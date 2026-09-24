!> GPU cuFFT/cuSPARSE Poisson solve. The fully periodic (x,y,z) branch is pencil-decomposed across MPI ranks/GPUs (multi-GPU capable); all other branches are single-GPU (nprocs==1 only)
Module poisson_gpu

  Use iso_fortran_env, Only : Int32, Int64
  Use cufft
  Use cusparse
  Use decomp, Only : decomp_poisson, transpose_x_to_y, transpose_y_to_x, transpose_y_to_z, transpose_z_to_y
  Use global, Only : nxp_global, nzp_global, mx, mz, nyg, rhs_p, rhs_p_hat, Dyy, kxx, kyy, kzz, x_bc_type, y_bc_type, pi, &
                     z_bc_type, nzm_global, Qz, sqrt_w_z, lambda_z, &
                     poisson_y_r, poisson_x_r, poisson_x_c, poisson_y_c, poisson_z_c

  Implicit None

  Logical :: plans_created = .False.
  Integer :: nslabs   ! number of interior y-planes batched per transform (nyg-2)


  ! Genuine 3D periodic cuFFT transform (y_bc_type==0 AND x_bc_type==0 only)
  ! (pencil-decomposed: batched 1-D cuFFT plans, one per pencil orientation, wired between
  ! 2decomp&fft's GPU (CUDA-aware) transposes exactly like projection.f90's CPU FFTW path)
  Integer :: plan_p3_x, plan_p3_y, plan_p3_z

  ! DCT-IV (x_bc_type==1) state: zero-padded 4N-point Z2Z FFT + twiddle, then z-FFT
  Integer :: plan_dct_L, plan_z, Lx
  Complex(Int64), Allocatable :: dct_ext(:,:,:), zwork(:,:,:)

  ! Batched cuSPARSE tridiagonal solve state
  Type(cusparseHandle) :: cusparse_h
  Logical :: gtsv_created = .False.
  Integer :: gtsv_m, gtsv_batch
  Complex(Int64), Allocatable :: gtsv_dl(:), gtsv_d(:), gtsv_du(:), gtsv_x(:)
  Character(1), Allocatable :: gtsv_buf(:)

  ! 4-wall duct (z_bc_type==1 .And. y_bc_type==1) state: x FFT (x_bc_type==0 only, so
  ! far) + device-resident z-eigenbasis (Qz/sqrt_w_z/lambda_z, built once on the host
  ! at init from Dzz -- see initialization.f90) applied per x-mode, then a batched
  ! y-tridiagonal cuSPARSE solve per (x-mode, z-eigenmode) pair -- mirrors
  ! gpu_solve_tridiagonal_batched but with lambda_z(m) standing in for kzz(k)
  Integer :: plan_duct_x_fwd, plan_duct_x_bwd
  Complex(Int64), Allocatable :: duct_c(:,:,:)   ! (nx1, nyg-2, nzm_global): x-mode, y, z (physical, then eigenmode)
  Real(Int64), Allocatable :: Qz_gpu(:,:), sqrt_w_z_gpu(:), lambda_z_gpu(:)

Contains

  !> Create the batched cuFFT plans once, on first use; branches on x_bc_type
  Subroutine gpu_poisson_init

    Integer :: ierr, nx1, nz1, ny_i, nyp_l

    nx1    = Int(nxp_global)
    nz1    = Int(nzp_global)
    nslabs = nyg - 2

    ! kxx/kzz are set once in initialization.f90 and never touched again, so a
    ! one-time upload here (rather than gpu_solve_periodic_3d/
    ! gpu_solve_tridiagonal_batched re-copying them from host on every call)
    ! is all they ever need.
    !$acc enter data create(kxx,kzz)
    !$acc update device(kxx,kzz)

    If ( z_bc_type == 1 ) Then

       ! Case (a), 4-wall duct only; case (b) (z wall alone, y periodic) is not yet
       ! GPU-ported -- guarded at init time already (initialization.f90), Stop here
       ! too as a defensive check against this module being reached any other way
       If ( y_bc_type /= 1 .Or. x_bc_type /= 0 ) &
            Stop 'ERROR: GPU_POISSON z_bc_type=1 only implemented for y_bc_type=1, x_bc_type=0 (4-wall duct, periodic x)'

       ny_i = nyg - 2
       nz1  = Int(nzm_global)

       Allocate( duct_c(nx1, ny_i, nz1) )
       !$acc enter data create(duct_c)

       ierr =         cufftPlanMany( plan_duct_x_fwd, 1, [nx1], [nx1], 1, nx1, [nx1], 1, nx1, CUFFT_Z2Z, ny_i*nz1 )
       ierr = ierr + cufftPlanMany( plan_duct_x_bwd, 1, [nx1], [nx1], 1, nx1, [nx1], 1, nx1, CUFFT_Z2Z, ny_i*nz1 )
       If ( ierr /= 0 ) Stop 'ERROR: cuFFT duct x-transform plan creation failed'

       ! Qz/sqrt_w_z/lambda_z are built once on the host in initialization.f90 (tiny
       ! LAPACK dstev call, nzm_global x nzm_global) -- just upload them here
       Allocate( Qz_gpu(nzm_global,nzm_global), sqrt_w_z_gpu(nzm_global), lambda_z_gpu(nzm_global) )
       Qz_gpu       = Qz
       sqrt_w_z_gpu = sqrt_w_z
       lambda_z_gpu = lambda_z
       !$acc enter data create(Qz_gpu,sqrt_w_z_gpu,lambda_z_gpu)
       !$acc update device(Qz_gpu,sqrt_w_z_gpu,lambda_z_gpu)

    Else If ( x_bc_type == 0 ) Then

       ! x,z periodic (y periodic or walls): pencil-decomposed, multi-GPU capable
       ! kyy only exists (is Allocated) in the periodic-y case
       If ( y_bc_type == 0 ) Then
          !$acc enter data create(kyy)
          !$acc update device(kyy)
       End If

       ! Periodic y transforms only the first nyp=ysz(2)-1 y-slots ("last interior cell is a
       ! redundant duplicate of the first", as for periodic x/z); the excluded last cell
       ! (rhs_p y-index nyg-1) is restored by a plain copy in projection.f90.
       nyp_l = decomp_poisson%ysz(2) - 1
       !$acc enter data create(poisson_y_r,poisson_x_r,poisson_x_c,poisson_y_c,poisson_z_c)

       ! x-pencil: contiguous length-nxp transforms, batch over (j,k)
       ierr =         cufftPlanMany( plan_p3_x, 1, [nx1], [nx1], 1, nx1, [nx1], 1, nx1, CUFFT_Z2Z, &
                                     decomp_poisson%xsz(2)*decomp_poisson%xsz(3) )
       ! z-pencil: length-nzp transforms at stride zsz1*zsz2, batch over (i,j)
       ierr = ierr + cufftPlanMany( plan_p3_z, 1, [nz1], [nz1], decomp_poisson%zsz(1)*decomp_poisson%zsz(2), 1, &
                                    [nz1], decomp_poisson%zsz(1)*decomp_poisson%zsz(2), 1, CUFFT_Z2Z, &
                                    decomp_poisson%zsz(1)*decomp_poisson%zsz(2) )
       ! y-pencil (periodic y only): length-nyp transforms at stride ysz1, batch over i (one plan execution per k-slice)
       If ( y_bc_type == 0 ) Then
          ierr = ierr + cufftPlanMany( plan_p3_y, 1, [nyp_l], [nyp_l], decomp_poisson%ysz(1), 1, &
                                       [nyp_l], decomp_poisson%ysz(1), 1, CUFFT_Z2Z, decomp_poisson%ysz(1) )
       End If
       If ( ierr /= 0 ) Stop 'ERROR: cuFFT pencil periodic plan creation failed'

    Else

       ! dct_ext is x-fastest (batch over z,yslab); zwork is z-fastest (batch over mode,yslab)
       Lx = 4*nx1
       Allocate( dct_ext(Lx,nz1,nslabs), zwork(nz1,nx1,nslabs) )
       !$acc enter data create(dct_ext,zwork)

       ierr = cufftPlanMany( plan_dct_L, 1, [Lx], [Lx], 1, Lx, [Lx], 1, Lx, &
                             CUFFT_Z2Z, nz1*nslabs )
       ierr = ierr + cufftPlanMany( plan_z, 1, [nz1], [nz1], 1, nz1, [nz1], 1, nz1, &
                                    CUFFT_Z2Z, nx1*nslabs )
       If ( ierr /= 0 ) Stop 'ERROR: cuFFT DCT-IV batched plan creation failed'

    End If

    plans_created = .True.

  End Subroutine gpu_poisson_init

  !> Pencil-decomposed forward transform, fully periodic case: y(real) -> x (cuFFT) -> y -> z (cuFFT) -> y, leaving poisson_y_c in (kx,y,kz) space; mirrors projection.f90's CPU chain with device-resident arrays and CUDA-aware transposes
  Subroutine gpu_forward_transform_3d

    Integer :: i, j, k, ierr, n1, n2, n3

    If ( .Not. plans_created ) Call gpu_poisson_init

    n1 = decomp_poisson%ysz(1); n2 = decomp_poisson%ysz(2); n3 = decomp_poisson%ysz(3)

    !$acc parallel loop collapse(3) present(rhs_p,poisson_y_r)
    Do k = 1, n3
       Do j = 1, n2
          Do i = 1, n1
             poisson_y_r(i,j,k) = rhs_p(i+1,j+1,k+1)
          End Do
       End Do
    End Do
    !$acc end parallel loop

    Call transpose_y_to_x( poisson_y_r, poisson_x_r, decomp_poisson )

    n1 = decomp_poisson%xsz(1); n2 = decomp_poisson%xsz(2); n3 = decomp_poisson%xsz(3)
    !$acc parallel loop collapse(3) present(poisson_x_r,poisson_x_c)
    Do k = 1, n3
       Do j = 1, n2
          Do i = 1, n1
             poisson_x_c(i,j,k) = dcmplx( poisson_x_r(i,j,k) )
          End Do
       End Do
    End Do
    !$acc end parallel loop
    !$acc host_data use_device(poisson_x_c)
    ierr = cufftExecZ2Z( plan_p3_x, poisson_x_c, poisson_x_c, CUFFT_FORWARD )
    !$acc end host_data
    If ( ierr /= 0 ) Stop 'ERROR: cuFFT pencil x forward exec failed'
    Call transpose_x_to_y( poisson_x_c, poisson_y_c, decomp_poisson )

    Call transpose_y_to_z( poisson_y_c, poisson_z_c, decomp_poisson )
    !$acc host_data use_device(poisson_z_c)
    ierr = cufftExecZ2Z( plan_p3_z, poisson_z_c, poisson_z_c, CUFFT_FORWARD )
    !$acc end host_data
    If ( ierr /= 0 ) Stop 'ERROR: cuFFT pencil z forward exec failed'
    Call transpose_z_to_y( poisson_z_c, poisson_y_c, decomp_poisson )

  End Subroutine gpu_forward_transform_3d

  !> y-FFT, elementwise divide by (kxx+kyy+kzz) in Fourier space, inverse y-FFT on this rank's y-pencil; (0,0,0) mode fixed to zero (periodic Poisson null-space), owned by whichever rank holds global mode (0,0,0)
  Subroutine gpu_solve_periodic_3d

    Integer :: i, j, k, k_global, i_global, ierr, nyp, n1, n3, yst1, yst3
    Real(Int64) :: inv_nyp

    n1   = decomp_poisson%ysz(1)
    n3   = decomp_poisson%ysz(3)
    nyp  = decomp_poisson%ysz(2) - 1
    yst1 = decomp_poisson%yst(1)
    yst3 = decomp_poisson%yst(3)
    inv_nyp = 1d0 / Real(nyp,Int64)

    Do k = 1, n3
       !$acc host_data use_device(poisson_y_c)
       ierr = cufftExecZ2Z( plan_p3_y, poisson_y_c(:,:,k), poisson_y_c(:,:,k), CUFFT_FORWARD )
       !$acc end host_data
       If ( ierr /= 0 ) Stop 'ERROR: cuFFT pencil y forward exec failed'
    End Do

    ! kxx/kyy/kzz are made persistently device-resident once in gpu_poisson_init
    !$acc parallel loop collapse(3) present(kxx,kyy,kzz,poisson_y_c) private(i_global,k_global)
    Do k = 1, n3
       Do j = 1, nyp
          Do i = 1, n1
             i_global = yst1 + i - 2
             k_global = yst3 + k - 2
             If ( i_global==0 .And. j==1 .And. k_global==0 ) Then
                poisson_y_c(i,j,k) = (0d0,0d0)
             Else
                poisson_y_c(i,j,k) = poisson_y_c(i,j,k) / ( kxx(i_global) + kyy(j-1) + kzz(k_global) )
             End If
          End Do
       End Do
    End Do
    !$acc end parallel loop

    Do k = 1, n3
       !$acc host_data use_device(poisson_y_c)
       ierr = cufftExecZ2Z( plan_p3_y, poisson_y_c(:,:,k), poisson_y_c(:,:,k), CUFFT_INVERSE )
       !$acc end host_data
       If ( ierr /= 0 ) Stop 'ERROR: cuFFT pencil y inverse exec failed'
    End Do

    !$acc parallel loop collapse(3) present(poisson_y_c)
    Do k = 1, n3
       Do j = 1, nyp
          Do i = 1, n1
             poisson_y_c(i,j,k) = poisson_y_c(i,j,k) * inv_nyp
          End Do
       End Do
    End Do
    !$acc end parallel loop

  End Subroutine gpu_solve_periodic_3d

  !> Pencil-decomposed inverse transform, fully periodic case: y -> z (inverse cuFFT) -> y -> x (inverse cuFFT) -> y(real) -> rhs_p, normalised by nxp*nzp (the y factor is applied in gpu_solve_periodic_3d)
  Subroutine gpu_inverse_transform_3d

    Integer :: i, j, k, ierr, n1, n2, n3
    Real(Int64) :: inv_norm

    Call transpose_y_to_z( poisson_y_c, poisson_z_c, decomp_poisson )
    !$acc host_data use_device(poisson_z_c)
    ierr = cufftExecZ2Z( plan_p3_z, poisson_z_c, poisson_z_c, CUFFT_INVERSE )
    !$acc end host_data
    If ( ierr /= 0 ) Stop 'ERROR: cuFFT pencil z inverse exec failed'
    Call transpose_z_to_y( poisson_z_c, poisson_y_c, decomp_poisson )

    Call transpose_y_to_x( poisson_y_c, poisson_x_c, decomp_poisson )
    !$acc host_data use_device(poisson_x_c)
    ierr = cufftExecZ2Z( plan_p3_x, poisson_x_c, poisson_x_c, CUFFT_INVERSE )
    !$acc end host_data
    If ( ierr /= 0 ) Stop 'ERROR: cuFFT pencil x inverse exec failed'

    inv_norm = 1d0 / ( Real(nxp_global,Int64) * Real(nzp_global,Int64) )
    n1 = decomp_poisson%xsz(1); n2 = decomp_poisson%xsz(2); n3 = decomp_poisson%xsz(3)
    !$acc parallel loop collapse(3) present(poisson_x_c,poisson_x_r)
    Do k = 1, n3
       Do j = 1, n2
          Do i = 1, n1
             poisson_x_r(i,j,k) = Real( poisson_x_c(i,j,k), Int64 ) * inv_norm
          End Do
       End Do
    End Do
    !$acc end parallel loop
    Call transpose_x_to_y( poisson_x_r, poisson_y_r, decomp_poisson )

    n1 = decomp_poisson%ysz(1); n2 = decomp_poisson%ysz(2); n3 = decomp_poisson%ysz(3)
    !$acc parallel loop collapse(3) present(poisson_y_r,rhs_p)
    Do k = 1, n3
       Do j = 1, n2
          Do i = 1, n1
             rhs_p(i+1,j+1,k+1) = poisson_y_r(i,j,k)
          End Do
       End Do
    End Do
    !$acc end parallel loop

  End Subroutine gpu_inverse_transform_3d

  !> Forward DCT-IV(x) then complex FFT(z) of every y-slab of rhs_p into rhs_p_hat, batched; x_bc_type==1 only, direct kx=i/kz=k indexing
  Subroutine gpu_forward_transform_dct_slabs

    Integer :: j, ix, iz, imode, ierr, nx1, nz1
    Real(Int64) :: theta
    Complex(Int64) :: w

    If ( .Not. plans_created ) Call gpu_poisson_init

    nx1 = Int(nxp_global)
    nz1 = Int(nzp_global)

    ! zero-padded load, x-fastest, batched over (z-line,yslab)
    !$acc parallel loop collapse(3) present(rhs_p,dct_ext)
    Do j = 2, nyg-1
       Do iz = 1, nz1
          Do ix = 1, Lx
             If ( ix <= nx1 ) Then
                dct_ext(ix,iz,j-1) = dcmplx( rhs_p(ix+1,j,iz+1) )
             Else
                dct_ext(ix,iz,j-1) = (0d0,0d0)
             End If
          End Do
       End Do
    End Do
    !$acc end parallel loop

    !$acc host_data use_device(dct_ext)
    ierr = cufftExecZ2Z( plan_dct_L, dct_ext, dct_ext, CUFFT_FORWARD )
    !$acc end host_data
    If ( ierr /= 0 ) Stop 'ERROR: cuFFT DCT-IV forward exec failed'

    ! extract DCT-IV coeffs via twiddle (S[k]=X4N(4N-2k-1), index Lx-2*imode), transpose into zwork
    !$acc parallel loop collapse(3) present(dct_ext,zwork)
    Do j = 2, nyg-1
       Do iz = 1, nz1
          Do imode = 0, mx
             theta = pi*Real(2*imode+1,Int64)/Real(4*nx1,Int64)
             w = dcmplx( dcos(theta), dsin(theta) )
             zwork(iz,imode+1,j-1) = 2d0*Real( w*dct_ext(Lx-2*imode,iz,j-1), Int64 )
          End Do
       End Do
    End Do
    !$acc end parallel loop

    !$acc host_data use_device(zwork)
    ierr = cufftExecZ2Z( plan_z, zwork, zwork, CUFFT_FORWARD )
    !$acc end host_data
    If ( ierr /= 0 ) Stop 'ERROR: cuFFT DCT-IV z-direction forward exec failed'

    !$acc parallel loop collapse(3) present(zwork,rhs_p_hat)
    Do j = 2, nyg-1
       Do imode = 0, mx
          Do iz = 1, nz1
             rhs_p_hat(j,imode,iz-1) = zwork(iz,imode+1,j-1)
          End Do
       End Do
    End Do
    !$acc end parallel loop

  End Subroutine gpu_forward_transform_dct_slabs

  !> Inverse of gpu_forward_transform_dct_slabs: inverse z-FFT then self-inverse DCT-IV(x), normalised by 2*nxp_global*nzp_global; x_bc_type==1 only
  Subroutine gpu_inverse_transform_dct_slabs

    Integer :: j, ix, iz, imode, ierr, nx1, nz1
    Real(Int64) :: theta, norm
    Complex(Int64) :: w

    nx1  = Int(nxp_global)
    nz1  = Int(nzp_global)
    norm = 2d0*Real(nxp_global,Int64)*Real(nzp_global,Int64)

    !$acc parallel loop collapse(3) present(rhs_p_hat,zwork)
    Do j = 2, nyg-1
       Do imode = 0, mx
          Do iz = 1, nz1
             zwork(iz,imode+1,j-1) = rhs_p_hat(j,imode,iz-1)
          End Do
       End Do
    End Do
    !$acc end parallel loop

    !$acc host_data use_device(zwork)
    ierr = cufftExecZ2Z( plan_z, zwork, zwork, CUFFT_INVERSE )
    !$acc end host_data
    If ( ierr /= 0 ) Stop 'ERROR: cuFFT DCT-IV z-direction inverse exec failed'

    ! zero-pad the (real, up to roundoff) inverse-z-FFT result back into x-fastest dct_ext
    !$acc parallel loop collapse(3) present(zwork,dct_ext)
    Do j = 2, nyg-1
       Do iz = 1, nz1
          Do ix = 1, Lx
             If ( ix <= nx1 ) Then
                dct_ext(ix,iz,j-1) = dcmplx( Real(zwork(iz,ix,j-1),Int64) )
             Else
                dct_ext(ix,iz,j-1) = (0d0,0d0)
             End If
          End Do
       End Do
    End Do
    !$acc end parallel loop

    !$acc host_data use_device(dct_ext)
    ierr = cufftExecZ2Z( plan_dct_L, dct_ext, dct_ext, CUFFT_FORWARD )
    !$acc end host_data
    If ( ierr /= 0 ) Stop 'ERROR: cuFFT DCT-IV inverse(self) exec failed'

    !$acc parallel loop collapse(3) present(dct_ext,rhs_p)
    Do j = 2, nyg-1
       Do iz = 1, nz1
          Do imode = 0, mx
             theta = pi*Real(2*imode+1,Int64)/Real(4*nx1,Int64)
             w = dcmplx( dcos(theta), dsin(theta) )
             rhs_p(imode+2,j,iz+1) = 2d0*Real( w*dct_ext(Lx-2*imode,iz,j-1), Int64 ) / norm
          End Do
       End Do
    End Do
    !$acc end parallel loop

  End Subroutine gpu_inverse_transform_dct_slabs

  !> Create the cuSPARSE handle and batched tridiagonal-solve arrays once, on first use
  Subroutine gpu_gtsv_init

    Integer :: ierr
    Integer(Int64) :: bufsize

    gtsv_m     = nyg - 2
    If ( x_bc_type == 0 .And. z_bc_type == 0 ) Then
       ! pencil path: one system per (kx,kz) mode this rank owns in its y-pencil
       gtsv_batch = decomp_poisson%ysz(1) * decomp_poisson%ysz(3)
    Else
       gtsv_batch = ( Int(mx,Int32) + 1 ) * ( Int(mz,Int32) + 1 )
    End If

    Allocate( gtsv_dl(gtsv_m*gtsv_batch), gtsv_d(gtsv_m*gtsv_batch), &
              gtsv_du(gtsv_m*gtsv_batch), gtsv_x(gtsv_m*gtsv_batch) )
    ! placeholder values for the sizing query below (bufferSize depends only on m/batchCount)
    gtsv_dl = (0d0,0d0)
    gtsv_d  = (1d0,0d0)
    gtsv_du = (0d0,0d0)
    gtsv_x  = (0d0,0d0)

    ierr = cusparseCreate( cusparse_h )
    If ( ierr /= 0 ) Stop 'ERROR: cusparseCreate failed'

    ! gtsv_dl/d/du/x are pure per-call scratch (fully overwritten before being
    ! read each call in gpu_solve_tridiagonal_batched), so a one-time device
    ! allocation here -- instead of create/destroy on every call -- is safe;
    ! the placeholder values only need to be on-device for this sizing query.
    !$acc enter data create(gtsv_dl,gtsv_d,gtsv_du,gtsv_x)
    !$acc update device(gtsv_dl,gtsv_d,gtsv_du,gtsv_x)
    !$acc host_data use_device(gtsv_dl,gtsv_d,gtsv_du,gtsv_x)
    ierr = cusparseZgtsvInterleavedBatch_bufferSizeExt( cusparse_h, CUSPARSE_ALG1, gtsv_m, &
                 gtsv_dl, gtsv_d, gtsv_du, gtsv_x, gtsv_batch, bufsize )
    !$acc end host_data
    If ( ierr /= 0 ) Stop 'ERROR: cusparseZgtsvInterleavedBatch_bufferSizeExt failed'

    Allocate( gtsv_buf(bufsize) )
    !$acc enter data create(gtsv_buf)

    gtsv_created = .True.

  End Subroutine gpu_gtsv_init

  !> Pack rhs_p_hat, batch-solve all y-tridiagonal systems via cuSPARSE, unpack back
  Subroutine gpu_solve_tridiagonal_batched

    Integer :: imode, k, b, ii, j, idx, ierr, mx_i, mz_i

    If ( .Not. gtsv_created ) Call gpu_gtsv_init

    mx_i = Int(mx,Int32)
    mz_i = Int(mz,Int32)

    ! rhs_p_hat/Dyy/kxx/kzz/gtsv_* are all persistently device-resident
    ! (initialization.f90 for rhs_p_hat/Dyy, gpu_poisson_init for kxx/kzz,
    ! gpu_gtsv_init for gtsv_*; Dyy's evolving Robin-BC elements are kept
    ! current via wallmodel.f90's own targeted !$acc update device) --
    ! present() fails loudly instead of silently re-transferring if that
    ! assumption ever breaks
    !$acc data present(Dyy,kxx,kzz,rhs_p_hat,gtsv_dl,gtsv_d,gtsv_du,gtsv_x,gtsv_buf)

    !$acc parallel loop collapse(2) present(Dyy,kxx,kzz,rhs_p_hat,gtsv_dl,gtsv_d,gtsv_du,gtsv_x)
    Do k = 0, mz_i
       Do imode = 0, mx_i
          b = k*(mx_i+1) + imode
          Do ii = 0, gtsv_m-1
             j   = ii + 2
             idx = ii*gtsv_batch + b + 1
             gtsv_d(idx) = Dyy(j,j) + kxx(imode) + kzz(k)
             If ( ii > 0 ) Then
                gtsv_dl(idx) = Dyy(j,j-1)
             Else
                gtsv_dl(idx) = (0d0,0d0)
             End If
             If ( ii < gtsv_m-1 ) Then
                gtsv_du(idx) = Dyy(j,j+1)
             Else
                gtsv_du(idx) = (0d0,0d0)
             End If
             gtsv_x(idx) = rhs_p_hat(j,imode,k)
          End Do
       End Do
    End Do
    !$acc end parallel loop

    ! Remove singularity 00 mode (periodic case only): ii=0,b=0 -> idx=1
    !$acc kernels present(gtsv_d)
    If ( x_bc_type == 0 ) gtsv_d(1) = 3d0/2d0*gtsv_d(1)
    !$acc end kernels

    !$acc host_data use_device(gtsv_dl,gtsv_d,gtsv_du,gtsv_x,gtsv_buf)
    ierr = cusparseZgtsvInterleavedBatch( cusparse_h, CUSPARSE_ALG1, gtsv_m, &
                 gtsv_dl, gtsv_d, gtsv_du, gtsv_x, gtsv_batch, gtsv_buf )
    !$acc end host_data
    If ( ierr /= 0 ) Stop 'ERROR: cusparseZgtsvInterleavedBatch failed'

    !$acc parallel loop collapse(2) present(rhs_p_hat,gtsv_x)
    Do k = 0, mz_i
       Do imode = 0, mx_i
          b = k*(mx_i+1) + imode
          Do ii = 0, gtsv_m-1
             j   = ii + 2
             idx = ii*gtsv_batch + b + 1
             rhs_p_hat(j,imode,k) = gtsv_x(idx)
          End Do
       End Do
    End Do
    !$acc end parallel loop

    !$acc end data

  End Subroutine gpu_solve_tridiagonal_batched

  !> Wall-normal (y) solve for the pencil path with x,z periodic and y walls: batched cuSPARSE tridiagonal solve, one system per (kx,kz) mode this rank owns, operating in place on the y-pencil poisson_y_c; global mode indices come from decomp_poisson%yst
  Subroutine gpu_solve_tridiagonal_pencil

    Integer :: i, k, ii, j, idx, b, ierr, n1, n3, yst1, yst3, i_global, k_global

    If ( .Not. gtsv_created ) Call gpu_gtsv_init

    n1   = decomp_poisson%ysz(1)
    n3   = decomp_poisson%ysz(3)
    yst1 = decomp_poisson%yst(1)
    yst3 = decomp_poisson%yst(3)

    !$acc data present(Dyy,kxx,kzz,poisson_y_c,gtsv_dl,gtsv_d,gtsv_du,gtsv_x,gtsv_buf)

    !$acc parallel loop collapse(2) present(Dyy,kxx,kzz,poisson_y_c,gtsv_dl,gtsv_d,gtsv_du,gtsv_x) private(i_global,k_global,b,ii,j,idx)
    Do k = 1, n3
       Do i = 1, n1
          i_global = yst1 + i - 2
          k_global = yst3 + k - 2
          b = (k-1)*n1 + (i-1)
          Do ii = 0, gtsv_m-1
             j   = ii + 2
             idx = ii*gtsv_batch + b + 1
             gtsv_d(idx) = Dyy(j,j) + kxx(i_global) + kzz(k_global)
             ! remove the periodic null-space singularity of the (0,0) mode, whichever rank owns it
             If ( ii == 0 .And. i_global == 0 .And. k_global == 0 ) gtsv_d(idx) = 3d0/2d0*gtsv_d(idx)
             If ( ii > 0 ) Then
                gtsv_dl(idx) = Dyy(j,j-1)
             Else
                gtsv_dl(idx) = (0d0,0d0)
             End If
             If ( ii < gtsv_m-1 ) Then
                gtsv_du(idx) = Dyy(j,j+1)
             Else
                gtsv_du(idx) = (0d0,0d0)
             End If
             gtsv_x(idx) = poisson_y_c(i,ii+1,k)
          End Do
       End Do
    End Do
    !$acc end parallel loop

    !$acc host_data use_device(gtsv_dl,gtsv_d,gtsv_du,gtsv_x,gtsv_buf)
    ierr = cusparseZgtsvInterleavedBatch( cusparse_h, CUSPARSE_ALG1, gtsv_m, &
                 gtsv_dl, gtsv_d, gtsv_du, gtsv_x, gtsv_batch, gtsv_buf )
    !$acc end host_data
    If ( ierr /= 0 ) Stop 'ERROR: cusparseZgtsvInterleavedBatch (pencil) failed'

    !$acc parallel loop collapse(2) present(poisson_y_c,gtsv_x) private(b,ii,idx)
    Do k = 1, n3
       Do i = 1, n1
          b = (k-1)*n1 + (i-1)
          Do ii = 0, gtsv_m-1
             idx = ii*gtsv_batch + b + 1
             poisson_y_c(i,ii+1,k) = gtsv_x(idx)
          End Do
       End Do
    End Do
    !$acc end parallel loop

    !$acc end data

  End Subroutine gpu_solve_tridiagonal_pencil

  !> Forward x-FFT (x_bc_type==0 only) then z-eigenmode transform of every interior point of rhs_p into rhs_p_hat, batched; z_bc_type==1 .And. y_bc_type==1 (4-wall duct) only. Fully device-resident.
  Subroutine gpu_forward_transform_duct

    Integer :: ix, jy, iz, i, m, k, ierr, nx1, ny_i, nz1
    Complex(Int64) :: acc

    If ( .Not. plans_created ) Call gpu_poisson_init

    nx1  = Int(nxp_global)
    ny_i = nyg - 2
    nz1  = Int(nzm_global)

    !$acc parallel loop collapse(3) present(rhs_p,duct_c)
    Do iz = 1, nz1
       Do jy = 1, ny_i
          Do ix = 1, nx1
             duct_c(ix,jy,iz) = dcmplx( rhs_p(ix+1,jy+1,iz+1) )
          End Do
       End Do
    End Do
    !$acc end parallel loop

    !$acc host_data use_device(duct_c)
    ierr = cufftExecZ2Z( plan_duct_x_fwd, duct_c, duct_c, CUFFT_FORWARD )
    !$acc end host_data
    If ( ierr /= 0 ) Stop 'ERROR: cuFFT duct x-forward exec failed'

    ! z-eigenmode transform per x-mode (mirrors the CPU path's Matmul(z_hat,Qz) in
    ! projection.f90, scaled first by sqrt_w_z to recover Dzz's own eigenbasis from
    ! Qz's orthonormal one -- see global.f90's lambda_z/Qz/sqrt_w_z comment)
    !$acc parallel loop collapse(3) private(acc) present(duct_c,Qz_gpu,sqrt_w_z_gpu,rhs_p_hat)
    Do i = 0, Int(mx,Int32)
       Do jy = 1, ny_i
          Do m = 0, nz1-1
             acc = (0d0,0d0)
             Do k = 1, nz1
                acc = acc + duct_c(i+1,jy,k) * sqrt_w_z_gpu(k) * Qz_gpu(k,m+1)
             End Do
             rhs_p_hat(jy+1, i, m) = acc
          End Do
       End Do
    End Do
    !$acc end parallel loop

  End Subroutine gpu_forward_transform_duct

  !> Inverse of gpu_forward_transform_duct: inverse z-eigenmode transform then inverse x-FFT, normalised by nxp_global; z_bc_type==1 .And. y_bc_type==1 only
  Subroutine gpu_inverse_transform_duct

    Integer :: ix, jy, iz, i, m, k, ierr, nx1, ny_i, nz1
    Real(Int64) :: norm
    Complex(Int64) :: acc

    nx1  = Int(nxp_global)
    ny_i = nyg - 2
    nz1  = Int(nzm_global)
    norm = Real(nx1,Int64)

    ! inverse z-eigenmode transform (Qz orthonormal -> inverse == transpose), then
    ! unscale by sqrt_w_z to recover Dzz's own eigenbasis
    !$acc parallel loop collapse(3) private(acc) present(duct_c,Qz_gpu,sqrt_w_z_gpu,rhs_p_hat)
    Do i = 0, Int(mx,Int32)
       Do jy = 1, ny_i
          Do k = 1, nz1
             acc = (0d0,0d0)
             Do m = 0, nz1-1
                acc = acc + rhs_p_hat(jy+1,i,m) * Qz_gpu(k,m+1)
             End Do
             duct_c(i+1,jy,k) = acc / sqrt_w_z_gpu(k)
          End Do
       End Do
    End Do
    !$acc end parallel loop

    !$acc host_data use_device(duct_c)
    ierr = cufftExecZ2Z( plan_duct_x_bwd, duct_c, duct_c, CUFFT_INVERSE )
    !$acc end host_data
    If ( ierr /= 0 ) Stop 'ERROR: cuFFT duct x-inverse exec failed'

    !$acc parallel loop collapse(3) present(duct_c,rhs_p)
    Do iz = 1, nz1
       Do jy = 1, ny_i
          Do ix = 1, nx1
             rhs_p(ix+1,jy+1,iz+1) = Real(duct_c(ix,jy,iz),Int64) / norm
          End Do
       End Do
    End Do
    !$acc end parallel loop

  End Subroutine gpu_inverse_transform_duct

  !> Pack rhs_p_hat (already forward z-eigenmode-transformed by gpu_forward_transform_duct), batch-solve one y-tridiagonal system per (x-mode, z-eigenmode) pair via cuSPARSE, unpack back; z_bc_type==1 .And. y_bc_type==1 (4-wall duct) only -- lambda_z(m) stands in for kzz(k) in gpu_solve_tridiagonal_batched
  Subroutine gpu_solve_duct_tridiagonal_batched

    Integer :: imode, m, b, ii, j, idx, ierr, mx_i, null_m

    If ( .Not. gtsv_created ) Call gpu_gtsv_init

    mx_i = Int(mx,Int32)

    ! null_m: 0-based index of the z-eigenmode nearest zero (Dzz's null mode, always
    ! the LAST/ascending-largest eigenvalue of lambda_z -- see initialization.f90);
    ! found on the host once per call, cheap (nzm_global small)
    null_m = Maxloc( lambda_z, dim=1 ) - 1

    !$acc data present(Dyy,kxx,lambda_z_gpu,rhs_p_hat,gtsv_dl,gtsv_d,gtsv_du,gtsv_x,gtsv_buf)

    !$acc parallel loop collapse(2) present(Dyy,kxx,lambda_z_gpu,rhs_p_hat,gtsv_dl,gtsv_d,gtsv_du,gtsv_x)
    Do m = 0, nzm_global-1
       Do imode = 0, mx_i
          b = m*(mx_i+1) + imode
          Do ii = 0, gtsv_m-1
             j   = ii + 2
             idx = ii*gtsv_batch + b + 1
             gtsv_d(idx) = Dyy(j,j) + kxx(imode) + lambda_z_gpu(m+1)
             If ( ii > 0 ) Then
                gtsv_dl(idx) = Dyy(j,j-1)
             Else
                gtsv_dl(idx) = (0d0,0d0)
             End If
             If ( ii < gtsv_m-1 ) Then
                gtsv_du(idx) = Dyy(j,j+1)
             Else
                gtsv_du(idx) = (0d0,0d0)
             End If
             gtsv_x(idx) = rhs_p_hat(j,imode,m)
          End Do
       End Do
    End Do
    !$acc end parallel loop

    ! Remove singularity of the 00 mode (x_bc_type==0's zero x-wavenumber paired with
    ! z's null eigenmode); ii=0,imode=0,m=null_m -> idx=null_m*(mx_i+1)+1
    !$acc kernels present(gtsv_d)
    If ( x_bc_type == 0 ) gtsv_d( null_m*(mx_i+1) + 1 ) = 3d0/2d0*gtsv_d( null_m*(mx_i+1) + 1 )
    !$acc end kernels

    !$acc host_data use_device(gtsv_dl,gtsv_d,gtsv_du,gtsv_x,gtsv_buf)
    ierr = cusparseZgtsvInterleavedBatch( cusparse_h, CUSPARSE_ALG1, gtsv_m, &
                 gtsv_dl, gtsv_d, gtsv_du, gtsv_x, gtsv_batch, gtsv_buf )
    !$acc end host_data
    If ( ierr /= 0 ) Stop 'ERROR: cusparseZgtsvInterleavedBatch (duct) failed'

    !$acc parallel loop collapse(2) present(rhs_p_hat,gtsv_x)
    Do m = 0, nzm_global-1
       Do imode = 0, mx_i
          b = m*(mx_i+1) + imode
          Do ii = 0, gtsv_m-1
             j   = ii + 2
             idx = ii*gtsv_batch + b + 1
             rhs_p_hat(j,imode,m) = gtsv_x(idx)
          End Do
       End Do
    End Do
    !$acc end parallel loop

    !$acc end data

  End Subroutine gpu_solve_duct_tridiagonal_batched

End Module poisson_gpu
