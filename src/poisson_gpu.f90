!> GPU cuFFT/cuSPARSE Poisson solve, pencil-decomposed across MPI ranks/GPUs (multi-GPU capable) for every supported BC combination: x periodic (FFT) or inflow/outflow (DCT-IV), z periodic (FFT) or walls (4-wall duct eigenmode solve), y periodic (FFT) or walls (tridiagonal)
Module poisson_gpu

  Use iso_fortran_env, Only : Int32, Int64
  Use cufft
  Use cusparse
  Use decomp, Only : decomp_poisson, transpose_x_to_y, transpose_y_to_x, transpose_y_to_z, transpose_z_to_y
  Use global, Only : nxp_global, nzp_global, nyg, rhs_p, Dyy, kxx, kyy, kzz, x_bc_type, y_bc_type, pi, &
                     z_bc_type, nzm_global, Qz, sqrt_w_z, lambda_z, &
                     poisson_y_r, poisson_x_r, poisson_x_c, poisson_y_c, poisson_z_c

  Implicit None

  Logical :: plans_created = .False.

  ! Genuine 3D periodic cuFFT transform (y_bc_type==0 AND x_bc_type==0 only)
  ! (pencil-decomposed: batched 1-D cuFFT plans, one per pencil orientation, wired between
  ! 2decomp&fft's GPU (CUDA-aware) transposes exactly like projection.f90's CPU FFTW path)
  Integer :: plan_p3_x, plan_p3_y, plan_p3_z

  ! DCT-IV (x_bc_type==1) state: zero-padded 4N-point Z2Z FFT of every x-line of this rank's x-pencil + twiddle
  Integer :: plan_dct_L, Lx
  Complex(Int64), Allocatable :: dct_ext(:,:,:)

  ! Batched cuSPARSE tridiagonal solve state
  Type(cusparseHandle) :: cusparse_h
  Logical :: gtsv_created = .False.
  Integer :: gtsv_m, gtsv_batch
  Complex(Int64), Allocatable :: gtsv_dl(:), gtsv_d(:), gtsv_du(:), gtsv_x(:)
  Character(1), Allocatable :: gtsv_buf(:)

  ! 4-wall duct (z_bc_type==1 .And. y_bc_type==1, x_bc_type==0) state: the x-FFT runs in the x-pencil like the
  ! periodic case; the z-eigenbasis (Qz/sqrt_w_z/lambda_z, built once on the host at init from Dzz -- see
  ! initialization.f90) is applied per local x-mode in the y-pencil (p_col==1: z fully local), then a batched
  ! y-tridiagonal cuSPARSE solve per (x-mode, z-eigenmode) pair
  Complex(Int64), Allocatable :: duct_c(:,:,:)   ! y-pencil shape: (local x-mode, y, z-eigenmode)
  Real(Int64), Allocatable :: Qz_gpu(:,:), sqrt_w_z_gpu(:), lambda_z_gpu(:)
  Real(Int64) :: lambda_z_tol

Contains

  !> Create the batched cuFFT plans once, on first use (pencil-decomposed: one batched 1-D plan per pencil orientation)
  Subroutine gpu_poisson_init

    Integer :: ierr, nx1, nz1, nyp_l

    nx1    = Int(nxp_global)
    nz1    = Int(nzp_global)

    ! kxx/kzz are set once in initialization.f90 and never touched again, so a
    ! one-time upload here (rather than re-copying them from host on every call)
    ! is all they ever need.
    !$acc enter data create(kxx,kzz)
    !$acc update device(kxx,kzz)

    If ( z_bc_type == 1 ) Then
       ! Case (a), 4-wall duct only; case (b) (z wall alone, y periodic) is not yet
       ! GPU-ported -- guarded at init time already (initialization.f90), Stop here
       ! too as a defensive check against this module being reached any other way
       If ( y_bc_type /= 1 .Or. x_bc_type /= 0 ) &
            Stop 'ERROR: GPU_POISSON z_bc_type=1 only implemented for y_bc_type=1, x_bc_type=0 (4-wall duct, periodic x)'
    End If

    ! kyy only exists (is Allocated) in the periodic-y case
    If ( y_bc_type == 0 ) Then
       !$acc enter data create(kyy)
       !$acc update device(kyy)
    End If

    ! Periodic y transforms only the first nyp=ysz(2)-1 y-slots ("last interior cell is a
    ! redundant duplicate of the first", as for periodic x/z); the excluded last cell
    ! (rhs_p y-index nyg-1) is restored by a plain copy in projection.f90.
    nyp_l = decomp_poisson%ysz(2) - 1
    !$acc enter data create(poisson_y_r,poisson_x_r,poisson_x_c,poisson_y_c)
    If ( z_bc_type == 0 ) Then
       !$acc enter data create(poisson_z_c)
    End If

    ierr = 0
    If ( x_bc_type == 0 ) Then
       ! x-pencil: contiguous length-nxp transforms, batch over (j,k)
       ierr = ierr + cufftPlanMany( plan_p3_x, 1, [nx1], [nx1], 1, nx1, [nx1], 1, nx1, CUFFT_Z2Z, &
                                    decomp_poisson%xsz(2)*decomp_poisson%xsz(3) )
    Else
       ! DCT-IV via zero-padded 4N FFT, batch over the x-pencil's (j,k) lines
       Lx = 4*nx1
       Allocate( dct_ext(Lx,decomp_poisson%xsz(2),decomp_poisson%xsz(3)) )
       !$acc enter data create(dct_ext)
       ierr = ierr + cufftPlanMany( plan_dct_L, 1, [Lx], [Lx], 1, Lx, [Lx], 1, Lx, CUFFT_Z2Z, &
                                    decomp_poisson%xsz(2)*decomp_poisson%xsz(3) )
    End If
    If ( z_bc_type == 0 ) Then
       ! z-pencil: length-nzp transforms at stride zsz1*zsz2, batch over (i,j)
       ierr = ierr + cufftPlanMany( plan_p3_z, 1, [nz1], [nz1], decomp_poisson%zsz(1)*decomp_poisson%zsz(2), 1, &
                                    [nz1], decomp_poisson%zsz(1)*decomp_poisson%zsz(2), 1, CUFFT_Z2Z, &
                                    decomp_poisson%zsz(1)*decomp_poisson%zsz(2) )
    End If
    ! y-pencil (periodic y only): length-nyp transforms at stride ysz1, batch over i (one plan execution per k-slice)
    If ( y_bc_type == 0 ) Then
       ierr = ierr + cufftPlanMany( plan_p3_y, 1, [nyp_l], [nyp_l], decomp_poisson%ysz(1), 1, &
                                    [nyp_l], decomp_poisson%ysz(1), 1, CUFFT_Z2Z, decomp_poisson%ysz(1) )
    End If
    If ( ierr /= 0 ) Stop 'ERROR: cuFFT pencil plan creation failed'

    If ( z_bc_type == 1 ) Then
       Allocate( duct_c(decomp_poisson%ysz(1), decomp_poisson%ysz(2), decomp_poisson%ysz(3)) )
       !$acc enter data create(duct_c)
       ! Qz/sqrt_w_z/lambda_z are built once on the host in initialization.f90 -- just upload them here
       Allocate( Qz_gpu(nzm_global,nzm_global), sqrt_w_z_gpu(nzm_global), lambda_z_gpu(nzm_global) )
       Qz_gpu       = Qz
       sqrt_w_z_gpu = sqrt_w_z
       lambda_z_gpu = lambda_z
       lambda_z_tol = 1d-8*Maxval(Abs(lambda_z))
       !$acc enter data create(Qz_gpu,sqrt_w_z_gpu,lambda_z_gpu)
       !$acc update device(Qz_gpu,sqrt_w_z_gpu,lambda_z_gpu)
    End If

    plans_created = .True.

  End Subroutine gpu_poisson_init

  !> Forward DCT-IV along x of every line of this rank's x-pencil: poisson_x_r -> poisson_x_c (real-valued), via a zero-padded 4N-point FFT + twiddle
  Subroutine gpu_dct_x_forward

    Integer :: i, j, k, ierr, nx1, n2, n3
    Real(Int64) :: theta
    Complex(Int64) :: w

    nx1 = Int(nxp_global)
    n2  = decomp_poisson%xsz(2); n3 = decomp_poisson%xsz(3)

    !$acc parallel loop collapse(3) present(poisson_x_r,dct_ext)
    Do k = 1, n3
       Do j = 1, n2
          Do i = 1, Lx
             If ( i <= nx1 ) Then
                dct_ext(i,j,k) = dcmplx( poisson_x_r(i,j,k) )
             Else
                dct_ext(i,j,k) = (0d0,0d0)
             End If
          End Do
       End Do
    End Do
    !$acc end parallel loop

    !$acc host_data use_device(dct_ext)
    ierr = cufftExecZ2Z( plan_dct_L, dct_ext, dct_ext, CUFFT_FORWARD )
    !$acc end host_data
    If ( ierr /= 0 ) Stop 'ERROR: cuFFT DCT-IV forward exec failed'

    ! extract DCT-IV coeffs via twiddle (S[k]=X4N(4N-2k-1), index Lx-2*imode)
    !$acc parallel loop collapse(3) present(dct_ext,poisson_x_c) private(theta,w)
    Do k = 1, n3
       Do j = 1, n2
          Do i = 0, nx1-1
             theta = pi*Real(2*i+1,Int64)/Real(4*nx1,Int64)
             w = dcmplx( dcos(theta), dsin(theta) )
             poisson_x_c(i+1,j,k) = 2d0*Real( w*dct_ext(Lx-2*i,j,k), Int64 )
          End Do
       End Do
    End Do
    !$acc end parallel loop

  End Subroutine gpu_dct_x_forward

  !> Inverse of gpu_dct_x_forward (DCT-IV is self-inverse up to 2N): Real(poisson_x_c) -> poisson_x_r, scaled by 1/(2*nxp*extra_norm) with extra_norm = nzp for the periodic-z FFT (1 otherwise)
  Subroutine gpu_dct_x_inverse(extra_norm)

    Real(Int64), Intent(In) :: extra_norm
    Integer :: i, j, k, ierr, nx1, n2, n3
    Real(Int64) :: theta, norm
    Complex(Int64) :: w

    nx1  = Int(nxp_global)
    n2   = decomp_poisson%xsz(2); n3 = decomp_poisson%xsz(3)
    norm = 2d0*Real(nxp_global,Int64)*extra_norm

    !$acc parallel loop collapse(3) present(poisson_x_c,dct_ext)
    Do k = 1, n3
       Do j = 1, n2
          Do i = 1, Lx
             If ( i <= nx1 ) Then
                dct_ext(i,j,k) = dcmplx( Real(poisson_x_c(i,j,k),Int64) )
             Else
                dct_ext(i,j,k) = (0d0,0d0)
             End If
          End Do
       End Do
    End Do
    !$acc end parallel loop

    !$acc host_data use_device(dct_ext)
    ierr = cufftExecZ2Z( plan_dct_L, dct_ext, dct_ext, CUFFT_FORWARD )
    !$acc end host_data
    If ( ierr /= 0 ) Stop 'ERROR: cuFFT DCT-IV inverse(self) exec failed'

    !$acc parallel loop collapse(3) present(dct_ext,poisson_x_r) private(theta,w)
    Do k = 1, n3
       Do j = 1, n2
          Do i = 0, nx1-1
             theta = pi*Real(2*i+1,Int64)/Real(4*nx1,Int64)
             w = dcmplx( dcos(theta), dsin(theta) )
             poisson_x_r(i+1,j,k) = 2d0*Real( w*dct_ext(Lx-2*i,j,k), Int64 ) / norm
          End Do
       End Do
    End Do
    !$acc end parallel loop

  End Subroutine gpu_dct_x_inverse

  !> Pencil-decomposed forward transform: y(real) -> x (cuFFT or DCT-IV) -> y -> z (cuFFT, periodic z only) -> y, leaving poisson_y_c in (kx,y,kz) space; mirrors projection.f90's CPU chain with device-resident arrays and CUDA-aware transposes
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

    If ( x_bc_type == 0 ) Then
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
    Else
       Call gpu_dct_x_forward
    End If
    Call transpose_x_to_y( poisson_x_c, poisson_y_c, decomp_poisson )

    ! z FFT for periodic z; the duct (z_bc_type==1) stays in physical z (eigenmode transform in the solve)
    If ( z_bc_type == 0 ) Then
       Call transpose_y_to_z( poisson_y_c, poisson_z_c, decomp_poisson )
       !$acc host_data use_device(poisson_z_c)
       ierr = cufftExecZ2Z( plan_p3_z, poisson_z_c, poisson_z_c, CUFFT_FORWARD )
       !$acc end host_data
       If ( ierr /= 0 ) Stop 'ERROR: cuFFT pencil z forward exec failed'
       Call transpose_z_to_y( poisson_z_c, poisson_y_c, decomp_poisson )
    End If

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

  !> Pencil-decomposed inverse transform: y -> z (inverse cuFFT, periodic z) -> y -> x (inverse cuFFT or DCT-IV) -> y(real) -> rhs_p, normalised here except for the y factor (applied in gpu_solve_periodic_3d)
  Subroutine gpu_inverse_transform_3d

    Integer :: i, j, k, ierr, n1, n2, n3
    Real(Int64) :: inv_norm, znorm

    If ( z_bc_type == 0 ) Then
       Call transpose_y_to_z( poisson_y_c, poisson_z_c, decomp_poisson )
       !$acc host_data use_device(poisson_z_c)
       ierr = cufftExecZ2Z( plan_p3_z, poisson_z_c, poisson_z_c, CUFFT_INVERSE )
       !$acc end host_data
       If ( ierr /= 0 ) Stop 'ERROR: cuFFT pencil z inverse exec failed'
       Call transpose_z_to_y( poisson_z_c, poisson_y_c, decomp_poisson )
    End If

    Call transpose_y_to_x( poisson_y_c, poisson_x_c, decomp_poisson )
    ! z solved via a real-space tridiagonal/eigenmode solve (duct) needs no z normalisation
    If ( z_bc_type == 0 ) Then
       znorm = Real(nzp_global,Int64)
    Else
       znorm = 1d0
    End If
    If ( x_bc_type == 0 ) Then
       !$acc host_data use_device(poisson_x_c)
       ierr = cufftExecZ2Z( plan_p3_x, poisson_x_c, poisson_x_c, CUFFT_INVERSE )
       !$acc end host_data
       If ( ierr /= 0 ) Stop 'ERROR: cuFFT pencil x inverse exec failed'

       inv_norm = 1d0 / ( Real(nxp_global,Int64) * znorm )
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
    Else
       Call gpu_dct_x_inverse(znorm)
    End If
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

  !> Create the cuSPARSE handle and batched tridiagonal-solve arrays once, on first use
  Subroutine gpu_gtsv_init

    Integer :: ierr
    Integer(Int64) :: bufsize

    gtsv_m     = nyg - 2
    ! one system per (kx,kz) mode (or (kx,z-eigenmode) for the duct) this rank owns in its y-pencil
    gtsv_batch = decomp_poisson%ysz(1) * decomp_poisson%ysz(3)

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

  !> Wall-normal (y) solve for the pencil path with periodic-FFT z and y walls (x periodic FFT or DCT-IV): batched cuSPARSE tridiagonal solve, one system per (kx,kz) mode this rank owns, operating in place on the y-pencil poisson_y_c; global mode indices come from decomp_poisson%yst
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
             If ( x_bc_type == 0 .And. ii == 0 .And. i_global == 0 .And. k_global == 0 ) gtsv_d(idx) = 3d0/2d0*gtsv_d(idx)
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

  !> 4-wall duct solve in the y-pencil (z fully local, p_col==1): z-eigenmode transform of poisson_y_c per local x-mode, batched cuSPARSE y-tridiagonal solve per (x-mode, z-eigenmode) pair with lambda_z(m) standing in for kzz, then the inverse eigenmode transform back into poisson_y_c (Qz orthonormal -> inverse == transpose; sqrt_w_z scaling recovers Dzz's own eigenbasis, see global.f90's lambda_z/Qz/sqrt_w_z comment)
  Subroutine gpu_solve_duct_pencil

    Integer :: i, jy, k, m, ii, j, idx, b, ierr, n1, ny_i, nz1, yst1, i_global
    Complex(Int64) :: acc

    If ( .Not. gtsv_created ) Call gpu_gtsv_init

    n1   = decomp_poisson%ysz(1)
    ny_i = decomp_poisson%ysz(2)
    nz1  = Int(nzm_global)
    yst1 = decomp_poisson%yst(1)

    !$acc data present(Dyy,kxx,lambda_z_gpu,Qz_gpu,sqrt_w_z_gpu,poisson_y_c,duct_c,gtsv_dl,gtsv_d,gtsv_du,gtsv_x,gtsv_buf)

    !$acc parallel loop collapse(3) private(acc)
    Do i = 1, n1
       Do jy = 1, ny_i
          Do m = 1, nz1
             acc = (0d0,0d0)
             Do k = 1, nz1
                acc = acc + poisson_y_c(i,jy,k) * sqrt_w_z_gpu(k) * Qz_gpu(k,m)
             End Do
             duct_c(i,jy,m) = acc
          End Do
       End Do
    End Do
    !$acc end parallel loop

    !$acc parallel loop collapse(2) private(i_global,b,ii,j,idx)
    Do m = 1, nz1
       Do i = 1, n1
          i_global = yst1 + i - 2
          b = (m-1)*n1 + (i-1)
          Do ii = 0, gtsv_m-1
             j   = ii + 2
             idx = ii*gtsv_batch + b + 1
             gtsv_d(idx) = Dyy(j,j) + kxx(i_global) + lambda_z_gpu(m)
             ! remove the null mode: zero x-wavenumber paired with z's zero eigenvalue (Dzz is
             ! negative-semi-definite, so found by magnitude rather than a hardcoded index)
             If ( x_bc_type == 0 .And. ii == 0 .And. i_global == 0 .And. Abs(lambda_z_gpu(m)) < lambda_z_tol ) &
                  gtsv_d(idx) = 3d0/2d0*gtsv_d(idx)
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
             gtsv_x(idx) = duct_c(i,ii+1,m)
          End Do
       End Do
    End Do
    !$acc end parallel loop

    !$acc host_data use_device(gtsv_dl,gtsv_d,gtsv_du,gtsv_x,gtsv_buf)
    ierr = cusparseZgtsvInterleavedBatch( cusparse_h, CUSPARSE_ALG1, gtsv_m, &
                 gtsv_dl, gtsv_d, gtsv_du, gtsv_x, gtsv_batch, gtsv_buf )
    !$acc end host_data
    If ( ierr /= 0 ) Stop 'ERROR: cusparseZgtsvInterleavedBatch (duct) failed'

    !$acc parallel loop collapse(2) private(b,ii,idx)
    Do m = 1, nz1
       Do i = 1, n1
          b = (m-1)*n1 + (i-1)
          Do ii = 0, gtsv_m-1
             idx = ii*gtsv_batch + b + 1
             duct_c(i,ii+1,m) = gtsv_x(idx)
          End Do
       End Do
    End Do
    !$acc end parallel loop

    !$acc parallel loop collapse(3) private(acc)
    Do i = 1, n1
       Do jy = 1, ny_i
          Do k = 1, nz1
             acc = (0d0,0d0)
             Do m = 1, nz1
                acc = acc + duct_c(i,jy,m) * Qz_gpu(k,m)
             End Do
             poisson_y_c(i,jy,k) = acc / sqrt_w_z_gpu(k)
          End Do
       End Do
    End Do
    !$acc end parallel loop

    !$acc end data

  End Subroutine gpu_solve_duct_pencil

End Module poisson_gpu
