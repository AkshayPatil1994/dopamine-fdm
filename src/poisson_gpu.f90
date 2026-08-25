!> Single-GPU cuFFT periodic Poisson transform (nprocs==1 only)
Module poisson_gpu

  Use iso_fortran_env, Only : Int32, Int64
  Use cufft
  Use cusparse
  Use global, Only : nxp_global, nzp_global, mx, mz, nyg, rhs_p, rhs_p_hat, Dyy, kxx, kzz, x_bc_type, pi

  Implicit None

  Integer :: plan_fwd, plan_bwd
  Logical :: plans_created = .False.
  Integer :: nslabs   ! number of interior y-planes batched per transform (nyg-2)

  Complex(Int64), Allocatable :: slab(:,:,:)

  ! DCT-IV (x_bc_type==1) state: zero-padded 4N-point Z2Z FFT + twiddle, then z-FFT
  Integer :: plan_dct_L, plan_z, Lx
  Complex(Int64), Allocatable :: dct_ext(:,:,:), zwork(:,:,:)

  ! Batched cuSPARSE tridiagonal solve state
  Type(cusparseHandle) :: cusparse_h
  Logical :: gtsv_created = .False.
  Integer :: gtsv_m, gtsv_batch
  Complex(Int64), Allocatable :: gtsv_dl(:), gtsv_d(:), gtsv_du(:), gtsv_x(:)
  Character(1), Allocatable :: gtsv_buf(:)

Contains

  !> Create the batched cuFFT plans once, on first use; branches on x_bc_type
  Subroutine gpu_poisson_init

    Integer :: ierr, nx1, nz1

    nx1    = Int(nxp_global)
    nz1    = Int(nzp_global)
    nslabs = nyg - 2

    If ( x_bc_type == 0 ) Then

       Allocate( slab( nx1, nz1, nslabs ) )
       !$acc enter data create(slab)

       ! cufftPlanMany dims use the reversed-relative-to-Fortran convention of cufftPlan2D
       ierr = cufftPlanMany( plan_fwd, 2, [nz1, nx1], [nz1, nx1], 1, nx1*nz1, &
                                          [nz1, nx1], 1, nx1*nz1, CUFFT_Z2Z, nslabs )
       ierr = ierr + cufftPlanMany( plan_bwd, 2, [nz1, nx1], [nz1, nx1], 1, nx1*nz1, &
                                                [nz1, nx1], 1, nx1*nz1, CUFFT_Z2Z, nslabs )
       If ( ierr /= 0 ) Stop 'ERROR: cuFFT batched plan creation failed'

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

  !> Forward-transform every y-slab of rhs_p into rhs_p_hat in one batched cuFFT call; direct kx=i, kz=k indexing (no MPI transpose, unlike the FFTW-MPI path). Fully device-resident: rhs_p/rhs_p_hat/slab never touch host here.
  Subroutine gpu_forward_transform_all_slabs

    Integer :: j, ix, iz, ierr, nx1, nz1

    If ( .Not. plans_created ) Call gpu_poisson_init

    nx1 = Int(nxp_global)
    nz1 = Int(nzp_global)

    !$acc parallel loop collapse(3) present(rhs_p,slab)
    Do j = 2, nyg-1
       Do iz = 1, nz1
          Do ix = 1, nx1
             slab(ix,iz,j-1) = dcmplx( rhs_p(ix+1,j,iz+1) )
          End Do
       End Do
    End Do
    !$acc end parallel loop

    !$acc host_data use_device(slab)
    ierr = cufftExecZ2Z( plan_fwd, slab, slab, CUFFT_FORWARD )
    !$acc end host_data
    If ( ierr /= 0 ) Stop 'ERROR: cuFFT batched forward exec failed'

    !$acc parallel loop collapse(3) present(slab,rhs_p_hat)
    Do j = 2, nyg-1
       Do iz = 1, nz1
          Do ix = 1, nx1
             rhs_p_hat(j,ix-1,iz-1) = slab(ix,iz,j-1)
          End Do
       End Do
    End Do
    !$acc end parallel loop

  End Subroutine gpu_forward_transform_all_slabs

  !> Inverse-transform every y-slab of rhs_p_hat back into rhs_p in one batched cuFFT call, normalised by nxp_global*nzp_global. Fully device-resident.
  Subroutine gpu_inverse_transform_all_slabs

    Integer :: j, ix, iz, ierr, nx1, nz1
    Real(Int64) :: norm

    nx1  = Int(nxp_global)
    nz1  = Int(nzp_global)
    norm = Real(nxp_global,Int64) * Real(nzp_global,Int64)

    !$acc parallel loop collapse(3) present(rhs_p_hat,slab)
    Do j = 2, nyg-1
       Do iz = 1, nz1
          Do ix = 1, nx1
             slab(ix,iz,j-1) = rhs_p_hat(j,ix-1,iz-1)
          End Do
       End Do
    End Do
    !$acc end parallel loop

    !$acc host_data use_device(slab)
    ierr = cufftExecZ2Z( plan_bwd, slab, slab, CUFFT_INVERSE )
    !$acc end host_data
    If ( ierr /= 0 ) Stop 'ERROR: cuFFT batched inverse exec failed'

    !$acc parallel loop collapse(3) present(slab,rhs_p)
    Do j = 2, nyg-1
       Do iz = 1, nz1
          Do ix = 1, nx1
             rhs_p(ix+1,j,iz+1) = Real(slab(ix,iz,j-1),Int64) / norm
          End Do
       End Do
    End Do
    !$acc end parallel loop

  End Subroutine gpu_inverse_transform_all_slabs

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
    gtsv_batch = ( Int(mx,Int32) + 1 ) * ( Int(mz,Int32) + 1 )

    Allocate( gtsv_dl(gtsv_m*gtsv_batch), gtsv_d(gtsv_m*gtsv_batch), &
              gtsv_du(gtsv_m*gtsv_batch), gtsv_x(gtsv_m*gtsv_batch) )
    ! placeholder values for the sizing query below (bufferSize depends only on m/batchCount)
    gtsv_dl = (0d0,0d0)
    gtsv_d  = (1d0,0d0)
    gtsv_du = (0d0,0d0)
    gtsv_x  = (0d0,0d0)

    ierr = cusparseCreate( cusparse_h )
    If ( ierr /= 0 ) Stop 'ERROR: cusparseCreate failed'

    !$acc data copy(gtsv_dl,gtsv_d,gtsv_du,gtsv_x)
    !$acc host_data use_device(gtsv_dl,gtsv_d,gtsv_du,gtsv_x)
    ierr = cusparseZgtsvInterleavedBatch_bufferSizeExt( cusparse_h, CUSPARSE_ALG1, gtsv_m, &
                 gtsv_dl, gtsv_d, gtsv_du, gtsv_x, gtsv_batch, bufsize )
    !$acc end host_data
    !$acc end data
    If ( ierr /= 0 ) Stop 'ERROR: cusparseZgtsvInterleavedBatch_bufferSizeExt failed'

    Allocate( gtsv_buf(bufsize) )

    gtsv_created = .True.

  End Subroutine gpu_gtsv_init

  !> Pack rhs_p_hat, batch-solve all y-tridiagonal systems via cuSPARSE, unpack back
  Subroutine gpu_solve_tridiagonal_batched

    Integer :: imode, k, b, ii, j, idx, ierr, mx_i, mz_i

    If ( .Not. gtsv_created ) Call gpu_gtsv_init

    mx_i = Int(mx,Int32)
    mz_i = Int(mz,Int32)

    ! rhs_p_hat is persistently device-resident (initialization.f90) -- present()
    ! fails loudly instead of silently re-transferring if that assumption ever breaks
    !$acc data copyin(Dyy,kxx,kzz) present(rhs_p_hat) &
    !$acc      create(gtsv_dl,gtsv_d,gtsv_du,gtsv_x,gtsv_buf)

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

End Module poisson_gpu
