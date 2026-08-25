! Module for fractional step method
Module projection

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use mpi
  Use profiler
#ifdef GPU_POISSON
  Use poisson_gpu
#endif

  ! prevent implicit typing
  Implicit None

Contains

  ! Compute incompressible velocity with fractional step method
  Subroutine compute_projection_step

    Call profiler_start(PROF_PROJECTION)
    Call compute_pseudo_pressure_rhs
    Call profiler_stop(PROF_PROJECTION)

    Call solve_poisson_equation

    Call profiler_start(PROF_PROJECTION)
    Call project_velocity
    Call profiler_stop(PROF_PROJECTION)

  End Subroutine compute_projection_step

  !> Compute right-hand side for pseudo-pressure equation Laplacian(p) = rhs_p = div(vel*)
  Subroutine compute_pseudo_pressure_rhs

    Integer(Int32) :: i, j, k
    Real   (Int64) :: inv_dx_p, inv_dz_p, inv_dyj

    inv_dx_p = 1d0 / dx
    inv_dz_p = 1d0 / dz

    ! rhs_p located at cell centers
    !$acc parallel loop collapse(2) present(rhs_p,U,V,W,y)
    Do k = 2, nzg-1
       Do j = 2, nyg-1
          inv_dyj = 1d0 / ( y(j) - y(j-1) )   ! varies with j
          Do i = 2, nxg-1
             rhs_p(i,j,k) = ( U(i,j,k) - U(i-1,j,k) ) * inv_dx_p + &
                             ( V(i,j,k) - V(i,j-1,k) ) * inv_dyj  + &
                             ( W(i,j,k) - W(i,j,k-1) ) * inv_dz_p
          End Do
       End Do
    End Do
    !$acc end parallel loop

  End Subroutine compute_pseudo_pressure_rhs

  !> Solve pseudo-pressure equation (fast Poisson solver)
  Subroutine solve_poisson_equation

    Integer(Int32) :: i, j, k, k_global, i_global, info
    Real   (Int64) :: dum, dumref, maxerr, wavenum_sum

    Call profiler_start(PROF_POISSON_FFT)
#ifdef GPU_POISSON
    ! rhs_p/rhs_p_hat stay device-resident through the whole Poisson stage; single-GPU cuFFT transform (nprocs==1 only)
    If ( x_bc_type == 0 ) Then
       Call gpu_forward_transform_all_slabs
    Else
       Call gpu_forward_transform_dct_slabs
    End If
#else
    ! Forward transform interior points
    Do j = 2, nyg-1
       If ( x_bc_type == 0 ) Then
          plane = dcmplx( rhs_p ( 2:nxp+1, j, 2:nzp+1 ) )
       Else
          xplane = rhs_p ( 2:nxp+1, j, 2:nzp+1 )
          Call fftw_execute_r2r(plan_dct, xplane, xplane)
          plane = dcmplx( xplane )
       End If
       Call fftw_mpi_execute_dft(plan_d, plane, plane_hat)
       rhs_p_hat (j, :, :) = plane_hat
    End Do
#endif
    Call profiler_stop(PROF_POISSON_FFT)

    Call profiler_start(PROF_POISSON_TRIDIAG)
#ifdef GPU_POISSON
    ! Batched cuSPARSE solve of all (mx+1)*(mz+1) y-tridiagonal systems in one call
    Call gpu_solve_tridiagonal_batched
#else
    ! Zgtsv overwrites DL, D, DU in-place, so re-initialise all three inside the loop before every call
    ! solve for each mode
    Do k = 0, mz
       Do i = 0, mx
          If ( x_bc_type == 0 ) Then
             ! Periodic mode-index mapping
             i_global = imode_map_fft( i, k )
             k_global = kmode_map_fft( i, k )
          Else
             ! Inflow/outflow mode-index mapping
             i_global = i
             k_global = Int(local_k_offset,Int32) + k
          End If
          ! Re-fill DL and DU (Zgtsv destroys them on exit)
          Do j = 2, nyg-2
             DL(j) = Dyy(j+1,j)   ! lower diagonal
             DU(j) = Dyy(j,j+1)   ! upper diagonal
          End Do
          ! kxx+kzz is constant for this (i,k) mode
          wavenum_sum = kxx(i_global) + kzz(k_global)
          ! form diagonal
          Do j = 2, nyg-1
             D(j) = Dyy(j,j) + wavenum_sum
          End Do
          ! Remove singularity 00 mode
          If ( x_bc_type==0 .And. i_global==0 .And. k_global==0 ) D(2) = 3d0/2d0*D(2)
          ! solve M*u = rhs; rhs_p_hat(2:nyg-1,i,k) is contiguous (y-first
          ! layout) so Zgtsv writes the solution directly in-place.
          Call Zgtsv( nr, nrhs, DL, D, DU, rhs_p_hat(2:nyg-1,i,k), nr, info)
       End Do
    End Do
#endif
    Call profiler_stop(PROF_POISSON_TRIDIAG)

    Call profiler_start(PROF_POISSON_FFT)
#ifdef GPU_POISSON
    If ( x_bc_type == 0 ) Then
       Call gpu_inverse_transform_all_slabs
    Else
       Call gpu_inverse_transform_dct_slabs
    End If
#else
    ! Inverse transform (mirror of the forward pass above)
    Do j = 2, nyg-1
       plane_hat = rhs_p_hat (j, :, :)
       Call fftw_mpi_execute_dft(plan_i, plane_hat, plane)
       If ( x_bc_type == 0 ) Then
          rhs_p ( 2:nxp+1, j, 2:nzp+1 ) = plane/Real( nxp_global*nzp_global, 8)
       Else
          xplane = Real( plane, 8 )
          Call fftw_execute_r2r(plan_dct, xplane, xplane)
          rhs_p ( 2:nxp+1, j, 2:nzp+1 ) = xplane/Real( 2*nxp_global*nzp_global, 8)
       End If
    End Do
#endif

    ! boundary conditions for pressure in x and z
    Call apply_periodic_xz_pressure

    ! update ghost interior planes
    Call update_ghost_interior_planes_pressure
    Call profiler_stop(PROF_POISSON_FFT)

    ! Save physical pressure at the end of the full time step
    If ( rk_step == 3 ) Then
#ifdef GPU_POISSON
       ! rhs_p is device-resident throughout; P is host-only, so pull rhs_p back
       ! just for this end-of-step save -- only 1 in 3 RK substages pays this transfer
       !$acc update host(rhs_p)
#endif
       P( 2:nxg-1, 2:nyg-1, 2:nzg-1 ) = rhs_p( 2:nxg-1, 2:nyg-1, 2:nzg-1 )/(dt*rk_coef(3,3))
       P(  1,:,:) = P(    2,:,:)
       P(nxg,:,:) = P(nxg-1,:,:)
       P(:,  1,:) = P(:,    2,:)
       P(:,nyg,:) = P(:,nyg-1,:)
       P(:,:,  1) = P(:,:,    2)
       P(:,:,nzg) = P(:,:,nzg-1)
    End If

  End Subroutine solve_poisson_equation

  ! Boundary conditions for pseudo-pressure in x and z
  Subroutine apply_periodic_xz_pressure

    ! x periodic ghost-fill (inflow/outflow needs none)
    ! GPU-resident: rhs_p never leaves the device for this fast path
    !$acc kernels present(rhs_p)
    If ( x_bc_type == 0 ) rhs_p ( nxg-1, :, : ) = rhs_p ( 2, :, : )
    !$acc end kernels

    ! apply periodicity in z (only first and last processor, MPI needed)
    ! nprocs==1: both ranks are the same, so copy locally to avoid a blocking self-send deadlock
    If ( nprocs == 1 ) Then
       !$acc kernels present(rhs_p)
       rhs_p ( 2:nxg-1, :, nzp+1+1 ) = rhs_p ( 2:nxg-1, :, 2 )
       !$acc end kernels
    Elseif ( myid==0 ) Then
       buffer_p = rhs_p ( 2:nxg-1, :, 2 )
       Call Mpi_send(buffer_p, (nxg-2)*(nyg-2), MPI_real8, nprocs-1, 0, &
            MPI_COMM_WORLD,ierr)
    Elseif ( myid==nprocs-1 ) Then
       Call Mpi_recv(buffer_p, (nxg-2)*(nyg-2), MPI_real8, 0, 0, &
            MPI_COMM_WORLD,istat,ierr)
       rhs_p ( 2:nxg-1, :, nzp+1+1 ) = buffer_p
    End If

  End Subroutine apply_periodic_xz_pressure

  ! Update ghost interior planes for pressure
  Subroutine update_ghost_interior_planes_pressure

    Integer(Int32) :: sendto, recvfrom
    Integer(Int32) :: tagto,  tagfrom
    
    ! update P
    ! send to bottom processor, receive from top one
    sendto   = myid - 1
    tagto    = myid - 1
    recvfrom = myid + 1
    tagfrom  = myid 
    If ( myid==0 ) Then
       sendto = MPI_PROC_NULL
       tagto  = 0
    End If
    If ( myid==nprocs-1 ) Then
       recvfrom = MPI_PROC_NULL
       tagfrom  = MPI_ANY_TAG
    End If
    buffer_ps = rhs_p(:,:,2)  ! send buffer
    Call Mpi_sendrecv(buffer_ps, (nxg-2)*(nyg-2), Mpi_real8, sendto, tagto,        &
         buffer_pr, (nxg-2)*(nyg-2), Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
         istat, ierr)   
    If ( myid/=nprocs-1 ) rhs_p(:,:,nzg) = buffer_pr ! received buffer 

  End Subroutine update_ghost_interior_planes_pressure

  !> Project velocity field: vel = vel* - gradient(p)
  Subroutine project_velocity

    Integer(Int32) :: i, j, k

    Real(Int64) :: inv_dx_proj, inv_dz_proj, inv_dygj
    inv_dx_proj = 1d0 / dx
    inv_dz_proj = 1d0 / dz

    !$acc kernels present(U,V,W,rhs_p,yg)
    ! project U (interior points)
    Do i = 2, nx-1
       U(i,2:nyg-1,2:nzg-1) = U(i,2:nyg-1,2:nzg-1) - &
            ( rhs_p(i+1,2:nyg-1,2:nzg-1) - rhs_p(i,2:nyg-1,2:nzg-1) ) * inv_dx_proj
    End Do

    ! Outflow face correction (x_bc_type==1)
    If ( x_bc_type == 1 ) Then
       U(nx,2:nyg-1,2:nzg-1) = U(nx,2:nyg-1,2:nzg-1) + &
            2d0 * rhs_p(nx,2:nyg-1,2:nzg-1) * inv_dx_proj
    End If

    ! project V (interior points); yg spacing varies with j
    Do j = 2, ny-1
       inv_dygj = 1d0 / ( yg(j+1) - yg(j) )
       V(2:nxg-1,j,2:nzg-1) = V(2:nxg-1,j,2:nzg-1) - &
            ( rhs_p(2:nxg-1,j+1,2:nzg-1) - rhs_p(2:nxg-1,j,2:nzg-1) ) * inv_dygj
    End Do

    ! project W (interior points)
    Do k = 2, nz-1
       W(2:nxg-1,2:nyg-1,k) = W(2:nxg-1,2:nyg-1,k) - &
            ( rhs_p(2:nxg-1,2:nyg-1,k+1) - rhs_p(2:nxg-1,2:nyg-1,k) ) * inv_dz_proj
    End Do
    !$acc end kernels

  End Subroutine project_velocity

  !> Compute maximum divergence for interior points
  Subroutine check_divergence(max_divergence)

    Real(Int64), Intent(Out) :: max_divergence

    Real   (Int64) :: max_divergence_local, div
    Integer(Int32) :: i, j, k    

    div = 0d0
    Do k = 2, nzg-1     
       Do j = 2, nyg-1
          Do i = 2, nxg-1
             If ( ibm_input_mode >= 1 ) Then
                ! Skip solid cells: ghost-cell face velocities make divergence inside the solid physically meaningless
                If ( phi(i,j,k) < 0d0 ) Cycle
                ! Skip IBM-adjacent fluid cells
                If ( phi(i-1,j,k) < 0d0 .Or. phi(i+1,j,k) < 0d0 .Or. &
                     phi(i,j-1,k) < 0d0 .Or. phi(i,j+1,k) < 0d0 .Or. &
                     phi(i,j,k-1) < 0d0 .Or. phi(i,j,k+1) < 0d0 ) Cycle
             End If
             div  = Max(div, Abs(( U(i,j,k) - U(i-1,j,k) )/( x(i)-x(i-1) ) + &
                                  ( V(i,j,k) - V(i,j-1,k) )/( y(j)-y(j-1) ) + &
                                  ( W(i,j,k) - W(i,j,k-1) )/( z(k)-z(k-1) )) )
          End Do
       End Do
    End Do

    max_divergence_local = div

    Call MPI_Reduce(max_divergence_local,max_divergence,1,MPI_real8,MPI_max,0,MPI_COMM_WORLD,ierr)

  End Subroutine check_divergence

End module projection
