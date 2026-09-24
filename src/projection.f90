! Module for fractional step method
Module projection

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use mpi
  Use decomp, Only : decomp_poisson, transpose_x_to_y, transpose_y_to_x, transpose_y_to_z, transpose_z_to_y, &
                      z_halo_neighbors, z_periodic_partner, x_halo_neighbors, x_periodic_partner
  Use boundary_conditions, Only : update_ghost_interior_planes, update_ghost_interior_planes_x, &
                      apply_periodic_bc_x, apply_periodic_bc_z
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
    Real   (Int64) :: inv_dx_p, inv_dzk, inv_dyj

    inv_dx_p = 1d0 / dx

    ! rhs_p located at cell centers; z spacing z(k)-z(k-1) varies with k when z_bc_type==1
    ! and alpha_grid_z>0 (spanwise stretching) -- constant dz only when z is uniform (periodic,
    ! or an unstretched wall), in which case z(k)-z(k-1) reduces to dz exactly.
    !$acc parallel loop collapse(2) present(rhs_p,U,V,W,y,z)
    Do k = 2, nzg-1
       Do j = 2, nyg-1
          ! computed per-j (not hoisted above the j loop) since collapse(2) requires
          ! the two collapsed Do statements to be immediately nested
          inv_dzk = 1d0 / ( z(k) - z(k-1) )   ! varies with k
          inv_dyj = 1d0 / ( y(j) - y(j-1) )   ! varies with j
          Do i = 2, nxg-1
             rhs_p(i,j,k) = ( U(i,j,k) - U(i-1,j,k) ) * inv_dx_p + &
                             ( V(i,j,k) - V(i,j-1,k) ) * inv_dyj  + &
                             ( W(i,j,k) - W(i,j,k-1) ) * inv_dzk
          End Do
       End Do
    End Do
    !$acc end parallel loop

  End Subroutine compute_pseudo_pressure_rhs

  !> Solve pseudo-pressure equation (fast Poisson solver)
  Subroutine solve_poisson_equation

    Integer(Int32) :: i, j, k, k_global, i_global, j_global, info, nyp
    Real   (Int64) :: dum, dumref, maxerr, wavenum_sum, inv_pdt
    Logical        :: is_first_p, is_last_p
    Integer(Int32) :: partner_p

    Call profiler_start(PROF_POISSON_FFT)
#ifdef GPU_POISSON
    ! rhs_p/rhs_p_hat stay device-resident through the whole Poisson stage; single-GPU cuFFT transform (nprocs==1 only)
    ! z_bc_type==1 (4-wall duct, y_bc_type==1 .And. x_bc_type==0 only -- other z-wall
    ! combinations are rejected for GPU_POISSON builds at init time, initialization.f90)
    ! pencil-decomposed transform chain (x FFT or DCT-IV, z FFT unless duct), multi-GPU capable for every BC combination
    Call gpu_forward_transform_3d
#else
    ! Forward chain, common part: y(real) -> x(local FFT/DCT) -> y(complex)
    poisson_y_r = rhs_p ( 2:decomp_poisson%ysz(1)+1, 2:decomp_poisson%ysz(2)+1, 2:decomp_poisson%ysz(3)+1 )
    Call transpose_y_to_x( poisson_y_r, poisson_x_r, decomp_poisson )
    If ( x_bc_type == 0 ) Then
       poisson_x_c = dcmplx( poisson_x_r )
       Call fftw_execute_dft(plan_fx_fwd, poisson_x_c, poisson_x_c)
    Else
       Call fftw_execute_r2r(plan_dct, poisson_x_r, poisson_x_r)
       poisson_x_c = dcmplx( poisson_x_r )
    End If
    Call transpose_x_to_y( poisson_x_c, poisson_y_c, decomp_poisson )

    If ( z_bc_type == 0 ) Then
       ! ---- default: z transformed first (FFT, always periodic today) -> back to y-pencil ----
       Call transpose_y_to_z( poisson_y_c, poisson_z_c, decomp_poisson )
       Call fftw_execute_dft(plan_fz_fwd, poisson_z_c, poisson_z_c)
       Call transpose_z_to_y( poisson_z_c, poisson_y_c, decomp_poisson )
    Else If ( y_bc_type == 0 ) Then
       ! ---- spanwise wall only (z_bc_type==1, y periodic): FFT y instead, while still fully
       ! local in the y-pencil, exactly mirroring the periodic-y FFT below but leaving the
       ! result in y-wavenumber space (no diagonal divide here -- z still needs to be solved,
       ! via Zgtsv/Dzz in the tridiagonal section below).
       nyp = decomp_poisson%ysz(2) - 1
       Do k = 1, decomp_poisson%ysz(3)
          Call fftw_execute_dft( plan_fy_fwd, poisson_y_c(:,1:nyp,k), poisson_y_c(:,1:nyp,k) )
       End Do
    Else
       ! ---- 4-wall duct (y_bc_type==1 .And. z_bc_type==1): neither direction separates via
       ! FFT, so nothing to do here -- stay in physical (y,z) space. p_col==1 is enforced for
       ! this combination (see decomp_auto_factorize / read_input_parameters), so the y-pencil
       ! is already fully local in z too; the z-eigenmode transform + per-mode y-tridiagonal
       ! solve happens entirely in the tridiagonal section below.
    End If
#endif
    Call profiler_stop(PROF_POISSON_FFT)

    Call profiler_start(PROF_POISSON_TRIDIAG)
#ifdef GPU_POISSON
    If ( z_bc_type == 1 ) Then
       ! 4-wall duct: z-eigenmode transform + batched cuSPARSE solve of one y-tridiagonal system per (x-mode, z-eigenmode) pair
       Call gpu_solve_duct_pencil
    Else If ( x_bc_type == 0 .And. y_bc_type == 0 ) Then
       ! Fully periodic: elementwise divide by kxx+kyy+kzz in Fourier space (no tridiagonal solve needed)
       Call gpu_solve_periodic_3d
    Else
       ! y walls: per-rank batched cuSPARSE tridiagonal solve over the local (kx,kz) modes
       Call gpu_solve_tridiagonal_pencil
    End If
#else
    If ( z_bc_type == 1 .And. y_bc_type == 1 ) Then
       ! ---- 4-wall duct: decouple z via its eigenbasis (Qz/lambda_z, built once at init from
       ! Dzz's interior tridiagonal part -- see initialization.f90), then solve one y-tridiagonal
       ! system per z-eigenmode, exactly like the periodic-z/wall-y branch below but with
       ! lambda_z(m) standing in for kzz(k_global). p_col==1 is enforced for this BC combination,
       ! so poisson_y_c is already fully local in both y and z here -- no z-pencil transpose needed.
       Block
         Integer(Int32) :: m
         ! Qz_c/sqrt_w_z_c/z_hat_duct are built once at init (initialization.f90) instead of
         ! Allocated/Deallocated on every call here (3x per step, one per RK substage)
         Do i = 1, decomp_poisson%ysz(1)
            i_global = decomp_poisson%yst(1) + i - 2

            ! forward z-eigenmode transform: scale by sqrt(cell width) first (recovers Dzz's own,
            ! generally non-orthogonal eigenvectors from S's orthonormal ones -- see global.f90's
            ! lambda_z/Qz/sqrt_w_z comment), then project onto S's eigenbasis
            Do k = 1, nzm_global
               z_hat_duct(:,k) = poisson_y_c(i,:,k) * sqrt_w_z_c(k)
            End Do
            z_hat_duct = Matmul( z_hat_duct, Qz_c )

            Do m = 1, nzm_global
               Do j = 2, nyg-2
                  DL(j) = Dyy(j+1,j)
                  DU(j) = Dyy(j,j+1)
               End Do
               wavenum_sum = kxx(i_global) + lambda_z(m)
               Do j = 2, nyg-1
                  D(j) = Dyy(j,j) + wavenum_sum
               End Do
               ! Remove singularity of the 00 mode: x_bc_type==0's zero x-wavenumber paired with
               ! z's zero eigenvalue (the constant-function null mode of the Neumann-pressure
               ! z-operator). Checked by magnitude, not a hardcoded index -- dstev returns
               ! eigenvalues ascending, and Dzz is negative-semi-definite, so the zero eigenvalue
               ! is the LARGEST (last, not first) of the ascending list.
               If ( x_bc_type==0 .And. i_global==0 .And. Abs(lambda_z(m)) < 1d-8*Maxval(Abs(lambda_z)) ) D(2) = 3d0/2d0*D(2)
               Call Zgtsv( nr, nrhs, DL, D, DU, z_hat_duct(:,m), nr, info)
            End Do

            ! inverse z-eigenmode transform: unproject (Qz orthonormal -> inverse == transpose),
            ! then unscale by sqrt(cell width) to recover Dzz's own eigenbasis
            poisson_y_c(i,:,:) = Matmul( z_hat_duct, Transpose(Qz_c) )
            Do k = 1, nzm_global
               poisson_y_c(i,:,k) = poisson_y_c(i,:,k) / sqrt_w_z_c(k)
            End Do
         End Do
       End Block
    ElseIf ( z_bc_type == 1 ) Then
       ! ---- spanwise wall only (y periodic): Zgtsv-in-z, one mode per (x,y)-wavenumber pair this rank owns, in the
       ! z-pencil (fully local in z there -- reuses the same transpose the z-FFT path uses). The
       ! z-pencil splits (x,y), unlike the y-pencil above which splits (x,z), so the y-mode's global
       ! index needs decomp_poisson%zst(2) (mirroring i_global's use of yst(1)/zst(1) elsewhere);
       ! poisson_y_c's untouched padding column (index nyp+1, beyond the nyp real y-DOF) rides along
       ! through the transpose and is simply skipped here (j_global >= nyp).
       Call transpose_y_to_z( poisson_y_c, poisson_z_c, decomp_poisson )
       Do j = 1, decomp_poisson%zsz(2)
          j_global = decomp_poisson%zst(2) + j - 2
          If ( j_global < 0 .Or. j_global >= nyp ) Cycle
          Do i = 1, decomp_poisson%zsz(1)
             i_global = decomp_poisson%zst(1) + i - 2
             ! Re-fill DL and DU (Zgtsv destroys them on exit)
             Do k = 2, nzg_global-2
                DL(k) = Dzz(k+1,k)   ! lower diagonal
                DU(k) = Dzz(k,k+1)   ! upper diagonal
             End Do
             ! kxx+kyy is constant for this (i,j) mode
             wavenum_sum = kxx(i_global) + kyy(j_global)
             Do k = 2, nzg_global-1
                D(k) = Dzz(k,k) + wavenum_sum
             End Do
             ! Remove singularity of the 00 mode (mirrors the y-wall branch's own null-space fix below)
             If ( x_bc_type==0 .And. i_global==0 .And. j_global==0 ) D(2) = 3d0/2d0*D(2)
             Call Zgtsv( nr, nrhs, DL, D, DU, poisson_z_c(i,j,:), nr, info)
          End Do
       End Do
       Call transpose_z_to_y( poisson_z_c, poisson_y_c, decomp_poisson )
    ElseIf ( y_bc_type == 0 ) Then
       ! Periodic y: local batched 1-D FFT along dim2 of poisson_y_c (already
       ! fully local in y for this rank), one k-slice at a time (fftw_plan_many_dft
       ! only batches over one stride/dist pair, here the stride-1 x-index),
       ! then a diagonal divide by kxx+kyy+kzz per mode, then inverse FFT --
       ! replaces the Zgtsv tridiagonal solve entirely for this branch.
       ! nyp = decomp_poisson%ysz(2)-1: same "last cell is a redundant duplicate
       ! of the first" convention as periodic x/z (see apply_periodic_bc_y and
       ! the reconstruction copy below) -- only the first nyp of ysz(2) y-slots
       ! are independent DOF for the FFT.
       nyp = decomp_poisson%ysz(2) - 1
       Do k = 1, decomp_poisson%ysz(3)
          Call fftw_execute_dft( plan_fy_fwd, poisson_y_c(:,1:nyp,k), poisson_y_c(:,1:nyp,k) )
       End Do
       Do k = 1, decomp_poisson%ysz(3)
          k_global = decomp_poisson%yst(3) + k - 2
          Do i = 1, decomp_poisson%ysz(1)
             i_global = decomp_poisson%yst(1) + i - 2
             Do j = 1, nyp
                wavenum_sum = kxx(i_global) + kyy(j-1) + kzz(k_global)
                If ( i_global==0 .And. j==1 .And. k_global==0 ) Then
                   poisson_y_c(i,j,k) = (0d0,0d0)   ! periodic Poisson null-space fix: reference pressure = 0
                Else
                   poisson_y_c(i,j,k) = poisson_y_c(i,j,k) / wavenum_sum
                End If
             End Do
          End Do
       End Do
       Do k = 1, decomp_poisson%ysz(3)
          Call fftw_execute_dft( plan_fy_inv, poisson_y_c(:,1:nyp,k), poisson_y_c(:,1:nyp,k) )
       End Do
       poisson_y_c(:,1:nyp,:) = poisson_y_c(:,1:nyp,:) / Real(nyp, Int64)
    Else
       ! Zgtsv overwrites DL, D, DU in-place, so re-initialise all three inside the loop before every call
       ! solve for each mode this rank owns (local index -> global mode index via decomp_poisson%yst)
       Do k = 1, decomp_poisson%ysz(3)
          k_global = decomp_poisson%yst(3) + k - 2
          Do i = 1, decomp_poisson%ysz(1)
             i_global = decomp_poisson%yst(1) + i - 2
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
             ! solve M*u = rhs; poisson_y_c(i,:,k) is a length-nym vector like D/DL/DU (Zgtsv only cares about sequence position, not the declared bounds), written in-place
             Call Zgtsv( nr, nrhs, DL, D, DU, poisson_y_c(i,:,k), nr, info)
          End Do
       End Do
    End If
#endif
    Call profiler_stop(PROF_POISSON_TRIDIAG)

    Call profiler_start(PROF_POISSON_FFT)
#ifdef GPU_POISSON
    Call gpu_inverse_transform_3d
#else
    If ( z_bc_type == 0 ) Then
       ! ---- default: z -> y (inverse local FFT) ----
       Call transpose_y_to_z( poisson_y_c, poisson_z_c, decomp_poisson )
       Call fftw_execute_dft(plan_fz_inv, poisson_z_c, poisson_z_c)
       Call transpose_z_to_y( poisson_z_c, poisson_y_c, decomp_poisson )
    Else If ( y_bc_type == 0 ) Then
       ! ---- spanwise wall only: inverse y-FFT (z was already inverted in physical space by the Zgtsv above) ----
       Do k = 1, decomp_poisson%ysz(3)
          Call fftw_execute_dft( plan_fy_inv, poisson_y_c(:,1:nyp,k), poisson_y_c(:,1:nyp,k) )
       End Do
       poisson_y_c(:,1:nyp,:) = poisson_y_c(:,1:nyp,:) / Real(nyp, Int64)
    Else
       ! ---- 4-wall duct: already back in physical (y,z) space -- the tridiagonal section's
       ! inverse z-eigenmode transform (Qz orthonormal) handled z, and the y-tridiagonal solve
       ! is a direct physical-space solve, no y transform to invert either.
    End If
    ! Inverse chain, common tail: y -> x(inverse local FFT/DCT) -> y(real) -> rhs_p
    Call transpose_y_to_x( poisson_y_c, poisson_x_c, decomp_poisson )
    If ( x_bc_type == 0 ) Then
       Call fftw_execute_dft(plan_fx_inv, poisson_x_c, poisson_x_c)
       If ( z_bc_type == 0 ) Then
          poisson_x_r = Real( poisson_x_c, 8 ) / Real( nxp_global*nzp_global, 8)
       Else
          ! z solved via a real-space tridiagonal (Zgtsv), not an FFT -- no z normalisation factor needed
          poisson_x_r = Real( poisson_x_c, 8 ) / Real( nxp_global, 8)
       End If
    Else
       poisson_x_r = Real( poisson_x_c, 8 )
       Call fftw_execute_r2r(plan_dct, poisson_x_r, poisson_x_r)
       If ( z_bc_type == 0 ) Then
          poisson_x_r = poisson_x_r / Real( 2*nxp_global*nzp_global, 8)
       Else
          poisson_x_r = poisson_x_r / Real( 2*nxp_global, 8)
       End If
    End If
    Call transpose_x_to_y( poisson_x_r, poisson_y_r, decomp_poisson )
    rhs_p ( 2:decomp_poisson%ysz(1)+1, 2:decomp_poisson%ysz(2)+1, 2:decomp_poisson%ysz(3)+1 ) = poisson_y_r
#endif

    ! boundary conditions for pressure in x and z
    Call apply_periodic_xz_pressure

    ! periodic y: reconstruct the redundant last interior cell (nyg-1), excluded
    ! from the FFT above, exactly mirroring apply_periodic_xz_pressure's x-periodic fixup
    If ( y_bc_type == 0 ) Then
       !$acc kernels present(rhs_p)
       rhs_p ( :, nyg-1, : ) = rhs_p ( :, 2, : )
       !$acc end kernels
    End If

    ! update ghost interior planes
    Call update_ghost_interior_planes_pressure
    Call update_ghost_interior_planes_pressure_x
    Call profiler_stop(PROF_POISSON_FFT)

    ! Save physical pressure at the end of the full time step
    If ( rk_step == 3 ) Then
#ifdef GPU_POISSON
       ! P is device-resident like rhs_p: computed and ghost-filled entirely on the device (no host round
       ! trips); the host copy is refreshed only when a host consumer needs it (see time_integration.f90)
       inv_pdt = 1d0/(dt*rk_coef(3,3))
       !$acc parallel loop collapse(3) present(P,rhs_p)
       Do k = 2, nzg-1
          Do j = 2, nyg-1
             Do i = 2, nxg-1
                P(i,j,k) = rhs_p(i,j,k)*inv_pdt
             End Do
          End Do
       End Do
       !$acc end parallel loop
       !$acc kernels present(P)
       P(:,  1,:) = P(:,    2,:)
       P(:,nyg,:) = P(:,nyg-1,:)
       !$acc end kernels
#else
       P( 2:nxg-1, 2:nyg-1, 2:nzg-1 ) = rhs_p( 2:nxg-1, 2:nyg-1, 2:nzg-1 )/(dt*rk_coef(3,3))
       P(:,  1,:) = P(:,    2,:)
       P(:,nyg,:) = P(:,nyg-1,:)
#endif
       ! x/z ghosts: cross-rank halo + domain-periodic wrap (id=4, generic scalar-shaped convention), not
       ! a same-rank self-copy -- a self-copy here duplicated the neighbour rank's own edge plane in the
       ! output snapshot at every interior MPI boundary (visible as spurious repeated planes along z/x)
       Call update_ghost_interior_planes_x(P,4)
       Call update_ghost_interior_planes(P,4)
       If ( x_bc_type == 0 ) Then
          Call apply_periodic_bc_x(P,4)
       Else
          Call x_periodic_partner(is_first_p, is_last_p, partner_p)
#ifdef GPU_POISSON
          !$acc kernels present(P)
          If ( is_first_p ) P(1,:,:) = P(2,:,:)
          If ( is_last_p  ) P(nxg,:,:) = P(nxg-1,:,:)
          !$acc end kernels
#else
          If ( is_first_p ) P(1,:,:) = P(2,:,:)
          If ( is_last_p  ) P(nxg,:,:) = P(nxg-1,:,:)
#endif
       End If
       If ( z_bc_type == 0 ) Then
          Call apply_periodic_bc_z(P,4)
       Else
          Call z_periodic_partner(is_first_p, is_last_p, partner_p)
#ifdef GPU_POISSON
          !$acc kernels present(P)
          If ( is_first_p ) P(:,:,1)   = P(:,:,2)
          If ( is_last_p  ) P(:,:,nzg) = P(:,:,nzg-1)
          !$acc end kernels
#else
          If ( is_first_p ) P(:,:,1)   = P(:,:,2)
          If ( is_last_p  ) P(:,:,nzg) = P(:,:,nzg-1)
#endif
       End If
    End If

  End Subroutine solve_poisson_equation

  ! Boundary conditions for pseudo-pressure in x and z
  Subroutine apply_periodic_xz_pressure

    Logical        :: is_first, is_last
    Integer(Int32) :: partner
    Logical        :: is_first_x, is_last_x
    Integer(Int32) :: partner_x

    ! x periodic ghost-fill (inflow/outflow needs none): only the row owning the
    ! domain's x tail is missing this point (nxp_global's periodic reduction), and
    ! only the row owning x=1 has the source value -- MPI needed unless p_row==1
    If ( x_bc_type == 0 ) Then
       Call x_periodic_partner(is_first_x, is_last_x, partner_x)
       If ( is_first_x .And. is_last_x ) Then
          !$acc kernels present(rhs_p)
          rhs_p ( nxg-1, :, : ) = rhs_p ( 2, :, : )
          !$acc end kernels
       Elseif ( is_first_x ) Then
#ifdef GPU_POISSON
          !$acc update host(rhs_p(2,:,:))
#endif
          buffer_px = rhs_p ( 2, :, : )
          Call Mpi_send(buffer_px, (nyg-2)*(nzg-1), MPI_real8, partner_x, 0, MPI_COMM_WORLD, ierr)
       Elseif ( is_last_x ) Then
          Call Mpi_recv(buffer_px, (nyg-2)*(nzg-1), MPI_real8, partner_x, 0, MPI_COMM_WORLD, istat, ierr)
          rhs_p ( nxg-1, :, : ) = buffer_px
#ifdef GPU_POISSON
          !$acc update device(rhs_p(nxg-1,:,:))
#endif
       End If
    End If

    ! apply periodicity in z (only the column's first and last rank, MPI needed); no
    ! equivalent needed for z_bc_type==1 (wall) -- like y_bc_type==1, project_velocity
    ! never reads rhs_p's outermost z ghost plane in that case (W is exactly zero at
    ! the wall face itself, so no pressure-gradient correction is needed there)
    If ( z_bc_type == 0 ) Then
       Call z_periodic_partner(is_first, is_last, partner)
       ! p_col==1: this rank self-pairs, so copy locally to avoid a blocking self-send deadlock
       If ( is_first .And. is_last ) Then
          !$acc kernels present(rhs_p)
          rhs_p ( 2:nxg-1, :, nzp+1+1 ) = rhs_p ( 2:nxg-1, :, 2 )
          !$acc end kernels
       Elseif ( is_first ) Then
#ifdef GPU_POISSON
          !$acc update host(rhs_p(2:nxg-1,:,2))
#endif
          buffer_p = rhs_p ( 2:nxg-1, :, 2 )
          Call Mpi_send(buffer_p, (nxg-2)*(nyg-2), MPI_real8, partner, 0, &
               MPI_COMM_WORLD,ierr)
       Elseif ( is_last ) Then
          Call Mpi_recv(buffer_p, (nxg-2)*(nyg-2), MPI_real8, partner, 0, &
               MPI_COMM_WORLD,istat,ierr)
          rhs_p ( 2:nxg-1, :, nzp+1+1 ) = buffer_p
#ifdef GPU_POISSON
          !$acc update device(rhs_p(2:nxg-1,:,nzp+1+1))
#endif
       End If
    End If

  End Subroutine apply_periodic_xz_pressure

  ! Update ghost interior planes for pressure
  Subroutine update_ghost_interior_planes_pressure

    Integer(Int32) :: up, down

    Call z_halo_neighbors(up, down)

    ! update P: send towards +z, receive from -z
#ifdef GPU_POISSON
    !$acc update host(rhs_p(2:nxg-1,:,2))
#endif
    buffer_ps = rhs_p(2:nxg-1,:,2)  ! send buffer
    Call Mpi_sendrecv(buffer_ps, (nxg-2)*(nyg-2), Mpi_real8, down, 0,             &
         buffer_pr, (nxg-2)*(nyg-2), Mpi_real8, up,   0, MPI_COMM_WORLD, &
         istat, ierr)
    If ( up /= MPI_PROC_NULL ) Then
       rhs_p(2:nxg-1,:,nzg) = buffer_pr ! received buffer
#ifdef GPU_POISSON
       !$acc update device(rhs_p(2:nxg-1,:,nzg))
#endif
    End If

  End Subroutine update_ghost_interior_planes_pressure

  !> Fill rhs_p's x row-seam ghost (nxg) from the next row's first interior point; needed by project_velocity's U update at the row-to-row boundary (distinct from the domain periodic wrap, handled in apply_periodic_xz_pressure)
  Subroutine update_ghost_interior_planes_pressure_x

    Integer(Int32) :: up, down

    Call x_halo_neighbors(up, down)

#ifdef GPU_POISSON
    !$acc update host(rhs_p(2,:,:))
#endif
    buffer_pgxs = rhs_p(2,:,:)
    Call Mpi_sendrecv(buffer_pgxs, (nyg-2)*(nzg-1), Mpi_real8, down, 0,             &
         buffer_pgxr, (nyg-2)*(nzg-1), Mpi_real8, up,   0, MPI_COMM_WORLD, istat, ierr)
    If ( up /= MPI_PROC_NULL ) Then
       rhs_p(nxg,:,:) = buffer_pgxr
#ifdef GPU_POISSON
       !$acc update device(rhs_p(nxg,:,:))
#endif
    End If

  End Subroutine update_ghost_interior_planes_pressure_x

  !> Project velocity field: vel = vel* - gradient(p)
  Subroutine project_velocity

    Integer(Int32) :: i, j, k
    Logical        :: is_first_x, is_last_x
    Integer(Int32) :: partner_x

    Real(Int64) :: inv_dx_proj, inv_dzgk, inv_dygj
    inv_dx_proj = 1d0 / dx

    Call x_periodic_partner(is_first_x, is_last_x, partner_x)

    !$acc kernels present(U,V,W,rhs_p,yg,zg)
    ! project U (interior points)
    Do i = 2, nx-1
       U(i,2:nyg-1,2:nzg-1) = U(i,2:nyg-1,2:nzg-1) - &
            ( rhs_p(i+1,2:nyg-1,2:nzg-1) - rhs_p(i,2:nyg-1,2:nzg-1) ) * inv_dx_proj
    End Do

    ! Row-to-row seam face (x_bc_type-independent): only the last row's face nx is
    ! the domain boundary (periodic wrap / outflow, handled separately below); every
    ! other row's face nx is a genuine interior face against the next row, needing
    ! the x-halo-filled ghost at nxg (update_ghost_interior_planes_pressure_x)
    If ( .Not. is_last_x ) Then
       U(nx,2:nyg-1,2:nzg-1) = U(nx,2:nyg-1,2:nzg-1) - &
            ( rhs_p(nxg,2:nyg-1,2:nzg-1) - rhs_p(nx,2:nyg-1,2:nzg-1) ) * inv_dx_proj
    End If

    ! Outflow face correction (x_bc_type==1, domain boundary only)
    If ( x_bc_type == 1 .And. is_last_x ) Then
       U(nx,2:nyg-1,2:nzg-1) = U(nx,2:nyg-1,2:nzg-1) + &
            2d0 * rhs_p(nx,2:nyg-1,2:nzg-1) * inv_dx_proj
    End If

    ! project V (interior points); yg spacing varies with j
    Do j = 2, ny-1
       inv_dygj = 1d0 / ( yg(j+1) - yg(j) )
       V(2:nxg-1,j,2:nzg-1) = V(2:nxg-1,j,2:nzg-1) - &
            ( rhs_p(2:nxg-1,j+1,2:nzg-1) - rhs_p(2:nxg-1,j,2:nzg-1) ) * inv_dygj
    End Do

    ! project W (interior points); zg spacing varies with k (z_bc_type==1, alpha_grid_z>0)
    Do k = 2, nz-1
       inv_dzgk = 1d0 / ( zg(k+1) - zg(k) )
       W(2:nxg-1,2:nyg-1,k) = W(2:nxg-1,2:nyg-1,k) - &
            ( rhs_p(2:nxg-1,2:nyg-1,k+1) - rhs_p(2:nxg-1,2:nyg-1,k) ) * inv_dzgk
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
