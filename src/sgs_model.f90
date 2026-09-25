!> SGS eddy-viscosity models for LES (sgs_model: 0=DNS, 1=Vreman)
Module sgs_models

  Use iso_fortran_env, Only : Int32, Int64
  Use global
  Use mpi
  Use decomp, Only : z_halo_neighbors, z_periodic_partner, x_periodic_partner
  Use boundary_conditions, Only : apply_periodic_bc_z, update_ghost_interior_planes_x, update_ghost_interior_planes
#ifdef GPU_POISSON
  Use boundary_conditions, Only : gpu_periodic_wrap
#endif

  Implicit None

  ! Module-level halo buffers for nu_t z-exchange (avoids per-call heap allocation)
  Real(Int64), Allocatable, Dimension(:,:,:) :: sgs_snd_lo, sgs_snd_hi
  Real(Int64), Allocatable, Dimension(:,:,:) :: sgs_rcv_lo, sgs_rcv_hi

Contains

  !> Dispatch to the selected SGS model (sgs_model==0: DNS no-op, ==1: Vreman)
  Subroutine compute_sgs_model(U_, V_, W_, nu_t_)

    Real(Int64), Dimension(nx,  nyg, nzg), Intent(In)    :: U_
    Real(Int64), Dimension(nxg, ny,  nzg), Intent(In)    :: V_
    Real(Int64), Dimension(nxg, nyg, nz ), Intent(In)    :: W_
    Real(Int64), Dimension(nxg, nyg, nzg), Intent(InOut) :: nu_t_

    Select Case (sgs_model)
    Case (0)
       ! DNS: nu_t stays zero (set during initialisation, never modified here)
       Return
    Case (1)
       Call compute_vreman(U_, V_, W_, nu_t_)
    Case Default
       If (myid == 0) Write(*,*) 'WARNING: unknown sgs_model =', sgs_model, ' — using DNS (nu_t=0).'
       Return
    End Select

  End Subroutine compute_sgs_model


  !> Vreman (2004) SGS model: nu_t=c_V*sqrt(B_beta/(alpha_ij alpha_ij))
  Subroutine compute_vreman(U_, V_, W_, nu_t_)

    Real(Int64), Dimension(nx,  nyg, nzg), Intent(In)    :: U_
    Real(Int64), Dimension(nxg, ny,  nzg), Intent(In)    :: V_
    Real(Int64), Dimension(nxg, nyg, nz ), Intent(In)    :: W_
    Real(Int64), Dimension(nxg, nyg, nzg), Intent(InOut) :: nu_t_

    ! local scalars
    Integer(Int32) :: i, j, k
    Logical        :: is_first, is_last
    Integer(Int32) :: partner
    Real   (Int64) :: xbuf2(nyg,nzg), xbuf3(nyg,nzg,2)
    Real   (Int64) :: dx_c, dy_c, dz_c          ! local filter widths (IBM pass only)
    Real   (Int64) :: dx2, dz2                  ! dx^2 (uniform) and dz_c^2 (per k, z may be stretched)
    Real   (Int64) :: inv_dx                    ! 1/dx
    Real   (Int64) :: inv_zg_km, inv_zg_kp      ! per-k z-grid reciprocals
    Real   (Int64) :: inv_dz_k                  ! 1/(z(k)-z(k-1)) — per k
    Real   (Int64) :: inv_yg_jm, inv_yg_jp      ! per-j y-grid reciprocals
    Real   (Int64) :: inv_dy_j                  ! 1/(y(j)-y(j-1)) — per j
    Real   (Int64) :: dy2_c                     ! dy_c^2 — per j (standard path)
    Real   (Int64) :: a11, a12, a13             ! alpha_1j = du_j/dx_1 (at cell ctr)
    Real   (Int64) :: a21, a22, a23             ! alpha_2j = du_j/dx_2
    Real   (Int64) :: a31, a32, a33             ! alpha_3j = du_j/dx_3
    Real   (Int64) :: b11, b12, b13, b22, b23, b33  ! beta_mn
    Real   (Int64) :: B_beta, alpha_sq, c_V, nu_t_loc
    Logical        :: ibm_active               ! Guard Umask_cc accesses

    c_V = 2.5d0 * Cs_vreman**2
    ibm_active = ibm_input_mode >= 1 .And. Allocated(Umask_cc)

    inv_dx = 1d0 / dx
    dx2    = dx * dx

    ! Pass 1: standard (non-IBM) stencil, unconditional so it stays GPU-offloadable without phi residency; ibm_active discards it via Pass 2/solid-cell zeroing below
    !$acc parallel loop collapse(2) present(U_,V_,W_,nu_t_,y,yg,z,zg)
    Do k = 2, nzg-1
       Do j = 2, nyg-1

          ! Per-(k,j) constants: y- and z-grid spacings (non-uniform y, and z when z_bc_type==1 with alpha_grid_z>0; x uniform)
          dz_c      = z(k) - z(k-1)
          dz2       = dz_c * dz_c
          inv_dz_k  = 1d0 / dz_c
          inv_zg_km = 1d0 / ( zg(k)   - zg(k-1) )
          inv_zg_kp = 1d0 / ( zg(k+1) - zg(k  ) )
          dy_c      = y(j) - y(j-1)
          dy2_c     = dy_c * dy_c
          inv_yg_jm = 1d0 / ( yg(j)   - yg(j-1) )
          inv_yg_jp = 1d0 / ( yg(j+1) - yg(j  ) )
          inv_dy_j  = 1d0 / dy_c

          Do i = 2, nxg-1

             ! ------ velocity gradient tensor (standard stencils) -------
             a11 = ( U_(i,j,k)   - U_(i-1,j,k)   ) * inv_dx
             ! Non-uniform y: use per-j precomputed inverses
             a21 = 0.5d0*( ( U_(i,j,k)   - U_(i,j-1,k)  ) * inv_yg_jm + &
                           ( U_(i,j+1,k) - U_(i,j,  k)  ) * inv_yg_jp )
             a31 = 0.5d0*( ( U_(i,j,k)   - U_(i,j,k-1)  ) * inv_zg_km + &
                           ( U_(i,j,k+1) - U_(i,j,k  )  ) * inv_zg_kp )
             a12 = 0.5d0*( ( V_(i,j,k)   - V_(i-1,j,k)  ) + &
                           ( V_(i+1,j,k) - V_(i,  j,k)  ) ) * inv_dx
             ! Non-uniform y: V on y-faces
             a22 = ( V_(i,j,k)   - V_(i,j-1,k)   ) * inv_dy_j
             a32 = 0.5d0*( ( V_(i,j,k)   - V_(i,j,k-1)  ) * inv_zg_km + &
                           ( V_(i,j,k+1) - V_(i,j,k  )  ) * inv_zg_kp )
             a13 = 0.5d0*( ( W_(i,j,k)   - W_(i-1,j,k)  ) + &
                           ( W_(i+1,j,k) - W_(i,  j,k)  ) ) * inv_dx
             ! Non-uniform y
             a23 = 0.5d0*( ( W_(i,j,k)   - W_(i,j-1,k)  ) * inv_yg_jm + &
                           ( W_(i,j+1,k) - W_(i,j,  k)  ) * inv_yg_jp )
             ! W on z-faces
             a33 = ( W_(i,j,k)   - W_(i,j,k-1)   ) * inv_dz_k

             ! ------ beta_mn (standard: dx_c=dx, dz_c=z(k)-z(k-1), dy_c=y(j)-y(j-1)) --
             b11 = dx2*(a11*a11) + dy2_c*(a21*a21) + dz2*(a31*a31)
             b22 = dx2*(a12*a12) + dy2_c*(a22*a22) + dz2*(a32*a32)
             b33 = dx2*(a13*a13) + dy2_c*(a23*a23) + dz2*(a33*a33)
             b12 = dx2*(a11*a12) + dy2_c*(a21*a22) + dz2*(a31*a32)
             b13 = dx2*(a11*a13) + dy2_c*(a21*a23) + dz2*(a31*a33)
             b23 = dx2*(a12*a13) + dy2_c*(a22*a23) + dz2*(a32*a33)

             ! ------ B_beta: second invariant of beta -------------------
             B_beta = b11*b22 - b12**2 &
                    + b11*b33 - b13**2 &
                    + b22*b33 - b23**2

             ! ------ alpha_ij alpha_ij = sum of squared gradients -------
             alpha_sq = a11**2 + a12**2 + a13**2 &
                      + a21**2 + a22**2 + a23**2 &
                      + a31**2 + a32**2 + a33**2

             ! ------ eddy viscosity -------------------------------------
             If ( alpha_sq > 1d-20 .And. B_beta > 0d0 ) Then
                nu_t_loc = c_V * Sqrt( B_beta / alpha_sq )
             Else
                nu_t_loc = 0d0
             End If

             nu_t_(i,j,k) = nu_t_loc

          End Do
       End Do
    End Do
    !$acc end parallel loop

#ifdef GPU_POISSON
    ! Without IBM there is no host-only Pass 2, so wall zeroing, ghost fills and MPI halos all stay on the device
    ! (previously nu_t made three full-array host<->device round trips per RK substage)
    If ( .Not. ibm_active ) Then
       Call vreman_ghosts_device(nu_t_)
       ! host copy only needed by host consumers (snapshot write syncs it itself, see input_output.f90)
       If ( rsb_active == 1 ) Then
          !$acc update host(nu_t_)
       End If
       Return
    End If
#endif

    ! Pass 2, wall zeroing, x-periodicity, and MPI halo exchange below all run on host; sync Pass 1's GPU result back first
    !$acc update host(nu_t_)

    ! Pass 2: IBM corrections (only when ibm_active); re-visits all fluid cells in the domain (O(volume), not O(surface)), does not affect Pass 1
    If ( ibm_active ) Then
       Do k = 2, nzg-1
          Do j = 2, nyg-1
             dz_c      = z(k) - z(k-1)
             inv_dz_k  = 1d0 / dz_c
             inv_zg_km = 1d0 / ( zg(k)   - zg(k-1) )
             inv_zg_kp = 1d0 / ( zg(k+1) - zg(k  ) )
             dy_c      = y(j) - y(j-1)
             inv_yg_jm = 1d0 / ( yg(j)   - yg(j-1) )
             inv_yg_jp = 1d0 / ( yg(j+1) - yg(j  ) )
             inv_dy_j  = 1d0 / dy_c

             Do i = 2, nxg-1
                If ( phi(i,j,k) <= 0d0 ) Cycle   ! solid cell: zeroed in post-processing

                ! Re-evaluate gradients; then apply one-sided corrections where needed
                a11 = ( U_(i,j,k)   - U_(i-1,j,k)   ) * inv_dx
                a21 = 0.5d0*( (U_(i,j,k)-U_(i,j-1,k))*inv_yg_jm + (U_(i,j+1,k)-U_(i,j,k))*inv_yg_jp )
                a31 = 0.5d0*( (U_(i,j,k)-U_(i,j,k-1))*inv_zg_km + (U_(i,j,k+1)-U_(i,j,k))*inv_zg_kp )
                a12 = 0.5d0*( (V_(i,j,k)-V_(i-1,j,k)) + (V_(i+1,j,k)-V_(i,j,k)) ) * inv_dx
                a22 = ( V_(i,j,k)   - V_(i,j-1,k)   ) * inv_dy_j
                a32 = 0.5d0*( (V_(i,j,k)-V_(i,j,k-1))*inv_zg_km + (V_(i,j,k+1)-V_(i,j,k))*inv_zg_kp )
                a13 = 0.5d0*( (W_(i,j,k)-W_(i-1,j,k)) + (W_(i+1,j,k)-W_(i,j,k)) ) * inv_dx
                a23 = 0.5d0*( (W_(i,j,k)-W_(i,j-1,k))*inv_yg_jm + (W_(i,j+1,k)-W_(i,j,k))*inv_yg_jp )
                a33 = ( W_(i,j,k)   - W_(i,j,k-1)   ) * inv_dz_k

                ! One-sided stencil corrections near IBM solid faces
                If ( Umask_cc(i,j-1,k) < 0.5d0 .And. Umask_cc(i,j+1,k) < 0.5d0 ) Then
                   a21 = 0d0;  a23 = 0d0
                Else If ( Umask_cc(i,j-1,k) < 0.5d0 ) Then
                   a21 = ( U_(i,j+1,k) - U_(i,j,k) ) * inv_yg_jp
                   a23 = ( W_(i,j+1,k) - W_(i,j,k) ) * inv_yg_jp
                Else If ( Umask_cc(i,j+1,k) < 0.5d0 ) Then
                   a21 = ( U_(i,j,k)   - U_(i,j-1,k) ) * inv_yg_jm
                   a23 = ( W_(i,j,k)   - W_(i,j-1,k) ) * inv_yg_jm
                End If

                If ( Umask_cc(i,j,k-1) < 0.5d0 .And. Umask_cc(i,j,k+1) < 0.5d0 ) Then
                   a31 = 0d0;  a32 = 0d0
                Else If ( Umask_cc(i,j,k-1) < 0.5d0 ) Then
                   a31 = ( U_(i,j,k+1) - U_(i,j,k) ) * inv_zg_kp
                   a32 = ( V_(i,j,k+1) - V_(i,j,k) ) * inv_zg_kp
                Else If ( Umask_cc(i,j,k+1) < 0.5d0 ) Then
                   a31 = ( U_(i,j,k)   - U_(i,j,k-1) ) * inv_zg_km
                   a32 = ( V_(i,j,k)   - V_(i,j,k-1) ) * inv_zg_km
                End If

                If ( Umask_cc(i-1,j,k) < 0.5d0 .And. Umask_cc(i+1,j,k) < 0.5d0 ) Then
                   a12 = 0d0;  a13 = 0d0
                Else If ( Umask_cc(i-1,j,k) < 0.5d0 ) Then
                   a12 = ( V_(i+1,j,k) - V_(i,j,k) ) * inv_dx
                   a13 = ( W_(i+1,j,k) - W_(i,j,k) ) * inv_dx
                Else If ( Umask_cc(i+1,j,k) < 0.5d0 ) Then
                   a12 = ( V_(i,j,k)   - V_(i-1,j,k) ) * inv_dx
                   a13 = ( W_(i,j,k)   - W_(i-1,j,k) ) * inv_dx
                End If

                ! IBM filter-width clamping: collapse filter to zero at the surface
                dx_c = Min(dx, 2d0*phi(i,j,k))
                dy_c = Min(y(j)-y(j-1), 2d0*phi(i,j,k))
                dz_c = Min(z(k)-z(k-1), 2d0*phi(i,j,k))

                b11 = dx_c*dx_c*(a11*a11) + dy_c*dy_c*(a21*a21) + dz_c*dz_c*(a31*a31)
                b22 = dx_c*dx_c*(a12*a12) + dy_c*dy_c*(a22*a22) + dz_c*dz_c*(a32*a32)
                b33 = dx_c*dx_c*(a13*a13) + dy_c*dy_c*(a23*a23) + dz_c*dz_c*(a33*a33)
                b12 = dx_c*dx_c*(a11*a12) + dy_c*dy_c*(a21*a22) + dz_c*dz_c*(a31*a32)
                b13 = dx_c*dx_c*(a11*a13) + dy_c*dy_c*(a21*a23) + dz_c*dz_c*(a31*a33)
                b23 = dx_c*dx_c*(a12*a13) + dy_c*dy_c*(a22*a23) + dz_c*dz_c*(a32*a33)

                B_beta = b11*b22 - b12**2 &
                       + b11*b33 - b13**2 &
                       + b22*b33 - b23**2

                alpha_sq = a11**2 + a12**2 + a13**2 &
                         + a21**2 + a22**2 + a23**2 &
                         + a31**2 + a32**2 + a33**2

                If ( alpha_sq > 1d-20 .And. B_beta > 0d0 ) Then
                   nu_t_(i,j,k) = c_V * Sqrt( B_beta / alpha_sq )
                Else
                   nu_t_(i,j,k) = 0d0
                End If

             End Do
          End Do
       End Do
    End If

    ! No-slip flat walls only (j=1, nyg): zero nu_t to avoid polluting the Robin BC; leave nu_t alone at a free-slip boundary (bc_face_y*==2), which has no molecular sublayer to damp it
    If ( y_bc_type == 1 .And. bc_face_ylo == 1 ) nu_t_(:,  1,:) = 0d0
    If ( y_bc_type == 1 .And. bc_face_yhi == 1 ) nu_t_(:,nyg,:) = 0d0

    ! y-periodicity (y_bc_type==0): fill ghost planes j=1,nyg (never written by
    ! Pass 1/2 above, but read by compute_rhs_v/w at the y boundaries)
    If ( y_bc_type == 0 ) Then
       nu_t_(:,    1,:) = nu_t_(:,nyg-2,:)
       nu_t_(:,nyg-1,:) = nu_t_(:,    2,:)
       nu_t_(:,nyg  ,:) = nu_t_(:,    3,:)
    End If

    ! IBM solid cells: suppress SGS stress to avoid polluting adjacent fluid
    If ( ibm_input_mode >= 1 .And. Allocated(Umask_cc) ) Then
       Where ( Umask_cc(2:nxg-1, 2:nyg-1, 2:nzg-1) < 0.5d0 )
          nu_t_(2:nxg-1, 2:nyg-1, 2:nzg-1) = 0d0
       End Where
    End If

    ! x-periodicity: fill ghost planes i=1,nxg (never written above, but read by compute_rhs_v/w at the x boundaries)
    ! (x-split: seam planes come from the x-neighbour, the wrap from the partner rank at the opposite domain edge)
#ifdef GPU_POISSON
    ! this (IBM, host-only Pass 2) path works on the host copy, but the x seam exchange now runs on device memory:
    ! push the just-corrected host nu_t to the device, exchange there, and bring the filled ghost planes back
    !$acc update device(nu_t_)
#endif
    Call update_ghost_interior_planes_x(nu_t_, 2)
#ifdef GPU_POISSON
    !$acc update host(nu_t_)
#endif
    Call x_periodic_partner(is_first, is_last, partner)
    If ( x_bc_type == 0 ) Then
       ! periodic: ghost 1 <- cell nxg-2, cells nxg-1,nxg <- 2,3 (cell nxg-1 duplicates cell 2), as apply_periodic_bc_x does for cell-centred fields
       If ( is_first .And. is_last ) Then
          nu_t_(1,    :,:) = nu_t_(nxg-2,:,:)
          nu_t_(nxg-1,:,:) = nu_t_(2,    :,:)
          nu_t_(nxg,  :,:) = nu_t_(3,    :,:)
       Else If ( is_first ) Then
          xbuf3(:,:,1) = nu_t_(2,:,:);  xbuf3(:,:,2) = nu_t_(3,:,:)
          Call Mpi_sendrecv(xbuf3, 2*nyg*nzg, Mpi_real8, partner, 13, xbuf2, nyg*nzg, Mpi_real8, partner, 14, &
                            MPI_COMM_WORLD, istat, ierr)
          nu_t_(1,:,:) = xbuf2
       Else If ( is_last ) Then
          xbuf2 = nu_t_(nxg-2,:,:)
          Call Mpi_sendrecv(xbuf2, nyg*nzg, Mpi_real8, partner, 14, xbuf3, 2*nyg*nzg, Mpi_real8, partner, 13, &
                            MPI_COMM_WORLD, istat, ierr)
          nu_t_(nxg-1,:,:) = xbuf3(:,:,1);  nu_t_(nxg,:,:) = xbuf3(:,:,2)
       End If
    Else
       ! inflow/outflow: zero-gradient at the domain edges (interior-seam planes were filled by the exchange above)
       If ( is_first ) nu_t_(1,  :,:) = nu_t_(2,    :,:)
       If ( is_last  ) nu_t_(nxg,:,:) = nu_t_(nxg-1,:,:)
    End If
    ! Ring exchange for intermediate ranks (host-only); rank-0/rank-(nprocs-1) wrap handled below.
    Call update_ghost_interior_planes_nut(nu_t_)

    ! Push host state to device before the z-boundary fill below, which runs device-resident
    ! at nprocs==1 and would otherwise have its result clobbered by a later blanket update device
    !$acc update device(nu_t_)
    If ( z_bc_type == 0 ) Then
       Call apply_periodic_bc_z(nu_t_, 4)
    Else
       ! No-slip/wall-model z walls: zero nu_t at the wall faces (mirrors the y-wall
       ! zeroing above) instead of periodically wrapping from the opposite domain
       ! edge -- z is domain-decomposed, so only the rank(s) owning the z=0/z=Lz
       ! physical boundary act (as in apply_Dirichlet_bc_z/apply_Robin_bc_z).
       Call z_periodic_partner(is_first, is_last, partner)
       !$acc kernels present(nu_t_)
       If ( is_first ) nu_t_(:,:,1)   = 0d0
       If ( is_last  ) nu_t_(:,:,nzg) = 0d0
       !$acc end kernels
    End If
    ! Sync back to host for output_data's snapshot write and compute_wall_model's host-only bits
    !$acc update host(nu_t_)

  End Subroutine compute_vreman


#ifdef GPU_POISSON
  !> Device-resident equivalent of compute_vreman's post-Pass-1 host block (non-IBM only): y wall/periodic ghosts, x seam + periodic wrap, z seam + z BC, all on nu_t_ in device memory, MPI on device buffers
  Subroutine vreman_ghosts_device(nu_t_)

    Real(Int64), Dimension(nxg, nyg, nzg), Intent(InOut) :: nu_t_
    Logical        :: is_first, is_last
    Integer(Int32) :: partner

    ! No-slip flat walls only (j=1, nyg): zero nu_t; leave it alone at a free-slip boundary
    If ( y_bc_type == 1 .And. bc_face_ylo == 1 ) Then
       !$acc kernels present(nu_t_)
       nu_t_(:,  1,:) = 0d0
       !$acc end kernels
    End If
    If ( y_bc_type == 1 .And. bc_face_yhi == 1 ) Then
       !$acc kernels present(nu_t_)
       nu_t_(:,nyg,:) = 0d0
       !$acc end kernels
    End If
    If ( y_bc_type == 0 ) Then
       !$acc kernels present(nu_t_)
       nu_t_(:,    1,:) = nu_t_(:,nyg-2,:)
       nu_t_(:,nyg-1,:) = nu_t_(:,    2,:)
       nu_t_(:,nyg  ,:) = nu_t_(:,    3,:)
       !$acc end kernels
    End If

    ! x seam planes from the x-neighbour, then the periodic wrap (single-plane, cell-centred)
    Call update_ghost_interior_planes_x(nu_t_, 2)
    Call x_periodic_partner(is_first, is_last, partner)
    If ( x_bc_type == 0 ) Then
       If ( is_first .And. is_last ) Then
          !$acc kernels present(nu_t_)
          nu_t_(1,    :,:) = nu_t_(nxg-2,:,:)
          nu_t_(nxg-1,:,:) = nu_t_(2,    :,:)
          nu_t_(nxg,  :,:) = nu_t_(3,    :,:)
          !$acc end kernels
       Else
          Call gpu_periodic_wrap(nu_t_, .True., 2, nxg, is_first, is_last, partner, 13)
       End If
    Else
       !$acc kernels present(nu_t_)
       If ( is_first ) nu_t_(1,  :,:) = nu_t_(2,    :,:)
       If ( is_last  ) nu_t_(nxg,:,:) = nu_t_(nxg-1,:,:)
       !$acc end kernels
    End If

    ! z seam planes, then z periodic wrap or wall zeroing
    Call update_ghost_interior_planes(nu_t_, 2)
    If ( z_bc_type == 0 ) Then
       Call apply_periodic_bc_z(nu_t_, 4)
    Else
       Call z_periodic_partner(is_first, is_last, partner)
       !$acc kernels present(nu_t_)
       If ( is_first ) nu_t_(:,:,1)   = 0d0
       If ( is_last  ) nu_t_(:,:,nzg) = 0d0
       !$acc end kernels
    End If

  End Subroutine vreman_ghosts_device
#endif

  !> MPI ring exchange for nu_t z-halo (cell-centred, 1 ghost/side)
  Subroutine update_ghost_interior_planes_nut(F_)

    Real(Int64), Dimension(nxg, nyg, nzg), Intent(InOut) :: F_

    Integer(Int32) :: up, down

    ! Allocate persistent single-plane buffers on first call
    If (.Not. Allocated(sgs_snd_lo)) Then
       Allocate( sgs_snd_lo(nxg,nyg,1), sgs_snd_hi(nxg,nyg,1) )
       Allocate( sgs_rcv_lo(nxg,nyg,1), sgs_rcv_hi(nxg,nyg,1) )
    End If

    Call z_halo_neighbors(up, down)

    !-- Exchange A: send nzg-1 towards +z; receive z=1 from -z -------
    sgs_snd_lo(:,:,1) = F_(:,:,nzg-1)
    Call Mpi_Sendrecv( sgs_snd_lo, nxg*nyg, Mpi_real8, up,   0, &
                       sgs_rcv_lo, nxg*nyg, Mpi_real8, down, 0, &
                       MPI_COMM_WORLD, istat, ierr )
    If (down /= MPI_PROC_NULL) F_(:,:,1) = sgs_rcv_lo(:,:,1)

    !-- Exchange B: send z=2 towards -z; receive z=nzg from +z -----
    sgs_snd_hi(:,:,1) = F_(:,:,2)
    Call Mpi_Sendrecv( sgs_snd_hi, nxg*nyg, Mpi_real8, down, 0, &
                       sgs_rcv_hi, nxg*nyg, Mpi_real8, up,   0, &
                       MPI_COMM_WORLD, istat, ierr )
    If (up /= MPI_PROC_NULL) F_(:,:,nzg) = sgs_rcv_hi(:,:,1)

  End Subroutine update_ghost_interior_planes_nut

End Module sgs_models
