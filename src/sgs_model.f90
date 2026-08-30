!> SGS eddy-viscosity models for LES (sgs_model: 0=DNS, 1=Vreman)
Module sgs_models

  Use iso_fortran_env, Only : Int32, Int64
  Use global
  Use mpi
  Use decomp, Only : z_halo_neighbors
  Use boundary_conditions, Only : apply_periodic_bc_z

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
    Real   (Int64) :: dx_c, dy_c, dz_c          ! local filter widths (IBM pass only)
    Real   (Int64) :: dx2, dz2                  ! uniform-grid squares
    Real   (Int64) :: inv_dx, inv_dz            ! 1/dx, 1/dz
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
    inv_dz = 1d0 / dz
    dx2    = dx * dx
    dz2    = dz * dz

    ! Pass 1: standard (non-IBM) stencil, unconditional so it stays GPU-offloadable without phi residency; ibm_active discards it via Pass 2/solid-cell zeroing below
    !$acc parallel loop collapse(2) present(U_,V_,W_,nu_t_,y,yg)
    Do k = 2, nzg-1
       Do j = 2, nyg-1

          ! Per-j constants: y-grid spacings (non-uniform y; x,z uniform)
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
             a31 = 0.5d0*( ( U_(i,j,k)   - U_(i,j,k-1)  ) + &
                           ( U_(i,j,k+1) - U_(i,j,k  )  ) ) * inv_dz
             a12 = 0.5d0*( ( V_(i,j,k)   - V_(i-1,j,k)  ) + &
                           ( V_(i+1,j,k) - V_(i,  j,k)  ) ) * inv_dx
             ! Non-uniform y: V on y-faces
             a22 = ( V_(i,j,k)   - V_(i,j-1,k)   ) * inv_dy_j
             a32 = 0.5d0*( ( V_(i,j,k)   - V_(i,j,k-1)  ) + &
                           ( V_(i,j,k+1) - V_(i,j,k  )  ) ) * inv_dz
             a13 = 0.5d0*( ( W_(i,j,k)   - W_(i-1,j,k)  ) + &
                           ( W_(i+1,j,k) - W_(i,  j,k)  ) ) * inv_dx
             ! Non-uniform y
             a23 = 0.5d0*( ( W_(i,j,k)   - W_(i,j-1,k)  ) * inv_yg_jm + &
                           ( W_(i,j+1,k) - W_(i,j,  k)  ) * inv_yg_jp )
             ! Uniform z: W on z-faces
             a33 = ( W_(i,j,k)   - W_(i,j,k-1)   ) * inv_dz

             ! ------ beta_mn (standard: dx_c=dx, dz_c=dz, dy_c=y(j)-y(j-1)) --
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

    ! Pass 2, wall zeroing, x-periodicity, and MPI halo exchange below all run on host; sync Pass 1's GPU result back first
    !$acc update host(nu_t_)

    ! Pass 2: IBM corrections (only when ibm_active); re-visits all fluid cells in the domain (O(volume), not O(surface)), does not affect Pass 1
    If ( ibm_active ) Then
       Do k = 2, nzg-1
          Do j = 2, nyg-1
             dy_c      = y(j) - y(j-1)
             inv_yg_jm = 1d0 / ( yg(j)   - yg(j-1) )
             inv_yg_jp = 1d0 / ( yg(j+1) - yg(j  ) )
             inv_dy_j  = 1d0 / dy_c

             Do i = 2, nxg-1
                If ( phi(i,j,k) <= 0d0 ) Cycle   ! solid cell: zeroed in post-processing

                ! Re-evaluate gradients; then apply one-sided corrections where needed
                a11 = ( U_(i,j,k)   - U_(i-1,j,k)   ) * inv_dx
                a21 = 0.5d0*( (U_(i,j,k)-U_(i,j-1,k))*inv_yg_jm + (U_(i,j+1,k)-U_(i,j,k))*inv_yg_jp )
                a31 = 0.5d0*( (U_(i,j,k)-U_(i,j,k-1)) + (U_(i,j,k+1)-U_(i,j,k)) ) * inv_dz
                a12 = 0.5d0*( (V_(i,j,k)-V_(i-1,j,k)) + (V_(i+1,j,k)-V_(i,j,k)) ) * inv_dx
                a22 = ( V_(i,j,k)   - V_(i,j-1,k)   ) * inv_dy_j
                a32 = 0.5d0*( (V_(i,j,k)-V_(i,j,k-1)) + (V_(i,j,k+1)-V_(i,j,k)) ) * inv_dz
                a13 = 0.5d0*( (W_(i,j,k)-W_(i-1,j,k)) + (W_(i+1,j,k)-W_(i,j,k)) ) * inv_dx
                a23 = 0.5d0*( (W_(i,j,k)-W_(i,j-1,k))*inv_yg_jm + (W_(i,j+1,k)-W_(i,j,k))*inv_yg_jp )
                a33 = ( W_(i,j,k)   - W_(i,j,k-1)   ) * inv_dz

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
                   a31 = ( U_(i,j,k+1) - U_(i,j,k) ) * inv_dz
                   a32 = ( V_(i,j,k+1) - V_(i,j,k) ) * inv_dz
                Else If ( Umask_cc(i,j,k+1) < 0.5d0 ) Then
                   a31 = ( U_(i,j,k)   - U_(i,j,k-1) ) * inv_dz
                   a32 = ( V_(i,j,k)   - V_(i,j,k-1) ) * inv_dz
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
                dz_c = Min(dz, 2d0*phi(i,j,k))

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
    nu_t_(1,  :,:) = nu_t_(nxg-1,:,:)
    nu_t_(nxg,:,:) = nu_t_(2,    :,:)

    ! Ring exchange for intermediate ranks (host-only); rank-0/rank-(nprocs-1) wrap handled below.
    Call update_ghost_interior_planes_nut(nu_t_)

    ! Push host state to device before apply_periodic_bc_z, which runs device-resident at nprocs==1 and would otherwise have its z-wrap fill clobbered by a later blanket update device
    !$acc update device(nu_t_)
    Call apply_periodic_bc_z(nu_t_, 4)
    ! Sync back to host for output_data's snapshot write and compute_wall_model's host-only bits
    !$acc update host(nu_t_)

  End Subroutine compute_vreman


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
