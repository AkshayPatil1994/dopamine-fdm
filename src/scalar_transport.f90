!> Passive suspended-sediment scalar transport (settling, MUSCL advection, diffusion)
Module scalar_transport

  Use iso_fortran_env, Only : Int32, Int64
  Use global
  Use mpi
  Use decomp, Only : z_halo_neighbors
  Use boundary_conditions, Only : apply_periodic_bc_z, apply_inflow_bc_scalar_x_C, outflow_convection_velocity

  Implicit None

  ! Module-level halo buffers for C z-exchange (avoids per-call heap allocation)
  Real(Int64), Allocatable, Dimension(:,:,:) :: sc_snd_lo, sc_snd_hi
  Real(Int64), Allocatable, Dimension(:,:,:) :: sc_rcv_lo, sc_rcv_hi

Contains

  !> Soulsby (1997) settling velocity, stored in ws (global)
  Subroutine compute_settling_velocity

    Real(Int64) :: Dstar

    Dstar = d_s * ( ((rho_s/rho_f) - 1d0) * grav / nu**2 )**(1d0/3d0)
    ws    = (nu / d_s) * ( Sqrt(10.36d0**2 + 1.049d0*Dstar**3) - 10.36d0 )

  End Subroutine compute_settling_velocity


  !> van Leer (1974) harmonic-mean slope limiter on a non-uniform grid
  Pure Function vanleer_slope(gf, gb) Result(sigma)
    !$acc routine seq

    Real(Int64), Intent(In) :: gf, gb
    Real(Int64)             :: sigma

    If ( gf*gb > 0d0 ) Then
       sigma = 2d0*gf*gb / (gf + gb)
    Else
       sigma = 0d0
    End If

  End Function vanleer_slope


  !> Compute RHS for scalar advection-diffusion (MUSCL advection + central diffusion), Sc/Sc_t diffusivity, with settling velocity
  Subroutine compute_rhs_scalar(C_, U_, V_, W_, Fc_)

    Real(Int64), Dimension(nxg, nyg, nzg), Intent(In)    :: C_
    Real(Int64), Dimension(nx,  nyg, nzg), Intent(In)    :: U_
    Real(Int64), Dimension(nxg, ny,  nzg), Intent(In)    :: V_
    Real(Int64), Dimension(nxg, nyg, nz ), Intent(In)    :: W_
    Real(Int64), Dimension(2:nxg-1, 2:nyg-1, 2:nzg-1), Intent(Out) :: Fc_

    Call compute_rhs_scalar_core(C_, U_, V_, W_, ws, nu/Sc, 1d0/Sc_t, Fc_)

  End Subroutine compute_rhs_scalar


  !> Shared MUSCL advection-diffusion core for any cell-centred scalar (sediment concentration, temperature, ...)
  Subroutine compute_rhs_scalar_core(C_, U_, V_, W_, w_settle, kappa_mol, kappa_t_inv, Fc_)

    Real(Int64), Dimension(nxg, nyg, nzg), Intent(In)    :: C_
    Real(Int64), Dimension(nx,  nyg, nzg), Intent(In)    :: U_
    Real(Int64), Dimension(nxg, ny,  nzg), Intent(In)    :: V_
    Real(Int64), Dimension(nxg, nyg, nz ), Intent(In)    :: W_
    Real(Int64), Intent(In) :: w_settle, kappa_mol, kappa_t_inv
    Real(Int64), Dimension(2:nxg-1, 2:nyg-1, 2:nzg-1), Intent(Out) :: Fc_

    Integer(Int32) :: i, j, k
    Real   (Int64) :: adv_x, adv_y, adv_z
    Real   (Int64) :: diff_x, diff_y, diff_z
    Real   (Int64) :: uf, vf, wf                   ! face velocity (incl. settling)
    Real   (Int64) :: kappa_lo, kappa_hi            ! face diffusivities
    Real   (Int64) :: C_lo, C_hi                    ! TVD reconstructed face values
    Real   (Int64) :: gf, gb, slp                   ! face gradients & limited slope
    Real   (Int64) :: dx_f, dy_f, dz_f

    !$acc parallel loop collapse(3) present(C_,U_,V_,W_,Fc_,nu_t,phi,x,xg,y,yg,z,zg)
    Do k = 2, nzg-1
       Do j = 2, nyg-1
          Do i = 2, nxg-1

             ! x-advection (U_ on x-faces): MUSCL reconstruction

             ! -- high face (i+1/2) at x(i): U_(i,j,k)
             uf = U_(i,j,k)
             If ( uf >= 0d0 ) Then
                ! upwind cell = i
                gf  = ( C_(i+1,j,k) - C_(i,j,k)   ) / ( xg(i+1) - xg(i)   )
                gb  = ( C_(i,j,k)   - C_(i-1,j,k) ) / ( xg(i)   - xg(i-1) )
                slp = vanleer_slope(gf, gb)
                C_hi = C_(i,j,k) + slp * ( x(i) - xg(i) )
             Else
                ! upwind cell = i+1; guard far-upwind C_(i+2) out of bounds at i=nxg-1
                If ( i == nxg-1 ) Then
                   C_hi = C_(i+1,j,k)
                Else
                   gf  = ( C_(i+2,j,k) - C_(i+1,j,k) ) / ( xg(i+2) - xg(i+1) )
                   gb  = ( C_(i+1,j,k) - C_(i,j,k)   ) / ( xg(i+1) - xg(i)   )
                   slp = vanleer_slope(gf, gb)
                   C_hi = C_(i+1,j,k) + slp * ( x(i) - xg(i+1) )
                End If
             End If

             ! -- low face (i-1/2) at x(i-1): U_(i-1,j,k)
             uf = U_(i-1,j,k)
             If ( uf >= 0d0 ) Then
                ! upwind cell = i-1; guard far-upwind C_(i-2) out of bounds at i=2
                If ( i == 2 ) Then
                   C_lo = C_(i-1,j,k)
                Else
                   gf  = ( C_(i,j,k)   - C_(i-1,j,k) ) / ( xg(i)   - xg(i-1) )
                   gb  = ( C_(i-1,j,k) - C_(i-2,j,k) ) / ( xg(i-1) - xg(i-2) )
                   slp = vanleer_slope(gf, gb)
                   C_lo = C_(i-1,j,k) + slp * ( x(i-1) - xg(i-1) )
                End If
             Else
                ! upwind cell = i
                gf  = ( C_(i+1,j,k) - C_(i,j,k)   ) / ( xg(i+1) - xg(i)   )
                gb  = ( C_(i,j,k)   - C_(i-1,j,k) ) / ( xg(i)   - xg(i-1) )
                slp = vanleer_slope(gf, gb)
                C_lo = C_(i,j,k) + slp * ( x(i-1) - xg(i) )
             End If

             dx_f  = x(i) - x(i-1)
             adv_x = ( U_(i,j,k)*C_hi - U_(i-1,j,k)*C_lo ) / dx_f

             !-------- y-advection (V_ - ws on y-faces) --------------------
             ! -- high face (j+1/2) at y(j): V_(i,j,k); subtract ws
             vf = V_(i,j,k) - w_settle
             If ( vf >= 0d0 ) Then
                ! upwind cell = j
                gf  = ( C_(i,j+1,k) - C_(i,j,k)   ) / ( yg(j+1) - yg(j)   )
                gb  = ( C_(i,j,k)   - C_(i,j-1,k) ) / ( yg(j)   - yg(j-1) )
                slp = vanleer_slope(gf, gb)
                C_hi = C_(i,j,k) + slp * ( y(j) - yg(j) )
             Else
                ! upwind cell = j+1; guard C_(i,j+2) out of bounds at j=nyg-1
                If ( j == nyg-1 ) Then
                   C_hi = C_(i,j+1,k)
                Else
                   gf  = ( C_(i,j+2,k) - C_(i,j+1,k) ) / ( yg(j+2) - yg(j+1) )
                   gb  = ( C_(i,j+1,k) - C_(i,j,k)   ) / ( yg(j+1) - yg(j)   )
                   slp = vanleer_slope(gf, gb)
                   C_hi = C_(i,j+1,k) + slp * ( y(j) - yg(j+1) )
                End If
             End If

             ! -- low face (j-1/2) at y(j-1): V_(i,j-1,k); subtract ws
             vf = V_(i,j-1,k) - w_settle
             If ( vf >= 0d0 ) Then
                ! upwind cell = j-1; guard C_(i,j-2) out of bounds at j=2
                If ( j == 2 ) Then
                   C_lo = C_(i,j-1,k)
                Else
                   gf  = ( C_(i,j,k)   - C_(i,j-1,k) ) / ( yg(j)   - yg(j-1) )
                   gb  = ( C_(i,j-1,k) - C_(i,j-2,k) ) / ( yg(j-1) - yg(j-2) )
                   slp = vanleer_slope(gf, gb)
                   C_lo = C_(i,j-1,k) + slp * ( y(j-1) - yg(j-1) )
                End If
             Else
                ! upwind cell = j
                gf  = ( C_(i,j+1,k) - C_(i,j,k)   ) / ( yg(j+1) - yg(j)   )
                gb  = ( C_(i,j,k)   - C_(i,j-1,k) ) / ( yg(j)   - yg(j-1) )
                slp = vanleer_slope(gf, gb)
                C_lo = C_(i,j,k) + slp * ( y(j-1) - yg(j) )
             End If

             dy_f  = y(j) - y(j-1)
             adv_y = ( (V_(i,j,k) - w_settle)*C_hi - (V_(i,j-1,k) - w_settle)*C_lo ) / dy_f

             !-------- z-advection (W_ on z-faces) -------------------------
             ! -- high face (k+1/2) at z(k): W_(i,j,k)
             wf = W_(i,j,k)
             If ( wf >= 0d0 ) Then
                ! upwind cell = k
                gf  = ( C_(i,j,k+1) - C_(i,j,k)   ) / ( zg(k+1) - zg(k)   )
                gb  = ( C_(i,j,k)   - C_(i,j,k-1) ) / ( zg(k)   - zg(k-1) )
                slp = vanleer_slope(gf, gb)
                C_hi = C_(i,j,k) + slp * ( z(k) - zg(k) )
             Else
                ! upwind cell = k+1; guard C_(i,j,k+2) out of bounds at k=nzg-1
                If ( k == nzg-1 ) Then
                   C_hi = C_(i,j,k+1)
                Else
                   gf  = ( C_(i,j,k+2) - C_(i,j,k+1) ) / ( zg(k+2) - zg(k+1) )
                   gb  = ( C_(i,j,k+1) - C_(i,j,k)   ) / ( zg(k+1) - zg(k)   )
                   slp = vanleer_slope(gf, gb)
                   C_hi = C_(i,j,k+1) + slp * ( z(k) - zg(k+1) )
                End If
             End If

             ! -- low face (k-1/2) at z(k-1): W_(i,j,k-1)
             wf = W_(i,j,k-1)
             If ( wf >= 0d0 ) Then
                ! upwind cell = k-1; guard C_(i,j,k-2) out of bounds at k=2
                If ( k == 2 ) Then
                   C_lo = C_(i,j,k-1)
                Else
                   gf  = ( C_(i,j,k)   - C_(i,j,k-1) ) / ( zg(k)   - zg(k-1) )
                   gb  = ( C_(i,j,k-1) - C_(i,j,k-2) ) / ( zg(k-1) - zg(k-2) )
                   slp = vanleer_slope(gf, gb)
                   C_lo = C_(i,j,k-1) + slp * ( z(k-1) - zg(k-1) )
                End If
             Else
                ! upwind cell = k
                gf  = ( C_(i,j,k+1) - C_(i,j,k)   ) / ( zg(k+1) - zg(k)   )
                gb  = ( C_(i,j,k)   - C_(i,j,k-1) ) / ( zg(k)   - zg(k-1) )
                slp = vanleer_slope(gf, gb)
                C_lo = C_(i,j,k) + slp * ( z(k-1) - zg(k) )
             End If

             dz_f  = z(k) - z(k-1)
             adv_z = ( W_(i,j,k)*C_hi - W_(i,j,k-1)*C_lo ) / dz_f

             !-------- x-diffusion -----------------------------------------
             ! kappa at high x-face: average of nu_t at (i,j,k) and (i+1,j,k)
             kappa_hi = kappa_mol + 0.5d0*(nu_t(i,j,k) + nu_t(i+1,j,k))*kappa_t_inv
             kappa_lo = kappa_mol + 0.5d0*(nu_t(i,j,k) + nu_t(i-1,j,k))*kappa_t_inv
             diff_x = ( kappa_hi*(C_(i+1,j,k) - C_(i,j,k))/(xg(i+1)-xg(i)) &
                      - kappa_lo*(C_(i,j,k) - C_(i-1,j,k))/(xg(i)-xg(i-1)) ) / dx_f

             !-------- y-diffusion -----------------------------------------
             kappa_hi = kappa_mol + 0.5d0*(nu_t(i,j,k) + nu_t(i,j+1,k))*kappa_t_inv
             kappa_lo = kappa_mol + 0.5d0*(nu_t(i,j,k) + nu_t(i,j-1,k))*kappa_t_inv
             diff_y = ( kappa_hi*(C_(i,j+1,k) - C_(i,j,k))/(yg(j+1)-yg(j)) &
                      - kappa_lo*(C_(i,j,k) - C_(i,j-1,k))/(yg(j)-yg(j-1)) ) / dy_f

             !-------- z-diffusion -----------------------------------------
             kappa_hi = kappa_mol + 0.5d0*(nu_t(i,j,k) + nu_t(i,j,k+1))*kappa_t_inv
             kappa_lo = kappa_mol + 0.5d0*(nu_t(i,j,k) + nu_t(i,j,k-1))*kappa_t_inv
             diff_z = ( kappa_hi*(C_(i,j,k+1) - C_(i,j,k))/(zg(k+1)-zg(k)) &
                      - kappa_lo*(C_(i,j,k) - C_(i,j,k-1))/(zg(k)-zg(k-1)) ) / dz_f

             !-------- Assemble RHS ----------------------------------------
             Fc_(i,j,k) = -adv_x - adv_y - adv_z + diff_x + diff_y + diff_z

             !-------- IBM gate: zero RHS inside solid --------------------
             If ( ibm_input_mode >= 1 ) Then
                If ( phi(i,j,k) < 0d0 ) Fc_(i,j,k) = 0d0
             End If

          End Do
       End Do
    End Do
    !$acc end parallel loop

  End Subroutine compute_rhs_scalar_core


  !> MPI ring exchange for C ghost planes (intermediate ranks only)
  Subroutine update_ghost_scalar(C_)

    Real(Int64), Dimension(nxg, nyg, nzg), Intent(InOut) :: C_

    Integer(Int32) :: up, down

    ! Allocate persistent single-plane buffers on first call
    If (.Not. Allocated(sc_snd_lo)) Then
       Allocate( sc_snd_lo(nxg,nyg,1), sc_snd_hi(nxg,nyg,1) )
       Allocate( sc_rcv_lo(nxg,nyg,1), sc_rcv_hi(nxg,nyg,1) )
    End If

    Call z_halo_neighbors(up, down)

    !-- Exchange A: send nzg-1 towards +z; receive z=1 from -z -------
    sc_snd_lo(:,:,1) = C_(:,:,nzg-1)
    Call Mpi_Sendrecv( sc_snd_lo, nxg*nyg, Mpi_real8, up,   0, &
                       sc_rcv_lo, nxg*nyg, Mpi_real8, down, 0, &
                       MPI_COMM_WORLD, istat, ierr )
    If (down /= MPI_PROC_NULL) C_(:,:,1) = sc_rcv_lo(:,:,1)

    !-- Exchange B: send z=2 towards -z; receive z=nzg from +z -----
    sc_snd_hi(:,:,1) = C_(:,:,2)
    Call Mpi_Sendrecv( sc_snd_hi, nxg*nyg, Mpi_real8, down, 0, &
                       sc_rcv_hi, nxg*nyg, Mpi_real8, up,   0, &
                       MPI_COMM_WORLD, istat, ierr )
    If (up /= MPI_PROC_NULL) C_(:,:,nzg) = sc_rcv_hi(:,:,1)

  End Subroutine update_ghost_scalar


  !> Apply boundary conditions to scalar C (x-periodic or Dirichlet inflow/convective outflow, z-MPI halo, y-wall fluxes)
  Subroutine apply_scalar_bc(C_)

    Real(Int64), Dimension(nxg, nyg, nzg), Intent(InOut) :: C_

    Real(Int64) :: Uc, courant

    ! x direction: periodic, or Dirichlet inflow / convective outflow
    If ( x_bc_type == 0 ) Then
       C_(1,   :,:) = C_(nxg-2,:,:)
       C_(nxg-1,:,:) = C_(2,   :,:)
       C_(nxg,  :,:) = C_(3,   :,:)
    Else
       Call apply_inflow_bc_scalar_x_C(C_)
       Uc = outflow_convection_velocity()
       courant = Min(Max(Uc,0d0)*dt/dx, 1d0)
       C_(nxg,:,:) = C_(nxg,:,:) - courant*( C_(nxg,:,:) - C_(nxg-1,:,:) )
    End If

    ! z-halo via MPI (ring exchange, non-periodic)
    Call update_ghost_scalar(C_)
    ! z-periodic wrap: rank 0 ↔ rank nprocs-1 (sets z=1 on rank 0 and
    ! fills the periodic-copy cell nzg-1 and ghost nzg on rank nprocs-1)
    ! Push host state to device first: apply_periodic_bc_z runs device-resident at nprocs==1, else its z-wrap fill is clobbered by the caller's later blanket update device
    !$acc update device(C_)
    Call apply_periodic_bc_z(C_, 4)
    !$acc update host(C_)

    ! y-bottom ghost
    If ( sed_bc_bot == 0 ) Then
       C_(:,1,:) = C_(:,2,:)
    Else
       C_(:,1,:) = 2d0*C_ref - C_(:,2,:)
    End If

    ! y-top ghost
    C_(:,nyg,:) = C_(:,nyg-1,:)

    ! Zero C inside IBM solid cells after every ghost/BC update to prevent
    ! the RK3 update leaking a non-zero IC value into adjacent fluid.
    If ( ibm_input_mode >= 1 .And. Allocated(phi) ) Then
       Where ( phi(2:nxg-1, 2:nyg-1, 2:nzg-1) <= 0d0 )
          C_(2:nxg-1, 2:nyg-1, 2:nzg-1) = 0d0
       End Where
    End If

  End Subroutine apply_scalar_bc


End Module scalar_transport
