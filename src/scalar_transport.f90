!> Passive suspended-sediment scalar transport (settling, MUSCL advection, diffusion)
Module scalar_transport

  Use iso_fortran_env, Only : Int32, Int64
  Use global
  Use mpi
  Use decomp, Only : z_halo_neighbors, x_halo_neighbors, x_periodic_partner, z_periodic_partner
  Use boundary_conditions, Only : apply_periodic_bc_z, apply_periodic_bc_x, update_ghost_interior_planes_x, &
                                  apply_inflow_bc_scalar_x_C, outflow_convection_velocity, outflow_stage_dt

  Implicit None

  ! Module-level halo buffers for C z-exchange (avoids per-call heap allocation)
  Real(Int64), Allocatable, Dimension(:,:,:) :: sc_snd_lo, sc_snd_hi
  Real(Int64), Allocatable, Dimension(:,:,:) :: sc_rcv_lo, sc_rcv_hi

  ! Scalar copy with a second x and z ghost layer (indices 0 and nxg+1 / nzg+1) for the MUSCL far-upwind cell, so the scheme
  ! stays second order across MPI seams and periodic wraps instead of dropping to first order there
  Real(Int64), Allocatable, Dimension(:,:,:) :: Cpad
  ! MPI staging planes for scalar_fill_pad (device-resident in GPU builds; the exchange itself runs on the host copies)
  Real(Int64), Allocatable, Dimension(:,:) :: xs, xr, zs, zr

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

    Call scalar_fill_pad(C_)

    !$acc parallel loop collapse(3) present(Cpad,U_,V_,W_,Fc_,nu_t,phi,x,xg,y,yg,z,zg,weight_y_0,weight_y_1,weight_z_0,weight_z_1)
    Do k = 2, nzg-1
       Do j = 2, nyg-1
          Do i = 2, nxg-1

             ! x-advection (U_ on x-faces): MUSCL reconstruction

             ! -- high face (i+1/2) at x(i): U_(i,j,k)
             uf = U_(i,j,k)
             If ( uf >= 0d0 ) Then
                ! upwind cell = i
                gf  = ( Cpad(i+1,j,k) - Cpad(i,j,k)   ) / ( xg(i+1) - xg(i)   )
                gb  = ( Cpad(i,j,k)   - Cpad(i-1,j,k) ) / ( xg(i)   - xg(i-1) )
                slp = vanleer_slope(gf, gb)
                C_hi = Cpad(i,j,k) + slp * ( x(i) - xg(i) )
             Else
                ! upwind cell = i+1; far-upwind C(i+2) comes from the padded halo (x is uniform, so dx stands in for xg(i+2)-xg(i+1))
                gf  = ( Cpad(i+2,j,k) - Cpad(i+1,j,k) ) / dx
                gb  = ( Cpad(i+1,j,k) - Cpad(i,j,k)   ) / ( xg(i+1) - xg(i)   )
                slp = vanleer_slope(gf, gb)
                C_hi = Cpad(i+1,j,k) + slp * ( x(i) - xg(i+1) )
             End If

             ! -- low face (i-1/2) at x(i-1): U_(i-1,j,k)
             uf = U_(i-1,j,k)
             If ( uf >= 0d0 ) Then
                ! upwind cell = i-1; far-upwind C(i-2) comes from the padded halo
                gf  = ( Cpad(i,j,k)   - Cpad(i-1,j,k) ) / ( xg(i)   - xg(i-1) )
                gb  = ( Cpad(i-1,j,k) - Cpad(i-2,j,k) ) / dx
                slp = vanleer_slope(gf, gb)
                C_lo = Cpad(i-1,j,k) + slp * ( x(i-1) - xg(i-1) )
             Else
                ! upwind cell = i
                gf  = ( Cpad(i+1,j,k) - Cpad(i,j,k)   ) / ( xg(i+1) - xg(i)   )
                gb  = ( Cpad(i,j,k)   - Cpad(i-1,j,k) ) / ( xg(i)   - xg(i-1) )
                slp = vanleer_slope(gf, gb)
                C_lo = Cpad(i,j,k) + slp * ( x(i-1) - xg(i) )
             End If

             dx_f  = x(i) - x(i-1)
             adv_x = ( U_(i,j,k)*C_hi - U_(i-1,j,k)*C_lo ) / dx_f

             !-------- y-advection (V_ - ws on y-faces) --------------------
             ! -- high face (j+1/2) at y(j): V_(i,j,k); subtract ws
             vf = V_(i,j,k) - w_settle
             If ( vf >= 0d0 ) Then
                ! upwind cell = j
                gf  = ( Cpad(i,j+1,k) - Cpad(i,j,k)   ) / ( yg(j+1) - yg(j)   )
                gb  = ( Cpad(i,j,k)   - Cpad(i,j-1,k) ) / ( yg(j)   - yg(j-1) )
                slp = vanleer_slope(gf, gb)
                C_hi = Cpad(i,j,k) + slp * ( y(j) - yg(j) )
             Else
                ! upwind cell = j+1; guard C(i,j+2) out of bounds at j=nyg-1
                If ( j == nyg-1 ) Then
                   C_hi = Cpad(i,j+1,k)
                Else
                   gf  = ( Cpad(i,j+2,k) - Cpad(i,j+1,k) ) / ( yg(j+2) - yg(j+1) )
                   gb  = ( Cpad(i,j+1,k) - Cpad(i,j,k)   ) / ( yg(j+1) - yg(j)   )
                   slp = vanleer_slope(gf, gb)
                   C_hi = Cpad(i,j+1,k) + slp * ( y(j) - yg(j+1) )
                End If
             End If

             ! -- low face (j-1/2) at y(j-1): V_(i,j-1,k); subtract ws
             vf = V_(i,j-1,k) - w_settle
             If ( vf >= 0d0 ) Then
                ! upwind cell = j-1; guard C(i,j-2) out of bounds at j=2
                If ( j == 2 ) Then
                   C_lo = Cpad(i,j-1,k)
                Else
                   gf  = ( Cpad(i,j,k)   - Cpad(i,j-1,k) ) / ( yg(j)   - yg(j-1) )
                   gb  = ( Cpad(i,j-1,k) - Cpad(i,j-2,k) ) / ( yg(j-1) - yg(j-2) )
                   slp = vanleer_slope(gf, gb)
                   C_lo = Cpad(i,j-1,k) + slp * ( y(j-1) - yg(j-1) )
                End If
             Else
                ! upwind cell = j
                gf  = ( Cpad(i,j+1,k) - Cpad(i,j,k)   ) / ( yg(j+1) - yg(j)   )
                gb  = ( Cpad(i,j,k)   - Cpad(i,j-1,k) ) / ( yg(j)   - yg(j-1) )
                slp = vanleer_slope(gf, gb)
                C_lo = Cpad(i,j,k) + slp * ( y(j-1) - yg(j) )
             End If

             dy_f  = y(j) - y(j-1)
             adv_y = ( (V_(i,j,k) - w_settle)*C_hi - (V_(i,j-1,k) - w_settle)*C_lo ) / dy_f

             !-------- z-advection (W_ on z-faces) -------------------------
             ! -- high face (k+1/2) at z(k): W_(i,j,k)
             wf = W_(i,j,k)
             If ( wf >= 0d0 ) Then
                ! upwind cell = k
                gf  = ( Cpad(i,j,k+1) - Cpad(i,j,k)   ) / ( zg(k+1) - zg(k)   )
                gb  = ( Cpad(i,j,k)   - Cpad(i,j,k-1) ) / ( zg(k)   - zg(k-1) )
                slp = vanleer_slope(gf, gb)
                C_hi = Cpad(i,j,k) + slp * ( z(k) - zg(k) )
             Else
                ! upwind cell = k+1; far-upwind C(k+2) comes from the padded halo (scalars need uniform z, so dz stands in for zg(k+2)-zg(k+1))
                gf  = ( Cpad(i,j,k+2) - Cpad(i,j,k+1) ) / dz
                gb  = ( Cpad(i,j,k+1) - Cpad(i,j,k)   ) / ( zg(k+1) - zg(k)   )
                slp = vanleer_slope(gf, gb)
                C_hi = Cpad(i,j,k+1) + slp * ( z(k) - zg(k+1) )
             End If

             ! -- low face (k-1/2) at z(k-1): W_(i,j,k-1)
             wf = W_(i,j,k-1)
             If ( wf >= 0d0 ) Then
                ! upwind cell = k-1; far-upwind C(k-2) comes from the padded halo
                gf  = ( Cpad(i,j,k)   - Cpad(i,j,k-1) ) / ( zg(k)   - zg(k-1) )
                gb  = ( Cpad(i,j,k-1) - Cpad(i,j,k-2) ) / dz
                slp = vanleer_slope(gf, gb)
                C_lo = Cpad(i,j,k-1) + slp * ( z(k-1) - zg(k-1) )
             Else
                ! upwind cell = k
                gf  = ( Cpad(i,j,k+1) - Cpad(i,j,k)   ) / ( zg(k+1) - zg(k)   )
                gb  = ( Cpad(i,j,k)   - Cpad(i,j,k-1) ) / ( zg(k)   - zg(k-1) )
                slp = vanleer_slope(gf, gb)
                C_lo = Cpad(i,j,k) + slp * ( z(k-1) - zg(k) )
             End If

             dz_f  = z(k) - z(k-1)
             adv_z = ( W_(i,j,k)*C_hi - W_(i,j,k-1)*C_lo ) / dz_f

             !-------- x-diffusion -----------------------------------------
             ! kappa at high x-face (x uniform: plain average; y, z use the interpolation weights): average of nu_t at (i,j,k) and (i+1,j,k)
             kappa_hi = kappa_mol + 0.5d0*(nu_t(i,j,k) + nu_t(i+1,j,k))*kappa_t_inv
             kappa_lo = kappa_mol + 0.5d0*(nu_t(i,j,k) + nu_t(i-1,j,k))*kappa_t_inv
             diff_x = ( kappa_hi*(Cpad(i+1,j,k) - Cpad(i,j,k))/(xg(i+1)-xg(i)) &
                      - kappa_lo*(Cpad(i,j,k) - Cpad(i-1,j,k))/(xg(i)-xg(i-1)) ) / dx_f

             !-------- y-diffusion -----------------------------------------
             kappa_hi = kappa_mol + ( weight_y_0(j  )*nu_t(i,j,k) + weight_y_1(j  )*nu_t(i,j+1,k) )*kappa_t_inv
             kappa_lo = kappa_mol + ( weight_y_0(j-1)*nu_t(i,j-1,k) + weight_y_1(j-1)*nu_t(i,j,k) )*kappa_t_inv
             diff_y = ( kappa_hi*(Cpad(i,j+1,k) - Cpad(i,j,k))/(yg(j+1)-yg(j)) &
                      - kappa_lo*(Cpad(i,j,k) - Cpad(i,j-1,k))/(yg(j)-yg(j-1)) ) / dy_f

             !-------- z-diffusion -----------------------------------------
             kappa_hi = kappa_mol + ( weight_z_0(k  )*nu_t(i,j,k) + weight_z_1(k  )*nu_t(i,j,k+1) )*kappa_t_inv
             kappa_lo = kappa_mol + ( weight_z_0(k-1)*nu_t(i,j,k-1) + weight_z_1(k-1)*nu_t(i,j,k) )*kappa_t_inv
             diff_z = ( kappa_hi*(Cpad(i,j,k+1) - Cpad(i,j,k))/(zg(k+1)-zg(k)) &
                      - kappa_lo*(Cpad(i,j,k) - Cpad(i,j,k-1))/(zg(k)-zg(k-1)) ) / dz_f

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


  !> Load C_ into the padded work array Cpad and fill its outer x/z ghost layer: neighbour-rank cells across MPI seams, the periodic partner's cells across a periodic wrap, else a copy of the adjacent ghost cell (zero slope, i.e. first-order upwind at true domain edges)
  Subroutine scalar_fill_pad(C_)

    Real(Int64), Dimension(nxg, nyg, nzg), Intent(In) :: C_

    Integer(Int32) :: i, j, k, up, down, xsend_up, xsend_dn, zsend_up, zsend_dn
    Logical        :: is_first_x, is_last_x, is_first_z, is_last_z
    Integer(Int32) :: partner_x, partner_z

    If ( .Not. Allocated(Cpad) ) Then
       Allocate( Cpad(0:nxg+1, 1:nyg, 0:nzg+1) )
       Allocate( xs(nyg,nzg), xr(nyg,nzg), zs(nxg+2,nyg), zr(nxg+2,nyg) )
       !$acc enter data create(Cpad,xs,xr,zs,zr)
    End If

    !$acc parallel loop collapse(3) present(C_,Cpad)
    Do k = 1, nzg
       Do j = 1, nyg
          Do i = 1, nxg
             Cpad(i,j,k) = C_(i,j,k)
          End Do
       End Do
    End Do
    !$acc end parallel loop

    Call x_periodic_partner(is_first_x, is_last_x, partner_x)
    Call z_periodic_partner(is_first_z, is_last_z, partner_z)

    If ( is_first_x .And. is_last_x .And. is_first_z .And. is_last_z ) Then
       ! single rank: every outer plane is a local wrap or edge copy (device-resident on GPU builds)
       !$acc kernels present(Cpad)
       If ( x_bc_type == 0 ) Then
          Cpad(0,   :,1:nzg) = Cpad(nxg-3,:,1:nzg)
          Cpad(nxg+1,:,1:nzg) = Cpad(4,   :,1:nzg)
       Else
          Cpad(0,   :,1:nzg) = Cpad(1,  :,1:nzg)
          Cpad(nxg+1,:,1:nzg) = Cpad(nxg,:,1:nzg)
       End If
       Cpad(:,:,0)     = Cpad(:,:,nzg-3)
       Cpad(:,:,nzg+1) = Cpad(:,:,4)
       !$acc end kernels
       Return
    End If

    ! multi-rank: x planes first for k=1..nzg, then z planes over the x-extended extent so corners are consistent. The MPI
    ! exchange runs on host copies; in GPU builds only the planes that are sent or received go through the host.
    Call x_halo_neighbors(up, down)
    xsend_up = nxg-2
    xsend_dn = 3
    If ( x_bc_type == 0 ) Then
       If ( is_last_x  ) Then
          up = partner_x
          xsend_up = nxg-3
       End If
       If ( is_first_x ) Then
          down = partner_x
          xsend_dn = 4
       End If
    End If

    ! x planes: pack on the device into contiguous buffers, stage those through the host for MPI, unpack on the device
    !$acc kernels present(Cpad,xs)
    xs = Cpad(xsend_up,:,1:nzg)
    !$acc end kernels
    !$acc update host(xs)
    Call Mpi_Sendrecv( xs, nyg*nzg, Mpi_real8, up,   20, &
                       xr, nyg*nzg, Mpi_real8, down, 20, MPI_COMM_WORLD, istat, ierr )
    !$acc update device(xr)
    If ( down /= MPI_PROC_NULL ) Then
       !$acc kernels present(Cpad,xr)
       Cpad(0,:,1:nzg) = xr
       !$acc end kernels
    Else
       !$acc kernels present(Cpad)
       Cpad(0,:,1:nzg) = Cpad(1,:,1:nzg)
       !$acc end kernels
    End If
    !$acc kernels present(Cpad,xs)
    xs = Cpad(xsend_dn,:,1:nzg)
    !$acc end kernels
    !$acc update host(xs)
    Call Mpi_Sendrecv( xs, nyg*nzg, Mpi_real8, down, 21, &
                       xr, nyg*nzg, Mpi_real8, up,   21, MPI_COMM_WORLD, istat, ierr )
    !$acc update device(xr)
    If ( up /= MPI_PROC_NULL ) Then
       !$acc kernels present(Cpad,xr)
       Cpad(nxg+1,:,1:nzg) = xr
       !$acc end kernels
    Else
       !$acc kernels present(Cpad)
       Cpad(nxg+1,:,1:nzg) = Cpad(nxg,:,1:nzg)
       !$acc end kernels
    End If

    ! z: scalars are z-periodic only, so first/last columns always wrap
    Call z_halo_neighbors(up, down)
    zsend_up = nzg-2
    zsend_dn = 3
    If ( is_last_z  ) Then
       up = partner_z
       zsend_up = nzg-3
    End If
    If ( is_first_z ) Then
       down = partner_z
       zsend_dn = 4
    End If

    !$acc kernels present(Cpad,zs)
    zs = Cpad(:,:,zsend_up)
    !$acc end kernels
    !$acc update host(zs)
    Call Mpi_Sendrecv( zs, (nxg+2)*nyg, Mpi_real8, up,   22, &
                       zr, (nxg+2)*nyg, Mpi_real8, down, 22, MPI_COMM_WORLD, istat, ierr )
    !$acc update device(zr)
    !$acc kernels present(Cpad,zr)
    Cpad(:,:,0) = zr
    !$acc end kernels
    !$acc kernels present(Cpad,zs)
    zs = Cpad(:,:,zsend_dn)
    !$acc end kernels
    !$acc update host(zs)
    Call Mpi_Sendrecv( zs, (nxg+2)*nyg, Mpi_real8, down, 23, &
                       zr, (nxg+2)*nyg, Mpi_real8, up,   23, MPI_COMM_WORLD, istat, ierr )
    !$acc update device(zr)
    !$acc kernels present(Cpad,zr)
    Cpad(:,:,nzg+1) = zr
    !$acc end kernels

  End Subroutine scalar_fill_pad


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


  !> z-halo exchange plus the periodic x/z wraps for a cell-centred scalar, after the caller's x seam exchange and inflow/outflow update; host-side (runs the wraps on the device via a round trip)
  Subroutine finish_scalar_halos(C_)

    Real(Int64), Dimension(nxg, nyg, nzg), Intent(InOut) :: C_

    ! z-halo via MPI (ring exchange, non-periodic)
    Call update_ghost_scalar(C_)
    ! Push host state to device first: the periodic wraps run device-resident at nprocs==1, else their fill is clobbered by the caller's later blanket update device
    !$acc update device(C_)
    If ( x_bc_type == 0 ) Call apply_periodic_bc_x(C_, 2)
    Call apply_periodic_bc_z(C_, 4)
    !$acc update host(C_)

  End Subroutine finish_scalar_halos


  !> Apply boundary conditions to scalar C (x-periodic or Dirichlet inflow/convective outflow, z-MPI halo, y-wall fluxes)
  Subroutine apply_scalar_bc(C_)

    Real(Int64), Dimension(nxg, nyg, nzg), Intent(InOut) :: C_

    Real(Int64) :: Uc, courant
    Logical :: is_first_x, is_last_x
    Integer(Int32) :: partner_x

    ! x direction: interior-rank seam planes, then Dirichlet inflow / convective outflow (periodic wrap comes in finish_scalar_halos)
    ! the seam exchange runs on device memory in multi-GPU builds: push the host state, exchange, bring the planes back
    !$acc update device(C_)
    Call update_ghost_interior_planes_x(C_, 4)
    !$acc update host(C_)
    If ( x_bc_type == 1 ) Then
       Call x_periodic_partner(is_first_x, is_last_x, partner_x)
       Call apply_inflow_bc_scalar_x_C(C_)
       Uc = outflow_convection_velocity()
       courant = Min(Max(Uc,0d0)*outflow_stage_dt()/dx, 1d0)
       If ( is_last_x ) C_(nxg,:,:) = C_(nxg,:,:) - courant*( C_(nxg,:,:) - C_(nxg-1,:,:) )
    End If

    Call finish_scalar_halos(C_)

    ! y-bottom ghost
    If ( sed_bc_bot == 0 ) Then
       C_(:,1,:) = C_(:,2,:)
    Else
       C_(:,1,:) = 2d0*C_ref - C_(:,2,:)
    End If

    ! y-top ghost
    C_(:,nyg,:) = C_(:,nyg-1,:)

    ! IBM solid cells are deliberately not zeroed: image-point interpolation and the MUSCL stencil read them, and zeros would drain nearby fluid; the no-flux ghost fill (apply_ghost_cell_ibm_scalar_noflux) supplies the boundary condition

  End Subroutine apply_scalar_bc


End Module scalar_transport
