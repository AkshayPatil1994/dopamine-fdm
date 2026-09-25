!> Passive suspended-sediment scalar transport (settling, MUSCL advection, diffusion)
Module scalar_transport

  Use iso_fortran_env, Only : Int32, Int64
  Use global
  Use mpi
  Use decomp, Only : z_halo_neighbors, x_halo_neighbors, z_periodic_partner, x_periodic_partner
  Use boundary_conditions, Only : apply_periodic_bc_x, apply_periodic_bc_z, update_ghost_interior_planes_x, &
                                  apply_inflow_bc_scalar_x_C, outflow_convection_velocity

  Implicit None

  ! Module-level halo buffers for C z-exchange (avoids per-call heap allocation)
  Real(Int64), Allocatable, Dimension(:,:,:) :: sc_snd_lo, sc_snd_hi
  Real(Int64), Allocatable, Dimension(:,:,:) :: sc_rcv_lo, sc_rcv_hi

  ! Far planes for the MUSCL reconstruction: the second cell beyond a face at the local array edge (local index 0 and
  ! n+1), taken from the neighbouring rank or the periodic partner, so the stencil does not depend on the decomposition.
  ! have_far_* = .False. at a true non-periodic domain boundary (first-order there); dfar_* = centre spacing from the
  ! ghost cell to the far cell.
  Real(Int64), Allocatable, Dimension(:,:) :: far_xlo, far_xhi, far_zlo, far_zhi
  Real(Int64), Allocatable, Dimension(:,:) :: fsend_xdn, fsend_xup, fsend_zdn, fsend_zup
  Logical        :: have_far_xlo = .False., have_far_xhi = .False., have_far_zlo = .False., have_far_zhi = .False.
  Real(Int64)    :: dfar_xlo = 1d0, dfar_xhi = 1d0, dfar_zlo = 1d0, dfar_zhi = 1d0
  Logical        :: far_ready = .False.

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


  !> Refresh the MUSCL far planes of C_ (device-resident in GPU builds; only four small planes go through the host)
  Subroutine update_far_planes(C_)

    Real(Int64), Dimension(nxg, nyg, nzg), Intent(In) :: C_

    Logical        :: is_first, is_last, per_x, per_z
    Integer(Int32) :: up, down, partner, dst_dn, dst_up, src_dn, src_up, p, nplane
    Integer(Int32) :: kg2, ig2

    per_x = ( x_bc_type == 0 )
    per_z = ( z_bc_type == 0 )

    If ( .Not. far_ready ) Then
       Allocate( far_xlo(nyg,nzg), far_xhi(nyg,nzg), far_zlo(nxg,nyg), far_zhi(nxg,nyg) )
       Allocate( fsend_xdn(nyg,nzg), fsend_xup(nyg,nzg), fsend_zdn(nxg,nyg), fsend_zup(nxg,nyg) )
       far_xlo = 0d0;  far_xhi = 0d0;  far_zlo = 0d0;  far_zhi = 0d0
       fsend_xdn = 0d0;  fsend_xup = 0d0;  fsend_zdn = 0d0;  fsend_zup = 0d0
       !$acc enter data copyin(far_xlo,far_xhi,far_zlo,far_zhi,fsend_xdn,fsend_xup,fsend_zdn,fsend_zup)

       ! availability and centre spacings, from the global grid (cells wrap as F(1)=F(n-2), F(n)=F(3), F(n+1)=F(4))
       Call z_halo_neighbors(up, down)
       have_far_zhi = ( up   /= MPI_PROC_NULL ) .Or. per_z
       have_far_zlo = ( down /= MPI_PROC_NULL ) .Or. per_z
       kg2 = kg2_global(myid)
       If ( up /= MPI_PROC_NULL ) Then
          dfar_zhi = zg_global(kg2+1) - zg_global(kg2)
       Else
          dfar_zhi = zg_global(4) - zg_global(3)
       End If
       If ( down /= MPI_PROC_NULL ) Then
          dfar_zlo = zg_global(kg1_global(myid)) - zg_global(kg1_global(myid)-1)
       Else
          dfar_zlo = zg_global(nzg_global-2) - zg_global(nzg_global-3)
       End If

       Call x_halo_neighbors(up, down)
       have_far_xhi = ( up   /= MPI_PROC_NULL ) .Or. per_x
       have_far_xlo = ( down /= MPI_PROC_NULL ) .Or. per_x
       ig2 = ig2_global(myid)
       If ( up /= MPI_PROC_NULL ) Then
          dfar_xhi = xg_global(ig2+1) - xg_global(ig2)
       Else
          dfar_xhi = xg_global(4) - xg_global(3)
       End If
       If ( down /= MPI_PROC_NULL ) Then
          dfar_xlo = xg_global(ig1_global(myid)) - xg_global(ig1_global(myid)-1)
       Else
          dfar_xlo = xg_global(nxg_global-2) - xg_global(nxg_global-3)
       End If
       far_ready = .True.
    End If

    !-- z ------------------------------------------------------------------
    Call z_halo_neighbors(up, down)
    Call z_periodic_partner(is_first, is_last, partner)
    If ( is_first .And. is_last ) Then
       If ( per_z ) Then
          !$acc kernels present(C_,far_zlo,far_zhi)
          far_zhi = C_(:,:,4)
          far_zlo = C_(:,:,nzg-3)
          !$acc end kernels
       End If
    Else
       dst_dn = down;  src_up = up;  dst_up = up;  src_dn = down
       If ( per_z .And. is_first ) Then
          dst_dn = partner;  src_dn = partner
       End If
       If ( per_z .And. is_last ) Then
          dst_up = partner;  src_up = partner
       End If
       ! plane sent to the lower neighbour becomes its high far plane; the periodic wrap sends plane 4 instead of 3
       p = 3;  If ( is_first ) p = 4
       !$acc kernels present(C_,fsend_zdn)
       fsend_zdn = C_(:,:,p)
       !$acc end kernels
       p = nzg-2;  If ( is_last ) p = nzg-3
       !$acc kernels present(C_,fsend_zup)
       fsend_zup = C_(:,:,p)
       !$acc end kernels
       !$acc update host(fsend_zdn,fsend_zup)
       nplane = nxg*nyg
       Call Mpi_sendrecv( fsend_zdn, nplane, Mpi_real8, dst_dn, 31, far_zhi, nplane, Mpi_real8, src_up, 31, &
                          MPI_COMM_WORLD, istat, ierr )
       Call Mpi_sendrecv( fsend_zup, nplane, Mpi_real8, dst_up, 32, far_zlo, nplane, Mpi_real8, src_dn, 32, &
                          MPI_COMM_WORLD, istat, ierr )
       !$acc update device(far_zlo,far_zhi)
    End If

    !-- x ------------------------------------------------------------------
    Call x_halo_neighbors(up, down)
    Call x_periodic_partner(is_first, is_last, partner)
    If ( is_first .And. is_last ) Then
       If ( per_x ) Then
          !$acc kernels present(C_,far_xlo,far_xhi)
          far_xhi = C_(4,:,:)
          far_xlo = C_(nxg-3,:,:)
          !$acc end kernels
       End If
    Else
       dst_dn = down;  src_up = up;  dst_up = up;  src_dn = down
       If ( per_x .And. is_first ) Then
          dst_dn = partner;  src_dn = partner
       End If
       If ( per_x .And. is_last ) Then
          dst_up = partner;  src_up = partner
       End If
       p = 3;  If ( is_first ) p = 4
       !$acc kernels present(C_,fsend_xdn)
       fsend_xdn = C_(p,:,:)
       !$acc end kernels
       p = nxg-2;  If ( is_last ) p = nxg-3
       !$acc kernels present(C_,fsend_xup)
       fsend_xup = C_(p,:,:)
       !$acc end kernels
       !$acc update host(fsend_xdn,fsend_xup)
       nplane = nyg*nzg
       Call Mpi_sendrecv( fsend_xdn, nplane, Mpi_real8, dst_dn, 33, far_xhi, nplane, Mpi_real8, src_up, 33, &
                          MPI_COMM_WORLD, istat, ierr )
       Call Mpi_sendrecv( fsend_xup, nplane, Mpi_real8, dst_up, 34, far_xlo, nplane, Mpi_real8, src_dn, 34, &
                          MPI_COMM_WORLD, istat, ierr )
       !$acc update device(far_xlo,far_xhi)
    End If

  End Subroutine update_far_planes


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
    Real   (Int64) :: c2, ddf                       ! far upwind cell value and its centre spacing
    Logical        :: hxlo, hxhi, hzlo, hzhi         ! far planes available (local copies for the device kernel)
    Real   (Int64) :: dxlo, dxhi, dzlo, dzhi

    Call update_far_planes(C_)
    hxlo = have_far_xlo;  hxhi = have_far_xhi;  hzlo = have_far_zlo;  hzhi = have_far_zhi
    dxlo = dfar_xlo;      dxhi = dfar_xhi;      dzlo = dfar_zlo;      dzhi = dfar_zhi

    !$acc parallel loop collapse(3) present(C_,U_,V_,W_,Fc_,nu_t,phi,x,xg,y,yg,z,zg,weight_y_0,weight_y_1,weight_z_0,weight_z_1,far_xlo,far_xhi,far_zlo,far_zhi)
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
                ! upwind cell = i+1; the far-upwind cell at i=nxg-1 is the far plane (first order where there is none)
                If ( i == nxg-1 .And. .Not. hxhi ) Then
                   C_hi = C_(i+1,j,k)
                Else
                   If ( i == nxg-1 ) Then
                      c2 = far_xhi(j,k);  ddf = dxhi
                   Else
                      c2 = C_(i+2,j,k);   ddf = xg(i+2) - xg(i+1)
                   End If
                   gf  = ( c2 - C_(i+1,j,k) ) / ddf
                   gb  = ( C_(i+1,j,k) - C_(i,j,k)   ) / ( xg(i+1) - xg(i)   )
                   slp = vanleer_slope(gf, gb)
                   C_hi = C_(i+1,j,k) + slp * ( x(i) - xg(i+1) )
                End If
             End If

             ! -- low face (i-1/2) at x(i-1): U_(i-1,j,k)
             uf = U_(i-1,j,k)
             If ( uf >= 0d0 ) Then
                ! upwind cell = i-1; the far-upwind cell at i=2 is the far plane (first order where there is none)
                If ( i == 2 .And. .Not. hxlo ) Then
                   C_lo = C_(i-1,j,k)
                Else
                   If ( i == 2 ) Then
                      c2 = far_xlo(j,k);  ddf = dxlo
                   Else
                      c2 = C_(i-2,j,k);   ddf = xg(i-1) - xg(i-2)
                   End If
                   gf  = ( C_(i,j,k)   - C_(i-1,j,k) ) / ( xg(i)   - xg(i-1) )
                   gb  = ( C_(i-1,j,k) - c2 ) / ddf
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
                ! upwind cell = k+1; the far-upwind cell at k=nzg-1 is the far plane (first order where there is none)
                If ( k == nzg-1 .And. .Not. hzhi ) Then
                   C_hi = C_(i,j,k+1)
                Else
                   If ( k == nzg-1 ) Then
                      c2 = far_zhi(i,j);  ddf = dzhi
                   Else
                      c2 = C_(i,j,k+2);   ddf = zg(k+2) - zg(k+1)
                   End If
                   gf  = ( c2 - C_(i,j,k+1) ) / ddf
                   gb  = ( C_(i,j,k+1) - C_(i,j,k)   ) / ( zg(k+1) - zg(k)   )
                   slp = vanleer_slope(gf, gb)
                   C_hi = C_(i,j,k+1) + slp * ( z(k) - zg(k+1) )
                End If
             End If

             ! -- low face (k-1/2) at z(k-1): W_(i,j,k-1)
             wf = W_(i,j,k-1)
             If ( wf >= 0d0 ) Then
                ! upwind cell = k-1; the far-upwind cell at k=2 is the far plane (first order where there is none)
                If ( k == 2 .And. .Not. hzlo ) Then
                   C_lo = C_(i,j,k-1)
                Else
                   If ( k == 2 ) Then
                      c2 = far_zlo(i,j);  ddf = dzlo
                   Else
                      c2 = C_(i,j,k-2);   ddf = zg(k-1) - zg(k-2)
                   End If
                   gf  = ( C_(i,j,k)   - C_(i,j,k-1) ) / ( zg(k)   - zg(k-1) )
                   gb  = ( C_(i,j,k-1) - c2 ) / ddf
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
             ! kappa at high x-face (x uniform: plain average; y, z use the interpolation weights): average of nu_t at (i,j,k) and (i+1,j,k)
             kappa_hi = kappa_mol + 0.5d0*(nu_t(i,j,k) + nu_t(i+1,j,k))*kappa_t_inv
             kappa_lo = kappa_mol + 0.5d0*(nu_t(i,j,k) + nu_t(i-1,j,k))*kappa_t_inv
             diff_x = ( kappa_hi*(C_(i+1,j,k) - C_(i,j,k))/(xg(i+1)-xg(i)) &
                      - kappa_lo*(C_(i,j,k) - C_(i-1,j,k))/(xg(i)-xg(i-1)) ) / dx_f

             !-------- y-diffusion -----------------------------------------
             kappa_hi = kappa_mol + ( weight_y_0(j  )*nu_t(i,j,k) + weight_y_1(j  )*nu_t(i,j+1,k) )*kappa_t_inv
             kappa_lo = kappa_mol + ( weight_y_0(j-1)*nu_t(i,j-1,k) + weight_y_1(j-1)*nu_t(i,j,k) )*kappa_t_inv
             diff_y = ( kappa_hi*(C_(i,j+1,k) - C_(i,j,k))/(yg(j+1)-yg(j)) &
                      - kappa_lo*(C_(i,j,k) - C_(i,j-1,k))/(yg(j)-yg(j-1)) ) / dy_f

             !-------- z-diffusion -----------------------------------------
             kappa_hi = kappa_mol + ( weight_z_0(k  )*nu_t(i,j,k) + weight_z_1(k  )*nu_t(i,j,k+1) )*kappa_t_inv
             kappa_lo = kappa_mol + ( weight_z_0(k-1)*nu_t(i,j,k-1) + weight_z_1(k-1)*nu_t(i,j,k) )*kappa_t_inv
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

    ! x direction: seam planes from the x-neighbour rank (p_row>1), then periodic wrap between the first and last
    ! x rank (a local wrap within each slab is only right when x is not split), or Dirichlet inflow / convective
    ! outflow. The exchange/wrap routines act on device memory in GPU builds, hence the update pair.
#ifdef GPU_POISSON
    !$acc update device(C_)
#endif
    Call update_ghost_interior_planes_x(C_, 4)
    If ( x_bc_type == 0 ) Call apply_periodic_bc_x(C_, 4)
#ifdef GPU_POISSON
    !$acc update host(C_)
#endif
    If ( x_bc_type /= 0 ) Then
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
