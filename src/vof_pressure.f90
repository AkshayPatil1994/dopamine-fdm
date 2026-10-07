!> Variable-density pressure operator of the two-fluid solver: face coefficients beta = 1/rho_face from the VOF field, the
!  operator A phi = div(beta grad phi) matching the discrete Laplacian of the fast Poisson solver when beta is constant, and a
!  preconditioned CG whose preconditioner is that fast solver (solve_poisson_equation) with the lightest-fluid coefficient
!  (tests/vof/spike_pressure_pcg.f90: 17/28/39 iterations to 1e-6 at density ratios 10/100/1000, independent of the grid).
!
!  Conventions follow the solver: cell-centred vectors have the shape (nxg,nyg,nzg) with ghost layers 1 and n; face coefficients
!  bu(1:nx,:,:), bv(:,1:ny,:), bw(:,:,1:nz) sit on the same faces as U, V, W. The periodic duplicate cell of the last rank and the
!  ghost cells carry zero weight in the inner products, so the CG sees exactly the distinct unknowns of the Poisson problem.
Module vof_pressure

  Use iso_fortran_env, Only : Int32, Int64
  Use global
  Use mpi
  Use decomp, Only : x_periodic_partner, z_periodic_partner
  Use boundary_conditions, Only : update_ghost_interior_planes_x
  Use scalar_transport, Only : finish_scalar_halos
  Use projection, Only : solve_poisson_equation, pois_layered, pois_bf, pois_bh, poisson_layered_supported
  Use vof_plic
  Use vof_advect, Only : vof_reconstruct, vof_mx, vof_my, vof_mz, vof_al

  Implicit None

  Real(Int64), Allocatable, Dimension(:,:,:) :: vp_rho, vp_bu, vp_bv, vp_bw, vp_w, vp_rfu, vp_rfv, vp_rfw, vp_rau, vp_rav, vp_raw
  Real(Int64), Allocatable, Dimension(:) :: vp_hy
  Real(Int64), Allocatable, Dimension(:,:,:) :: vp_r, vp_z, vp_d, vp_ap
  Real(Int64) :: vp_beta0 = 1d0, vp_wsum = 1d0
  Integer(Int32) :: vp_iters_last = 0
  Integer(Int64) :: vp_its_total = 0   ! PCG iterations since the start of the run
  Real(Int64) :: vp_res_last = 0d0

Contains

  Subroutine vp_init

    Integer(Int32) :: i, j, k, ihi, jhi, khi
    Logical :: is_first, is_last
    Integer(Int32) :: partner
    Real(Int64) :: wl

    Allocate( vp_rho(nxg,nyg,nzg), vp_bu(nx,nyg,nzg), vp_bv(nxg,ny,nzg), vp_bw(nxg,nyg,nz), vp_w(nxg,nyg,nzg) )
    Allocate( vp_r(nxg,nyg,nzg), vp_z(nxg,nyg,nzg), vp_d(nxg,nyg,nzg), vp_ap(nxg,nyg,nzg) )
    Allocate( vp_rfu(nx,nyg,nzg), vp_rfv(nxg,ny,nzg), vp_rfw(nxg,nyg,nz), vp_hy(nyg) )
    Allocate( vp_rau(nx,nyg,nzg), vp_rav(nxg,ny,nzg), vp_raw(nxg,nyg,nz) )
    vp_rau = vof_rho_g;  vp_rav = vof_rho_g;  vp_raw = vof_rho_g
    vp_rfu = vof_rho_g;  vp_rfv = vof_rho_g;  vp_rfw = vof_rho_g
    Allocate( pois_bf(nyg), pois_bh(nyg) )
    pois_bf = 1d0/vof_rho_g;  pois_bh = 1d0/vof_rho_g
    Do j = 2, nyg-1
       vp_hy(j) = y(j) - y(j-1)
    End Do
    vp_hy(1) = vp_hy(2);  vp_hy(nyg) = vp_hy(nyg-1)
    vp_rho = vof_rho_g;  vp_bu = 1d0/vof_rho_g;  vp_bv = 1d0/vof_rho_g;  vp_bw = 1d0/vof_rho_g
    vp_r = 0d0;  vp_z = 0d0;  vp_d = 0d0;  vp_ap = 0d0
    vp_beta0 = 1d0/Min(vof_rho_l, vof_rho_g)

    ihi = nxg-1;  jhi = nyg-1;  khi = nzg-1
    If ( x_bc_type == 0 ) Then
       Call x_periodic_partner(is_first, is_last, partner)
       If ( is_last ) ihi = nxg-2
    End If
    If ( z_bc_type == 0 ) Then
       Call z_periodic_partner(is_first, is_last, partner)
       If ( is_last ) khi = nzg-2
    End If
    If ( y_bc_type == 0 ) jhi = nyg-2

    vp_w = 0d0
    wl = 0d0
    Do k = 2, khi
       Do j = 2, jhi
          Do i = 2, ihi
             vp_w(i,j,k) = dx*(y(j) - y(j-1))*(z(k) - z(k-1))
             wl = wl + vp_w(i,j,k)
          End Do
       End Do
    End Do
    Call MPI_Allreduce(wl, vp_wsum, 1, MPI_real8, MPI_SUM, MPI_COMM_WORLD, ierr)

  End Subroutine vp_init


  !> Cell densities (ghosts included), face densities (arithmetic in x and z, volume-weighted across the stretched y cells so
  !  that rho_face * u is the mass-consistent momentum of the staggered control volume) and beta = 1/rho_face, from the padded C
  Subroutine vp_set_density(Cp)

    Real(Int64), Intent(In) :: Cp(0:nxg+1,0:nyg+1,0:nzg+1)
    Integer(Int32) :: i, j, k
    Real(Int64) :: dr, hl, hu

    Do k = 1, nzg
       Do j = 1, nyg
          Do i = 1, nxg
             vp_rho(i,j,k) = vof_rho_g + (vof_rho_l - vof_rho_g)*Cp(i,j,k)
          End Do
       End Do
    End Do
    Do k = 1, nzg
       Do j = 1, nyg
          Do i = 1, Min(nx, nxg-1)
             vp_rau(i,j,k) = 0.5d0*( vp_rho(i,j,k) + vp_rho(i+1,j,k) )
          End Do
       End Do
    End Do
    ! on a rank with a +x neighbour the last face is the duplicate of the neighbour's first one: its far cell is the pad layer
    If ( nx == nxg ) vp_rau(nx,:,:) = 0.5d0*( vp_rho(nx,:,:) + vof_rho_g + (vof_rho_l - vof_rho_g)*Cp(nxg+1,1:nyg,1:nzg) )
    Do k = 1, nzg
       Do j = 1, ny
          Do i = 1, nxg
             vp_rav(i,j,k) = ( vp_rho(i,j,k)*vp_hy(j) + vp_rho(i,j+1,k)*vp_hy(j+1) )/( vp_hy(j) + vp_hy(j+1) )
          End Do
       End Do
    End Do
    Do k = 1, nz
       Do j = 1, nyg
          Do i = 1, nxg
             vp_raw(i,j,k) = 0.5d0*( vp_rho(i,j,k) + vp_rho(i,j,Min(k+1,nzg)) )
          End Do
       End Do
    End Do
    If ( vof_geo_density >= 1 ) Then
       ! geometric staggered-cell density (Fuster/Arrufat): the liquid in the half cells that make up the control volume of the
       ! face is taken from each cell's reconstructed plane, so that a flat interface is hydrostatically exact whatever C is
       Call vof_reconstruct(Cp, nxg, nyg, nzg, vof_normal_scheme)
       dr = vof_rho_l - vof_rho_g
       Do k = 1, nzg
          Do j = 1, nyg
             Do i = 1, Min(nx, nxg-1)
                vp_rfu(i,j,k) = vof_rho_g + dr*0.5d0*( half_frac(Cp(i,j,k), i,j,k, 1, .True.) &
                                                     + half_frac(Cp(i+1,j,k), i+1,j,k, 1, .False.) )
             End Do
          End Do
       End Do
       If ( nx == nxg ) vp_rfu(nx,:,:) = vp_rau(nx,:,:)
       Do k = 1, nzg
          Do j = 1, ny
             hl = vp_hy(j);  hu = vp_hy(j+1)
             Do i = 1, nxg
                vp_rfv(i,j,k) = vof_rho_g + dr*( hl*half_frac(Cp(i,j,k), i,j,k, 2, .True.) &
                                               + hu*half_frac(Cp(i,j+1,k), i,j+1,k, 2, .False.) )/( hl + hu )
             End Do
          End Do
       End Do
       Do k = 1, nz
          Do j = 1, nyg
             Do i = 1, nxg
                vp_rfw(i,j,k) = vof_rho_g + dr*0.5d0*( half_frac(Cp(i,j,k), i,j,k, 3, .True.) &
                                                     + half_frac(Cp(i,j,Min(k+1,nzg)), i,j,Min(k+1,nzg), 3, .False.) )
             End Do
          End Do
       End Do
    Else
       vp_rfu = vp_rau;  vp_rfv = vp_rav;  vp_rfw = vp_raw
    End If
    vp_bu = 1d0/vp_rfu;  vp_bv = 1d0/vp_rfv;  vp_bw = 1d0/vp_rfw
    If ( vof_layered_precond >= 1 .And. poisson_layered_supported() ) Call layer_coefficients

  End Subroutine vp_set_density


  !> Row-wise effective coefficients of the layered preconditioner: arithmetic horizontal-mean beta in the plane (parallel
  !  conduction), harmonic mean of beta across the rows (series), weighted by the cell volumes
  Subroutine layer_coefficients

    Integer(Int32) :: i, j, k
    Real(Int64) :: acc(3,nyg), glb(3,nyg), wk

    acc = 0d0
    Do k = 1, nzg
       Do j = 1, nyg
          Do i = 1, nxg
             wk = vp_w(i,j,k)
             If ( wk == 0d0 ) Cycle
             acc(1,j) = acc(1,j) + wk
             acc(2,j) = acc(2,j) + wk*0.5d0*( vp_bu(Min(i,nx),j,k) + vp_bw(i,j,Min(k,nz)) )
             If ( j <= ny ) acc(3,j) = acc(3,j) + wk*vp_rfv(i,j,k)
          End Do
       End Do
    End Do
    Call MPI_Allreduce(acc, glb, 3*nyg, MPI_real8, MPI_SUM, MPI_COMM_WORLD, ierr)
    Do j = 1, nyg
       If ( glb(1,j) > 0d0 ) Then
          pois_bh(j) = glb(2,j)/glb(1,j)
          If ( j <= ny ) pois_bf(j) = glb(1,j)/glb(3,j)
       End If
    End Do

  End Subroutine layer_coefficients


  !> Liquid fraction of the upper (upper=.True.) or lower half of the cell (i,j,k) along direction dir, from its plane
  Function half_frac(c, i, j, k, dir, upper) Result(f)

    Real(Int64),    Intent(In) :: c
    Integer(Int32), Intent(In) :: i, j, k, dir
    Logical,        Intent(In) :: upper
    Real(Int64) :: f, xlo(3), w(3)

    If ( c <= vof_eps ) Then
       f = 0d0
    Else If ( c >= 1d0 - vof_eps ) Then
       f = 1d0
    Else
       xlo = -0.5d0;  w = 1d0
       w(dir) = 0.5d0
       If ( upper ) xlo(dir) = 0d0
       f = plic_subbox_fraction( vof_mx(i,j,k), vof_my(i,j,k), vof_mz(i,j,k), vof_al(i,j,k), &
                                 xlo(1), xlo(2), xlo(3), w(1), w(2), w(3) )
    End If

  End Function half_frac


  !> Ghost layers of a cell-centred field: rank seams and periodic wraps, zero gradient at y walls (periodic wrap for y_bc_type=0),
  !  zero gradient at the x inlet and, at the x outlet, mirror (antisym=.False.) or the Dirichlet-zero reflection of the pressure
  !  problem (antisym=.True.)
  Subroutine vp_halo(F, antisym)

    Real(Int64), Intent(InOut) :: F(nxg,nyg,nzg)
    Logical,     Intent(In)    :: antisym

    Logical :: is_first, is_last
    Integer(Int32) :: partner

    Call update_ghost_interior_planes_x(F, 4)
    If ( x_bc_type == 1 ) Then
       Call x_periodic_partner(is_first, is_last, partner)
       If ( is_first ) F(1,:,:) = F(2,:,:)
       If ( is_last ) Then
          If ( antisym ) Then
             F(nxg,:,:) = -F(nxg-1,:,:)
          Else
             F(nxg,:,:) = F(nxg-1,:,:)
          End If
       End If
    End If
    Call finish_scalar_halos(F)
    If ( y_bc_type == 0 ) Then
       F(:,1,:)       = F(:,nyg-2,:)
       F(:,nyg-1,:)   = F(:,2,:)
       F(:,nyg,:)     = F(:,3,:)
    Else
       F(:,1,:)   = F(:,2,:)
       F(:,nyg,:) = F(:,nyg-1,:)
    End If

  End Subroutine vp_halo


  !> out = div(beta grad phi) at the interior cells; phi needs valid ghost layers (vp_halo with antisym=.True.)
  Subroutine vp_apply_A(phi, out)

    Real(Int64), Intent(In)  :: phi(nxg,nyg,nzg)
    Real(Int64), Intent(Out) :: out(nxg,nyg,nzg)

    Integer(Int32) :: i, j, k
    Real(Int64) :: inv_dx2, hy, hz, ihy_hi, ihy_lo, ihz_hi, ihz_lo, fxh, fxl, fyh, fyl, fzh, fzl

    inv_dx2 = 1d0/(dx*dx)
    out = 0d0
    Do k = 2, nzg-1
       hz = z(k) - z(k-1)
       ihz_hi = 1d0/(zg(k+1) - zg(k));  ihz_lo = 1d0/(zg(k) - zg(k-1))
       Do j = 2, nyg-1
          hy = y(j) - y(j-1)
          ihy_hi = 1d0/(yg(j+1) - yg(j));  ihy_lo = 1d0/(yg(j) - yg(j-1))
          Do i = 2, nxg-1
             fxh = vp_bu(i,j,k)*( phi(i+1,j,k) - phi(i,j,k) )
             fxl = vp_bu(i-1,j,k)*( phi(i,j,k) - phi(i-1,j,k) )
             fyh = vp_bv(i,j,k)*( phi(i,j+1,k) - phi(i,j,k) )*ihy_hi
             fyl = vp_bv(i,j-1,k)*( phi(i,j,k) - phi(i,j-1,k) )*ihy_lo
             fzh = vp_bw(i,j,k)*( phi(i,j,k+1) - phi(i,j,k) )*ihz_hi
             fzl = vp_bw(i,j,k-1)*( phi(i,j,k) - phi(i,j,k-1) )*ihz_lo
             out(i,j,k) = (fxh - fxl)*inv_dx2 + (fyh - fyl)/hy + (fzh - fzl)/hz
          End Do
       End Do
    End Do

  End Subroutine vp_apply_A


  !> Face fields beta grad(phi): gu(1:nx), gv(1:ny), gw(1:nz) (ghost layers of phi must be valid)
  Subroutine vp_grad(phi, gu, gv, gw)

    Real(Int64), Intent(In)  :: phi(nxg,nyg,nzg)
    Real(Int64), Intent(Out) :: gu(nx,nyg,nzg), gv(nxg,ny,nzg), gw(nxg,nyg,nz)

    Integer(Int32) :: i, j, k
    Real(Int64) :: inv_dx

    inv_dx = 1d0/dx
    gu = 0d0;  gv = 0d0;  gw = 0d0
    Do k = 1, nzg
       Do j = 1, nyg
          Do i = 1, Min(nx, nxg-1)
             gu(i,j,k) = vp_bu(i,j,k)*( phi(i+1,j,k) - phi(i,j,k) )*inv_dx
          End Do
       End Do
    End Do
    Do k = 1, nzg
       Do j = 1, ny
          Do i = 1, nxg
             gv(i,j,k) = vp_bv(i,j,k)*( phi(i,j+1,k) - phi(i,j,k) )/( yg(j+1) - yg(j) )
          End Do
       End Do
    End Do
    Do k = 1, Min(nz, nzg-1)
       Do j = 1, nyg
          Do i = 1, nxg
             gw(i,j,k) = vp_bw(i,j,k)*( phi(i,j,k+1) - phi(i,j,k) )/( zg(k+1) - zg(k) )
          End Do
       End Do
    End Do

  End Subroutine vp_grad


  !> d = div(F) at the interior cells for face fields with the shapes of U, V, W
  Subroutine vp_div(Fu, Fv, Fw, d)

    Real(Int64), Intent(In)  :: Fu(nx,nyg,nzg), Fv(nxg,ny,nzg), Fw(nxg,nyg,nz)
    Real(Int64), Intent(Out) :: d(nxg,nyg,nzg)

    Integer(Int32) :: i, j, k
    Real(Int64) :: inv_dx, inv_hy, inv_hz

    inv_dx = 1d0/dx
    d = 0d0
    Do k = 2, nzg-1
       inv_hz = 1d0/(z(k) - z(k-1))
       Do j = 2, nyg-1
          inv_hy = 1d0/(y(j) - y(j-1))
          Do i = 2, nxg-1
             d(i,j,k) = ( Fu(i,j,k) - Fu(i-1,j,k) )*inv_dx + ( Fv(i,j,k) - Fv(i,j-1,k) )*inv_hy &
                      + ( Fw(i,j,k) - Fw(i,j,k-1) )*inv_hz
          End Do
       End Do
    End Do

  End Subroutine vp_div


  Function vp_dot(a, b) Result(s)

    Real(Int64), Intent(In) :: a(nxg,nyg,nzg), b(nxg,nyg,nzg)
    Real(Int64) :: s, sl

    sl = Sum( vp_w*a*b )
    Call MPI_Allreduce(sl, s, 1, MPI_real8, MPI_SUM, MPI_COMM_WORLD, ierr)

  End Function vp_dot


  !> max|a| over the distinct unknowns (zero-weight ghost and duplicate cells excluded)
  Function vp_maxabs(a) Result(s)

    Real(Int64), Intent(In) :: a(nxg,nyg,nzg)
    Real(Int64) :: s, sl

    sl = MaxVal( Abs(a), mask=( vp_w > 0d0 ) )
    Call MPI_Allreduce(sl, s, 1, MPI_real8, MPI_MAX, MPI_COMM_WORLD, ierr)

  End Function vp_maxabs


  Subroutine vp_remove_mean(a)

    Real(Int64), Intent(InOut) :: a(nxg,nyg,nzg)
    Real(Int64) :: sl, s

    ! the pressure is pinned by the Dirichlet outlet when x is not periodic, so the operator is not singular: a mean removal
    ! would perturb the system and stall the PCG at a finite residual
    If ( x_bc_type == 1 ) Return
    sl = Sum( vp_w*a )
    Call MPI_Allreduce(sl, s, 1, MPI_real8, MPI_SUM, MPI_COMM_WORLD, ierr)
    a = a - s/vp_wsum

  End Subroutine vp_remove_mean


  !> z = M^-1 r: constant-coefficient fast solve with beta0, zero weighted mean (keeps the preconditioner symmetric)
  Subroutine vp_precond(r, z)

    Real(Int64), Intent(In)  :: r(nxg,nyg,nzg)
    Real(Int64), Intent(Out) :: z(nxg,nyg,nzg)

    z = 0d0
    If ( vof_layered_precond >= 1 .And. poisson_layered_supported() ) Then
       rhs_p(2:nxg,2:nyg-1,2:nzg) = r(2:nxg,2:nyg-1,2:nzg)
       pois_layered = .True.
       Call solve_poisson_equation(skip_p_save=.True.)
       pois_layered = .False.
    Else
       rhs_p(2:nxg,2:nyg-1,2:nzg) = r(2:nxg,2:nyg-1,2:nzg)/vp_beta0
       Call solve_poisson_equation(skip_p_save=.True.)
    End If
    z(2:nxg,2:nyg-1,2:nzg) = rhs_p(2:nxg,2:nyg-1,2:nzg)
    Call vp_remove_mean(z)

  End Subroutine vp_precond


  !> Up to nit PCG iterations on A x = f from x = 0 (f need not have zero mean: it is removed), stopping early once the residual
  !  falls below tol times its initial norm (tol = 0: always nit iterations). With vof_div_tol > 0 and rscale given, the stop is
  !  instead rscale*max|r| < vof_div_tol: r is the divergence left after the correction, rscale the time scale of f.
  !  Sets vp_iters_last and the final relative residual.
  Subroutine vp_pcg(f, nit, tol, x, rscale)

    Real(Int64),    Intent(In)  :: f(nxg,nyg,nzg), tol
    Integer(Int32), Intent(In)  :: nit
    Real(Int64),    Intent(Out) :: x(nxg,nyg,nzg)
    Real(Int64),    Intent(In), Optional :: rscale

    Integer(Int32) :: it
    Real(Int64) :: rz, rzn, alpha, r0, rn
    Logical :: abs_stop

    abs_stop = ( vof_div_tol > 0d0 .And. Present(rscale) )
    x = 0d0
    vp_iters_last = 0
    vp_res_last = 0d0
    vp_r = f
    Call vp_remove_mean(vp_r)
    r0 = Sqrt(vp_dot(vp_r, vp_r))
    If ( r0 == 0d0 ) Return
    If ( abs_stop ) Then
       If ( rscale*vp_maxabs(vp_r) < vof_div_tol ) Return
    End If
    Call vp_precond(vp_r, vp_z)
    vp_d = vp_z
    rz = vp_dot(vp_r, vp_z)
    Do it = 1, nit
       Call vp_halo(vp_d, .True.)
       Call vp_apply_A(vp_d, vp_ap)
       alpha = rz/vp_dot(vp_d, vp_ap)
       x = x + alpha*vp_d
       vp_r = vp_r - alpha*vp_ap
       Call vp_remove_mean(vp_r)
       vp_iters_last = it
       vp_its_total = vp_its_total + 1
       rn = Sqrt(vp_dot(vp_r, vp_r))
       vp_res_last = rn/r0
       If ( it == nit ) Exit
       If ( abs_stop ) Then
          If ( rscale*vp_maxabs(vp_r) < vof_div_tol ) Exit
       Else If ( rn < tol*r0 ) Then
          Exit
       End If
       Call vp_precond(vp_r, vp_z)
       rzn = vp_dot(vp_r, vp_z)
       vp_d = vp_z + (rzn/rz)*vp_d
       rz = rzn
    End Do

  End Subroutine vp_pcg


  !> Development check, run from vof_init when vof_selftest=1: PCG iteration counts to 1e-6 and 1e-10 for a rough zero-mean right-hand
  !  side on the current density field, and the symmetry of A (a.Ab - b.Aa)
  Subroutine vp_selftest

    Real(Int64), Allocatable :: f(:,:,:), x(:,:,:), a(:,:,:), b(:,:,:), Aa(:,:,:), Ab(:,:,:)
    Integer(Int32) :: i, j, k
    Real(Int64) :: sab, sba, tol_list(2)
    Integer(Int32) :: it

    Allocate( f(nxg,nyg,nzg), x(nxg,nyg,nzg), a(nxg,nyg,nzg), b(nxg,nyg,nzg), Aa(nxg,nyg,nzg), Ab(nxg,nyg,nzg) )
    f = 0d0
    Do k = 2, nzg-1
       Do j = 2, nyg-1
          Do i = 2, nxg-1
             f(i,j,k) = Sin(6.283d0*xg(i) + 3.14d0*yg(j) + 4d0*zg(k)) + Cos(2d0*xg(i) - 9.42d0*yg(j) + 6d0*zg(k))
             a(i,j,k) = Sin(3d0*xg(i) + 2d0*yg(j)) + Cos(5d0*zg(k))
             b(i,j,k) = Cos(2d0*xg(i) - 3d0*yg(j) + zg(k))
          End Do
       End Do
    End Do
    f = f*vp_w/Max(vp_w, 1d-300)
    a = a*vp_w/Max(vp_w, 1d-300)
    b = b*vp_w/Max(vp_w, 1d-300)
    Call vp_halo(a, .True.);  Call vp_halo(b, .True.)
    Call vp_apply_A(a, Aa);  Call vp_apply_A(b, Ab)
    sab = vp_dot(a, Ab);  sba = vp_dot(b, Aa)
    If ( myid == 0 ) Write(*,'(A,3ES14.5)') '   vp_selftest symmetry a.Ab, b.Aa, rel diff = ', sab, sba, &
         Abs(sab - sba)/Max(Abs(sab), 1d-300)

    tol_list = (/ 1d-6, 1d-10 /)
    Do it = 1, 2
       Call vp_pcg(f, 200, tol_list(it), x)
       If ( myid == 0 ) Write(*,'(A,ES9.1,A,I4,A,ES10.2)') '   vp_selftest PCG tol', tol_list(it), ' iterations', &
            vp_iters_last, ' final residual', vp_res_last
       Call vp_remove_mean(x)
       sab = vp_dot(x, x)
       If ( myid == 0 ) Write(*,'(A,2ES20.12)') '   vp_selftest solution norm^2, f.x = ', sab, vp_dot(f, x)
    End Do
    Deallocate( f, x, a, b, Aa, Ab )

  End Subroutine vp_selftest

End Module vof_pressure
