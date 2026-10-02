!> Spike C: momentum form of a one-fluid VOF solver at large density ratio, in a serial 2-D closed-box mini-solver that combines the
!  validated pieces (PLIC sweeps, K=0 pressure predictor + PCG increment with the fast solver as preconditioner, still-water
!  hydrostatic subtraction). Per step (Strang): advect(dt/2) - project - RK3 forces with frozen rho - advect(dt/2) - project.
!    form B (conservative): rho*u is advected in the same sweeps as C with the sweeps' own mass fluxes; the Weymouth-Yue
!    c-tilde*div(u)
!       correction enters the momentum update with rho-tilde = rho_g + (rho_l-rho_g) c-tilde so mass and momentum stay consistent
!       (uniform u is reproduced exactly at any density ratio);
!    form A (velocity): u advected inside the RK3 stages with the skew-symmetric central operator, no density in the convection.
!
!  Results (N = 48, sloshing amplitude 2.4 cells, g such that T ~ 1, half-period error versus two-layer linear theory):
!   - Hydrostatic subtraction with the still-water profile keeps a flat interface inside a cell at rest to 1e-32 (T1).
!   - Gas-phase spurious velocities: at ratio 1000 the lateral hydrostatic pressure of sliver cells (a few % liquid in a gas cell)
!     acts on faces of 10-1000x gas density: |u_gas| ~ 0.1-0.2 after one dt, independent of h and of the amplitude, growing with the
!     ratio (clean at ratio 10). They are a physical-looking but non-convergent artifact of cell-averaged mixing, and they set the
!     time step: with a fixed dt the Courant number exceeds 1 and the period error appears to grow with refinement (13-20 %);
!     with Co <= 0.4 the errors are small and bounded. The step count then grows (349 vs 64 steps at ratio 1000, central).
!   - Conservative rho*u transport with the c-tilde-consistent update (form 1) is stable and accurate once Co is controlled:
!     central interpolation +0.9 % (ratio 100) and +0.5 % (1000); first-order upwind +2.1 % and damped/stalled at 1000;
!     van Leer +2.0 % and +6.6 %. The no-transport control (form 0) gives -1.0 % and -1.6 %.
!   - Density regularisation (snap thresholds 0.05-0.2) does not cure the gas currents.
Module spike_c_mod

  Use iso_fortran_env, Only : Int32, Int64
  Use vof_plic
  Use vof_normals
  Use vof_advect
  Use spike_proj_lib, Only : proj_init, n, h, pi, grad_beta, divergence, fast_solve, exact_solve, cg_iterations
  Implicit None

  Integer(Int32), Parameter :: nz = 5
  Integer(Int32) :: n1, n2, n3 = nz
  Real(Int64) :: rho_l, rho_g, gacc
  Real(Int64), Allocatable :: Cp(:,:,:), ux(:,:), uy(:,:), qx(:,:), qy(:,:), ppre(:,:), rhos(:), hh(:)
  Integer(Int32) :: mcg = 3
  Real(Int64) :: cfl_lim = 0d0      ! > 0: adaptive dt = min(T/64, cfl_lim h / max|u|)
  Integer(Int32) :: nsteps_used = 0
  Real(Int64) :: snap = 0d0         ! fractions below snap (above 1-snap) count as pure gas (liquid) in the dynamics' density
  Integer(Int32) :: mom_scheme = 0    ! 0 central, 1 first-order upwind, 2 van Leer limited upwind-biased
  Real(Int64), Parameter :: acoef(3,3) = Reshape( (/ 8d0/15d0, 1d0/4d0, 1d0/4d0, &
                                                     0d0,      5d0/12d0, 0d0,    &
                                                     0d0,      0d0,      3d0/4d0 /), (/3,3/) )

Contains

  Subroutine setup(nn, rl, rg, g)

    Integer(Int32), Intent(In) :: nn
    Real(Int64),    Intent(In) :: rl, rg, g

    Call proj_init(nn)
    n1 = n + 2;  n2 = n + 2
    rho_l = rl;  rho_g = rg;  gacc = g
    If ( Allocated(Cp) ) Deallocate( Cp, ux, uy, qx, qy, ppre, rhos, hh )
    Allocate( Cp(0:n1+1,0:n2+1,0:n3+1), ux(n-1,n), uy(n,n-1), qx(n-1,n), qy(n,n-1), ppre(n,n), rhos(n), hh(Max(n1,n3)) )
    hh = h
    Call vof_advect_init(n1, n2, n3)
    Cp = 0d0;  ux = 0d0;  uy = 0d0;  qx = 0d0;  qy = 0d0;  ppre = 0d0

  End Subroutine setup


  !> Face value of the transported velocity between the two middle samples (u_m1, u_0 | u_1, u_2 are the four consecutive samples
  !  around the cv face; the face sits between u_0 and u_1, sample order along the flux direction), given the sign of the mass flux
  Pure Function pick(um1, u0, u1, u2, mflux) Result(v)

    Real(Int64), Intent(In) :: um1, u0, u1, u2, mflux
    Real(Int64) :: v, r, phi, d

    ! the cv-face value for a flux through the point halfway between the samples u0 and u1
    Select Case(mom_scheme)
    Case(0)
       v = 0.5d0*( u0 + u1 )
    Case(1)
       v = Merge(u0, u1, mflux >= 0d0)
    Case Default
       If ( mflux >= 0d0 ) Then
          d = u1 - u0
          If ( Abs(d) < 1d-300 ) Then
             v = u0
          Else
             r = (u0 - um1)/d
             phi = (r + Abs(r))/(1d0 + Abs(r))
             v = u0 + 0.5d0*phi*d
          End If
       Else
          d = u0 - u1
          If ( Abs(d) < 1d-300 ) Then
             v = u1
          Else
             r = (u1 - u2)/d
             phi = (r + Abs(r))/(1d0 + Abs(r))
             v = u1 + 0.5d0*phi*d
          End If
       End If
    End Select

  End Function pick


  !> x-face array value with zero outside the interior faces (closed walls)
  Pure Function xf(arr, a, b) Result(v)

    Real(Int64), Intent(In) :: arr(n-1,n)
    Integer(Int32), Intent(In) :: a, b
    Real(Int64) :: v

    v = 0d0
    If ( a >= 1 .And. a <= n-1 .And. b >= 1 .And. b <= n ) v = arr(a,b)

  End Function xf


  Pure Function yf(arr, a, b) Result(v)

    Real(Int64), Intent(In) :: arr(n,n-1)
    Integer(Int32), Intent(In) :: a, b
    Real(Int64) :: v

    v = 0d0
    If ( a >= 1 .And. a <= n .And. b >= 1 .And. b <= n-1 ) v = arr(a,b)

  End Function yf


  Pure Function rho_of(c) Result(r)

    Real(Int64), Intent(In) :: c
    Real(Int64) :: r

    r = rho_g + (rho_l - rho_g)*Min(1d0, Max(0d0, (c - snap)/(1d0 - 2d0*snap)))

  End Function rho_of


  Subroutine rho_cells(rho)

    Real(Int64), Intent(Out) :: rho(n,n)
    Integer(Int32) :: a, b

    Do b = 1, n
       Do a = 1, n
          rho(a,b) = rho_of(Cp(a+1,b+1,3))
       End Do
    End Do

  End Subroutine rho_cells


  Subroutine face_rho(rho, rfx, rfy)

    Real(Int64), Intent(In)  :: rho(n,n)
    Real(Int64), Intent(Out) :: rfx(n-1,n), rfy(n,n-1)

    rfx = 0.5d0*( rho(1:n-1,:) + rho(2:n,:) )
    rfy = 0.5d0*( rho(:,1:n-1) + rho(:,2:n) )

  End Subroutine face_rho


  Subroutine fill_pad_box(C, m1, m2, m3)

    Integer(Int32), Intent(In)    :: m1, m2, m3
    Real(Int64),    Intent(InOut) :: C(0:m1+1,0:m2+1,0:m3+1)

    C(1,:,:) = C(2,:,:);        C(0,:,:) = C(3,:,:)
    C(m1,:,:) = C(m1-1,:,:);    C(m1+1,:,:) = C(m1-2,:,:)
    C(:,1,:) = C(:,2,:);        C(:,0,:) = C(:,3,:)
    C(:,m2,:) = C(:,m2-1,:);    C(:,m2+1,:) = C(:,m2-2,:)
    C(:,:,1) = C(:,:,m3-1);     C(:,:,m3) = C(:,:,2)
    C(:,:,0) = C(:,:,m3-2);     C(:,:,m3+1) = C(:,:,3)

  End Subroutine fill_pad_box


  !> Weymouth-Yue sweeps for C and the conservative momentum update; (ut,vt) is the div-free transporting velocity
  Subroutine advect_half(dtau, istep, ut, vt)

    Real(Int64),    Intent(In) :: dtau, ut(n-1,n), vt(n,n-1)
    Integer(Int32), Intent(In) :: istep
    Real(Int64) :: U3(n1-1,n2,n3), V3(n1,n2-1,n3), W3(n1,n2,n3-1)
    Real(Int64) :: Mx(0:n,n), My(n,0:n), rt(n,n), Dd(n,n), dq(n-1,n), dqy(n,n-1)
    Real(Int64) :: rho(n,n), rfx(n-1,n), rfy(n,n-1), vcv, cl, co, area
    Integer(Int32) :: isw, d, a, b, order(2), i, j, k
    Real(Int64), Allocatable :: uqx(:,:), uqy(:,:)

    vcv = h**3
    area = h*h
    U3 = 0d0;  V3 = 0d0;  W3 = 0d0
    Do b = 1, n
       Do a = 1, n-1
          U3(a+1,b+1,:) = ut(a,b)
       End Do
    End Do
    Do b = 1, n-1
       Do a = 1, n
          V3(a+1,b+1,:) = vt(a,b)
       End Do
    End Do
    Allocate( uqx(n-1,n), uqy(n,n-1) )
    Do k = 1, n3
       Do j = 1, n2
          Do i = 1, n1
             vof_cc(i,j,k) = Merge(1d0, 0d0, Cp(i,j,k) > 0.5d0)
          End Do
       End Do
    End Do
    If ( Mod(istep,2) == 0 ) Then
       order = (/ 1, 2 /)
    Else
       order = (/ 2, 1 /)
    End If

    Do isw = 1, 2
       d = order(isw)
       Call fill_pad_box(Cp, n1, n2, n3)
       Call vof_reconstruct(Cp, n1, n2, n3, VOF_SCHEME_YOUNGS)
       Call rho_cells(rho)
       Call face_rho(rho, rfx, rfy)
       uqx = qx/rfx;  uqy = qy/rfy
       Mx = 0d0;  My = 0d0
       If ( d == 1 ) Then
          Call vof_sweep(1, Cp, n1, n2, n3, U3, hh(1:n1), hh(1:n2), hh(1:n3), dtau, co, cl)
          Do b = 1, n
             Do a = 1, n-1
                Mx(a,b) = rho_g*ut(a,b)*dtau*area + (rho_l - rho_g)*vof_flux(a+1,b+1,3)
             End Do
          End Do
       Else
          Call vof_sweep(2, Cp, n1, n2, n3, V3, hh(1:n1), hh(1:n2), hh(1:n3), dtau, co, cl)
          Do b = 1, n-1
             Do a = 1, n
                My(a,b) = rho_g*vt(a,b)*dtau*area + (rho_l - rho_g)*vof_flux(a+1,b+1,3)
             End Do
          End Do
       End If
       ! rho-tilde * (volume-flux divergence) per cell, from the frozen c-tilde
       Do b = 1, n
          Do a = 1, n
             rt(a,b) = rho_g + (rho_l - rho_g)*vof_cc(a+1,b+1,3)
             If ( d == 1 ) Then
                Dd(a,b) = dtau*area*( xf(ut,a,b) - xf(ut,a-1,b) )
             Else
                Dd(a,b) = dtau*area*( yf(vt,a,b) - yf(vt,a,b-1) )
             End If
          End Do
       End Do
       dq = 0d0;  dqy = 0d0
       If ( d == 1 ) Then
          Do b = 1, n
             Do a = 1, n-1
                dq(a,b) = -( Qc_x(a+1,b) - Qc_x(a,b) ) + uqx(a,b)*0.5d0*( rt(a,b)*Dd(a,b) + rt(a+1,b)*Dd(a+1,b) )
             End Do
          End Do
          Do b = 1, n-1
             Do a = 1, n
                dqy(a,b) = -( Qcorner_x(a,b,1) - Qcorner_x(a,b,0) ) + uqy(a,b)*0.5d0*( rt(a,b)*Dd(a,b) + rt(a,b+1)*Dd(a,b+1) )
             End Do
          End Do
       Else
          Do b = 1, n-1
             Do a = 1, n
                dqy(a,b) = -( Qc_y(a,b+1) - Qc_y(a,b) ) + uqy(a,b)*0.5d0*( rt(a,b)*Dd(a,b) + rt(a,b+1)*Dd(a,b+1) )
             End Do
          End Do
          Do b = 1, n
             Do a = 1, n-1
                dq(a,b) = -( Qcorner_y(a,b,1) - Qcorner_y(a,b,0) ) + uqx(a,b)*0.5d0*( rt(a,b)*Dd(a,b) + rt(a+1,b)*Dd(a+1,b) )
             End Do
          End Do
       End If
       qx = qx + dq/vcv
       qy = qy + dqy/vcv
    End Do
    Call fill_pad_box(Cp, n1, n2, n3)
    Deallocate( uqx, uqy )

  Contains

    Function Qc_x(a, b) Result(q)

      Integer(Int32), Intent(In) :: a, b
      Real(Int64) :: q, mc, ub

      mc = 0.5d0*( Mx(a-1,b) + Mx(a,b) )
      ub = pick(xf(uqx,a-2,b), xf(uqx,a-1,b), xf(uqx,a,b), xf(uqx,a+1,b), mc)
      q = mc*ub

    End Function Qc_x

    Function Qc_y(a, b) Result(q)

      Integer(Int32), Intent(In) :: a, b
      Real(Int64) :: q, mc, ub

      mc = 0.5d0*( My(a,b-1) + My(a,b) )
      ub = pick(yf(uqy,a,b-2), yf(uqy,a,b-1), yf(uqy,a,b), yf(uqy,a,b+1), mc)
      q = mc*ub

    End Function Qc_y

    !> x-sweep momentum flux of qy through the vertical cv face at the corner left (side=0) / right (side=1) of column a, row b+1/2
    Function Qcorner_x(a, b, side) Result(q)

      Integer(Int32), Intent(In) :: a, b, side
      Real(Int64) :: q, mm_, ub
      Integer(Int32) :: ia

      ia = a - 1 + side
      mm_ = 0.5d0*( Mx(ia,b) + Mx(ia,b+1) )
      If ( ia >= 1 .And. ia <= n-1 ) Then
         ub = pick(yf(uqy,ia-1,b), yf(uqy,ia,b), yf(uqy,ia+1,b), yf(uqy,ia+2,b), mm_)
      Else
         ub = 0d0
      End If
      q = mm_*ub

    End Function Qcorner_x

    !> y-sweep momentum flux of qx through the horizontal cv face below (side=0) / above (side=1) row b, at x-face a+1/2
    Function Qcorner_y(a, b, side) Result(q)

      Integer(Int32), Intent(In) :: a, b, side
      Real(Int64) :: q, mm_, ub
      Integer(Int32) :: ib

      ib = b - 1 + side
      mm_ = 0.5d0*( My(a,ib) + My(a+1,ib) )
      If ( ib >= 1 .And. ib <= n-1 ) Then
         ub = pick(xf(uqx,a,ib-1), xf(uqx,a,ib), xf(uqx,a,ib+1), xf(uqx,a,ib+2), mm_)
      Else
         ub = 0d0
      End If
      q = mm_*ub

    End Function Qcorner_y

  End Subroutine advect_half


  !> Project the not-divergence-free velocity (qx,qy)/rho_f with the current density (exact variable-coefficient solve)
  Subroutine project_now(uxo, uyo, ux_out, uy_out)

    Real(Int64), Intent(In)  :: uxo(n-1,n), uyo(n,n-1)
    Real(Int64), Intent(Out) :: ux_out(n-1,n), uy_out(n,n-1)
    Real(Int64) :: rho(n,n), f(n,n), phi(n,n), gx(n-1,n), gy(n,n-1)

    Call rho_cells(rho)
    Call divergence(uxo, uyo, f)
    Call exact_solve(rho, f, phi)
    Call grad_beta(rho, phi, gx, gy)
    ux_out = uxo - gx
    uy_out = uyo - gy

  End Subroutine project_now


  !> RK3 (Wray) force sub-step with frozen density: gravity (hydrostatic part removed) + pressure, V2 + mcg PCG iterations, K=0
  Subroutine forces_rk3(dt, ux_io, uy_io)

    Real(Int64), Intent(In)    :: dt
    Real(Int64), Intent(InOut) :: ux_io(n-1,n), uy_io(n,n-1)
    Real(Int64) :: rho(n,n), rfx(n-1,n), rfy(n,n-1), gv(n,n-1), gxs(n-1,n,3), gys(n,n-1,3)
    Real(Int64) :: usx(n-1,n), usy(n,n-1), bx(n-1,n), by(n,n-1), f(n,n), dl(n,n), dx_(n-1,n), dy_(n,n-1)
    Real(Int64) :: uox(n-1,n), uoy(n,n-1), beta0, ass, rsf, p
    Integer(Int32) :: s, j, a, b

    Call rho_cells(rho)
    Call face_rho(rho, rfx, rfy)
    beta0 = 1d0/rho_g
    Do b = 1, n-1
       rsf = 0.5d0*( rhos(b) + rhos(b+1) )
       Do a = 1, n
          gv(a,b) = -gacc*( rfy(a,b) - rsf )/rfy(a,b)
       End Do
    End Do
    uox = ux_io;  uoy = uy_io
    Do s = 1, 3
       ass = acoef(s,s)
       usx = uox
       usy = uoy + dt*ass*gv
       Do j = 1, s-1
          usx = usx - dt*acoef(s,j)*gxs(:,:,j)
          usy = usy + dt*acoef(s,j)*gv - dt*acoef(s,j)*gys(:,:,j)
       End Do
       Call grad_beta(rho, ppre, bx, by)
       Call divergence(usx - dt*ass*bx, usy - dt*ass*by, f)
       If ( mcg == 0 ) Then
          Call fast_solve(beta0, f/(dt*ass), dl)
          Call grad_beta(rho, dl, dx_, dy_, beta0)
       Else
          Call cg_iterations(rho, f/(dt*ass), beta0, mcg, dl)
          Call grad_beta(rho, dl, dx_, dy_)
       End If
       gxs(:,:,s) = bx + dx_
       gys(:,:,s) = by + dy_
       ppre = ppre + dl
       p = 0d0
       ux_io = usx - dt*ass*gxs(:,:,s)
       uy_io = usy - dt*ass*gys(:,:,s)
    End Do

  End Subroutine forces_rk3


  !> Exact dynamic pressure of the initial state (u = 0): D(beta G p) = D(g_v), used as the first pressure predictor
  Subroutine init_pressure()

    Real(Int64) :: rho(n,n), rfx(n-1,n), rfy(n,n-1), gv(n,n-1), f(n,n), zx(n-1,n)
    Real(Int64) :: rsf
    Integer(Int32) :: a, b

    Call rho_cells(rho)
    Call face_rho(rho, rfx, rfy)
    Do b = 1, n-1
       rsf = 0.5d0*( rhos(b) + rhos(b+1) )
       Do a = 1, n
          gv(a,b) = -gacc*( rfy(a,b) - rsf )/rfy(a,b)
       End Do
    End Do
    zx = 0d0
    Call divergence(zx, gv, f)
    Call exact_solve(rho, f, ppre)

  End Subroutine init_pressure


  Function rho_row(b) Result(r)

    Integer(Int32), Intent(In) :: b
    Real(Int64) :: r(n)
    Integer(Int32) :: a

    Do a = 1, n
       r(a) = rho_of(Cp(a+1,b+1,3))
    End Do

  End Function rho_row


  Subroutine hydrostatic_profile()

    Integer(Int32) :: b

    Do b = 1, n
       rhos(b) = Sum(rho_row(b))/n
    End Do

  End Subroutine hydrostatic_profile


  Subroutine momentum_from_velocity()

    Real(Int64) :: rho(n,n), rfx(n-1,n), rfy(n,n-1)

    Call rho_cells(rho)
    Call face_rho(rho, rfx, rfy)
    qx = rfx*ux
    qy = rfy*uy

  End Subroutine momentum_from_velocity


  !> control (form 0): C is transported with the same Strang sweeps but momentum is not (velocity kept, re-projected with the new
  !  density), i.e. no convection of momentum -- equivalent at linear order to a velocity-form scheme for a small-amplitude wave
  Subroutine step_form0(dt, istep)

    Real(Int64), Intent(In) :: dt
    Integer(Int32), Intent(In) :: istep
    Real(Int64) :: uxn(n-1,n), uyn(n,n-1), qs_x(n-1,n), qs_y(n,n-1)

    Call momentum_from_velocity()
    qs_x = qx;  qs_y = qy
    Call advect_half(0.5d0*dt, istep, ux, uy)
    qx = qs_x;  qy = qs_y
    uxn = ux;  uyn = uy
    Call project_now(uxn, uyn, ux, uy)
    Call forces_rk3(dt, ux, uy)
    Call momentum_from_velocity()
    Call advect_half(0.5d0*dt, istep+1, ux, uy)
    uxn = ux;  uyn = uy
    Call project_now(uxn, uyn, ux, uy)

  End Subroutine step_form0


  !> one full Strang step with conservative (form 1) momentum transport
  Subroutine step_form1(dt, istep)

    Real(Int64), Intent(In) :: dt
    Integer(Int32), Intent(In) :: istep
    Real(Int64) :: rho(n,n), rfx(n-1,n), rfy(n,n-1), uxn(n-1,n), uyn(n,n-1)

    Call momentum_from_velocity()
    Call advect_half(0.5d0*dt, istep, ux, uy)
    Call rho_cells(rho);  Call face_rho(rho, rfx, rfy)
    uxn = qx/rfx;  uyn = qy/rfy
    Call project_now(uxn, uyn, ux, uy)
    Call forces_rk3(dt, ux, uy)
    Call momentum_from_velocity()
    Call advect_half(0.5d0*dt, istep+1, ux, uy)
    Call rho_cells(rho);  Call face_rho(rho, rfx, rfy)
    uxn = qx/rfx;  uyn = qy/rfy
    Call project_now(uxn, uyn, ux, uy)

  End Subroutine step_form1


  Subroutine energies(ke, pe)

    Real(Int64), Intent(Out) :: ke, pe
    Real(Int64) :: rho(n,n), rfx(n-1,n), rfy(n,n-1)
    Integer(Int32) :: a, b

    Call rho_cells(rho);  Call face_rho(rho, rfx, rfy)
    ke = 0.5d0*h**3*( Sum(rfx*ux**2) + Sum(rfy*uy**2) )
    pe = 0d0
    Do b = 1, n
       Do a = 1, n
          pe = pe + rho(a,b)*gacc*(b-0.5d0)*h*h**3
       End Do
    End Do

  End Subroutine energies


  Subroutine init_wave(amp)

    Real(Int64), Intent(In) :: amp
    Integer(Int32) :: a, b, k
    Real(Int64) :: xc, eta, s, x0, y0, beta

    Do b = 1, n
       Do a = 1, n
          xc = (a-0.5d0)*h
          eta = amp*Cos(pi*xc)
          s = -amp*pi*Sin(pi*xc)
          x0 = (a-1)*h;  y0 = (b-1)*h
          beta = 0.5d0 + eta - s*xc - y0 + s*x0
          Do k = 1, n3
             Cp(a+1,b+1,k) = plic_box_fraction(-s*h, h, 0d0, beta)
          End Do
       End Do
    End Do
    Call fill_pad_box(Cp, n1, n2, n3)

  End Subroutine init_wave


  !> sloshing wave (amplitude 2.4 cells, N = 48): half-period error versus the analytic two-layer linear theory
  Subroutine slosh(rr, fm, ms, nn, per_err, maxu, nc)

    Real(Int64),    Intent(In)  :: rr
    Integer(Int32), Intent(In)  :: fm, ms, nn
    Real(Int64),    Intent(Out) :: per_err, maxu
    Integer(Int32), Intent(Out) :: nc
    Real(Int64) :: dt, per, wth, t, vl, vlold, tc(20), dt0, umx
    Integer(Int32) :: istep

    mom_scheme = ms
    Call setup(nn, 1d0, 1d0/rr, 13.7d0)
    Call init_wave(0.05d0)
    Call hydrostatic_profile()
    Call init_pressure()
    ux = 0d0;  uy = 0d0
    wth = Sqrt( gacc*pi*(rho_l - rho_g)/( (rho_l + rho_g)/Tanh(pi*0.5d0) ) )
    per = 2d0*pi/wth
    dt = per/64d0
    nc = 0;  vlold = 0d0;  maxu = 0d0
    per_err = 1d30
    t = 0d0;  istep = 0;  dt0 = dt
    Do While ( t < 0.8d0*per .And. istep < 200000 )
       If ( cfl_lim > 0d0 ) Then
          umx = Max(Maxval(Abs(ux)), Maxval(Abs(uy)), 1d-3)
          dt = Min(dt0, cfl_lim*h/umx)
       End If
       istep = istep + 1
       If ( fm == 0 ) Then
          Call step_form0(dt, 2*istep)
       Else
          Call step_form1(dt, 2*istep)
       End If
       t = t + dt
       vl = h*Sum(Cp(2,2:n+1,3)) - 0.5d0
       If ( vl /= vl ) Exit
       maxu = Max(maxu, Maxval(Abs(ux)), Maxval(Abs(uy)))
       If ( istep > 1 .And. vl*vlold < 0d0 .And. nc < 20 ) Then
          nc = nc + 1;  tc(nc) = t - dt*vl/(vl - vlold)
       End If
       vlold = vl
    End Do
    nsteps_used = istep
    If ( nc >= 2 ) per_err = (tc(2) - tc(1))/(0.5d0*per) - 1d0

  End Subroutine slosh


  Subroutine run_all

    Real(Int64) :: dt, e, mu, rr(5)
    Integer(Int32) :: istep, b, nc, ms, ir, fm, ig, nlist(3) = (/ 24, 48, 96 /)
    Real(Int64), Parameter :: snaps(4) = (/ 0d0, 0.05d0, 0.1d0, 0.2d0 /)
    Character(len=8) :: arg

    ! --- T1: still water, interface inside a cell
    Call setup(48, 1d0, 1d-3, 10d0)
    Do b = 1, n
       Cp(:,b+1,:) = Min(1d0, Max(0d0, (0.5d0 + 0.3d0*h - (b-1)*h)/h))
    End Do
    Call fill_pad_box(Cp, n1, n2, n3)
    Call hydrostatic_profile()
    Call init_pressure()
    ux = 0d0;  uy = 0d0
    dt = 0.01d0
    Do istep = 0, 99
       Call step_form1(dt, 2*istep)
    End Do
    Write(*,'(a,es10.2)') 'T1 still water (ratio 1000), max|u| after 100 steps:', Max(Maxval(Abs(ux)), Maxval(Abs(uy)))

    ! --- T3: sloshing wave, mode k = pi in a unit box, interface at y = 0.5, 2 periods
    rr = (/ 10d0, 100d0, 250d0, 500d0, 1000d0 /)
    If ( Command_argument_count() == 0 ) Write(*, &
         '(a)') 'T3 sloshing: half-period error (measured/analytic - 1) and max|u|; 99 = unstable / no oscillation'
    Do fm = 0, Merge(-1, 1, Command_argument_count() > 0)
       Do ms = 0, Merge(2, 0, fm == 1)
          Write(*,'(a,i2,a,i2)') ' momentum form', fm, '   interpolation scheme', ms
          Do ir = 1, 5
             Call slosh(rr(ir), fm, ms, 48, e, mu, nc)
             If ( nc < 2 ) Then
                Write(*,'(a,f7.0,a)') '   ratio', rr(ir), '   99'
             Else
                Write(*,'(a,f7.0,a,f9.4,a,es9.2)') '   ratio', rr(ir), '   period err', e, '   max|u|', mu
             End If
          End Do
       End Do
    End Do

    arg = ' '
    If ( Command_argument_count() > 0 ) Call Get_command_argument(1, arg)
    If ( arg == 'snap' ) Then
       Write(*,'(a)') 'T3 density regularisation (form 0 control): period error versus snap threshold and N'
       Do ir = 4, 5
          Do ig = 1, 4
             Do ms = 1, 2
                snap = snaps(ig)
                Call slosh(rr(ir), 0, 0, 24*2**ms, e, mu, nc)
                Write(*,'(a,f7.0,a,f5.2,a,i4,a,f9.4,a,es9.2)') '   ratio', rr(ir), '  snap', snap, '  N', 24*2**ms, &
                     '   period err', e, '   max|u|', mu
             End Do
          End Do
       End Do
       snap = 0d0
    End If
    If ( arg == 'cfl' ) Then
       Write(*,'(a)') 'T3 with adaptive dt (Co <= 0.4), N = 48: half-period error and max|u| versus ratio and momentum treatment'
       cfl_lim = 0.4d0
       Do ir = 2, 5, 3
          Call slosh(rr(ir), 0, 0, 48, e, mu, nc)
          Write(*,'(a,f7.0,a,f9.4,a,es9.2,a,i7)') '   ratio', rr(ir), '  form 0 (no momentum transport) period err', e, &
               '  max|u|', mu, '  steps', nsteps_used
          Do ms = 0, 2
             Call slosh(rr(ir), 1, ms, 48, e, mu, nc)
             If ( nc < 2 ) Then
                Write(*,'(a,f7.0,a,i2,a,es9.2,a,i7)') '   ratio', rr(ir), '  form 1 scheme', ms, &
                     '  no 2nd zero-crossing (stalled/damped)  max|u|', mu, '  steps', nsteps_used
             Else
                Write(*,'(a,f7.0,a,i2,a,f9.4,a,es9.2,a,i7)') '   ratio', rr(ir), '  form 1 scheme', ms, '  period err', &
                     e, '  max|u|', mu, '  steps', nsteps_used
             End If
          End Do
       End Do
       cfl_lim = 0d0
    End If
    If ( arg == 'conv' ) Then
       Write(*,'(a)') 'T3 convergence (form 1, upwind), amplitude 0.05: period error versus N'
       Do ir = 2, 5, 3
          Do ig = 1, 3
             Call slosh(rr(ir), 1, 1, nlist(ig), e, mu, nc)
             Write(*,'(a,f7.0,a,i4,a,f9.4,a,es9.2)') '   ratio', rr(ir), '  N', nlist(ig), '   period err', e, '   max|u|', mu
          End Do
       End Do
    End If

  End Subroutine run_all


End Module spike_c_mod


Program spike_c_momentum
  Use spike_c_mod, Only : run_all
  Implicit None
  Call run_all
End Program spike_c_momentum
