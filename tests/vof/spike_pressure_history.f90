!> Spike A'' (closed loop): the solver's RK3 (Wray) with a projection at every stage, a frozen two-fluid density field and an
!  oscillating non-gradient body force, comparing the production-compatible "explicit extrapolated pressure force + constant-
!  coefficient projection" (V2) against the exact variable-density projection at every stage.
!
!  Stage s of the low-storage scheme:  u^s = u^o + dt * sum_{j<=s} a(s,j) (F_j - G_j),  G_j = applied pressure force at stage j.
!  Reference: G_s = beta grad p_s with D(beta G p_s) solved exactly. V2: p_s = ptilde_s + delta_s, G_s = beta grad ptilde_s +
!  beta0 grad delta_s, where ptilde_s is a Lagrange extrapolation (order K) of the stored stage pressures to the stage time and delta_s
!  comes from the fast constant-coefficient solve, so u^s is divergence-free to round-off. Two stage-time conventions for the
!  stored pressures are compared (A: the time of u^s, rk_t; B: the time at which F_s is evaluated).
!
!  Results (N=48, ratio 10 and 1000, w*dt = 0.05..0.4, liquid-forced frozen two-fluid field, errors vs exact projection at every stage):
!   - No pressure predictor (the existing projection used as is): 42-49 % velocity error in both phases at every dt.
!   - beta0 must be the largest 1/rho (rho0 = lighter fluid); beta0 = 1/rho_l or the mean makes the pressure history diverge,
!     because the stored pressure receives only beta/beta0 of the true correction (contraction 1 - beta/beta0 ~ 0.999 in the liquid).
!   - Stage-pressure time stamps: convention B (the time at which F_s is evaluated) is exact: p_s = p*(tau_B) to round-off for this
!     linear problem, so a-priori extrapolation errors are 2e-4 (K=1) and 1e-5 (K=2) at w*dt = 0.05. Convention A is wrong.
!   - Closed loop, single fast solve (V2): K=0 gives 2.6 % (w*dt = 0.1) and K=1 1.1 %; first order in dt, independent of the density
!     ratio, a steady bias (identical every period over 10 periods), not secular drift. K=2 is unstable.
!   - V2 plus m PCG iterations on the increment (fast-solve preconditioner), K=0: m=1 costs the same number of fast solves as V2 and
!     is ~20x more accurate (1.2e-3) but no longer exactly divergence-free (~2e-3 of the force scale); m=2: 1.3e-4 liquid;
!     m=3: 3.6e-5 liquid / 6.7e-4 gas. K>=1 with m>=1 is unstable at ratio 1000: use K=0.
Program spike_pressure_history

  Use iso_fortran_env, Only : Int32, Int64
  Use spike_proj_lib
  Implicit None

  Real(Int64), Parameter :: acoef(3,3) = Reshape( (/ 8d0/15d0, 1d0/4d0, 1d0/4d0, &
                                                     0d0,      5d0/12d0, 0d0,    &
                                                     0d0,      0d0,      3d0/4d0 /), (/3,3/) )
  Real(Int64), Parameter :: rkt(3) = (/ 8d0/15d0, 2d0/3d0, 1d0 /), tfr(3) = (/ 0d0, 8d0/15d0, 2d0/3d0 /)
  Real(Int64), Parameter :: omega = 2d0*pi, amp = 1d0
  Integer(Int32) :: nper = 2
  Real(Int64) :: ratios(2) = (/ 10d0, 1000d0 /), wdts(4) = (/ 0.4d0, 0.2d0, 0.1d0, 0.05d0 /)
  Real(Int64), Allocatable :: rho(:,:), fsx(:,:), fsy(:,:)
  Real(Int64), Allocatable :: uxr(:,:,:), uyr(:,:,:), pr(:,:,:), uxv(:,:,:), uyv(:,:,:), pv(:,:,:)
  Real(Int64) :: rho_l, rho_g, dt, beta0, b0s(3), eps_ap, el, eg, divmax
  Integer(Int32) :: ir, iw, nsteps, K, mode, iters

  Call proj_init(48)
  Allocate( rho(n,n), fsx(n-1,n), fsy(n,n-1) )

  If ( Command_argument_count() > 0 ) Then
     Call long_run
     Stop
  End If

  Write(*,'(a)') 'ratio  w*dt  tau  K nCG   eps_apriori   err_liq    err_gas    max|div|*h/|f|   (nCG=0: V2 fast solve only)'
  Write(*,'(a)') '(velocity errors are relative to the reference u of the same phase)'
  Do ir = 1, 2
     rho_l = 1d0;  rho_g = rho_l/ratios(ir)
     Call build_density(rho_l, rho_g, rho)
     Call make_force(rho, fsx, fsy)
     b0s = (/ 1d0/rho_l, 1d0/rho_g, Sum(1d0/rho)/(n*n) /)
     Do iw = 1, 4
        nsteps = Nint(nper*2d0*pi/wdts(iw))
        dt = Real(nper,Int64)/nsteps
        Allocate( uxr(n-1,n,nsteps), uyr(n,n-1,nsteps), pr(n,n,3*nsteps), uxv(n-1,n,nsteps), uyv(n,n-1,nsteps), pv(n,n,3*nsteps) )
        Call simulate(0, nsteps, dt, rho, 0d0, 0, 1, uxr, uyr, pr, divmax, iters, 0)
        Do K = -1, 2
           Call apriori_eps(pr, nsteps, dt, K, .False., eps_ap)
           Do mode = 1, 4
              ! mode 1: V2 (fast solve, exactly div-free); modes 2-4: 1-3 PCG iterations on the increment (hybrid)
              beta0 = b0s(2)
              If ( mode == 1 ) Then
                 Call simulate(1, nsteps, dt, rho, beta0, K, 0, uxv, uyv, pv, divmax, iters, 0)
              Else
                 Call simulate(2, nsteps, dt, rho, beta0, K, 0, uxv, uyv, pv, divmax, iters, mode-1)
              End If
              Call vel_error(rho, rho_l, uxr, uyr, uxv, uyv, nsteps/2 + 1, nsteps, el, eg)
              Write(*,'(f6.0,f6.2,2x,a,i3,i4,es12.3,3es11.3)') ratios(ir), wdts(iw), 'B', K, mode-1, eps_ap, el, eg, divmax
           End Do
        End Do
        Deallocate( uxr, uyr, pr, uxv, uyv, pv )
     End Do
  End Do

Contains

  !> Error per period over 10 periods (ratio 1000, w*dt = 0.1): does the single-solve error accumulate?
  Subroutine long_run

    Integer(Int32) :: ip, kk, md, ns, spp
    Real(Int64) :: rl, rg, dtl, b0, e1, e2, dv
    Real(Int64), Allocatable :: a1(:,:,:), a2(:,:,:), p1(:,:,:), b1(:,:,:), b2(:,:,:), p2(:,:,:)
    Integer(Int32) :: it

    nper = 10
    rl = 1d0;  rg = 1d-3
    Call build_density(rl, rg, rho)
    Call make_force(rho, fsx, fsy)
    ns = Nint(nper*2d0*pi/0.1d0);  dtl = Real(nper,Int64)/ns
    spp = ns/nper
    Allocate( a1(n-1,n,ns), a2(n,n-1,ns), p1(n,n,3*ns), b1(n-1,n,ns), b2(n,n-1,ns), p2(n,n,3*ns) )
    b0 = 1d0/rg
    Call simulate(0, ns, dtl, rho, 0d0, 0, 0, a1, a2, p1, dv, it, 0)
    Write(*,'(a)') 'ratio 1000, w*dt=0.1, 10 periods: liquid / gas velocity error per period'
    Do kk = 0, 1
       Do md = 0, 3
          If ( md == 0 ) Then
             Call simulate(1, ns, dtl, rho, b0, kk, 0, b1, b2, p2, dv, it, 0)
          Else
             Call simulate(2, ns, dtl, rho, b0, kk, 0, b1, b2, p2, dv, it, md)
          End If
          Write(*,'(a,i2,a,i2)') 'K=', kk, '  nCG=', md
          Do ip = 1, nper
             Call vel_error(rho, rl, a1, a2, b1, b2, (ip-1)*spp + 1, ip*spp, e1, e2)
             Write(*,'(i4,2es12.3)') ip, e1, e2
          End Do
       End Do
    End Do

  End Subroutine long_run


  Subroutine build_density(rl, rg, rho)

    Real(Int64), Intent(In)  :: rl, rg
    Real(Int64), Intent(Out) :: rho(n,n)
    Integer(Int32) :: i, j
    Real(Int64) :: x, y

    Do j = 1, n
       y = (j-0.5d0)*h
       Do i = 1, n
          x = (i-0.5d0)*h
          rho(i,j) = Merge(rl, rg, y < 0.5d0 + 0.05d0*Cos(2d0*pi*x))
       End Do
    End Do

  End Subroutine build_density


  !> Acceleration beta_f * S with the solenoidal per-volume force S = liquid fraction * (psi_y, -psi_x), psi = sin(pi x) sin(pi y): a force per unit
  !  volume acting on a density jump is not a gradient per unit mass, so a pressure (baroclinic-like) is genuinely required
  Subroutine make_force(rho, fx, fy)

    Real(Int64), Intent(In)  :: rho(n,n)
    Real(Int64), Intent(Out) :: fx(n-1,n), fy(n,n-1)
    Integer(Int32) :: i, j

    Do j = 1, n
       Do i = 1, n-1
          fx(i,j) = 2d0/(rho(i,j) + rho(i+1,j))*cliq(rho(i,j), rho(i+1,j))*pi*Sin(pi*i*h)*Cos(pi*(j-0.5d0)*h)
       End Do
    End Do
    Do j = 1, n-1
       Do i = 1, n
          fy(i,j) = -2d0/(rho(i,j) + rho(i,j+1))*cliq(rho(i,j), rho(i,j+1))*pi*Cos(pi*(i-0.5d0)*h)*Sin(pi*j*h)
       End Do
    End Do

  End Subroutine make_force


  !> liquid fraction of a face (mean of its two cells; the liquid has the larger density)
  Pure Function cliq(ra, rb) Result(c)

    Real(Int64), Intent(In) :: ra, rb
    Real(Int64) :: c

    c = 0.5d0*( Merge(1d0, 0d0, ra >= 1d0) + Merge(1d0, 0d0, rb >= 1d0) )

  End Function cliq


  Pure Function stage_time(step, s, optA, dt) Result(t)

    Integer(Int32), Intent(In) :: step, s
    Logical, Intent(In) :: optA
    Real(Int64), Intent(In) :: dt
    Real(Int64) :: t

    If ( optA ) Then
       t = (step-1)*dt + rkt(s)*dt
    Else
       t = (step-1)*dt + tfr(s)*dt
    End If

  End Function stage_time


  !> Lagrange extrapolation of order K to time tt from the (up to K+1) stored pressures preceding entry idx (idx counts steps*3+stage)
  Subroutine extrapolate(hist, idx, dt, K, optA, tt, pt)

    Real(Int64), Intent(In)  :: hist(n,n,*), dt, tt
    Integer(Int32), Intent(In) :: idx, K
    Logical, Intent(In) :: optA
    Real(Int64), Intent(Out) :: pt(n,n)
    Integer(Int32) :: m, i, j, ii, jj
    Real(Int64) :: w, ti, tj

    pt = 0d0
    If ( K < 0 ) Return
    m = Min(K+1, idx-1)
    Do i = 1, m
       ii = idx - i
       ti = stage_time((ii-1)/3 + 1, Mod(ii-1,3) + 1, optA, dt)
       w = 1d0
       Do j = 1, m
          If ( j == i ) Cycle
          jj = idx - j
          tj = stage_time((jj-1)/3 + 1, Mod(jj-1,3) + 1, optA, dt)
          w = w*(tt - tj)/(ti - tj)
       End Do
       pt = pt + w*hist(:,:,ii)
    End Do

  End Subroutine extrapolate


  Subroutine simulate(mode, nsteps, dt, rho, beta0, K, optAi, uxh, uyh, phist, divmax, iters, nit_cg)

    Integer(Int32), Intent(In) :: mode, nsteps, K, optAi, nit_cg
    Real(Int64), Intent(In) :: dt, beta0, rho(n,n)
    Real(Int64), Intent(Out) :: uxh(n-1,n,nsteps), uyh(n,n-1,nsteps), phist(n,n,*), divmax
    Integer(Int32), Intent(Out) :: iters
    Real(Int64) :: ux(n-1,n), uy(n,n-1), uox(n-1,n), uoy(n,n-1), usx(n-1,n), usy(n,n-1), gx(n-1,n,3), gy(n,n-1,3), cf(3)
    Real(Int64) :: pt(n,n), ps(n,n), f(n,n), dl(n,n), bx(n-1,n), by(n,n-1), dx_(n-1,n), dy_(n,n-1), d(n,n), t0, tt, ass
    Integer(Int32) :: step, s, j, idx, it
    Logical :: optA

    optA = optAi == 1
    ux = 0d0;  uy = 0d0
    divmax = 0d0;  iters = 0
    Do step = 1, nsteps
       t0 = (step-1)*dt
       uox = ux;  uoy = uy
       Do s = 1, 3
          cf(s) = amp*Cos(omega*(t0 + tfr(s)*dt))
          ass = acoef(s,s)
          usx = uox + dt*ass*cf(s)*fsx
          usy = uoy + dt*ass*cf(s)*fsy
          Do j = 1, s-1
             usx = usx + dt*acoef(s,j)*(cf(j)*fsx - gx(:,:,j))
             usy = usy + dt*acoef(s,j)*(cf(j)*fsy - gy(:,:,j))
          End Do
          idx = (step-1)*3 + s
          tt = stage_time(step, s, optA, dt)
          If ( mode == 0 ) Then
             Call divergence(usx, usy, f)
             Call exact_solve(rho, f/(dt*ass), ps, it)
             iters = iters + it
             Call grad_beta(rho, ps, gx(:,:,s), gy(:,:,s))
          Else
             Call extrapolate(phist, idx, dt, K, optA, tt, pt)
             Call grad_beta(rho, pt, bx, by)
             Call divergence(usx - dt*ass*bx, usy - dt*ass*by, f)
             If ( mode == 2 ) Then
                Call cg_iterations(rho, f/(dt*ass), beta0, nit_cg, dl)
                Call grad_beta(rho, dl, dx_, dy_)
             Else
                Call fast_solve(beta0, f/(dt*ass), dl)
                Call grad_beta(rho, dl, dx_, dy_, beta0)
             End If
             gx(:,:,s) = bx + dx_
             gy(:,:,s) = by + dy_
             ps = pt + dl
          End If
          phist(:,:,idx) = ps
          If ( s == 3 ) Then
             ux = usx - dt*ass*gx(:,:,s)
             uy = usy - dt*ass*gy(:,:,s)
             Call divergence(ux, uy, d)
             divmax = Max(divmax, Maxval(Abs(d))*h/Max(Maxval(Abs(fsx)), 1d-300))
          End If
       End Do
       uxh(:,:,step) = ux
       uyh(:,:,step) = uy
    End Do

  End Subroutine simulate


  !> max over the last period of |p - ptilde| / max|p| using the reference pressure history itself
  Subroutine apriori_eps(pr, nsteps, dt, K, optA, eps)

    Real(Int64), Intent(In) :: pr(n,n,*), dt
    Integer(Int32), Intent(In) :: nsteps, K
    Logical, Intent(In) :: optA
    Real(Int64), Intent(Out) :: eps
    Real(Int64) :: pt(n,n), pmax
    Integer(Int32) :: step, s, idx

    eps = 0d0;  pmax = 0d0
    Do step = nsteps/2 + 1, nsteps
       Do s = 1, 3
          idx = (step-1)*3 + s
          Call extrapolate(pr, idx, dt, K, optA, stage_time(step, s, optA, dt), pt)
          eps = Max(eps, Maxval(Abs(pr(:,:,idx) - pt)))
          pmax = Max(pmax, Maxval(Abs(pr(:,:,idx))))
       End Do
    End Do
    eps = eps/pmax

  End Subroutine apriori_eps


  Subroutine vel_error(rho, rho_l, uxr, uyr, uxv, uyv, i1, i2, el, eg)

    Real(Int64), Intent(In) :: rho(n,n), rho_l, uxr(n-1,n,*), uyr(n,n-1,*), uxv(n-1,n,*), uyv(n,n-1,*)
    Integer(Int32), Intent(In) :: i1, i2
    Real(Int64), Intent(Out) :: el, eg
    Real(Int64) :: nl, ng
    Integer(Int32) :: step, i, j
    Logical :: liq

    el = 0d0;  eg = 0d0;  nl = 0d0;  ng = 0d0
    Do step = i1, i2
       Do j = 1, n
          Do i = 1, n-1
             If ( rho(i,j) > 0.5d0*rho_l .And. rho(i+1,j) > 0.5d0*rho_l ) Then
                el = Max(el, Abs(uxv(i,j,step) - uxr(i,j,step)));  nl = Max(nl, Abs(uxr(i,j,step)))
             Else If ( rho(i,j) < 0.5d0*rho_l .And. rho(i+1,j) < 0.5d0*rho_l ) Then
                eg = Max(eg, Abs(uxv(i,j,step) - uxr(i,j,step)));  ng = Max(ng, Abs(uxr(i,j,step)))
             End If
          End Do
       End Do
       Do j = 1, n-1
          Do i = 1, n
             liq = rho(i,j) > 0.5d0*rho_l .And. rho(i,j+1) > 0.5d0*rho_l
             If ( liq ) Then
                el = Max(el, Abs(uyv(i,j,step) - uyr(i,j,step)));  nl = Max(nl, Abs(uyr(i,j,step)))
             Else If ( rho(i,j) < 0.5d0*rho_l .And. rho(i,j+1) < 0.5d0*rho_l ) Then
                eg = Max(eg, Abs(uyv(i,j,step) - uyr(i,j,step)));  ng = Max(ng, Abs(uyr(i,j,step)))
             End If
          End Do
       End Do
    End Do
    el = el/Max(nl, 1d-300)
    eg = eg/Max(ng, 1d-300)

  End Subroutine vel_error

End Program spike_pressure_history
