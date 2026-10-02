!> Wave generation for the two-fluid solver: target free-surface elevation and velocity of a long-crested wave travelling in +x over
!  the still-water depth d = vof_level (the bed is y = 0), used for the Dirichlet inlet and the relaxation zones.
!
!    wave_type 1  linear (Airy) wave of height wave_height and period wave_period
!    wave_type 2  Rienecker-Fenton stream-function wave of the same height and period (nonlinear)
!    wave_type 3  JONSWAP spectrum (wave_Hs, wave_Tp, wave_gamma) as a sum of wave_nfreq linear components with seeded random
!                 phases, finite-depth dispersion and Wheeler stretching of the velocity profile
!
!  Above the surface the water velocity is not defined; the air carries the uniform return flow that makes the inlet flux zero (a
!  closed flume): u_air = -Q/(Ly - d - eta), Q the volume flux of the water column. All setup (dispersion roots, spectrum, phases,
!  the stream-function solve) runs on every rank from the same inputs and the same seed, so no communication is needed.
Module waves

  Use iso_fortran_env, Only : Int32, Int64
  Use global, Only : nyg, ny, y, yg, Ly_i, vof_level, vof_grav, wave_type, wave_height, wave_period, wave_phase, wave_sf_n, &
                     wave_current_mode, wave_Hs, wave_Tp, wave_gamma, wave_nfreq, wave_seed, wave_ramp_time

  Implicit None

  Integer(Int32), Parameter :: WAVE_MAX_COMP = 4096, SF_MAX_N = 64
  Integer(Int32) :: wv_n = 0                   ! number of linear components
  Real(Int64), Allocatable, Dimension(:) :: wv_a, wv_k, wv_w, wv_ph
  Real(Int64) :: wv_d = 0d0, wv_g = 9.81d0
  ! stream-function solution (wave_type 2): wave number, phase speed, Eulerian mean-current definition, coefficients
  Integer(Int32) :: sf_n = 0
  Real(Int64) :: sf_k = 0d0, sf_c = 0d0, sf_ubar = 0d0, sf_q = 0d0, sf_r = 0d0
  Real(Int64) :: sf_b(0:SF_MAX_N) = 0d0, sf_e(0:SF_MAX_N) = 0d0
  ! inlet profiles on the solver's y grid, refreshed by wave_inlet_profile
  Real(Int64), Allocatable, Dimension(:) :: wv_in_u, wv_in_v
  Real(Int64) :: wv_in_eta = 0d0

Contains

  !> Wave number from the finite-depth dispersion relation w^2 = g k tanh(k d)
  Function wave_number(om, d, g) Result(k)

    Real(Int64), Intent(In) :: om, d, g
    Real(Int64) :: k, f, df
    Integer(Int32) :: it

    k = om*om/g
    If ( k*d < 1d0 ) k = om/Sqrt(g*d)
    Do it = 1, 100
       f  = g*k*Tanh(k*d) - om*om
       df = g*Tanh(k*d) + g*k*d/Cosh(k*d)**2
       k = k - f/df
       If ( Abs(f) < 1d-14*om*om ) Exit
    End Do

  End Function wave_number


  Subroutine wave_init

    Integer(Int32) :: i, nf
    Real(Int64) :: om, fp, fmin, fmax, df, f, sigma, r, alpha_j, m0, sj, hs, pi_
    Real(Int64), Allocatable :: sp(:)
    Integer(Int32), Allocatable :: seed(:)
    Integer(Int32) :: nseed

    pi_ = 4d0*Atan(1d0)
    wv_d = vof_level
    wv_g = vof_grav
    Allocate( wv_in_u(nyg), wv_in_v(ny) )
    wv_in_u = 0d0;  wv_in_v = 0d0

    Select Case(wave_type)
    Case(1)
       wv_n = 1
       Allocate( wv_a(1), wv_k(1), wv_w(1), wv_ph(1) )
       om = 2d0*pi_/wave_period
       wv_w(1) = om;  wv_k(1) = wave_number(om, wv_d, wv_g);  wv_a(1) = 0.5d0*wave_height;  wv_ph(1) = wave_phase
    Case(2)
       Call stream_function_solve
    Case(3)
       nf = Min(Max(wave_nfreq, 8), WAVE_MAX_COMP)
       wv_n = nf
       Allocate( wv_a(nf), wv_k(nf), wv_w(nf), wv_ph(nf), sp(nf) )
       fp = 1d0/wave_Tp
       fmin = 0.5d0*fp;  fmax = 3.5d0*fp
       df = (fmax - fmin)/nf
       Call Random_seed(size=nseed)
       Allocate( seed(nseed) )
       seed = wave_seed + 37*[(i, i=1,nseed)]
       Call Random_seed(put=seed)
       m0 = 0d0
       Do i = 1, nf
          f = fmin + (Real(i,Int64) - 0.5d0)*df
          sigma = Merge(0.07d0, 0.09d0, f <= fp)
          r = Exp( -(f - fp)**2/(2d0*sigma*sigma*fp*fp) )
          sj = f**(-5)*Exp(-1.25d0*(fp/f)**4)*wave_gamma**r
          sp(i) = sj
          m0 = m0 + sj*df
       End Do
       ! scale the spectrum so that 4 sqrt(m0) = Hs
       alpha_j = (wave_Hs/4d0)**2/m0
       Do i = 1, nf
          f = fmin + (Real(i,Int64) - 0.5d0)*df
          Call Random_number(r)
          wv_ph(i) = 2d0*pi_*r
          wv_w(i) = 2d0*pi_*f
          wv_k(i) = wave_number(wv_w(i), wv_d, wv_g)
          wv_a(i) = Sqrt(2d0*alpha_j*sp(i)*df)
       End Do
       Deallocate( sp, seed )
    Case Default
       wv_n = 0
    End Select

  End Subroutine wave_init


  !> Start-up amplitude ramp 0 -> 1 over wave_ramp_time (raised cosine)
  Function wave_ramp(t) Result(r)

    Real(Int64), Intent(In) :: t
    Real(Int64) :: r

    r = 1d0
    If ( wave_ramp_time > 0d0 .And. t < wave_ramp_time ) r = 0.5d0*( 1d0 - Cos(4d0*Atan(1d0)*Max(t, 0d0)/wave_ramp_time) )

  End Function wave_ramp


  Function wave_eta(x, t) Result(eta)

    Real(Int64), Intent(In) :: x, t
    Real(Int64) :: eta

    eta = wave_ramp(t)*wave_eta_raw(x, t)

  End Function wave_eta


  Subroutine wave_vel(x, y, t, u, v)

    Real(Int64), Intent(In)  :: x, y, t
    Real(Int64), Intent(Out) :: u, v

    Call wave_vel_raw(x, y, t, u, v)
    u = wave_ramp(t)*u
    v = wave_ramp(t)*v

  End Subroutine wave_vel


  !> Surface elevation above the still level at (x, t), without the start-up ramp
  Function wave_eta_raw(x, t) Result(eta)

    Real(Int64), Intent(In) :: x, t
    Real(Int64) :: eta
    Integer(Int32) :: i, j

    eta = 0d0
    If ( wave_type == 2 ) Then
       Do j = 0, sf_n
          eta = eta + sf_e(j)*Cos(Real(j,Int64)*(sf_k*x - 2d0*4d0*Atan(1d0)*t/wave_period + wave_phase))
       End Do
       Return
    End If
    Do i = 1, wv_n
       eta = eta + wv_a(i)*Cos(wv_k(i)*x - wv_w(i)*t + wv_ph(i))
    End Do

  End Function wave_eta_raw


  !> Water velocity (u, v) at the point (x, y), y from the bed, for y up to the surface (the caller handles the air)
  Subroutine wave_vel_raw(x, y, t, u, v)

    Real(Int64), Intent(In)  :: x, y, t
    Real(Int64), Intent(Out) :: u, v
    Integer(Int32) :: i, j
    Real(Int64) :: th, ye, eta, ch, sh, kd, jj

    u = 0d0;  v = 0d0
    If ( wave_type == 2 ) Then
       th = sf_k*x - 2d0*4d0*Atan(1d0)*t/wave_period + wave_phase
       ye = Min(y, wv_d + wave_eta_raw(x, t))
       u = sf_c + sf_b(0)
       Do j = 1, sf_n
          jj = Real(j,Int64)
          ch = ratio_cosh(jj*sf_k*ye, jj*sf_k*wv_d)
          sh = ratio_sinh(jj*sf_k*ye, jj*sf_k*wv_d)
          u = u + jj*sf_k*sf_b(j)*ch*Cos(jj*th)
          v = v + jj*sf_k*sf_b(j)*sh*Sin(jj*th)
       End Do
       Return
    End If
    eta = wave_eta_raw(x, t)
    Do i = 1, wv_n
       th = wv_k(i)*x - wv_w(i)*t + wv_ph(i)
       ye = y
       If ( wave_type == 3 ) Then
          ! Wheeler stretching: map [0, d+eta] onto [0, d]; above the surface keep the surface value
          ye = Min(y, wv_d + eta)*wv_d/(wv_d + eta)
       Else
          ye = Min(y, wv_d)
       End If
       kd = wv_k(i)*wv_d
       ! ratios written in exponentials-safe form
       ch = Cosh(wv_k(i)*ye)/Sinh(kd)
       sh = Sinh(wv_k(i)*ye)/Sinh(kd)
       u = u + wv_a(i)*wv_w(i)*ch*Cos(th)
       v = v + wv_a(i)*wv_w(i)*sh*Sin(th)
    End Do

  End Subroutine wave_vel_raw


  !> cosh(a)/cosh(b) and sinh(a)/cosh(b) for 0 <= a <= b without overflow
  Function ratio_cosh(a, b) Result(r)

    Real(Int64), Intent(In) :: a, b
    Real(Int64) :: r

    r = ( Exp(a - b) + Exp(-a - b) )/( 1d0 + Exp(-2d0*b) )

  End Function ratio_cosh


  Function ratio_sinh(a, b) Result(r)

    Real(Int64), Intent(In) :: a, b
    Real(Int64) :: r

    r = ( Exp(a - b) - Exp(-a - b) )/( 1d0 + Exp(-2d0*b) )

  End Function ratio_sinh


  !> Residuals of the Fenton Fourier stream-function problem for z = (k, c, B0, psi_s, R, B_1..B_N, s_0..s_N), height h
  Subroutine sf_residual(z, n, h, res)

    Integer(Int32), Intent(In)  :: n
    Real(Int64),    Intent(In)  :: z(2*n+6), h
    Real(Int64),    Intent(Out) :: res(2*n+6)

    Integer(Int32) :: m, j
    Real(Int64) :: k, c, b0, psis, rr, bj(n), sm(0:n), pi_, psi, um, vm, sk, sjm, ch, sh, ang, tsum

    pi_ = 4d0*Atan(1d0)
    k = z(1);  c = z(2);  b0 = z(3);  psis = z(4);  rr = z(5)
    bj = z(6:5+n)
    sm = z(6+n:6+2*n)
    Do m = 0, n
       psi = b0*sm(m)
       um = b0
       vm = 0d0
       Do j = 1, n
          ang = Real(j*m,Int64)*pi_/n
          ch = ratio_cosh(j*k*sm(m), j*k*wv_d)
          sh = ratio_sinh(j*k*sm(m), j*k*wv_d)
          psi = psi + bj(j)*sh*Cos(ang)
          um = um + j*k*bj(j)*ch*Cos(ang)
          vm = vm + j*k*bj(j)*sh*Sin(ang)
       End Do
       res(m+1) = psi - psis
       res(n+2+m) = 0.5d0*(um*um + vm*vm) + wv_g*sm(m) - rr
    End Do
    tsum = 0.5d0*(sm(0) + sm(n))
    Do m = 1, n-1
       tsum = tsum + sm(m)
    End Do
    res(2*n+3) = tsum/n - wv_d
    res(2*n+4) = sm(0) - sm(n) - h
    res(2*n+5) = k*c - 2d0*pi_/wave_period
    If ( wave_current_mode == 2 ) Then
       res(2*n+6) = b0 + c
    Else
       res(2*n+6) = psis + c*wv_d
    End If

  End Subroutine sf_residual


  !> Newton solve with a finite-difference Jacobian, continuing in wave height from a small wave to wave_height
  Subroutine stream_function_solve

    Integer(Int32) :: n, nu, it, ih, nh, i, info
    Real(Int64), Allocatable :: z(:), res(:), res2(:), jac(:,:), dz(:)
    Integer(Int32), Allocatable :: piv(:)
    Real(Int64) :: h, k0, om, pi_, a, eps, zs, rn
    Integer(Int32) :: j, m

    pi_ = 4d0*Atan(1d0)
    n = Min(Max(wave_sf_n, 4), SF_MAX_N)
    nu = 2*n + 6
    Allocate( z(nu), res(nu), res2(nu), jac(nu,nu), dz(nu), piv(nu) )
    om = 2d0*pi_/wave_period
    k0 = wave_number(om, wv_d, wv_g)
    nh = 8
    a = 0.5d0*wave_height/nh
    ! linear-wave start for the smallest height
    z(1) = k0;  z(2) = om/k0;  z(3) = -z(2);  z(4) = -z(2)*wv_d;  z(5) = 0.5d0*z(2)**2 + wv_g*wv_d
    z(6:5+n) = 0d0
    z(6) = a*om/(k0*Tanh(k0*wv_d))
    Do m = 0, n
       z(6+n+m) = wv_d + a*Cos(m*pi_/n)
    End Do
    Do ih = 1, nh
       h = wave_height*ih/nh
       Do it = 1, 60
          Call sf_residual(z, n, h, res)
          rn = Maxval(Abs(res))
          If ( rn < 1d-12*Max(1d0, wv_g*wv_d) ) Exit
          Do i = 1, nu
             eps = 1d-7*Max(Abs(z(i)), 1d-2)
             zs = z(i);  z(i) = zs + eps
             Call sf_residual(z, n, h, res2)
             z(i) = zs
             jac(:,i) = (res2 - res)/eps
          End Do
          dz = -res
          Call dgesv(nu, 1, jac, nu, piv, dz, nu, info)
          If ( info /= 0 ) Exit
          z = z + dz
       End Do
    End Do
    Call sf_residual(z, n, wave_height, res)
    If ( Maxval(Abs(res)) > 1d-8*Max(1d0, wv_g*wv_d) ) Then
       Write(*,'(A,ES10.2)') ' WARNING: stream-function solve did not converge, residual ', Maxval(Abs(res))
    End If

    sf_n = n
    sf_k = z(1);  sf_c = z(2);  sf_b(0) = z(3);  sf_q = z(4);  sf_r = z(5)
    sf_b(1:n) = z(6:5+n)
    ! surface cosine series, relative to the still level for j = 0
    Do j = 0, n
       a = 0.5d0*z(6+n) + 0.5d0*(-1d0)**j*z(6+2*n)
       Do m = 1, n-1
          a = a + z(6+n+m)*Cos(j*m*pi_/n)
       End Do
       sf_e(j) = a*Merge(1d0, 2d0, j == 0)/n
    End Do
    sf_e(0) = sf_e(0) - wv_d
    wv_n = 0
    Deallocate( z, res, res2, jac, dz, piv )

  End Subroutine stream_function_solve


  !> Full target at the inlet plane x = 0 on the solver's grid: eta, u at the cell-centre heights yg, v at the face heights y.
  !  Water below the surface follows the wave; above it the air returns the water flux uniformly.
  Subroutine wave_inlet_profile(t)

    Real(Int64), Intent(In) :: t
    Call wave_profile_at(0d0, t, wv_in_eta, wv_in_u, wv_in_v)

  End Subroutine wave_inlet_profile


  !> Same as the inlet profile at an arbitrary station x
  Subroutine wave_profile_at(x, t, eta, uc, vf)

    Real(Int64), Intent(In)  :: x, t
    Real(Int64), Intent(Out) :: eta, uc(nyg), vf(ny)
    Integer(Int32), Parameter :: nq = 80
    Integer(Int32) :: j, m
    Real(Int64) :: q, ys, u, v, hs, dys, ua, cf
    Integer(Int32) :: jj

    eta = wave_eta(x, t)
    hs = wv_d + eta
    ! water volume flux per unit span, trapezoid in y
    q = 0d0
    dys = hs/nq
    Do m = 0, nq
       ys = m*dys
       Call wave_vel(x, ys, t, u, v)
       q = q + Merge(0.5d0, 1d0, m == 0 .Or. m == nq)*u*dys
    End Do
    ua = -q/Max(Ly_i - hs, 1d-12)
    ! cell-averaged inlet velocity: the water part of each cell moves with the wave, the air part with the return flow, so that
    ! the face velocities integrate to the same zero net flux as the analytic profile (no flux mismatch in the surface cell)
    Do j = 1, nyg
       jj = Min(Max(j, 2), nyg-1)
       cf = Min(1d0, Max(0d0, (hs - y(jj-1))/(y(jj) - y(jj-1))))
       Call wave_vel(x, Max(Min(yg(jj), hs), 0d0), t, u, v)
       uc(j) = cf*u + (1d0 - cf)*ua
    End Do
    Do j = 1, ny
       If ( y(j) <= hs ) Then
          Call wave_vel(x, Max(y(j), 0d0), t, u, v)
          vf(j) = v
       Else
          vf(j) = 0d0
       End If
    End Do

  End Subroutine wave_profile_at

End Module waves
