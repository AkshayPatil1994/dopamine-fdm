!> Wave generation for the two-fluid solver: target free-surface elevation and velocity of a long-crested wave
!  in +x over the still-water depth d = vof_level (bed at y = 0), for the Dirichlet inlet and the relaxation
!  zones. wave_type 1 linear (Airy), 2 Rienecker-Fenton stream function, 3 JONSWAP (wave_Hs, wave_Tp,
!  wave_gamma) as wave_nfreq linear components with seeded random phases, finite-depth dispersion and Wheeler
!  stretching; 4 an elevation record eta(t) (wave_series_file) FFT-decomposed into linear components
!  with finite-depth dispersion and Wheeler stretching, the record time being the simulation time. A mean
!  current (current_type) adds to the water velocity. Above the surface the water velocity is undefined;
!  the air carries the uniform return flow that zeroes the inlet flux (closed flume): u_air = -Q/(Ly - d - eta), Q
!  the water column's volume flux. Setup runs on every rank from the same inputs and
!  seed, so no communication is needed.
Module waves

  Use iso_fortran_env, Only : Int32, Int64
  Use global, Only : nyg, ny, y, yg, Ly_i, vof_level, grav, wave_type, wave_height, wave_period, wave_phase, wave_sf_n, &
                     wave_current_mode, wave_Hs, wave_Tp, wave_gamma, wave_nfreq, wave_seed, wave_ramp_time, dx, &
                     wave_series_file, wave_series_pad, wave_series_taper, wave_series_fmin, wave_series_fmax, &
                     current_type, current_U, current_z0, current_n

  Implicit None

  Integer(Int32), Parameter :: WAVE_MAX_COMP = 4096, SF_MAX_N = 64
  Integer(Int32) :: wv_n = 0                   ! number of linear components
  Real(Int64), Allocatable, Dimension(:) :: wv_a, wv_k, wv_w, wv_ph
  Real(Int64) :: wv_d = 0d0, wv_g = 9.81d0
  ! stream-function solution (wave_type 2): wave number, phase speed, Eulerian mean-current definition,
  ! coefficients
  Integer(Int32) :: sf_n = 0
  Real(Int64) :: sf_k = 0d0, sf_c = 0d0, sf_ubar = 0d0, sf_q = 0d0, sf_r = 0d0
  Real(Int64) :: sf_b(0:SF_MAX_N) = 0d0, sf_e(0:SF_MAX_N) = 0d0
  ! m, mv: water fraction (x ramp) at the u/w and v heights
  Real(Int64), Allocatable, Dimension(:) :: wv_in_u, wv_in_v, wv_in_m, wv_in_mv
  Real(Int64) :: wv_in_eta = 0d0
  ! time of the stored inlet profile: it is a pure function of t, so equal-t calls reuse it
  Real(Int64) :: wv_in_t = -1d300

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
    Real(Int64) :: om, fp, fmin, fmax, df, f, sigma, r, alpha_j, m0, sj, pi_
    Real(Int64), Allocatable :: sp(:)
    Integer(Int32), Allocatable :: seed(:)
    Integer(Int32) :: nseed

    pi_ = 4d0*Atan(1d0)
    wv_d = vof_level
    wv_g = grav
    Allocate( wv_in_u(nyg), wv_in_v(ny), wv_in_m(nyg), wv_in_mv(ny) )
    wv_in_u = 0d0;  wv_in_v = 0d0;  wv_in_m = 0d0;  wv_in_mv = 0d0
    !$acc enter data copyin(wv_in_u,wv_in_v,wv_in_m,wv_in_mv)

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
    Case(4)
       Call wave_series_init
    Case Default
       wv_n = 0
    End Select

  End Subroutine wave_init


  !> In-place radix-2 FFT (forward, e^{-i...}) of z(0:m-1), m a power of two
  Subroutine wave_fft(z, m)

    Integer(Int32), Intent(In) :: m
    Complex(Int64), Intent(InOut) :: z(0:m-1)
    Integer(Int32) :: i, j, l, bit, len
    Complex(Int64) :: w, wl, a, c
    Real(Int64) :: ang

    j = 0
    Do i = 1, m-1
       bit = m/2
       Do While ( IAnd(j, bit) /= 0 )
          j = IEor(j, bit);  bit = bit/2
       End Do
       j = IEor(j, bit)
       If ( i < j ) Then
          a = z(i);  z(i) = z(j);  z(j) = a
       End If
    End Do
    len = 2
    Do While ( len <= m )
       ang = -8d0*Atan(1d0)/len
       wl = Cmplx(Cos(ang), Sin(ang), Int64)
       Do i = 0, m-1, len
          w = (1d0, 0d0)
          Do l = 0, len/2-1
             a = z(i+l);  c = w*z(i+l+len/2)
             z(i+l) = a + c;  z(i+l+len/2) = a - c
             w = w*wl
          End Do
       End Do
       len = 2*len
    End Do

  End Subroutine wave_fft


  !> wave_type 4: decompose the elevation record eta(t) at x=0 by FFT into linear components
  !  (w = 2 pi f, k from the finite-depth dispersion relation, a cos(k x - w t + ph)), record time =
  !  simulation time; then write fields/wave_series_check.out comparing the
  !  production wave_eta_raw(0,t) with the input record. The mean of the record is removed.
  Subroutine wave_series_init

    Use mpi, Only : myid, MPI_COMM_WORLD

    Integer(Int32) :: n, ns, m, i, j, ios, nt, nk, iu, ipk
    Real(Int64) :: pi_, kg, t0, dts, rl, fr, mean, var_in, fmx, amax, lo, hi, thr, tot2, e, rms_e, rms_i, max_e, max_i, slope
    Real(Int64), Allocatable :: tr(:), er(:), xs(:), amp(:), phs(:)
    Complex(Int64), Allocatable :: z(:)
    Character(300) :: line
    Real(Int64) :: col(2)

    pi_ = 4d0*Atan(1d0)
    Open(newunit=iu, file=Trim(wave_series_file), status='old', action='read', iostat=ios)
    If ( ios /= 0 ) Call series_abort( 'cannot open wave_series_file' )
    n = 0
    Do
       Read(iu,'(A)',iostat=ios) line
       If ( ios /= 0 ) Exit
       line = Adjustl(line)
       If ( Len_trim(line) == 0 .Or. line(1:1) == '#' ) Cycle
       n = n + 1
    End Do
    If ( n < 8 ) Call series_abort( 'wave_series_file needs at least 8 data rows' )
    Allocate( tr(n), er(n) )
    Rewind(iu)
    i = 0
    Do
       Read(iu,'(A)',iostat=ios) line
       If ( ios /= 0 ) Exit
       line = Adjustl(line)
       If ( Len_trim(line) == 0 .Or. line(1:1) == '#' ) Cycle
       Read(line,*,iostat=ios) col
       If ( ios /= 0 ) Call series_abort( 'wave_series_file: a data row does not parse as two columns (t eta)' )
       i = i + 1
       tr(i) = col(1);  er(i) = col(2)
    End Do
    Close(iu)
    If ( Any(Abs(tr) > Huge(1d0)) .Or. Any(tr /= tr) .Or. Any(Abs(er) > Huge(1d0)) .Or. Any(er /= er) ) &
         Call series_abort( 'wave_series_file: non-finite time or elevation' )
    If ( Any(tr(2:n) <= tr(1:n-1)) ) Call series_abort( 'wave_series_file: t must be strictly increasing' )

    ! uniform resampling at the finest input step
    t0 = tr(1)
    rl = tr(n) - t0
    ns = Min( Nint(rl/Minval(tr(2:n) - tr(1:n-1))) + 1, 2**22 )
    dts = rl/(ns - 1)
    Allocate( xs(ns) )
    j = 1
    Do i = 1, ns
       Do While ( j < n-1 .And. tr(j+1) < t0 + (i-1)*dts )
          j = j + 1
       End Do
       fr = Min(1d0, Max(0d0, (t0 + (i-1)*dts - tr(j))/(tr(j+1) - tr(j))))
       xs(i) = er(j) + fr*(er(j+1) - er(j))
    End Do
    mean = Sum(xs)/ns
    xs = xs - mean
    var_in = Sum(xs**2)/ns
    nt = Nint(wave_series_taper*ns)
    Do i = 1, nt
       fr = 0.5d0*(1d0 - Cos(pi_*(i - 1)/nt))
       xs(i) = fr*xs(i);  xs(ns+1-i) = fr*xs(ns+1-i)
    End Do

    If ( ns*(1d0 + wave_series_pad) > 2d0**26 ) Call series_abort( 'wave_series_pad makes the padded record exceed 2^26 samples' )
    m = 1
    Do While ( m < ns*(1d0 + wave_series_pad) )
       m = 2*m
    End Do
    Allocate( z(0:m-1) )
    z = (0d0, 0d0)
    z(0:ns-1) = xs
    Call wave_fft(z, m)

    ! bins k = 1..m/2-1 at f = k/(m dts); band, amplitude floor (1e-3 of the peak) and component cap
    Allocate( amp(m/2-1), phs(m/2-1) )
    amp = 2d0*Abs(z(1:m/2-1))/m
    phs = [( Atan2(Aimag(z(i)), Real(z(i))), i = 1, m/2-1 )]
    tot2 = Sum(amp**2)
    fmx = wave_series_fmax
    If ( fmx <= 0d0 ) Then
       ! shortest wavelength of 20 cells, finite-depth dispersion (same relation as the components)
       kg = 2d0*pi_/(20d0*Max(dx, 1d-300))
       fmx = Merge( Sqrt(wv_g*kg*Tanh(kg*wv_d))/(2d0*pi_), Huge(1d0), dx > 0d0 )
    End If
    Do i = 1, m/2-1
       fr = i/(m*dts)
       If ( fr < wave_series_fmin .Or. fr > fmx ) amp(i) = 0d0
    End Do
    ipk = Maxloc(amp, 1)
    amax = Maxval(amp)
    If ( amax <= 0d0 ) Call series_abort( 'no spectral content inside [wave_series_fmin, wave_series_fmax]' )
    lo = 1d-3*amax;  hi = amax
    If ( Count(amp >= lo) > WAVE_MAX_COMP ) Then
       Do i = 1, 60
          thr = 0.5d0*(lo + hi)
          If ( Count(amp >= thr) > WAVE_MAX_COMP ) Then
             lo = thr
          Else
             hi = thr
          End If
       End Do
       lo = hi
    End If
    nk = Count(amp >= lo)
    wv_n = nk
    Allocate( wv_a(nk), wv_k(nk), wv_w(nk), wv_ph(nk) )
    j = 0
    Do i = 1, m/2-1
       If ( amp(i) < lo ) Cycle
       j = j + 1
       wv_w(j)  = 2d0*pi_*i/(m*dts)
       wv_k(j)  = wave_number(wv_w(j), wv_d, wv_g)
       wv_a(j)  = amp(i)
       wv_ph(j) = wv_w(j)*t0 - phs(i)
    End Do

    ! largest surface slope |d eta/dx| of the linear sum over the record, which no wave of this form survives beyond ~1
    slope = 0d0
    Do i = 1, n
       slope = Max(slope, Abs(Sum( wv_a*wv_k*Sin(wv_ph - wv_w*tr(i)) )))
    End Do
    If ( slope >= 1d0 ) Then
       Write(*,'(A,F6.2)') ' ERROR: wave_series max surface slope = ', slope
       Write(*,'(A)') ' ERROR: the record is steeper than a linear sum of waves can represent (slope >= 1), breaking at the inlet.'
       Write(*,'(A)') ' ERROR: a sharp record start/end leaks high frequencies: use wave_series_taper or limit the band.'
       Write(*,'(A)') ' ERROR: reduce the amplitudes or restrict the band (wave_series_fmax).'
       Call MPI_Abort(MPI_COMM_WORLD, 1, ios)
    End If
    If ( vof_level + Maxval(er) - mean >= Ly_i .Or. vof_level + Minval(er) - mean <= 0d0 ) &
         Call series_abort( 'the record leaves the domain: vof_level + eta must stay inside (0, Ly)' )
    Deallocate( xs, z, amp, phs )
    If ( myid /= 0 ) Then
       Deallocate( tr, er )
       Return
    End If
    If ( nk > 1000 ) Write(*,'(A)') ' WARNING: wave_series has > 1000 components: the relaxation zones cost O(components) per cell'

    Write(*,'(A,I7,A,I7,A,ES10.3,A,F9.3,A)') ' wave_series: ', n, ' samples -> ', nk, ' components, dt =', dts, &
         ', record ', rl, ' s'
    Write(*,'(A,F8.4,A,F8.4,A,F6.2,A)') '   peak f =', ipk/(m*dts), ' Hz, band up to', Min(fmx, 1d6), ' Hz, max slope =', slope, ''
    Write(*,'(A,F6.3)') '   variance fraction kept in components = ', Sum(wv_a**2)/tot2
    Write(*,'(A,F6.3)') '   variance fraction with k*vof_level < pi (shallower than deep water) = ', &
         Sum( Merge(wv_a**2, 0d0, wv_k*wv_d < pi_) )/Sum(wv_a**2)
    If ( slope > 0.5d0 ) Write(*,'(A)') ' WARNING: wave_series is steep (max slope > 0.5); linear superposition is poor'
    If ( wave_ramp_time <= 0d0 .And. Abs(er(1) - mean) > 0.05d0*Maxval(Abs(er - mean)) ) &
         Write(*,'(A)') ' WARNING: wave_series starts away from rest and wave_ramp_time = 0: an abrupt start-up'

    Call execute_command_line('mkdir -p fields')
    Open(newunit=iu, file='fields/wave_series_check.out', status='replace')
    rms_e = 0d0;  rms_i = 0d0;  max_e = 0d0;  max_i = 0d0
    Do i = 1, n
       e = wave_eta_raw(0d0, tr(i)) - (er(i) - mean)
       rms_e = rms_e + e*e;  rms_i = rms_i + (er(i) - mean)**2
       max_e = Max(max_e, Abs(e));  max_i = Max(max_i, Abs(er(i) - mean))
    End Do
    rms_e = Sqrt(rms_e/n);  rms_i = Sqrt(rms_i/n)
    Write(iu,'(A)') '# wave_type 4: input record vs the linear components imposed at the inlet (x=0, no start-up ramp)'
    Write(iu,'(A,I8,A,ES12.4,A,ES12.4,A,I6,A,ES12.4)') '# samples', n, '  dt_resampled', dts, '  padded length', m*dts, &
         '  components', nk, '  fmax', Min(fmx, 1d6)
    Write(iu,'(A,ES12.4,A,ES12.4)') '# Hm0 input (4 sigma)', 4d0*Sqrt(var_in), '  Hm0 reconstruction', &
         4d0*Sqrt(0.5d0*Sum(wv_a**2)*m/ns)
    Write(iu,'(A,ES12.4,A,ES12.4)') '# rms error / rms input', rms_e/rms_i, '  max error / max input', max_e/max_i
    Write(iu,'(A)') '# columns: t  eta_input_minus_mean  eta_reconstructed  error'
    Do i = 1, n
       e = wave_eta_raw(0d0, tr(i))
       Write(iu,'(4ES16.8)') tr(i), er(i) - mean, e, e - (er(i) - mean)
    End Do
    Close(iu)
    Write(*,'(A,ES10.3,A,ES10.3,A)') '   reconstruction error: rms/rms_in =', rms_e/rms_i, ', max/max_in =', max_e/max_i, &
         '  (fields/wave_series_check.out)'
    If ( rms_e/rms_i > 0.05d0 ) Write(*,'(A)') ' WARNING: wave_series reconstruction error > 5%: band, component cap or taper '// &
         'removes content (see fields/wave_series_check.out)'
    Deallocate( tr, er )

  Contains

    Subroutine series_abort(msg)

      Character(*), Intent(In) :: msg

      Write(*,'(A,A)') ' ERROR: &WAVES ', msg
      Call MPI_Abort(MPI_COMM_WORLD, 1, ios)

    End Subroutine series_abort

  End Subroutine wave_series_init


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


  !> Water velocity (u, v) at the point (x, y), y from the bed, for y up to the surface (the caller handles
  !  the air)
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
       If ( wave_type >= 3 ) Then
          ! Wheeler stretching: map [0, d+eta] onto [0, d]; above the surface keep the surface value
          ye = Min(y, wv_d + eta)*wv_d/(wv_d + eta)
       Else
          ye = Min(y, wv_d)
       End If
       kd = wv_k(i)*wv_d
       ! cosh(k ye)/sinh(kd) and sinh(k ye)/sinh(kd) without overflow for large kd
       ch = ratio_cosh(wv_k(i)*ye, kd)/Tanh(kd)
       sh = ratio_sinh(wv_k(i)*ye, kd)/Tanh(kd)
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


  !> Residuals of the Fenton Fourier stream-function problem for z = (k, c, B0, psi_s, R, B_1..B_N, s_0..s_N),
  !  height h
  Subroutine sf_residual(z, n, h, res)

    Integer(Int32), Intent(In)  :: n
    Real(Int64),    Intent(In)  :: z(2*n+6), h
    Real(Int64),    Intent(Out) :: res(2*n+6)

    Integer(Int32) :: m, j
    Real(Int64) :: k, c, b0, psis, rr, bj(n), sm(0:n), pi_, psi, um, vm, ch, sh, ang, tsum

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


  !> Newton solve with a finite-difference Jacobian, continuing in wave height from a small wave to
  !  wave_height
  Subroutine stream_function_solve

    Use mpi
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
       Write(*,'(A,ES10.2)') ' ERROR: stream-function solve did not converge, residual ', Maxval(Abs(res))
       Write(*,'(A,F7.4,A,F7.4,A)') ' ERROR: wave steepness H/L = ', wave_height*k0/(2d0*pi_), &
            ' vs breaking limit ~', 0.142d0*Tanh(k0*wv_d), ' (0.142 tanh(kd), linear k)'
       Write(*,'(A)') ' ERROR: the requested wave (wave_height, wave_period, vof_level) is at or beyond the breaking limit, '// &
            'so no periodic stream-function wave exists. Reduce wave_height, increase wave_period or increase vof_level.'
       Call MPI_Abort(MPI_COMM_WORLD, 1, i)
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


  !> Full target at the inlet plane x = 0 on the solver's grid: eta, u at the cell-centre heights yg, v at the
  !  face heights y. Water below the surface follows the wave; above it the air returns the water flux
  !  uniformly.
  Subroutine wave_inlet_profile(t)

    Real(Int64), Intent(In) :: t
    Real(Int64) :: hs, rm
    Integer(Int32) :: j, jj
    If ( t == wv_in_t ) Return
    Call wave_profile_at(0d0, t, wv_in_eta, wv_in_u, wv_in_v)
    wv_in_t = t
    hs = wv_d + wv_in_eta
    rm = wave_ramp(t)
    Do j = 1, nyg
       jj = Min(Max(j, 2), nyg-1)
       wv_in_m(j) = rm*Min(1d0, Max(0d0, (hs - y(jj-1))/(y(jj) - y(jj-1))))
    End Do
    Do j = 1, ny
       wv_in_mv(j) = Merge(rm, 0d0, y(j) <= hs)
    End Do
    !$acc update device(wv_in_u,wv_in_v,wv_in_m,wv_in_mv)

  End Subroutine wave_inlet_profile


  !> Mean current at height ys (clamped to the still-water depth, so the water above it keeps the surface
  !  value); depth-mean current_U
  Function current_vel(ys, t) Result(u)

    Real(Int64), Intent(In) :: ys, t
    Real(Int64) :: u, yc

    u = 0d0
    yc = Min(Max(ys, 0d0), wv_d)
    Select Case(current_type)
    Case(1)
       u = current_U
    Case(2)
       u = current_U*Log(Max(yc, current_z0)/current_z0)/( Log(wv_d/current_z0) - 1d0 + current_z0/wv_d )
    Case(3)
       u = current_U*(current_n + 1d0)/current_n*(yc/wv_d)**(1d0/current_n)
    End Select
    u = wave_ramp(t)*u

  End Function current_vel


  !> Cell-centre velocity of the current alone (water below the still level, nothing above): the target of the
  !  absorption zone
  Subroutine current_profile(t, uc)

    Real(Int64), Intent(In)  :: t
    Real(Int64), Intent(Out) :: uc(nyg)
    Integer(Int32) :: j, jj
    Real(Int64) :: cf

    Do j = 1, nyg
       jj = Min(Max(j, 2), nyg-1)
       cf = Min(1d0, Max(0d0, (wv_d - y(jj-1))/(y(jj) - y(jj-1))))
       uc(j) = cf*current_vel(yg(jj), t)
    End Do

  End Subroutine current_profile


  !> Same as the inlet profile at an arbitrary station x
  Subroutine wave_profile_at(x, t, eta, uc, vf)

    Real(Int64), Intent(In)  :: x, t
    Real(Int64), Intent(Out) :: eta, uc(nyg), vf(ny)
    Integer(Int32), Parameter :: nq = 80
    Integer(Int32) :: j, m
    Real(Int64) :: q, ys, u, v, hs, dys, ua, cf
    Integer(Int32) :: jj
    Real(Int64) :: eraw, cj(SF_MAX_N), sj(SF_MAX_N), th

    eta = wave_eta(x, t)
    hs = wv_d + eta
    ! stream function: the harmonics in x do not depend on the height, evaluate them once per station
    If ( wave_type == 2 ) Then
       eraw = wave_eta_raw(x, t)
       th = sf_k*x - 2d0*4d0*Atan(1d0)*t/wave_period + wave_phase
       Do j = 1, sf_n
          cj(j) = Cos(Real(j,Int64)*th);  sj(j) = Sin(Real(j,Int64)*th)
       End Do
    End If
    ! water volume flux per unit span, trapezoid in y
    q = 0d0
    dys = hs/nq
    Do m = 0, nq
       ys = m*dys
       Call vel_at(ys)
       q = q + Merge(0.5d0, 1d0, m == 0 .Or. m == nq)*u*dys
    End Do
    ua = -q/Max(Ly_i - hs, 1d-12)
    ! cell-averaged inlet velocity: the water part of each cell moves with the wave, the air part with the
    ! return flow, so that the face velocities integrate to the same zero net flux as the analytic profile
    ! (no flux mismatch in the surface cell)
    Do j = 1, nyg
       jj = Min(Max(j, 2), nyg-1)
       cf = Min(1d0, Max(0d0, (hs - y(jj-1))/(y(jj) - y(jj-1))))
       Call vel_at(Max(Min(yg(jj), hs), 0d0))
       uc(j) = cf*( u + current_vel(yg(jj), t) ) + (1d0 - cf)*ua
    End Do
    Do j = 1, ny
       If ( y(j) <= hs ) Then
          Call vel_at(Max(y(j), 0d0))
          vf(j) = v
       Else
          vf(j) = 0d0
       End If
    End Do

  Contains

    !> wave_vel at height yy; the stream-function exponentials are powers of one exponential per height
    !  (equal to the direct evaluation to round-off)
    Subroutine vel_at(yy)

      Real(Int64), Intent(In) :: yy
      Integer(Int32) :: m2
      Real(Int64) :: ye, pe, qe, re, pj, qj, rj, w

      If ( wave_type /= 2 ) Then
         Call wave_vel(x, yy, t, u, v)
         Return
      End If
      ye = Min(yy, wv_d + eraw)
      pe = Exp(sf_k*(ye - wv_d));  qe = Exp(-sf_k*(ye + wv_d));  re = Exp(-2d0*sf_k*wv_d)
      pj = 1d0;  qj = 1d0;  rj = 1d0
      u = sf_c + sf_b(0);  v = 0d0
      Do m2 = 1, sf_n
         pj = pj*pe;  qj = qj*qe;  rj = rj*re
         w = Real(m2,Int64)*sf_k*sf_b(m2)/( 1d0 + rj )
         u = u + w*(pj + qj)*cj(m2)
         v = v + w*(pj - qj)*sj(m2)
      End Do
      u = wave_ramp(t)*u;  v = wave_ramp(t)*v

    End Subroutine vel_at

  End Subroutine wave_profile_at

End Module waves
