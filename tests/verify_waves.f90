!> Standalone check of src/waves.f90: linear dispersion, the stream-function solution (kinematic and dynamic surface conditions,
!  height, mean level, zero mean flux) and the JONSWAP significant wave height. Exit code 1 on failure.
Program verify_waves

  Use iso_fortran_env, Only : Int32, Int64
  Use global, Only : nyg, ny, vof_level, grav, wave_type, wave_height, wave_period, wave_phase, wave_sf_n, wave_current_mode, &
                     wave_Hs, wave_Tp, wave_gamma, wave_nfreq, wave_seed, &
                     current_type, current_U, current_z0, current_n, Ly_i, wave_series_file, wave_series_taper, wave_series_fmax
  Use mpi, Only : myid
  Use waves

  Implicit None

  Real(Int64), Parameter :: pi = 3.14159265358979323846d0
  Real(Int64) :: d, g, x, eta, u, v, etax, h, res_k, res_b, bmin, bmax, bern, q, tt, ys, mean, var
  Real(Int64) :: ex, em, ux, uy, ts, uref, ys2, k1, k2, w1, w2, errmax, vbed
  Integer(Int32) :: iu
  Integer(Int32) :: i, j, nfail, nt, nx_s
  Real(Int64), Allocatable :: ser(:)

  nfail = 0
  nyg = 1;  ny = 1
  vof_level = 0.5d0;  grav = 9.81d0
  d = vof_level;  g = grav

  ! ---- linear wave
  wave_type = 1;  wave_height = 0.05d0;  wave_period = 1.5d0;  wave_phase = 0d0
  Call wave_init
  Write(*,'(A,2ES12.4)') 'linear: k, omega^2/(g k tanh(kd)) - 1 =', wv_k(1), wv_w(1)**2/(g*wv_k(1)*Tanh(wv_k(1)*d)) - 1d0
  If ( Abs(wave_eta(0d0, 0d0) - 0.025d0) > 1d-14 ) Then
     Write(*,*) 'FAIL linear crest';  nfail = nfail + 1
  End If
  Deallocate( wv_a, wv_k, wv_w, wv_ph, wv_in_u, wv_in_v, wv_in_m, wv_in_mv )

  ! ---- stream function, H/d = 0.4
  wave_type = 2;  wave_height = 0.2d0;  wave_period = 2.2d0;  wave_sf_n = 20;  wave_current_mode = 1
  Call wave_init
  Write(*,'(A,3ES14.6)') 'stream function: k, c, kd =', sf_k, sf_c, sf_k*d
  Write(*,'(A,2ES14.6)') '  crest, trough elevation =', wave_eta(0d0, 0d0), wave_eta(pi/sf_k, 0d0)
  If ( Abs(wave_eta(0d0,0d0) - wave_eta(pi/sf_k,0d0) - wave_height) > 1d-9 ) Then
     Write(*,*) 'FAIL height';  nfail = nfail + 1
  End If
  ! mean level over a wavelength
  mean = 0d0
  nx_s = 400
  Do i = 0, nx_s-1
     mean = mean + wave_eta(2d0*pi/sf_k*(i+0.5d0)/nx_s, 0d0)/nx_s
  End Do
  Write(*,'(A,ES12.3)') '  mean elevation =', mean
  If ( Abs(mean) > 1d-9 ) Then
     Write(*,*) 'FAIL mean level';  nfail = nfail + 1
  End If
  ! kinematic and dynamic conditions on the surface
  bmin = 1d300;  bmax = -1d300;  res_k = 0d0
  Do i = 0, 40
     x = 2d0*pi/sf_k*i/40d0
     eta = wave_eta(x, 0d0)
     etax = (wave_eta(x + 1d-6, 0d0) - wave_eta(x - 1d-6, 0d0))/2d-6
     Call wave_vel(x, d + eta, 0d0, u, v)
     res_k = Max(res_k, Abs((u - sf_c)*etax - v)/sf_c)
     bern = 0.5d0*((u - sf_c)**2 + v**2) + g*(d + eta)
     bmin = Min(bmin, bern);  bmax = Max(bmax, bern)
  End Do
  Write(*,'(A,2ES12.3)') '  kinematic residual/c, Bernoulli spread =', res_k, bmax - bmin
  If ( res_k > 1d-6 .Or. bmax - bmin > 1d-6 ) Then
     Write(*,*) 'FAIL surface conditions';  nfail = nfail + 1
  End If
  ! zero mean volume flux through a fixed section over one period
  q = 0d0
  nt = 400
  Do i = 0, nt-1
     tt = wave_period*(i + 0.5d0)/nt
     eta = wave_eta(0.3d0, tt)
     Do j = 0, 200
        ys = (d + eta)*(j + 0.5d0)/201d0
        Call wave_vel(0.3d0, ys, tt, u, v)
        q = q + u*(d + eta)/201d0/nt
     End Do
  End Do
  Write(*,'(A,ES12.3)') '  time-mean flux =', q
  If ( Abs(q) > 1d-5 ) Then
     Write(*,*) 'FAIL mean flux';  nfail = nfail + 1
  End If
  Deallocate( wv_in_u, wv_in_v, wv_in_m, wv_in_mv )

  ! ---- elevation series (wave_type 4): two tones at t0 /= 0 recovered by the FFT decomposition
  myid = 0;  Ly_i = 1d0
  w1 = 2d0*pi*0.8d0;  w2 = 2d0*pi*1.3d0;  k1 = wave_number(w1, d, g);  k2 = wave_number(w2, d, g)
  Open(newunit=iu, file='verify_waves_series.txt', status='replace')
  Write(iu,'(A)') '# t eta'
  Do i = 0, 799
     ts = 2d0 + 0.05d0*i
     Write(iu,'(2ES22.14)') ts, 0.02d0*Cos(w1*ts + 0.7d0) + 0.01d0*Cos(w2*ts - 0.4d0) + 0.003d0
  End Do
  Close(iu)
  wave_type = 4;  wave_series_file = 'verify_waves_series.txt'
  wave_ramp_time = 0d0;  wave_series_taper = 0d0;  wave_series_fmax = 4d0
  Call wave_init
  errmax = 0d0
  ! the band limit rings at the abrupt record ends only: check the interior
  Do i = 40, 760
     ts = 2d0 + 0.05d0*i
     errmax = Max(errmax, Abs(wave_eta(0d0, ts) - (0.02d0*Cos(w1*ts + 0.7d0) + 0.01d0*Cos(w2*ts - 0.4d0))))
  End Do
  Write(*,'(A,I5,A,ES12.4)') 'series: components', wv_n, ', max |eta_rec - eta_in| / H =', errmax/0.06d0
  If ( errmax > 0.03d0*0.06d0 ) Then
     Write(*,*) 'FAIL series reconstruction';  nfail = nfail + 1
  End If
  ! finite-depth velocity of the two tones (Wheeler stretching maps [0, d+eta] onto [0, d]): the surface (coth) and half depth
  ! (cosh ratio) horizontal velocity, and zero vertical velocity at the bed
  Do i = 1, 3
     ts = 10d0 + 0.37d0*i
     eta = wave_eta(0d0, ts)
     Call wave_vel(0d0, d + eta, ts, u, v)
     uref = 0.02d0*w1/Tanh(k1*d)*Cos(w1*ts + 0.7d0) + 0.01d0*w2/Tanh(k2*d)*Cos(w2*ts - 0.4d0)
     Write(*,'(A,2ES12.4)') '  surface u, analytic =', u, uref
     If ( Abs(u - uref) > 0.05d0*(0.02d0*w1 + 0.01d0*w2) ) Then
        Write(*,*) 'FAIL series surface velocity';  nfail = nfail + 1
     End If
     ys2 = 0.5d0*d*d/(d + eta)
     Call wave_vel(0d0, 0.5d0*d, ts, u, v)
     uref = 0.02d0*w1*Cosh(k1*ys2)/Sinh(k1*d)*Cos(w1*ts + 0.7d0) + 0.01d0*w2*Cosh(k2*ys2)/Sinh(k2*d)*Cos(w2*ts - 0.4d0)
     If ( Abs(u - uref) > 0.05d0*(0.02d0*w1 + 0.01d0*w2) ) Then
        Write(*,*) 'FAIL series half-depth velocity';  nfail = nfail + 1
     End If
     Call wave_vel(0d0, 0d0, ts, ux, vbed)
     If ( Abs(vbed) > 1d-12 ) Then
        Write(*,*) 'FAIL series v at the bed';  nfail = nfail + 1
     End If
  End Do
  Deallocate( wv_a, wv_k, wv_w, wv_ph, wv_in_u, wv_in_v, wv_in_m, wv_in_mv )

  ! ---- JONSWAP: significant height from the variance of a long record
  wave_type = 3;  wave_Hs = 0.1d0;  wave_Tp = 1.6d0;  wave_gamma = 3.3d0;  wave_nfreq = 400;  wave_seed = 7
  Call wave_init
  nt = 20000
  Allocate( ser(nt) )
  Do i = 1, nt
     ser(i) = wave_eta(0d0, 0.05d0*i)
  End Do
  mean = Sum(ser)/nt
  var = Sum((ser - mean)**2)/nt
  Write(*,'(A,2ES12.4)') 'JONSWAP: 4 sigma, Hs =', 4d0*Sqrt(var), wave_Hs
  If ( Abs(4d0*Sqrt(var)/wave_Hs - 1d0) > 0.1d0 ) Then
     Write(*,*) 'FAIL Hs';  nfail = nfail + 1
  End If

  ! ---- mean current: depth-mean equals current_U for every profile, zero at the bed for the log law
  current_U = 0.3d0;  current_z0 = 1d-3;  current_n = 7d0
  Do i = 1, 3
     current_type = i
     mean = 0d0
     nt = 20000
     Do j = 1, nt
        mean = mean + current_vel(d*(j - 0.5d0)/nt, 1d9)/nt
     End Do
     Write(*,'(A,I2,A,ES12.4)') 'current type', i, ': depth mean / current_U - 1 =', mean/current_U - 1d0
     If ( Abs(mean/current_U - 1d0) > 2d-3 ) Then
        Write(*,*) 'FAIL current mean';  nfail = nfail + 1
     End If
  End Do
  current_type = 2
  If ( Abs(current_vel(0d0, 1d9)) > 0d0 .Or. current_vel(d, 1d9) <= current_U ) Then
     Write(*,*) 'FAIL log profile shape';  nfail = nfail + 1
  End If

  If ( nfail > 0 ) Then
     Write(*,'(A,I3,A)') 'verify_waves: ', nfail, ' check(s) FAILED'
     Stop 1
  End If
  Write(*,'(A)') 'verify_waves: all checks passed'

End Program verify_waves
