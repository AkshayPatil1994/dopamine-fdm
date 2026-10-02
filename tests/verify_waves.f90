!> Standalone check of src/waves.f90: linear dispersion, the stream-function solution (kinematic and dynamic surface conditions,
!  height, mean level, zero mean flux) and the JONSWAP significant wave height. Exit code 1 on failure.
Program verify_waves

  Use iso_fortran_env, Only : Int32, Int64
  Use global, Only : nyg, ny, vof_level, vof_grav, wave_type, wave_height, wave_period, wave_phase, wave_sf_n, wave_current_mode, &
                     wave_Hs, wave_Tp, wave_gamma, wave_nfreq, wave_seed
  Use waves

  Implicit None

  Real(Int64), Parameter :: pi = 3.14159265358979323846d0
  Real(Int64) :: d, g, x, eta, u, v, etax, h, res_k, res_b, bmin, bmax, bern, q, tt, ys, mean, var
  Real(Int64) :: ex, em, ux, uy
  Integer(Int32) :: i, j, nfail, nt, nx_s
  Real(Int64), Allocatable :: ser(:)

  nfail = 0
  nyg = 1;  ny = 1
  vof_level = 0.5d0;  vof_grav = 9.81d0
  d = vof_level;  g = vof_grav

  ! ---- linear wave
  wave_type = 1;  wave_height = 0.05d0;  wave_period = 1.5d0;  wave_phase = 0d0
  Call wave_init
  Write(*,'(A,2ES12.4)') 'linear: k, omega^2/(g k tanh(kd)) - 1 =', wv_k(1), wv_w(1)**2/(g*wv_k(1)*Tanh(wv_k(1)*d)) - 1d0
  If ( Abs(wave_eta(0d0, 0d0) - 0.025d0) > 1d-14 ) Then
     Write(*,*) 'FAIL linear crest';  nfail = nfail + 1
  End If
  Deallocate( wv_a, wv_k, wv_w, wv_ph, wv_in_u, wv_in_v )

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
  Deallocate( wv_in_u, wv_in_v )

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

  If ( nfail > 0 ) Then
     Write(*,'(A,I3,A)') 'verify_waves: ', nfail, ' check(s) FAILED'
     Stop 1
  End If
  Write(*,'(A)') 'verify_waves: all checks passed'

End Program verify_waves
