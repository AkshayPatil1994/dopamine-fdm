!> Standalone check of the momentum face-reconstruction kernels (src/mom_recon.f90): periodic linear advection u_t + u_x = 0 of point
!  values, SSP-RK3 at CFL 0.1 so the error is the spatial one. Reports (A) amplitude and phase error after one period per
!  wavenumber k dx, (B) observed order on a smooth profile, (C) overshoot / undershoot on a step.
Program test_mom_kernels

  Use iso_fortran_env, Only : Int32, Int64
  Use mom_recon
  Implicit None

  Real(Int64), Parameter :: pi = 3.14159265358979323846d0
  Integer(Int32), Parameter :: nsch = 9
  Integer(Int32) :: schemes(nsch) = (/ 0, 5, 4, 8, 2, 3, 6, 7, 1 /)
  Character(len=12) :: names(nsch) = (/ 'central2    ', 'central4    ', 'QUICK       ', 'QUICK-clip  ', 'Koren       ', &
                                         'WENO3-Z-clip', 'WENO5-Z     ', 'WENO5-Z-clip', 'upwind1     ' /)
  Integer(Int32) :: is, im, n, ms(7), k
  Real(Int64) :: amp, ph, e32, e64, e128, mn, mx, l1
  Real(Int64), Allocatable :: u(:)

  ms = (/ 2, 4, 8, 12, 16, 20, 24 /)
  Write(*,'(A)') '(A) one period, N=64, CFL 0.1: amplitude ratio / phase error [deg] per k dx = 2 pi m / 64'
  Write(*,'(A12,7(5X,A,F4.2))') 'scheme', ('kdx=', 2d0*pi*ms(k)/64d0, k=1,7)
  Do is = 1, nsch
     Write(*,'(A12)',advance='no') names(is)
     Do im = 1, 7
        n = 64
        Allocate( u(n) )
        Do k = 1, n
           u(k) = 0.01d0*Sin( 2d0*pi*ms(im)*(k-1)/n )
        End Do
        Call advect(u, n, 1d0, schemes(is), 0.1d0)
        Call project(u, n, ms(im), amp, ph)
        Write(*,'(F8.4,A,F7.2)',advance='no') amp/0.01d0, '/', ph
        Deallocate( u )
     End Do
     Write(*,*)
  End Do

  Write(*,'(/,A)') '(B) smooth profile sin(2 pi x)+0.5 sin(6 pi x): L2 error after one period at N=32/64/128 and observed order'
  Do is = 1, nsch
     e32 = smooth_err(32, schemes(is));  e64 = smooth_err(64, schemes(is));  e128 = smooth_err(128, schemes(is))
     Write(*,'(A12,3ES11.3,2F7.2)') names(is), e32, e64, e128, Log(e32/e64)/Log(2d0), Log(e64/e128)/Log(2d0)
  End Do

  Write(*,'(/,A)') '(C) square wave (N=100, one period): min, max, L1 error'
  Do is = 1, nsch
     n = 100
     Allocate( u(n) )
     Do k = 1, n
        u(k) = Merge(1d0, 0d0, (k-1) >= 25 .And. (k-1) < 75)
     End Do
     Call advect(u, n, 1d0, schemes(is), 0.1d0)
     mn = Minval(u);  mx = Maxval(u)
     l1 = 0d0
     Do k = 1, n
        l1 = l1 + Abs( u(k) - Merge(1d0, 0d0, (k-1) >= 25 .And. (k-1) < 75) )/n
     End Do
     Write(*,'(A12,3F10.5)') names(is), mn, mx, l1
     Deallocate( u )
  End Do

  ! (D) the vectorizable row kernel against mom_face + the refill-Courant blend, both flux signs and blend weights 0, (0,1), 1
  Block
    Integer(Int32), Parameter :: nr = 4000
    Real(Int64) :: pad(nr+5), mda(nr), r1(nr), r2(nr), f(nr), sm(-2:3), mf, cc, up, tb, ref, err, cm0, cm1, ivcv
    cm0 = 2d-3;  cm1 = 1d-2;  ivcv = 1d3
    Call random_number(pad);  Call random_number(mda);  Call random_number(r1);  Call random_number(r2)
    pad = pad + Sin( 0.05d0*[(Real(k,Int64), k = 1, nr+5)] )
    mda = mda - 0.5d0
    r1 = 1d0 + 999d0*Merge(1d0, 0d0, r1 > 0.9d0);  r2 = r1 + Merge(0d0, 998d0, r2 < 0.95d0)
    Call mom_flux_weno5z_row(nr, pad(1:nr), pad(2:nr+1), pad(3:nr+2), pad(4:nr+3), pad(5:nr+4), pad(6:nr+5), &
                             mda, mda, r1, r2, ivcv, cm0, cm1, f)
    err = 0d0
    Do k = 1, nr
       sm = pad(k:k+5);  mf = mda(k)
       cc = Abs(mf)*Abs(1d0/r1(k) - 1d0/r2(k))*ivcv
       up = Merge(sm(0), sm(1), mf >= 0d0)
       tb = Min(1d0, Max(0d0, (cc - cm0)/(cm1 - cm0)))
       ref = mf*( up + (1d0 - tb*tb*(3d0 - 2d0*tb))*( mom_face(sm, mf, MOM_WENO5) - up ) )
       err = Max(err, Abs(f(k) - ref)/Max(Abs(ref), 1d-3*Abs(mf)))
    End Do
    Write(*,'(/,A,ES10.2)') '(D) row kernel vs mom_face + blend, max relative difference: ', err
    If ( err > 1d-12 ) Then
       Write(*,'(A)') 'FAIL'
       Stop 1
    End If
  End Block

Contains

  Subroutine advect(u, n, tend, sch, cfl)

    Integer(Int32), Intent(In)    :: n, sch
    Real(Int64),    Intent(InOut):: u(n)
    Real(Int64),    Intent(In)    :: tend, cfl
    Real(Int64) :: dx, dt, t, u1(n), u2(n)
    Integer(Int32) :: nst, s

    dx = 1d0/n;  dt = cfl*dx;  nst = Nint(tend/dt);  dt = tend/nst
    Do s = 1, nst
       u1 = u + dt*rhs(u, n, sch)
       u2 = 0.75d0*u + 0.25d0*( u1 + dt*rhs(u1, n, sch) )
       u  = u/3d0 + 2d0/3d0*( u2 + dt*rhs(u2, n, sch) )
    End Do

  End Subroutine advect

  Function rhs(u, n, sch) Result(r)

    Integer(Int32), Intent(In) :: n, sch
    Real(Int64),    Intent(In) :: u(n)
    Real(Int64) :: r(n), f(0:n), st(-2:3)
    Integer(Int32) :: i, m, j

    Do i = 0, n          ! face between cell i and i+1 (periodic)
       Do m = -2, 3
          j = Modulo( i + m - 1, n ) + 1
          st(m) = u(j)
       End Do
       f(i) = mom_face(st, 1d0, sch)
    End Do
    Do i = 1, n
       r(i) = -( f(i) - f(i-1) )*n
    End Do

  End Function rhs

  Subroutine project(u, n, m, amp, ph)

    Integer(Int32), Intent(In)  :: n, m
    Real(Int64),    Intent(In)  :: u(n)
    Real(Int64),    Intent(Out) :: amp, ph
    Real(Int64) :: a, b
    Integer(Int32) :: k

    a = 0d0;  b = 0d0
    Do k = 1, n
       a = a + 2d0/n*u(k)*Sin( 2d0*pi*m*(k-1)/n )
       b = b + 2d0/n*u(k)*Cos( 2d0*pi*m*(k-1)/n )
    End Do
    amp = Sqrt(a*a + b*b)
    ph = Atan2(b, a)*180d0/pi      ! exact solution after one period: phase 0
  End Subroutine project

  Function smooth_err(n, sch) Result(e)

    Integer(Int32), Intent(In) :: n, sch
    Real(Int64) :: e
    Real(Int64), Allocatable :: u(:), u0(:)
    Integer(Int32) :: k

    Allocate( u(n), u0(n) )
    Do k = 1, n
       u0(k) = Sin( 2d0*pi*(k-1)/n ) + 0.5d0*Sin( 6d0*pi*(k-1)/n )
    End Do
    u = u0
    Call advect(u, n, 1d0, sch, 0.1d0)
    e = Sqrt( Sum( (u - u0)**2 )/n )
    Deallocate( u, u0 )
  End Function smooth_err

End Program test_mom_kernels
