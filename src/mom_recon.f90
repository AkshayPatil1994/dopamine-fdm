!> Face reconstruction kernels for the convective flux of the transported velocity. Pure functions with no solver dependencies so they
!  can be unit-tested (tests/vof/test_mom_kernels.f90) and shared between the staggered momentum update and any scalar transport.
!
!  Face between cells b1 (lower index) and c1 (upper index); stencil  a3 a2 | b1 c1 | c2 c3  ->  values (s(-2:3)) in index order
!  s(-2)=a3 s(-1)=a2 s(0)=b1 | s(1)=c1 s(2)=c2 s(3)=c3.  mf is the mass flux through the face: its sign selects the upwind side.
!
!  Schemes (mom_scheme):
!    0 central 2nd order            1 first-order upwind             2 Koren (bounded, TVD, 3rd order)
!    3 WENO3-Z clipped to stencil   4 QUICK (3rd order, unbounded)   5 central 4th order
!    6 WENO5-Z (5th order)          7 WENO5-Z clipped                8 QUICK clipped to [min,max] of the three nearest cells
!  "clipped" = result limited to the range of the three cells nearest the face (up, down, far-up): max-principle for the velocity.
Module mom_recon

  Use iso_fortran_env, Only : Int32, Int64

  Implicit None

  Integer(Int32), Parameter :: MOM_CENTRAL = 0, MOM_UPWIND = 1, MOM_KOREN = 2, MOM_WENO3 = 3, MOM_QUICK = 4, MOM_CENTRAL4 = 5, &
                               MOM_WENO5 = 6, MOM_WENO5C = 7, MOM_QUICKC = 8

Contains

  !> Face value for scheme sch from the six-point stencil s(-2:3) and the sign of the face flux mf
  Pure Function mom_face(s, mf, sch) Result(r)
    !$acc routine seq

    Real(Int64),    Intent(In) :: s(-2:3), mf
    Integer(Int32), Intent(In) :: sch
    Real(Int64) :: r

    Real(Int64) :: v(-2:2)   ! upwind-oriented stencil: v(0) = upwind cell, v(1) = downwind cell, v(-1), v(-2) further upwind, v(2) further downwind

    If ( sch == MOM_CENTRAL ) Then
       r = 0.5d0*( s(0) + s(1) );  Return
    End If
    If ( sch == MOM_CENTRAL4 ) Then
       r = ( -s(-1) + 7d0*s(0) + 7d0*s(1) - s(2) )/12d0;  Return
    End If
    If ( mf >= 0d0 ) Then
       v(-2) = s(-2);  v(-1) = s(-1);  v(0) = s(0);  v(1) = s(1);  v(2) = s(2)
    Else
       v(-2) = s(3);   v(-1) = s(2);   v(0) = s(1);  v(1) = s(0);  v(2) = s(-1)
    End If
    Select Case(sch)
    Case(MOM_UPWIND)
       r = v(0)
    Case(MOM_KOREN)
       r = koren(v(-1), v(0), v(1))
    Case(MOM_WENO3)
       r = clip3( weno3z(v(-1), v(0), v(1)), v(-1), v(0), v(1) )
    Case(MOM_QUICK)
       r = quick(v(-1), v(0), v(1))
    Case(MOM_QUICKC)
       r = clip3( quick(v(-1), v(0), v(1)), v(-1), v(0), v(1) )
    Case(MOM_WENO5)
       r = weno5z(v(-2), v(-1), v(0), v(1), v(2))
    Case(MOM_WENO5C)
       r = clip3( weno5z(v(-2), v(-1), v(0), v(1), v(2)), v(-1), v(0), v(1) )
    Case Default
       r = 0.5d0*( s(0) + s(1) )
    End Select

  End Function mom_face


  Pure Function clip3(r, a, b, c) Result(q)
    !$acc routine seq
    Real(Int64), Intent(In) :: r, a, b, c
    Real(Int64) :: q
    q = Min( Max(r, Min(a, b, c)), Max(a, b, c) )
  End Function clip3


  !> Koren limiter, far upwind f, upwind u, downwind d
  Pure Function koren(f, u, d) Result(r)
    !$acc routine seq
    Real(Int64), Intent(In) :: f, u, d
    Real(Int64) :: r, rr, ph

    r = u
    If ( d == u ) Return
    rr = ( u - f )/( d - u )
    ph = Max( 0d0, Min( 2d0*rr, ( 1d0 + 2d0*rr )/3d0, 2d0 ) )
    r = u + 0.5d0*ph*( d - u )
  End Function koren


  Pure Function quick(f, u, d) Result(r)
    !$acc routine seq
    Real(Int64), Intent(In) :: f, u, d
    Real(Int64) :: r
    r = ( -f + 6d0*u + 3d0*d )/8d0
  End Function quick


  Pure Function weno3z(f, u, d) Result(r)
    !$acc routine seq
    Real(Int64), Intent(In) :: f, u, d
    Real(Int64) :: r, s0, s1, tz, sc, w0, w1

    s0 = ( u - f )**2;  s1 = ( d - u )**2;  tz = Abs(s0 - s1)
    sc = 1d-12 + 1d-6*( f**2 + u**2 + d**2 )
    w0 = ( 1d0/3d0 )*( 1d0 + tz/(s0 + sc) );  w1 = ( 2d0/3d0 )*( 1d0 + tz/(s1 + sc) )
    r = ( w0*( -f + 3d0*u )*0.5d0 + w1*( u + d )*0.5d0 )/( w0 + w1 )
  End Function weno3z


  !> WENO5-Z (Borges et al. 2008), face value to the downwind side of cell v(0)
  Pure Function weno5z(vm2, vm1, v0, vp1, vp2) Result(r)
    !$acc routine seq
    Real(Int64), Intent(In) :: vm2, vm1, v0, vp1, vp2
    Real(Int64) :: r, p0, p1, p2, b0, b1, b2, t5, a0, a1, a2, eps

    p0 = ( 2d0*vm2 - 7d0*vm1 + 11d0*v0 )/6d0
    p1 = ( -vm1 + 5d0*v0 + 2d0*vp1 )/6d0
    p2 = ( 2d0*v0 + 5d0*vp1 - vp2 )/6d0
    b0 = 13d0/12d0*( vm2 - 2d0*vm1 + v0 )**2  + 0.25d0*( vm2 - 4d0*vm1 + 3d0*v0 )**2
    b1 = 13d0/12d0*( vm1 - 2d0*v0 + vp1 )**2  + 0.25d0*( vm1 - vp1 )**2
    b2 = 13d0/12d0*( v0 - 2d0*vp1 + vp2 )**2  + 0.25d0*( 3d0*v0 - 4d0*vp1 + vp2 )**2
    eps = 1d-30 + 1d-12*( vm2**2 + vm1**2 + v0**2 + vp1**2 + vp2**2 )
    t5 = Abs(b0 - b2)
    a0 = 0.1d0*( 1d0 + ( t5/(b0 + eps) )**2 )
    a1 = 0.6d0*( 1d0 + ( t5/(b1 + eps) )**2 )
    a2 = 0.3d0*( 1d0 + ( t5/(b2 + eps) )**2 )
    r = ( a0*p0 + a1*p1 + a2*p2 )/( a0 + a1 + a2 )
  End Function weno5z


  !> Momentum flux mfc*face value of the WENO5-Z scheme (mom_face, scheme 6, blended to upwind by the refill Courant number as in
  !  vof_twofluid fmom) for a row of n faces. The stencil rows sm2..s3 hold s(-2:3) of every face, mda/mdb the two cell mass fluxes
  !  whose mean is the face flux, r1/r2 the control-volume densities each side, ivcv = 1/(control-volume volume). Unit-stride rows
  !  without branches so that the loop vectorizes; the three weights share one division (rational form of the WENO-Z weights).
  Subroutine mom_flux_weno5z_row(n, sm2, sm1, s0, s1, s2, s3, mda, mdb, r1, r2, ivcv, cm0, cm1, f)

    Integer(Int32), Intent(In)  :: n
    Real(Int64),    Intent(In)  :: sm2(n), sm1(n), s0(n), s1(n), s2(n), s3(n), mda(n), mdb(n), r1(n), r2(n), ivcv, cm0, cm1
    Real(Int64),    Intent(Out) :: f(n)

    Integer(Int32) :: i
    Real(Int64) :: mf, vm2, vm1, v0, vp1, vp2, p0, p1, p2, b0, b1, b2, t5, eps, d0, d1, d2, n0, n1, n2, tt, r, up, t, icm
    Logical :: fwd

    icm = 1d0/(cm1 - cm0)
    !GCC$ ivdep
    Do i = 1, n
       mf = 0.5d0*( mda(i) + mdb(i) )
       fwd = ( mf >= 0d0 )
       vm2 = Merge(sm2(i), s3(i), fwd);  vm1 = Merge(sm1(i), s2(i), fwd);  v0 = Merge(s0(i), s1(i), fwd)
       vp1 = Merge(s1(i), s0(i), fwd);   vp2 = Merge(s2(i), sm1(i), fwd)
       p0 = ( 2d0*vm2 - 7d0*vm1 + 11d0*v0 )/6d0
       p1 = ( -vm1 + 5d0*v0 + 2d0*vp1 )/6d0
       p2 = ( 2d0*v0 + 5d0*vp1 - vp2 )/6d0
       b0 = 13d0/12d0*( vm2 - 2d0*vm1 + v0 )**2  + 0.25d0*( vm2 - 4d0*vm1 + 3d0*v0 )**2
       b1 = 13d0/12d0*( vm1 - 2d0*v0 + vp1 )**2  + 0.25d0*( vm1 - vp1 )**2
       b2 = 13d0/12d0*( v0 - 2d0*vp1 + vp2 )**2  + 0.25d0*( 3d0*v0 - 4d0*vp1 + vp2 )**2
       eps = 1d-30 + 1d-12*( vm2**2 + vm1**2 + v0**2 + vp1**2 + vp2**2 )
       t5 = Abs(b0 - b2)
       tt = t5*t5
       b0 = b0 + eps;  b1 = b1 + eps;  b2 = b2 + eps
       d0 = b0*b0;  d1 = b1*b1;  d2 = b2*b2
       n0 = 0.1d0*( d0 + tt )*d1*d2;  n1 = 0.6d0*( d1 + tt )*d0*d2;  n2 = 0.3d0*( d2 + tt )*d0*d1
       r = ( n0*p0 + n1*p1 + n2*p2 )/( n0 + n1 + n2 )
       up = v0
       t = Min(1d0, Max(0d0, ( Abs(mf)*Abs(r2(i) - r1(i))/( r1(i)*r2(i) )*ivcv - cm0 )*icm))
       f(i) = mf*( up + ( 1d0 - t*t*( 3d0 - 2d0*t ) )*( r - up ) )
    End Do

  End Subroutine mom_flux_weno5z_row

End Module mom_recon
