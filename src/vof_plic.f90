!> Geometry of a planar interface inside an axis-aligned box (PLIC): volume under a plane and the inverse (plane position for a given
!  volume fraction). Pure functions with no solver dependencies, so they can be unit-tested and reused for wave initialisation.
!
!  Conventions: cell coordinates are scaled to the unit cube centred on the origin (x_i in [-1/2,1/2]); the interface is m.x = alpha
!  with m the L1-normalised normal pointing from liquid to gas (m = -grad C), liquid on the side m.x < alpha.
Module vof_plic

  Use iso_fortran_env, Only : Int64

  Implicit None

  Real(Int64), Parameter :: vof_eps = 1d-12   ! C below this is gas, above 1-vof_eps is liquid

Contains

  !> Fraction of the unit cube [0,1]^3 satisfying a1*xi1 + a2*xi2 + a3*xi3 <= beta (any signs of a)
  Pure Function plic_box_fraction(a1, a2, a3, beta) Result(f)
    !$acc routine seq

    Real(Int64), Intent(In) :: a1, a2, a3, beta
    Real(Int64) :: f

    Real(Int64) :: s, b, n1, n2, n3, al, al0, b1, b2, b3, b12, bm, pr, tmp

    b = beta
    If ( a1 < 0d0 ) b = b - a1
    If ( a2 < 0d0 ) b = b - a2
    If ( a3 < 0d0 ) b = b - a3
    s = Abs(a1) + Abs(a2) + Abs(a3)

    If ( s < 1d-300 ) Then
       f = Merge(1d0, 0d0, beta > 0d0)
       Return
    End If
    If ( b <= 0d0 ) Then
       f = 0d0
       Return
    End If
    If ( b >= s ) Then
       f = 1d0
       Return
    End If

    n1 = Abs(a1)/s
    n2 = Abs(a2)/s
    n3 = Abs(a3)/s
    al  = b/s
    al0 = Min(al, 1d0 - al)

    b1 = Min(n1, n2)
    b3 = Max(n1, n2)
    b2 = n3
    If ( b2 < b1 ) Then
       tmp = b1;  b1 = b2;  b2 = tmp
    Else If ( b2 > b3 ) Then
       tmp = b3;  b3 = b2;  b2 = tmp
    End If
    b12 = b1 + b2
    bm  = Min(b12, b3)
    pr  = Max(6d0*b1*b2*b3, 1d-50)

    If ( al0 < b1 ) Then
       tmp = al0*al0*al0/pr
    Else If ( al0 < b2 ) Then
       tmp = 0.5d0*al0*(al0 - b1)/(b2*b3) + b1*b1*b1/pr
    Else If ( al0 < bm ) Then
       tmp = ( al0*al0*(3d0*b12 - al0) + b1*b1*(b1 - 3d0*al0) + b2*b2*(b2 - 3d0*al0) )/pr
    Else If ( b12 < b3 ) Then
       tmp = (al0 - 0.5d0*bm)/b3
    Else
       tmp = ( al0*al0*(3d0 - 2d0*al0) + b1*b1*(b1 - 3d0*al0) + b2*b2*(b2 - 3d0*al0) + b3*b3*(b3 - 3d0*al0) )/pr
    End If

    If ( al <= 0.5d0 ) Then
       f = tmp
    Else
       f = 1d0 - tmp
    End If

  End Function plic_box_fraction


  !> Plane constant alpha (centred cell coordinates, normal L1-normalised on return) for liquid fraction c in 0<c<1
  Pure Function plic_alpha(c, m1, m2, m3) Result(alpha)
    !$acc routine seq

    Real(Int64), Intent(In) :: c, m1, m2, m3
    Real(Int64) :: alpha

    Real(Int64) :: s, n1, n2, n3, ma, mb, mc, tmp, m12, pr, v1, v2, v3, mm, ch, p12, q, teta, cs, p

    s = Abs(m1) + Abs(m2) + Abs(m3)
    n1 = Abs(m1)/s
    n2 = Abs(m2)/s
    n3 = Abs(m3)/s

    ma = Min(n1, n2)
    mc = Max(n1, n2)
    mb = n3
    If ( mb < ma ) Then
       tmp = ma;  ma = mb;  mb = tmp
    Else If ( mb > mc ) Then
       tmp = mc;  mc = mb;  mb = tmp
    End If
    m12 = ma + mb
    pr  = Max(6d0*ma*mb*mc, 1d-50)
    v1  = ma*ma*ma/pr
    v2  = v1 + (mb - ma)/(2d0*mc)
    If ( mc < m12 ) Then
       mm = mc
       v3 = ( mc*mc*(3d0*m12 - mc) + ma*ma*(ma - 3d0*mc) + mb*mb*(mb - 3d0*mc) )/pr
    Else
       mm = m12
       v3 = mm/(2d0*mc)
    End If

    ch = Min(c, 1d0 - c)
    If ( ch < v1 ) Then
       alpha = (pr*ch)**(1d0/3d0)
    Else If ( ch < v2 ) Then
       alpha = ( ma + Sqrt(ma*ma + 8d0*mb*mc*(ch - v1)) )/2d0
    Else If ( ch < v3 ) Then
       p12  = Sqrt(2d0*ma*mb)
       q    = 3d0*(m12 - 2d0*mc*ch)/(4d0*p12)
       teta = ACos(Max(-1d0, Min(1d0, q)))/3d0
       cs   = Cos(teta)
       alpha = p12*( Sqrt(3d0*(1d0 - cs*cs)) - cs ) + m12
    Else If ( m12 < mc ) Then
       alpha = mc*ch + mm/2d0
    Else
       p    = ma*(mb + mc) + mb*mc - 0.25d0
       p12  = Sqrt(p)
       q    = 3d0*ma*mb*mc*(0.5d0 - ch)/(2d0*p*p12)
       teta = ACos(Max(-1d0, Min(1d0, q)))/3d0
       cs   = Cos(teta)
       alpha = p12*( Sqrt(3d0*(1d0 - cs*cs)) - cs ) + 0.5d0
    End If
    If ( c > 0.5d0 ) alpha = 1d0 - alpha

    alpha = alpha - 0.5d0

  End Function plic_alpha


  !> Liquid fraction of a sub-box of the cell for the plane m.x = alpha (m L1-normalised): the box spans x_i in [xlo_i, xlo_i + w_i]
  Pure Function plic_subbox_fraction(m1, m2, m3, alpha, xlo1, xlo2, xlo3, w1, w2, w3) Result(f)
    !$acc routine seq

    Real(Int64), Intent(In) :: m1, m2, m3, alpha, xlo1, xlo2, xlo3, w1, w2, w3
    Real(Int64) :: f

    f = plic_box_fraction( m1*w1, m2*w2, m3*w3, alpha - m1*xlo1 - m2*xlo2 - m3*xlo3 )

  End Function plic_subbox_fraction

End Module vof_plic
