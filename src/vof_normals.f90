!> Interface-normal estimators for PLIC, independent of the flux and update code. Each takes the 3x3x3 neighbourhood of C (in index
!  space, so cell aspect ratio is carried by the caller's scaling of the plane) and returns the L1-normalised normal m = -grad C
!  pointing from liquid to gas.
Module vof_normals

  Use iso_fortran_env, Only : Int32, Int64
  Use vof_plic, Only : vof_eps

  Implicit None

  Integer(Int32), Parameter :: VOF_SCHEME_YOUNGS = 1, VOF_SCHEME_CC = 2

Contains

  !> Youngs finite-difference gradient (27-point, weights 4/2/1)
  Pure Subroutine vof_normal_youngs(c, m)
    !$acc routine seq

    Real(Int64), Intent(In)  :: c(-1:1,-1:1,-1:1)
    Real(Int64), Intent(Out) :: m(3)

    Integer(Int32) :: a, b
    Real(Int64) :: g(3), wt, s

    g = 0d0
    Do b = -1, 1
       Do a = -1, 1
          wt = 1d0
          If ( a == 0 .And. b == 0 ) Then
             wt = 4d0
          Else If ( a == 0 .Or. b == 0 ) Then
             wt = 2d0
          End If
          g(1) = g(1) + wt*( c( 1,a,b) - c(-1,a,b) )
          g(2) = g(2) + wt*( c(a, 1,b) - c(a,-1,b) )
          g(3) = g(3) + wt*( c(a,b, 1) - c(a,b,-1) )
       End Do
    End Do

    s = Abs(g(1)) + Abs(g(2)) + Abs(g(3))
    If ( s < 1d-300 ) Then
       m = (/ 0d0, 0d0, 1d0 /)
    Else
       m = -g/s
    End If

  End Subroutine vof_normal_youngs


  !> Centred-columns height-function normal along the Youngs-dominant axis, falling back to Youngs when the 5 columns used are not
  !  all saturated at both ends (the interface does not cross the whole 3-cell column)
  Pure Subroutine vof_normal_cc(c, m)
    !$acc routine seq

    Real(Int64), Intent(In)  :: c(-1:1,-1:1,-1:1)
    Real(Int64), Intent(Out) :: m(3)

    Real(Int64) :: my(3), h(-1:1,-1:1), sgn, s, lo, hi
    Integer(Int32) :: d, a, b, p, q
    Logical :: complete

    Call vof_normal_youngs(c, my)
    d = 1
    If ( Abs(my(2)) > Abs(my(d)) ) d = 2
    If ( Abs(my(3)) > Abs(my(d)) ) d = 3
    sgn = Sign(1d0, my(d))

    ! columns along d at transverse offsets (a,b) in the two other axes
    Do b = -1, 1
       Do a = -1, 1
          h(a,b) = 0d0
          Do p = -1, 1
             If ( d == 1 ) Then
                h(a,b) = h(a,b) + c(p,a,b)
             Else If ( d == 2 ) Then
                h(a,b) = h(a,b) + c(a,p,b)
             Else
                h(a,b) = h(a,b) + c(a,b,p)
             End If
          End Do
       End Do
    End Do

    complete = .True.
    Do q = 1, 5
       If ( q == 1 ) Then
          a = 0;  b = 0
       Else If ( q == 2 ) Then
          a = 1;  b = 0
       Else If ( q == 3 ) Then
          a = -1; b = 0
       Else If ( q == 4 ) Then
          a = 0;  b = 1
       Else
          a = 0;  b = -1
       End If
       If ( d == 1 ) Then
          lo = c(-1,a,b);  hi = c(1,a,b)
       Else If ( d == 2 ) Then
          lo = c(a,-1,b);  hi = c(a,1,b)
       Else
          lo = c(a,b,-1);  hi = c(a,b,1)
       End If
       ! gas on the +d side of the interface (sgn>0): liquid at the low end, gas at the high end, and vice versa
       If ( sgn > 0d0 ) Then
          If ( lo < 1d0 - vof_eps .Or. hi > vof_eps ) complete = .False.
       Else
          If ( hi < 1d0 - vof_eps .Or. lo > vof_eps ) complete = .False.
       End If
    End Do

    If ( .Not. complete ) Then
       m = my
       Return
    End If

    ! m_a = -dH/da for both orientations (liquid below: h=H, m_d=+1; liquid above: h=top-H, m_d=-1); the two transverse
    ! axes of h(a,b) are in increasing axis order
    m(d) = sgn
    If ( d == 1 ) Then
       m(2) = -0.5d0*( h(1,0) - h(-1,0) )
       m(3) = -0.5d0*( h(0,1) - h(0,-1) )
    Else If ( d == 2 ) Then
       m(1) = -0.5d0*( h(1,0) - h(-1,0) )
       m(3) = -0.5d0*( h(0,1) - h(0,-1) )
    Else
       m(1) = -0.5d0*( h(1,0) - h(-1,0) )
       m(2) = -0.5d0*( h(0,1) - h(0,-1) )
    End If
    s = Abs(m(1)) + Abs(m(2)) + Abs(m(3))
    m = m/s

  End Subroutine vof_normal_cc

End Module vof_normals
