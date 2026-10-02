!> Spike A': single-solve pressure splitting versus the exact variable-density projection (serial 2-D, Neumann walls, cell-centred p,
!  face velocities; the constant-coefficient solve is a cosine-transform solver, the same structure as the production FFT chain).
!
!  Exact problem: D(beta G p) = D u*/dt with beta = 1/rho (arithmetic-rho face average), u = u* - dt beta G p. The exact p is found
!  by preconditioned CG to round-off and then perturbed to ptilde = p + delta to mimic an extrapolated-pressure error. Variants:
!    V1 (Dodd-Ferrante): beta0 L p1 = D u*/dt + D((beta0-beta) G ptilde), u = u* - dt beta G p1       -> residual divergence
!    V2 (incremental)  : u** = u* - dt beta G ptilde, beta0 L dp = D u**/dt, u = u** - dt beta0 G dp  -> exactly div-free, force error
!  Errors are reported per phase (liquid / gas / mixed faces) relative to the largest exact pressure force in that class, and the
!  divergence relative to the largest exact pressure force: both are the fraction of the pressure force that is wrong.
!
!  Result (N=64, ratio 1..1000): V2 reduces to "explicit pressure predictor + ordinary constant-coefficient projection" (its result
!  is independent of beta0). It is divergence-free to round-off at every ratio, needs no variable-coefficient solve, and its error
!  is the solenoidal part of beta*grad(p - ptilde) (~ grad(beta) x grad(delta), confined to the interface): ~0.5*eps of the local
!  pressure force when delta = eps*p, independent of the density ratio up to 1000. V1 leaves a residual divergence ~eps and is
!  only usable with beta0 = max(beta) (1/rho_g); with beta0 = 1/rho_l its gas-phase error is ~500*eps. When delta is spatially
!  unrelated to p (shape 2) the gas-phase relative error is large for both: the extrapolation error must share p's structure.
Program spike_pressure_split

  Use iso_fortran_env, Only : Int32, Int64
  Implicit None

  Real(Int64), Parameter :: pi = 3.14159265358979323846d0
  Integer(Int32), Parameter :: n = 64
  Real(Int64) :: h, ratio, rho_l, rho_g, beta0, eps, beta0s(3), epss(3), ratios(4)
  Real(Int64), Allocatable :: rho(:,:), qx(:,:), lam(:), pex(:,:), uxs(:,:), uys(:,:), uxe(:,:), uye(:,:)
  Integer(Int32) :: ir, ib, ie, ishape
  Character(len=14) :: bnames(3) = (/ Character(len=14) :: 'beta0=1/rho_l  ', 'beta0=1/rho_g  ', 'beta0=mean     ' /)

  h = 1d0/n
  ratios = (/ 1d0, 10d0, 100d0, 1000d0 /)
  epss = (/ 1d-1, 1d-2, 1d-3 /)
  Allocate( rho(n,n), qx(n,n), lam(n), pex(n,n), uxs(n-1,n), uys(n,n-1), uxe(n-1,n), uye(n,n-1) )
  Call build_dct(qx, lam)

  Write(*,'(a)') 'columns: V1 [liq gas mix | div]   V2 [liq gas mix | div]   (relative to the exact pressure force in that class)'
  Do ir = 1, 4
     ratio = ratios(ir)
     rho_l = 1d0;  rho_g = rho_l/ratio
     Call build_density(rho_l, rho_g, rho)
     Call make_ustar(uxs, uys)
     Call exact_solution(rho, qx, lam, uxs, uys, pex, uxe, uye)
     beta0s = (/ 1d0/rho_l, 1d0/rho_g, Sum(1d0/rho)/(n*n) /)
     Write(*,'(a,es8.1)') 'rho_l/rho_g = ', ratio
     Do ib = 1, 3
        beta0 = beta0s(ib)
        Do ishape = 1, 2
           Do ie = 1, 3
              eps = epss(ie)
              Call run_variants(rho, qx, lam, beta0, eps, ishape, bnames(ib))
           End Do
        End Do
     End Do
  End Do

Contains

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


  Subroutine make_ustar(ux, uy)

    Real(Int64), Intent(Out) :: ux(n-1,n), uy(n,n-1)
    Integer(Int32) :: i, j

    Do j = 1, n
       Do i = 1, n-1
          ux(i,j) = Sin(2d0*pi*i*h)*Cos(2d0*pi*(j-0.5d0)*h)
       End Do
    End Do
    Do j = 1, n-1
       Do i = 1, n
          uy(i,j) = Sin(pi*j*h)*Cos(2d0*pi*(i-0.5d0)*h)
       End Do
    End Do

  End Subroutine make_ustar


  Subroutine build_dct(q, lam)

    Real(Int64), Intent(Out) :: q(n,n), lam(n)
    Integer(Int32) :: i, k
    Real(Int64) :: nrm

    Do k = 0, n-1
       nrm = Sqrt(2d0/n);  If ( k == 0 ) nrm = Sqrt(1d0/n)
       Do i = 1, n
          q(i,k+1) = nrm*Cos(pi*k*(i-0.5d0)/n)
       End Do
       lam(k+1) = -(2d0 - 2d0*Cos(pi*k/n))/h**2
    End Do

  End Subroutine build_dct


  !> Face gradient times face coefficient (beta=1/rho_f, rho_f arithmetic) : fx = beta*(p(i+1)-p(i))/h, interior faces only
  Subroutine grad_beta(rho, p, fx, fy, constant_beta)

    Real(Int64), Intent(In)  :: rho(n,n), p(n,n)
    Real(Int64), Intent(Out) :: fx(n-1,n), fy(n,n-1)
    Real(Int64), Intent(In), Optional :: constant_beta
    Integer(Int32) :: i, j
    Real(Int64) :: bx, by

    Do j = 1, n
       Do i = 1, n-1
          bx = 2d0/(rho(i,j) + rho(i+1,j))
          If ( Present(constant_beta) ) bx = constant_beta
          fx(i,j) = bx*(p(i+1,j) - p(i,j))/h
       End Do
    End Do
    Do j = 1, n-1
       Do i = 1, n
          by = 2d0/(rho(i,j) + rho(i,j+1))
          If ( Present(constant_beta) ) by = constant_beta
          fy(i,j) = by*(p(i,j+1) - p(i,j))/h
       End Do
    End Do

  End Subroutine grad_beta


  Subroutine divergence(fx, fy, d)

    Real(Int64), Intent(In)  :: fx(n-1,n), fy(n,n-1)
    Real(Int64), Intent(Out) :: d(n,n)
    Integer(Int32) :: i, j
    Real(Int64) :: e, w, no, so

    Do j = 1, n
       Do i = 1, n
          e = 0d0;  w = 0d0;  no = 0d0;  so = 0d0
          If ( i < n ) e = fx(i,j)
          If ( i > 1 ) w = fx(i-1,j)
          If ( j < n ) no = fy(i,j)
          If ( j > 1 ) so = fy(i,j-1)
          d(i,j) = (e - w + no - so)/h
       End Do
    End Do

  End Subroutine divergence


  Subroutine precond(beta0, r, z)

    Real(Int64), Intent(In)  :: beta0, r(n,n)
    Real(Int64), Intent(Out) :: z(n,n)
    Real(Int64) :: t(n,n)
    Integer(Int32) :: a, b

    t = MatMul(Transpose(qx), MatMul(r, qx))
    Do b = 1, n
       Do a = 1, n
          If ( a == 1 .And. b == 1 ) Then
             t(a,b) = 0d0
          Else
             t(a,b) = t(a,b)/(beta0*(lam(a) + lam(b)))
          End If
       End Do
    End Do
    z = MatMul(qx, MatMul(t, Transpose(qx)))

  End Subroutine precond


  !> Constant-coefficient fast solve L p = f/beta0 (zero-mean p)
  Subroutine fast_solve(beta0, f, p)

    Real(Int64), Intent(In)  :: beta0, f(n,n)
    Real(Int64), Intent(Out) :: p(n,n)

    Call precond(beta0, f - Sum(f)/(n*n), p)

  End Subroutine fast_solve


  !> Exact variable-density projection: D(beta G p) = D u* (dt = 1) by CG to round-off; returns p and the projected face velocity
  Subroutine exact_solution(rho, qx_, lam_, uxs, uys, p, ux, uy)

    Real(Int64), Intent(In)  :: rho(n,n), qx_(n,n), lam_(n), uxs(n-1,n), uys(n,n-1)
    Real(Int64), Intent(Out) :: p(n,n), ux(n-1,n), uy(n,n-1)
    Real(Int64) :: f(n,n), r(n,n), z(n,n), d(n,n), ap(n,n), fx(n-1,n), fy(n,n-1), rz, rzn, alpha, r0
    Integer(Int32) :: it
    Real(Int64) :: b0

    Call divergence(uxs, uys, f)
    f = f - Sum(f)/(n*n)
    b0 = Sum(1d0/rho)/(n*n)
    p = 0d0
    r = f
    r0 = Sqrt(Sum(r*r))
    Call precond(b0, r, z)
    d = z
    rz = Sum(r*z)
    Do it = 1, 500
       Call grad_beta(rho, d, fx, fy)
       Call divergence(fx, fy, ap)
       alpha = rz/Sum(d*ap)
       p = p + alpha*d
       r = r - alpha*ap
       r = r - Sum(r)/(n*n)
       If ( Sqrt(Sum(r*r)) < 1d-13*r0 ) Exit
       Call precond(b0, r, z)
       rzn = Sum(r*z)
       d = z + (rzn/rz)*d
       rz = rzn
    End Do
    Call grad_beta(rho, p, fx, fy)
    ux = uxs - fx
    uy = uys - fy

  End Subroutine exact_solution


  Subroutine run_variants(rho, qx_, lam_, beta0, eps, ishape, label)

    Real(Int64), Intent(In) :: rho(n,n), qx_(n,n), lam_(n), beta0, eps
    Integer(Int32), Intent(In) :: ishape
    Character(len=*), Intent(In) :: label
    Real(Int64) :: pt(n,n), delta(n,n), f(n,n), p1(n,n), dp(n,n), fx(n-1,n), fy(n,n-1), cx(n-1,n), cy(n,n-1), d(n,n)
    Real(Int64) :: u1x(n-1,n), u1y(n,n-1), u2x(n-1,n), u2y(n,n-1), f0x(n-1,n), f0y(n,n-1), pmax
    Real(Int64) :: e1(3), e2(3), f0c(3), div1, div2, f0max
    Integer(Int32) :: i, j, cls
    Real(Int64) :: x, y

    pmax = Maxval(Abs(pex))
    Do j = 1, n
       Do i = 1, n
          x = (i-0.5d0)*h;  y = (j-0.5d0)*h
          If ( ishape == 1 ) Then
             delta(i,j) = eps*pex(i,j)
          Else
             delta(i,j) = eps*pmax*Cos(3d0*pi*x)*Cos(2d0*pi*y)
          End If
       End Do
    End Do
    pt = pex + delta

    ! V1
    Call grad_beta(rho, pt, fx, fy)
    Call grad_beta(rho, pt, cx, cy, beta0)
    Call divergence(cx - fx, cy - fy, d)
    Call divergence(uxs, uys, f)
    f = f + d
    Call fast_solve(beta0, f, p1)
    Call grad_beta(rho, p1, fx, fy)
    u1x = uxs - fx
    u1y = uys - fy
    Call divergence(u1x, u1y, d)
    div1 = Maxval(Abs(d))*h

    ! V2
    Call grad_beta(rho, pt, fx, fy)
    Call divergence(uxs - fx, uys - fy, f)
    Call fast_solve(beta0, f, dp)
    Call grad_beta(rho, dp, fx, fy, beta0)
    Call grad_beta(rho, pt, cx, cy)
    u2x = uxs - cx - fx
    u2y = uys - cy - fy
    Call divergence(u2x, u2y, d)
    div2 = Maxval(Abs(d))*h

    f0x = uxs - uxe
    f0y = uys - uye
    f0max = Max(Maxval(Abs(f0x)), Maxval(Abs(f0y)))
    e1 = 0d0;  e2 = 0d0;  f0c = 0d0
    Do j = 1, n
       Do i = 1, n-1
          cls = face_class(rho(i,j), rho(i+1,j))
          f0c(cls) = Max(f0c(cls), Abs(f0x(i,j)))
          e1(cls) = Max(e1(cls), Abs(u1x(i,j) - uxe(i,j)))
          e2(cls) = Max(e2(cls), Abs(u2x(i,j) - uxe(i,j)))
       End Do
    End Do
    Do j = 1, n-1
       Do i = 1, n
          cls = face_class(rho(i,j), rho(i,j+1))
          f0c(cls) = Max(f0c(cls), Abs(f0y(i,j)))
          e1(cls) = Max(e1(cls), Abs(u1y(i,j) - uye(i,j)))
          e2(cls) = Max(e2(cls), Abs(u2y(i,j) - uye(i,j)))
       End Do
    End Do
    Write(*,'(2x,a,a,i1,a,es8.1,a,3es9.2,a,es9.2,a,3es9.2,a,es9.2)') label, ' shape', ishape, ' eps', eps, '  V1', &
         e1/Max(f0c, 1d-300), ' |', div1/f0max, '  V2', e2/Max(f0c, 1d-300), ' |', div2/f0max

  End Subroutine run_variants


  Pure Function face_class(ra, rb) Result(c)

    Real(Int64), Intent(In) :: ra, rb
    Integer(Int32) :: c

    If ( ra == rb ) Then
       c = Merge(1, 2, ra == Maxval((/ ra, rb /)) .And. ra >= 0.5d0)
    Else
       c = 3
    End If

  End Function face_class

End Program spike_pressure_split
