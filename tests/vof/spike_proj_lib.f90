!> Shared helpers for the pressure spikes: serial 2-D, Neumann walls, cell-centred p and face velocities on an n x n grid of the unit
!  square. Constant-coefficient fast solve (cosine-transform basis, the structure of the production FFT chain) and an exact
!  variable-coefficient solve D(beta G p) = f by CG preconditioned with the fast solve.
Module spike_proj_lib

  Use iso_fortran_env, Only : Int32, Int64
  Implicit None

  Real(Int64), Parameter :: pi = 3.14159265358979323846d0
  Integer(Int32) :: n = 0
  Real(Int64) :: h = 0d0
  Real(Int64), Allocatable :: qx(:,:), lam(:)

Contains

  Subroutine proj_init(nn)

    Integer(Int32), Intent(In) :: nn
    Integer(Int32) :: i, k
    Real(Int64) :: nrm

    n = nn;  h = 1d0/n
    If ( Allocated(qx) ) Deallocate( qx, lam )
    Allocate( qx(n,n), lam(n) )
    Do k = 0, n-1
       nrm = Sqrt(2d0/n);  If ( k == 0 ) nrm = Sqrt(1d0/n)
       Do i = 1, n
          qx(i,k+1) = nrm*Cos(pi*k*(i-0.5d0)/n)
       End Do
       lam(k+1) = -(2d0 - 2d0*Cos(pi*k/n))/h**2
    End Do

  End Subroutine proj_init


  !> Face gradient times face coefficient beta = 2/(rho_i+rho_j) (interior faces); constant_beta overrides the coefficient
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


  !> Constant-coefficient solve beta0 * L p = f (zero-mean f and p)
  Subroutine fast_solve(beta0, f, p)

    Real(Int64), Intent(In)  :: beta0, f(n,n)
    Real(Int64), Intent(Out) :: p(n,n)
    Real(Int64) :: t(n,n)
    Integer(Int32) :: a, b

    t = MatMul(Transpose(qx), MatMul(f - Sum(f)/(n*n), qx))
    Do b = 1, n
       Do a = 1, n
          If ( a == 1 .And. b == 1 ) Then
             t(a,b) = 0d0
          Else
             t(a,b) = t(a,b)/(beta0*(lam(a) + lam(b)))
          End If
       End Do
    End Do
    p = MatMul(qx, MatMul(t, Transpose(qx)))

  End Subroutine fast_solve


  !> Exact solve of D(beta G p) = f by PCG (fast solve with the mean beta as preconditioner) to a relative residual of 1e-13
  Subroutine exact_solve(rho, f, p, iters)

    Real(Int64), Intent(In)  :: rho(n,n), f(n,n)
    Real(Int64), Intent(Out) :: p(n,n)
    Integer(Int32), Intent(Out), Optional :: iters
    Real(Int64) :: r(n,n), z(n,n), d(n,n), ap(n,n), fx(n-1,n), fy(n,n-1), rz, rzn, alpha, r0, b0
    Integer(Int32) :: it

    b0 = Sum(1d0/rho)/(n*n)
    p = 0d0
    r = f - Sum(f)/(n*n)
    r0 = Sqrt(Sum(r*r))
    If ( r0 == 0d0 ) Then
       If ( Present(iters) ) iters = 0
       Return
    End If
    Call fast_solve(b0, r, z)
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
       Call fast_solve(b0, r, z)
       rzn = Sum(r*z)
       d = z + (rzn/rz)*d
       rz = rzn
    End Do
    If ( Present(iters) ) iters = it

  End Subroutine exact_solve


  !> nit PCG iterations (fast-solve preconditioner with beta0) on D(beta G x) = f from x = 0
  Subroutine cg_iterations(rho, f, beta0, nit, x)

    Real(Int64), Intent(In)  :: rho(n,n), f(n,n), beta0
    Integer(Int32), Intent(In) :: nit
    Real(Int64), Intent(Out) :: x(n,n)
    Real(Int64) :: r(n,n), z(n,n), d(n,n), ap(n,n), fx(n-1,n), fy(n,n-1), rz, rzn, alpha
    Integer(Int32) :: it

    x = 0d0
    r = f - Sum(f)/(n*n)
    If ( Maxval(Abs(r)) == 0d0 ) Return
    Call fast_solve(beta0, r, z)
    d = z
    rz = Sum(r*z)
    Do it = 1, nit
       Call grad_beta(rho, d, fx, fy)
       Call divergence(fx, fy, ap)
       alpha = rz/Sum(d*ap)
       x = x + alpha*d
       r = r - alpha*ap
       r = r - Sum(r)/(n*n)
       If ( it == nit ) Exit
       Call fast_solve(beta0, r, z)
       rzn = Sum(r*z)
       d = z + (rzn/rz)*d
       rz = rzn
    End Do

  End Subroutine cg_iterations

End Module spike_proj_lib
