!> Spike A: how well does the constant-coefficient fast Poisson solver (cosine transforms in x and y, the structure of the solver's
!  FFT/tridiagonal chain) precondition the variable-density pressure equation div( (1/rho) grad p ) = f across a sharp interface?
!  Serial 2-D, Neumann walls; reports preconditioned-CG iteration counts versus density ratio. The constant coefficient of the
!  preconditioner is irrelevant to CG (a scalar multiple of M gives identical iterates), so only the ratio matters.
!  Result (white-noise rhs, cold start): iterations to 1e-6 are 17 / 28 / 37-39 at ratio 10 / 100 / 1000, independent of N (64, 128).
!  The fast solver therefore works as a preconditioner but costs ~30-60 solves per projection at ratio 1e3 -- 3x that per RK step.
Program spike_pressure_pcg

  Use iso_fortran_env, Only : Int32, Int64
  Implicit None

  Real(Int64), Parameter :: pi = 3.14159265358979323846d0
  Integer(Int32) :: n, ir, ip, itols(5), nratio
  Real(Int64) :: ratio, rho_l, rho_g, beta0
  Real(Int64), Allocatable :: rho(:,:), qx(:,:), lam(:)
  Real(Int64) :: ratios(4) = (/ 1d0, 10d0, 100d0, 1000d0 /)
  Character(len=24) :: names(4) = (/ Character(len=24) :: 'beta0 = 1/rho_l (min)', 'beta0 = 1/rho_g (max)', &
                                      'beta0 = vol-mean beta', 'beta0 = 1/mean(rho)' /)

  nratio = 4
  Do n = 64, 128, 64
     Write(*,'(a,i0)') 'N = ', n
     Allocate( rho(n,n), qx(n,n), lam(n) )
     Call build_dct(n, qx, lam)
     Do ir = 1, nratio
        ratio = ratios(ir)
        rho_l = 1d0;  rho_g = rho_l/ratio
        Call build_density(n, rho_l, rho_g, rho)
        Write(*,'(a,es8.1)') '  rho_l/rho_g = ', ratio
        Do ip = 3, 3
           Select Case(ip)
           Case(1); beta0 = 1d0/rho_l
           Case(2); beta0 = 1d0/rho_g
           Case(3); beta0 = Sum(1d0/rho)/Real(n*n,Int64)
           Case(4); beta0 = 1d0/(Sum(rho)/Real(n*n,Int64))
           End Select
           Call pcg(n, rho, qx, lam, beta0, itols)
           Write(*,'(4x,a,5i6)') 'iterations to 1e-2,1e-4,1e-6,1e-8,1e-10 :', itols
        End Do
     End Do
     Deallocate( rho, qx, lam )
  End Do

Contains

  Subroutine build_density(n, rl, rg, rho)

    Integer(Int32), Intent(In)  :: n
    Real(Int64),    Intent(In)  :: rl, rg
    Real(Int64),    Intent(Out) :: rho(n,n)
    Integer(Int32) :: i, j
    Real(Int64) :: x, y, h

    h = 1d0/n
    Do j = 1, n
       y = (j-0.5d0)*h
       Do i = 1, n
          x = (i-0.5d0)*h
          rho(i,j) = Merge(rl, rg, y < 0.5d0 + 0.05d0*Cos(2d0*pi*x))
       End Do
    End Do

  End Subroutine build_density


  !> Orthonormal DCT-II basis of the 1-D Neumann second difference on n cells (h=1/n) and its eigenvalues (<=0)
  Subroutine build_dct(n, q, lam)

    Integer(Int32), Intent(In)  :: n
    Real(Int64),    Intent(Out) :: q(n,n), lam(n)
    Integer(Int32) :: i, k
    Real(Int64) :: h, nrm

    h = 1d0/n
    Do k = 0, n-1
       nrm = Sqrt(2d0/n);  If ( k == 0 ) nrm = Sqrt(1d0/n)
       Do i = 1, n
          q(i,k+1) = nrm*Cos(pi*k*(i-0.5d0)/n)
       End Do
       lam(k+1) = -(2d0 - 2d0*Cos(pi*k/n))/h**2
    End Do

  End Subroutine build_dct


  Subroutine apply_a(n, rho, p, ap)

    Integer(Int32), Intent(In)  :: n
    Real(Int64),    Intent(In)  :: rho(n,n), p(n,n)
    Real(Int64),    Intent(Out) :: ap(n,n)
    Integer(Int32) :: i, j
    Real(Int64) :: h2, fe, fw, fn, fs

    h2 = (1d0/n)**2
    Do j = 1, n
       Do i = 1, n
          fe = 0d0;  fw = 0d0;  fn = 0d0;  fs = 0d0
          If ( i < n ) fe = (p(i+1,j) - p(i,j)) * 2d0/(rho(i,j) + rho(i+1,j))
          If ( i > 1 ) fw = (p(i,j) - p(i-1,j)) * 2d0/(rho(i,j) + rho(i-1,j))
          If ( j < n ) fn = (p(i,j+1) - p(i,j)) * 2d0/(rho(i,j) + rho(i,j+1))
          If ( j > 1 ) fs = (p(i,j) - p(i,j-1)) * 2d0/(rho(i,j) + rho(i,j-1))
          ap(i,j) = (fe - fw + fn - fs)/h2
       End Do
    End Do

  End Subroutine apply_a


  Subroutine precond(n, qx, lam, beta0, r, z)

    Integer(Int32), Intent(In)  :: n
    Real(Int64),    Intent(In)  :: qx(n,n), lam(n), beta0, r(n,n)
    Real(Int64),    Intent(Out) :: z(n,n)
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


  Subroutine pcg(n, rho, qx, lam, beta0, itols)

    Integer(Int32), Intent(In)  :: n
    Real(Int64),    Intent(In)  :: rho(n,n), qx(n,n), lam(n), beta0
    Integer(Int32), Intent(Out) :: itols(5)
    Real(Int64), Allocatable :: b(:,:), x(:,:), r(:,:), z(:,:), p(:,:), ap(:,:)
    Real(Int64) :: rz, rzn, alpha, beta, r0, tol
    Integer(Int32) :: it, i, j
    Real(Int64) :: seed

    Allocate( b(n,n), x(n,n), r(n,n), z(n,n), p(n,n), ap(n,n) )
    seed = 0.37d0
    Do j = 1, n
       Do i = 1, n
          seed = Mod(seed*9301d0 + 0.49297d0, 1d0)
          b(i,j) = seed - 0.5d0
       End Do
    End Do
    b = b - Sum(b)/Real(n*n,Int64)
    x = 0d0
    r = b
    r0 = Sqrt(Sum(r*r))
    tol = 1d-8
    Call precond(n, qx, lam, beta0, r, z)
    p = z
    rz = Sum(r*z)
    itols = 2000
    Do it = 1, 2000
       Call apply_a(n, rho, p, ap)
       alpha = rz/Sum(p*ap)
       x = x + alpha*p
       r = r - alpha*ap
       r = r - Sum(r)/Real(n*n,Int64)
       Do j = 1, 5
          If ( itols(j) == 2000 .And. Sqrt(Sum(r*r)) < 10d0**(-2*j)*r0 ) itols(j) = it
       End Do
       If ( itols(5) < 2000 ) Exit
       Call precond(n, qx, lam, beta0, r, z)
       rzn = Sum(r*z)
       beta = rzn/rz
       rz = rzn
       p = z + beta*p
    End Do

  End Subroutine pcg

End Program spike_pressure_pcg
