!> Spike B: temporal treatment of the VOF advance inside the solver's RK3 (Wray) step, for a prescribed, discretely divergence-free,
!  time-dependent progressive-wave velocity field (streamfunction psi = A cos(kx - w t) sin(pi y), closed top/bottom, periodic in x)
!  with a flat interface that is displaced by a few cells -- the regime wave simulations sit in (Co <= 0.5 gives w*dt ~ 0.05-0.3).
!  Each strategy is compared with a converged reference of the same spatial scheme (dt/8, exact midpoint velocity), which isolates the
!  temporal error. Strategies (stage velocities are taken exact, i.e. the best case for the stage-based ones):
!    S1 one step, exact midpoint velocity (unreachable upper bound)     S4 one step, AB2-extrapolated midpoint 1.5u^n - 0.5u^(n-1)
!    S2 one step, start-of-step velocity u^n                            S5 three sub-steps (8/15, 2/15, 1/3 dt) with u^n, u^(1), u^(2)
!    S3 one step, stage-1 velocity u^(1) (t + 8/15 dt)                  S6 three sub-steps with the exact velocity at each sub-interval midpoint
!    S7 three sub-steps with linearly extrapolated midpoint velocities built only from stage velocities already available when the
!       sub-step is needed (u^n and the previous step's u^(2); u^n, u^(1); u^(1), u^(2)) -- the variant a stage-consistent solver can run
!
!  Result (N=64, 5 periods, displacement amplitude 3 cells, w*dt = 0.3/0.15/0.075; L1 vs converged reference at w*dt = 0.15):
!    S2 u^n 3.5e-5 (1st order) | S5 stage-left-point sub-steps 1.4e-5 (1st order, 3x cost) | S4 AB2 9e-6 (2nd) | S3 u^(1) 2.9e-6
!    S7 stage-consistent sub-steps with extrapolated midpoint velocities 1.5e-6 (2nd order, observed 1.97) | S1 exact midpoint 1.2e-6.
!  S7 also yields C at every stage time t_s, which a mass-consistent momentum update needs for rho^(s); S3 only gives C^(n+1).
Program spike_c_timeint

  Use iso_fortran_env, Only : Int32, Int64
  Use vof_plic
  Use vof_normals
  Use vof_advect
  Implicit None

  Real(Int64), Parameter :: pi = 3.14159265358979323846d0, omega = 2d0*pi, kwav = 2d0*pi
  Integer(Int32), Parameter :: n = 64, n1 = n + 2, n2 = n + 2, n3 = 5, nper = 5
  Real(Int64), Parameter :: h = 1d0/n, amp_cells = 3d0
  Real(Int64) :: hh(n1)
  Real(Int64) :: wdts(3) = (/ 0.3d0, 0.15d0, 0.075d0 /)
  Real(Int64), Allocatable :: cref(:,:,:), cs(:,:,:)
  Real(Int64) :: dt, e(7), prev(7), v0, v1, cmin, cmax, co, cl, cotot(7), rate
  Integer(Int32) :: iw, ns, is
  Integer(Int64) :: nint_

  hh = h
  Allocate( cref(0:n1+1,0:n2+1,0:n3+1), cs(0:n1+1,0:n2+1,0:n3+1) )
  prev = 0d0
  Write(*,'(a)') ' w*dt   Co     L1 error vs converged reference:'
  Write(*,'(a)') '                      S1        S2        S3        S4        S5        S6        S7'
  Do iw = 1, 3
     ns = Nint(nper*2d0*pi/wdts(iw))
     dt = Real(nper,Int64)/ns
     Call init_c(cref)
     Call run(8, ns, dt, cref, co)
     Do is = 1, 7
        Call init_c(cs)
        Call run(is, ns, dt, cs, cotot(is))
        e(is) = l1(cs, cref)
     End Do
     Write(*,'(f6.3,f6.2,12x,7es10.2)') wdts(iw), Maxval(cotot), e
     If ( iw > 1 ) Write(*,'(a,3f7.2)') ' order S3,S5,S7:', Log(prev(3)/e(3))/Log(2d0), &
          Log(prev(5)/e(5))/Log(2d0), Log(prev(7)/e(7))/Log(2d0)
     prev = e
  End Do
  Call vof_local_stats(cs, n1, n2, n3, hh, hh, hh(1:n3), v1, cmin, cmax, nint_)
  Write(*,'(a,es10.2,a,es10.2)') 'last run: C range excess ', Max(-cmin, cmax-1d0), '   interface cells ', Real(nint_,Int64)

Contains

  Function psi(x, y, t) Result(p)

    Real(Int64), Intent(In) :: x, y, t
    Real(Int64) :: p, a

    a = amp_cells*h*omega/kwav
    p = a*Cos(kwav*x - omega*t)*Sin(pi*y)

  End Function psi


  Subroutine vel(t, U, V, W)

    Real(Int64), Intent(In)  :: t
    Real(Int64), Intent(Out) :: U(n1-1,n2,n3), V(n1,n2-1,n3), W(n1,n2,n3-1)
    Integer(Int32) :: i, j

    W = 0d0
    Do j = 1, n2
       Do i = 1, n1-1
          U(i,j,:) = ( psi((i-1)*h, (j-1)*h, t) - psi((i-1)*h, (j-2)*h, t) )/h
       End Do
    End Do
    Do j = 1, n2-1
       Do i = 1, n1
          V(i,j,:) = -( psi((i-1)*h, (j-1)*h, t) - psi((i-2)*h, (j-1)*h, t) )/h
       End Do
    End Do

  End Subroutine vel


  Subroutine fill_pad_wall(Cp, m1, m2, m3)

    Integer(Int32), Intent(In)    :: m1, m2, m3
    Real(Int64),    Intent(InOut) :: Cp(0:m1+1,0:m2+1,0:m3+1)

    Cp(1,:,:)    = Cp(m1-1,:,:)
    Cp(m1,:,:)   = Cp(2,:,:)
    Cp(0,:,:)    = Cp(m1-2,:,:)
    Cp(m1+1,:,:) = Cp(3,:,:)
    Cp(:,1,:)    = Cp(:,2,:)
    Cp(:,0,:)    = Cp(:,3,:)
    Cp(:,m2,:)   = Cp(:,m2-1,:)
    Cp(:,m2+1,:) = Cp(:,m2-2,:)
    Cp(:,:,1)    = Cp(:,:,m3-1)
    Cp(:,:,m3)   = Cp(:,:,2)
    Cp(:,:,0)    = Cp(:,:,m3-2)
    Cp(:,:,m3+1) = Cp(:,:,3)

  End Subroutine fill_pad_wall


  !> flat interface y = 0.5 (liquid below), exact cell fractions
  Subroutine init_c(Cp)

    Real(Int64), Intent(Out) :: Cp(0:n1+1,0:n2+1,0:n3+1)
    Integer(Int32) :: j

    Do j = 0, n2+1
       Cp(:,j,:) = Min(1d0, Max(0d0, (0.5d0 - (j-2)*h)/h))
    End Do

  End Subroutine init_c


  Function l1(a, b) Result(e)

    Real(Int64), Intent(In) :: a(0:n1+1,0:n2+1,0:n3+1), b(0:n1+1,0:n2+1,0:n3+1)
    Real(Int64) :: e

    e = Sum(Abs(a(2:n1-1,2:n2-1,2:n3-1) - b(2:n1-1,2:n2-1,2:n3-1)))*h**3

  End Function l1


  !> strategy 1-6 as described above; 7 = reference (8 sub-steps of dt/8, exact midpoint)
  Subroutine run(strat, ns, dt, Cp, comax)

    Integer(Int32), Intent(In)    :: strat, ns
    Real(Int64),    Intent(In)    :: dt
    Real(Int64),    Intent(InOut) :: Cp(0:n1+1,0:n2+1,0:n3+1)
    Real(Int64),    Intent(Out)   :: comax
    Real(Int64) :: U(n1-1,n2,n3), V(n1,n2-1,n3), W(n1,n2,n3-1), U0(n1-1,n2,n3), V0(n1,n2-1,n3)
    Real(Int64) :: Um(n1-1,n2,n3), Vm(n1,n2-1,n3), U1(n1-1,n2,n3), V1(n1,n2-1,n3), U2(n1-1,n2,n3), V2(n1,n2-1,n3)
    Real(Int64) :: t0, sub
    Integer(Int32) :: step, s, m, istep
    Real(Int64), Parameter :: tt(0:3) = (/ 0d0, 8d0/15d0, 2d0/3d0, 1d0 /)

    comax = 0d0
    istep = 0
    Call fill_pad_wall(Cp, n1, n2, n3)
    Do step = 1, ns
       t0 = (step-1)*dt
       Select Case(strat)
       Case(1)
          Call vel(t0 + 0.5d0*dt, U, V, W);  Call adv(Cp, U, V, W, dt, istep, comax)
       Case(2)
          Call vel(t0, U, V, W);  Call adv(Cp, U, V, W, dt, istep, comax)
       Case(3)
          Call vel(t0 + tt(1)*dt, U, V, W);  Call adv(Cp, U, V, W, dt, istep, comax)
       Case(4)
          Call vel(t0, U, V, W)
          Call vel(t0 - dt, U0, V0, W)
          If ( step == 1 ) Then
             Call vel(t0 + 0.5d0*dt, U, V, W)
          Else
             U = 1.5d0*U - 0.5d0*U0;  V = 1.5d0*V - 0.5d0*V0
          End If
          Call adv(Cp, U, V, W, dt, istep, comax)
       Case(5)
          Do s = 1, 3
             Call vel(t0 + tt(s-1)*dt, U, V, W);  Call adv(Cp, U, V, W, dt*(tt(s) - tt(s-1)), istep, comax)
          End Do
       Case(6)
          Do s = 1, 3
             Call vel(t0 + 0.5d0*(tt(s-1) + tt(s))*dt, U, V, W);  Call adv(Cp, U, V, W, dt*(tt(s) - tt(s-1)), istep, comax)
          End Do
       Case(7)
          Call vel(t0, U0, V0, W)
          Call vel(t0 - dt/3d0, Um, Vm, W)
          Call vel(t0 + tt(1)*dt, U1, V1, W)
          Call vel(t0 + tt(2)*dt, U2, V2, W)
          ! sub-interval midpoints 4/15, 3/5, 5/6 from (u^n, u^(2)_prev), (u^n, u^(1)), (u^(1), u^(2))
          U = U0 + (4d0/15d0)*(U0 - Um)/(1d0/3d0);  V = V0 + (4d0/15d0)*(V0 - Vm)/(1d0/3d0)
          Call adv(Cp, U, V, W, dt*tt(1), istep, comax)
          U = U1 + (0.6d0 - tt(1))*(U1 - U0)/tt(1);  V = V1 + (0.6d0 - tt(1))*(V1 - V0)/tt(1)
          Call adv(Cp, U, V, W, dt*(tt(2) - tt(1)), istep, comax)
          U = U2 + (5d0/6d0 - tt(2))*(U2 - U1)/(tt(2) - tt(1));  V = V2 + (5d0/6d0 - tt(2))*(V2 - V1)/(tt(2) - tt(1))
          Call adv(Cp, U, V, W, dt*(tt(3) - tt(2)), istep, comax)
       Case(8)
          sub = dt/8d0
          Do m = 1, 8
             Call vel(t0 + (m-0.5d0)*sub, U, V, W);  Call adv(Cp, U, V, W, sub, istep, comax)
          End Do
       End Select
    End Do

  End Subroutine run


  Subroutine adv(Cp, U, V, W, dts, istep, comax)

    Real(Int64),    Intent(InOut) :: Cp(0:n1+1,0:n2+1,0:n3+1)
    Real(Int64),    Intent(In)    :: U(n1-1,n2,n3), V(n1,n2-1,n3), W(n1,n2,n3-1), dts
    Integer(Int32), Intent(InOut) :: istep
    Real(Int64),    Intent(InOut) :: comax
    Real(Int64) :: co, cl

    Call vof_advect_step(Cp, n1, n2, n3, U, V, W, hh, hh, hh(1:n3), dts, istep, VOF_SCHEME_YOUNGS, fill_pad_wall, co, cl)
    istep = istep + 1
    comax = Max(comax, co)

  End Subroutine adv

End Program spike_c_timeint
