!> Prescribed-velocity validation of src/vof_advect.f90 (serial, periodic padding, no solver): oblique-band translation against the
!  exact geometric solution, Zalesak slotted-disk rotation and the LeVeque vortex-reversal test (forward/reverse distortion), with
!  volume conservation, boundedness and convergence reported.
Program test_vof_advect

  Use iso_fortran_env, Only : Int32, Int64
  Use vof_plic
  Use vof_normals
  Use vof_advect

  Implicit None

  Real(Int64), Parameter :: pi = 3.14159265358979323846d0
  Integer(Int32) :: n, iscm, ic, nlev
  Real(Int64) :: max_drift = 0d0, max_excess = 0d0, min_order = 1d3
  Real(Int64) :: e_prev, e_now, rate
  Real(Int64), Parameter :: cos_list(5) = (/ 0.05d0, 0.2d0, 0.5d0, 0.8d0, 0.95d0 /)
  Logical :: ok

  ok = .True.
  nlev = 2
  If ( Command_argument_count() > 0 ) nlev = 3

  Write(*,'(a)') '--- Case 1: oblique band translation (exact reference), 3-D velocity, Co ~ 0.4'
  Do iscm = 1, 2
     e_prev = 0d0
     Write(*,'(a,i0)') 'normal scheme ', iscm
     Do ic = 1, nlev
        n = 16*2**ic
        Call run_translation(n, iscm, e_now)
        rate = 0d0
        If ( ic > 1 ) rate = Log(e_prev/e_now)/Log(2d0)
        If ( ic > 1 ) min_order = Min(min_order, rate)
        Write(*,'(a,i4,a,es10.3,a,f6.2)') '  N=', n, '  L1 err=', e_now, '  order=', rate
        e_prev = e_now
     End Do
  End Do

  Write(*,'(a)') '--- Case 2: Zalesak slotted disk, one revolution'
  Do iscm = 1, 2
     e_prev = 0d0
     Write(*,'(a,i0)') 'normal scheme ', iscm
     Do ic = 1, nlev
        n = 50*2**(ic-1)
        Call run_zalesak(n, iscm, e_now)
        rate = 0d0
        If ( ic > 1 ) rate = Log(e_prev/e_now)/Log(2d0)
        If ( ic > 1 ) min_order = Min(min_order, rate)
        Write(*,'(a,i4,a,es10.3,a,f6.2)') '  N=', n, '  L1 err=', e_now, '  order=', rate
        e_prev = e_now
     End Do
  End Do

  Write(*,'(a)') '--- Case 3: vortex reversal (forward/reverse), T=2 and T=8'
  Do iscm = 1, 2
     Write(*,'(a,i0)') 'normal scheme ', iscm
     Do ic = 1, nlev-1
        n = 64*2**(ic-1)
        Call run_vortex(n, iscm, 2d0)
        Call run_vortex(n, iscm, 8d0)
     End Do
  End Do

  Write(*,'(a)') '--- Case 4: disk translation, grid-orientation sensitivity (N=64, Co=0.4) and CFL sensitivity (45 deg)'
  Do ic = 0, 4
     Call run_disk_translate(64, 1, Real(ic,Int64)*11.25d0, 0.4d0, e_now)
  End Do
  Do ic = 1, 5
     Call run_disk_translate(64, 1, 45d0, cos_list(ic), e_now)
  End Do

  Write(*,'(a)') '--- Case 5: THINC (diffuse-interface) fluxes, vortex reversal and disk translation: conservation and boundedness'
  vof_flux_scheme = 1
  Call run_vortex(64, 1, 2d0)
  Call run_disk_translate(64, 1, 22.5d0, 0.4d0, e_now)
  vof_flux_scheme = 0

  If ( max_drift > 1d-12 .Or. max_excess > 1d-12 .Or. min_order < 1.8d0 ) ok = .False.
  Write(*,'(a,es9.2,a,es9.2,a,f5.2)') 'max |volume drift|=', max_drift, '  max C excess=', max_excess, '  min order=', min_order
  If ( ok ) Then
     Write(*,'(a)') 'PASS'
  Else
     Write(*,'(a)') 'FAIL'
     Stop 1
  End If

Contains

  Subroutine fill_pad_periodic(Cp, n1, n2, n3)

    Integer(Int32), Intent(In)    :: n1, n2, n3
    Real(Int64),    Intent(InOut) :: Cp(0:n1+1,0:n2+1,0:n3+1)

    Cp(1,:,:)    = Cp(n1-1,:,:)
    Cp(n1,:,:)   = Cp(2,:,:)
    Cp(0,:,:)    = Cp(n1-2,:,:)
    Cp(n1+1,:,:) = Cp(3,:,:)
    Cp(:,1,:)    = Cp(:,n2-1,:)
    Cp(:,n2,:)   = Cp(:,2,:)
    Cp(:,0,:)    = Cp(:,n2-2,:)
    Cp(:,n2+1,:) = Cp(:,3,:)
    Cp(:,:,1)    = Cp(:,:,n3-1)
    Cp(:,:,n3)   = Cp(:,:,2)
    Cp(:,:,0)    = Cp(:,:,n3-2)
    Cp(:,:,n3+1) = Cp(:,:,3)

  End Subroutine fill_pad_periodic


  !> Exact cell fraction of the periodic band a <= (x+y) mod 1 < b (normal (1,1,0)), cell [x0,x0+h]x[y0,y0+h]
  Function band_fraction(x0, y0, h, a, b) Result(f)

    Real(Int64), Intent(In) :: x0, y0, h, a, b
    Real(Int64) :: f
    Integer(Int32) :: m

    f = 0d0
    Do m = -3, 3
       f = f + plic_box_fraction(h, h, 0d0, b + m - x0 - y0) - plic_box_fraction(h, h, 0d0, a + m - x0 - y0)
    End Do

  End Function band_fraction


  Subroutine run_translation(n, scheme, err)

    Integer(Int32), Intent(In)  :: n, scheme
    Real(Int64),    Intent(Out) :: err

    Integer(Int32) :: n1, n2, n3, i, j, k, istep, nsteps
    Real(Int64) :: h, dt, uu, vv, ww, tend, co, cl, cotot, cltot, x0, y0, vol, v0, v1, cmin, cmax, vel(3)
    Real(Int64), Allocatable :: Cp(:,:,:), U(:,:,:), V(:,:,:), W(:,:,:), hh(:)
    Integer(Int64) :: nint

    n1 = n + 2;  n2 = n + 2;  n3 = 6
    h = 1d0/n
    Allocate( Cp(0:n1+1,0:n2+1,0:n3+1), U(n1-1,n2,n3), V(n1,n2-1,n3), W(n1,n2,n3-1), hh(Max(n1,n3)) )
    hh = h
    ! band (x+y) mod 1 in [0.25, 0.6): velocity (0.7,0.3,-0.4) shifts x+y by 1.0 per unit time -> after t=0.5 by 0.5
    uu = 0.7d0;  vv = 0.3d0;  ww = -0.4d0
    U = uu;  V = vv;  W = ww
    vel = (/ uu, vv, ww /)
    tend = 0.5d0
    dt = 0.4d0*h/Maxval(Abs(vel))
    nsteps = Ceiling(tend/dt)
    dt = tend/nsteps

    Cp = 0d0
    Do k = 1, n3
       Do j = 1, n2
          Do i = 1, n1
             Cp(i,j,k) = band_fraction((i-2)*h, (j-2)*h, h, 0.25d0, 0.6d0)
          End Do
       End Do
    End Do
    Call fill_pad_periodic(Cp, n1, n2, n3)
    Call vof_local_stats(Cp, n1, n2, n3, hh(1:n1), hh(1:n2), hh(1:n3), v0, cmin, cmax, nint)
    cotot = 0d0;  cltot = 0d0
    Do istep = 0, nsteps-1
       Call vof_advect_step(Cp, n1, n2, n3, U, V, W, hh(1:n1), hh(1:n2), hh(1:n3), dt, istep, scheme, fill_pad_periodic, co, cl)
       cotot = Max(cotot, co);  cltot = cltot + cl
    End Do
    Call vof_local_stats(Cp, n1, n2, n3, hh(1:n1), hh(1:n2), hh(1:n3), v1, cmin, cmax, nint)

    err = 0d0
    vol = h**3
    Do k = 2, n3-1
       Do j = 2, n2-1
          Do i = 2, n1-1
             x0 = (i-2)*h - uu*tend;  y0 = (j-2)*h - vv*tend
             err = err + Abs( Cp(i,j,k) - band_fraction(x0, y0, h, 0.25d0, 0.6d0) )*vol
          End Do
       End Do
    End Do
    max_drift = Max(max_drift, Abs((v1-v0)/v0));  max_excess = Max(max_excess, -cmin, cmax-1d0)
    If ( n == 64 ) Write(*,'(a,es9.2,a,es9.2,a,es9.2,a,f5.2)') '      vol drift=', (v1-v0)/v0, '  clip=', cltot, &
         '  C range excess=', Max(-cmin, cmax-1d0), '  Co=', cotot

  End Subroutine run_translation


  !> Disk (centre (0.5,0.75), R=0.15) minus slot (width 0.05, down to y=0.7): sub-sampled volume fraction of cell (x0,y0,h)
  Function zalesak_fraction(x0, y0, h) Result(f)

    Real(Int64), Intent(In) :: x0, y0, h
    Real(Int64) :: f, x, y
    Integer(Int32) :: a, b
    Integer(Int32), Parameter :: ns = 16

    f = 0d0
    Do b = 1, ns
       y = y0 + (b-0.5d0)*h/ns
       Do a = 1, ns
          x = x0 + (a-0.5d0)*h/ns
          If ( (x-0.5d0)**2 + (y-0.75d0)**2 <= 0.15d0**2 ) Then
             If ( .Not. ( Abs(x-0.5d0) <= 0.025d0 .And. y <= 0.8250d0 ) ) f = f + 1d0
          End If
       End Do
    End Do
    f = f/(ns*ns)

  End Function zalesak_fraction


  Subroutine run_zalesak(n, scheme, err)

    Integer(Int32), Intent(In)  :: n, scheme
    Real(Int64),    Intent(Out) :: err

    Integer(Int32) :: n1, n2, n3, i, j, k, istep, nsteps
    Real(Int64) :: h, dt, om, co, cl, cltot, cotot, v0, v1, cmin, cmax, xf, yl, yh, xn, yf
    Real(Int64), Allocatable :: Cp(:,:,:), C0(:,:,:), U(:,:,:), V(:,:,:), W(:,:,:), hh(:)
    Integer(Int64) :: nint

    n1 = n + 2;  n2 = n + 2;  n3 = 5
    h = 1d0/n
    om = 2d0*pi
    Allocate( Cp(0:n1+1,0:n2+1,0:n3+1), C0(0:n1+1,0:n2+1,0:n3+1), U(n1-1,n2,n3), V(n1,n2-1,n3), W(n1,n2,n3-1), hh(Max(n1,n3)) )
    hh = h
    W = 0d0
    ! psi = -om/2 ((x-.5)^2+(y-.5)^2): u = dpsi/dy at corners, v = -dpsi/dx
    Do k = 1, n3
       Do j = 1, n2
          Do i = 1, n1-1
             xf = (i-1)*h;  yl = (j-2)*h;  yh = (j-1)*h
             U(i,j,k) = ( psi(xf,yh) - psi(xf,yl) )/h
          End Do
       End Do
       Do j = 1, n2-1
          Do i = 1, n1
             yf = (j-1)*h
             xn = (i-2)*h
             V(i,j,k) = -( psi(xn+h,yf) - psi(xn,yf) )/h
          End Do
       End Do
    End Do

    Cp = 0d0
    Do k = 1, n3
       Do j = 1, n2
          Do i = 1, n1
             Cp(i,j,k) = zalesak_fraction((i-2)*h, (j-2)*h, h)
          End Do
       End Do
    End Do
    Call fill_pad_periodic(Cp, n1, n2, n3)
    C0 = Cp
    Call vof_local_stats(Cp, n1, n2, n3, hh(1:n1), hh(1:n2), hh(1:n3), v0, cmin, cmax, nint)

    dt = 0.4d0*h/(om*0.7072d0)
    nsteps = Ceiling(1d0/dt)
    dt = 1d0/nsteps
    cotot = 0d0;  cltot = 0d0
    Do istep = 0, nsteps-1
       Call vof_advect_step(Cp, n1, n2, n3, U, V, W, hh(1:n1), hh(1:n2), hh(1:n3), dt, istep, scheme, fill_pad_periodic, co, cl)
       cotot = Max(cotot, co);  cltot = cltot + cl
    End Do
    Call vof_local_stats(Cp, n1, n2, n3, hh(1:n1), hh(1:n2), hh(1:n3), v1, cmin, cmax, nint)

    err = 0d0
    Do k = 2, n3-1
       Do j = 2, n2-1
          Do i = 2, n1-1
             err = err + Abs(Cp(i,j,k) - C0(i,j,k))*h**3
          End Do
       End Do
    End Do
    max_drift = Max(max_drift, Abs((v1-v0)/v0));  max_excess = Max(max_excess, -cmin, cmax-1d0)
    If ( n == 100 ) Write(*,'(a,es9.2,a,es9.2,a,es9.2,a,f5.2,a,i0)') '      vol drift=', (v1-v0)/v0, '  clip=', cltot, &
         '  C range excess=', Max(-cmin, cmax-1d0), '  Co=', cotot, '  interface cells=', nint

  End Subroutine run_zalesak


  Function psi(x, y) Result(p)

    Real(Int64), Intent(In) :: x, y
    Real(Int64) :: p

    p = -pi*( (x-0.5d0)**2 + (y-0.5d0)**2 )

  End Function psi


  Subroutine run_vortex(n, scheme, tper)

    Integer(Int32), Intent(In) :: n, scheme
    Real(Int64),    Intent(In) :: tper

    Integer(Int32) :: n1, n2, n3, i, j, k, istep, nsteps
    Real(Int64) :: h, dt, t, co, cl, cltot, cotot, v0, v1, cmin, cmax, err, xf, yl, yh, xn, yf, cs
    Real(Int64), Allocatable :: Cp(:,:,:), C0(:,:,:), U(:,:,:), V(:,:,:), W(:,:,:), hh(:)
    Integer(Int64) :: nint

    n1 = n + 2;  n2 = n + 2;  n3 = 5
    h = 1d0/n
    Allocate( Cp(0:n1+1,0:n2+1,0:n3+1), C0(0:n1+1,0:n2+1,0:n3+1), U(n1-1,n2,n3), V(n1,n2-1,n3), W(n1,n2,n3-1), hh(Max(n1,n3)) )
    hh = h
    W = 0d0
    Cp = 0d0
    Do k = 1, n3
       Do j = 1, n2
          Do i = 1, n1
             Cp(i,j,k) = disk_fraction((i-2)*h, (j-2)*h, h)
          End Do
       End Do
    End Do
    Call fill_pad_periodic(Cp, n1, n2, n3)
    C0 = Cp
    Call vof_local_stats(Cp, n1, n2, n3, hh(1:n1), hh(1:n2), hh(1:n3), v0, cmin, cmax, nint)

    dt = 0.4d0*h
    nsteps = Ceiling(tper/dt)
    dt = tper/nsteps
    cotot = 0d0;  cltot = 0d0
    Do istep = 0, nsteps-1
       t = (istep + 0.5d0)*dt
       cs = Cos(pi*t/tper)
       Do k = 1, n3
          Do j = 1, n2
             Do i = 1, n1-1
                xf = (i-1)*h;  yl = (j-2)*h;  yh = (j-1)*h
                U(i,j,k) = cs*( vpsi(xf,yh) - vpsi(xf,yl) )/h
             End Do
          End Do
          Do j = 1, n2-1
             Do i = 1, n1
                yf = (j-1)*h
                xn = (i-2)*h
                V(i,j,k) = -cs*( vpsi(xn+h,yf) - vpsi(xn,yf) )/h
             End Do
          End Do
       End Do
       Call vof_advect_step(Cp, n1, n2, n3, U, V, W, hh(1:n1), hh(1:n2), hh(1:n3), dt, istep, scheme, fill_pad_periodic, co, cl)
       cotot = Max(cotot, co);  cltot = cltot + cl
    End Do
    Call vof_local_stats(Cp, n1, n2, n3, hh(1:n1), hh(1:n2), hh(1:n3), v1, cmin, cmax, nint)
    err = 0d0
    Do k = 2, n3-1
       Do j = 2, n2-1
          Do i = 2, n1-1
             err = err + Abs(Cp(i,j,k) - C0(i,j,k))*h**3
          End Do
       End Do
    End Do
    max_drift = Max(max_drift, Abs((v1-v0)/v0));  max_excess = Max(max_excess, -cmin, cmax-1d0)
    Write(*,'(a,i4,a,f4.1,a,es10.3,a,es9.2,a,es9.2,a,f5.2)') '  N=', n, ' T=', tper, '  L1 err=', err, '  vol drift=', &
         (v1-v0)/v0, '  C excess=', Max(-cmin, cmax-1d0), '  Co=', cotot

  End Subroutine run_vortex


  !> psi = sin^2(pi x) sin^2(pi y)/pi (the vortex-reversal streamfunction at unit amplitude)
  Function vpsi(x, y) Result(p)

    Real(Int64), Intent(In) :: x, y
    Real(Int64) :: p

    p = Sin(pi*x)**2*Sin(pi*y)**2/pi

  End Function vpsi


  !> Disk of radius 0.15 translated at speed 1 along angle theta (deg) for t=0.5, compared with the sub-sampled disk at its exact
  !  final position; co is the target Courant number of the largest velocity component
  Subroutine run_disk_translate(n, scheme, theta, co, err)

    Integer(Int32), Intent(In)  :: n, scheme
    Real(Int64),    Intent(In)  :: theta, co
    Real(Int64),    Intent(Out) :: err

    Integer(Int32) :: n1, n2, n3, i, j, k, istep, nsteps
    Real(Int64) :: h, dt, uu, vv, cobs, cl, cltot, cotot, v0, v1, cmin, cmax, tend, cx, cy
    Real(Int64), Allocatable :: Cp(:,:,:), U(:,:,:), V(:,:,:), W(:,:,:), hh(:)
    Integer(Int64) :: nint
    Integer(Int32) :: mx, my
    Real(Int64) :: cref

    n1 = n + 2;  n2 = n + 2;  n3 = 5
    h = 1d0/n
    Allocate( Cp(0:n1+1,0:n2+1,0:n3+1), U(n1-1,n2,n3), V(n1,n2-1,n3), W(n1,n2,n3-1), hh(Max(n1,n3)) )
    hh = h
    uu = Cos(theta*pi/180d0);  vv = Sin(theta*pi/180d0)
    U = uu;  V = vv;  W = 0d0
    tend = 0.5d0
    dt = co*h/Max(Abs(uu), Abs(vv))
    nsteps = Ceiling(tend/dt)
    dt = tend/nsteps
    Do k = 1, n3
       Do j = 1, n2
          Do i = 1, n1
             Cp(i,j,k) = disk_fraction((i-2)*h, (j-2)*h, h)
          End Do
       End Do
    End Do
    Call fill_pad_periodic(Cp, n1, n2, n3)
    Call vof_local_stats(Cp, n1, n2, n3, hh(1:n1), hh(1:n2), hh(1:n3), v0, cmin, cmax, nint)
    cotot = 0d0;  cltot = 0d0
    Do istep = 0, nsteps-1
       Call vof_advect_step(Cp, n1, n2, n3, U, V, W, hh(1:n1), hh(1:n2), hh(1:n3), dt, istep, scheme, fill_pad_periodic, cobs, cl)
       cotot = Max(cotot, cobs);  cltot = cltot + cl
    End Do
    Call vof_local_stats(Cp, n1, n2, n3, hh(1:n1), hh(1:n2), hh(1:n3), v1, cmin, cmax, nint)
    cx = uu*tend;  cy = vv*tend
    err = 0d0
    Do k = 2, n3-1
       Do j = 2, n2-1
          Do i = 2, n1-1
             cref = 0d0
             Do mx = -1, 1
                Do my = -1, 1
                   cref = cref + disk_fraction((i-2)*h - cx + mx, (j-2)*h - cy + my, h)
                End Do
             End Do
             err = err + Abs( Cp(i,j,k) - cref )*h**3
          End Do
       End Do
    End Do
    max_drift = Max(max_drift, Abs((v1-v0)/v0));  max_excess = Max(max_excess, -cmin, cmax-1d0)
    Write(*,'(a,f6.2,a,f5.2,a,es10.3,a,es9.2)') '  theta=', theta, ' Co=', cotot, '  L1 err=', err, '  vol drift=', (v1-v0)/v0

  End Subroutine run_disk_translate


  Function disk_fraction(x0, y0, h) Result(f)

    Real(Int64), Intent(In) :: x0, y0, h
    Real(Int64) :: f, x, y
    Integer(Int32) :: a, b
    Integer(Int32), Parameter :: ns = 16

    f = 0d0
    Do b = 1, ns
       y = y0 + (b-0.5d0)*h/ns
       Do a = 1, ns
          x = x0 + (a-0.5d0)*h/ns
          If ( (x-0.5d0)**2 + (y-0.75d0)**2 <= 0.15d0**2 ) f = f + 1d0
       End Do
    End Do
    f = f/(ns*ns)

  End Function disk_fraction

End Program test_vof_advect
