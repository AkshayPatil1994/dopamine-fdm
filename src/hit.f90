! Homogeneous isotropic turbulence in a triply periodic box: random-phase solenoidal initial condition
! (ic_type=8), band-limited deterministic forcing (hit_forcing=1) and TKE/dissipation diagnostics (hit_stats=1).
!
! Every periodic direction has n_global-2 distinct points and period (n_global-2)*d (see apply_periodic_bc_*).
! Fourier coefficients are taken at the true staggered positions of each velocity component, and the
! solenoidal projection uses the modified wavenumber 2*sin(k*d/2)/d, so a projected field has zero discrete
! divergence on the MAC grid. The transforms are separable direct sums over the few forced modes
Module hit

  Use iso_fortran_env, Only : Int32, Int64
  Use global
  Use mpi
  Use decomp, Only : x_periodic_partner, z_periodic_partner

  Implicit None

  Real(Int64), Allocatable, Dimension(:,:,:) :: fU, fV, fW   ! forcing at the U, V, W points

  ! phase tables e^{+i k_m x} at the face (f) and centre (c) points of each direction, masks of the
  ! distinct points this rank owns, and modified wavenumbers km(1:3, m)
  Type ptabs
     Complex(Int64), Allocatable :: ex_f(:,:), ex_c(:,:), ey_f(:,:), ey_c(:,:), ez_f(:,:), ez_c(:,:)
     Real   (Int64), Allocatable :: mx_f(:), mx_c(:), my_f(:), my_c(:), mz_f(:), mz_c(:), km(:,:)
     Integer(Int32) :: nmax
  End Type

  ! forcing tables as separate real cos/sin arrays, device-resident under OpenACC
  Type rtabs
     Real(Int64), Allocatable :: cx_f(:,:), sx_f(:,:), cx_c(:,:), sx_c(:,:), cy_f(:,:), sy_f(:,:), cy_c(:,:), sy_c(:,:), &
                                 cz_f(:,:), sz_f(:,:), cz_c(:,:), sz_c(:,:)
  End Type

  Type(ptabs), Save :: TF, TM
  Type(rtabs), Save :: RT
  Real(Int64), Allocatable, Save :: Wa_r(:,:,:), Wa_i(:,:,:), Wb_r(:,:,:), Wb_i(:,:,:), Cc_r(:,:,:), Cc_i(:,:,:)
  Logical,     Save :: hit_ready = .False.
  Logical,     Allocatable, Save :: inband(:,:,:)
  Real(Int64), Save :: n_tot = 1d0
  Logical,     Save :: stats_started = .False.
  Integer(Int64), Save :: rng_state = 1_Int64

Contains

  Subroutine build_tables(nmax, T)

    Integer(Int32), Intent(In)  :: nmax
    Type(ptabs),    Intent(Out) :: T

    Logical        :: is_first, is_last
    Integer(Int32) :: partner, i, m
    Real(Int64)    :: d(3), P(3), kk, cnt(3), cnt_g(3)

    d(1) = x_global(2) - x_global(1)
    d(2) = y_global(2) - y_global(1)
    d(3) = z_global(2) - z_global(1)
    P(1) = Real(nx_global-2, Int64)*d(1)
    P(2) = Real(ny_global-2, Int64)*d(2)
    P(3) = Real(nz_global-2, Int64)*d(3)
    T%nmax = nmax

    Allocate ( T%ex_f(nx, -nmax:nmax), T%ex_c(nxg, -nmax:nmax), T%ey_f(ny, -nmax:nmax), T%ey_c(nyg, -nmax:nmax), &
               T%ez_f(nz, -nmax:nmax), T%ez_c(nzg, -nmax:nmax) )
    Allocate ( T%mx_f(nx), T%mx_c(nxg), T%my_f(ny), T%my_c(nyg), T%mz_f(nz), T%mz_c(nzg), T%km(3,-nmax:nmax) )

    Do m = -nmax, nmax
       Do i = 1, 3
          kk = 2d0*pi*Real(m, Int64)/P(i)
          T%km(i,m) = 2d0*Sin(0.5d0*kk*d(i))/d(i)
       End Do
       kk = 2d0*pi*Real(m, Int64)/P(1)
       Do i = 1, nx
          T%ex_f(i,m) = Exp( Cmplx(0d0, kk*Real(i1_global(myid)+i-2, Int64)*d(1), Int64) )
       End Do
       Do i = 1, nxg
          T%ex_c(i,m) = Exp( Cmplx(0d0, kk*(Real(ig1_global(myid)+i-2, Int64)-1.5d0)*d(1), Int64) )
       End Do
       kk = 2d0*pi*Real(m, Int64)/P(2)
       Do i = 1, ny
          T%ey_f(i,m) = Exp( Cmplx(0d0, kk*Real(i-1, Int64)*d(2), Int64) )
       End Do
       Do i = 1, nyg
          T%ey_c(i,m) = Exp( Cmplx(0d0, kk*(Real(i-1, Int64)-0.5d0)*d(2), Int64) )
       End Do
       kk = 2d0*pi*Real(m, Int64)/P(3)
       Do i = 1, nz
          T%ez_f(i,m) = Exp( Cmplx(0d0, kk*Real(k1_global(myid)+i-2, Int64)*d(3), Int64) )
       End Do
       Do i = 1, nzg
          T%ez_c(i,m) = Exp( Cmplx(0d0, kk*(Real(k1_global(myid)+i-1, Int64)-1.5d0)*d(3), Int64) )
       End Do
    End Do

    ! distinct points: faces 2..n-1; centres 2..ng-1, except the last rank's last centre (duplicate of centre 2)
    T%mx_f = 0d0;  T%mx_f(2:nx-1)  = 1d0
    T%my_f = 0d0;  T%my_f(2:ny-1)  = 1d0
    T%mz_f = 0d0;  T%mz_f(2:nz-1)  = 1d0
    T%mx_c = 0d0;  T%mx_c(2:nxg-1) = 1d0
    T%my_c = 0d0;  T%my_c(2:nyg-2) = 1d0
    T%mz_c = 0d0;  T%mz_c(2:nzg-1) = 1d0
    Call x_periodic_partner(is_first, is_last, partner)
    If ( is_last ) T%mx_c(nxg-1) = 0d0
    Call z_periodic_partner(is_first, is_last, partner)
    If ( is_last ) T%mz_c(nzg-1) = 0d0

    ! the ranks must tile the distinct points exactly once
    cnt(1) = Sum(T%mx_f)*Sum(T%my_c)*Sum(T%mz_c)
    cnt(2) = Sum(T%mx_c)*Sum(T%my_f)*Sum(T%mz_c)
    cnt(3) = Sum(T%mx_c)*Sum(T%my_c)*Sum(T%mz_f)
    Call MPI_Allreduce(cnt, cnt_g, 3, MPI_real8, MPI_SUM, MPI_COMM_WORLD, ierr)
    n_tot = Real(nx_global-2, Int64)*Real(ny_global-2, Int64)*Real(nz_global-2, Int64)
    If ( Any(Abs(cnt_g - n_tot) > 0.5d0) ) Then
       If ( myid == 0 ) Write(*,'(A,3F14.0,A,F14.0)') ' ERROR: HIT point count ', cnt_g, ' /= ', n_tot
       Call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    End If

  End Subroutine build_tables

  ! c(mx,my,mz) = (1/n_tot) sum_points F e^{-i k.x} over the distinct points, (cr,ci) real and imaginary parts.
  ! Three separable stages on the device; only the K^3 result is reduced over the ranks on the host
  Subroutine force_analyse(F, n1, n2, n3, c1, s1, m1, c2, s2, m2, c3, s3, m3, nm, cr, ci)

    Integer(Int32), Intent(In)  :: n1, n2, n3, nm
    Real   (Int64), Intent(In)  :: F(n1,n2,n3), m1(n1), m2(n2), m3(n3)
    Real   (Int64), Intent(In)  :: c1(n1,-nm:nm), s1(n1,-nm:nm), c2(n2,-nm:nm), s2(n2,-nm:nm), c3(n3,-nm:nm), s3(n3,-nm:nm)
    Real   (Int64), Intent(Out) :: cr(-nm:nm,-nm:nm,-nm:nm), ci(-nm:nm,-nm:nm,-nm:nm)

    Integer(Int32) :: i, j, k, mx, my, mz
    Real   (Int64) :: sr, si, w, ar, ai
    Real   (Int64) :: cl(-nm:nm,-nm:nm,-nm:nm,2), cg(-nm:nm,-nm:nm,-nm:nm,2)

    !$acc parallel loop gang collapse(3) present(F,m1,c1,s1,Wa_r,Wa_i)
    Do k = 1, n3
       Do j = 1, n2
          Do mx = -nm, nm
             sr = 0d0;  si = 0d0
             !$acc loop vector reduction(+:sr,si)
             Do i = 1, n1
                w  = m1(i)*F(i,j,k)
                sr = sr + w*c1(i,mx)
                si = si - w*s1(i,mx)
             End Do
             Wa_r(mx,j,k) = sr;  Wa_i(mx,j,k) = si
          End Do
       End Do
    End Do
    !$acc parallel loop collapse(3) present(m2,c2,s2,Wa_r,Wa_i,Wb_r,Wb_i) private(sr,si,ar,ai)
    Do k = 1, n3
       Do my = -nm, nm
          Do mx = -nm, nm
             sr = 0d0;  si = 0d0
             !$acc loop seq
             Do j = 1, n2
                ar = m2(j)*Wa_r(mx,j,k);  ai = m2(j)*Wa_i(mx,j,k)
                sr = sr + ar*c2(j,my) + ai*s2(j,my)
                si = si + ai*c2(j,my) - ar*s2(j,my)
             End Do
             Wb_r(mx,my,k) = sr;  Wb_i(mx,my,k) = si
          End Do
       End Do
    End Do
    !$acc parallel loop collapse(3) present(m3,c3,s3,Wb_r,Wb_i,Cc_r,Cc_i) private(sr,si,ar,ai)
    Do mz = -nm, nm
       Do my = -nm, nm
          Do mx = -nm, nm
             sr = 0d0;  si = 0d0
             !$acc loop seq
             Do k = 1, n3
                ar = m3(k)*Wb_r(mx,my,k);  ai = m3(k)*Wb_i(mx,my,k)
                sr = sr + ar*c3(k,mz) + ai*s3(k,mz)
                si = si + ai*c3(k,mz) - ar*s3(k,mz)
             End Do
             Cc_r(mx,my,mz) = sr;  Cc_i(mx,my,mz) = si
          End Do
       End Do
    End Do
    !$acc update host(Cc_r,Cc_i)
    cl(:,:,:,1) = Cc_r;  cl(:,:,:,2) = Cc_i
    Call MPI_Allreduce(cl, cg, Size(cl), MPI_real8, MPI_SUM, MPI_COMM_WORLD, ierr)
    cr = cg(:,:,:,1)/n_tot;  ci = cg(:,:,:,2)/n_tot

  End Subroutine force_analyse

  ! F(i,j,k) = Re sum_m c(m) e^{i k.x} on the device; (cr,ci) are host arrays copied over here
  Subroutine force_synth(F, n1, n2, n3, c1, s1, c2, s2, c3, s3, nm, cr, ci)

    Integer(Int32), Intent(In)    :: n1, n2, n3, nm
    Real   (Int64), Intent(InOut) :: F(n1,n2,n3)
    Real   (Int64), Intent(In)    :: c1(n1,-nm:nm), s1(n1,-nm:nm), c2(n2,-nm:nm), s2(n2,-nm:nm), c3(n3,-nm:nm), s3(n3,-nm:nm)
    Real   (Int64), Intent(In)    :: cr(-nm:nm,-nm:nm,-nm:nm), ci(-nm:nm,-nm:nm,-nm:nm)

    Integer(Int32) :: i, j, k, mx, my, mz
    Real   (Int64) :: sr, si, tr, ti, s

    Cc_r = cr;  Cc_i = ci
    !$acc update device(Cc_r,Cc_i)
    !$acc parallel loop collapse(3) present(c3,s3,Cc_r,Cc_i,Wb_r,Wb_i) private(sr,si)
    Do k = 1, n3
       Do my = -nm, nm
          Do mx = -nm, nm
             sr = 0d0;  si = 0d0
             !$acc loop seq
             Do mz = -nm, nm
                sr = sr + Cc_r(mx,my,mz)*c3(k,mz) - Cc_i(mx,my,mz)*s3(k,mz)
                si = si + Cc_r(mx,my,mz)*s3(k,mz) + Cc_i(mx,my,mz)*c3(k,mz)
             End Do
             Wb_r(mx,my,k) = sr;  Wb_i(mx,my,k) = si
          End Do
       End Do
    End Do
    !$acc parallel loop collapse(3) present(c2,s2,Wb_r,Wb_i,Wa_r,Wa_i) private(tr,ti)
    Do k = 1, n3
       Do j = 1, n2
          Do mx = -nm, nm
             tr = 0d0;  ti = 0d0
             !$acc loop seq
             Do my = -nm, nm
                tr = tr + Wb_r(mx,my,k)*c2(j,my) - Wb_i(mx,my,k)*s2(j,my)
                ti = ti + Wb_r(mx,my,k)*s2(j,my) + Wb_i(mx,my,k)*c2(j,my)
             End Do
             Wa_r(mx,j,k) = tr;  Wa_i(mx,j,k) = ti
          End Do
       End Do
    End Do
    !$acc parallel loop collapse(3) present(F,c1,s1,Wa_r,Wa_i) private(s)
    Do k = 1, n3
       Do j = 1, n2
          Do i = 1, n1
             s = 0d0
             !$acc loop seq
             Do mx = -nm, nm
                s = s + Wa_r(mx,j,k)*c1(i,mx) - Wa_i(mx,j,k)*s1(i,mx)
             End Do
             F(i,j,k) = s
          End Do
       End Do
    End Do

  End Subroutine force_synth

  ! F(i,j,k) = Re sum_m c(m) e^{i k.x} at all local points (ghost points included)
  Subroutine synthesise(F, n1, n2, n3, e1, e2, e3, nmax, c)

    Integer(Int32), Intent(In)  :: n1, n2, n3, nmax
    Real   (Int64), Intent(Out) :: F(n1,n2,n3)
    Complex(Int64), Intent(In)  :: e1(n1,-nmax:nmax), e2(n2,-nmax:nmax), e3(n3,-nmax:nmax)
    Complex(Int64), Intent(In)  :: c(-nmax:nmax,-nmax:nmax,-nmax:nmax)

    Complex(Int64), Allocatable :: A(:,:,:), B(:,:)
    Complex(Int64) :: s
    Integer(Int32) :: i, j, k, m

    Allocate ( A(-nmax:nmax,-nmax:nmax,n3), B(-nmax:nmax,n2) )
    Do k = 1, n3
       A(:,:,k) = (0d0,0d0)
       Do m = -nmax, nmax
          A(:,:,k) = A(:,:,k) + c(:,:,m)*e3(k,m)
       End Do
    End Do
    Do k = 1, n3
       B = (0d0,0d0)
       Do j = 1, n2
          Do m = -nmax, nmax
             B(:,j) = B(:,j) + A(:,m,k)*e2(j,m)
          End Do
       End Do
       Do j = 1, n2
          Do i = 1, n1
             s = (0d0,0d0)
             Do m = -nmax, nmax
                s = s + B(m,j)*e1(i,m)
             End Do
             F(i,j,k) = Real(s, Int64)
          End Do
       End Do
    End Do
    Deallocate(A, B)

  End Subroutine synthesise

  ! Solenoidal projection (modified wavenumbers) of a mode (cu,cv,cw)
  Subroutine project_mode(T, mx, my, mz, cu, cv, cw)

    Type(ptabs),    Intent(In)    :: T
    Integer(Int32), Intent(In)    :: mx, my, mz
    Complex(Int64), Intent(InOut) :: cu, cv, cw
    Real   (Int64) :: kx, ky, kz, k2
    Complex(Int64) :: dot

    kx = T%km(1,mx);  ky = T%km(2,my);  kz = T%km(3,mz)
    k2 = kx*kx + ky*ky + kz*kz
    If ( k2 <= 0d0 ) Return
    dot = (kx*cu + ky*cv + kz*cw)/k2
    cu = cu - kx*dot
    cv = cv - ky*dot
    cw = cw - kz*dot

  End Subroutine project_mode

  Subroutine hit_setup

    Integer(Int32) :: nmax, mx, my, mz

    nmax = Floor(hit_kf)
    Call build_tables(nmax, TF)
    Call build_tables(0, TM)
    Allocate ( inband(-nmax:nmax,-nmax:nmax,-nmax:nmax) )
    inband = .False.
    Do mz = -nmax, nmax
       Do my = -nmax, nmax
          Do mx = -nmax, nmax
             inband(mx,my,mz) = ( mx*mx + my*my + mz*mz > 0 ) .And. &
                                ( Real(mx*mx + my*my + mz*mz, Int64) <= hit_kf*hit_kf )
          End Do
       End Do
    End Do
    Allocate ( fU(nx,nyg,nzg), fV(nxg,ny,nzg), fW(nxg,nyg,nz) )
    fU = 0d0;  fV = 0d0;  fW = 0d0
    Allocate ( RT%cx_f(nx,-nmax:nmax),  RT%sx_f(nx,-nmax:nmax),  RT%cx_c(nxg,-nmax:nmax), RT%sx_c(nxg,-nmax:nmax), &
               RT%cy_f(ny,-nmax:nmax),  RT%sy_f(ny,-nmax:nmax),  RT%cy_c(nyg,-nmax:nmax), RT%sy_c(nyg,-nmax:nmax), &
               RT%cz_f(nz,-nmax:nmax),  RT%sz_f(nz,-nmax:nmax),  RT%cz_c(nzg,-nmax:nmax), RT%sz_c(nzg,-nmax:nmax) )
    RT%cx_f = Real(TF%ex_f, Int64);  RT%sx_f = Aimag(TF%ex_f)
    RT%cx_c = Real(TF%ex_c, Int64);  RT%sx_c = Aimag(TF%ex_c)
    RT%cy_f = Real(TF%ey_f, Int64);  RT%sy_f = Aimag(TF%ey_f)
    RT%cy_c = Real(TF%ey_c, Int64);  RT%sy_c = Aimag(TF%ey_c)
    RT%cz_f = Real(TF%ez_f, Int64);  RT%sz_f = Aimag(TF%ez_f)
    RT%cz_c = Real(TF%ez_c, Int64);  RT%sz_c = Aimag(TF%ez_c)
    Allocate ( Wa_r(-nmax:nmax,nyg,nzg), Wa_i(-nmax:nmax,nyg,nzg), Wb_r(-nmax:nmax,-nmax:nmax,nzg), &
               Wb_i(-nmax:nmax,-nmax:nmax,nzg), Cc_r(-nmax:nmax,-nmax:nmax,-nmax:nmax), Cc_i(-nmax:nmax,-nmax:nmax,-nmax:nmax) )
    !$acc enter data copyin(fU,fV,fW,RT%cx_f,RT%sx_f,RT%cx_c,RT%sx_c,RT%cy_f,RT%sy_f,RT%cy_c,RT%sy_c) &
    !$acc& copyin(RT%cz_f,RT%sz_f,RT%cz_c,RT%sz_c) copyin(TF%mx_f,TF%mx_c,TF%my_f,TF%my_c,TF%mz_f,TF%mz_c) &
    !$acc& create(Wa_r,Wa_i,Wb_r,Wb_i,Cc_r,Cc_i)
    hit_ready = .True.
    If ( myid == 0 ) Write(*,'(A,I0,A,F5.2,A,ES10.3)') '   HIT forcing: ', Count(inband), &
         ' modes with |n| <= ', hit_kf, ', power eps_in = ', hit_eps

  End Subroutine hit_setup

  ! f_hat = eps_in/(2 E_band) u_hat on the solenoidal modes of the forced band, so the forcing injects exactly
  ! eps_in per unit mass; evaluated from the stage velocity, never touches the mean
  Subroutine hit_update_forcing(U_, V_, W_)

    Real(Int64), Intent(In) :: U_(nx,nyg,nzg), V_(nxg,ny,nzg), W_(nxg,nyg,nz)

    Real(Int64), Allocatable :: ur(:,:,:), ui(:,:,:), vr(:,:,:), vi(:,:,:), wr(:,:,:), wi(:,:,:)
    Complex(Int64) :: au, av, aw
    Integer(Int32) :: nm, mx, my, mz
    Real(Int64)    :: E, sc

    If ( .Not. hit_ready ) Call hit_setup
    nm = TF%nmax
    Allocate ( ur(-nm:nm,-nm:nm,-nm:nm), ui(-nm:nm,-nm:nm,-nm:nm), vr(-nm:nm,-nm:nm,-nm:nm), vi(-nm:nm,-nm:nm,-nm:nm), &
               wr(-nm:nm,-nm:nm,-nm:nm), wi(-nm:nm,-nm:nm,-nm:nm) )

    Call force_analyse(U_, nx, nyg, nzg, RT%cx_f, RT%sx_f, TF%mx_f, RT%cy_c, RT%sy_c, TF%my_c, &
                       RT%cz_c, RT%sz_c, TF%mz_c, nm, ur, ui)
    Call force_analyse(V_, nxg, ny, nzg, RT%cx_c, RT%sx_c, TF%mx_c, RT%cy_f, RT%sy_f, TF%my_f, &
                       RT%cz_c, RT%sz_c, TF%mz_c, nm, vr, vi)
    Call force_analyse(W_, nxg, nyg, nz, RT%cx_c, RT%sx_c, TF%mx_c, RT%cy_c, RT%sy_c, TF%my_c, &
                       RT%cz_f, RT%sz_f, TF%mz_f, nm, wr, wi)

    E = 0d0
    Do mz = -nm, nm
       Do my = -nm, nm
          Do mx = -nm, nm
             If ( inband(mx,my,mz) ) Then
                au = Cmplx(ur(mx,my,mz), ui(mx,my,mz), Int64)
                av = Cmplx(vr(mx,my,mz), vi(mx,my,mz), Int64)
                aw = Cmplx(wr(mx,my,mz), wi(mx,my,mz), Int64)
                Call project_mode(TF, mx, my, mz, au, av, aw)
                E = E + 0.5d0*( Abs(au)**2 + Abs(av)**2 + Abs(aw)**2 )
                ur(mx,my,mz) = Real(au, Int64);  ui(mx,my,mz) = Aimag(au)
                vr(mx,my,mz) = Real(av, Int64);  vi(mx,my,mz) = Aimag(av)
                wr(mx,my,mz) = Real(aw, Int64);  wi(mx,my,mz) = Aimag(aw)
             Else
                ur(mx,my,mz) = 0d0;  ui(mx,my,mz) = 0d0;  vr(mx,my,mz) = 0d0
                vi(mx,my,mz) = 0d0;  wr(mx,my,mz) = 0d0;  wi(mx,my,mz) = 0d0
             End If
          End Do
       End Do
    End Do
    sc = 0d0
    If ( E > 1d-300 ) sc = hit_eps/(2d0*E)
    ur = sc*ur;  ui = sc*ui;  vr = sc*vr;  vi = sc*vi;  wr = sc*wr;  wi = sc*wi

    Call force_synth(fU, nx, nyg, nzg, RT%cx_f, RT%sx_f, RT%cy_c, RT%sy_c, RT%cz_c, RT%sz_c, nm, ur, ui)
    Call force_synth(fV, nxg, ny, nzg, RT%cx_c, RT%sx_c, RT%cy_f, RT%sy_f, RT%cz_c, RT%sz_c, nm, vr, vi)
    Call force_synth(fW, nxg, nyg, nz, RT%cx_c, RT%sx_c, RT%cy_c, RT%sy_c, RT%cz_f, RT%sz_f, nm, wr, wi)
    Deallocate(ur, ui, vr, vi, wr, wi)

  End Subroutine hit_update_forcing

  ! ic_type=8: random-phase solenoidal field with a Passot-Pouquet spectrum E(k) ~ k^4 exp(-2 (k/kp)^2)
  ! (k in units of the fundamental wavenumber), rescaled to u_rms = Utarget per component. Every rank draws
  ! the same coefficients, so the field does not depend on the decomposition
  Subroutine hit_generate_ic(U_, V_, W_)

    Real(Int64), Intent(Out) :: U_(nx,nyg,nzg), V_(nxg,ny,nzg), W_(nxg,nyg,nz)

    Type(ptabs) :: T
    Complex(Int64), Allocatable :: cu(:,:,:), cv(:,:,:), cw(:,:,:)
    Integer(Int32) :: nmax, mx, my, mz, i, nmode
    Real(Int64)    :: n2, amp, E, Etot, r(6), g(6), norm, ang
    Complex(Int64) :: au, av, aw

    nmax = Min( Ceiling(2.5d0*ic_kp), (Min(nx_global, ny_global, nz_global) - 2)/3 )
    Call build_tables(nmax, T)
    Allocate ( cu(-nmax:nmax,-nmax:nmax,-nmax:nmax), cv(-nmax:nmax,-nmax:nmax,-nmax:nmax), &
               cw(-nmax:nmax,-nmax:nmax,-nmax:nmax) )
    cu = (0d0,0d0);  cv = (0d0,0d0);  cw = (0d0,0d0)

    rng_state = 1_Int64 + Mod(Int(ic_seed, Int64)*104729_Int64, 2147483646_Int64)

    nmode = 0
    Etot  = 0d0
    Do mz = -nmax, nmax
       Do my = -nmax, nmax
          Do mx = -nmax, nmax
             ! one of each +-n pair, the other is its conjugate
             If ( .Not. ( mx > 0 .Or. ( mx == 0 .And. my > 0 ) .Or. ( mx == 0 .And. my == 0 .And. mz > 0 ) ) ) Cycle
             n2 = Real(mx*mx + my*my + mz*mz, Int64)
             If ( n2 > Real(nmax*nmax, Int64) ) Cycle
             Do i = 1, 6
                r(i) = urand()
             End Do
             Do i = 1, 3
                norm = Sqrt(-2d0*Log(r(2*i-1)))
                ang  = 2d0*pi*r(2*i)
                g(2*i-1) = norm*Cos(ang)
                g(2*i)   = norm*Sin(ang)
             End Do
             au = Cmplx(g(1), g(2), Int64);  av = Cmplx(g(3), g(4), Int64);  aw = Cmplx(g(5), g(6), Int64)
             Call project_mode(T, mx, my, mz, au, av, aw)
             norm = Sqrt(Abs(au)**2 + Abs(av)**2 + Abs(aw)**2)
             If ( norm <= 0d0 ) Cycle
             ! spectral energy density per mode: E(k)/(4 pi k^2)
             E   = n2**2*Exp(-2d0*n2/ic_kp**2)/(4d0*pi*n2)
             amp = Sqrt(2d0*E)/norm
             cu(mx,my,mz) = amp*au;  cv(mx,my,mz) = amp*av;  cw(mx,my,mz) = amp*aw
             cu(-mx,-my,-mz) = Conjg(cu(mx,my,mz))
             cv(-mx,-my,-mz) = Conjg(cv(mx,my,mz))
             cw(-mx,-my,-mz) = Conjg(cw(mx,my,mz))
             Etot  = Etot + 2d0*E
             nmode = nmode + 1
          End Do
       End Do
    End Do
    ! TKE = 1.5 u_rms^2 = 0.5 sum |c|^2 over all modes
    amp = Sqrt(1.5d0*Utarget**2/Etot)
    cu = amp*cu;  cv = amp*cv;  cw = amp*cw

    Call synthesise(U_, nx,  nyg, nzg, T%ex_f, T%ey_c, T%ez_c, nmax, cu)
    Call synthesise(V_, nxg, ny,  nzg, T%ex_c, T%ey_f, T%ez_c, nmax, cv)
    Call synthesise(W_, nxg, nyg, nz,  T%ex_c, T%ey_c, T%ez_f, nmax, cw)

    If ( myid == 0 ) Write(*,'(A,I0,A,F6.2,A,I0,A,ES10.3)') '   IC = HIT random solenoidal field (ic_type=8): ', nmode, &
         ' mode pairs, kp = ', ic_kp, ', nmax = ', nmax, ', u_rms = ', Utarget
    Deallocate(cu, cv, cw)

  End Subroutine hit_generate_ic

  ! Lehmer generator (48271 mod 2^31-1) in (0,1): the same field with every compiler, unlike Random_number
  Real(Int64) Function urand()

    rng_state = Mod(48271_Int64*rng_state, 2147483647_Int64)
    urand = Real(rng_state, Int64)/2147483647d0

  End Function urand

  ! sum over the owned points of ((F(p+1)-F(p))/h)^2 along direction d; the masks are indexed by the lower
  ! point p along d
  Real(Int64) Function grad_sq(F, n1, n2, n3, d, h, m1, m2, m3)

    Integer(Int32), Intent(In) :: n1, n2, n3, d
    Real   (Int64), Intent(In) :: F(n1,n2,n3), h, m1(n1), m2(n2), m3(n3)
    Integer(Int32) :: i, j, k, e1, e2, e3
    Real   (Int64) :: s, g

    e1 = n1;  e2 = n2;  e3 = n3
    If ( d == 1 ) e1 = n1 - 1
    If ( d == 2 ) e2 = n2 - 1
    If ( d == 3 ) e3 = n3 - 1
    s = 0d0
    Do k = 1, e3
       Do j = 1, e2
          Do i = 1, e1
             If ( m1(i)*m2(j)*m3(k) == 0d0 ) Cycle
             If ( d == 1 ) g = F(i+1,j,k) - F(i,j,k)
             If ( d == 2 ) g = F(i,j+1,k) - F(i,j,k)
             If ( d == 3 ) g = F(i,j,k+1) - F(i,j,k)
             s = s + (g/h)**2
          End Do
       End Do
    End Do
    grad_sq = s

  End Function grad_sq

  ! mask padded with zeros to length n, or shifted down by one (m(p+1)) when sh is set
  Function mask_for(m, n, sh) Result(r)

    Real   (Int64), Intent(In) :: m(:)
    Integer(Int32), Intent(In) :: n
    Logical,        Intent(In) :: sh
    Real   (Int64) :: r(n)
    Integer(Int32) :: p, q

    r = 0d0
    Do p = 1, n
       q = p
       If ( sh ) q = p + 1
       If ( q <= Size(m) ) r(p) = m(q)
    End Do

  End Function mask_for

  ! TKE, dissipation eps = nu <du_i/dx_j du_i/dx_j> and the derived scales, appended to hit_stats.dat
  Subroutine hit_write_stats

    Real(Int64) :: s(15), g(15), d(3), u2, tke, eps, urms, lam, rel, eta, Lint
    Integer(Int32) :: i

    If ( .Not. hit_ready ) Then
       Call build_tables(0, TM)
       hit_ready = .True.
    End If
    d(1) = x_global(2) - x_global(1)
    d(2) = y_global(2) - y_global(1)
    d(3) = z_global(2) - z_global(1)
    !$acc update host(U,V,W)

    s = 0d0
    s(1) = Sum_masked(U, nx,  nyg, nzg, TM%mx_f, TM%my_c, TM%mz_c, 1)
    s(2) = Sum_masked(V, nxg, ny,  nzg, TM%mx_c, TM%my_f, TM%mz_c, 1)
    s(3) = Sum_masked(W, nxg, nyg, nz,  TM%mx_c, TM%my_c, TM%mz_f, 1)
    s(4) = Sum_masked(U, nx,  nyg, nzg, TM%mx_f, TM%my_c, TM%mz_c, 2)
    s(5) = Sum_masked(V, nxg, ny,  nzg, TM%mx_c, TM%my_f, TM%mz_c, 2)
    s(6) = Sum_masked(W, nxg, nyg, nz,  TM%mx_c, TM%my_c, TM%mz_f, 2)
    s(7)  = grad_sq(U, nx,  nyg, nzg, 1, d(1), mask_for(TM%mx_c,nx,.True.),  mask_for(TM%my_c,nyg,.False.), &
                    mask_for(TM%mz_c,nzg,.False.))
    s(8)  = grad_sq(V, nxg, ny,  nzg, 2, d(2), mask_for(TM%mx_c,nxg,.False.), mask_for(TM%my_c,ny,.True.),  &
                    mask_for(TM%mz_c,nzg,.False.))
    s(9)  = grad_sq(W, nxg, nyg, nz,  3, d(3), mask_for(TM%mx_c,nxg,.False.), mask_for(TM%my_c,nyg,.False.), &
                    mask_for(TM%mz_c,nz,.True.))
    s(10) = grad_sq(U, nx,  nyg, nzg, 2, d(2), mask_for(TM%mx_f,nx,.False.),  mask_for(TM%my_f,nyg,.False.), &
                    mask_for(TM%mz_c,nzg,.False.))
    s(11) = grad_sq(V, nxg, ny,  nzg, 1, d(1), mask_for(TM%mx_f,nxg,.False.), mask_for(TM%my_f,ny,.False.),  &
                    mask_for(TM%mz_c,nzg,.False.))
    s(12) = grad_sq(U, nx,  nyg, nzg, 3, d(3), mask_for(TM%mx_f,nx,.False.),  mask_for(TM%my_c,nyg,.False.), &
                    mask_for(TM%mz_f,nzg,.False.))
    s(13) = grad_sq(W, nxg, nyg, nz,  1, d(1), mask_for(TM%mx_f,nxg,.False.), mask_for(TM%my_c,nyg,.False.), &
                    mask_for(TM%mz_f,nz,.False.))
    s(14) = grad_sq(V, nxg, ny,  nzg, 3, d(3), mask_for(TM%mx_c,nxg,.False.), mask_for(TM%my_f,ny,.False.),  &
                    mask_for(TM%mz_f,nzg,.False.))
    s(15) = grad_sq(W, nxg, nyg, nz,  2, d(2), mask_for(TM%mx_c,nxg,.False.), mask_for(TM%my_f,nyg,.False.), &
                    mask_for(TM%mz_f,nz,.False.))
    Call MPI_Allreduce(s, g, 15, MPI_real8, MPI_SUM, MPI_COMM_WORLD, ierr)
    g = g/n_tot

    If ( myid /= 0 ) Return
    u2   = 0d0
    Do i = 1, 3
       u2 = u2 + g(3+i) - g(i)**2
    End Do
    tke  = 0.5d0*u2
    eps  = nu*Sum(g(7:15))
    urms = Sqrt(2d0*tke/3d0)
    lam  = Sqrt(15d0*nu*urms**2/eps)
    rel  = urms*lam/nu
    eta  = (nu**3/eps)**0.25d0
    Lint = urms**3/eps
    If ( .Not. stats_started ) Then
       stats_started = .True.
       If ( restart == 0 ) Then
          Open(unit=88, file='hit_stats.dat', status='replace')
          Write(88,'(A)') '# step t TKE eps u_rms lambda Re_lambda eta eta/dx L=u^3/eps T_e=L/u_rms'
          Close(88)
       End If
    End If
    Open(unit=88, file='hit_stats.dat', status='unknown', position='append')
    Write(88,'(I10,10ES15.7)') istep + nstep_init, t, tke, eps, urms, lam, rel, eta, eta/d(1), Lint, Lint/urms
    Close(88)

  End Subroutine hit_write_stats

  Real(Int64) Function Sum_masked(F, n1, n2, n3, m1, m2, m3, pw)

    Integer(Int32), Intent(In) :: n1, n2, n3, pw
    Real   (Int64), Intent(In) :: F(n1,n2,n3), m1(n1), m2(n2), m3(n3)
    Integer(Int32) :: i, j, k
    Real   (Int64) :: s

    s = 0d0
    Do k = 1, n3
       Do j = 1, n2
          Do i = 1, n1
             s = s + m1(i)*m2(j)*m3(k)*F(i,j,k)**pw
          End Do
       End Do
    End Do
    Sum_masked = s

  End Function Sum_masked

End Module hit
