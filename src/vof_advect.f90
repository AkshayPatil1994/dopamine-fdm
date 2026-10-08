!> Geometric PLIC advection of the liquid volume fraction on the staggered Cartesian grid: Weymouth & Yue (2010) directional sweeps,
!  each face flux being the exact liquid volume in the slab |u| dt upstream of the face.
!
!  Array conventions follow the solver: cell-centred C lives on a padded array Cp(0:n1+1,0:n2+1,0:n3+1) whose interior cells are
!  2..n-1 (1 and n are ghost cells, 0 and n+1 the second layer needed for the donor cell's normal); face velocities U(1:n1-1,:,:),
!  V(:,1:n2-1,:), W(:,:,1:n3-1) sit on the upper face of cell i (U(i) between cells i and i+1). Cell widths are h1(1:n1), h2(1:n2),
!  h3(1:n3), so stretched directions are allowed. The velocity must be discretely divergence-free for the update to stay bounded
!  and conservative; no module here depends on MPI, the solver globals or the I/O.
Module vof_advect

  Use iso_fortran_env, Only : Int32, Int64
  Use vof_plic
  Use vof_normals

  Implicit None

  Real(Int64), Allocatable, Dimension(:,:,:) :: vof_mx, vof_my, vof_mz, vof_al, vof_flux, vof_cc
  ! face-flux model: 0 = geometric PLIC, 1 = THINC (tanh profile, diffuse interface over ~2-3 cells; beta sets the sharpness)
  Integer(Int32) :: vof_flux_scheme = 0
  Real(Int64)    :: vof_thinc_beta = 2d0

  Abstract Interface
     Subroutine fill_pad_iface(Cp, n1, n2, n3)
       Import :: Int32, Int64
       Integer(Int32), Intent(In)    :: n1, n2, n3
       Real(Int64),    Intent(InOut) :: Cp(0:n1+1,0:n2+1,0:n3+1)
     End Subroutine fill_pad_iface
  End Interface

Contains

  Subroutine vof_advect_init(n1, n2, n3)

    Integer(Int32), Intent(In) :: n1, n2, n3

    If ( Allocated(vof_mx) ) Then
       If ( Size(vof_mx,1) == n1 .And. Size(vof_mx,2) == n2 .And. Size(vof_mx,3) == n3 ) Return
       !$acc exit data delete(vof_mx,vof_my,vof_mz,vof_al,vof_flux,vof_cc)
       Deallocate( vof_mx, vof_my, vof_mz, vof_al, vof_flux, vof_cc )
    End If
    Allocate( vof_mx(n1,n2,n3), vof_my(n1,n2,n3), vof_mz(n1,n2,n3), vof_al(n1,n2,n3), vof_flux(n1,n2,n3), vof_cc(n1,n2,n3) )
    vof_mx = 0d0;  vof_my = 0d0;  vof_mz = 0d0;  vof_al = 0d0;  vof_flux = 0d0;  vof_cc = 0d0
    !$acc enter data copyin(vof_mx,vof_my,vof_mz,vof_al,vof_flux,vof_cc)

  End Subroutine vof_advect_init


  !> Plane normal and constant in every cell holding an interface (0 < C < 1), from the padded C
  Subroutine vof_reconstruct(Cp, n1, n2, n3, scheme)

    Integer(Int32), Intent(In) :: n1, n2, n3, scheme
    Real(Int64),    Intent(In) :: Cp(0:n1+1,0:n2+1,0:n3+1)

    Integer(Int32) :: i, j, k
    Real(Int64) :: st(-1:1,-1:1,-1:1), m(3), c

    !$acc parallel loop collapse(3) present(Cp,vof_mx,vof_my,vof_mz,vof_al) private(st,m,c)
    Do k = 1, n3
       Do j = 1, n2
          Do i = 1, n1
             c = Cp(i,j,k)
             If ( c > vof_eps .And. c < 1d0 - vof_eps ) Then
                st = Cp(i-1:i+1,j-1:j+1,k-1:k+1)
                If ( scheme == VOF_SCHEME_CC ) Then
                   Call vof_normal_cc(st, m)
                Else
                   Call vof_normal_youngs(st, m)
                End If
                vof_mx(i,j,k) = m(1)
                vof_my(i,j,k) = m(2)
                vof_mz(i,j,k) = m(3)
                vof_al(i,j,k) = plic_alpha(c, m(1), m(2), m(3))
             Else
                vof_mx(i,j,k) = 0d0
                vof_my(i,j,k) = 0d0
                vof_mz(i,j,k) = 0d0
                vof_al(i,j,k) = 0d0
             End If
          End Do
       End Do
    End Do
    !$acc end parallel loop

  End Subroutine vof_reconstruct


  Pure Function lncosh(x) Result(r)
    !$acc routine seq
    Real(Int64), Intent(In) :: x
    Real(Int64) :: r, a
    a = Abs(x)
    r = a + Log(1d0 + Exp(-2d0*a)) - 0.6931471805599453d0
  End Function lncosh


  !> One directional sweep (dir = 1,2,3): face volume fluxes from the reconstructed planes, then the conservative update with the
  !  c-tilde * div(u) correction. co_max = max donor-cell Courant number seen, clip_loss = liquid volume removed by clamping C to [0,1]
  Subroutine vof_sweep(dir, Cp, n1, n2, n3, uf, h1, h2, h3, dt, co_max, clip_loss)

    Integer(Int32), Intent(In)    :: dir, n1, n2, n3
    Real(Int64),    Intent(InOut) :: Cp(0:n1+1,0:n2+1,0:n3+1)
    Real(Int64),    Intent(In)    :: uf(:,:,:)
    Real(Int64),    Intent(In)    :: h1(n1), h2(n2), h3(n3), dt
    Real(Int64),    Intent(InOut) :: co_max, clip_loss

    Integer(Int32) :: i, j, k, di, dj, dk, id, jd, kd, i0, i1, j0, j1, k0, k1
    Real(Int64) :: u, hd, vol, s, c, frac, fl, xlo(3), w(3), cn, cnc, vc, area, comax, closs, mdir, mn, gam, beta, x0, xa, xb

    di = 0;  dj = 0;  dk = 0
    i0 = 2;  i1 = n1-1;  j0 = 2;  j1 = n2-1;  k0 = 2;  k1 = n3-1
    If ( dir == 1 ) Then
       di = 1;  i0 = 1
    Else If ( dir == 2 ) Then
       dj = 1;  j0 = 1
    Else
       dk = 1;  k0 = 1
    End If

    comax = co_max
    !$acc parallel loop collapse(3) present(Cp,uf,h1,h2,h3,vof_mx,vof_my,vof_mz,vof_al,vof_flux) &
    !$acc& private(xlo,w) reduction(max:comax)
    Do k = k0, k1
       Do j = j0, j1
          Do i = i0, i1
             u = uf(i,j,k)
             If ( u >= 0d0 ) Then
                id = i;  jd = j;  kd = k
             Else
                id = i+di;  jd = j+dj;  kd = k+dk
             End If
             If ( dir == 1 ) Then
                hd = h1(id)
             Else If ( dir == 2 ) Then
                hd = h2(jd)
             Else
                hd = h3(kd)
             End If
             vol = h1(id)*h2(jd)*h3(kd)
             s = Abs(u)*dt/hd
             comax = Max(comax, s)
             s = Min(s, 1d0)
             c = Cp(id,jd,kd)
             If ( c <= vof_eps ) Then
                frac = 0d0
             Else If ( c >= 1d0 - vof_eps ) Then
                frac = 1d0
             Else If ( vof_flux_scheme == 1 ) Then
                If ( dir == 1 ) Then
                   mdir = vof_mx(id,jd,kd)
                Else If ( dir == 2 ) Then
                   mdir = vof_my(id,jd,kd)
                Else
                   mdir = vof_mz(id,jd,kd)
                End If
                mn = Sqrt( vof_mx(id,jd,kd)**2 + vof_my(id,jd,kd)**2 + vof_mz(id,jd,kd)**2 )
                ! the profile steepens only across the interface: beta scales with the normal's component along the sweep
                beta = vof_thinc_beta*Abs(mdir)/Max(mn, 1d-300)
                If ( beta < 1d-3 .Or. s < 1d-14 ) Then
                   frac = c
                Else
                   gam = Merge(1d0, -1d0, mdir < 0d0)
                   x0 = Atanh( -Tanh(beta*gam*(c - 0.5d0))/Tanh(0.5d0*beta) )/beta
                   If ( u >= 0d0 ) Then
                      xa = 0.5d0 - s;  xb = 0.5d0
                   Else
                      xa = -0.5d0;  xb = -0.5d0 + s
                   End If
                   frac = 0.5d0 + gam/(2d0*beta*s)*( lncosh(beta*(xb - x0)) - lncosh(beta*(xa - x0)) )
                   frac = Min(1d0, Max(0d0, frac))
                End If
             Else
                xlo = -0.5d0
                w   = 1d0
                w(dir) = s
                If ( u >= 0d0 ) xlo(dir) = 0.5d0 - s
                frac = plic_subbox_fraction( vof_mx(id,jd,kd), vof_my(id,jd,kd), vof_mz(id,jd,kd), vof_al(id,jd,kd), &
                                             xlo(1), xlo(2), xlo(3), w(1), w(2), w(3) )
             End If
             vof_flux(i,j,k) = Sign(1d0, u)*frac*s*vol
          End Do
       End Do
    End Do
    !$acc end parallel loop
    co_max = comax

    closs = clip_loss
    !$acc parallel loop collapse(3) present(Cp,uf,h1,h2,h3,vof_flux,vof_cc) reduction(+:closs)
    Do k = 2, n3-1
       Do j = 2, n2-1
          Do i = 2, n1-1
             vc = h1(i)*h2(j)*h3(k)
             If ( dir == 1 ) Then
                area = vc/h1(i)
             Else If ( dir == 2 ) Then
                area = vc/h2(j)
             Else
                area = vc/h3(k)
             End If
             fl = vof_flux(i-di,j-dj,k-dk) - vof_flux(i,j,k)
             cn = Cp(i,j,k) + ( fl + vof_cc(i,j,k)*dt*area*( uf(i,j,k) - uf(i-di,j-dj,k-dk) ) )/vc
             cnc = Min(1d0, Max(0d0, cn))
             closs = closs + (cn - cnc)*vc
             Cp(i,j,k) = cnc
          End Do
       End Do
    End Do
    !$acc end parallel loop
    clip_loss = closs

  End Subroutine vof_sweep


  !> Full advection step: normals/planes are rebuilt before each sweep; the sweep order alternates with istep to cancel the
  !  directional-splitting bias. c-tilde (the c > 1/2 indicator) is frozen for the whole step so the three corrections sum to
  !  c-tilde * div(u) = 0.
  Subroutine vof_advect_step(Cp, n1, n2, n3, U, V, W, h1, h2, h3, dt, istep, scheme, fill_pad, co_max, clip_loss)

    Integer(Int32), Intent(In)    :: n1, n2, n3, istep, scheme
    Real(Int64),    Intent(InOut) :: Cp(0:n1+1,0:n2+1,0:n3+1)
    Real(Int64),    Intent(In)    :: U(:,:,:), V(:,:,:), W(:,:,:)
    Real(Int64),    Intent(In)    :: h1(n1), h2(n2), h3(n3), dt
    Procedure(fill_pad_iface)     :: fill_pad
    Real(Int64),    Intent(Out)   :: co_max, clip_loss

    Integer(Int32) :: isw, d, order(3)
    Integer(Int32) :: i, j, k

    co_max = 0d0
    clip_loss = 0d0
    Call vof_advect_init(n1, n2, n3)

    !$acc parallel loop collapse(3) present(Cp,vof_cc)
    Do k = 1, n3
       Do j = 1, n2
          Do i = 1, n1
             vof_cc(i,j,k) = Merge(1d0, 0d0, Cp(i,j,k) > 0.5d0)
          End Do
       End Do
    End Do
    !$acc end parallel loop

    If ( Mod(istep, 2) == 0 ) Then
       order = (/ 1, 2, 3 /)
    Else
       order = (/ 3, 2, 1 /)
    End If

    Do isw = 1, 3
       d = order(isw)
       Call fill_pad(Cp, n1, n2, n3)
       Call vof_reconstruct(Cp, n1, n2, n3, scheme)
       If ( d == 1 ) Then
          Call vof_sweep(1, Cp, n1, n2, n3, U, h1, h2, h3, dt, co_max, clip_loss)
       Else If ( d == 2 ) Then
          Call vof_sweep(2, Cp, n1, n2, n3, V, h1, h2, h3, dt, co_max, clip_loss)
       Else
          Call vof_sweep(3, Cp, n1, n2, n3, W, h1, h2, h3, dt, co_max, clip_loss)
       End If
    End Do
    Call fill_pad(Cp, n1, n2, n3)

  End Subroutine vof_advect_step


  !> Local-rank diagnostics (reduce across ranks in the caller): liquid volume, min/max C, interface-cell count, volume of the
  !  interface cells' mixed liquid. Ghost cells are excluded.
  Subroutine vof_local_stats(Cp, n1, n2, n3, h1, h2, h3, vliq, cmin, cmax, nint)

    Integer(Int32), Intent(In)  :: n1, n2, n3
    Real(Int64),    Intent(In)  :: Cp(0:n1+1,0:n2+1,0:n3+1), h1(n1), h2(n2), h3(n3)
    Real(Int64),    Intent(Out) :: vliq, cmin, cmax
    Integer(Int64), Intent(Out) :: nint

    Integer(Int32) :: i, j, k
    Real(Int64) :: c, v, cmn, cmx
    Integer(Int64) :: ni

    v = 0d0;  cmn = 1d300;  cmx = -1d300;  ni = 0
    !$acc parallel loop collapse(3) present(Cp,h1,h2,h3) reduction(+:v,ni) reduction(min:cmn) reduction(max:cmx)
    Do k = 2, n3-1
       Do j = 2, n2-1
          Do i = 2, n1-1
             c = Cp(i,j,k)
             v = v + c*h1(i)*h2(j)*h3(k)
             cmn = Min(cmn, c)
             cmx = Max(cmx, c)
             If ( c > vof_eps .And. c < 1d0 - vof_eps ) ni = ni + 1
          End Do
       End Do
    End Do
    !$acc end parallel loop
    vliq = v;  cmin = cmn;  cmax = cmx;  nint = ni

  End Subroutine vof_local_stats

End Module vof_advect
