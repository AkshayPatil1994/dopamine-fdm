!> Mask-aware geometric multigrid preconditioner for the masked PCG of vof_pressure (pcg_precond = 1).
!
!  The operator is the symmetric conductance form S = -V A of A = div(beta grad) (V the cell volumes), so that z = A^-1 r is
!  computed as z = -S^-1 (V r) and the preconditioner is symmetric negative definite like A. Coarse cells are aggregates of
!  2x2x2 (or fewer, see mg_coarsen) rank-local fine cells, the transfers are piecewise constant (restriction = sum = transpose
!  of the prolongation) and the coarse conductances are the sums of the fine conductances crossing a coarse face, scaled by
!  (fine centre distance)/(coarse centre distance) in the coarsened direction. Closed faces therefore stay closed and a coarse
!  face is open if any fine face under it is open; solid cells and sealed pockets give empty or singular rows that the
!  smoother leaves alone. The smoother is a Chebyshev polynomial of the y-line-Jacobi iteration (exact tridiagonal solves along
!  the stretched y direction, which carries the grid anisotropy), the same polynomial before and after the coarse correction,
!  so the V-cycle from a zero guess is a fixed symmetric operator.
!
!  Level arrays are flat 1D storage with one offset per level, indexed (0:ni+1,0:nj+1,0:nk+1) through explicit-shape dummies:
!  cells 1:n are the distinct unknowns (the periodic duplicate cell of the last rank is not one), 0 and n+1 the ghost layer
!  (neighbour rank, periodic partner, Neumann copy or the Dirichlet reflection of the x outlet). cx(i) joins cells i and i+1.
Module vof_mg

  Use iso_fortran_env, Only : Int32, Int64
  Use global
  Use mpi
  Use decomp, Only : x_halo_neighbors, z_halo_neighbors, x_periodic_partner, z_periodic_partner

  Implicit None

  Integer(Int32), Parameter :: mg_maxlev = 24
  Integer(Int32) :: mg_nlev = 0, mg_deg = 3, mg_deg_coarse = 8
  Real(Int64) :: mg_ratio = 6d0, mg_ratio_coarse = 60d0, mg_cscale = 1d0, mg_eps = 1d-4
  Logical :: mg_coarsen_y = .True.
  Integer(Int32) :: mg_ni(mg_maxlev), mg_nj(mg_maxlev), mg_nk(mg_maxlev)
  Integer(Int32) :: mg_sx(mg_maxlev), mg_sy(mg_maxlev), mg_sz(mg_maxlev)   ! coarsening stride from level l to l+1
  Integer(Int32) :: mg_off(mg_maxlev+1)
  Real(Int64), Allocatable, Dimension(:) :: mg_cx, mg_cy, mg_cz, mg_dd, mg_x, mg_b, mg_t, mg_d, mg_m, mg_iw
  Real(Int64), Allocatable, Dimension(:,:) :: mg_fx, mg_fy, mg_fz   ! centre-distance ratios of the coarse faces, level l -> l+1
  Real(Int64), Allocatable, Dimension(:) :: mg_hy1, mg_hz1, mg_dcy1, mg_dcz1   ! level-1 cell sizes and centre distances
  Real(Int64), Allocatable, Dimension(:,:) :: mg_sbx, mg_rbx, mg_sbz, mg_rbz
  Integer(Int32) :: mg_xup, mg_xdn, mg_zup, mg_zdn
  Logical :: mg_xself, mg_zself
  Real(Int64) :: mg_xlo_gs, mg_xhi_gs   ! ghost = gs * boundary cell where there is no neighbour
  Logical :: mg_ready = .False.

Contains


  !> Hierarchy, neighbour topology and level-1 geometry; ihi, jhi, khi are the last unknown fine cells
  Subroutine mg_init(ihi, jhi, khi)

    Integer(Int32), Intent(In) :: ihi, jhi, khi

    Integer(Int32) :: l, j, k, nmin(2), nmin_g(2), ntot, nmax, partner
    Integer(Int32) :: ni, nj, nk, nic, njc, nkc, sx, sy, sz
    Logical :: is_first, is_last, cx_, cy_, cz_
    Real(Int64), Allocatable :: hx(:), hy(:), hz(:), dcx(:), dcy(:), dcz(:)
    Real(Int64), Allocatable :: hxn(:), hyn(:), hzn(:), dcxn(:), dcyn(:), dczn(:)

    If ( y_bc_type == 0 ) Then
       If ( myid == 0 ) Write(*,'(A)') ' ERROR: pcg_precond = 1 needs y walls (y_bc_type = 1)'
       Call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    End If

    Call x_halo_neighbors(mg_xup, mg_xdn)
    Call z_halo_neighbors(mg_zup, mg_zdn)
    mg_xself = .False.;  mg_zself = .False.
    If ( x_bc_type == 0 ) Then
       Call x_periodic_partner(is_first, is_last, partner)
       If ( is_first ) mg_xdn = partner
       If ( is_last ) mg_xup = partner
       mg_xself = ( is_first .And. is_last )
    End If
    If ( z_bc_type == 0 ) Then
       Call z_periodic_partner(is_first, is_last, partner)
       If ( is_first ) mg_zdn = partner
       If ( is_last ) mg_zup = partner
       mg_zself = ( is_first .And. is_last )
    End If
    mg_xlo_gs = 1d0;  mg_xhi_gs = 1d0
    If ( x_bc_type == 1 ) mg_xhi_gs = -1d0

    ni = ihi - 1;  nj = jhi - 1;  nk = khi - 1
    nmax = Max(ni, nj, nk) + 1
    Allocate( mg_fx(0:nmax,mg_maxlev), mg_fy(0:nmax,mg_maxlev), mg_fz(0:nmax,mg_maxlev) )
    mg_fx = 1d0;  mg_fy = 1d0;  mg_fz = 1d0
    Allocate( hx(0:ni+1), hy(0:nj+1), hz(0:nk+1), dcx(0:ni), dcy(0:nj), dcz(0:nk) )
    hx = dx;  dcx = dx
    Do j = 1, nj
       hy(j) = y(j+1) - y(j)
    End Do
    hy(0) = hy(1);  hy(nj+1) = hy(nj)
    Do j = 0, nj
       dcy(j) = yg(j+2) - yg(j+1)
    End Do
    Do k = 1, nk
       hz(k) = z(k+1) - z(k)
    End Do
    hz(0) = hz(1);  hz(nk+1) = hz(nk)
    Call mg_ghost_size(hz, nk, mg_zup, mg_zdn, mg_zself, 401)
    Do k = 0, nk
       dcz(k) = zg(k+2) - zg(k+1)
    End Do
    Allocate( mg_hy1(0:nj+1), mg_hz1(0:nk+1) )
    mg_hy1 = hy;  mg_hz1 = hz

    mg_ni(1) = ni;  mg_nj(1) = nj;  mg_nk(1) = nk
    mg_off(1) = 0
    l = 1
    Do
       ! a direction is coarsened if every rank holds at least two of its cells
       nmin = (/ mg_ni(l), mg_nk(l) /)
       Call MPI_Allreduce(nmin, nmin_g, 2, MPI_integer, MPI_MIN, MPI_COMM_WORLD, ierr)
       cx_ = ( nmin_g(1) >= 2 );  cz_ = ( nmin_g(2) >= 2 );  cy_ = ( mg_coarsen_y .And. mg_nj(l) >= 2 )
       If ( .Not. ( cx_ .Or. cy_ .Or. cz_ ) .Or. l == mg_maxlev ) Exit
       sx = Merge(2, 1, cx_);  sy = Merge(2, 1, cy_);  sz = Merge(2, 1, cz_)
       mg_sx(l) = sx;  mg_sy(l) = sy;  mg_sz(l) = sz
       ni = mg_ni(l);  nj = mg_nj(l);  nk = mg_nk(l)
       nic = (ni + sx - 1)/sx;  njc = (nj + sy - 1)/sy;  nkc = (nk + sz - 1)/sz
       mg_ni(l+1) = nic;  mg_nj(l+1) = njc;  mg_nk(l+1) = nkc
       mg_off(l+1) = mg_off(l) + (ni+2)*(nj+2)*(nk+2)
       Allocate( hxn(0:nic+1), hyn(0:njc+1), hzn(0:nkc+1), dcxn(0:nic), dcyn(0:njc), dczn(0:nkc) )
       Call mg_coarse_geometry(hx, dcx, ni, sx, nic, mg_xup, mg_xdn, mg_xself, 402, hxn, dcxn, mg_fx(0,l))
       Call mg_coarse_geometry(hy, dcy, nj, sy, njc, MPI_PROC_NULL, MPI_PROC_NULL, .False., 403, hyn, dcyn, mg_fy(0,l))
       Call mg_coarse_geometry(hz, dcz, nk, sz, nkc, mg_zup, mg_zdn, mg_zself, 404, hzn, dczn, mg_fz(0,l))
       Call move_alloc(hxn, hx);  Call move_alloc(hyn, hy);  Call move_alloc(hzn, hz)
       Call move_alloc(dcxn, dcx);  Call move_alloc(dcyn, dcy);  Call move_alloc(dczn, dcz)
       l = l + 1
    End Do
    mg_nlev = l
    mg_off(l+1) = mg_off(l) + (mg_ni(l)+2)*(mg_nj(l)+2)*(mg_nk(l)+2)
    ntot = mg_off(mg_nlev+1)
    Deallocate( hx, hy, hz, dcx, dcy, dcz )

    Allocate( mg_cx(ntot), mg_cy(ntot), mg_cz(ntot), mg_dd(ntot), mg_x(ntot), mg_b(ntot), mg_t(ntot), mg_d(ntot), &
              mg_m(ntot), mg_iw(ntot) )
    mg_cx = 0d0;  mg_cy = 0d0;  mg_cz = 0d0;  mg_dd = 0d0;  mg_x = 0d0;  mg_b = 0d0;  mg_t = 0d0;  mg_d = 0d0
    mg_m = 0d0;  mg_iw = 0d0
    Allocate( mg_sbx((mg_nj(1)+2)*(mg_nk(1)+2),2), mg_rbx((mg_nj(1)+2)*(mg_nk(1)+2),2) )
    Allocate( mg_sbz((mg_ni(1)+2)*(mg_nj(1)+2),2), mg_rbz((mg_ni(1)+2)*(mg_nj(1)+2),2) )
    !$acc enter data copyin(mg_cx,mg_cy,mg_cz,mg_dd,mg_x,mg_b,mg_t,mg_d,mg_m,mg_iw)
    !$acc enter data create(mg_sbx,mg_rbx,mg_sbz,mg_rbz)
    !$acc enter data copyin(mg_hy1,mg_hz1)
    mg_ready = .True.

    If ( myid == 0 ) Then
       Write(*,'(A,I0,A)') '   GMG preconditioner: ', mg_nlev, ' levels (rank-0 sizes ni,nj,nk):'
       Do l = 1, mg_nlev
          Write(*,'(A,I3,3I6)') '     level', l, mg_ni(l), mg_nj(l), mg_nk(l)
       End Do
    End If

  End Subroutine mg_init




  !> Ghost sizes of a 1D size array h(1:n): the neighbour's end cell, or the own cell where there is no neighbour
  Subroutine mg_ghost_size(h, n, up, dn, self, tag)

    Integer(Int32), Intent(In) :: n, up, dn, tag
    Logical, Intent(In) :: self
    Real(Int64), Intent(InOut) :: h(0:n+1)
    Real(Int64) :: sl, sh

    h(0) = h(1);  h(n+1) = h(n)
    If ( self ) Then
       h(0) = h(n);  h(n+1) = h(1)
       Return
    End If
    sl = h(n);  sh = h(1)
    Call MPI_Sendrecv(sl, 1, MPI_real8, up, tag, h(0), 1, MPI_real8, dn, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    Call MPI_Sendrecv(sh, 1, MPI_real8, dn, tag+10, h(n+1), 1, MPI_real8, up, tag+10, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

  End Subroutine mg_ghost_size


  !> Sizes and centre distances of the coarse cells of one direction, and the ratio fr(I) of the centre distance of the fine face
  !  under each coarse face to that of the coarse face (conductance of a coarsened direction scales with it)
  Subroutine mg_coarse_geometry(h, dc, n, s, nc, up, dn, self, tag, hc, dcc, fr)

    Integer(Int32), Intent(In) :: n, s, nc, up, dn, tag
    Logical, Intent(In) :: self
    Real(Int64), Intent(In) :: h(0:n+1), dc(0:n)
    Real(Int64), Intent(Out) :: hc(0:nc+1), dcc(0:nc), fr(0:nc)
    Integer(Int32) :: ic, ii

    hc = 0d0
    Do ic = 1, nc
       Do ii = s*ic - (s-1), Min(s*ic, n)
          hc(ic) = hc(ic) + h(ii)
       End Do
    End Do
    Call mg_ghost_size(hc, nc, up, dn, self, tag)
    Do ic = 0, nc
       dcc(ic) = 0.5d0*( hc(ic) + hc(ic+1) )
       fr(ic) = dc(Min(s*ic, n))/dcc(ic)
    End Do

  End Subroutine mg_coarse_geometry


  !> Conductances of all levels from the fine face coefficients (beta on the faces of U, V, W; zero on closed faces), the
  !  diagonals and the tridiagonal factors of the line smoother; call whenever the coefficients change
  Subroutine mg_set_coef(bu, bv, bw)

    Real(Int64), Intent(In) :: bu(nx,nyg,nzg), bv(nxg,ny,nzg), bw(nxg,nyg,nz)
    Integer(Int32) :: l, o, oc

    Call mg_fine_coef(mg_ni(1), mg_nj(1), mg_nk(1), bu, bv, bw, mg_cx, mg_cy, mg_cz)
    Do l = 1, mg_nlev-1
       o = mg_off(l) + 1;  oc = mg_off(l+1) + 1
       Call mg_coarsen_coef(mg_ni(l), mg_nj(l), mg_nk(l), mg_cx(o), mg_cy(o), mg_cz(o), &
                            mg_ni(l+1), mg_nj(l+1), mg_nk(l+1), mg_sx(l), mg_sy(l), mg_sz(l), &
                            mg_fx(0,l), mg_fy(0,l), mg_fz(0,l), mg_cx(oc), mg_cy(oc), mg_cz(oc))
    End Do
    Do l = 1, mg_nlev
       o = mg_off(l) + 1
       Call mg_factor(mg_ni(l), mg_nj(l), mg_nk(l), mg_cx(o), mg_cy(o), mg_cz(o), mg_dd(o), mg_m(o), mg_iw(o))
    End Do
    !$acc update device(mg_cx,mg_cy,mg_cz,mg_dd,mg_m,mg_iw)

  End Subroutine mg_set_coef


  !> Level-1 conductances beta*area/distance on the faces between unknown cells; physical Neumann faces carry nothing
  Subroutine mg_fine_coef(ni, nj, nk, bu, bv, bw, cx, cy, cz)

    Integer(Int32), Intent(In) :: ni, nj, nk
    Real(Int64), Intent(In) :: bu(nx,nyg,nzg), bv(nxg,ny,nzg), bw(nxg,nyg,nz)
    Real(Int64), Intent(Out) :: cx(0:ni+1,0:nj+1,0:nk+1), cy(0:ni+1,0:nj+1,0:nk+1), cz(0:ni+1,0:nj+1,0:nk+1)
    Integer(Int32) :: i, j, k, partner
    Logical :: is_first, is_last

    cx = 0d0;  cy = 0d0;  cz = 0d0
    Do k = 1, nk
       Do j = 1, nj
          Do i = 0, ni
             cx(i,j,k) = bu(i+1,j+1,k+1)*mg_hy1(j)*mg_hz1(k)/dx
          End Do
       End Do
    End Do
    Do k = 1, nk
       Do j = 1, nj-1
          Do i = 1, ni
             cy(i,j,k) = bv(i+1,j+1,k+1)*dx*mg_hz1(k)/( yg(j+2) - yg(j+1) )
          End Do
       End Do
    End Do
    Do k = 0, nk
       Do j = 1, nj
          Do i = 1, ni
             cz(i,j,k) = bw(i+1,j+1,k+1)*dx*mg_hy1(j)/( zg(k+2) - zg(k+1) )
          End Do
       End Do
    End Do
    If ( x_bc_type == 1 ) Then
       Call x_periodic_partner(is_first, is_last, partner)
       If ( is_first ) cx(0,:,:) = 0d0
    End If
    If ( z_bc_type == 1 ) Then
       Call z_periodic_partner(is_first, is_last, partner)
       If ( is_first ) cz(:,:,0) = 0d0
       If ( is_last ) cz(:,:,nk) = 0d0
    End If

  End Subroutine mg_fine_coef


  !> Coarse conductances: sum of the fine conductances through a coarse face, times the centre-distance ratio fr (and mg_cscale)
  Subroutine mg_coarsen_coef(nif, njf, nkf, cxf, cyf, czf, nic, njc, nkc, sx, sy, sz, fxr, fyr, fzr, cxc, cyc, czc)

    Integer(Int32), Intent(In) :: nif, njf, nkf, nic, njc, nkc, sx, sy, sz
    Real(Int64), Intent(In) :: cxf(0:nif+1,0:njf+1,0:nkf+1), cyf(0:nif+1,0:njf+1,0:nkf+1), czf(0:nif+1,0:njf+1,0:nkf+1)
    Real(Int64), Intent(In) :: fxr(0:nic), fyr(0:njc), fzr(0:nkc)
    Real(Int64), Intent(Out) :: cxc(0:nic+1,0:njc+1,0:nkc+1), cyc(0:nic+1,0:njc+1,0:nkc+1), czc(0:nic+1,0:njc+1,0:nkc+1)
    Integer(Int32) :: I, J, K, ii, jj, kk, lc
    Real(Int64) :: s

    cxc = 0d0;  cyc = 0d0;  czc = 0d0
    Do K = 1, nkc
       Do J = 1, njc
          Do I = 0, nic
             lc = Min(sx*I, nif);  s = 0d0
             Do kk = sz*K - (sz-1), Min(sz*K, nkf)
                Do jj = sy*J - (sy-1), Min(sy*J, njf)
                   s = s + cxf(lc,jj,kk)
                End Do
             End Do
             cxc(I,J,K) = mg_cscale*fxr(I)*s
          End Do
       End Do
    End Do
    Do K = 1, nkc
       Do J = 0, njc
          Do I = 1, nic
             lc = Min(sy*J, njf);  s = 0d0
             Do kk = sz*K - (sz-1), Min(sz*K, nkf)
                Do ii = sx*I - (sx-1), Min(sx*I, nif)
                   s = s + cyf(ii,lc,kk)
                End Do
             End Do
             cyc(I,J,K) = mg_cscale*fyr(J)*s
          End Do
       End Do
    End Do
    Do K = 0, nkc
       Do J = 1, njc
          Do I = 1, nic
             lc = Min(sz*K, nkf);  s = 0d0
             Do jj = sy*J - (sy-1), Min(sy*J, njf)
                Do ii = sx*I - (sx-1), Min(sx*I, nif)
                   s = s + czf(ii,jj,lc)
                End Do
             End Do
             czc(I,J,K) = mg_cscale*fzr(K)*s
          End Do
       End Do
    End Do

  End Subroutine mg_coarsen_coef


  !> Diagonal of S and the factors of the y-line tridiagonal T = (1+eps) diag - offdiag_y (Dirichlet reflection of the x outlet
  !  folded into the diagonal); a zero diagonal (solid cell) gives a zero pivot inverse, so the cell is left alone
  Subroutine mg_factor(ni, nj, nk, cx, cy, cz, dd, m, iw)

    Integer(Int32), Intent(In) :: ni, nj, nk
    Real(Int64), Intent(In) :: cx(0:ni+1,0:nj+1,0:nk+1), cy(0:ni+1,0:nj+1,0:nk+1), cz(0:ni+1,0:nj+1,0:nk+1)
    Real(Int64), Intent(Out) :: dd(0:ni+1,0:nj+1,0:nk+1), m(0:ni+1,0:nj+1,0:nk+1), iw(0:ni+1,0:nj+1,0:nk+1)
    Integer(Int32) :: i, j, k
    Real(Int64) :: dt, w, mm
    Logical :: dirichlet_hi

    dirichlet_hi = ( x_bc_type == 1 .And. mg_xup == MPI_PROC_NULL )
    dd = 0d0;  m = 0d0;  iw = 0d0
    Do k = 1, nk
       Do j = 1, nj
          Do i = 1, ni
             dd(i,j,k) = cx(i,j,k) + cx(i-1,j,k) + cy(i,j,k) + cy(i,j-1,k) + cz(i,j,k) + cz(i,j,k-1)
          End Do
       End Do
    End Do
    Do k = 1, nk
       Do i = 1, ni
          w = 0d0
          Do j = 1, nj
             dt = (1d0 + mg_eps)*dd(i,j,k)
             If ( dirichlet_hi .And. i == ni ) dt = dt + (1d0 + mg_eps)*cx(i,j,k)
             mm = 0d0
             If ( j > 1 ) Then
                If ( cy(i,j-1,k) > 0d0 .And. w > 0d0 ) mm = cy(i,j-1,k)/w
                dt = dt - mm*cy(i,j-1,k)
             End If
             w = dt
             m(i,j,k) = mm
             If ( w > 0d0 ) iw(i,j,k) = 1d0/w
          End Do
       End Do
    End Do

  End Subroutine mg_factor


  !> Ghost layers of a level array: neighbour ranks and periodic partners by MPI, Neumann copies at the y walls and at physical
  !  x/z boundaries, the reflection mg_xhi_gs at a Dirichlet x outlet
  Subroutine mg_halo(a, ni, nj, nk)

    Integer(Int32), Intent(In) :: ni, nj, nk
    Real(Int64), Intent(InOut) :: a(0:ni+1,0:nj+1,0:nk+1)
    Integer(Int32) :: i, j, k, n
    Real(Int64) :: glo, ghi

    !$acc kernels present(a) async(1)
    a(:,0,:) = a(:,1,:)
    a(:,nj+1,:) = a(:,nj,:)
    !$acc end kernels
    If ( mg_xself ) Then
       !$acc kernels present(a) async(1)
       a(0,:,:) = a(ni,:,:)
       a(ni+1,:,:) = a(1,:,:)
       !$acc end kernels
    Else
       n = (nj+2)*(nk+2)
       !$acc parallel loop collapse(2) present(a,mg_sbx) async(1)
       Do k = 0, nk+1
          Do j = 0, nj+1
             mg_sbx(1+j+(nj+2)*k,1) = a(ni,j,k)
             mg_sbx(1+j+(nj+2)*k,2) = a(1,j,k)
          End Do
       End Do
       !$acc wait(1)
       !$acc update host(mg_sbx(1:n,1:2))
       Call MPI_Sendrecv(mg_sbx(1,1), n, MPI_real8, mg_xup, 301, mg_rbx(1,1), n, MPI_real8, mg_xdn, 301, &
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
       Call MPI_Sendrecv(mg_sbx(1,2), n, MPI_real8, mg_xdn, 302, mg_rbx(1,2), n, MPI_real8, mg_xup, 302, &
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
       !$acc update device(mg_rbx(1:n,1:2))
       glo = mg_xlo_gs;  ghi = mg_xhi_gs
       !$acc parallel loop collapse(2) present(a,mg_rbx) firstprivate(glo,ghi) async(1)
       Do k = 0, nk+1
          Do j = 0, nj+1
             If ( mg_xdn /= MPI_PROC_NULL ) Then
                a(0,j,k) = mg_rbx(1+j+(nj+2)*k,1)
             Else
                a(0,j,k) = glo*a(1,j,k)
             End If
             If ( mg_xup /= MPI_PROC_NULL ) Then
                a(ni+1,j,k) = mg_rbx(1+j+(nj+2)*k,2)
             Else
                a(ni+1,j,k) = ghi*a(ni,j,k)
             End If
          End Do
       End Do
    End If
    If ( mg_zself ) Then
       !$acc kernels present(a) async(1)
       a(:,:,0) = a(:,:,nk)
       a(:,:,nk+1) = a(:,:,1)
       !$acc end kernels
    Else
       n = (ni+2)*(nj+2)
       !$acc parallel loop collapse(2) present(a,mg_sbz) async(1)
       Do j = 0, nj+1
          Do i = 0, ni+1
             mg_sbz(1+i+(ni+2)*j,1) = a(i,j,nk)
             mg_sbz(1+i+(ni+2)*j,2) = a(i,j,1)
          End Do
       End Do
       !$acc wait(1)
       !$acc update host(mg_sbz(1:n,1:2))
       Call MPI_Sendrecv(mg_sbz(1,1), n, MPI_real8, mg_zup, 303, mg_rbz(1,1), n, MPI_real8, mg_zdn, 303, &
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
       Call MPI_Sendrecv(mg_sbz(1,2), n, MPI_real8, mg_zdn, 304, mg_rbz(1,2), n, MPI_real8, mg_zup, 304, &
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
       !$acc update device(mg_rbz(1:n,1:2))
       !$acc parallel loop collapse(2) present(a,mg_rbz) async(1)
       Do j = 0, nj+1
          Do i = 0, ni+1
             If ( mg_zdn /= MPI_PROC_NULL ) Then
                a(i,j,0) = mg_rbz(1+i+(ni+2)*j,1)
             Else
                a(i,j,0) = a(i,j,1)
             End If
             If ( mg_zup /= MPI_PROC_NULL ) Then
                a(i,j,nk+1) = mg_rbz(1+i+(ni+2)*j,2)
             Else
                a(i,j,nk+1) = a(i,j,nk)
             End If
          End Do
       End Do
    End If

  End Subroutine mg_halo


  !> t = b - S x on the unknown cells (x with valid ghost layers)
  Subroutine mg_resid(ni, nj, nk, x, b, cx, cy, cz, dd, t)

    Integer(Int32), Intent(In) :: ni, nj, nk
    Real(Int64), Intent(In) :: x(0:ni+1,0:nj+1,0:nk+1), b(0:ni+1,0:nj+1,0:nk+1), dd(0:ni+1,0:nj+1,0:nk+1)
    Real(Int64), Intent(In) :: cx(0:ni+1,0:nj+1,0:nk+1), cy(0:ni+1,0:nj+1,0:nk+1), cz(0:ni+1,0:nj+1,0:nk+1)
    Real(Int64), Intent(InOut) :: t(0:ni+1,0:nj+1,0:nk+1)
    Integer(Int32) :: i, j, k

    !$acc parallel loop collapse(3) present(x,b,cx,cy,cz,dd,t) async(1)
    Do k = 1, nk
       Do j = 1, nj
          Do i = 1, ni
             t(i,j,k) = b(i,j,k) - ( dd(i,j,k)*x(i,j,k) - cx(i,j,k)*x(i+1,j,k) - cx(i-1,j,k)*x(i-1,j,k) &
                                   - cy(i,j,k)*x(i,j+1,k) - cy(i,j-1,k)*x(i,j-1,k) &
                                   - cz(i,j,k)*x(i,j,k+1) - cz(i,j,k-1)*x(i,j,k-1) )
          End Do
       End Do
    End Do

  End Subroutine mg_resid


  !> t <- T^-1 t: one tridiagonal solve along y per column
  Subroutine mg_line(ni, nj, nk, t, cy, m, iw)

    Integer(Int32), Intent(In) :: ni, nj, nk
    Real(Int64), Intent(InOut) :: t(0:ni+1,0:nj+1,0:nk+1)
    Real(Int64), Intent(In) :: cy(0:ni+1,0:nj+1,0:nk+1), m(0:ni+1,0:nj+1,0:nk+1), iw(0:ni+1,0:nj+1,0:nk+1)
    Integer(Int32) :: i, j, k

#ifdef GPU_POISSON
    !$acc parallel loop collapse(2) gang vector present(t,cy,m,iw) async(1)
    Do k = 1, nk
       Do i = 1, ni
          !$acc loop seq
          Do j = 2, nj
             t(i,j,k) = t(i,j,k) + m(i,j,k)*t(i,j-1,k)
          End Do
          t(i,nj,k) = t(i,nj,k)*iw(i,nj,k)
          !$acc loop seq
          Do j = nj-1, 1, -1
             t(i,j,k) = ( t(i,j,k) + cy(i,j,k)*t(i,j+1,k) )*iw(i,j,k)
          End Do
       End Do
    End Do
#else
    Do k = 1, nk
       Do j = 2, nj
          Do i = 1, ni
             t(i,j,k) = t(i,j,k) + m(i,j,k)*t(i,j-1,k)
          End Do
       End Do
       Do i = 1, ni
          t(i,nj,k) = t(i,nj,k)*iw(i,nj,k)
       End Do
       Do j = nj-1, 1, -1
          Do i = 1, ni
             t(i,j,k) = ( t(i,j,k) + cy(i,j,k)*t(i,j+1,k) )*iw(i,j,k)
          End Do
       End Do
    End Do
#endif

  End Subroutine mg_line


  !> d = c1 d + c2 t; x = x + d (x = d when zero)
  Subroutine mg_upd(ni, nj, nk, x, d, t, c1, c2, zero)

    Integer(Int32), Intent(In) :: ni, nj, nk
    Real(Int64), Intent(InOut) :: x(0:ni+1,0:nj+1,0:nk+1), d(0:ni+1,0:nj+1,0:nk+1)
    Real(Int64), Intent(In) :: t(0:ni+1,0:nj+1,0:nk+1), c1, c2
    Logical, Intent(In) :: zero
    Integer(Int32) :: i, j, k

    !$acc parallel loop collapse(3) present(x,d,t) async(1)
    Do k = 1, nk
       Do j = 1, nj
          Do i = 1, ni
             d(i,j,k) = c1*d(i,j,k) + c2*t(i,j,k)
             If ( zero ) Then
                x(i,j,k) = d(i,j,k)
             Else
                x(i,j,k) = x(i,j,k) + d(i,j,k)
             End If
          End Do
       End Do
    End Do

  End Subroutine mg_upd


  !> Chebyshev polynomial of degree deg in T^-1 S on [lo, hi] (hi = 2 bounds the line-Jacobi spectrum), from x = 0 when zero
  Subroutine mg_cheb(ni, nj, nk, x, b, t, d, cx, cy, cz, dd, m, iw, deg, lo, hi, zero)

    Integer(Int32), Intent(In) :: ni, nj, nk, deg
    Real(Int64), Intent(InOut) :: x(0:ni+1,0:nj+1,0:nk+1), t(0:ni+1,0:nj+1,0:nk+1), d(0:ni+1,0:nj+1,0:nk+1)
    Real(Int64), Intent(In) :: b(0:ni+1,0:nj+1,0:nk+1), cx(0:ni+1,0:nj+1,0:nk+1), cy(0:ni+1,0:nj+1,0:nk+1)
    Real(Int64), Intent(In) :: cz(0:ni+1,0:nj+1,0:nk+1), dd(0:ni+1,0:nj+1,0:nk+1), m(0:ni+1,0:nj+1,0:nk+1)
    Real(Int64), Intent(In) :: iw(0:ni+1,0:nj+1,0:nk+1), lo, hi
    Logical, Intent(In) :: zero
    Integer(Int32) :: it
    Real(Int64) :: theta, delta, sigma, rho0, rho1

    theta = 0.5d0*(hi + lo);  delta = 0.5d0*(hi - lo);  sigma = theta/delta;  rho0 = 1d0/sigma
    Do it = 1, deg
       If ( it == 1 .And. zero ) Then
          Call mg_copy(ni, nj, nk, b, t)
       Else
          Call mg_halo(x, ni, nj, nk)
          Call mg_resid(ni, nj, nk, x, b, cx, cy, cz, dd, t)
       End If
       Call mg_line(ni, nj, nk, t, cy, m, iw)
       If ( it == 1 ) Then
          Call mg_upd(ni, nj, nk, x, d, t, 0d0, 1d0/theta, zero)
       Else
          rho1 = 1d0/(2d0*sigma - rho0)
          Call mg_upd(ni, nj, nk, x, d, t, rho1*rho0, 2d0*rho1/delta, .False.)
          rho0 = rho1
       End If
    End Do

  End Subroutine mg_cheb


  Subroutine mg_copy(ni, nj, nk, a, c)

    Integer(Int32), Intent(In) :: ni, nj, nk
    Real(Int64), Intent(In) :: a(0:ni+1,0:nj+1,0:nk+1)
    Real(Int64), Intent(InOut) :: c(0:ni+1,0:nj+1,0:nk+1)
    Integer(Int32) :: i, j, k

    !$acc parallel loop collapse(3) present(a,c) async(1)
    Do k = 1, nk
       Do j = 1, nj
          Do i = 1, ni
             c(i,j,k) = a(i,j,k)
          End Do
       End Do
    End Do

  End Subroutine mg_copy


  !> bc = sum of the fine values under each coarse cell (transpose of the piecewise-constant prolongation)
  Subroutine mg_restrict(nif, njf, nkf, tf, nic, njc, nkc, sx, sy, sz, bc)

    Integer(Int32), Intent(In) :: nif, njf, nkf, nic, njc, nkc, sx, sy, sz
    Real(Int64), Intent(In) :: tf(0:nif+1,0:njf+1,0:nkf+1)
    Real(Int64), Intent(InOut) :: bc(0:nic+1,0:njc+1,0:nkc+1)
    Integer(Int32) :: I, J, K, ii, jj, kk
    Real(Int64) :: s

    !$acc parallel loop collapse(3) gang vector present(tf,bc) async(1)
    Do K = 1, nkc
       Do J = 1, njc
          Do I = 1, nic
             s = 0d0
             Do kk = sz*K - (sz-1), Min(sz*K, nkf)
                Do jj = sy*J - (sy-1), Min(sy*J, njf)
                   Do ii = sx*I - (sx-1), Min(sx*I, nif)
                      s = s + tf(ii,jj,kk)
                   End Do
                End Do
             End Do
             bc(I,J,K) = s
          End Do
       End Do
    End Do

  End Subroutine mg_restrict


  Subroutine mg_prolong(nif, njf, nkf, xf, nic, njc, nkc, sx, sy, sz, xc)

    Integer(Int32), Intent(In) :: nif, njf, nkf, nic, njc, nkc, sx, sy, sz
    Real(Int64), Intent(InOut) :: xf(0:nif+1,0:njf+1,0:nkf+1)
    Real(Int64), Intent(In) :: xc(0:nic+1,0:njc+1,0:nkc+1)
    Integer(Int32) :: i, j, k

    !$acc parallel loop collapse(3) present(xf,xc) async(1)
    Do k = 1, nkf
       Do j = 1, njf
          Do i = 1, nif
             xf(i,j,k) = xf(i,j,k) + xc((i+sx-1)/sx,(j+sy-1)/sy,(k+sz-1)/sz)
          End Do
       End Do
    End Do

  End Subroutine mg_prolong


  !> One symmetric V-cycle from a zero guess on the level-1 right-hand side mg_b
  Subroutine mg_vcycle

    Integer(Int32) :: l, o, oc
    Real(Int64), Parameter :: hi = 2d0

    Do l = 1, mg_nlev-1
       o = mg_off(l) + 1;  oc = mg_off(l+1) + 1
       Call mg_cheb(mg_ni(l), mg_nj(l), mg_nk(l), mg_x(o), mg_b(o), mg_t(o), mg_d(o), mg_cx(o), mg_cy(o), mg_cz(o), &
                    mg_dd(o), mg_m(o), mg_iw(o), mg_deg, hi/mg_ratio, hi, .True.)
       Call mg_halo(mg_x(o), mg_ni(l), mg_nj(l), mg_nk(l))
       Call mg_resid(mg_ni(l), mg_nj(l), mg_nk(l), mg_x(o), mg_b(o), mg_cx(o), mg_cy(o), mg_cz(o), mg_dd(o), mg_t(o))
       Call mg_restrict(mg_ni(l), mg_nj(l), mg_nk(l), mg_t(o), mg_ni(l+1), mg_nj(l+1), mg_nk(l+1), &
                        mg_sx(l), mg_sy(l), mg_sz(l), mg_b(oc))
    End Do
    l = mg_nlev;  o = mg_off(l) + 1
    Call mg_cheb(mg_ni(l), mg_nj(l), mg_nk(l), mg_x(o), mg_b(o), mg_t(o), mg_d(o), mg_cx(o), mg_cy(o), mg_cz(o), &
                 mg_dd(o), mg_m(o), mg_iw(o), mg_deg_coarse, hi/mg_ratio_coarse, hi, .True.)
    Do l = mg_nlev-1, 1, -1
       o = mg_off(l) + 1;  oc = mg_off(l+1) + 1
       Call mg_prolong(mg_ni(l), mg_nj(l), mg_nk(l), mg_x(o), mg_ni(l+1), mg_nj(l+1), mg_nk(l+1), &
                       mg_sx(l), mg_sy(l), mg_sz(l), mg_x(oc))
       Call mg_cheb(mg_ni(l), mg_nj(l), mg_nk(l), mg_x(o), mg_b(o), mg_t(o), mg_d(o), mg_cx(o), mg_cy(o), mg_cz(o), &
                    mg_dd(o), mg_m(o), mg_iw(o), mg_deg, hi/mg_ratio, hi, .False.)
    End Do

  End Subroutine mg_vcycle


  !> z = M^-1 r = -S^-1 (V r) for fine cell-centred vectors r, z of the shape (nxg,nyg,nzg)
  Subroutine mg_precond(r, z)

    Real(Int64), Intent(In)  :: r(nxg,nyg,nzg)
    Real(Int64), Intent(Out) :: z(nxg,nyg,nzg)

    Call mg_load(mg_ni(1), mg_nj(1), mg_nk(1), r, mg_hy1, mg_hz1, mg_b)
    Call mg_vcycle
    Call mg_store(mg_ni(1), mg_nj(1), mg_nk(1), mg_x, z)
    !$acc wait(1)

  End Subroutine mg_precond


  Subroutine mg_load(ni, nj, nk, r, hy, hz, b)

    Integer(Int32), Intent(In) :: ni, nj, nk
    Real(Int64), Intent(In) :: r(nxg,nyg,nzg), hy(0:nj+1), hz(0:nk+1)
    Real(Int64), Intent(InOut) :: b(0:ni+1,0:nj+1,0:nk+1)
    Integer(Int32) :: i, j, k

    !$acc parallel loop collapse(3) present(r,hy,hz,b) async(1)
    Do k = 1, nk
       Do j = 1, nj
          Do i = 1, ni
             b(i,j,k) = -dx*hy(j)*hz(k)*r(i+1,j+1,k+1)
          End Do
       End Do
    End Do

  End Subroutine mg_load


  Subroutine mg_store(ni, nj, nk, x, z)

    Integer(Int32), Intent(In) :: ni, nj, nk
    Real(Int64), Intent(In) :: x(0:ni+1,0:nj+1,0:nk+1)
    Real(Int64), Intent(InOut) :: z(nxg,nyg,nzg)
    Integer(Int32) :: i, j, k

    !$acc parallel loop collapse(3) present(x,z) async(1)
    Do k = 1, nzg
       Do j = 1, nyg
          Do i = 1, nxg
             z(i,j,k) = 0d0
          End Do
       End Do
    End Do
    !$acc parallel loop collapse(3) present(x,z) async(1)
    Do k = 1, nk
       Do j = 1, nj
          Do i = 1, ni
             z(i+1,j+1,k+1) = x(i,j,k)
          End Do
       End Do
    End Do

  End Subroutine mg_store

End Module vof_mg
