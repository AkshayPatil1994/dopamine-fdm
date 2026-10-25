!> Mask-aware geometric multigrid preconditioner for the masked PCG of vof_pressure (pcg_precond = 1). The
!  operator is the symmetric conductance form S = -V A of A = div(beta grad) (V the cell volumes), so z = A^-1
!  r = -S^-1 (V r) and the preconditioner is symmetric negative definite like A. Coarse cells aggregate up to
!  2x2x2 rank-local fine cells (strides mg_sx/mg_sy/mg_sz); transfers are piecewise constant; coarse
!  conductances sum the fine ones crossing a coarse face, scaled by the centre-distance ratio. A coarse face
!  is open if any fine face under it is; solid cells and sealed pockets give empty or singular rows the
!  smoother leaves alone. The smoother is a Chebyshev polynomial of the y-line-Jacobi iteration (exact
!  tridiagonal solves along the stretched y direction, which carries the anisotropy), identical before and
!  after the coarse correction, so the V-cycle from a zero guess is a fixed symmetric operator.
Module vof_mg

  Use iso_fortran_env, Only : Int32, Int64, Real32
  Use global
  Use mpi
  Use decomp, Only : x_halo_neighbors, z_halo_neighbors, x_periodic_partner, z_periodic_partner

  Implicit None

  ! storage kind of the V-cycle (the outer PCG stays in double): the preconditioner only needs ~1e-7 relative
  ! accuracy and is memory-bound
  Integer(Int32), Parameter :: mgk = Real32
  Integer(Int32) :: mg_mpi_t
  Integer(Int32), Parameter :: mg_maxlev = 24
  Integer(Int32) :: mg_nlev = 0, mg_deg = 3, mg_deg_coarse = 8
  Real(Int64) :: mg_ratio = 6d0, mg_ratio_coarse = 60d0, mg_cscale = 1d0, mg_eps = 1d-4
  Logical :: mg_coarsen_y = .True.
  Integer(Int32) :: mg_ni(mg_maxlev), mg_nj(mg_maxlev), mg_nk(mg_maxlev)
  Integer(Int32) :: mg_sx(mg_maxlev), mg_sy(mg_maxlev), mg_sz(mg_maxlev)   ! coarsening stride from level l to l+1
  Integer(Int32) :: mg_off(mg_maxlev+1)
  ! flat storage, one offset (mg_off) per level, indexed (0:ni+1,0:nj+1,0:nk+1): cells 1:n are the unknowns,
  ! 0 and n+1 the ghost layer (neighbour rank, periodic partner, Neumann copy or Dirichlet reflection of the
  ! x outlet); cx(i) joins cells i, i+1
  Real(mgk), Allocatable, Dimension(:) :: mg_cx, mg_cy, mg_cz, mg_dd, mg_x, mg_b, mg_t, mg_d, mg_m, mg_iw, mg_s
  ! centre-distance ratios of the coarse faces, level l -> l+1
  Real(Int64), Allocatable, Dimension(:,:) :: mg_fx, mg_fy, mg_fz
  Real(Int64), Allocatable, Dimension(:) :: mg_hy1, mg_hz1, mg_dcy1, mg_dcz1   ! level-1 cell sizes and centre distances
  Real(mgk), Allocatable, Dimension(:,:) :: mg_sbx, mg_rbx, mg_sbz, mg_rbz
  Integer(Int32) :: mg_xup, mg_xdn, mg_zup, mg_zdn
  Logical :: mg_xself, mg_zself
  Real(Int64) :: mg_xlo_gs, mg_xhi_gs   ! ghost = gs * boundary cell where there is no neighbour
  Logical :: mg_ready = .False.
  ! direct coarsest level: the coarsest level (all ranks together, nl own cells of ng) is solved exactly by the
  ! explicit inverse rows of its banded matrix, so the cycle does not depend on how many sweeps it can afford
  Logical :: mg_direct = .False.
  Integer(Int32) :: mg_dnl = 0, mg_dng = 0, mg_db = 1, mg_dnx = 1, mg_dnz = 1
  Integer(Int32), Allocatable :: mg_dcnt(:), mg_ddsp(:), mg_dpos(:), mg_dgi(:,:,:)
  Real(mgk), Allocatable :: mg_dinv(:,:), mg_drecv(:), mg_dout(:), mg_dpk(:,:), mg_dcg(:,:)

Contains


  !> Hierarchy, neighbour topology and level-1 geometry; ihi, jhi, khi are the last unknown fine cells
  Subroutine mg_init(ihi, jhi, khi)

    Integer(Int32), Intent(In) :: ihi, jhi, khi

    Integer(Int32) :: l, j, k, nmin(2), nmin_g(2), nsum(3), nsum_g(3), ntot, nmax, partner
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
    mg_mpi_t = Merge(MPI_REAL4, MPI_REAL8, mgk == Real32)
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
       nsum = (/ mg_ni(l), mg_nk(l), mg_ni(l)*mg_nj(l)*mg_nk(l) /)
       Call MPI_Allreduce(nsum, nsum_g, 3, MPI_integer, MPI_SUM, MPI_COMM_WORLD, ierr)
       mg_dnx = nsum_g(1)/p_col;  mg_dnz = nsum_g(2)/p_row;  mg_dng = nsum_g(3)
       mg_db = mg_direct_band(mg_nj(l), mg_dnx, mg_dnz)
       ! band Cholesky of the n x n matrix plus nl solves for the own rows of its inverse
       If ( mg_dng <= 1024 .And. Real(mg_dng,Int64)*Real(mg_db,Int64)*( Real(mg_db,Int64) &
            + 2d0*Real((mg_dng+p_row*p_col-1)/(p_row*p_col),Int64) ) <= 6d6 ) Then
          mg_direct = .True.
          Exit
       End If
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

    If ( mg_direct ) Call mg_direct_setup
    Allocate( mg_cx(ntot), mg_cy(ntot), mg_cz(ntot), mg_dd(ntot), mg_x(ntot), mg_b(ntot), mg_t(ntot), mg_d(ntot), &
              mg_m(ntot), mg_iw(ntot), mg_s(ntot) )
    mg_cx = 0d0;  mg_cy = 0d0;  mg_cz = 0d0;  mg_dd = 0d0;  mg_x = 0d0;  mg_b = 0d0;  mg_t = 0d0;  mg_d = 0d0
    mg_m = 0d0;  mg_iw = 0d0;  mg_s = 0d0
    mg_deg = mg_smooth_deg
    Allocate( mg_sbx((mg_nj(1)+2)*(mg_nk(1)+2),2), mg_rbx((mg_nj(1)+2)*(mg_nk(1)+2),2) )
    Allocate( mg_sbz((mg_ni(1)+2)*(mg_nj(1)+2),2), mg_rbz((mg_ni(1)+2)*(mg_nj(1)+2),2) )
    !$acc enter data copyin(mg_cx,mg_cy,mg_cz,mg_dd,mg_x,mg_b,mg_t,mg_d,mg_m,mg_iw,mg_s)
    !$acc enter data create(mg_sbx,mg_rbx,mg_sbz,mg_rbz)
    !$acc enter data copyin(mg_hy1,mg_hz1,mg_fx,mg_fy,mg_fz)
    mg_ready = .True.

    If ( myid == 0 ) Then
       Write(*,'(A,I0,A)') '   GMG preconditioner: ', mg_nlev, ' levels (rank-0 sizes ni,nj,nk):'
       Do l = 1, mg_nlev
          Write(*,'(A,I3,3I6)') '     level', l, mg_ni(l), mg_nj(l), mg_nk(l)
       End Do
       If ( mg_direct ) Write(*,'(A,I0,A,I0)') '   coarsest level solved directly: ', mg_dng, ' cells, half-bandwidth ', mg_db
    End If

  End Subroutine mg_init



  !> Half-bandwidth of the coarsest-level matrix, ordered j + nj*(fold(i) - 1 + nx*(fold(k) - 1))
  Function mg_direct_band(nj, nx, nz) Result(b)

    Integer(Int32), Intent(In) :: nj, nx, nz
    Integer(Int32) :: b

    b = 1
    If ( nx > 1 ) b = Max(b, nj*Merge(2, 1, x_bc_type == 0 .And. nx > 2))
    If ( nz > 1 ) b = Max(b, nj*nx*Merge(2, 1, z_bc_type == 0 .And. nz > 2))

  End Function mg_direct_band


  !> Position of cell i of n: periodic directions are folded (1, n, 2, n-1, ...) so that the wrap-around
  !  neighbours stay within two positions
  Function mg_fold(i, n, per) Result(p)

    Integer(Int32), Intent(In) :: i, n
    Logical, Intent(In) :: per
    Integer(Int32) :: p

    p = i
    If ( per .And. n > 2 ) p = Merge(2*i - 1, 2*(n - i + 1), i <= (n + 1)/2)

  End Function mg_fold


  !> Gather layout of the coarsest level over all ranks (row-major rank = x index * p_col + z index) and the
  !  band position of every gathered cell
  Subroutine mg_direct_setup

    Integer(Int32) :: l, np, r, row, col, rr, i, j, k, ni, nj, nk, xoff, zoff, idx
    Integer(Int32), Allocatable :: nia(:), nka(:)

    np = p_row*p_col
    l = mg_nlev
    ni = mg_ni(l);  nj = mg_nj(l);  nk = mg_nk(l)
    Allocate( nia(0:np-1), nka(0:np-1), mg_dcnt(0:np-1), mg_ddsp(0:np-1) )
    Call MPI_Allgather(ni, 1, MPI_integer, nia, 1, MPI_integer, MPI_COMM_WORLD, ierr)
    Call MPI_Allgather(nk, 1, MPI_integer, nka, 1, MPI_integer, MPI_COMM_WORLD, ierr)
    mg_dcnt = nia*nj*nka
    mg_ddsp(0) = 0
    Do r = 1, np-1
       mg_ddsp(r) = mg_ddsp(r-1) + mg_dcnt(r-1)
    End Do
    Allocate( mg_dpos(mg_dng), mg_dgi(mg_dnx,nj,mg_dnz) )
    Do r = 0, np-1
       row = r/p_col;  col = Mod(r, p_col)
       xoff = 0
       Do rr = 0, row-1
          xoff = xoff + nia(rr*p_col)
       End Do
       zoff = Sum(nka(row*p_col:row*p_col+col-1))
       Do k = 1, nka(r)
          Do j = 1, nj
             Do i = 1, nia(r)
                idx = mg_ddsp(r) + i + nia(r)*((j-1) + nj*(k-1))
                mg_dgi(xoff+i,j,zoff+k) = idx
                mg_dpos(idx) = j + nj*( mg_fold(xoff+i, mg_dnx, x_bc_type == 0) - 1 &
                                      + mg_dnx*(mg_fold(zoff+k, mg_dnz, z_bc_type == 0) - 1) )
             End Do
          End Do
       End Do
    End Do
    mg_dnl = ni*nj*nk
    Allocate( mg_dinv(mg_dng,mg_dnl), mg_drecv(mg_dng), mg_dout(mg_dnl), mg_dpk(3,mg_dnl), mg_dcg(3,mg_dng) )
    mg_dinv = 0d0;  mg_drecv = 0d0;  mg_dout = 0d0
    !$acc enter data create(mg_dinv,mg_drecv,mg_dout,mg_dpk)
    Deallocate( nia, nka )

  End Subroutine mg_direct_setup


  !> Adds the face (p, q) of conductance c to the banded lower storage ab (off-diagonal and, when p == q, the
  !  diagonal part) and to the face sum ds
  Subroutine mg_dpair(ab, ds, b, n, p, q, c)

    Integer(Int32), Intent(In) :: b, n, p, q
    Real(Int64), Intent(InOut) :: ab(0:b,n), ds(n)
    Real(Int64), Intent(In) :: c

    ds(p) = ds(p) + c;  ds(q) = ds(q) + c
    If ( p == q ) Then
       ab(0,p) = ab(0,p) - 2d0*c
    Else
       ab(Abs(p-q),Min(p,q)) = ab(Abs(p-q),Min(p,q)) - c
    End If

  End Subroutine mg_dpair


  !> Gathers the coarsest-level conductances, assembles K = S + eps*diag(S) (the Dirichlet x outlet folded
  !  into the diagonal) in a band, factors it and stores the rows of K^-1 of the own cells (K is symmetric).
  !  Cells without an open face stay zero.
  Subroutine mg_direct_build

    Integer(Int32) :: l, o, np, n, b, nl, ni, nj, nk, ix, j, kz, ix2, kz2, idx, p, q, pc, c, g, r
    Real(Int64), Allocatable :: ab(:,:), ds(:), dr(:), y(:), xs(:)
    Real(Int64) :: s, cf

    np = p_row*p_col
    l = mg_nlev;  o = mg_off(l) + 1
    ni = mg_ni(l);  nj = mg_nj(l);  nk = mg_nk(l)
    n = mg_dng;  b = mg_db;  nl = mg_dnl
    Call mg_dpack(ni, nj, nk, mg_cx(o), mg_cy(o), mg_cz(o), mg_dpk)
    !$acc wait(1)
    !$acc update host(mg_dpk)
    If ( np == 1 ) Then
       mg_dcg = mg_dpk
    Else
       Call MPI_Allgatherv(mg_dpk, 3*nl, mg_mpi_t, mg_dcg, 3*mg_dcnt, 3*mg_ddsp, mg_mpi_t, MPI_COMM_WORLD, ierr)
    End If

    Allocate( ab(0:b,n), ds(n), dr(n), y(n), xs(n) )
    ab = 0d0;  ds = 0d0;  dr = 0d0
    Do kz = 1, mg_dnz
       Do j = 1, nj
          Do ix = 1, mg_dnx
             idx = mg_dgi(ix,j,kz);  p = mg_dpos(idx)
             cf = mg_dcg(1,idx)
             If ( cf > 0d0 ) Then
                ix2 = ix + 1
                If ( ix2 > mg_dnx .And. x_bc_type == 0 ) ix2 = 1
                If ( ix2 > mg_dnx ) Then
                   ds(p) = ds(p) + cf;  dr(p) = cf
                Else
                   Call mg_dpair(ab, ds, b, n, p, mg_dpos(mg_dgi(ix2,j,kz)), cf)
                End If
             End If
             cf = mg_dcg(2,idx)
             If ( cf > 0d0 .And. j < nj ) Call mg_dpair(ab, ds, b, n, p, mg_dpos(mg_dgi(ix,j+1,kz)), cf)
             cf = mg_dcg(3,idx)
             If ( cf > 0d0 ) Then
                kz2 = kz + 1
                If ( kz2 > mg_dnz .And. z_bc_type == 0 ) kz2 = 1
                If ( kz2 <= mg_dnz ) Call mg_dpair(ab, ds, b, n, p, mg_dpos(mg_dgi(ix,j,kz2)), cf)
             End If
          End Do
       End Do
    End Do
    Do p = 1, n
       If ( ds(p) > 0d0 ) Then
          ab(0,p) = ab(0,p) + (1d0 + mg_eps)*ds(p) + dr(p)
       Else
          ab(0,p) = 1d0
       End If
    End Do
    Do p = 1, n
       ab(0,p) = Sqrt(ab(0,p))
       Do r = p+1, Min(n, p+b)
          ab(r-p,p) = ab(r-p,p)/ab(0,p)
       End Do
       Do q = p+1, Min(n, p+b)
          Do r = q, Min(n, p+b)
             ab(r-q,q) = ab(r-q,q) - ab(r-p,p)*ab(q-p,p)
          End Do
       End Do
    End Do
    mg_dinv = 0d0
    Do c = 1, nl
       pc = mg_dpos(mg_ddsp(myid) + c)
       If ( ds(pc) <= 0d0 ) Cycle
       y = 0d0
       y(pc) = 1d0/ab(0,pc)
       Do r = pc+1, n
          s = 0d0
          Do q = Max(pc, r-b), r-1
             s = s + ab(r-q,q)*y(q)
          End Do
          y(r) = -s/ab(0,r)
       End Do
       Do p = n, 1, -1
          s = y(p)
          Do q = p+1, Min(n, p+b)
             s = s - ab(q-p,p)*xs(q)
          End Do
          xs(p) = s/ab(0,p)
       End Do
       Do g = 1, n
          mg_dinv(g,c) = Real(xs(mg_dpos(g)), mgk)
       End Do
    End Do
    !$acc update device(mg_dinv)
    Deallocate( ab, ds, dr, y, xs )

  End Subroutine mg_direct_build


  Subroutine mg_dpack(ni, nj, nk, cx, cy, cz, pk)

    Integer(Int32), Intent(In) :: ni, nj, nk
    Real(mgk), Intent(In) :: cx(0:ni+1,0:nj+1,0:nk+1), cy(0:ni+1,0:nj+1,0:nk+1), cz(0:ni+1,0:nj+1,0:nk+1)
    Real(mgk), Intent(InOut) :: pk(3,ni*nj*nk)
    Integer(Int32) :: i, j, k, c

    !$acc parallel loop collapse(3) private(c) present(cx,cy,cz,pk) async(1)
    Do k = 1, nk
       Do j = 1, nj
          Do i = 1, ni
             c = i + ni*((j-1) + nj*(k-1))
             pk(1,c) = cx(i,j,k);  pk(2,c) = cy(i,j,k);  pk(3,c) = cz(i,j,k)
          End Do
       End Do
    End Do

  End Subroutine mg_dpack


  !> x = K^-1 b on the coarsest level: gather b over the ranks, one dot product with the stored inverse row per
  !  own cell
  Subroutine mg_direct_solve(ni, nj, nk, b, x)

    Integer(Int32), Intent(In) :: ni, nj, nk
    Real(mgk), Intent(In) :: b(0:ni+1,0:nj+1,0:nk+1)
    Real(mgk), Intent(InOut) :: x(0:ni+1,0:nj+1,0:nk+1)
    Integer(Int32) :: i, j, k, c, g, d0, nl
    Real(mgk) :: s

    nl = mg_dnl
    d0 = mg_ddsp(myid)
    !$acc parallel loop collapse(3) private(c) present(b,mg_drecv) async(1)
    Do k = 1, nk
       Do j = 1, nj
          Do i = 1, ni
             c = i + ni*((j-1) + nj*(k-1))
             mg_drecv(d0+c) = b(i,j,k)
          End Do
       End Do
    End Do
    If ( p_row*p_col > 1 ) Then
       !$acc wait(1)
       !$acc update host(mg_drecv(d0+1:d0+nl))
       Call MPI_Allgatherv(MPI_IN_PLACE, 0, mg_mpi_t, mg_drecv, mg_dcnt, mg_ddsp, mg_mpi_t, MPI_COMM_WORLD, ierr)
       !$acc update device(mg_drecv) async(1)
    End If
    !$acc parallel loop gang private(s) present(mg_dinv,mg_drecv,mg_dout) async(1)
    Do c = 1, nl
       s = 0.0
       !$acc loop vector reduction(+:s)
       Do g = 1, mg_dng
          s = s + mg_dinv(g,c)*mg_drecv(g)
       End Do
       mg_dout(c) = s
    End Do
    !$acc parallel loop collapse(3) private(c) present(x,mg_dout) async(1)
    Do k = 1, nk
       Do j = 1, nj
          Do i = 1, ni
             c = i + ni*((j-1) + nj*(k-1))
             x(i,j,k) = mg_dout(c)
          End Do
       End Do
    End Do

  End Subroutine mg_direct_solve




  !> Ghost sizes of a 1D size array h(1:n): the neighbour's end cell, or the own cell where there is no
  !  neighbour
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


  !> Sizes and centre distances of the coarse cells of one direction, and the ratio fr(I) of the centre
  !  distance of the fine face under each coarse face to that of the coarse face (conductance of a coarsened
  !  direction scales with it)
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


  !> Conductances of all levels from the fine face coefficients (beta on the faces of U, V, W; zero on closed
  !  faces), the diagonals and the tridiagonal factors of the line smoother; call whenever the coefficients
  !  change. Device-resident: bu, bv, bw must be present on the device in a GPU build.
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
    If ( mg_direct ) Call mg_direct_build

  End Subroutine mg_set_coef


  !> Level-1 conductances beta*area/distance on the faces between unknown cells; physical Neumann faces carry
  !  nothing
  Subroutine mg_fine_coef(ni, nj, nk, bu, bv, bw, cx, cy, cz)

    Integer(Int32), Intent(In) :: ni, nj, nk
    Real(Int64), Intent(In) :: bu(nx,nyg,nzg), bv(nxg,ny,nzg), bw(nxg,nyg,nz)
    Real(mgk), Intent(Out) :: cx(0:ni+1,0:nj+1,0:nk+1), cy(0:ni+1,0:nj+1,0:nk+1), cz(0:ni+1,0:nj+1,0:nk+1)
    Integer(Int32) :: i, j, k, partner
    Logical :: is_first, is_last

    !$acc kernels present(cx,cy,cz) async(1)
    cx = 0d0;  cy = 0d0;  cz = 0d0
    !$acc end kernels
    !$acc parallel loop collapse(3) present(bu,cx,mg_hy1,mg_hz1) async(1)
    Do k = 1, nk
       Do j = 1, nj
          Do i = 0, ni
             cx(i,j,k) = bu(i+1,j+1,k+1)*mg_hy1(j)*mg_hz1(k)/dx
          End Do
       End Do
    End Do
    !$acc parallel loop collapse(3) present(bv,cy,mg_hz1,yg) async(1)
    Do k = 1, nk
       Do j = 1, nj-1
          Do i = 1, ni
             cy(i,j,k) = bv(i+1,j+1,k+1)*dx*mg_hz1(k)/( yg(j+2) - yg(j+1) )
          End Do
       End Do
    End Do
    !$acc parallel loop collapse(3) present(bw,cz,mg_hy1,zg) async(1)
    Do k = 0, nk
       Do j = 1, nj
          Do i = 1, ni
             cz(i,j,k) = bw(i+1,j+1,k+1)*dx*mg_hy1(j)/( zg(k+2) - zg(k+1) )
          End Do
       End Do
    End Do
    If ( x_bc_type == 1 ) Then
       Call x_periodic_partner(is_first, is_last, partner)
       !$acc kernels present(cx) async(1)
       If ( is_first ) cx(0,:,:) = 0d0
       !$acc end kernels
    End If
    If ( z_bc_type == 1 ) Then
       Call z_periodic_partner(is_first, is_last, partner)
       !$acc kernels present(cz) async(1)
       If ( is_first ) cz(:,:,0) = 0d0
       If ( is_last ) cz(:,:,nk) = 0d0
       !$acc end kernels
    End If

  End Subroutine mg_fine_coef


  !> Coarse conductances: sum of the fine conductances through a coarse face, times the centre-distance ratio
  !  fr (and mg_cscale)
  Subroutine mg_coarsen_coef(nif, njf, nkf, cxf, cyf, czf, nic, njc, nkc, sx, sy, sz, fxr, fyr, fzr, cxc, cyc, czc)

    Integer(Int32), Intent(In) :: nif, njf, nkf, nic, njc, nkc, sx, sy, sz
    Real(mgk), Intent(In) :: cxf(0:nif+1,0:njf+1,0:nkf+1), cyf(0:nif+1,0:njf+1,0:nkf+1), czf(0:nif+1,0:njf+1,0:nkf+1)
    Real(Int64), Intent(In) :: fxr(0:nic), fyr(0:njc), fzr(0:nkc)
    Real(mgk), Intent(Out) :: cxc(0:nic+1,0:njc+1,0:nkc+1), cyc(0:nic+1,0:njc+1,0:nkc+1), czc(0:nic+1,0:njc+1,0:nkc+1)
    Integer(Int32) :: I, J, K, ii, jj, kk, lc
    Real(Int64) :: s

    !$acc kernels present(cxc,cyc,czc) async(1)
    cxc = 0d0;  cyc = 0d0;  czc = 0d0
    !$acc end kernels
    !$acc parallel loop collapse(3) private(lc,s) present(cxf,cxc,fxr) async(1)
    Do K = 1, nkc
       Do J = 1, njc
          Do I = 0, nic
             lc = Min(sx*I, nif);  s = 0d0
             !$acc loop seq
             Do kk = sz*K - (sz-1), Min(sz*K, nkf)
                !$acc loop seq
                Do jj = sy*J - (sy-1), Min(sy*J, njf)
                   s = s + cxf(lc,jj,kk)
                End Do
             End Do
             cxc(I,J,K) = mg_cscale*fxr(I)*s
          End Do
       End Do
    End Do
    !$acc parallel loop collapse(3) private(lc,s) present(cyf,cyc,fyr) async(1)
    Do K = 1, nkc
       Do J = 0, njc
          Do I = 1, nic
             lc = Min(sy*J, njf);  s = 0d0
             !$acc loop seq
             Do kk = sz*K - (sz-1), Min(sz*K, nkf)
                !$acc loop seq
                Do ii = sx*I - (sx-1), Min(sx*I, nif)
                   s = s + cyf(ii,lc,kk)
                End Do
             End Do
             cyc(I,J,K) = mg_cscale*fyr(J)*s
          End Do
       End Do
    End Do
    !$acc parallel loop collapse(3) private(lc,s) present(czf,czc,fzr) async(1)
    Do K = 0, nkc
       Do J = 1, njc
          Do I = 1, nic
             lc = Min(sz*K, nkf);  s = 0d0
             !$acc loop seq
             Do jj = sy*J - (sy-1), Min(sy*J, njf)
                !$acc loop seq
                Do ii = sx*I - (sx-1), Min(sx*I, nif)
                   s = s + czf(ii,jj,lc)
                End Do
             End Do
             czc(I,J,K) = mg_cscale*fzr(K)*s
          End Do
       End Do
    End Do

  End Subroutine mg_coarsen_coef


  !> Diagonal of S and the factors of the y-line tridiagonal T = (1+eps) diag - offdiag_y (Dirichlet
  !  reflection of the x outlet folded into the diagonal); a zero diagonal (solid cell) gives a zero pivot
  !  inverse, so the cell is left alone
  Subroutine mg_factor(ni, nj, nk, cx, cy, cz, dd, m, iw)

    Integer(Int32), Intent(In) :: ni, nj, nk
    Real(mgk), Intent(In) :: cx(0:ni+1,0:nj+1,0:nk+1), cy(0:ni+1,0:nj+1,0:nk+1), cz(0:ni+1,0:nj+1,0:nk+1)
    Real(mgk), Intent(Out) :: dd(0:ni+1,0:nj+1,0:nk+1), m(0:ni+1,0:nj+1,0:nk+1), iw(0:ni+1,0:nj+1,0:nk+1)
    Integer(Int32) :: i, j, k
    Real(Int64) :: dt, w, mm
    Logical :: dirichlet_hi

    dirichlet_hi = ( x_bc_type == 1 .And. mg_xup == MPI_PROC_NULL )
    !$acc kernels present(dd,m,iw) async(1)
    dd = 0d0;  m = 0d0;  iw = 0d0
    !$acc end kernels
    !$acc parallel loop collapse(3) present(cx,cy,cz,dd) async(1)
    Do k = 1, nk
       Do j = 1, nj
          Do i = 1, ni
             dd(i,j,k) = 0d0 + cx(i,j,k) + cx(i-1,j,k) + cy(i,j,k) + cy(i,j-1,k) + cz(i,j,k) + cz(i,j,k-1)
          End Do
       End Do
    End Do
    !$acc parallel loop collapse(2) private(w,dt,mm) present(cx,cy,dd,m,iw) async(1)
    Do k = 1, nk
       Do i = 1, ni
          w = 0d0
          !$acc loop seq
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


  !> Ghost layers of a level array: neighbour ranks and periodic partners by MPI, Neumann copies at the y
  !  walls and at physical x/z boundaries, the reflection mg_xhi_gs at a Dirichlet x outlet
  Subroutine mg_halo(a, ni, nj, nk)

    Integer(Int32), Intent(In) :: ni, nj, nk
    Real(mgk), Intent(InOut) :: a(0:ni+1,0:nj+1,0:nk+1)
    Integer(Int32) :: i, j, k, nx_, nz_, nreq, req(8)
    Logical :: xmpi, zmpi
    Real(Int64) :: glo, ghi

    !$acc kernels present(a) async(1)
    a(:,0,:) = a(:,1,:)
    a(:,nj+1,:) = a(:,nj,:)
    !$acc end kernels
    xmpi = .Not. ( mg_xself .Or. ( mg_xdn == MPI_PROC_NULL .And. mg_xup == MPI_PROC_NULL ) )
    zmpi = .Not. ( mg_zself .Or. ( mg_zdn == MPI_PROC_NULL .And. mg_zup == MPI_PROC_NULL ) )
    glo = mg_xlo_gs;  ghi = mg_xhi_gs
    If ( mg_xself ) Then
       !$acc kernels present(a) async(1)
       a(0,:,:) = a(ni,:,:)
       a(ni+1,:,:) = a(1,:,:)
       !$acc end kernels
    Else If ( .Not. xmpi ) Then
       ! no x neighbour: fill on the device, no host staging of the send/receive planes
       !$acc kernels present(a) async(1)
       a(0,:,:) = glo*a(1,:,:)
       a(ni+1,:,:) = ghi*a(ni,:,:)
       !$acc end kernels
    End If
    If ( mg_zself ) Then
       !$acc kernels present(a) async(1)
       a(:,:,0) = a(:,:,nk)
       a(:,:,nk+1) = a(:,:,1)
       !$acc end kernels
    Else If ( .Not. zmpi ) Then
       !$acc kernels present(a) async(1)
       a(:,:,0) = a(:,:,1)
       a(:,:,nk+1) = a(:,:,nk)
       !$acc end kernels
    End If
    If ( .Not. ( xmpi .Or. zmpi ) ) Return

    ! both directions in one exchange (the stencil does not read the corners) with a single host synchronisation:
    ! every wait on a GPU shared by several ranks costs far more than the copy
    nx_ = (nj+2)*(nk+2);  nz_ = (ni+2)*(nj+2)
    If ( xmpi ) Then
       !$acc parallel loop collapse(2) present(a,mg_sbx) async(1)
       Do k = 0, nk+1
          Do j = 0, nj+1
             mg_sbx(1+j+(nj+2)*k,1) = a(ni,j,k)
             mg_sbx(1+j+(nj+2)*k,2) = a(1,j,k)
          End Do
       End Do
       !$acc update host(mg_sbx(1:nx_,1:2)) async(1)
    End If
    If ( zmpi ) Then
       !$acc parallel loop collapse(2) present(a,mg_sbz) async(1)
       Do j = 0, nj+1
          Do i = 0, ni+1
             mg_sbz(1+i+(ni+2)*j,1) = a(i,j,nk)
             mg_sbz(1+i+(ni+2)*j,2) = a(i,j,1)
          End Do
       End Do
       !$acc update host(mg_sbz(1:nz_,1:2)) async(1)
    End If
    !$acc wait(1)
    nreq = 0
    If ( xmpi ) Then
       Call MPI_Irecv(mg_rbx(1,1), nx_, mg_mpi_t, mg_xdn, 301, MPI_COMM_WORLD, req(nreq+1), ierr)
       Call MPI_Irecv(mg_rbx(1,2), nx_, mg_mpi_t, mg_xup, 302, MPI_COMM_WORLD, req(nreq+2), ierr)
       Call MPI_Isend(mg_sbx(1,1), nx_, mg_mpi_t, mg_xup, 301, MPI_COMM_WORLD, req(nreq+3), ierr)
       Call MPI_Isend(mg_sbx(1,2), nx_, mg_mpi_t, mg_xdn, 302, MPI_COMM_WORLD, req(nreq+4), ierr)
       nreq = nreq + 4
    End If
    If ( zmpi ) Then
       Call MPI_Irecv(mg_rbz(1,1), nz_, mg_mpi_t, mg_zdn, 303, MPI_COMM_WORLD, req(nreq+1), ierr)
       Call MPI_Irecv(mg_rbz(1,2), nz_, mg_mpi_t, mg_zup, 304, MPI_COMM_WORLD, req(nreq+2), ierr)
       Call MPI_Isend(mg_sbz(1,1), nz_, mg_mpi_t, mg_zup, 303, MPI_COMM_WORLD, req(nreq+3), ierr)
       Call MPI_Isend(mg_sbz(1,2), nz_, mg_mpi_t, mg_zdn, 304, MPI_COMM_WORLD, req(nreq+4), ierr)
       nreq = nreq + 4
    End If
    Call MPI_Waitall(nreq, req, MPI_STATUSES_IGNORE, ierr)
    If ( xmpi ) Then
       !$acc update device(mg_rbx(1:nx_,1:2)) async(1)
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
    If ( zmpi ) Then
       !$acc update device(mg_rbz(1:nz_,1:2)) async(1)
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
    Real(mgk), Intent(In) :: x(0:ni+1,0:nj+1,0:nk+1), b(0:ni+1,0:nj+1,0:nk+1), dd(0:ni+1,0:nj+1,0:nk+1)
    Real(mgk), Intent(In) :: cx(0:ni+1,0:nj+1,0:nk+1), cy(0:ni+1,0:nj+1,0:nk+1), cz(0:ni+1,0:nj+1,0:nk+1)
    Real(mgk), Intent(InOut) :: t(0:ni+1,0:nj+1,0:nk+1)
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


#ifndef GPU_POISSON
  !> t <- T^-1 t: one tridiagonal solve along y per column
  Subroutine mg_line(ni, nj, nk, t, cy, m, iw)

    Integer(Int32), Intent(In) :: ni, nj, nk
    Real(mgk), Intent(InOut) :: t(0:ni+1,0:nj+1,0:nk+1)
    Real(mgk), Intent(In) :: cy(0:ni+1,0:nj+1,0:nk+1), m(0:ni+1,0:nj+1,0:nk+1), iw(0:ni+1,0:nj+1,0:nk+1)
    Integer(Int32) :: i, j, k

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

  End Subroutine mg_line


  !> d = c1 d + c2 t; x = x + d (x = d when zero)
  Subroutine mg_upd(ni, nj, nk, x, d, t, c1, c2, zero)

    Integer(Int32), Intent(In) :: ni, nj, nk
    Real(mgk), Intent(InOut) :: x(0:ni+1,0:nj+1,0:nk+1), d(0:ni+1,0:nj+1,0:nk+1)
    Real(mgk), Intent(In) :: t(0:ni+1,0:nj+1,0:nk+1)
    Real(Int64), Intent(In) :: c1, c2
    Logical, Intent(In) :: zero
    Integer(Int32) :: i, j, k
    Real(mgk) :: c1s, c2s

    c1s = Real(c1, mgk);  c2s = Real(c2, mgk)
    Do k = 1, nk
       Do j = 1, nj
          Do i = 1, ni
             d(i,j,k) = c1s*d(i,j,k) + c2s*t(i,j,k)
             If ( zero ) Then
                x(i,j,k) = d(i,j,k)
             Else
                x(i,j,k) = x(i,j,k) + d(i,j,k)
             End If
          End Do
       End Do
    End Do

  End Subroutine mg_upd
#endif

#ifdef GPU_POISSON
  !> Forward elimination of the y-line solve applied to the residual b - S x (to b when zero) in one pass
  Subroutine mg_rfwd(ni, nj, nk, x, b, cx, cy, cz, dd, m, t, zero)

    Integer(Int32), Intent(In) :: ni, nj, nk
    Real(mgk), Intent(In) :: x(0:ni+1,0:nj+1,0:nk+1), b(0:ni+1,0:nj+1,0:nk+1), dd(0:ni+1,0:nj+1,0:nk+1)
    Real(mgk), Intent(In) :: cx(0:ni+1,0:nj+1,0:nk+1), cy(0:ni+1,0:nj+1,0:nk+1), cz(0:ni+1,0:nj+1,0:nk+1)
    Real(mgk), Intent(In) :: m(0:ni+1,0:nj+1,0:nk+1)
    Real(mgk), Intent(InOut) :: t(0:ni+1,0:nj+1,0:nk+1)
    Logical, Intent(In) :: zero
    Integer(Int32) :: i, j, k
    Real(mgk) :: r, tp

    !$acc parallel loop collapse(2) gang vector private(r,tp) present(x,b,cx,cy,cz,dd,m,t) async(1)
    Do k = 1, nk
       Do i = 1, ni
          tp = 0.0
          !$acc loop seq
          !dir$ unroll 4
          Do j = 1, nj
             If ( zero ) Then
                r = b(i,j,k)
             Else
                r = b(i,j,k) - ( dd(i,j,k)*x(i,j,k) - cx(i,j,k)*x(i+1,j,k) - cx(i-1,j,k)*x(i-1,j,k) &
                               - cy(i,j,k)*x(i,j+1,k) - cy(i,j-1,k)*x(i,j-1,k) &
                               - cz(i,j,k)*x(i,j,k+1) - cz(i,j,k-1)*x(i,j,k-1) )
             End If
             r = r + m(i,j,k)*tp
             t(i,j,k) = r
             tp = r
          End Do
       End Do
    End Do

  End Subroutine mg_rfwd


  !> Back substitution of the y-line solve and the Chebyshev update d = c1 d + c2 t, x = x + d (x = d when
  !  zero) in one pass
  Subroutine mg_bwdupd(ni, nj, nk, x, d, t, cy, iw, c1, c2, zero)

    Integer(Int32), Intent(In) :: ni, nj, nk
    Real(mgk), Intent(InOut) :: x(0:ni+1,0:nj+1,0:nk+1), d(0:ni+1,0:nj+1,0:nk+1)
    Real(mgk), Intent(In) :: t(0:ni+1,0:nj+1,0:nk+1), cy(0:ni+1,0:nj+1,0:nk+1), iw(0:ni+1,0:nj+1,0:nk+1)
    Real(Int64), Intent(In) :: c1, c2
    Logical, Intent(In) :: zero
    Integer(Int32) :: i, j, k
    Real(mgk) :: c1s, c2s, tn, dn

    c1s = Real(c1, mgk);  c2s = Real(c2, mgk)
    !$acc parallel loop collapse(2) gang vector private(tn,dn) present(x,d,t,cy,iw) async(1)
    Do k = 1, nk
       Do i = 1, ni
          tn = 0.0
          !$acc loop seq
          !dir$ unroll 4
          Do j = nj, 1, -1
             If ( j == nj ) Then
                tn = t(i,j,k)*iw(i,j,k)
             Else
                tn = ( t(i,j,k) + cy(i,j,k)*tn )*iw(i,j,k)
             End If
             dn = c1s*d(i,j,k) + c2s*tn
             d(i,j,k) = dn
             If ( zero ) Then
                x(i,j,k) = dn
             Else
                x(i,j,k) = x(i,j,k) + dn
             End If
          End Do
       End Do
    End Do

  End Subroutine mg_bwdupd
#endif


  !> Chebyshev polynomial of degree deg in T^-1 S on [lo, hi] (hi = 2 bounds the line-Jacobi spectrum), from x
  !  = 0 when zero
  Subroutine mg_cheb(ni, nj, nk, x, b, t, d, cx, cy, cz, dd, m, iw, deg, lo, hi, zero)

    Integer(Int32), Intent(In) :: ni, nj, nk, deg
    Real(mgk), Intent(InOut) :: x(0:ni+1,0:nj+1,0:nk+1), t(0:ni+1,0:nj+1,0:nk+1), d(0:ni+1,0:nj+1,0:nk+1)
    Real(mgk), Intent(In) :: b(0:ni+1,0:nj+1,0:nk+1), cx(0:ni+1,0:nj+1,0:nk+1), cy(0:ni+1,0:nj+1,0:nk+1)
    Real(mgk), Intent(In) :: cz(0:ni+1,0:nj+1,0:nk+1), dd(0:ni+1,0:nj+1,0:nk+1), m(0:ni+1,0:nj+1,0:nk+1)
    Real(mgk), Intent(In) :: iw(0:ni+1,0:nj+1,0:nk+1)
    Real(Int64), Intent(In) :: lo, hi
    Logical, Intent(In) :: zero
    Integer(Int32) :: it
    Logical :: first
    Real(Int64) :: theta, delta, sigma, rho0, rho1, c1, c2

    theta = 0.5d0*(hi + lo);  delta = 0.5d0*(hi - lo);  sigma = theta/delta;  rho0 = 1d0/sigma
    Do it = 1, deg
       first = ( it == 1 .And. zero )
       If ( .Not. first ) Call mg_halo(x, ni, nj, nk)
       If ( it == 1 ) Then
          c1 = 0d0;  c2 = 1d0/theta
       Else
          rho1 = 1d0/(2d0*sigma - rho0)
          c1 = rho1*rho0;  c2 = 2d0*rho1/delta
          rho0 = rho1
       End If
#ifdef GPU_POISSON
       Call mg_rfwd(ni, nj, nk, x, b, cx, cy, cz, dd, m, t, first)
       Call mg_bwdupd(ni, nj, nk, x, d, t, cy, iw, c1, c2, first)
#else
       If ( first ) Then
          Call mg_copy(ni, nj, nk, b, t)
       Else
          Call mg_resid(ni, nj, nk, x, b, cx, cy, cz, dd, t)
       End If
       Call mg_line(ni, nj, nk, t, cy, m, iw)
       Call mg_upd(ni, nj, nk, x, d, t, c1, c2, first)
#endif
    End Do

  End Subroutine mg_cheb


  Subroutine mg_copy(ni, nj, nk, a, c)

    Integer(Int32), Intent(In) :: ni, nj, nk
    Real(mgk), Intent(In) :: a(0:ni+1,0:nj+1,0:nk+1)
    Real(mgk), Intent(InOut) :: c(0:ni+1,0:nj+1,0:nk+1)
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
    Real(mgk), Intent(In) :: tf(0:nif+1,0:njf+1,0:nkf+1)
    Real(mgk), Intent(InOut) :: bc(0:nic+1,0:njc+1,0:nkc+1)
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
    Real(mgk), Intent(InOut) :: xf(0:nif+1,0:njf+1,0:nkf+1)
    Real(mgk), Intent(In) :: xc(0:nic+1,0:njc+1,0:nkc+1)
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


  !> One symmetric cycle from a zero guess on the level-1 right-hand side mg_b (V-cycle for mg_cycle = 1,
  !  W-cycle for 2)
  Subroutine mg_vcycle

    Call mg_level(1)

  End Subroutine mg_vcycle


  !> Cycle on level l from a zero guess: smooth, restrict the residual, mg_cycle visits of the coarse level
  !  (each one on the residual left by the previous ones, their corrections added), prolong, smooth. Every
  !  visit is a fixed symmetric operator, so the result is too
  Recursive Subroutine mg_level(l)

    Integer(Int32), Intent(In) :: l
    Integer(Int32) :: o, oc, g, ni, nj, nk
    Real(Int64), Parameter :: hi = 2d0

    o = mg_off(l) + 1
    If ( l == mg_nlev ) Then
       If ( mg_direct ) Then
          Call mg_direct_solve(mg_ni(l), mg_nj(l), mg_nk(l), mg_b(o), mg_x(o))
       Else
          Call mg_cheb(mg_ni(l), mg_nj(l), mg_nk(l), mg_x(o), mg_b(o), mg_t(o), mg_d(o), mg_cx(o), mg_cy(o), mg_cz(o), &
                       mg_dd(o), mg_m(o), mg_iw(o), mg_deg_coarse, hi/mg_ratio_coarse, hi, .True.)
       End If
       Return
    End If
    oc = mg_off(l+1) + 1
    ni = mg_ni(l+1);  nj = mg_nj(l+1);  nk = mg_nk(l+1)
    Call mg_cheb(mg_ni(l), mg_nj(l), mg_nk(l), mg_x(o), mg_b(o), mg_t(o), mg_d(o), mg_cx(o), mg_cy(o), mg_cz(o), &
                 mg_dd(o), mg_m(o), mg_iw(o), mg_deg, hi/mg_ratio, hi, .True.)
    Call mg_halo(mg_x(o), mg_ni(l), mg_nj(l), mg_nk(l))
    Call mg_resid(mg_ni(l), mg_nj(l), mg_nk(l), mg_x(o), mg_b(o), mg_cx(o), mg_cy(o), mg_cz(o), mg_dd(o), mg_t(o))
    Call mg_restrict(mg_ni(l), mg_nj(l), mg_nk(l), mg_t(o), ni, nj, nk, mg_sx(l), mg_sy(l), mg_sz(l), mg_b(oc))
    Do g = 1, mg_cycle
       Call mg_level(l+1)
       If ( mg_cycle == 1 ) Exit
       If ( g == 1 ) Then
          Call mg_copy(ni, nj, nk, mg_x(oc), mg_s(oc))
       Else
          Call mg_add(ni, nj, nk, mg_x(oc), mg_s(oc))
       End If
       If ( g == mg_cycle ) Exit
       Call mg_halo(mg_x(oc), ni, nj, nk)
       Call mg_resid(ni, nj, nk, mg_x(oc), mg_b(oc), mg_cx(oc), mg_cy(oc), mg_cz(oc), mg_dd(oc), mg_t(oc))
       Call mg_copy(ni, nj, nk, mg_t(oc), mg_b(oc))
    End Do
    If ( mg_cycle > 1 ) Call mg_copy(ni, nj, nk, mg_s(oc), mg_x(oc))
    Call mg_prolong(mg_ni(l), mg_nj(l), mg_nk(l), mg_x(o), ni, nj, nk, mg_sx(l), mg_sy(l), mg_sz(l), mg_x(oc))
    Call mg_cheb(mg_ni(l), mg_nj(l), mg_nk(l), mg_x(o), mg_b(o), mg_t(o), mg_d(o), mg_cx(o), mg_cy(o), mg_cz(o), &
                 mg_dd(o), mg_m(o), mg_iw(o), mg_deg, hi/mg_ratio, hi, .False.)

  End Subroutine mg_level


  Subroutine mg_add(ni, nj, nk, a, c)

    Integer(Int32), Intent(In) :: ni, nj, nk
    Real(mgk), Intent(In) :: a(0:ni+1,0:nj+1,0:nk+1)
    Real(mgk), Intent(InOut) :: c(0:ni+1,0:nj+1,0:nk+1)
    Integer(Int32) :: i, j, k

    !$acc parallel loop collapse(3) present(a,c) async(1)
    Do k = 1, nk
       Do j = 1, nj
          Do i = 1, ni
             c(i,j,k) = c(i,j,k) + a(i,j,k)
          End Do
       End Do
    End Do

  End Subroutine mg_add


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
    Real(mgk), Intent(InOut) :: b(0:ni+1,0:nj+1,0:nk+1)
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
    Real(mgk), Intent(In) :: x(0:ni+1,0:nj+1,0:nk+1)
    Real(Int64), Intent(InOut) :: z(nxg,nyg,nzg)
    Integer(Int32) :: i, j, k

    ! only the cells outside the unknowns (ghost layers and the periodic duplicate) are not overwritten below
    !$acc kernels present(z) async(1)
    z(1,:,:) = 0d0;  z(nxg,:,:) = 0d0;  z(:,1,:) = 0d0;  z(:,nyg,:) = 0d0;  z(:,:,1) = 0d0;  z(:,:,nzg) = 0d0
    z(ni+2:nxg,:,:) = 0d0;  z(:,nj+2:nyg,:) = 0d0;  z(:,:,nk+2:nzg) = 0d0
    !$acc end kernels
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
