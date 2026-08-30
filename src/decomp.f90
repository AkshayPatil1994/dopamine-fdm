!> 2decomp&fft pencil decomposition bookkeeping: the main (face-count) pencil grid, plus a dedicated Poisson pencil sized to the actual Fourier-transform grid
Module decomp

  Use iso_fortran_env, Only : Int32, Int64
  Use global,   Only : nx_global, ny_global, nz_global, p_row, p_col, x_bc_type
  Use mpi,      Only : nprocs, myid, MPI_PROC_NULL, k1_global, k2_global, kg1_global, kg2_global, &
                        i1_global, i2_global, ig1_global, ig2_global
  Use decomp_2d, Only : decomp_2d_init, decomp_2d_finalize, decomp_info, decomp_info_init, &
       d2d_xstart => xstart, d2d_xend => xend, d2d_xsize => xsize, &
       d2d_ystart => ystart, d2d_yend => yend, d2d_ysize => ysize, &
       d2d_zstart => zstart, d2d_zend => zend, d2d_zsize => zsize, &
       transpose_x_to_y, transpose_y_to_x, transpose_y_to_z, transpose_z_to_y
  Use m_decomp_pool, Only : decomp_pool
  Use decomp_2d_constants, Only : complex_type

  ! prevent implicit typing
  Implicit None

  ! dedicated pencil grid for the pressure Poisson solve, sized to the
  ! Fourier-transform grid (nxp_global x nym_global x nzp_global), which
  ! differs from the main (face-count) pencil grid above
  Type(decomp_info) :: decomp_poisson

Contains

  !> Initialise the main 2decomp&fft pencil grid; nx/ny/nz here are face-point counts, matching this code's existing (not cell-count) convention
  Subroutine decomp_init_pencil

    If ( p_row == 0 .And. p_col == 0 ) Call decomp_auto_factorize(p_row, p_col)

    ! complex_pool=.true.: the Poisson solve's transpose_* calls move complex
    ! data, and need the shared 2decomp&fft memory pool to have a complex-type
    ! shape registered (real-only by default) or they abort at run time.
    Call decomp_2d_init(nx_global, ny_global, nz_global, p_row, p_col, complex_pool=.true.)

  End Subroutine decomp_init_pencil

  !> Auto-pick p_row/p_col when both are left 0: prefer pure z-split (p_row=1) since it skips x-transpose communication entirely; only fall back to a true 2D split if z-only can't give every rank >=2 interior z-cells, then pick the nprocs factor pair closest to the grid's x:z aspect ratio among pairs leaving >=2 interior cells per rank in both directions
  Subroutine decomp_auto_factorize(p_row_out, p_col_out)

    Integer(Int32), Intent(Out) :: p_row_out, p_col_out

    Integer(Int32) :: nxm_g, nzm_g, row_cap, col_cap, trial_row, trial_col, best_row, best_col
    Real   (Int64) :: target_ratio, score, best_score
    Logical        :: found

    nxm_g = nx_global - 1
    nzm_g = nz_global - 1
    ! 2decomp&fft's other pencil views (x-pencil, z-pencil) also split y, so p_row<=min(nx,ny) and p_col<=min(ny,nz) is a hard library requirement, not just this code's own >=2-interior-cell floor
    row_cap = Min(nx_global, ny_global)
    col_cap = Min(ny_global, nz_global)

    If ( nzm_g/nprocs >= 2 .And. nprocs <= col_cap ) Then
       p_row_out = 1
       p_col_out = nprocs
       Return
    End If

    target_ratio = Real(nxm_g, Int64) / Real(nzm_g, Int64)
    best_score = Huge(best_score)
    found = .False.
    Do trial_row = 1, nprocs
       If ( Mod(nprocs, trial_row) /= 0 ) Cycle
       trial_col = nprocs / trial_row
       If ( nxm_g/trial_row < 2 .Or. nzm_g/trial_col < 2 ) Cycle
       If ( trial_row > row_cap .Or. trial_col > col_cap ) Cycle
       score = Abs( Log( Real(trial_row,Int64)/Real(trial_col,Int64) / target_ratio ) )
       If ( score < best_score ) Then
          best_score = score
          best_row = trial_row
          best_col = trial_col
          found = .True.
       End If
    End Do

    If ( .Not. found ) Stop 'ERROR: no valid p_row/p_col auto-split for nprocs; set explicitly'

    p_row_out = best_row
    p_col_out = best_col

  End Subroutine decomp_auto_factorize

  !> Register the Poisson-solve pencil grid (nxp_g x nyp_g x nzp_g: the Fourier-transform grid, not the face-point grid above)
  Subroutine decomp_init_poisson_pencil(nxp_g, nyp_g, nzp_g)

    Integer(Int32), Intent(In) :: nxp_g, nyp_g, nzp_g

    Call decomp_info_init(nxp_g, nyp_g, nzp_g, decomp_poisson)
    ! register decomp_poisson's complex shape with the shared memory pool too (see note above)
    Call decomp_pool%new_shape(complex_type, decomp_poisson)

  End Subroutine decomp_init_poisson_pencil

  !> Build flat-rank (0:nprocs-1) x/z ownership arrays: i1_global/i2_global/
  !> ig1_global/ig2_global (x-split via p_row) and k1_global/k2_global/
  !> kg1_global/kg2_global (z-split via p_col). Assumes MPI's standard
  !> row-major Cartesian rank ordering (rank = row*p_col + col, the same
  !> convention MPI_CART_CREATE uses with reorder=.false.) -- self-checked
  !> against this rank's own decomp_main partition before trusting it.
  Subroutine decomp_build_xz_ranges

    Integer(Int32) :: r, row, col
    Integer(Int32), Allocatable :: xst(:), xen(:), zst(:), zen(:)
    Integer(Int32), Allocatable :: xst_chk(:), xen_chk(:)
    Integer(Int32) :: chk_i1, chk_i2, chk_j1, chk_j2, chk_k1, chk_k2

    ! Self-check the row-major rank<->(row,col) mapping assumption against
    ! this rank's own decomp_main partition (same algorithm, same nx_global
    ! input -- an apples-to-apples check of the mapping alone) before
    ! trusting it to build every rank's x/z ownership below.
    Allocate ( xst_chk(0:p_row-1), xen_chk(0:p_row-1) )
    Call distribute_1d( nx_global, p_row, xst_chk, xen_chk )
    Call y_pencil_local_range(chk_i1, chk_i2, chk_j1, chk_j2, chk_k1, chk_k2)
    row = myid / p_col
    If ( xst_chk(row) /= chk_i1 .Or. xen_chk(row) /= chk_i2 ) &
         Stop 'ERROR: rank<->(row,col) mapping assumption wrong -- cannot safely build x/z ownership'
    Deallocate ( xst_chk, xen_chk )

    Allocate ( xst(0:p_row-1), xen(0:p_row-1) )
    Allocate ( zst(0:p_col-1), zen(0:p_col-1) )
    ! Distribute decomp_poisson's own (periodic-reduced) point count, not nx_global-1/nz_global-1, so the two independent "remainder on last ranks" splits can't disagree on a middle rank and silently truncate its rhs_p write-back; z is always periodic, x only for x_bc_type==0 (DCT-IV keeps the full count); the dropped point is restored on the true last rank/column below
    If ( x_bc_type == 0 ) Then
       Call distribute_1d( nx_global-2, p_row, xst, xen )
       xen(p_row-1) = xen(p_row-1) + 1
    Else
       Call distribute_1d( nx_global-1, p_row, xst, xen )
    End If
    Call distribute_1d( nz_global-2, p_col, zst, zen )
    zen(p_col-1) = zen(p_col-1) + 1

    Allocate (  i1_global(0:nprocs-1),  i2_global(0:nprocs-1) )
    Allocate ( ig1_global(0:nprocs-1), ig2_global(0:nprocs-1) )
    Allocate (  k1_global(0:nprocs-1),  k2_global(0:nprocs-1) )
    Allocate ( kg1_global(0:nprocs-1), kg2_global(0:nprocs-1) )

    Do r = 0, nprocs-1
       row = r / p_col
       col = Mod(r, p_col)

       ig1_global(r) = xst(row)
       ig2_global(r) = xen(row) + 2
       i1_global(r)  = ig1_global(r)
       i2_global(r)  = ig2_global(r)
       If ( row == p_row-1 ) i2_global(r) = nx_global

       kg1_global(r) = zst(col)
       kg2_global(r) = zen(col) + 2
       k1_global(r)  = kg1_global(r)
       k2_global(r)  = kg2_global(r)
       If ( col == p_col-1 ) k2_global(r) = nz_global
    End Do

    Deallocate ( xst, xen, zst, zen )

  End Subroutine decomp_build_xz_ranges

  !> z-halo neighbor ranks (same row, adjacent column); MPI_PROC_NULL at the row's column edges (interior halo exchange -- not periodic wraparound, see z_periodic_partner)
  Subroutine z_halo_neighbors(up, down)

    Integer(Int32), Intent(Out) :: up, down
    Integer(Int32) :: row, col

    row = myid / p_col
    col = Mod(myid, p_col)

    If ( col == p_col-1 ) Then
       up = MPI_PROC_NULL
    Else
       up = myid + 1
    End If

    If ( col == 0 ) Then
       down = MPI_PROC_NULL
    Else
       down = myid - 1
    End If

  End Subroutine z_halo_neighbors

  !> x-halo neighbor ranks (same column, adjacent row); MPI_PROC_NULL at the column's row edges
  Subroutine x_halo_neighbors(up, down)

    Integer(Int32), Intent(Out) :: up, down
    Integer(Int32) :: row, col

    row = myid / p_col
    col = Mod(myid, p_col)

    If ( row == p_row-1 ) Then
       up = MPI_PROC_NULL
    Else
       up = myid + p_col
    End If

    If ( row == 0 ) Then
       down = MPI_PROC_NULL
    Else
       down = myid - p_col
    End If

  End Subroutine x_halo_neighbors

  !> This rank's periodic-wraparound partner in z (only meaningful if is_first or is_last)
  Subroutine z_periodic_partner(is_first, is_last, partner)

    Logical,        Intent(Out) :: is_first, is_last
    Integer(Int32), Intent(Out) :: partner
    Integer(Int32) :: row, col

    row = myid / p_col
    col = Mod(myid, p_col)
    is_first = ( col == 0 )
    is_last  = ( col == p_col-1 )
    partner  = -1
    If ( is_first ) partner = row*p_col + (p_col-1)
    If ( is_last  ) partner = row*p_col

  End Subroutine z_periodic_partner

  !> This rank's periodic-wraparound partner in x (only meaningful if is_first or is_last)
  Subroutine x_periodic_partner(is_first, is_last, partner)

    Logical,        Intent(Out) :: is_first, is_last
    Integer(Int32), Intent(Out) :: partner
    Integer(Int32) :: row, col

    row = myid / p_col
    col = Mod(myid, p_col)
    is_first = ( row == 0 )
    is_last  = ( row == p_row-1 )
    partner  = -1
    If ( is_first ) partner = (p_row-1)*p_col + col
    If ( is_last  ) partner = col

  End Subroutine x_periodic_partner

  !> Split n points across nprocs ranks as evenly as possible (remainder on the last ranks); replicates 2decomp&fft's own distribute() algorithm locally (no MPI/decomp_info needed), so callers aren't restricted to evenly-divisible sizes
  Subroutine distribute_1d(n, nprocs, st, en)

    Integer(Int32), Intent(In)  :: n, nprocs
    Integer(Int32), Intent(Out) :: st(0:nprocs-1), en(0:nprocs-1)

    Integer(Int32) :: i, size1, nu, nl

    size1 = n / nprocs
    nu = n - size1*nprocs
    nl = nprocs - nu
    st(0) = 1
    en(0) = size1
    Do i = 1, nl-1
       st(i) = st(i-1) + size1
       en(i) = en(i-1) + size1
    End Do
    size1 = size1 + 1
    Do i = nl, nprocs-1
       st(i) = en(i-1) + 1
       en(i) = en(i-1) + size1
    End Do

  End Subroutine distribute_1d

  !> Local face-point index range (1-based, inclusive) owned by this rank in the y-pencil, mirroring src/mpi.f90's k1_global/k2_global convention for z
  Subroutine y_pencil_local_range(i1, i2, j1, j2, k1, k2)

    Integer(Int32), Intent(Out) :: i1, i2, j1, j2, k1, k2

    i1 = d2d_ystart(1); i2 = d2d_yend(1)
    j1 = d2d_ystart(2); j2 = d2d_yend(2)
    k1 = d2d_ystart(3); k2 = d2d_yend(3)

  End Subroutine y_pencil_local_range

End Module decomp
