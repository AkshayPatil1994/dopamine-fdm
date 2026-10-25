module raycast_mod
    !
    ! raycast_mod — point-in-solid classification of grid points against a triangle surface mesh.
    !
    ! sign_method = 1 (default):
    !   Tier 1  From every grid point six half-line rays are cast (+-x, +-y, +-z). Along a ray the
    !           ORIENTATION-SIGNED crossings of the solid's triangles are summed (+1 where the surface normal
    !           points along the ray, -1 against it); that sum is the point's winding number (0 outside, 1
    !           inside, 2 where two shells overlap or a patch is doubled) -- provided the ray does not meet a
    !           defect. A clean ray gives the right answer even if the mesh is open elsewhere (a missing top
    !             face
    !           does not matter to a ray going down), so each ray is an independent witness. On an open shell
    !             a
    !           ray that escapes without meeting any triangle abstains (it may simply have left through a
    !             hole).
    !           A point is decided here when at least three rays vote and they are unanimous. (Plain crossing
    !           parity cannot do this: after a missed crossing it is silently wrong for the rest of the line,
    !           and doubled or overlapping shells cancel it.)
    !   Tier 2  Every other point (rays that disagree or that all abstain) gets the generalized winding number
    !           of the solid evaluated with a BVH, thresholded at 1/2; see winding_mod. This degrades
    !           gracefully with holes (w stays above 1/2 until the hole subtends more than half the sphere)
    !           and needs no closed shell.
    !   Orientation, duplicate triangles, degenerate slivers and (planar) hole lids come from meshprep_mod.
    ! sign_method = 0 : the original unsigned 3-axis parity majority vote (needs a watertight shell).
    ! sign_method = 2 : winding number at every point, plus a count of how often tier 1 disagrees
    !   (validation).
    ! A solid flagged `ground` in the manifest is an open sheet classified as "below its topmost surface".
    !
    ! The work is split across MPI ranks: tier-1 lines round-robin over ranks (their votes are combined with
    !   one
    ! all-reduce per solid), tier 2 and the final write-back on each rank's own x-slab.
    !
    use mpi
    use utils_io,          only : dp, solid_is_ground, sign_method
    use spatial_hash_mod,  only : hcells, csr_offset, csr_triid, cell_flat, hash_ix
    use winding_mod,       only : bvh_t, bvh_build, bvh_free, winding_number
    implicit none

    ! Repaired mesh view used by the robust classifier (set once by raycast_set_mesh)
    integer, allocatable, private :: rc_faces(:,:)        ! faces with the orientation repair applied
    logical, allocatable, private :: rc_active(:)         ! false: degenerate or cancelled duplicate triangle
    type(bvh_t), allocatable, private :: rc_bvh(:)        ! one winding-number hierarchy per solid
    logical, allocatable, private :: rc_open(:)           ! solid has boundary / non-manifold edges (not a closed shell)
    real(dp), allocatable, private :: rc_lo(:,:), rc_hi(:,:)   ! per-solid bounding box of the active triangles

contains

    ! -----------------------------------------------------------------------
    ! Gather every triangle ID whose hash bin lies in the 3x3 neighbourhood
    ! (in the two axes perpendicular to `axis`) of (coord_b, coord_c), across
    ! the full hash-cell range along `axis` — i.e. every candidate that could
    ! be crossed by the infinite line through (coord_b, coord_c) running
    ! along `axis`. axis/dim_b/dim_c follow the 1=x, 2=y, 3=z convention.
    subroutine query_hash_line(axis, coord_b, coord_c, candidate_tris, n_candidates)
        implicit none
        integer,  intent(in)    :: axis
        real(dp), intent(in)    :: coord_b, coord_c
        integer,  allocatable, intent(inout) :: candidate_tris(:)
        integer,  intent(out)   :: n_candidates
        integer :: dim_b, dim_c, hb, hc, cb, cc, i_along, cflat, lo, hi, n_tmp, need

        call perp_dims(axis, dim_b, dim_c)
        hb = hash_ix(coord_b, dim_b)
        hc = hash_ix(coord_c, dim_c)

        if (.not. allocated(candidate_tris)) allocate(candidate_tris(256))
        n_tmp = 0
        do cb = max(1, hb-1), min(hcells(dim_b), hb+1)
            do cc = max(1, hc-1), min(hcells(dim_c), hc+1)
                do i_along = 1, hcells(axis)
                    cflat = hcell_flat(axis, i_along, cb, cc)
                    lo = csr_offset(cflat); hi = csr_offset(cflat+1) - 1
                    if (hi >= lo) then
                        need = n_tmp + (hi - lo + 1)
                        if (need > size(candidate_tris)) call grow_int(candidate_tris, need)
                        candidate_tris(n_tmp+1:need) = csr_triid(lo:hi)
                        n_tmp = need
                    end if
                end do
            end do
        end do
        n_candidates = n_tmp
    end subroutine query_hash_line

    pure subroutine perp_dims(axis, dim_b, dim_c)
        integer, intent(in)  :: axis
        integer, intent(out) :: dim_b, dim_c
        select case (axis)
        case (1); dim_b = 2; dim_c = 3
        case (2); dim_b = 1; dim_c = 3
        case default; dim_b = 1; dim_c = 2
        end select
    end subroutine perp_dims

    pure integer function hcell_flat(axis, i_along, jb, jc)
        integer, intent(in) :: axis, i_along, jb, jc
        select case (axis)
        case (1); hcell_flat = cell_flat(i_along, jb, jc)
        case (2); hcell_flat = cell_flat(jb, i_along, jc)
        case default; hcell_flat = cell_flat(jb, jc, i_along)
        end select
    end function hcell_flat

    subroutine grow_int(arr, need)
        integer, allocatable, intent(inout) :: arr(:)
        integer, intent(in) :: need
        integer, allocatable :: tmp(:)
        allocate(tmp(max(need, 2*size(arr))))
        tmp(1:size(arr)) = arr
        call move_alloc(tmp, arr)
    end subroutine grow_int

    ! -----------------------------------------------------------------------
    ! Build the repaired mesh view for the robust classifier: faces with flipped triangles re-wound, the
    !   active mask,
    ! and one winding-number hierarchy per (non-ground) solid. Call once, on every rank, after meshprep's
    !   analysis.
    subroutine raycast_set_mesh(vertices, faces, face_solid_id, nfaces, nsolids, flip, active, is_open)
        implicit none
        real(dp), intent(in) :: vertices(:,:)
        integer,  intent(in) :: faces(:,:), face_solid_id(:), nfaces, nsolids
        logical,  intent(in) :: flip(:), active(:), is_open(:)
        integer, allocatable :: list(:)
        integer :: f, s, n, k

        if (allocated(rc_faces)) deallocate(rc_faces)
        if (allocated(rc_active)) deallocate(rc_active)
        if (allocated(rc_bvh)) then
            do s = 1, size(rc_bvh)
                call bvh_free(rc_bvh(s))
            end do
            deallocate(rc_bvh)
        end if
        if (allocated(rc_open)) deallocate(rc_open, rc_lo, rc_hi)
        allocate(rc_faces(3,nfaces), rc_active(nfaces), rc_bvh(nsolids), rc_open(nsolids), rc_lo(3,nsolids), rc_hi(3,nsolids))
        rc_faces  = faces(:,1:nfaces)
        rc_active = active(1:nfaces)
        rc_open   = is_open(1:nsolids)
        rc_lo = huge(1.0_dp); rc_hi = -huge(1.0_dp)
        do f = 1, nfaces
            if (.not. rc_active(f)) cycle
            do k = 1, 3
                rc_lo(:,face_solid_id(f)) = min(rc_lo(:,face_solid_id(f)), vertices(:,faces(k,f)))
                rc_hi(:,face_solid_id(f)) = max(rc_hi(:,face_solid_id(f)), vertices(:,faces(k,f)))
            end do
        end do
        do f = 1, nfaces
            if (flip(f)) then
                rc_faces(2,f) = faces(3,f); rc_faces(3,f) = faces(2,f)
            end if
        end do
        if (sign_method == 0) return
        do s = 1, nsolids
            if (solid_is_ground(s)) cycle
            n = count(rc_active .and. face_solid_id(1:nfaces) == s)
            allocate(list(n))
            n = 0
            do f = 1, nfaces
                if (rc_active(f) .and. face_solid_id(f) == s) then
                    n = n + 1; list(n) = f
                end if
            end do
            call bvh_build(rc_bvh(s), vertices, rc_faces, list)
            deallocate(list)
        end do
    end subroutine raycast_set_mesh

    ! Hash column through (coord_b, coord_c): the triangles whose bounding box contains the line. Exact — a
    !   crossing
    ! implies the point lies inside the triangle's projected box, hence in the same hash column.
    subroutine query_hash_column(axis, coord_b, coord_c, candidate_tris, n_candidates)
        implicit none
        integer,  intent(in)    :: axis
        real(dp), intent(in)    :: coord_b, coord_c
        integer,  allocatable, intent(inout) :: candidate_tris(:)
        integer,  intent(out)   :: n_candidates
        integer :: dim_b, dim_c, hb, hc, i_along, cflat, lo, hi, need

        call perp_dims(axis, dim_b, dim_c)
        hb = hash_ix(coord_b, dim_b)
        hc = hash_ix(coord_c, dim_c)
        if (.not. allocated(candidate_tris)) allocate(candidate_tris(256))
        n_candidates = 0
        do i_along = 1, hcells(axis)
            cflat = hcell_flat(axis, i_along, hb, hc)
            lo = csr_offset(cflat); hi = csr_offset(cflat+1) - 1
            if (hi >= lo) then
                need = n_candidates + (hi - lo + 1)
                if (need > size(candidate_tris)) call grow_int(candidate_tris, need)
                candidate_tris(n_candidates+1:need) = csr_triid(lo:hi)
                n_candidates = need
            end if
        end do
    end subroutine query_hash_column

    ! -----------------------------------------------------------------------
    ! Ray/triangle crossing test along `axis`, through (coord_b, coord_c) in
    ! the two perpendicular axes. Winding-independent unless `sgn` is requested: uses only vertex
    ! positions, never the stored/averaged triangle normal. Returns whether
    ! the infinite line crosses the triangle and, if so, the coordinate along
    ! `axis` at the crossing. `sgn` is the sign of the triangle normal's component along the ray (the normal
    !   of
    ! the vertex order v0,v1,v2), which is what the signed winding count needs.
    pure subroutine ray_triangle_crossing(axis, coord_b, coord_c, v0, v1, v2, hit, coord_a, sgn)
        implicit none
        integer,  intent(in)  :: axis
        real(dp), intent(in)  :: coord_b, coord_c
        real(dp), intent(in)  :: v0(3), v1(3), v2(3)
        logical,  intent(out) :: hit
        real(dp), intent(out) :: coord_a
        integer,  intent(out), optional :: sgn
        integer  :: dim_b, dim_c
        real(dp) :: b0, c0, e1b, e1c, e2b, e2c, pb, pc
        real(dp) :: denom, s, t, u
        real(dp) :: edge1(3), edge2(3), nrm(3)

        hit = .false.
        coord_a = 0.0_dp
        if (present(sgn)) sgn = 0
        call perp_dims(axis, dim_b, dim_c)

        b0 = v0(dim_b); c0 = v0(dim_c)
        e1b = v1(dim_b) - b0; e1c = v1(dim_c) - c0
        e2b = v2(dim_b) - b0; e2c = v2(dim_c) - c0
        pb  = coord_b - b0;   pc  = coord_c - c0

        ! 2-D barycentric test of (coord_b,coord_c) against the projected triangle
        denom = e1b*e2c - e2b*e1c
        ! projects to a degenerate line/point (relative test: independent of the length unit of the geometry)
        if (abs(denom) <= 1.0e-14_dp * (e1b*e1b + e1c*e1c + e2b*e2b + e2c*e2c) .or. denom == 0.0_dp) return

        s = (pb*e2c - e2b*pc) / denom
        t = (e1b*pc - pb*e1c) / denom
        u = 1.0_dp - s - t
        if (s < 0.0_dp .or. t < 0.0_dp .or. u < 0.0_dp) return

        edge1 = v1 - v0
        edge2 = v2 - v0
        nrm(1) = edge1(2)*edge2(3) - edge1(3)*edge2(2)
        nrm(2) = edge1(3)*edge2(1) - edge1(1)*edge2(3)
        nrm(3) = edge1(1)*edge2(2) - edge1(2)*edge2(1)
        if (nrm(axis) == 0.0_dp) return   ! triangle parallel to ray axis

        coord_a = v0(axis) - (nrm(dim_b)*(coord_b - v0(dim_b)) + nrm(dim_c)*(coord_c - v0(dim_c))) / nrm(axis)
        hit = .true.
        if (present(sgn)) sgn = merge(1, -1, nrm(axis) > 0.0_dp)
    end subroutine ray_triangle_crossing

    ! -----------------------------------------------------------------------
    ! sign_method = 0: classify every point (xin(xlo:xhi), yin(sy:ey), zin(sz:ez)) as inside or
    !   outside
    ! solid `sid` via a 3-axis ray-crossing-PARITY majority vote. Parity ignores orientation, but it is
    !   flipped for the
    ! rest of any ray through a hole, and cancelled by doubled or overlapping shells. (sx:ex) only fixes the
    !   jitter
    ! scale so every rank perturbs the rays identically; xlo:xhi is this rank's slab.
    subroutine classify_solid(vertices, faces, face_solid_id, nfaces, sid, &
                               xin, yin, zin, sx, ex, sy, ey, sz, ez, xlo, xhi, inside)
        implicit none
        real(dp), intent(in) :: vertices(:,:)
        integer,  intent(in) :: faces(:,:)
        integer,  intent(in) :: face_solid_id(:)
        integer,  intent(in) :: nfaces, sid
        real(dp), intent(in) :: xin(:), yin(:), zin(:)
        integer,  intent(in) :: sx, ex, sy, ey, sz, ez, xlo, xhi
        integer,  intent(out) :: inside(xlo:xhi, sy:ey, sz:ez)   ! vote count, 0-3

        integer, allocatable :: cands(:)
        real(dp), allocatable :: crossings(:)
        integer :: n_cand, n_cross, ii, jj, kk, fc, fid, cap
        real(dp) :: v0(3), v1(3), v2(3)
        logical  :: hit
        real(dp) :: coord_a
        real(dp) :: scale, eps_x, eps_y, eps_z
        integer, allocatable :: last_seen(:)
        integer :: query_gen

        inside = 0
        if (xlo > xhi) return
        allocate(cands(256))

        ! build_hash bins a triangle into every hash cell its AABB overlaps. A crossing-parity count is not
        ! idempotent under that duplication, so each triangle is tagged with the query that last processed it
        !   and
        ! every candidate list is deduplicated in O(1) per candidate, without an O(nfaces) reset between
        !   columns.
        allocate(last_seen(nfaces))
        last_seen = 0
        query_gen = 0

        ! Deterministic per-axis jitter: a query point that lands exactly on the shared edge between two
        !   triangles
        ! (the diagonal of a quad face split into 2 triangles, common whenever a grid point coincides with a
        !   face's
        ! own centre) would register as a crossing of BOTH triangles, doubling the count and flipping parity.
        scale = max(abs(xin(ex)-xin(sx)), abs(yin(ey)-yin(sy)), abs(zin(ez)-zin(sz)), 1.0_dp)
        eps_x = scale * 7.5307258e-8_dp
        eps_y = scale * 3.1415927e-8_dp
        eps_z = scale * 4.1421356e-8_dp

        ! ---- axis = 3 (z): one line per (ii,jj), varying kk ----
        do ii = xlo, xhi
            do jj = sy, ey
                call query_hash_line(3, xin(ii)+eps_x, yin(jj)+eps_y, cands, n_cand)
                query_gen = query_gen + 1
                cap = max(n_cand, 1)
                allocate(crossings(cap))
                n_cross = 0
                do fc = 1, n_cand
                    fid = cands(fc)
                    if (last_seen(fid) == query_gen) cycle
                    last_seen(fid) = query_gen
                    if (face_solid_id(fid) /= sid) cycle
                    v0 = vertices(:, faces(1,fid)); v1 = vertices(:, faces(2,fid)); v2 = vertices(:, faces(3,fid))
                    call ray_triangle_crossing(3, xin(ii)+eps_x, yin(jj)+eps_y, v0, v1, v2, hit, coord_a)
                    if (hit) then
                        n_cross = n_cross + 1
                        crossings(n_cross) = coord_a
                    end if
                end do
                call sort_real(crossings, n_cross)
                call vote_along(crossings, n_cross, zin, sz, ez, inside(ii, jj, :))
                deallocate(crossings)
            end do
        end do

        ! ---- axis = 1 (x): one line per (jj,kk), varying ii ----
        do jj = sy, ey
            do kk = sz, ez
                call query_hash_line(1, yin(jj)+eps_y, zin(kk)+eps_z, cands, n_cand)
                query_gen = query_gen + 1
                cap = max(n_cand, 1)
                allocate(crossings(cap))
                n_cross = 0
                do fc = 1, n_cand
                    fid = cands(fc)
                    if (last_seen(fid) == query_gen) cycle
                    last_seen(fid) = query_gen
                    if (face_solid_id(fid) /= sid) cycle
                    v0 = vertices(:, faces(1,fid)); v1 = vertices(:, faces(2,fid)); v2 = vertices(:, faces(3,fid))
                    call ray_triangle_crossing(1, yin(jj)+eps_y, zin(kk)+eps_z, v0, v1, v2, hit, coord_a)
                    if (hit) then
                        n_cross = n_cross + 1
                        crossings(n_cross) = coord_a
                    end if
                end do
                call sort_real(crossings, n_cross)
                call vote_along(crossings, n_cross, xin, xlo, xhi, inside(:, jj, kk))
                deallocate(crossings)
            end do
        end do

        ! ---- axis = 2 (y): one line per (ii,kk), varying jj ----
        do ii = xlo, xhi
            do kk = sz, ez
                call query_hash_line(2, xin(ii)+eps_x, zin(kk)+eps_z, cands, n_cand)
                query_gen = query_gen + 1
                cap = max(n_cand, 1)
                allocate(crossings(cap))
                n_cross = 0
                do fc = 1, n_cand
                    fid = cands(fc)
                    if (last_seen(fid) == query_gen) cycle
                    last_seen(fid) = query_gen
                    if (face_solid_id(fid) /= sid) cycle
                    v0 = vertices(:, faces(1,fid)); v1 = vertices(:, faces(2,fid)); v2 = vertices(:, faces(3,fid))
                    call ray_triangle_crossing(2, xin(ii)+eps_x, zin(kk)+eps_z, v0, v1, v2, hit, coord_a)
                    if (hit) then
                        n_cross = n_cross + 1
                        crossings(n_cross) = coord_a
                    end if
                end do
                call sort_real(crossings, n_cross)
                call vote_along(crossings, n_cross, yin, sy, ey, inside(ii, :, kk))
                deallocate(crossings)
            end do
        end do

        deallocate(cands, last_seen)
    end subroutine classify_solid

    ! -----------------------------------------------------------------------
    ! Ground-surface classification for an open sheet (terrain/bathymetry) that encloses no volume:
    ! a point is solid when any of the sheet's triangles lie above it along the vertical (axis 3),
    ! i.e. it is below the topmost surface. Needs no closed shell, so open edges and missing caps are
    ! harmless, and unlike parity it is immune to overlapping/duplicated sheets (e.g. touching roofs).
    ! Covers xlo:xhi (this rank's slab) over sy:ey and szc:ez, where szc reaches down to the domain floor.
    subroutine classify_ground(vertices, faces, face_solid_id, nfaces, sid, &
                               xin, yin, zin, sx, ex, sy, ey, sz, ez, xlo, xhi, szc, inside)
        implicit none
        real(dp), intent(in) :: vertices(:,:)
        integer,  intent(in) :: faces(:,:)
        integer,  intent(in) :: face_solid_id(:)
        integer,  intent(in) :: nfaces, sid
        real(dp), intent(in) :: xin(:), yin(:), zin(:)
        integer,  intent(in) :: sx, ex, sy, ey, sz, ez, xlo, xhi, szc
        logical,  intent(out) :: inside(xlo:xhi, sy:ey, szc:ez)

        integer, allocatable :: cands(:), last_seen(:)
        real(dp), allocatable :: crossings(:)
        integer :: n_cand, n_cross, ii, jj, kk, fc, fid, query_gen, nbelow
        logical :: hit
        real(dp) :: coord_a, scale, eps_x, eps_y

        inside = .false.
        if (xlo > xhi) return
        allocate(cands(256), last_seen(nfaces))
        last_seen = 0
        query_gen = 0
        ! Same transverse jitter as classify_solid (avoids exact shared-edge double counts)
        scale = max(abs(xin(ex)-xin(sx)), abs(yin(ey)-yin(sy)), abs(zin(ez)-zin(sz)), 1.0_dp)
        eps_x = scale * 7.5307258e-8_dp
        eps_y = scale * 3.1415927e-8_dp

        do ii = xlo, xhi
            do jj = sy, ey
                call query_hash_column(3, xin(ii)+eps_x, yin(jj)+eps_y, cands, n_cand)
                query_gen = query_gen + 1
                allocate(crossings(max(n_cand, 1)))
                n_cross = 0
                do fc = 1, n_cand
                    fid = cands(fc)
                    if (last_seen(fid) == query_gen) cycle
                    last_seen(fid) = query_gen
                    if (face_solid_id(fid) /= sid) cycle
                    if (allocated(rc_active)) then
                        if (.not. rc_active(fid)) cycle
                    end if
                    call ray_triangle_crossing(3, xin(ii)+eps_x, yin(jj)+eps_y, vertices(:,faces(1,fid)), &
                                               vertices(:,faces(2,fid)), vertices(:,faces(3,fid)), hit, coord_a)
                    if (hit) then
                        n_cross = n_cross + 1
                        crossings(n_cross) = coord_a
                    end if
                end do
                call sort_real(crossings, n_cross)
                nbelow = 0
                do kk = szc, ez
                    do while (nbelow < n_cross)
                        if (crossings(nbelow+1) >= zin(kk)) exit
                        nbelow = nbelow + 1
                    end do
                    if (n_cross > nbelow) inside(ii, jj, kk) = .true.
                end do
                deallocate(crossings)
            end do
        end do

        deallocate(cands, last_seen)
    end subroutine classify_ground

    ! Given sorted crossing coordinates, sweep the (already-sorted) grid
    ! coordinates line(lo:hi) and accumulate a +1 vote wherever the parity of
    ! crossings strictly below that coordinate is odd (inside).
    subroutine vote_along(crossings, n_cross, line, lo, hi, votes)
        implicit none
        real(dp), intent(in)    :: crossings(:)
        integer,  intent(in)    :: n_cross, lo, hi
        real(dp), intent(in)    :: line(:)
        integer,  intent(inout) :: votes(lo:hi)
        integer :: idx, ptr

        ptr = 0
        do idx = lo, hi
            do while (ptr < n_cross)
                if (crossings(ptr+1) >= line(idx)) exit
                ptr = ptr + 1
            end do
            if (mod(ptr, 2) == 1) votes(idx) = votes(idx) + 1
        end do
    end subroutine vote_along

    ! -----------------------------------------------------------------------
    ! Tier 1 of the robust classifier for one grid line: signed crossings of solid `sid` along `axis` through
    !   (cb,cc).
    ! For each point line(i) two half-line rays are evaluated (towards -axis and +axis); each contributes to
    !   lcode(i)
    !   1  if its winding number is non-zero (inside),   8  if it is zero (outside),   0  if the ray abstains.
    ! A ray with no crossing at all is a clean "outside" for a closed shell; for an open shell it abstains
    !   while the
    ! point is inside the solid's bounding box (it may have left through a hole).
    subroutine line_codes(axis, cb, cc, line, lo, hi, vertices, sid, face_solid_id, nfaces, &
                          cands, cross_c, cross_s, last_seen, query_gen, lcode)
        implicit none
        integer,  intent(in)    :: axis, lo, hi, sid, nfaces
        real(dp), intent(in)    :: cb, cc, line(:), vertices(:,:)
        integer,  intent(in)    :: face_solid_id(:)
        integer,  allocatable, intent(inout) :: cands(:), cross_s(:)
        real(dp), allocatable, intent(inout) :: cross_c(:)
        integer,  intent(inout) :: last_seen(:), query_gen
        integer(1), intent(out) :: lcode(lo:hi)
        integer :: n_cand, n_cross, fc, fid, sg, idx, ptr, sbelow, total, dim_b, dim_c, code
        real(dp) :: ca
        logical :: hit, transverse_in_box, may_abstain

        call query_hash_column(axis, cb, cc, cands, n_cand)
        query_gen = query_gen + 1
        if (n_cand > size(cross_c)) then
            deallocate(cross_c, cross_s); allocate(cross_c(2*n_cand), cross_s(2*n_cand))
        end if
        n_cross = 0
        do fc = 1, n_cand
            fid = cands(fc)
            if (last_seen(fid) == query_gen) cycle
            last_seen(fid) = query_gen
            if (face_solid_id(fid) /= sid) cycle
            if (.not. rc_active(fid)) cycle
            call ray_triangle_crossing(axis, cb, cc, vertices(:,rc_faces(1,fid)), vertices(:,rc_faces(2,fid)), &
                                       vertices(:,rc_faces(3,fid)), hit, ca, sg)
            if (hit) then
                n_cross = n_cross + 1
                cross_c(n_cross) = ca; cross_s(n_cross) = sg
            end if
        end do
        if (n_cross > 1) call sort_pairs(cross_c, cross_s, n_cross)
        total = sum(cross_s(1:n_cross))

        call perp_dims(axis, dim_b, dim_c)
        transverse_in_box = rc_open(sid) .and. cb >= rc_lo(dim_b,sid) .and. cb <= rc_hi(dim_b,sid) &
                                         .and. cc >= rc_lo(dim_c,sid) .and. cc <= rc_hi(dim_c,sid)
        ptr = 0; sbelow = 0
        do idx = lo, hi
            do while (ptr < n_cross)
                if (cross_c(ptr+1) >= line(idx)) exit
                ptr = ptr + 1
                sbelow = sbelow + cross_s(ptr)
            end do
            may_abstain = transverse_in_box
            if (may_abstain) may_abstain = (line(idx) >= rc_lo(axis,sid) .and. line(idx) <= rc_hi(axis,sid))
            code = 0
            ! ray towards -axis: winding = minus the signed crossings below
            if (ptr > 0) then
                code = code + merge(1, 8, sbelow /= 0)
            else if (.not. may_abstain) then
                code = code + 8
            end if
            ! ray towards +axis: winding = signed crossings above
            if (ptr < n_cross) then
                code = code + merge(1, 8, total - sbelow /= 0)
            else if (.not. may_abstain) then
                code = code + 8
            end if
            lcode(idx) = int(code, 1)
        end do
    end subroutine line_codes

    ! -----------------------------------------------------------------------
    ! Robust classification of solid `sid` (sign_method 1 or 2): tier-1 signed rays + tier-2 winding number.
    ! Returns inside(xlo:xhi, sy:ey, szc:ez) for this rank's slab; all ranks must call (collective).
    subroutine classify_solid_robust(vertices, face_solid_id, nfaces, sid, xin, yin, zin, &
                                     sx, ex, sy, ey, sz, ez, xlo, xhi, szc, myid, nproc, comm, inside)
        implicit none
        real(dp), intent(in) :: vertices(:,:)
        integer,  intent(in) :: face_solid_id(:), nfaces, sid
        real(dp), intent(in) :: xin(:), yin(:), zin(:)
        integer,  intent(in) :: sx, ex, sy, ey, sz, ez, xlo, xhi, szc, myid, nproc, comm
        logical,  intent(out) :: inside(xlo:xhi, sy:ey, szc:ez)

        integer(1), allocatable :: code(:,:,:)
        integer(1), allocatable :: lcode(:)
        integer, allocatable :: cands(:), cross_s(:), last_seen(:)
        real(dp), allocatable :: cross_c(:)
        integer :: ii, jj, kk, query_gen, lc, n_in, n_out, ierr
        integer :: counts(5), gcounts(5)
        real(dp) :: scale, eps_x, eps_y, eps_z, w, p(3)
        logical :: ray_in, ray_out, decided

        allocate(code(sx:ex, sy:ey, szc:ez)); code = 0_1
        allocate(cands(256), cross_c(256), cross_s(256), last_seen(nfaces))
        last_seen = 0; query_gen = 0
        scale = max(abs(xin(ex)-xin(sx)), abs(yin(ey)-yin(sy)), abs(zin(ez)-zin(sz)), 1.0_dp)
        eps_x = scale * 7.5307258e-8_dp
        eps_y = scale * 3.1415927e-8_dp
        eps_z = scale * 4.1421356e-8_dp

        if (sign_method >= 1) then
            ! ---- tier 1: lines round-robin over ranks ----
            lc = 0
            allocate(lcode(szc:ez))
            do ii = sx, ex
                do jj = sy, ey
                    lc = lc + 1
                    if (mod(lc, nproc) /= myid) cycle
                    call line_codes(3, xin(ii)+eps_x, yin(jj)+eps_y, zin, szc, ez, vertices, sid, face_solid_id, nfaces, &
                                    cands, cross_c, cross_s, last_seen, query_gen, lcode)
                    code(ii,jj,:) = code(ii,jj,:) + lcode
                end do
            end do
            deallocate(lcode)
            allocate(lcode(sx:ex))
            do jj = sy, ey
                do kk = szc, ez
                    lc = lc + 1
                    if (mod(lc, nproc) /= myid) cycle
                    call line_codes(1, yin(jj)+eps_y, zin(kk)+eps_z, xin, sx, ex, vertices, sid, face_solid_id, nfaces, &
                                    cands, cross_c, cross_s, last_seen, query_gen, lcode)
                    code(:,jj,kk) = code(:,jj,kk) + lcode
                end do
            end do
            deallocate(lcode)
            allocate(lcode(sy:ey))
            do ii = sx, ex
                do kk = szc, ez
                    lc = lc + 1
                    if (mod(lc, nproc) /= myid) cycle
                    call line_codes(2, xin(ii)+eps_x, zin(kk)+eps_z, yin, sy, ey, vertices, sid, face_solid_id, nfaces, &
                                    cands, cross_c, cross_s, last_seen, query_gen, lcode)
                    code(ii,:,kk) = code(ii,:,kk) + lcode
                end do
            end do
            deallocate(lcode)
            if (nproc > 1) call allreduce_i1(code, size(code), comm)
        end if

        ! ---- decide: rays where unanimous, winding number elsewhere ----
        counts = 0      ! 1: ray inside, 2: ray outside, 3: winding number, 4: winding number near 1/2, 5: tier-1 vs winding mismatch
        inside = .false.
        do kk = szc, ez
            do jj = sy, ey
                do ii = max(xlo, sx), min(xhi, ex)
                    n_in  = mod(int(code(ii,jj,kk)), 8)
                    n_out = int(code(ii,jj,kk)) / 8
                    ray_in  = (n_in >= 3 .and. n_out == 0)
                    ray_out = (n_out >= 3 .and. n_in == 0)
                    decided = ray_in .or. ray_out
                    if (decided .and. sign_method == 1) then
                        inside(ii,jj,kk) = ray_in
                        if (ray_in) then; counts(1) = counts(1) + 1; else; counts(2) = counts(2) + 1; end if
                    else
                        p = [xin(ii), yin(jj), zin(kk)]
                        w = abs(winding_number(rc_bvh(sid), vertices, rc_faces, p))
                        inside(ii,jj,kk) = (w > 0.5_dp)
                        counts(3) = counts(3) + 1
                        if (abs(w - 0.5_dp) < 0.1_dp) counts(4) = counts(4) + 1
                        if (decided .and. (ray_in .neqv. inside(ii,jj,kk))) counts(5) = counts(5) + 1
                    end if
                end do
            end do
        end do
        call MPI_REDUCE(counts, gcounts, 5, MPI_INTEGER, MPI_SUM, 0, comm, ierr)
        if (myid == 0) then
            if (sign_method == 1) then
                print '(A,I0,A,I0,A,I0,A,I0,A,I0,A)', "    solid ", sid, ": ", gcounts(1)+gcounts(2), &
                    " points decided by signed rays, ", gcounts(3), " by winding number (", gcounts(4), " near the 1/2 threshold)"
            else
                print '(A,I0,A,I0,A,I0,A)', "    solid ", sid, ": winding number at all ", gcounts(3), &
                    " points; ", gcounts(5), " disagree with an unanimous signed-ray decision"
            end if
        end if
        deallocate(code, cands, cross_c, cross_s, last_seen)
    end subroutine classify_solid_robust

    ! MPI_SUM of a (large) int8 array, chunked to stay within the default-integer count limit
    subroutine allreduce_i1(buf, n, comm)
        integer(1), intent(inout) :: buf(*)
        integer,    intent(in)    :: n, comm
        integer, parameter :: chunk = 500000000
        integer :: off, cnt, ierr
        off = 0
        do while (off < n)
            cnt = min(chunk, n - off)
            call MPI_ALLREDUCE(MPI_IN_PLACE, buf(off+1), cnt, MPI_INTEGER1, MPI_SUM, comm, ierr)
            off = off + cnt
        end do
    end subroutine allreduce_i1

    ! -----------------------------------------------------------------------
    ! Top-level entry point: classify every point of this rank's x-slab xlo:xhi (intersected with the tagged
    !   box
    ! sx:ex, sy:ey, sz:ez; for `ground` solids the box reaches down to the domain floor) as solid or fluid,
    !   one
    ! solid at a time (first solid to claim a point wins), and write the result into `grid`/`objid_grid`,
    !   which hold
    ! this slab (xlo:xhi, 1:ny, 1:nz). Cells that already carry a real near-surface distance (from
    ! compute_narrowband_sdf / compute_scalar_distance_face) keep their magnitude, with the sign forced to
    !   match the
    ! classification; cells still at the sentinel value keep the sentinel magnitude (fast_sweep_3d fills those
    !   in
    ! afterwards). Collective over comm: every rank must call it, including ranks whose slab is empty (xlo >
    !   xhi).
    subroutine raycast_classify_domain(vertices, faces, face_solid_id, nfaces, nsolids, &
                                        xin, yin, zin, sx, ex, sy, ey, sz, ez, xlo, xhi, &
                                        scalarvalue_in, grid, objid_grid, myid, nproc, comm)
        implicit none
        real(dp), intent(in)    :: vertices(:,:)
        integer,  intent(in)    :: faces(:,:)
        integer,  intent(in)    :: face_solid_id(:)
        integer,  intent(in)    :: nfaces, nsolids
        real(dp), intent(in)    :: xin(:), yin(:), zin(:)
        integer,  intent(in)    :: sx, ex, sy, ey, sz, ez, xlo, xhi, myid, nproc, comm
        real(dp), intent(in)    :: scalarvalue_in
        real(dp), intent(inout) :: grid(xlo:xhi, 1:size(yin), 1:size(zin))
        real(dp), intent(inout) :: objid_grid(xlo:xhi, 1:size(yin), 1:size(zin))

        integer, allocatable :: votes(:,:,:)
        logical, allocatable :: any_inside(:,:,:), inside_s(:,:,:)
        real(dp), allocatable :: winning_id(:,:,:)
        integer :: sid, ii, jj, kk, szc, bx_lo, bx_hi

        szc = sz
        if (any(solid_is_ground(1:nsolids))) szc = 1      ! an open ground sheet has solid all the way down to the floor
        bx_lo = max(xlo, sx); bx_hi = min(xhi, ex)

        allocate(inside_s(bx_lo:bx_hi, sy:ey, szc:ez))
        allocate(any_inside(bx_lo:bx_hi, sy:ey, szc:ez))
        allocate(winning_id(bx_lo:bx_hi, sy:ey, szc:ez))
        any_inside = .false.
        winning_id = 0.0_dp

        if (myid == 0 .and. sign_method /= 0) print *, "*** Solid/fluid classification (sign_method =", sign_method, ") ***"
        do sid = 1, nsolids
            if (solid_is_ground(sid)) then
                call classify_ground(vertices, faces, face_solid_id, nfaces, sid, &
                                     xin, yin, zin, sx, ex, sy, ey, sz, ez, bx_lo, bx_hi, szc, inside_s)
            else if (sign_method == 0) then
                allocate(votes(bx_lo:bx_hi, sy:ey, sz:ez))
                call classify_solid(vertices, faces, face_solid_id, nfaces, sid, &
                                    xin, yin, zin, sx, ex, sy, ey, sz, ez, bx_lo, bx_hi, votes)
                inside_s = .false.
                inside_s(:,:,sz:ez) = (votes >= 2)
                deallocate(votes)
            else
                call classify_solid_robust(vertices, face_solid_id, nfaces, sid, xin, yin, zin, &
                                           sx, ex, sy, ey, sz, ez, bx_lo, bx_hi, szc, myid, nproc, comm, inside_s)
            end if
            do kk = szc, ez
                do jj = sy, ey
                    do ii = bx_lo, bx_hi
                        if (.not. any_inside(ii,jj,kk) .and. inside_s(ii,jj,kk)) then
                            any_inside(ii,jj,kk) = .true.
                            winning_id(ii,jj,kk) = real(sid, dp)
                        end if
                    end do
                end do
            end do
        end do

        do kk = szc, ez
            do jj = sy, ey
                do ii = bx_lo, bx_hi
                    if (any_inside(ii,jj,kk)) then
                        if (abs(grid(ii,jj,kk)) < scalarvalue_in * 0.5_dp) then
                            grid(ii,jj,kk) = -abs(grid(ii,jj,kk))
                        else
                            grid(ii,jj,kk) = -scalarvalue_in
                        end if
                        objid_grid(ii,jj,kk) = winning_id(ii,jj,kk)
                    else
                        if (abs(grid(ii,jj,kk)) < scalarvalue_in * 0.5_dp) then
                            grid(ii,jj,kk) = abs(grid(ii,jj,kk))
                        else
                            grid(ii,jj,kk) = scalarvalue_in
                        end if
                        objid_grid(ii,jj,kk) = 0.0_dp
                    end if
                end do
            end do
        end do

        deallocate(inside_s, any_inside, winning_id)
    end subroutine raycast_classify_domain

    ! Small-n insertion sort (crossing counts per line are typically tiny)
    pure subroutine sort_real(arr, n)
        real(dp), intent(inout) :: arr(:)
        integer,  intent(in)    :: n
        integer  :: i, j
        real(dp) :: key
        do i = 2, n
            key = arr(i); j = i - 1
            do while (j >= 1)
                if (arr(j) <= key) exit
                arr(j+1) = arr(j); j = j - 1
            end do
            arr(j+1) = key
        end do
    end subroutine sort_real

    ! Sort crossing coordinates ascending, carrying the crossing sign along (median-of-3 quicksort, insertion
    !   for small
    ! ranges; crossing counts per line are tiny for typical solids but reach thousands for finely tessellated
    !   terrain).
    recursive subroutine sort_pairs(c, s, n)
        real(dp), intent(inout) :: c(:)
        integer,  intent(inout) :: s(:)
        integer,  intent(in)    :: n
        call qsort_pairs(c, s, 1, n)
    end subroutine sort_pairs

    recursive subroutine qsort_pairs(c, s, lo, hi)
        real(dp), intent(inout) :: c(:)
        integer,  intent(inout) :: s(:)
        integer,  intent(in)    :: lo, hi
        integer :: i, j, m, si
        real(dp) :: pivot, ci
        if (hi - lo < 24) then
            do i = lo + 1, hi
                ci = c(i); si = s(i); j = i - 1
                do while (j >= lo)
                    if (c(j) <= ci) exit
                    c(j+1) = c(j); s(j+1) = s(j); j = j - 1
                end do
                c(j+1) = ci; s(j+1) = si
            end do
            return
        end if
        m = (lo + hi) / 2
        pivot = max(min(c(lo), c(m)), min(max(c(lo), c(m)), c(hi)))
        i = lo; j = hi
        do
            do while (c(i) < pivot); i = i + 1; end do
            do while (c(j) > pivot); j = j - 1; end do
            if (i <= j) then
                ci = c(i); c(i) = c(j); c(j) = ci
                si = s(i); s(i) = s(j); s(j) = si
                i = i + 1; j = j - 1
            end if
            if (i > j) exit
        end do
        if (lo < j) call qsort_pairs(c, s, lo, j)
        if (i < hi) call qsort_pairs(c, s, i, hi)
    end subroutine qsort_pairs

end module raycast_mod
