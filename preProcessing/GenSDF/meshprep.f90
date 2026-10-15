module meshprep_mod
    !
    ! meshprep_mod — mesh health analysis and repair ahead of the solid/fluid classification.
    !
    ! Real-world surface meshes are rarely clean: STL files carry three private copies of every vertex, exporters
    ! leave pinholes and cracks, patches are modelled twice, and triangle orientation is often mixed. The ray and
    ! winding-number classifiers in raycast_mod need to know what they are dealing with, so this module
    !
    !   1. welds coincident vertices (per solid, within a small tolerance, 27-cell neighbour search so that points
    !      straddling a quantisation boundary still merge),
    !   2. drops degenerate triangles and cancels duplicated ones (same welded vertices; copies of opposite
    !      orientation annihilate, copies of equal orientation collapse to one),
    !   3. builds the edge graph and reports boundary edges (holes), non-manifold edges and inconsistently oriented
    !      edges, with the location and size of every hole loop (and an OBJ of the open edges to look at),
    !   4. orients each connected patch consistently by flood-filling across manifold edges. A patch whose file
    !      orientation contradicts itself is made outward-facing (flux of the position vector about the patch
    !      centroid); a patch that already was consistent is left as the file gave it, so a hollow body whose
    !      cavity surface faces inward keeps that meaning.
    !
    ! Output per triangle: active (kept) and flip (swap vertices 2 and 3). The geometry arrays are not modified.
    !
    use utils_io, only : dp, vertical_axis
    implicit none
    private
    public :: prepare_mesh

    integer, parameter :: i8 = selected_int_kind(18)

contains

    ! ------------------------------------------------------------------------------------------------------------------
    subroutine prepare_mesh(procid, vertices, faces, face_solid_id, nfaces, nsolids, is_ground, weld_tol, crack_dist, &
                            do_cap, flip, active, is_open, cap_faces, cap_solid, report_dir)
        implicit none
        integer,  intent(in)  :: procid
        real(dp), intent(in)  :: vertices(:,:)
        integer,  intent(in)  :: faces(:,:), face_solid_id(:), nfaces, nsolids
        logical,  intent(in)  :: is_ground(:), do_cap
        real(dp), intent(in)  :: weld_tol, crack_dist     ! boundary loops closer than crack_dist are seams, not holes
        logical,  allocatable, intent(out) :: flip(:), active(:)
        logical,  intent(out) :: is_open(:)               ! solid is not a closed manifold shell (after capping)
        integer,  allocatable, intent(out) :: cap_faces(:,:), cap_solid(:)   ! lid triangles (vertex ids) and their solid
        character(len=*), intent(in) :: report_dir        ! directory for open_edges.obj ('' = do not write)

        integer,  allocatable :: rep(:), wf(:,:), order(:), nbr(:,:), comp(:), queue(:)
        integer,  allocatable :: rec_tri(:), rec_slot(:), edge_a(:), edge_b(:), edge_s(:), edge_f(:), edge_state(:)
        integer,  allocatable :: n_deg(:), n_dup(:), n_bnd(:), n_nonman(:), n_incons(:), n_flip(:), n_patch(:), n_fixed(:)
        integer,  allocatable :: n_left(:), n_lids(:), n_lidtri(:)
        integer(i8), allocatable :: key(:)
        logical,  allocatable :: nrel(:,:), fl(:), comp_bad(:), cavity(:)
        real(dp), allocatable :: cs(:,:), wsum(:), outw(:)
        real(dp) :: ext, minarea2, av(3), ce(3), w
        integer(i8) :: nvp
        integer :: nv, f, k, a, b, c, ne, nact, nbe, n_conflict
        integer :: i, j, grp, ncomp, head, tail, u, v, slot, t1, t2, s1, s2
        logical :: d1, d2

        nv = size(vertices, 2)
        allocate(flip(nfaces), active(nfaces))
        flip = .false.; active = .true.; is_open = .false.
        if (nfaces == 0) return
        allocate(n_deg(nsolids), n_dup(nsolids), n_bnd(nsolids), n_nonman(nsolids), n_incons(nsolids), &
                 n_flip(nsolids), n_patch(nsolids), n_fixed(nsolids), n_left(nsolids), n_lids(nsolids), n_lidtri(nsolids))
        n_deg = 0; n_dup = 0; n_bnd = 0; n_nonman = 0; n_incons = 0; n_flip = 0; n_patch = 0; n_fixed = 0
        n_left = 0; n_lids = 0; n_lidtri = 0
        allocate(cap_faces(3,0), cap_solid(0))

        ! ---- 1. weld vertices (vertex ids are solid-specific because each solid's vertices are appended separately) ----
        allocate(rep(nv), wf(3,nfaces))
        call weld_vertices(vertices, faces, face_solid_id, nfaces, weld_tol, rep)
        do f = 1, nfaces
            wf(:,f) = rep(faces(:,f))
        end do
        nvp = int(nv, i8) + 1_i8

        ! ---- 2. degenerate and duplicate triangles ----
        ext = maxval(vertices) - minval(vertices)
        minarea2 = (1.0e-12_dp * max(ext, tiny(1.0_dp)))**4          ! |e1 x e2|^2 below this: zero-area sliver
        do f = 1, nfaces
            s1 = face_solid_id(f)
            if (wf(1,f) == wf(2,f) .or. wf(2,f) == wf(3,f) .or. wf(1,f) == wf(3,f)) then
                active(f) = .false.; n_deg(s1) = n_deg(s1) + 1; cycle
            end if
            call tri_area_vec_centre(vertices, faces(:,f), av, ce)
            if (4.0_dp*dot_product(av,av) < minarea2) then
                active(f) = .false.; n_deg(s1) = n_deg(s1) + 1
            end if
        end do
        call cancel_duplicates(wf, face_solid_id, nfaces, nvp, active, n_dup)

        ! ---- 3. edge graph: one record per active-triangle edge, sorted by (low vertex, high vertex) ----
        nact = count(active)
        ne = 3 * nact
        allocate(key(ne), order(ne), rec_tri(ne), rec_slot(ne))
        i = 0
        do f = 1, nfaces
            if (.not. active(f)) cycle
            do k = 1, 3
                a = wf(k,f); b = wf(mod(k,3)+1, f)
                i = i + 1
                key(i) = int(min(a,b), i8) * nvp + int(max(a,b), i8)
                rec_tri(i) = f; rec_slot(i) = k; order(i) = i
            end do
        end do
        call sort_keys(key, order, ne)

        allocate(nbr(3,nfaces), nrel(3,nfaces), edge_a(ne), edge_b(ne), edge_s(ne), edge_f(ne), edge_state(ne))
        nbr = 0; nrel = .false.; nbe = 0
        i = 1
        do while (i <= ne)
            j = i
            do while (j < ne)
                if (key(j+1) /= key(i)) exit
                j = j + 1
            end do
            grp = j - i + 1
            t1 = rec_tri(order(i)); slot = rec_slot(order(i))
            s1 = face_solid_id(t1)
            if (grp == 1) then
                n_bnd(s1) = n_bnd(s1) + 1
                nbe = nbe + 1
                edge_a(nbe) = wf(slot,t1); edge_b(nbe) = wf(mod(slot,3)+1,t1); edge_s(nbe) = s1; edge_f(nbe) = t1
            else if (grp == 2) then
                t2 = rec_tri(order(i+1)); k = rec_slot(order(i+1))
                d1 = wf(slot,t1) < wf(mod(slot,3)+1,t1)
                d2 = wf(k,t2)    < wf(mod(k,3)+1,t2)
                if (t1 /= t2) then
                    ! consistently oriented neighbours traverse their shared edge in opposite directions
                    nbr(slot,t1) = t2; nrel(slot,t1) = (d1 .eqv. d2)
                    nbr(k,t2)    = t1; nrel(k,t2)    = (d1 .eqv. d2)
                    if (d1 .eqv. d2) n_incons(s1) = n_incons(s1) + 1
                end if
            else
                n_nonman(s1) = n_nonman(s1) + 1
            end if
            i = j + 1
        end do

        ! ---- 4. orient patches by flood fill across manifold edges ----
        allocate(comp(nfaces), fl(nfaces), queue(nfaces))
        comp = 0; fl = .false.
        ncomp = 0; n_conflict = 0
        do f = 1, nfaces
            if (.not. active(f) .or. comp(f) /= 0) cycle
            ncomp = ncomp + 1
            comp(f) = ncomp; fl(f) = .false.
            n_patch(face_solid_id(f)) = n_patch(face_solid_id(f)) + 1
            head = 1; tail = 1; queue(1) = f
            do while (head <= tail)
                u = queue(head); head = head + 1
                do slot = 1, 3
                    v = nbr(slot,u)
                    if (v == 0) cycle
                    ! v must be flipped (relative to its file order) iff u is flipped XOR the file orientations disagree
                    if (comp(v) == 0) then
                        comp(v) = ncomp
                        fl(v) = (fl(u) .neqv. nrel(slot,u))
                        tail = tail + 1; queue(tail) = v
                    else if (fl(v) .neqv. (fl(u) .neqv. nrel(slot,u))) then
                        n_conflict = n_conflict + 1          ! each such edge is seen from both sides
                    end if
                end do
            end do
        end do
        n_conflict = n_conflict / 2

        ! a patch needs repair if the file orientation contradicts itself somewhere on it
        allocate(comp_bad(ncomp))
        comp_bad = .false.
        do f = 1, nfaces
            if (active(f)) then
                if (any(nrel(:,f) .and. nbr(:,f) /= 0)) comp_bad(comp(f)) = .true.
            end if
        end do

        ! repaired patches are made outward-facing: flux of (triangle centre - patch centroid) through the patch
        allocate(cs(3,ncomp), wsum(ncomp), outw(ncomp))
        cs = 0.0_dp; wsum = 0.0_dp; outw = 0.0_dp
        do f = 1, nfaces
            if (.not. active(f)) cycle
            c = comp(f)
            if (.not. comp_bad(c)) cycle
            call tri_area_vec_centre(vertices, faces(:,f), av, ce)
            w = sqrt(dot_product(av,av))
            cs(:,c) = cs(:,c) + w * ce
            wsum(c) = wsum(c) + w
        end do
        do c = 1, ncomp
            if (wsum(c) > 0.0_dp) cs(:,c) = cs(:,c) / wsum(c)
        end do
        do f = 1, nfaces
            if (.not. active(f)) cycle
            c = comp(f)
            if (.not. comp_bad(c)) cycle
            call tri_area_vec_centre(vertices, faces(:,f), av, ce)
            if (fl(f)) av = -av
            outw(c) = outw(c) + dot_product(av, ce - cs(:,c))
        end do
        ! A repaired patch that sits inside an odd number of other patches of its solid is a cavity wall: its surface
        ! must face into the cavity, i.e. the opposite of the "outward" choice (nesting parity from a ray).
        allocate(cavity(ncomp))
        cavity = .false.
        call nesting_parity(vertices, faces, face_solid_id, nfaces, active, comp, comp_bad, ncomp, cavity)
        do c = 1, ncomp
            if (comp_bad(c) .and. cavity(c)) outw(c) = -outw(c)
        end do
        do f = 1, nfaces
            if (.not. active(f)) cycle
            c = comp(f)
            if (.not. comp_bad(c)) cycle
            flip(f) = fl(f)
            if (outw(c) < 0.0_dp) flip(f) = .not. flip(f)
            if (flip(f)) n_flip(face_solid_id(f)) = n_flip(face_solid_id(f)) + 1
        end do
        do c = 1, ncomp
            if (comp_bad(c)) then
                do f = 1, nfaces
                    if (active(f) .and. comp(f) == c) then
                        n_fixed(face_solid_id(f)) = n_fixed(face_solid_id(f)) + 1; exit
                    end if
                end do
            end if
        end do

        ! ---- 5. lids over planar boundary loops ----
        edge_state = 0
        n_left = n_bnd
        if (do_cap .and. nbe > 0) then
            call cap_boundary_loops(vertices, nbe, edge_a, edge_b, edge_s, edge_f, flip, is_ground, nsolids, nv, crack_dist, &
                                    cap_faces, cap_solid, edge_state, n_lids, n_lidtri)
            do i = 1, nbe
                if (edge_state(i) == 1) n_left(edge_s(i)) = n_left(edge_s(i)) - 1
            end do
        end if
        is_open = (n_left > 0 .or. n_nonman > 0) .and. .not. is_ground(1:nsolids)
        if (procid == 0) &
            call print_report(vertices, nfaces, nsolids, is_ground, face_solid_id, active, n_deg, n_dup, n_bnd, n_nonman, &
                              n_incons, n_flip, n_patch, n_fixed, n_conflict, nbe, edge_a, edge_b, edge_s, edge_state, &
                              n_lids, n_lidtri, n_left, do_cap, report_dir)
    end subroutine prepare_mesh


    ! For every patch flagged in `want`, decide whether it lies inside an odd number of other patches of its solid.
    ! Cast a ray in +x from the patch's extreme vertex (so the patch itself cannot be met) and count crossings of the
    ! other patches' triangles; the transverse position gets a tiny irrational offset to avoid grazing an edge.
    subroutine nesting_parity(vertices, faces, face_solid_id, nfaces, active, comp, want, ncomp, odd)
        real(dp), intent(in)  :: vertices(:,:)
        integer,  intent(in)  :: faces(:,:), face_solid_id(:), nfaces, comp(:), ncomp
        logical,  intent(in)  :: active(:), want(:)
        logical,  intent(out) :: odd(:)
        integer,  allocatable :: apex(:)
        real(dp), allocatable :: apx(:)
        integer :: f, k, c, v, ncross, sc
        real(dp) :: py, pz, p0(3), v1(3), v2(3), v3(3), s, t, u, x
        real(dp) :: e1y, e1z, e2y, e2z, denom, e1(3), e2(3), nr(3)

        odd = .false.
        allocate(apex(ncomp), apx(ncomp))
        apex = 0; apx = -huge(1.0_dp)
        do f = 1, nfaces
            if (.not. active(f)) cycle
            c = comp(f)
            if (.not. want(c)) cycle
            do k = 1, 3
                v = faces(k,f)
                if (vertices(1,v) > apx(c)) then; apx(c) = vertices(1,v); apex(c) = v; end if
            end do
        end do
        do c = 1, ncomp
            if (.not. want(c) .or. apex(c) == 0) cycle
            p0 = vertices(:,apex(c))
            sc = face_solid_id(first_face_of(c, comp, active, nfaces))
            py = p0(2) + 1.234567e-9_dp * (1.0_dp + abs(p0(2)))
            pz = p0(3) + 2.718281e-9_dp * (1.0_dp + abs(p0(3)))
            ncross = 0
            do f = 1, nfaces
                if (.not. active(f)) cycle
                if (comp(f) == c) cycle
                if (face_solid_id(f) /= sc) cycle
                v1 = vertices(:,faces(1,f)); v2 = vertices(:,faces(2,f)); v3 = vertices(:,faces(3,f))
                if (max(v1(1),v2(1),v3(1)) < p0(1)) cycle
                e1y = v2(2)-v1(2); e1z = v2(3)-v1(3); e2y = v3(2)-v1(2); e2z = v3(3)-v1(3)
                denom = e1y*e2z - e2y*e1z
                if (denom == 0.0_dp) cycle
                s = ((py-v1(2))*e2z - e2y*(pz-v1(3))) / denom
                t = (e1y*(pz-v1(3)) - (py-v1(2))*e1z) / denom
                u = 1.0_dp - s - t
                if (s < 0.0_dp .or. t < 0.0_dp .or. u < 0.0_dp) cycle
                e1 = v2 - v1; e2 = v3 - v1
                nr = [e1(2)*e2(3)-e1(3)*e2(2), e1(3)*e2(1)-e1(1)*e2(3), e1(1)*e2(2)-e1(2)*e2(1)]
                if (nr(1) == 0.0_dp) cycle
                x = v1(1) - (nr(2)*(py-v1(2)) + nr(3)*(pz-v1(3))) / nr(1)
                if (x > p0(1)) ncross = ncross + 1
            end do
            odd(c) = (mod(ncross, 2) == 1)
        end do
        deallocate(apex, apx)
    end subroutine nesting_parity

    integer function first_face_of(c, comp, active, nfaces)
        integer, intent(in) :: c, comp(:), nfaces
        logical, intent(in) :: active(:)
        integer :: f
        first_face_of = 1
        do f = 1, nfaces
            if (active(f) .and. comp(f) == c) then; first_face_of = f; return; end if
        end do
    end function first_face_of

    ! ------------------------------------------------------------------------------------------------------------------
    subroutine tri_area_vec_centre(vertices, tri, av, ce)
        real(dp), intent(in)  :: vertices(:,:)
        integer,  intent(in)  :: tri(3)
        real(dp), intent(out) :: av(3), ce(3)
        real(dp) :: e1(3), e2(3)
        e1 = vertices(:,tri(2)) - vertices(:,tri(1))
        e2 = vertices(:,tri(3)) - vertices(:,tri(1))
        av = 0.5_dp * [e1(2)*e2(3)-e1(3)*e2(2), e1(3)*e2(1)-e1(1)*e2(3), e1(1)*e2(2)-e1(2)*e2(1)]
        ce = (vertices(:,tri(1)) + vertices(:,tri(2)) + vertices(:,tri(3))) / 3.0_dp
    end subroutine tri_area_vec_centre

    ! ------------------------------------------------------------------------------------------------------------------
    ! Merge vertices closer than tol (Chebyshev). Chained hash of tol-sized cells; each vertex looks into its 27
    ! neighbour cells so two points a hair apart on either side of a cell boundary still merge.
    ! rep(v) = the earliest vertex it merged into.
    subroutine weld_vertices(vertices, faces, face_solid_id, nfaces, tol, rep)
        real(dp), intent(in)  :: vertices(:,:)
        integer,  intent(in)  :: faces(:,:), face_solid_id(:), nfaces
        real(dp), intent(in)  :: tol
        integer,  intent(out) :: rep(:)
        integer,  allocatable :: head(:), nxt(:), vsolid(:)
        integer(i8) :: cx, cy, cz
        integer :: nv, v, u, f, k, ddx, ddy, ddz, h, msize
        real(dp) :: lo(3), inv
        logical :: found

        nv = size(vertices, 2)
        allocate(vsolid(nv)); vsolid = 0
        do f = 1, nfaces
            do k = 1, 3
                vsolid(faces(k,f)) = face_solid_id(f)
            end do
        end do
        lo = [minval(vertices(1,:)), minval(vertices(2,:)), minval(vertices(3,:))]
        inv = 1.0_dp / max(tol, tiny(1.0_dp))
        msize = 1
        do while (msize < 2*nv)
            msize = msize * 2
        end do
        allocate(head(0:msize-1), nxt(nv))
        head = 0; nxt = 0
        do v = 1, nv
            rep(v) = v
            cx = int(floor((vertices(1,v)-lo(1))*inv), i8)
            cy = int(floor((vertices(2,v)-lo(2))*inv), i8)
            cz = int(floor((vertices(3,v)-lo(3))*inv), i8)
            found = .false.
            search: do ddx = -1, 1
                do ddy = -1, 1
                    do ddz = -1, 1
                        h = cell_hash(cx+ddx, cy+ddy, cz+ddz, msize)
                        u = head(h)
                        do while (u /= 0)
                            if (vsolid(u) == vsolid(v)) then
                                if (maxval(abs(vertices(:,u) - vertices(:,v))) <= tol) then
                                    rep(v) = rep(u); found = .true.; exit search
                                end if
                            end if
                            u = nxt(u)
                        end do
                    end do
                end do
            end do search
            if (.not. found) then
                h = cell_hash(cx, cy, cz, msize)
                nxt(v) = head(h); head(h) = v
            end if
        end do
        deallocate(head, nxt, vsolid)
    end subroutine weld_vertices

    pure integer function cell_hash(ix, iy, iz, msize)
        integer(i8), intent(in) :: ix, iy, iz
        integer,     intent(in) :: msize
        integer(i8) :: h
        h = ieor(ieor(ix * 73856093_i8, iy * 19349663_i8), iz * 83492791_i8)
        cell_hash = int(iand(h, int(msize-1, i8)))
    end function cell_hash

    ! ------------------------------------------------------------------------------------------------------------------
    ! Triangles sharing all three welded vertices (within one solid): the net orientation decides what survives.
    ! n copies in agreement keep one; k copies one way and k the other cancel completely (an internal wall).
    subroutine cancel_duplicates(wf, face_solid_id, nfaces, nvp, active, n_dup)
        integer,     intent(in)    :: wf(:,:), face_solid_id(:), nfaces
        integer(i8), intent(in)    :: nvp
        logical,     intent(inout) :: active(:)
        integer,     intent(inout) :: n_dup(:)
        integer(i8), allocatable :: key(:)
        integer,     allocatable :: ord(:), tri(:), srt(:,:)
        integer :: f, m, i, j, p, q, fp, fq, net, even_f, odd_f, cnt
        logical, allocatable :: done(:)

        m = count(active)
        if (m == 0) return
        allocate(key(m), ord(m), tri(m), srt(3,m), done(m))
        i = 0
        do f = 1, nfaces
            if (.not. active(f)) cycle
            i = i + 1
            tri(i) = f
            srt(1,i) = min(wf(1,f), wf(2,f), wf(3,f))
            srt(3,i) = max(wf(1,f), wf(2,f), wf(3,f))
            srt(2,i) = wf(1,f) + wf(2,f) + wf(3,f) - srt(1,i) - srt(3,i)
            key(i)   = int(srt(1,i), i8) * nvp + int(srt(2,i), i8)
            ord(i)   = i
        end do
        done = .false.
        call sort_keys(key, ord, m)
        i = 1
        do while (i <= m)
            j = i
            do while (j < m)
                if (key(j+1) /= key(i)) exit
                j = j + 1
            end do
            if (j > i) then
                do p = i, j
                    if (done(p)) cycle
                    fp = tri(ord(p))
                    net = 0; cnt = 0; even_f = 0; odd_f = 0
                    do q = p, j
                        if (done(q)) cycle
                        fq = tri(ord(q))
                        if (srt(3,ord(q)) /= srt(3,ord(p)) .or. face_solid_id(fq) /= face_solid_id(fp)) cycle
                        done(q) = .true.
                        cnt = cnt + 1
                        if (is_cyclic(wf(:,fq), srt(:,ord(q)))) then
                            net = net + 1; if (even_f == 0) even_f = fq
                        else
                            net = net - 1; if (odd_f == 0) odd_f = fq
                        end if
                        active(fq) = .false.
                    end do
                    if (net > 0) active(even_f) = .true.
                    if (net < 0) active(odd_f)  = .true.
                    n_dup(face_solid_id(fp)) = n_dup(face_solid_id(fp)) + cnt - merge(1, 0, net /= 0)
                end do
            end if
            i = j + 1
        end do
        deallocate(key, ord, tri, srt, done)
    end subroutine cancel_duplicates

    ! stored vertex order is a cyclic rotation of the sorted order
    pure logical function is_cyclic(w, s)
        integer, intent(in) :: w(3), s(3)
        is_cyclic = (w(1)==s(1) .and. w(2)==s(2)) .or. (w(1)==s(2) .and. w(2)==s(3)) .or. (w(1)==s(3) .and. w(2)==s(1))
    end function is_cyclic

    ! ------------------------------------------------------------------------------------------------------------------
    subroutine print_report(vertices, nfaces, nsolids, is_ground, face_solid_id, active, n_deg, n_dup, n_bnd, n_nonman, &
                            n_incons, n_flip, n_patch, n_fixed, n_conflict, nbe, edge_a, edge_b, edge_s, edge_state, &
                            n_lids, n_lidtri, n_left, do_cap, report_dir)
        real(dp), intent(in) :: vertices(:,:)
        integer,  intent(in) :: nfaces, nsolids, face_solid_id(:), n_conflict, nbe
        logical,  intent(in) :: is_ground(:), active(:), do_cap
        integer,  intent(in) :: n_deg(:), n_dup(:), n_bnd(:), n_nonman(:), n_incons(:), n_flip(:), n_patch(:), n_fixed(:)
        integer,  intent(in) :: edge_a(:), edge_b(:), edge_s(:), edge_state(:), n_lids(:), n_lidtri(:), n_left(:)
        character(len=*), intent(in) :: report_dir
        integer, allocatable :: uf(:), eroot(:), eord(:), loop_n(:), uf_solid(:)
        integer(i8), allocatable :: ekey(:)
        real(dp), allocatable :: loop_lo(:,:), loop_hi(:,:), loop_len(:)
        integer, allocatable :: loop_state(:)
        integer :: s, e, i, j, nloop, best, shown, cnt_tri, unit, ios
        real(dp) :: len_e
        character(len=1024) :: line
        logical :: any_open

        print *, "--- Mesh check (welded, duplicates cancelled, orientation analysed) ---"
        any_open = .false.
        do s = 1, nsolids
            cnt_tri = count(face_solid_id(1:nfaces) == s)
            write(line,'(A,I0,A,I0,A,I0,A,I0,A,I0,A)') " solid ", s, ": ", cnt_tri, " triangles (", n_deg(s), &
                " degenerate, ", n_dup(s), " duplicate dropped), ", n_patch(s), " connected patch(es)"
            print '(A)', trim(line)
            if (is_ground(s)) then
                print '(A)', "    ground surface: open by design, classified as 'everything below the topmost sheet'"
                cycle
            end if
            if (n_bnd(s) == 0 .and. n_nonman(s) == 0 .and. n_incons(s) == 0) then
                print '(A)', "    watertight, consistently oriented"
            else
                if (n_bnd(s) > 0) then
                    write(line,'(A,I0,A)') "    OPEN: ", n_bnd(s), " boundary edge(s) (holes / cracks / missing faces)"
                    print '(A)', trim(line)
                    if (n_lids(s) > 0) then
                        write(line,'(A,I0,A,I0,A)') "    -> closed ", n_lids(s), " planar hole loop(s) with a lid of ", n_lidtri(s), " triangle(s)"
                        print '(A)', trim(line)
                    end if
                    if (n_left(s) > 0) then
                        any_open = .true.
                        write(line,'(A,I0,A)') "    -> ", n_left(s), " boundary edge(s) left open (seam/crack, non-planar or branching loop)"
                        print '(A)', trim(line)
                    else if (.not. do_cap) then
                        any_open = .true.
                    end if
                end if
                if (n_nonman(s) > 0) then
                    write(line,'(A,I0,A)') "    ", n_nonman(s), " non-manifold edge(s) (more than two triangles share an edge)"
                    print '(A)', trim(line)
                end if
                if (n_incons(s) > 0) then
                    write(line,'(A,I0,A,I0,A,I0,A)') "    ", n_incons(s), " inconsistently oriented edge(s) -> repaired in ", &
                        n_fixed(s), " patch(es), ", n_flip(s), " triangle(s) re-oriented"
                    print '(A)', trim(line)
                end if
            end if
        end do
        if (n_conflict > 0) print '(A,I0,A)', "    WARNING: ", n_conflict, " edge(s) where orientation cannot be made consistent (non-orientable)"

        ! ---- locate the holes: group boundary edges into loops (union-find on welded vertex ids) ----
        if (nbe > 0) then
            allocate(uf(size(vertices,2)))
            do i = 1, size(uf)
                uf(i) = i
            end do
            do e = 1, nbe
                call uf_union(uf, edge_a(e), edge_b(e))
            end do
            allocate(eroot(nbe), eord(nbe), ekey(nbe))
            do e = 1, nbe
                eroot(e) = uf_find(uf, edge_a(e))
                ekey(e)  = int(edge_s(e), i8) * (int(size(uf), i8) + 1_i8) + int(eroot(e), i8)
                eord(e)  = e
            end do
            call sort_keys(ekey, eord, nbe)
            nloop = 1
            do e = 2, nbe
                if (ekey(e) /= ekey(e-1)) nloop = nloop + 1
            end do
            allocate(loop_n(nloop), loop_lo(3,nloop), loop_hi(3,nloop), loop_len(nloop))
            allocate(uf_solid(nloop), loop_state(nloop))
            loop_state = 0
            loop_n = 0; loop_len = 0.0_dp; loop_lo = huge(1.0_dp); loop_hi = -huge(1.0_dp)
            j = 0
            do e = 1, nbe
                if (e == 1) then
                    j = 1
                else if (ekey(e) /= ekey(e-1)) then
                    j = j + 1
                end if
                i = eord(e)
                uf_solid(j) = edge_s(i); loop_state(j) = edge_state(i)
                loop_n(j) = loop_n(j) + 1
                loop_lo(:,j) = min(loop_lo(:,j), min(vertices(:,edge_a(i)), vertices(:,edge_b(i))))
                loop_hi(:,j) = max(loop_hi(:,j), max(vertices(:,edge_a(i)), vertices(:,edge_b(i))))
                len_e = sqrt(sum((vertices(:,edge_a(i)) - vertices(:,edge_b(i)))**2))
                loop_len(j) = loop_len(j) + len_e
            end do
            do s = 1, nsolids
                if (is_ground(s) .or. n_bnd(s) == 0) cycle
                shown = 0
                do
                    best = 0
                    do j = 1, nloop
                        if (uf_solid(j) == s .and. loop_n(j) > 0) then
                            if (best == 0) then
                                best = j
                            else if (loop_n(j) > loop_n(best)) then
                                best = j
                            end if
                        end if
                    end do
                    if (best == 0 .or. shown >= 6) exit
                    shown = shown + 1
                    call to_user_axes(loop_lo(:,best)); call to_user_axes(loop_hi(:,best))
                    print '(A,I0,A,I0,A,ES9.2,A,3ES11.3,A,3ES11.3,A)', "    solid ", s, " hole/crack: ", loop_n(best), &
                        " edges, perimeter ", loop_len(best), ", bbox ", min(loop_lo(:,best),loop_hi(:,best)), " ..", &
                        max(loop_lo(:,best),loop_hi(:,best)), merge("  [lid]      ", merge("  [seam]     ", "  [left open]", loop_state(best) == 3), loop_state(best) == 1)
                    loop_n(best) = 0
                end do
                shown = count(uf_solid == s .and. loop_n > 0)
                if (shown > 0) print '(A,I0,A)', "    ... and ", shown, " more open loop(s)"
            end do
            if (len_trim(report_dir) > 0) then
                open(newunit=unit, file=trim(report_dir)//'open_edges.obj', status='replace', action='write', iostat=ios)
                if (ios == 0) then
                    write(unit,'(A)') '# boundary edges of the input surface meshes (holes, cracks, missing faces)'
                    j = 0
                    do e = 1, min(nbe, 400000)
                        if (is_ground(edge_s(e))) cycle
                        write(unit,'(A,3ES16.8)') 'v ', user_axes(vertices(:,edge_a(e)))
                        write(unit,'(A,3ES16.8)') 'v ', user_axes(vertices(:,edge_b(e)))
                        write(unit,'(A,I0,A,I0)') 'l ', 2*j+1, ' ', 2*j+2
                        j = j + 1
                    end do
                    close(unit)
                    if (j > 0) print '(A)', "    open edges written to "//trim(report_dir)//"open_edges.obj"
                end if
            end if
            deallocate(uf, eroot, eord, ekey, loop_n, loop_lo, loop_hi, loop_len, uf_solid, loop_state)
        end if
        if (any_open) print '(A)', "    open shell(s) remain: the signed-ray integrity test and winding-number fallback decide solid/fluid."
        print *, "-------------------------------------------------------------------------"
    end subroutine print_report

    ! mesh coordinates are GenSDF-internal (vertical = axis 3); show them in the user's axes
    pure function user_axes(p) result(q)
        real(dp), intent(in) :: p(3)
        real(dp) :: q(3)
        q = p
        if (vertical_axis == 2) q = [p(1), p(3), p(2)]
    end function user_axes

    subroutine to_user_axes(p)
        real(dp), intent(inout) :: p(3)
        p = user_axes(p)
    end subroutine to_user_axes

    integer function uf_find(uf, x)
        integer, intent(inout) :: uf(:)
        integer, intent(in) :: x
        integer :: r, y, nxt
        r = x
        do while (uf(r) /= r)
            r = uf(r)
        end do
        y = x
        do while (uf(y) /= r)
            nxt = uf(y); uf(y) = r; y = nxt
        end do
        uf_find = r
    end function uf_find

    subroutine uf_union(uf, a, b)
        integer, intent(inout) :: uf(:)
        integer, intent(in) :: a, b
        integer :: ra, rb
        ra = uf_find(uf, a); rb = uf_find(uf, b)
        if (ra /= rb) uf(max(ra,rb)) = min(ra,rb)
    end subroutine uf_union


    ! ------------------------------------------------------------------------------------------------------------------
    ! Lids over boundary loops. Boundary edge e is a directed edge a->b of triangle edge_f(e) (after the orientation
    ! repair, reversed if that triangle is flipped); a lid must traverse it b->a. Lid edges chain head-to-tail into loops.
    ! A loop is capped when it is simple (no vertex repeated, no branching), has >= 3 vertices, is nearly planar
    ! (max out-of-plane distance <= planar_tol * rms radius) and is not part of a seam. Planar loops are triangulated by
    ! ear clipping, so concave outlines (L-shaped floors, window frames) get a lid that stays inside the outline.
    ! Seam: a loop with a vertex within crack_dist of another loop of the same solid is the edge of a crack or of two
    ! patches that merely do not share vertices. Lids on both sides of such a crack would be two independent plates
    ! that cross each other; the winding number of the open surface handles a narrow crack correctly, so both stay open.
    subroutine cap_boundary_loops(vertices, nbe, edge_a, edge_b, edge_s, edge_f, flip, is_ground, nsolids, nv, crack_dist, &
                                  cap_faces, cap_solid, edge_state, n_lids, n_lidtri)
        real(dp), intent(in)  :: vertices(:,:), crack_dist
        integer,  intent(in)  :: nbe, edge_a(:), edge_b(:), edge_s(:), edge_f(:), nsolids, nv
        logical,  intent(in)  :: flip(:), is_ground(:)
        integer,  allocatable, intent(inout) :: cap_faces(:,:), cap_solid(:)
        integer,  intent(inout) :: edge_state(:), n_lids(:), n_lidtri(:)
        real(dp), parameter :: planar_tol = 0.2_dp
        integer,  parameter :: ear_clip_max = 2500
        integer, allocatable :: c_start(:), c_end(:), ord(:), tri(:,:), stamp(:), tmp_f(:,:), tmp_s(:)
        integer, allocatable :: lv(:), le(:), lstart(:), lcount(:), lsolid(:), vloop(:), nxt_hash(:), head_hash(:)
        integer(i8), allocatable :: skey(:)
        logical, allocatable :: visited(:), lok(:), lcrack(:)
        real(dp), allocatable :: pts(:,:)
        integer :: e, i, j, k, n, lo, hi, cur, nxt, v, u, ncap, ntri, s, stampc, nloop, nlv, l, msize, h
        integer :: ddx, ddy, ddz
        integer(i8) :: cx, cy, cz
        logical :: ok, planar
        real(dp) :: nrm(3), rms, dev, ctr(3), lo3(3), inv

        allocate(c_start(nbe), c_end(nbe), ord(nbe), skey(nbe), visited(nbe), stamp(nv))
        allocate(lv(nbe), le(nbe), lstart(nbe+1), lcount(nbe), lsolid(nbe), lok(nbe), lcrack(nbe))
        do e = 1, nbe
            if (flip(edge_f(e))) then
                c_start(e) = edge_a(e); c_end(e) = edge_b(e)
            else
                c_start(e) = edge_b(e); c_end(e) = edge_a(e)
            end if
            skey(e) = int(c_start(e), i8); ord(e) = e
        end do
        call sort_keys(skey, ord, nbe)
        visited = .false.; stamp = 0; stampc = 0
        ncap = size(cap_faces, 2)

        ! ---- pass 1: trace the loops (vertex list lv, edge list le, per-loop offsets) ----
        nloop = 0; nlv = 0
        do e = 1, nbe
            if (visited(e) .or. is_ground(edge_s(e))) cycle
            nloop = nloop + 1
            lstart(nloop) = nlv + 1
            n = 0; cur = e; ok = .true.
            do
                visited(cur) = .true.
                n = n + 1; lv(nlv+n) = c_start(cur); le(nlv+n) = cur
                v = c_end(cur)
                if (v == c_start(e)) exit
                call edges_from(v, nbe, skey, ord, lo, hi)
                if (hi < lo) then
                    ok = .false.; exit                      ! open chain (should not happen)
                end if
                if (hi > lo) ok = .false.                   ! branching vertex: not a simple loop
                nxt = 0
                do k = lo, hi
                    if (.not. visited(ord(k))) then; nxt = ord(k); exit; end if
                end do
                if (nxt == 0) then
                    ok = .false.; exit
                end if
                cur = nxt
                if (n > nbe) then; ok = .false.; exit; end if
            end do
            if (ok .and. n < 3) ok = .false.
            if (ok) then       ! vertices must be distinct
                stampc = stampc + 1
                do i = 1, n
                    if (stamp(lv(nlv+i)) == stampc) then; ok = .false.; exit; end if
                    stamp(lv(nlv+i)) = stampc
                end do
            end if
            lcount(nloop) = n; lsolid(nloop) = edge_s(e); lok(nloop) = ok
            nlv = nlv + n
        end do
        lstart(nloop+1) = nlv + 1

        ! ---- pass 2: seams. Hash every loop vertex; a vertex of another loop of the same solid within crack_dist
        !      marks both loops ----
        lcrack = .false.
        if (nloop > 1 .and. crack_dist > 0.0_dp) then
            allocate(vloop(nlv), nxt_hash(nlv))
            msize = 1
            do while (msize < 2*nlv); msize = msize * 2; end do
            allocate(head_hash(0:msize-1))
            head_hash = 0; nxt_hash = 0
            lo3 = [minval(vertices(1,:)), minval(vertices(2,:)), minval(vertices(3,:))]
            inv = 1.0_dp / crack_dist
            do l = 1, nloop
                do i = lstart(l), lstart(l+1) - 1
                    vloop(i) = l
                end do
            end do
            do i = 1, nlv
                v = lv(i)
                cx = int(floor((vertices(1,v)-lo3(1))*inv), i8)
                cy = int(floor((vertices(2,v)-lo3(2))*inv), i8)
                cz = int(floor((vertices(3,v)-lo3(3))*inv), i8)
                do ddx = -1, 1
                    do ddy = -1, 1
                        do ddz = -1, 1
                            h = cell_hash(cx+ddx, cy+ddy, cz+ddz, msize)
                            j = head_hash(h)
                            do while (j /= 0)
                                if (vloop(j) /= vloop(i) .and. lsolid(vloop(j)) == lsolid(vloop(i))) then
                                    u = lv(j)
                                    if (maxval(abs(vertices(:,u) - vertices(:,v))) <= crack_dist) then
                                        lcrack(vloop(i)) = .true.; lcrack(vloop(j)) = .true.
                                    end if
                                end if
                                j = nxt_hash(j)
                            end do
                        end do
                    end do
                end do
                h = cell_hash(cx, cy, cz, msize)
                nxt_hash(i) = head_hash(h); head_hash(h) = i
            end do
            deallocate(vloop, nxt_hash, head_hash)
        end if

        ! ---- pass 3: cap what is simple, planar and not a seam ----
        do l = 1, nloop
            n = lcount(l)
            ok = lok(l) .and. .not. lcrack(l)
            if (ok) then
                allocate(pts(3,n))
                do i = 1, n
                    pts(:,i) = vertices(:,lv(lstart(l)+i-1))
                end do
                call plane_fit(pts, n, ctr, nrm, rms, dev, planar)
                ok = planar .and. dev <= planar_tol * rms
                if (ok) then
                    allocate(tri(3, n-2))
                    call triangulate_polygon(pts, n, ctr, nrm, ear_clip_max, tri, ntri)
                    s = lsolid(l)
                    if (ncap + ntri > size(cap_faces, 2)) then
                        allocate(tmp_f(3, max(2*(ncap+ntri), 64)), tmp_s(max(2*(ncap+ntri), 64)))
                        tmp_f(:,1:ncap) = cap_faces(:,1:ncap); tmp_s(1:ncap) = cap_solid(1:ncap)
                        call move_alloc(tmp_f, cap_faces); call move_alloc(tmp_s, cap_solid)
                    end if
                    do i = 1, ntri
                        cap_faces(:, ncap+i) = [lv(lstart(l)+tri(1,i)-1), lv(lstart(l)+tri(2,i)-1), lv(lstart(l)+tri(3,i)-1)]
                        cap_solid(ncap+i) = s
                    end do
                    ncap = ncap + ntri
                    n_lids(s) = n_lids(s) + 1; n_lidtri(s) = n_lidtri(s) + ntri
                    deallocate(tri)
                end if
                deallocate(pts)
            end if
            do i = 1, n
                edge_state(le(lstart(l)+i-1)) = merge(1, merge(3, 2, lcrack(l)), ok)
            end do
        end do
        ! trim the lid arrays to the used size
        allocate(tmp_f(3, ncap), tmp_s(ncap))
        if (ncap > 0) then
            tmp_f = cap_faces(:,1:ncap); tmp_s = cap_solid(1:ncap)
        end if
        call move_alloc(tmp_f, cap_faces); call move_alloc(tmp_s, cap_solid)
        deallocate(c_start, c_end, ord, skey, visited, stamp, lv, le, lstart, lcount, lsolid, lok, lcrack)
    end subroutine cap_boundary_loops

    ! index range [lo,hi] in the start-vertex-sorted order of the lid edges that begin at vertex v
    subroutine edges_from(v, nbe, skey, ord, lo, hi)
        integer,     intent(in)  :: v, nbe
        integer(i8), intent(in)  :: skey(:)
        integer,     intent(in)  :: ord(:)
        integer,     intent(out) :: lo, hi
        integer :: a, b, mid
        a = 1; b = nbe + 1                 ! lower bound
        do while (a < b)
            mid = (a + b) / 2
            if (skey(mid) < int(v, i8)) then; a = mid + 1; else; b = mid; end if
        end do
        lo = a
        a = lo; b = nbe + 1                ! upper bound
        do while (a < b)
            mid = (a + b) / 2
            if (skey(mid) <= int(v, i8)) then; a = mid + 1; else; b = mid; end if
        end do
        hi = a - 1
    end subroutine edges_from

    ! Newell normal (loop orientation), centroid, rms radius and max out-of-plane deviation of a closed polyline
    subroutine plane_fit(p, n, ctr, nrm, rms, dev, ok)
        real(dp), intent(in)  :: p(:,:)
        integer,  intent(in)  :: n
        real(dp), intent(out) :: ctr(3), nrm(3), rms, dev
        logical,  intent(out) :: ok
        integer :: i, j
        real(dp) :: q(3), r(3), m
        ctr = sum(p(:,1:n), dim=2) / real(n, dp)
        nrm = 0.0_dp
        do i = 1, n
            j = mod(i, n) + 1
            q = p(:,i) - ctr; r = p(:,j) - ctr
            nrm = nrm + [q(2)*r(3)-q(3)*r(2), q(3)*r(1)-q(1)*r(3), q(1)*r(2)-q(2)*r(1)]
        end do
        m = sqrt(dot_product(nrm, nrm))
        ok = (m > 0.0_dp)
        rms = 0.0_dp; dev = 0.0_dp
        if (.not. ok) return
        nrm = nrm / m
        do i = 1, n
            q = p(:,i) - ctr
            rms = rms + dot_product(q, q)
            dev = max(dev, abs(dot_product(q, nrm)))
        end do
        rms = sqrt(rms / real(n, dp))
    end subroutine plane_fit

    ! Ear-clipping triangulation of a planar polygon given in counter-clockwise order about nrm. tri holds local vertex
    ! indices; the triangles inherit the polygon's orientation. For very long loops (> ear_max) a fan from vertex 1 is used.
    subroutine triangulate_polygon(p, n, ctr, nrm, ear_max, tri, ntri)
        real(dp), intent(in)  :: p(:,:), ctr(3), nrm(3)
        integer,  intent(in)  :: n, ear_max
        integer,  intent(out) :: tri(:,:), ntri
        real(dp), allocatable :: x(:), y(:)
        integer,  allocatable :: prv(:), nxt(:)
        real(dp) :: u(3), v(3), a(3), eps, best, cr
        integer :: i, j, t, nrem, ia, ib, ic, start, pick, axis
        logical :: found, blocked

        ntri = 0
        if (n < 3) return
        if (n > ear_max) then
            do i = 2, n - 1
                ntri = ntri + 1; tri(:,ntri) = [1, i, i+1]
            end do
            return
        end if
        axis = minloc(abs(nrm), dim=1)
        a = 0.0_dp; a(axis) = 1.0_dp
        u = [nrm(2)*a(3)-nrm(3)*a(2), nrm(3)*a(1)-nrm(1)*a(3), nrm(1)*a(2)-nrm(2)*a(1)]
        u = u / sqrt(dot_product(u,u))
        v = [nrm(2)*u(3)-nrm(3)*u(2), nrm(3)*u(1)-nrm(1)*u(3), nrm(1)*u(2)-nrm(2)*u(1)]
        allocate(x(n), y(n), prv(n), nxt(n))
        do i = 1, n
            x(i) = dot_product(p(:,i) - ctr, u); y(i) = dot_product(p(:,i) - ctr, v)
            prv(i) = mod(i + n - 2, n) + 1; nxt(i) = mod(i, n) + 1
        end do
        eps = 1.0e-12_dp * max(maxval(x*x + y*y), tiny(1.0_dp))
        nrem = n; start = 1
        do while (nrem > 3)
            found = .false.; best = -huge(1.0_dp); pick = start
            i = start
            do t = 1, nrem
                ia = prv(i); ib = i; ic = nxt(i)
                cr = (x(ib)-x(ia))*(y(ic)-y(ia)) - (y(ib)-y(ia))*(x(ic)-x(ia))
                if (cr > best) then; best = cr; pick = i; end if
                if (cr > eps) then
                    blocked = .false.
                    j = nxt(ic)
                    do while (j /= ia)
                        if (in_triangle(x(j), y(j), x(ia), y(ia), x(ib), y(ib), x(ic), y(ic))) then
                            blocked = .true.; exit
                        end if
                        j = nxt(j)
                    end do
                    if (.not. blocked) then; found = .true.; pick = i; exit; end if
                end if
                i = nxt(i)
            end do
            ia = prv(pick); ib = pick; ic = nxt(pick)
            cr = (x(ib)-x(ia))*(y(ic)-y(ia)) - (y(ib)-y(ia))*(x(ic)-x(ia))
            if (found .or. cr > eps) then
                ntri = ntri + 1; tri(:,ntri) = [ia, ib, ic]
            end if
            nxt(ia) = ic; prv(ic) = ia
            start = ic
            nrem = nrem - 1
        end do
        ia = start; ib = nxt(ia); ic = nxt(ib)
        ntri = ntri + 1; tri(:,ntri) = [ia, ib, ic]
        deallocate(x, y, prv, nxt)
    end subroutine triangulate_polygon

    ! point inside or on the boundary of the counter-clockwise triangle (a,b,c)
    pure logical function in_triangle(px, py, ax, ay, bx, by, cx, cy)
        real(dp), intent(in) :: px, py, ax, ay, bx, by, cx, cy
        real(dp) :: d1, d2, d3
        d1 = (bx-ax)*(py-ay) - (by-ay)*(px-ax)
        d2 = (cx-bx)*(py-by) - (cy-by)*(px-bx)
        d3 = (ax-cx)*(py-cy) - (ay-cy)*(px-cx)
        in_triangle = (d1 >= 0.0_dp .and. d2 >= 0.0_dp .and. d3 >= 0.0_dp)
    end function in_triangle

    ! ------------------------------------------------------------------------------------------------------------------
    ! Bottom-up merge sort of int64 keys carrying an integer payload.
    subroutine sort_keys(key, idx, n)
        integer(i8), intent(inout) :: key(:)
        integer,     intent(inout) :: idx(:)
        integer,     intent(in)    :: n
        integer(i8), allocatable :: k2(:)
        integer,     allocatable :: i2(:)
        integer :: width, lo, mid, hi, i, j, o
        logical :: in_main
        if (n < 2) return
        allocate(k2(n), i2(n))
        width = 1; in_main = .true.
        do while (width < n)
            do lo = 1, n, 2*width
                mid = min(lo + width - 1, n); hi = min(lo + 2*width - 1, n)
                i = lo; j = mid + 1; o = lo
                if (in_main) then
                    do while (i <= mid .and. j <= hi)
                        if (key(i) <= key(j)) then
                            k2(o) = key(i); i2(o) = idx(i); i = i + 1
                        else
                            k2(o) = key(j); i2(o) = idx(j); j = j + 1
                        end if
                        o = o + 1
                    end do
                    do while (i <= mid); k2(o) = key(i); i2(o) = idx(i); i = i + 1; o = o + 1; end do
                    do while (j <= hi);  k2(o) = key(j); i2(o) = idx(j); j = j + 1; o = o + 1; end do
                else
                    do while (i <= mid .and. j <= hi)
                        if (k2(i) <= k2(j)) then
                            key(o) = k2(i); idx(o) = i2(i); i = i + 1
                        else
                            key(o) = k2(j); idx(o) = i2(j); j = j + 1
                        end if
                        o = o + 1
                    end do
                    do while (i <= mid); key(o) = k2(i); idx(o) = i2(i); i = i + 1; o = o + 1; end do
                    do while (j <= hi);  key(o) = k2(j); idx(o) = i2(j); j = j + 1; o = o + 1; end do
                end if
            end do
            in_main = .not. in_main
            width = 2 * width
        end do
        if (.not. in_main) then
            key(1:n) = k2(1:n); idx(1:n) = i2(1:n)
        end if
        deallocate(k2, i2)
    end subroutine sort_keys

end module meshprep_mod
