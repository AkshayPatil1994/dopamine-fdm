module winding_mod
    !
    ! winding_mod — generalized winding number (Jacobson, Kavan & Sorkine-Hornung 2013) of a triangle mesh, evaluated
    ! with a bounding-volume hierarchy and a dipole far-field approximation (Barill et al. 2018, zeroth order).
    !
    !   w(p) = (1/4pi) * sum_triangles  SolidAngle(triangle seen from p)
    !
    ! For a closed, consistently oriented surface w is 1 inside and 0 outside (-1/0 if the orientation is inverted;
    ! 2 where two shells overlap). For an OPEN surface w is a harmonic function that varies smoothly across the holes
    ! (inside a box with one face missing it is 1 - Omega_hole/4pi), so thresholding at 1/2 still separates inside
    ! from outside — unlike ray-crossing parity, which is flipped for the whole remainder of any ray through a hole.
    ! Triangle soups, cracks, T-junctions and duplicated or overlapping patches degrade gracefully instead of
    ! catastrophically. Cost per point is O(log N) with the hierarchy.
    !
    use utils_io, only : dp
    implicit none
    private
    public :: bvh_t, bvh_build, bvh_free, winding_number

    real(dp), parameter :: pi   = 3.14159265358979323846_dp
    ! A node is replaced by its dipole when the query is farther than beta * (node radius). Leaves hold <= leaf_size triangles.
    real(dp), parameter :: beta = 3.0_dp
    integer,  parameter :: leaf_size = 6

    type :: bvh_t
        integer :: nnode = 0, ntri = 0
        integer,  allocatable :: left(:), right(:)      ! children (0 for a leaf)
        integer,  allocatable :: first(:), cnt(:)       ! triangle range in tid
        integer,  allocatable :: tid(:)                 ! triangle ids, permuted so every node is a contiguous range
        real(dp), allocatable :: cen(:,:)               ! node centre (area-weighted centroid)
        real(dp), allocatable :: avec(:,:)              ! node area vector: sum of 0.5*(e1 x e2)
        real(dp), allocatable :: rad(:)                 ! radius of the sphere about cen enclosing the node's vertices
    end type bvh_t

contains

    subroutine bvh_free(b)
        type(bvh_t), intent(inout) :: b
        if (allocated(b%left)) deallocate(b%left, b%right, b%first, b%cnt, b%tid, b%cen, b%avec, b%rad)
        b%nnode = 0; b%ntri = 0
    end subroutine bvh_free

    ! Build over the triangles listed in tri_list (indices into faces). faces must already carry the repaired orientation.
    subroutine bvh_build(b, vertices, faces, tri_list)
        type(bvh_t), intent(out) :: b
        real(dp), intent(in) :: vertices(:,:)
        integer,  intent(in) :: faces(:,:), tri_list(:)
        real(dp), allocatable :: tc(:,:)
        integer :: n, i, node, lo, hi, mid, axis, maxnode, sp
        integer :: stk_node(128), stk_lo(128), stk_hi(128)
        real(dp) :: bmin(3), bmax(3), ext(3)

        n = size(tri_list)
        b%ntri = n
        if (n == 0) return
        maxnode = max(1, 4 * ((n + leaf_size - 1) / leaf_size) + 4)
        allocate(b%left(maxnode), b%right(maxnode), b%first(maxnode), b%cnt(maxnode))
        allocate(b%cen(3,maxnode), b%avec(3,maxnode), b%rad(maxnode), b%tid(n), tc(3,n))
        b%left = 0; b%right = 0
        b%tid = tri_list
        do i = 1, n
            tc(:,i) = (vertices(:,faces(1,b%tid(i))) + vertices(:,faces(2,b%tid(i))) + vertices(:,faces(3,b%tid(i)))) / 3.0_dp
        end do

        b%nnode = 1
        sp = 1; stk_node(1) = 1; stk_lo(1) = 1; stk_hi(1) = n
        do while (sp > 0)
            node = stk_node(sp); lo = stk_lo(sp); hi = stk_hi(sp); sp = sp - 1
            b%first(node) = lo; b%cnt(node) = hi - lo + 1
            call node_moments(b, vertices, faces, node, lo, hi)
            if (hi - lo + 1 <= leaf_size) cycle
            bmin = minval(tc(:,lo:hi), dim=2); bmax = maxval(tc(:,lo:hi), dim=2)
            ext = bmax - bmin
            axis = maxloc(ext, dim=1)
            if (ext(axis) <= 0.0_dp) cycle                 ! coincident centroids: keep as a (large) leaf
            mid = (lo + hi) / 2
            call select_nth(b%tid, tc, lo, hi, mid, axis)
            b%left(node)  = b%nnode + 1
            b%right(node) = b%nnode + 2
            b%nnode = b%nnode + 2
            sp = sp + 1; stk_node(sp) = b%right(node); stk_lo(sp) = mid + 1; stk_hi(sp) = hi
            sp = sp + 1; stk_node(sp) = b%left(node);  stk_lo(sp) = lo;      stk_hi(sp) = mid
        end do
        deallocate(tc)
    end subroutine bvh_build

    ! area vector, area-weighted centroid and enclosing radius of the triangles tid(lo:hi)
    subroutine node_moments(b, vertices, faces, node, lo, hi)
        type(bvh_t), intent(inout) :: b
        real(dp), intent(in) :: vertices(:,:)
        integer,  intent(in) :: faces(:,:), node, lo, hi
        real(dp) :: av(3), ce(3), e1(3), e2(3), wsum, w, c0(3), r2, d(3)
        integer :: i, t, k

        b%avec(:,node) = 0.0_dp; c0 = 0.0_dp; wsum = 0.0_dp
        do i = lo, hi
            t = b%tid(i)
            e1 = vertices(:,faces(2,t)) - vertices(:,faces(1,t))
            e2 = vertices(:,faces(3,t)) - vertices(:,faces(1,t))
            av = 0.5_dp * [e1(2)*e2(3)-e1(3)*e2(2), e1(3)*e2(1)-e1(1)*e2(3), e1(1)*e2(2)-e1(2)*e2(1)]
            ce = (vertices(:,faces(1,t)) + vertices(:,faces(2,t)) + vertices(:,faces(3,t))) / 3.0_dp
            w  = sqrt(dot_product(av, av))
            b%avec(:,node) = b%avec(:,node) + av
            c0 = c0 + w * ce
            wsum = wsum + w
        end do
        if (wsum > 0.0_dp) then
            c0 = c0 / wsum
        else
            c0 = (vertices(:,faces(1,b%tid(lo))) + vertices(:,faces(2,b%tid(lo))) + vertices(:,faces(3,b%tid(lo)))) / 3.0_dp
        end if
        b%cen(:,node) = c0
        r2 = 0.0_dp
        do i = lo, hi
            t = b%tid(i)
            do k = 1, 3
                d = vertices(:,faces(k,t)) - c0
                r2 = max(r2, dot_product(d, d))
            end do
        end do
        b%rad(node) = sqrt(r2)
    end subroutine node_moments

    ! Quickselect on tc(axis,lo:hi): afterwards position nth holds the median, smaller on its left, larger on its right.
    subroutine select_nth(tid, tc, lo0, hi0, nth, axis)
        integer,  intent(inout) :: tid(:)
        real(dp), intent(inout) :: tc(:,:)
        integer,  intent(in) :: lo0, hi0, nth, axis
        integer :: lo, hi, i, j, m
        real(dp) :: pivot
        lo = lo0; hi = hi0
        do while (lo < hi)
            m = (lo + hi) / 2
            pivot = median3(tc(axis,lo), tc(axis,m), tc(axis,hi))
            i = lo; j = hi
            do
                do while (tc(axis,i) < pivot); i = i + 1; end do
                do while (tc(axis,j) > pivot); j = j - 1; end do
                if (i <= j) then
                    call swap_item(tid, tc, i, j)
                    i = i + 1; j = j - 1
                end if
                if (i > j) exit
            end do
            if (nth <= j) then
                hi = j
            else if (nth >= i) then
                lo = i
            else
                exit
            end if
        end do
    end subroutine select_nth

    pure real(dp) function median3(a, b, c)
        real(dp), intent(in) :: a, b, c
        median3 = max(min(a,b), min(max(a,b), c))
    end function median3

    subroutine swap_item(tid, tc, i, j)
        integer,  intent(inout) :: tid(:)
        real(dp), intent(inout) :: tc(:,:)
        integer,  intent(in) :: i, j
        integer  :: it
        real(dp) :: tmp(3)
        it = tid(i); tid(i) = tid(j); tid(j) = it
        tmp = tc(:,i); tc(:,i) = tc(:,j); tc(:,j) = tmp
    end subroutine swap_item

    ! ------------------------------------------------------------------------------------------------------------------
    ! Winding number at p. Positive for points inside a surface whose triangles are wound counter-clockwise seen from
    ! outside.
    function winding_number(b, vertices, faces, p) result(w)
        type(bvh_t), intent(in) :: b
        real(dp), intent(in) :: vertices(:,:), p(3)
        integer,  intent(in) :: faces(:,:)
        real(dp) :: w, d(3), d2, r
        integer :: stk(128), sp, node, i, t

        w = 0.0_dp
        if (b%nnode == 0) return
        sp = 1; stk(1) = 1
        do while (sp > 0)
            node = stk(sp); sp = sp - 1
            d  = b%cen(:,node) - p
            d2 = dot_product(d, d)
            r  = beta * b%rad(node)
            if (d2 > r*r) then
                ! far field: dipole  A.(c-p)/|c-p|^3
                w = w + dot_product(b%avec(:,node), d) / (d2 * sqrt(d2))
            else if (b%left(node) == 0) then
                do i = b%first(node), b%first(node) + b%cnt(node) - 1
                    t = b%tid(i)
                    w = w + 2.0_dp * solid_angle_half(vertices(:,faces(1,t)) - p, vertices(:,faces(2,t)) - p, &
                                                      vertices(:,faces(3,t)) - p)
                end do
            else
                sp = sp + 1; stk(sp) = b%right(node)
                sp = sp + 1; stk(sp) = b%left(node)
            end if
        end do
        w = w / (4.0_dp * pi)
    end function winding_number

    ! Half the signed solid angle of the triangle (a,b,c) (vertex vectors relative to the query point), Van Oosterom &
    ! Strackee: tan(Omega/2) = a.(b x c) / (|a||b||c| + (a.b)|c| + (a.c)|b| + (b.c)|a|).
    ! For a triangle wound counter-clockwise seen from outside the solid and a query point inside, a.(b x c) > 0.
    pure function solid_angle_half(a, b, c) result(om2)
        real(dp), intent(in) :: a(3), b(3), c(3)
        real(dp) :: om2, la, lb, lc, num, den
        la = sqrt(dot_product(a,a)); lb = sqrt(dot_product(b,b)); lc = sqrt(dot_product(c,c))
        num = a(1)*(b(2)*c(3)-b(3)*c(2)) + a(2)*(b(3)*c(1)-b(1)*c(3)) + a(3)*(b(1)*c(2)-b(2)*c(1))
        den = la*lb*lc + dot_product(a,b)*lc + dot_product(a,c)*lb + dot_product(b,c)*la
        om2 = atan2(num, den)
    end function solid_angle_half

end module winding_mod
