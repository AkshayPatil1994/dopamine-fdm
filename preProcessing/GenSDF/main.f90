program generatesdf
    !
    ! GenSDF — Signed-Distance-Field Generator
    !
    use utils_io
    use narrowband_mod,   only : compute_scalar_distance_face, compute_narrowband_sdf
    use raycast_mod,      only : raycast_classify_domain, raycast_set_mesh
    use meshprep_mod,     only : prepare_mesh
    use spatial_hash_mod, only : build_hash, free_hash
    use fast_sweep_mod,   only : fast_sweep_3d
    use mpi
#ifdef GPU_SDF
    use openacc
#endif

    implicit none

    !  geometry
    integer :: nfaces, nvertices, nnormals, nsolids
    real(dp), allocatable :: vertices(:,:), normals(:,:)
    integer,  allocatable :: faces(:,:), face_normals(:,:)
    integer,  allocatable :: face_solid_id(:)
    real(dp) :: bbox_min(3), bbox_max(3)
    integer  :: sx, ex, sy, ey, sz, ez
    integer  :: cx_lo, cx_hi   ! per-rank x-slab intersected with the AABB+buffer [sx,ex]

    !  SDF arrays
    real(dp), allocatable :: sdf(:,:,:)
    real(dp), allocatable :: objid(:,:,:)
    real(dp), allocatable :: flood_fill_arr(:,:,:)   ! assembled field, root rank only (for the write)
    real(dp), allocatable :: objid_fill_arr(:,:,:)

    !  MPI neighbours for FSM 
    integer :: left_rank, right_rank

    !  MPI 
    integer :: myid, ierror, mpi_dx
    integer, allocatable :: decomp_x_start(:), decomp_x_end(:), decomp_size(:)

    !  timing 
    real(dp) :: startTime, endTime, totalTime, time1, time2

    !  misc
    integer :: ii
    logical :: debug = .false.

    !  mesh analysis (watertightness, orientation repair), broadcast from rank 0
    logical,  allocatable :: mesh_flip(:), mesh_active(:), mesh_open(:)
    integer,  allocatable :: cap_faces(:,:), cap_solid(:)
    integer :: ncap
    real(dp) :: weld_tol

#ifdef GPU_SDF
    !  multi-GPU device binding
    integer :: node_comm, node_rank, num_gpus
#endif

    ! MPI init and domain decomposition
    call MPI_INIT(ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid,   ierror)
    allocate(decomp_x_start(0:nprocs-1), decomp_x_end(0:nprocs-1), decomp_size(0:nprocs-1))

#ifdef GPU_SDF
    ! One GPU per rank (node-local rank mod GPU count).
    call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, node_comm, ierror)
    call MPI_COMM_RANK(node_comm, node_rank, ierror)
    num_gpus = acc_get_num_devices(acc_device_nvidia)
    if (num_gpus < 1) error stop "GPU_SDF: no NVIDIA GPU visible to this rank"
    call acc_set_device_num(mod(node_rank, num_gpus), acc_device_nvidia)
    call MPI_COMM_FREE(node_comm, ierror)
#endif

    if (myid == 0) then
        startTime = 0.0_dp
        call cpu_time(startTime)
        call printlogo()
        print *, "*** Starting with", nprocs, "MPI ranks ***"
    end if
    call MPI_BARRIER(MPI_COMM_WORLD, ierror)

    call read_inputfile(myid)
#ifdef GPU_SDF
    !$acc update device(scalarvalue)
#endif
    call read_cans_grid(myid, 'data/', dp, ng, r0, non_uniform_grid, xp, yp, zp, xf, yf, zf, dz, lx, ly, lz)
    allocate(dz_inverse(ng(3)))
    call setup_grid_spacing(myid, xf, yf, zf, ng(3), dx, dy, dz, dx_inverse, dy_inverse, dz_inverse)
    ! dy above is the uniform-grid spacing; with a stretched spanwise grid keep the smallest cell width as the
    !   scalar dy
    if (y_stretched) then
        dy = minval(dyv)
        dy_inverse = 1.0_dp / dy
    end if
#ifdef GPU_SDF
    ! dz is allocated via dummy-argument aliasing in read_cans_grid, so declare create's automatic on-allocate
    !   hook never fires for it.
    !$acc enter data create(dz)
    !$acc enter data create(dyv)
    !$acc update device(dx, dy, dz, dyv)
#endif

    ! Load geometry: a single .stl/.obj file (one solid, ID 1) or a .list manifest
    ! of multiple .stl/.obj files, each assigned a solid ID by its manifest order.
    call load_multi_object_geometry(myid, trim(inputfilename), vertices, normals, faces, face_normals, &
                                    face_solid_id, nvertices, nnormals, nfaces, nsolids)

    ! Solver clamps solid IDs into src/global.f90's max_ibm_objects-sized BC arrays; keep both in sync.
    if (nsolids > 15) then
        print *, 'FATAL: manifest has ', nsolids, ' solids, exceeds solver max_ibm_objects=15'
        error stop
    end if

    call check_object_resolution(myid, vertices, faces, face_solid_id, nfaces, nsolids)

    ! Mesh health: weld, drop duplicates, locate holes, repair triangle orientation, lid planar holes. Every
    !   rank runs
    ! the same deterministic analysis (only rank 0 reports), so no communication is needed.
    weld_tol = max(1.0e-3_dp * min(dx, dy, minval(dz)), 1.0e-9_dp * (maxval(vertices) - minval(vertices)))
    allocate(mesh_open(nsolids))
    call prepare_mesh(myid, vertices, faces, face_solid_id, nfaces, nsolids, solid_is_ground, weld_tol, &
                      2.0_dp * max(dx, dy, maxval(dz)), cap_holes /= 0, &
                      mesh_flip, mesh_active, mesh_open, cap_faces, cap_solid, 'data/')
    ncap = size(cap_solid)
    if (ncap > 0) call append_lids()
    if (sign_method == 0 .and. myid == 0 .and. any(mesh_open)) &
        print *, "WARNING: sign_method = 0 (parity) on an open mesh: expect misclassified regions behind the holes."

#ifdef GPU_SDF
    ! Geometry is read-only and reused unchanged across the p/u/v/w passes below.
    !$acc enter data copyin(vertices, normals, faces, face_normals, face_solid_id)
#endif

    call getbbox(myid, vertices, nvertices, bbox_min, bbox_max)
    call tagminmax(myid, xf, yp, zp, bbox_min, bbox_max, nx, ny, nz, dx, dy, dz(2), buffer_points, &
                   sx, ex, sy, ey, sz, ez, dyv, dz)

    !  domain decomposition: x-slabs over the FULL domain, not just the AABB, so fast_sweep_3d can reach
    !    x-cells outside the
    !  AABB+buffer box (otherwise they stay at the sentinel value and the field has a cliff at the AABB
    !    boundary). Only the
    !  near-surface computation is restricted to the AABB (cx_lo/cx_hi below).
    mpi_dx = ceiling(real(nx, dp) / real(nprocs, dp))
    decomp_x_start(0) = 1
    decomp_x_end(0)   = min(mpi_dx, nx)
    decomp_size(0)     = decomp_x_end(0) - decomp_x_start(0) + 1
    do ii = 1, nprocs-2
        decomp_x_start(ii) = decomp_x_end(ii-1) + 1
        decomp_x_end(ii)   = min(decomp_x_start(ii) + mpi_dx - 1, nx)
        decomp_size(ii)     = decomp_x_end(ii) - decomp_x_start(ii) + 1
    end do
    if (nprocs > 1) then
        decomp_x_start(nprocs-1) = decomp_x_end(nprocs-2) + 1
        decomp_x_end(nprocs-1)   = nx
        decomp_size(nprocs-1)     = decomp_x_end(nprocs-1) - decomp_x_start(nprocs-1) + 1
    end if

    ! This rank's x-slab intersected with the AABB+buffer box: the only range
    ! where near-surface distance computation actually needs to run. Empty
    ! (cx_lo > cx_hi) when this rank's slab lies entirely outside the AABB.
    cx_lo = max(decomp_x_start(myid), sx)
    cx_hi = min(decomp_x_end(myid),   ex)

    if (debug .and. myid == 0) then
        do ii = 0, nprocs-1
            print *, "PID:", ii, "xstart:", decomp_x_start(ii), "xend:", decomp_x_end(ii)
        end do
    end if

    call raycast_set_mesh(vertices, faces, face_solid_id, nfaces, nsolids, mesh_flip, mesh_active, mesh_open)

    if (myid == 0) then
        call cpu_time(time1)
        print *, "-- Pre-processing done in", time1-startTime, "s"
        call estimated_memoryusage(nprocs, nfaces, nvertices, nx, ny, nz)
        print *, "*** Calculating SDF | cell-centres ***"
    end if
    call MPI_BARRIER(MPI_COMM_WORLD, ierror)

    left_rank  = myid - 1
    right_rank = myid + 1
    if (myid == 0)        left_rank  = -1      ! Neumann BC at left boundary
    if (myid == nprocs-1) right_rank = nprocs  ! Neumann BC at right boundary

    ! Build spatial hash once (shared geometry, reused for all stagger locations).
    ! Needed unconditionally: raycast_classify_domain (the primary sign/objid
    ! source) relies on it too, not just the fast-sweep narrow-band path.
    call build_hash(vertices, faces, nfaces, bbox_min, bbox_max, dx, dy, minval(dz))

    ! Per-rank slab of the field (and of the solid-object ID); the assembled global array exists on rank 0
    !   only.
    allocate(sdf(decomp_x_start(myid):decomp_x_end(myid), ny, nz))
    allocate(objid(decomp_x_start(myid):decomp_x_end(myid), ny, nz))
    if (myid == 0) then
        allocate(flood_fill_arr(nx, ny, nz))
        allocate(objid_fill_arr(nx, ny, nz))
    end if

    ! Cell-centre SDF is always written; face-staggered SDFs (u, v, w) only when compute_face_sdf = .true.
    ! (objid is not written for the face passes: the solver's thermal ghost list is cell-centre-only.)
    call compute_field(xp, yp, zp, 'data/sdfp.bin', 'cell-centres', .true.)
    if (compute_face_sdf) then
        call compute_field(xf, yp, zp, 'data/sdfu.bin', 'u-faces', .false.)
        call compute_field(xp, yf, zp, 'data/sdfv.bin', 'v-faces', .false.)
        call compute_field(xp, yp, zf, 'data/sdfw.bin', 'w-faces', .false.)
    end if

    if (myid == 0) then
        call cpu_time(endTime)
        totalTime = endTime - startTime
        print *, "*** SDF complete in", totalTime, "s ***"
    end if

    call free_hash()
    call MPI_FINALIZE(ierror)
contains

    ! Append the hole lids made by prepare_mesh to the geometry arrays (already wound consistently with their
    !   neighbours,
    ! so they carry no flip; per-face geometric normals keep the narrow-band normal bookkeeping valid).
    subroutine append_lids()
        real(dp), allocatable :: new_normals(:,:)
        integer,  allocatable :: tmp_i(:,:)
        logical,  allocatable :: tmp_l(:)
        real(dp) :: e1(3), e2(3), nn(3), mag
        integer :: i, nf0, nn0

        nf0 = nfaces; nn0 = nnormals
        allocate(new_normals(3, ncap))
        do i = 1, ncap
            e1 = vertices(:,cap_faces(2,i)) - vertices(:,cap_faces(1,i))
            e2 = vertices(:,cap_faces(3,i)) - vertices(:,cap_faces(1,i))
            nn = [e1(2)*e2(3)-e1(3)*e2(2), e1(3)*e2(1)-e1(1)*e2(3), e1(1)*e2(2)-e1(2)*e2(1)]
            mag = sqrt(dot_product(nn, nn))
            if (mag > 0.0_dp) nn = nn / mag
            new_normals(:,i) = nn
        end do
        normals = reshape([normals(:,1:nn0), new_normals], [3, nn0 + ncap])
        allocate(tmp_i(3, ncap))
        tmp_i = faces(:,1:nf0)
        faces = reshape([faces(:,1:nf0), cap_faces], [3, nf0 + ncap])
        do i = 1, ncap
            tmp_i(:,i) = nn0 + i
        end do
        face_normals = reshape([face_normals(:,1:nf0), tmp_i], [3, nf0 + ncap])
        face_solid_id = [face_solid_id(1:nf0), cap_solid]
        allocate(tmp_l(nf0 + ncap))
        tmp_l(1:nf0) = mesh_flip; tmp_l(nf0+1:) = .false.
        call move_alloc(tmp_l, mesh_flip)
        allocate(tmp_l(nf0 + ncap))
        tmp_l(1:nf0) = mesh_active; tmp_l(nf0+1:) = .true.
        call move_alloc(tmp_l, mesh_active)
        nfaces = nf0 + ncap; nnormals = nn0 + ncap
        deallocate(new_normals, tmp_i)
    end subroutine append_lids

    ! One staggering of the field: unsigned near-surface distance per rank slab -> sign from the solid/fluid
    ! classification (each rank classifies its own slab) -> optional fast-sweep fill -> gather on rank 0 and
    !   write.
    subroutine compute_field(xs, ys, zs, fname, label, write_objid)
        implicit none
        real(dp), intent(in) :: xs(:), ys(:), zs(:)
        character(len=*), intent(in) :: fname, label
        logical, intent(in) :: write_objid

        sdf = scalarvalue
        objid = 0.0_dp
        if (cx_lo <= cx_hi) then
            if (use_fast_sweep) then
                call compute_narrowband_sdf(myid, cx_lo, cx_hi, &
                                            sy, ey, sz, ez, xs, ys, zs, nfaces, faces, face_normals, &
                                            vertices, normals, face_solid_id, narrow_band_width, &
                                            sdf(cx_lo:cx_hi,:,:), objid(cx_lo:cx_hi,:,:))
            else
                call compute_scalar_distance_face(myid, cx_lo, cx_hi, &
                                                  sy, ey, sz, ez, xs, ys, zs, nfaces, faces, face_normals, &
                                                  vertices, normals, face_solid_id, buffer_points, &
                                                  sdf(cx_lo:cx_hi,:,:), objid(cx_lo:cx_hi,:,:))
            end if
        end if
        call MPI_BARRIER(MPI_COMM_WORLD, ierror)

        ! Collective: ranks whose slab misses the box pass an empty range but still take part in the
        !   reductions.
        call raycast_classify_domain(vertices, faces, face_solid_id, nfaces, nsolids, &
                                     xs, ys, zs, sx, ex, sy, ey, sz, ez, cx_lo, cx_hi, scalarvalue, &
                                     sdf(cx_lo:cx_hi,:,:), objid(cx_lo:cx_hi,:,:), myid, nprocs, MPI_COMM_WORLD)

        if (use_fast_sweep) then
            ! objid is not touched by the fast-sweep propagation — it only needs to be correct in the
            !   narrow-band
            ! interface cells where ghost points actually live.
            call fast_sweep_3d(sdf, decomp_x_start(myid), decomp_x_end(myid), 1, ny, 1, nz, &
                               xs, ys, zs, myid, left_rank, right_rank)
        end if
        call gather_array(myid, nprocs, sdf, nx, ny, nz, decomp_x_start, decomp_x_end, decomp_size, &
                          1, ny, 1, nz, scalarvalue, flood_fill_arr, ierror)
        if (write_objid) &
            call gather_array(myid, nprocs, objid, nx, ny, nz, decomp_x_start, decomp_x_end, decomp_size, &
                              1, ny, 1, nz, 0.0_dp, objid_fill_arr, ierror)
        if (myid == 0) then
            call cpu_time(time1)
            print *, "*** Writing | ", label, " ***"
            call write_sdf_padded(fname, flood_fill_arr)
            if (write_objid) call write_sdf_padded('data/sdfp_objid.bin', objid_fill_arr)
            call cpu_time(time2)
            print *, "-- Write done in", time2-time1, "s |", label
            if (label == 'cell-centres' .and. compute_face_sdf) print *, "*** Calculating SDF | face-staggered ***"
        end if
        call MPI_BARRIER(MPI_COMM_WORLD, ierror)
    end subroutine compute_field

end program generatesdf
