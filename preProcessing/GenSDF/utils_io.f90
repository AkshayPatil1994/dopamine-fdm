module utils_io
    !
    ! utils_io — data, grid I/O, geometry I/O, and output routines for GenSDF v2.
    !
    implicit none

    !  precision 
    integer, parameter, public :: dp = selected_real_kind(15,307), &
                                  sp = selected_real_kind(6 , 37)

    !  geometry 
    character(len=512) :: inputfilename

    !  Cartesian grid 
    integer, protected :: nx, ny, nz
    real(dp), dimension(3), protected :: r0
    integer, dimension(3) :: ng
    logical, protected :: non_uniform_grid
    real(dp), allocatable, dimension(:) :: xp, yp, zp, xf, yf, zf
    real(dp) :: lx, ly, lz, dx, dy, dx_inverse, dy_inverse
    real(dp), allocatable, dimension(:) :: dz, dz_inverse
    ! per-cell width along GenSDF's internal y (= solver z when vertical_axis==2); uniform unless
    !   data/grid_z.out is present
    real(dp), allocatable, dimension(:) :: dyv
    logical :: y_stretched = .false.
    ! Per-solid flag set from a trailing `ground` keyword in the .list manifest: the solid is an open surface
    ! (terrain/bathymetry) and everything below it along the vertical axis is solid. Sized by
    !   load_multi_object_geometry.
    logical, allocatable :: solid_is_ground(:)

    !  SDF / algorithm control 
    real(dp), protected :: scalarvalue
    integer,  protected :: buffer_points
    logical,  protected :: use_fast_sweep
    integer,  protected :: narrow_band_width
    !  2 = solver wall-normal is y;  3 = solver wall-normal is z (default, z is wall-normal)
    integer,  protected :: vertical_axis = 3
    ! Compute face SDFs (sdfu, sdfv, sdfw) in addition to cell-centre (sdfp).
    ! Default .false. — only sdfp is written.
    logical,  protected :: compute_face_sdf = .false.
    ! Solid/fluid classification of the grid points (optional last line of parameters.in):
    !   0 = 3-axis ray-parity majority vote (needs a watertight shell)
    !   1 = default: orientation-signed ray crossings with per-line integrity checks, undecided
    !       points
    !       resolved by the generalized winding number (tolerates holes, cracks,
    !         flipped/duplicated/overlapping triangles)
    !   2 = winding number at every point (slow reference; also prints how often the ray tier disagrees)
    integer,  protected :: sign_method = 1
    ! Close simple, nearly planar boundary loops of non-ground solids with a triangulated lid before
    !   classifying (1 = yes).
    ! Optional second trailing line of parameters.in.
    integer,  protected :: cap_holes = 1

    !  MPI 
    integer :: nprocs

    !  Auxiliary
    integer, protected :: pbarwidth

    !$acc declare create(dx, dy, scalarvalue)
    ! dz is referenced inside godunov_update_cell (!$acc routine seq), which requires
    ! module allocatables it touches to be declared, not just enter-data'd.
    !$acc declare create(dz)
    !$acc declare create(dyv)

contains

    ! Read input file 
    subroutine read_inputfile(procid)
        implicit none
        integer, intent(in) :: procid
        integer :: iunit, ierr
        character(len=512) :: dummyline

        iunit = 10
        open(newunit=iunit, file='parameters.in', status='old', action='read', iostat=ierr)
        if (ierr == 0) then
            call skip_and_read(iunit, dummyline, ierr); read(dummyline,*) inputfilename
            call skip_and_read(iunit, dummyline, ierr); read(dummyline,*) scalarvalue, buffer_points, pbarwidth
            call skip_and_read(iunit, dummyline, ierr); read(dummyline,*) nx, ny, nz
            call skip_and_read(iunit, dummyline, ierr); read(dummyline,*) r0(1), r0(2), r0(3)
            call skip_and_read(iunit, dummyline, ierr); read(dummyline,*) non_uniform_grid
            call skip_and_read(iunit, dummyline, ierr); read(dummyline,*) use_fast_sweep
            call skip_and_read(iunit, dummyline, ierr); read(dummyline,*) narrow_band_width
            call skip_and_read(iunit, dummyline, ierr); read(dummyline,*) vertical_axis
            call skip_and_read(iunit, dummyline, ierr); read(dummyline,*) compute_face_sdf
            ! Optional trailing line: sign_method (absent => default 1)
            call skip_and_read(iunit, dummyline, ierr)
            if (ierr == 0) read(dummyline,*,iostat=ierr) sign_method
            if (sign_method < 0 .or. sign_method > 2) then
                if (procid == 0) print *, "ERROR: sign_method must be 0, 1 or 2, got:", sign_method
                error stop
            end if
            ! Optional: cap_holes (1 = close planar holes with lids, 0 = leave the mesh as given)
            call skip_and_read(iunit, dummyline, ierr)
            if (ierr == 0) read(dummyline,*,iostat=ierr) cap_holes
            if (vertical_axis /= 2 .and. vertical_axis /= 3) then
                if (procid == 0) &
                    print *, "ERROR: vertical_axis must be 2 or 3, got:", vertical_axis
                error stop
            end if
            ! Permute so GenSDF's internal dir 3 (z) always = wall-normal direction
            if (vertical_axis == 2) then
                ng = [nx, nz, ny]
                nx = ng(1); ny = ng(2); nz = ng(3)
            end if
            ng(1) = nx
            ng(2) = ny
            ng(3) = nz
            if (procid == 0) print *, "*** Input file successfully read ***"
        else
            if (procid == 0) print *, "ERROR: parameters.in encountered a problem!"
            error stop
        end if
        close(iunit)

    contains
        ! Skip blank lines and lines whose first non-blank character is '#', then
        ! return the first data line in `line`.
        subroutine skip_and_read(u, line, ios)
            integer,          intent(in)  :: u
            character(len=*), intent(out) :: line
            integer,          intent(out) :: ios
            character(len=1) :: fc
            integer :: pos
            do
                read(u,'(A)',iostat=ios) line
                if (ios /= 0) return
                pos = verify(line, ' ')           ! index of first non-space
                if (pos == 0) cycle               ! blank line → skip
                fc = line(pos:pos)
                if (fc /= '#') return             ! found a data line
                ! else: comment line → skip
            end do
        end subroutine skip_and_read
    end subroutine read_inputfile

    ! CaNS grid still default as it is easy and common 
    ! Read CaNS-style grid files: geometry.out (grid dimensions and lengths) and optionally grid.out
    !   (non-uniform z).
    subroutine read_cans_grid(procid, loc, iprecision, npoints, origin, non_uni_grid, &
                               xin_p, yin_p, zin_p, xin_f, yin_f, zin_f, dzin, &
                               inputlx, inputly, inputlz)
        implicit none
        integer,  intent(in)  :: procid, iprecision
        character(len=*), intent(in) :: loc
        real(dp), dimension(3), intent(in) :: origin
        logical,  intent(in)  :: non_uni_grid
        integer,  dimension(3), intent(out) :: npoints
        real(dp), intent(out), allocatable, dimension(:) :: xin_p, yin_p, zin_p
        real(dp), intent(out), allocatable, dimension(:) :: xin_f, yin_f, zin_f
        real(dp), intent(out), allocatable, dimension(:) :: dzin
        real(dp), intent(out) :: inputlx, inputly, inputlz
        logical :: fexists
        integer :: iter, ii, jj, unit
        real(dp) :: dl(3), l(3)
        real(sp), allocatable :: grid_z4(:,:)
        real(dp), allocatable :: grid_z8(:,:)
        character(len=512) :: geofile, grdfile, linebuf
        integer :: filerr

        inquire(file=trim(loc), exist=fexists)
        if (.not. fexists) then
            print *, "Input directory does not exist: ", trim(loc)
            stop
        end if

        geofile = trim(loc) // "geometry.out"
        open(newunit=unit, file=geofile, status='old')
        ! Skip comment lines starting with '#'
        do
            read(unit,'(A)',iostat=filerr) linebuf
            if (filerr /= 0) exit
            if (len_trim(adjustl(linebuf)) == 0) cycle
            if (linebuf(verify(linebuf,' '):verify(linebuf,' ')) /= '#') then
                read(linebuf, *) npoints; exit
            end if
        end do
        do
            read(unit,'(A)',iostat=filerr) linebuf
            if (filerr /= 0) exit
            if (len_trim(adjustl(linebuf)) == 0) cycle
            if (linebuf(verify(linebuf,' '):verify(linebuf,' ')) /= '#') then
                read(linebuf, *) l; exit
            end if
        end do
        close(unit)
        ! Permute so GenSDF's internal dir 3 (z) = wall-normal direction
        if (vertical_axis == 2) then
            npoints = [npoints(1), npoints(3), npoints(2)]
            l       = [l(1),       l(3),       l(2)      ]
        end if
        inputlx = l(1); inputly = l(2); inputlz = l(3)
        dl = l / real(npoints, dp)

        allocate(xin_p(npoints(1)), yin_p(npoints(2)), zin_p(npoints(3)))
        allocate(xin_f(npoints(1)), yin_f(npoints(2)), zin_f(npoints(3)))
        allocate(dzin(npoints(3)))
        allocate(dyv(npoints(2)))
        dyv = dl(2)

        do iter = 1, npoints(1)
            xin_p(iter) = origin(1) + (iter - 0.5_dp) * dl(1)
            xin_f(iter) = xin_p(iter) + dl(1) / 2.0_dp
        end do
        do iter = 1, npoints(2)
            yin_p(iter) = origin(2) + (iter - 0.5_dp) * dl(2)
            yin_f(iter) = yin_p(iter) + dl(2) / 2.0_dp
        end do
        do iter = 1, npoints(3)
            zin_p(iter) = origin(3) + (iter - 0.5_dp) * dl(3)
            zin_f(iter) = zin_p(iter) + dl(3) / 2.0_dp
        end do

        if (non_uni_grid) then
            grdfile = trim(loc) // "grid.out"
            if (iprecision == 4) then
                open(newunit=unit, file=grdfile, status='old', action='read', iostat=filerr)
                if (filerr /= 0) error stop "Error opening grid.out"
                allocate(grid_z4(npoints(3),5))
                do ii = 1, npoints(3)
                    read(unit,*,iostat=filerr) (grid_z4(ii,jj), jj=1,5)
                    if (filerr /= 0) error stop "Error reading grid.out"
                end do
                close(unit)
                zin_p = origin(3) + grid_z4(:,2)
                zin_f = origin(3) + grid_z4(:,3)
            else
                open(newunit=unit, file=grdfile, status='old', action='read', iostat=filerr)
                if (filerr /= 0) error stop "Error opening grid.out"
                allocate(grid_z8(npoints(3),5))
                do ii = 1, npoints(3)
                    read(unit,*,iostat=filerr) (grid_z8(ii,jj), jj=1,5)
                    if (filerr /= 0) error stop "Error reading grid.out"
                end do
                close(unit)
                zin_p = origin(3) + grid_z8(:,2)
                zin_f = origin(3) + grid_z8(:,3)
            end if
        end if

        ! Optional stretched spanwise grid (solver z, GenSDF's internal y when vertical_axis==2): same
        !   5-column
        ! format as grid.out, written by the solver's genGridandIC when alpha_grid_z>0; absent => uniform y.
        if (vertical_axis == 2) then
            grdfile = trim(loc) // "grid_z.out"
            inquire(file=trim(grdfile), exist=fexists)
            if (fexists) then
                if (allocated(grid_z8)) deallocate(grid_z8)
                allocate(grid_z8(npoints(2),5))
                open(newunit=unit, file=grdfile, status='old', action='read', iostat=filerr)
                if (filerr /= 0) error stop "Error opening grid_z.out"
                do ii = 1, npoints(2)
                    read(unit,*,iostat=filerr) (grid_z8(ii,jj), jj=1,5)
                    if (filerr /= 0) error stop "Error reading grid_z.out"
                end do
                close(unit)
                yin_p = origin(2) + grid_z8(:,2)
                yin_f = origin(2) + grid_z8(:,3)
                dyv   = 2.0_dp * (yin_f - yin_p)     ! cell width (centres are face midpoints)
                y_stretched = .true.
                deallocate(grid_z8)
                if (procid == 0) print *, "*** Read stretched spanwise grid from grid_z.out ***"
            end if
        end if

        if (procid == 0) print *, "*** Successfully read the grid ***"
    end subroutine read_cans_grid

    ! --
    subroutine setup_grid_spacing(procid, xin, yin, zin, nz_in, &
                                   dx_out, dy_out, dz_out, dx_inv_out, dy_inv_out, dz_inv_out)
        implicit none
        integer,  intent(in)  :: procid, nz_in
        real(dp), intent(in),  dimension(:) :: xin, yin, zin
        real(dp), intent(out) :: dx_out, dy_out, dx_inv_out, dy_inv_out
        real(dp), intent(out), dimension(:) :: dz_out, dz_inv_out
        integer :: ii

        dx_out = xin(2) - xin(1)
        dy_out = yin(2) - yin(1)
        dx_inv_out = 1.0_dp / dx_out
        dy_inv_out = 1.0_dp / dy_out
        ! z spacing: face(ii+1) - face(ii) gives exact width of cell ii
        ! Loop only to nz_in-1 to avoid reading zin(nz_in+1) past the array end;
        ! the last cell reuses the second-to-last spacing.
        do ii = 1, nz_in - 1
            dz_out(ii)     = zin(ii+1) - zin(ii)
            dz_inv_out(ii) = 1.0_dp    / dz_out(ii)
        end do
        dz_out(nz_in)     = dz_out(nz_in - 1)
        dz_inv_out(nz_in) = dz_inv_out(nz_in - 1)

        if (procid == 0) print *, "*** Successfully set up grid spacing ***"
    end subroutine setup_grid_spacing

    ! Read OBJ file (ASCII). Polygons are fan-triangulated (a quad exported as 'f 1 2 3 4' must not silently
    !   lose half its
    ! area: that opens a hole in the shell). Face tokens may be v, v/vt, v//vn or v/vt/vn, with negative
    !   (relative) indices.
    ! Per-face geometric normals are synthesized when the file carries no (or incomplete) vertex normals.
    subroutine read_obj(procid, filename, vertices, normals, faces, face_normals, &
                        num_vertices, num_normals, num_faces)
        implicit none
        integer,  intent(in)  :: procid
        character(len=*), intent(in) :: filename
        real(dp), allocatable, intent(out) :: vertices(:,:), normals(:,:)
        integer,  allocatable, intent(out) :: faces(:,:), face_normals(:,:)
        integer,  intent(out) :: num_vertices, num_normals, num_faces
        integer, parameter :: max_tok = 256
        integer :: unit, ios, v_count, vn_count, f_count, ntok, itok
        integer :: idx_v(max_tok), idx_n(max_tok)
        character(len=2048) :: line
        character(len=64)   :: tok(max_tok)
        real(dp) :: p(3)
        logical  :: have_all_normals

        open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, 'Error opening file: ', trim(filename); error stop
        end if

        v_count = 0; vn_count = 0; f_count = 0
        do
            read(unit,'(A)',iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (line(1:2) == 'v ') then
                v_count = v_count + 1
            else if (line(1:3) == 'vn ') then
                vn_count = vn_count + 1
            else if (line(1:2) == 'f ') then
                call split_tokens(line(3:), tok, ntok, max_tok)
                f_count = f_count + max(0, ntok - 2)
            end if
        end do

        allocate(vertices(3,v_count), normals(3,max(vn_count,1)))
        allocate(faces(3,f_count), face_normals(3,f_count))
        num_vertices = v_count; num_faces = f_count
        have_all_normals = (vn_count > 0)

        rewind(unit)
        v_count = 0; vn_count = 0; f_count = 0
        do
            read(unit,'(A)',iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (line(1:2) == 'v ') then
                read(line(3:),*) p
                v_count = v_count + 1
                vertices(:, v_count) = p
            else if (line(1:3) == 'vn ') then
                read(line(4:),*) p
                vn_count = vn_count + 1
                normals(:, vn_count) = p
            else if (line(1:2) == 'f ') then
                call split_tokens(line(3:), tok, ntok, max_tok)
                if (ntok < 3) cycle
                do itok = 1, ntok
                    call parse_face_token(tok(itok), v_count, vn_count, idx_v(itok), idx_n(itok))
                    if (idx_n(itok) == 0) have_all_normals = .false.
                end do
                do itok = 2, ntok - 1
                    f_count = f_count + 1
                    faces(:, f_count)        = [idx_v(1), idx_v(itok), idx_v(itok+1)]
                    face_normals(:, f_count) = [idx_n(1), idx_n(itok), idx_n(itok+1)]
                end do
            end if
        end do
        close(unit)

        if (num_faces > 0) then
            if (minval(faces) < 1 .or. maxval(faces) > num_vertices) then
                print *, 'FATAL: OBJ face references a vertex that does not exist: ', trim(filename); error stop
            end if
        end if

        num_normals = vn_count
        if (have_all_normals .and. num_faces > 0) then
            if (minval(face_normals) < 1 .or. maxval(face_normals) > num_normals) have_all_normals = .false.
        end if
        if (.not. have_all_normals) then
            ! Winding-derived per-face normals (the SDF sign is decided by the ray/winding classifier, not by
            !   these)
            deallocate(normals); allocate(normals(3, max(num_faces,1)))
            call synthesize_face_normals(vertices, faces, num_faces, normals, face_normals)
            num_normals = num_faces
            if (procid == 0) print *, 'Note: OBJ has no usable vertex normals; using per-face geometric normals.'
        end if

        call permute_geometry_axes(vertices, normals, faces, face_normals, num_vertices, num_normals, num_faces)

        if (procid == 0) then
            print *, 'Successfully read OBJ: ', trim(filename)
            print *, 'Vertices: ', num_vertices, '  Normals: ', num_normals, '  Triangles: ', num_faces
        end if
    end subroutine read_obj

    ! Split a text line into blank-separated tokens (tabs count as blanks); tokens past max_tok are dropped.
    subroutine split_tokens(str, tok, ntok, max_tok)
        implicit none
        character(len=*), intent(in)  :: str
        integer,          intent(in)  :: max_tok
        character(len=*), intent(out) :: tok(max_tok)
        integer,          intent(out) :: ntok
        integer :: i, n, start
        logical :: in_tok
        ntok = 0; in_tok = .false.; start = 1
        n = len_trim(str)
        do i = 1, n + 1
            if (i <= n) then
                if (str(i:i) == ' ' .or. str(i:i) == char(9) .or. str(i:i) == char(13)) then
                    if (in_tok) then
                        if (ntok < max_tok) then; ntok = ntok + 1; tok(ntok) = str(start:i-1); end if
                        in_tok = .false.
                    end if
                else if (.not. in_tok) then
                    in_tok = .true.; start = i
                end if
            else if (in_tok) then
                if (ntok < max_tok) then; ntok = ntok + 1; tok(ntok) = str(start:n); end if
            end if
        end do
    end subroutine split_tokens

    ! One OBJ face token: v, v/vt, v//vn or v/vt/vn. Negative indices count back from the last element read so
    !   far.
    ! vn_idx = 0 when the token carries no normal.
    subroutine parse_face_token(tok, v_count, vn_count, v_idx, vn_idx)
        implicit none
        character(len=*), intent(in)  :: tok
        integer,          intent(in)  :: v_count, vn_count
        integer,          intent(out) :: v_idx, vn_idx
        integer :: s1, s2, n
        n  = len_trim(tok)
        s1 = index(tok(1:n), '/')
        v_idx = 0; vn_idx = 0
        if (s1 == 0) then
            read(tok(1:n), *) v_idx
        else
            read(tok(1:s1-1), *) v_idx
            s2 = index(tok(s1+1:n), '/')
            if (s2 > 0) then
                s2 = s1 + s2
                if (s2 < n) read(tok(s2+1:n), *) vn_idx
            end if
        end if
        if (v_idx  < 0) v_idx  = v_count  + 1 + v_idx
        if (vn_idx < 0) vn_idx = vn_count + 1 + vn_idx
    end subroutine parse_face_token

    ! One unit normal per face from its winding; face_normals(:,f) all point at that face's entry
    subroutine synthesize_face_normals(vertices, faces, nf, normals, face_normals)
        implicit none
        real(dp), intent(in)  :: vertices(:,:)
        integer,  intent(in)  :: faces(:,:), nf
        real(dp), intent(out) :: normals(:,:)
        integer,  intent(out) :: face_normals(:,:)
        integer  :: f
        real(dp) :: e1(3), e2(3), n(3), m
        do f = 1, nf
            e1 = vertices(:,faces(2,f)) - vertices(:,faces(1,f))
            e2 = vertices(:,faces(3,f)) - vertices(:,faces(1,f))
            n  = [e1(2)*e2(3)-e1(3)*e2(2), e1(3)*e2(1)-e1(1)*e2(3), e1(1)*e2(2)-e1(2)*e2(1)]
            m  = sqrt(dot_product(n,n))
            if (m > 0.0_dp) n = n / m
            normals(:,f) = n
            face_normals(:,f) = f
        end do
    end subroutine synthesize_face_normals

    ! GenSDF's internal direction 3 is the wall-normal axis; with vertical_axis == 2 swap geometry y<->z.
    subroutine permute_geometry_axes(vertices, normals, faces, face_normals, nv, nn, nf)
        implicit none
        real(dp), intent(inout) :: vertices(:,:), normals(:,:)
        integer,  intent(inout) :: faces(:,:), face_normals(:,:)
        integer,  intent(in)    :: nv, nn, nf
        real(dp), allocatable :: tmp_row(:)
        integer,  allocatable :: tmp_irow(:)
        if (vertical_axis /= 2) return
        allocate(tmp_row(nv))
        tmp_row = vertices(2,:nv); vertices(2,:nv) = vertices(3,:nv); vertices(3,:nv) = tmp_row
        deallocate(tmp_row)
        allocate(tmp_row(nn))
        tmp_row = normals(2,:nn); normals(2,:nn) = normals(3,:nn); normals(3,:nn) = tmp_row
        deallocate(tmp_row)
        ! y<->z is a reflection (det=-1): reverse winding so cross-product normals stay outward
        allocate(tmp_irow(nf))
        tmp_irow = faces(2,:nf); faces(2,:nf) = faces(3,:nf); faces(3,:nf) = tmp_irow
        tmp_irow = face_normals(2,:nf); face_normals(2,:nf) = face_normals(3,:nf); face_normals(3,:nf) = tmp_irow
        deallocate(tmp_irow)
    end subroutine permute_geometry_axes

    ! Read STL file: binary (size-checked) or ASCII ('solid ... facet normal ... vertex ...'), detected from
    !   the file itself.
    subroutine read_stl_binary(procid, filename, vertices, normals, faces, face_normals, &
                               num_vertices, num_normals, num_faces)
        !
        ! Binary layout (little-endian): 80-byte header, uint32 triangle count, then per triangle 12 float32
        ! (normal, 3 vertices) and 2 attribute bytes. Output arrays:
        !   vertices(3, 3*nf) — 3 vertices per triangle (no sharing), normals(3, nf) — per-face, faces packed
        !     sequentially.
        !
        use iso_fortran_env, only : int32, int64
        implicit none
        integer,  intent(in)  :: procid
        character(len=*), intent(in) :: filename
        real(dp), allocatable, intent(out) :: vertices(:,:), normals(:,:)
        integer,  allocatable, intent(out) :: faces(:,:), face_normals(:,:)
        integer,  intent(out) :: num_vertices, num_normals, num_faces
        integer            :: unit, ios, fi
        integer(int32)     :: nf_i32
        integer(int64)     :: fsize
        character(len=80)  :: header
        integer(2)         :: attr
        real(4)            :: fn(3), v0(3), v1(3), v2(3)  ! float32 reads
        logical            :: is_binary

        open(newunit=unit, file=trim(filename), access='stream', form='unformatted', &
             status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, "Error opening STL file: ", trim(filename); error stop
        end if
        inquire(unit=unit, size=fsize)

        is_binary = .false.
        if (fsize >= 84) then
            read(unit, pos=1) header
            read(unit) nf_i32
            is_binary = (fsize == 84_int64 + 50_int64 * int(max(nf_i32,0_int32), int64))
        end if
        if (.not. is_binary) then
            close(unit)
            if (fsize >= 5) then
                call read_stl_ascii(procid, filename, vertices, normals, faces, face_normals, &
                                    num_vertices, num_normals, num_faces)
                return
            end if
            print *, "FATAL: not a valid STL file: ", trim(filename); error stop
        end if

        num_faces    = int(nf_i32)
        num_vertices = 3 * num_faces
        num_normals  = num_faces

        allocate(vertices(3, num_vertices))
        allocate(normals (3, num_normals))
        allocate(faces      (3, num_faces))
        allocate(face_normals(3, num_faces))

        do fi = 1, num_faces
            read(unit) fn(1), fn(2), fn(3)    ! face normal (float32)
            read(unit) v0(1), v0(2), v0(3)    ! vertex 0
            read(unit) v1(1), v1(2), v1(3)    ! vertex 1
            read(unit) v2(1), v2(2), v2(3)    ! vertex 2
            read(unit) attr                   ! attribute bytes (ignored)

            normals(:,fi) = real(fn,dp)
            vertices(:, 3*fi-2) = real(v0,dp)
            vertices(:, 3*fi-1) = real(v1,dp)
            vertices(:, 3*fi  ) = real(v2,dp)
            faces(:,fi)        = [3*fi-2, 3*fi-1, 3*fi]
            face_normals(:,fi) = [fi, fi, fi]
        end do
        close(unit)

        call permute_geometry_axes(vertices, normals, faces, face_normals, num_vertices, num_normals, num_faces)

        if (procid == 0) then
            print *, "Successfully read binary STL: ", trim(filename)
            print *, "Faces: ", num_faces
        end if
    end subroutine read_stl_binary

    subroutine read_stl_ascii(procid, filename, vertices, normals, faces, face_normals, &
                              num_vertices, num_normals, num_faces)
        implicit none
        integer,  intent(in)  :: procid
        character(len=*), intent(in) :: filename
        real(dp), allocatable, intent(out) :: vertices(:,:), normals(:,:)
        integer,  allocatable, intent(out) :: faces(:,:), face_normals(:,:)
        integer,  intent(out) :: num_vertices, num_normals, num_faces
        integer :: unit, ios, fi, nvert_in_facet, p
        character(len=512) :: line
        real(dp) :: q(3)

        open(newunit=unit, file=trim(filename), form='formatted', access='sequential', status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, "Error opening STL file: ", trim(filename); error stop
        end if
        num_faces = 0
        do
            read(unit,'(A)',iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (line(1:5) == 'facet' .or. line(1:5) == 'FACET') num_faces = num_faces + 1
        end do
        if (num_faces == 0) then
            print *, "FATAL: STL file has no facets (neither valid binary nor ASCII): ", trim(filename); error stop
        end if
        num_vertices = 3 * num_faces; num_normals = num_faces
        allocate(vertices(3,num_vertices), normals(3,num_normals), faces(3,num_faces), face_normals(3,num_faces))
        normals = 0.0_dp

        rewind(unit)
        fi = 0; nvert_in_facet = 0
        do
            read(unit,'(A)',iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (line(1:5) == 'facet' .or. line(1:5) == 'FACET') then
                fi = fi + 1; nvert_in_facet = 0
                p = index(line, 'normal'); if (p == 0) p = index(line, 'NORMAL')
                if (p > 0) then
                    read(line(p+6:),*,iostat=ios) q
                    if (ios == 0) normals(:,fi) = q
                end if
            else if (line(1:6) == 'vertex' .or. line(1:6) == 'VERTEX') then
                read(line(7:),*) q
                nvert_in_facet = nvert_in_facet + 1
                if (nvert_in_facet <= 3) vertices(:, 3*(fi-1)+nvert_in_facet) = q
            end if
        end do
        close(unit)
        do fi = 1, num_faces
            faces(:,fi)        = [3*fi-2, 3*fi-1, 3*fi]
            face_normals(:,fi) = [fi, fi, fi]
        end do

        call permute_geometry_axes(vertices, normals, faces, face_normals, num_vertices, num_normals, num_faces)
        if (procid == 0) then
            print *, "Successfully read ASCII STL: ", trim(filename)
            print *, "Faces: ", num_faces
        end if
    end subroutine read_stl_ascii

    ! Load geometry for one or more solids, tagging every face with its owning solid ID.
    subroutine load_multi_object_geometry(procid, filename, vertices, normals, faces, face_normals, &
                                          face_solid_id, num_vertices, num_normals, num_faces, num_solids)
        ! `filename` is a single .stl/.obj file (solid ID 1) or a `.list` manifest (one path per line, 1-based
        !   order is the solid ID).
        ! A manifest line ending in ` ground` marks that solid as an open ground surface (see
        !   solid_is_ground).
        implicit none
        integer,  intent(in)  :: procid
        character(len=*), intent(in) :: filename
        real(dp), allocatable, intent(out) :: vertices(:,:), normals(:,:)
        integer,  allocatable, intent(out) :: faces(:,:), face_normals(:,:)
        integer,  allocatable, intent(out) :: face_solid_id(:)
        integer,  intent(out) :: num_vertices, num_normals, num_faces, num_solids

        character(len=512), allocatable :: manifest(:)
        real(dp), allocatable :: v_i(:,:), n_i(:,:)
        integer,  allocatable :: f_i(:,:), fn_i(:,:)
        integer :: nv_i, nn_i, nf_i
        integer :: solid, unit, ios, nlines
        character(len=512) :: line, aline
        character(len=4)   :: ext
        integer :: ii, pos, vbase, nbase, fbase

        ii = len_trim(filename)
        ext = filename(max(1,ii-3):ii)

        if (ext == '.lst' .or. filename(max(1,ii-4):ii) == '.list') then
            ! Count manifest entries first
            open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
            if (ios /= 0) then
                print *, 'Error opening manifest file: ', trim(filename); stop
            end if
            nlines = 0
            do
                read(unit,'(A)',iostat=ios) line
                if (ios /= 0) exit
                aline = adjustl(line)
                pos = verify(aline, ' ')
                if (pos == 0 .or. aline(pos:pos) == '#') cycle
                nlines = nlines + 1
            end do
            rewind(unit)
            allocate(manifest(nlines))
            nlines = 0
            do
                read(unit,'(A)',iostat=ios) line
                if (ios /= 0) exit
                aline = adjustl(line)
                pos = verify(aline, ' ')
                if (pos == 0 .or. aline(pos:pos) == '#') cycle
                nlines = nlines + 1
                manifest(nlines) = trim(aline)
            end do
            close(unit)
        else
            allocate(manifest(1))
            manifest(1) = trim(filename)
        end if

        num_solids   = size(manifest)
        allocate(solid_is_ground(num_solids))
        solid_is_ground = .false.
        do solid = 1, num_solids
            ii = len_trim(manifest(solid))
            if (ii > 7) then
                if (manifest(solid)(ii-6:ii) == ' ground') then
                    solid_is_ground(solid) = .true.
                    manifest(solid) = manifest(solid)(1:ii-7)
                end if
            end if
        end do
        num_vertices = 0; num_normals = 0; num_faces = 0
        allocate(vertices(3,0), normals(3,0), faces(3,0), face_normals(3,0), face_solid_id(0))

        do solid = 1, num_solids
            ii = len_trim(manifest(solid))
            ext = manifest(solid)(max(1,ii-3):ii)
            if (ext == '.stl' .or. ext == '.STL') then
                call read_stl_binary(procid, trim(manifest(solid)), v_i, n_i, f_i, fn_i, nv_i, nn_i, nf_i)
            else
                call read_obj(procid, trim(manifest(solid)), v_i, n_i, f_i, fn_i, nv_i, nn_i, nf_i)
            end if

            vbase = num_vertices; nbase = num_normals; fbase = num_faces

            vertices     = reshape([vertices, v_i],           [3, vbase + nv_i])
            normals      = reshape([normals,  n_i],            [3, nbase + nn_i])
            faces        = reshape([faces,    f_i + vbase],    [3, fbase + nf_i])
            face_normals = reshape([face_normals, fn_i + nbase], [3, fbase + nf_i])
            face_solid_id = [face_solid_id, spread(solid, 1, nf_i)]

            num_vertices = vbase + nv_i
            num_normals  = nbase + nn_i
            num_faces    = fbase + nf_i

            deallocate(v_i, n_i, f_i, fn_i)
        end do

        if (procid == 0) then
            print *, 'Loaded ', num_solids, ' solid(s), total faces: ', num_faces
        end if
    end subroutine load_multi_object_geometry

    ! Get axis-aligned bounding box of the geometry
    subroutine getbbox(procid, vertices_in, num_vertices_in, bbox_min, bbox_max)
        implicit none
        integer,  intent(in)  :: procid, num_vertices_in
        real(dp), intent(in)  :: vertices_in(3, num_vertices_in)
        real(dp), intent(out) :: bbox_min(3), bbox_max(3)
        integer :: iter
        bbox_min = vertices_in(:,1); bbox_max = vertices_in(:,1)
        do iter = 2, num_vertices_in
            where (vertices_in(:,iter) < bbox_min) bbox_min = vertices_in(:,iter)
            where (vertices_in(:,iter) > bbox_max) bbox_max = vertices_in(:,iter)
        end do
        if (procid == 0) then
            print *, "Geometry bbox min:", bbox_min
            print *, "Geometry bbox max:", bbox_max
        end if
    end subroutine getbbox

    ! Warn about solids that are too small, in any direction, relative to the local grid
    ! spacing to be reliably captured: a solid narrower than ~2 cells can fail to enclose
    ! any grid node at all, at which point no per-point sampling strategy (narrow-band or
    ! brute-force) can detect it, and fill_internal's topological classification has no
    ! shell of near-surface cells to close around it.
    subroutine check_object_resolution(procid, vertices, faces, face_solid_id, nfaces, nsolids)
        implicit none
        integer,  intent(in) :: procid, nfaces, nsolids
        real(dp), intent(in) :: vertices(:,:)
        integer,  intent(in) :: faces(:,:)
        integer,  intent(in) :: face_solid_id(:)
        real(dp), allocatable :: smin(:,:), smax(:,:)
        logical,  allocatable :: seen(:)
        integer  :: fi, sid, vi, vidx, d, kk
        real(dp) :: extent(3), local_dz, local_dy, min_cells
        real(dp), parameter :: MIN_CELLS_OK = 2.0_dp

        if (procid /= 0) return

        allocate(smin(3,nsolids), smax(3,nsolids), seen(nsolids))
        seen = .false.

        do fi = 1, nfaces
            sid = face_solid_id(fi)
            do vi = 1, 3
                vidx = faces(vi,fi)
                if (.not. seen(sid)) then
                    smin(:,sid) = vertices(:,vidx)
                    smax(:,sid) = vertices(:,vidx)
                    seen(sid)   = .true.
                else
                    do d = 1, 3
                        if (vertices(d,vidx) < smin(d,sid)) smin(d,sid) = vertices(d,vidx)
                        if (vertices(d,vidx) > smax(d,sid)) smax(d,sid) = vertices(d,vidx)
                    end do
                end if
            end do
        end do

        do sid = 1, nsolids
            if (.not. seen(sid)) cycle
            extent = smax(:,sid) - smin(:,sid)

            ! Coarsest local z spacing spanned by this solid (dz may be non-uniform);
            ! conservative on purpose so a stretched region isn't reported as resolved
            ! just because the grid is fine somewhere else in the domain.
            local_dz = dz(1)
            do kk = 1, size(zp)
                if (zp(kk) >= smin(3,sid) - dz(kk) .and. zp(kk) <= smax(3,sid) + dz(kk)) &
                    local_dz = max(local_dz, dz(kk))
            end do

            local_dy = dyv(1)
            do kk = 1, size(yp)
                if (yp(kk) >= smin(2,sid) - dyv(kk) .and. yp(kk) <= smax(2,sid) + dyv(kk)) &
                    local_dy = max(local_dy, dyv(kk))
            end do

            min_cells = min(extent(1)/dx, extent(2)/local_dy, extent(3)/local_dz)
            if (min_cells < MIN_CELLS_OK) then
                print '(A,I0,A,F6.2,A)', 'WARNING: solid ', sid, &
                    ' spans only ', min_cells, &
                    ' grid cells in its narrowest direction — below reliable SDF ' // &
                    'resolution (need >= 2); it may be dropped or misrepresented.'
            end if
        end do

        deallocate(smin, smax, seen)
    end subroutine check_object_resolution

    ! Tag grid points inside bounding box (with buffer) for focused SDF computation.
    ! The buffer is `sdfresolution` cells. With a stretched grid pass the per-cell widths (dy_cell/dz_cell) so
    !   the buffer
    ! is that many *local* cells rather than that many cells of the smallest spacing (which can be a small
    !   fraction of a
    ! cell where the grid is coarse, leaving the near-surface band outside the tagged box).
    subroutine tagminmax(procid, xin, yin, zin, bbox_min, bbox_max, &
                          nx_in, ny_in, nz_in, dx_in, dy_in, dz_in, sdfresolution, &
                          sx, ex, sy, ey, sz, ez, dy_cell, dz_cell)
        implicit none
        integer,  intent(in)  :: procid, nx_in, ny_in, nz_in, sdfresolution
        real(dp), dimension(:), intent(in) :: xin, yin, zin
        real(dp), intent(in)  :: dx_in, dy_in, dz_in
        real(dp), dimension(:), intent(in) :: bbox_min, bbox_max
        integer,  intent(out) :: sx, ex, sy, ey, sz, ez
        real(dp), intent(in), optional :: dy_cell(:), dz_cell(:)
        integer  :: iteri, iterj, iterk
        real(dp) :: buf, wy, wz

        sx = 1; ex = nx_in; sy = 1; ey = ny_in; sz = 1; ez = nz_in
        buf = real(sdfresolution, dp)

        do iteri = 1, nx_in
            if (xin(iteri) <= bbox_min(1) - buf*dx_in) then; sx = iteri
            else if (xin(iteri) <= bbox_max(1) + buf*dx_in) then; ex = iteri
            else; exit; end if
        end do
        do iterj = 1, ny_in
            wy = dy_in; if (present(dy_cell)) wy = dy_cell(iterj)
            if (yin(iterj) <= bbox_min(2) - buf*wy) then; sy = iterj
            else if (yin(iterj) <= bbox_max(2) + buf*wy) then; ey = iterj
            else; exit; end if
        end do
        do iterk = 1, nz_in
            wz = dz_in; if (present(dz_cell)) wz = dz_cell(iterk)
            if (zin(iterk) <= bbox_min(3) - buf*wz) then; sz = iterk
            else if (zin(iterk) <= bbox_max(3) + buf*wz) then; ez = iterk
            else; exit; end if
        end do

        if (procid == 0) then
            print *, "AABB x:", sx, xin(sx), "|", ex, xin(ex)
            print *, "AABB y:", sy, yin(sy), "|", ey, yin(ey)
            print *, "AABB z:", sz, zin(sz), "|", ez, zin(ez)
        end if
    end subroutine tagminmax

    ! Gather local SDF arrays from all ranks into a global array on the root rank
    subroutine gather_array(proc_id, num_procs, local_sdf, global_nx, global_ny, global_nz, &
                             start_indices, end_indices, chunk_sizes, s_y, e_y, s_z, e_z, &
                             input_scalarvalue, assembled_array, ierr_mpi)
        use mpi
        implicit none
        integer,  intent(in)  :: proc_id, num_procs, global_nx, global_ny, global_nz
        integer,  intent(in)  :: s_y, e_y, s_z, e_z
        integer,  intent(in),  dimension(0:) :: start_indices, end_indices, chunk_sizes
        real(dp), intent(in),  dimension(:,:,:) :: local_sdf
        real(dp), intent(in)  :: input_scalarvalue
        real(dp), allocatable, intent(inout), dimension(:,:,:) :: assembled_array
        integer,  intent(out) :: ierr_mpi
        real(dp), allocatable, dimension(:) :: sendbuf, recvbuf
        integer,  allocatable, dimension(:)  :: recvcounts, displs
        integer :: rank, local_count

        ! Build per-rank send counts and receive displacements
        allocate(recvcounts(0:num_procs-1), displs(0:num_procs-1))
        do rank = 0, num_procs-1
            recvcounts(rank) = chunk_sizes(rank) * global_ny * global_nz
        end do
        displs(0) = 0
        do rank = 1, num_procs-1
            displs(rank) = displs(rank-1) + recvcounts(rank-1)
        end do

        ! Pack local slab into a contiguous send buffer
        local_count = chunk_sizes(proc_id) * global_ny * global_nz
        allocate(sendbuf(local_count))
        sendbuf = reshape(local_sdf, [local_count])

        ! Root receives all slabs; allocate assembled_array if not already done
        if (proc_id == 0) then
            if (.not. allocated(assembled_array)) then
                allocate(assembled_array(global_nx, global_ny, global_nz))
            end if
            assembled_array = input_scalarvalue
            allocate(recvbuf(sum(recvcounts)))
        else
            allocate(recvbuf(1))   ! unused on non-root ranks
        end if

        call MPI_GATHERV(sendbuf, local_count, MPI_DOUBLE_PRECISION, &
                         recvbuf, recvcounts, displs, MPI_DOUBLE_PRECISION, &
                         0, MPI_COMM_WORLD, ierr_mpi)

        if (proc_id == 0) then
            do rank = 0, num_procs-1
                assembled_array(start_indices(rank):end_indices(rank), :, :) = &
                    reshape(recvbuf(displs(rank)+1 : displs(rank)+recvcounts(rank)), &
                            [chunk_sizes(rank), global_ny, global_nz])
            end do
            deallocate(recvbuf)
        end if

        deallocate(sendbuf, recvcounts, displs)
    end subroutine gather_array

    ! Write SDF to file, with ghost padding in x and wall-normal directions as expected by the solver's IBM
    !   reader.
    subroutine write_sdf_padded(outputfilename, interior)
        !
        ! Write SDF ghost-padded in x and the wall-normal direction, matching
        ! the solver IBM reader expectation of (nxm+2) x (n_wall+2) x n_span
        ! pages in solver coordinates.
        !
        implicit none
        character(len=*), intent(in) :: outputfilename
        real(dp), intent(in), dimension(:,:,:) :: interior
        integer :: d1, d2, d3, unit, j, k
        real(dp), allocatable :: padded(:,:,:)

        d1 = size(interior,1); d2 = size(interior,2); d3 = size(interior,3)

        ! Solver based padding on interior SDF
        allocate(padded(d1+2, d3+2, d2))
        padded = 0.0_dp
        do k = 1, d2           ! k loops spanwise = GenSDF dim 2
            do j = 1, d3       ! j loops wall-normal = GenSDF dim 3
                padded(2:d1+1, j+1, k) = interior(:, k, j)
            end do
        end do
        ! x ghost (dim 1)
        padded(1,    2:d3+1, :) = padded(2,    2:d3+1, :)
        padded(d1+2, 2:d3+1, :) = padded(d1+1, 2:d3+1, :)
        ! wall-normal ghost (dim 2)
        padded(:, 1,      :) = padded(:, 2,     :)
        padded(:, d3+2,   :) = padded(:, d3+1,  :)

        open(newunit=unit, file=trim(outputfilename), access='stream', status='replace', form='unformatted', &
             convert='big_endian')
        write(unit) padded
        close(unit)
        deallocate(padded)
    end subroutine write_sdf_padded

    ! Estimate memory usage based on input parameters
    subroutine estimated_memoryusage(n_procs, input_nfaces, input_nvertices, in_nx, in_ny, in_nz)
        implicit none
        integer, intent(in) :: n_procs, input_nfaces, input_nvertices, in_nx, in_ny, in_nz
        real(dp) :: total
        character(len=256) :: str
        total  = 6.0_dp*real(input_nfaces,    dp)*real(sizeof(in_nx),   dp)
        total  = total + 6.0_dp*real(input_nvertices, dp)*real(sizeof(total), dp)
        total  = total + 2.0_dp*real(in_nx*in_ny*in_nz, dp)*real(sizeof(total), dp)
        write(str,'(F5.2)') (n_procs*total)/1e9_dp
        print *, "-- Estimated minimum memory:", trim(str), " GiB(s)"
    end subroutine estimated_memoryusage

    ! Progress bar
    subroutine show_progress(current, total, width, start_time)
        use mpi, only : MPI_WTIME
        implicit none
        integer,  intent(in) :: current, total, width
        real(dp), intent(in) :: start_time
        integer  :: progress
        real(dp) :: percent, elapsed, remaining, now
        character(len=:), allocatable :: bar

        now = MPI_WTIME()
        elapsed = now - start_time
        if (current > 0) then
            remaining = elapsed * (real(total,dp) - real(current,dp)) / real(current,dp)
        else
            remaining = 0.0_dp
        end if
        if (total > 0) then
            percent  = real(current,dp) / real(total,dp) * 100.0_dp
            progress = int(real(current,dp) / real(total,dp) * width)
        else
            percent  = 0.0_dp
            progress = 0
        end if
        allocate(character(len=width) :: bar)
        bar = repeat('|', progress) // repeat(' ', width-progress)
        write(*,'(a,"|",a,"|",f6.2,"%",1x,i8,"/",i8," Elapsed:",f8.2,"s Remaining:",f8.2,"s")', &
              advance='no') char(13), bar, percent, current, total, elapsed, remaining
        flush(6)
        if (current == total) print *
        deallocate(bar)
    end subroutine show_progress

    ! LOGOOOO 
    subroutine printlogo()
        implicit none
        print *, " ░▒▓██████▓▒░░▒▓████████▓▒░▒▓███████▓▒░ ░▒▓███████▓▒░▒▓███████▓▒░░▒▓████████▓▒░ "
        print *, "░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░        "
        print *, "░▒▓█▓▒░      ░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░        "
        print *, "░▒▓█▓▒▒▓███▓▒░▒▓██████▓▒░ ░▒▓█▓▒░░▒▓█▓▒░░▒▓██████▓▒░░▒▓█▓▒░░▒▓█▓▒░▒▓██████▓▒░   "
        print *, "░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░      ░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░        "
        print *, "░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░      ░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░        "
        print *, " ░▒▓██████▓▒░░▒▓████████▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓███████▓▒░░▒▓███████▓▒░░▒▓█▓▒░        "
    end subroutine printlogo

end module utils_io
