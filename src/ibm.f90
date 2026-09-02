!> Ghost-Cell Immersed Boundary Method (Tseng & Ferziger 2003)
Module ibm

  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use mpi
  Use decomp, Only : z_halo_neighbors, x_periodic_partner

  Implicit None

  ! Wall velocity (no-slip by default; set non-zero for moving walls)
  Real(Int64) :: U_wall = 0d0
  Real(Int64) :: V_wall = 0d0
  Real(Int64) :: W_wall = 0d0

  ! Minimum number of fluid-side cells to look ahead for image point
  Integer(Int32), Parameter :: n_image_layers = 2

  ! Accumulators for Method 1 IBM force (summed over 6 IBM applications/step)
  Real(Int64) :: ibm_Fx_acc = 0d0
  Real(Int64) :: ibm_Fy_acc = 0d0
  Real(Int64) :: ibm_Fz_acc = 0d0

Contains

  !  Read precomputed SDF and populate ghost-cell lists for U, V, W.
  !  Must be called after grid setup (xg/yg/zg must be initialised).
  Subroutine setup_ibm

    Integer(Int32) :: n_ghost_u_global, n_ghost_v_global, n_ghost_w_global

    If ( myid==0 ) Write(*,*) 'IBM: reading cell-centre SDF from ', Trim(ibm_sdf_file), '...'
    Call read_phi_from_sdf_file

    If ( myid==0 ) Write(*,*) 'IBM: building ghost-cell lists for U, V, W...'
    Call build_ghost_list_u
    Call build_ghost_list_v
    Call build_ghost_list_w

    If ( myid==0 ) Write(*,*) 'IBM: building cell-centre ghost list for pressure integration...'
    Call build_ghost_list_cc

    If ( myid==0 ) Then
       Write(*,*) 'IBM: ghost cells (rank 0) — U:', n_ghost_u, ' V:', n_ghost_v, ' W:', n_ghost_w, &
                  ' CC:', n_ghost_cc
    End If

    ! Global sum across all ranks for diagnostic
    Call MPI_Reduce(n_ghost_u, n_ghost_u_global, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    Call MPI_Reduce(n_ghost_v, n_ghost_v_global, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    Call MPI_Reduce(n_ghost_w, n_ghost_w_global, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    If ( myid==0 ) Then
       Write(*,*) 'IBM: ghost cells (GLOBAL) — U:', n_ghost_u_global, ' V:', n_ghost_v_global, ' W:', n_ghost_w_global
    End If

    ! Device-resident IBM data: phi (read once above) and the ghost-cell lists (built once above), never modified afterward
    !$acc enter data copyin(phi)
    !$acc enter data copyin(ghost_u_idx,ghost_u_wgt,ghost_u_img,ghost_u_ref,ghost_u_nrm,ghost_u_yref,ghost_u_dGB,ghost_u_objid)
    !$acc enter data copyin(ghost_v_idx,ghost_v_wgt,ghost_v_img,ghost_v_ref,ghost_v_nrm,ghost_v_yref,ghost_v_dGB,ghost_v_objid)
    !$acc enter data copyin(ghost_w_idx,ghost_w_wgt,ghost_w_img,ghost_w_ref,ghost_w_nrm,ghost_w_yref,ghost_w_dGB,ghost_w_objid)
    ! Cell-centre ghost list: used both by the host-only Method-2 force diagnostic and (when boussinesq_flag>=1) by the device-resident apply_ghost_cell_ibm_scalar
    !$acc enter data copyin(ghost_cc_idx,ghost_cc_wgt_cc,ghost_cc_img_cc,ghost_cc_objid)

  End Subroutine setup_ibm

  !> Exchange phi z-ghost planes between row-adjacent MPI ranks; domain-boundary Neumann values from the caller are left untouched
  Subroutine exchange_phi_ghost_planes

    Real   (Int64) :: buf_s(nxg,nyg), buf_r(nxg,nyg)
    Integer(Int32) :: up, down

    Call z_halo_neighbors(up, down)

    ! Pass 1: send k=nzg-1 towards +z; receive from -z into k=1.
    buf_s = phi(:,:,nzg-1)
    Call MPI_Sendrecv( buf_s, nxg*nyg, MPI_real8, up,   0, &
                       buf_r, nxg*nyg, MPI_real8, down, 0, &
                       MPI_COMM_WORLD, istat, ierr )
    If ( down /= MPI_PROC_NULL ) phi(:,:,1) = buf_r

    ! Pass 2: send k=2 towards -z; receive from +z into k=nzg.
    buf_s = phi(:,:,2)
    Call MPI_Sendrecv( buf_s, nxg*nyg, MPI_real8, down, 0, &
                       buf_r, nxg*nyg, MPI_real8, up,   0, &
                       MPI_COMM_WORLD, istat, ierr )
    If ( up /= MPI_PROC_NULL ) phi(:,:,nzg) = buf_r

  End Subroutine exchange_phi_ghost_planes

  !> Read a distributed cell-centre scalar field from file: (nxg_global,nyg_global,nzm_global) Real(8), column-major, big-endian (x,y already ghosted in-file, z interior-only, same convention as xg_global/kg-based reads elsewhere)
  Subroutine read_distributed_scalar_field(filename, field)

    Character(*), Intent(In) :: filename
    Real(Int64), Dimension(nxg,nyg,nzg), Intent(InOut) :: field

    Integer(Int32) :: iproc, nxge_r, nzge_r, n_interior
    Integer(Int32) :: sdf_unit
    Real   (Int64), Allocatable :: global_field(:,:,:), tmp_read(:), send_buf(:,:,:)

    ! Rank iproc owns global x-columns ig1_global(iproc):ig2_global(iproc) (already ghosted in-file)
    ! and interior z-planes kg1_global(iproc):kg2_global(iproc)-2.
    If ( myid==0 ) Then

       Open(newunit=sdf_unit, file=Trim(filename), access='stream', form='unformatted', action='read', convert='big_endian')
       Allocate( tmp_read(Int(nxg_global, Int64) * Int(nyg_global, Int64) * Int(nzm_global, Int64)) )
       Read(sdf_unit) tmp_read
       Close(sdf_unit)
       Allocate( global_field(nxg_global, nyg_global, nzm_global) )
       global_field = Reshape(tmp_read, [nxg_global, nyg_global, nzm_global])
       Deallocate(tmp_read)

       Do iproc = 0, nprocs-1
          nxge_r     = ig2_global(iproc) - ig1_global(iproc) + 1
          nzge_r     = kg2_global(iproc) - kg1_global(iproc) + 1
          n_interior = nzge_r - 2
          Allocate( send_buf(nxge_r, nyg_global, n_interior) )
          send_buf = global_field( ig1_global(iproc):ig2_global(iproc), :, kg1_global(iproc):kg2_global(iproc)-2 )
          If ( iproc == 0 ) Then
             field(:,:,2:nzge_r-1) = send_buf
          Else
             Call MPI_Send(send_buf, nxge_r*nyg_global*n_interior, MPI_real8, iproc, iproc, MPI_COMM_WORLD, ierr)
          End If
          Deallocate(send_buf)
       End Do
       Deallocate(global_field)

    Else
       Call MPI_Recv(field(:,:,2:nzg-1), nxg*nyg*(nzg-2), MPI_real8, 0, myid, MPI_COMM_WORLD, istat, ierr)
    End If

  End Subroutine read_distributed_scalar_field

  !> Read cell-centre SDF from ibm_sdf_file and derive Umask_cc from sign(phi); optionally reads the per-solid object-ID field from ibm_objid_file
  Subroutine read_phi_from_sdf_file

    Integer(Int32) :: i, j, k
    Logical        :: is_first_x, is_last_x
    Integer(Int32) :: partner_x

    If (Allocated(phi))       Deallocate(phi)
    If (Allocated(Umask_cc))  Deallocate(Umask_cc)
    If (Allocated(ibm_obj_id)) Deallocate(ibm_obj_id)
    Allocate ( phi       (nxg, nyg, nzg) )
    Allocate ( Umask_cc  (nxg, nyg, nzg) )
    Allocate ( ibm_obj_id(nxg, nyg, nzg) )
    phi        = 0d0
    Umask_cc   = 0d0
    ibm_obj_id = 0d0

    Call read_distributed_scalar_field(ibm_sdf_file, phi)

    ! Apply Neumann (zero-gradient) ghost BCs for phi at domain boundaries; x=1/nxg is an
    ! inter-rank seam (already correctly populated straight from the global file slice) on
    ! any rank that doesn't own the true global x edge, so skip it there
    Call x_periodic_partner(is_first_x, is_last_x, partner_x)
    If ( is_first_x ) phi(1,:,:)   = phi(2,:,:)
    If ( is_last_x  ) phi(nxg,:,:) = phi(nxg-1,:,:)
    phi(:,1,:)   = phi(:,2,:)
    phi(:,nyg,:) = phi(:,nyg-1,:)
    phi(:,:,1)   = phi(:,:,2)
    phi(:,:,nzg) = phi(:,:,nzg-1)
    ! Overwrite interior-rank z ghost planes with actual neighbour values
    Call exchange_phi_ghost_planes

    ! Derive Umask_cc from sign of phi (positive = fluid)
    Do k = 1, nzg
       Do j = 1, nyg
          Do i = 1, nxg
             If ( phi(i,j,k) > 0d0 ) Then
                Umask_cc(i,j,k) = 1d0
             Else
                Umask_cc(i,j,k) = 0d0
             End If
          End Do
       End Do
    End Do

    ! Optional per-solid object-ID field (companion to phi, written by GenSDF's .list manifest mode)
    If ( Len_Trim(ibm_objid_file) > 0 ) Then
       If ( myid==0 ) Write(*,*) 'IBM: reading per-object ID field from ', Trim(ibm_objid_file), '...'
       Call read_distributed_scalar_field(ibm_objid_file, ibm_obj_id)
    End If

  End Subroutine read_phi_from_sdf_file

  !  Build ghost-cell list for U (defined at x-faces)
  Subroutine build_ghost_list_u

    Integer(Int32) :: i, j, k, ng, nd
    Integer(Int32) :: ii, jj, kk, ir, jr, kr
    Real   (Int64) :: nx_, ny_, nz_
    Real   (Int64) :: xGc, yGc, zGc, xB, yB, zB, xI, yI, zI
    Real   (Int64) :: dGB, dGI
    Real   (Int64) :: yref_

    ! Pass 1: count ghost cells (no storage) — avoids any fixed-size estimate
    ng = 0;  nd = 0
    Do k = 2, nzg-1
       Do j = 2, nyg-1
          Do i = 2, nx-1
             If ( 0.5d0*(phi(i,j,k)+phi(i+1,j,k)) < 0d0 ) Then
                If ( 0.5d0*(phi(i-1,j,k)+phi(i,  j,  k  )) >= 0d0 .Or. &
                     0.5d0*(phi(i+1,j,k)+phi(i+2,j,  k  )) >= 0d0 .Or. &
                     0.5d0*(phi(i,  j-1,k)+phi(i+1,j-1,k)) >= 0d0 .Or. &
                     0.5d0*(phi(i,  j+1,k)+phi(i+1,j+1,k)) >= 0d0 .Or. &
                     0.5d0*(phi(i,  j,  k-1)+phi(i+1,j,k-1)) >= 0d0 .Or. &
                     0.5d0*(phi(i,  j,  k+1)+phi(i+1,j,k+1)) >= 0d0 ) Then
                   Call compute_normal_at_face_u(i,j,k, nx_,ny_,nz_)
                   dGB = Abs( 0.5d0*(phi(i,j,k)+phi(i+1,j,k)) )
                   xI  = x(i)  + 2d0*dGB*nx_
                   yI  = yg(j) + 2d0*dGB*ny_
                   zI  = zg(k) + 2d0*dGB*nz_
                   Call find_stencil_u(xI, yI, zI, ii, jj, kk)
                   If ( ii>=2 .And. ii<=nx-1 .And. jj>=2 .And. jj<=nyg-1 .And. &
                        kk>=2 .And. kk<=nzg-1 ) Then
                      If ( Umask_cc(ii,  jj,  kk  ) > 0.5d0 .And. Umask_cc(ii+1,jj,  kk  ) > 0.5d0 .And. &
                           Umask_cc(ii,  jj+1,kk  ) > 0.5d0 .And. Umask_cc(ii+1,jj+1,kk  ) > 0.5d0 .And. &
                           Umask_cc(ii,  jj,  kk+1) > 0.5d0 .And. Umask_cc(ii+1,jj,  kk+1) > 0.5d0 .And. &
                           Umask_cc(ii,  jj+1,kk+1) > 0.5d0 .And. Umask_cc(ii+1,jj+1,kk+1) > 0.5d0 ) Then
                         ng = ng + 1
                      Else
                         nd = nd + 1   ! image stencil clips solid — concave corner
                      End If
                   Else
                      nd = nd + 1      ! image outside domain bounds
                   End If
                End If
             End If
          End Do
       End Do
    End Do
    If ( nd > 0 .And. myid == 0 ) &
       Write(*,'(A,I0,A)') '[IBM] build_ghost_list_u: ', nd, &
          ' ghost cell(s) dropped (image in solid or outside domain)'

    ! Allocate final arrays with exact count
    n_ghost_u = ng
    Allocate ( ghost_u_idx(3,n_ghost_u), ghost_u_img(3,n_ghost_u), ghost_u_wgt(8,n_ghost_u) )
    Allocate ( ghost_u_nrm(3,n_ghost_u), ghost_u_yref(n_ghost_u)  )
    Allocate ( ghost_u_ref(3,n_ghost_u) )
    Allocate ( ghost_u_dGB(  n_ghost_u) )
    Allocate ( ghost_u_xB(3, n_ghost_u) )
    Allocate ( ghost_u_img_cc(3, n_ghost_u) )
    Allocate ( ghost_u_wgt_cc(8, n_ghost_u) )
    Allocate ( ghost_u_objid(n_ghost_u) )

    ! Pass 2: identical traversal, fill ghost arrays directly
    ng = 0
    Do k = 2, nzg-1
       Do j = 2, nyg-1
          Do i = 2, nx-1
             If ( 0.5d0*(phi(i,j,k)+phi(i+1,j,k)) < 0d0 ) Then
                If ( 0.5d0*(phi(i-1,j,k)+phi(i,  j,  k  )) >= 0d0 .Or. &
                     0.5d0*(phi(i+1,j,k)+phi(i+2,j,  k  )) >= 0d0 .Or. &
                     0.5d0*(phi(i,  j-1,k)+phi(i+1,j-1,k)) >= 0d0 .Or. &
                     0.5d0*(phi(i,  j+1,k)+phi(i+1,j+1,k)) >= 0d0 .Or. &
                     0.5d0*(phi(i,  j,  k-1)+phi(i+1,j,k-1)) >= 0d0 .Or. &
                     0.5d0*(phi(i,  j,  k+1)+phi(i+1,j,k+1)) >= 0d0 ) Then

                   xGc = x(i)
                   yGc = yg(j)
                   zGc = zg(k)

                   Call compute_normal_at_face_u(i,j,k, nx_,ny_,nz_)

                   dGB = Abs( 0.5d0*(phi(i,j,k)+phi(i+1,j,k)) )

                   xB = xGc + dGB*nx_
                   yB = yGc + dGB*ny_
                   zB = zGc + dGB*nz_

                   dGI = 2d0*dGB
                   xI  = xGc + dGI*nx_
                   yI  = yGc + dGI*ny_
                   zI  = zGc + dGI*nz_

                   Call find_stencil_u(xI, yI, zI, ii, jj, kk)

                   If ( ii>=2 .And. ii<=nx-1 .And. jj>=2 .And. jj<=nyg-1 .And. &
                        kk>=2 .And. kk<=nzg-1 ) Then
                      If ( Umask_cc(ii,  jj,  kk  ) > 0.5d0 .And. Umask_cc(ii+1,jj,  kk  ) > 0.5d0 .And. &
                           Umask_cc(ii,  jj+1,kk  ) > 0.5d0 .And. Umask_cc(ii+1,jj+1,kk  ) > 0.5d0 .And. &
                           Umask_cc(ii,  jj,  kk+1) > 0.5d0 .And. Umask_cc(ii+1,jj,  kk+1) > 0.5d0 .And. &
                           Umask_cc(ii,  jj+1,kk+1) > 0.5d0 .And. Umask_cc(ii+1,jj+1,kk+1) > 0.5d0 ) Then

                         ng = ng + 1
                         ghost_u_idx(1,ng) = i
                         ghost_u_idx(2,ng) = j
                         ghost_u_idx(3,ng) = k
                         ghost_u_img(1,ng) = ii
                         ghost_u_img(2,ng) = jj
                         ghost_u_img(3,ng) = kk
                         Call trilinear_weights_u(xI, yI, zI, ii, jj, kk, ghost_u_wgt(1:8,ng))
                         ghost_u_nrm(1,ng)  = nx_
                         ghost_u_nrm(2,ng)  = ny_
                         ghost_u_nrm(3,ng)  = nz_
                         Call find_stencil_centre(xB + Real(n_image_layers,8)*dymin*nx_, &
                                                  yB + Real(n_image_layers,8)*dymin*ny_, &
                                                  zB + Real(n_image_layers,8)*dymin*nz_, &
                                                  ir, jr, kr)
                         yref_ = Max( Abs( (xg(ir)-xB)*nx_ + (yg(jr)-yB)*ny_ + (zg(kr)-zB)*nz_ ), 1d-14 )
                         ghost_u_yref(ng)   = yref_
                         ghost_u_ref(1,ng) = ir
                         ghost_u_ref(2,ng) = jr
                         ghost_u_ref(3,ng) = kr

                      End If
                   End If

                End If
             End If
          End Do
       End Do
    End Do

    ! Store dGB, boundary-point coords, and precomputed pressure stencil
    Do ng = 1, n_ghost_u
       i = ghost_u_idx(1,ng);  j = ghost_u_idx(2,ng);  k = ghost_u_idx(3,ng)
       ghost_u_objid(ng) = Min(Max(Nint(ibm_obj_id(i,j,k)), 0), max_ibm_objects)
       ghost_u_dGB(ng)   = Abs( 0.5d0*(phi(i,j,k)+phi(i+1,j,k)) )
       ghost_u_xB(1,ng)  = x (i)  + ghost_u_dGB(ng)*ghost_u_nrm(1,ng)
       ghost_u_xB(2,ng)  = yg(j)  + ghost_u_dGB(ng)*ghost_u_nrm(2,ng)
       ghost_u_xB(3,ng)  = zg(k)  + ghost_u_dGB(ng)*ghost_u_nrm(3,ng)
       ! Cell-centre stencil at image point I = xB + dGB*nrm = xG + 2*dGB*nrm
       ! (same xI as velocity image point, on cell-centre grid for P interpolation)
       Call find_stencil_centre( ghost_u_xB(1,ng) + ghost_u_dGB(ng)*ghost_u_nrm(1,ng), &
                                 ghost_u_xB(2,ng) + ghost_u_dGB(ng)*ghost_u_nrm(2,ng), &
                                 ghost_u_xB(3,ng) + ghost_u_dGB(ng)*ghost_u_nrm(3,ng), &
                                 ii, jj, kk )
       Call trilinear_weights( ghost_u_xB(1,ng) + ghost_u_dGB(ng)*ghost_u_nrm(1,ng), &
                               ghost_u_xB(2,ng) + ghost_u_dGB(ng)*ghost_u_nrm(2,ng), &
                               ghost_u_xB(3,ng) + ghost_u_dGB(ng)*ghost_u_nrm(3,ng), &
                               ii, jj, kk, ghost_u_wgt_cc(1:8,ng) )
       ghost_u_img_cc(1,ng) = ii;  ghost_u_img_cc(2,ng) = jj;  ghost_u_img_cc(3,ng) = kk
    End Do

  End Subroutine build_ghost_list_u

  !  Build ghost-cell list for V (y-faces)
  Subroutine build_ghost_list_v

    Integer(Int32) :: i, j, k, ng, nd
    Integer(Int32) :: ii, jj, kk, ir, jr, kr
    Real   (Int64) :: nx_, ny_, nz_
    Real   (Int64) :: xGc, yGc, zGc, dGB, xB, yB, zB, xI, yI, zI, yref_

    ! Pass 1: count ghost cells (no storage)
    ng = 0;  nd = 0
    Do k = 2, nzg-1
       Do j = 2, ny-1
          Do i = 2, nxg-1
             If ( 0.5d0*(phi(i,j,k)+phi(i,j+1,k)) < 0d0 ) Then
                If ( 0.5d0*(phi(i-1,j,k)+phi(i-1,j+1,k)) >= 0d0 .Or. &
                     0.5d0*(phi(i+1,j,k)+phi(i+1,j+1,k)) >= 0d0 .Or. &
                     0.5d0*(phi(i,j-1,k)+phi(i,j,  k))   >= 0d0 .Or. &
                     0.5d0*(phi(i,j+1,k)+phi(i,j+2,k))   >= 0d0 .Or. &
                     0.5d0*(phi(i,j,k-1)+phi(i,j+1,k-1)) >= 0d0 .Or. &
                     0.5d0*(phi(i,j,k+1)+phi(i,j+1,k+1)) >= 0d0 ) Then
                   Call compute_normal_at_face_v(i,j,k, nx_,ny_,nz_)
                   dGB = Abs( 0.5d0*(phi(i,j,k)+phi(i,j+1,k)) )
                   xI = xg(i) + 2d0*dGB*nx_
                   yI = y(j)  + 2d0*dGB*ny_
                   zI = zg(k) + 2d0*dGB*nz_
                   Call find_stencil_v(xI, yI, zI, ii, jj, kk)
                   If ( ii>=2 .And. ii<=nxg-1 .And. jj>=2 .And. jj<=ny-1 .And. &
                        kk>=2 .And. kk<=nzg-1 ) Then
                      If ( Umask_cc(ii,  jj,  kk  ) > 0.5d0 .And. Umask_cc(ii+1,jj,  kk  ) > 0.5d0 .And. &
                           Umask_cc(ii,  jj+1,kk  ) > 0.5d0 .And. Umask_cc(ii+1,jj+1,kk  ) > 0.5d0 .And. &
                           Umask_cc(ii,  jj,  kk+1) > 0.5d0 .And. Umask_cc(ii+1,jj,  kk+1) > 0.5d0 .And. &
                           Umask_cc(ii,  jj+1,kk+1) > 0.5d0 .And. Umask_cc(ii+1,jj+1,kk+1) > 0.5d0 ) Then
                         ng = ng + 1
                      Else
                         nd = nd + 1   ! image stencil clips solid — concave corner
                      End If
                   Else
                      nd = nd + 1      ! image outside domain bounds
                   End If
                End If
             End If
          End Do
       End Do
    End Do
    If ( nd > 0 .And. myid == 0 ) &
       Write(*,'(A,I0,A)') '[IBM] build_ghost_list_v: ', nd, &
          ' ghost cell(s) dropped (image in solid or outside domain)'

    ! Allocate final arrays with exact count
    n_ghost_v = ng
    Allocate ( ghost_v_idx(3,n_ghost_v), ghost_v_img(3,n_ghost_v), ghost_v_wgt(8,n_ghost_v) )
    Allocate ( ghost_v_nrm(3,n_ghost_v), ghost_v_yref(n_ghost_v)  )
    Allocate ( ghost_v_ref(3,n_ghost_v) )
    Allocate ( ghost_v_dGB(  n_ghost_v) )
    Allocate ( ghost_v_xB(3, n_ghost_v) )
    Allocate ( ghost_v_img_cc(3, n_ghost_v) )
    Allocate ( ghost_v_wgt_cc(8, n_ghost_v) )
    Allocate ( ghost_v_objid(n_ghost_v) )

    ! Pass 2: identical traversal, fill ghost arrays directly
    ng = 0
    Do k = 2, nzg-1
       Do j = 2, ny-1
          Do i = 2, nxg-1
             If ( 0.5d0*(phi(i,j,k)+phi(i,j+1,k)) < 0d0 ) Then
                If ( 0.5d0*(phi(i-1,j,k)+phi(i-1,j+1,k)) >= 0d0 .Or. &
                     0.5d0*(phi(i+1,j,k)+phi(i+1,j+1,k)) >= 0d0 .Or. &
                     0.5d0*(phi(i,j-1,k)+phi(i,j,  k))   >= 0d0 .Or. &
                     0.5d0*(phi(i,j+1,k)+phi(i,j+2,k))   >= 0d0 .Or. &
                     0.5d0*(phi(i,j,k-1)+phi(i,j+1,k-1)) >= 0d0 .Or. &
                     0.5d0*(phi(i,j,k+1)+phi(i,j+1,k+1)) >= 0d0 ) Then

                   xGc = xg(i)
                   yGc = y(j)
                   zGc = zg(k)

                   Call compute_normal_at_face_v(i,j,k, nx_,ny_,nz_)

                   dGB = Abs( 0.5d0*(phi(i,j,k)+phi(i,j+1,k)) )
                   xB = xGc + dGB*nx_;  yB = yGc + dGB*ny_;  zB = zGc + dGB*nz_
                   xI = xGc + 2d0*dGB*nx_;  yI = yGc + 2d0*dGB*ny_;  zI = zGc + 2d0*dGB*nz_

                   Call find_stencil_v(xI, yI, zI, ii, jj, kk)
                   If ( ii>=2 .And. ii<=nxg-1 .And. jj>=2 .And. jj<=ny-1 .And. &
                        kk>=2 .And. kk<=nzg-1 ) Then
                      If ( Umask_cc(ii,  jj,  kk  ) > 0.5d0 .And. Umask_cc(ii+1,jj,  kk  ) > 0.5d0 .And. &
                           Umask_cc(ii,  jj+1,kk  ) > 0.5d0 .And. Umask_cc(ii+1,jj+1,kk  ) > 0.5d0 .And. &
                           Umask_cc(ii,  jj,  kk+1) > 0.5d0 .And. Umask_cc(ii+1,jj,  kk+1) > 0.5d0 .And. &
                           Umask_cc(ii,  jj+1,kk+1) > 0.5d0 .And. Umask_cc(ii+1,jj+1,kk+1) > 0.5d0 ) Then
                         ng = ng + 1
                         ghost_v_idx(1,ng)=i; ghost_v_idx(2,ng)=j; ghost_v_idx(3,ng)=k
                         ghost_v_img(1,ng)=ii; ghost_v_img(2,ng)=jj; ghost_v_img(3,ng)=kk
                         Call trilinear_weights_v(xI, yI, zI, ii, jj, kk, ghost_v_wgt(1:8,ng))
                         ghost_v_nrm(:,ng)=[nx_,ny_,nz_]
                         Call find_stencil_centre(xB+Real(n_image_layers,8)*dymin*nx_, &
                                                  yB+Real(n_image_layers,8)*dymin*ny_, &
                                                  zB+Real(n_image_layers,8)*dymin*nz_, &
                                                  ir,jr,kr)
                         ghost_v_yref(ng) = Max( Abs( (xg(ir)-xB)*nx_ + (yg(jr)-yB)*ny_ + (zg(kr)-zB)*nz_ ), 1d-14 )
                         ghost_v_ref(:,ng)=[ir,jr,kr]
                      End If
                   End If
                End If
             End If
          End Do
       End Do
    End Do

    Do ng = 1, n_ghost_v
       i = ghost_v_idx(1,ng);  j = ghost_v_idx(2,ng);  k = ghost_v_idx(3,ng)
       ghost_v_objid(ng) = Min(Max(Nint(ibm_obj_id(i,j,k)), 0), max_ibm_objects)
       ghost_v_dGB(ng)   = Abs( 0.5d0*(phi(i,j,k)+phi(i,j+1,k)) )
       ghost_v_xB(1,ng)  = xg(i)  + ghost_v_dGB(ng)*ghost_v_nrm(1,ng)
       ghost_v_xB(2,ng)  = y (j)  + ghost_v_dGB(ng)*ghost_v_nrm(2,ng)
       ghost_v_xB(3,ng)  = zg(k)  + ghost_v_dGB(ng)*ghost_v_nrm(3,ng)
       Call find_stencil_centre( ghost_v_xB(1,ng) + ghost_v_dGB(ng)*ghost_v_nrm(1,ng), &
                                 ghost_v_xB(2,ng) + ghost_v_dGB(ng)*ghost_v_nrm(2,ng), &
                                 ghost_v_xB(3,ng) + ghost_v_dGB(ng)*ghost_v_nrm(3,ng), &
                                 ii, jj, kk )
       Call trilinear_weights( ghost_v_xB(1,ng) + ghost_v_dGB(ng)*ghost_v_nrm(1,ng), &
                               ghost_v_xB(2,ng) + ghost_v_dGB(ng)*ghost_v_nrm(2,ng), &
                               ghost_v_xB(3,ng) + ghost_v_dGB(ng)*ghost_v_nrm(3,ng), &
                               ii, jj, kk, ghost_v_wgt_cc(1:8,ng) )
       ghost_v_img_cc(1,ng) = ii;  ghost_v_img_cc(2,ng) = jj;  ghost_v_img_cc(3,ng) = kk
    End Do

  End Subroutine build_ghost_list_v

  !  Build ghost-cell list for W (z-faces)
  Subroutine build_ghost_list_w

    Integer(Int32) :: i, j, k, ng, nd
    Integer(Int32) :: ii, jj, kk, ir, jr, kr
    Real   (Int64) :: nx_, ny_, nz_
    Real   (Int64) :: xGc, yGc, zGc, dGB, xB, yB, zB, xI, yI, zI, yref_

    ! Pass 1: count ghost cells (no storage)
    ng = 0;  nd = 0
    Do k = 2, nz-1
       Do j = 2, nyg-1
          Do i = 2, nxg-1
             If ( 0.5d0*(phi(i,j,k)+phi(i,j,k+1)) < 0d0 ) Then
                If ( 0.5d0*(phi(i-1,j,k)+phi(i-1,j,k+1)) >= 0d0 .Or. &
                     0.5d0*(phi(i+1,j,k)+phi(i+1,j,k+1)) >= 0d0 .Or. &
                     0.5d0*(phi(i,j-1,k)+phi(i,j-1,k+1)) >= 0d0 .Or. &
                     0.5d0*(phi(i,j+1,k)+phi(i,j+1,k+1)) >= 0d0 .Or. &
                     0.5d0*(phi(i,j,k-1)+phi(i,j,k))               >= 0d0 .Or. &
                     0.5d0*(phi(i,j,k+1)+phi(i,j,min(k+2,nzg)))   >= 0d0 ) Then
                   Call compute_normal_at_face_w(i,j,k, nx_,ny_,nz_)
                   dGB = Abs( 0.5d0*(phi(i,j,k)+phi(i,j,k+1)) )
                   xI = xg(i) + 2d0*dGB*nx_
                   yI = yg(j) + 2d0*dGB*ny_
                   zI = z(k)  + 2d0*dGB*nz_
                   Call find_stencil_w(xI, yI, zI, ii, jj, kk)
                   If ( ii>=2 .And. ii<=nxg-1 .And. jj>=2 .And. jj<=nyg-1 .And. &
                        kk>=2 .And. kk<=nz-1 ) Then
                      If ( Umask_cc(ii,  jj,  kk  ) > 0.5d0 .And. Umask_cc(ii+1,jj,  kk  ) > 0.5d0 .And. &
                           Umask_cc(ii,  jj+1,kk  ) > 0.5d0 .And. Umask_cc(ii+1,jj+1,kk  ) > 0.5d0 .And. &
                           Umask_cc(ii,  jj,  kk+1) > 0.5d0 .And. Umask_cc(ii+1,jj,  kk+1) > 0.5d0 .And. &
                           Umask_cc(ii,  jj+1,kk+1) > 0.5d0 .And. Umask_cc(ii+1,jj+1,kk+1) > 0.5d0 ) Then
                         ng = ng + 1
                      Else
                         nd = nd + 1   ! image stencil clips solid — concave corner
                      End If
                   Else
                      nd = nd + 1      ! image outside domain bounds
                   End If
                End If
             End If
          End Do
       End Do
    End Do
    If ( nd > 0 .And. myid == 0 ) &
       Write(*,'(A,I0,A)') '[IBM] build_ghost_list_w: ', nd, &
          ' ghost cell(s) dropped (image in solid or outside domain)'

    ! Allocate final arrays with exact count
    n_ghost_w = ng
    Allocate ( ghost_w_idx(3,n_ghost_w), ghost_w_img(3,n_ghost_w), ghost_w_wgt(8,n_ghost_w) )
    Allocate ( ghost_w_nrm(3,n_ghost_w), ghost_w_yref(n_ghost_w)  )
    Allocate ( ghost_w_ref(3,n_ghost_w) )
    Allocate ( ghost_w_dGB(  n_ghost_w) )
    Allocate ( ghost_w_xB(3, n_ghost_w) )
    Allocate ( ghost_w_img_cc(3, n_ghost_w) )
    Allocate ( ghost_w_wgt_cc(8, n_ghost_w) )
    Allocate ( ghost_w_objid(n_ghost_w) )

    ! Pass 2: identical traversal, fill ghost arrays directly
    ng = 0
    Do k = 2, nz-1
       Do j = 2, nyg-1
          Do i = 2, nxg-1
             If ( 0.5d0*(phi(i,j,k)+phi(i,j,k+1)) < 0d0 ) Then
                If ( 0.5d0*(phi(i-1,j,k)+phi(i-1,j,k+1)) >= 0d0 .Or. &
                     0.5d0*(phi(i+1,j,k)+phi(i+1,j,k+1)) >= 0d0 .Or. &
                     0.5d0*(phi(i,j-1,k)+phi(i,j-1,k+1)) >= 0d0 .Or. &
                     0.5d0*(phi(i,j+1,k)+phi(i,j+1,k+1)) >= 0d0 .Or. &
                     0.5d0*(phi(i,j,k-1)+phi(i,j,k))               >= 0d0 .Or. &
                     0.5d0*(phi(i,j,k+1)+phi(i,j,min(k+2,nzg)))   >= 0d0 ) Then

                   xGc = xg(i)
                   yGc = yg(j)
                   zGc = z(k)

                   Call compute_normal_at_face_w(i,j,k, nx_,ny_,nz_)

                   dGB = Abs( 0.5d0*(phi(i,j,k)+phi(i,j,k+1)) )
                   xB = xGc+dGB*nx_;  yB = yGc+dGB*ny_;  zB = zGc+dGB*nz_
                   xI = xGc+2d0*dGB*nx_;  yI = yGc+2d0*dGB*ny_;  zI = zGc+2d0*dGB*nz_

                   Call find_stencil_w(xI, yI, zI, ii, jj, kk)
                   If ( ii>=2 .And. ii<=nxg-1 .And. jj>=2 .And. jj<=nyg-1 .And. &
                        kk>=2 .And. kk<=nz-1 ) Then
                      If ( Umask_cc(ii,  jj,  kk  ) > 0.5d0 .And. Umask_cc(ii+1,jj,  kk  ) > 0.5d0 .And. &
                           Umask_cc(ii,  jj+1,kk  ) > 0.5d0 .And. Umask_cc(ii+1,jj+1,kk  ) > 0.5d0 .And. &
                           Umask_cc(ii,  jj,  kk+1) > 0.5d0 .And. Umask_cc(ii+1,jj,  kk+1) > 0.5d0 .And. &
                           Umask_cc(ii,  jj+1,kk+1) > 0.5d0 .And. Umask_cc(ii+1,jj+1,kk+1) > 0.5d0 ) Then
                         ng = ng + 1
                         ghost_w_idx(1,ng)=i; ghost_w_idx(2,ng)=j; ghost_w_idx(3,ng)=k
                         ghost_w_img(1,ng)=ii; ghost_w_img(2,ng)=jj; ghost_w_img(3,ng)=kk
                         Call trilinear_weights_w(xI, yI, zI, ii, jj, kk, ghost_w_wgt(1:8,ng))
                         ghost_w_nrm(:,ng)=[nx_,ny_,nz_]
                         Call find_stencil_centre(xB+Real(n_image_layers,8)*dymin*nx_, &
                                                  yB+Real(n_image_layers,8)*dymin*ny_, &
                                                  zB+Real(n_image_layers,8)*dymin*nz_, &
                                                  ir,jr,kr)
                         ghost_w_yref(ng) = Max( Abs( (xg(ir)-xB)*nx_ + (yg(jr)-yB)*ny_ + (zg(kr)-zB)*nz_ ), 1d-14 )
                         ghost_w_ref(:,ng)=[ir,jr,kr]
                      End If
                   End If
                End If
             End If
          End Do
       End Do
    End Do

    Do ng = 1, n_ghost_w
       i = ghost_w_idx(1,ng);  j = ghost_w_idx(2,ng);  k = ghost_w_idx(3,ng)
       ghost_w_objid(ng) = Min(Max(Nint(ibm_obj_id(i,j,k)), 0), max_ibm_objects)
       ghost_w_dGB(ng)   = Abs( 0.5d0*(phi(i,j,k)+phi(i,j,k+1)) )
       ghost_w_xB(1,ng)  = xg(i)  + ghost_w_dGB(ng)*ghost_w_nrm(1,ng)
       ghost_w_xB(2,ng)  = yg(j)  + ghost_w_dGB(ng)*ghost_w_nrm(2,ng)
       ghost_w_xB(3,ng)  = z (k)  + ghost_w_dGB(ng)*ghost_w_nrm(3,ng)
       Call find_stencil_centre( ghost_w_xB(1,ng) + ghost_w_dGB(ng)*ghost_w_nrm(1,ng), &
                                 ghost_w_xB(2,ng) + ghost_w_dGB(ng)*ghost_w_nrm(2,ng), &
                                 ghost_w_xB(3,ng) + ghost_w_dGB(ng)*ghost_w_nrm(3,ng), &
                                 ii, jj, kk )
       Call trilinear_weights( ghost_w_xB(1,ng) + ghost_w_dGB(ng)*ghost_w_nrm(1,ng), &
                               ghost_w_xB(2,ng) + ghost_w_dGB(ng)*ghost_w_nrm(2,ng), &
                               ghost_w_xB(3,ng) + ghost_w_dGB(ng)*ghost_w_nrm(3,ng), &
                               ii, jj, kk, ghost_w_wgt_cc(1:8,ng) )
       ghost_w_img_cc(1,ng) = ii;  ghost_w_img_cc(2,ng) = jj;  ghost_w_img_cc(3,ng) = kk
    End Do

  End Subroutine build_ghost_list_w

  !> Apply ghost-cell IBM every RK sub-step in place of volume-penalisation
  Subroutine apply_ghost_cell_ibm(U_,V_,W_)

    Real(Int64), Dimension(nx, nyg,nzg), Intent(InOut) :: U_
    Real(Int64), Dimension(nxg, ny,nzg), Intent(InOut) :: V_
    Real(Int64), Dimension(nxg,nyg, nz), Intent(InOut) :: W_

    Integer(Int32) :: n, i, j, k
    Real   (Int64) :: U_I

    ! Zero fully-solid faces using phi-based averages (face-centred; Umask_cc is cell-centred).
    ! Ghost cells are also zeroed here and corrected below from fluid image points.
    !$acc parallel loop collapse(3) present(phi,U_)
    Do k = 2, nzg-1
       Do j = 2, nyg-1
          Do i = 2, nx-1
             If ( 0.5d0*(phi(i,j,k)+phi(i+1,j,k)) < 0d0 ) U_(i,j,k) = 0d0
          End Do
       End Do
    End Do
    !$acc end parallel loop
    !$acc parallel loop collapse(3) present(phi,V_)
    Do k = 2, nzg-1
       Do j = 2, ny-1
          Do i = 2, nxg-1
             If ( 0.5d0*(phi(i,j,k)+phi(i,j+1,k)) < 0d0 ) V_(i,j,k) = 0d0
          End Do
       End Do
    End Do
    !$acc end parallel loop
    !$acc parallel loop collapse(3) present(phi,W_)
    Do k = 2, nz-1
       Do j = 2, nyg-1
          Do i = 2, nxg-1
             If ( 0.5d0*(phi(i,j,k)+phi(i,j,k+1)) < 0d0 ) W_(i,j,k) = 0d0
          End Do
       End Do
    End Do
    !$acc end parallel loop

    !--- U ghost cells: U_G = 2*U_wall - U_I  (image at 2*dGB, so r = 0.5 always) ---
    !$acc parallel loop present(U_,ghost_u_idx,ghost_u_wgt,ghost_u_img) private(U_I)
    Do n = 1, n_ghost_u
       i = ghost_u_idx(1,n);  j = ghost_u_idx(2,n);  k = ghost_u_idx(3,n)
       U_I = trilinear_interp_u(U_, ghost_u_wgt(1:8,n), ghost_u_img(:,n))
       U_(i,j,k) = 2d0*U_wall - U_I
    End Do
    !$acc end parallel loop

    !--- V ghost cells ---
    !$acc parallel loop present(V_,ghost_v_idx,ghost_v_wgt,ghost_v_img) private(U_I)
    Do n = 1, n_ghost_v
       i = ghost_v_idx(1,n);  j = ghost_v_idx(2,n);  k = ghost_v_idx(3,n)
       U_I = trilinear_interp_v(V_, ghost_v_wgt(1:8,n), ghost_v_img(:,n))
       V_(i,j,k) = 2d0*V_wall - U_I
    End Do
    !$acc end parallel loop

    !--- W ghost cells ---
    !$acc parallel loop present(W_,ghost_w_idx,ghost_w_wgt,ghost_w_img) private(U_I)
    Do n = 1, n_ghost_w
       i = ghost_w_idx(1,n);  j = ghost_w_idx(2,n);  k = ghost_w_idx(3,n)
       U_I = trilinear_interp_w(W_, ghost_w_wgt(1:8,n), ghost_w_img(:,n))
       W_(i,j,k) = 2d0*W_wall - U_I
    End Do
    !$acc end parallel loop

  End Subroutine apply_ghost_cell_ibm

  !> Ghost-cell thermal condition at the immersed boundary, per solid ID: ibm_T_bc_type(id) 0=adiabatic (T_ghost=T_image), 1=isothermal (mirror against ibm_T_wall(id))
  Subroutine apply_ghost_cell_ibm_scalar(T_)

    Real(Int64), Dimension(nxg,nyg,nzg), Intent(InOut) :: T_

    Integer(Int32) :: n, i, j, k, oid
    Real   (Int64) :: T_I

    !$acc parallel loop present(T_,ghost_cc_idx,ghost_cc_wgt_cc,ghost_cc_img_cc,ghost_cc_objid,ibm_T_bc_type,ibm_T_wall) private(T_I,oid)
    Do n = 1, n_ghost_cc
       i = ghost_cc_idx(1,n);  j = ghost_cc_idx(2,n);  k = ghost_cc_idx(3,n)
       oid = ghost_cc_objid(n)
       T_I = trilinear_interp_p(T_, ghost_cc_wgt_cc(1:8,n), ghost_cc_img_cc(:,n))
       If ( ibm_T_bc_type(oid) == 1 ) Then
          T_(i,j,k) = 2d0*ibm_T_wall(oid) - T_I   ! isothermal: Dirichlet mirror
       Else
          T_(i,j,k) = T_I                         ! adiabatic: zero-gradient
       End If
    End Do
    !$acc end parallel loop

  End Subroutine apply_ghost_cell_ibm_scalar

  !                     Helper routines

  !  Zero the per-step IBM force accumulators.
  !  Call once per time step before the first RK IBM application.
  Subroutine zero_ibm_stage_accumulators
    ibm_Fx_acc = 0d0
    ibm_Fy_acc = 0d0
    ibm_Fz_acc = 0d0
  End Subroutine zero_ibm_stage_accumulators

  !> Accumulate IBM momentum exchange from one IBM application (pre/post state)
  Subroutine accumulate_ibm_stage_impulse(U_pre_, V_pre_, W_pre_, U_, V_, W_)

    Real(Int64), Dimension(nx, nyg,nzg), Intent(In) :: U_pre_, U_
    Real(Int64), Dimension(nxg, ny,nzg), Intent(In) :: V_pre_, V_
    Real(Int64), Dimension(nxg,nyg, nz), Intent(In) :: W_pre_, W_

    Integer(Int32) :: n, i, j, k
    Real   (Int64) :: dV

    Do n = 1, n_ghost_u
       i = ghost_u_idx(1,n);  j = ghost_u_idx(2,n);  k = ghost_u_idx(3,n)
       dV = (x(i+1)-x(i-1))*0.5d0 * (y(j)-y(j-1)) * (z(k)-z(k-1))
       ibm_Fx_acc = ibm_Fx_acc - (U_(i,j,k) - U_pre_(i,j,k)) * dV
    End Do

    Do n = 1, n_ghost_v
       i = ghost_v_idx(1,n);  j = ghost_v_idx(2,n);  k = ghost_v_idx(3,n)
       dV = (xg(i+1)-xg(i)) * (y(j+1)-y(j-1))*0.5d0 * (z(k)-z(k-1))
       ibm_Fy_acc = ibm_Fy_acc - (V_(i,j,k) - V_pre_(i,j,k)) * dV
    End Do

    Do n = 1, n_ghost_w
       i = ghost_w_idx(1,n);  j = ghost_w_idx(2,n);  k = ghost_w_idx(3,n)
       dV = (xg(i+1)-xg(i)) * (y(j)-y(j-1)) * (z(k+1)-z(k-1))*0.5d0
       ibm_Fz_acc = ibm_Fz_acc - (W_(i,j,k) - W_pre_(i,j,k)) * dV
    End Do

  End Subroutine accumulate_ibm_stage_impulse

  !> Build cell-centre ghost list for rigorous pressure integration
  Subroutine build_ghost_list_cc

    Integer(Int32) :: i, j, k, ng, nd
    Integer(Int32) :: ii, jj, kk
    Integer(Int32) :: iu, ju, ku, iv, jv, kv, iw, jw, kw
    Real   (Int64) :: nx_, ny_, nz_, dGB, xI, yI, zI

    ! Pass 1: count qualifying cells
    ng = 0;  nd = 0
    Do k = 2, nzg-1
       Do j = 2, nyg-1
          Do i = 2, nxg-1
             If ( phi(i,j,k) < 0d0 ) Then
                If ( phi(i-1,j,k) >= 0d0 .Or. phi(i+1,j,k) >= 0d0 .Or. &
                     phi(i,j-1,k) >= 0d0 .Or. phi(i,j+1,k) >= 0d0 .Or. &
                     phi(i,j,k-1) >= 0d0 .Or. phi(i,j,k+1) >= 0d0 ) Then
                   Call compute_normal_at_cc(i,j,k, nx_,ny_,nz_)
                   dGB = Abs(phi(i,j,k))
                   xI = xg(i) + 2d0*dGB*nx_
                   yI = yg(j) + 2d0*dGB*ny_
                   zI = zg(k) + 2d0*dGB*nz_
                   Call find_stencil_centre(xI, yI, zI, ii, jj, kk)
                   If ( ii>=2 .And. ii<=nxg-1 .And. jj>=2 .And. jj<=nyg-1 .And. &
                        kk>=2 .And. kk<=nzg-1 ) Then
                      If ( Umask_cc(ii,jj,kk) > 0.5d0 ) Then
                         ng = ng + 1
                      Else
                         nd = nd + 1   ! image in solid — concave corner
                      End If
                   Else
                      nd = nd + 1      ! image outside domain bounds
                   End If
                End If
             End If
          End Do
       End Do
    End Do
    If ( nd > 0 .And. myid == 0 ) &
       Write(*,'(A,I0,A)') '[IBM] build_ghost_list_cc: ', nd, &
          ' ghost cell(s) dropped (image in solid or outside domain)'

    n_ghost_cc = ng
    Allocate ( ghost_cc_idx   (3, n_ghost_cc) )
    Allocate ( ghost_cc_nrm   (3, n_ghost_cc) )
    Allocate ( ghost_cc_dGB   (   n_ghost_cc) )
    Allocate ( ghost_cc_img_cc(3, n_ghost_cc) )
    Allocate ( ghost_cc_wgt_cc(8, n_ghost_cc) )
    Allocate ( ghost_cc_objid (   n_ghost_cc) )
    ! Staggered image-point stencils for surface-sampling viscous traction.
    Allocate ( ghost_cc_img_u (3, n_ghost_cc) )
    Allocate ( ghost_cc_img_v (3, n_ghost_cc) )
    Allocate ( ghost_cc_img_w (3, n_ghost_cc) )
    Allocate ( ghost_cc_wgt_u (8, n_ghost_cc) )
    Allocate ( ghost_cc_wgt_v (8, n_ghost_cc) )
    Allocate ( ghost_cc_wgt_w (8, n_ghost_cc) )

    ! Pass 2: identical traversal, fill arrays
    ng = 0
    Do k = 2, nzg-1
       Do j = 2, nyg-1
          Do i = 2, nxg-1
             If ( phi(i,j,k) < 0d0 ) Then
                If ( phi(i-1,j,k) >= 0d0 .Or. phi(i+1,j,k) >= 0d0 .Or. &
                     phi(i,j-1,k) >= 0d0 .Or. phi(i,j+1,k) >= 0d0 .Or. &
                     phi(i,j,k-1) >= 0d0 .Or. phi(i,j,k+1) >= 0d0 ) Then
                   Call compute_normal_at_cc(i,j,k, nx_,ny_,nz_)
                   dGB = Abs(phi(i,j,k))
                   xI = xg(i) + 2d0*dGB*nx_
                   yI = yg(j) + 2d0*dGB*ny_
                   zI = zg(k) + 2d0*dGB*nz_
                   Call find_stencil_centre(xI, yI, zI, ii, jj, kk)
                   If ( ii>=2 .And. ii<=nxg-1 .And. jj>=2 .And. jj<=nyg-1 .And. &
                        kk>=2 .And. kk<=nzg-1 ) Then
                      If ( Umask_cc(ii,jj,kk) > 0.5d0 ) Then
                         ng = ng + 1
                         ghost_cc_idx(1,ng) = i
                         ghost_cc_idx(2,ng) = j
                         ghost_cc_idx(3,ng) = k
                         ghost_cc_objid(ng) = Min(Max(Nint(ibm_obj_id(i,j,k)), 0), max_ibm_objects)
                         ghost_cc_nrm(1,ng) = nx_
                         ghost_cc_nrm(2,ng) = ny_
                         ghost_cc_nrm(3,ng) = nz_
                         ghost_cc_dGB(ng)   = dGB
                         ghost_cc_img_cc(1,ng) = ii
                         ghost_cc_img_cc(2,ng) = jj
                         ghost_cc_img_cc(3,ng) = kk
                         Call trilinear_weights(xI, yI, zI, ii, jj, kk, &
                                                ghost_cc_wgt_cc(1:8,ng))
                         ! Staggered velocity stencils at the same image point I,
                         ! precomputed for sample_ibm_surface viscous traction.
                         Call find_stencil_u(xI, yI, zI, iu, ju, ku)
                         Call trilinear_weights_u(xI, yI, zI, iu, ju, ku, &
                                                  ghost_cc_wgt_u(1:8,ng))
                         ghost_cc_img_u(1,ng) = iu
                         ghost_cc_img_u(2,ng) = ju
                         ghost_cc_img_u(3,ng) = ku
                         Call find_stencil_v(xI, yI, zI, iv, jv, kv)
                         Call trilinear_weights_v(xI, yI, zI, iv, jv, kv, &
                                                  ghost_cc_wgt_v(1:8,ng))
                         ghost_cc_img_v(1,ng) = iv
                         ghost_cc_img_v(2,ng) = jv
                         ghost_cc_img_v(3,ng) = kv
                         Call find_stencil_w(xI, yI, zI, iw, jw, kw)
                         Call trilinear_weights_w(xI, yI, zI, iw, jw, kw, &
                                                  ghost_cc_wgt_w(1:8,ng))
                         ghost_cc_img_w(1,ng) = iw
                         ghost_cc_img_w(2,ng) = jw
                         ghost_cc_img_w(3,ng) = kw
                      End If
                   End If
                End If
             End If
          End Do
       End Do
    End Do

  End Subroutine build_ghost_list_cc

  !  Wall-normal at a cell centre (i,j,k): central-difference
  !  gradient of phi on the cell-centre grid, normalised.
  Subroutine compute_normal_at_cc(i,j,k, nx_,ny_,nz_)
    Integer(Int32), Intent(In)  :: i, j, k
    Real   (Int64), Intent(Out) :: nx_, ny_, nz_
    Real   (Int64)              :: nmag, h_up, h_dn
    nx_  = (phi(i+1,j,k) - phi(i-1,j,k)) / (xg(i+1) - xg(i-1))
    ! Non-uniform central difference in y (2nd-order on stretched meshes)
    h_up = yg(j+1) - yg(j);  h_dn = yg(j) - yg(j-1)
    ny_  = ( h_dn**2*phi(i,j+1,k) + (h_up**2-h_dn**2)*phi(i,j,k) - h_up**2*phi(i,j-1,k) ) &
           / ( h_up * h_dn * (h_up + h_dn) )
    nz_  = (phi(i,j,k+1) - phi(i,j,k-1)) / (zg(k+1) - zg(k-1))
    nmag = Sqrt(nx_**2 + ny_**2 + nz_**2)
    If (nmag > 1d-14) Then
       nx_ = nx_/nmag;  ny_ = ny_/nmag;  nz_ = nz_/nmag
    Else
       nx_ = 0d0;  ny_ = 1d0;  nz_ = 0d0
    End If
  End Subroutine compute_normal_at_cc

  !  Wall-normal at a U-face (i,j,k): gradient of phi at surrounding
  !  cell centres, normalised.
  Subroutine compute_normal_at_face_u(i,j,k, nx_,ny_,nz_)
    Integer(Int32), Intent(In)  :: i, j, k
    Real   (Int64), Intent(Out) :: nx_, ny_, nz_
    Real   (Int64)              :: nmag, h_up, h_dn

    nx_  = ( phi(i+1,j,k) - phi(i,j,k) ) / ( xg(i+1) - xg(i) )
    ! Non-uniform central difference in y (2nd-order on stretched meshes)
    h_up = yg(j+1) - yg(j);  h_dn = yg(j) - yg(j-1)
    ny_  = ( h_dn**2*phi(i,j+1,k) + (h_up**2-h_dn**2)*phi(i,j,k) - h_up**2*phi(i,j-1,k) ) &
           / ( h_up * h_dn * (h_up + h_dn) )
    nz_  = ( phi(i,j,k+1) - phi(i,j,k-1) ) / ( zg(k+1) - zg(k-1) )
    nmag = Sqrt(nx_**2 + ny_**2 + nz_**2)
    If (nmag > 1d-14) Then
       nx_ = nx_/nmag;  ny_ = ny_/nmag;  nz_ = nz_/nmag
    Else
       nx_ = 0d0;  ny_ = 1d0;  nz_ = 0d0   ! fallback: wall-normal in y
    End If
  End Subroutine compute_normal_at_face_u

  Subroutine compute_normal_at_face_v(i,j,k, nx_,ny_,nz_)
    Integer(Int32), Intent(In)  :: i, j, k
    Real   (Int64), Intent(Out) :: nx_, ny_, nz_
    Real   (Int64)              :: nmag
    nx_ = ( phi(i+1,j,k) - phi(i-1,j,k) ) / ( xg(i+1) - xg(i-1) )
    ny_ = ( phi(i,j+1,k) - phi(i,j,k) )   / ( yg(j+1) - yg(j) )
    nz_ = ( phi(i,j,k+1) - phi(i,j,k-1) ) / ( zg(k+1) - zg(k-1) )
    nmag = Sqrt(nx_**2+ny_**2+nz_**2)
    If (nmag>1d-14) Then; nx_=nx_/nmag; ny_=ny_/nmag; nz_=nz_/nmag
    Else; nx_=0d0; ny_=1d0; nz_=0d0; End If
  End Subroutine compute_normal_at_face_v

  Subroutine compute_normal_at_face_w(i,j,k, nx_,ny_,nz_)
    Integer(Int32), Intent(In)  :: i, j, k
    Real   (Int64), Intent(Out) :: nx_, ny_, nz_
    Real   (Int64)              :: nmag, h_up, h_dn
    nx_  = ( phi(i+1,j,k) - phi(i-1,j,k) ) / ( xg(i+1) - xg(i-1) )
    ! Non-uniform central difference in y (2nd-order on stretched meshes)
    h_up = yg(j+1) - yg(j);  h_dn = yg(j) - yg(j-1)
    ny_  = ( h_dn**2*phi(i,j+1,k) + (h_up**2-h_dn**2)*phi(i,j,k) - h_up**2*phi(i,j-1,k) ) &
           / ( h_up * h_dn * (h_up + h_dn) )
    nz_  = ( phi(i,j,k+1) - phi(i,j,k) )   / ( zg(k+1) - zg(k) )
    nmag = Sqrt(nx_**2+ny_**2+nz_**2)
    If (nmag>1d-14) Then; nx_=nx_/nmag; ny_=ny_/nmag; nz_=nz_/nmag
    Else; nx_=0d0; ny_=1d0; nz_=0d0; End If
  End Subroutine compute_normal_at_face_w

  !  Find lower-left cell-centre index (ii,jj,kk) such that
  !  xg(ii) <= xp < xg(ii+1) etc. — trilinear stencil anchor.
  Subroutine find_stencil_centre(xp, yp, zp, ii, jj, kk)
    Real   (Int64), Intent(In)  :: xp, yp, zp
    Integer(Int32), Intent(Out) :: ii, jj, kk
    Integer(Int32) :: n

    ii = 2
    Do n = 2, nxg-2
       If ( xg(n) <= xp ) ii = n
    End Do

    jj = 2
    Do n = 2, nyg-2
       If ( yg(n) <= yp ) jj = n
    End Do

    kk = 2
    Do n = 2, nzg-2
       If ( zg(n) <= zp ) kk = n
    End Do
  End Subroutine find_stencil_centre

  ! Staggered stencil finders, one per velocity component: U uses x-faces x cell-centres, V cell-centres x y-faces, W cell-centres x z-faces
  Subroutine find_stencil_u(xp, yp, zp, ii, jj, kk)
    Real   (Int64), Intent(In)  :: xp, yp, zp
    Integer(Int32), Intent(Out) :: ii, jj, kk
    Integer(Int32) :: n
    ii = 2
    Do n = 2, nx-1
       If ( x(n) <= xp ) ii = n
    End Do
    jj = 2
    Do n = 2, nyg-2
       If ( yg(n) <= yp ) jj = n
    End Do
    kk = 2
    Do n = 2, nzg-2
       If ( zg(n) <= zp ) kk = n
    End Do
  End Subroutine find_stencil_u

  Subroutine find_stencil_v(xp, yp, zp, ii, jj, kk)
    Real   (Int64), Intent(In)  :: xp, yp, zp
    Integer(Int32), Intent(Out) :: ii, jj, kk
    Integer(Int32) :: n
    ii = 2
    Do n = 2, nxg-2
       If ( xg(n) <= xp ) ii = n
    End Do
    jj = 2
    Do n = 2, ny-1
       If ( y(n) <= yp ) jj = n
    End Do
    kk = 2
    Do n = 2, nzg-2
       If ( zg(n) <= zp ) kk = n
    End Do
  End Subroutine find_stencil_v

  Subroutine find_stencil_w(xp, yp, zp, ii, jj, kk)
    Real   (Int64), Intent(In)  :: xp, yp, zp
    Integer(Int32), Intent(Out) :: ii, jj, kk
    Integer(Int32) :: n
    ii = 2
    Do n = 2, nxg-2
       If ( xg(n) <= xp ) ii = n
    End Do
    jj = 2
    Do n = 2, nyg-2
       If ( yg(n) <= yp ) jj = n
    End Do
    kk = 2
    Do n = 2, nz-1
       If ( z(n) <= zp ) kk = n
    End Do
  End Subroutine find_stencil_w

  !> Compute 8 trilinear weights (unit-cube corner order) for point (xp,yp,zp) anchored at centre index (ii,jj,kk)
  Subroutine trilinear_weights(xp, yp, zp, ii, jj, kk, w)
    Real   (Int64), Intent(In)  :: xp, yp, zp
    Integer(Int32), Intent(In)  :: ii, jj, kk
    Real   (Int64), Intent(Out) :: w(8)
    Real   (Int64) :: tx, ty, tz

    tx = (xp - xg(ii)) / Max(xg(ii+1)-xg(ii), 1d-14)
    ty = (yp - yg(jj)) / Max(yg(jj+1)-yg(jj), 1d-14)
    tz = (zp - zg(kk)) / Max(zg(kk+1)-zg(kk), 1d-14)

    tx = Max(0d0, Min(1d0, tx))
    ty = Max(0d0, Min(1d0, ty))
    tz = Max(0d0, Min(1d0, tz))

    w(1) = (1d0-tx)*(1d0-ty)*(1d0-tz)
    w(2) =      tx *(1d0-ty)*(1d0-tz)
    w(3) = (1d0-tx)*     ty *(1d0-tz)
    w(4) =      tx *     ty *(1d0-tz)
    w(5) = (1d0-tx)*(1d0-ty)*     tz
    w(6) =      tx *(1d0-ty)*     tz
    w(7) = (1d0-tx)*     ty *     tz
    w(8) =      tx *     ty *     tz
  End Subroutine trilinear_weights

  ! Staggered trilinear weight routines, one per velocity component's grid, matching find_stencil_u/v/w
  Subroutine trilinear_weights_u(xp, yp, zp, ii, jj, kk, w)
    Real   (Int64), Intent(In)  :: xp, yp, zp
    Integer(Int32), Intent(In)  :: ii, jj, kk
    Real   (Int64), Intent(Out) :: w(8)
    Real   (Int64) :: tx, ty, tz
    tx = (xp - x(ii))  / Max(x(ii+1)  - x(ii),  1d-14)
    ty = (yp - yg(jj)) / Max(yg(jj+1) - yg(jj), 1d-14)
    tz = (zp - zg(kk)) / Max(zg(kk+1) - zg(kk), 1d-14)
    tx = Max(0d0, Min(1d0, tx));  ty = Max(0d0, Min(1d0, ty));  tz = Max(0d0, Min(1d0, tz))
    w(1) = (1d0-tx)*(1d0-ty)*(1d0-tz);  w(2) =      tx *(1d0-ty)*(1d0-tz)
    w(3) = (1d0-tx)*     ty *(1d0-tz);  w(4) =      tx *     ty *(1d0-tz)
    w(5) = (1d0-tx)*(1d0-ty)*     tz ;  w(6) =      tx *(1d0-ty)*     tz
    w(7) = (1d0-tx)*     ty *     tz ;  w(8) =      tx *     ty *     tz
  End Subroutine trilinear_weights_u

  Subroutine trilinear_weights_v(xp, yp, zp, ii, jj, kk, w)
    Real   (Int64), Intent(In)  :: xp, yp, zp
    Integer(Int32), Intent(In)  :: ii, jj, kk
    Real   (Int64), Intent(Out) :: w(8)
    Real   (Int64) :: tx, ty, tz
    tx = (xp - xg(ii)) / Max(xg(ii+1) - xg(ii), 1d-14)
    ty = (yp - y(jj))  / Max(y(jj+1)  - y(jj),  1d-14)
    tz = (zp - zg(kk)) / Max(zg(kk+1) - zg(kk), 1d-14)
    tx = Max(0d0, Min(1d0, tx));  ty = Max(0d0, Min(1d0, ty));  tz = Max(0d0, Min(1d0, tz))
    w(1) = (1d0-tx)*(1d0-ty)*(1d0-tz);  w(2) =      tx *(1d0-ty)*(1d0-tz)
    w(3) = (1d0-tx)*     ty *(1d0-tz);  w(4) =      tx *     ty *(1d0-tz)
    w(5) = (1d0-tx)*(1d0-ty)*     tz ;  w(6) =      tx *(1d0-ty)*     tz
    w(7) = (1d0-tx)*     ty *     tz ;  w(8) =      tx *     ty *     tz
  End Subroutine trilinear_weights_v

  Subroutine trilinear_weights_w(xp, yp, zp, ii, jj, kk, w)
    Real   (Int64), Intent(In)  :: xp, yp, zp
    Integer(Int32), Intent(In)  :: ii, jj, kk
    Real   (Int64), Intent(Out) :: w(8)
    Real   (Int64) :: tx, ty, tz
    tx = (xp - xg(ii)) / Max(xg(ii+1) - xg(ii), 1d-14)
    ty = (yp - yg(jj)) / Max(yg(jj+1) - yg(jj), 1d-14)
    tz = (zp - z(kk))  / Max(z(kk+1)  - z(kk),  1d-14)
    tx = Max(0d0, Min(1d0, tx));  ty = Max(0d0, Min(1d0, ty));  tz = Max(0d0, Min(1d0, tz))
    w(1) = (1d0-tx)*(1d0-ty)*(1d0-tz);  w(2) =      tx *(1d0-ty)*(1d0-tz)
    w(3) = (1d0-tx)*     ty *(1d0-tz);  w(4) =      tx *     ty *(1d0-tz)
    w(5) = (1d0-tx)*(1d0-ty)*     tz ;  w(6) =      tx *(1d0-ty)*     tz
    w(7) = (1d0-tx)*     ty *     tz ;  w(8) =      tx *     ty *     tz
  End Subroutine trilinear_weights_w

  !> Trilinear interpolation at image point using stencil anchor (img) and weights (w) from setup_ibm, on the component's staggered grid
  Real(Int64) Function trilinear_interp_u(U_, w, img)
    !$acc routine seq
    Real   (Int64), Dimension(nx,nyg,nzg), Intent(In) :: U_
    Real   (Int64), Dimension(8),          Intent(In) :: w
    Integer(Int32), Dimension(3),          Intent(In) :: img   ! image stencil anchor on U grid
    Integer(Int32) :: i0, j0, k0
    i0 = img(1);  j0 = img(2);  k0 = img(3)
    trilinear_interp_u = &
         w(1)*U_(i0  ,j0  ,k0  ) + w(2)*U_(i0+1,j0  ,k0  ) + &
         w(3)*U_(i0  ,j0+1,k0  ) + w(4)*U_(i0+1,j0+1,k0  ) + &
         w(5)*U_(i0  ,j0  ,k0+1) + w(6)*U_(i0+1,j0  ,k0+1) + &
         w(7)*U_(i0  ,j0+1,k0+1) + w(8)*U_(i0+1,j0+1,k0+1)
  End Function trilinear_interp_u

  Real(Int64) Function trilinear_interp_v(V_, w, img)
    !$acc routine seq
    Real   (Int64), Dimension(nxg,ny,nzg), Intent(In) :: V_
    Real   (Int64), Dimension(8),          Intent(In) :: w
    Integer(Int32), Dimension(3),          Intent(In) :: img   ! image stencil anchor on V grid
    Integer(Int32) :: i0, j0, k0
    i0 = img(1);  j0 = img(2);  k0 = img(3)
    trilinear_interp_v = &
         w(1)*V_(i0  ,j0  ,k0  ) + w(2)*V_(i0+1,j0  ,k0  ) + &
         w(3)*V_(i0  ,j0+1,k0  ) + w(4)*V_(i0+1,j0+1,k0  ) + &
         w(5)*V_(i0  ,j0  ,k0+1) + w(6)*V_(i0+1,j0  ,k0+1) + &
         w(7)*V_(i0  ,j0+1,k0+1) + w(8)*V_(i0+1,j0+1,k0+1)
  End Function trilinear_interp_v

  Real(Int64) Function trilinear_interp_w(W_, w, img)
    !$acc routine seq
    Real   (Int64), Dimension(nxg,nyg,nz), Intent(In) :: W_
    Real   (Int64), Dimension(8),          Intent(In) :: w
    Integer(Int32), Dimension(3),          Intent(In) :: img   ! image stencil anchor on W grid
    Integer(Int32) :: i0, j0, k0
    i0 = img(1);  j0 = img(2);  k0 = img(3)
    trilinear_interp_w = &
         w(1)*W_(i0  ,j0  ,k0  ) + w(2)*W_(i0+1,j0  ,k0  ) + &
         w(3)*W_(i0  ,j0+1,k0  ) + w(4)*W_(i0+1,j0+1,k0  ) + &
         w(5)*W_(i0  ,j0  ,k0+1) + w(6)*W_(i0+1,j0  ,k0+1) + &
         w(7)*W_(i0  ,j0+1,k0+1) + w(8)*W_(i0+1,j0+1,k0+1)
  End Function trilinear_interp_w

  !  Trilinear interpolation of the cell-centred pressure field P
  !  to an arbitrary point, using the centre-grid stencil.
  Real(Int64) Function trilinear_interp_p(P_, w, img)
    !$acc routine seq
    Real   (Int64), Dimension(nxg,nyg,nzg), Intent(In) :: P_
    Real   (Int64), Dimension(8),           Intent(In) :: w
    Integer(Int32), Dimension(3),           Intent(In) :: img
    Integer(Int32) :: i0, j0, k0
    i0 = img(1);  j0 = img(2);  k0 = img(3)
    trilinear_interp_p = &
         w(1)*P_(i0  ,j0  ,k0  ) + w(2)*P_(i0+1,j0  ,k0  ) + &
         w(3)*P_(i0  ,j0+1,k0  ) + w(4)*P_(i0+1,j0+1,k0  ) + &
         w(5)*P_(i0  ,j0  ,k0+1) + w(6)*P_(i0+1,j0  ,k0+1) + &
         w(7)*P_(i0  ,j0+1,k0+1) + w(8)*P_(i0+1,j0+1,k0+1)
  End Function trilinear_interp_p

  !> Compute IBM forces every nsampling steps via Method 1 (momentum exchange) and Method 2 (surface integral)
  Subroutine compute_ibm_forces(U_, V_, W_,                           &
       Fx_ibm,  Fy_ibm,  Fz_ibm,                                      &
       Fx_pres, Fy_pres, Fz_pres,                                     &
       Fx_visc, Fy_visc, Fz_visc)

    Real(Int64), Dimension(nx, nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg, ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg, nz), Intent(In) :: W_
    Real(Int64), Intent(Out) :: Fx_ibm, Fy_ibm, Fz_ibm
    Real(Int64), Intent(Out) :: Fx_pres, Fy_pres, Fz_pres
    Real(Int64), Intent(Out) :: Fx_visc, Fy_visc, Fz_visc

    Integer(Int32) :: n, i, j, k
    Real   (Int64) :: dV, dA, U_I, p_I, nu_t_B, dGB
    Real   (Int64) :: lFx_pres, lFy_pres, lFz_pres
    Real   (Int64) :: lFx_visc, lFy_visc, lFz_visc

    lFx_pres=0d0; lFy_pres=0d0; lFz_pres=0d0
    lFx_visc=0d0; lFy_visc=0d0; lFz_visc=0d0

    !--- Method 1: use accumulated impulse from all 3 RK stages ---
    ! ibm_F?_acc holds the rank-local sum of (U_after_IBM - U_before_IBM)*dV; dividing by dt gives force per unit time
    Call MPI_Allreduce(ibm_Fx_acc/dt, Fx_ibm, 1, MPI_real8, MPI_SUM, MPI_COMM_WORLD, ierr)
    Call MPI_Allreduce(ibm_Fy_acc/dt, Fy_ibm, 1, MPI_real8, MPI_SUM, MPI_COMM_WORLD, ierr)
    Call MPI_Allreduce(ibm_Fz_acc/dt, Fz_ibm, 1, MPI_real8, MPI_SUM, MPI_COMM_WORLD, ierr)

    !--- Method 2 pressure: cell-centre ghost list (one sample per interface cell) ---
    Do n = 1, n_ghost_cc
       i = ghost_cc_idx(1,n);  j = ghost_cc_idx(2,n);  k = ghost_cc_idx(3,n)
       dGB = ghost_cc_dGB(n)
       dV  = (xg(i+1)-xg(i-1))*0.5d0 * (yg(j+1)-yg(j-1))*0.5d0 * (zg(k+1)-zg(k-1))*0.5d0
       dA  = dV / Max(dGB, 1d-14)
       p_I = trilinear_interp_p(P, ghost_cc_wgt_cc(1:8,n), ghost_cc_img_cc(:,n))
       lFx_pres = lFx_pres - p_I * ghost_cc_nrm(1,n) * dA
       lFy_pres = lFy_pres - p_I * ghost_cc_nrm(2,n) * dA
       lFz_pres = lFz_pres - p_I * ghost_cc_nrm(3,n) * dA
    End Do

    !--- Method 2 viscous: staggered ghost lists, one component each ---
    ! dA is the perpendicular face area; using dV/dGB instead would give a 1/dGB^2 singularity since the gradient already carries 1/dGB
    Do n = 1, n_ghost_u
       i = ghost_u_idx(1,n);  j = ghost_u_idx(2,n);  k = ghost_u_idx(3,n)
       dGB = ghost_u_dGB(n)
       dA  = (y(j)-y(j-1)) * (z(k)-z(k-1))             ! y-z face area
       U_I = trilinear_interp_u(U_, ghost_u_wgt(1:8,n), ghost_u_img(:,n))
       nu_t_B = 0.5d0*(nu_t(i,j,k) + nu_t(Min(i+1,nxg),j,k))
       lFx_visc = lFx_visc + (nu + nu_t_B) * (U_I - U_wall) / Max(dGB, 1d-14) * dA
    End Do

    Do n = 1, n_ghost_v
       i = ghost_v_idx(1,n);  j = ghost_v_idx(2,n);  k = ghost_v_idx(3,n)
       dGB = ghost_v_dGB(n)
       dA  = (xg(i+1)-xg(i)) * (z(k)-z(k-1))           ! x-z face area
       U_I = trilinear_interp_v(V_, ghost_v_wgt(1:8,n), ghost_v_img(:,n))
       nu_t_B = 0.5d0*(nu_t(i,j,k) + nu_t(i,Min(j+1,nyg),k))
       lFy_visc = lFy_visc + (nu + nu_t_B) * (U_I - V_wall) / Max(dGB, 1d-14) * dA
    End Do

    Do n = 1, n_ghost_w
       i = ghost_w_idx(1,n);  j = ghost_w_idx(2,n);  k = ghost_w_idx(3,n)
       dGB = ghost_w_dGB(n)
       dA  = (xg(i+1)-xg(i)) * (y(j)-y(j-1))           ! x-y face area
       U_I = trilinear_interp_w(W_, ghost_w_wgt(1:8,n), ghost_w_img(:,n))
       nu_t_B = 0.5d0*(nu_t(i,j,k) + nu_t(i,j,Min(k+1,nzg)))
       lFz_visc = lFz_visc + (nu + nu_t_B) * (U_I - W_wall) / Max(dGB, 1d-14) * dA
    End Do

    !--- Global reduction for Method 2 ---
    Call MPI_Allreduce(lFx_pres, Fx_pres, 1, MPI_real8, MPI_SUM, MPI_COMM_WORLD, ierr)
    Call MPI_Allreduce(lFy_pres, Fy_pres, 1, MPI_real8, MPI_SUM, MPI_COMM_WORLD, ierr)
    Call MPI_Allreduce(lFz_pres, Fz_pres, 1, MPI_real8, MPI_SUM, MPI_COMM_WORLD, ierr)
    Call MPI_Allreduce(lFx_visc, Fx_visc, 1, MPI_real8, MPI_SUM, MPI_COMM_WORLD, ierr)
    Call MPI_Allreduce(lFy_visc, Fy_visc, 1, MPI_real8, MPI_SUM, MPI_COMM_WORLD, ierr)
    Call MPI_Allreduce(lFz_visc, Fz_visc, 1, MPI_real8, MPI_SUM, MPI_COMM_WORLD, ierr)

  End Subroutine compute_ibm_forces

  !> Export the immersed-surface (phi=0) sample as a ParaView point cloud
  Subroutine sample_ibm_surface(U_, V_, W_)

    Real(Int64), Dimension(nx, nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg, ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg, nz), Intent(In) :: W_

    Integer(Int32), Parameter :: NF = 13     ! fields per point (see layout below)
    Integer(Int32) :: n, i, j, k, nl, ntot, iproc, off, funit
    Real   (Int64) :: dV, dA, dGB, invd, pB, nutB, uI, vI, wI
    Real   (Int64), Allocatable :: lbuf(:,:), gbuf(:,:)
    Integer(Int32), Allocatable :: counts(:), recvc(:), displs(:)
    Character(256) :: fname
    Character(16)  :: ext
    Logical        :: dir_exists

    nl = n_ghost_cc

    ! Fill the rank-local per-point buffer; row layout (NF=13): 1-3 x,y,z 4-6 nx,ny,nz 7 p 8-10 fp_x,fp_y,fp_z 11-13 fv_x,fv_y,fv_z
    Allocate ( lbuf(NF, Max(nl,1)) )
    Do n = 1, nl
       i = ghost_cc_idx(1,n);  j = ghost_cc_idx(2,n);  k = ghost_cc_idx(3,n)
       dGB  = ghost_cc_dGB(n)
       invd = 1d0 / Max(dGB, 1d-14)
       dV   = (xg(i+1)-xg(i-1))*0.5d0 * (yg(j+1)-yg(j-1))*0.5d0 * (zg(k+1)-zg(k-1))*0.5d0
       dA   = dV * invd

       ! interpolate pressure, turbulent viscosity and velocity at image point I
       pB   = trilinear_interp_p(P,    ghost_cc_wgt_cc(1:8,n), ghost_cc_img_cc(:,n))
       nutB = trilinear_interp_p(nu_t, ghost_cc_wgt_cc(1:8,n), ghost_cc_img_cc(:,n))
       uI   = trilinear_interp_u(U_,   ghost_cc_wgt_u (1:8,n), ghost_cc_img_u (:,n))
       vI   = trilinear_interp_v(V_,   ghost_cc_wgt_v (1:8,n), ghost_cc_img_v (:,n))
       wI   = trilinear_interp_w(W_,   ghost_cc_wgt_w (1:8,n), ghost_cc_img_w (:,n))

       ! boundary-point coordinates B = cell-centre + dGB*nrm
       lbuf(1,n) = xg(i) + dGB*ghost_cc_nrm(1,n)
       lbuf(2,n) = yg(j) + dGB*ghost_cc_nrm(2,n)
       lbuf(3,n) = zg(k) + dGB*ghost_cc_nrm(3,n)
       lbuf(4,n) = ghost_cc_nrm(1,n)
       lbuf(5,n) = ghost_cc_nrm(2,n)
       lbuf(6,n) = ghost_cc_nrm(3,n)
       lbuf(7,n) = pB
       ! pressure force per point  (-p n dA)
       lbuf(8,n)  = -pB * ghost_cc_nrm(1,n) * dA
       lbuf(9,n)  = -pB * ghost_cc_nrm(2,n) * dA
       lbuf(10,n) = -pB * ghost_cc_nrm(3,n) * dA
       ! viscous force per point  ((nu+nu_t) dU/dn dA)
       lbuf(11,n) = (nu + nutB) * (uI - U_wall) * invd * dA
       lbuf(12,n) = (nu + nutB) * (vI - V_wall) * invd * dA
       lbuf(13,n) = (nu + nutB) * (wI - W_wall) * invd * dA
    End Do

    ! ── Gather point counts, then the buffers, onto rank 0 ───────────
    Allocate ( counts(nprocs), recvc(nprocs), displs(nprocs) )
    Call MPI_Gather(nl, 1, MPI_integer, counts, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)

    ntot = 0
    If ( myid == 0 ) Then
       ntot = Sum(counts)
       off = 0
       Do iproc = 1, nprocs
          recvc(iproc)  = counts(iproc) * NF
          displs(iproc) = off
          off = off + recvc(iproc)
       End Do
    Else
       recvc  = 0
       displs = 0
    End If

    ! gbuf only needs to hold data on rank 0; a size-1 stub elsewhere keeps
    ! the MPI_Gatherv recvbuf argument valid (associated) on every rank.
    Allocate ( gbuf(NF, Max(ntot,1)) )
    Call MPI_Gatherv(lbuf, nl*NF, MPI_real8, &
                     gbuf, recvc, displs, MPI_real8, 0, MPI_COMM_WORLD, ierr)

    ! ── Rank 0 writes the binary point cloud ─────────────────────────
    If ( myid == 0 ) Then
       Inquire(file='ibm_surface/.', exist=dir_exists)
       If ( .Not. dir_exists ) &
          Call execute_command_line('mkdir -p ibm_surface', wait=.True.)

       Write(ext,'(I8.8)') istep + nstep_init
       fname = 'ibm_surface/surface.' // Trim(Adjustl(ext)) // '.bin'
       Open(newunit=funit, file=Trim(fname), access='stream', &
            form='unformatted', action='write', status='replace')
       Write(funit) ntot                    ! Int32  : number of surface points
       Write(funit) t                       ! float64: simulation time
       Write(funit) gbuf(1,1:ntot)          ! x
       Write(funit) gbuf(2,1:ntot)          ! y
       Write(funit) gbuf(3,1:ntot)          ! z
       Write(funit) gbuf(4,1:ntot)          ! nx
       Write(funit) gbuf(5,1:ntot)          ! ny
       Write(funit) gbuf(6,1:ntot)          ! nz
       Write(funit) gbuf(7,1:ntot)          ! pressure
       Write(funit) gbuf(8,1:ntot)          ! fp_x
       Write(funit) gbuf(9,1:ntot)          ! fp_y
       Write(funit) gbuf(10,1:ntot)         ! fp_z
       Write(funit) gbuf(11,1:ntot)         ! fv_x
       Write(funit) gbuf(12,1:ntot)         ! fv_y
       Write(funit) gbuf(13,1:ntot)         ! fv_z
       Close(funit)
    End If

    Deallocate ( lbuf, counts, recvc, displs, gbuf )

  End Subroutine sample_ibm_surface

End Module ibm

