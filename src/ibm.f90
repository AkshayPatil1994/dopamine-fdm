!> Ghost-Cell Immersed Boundary Method (Tseng & Ferziger 2003)
Module ibm

  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use mpi
  Use decomp, Only : z_halo_neighbors, x_halo_neighbors, x_periodic_partner, z_periodic_partner
  Use boundary_conditions, Only : exchange_velocity_halos
  Use halo_pad, Only : pad_axis, pad_field

  Implicit None

  ! Wall velocity (no-slip by default; set non-zero for moving walls)
  Real(Int64) :: U_wall = 0d0
  Real(Int64) :: V_wall = 0d0
  Real(Int64) :: W_wall = 0d0

  ! Minimum number of fluid-side cells to look ahead for image point
  Integer(Int32), Parameter :: n_image_layers = 2

  ! Ghost cells whose image stencil is not usable, counted while building the lists (1=U, 2=V, 3=W, 4=cell-centre):
  ! stencil clipping solid (concave corner) vs image outside this rank's local array (domain edge, or across a rank seam)
  Integer(Int32) :: ibm_drop_solid(4) = 0, ibm_drop_outside(4) = 0

  ! Extended halo for everything that samples the fields at image points. The solver arrays carry ONE ghost plane per side in x
  ! and z; the image point of a ghost cell lies up to ~2 cells from it, so near a rank seam its interpolation stencil (and the
  ! SDF around it) reaches past that plane. ibm_E extra planes (halo_pad.f90) make the stencil, and hence the ghost-cell
  ! boundary condition, independent of the rank layout. All setup-time lookups use the xge/xe/zge/ze axes and the phie/Umaske
  ! fields (bounds 1-ibm_E : n+ibm_E); the runtime kernels interpolate from Uext/Vext/Wext/Cext, refreshed before use.
  Integer(Int32), Save :: ibm_E = 2
  Logical,        Save :: ibm_lo_x = .False., ibm_hi_x = .False., ibm_lo_z = .False., ibm_hi_z = .False.
  Integer(Int32), Save :: ext_cx_lo, ext_cx_hi, ext_fx_lo, ext_fx_hi, ext_cz_lo, ext_cz_hi, ext_fz_lo, ext_fz_hi
  Logical,        Save :: stencil_out = .False.
  Real   (Int64), Allocatable, Dimension(:)     :: xge, xe, zge, ze
  Real   (Int64), Allocatable, Dimension(:,:,:) :: phie, Umaske
  Real   (Int64), Allocatable, Dimension(:,:,:) :: Uext, Vext, Wext, Cext
  Real   (Int64), Allocatable, Dimension(:)     :: ext_bs, ext_br

  ! Accumulators for Method 1 IBM force (summed over 6 IBM applications/step)
  Real(Int64) :: ibm_Fx_acc = 0d0
  Real(Int64) :: ibm_Fy_acc = 0d0
  Real(Int64) :: ibm_Fz_acc = 0d0

  ! Image-point values gathered from the pre-update field before any ghost cell is overwritten
  Real(Int64), Allocatable :: ghost_img_val(:)

Contains

  !  Read precomputed SDF and populate ghost-cell lists for U, V, W.
  !  Must be called after grid setup (xg/yg/zg must be initialised).
  Subroutine setup_ibm

    Integer(Int32) :: n_ghost_u_global, n_ghost_v_global, n_ghost_w_global

    If ( myid==0 ) Write(*,*) 'IBM: reading cell-centre SDF from ', Trim(ibm_sdf_file), '...'
    Call read_phi_from_sdf_file

    Call setup_ibm_ext

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

    ! Dropped ghost cells (no boundary condition applied there), summed over ranks; 'outside' > 0 for np>1 beyond what np=1 shows is a rank-seam loss
    Block
      Integer(Int32) :: ds(4), dout(4)
      Call MPI_Reduce(ibm_drop_solid,   ds,   4, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      Call MPI_Reduce(ibm_drop_outside, dout, 4, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      If ( myid==0 ) Then
         Write(*,'(A,4(I0,1X))') ' IBM: dropped ghost cells (GLOBAL) U,V,W,CC, image clips solid: ', ds
         Write(*,'(A,4(I0,1X))') ' IBM: dropped ghost cells (GLOBAL) U,V,W,CC, image outside local array: ', dout
      End If
    End Block

    Call trace_ghost_lists

    ! Global sum across all ranks for diagnostic
    Call MPI_Reduce(n_ghost_u, n_ghost_u_global, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    Call MPI_Reduce(n_ghost_v, n_ghost_v_global, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    Call MPI_Reduce(n_ghost_w, n_ghost_w_global, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    If ( myid==0 ) Then
       Write(*,*) 'IBM: ghost cells (GLOBAL) — U:', n_ghost_u_global, ' V:', n_ghost_v_global, ' W:', n_ghost_w_global
    End If

    Allocate ( ghost_img_val(Max(n_ghost_u, n_ghost_v, n_ghost_w, 1)) )
    !$acc enter data create(ghost_img_val)

    ! Device-resident IBM data: phi (read once above) and the ghost-cell lists (built once above), never modified afterward
    !$acc enter data copyin(phi)
    !$acc enter data copyin(ghost_u_idx,ghost_u_wgt,ghost_u_img,ghost_u_ref,ghost_u_nrm,ghost_u_yref,ghost_u_dGB,ghost_u_dGI,ghost_u_objid)
    !$acc enter data copyin(ghost_v_idx,ghost_v_wgt,ghost_v_img,ghost_v_ref,ghost_v_nrm,ghost_v_yref,ghost_v_dGB,ghost_v_dGI,ghost_v_objid)
    !$acc enter data copyin(ghost_w_idx,ghost_w_wgt,ghost_w_img,ghost_w_ref,ghost_w_nrm,ghost_w_yref,ghost_w_dGB,ghost_w_dGI,ghost_w_objid)
    ! Cell-centre ghost list: used both by the host-only Method-2 force diagnostic and (when boussinesq_flag>=1) by the device-resident apply_ghost_cell_ibm_scalar
    !$acc enter data copyin(ghost_cc_idx,ghost_cc_wgt_cc,ghost_cc_img_cc,ghost_cc_objid,ghost_cc_dGB,ghost_cc_dGI)

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

  !> Exchange phi x-ghost planes between x-neighbour ranks (interior-rank seams only; true
  !> domain x-edges are handled by the caller's Neumann assignment, same split as z above).
  !> Local buffers rather than the shared buffer_bcx* set (boundary_conditions module):
  !> smooth_ibm_corners runs from readSDF, before those buffers are allocated (initialization.f90).
  Subroutine exchange_phi_x_ghost_planes

    Real   (Int64) :: buf_s(nyg,nzg), buf_r(nyg,nzg)
    Integer(Int32) :: up, down

    Call x_halo_neighbors(up, down)

    ! Pass 1: send i=nxg-1 towards +x; receive from -x into i=1.
    buf_s = phi(nxg-1,:,:)
    Call MPI_Sendrecv( buf_s, nyg*nzg, MPI_real8, up,   0, &
                       buf_r, nyg*nzg, MPI_real8, down, 0, &
                       MPI_COMM_WORLD, istat, ierr )
    If ( down /= MPI_PROC_NULL ) phi(1,:,:) = buf_r

    ! Pass 2: send i=2 towards -x; receive from +x into i=nxg.
    buf_s = phi(2,:,:)
    Call MPI_Sendrecv( buf_s, nyg*nzg, MPI_real8, down, 0, &
                       buf_r, nyg*nzg, MPI_real8, up,   0, &
                       MPI_COMM_WORLD, istat, ierr )
    If ( up /= MPI_PROC_NULL ) phi(nxg,:,:) = buf_r

  End Subroutine exchange_phi_x_ghost_planes

  !> Periodic wrap of the cell-centred phi in x, host-only (setup runs before any device data or device scalars exist, so the
  !  device-resident apply_periodic_bc_x cannot be used here): ghost 1 <- cell nxg-2, cells nxg-1, nxg <- 2, 3 (the periodic
  !  duplicate cell convention of apply_periodic_bc_x), across the partner rank when x is split
  Subroutine phi_wrap_x_host

    Logical        :: is_first, is_last
    Integer(Int32) :: partner
    Real   (Int64), Allocatable :: b2(:,:,:), b1(:,:)

    Call x_periodic_partner(is_first, is_last, partner)
    If ( is_first .And. is_last ) Then
       phi(1,    :,:) = phi(nxg-2,:,:)
       phi(nxg-1,:,:) = phi(2,    :,:)
       phi(nxg,  :,:) = phi(3,    :,:)
    Else If ( is_first .Or. is_last ) Then
       Allocate( b2(2,nyg,nzg), b1(nyg,nzg) )
       If ( is_first ) Then
          b2(1,:,:) = phi(2,:,:);  b2(2,:,:) = phi(3,:,:)
          Call MPI_Sendrecv( b2, 2*nyg*nzg, MPI_real8, partner, 41, b1, nyg*nzg, MPI_real8, partner, 42, &
                             MPI_COMM_WORLD, istat, ierr )
          phi(1,:,:) = b1
       Else
          b1 = phi(nxg-2,:,:)
          Call MPI_Sendrecv( b1, nyg*nzg, MPI_real8, partner, 42, b2, 2*nyg*nzg, MPI_real8, partner, 41, &
                             MPI_COMM_WORLD, istat, ierr )
          phi(nxg-1,:,:) = b2(1,:,:);  phi(nxg,:,:) = b2(2,:,:)
       End If
       Deallocate( b2, b1 )
    End If

  End Subroutine phi_wrap_x_host

  !> Periodic wrap of phi in z, host-only (see phi_wrap_x_host)
  Subroutine phi_wrap_z_host

    Logical        :: is_first, is_last
    Integer(Int32) :: partner
    Real   (Int64), Allocatable :: b2(:,:,:), b1(:,:)

    Call z_periodic_partner(is_first, is_last, partner)
    If ( is_first .And. is_last ) Then
       phi(:,:,1    ) = phi(:,:,nzg-2)
       phi(:,:,nzg-1) = phi(:,:,2    )
       phi(:,:,nzg  ) = phi(:,:,3    )
    Else If ( is_first .Or. is_last ) Then
       Allocate( b2(nxg,nyg,2), b1(nxg,nyg) )
       If ( is_first ) Then
          b2(:,:,1) = phi(:,:,2);  b2(:,:,2) = phi(:,:,3)
          Call MPI_Sendrecv( b2, 2*nxg*nyg, MPI_real8, partner, 43, b1, nxg*nyg, MPI_real8, partner, 44, &
                             MPI_COMM_WORLD, istat, ierr )
          phi(:,:,1) = b1
       Else
          b1 = phi(:,:,nzg-2)
          Call MPI_Sendrecv( b1, nxg*nyg, MPI_real8, partner, 44, b2, 2*nxg*nyg, MPI_real8, partner, 43, &
                             MPI_COMM_WORLD, istat, ierr )
          phi(:,:,nzg-1) = b2(:,:,1);  phi(:,:,nzg) = b2(:,:,2)
       End If
       Deallocate( b2, b1 )
    End If

  End Subroutine phi_wrap_z_host

  !> Read a distributed cell-centre scalar field from file: (nxg_global,nyg_global,nzm_global) Real(8), column-major, big-endian (x,y already ghosted in-file, z interior-only, same convention as xg_global/kg-based reads elsewhere)
  Subroutine read_distributed_scalar_field(filename, field)

    Character(*), Intent(In) :: filename
    Real(Int64), Dimension(nxg,nyg,nzg), Intent(InOut) :: field

    Integer(Int32) :: iproc, nxge_r, nzge_r, n_interior
    Integer(Int32) :: sdf_unit
    Integer(Int64) :: file_size, expected_size, expected_elements
    Real   (Int64), Allocatable :: global_field(:,:,:), tmp_read(:), send_buf(:,:,:)

    ! Rank iproc owns global x-columns ig1_global(iproc):ig2_global(iproc) (already ghosted in-file)
    ! and interior z-planes kg1_global(iproc):kg2_global(iproc)-2.
    If ( myid==0 ) Then

       expected_elements = Int(nxg_global, Int64) * Int(nyg_global, Int64) * Int(nzm_global, Int64)
       expected_size = expected_elements * Int(storage_size(1d0)/8, Int64)
       Inquire(file=Trim(filename), size=file_size)
       If ( file_size /= expected_size ) Then
          Write(error_unit,*) 'IBM: file size mismatch reading ', Trim(filename)
          Write(error_unit,*) '  expected ', expected_elements, ' elements (', expected_size, ' bytes),', &
               ' found ', file_size, ' bytes'
          Call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
       End If

       Open(newunit=sdf_unit, file=Trim(filename), access='stream', form='unformatted', action='read', convert='big_endian')
       Allocate( tmp_read(expected_elements) )
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
    Integer(Int32) :: objid_max_local, objid_max_global

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

    ! Ghost BCs for phi at domain boundaries: true periodic wrap when x_bc_type/z_bc_type
    ! select periodic (matches solve_poisson_equation's treatment of P, projection.f90),
    ! Neumann (zero-gradient) otherwise. x=1/nxg (resp. z=1/nzg) is an inter-rank seam
    ! (already correctly populated straight from the global file slice) on any rank that
    ! doesn't own the true global domain edge, so the Neumann branch skips it there.
    ! is_first_x/is_last_x are needed below regardless of branch (passed to
    ! smooth_ibm_corners), so resolve them unconditionally -- cheap, no communication.
    Call x_periodic_partner(is_first_x, is_last_x, partner_x)
    If ( x_bc_type == 0 ) Then
       Call phi_wrap_x_host
    Else
       If ( is_first_x ) phi(1,:,:)   = phi(2,:,:)
       If ( is_last_x  ) phi(nxg,:,:) = phi(nxg-1,:,:)
    End If
    phi(:,1,:)   = phi(:,2,:)
    phi(:,nyg,:) = phi(:,nyg-1,:)
    If ( z_bc_type == 0 ) Then
       Call phi_wrap_z_host
    Else
       phi(:,:,1)   = phi(:,:,2)
       phi(:,:,nzg) = phi(:,:,nzg-1)
    End If
    ! Overwrite interior-rank z ghost planes with actual neighbour values
    Call exchange_phi_ghost_planes

    ! Optional SDF corner-rounding (see smooth_ibm): must run before Umask_cc/ghost lists are
    ! derived from phi, so the rounded geometry is what the ghost-cell lists actually see
    Call smooth_ibm_corners(is_first_x, is_last_x)

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

       ! every ghost_*_objid assignment below clamps into [0,max_ibm_objects]; catch an
       ! out-of-range ID here instead of silently mis-mapping it onto the wrong object's
       ! per-object BC/roughness settings
       objid_max_local = Nint(Maxval(ibm_obj_id))
       Call MPI_Allreduce(objid_max_local, objid_max_global, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
       If ( objid_max_global > max_ibm_objects ) Then
          If ( myid==0 ) Write(*,'(A,I0,A,I0,A)') ' ERROR: ibm_objid_file contains object ID ', &
               objid_max_global, ', but max_ibm_objects = ', max_ibm_objects, &
               ' (per-object BC/roughness arrays only go that far)'
          Call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
       End If
    End If

  End Subroutine read_phi_from_sdf_file

  !> Round sharp SDF corners by smooth_ibm passes of 6-point Jacobi averaging on phi (no-op
  !> when smooth_ibm<=0, the default -- exact original sharp geometry). At a sharp corner the
  !> ghost-cell surface normal (compute_normal_at_cc/_face_*) is discontinuous, and the
  !> resulting rough near-body reconstruction has nothing to damp it once injected into the
  !> (non-dissipative) skew-symmetric advection scheme, which is what produces Gibbs-like
  !> ringing there. Smoothing phi directly blunts the corner by ~smooth_ibm grid cells, giving
  !> the ghost-cell method a continuous normal to work with at the cost of that much geometric
  !> fidelity right at the corner.
  Subroutine smooth_ibm_corners(is_first_x, is_last_x)

    Logical, Intent(In) :: is_first_x, is_last_x

    Real   (Int64), Allocatable :: phi_new(:,:,:)
    Integer(Int32) :: i, j, k, n

    If ( smooth_ibm <= 0 ) Return

    If ( myid==0 ) Write(*,'(A,I0,A)') ' IBM: rounding SDF corners (smooth_ibm = ', smooth_ibm, ' passes)...'

    Allocate ( phi_new(nxg,nyg,nzg) )

    Do n = 1, smooth_ibm
       Do k = 2, nzg-1
          Do j = 2, nyg-1
             Do i = 2, nxg-1
                phi_new(i,j,k) = ( phi(i,j,k) + phi(i-1,j,k) + phi(i+1,j,k) &
                                  + phi(i,j-1,k) + phi(i,j+1,k) &
                                  + phi(i,j,k-1) + phi(i,j,k+1) ) / 7d0
             End Do
          End Do
       End Do
       phi(2:nxg-1,2:nyg-1,2:nzg-1) = phi_new(2:nxg-1,2:nyg-1,2:nzg-1)

       ! Refresh ghosts before the next pass: true periodic wrap when x_bc_type/z_bc_type
       ! select periodic, Neumann at true domain boundaries otherwise, MPI exchange at
       ! interior-rank seams (both x and z; y is never decomposed)
       If ( x_bc_type == 0 ) Then
          Call phi_wrap_x_host
       Else
          If ( is_first_x ) phi(1,:,:)   = phi(2,:,:)
          If ( is_last_x  ) phi(nxg,:,:) = phi(nxg-1,:,:)
       End If
       phi(:,1,:)   = phi(:,2,:)
       phi(:,nyg,:) = phi(:,nyg-1,:)
       If ( z_bc_type == 0 ) Then
          Call phi_wrap_z_host
       Else
          phi(:,:,1)   = phi(:,:,2)
          phi(:,:,nzg) = phi(:,:,nzg-1)
       End If
       Call exchange_phi_x_ghost_planes
       Call exchange_phi_ghost_planes
    End Do

    Deallocate ( phi_new )

  End Subroutine smooth_ibm_corners

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
             If ( 0.5d0*(phie(i,j,k)+phie(i+1,j,k)) < 0d0 ) Then
                If ( 0.5d0*(phie(i-1,j,k)+phie(i,  j,  k  )) >= 0d0 .Or. &
                     0.5d0*(phie(i+1,j,k)+phie(i+2,j,  k  )) >= 0d0 .Or. &
                     0.5d0*(phie(i,  j-1,k)+phie(i+1,j-1,k)) >= 0d0 .Or. &
                     0.5d0*(phie(i,  j+1,k)+phie(i+1,j+1,k)) >= 0d0 .Or. &
                     0.5d0*(phie(i,  j,  k-1)+phie(i+1,j,k-1)) >= 0d0 .Or. &
                     0.5d0*(phie(i,  j,  k+1)+phie(i+1,j,k+1)) >= 0d0 ) Then
                   Call compute_normal_at_face_u(i,j,k, nx_,ny_,nz_)
                   dGB = Abs( 0.5d0*(phie(i,j,k)+phie(i+1,j,k)) )
                   dGI = Max( 2d0*dGB, Real(n_image_layers,8)*dymin )
                   xI  = xe(i)  + dGI*nx_
                   yI  = yg(j) + dGI*ny_
                   zI  = zge(k) + dGI*nz_
                   Call find_stencil_u(xI, yI, zI, ii, jj, kk)
                   If ( .Not. stencil_out .And. ii>=ext_fx_lo .And. ii<=ext_fx_hi .And. jj>=1 .And. jj<=nyg-1 .And. &
                        kk>=ext_cz_lo .And. kk<=ext_cz_hi ) Then
                      If ( Umaske(ii,  jj,  kk  ) > 0.5d0 .And. Umaske(ii+1,jj,  kk  ) > 0.5d0 .And. &
                           Umaske(ii,  jj+1,kk  ) > 0.5d0 .And. Umaske(ii+1,jj+1,kk  ) > 0.5d0 .And. &
                           Umaske(ii,  jj,  kk+1) > 0.5d0 .And. Umaske(ii+1,jj,  kk+1) > 0.5d0 .And. &
                           Umaske(ii,  jj+1,kk+1) > 0.5d0 .And. Umaske(ii+1,jj+1,kk+1) > 0.5d0 ) Then
                         ng = ng + 1
                      Else
                         nd = nd + 1   ! image stencil clips solid — concave corner
                         ibm_drop_solid(1) = ibm_drop_solid(1) + 1
                      End If
                   Else
                      nd = nd + 1      ! image outside domain bounds
                      ibm_drop_outside(1) = ibm_drop_outside(1) + 1
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
    Allocate ( ghost_u_dGI(  n_ghost_u) )
    Allocate ( ghost_u_xB(3, n_ghost_u) )
    Allocate ( ghost_u_img_cc(3, n_ghost_u) )
    Allocate ( ghost_u_wgt_cc(8, n_ghost_u) )
    Allocate ( ghost_u_objid(n_ghost_u) )

    ! Pass 2: identical traversal, fill ghost arrays directly
    ng = 0
    Do k = 2, nzg-1
       Do j = 2, nyg-1
          Do i = 2, nx-1
             If ( 0.5d0*(phie(i,j,k)+phie(i+1,j,k)) < 0d0 ) Then
                If ( 0.5d0*(phie(i-1,j,k)+phie(i,  j,  k  )) >= 0d0 .Or. &
                     0.5d0*(phie(i+1,j,k)+phie(i+2,j,  k  )) >= 0d0 .Or. &
                     0.5d0*(phie(i,  j-1,k)+phie(i+1,j-1,k)) >= 0d0 .Or. &
                     0.5d0*(phie(i,  j+1,k)+phie(i+1,j+1,k)) >= 0d0 .Or. &
                     0.5d0*(phie(i,  j,  k-1)+phie(i+1,j,k-1)) >= 0d0 .Or. &
                     0.5d0*(phie(i,  j,  k+1)+phie(i+1,j,k+1)) >= 0d0 ) Then

                   xGc = xe(i)
                   yGc = yg(j)
                   zGc = zge(k)

                   Call compute_normal_at_face_u(i,j,k, nx_,ny_,nz_)

                   dGB = Abs( 0.5d0*(phie(i,j,k)+phie(i+1,j,k)) )

                   xB = xGc + dGB*nx_
                   yB = yGc + dGB*ny_
                   zB = zGc + dGB*nz_

                   dGI = Max( 2d0*dGB, Real(n_image_layers,8)*dymin )
                   xI  = xGc + dGI*nx_
                   yI  = yGc + dGI*ny_
                   zI  = zGc + dGI*nz_

                   Call find_stencil_u(xI, yI, zI, ii, jj, kk)

                   If ( .Not. stencil_out .And. ii>=ext_fx_lo .And. ii<=ext_fx_hi .And. jj>=1 .And. jj<=nyg-1 .And. &
                        kk>=ext_cz_lo .And. kk<=ext_cz_hi ) Then
                      If ( Umaske(ii,  jj,  kk  ) > 0.5d0 .And. Umaske(ii+1,jj,  kk  ) > 0.5d0 .And. &
                           Umaske(ii,  jj+1,kk  ) > 0.5d0 .And. Umaske(ii+1,jj+1,kk  ) > 0.5d0 .And. &
                           Umaske(ii,  jj,  kk+1) > 0.5d0 .And. Umaske(ii+1,jj,  kk+1) > 0.5d0 .And. &
                           Umaske(ii,  jj+1,kk+1) > 0.5d0 .And. Umaske(ii+1,jj+1,kk+1) > 0.5d0 ) Then

                         ng = ng + 1
                         ghost_u_idx(1,ng) = i
                         ghost_u_idx(2,ng) = j
                         ghost_u_idx(3,ng) = k
                         ghost_u_dGI(ng)   = dGI
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
                         yref_ = Max( Abs( (xge(ir)-xB)*nx_ + (yg(jr)-yB)*ny_ + (zge(kr)-zB)*nz_ ), 1d-14 )
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
       ghost_u_dGB(ng)   = Abs( 0.5d0*(phie(i,j,k)+phie(i+1,j,k)) )
       ghost_u_xB(1,ng)  = xe(i)  + ghost_u_dGB(ng)*ghost_u_nrm(1,ng)
       ghost_u_xB(2,ng)  = yg(j)  + ghost_u_dGB(ng)*ghost_u_nrm(2,ng)
       ghost_u_xB(3,ng)  = zge(k)  + ghost_u_dGB(ng)*ghost_u_nrm(3,ng)
       ! Cell-centre stencil at image point I = xB + (dGI-dGB)*nrm = xG + dGI*nrm
       ! (same xI as velocity image point, on cell-centre grid for P interpolation)
       Call find_stencil_centre( ghost_u_xB(1,ng) + (ghost_u_dGI(ng)-ghost_u_dGB(ng))*ghost_u_nrm(1,ng), &
                                 ghost_u_xB(2,ng) + (ghost_u_dGI(ng)-ghost_u_dGB(ng))*ghost_u_nrm(2,ng), &
                                 ghost_u_xB(3,ng) + (ghost_u_dGI(ng)-ghost_u_dGB(ng))*ghost_u_nrm(3,ng), &
                                 ii, jj, kk )
       Call trilinear_weights( ghost_u_xB(1,ng) + (ghost_u_dGI(ng)-ghost_u_dGB(ng))*ghost_u_nrm(1,ng), &
                               ghost_u_xB(2,ng) + (ghost_u_dGI(ng)-ghost_u_dGB(ng))*ghost_u_nrm(2,ng), &
                               ghost_u_xB(3,ng) + (ghost_u_dGI(ng)-ghost_u_dGB(ng))*ghost_u_nrm(3,ng), &
                               ii, jj, kk, ghost_u_wgt_cc(1:8,ng) )
       ghost_u_img_cc(1,ng) = ii;  ghost_u_img_cc(2,ng) = jj;  ghost_u_img_cc(3,ng) = kk
    End Do

  End Subroutine build_ghost_list_u

  !  Build ghost-cell list for V (y-faces)
  Subroutine build_ghost_list_v

    Integer(Int32) :: i, j, k, ng, nd
    Integer(Int32) :: ii, jj, kk, ir, jr, kr
    Real   (Int64) :: nx_, ny_, nz_
    Real   (Int64) :: xGc, yGc, zGc, dGB, dGI, xB, yB, zB, xI, yI, zI, yref_

    ! Pass 1: count ghost cells (no storage)
    ng = 0;  nd = 0
    Do k = 2, nzg-1
       Do j = 2, ny-1
          Do i = 2, nxg-1
             If ( phi_v(i,j,k) < 0d0 ) Then
                If ( phi_v(i-1,j,k) >= 0d0 .Or. &
                     phi_v(i+1,j,k) >= 0d0 .Or. &
                     phi_v(i,j-1,k)   >= 0d0 .Or. &
                     phi_v(i,j+1,k)   >= 0d0 .Or. &
                     phi_v(i,j,k-1) >= 0d0 .Or. &
                     phi_v(i,j,k+1) >= 0d0 ) Then
                   Call compute_normal_at_face_v(i,j,k, nx_,ny_,nz_)
                   dGB = Abs( phi_v(i,j,k) )
                   dGI = Max( 2d0*dGB, Real(n_image_layers,8)*dymin )
                   xI = xge(i) + dGI*nx_
                   yI = y(j)  + dGI*ny_
                   zI = zge(k) + dGI*nz_
                   Call find_stencil_v(xI, yI, zI, ii, jj, kk)
                   If ( .Not. stencil_out .And. ii>=ext_cx_lo .And. ii<=ext_cx_hi .And. jj>=2 .And. jj<=ny-1 .And. &
                        kk>=ext_cz_lo .And. kk<=ext_cz_hi ) Then
                      If ( Umaske(ii,  jj,  kk  ) > 0.5d0 .And. Umaske(ii+1,jj,  kk  ) > 0.5d0 .And. &
                           Umaske(ii,  jj+1,kk  ) > 0.5d0 .And. Umaske(ii+1,jj+1,kk  ) > 0.5d0 .And. &
                           Umaske(ii,  jj,  kk+1) > 0.5d0 .And. Umaske(ii+1,jj,  kk+1) > 0.5d0 .And. &
                           Umaske(ii,  jj+1,kk+1) > 0.5d0 .And. Umaske(ii+1,jj+1,kk+1) > 0.5d0 ) Then
                         ng = ng + 1
                      Else
                         nd = nd + 1   ! image stencil clips solid — concave corner
                         ibm_drop_solid(2) = ibm_drop_solid(2) + 1
                      End If
                   Else
                      nd = nd + 1      ! image outside domain bounds
                      ibm_drop_outside(2) = ibm_drop_outside(2) + 1
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
    Allocate ( ghost_v_dGI(  n_ghost_v) )
    Allocate ( ghost_v_xB(3, n_ghost_v) )
    Allocate ( ghost_v_img_cc(3, n_ghost_v) )
    Allocate ( ghost_v_wgt_cc(8, n_ghost_v) )
    Allocate ( ghost_v_objid(n_ghost_v) )

    ! Pass 2: identical traversal, fill ghost arrays directly
    ng = 0
    Do k = 2, nzg-1
       Do j = 2, ny-1
          Do i = 2, nxg-1
             If ( phi_v(i,j,k) < 0d0 ) Then
                If ( phi_v(i-1,j,k) >= 0d0 .Or. &
                     phi_v(i+1,j,k) >= 0d0 .Or. &
                     phi_v(i,j-1,k)   >= 0d0 .Or. &
                     phi_v(i,j+1,k)   >= 0d0 .Or. &
                     phi_v(i,j,k-1) >= 0d0 .Or. &
                     phi_v(i,j,k+1) >= 0d0 ) Then

                   xGc = xge(i)
                   yGc = y(j)
                   zGc = zge(k)

                   Call compute_normal_at_face_v(i,j,k, nx_,ny_,nz_)

                   dGB = Abs( phi_v(i,j,k) )
                   dGI = Max( 2d0*dGB, Real(n_image_layers,8)*dymin )
                   xB = xGc + dGB*nx_;  yB = yGc + dGB*ny_;  zB = zGc + dGB*nz_
                   xI = xGc + dGI*nx_;  yI = yGc + dGI*ny_;  zI = zGc + dGI*nz_

                   Call find_stencil_v(xI, yI, zI, ii, jj, kk)
                   If ( .Not. stencil_out .And. ii>=ext_cx_lo .And. ii<=ext_cx_hi .And. jj>=2 .And. jj<=ny-1 .And. &
                        kk>=ext_cz_lo .And. kk<=ext_cz_hi ) Then
                      If ( Umaske(ii,  jj,  kk  ) > 0.5d0 .And. Umaske(ii+1,jj,  kk  ) > 0.5d0 .And. &
                           Umaske(ii,  jj+1,kk  ) > 0.5d0 .And. Umaske(ii+1,jj+1,kk  ) > 0.5d0 .And. &
                           Umaske(ii,  jj,  kk+1) > 0.5d0 .And. Umaske(ii+1,jj,  kk+1) > 0.5d0 .And. &
                           Umaske(ii,  jj+1,kk+1) > 0.5d0 .And. Umaske(ii+1,jj+1,kk+1) > 0.5d0 ) Then
                         ng = ng + 1
                         ghost_v_idx(1,ng)=i; ghost_v_idx(2,ng)=j; ghost_v_idx(3,ng)=k
                         ghost_v_dGI(ng) = dGI
                         ghost_v_img(1,ng)=ii; ghost_v_img(2,ng)=jj; ghost_v_img(3,ng)=kk
                         Call trilinear_weights_v(xI, yI, zI, ii, jj, kk, ghost_v_wgt(1:8,ng))
                         ghost_v_nrm(:,ng)=[nx_,ny_,nz_]
                         Call find_stencil_centre(xB+Real(n_image_layers,8)*dymin*nx_, &
                                                  yB+Real(n_image_layers,8)*dymin*ny_, &
                                                  zB+Real(n_image_layers,8)*dymin*nz_, &
                                                  ir,jr,kr)
                         ghost_v_yref(ng) = Max( Abs( (xge(ir)-xB)*nx_ + (yg(jr)-yB)*ny_ + (zge(kr)-zB)*nz_ ), 1d-14 )
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
       ghost_v_dGB(ng)   = Abs( phi_v(i,j,k) )
       ghost_v_xB(1,ng)  = xge(i)  + ghost_v_dGB(ng)*ghost_v_nrm(1,ng)
       ghost_v_xB(2,ng)  = y (j)  + ghost_v_dGB(ng)*ghost_v_nrm(2,ng)
       ghost_v_xB(3,ng)  = zge(k)  + ghost_v_dGB(ng)*ghost_v_nrm(3,ng)
       Call find_stencil_centre( ghost_v_xB(1,ng) + (ghost_v_dGI(ng)-ghost_v_dGB(ng))*ghost_v_nrm(1,ng), &
                                 ghost_v_xB(2,ng) + (ghost_v_dGI(ng)-ghost_v_dGB(ng))*ghost_v_nrm(2,ng), &
                                 ghost_v_xB(3,ng) + (ghost_v_dGI(ng)-ghost_v_dGB(ng))*ghost_v_nrm(3,ng), &
                                 ii, jj, kk )
       Call trilinear_weights( ghost_v_xB(1,ng) + (ghost_v_dGI(ng)-ghost_v_dGB(ng))*ghost_v_nrm(1,ng), &
                               ghost_v_xB(2,ng) + (ghost_v_dGI(ng)-ghost_v_dGB(ng))*ghost_v_nrm(2,ng), &
                               ghost_v_xB(3,ng) + (ghost_v_dGI(ng)-ghost_v_dGB(ng))*ghost_v_nrm(3,ng), &
                               ii, jj, kk, ghost_v_wgt_cc(1:8,ng) )
       ghost_v_img_cc(1,ng) = ii;  ghost_v_img_cc(2,ng) = jj;  ghost_v_img_cc(3,ng) = kk
    End Do

  End Subroutine build_ghost_list_v

  !  Build ghost-cell list for W (z-faces)
  Subroutine build_ghost_list_w

    Integer(Int32) :: i, j, k, ng, nd
    Integer(Int32) :: ii, jj, kk, ir, jr, kr
    Real   (Int64) :: nx_, ny_, nz_
    Real   (Int64) :: xGc, yGc, zGc, dGB, dGI, xB, yB, zB, xI, yI, zI, yref_

    ! Pass 1: count ghost cells (no storage)
    ng = 0;  nd = 0
    Do k = 2, nz-1
       Do j = 2, nyg-1
          Do i = 2, nxg-1
             If ( phi_w(i,j,k) < 0d0 ) Then
                If ( 0.5d0*(phie(i-1,j,k)+phie(i-1,j,k+1)) >= 0d0 .Or. &
                     0.5d0*(phie(i+1,j,k)+phie(i+1,j,k+1)) >= 0d0 .Or. &
                     0.5d0*(phie(i,j-1,k)+phie(i,j-1,k+1)) >= 0d0 .Or. &
                     0.5d0*(phie(i,j+1,k)+phie(i,j+1,k+1)) >= 0d0 .Or. &
                     phi_w(i,j,k-1) >= 0d0 .Or. &
                     phi_w(i,j,k+1) >= 0d0 ) Then
                   Call compute_normal_at_face_w(i,j,k, nx_,ny_,nz_)
                   dGB = Abs( phi_w(i,j,k) )
                   dGI = Max( 2d0*dGB, Real(n_image_layers,8)*dymin )
                   xI = xge(i) + dGI*nx_
                   yI = yg(j) + dGI*ny_
                   zI = ze(k)  + dGI*nz_
                   Call find_stencil_w(xI, yI, zI, ii, jj, kk)
                   If ( .Not. stencil_out .And. ii>=ext_cx_lo .And. ii<=ext_cx_hi .And. jj>=1 .And. jj<=nyg-1 .And. &
                        kk>=ext_fz_lo .And. kk<=ext_fz_hi ) Then
                      If ( Umaske(ii,  jj,  kk  ) > 0.5d0 .And. Umaske(ii+1,jj,  kk  ) > 0.5d0 .And. &
                           Umaske(ii,  jj+1,kk  ) > 0.5d0 .And. Umaske(ii+1,jj+1,kk  ) > 0.5d0 .And. &
                           Umaske(ii,  jj,  kk+1) > 0.5d0 .And. Umaske(ii+1,jj,  kk+1) > 0.5d0 .And. &
                           Umaske(ii,  jj+1,kk+1) > 0.5d0 .And. Umaske(ii+1,jj+1,kk+1) > 0.5d0 ) Then
                         ng = ng + 1
                      Else
                         nd = nd + 1   ! image stencil clips solid — concave corner
                         ibm_drop_solid(3) = ibm_drop_solid(3) + 1
                      End If
                   Else
                      nd = nd + 1      ! image outside domain bounds
                      ibm_drop_outside(3) = ibm_drop_outside(3) + 1
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
    Allocate ( ghost_w_dGI(  n_ghost_w) )
    Allocate ( ghost_w_xB(3, n_ghost_w) )
    Allocate ( ghost_w_img_cc(3, n_ghost_w) )
    Allocate ( ghost_w_wgt_cc(8, n_ghost_w) )
    Allocate ( ghost_w_objid(n_ghost_w) )

    ! Pass 2: identical traversal, fill ghost arrays directly
    ng = 0
    Do k = 2, nz-1
       Do j = 2, nyg-1
          Do i = 2, nxg-1
             If ( phi_w(i,j,k) < 0d0 ) Then
                If ( 0.5d0*(phie(i-1,j,k)+phie(i-1,j,k+1)) >= 0d0 .Or. &
                     0.5d0*(phie(i+1,j,k)+phie(i+1,j,k+1)) >= 0d0 .Or. &
                     0.5d0*(phie(i,j-1,k)+phie(i,j-1,k+1)) >= 0d0 .Or. &
                     0.5d0*(phie(i,j+1,k)+phie(i,j+1,k+1)) >= 0d0 .Or. &
                     phi_w(i,j,k-1) >= 0d0 .Or. &
                     phi_w(i,j,k+1) >= 0d0 ) Then

                   xGc = xge(i)
                   yGc = yg(j)
                   zGc = ze(k)

                   Call compute_normal_at_face_w(i,j,k, nx_,ny_,nz_)

                   dGB = Abs( phi_w(i,j,k) )
                   dGI = Max( 2d0*dGB, Real(n_image_layers,8)*dymin )
                   xB = xGc+dGB*nx_;  yB = yGc+dGB*ny_;  zB = zGc+dGB*nz_
                   xI = xGc+dGI*nx_;  yI = yGc+dGI*ny_;  zI = zGc+dGI*nz_

                   Call find_stencil_w(xI, yI, zI, ii, jj, kk)
                   If ( .Not. stencil_out .And. ii>=ext_cx_lo .And. ii<=ext_cx_hi .And. jj>=1 .And. jj<=nyg-1 .And. &
                        kk>=ext_fz_lo .And. kk<=ext_fz_hi ) Then
                      If ( Umaske(ii,  jj,  kk  ) > 0.5d0 .And. Umaske(ii+1,jj,  kk  ) > 0.5d0 .And. &
                           Umaske(ii,  jj+1,kk  ) > 0.5d0 .And. Umaske(ii+1,jj+1,kk  ) > 0.5d0 .And. &
                           Umaske(ii,  jj,  kk+1) > 0.5d0 .And. Umaske(ii+1,jj,  kk+1) > 0.5d0 .And. &
                           Umaske(ii,  jj+1,kk+1) > 0.5d0 .And. Umaske(ii+1,jj+1,kk+1) > 0.5d0 ) Then
                         ng = ng + 1
                         ghost_w_idx(1,ng)=i; ghost_w_idx(2,ng)=j; ghost_w_idx(3,ng)=k
                         ghost_w_dGI(ng) = dGI
                         ghost_w_img(1,ng)=ii; ghost_w_img(2,ng)=jj; ghost_w_img(3,ng)=kk
                         Call trilinear_weights_w(xI, yI, zI, ii, jj, kk, ghost_w_wgt(1:8,ng))
                         ghost_w_nrm(:,ng)=[nx_,ny_,nz_]
                         Call find_stencil_centre(xB+Real(n_image_layers,8)*dymin*nx_, &
                                                  yB+Real(n_image_layers,8)*dymin*ny_, &
                                                  zB+Real(n_image_layers,8)*dymin*nz_, &
                                                  ir,jr,kr)
                         ghost_w_yref(ng) = Max( Abs( (xge(ir)-xB)*nx_ + (yg(jr)-yB)*ny_ + (zge(kr)-zB)*nz_ ), 1d-14 )
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
       ghost_w_dGB(ng)   = Abs( phi_w(i,j,k) )
       ghost_w_xB(1,ng)  = xge(i)  + ghost_w_dGB(ng)*ghost_w_nrm(1,ng)
       ghost_w_xB(2,ng)  = yg(j)  + ghost_w_dGB(ng)*ghost_w_nrm(2,ng)
       ghost_w_xB(3,ng)  = ze(k)  + ghost_w_dGB(ng)*ghost_w_nrm(3,ng)
       Call find_stencil_centre( ghost_w_xB(1,ng) + (ghost_w_dGI(ng)-ghost_w_dGB(ng))*ghost_w_nrm(1,ng), &
                                 ghost_w_xB(2,ng) + (ghost_w_dGI(ng)-ghost_w_dGB(ng))*ghost_w_nrm(2,ng), &
                                 ghost_w_xB(3,ng) + (ghost_w_dGI(ng)-ghost_w_dGB(ng))*ghost_w_nrm(3,ng), &
                                 ii, jj, kk )
       Call trilinear_weights( ghost_w_xB(1,ng) + (ghost_w_dGI(ng)-ghost_w_dGB(ng))*ghost_w_nrm(1,ng), &
                               ghost_w_xB(2,ng) + (ghost_w_dGI(ng)-ghost_w_dGB(ng))*ghost_w_nrm(2,ng), &
                               ghost_w_xB(3,ng) + (ghost_w_dGI(ng)-ghost_w_dGB(ng))*ghost_w_nrm(3,ng), &
                               ii, jj, kk, ghost_w_wgt_cc(1:8,ng) )
       ghost_w_img_cc(1,ng) = ii;  ghost_w_img_cc(2,ng) = jj;  ghost_w_img_cc(3,ng) = kk
    End Do

  End Subroutine build_ghost_list_w

  !> Debug (DOPAMINE_TRACE_DIR set, see debug_trace.f90): write each rank's ghost-cell lists as global (i,j,k) so that
  !  layouts can be compared cell by cell (tests/regression/ibm_ghost_diff.py)
  Subroutine trace_ghost_lists

    Character(Len=1024) :: dir, fname
    Integer(Int32)      :: length, stat, n, u
    Integer(Int32)      :: gi, gk

    Call Get_Environment_Variable('DOPAMINE_TRACE_DIR', dir, length, stat)
    If ( stat /= 0 .Or. length == 0 ) Return

    Write(fname,'(A,A,I4.4)') dir(1:length), '/ibm_ghosts.rank', myid
    Open(newunit=u, file=Trim(fname), status='replace', action='write')
    Do n = 1, n_ghost_u
       gi = i1_global(myid) + ghost_u_idx(1,n) - 1;  gk = kg1_global(myid) + ghost_u_idx(3,n) - 1
       Write(u,'(A,3(1X,I0))') 'U', gi, ghost_u_idx(2,n), gk
    End Do
    Do n = 1, n_ghost_v
       gi = ig1_global(myid) + ghost_v_idx(1,n) - 1;  gk = kg1_global(myid) + ghost_v_idx(3,n) - 1
       Write(u,'(A,3(1X,I0))') 'V', gi, ghost_v_idx(2,n), gk
    End Do
    Do n = 1, n_ghost_w
       gi = ig1_global(myid) + ghost_w_idx(1,n) - 1;  gk = k1_global(myid) + ghost_w_idx(3,n) - 1
       Write(u,'(A,3(1X,I0))') 'W', gi, ghost_w_idx(2,n), gk
    End Do
    Close(u)

    ! the local SDF with its global offsets, to compare layouts cell by cell (ibm_ghost_diff.py --phi)
    Write(fname,'(A,A,I4.4)') dir(1:length), '/ibm_phi.rank', myid
    Open(newunit=u, file=Trim(fname), status='replace', access='stream', form='unformatted', action='write')
    Write(u) nxg, nyg, nzg, ig1_global(myid), kg1_global(myid)
    Write(u) phi
    Close(u)

  End Subroutine trace_ghost_lists

  !> Fe(1-E:n1+E, :, 1-E:n3+E) = F plus E extra planes in x and z from the neighbour rank / periodic partner / edge (the rules of
  !  halo_pad.pad_field). Device-aware: the copy and the pack/unpack run on the device in GPU builds, only the small exchanged
  !  planes go through the host for MPI.
  Subroutine ibm_fill_ext(F, Fe, n1, n2, n3, xface, zface)

    Integer(Int32), Intent(In)    :: n1, n2, n3
    Real   (Int64), Intent(In)    :: F(n1,n2,n3)
    Real   (Int64), Intent(InOut) :: Fe(1-ibm_E:n1+ibm_E, n2, 1-ibm_E:n3+ibm_E)
    Logical,        Intent(In)    :: xface, zface

    Logical        :: is_first, is_last, per
    Integer(Int32) :: up, down, partner, dst_dn, dst_up, src_dn, src_up, sdn, sup, m, i, j, k, E, cnt, idx

    E = ibm_E

    !$acc parallel loop collapse(3) present(F,Fe)
    Do k = 1, n3
       Do j = 1, n2
          Do i = 1, n1
             Fe(i,j,k) = F(i,j,k)
          End Do
       End Do
    End Do
    !$acc end parallel loop

    !-- x --
    per = ( x_bc_type == 0 )
    Call x_halo_neighbors(up, down)
    Call x_periodic_partner(is_first, is_last, partner)
    If ( is_first .And. is_last ) Then
       If ( per ) Then
          sdn = Merge(3, 4, xface);  sup = Merge(n1-2, n1-3, xface)
          !$acc parallel loop collapse(3) present(Fe)
          Do k = 1, n3
             Do j = 1, n2
                Do m = 1, E
                   Fe(n1+m,j,k) = Fe(sdn+m-1,j,k)
                   Fe(1-m, j,k) = Fe(sup-m+1,j,k)
                End Do
             End Do
          End Do
          !$acc end parallel loop
       Else
          !$acc parallel loop collapse(3) present(Fe)
          Do k = 1, n3
             Do j = 1, n2
                Do m = 1, E
                   Fe(n1+m,j,k) = Fe(n1,j,k)
                   Fe(1-m, j,k) = Fe(1, j,k)
                End Do
             End Do
          End Do
          !$acc end parallel loop
       End If
    Else
       dst_dn = down;  src_up = up;  dst_up = up;  src_dn = down
       sdn = 3;  sup = n1-2
       If ( per .And. is_first ) Then
          dst_dn = partner;  src_dn = partner;  sdn = Merge(3, 4, xface)
       End If
       If ( per .And. is_last ) Then
          dst_up = partner;  src_up = partner;  sup = Merge(n1-2, n1-3, xface)
       End If
       cnt = E*n2*n3
       !$acc parallel loop collapse(3) present(Fe,ext_bs)
       Do k = 1, n3
          Do j = 1, n2
             Do m = 1, E
                ext_bs(m + E*((j-1) + n2*(k-1))) = Fe(sdn+m-1,j,k)
             End Do
          End Do
       End Do
       !$acc end parallel loop
       !$acc update host(ext_bs(1:cnt))
       Call MPI_Sendrecv( ext_bs, cnt, MPI_real8, dst_dn, 71, ext_br, cnt, MPI_real8, src_up, 71, MPI_COMM_WORLD, istat, ierr )
       !$acc update device(ext_br(1:cnt))
       !$acc parallel loop collapse(3) present(Fe,ext_br)
       Do k = 1, n3
          Do j = 1, n2
             Do m = 1, E
                If ( src_up /= MPI_PROC_NULL ) Then
                   Fe(n1+m,j,k) = ext_br(m + E*((j-1) + n2*(k-1)))
                Else
                   Fe(n1+m,j,k) = Fe(n1,j,k)
                End If
             End Do
          End Do
       End Do
       !$acc end parallel loop
       !$acc parallel loop collapse(3) present(Fe,ext_bs)
       Do k = 1, n3
          Do j = 1, n2
             Do m = 1, E
                ext_bs(m + E*((j-1) + n2*(k-1))) = Fe(sup-E+m,j,k)
             End Do
          End Do
       End Do
       !$acc end parallel loop
       !$acc update host(ext_bs(1:cnt))
       Call MPI_Sendrecv( ext_bs, cnt, MPI_real8, dst_up, 72, ext_br, cnt, MPI_real8, src_dn, 72, MPI_COMM_WORLD, istat, ierr )
       !$acc update device(ext_br(1:cnt))
       !$acc parallel loop collapse(3) present(Fe,ext_br)
       Do k = 1, n3
          Do j = 1, n2
             Do m = 1, E
                If ( src_dn /= MPI_PROC_NULL ) Then
                   Fe(m-E,j,k) = ext_br(m + E*((j-1) + n2*(k-1)))
                Else
                   Fe(m-E,j,k) = Fe(1,j,k)
                End If
             End Do
          End Do
       End Do
       !$acc end parallel loop
    End If

    !-- z (over the x-extended extent) --
    per = ( z_bc_type == 0 )
    Call z_halo_neighbors(up, down)
    Call z_periodic_partner(is_first, is_last, partner)
    If ( is_first .And. is_last ) Then
       If ( per ) Then
          sdn = Merge(3, 4, zface);  sup = Merge(n3-2, n3-3, zface)
          !$acc parallel loop collapse(3) present(Fe)
          Do m = 1, E
             Do j = 1, n2
                Do i = 1-E, n1+E
                   Fe(i,j,n3+m) = Fe(i,j,sdn+m-1)
                   Fe(i,j,1-m ) = Fe(i,j,sup-m+1)
                End Do
             End Do
          End Do
          !$acc end parallel loop
       Else
          !$acc parallel loop collapse(3) present(Fe)
          Do m = 1, E
             Do j = 1, n2
                Do i = 1-E, n1+E
                   Fe(i,j,n3+m) = Fe(i,j,n3)
                   Fe(i,j,1-m ) = Fe(i,j,1)
                End Do
             End Do
          End Do
          !$acc end parallel loop
       End If
    Else
       dst_dn = down;  src_up = up;  dst_up = up;  src_dn = down
       sdn = 3;  sup = n3-2
       If ( per .And. is_first ) Then
          dst_dn = partner;  src_dn = partner;  sdn = Merge(3, 4, zface)
       End If
       If ( per .And. is_last ) Then
          dst_up = partner;  src_up = partner;  sup = Merge(n3-2, n3-3, zface)
       End If
       cnt = (n1+2*E)*n2*E
       !$acc parallel loop collapse(3) present(Fe,ext_bs)
       Do m = 1, E
          Do j = 1, n2
             Do i = 1, n1+2*E
                ext_bs(i + (n1+2*E)*((j-1) + n2*(m-1))) = Fe(i-E,j,sdn+m-1)
             End Do
          End Do
       End Do
       !$acc end parallel loop
       !$acc update host(ext_bs(1:cnt))
       Call MPI_Sendrecv( ext_bs, cnt, MPI_real8, dst_dn, 73, ext_br, cnt, MPI_real8, src_up, 73, MPI_COMM_WORLD, istat, ierr )
       !$acc update device(ext_br(1:cnt))
       !$acc parallel loop collapse(3) present(Fe,ext_br)
       Do m = 1, E
          Do j = 1, n2
             Do i = 1, n1+2*E
                If ( src_up /= MPI_PROC_NULL ) Then
                   Fe(i-E,j,n3+m) = ext_br(i + (n1+2*E)*((j-1) + n2*(m-1)))
                Else
                   Fe(i-E,j,n3+m) = Fe(i-E,j,n3)
                End If
             End Do
          End Do
       End Do
       !$acc end parallel loop
       !$acc parallel loop collapse(3) present(Fe,ext_bs)
       Do m = 1, E
          Do j = 1, n2
             Do i = 1, n1+2*E
                ext_bs(i + (n1+2*E)*((j-1) + n2*(m-1))) = Fe(i-E,j,sup-E+m)
             End Do
          End Do
       End Do
       !$acc end parallel loop
       !$acc update host(ext_bs(1:cnt))
       Call MPI_Sendrecv( ext_bs, cnt, MPI_real8, dst_up, 74, ext_br, cnt, MPI_real8, src_dn, 74, MPI_COMM_WORLD, istat, ierr )
       !$acc update device(ext_br(1:cnt))
       !$acc parallel loop collapse(3) present(Fe,ext_br)
       Do m = 1, E
          Do j = 1, n2
             Do i = 1, n1+2*E
                If ( src_dn /= MPI_PROC_NULL ) Then
                   Fe(i-E,j,m-E) = ext_br(i + (n1+2*E)*((j-1) + n2*(m-1)))
                Else
                   Fe(i-E,j,m-E) = Fe(i-E,j,1)
                End If
             End Do
          End Do
       End Do
       !$acc end parallel loop
    End If

  End Subroutine ibm_fill_ext

  !> Refresh Uext/Vext/Wext from the current velocity (device-resident in GPU builds)
  Subroutine ibm_fill_ext_uvw(U_, V_, W_)
    Real(Int64), Dimension(nx, nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg, ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg, nz), Intent(In) :: W_
    Call ibm_fill_ext( U_, Uext, nx,  nyg, nzg, .True.,  .False. )
    Call ibm_fill_ext( V_, Vext, nxg, ny,  nzg, .False., .False. )
    Call ibm_fill_ext( W_, Wext, nxg, nyg, nz,  .False., .True.  )
  End Subroutine ibm_fill_ext_uvw

  !> Refresh Cext from a cell-centred scalar (device-resident in GPU builds)
  Subroutine ibm_fill_ext_cc(C_)
    Real(Int64), Dimension(nxg,nyg,nzg), Intent(In) :: C_
    Call ibm_fill_ext( C_, Cext, nxg, nyg, nzg, .False., .False. )
  End Subroutine ibm_fill_ext_cc

  !> Apply ghost-cell IBM every RK sub-step in place of volume-penalisation
  Subroutine apply_ghost_cell_ibm(U_,V_,W_)

    Real(Int64), Dimension(nx, nyg,nzg), Intent(InOut) :: U_
    Real(Int64), Dimension(nxg, ny,nzg), Intent(InOut) :: V_
    Real(Int64), Dimension(nxg,nyg, nz), Intent(InOut) :: W_

    Integer(Int32) :: n, i, j, k, ee
    Real   (Int64) :: r

    ee = ibm_E
    ! Image-point stencils can reach into the seam halo planes, which are stale after the RK update / projection that precedes this call
    If ( nprocs > 1 ) Call exchange_velocity_halos

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
    !$acc parallel loop collapse(3) present(phi,V_,y,yg)
    Do k = 2, nzg-1
       Do j = 2, ny-1
          Do i = 2, nxg-1
             If ( ( ( yg(j+1) - y(j) )*phi(i,j,k) + ( y(j) - yg(j) )*phi(i,j+1,k) ) / ( yg(j+1) - yg(j) ) < 0d0 ) V_(i,j,k) = 0d0
          End Do
       End Do
    End Do
    !$acc end parallel loop
    !$acc parallel loop collapse(3) present(phi,W_,z,zg)
    Do k = 2, nz-1
       Do j = 2, nyg-1
          Do i = 2, nxg-1
             If ( ( ( zg(k+1) - z(k) )*phi(i,j,k) + ( z(k) - zg(k) )*phi(i,j,k+1) ) / ( zg(k+1) - zg(k) ) < 0d0 ) W_(i,j,k) = 0d0
          End Do
       End Do
    End Do
    !$acc end parallel loop

    ! image-point stencils read the extended copies (extra planes beyond the seam ghost plane)
    Call ibm_fill_ext_uvw(U_, V_, W_)

    !--- U ghost cells: general two-point mirror U_G = (U_wall - r*U_I)/(1-r), r = dGB/dGI
    !    (reduces to the textbook 2*U_wall - U_I when the image sits at the unclamped 2*dGB, r=0.5) ---
    !$acc parallel loop present(Uext,ghost_u_wgt,ghost_u_img,ghost_img_val)
    Do n = 1, n_ghost_u
       ghost_img_val(n) = trilinear_interp_u(Uext, ghost_u_wgt(1:8,n), ghost_u_img(:,n), ee)
    End Do
    !$acc end parallel loop
    !$acc parallel loop present(U_,ghost_u_idx,ghost_u_dGB,ghost_u_dGI,ghost_img_val) private(r)
    Do n = 1, n_ghost_u
       i = ghost_u_idx(1,n);  j = ghost_u_idx(2,n);  k = ghost_u_idx(3,n)
       r = ghost_u_dGB(n) / ghost_u_dGI(n)
       U_(i,j,k) = ( U_wall - r*ghost_img_val(n) ) / ( 1d0 - r )
    End Do
    !$acc end parallel loop

    !--- V ghost cells ---
    !$acc parallel loop present(Vext,ghost_v_wgt,ghost_v_img,ghost_img_val)
    Do n = 1, n_ghost_v
       ghost_img_val(n) = trilinear_interp_v(Vext, ghost_v_wgt(1:8,n), ghost_v_img(:,n), ee)
    End Do
    !$acc end parallel loop
    !$acc parallel loop present(V_,ghost_v_idx,ghost_v_dGB,ghost_v_dGI,ghost_img_val) private(r)
    Do n = 1, n_ghost_v
       i = ghost_v_idx(1,n);  j = ghost_v_idx(2,n);  k = ghost_v_idx(3,n)
       r = ghost_v_dGB(n) / ghost_v_dGI(n)
       V_(i,j,k) = ( V_wall - r*ghost_img_val(n) ) / ( 1d0 - r )
    End Do
    !$acc end parallel loop

    !--- W ghost cells ---
    !$acc parallel loop present(Wext,ghost_w_wgt,ghost_w_img,ghost_img_val)
    Do n = 1, n_ghost_w
       ghost_img_val(n) = trilinear_interp_w(Wext, ghost_w_wgt(1:8,n), ghost_w_img(:,n), ee)
    End Do
    !$acc end parallel loop
    !$acc parallel loop present(W_,ghost_w_idx,ghost_w_dGB,ghost_w_dGI,ghost_img_val) private(r)
    Do n = 1, n_ghost_w
       i = ghost_w_idx(1,n);  j = ghost_w_idx(2,n);  k = ghost_w_idx(3,n)
       r = ghost_w_dGB(n) / ghost_w_dGI(n)
       W_(i,j,k) = ( W_wall - r*ghost_img_val(n) ) / ( 1d0 - r )
    End Do
    !$acc end parallel loop

    ! The neighbour's halo copies of the cells just modified are stale until refreshed; the next SGS/RHS reads them
    If ( nprocs > 1 ) Call exchange_velocity_halos

  End Subroutine apply_ghost_cell_ibm

  !> Ghost-cell thermal condition at the immersed boundary, per solid ID: ibm_T_bc_type(id) 0=adiabatic (T_ghost=T_image), 1=isothermal (mirror against ibm_T_wall(id))
  Subroutine apply_ghost_cell_ibm_scalar(T_)

    Real(Int64), Dimension(nxg,nyg,nzg), Intent(InOut) :: T_

    Integer(Int32) :: n, i, j, k, oid, ee
    Real   (Int64) :: T_I, r

    ee = ibm_E
    Call ibm_fill_ext_cc(T_)
    !$acc parallel loop present(T_,Cext,ghost_cc_idx,ghost_cc_wgt_cc,ghost_cc_img_cc,ghost_cc_objid,ibm_T_bc_type,ibm_T_wall,ghost_cc_dGB,ghost_cc_dGI) private(T_I,oid,r)
    Do n = 1, n_ghost_cc
       i = ghost_cc_idx(1,n);  j = ghost_cc_idx(2,n);  k = ghost_cc_idx(3,n)
       oid = ghost_cc_objid(n)
       T_I = trilinear_interp_p(Cext, ghost_cc_wgt_cc(1:8,n), ghost_cc_img_cc(:,n), ee)
       If ( ibm_T_bc_type(oid) == 1 ) Then
          ! isothermal: general two-point mirror, r = dGB/dGI (reduces to 2*T_wall-T_I at r=0.5)
          r = ghost_cc_dGB(n) / ghost_cc_dGI(n)
          T_(i,j,k) = ( ibm_T_wall(oid) - r*T_I ) / ( 1d0 - r )
       Else
          T_(i,j,k) = T_I                         ! adiabatic: zero-gradient
       End If
    End Do
    !$acc end parallel loop
    ! host consumers (snapshots, restart, probes, particles) read the host copy
    !$acc update host(T_)

  End Subroutine apply_ghost_cell_ibm_scalar

  !> No-flux (zero-gradient) ghost-cell condition for a cell-centred scalar at the immersed boundary: each ghost cell takes its image-point value (used for sediment; solid cells are otherwise held at C=0, which would make the body an absorbing surface)
  Subroutine apply_ghost_cell_ibm_scalar_noflux(C_)

    Real(Int64), Dimension(nxg,nyg,nzg), Intent(InOut) :: C_

    Integer(Int32) :: n, i, j, k, ee

    ee = ibm_E
    Call ibm_fill_ext_cc(C_)
    !$acc parallel loop present(C_,Cext,ghost_cc_idx,ghost_cc_wgt_cc,ghost_cc_img_cc)
    Do n = 1, n_ghost_cc
       i = ghost_cc_idx(1,n);  j = ghost_cc_idx(2,n);  k = ghost_cc_idx(3,n)
       C_(i,j,k) = trilinear_interp_p(Cext, ghost_cc_wgt_cc(1:8,n), ghost_cc_img_cc(:,n), ee)
    End Do
    !$acc end parallel loop
    !$acc update host(C_)

  End Subroutine apply_ghost_cell_ibm_scalar_noflux

  !                     Helper routines

  !  Zero the per-step IBM force accumulators.
  !  Call once per time step before the first RK IBM application.
  Subroutine zero_ibm_stage_accumulators
    ibm_Fx_acc = 0d0
    ibm_Fy_acc = 0d0
    ibm_Fz_acc = 0d0
  End Subroutine zero_ibm_stage_accumulators

  !> Flags for the periodic duplicate cell (cell nxg-1 / nzg-1 on the last rank repeats the first cell): force sums skip it so
  !  the wrapped image of a boundary cell is not counted twice
  Subroutine dup_cell_flags(skip_x, skip_z)
    Logical, Intent(Out) :: skip_x, skip_z
    Logical        :: is_first, is_last
    Integer(Int32) :: partner
    Call x_periodic_partner(is_first, is_last, partner)
    skip_x = ( x_bc_type == 0 .And. is_last )
    Call z_periodic_partner(is_first, is_last, partner)
    skip_z = ( z_bc_type == 0 .And. is_last )
  End Subroutine dup_cell_flags

  !> Accumulate IBM momentum exchange from one IBM application (pre/post state)
  Subroutine accumulate_ibm_stage_impulse(U_pre_, V_pre_, W_pre_, U_, V_, W_)

    Real(Int64), Dimension(nx, nyg,nzg), Intent(In) :: U_pre_, U_
    Real(Int64), Dimension(nxg, ny,nzg), Intent(In) :: V_pre_, V_
    Real(Int64), Dimension(nxg,nyg, nz), Intent(In) :: W_pre_, W_

    Integer(Int32) :: n, i, j, k
    Real   (Int64) :: dV
    Logical        :: skip_x, skip_z

    Call dup_cell_flags(skip_x, skip_z)

    Do n = 1, n_ghost_u
       i = ghost_u_idx(1,n);  j = ghost_u_idx(2,n);  k = ghost_u_idx(3,n)
       If ( skip_z .And. k == nzg-1 ) Cycle
       dV = (x(i+1)-x(i-1))*0.5d0 * (y(j)-y(j-1)) * (z(k)-z(k-1))
       ibm_Fx_acc = ibm_Fx_acc - (U_(i,j,k) - U_pre_(i,j,k)) * dV
    End Do

    Do n = 1, n_ghost_v
       i = ghost_v_idx(1,n);  j = ghost_v_idx(2,n);  k = ghost_v_idx(3,n)
       If ( skip_x .And. i == nxg-1 ) Cycle
       If ( skip_z .And. k == nzg-1 ) Cycle
       dV = (xg(i+1)-xg(i)) * (y(j+1)-y(j-1))*0.5d0 * (z(k)-z(k-1))
       ibm_Fy_acc = ibm_Fy_acc - (V_(i,j,k) - V_pre_(i,j,k)) * dV
    End Do

    Do n = 1, n_ghost_w
       i = ghost_w_idx(1,n);  j = ghost_w_idx(2,n);  k = ghost_w_idx(3,n)
       If ( skip_x .And. i == nxg-1 ) Cycle
       dV = (xg(i+1)-xg(i)) * (y(j)-y(j-1)) * (z(k+1)-z(k-1))*0.5d0
       ibm_Fz_acc = ibm_Fz_acc - (W_(i,j,k) - W_pre_(i,j,k)) * dV
    End Do

  End Subroutine accumulate_ibm_stage_impulse

  !> Build cell-centre ghost list for rigorous pressure integration
  Subroutine build_ghost_list_cc

    Integer(Int32) :: i, j, k, ng, nd
    Integer(Int32) :: ii, jj, kk
    Integer(Int32) :: iu, ju, ku, iv, jv, kv, iw, jw, kw
    Real   (Int64) :: nx_, ny_, nz_, dGB, dGI, xI, yI, zI

    ! Pass 1: count qualifying cells
    ng = 0;  nd = 0
    Do k = 2, nzg-1
       Do j = 2, nyg-1
          Do i = 2, nxg-1
             If ( phie(i,j,k) < 0d0 ) Then
                If ( phie(i-1,j,k) >= 0d0 .Or. phie(i+1,j,k) >= 0d0 .Or. &
                     phie(i,j-1,k) >= 0d0 .Or. phie(i,j+1,k) >= 0d0 .Or. &
                     phie(i,j,k-1) >= 0d0 .Or. phie(i,j,k+1) >= 0d0 ) Then
                   Call compute_normal_at_cc(i,j,k, nx_,ny_,nz_)
                   dGB = Abs(phie(i,j,k))
                   dGI = Max( 2d0*dGB, Real(n_image_layers,8)*dymin )
                   xI = xge(i) + dGI*nx_
                   yI = yg(j) + dGI*ny_
                   zI = zge(k) + dGI*nz_
                   Call find_stencil_centre(xI, yI, zI, ii, jj, kk)
                   If ( .Not. stencil_out .And. ii>=ext_cx_lo .And. ii<=ext_cx_hi .And. jj>=1 .And. jj<=nyg-1 .And. &
                        kk>=ext_cz_lo .And. kk<=ext_cz_hi ) Then
                      If ( Umaske(ii,jj,kk) > 0.5d0 ) Then
                         ng = ng + 1
                      Else
                         nd = nd + 1   ! image in solid — concave corner
                         ibm_drop_solid(4) = ibm_drop_solid(4) + 1
                      End If
                   Else
                      nd = nd + 1      ! image outside domain bounds
                      ibm_drop_outside(4) = ibm_drop_outside(4) + 1
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
    Allocate ( ghost_cc_dGI   (   n_ghost_cc) )
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
             If ( phie(i,j,k) < 0d0 ) Then
                If ( phie(i-1,j,k) >= 0d0 .Or. phie(i+1,j,k) >= 0d0 .Or. &
                     phie(i,j-1,k) >= 0d0 .Or. phie(i,j+1,k) >= 0d0 .Or. &
                     phie(i,j,k-1) >= 0d0 .Or. phie(i,j,k+1) >= 0d0 ) Then
                   Call compute_normal_at_cc(i,j,k, nx_,ny_,nz_)
                   dGB = Abs(phie(i,j,k))
                   dGI = Max( 2d0*dGB, Real(n_image_layers,8)*dymin )
                   xI = xge(i) + dGI*nx_
                   yI = yg(j) + dGI*ny_
                   zI = zge(k) + dGI*nz_
                   Call find_stencil_centre(xI, yI, zI, ii, jj, kk)
                   If ( .Not. stencil_out .And. ii>=ext_cx_lo .And. ii<=ext_cx_hi .And. jj>=1 .And. jj<=nyg-1 .And. &
                        kk>=ext_cz_lo .And. kk<=ext_cz_hi ) Then
                      If ( Umaske(ii,jj,kk) > 0.5d0 ) Then
                         ng = ng + 1
                         ghost_cc_idx(1,ng) = i
                         ghost_cc_idx(2,ng) = j
                         ghost_cc_idx(3,ng) = k
                         ghost_cc_objid(ng) = Min(Max(Nint(ibm_obj_id(i,j,k)), 0), max_ibm_objects)
                         ghost_cc_nrm(1,ng) = nx_
                         ghost_cc_nrm(2,ng) = ny_
                         ghost_cc_nrm(3,ng) = nz_
                         ghost_cc_dGB(ng)   = dGB
                         ghost_cc_dGI(ng)   = dGI
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
  !> SDF at the W (z-face) location (i,j,k), interpolated from the two adjacent cell centres with the local
  !  z weights (centres are face midpoints, so this is a plain average only on a uniform z grid); built from
  !  z/zg directly since setup_ibm runs before the interpolation weight arrays are allocated (host-only:
  !  the device W-zeroing loop below inlines the same expression, module arrays can't be used in an acc routine)
  Pure Function phi_w(i,j,k) Result(v)
    Integer(Int32), Intent(In) :: i, j, k
    Real   (Int64) :: v, w0
    w0 = ( zge(k+1) - ze(k) ) / ( zge(k+1) - zge(k) )
    v  = w0*phie(i,j,k) + ( 1d0 - w0 )*phie(i,j,k+1)
  End Function phi_w

  !> SDF at the V (y-face) location (i,j,k), interpolated from the two adjacent cell centres with the local
  !  y weights (plain average only on a uniform y grid); host-only, the device V-zeroing loop inlines the same expression
  Pure Function phi_v(i,j,k) Result(v)
    Integer(Int32), Intent(In) :: i, j, k
    Real   (Int64) :: v, w0
    w0 = ( yg(j+1) - y(j) ) / ( yg(j+1) - yg(j) )
    v  = w0*phie(i,j,k) + ( 1d0 - w0 )*phie(i,j+1,k)
  End Function phi_v

  Subroutine compute_normal_at_cc(i,j,k, nx_,ny_,nz_)
    Integer(Int32), Intent(In)  :: i, j, k
    Real   (Int64), Intent(Out) :: nx_, ny_, nz_
    Real   (Int64)              :: nmag, h_up, h_dn
    nx_  = (phie(i+1,j,k) - phie(i-1,j,k)) / (xge(i+1) - xge(i-1))
    ! Non-uniform central difference in y (2nd-order on stretched meshes)
    h_up = yg(j+1) - yg(j);  h_dn = yg(j) - yg(j-1)
    ny_  = ( h_dn**2*phie(i,j+1,k) + (h_up**2-h_dn**2)*phie(i,j,k) - h_up**2*phie(i,j-1,k) ) &
           / ( h_up * h_dn * (h_up + h_dn) )
    ! Non-uniform central difference in ze(2nd-order on stretched meshes)
    h_up = zge(k+1) - zge(k);  h_dn = zge(k) - zge(k-1)
    nz_  = ( h_dn**2*phie(i,j,k+1) + (h_up**2-h_dn**2)*phie(i,j,k) - h_up**2*phie(i,j,k-1) ) &
           / ( h_up * h_dn * (h_up + h_dn) )
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

    nx_  = ( phie(i+1,j,k) - phie(i,j,k) ) / ( xge(i+1) - xge(i) )
    ! Non-uniform central difference in y (2nd-order on stretched meshes)
    h_up = yg(j+1) - yg(j);  h_dn = yg(j) - yg(j-1)
    ny_  = ( h_dn**2*phie(i,j+1,k) + (h_up**2-h_dn**2)*phie(i,j,k) - h_up**2*phie(i,j-1,k) ) &
           / ( h_up * h_dn * (h_up + h_dn) )
    ! Non-uniform central difference in ze(2nd-order on stretched meshes)
    h_up = zge(k+1) - zge(k);  h_dn = zge(k) - zge(k-1)
    nz_  = ( h_dn**2*phie(i,j,k+1) + (h_up**2-h_dn**2)*phie(i,j,k) - h_up**2*phie(i,j,k-1) ) &
           / ( h_up * h_dn * (h_up + h_dn) )
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
    Real   (Int64)              :: nmag, h_up, h_dn
    nx_ = ( phie(i+1,j,k) - phie(i-1,j,k) ) / ( xge(i+1) - xge(i-1) )
    ny_ = ( phie(i,j+1,k) - phie(i,j,k) )   / ( yg(j+1) - yg(j) )
    ! Non-uniform central difference in ze(2nd-order on stretched meshes)
    h_up = zge(k+1) - zge(k);  h_dn = zge(k) - zge(k-1)
    nz_ = ( h_dn**2*phie(i,j,k+1) + (h_up**2-h_dn**2)*phie(i,j,k) - h_up**2*phie(i,j,k-1) ) &
           / ( h_up * h_dn * (h_up + h_dn) )
    nmag = Sqrt(nx_**2+ny_**2+nz_**2)
    If (nmag>1d-14) Then; nx_=nx_/nmag; ny_=ny_/nmag; nz_=nz_/nmag
    Else; nx_=0d0; ny_=1d0; nz_=0d0; End If
  End Subroutine compute_normal_at_face_v

  Subroutine compute_normal_at_face_w(i,j,k, nx_,ny_,nz_)
    Integer(Int32), Intent(In)  :: i, j, k
    Real   (Int64), Intent(Out) :: nx_, ny_, nz_
    Real   (Int64)              :: nmag, h_up, h_dn
    nx_  = ( phie(i+1,j,k) - phie(i-1,j,k) ) / ( xge(i+1) - xge(i-1) )
    ! Non-uniform central difference in y (2nd-order on stretched meshes)
    h_up = yg(j+1) - yg(j);  h_dn = yg(j) - yg(j-1)
    ny_  = ( h_dn**2*phie(i,j+1,k) + (h_up**2-h_dn**2)*phie(i,j,k) - h_up**2*phie(i,j-1,k) ) &
           / ( h_up * h_dn * (h_up + h_dn) )
    nz_  = ( phie(i,j,k+1) - phie(i,j,k) )   / ( zge(k+1) - zge(k) )
    nmag = Sqrt(nx_**2+ny_**2+nz_**2)
    If (nmag>1d-14) Then; nx_=nx_/nmag; ny_=ny_/nmag; nz_=nz_/nmag
    Else; nx_=0d0; ny_=1d0; nz_=0d0; End If
  End Subroutine compute_normal_at_face_w

  !> Set stencil_out when a point lies beyond the extended range on a side where the extended planes exist (data would be
  !  missing); on a non-periodic domain edge the legacy clamp to the edge stencil is kept
  Subroutine flag_out(p, lo_coord, hi_coord, has_lo, has_hi)
    Real   (Int64), Intent(In) :: p, lo_coord, hi_coord
    Logical,        Intent(In) :: has_lo, has_hi
    If ( has_lo .And. p <  lo_coord ) stencil_out = .True.
    If ( has_hi .And. p >= hi_coord ) stencil_out = .True.
  End Subroutine flag_out

  !> Extended-halo setup (host): depth ibm_E from the image distance, padded axes, SDF and fluid mask, anchor ranges, staging buffers
  Subroutine setup_ibm_ext

    Logical        :: is_first, is_last, per
    Integer(Int32) :: up, down, partner, Eloc, Emax, wmin, wmin_g, cnt
    Real   (Int64) :: dmin_h

    ! depth: image distance ~ max(2 dGB, n_image_layers*dymin) plus the trilinear corner and one cell of margin
    dmin_h = Min( x_global(2)-x_global(1), z_global(2)-z_global(1) )
    Eloc   = 2 + Ceiling( Real(n_image_layers,8)*dymin / dmin_h )
    If ( Eloc > 8 .And. myid == 0 ) Write(*,'(A,I0,A)') ' WARNING: IBM image-point halo depth ', Eloc, &
         ' exceeds the cap of 8 cells; farther image points will be dropped (reduce n_image_layers)'
    Eloc   = Min( Max(Eloc,2), 8 )
    Call MPI_Allreduce(Eloc, Emax, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    ibm_E = Emax
    wmin  = Min( nxg-2, nzg-2 )
    Call MPI_Allreduce(wmin, wmin_g, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
    If ( nprocs > 1 .And. wmin_g < ibm_E ) Then
       If ( myid == 0 ) Write(*,'(A,I0,A,I0,A)') ' ERROR: IBM needs at least ', ibm_E, ' interior cells per rank in x and z ', &
            'for its image-point stencils, but the thinnest slab has ', wmin_g, '. Use fewer ranks (or a different p_row/p_col).'
       Call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    End If
    If ( nprocs == 1 .And. wmin_g < ibm_E .And. myid == 0 ) Write(*,'(A,I0,A,I0,A)') &
         ' WARNING: IBM halo depth ', ibm_E, ' exceeds the domain interior (', wmin_g, ' cells)'

    Call x_halo_neighbors(up, down)
    Call x_periodic_partner(is_first, is_last, partner)
    per = ( x_bc_type == 0 )
    ibm_lo_x = ( down /= MPI_PROC_NULL ) .Or. per
    ibm_hi_x = ( up   /= MPI_PROC_NULL ) .Or. per
    Call z_halo_neighbors(up, down)
    per = ( z_bc_type == 0 )
    ibm_lo_z = ( down /= MPI_PROC_NULL ) .Or. per
    ibm_hi_z = ( up   /= MPI_PROC_NULL ) .Or. per

    ext_cx_lo = 1 - Merge(ibm_E, 0, ibm_lo_x);  ext_cx_hi = nxg - 1 + Merge(ibm_E, 0, ibm_hi_x)
    ext_fx_lo = Merge(1 - ibm_E, 2, ibm_lo_x);  ext_fx_hi = nx  - 1 + Merge(ibm_E, 0, ibm_hi_x)
    ext_cz_lo = 1 - Merge(ibm_E, 0, ibm_lo_z);  ext_cz_hi = nzg - 1 + Merge(ibm_E, 0, ibm_hi_z)
    ext_fz_lo = Merge(1 - ibm_E, 2, ibm_lo_z);  ext_fz_hi = nz  - 1 + Merge(ibm_E, 0, ibm_hi_z)

    Allocate( xe (1-ibm_E:nx +ibm_E), xge(1-ibm_E:nxg+ibm_E) )
    Allocate( ze (1-ibm_E:nz +ibm_E), zge(1-ibm_E:nzg+ibm_E) )
    Call pad_axis( x,  nx,  x_global,  nx_global,  i1_global(myid),  ibm_E, xe  )
    Call pad_axis( xg, nxg, xg_global, nxg_global, ig1_global(myid), ibm_E, xge )
    Call pad_axis( z,  nz,  z_global,  nz_global,  k1_global(myid),  ibm_E, ze  )
    Call pad_axis( zg, nzg, zg_global, nzg_global, kg1_global(myid), ibm_E, zge )

    Allocate( phie(1-ibm_E:nxg+ibm_E, nyg, 1-ibm_E:nzg+ibm_E), Umaske(1-ibm_E:nxg+ibm_E, nyg, 1-ibm_E:nzg+ibm_E) )
    Call pad_field( phi, nxg, nyg, nzg, .False., .False., ibm_E, phie )
    Where ( phie > 0d0 )
       Umaske = 1d0
    Elsewhere
       Umaske = 0d0
    End Where

    Allocate( Uext(1-ibm_E:nx +ibm_E, nyg, 1-ibm_E:nzg+ibm_E), Vext(1-ibm_E:nxg+ibm_E, ny,  1-ibm_E:nzg+ibm_E) )
    Allocate( Wext(1-ibm_E:nxg+ibm_E, nyg, 1-ibm_E:nz +ibm_E), Cext(1-ibm_E:nxg+ibm_E, nyg, 1-ibm_E:nzg+ibm_E) )
    Uext = 0d0;  Vext = 0d0;  Wext = 0d0;  Cext = 0d0
    cnt = Max( ibm_E*nyg*nzg, (nxg+2*ibm_E)*nyg*ibm_E ) + 8
    Allocate( ext_bs(cnt), ext_br(cnt) )
    ext_bs = 0d0;  ext_br = 0d0
    !$acc enter data create(Uext,Vext,Wext,Cext,ext_bs,ext_br)

    If ( myid == 0 ) Write(*,'(A,I0,A)') ' IBM: extended halo of ', ibm_E, ' planes for image-point stencils'

  End Subroutine setup_ibm_ext

  !  Find lower-left cell-centre index (ii,jj,kk) such that
  !  xg(ii) <= xp < xg(ii+1) etc. — trilinear stencil anchor.
  !  Search range is the full [1,Ng-1] the single MPI ghost/halo plane on each
  !  side actually supports (not [2,Ng-2]): an image point near a rank's own
  !  z-seam can legitimately anchor in that halo plane, and restricting the
  !  search to the interior silently mis-resolves it to the wrong index
  !  (previously observed as ghost cells dropped/misclassified differently
  !  depending on MPI rank count for cells near a z-decomposition seam).
  Subroutine find_stencil_centre(xp, yp, zp, ii, jj, kk)
    Real   (Int64), Intent(In)  :: xp, yp, zp
    Integer(Int32), Intent(Out) :: ii, jj, kk
    Integer(Int32) :: n
    stencil_out = .False.
    ii = ext_cx_lo
    Do n = ext_cx_lo, ext_cx_hi
       If ( xge(n) <= xp ) ii = n
    End Do
    jj = 1
    Do n = 1, nyg-1
       If ( yg(n) <= yp ) jj = n
    End Do
    kk = ext_cz_lo
    Do n = ext_cz_lo, ext_cz_hi
       If ( zge(n) <= zp ) kk = n
    End Do
    Call flag_out(xp, xge(ext_cx_lo), xge(ext_cx_hi+1), ibm_lo_x, ibm_hi_x)
    Call flag_out(zp, zge(ext_cz_lo), zge(ext_cz_hi+1), ibm_lo_z, ibm_hi_z)
  End Subroutine find_stencil_centre

  ! Staggered stencil finders, one per velocity component: U uses x-faces x cell-centres, V cell-centres x y-faces, W cell-centres x z-faces
  ! Own-direction (x/y/z) loops below are left at [2,N-1]: those arrays aren't
  ! MPI z-halo-limited to 1 plane the way xg/yg/zg are. The xg/yg/zg loops use
  ! the full [1,Ng-1] range for the same reason as find_stencil_centre above.
  Subroutine find_stencil_u(xp, yp, zp, ii, jj, kk)
    Real   (Int64), Intent(In)  :: xp, yp, zp
    Integer(Int32), Intent(Out) :: ii, jj, kk
    Integer(Int32) :: n
    stencil_out = .False.
    ii = ext_fx_lo
    Do n = ext_fx_lo, ext_fx_hi
       If ( xe(n) <= xp ) ii = n
    End Do
    jj = 1
    Do n = 1, nyg-1
       If ( yg(n) <= yp ) jj = n
    End Do
    kk = ext_cz_lo
    Do n = ext_cz_lo, ext_cz_hi
       If ( zge(n) <= zp ) kk = n
    End Do
    Call flag_out(xp, xe(ext_fx_lo), xe(ext_fx_hi+1), ibm_lo_x, ibm_hi_x)
    Call flag_out(zp, zge(ext_cz_lo), zge(ext_cz_hi+1), ibm_lo_z, ibm_hi_z)
  End Subroutine find_stencil_u

  Subroutine find_stencil_v(xp, yp, zp, ii, jj, kk)
    Real   (Int64), Intent(In)  :: xp, yp, zp
    Integer(Int32), Intent(Out) :: ii, jj, kk
    Integer(Int32) :: n
    stencil_out = .False.
    ii = ext_cx_lo
    Do n = ext_cx_lo, ext_cx_hi
       If ( xge(n) <= xp ) ii = n
    End Do
    jj = 2
    Do n = 2, ny-1
       If ( y(n) <= yp ) jj = n
    End Do
    kk = ext_cz_lo
    Do n = ext_cz_lo, ext_cz_hi
       If ( zge(n) <= zp ) kk = n
    End Do
    Call flag_out(xp, xge(ext_cx_lo), xge(ext_cx_hi+1), ibm_lo_x, ibm_hi_x)
    Call flag_out(zp, zge(ext_cz_lo), zge(ext_cz_hi+1), ibm_lo_z, ibm_hi_z)
  End Subroutine find_stencil_v

  Subroutine find_stencil_w(xp, yp, zp, ii, jj, kk)
    Real   (Int64), Intent(In)  :: xp, yp, zp
    Integer(Int32), Intent(Out) :: ii, jj, kk
    Integer(Int32) :: n
    stencil_out = .False.
    ii = ext_cx_lo
    Do n = ext_cx_lo, ext_cx_hi
       If ( xge(n) <= xp ) ii = n
    End Do
    jj = 1
    Do n = 1, nyg-1
       If ( yg(n) <= yp ) jj = n
    End Do
    kk = ext_fz_lo
    Do n = ext_fz_lo, ext_fz_hi
       If ( ze(n) <= zp ) kk = n
    End Do
    Call flag_out(xp, xge(ext_cx_lo), xge(ext_cx_hi+1), ibm_lo_x, ibm_hi_x)
    Call flag_out(zp, ze(ext_fz_lo), ze(ext_fz_hi+1), ibm_lo_z, ibm_hi_z)
  End Subroutine find_stencil_w

  !> Compute 8 trilinear weights (unit-cube corner order) for point (xp,yp,zp) anchored at centre index (ii,jj,kk)
  Subroutine trilinear_weights(xp, yp, zp, ii, jj, kk, w)
    Real   (Int64), Intent(In)  :: xp, yp, zp
    Integer(Int32), Intent(In)  :: ii, jj, kk
    Real   (Int64), Intent(Out) :: w(8)
    Real   (Int64) :: tx, ty, tz

    tx = (xp - xge(ii)) / Max(xge(ii+1)-xge(ii), 1d-14)
    ty = (yp - yg(jj)) / Max(yg(jj+1)-yg(jj), 1d-14)
    tz = (zp - zge(kk)) / Max(zge(kk+1)-zge(kk), 1d-14)

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
    tx = (xp - xe(ii))  / Max(xe(ii+1)  - xe(ii),  1d-14)
    ty = (yp - yg(jj)) / Max(yg(jj+1) - yg(jj), 1d-14)
    tz = (zp - zge(kk)) / Max(zge(kk+1) - zge(kk), 1d-14)
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
    tx = (xp - xge(ii)) / Max(xge(ii+1) - xge(ii), 1d-14)
    ty = (yp - y(jj))  / Max(y(jj+1)  - y(jj),  1d-14)
    tz = (zp - zge(kk)) / Max(zge(kk+1) - zge(kk), 1d-14)
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
    tx = (xp - xge(ii)) / Max(xge(ii+1) - xge(ii), 1d-14)
    ty = (yp - yg(jj)) / Max(yg(jj+1) - yg(jj), 1d-14)
    tz = (zp - ze(kk))  / Max(ze(kk+1)  - ze(kk),  1d-14)
    tx = Max(0d0, Min(1d0, tx));  ty = Max(0d0, Min(1d0, ty));  tz = Max(0d0, Min(1d0, tz))
    w(1) = (1d0-tx)*(1d0-ty)*(1d0-tz);  w(2) =      tx *(1d0-ty)*(1d0-tz)
    w(3) = (1d0-tx)*     ty *(1d0-tz);  w(4) =      tx *     ty *(1d0-tz)
    w(5) = (1d0-tx)*(1d0-ty)*     tz ;  w(6) =      tx *(1d0-ty)*     tz
    w(7) = (1d0-tx)*     ty *     tz ;  w(8) =      tx *     ty *     tz
  End Subroutine trilinear_weights_w

  !> Trilinear interpolation at image point using stencil anchor (img) and weights (w) from setup_ibm, on the component's staggered grid
  Real(Int64) Function trilinear_interp_u(U_, w, img, e)
    !$acc routine seq
    Integer(Int32), Intent(In) :: e   ! extended-halo depth: the array is dimensioned (1-e:n+e, :, 1-e:n+e)
    Real   (Int64), Dimension(1-e:,:,1-e:), Intent(In) :: U_
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

  Real(Int64) Function trilinear_interp_v(V_, w, img, e)
    !$acc routine seq
    Integer(Int32), Intent(In) :: e   ! extended-halo depth: the array is dimensioned (1-e:n+e, :, 1-e:n+e)
    Real   (Int64), Dimension(1-e:,:,1-e:), Intent(In) :: V_
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

  Real(Int64) Function trilinear_interp_w(W_, w, img, e)
    !$acc routine seq
    Integer(Int32), Intent(In) :: e   ! extended-halo depth: the array is dimensioned (1-e:n+e, :, 1-e:n+e)
    Real   (Int64), Dimension(1-e:,:,1-e:), Intent(In) :: W_
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
  Real(Int64) Function trilinear_interp_p(P_, w, img, e)
    !$acc routine seq
    Integer(Int32), Intent(In) :: e   ! extended-halo depth: the array is dimensioned (1-e:n+e, :, 1-e:n+e)
    Real   (Int64), Dimension(1-e:,:,1-e:), Intent(In) :: P_
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
    Real   (Int64) :: dV, dA, U_I, p_I, nu_t_B, dGB, dIB
    Real   (Int64) :: lFx_pres, lFy_pres, lFz_pres
    Real   (Int64) :: lFx_visc, lFy_visc, lFz_visc
    Logical        :: skip_x, skip_z

    Real   (Int64), Allocatable :: Ue(:,:,:), Ve(:,:,:), We(:,:,:), Pe(:,:,:)

    ! host copies with the extended halo, for the image-point stencils
    Allocate( Ue(1-ibm_E:nx +ibm_E, nyg, 1-ibm_E:nzg+ibm_E), Ve(1-ibm_E:nxg+ibm_E, ny,  1-ibm_E:nzg+ibm_E), &
              We(1-ibm_E:nxg+ibm_E, nyg, 1-ibm_E:nz +ibm_E), Pe(1-ibm_E:nxg+ibm_E, nyg, 1-ibm_E:nzg+ibm_E) )
    Call pad_field( U_, nx,  nyg, nzg, .True.,  .False., ibm_E, Ue )
    Call pad_field( V_, nxg, ny,  nzg, .False., .False., ibm_E, Ve )
    Call pad_field( W_, nxg, nyg, nz,  .False., .True.,  ibm_E, We )
    Call pad_field( P,  nxg, nyg, nzg, .False., .False., ibm_E, Pe )

    Call dup_cell_flags(skip_x, skip_z)
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
       If ( skip_x .And. i == nxg-1 ) Cycle
       If ( skip_z .And. k == nzg-1 ) Cycle
       dGB = ghost_cc_dGB(n)
       dV  = (xg(i+1)-xg(i-1))*0.5d0 * (yg(j+1)-yg(j-1))*0.5d0 * (zg(k+1)-zg(k-1))*0.5d0
       dA  = dV / Max(dGB, 1d-14)
       p_I = trilinear_interp_p(Pe, ghost_cc_wgt_cc(1:8,n), ghost_cc_img_cc(:,n), ibm_E)
       lFx_pres = lFx_pres - p_I * ghost_cc_nrm(1,n) * dA
       lFy_pres = lFy_pres - p_I * ghost_cc_nrm(2,n) * dA
       lFz_pres = lFz_pres - p_I * ghost_cc_nrm(3,n) * dA
    End Do

    !--- Method 2 viscous: staggered ghost lists, one component each ---
    ! dA is the perpendicular face area; using dV/dGB instead would give a 1/dGB^2 singularity since the gradient already carries 1/dGB
    ! dIB = actual boundary-to-image distance (dGI-dGB): equals dGB only when dGI is the
    ! unclamped mirror 2*dGB; the gradient must divide by the real B-to-I distance (see ghost_*_dGI)
    Do n = 1, n_ghost_u
       i = ghost_u_idx(1,n);  j = ghost_u_idx(2,n);  k = ghost_u_idx(3,n)
       If ( skip_z .And. k == nzg-1 ) Cycle
       dGB = ghost_u_dGB(n)
       dIB = ghost_u_dGI(n) - dGB
       dA  = (y(j)-y(j-1)) * (z(k)-z(k-1))             ! y-z face area
       U_I = trilinear_interp_u(Ue, ghost_u_wgt(1:8,n), ghost_u_img(:,n), ibm_E)
       nu_t_B = 0.5d0*(nu_t(i,j,k) + nu_t(Min(i+1,nxg),j,k))
       lFx_visc = lFx_visc + (nu + nu_t_B) * (U_I - U_wall) / Max(dIB, 1d-14) * dA
    End Do

    Do n = 1, n_ghost_v
       i = ghost_v_idx(1,n);  j = ghost_v_idx(2,n);  k = ghost_v_idx(3,n)
       If ( skip_x .And. i == nxg-1 ) Cycle
       If ( skip_z .And. k == nzg-1 ) Cycle
       dGB = ghost_v_dGB(n)
       dIB = ghost_v_dGI(n) - dGB
       dA  = (xg(i+1)-xg(i)) * (z(k)-z(k-1))           ! x-z face area
       U_I = trilinear_interp_v(Ve, ghost_v_wgt(1:8,n), ghost_v_img(:,n), ibm_E)
       nu_t_B = 0.5d0*(nu_t(i,j,k) + nu_t(i,Min(j+1,nyg),k))
       lFy_visc = lFy_visc + (nu + nu_t_B) * (U_I - V_wall) / Max(dIB, 1d-14) * dA
    End Do

    Do n = 1, n_ghost_w
       i = ghost_w_idx(1,n);  j = ghost_w_idx(2,n);  k = ghost_w_idx(3,n)
       If ( skip_x .And. i == nxg-1 ) Cycle
       dGB = ghost_w_dGB(n)
       dIB = ghost_w_dGI(n) - dGB
       dA  = (xg(i+1)-xg(i)) * (y(j)-y(j-1))           ! x-y face area
       U_I = trilinear_interp_w(We, ghost_w_wgt(1:8,n), ghost_w_img(:,n), ibm_E)
       nu_t_B = 0.5d0*(nu_t(i,j,k) + nu_t(i,j,Min(k+1,nzg)))
       lFz_visc = lFz_visc + (nu + nu_t_B) * (U_I - W_wall) / Max(dIB, 1d-14) * dA
    End Do

    Deallocate( Ue, Ve, We, Pe )

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
    Real   (Int64) :: dV, dA, dGB, invd, invdI, pB, nutB, uI, vI, wI
    Real   (Int64), Allocatable :: lbuf(:,:), gbuf(:,:)
    Integer(Int32), Allocatable :: counts(:), recvc(:), displs(:)
    Character(256) :: fname
    Character(16)  :: ext
    Logical        :: dir_exists

    Real   (Int64), Allocatable :: Ue(:,:,:), Ve(:,:,:), We(:,:,:), Pe(:,:,:), Ne(:,:,:)

    nl = n_ghost_cc

    ! host copies with the extended halo, for the image-point stencils
    Allocate( Ue(1-ibm_E:nx +ibm_E, nyg, 1-ibm_E:nzg+ibm_E), Ve(1-ibm_E:nxg+ibm_E, ny,  1-ibm_E:nzg+ibm_E), &
              We(1-ibm_E:nxg+ibm_E, nyg, 1-ibm_E:nz +ibm_E), Pe(1-ibm_E:nxg+ibm_E, nyg, 1-ibm_E:nzg+ibm_E), &
              Ne(1-ibm_E:nxg+ibm_E, nyg, 1-ibm_E:nzg+ibm_E) )
    Call pad_field( U_,   nx,  nyg, nzg, .True.,  .False., ibm_E, Ue )
    Call pad_field( V_,   nxg, ny,  nzg, .False., .False., ibm_E, Ve )
    Call pad_field( W_,   nxg, nyg, nz,  .False., .True.,  ibm_E, We )
    Call pad_field( P,    nxg, nyg, nzg, .False., .False., ibm_E, Pe )
    Call pad_field( nu_t, nxg, nyg, nzg, .False., .False., ibm_E, Ne )

    ! Fill the rank-local per-point buffer; row layout (NF=13): 1-3 x,y,z 4-6 nx,ny,nz 7 p 8-10 fp_x,fp_y,fp_z 11-13 fv_x,fv_y,fv_z
    Allocate ( lbuf(NF, Max(nl,1)) )
    Do n = 1, nl
       i = ghost_cc_idx(1,n);  j = ghost_cc_idx(2,n);  k = ghost_cc_idx(3,n)
       dGB  = ghost_cc_dGB(n)
       invd = 1d0 / Max(dGB, 1d-14)
       ! Distance from boundary B to image point I: dGB when I is at the unclamped mirror
       ! location (dGI=2*dGB), but (dGI-dGB) when dGI was clamped up (see ghost_cc_dGI) --
       ! the velocity gradient below must divide by the actual B-to-I distance, not dGB.
       invdI = 1d0 / Max(ghost_cc_dGI(n) - dGB, 1d-14)
       dV   = (xg(i+1)-xg(i-1))*0.5d0 * (yg(j+1)-yg(j-1))*0.5d0 * (zg(k+1)-zg(k-1))*0.5d0
       dA   = dV * invd

       ! interpolate pressure, turbulent viscosity and velocity at image point I
       pB   = trilinear_interp_p(Pe,   ghost_cc_wgt_cc(1:8,n), ghost_cc_img_cc(:,n), ibm_E)
       nutB = trilinear_interp_p(Ne,   ghost_cc_wgt_cc(1:8,n), ghost_cc_img_cc(:,n), ibm_E)
       uI   = trilinear_interp_u(Ue,   ghost_cc_wgt_u (1:8,n), ghost_cc_img_u (:,n), ibm_E)
       vI   = trilinear_interp_v(Ve,   ghost_cc_wgt_v (1:8,n), ghost_cc_img_v (:,n), ibm_E)
       wI   = trilinear_interp_w(We,   ghost_cc_wgt_w (1:8,n), ghost_cc_img_w (:,n), ibm_E)

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
       lbuf(11,n) = (nu + nutB) * (uI - U_wall) * invdI * dA
       lbuf(12,n) = (nu + nutB) * (vI - V_wall) * invdI * dA
       lbuf(13,n) = (nu + nutB) * (wI - W_wall) * invdI * dA
    End Do

    Deallocate( Ue, Ve, We, Pe, Ne )

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

