!   Module for boundary conditions
Module boundary_conditions

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use mpi
  Use decomp, Only : z_halo_neighbors, x_halo_neighbors, z_periodic_partner, x_periodic_partner
  Use synthetic_eddy_method

  ! prevent implicit typing
  Implicit None

Contains

  ! Apply boundary conditions to velocities in the 3 directions; after_projection semantics
  Subroutine apply_boundary_conditions(after_projection)

    Logical, Intent(In), Optional :: after_projection
    Logical :: skip_outflow_u
    Real   (Int64) :: Uc

    skip_outflow_u = .False.
    If ( Present(after_projection) ) skip_outflow_u = after_projection

    ! interior region (z-direction)
    Call update_ghost_interior_planes(U,1)
    Call update_ghost_interior_planes(V,2)
    Call update_ghost_interior_planes(W,3)

    ! interior region (x-direction: needed even under x_bc_type==1, for
    ! interior row boundaries between the domain's x=1/x=nx_global edges)
    Call update_ghost_interior_planes_x(U,1)
    Call update_ghost_interior_planes_x(V,2)
    Call update_ghost_interior_planes_x(W,2)

    ! x direction: periodic, or Dirichlet-SEM inflow / convective outflow
    If ( x_bc_type == 0 ) Then
       Call apply_periodic_bc_x(U,1)
       Call apply_periodic_bc_x(V,2)
       Call apply_periodic_bc_x(W,2)
    Else
       Call update_inflow_recycle(t)   ! no-op unless inflow_type==2 (recycled precursor inflow); advances the bracketing donor frames to the current stage time
       Call apply_inflow_bc_x (U,1)
       Call apply_inflow_bc_x (V,2)
       Call apply_inflow_bc_x (W,3)
       Uc = outflow_convection_velocity()
       If ( .Not. skip_outflow_u ) Call apply_outflow_bc_x(U,Uc)
       ! V/W's convective relaxation isn't idempotent (unlike the zero-gradient BC it replaced), so apply once per RK stage, on the post-projection refresh, not twice
       If ( skip_outflow_u ) Call apply_outflow_bc_x(V,Uc)
       If ( skip_outflow_u ) Call apply_outflow_bc_x(W,Uc)
    End If

    ! apply periodicity in z 
    Call apply_periodic_bc_z(U,1)
    Call apply_periodic_bc_z(V,2)
    Call apply_periodic_bc_z(W,3)

    ! y direction: periodic, or Robin (wall/slip-length)
    If ( y_bc_type == 0 ) Then
       Call apply_periodic_bc_y(U,2)
       Call apply_periodic_bc_y(V,1)
       Call apply_periodic_bc_y(W,2)
    Else
       ! U boundary condition at the wall
       Call apply_Robin_bc_y(U,alpha_x,2)

       ! V boundary condition at the wall
       Call apply_Robin_bc_y(V,alpha_y,1)

       ! W boundary condition at the wall
       Call apply_Robin_bc_y(W,alpha_z,2)
    End If

    ! compute boundary conditions for pseudo-pressure
  End Subroutine apply_boundary_conditions
  
  ! Periodicity in x, MPI communication required; id: 1=x-faces, 2=x-centres; F is in/out
  Subroutine apply_periodic_bc_x(F,id)

    Real   (Int64), Intent(InOut) :: F(:,:,:)
    Integer(Int32), Intent(In)    :: id
    Logical :: is_first, is_last
    Integer(Int32) :: partner
    Integer(Int32) :: n2, n3
    Real(Int64), Allocatable :: b1(:,:), be(:,:)
    Real(Int64), Allocatable :: bpair(:,:,:)

    Call x_periodic_partner(is_first, is_last, partner)

    ! p_row==1: this rank self-pairs (its column has no other row), so wrap
    ! around via direct local copy instead of Mpi_sendrecv
    If ( is_first .And. is_last ) Then
       !$acc kernels present(F)
       If ( id==1 ) Then
          F( 1,:,:) = F(nx-1,:,:)
          F(nx,:,:) = F(   2,:,:)
       Else
          F(    1,:,:) = F(nxg-2,:,:)
          F(nxg-1,:,:) = F(    2,:,:) ! see note*
          F(nxg  ,:,:) = F(    3,:,:)
       End If
       !$acc end kernels
       Return
    End If

    n2 = Size(F,2)
    n3 = Size(F,3)

    If ( id == 1 ) Then
       ! F at x faces: exactly one ghost plane needed on each side, one message each way
       Allocate ( b1(n2,n3), be(n2,n3) )
       If ( is_first ) Then
          b1 = F(2,:,:)
          Call Mpi_sendrecv(b1, n2*n3, Mpi_real8, partner, 10, &
               be, n2*n3, Mpi_real8, partner, 10, MPI_COMM_WORLD, istat, ierr)
          F(1,:,:) = be
       Else If ( is_last ) Then
          be = F(nx-1,:,:)
          Call Mpi_sendrecv(be, n2*n3, Mpi_real8, partner, 10, &
               b1, n2*n3, Mpi_real8, partner, 10, MPI_COMM_WORLD, istat, ierr)
          F(nx,:,:) = b1
       End If
       Deallocate(b1,be)
    Else
       ! F at x centres: is_first needs 1 ghost plane but must supply 2 (is_last needs both);
       ! send/receive the 2-plane side as a single combined message instead of two sequential ones
       Allocate ( bpair(n2,n3,2), be(n2,n3) )
       If ( is_first ) Then
          bpair(:,:,1) = F(2,:,:)
          bpair(:,:,2) = F(3,:,:)
          Call Mpi_sendrecv(bpair, n2*n3*2, Mpi_real8, partner, 11, &
               be, n2*n3, Mpi_real8, partner, 12, MPI_COMM_WORLD, istat, ierr)
          F(1,:,:) = be
       Else If ( is_last ) Then
          be = F(nxg-2,:,:)
          Call Mpi_sendrecv(be, n2*n3, Mpi_real8, partner, 12, &
               bpair, n2*n3*2, Mpi_real8, partner, 11, MPI_COMM_WORLD, istat, ierr)
          F(nxg-1,:,:) = bpair(:,:,1) ! see note*
          F(nxg  ,:,:) = bpair(:,:,2)
       End If
       Deallocate(bpair,be)
    End If

    ! Note: handles non-periodic ICs; redundant after the first time step.
  End Subroutine apply_periodic_bc_x

  ! Dirichlet inflow (x_bc_type==1): F face/ghost = mean profile + SEM fluctuation (sem.f90); comp: 1=U,2=V,3=W; F in/out; no-op except on the rank owning the x=1 boundary
  Subroutine apply_inflow_bc_x(F,comp)

    Real   (Int64), Intent(InOut) :: F(:,:,:)
    Integer(Int32), Intent(In)    :: comp

    Integer(Int32) :: j, k, n2, n3
    Real   (Int64) :: up, vp, wp, target_val
    Logical :: is_first, is_last
    Integer(Int32) :: partner

    Call x_periodic_partner(is_first, is_last, partner)
    If ( .Not. is_first ) Return

    n2 = Size(F,2)
    n3 = Size(F,3)

    ! F's y/z coords depend on staggered-grid placement: U=(yc,zc), V=(yf,zc), W=(yc,zf); every (j,k) is independent -- sem_fluctuation/mean_profile_U/recycle_value are !$acc routine seq
    !$acc parallel loop collapse(2) present(F,y,yg,z,zg) private(up,vp,wp,target_val)
    Do k = 1, n3
       Do j = 1, n2
          If ( inflow_type == 2 ) Then
             ! recycled precursor inflow: donor data is already the full field (mean+fluctuation), no mean_profile/SEM decomposition needed
             If ( comp == 1 ) Then
                F(1,j,k) = recycle_value(1, j, k)               ! Dirichlet face value, exact
             Else If ( comp == 2 ) Then
                target_val = recycle_value(2, j, k)
                F(1,j,k) = 2d0*target_val - F(2,j,k)             ! ghost-cell mirror
             Else
                target_val = recycle_value(3, j, k)
                F(1,j,k) = 2d0*target_val - F(2,j,k)             ! ghost-cell mirror
             End If
          Else If ( comp == 1 ) Then
             Call sem_fluctuation( 1, j, k, yg(j), zg(k), t, up, vp, wp )
             target_val = mean_profile_U( yg(j) ) + up
             F(1,j,k) = target_val                  ! Dirichlet face value, exact
          Else If ( comp == 2 ) Then
             Call sem_fluctuation( 2, j, k, y(j), zg(k), t, up, vp, wp )
             target_val = vp
             F(1,j,k) = 2d0*target_val - F(2,j,k)    ! ghost-cell mirror
          Else
             Call sem_fluctuation( 3, j, k, yg(j), z(k), t, up, vp, wp )
             target_val = wp
             F(1,j,k) = 2d0*target_val - F(2,j,k)    ! ghost-cell mirror
          End If
       End Do
    End Do
    !$acc end parallel loop

  End Subroutine apply_inflow_bc_x

  ! Dirichlet inflow for the cell-centred temperature field (x_bc_type==1, host-only: apply_temperature_bc runs host-side, see thermal_transport.f90); no fluctuation coupling; no-op except on the rank owning the x=1 boundary
  Subroutine apply_inflow_bc_scalar_x(T_)

    Real(Int64), Intent(InOut) :: T_(:,:,:)

    Integer(Int32) :: j, k, n2, n3
    Logical :: is_first, is_last
    Integer(Int32) :: partner

    Call x_periodic_partner(is_first, is_last, partner)
    If ( .Not. is_first ) Return

    n2 = Size(T_,2)
    n3 = Size(T_,3)

    Call update_inflow_recycle(t)   ! no-op unless inflow_type==2; self-contained here in case this runs before apply_boundary_conditions in some call ordering

    If ( inflow_type == 2 ) Then
       Do k = 1, n3
          Do j = 1, n2
             T_(1,j,k) = recycle_value(4, j, k)
          End Do
       End Do
    Else
       Do k = 1, n3
          Do j = 1, n2
             T_(1,j,k) = mean_profile_T( yg(j) )
          End Do
       End Do
    End If

  End Subroutine apply_inflow_bc_scalar_x

  ! Dirichlet inflow for the cell-centred sediment concentration field (x_bc_type==1, sediment_flag==1; host-only, mirrors apply_inflow_bc_scalar_x); no mean-profile file support (unlike T) -- constant C_ref fallback, or recycled donor value under inflow_type==2; no-op except on the rank owning the x=1 boundary
  Subroutine apply_inflow_bc_scalar_x_C(C_)

    Real(Int64), Intent(InOut) :: C_(:,:,:)

    Integer(Int32) :: j, k, n2, n3
    Logical :: is_first, is_last
    Integer(Int32) :: partner

    Call x_periodic_partner(is_first, is_last, partner)
    If ( .Not. is_first ) Return

    n2 = Size(C_,2)
    n3 = Size(C_,3)

    Call update_inflow_recycle(t)   ! no-op unless inflow_type==2; self-contained here in case this runs before apply_boundary_conditions in some call ordering

    If ( inflow_type == 2 ) Then
       Do k = 1, n3
          Do j = 1, n2
             C_(1,j,k) = recycle_value(5, j, k)
          End Do
       End Do
    Else
       C_(1,:,:) = C_ref
    End If

  End Subroutine apply_inflow_bc_scalar_x_C

  ! Plane-averaged streamwise velocity at the outlet (x_bc_type==1, no MPI reduction: local-rank average, consistent with the other x_bc_type==1 routines' no-MPI scope), used as the convection speed for apply_outflow_bc_x
  Function outflow_convection_velocity() Result(Uc)

    Real(Int64) :: Uc
    Integer(Int32) :: n2, n3

    n2 = Size(U,2)
    n3 = Size(U,3)

    !$acc kernels present(U)
    Uc = Sum(U(nx-1,2:n2-1,2:n3-1)) / Real((n2-2)*(n3-2),Int64)
    !$acc end kernels

  End Function outflow_convection_velocity

  ! Convective outflow (x_bc_type==1): dF/dt + Uc*dF/dx = 0, explicit upwind from interior; pairs with homogeneous Dirichlet pressure BC in the DCT-IV solve; F (U,V, or W) in/out; no-op except on the rank owning the x=nx_global boundary
  Subroutine apply_outflow_bc_x(F,Uc)

    Real(Int64), Intent(InOut) :: F(:,:,:)
    Real(Int64), Intent(In)    :: Uc
    Integer(Int32) :: nlast
    Real   (Int64) :: courant
    Logical :: is_first, is_last
    Integer(Int32) :: partner

    Call x_periodic_partner(is_first, is_last, partner)
    If ( .Not. is_last ) Return

    nlast = Size(F,1)
    ! clip to [0,1]: negative Uc (local backflow) would pull the outlet value from outside the domain, and Courant>1 is unconditionally unstable
    courant = Min(Max(Uc,0d0)*dt/dx,1d0)

    !$acc kernels present(F)
    F(nlast,:,:) = F(nlast,:,:) - courant*( F(nlast,:,:) - F(nlast-1,:,:) )
    !$acc end kernels

  End Subroutine apply_outflow_bc_x

  ! Periodicity in z, MPI communication required
  Subroutine apply_periodic_bc_z(F,id)

    Real   (Int64), Intent(InOut) :: F(:,:,:)
    Integer(Int32), Intent(In)    :: id
    Logical :: is_first, is_last
    Integer(Int32) :: partner

    Call z_periodic_partner(is_first, is_last, partner)

    ! p_col==1: this rank self-pairs (its row has no other column), so wrap
    ! around via direct local copy instead of Mpi_sendrecv (as in apply_periodic_bc_x)
    If ( is_first .And. is_last ) Then
       !$acc kernels present(F)
       If ( id == 3 ) Then
          F(:,:,1)  = F(:,:,nz-1)
          F(:,:,nz) = F(:,:,2)
       Else
          F(:,:,1)     = F(:,:,nzg-2)
          F(:,:,nzg-1) = F(:,:,2)
          F(:,:,nzg)   = F(:,:,3)
       End If
       !$acc end kernels
       Return
    End If

    ! save planes
    If ( is_first ) Then
      ! begin planes
      If (id == 3) Then
        ! F defined at z faces
        buffer_wi(:,:)   = F(:,:,2)
      Elseif (id == 1) Then 
        ! F defined at z centers
         buffer_ui(:,:,2) = F(:,:,2) 
         buffer_ui(:,:,3) = F(:,:,3)
      Elseif (id == 2) Then
        ! F defined at z centers
         buffer_vi(:,:,2) = F(:,:,2) 
         buffer_vi(:,:,3) = F(:,:,3)
      Elseif (id == 4) Then
        ! F defined at z centers
         buffer_ci(:,:,2) = F(:,:,2) 
         buffer_ci(:,:,3) = F(:,:,3)
      End If      

    End If

    If ( is_last ) Then
       ! end planes
      If (id == 3) Then
        ! F defined at z faces
        buffer_we(:,:) = F(:,:, nz-1)
      Elseif (id == 1) Then
        ! F defined at z centers
        buffer_ue(:,:) = F(:,:,nzg-2)
      Elseif (id == 2) Then
        ! F defined at z centers
        buffer_ve(:,:) = F(:,:,nzg-2)
      Elseif (id == 4) Then
        ! F defined at z centers
        buffer_ce(:,:) = F(:,:,nzg-2)
      End If
    End If

    ! communicate planes
    If ( is_first ) Then
      If (id == 3) Then
        ! Send/receive W
        Call Mpi_sendrecv(buffer_wi, nxg*nyg, Mpi_real8, partner, 3,  &
        buffer_we, nxg*nyg, Mpi_real8, partner, 3, MPI_COMM_WORLD,    &
        istat, ierr)
      Elseif (id == 1) Then
        ! Send/receive U
        Call Mpi_sendrecv(buffer_ui, nx*nyg*2, Mpi_real8, partner, 1, &
        buffer_ue, nx*nyg, Mpi_real8, partner, 1, MPI_COMM_WORLD,     &
        istat, ierr)
      Elseif (id == 2) Then
        ! Send/receive V
        Call Mpi_sendrecv(buffer_vi, nxg*ny*2, Mpi_real8, partner, 2, &
        buffer_ve, nxg*ny, Mpi_real8, partner, 2, MPI_COMM_WORLD,     &
        istat, ierr)
      Elseif (id == 4) Then
        ! Send/receive C
        Call Mpi_sendrecv(buffer_ci, nxg*nyg*2, Mpi_real8, partner, 4, &
        buffer_ce, nxg*nyg, Mpi_real8, partner, 4, MPI_COMM_WORLD,     &
        istat, ierr)
      End If
    End If

    If ( is_last ) Then
      If (id == 3) Then
        ! Send/receive W
        Call Mpi_sendrecv(buffer_we, nxg*nyg, Mpi_real8, partner, 3, &
        buffer_wi, nxg*nyg, Mpi_real8, partner, 3, MPI_COMM_WORLD,   &
        istat, ierr)
      Elseif (id == 1) Then
        ! Send/receive U
        Call Mpi_sendrecv(buffer_ue, nx*nyg, Mpi_real8, partner, 1,  &
        buffer_ui, nx*nyg*2, Mpi_real8, partner, 1, MPI_COMM_WORLD,  &
        istat, ierr)
      Elseif (id == 2) Then
        ! Send/receive V
        Call Mpi_sendrecv(buffer_ve, nxg*ny, Mpi_real8, partner, 2,  &
        buffer_vi, nxg*ny*2, Mpi_real8, partner, 2, MPI_COMM_WORLD,  &
        istat, ierr)
      Elseif (id == 4) Then
        ! Send/receive C
        Call Mpi_sendrecv(buffer_ce, nxg*nyg, Mpi_real8, partner, 4,  &
        buffer_ci, nxg*nyg*2, Mpi_real8, partner, 4, MPI_COMM_WORLD,  &
        istat, ierr)
      End If
    End If

    ! apply conditions
    If ( is_first ) Then
      If (id == 3) Then
        F(:,:,1) = buffer_we(:,:)       ! W_global(:,:,nz_global-1)
      Elseif (id == 1) Then
        F(:,:,1) = buffer_ue(:,:)       ! U_global(:,:,nzg_global-2)
      Elseif (id == 2) Then
        F(:,:,1) = buffer_ve(:,:)       ! U_global(:,:,nzg_global-2)
      Elseif (id == 4) Then
        F(:,:,1) = buffer_ce(:,:)       ! U_global(:,:,nzg_global-2)
      End If
    End If
    If ( is_last ) Then
      If (id == 3) Then
        F(:,:,nz   ) = buffer_wi(:,:)   ! W_global(:,:,2) 
      Elseif (id == 1) Then
        F(:,:,nzg-1) = buffer_ui(:,:,2) ! U_global(:,:,2)
        F(:,:,nzg  ) = buffer_ui(:,:,3) ! U_global(:,:,3)
      Elseif (id == 2) Then
        F(:,:,nzg-1) = buffer_vi(:,:,2) ! V_global(:,:,2)
        F(:,:,nzg  ) = buffer_vi(:,:,3) ! V_global(:,:,3)
      Elseif (id == 4) Then
        F(:,:,nzg-1) = buffer_ci(:,:,2) ! C_global(:,:,2)
        F(:,:,nzg  ) = buffer_ci(:,:,3) ! C_global(:,:,3)
      End If    
    End If   
    
  End Subroutine apply_periodic_bc_z

  ! Periodicity in y, no MPI needed (y is never domain-decomposed, every rank
  ! already owns the full y-extent); id: 1=y-faces (V), 2=y-centres (U,W); F in/out
  Subroutine apply_periodic_bc_y(F,id)

    Real   (Int64), Intent(InOut) :: F(:,:,:)
    Integer(Int32), Intent(In)    :: id

    !$acc kernels present(F)
    If ( id==1 ) Then
       ! F defined at y faces (V): face 1 and face ny are the same physical
       ! periodic location (y=0 == y=Ly); each takes the opposite side's
       ! nearest interior face
       F(:, 1,:) = F(:,ny-1,:)
       F(:,ny,:) = F(:,   2,:)
    Else
       ! F defined at y centres (U,W): 3-point wraparound
       F(:,    1,:) = F(:,nyg-2,:)
       F(:,nyg-1,:) = F(:,    2,:)
       F(:,nyg  ,:) = F(:,    3,:)
    End If
    !$acc end kernels

  End Subroutine apply_periodic_bc_y

  ! Dirichlet boundary condition in y (no MPI): F; id: 1=y-faces, 2=y-centres; F in/out
  Subroutine apply_Dirichlet_bc_y(F,id)

    Real   (Int64), Intent(InOut) :: F(:,:,:)
    Integer(Int32), Intent(In)    :: id

    If ( id==1 ) Then
       ! F defined at y faces
       F(:, 1,:) = 0d0
       F(:,ny,:) = 0d0
    Else
       ! F defined at y centers
       F(:,  1,:) = -F(:,    2,:)
       F(:,nyg,:) = -F(:,nyg-1,:)
    End If

  End Subroutine apply_Dirichlet_bc_y

  ! Neumann boundary condition in y (no MPI): F; id: 1=y-faces, 2=y-centres; F in/out
  Subroutine apply_Neumann_bc_y(F,id)

    Real   (Int64), Intent(InOut) :: F(:,:,:)
    Integer(Int32), Intent(In)    :: id

    If ( id==1 ) Then
       ! F defined at y faces (first order)
       F(:, 1,:) = F(:,   2,:) 
       F(:,ny,:) = F(:,ny-1,:) 
    Else
       ! F defined at y centers (second order)
       F(:,  1,:) = F(:,    2,:) 
       F(:,nyg,:) = F(:,nyg-1,:) 
    End If

  End Subroutine apply_Neumann_bc_y

  ! Non-zero Neumann boundary condition in y (no MPI): F, alpha (slip length); id: 1=y-faces, 2=y-centres; F in/out
  Subroutine apply_nonzero_Neumann_bc_y(F,alpha,id)

    Real   (Int64), Intent(InOut) :: F(:,:,:)
    Real   (Int64), Intent(In)    :: alpha(:,:,:)
    Integer(Int32), Intent(In)    :: id

    If ( id==1 ) Then
       ! F defined at y faces (first order)
       F(:, 1,:) = F(:,   2,:) - (y (2) - y(   1))*alpha(:,1,:)
       F(:,ny,:) = F(:,ny-1,:) - (y(ny) - y(ny-1))*alpha(:,2,:)
    Else
       ! F defined at y centers (second order)
       F(:,  1,:) = F(:,    2,:) - (yg  (2) - yg(    1))*alpha(:,1,:)
       F(:,nyg,:) = F(:,nyg-1,:) - (yg(nyg) - yg(nyg-1))*alpha(:,2,:)
    End If

  End Subroutine apply_nonzero_Neumann_bc_y

  ! Robin (slip-length) BC in y; formulas
  Subroutine apply_Robin_bc_y(F,alpha,id)

    Real   (Int64), Intent(InOut) :: F(:,:,:)
    Real   (Int64), Intent(In)    :: alpha(:,:,:)
    Integer(Int32), Intent(In)    :: id

    !$acc kernels present(F,alpha,y,yg)
    If ( id==1 ) Then
       ! F defined at y faces (this is first order, could move to second in future)
       F(:, 1,:) = alpha(:,1,:)*F(:,   2,:) / (y( 2)-y(   1)) / ( alpha(:,1,:)/(y( 2)-y(   1)) + 1d0 )
       ! Version for same alpha at each wall:
       F(:,ny,:) = alpha(:,2,:)*F(:,ny-1,:) / (y(ny)-y(ny-1)) / ( alpha(:,2,:)/(y(ny)-y(ny-1)) + 1d0 )
    Else
       ! F defined at y centers (this is second order)
       F(:,  1,:) = ( 2d0*alpha(:,1,:)/(yg(  2) - yg(    1)) - 1d0 )*F(:,    2,:) / ( 2d0*alpha(:,1,:)/(yg(  2) - yg(    1)) + 1d0 )
       ! Version for same alpha at each wall:
       F(:,nyg,:) = ( 2d0*alpha(:,2,:)/(yg(nyg) - yg(nyg-1)) - 1d0 )*F(:,nyg-1,:) / ( 2d0*alpha(:,2,:)/(yg(nyg) - yg(nyg-1)) + 1d0 )
    End If
    !$acc end kernels
    
  End Subroutine apply_Robin_bc_y

  !          Update ghost interior planes (z-direction, same row/adjacent column)
  !          Both directions (+z,-z) are exchanged concurrently via non-blocking MPI instead of two sequential blocking round trips
  Subroutine update_ghost_interior_planes(F,id)

    Real   (Int64), Intent(InOut) :: F(:,:,:)
    Integer(Int32), Intent(In)    :: id

    Integer(Int32) :: up, down
    Integer(Int32) :: reqs(4)

    Call z_halo_neighbors(up, down)

    If (id == 1) Then
      ! update U
      buffer_us(:,:,1) = F(:,:,nzg-1) ! send buffer, towards +z
      buffer_us(:,:,2) = F(:,:,2)     ! send buffer, towards -z
      Call Mpi_irecv(buffer_ur(:,:,1), nx*nyg, Mpi_real8, down, 0, MPI_COMM_WORLD, reqs(1), ierr)
      Call Mpi_irecv(buffer_ur(:,:,2), nx*nyg, Mpi_real8, up,   0, MPI_COMM_WORLD, reqs(2), ierr)
      Call Mpi_isend(buffer_us(:,:,1), nx*nyg, Mpi_real8, up,   0, MPI_COMM_WORLD, reqs(3), ierr)
      Call Mpi_isend(buffer_us(:,:,2), nx*nyg, Mpi_real8, down, 0, MPI_COMM_WORLD, reqs(4), ierr)
      Call Mpi_waitall(4, reqs, MPI_STATUSES_IGNORE, ierr)
      If ( down /= MPI_PROC_NULL ) F(:,:,1)   = buffer_ur(:,:,1) ! received from -z neighbour
      If ( up   /= MPI_PROC_NULL ) F(:,:,nzg) = buffer_ur(:,:,2) ! received from +z neighbour

    Elseif (id == 2) Then
      ! update V
      buffer_vs(:,:,1) = F(:,:,nzg-1) ! send buffer, towards +z
      buffer_vs(:,:,2) = F(:,:,2)     ! send buffer, towards -z
      Call Mpi_irecv(buffer_vr(:,:,1), nxg*ny, Mpi_real8, down, 0, MPI_COMM_WORLD, reqs(1), ierr)
      Call Mpi_irecv(buffer_vr(:,:,2), nxg*ny, Mpi_real8, up,   0, MPI_COMM_WORLD, reqs(2), ierr)
      Call Mpi_isend(buffer_vs(:,:,1), nxg*ny, Mpi_real8, up,   0, MPI_COMM_WORLD, reqs(3), ierr)
      Call Mpi_isend(buffer_vs(:,:,2), nxg*ny, Mpi_real8, down, 0, MPI_COMM_WORLD, reqs(4), ierr)
      Call Mpi_waitall(4, reqs, MPI_STATUSES_IGNORE, ierr)
      If ( down /= MPI_PROC_NULL ) F(:,:,1)   = buffer_vr(:,:,1) ! received from -z neighbour
      If ( up   /= MPI_PROC_NULL ) F(:,:,nzg) = buffer_vr(:,:,2) ! received from +z neighbour

    Elseif (id == 3) Then
      ! update W
      buffer_ws(:,:,1) = F(:,:,nz-1) ! send buffer, towards +z
      buffer_ws(:,:,2) = F(:,:,2)    ! send buffer, towards -z
      Call Mpi_irecv(buffer_wr(:,:,1), nxg*nyg, Mpi_real8, down, 0, MPI_COMM_WORLD, reqs(1), ierr)
      Call Mpi_irecv(buffer_wr(:,:,2), nxg*nyg, Mpi_real8, up,   0, MPI_COMM_WORLD, reqs(2), ierr)
      Call Mpi_isend(buffer_ws(:,:,1), nxg*nyg, Mpi_real8, up,   0, MPI_COMM_WORLD, reqs(3), ierr)
      Call Mpi_isend(buffer_ws(:,:,2), nxg*nyg, Mpi_real8, down, 0, MPI_COMM_WORLD, reqs(4), ierr)
      Call Mpi_waitall(4, reqs, MPI_STATUSES_IGNORE, ierr)
      If ( down /= MPI_PROC_NULL ) F(:,:,1)  = buffer_wr(:,:,1) ! received from -z neighbour
      If ( up   /= MPI_PROC_NULL ) F(:,:,nz) = buffer_wr(:,:,2) ! received from +z neighbour

    Elseif (id == 4) Then
      ! update term (scalar)
      buffer_ws(:,:,1) = F(:,:,nzg-1) ! send buffer, towards +z
      buffer_ws(:,:,2) = F(:,:,2)     ! send buffer, towards -z
      Call Mpi_irecv(buffer_wr(:,:,1), nxg*nyg, Mpi_real8, down, 0, MPI_COMM_WORLD, reqs(1), ierr)
      Call Mpi_irecv(buffer_wr(:,:,2), nxg*nyg, Mpi_real8, up,   0, MPI_COMM_WORLD, reqs(2), ierr)
      Call Mpi_isend(buffer_ws(:,:,1), nxg*nyg, Mpi_real8, up,   0, MPI_COMM_WORLD, reqs(3), ierr)
      Call Mpi_isend(buffer_ws(:,:,2), nxg*nyg, Mpi_real8, down, 0, MPI_COMM_WORLD, reqs(4), ierr)
      Call Mpi_waitall(4, reqs, MPI_STATUSES_IGNORE, ierr)
      If ( down /= MPI_PROC_NULL ) F(:,:,1)   = buffer_wr(:,:,1) ! received from -z neighbour
      If ( up   /= MPI_PROC_NULL ) F(:,:,nzg) = buffer_wr(:,:,2) ! received from +z neighbour
    End if

  End Subroutine update_ghost_interior_planes

  !          Update ghost interior planes (x-direction, same column/adjacent row)
  ! Both directions (+x,-x) are exchanged concurrently via non-blocking MPI instead of two sequential blocking round trips
  Subroutine update_ghost_interior_planes_x(F,id)

    Real   (Int64), Intent(InOut) :: F(:,:,:)
    Integer(Int32), Intent(In)    :: id

    Integer(Int32) :: up, down
    Integer(Int32) :: n2, n3
    Integer(Int32) :: reqs(4)
    Real(Int64), Allocatable :: bs1(:,:), bs2(:,:), br1(:,:), br2(:,:)

    Call x_halo_neighbors(up, down)

    n2 = Size(F,2)
    n3 = Size(F,3)
    Allocate ( bs1(n2,n3), bs2(n2,n3), br1(n2,n3), br2(n2,n3) )

    If ( id == 1 ) Then
       ! F defined at x faces (U): interior faces are 2..nx-1
       bs1 = F(nx-1,:,:) ! send buffer, towards +x
       bs2 = F(2,:,:)    ! send buffer, towards -x
    Else
       ! F defined at x centres (V,W,scalars): interior centres are 2..nxg-1
       bs1 = F(nxg-1,:,:) ! send buffer, towards +x
       bs2 = F(2,:,:)     ! send buffer, towards -x
    End If

    Call Mpi_irecv(br1, n2*n3, Mpi_real8, down, 0, MPI_COMM_WORLD, reqs(1), ierr)
    Call Mpi_irecv(br2, n2*n3, Mpi_real8, up,   0, MPI_COMM_WORLD, reqs(2), ierr)
    Call Mpi_isend(bs1, n2*n3, Mpi_real8, up,   0, MPI_COMM_WORLD, reqs(3), ierr)
    Call Mpi_isend(bs2, n2*n3, Mpi_real8, down, 0, MPI_COMM_WORLD, reqs(4), ierr)
    Call Mpi_waitall(4, reqs, MPI_STATUSES_IGNORE, ierr)

    If ( id == 1 ) Then
       If ( down /= MPI_PROC_NULL ) F(1,:,:)  = br1 ! received from -x neighbour
       If ( up   /= MPI_PROC_NULL ) F(nx,:,:) = br2 ! received from +x neighbour
    Else
       If ( down /= MPI_PROC_NULL ) F(1,:,:)   = br1 ! received from -x neighbour
       If ( up   /= MPI_PROC_NULL ) F(nxg,:,:) = br2 ! received from +x neighbour
    End If

    Deallocate(bs1, bs2, br1, br2)

  End Subroutine update_ghost_interior_planes_x

End Module boundary_conditions
