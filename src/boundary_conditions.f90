!   Module for boundary conditions
Module boundary_conditions

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use mpi
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

    ! interior region
    Call update_ghost_interior_planes(U,1)
    Call update_ghost_interior_planes(V,2)
    Call update_ghost_interior_planes(W,3)

    ! x direction: periodic, or Dirichlet-SEM inflow / convective outflow
    If ( x_bc_type == 0 ) Then
       Call apply_periodic_bc_x(U,1)
       Call apply_periodic_bc_x(V,2)
       Call apply_periodic_bc_x(W,2)
    Else
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

    ! U boundary condition at the wall
    Call apply_Robin_bc_y(U,alpha_x,2)

    ! V boundary condition at the wall
    Call apply_Robin_bc_y(V,alpha_y,1)

    ! W boundary condition at the wall
    Call apply_Robin_bc_y(W,alpha_z,2)

    ! compute boundary conditions for pseudo-pressure
  End Subroutine apply_boundary_conditions
  
  ! Periodicity in x (no MPI): F is the array; id: 1=x-faces, 2=x-centres; F is in/out
  Subroutine apply_periodic_bc_x(F,id)

    Real   (Int64), Intent(InOut) :: F(:,:,:)
    Integer(Int32), Intent(In)    :: id

    !$acc kernels present(F)
    If ( id==1 ) Then
       ! F defined at x faces
       F( 1,:,:) = F(nx-1,:,:)
       F(nx,:,:) = F(   2,:,:)
    Else
       ! F defined at x centers
       F(    1,:,:) = F(nxg-2,:,:)
       F(nxg-1,:,:) = F(    2,:,:) ! see note*
       F(nxg  ,:,:) = F(    3,:,:)
    End If
    !$acc end kernels

    ! Note: handles non-periodic ICs; redundant after the first time step.
  End Subroutine apply_periodic_bc_x

  ! Dirichlet inflow (x_bc_type==1, no MPI): F face/ghost = mean profile + SEM fluctuation (sem.f90); comp: 1=U,2=V,3=W; F in/out
  Subroutine apply_inflow_bc_x(F,comp)

    Real   (Int64), Intent(InOut) :: F(:,:,:)
    Integer(Int32), Intent(In)    :: comp

    Integer(Int32) :: j, k, n2, n3
    Real   (Int64) :: up, vp, wp, target_val

    n2 = Size(F,2)
    n3 = Size(F,3)

    ! F's y/z coords depend on staggered-grid placement: U=(yc,zc), V=(yf,zc), W=(yc,zf); every (j,k) is independent -- sem_fluctuation/mean_profile_U are !$acc routine seq
    !$acc parallel loop collapse(2) present(F,y,yg,z,zg) private(up,vp,wp,target_val)
    Do k = 1, n3
       Do j = 1, n2
          If ( comp == 1 ) Then
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

  ! Convective outflow (x_bc_type==1, no MPI): dF/dt + Uc*dF/dx = 0, explicit upwind from interior; pairs with homogeneous Dirichlet pressure BC in the DCT-IV solve; F (U,V, or W) in/out
  Subroutine apply_outflow_bc_x(F,Uc)

    Real(Int64), Intent(InOut) :: F(:,:,:)
    Real(Int64), Intent(In)    :: Uc
    Integer(Int32) :: nlast
    Real   (Int64) :: courant

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
    Integer(Int64) :: n(3)

    ! nprocs==1: myid self-pairs, so wrap around via direct local copy instead of Mpi_sendrecv (as in apply_periodic_bc_x)
    If ( nprocs == 1 ) Then
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
    If ( myid==0 ) Then
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

    If ( myid==nprocs-1 ) Then
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
    If ( myid==0 ) Then
      If (id == 3) Then
        ! Send/receive W
        Call Mpi_sendrecv(buffer_wi, nxg*nyg, Mpi_real8, nprocs-1, 3,  &
        buffer_we, nxg*nyg, Mpi_real8, nprocs-1, 3, MPI_COMM_WORLD,    &
        istat, ierr)
      Elseif (id == 1) Then
        ! Send/receive U 
        Call Mpi_sendrecv(buffer_ui, nx*nyg*2, Mpi_real8, nprocs-1, 1, &
        buffer_ue, nx*nyg, Mpi_real8, nprocs-1, 1, MPI_COMM_WORLD,     &
        istat, ierr)
      Elseif (id == 2) Then
        ! Send/receive V 
        Call Mpi_sendrecv(buffer_vi, nxg*ny*2, Mpi_real8, nprocs-1, 2, &
        buffer_ve, nxg*ny, Mpi_real8, nprocs-1, 2, MPI_COMM_WORLD,     &
        istat, ierr)
      Elseif (id == 4) Then
        ! Send/receive C 
        Call Mpi_sendrecv(buffer_ci, nxg*nyg*2, Mpi_real8, nprocs-1, 4, &
        buffer_ce, nxg*nyg, Mpi_real8, nprocs-1, 4, MPI_COMM_WORLD,     &
        istat, ierr)
      End If
    End If

    If ( myid==nprocs-1 ) Then
      If (id == 3) Then
        ! Send/receive W
        Call Mpi_sendrecv(buffer_we, nxg*nyg, Mpi_real8, 0, 3, &
        buffer_wi, nxg*nyg, Mpi_real8, 0, 3, MPI_COMM_WORLD,   &
        istat, ierr)
      Elseif (id == 1) Then
        ! Send/receive U 
        Call Mpi_sendrecv(buffer_ue, nx*nyg, Mpi_real8, 0, 1,  &
        buffer_ui, nx*nyg*2, Mpi_real8, 0, 1, MPI_COMM_WORLD,  &
        istat, ierr)
      Elseif (id == 2) Then
        ! Send/receive V 
        Call Mpi_sendrecv(buffer_ve, nxg*ny, Mpi_real8, 0, 2,  &
        buffer_vi, nxg*ny*2, Mpi_real8, 0, 2, MPI_COMM_WORLD,  &
        istat, ierr)
      Elseif (id == 4) Then
        ! Send/receive C 
        Call Mpi_sendrecv(buffer_ce, nxg*nyg, Mpi_real8, 0, 4,  &
        buffer_ci, nxg*nyg*2, Mpi_real8, 0, 4, MPI_COMM_WORLD,  &
        istat, ierr)
      End If
    End If

    ! apply conditions
    If ( myid==0 ) Then
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
    If ( myid==nprocs-1 ) Then
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

  !          Update ghost interior planes
  Subroutine update_ghost_interior_planes(F,id)

    Real   (Int64), Intent(InOut) :: F(:,:,:)    
    Integer(Int32), Intent(In)    :: id

    Integer(Int32) :: sendto, recvfrom
    Integer(Int32) :: tagto,  tagfrom
    
    If (id == 1) Then
      ! update U
      ! send to top processor, receive from bottom one
      sendto   = myid + 1
      tagto    = myid + 1
      recvfrom = myid - 1
      tagfrom  = myid 
      If ( myid==0 ) Then 
         recvfrom = MPI_PROC_NULL
         tagfrom  = MPI_ANY_TAG
      End If
      If ( myid==nprocs-1 ) Then
         sendto = MPI_PROC_NULL
         tagto  = 0
      End If
      buffer_us = F(:,:,nzg-1) ! send buffer
      Call Mpi_sendrecv(buffer_us, nx*nyg, Mpi_real8, sendto, tagto,        &
           buffer_ur, nx*nyg, Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
           istat, ierr)   
      If ( myid/=0 ) F(:,:,1) = buffer_ur ! received buffer
      
      ! send to bottom processor, receive from top one
      sendto   = myid - 1
      tagto    = myid - 1
      recvfrom = myid + 1
      tagfrom  = myid 
      If ( myid==0 ) Then
         sendto = MPI_PROC_NULL
         tagto  = 0
      End If
      If ( myid==nprocs-1 ) Then
         recvfrom = MPI_PROC_NULL
         tagfrom  = MPI_ANY_TAG
      End If
      buffer_us = F(:,:,2)  ! send buffer
      Call Mpi_sendrecv(buffer_us, nx*nyg, Mpi_real8, sendto, tagto,        &
           buffer_ur, nx*nyg, Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
           istat, ierr)   
      If ( myid/=nprocs-1 ) F(:,:,nzg) = buffer_ur ! received buffer

    Elseif (id == 2) Then
      ! update V
      ! send to top processor, receive from bottom one
      sendto   = myid + 1
      tagto    = myid + 1
      recvfrom = myid - 1
      tagfrom  = myid 
      If ( myid==0 ) Then
         recvfrom = MPI_PROC_NULL
         tagfrom  = MPI_ANY_TAG
      End If
      If ( myid==nprocs-1 ) Then
         sendto = MPI_PROC_NULL
         tagto  = 0
      End If
      buffer_vs = F(:,:,nzg-1) ! send buffer
      Call Mpi_sendrecv(buffer_vs, nxg*ny, Mpi_real8, sendto, tagto,        &
           buffer_vr, nxg*ny, Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
           istat, ierr)   
      If ( myid/=0 ) F(:,:,1) = buffer_vr ! received buffer
      
      ! send to bottom processor, receive from top one
      sendto   = myid - 1
      tagto    = myid - 1
      recvfrom = myid + 1
      tagfrom  = myid 
      If ( myid==0 ) Then
         sendto = MPI_PROC_NULL
         tagto  = 0
      End If
      If ( myid==nprocs-1 ) Then
         recvfrom = MPI_PROC_NULL
         tagfrom  = MPI_ANY_TAG
      End If
      buffer_vs = F(:,:,2)  ! send buffer
      Call Mpi_sendrecv(buffer_vs, nxg*ny, Mpi_real8, sendto, tagto,        &
           buffer_vr, nxg*ny, Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
           istat, ierr)   
      If ( myid/=nprocs-1 ) F(:,:,nzg) = buffer_vr ! received buffer
      
    Elseif (id == 3) Then
      ! update W
      ! send to top processor, receive from bottom one
      sendto   = myid + 1
      tagto    = myid + 1
      recvfrom = myid - 1
      tagfrom  = myid 
      If ( myid==0 ) Then
         recvfrom = MPI_PROC_NULL
         tagfrom  = MPI_ANY_TAG
      End If
      If ( myid==nprocs-1 ) Then
         sendto = MPI_PROC_NULL
         tagto  = 0
      End If
      buffer_ws = F(:,:,nz-1)    ! send buffer
      Call Mpi_sendrecv(buffer_ws, nxg*nyg, Mpi_real8, sendto, tagto,        &
           buffer_wr, nxg*nyg, Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
           istat, ierr)   
      If ( myid/=0 ) F(:,:,1) = buffer_wr ! received buffer
      
      ! send to bottom processor, receive from top one
      sendto   = myid - 1
      tagto    = myid - 1
      recvfrom = myid + 1
      tagfrom  = myid 
      If ( myid==0 ) Then
         sendto = MPI_PROC_NULL
         tagto  = 0
      End If
      If ( myid==nprocs-1 ) Then
         recvfrom = MPI_PROC_NULL
         tagfrom  = MPI_ANY_TAG
      End If
      buffer_ws = F(:,:,2)  ! send buffer
      Call Mpi_sendrecv(buffer_ws, nxg*nyg, Mpi_real8, sendto, tagto,        &
           buffer_wr, nxg*nyg, Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
           istat, ierr)   
      If ( myid/=nprocs-1 ) F(:,:,nz) = buffer_wr ! received buffer     

    Elseif (id == 4) Then
      ! update term
      ! send to top processor, receive from bottom one
      sendto   = myid + 1
      tagto    = myid + 1
      recvfrom = myid - 1
      tagfrom  = myid 
      If ( myid==0 ) Then
         recvfrom = MPI_PROC_NULL
         tagfrom  = MPI_ANY_TAG
      End If
      If ( myid==nprocs-1 ) Then
         sendto = MPI_PROC_NULL
         tagto  = 0
      End If
      buffer_ws = F(:,:,nzg-1)    ! send buffer
      Call Mpi_sendrecv(buffer_ws, nxg*nyg, Mpi_real8, sendto, tagto,        &
           buffer_wr, nxg*nyg, Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
           istat, ierr)   
      If ( myid/=0 ) F(:,:,1) = buffer_wr ! received buffer
      
      ! send to bottom processor, receive from top one
      sendto   = myid - 1
      tagto    = myid - 1
      recvfrom = myid + 1
      tagfrom  = myid 
      If ( myid==0 ) Then
         sendto = MPI_PROC_NULL
         tagto  = 0
      End If
      If ( myid==nprocs-1 ) Then
         recvfrom = MPI_PROC_NULL
         tagfrom  = MPI_ANY_TAG
      End If
      buffer_ws = F(:,:,2)  ! send buffer
      Call Mpi_sendrecv(buffer_ws, nxg*nyg, Mpi_real8, sendto, tagto,        &
           buffer_wr, nxg*nyg, Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
           istat, ierr)   
      If ( myid/=nprocs-1 ) F(:,:,nzg) = buffer_wr ! received buffer     
    End if  
    
  End Subroutine update_ghost_interior_planes

End Module boundary_conditions
