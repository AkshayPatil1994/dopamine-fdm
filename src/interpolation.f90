! Module for interpolation
Module interpolation

  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use mpi,             Only : ierr, MPI_COMM_WORLD
  Use global,          Only : weight_y_0, weight_y_1, weight_z_0, weight_z_1

  Implicit None

Contains

  !> Linear interpolation of u in x (uniform mesh); di: 1=faces to center, 2=centers to faces
  Pure Subroutine interpolate_x(u,ui,di)

    Real    (Int64), Intent(In)  :: u(:,:,:)
    Real    (Int64), Intent(Out) :: ui(:,:,:)
    Integer (Int32), Intent(In)  :: di

    Integer (Int32) :: n1

    ! Uniform x-grid: simple 0.5 average; n1 computed on host since Size() on an assumed-shape dummy inside
    ! an ACC region crashed under nvfortran
    n1 = Size(u,1)
    !$acc kernels present(u,ui)
    ui(1:n1-1, :, :) = 0.5d0*( u(1:n1-1,:,:) + u(2:n1,:,:) )
    !$acc end kernels

  End Subroutine interpolate_x

  !> Linear interpolation of u in y; di: 1=faces to center (normal average), 2=centers to faces (weighted
  !  average)
  Subroutine interpolate_y(u,ui,di)

    Real    (Int64), Intent(In)  :: u(:,:,:)
    Real    (Int64), Intent(Out) :: ui(:,:,:)
    Integer (Int32), Intent(In)  :: di

    Integer (Int32) :: n2, n3, i2, i3

    n2 = Size(u, 2)
    n3 = Size(u, 3)

    If ( di==1 ) Then
       ! faces to centers (normal average)
       !$acc kernels present(u,ui)
       ui(:, 1:n2-1, :) = 0.5d0*( u(:,1:n2-1,:) + u(:,2:n2,:) )
       !$acc end kernels
    Elseif ( di==2 ) Then
       ! centers to faces (weighted average)
       !$acc kernels present(u,ui,weight_y_0,weight_y_1)
       Do i3 = 1, n3
          Do i2 = 1, n2-1
             ui(:, i2, i3) = weight_y_0(i2)*u(:,i2,i3) + weight_y_1(i2)*u(:,i2+1,i3)
          End Do
       End Do
       !$acc end kernels
    Else
       Write(*,'(A)') 'Error: invalid interpolation'
       Call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    End If
   
  End Subroutine interpolate_y

  !> Linear interpolation of u in z; di: 1=faces to center (centres are face midpoints, plain average),
  !  2=centers to faces (weighted, z may be stretched)
  Pure Subroutine interpolate_z(u,ui,di)

    Real    (Int64), Intent(In)  :: u(:,:,:)
    Real    (Int64), Intent(Out) :: ui(:,:,:)
    Integer (Int32), Intent(In)  :: di

    Integer (Int32) :: n3, i3

    n3 = Size(u,3)
    If ( di==2 ) Then
       !$acc kernels present(u,ui,weight_z_0,weight_z_1)
       Do i3 = 1, n3-1
          ui(:, :, i3) = weight_z_0(i3)*u(:,:,i3) + weight_z_1(i3)*u(:,:,i3+1)
       End Do
       !$acc end kernels
    Else
       ! n3 computed on host, see interpolate_x
       !$acc kernels present(u,ui)
       ui(:, :, 1:n3-1) = 0.5d0*( u(:,:,1:n3-1) + u(:,:,2:n3) )
       !$acc end kernels
    End If

  End Subroutine interpolate_z

End Module interpolation
