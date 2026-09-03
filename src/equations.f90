! Module to compute right-hand side of Navier-Stokes eq.
Module equations

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global,          Only : x, xm, xg, y, ym, yg, z, zm, zg, term_1, &
                              term_2, term, nx, nxg, ny, nyg, nz, nzg, &
                              nu, dPdx, dPdz, yg_m, nu_t, in1, in2,    &
                              weight_y_0, weight_y_1, dx, dz,         &
                              boussinesq_flag, beta_T, grav, T_ref, Tscal, &
                              advection_scheme, uav_active
  Use interpolation
  Use uav_actuator, Only : apply_uav_forcing
  
  ! prevent implicit typing
  Implicit None
  
Contains

  ! Compute RHS for du/dt: du/dt = -du^2/dx - duv/dy - duw/dz + div(nu grad(u)) + 2*Omega_z*v
  Subroutine compute_rhs_u(U_,V_,W_,rhs_u)

    Real(Int64), Dimension(nx,nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg,ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg,nz), Intent(In) :: W_
    Real(Int64), Dimension(2:nx-1,2:nyg-1,2:nzg-1), Intent(Out) :: rhs_u

    ! local variables
    Integer(Int32) :: i, j, k
    Real   (Int64) :: dy_3
    Real   (Int64) :: nu_x1, nu_x2, nu_y1, nu_y2, nu_z1, nu_z2
    ! Uniform x/z grid: all face/centre spacings equal dx and dz.
    Real   (Int64) :: inv_dx, inv_dz, inv_dx2
    Real   (Int64) :: inv_dy_j   ! 1/(y(j)-y(j-1))
    Real   (Int64) :: inv_dy3    ! 1/dy_3
    Real   (Int64) :: w_adv, w_div   ! convective blend weights, see advection_scheme

    inv_dx  = 1d0 / dx
    inv_dz  = 1d0 / dz
    inv_dx2 = inv_dx * inv_dx   ! used for second derivative in x: 1/dx^2

    ! advective-form weight: 0.5 gives skew-symmetric (default), 0 gives pure divergence-form central
    w_adv = 0.5d0
    If ( advection_scheme == 1 ) w_adv = 0d0
    w_div = 1d0 - w_adv

    ! compute convective terms

    ! 1)---------compute -du^2/dx--------------

    ! interpolate u in x (faces to centers)
    Call interpolate_x(U_,term_1,in1)

    ! -du^2/dx: w_div*(divergence form) + w_adv*(advective form u*du/dx), skew-symmetric suppresses aliasing near sharp gradients
    !$acc parallel loop collapse(3) present(term_1,rhs_u,U_)
    Do k=2,nzg-1
       Do j=2,nyg-1
          Do i=2,nx-1
             rhs_u(i,j,k) = -( w_div*( term_1(i,j,k)*term_1(i,j,k) - term_1(i-1,j,k)*term_1(i-1,j,k) ) + &
                                w_adv*U_(i,j,k)*( term_1(i,j,k) - term_1(i-1,j,k) ) )*inv_dx
          End Do
       End Do
    End Do
    !$acc end parallel loop

    ! 2)-----------compute -duv/dy--------------

    ! interpolate u in y (centers to faces)
    Call interpolate_y(U_,term_1(1:nx,:,:),in2)

    ! interpolate v in x (centers to faces)
    Call interpolate_x(V_,term_2(:,1:ny,:),in2)

    ! -duv/dy: w_div*(divergence form) + w_adv*(advective form v*du/dy)
    !$acc parallel loop collapse(3) present(term_1,term_2,rhs_u,y)
    Do k=2,nzg-1
       Do j=2,nyg-1
          Do i=2,nx-1
             rhs_u(i,j,k) = rhs_u(i,j,k) - ( w_div*( term_1(i,j,k)*term_2(i,j,k) - term_1(i,j-1,k)*term_2(i,j-1,k) ) + &
                  w_adv*0.5d0*(term_2(i,j-1,k)+term_2(i,j,k))*( term_1(i,j,k) - term_1(i,j-1,k) ) ) &
                  / ( y(j) - y(j-1) )
          End Do
       End Do
    End Do
    !$acc end parallel loop

    ! 3)--------------compute -duw/dz--------------

    ! interpolate u in z (centers to faces)
    Call interpolate_z(U_,term_1(1:nx,:,:),in2)

    ! interpolate w in x (centers to faces)
    Call interpolate_x(W_,term_2(:,:,1:nz),in2)

    ! -duw/dz: w_div*(divergence form) + w_adv*(advective form w*du/dz)
    !$acc parallel loop collapse(3) present(term_1,term_2,rhs_u)
    Do k=2,nzg-1
       Do j=2,nyg-1
          Do i=2,nx-1
             rhs_u(i,j,k) = rhs_u(i,j,k) - ( w_div*( term_1(i,j,k)*term_2(i,j,k) - term_1(i,j,k-1)*term_2(i,j,k-1) ) + &
                  w_adv*0.5d0*(term_2(i,j,k-1)+term_2(i,j,k))*( term_1(i,j,k) - term_1(i,j,k-1) ) )*inv_dz
          End Do
       End Do
    End Do
    !$acc end parallel loop

    ! compute viscous terms

    ! 4)--------compute d( 2*(nu_t+nu)*S1j )/dxj----

    ! first derivation in y: du/dy + dv/dx goes to yg_m(1:ny)
    !$acc kernels present(term_1,U_,yg)
    Do i = 1, nx
       Do k = 1, nzg
          term_1(i,1:nyg-1,k) = ( U_(i,2:nyg,k) - U_(i,1:nyg-1,k) )/( yg(2:nyg) - yg(1:nyg-1) )
       End Do
    End Do
    !$acc end kernels
    ! interpolate du/dy from yg_m(1:ny) to faces y(1:ny): interpolate_y_2nd disabled, term_1 used directly below (no interpolation)

    !$acc parallel loop collapse(2) present(term_1,rhs_u,nu_t,weight_y_0,weight_y_1,U_,V_,W_,y)
    Do k=2,nzg-1

       Do j=2,nyg-1

          ! second derivation in y (U at faces) — non-uniform y only
          dy_3   = y(j) - y(j-1)
          inv_dy3 = 1d0 / dy_3

          Do i=2,nx-1

             ! total viscosity at x locations
             nu_x1  = nu + nu_t(i  ,j,k)
             nu_x2  = nu + nu_t(i+1,j,k)
              
             ! Wall-face viscosity forced to nu only
             If ( j-1 == 1 ) Then
                nu_y1 = nu
             Else
                nu_y1 = nu + 0.5d0*( weight_y_0(j-1)*nu_t(i,  j-1,k) + weight_y_1(j-1)*nu_t(i,  j  ,k) + &
                                     weight_y_0(j-1)*nu_t(i+1,j-1,k) + weight_y_1(j-1)*nu_t(i+1,j  ,k) )
             End If
             If ( j+1 == nyg ) Then
                nu_y2 = nu
             Else
                nu_y2 = nu + 0.5d0*( weight_y_0(j  )*nu_t(i,  j,  k) + weight_y_1(j  )*nu_t(i,  j+1,k) + &
                                     weight_y_0(j  )*nu_t(i+1,j,  k) + weight_y_1(j  )*nu_t(i+1,j+1,k) )
             End If

             ! total viscosity at z locations
             nu_z1 = nu + 0.25d0*(nu_t(i,j,k)+nu_t(i,j,k-1)+nu_t(i+1,j,k)+nu_t(i+1,j,k-1))
             nu_z2 = nu + 0.25d0*(nu_t(i,j,k)+nu_t(i,j,k+1)+nu_t(i+1,j,k)+nu_t(i+1,j,k+1))

             ! viscous term, fused directly into rhs_u (was written to scratch `term` then accumulated separately)
             rhs_u(i,j,k) = rhs_u(i,j,k) +                                                                  &
                           2d0*inv_dx2*(nu_x2*(U_(i+1,j,k)-U_(i,j,k)) - nu_x1*(U_(i,j,k)-U_(i-1,j,k)) )    + & !d(2(nu+nu_t)*du/dx)/dx
                           inv_dy3*(nu_y2*( term_1(i,j  ,k) + (V_(i+1,j  ,k)-V_(i,j  ,k))*inv_dx )         - &
                                     nu_y1*( term_1(i,j-1,k) + (V_(i+1,j-1,k)-V_(i,j-1,k))*inv_dx ) )       + & !d((nu+nu_t)*(du/dy+dv/dx))/dy
                           inv_dz*(nu_z2*( (U_(i,j,k+1)-U_(i,j,k  ))*inv_dz + (W_(i+1,j,k  )-W_(i,j,k  ))*inv_dx ) - &
                                   nu_z1*( (U_(i,j,k  )-U_(i,j,k-1))*inv_dz + (W_(i+1,j,k-1)-W_(i,j,k-1))*inv_dx ) ) + & !d((nu+nu_t)*(du/dz+dw/dx))/dz
                           dPdx ! constant pressure gradient forcing

          End Do
       End Do
    End Do
    !$acc end parallel loop

  End Subroutine compute_rhs_u

  !                       Compute RHS for dv/dt
  ! dv/dt = -duv/dx - dv^2/dy - dvw/dz + div(nu grad(v)) - 2*Omega_z*u
  Subroutine compute_rhs_v(U_,V_,W_,rhs_v)

    Real(Int64), Dimension(nx,nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg,ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg,nz), Intent(In) :: W_
    Real(Int64), Dimension(2:nxg-1,2:ny-1,2:nzg-1), Intent(Out) :: rhs_v

    ! local variables
    Integer(Int32) :: i, j, k
    Real   (Int64) :: dx_1, dx_2, dx_3, maxerr
    Real   (Int64) :: dy_1, dy_2, dy_3
    Real   (Int64) :: dz_1, dz_2, dz_3
    Real   (Int64) :: nu_x1, nu_x2, nu_y1, nu_y2, nu_z1, nu_z2
    Real   (Int64) :: inv_dx, inv_dz, inv_dx2
    Real   (Int64) :: inv_dyg_j               ! 1/(yg(j+1)-yg(j))
    Real   (Int64) :: inv_dy1, inv_dy2, inv_dy3, two_inv_dy3  ! y-spacing inverses
    Real   (Int64) :: w_adv, w_div   ! convective blend weights, see advection_scheme

    inv_dx  = 1d0 / dx
    inv_dz  = 1d0 / dz
    inv_dx2 = inv_dx * inv_dx

    ! advective-form weight: 0.5 gives skew-symmetric (default), 0 gives pure divergence-form central
    w_adv = 0.5d0
    If ( advection_scheme == 1 ) w_adv = 0d0
    w_div = 1d0 - w_adv

    ! compute convective terms

    ! 1)---------compute -dv^2/dy--------------

    ! interpolate v in y (faces to centers)
    Call interpolate_y(V_,term_1,in1)

    ! -dv^2/dy: w_div*(divergence form) + w_adv*(advective form v*dv/dy)
    !$acc parallel loop collapse(3) present(term_1,rhs_v,yg,V_)
    Do k=2,nzg-1
       Do j=2,ny-1
          Do i=2,nxg-1
             rhs_v(i,j,k) = -( w_div*( term_1(i,j,k)*term_1(i,j,k) - term_1(i,j-1,k)*term_1(i,j-1,k) ) + &
                                w_adv*V_(i,j,k)*( term_1(i,j,k) - term_1(i,j-1,k) ) ) / ( yg(j+1) - yg(j) )
          End Do
       End Do
    End Do
    !$acc end parallel loop

    ! interpolate -dv^2/dy from yg_m(2:ny-1) to y(2:end-1)
    !Call interpolate_y_2nd(yg_m(2:ny-1),term_2(:,2:ny-1,:),y(2:ny-1),term(:,2:ny-1,:))

    ! 2)-----------compute -duv/dx--------------

    ! interpolate u in y (centers to faces)
    Call interpolate_y(U_,term_1(1:nx,:,:),in2)

    ! interpolate v in x (centers to faces)
    Call interpolate_x(V_,term_2(:,1:ny,:),in2)

    ! -duv/dx: w_div*(divergence form) + w_adv*(advective form u*dv/dx)
    !$acc parallel loop collapse(3) present(term_1,term_2,rhs_v)
    Do k=2,nzg-1
       Do j=2,ny-1
          Do i=2,nxg-1
             rhs_v(i,j,k) = rhs_v(i,j,k) - ( w_div*( term_1(i,j,k)*term_2(i,j,k) - term_1(i-1,j,k)*term_2(i-1,j,k) ) + &
                  w_adv*0.5d0*(term_1(i-1,j,k)+term_1(i,j,k))*( term_2(i,j,k) - term_2(i-1,j,k) ) )*inv_dx
          End Do
       End Do
    End Do
    !$acc end parallel loop

    ! 3)--------------compute -dvw/dz--------------

    ! interpolate v in z (centers to faces)
    Call interpolate_z(V_,term_1(:,1:ny,:),in2)

    ! interpolate w in y (centers to faces)
    Call interpolate_y(W_,term_2(:,:,1:nz),in2)

    ! -dvw/dz: w_div*(divergence form) + w_adv*(advective form w*dv/dz)
    !$acc parallel loop collapse(3) present(term_1,term_2,rhs_v)
    Do k=2,nzg-1
       Do j=2,ny-1
          Do i=2,nxg-1
             rhs_v(i,j,k) = rhs_v(i,j,k) - ( w_div*( term_1(i,j,k)*term_2(i,j,k) - term_1(i,j,k-1)*term_2(i,j,k-1) ) + &
                  w_adv*0.5d0*(term_2(i,j,k-1)+term_2(i,j,k))*( term_1(i,j,k) - term_1(i,j,k-1) ) )*inv_dz
          End Do
       End Do
    End Do
    !$acc end parallel loop

    ! compute viscous terms

    ! 4)-----------compute d( 2*(nu+nu_t)*S2j )/dxj--

    ! interpolate eddy viscosity to faces

    ! second order remain, no need to interpolate
    !$acc parallel loop collapse(2) present(rhs_v,V_,U_,W_,nu_t,weight_y_0,weight_y_1,y,yg)
    Do k=2,nzg-1

       Do j=2,ny-1

          ! y spacings at cell-faces (dy_1,dy_2) and cell-centres (dy_3); non-uniform y.
          dy_1 = y(  j) - y(j-1)
          dy_2 = y(j+1) - y(j  )
          dy_3 = yg(j+1) - yg(j)
          inv_dy1     = 1d0 / dy_1
          inv_dy2     = 1d0 / dy_2
          inv_dy3     = 1d0 / dy_3
          two_inv_dy3 = 2d0 * inv_dy3

          Do i=2,nxg-1

             ! eddy viscosity at x locations
             nu_x1 = nu + 0.5d0*( weight_y_0(j)*nu_t(i-1,j,k) + weight_y_1(j)*nu_t(i-1,j+1,k) + & 
                                  weight_y_0(j)*nu_t(i  ,j,k) + weight_y_1(j)*nu_t(i  ,j+1,k) )
             nu_x2 = nu + 0.5d0*( weight_y_0(j)*nu_t(i,  j,k) + weight_y_1(j)*nu_t(i,  j+1,k) + & 
                                  weight_y_0(j)*nu_t(i+1,j,k) + weight_y_1(j)*nu_t(i+1,j+1,k) )

             ! eddy viscosity at the centers
             nu_y1 = nu + nu_t(i,j  ,k) 
             nu_y2 = nu + nu_t(i,j+1,k) 

             ! eddy viscosity at z locations
             nu_z1 = nu + 0.5d0*( weight_y_0(j)*nu_t(i,j,k-1) + weight_y_1(j)*nu_t(i,j+1,k-1) + &
                                  weight_y_0(j)*nu_t(i,j,k  ) + weight_y_1(j)*nu_t(i,j+1,k  ) )
             nu_z2 = nu + 0.5d0*( weight_y_0(j)*nu_t(i,j,  k) + weight_y_1(j)*nu_t(i,j+1,  k) + &
                                  weight_y_0(j)*nu_t(i,j,k+1) + weight_y_1(j)*nu_t(i,j+1,k+1) )

             ! viscous term, fused directly into rhs_v (was written to scratch term_2 then accumulated separately)
             rhs_v(i,j,k) = rhs_v(i,j,k) +                                                                                     &
                             inv_dx*(nu_x2*( (V_(i+1,j,k) - V_(i  ,j,k))*inv_dx + (U_(i  ,j+1,k) - U_(i  ,j,k))*inv_dy3 )   - &
                                     nu_x1*( (V_(i,  j,k) - V_(i-1,j,k))*inv_dx + (U_(i-1,j+1,k) - U_(i-1,j,k))*inv_dy3 ) ) + & !d((nu+nu_t)(dv/dx+du/dy))/dx
                             two_inv_dy3*(nu_y2*(V_(i,j+1,k) - V_(i,j  ,k))*inv_dy2 - &
                             nu_y1*(V_(i,j  ,k) - V_(i,j-1,k))*inv_dy1 ) + & !d(2(nu+nu_t)dv/dy)/dy
                             inv_dz*(nu_z2*( (V_(i,j,k+1) - V_(i,j,k  ))*inv_dz + (W_(i,j+1,k  ) - W_(i,j,k  ))*inv_dy3 )     - &
                                     nu_z1*( (V_(i,j,k  ) - V_(i,j,k-1))*inv_dz + (W_(i,j+1,k-1) - W_(i,j,k-1))*inv_dy3 ) )       !d((nu+nu_t)(dv/dz+dw/dy))/dz

          End Do
       End Do
    End Do
    !$acc end parallel loop

    ! Boussinesq buoyancy: gravity acts along -y; T_face is the stretched-grid-aware cell-centre-to-v-face interpolation of Tscal (same weight_y_0/weight_y_1 used for nu_t elsewhere in this file)
    If ( boussinesq_flag >= 1 ) Then
       !$acc parallel loop collapse(3) present(rhs_v,Tscal,weight_y_0,weight_y_1)
       Do k=2,nzg-1
          Do j=2,ny-1
             Do i=2,nxg-1
                rhs_v(i,j,k) = rhs_v(i,j,k) + beta_T*grav*( weight_y_0(j)*Tscal(i,j,k) + weight_y_1(j)*Tscal(i,j+1,k) - T_ref )
             End Do
          End Do
       End Do
       !$acc end parallel loop
    End If

    ! UAV actuator disk: static-disk vertical reaction force (Phase 1, see docs/UAV_ActuatorDisk_Design.md)
    ! apply_uav_forcing is host-only (no OpenACC); on GPU builds rhs_v (Fv1/Fv2/Fv3) is a
    ! separate device allocation from the terms computed above, so it must be pulled to the
    ! host before the host loop touches it and pushed back before the RK update reads it on device.
    If ( uav_active >= 1 ) Then
       !$acc update host(rhs_v)
       Call apply_uav_forcing(rhs_v)
       !$acc update device(rhs_v)
    End If

  End Subroutine compute_rhs_v

  !                Compute RHS for dw/dt
  !  dw/dt = -duw/dx - dvw/dy - dw^2/dz + div(nu grad(w))
  Subroutine compute_rhs_w(U_,V_,W_,rhs_w)

    Real(Int64), Dimension(nx,nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg,ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg,nz), Intent(In) :: W_
    Real(Int64), Dimension(2:nxg-1,2:nyg-1,2:nz-1), Intent(Out) :: rhs_w

    ! local variables
    Integer(Int32) :: i, j, k
    Real   (Int64) :: dy_3
    Real   (Int64) :: nu_x1, nu_x2, nu_y1, nu_y2, nu_z1, nu_z2
    Real   (Int64) :: inv_dx, inv_dz, inv_dz2
    Real   (Int64) :: inv_dy_j   ! 1/(y(j)-y(j-1))
    Real   (Int64) :: inv_dy3    ! 1/dy_3
    Real   (Int64) :: w_adv, w_div   ! convective blend weights, see advection_scheme

    inv_dx  = 1d0 / dx
    inv_dz  = 1d0 / dz
    inv_dz2 = inv_dz * inv_dz

    ! advective-form weight: 0.5 gives skew-symmetric (default), 0 gives pure divergence-form central
    w_adv = 0.5d0
    If ( advection_scheme == 1 ) w_adv = 0d0
    w_div = 1d0 - w_adv

    ! compute convective terms

    ! 1)---------compute -dw^2/dz--------------

    ! interpolate w in z (faces to centers)
    Call interpolate_z(W_,term_1,in1)

    ! -dw^2/dz: w_div*(divergence form) + w_adv*(advective form w*dw/dz)
    !$acc parallel loop collapse(3) present(term_1,rhs_w,W_)
    Do k=2,nz-1
       Do j=2,nyg-1
          Do i=2,nxg-1
             rhs_w(i,j,k) = -( w_div*( term_1(i,j,k)*term_1(i,j,k) - term_1(i,j,k-1)*term_1(i,j,k-1) ) + &
                                w_adv*W_(i,j,k)*( term_1(i,j,k) - term_1(i,j,k-1) ) )*inv_dz
          End Do
       End Do
    End Do
    !$acc end parallel loop

    ! 2)-----------compute -duw/dx--------------

    ! interpolate u in z (centers to faces)
    Call interpolate_z(U_,term_1(1:nx,:,:),in2)

    ! interpolate w in x (centers to faces)
    Call interpolate_x(W_,term_2(:,:,1:nz),in2)

    ! -duw/dx: w_div*(divergence form) + w_adv*(advective form u*dw/dx)
    !$acc parallel loop collapse(3) present(term_1,term_2,rhs_w)
    Do k=2,nz-1
       Do j=2,nyg-1
          Do i=2,nxg-1
             rhs_w(i,j,k) = rhs_w(i,j,k) - ( w_div*( term_1(i,j,k)*term_2(i,j,k) - term_1(i-1,j,k)*term_2(i-1,j,k) ) + &
                  w_adv*0.5d0*(term_1(i-1,j,k)+term_1(i,j,k))*( term_2(i,j,k) - term_2(i-1,j,k) ) )*inv_dx
          End Do
       End Do
    End Do
    !$acc end parallel loop

    ! 3)--------------compute -dvw/dy--------------

    ! interpolate v in z (centers to faces)
    Call interpolate_z(V_,term_1(:,1:ny,:),in2)

    ! interpolate w in y (centers to faces)
    Call interpolate_y(W_,term_2(:,:,1:nz),in2)

    ! -dvw/dy: w_div*(divergence form) + w_adv*(advective form v*dw/dy)
    !$acc parallel loop collapse(3) present(term_1,term_2,rhs_w,y)
    Do k=2,nz-1
       Do j=2,nyg-1
          Do i=2,nxg-1
             rhs_w(i,j,k) = rhs_w(i,j,k) - ( w_div*( term_1(i,j,k)*term_2(i,j,k) - term_1(i,j-1,k)*term_2(i,j-1,k) ) + &
                  w_adv*0.5d0*(term_1(i,j-1,k)+term_1(i,j,k))*( term_2(i,j,k) - term_2(i,j-1,k) ) ) &
                  / ( y(j) - y(j-1) )
          End Do
       End Do
    End Do
    !$acc end parallel loop

    ! compute viscous terms

    ! 4)----------compute d( 2(nu+nu_t)*S3j )/dxj----

    ! first derivation in y: dw/dy goes to yg_m(1:ny)
    !$acc kernels present(term_1,W_,yg)
    Do i = 1, nxg
       Do k = 1, nz
          term_1(i,1:nyg-1,k) = ( W_(i,2:nyg,k) - W_(i,1:nyg-1,k) )/( yg(2:nyg) - yg(1:nyg-1) )
       End Do
    End Do
    !$acc end kernels
    ! interpolate dw/dy from yg_m(1:ny) to faces y(1:ny): interpolate_y_2nd disabled, term_1 used directly below (no interpolation)

    ! nu_t at the required face locations is read directly from the global array
    ! in the loop below; a separate interpolate_y pass is not needed.

    !$acc parallel loop collapse(2) present(term_1,rhs_w,nu_t,weight_y_0,weight_y_1,U_,V_,W_,y)
    Do k=2,nz-1

       Do j=2,nyg-1

          ! second derivation in y (W at faces) — non-uniform y only
          dy_3   = y(j) - y(j-1)
          inv_dy3 = 1d0 / dy_3

          Do i=2,nxg-1

             ! eddy viscosity at x locations
             nu_x1 = nu + 0.25d0*(nu_t(i,j,k)+nu_t(i,j,k+1)+nu_t(i-1,j,k)+nu_t(i-1,j,k+1))
             nu_x2 = nu + 0.25d0*(nu_t(i,j,k)+nu_t(i,j,k+1)+nu_t(i+1,j,k)+nu_t(i+1,j,k+1))

             ! Wall-face viscosity forced to nu only
             If ( j-1 == 1 ) Then
                nu_y1 = nu
             Else
                nu_y1 = nu + 0.5d0*( weight_y_0(j-1)*nu_t(i,j-1,  k) + weight_y_1(j-1)*nu_t(i,j  ,  k) + &
                                     weight_y_0(j-1)*nu_t(i,j-1,k+1) + weight_y_1(j-1)*nu_t(i,j  ,k+1) )
             End If
             If ( j+1 == nyg ) Then
                nu_y2 = nu
             Else
                nu_y2 = nu + 0.5d0*( weight_y_0(j  )*nu_t(i,j  ,  k) + weight_y_1(j  )*nu_t(i,j+1,  k) + &
                                     weight_y_0(j  )*nu_t(i,j  ,k+1) + weight_y_1(j  )*nu_t(i,j+1,k+1) )
             End If

             ! eddy viscosity at z locations
             nu_z1 = nu + nu_t(i,j,k)
             nu_z2 = nu + nu_t(i,j,k+1)

             ! viscous term, fused directly into rhs_w (was written to scratch `term` then accumulated separately)
             rhs_w(i,j,k) = rhs_w(i,j,k) +                                                                                &
                           inv_dx*(nu_x2*( (W_(i+1,j,k) - W_(i  ,j,k))*inv_dx + (U_(i  ,j,k+1) - U_(i  ,j,k))*inv_dz )   - &
                                   nu_x1*( (W_(i  ,j,k) - W_(i-1,j,k))*inv_dx + (U_(i-1,j,k+1) - U_(i-1,j,k))*inv_dz ) ) + & !d((nu+nu_t)*(dw/dx+du/dz))/dx
                           inv_dy3*(nu_y2*( term_1(i,j  ,k) + (V_(i,j  ,k+1) - V_(i,j  ,k))*inv_dz )                   - &
                                     nu_y1*( term_1(i,j-1,k) + (V_(i,j-1,k+1) - V_(i,j-1,k))*inv_dz ) )                 + & !d((nu+nu_t)*(dw/dy+dv/dz))/dy
                           2d0*inv_dz*(nu_z2*(W_(i,j,k+1)-W_(i,j,k))*inv_dz - nu_z1*(W_(i,j,k)-W_(i,j,k-1))*inv_dz )     + & !d(2(nu+nu_t)*dw/dz)/dz
                           dPdz ! constant pressure gradient forcing

          End Do
       End Do
    End Do
    !$acc end parallel loop

  End Subroutine compute_rhs_w

End Module equations
