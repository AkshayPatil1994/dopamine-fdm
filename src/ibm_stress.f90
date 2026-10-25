!     Wall-model stress on the staircase faces of an immersed body
Module ibm_stress

  Use iso_fortran_env, Only : Int32, Int64
  Use global
  Use mpi
  Use ieee_arithmetic, Only : ieee_value, ieee_quiet_nan
  Use wallmodel, Only : solve_u_tau_reichardt, solve_u_tau_rough
  Use ibm, Only : dup_cell_flags

  Implicit None

  ! Kinematic momentum flux through an edge between a fluid and a closed face, replacing the viscous stencil
  ! flux. Last index: U faces (y-lo, y-hi, z-lo, z-hi), V (x-lo, x-hi, z-lo, z-hi), W (x-lo, x-hi, y-lo,
  ! y-hi); ovr_none = no override
  Real(Int64), Parameter :: ovr_none = Huge(1d0)
  Real(Int64), Allocatable, Dimension(:,:,:,:) :: ovr_u, ovr_v, ovr_w
  Real(Int64), Allocatable, Dimension(:,:,:)   :: nu_const

Contains

  !> Allocate the override arrays: full size when the stress is used, one element otherwise (the viscous
  !  kernels always list them); nu_const is the single-phase kinematic viscosity
  Subroutine ibm_stress_init(full)

    Logical, Intent(In) :: full

    If ( Allocated(ovr_u) ) Then
       !$acc exit data delete(ovr_u,ovr_v,ovr_w,nu_const)
       Deallocate( ovr_u, ovr_v, ovr_w, nu_const )
    End If
    If ( full ) Then
       Allocate( ovr_u(nx,nyg,nzg,4), ovr_v(nxg,ny,nzg,4), ovr_w(nxg,nyg,nz,4), nu_const(nxg,nyg,nzg) )
    Else
       Allocate( ovr_u(1,1,1,4), ovr_v(1,1,1,4), ovr_w(1,1,1,4), nu_const(1,1,1) )
    End If
    ovr_u = ovr_none;  ovr_v = ovr_none;  ovr_w = ovr_none;  nu_const = nu
    !$acc enter data copyin(ovr_u,ovr_v,ovr_w,nu_const)

  End Subroutine ibm_stress_init


  !> Log-law wall shear at distance d from the closed face along the tangential velocity (ua along the
  !  component, ub the other); sg = +1 if the closed face is on the low side, -1 on the high side
  Function wall_flux(ua, ub, d, nul, oid, sg) Result(f)
    !$acc routine seq

    Real(Int64),    Intent(In) :: ua, ub, d, nul, sg
    Integer(Int32), Intent(In) :: oid
    Real(Int64) :: f, ut, um

    um = Sqrt(ua*ua + ub*ub)
    f = 0d0
    If ( um == 0d0 ) Return
    If ( ibm_z0(oid) > 0d0 ) Then
       Call solve_u_tau_rough(um, d, ibm_z0(oid), ut)
    Else
       Call solve_u_tau_reichardt(um, d, nul, ut)
    End If
    f = sg*ut*ut*ua/um

  End Function wall_flux


  Pure Function solid_id(ob, i, j, k) Result(oid)
    !$acc routine seq

    Real(Int64),    Intent(In) :: ob(nxg,nyg,nzg)
    Integer(Int32), Intent(In) :: i, j, k
    Integer(Int32) :: oid

    oid = Min(Max(Nint(ob(i,j,k)), 0), max_ibm_objects)

  End Function solid_id


  !> Mean of the open faces among four (mask 1 = open); 0 when all four are closed
  Pure Function open_mean4(f1, f2, f3, f4, m1, m2, m3, m4) Result(f)
    !$acc routine seq

    Real(Int64), Intent(In) :: f1, f2, f3, f4, m1, m2, m3, m4
    Real(Int64) :: f, n

    n = m1 + m2 + m3 + m4
    f = 0d0
    If ( n > 0.5d0 ) f = ( m1*f1 + m2*f2 + m3*f3 + m4*f4 )/n

  End Function open_mean4


  !> Roughness id of the closed face between two cells: that of the solid cell (phi < 0), the larger one when
  !  both are solid
  Pure Function face_id(ph, ob, i1, j1, k1, i2, j2, k2) Result(oid)
    !$acc routine seq

    Real(Int64),    Intent(In) :: ph(nxg,nyg,nzg), ob(nxg,nyg,nzg)
    Integer(Int32), Intent(In) :: i1, j1, k1, i2, j2, k2
    Integer(Int32) :: oid

    oid = 0
    If ( ph(i1,j1,k1) < 0d0 ) oid = solid_id(ob, i1, j1, k1)
    If ( ph(i2,j2,k2) < 0d0 ) oid = Max( oid, solid_id(ob, i2, j2, k2) )

  End Function face_id


  !> Set the flux overrides from the current velocity. mu, mv, mw: 1 (open) / 0 (closed) face masks; nu_cell:
  !  kinematic viscosity at cell centres, <= 0 where the wall model is off. Domain-wall rows stay with the
  !  flat-wall model.
  Subroutine ibm_stress_update(U, V, W, mu, mv, mw, nu_cell)

    Real(Int64), Intent(In) :: U(nx,nyg,nzg), V(nxg,ny,nzg), W(nxg,nyg,nz)
    Real(Int64), Intent(In) :: mu(nx,nyg,nzg), mv(nxg,ny,nzg), mw(nxg,nyg,nz), nu_cell(nxg,nyg,nzg)

    Integer(Int32) :: i, j, k
    Real(Int64) :: nul, ua, ub

    !$acc kernels present(ovr_u,ovr_v,ovr_w)
    ovr_u = ovr_none;  ovr_v = ovr_none;  ovr_w = ovr_none
    !$acc end kernels

    !$acc parallel loop collapse(3) private(nul,ua,ub) present(U,V,W,mu,mv,mw,nu_cell,ovr_u,phi,ibm_obj_id,y,yg,z,zg)
    Do k = 2, nzg-1
       Do j = 2, nyg-1
          Do i = 2, nx-1
             If ( mu(i,j,k) < 0.5d0 ) Cycle
             nul = Min(nu_cell(i,j,k), nu_cell(i+1,j,k))
             If ( nul <= 0d0 ) Cycle
             ua = U(i,j,k)
             ub = open_mean4( W(i,j,k-1), W(i,j,k), W(i+1,j,k-1), W(i+1,j,k), mw(i,j,k-1), mw(i,j,k), mw(i+1,j,k-1), mw(i+1,j,k) )
             If ( j-1 >= 2 ) Then;  If ( mu(i,j-1,k) < 0.5d0 ) &
                ovr_u(i,j,k,1) = wall_flux(ua, ub, yg(j) - y(j-1), nul, face_id(phi, ibm_obj_id, i,j-1,k, i+1,j-1,k), 1d0);  End If
             If ( j+1 <= nyg-1 ) Then;  If ( mu(i,j+1,k) < 0.5d0 ) &
                ovr_u(i,j,k,2) = wall_flux(ua, ub, y(j) - yg(j), nul, face_id(phi, ibm_obj_id, i,j+1,k, i+1,j+1,k), -1d0);  End If
             ub = open_mean4( V(i,j-1,k), V(i,j,k), V(i+1,j-1,k), V(i+1,j,k), mv(i,j-1,k), mv(i,j,k), mv(i+1,j-1,k), mv(i+1,j,k) )
             If ( mu(i,j,k-1) < 0.5d0 ) ovr_u(i,j,k,3) = wall_flux(ua, ub, zg(k) - z(k-1), nul, &
                face_id(phi, ibm_obj_id, i,j,k-1, i+1,j,k-1), 1d0)
             If ( mu(i,j,k+1) < 0.5d0 ) ovr_u(i,j,k,4) = wall_flux(ua, ub, z(k) - zg(k), nul, &
                face_id(phi, ibm_obj_id, i,j,k+1, i+1,j,k+1), -1d0)
          End Do
       End Do
    End Do

    !$acc parallel loop collapse(3) private(nul,ua,ub) present(U,V,W,mu,mv,mw,nu_cell,ovr_v,phi,ibm_obj_id,x,xg,z,zg)
    Do k = 2, nzg-1
       Do j = 2, ny-1
          Do i = 2, nxg-1
             If ( mv(i,j,k) < 0.5d0 ) Cycle
             nul = Min(nu_cell(i,j,k), nu_cell(i,j+1,k))
             If ( nul <= 0d0 ) Cycle
             ua = V(i,j,k)
             ub = open_mean4( W(i,j,k-1), W(i,j,k), W(i,j+1,k-1), W(i,j+1,k), mw(i,j,k-1), mw(i,j,k), mw(i,j+1,k-1), mw(i,j+1,k) )
             If ( mv(i-1,j,k) < 0.5d0 ) ovr_v(i,j,k,1) = wall_flux(ua, ub, xg(i) - x(i-1), nul, &
                face_id(phi, ibm_obj_id, i-1,j,k, i-1,j+1,k), 1d0)
             If ( mv(i+1,j,k) < 0.5d0 ) ovr_v(i,j,k,2) = wall_flux(ua, ub, x(i) - xg(i), nul, &
                face_id(phi, ibm_obj_id, i+1,j,k, i+1,j+1,k), -1d0)
             ub = open_mean4( U(i-1,j,k), U(i,j,k), U(i-1,j+1,k), U(i,j+1,k), mu(i-1,j,k), mu(i,j,k), mu(i-1,j+1,k), mu(i,j+1,k) )
             If ( mv(i,j,k-1) < 0.5d0 ) ovr_v(i,j,k,3) = wall_flux(ua, ub, zg(k) - z(k-1), nul, &
                face_id(phi, ibm_obj_id, i,j,k-1, i,j+1,k-1), 1d0)
             If ( mv(i,j,k+1) < 0.5d0 ) ovr_v(i,j,k,4) = wall_flux(ua, ub, z(k) - zg(k), nul, &
                face_id(phi, ibm_obj_id, i,j,k+1, i,j+1,k+1), -1d0)
          End Do
       End Do
    End Do

    !$acc parallel loop collapse(3) private(nul,ua,ub) present(U,V,W,mu,mv,mw,nu_cell,ovr_w,phi,ibm_obj_id,x,xg,y,yg)
    Do k = 2, nz-1
       Do j = 2, nyg-1
          Do i = 2, nxg-1
             If ( mw(i,j,k) < 0.5d0 ) Cycle
             nul = Min(nu_cell(i,j,k), nu_cell(i,j,k+1))
             If ( nul <= 0d0 ) Cycle
             ua = W(i,j,k)
             ub = open_mean4( V(i,j-1,k), V(i,j,k), V(i,j-1,k+1), V(i,j,k+1), mv(i,j-1,k), mv(i,j,k), mv(i,j-1,k+1), mv(i,j,k+1) )
             If ( mw(i-1,j,k) < 0.5d0 ) ovr_w(i,j,k,1) = wall_flux(ua, ub, xg(i) - x(i-1), nul, &
                face_id(phi, ibm_obj_id, i-1,j,k, i-1,j,k+1), 1d0)
             If ( mw(i+1,j,k) < 0.5d0 ) ovr_w(i,j,k,2) = wall_flux(ua, ub, x(i) - xg(i), nul, &
                face_id(phi, ibm_obj_id, i+1,j,k, i+1,j,k+1), -1d0)
             ub = open_mean4( U(i-1,j,k), U(i,j,k), U(i-1,j,k+1), U(i,j,k+1), mu(i-1,j,k), mu(i,j,k), mu(i-1,j,k+1), mu(i,j,k+1) )
             If ( j-1 >= 2 ) Then;  If ( mw(i,j-1,k) < 0.5d0 ) &
                ovr_w(i,j,k,3) = wall_flux(ua, ub, yg(j) - y(j-1), nul, face_id(phi, ibm_obj_id, i,j-1,k, i,j-1,k+1), 1d0);  End If
             If ( j+1 <= nyg-1 ) Then;  If ( mw(i,j+1,k) < 0.5d0 ) &
                ovr_w(i,j,k,4) = wall_flux(ua, ub, y(j) - yg(j), nul, face_id(phi, ibm_obj_id, i,j+1,k, i,j+1,k+1), -1d0);  End If
          End Do
       End Do
    End Do

  End Subroutine ibm_stress_update


  !> Loads on the staircase body (ibm_method = 1): pressure on closed faces (linear extrapolation of the two
  !  fluid cells normal to the face) and wall shear: log-law stress if ibm_wall_model_flag = 1, else (nu +
  !  nu_t) u/(h/2); the impulse column is NaN
  Subroutine ibm_stair_forces(U, V, W, Fx_ibm, Fy_ibm, Fz_ibm, Fx_pres, Fy_pres, Fz_pres, Fx_visc, Fy_visc, Fz_visc)

    Real(Int64), Intent(In)  :: U(nx,nyg,nzg), V(nxg,ny,nzg), W(nxg,nyg,nz)
    Real(Int64), Intent(Out) :: Fx_ibm, Fy_ibm, Fz_ibm, Fx_pres, Fy_pres, Fz_pres, Fx_visc, Fy_visc, Fz_visc

    Integer(Int32) :: i, j, k, d, e, e2, sgn, ii, jj, kk, ihi, khi, de(3)
    Real(Int64) :: lpres(3), lvisc(3), gp(3), gv(3), uc(3), area, hn, pface, p2
    Logical :: skip_x, skip_z

    Call dup_cell_flags(skip_x, skip_z)
    ihi = nxg-1;  khi = nzg-1
    If ( skip_x ) ihi = nxg-2
    If ( skip_z ) khi = nzg-2
    lpres = 0d0;  lvisc = 0d0

    Do k = 2, khi
       Do j = 2, nyg-1
          Do i = 2, ihi
             If ( phi(i,j,k) < 0d0 ) Cycle
             uc(1) = 0.5d0*( U(i-1,j,k) + U(i,j,k) )
             uc(2) = 0.5d0*( V(i,j-1,k) + V(i,j,k) )
             uc(3) = 0.5d0*( W(i,j,k-1) + W(i,j,k) )
             Do d = 1, 3
                Do sgn = -1, 1, 2
                   de = 0;  de(d) = sgn
                   ii = i + de(1);  jj = j + de(2);  kk = k + de(3)
                   If ( phi(ii,jj,kk) >= 0d0 ) Cycle
                   If ( d == 1 ) Then
                      area = (y(j)-y(j-1))*(z(k)-z(k-1));  hn = dx
                   Else If ( d == 2 ) Then
                      area = dx*(z(k)-z(k-1));  hn = y(j)-y(j-1)
                   Else
                      area = dx*(y(j)-y(j-1));  hn = z(k)-z(k-1)
                   End If
                   pface = P(i,j,k)
                   ii = i - de(1);  jj = j - de(2);  kk = k - de(3)
                   If ( phi(ii,jj,kk) >= 0d0 ) Then
                      p2 = P(ii,jj,kk);  pface = pface + 0.5d0*( pface - p2 )
                   End If
                   lpres(d) = lpres(d) + Real(sgn,Int64)*pface*area
                   Do e = 1, 3
                      If ( e == d ) Cycle
                      e2 = 6 - d - e
                      If ( ibm_wall_model_flag == 1 ) Then
                         lvisc(e) = lvisc(e) + area*wall_flux( uc(e), uc(e2), 0.5d0*hn, nu, &
                              solid_id(ibm_obj_id, i+de(1), j+de(2), k+de(3)), 1d0 )
                      Else
                         lvisc(e) = lvisc(e) + (nu + nu_t(i,j,k))*uc(e)/(0.5d0*hn)*area
                      End If
                   End Do
                End Do
             End Do
          End Do
       End Do
    End Do

    Fx_ibm = ieee_value(1d0, ieee_quiet_nan);  Fy_ibm = Fx_ibm;  Fz_ibm = Fx_ibm
    Call MPI_Allreduce(lpres, gp, 3, MPI_real8, MPI_SUM, MPI_COMM_WORLD, ierr)
    Call MPI_Allreduce(lvisc, gv, 3, MPI_real8, MPI_SUM, MPI_COMM_WORLD, ierr)
    Fx_pres = gp(1);  Fy_pres = gp(2);  Fz_pres = gp(3)
    Fx_visc = gv(1);  Fy_visc = gv(2);  Fz_visc = gv(3)

  End Subroutine ibm_stair_forces

End Module ibm_stress
