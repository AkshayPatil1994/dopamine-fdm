!     Wall-model stress on the staircase faces of an immersed body
Module ibm_stress

  Use iso_fortran_env, Only : Int32, Int64
  Use global
  Use wallmodel, Only : solve_u_tau_reichardt, solve_u_tau_rough

  Implicit None

  ! Kinematic momentum flux through an edge between a fluid face and a closed face, replacing the viscous flux the stencil carries there.
  ! Last index: U faces (y-lo, y-hi, z-lo, z-hi), V faces (x-lo, x-hi, z-lo, z-hi), W faces (x-lo, x-hi, y-lo, y-hi); ovr_none = no override
  Real(Int64), Parameter :: ovr_none = Huge(1d0)
  Real(Int64), Allocatable, Dimension(:,:,:,:) :: ovr_u, ovr_v, ovr_w

Contains

  Subroutine ibm_stress_init

    Allocate( ovr_u(nx,nyg,nzg,4), ovr_v(nxg,ny,nzg,4), ovr_w(nxg,nyg,nz,4) )
    ovr_u = ovr_none;  ovr_v = ovr_none;  ovr_w = ovr_none

  End Subroutine ibm_stress_init


  !> Wall shear of the log law at distance d from the closed face, along the tangential velocity (ua along the component, ub the other
  !  tangential one); sg = +1 when the closed face is on the low side of the fluid face, -1 on the high side
  Function wall_flux(ua, ub, d, nul, oid, sg) Result(f)

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


  Pure Function solid_id(i, j, k) Result(oid)

    Integer(Int32), Intent(In) :: i, j, k
    Integer(Int32) :: oid

    oid = 0
    If ( Allocated(ibm_obj_id) ) oid = Min(Max(Nint(ibm_obj_id(i,j,k)), 0), max_ibm_objects)

  End Function solid_id


  !> Set the flux overrides from the current velocity. mu, mv, mw are the 1 (open) / 0 (closed) face masks; nu_cell is the kinematic
  !  viscosity at the cell centres, <= 0 where the wall model is off (the stencil flux stays). Domain-wall rows are left to the flat-wall model.
  Subroutine ibm_stress_update(U, V, W, mu, mv, mw, nu_cell)

    Real(Int64), Intent(In) :: U(nx,nyg,nzg), V(nxg,ny,nzg), W(nxg,nyg,nz)
    Real(Int64), Intent(In) :: mu(nx,nyg,nzg), mv(nxg,ny,nzg), mw(nxg,nyg,nz), nu_cell(nxg,nyg,nzg)

    Integer(Int32) :: i, j, k
    Real(Int64) :: nul, ua, ub

    ovr_u = ovr_none;  ovr_v = ovr_none;  ovr_w = ovr_none

    Do k = 2, nzg-1
       Do j = 2, nyg-1
          Do i = 2, nx-1
             If ( mu(i,j,k) < 0.5d0 ) Cycle
             nul = Min(nu_cell(i,j,k), nu_cell(i+1,j,k))
             If ( nul <= 0d0 ) Cycle
             ua = U(i,j,k)
             ub = 0.25d0*( W(i,j,k-1) + W(i,j,k) + W(i+1,j,k-1) + W(i+1,j,k) )
             If ( j-1 >= 2 ) Then;  If ( mu(i,j-1,k) < 0.5d0 ) &
                ovr_u(i,j,k,1) = wall_flux(ua, ub, yg(j) - y(j-1), nul, solid_id(i,j-1,k), 1d0);  End If
             If ( j+1 <= nyg-1 ) Then;  If ( mu(i,j+1,k) < 0.5d0 ) &
                ovr_u(i,j,k,2) = wall_flux(ua, ub, y(j) - yg(j), nul, solid_id(i,j+1,k), -1d0);  End If
             ub = 0.25d0*( V(i,j-1,k) + V(i,j,k) + V(i+1,j-1,k) + V(i+1,j,k) )
             If ( mu(i,j,k-1) < 0.5d0 ) ovr_u(i,j,k,3) = wall_flux(ua, ub, zg(k) - z(k-1), nul, solid_id(i,j,k-1), 1d0)
             If ( mu(i,j,k+1) < 0.5d0 ) ovr_u(i,j,k,4) = wall_flux(ua, ub, z(k) - zg(k), nul, solid_id(i,j,k+1), -1d0)
          End Do
       End Do
    End Do

    Do k = 2, nzg-1
       Do j = 2, ny-1
          Do i = 2, nxg-1
             If ( mv(i,j,k) < 0.5d0 ) Cycle
             nul = Min(nu_cell(i,j,k), nu_cell(i,j+1,k))
             If ( nul <= 0d0 ) Cycle
             ua = V(i,j,k)
             ub = 0.25d0*( W(i,j,k-1) + W(i,j,k) + W(i,j+1,k-1) + W(i,j+1,k) )
             If ( mv(i-1,j,k) < 0.5d0 ) ovr_v(i,j,k,1) = wall_flux(ua, ub, xg(i) - x(i-1), nul, solid_id(i-1,j,k), 1d0)
             If ( mv(i+1,j,k) < 0.5d0 ) ovr_v(i,j,k,2) = wall_flux(ua, ub, x(i) - xg(i), nul, solid_id(i+1,j,k), -1d0)
             ub = 0.25d0*( U(i-1,j,k) + U(i,j,k) + U(i-1,j+1,k) + U(i,j+1,k) )
             If ( mv(i,j,k-1) < 0.5d0 ) ovr_v(i,j,k,3) = wall_flux(ua, ub, zg(k) - z(k-1), nul, solid_id(i,j,k-1), 1d0)
             If ( mv(i,j,k+1) < 0.5d0 ) ovr_v(i,j,k,4) = wall_flux(ua, ub, z(k) - zg(k), nul, solid_id(i,j,k+1), -1d0)
          End Do
       End Do
    End Do

    Do k = 2, nz-1
       Do j = 2, nyg-1
          Do i = 2, nxg-1
             If ( mw(i,j,k) < 0.5d0 ) Cycle
             nul = Min(nu_cell(i,j,k), nu_cell(i,j,k+1))
             If ( nul <= 0d0 ) Cycle
             ua = W(i,j,k)
             ub = 0.25d0*( V(i,j-1,k) + V(i,j,k) + V(i,j-1,k+1) + V(i,j,k+1) )
             If ( mw(i-1,j,k) < 0.5d0 ) ovr_w(i,j,k,1) = wall_flux(ua, ub, xg(i) - x(i-1), nul, solid_id(i-1,j,k), 1d0)
             If ( mw(i+1,j,k) < 0.5d0 ) ovr_w(i,j,k,2) = wall_flux(ua, ub, x(i) - xg(i), nul, solid_id(i+1,j,k), -1d0)
             ub = 0.25d0*( U(i-1,j,k) + U(i,j,k) + U(i-1,j,k+1) + U(i,j,k+1) )
             If ( j-1 >= 2 ) Then;  If ( mw(i,j-1,k) < 0.5d0 ) &
                ovr_w(i,j,k,3) = wall_flux(ua, ub, yg(j) - y(j-1), nul, solid_id(i,j-1,k), 1d0);  End If
             If ( j+1 <= nyg-1 ) Then;  If ( mw(i,j+1,k) < 0.5d0 ) &
                ovr_w(i,j,k,4) = wall_flux(ua, ub, y(j) - yg(j), nul, solid_id(i,j+1,k), -1d0);  End If
          End Do
       End Do
    End Do

  End Subroutine ibm_stress_update

End Module ibm_stress
