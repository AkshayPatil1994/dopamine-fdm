!> Hydrodynamic loads on immersed bodies in the two-fluid solver (the single-fluid routines in ibm.f90 assume unit density and a
!  kinematic pressure), written to the same ibm_forces.csv: pressure and viscous traction summed over the faces of the staircase
!  body (the faces on which the masked pressure operator has its Neumann condition). The IBM-impulse columns are NaN.
!  The pressure of the two-fluid solver is the dynamic part with the still-water row-mean hydrostatic pressure removed; the load
!  uses the total pressure P + p_s(y), p_s(y) = g * int_y^top rho_ref dy (zero at the top), so a body in still water feels
!  buoyancy. The viscous stress uses the mixture viscosity at the image point plus the SGS part rho * nu_t.
Module vof_ibm

  Use iso_fortran_env, Only : Int32, Int64
  Use global
  Use mpi
  Use halo_pad, Only : pad_field
  Use ieee_arithmetic, Only : ieee_value, ieee_quiet_nan
  Use ibm, Only : U_wall, V_wall, W_wall, dup_cell_flags
  Use vof_state, Only : Cv, hy, vof_rho_ref
  Use ibm_stress, Only : wall_flux, solid_id

  Implicit None

Contains

  !> Hydrostatic reference pressure of the still-water row densities at the cell centres, zero at the top row
  Subroutine hydrostatic_reference(ps)

    Real(Int64), Intent(Out) :: ps(nyg)
    Integer(Int32) :: j
    Real(Int64) :: acc

    ps = 0d0
    If ( vof_hsplit /= 0 ) Return   ! the well-balanced split carries a different reference: loads are then the dynamic part only
    acc = 0d0
    Do j = nyg-1, 2, -1
       ps(j) = vof_grav*( acc + 0.5d0*vof_rho_ref(j)*hy(j) )
       acc = acc + vof_rho_ref(j)*hy(j)
    End Do
    ps(1) = ps(2);  ps(nyg) = 0d0

  End Subroutine hydrostatic_reference


  !> Forces on the bodies. Method 1: IBM impulse over dt (density weighted). Method 2: pressure and viscous traction summed over
  !  the faces of the staircase body (closed faces between a fluid cell, phi >= 0, and a solid cell, phi < 0), the same faces on
  !  which the pressure operator has its Neumann condition, so a uniform pressure gives exactly zero and the still-water
  !  hydrostatic pressure gives the buoyancy of the staircase body. The face pressure is the linear extrapolation of the two
  !  fluid cells normal to the face (exact for the hydrostatic gradient); the shear uses the cell-centre velocity at half a cell,
  !  or the log-law wall stress of the staircase wall model when ibm_wall_model_flag = 1 (the stress the flow receives).
  Subroutine vof_compute_ibm_forces(Fx_ibm, Fy_ibm, Fz_ibm, Fx_pres, Fy_pres, Fz_pres, Fx_visc, Fy_visc, Fz_visc)

    Real(Int64), Intent(Out) :: Fx_ibm, Fy_ibm, Fz_ibm, Fx_pres, Fy_pres, Fz_pres, Fx_visc, Fy_visc, Fz_visc

    Integer(Int32), Parameter :: E2 = 2
    Integer(Int32) :: i, j, k, d, sgn, ii, jj, kk, i2, j2, k2, ihi, khi, e, eo, de(3)
    Real(Int64) :: lpres(3), lvisc(3), mul, mug, area, pF, pF2, pface, mu_c, hn, uc(3), cw(3), rho_c
    Logical :: skip_x, skip_z, ibm_wm
    Real(Int64), Allocatable :: ps(:), Ptot(:,:,:), Pp(:,:,:), Php(:,:,:), Ue(:,:,:), Ve(:,:,:), We(:,:,:)

    mul = vof_rho_l*vof_nu_l;  mug = vof_rho_g*vof_nu_g
    cw = (/ U_wall, V_wall, W_wall /)
    ibm_wm = ( ibm_wall_model_flag == 1 )

    Allocate( ps(nyg), Ptot(nxg,nyg,nzg) )
    Call hydrostatic_reference(ps)
    Do j = 1, nyg
       Ptot(:,j,:) = P(:,j,:) + ps(j)
    End Do
    Allocate( Pp(1-E2:nxg+E2, nyg, 1-E2:nzg+E2), Php(1-E2:nxg+E2, nyg, 1-E2:nzg+E2) )
    Allocate( Ue(1-E2:nx+E2, nyg, 1-E2:nzg+E2), Ve(1-E2:nxg+E2, ny, 1-E2:nzg+E2), We(1-E2:nxg+E2, nyg, 1-E2:nz+E2) )
    Call pad_field( Ptot, nxg, nyg, nzg, .False., .False., E2, Pp )
    Call pad_field( phi,  nxg, nyg, nzg, .False., .False., E2, Php )
    Call pad_field( U, nx,  nyg, nzg, .True.,  .False., E2, Ue )
    Call pad_field( V, nxg, ny,  nzg, .False., .False., E2, Ve )
    Call pad_field( W, nxg, nyg, nz,  .False., .True.,  E2, We )

    Call dup_cell_flags(skip_x, skip_z)
    ihi = nxg-1;  khi = nzg-1
    If ( skip_x ) ihi = nxg-2
    If ( skip_z ) khi = nzg-2
    lpres = 0d0;  lvisc = 0d0

    Do k = 2, khi
       Do j = 2, nyg-1
          Do i = 2, ihi
             If ( Php(i,j,k) < 0d0 ) Cycle
             ! cell-centre velocity of the fluid cell
             uc(1) = 0.5d0*( Ue(i-1,j,k) + Ue(i,j,k) )
             uc(2) = 0.5d0*( Ve(i,j-1,k) + Ve(i,j,k) )
             uc(3) = 0.5d0*( We(i,j,k-1) + We(i,j,k) )
             mu_c = mug + (mul - mug)*Cv(i,j,k) + ( vof_rho_g + (vof_rho_l - vof_rho_g)*Cv(i,j,k) )*nu_t(i,j,k)
             Do d = 1, 3
                Do sgn = -1, 1, 2
                   de = 0;  de(d) = sgn
                   ii = i + de(1);  jj = j + de(2);  kk = k + de(3)
                   If ( jj < 1 .Or. jj > nyg ) Cycle
                   If ( Php(ii,jj,kk) >= 0d0 ) Cycle
                   ! closed face between this fluid cell and a solid cell in direction sgn*d
                   If ( d == 1 ) Then
                      area = (y(j)-y(j-1))*(z(k)-z(k-1));  hn = dx
                   Else If ( d == 2 ) Then
                      area = dx*(z(k)-z(k-1));  hn = y(j)-y(j-1)
                   Else
                      area = dx*(y(j)-y(j-1));  hn = z(k)-z(k-1)
                   End If
                   ! pressure at the face: linear extrapolation from this cell and the next fluid cell away from the body
                   i2 = i - de(1);  j2 = j - de(2);  k2 = k - de(3)
                   pF = Pp(i,j,k)
                   pface = pF
                   If ( j2 >= 1 .And. j2 <= nyg ) Then
                      If ( Php(i2,j2,k2) >= 0d0 ) Then
                         pF2 = Pp(i2,j2,k2)
                         pface = pF + 0.5d0*( pF - pF2 )
                      End If
                   End If
                   lpres(d) = lpres(d) + Real(sgn,Int64)*pface*area
                   ! shear: the tangential velocity of the cell relative to the wall over half a cell
                   rho_c = vof_rho_g + (vof_rho_l - vof_rho_g)*Cv(i,j,k)
                   Do e = 1, 3
                      If ( e == d ) Cycle
                      eo = 6 - d - e
                      If ( ibm_wm ) Then
                         lvisc(e) = lvisc(e) + rho_c*area*wall_flux( uc(e) - cw(e), uc(eo) - cw(eo), 0.5d0*hn, &
                              ( mug + (mul - mug)*Cv(i,j,k) )/rho_c, solid_id(ii,jj,kk), 1d0 )
                      Else
                         lvisc(e) = lvisc(e) + mu_c*( uc(e) - cw(e) )/(0.5d0*hn)*area
                      End If
                   End Do
                End Do
             End Do
          End Do
       End Do
    End Do

    ! Method 1 (the IBM impulse) is not available: the pressure acts through the closed faces of the masked operator, so the
    ! momentum the ghost-cell IBM gives to the fluid no longer carries the load. NaN in the CSV rather than a misleading zero.
    Fx_ibm = ieee_value(1d0, ieee_quiet_nan);  Fy_ibm = Fx_ibm;  Fz_ibm = Fx_ibm
    Call MPI_Allreduce(lpres(1), Fx_pres, 1, MPI_real8, MPI_SUM, MPI_COMM_WORLD, ierr)
    Call MPI_Allreduce(lpres(2), Fy_pres, 1, MPI_real8, MPI_SUM, MPI_COMM_WORLD, ierr)
    Call MPI_Allreduce(lpres(3), Fz_pres, 1, MPI_real8, MPI_SUM, MPI_COMM_WORLD, ierr)
    Call MPI_Allreduce(lvisc(1), Fx_visc, 1, MPI_real8, MPI_SUM, MPI_COMM_WORLD, ierr)
    Call MPI_Allreduce(lvisc(2), Fy_visc, 1, MPI_real8, MPI_SUM, MPI_COMM_WORLD, ierr)
    Call MPI_Allreduce(lvisc(3), Fz_visc, 1, MPI_real8, MPI_SUM, MPI_COMM_WORLD, ierr)

    Deallocate( ps, Ptot, Pp, Php, Ue, Ve, We )

  End Subroutine vof_compute_ibm_forces

End Module vof_ibm
