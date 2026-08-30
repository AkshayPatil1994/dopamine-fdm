!> Boussinesq temperature transport (advection-diffusion) and its boundary conditions
Module thermal_transport

  Use iso_fortran_env,  Only : Int32, Int64
  Use global
  Use mpi
  Use boundary_conditions, Only : apply_periodic_bc_z, apply_inflow_bc_scalar_x, outflow_convection_velocity
  Use scalar_transport,    Only : compute_rhs_scalar_core, update_ghost_scalar

  Implicit None

Contains

  !> Compute RHS for temperature advection-diffusion; molecular+turbulent diffusivity from Pr/Pr_t, no settling term
  Subroutine compute_rhs_temperature(T_, U_, V_, W_, Ft_)

    Real(Int64), Dimension(nxg, nyg, nzg), Intent(In)    :: T_
    Real(Int64), Dimension(nx,  nyg, nzg), Intent(In)    :: U_
    Real(Int64), Dimension(nxg, ny,  nzg), Intent(In)    :: V_
    Real(Int64), Dimension(nxg, nyg, nz ), Intent(In)    :: W_
    Real(Int64), Dimension(2:nxg-1, 2:nyg-1, 2:nzg-1), Intent(Out) :: Ft_

    Call compute_rhs_scalar_core(T_, U_, V_, W_, 0d0, nu/Pr, 1d0/Pr_t, Ft_)

  End Subroutine compute_rhs_temperature


  !> Apply boundary conditions to temperature T (x-periodic or SEM-inflow/convective-outflow, z-MPI halo, y-wall adiabatic/isothermal); runs host-side (caller syncs Tscal host/device around this call)
  Subroutine apply_temperature_bc(T_)

    Real(Int64), Dimension(nxg, nyg, nzg), Intent(InOut) :: T_

    Real(Int64) :: Uc, courant

    ! x direction: periodic, or Dirichlet-SEM inflow / convective outflow
    If ( x_bc_type == 0 ) Then
       T_(1,   :,:) = T_(nxg-2,:,:)
       T_(nxg-1,:,:) = T_(2,   :,:)
       T_(nxg,  :,:) = T_(3,   :,:)
    Else
       Call apply_inflow_bc_scalar_x(T_)
       Uc = outflow_convection_velocity()   ! scalar-only, safe to call from host code (see boundary_conditions.f90)
       courant = Min(Max(Uc,0d0)*dt/dx, 1d0)
       T_(nxg,:,:) = T_(nxg,:,:) - courant*( T_(nxg,:,:) - T_(nxg-1,:,:) )
    End If

    ! z-halo via MPI (ring exchange, non-periodic) + z-periodic wrap
    Call update_ghost_scalar(T_)
    ! Push host state to device first: apply_periodic_bc_z runs device-resident at nprocs==1, else its z-wrap fill is clobbered by the caller's later blanket update device
    !$acc update device(T_)
    Call apply_periodic_bc_z(T_, 4)
    !$acc update host(T_)

    ! y-bottom ghost: 0=adiabatic (zero-gradient), 1=isothermal (Dirichlet mirror)
    If ( T_bc_bot == 0 ) Then
       T_(:,1,:) = T_(:,2,:)
    Else
       T_(:,1,:) = 2d0*T_wall_bot - T_(:,2,:)
    End If

    ! y-top ghost
    If ( T_bc_top == 0 ) Then
       T_(:,nyg,:) = T_(:,nyg-1,:)
    Else
       T_(:,nyg,:) = 2d0*T_wall_top - T_(:,nyg-1,:)
    End If

    ! Fallback safety net inside IBM solid cells; apply_ghost_cell_ibm_scalar (ibm.f90) enforces the actual physical wall condition
    If ( ibm_input_mode >= 1 .And. Allocated(phi) ) Then
       Where ( phi(2:nxg-1, 2:nyg-1, 2:nzg-1) <= 0d0 )
          T_(2:nxg-1, 2:nyg-1, 2:nzg-1) = T_ref
       End Where
    End If

  End Subroutine apply_temperature_bc


End Module thermal_transport
