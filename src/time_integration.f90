! Module for temporal integration
Module time_integration

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use mpi
  Use equations
  Use projection
  Use boundary_conditions
  Use wallmodel
  Use sgs_models
  Use ibm
  Use scalar_transport,  Only : compute_rhs_scalar, apply_scalar_bc
  Use thermal_transport, Only : compute_rhs_temperature, apply_temperature_bc
  Use monitor,          Only : compute_cfl, write_force_csv, compute_bulk_velocity
  Use profiler

  ! prevent implicit typing
  Implicit None

  ! Persistent velocity snapshot buffers for IBM force computation.
  ! Allocated once on the first qualifying step; never deallocated until finalise.
  Real(Int64), Allocatable :: U_pre(:,:,:), V_pre(:,:,:), W_pre(:,:,:)

Contains

  !> Update imposed pressure gradients for oscillatory/pulsatile forcing: dPdx=dPdx_t+Ub_x*waveOmega_x*cos(waveOmega_x*t+phi_wave_x), analogous in z
  Subroutine update_pressure_forcing

    If ( waveOmega_x /= 0d0 ) Then
       dPdx = dPdx_t + Ub_x * waveOmega_x * Cos(waveOmega_x * t + phi_wave_x)
    Else
       dPdx = dPdx_t
    End If

    If ( waveOmega_z /= 0d0 ) Then
       dPdz = dPdz_t + Ub_z * waveOmega_z * Cos(waveOmega_z * t + phi_wave_z)
    Else
       dPdz = dPdz_t
    End If

  End Subroutine update_pressure_forcing

  ! Explicit Runge-Kutta 3 steps
  Subroutine compute_time_step_RK3

    Real(Int64) :: to
    Real(Int64) :: Fx_ibm,  Fy_ibm,  Fz_ibm
    Real(Int64) :: Fx_pres, Fy_pres, Fz_pres
    Real(Int64) :: Fx_visc, Fy_visc, Fz_visc
    Real(Int64) :: cfl_conv, cfl_visc, dt_new, dt_presnap
    Real(Int64) :: Ub_now, dU_cmfr
    Logical     :: needs_final_sync

    ! Enforce IBM before saving old state (zeroes solid cells on step 1, no-op thereafter); apply_ghost_cell_ibm is device-resident and U,V,W are already device-current, so no sync needed
    If ( ibm_input_mode >= 1 ) Then
       Call profiler_start(PROF_IBM)
       Call apply_ghost_cell_ibm(U, V, W)
       ! compute_sgs_model's Pass 2 below is host-only when IBM is active and needs current U,V,W
       !$acc update host(U,V,W)
       Call profiler_stop(PROF_IBM)
    End If

    ! Compute SGS eddy viscosity before the CFL check so viscous CFL uses current-step nu_t, not the previous step's (DNS no-op when sgs_model==0)
    ! U,V,W already device-resident (see above); compute_sgs_model leaves nu_t device-resident too
    Call profiler_start(PROF_SGS)
    Call compute_sgs_model(U,V,W,nu_t)
    Call profiler_stop(PROF_SGS)

    ! CFL check and adaptive dt, evaluated from the start-of-step velocity so dt is enforced on the current step, not lagged by one
    Call profiler_start(PROF_CFL)
    Call compute_cfl(cfl_conv, cfl_visc)
    Call profiler_stop(PROF_CFL)
    cfl_conv_last = cfl_conv
    cfl_visc_last = cfl_visc
    cfl_current   = Max(cfl_conv, cfl_visc)
    If ( cfl_adaptive == 1 .And. cfl_current > 0d0 ) Then
       dt_new = dt * cfl_target / cfl_current * cfl_safety
       dt = Max(dt_min, Min(dt_max, dt_new))
    End If

    ! nsave<0: shrink dt (even below dt_min) so t lands exactly on tsave_next, keeping physical-time-based snapshots uniformly spaced under adaptive dt
    dt_presnap = dt
    If ( nsave < 0 .And. t + dt > tsave_next ) Then
       dt = tsave_next - t
    End If

    ! save previous state
    to = t
    Call profiler_start(PROF_RK_UPDATE)
    ! U,V,W are already device-resident; copy directly on-device instead of bouncing through host
    !$acc kernels present(U,V,W,Uo,Vo,Wo)
    Uo = U
    Vo = V
    Wo = W
    !$acc end kernels
    Call profiler_stop(PROF_RK_UPDATE)
    If ( sediment_flag >= 1 ) Then
       !$acc kernels present(Cscal,Cscal_o)
       Cscal_o = Cscal
       !$acc end kernels
    End If
    If ( boussinesq_flag >= 1 ) Then
       !$acc kernels present(Tscal,Tscal_o)
       Tscal_o = Tscal
       !$acc end kernels
    End If

    ! Update imposed pressure gradient (oscillatory / steady forcing); under constant-mass-flux forcing dPdx carries no direct RHS forcing (the CMFR correction below is the sole forcing mechanism)
    If ( flow_forcing_mode == 0 ) Then
       Call update_pressure_forcing
    Else
       dPdx = 0d0
    End If

    ! step 1
    rk_step = 1
    ! nu_t, U,V,W stay device-resident throughout; compute_wall_model is device-resident too (except the host-only compute_pseudo_pressure_bc_for_robin_bc), so no update device needed here
    Call profiler_start(PROF_WALLMODEL)
    Call compute_wall_model(U,V,W,nu_t) ! sets alpha + EQWM ghost cells on old velocity
    Call profiler_stop(PROF_WALLMODEL)
    Call profiler_start(PROF_RHS)
    Call compute_rhs_u(U,V,W,Fu1)
    Call compute_rhs_v(U,V,W,Fv1)
    Call compute_rhs_w(U,V,W,Fw1)
    Call profiler_stop(PROF_RHS)

    ! RK-stage velocity update, GPU-resident (Fu1/Fw1 never leave device; Fv1 round-trips to host inside apply_uav_forcing when uav_active>=1)
    Call profiler_start(PROF_RK_UPDATE)
    !$acc kernels present(U,V,W,Uo,Vo,Wo,Fu1,Fv1,Fw1)
    U(2:nx-1,2:nyg-1,2:nzg-1) = Uo(2:nx-1,2:nyg-1,2:nzg-1) + dt*rk_coef(1,1)*Fu1
    V(2:nxg-1,2:ny-1,2:nzg-1) = Vo(2:nxg-1,2:ny-1,2:nzg-1) + dt*rk_coef(1,1)*Fv1
    W(2:nxg-1,2:nyg-1,2:nz-1) = Wo(2:nxg-1,2:nyg-1,2:nz-1) + dt*rk_coef(1,1)*Fw1
    !$acc end kernels
    Call profiler_stop(PROF_RK_UPDATE)

    ! Apply ghost-cell IBM on updated velocity (2nd-order; no-slip or EQWM)
    ! Always zero solid cells first; then overwrite ghost cells with EQWM if active.
    If ( ibm_input_mode >= 1 ) Then
       Call profiler_start(PROF_IBM)
       ! Stage 1: zero accumulators, capture first IBM impulse; only the host-only impulse snapshot needs an explicit sync, when active
       If ( nsampling > 0 .And. Mod(istep, nsampling) == 0 ) Then
          If ( .Not. Allocated(U_pre) ) Then
             Allocate( U_pre(nx, nyg,nzg), V_pre(nxg, ny,nzg), W_pre(nxg,nyg, nz) )
          End If
          Call zero_ibm_stage_accumulators()
          !$acc update host(U,V,W)
          U_pre = U;  V_pre = V;  W_pre = W
       End If
       Call apply_ghost_cell_ibm(U, V, W)
       If ( nsampling > 0 .And. Mod(istep, nsampling) == 0 ) Then
          ! apply_ghost_cell_ibm runs device-resident; sync before this host-only impulse read
          !$acc update host(U,V,W)
          Call accumulate_ibm_stage_impulse(U_pre, V_pre, W_pre, U, V, W)
       End If
       If ( ibm_wall_model_flag == 1 ) Then
          Call compute_ibm_wall_model(U, V, W, nu_t)
       End If
       ! No update device needed: apply_ghost_cell_ibm/compute_ibm_wall_model both write device
       ! directly now, and the host-only reads above never modify U,V,W (Intent(In) throughout)
       Call profiler_stop(PROF_IBM)
    End If

    t = to + rk_t(rk_step)*dt

    Call profiler_start(PROF_BC)
    Call apply_boundary_conditions
    Call profiler_stop(PROF_BC)
    Call compute_projection_step
    Call profiler_start(PROF_BC)
    Call apply_boundary_conditions(after_projection=.True.)
    Call profiler_stop(PROF_BC)
    ! U,V,W now correct+resident on device; only the IBM re-enforce block below (if active) needs a host mirror here
    If ( ibm_input_mode >= 1 ) Then
       Call profiler_start(PROF_RK_UPDATE)
       !$acc update host(U,V,W)
       Call profiler_stop(PROF_RK_UPDATE)
       ! Re-enforce IBM: projection corrupts ghost-cell velocities via grad(p) correction.
       ! Accumulate the grad(P)-correction impulse so Method 1 captures the full force.
       Call profiler_start(PROF_IBM)
       If ( nsampling > 0 .And. Mod(istep, nsampling) == 0 ) Then
          U_pre = U;  V_pre = V;  W_pre = W
       End If
       Call apply_ghost_cell_ibm(U, V, W)
       If ( nsampling > 0 .And. Mod(istep, nsampling) == 0 ) Then
          ! apply_ghost_cell_ibm runs device-resident; sync before this host-only impulse read
          !$acc update host(U,V,W)
          Call accumulate_ibm_stage_impulse(U_pre, V_pre, W_pre, U, V, W)
       End If
       If ( ibm_wall_model_flag == 1 ) Call compute_ibm_wall_model(U, V, W, nu_t)
       ! Refresh host mirror: next substage's SGS Pass 2 reads U,V,W host-only; see GPU_porting.md
       !$acc update host(U,V,W)
       Call profiler_stop(PROF_IBM)
    End If

    ! Scalar step 1
    If ( sediment_flag >= 1 ) Then
       Call profiler_start(PROF_SCALAR)
       Call compute_rhs_scalar(Cscal, U, V, W, Fcs1)
       !$acc kernels present(Cscal,Cscal_o,Fcs1)
       Cscal(2:nxg-1,2:nyg-1,2:nzg-1) = Cscal_o(2:nxg-1,2:nyg-1,2:nzg-1) + dt*rk_coef(1,1)*Fcs1
       !$acc end kernels
       !$acc update host(Cscal)
       Call apply_scalar_bc(Cscal)
       !$acc update device(Cscal)
       Call profiler_stop(PROF_SCALAR)
    End If

    ! Temperature step 1 (buoyancy in stage n's compute_rhs_v reads Tscal as finalized at the end of stage n-1)
    If ( boussinesq_flag >= 1 ) Then
       Call profiler_start(PROF_SCALAR)
       Call compute_rhs_temperature(Tscal, U, V, W, Ft1)
       !$acc kernels present(Tscal,Tscal_o,Ft1)
       Tscal(2:nxg-1,2:nyg-1,2:nzg-1) = Tscal_o(2:nxg-1,2:nyg-1,2:nzg-1) + dt*rk_coef(1,1)*Ft1
       !$acc end kernels
       !$acc update host(Tscal)
       Call apply_temperature_bc(Tscal)
       !$acc update device(Tscal)
       If ( ibm_input_mode >= 1 ) Call apply_ghost_cell_ibm_scalar(Tscal)
       Call profiler_stop(PROF_SCALAR)
    End If

    ! step 2
    rk_step = 2
    Call profiler_start(PROF_SGS)
    ! U,V,W already device-resident; nothing in this design ever writes them from the host
    Call compute_sgs_model(U,V,W,nu_t)
    Call profiler_stop(PROF_SGS)
    Call profiler_start(PROF_WALLMODEL)
    Call compute_wall_model(U,V,W,nu_t)
    Call profiler_stop(PROF_WALLMODEL)
    Call profiler_start(PROF_RHS)
    Call compute_rhs_u(U,V,W,Fu2)
    Call compute_rhs_v(U,V,W,Fv2)
    Call compute_rhs_w(U,V,W,Fw2)
    Call profiler_stop(PROF_RHS)

    Call profiler_start(PROF_RK_UPDATE)
    !$acc kernels present(U,V,W,Uo,Vo,Wo,Fu1,Fv1,Fw1,Fu2,Fv2,Fw2)
    U(2:nx-1,2:nyg-1,2:nzg-1) = Uo(2:nx-1,2:nyg-1,2:nzg-1) + dt*( rk_coef(2,1)*Fu1 + rk_coef(2,2)*Fu2 )
    V(2:nxg-1,2:ny-1,2:nzg-1) = Vo(2:nxg-1,2:ny-1,2:nzg-1) + dt*( rk_coef(2,1)*Fv1 + rk_coef(2,2)*Fv2 )
    W(2:nxg-1,2:nyg-1,2:nz-1) = Wo(2:nxg-1,2:nyg-1,2:nz-1) + dt*( rk_coef(2,1)*Fw1 + rk_coef(2,2)*Fw2 )
    !$acc end kernels
    Call profiler_stop(PROF_RK_UPDATE)

    If ( ibm_input_mode >= 1 ) Then
       Call profiler_start(PROF_IBM)
       ! Stage 2: capture IBM impulse; only the host-only impulse snapshot below needs an explicit sync, when active
       If ( nsampling > 0 .And. Mod(istep, nsampling) == 0 ) Then
          !$acc update host(U,V,W)
          U_pre = U;  V_pre = V;  W_pre = W
       End If
       Call apply_ghost_cell_ibm(U, V, W)
       If ( nsampling > 0 .And. Mod(istep, nsampling) == 0 ) Then
          ! apply_ghost_cell_ibm runs device-resident; sync before this host-only impulse read
          !$acc update host(U,V,W)
          Call accumulate_ibm_stage_impulse(U_pre, V_pre, W_pre, U, V, W)
       End If
       If ( ibm_wall_model_flag == 1 ) Then
          Call compute_ibm_wall_model(U, V, W, nu_t)
       End If
       ! No update device needed: apply_ghost_cell_ibm/compute_ibm_wall_model both write device
       ! directly now, and the host-only reads above never modify U,V,W (Intent(In) throughout)
       Call profiler_stop(PROF_IBM)
    End If

    t = to + rk_t(rk_step)*dt

    Call profiler_start(PROF_BC)
    Call apply_boundary_conditions
    Call profiler_stop(PROF_BC)
    Call compute_projection_step
    Call profiler_start(PROF_BC)
    Call apply_boundary_conditions(after_projection=.True.)
    Call profiler_stop(PROF_BC)
    ! U,V,W now correct+resident on device; only the IBM re-enforce block below (if active) needs a host mirror here
    If ( ibm_input_mode >= 1 ) Then
       Call profiler_start(PROF_RK_UPDATE)
       !$acc update host(U,V,W)
       Call profiler_stop(PROF_RK_UPDATE)
       ! Re-enforce IBM: projection corrupts ghost-cell velocities via grad(p) correction.
       ! Accumulate the grad(P)-correction impulse so Method 1 captures the full force.
       Call profiler_start(PROF_IBM)
       If ( nsampling > 0 .And. Mod(istep, nsampling) == 0 ) Then
          U_pre = U;  V_pre = V;  W_pre = W
       End If
       Call apply_ghost_cell_ibm(U, V, W)
       If ( nsampling > 0 .And. Mod(istep, nsampling) == 0 ) Then
          ! apply_ghost_cell_ibm runs device-resident; sync before this host-only impulse read
          !$acc update host(U,V,W)
          Call accumulate_ibm_stage_impulse(U_pre, V_pre, W_pre, U, V, W)
       End If
       If ( ibm_wall_model_flag == 1 ) Call compute_ibm_wall_model(U, V, W, nu_t)
       ! Refresh host mirror: see the matching comment on substage 1's re-enforce-IBM block above.
       !$acc update host(U,V,W)
       Call profiler_stop(PROF_IBM)
    End If

    ! Scalar step 2
    If ( sediment_flag >= 1 ) Then
       Call profiler_start(PROF_SCALAR)
       Call compute_rhs_scalar(Cscal, U, V, W, Fcs2)
       !$acc kernels present(Cscal,Cscal_o,Fcs1,Fcs2)
       Cscal(2:nxg-1,2:nyg-1,2:nzg-1) = Cscal_o(2:nxg-1,2:nyg-1,2:nzg-1) + &
            dt*( rk_coef(2,1)*Fcs1 + rk_coef(2,2)*Fcs2 )
       !$acc end kernels
       !$acc update host(Cscal)
       Call apply_scalar_bc(Cscal)
       !$acc update device(Cscal)
       Call profiler_stop(PROF_SCALAR)
    End If

    ! Temperature step 2
    If ( boussinesq_flag >= 1 ) Then
       Call profiler_start(PROF_SCALAR)
       Call compute_rhs_temperature(Tscal, U, V, W, Ft2)
       !$acc kernels present(Tscal,Tscal_o,Ft1,Ft2)
       Tscal(2:nxg-1,2:nyg-1,2:nzg-1) = Tscal_o(2:nxg-1,2:nyg-1,2:nzg-1) + &
            dt*( rk_coef(2,1)*Ft1 + rk_coef(2,2)*Ft2 )
       !$acc end kernels
       !$acc update host(Tscal)
       Call apply_temperature_bc(Tscal)
       !$acc update device(Tscal)
       If ( ibm_input_mode >= 1 ) Call apply_ghost_cell_ibm_scalar(Tscal)
       Call profiler_stop(PROF_SCALAR)
    End If

    ! step 3
    rk_step = 3
    Call profiler_start(PROF_SGS)
    ! U,V,W already device-resident; nothing in this design ever writes them from the host
    Call compute_sgs_model(U,V,W,nu_t)
    Call profiler_stop(PROF_SGS)
    Call profiler_start(PROF_WALLMODEL)
    Call compute_wall_model(U,V,W,nu_t)
    Call profiler_stop(PROF_WALLMODEL)
    Call profiler_start(PROF_RHS)
    Call compute_rhs_u(U,V,W,Fu3)
    Call compute_rhs_v(U,V,W,Fv3)
    Call compute_rhs_w(U,V,W,Fw3)
    Call profiler_stop(PROF_RHS)

    Call profiler_start(PROF_RK_UPDATE)
    !$acc kernels present(U,V,W,Uo,Vo,Wo,Fu1,Fv1,Fw1,Fu2,Fv2,Fw2,Fu3,Fv3,Fw3)
    U(2:nx-1,2:nyg-1,2:nzg-1) = Uo(2:nx-1,2:nyg-1,2:nzg-1) + &
         dt*( rk_coef(3,1)*Fu1 + rk_coef(3,2)*Fu2 + rk_coef(3,3)*Fu3 )
    V(2:nxg-1,2:ny-1,2:nzg-1) = Vo(2:nxg-1,2:ny-1,2:nzg-1) + &
         dt*( rk_coef(3,1)*Fv1 + rk_coef(3,2)*Fv2 + rk_coef(3,3)*Fv3 )
    W(2:nxg-1,2:nyg-1,2:nz-1) = Wo(2:nxg-1,2:nyg-1,2:nz-1) + &
         dt*( rk_coef(3,1)*Fw1 + rk_coef(3,2)*Fw2 + rk_coef(3,3)*Fw3 )
    !$acc end kernels
    Call profiler_stop(PROF_RK_UPDATE)

    If ( ibm_input_mode >= 1 ) Then
       Call profiler_start(PROF_IBM)
       ! Stage 3: capture final IBM impulse (accumulators zeroed at stage 1); only the host-only impulse snapshot below needs an explicit sync
       If ( nsampling > 0 .And. Mod(istep, nsampling) == 0 ) Then
          !$acc update host(U,V,W)
          U_pre = U;  V_pre = V;  W_pre = W
       End If
       Call apply_ghost_cell_ibm(U, V, W)
       If ( nsampling > 0 .And. Mod(istep, nsampling) == 0 ) Then
          ! apply_ghost_cell_ibm runs device-resident; sync before this host-only impulse read
          !$acc update host(U,V,W)
          Call accumulate_ibm_stage_impulse(U_pre, V_pre, W_pre, U, V, W)
       End If
       If ( ibm_wall_model_flag == 1 ) Then
          Call compute_ibm_wall_model(U, V, W, nu_t)
       End If
       ! No update device needed: apply_ghost_cell_ibm/compute_ibm_wall_model both write device
       ! directly now, and the host-only reads above never modify U,V,W (Intent(In) throughout)
       Call profiler_stop(PROF_IBM)
    End If

    t = to + rk_t(rk_step)*dt

    Call profiler_start(PROF_BC)
    Call apply_boundary_conditions
    Call profiler_stop(PROF_BC)
    Call compute_projection_step
    Call profiler_start(PROF_BC)
    Call apply_boundary_conditions(after_projection=.True.)
    Call profiler_stop(PROF_BC)

    ! Constant-mass-flux forcing: shift U by a uniform constant so the volume-averaged bulk velocity hits Ub_target exactly (div-free preserved since the shift is spatially uniform); dPdx tracks the equivalent forcing for diagnostics
    If ( flow_forcing_mode == 1 ) Then
       Call profiler_start(PROF_CMFR)
       Call compute_bulk_velocity(Ub_now)
       dU_cmfr = Ub_target - Ub_now
       !$acc kernels present(U)
       U(2:nx-1,2:nyg-1,2:nzg-1) = U(2:nx-1,2:nyg-1,2:nzg-1) + dU_cmfr
       !$acc end kernels
       dPdx_cmfr = -dU_cmfr/dt   ! diagnostic only; dPdx itself stays 0 (see update_pressure_forcing guard above)
       Call profiler_stop(PROF_CMFR)
       Call profiler_start(PROF_BC)
       Call apply_boundary_conditions(after_projection=.True.)
       Call profiler_stop(PROF_BC)
    End If

    ! Final sync: host U,V,W only needed this step if a host-only consumer will actually run (IBM re-enforce below, RSB/TI-rescale accumulation, monitor/divergence check, a field snapshot, or a slice/line probe)
    needs_final_sync = ( ibm_input_mode >= 1 ) .Or. &
         ( rsb_active == 1 .And. istep >= rsb_nstart ) .Or. &
         ( ti_rescale_active == 1 .And. istep >= ti_rescale_nstart ) .Or. &
         ( Mod(istep, nmonitor) == 0 ) .Or. &
         ( nsave > 0 .And. Mod(istep, nsave) == 0 ) .Or. &
         ( nsave < 0 .And. t >= tsave_next - 1d-10 ) .Or. &
         ( n_slices > 0 .And. slice_freq > 0 .And. Mod(istep, slice_freq) == 0 ) .Or. &
         ( n_lines  > 0 .And. line_freq  > 0 .And. Mod(istep, line_freq)  == 0 )
    If ( needs_final_sync ) Then
       Call profiler_start(PROF_RK_UPDATE)
       !$acc update host(U,V,W)
       Call profiler_stop(PROF_RK_UPDATE)
    End If
    ! Re-enforce IBM: projection corrupts ghost-cell velocities via grad(p) correction.
    ! Accumulate the grad(P)-correction impulse so Method 1 captures the full force.
    If ( ibm_input_mode >= 1 ) Then
       Call profiler_start(PROF_IBM)
       If ( nsampling > 0 .And. Mod(istep, nsampling) == 0 ) Then
          U_pre = U;  V_pre = V;  W_pre = W
       End If
       Call apply_ghost_cell_ibm(U, V, W)
       If ( nsampling > 0 .And. Mod(istep, nsampling) == 0 ) Then
          ! apply_ghost_cell_ibm runs device-resident; sync before this host-only impulse read
          !$acc update host(U,V,W)
          Call accumulate_ibm_stage_impulse(U_pre, V_pre, W_pre, U, V, W)
       End If
       If ( ibm_wall_model_flag == 1 ) Call compute_ibm_wall_model(U, V, W, nu_t)
       ! Refresh host mirror: output_monitor/output_data/check_divergence read U,V,W host-only
       !$acc update host(U,V,W)
       Call profiler_stop(PROF_IBM)
    End If

    ! Scalar step 3
    If ( sediment_flag >= 1 ) Then
       Call profiler_start(PROF_SCALAR)
       Call compute_rhs_scalar(Cscal, U, V, W, Fcs3)
       !$acc kernels present(Cscal,Cscal_o,Fcs1,Fcs2,Fcs3)
       Cscal(2:nxg-1,2:nyg-1,2:nzg-1) = Cscal_o(2:nxg-1,2:nyg-1,2:nzg-1) + &
            dt*( rk_coef(3,1)*Fcs1 + rk_coef(3,2)*Fcs2 + rk_coef(3,3)*Fcs3 )
       !$acc end kernels
       !$acc update host(Cscal)
       Call apply_scalar_bc(Cscal)
       !$acc update device(Cscal)
       Call profiler_stop(PROF_SCALAR)
    End If

    ! Temperature step 3
    If ( boussinesq_flag >= 1 ) Then
       Call profiler_start(PROF_SCALAR)
       Call compute_rhs_temperature(Tscal, U, V, W, Ft3)
       !$acc kernels present(Tscal,Tscal_o,Ft1,Ft2,Ft3)
       Tscal(2:nxg-1,2:nyg-1,2:nzg-1) = Tscal_o(2:nxg-1,2:nyg-1,2:nzg-1) + &
            dt*( rk_coef(3,1)*Ft1 + rk_coef(3,2)*Ft2 + rk_coef(3,3)*Ft3 )
       !$acc end kernels
       !$acc update host(Tscal)
       Call apply_temperature_bc(Tscal)
       !$acc update device(Tscal)
       If ( ibm_input_mode >= 1 ) Call apply_ghost_cell_ibm_scalar(Tscal)
       Call profiler_stop(PROF_SCALAR)
    End If

    ! ── IBM force output ─────────────────────────────────────────────
    If ( ibm_input_mode >= 1 .And. nsampling > 0 .And. Mod(istep, nsampling) == 0 ) Then
       Call profiler_start(PROF_IBM)
       Call compute_ibm_forces(U, V, W, &
            Fx_ibm,  Fy_ibm,  Fz_ibm, &
            Fx_pres, Fy_pres, Fz_pres, &
            Fx_visc, Fy_visc, Fz_visc)
       Call write_force_csv(Fx_ibm,  Fy_ibm,  Fz_ibm, &
                            Fx_pres, Fy_pres, Fz_pres, &
                            Fx_visc, Fy_visc, Fz_visc)
       ! U_pre/V_pre/W_pre are persistent module-level allocatables reused
       ! at all 3 RK stages each qualifying step; freed only at finalise.
       Call profiler_stop(PROF_IBM)
    End If

    ! ── IBM surface field output (per-point pressure / viscous force) ─
    If ( ibm_input_mode >= 1 .And. ibm_surface_nsampling > 0 .And. &
         Mod(istep, ibm_surface_nsampling) == 0 ) Then
       Call profiler_start(PROF_IBM)
       Call sample_ibm_surface(U, V, W)
       Call profiler_stop(PROF_IBM)
    End If

    ! restore the pre-snap dt so next step's CFL-based scaling isn't anchored to the output-snapped value
    dt = dt_presnap

  End Subroutine compute_time_step_RK3

End Module time_integration
