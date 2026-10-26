! Module for temporal integration
Module time_integration

  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use mpi
  Use equations
  Use adv_hybrid, Only : adv_pad
  Use projection
  Use boundary_conditions
  Use wallmodel
  Use sgs_models
  Use ibm
  Use ibm_stress,       Only : ibm_stress_update, nu_const, ibm_stair_forces, ovr_u, ovr_v, ovr_w
  Use vof_pressure,     Only : vp_project_staircase, vp_mu, vp_mv, vp_mw
  Use scalar_transport,  Only : compute_rhs_scalar, apply_scalar_bc
  Use thermal_transport, Only : compute_rhs_temperature, apply_temperature_bc
  Use monitor,          Only : compute_cfl, select_dt, write_force_csv, compute_bulk_velocity
  Use particles,        Only : advance_particles
  Use vof_state,        Only : vof_advance_substep, vof_save_stage, vof_end_step, vof_output_monitor
  Use profiler
  Use debug_trace,      Only : trace_stage, trace_interior, trace_active

  Implicit None

  ! Velocity snapshot buffers for the IBM force computation, allocated on the first qualifying step and kept
  ! until finalise.
  Real(Int64), Allocatable :: U_pre(:,:,:), V_pre(:,:,:), W_pre(:,:,:)

Contains

  !> Projection: the fast Poisson solver; with the staircase IBM (ibm_method = 1) the closed faces are zeroed
  !  around it
  Subroutine project_velocity_field

    If ( ibm_method == 1 ) Then
       Call vp_project_staircase
    Else
       Call compute_projection_step
    End If

  End Subroutine project_velocity_field

  !> Predictor with the last stage pressure: U -= dt*c*grad(Pst), via project_velocity (same seam/outflow face
  !  handling as the projection). The IBM block that follows zeroes the body faces, so the right-hand side of
  !  the whole-box solve has no flux through the body and the solve supplies only the increment (the
  !  pressure-free predictor of ibm_pressure = 0 discards the full-gradient flux).
  Subroutine add_pressure_predictor(c)

    Real(Int64), Intent(In) :: c

    !$acc kernels present(rhs_p,Pst)
    rhs_p = dt*c*Pst
    !$acc end kernels
    If ( outflow_backflow_beta > 0d0 ) Then
       !$acc kernels present(phib,Pb)
       phib = dt*c*Pb
       !$acc end kernels
    End If
    Call apply_periodic_xz_pressure
    Call update_ghost_interior_planes_pressure
    Call update_ghost_interior_planes_pressure_x
    Call project_velocity

  End Subroutine add_pressure_predictor

  !> Stage pressure after the projection: Pst += (solved pseudo-pressure)/(dt*c). After the last stage Pst is
  !  copied into P (snapshots and restart files); read_pressure_restart restores Pst so a hot start continues
  !  the incremental form.
  Subroutine update_stage_pressure(c)

    Real(Int64), Intent(In) :: c

    ! interior only: the ghost planes of rhs_p still hold the predictor's dt*c*Pst, so a whole-array update
    ! would double them every stage
    !$acc kernels present(rhs_p,Pst)
    Pst(2:nxg-1,:,2:nzg-1) = Pst(2:nxg-1,:,2:nzg-1) + rhs_p(2:nxg-1,:,2:nzg-1)/(dt*c)
    !$acc end kernels
    If ( outflow_backflow_beta > 0d0 ) Then
       !$acc kernels present(phib,Pb)
       Pb = Pb + phib/(dt*c)
       !$acc end kernels
    End If
    If ( rk_step == 3 ) Then
       !$acc kernels present(P,Pst)
       P(2:nxg-1,2:nyg-1,2:nzg-1) = Pst(2:nxg-1,2:nyg-1,2:nzg-1)
       !$acc end kernels
       Call finish_pressure_ghosts(P)
    End If

  End Subroutine update_stage_pressure

  !> Wall-model stress on the staircase faces of the immersed body (ibm_method = 1 with the IBM wall model)
  Subroutine update_ibm_stress

    If ( .Not. ibm_stress_on ) Return
    Call ibm_stress_update(U, V, W, vp_mu, vp_mv, vp_mw, nu_const)

  End Subroutine update_ibm_stress

  !> Update imposed pressure gradients for oscillatory forcing:
  !  dPdx=dPdx_t+Ub_x*waveOmega_x*cos(waveOmega_x*t+phi_wave_x), analogous in z
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
    Integer(Int32) :: i, j, k
    Real(Int64) :: Fx_ibm,  Fy_ibm,  Fz_ibm
    Real(Int64) :: Fx_pres, Fy_pres, Fz_pres
    Real(Int64) :: Fx_visc, Fy_visc, Fz_visc
    Real(Int64) :: cfl_conv, cfl_visc, cfl_accel, dt_presnap
    Real(Int64) :: Ub_now, dU_cmfr
    Logical     :: needs_final_sync

    Call trace_stage('step_start')

    ! Enforce IBM before saving the old state (zeroes solid cells on step 1, no-op afterwards); U,V,W are
    ! device-current, so no sync is needed
    If ( ibm_input_mode >= 1 ) Then
       Call profiler_start(PROF_IBM)
       Call apply_ghost_cell_ibm(U, V, W)
       Call profiler_stop(PROF_IBM)
    End If

    ! SGS eddy viscosity before the CFL check so the viscous CFL uses the current nu_t (no-op for DNS);
    ! compute_sgs_model leaves nu_t device-resident
    Call profiler_start(PROF_SGS)
    Call compute_sgs_model(U,V,W,nu_t)
    Call profiler_stop(PROF_SGS)

    ! CFL check and adaptive dt from the start-of-step velocity, so dt is enforced on the current step; a
    ! fixed dt needs it only for the monitor line
    If ( cfl_adaptive == 1 .Or. Mod(istep, nmonitor) == 0 ) Then
       Call profiler_start(PROF_CFL)
       Call compute_cfl(cfl_conv, cfl_visc, cfl_accel)
       Call profiler_stop(PROF_CFL)
       cfl_conv_last  = cfl_conv
       cfl_visc_last  = cfl_visc
       cfl_accel_last = cfl_accel
       cfl_current    = Max(cfl_conv, cfl_visc, cfl_accel)
    End If
    ! nsave<0 / nsteps<0: dt is shrunk so t lands exactly on tsave_next / sim_end_time (uniform physical-time
    ! snapshots under adaptive dt)
    Call select_dt(cfl_current, dt_presnap)

    to = t
    Call profiler_start(PROF_RK_UPDATE)
    ! device-to-device copy as explicit parallel loops: as whole-array assignments under `kernels` nvfortran
    ! ran them in scalar mode (one thread, ~80% of the step at 8.6M points)
    !$acc parallel loop collapse(3) present(U,Uo)
    Do k = 1, Size(U,3)
       Do j = 1, Size(U,2)
          Do i = 1, Size(U,1)
             Uo(i,j,k) = U(i,j,k)
          End Do
       End Do
    End Do
    !$acc end parallel loop
    !$acc parallel loop collapse(3) present(V,Vo)
    Do k = 1, Size(V,3)
       Do j = 1, Size(V,2)
          Do i = 1, Size(V,1)
             Vo(i,j,k) = V(i,j,k)
          End Do
       End Do
    End Do
    !$acc end parallel loop
    !$acc parallel loop collapse(3) present(W,Wo)
    Do k = 1, Size(W,3)
       Do j = 1, Size(W,2)
          Do i = 1, Size(W,1)
             Wo(i,j,k) = W(i,j,k)
          End Do
       End Do
    End Do
    !$acc end parallel loop
    ! particles.f90's added-mass term reads Uo/Vo/Wo host-side; on a GPU build they stay device-only until
    ! synced
    If ( particles_active >= 1 .And. particle_mode == 1 .And. particle_added_mass == 1 ) Then
       !$acc update host(Uo,Vo,Wo)
    End If
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

    ! Imposed pressure gradient (oscillatory / steady); under constant-mass-flux forcing dPdx has no direct
    ! RHS forcing (the CMFR correction below is the only forcing)
    If ( flow_forcing_mode == 0 ) Then
       Call update_pressure_forcing
    Else
       dPdx = 0d0
    End If

    ! RK stage 1
    rk_step = 1
    ! compute_wall_model is device-resident (except the host-only compute_pseudo_pressure_bc_for_robin_bc), so
    ! no update device is needed
    Call profiler_start(PROF_WALLMODEL)
    Call update_wm_filter(U,W)
    Call compute_wall_model(U,V,W) ! sets the wall stress + EQWM ghost cells on old velocity
    Call update_ibm_stress
    Call profiler_stop(PROF_WALLMODEL)
    Call profiler_start(PROF_RHS)
    If ( advection_scheme == 2 ) Call adv_pad(U,V,W)
    Call compute_rhs_u(U,V,W,Fu1)
    Call compute_rhs_v(U,V,W,Fv1)
    Call compute_rhs_w(U,V,W,Fw1)
    Call profiler_stop(PROF_RHS)

    ! Scalar RHS from the same stage-start state (C, T, U,V,W) as the momentum RHS so the scalar RK3 sees
    ! consistent stage times; its update follows the velocity stage below
    If ( sediment_flag >= 1 .Or. boussinesq_flag >= 1 ) Call profiler_start(PROF_SCALAR)
    If ( sediment_flag >= 1 ) Call compute_rhs_scalar(Cscal, U, V, W, Fcs1)
    If ( boussinesq_flag >= 1 ) Call compute_rhs_temperature(Tscal, U, V, W, Ft1)
    If ( trace_active() ) Then
       If ( sediment_flag   >= 1 ) Then
          !$acc update host(Fcs1)
          Call trace_interior('s1_rhsS', 'Fc', Fcs1, 1)
       End If
       If ( boussinesq_flag >= 1 ) Then
          !$acc update host(Ft1)
          Call trace_interior('s1_rhsS', 'Ft', Ft1, 1)
       End If
       Call trace_stage('s1_rhsS')
    End If
    If ( sediment_flag >= 1 .Or. boussinesq_flag >= 1 ) Call profiler_stop(PROF_SCALAR)

    ! RK-stage velocity update, GPU-resident (Fu1/Fw1 never leave the device; Fv1 round-trips
    ! to the host inside apply_uav_forcing when uav_active>=1)
    Call profiler_start(PROF_RK_UPDATE)
    !$acc kernels present(U,V,W,Uo,Vo,Wo,Fu1,Fv1,Fw1)
    U(2:nx-1,2:nyg-1,2:nzg-1) = Uo(2:nx-1,2:nyg-1,2:nzg-1) + dt*rk_coef(1,1)*Fu1
    V(2:nxg-1,2:ny-1,2:nzg-1) = Vo(2:nxg-1,2:ny-1,2:nzg-1) + dt*rk_coef(1,1)*Fv1
    W(2:nxg-1,2:nyg-1,2:nz-1) = Wo(2:nxg-1,2:nyg-1,2:nz-1) + dt*rk_coef(1,1)*Fw1
    !$acc end kernels
    Call profiler_stop(PROF_RK_UPDATE)
    ! ibm_pressure: the stage pressure of the last solve enters the predictor before the body faces are zeroed
    If ( ibm_pressure > 0 ) Call add_pressure_predictor(rk_coef(1,1))

    ! Ghost-cell IBM on the updated velocity (2nd order; no-slip or EQWM): zero solid cells first, then
    ! overwrite ghost cells with EQWM if active
    If ( ibm_input_mode >= 1 ) Then
       Call profiler_start(PROF_IBM)
       ! Stage 1: zero accumulators, capture the first IBM impulse; only the host-only impulse snapshot needs
       ! a sync, when active
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
          Call compute_ibm_wall_model(U, V, W)
       End If
       ! No update device needed: apply_ghost_cell_ibm/compute_ibm_wall_model write on the device and the
       ! host-only reads above never modify U,V,W
       Call profiler_stop(PROF_IBM)
    End If

    t = to + rk_t(rk_step)*dt

    If ( vof_active >= 1 ) Call vof_advance_substep(1)

    Call profiler_start(PROF_BC)
    Call apply_boundary_conditions
    Call trace_stage('s1_bc_pre')
    Call profiler_stop(PROF_BC)
    Call project_velocity_field
    If ( ibm_pressure > 0 ) Call update_stage_pressure(rk_coef(1,1))   ! the solve returned the increment
    Call trace_stage('s1_proj')
    Call profiler_start(PROF_BC)
    Call apply_boundary_conditions(after_projection=.True.)
    Call trace_stage('s1_bc_post1')
    Call profiler_stop(PROF_BC)
    ! U,V,W are correct on the device; the IBM re-enforce below needs a host mirror only for the impulse
    ! snapshot. ibm_pressure = 2 skips the re-enforcement: it would discard the flux the solve put through the
    ! body faces.
    If ( ibm_input_mode >= 1 .And. ibm_pressure < 2 ) Then
       ! Re-enforce IBM: the projection's grad(p) correction corrupts ghost-cell velocities; accumulate that
       ! impulse so Method 1 captures the full force.
       Call profiler_start(PROF_IBM)
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
       If ( ibm_wall_model_flag == 1 ) Call compute_ibm_wall_model(U, V, W)
       Call profiler_stop(PROF_IBM)
    End If

    If ( vof_active >= 1 ) Call vof_save_stage(1)

    ! Scalar RK stage 1
    If ( sediment_flag >= 1 ) Then
       Call profiler_start(PROF_SCALAR)
       !$acc kernels present(Cscal,Cscal_o,Fcs1)
       Cscal(2:nxg-1,2:nyg-1,2:nzg-1) = Cscal_o(2:nxg-1,2:nyg-1,2:nzg-1) + dt*rk_coef(1,1)*Fcs1
       !$acc end kernels
       Call apply_scalar_bc(Cscal)
       If ( ibm_input_mode >= 1 ) Call apply_ghost_cell_ibm_scalar_noflux(Cscal)
       Call trace_stage('s1_C')
       Call profiler_stop(PROF_SCALAR)
    End If

    ! Temperature RK stage 1 (buoyancy in stage n's compute_rhs_v reads Tscal as finalized at the end of stage
    ! n-1)
    If ( boussinesq_flag >= 1 ) Then
       Call profiler_start(PROF_SCALAR)
       !$acc kernels present(Tscal,Tscal_o,Ft1)
       Tscal(2:nxg-1,2:nyg-1,2:nzg-1) = Tscal_o(2:nxg-1,2:nyg-1,2:nzg-1) + dt*rk_coef(1,1)*Ft1
       !$acc end kernels
       Call apply_temperature_bc(Tscal)
       If ( ibm_input_mode >= 1 ) Call apply_ghost_cell_ibm_scalar(Tscal)
       Call trace_stage('s1_T')
       Call profiler_stop(PROF_SCALAR)
    End If

    ! RK stage 2
    rk_step = 2
    ! t is now this stage's time, so time-dependent forcing is re-evaluated per stage
    If ( flow_forcing_mode == 0 ) Call update_pressure_forcing
    Call profiler_start(PROF_SGS)
    ! U,V,W are device-resident; nothing writes them from the host
    Call compute_sgs_model(U,V,W,nu_t)
    Call profiler_stop(PROF_SGS)
    Call profiler_start(PROF_WALLMODEL)
    Call compute_wall_model(U,V,W)
    Call update_ibm_stress
    Call profiler_stop(PROF_WALLMODEL)
    Call profiler_start(PROF_RHS)
    If ( advection_scheme == 2 ) Call adv_pad(U,V,W)
    Call compute_rhs_u(U,V,W,Fu2)
    Call compute_rhs_v(U,V,W,Fv2)
    Call compute_rhs_w(U,V,W,Fw2)
    Call profiler_stop(PROF_RHS)

    ! Scalar RHS from the same stage-start state (C, T, U,V,W) as the momentum RHS so the scalar RK3 sees
    ! consistent stage times; its update follows the velocity stage below
    If ( sediment_flag >= 1 .Or. boussinesq_flag >= 1 ) Call profiler_start(PROF_SCALAR)
    If ( sediment_flag >= 1 ) Call compute_rhs_scalar(Cscal, U, V, W, Fcs2)
    If ( boussinesq_flag >= 1 ) Call compute_rhs_temperature(Tscal, U, V, W, Ft2)
    If ( trace_active() ) Then
       If ( sediment_flag   >= 1 ) Then
          !$acc update host(Fcs2)
          Call trace_interior('s2_rhsS', 'Fc', Fcs2, 1)
       End If
       If ( boussinesq_flag >= 1 ) Then
          !$acc update host(Ft2)
          Call trace_interior('s2_rhsS', 'Ft', Ft2, 1)
       End If
       Call trace_stage('s2_rhsS')
    End If
    If ( sediment_flag >= 1 .Or. boussinesq_flag >= 1 ) Call profiler_stop(PROF_SCALAR)

    Call profiler_start(PROF_RK_UPDATE)
    !$acc kernels present(U,V,W,Uo,Vo,Wo,Fu1,Fv1,Fw1,Fu2,Fv2,Fw2)
    U(2:nx-1,2:nyg-1,2:nzg-1) = Uo(2:nx-1,2:nyg-1,2:nzg-1) + dt*( rk_coef(2,1)*Fu1 + rk_coef(2,2)*Fu2 )
    V(2:nxg-1,2:ny-1,2:nzg-1) = Vo(2:nxg-1,2:ny-1,2:nzg-1) + dt*( rk_coef(2,1)*Fv1 + rk_coef(2,2)*Fv2 )
    W(2:nxg-1,2:nyg-1,2:nz-1) = Wo(2:nxg-1,2:nyg-1,2:nz-1) + dt*( rk_coef(2,1)*Fw1 + rk_coef(2,2)*Fw2 )
    !$acc end kernels
    ! ibm_pressure: the stage pressure of the last solve enters the predictor before the body faces are zeroed
    If ( ibm_pressure > 0 ) Call add_pressure_predictor(rk_coef(2,2))
    Call profiler_stop(PROF_RK_UPDATE)

    If ( ibm_input_mode >= 1 ) Then
       Call profiler_start(PROF_IBM)
       ! Stage 2: capture the IBM impulse; only the host-only impulse snapshot below needs a sync, when active
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
          Call compute_ibm_wall_model(U, V, W)
       End If
       ! No update device needed: apply_ghost_cell_ibm/compute_ibm_wall_model write on the device and the
       ! host-only reads above never modify U,V,W
       Call profiler_stop(PROF_IBM)
    End If

    t = to + rk_t(rk_step)*dt

    If ( vof_active >= 1 ) Call vof_advance_substep(2)

    Call profiler_start(PROF_BC)
    Call apply_boundary_conditions
    Call trace_stage('s2_bc_pre')
    Call profiler_stop(PROF_BC)
    Call project_velocity_field
    If ( ibm_pressure > 0 ) Call update_stage_pressure(rk_coef(2,2))   ! the solve returned the increment
    Call trace_stage('s2_proj')
    Call profiler_start(PROF_BC)
    Call apply_boundary_conditions(after_projection=.True.)
    Call trace_stage('s2_bc_post1')
    Call profiler_stop(PROF_BC)
    ! U,V,W are correct on the device; the IBM re-enforce below needs a host mirror only for the impulse
    ! snapshot. ibm_pressure = 2 skips the re-enforcement: it would discard the flux the solve put through the
    ! body faces.
    If ( ibm_input_mode >= 1 .And. ibm_pressure < 2 ) Then
       ! Re-enforce IBM: the projection's grad(p) correction corrupts ghost-cell velocities; accumulate that
       ! impulse so Method 1 captures the full force.
       Call profiler_start(PROF_IBM)
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
       If ( ibm_wall_model_flag == 1 ) Call compute_ibm_wall_model(U, V, W)
       Call profiler_stop(PROF_IBM)
    End If

    If ( vof_active >= 1 ) Call vof_save_stage(2)

    ! Scalar RK stage 2
    If ( sediment_flag >= 1 ) Then
       Call profiler_start(PROF_SCALAR)
       !$acc kernels present(Cscal,Cscal_o,Fcs1,Fcs2)
       Cscal(2:nxg-1,2:nyg-1,2:nzg-1) = Cscal_o(2:nxg-1,2:nyg-1,2:nzg-1) + &
            dt*( rk_coef(2,1)*Fcs1 + rk_coef(2,2)*Fcs2 )
       !$acc end kernels
       Call apply_scalar_bc(Cscal)
       If ( ibm_input_mode >= 1 ) Call apply_ghost_cell_ibm_scalar_noflux(Cscal)
       Call trace_stage('s2_C')
       Call profiler_stop(PROF_SCALAR)
    End If

    ! Temperature RK stage 2
    If ( boussinesq_flag >= 1 ) Then
       Call profiler_start(PROF_SCALAR)
       !$acc kernels present(Tscal,Tscal_o,Ft1,Ft2)
       Tscal(2:nxg-1,2:nyg-1,2:nzg-1) = Tscal_o(2:nxg-1,2:nyg-1,2:nzg-1) + &
            dt*( rk_coef(2,1)*Ft1 + rk_coef(2,2)*Ft2 )
       !$acc end kernels
       Call apply_temperature_bc(Tscal)
       If ( ibm_input_mode >= 1 ) Call apply_ghost_cell_ibm_scalar(Tscal)
       Call trace_stage('s2_T')
       Call profiler_stop(PROF_SCALAR)
    End If

    ! RK stage 3
    rk_step = 3
    ! t is now this stage's time, so time-dependent forcing is re-evaluated per stage
    If ( flow_forcing_mode == 0 ) Call update_pressure_forcing
    Call profiler_start(PROF_SGS)
    ! U,V,W are device-resident; nothing writes them from the host
    Call compute_sgs_model(U,V,W,nu_t)
    Call profiler_stop(PROF_SGS)
    Call profiler_start(PROF_WALLMODEL)
    Call compute_wall_model(U,V,W)
    Call update_ibm_stress
    Call profiler_stop(PROF_WALLMODEL)
    Call profiler_start(PROF_RHS)
    If ( advection_scheme == 2 ) Call adv_pad(U,V,W)
    Call compute_rhs_u(U,V,W,Fu3)
    Call compute_rhs_v(U,V,W,Fv3)
    Call compute_rhs_w(U,V,W,Fw3)
    Call profiler_stop(PROF_RHS)

    ! Scalar RHS from the same stage-start state (C, T, U,V,W) as the momentum RHS so the scalar RK3 sees
    ! consistent stage times; its update follows the velocity stage below
    If ( sediment_flag >= 1 .Or. boussinesq_flag >= 1 ) Call profiler_start(PROF_SCALAR)
    If ( sediment_flag >= 1 ) Call compute_rhs_scalar(Cscal, U, V, W, Fcs3)
    If ( boussinesq_flag >= 1 ) Call compute_rhs_temperature(Tscal, U, V, W, Ft3)
    If ( trace_active() ) Then
       If ( sediment_flag   >= 1 ) Then
          !$acc update host(Fcs3)
          Call trace_interior('s3_rhsS', 'Fc', Fcs3, 1)
       End If
       If ( boussinesq_flag >= 1 ) Then
          !$acc update host(Ft3)
          Call trace_interior('s3_rhsS', 'Ft', Ft3, 1)
       End If
       Call trace_stage('s3_rhsS')
    End If
    If ( sediment_flag >= 1 .Or. boussinesq_flag >= 1 ) Call profiler_stop(PROF_SCALAR)

    Call profiler_start(PROF_RK_UPDATE)
    !$acc kernels present(U,V,W,Uo,Vo,Wo,Fu1,Fv1,Fw1,Fu2,Fv2,Fw2,Fu3,Fv3,Fw3)
    U(2:nx-1,2:nyg-1,2:nzg-1) = Uo(2:nx-1,2:nyg-1,2:nzg-1) + &
         dt*( rk_coef(3,1)*Fu1 + rk_coef(3,2)*Fu2 + rk_coef(3,3)*Fu3 )
    V(2:nxg-1,2:ny-1,2:nzg-1) = Vo(2:nxg-1,2:ny-1,2:nzg-1) + &
         dt*( rk_coef(3,1)*Fv1 + rk_coef(3,2)*Fv2 + rk_coef(3,3)*Fv3 )
    W(2:nxg-1,2:nyg-1,2:nz-1) = Wo(2:nxg-1,2:nyg-1,2:nz-1) + &
         dt*( rk_coef(3,1)*Fw1 + rk_coef(3,2)*Fw2 + rk_coef(3,3)*Fw3 )
    !$acc end kernels
    ! ibm_pressure: the stage pressure of the last solve enters the predictor before the body faces are zeroed
    If ( ibm_pressure > 0 ) Call add_pressure_predictor(rk_coef(3,3))
    Call profiler_stop(PROF_RK_UPDATE)

    If ( ibm_input_mode >= 1 ) Then
       Call profiler_start(PROF_IBM)
       ! Stage 3: capture the final IBM impulse (accumulators zeroed at stage 1); only the host-only impulse
       ! snapshot needs a sync
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
          Call compute_ibm_wall_model(U, V, W)
       End If
       ! No update device needed: apply_ghost_cell_ibm/compute_ibm_wall_model write on the device and the
       ! host-only reads above never modify U,V,W
       Call profiler_stop(PROF_IBM)
    End If

    t = to + rk_t(rk_step)*dt

    If ( vof_active >= 1 ) Call vof_advance_substep(3)

    Call profiler_start(PROF_BC)
    Call apply_boundary_conditions
    Call trace_stage('s3_bc_pre')
    Call profiler_stop(PROF_BC)
    Call project_velocity_field
    If ( ibm_pressure > 0 ) Call update_stage_pressure(rk_coef(3,3))   ! the solve returned the increment
    Call trace_stage('s3_proj')
    Call profiler_start(PROF_BC)
    Call apply_boundary_conditions(after_projection=.True.)
    Call trace_stage('s3_bc_post1')
    Call profiler_stop(PROF_BC)

    ! Constant-mass-flux forcing: shift U by a uniform constant so the bulk velocity hits Ub_target exactly
    ! (div-free is preserved by a uniform shift)
    If ( flow_forcing_mode == 1 ) Then
       Call profiler_start(PROF_CMFR)
       Call compute_bulk_velocity(Ub_now)
       dU_cmfr = Ub_target - Ub_now
       !$acc kernels present(U)
       U(2:nx-1,2:nyg-1,2:nzg-1) = U(2:nx-1,2:nyg-1,2:nzg-1) + dU_cmfr
       !$acc end kernels
       Call profiler_stop(PROF_CMFR)
       Call profiler_start(PROF_BC)
       Call apply_boundary_conditions(after_projection=.True.)
       Call trace_stage('s3_bc_post2')
       Call profiler_stop(PROF_BC)
    End If

    ! Final sync: host U,V,W are needed only if a host-only consumer runs (RSB, monitor/divergence check,
    ! snapshot, slice/line probe)
    needs_final_sync = ( particles_active >= 1 ) .Or. &
         ( rsb_active == 1 .And. istep >= rsb_nstart ) .Or. &
         ( Mod(istep, nmonitor) == 0 ) .Or. &
         ( nsave > 0 .And. Mod(istep, nsave) == 0 ) .Or. &
         ( nsave < 0 .And. t >= tsave_next - 1d-10 ) .Or. &
         ( chunked_snapshot .And. chunk_nsave > 0 .And. Mod(istep, Max(chunk_nsave,1)) == 0 ) .Or. &
         ( chunked_snapshot .And. chunk_nsave < 0 .And. t >= chunk_tsave_next - 1d-10 ) .Or. &
         ( n_slices > 0 .And. slice_freq > 0 .And. Mod(istep, slice_freq) == 0 ) .Or. &
         ( n_slices > 0 .And. slice_freq < 0 .And. t >= slice_tnext - 1d-10 ) .Or. &
         ( n_lines  > 0 .And. line_freq  > 0 .And. Mod(istep, line_freq)  == 0 ) .Or. &
         ( n_lines  > 0 .And. line_freq  < 0 .And. t >= line_tnext  - 1d-10 )
    If ( needs_final_sync ) Then
       Call profiler_start(PROF_RK_UPDATE)
       If ( ibm_input_mode == 0 ) Then
          !$acc update host(U,V,W)
       End If
       !$acc update host(P)
       ! U,V,W with the IBM are refreshed after the final re-enforce below. P (and nu_t in LES) are
       ! device-resident too: refresh their host copies for the same consumers; snapshots also sync them in
       ! input_output.f90
       If ( sgs_model /= 0 ) Then
          !$acc update host(nu_t)
       End If
       Call profiler_stop(PROF_RK_UPDATE)
    End If
    ! Re-enforce IBM: the projection's grad(p) correction corrupts ghost-cell velocities; accumulate that
    ! impulse so Method 1 captures the full force. ibm_pressure = 2 skips this: it would discard the flux the
    ! solve put through the body faces.
    If ( ibm_input_mode >= 1 .And. ibm_pressure < 2 ) Then
       Call profiler_start(PROF_IBM)
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
       If ( ibm_wall_model_flag == 1 ) Call compute_ibm_wall_model(U, V, W)
       ! Refresh host mirror for the host-only consumers (monitor/output/check_divergence, particles, probes)
       If ( needs_final_sync ) Then
          !$acc update host(U,V,W)
       End If
       Call profiler_stop(PROF_IBM)
    Else If ( ibm_input_mode >= 1 .And. needs_final_sync ) Then
       ! ibm_pressure = 2: no re-enforcement, but the host mirror is still needed by the host-only consumers
       !$acc update host(U,V,W)
    End If

    ! Advance particles once per full step against the final host-resident U,V,W; needs_final_sync guarantees
    ! the host mirror is current whenever particles_active>=1
    If ( particles_active >= 1 ) Call advance_particles

    ! Scalar RK stage 3
    If ( sediment_flag >= 1 ) Then
       Call profiler_start(PROF_SCALAR)
       !$acc kernels present(Cscal,Cscal_o,Fcs1,Fcs2,Fcs3)
       Cscal(2:nxg-1,2:nyg-1,2:nzg-1) = Cscal_o(2:nxg-1,2:nyg-1,2:nzg-1) + &
            dt*( rk_coef(3,1)*Fcs1 + rk_coef(3,2)*Fcs2 + rk_coef(3,3)*Fcs3 )
       !$acc end kernels
       Call apply_scalar_bc(Cscal)
       If ( ibm_input_mode >= 1 ) Call apply_ghost_cell_ibm_scalar_noflux(Cscal)
       Call trace_stage('s3_C')
       Call profiler_stop(PROF_SCALAR)
    End If

    ! Temperature RK stage 3
    If ( boussinesq_flag >= 1 ) Then
       Call profiler_start(PROF_SCALAR)
       !$acc kernels present(Tscal,Tscal_o,Ft1,Ft2,Ft3)
       Tscal(2:nxg-1,2:nyg-1,2:nzg-1) = Tscal_o(2:nxg-1,2:nyg-1,2:nzg-1) + &
            dt*( rk_coef(3,1)*Ft1 + rk_coef(3,2)*Ft2 + rk_coef(3,3)*Ft3 )
       !$acc end kernels
       Call apply_temperature_bc(Tscal)
       If ( ibm_input_mode >= 1 ) Call apply_ghost_cell_ibm_scalar(Tscal)
       Call trace_stage('s3_T')
       Call profiler_stop(PROF_SCALAR)
    End If

    ! ── IBM force output ─────────────────────────────────────────────
    If ( ibm_input_mode >= 1 .And. nsampling > 0 .And. Mod(istep, nsampling) == 0 ) Then
       Call profiler_start(PROF_IBM)
       If ( ibm_method == 1 ) Then
          !$acc update host(U,V,W,P,nu_t)
          Call ibm_stair_forces(U, V, W, Fx_ibm, Fy_ibm, Fz_ibm, Fx_pres, Fy_pres, Fz_pres, Fx_visc, Fy_visc, Fz_visc)
       Else
          Call compute_ibm_forces(U, V, W, &
               Fx_ibm,  Fy_ibm,  Fz_ibm, &
               Fx_pres, Fy_pres, Fz_pres, &
               Fx_visc, Fy_visc, Fz_visc)
       End If
       Call write_force_csv(Fx_ibm,  Fy_ibm,  Fz_ibm, &
                            Fx_pres, Fy_pres, Fz_pres, &
                            Fx_visc, Fy_visc, Fz_visc)
       ! U_pre/V_pre/W_pre are module-level allocatables reused at all 3 RK stages of each qualifying step;
       ! freed at finalise
       Call profiler_stop(PROF_IBM)
    End If

    ! ── IBM surface field output (per-point pressure / viscous force) ─
    If ( ibm_input_mode >= 1 .And. ibm_surface_nsampling > 0 .And. &
         Mod(istep, ibm_surface_nsampling) == 0 ) Then
       Call profiler_start(PROF_IBM)
       Call sample_ibm_surface(U, V, W)
       Call profiler_stop(PROF_IBM)
    End If

    If ( vof_active >= 1 ) Then
       Call vof_end_step
       If ( Mod(istep, nmonitor) == 0 ) Call vof_output_monitor
    End If

    Call trace_stage('step_end')
    dt_step = dt

    ! restore the pre-snap dt so the next step's CFL scaling is not anchored to the output-snapped value
    dt = dt_presnap

  End Subroutine compute_time_step_RK3

End Module time_integration
