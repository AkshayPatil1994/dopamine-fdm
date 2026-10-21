! Module for I/O
Module input_output

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use mpi
  Use genGridAndIC
  Use particles, Only : write_particle_restart, write_particle_snapshot
  Use time_stats, Only : write_tavg

  ! prevent implicit typing
  Implicit None

  Integer(Int32), Parameter :: snap_magic = -4242   ! first word of a box snapshot; a restart file starts with nx_global > 0

Contains

  !  Read input parameters from a Fortran-namelist file. Defaults to
  !  'input_parameters' in the working directory (the main solver's usual
  !  call, see initialize() in initialization.f90); pass input_file to read
  !  a different path instead (e.g. dopamine_esem_main.f90's --input= option).
  Subroutine read_input_parameters(input_file)

    Character(*), Intent(In), Optional :: input_file
    Character(:), Allocatable :: fname
    Integer(Int32) :: ios, unit_in
    Real(Int64)    :: grav_sediment
    ! Local aliases so the namelist uses short, user-friendly names
    Integer(Int32) :: nx, ny, nz, nks
    Real   (Int64) :: Lx, Ly, Lz, alpha_grid
    Real   (Int64) :: T_wave_x, T_wave_z

    ! ---- Namelist group declarations --------------------------------
    Namelist /DOMAIN/ nx, ny, nz, Lx, Ly, Lz, alpha_grid, grid_type, p_row, p_col, alpha_grid_z

    Namelist /PHYSICS/ nu, dPdx, dPdz, Ub_x, Ub_z, T_wave_x, T_wave_z, &
                       phi_wave_x, phi_wave_z, sgs_model, Cs_vreman, flat_wall_model_flag, &
                       z0_ylo, z0_yhi, z0h_ylo, z0h_yhi, wm_filter_mode, wm_filter_Tf_fac, &
                       flow_forcing_mode, Ub_target, advection_scheme, adv_order, adv_beta0, adv_beta_ibm, adv_band, &
                       rotation_active, Omega_x

    Namelist /NUMERICS/ dt, nsteps, nsave, nmonitor, sim_end_time, tsave, &
                        cfl_adaptive, cfl_target, cfl_safety, dt_min, dt_max, pcg_precond, mg_cycle, mg_smooth_deg

    Namelist /BOUNDARY_CONDITIONS/ bc_face_ylo, bc_face_yhi, x_bc_type, y_bc_type, z_bc_type, outflow_backflow_beta

    Namelist /INFLOW/ inflow_type, inflow_Uconst, inflow_profile_file, inflow_temperature_file, &
                       sem_profile_format, sem_Lscale_ratio_y, sem_Lscale_ratio_z, &
                       sem_n_eddies, sem_length_scale, sem_seed, &
                       sem_ensemble_samples, sem_ensemble_periods, &
                       sem_sigma_file, sem_eddy_placement, sem_use_esem, &
                       sem_divergence_free, sem_wall_damping, sem_wall_damping_Aplus, &
                       inflow_recycle_file, inflow_recycle_loop, inflow_recycle_t_offset, &
                       inflow_recycle_shift_z, inflow_recycle_seed

    Namelist /IBM/ ibm_input_mode, ibm_method, ibm_pressure, ibm_wall_model_flag, &
                   ibm_sdf_file, ibm_objid_file, ks, nks, nsampling, ibm_surface_nsampling, &
                   ibm_T_bc_type, ibm_T_wall, ibm_z0, smooth_ibm

    Namelist /INITIAL_CONDITIONS/ Utarget, nstep_init, restart, t_start, &
                                   scalar_restart, ic_type, noise_percent

    Namelist /IO/ filein, fileout, full_snapshot, chunked_snapshot, snapshot_location, &
                  x_start, x_end, y_start, y_end, z_start, z_end

    Namelist /SEDIMENT/ sediment_flag, d_s, rho_s, rho_f, grav, Sc, Sc_t, C_ref, sed_bc_bot, &
                         C_ic_type, C_ic_height

    Namelist /BOUSSINESQ/ boussinesq_flag, beta_T, T_ref, grav, Pr, Pr_t, &
                           T_bc_bot, T_bc_top, T_wall_bot, T_wall_top, &
                           T_ic_type, T_ic_grad

    Namelist /STATISTICS/ rsb_active, rsb_freq, rsb_nstart, rsb_hom_dir, rsb_fileout, &
                          tavg_active, tavg_nstart, tavg_reset, tavg_T, tavg_C, &
                          n_slices, slice_freq, slice_interval, slice_dir, slice_pos, slice_comps, slice_fileout, &
                          n_lines,  line_freq,  line_interval, line_dir,  line_pos1,  line_pos2,   &
                          line_start, line_end, line_comps, line_fileout

    Namelist /UAV/ uav_active, uav_xc, uav_yc, uav_zc, uav_disk_radius, &
                   uav_n_r, uav_n_theta, uav_hover_thrust, uav_kernel_ncell, &
                   uav_path_active, uav_path_file, uav_thrust_active, uav_thrust_file, &
                   uav_load_profile, uav_tilt_active, uav_tilt_tau, uav_grav, uav_swirl_frac

    Namelist /PARTICLES/ particles_active, n_particles_init, &
                         particle_seed_xmin, particle_seed_xmax, particle_seed_ymin, particle_seed_ymax, &
                         particle_seed_zmin, particle_seed_zmax, particle_seed_seed, &
                         bc_particle_x, bc_particle_y, bc_particle_z, &
                         particle_reinit_on_exit, particle_max_age, particle_restart_file, &
                         particle_restart_load, &
                         particle_mode, particle_diam, particle_rho, particle_rho_f, &
                         particle_added_mass, particle_brownian, particle_temp_abs, &
                         particle_ibm_bc, particle_ibm_tau_crit, particle_resuspend_ucrit, &
                         particle_boussinesq_coupling, particle_deposit_file, particle_deposit_freq, &
                         sgs_particle_model, particle_langevin_C0

    Namelist /WAVES/ wave_type, wave_height, wave_period, wave_phase, wave_sf_n, wave_current_mode, wave_Hs, wave_Tp, &
                     wave_gamma, wave_nfreq, wave_seed, wave_gen_len, wave_abs_len, wave_relax_rate, wave_gauge_x, wave_ramp_time, &
                     current_type, current_U, current_z0, current_n, wave_turb

    Namelist /VOF/ vof_active, vof_ic_type, vof_level, vof_center, vof_radius, vof_normal_scheme, vof_method, vof_beta, &
                   vof_geo_density, vof_inertia, &
                   vof_flow, vof_rho_l, vof_rho_g, vof_nu_l, vof_nu_g, vof_grav, vof_sigma, &
                   vof_pcg_iters, vof_pcg_tol, vof_mom_scheme, vof_mom_cm0, vof_mom_cm1, &
                   vof_rk_mom, vof_rk_nth, vof_co_sub, vof_nsub, vof_freeze_ut, vof_adv_iters, vof_adv_tol, &
                   vof_tr_floor, &
                   vof_div_tol, &
                   vof_u0, vof_wave_amp, vof_wave_lambda, vof_wave_stokes, vof_tgv, vof_prescribed, vof_presc_T, &
                   vof_selftest, vof_smooth_w, vof_hsplit, vof_debug, vof_snap_dt

    ! ---- Defaults (variables not in the file keep these values) ------
    nx = 4; ny = 4; nz = 4
    Lx = 1d0; Ly = 1d0; Lz = 1d0
    alpha_grid = 1d0
    grid_type  = 1
    p_row = 0; p_col = 0
    Ub_x = 0d0;  Ub_z = 0d0
    T_wave_x = 0d0;  T_wave_z = 0d0
    phi_wave_x = 0d0;  phi_wave_z = 0d0
    flow_forcing_mode = 0
    Ub_target = 0d0
    nks = 0
    ic_type       = 1
    noise_percent = 5.0d0
    advection_scheme = 0
    adv_order = 3;  adv_beta0 = 0.05d0;  adv_beta_ibm = 1d0;  adv_band = 2
    nu = -1d0;  dt = 0d0;  dPdx = 0d0;  dPdz = 0d0;  Utarget = 0d0   ! nu and dt must be set (checked after the read)
    nsteps = 0; nsave = 0; nmonitor = 0   ! no valid default: must be set in &NUMERICS (checked after the read)

    fname = 'input_parameters'
    If ( Present(input_file) ) fname = input_file

    ! ---- Only rank 0 reads the file ----------------------------------
    If ( myid==0 ) Then

       Open(newunit=unit_in, file=fname, status='old',    &
            action='read', iostat=ios)
       If ( ios /= 0 ) Call abort_input( 'ERROR: cannot open input parameters file: '//fname )

       Rewind(unit_in)
       Read(unit_in, nml=DOMAIN,              iostat=ios)
       If (ios /= 0) Call abort_input( 'ERROR: &DOMAIN missing or failed to parse (check for unrecognized variable names)' )

       Rewind(unit_in)
       Read(unit_in, nml=PHYSICS,             iostat=ios)
       If (ios /= 0) Call abort_input( 'ERROR: &PHYSICS missing or failed to parse (check for unrecognized variable names)' )

       Rewind(unit_in)
       Read(unit_in, nml=NUMERICS,            iostat=ios)
       If (ios /= 0) Call abort_input( 'ERROR: &NUMERICS missing or failed to parse (check for unrecognized variable names)' )

       If ( namelist_group_present(unit_in, 'BOUNDARY_CONDITIONS') ) Then
          Rewind(unit_in)
          Read(unit_in, nml=BOUNDARY_CONDITIONS, iostat=ios)
          If (ios /= 0) Call abort_input( 'ERROR: &BOUNDARY_CONDITIONS present but failed to parse (check variable names)' )
       Else
          Write(*,'(A)') ' INFO: no &BOUNDARY_CONDITIONS found, using defaults'
       End If

       If ( namelist_group_present(unit_in, 'INFLOW') ) Then
          Rewind(unit_in)
          Read(unit_in, nml=INFLOW,              iostat=ios)
          If (ios /= 0) Call abort_input( 'ERROR: &INFLOW present but failed to parse (check variable names)' )
       Else If ( x_bc_type == 1 ) Then
          Write(*,'(A)') ' INFO: no &INFLOW found, using defaults (inflow_type=0, constant uniform flow)'
       End If

       If ( namelist_group_present(unit_in, 'IBM') ) Then
          Rewind(unit_in)
          Read(unit_in, nml=IBM,                 iostat=ios)
          If (ios /= 0) Call abort_input( 'ERROR: &IBM present but failed to parse (check variable names)' )
       Else
          Write(*,'(A)') ' INFO: no &IBM found, using defaults (no IBM body)'
          ibm_input_mode      = 0
          ibm_method          = 0
          ibm_pressure        = 0
          ibm_wall_model_flag = 0
          ibm_sdf_file        = 'SDF_in'
          ks                  = 0d0
          nks                 = 0
          nsampling           = 0
          ibm_surface_nsampling = 0
          smooth_ibm          = 0
       End If

       Rewind(unit_in)
       Read(unit_in, nml=INITIAL_CONDITIONS,  iostat=ios)
       If (ios /= 0) Call abort_input( 'ERROR: &INITIAL_CONDITIONS missing or failed to parse ' // &
            '(check for unrecognized variable names)' )

       Rewind(unit_in)
       Read(unit_in, nml=IO,                  iostat=ios)
       If (ios /= 0) Call abort_input( 'ERROR: &IO missing or failed to parse (check for unrecognized variable names)' )

       If ( namelist_group_present(unit_in, 'SEDIMENT') ) Then
          Rewind(unit_in)
          Read(unit_in, nml=SEDIMENT,            iostat=ios)
          If (ios /= 0) Call abort_input( 'ERROR: &SEDIMENT present but failed to parse (check variable names)' )
       Else
          Write(*,'(A)') ' INFO: no &SEDIMENT found, scalar transport disabled'
       End If
       grav_sediment = grav

       If ( namelist_group_present(unit_in, 'BOUSSINESQ') ) Then
          Rewind(unit_in)
          Read(unit_in, nml=BOUSSINESQ,          iostat=ios)
          If (ios /= 0) Call abort_input( 'ERROR: &BOUSSINESQ present but failed to parse (check variable names)' )
          If ( sediment_flag >= 1 .And. grav /= grav_sediment ) Then
             Write(*,'(A)') ' WARNING: &SEDIMENT and &BOUSSINESQ specify different grav, using &BOUSSINESQ value'
          End If
       Else
          Write(*,'(A)') ' INFO: no &BOUSSINESQ found, thermal stratification disabled'
       End If

       If ( namelist_group_present(unit_in, 'STATISTICS') ) Then
          Rewind(unit_in)
          Read(unit_in, nml=STATISTICS,          iostat=ios)
          If (ios /= 0) Call abort_input( 'ERROR: &STATISTICS present but failed to parse (check variable names)' )
       Else
          Write(*,'(A)') ' INFO: no &STATISTICS found, Reynolds stress budget disabled'
       End If

       If ( namelist_group_present(unit_in, 'UAV') ) Then
          Rewind(unit_in)
          Read(unit_in, nml=UAV,                 iostat=ios)
          If (ios /= 0) Call abort_input( 'ERROR: &UAV present but failed to parse (check variable names)' )
       Else
          Write(*,'(A)') ' INFO: no &UAV found, UAV actuator disk disabled'
       End If
       ! particle_seed_*max defaults to the full domain; resolved here (using the DOMAIN
       ! group's local Lx/Ly/Lz aliases, already read above) rather than in the declaration
       ! block, since the domain size isn't known until DOMAIN is read.
       If ( namelist_group_present(unit_in, 'PARTICLES') ) Then
          Rewind(unit_in)
          Read(unit_in, nml=PARTICLES,           iostat=ios)
          If (ios /= 0) Call abort_input( 'ERROR: &PARTICLES present but failed to parse (check variable names)' )
          If ( particle_seed_xmax < 0d0 ) particle_seed_xmax = Lx
          If ( particle_seed_ymax < 0d0 ) particle_seed_ymax = Ly
          If ( particle_seed_zmax < 0d0 ) particle_seed_zmax = Lz
          If ( bc_particle_x < 0 ) bc_particle_x = Merge(0, 1, x_bc_type == 0)
          If ( bc_particle_y < 0 ) bc_particle_y = Merge(0, 2, y_bc_type == 0)
          If ( bc_particle_z < 0 ) bc_particle_z = Merge(0, 2, z_bc_type == 0)
          If ( particle_reinit_on_exit < 0 .Or. particle_reinit_on_exit > 1 ) Then
             Call abort_input( 'ERROR: &PARTICLES particle_reinit_on_exit must be 0 (none) or 1 (inflow)' )
          End If
       Else
          Write(*,'(A)') ' INFO: no &PARTICLES found, point-particle tracking disabled'
       End If

       If ( namelist_group_present(unit_in, 'WAVES') ) Then
          Rewind(unit_in)
          Read(unit_in, nml=WAVES,               iostat=ios)
          If (ios /= 0) Call abort_input( 'ERROR: &WAVES present but failed to parse (check variable names)' )
       End If

       If ( namelist_group_present(unit_in, 'VOF') ) Then
          Rewind(unit_in)
          Read(unit_in, nml=VOF,                 iostat=ios)
          If (ios /= 0) Call abort_input( 'ERROR: &VOF present but failed to parse (check variable names)' )
       Else
          Write(*,'(A)') ' INFO: no &VOF found, free-surface VOF disabled'
       End If

       Close(unit_in)

       If ( nsteps == 0 ) Call abort_input( 'ERROR: &NUMERICS nsteps must be set (>0 fixed step count, <0 run until sim_end_time)' )
       If ( n_slices > 0 .And. slice_freq < 0 .And. slice_interval <= 0d0 ) &
          Call abort_input( 'ERROR: &STATISTICS slice_freq<0 requires slice_interval > 0' )
       If ( n_lines > 0 .And. line_freq < 0 .And. line_interval <= 0d0 ) &
          Call abort_input( 'ERROR: &STATISTICS line_freq<0 requires line_interval > 0' )
       If ( nsave == 0 ) Call abort_input( 'ERROR: &NUMERICS nsave must be set (>0 every nsave steps, <0 every tsave time units)' )
       If ( nmonitor <= 0 ) Call abort_input( 'ERROR: &NUMERICS nmonitor must be set to a positive step interval' )
       ! Map short namelist names to the global variable names
       nx_global  = nx;  ny_global  = ny;  nz_global  = nz
       Lx_i       = Lx;  Ly_i       = Ly;  Lz_i       = Lz
       Call check_vof_inputs
       Call check_run_inputs
       If ( pcg_precond < -1 .Or. pcg_precond > 1 ) &
            Call abort_input( 'ERROR: pcg_precond must be 0 (fast Poisson solver), 1 (multigrid) or -1 (automatic)' )
       If ( mg_cycle == 0 .Or. mg_cycle < -1 .Or. mg_cycle > 3 ) &
            Call abort_input( 'ERROR: mg_cycle must be 1 (V), 2 (W), 3 or -1 (automatic)' )
       ! the W cycle halves the PCG iterations but visits the coarse levels 2^L times: launch-bound tiny kernels on a GPU
       If ( mg_cycle < 0 ) Then
          mg_cycle = 2
#ifdef GPU_POISSON
          mg_cycle = 1
#endif
       End If
       If ( mg_smooth_deg < 1 ) Call abort_input( 'ERROR: mg_smooth_deg must be >= 1' )
       If ( pcg_precond < 0 ) Then
          pcg_precond = 0
          If ( y_bc_type /= 0 ) pcg_precond = 1
       End If
       If ( ibm_method < 0 .Or. ibm_method > 1 ) &
            Call abort_input( 'ERROR: ibm_method must be 0 (ghost-cell) or 1 (staircase)' )
       If ( ibm_pressure < 0 .Or. ibm_pressure > 2 ) Call abort_input( 'ERROR: ibm_pressure must be 0, 1 or 2' )
       If ( ibm_pressure > 0 .And. ibm_input_mode < 1 ) &
            Call abort_input( 'ERROR: ibm_pressure > 0 needs an immersed body (ibm_input_mode >= 1)' )
       If ( outflow_backflow_beta < 0d0 ) Call abort_input( 'ERROR: outflow_backflow_beta must be >= 0' )
       If ( outflow_backflow_beta > 0d0 .And. ( x_bc_type /= 1 .Or. ibm_pressure < 1 ) ) &
            Call abort_input( 'ERROR: outflow_backflow_beta > 0 needs x_bc_type = 1 and ibm_pressure >= 1' )
       If ( ibm_pressure == 2 .And. ibm_method /= 1 ) &
            Call abort_input( 'ERROR: ibm_pressure = 2 needs the staircase method (ibm_method = 1)' )
       If ( ibm_method == 1 .And. ibm_input_mode < 1 ) &
            Call abort_input( 'ERROR: ibm_method=1 needs an immersed body (ibm_input_mode >= 1)' )
#ifdef GPU_POISSON
       If ( pcg_precond < 1 .And. vof_flow >= 1 ) &
            Call abort_input( 'ERROR: the two-fluid solver in a GPU build needs the multigrid preconditioner (pcg_precond = 1)' )
#endif
       alphaGrid  = alpha_grid
       ! grid_type is used directly from the global module (no alias needed)
       nks_global = nks

       ! Derive angular frequencies from wave periods (omega = 2*pi/T; 0 if T=0)
       If ( T_wave_x > 0d0 ) Then
          waveOmega_x = 8d0*Atan(1d0) / T_wave_x   ! 2*pi / T_wave_x
       Else
          waveOmega_x = 0d0
       End If
       If ( T_wave_z > 0d0 ) Then
          waveOmega_z = 8d0*Atan(1d0) / T_wave_z   ! 2*pi / T_wave_z
       Else
          waveOmega_z = 0d0
       End If

       ! Time-dependent forcing starts from the input gradients
       dPdx_t   = dPdx
       dPdz_t   = dPdz

       ! Rotation axis (streamwise, x) passes through the domain centerline
       y0_rot = 0.5d0 * Ly_i
       z0_rot = 0.5d0 * Lz_i

       If ( flow_forcing_mode == 1 .And. T_wave_x > 0d0 ) Then
          Call abort_input( 'ERROR: flow_forcing_mode=1 (constant mass flux) is incompatible with oscillatory forcing (T_wave_x)' )
       End If
       If ( flow_forcing_mode == 1 .And. x_bc_type /= 0 ) Then
          Call abort_input( 'ERROR: flow_forcing_mode=1 (constant mass flux) requires periodic streamwise BC (x_bc_type=0)' )
       End If
       ! 4-wall duct: the coupled 2D (y,z) pressure solve (solve_poisson_equation) needs the full
       ! z-extent locally available within the y-pencil for every rank, which only holds when z is
       ! not MPI-split -- i.e. p_col==1 (x-only decomposition, via p_row). Auto mode (p_row=p_col=0)
       ! handles this itself (decomp_auto_factorize forces p_col=1 for this BC combination); an
       ! explicit p_col/=1 request is rejected here rather than silently doing the wrong thing.
       If ( alpha_grid_z > 0d0 .And. z_bc_type == 0 ) Then
          Call abort_input( 'ERROR: alpha_grid_z>0 (spanwise grid stretching) requires z_bc_type=1 -- periodic z ' // &
               '(z_bc_type=0) is FFT-based and needs uniform spacing' )
       End If
       If ( z_bc_type == 1 .And. flat_wall_model_flag == 2 ) Then
          Call abort_input( 'ERROR: flat_wall_model_flag=2 (rough EQWM) is not yet supported for z walls (z_bc_type=1) ' // &
               '-- use flat_wall_model_flag=0 (DNS no-slip) or 1 (smooth Reichardt EQWM) with a spanwise wall' )
       End If
       If ( y_bc_type == 1 .And. z_bc_type == 1 .And. p_col > 1 ) Then
          Call abort_input( 'ERROR: y_bc_type=1 and z_bc_type=1 (4-wall duct) requires p_col=1 (decompose only in ' // &
               'x, via p_row) -- the 2D (y,z) pressure solve needs the full z-extent local to every ' // &
               'rank; leave p_row/p_col=0 for auto, or set p_col=1 explicitly' )
       End If

       Write(*,'(A)') ' Input parameters read from namelist file.'
       Write(*,'(A,3I6)')    '   Grid      (nx,ny,nz)        = ', nx, ny, nz
       Write(*,'(A,3F10.4)') '   Domain    (Lx,Ly,Lz)        = ', Lx, Ly, Lz
       Write(*,'(A,I2)')     '   grid_type                   = ', grid_type
       Write(*,'(A,2I4)')    '   p_row, p_col (0,0=auto)     = ', p_row, p_col
       If ( vof_flow >= 1 ) Then
          Write(*,'(A)')     '   nu                          not used (two-fluid viscosities vof_nu_l, vof_nu_g)'
       Else
          Write(*,'(A,E12.4)')  '   nu                          = ', nu
       End If
       Write(*,'(A,I2)')     '   advection_scheme (0=skew-symmetric,1=pure central,2=skew+upwind) = ', advection_scheme
       If ( advection_scheme == 2 ) Write(*,'(A,I2,A,F6.3,A,F6.3,A,I2)') '   adv_order ', adv_order, '  adv_beta0 ', adv_beta0, &
            '  adv_beta_ibm ', adv_beta_ibm, '  adv_band ', adv_band
       Write(*,'(A,I2)')     '   flow_forcing_mode (0=dPdx,1=const mass flux) = ', flow_forcing_mode
       If ( flow_forcing_mode == 1 ) Then
          Write(*,'(A,E12.4)') '   Ub_target                   = ', Ub_target
       End If
       Write(*,'(A,E12.4)')  '   dPdx                        = ', dPdx
       Write(*,'(A,E12.4)')  '   dPdz                        = ', dPdz
       If ( rotation_active >= 1 ) Then
          Write(*,'(A,E12.4)') '   Omega_x (streamwise rotation) = ', Omega_x
          Write(*,'(A,2E12.4)') '   rotation axis (y0_rot,z0_rot) = ', y0_rot, z0_rot
       End If
       If ( T_wave_x > 0d0 ) Then
          Write(*,'(A,E12.4)') '   Ub_x                        = ', Ub_x
          Write(*,'(A,E12.4)') '   T_wave_x                    = ', T_wave_x
          Write(*,'(A,E12.4)') '   wave_omega_x (derived)      = ', waveOmega_x
          Write(*,'(A,E12.4)') '   phi_wave_x                  = ', phi_wave_x
       End If
       If ( T_wave_z > 0d0 ) Then
          Write(*,'(A,E12.4)') '   Ub_z                        = ', Ub_z
          Write(*,'(A,E12.4)') '   T_wave_z                    = ', T_wave_z
          Write(*,'(A,E12.4)') '   wave_omega_z (derived)      = ', waveOmega_z
          Write(*,'(A,E12.4)') '   phi_wave_z                  = ', phi_wave_z
       End If
       Write(*,'(A,E12.4)')  '   dt                          = ', dt
       If ( nsteps < 0 ) Then
          Write(*,'(A,E12.4)') '   sim_end_time (nsteps<0)     = ', sim_end_time
       Else
          Write(*,'(A,I10)')   '   nsteps                      = ', nsteps
       End If
       If ( nsave < 0 ) Then
          Write(*,'(A,E12.4)') '   tsave (nsave<0)             = ', tsave
       Else
          Write(*,'(A,I10)')   '   nsave                       = ', nsave
       End If
       Write(*,'(A,L2,A,L2,A,I2)') '   full_snapshot, chunked, location = ', full_snapshot, ',', chunked_snapshot, ',', &
            snapshot_location
       If ( .Not. full_snapshot ) Write(*,'(A)') '   NOTE: full_snapshot = F, no restart file will be written'
       If ( full_snapshot .And. snapshot_location == 1 ) &
            Write(*,'(A)') '   NOTE: snapshot_location = 1 writes a cell-centred full snapshot that cannot be used for restart'
       Write(*,'(A,I2)')     '   sgs_model                   = ', sgs_model
       Write(*,'(A,E12.4)')  '   Cs_vreman                   = ', Cs_vreman
       Write(*,'(A,I2)')     '   flat_wall_model_flag        = ', flat_wall_model_flag
       If ( flat_wall_model_flag == 2 ) Then
          Write(*,'(A,E12.4,A,E12.4)') '   z0_ylo, z0_yhi              = ', z0_ylo,  ', ', z0_yhi
          Write(*,'(A,E12.4,A,E12.4)') '   z0h_ylo, z0h_yhi            = ', z0h_ylo, ', ', z0h_yhi
       End If
       If ( wm_filter_mode > 0 ) Write(*,'(A,I2,A,F8.3)') '   wm_filter_mode, Tf_fac      = ', &
          wm_filter_mode, ', ', wm_filter_Tf_fac
       Write(*,'(A,2I3)')    '   BC y-walls (ylo/yhi)        = ', bc_face_ylo, bc_face_yhi
       Write(*,'(A,I2)')     '   x_bc_type (0=periodic,1=inflow/outflow) = ', x_bc_type
       If ( outflow_backflow_beta > 0d0 ) Write(*,'(A,E12.4)') '   outflow_backflow_beta       = ', outflow_backflow_beta
       Write(*,'(A,I2)')     '   y_bc_type (0=periodic,1=wall)           = ', y_bc_type
       Write(*,'(A,I2)')     '   z_bc_type (0=periodic,1=wall DNS no-slip) = ', z_bc_type
       If ( z_bc_type == 1 .And. alpha_grid_z > 0d0 ) Then
          Write(*,'(A,E12.4)') '   alpha_grid_z (spanwise stretching)      = ', alpha_grid_z
       End If
       If ( x_bc_type == 1 ) Then
          If ( inflow_type == 0 ) Then
             Write(*,'(A,E12.4)') '   inflow_type=0 (uniform), inflow_Uconst = ', inflow_Uconst
          Else If ( inflow_type == 2 ) Then
             Write(*,'(A,A)')  '   inflow_type=2 (recycled precursor slice), inflow_recycle_file = ', Trim(inflow_recycle_file)
             Write(*,'(A,I2)') '   inflow_recycle_loop         = ', inflow_recycle_loop
             Write(*,'(A,E12.4)') '   inflow_recycle_t_offset  = ', inflow_recycle_t_offset
             If ( inflow_recycle_loop == 1 ) Then
                Write(*,'(A,I2)') '   inflow_recycle_shift_z      = ', inflow_recycle_shift_z
                If ( inflow_recycle_shift_z == 1 ) Write(*,'(A,I10)') '   inflow_recycle_seed          = ', inflow_recycle_seed
             End If
          Else If ( inflow_type == 3 ) Then
             Write(*,'(A,I2)') '   inflow_type=3 (wave inlet), wave_type = ', wave_type
          Else If ( inflow_type == 1 ) Then
             Write(*,'(A,A)')     '   inflow_type=1 (SEM), inflow_profile_file = ', Trim(inflow_profile_file)
             If ( sem_profile_format == 0 ) Then
                Write(*,'(A)')    '   sem_profile_format=0 (y U V W uu vv ww uv Reynolds-stress file)'
             Else
                Write(*,'(A)')    '   sem_profile_format=1 (wind-tunnel z U Iu Iv Iw [+length scales] file)'
                Write(*,'(A,F6.3)') '   sem_Lscale_ratio_y (Loy/Lox fallback)  = ', sem_Lscale_ratio_y
                Write(*,'(A,F6.3)') '   sem_Lscale_ratio_z (Loz/Lox fallback)  = ', sem_Lscale_ratio_z
             End If
             If ( sem_use_esem == 1 ) Then
                Write(*,'(A)')    '   sem_use_esem=1 -> Ensemble SEM (empirical/exact normalisation)'
             Else
                Write(*,'(A)')    '   sem_use_esem=0 -> classical SEM (Jarrin 2006, analytical normalisation)'
             End If
             Write(*,'(A,I6)')    '   sem_n_eddies                = ', sem_n_eddies
             Write(*,'(A,E12.4)') '   sem_length_scale            = ', sem_length_scale
             Write(*,'(A,I10)')   '   sem_seed                    = ', sem_seed
             If ( sem_use_esem == 1 ) Then
                Write(*,'(A,I6)') '   sem_ensemble_samples        = ', sem_ensemble_samples
                Write(*,'(A,I6)') '   sem_ensemble_periods        = ', sem_ensemble_periods
             End If
             If ( Len_trim(sem_sigma_file) > 0 ) Then
                Write(*,'(A,A)')  '   sem_sigma_file (inhomogeneous) = ', Trim(sem_sigma_file)
                Write(*,'(A,I2)') '   sem_eddy_placement (0=uniform,1=PDF) = ', sem_eddy_placement
             Else
                Write(*,'(A)')    '   sem_sigma_file not set -> homogeneous length scale (backward compatible)'
             End If
             If ( sem_wall_damping == 1 ) Then
                Write(*,'(A,F6.2)') '   sem_wall_damping=1, sem_wall_damping_Aplus = ', sem_wall_damping_Aplus
             End If
          End If
       End If
       Write(*,'(A,I2)')     '   ibm_input_mode              = ', ibm_input_mode
       If ( ibm_input_mode >= 1 ) Then
          Write(*,'(A,I2)')     '   ibm_method (0 ghost-cell, 1 staircase) = ', ibm_method
          Write(*,'(A,I2)')     '   ibm_pressure (0 off, 1 incr., 2 +flux) = ', ibm_pressure
          Write(*,'(A,I2)')     '   ibm_wall_model_flag         = ', ibm_wall_model_flag
          Write(*,'(A,I4)')     '   smooth_ibm (SDF corner-rounding passes, 0=off) = ', smooth_ibm
       End If
       Write(*,'(A,I8)')     '   nsampling                   = ', nsampling
       Write(*,'(A,I8)')     '   ibm_surface_nsampling       = ', ibm_surface_nsampling
       Write(*,'(A,I2)')     '   cfl_adaptive                = ', cfl_adaptive
       If ( cfl_adaptive == 1 ) Then
          Write(*,'(A,F8.4)') '   cfl_target                  = ', cfl_target
          Write(*,'(A,F8.4)') '   cfl_safety                  = ', cfl_safety
       End If
       Write(*,'(A,I2)')     '   restart                     = ', restart
       If ( restart == 1 ) Then
          Write(*,'(A,I8)')    '   nstep_init                  = ', nstep_init
          If ( t_start >= 0d0 ) &
             Write(*,'(A,F12.6)') '   t_start                     = ', t_start
       End If
       Write(*,'(A,I2)')     '   scalar_restart               = ', scalar_restart
       Write(*,'(A,I2)')     '   ic_type                     = ', ic_type
       Write(*,'(A,F7.2)')   '   noise_percent               = ', noise_percent
       If ( rsb_active == 1 ) Then
          Write(*,'(A,I2)')     '   rsb_active                  = ', rsb_active
          Write(*,'(A,I8)')     '   rsb_freq                    = ', rsb_freq
          Write(*,'(A,I8)')     '   rsb_nstart                  = ', rsb_nstart
          Write(*,'(A,A)')      '   rsb_hom_dir                 = ', Trim(rsb_hom_dir)
          Write(*,'(A,A)')      '   rsb_fileout                 = ', Trim(rsb_fileout)
       End If
       If ( tavg_active == 1 ) Then
          Write(*,'(A,I8)')     '   tavg_nstart                 = ', tavg_nstart
          Write(*,'(A,I2)')     '   tavg_reset                  = ', tavg_reset
          Write(*,'(A,2I2)')    '   tavg_T, tavg_C              = ', tavg_T, tavg_C
       End If
       If ( boussinesq_flag >= 1 ) Then
          Write(*,'(A,I2)')    '   boussinesq_flag              = ', boussinesq_flag
          Write(*,'(A,E12.4)') '   beta_T                       = ', beta_T
          Write(*,'(A,E12.4)') '   T_ref                        = ', T_ref
          Write(*,'(A,E12.4)') '   grav                         = ', grav
          Write(*,'(A,F6.3)')  '   Pr                           = ', Pr
          Write(*,'(A,F6.3)')  '   Pr_t                         = ', Pr_t
          Write(*,'(A,2I3)')   '   T_bc (bot/top, 0=adiabatic,1=isothermal) = ', T_bc_bot, T_bc_top
       End If
       If ( sediment_flag >= 1 ) Then
          Write(*,'(A,I2)')    '   sediment_flag                = ', sediment_flag
          Write(*,'(A,2E12.4)') '   d_s, rho_s                   = ', d_s, rho_s
          Write(*,'(A,E12.4)') '   grav                         = ', grav
          Write(*,'(A,2F6.3)') '   Sc, Sc_t                     = ', Sc, Sc_t
       End If
       If ( uav_active >= 1 ) Write(*,'(A,I2)') '   uav_active                   = ', uav_active
       If ( wave_type > 0 ) Then
          Write(*,'(A,I2)')    '   wave_type (1 Airy,2 stream function,3 JONSWAP) = ', wave_type
          Write(*,'(A,2E12.4)') '   wave_gen_len, wave_abs_len   = ', wave_gen_len, wave_abs_len
          If ( current_type > 0 ) Write(*,'(A,I2)') '   current_type                 = ', current_type
          If ( wave_turb == 1 ) Write(*,'(A)') '   wave_turb=1 (precursor fluctuation on the wave inlet)'
       End If
       If ( n_slices > 0 ) Then
          Write(*,'(A,I4)')  '   n_slices                    = ', n_slices
          Write(*,'(A,I8)')  '   slice_freq                  = ', slice_freq
          If ( slice_freq < 0 ) Write(*,'(A,E12.4)') '   slice_interval              = ', slice_interval
       End If
       If ( n_lines > 0 ) Then
          Write(*,'(A,I4)')  '   n_lines                     = ', n_lines
          Write(*,'(A,I8)')  '   line_freq                   = ', line_freq
          If ( line_freq < 0 ) Write(*,'(A,E12.4)') '   line_interval               = ', line_interval
       End If
       If ( particles_active >= 1 ) Then
          Write(*,'(A,I2)')    '   particles_active            = ', particles_active
          Write(*,'(A,I8)')    '   n_particles_init            = ', n_particles_init
          Write(*,'(A,3I3)')   '   bc_particle (x,y,z) (0=periodic,1=exit,2=reflect,3=absorb) = ', &
               bc_particle_x, bc_particle_y, bc_particle_z
          Write(*,'(A,I2)')    '   particle_reinit_on_exit (0=none,1=inflow) = ', particle_reinit_on_exit
          If ( restart == 1 ) Write(*,'(A,I2)') '   particle_restart_load (0=fresh seed,1=from file) = ', particle_restart_load
          Write(*,'(A,I2)')    '   particle_mode (0=tracer,1=inertial)       = ', particle_mode
          If ( particle_mode >= 1 ) Then
             Write(*,'(A,E12.4)') '   particle_diam                = ', particle_diam
             Write(*,'(A,E12.4)') '   particle_rho                 = ', particle_rho
             Write(*,'(A,E12.4)') '   particle_rho_f               = ', particle_rho_f
             Write(*,'(A,I2)')    '   particle_added_mass          = ', particle_added_mass
             Write(*,'(A,I2)')    '   particle_brownian            = ', particle_brownian
          End If
       End If
       If ( vof_active >= 1 ) Then
          Write(*,'(A,I2)')    '   vof_active                  = ', vof_active
          Write(*,'(A,I2)')    '   vof_ic_type (1 plane,2 drop,3 bubble,4 wave,5 column,6 disk,7 Zalesak) = ', vof_ic_type
          Write(*,'(A,I2)')    '   vof_normal_scheme (1=Youngs,2=height function) = ', vof_normal_scheme
          Write(*,'(A,I2,A,F6.2)') '   vof_method (1=PLIC,2=THINC) = ', vof_method, '  vof_beta = ', vof_beta
          Write(*,'(A,I2)')    '   vof_flow (0=passive,1=two-fluid)       = ', vof_flow
          If ( vof_flow >= 1 ) Then
             Write(*,'(A,2E12.4)') '   vof_rho_l, vof_rho_g        = ', vof_rho_l, vof_rho_g
             Write(*,'(A,2E12.4)') '   vof_nu_l,  vof_nu_g         = ', vof_nu_l, vof_nu_g
             Write(*,'(A,E12.4)')  '   vof_grav                    = ', vof_grav
             Write(*,'(A,I3)')     '   vof_pcg_iters               = ', vof_pcg_iters
             Write(*,'(A,I2)')     '   vof_mom_scheme (0 central,1 upwind,2 Koren,4 QUICK,6 WENO5-Z) = ', vof_mom_scheme
          End If
       End If
    End If

    ! ---- Broadcast everything to all ranks ---------------------------

    Call Mpi_bcast ( nx_global,            1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( ny_global,            1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( nz_global,            1, MPI_integer,   0, MPI_COMM_WORLD, ierr )

    Call Mpi_bcast ( Lx_i,                 1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( Ly_i,                 1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( Lz_i,                 1, MPI_real8,     0, MPI_COMM_WORLD, ierr )

    Call Mpi_bcast ( alphaGrid,            1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( alpha_grid_z,         1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( grid_type,            1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( Utarget,              1, MPI_real8,     0, MPI_COMM_WORLD, ierr )

    Call Mpi_bcast ( p_row,                1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( p_col,                1, MPI_integer,   0, MPI_COMM_WORLD, ierr )

    Call Mpi_bcast ( ibm_input_mode,       1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( ibm_method,           1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( ibm_pressure,         1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( pcg_precond,          1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( mg_cycle,             1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( mg_smooth_deg,        1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( ibm_sdf_file,  Len(ibm_sdf_file),  MPI_character, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( ibm_objid_file, Len(ibm_objid_file), MPI_character, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( ibm_wall_model_flag,  1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( smooth_ibm,           1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( nks_global,           1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( ks,                   1, MPI_real8,     0, MPI_COMM_WORLD, ierr )

    Call Mpi_bcast ( dt,                   1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( nu,                   1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( advection_scheme,     1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( adv_order,            1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( adv_band,             1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( adv_beta0,            1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( adv_beta_ibm,         1, MPI_real8,     0, MPI_COMM_WORLD, ierr )

    Call Mpi_bcast ( dPdx,                 1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( dPdx_t,               1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( dPdz,                 1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( dPdz_t,               1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( Ub_x,                 1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( Ub_z,                 1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( waveOmega_x,          1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( waveOmega_z,          1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( phi_wave_x,           1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( phi_wave_z,           1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( flow_forcing_mode,    1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( Ub_target,            1, MPI_real8,     0, MPI_COMM_WORLD, ierr )

    Call Mpi_bcast ( rotation_active,      1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( Omega_x,              1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( y0_rot,               1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( z0_rot,               1, MPI_real8,     0, MPI_COMM_WORLD, ierr )

    Call Mpi_bcast ( nstep_init,           1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( t_start,              1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( nsteps,               1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( nsave,                1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( full_snapshot,        1, MPI_logical,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( chunked_snapshot,     1, MPI_logical,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( snapshot_location,    1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( x_start,              1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( x_end,                1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( y_start,              1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( y_end,                1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( z_start,              1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( z_end,                1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( nmonitor,             1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( sim_end_time,         1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( tsave,                1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( nsampling,            1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( ibm_surface_nsampling,1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( cfl_adaptive,         1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( cfl_target,           1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( cfl_safety,           1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( dt_min,               1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( dt_max,               1, MPI_real8,     0, MPI_COMM_WORLD, ierr )

    Call Mpi_bcast ( restart,              1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( scalar_restart,       1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( ic_type,              1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( noise_percent,        1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( bc_face_ylo,          1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( bc_face_yhi,          1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( x_bc_type,            1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( outflow_backflow_beta, 1, MPI_real8,    0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( y_bc_type,            1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( z_bc_type,            1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( inflow_type,          1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( inflow_Uconst,        1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( inflow_profile_file,  Len(inflow_profile_file), MPI_character, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( inflow_temperature_file, Len(inflow_temperature_file), MPI_character, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( sem_profile_format,   1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( sem_Lscale_ratio_y,   1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( sem_Lscale_ratio_z,   1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( sem_n_eddies,         1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( sem_length_scale,     1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( sem_seed,             1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( sem_ensemble_samples, 1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( sem_ensemble_periods, 1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( sem_sigma_file,       Len(sem_sigma_file), MPI_character, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( sem_eddy_placement,   1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( sem_use_esem,         1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( sem_divergence_free,  1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( sem_wall_damping,       1, MPI_integer, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( sem_wall_damping_Aplus, 1, MPI_real8,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( inflow_recycle_file,    Len(inflow_recycle_file), MPI_character, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( inflow_recycle_loop,    1, MPI_integer, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( inflow_recycle_t_offset,1, MPI_real8,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( inflow_recycle_shift_z, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( inflow_recycle_seed,    1, MPI_integer, 0, MPI_COMM_WORLD, ierr )

    Call Mpi_bcast ( sgs_model,            1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( Cs_vreman,            1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( flat_wall_model_flag, 1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( z0_ylo,               1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( z0_yhi,               1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( z0h_ylo,              1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( z0h_yhi,              1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( wm_filter_mode,       1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( wm_filter_Tf_fac,     1, MPI_real8,     0, MPI_COMM_WORLD, ierr )

    Call Mpi_bcast ( sediment_flag,        1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( sed_bc_bot,           1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( d_s,                  1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( rho_s,                1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( rho_f,                1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( grav,                 1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( Sc,                   1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( Sc_t,                 1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( C_ref,                1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( C_ic_type,            1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( C_ic_height,          1, MPI_real8,     0, MPI_COMM_WORLD, ierr )

    Call Mpi_bcast ( uav_active,           1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( uav_xc,               1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( uav_yc,               1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( uav_zc,               1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( uav_disk_radius,      1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( uav_n_r,              1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( uav_n_theta,          1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( uav_hover_thrust,     1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( uav_kernel_ncell,     1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( uav_path_active,      1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( uav_path_file, Len(uav_path_file), MPI_character, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( uav_thrust_active,    1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( uav_thrust_file, Len(uav_thrust_file), MPI_character, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( uav_load_profile,     1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( uav_tilt_active,      1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( uav_tilt_tau,         1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( uav_grav,             1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( uav_swirl_frac,       1, MPI_real8,     0, MPI_COMM_WORLD, ierr )

    Call Mpi_bcast ( boussinesq_flag,      1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( beta_T,               1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( T_ref,                1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( Pr,                   1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( Pr_t,                 1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( T_bc_bot,             1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( T_bc_top,             1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( T_wall_bot,           1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( T_wall_top,           1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( T_ic_type,            1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( T_ic_grad,            1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( ibm_T_bc_type, Size(ibm_T_bc_type), MPI_integer, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( ibm_T_wall,    Size(ibm_T_wall),    MPI_real8,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( ibm_z0,        Size(ibm_z0),        MPI_real8,   0, MPI_COMM_WORLD, ierr )

    Call Mpi_bcast ( rsb_active,           1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( rsb_freq,             1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( rsb_nstart,           1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( rsb_hom_dir, Len(rsb_hom_dir), MPI_character, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( rsb_fileout, Len(rsb_fileout), MPI_character, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( tavg_active,          1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( tavg_nstart,          1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( tavg_reset,           1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( tavg_T,               1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( tavg_C,               1, MPI_integer,   0, MPI_COMM_WORLD, ierr )

    Call Mpi_bcast ( n_slices,             1,              MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( slice_freq,           1,              MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( slice_interval,       1,              MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( slice_dir,    Len(slice_dir(1))*MAX_PROBES,   MPI_character, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( slice_pos,            MAX_PROBES,     MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( slice_comps,  Len(slice_comps(1))*MAX_PROBES, MPI_character, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( slice_fileout,Len(slice_fileout(1))*MAX_PROBES,MPI_character,0, MPI_COMM_WORLD, ierr )

    Call Mpi_bcast ( n_lines,              1,              MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( line_freq,            1,              MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( line_interval,        1,              MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( line_dir,     Len(line_dir(1))*MAX_PROBES,    MPI_character, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( line_pos1,            MAX_PROBES,     MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( line_pos2,            MAX_PROBES,     MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( line_start,           MAX_PROBES,     MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( line_end,             MAX_PROBES,     MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( line_comps,   Len(line_comps(1))*MAX_PROBES,  MPI_character, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( line_fileout, Len(line_fileout(1))*MAX_PROBES,MPI_character, 0, MPI_COMM_WORLD, ierr )

    Call Mpi_bcast ( particles_active,      1, MPI_integer, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( n_particles_init,      1, MPI_integer, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( particle_seed_xmin,    1, MPI_real8,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( particle_seed_xmax,    1, MPI_real8,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( particle_seed_ymin,    1, MPI_real8,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( particle_seed_ymax,    1, MPI_real8,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( particle_seed_zmin,    1, MPI_real8,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( particle_seed_zmax,    1, MPI_real8,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( particle_seed_seed,    1, MPI_integer, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( bc_particle_x,         1, MPI_integer, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( bc_particle_y,         1, MPI_integer, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( bc_particle_z,         1, MPI_integer, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( particle_reinit_on_exit, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( particle_max_age,      1, MPI_real8,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( particle_restart_file, Len(particle_restart_file), MPI_character, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( particle_restart_load, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( particle_mode,         1, MPI_integer, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( particle_diam,         1, MPI_real8,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( particle_rho,          1, MPI_real8,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( particle_rho_f,        1, MPI_real8,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( particle_added_mass,   1, MPI_integer, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( particle_brownian,     1, MPI_integer, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( particle_temp_abs,     1, MPI_real8,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( particle_ibm_bc, Size(particle_ibm_bc), MPI_integer, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( particle_ibm_tau_crit,    1, MPI_real8, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( particle_resuspend_ucrit, 1, MPI_real8, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( particle_boussinesq_coupling, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( particle_deposit_file, Len(particle_deposit_file), MPI_character, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( particle_deposit_freq, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( sgs_particle_model,    1, MPI_integer, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( particle_langevin_C0,  1, MPI_real8,   0, MPI_COMM_WORLD, ierr )

    Call Mpi_bcast ( wave_type,         1, MPI_integer, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( wave_height,       1, MPI_real8, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( wave_period,       1, MPI_real8, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( wave_phase,        1, MPI_real8, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( wave_sf_n,         1, MPI_integer, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( wave_current_mode, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( wave_Hs,           1, MPI_real8, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( wave_Tp,           1, MPI_real8, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( wave_gamma,        1, MPI_real8, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( wave_nfreq,        1, MPI_integer, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( wave_seed,         1, MPI_integer, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( wave_gen_len,      1, MPI_real8, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( wave_abs_len,      1, MPI_real8, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( wave_relax_rate,   1, MPI_real8, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( wave_gauge_x,      MAX_WAVEGAUGES, MPI_real8, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( wave_ramp_time,    1, MPI_real8, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( current_type,      1, MPI_integer, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( wave_turb,         1, MPI_integer, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( current_U,         1, MPI_real8, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( current_z0,        1, MPI_real8, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( current_n,         1, MPI_real8, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( vof_active,            1, MPI_integer, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( vof_ic_type,           1, MPI_integer, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( vof_level,             1, MPI_real8,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( vof_center,            3, MPI_real8,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( vof_radius,            1, MPI_real8,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( vof_normal_scheme,     1, MPI_integer, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( vof_method,            1, MPI_integer, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( vof_geo_density,       1, MPI_integer, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( vof_inertia,           1, MPI_integer, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( vof_beta,              1, MPI_real8,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( vof_flow,              1, MPI_integer, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( vof_rho_l,             1, MPI_real8,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( vof_rho_g,             1, MPI_real8,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( vof_nu_l,              1, MPI_real8,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( vof_nu_g,              1, MPI_real8,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( vof_grav,              1, MPI_real8,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( vof_pcg_iters,         1, MPI_integer, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( vof_pcg_tol,           1, MPI_real8,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( vof_mom_scheme,        1, MPI_integer, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( vof_selftest,          1, MPI_integer, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( vof_u0,                1, MPI_real8,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( vof_adv_iters,         1, MPI_integer, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( vof_adv_tol,           1, MPI_real8,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( vof_tr_floor,          1, MPI_real8,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( vof_div_tol,           1, MPI_real8,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( vof_wave_amp,          1, MPI_real8,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( vof_wave_lambda,       1, MPI_real8,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( vof_nsub,             1, MPI_integer,0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( vof_rk_mom,           1, MPI_integer,0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( vof_rk_nth,           1, MPI_integer,0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( vof_sigma,            1, MPI_real8,  0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( vof_wave_stokes,      1, MPI_integer,0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( vof_snap_dt,          1, MPI_real8,  0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( vof_prescribed,       1, MPI_integer,0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( vof_presc_T,          1, MPI_real8,  0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( vof_debug,            1, MPI_integer,0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( vof_tgv,              1, MPI_integer,0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( vof_freeze_ut,        1, MPI_integer,0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( vof_hsplit,           1, MPI_integer,0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( vof_co_sub,           1, MPI_real8,  0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( vof_mom_cm0,          1, MPI_real8,  0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( vof_mom_cm1,          1, MPI_real8,  0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( vof_smooth_w,         1, MPI_real8,  0, MPI_COMM_WORLD, ierr )

  End Subroutine read_input_parameters

  !> Range checks of the single-phase inputs that had none (an unknown value used to act as a silent default, e.g. sgs_model = 7 ran as DNS)
  !  and warnings for combinations that are accepted but weak or contradictory (see docs/Input-Parameters.md, docs/Numerics.md 6.4)
  Subroutine check_run_inputs

    Call check_numerics
    Call check_snapshot_inputs
    Call check_physics_and_bc
    Call check_combination_warnings

  End Subroutine check_run_inputs


  !> snapshot_location range and the sign convention of the chunk box (all six negative = face indices, none negative = lengths)
  Subroutine check_snapshot_inputs

    Integer(Int32) :: nneg

    If ( snapshot_location < 0 .Or. snapshot_location > 1 ) &
         Call abort_input( 'ERROR: &IO snapshot_location must be 0 (staggered) or 1 (cell-centred)' )
    If ( .Not. chunked_snapshot ) Return
    nneg = Count( [ x_start, x_end, y_start, y_end, z_start, z_end ] < 0d0 )
    If ( nneg /= 0 .And. nneg /= 6 ) &
         Call abort_input( 'ERROR: &IO chunk box x/y/z_start/end must be all lengths (>= 0) or all negative face indices' )

  End Subroutine check_snapshot_inputs


  !> Required values, ranges of the time stepping and the momentum convection
  Subroutine check_numerics

    If ( nu < 0d0 .And. vof_flow < 1 ) Call abort_input( 'ERROR: &PHYSICS nu (kinematic viscosity) must be set' )
    If ( nu < 0d0 ) nu = 0d0   ! two-fluid runs take the viscosities from &VOF
    If ( sgs_model < 0 .Or. sgs_model > 1 ) Call abort_input( 'ERROR: sgs_model must be 0 (DNS) or 1 (Vreman)' )
    If ( advection_scheme < 0 .Or. advection_scheme > 2 ) Call abort_input( 'ERROR: advection_scheme must be 0 ' // &
         '(skew-symmetric), 1 (divergence-form central) or 2 (skew-symmetric + upwind-biased dissipation)' )
    If ( advection_scheme == 2 ) Then
       If ( adv_order /= 3 .And. adv_order /= 5 ) Call abort_input( 'ERROR: adv_order must be 3 or 5' )
       If ( adv_beta0 < 0d0 .Or. adv_beta_ibm < 0d0 ) Call abort_input( 'ERROR: adv_beta0 and adv_beta_ibm must be >= 0' )
       If ( adv_band < 0 .Or. adv_band > 8 ) Call abort_input( 'ERROR: adv_band must be in 0..8' )
       If ( vof_flow >= 1 ) Call abort_input( 'ERROR: advection_scheme = 2 is for the single-phase solver (vof_flow = 0)' )
    End If
    If ( cfl_adaptive < 0 .Or. cfl_adaptive > 1 ) Call abort_input( 'ERROR: cfl_adaptive must be 0 (fixed dt) or 1 (adaptive dt)' )
    If ( dt <= 0d0 ) Call abort_input( 'ERROR: &NUMERICS dt must be > 0 (it is the initial dt when cfl_adaptive = 1)' )
    If ( cfl_adaptive == 1 ) Then
       If ( cfl_target <= 0d0 .Or. cfl_safety <= 0d0 ) Call abort_input( 'ERROR: cfl_target and cfl_safety must be > 0' )
       If ( dt_min > dt_max ) Call abort_input( 'ERROR: dt_min must not exceed dt_max' )
    End If

  End Subroutine check_numerics


  !> Ranges of the wall-model, immersed-body, boundary-condition and feature flags
  Subroutine check_physics_and_bc

    If ( flat_wall_model_flag < 0 .Or. flat_wall_model_flag > 2 ) &
         Call abort_input( 'ERROR: flat_wall_model_flag must be 0 (no-slip), 1 (smooth log law) or 2 (rough log law)' )
    If ( wm_filter_mode < 0 .Or. wm_filter_mode > 1 ) &
         Call abort_input( 'ERROR: wm_filter_mode must be 0 (raw wall-model input) or 1 (temporal filter)' )
    If ( wm_filter_mode == 1 .And. ( flat_wall_model_flag == 0 .Or. wm_filter_Tf_fac <= 0d0 .Or. z_bc_type == 1 ) ) &
         Call abort_input( 'ERROR: wm_filter_mode=1 needs flat_wall_model_flag>0, wm_filter_Tf_fac>0 and no z walls (z_bc_type=0)' )
    If ( ibm_wall_model_flag < 0 .Or. ibm_wall_model_flag > 1 ) &
         Call abort_input( 'ERROR: ibm_wall_model_flag must be 0 (no-slip) or 1 (log-law EQWM)' )
    If ( ibm_input_mode < 0 .Or. ibm_input_mode > 1 ) &
         Call abort_input( 'ERROR: ibm_input_mode must be 0 (no body) or 1 (SDF immersed boundary)' )
    If ( x_bc_type < 0 .Or. x_bc_type > 1 ) Call abort_input( 'ERROR: x_bc_type must be 0 (periodic) or 1 (inflow/outflow)' )
    If ( y_bc_type < 0 .Or. y_bc_type > 1 ) Call abort_input( 'ERROR: y_bc_type must be 0 (periodic) or 1 (walls)' )
    If ( z_bc_type < 0 .Or. z_bc_type > 1 ) Call abort_input( 'ERROR: z_bc_type must be 0 (periodic) or 1 (walls)' )
    If ( sediment_flag < 0 ) Call abort_input( 'ERROR: sediment_flag must be 0 or 1' )
    If ( boussinesq_flag < 0 ) Call abort_input( 'ERROR: boussinesq_flag must be 0 or 1' )

  End Subroutine check_physics_and_bc


  !> Accepted but weak or contradictory combinations: warnings only
  Subroutine check_combination_warnings

    Integer(Int32) :: n

    ! immersed body: weak or contradictory combinations
    If ( ibm_input_mode >= 1 .And. vof_flow < 1 ) Then
       If ( ibm_method == 0 ) Then
          Write(*,'(A)') ' WARNING: ghost-cell IBM (ibm_method = 0) needs walls resolved by >= 2 cells; thin or sharp ' // &
               'geometry (buildings, terrain) is better served by ibm_method = 1, ibm_pressure = 2'
          If ( ibm_wall_model_flag == 1 ) Write(*,'(A)') ' WARNING: ibm_wall_model_flag = 1 with ibm_method = 0 only ' // &
               'rewrites the ghost velocities and carries almost no stress'
       Else If ( ibm_pressure == 0 ) Then
          Write(*,'(A)') ' WARNING: ibm_method = 1 with ibm_pressure = 0 discards the pressure-gradient flux through ' // &
               'the body faces (divergence next to the body, velocity can grow at sharp geometry): ibm_pressure = 2 is recommended'
       End If
       If ( flat_wall_model_flag > 0 ) Write(*,'(A)') &
            ' INFO: flat_wall_model_flag > 0 together with an immersed body: the y-wall model and the IBM wall stress both act'
    End If
    If ( ibm_input_mode < 1 .And. ibm_wall_model_flag == 1 ) &
         Write(*,'(A)') ' WARNING: ibm_wall_model_flag = 1 is ignored without an immersed body (ibm_input_mode = 0)'
    If ( boussinesq_flag < 1 ) Then
       Do n = 0, Size(ibm_T_bc_type) - 1
          If ( ibm_T_bc_type(n) == 1 ) Then
             Write(*,'(A,I0,A)') ' WARNING: ibm_T_bc_type(', n, ') = 1 (isothermal) is ignored without boussinesq_flag = 1'
             Exit
          End If
       End Do
    End If
    If ( x_bc_type == 0 .And. inflow_type /= 0 ) &
         Write(*,'(A)') ' WARNING: &INFLOW inflow_type is ignored with a periodic streamwise direction (x_bc_type = 0)'
    If ( restart == 1 .And. cfl_adaptive == 1 .And. t_start < 0d0 .And. vof_flow < 1 ) &
         Write(*,'(A)') ' WARNING: restart with adaptive dt and no t_start: the time restarts from nstep_init*dt, not the true time'

  End Subroutine check_combination_warnings


  !> Reject feature combinations the free-surface solver does not support (rank 0, right after the namelists are read). UAV is
  !  excluded in both VOF modes; the two-fluid solver (vof_flow>=1) additionally excludes the features listed below.
  Subroutine check_vof_inputs

    If ( vof_active >= 1 .And. uav_active >= 1 ) &
         Call abort_input( 'ERROR: &VOF (vof_active>=1) cannot be combined with the UAV actuator (uav_active>=1)' )
    If ( vof_flow >= 1 .And. vof_active < 1 ) &
         Call abort_input( 'ERROR: vof_flow>=1 requires vof_active=1' )
    If ( wave_type <= 0 .And. ( wave_height > 0d0 .Or. wave_Hs > 0d0 ) ) &
         Write(*,'(A)') ' WARNING: &WAVES wave_height/wave_Hs is set but wave_type = 0: no waves are generated'
    If ( particles_active < 1 .And. n_particles_init > 0 ) &
         Write(*,'(A)') ' WARNING: &PARTICLES n_particles_init > 0 but particles_active = 0: the particles are disabled'
    If ( inflow_type < 0 .Or. inflow_type > 3 ) Call abort_input( 'ERROR: &INFLOW inflow_type must be 0..3' )
    If ( wave_turb /= 0 .And. wave_turb /= 1 ) Call abort_input( 'ERROR: &WAVES wave_turb must be 0 or 1' )
    If ( wave_turb == 1 .And. ( wave_type <= 0 .Or. Len_trim(inflow_recycle_file) == 0 ) ) &
         Call abort_input( 'ERROR: &WAVES wave_turb=1 needs wave_type>0 and inflow_recycle_file (the precursor slice)' )
    If ( wave_type > 0 ) Then
       If ( vof_active < 1 .Or. vof_flow < 1 ) &
            Call abort_input( 'ERROR: &WAVES wave_type>0 requires vof_active=1 and vof_flow=1' )
       If ( x_bc_type /= 1 .Or. inflow_type /= 3 ) &
            Call abort_input( 'ERROR: &WAVES wave_type>0 requires x_bc_type=1 and inflow_type=3 (wave inlet)' )
    End If
    If ( inflow_type == 3 .And. wave_type <= 0 ) &
         Call abort_input( 'ERROR: inflow_type=3 (wave inlet) requires wave_type>0 in &WAVES' )
    If ( ( wave_gen_len > 0d0 .Or. wave_abs_len > 0d0 ) .And. ( wave_type <= 0 .Or. x_bc_type /= 1 ) ) &
         Call abort_input( 'ERROR: wave_gen_len/wave_abs_len (relaxation zones) require wave_type>0 and x_bc_type=1' )
    If ( current_type < 0 .Or. current_type > 3 ) Call abort_input( 'ERROR: &WAVES current_type must be 0..3' )
    If ( current_type > 0 .And. wave_type <= 0 ) &
         Call abort_input( 'ERROR: &WAVES current_type>0 requires wave_type>0 (wave inlet)' )
    If ( current_type == 2 .And. ( current_z0 <= 0d0 .Or. current_z0 >= 0.5d0*vof_level ) ) &
         Call abort_input( 'ERROR: &WAVES current_type=2 needs 0 < current_z0 < vof_level/2' )
    If ( current_type == 3 .And. current_n < 1d0 ) Call abort_input( 'ERROR: &WAVES current_type=3 needs current_n >= 1' )
    If ( wave_abs_len > 0d0 .And. wave_current_mode == 2 ) &
         Write(*,'(A)') ' WARNING: the absorption zone relaxes the velocity to rest: wrong for a mean current (wave_current_mode=2)'
    If ( wave_type > 0 ) Then
       ! gauges: the writer takes the leading entries >= 0, so a hole would write the wrong stations; all must lie inside the flume
       Block
          Integer(Int32) :: ig, ng_lead
          ng_lead = 0
          Do ig = 1, Size(wave_gauge_x)
             If ( wave_gauge_x(ig) < 0d0 ) Exit
             ng_lead = ig
          End Do
          If ( Count(wave_gauge_x >= 0d0) /= ng_lead ) &
               Call abort_input( 'ERROR: &WAVES wave_gauge_x entries must be contiguous from the first (negative = unused)' )
          Do ig = 1, ng_lead
             If ( wave_gauge_x(ig) >= Lx_i ) Call abort_input( 'ERROR: &WAVES wave_gauge_x lies outside the domain (x >= Lx)' )
          End Do
       End Block
       If ( wave_type <= 2 ) Then
          If ( wave_period <= 0d0 .Or. wave_height <= 0d0 ) &
               Call abort_input( 'ERROR: &WAVES wave_type 1/2 need wave_period > 0 and wave_height > 0' )
       Else
          If ( wave_Tp <= 0d0 .Or. wave_Hs <= 0d0 ) &
               Call abort_input( 'ERROR: &WAVES wave_type 3 needs wave_Tp > 0 and wave_Hs > 0' )
       End If
       If ( wave_gen_len > 0d0 .And. wave_abs_len > 0d0 .And. wave_gen_len + wave_abs_len > Lx_i ) &
            Call abort_input( 'ERROR: &WAVES the generation and absorption zones overlap (wave_gen_len + wave_abs_len > Lx)' )
       If ( vof_level <= 0d0 .Or. vof_level >= Ly_i ) Call abort_input( 'ERROR: &VOF vof_level must lie inside (0, Ly) for waves' )
       If ( wave_type <= 2 .And. 0.5d0*wave_height >= vof_level ) &
            Write(*,'(A)') ' WARNING: wave amplitude is not smaller than the still-water depth vof_level'
       If ( wave_type <= 2 .And. vof_level + 0.5d0*wave_height >= Ly_i ) &
            Call abort_input( 'ERROR: &WAVES the wave crest reaches the top of the domain (vof_level + H/2 >= Ly)' )
    End If
    If ( vof_flow >= 1 .And. restart == 1 .And. t_start < 0d0 .And. cfl_adaptive == 1 ) Then
       ! without t_start the restart clock is nstep_init*dt, which is wrong under adaptive dt
       If ( wave_type > 0 .Or. nsave < 0 .Or. nsteps < 0 ) &
            Call abort_input( 'ERROR: restart with adaptive dt needs t_start (the physical time of the snapshot): the clock ' // &
            'nstep_init*dt would reset the wave phase, ramp and output cadence' )
       Write(*,'(A)') ' WARNING: restart with adaptive dt and no t_start: t = nstep_init*dt is not the physical time'
    End If
    If ( vof_flow >= 1 ) Then
       If ( rotation_active >= 1 ) Call abort_input( 'ERROR: vof_flow=1 does not support rotation' )
       If ( particles_active >= 1 .And. Abs(particle_rho_f - vof_rho_l) > 1d-12*vof_rho_l ) &
            Call abort_input( 'ERROR: vof_flow=1 particles move in the liquid: particle_rho_f must equal vof_rho_l' )
       If ( particles_active >= 1 .And. particle_mode >= 1 .And. Abs(grav - vof_grav) > 1d-12*Max(1d0, Abs(vof_grav)) ) &
            Call abort_input( 'ERROR: vof_flow=1 inertial particles need grav (&SEDIMENT/&BOUSSINESQ) equal to vof_grav (&VOF)' )
       If ( boussinesq_flag >= 1 .And. Abs(grav - vof_grav) > 1d-12*Max(1d0, Abs(vof_grav)) ) &
            Call abort_input( 'ERROR: vof_flow=1 with Boussinesq needs grav (&BOUSSINESQ) equal to vof_grav (&VOF)' )
       If ( boussinesq_flag >= 1 .And. vof_hsplit /= 0 ) &
            Call abort_input( 'ERROR: vof_flow=1 Boussinesq buoyancy needs vof_hsplit = 0' )
       If ( dPdx /= 0d0 .Or. dPdz /= 0d0 .Or. flow_forcing_mode /= 0 ) &
            Call abort_input( 'ERROR: vof_flow=1 does not apply dPdx, dPdz or mass-flux forcing (set them to zero)' )
       If ( y_bc_type /= 1 ) Call abort_input( 'ERROR: vof_flow=1 supports wall-bounded y only (y_bc_type=1)' )
       If ( z_bc_type /= 0 ) Call abort_input( 'ERROR: vof_flow=1 supports periodic z only (z_bc_type=0)' )
       If ( x_bc_type == 1 .And. inflow_type /= 3 .And. inflow_type /= 0 ) &
            Call abort_input( 'ERROR: vof_flow=1 with an x inlet supports only inflow_type=3 (wave inlet) or 0 (uniform)' )
       If ( rsb_active == 1 ) &
            Call abort_input( 'ERROR: vof_flow=1 does not support the Reynolds-stress budget' )
       If ( tavg_active == 1 ) &
            Call abort_input( 'ERROR: vof_flow=1 does not support the time-averaged statistics (tavg_active)' )
       If ( ibm_method == 1 ) &
            Call abort_input( 'ERROR: vof_flow=1 always uses the staircase wall stress; leave ibm_method = 0' )
       If ( ibm_input_mode >= 1 .And. ibm_surface_nsampling > 0 ) &
            Call abort_input( 'ERROR: the IBM surface field dump (ibm_surface_nsampling>0) is not available with vof_flow=1' )
    End If

  End Subroutine check_vof_inputs

  !> Report an input-file error and abort every rank (a bare Stop on rank 0 alone leaves the others blocked in the bcast below)
  Subroutine abort_input(msg)

    Character(*), Intent(In) :: msg

    Write(*,'(A)') Trim(msg)
    Call MPI_Abort(MPI_COMM_WORLD, 1, ierr)

  End Subroutine abort_input

  ! Scan the input file for an uncommented "&<group_name>" header, without consuming it; iostat-ambiguity rationale
  Function namelist_group_present(unit_in, group_name) result(found)

    Integer(Int32),   Intent(In) :: unit_in
    Character(Len=*), Intent(In) :: group_name
    Logical :: found

    Character(Len=512) :: line
    Character(Len=Len(group_name)) :: tag_upper
    Integer(Int32) :: ios_local, glen, ic
    Character(Len=1) :: ch

    found = .False.
    glen  = Len_trim(group_name)
    Do ic = 1, glen
       ch = group_name(ic:ic)
       If ( ch >= 'a' .And. ch <= 'z' ) ch = Achar(Iachar(ch) - 32)
       tag_upper(ic:ic) = ch
    End Do

    Rewind(unit_in)
    Do
       Read(unit_in, '(A)', iostat=ios_local) line
       If ( ios_local /= 0 ) Exit   ! EOF: group not found
       line = Adjustl(line)
       If ( line(1:1) == '!' ) Cycle           ! comment line
       If ( line(1:1) /= '&' ) Cycle           ! not a group header
       ! Upper-case just the slice we're about to compare
       Do ic = 1, glen
          ch = line(ic+1:ic+1)
          If ( ch >= 'a' .And. ch <= 'z' ) ch = Achar(Iachar(ch) - 32)
          line(ic+1:ic+1) = ch
       End Do
       If ( line(2:glen+1) /= tag_upper ) Cycle
       ! Require the tag to end here (whitespace/EOL), not be a longer name
       If ( Len_trim(line) > glen+1 ) Then
          ch = line(glen+2:glen+2)
          If ( ch /= ' ' ) Cycle
       End If
       found = .True.
       Exit
    End Do

  End Function namelist_group_present

  !    Generates an initial condition for channel
  ! Output: U,V,W
  Subroutine init_flow

	! Generate the grid
	Call generateGrid

	! Generate the initial conditions
	Call generateIc

	! NOTE: readSDF (IBM setup) is deferred to initialize() so that
	! xg/yg/zg are populated before compute_normal_at_face_* is called.

  End Subroutine init_flow

  ! Read binary snapshot: mesh, U,V and W; Input: filein; Output: U,V,W,x,y,z
  Subroutine read_input_data

    Integer(Int32) ::  nx_global_f,  ny_global_f,  nz_global_f
    Integer(Int32) :: nxm_global_f, nym_global_f, nzm_global_f
    Integer(Int64) :: file_size_actual, file_size_expected

    ! processor 0 Reads the all the data
    If ( myid==0 ) Then

       Write(*,'(A,A,A)') '  Reading restart file: ', Trim(Adjustl(filein)), ' ...'

       Open(1,file=filein,access='stream',form='unformatted',action='Read',convert='little_endian')
       Read(1) nx_global_f
       Close(1)
       If ( nx_global_f == snap_magic ) Then
          Write(*,'(A,A)') ' ERROR: filein is a cell-centred/chunked box snapshot and cannot be used for restart: ', &
               Trim(Adjustl(filein))
          Call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
       End If

       ! ── file-size sanity check: verify expected byte layout ──────────────
       file_size_expected = &
            Int(6,  Int64)*4_Int64 + &
            Int(nx_global + ny_global + nz_global + &
                nxm_global + nym_global + nzm_global, Int64)*8_Int64 + &
            Int(3,  Int64)*4_Int64 + &
            Int(nx_global,  Int64)*Int(nyg_global, Int64)*Int(nzg_global, Int64)*8_Int64 + &
            Int(3,  Int64)*4_Int64 + &
            Int(nxg_global, Int64)*Int(ny_global,  Int64)*Int(nzg_global, Int64)*8_Int64 + &
            Int(3,  Int64)*4_Int64 + &
            Int(nxg_global, Int64)*Int(nyg_global, Int64)*Int(nz_global,  Int64)*8_Int64 + &
            Int(3,  Int64)*4_Int64 + &
            Int(nxg_global, Int64)*Int(nyg_global, Int64)*Int(nzg_global, Int64)*8_Int64
       Inquire(file=Trim(Adjustl(filein)), size=file_size_actual)
       If (file_size_actual < 0_Int64) Then
          Write(*,'(A)') '  WARNING: cannot determine restart file size (skipping check)'
       Else If (file_size_actual < file_size_expected) Then
          Write(*,'(A)') ' '
          Write(*,'(A)')    ' ERROR: restart file is smaller than expected for current grid dimensions.'
          Write(*,'(A,A)')  '   File      : ', Trim(Adjustl(filein))
          Write(*,'(A,I0)') '   Expected  >= ', file_size_expected, ' bytes'
          Write(*,'(A,I0)') '   Found     : ', file_size_actual,   ' bytes'
          Write(*,'(A)')    '   Check that nx/ny/nz in input_parameters match the saved file.'
          Call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
       Else
          Write(*,'(A,I0,A)') '  File size check passed (', file_size_actual, ' bytes)'
       End If

       Open(1,file=filein,access='stream',form='unformatted',action='Read')

       ! mesh
       Read(1) nx_global_f
       If ( nx_global_f/=nx_global ) Then
          Write(*,'(A,I0,A,I0)') ' ERROR: restart file nx = ', nx_global_f, ' but input_parameters gives ', nx_global
          Call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
       End If
       Read(1) x_global

       Read(1) ny_global_f
       If ( ny_global_f/=ny_global ) Then
          Write(*,'(A,I0,A,I0)') ' ERROR: restart file ny = ', ny_global_f, ' but input_parameters gives ', ny_global
          Call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
       End If
       Read(1) y_global

       Read(1) nz_global_f
       If ( nz_global_f/=nz_global ) Then
          Write(*,'(A,I0,A,I0)') ' ERROR: restart file nz = ', nz_global_f, ' but input_parameters gives ', nz_global
          Call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
       End If
       Read(1) z_global

       Read(1) nxm_global_f
       If ( nxm_global_f/=nxm_global ) Then
          Write(*,'(A,I0,A,I0)') ' ERROR: restart file nxm = ', nxm_global_f, ' but input_parameters gives ', nxm_global
          Call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
       End If
       Read(1) xm_global

       Read(1) nym_global_f
       If ( nym_global_f/=nym_global ) Then
          Write(*,'(A,I0,A,I0)') ' ERROR: restart file nym = ', nym_global_f, ' but input_parameters gives ', nym_global
          Call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
       End If
       Read(1) ym_global

       Read(1) nzm_global_f
       If ( nzm_global_f/=nzm_global ) Then
          Write(*,'(A,I0,A,I0)') ' ERROR: restart file nzm = ', nzm_global_f, ' but input_parameters gives ', nzm_global
          Call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
       End If
       Read(1) zm_global

    End If

    ! U, V, W: x/z-decomposed field blocks, gathered/scattered via rank 0 (see read_distributed_field_block)
    Call read_distributed_field_block(1, U, .True.,  .False.)
    Call read_distributed_field_block(1, V, .False., .False.)
    Call read_distributed_field_block(1, W, .False., .True.)

    ! close file
    If (myid==0) Then
       Close(1)
    End If

    ! send data to all other processors
    ! mesh
    Call Mpi_bcast ( x_global,nx_global,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( y_global,ny_global,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( z_global,nz_global,MPI_real8,0,MPI_COMM_WORLD,ierr )

    Call Mpi_bcast ( xm_global,nxm_global,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( ym_global,nym_global,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( zm_global,nzm_global,MPI_real8,0,MPI_COMM_WORLD,ierr )

    ! set solution for zero step
    Uo = U
    Vo = V
    Wo = W

! NOTE: readSDF (IBM setup) is deferred to initialize() so that
	! xg/yg/zg are populated before compute_normal_at_face_* is called.

  End Subroutine read_input_data

  !> Per-rank block shape/placement within a gathered/scattered global array; ix1:ix2,iz1:iz2 as in read/write_distributed_field_block.
  Subroutine distributed_block_range(iproc, is_x_face, is_z_face, ix1, ix2, iz1, iz2)

    Integer(Int32), Intent(In)  :: iproc
    Logical,        Intent(In)  :: is_x_face, is_z_face
    Integer(Int32), Intent(Out) :: ix1, ix2, iz1, iz2

    If ( is_x_face ) Then
       ix1 = i1_global(iproc); ix2 = i2_global(iproc)
    Else
       ix1 = ig1_global(iproc); ix2 = ig2_global(iproc)
    End If
    If ( is_z_face ) Then
       iz1 = k1_global(iproc); iz2 = k2_global(iproc)
    Else
       iz1 = kg1_global(iproc); iz2 = kg2_global(iproc)
    End If

  End Subroutine distributed_block_range

  !> Read one x/z-decomposed field block from a stream unit open on rank 0 only, scattering to every rank's local array. Every non-zero rank's receive is posted up front (non-blocking) so rank 0's sends aren't serialized behind a one-at-a-time round trip.
  Subroutine read_distributed_field_block(unit_no, field_local, is_x_face, is_z_face)

    Integer(Int32), Intent(In) :: unit_no
    Real(Int64), Dimension(:,:,:), Intent(InOut) :: field_local
    Logical, Intent(In) :: is_x_face, is_z_face

    Integer(Int32) :: nn(3), iproc, ix1, ix2, iz1, iz2, nxl, nzl
    Integer(Int32) :: expect_n1, expect_n3
    Real(Int64), Allocatable :: global_buf(:,:,:)
    Type :: send_block_t
       Real(Int64), Allocatable :: d(:,:,:)
    End Type send_block_t
    Type(send_block_t), Allocatable :: send_blocks(:)
    Integer(Int32), Allocatable :: reqs(:)

    If ( myid == 0 ) Then
       Read(unit_no) nn

       ! nn is the field's global (nx1,nx2,nx3) block shape as written; cross-check it against what this call expects (x/z face vs
       ! ghost count per is_x_face/is_z_face, y the caller's full local extent) rather than risk a shape-mismatch crash or silently
       ! wrong data in the ix1:ix2,:,iz1:iz2 slice below.
       expect_n1 = Merge(nx_global, nxg_global, is_x_face)
       expect_n3 = Merge(nz_global, nzg_global, is_z_face)
       If ( nn(1) /= expect_n1 .Or. nn(2) /= Size(field_local,2) .Or. nn(3) /= expect_n3 ) Then
          Write(*,'(A,3(I0,1X),A,3(I0,1X))') ' ERROR: restart file field block size mismatch: found ', &
               nn, ' expected ', expect_n1, Size(field_local,2), expect_n3
          Call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
       End If

       Allocate( global_buf(nn(1), nn(2), nn(3)) )
       Read(unit_no) global_buf

       Call distributed_block_range(0, is_x_face, is_z_face, ix1, ix2, iz1, iz2)
       field_local = global_buf(ix1:ix2, :, iz1:iz2)

       If ( nprocs > 1 ) Then
          Allocate( send_blocks(nprocs-1), reqs(nprocs-1) )
          Do iproc = 1, nprocs-1
             Call distributed_block_range(iproc, is_x_face, is_z_face, ix1, ix2, iz1, iz2)
             nxl = ix2 - ix1 + 1
             nzl = iz2 - iz1 + 1
             Allocate( send_blocks(iproc)%d(nxl, nn(2), nzl) )
             send_blocks(iproc)%d = global_buf(ix1:ix2, :, iz1:iz2)
             Call Mpi_isend(send_blocks(iproc)%d, nxl*nn(2)*nzl, Mpi_real8, iproc, iproc, MPI_COMM_WORLD, reqs(iproc), ierr)
          End Do
          Call Mpi_waitall(nprocs-1, reqs, MPI_STATUSES_IGNORE, ierr)
          Deallocate(send_blocks, reqs)
       End If
       Deallocate(global_buf)
    Else
       Call Mpi_recv(field_local, Size(field_local), Mpi_real8, 0, myid, MPI_COMM_WORLD, istat, ierr)
    End If

  End Subroutine read_distributed_field_block

  !> Skip one x/z-decomposed field block on a stream unit open on rank 0 only, without reading or distributing the payload; used
  !  by read_scalar_restart/read_temperature_restart to pass the U/V/W/P blocks already consumed by read_input_data. Arguments
  !  as in read_distributed_field_block (dim2_size = its Size(field_local,2)).
  Subroutine skip_distributed_field_block(unit_no, dim2_size, is_x_face, is_z_face)

    Integer(Int32), Intent(In) :: unit_no, dim2_size
    Logical, Intent(In) :: is_x_face, is_z_face

    Integer(Int32) :: nn(3), expect_n1, expect_n3
    Integer(Int64) :: cur_pos, payload_bytes

    If ( myid == 0 ) Then
       Read(unit_no) nn

       expect_n1 = Merge(nx_global, nxg_global, is_x_face)
       expect_n3 = Merge(nz_global, nzg_global, is_z_face)
       If ( nn(1) /= expect_n1 .Or. nn(2) /= dim2_size .Or. nn(3) /= expect_n3 ) Then
          Write(*,'(A,3(I0,1X),A,3(I0,1X))') ' ERROR: restart file field block size mismatch while skipping: found ', &
               nn, ' expected ', expect_n1, dim2_size, expect_n3
          Call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
       End If

       payload_bytes = Int(nn(1),Int64)*Int(nn(2),Int64)*Int(nn(3),Int64)*8_Int64
       Inquire(unit=unit_no, POS=cur_pos)
       Read(unit_no, POS=cur_pos+payload_bytes)   ! pure seek: no data item, so nothing is transferred
    End If

  End Subroutine skip_distributed_field_block

  !> Write one x/z-decomposed field block (3-int header + global-sized data) to a stream unit open on rank 0 only, gathering from every rank's local array. Every non-zero rank's contribution is received into its own buffer via a single posted-up-front Waitall instead of one-at-a-time blocking Mpi_recv, so the gather isn't serialized behind rank 0's receive order. is_x_face/is_z_face as in read_distributed_field_block.
  Subroutine write_distributed_field_block(unit_no, field_local, nxf_global, nyf_global, nzf_global, is_x_face, is_z_face)

    Integer(Int32), Intent(In) :: unit_no, nxf_global, nyf_global, nzf_global
    Real(Int64), Dimension(:,:,:), Intent(In) :: field_local
    Logical, Intent(In) :: is_x_face, is_z_face

    Integer(Int32) :: iproc, ix1, ix2, iz1, iz2, nxl, nzl
    Real(Int64), Allocatable :: global_buf(:,:,:)
    Type :: recv_block_t
       Real(Int64), Allocatable :: d(:,:,:)
    End Type recv_block_t
    Type(recv_block_t), Allocatable :: recv_blocks(:)
    Integer(Int32), Allocatable :: reqs(:)

    If ( myid == 0 ) Then
       Allocate( global_buf(nxf_global, nyf_global, nzf_global) )

       Call distributed_block_range(0, is_x_face, is_z_face, ix1, ix2, iz1, iz2)
       global_buf(ix1:ix2, :, iz1:iz2) = field_local

       If ( nprocs > 1 ) Then
          Allocate( recv_blocks(nprocs-1), reqs(nprocs-1) )
          Do iproc = 1, nprocs-1
             Call distributed_block_range(iproc, is_x_face, is_z_face, ix1, ix2, iz1, iz2)
             nxl = ix2 - ix1 + 1
             nzl = iz2 - iz1 + 1
             Allocate( recv_blocks(iproc)%d(nxl, nyf_global, nzl) )
             Call Mpi_irecv(recv_blocks(iproc)%d, nxl*nyf_global*nzl, Mpi_real8, iproc, iproc, MPI_COMM_WORLD, reqs(iproc), ierr)
          End Do
          Call Mpi_waitall(nprocs-1, reqs, MPI_STATUSES_IGNORE, ierr)
          Do iproc = 1, nprocs-1
             Call distributed_block_range(iproc, is_x_face, is_z_face, ix1, ix2, iz1, iz2)
             global_buf(ix1:ix2, :, iz1:iz2) = recv_blocks(iproc)%d
          End Do
          Deallocate(recv_blocks, reqs)
       End If

       Write(unit_no) nxf_global, nyf_global, nzf_global
       Write(unit_no) global_buf
       Deallocate(global_buf)
    Else
       Call Mpi_send(field_local, Size(field_local), Mpi_real8, 0, myid, MPI_COMM_WORLD, ierr)
    End If

  End Subroutine write_distributed_field_block

  !> Coordinate-mode face range f1:f2 of the box s:e on the faces xf(1:nf): nearest faces, clipped to the domain, at least one cell
  !  wide; ok = .False. if s > e or the box lies outside the domain
  Subroutine coordinate_face_range(xf, nf, s, e, f1, f2, ok)

    Integer(Int32), Intent(In)  :: nf
    Real(Int64),    Intent(In)  :: xf(nf), s, e
    Integer(Int32), Intent(Out) :: f1, f2
    Logical,        Intent(Out) :: ok

    ok = ( s <= e .And. e >= xf(1) .And. s <= xf(nf) )
    f1 = Minloc( Abs(xf - Min(Max(s, xf(1)), xf(nf))), 1 )
    f2 = Minloc( Abs(xf - Min(Max(e, xf(1)), xf(nf))), 1 )
    If ( f2 <= f1 ) Then
       If ( f1 < nf ) Then
          f2 = f1 + 1
       Else
          f1 = f2 - 1
       End If
    End If

  End Subroutine coordinate_face_range


  !> Resolve the chunk box inputs into global 1-based face indices chunk_f1/chunk_f2 (cells f1..f2-1 lie inside), aborting if the
  !  box is empty or outside the domain. Index mode (all six inputs negative) maps v to the 0-based face index -v-1.
  Subroutine resolve_chunk_box

    Logical        :: index_mode, ok
    Integer(Int32) :: a, nf(3), n1, n2
    Real(Int64)    :: s(3), e(3)

    If ( .Not. chunked_snapshot ) Return

    index_mode = ( x_start < 0d0 )
    s  = [ x_start, y_start, z_start ]
    e  = [ x_end,   y_end,   z_end   ]
    nf = [ nx_global, ny_global, nz_global ]

    Do a = 1, 3
       If ( index_mode ) Then
          n1 = Nint(-s(a)) - 1
          n2 = Nint(-e(a)) - 1
          ok = ( n1 >= 0 .And. n2 <= nf(a)-1 .And. n1 < n2 )
          chunk_f1(a) = n1 + 1
          chunk_f2(a) = n2 + 1
       Else
          Select Case (a)
          Case (1);  Call coordinate_face_range(x_global, nf(a), s(a), e(a), chunk_f1(a), chunk_f2(a), ok)
          Case (2);  Call coordinate_face_range(y_global, nf(a), s(a), e(a), chunk_f1(a), chunk_f2(a), ok)
          Case (3);  Call coordinate_face_range(z_global, nf(a), s(a), e(a), chunk_f1(a), chunk_f2(a), ok)
          End Select
       End If
       If ( .Not. ok ) Then
          If ( myid == 0 ) Write(*,'(A,I0,A)') ' ERROR: chunk box is empty or outside the domain in direction ', a, &
               ' (index mode needs 0 <= start < end <= n_faces-1)'
          Call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
       End If
    End Do

    If ( myid == 0 ) Then
       Write(*,'(A,3(1X,I0,A,I0))') '   chunk box face indices (1-based)  :', &
            chunk_f1(1), ':', chunk_f2(1), chunk_f1(2), ':', chunk_f2(2), chunk_f1(3), ':', chunk_f2(3)
       Write(*,'(A,3(1X,ES11.4,A,ES11.4))') '   chunk box extent                  :', &
            x_global(chunk_f1(1)), ':', x_global(chunk_f2(1)), y_global(chunk_f1(2)), ':', y_global(chunk_f2(2)), &
            z_global(chunk_f1(3)), ':', z_global(chunk_f2(3))
    End If

  End Subroutine resolve_chunk_box


  !> Stagger flags st and global index range gl:gh of field kind ic (0 scalar, 1/2/3 = U/V/W) over the face box f1:f2:
  !  a staggered direction holds faces f1..f2, a centred one holds cells f1..f2-1
  Subroutine field_box(ic, f1, f2, st, gl, gh)

    Integer(Int32), Intent(In)  :: ic, f1(3), f2(3)
    Integer(Int32), Intent(Out) :: st(3), gl(3), gh(3)

    st = 0
    If ( snapshot_location == 0 .And. ic > 0 ) st(ic) = 1
    gl = f1
    gh = f2 - 1 + st

  End Subroutine field_box


  !> Part of the field box gl:gh owned by rank r (first global index lo, counts n; n = 0 when it owns none): cells ig1..ig2-2 in x,
  !  kg1..kg2-2 in z, all of y; the closing face of a staggered direction belongs to the last rank
  Subroutine box_extent(r, st, gl, gh, lo, n)

    Integer(Int32), Intent(In)  :: r, st(3), gl(3), gh(3)
    Integer(Int32), Intent(Out) :: lo(3), n(3)
    Integer(Int32) :: own_lo(3), own_hi(3)

    own_lo = [ ig1_global(r), 1, kg1_global(r) ]
    own_hi = [ ig2_global(r)-2, nym_global, kg2_global(r)-2 ]
    If ( st(1) == 1 .And. own_hi(1) == nxm_global ) own_hi(1) = nx_global
    If ( st(2) == 1 ) own_hi(2) = ny_global
    If ( st(3) == 1 .And. own_hi(3) == nzm_global ) own_hi(3) = nz_global

    lo = Max(gl, own_lo)
    n  = Max(0, Min(gh, own_hi) - lo + 1)

  End Subroutine box_extent


  !> Fill buf with this rank's share lo:lo+n-1 of field A; velocities are averaged to cell centres when snapshot_location = 1
  Subroutine fill_box_block(A, ic, st, lo, n, buf)

    Real(Int64),    Intent(In)  :: A(:,:,:)
    Integer(Int32), Intent(In)  :: ic, st(3), lo(3), n(3)
    Real(Int64),    Intent(Out) :: buf(:,:,:)

    Integer(Int32) :: i, j, k, ia, ja, ka, ox, oy, oz
    Logical        :: centred

    centred = ( snapshot_location == 1 )
    ox = 2 - ig1_global(myid) - st(1)
    oy = 1 - st(2)
    oz = 2 - kg1_global(myid) - st(3)

    Do k = 1, n(3)
       ka = lo(3) + k - 1 + oz
       Do j = 1, n(2)
          ja = lo(2) + j - 1 + oy
          Do i = 1, n(1)
             ia = lo(1) + i - 1 + ox
             If ( centred .And. ic == 1 ) Then
                buf(i,j,k) = 0.5d0*( A(ia-1,ja,ka) + A(ia,ja,ka) )
             Else If ( centred .And. ic == 2 ) Then
                buf(i,j,k) = 0.5d0*( A(ia,ja-1,ka) + A(ia,ja,ka) )
             Else If ( centred .And. ic == 3 ) Then
                buf(i,j,k) = 0.5d0*( A(ia,ja,ka-1) + A(ia,ja,ka) )
             Else
                buf(i,j,k) = A(ia,ja,ka)
             End If
          End Do
       End Do
    End Do

  End Subroutine fill_box_block


  !> Collective: gather one field of the box f1:f2 onto rank 0 and append it to the stream unit_no open there
  Subroutine write_box_field(unit_no, A, ic, f1, f2)

    Integer(Int32), Intent(In) :: unit_no, ic, f1(3), f2(3)
    Real(Int64),    Intent(In) :: A(:,:,:)

    Integer(Int32) :: st(3), gl(3), gh(3), lo(3), n(3), r
    Real(Int64), Allocatable :: buf(:,:,:), gbuf(:,:,:), rbuf(:,:,:)

    Call field_box(ic, f1, f2, st, gl, gh)
    Call box_extent(myid, st, gl, gh, lo, n)
    If ( All(n > 0) ) Then
       Allocate( buf(n(1),n(2),n(3)) )
       Call fill_box_block(A, ic, st, lo, n, buf)
    End If

    If ( myid == 0 ) Then
       Allocate( gbuf(gh(1)-gl(1)+1, gh(2)-gl(2)+1, gh(3)-gl(3)+1) )
       Do r = 0, nprocs-1
          Call box_extent(r, st, gl, gh, lo, n)
          If ( Any(n <= 0) ) Cycle
          If ( r == 0 ) Then
             gbuf(lo(1)-gl(1)+1:lo(1)-gl(1)+n(1), lo(2)-gl(2)+1:lo(2)-gl(2)+n(2), lo(3)-gl(3)+1:lo(3)-gl(3)+n(3)) = buf
          Else
             Allocate( rbuf(n(1),n(2),n(3)) )
             Call Mpi_recv(rbuf, n(1)*n(2)*n(3), Mpi_real8, r, r, MPI_COMM_WORLD, istat, ierr)
             gbuf(lo(1)-gl(1)+1:lo(1)-gl(1)+n(1), lo(2)-gl(2)+1:lo(2)-gl(2)+n(2), lo(3)-gl(3)+1:lo(3)-gl(3)+n(3)) = rbuf
             Deallocate( rbuf )
          End If
       End Do
       Write(unit_no) gbuf
       Deallocate( gbuf )
    Else If ( All(n > 0) ) Then
       Call Mpi_send(buf, n(1)*n(2)*n(3), Mpi_real8, 0, myid, MPI_COMM_WORLD, ierr)
    End If

  End Subroutine write_box_field


  !> Collective: write the face box f1:f2 of every field to fname. Layout (little-endian stream): snap_magic, snapshot_location,
  !  step, nfields (Int32); t; f1, f2; per direction x_global(f1:f2) then xm_global(f1:f2-1); per field name(8), st, gl, gh;
  !  then the field blocks in the same order, each (gh-gl+1) column-major
  Subroutine write_box_snapshot(fname, f1, f2)

    Character(*),   Intent(In) :: fname
    Integer(Int32), Intent(In) :: f1(3), f2(3)

    Character(8)   :: names(8)
    Integer(Int32) :: fid(8), nf, q, st(3), gl(3), gh(3), unit_no, ios
    Integer(Int64) :: fsize

    names(1) = 'U';  names(2) = 'V';  names(3) = 'W';  names(4) = 'P'
    fid = [ 1, 2, 3, 4, 0, 0, 0, 0 ]
    nf  = 4
    If ( sediment_flag >= 1 )   Then;  nf = nf + 1;  names(nf) = 'C';     fid(nf) = 5;  End If
    If ( sgs_model /= 0 )       Then;  nf = nf + 1;  names(nf) = 'nu_t';  fid(nf) = 6;  End If
    If ( boussinesq_flag >= 1 ) Then;  nf = nf + 1;  names(nf) = 'T';     fid(nf) = 7;  End If
    If ( vof_active >= 1 )      Then;  nf = nf + 1;  names(nf) = 'Cvof';  fid(nf) = 8;  End If

    If ( myid == 0 ) Then
       Write(*,'(A,A,A)') '  Writing snapshot: ', Trim(Adjustl(fname)), ' ...'
       Open(newunit=unit_no, file=Trim(Adjustl(fname)), access='stream', form='unformatted', action='write', &
            status='replace', convert='little_endian', iostat=ios)
       If ( ios /= 0 ) Then
          Write(*,'(A,A)') ' ERROR: cannot open snapshot file ', Trim(Adjustl(fname))
          Call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
       End If
       Write(unit_no) snap_magic, snapshot_location, istep + nstep_init, nf
       Write(unit_no) t
       Write(unit_no) f1, f2
       Write(unit_no) x_global(f1(1):f2(1)), xm_global(f1(1):f2(1)-1)
       Write(unit_no) y_global(f1(2):f2(2)), ym_global(f1(2):f2(2)-1)
       Write(unit_no) z_global(f1(3):f2(3)), zm_global(f1(3):f2(3)-1)
       Do q = 1, nf
          Call field_box(Merge(fid(q), 0, fid(q) <= 3), f1, f2, st, gl, gh)
          Write(unit_no) names(q), st, gl, gh
       End Do
    End If

    Do q = 1, nf
       Select Case ( fid(q) )
       Case (1);  Call write_box_field(unit_no, U,       1, f1, f2)
       Case (2);  Call write_box_field(unit_no, V,       2, f1, f2)
       Case (3);  Call write_box_field(unit_no, W,       3, f1, f2)
       Case (4);  Call write_box_field(unit_no, P,       0, f1, f2)
       Case (5);  Call write_box_field(unit_no, Cscal,   0, f1, f2)
       Case (6);  Call write_box_field(unit_no, nu_t,    0, f1, f2)
       Case (7);  Call write_box_field(unit_no, Tscal,   0, f1, f2)
       Case (8);  Call write_box_field(unit_no, Cvof_io, 0, f1, f2)
       End Select
    End Do

    If ( myid == 0 ) Then
       Close(unit_no)
       Inquire(file=Trim(Adjustl(fname)), size=fsize)
       Write(*,'(A,A,A,F8.3,A)') '  Wrote ', Trim(Adjustl(fname)), '  (', fsize/1e9, ' GB)'
    End If

  End Subroutine write_box_snapshot


  ! Input: U,V,W,x,y,z,xm,ym,zm
  ! Output: fileout
  Subroutine output_data

    Character(200)   :: fname, pfname
    Character(8)     :: ext, pext
    Integer  (Int64) :: fsize
    Logical          :: dirExists
    Logical          :: save_now

    ! nsave>0: step-count cadence; nsave<0: physical-time cadence, fires once t reaches tsave_next (kept on tsave multiples by dt snapping in compute_time_step_RK3)
    If ( nsave > 0 ) Then
       save_now = ( Mod(istep,nsave) == 0 )
    Else
       ! tolerance guards against t landing one ULP short of tsave_next after dt-snapping, which would otherwise stall tsave_next forever
       save_now = ( t >= tsave_next - 1d-10 )
       If ( save_now ) tsave_next = tsave_next + tsave
    End If

    If ( save_now ) then

       If ( vof_active >= 1 ) Call vof_sync_io_hook()
       If ( sediment_flag >= 1 ) Then
          !$acc update host(Cscal)
       End If
       If ( boussinesq_flag >= 1 ) Then
          !$acc update host(Tscal)
       End If
#ifdef GPU_POISSON
       !$acc update host(U,V,W,P)
       If ( sgs_model /= 0 ) Then
          !$acc update host(nu_t)
       End If
#endif

       Write(ext,'(I8)') istep + nstep_init

       If ( myid==0 ) Then
          !Create a directory to store the data, if needed
          inquire(file=trim("fields")//'/.', exist=dirExists)   !Works for gfortran
          If ( .not. dirExists ) call system('mkdir -p fields/')
       End If

       ! only the staggered full snapshot is the restart format; the centred full snapshot reuses the box writer on the whole domain
       If ( full_snapshot .And. snapshot_location == 1 ) Then
          fname = 'fields/'//Trim(Adjustl(fileout))//'.'//Trim(Adjustl(ext))
          Call write_box_snapshot(fname, [1, 1, 1], [nx_global, ny_global, nz_global])
       End If

       ! with both snapshots on, the chunk gets the Block suffix so it cannot clobber the full file
       If ( chunked_snapshot ) Then
          fname = 'fields/'//Trim(Adjustl(fileout))//'.'//Trim(Adjustl(ext))
          If ( full_snapshot ) fname = 'fields/'//Trim(Adjustl(fileout))//'Block.'//Trim(Adjustl(ext))
          Call write_box_snapshot(fname, chunk_f1, chunk_f2)
       End If

       If ( full_snapshot .And. snapshot_location == 0 ) Then

       ! processor 0 writes the data
       If ( myid==0 ) Then

          fname = 'fields/'//Trim(Adjustl(fileout))//'.'//Trim(Adjustl(ext))

          Write(*,'(A,A,A)') '  Writing snapshot: ', Trim(Adjustl(fname)), ' ...'
          Open(1,file=fname,access='stream',form='unformatted',action='write')

          ! mesh
          Write(1) Shape(x_global), x_global
          Write(1) Shape(y_global), y_global
          Write(1) Shape(z_global), z_global

          Write(1) Shape(xm_global), xm_global
          Write(1) Shape(ym_global), ym_global
          Write(1) Shape(zm_global), zm_global

       End If

       Call write_distributed_field_block(1, U, nx_global,  nyg_global, nzg_global, .True.,  .False.)
       Call write_distributed_field_block(1, V, nxg_global, ny_global,  nzg_global, .False., .False.)
       Call write_distributed_field_block(1, W, nxg_global, nyg_global, nz_global,  .False., .True.)
       Call write_distributed_field_block(1, P, nxg_global, nyg_global, nzg_global, .False., .False.)

       ! C (scalar concentration, cell-centred)
       If ( sediment_flag >= 1 ) Then
          Call write_distributed_field_block(1, Cscal, nxg_global, nyg_global, nzg_global, .False., .False.)
       End If

       ! nu_t (SGS turbulent viscosity, cell-centred) — only written when LES is active
       If ( sgs_model /= 0 ) Then
          Call write_distributed_field_block(1, nu_t, nxg_global, nyg_global, nzg_global, .False., .False.)
       End If

       ! T (Boussinesq temperature, cell-centred)
       If ( boussinesq_flag >= 1 ) Then
          Call write_distributed_field_block(1, Tscal, nxg_global, nyg_global, nzg_global, .False., .False.)
       End If

       ! VOF liquid fraction (cell-centred), always the last field block
       If ( vof_active >= 1 ) Then
          Call write_distributed_field_block(1, Cvof_io, nxg_global, nyg_global, nzg_global, .False., .False.)
       End If

       ! close file and report size
       If (myid==0) Then
          Close(1)
          Inquire(file=Trim(Adjustl(fname)), size=fsize)
          If (fsize > 0_Int64) Then
             Write(*,'(A,A,A,F6.2,A)') '  Wrote ', Trim(Adjustl(fname)), &
                  '  (', fsize/1e9, ' GiB)'
          End If
       End If

       End If

       ! Particle state (src/particles.f90): separate latest-state restart file plus a timestamped snapshot for ParaView
       ! (postProcessing/generate_particles_xmf.py). write_particle_snapshot is collective and all ranks must agree on the filename.
       If ( particles_active >= 1 ) Then
          Call write_particle_restart
          Write(pext,'(I8)') istep + nstep_init
          pfname = 'fields/'//Trim(Adjustl(fileout))//'_particles.'//Trim(Adjustl(pext))
          Call write_particle_snapshot(pfname)
       End If

       Call write_tavg(istep + nstep_init)

    End If

  End Subroutine output_data


  !> Open the restart file as unit 2 on rank 0, positioned after the six mesh blocks (size header, then the coordinates)
  Subroutine open_restart_skip_mesh

    Integer(Int32) :: m, n, nn(6)
    Real(Int64), Allocatable :: tmp(:)

    If ( myid /= 0 ) Return
    Open(2, file=Trim(Adjustl(filein)), access='stream', form='unformatted', action='read')
    nn = (/ nx_global, ny_global, nz_global, nxm_global, nym_global, nzm_global /)
    Do m = 1, 6
       Allocate( tmp(nn(m)) )
       Read(2) n
       Read(2) tmp
       Deallocate( tmp )
    End Do

  End Subroutine open_restart_skip_mesh


  !> Read one block of a restart file into F, skipping the mesh and every block before it. Layout of the file (output_data): U, V, W,
  !  P, then C (sediment_flag >= 1), nu_t (LES), T (boussinesq_flag >= 1) and the VOF fraction. The optional blocks are located from THIS
  !  run's flags, not from what the file was written with: restarting with a flag flipped relative to the writing run desyncs the skip
  !  and misreads everything after it. blk: 1 = P, 2 = C, 3 = T, 4 = VOF fraction.
  Subroutine read_restart_block(blk, F)

    Integer(Int32), Intent(In) :: blk
    Real(Int64), Dimension(:,:,:), Intent(InOut) :: F

    Logical :: present_blk(4)
    Integer(Int32) :: m

    Call open_restart_skip_mesh
    ! U, V, W (already consumed by read_input_data)
    Call skip_distributed_field_block(2, nyg, .True.,  .False.)
    Call skip_distributed_field_block(2, ny,  .False., .False.)
    Call skip_distributed_field_block(2, nyg, .False., .True.)
    present_blk = (/ .True., sediment_flag >= 1, boussinesq_flag >= 1, .True. /)
    Do m = 1, blk
       If ( m == 3 .And. sgs_model /= 0 ) Call skip_distributed_field_block(2, nyg, .False., .False.)   ! nu_t sits between C and T
       If ( .Not. present_blk(m) ) Cycle
       If ( m == blk ) Then
          Call read_distributed_field_block(2, F, .False., .False.)
       Else
          Call skip_distributed_field_block(2, nyg, .False., .False.)
       End If
    End Do
    If ( myid == 0 ) Close(2)

  End Subroutine read_restart_block


  !> Read scalar C from a restart file
  Subroutine read_scalar_restart

    Call read_restart_block(2, Cscal)
    Cscal_o = Cscal

  End Subroutine read_scalar_restart


  !> Read the pressure block of a restart file into P (written there from Pst when ibm_pressure > 0)
  Subroutine read_pressure_restart

    Call read_restart_block(1, P)

  End Subroutine read_pressure_restart


  !> Read the VOF fraction into Cvof_io and the pressure history of the two-fluid step into P from a restart file
  Subroutine read_vof_restart

    Call read_restart_block(1, P)
    Call read_restart_block(4, Cvof_io)

  End Subroutine read_vof_restart


  !> Read temperature T from a restart file
  Subroutine read_temperature_restart

    Call read_restart_block(3, Tscal)
    Tscal_o = Tscal
    ! Device copy deferred: Tscal/Tscal_o are not yet `!$acc enter data create`d at this point in startup (that happens later, unconditionally, in initialization.f90), so updating the device here would be undefined behavior.

  End Subroutine read_temperature_restart


End Module input_output
