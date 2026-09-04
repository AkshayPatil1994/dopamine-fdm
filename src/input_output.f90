! Module for I/O
Module input_output

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use mpi
  Use genGridAndIC    

  ! prevent implicit typing
  Implicit None

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
    Namelist /DOMAIN/ nx, ny, nz, Lx, Ly, Lz, alpha_grid, grid_type, p_row, p_col

    Namelist /PHYSICS/ nu, dPdx, dPdz, Ub_x, Ub_z, T_wave_x, T_wave_z, &
                       phi_wave_x, phi_wave_z, sgs_model, Cs_vreman, flat_wall_model_flag, &
                       z0_ylo, z0_yhi, z0h_ylo, z0h_yhi, &
                       flow_forcing_mode, Ub_target, advection_scheme

    Namelist /NUMERICS/ dt, nsteps, nsave, nmonitor, sim_end_time, tsave, &
                        cfl_adaptive, cfl_target, cfl_safety, dt_min, dt_max

    Namelist /BOUNDARY_CONDITIONS/ bc_face_ylo, bc_face_yhi, x_bc_type, y_bc_type

    Namelist /INFLOW/ inflow_type, inflow_Uconst, inflow_profile_file, inflow_temperature_file, &
                       sem_profile_format, sem_Lscale_ratio_y, sem_Lscale_ratio_z, &
                       sem_n_eddies, sem_length_scale, sem_seed, &
                       sem_ensemble_samples, sem_ensemble_periods, &
                       sem_sigma_file, sem_eddy_placement, sem_use_esem, &
                       sem_divergence_free, sem_wall_damping, sem_wall_damping_Aplus, &
                       inflow_recycle_file, inflow_recycle_loop, inflow_recycle_t_offset, &
                       inflow_recycle_shift_z, inflow_recycle_seed

    Namelist /IBM/ ibm_input_mode, ibm_wall_model_flag, &
                   ibm_sdf_file, ibm_objid_file, ks, nks, nsampling, ibm_surface_nsampling, &
                   ibm_T_bc_type, ibm_T_wall, ibm_z0

    Namelist /INITIAL_CONDITIONS/ Utarget, nstep_init, restart, t_start, &
                                   scalar_restart, ic_type, noise_percent

    Namelist /IO/ filein, fileout

    Namelist /SEDIMENT/ sediment_flag, d_s, rho_s, rho_f, grav, Sc, Sc_t, C_ref, sed_bc_bot, &
                         C_ic_type, C_ic_height

    Namelist /BOUSSINESQ/ boussinesq_flag, beta_T, T_ref, grav, Pr, Pr_t, &
                           T_bc_bot, T_bc_top, T_wall_bot, T_wall_top, &
                           T_ic_type, T_ic_grad, ibm_T_bc_type, ibm_T_wall

    Namelist /STATISTICS/ rsb_active, rsb_freq, rsb_nstart, rsb_hom_dir, rsb_fileout, &
                          n_slices, slice_freq, slice_dir, slice_pos, slice_comps, slice_fileout, &
                          n_lines,  line_freq,  line_dir,  line_pos1,  line_pos2,   &
                          line_start, line_end, line_comps, line_fileout

    Namelist /UAV/ uav_active, uav_xc, uav_yc, uav_zc, uav_disk_radius, &
                   uav_n_r, uav_n_theta, uav_hover_thrust, uav_kernel_ncell, &
                   uav_path_active, uav_path_file, uav_thrust_active, uav_thrust_file

    Namelist /TI_RESCALE/ ti_rescale_active, ti_rescale_x, ti_rescale_nstart, &
                          ti_rescale_freq, ti_rescale_relax, ti_rescale_clip, ti_rescale_abs_clip, &
                          ti_rescale_filter_alpha, ti_rescale_deadband, ti_rescale_relax_min, &
                          ti_rescale_u_active, ti_rescale_u_relax, ti_rescale_u_relax_min, &
                          ti_rescale_u_clip, ti_rescale_u_abs_clip, ti_rescale_u_deadband

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

    fname = 'input_parameters'
    If ( Present(input_file) ) fname = input_file

    ! ---- Only rank 0 reads the file ----------------------------------
    If ( myid==0 ) Then

       Open(newunit=unit_in, file=fname, status='old',    &
            action='read', iostat=ios)
       If ( ios /= 0 ) Stop 'ERROR: cannot open input parameters file: '//fname

       Rewind(unit_in)
       Read(unit_in, nml=DOMAIN,              iostat=ios)
       If (ios /= 0) Stop 'ERROR: &DOMAIN missing or failed to parse (check for unrecognized variable names)'

       Rewind(unit_in)
       Read(unit_in, nml=PHYSICS,             iostat=ios)
       If (ios /= 0) Stop 'ERROR: &PHYSICS missing or failed to parse (check for unrecognized variable names)'

       Rewind(unit_in)
       Read(unit_in, nml=NUMERICS,            iostat=ios)
       If (ios /= 0) Stop 'ERROR: &NUMERICS missing or failed to parse (check for unrecognized variable names)'

       If ( namelist_group_present(unit_in, 'BOUNDARY_CONDITIONS') ) Then
          Rewind(unit_in)
          Read(unit_in, nml=BOUNDARY_CONDITIONS, iostat=ios)
          If (ios /= 0) Stop 'ERROR: &BOUNDARY_CONDITIONS present but failed to parse (check variable names)'
       Else
          Write(*,'(A)') ' INFO: no &BOUNDARY_CONDITIONS found, using defaults'
       End If

       If ( namelist_group_present(unit_in, 'INFLOW') ) Then
          Rewind(unit_in)
          Read(unit_in, nml=INFLOW,              iostat=ios)
          If (ios /= 0) Stop 'ERROR: &INFLOW present but failed to parse (check variable names)'
       Else If ( x_bc_type == 1 ) Then
          Write(*,'(A)') ' INFO: no &INFLOW found, using defaults (inflow_type=0, constant uniform flow)'
       End If

       If ( namelist_group_present(unit_in, 'IBM') ) Then
          Rewind(unit_in)
          Read(unit_in, nml=IBM,                 iostat=ios)
          If (ios /= 0) Stop 'ERROR: &IBM present but failed to parse (check variable names)'
       Else
          Write(*,'(A)') ' INFO: no &IBM found, using defaults (no IBM body)'
          ibm_input_mode      = 0
          ibm_wall_model_flag = 0
          ibm_sdf_file        = 'SDF_in'
          ks                  = 0d0
          nks                 = 0
          nsampling           = 0
          ibm_surface_nsampling = 0
       End If

       Rewind(unit_in)
       Read(unit_in, nml=INITIAL_CONDITIONS,  iostat=ios)
       If (ios /= 0) Stop 'ERROR: &INITIAL_CONDITIONS missing or failed to parse (check for unrecognized variable names)'

       Rewind(unit_in)
       Read(unit_in, nml=IO,                  iostat=ios)
       If (ios /= 0) Stop 'ERROR: &IO missing or failed to parse (check for unrecognized variable names)'

       If ( namelist_group_present(unit_in, 'SEDIMENT') ) Then
          Rewind(unit_in)
          Read(unit_in, nml=SEDIMENT,            iostat=ios)
          If (ios /= 0) Stop 'ERROR: &SEDIMENT present but failed to parse (check variable names)'
       Else
          Write(*,'(A)') ' INFO: no &SEDIMENT found, scalar transport disabled'
       End If
       grav_sediment = grav

       If ( namelist_group_present(unit_in, 'BOUSSINESQ') ) Then
          Rewind(unit_in)
          Read(unit_in, nml=BOUSSINESQ,          iostat=ios)
          If (ios /= 0) Stop 'ERROR: &BOUSSINESQ present but failed to parse (check variable names)'
          If ( sediment_flag >= 1 .And. grav /= grav_sediment ) Then
             Write(*,'(A)') ' WARNING: &SEDIMENT and &BOUSSINESQ specify different grav, using &BOUSSINESQ value'
          End If
       Else
          Write(*,'(A)') ' INFO: no &BOUSSINESQ found, thermal stratification disabled'
       End If

       If ( namelist_group_present(unit_in, 'STATISTICS') ) Then
          Rewind(unit_in)
          Read(unit_in, nml=STATISTICS,          iostat=ios)
          If (ios /= 0) Stop 'ERROR: &STATISTICS present but failed to parse (check variable names)'
       Else
          Write(*,'(A)') ' INFO: no &STATISTICS found, Reynolds stress budget disabled'
       End If

       If ( namelist_group_present(unit_in, 'UAV') ) Then
          Rewind(unit_in)
          Read(unit_in, nml=UAV,                 iostat=ios)
          If (ios /= 0) Stop 'ERROR: &UAV present but failed to parse (check variable names)'
       Else
          Write(*,'(A)') ' INFO: no &UAV found, UAV actuator disk disabled'
       End If

       If ( namelist_group_present(unit_in, 'TI_RESCALE') ) Then
          Rewind(unit_in)
          Read(unit_in, nml=TI_RESCALE,          iostat=ios)
          If (ios /= 0) Stop 'ERROR: &TI_RESCALE present but failed to parse (check variable names)'
       End If

       Close(unit_in)

       ! Map short namelist names to the global variable names
       nx_global  = nx;  ny_global  = ny;  nz_global  = nz
       Lx_i       = Lx;  Ly_i       = Ly;  Lz_i       = Lz
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

       ! Derived pressure-gradient references
       dPdx0    = dPdx
       dPdx_t   = dPdx
       dPdx_ref = dPdx
       dPdz_t   = dPdz

       If ( flow_forcing_mode == 1 .And. T_wave_x > 0d0 ) Then
          Stop 'ERROR: flow_forcing_mode=1 (constant mass flux) is incompatible with oscillatory forcing (T_wave_x)'
       End If
       If ( flow_forcing_mode == 1 .And. x_bc_type /= 0 ) Then
          Stop 'ERROR: flow_forcing_mode=1 (constant mass flux) requires periodic streamwise BC (x_bc_type=0)'
       End If

       Write(*,'(A)') ' Input parameters read from namelist file.'
       Write(*,'(A,3I6)')    '   Grid      (nx,ny,nz)        = ', nx, ny, nz
       Write(*,'(A,3F10.4)') '   Domain    (Lx,Ly,Lz)        = ', Lx, Ly, Lz
       Write(*,'(A,I2)')     '   grid_type                   = ', grid_type
       Write(*,'(A,2I4)')    '   p_row, p_col (0,0=auto)     = ', p_row, p_col
       Write(*,'(A,E12.4)')  '   nu                          = ', nu
       Write(*,'(A,I2)')     '   advection_scheme (0=skew-symmetric,1=pure central) = ', advection_scheme
       Write(*,'(A,I2)')     '   flow_forcing_mode (0=dPdx,1=const mass flux) = ', flow_forcing_mode
       If ( flow_forcing_mode == 1 ) Then
          Write(*,'(A,E12.4)') '   Ub_target                   = ', Ub_target
       End If
       Write(*,'(A,E12.4)')  '   dPdx                        = ', dPdx
       Write(*,'(A,E12.4)')  '   dPdz                        = ', dPdz
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
       Write(*,'(A,I2)')     '   sgs_model                   = ', sgs_model
       Write(*,'(A,E12.4)')  '   Cs_vreman                   = ', Cs_vreman
       Write(*,'(A,I2)')     '   flat_wall_model_flag        = ', flat_wall_model_flag
       If ( flat_wall_model_flag == 2 ) Then
          Write(*,'(A,E12.4,A,E12.4)') '   z0_ylo, z0_yhi              = ', z0_ylo,  ', ', z0_yhi
          Write(*,'(A,E12.4,A,E12.4)') '   z0h_ylo, z0h_yhi            = ', z0h_ylo, ', ', z0h_yhi
       End If
       Write(*,'(A,2I3)')    '   BC y-walls (ylo/yhi)        = ', bc_face_ylo, bc_face_yhi
       Write(*,'(A,I2)')     '   x_bc_type (0=periodic,1=inflow/outflow) = ', x_bc_type
       Write(*,'(A,I2)')     '   y_bc_type (0=periodic,1=wall)           = ', y_bc_type
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
          Else
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
             If ( ti_rescale_active == 1 ) Then
                Write(*,'(A,F10.4)') '   ti_rescale_active=1, ti_rescale_x   = ', ti_rescale_x
                Write(*,'(A,I8)')    '   ti_rescale_nstart           = ', ti_rescale_nstart
                Write(*,'(A,I8)')    '   ti_rescale_freq             = ', ti_rescale_freq
                Write(*,'(A,F6.3)')  '   ti_rescale_relax            = ', ti_rescale_relax
                Write(*,'(A,F6.3)')  '   ti_rescale_clip             = ', ti_rescale_clip
                Write(*,'(A,F6.3)')  '   ti_rescale_abs_clip         = ', ti_rescale_abs_clip
                Write(*,'(A,F6.3)')  '   ti_rescale_filter_alpha     = ', ti_rescale_filter_alpha
                Write(*,'(A,F6.3)')  '   ti_rescale_deadband         = ', ti_rescale_deadband
                Write(*,'(A,F6.3)')  '   ti_rescale_relax_min        = ', ti_rescale_relax_min
                If ( ti_rescale_u_active == 1 ) Then
                   Write(*,'(A)')       '   ti_rescale_u_active=1 (mean-profile rescaling enabled)'
                   Write(*,'(A,F6.3)')  '   ti_rescale_u_relax          = ', ti_rescale_u_relax
                   Write(*,'(A,F6.3)')  '   ti_rescale_u_relax_min      = ', ti_rescale_u_relax_min
                   Write(*,'(A,F6.3)')  '   ti_rescale_u_clip           = ', ti_rescale_u_clip
                   Write(*,'(A,F6.3)')  '   ti_rescale_u_abs_clip       = ', ti_rescale_u_abs_clip
                   Write(*,'(A,F6.3)')  '   ti_rescale_u_deadband       = ', ti_rescale_u_deadband
                End If
             End If
          End If
       End If
       Write(*,'(A,I2)')     '   ibm_input_mode              = ', ibm_input_mode
       Write(*,'(A,I2)')     '   ibm_wall_model_flag         = ', ibm_wall_model_flag
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
       If ( boussinesq_flag >= 1 ) Then
          Write(*,'(A,I2)')    '   boussinesq_flag              = ', boussinesq_flag
          Write(*,'(A,E12.4)') '   beta_T                       = ', beta_T
          Write(*,'(A,E12.4)') '   T_ref                        = ', T_ref
          Write(*,'(A,E12.4)') '   grav                         = ', grav
          Write(*,'(A,F6.3)')  '   Pr                           = ', Pr
          Write(*,'(A,F6.3)')  '   Pr_t                         = ', Pr_t
          Write(*,'(A,2I3)')   '   T_bc (bot/top, 0=adiabatic,1=isothermal) = ', T_bc_bot, T_bc_top
       End If
       If ( n_slices > 0 ) Then
          Write(*,'(A,I4)')  '   n_slices                    = ', n_slices
          Write(*,'(A,I8)')  '   slice_freq                  = ', slice_freq
       End If
       If ( n_lines > 0 ) Then
          Write(*,'(A,I4)')  '   n_lines                     = ', n_lines
          Write(*,'(A,I8)')  '   line_freq                   = ', line_freq
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
    Call Mpi_bcast ( grid_type,            1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( Utarget,              1, MPI_real8,     0, MPI_COMM_WORLD, ierr )

    Call Mpi_bcast ( p_row,                1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( p_col,                1, MPI_integer,   0, MPI_COMM_WORLD, ierr )

    Call Mpi_bcast ( ibm_input_mode,       1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( ibm_sdf_file,  Len(ibm_sdf_file),  MPI_character, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( ibm_objid_file, Len(ibm_objid_file), MPI_character, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( ibm_wall_model_flag,  1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( nks_global,           1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( ks,                   1, MPI_real8,     0, MPI_COMM_WORLD, ierr )

    Call Mpi_bcast ( dt,                   1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( nu,                   1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( advection_scheme,     1, MPI_integer,   0, MPI_COMM_WORLD, ierr )

    Call Mpi_bcast ( dPdx,                 1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( dPdx0,                1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( dPdx_t,               1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( dPdz,                 1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( dPdz_t,               1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( dPdx_ref,             1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( Ub_x,                 1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( Ub_z,                 1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( waveOmega_x,          1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( waveOmega_z,          1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( phi_wave_x,           1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( phi_wave_z,           1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( flow_forcing_mode,    1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( Ub_target,            1, MPI_real8,     0, MPI_COMM_WORLD, ierr )

    Call Mpi_bcast ( nstep_init,           1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( t_start,              1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( nsteps,               1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( nsave,                1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
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
    Call Mpi_bcast ( y_bc_type,            1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
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

    Call Mpi_bcast ( n_slices,             1,              MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( slice_freq,           1,              MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( slice_dir,    Len(slice_dir(1))*MAX_PROBES,   MPI_character, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( slice_pos,            MAX_PROBES,     MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( slice_comps,  Len(slice_comps(1))*MAX_PROBES, MPI_character, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( slice_fileout,Len(slice_fileout(1))*MAX_PROBES,MPI_character,0, MPI_COMM_WORLD, ierr )

    Call Mpi_bcast ( n_lines,              1,              MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( line_freq,            1,              MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( line_dir,     Len(line_dir(1))*MAX_PROBES,    MPI_character, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( line_pos1,            MAX_PROBES,     MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( line_pos2,            MAX_PROBES,     MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( line_start,           MAX_PROBES,     MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( line_end,             MAX_PROBES,     MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( line_comps,   Len(line_comps(1))*MAX_PROBES,  MPI_character, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( line_fileout, Len(line_fileout(1))*MAX_PROBES,MPI_character, 0, MPI_COMM_WORLD, ierr )

    Call Mpi_bcast ( ti_rescale_active,    1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( ti_rescale_x,         1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( ti_rescale_nstart,    1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( ti_rescale_freq,      1, MPI_integer,   0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( ti_rescale_relax,     1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( ti_rescale_clip,      1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( ti_rescale_abs_clip,  1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( ti_rescale_filter_alpha, 1, MPI_real8,  0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( ti_rescale_deadband,  1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( ti_rescale_relax_min, 1, MPI_real8,     0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( ti_rescale_u_active,   1, MPI_integer,  0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( ti_rescale_u_relax,    1, MPI_real8,    0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( ti_rescale_u_relax_min,1, MPI_real8,    0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( ti_rescale_u_clip,     1, MPI_real8,    0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( ti_rescale_u_abs_clip, 1, MPI_real8,    0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast ( ti_rescale_u_deadband, 1, MPI_real8,    0, MPI_COMM_WORLD, ierr )

  End Subroutine read_input_parameters

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
       If ( nx_global_f/=nx_global ) Stop 'nx_f/=nx'
       Read(1) x_global

       Read(1) ny_global_f
       If ( ny_global_f/=ny_global ) Stop 'ny_f/=ny'
       Read(1) y_global

       Read(1) nz_global_f
       If ( nz_global_f/=nz_global ) Stop 'nz_f/=nz'
       Read(1) z_global

       Read(1) nxm_global_f
       If ( nxm_global_f/=nxm_global ) Stop 'nxm_f/=nxm'
       Read(1) xm_global

       Read(1) nym_global_f
       If ( nym_global_f/=nym_global ) Stop 'nym_f/=nym'
       Read(1) ym_global

       Read(1) nzm_global_f
       If ( nzm_global_f/=nzm_global ) Stop 'nzm_f/=nzm'
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

       ! nn is the field's global (nx1,nx2,nx3) block shape as written; cross-check
       ! against what this call expects (x/z face-vs-ghost count per is_x_face/
       ! is_z_face, y always the caller's own full local extent since y is never
       ! MPI-decomposed) instead of trusting it blindly -- a mismatch here would
       ! otherwise only surface later as an out-of-bounds/shape-mismatch crash (or,
       ! worse, silently wrong data) inside the ix1:ix2,:,iz1:iz2 slice below.
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

  ! Input: U,V,W,x,y,z,xm,ym,zm
  ! Output: fileout
  Subroutine output_data

    Character(200)   :: fname
    Character(8)     :: ext
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

       ! processor 0 writes the data
       If ( myid==0 ) Then

          !Create a directory to store the data, if needed
          inquire(file=trim("fields")//'/.', exist=dirExists)   !Works for gfortran
          If ( .not. dirExists ) call system('mkdir -p fields/')

          Write(ext,'(I8)') istep + nstep_init

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

  End Subroutine output_data


  ! Read scalar C from a restart file; caller/layout details
  Subroutine read_scalar_restart

    Real(Int64), Allocatable :: skip_U(:,:,:), skip_V(:,:,:), skip_W(:,:,:), skip_P(:,:,:)

    If ( myid == 0 ) Then
       Open(2, file=Trim(Adjustl(filein)), access='stream', &
            form='unformatted', action='read')
       ! Skip mesh arrays
       Block
          Integer(Int32) :: n1, n2, n3, n4, n5, n6
          Real   (Int64), Allocatable :: tmp1(:), tmp2(:), tmp3(:), tmp4(:), &
                                         tmp5(:), tmp6(:)
          Allocate(tmp1(nx_global)); Read(2) n1; Read(2) tmp1
          Allocate(tmp2(ny_global)); Read(2) n2; Read(2) tmp2
          Allocate(tmp3(nz_global)); Read(2) n3; Read(2) tmp3
          Allocate(tmp4(nxm_global)); Read(2) n4; Read(2) tmp4
          Allocate(tmp5(nym_global)); Read(2) n5; Read(2) tmp5
          Allocate(tmp6(nzm_global)); Read(2) n6; Read(2) tmp6
          Deallocate(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)
       End Block
    End If

    ! Skip U, V, W, P blocks (discarded — this restart file's U/V/W/P were already consumed by read_input_data)
    Allocate( skip_U(nx,nyg,nzg), skip_V(nxg,ny,nzg), skip_W(nxg,nyg,nz), skip_P(nxg,nyg,nzg) )
    Call read_distributed_field_block(2, skip_U, .True.,  .False.)
    Call read_distributed_field_block(2, skip_V, .False., .False.)
    Call read_distributed_field_block(2, skip_W, .False., .True.)
    Call read_distributed_field_block(2, skip_P, .False., .False.)
    Deallocate( skip_U, skip_V, skip_W, skip_P )

    ! C block
    Call read_distributed_field_block(2, Cscal, .False., .False.)

    If ( myid == 0 ) Close(2)
    Cscal_o = Cscal

  End Subroutine read_scalar_restart


  ! Read temperature T from a restart file; skips U,V,W,P and any C/nu_t blocks ahead of it
  Subroutine read_temperature_restart

    Real(Int64), Allocatable :: skip_U(:,:,:), skip_V(:,:,:), skip_W(:,:,:), skip_P(:,:,:)

    If ( myid == 0 ) Then
       Open(2, file=Trim(Adjustl(filein)), access='stream', &
            form='unformatted', action='read')
       Block
          Integer(Int32) :: n1, n2, n3, n4, n5, n6
          Real   (Int64), Allocatable :: tmp1(:), tmp2(:), tmp3(:), tmp4(:), &
                                         tmp5(:), tmp6(:)
          Allocate(tmp1(nx_global)); Read(2) n1; Read(2) tmp1
          Allocate(tmp2(ny_global)); Read(2) n2; Read(2) tmp2
          Allocate(tmp3(nz_global)); Read(2) n3; Read(2) tmp3
          Allocate(tmp4(nxm_global)); Read(2) n4; Read(2) tmp4
          Allocate(tmp5(nym_global)); Read(2) n5; Read(2) tmp5
          Allocate(tmp6(nzm_global)); Read(2) n6; Read(2) tmp6
          Deallocate(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)
       End Block
    End If

    ! U,V,W,P blocks (discarded — already consumed by read_input_data)
    Allocate( skip_U(nx,nyg,nzg), skip_V(nxg,ny,nzg), skip_W(nxg,nyg,nz), skip_P(nxg,nyg,nzg) )
    Call read_distributed_field_block(2, skip_U, .True.,  .False.)
    Call read_distributed_field_block(2, skip_V, .False., .False.)
    Call read_distributed_field_block(2, skip_W, .False., .True.)
    Call read_distributed_field_block(2, skip_P, .False., .False.)

    ! C block (only present if sediment was active when the file was written) — discarded
    If ( sediment_flag >= 1 ) Call read_distributed_field_block(2, skip_P, .False., .False.)

    ! nu_t block (only present if LES was active when the file was written) — discarded
    If ( sgs_model /= 0 ) Call read_distributed_field_block(2, skip_P, .False., .False.)

    Deallocate( skip_U, skip_V, skip_W, skip_P )

    ! T block
    Call read_distributed_field_block(2, Tscal, .False., .False.)

    If ( myid == 0 ) Close(2)
    Tscal_o = Tscal
    ! Device copy deferred: Tscal/Tscal_o are not yet `!$acc enter data create`d at this point in startup (that happens later, unconditionally, in initialization.f90), so updating the device here would be undefined behavior.

  End Subroutine read_temperature_restart


End Module input_output
