! Module with all shared global variables
Module global

  ! General Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use, Intrinsic :: iso_c_binding

  ! prevent implicit typing
  Implicit None

  ! FFTW
  Include 'fftw3.f03'

  ! Declarations

  ! step number
  Integer(Int32) :: istep, rk_step
  Real   (Int64) :: time1, time2, time_wall_start

  ! constants
  Real(Int64) :: pi = 4d0*datan(1d0)
  !$acc declare create(pi)

  ! files
  Character(200) :: filein, fileout
  Integer(Int32) :: nsave, nmonitor

  ! domain size
  Real(Int64) :: Lx, Lz, Ly, Lxp, Lzp

  ! steps
  Integer(Int32) :: nsteps, nstep_init
  Real   (Int64) :: dt, t

  ! Time-based stopping/save control: sim_end_time/tsave override nsteps/nsave when negative; tsave_next is the next due save time
  Real   (Int64) :: sim_end_time = 1d30
  Real   (Int64) :: tsave        = 1d30
  Real   (Int64) :: tsave_next   = 0d0

  ! viscosity
  Real(Int64) :: nu
  !$acc declare create(nu)

  ! convective-term discretization: 0=skew-symmetric (default, recommended -- suppresses
  ! aliasing-driven ringing near sharp/immersed-boundary gradients), 1=pure divergence-form central
  ! Host-only: consumed once per RHS call to derive a blend weight, never referenced on device.
  Integer(Int32) :: advection_scheme = 0

  ! global face points
  Integer(Int32) :: nx_global, ny_global, nz_global

  ! 2decomp&fft process grid (0,0 = let the library auto-factorize nprocs)
  Integer(Int32) :: p_row = 0, p_col = 0

  ! global roughness points
  Real(Int64) :: ks
  Integer(Int32) :: nks_global

  ! restart flag: 0=fresh run, 1=hot-start from filein
  Integer(Int32) :: restart = 0
  ! scalar_restart: 1=read C from restart file (default), 0=use IC even when restart=1
  Integer(Int32) :: scalar_restart = 1

  ! IBM control flags: ibm_input_mode 0=no body,1=SDF ghost-cell IBM; ibm_wall_model_flag 0=DNS no-slip,1=log-law EQWM
  Integer(Int32) :: ibm_input_mode      = 0            ! default: no IBM body
  Character(200) :: ibm_sdf_file        = 'SDF_in'     ! cell-centre SDF file
  Character(200) :: ibm_objid_file      = ''           ! optional per-solid ID field (GenSDF sdfp_objid.bin); '' = single uniform IBM condition
  Integer(Int32) :: ibm_wall_model_flag = 0             ! default: DNS no-slip

  !  y-wall boundary condition type
  !    1 = Dirichlet (no-slip)  2 = Neumann (free-slip)
  Integer(Int32) :: bc_face_ylo = 1, bc_face_yhi = 1

  ! local face points
  Integer(Int32) :: nx, ny, nz
  !$acc declare create(nx,ny,nz)

  ! global center points
  Integer(Int32) :: nxm_global, nym_global, nzm_global

  ! local center points
  Integer(Int32) :: nxm, nym, nzm

  ! global center points + ghost cells
  Integer(Int32) :: nxg_global, nyg_global, nzg_global

  ! local center points + ghost cells
  Integer(Int32) :: nxg, nyg, nzg
  !$acc declare create(nxg,nyg,nzg)

  ! global grid at face points
  Real(Int64), Allocatable, Dimension(:) :: x_global, y_global, z_global

  ! local grid at face points
  Real(Int64), Allocatable, Dimension(:) :: x, y, z

  ! global grid at middle points
  Real(Int64), Allocatable, Dimension(:) :: xm_global, ym_global, zm_global

  ! local grid at middle points
  Real(Int64), Allocatable, Dimension(:) :: xm, ym, zm

  ! global grid at middle points + ghost cells
  Real(Int64), Allocatable, Dimension(:) :: xg_global, yg_global, zg_global

  ! local grid at middle points + ghost cells
  Real(Int64), Allocatable, Dimension(:) :: xg, yg, zg

  ! middle points for yg->yg_m and yg_m->yg_mm
  Real(Int64), Allocatable, Dimension(:) :: yg_m, yg_mm

  ! local velocities and pressure
  Real(Int64), Allocatable, Dimension(:,:,:) :: U,V,W,P
  Real(Int64), Allocatable, Dimension(:,:,:) :: Uo,Vo,Wo,Po

  ! local auxiliary 
  Real(Int64), Allocatable, Dimension(:,:,:) :: term_1, term_2, term

  ! local rhs for velocities and pressure
  Real(Int64), Allocatable, Dimension(:,:,:) :: rhs_p
  Real(Int64), Allocatable, Dimension(:,:,:) :: rhs_uo, rhs_vo, rhs_wo

  ! local rhs for pressure in Fourier
  ! rhs_p_hat storage order
  Complex(Int64), Dimension(:,:,:), Allocatable :: rhs_p_hat

  ! local auxiliary arrays for MPI_sendrev boundary conditions
  Real(Int64), Allocatable, Dimension(:,:,:) :: buffer_ui, buffer_vi, buffer_ci
  Real(Int64), Allocatable, Dimension(:,:)   :: buffer_ue, buffer_ve, buffer_we, buffer_wi, buffer_ce
  Real(Int64), Allocatable, Dimension(:,:)   :: buffer_p

  ! local auxiliary arrays for MPI_isend/irecv interior planes; 3rd dim: 1=+z exchange, 2=-z exchange, issued concurrently
  Real(Int64), Allocatable, Dimension(:,:,:) :: buffer_us, buffer_ur
  Real(Int64), Allocatable, Dimension(:,:,:) :: buffer_vs, buffer_vr
  Real(Int64), Allocatable, Dimension(:,:,:) :: buffer_ws, buffer_wr
  Real(Int64), Allocatable, Dimension(:,:)   :: buffer_ps, buffer_pr
  
  ! local pencil work arrays for the Poisson pencil-transpose chain (2decomp&fft); rhs_p_hat below (y-pencil, post z-FFT) is shared with the GPU_POISSON path
  Real   (Int64), Allocatable, Dimension(:,:,:) :: poisson_y_r   ! y-pencil, real: interfaces with rhs_p
  Real   (Int64), Allocatable, Dimension(:,:,:) :: poisson_x_r   ! x-pencil, real: DCT-IV path only
  Complex(Int64), Allocatable, Dimension(:,:,:) :: poisson_x_c   ! x-pencil, complex
  Complex(Int64), Allocatable, Dimension(:,:,:) :: poisson_y_c   ! y-pencil, complex: also the Zgtsv operand
  Complex(Int64), Allocatable, Dimension(:,:,:) :: poisson_z_c   ! z-pencil, complex

  ! Fourier points and wave numbers
  Integer(C_INTPTR_T) :: nxp_global, nzp_global
  Integer(Int32)      :: nzp   ! local z-slab count for physical-space rhs_p ghost fill (old 1D z-slab convention; unrelated to the pencil-transpose chain above, unchanged by the 2D-pencil Poisson port until physical-space arrays move to the pencil layout too)
  Integer(C_INTPTR_T) :: mx_global, mz_global
  Integer(C_INTPTR_T) :: mx, mz   ! this rank's local mode-count range in the y-pencil after the transpose chain (== mx_global/mz_global whenever x/z aren't split, e.g. nprocs==1)
  Real   (Int64)      :: dx, dz
  Real   (Int64), Dimension(:), Allocatable :: kxx, kyy, kzz

  ! local (non-MPI) FFTW plans for the post-transpose 1-D transforms
  Type(C_PTR) :: plan_fx_fwd, plan_fx_inv   ! complex 1-D FFT in x (periodic path)
  Type(C_PTR) :: plan_fy_fwd, plan_fy_inv   ! complex 1-D FFT in y (periodic-y path only, y_bc_type==0)
  Type(C_PTR) :: plan_fz_fwd, plan_fz_inv   ! complex 1-D FFT in z (both paths)
  Type(C_PTR) :: plan_dct                    ! real 1-D DCT-IV in x (inflow/outflow path)

  ! streamwise (x) pressure/velocity BC selector: 0=periodic, 1=inflow/outflow
  Integer(Int32) :: x_bc_type  = 0

  ! wall-normal (y) pressure/velocity BC selector: 0=periodic, 1=wall (uses bc_face_ylo/yhi as today)
  Integer(Int32) :: y_bc_type  = 1

  ! &INFLOW streamwise inflow condition (x_bc_type==1 only): inflow_type 0=constant, 1=SEM, 2=recycled precursor slice
  Integer(Int32) :: inflow_type        = 0
  Real   (Int64) :: inflow_Uconst      = 0d0
  Character(200) :: inflow_profile_file = 'inflow_profile.dat'

  ! recycled precursor inflow (inflow_type==2): reads an x-normal slice produced by
  ! probe_output's slice writer (dir='x') from an independent donor simulation with
  ! identical ny/nz, and imposes it as a time-interpolated Dirichlet inflow plane;
  ! inflow_recycle_file is the donor's slice_fileout base name (reads <file>.bin,
  ! <file>_meta.txt, <file>_times.bin)
  Character(200) :: inflow_recycle_file    = ''
  Integer(Int32) :: inflow_recycle_loop    = 1     ! 1: wrap the donor time series when this run's t exceeds its range; 0: clamp to the last frame
  Real   (Int64) :: inflow_recycle_t_offset = 0d0  ! donor_time = t + inflow_recycle_t_offset, before looping/clamping
  ! spanwise (z) shift on every donor-loop wrap (inflow_recycle_loop==1 only), to avoid the
  ! spurious phase-locked streamwise structure spacing (Uconv*donor_duration) that a plain
  ! repeated loop would inject -- same purpose as sem.f90's per-recycle eddy re-randomization
  Integer(Int32) :: inflow_recycle_shift_z = 1
  Integer(Int32) :: inflow_recycle_seed    = 12345
  ! optional SEM inflow mean-temperature profile (y T columns); unset -> mean_profile_T falls back to T_ref
  Character(200) :: inflow_temperature_file = ''
  Integer(Int32) :: sem_profile_format = 0   ! 0=Reynolds-stress (y U V W uu vv ww uv), 1=wind-tunnel TI (see sem.md)
  Real   (Int64) :: sem_Lscale_ratio_y = 0.3d0   ! fallback Loy/Lox when only a length scale's x-component is given (sem_profile_format==1)
  Real   (Int64) :: sem_Lscale_ratio_z = 0.2d0   ! fallback Loz/Lox, same as above
  Integer(Int32) :: sem_n_eddies       = 200      ! <=0: auto-tuned from sem_length_scale and domain geometry
  Real   (Int64) :: sem_length_scale   = 0.01d0   ! <=0: auto-tuned from grid resolution and the inflow profile
  Integer(Int32) :: sem_seed           = 12345

  ! Ensemble SEM (ESEM) config
  Integer(Int32) :: sem_ensemble_samples = 100
  Integer(Int32) :: sem_ensemble_periods = 8
  Character(200) :: sem_sigma_file       = ''
  Integer(Int32) :: sem_eddy_placement   = 0
  Integer(Int32) :: sem_use_esem         = 1
  Integer(Int32) :: sem_divergence_free  = 0

  ! near-wall eddy-size damping (Van Driest form)
  Integer(Int32) :: sem_wall_damping        = 0
  Real   (Int64) :: sem_wall_damping_Aplus  = 25d0

  ! device residency for the scalars sem.f90's per-step (!$acc routine seq) call chain reads directly
  !$acc declare create(inflow_type, inflow_Uconst, sem_n_eddies, sem_length_scale, sem_seed, sem_eddy_placement, sem_use_esem, sem_divergence_free)

  ! finite differences (second derivative)
  Real(Int64) :: ddx1, ddx2, ddx3
  Real(Int64) :: ddy1, ddy2, ddy3
  Real(Int64) :: ddz1, ddz2, ddz3

  ! linear solver
  Integer (Int32) :: nr, nrhs
  Integer (Int32), Dimension(:),   Allocatable :: pivot  
  Complex (Int64), Dimension(:),   Allocatable :: D, DL, DU
  Complex (Int64), Dimension(:,:), Allocatable :: Dyy

  ! pressure gradients
  Real(Int64) :: dPdx, dPdy, dPdz, dPdx_ref, dPdx0

  ! Oscillatory pressure gradient (x and z)
  Real(Int64) :: dPdx_t, dPdz_t
  Real(Int64) :: Ub_x, Ub_z, waveOmega_x, waveOmega_z, phi_wave_x, phi_wave_z

  ! Constant mass-flux (CMFR) forcing: mode 0 = prescribed dPdx (default), 1 = hold bulk velocity at Ub_target
  Integer(Int32) :: flow_forcing_mode
  Real(Int64) :: Ub_target
  Real(Int64) :: dPdx_cmfr   ! diagnostic-only equivalent forcing under CMFR; never fed back into compute_rhs_u

  ! interpolation weights 
  Integer(Int32) :: in1, in2
  Real(Int64), Dimension(:), Allocatable :: weight_y_0, weight_y_1

  ! actual pressure boundary conditions
  Real   (Int64) :: coef_bc_1, coef_bc_2
  Real   (Int64), Dimension(:,:), Allocatable :: bc_1,     bc_2
  Complex(Int64), Dimension(:,:), Allocatable :: bc_1_hat, bc_2_hat
  Logical(Int32) :: pressure_computed

  ! Runge-Kutta 3 coefficients and buffers
  Real(Int64), Dimension(:),     Allocatable :: rk_t
  Real(Int64), Dimension(:,:),   Allocatable :: rk_coef
  Real(Int64), Dimension(:,:,:), Allocatable :: Fu1, Fu2, Fu3
  Real(Int64), Dimension(:,:,:), Allocatable :: Fv1, Fv2, Fv3
  Real(Int64), Dimension(:,:,:), Allocatable :: Fw1, Fw2, Fw3

  !	Eddy Viscosity
  Real   (Int64), Allocatable, Dimension(:,:,:)   :: nu_t

  ! SGS model control: sgs_model 0=DNS,1=Vreman; Cs_vreman is Smagorinsky-equivalent constant (c_V = 2.5*Cs_vreman^2)
  Integer(Int32) :: sgs_model  = 0       ! default: DNS
  Real   (Int64) :: Cs_vreman  = 0.17d0   ! default Vreman constant

  ! Flat-wall equilibrium wall model flag: 0=DNS no-slip (default), 1=smooth log-law EQWM, 2=rough z0 EQWM
  Integer(Int32) :: flat_wall_model_flag = 0

  ! Rough-wall EQWM roughness lengths [m] (used when flat_wall_model_flag==2); momentum (z0)
  ! and thermal (z0h) roughness are independent, per wall
  Real   (Int64) :: z0_ylo  = 0d0, z0_yhi  = 0d0
  Real   (Int64) :: z0h_ylo = 0d0, z0h_yhi = 0d0

  ! Rough-EQWM matching-height grid index (flat_wall_model_flag==2 only): the
  ! u_tau/theta_tau log-law solves sample U/W/T here instead of the literal
  ! first interior cell (j=2), so a fine near-wall grid doesn't put the sample
  ! point inside the roughness sublayer; defaults to 2 (no shift) otherwise.
  ! Computed once in initialization.f90 after the grid is built.
  Integer(Int32) :: j_match_ylo = 2, j_match_yhi = 2
  !$acc declare create(flat_wall_model_flag,z0_ylo,z0_yhi,z0h_ylo,z0h_yhi,j_match_ylo,j_match_yhi)

  ! wall-model Robin BC coefficient arrays
  Real   (Int64), Allocatable, Dimension(:,:,:) :: alpha_x, alpha_y, alpha_z

  ! Thermal Robin-BC coefficient (flat-wall rough EQWM, T_bc_bot/top==2); cell-centred in x,z like alpha_z
  Real   (Int64), Allocatable, Dimension(:,:,:) :: alpha_T

  ! Persisted Obukhov length (nxg,nzg, cell-centred -- matches alpha_T's grid), seeded
  ! neutral and iterated in place each call by solve_most (compute_flat_wall_thermal_eqwm);
  ! carrying it across calls/timesteps keeps the fixed-point iteration's cost low once
  ! the flow is quasi-steady. Only meaningful where T_bc_bot/top==2.
  Real   (Int64), Allocatable, Dimension(:,:) :: L_obukhov_ylo, L_obukhov_yhi
  
  ! Auxillary data variables for roughness
  Real   (Int64) :: Utarget
  Real   (Int64) :: Lx_i, Ly_i, Lz_i
  Real   (Int64) :: alphaGrid

  ! Initial condition type (ic_type 1-6; 6=Taylor-Green Vortex, requires x_bc_type=0 and y_bc_type=0) and noise_percent
  Integer(Int32) :: ic_type      = 1
  Real   (Int64) :: noise_percent = 5.0d0

  ! Vertical grid type (grid_type 1-7)
  Integer(Int32) :: grid_type = 1

  !	Grid sizes for fft
  Real	 (Int64) :: dxmin, dymin, dzmin, Delta

  ! Ghost-cell IBM data structures (phi, Umask_cc, ghost_?_* lists)
  Real   (Int64), Allocatable, Dimension(:,:,:) :: phi
  Real   (Int64), Allocatable, Dimension(:,:,:) :: Umask_cc

  ! Per-cell solid ID from ibm_objid_file (0 = unset/legacy single-object); rounded to Integer when consumed
  Real   (Int64), Allocatable, Dimension(:,:,:) :: ibm_obj_id

  ! ghost-cell lists for U (x-faces), V (y-faces), W (z-faces)
  Integer(Int32) :: n_ghost_u, n_ghost_v, n_ghost_w

  Integer(Int32), Allocatable, Dimension(:,:) :: ghost_u_idx   ! (3, n_ghost_u)
  Real   (Int64), Allocatable, Dimension(:,:) :: ghost_u_wgt   ! (9, n_ghost_u)
  Real   (Int64), Allocatable, Dimension(:,:) :: ghost_u_nrm   ! (3, n_ghost_u)
  Real   (Int64), Allocatable, Dimension(:)   :: ghost_u_yref  ! (   n_ghost_u)
  Integer(Int32), Allocatable, Dimension(:,:) :: ghost_u_ref   ! (3, n_ghost_u)

  Integer(Int32), Allocatable, Dimension(:,:) :: ghost_v_idx
  Real   (Int64), Allocatable, Dimension(:,:) :: ghost_v_wgt
  Real   (Int64), Allocatable, Dimension(:,:) :: ghost_v_nrm
  Real   (Int64), Allocatable, Dimension(:)   :: ghost_v_yref
  Integer(Int32), Allocatable, Dimension(:,:) :: ghost_v_ref

  Integer(Int32), Allocatable, Dimension(:,:) :: ghost_w_idx
  Real   (Int64), Allocatable, Dimension(:,:) :: ghost_w_wgt
  Real   (Int64), Allocatable, Dimension(:,:) :: ghost_w_nrm
  Real   (Int64), Allocatable, Dimension(:)   :: ghost_w_yref
  Integer(Int32), Allocatable, Dimension(:,:) :: ghost_w_ref

  ! Image-point stencil anchor: ghost_?_img(1:3,n) = lower-left (ii,jj,kk) of the enclosing 2x2x2 cube, staggered grid
  Integer(Int32), Allocatable, Dimension(:,:) :: ghost_u_img   ! (3, n_ghost_u)
  Integer(Int32), Allocatable, Dimension(:,:) :: ghost_v_img   ! (3, n_ghost_v)
  Integer(Int32), Allocatable, Dimension(:,:) :: ghost_w_img   ! (3, n_ghost_w)

  ! Per-ghost-point solid ID, resolved from ibm_obj_id at the ghost cell's own index
  ! (mirrors ghost_cc_objid); used to look up per-object roughness (ibm_z0)
  Integer(Int32), Allocatable, Dimension(:)   :: ghost_u_objid   ! (n_ghost_u)
  Integer(Int32), Allocatable, Dimension(:)   :: ghost_v_objid   ! (n_ghost_v)
  Integer(Int32), Allocatable, Dimension(:)   :: ghost_w_objid   ! (n_ghost_w)

  ! Distance from ghost cell G to boundary point B along the wall normal.
  ! Used for surface-integral force Method 2: dA = dV / dGB.
  Real(Int64), Allocatable, Dimension(:) :: ghost_u_dGB   ! (n_ghost_u)
  Real(Int64), Allocatable, Dimension(:) :: ghost_v_dGB   ! (n_ghost_v)
  Real(Int64), Allocatable, Dimension(:) :: ghost_w_dGB   ! (n_ghost_w)

  ! Physical coordinates of boundary point B = G + dGB * nrm.
  ! Used for pressure interpolation in force Method 2.
  Real(Int64), Allocatable, Dimension(:,:) :: ghost_u_xB  ! (3, n_ghost_u)
  Real(Int64), Allocatable, Dimension(:,:) :: ghost_v_xB  ! (3, n_ghost_v)
  Real(Int64), Allocatable, Dimension(:,:) :: ghost_w_xB  ! (3, n_ghost_w)

  ! Cell-centre trilinear stencil at image point I = G+2*dGB*nrm, precomputed for pressure interpolation (force Method 2)
  Integer(Int32), Allocatable, Dimension(:,:) :: ghost_u_img_cc  ! (3, n_ghost_u)
  Integer(Int32), Allocatable, Dimension(:,:) :: ghost_v_img_cc  ! (3, n_ghost_v)
  Integer(Int32), Allocatable, Dimension(:,:) :: ghost_w_img_cc  ! (3, n_ghost_w)
  Real   (Int64), Allocatable, Dimension(:,:) :: ghost_u_wgt_cc  ! (8, n_ghost_u)
  Real   (Int64), Allocatable, Dimension(:,:) :: ghost_v_wgt_cc  ! (8, n_ghost_v)
  Real   (Int64), Allocatable, Dimension(:,:) :: ghost_w_wgt_cc  ! (8, n_ghost_w)

  ! Cell-centre ghost list for rigorous pressure-force integration (Method 2): one entry per solid/fluid interface cell, no overcounting
  Integer(Int32) :: n_ghost_cc = 0
  Integer(Int32), Allocatable, Dimension(:,:) :: ghost_cc_idx     ! (3, n_ghost_cc)
  Real   (Int64), Allocatable, Dimension(:,:) :: ghost_cc_nrm     ! (3, n_ghost_cc)
  Real   (Int64), Allocatable, Dimension(:)   :: ghost_cc_dGB     ! (   n_ghost_cc)
  Integer(Int32), Allocatable, Dimension(:,:) :: ghost_cc_img_cc  ! (3, n_ghost_cc)
  Real   (Int64), Allocatable, Dimension(:,:) :: ghost_cc_wgt_cc  ! (8, n_ghost_cc)
  Integer(Int32), Allocatable, Dimension(:)   :: ghost_cc_objid   ! (   n_ghost_cc) per-ghost-point solid ID, resolved from ibm_obj_id

  ! Precomputed staggered-velocity image-point stencils for viscous-traction export in sample_ibm_surface (built in build_ghost_list_cc)
  Integer(Int32), Allocatable, Dimension(:,:) :: ghost_cc_img_u  ! (3, n_ghost_cc)
  Integer(Int32), Allocatable, Dimension(:,:) :: ghost_cc_img_v  ! (3, n_ghost_cc)
  Integer(Int32), Allocatable, Dimension(:,:) :: ghost_cc_img_w  ! (3, n_ghost_cc)
  Real   (Int64), Allocatable, Dimension(:,:) :: ghost_cc_wgt_u  ! (8, n_ghost_cc)
  Real   (Int64), Allocatable, Dimension(:,:) :: ghost_cc_wgt_v  ! (8, n_ghost_cc)
  Real   (Int64), Allocatable, Dimension(:,:) :: ghost_cc_wgt_w  ! (8, n_ghost_cc)

  ! IBM force monitoring: nsampling writes ibm_forces.csv every N steps, ibm_surface_nsampling dumps per-point surface field
  Integer(Int32) :: nsampling             = 0
  Integer(Int32) :: ibm_force_unit        = -1
  Integer(Int32) :: ibm_surface_nsampling = 0

  ! CFL monitoring/adaptive dt: cfl_adaptive 0=fixed,1=adaptive; cfl_target/cfl_safety drive dt within [dt_min,dt_max]; cfl_current is global max
  Integer(Int32) :: cfl_adaptive  = 0
  Real   (Int64) :: cfl_target    = 0.5d0
  Real   (Int64) :: cfl_safety    = 0.9d0
  Real   (Int64) :: dt_min        = 1d-10
  Real   (Int64) :: dt_max        = 1d10
  Real   (Int64) :: cfl_current   = 0d0
  Real   (Int64) :: cfl_conv_last = 0d0   ! convective CFL from last step
  Real   (Int64) :: cfl_visc_last = 0d0   ! viscous CFL from last step

  ! Suspended sediment transport: sediment_flag 0=off,1=on; sed_bc_bot 0=flux,1=equilibrium; C_ic_type 0=uniform,1=Rouse,2=ramp,3=slab
  Integer(Int32) :: sediment_flag = 0
  Integer(Int32) :: sed_bc_bot    = 0
  Integer(Int32) :: C_ic_type     = 0
  Real   (Int64) :: d_s           = 1d-4    ! particle diameter [m]
  Real   (Int64) :: rho_s         = 2650d0  ! particle density [kg/m^3]
  Real   (Int64) :: rho_f         = 1000d0  ! fluid density [kg/m^3]
  Real   (Int64) :: grav          = 9.81d0  ! gravitational acceleration [m/s^2]
  Real   (Int64) :: Sc            = 1d0     ! molecular Schmidt number
  Real   (Int64) :: Sc_t          = 0.7d0   ! turbulent Schmidt number
  Real   (Int64) :: ws            = 0d0     ! settling velocity (computed at init)
  Real   (Int64) :: C_ref         = 0d0     ! reference near-bed concentration
  Real   (Int64) :: C_ic_height   = 0d0     ! slab height for C_ic_type==3 [m]
  !$acc declare create(grav)

  ! Scalar concentration field (cell-centred: nxg x nyg x nzg)
  Real(Int64), Allocatable, Dimension(:,:,:) :: Cscal, Cscal_o
  Real(Int64), Allocatable, Dimension(:,:,:) :: Fcs1, Fcs2, Fcs3  ! RK3 stage RHS

  ! Boussinesq thermal stratification: boussinesq_flag 0=off,1=on; gravity acts along -y
  ! T_bc_bot/top: 0=adiabatic (zero-gradient), 1=isothermal (Dirichlet against T_wall_bot/top),
  ! 2=rough EQWM flux BC (isothermal target T_wall_bot/top, flux set via z0h_ylo/yhi log law;
  ! requires flat_wall_model_flag=2)
  Integer(Int32) :: boussinesq_flag = 0
  Real   (Int64) :: beta_T   = 0d0    ! thermal expansion coefficient [1/K]
  Real   (Int64) :: T_ref    = 0d0    ! reference temperature for buoyancy [K]
  Real   (Int64) :: Pr       = 0.7d0  ! molecular Prandtl number
  Real   (Int64) :: Pr_t     = 0.85d0 ! turbulent Prandtl number
  Integer(Int32) :: T_bc_bot = 0, T_bc_top = 0
  Real   (Int64) :: T_wall_bot = 0d0, T_wall_top = 0d0
  Integer(Int32) :: T_ic_type = 0     ! 0=uniform (T_ref), 1=linear gradient in y
  Real   (Int64) :: T_ic_grad = 0d0   ! [K/m], used when T_ic_type==1
  ! Per-object IBM thermal BC by solid ID (0=adiabatic, 1=isothermal against ibm_T_wall(id); 0 is the default/legacy slot).
  Integer(Int32), Parameter :: max_ibm_objects = 15
  Integer(Int32) :: ibm_T_bc_type(0:max_ibm_objects) = 0
  Real   (Int64) :: ibm_T_wall   (0:max_ibm_objects) = 0d0
  ! Per-object IBM momentum roughness length [m] (0 = smooth Reichardt EQWM, the
  ! default); only consumed when ibm_wall_model_flag=1. Independent thermal
  ! roughness (ibm_z0h) is not yet wired into the IBM thermal ghost-cell BC --
  ! see apply_ghost_cell_ibm_scalar's plain adiabatic/isothermal mirror, which
  ! would need its own reference-cell EQWM infrastructure (like ghost_u_ref/yref)
  ! to support a flux-consistent rough BC; deferred.
  Real   (Int64) :: ibm_z0(0:max_ibm_objects) = 0d0
  !$acc declare create(boussinesq_flag,beta_T,T_ref,Pr,Pr_t,ibm_T_bc_type,ibm_T_wall,ibm_z0)

  ! UAV actuator disk (src/uav_actuator.f90): Phase 1 scope only -- a STATIC
  ! disk centred at (uav_xc,uav_yc,uav_zc) with uniform loading, applying a
  ! purely vertical (y) reaction force to the fluid. Path/time-dependence and
  ! horizontal force components are later-phase extensions (see
  ! docs/UAV_ActuatorDisk_Design.md).
  ! uav_active:        0=off (default), 1=on
  ! uav_hover_thrust:  disk thrust in this solver's KINEMATIC convention,
  !                    i.e. (physical thrust)/(fluid density) [m^4/s^2] --
  !                    matches dPdx's kinematic convention (this solver
  !                    tracks P/rho, not P; there is no explicit rho anywhere)
  ! uav_kernel_ncell:  regularized-delta (Gaussian) kernel support radius, in
  !                    grid cells, used to spread each marker's force
  Integer(Int32) :: uav_active        = 0
  Real   (Int64) :: uav_xc            = 0d0
  Real   (Int64) :: uav_yc            = 0d0
  Real   (Int64) :: uav_zc            = 0d0
  Real   (Int64) :: uav_disk_radius   = 0.15d0
  Integer(Int32) :: uav_n_r           = 15
  Integer(Int32) :: uav_n_theta       = 24
  Real   (Int64) :: uav_hover_thrust  = 0d0
  Integer(Int32) :: uav_kernel_ncell  = 2
  ! Phase 2: path-following disk (still horizontal/untilted -- orientation
  ! tilt for cruise segments is a later phase, see design doc). When
  ! uav_path_active=0 (default) the disk stays at the fixed (uav_xc,uav_yc,
  ! uav_zc) above; when 1, its centre instead follows uav_path_file (rows
  ! "t x y z", monotonically increasing t) via cubic Hermite (Catmull-Rom
  ! tangent) interpolation, clamped to the first/last waypoint outside the
  ! file's time range.
  Integer(Int32) :: uav_path_active   = 0
  Character(200) :: uav_path_file     = ''
  !$acc declare create(T_bc_bot,T_bc_top,T_wall_bot,T_wall_top)

  ! Temperature field (cell-centred: nxg x nyg x nzg)
  Real(Int64), Allocatable, Dimension(:,:,:) :: Tscal, Tscal_o
  Real(Int64), Allocatable, Dimension(:,:,:) :: Ft1, Ft2, Ft3     ! RK3 stage RHS

  ! Reynolds stress budget (RSB) control and output file layout
  Integer(Int32) :: rsb_active  = 0
  Integer(Int32) :: rsb_freq    = 10
  Integer(Int32) :: rsb_nstart  = 0
  Character(200) :: rsb_hom_dir = 'x,z'
  Character(200) :: rsb_fileout = 'rsb'

  ! in-situ SEM inflow TI-profile rescaling
  Integer(Int32) :: ti_rescale_active = 0
  Real   (Int64) :: ti_rescale_x      = 0d0
  Integer(Int32) :: ti_rescale_nstart = 0      ! <0: auto-tuned from advection time and wall-shear timescale
  Integer(Int32) :: ti_rescale_freq   = 1000   ! <=0: auto-tuned from the SEM eddy turnover time
  Real   (Int64) :: ti_rescale_relax  = 0.3d0
  Real   (Int64) :: ti_rescale_clip   = 1.5d0
  Real   (Int64) :: ti_rescale_abs_clip = 2d0
  Real   (Int64) :: ti_rescale_filter_alpha = 0.25d0   ! EMA smoothing of the measured variance across windows
  Real   (Int64) :: ti_rescale_deadband = 0.03d0        ! skip the update where |ratio-1| is below this
  Real   (Int64) :: ti_rescale_relax_min = 0.05d0       ! floor for the Robbins-Monro-decayed gain

  ! mean-profile (U) companion to TI_RESCALE above: closes the loop on prof_U the same way TI_RESCALE closes it on prof_R11/22/33
  Integer(Int32) :: ti_rescale_u_active   = 0
  Real   (Int64) :: ti_rescale_u_relax    = 0.3d0
  Real   (Int64) :: ti_rescale_u_relax_min = 0.05d0
  Real   (Int64) :: ti_rescale_u_clip     = 0.5d0    ! max |correction| per window, as a fraction of Uconv_sem
  Real   (Int64) :: ti_rescale_u_abs_clip = 1.0d0    ! anti-windup: max cumulative |prof_U-prof_U_target|, as a fraction of Uconv_sem
  Real   (Int64) :: ti_rescale_u_deadband = 0.01d0   ! skip the update where |bias|/Uconv_sem is below this

  ! 2-D planar slice probes: config and output file layout
  Integer(Int32), Parameter :: MAX_PROBES = 8

  Integer(Int32) :: n_slices   = 0
  Integer(Int32) :: slice_freq = 100
  Character(4)   :: slice_dir    (MAX_PROBES) = 'z'
  Real   (Int64) :: slice_pos    (MAX_PROBES) = 0d0
  Character(8)   :: slice_comps  (MAX_PROBES) = 'UVW'
  Character(200) :: slice_fileout(MAX_PROBES) = 'slice'

  ! 1-D line probes: config and output file layout
  Integer(Int32) :: n_lines   = 0
  Integer(Int32) :: line_freq = 100
  Character(4)   :: line_dir    (MAX_PROBES) = 'y'
  Real   (Int64) :: line_pos1   (MAX_PROBES) = 0d0
  Real   (Int64) :: line_pos2   (MAX_PROBES) = 0d0
  Real   (Int64) :: line_start  (MAX_PROBES) = 0d0
  Real   (Int64) :: line_end    (MAX_PROBES) = 1d30
  Character(8)   :: line_comps  (MAX_PROBES) = 'UVW'
  Character(200) :: line_fileout(MAX_PROBES) = 'line'

End Module global
