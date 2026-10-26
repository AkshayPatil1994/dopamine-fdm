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

  ! full_snapshot: restartable full-domain file; chunked_snapshot: sub-box of it. snapshot_location 0 =
  ! staggered, 1 = cell-centred. Box bounds are lengths, or (all six negative) 0-based face indices as -(n+1);
  ! chunk_f1/f2 are the resolved 1-based global faces.
  Logical        :: full_snapshot = .True., chunked_snapshot = .False.
  Integer(Int32) :: snapshot_location = 0
  Real   (Int64) :: x_start = 0d0, x_end = 1d30, y_start = 0d0, y_end = 1d30, z_start = 0d0, z_end = 1d30
  Integer(Int32) :: chunk_f1(3) = 1, chunk_f2(3) = 2
  Integer(Int32) :: nsave, nmonitor

  ! domain size
  Real(Int64) :: Lx, Lz, Ly, Lxp, Lzp

  ! steps
  Integer(Int32) :: nsteps
  Integer(Int32) :: nstep_init = 0
  Real   (Int64) :: dt, t
  Real   (Int64) :: dt_step = 0d0   ! step just completed (dt itself is restored to its unsnapped value)
  ! start time overriding nstep_init*dt when >= 0 (restarts under adaptive dt); -1 = use nstep_init*dt
  Real   (Int64) :: t_start = -1d0

  ! sim_end_time/tsave are used only when nsteps<0 / nsave<0 (absolute end time / save interval); tsave_next
  ! is the next due save time
  Real   (Int64) :: sim_end_time = 1d30
  Real   (Int64) :: tsave        = 1d30
  Real   (Int64) :: tsave_next   = 0d0

  ! chunk-box write cadence, same convention as nsave/tsave; chunk_nsave = 0 follows nsave/tsave
  Integer(Int32) :: chunk_nsave      = 0
  Real   (Int64) :: chunk_tsave      = 1d30
  Real   (Int64) :: chunk_tsave_next = 0d0

  ! viscosity
  Real(Int64) :: nu
  !$acc declare create(nu)

  ! convective term: 0 skew-symmetric (default; suppresses aliasing ringing near sharp/IBM gradients), 1
  ! divergence-form central, 2 skew-symmetric plus upwind dissipation (below). Host-only: consumed once per
  ! RHS call to derive a blend weight.
  Integer(Int32) :: advection_scheme = 0
  ! cells next to a wall-modelled wall (flat_wall_model_flag > 0) that use the divergence form instead of the
  ! skew-symmetric blend
  Integer(Int32) :: adv_wall_layers = -1   ! -1: 8 with a flat-wall model and advection_scheme /= 1, else 0
  ! advection_scheme = 2 (adv_hybrid.f90): upwind-biased dissipation of order adv_order (3 or 5), coefficient
  ! adv_beta0, adv_beta_ibm within adv_band cells of a solid cell
  Integer(Int32) :: adv_order = 3, adv_band = 2
  Real   (Int64) :: adv_beta0 = 0.05d0, adv_beta_ibm = 1d0

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

  ! ibm_input_mode 0 = no body, 1 = SDF ghost-cell IBM; ibm_wall_model_flag 0 = no-slip, 1 = log-law EQWM
  Integer(Int32) :: ibm_input_mode      = 0            ! default: no IBM body
  Character(200) :: ibm_sdf_file        = 'SDF_in'     ! cell-centre SDF file
  ! per-solid ID field (GenSDF sdfp_objid.bin); '' = one uniform IBM condition
  Character(200) :: ibm_objid_file      = ''
  ! smooth_ibm: Jacobi smoothing passes of phi before the ghost-cell lists are built (0 = off). Rounds SDF
  ! corners (~smooth_ibm cells) so the ghost normal stays continuous; a jump there leaves rough fields the
  ! non-dissipative advection cannot damp (Gibbs ringing).
  Integer(Int32) :: smooth_ibm          = 0
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

  ! global grid at middle points + ghost cells
  Real(Int64), Allocatable, Dimension(:) :: xg_global, yg_global, zg_global

  ! local grid at middle points + ghost cells
  Real(Int64), Allocatable, Dimension(:) :: xg, yg, zg

  ! local velocities and pressure
  Real(Int64), Allocatable, Dimension(:,:,:) :: U,V,W,P
  Real(Int64), Allocatable, Dimension(:,:,:) :: Uo,Vo,Wo

  ! local auxiliary
  Real(Int64), Allocatable, Dimension(:,:,:) :: term_1, term_2

  ! pressure right-hand side, then the pseudo-pressure
  Real(Int64), Allocatable, Dimension(:,:,:) :: rhs_p

  ! local auxiliary arrays for MPI_sendrev boundary conditions
  Real(Int64), Allocatable, Dimension(:,:,:) :: buffer_ui, buffer_vi, buffer_ci
  Real(Int64), Allocatable, Dimension(:,:)   :: buffer_ue, buffer_ve, buffer_we, buffer_wi, buffer_ce
  Real(Int64), Allocatable, Dimension(:,:)   :: buffer_p

  ! MPI_isend/irecv interior planes; 3rd dim: 1 = +z exchange, 2 = -z exchange, issued concurrently
  Real(Int64), Allocatable, Dimension(:,:,:) :: buffer_us, buffer_ur
  Real(Int64), Allocatable, Dimension(:,:,:) :: buffer_vs, buffer_vr
  Real(Int64), Allocatable, Dimension(:,:,:) :: buffer_ws, buffer_wr
  ! x ghost-plane exchange (update_ghost_interior_planes_x): sized for the largest (nyg,nzg) field, smaller
  ! fields use a leading subslice
  Real(Int64), Allocatable, Dimension(:,:)   :: buffer_bcxs1, buffer_bcxs2, buffer_bcxr1, buffer_bcxr2
  Real(Int64), Allocatable, Dimension(:,:)   :: buffer_ps, buffer_pr
  Real(Int64), Allocatable, Dimension(:,:)   :: buffer_px
  Real(Int64), Allocatable, Dimension(:,:)   :: buffer_pgxs, buffer_pgxr
  
  ! local pencil work arrays for the Poisson pencil-transpose chain (2decomp&fft)
  Real   (Int64), Allocatable, Dimension(:,:,:) :: poisson_y_r   ! y-pencil, real: interfaces with rhs_p
  Real   (Int64), Allocatable, Dimension(:,:,:) :: poisson_x_r   ! x-pencil, real: DCT-IV path only
  Complex(Int64), Allocatable, Dimension(:,:,:) :: poisson_x_c   ! x-pencil, complex
  Complex(Int64), Allocatable, Dimension(:,:,:) :: poisson_y_c   ! y-pencil, complex: also the Zgtsv operand
  Complex(Int64), Allocatable, Dimension(:,:,:) :: poisson_z_c   ! z-pencil, complex

  ! Fourier points and wave numbers
  Integer(C_INTPTR_T) :: nxp_global, nzp_global
  Integer(Int32)      :: nzp   ! local z-slab count for the physical-space rhs_p ghost fill
  Integer(C_INTPTR_T) :: mx_global, mz_global
  ! local mode range in the y-pencil after the transposes (== global when x/z are not split)
  Integer(C_INTPTR_T) :: mx, mz
  Real   (Int64)      :: dx, dz
  Real   (Int64), Dimension(:), Allocatable :: kxx, kyy, kzz

  ! local (non-MPI) FFTW plans for the post-transpose 1-D transforms
  Type(C_PTR) :: plan_fx_fwd, plan_fx_inv   ! complex 1-D FFT in x (periodic path)
  Type(C_PTR) :: plan_fy_fwd, plan_fy_inv   ! complex 1-D FFT in y (periodic-y path only, y_bc_type==0)
  Type(C_PTR) :: plan_fz_fwd, plan_fz_inv   ! complex 1-D FFT in z (both paths)
  Type(C_PTR) :: plan_dct                    ! real 1-D DCT-IV in x (inflow/outflow path)

  ! streamwise (x) pressure/velocity BC selector: 0=periodic, 1=inflow/outflow
  Integer(Int32) :: x_bc_type  = 0

  ! Outlet backflow stabilisation (x_bc_type==1, ibm_pressure>0): static pressure pb = -beta*|u|^2/2*H(-u_n)
  ! on the outlet face of the Poisson solve, so re-entering fluid carries no net energy. 0 = off (pb = 0,
  ! bit-identical)
  Real(Int64) :: outflow_backflow_beta = 0d0
  ! outflow_Uc: plane-mean outlet speed of the last apply_boundary_conditions; Pb: outlet-face pressure of the
  ! last stage; phib: face value in rhs_p units (dt*c*Pb for the predictor, the increment of the solve for the
  ! projection)
  Real(Int64) :: outflow_Uc = 0d0
  Real(Int64), Allocatable, Dimension(:,:) :: Pb, phib

  ! wall-normal (y) BC selector: 0 = periodic, 1 = wall (bc_face_ylo/yhi)
  Integer(Int32) :: y_bc_type  = 1

  ! spanwise (z) BC selector: 0 = periodic, 1 = wall (DNS no-slip only). A 4-wall duct (y and z walls) uses a
  ! coupled per-mode 2D solve (decomp.f90, projection.f90) and needs p_col==1 (checked in
  ! read_input_parameters).
  Integer(Int32) :: z_bc_type  = 0

  ! &INFLOW (x_bc_type==1 only): inflow_type 0 constant, 1 SEM, 2 recycled precursor slice, 3 wave inlet
  Integer(Int32) :: inflow_type        = 0
  Real   (Int64) :: inflow_Uconst      = 0d0
  Character(200) :: inflow_profile_file = 'inflow_profile.dat'

  ! inflow_type 2: x-normal slice from probe_output's slice writer (dir='x') of a donor run with the same
  ! ny/nz, imposed as a time-interpolated Dirichlet plane. inflow_recycle_file is the donor's slice_fileout
  ! base (<file>.bin, _meta.txt, _times.bin)
  Character(200) :: inflow_recycle_file    = ''
  ! 1: wrap the donor series past its end; 0: clamp to the last frame
  Integer(Int32) :: inflow_recycle_loop    = 1
  Real   (Int64) :: inflow_recycle_t_offset = 0d0  ! donor_time = t + offset, before wrap/clamp
  ! z shift at every donor-loop wrap (inflow_recycle_loop==1), so a plain repeat does not imprint a
  ! phase-locked streamwise structure spacing (same purpose as the per-recycle eddy re-randomization in
  ! sem.f90)
  Integer(Int32) :: inflow_recycle_shift_z = 1
  Integer(Int32) :: inflow_recycle_seed    = 12345
  ! optional SEM inflow mean-temperature profile (y T columns); unset -> mean_profile_T falls back to T_ref
  Character(200) :: inflow_temperature_file = ''
  ! 0=Reynolds-stress (y U V W uu vv ww uv), 1=wind-tunnel TI (see sem.md)
  Integer(Int32) :: sem_profile_format = 0
  ! fallback Loy/Lox when only the x length scale is given (sem_profile_format==1)
  Real   (Int64) :: sem_Lscale_ratio_y = 0.3d0
  Real   (Int64) :: sem_Lscale_ratio_z = 0.2d0   ! fallback Loz/Lox, same as above
  Integer(Int32) :: sem_n_eddies       = 200      ! <=0: auto-tuned from sem_length_scale and domain geometry
  Real   (Int64) :: sem_length_scale   = 0.01d0   ! <=0: auto-tuned from grid resolution and inflow profile
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
  !$acc declare create(inflow_type, inflow_Uconst, sem_n_eddies, sem_length_scale, sem_seed, sem_eddy_placement, sem_use_esem, sem_divergence_free, sem_wall_damping)

  ! linear solver
  Integer (Int32) :: nr, nrhs
  Complex (Int64), Dimension(:),   Allocatable :: D, DL, DU
  Complex (Int64), Dimension(:,:), Allocatable :: Dyy
  ! spanwise wall pressure operator (z_bc_type==1 only), global extent, see Dyy
  Complex (Int64), Dimension(:,:), Allocatable :: Dzz

  ! 4-wall duct (y_bc_type==1 .And. z_bc_type==1): eigendecomposition of the z-wall operator decouples the
  ! (y,z) pressure problem into nzm_global y-tridiagonal solves. Dzz is symmetric only on a uniform z grid, so
  ! S = W^(1/2) Dzz W^(-1/2) (W = diag cell widths) is diagonalised; Dzz's eigenvectors are W^(-1/2) Qz, hence
  ! the sqrt_w_z scaling in solve_poisson_equation.
  Real (Int64), Dimension(:),   Allocatable :: lambda_z    ! eigenvalues, ascending, length nzm_global
  Real (Int64), Dimension(:,:), Allocatable :: Qz          ! orthonormal eigenvectors of S, nzm_global x nzm_global
  Real (Int64), Dimension(:),   Allocatable :: sqrt_w_z    ! sqrt(cell width), length nzm_global
  ! Complex copies of Qz, sqrt_w_z and the z_hat scratch, built at init (Zgtsv/Matmul in
  ! solve_poisson_equation work in Complex)
  Complex (Int64), Dimension(:,:), Allocatable :: Qz_c
  Complex (Int64), Dimension(:),   Allocatable :: sqrt_w_z_c
  Complex (Int64), Dimension(:,:), Allocatable :: z_hat_duct

  ! pressure gradients
  Real(Int64) :: dPdx, dPdz

  ! Oscillatory pressure gradient (x and z)
  Real(Int64) :: dPdx_t, dPdz_t
  Real(Int64) :: Ub_x, Ub_z, waveOmega_x, waveOmega_z, phi_wave_x, phi_wave_z

  ! Constant mass-flux (CMFR) forcing: mode 0 = prescribed dPdx (default), 1 = hold bulk velocity at Ub_target
  Integer(Int32) :: flow_forcing_mode
  Real(Int64) :: Ub_target

  ! Rigid-body rotation about the streamwise (x) axis: rotation_active 0=off,1=on.
  ! Adds Coriolis + centrifugal forcing about the duct centerline (y0_rot,z0_rot),
  ! set to the domain centerline (Ly/2,Lz/2) once Ly_i/Lz_i are known.
  Integer(Int32) :: rotation_active = 0
  Real(Int64) :: Omega_x = 0d0
  Real(Int64) :: y0_rot, z0_rot

  ! interpolation weights 
  Integer(Int32) :: in1, in2
  Real(Int64), Dimension(:), Allocatable :: weight_y_0, weight_y_1, weight_z_0, weight_z_1

  ! Runge-Kutta 3 coefficients and buffers
  Real(Int64), Dimension(:),     Allocatable :: rk_t
  Real(Int64), Dimension(:,:),   Allocatable :: rk_coef
  Real(Int64), Dimension(:,:,:), Allocatable :: Fu1, Fu2, Fu3
  Real(Int64), Dimension(:,:,:), Allocatable :: Fv1, Fv2, Fv3
  Real(Int64), Dimension(:,:,:), Allocatable :: Fw1, Fw2, Fw3

  !	Eddy Viscosity
  Real   (Int64), Allocatable, Dimension(:,:,:)   :: nu_t
  ! turbulent scalar diffusivity (compute_scalar_diffusivity); not in snapshots
  Real   (Int64), Allocatable, Dimension(:,:,:)   :: kappa_t

  ! sgs_model 0=DNS,1=Vreman,2=AMD; Cs_vreman is Smagorinsky-equivalent (c_V = 2.5*Cs_vreman^2)
  Integer(Int32) :: sgs_model  = 2       ! default: AMD
  Real   (Int64) :: Cs_vreman  = 0.17d0   ! default Vreman constant
  Real   (Int64) :: C_amd      = 1d0/3d0  ! AMD Poincare constant, 1/3 for second-order central differences
  ! 1: add the Abkar & Moin (2017) buoyancy production to the AMD numerator (boussinesq_flag >= 1)
  Integer(Int32) :: amd_buoyancy = 1

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

  ! flat-wall EQWM kinematic wall shear along x (U faces) and z (W faces) at the bottom (2nd index 1) and top
  ! (2) y wall; replaces the molecular wall-edge flux of compute_rhs_u / compute_rhs_w when
  ! flat_wall_model_flag > 0
  Real   (Int64), Allocatable, Dimension(:,:,:) :: tau_x, tau_z

  ! Temporal filter of the EQWM input (Yang, Park & Moin 2017): wm_filter_mode 0 raw matching-height velocity,
  ! 1 first-order u_f = (1-eps) u_f + eps u, eps = dt/Tf, Tf = wm_filter_Tf_fac * y_match/(kappa u_tau).
  ! wmf_u*/wmf_w*: filtered tangential vector at U faces (*x) and W faces (*z) of both walls; wmf_t*: the
  ! u_tau setting Tf. Seeded on the first step, also after a restart
  ! -1: 1 on a flat-wall model with no z walls and vof_flow = 0, else 0
  Integer(Int32) :: wm_filter_mode = -1
  Real   (Int64) :: wm_filter_Tf_fac = 1d0
  Logical        :: wmf_seeded = .False.
  Real   (Int64), Allocatable, Dimension(:,:,:) :: wmf_ux, wmf_wx, wmf_tx, wmf_uz, wmf_wz, wmf_tz

  ! IBM method: 0 ghost-cell (image-point mirror), 1 staircase (zero velocity in the solid, masked projection,
  ! wall-model stress as flux through edges next to closed faces, see ibm_stress); ibm_stress_on: viscous
  ! terms read the flux overrides of ibm_stress
  Integer(Int32) :: ibm_method = 0
  ! ibm_pressure (whole-box fast Poisson projection used with every immersed body): 0 pressure-free predictor,
  ! the solve removes the full gradient and body faces are zeroed afterwards, discarding that flux (mass
  ! defect at the body); 1 predictor carries -grad(Pst), Pst = pressure of the last stage, so the solve
  ! supplies only the increment; 2 as 1, and body faces keep the flux of the solve so fluid cells stay
  ! divergence-free (staircase, ibm_method = 1). Pst is updated per stage in time_integration.f90; see
  ! docs/Numerics.md 6.4
  Integer(Int32) :: ibm_pressure = 0
  Real(Int64), Allocatable, Dimension(:,:,:) :: Pst
  ! PCG preconditioner of vof_pressure (masked, variable-density): 0 fast Poisson solver, 1 geometric
  ! multigrid (vof_mg), -1 picks 1 with y walls when the build supports it, else 0
  Integer(Int32) :: pcg_precond = -1
  ! PCG multigrid cycle: 1 V, 2 W, 3 three coarse visits per level, -1 V on GPU, W on CPU
  Integer(Int32) :: mg_cycle = -1
  Integer(Int32) :: mg_smooth_deg = 4  ! degree of the Chebyshev smoother on every level but the coarsest
  Logical        :: ibm_stress_on = .False.
  !$acc declare create(ibm_stress_on)

  ! spanwise (z) wall-model shear stress on U and V at the low/high z wall (index 1/2); z_bc_type==1 and
  ! flat_wall_model_flag==1 only (smooth Reichardt; rough is rejected at input). W is wall-normal and stays
  ! exact no-penetration.
  Real   (Int64), Allocatable, Dimension(:,:,:) :: tau_zu, tau_zv

  ! Wall-edge temperature flux kappa*dT/dy at the low (1) / high (2) y wall from the MOST law
  ! (T_bc_bot/top==2), replacing the diffusive flux in compute_rhs_scalar_core; cell-centred in x,z
  Real   (Int64), Allocatable, Dimension(:,:,:) :: qT_wall

  ! Obukhov length (nxg,nzg, cell-centred), iterated in place by solve_most; carrying it across steps keeps
  ! the fixed-point iteration cheap. Only meaningful where T_bc_bot/top==2.
  Real   (Int64), Allocatable, Dimension(:,:) :: L_obukhov_ylo, L_obukhov_yhi
  
  ! Auxillary data variables for roughness
  Real   (Int64) :: Utarget
  Real   (Int64) :: Lx_i, Ly_i, Lz_i
  Real   (Int64) :: alphaGrid

  ! z grid stretching (z_bc_type==1 only; periodic z is FFT-based and needs uniform spacing): 0 = uniform, >0
  ! = symmetric tanh clustering at both z walls, same convention as alphaGrid of grid_type=2
  Real   (Int64) :: alpha_grid_z = 0d0

  ! ic_type 1-6 (6 = Taylor-Green, needs x_bc_type=0 and y_bc_type=0) and noise_percent
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

  ! Image-point stencil anchor: lower-left (ii,jj,kk) of the enclosing 2x2x2 cube on the staggered grid
  Integer(Int32), Allocatable, Dimension(:,:) :: ghost_u_img   ! (3, n_ghost_u)
  Integer(Int32), Allocatable, Dimension(:,:) :: ghost_v_img   ! (3, n_ghost_v)
  Integer(Int32), Allocatable, Dimension(:,:) :: ghost_w_img   ! (3, n_ghost_w)

  ! Solid ID per ghost point, from ibm_obj_id at the ghost cell (per-object roughness ibm_z0)
  Integer(Int32), Allocatable, Dimension(:)   :: ghost_u_objid   ! (n_ghost_u)
  Integer(Int32), Allocatable, Dimension(:)   :: ghost_v_objid   ! (n_ghost_v)
  Integer(Int32), Allocatable, Dimension(:)   :: ghost_w_objid   ! (n_ghost_w)

  ! Ghost-to-boundary distance along the wall normal (surface-integral force, dA = dV / dGB)
  Real(Int64), Allocatable, Dimension(:) :: ghost_u_dGB   ! (n_ghost_u)
  Real(Int64), Allocatable, Dimension(:) :: ghost_v_dGB   ! (n_ghost_v)
  Real(Int64), Allocatable, Dimension(:) :: ghost_w_dGB   ! (n_ghost_w)

  ! Ghost-to-image distance along the normal: Max(2*dGB, n_image_layers*dymin) keeps I outside G's cell; the
  ! reconstruction uses r = dGB/dGI (r = 0.5 mirror when unclamped)
  Real(Int64), Allocatable, Dimension(:) :: ghost_u_dGI   ! (n_ghost_u)
  Real(Int64), Allocatable, Dimension(:) :: ghost_v_dGI   ! (n_ghost_v)
  Real(Int64), Allocatable, Dimension(:) :: ghost_w_dGI   ! (n_ghost_w)

  ! Cell-centre ghost list for the pressure-force integration: one entry per solid/fluid interface cell
  Integer(Int32) :: n_ghost_cc = 0
  Integer(Int32), Allocatable, Dimension(:,:) :: ghost_cc_idx     ! (3, n_ghost_cc)
  Real   (Int64), Allocatable, Dimension(:,:) :: ghost_cc_nrm     ! (3, n_ghost_cc)
  Real   (Int64), Allocatable, Dimension(:)   :: ghost_cc_dGB     ! (   n_ghost_cc)
  Real   (Int64), Allocatable, Dimension(:)   :: ghost_cc_dGI     ! (n_ghost_cc) image distance, see ghost_u_dGI
  Integer(Int32), Allocatable, Dimension(:,:) :: ghost_cc_img_cc  ! (3, n_ghost_cc)
  Real   (Int64), Allocatable, Dimension(:,:) :: ghost_cc_wgt_cc  ! (8, n_ghost_cc)
  Integer(Int32), Allocatable, Dimension(:)   :: ghost_cc_objid   ! (n_ghost_cc) solid ID from ibm_obj_id

  ! Staggered-velocity image-point stencils for the viscous traction in sample_ibm_surface (built in
  ! build_ghost_list_cc)
  Integer(Int32), Allocatable, Dimension(:,:) :: ghost_cc_img_u  ! (3, n_ghost_cc)
  Integer(Int32), Allocatable, Dimension(:,:) :: ghost_cc_img_v  ! (3, n_ghost_cc)
  Integer(Int32), Allocatable, Dimension(:,:) :: ghost_cc_img_w  ! (3, n_ghost_cc)
  Real   (Int64), Allocatable, Dimension(:,:) :: ghost_cc_wgt_u  ! (8, n_ghost_cc)
  Real   (Int64), Allocatable, Dimension(:,:) :: ghost_cc_wgt_v  ! (8, n_ghost_cc)
  Real   (Int64), Allocatable, Dimension(:,:) :: ghost_cc_wgt_w  ! (8, n_ghost_cc)

  ! IBM forces: nsampling = interval (steps) of ibm_forces.csv; ibm_surface_nsampling = interval of per-point
  ! surface dumps
  Integer(Int32) :: nsampling             = 0
  Integer(Int32) :: ibm_force_unit        = -1
  Integer(Int32) :: ibm_surface_nsampling = 0

  ! cfl_adaptive 0=fixed,1=adaptive: dt follows cfl_target/cfl_safety within [dt_min,dt_max]; cfl_current is
  ! the global max
  Integer(Int32) :: cfl_adaptive  = 0
  Real   (Int64) :: cfl_target    = 0.5d0
  Real   (Int64) :: cfl_safety    = 0.9d0
  Real   (Int64) :: dt_min        = 1d-10
  Real   (Int64) :: dt_max        = 1d10
  Real   (Int64) :: cfl_current   = 0d0
  Real   (Int64) :: cfl_conv_last  = 0d0   ! convective CFL from last step
  Real   (Int64) :: cfl_visc_last  = 0d0   ! viscous CFL from last step
  Real   (Int64) :: cfl_accel_last = 0d0   ! source-term (buoyancy/rotation/UAV) accel-CFL from last step

  ! sediment_flag 0=off,1=on; sed_bc_bot 0=flux,1=equilibrium; C_ic_type 0=uniform,1=Rouse,2=ramp,3=slab
  Integer(Int32) :: sediment_flag = 0
  Integer(Int32) :: sed_bc_bot    = 0
  Integer(Int32) :: C_ic_type     = 0
  Real   (Int64) :: d_s           = 1d-4    ! particle diameter [m]
  Real   (Int64) :: rho_s         = 2650d0  ! particle density [kg/m^3]
  Real   (Int64) :: rho_f         = 1000d0  ! fluid density [kg/m^3]
  ! gravitational acceleration [m/s^2] along -y, the only gravity input (&PHYSICS): sediment settling, Boussinesq,
  ! particles, two-fluid VOF, waves, UAV tilt
  Real   (Int64) :: grav          = 9.81d0
  Real   (Int64) :: Sc            = 1d0     ! molecular Schmidt number
  Real   (Int64) :: Sc_t          = 0.7d0   ! turbulent Schmidt number
  Real   (Int64) :: ws            = 0d0     ! settling velocity (computed at init)
  Real   (Int64) :: C_ref         = 0d0     ! reference near-bed concentration
  Real   (Int64) :: C_ic_height   = 0d0     ! slab height for C_ic_type==3 [m]
  !$acc declare create(grav)

  ! Scalar concentration field (cell-centred: nxg x nyg x nzg)
  Real(Int64), Allocatable, Dimension(:,:,:) :: Cscal, Cscal_o
  Real(Int64), Allocatable, Dimension(:,:,:) :: Fcs1, Fcs2, Fcs3  ! RK3 stage RHS

  ! Boussinesq stratification: boussinesq_flag 0=off,1=on; gravity acts along -y. T_bc_bot/top: 0 adiabatic, 1
  ! isothermal (T_wall_bot/top), 2 rough EQWM flux BC (isothermal target, flux from the z0h_ylo/yhi log law;
  ! needs flat_wall_model_flag=2)
  Integer(Int32) :: boussinesq_flag = 0
  Real   (Int64) :: beta_T   = 0d0    ! thermal expansion coefficient [1/K]
  Real   (Int64) :: T_ref    = 0d0    ! reference temperature for buoyancy [K]
  Real   (Int64) :: Pr       = 0.7d0  ! molecular Prandtl number
  Real   (Int64) :: Pr_t     = 0.85d0 ! turbulent Prandtl number
  Integer(Int32) :: T_bc_bot = 0, T_bc_top = 0
  Real   (Int64) :: T_wall_bot = 0d0, T_wall_top = 0d0
  Integer(Int32) :: T_ic_type = 0     ! 0=uniform (T_ref), 1=linear gradient in y
  Real   (Int64) :: T_ic_grad = 0d0   ! [K/m], used when T_ic_type==1
  ! Per-object IBM thermal BC by solid ID (0 adiabatic, 1 isothermal at ibm_T_wall(id)); slot 0 is the default
  Integer(Int32), Parameter :: max_ibm_objects = 15
  Integer(Int32) :: ibm_T_bc_type(0:max_ibm_objects) = 0
  Real   (Int64) :: ibm_T_wall   (0:max_ibm_objects) = 0d0
  ! Per-object IBM momentum roughness length [m] (0 = smooth Reichardt EQWM), ibm_wall_model_flag=1 only. No
  ! thermal roughness: apply_ghost_cell_ibm_scalar uses a plain adiabatic/isothermal mirror
  Real   (Int64) :: ibm_z0(0:max_ibm_objects) = 0d0
  !$acc declare create(boussinesq_flag,beta_T,T_ref,Pr,Pr_t,ibm_T_bc_type,ibm_T_wall,ibm_z0)

  ! UAV actuator disk (src/uav_actuator.f90): uniform loading, vertical (y) reaction force, static or
  ! path-following, fixed or scheduled thrust.
  ! uav_active:       0=off (default), 1=on
  ! uav_hover_thrust: kinematic thrust (thrust/density) [m^4/s^2], as for dPdx (the solver tracks P/rho)
  ! uav_kernel_ncell: support radius of the Gaussian delta kernel spreading each marker's force, in grid cells
  Integer(Int32) :: uav_active        = 0
  Real   (Int64) :: uav_xc            = 0d0
  Real   (Int64) :: uav_yc            = 0d0
  Real   (Int64) :: uav_zc            = 0d0
  Real   (Int64) :: uav_disk_radius   = 0.15d0
  Integer(Int32) :: uav_n_r           = 15
  Integer(Int32) :: uav_n_theta       = 24
  Real   (Int64) :: uav_hover_thrust  = 0d0
  Integer(Int32) :: uav_kernel_ncell  = 2
  ! Path-following disk (untilted): with uav_path_active=1 the centre follows uav_path_file (rows 't x y z')
  ! by cubic Hermite (Catmull-Rom) interpolation, clamped at the end waypoints; otherwise it stays at
  ! (uav_xc,uav_yc,uav_zc).
  Integer(Int32) :: uav_path_active   = 0
  Character(200) :: uav_path_file     = ''
  ! uav_thrust_active=1: thrust follows uav_thrust_file (rows 't T', interpolated like the path) instead of
  ! uav_hover_thrust.
  Integer(Int32) :: uav_thrust_active = 0
  Character(200) :: uav_thrust_file   = ''
  ! uav_tilt_active: 1 = disk normal from the path acceleration (differentially-flat point mass, n ~ (ax,
  ! grav+ay, az)), low-pass filtered with uav_tilt_tau against Catmull-Rom knot jumps; needs
  ! uav_path_active=1. 0 = horizontal uav_load_profile: 0=uniform, 1=parabolic tip taper (marker share
  ! 1-(r/R)^2, renormalised); uav_tilt_tau: filter time constant [s]; uav_swirl_frac: tangential in-plane
  ! force per marker as a fraction of its thrust share (rotor torque reaction; sense arbitrary)
  Integer(Int32) :: uav_load_profile  = 0
  Integer(Int32) :: uav_tilt_active   = 0
  Real   (Int64) :: uav_tilt_tau      = 0.2d0
  Real   (Int64) :: uav_swirl_frac    = 0d0
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

  ! On-the-fly time-averaged statistics (time_stats.f90): written next to every field snapshot
  Integer(Int32) :: tavg_active = 0
  Integer(Int32) :: tavg_nstart = 0     ! first step (counting nstep_init) that is accumulated
  Integer(Int32) :: tavg_reset  = 0     ! 1: sums restart after every write (independent windows)
  Integer(Int32) :: tavg_T      = 0     ! 1: also T mean, variance and fluxes (boussinesq_flag >= 1)
  Integer(Int32) :: tavg_C      = 0     ! 1: also C mean, variance and fluxes (sediment_flag >= 1)

  ! 2-D planar slice probes: config and output file layout
  Integer(Int32), Parameter :: MAX_PROBES = 30
  Integer(Int32), Parameter :: MAX_WAVEGAUGES = 32

  Integer(Int32) :: n_slices   = 0
  Integer(Int32) :: slice_freq = 100
  Real   (Int64) :: slice_interval = 0d0   ! slice_freq<0: write every slice_interval time units
  Real   (Int64) :: slice_tnext    = 0d0
  Character(4)   :: slice_dir    (MAX_PROBES) = 'z'
  Real   (Int64) :: slice_pos    (MAX_PROBES) = 0d0
  Character(8)   :: slice_comps  (MAX_PROBES) = 'UVW'
  Character(200) :: slice_fileout(MAX_PROBES) = 'slice'

  ! Lagrangian particles (src/particles.f90): particles_active 0=off, 1=on. bc_particle_x/y/z: -1=auto
  ! (periodic if the fluid BC is, else exit(x)/reflect(y,z)), 0=periodic, 1=exit, 2=reflect, 3=absorb (no
  ! bounce). particle_reinit_on_exit: 0=none, 1=re-seed an exited particle at the inflow plane.
  Integer(Int32) :: particles_active       = 0
  Integer(Int32) :: n_particles_init       = 0
  ! xmax<0: full domain once Lx is known
  Real   (Int64) :: particle_seed_xmin = 0d0, particle_seed_xmax = -1d0
  Real   (Int64) :: particle_seed_ymin = 0d0, particle_seed_ymax = -1d0
  Real   (Int64) :: particle_seed_zmin = 0d0, particle_seed_zmax = -1d0
  Integer(Int32) :: particle_seed_seed     = 987654
  Integer(Int32) :: bc_particle_x = -1, bc_particle_y = -1, bc_particle_z = -1
  Integer(Int32) :: particle_reinit_on_exit = 0
  Real   (Int64) :: particle_max_age       = 1d30
  Character(200) :: particle_restart_file  = 'particles_restart'
  ! particle_restart_load: 1=read particles from particle_restart_file when restart==1; 0=always seed fresh
  ! (as scalar_restart for C)
  Integer(Int32) :: particle_restart_load  = 1

  ! particle_mode: 0=tracer, 1=inertial (reduced Maxey-Riley: Schiller-Naumann drag and gravity, always on).
  ! No Saffman-Mei lift: its near-wall gating needs the IBM wall distance.
  Integer(Int32) :: particle_mode        = 0
  Real   (Int64) :: particle_diam        = 1d-4    ! d_p [m]
  Real   (Int64) :: particle_rho         = 2650d0  ! rho_p [kg/m^3]
  ! rho_f of the particle force balance, independent of the sediment/Boussinesq rho_f
  Real   (Int64) :: particle_rho_f       = 1000d0
  ! 1: local Eulerian dU/dt added-mass approximation (advance_particles)
  Integer(Int32) :: particle_added_mass  = 0
  Integer(Int32) :: particle_brownian    = 0       ! 1: isotropic Stokes-Einstein Brownian kick
  Real   (Int64) :: particle_temp_abs    = 293d0   ! [K], particle_brownian==1 only

  ! IBM/SDF collision (ibm_input_mode>=1), BC per solid ID: 1=absorb, 2=reflect (default), 3=deposit_resuspend
  ! (reflect above relative speed particle_resuspend_ucrit, else absorb; a proxy for a pickup function).
  ! Restitution e=Min(1,tau_p/particle_ibm_tau_crit) is a heuristic (apply_ibm_collision in particles.f90).
  Integer(Int32) :: particle_ibm_bc(0:max_ibm_objects) = 2
  Real   (Int64) :: particle_ibm_tau_crit    = 1d-3
  Real   (Int64) :: particle_resuspend_ucrit = 1d30   ! 1d30: always deposits until set

  ! boussinesq_flag>=1 only: 0 keeps particle_rho_f; 1 uses the local rho_f*(1-beta_T*(T-T_ref)) in the
  ! buoyancy term only. One-way: particles never feed back into Tscal or the momentum.
  Integer(Int32) :: particle_boussinesq_coupling = 0
  ! Deposition diagnostic (absorb and deposit_resuspend): streamwise deposition rate written to
  ! particle_deposit_file every particle_deposit_freq monitor reports.
  Character(200) :: particle_deposit_file = 'particle_deposit_x.csv'
  Integer(Int32) :: particle_deposit_freq = 10

  ! LES sub-grid dispersion: 0=none, 1=isotropic Langevin with k_sgs/eps_sgs from nu_t by a Deardorff mixing
  ! length (nu_t=C_k*sqrt(k_sgs)*Delta), see compute_sgs_stats in particles.f90
  Integer(Int32) :: sgs_particle_model  = 0
  Real   (Int64) :: particle_langevin_C0 = 2.1d0   ! Kolmogorov constant

  ! ---- Geometric PLIC VOF liquid fraction and two-fluid solver (src/vof_*.f90, docs/Two-Phase-VOF.md) ----
  ! vof_active 1: C is transported; vof_flow 0: C only advected by the flow, 1: two-fluid Navier-Stokes
  ! (density/viscosity from C, conservative momentum transport with the VOF mass fluxes)
  Integer(Int32) :: vof_active        = 0
  Integer(Int32) :: vof_flow          = 0
  ! initial interface: 1 plane (liquid below y = vof_level), 2 sphere (vof_center, vof_radius), 3 gas bubble
  ! in liquid, 4 standing wave (vof_level, vof_wave_amp, vof_wave_lambda), 5 box x < vof_center(1), y <
  ! vof_center(2) (vof_radius > 0: also the periodic mirror image), 6 disk in the x-y plane (vof_center(1:2),
  ! vof_radius)
  Integer(Int32) :: vof_ic_type       = 1
  Real   (Int64) :: vof_level         = 0d0
  Real   (Int64) :: vof_center(3)     = 0d0
  Real   (Int64) :: vof_radius        = 0d0
  Real   (Int64) :: vof_wave_amp      = 0d0
  Real   (Int64) :: vof_wave_lambda   = 0d0
  ! 1: second-order Stokes surface and deep-water velocities of a +x wave
  Integer(Int32) :: vof_wave_stokes   = 0
  ! vof_normal_scheme 1 Youngs, 2 centred-column height function (Youngs fallback); vof_method 1 PLIC (sharp),
  ! 2 THINC (diffuse over 2-3 cells, vof_beta = sharpness)
  Integer(Int32) :: vof_normal_scheme = 1
  Integer(Int32) :: vof_method        = 1
  Real   (Int64) :: vof_beta          = 2d0
  ! fluids: densities, kinematic viscosities (mu = rho nu), gravity along -y, surface tension (0 = off)
  Real   (Int64) :: vof_rho_l         = 1000d0
  Real   (Int64) :: vof_rho_g         = 1d0
  Real   (Int64) :: vof_nu_l          = 1d-6
  Real   (Int64) :: vof_nu_g          = 1.5d-5
  Real   (Int64) :: vof_sigma         = 0d0
  Real   (Int64) :: vof_u0            = 0d0   ! uniform streamwise velocity at start-up (fresh start only)
  ! vof_mom_scheme 0 central, 1 upwind, 2 Koren, 4 QUICK, 5 central 4th order, 6 WENO5-Z; the face value
  ! blends to upwind where the staggered control-volume mass turns over (smoothstep over refill numbers
  ! vof_mom_cm0..cm1)
  Integer(Int32) :: vof_mom_scheme    = 6
  Real   (Int64) :: vof_mom_cm0       = 2d-3
  Real   (Int64) :: vof_mom_cm1       = 1d-2
  ! vof_geo_density 1: face densities from the reconstructed planes in the half cells of each staggered
  ! control volume; 0 = arithmetic mean (hydrostatically wrong at 1000:1)
  Integer(Int32) :: vof_geo_density   = 1
  ! vof_inertia 1: the pressure operator and every force act on the momentum of the staggered control volume
  ! with the transported mass (arithmetic, vp_ra*); the weight and stresses keep their geometric values
  Integer(Int32) :: vof_inertia       = 0
  ! vof_rk_mom 1: pseudo-time SSP-RK3 per geometric sweep with frozen exact mass fluxes (stable to a full-step
  ! Courant number of about 0.3; switched to 0 above density ratio 2e3), 0 forward Euler; vof_rk_nth
  ! pseudo-time segments (0 = 1); vof_nsub splits the advection half-step (0 = sub-step Courant number below
  ! vof_co_sub); vof_freeze_ut 1 projects the transporting velocity once per half-step
  Integer(Int32) :: vof_rk_mom        = 1
  Integer(Int32) :: vof_rk_nth        = 0
  Integer(Int32) :: vof_nsub          = 0
  Real   (Int64) :: vof_co_sub        = 1.2d-1
  Integer(Int32) :: vof_freeze_ut     = 1
  ! pressure: PCG on the force-stage increment to relative residual vof_pcg_tol (at most vof_pcg_iters
  ! iterations, 0 = always vof_pcg_iters); projection after each advection half-step: vof_adv_iters
  ! iterations, relative residual vof_adv_tol
  Integer(Int32) :: vof_pcg_iters     = 30
  Real   (Int64) :: vof_pcg_tol       = 2d-1
  Integer(Int32) :: vof_adv_iters     = 3
  Real   (Int64) :: vof_adv_tol       = 1d-8
  ! transporting-velocity projection (once per half-step) stops at residual vof_tr_floor * (largest
  ! velocity); the CG residual is the divergence left in the liquid volume
  Real   (Int64) :: vof_tr_floor      = 1d-14
  ! > 0: every pressure solve stops once dt*max|div u| is below this value, instead of the relative tests
  ! above; the caps still apply
  Real   (Int64) :: vof_div_tol       = 0d0
  ! passive tests: vof_prescribed 1 LeVeque vortex reversal, 2 Enright deformation (reversing at t =
  ! vof_presc_T/2), 3 steady rigid rotation with period vof_presc_T; vof_tgv 1 Taylor-Green initial velocity
  ! (amplitude vof_u0)
  Integer(Int32) :: vof_prescribed    = 0
  Real   (Int64) :: vof_presc_T       = 8d0
  Integer(Int32) :: vof_tgv           = 0
  ! diagnostics and development switches
  ! 1: vof_dev.dat / vof_prof.dat / vof_front.dat / vof_vmax.dat; with vof_snap_dt also slices
  Integer(Int32) :: vof_debug         = 0
  ! (vof_debug = 1, one rank) dump vof_snap_NNNN.dat every vof_snap_dt
  Real   (Int64) :: vof_snap_dt       = 0d0
  Integer(Int32) :: vof_selftest      = 0
  Integer(Int32) :: vof_hsplit        = 0     ! 1: well-balanced split p = p' - rho g (y - vof_level)
  Real   (Int64) :: vof_smooth_w      = 0d0   ! vof_ic_type 4: smoothed Heaviside of this half-width (cells)
  ! set by vof_init: refreshes Cvof_io from the device before a snapshot
  Procedure(), Pointer :: vof_sync_io_hook => Null()
  ! C without the second ghost layer, for snapshots and restart
  Real(Int64), Allocatable, Dimension(:,:,:) :: Cvof_io
  ! Waves (src/waves.f90): wave_type 0=off, 1=linear, 2=stream function, 3=JONSWAP; the inlet needs
  ! inflow_type=3. wave_current_mode: 1 zero mean volume flux (closed flume), 2 zero mean Eulerian velocity.
  ! Relaxation zones (waves2Foam weights): generation length from the inlet, absorption length at the outlet,
  ! relaxation rate [1/s] at full strength.
  Integer(Int32) :: wave_type         = 0
  Real   (Int64) :: wave_height       = 0d0
  Real   (Int64) :: wave_period       = 1d0
  Real   (Int64) :: wave_phase        = 0d0
  Integer(Int32) :: wave_sf_n         = 16
  Integer(Int32) :: wave_current_mode = 1
  Real   (Int64) :: wave_Hs           = 0d0
  Real   (Int64) :: wave_Tp           = 1d0
  Real   (Int64) :: wave_gamma        = 3.3d0
  Integer(Int32) :: wave_nfreq        = 200
  Integer(Int32) :: wave_seed         = 12345
  Real   (Int64) :: wave_gen_len      = 0d0
  Real   (Int64) :: wave_abs_len      = 0d0
  Real   (Int64) :: wave_relax_rate   = 20d0
  Real   (Int64) :: wave_ramp_time    = 0d0    ! smooth start-up ramp of the wave amplitude [s] (0 = none)
  ! mean current (throughflow) added to the wave in the water at the inlet and to the absorption target:
  ! current_type 0 none, 1 uniform, 2 log law (roughness current_z0), 3 power law (exponent 1/current_n);
  ! current_U is the depth mean; air carries only the return flow; ramped with wave_ramp_time. wave_turb = 1:
  ! recycled-precursor turbulence (inflow_recycle_file) at the inlet, its fluctuation (donor minus its
  ! z-averaged first frame) faded with the water fraction so air carries none
  Integer(Int32) :: wave_turb         = 0
  Integer(Int32) :: current_type      = 0
  Real   (Int64) :: current_U         = 0d0
  Real   (Int64) :: current_z0        = 1d-3
  Real   (Int64) :: current_n         = 7d0
  ! x of surface-elevation gauges written to vof_gauges.dat (<0 unused)
  Real   (Int64) :: wave_gauge_x(MAX_WAVEGAUGES) = -1d0

  ! 1-D line probes: config and output file layout
  Integer(Int32) :: n_lines   = 0
  Integer(Int32) :: line_freq = 100
  Real   (Int64) :: line_interval = 0d0    ! line_freq<0: write every line_interval time units
  Real   (Int64) :: line_tnext    = 0d0
  Character(4)   :: line_dir    (MAX_PROBES) = 'y'
  Real   (Int64) :: line_pos1   (MAX_PROBES) = 0d0
  Real   (Int64) :: line_pos2   (MAX_PROBES) = 0d0
  Real   (Int64) :: line_start  (MAX_PROBES) = 0d0
  Real   (Int64) :: line_end    (MAX_PROBES) = 1d30
  Character(8)   :: line_comps  (MAX_PROBES) = 'UVW'
  Character(200) :: line_fileout(MAX_PROBES) = 'line'

End Module global
