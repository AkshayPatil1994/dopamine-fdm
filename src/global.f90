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
  Integer(Int32) :: nsteps
  Integer(Int32) :: nstep_init = 0
  Real   (Int64) :: dt, t
  Real   (Int64) :: dt_step = 0d0   ! size of the step just completed (dt itself is restored to its pre-snap value at the end of each step)
  ! explicit restart start time (overrides nstep_init*dt when >=0d0 -- needed
  ! for restarts under adaptive dt, where step count no longer maps to a fixed
  ! dt*nstep_init); default -1d0 means "not given, fall back to nstep_init*dt"
  Real   (Int64) :: t_start = -1d0

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
  ! smooth_ibm: number of 6-point Jacobi smoothing passes applied to phi before the ghost-cell
  ! lists are built (0 = off, default -- exact original sharp SDF). Rounds sharp SDF corners by
  ! ~smooth_ibm grid cells so the ghost-cell reconstruction's local surface normal stays
  ! continuous there; a discontinuous corner normal otherwise injects a rough field that the
  ! (non-dissipative) skew-symmetric advection scheme can't damp, producing Gibbs-like ringing.
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
  ! x-direction ghost-plane exchange (update_ghost_interior_planes_x): sized to the
  ! largest (dim2,dim3) extent across the U/V/W/P fields it's called for (nyg,nzg);
  ! calls for a smaller field just use a leading (1:n2,1:n3) subslice of the same buffer
  Real(Int64), Allocatable, Dimension(:,:)   :: buffer_bcxs1, buffer_bcxs2, buffer_bcxr1, buffer_bcxr2
  Real(Int64), Allocatable, Dimension(:,:)   :: buffer_ps, buffer_pr
  Real(Int64), Allocatable, Dimension(:,:)   :: buffer_px
  Real(Int64), Allocatable, Dimension(:,:)   :: buffer_pgxs, buffer_pgxr
  
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

  ! spanwise (z) pressure/velocity BC selector: 0=periodic (default), 1=wall (DNS no-slip only --
  ! no wall model yet). y_bc_type==1 .And. z_bc_type==1 (4-wall duct) is supported via a coupled
  ! per-eigenmode 2D (y-z) solve (decomp.f90, initialization.f90, projection.f90), but requires
  ! p_col==1 so the full z-extent is local to every rank (enforced at input read time, see
  ! read_input_parameters).
  Integer(Int32) :: z_bc_type  = 0

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
  !$acc declare create(inflow_type, inflow_Uconst, sem_n_eddies, sem_length_scale, sem_seed, sem_eddy_placement, sem_use_esem, sem_divergence_free, sem_wall_damping)

  ! finite differences (second derivative)
  Real(Int64) :: ddx1, ddx2, ddx3
  Real(Int64) :: ddy1, ddy2, ddy3
  Real(Int64) :: ddz1, ddz2, ddz3

  ! linear solver
  Integer (Int32) :: nr, nrhs
  Integer (Int32), Dimension(:),   Allocatable :: pivot  
  Complex (Int64), Dimension(:),   Allocatable :: D, DL, DU
  Complex (Int64), Dimension(:,:), Allocatable :: Dyy
  Complex (Int64), Dimension(:,:), Allocatable :: Dzz   ! spanwise (z) wall pressure operator (z_bc_type==1 only), global extent -- see Dyy

  ! 4-wall duct (y_bc_type==1 .And. z_bc_type==1 only): eigendecomposition of the z-wall operator's
  ! interior tridiagonal part, used to decouple z from the coupled 2D (y,z) pressure Poisson
  ! problem into nzm_global independent 1D y-tridiagonal solves, one per z-eigenmode, each
  ! reusing the same Zgtsv-in-y machinery as the y_bc_type==1/z periodic case (just adding
  ! lambda_z(m) to the diagonal instead of kzz). See solve_poisson_equation.
  !
  ! Dzz itself is symmetric only for a UNIFORM z grid; a stretched z (alpha_grid_z>0) makes it
  ! non-symmetric, so we diagonalise the similarity-transformed S = W^(1/2) Dzz W^(-1/2) instead
  ! (W = diag(cell widths) -- a standard finite-volume trick: W*Dzz is exactly symmetric by
  ! construction, so S is too), which shares Dzz's eigenvalues and has orthonormal eigenvectors
  ! Qz. Dzz's own (non-orthogonal in general) eigenvectors are V = W^(-1/2)*Qz, so the forward/
  ! inverse transforms in solve_poisson_equation pre/post-multiply by sqrt_w_z = sqrt(diag(W)).
  ! Reduces exactly to the plain uniform-grid case when sqrt_w_z is constant (cancels out).
  Real (Int64), Dimension(:),   Allocatable :: lambda_z    ! eigenvalues, ascending, length nzm_global
  Real (Int64), Dimension(:,:), Allocatable :: Qz          ! orthonormal eigenvectors of S, nzm_global x nzm_global
  Real (Int64), Dimension(:),   Allocatable :: sqrt_w_z    ! sqrt(cell width), length nzm_global
  ! Complex copies of Qz/sqrt_w_z (Zgtsv/Matmul in solve_poisson_equation both work in Complex),
  ! and the z_hat scratch buffer -- all built once at init instead of every solve_poisson_equation
  ! call (was Allocate/Deallocate + Dcmplx conversion on every RK substage, i.e. 3x per step)
  Complex (Int64), Dimension(:,:), Allocatable :: Qz_c
  Complex (Int64), Dimension(:),   Allocatable :: sqrt_w_z_c
  Complex (Int64), Dimension(:,:), Allocatable :: z_hat_duct

  ! pressure gradients
  Real(Int64) :: dPdx, dPdy, dPdz, dPdx_ref, dPdx0

  ! Oscillatory pressure gradient (x and z)
  Real(Int64) :: dPdx_t, dPdz_t
  Real(Int64) :: Ub_x, Ub_z, waveOmega_x, waveOmega_z, phi_wave_x, phi_wave_z

  ! Constant mass-flux (CMFR) forcing: mode 0 = prescribed dPdx (default), 1 = hold bulk velocity at Ub_target
  Integer(Int32) :: flow_forcing_mode
  Real(Int64) :: Ub_target
  Real(Int64) :: dPdx_cmfr   ! diagnostic-only equivalent forcing under CMFR; never fed back into compute_rhs_u

  ! Rigid-body rotation about the streamwise (x) axis: rotation_active 0=off,1=on.
  ! Adds Coriolis + centrifugal forcing about the duct centerline (y0_rot,z0_rot),
  ! set to the domain centerline (Ly/2,Lz/2) once Ly_i/Lz_i are known.
  Integer(Int32) :: rotation_active = 0
  Real(Int64) :: Omega_x = 0d0
  Real(Int64) :: y0_rot, z0_rot

  ! interpolation weights 
  Integer(Int32) :: in1, in2
  Real(Int64), Dimension(:), Allocatable :: weight_y_0, weight_y_1, weight_z_0, weight_z_1

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

  ! spanwise (z) wall-model Robin BC coefficients (z_bc_type==1, flat_wall_model_flag==1 only --
  ! smooth Reichardt EQWM; rough (flag==2) is not yet supported for z walls). U,V are tangential
  ! to a z wall (Robin); W is the wall-normal component and stays exactly no-penetration, so
  ! there is no alpha_z_w -- mirrors how alpha_y (V, wall-normal at a y wall) is always 0.
  Real   (Int64), Allocatable, Dimension(:,:,:) :: alpha_z_u, alpha_z_v

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

  ! Spanwise (z) grid stretching (z_bc_type==1 only -- periodic z is FFT-based and needs
  ! uniform spacing): 0 (default) = uniform, matching all existing behaviour; >0 = symmetric
  ! tanh clustering at both z walls, same formula/parameter convention as grid_type=2's alphaGrid.
  Real   (Int64) :: alpha_grid_z = 0d0

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

  ! Distance from ghost cell G to image point I along the wall normal:
  ! Max(2*dGB, n_image_layers*dymin), clamped up from the mirror distance 2*dGB
  ! when the true dGB is small compared to the local grid spacing (small/thin
  ! features on a coarse grid) so I reliably lands outside G's own grid cell.
  ! Ghost reconstruction uses the general two-point form with r = dGB/dGI
  ! (reduces to the textbook r=0.5 mirror when dGI is not clamped).
  Real(Int64), Allocatable, Dimension(:) :: ghost_u_dGI   ! (n_ghost_u)
  Real(Int64), Allocatable, Dimension(:) :: ghost_v_dGI   ! (n_ghost_v)
  Real(Int64), Allocatable, Dimension(:) :: ghost_w_dGI   ! (n_ghost_w)

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
  Real   (Int64), Allocatable, Dimension(:)   :: ghost_cc_dGI     ! (   n_ghost_cc) image distance, see ghost_u_dGI
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
  Real   (Int64) :: cfl_conv_last  = 0d0   ! convective CFL from last step
  Real   (Int64) :: cfl_visc_last  = 0d0   ! viscous CFL from last step
  Real   (Int64) :: cfl_accel_last = 0d0   ! source-term (buoyancy/rotation/UAV) accel-CFL from last step

  ! Suspended sediment transport: sediment_flag 0=off,1=on; sed_bc_bot 0=flux,1=equilibrium; C_ic_type 0=uniform,1=Rouse,2=ramp,3=slab
  Integer(Int32) :: sediment_flag = 0
  Integer(Int32) :: sed_bc_bot    = 0
  Integer(Int32) :: C_ic_type     = 0
  Real   (Int64) :: d_s           = 1d-4    ! particle diameter [m]
  Real   (Int64) :: rho_s         = 2650d0  ! particle density [kg/m^3]
  Real   (Int64) :: rho_f         = 1000d0  ! fluid density [kg/m^3]
  Real   (Int64) :: grav          = 0d0     ! gravitational acceleration [m/s^2]; off unless set via &SEDIMENT or &BOUSSINESQ
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

  ! UAV actuator disk (src/uav_actuator.f90): a disk with uniform loading,
  ! applying a purely vertical (y) reaction force to the fluid; static or
  ! path-following, with a fixed or scheduled thrust. Horizontal force
  ! components are not modelled.
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
  ! Time-varying thrust schedule: when uav_thrust_active=0 (default), the
  ! disk uses the fixed uav_hover_thrust above for its whole run; when 1,
  ! it instead uses uav_thrust_file (rows "t T", same interpolation as the
  ! path) -- e.g. a takeoff surge above hover thrust, a reduced-thrust
  ! controlled descent, and a landing flare, all as a function of time
  ! (equivalently of position, since position is itself a function of time
  ! along uav_path_file).
  Integer(Int32) :: uav_thrust_active = 0
  Character(200) :: uav_thrust_file   = ''
  ! uav_load_profile:  0=uniform disk loading (default), 1=parabolic tip-taper
  !                    (marker share weighted by 1-(r/R)^2, renormalized to
  !                    still sum to 1) -- a drop-in reweighting of the marker
  !                    table, no change to the force-application code path.
  ! uav_tilt_active:   0=disk stays horizontal (default, i.e. identical to
  !                    the untilted model above); 1=disk normal is derived
  !                    automatically each step from the path's own kinematic
  !                    acceleration (differentially-flat point-mass tilt:
  !                    n ~ (ax, grav+ay, az)), low-pass filtered with time
  !                    constant uav_tilt_tau to tame the Catmull-Rom path's
  !                    knot-to-knot acceleration discontinuities. Requires
  !                    uav_path_active=1 to have any effect (a static disk's
  !                    path acceleration is identically zero).
  ! uav_tilt_tau:      low-pass time constant [s] for the tilt filter above.
  ! uav_swirl_frac:    0=no swirl (default); tangential (in-plane) reaction
  !                    force per marker as a fraction of that marker's own
  !                    thrust share, representing rotor torque reaction.
  !                    Rotation sense is an arbitrary modelling choice, not
  !                    derived from any tracked rotor RPM/direction.
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

  ! Bezier-parametrized SEM inflow Reynolds-stress optimization: a single-run (online) realization
  ! of Lamberti et al. 2018 (JWEIA 177:32-44) Sections 5-6.1. Matches the downstream v'^2/w'^2
  ! profiles at a station to the wind-tunnel target by fitting Bezier control points: step0
  ! (baseline) plus step1 (v'^2 AND w'^2 doubled together at the inflow, their Section 6.1 combined
  ! perturbation) give a per-control-point scalar secant slope for each decision variable, then ONE
  ! corrected profile is applied and verified. v'^2 and w'^2 are corrected independently (a real
  ! run showed the paper's coupled 2-variable weighted least-squares fit, Eq. 5, can be dangerously
  ! ill-conditioned: the two decision variables tend to move all three downstream stats in the same
  ! direction, so their Jacobian columns can be nearly collinear at some heights, and even Tikhonov
  ! regularization on the 2x2 solve wasn't enough to tame the resulting instability under real
  ! turbulent measurement noise -- a decoupled scalar secant per component needs no matrix
  ! inversion at all, so there is no collinearity to be unstable about). By default
  ! (inflow_opt_max_iter=1) the algorithm stops after that one corrected step: the paper itself
  ! only ever validates a single corrected step (Section 6.1) and explicitly lists further
  ! iteration and an automatic stopping criterion as unsolved future work (Section 7) -- so this is
  ! not an arbitrary simplification, it's matching what was actually shown to work. Continuing past
  ! that single step (inflow_opt_max_iter>1) is this codebase's own, unvalidated-by-the-paper
  ! extension: it secant-refines the slopes and keeps correcting with a Robbins-Monro-style step
  ! size (inflow_opt_relax/iter, decaying so noisy sequential online measurements average out
  ! instead of being chased), guarded by the same best-iterate/stall safety net either way. u'^2,
  ! shear stress and the mean profile are not decision variables (per the paper) and are left
  ! untouched. See docs/design-notes/sem.md.
  Integer(Int32) :: inflow_opt_active = 0
  Real   (Int64) :: inflow_opt_x      = 0d0
  Integer(Int32) :: inflow_opt_nstart = -1     ! <0: auto-tuned from advection time and wall-shear timescale
  Integer(Int32) :: inflow_opt_window = -1     ! <=0: auto-tuned from the SEM eddy turnover time; steps averaged per measurement phase
  Integer(Int32) :: n_bezier          = 8      ! Bezier control points spanning prof_y; endpoints fixed to the target
  ! A control point within wall_exclude_factor*wall_Ltaper_{lo,hi} of an active no-slip wall sits
  ! inside sem_fluctuation's own taper zone (see sem.f90), which deterministically suppresses the
  ! injected Reynolds stress toward zero there regardless of the Bezier target -- no correction can
  ! close that gap since it isn't a response-model error, it's the no-slip enforcement working as
  ! designed. Such control points are excluded from both the correction and the residual check.
  ! This has no counterpart in the paper (its offline runs were reviewed by eye), but is needed
  ! for an unattended online run to avoid chasing a structurally unfixable control point.
  Real   (Int64) :: inflow_opt_wall_exclude = 1d0
  ! A secant slope estimated from ONE fast online measurement window can be small/noisy at some
  ! control points (near-zero measured response between step0 and step1, unlike the paper's fully
  ! time-converged, independently-restarted offline perturbation runs), which blows up the Newton
  ! step dx=-(measured-target)/slope regardless of whether v'^2/w'^2 are solved jointly or (as
  ! here) independently -- a real run showed corrections up to ~30x target from this alone. This
  ! is a numerical-robustness safeguard, not a paper-fidelity choice: it caps the applied step to
  ! +-inflow_opt_trust of the target value at each control point.
  Real   (Int64) :: inflow_opt_trust    = 0.5d0
  Integer(Int32) :: inflow_opt_max_iter = 1    ! 1 (default): the paper's validated single corrected step. >1: this codebase's own experimental extension (see above)
  Real   (Int64) :: inflow_opt_relax    = 0.7d0 ! base step-size scale for the experimental iter>1 extension only (effective scale = inflow_opt_relax/iter); unused when inflow_opt_max_iter=1
  Real   (Int64) :: inflow_opt_tol      = 0.1d0 ! experimental extension only: stop iterating once the worst-case relative residual (|measured-target|/target, over u'^2/v'^2/w'^2 and all interior control points) drops below this

  ! 2-D planar slice probes: config and output file layout
  Integer(Int32), Parameter :: MAX_PROBES = 8

  Integer(Int32) :: n_slices   = 0
  Integer(Int32) :: slice_freq = 100
  Character(4)   :: slice_dir    (MAX_PROBES) = 'z'
  Real   (Int64) :: slice_pos    (MAX_PROBES) = 0d0
  Character(8)   :: slice_comps  (MAX_PROBES) = 'UVW'
  Character(200) :: slice_fileout(MAX_PROBES) = 'slice'

  ! Lagrangian point-particle tracking (src/particles.f90): particles_active 0=off,1=on.
  ! Per-direction particle BC codes (bc_particle_x/y/z): -1=auto (periodic if the matching
  ! fluid x_bc_type/y_bc_type/z_bc_type==0, else exit(x)/reflect(y,z)); explicit override
  ! 0=periodic,1=exit(inflow/outflow),2=reflect(wall),3=absorb(wall, no bounce).
  ! particle_reinit_on_exit: 0=none (population decays as particles exit), 1=inflow
  ! (replace an outflow exit with a new particle at the inflow plane).
  Integer(Int32) :: particles_active       = 0
  Integer(Int32) :: n_particles_init       = 0
  Real   (Int64) :: particle_seed_xmin = 0d0, particle_seed_xmax = -1d0   ! <0 (xmax): resolved to the full domain once Lx is known
  Real   (Int64) :: particle_seed_ymin = 0d0, particle_seed_ymax = -1d0
  Real   (Int64) :: particle_seed_zmin = 0d0, particle_seed_zmax = -1d0
  Integer(Int32) :: particle_seed_seed     = 987654
  Integer(Int32) :: bc_particle_x = -1, bc_particle_y = -1, bc_particle_z = -1
  Integer(Int32) :: particle_reinit_on_exit = 0
  Real   (Int64) :: particle_max_age       = 1d30
  Character(200) :: particle_restart_file  = 'particles_restart'
  ! particle_restart_load: 1=read particles from particle_restart_file when restart==1
  ! (default, mirrors the main fluid restart); 0=always seed fresh particles even when
  ! restart==1 (e.g. hot-starting the flow field but starting a new particle release) --
  ! same idea as scalar_restart for C.
  Integer(Int32) :: particle_restart_load  = 1

  ! Phase 2: inertial force model. particle_mode: 0=tracer (dx/dt=u_fluid, Phase 1
  ! behaviour, default), 1=inertial (independent particle velocity, Maxey-Riley-reduced
  ! ODE: nonlinear (Schiller-Naumann) drag + gravity, both always on in this mode; a
  ! small-Stokes-number inertial particle already reproduces settling/tracer-like
  ! behaviour on its own, so there is no separate "settling_tracer" mode). Saffman-Mei
  ! lift is deferred to Phase 3 (its usual near-wall gating needs the IBM SDF wall
  ! distance that Phase 3 introduces).
  Integer(Int32) :: particle_mode        = 0
  Real   (Int64) :: particle_diam        = 1d-4    ! d_p [m]
  Real   (Int64) :: particle_rho         = 2650d0  ! rho_p [kg/m^3]
  Real   (Int64) :: particle_rho_f       = 1000d0  ! rho_f [kg/m^3] for the particle force balance -- independent of sediment/Boussinesq rho_f by default
  Integer(Int32) :: particle_added_mass  = 0       ! 0=off (default), 1=on: local Eulerian dU/dt added-mass approximation (see advance_particles)
  Integer(Int32) :: particle_brownian    = 0       ! 0=off (default), 1=on: isotropic Stokes-Einstein Brownian kick
  Real   (Int64) :: particle_temp_abs    = 293d0   ! [K], particle_brownian==1 only

  ! Phase 3: IBM/SDF collision (ibm_input_mode>=1 only). Per-object BC by solid ID (same
  ! 0..max_ibm_objects convention as ibm_T_bc_type/ibm_z0): 1=absorb (deposit, removed),
  ! 2=reflect (default), 3=deposit_resuspend (reflect if the local relative speed exceeds
  ! particle_resuspend_ucrit, else absorb -- a simplified proxy for a full Shields/van Rijn
  ! pickup function). Reflection uses a simple Stokes-number-dependent restitution heuristic
  ! e=Min(1,tau_p/particle_ibm_tau_crit) (documented simplification, not a validated closed
  ! form) -- see particles.f90's apply_ibm_collision.
  Integer(Int32) :: particle_ibm_bc(0:max_ibm_objects) = 2
  Real   (Int64) :: particle_ibm_tau_crit    = 1d-3
  Real   (Int64) :: particle_resuspend_ucrit = 1d30   ! effectively "always deposits" until set

  ! Phase 4: Boussinesq coupling (boussinesq_flag>=1 only). 0=off (default): particle_rho_f
  ! stays the fixed value the user set. 1=on: the buoyancy term's fluid density is instead
  ! the local Boussinesq value rho_f*(1-beta_T*(T-T_ref)) interpolated at the particle -- only
  ! the buoyancy term uses this local value (drag/Re_p keep the constant particle_rho_f, a
  ! scoped simplification). Ships one-way (particles never feed back into Tscal/momentum).
  Integer(Int32) :: particle_boussinesq_coupling = 0
  ! Deposition diagnostic (any ibm_input_mode>=1 collision that removes a particle, both
  ! Phase 3's absorb and deposit_resuspend outcomes): streamwise deposition-rate accumulator,
  ! written to particle_deposit_file every particle_deposit_freq monitor reports. Structured
  ! so a later two-way concentration/deposition feedback into the flow is a small increment.
  Character(200) :: particle_deposit_file = 'particle_deposit_x.csv'
  Integer(Int32) :: particle_deposit_freq = 10

  ! Phase 5: LES sub-grid dispersion. 0=none (default, correct for DNS resolution); 1=langevin:
  ! a simplified (isotropic) Thomson/Weil-Sullivan-Moeng well-mixed Langevin model, diagnosing
  ! k_sgs/eps_sgs from the existing eddy-viscosity SGS model's nu_t via a Deardorff-style
  ! mixing-length closure (nu_t=C_k*sqrt(k_sgs)*Delta) rather than a transported k_sgs
  ! equation -- a documented simplification, see particles.f90's compute_sgs_stats.
  Integer(Int32) :: sgs_particle_model  = 0
  Real   (Int64) :: particle_langevin_C0 = 2.1d0   ! Kolmogorov constant

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
