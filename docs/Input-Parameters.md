# Input Parameters

All settings live in the Fortran namelist file `input_parameters` (read from the working directory). Two ready-made copies are in the repository root:

- `input_parameters`: the bare runnable file, no comments. Copy it and edit the numbers.
- `input_parameters_with_comments`: the same file with every option explained inline, including the optional groups (`&UAV`, `&PARTICLES`, `&WAVES`, `&VOF`, ...) as commented-out blocks.

**Rules**

- Each group is `&NAME ... /`; `!` starts a comment. Unknown variable names abort the run at start-up.
- **Required** groups: `&DOMAIN`, `&PHYSICS`, `&NUMERICS`, `&INITIAL_CONDITIONS`, `&IO`. All others are optional; omit a group to use its defaults (feature off).
- A default of **required** means there is none: set it. Lengths are in metres, times in seconds, temperatures in K. There is no explicit density: pressure gradients and thrusts are kinematic (divided by density).
- Quick start: the "Recipes" section at the end lists the few parameters to change for common cases.

Contents: [DOMAIN](#domain) · [PHYSICS](#physics) · [NUMERICS](#numerics) · [BOUNDARY_CONDITIONS](#boundary_conditions) · [INFLOW](#inflow) · [INFLOW_OPT](#inflow_opt) · [IBM](#ibm) · [INITIAL_CONDITIONS](#initial_conditions) · [IO](#io) · [SEDIMENT](#sediment) · [BOUSSINESQ](#boussinesq) · [STATISTICS](#statistics) · [UAV](#uav) · [PARTICLES](#particles) · [VOF](#vof) · [WAVES](#waves) · [Recipes](#recipes) · [RSB output](#rsb-output-files)

## DOMAIN

| Parameter | Default | Meaning / how to use |
|---|---|---|
| `nx, ny, nz` | required | Face points per direction (interior cells = n−1). |
| `Lx, Ly, Lz` | required | Domain lengths [m]. |
| `grid_type` | `1` | Wall-normal (y) grid: 1 uniform; 2 tanh, fine at both walls; 3 tanh, fine at bottom; 4 tanh, fine at top; 5 uniform sublayer `[0,ks]` (`nks` cells, set in `&IBM`) then tanh fine at the interface and the top wall; 6 same, coarse at the top; 7 same, coarse at the interface and fine at the top. |
| `alpha_grid` | `1.0` | Stretching strength for `grid_type` 2–7 (larger = stronger clustering). Ignored for type 1. |
| `alpha_grid_z` | `0.0` | Spanwise stretching: `0` uniform, `>0` symmetric tanh clustering at both z walls. Needs `z_bc_type=1` (aborts with periodic z). |
| `p_row, p_col` | `0, 0` | MPI pencil grid (`p_row` splits x, `p_col` splits z, y always local; product = rank count). `0,0` = automatic: pure z-slab if every rank keeps ≥2 z-cells, otherwise the factor pair closest to the x:z aspect ratio. With `y_bc_type=1` and `z_bc_type=1` (duct) auto forces `p_row=nprocs, p_col=1`; an explicit `p_col>1` there aborts. |

## PHYSICS

| Parameter | Default | Meaning / how to use |
|---|---|---|
| `nu` | required | Kinematic viscosity [m²/s]. With `vof_flow=1` the larger of `vof_nu_l/g` replaces it for the viscous CFL limit. |
| `dPdx, dPdz` | required | Mean streamwise / spanwise pressure gradient (kinematic, m/s²) that drives the flow. `dPdx` is ignored if `flow_forcing_mode=1`; `dPdx` also sets `u_tau` for `ic_type` 1 and 4. |
| `flow_forcing_mode` | `0` | `0` drive with `dPdx`/`dPdz`; `1` constant mass flux: a uniform shift of U each step holds the bulk velocity at `Ub_target` (needs `x_bc_type=0`, no `T_wave_x`; aborts otherwise). `dPdx` then only reports the equivalent forcing. |
| `Ub_target` | `0.0` | Target bulk streamwise velocity [m/s] for `flow_forcing_mode=1`. |
| `Ub_x, Ub_z` | `0.0` | Amplitude of oscillatory (wave) forcing: `dP/dx(t) = dPdx + Ub_x·ω_x·cos(ω_x t + φ_x)` (same in z). |
| `T_wave_x, T_wave_z` | `0.0` | Oscillation period [s]; `0` = steady. ω = 2π/T. |
| `phi_wave_x, phi_wave_z` | `0.0` | Phase offset [rad] (π/2 gives sine forcing; φ_z−φ_x sets the cross-wave lag). |
| `sgs_model` | `0` | `0` DNS, `1` Vreman LES. |
| `Cs_vreman` | `0.17` | Vreman constant (c_V = 2.5·Cs²). |
| `flat_wall_model_flag` | `0` | Wall model on flat y walls (and on z walls if `z_bc_type=1`): `0` no-slip DNS, `1` smooth log-law EQWM, `2` rough-z0 EQWM (y walls only; aborts with z walls). |
| `z0_ylo, z0_yhi` | `0.0` | Momentum roughness length [m] of the bottom / top wall; mode 2 only, must be `>0` on every no-slip wall (aborts otherwise). |
| `z0h_ylo, z0h_yhi` | `0.0` | Thermal roughness length [m]; only for `T_bc_bot/top=2`. |
| `advection_scheme` | `0` | `0` skew-symmetric (energy-conserving, less ringing at sharp IBM gradients; recommended), `1` pure divergence-form central. |
| `rotation_active` | `0` | `1` rigid-body rotation about x (Coriolis + centrifugal in the v/w equations); the axis is fixed at the centreline (Ly/2, Lz/2). |
| `Omega_x` | `0.0` | Rotation rate [rad/s], right-handed about +x (y turning towards z). |

## NUMERICS

| Parameter | Default | Meaning / how to use |
|---|---|---|
| `dt` | required | Time step [s] (initial step if `cfl_adaptive=1`). |
| `nsteps` | required | `>0` run exactly this many steps; `<0` ignore the count and run until `t ≥ sim_end_time`. |
| `nsave` | required | `>0` write a snapshot every `nsave` steps; `<0` every `tsave` time units (dt is shortened to land on multiples of `tsave`). |
| `nmonitor` | required | Print a monitor line every `nmonitor` steps (must be `>0`). |
| `sim_end_time` | `1e30` | End time; used only when `nsteps<0`. |
| `tsave` | `1e30` | Save interval; used only when `nsave<0`. |
| `cfl_adaptive` | `0` | `1` adapts `dt` to `cfl_target`. |
| `cfl_target` | `0.5` | Target CFL (adaptive). |
| `cfl_safety` | `0.9` | Safety factor on the `dt` update (adaptive). |
| `dt_min, dt_max` | `1e-10`, `1e10` | Bounds on `dt` (adaptive). |
| `pcg_precond` | `-1` | Preconditioner of the masked / variable-density PCG (`ibm_method=1`, `vof_flow≥1`): `0` fast Poisson solver, `1` mask-aware multigrid (needs y walls), `-1` automatic (1 with y walls, else 0; 0 for `vof_flow≥1` on GPU). GPU: `ibm_method=1` needs 1; `vof_flow≥1` cannot use 1. |

## BOUNDARY_CONDITIONS

| Parameter | Default | Meaning / how to use |
|---|---|---|
| `bc_face_ylo, bc_face_yhi` | `1, 1` | Bottom / top wall: `1` no-slip, `2` free-slip. Ignored if `y_bc_type=0`. |
| `x_bc_type` | `0` | Streamwise: `0` periodic (FFT pressure solve); `1` inflow/outflow (Dirichlet inflow from `&INFLOW`, convective outflow `dF/dt+Uc dF/dx=0` with `Uc` the outlet mean clipped to `[0,dx/dt]`, DCT-IV pressure solve; the outflow face is divergence-corrected so mass is conserved). |
| `y_bc_type` | `1` | Wall-normal: `0` periodic (needs `grid_type=1`; required for Taylor–Green `ic_type=6`), `1` walls. |
| `z_bc_type` | `0` | Spanwise: `0` periodic; `1` wall (no-slip or smooth EQWM only; no free-slip, no rough EQWM). With `y_bc_type=1` this is a 4-wall duct (coupled 2-D y–z pressure solve, needs `p_col=1`). |

**GPU build** (`ENABLE_GPU=ON`): the Poisson solve supports `x_bc_type=0|1`, multi-GPU, `y_bc_type=0` only with `x_bc_type=0`, and `z_bc_type=1` with `x_bc_type=0` (duct, or a z wall alone with periodic y). A duct with `x_bc_type=1` aborts: use the CPU build.

## INFLOW

Used when `x_bc_type=1`; omitted = `inflow_type=0`.

| Parameter | Default | Meaning / how to use |
|---|---|---|
| `inflow_type` | `0` | `0` uniform flow `U=inflow_Uconst`; `1` synthetic eddy method (SEM) from `inflow_profile_file`; `2` recycled precursor slice; `3` wave inlet of the two-fluid solver (see `&WAVES`). |
| `inflow_Uconst` | `0.0` | Uniform velocity [m/s] (type 0). |
| `inflow_profile_file` | `'inflow_profile.dat'` | Type 1 profile, rows sorted by increasing y (`#` and blank lines skipped, values outside the range clamped): `y U V W uu vv ww uv [uw vw]` [m, m/s, m²/s²]; only y, U, uu, vv, ww, uv are used (uw=vw=0). |
| `inflow_temperature_file` | `''` | Optional `y T` mean-temperature profile (needs `boussinesq_flag=1`); empty = uniform `T_ref`. |
| `sem_profile_format` | `0` | `0` Reynolds-stress file above; `1` wind-tunnel file `z U Iu Iv Iw [length scales]` (length scales from mixing-length theory if absent; ESEM only). |
| `sem_Lscale_ratio_y, sem_Lscale_ratio_z` | `0.3`, `0.2` | Fallback Loy/Lox, Loz/Lox when a format-1 length-scale column lacks them. |
| `sem_n_eddies` | `200` | Eddies in the virtual box: more = richer turbulence, proportionally more cost; not needed for stability. |
| `sem_length_scale` | `0.01` | Eddy radius [m] in homogeneous mode (no `sem_sigma_file`); also the floor of the mixing-length fallback. |
| `sem_seed` | `12345` | RNG seed (identical on all ranks). |
| `sem_ensemble_samples` | `100` | Samples per recycling period for the ESEM normalisation (one-time start-up cost; lower it for very fine inlets). |
| `sem_ensemble_periods` | `8` | Recycling periods spanned by the ensemble window (more = better estimate, linear cost). |
| `sem_sigma_file` | `''` | Optional inhomogeneous length scales, rows `y sigma_ux sigma_uy sigma_uz sigma_vx ... sigma_wz` (10 columns, all `>0`, sorted by y). Empty = one scale. |
| `sem_eddy_placement` | `0` | `0` uniform in y; `1` PDF-weighted towards small-eddy regions (needs `sem_sigma_file`). |
| `sem_use_esem` | `1` | `1` Ensemble SEM (empirical normalisation, recommended); `0` classical Jarrin SEM (comparison only). |
| `sem_divergence_free` | `0` | `1` divergence-free curl construction (Poletto et al. 2013); needs `sem_use_esem=1`. |
| `sem_wall_damping` | `0` | `1` Van Driest shrinking of near-wall eddies (needs `sem_use_esem=1`); use for wall-bounded profiles that show near-wall mean-velocity overshoot. |
| `sem_wall_damping_Aplus` | `25.0` | Van Driest A⁺ (larger = damping reaches further). |
| `inflow_recycle_file` | `''` | Type 2: base name of a donor run's x-normal slice (`<file>.bin`, `_meta.txt`, `_times.bin`). Required for type 2. Donor must have the same `ny, nz` (and Ly, Lz). |
| `inflow_recycle_loop` | `1` | `1` wrap around the donor's time range; `0` clamp to the last snapshot. |
| `inflow_recycle_t_offset` | `0.0` | Donor time = `t` + offset (start partway into the donor record). |
| `inflow_recycle_shift_z` | `1` | At each wrap apply a random circular z-shift (avoids a phase-locked periodic mode). Needs `inflow_recycle_loop=1`. |
| `inflow_recycle_seed` | `12345` | Seed of that shift (restart-safe). |

Producing a donor slice: in the donor run set in `&STATISTICS` `n_slices=1, slice_dir(1)='x', slice_pos(1)=<x>, slice_comps(1)='UVW', slice_fileout(1)='fields/recycle_slice'`. The slice holds cell-centre values, so V and W are imposed slightly off their staggered faces. If `boussinesq_flag` or `sediment_flag` is on, `slice_comps` must also contain `T` / `C` or the run aborts (`T` otherwise falls back to `inflow_temperature_file` / `T_ref`; `C` to `C_ref`). The `sem_*` settings are unused for type 2. Details: [Numerics §11](Numerics.md#11-streamwise-inflowoutflow-bc-x_bc_type--1).

## INFLOW_OPT

Optional, with `inflow_type=1`: tunes the inflow v′² and w′² profiles (Lamberti et al. 2018) so the statistics at a downstream station match the target. Phases: baseline measurement, doubled-v′²/w′² measurement, one secant correction, verification (reverts to the baseline if the correction is worse). Method and rationale: [Numerics §11.3](Numerics.md#113-inflow-reynolds-stress-optimisation-inflow_opt).

| Parameter | Default | Meaning / how to use |
|---|---|---|
| `inflow_opt_active` | `0` | `1` enable. |
| `inflow_opt_x` | `0.0` | Station [m] where statistics are matched (nearest cell centre). |
| `inflow_opt_nstart` | `-1` | Step at which sampling begins; `<0` auto from the advection and wall-shear times. |
| `inflow_opt_window` | `-1` | Steps per measurement; `<=0` auto from the SEM eddy turnover time. |
| `n_bezier` | `8` | Bezier control points over the height (end points fixed; n−2 optimised). |
| `inflow_opt_wall_exclude` | `1.0` | Control points within this factor of the wall taper length are excluded (their variance is suppressed by construction). |
| `inflow_opt_trust` | `0.5` | Caps each applied step to ±this fraction of the target (guards against noisy slopes). |
| `inflow_opt_max_iter` | `1` | `1` is the validated single corrected step; `>1` is an experimental iterative extension (each costs another window). |
| `inflow_opt_relax` | `0.7` | Step scale / iteration for `max_iter>1`. |
| `inflow_opt_tol` | `0.1` | Stop when the worst relative residual falls below this (`max_iter>1`). |

## IBM

Ghost-cell or staircase immersed boundaries from a precomputed signed-distance field (generate it with [GenSDF](Tools.md#gensdf)). Omit the group if there is no body.

| Parameter | Default | Meaning / how to use |
|---|---|---|
| `ibm_input_mode` | `0` | `0` no body; `1` read the cell-centre SDF `ibm_sdf_file` (solid where `phi<0`). |
| `ibm_method` | `0` | `0` ghost-cell (image-point mirror; resolve walls with ≥2 cells); `1` staircase (zero velocity in solid, masked PCG projection keeping closed faces at zero flux, log-law stress as viscous flux with `ibm_wall_model_flag=1`; needs `ibm_input_mode≥1`, GPU needs `pcg_precond=1`; not for `vof_flow=1`, which always uses the staircase stress). |
| `ibm_pcg_tol` | `1e-8` | Relative PCG residual of the masked projection of `ibm_method=1` (the divergence left after the projection is of this order relative to its pre-projection value; `1e-13` reaches round-off at about 1.3× the cost). |
| `ibm_projection` | `0` | Projection of `ibm_method=1`: `0` masked PCG (closed faces at exactly zero flux, divergence-free to `ibm_pcg_tol`); `1` fast whole-box Poisson solve that ignores the body, after which the closed faces are zeroed (as in direct-forcing/ghost-cell codes). `1` is 3–4× cheaper per run on the cases tested but is not divergence-free next to the body (see [Numerics §6.4](Numerics.md#64-staircase-method-ibm_method--1)). Not available with `vof_flow=1`. |
| `ibm_wall_model_flag` | `0` | `0` no-slip; `1` log-law EQWM on IBM surfaces. With `ibm_method=0` it only rewrites ghost velocities and carries almost no stress in near-1-D shear: use `ibm_method=1` (or `vof_flow=1`) for a wall-modelled IBM. |
| `ibm_sdf_file` | `'SDF_in'` | Path of the SDF file. |
| `ibm_objid_file` | `''` | Optional per-solid ID field (GenSDF `sdfp_objid.bin`) enabling per-object settings below; empty = one condition for all solids. |
| `ibm_z0(0:15)` | `0.0` | Per-object momentum roughness [m]; `0` = smooth Reichardt law. `ibm_wall_model_flag=1` only. No thermal-roughness counterpart: the thermal BC is always smooth. |
| `ibm_T_bc_type(0:15)` | `0` | Per-object thermal BC: `0` adiabatic, `1` isothermal (needs `boussinesq_flag=1`). Plain Dirichlet mirror, no wall model. |
| `ibm_T_wall(0:15)` | `0.0` | Wall temperature [K] where `ibm_T_bc_type=1`. |
| `ks, nks` | required for `grid_type` 5–7 | Roughness sublayer height [m] and number of uniform cells in it. |
| `nsampling` | `0` | IBM force output interval (`0` off). |
| `ibm_surface_nsampling` | `0` | Interval of surface-field dumps (pressure, pressure/viscous force) to `ibm_surface/surface.<step>.bin` (`0` off). |
| `smooth_ibm` | `0` | Jacobi smoothing passes on `phi` before the ghost lists are built; rounds sharp SDF corners to avoid Gibbs ringing. |

## INITIAL_CONDITIONS

| Parameter | Default | Meaning / how to use |
|---|---|---|
| `ic_type` | `1` | `1` log-law + noise; `2` linear/tent + noise; `3` zero mean + noise; `4` Reichardt channel profile (u_τ=√(\|dPdx\|Ly/2)) + structured waves; `5` inverse-linear (anti-tent, `Utarget` at walls) + noise, for fast transition; `6` Taylor–Green vortex of amplitude `Utarget` (analytic; needs `x_bc_type=0`, `y_bc_type=0`); `7` deterministic two-mode perturbation on uniform `Utarget` (ε₁=`noise_percent`/100, ε₂=0.04ε₁). |
| `Utarget` | required | Target bulk or centreline velocity [m/s] (vortex amplitude for type 6). |
| `noise_percent` | `5.0` | Noise amplitude as % of `Utarget` (types 1–3, 5; ignored for 4 and 6). |
| `restart` | `0` | `1` hot-start from `filein`. |
| `nstep_init` | `0` | Starting step number for a restart. |
| `t_start` | `-1.0` | Restart time [s]; `-1` uses `nstep_init·dt`. Needed with `cfl_adaptive=1`. |
| `scalar_restart` | `1` | `1` read C from the restart file; `0` re-initialise C even when `restart=1` (file has no scalar block). |

## IO

| Parameter | Default | Meaning / how to use |
|---|---|---|
| `filein` | required for restart | Restart file read when `restart=1`. |
| `fileout` | required | Base name of the output snapshots. |

## SEDIMENT

Suspended-sediment passive scalar with Soulsby (1997) settling. Omit to disable.

| Parameter | Default | Meaning / how to use |
|---|---|---|
| `sediment_flag` | `0` | `1` enable. |
| `d_s` | `1e-4` | Particle diameter [m]. |
| `rho_s, rho_f` | `2650`, `1000` | Particle and fluid density [kg/m³]. |
| `grav` | `0.0` | Gravity [m/s²], shared with `&BOUSSINESQ` and inertial `&PARTICLES` (warning if the groups disagree; `&BOUSSINESQ` wins). Set it for any of them. |
| `Sc, Sc_t` | `1.0`, `0.7` | Molecular / turbulent Schmidt numbers. |
| `sed_bc_bot` | `0` | Bottom BC: `0` zero flux, `1` fixed `C_ref`. |
| `C_ref` | `0.0` | Near-bed reference concentration (also the inflow value for `x_bc_type=1`). |
| `C_ic_type` | `0` | Initial C: `0` uniform `C_ref`, `1` Rouse, `2` linear ramp `C_ref`→0, `3` slab. |
| `C_ic_height` | `0.0` | Slab height [m] (`C_ic_type=3`). |

## BOUSSINESQ

Temperature transport with buoyancy in the y-momentum equation (gravity along −y). Omit to disable.

| Parameter | Default | Meaning / how to use |
|---|---|---|
| `boussinesq_flag` | `0` | `1` enable. |
| `beta_T` | `0.0` | Thermal expansion coefficient [1/K]. |
| `T_ref` | `0.0` | Reference temperature [K] for buoyancy and the uniform IC. |
| `grav` | `0.0` | Gravity [m/s²] (see `&SEDIMENT`). |
| `Pr, Pr_t` | `0.7`, `0.85` | Molecular / turbulent Prandtl numbers. |
| `T_bc_bot, T_bc_top` | `0, 0` | Wall BC: `0` adiabatic, `1` isothermal at `T_wall_*`, `2` rough-EQWM flux with Monin–Obukhov stability coupling (needs `flat_wall_model_flag=2` and `z0h_*>0`; the momentum stress uses the previous step's Obukhov length). |
| `T_wall_bot, T_wall_top` | `0.0` | Wall temperature [K] for BC types 1 and 2. |
| `T_ic_type` | `0` | `0` uniform `T_ref`; `1` linear gradient. |
| `T_ic_grad` | `0.0` | Gradient [K/m] for `T_ic_type=1`. |

Immersed-body and inflow thermal conditions are set in `&IBM` and `&INFLOW`.

## STATISTICS

Reynolds-stress budget and probes. Omit to disable.

| Parameter | Default | Meaning / how to use |
|---|---|---|
| `rsb_active` | `0` | `1` accumulate the Pope §7.4 budget. |
| `rsb_freq` | `10` | Window length in steps; one sample is written (and the accumulators reset) per window. |
| `rsb_nstart` | `0` | Step at which accumulation starts (skip spin-up). |
| `rsb_hom_dir` | `'x,z'` | Homogeneous averaging directions: `''` none, `'x'`, `'z'`, `'x,z'` (channel). |
| `rsb_fileout` | `'rsb'` | Output base path, e.g. `'stats/rsb'`. |
| `n_slices` | `0` | Number of 2-D slice probes (max 8). |
| `slice_freq` | `100` | Output every this many steps (shared). |
| `slice_dir(n)` | `'z'` | Normal direction `'x'`, `'y'` or `'z'`. |
| `slice_pos(n)` | `0.0` | Position along the normal [m]. |
| `slice_comps(n)` | `'UVW'` | Any of `U V W P T C` concatenated (`T` needs Boussinesq, `C` needs sediment; dropped with a warning otherwise). |
| `slice_fileout(n)` | `'slice'` | Base path. Writes `.bin` (appended `(ncomp,n1,n2)` big-endian float64 per output), `_meta.txt` and `_times.bin` (exact time of each snapshot, used by `inflow_type=2`). Cell-centre values; (n1,n2) = (ny,nz) for x-normal, (nx,nz) for y-normal, (nx,ny) for z-normal. |
| `n_lines` | `0` | Number of 1-D line probes (max 8). |
| `line_freq` | `100` | Output interval (shared). |
| `line_dir(n)` | `'y'` | Direction along the line. |
| `line_pos1(n), line_pos2(n)` | `0.0` | Transverse coordinates: x-line (y, z); y-line (x, z); z-line (x, y). |
| `line_start(n), line_end(n)` | `0.0`, `1e30` | Extent along the line (default full). |
| `line_comps(n)` | `'UVW'` | Components, as for slices. |
| `line_fileout(n)` | `'line'` | Base path; `.bin` appends `(ncomp,npts)` float64 per output, plus `_meta.txt`. |

Example: `n_slices=1, slice_freq=50, slice_dir(1)='z', slice_pos(1)=0.5, slice_comps(1)='UVW', slice_fileout(1)='slices/z05'`.

## UAV

Actuator disk whose thrust is spread as a reaction force over `rhs_u/v/w`. Omit to disable. Working cases: `examples/uav_*` ([Examples](Examples.md#uav-actuator-disk-examplesuav_)).

| Parameter | Default | Meaning / how to use |
|---|---|---|
| `uav_active` | `0` | `1` enable. |
| `uav_xc, uav_yc, uav_zc` | `0.0` | Fixed disk centre [m] (when `uav_path_active=0`). |
| `uav_disk_radius` | `0.15` | Disk radius [m]. |
| `uav_n_r, uav_n_theta` | `15`, `24` | Radial bands / azimuthal sectors of force markers. |
| `uav_hover_thrust` | `0.0` | Thrust / fluid density [m⁴/s²] (same kinematic convention as `dPdx`). |
| `uav_kernel_ncell` | `2` | Support radius [cells] of the Gaussian delta kernel. Keep the disk ≥ this·max(dx,dz) from a periodic x/z edge (the kernel does not wrap). |
| `uav_path_active` | `0` | `1` follow `uav_path_file`: rows `t x y z` sorted by strictly increasing t; Catmull-Rom cubic-Hermite interpolation, clamped outside the range (can under/overshoot a few % at a flat-hold→steep-change transition). |
| `uav_thrust_active` | `0` | `1` follow `uav_thrust_file` (rows `t T`, same interpolation) instead of the constant thrust, e.g. takeoff surge or landing flare. |
| `uav_load_profile` | `0` | `0` uniform loading; `1` parabolic tip taper (weight 1−(r/R)², renormalised). |
| `uav_tilt_active` | `0` | `1` tilt the disk normal along `(ax, uav_grav+ay, az)` from the path acceleration (needs `uav_path_active=1`). Direction only: the thrust magnitude is not trimmed. |
| `uav_tilt_tau` | `0.2` | Low-pass time constant [s] of the tilt (smooths the path's knot-to-knot acceleration jumps). |
| `uav_grav` | `9.81` | Gravity used only by the tilt model. |
| `uav_swirl_frac` | `0.0` | In-plane tangential reaction force as a fraction of each marker's thrust (rotor-torque reaction; net torque only, rotation sense arbitrary). |

All ranks read the path and thrust files independently (`#` and blank lines skipped).

## PARTICLES

Lagrangian point particles ([particles.f90](../src/particles.f90)). Omit to disable. Not supported with `vof_flow=1`.

| Parameter | Default | Meaning / how to use |
|---|---|---|
| `particles_active` | `0` | `1` enable. |
| `n_particles_init` | `0` | Particles seeded at start. |
| `particle_seed_xmin, particle_seed_ymin, particle_seed_zmin` | `0.0` | Lower bound of the seeding box [m]. |
| `particle_seed_xmax, particle_seed_ymax, particle_seed_zmax` | `-1` | Upper bound; `<0` = full domain. |
| `particle_seed_seed` | `987654` | RNG seed. |
| `bc_particle_x, bc_particle_y, bc_particle_z` | `-1` | Per-direction BC: `-1` auto (periodic where the fluid BC is periodic, else exit in x / reflect in y,z); `0` periodic, `1` exit, `2` reflect, `3` absorb. |
| `particle_reinit_on_exit` | `0` | `1` replace a particle leaving the outflow by a new one at the inflow plane; `0` population decays. |
| `particle_max_age` | `1e30` | Maximum particle age [s]. |
| `particle_restart_file` | `'particles_restart'` | Particle restart file. |
| `particle_restart_load` | `1` | With `restart=1`: `1` read particles, `0` seed fresh. |
| `particle_mode` | `0` | `0` tracer; `1` inertial (Schiller–Naumann drag + gravity from `grav`). |
| `particle_diam, particle_rho, particle_rho_f` | `1e-4`, `2650`, `1000` | Diameter [m], particle and fluid density [kg/m³] (inertial). |
| `particle_added_mass` | `0` | `1` local added-mass approximation. |
| `particle_brownian` | `0` | `1` Stokes–Einstein Brownian kick. |
| `particle_temp_abs` | `293` | Temperature [K] for Brownian motion. |
| `particle_ibm_bc(0:15)` | `2` | Per-solid-ID IBM collision: `1` absorb, `2` reflect (restitution `min(1, τ_p/particle_ibm_tau_crit)`), `3` reflect if the relative speed exceeds `particle_resuspend_ucrit`, else absorb. |
| `particle_ibm_tau_crit` | `1e-3` | Response-time scale of the restitution. |
| `particle_resuspend_ucrit` | `1e30` | Resuspension speed threshold [m/s] (effectively always deposits until set). |
| `particle_boussinesq_coupling` | `0` | `1` buoyancy with the local Boussinesq density (one-way). |
| `particle_deposit_file`, `particle_deposit_freq` | `'particle_deposit_x.csv'`, `10` | Deposition-rate CSV and its write interval in monitor reports. |
| `sgs_particle_model` | `0` | `1` simplified isotropic Langevin SGS dispersion for LES. |
| `particle_langevin_C0` | `2.1` | Kolmogorov constant of that model. |

## VOF

Geometric PLIC volume-of-fluid and the two-fluid solver. Method, validation, limits: [Two-Phase VOF](Two-Phase-VOF.md). With `vof_active=0` none of it runs.

| Parameter | Default | Meaning / how to use |
|---|---|---|
| `vof_active` | `0` | `1` advance the liquid fraction `C`. |
| `vof_flow` | `0` | `0` passive `C`; `1` two-fluid solver (density/viscosity from `C`, conservative momentum transport, variable-density projection, gravity). |
| `vof_ic_type` | `1` | Initial interface: `1` plane below `y=vof_level`; `2` sphere; `3` gas bubble; `4` standing/progressive wave; `5` dam-break box `x<vof_center(1), y<vof_center(2)` (`vof_radius>0` adds the periodic mirror); `6` disk in x–y; `7` Zalesak slotted disk (slot width `vof_radius/3`). |
| `vof_level` | `0.0` | Interface height. |
| `vof_center(3), vof_radius` | `0.0` | Centre and radius of the initial shape. |
| `vof_wave_amp, vof_wave_lambda` | `0.0` | `ic_type=4`: amplitude and wavelength (periodic length). |
| `vof_wave_stokes` | `0` | `1` adds the 2nd-order Stokes harmonic and deep-water velocity (breaking-wave start, steepness `ka=2π·amp/λ`). |
| `vof_smooth_w` | `0.0` | `ic_type=4`: smoothed Heaviside half-width [cells]. |
| `vof_normal_scheme` | `1` | Interface normal: `1` Youngs, `2` height function. |
| `vof_method` | `1` | `C` fluxes: `1` PLIC (sharp), `2` THINC (diffuse over 2–3 cells). |
| `vof_beta` | `2.0` | THINC steepness. |
| `vof_rho_l, vof_rho_g` | `1000`, `1` | Liquid / gas density. |
| `vof_nu_l, vof_nu_g` | `1e-6`, `1.5e-5` | Liquid / gas kinematic viscosity. |
| `vof_grav` | `9.81` | Gravity magnitude along −y. |
| `vof_sigma` | `0.0` | Surface tension; `>0` adds balanced-force `σκ∇C` and the capillary dt limit. |
| `vof_u0` | `0.0` | Uniform streamwise velocity at a fresh start (amplitude for `vof_tgv`/`vof_shear`). |
| `vof_mom_scheme` | `6` | Momentum face value: `0` central, `1` upwind, `2` Koren, `4` QUICK, `5` central 4th order, `6` WENO5-Z. |
| `vof_mom_cm0, vof_mom_cm1` | `2e-3`, `1e-2` | Refill-Courant range over which the face value blends from high order to upwind. |
| `vof_geo_density` | `1` | `1` face density from the reconstructed planes (hydrostatically exact); `0` arithmetic mean. |
| `vof_rk_mom` | `1` | `1` pseudo-time SSP-RK3 with frozen mass fluxes (auto 0 above density ratio 2000); `0` forward Euler. |
| `vof_rk_nth` | `0` | RK3 pseudo-time segments (`0` = 1). |
| `vof_nsub, vof_co_sub` | `0`, `0.12` | Advection sub-steps; `0` picks them so the sub-step Courant number stays below `vof_co_sub` (0.03 for Euler). |
| `vof_freeze_ut` | `1` | Project the transporting velocity once per half-step. |
| `vof_pcg_iters, vof_pcg_tol` | `30`, `0.2` | Force-stage pressure PCG: iteration cap and relative residual (`0` = always the cap). |
| `vof_adv_iters, vof_adv_tol` | `3`, `1e-8` | Projection after each advection half-step. |
| `vof_div_tol` | `0.0` | `>0`: stop every pressure solve once `dt·max\|div u\|` is below this. |
| `vof_layered_precond` | `0` | `1` row-mean-density PCG preconditioner (CPU, y walls). |
| `vof_cfl_max` | `0.4` | Courant limit of the advection step. |
| `vof_prescribed, vof_presc_T` | `0`, `8.0` | Passive tests (`vof_flow=0`): `1` LeVeque vortex, `2` Enright deformation (reversing at `T/2`, unit box; `vof_diag.dat` column 28 = L1 error), `3` rigid rotation about (0.5,0.5) of period `T` (Zalesak). |
| `vof_tgv` | `0` | `1` Taylor–Green initial velocity of amplitude `vof_u0`. |
| `vof_shear` | `0` | `1` gas above `vof_level` starts at `vof_u0`. |
| `vof_hsplit` | `0` | `1` well-balanced split `p = p' − ρg(y−vof_level)`. |
| `vof_frozen` | `0` | `1` hold density fixed (no interface transport; diagnostic). |
| `vof_selftest` | `0` | `1` run the PCG iteration-count self-test at start-up (development). |
| `vof_debug` | `0` | `1` write `vof_dev/prof/front/vmax.dat`. |
| `vof_snap_dt` | `0.0` | With `vof_debug=1`: interval of mid-plane slices `vof_snap_NNNN.dat` of C, u, v. |

Not supported with `vof_flow=1` (aborts): sediment, Boussinesq, particles, UAV, rotation; host-only (a GPU build aborts in `vof_init`). `dPdx` / `Ub_target` forcing is not applied; only periodic z. The ghost-cell IBM (`ibm_input_mode=1`) is supported.

## WAVES

Numerical wave flume: Dirichlet wave inlet plus relaxation zones. Needs `vof_flow=1`, `x_bc_type=1`, `inflow_type=3` ([details](Two-Phase-VOF.md#17-waves-inlet-and-relaxation-zones)).

| Parameter | Default | Meaning / how to use |
|---|---|---|
| `wave_type` | `0` | `0` off; `1` linear; `2` Rienecker–Fenton stream function; `3` JONSWAP. |
| `wave_height, wave_period, wave_phase` | `0`, `1`, `0` | Types 1–2. |
| `wave_sf_n` | `16` | Fourier order (type 2). |
| `wave_current_mode` | `1` | `1` zero mean volume flux (closed flume); `2` zero mean Eulerian velocity. |
| `wave_Hs, wave_Tp, wave_gamma` | `0`, `1`, `3.3` | JONSWAP significant height, peak period, peak enhancement. |
| `wave_nfreq, wave_seed` | `200`, `12345` | JONSWAP components and phase seed. |
| `wave_gen_len, wave_abs_len` | `0` | Generation (inlet) and absorption (outlet) zone lengths [m]. |
| `wave_relax_rate` | `20` | Relaxation rate [1/s] at full strength. |
| `wave_ramp_time` | `0` | Start-up ramp [s]. |
| `current_type` | `0` | Mean throughflow added in the water (and as absorption target): `0` none, `1` uniform, `2` log law, `3` power law. Needs `wave_type>0`; ramped with `wave_ramp_time`. |
| `current_U` | `0` | Depth-mean current [m/s]. |
| `current_z0` | `1e-3` | Bed roughness of the log-law current [m]. |
| `current_n` | `7` | Power-law exponent denominator, `u ~ (y/d)^(1/n)`. |
| `wave_turb` | `0` | `1` add the turbulence of a recycled precursor slice (`inflow_recycle_file`, same `ny, nz` as `inflow_type=2`): donor fluctuation × water fraction × ramp, V and W from the donor. Run the precursor with the same bed wall model and set `current_U`/`current_z0` to its mean profile. |
| `wave_gauge_x(8)` | `-1` | Surface-elevation gauge positions, written to `vof_gauges.dat` (`<0` unused). |

## Recipes

| Case | Set these |
|---|---|
| DNS periodic channel | `nx,ny,nz,Lx,Ly,Lz`, `grid_type=2`+`alpha_grid`, `nu`, `dPdx` (or `flow_forcing_mode=1`+`Ub_target`), `dt` or `cfl_adaptive=1`, `nsteps`, `nsave`, `nmonitor`, `Utarget`, `fileout` |
| Wall-modelled LES | as above plus `sgs_model=1`, `flat_wall_model_flag=1` (or 2 with `z0_ylo/yhi`), coarser grid |
| Rough wall (IBM) | `&IBM`: `ibm_input_mode=1`, `ibm_sdf_file`; for a wall model add `ibm_wall_model_flag=1`, `ibm_method=1`, `ibm_z0` |
| Restart | `restart=1`, `filein`, `nstep_init` (and `t_start` if adaptive dt) |
| Inflow/outflow | `x_bc_type=1`, `&INFLOW inflow_type`, `Utarget` |
| Sediment / thermal | `sediment_flag=1` + `grav`, or `boussinesq_flag=1` + `beta_T`, `T_ref`, `grav` |
| Two-phase wave flume | `&VOF vof_active=1, vof_flow=1, vof_ic_type=1, vof_level, vof_sigma`; `&WAVES wave_type, wave_height, wave_period, wave_gen_len, wave_abs_len`; `x_bc_type=1`, `inflow_type=3` |
| Profiles and budgets | `&STATISTICS rsb_active=1` / `n_slices` / `n_lines` |

## RSB output files

Written by rank 0 as little-endian float64 without record markers; each window appends one `(nc, out_nx, out_ny, out_nz)` slice, so the file layout is `(nc, out_nx, out_ny, out_nz, nsamples)`. `out_*` is 1 along every homogeneous direction, the full cell count otherwise. Component order of 6-tensors: 11, 22, 33, 12, 13, 23. On `restart=1` the module reads `nsamples` from `<base>.meta` and keeps appending; a fresh run truncates the files.

| File | nc | Content |
|---|---|---|
| `<base>_Umean.bin` | 3 | mean velocity ⟨U⟩, ⟨V⟩, ⟨W⟩ |
| `<base>_Rij.bin` | 6 | Reynolds stress ⟨u_i′u_j′⟩ |
| `<base>_Pij.bin` | 6 | production |
| `<base>_epsRes.bin` | 6 | resolved dissipation 2ν⟨∂u_i′/∂x_k ∂u_j′/∂x_k⟩ |
| `<base>_epsSGS.bin` | 6 | SGS dissipation 2⟨ν_t s_ij′⟩ |
| `<base>_PiStrain.bin` | 6 | pressure-strain |
| `<base>_DTij.bin` | 6 | turbulent diffusion |
| `<base>_Dnuij.bin` | 6 | viscous diffusion ν∇²R_ij |
| `<base>_PhiPij.bin` | 6 | pressure diffusion |
| `<base>_Resid.bin` | 6 | budget residual (closure check) |
| `<base>.meta` | – | text: grid shape, `nsamples`, endian, dtype, component order |

```python
import numpy as np

def load_rsb(base, term, nc, out_nx, out_ny, out_nz):
    raw = np.fromfile(f"{base}_{term}.bin", dtype="<f8")
    n = nc * out_nx * out_ny * out_nz
    return raw.reshape(nc, out_nx, out_ny, out_nz, len(raw) // n, order='F')

Rij = load_rsb("stats/rsb", "Rij", 6, 1, 128, 1)   # rsb_hom_dir='x,z'
# Rij[0, 0, :, 0, -1] is the final-sample R_11 profile
```

`dopamine_post.rsb.RSBStats` (or `dopamine-post rsb plot`) reads these files: see [Tools](Tools.md#dopamine_postrsb). Theory: [Numerics §9](Numerics.md#9-reynolds-stress-budget).
