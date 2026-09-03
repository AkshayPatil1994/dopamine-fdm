# UAV Actuator-Disk Design Notes (branch `uav`)

Status: Phase 1 (static disk, uniform loading, vertical-only force)
implemented and ground-effect-checked -- `src/uav_actuator.f90`, `&UAV`
namelist (see `docs/Input-Parameters.md`), a hover smoke test
(`examples/uav_hover_disk`), and an IGE/OGE ground-effect comparison
(`examples/uav_ground_effect`) showing the correct qualitative trend
(reduced induced velocity near the ground at fixed thrust). Phase 2
(path-following disk) is also implemented and smoke-tested
(`examples/uav_path_takeoff`): a vertical climb from h/R=2 to h/R=10
stays stable (divergence at machine precision throughout) and the
downwash signature in a fixed centreline probe visibly translates from
near the ground up to the new altitude as the disk climbs -- with a
noticeable lag behind the disk's instantaneous position, which is the
flow's own inertia/recirculation, not a bug in the force placement
(the force itself is applied at the disk's exact current position every
RK sub-stage; what's lagging is the fluid's response). The MPI rank-
boundary-crossing case is also checked (`examples/uav_path_cross_rank`,
2 ranks): the downwash footprint tracks the prescribed horizontal path
smoothly straight through the rank boundary, confirming
`apply_uav_forcing`'s per-rank local-cell contribution logic is correct,
not just untested. Phases 3+ (mean wind, ABL turbulence) are still
design-only, sketched below.

Goal: study how ambient turbulence and mean wind affect (and are affected
by) a small rotorcraft UAV during takeoff, landing, and a prescribed flight
trajectory close to the ground, using fdm-dopamine's existing LES/ABL
machinery.

## 1. Why an actuator disk, not a moving IBM body

The existing IBM (`src/ibm.f90`) is a static, SDF-based, no-slip ghost-cell
method: geometry is baked into `phi` and the `ghost_u/v/w` index lists once
at `setup_ibm`, then reused every step. Rebuilding those lists every step for
a body translating through the domain would mean re-ray-casting an SDF every
RK sub-step — expensive, and physically the wrong model anyway: we don't
need boundary-layer-resolved flow over the airframe, we need the rotor's
momentum footprint (induced velocity, wake, ground-effect recirculation) and
how ambient turbulence perturbs the disk loading. That is exactly what an
actuator-disk **body-force source term** gives, and it costs almost nothing
to move.

Reference practice: Sørensen & Shen (2002), Martinez-Tossas & Meneveau
(2017), and the SOWFA/PALM wind-turbine actuator-disk implementations all do
this — a rotor plane is a Lagrangian marker set that projects a momentum
sink/source into the momentum RHS via a regularized delta (Gaussian)
kernel; no solid boundary condition is imposed.

## 2. Model

### 2.1 Disk kinematics

- Input: a list of waypoints `(t_i, x_i, y_i, z_i)` — or `(s_i, x_i,y_i,z_i)`
  with a separately specified speed profile for takeoff/hover/landing
  segments (hover has zero path speed but the disk still needs a defined
  orientation).
- Interpolate position with a cubic spline (or PCHIP, to avoid overshoot near
  waypoints — important for a hover-to-climb transition) evaluated at every
  RK3 sub-stage time `to` used in `compute_time_step_RK3`.
- Orientation: default the thrust vector to `+y` (vertical) for
  takeoff/landing/hover segments, tilted toward the path tangent for cruise
  segments — simplest is a prescribed schedule (waypoint metadata carries a
  disk-tilt angle) rather than deriving it from a flight-dynamics model, at
  least in phase 1 (see §4).

### 2.2 Force

Two loading closures, selectable:

1. **Prescribed thrust** `T(t)` (a lookup table or constant hover thrust
   trimmed to weight) — decouples the study from a rotor aerodynamic model
   entirely; good for the first pass at "how does ambient turbulence disturb
   a UAV with a fixed control input".
2. **Momentum-theory closure** `T = 2*rho*A*a*(1-a)*|U_local + U_induced|^2`
   type relation (Rankine-Froude / Glauert), sampling `U_local` at the disk
   from the resolved+SGS field — lets ambient gusts feed back into
   instantaneous thrust/torque, which is the more interesting case for a
   "turbulence impact on the UAV" study but needs a stable iteration for the
   induction factor each step (cheap: a handful of fixed-point iterations
   over a few hundred markers).

Loading distribution across the disk: uniform to start; a simple radial
weighting (e.g. Prandtl tip correction or a specified `C_T(r)` profile) is a
drop-in refinement of the marker-weight table, not a structural change.

### 2.3 Projection onto the grid

- Marker set: ~15 radial x ~24 azimuthal points (~350) on the instantaneous
  disk, generated once in a disk-local frame at start-up and rotated +
  translated every step (cheap: a 3x3 rotation + translate per marker).
- Regularized delta: 3-point (Roma-Peskin) or truncated-Gaussian kernel,
  support radius ~2 grid cells in each direction, matching the local grid
  spacing at the marker's position (grid is stretched in `y`, see
  `genGridandIC.f90`, so the kernel width must be evaluated locally, not
  assumed uniform).
- Force is added into `rhs_u`, `rhs_v`, `rhs_w` in `equations.f90`, the same
  place the Boussinesq buoyancy term is added in `compute_rhs_v` (line ~314):
  a per-cell increment computed from the marker forces that fall in that
  cell's kernel support, accumulated as a source array (or, cheaper for a
  disk this small, looped marker-by-marker with atomic/reduction adds).
- Rank ownership: before projecting, test each marker's support box against
  this rank's `(ig1_global,ig2_global,kg1_global,kg2_global)` bounds (same
  pattern `ibm.f90`'s distributed-field readers already use) and skip
  markers with no overlap — O(n_markers) per rank, negligible next to the
  Poisson solve.

### 2.4 Where this plugs into the RK3 loop

`compute_time_step_RK3` (`time_integration.f90`) currently calls
`apply_ghost_cell_ibm` at each of the 3 RK stages plus the projection
re-enforce step. The UAV force term is different in kind — it's an explicit
RHS contribution, not a velocity constraint — so it belongs inside
`compute_rhs_u/v/w` (called once per stage, before the pressure projection),
evaluated at that stage's `to`. No re-enforcement step is needed afterward:
unlike a Dirichlet ghost-cell BC, a momentum source is not corrupted by the
projection's pressure-gradient correction (only its divergence content is
absorbed, exactly as with any other explicit force).

## 3. Sketch of a new module `src/uav_actuator.f90`

```
Module uav_actuator
  ! path_x(:), path_y(:), path_z(:), path_t(:): waypoints, read from a plain
  !   text file (mirrors inflow_profile_file's format conventions)
  ! spline coefficients built once at setup_uav (cubic natural spline in t)
  ! disk_radius, n_r, n_theta, marker arrays (local-frame + current world-frame)
  ! thrust_mode: 0=prescribed T(t) table, 1=momentum-theory closure
  !
  ! Subroutine setup_uav              -- read path/params, build markers, !$acc enter data
  ! Subroutine uav_disk_state(t, xc,yc,zc, that)   -- spline eval: centre + thrust unit vector
  ! Subroutine uav_place_markers(t)   -- rotate/translate local markers into world frame this stage
  ! Subroutine uav_compute_forces(U_,V_,W_, t)  -- momentum-theory sampling (mode 1 only)
  ! Subroutine uav_add_rhs(rhs_u,rhs_v,rhs_w)   -- project marker forces via regularized delta
End Module
```

Namelist sketch (`&UAV`, mirroring the existing `&IBM`/`&INFLOW` style):

```
&UAV
  uav_active       = 0,          ! 0=off (default), 1=on
  uav_path_file    = 'uav_path.dat',   ! columns: t x y z [tilt_deg]
  uav_disk_radius  = 0.15,       ! [m]
  uav_n_r          = 15,
  uav_n_theta      = 24,
  uav_thrust_mode  = 0,          ! 0=prescribed T(t), 1=momentum-theory
  uav_thrust_file  = 'uav_thrust.dat', ! (mode 0) columns: t T
  uav_hover_thrust = 0.0,        ! (mode 1) trim value, N -- weight estimate
  uav_kernel_ncell = 2,          ! regularized-delta support radius, grid cells
/
```

## 4. Turbulence / wind environment — reuse, don't rebuild

This is the part that's already almost entirely in place:

- **Mean wind + shear**: `flat_wall_model_flag=1` (smooth) or the rough-wall
  `ibm_z0`/MOST path (`wallmodel.f90`) already gives a log-law near-ground
  profile; `dPdx`/`Ub_target` (constant-mass-flux forcing) sets the mean
  wind speed aloft.
- **Turbulence injection**: `&INFLOW inflow_type=1` (ESEM synthetic-eddy
  method) lets us set an atmospheric turbulence-intensity/length-scale
  profile directly (`sem_profile_format=1`, wind-tunnel-style TI input) —
  good for a parametric "vary ambient TI, see UAV disk-load variance" sweep
  without running a separate precursor.
- **More realistic ABL turbulence**: `inflow_type=2` (recycled precursor)
  — run a precursor ABL/rough-wall channel once, then feed its recorded
  x-normal slice into the UAV domain as time-resolved inflow (this is
  exactly the `examples/precursor_successor` pattern already in the repo).
  This is the higher-fidelity option for gust-response studies, at the cost
  of a separate precursor run.
- **Ground**: reuse the existing rough-wall EQWM/MOST ground BC as-is (no
  changes needed) — `ibm_z0` or `flat_wall_model_flag` already models the
  surface the UAV takes off from/lands on.
- **Ground-effect / dust lofting**: `&SEDIMENT` (Soulsby settling scalar) is
  a ready-made proxy for dust/debris lofted by rotor downwash near the
  ground during takeoff/landing — set `sed_bc_bot` and seed `C_ref` near the
  disk's touchdown footprint rather than adding a new scalar.

So "vary turbulence and wind levels" in the study becomes: sweep
`sem_length_scale`/TI profile (or precursor `Ub_target`/roughness) and
`dPdx`/`Ub_target`, holding the UAV path and thrust schedule fixed, and
compare disk-load time series and near-ground wake statistics across runs.

## 5. Diagnostics (all pre-existing)

- `&STATISTICS` line/slice probes centered on the disk's footprint and
  along the flight path, for velocity/pressure time series feeding a
  disk-load or induced-velocity power spectrum.
- `reynolds_stress_budget.f90` (`rsb_active=1`) around the ground-effect
  region during hover/landing, to quantify how the rotor wake modifies
  local turbulence production/dissipation vs. the undisturbed ABL.
- A UAV analogue of `compute_ibm_forces`'s force accumulator: log
  instantaneous thrust, disk-integrated induced velocity, and (mode 1) the
  induction factor per step to a small CSV, the same way `write_force_csv`
  already does for IBM.

## 6. Phased implementation plan

1. **Static hover validation**: fixed disk in quiescent flow -- **done**:
   `src/uav_actuator.f90` + `examples/uav_hover_disk` (free-slip box, no
   ground). Ground-effect check -- **done**:
   `examples/uav_ground_effect/{near_ground,far_ground}` place the same
   disk at h/R=2 and h/R=8 above the domain's own no-slip y=0 wall (no IBM
   needed -- the ground is already grid-conforming). Both runs reach a
   quasi-steady state by ~step 400 (std < 1% of the mean over the last 10
   line-probe snapshots). Comparing the downward velocity just below the
   disk: near_ground/far_ground = 0.94, i.e. ~6% *less* induced velocity in
   ground effect at the same prescribed thrust -- the correct qualitative
   direction per classical actuator-disk ground-effect theory (Cheeseman &
   Bennett's image-disk model predicts a smaller ~1.6% reduction at this
   h/R for an idealized inviscid disk; our larger reduction is plausible
   given a real, compact, turbulent LES box rather than an idealized
   infinite-domain potential-flow disk, but this is a qualitative match, not
   a quantitative validation against the closed-form theory). See
   `examples/uav_ground_effect/analyze_ground_effect.py`. Remaining honest
   gap: a real quantitative comparison would need a longer time-averaging
   window, an ensemble of realizations (turbulent flow, single run), and
   ideally a comparison at several h/R rather than two points.
2. **Prescribed path, quiescent ambient** -- **done**, including the
   cross-rank case: `examples/uav_path_takeoff` (vertical climb, single
   rank) and `examples/uav_path_cross_rank` (horizontal path, `p_row=2`,
   `mpirun -np 2`, disk travels through x=0.5, the rank boundary). Both
   stay stable (divergence at machine precision throughout, including the
   crossing window) and an x-line probe through the disk's height/z shows
   the downwash footprint tracking the prescribed x(t) smoothly straight
   through the rank boundary -- no discontinuity, jump, or missing/doubled
   force at x=0.5. This exercises exactly the code path a single-rank run
   cannot: a marker whose kernel support spans two ranks, each
   contributing only the local cells it actually owns. Still open: a full
   takeoff -> cruise -> landing profile combining vertical and horizontal
   motion in one path (each has only been exercised separately so far).
3. **Add mean wind**: steady crosswind via `dPdx`/`Ub_target`; check wake
   deflection and asymmetric loading make physical sense.
4. **Add ABL turbulence**: SEM (`inflow_type=1`) first (fast, parametric),
   then precursor-recycled (`inflow_type=2`) for a higher-fidelity gust
   case; sweep turbulence intensity/length scale and quantify disk-load
   variance and near-ground TKE changes during takeoff/landing.
5. **Stretch — momentum-theory closure + light dynamics**: switch to
   `uav_thrust_mode=1` so loading responds to local gusts, and optionally
   let position deviate from the prescribed path via a simple point-mass
   force-balance (thrust vector + drag + weight) instead of pure
   kinematics — this is a genuine scope increase (needs a small ODE
   integrator coupled into the RK3 loop) and should only be taken on after
   phases 1-4 are validated.

## 7. Cost estimate

Per-step added cost is `O(n_markers x kernel_support)` — a few thousand
flops — plus an `O(n_markers)` rank-ownership test; both are negligible next
to the Poisson solve that already dominates runtime. The real added cost of
this study is not the actuator disk itself but:
- resolution needed near the ground/disk to resolve the wake and
  ground-effect recirculation (may need local grid stretching/clustering
  along the flight path, not just at the wall — a domain-size/grid-design
  decision, not a solver-cost one), and
- the ABL turbulence generation path chosen: SEM is essentially free
  (already amortized), whereas a full precursor run doubles the total
  simulation cost (a separate run) for higher-fidelity inflow.

## 8. Open questions

- Disk-tilt schedule for cruise segments: prescribed per-waypoint vs.
  derived from path tangent — affects how the marker-rotation code is
  parametrized.
- Whether ground-effect fidelity requires resolving the disk boundary layer
  at all (probably not, for a first pass — momentum-source models are known
  to under-predict very-near-wake details but capture the induced-velocity
  field and far-wake correctly).
- How far to take the two-way coupling (phase 5) before it's really a
  flight-dynamics + CFD co-simulation rather than a CFD study with a moving
  forcing term.
