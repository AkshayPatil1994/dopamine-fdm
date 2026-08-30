# Input Parameters Reference

← [Home](Home.md)

Parameters are grouped into Fortran namelists in `input_parameters`. Any namelist group
marked *optional* may be omitted entirely; the solver will use the documented defaults.
See also the runnable examples in `examples/` ([Examples](Examples.md)).

## `&DOMAIN`
| Parameter | Default | Description |
|-----------|---------|-------------|
| `nx, ny, nz` | — | Total face-point counts (interior cells = n−1) |
| `Lx, Ly, Lz` | — | Domain lengths |
| `grid_type` | `1` | Vertical grid: 1=uniform, 2=symmetric tanh, 3=tanh bottom, 4=tanh top, 5–7=sublayer+tanh variants |
| `alpha_grid` | `1.0` | Grid-stretching intensity (larger → stronger clustering) |
| `p_row, p_col` | `0, 0` | 2decomp&fft MPI pencil grid: `p_row` splits $x$, `p_col` splits $z$ ($y$ always local); `p_row * p_col` must equal the rank count. `0, 0` = auto-pick (prefer a pure z-slab split, `p_row=1`; fall back to a 2-D split matching the grid's aspect ratio). Set both explicitly to override |

See [Numerics § MPI parallelism](Numerics.md#10-mpi-parallelism).

## `&PHYSICS`
| Parameter | Default | Description |
|-----------|---------|-------------|
| `nu` | — | Kinematic viscosity |
| `dPdx` | — | Mean (steady) streamwise pressure gradient |
| `dPdz` | — | Mean (steady) spanwise pressure gradient |
| `Ub_x` | `0.0` | Streamwise oscillatory velocity amplitude |
| `Ub_z` | `0.0` | Spanwise oscillatory velocity amplitude |
| `T_wave_x` | `0.0` | Wave period $T_x$ [s] for streamwise forcing; `0` = steady; $\omega_x = 2\pi/T_x$ derived internally |
| `T_wave_z` | `0.0` | Wave period $T_z$ [s] for spanwise forcing; `0` = steady; $\omega_z = 2\pi/T_z$ derived internally |
| `phi_wave_x` | `0.0` | Phase offset $\varphi_x$ [rad] for streamwise oscillation (e.g. `1.5708` ≈ π/2 gives sine forcing) |
| `phi_wave_z` | `0.0` | Phase offset $\varphi_z$ [rad] for spanwise oscillation; difference $\varphi_z - \varphi_x$ sets the cross-wave phase lag |
| `sgs_model` | `0` | 0=DNS (ν_t=0), 1=Vreman SGS |
| `Cs_vreman` | `0.17` | Smagorinsky-equivalent constant for Vreman model ($c_V = 2.5\,C_s^2$) |
| `flat_wall_model_flag` | `0` | 0=no-slip, 1=log-law EQWM on flat walls |

See [Numerics § Forcing](Numerics.md#1-governing-equations) and
[§ SGS model](Numerics.md#5-sub-grid-scale-model).

## `&NUMERICS`
| Parameter | Default | Description |
|-----------|---------|-------------|
| `dt` | — | Time-step size |
| `nsteps` | — | Total number of time steps. If `nsteps < 0`, the step count is ignored and the run instead stops once physical time `t >= sim_end_time` |
| `nsave` | — | Write field snapshot every `nsave` steps. If `nsave < 0`, the step count is ignored and a snapshot is instead written every `tsave` physical time units |
| `nmonitor` | — | Print monitor line every `nmonitor` steps |
| `sim_end_time` | `1e30` | Physical end time; only used when `nsteps < 0` |
| `tsave` | `1e30` | Physical save interval; only used when `nsave < 0`. When `cfl_adaptive=1` (or `dt` otherwise doesn't divide evenly into `tsave`), `dt` is automatically shrunk on the step just before a save so `t` lands exactly on multiples of `tsave` |
| `cfl_adaptive` | `0` | 0=fixed `dt`, 1=CFL-adaptive time step |
| `cfl_target` | `0.5` | Target CFL number (adaptive mode) |
| `cfl_safety` | `0.9` | Safety factor on adaptive `dt` adjustment |
| `dt_min` | `1e-10` | Minimum allowed `dt` (adaptive mode) |
| `dt_max` | `1e10` | Maximum allowed `dt` (adaptive mode) |

See [Numerics § Time integration](Numerics.md#3-time-integration).

## `&BOUNDARY_CONDITIONS`
| Parameter | Default | Description |
|-----------|---------|-------------|
| `bc_face_ylo` | `1` | Bottom wall: 1=no-slip (Dirichlet), 2=free-slip (Neumann) |
| `bc_face_yhi` | `1` | Top wall: 1=no-slip (Dirichlet), 2=free-slip (Neumann) |
| `x_bc_type` | `0` | Streamwise BC: 0=periodic (spectral FFT pressure solve); 1=inflow/outflow (Dirichlet velocity at inflow, see `&INFLOW`; convective outflow; DCT-IV pressure solve). GPU build supports both, `nprocs=1` only |

## `&INFLOW` *(optional — used only when `x_bc_type = 1`)*
| Parameter | Default | Description |
|-----------|---------|-------------|
| `inflow_type` | `0` | 0=constant uniform flow (`U=inflow_Uconst`, `V=W=0`, no fluctuation); 1=synthetic eddy method (SEM) — mean `U(y)` and Reynolds-stress profiles from `inflow_profile_file` |
| `inflow_Uconst` | `0.0` | Uniform streamwise velocity (`inflow_type=0` only) |
| `inflow_profile_file` | `'inflow_profile.dat'` | Reference profile file (`inflow_type=1`): rows `y U V W uu vv ww uv [uw vw]`; only `y,U,uu,vv,ww,uv` are used (`uw=vw=0` assumed) |
| `sem_n_eddies` | `200` | Number of eddies in the virtual inflow box |
| `sem_length_scale` | `0.01` | Eddy shape-function radius [m] (homogeneous mode, i.e. `sem_sigma_file` unset) |
| `sem_seed` | `12345` | RNG seed for the eddy population (identical on every MPI rank) |
| `sem_ensemble_samples` | `100` | Samples per eddy-recycling period used to build the ESEM empirical normalisation (one-time start-up cost) |
| `sem_ensemble_periods` | `8` | Number of eddy-recycling periods the ensemble window spans |
| `sem_sigma_file` | `''` | Optional 10-column file (`y` + 9 `sigma_ij` components) for inhomogeneous eddy length scales (ESEM §4.1); empty = homogeneous |
| `sem_eddy_placement` | `0` | 0=uniform random placement; 1=PDF-weighted placement favouring smaller-eddy regions (ESEM §4.2, needs `sem_sigma_file`) |
| `sem_use_esem` | `1` | 1=Ensemble SEM (Schau et al. 2022) — empirical/exact normalisation, recommended; 0=classical SEM (Jarrin 2006) — analytical normalisation, comparison only |
| `sem_divergence_free` | `0` | 0=independent scalar eddy sums; 1=divergence-free construction (Poletto, Craft & Revell 2013) — curl of a vector potential, guarantees `div(u')=0` at the inflow plane. Requires `sem_use_esem=1` |

See [Numerics § Streamwise inflow/outflow BC](Numerics.md#11-streamwise-inflowoutflow-bc-x_bc_type--1)
and the [precursor/successor example](Examples.md#precursor_successor).

## `&IBM` *(optional)*
| Parameter | Default | Description |
|-----------|---------|-------------|
| `ibm_input_mode` | `0` | 0=no IBM body (smooth walls); 1=SDF-based ghost-cell IBM — reads the precomputed cell-centre signed-distance field `ibm_sdf_file`; solid/fluid mask `Umask_cc` is derived from `sign(phi)` |
| `ibm_wall_model_flag` | `0` | 0=no-slip, 1=log-law EQWM on IBM surfaces |
| `ibm_sdf_file` | `'SDF_in'` | Path to precomputed cell-centre SDF |
| `ibm_objid_file` | `''` | Optional per-solid object-ID field (e.g. `GenSDF`'s `sdfp_objid.bin`), for per-object boundary conditions; `''` = a single uniform IBM condition applied to every solid |
| `ks` | — | Roughness sublayer height (grid types 5–7) |
| `nks` | `0` | Number of uniform cells in roughness sublayer |
| `nsampling` | `0` | IBM force output interval (0 = disabled) |
| `ibm_surface_nsampling` | `0` | Per-point surface field (pressure + pressure/viscous force) dump interval, written to `ibm_surface/surface.<step>.bin` (0 = disabled) |

> **Note:** Single IBM path: `ibm_input_mode = 1` always
> reads a precomputed cell-centre SDF from `ibm_sdf_file`. Generate that SDF ahead of
> time with the `GenSDF` tool — see
> [Pre- and Post-Processing Tools § GenSDF](Tools.md#gensdf).

## `&INITIAL_CONDITIONS`
| Parameter | Default | Description |
|-----------|---------|-------------|
| `ic_type` | `1` | 1=log-law + noise; 2=linear/tent + noise; 3=zero mean + noise; 4=Reichardt turbulent channel profile + structured perturbation; 5=inverse-linear/anti-tent + noise |
| `noise_percent` | `5.0` | White-noise amplitude as % of `Utarget` (applied to U, V, W for types 1–3) |
| `Utarget` | — | Target bulk or centreline velocity |
| `nstep_init` | — | Starting step number (non-zero for hot-start logging) |
| `restart` | `0` | 0=fresh run, 1=hot-start from `filein` |

**IC type 4 — Reichardt profile** (`ic_type = 4`):
Builds the Reichardt (1951) composite law-of-the-wall profile, scaled by $u_\tau$ derived
from `dPdx` and the channel half-height, mirrored symmetrically for a two-wall channel.
The peak centreline velocity is $U_{cl}^+ \cdot u_\tau \approx 18\,\text{m/s}$ for
$Re_\tau \approx 395$. Perturbations are applied as structured streamwise/spanwise waves
rather than white noise.

**IC type 5 — inverse-linear profile** (`ic_type = 5`):
Sets a mean streamwise velocity that is large near the walls and decreases linearly to
zero at the channel centre (or top boundary), imposing a strong destabilising shear
throughout the domain to promote rapid laminar-to-turbulent transition:
$$U(y) = \begin{cases} U_\text{target}\left(1 - \dfrac{2\min(y,\,L_y-y)}{L_y}\right) & \text{closed channel (no-slip top)} \\[6pt] U_\text{target}\left(1 - \dfrac{y}{L_y}\right) & \text{open channel (free-slip top)} \end{cases}$$
White noise of amplitude `noise_percent`% is superimposed on U, V, and W.

See [Numerics § Initial conditions](Numerics.md#12-initial-conditions).

## `&IO`
| Parameter | Default | Description |
|-----------|---------|-------------|
| `filein` | — | Path to restart field file (used when `restart = 1`) |
| `fileout` | — | Base name for output snapshot files |

## `&SEDIMENT` *(optional — omit to disable)*
| Parameter | Default | Description |
|-----------|---------|-------------|
| `sediment_flag` | `0` | 0=off, 1=passive scalar with Soulsby settling |
| `d_s` | `1e-4` | Particle diameter [m] |
| `rho_s` | `2650.0` | Particle density [kg m⁻³] |
| `rho_f` | `1000.0` | Fluid density [kg m⁻³] |
| `grav` | `9.81` | Gravitational acceleration [m s⁻²] |
| `Sc` | `1.0` | Molecular Schmidt number |
| `Sc_t` | `0.7` | Turbulent Schmidt number |
| `C_ref` | `0.0` | Reference near-bed concentration |
| `sed_bc_bot` | `0` | Bottom BC for C: 0=zero normal flux, 1=equilibrium at `C_ref` |
| `C_ic_type` | `0` | Initial C profile: 0=uniform `C_ref`, 1=Rouse, 2=linear ramp, 3=slab |
| `C_ic_height` | `0.0` | Slab height [m] (used when `C_ic_type = 3`) |

See [Numerics § Scalar transport](Numerics.md#8-scalar-transport-suspended-sediment).

## `&STATISTICS` *(optional — omit to disable)*

Controls the Reynolds stress budget (RSB) module. When enabled, the module accumulates
windowed time-averages of all moments needed for the Pope §7.4 budget every `rsb_freq`
steps, then writes one slice to each output file at the end of each window.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `rsb_active` | `0` | 0=disabled, 1=enable RSB accumulation and output |
| `rsb_freq` | `10` | Accumulate every `rsb_freq` steps; write (and reset) after each accumulation |
| `rsb_nstart` | `0` | Step at which to begin accumulating (skip transient spin-up) |
| `rsb_hom_dir` | `'x,z'` | Comma-separated homogeneous directions for spatial averaging: `''`=none, `'x'`, `'z'`, `'x,z'` |
| `rsb_fileout` | `'rsb'` | Base path for output files (e.g. `'stats/rsb'`) |

### RSB output files

All files are little-endian float64 (no Fortran record markers), produced by rank 0
only. Each call to `output_rsb` appends one slice of shape `(nc, out_nx, out_ny, out_nz)`
to each `.bin` file, so the on-disk layout is `(nc, out_nx, out_ny, out_nz, nsamples)`.

| File | `nc` | Description |
|------|------|-------------|
| `<base>_Umean.bin` | 3 | Mean velocity vector $(\langle U \rangle, \langle V \rangle, \langle W \rangle)$ |
| `<base>_Rij.bin` | 6 | Reynolds stress tensor $R_{ij} = \langle u_i' u_j' \rangle$ |
| `<base>_Pij.bin` | 6 | Production $P_{ij}$ |
| `<base>_epsRes.bin` | 6 | Resolved dissipation $\varepsilon_{ij}^\text{res} = 2\nu\langle(\partial u_i'/\partial x_k)(\partial u_j'/\partial x_k)\rangle$ |
| `<base>_epsSGS.bin` | 6 | SGS dissipation $\varepsilon_{ij}^\text{sgs} = 2\langle\nu_t s_{ij}'\rangle$ |
| `<base>_PiStrain.bin` | 6 | Pressure-strain $\Pi_{ij}$ |
| `<base>_DTij.bin` | 6 | Turbulent diffusion $D^T_{ij}$ |
| `<base>_Dnuij.bin` | 6 | Viscous diffusion $D^\nu_{ij} = \nu\,\nabla^2 R_{ij}$ |
| `<base>_PhiPij.bin` | 6 | Pressure diffusion $\Phi^P_{ij}$ |
| `<base>_Resid.bin` | 6 | Budget residual (closure check) |
| `<base>.meta` | — | Text file: grid shape, `nsamples`, endian, dtype, component order |

The 6-component order for all tensor quantities is: **11, 22, 33, 12, 13, 23**.

**`out_nx`, `out_ny`, `out_nz`** are the output grid sizes after homogeneous averaging: 1
along any homogeneous direction, full cell-count otherwise.

**Restart continuity**: on a hot-start (`restart = 1`) the module reads `nsamples` from
the existing `.meta` file and continues appending. A fresh run (`restart = 0`) truncates
all pre-existing `.bin` files.

### Reading RSB output in Python

```python
import numpy as np

def load_rsb(base, term, nc, out_nx, out_ny, out_nz):
    """Load all samples for one RSB term.
    Returns array of shape (nc, out_nx, out_ny, out_nz, nsamples)."""
    raw = np.fromfile(f"{base}_{term}.bin", dtype="<f8")
    n = nc * out_nx * out_ny * out_nz
    nsamples = len(raw) // n
    return raw.reshape(nc, out_nx, out_ny, out_nz, nsamples, order='F')

# Example: 1-D channel profile (rsb_hom_dir='x,z'), nc=6 components
Rij = load_rsb("stats/rsb", "Rij", 6, 1, 128, 1)
# Rij[0, 0, :, 0, -1] → final-sample R_11 profile
```

Or use `postProcessing/read_RSBstats.py` — see
[Pre- and Post-Processing Tools](Tools.md#read_rsbstatspy).

See [Numerics § Reynolds stress budget](Numerics.md#9-reynolds-stress-budget).
