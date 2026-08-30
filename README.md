# dopamine-fdm

A finite-difference Navier–Stokes solver for turbulent channel and open-channel flows, with support for rough-wall immersed boundary methods (IBM), wall models, sub-grid scale (SGS) turbulence modelling, and an exact Reynolds stress budget analysis module.

## Features

- **3-D incompressible Navier–Stokes** — fractional-step projection method with a spectral (FFTW3-MPI) Poisson solver for the pressure correction
- **MPI parallelism** — z-slab domain decomposition via MPI; scales to large core counts
- **Staggered MAC grid** — U on x-faces, V on y-faces, W on z-faces, pressure at cell centres
- **Flexible vertical grid** — 7 stretching options (uniform, symmetric tanh, single-sided tanh, uniform roughness sublayer + tanh blend)
- **Immersed boundary method (IBM)** — ghost-cell method (Tseng & Ferziger 2003) for complex rough surfaces; accepts a binary face-point mask or a precomputed signed-distance field (SDF)
- **Wall models** — flat-wall and IBM-surface equilibrium wall models (log-law EQWM); selectable per wall
- **SGS turbulence model** — Vreman (2004) eddy-viscosity model; optional (set `sgs_model = 0` for DNS)
- **Flexible boundary conditions** — Dirichlet (no-slip) or Neumann (free-slip) independently on top and bottom y-walls
- **Initial conditions** — log-law profile, linear tent/ramp, inverse-linear anti-tent/inverted-ramp, zero mean, or Reichardt turbulent channel profile; all with configurable white-noise perturbation
- **Oscillatory / pulsatile forcing** — independent time-varying pressure gradients in $x$ and $z$: $f_x(t) = \mathrm{d}P/\mathrm{d}x_0 + U_{b,x}\,\omega_x\cos(\omega_x t)$, $f_z(t) = \mathrm{d}P/\mathrm{d}z_0 + U_{b,z}\,\omega_z\cos(\omega_z t)$; wave periods `T_wave_x` / `T_wave_z` are supplied directly and converted to $\omega = 2\pi/T$ internally; set to `0` for steady forcing; enables cross-wave boundary layer studies
- **Suspended-sediment scalar transport** — passive scalar *C* (cell-centred) advected with a van Leer (1974) slope-limited MUSCL scheme (TVD, second-order on non-uniform grids); Soulsby (1997) settling velocity $w_s$; effective diffusivity $\kappa = \nu/Sc + \nu_t/Sc_t$; conservative cell-width flux normalisation; IBM-aware RHS masking; RK3 time-integration
- **Reynolds stress budget** — windowed ensemble averaging of all six Pope §7.4 budget terms (production, pressure-strain, viscous diffusion, turbulent diffusion, pressure diffusion, resolved and SGS dissipation); growing binary output with restart-safe sample count
- **Binary I/O** — Fortran stream unformatted, big-endian float64 for field snapshots; hot-restart supported; file-size sanity check on restart
- **Post-processing** — `postProcessing/generateXMF.py` writes XDMF files for ParaView, including a combined velocity vector attribute via XDMF `JOIN`

## Repository layout

```
fdm-dopamine/
├── CMakeLists.txt           # CMake build definition
├── cmake/
│   └── FindFFTW.cmake       # FFTW3 finder module
├── src/                     # Fortran 90 source files
│   ├── mpi.f90              # MPI setup and domain decomposition
│   ├── profiler.f90         # Per-stage MPI_Wtime() profiler (printed at shutdown)
│   ├── global.f90           # Shared global variables
│   ├── interpolation.f90    # Velocity interpolation utilities
│   ├── sgs_model.f90        # Vreman SGS model
│   ├── equations.f90        # RHS convection + diffusion
│   ├── boundary_conditions.f90
│   ├── readMask.f90         # IBM mask / SDF reader
│   ├── ibm.f90              # Ghost-cell IBM module
│   ├── genGridandIC.f90     # Grid generation and initial conditions
│   ├── input_output.f90     # Namelist reader and field I/O
│   ├── initialization.f90
│   ├── scalar_transport.f90 # Suspended sediment: van Leer MUSCL advection, settling
│   ├── wallmodel.f90        # Flat-wall and IBM EQWM
│   ├── poisson_gpu.f90      # [ENABLE_GPU only] cuFFT + cuSPARSE GPU Poisson solve
│   ├── projection.f90       # Pressure projection (FFTW Poisson, or GPU when ENABLE_GPU=ON)
│   ├── time_integration.f90 # RK3 time stepping
│   ├── monitor.f90          # Runtime statistics
│   ├── reynolds_stress_budget.f90  # Pope §7.4 budget (all 6 components)
│   ├── finalization.f90
│   └── main.f90             # Entry point
├── postProcessing/
│   └── generateXMF.py       # Write XDMF metadata for ParaView
├── preProcessing/
│   └── GenSDF/              # Signed-distance field generator for IBM
├── input_parameters         # Runtime namelist (edit before running)
└── .gitignore
```

## Dependencies

| Library | Minimum version | Purpose |
|---------|----------------|---------|
| gfortran or ifort | GFortran ≥ 9 / Intel ≥ 2019 | Fortran compiler |
| MPI | any standard MPI-3 | Domain decomposition |
| FFTW3 | 3.3 with MPI support | Spectral Poisson solver |
| LAPACK / BLAS | any | Linear algebra |
| CMake | 3.18 | Build system |
| NVHPC SDK *(optional, GPU build only)* | 23.3+ | `nvfortran` + OpenACC + cuFFT + cuSPARSE |

On Debian/Ubuntu the CPU-build dependencies can be installed with:

```bash
sudo apt install gfortran libopenmpi-dev libfftw3-dev libfftw3-mpi-dev liblapack-dev libblas-dev cmake
```

## Building (CPU)

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

The executable `dopamine` is placed in `build/`.

> **Note**: The solver uses `-fconvert=big-endian` globally so that field snapshots are big-endian.
> Reynolds stress budget output uses `CONVERT='little_endian'` on the file `OPEN` to produce
> little-endian float64 (consistent with the coordinate `.bin` files used by `generateXMF.py`).

## Building (single-GPU, OpenACC + cuFFT + cuSPARSE) - Not tested extensively - Watch for bugs!

An optional `ENABLE_GPU` CMake target offloads the RHS/SGS/boundary-condition/projection
kernels and the pressure Poisson solve (cuFFT transform + cuSPARSE batched tridiagonal solve)
to a single NVIDIA GPU via OpenACC. The CPU build above is untouched and remains the default;
this is a separate, opt-in build.

**Requirements**:
- [NVIDIA HPC SDK](https://developer.nvidia.com/hpc-sdk) 23.3 or later (provides `nvfortran`,
  its bundled OpenMPI, and the `cufft`/`cusparse` device libraries) — e.g. installed under
  `/opt/nvidia/hpc_sdk`.
- A CUDA-capable NVIDIA GPU. `CMakeLists.txt` targets compute capability `cc86` (Ampere,
  e.g. RTX 30-series); edit the `-gpu=cc86` flag in `CMakeLists.txt` to match a different GPU
  before building for another architecture.

**Scope — read before using**: the GPU Poisson solve currently only supports **`nprocs=1`**
(single GPU, no MPI domain decomposition), with either **`x_bc_type=0`** (periodic, spectral FFT)
or **`x_bc_type=1`** (inflow/outflow, DCT-IV) streamwise BC. A runtime guard `Stop`s immediately on
startup if `nprocs/=1` or `x_bc_type` is anything other than 0/1 — Every other solver feature (IBM, wall models,
SGS, sediment transport, RSB statistics, etc.) works normally in the GPU build.

Put the NVHPC SDK's `nvfortran` and bundled MPI on your `PATH`/`LD_LIBRARY_PATH` first, e.g.:

```bash
export PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/23.3/comm_libs/openmpi4/openmpi-4.0.5/bin:/opt/nvidia/hpc_sdk/Linux_x86_64/23.3/compilers/bin:$PATH
export LD_LIBRARY_PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/23.3/comm_libs/openmpi4/openmpi-4.0.5/lib:/opt/nvidia/hpc_sdk/Linux_x86_64/23.3/compilers/lib:$LD_LIBRARY_PATH
```

Then configure and build with `nvfortran` and `-DENABLE_GPU=ON`:

```bash
cmake -S . -B build_gpu -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_Fortran_COMPILER=nvfortran -DENABLE_GPU=ON
cmake --build build_gpu -j$(nproc)
```

The executable `dopamine` is placed in `build_gpu/`. `-Minfo=accel` output during the build
lists every OpenACC region the compiler offloaded — useful to confirm a kernel you touched is
still being generated for the GPU.

## Running

1. Edit `input_parameters` to set domain size, grid type, physics, and IC options (see comments
   in that file). For the GPU build, also make sure `nprocs=1` (i.e. launch with `-np 1`) and
   `x_bc_type` is 0 or 1 in `&BOUNDARY_CONDITIONS` — see "Building (single-GPU...)" above.
2. Create the required output directories (adjust paths to match your `fileout` and `rsb_fileout` settings):
   ```bash
   mkdir -p fields restart stats
   ```
3. Launch with MPI (CPU build, any rank count):
   ```bash
   mpirun -np <N> ./build/dopamine
   ```
   or, for the GPU build (single rank only):
   ```bash
   mpirun -np 1 ./build_gpu/dopamine
   ```
   The solver reads `input_parameters` from the working directory in both cases.

Every run prints a per-stage profiler summary (`src/profiler.f90`) at shutdown — wall time and
percent-of-tracked-time for SGS, wall model, RHS, boundary conditions, Poisson FFT, Poisson
tridiagonal solve, projection, and CFL check — useful for spotting where time is going on either
build.

## Input parameters overview

Parameters are grouped into Fortran namelists in `input_parameters`.  Any namelist group marked *optional* may be omitted entirely; the solver will use the documented defaults.

### `&DOMAIN`
| Parameter | Default | Description |
|-----------|---------|-------------|
| `nx, ny, nz` | — | Total face-point counts (interior cells = n−1) |
| `Lx, Ly, Lz` | — | Domain lengths |
| `grid_type` | `1` | Vertical grid: 1=uniform, 2=symmetric tanh, 3=tanh bottom, 4=tanh top, 5–7=sublayer+tanh variants |
| `alpha_grid` | `1.0` | Grid-stretching intensity (larger → stronger clustering) |

### `&PHYSICS`
| Parameter | Default | Description |
|-----------|---------|-------------|
| `nu` | — | Kinematic viscosity |
| `dPdx` | — | Mean (steady) streamwise pressure gradient |
| `dPdz` | `0.0` | Mean (steady) spanwise pressure gradient |
| `Ub_x` | `0.0` | Streamwise oscillatory velocity amplitude |
| `Ub_z` | `0.0` | Spanwise oscillatory velocity amplitude |
| `T_wave_x` | `0.0` | Wave period $T_x$ [s] for streamwise forcing; `0` = steady; $\omega_x = 2\pi/T_x$ derived internally |
| `T_wave_z` | `0.0` | Wave period $T_z$ [s] for spanwise forcing; `0` = steady; $\omega_z = 2\pi/T_z$ derived internally |
| `phi_wave_x` | `0.0` | Phase offset $\varphi_x$ [rad] for streamwise oscillation (e.g. `1.5708` ≈ π/2 gives sine forcing) |
| `phi_wave_z` | `0.0` | Phase offset $\varphi_z$ [rad] for spanwise oscillation; difference $\varphi_z - \varphi_x$ sets the cross-wave phase lag |
| `sgs_model` | `0` | 0=DNS (ν_t=0), 1=Vreman SGS |
| `Cs_vreman` | `0.1` | Smagorinsky-equivalent constant for Vreman model ($c_V = 2.5\,C_s^2$) |
| `flat_wall_model_flag` | `0` | 0=no-slip, 1=log-law EQWM on flat walls |

### `&NUMERICS`
| Parameter | Default | Description |
|-----------|---------|-------------|
| `dt` | — | Time-step size |
| `nsteps` | — | Total number of time steps. If `nsteps < 0`, the step count is ignored and the run instead stops once physical time `t >= sim_end_time` |
| `nsave` | — | Write field snapshot every `nsave` steps. If `nsave < 0`, the step count is ignored and a snapshot is instead written every `tsave` physical time units |
| `nmonitor` | `1` | Print monitor line every `nmonitor` steps |
| `sim_end_time` | `1e30` | Physical end time; only used when `nsteps < 0` |
| `tsave` | `1e30` | Physical save interval; only used when `nsave < 0`. When `cfl_adaptive=1` (or `dt` otherwise doesn't divide evenly into `tsave`), `dt` is automatically shrunk on the step just before a save so `t` lands exactly on multiples of `tsave` |
| `cfl_adaptive` | `0` | 0=fixed `dt`, 1=CFL-adaptive time step |
| `cfl_target` | `0.5` | Target CFL number (adaptive mode) |
| `cfl_safety` | `0.9` | Safety factor on adaptive `dt` adjustment |
| `dt_min` | `1e-8` | Minimum allowed `dt` (adaptive mode) |
| `dt_max` | `1e10` | Maximum allowed `dt` (adaptive mode) |

### `&BOUNDARY_CONDITIONS`
| Parameter | Default | Description |
|-----------|---------|-------------|
| `bc_face_ylo` | `1` | Bottom wall: 1=no-slip (Dirichlet), 2=free-slip (Neumann) |
| `bc_face_yhi` | `1` | Top wall: 1=no-slip (Dirichlet), 2=free-slip (Neumann) |
| `x_bc_type` | `0` | Streamwise BC: 0=periodic (spectral FFT pressure solve); 1=inflow/outflow (Dirichlet velocity at inflow, see `&INFLOW`; convective outflow; DCT-IV pressure solve). GPU build supports both, `nprocs=1` only |

### `&INFLOW` *(optional — used only when `x_bc_type = 1`)*
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


### `&IBM` *(optional)*
| Parameter | Default | Description |
|-----------|---------|-------------|
| `ibm_input_mode` | `0` | 0=no body; 1=read `Umask_in` (face-point binary mask, φ computed via fast-sweep); 2=read `SDF_in` only (precomputed cell-centre signed-distance, solid/fluid from sign(φ)) |
| `ibm_wall_model_flag` | `0` | 0=no-slip, 1=log-law EQWM on IBM surfaces |
| `ibm_mask_file` | `'Umask_in'` | Path to face-point binary mask (mode 1) |
| `ibm_sdf_file` | `'SDF_in'` | Path to precomputed cell-centre SDF (mode 2) |
| `ks` | — | Roughness sublayer height (grid types 5–7) |
| `nks` | — | Number of uniform cells in roughness sublayer |
| `nsampling` | `0` | IBM force output interval (0 = disabled) |
| `ibm_surface_nsampling` | `0` | Per-point surface field (pressure + pressure/viscous force) dump interval, written to `ibm_surface/surface.<step>.bin` (0 = disabled) |

### `&INITIAL_CONDITIONS`
| Parameter | Default | Description |
|-----------|---------|-------------|
| `ic_type` | `1` | 1=log-law + noise; 2=linear/tent + noise; 3=zero mean + noise; 4=Reichardt turbulent channel profile + structured perturbation; 5=inverse-linear/anti-tent + noise |
| `noise_percent` | `5.0` | White-noise amplitude as % of `Utarget` (applied to U, V, W for types 1–3) |
| `Utarget` | — | Target bulk or centreline velocity |
| `nstep_init` | `0` | Starting step number (non-zero for hot-start logging) |
| `restart` | `0` | 0=fresh run, 1=hot-start from `filein` |

**IC type 4 — Reichardt profile** (`ic_type = 4`):  
Builds the Reichardt (1951) composite law-of-the-wall profile, scaled by $u_\tau$ derived from `dPdx` and the channel half-height, mirrored symmetrically for a two-wall channel.  The peak centreline velocity is $U_{cl}^+ \cdot u_\tau \approx 18\,\text{m/s}$ for $Re_\tau \approx 395$.  Perturbations are applied as structured streamwise/spanwise waves rather than white noise.

**IC type 5 — inverse-linear profile** (`ic_type = 5`):  
Sets a mean streamwise velocity that is large near the walls and decreases linearly to zero at the channel centre (or top boundary), imposing a strong destabilising shear throughout the domain to promote rapid laminar-to-turbulent transition:
$$U(y) = \begin{cases} U_\text{target}\left(1 - \dfrac{2\min(y,\,L_y-y)}{L_y}\right) & \text{closed channel (no-slip top)} \\[6pt] U_\text{target}\left(1 - \dfrac{y}{L_y}\right) & \text{open channel (free-slip top)} \end{cases}$$
White noise of amplitude `noise_percent`% is superimposed on U, V, and W.

### `&IO`
| Parameter | Default | Description |
|-----------|---------|-------------|
| `filein` | — | Path to restart field file (used when `restart = 1`) |
| `fileout` | — | Base name for output snapshot files |

### `&SEDIMENT` *(optional — omit to disable)*
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

### `&STATISTICS` *(optional — omit to disable)*

Controls the Reynolds stress budget (RSB) module.  When enabled, the module accumulates windowed time-averages of all moments needed for the Pope §7.4 budget every `rsb_freq` steps, then writes one slice to each output file at the end of each window.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `rsb_active` | `0` | 0=disabled, 1=enable RSB accumulation and output |
| `rsb_freq` | `10` | Accumulate every `rsb_freq` steps; write (and reset) after each accumulation |
| `rsb_nstart` | `0` | Step at which to begin accumulating (skip transient spin-up) |
| `rsb_hom_dir` | `'x,z'` | Comma-separated homogeneous directions for spatial averaging: `''`=none, `'x'`, `'z'`, `'x,z'` |
| `rsb_fileout` | `'rsb'` | Base path for output files (e.g. `'stats/rsb'`) |

#### RSB output files

All files are little-endian float64 (no Fortran record markers), produced by rank 0 only.  Each call to `output_rsb` appends one slice of shape `(nc, out_nx, out_ny, out_nz)` to each `.bin` file, so the on-disk layout is `(nc, out_nx, out_ny, out_nz, nsamples)`.

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

**`out_nx`, `out_ny`, `out_nz`** are the output grid sizes after homogeneous averaging: 1 along any homogeneous direction, full cell-count otherwise.

**Restart continuity**: on a hot-start (`restart = 1`) the module reads `nsamples` from the existing `.meta` file and continues appending.  A fresh run (`restart = 0`) truncates all pre-existing `.bin` files.

#### Reading RSB output in Python

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

## Output files

| Location | Content | Format |
|----------|---------|--------|
| `fields/` | Velocity (U, V, W), pressure (P), and — when `sgs_model /= 0` — SGS turbulent viscosity (ν_t) snapshots, written every `nsave` steps (or every `tsave` time units if `nsave < 0`) | Big-endian float64 stream |
| `restart/` | Hot-restart fields | Big-endian float64 stream |
| `stats/` | Monitor statistics (text) and, if enabled, RSB budget files (little-endian float64) | See above |

Field snapshots are Fortran stream unformatted, big-endian float64 (no record markers).  Each field block is preceded by a 3-integer size header.  The layout is:

| Block | Dimensions | Present |
|-------|-----------|---------|
| Grid arrays (x, y, z, xm, ym, zm) | 1-D, each preceded by 1 Int32 | always |
| U | `nx × nyg × nzg` | always |
| V | `nxg × ny × nzg` | always |
| W | `nxg × nyg × nz` | always |
| P | `nxg × nyg × nzg` | always |
| C | `nxg × nyg × nzg` | only when `sediment_flag >= 1` |
| ν_t | `nxg × nyg × nzg` | only when `sgs_model /= 0` |

```python
import numpy as np
data = np.fromfile("fields/channel_test.1", dtype=">f8")
```

### XDMF / ParaView post-processing

`postProcessing/generateXMF.py` generates XDMF metadata files that let ParaView open the binary snapshots directly.  It produces a single `Velocity` vector attribute using the XDMF `JOIN` function rather than three separate scalars:

```bash
python postProcessing/generateXMF.py --case channel_test --nx 513 --ny 128 --nz 257
```

## References

The methods implemented in this solver draw on the following published works.  Please cite the relevant papers when publishing results obtained with `fdm-dopamine`.

| Method | Citation | DOI |
|--------|----------|-----|
| Time integration (RK3) | Wray, A.A. (1990). *Minimal storage time advancement schemes for spectral methods*. NASA Ames Report. | — |
| Fractional-step projection | Kim, J. & Moin, P. (1985). J. Comput. Phys. **59**, 308–323. | [10.1016/0021-9991(85)90148-2](https://doi.org/10.1016/0021-9991(85)90148-2) |
| Spectral Poisson solver (FFTW3) | Frigo, M. & Johnson, S.G. (2005). Proc. IEEE **93**(2), 216–231. | [10.1109/JPROC.2004.840301](https://doi.org/10.1109/JPROC.2004.840301) |
| Ghost-cell IBM | Tseng, Y.-H. & Ferziger, J.H. (2003). J. Comput. Phys. **192**(2), 593–623. | [doi:10.1016/j.jcp.2003.07.024](https://doi.org/10.1016/j.jcp.2003.07.024) |
| Vreman SGS model | Vreman, A.W. (2004). Phys. Fluids **16**(10), 3670–3681. | [10.1063/1.1785131](https://doi.org/10.1063/1.1785131) |
| van Leer MUSCL limiter (scalar) | van Leer, B. (1974). J. Comput. Phys. **14**(4), 361–370. | [10.1016/0021-9991(74)90019-9](https://doi.org/10.1016/0021-9991(74)90019-9) |
| Settling velocity (sediment) | Soulsby, R.L. (1997). *Dynamics of Marine Sands*. Thomas Telford. | ISBN 978-0-7277-2584-5 |
| Reynolds stress budget | Pope, S.B. (2000). *Turbulent Flows*. Cambridge University Press. | [10.1017/CBO9780511840531](https://doi.org/10.1017/CBO9780511840531) |
| Reichardt IC profile | Reichardt, H. (1951). Z. Angew. Math. Mech. **31**(7), 208–219. | [10.1002/zamm.19510310704](https://doi.org/10.1002/zamm.19510310704) |
| Staggered MAC grid | Harlow, F.H. & Welch, J.E. (1965). Phys. Fluids **8**(12), 2182–2189. | [10.1063/1.1761178](https://doi.org/10.1063/1.1761178) |


## License

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

See [LICENSE](LICENSE) for the full text.

## Acknowledgments

The original form of the channel flow solver was shared by Adrian Lozano-Duran during my PhD and I am grateful for his generous help and support with the code.

## Citing this solver

If you use `dopamine-fdm` in your research, please consider citing the following publications.

1. Lozano-Durán, A., & Bae, H. J. (2019). Characteristic scales of Townsend's wall-attached eddies. *Journal of Fluid Mechanics*, 868, 698–725. doi:[10.1017/jfm.2019.209](https://doi.org/10.1017/jfm.2019.209) — *original DNS solver*
2. Patil, A., & Fringer, O. (2022). Drag enhancement by the addition of weak waves to a wave-current boundary layer over bumpy walls. *Journal of Fluid Mechanics*, 947, A3. doi:[10.1017/jfm.2022.628](https://doi.org/10.1017/jfm.2022.628) — *IBM + DNS*
3. Patil, A., & García-Sánchez, C. (2025). Should we care about the spatial heterogeneity in coral reefs under unidirectional turbulent flows? arXiv:[2506.03021](https://arxiv.org/abs/2506.03021) [physics.flu-dyn] — *improved IBM*
4. Patil, A., Paranjothi, U. C. K., & García-Sánchez, C. (2025). GenSDF: An MPI-Fortran based signed-distance-field generator for computational fluid dynamics applications. *SoftwareX*, 30, 102117. — *GenSDF*

<details>
<summary>BibTeX</summary>

```bibtex
@article{lozanoduran2019characteristic,
  title={Characteristic scales of {T}ownsend's wall-attached eddies},
  author={Lozano-Dur{\'a}n, Adri{\'a}n and Bae, Hyunji Jane},
  journal={Journal of Fluid Mechanics},
  volume={868},
  pages={698--725},
  year={2019},
  publisher={Cambridge University Press},
  doi={10.1017/jfm.2019.209}
}

@article{patil2022drag,
  title={Drag enhancement by the addition of weak waves to a wave-current boundary layer over bumpy walls},
  author={Patil, Akshay and Fringer, Oliver},
  journal={Journal of Fluid Mechanics},
  volume={947},
  pages={A3},
  year={2022},
  publisher={Cambridge University Press},
  doi={10.1017/jfm.2022.628}
}

@misc{patil2025carespatialheterogeneitycoral,
  title={Should we care about the spatial heterogeneity in coral reefs under unidirectional turbulent flows?},
  author={Patil, Akshay and Garc{\'i}a-S{\'a}nchez, Clara},
  year={2025},
  eprint={2506.03021},
  archivePrefix={arXiv},
  primaryClass={physics.flu-dyn},
  url={https://arxiv.org/abs/2506.03021}
}

@article{patil2025gensdf,
  title={{GenSDF}: An {MPI-Fortran} based signed-distance-field generator for computational fluid dynamics applications},
  author={Patil, Akshay and Paranjothi, Uma Chandrika Karrothu and Garc{\'i}a-S{\'a}nchez, Clara},
  journal={SoftwareX},
  volume={30},
  pages={102117},
  year={2025},
  publisher={Elsevier}
}
```

</details>
