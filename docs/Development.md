# Development

For maintainers and contributors: repository layout, the CPU/GPU consistency guarantees and debugging tools, and how to contribute. Contents: [Code structure](#code-structure) · [Decomposition consistency](#decomposition-consistency-cpu--multi-cpu--gpu--multi-gpu) · [Contributing](#contributing)

## Code structure

Repository layout and the role of each Fortran module.

```
fdm-dopamine/
├── CMakeLists.txt           # CMake build definition
├── cmake/
│   └── FindFFTW.cmake       # FFTW3 finder module
├── src/                     # Fortran 90 source files
│   ├── mpi.f90              # MPI setup and domain decomposition
│   ├── profiler.f90         # Per-stage MPI_Wtime() profiler (printed at shutdown)
│   ├── global.f90           # Shared global variables
│   ├── decomp.f90           # 2decomp&fft pencil-decomposition bookkeeping (main + Poisson pencil grids)
│   ├── interpolation.f90    # Velocity interpolation utilities
│   ├── sem.f90              # Synthetic Eddy Method inflow (ESEM/SEM, recycling)
│   ├── sgs_model.f90        # Vreman SGS model
│   ├── equations.f90        # RHS convection + diffusion
│   ├── boundary_conditions.f90
│   ├── readMask.f90         # IBM setup driver (reads precomputed SDF)
│   ├── ibm.f90              # Ghost-cell IBM module (extended halo for image-point stencils)
│   ├── halo_pad.f90         # Host-side padded copies (extra x/z planes) of decomposed fields, used by the IBM and particles
│   ├── particles.f90        # Lagrangian point particles (tracers, inertial, IBM collisions, reinjection)
│   ├── debug_trace.f90      # DOPAMINE_TRACE_DIR stage tracing for layout-consistency debugging
│   ├── genGridandIC.f90     # Grid generation and initial conditions
│   ├── input_output.f90     # Namelist reader and field I/O
│   ├── initialization.f90
│   ├── scalar_transport.f90 # Suspended sediment: van Leer MUSCL advection, settling
│   ├── thermal_transport.f90 # Passive thermal scalar transport
│   ├── wallmodel.f90        # Flat-wall and IBM EQWM
│   ├── ibm_stress.f90       # Wall-model stress on the staircase faces of an immersed body (ibm_method=1)
│   ├── poisson_gpu.f90      # [ENABLE_GPU only] cuFFT + cuSPARSE GPU Poisson solve
│   ├── projection.f90       # Pressure projection (FFTW + 2decomp&fft pencil transposes, or GPU when ENABLE_GPU=ON)
│   ├── time_integration.f90 # RK3 time stepping
│   ├── monitor.f90          # Runtime statistics
│   ├── reynolds_stress_budget.f90  # Pope §7.4 budget (all 6 components)
│   ├── probe_output.f90     # Line/slice probe output
│   ├── uav_actuator.f90     # UAV actuator-disk rotor forcing (moving marker ring)
│   ├── waves.f90            # Target waves of the numerical wave flume: Airy, Rienecker-Fenton stream function, JONSWAP (inlet + relaxation zones)
│   ├── vof_plic.f90         # PLIC geometry: plane constants, volume below a plane in a box (no solver dependencies)
│   ├── vof_normals.f90      # Interface normals (Youngs, centred-column height function)
│   ├── vof_advect.f90       # Split-sweep Weymouth-Yue advection of C with exact PLIC fluxes (or THINC); serial kernels, unit-tested
│   ├── vof_state.f90        # C on the decomposed grid: halos, initial shapes, prescribed test fields, diagnostics (vof_diag.dat)
│   ├── vof_pressure.f90     # Variable-density pressure operator, geometric face density, PCG with the fast Poisson solver as preconditioner
│   ├── vof_mg.f90           # Mask-aware geometric multigrid preconditioner of the masked PCG (pcg_precond=1)
│   ├── vof_ibm.f90          # Hydrodynamic loads on immersed bodies in the two-fluid solver (ibm_forces.csv)
│   ├── vof_curv.f90         # Height-function curvature for surface tension
│   ├── mom_recon.f90        # Momentum face-value kernels (central, Koren, QUICK, WENO3-Z/WENO5-Z), pure functions
│   ├── vof_twofluid.f90     # Two-fluid time step: consistent momentum transport, forces (viscous, gravity, surface tension), projections
│   ├── gpu_device.f90       # [ENABLE_GPU only] rank-to-GPU binding for multi-GPU runs
│   ├── dopamine_esem_main.f90  # Standalone ESEM inflow precursor generator (dopamine-ESEM executable)
│   ├── finalization.f90
│   └── main.f90             # Entry point
├── tests/
│   ├── regression/          # Small deterministic cases: np / layout parity, restart, VOF cases (vof_check.py)
│   └── vof/                 # Standalone VOF unit tests: PLIC geometry, prescribed-velocity advection, momentum kernels
├── postProcessing/
│   ├── pyproject.toml       # dopamine-fdm-post: pip install -e postProcessing/
│   ├── dopamine_post/       # Post-processing library (see Tools page) + `dopamine-post` CLI
│   ├── examples/            # Short demonstration scripts, one per dopamine_post module
│   └── animations/          # Showcase GIF scripts (see Tools page)
├── preProcessing/
│   └── GenSDF/              # Signed-distance field generator for IBM
├── docs/                    # Documentation (this folder) and README animations
├── input_parameters         # Runtime namelist, no comments (edit before running)
├── input_parameters_with_comments  # Same file, annotated, with all optional groups
└── .gitignore
```

## Decomposition consistency (CPU / multi-CPU / GPU / multi-GPU)

Goal: the same case gives the same answer (to roundoff) at every rank count and layout, on CPU and GPU, for every feature
combination (LES, IBM, Boussinesq, scalars, UAV, inflow/outflow, point particles). This is enforced by the regression suite:
each case is compared at np=1 vs 2, 2x2 and 4x1 (particles also 3), restart vs uninterrupted, and CPU vs GPU.

### Tools

| Tool | Use |
|---|---|
| `DOPAMINE_TRACE_DIR=<dir>` (+ `DOPAMINE_TRACE_STEPS`, `DOPAMINE_TRACE_START`) | `src/debug_trace.f90`: dump U,V,W,P,nu_t,T,C (and the scalar RHS) at named points of each RK stage, assembled by global index; `seam.log` lists every ghost cell that disagrees with the owning rank. Also writes the IBM ghost lists and local `phi`. Off unless the variable is set. |
| `tests/regression/trace_case.sh <case> <exe> [layout_a] [layout_b] [steps] [tol]` | Run a case at two layouts (`2`, `2x2`, ...) and report the first trace point / field / cell that differs. |
| `tests/regression/trace_restart.sh <case> <exe> [layout]` | Same for hot-start restart vs an uninterrupted run. |
| `tests/regression/compare_particles.py` | Final particle state by ID between two runs (layouts or CPU/GPU). |
| `tests/regression/ibm_ghost_counts.py`, `ibm_ghost_diff.py` | IBM ghost-cell bookkeeping between layouts (counts / cell by cell). |

`trace_compare.py` gates on the *interior* of the global array and scales by the field spread (with a floor of 1e-8 of the
field magnitude). The older `compare_monitor.py --field-tol` scales by the field maximum, which hides differences in a
fluctuation on a large offset (Boussinesq `T` around `T_ref`).

### Notes

- Multi-GPU needs `ibm_E` interior cells per rank in x and z (checked at setup, with a message).
- CPU-vs-GPU comparisons use a gfortran CPU build and an nvfortran GPU build; results agree to roundoff.
- Running several 4-rank GPU tests at once (`ctest -j`) can fail spuriously on 2 GPUs; run them sequentially.

## Contributing

License: AGPL-3.0-or-later, see [LICENSE](../LICENSE). Citation information is in the [README](../README.md#citing-this-solver) and [`CITATION.cff`](../CITATION.cff); the methods implemented are listed in [References](References.md).

- If you spot a mistake anywhere in the docs (including this wiki), please open an issue
  on GitHub.
- Bug reports and pull requests are welcome; there is no separate CONTRIBUTING template
  at present, so a descriptive issue or PR is the way to start.
