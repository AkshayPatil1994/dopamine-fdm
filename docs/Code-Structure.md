# Code Structure

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
│   ├── vof_curv.f90         # Height-function curvature for surface tension
│   ├── mom_recon.f90        # Momentum face-value kernels (central, Koren, QUICK, WENO3-Z/WENO5-Z), pure functions
│   ├── vof_twofluid.f90     # Two-fluid time step: consistent momentum transport, forces (viscous, gravity, surface tension), projections
│   ├── finalization.f90
│   └── main.f90             # Entry point
├── tests/
│   ├── regression/          # Small deterministic cases: np / layout parity, restart, VOF cases (vof_check.py)
│   └── vof/                 # Standalone VOF unit tests: PLIC geometry, prescribed-velocity advection, momentum kernels
├── postProcessing/
│   ├── pyproject.toml       # dopamine-fdm-post: pip install -e postProcessing/
│   ├── dopamine_post/       # Post-processing library (see Tools page) + `dopamine-post` CLI
│   ├── examples/            # Short demonstration scripts, one per dopamine_post module
│   └── animations/          # Showcase GIF scripts (see Animations page)
├── preProcessing/
│   └── GenSDF/              # Signed-distance field generator for IBM
├── input_parameters         # Runtime namelist (edit before running)
└── .gitignore
```
