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
│   ├── ibm.f90              # Ghost-cell IBM module
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
│   ├── finalization.f90
│   └── main.f90             # Entry point
├── postProcessing/
│   ├── generateXMF.py       # Write XDMF metadata for ParaView
│   ├── generate_UAVpath.py  # Render a uav_path_file as a moving-disk ParaView animation (see Tools page for the rest)
│   └── animations/          # Showcase GIF scripts (see Animations page)
├── preProcessing/
│   └── GenSDF/              # Signed-distance field generator for IBM
├── input_parameters         # Runtime namelist (edit before running)
└── .gitignore
```
