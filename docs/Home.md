# dopamine-fdm Wiki

`dopamine-fdm` is a parallel, finite-difference solver for the 3-D incompressible
Navier–Stokes equations, targeting turbulent channel and open-channel flows. It supports
rough-wall immersed boundary methods (IBM), equilibrium wall models, sub-grid scale (SGS)
turbulence modelling, suspended-sediment scalar transport, and an exact Reynolds stress
budget analysis module.

This wiki expands on the top-level [README](../README.md) with a page per topic. Start
here, then jump to whichever page matches what you're doing.

## Pages

| Page | Use it for |
|------|------------|
| [Installation & Running](Installation.md) | Dependencies, CPU/GPU build instructions, launching a run, output layout |
| [Input Parameters Reference](Input-Parameters.md) | Every namelist and parameter in `input_parameters`, grouped by section |
| [Numerics & Governing Equations](Numerics.md) | The discretisation, time integration, Poisson solver, SGS model, IBM, wall models, scalar transport, Reynolds stress budget, inflow/outflow BCs |
| [Examples](Examples.md) | Walkthroughs of the three bundled example cases |
| [Pre- and Post-Processing Tools](Tools.md) | `GenSDF` mesh/SDF generator and the Python scripts in `postProcessing/` |
| [Contributing & Citing](Contributing.md) | License, how to cite, acknowledgments, opening issues |

## At a glance

- **Discretisation** — staggered MAC grid, second-order central differences, low-storage
  RK3 time integration, fractional-step projection with a spectral (local FFTW3 +
  2decomp&fft pencil transposes) Poisson solver.
- **Parallelism** — [2decomp&fft](https://github.com/2decomp-fft/2decomp-fft) 2-D pencil
  MPI decomposition (auto-prefers a cheap 1-D z-slab split when the rank count allows);
  optional single-GPU (OpenACC + cuFFT + cuSPARSE) build for `nprocs=1`.
- **Physics** — DNS or Vreman-model LES, flat-wall/IBM log-law wall models, rough-wall
  ghost-cell IBM, oscillatory/pulsatile forcing, suspended-sediment transport, periodic or
  inflow/outflow (synthetic-eddy-method) streamwise boundary conditions.
- **Diagnostics** — full Pope §7.4 Reynolds stress budget, line/slice probes, per-stage
  profiler summary at shutdown.

See the repository layout in the [README](../README.md#repository-layout) for where each
piece of source lives.
