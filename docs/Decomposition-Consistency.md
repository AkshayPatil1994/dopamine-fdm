# Decomposition consistency (CPU / multi-CPU / GPU / multi-GPU)

Goal: the same case gives the same answer (to roundoff) at every rank count and layout, on CPU and GPU, for every feature
combination (LES, IBM, Boussinesq, scalars, UAV, inflow/outflow, point particles). This is enforced by the regression suite:
each case is compared at np=1 vs 2, 2x2 and 4x1 (particles also 3), restart vs uninterrupted, and CPU vs GPU.

## Tools

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

## What was wrong and how it was fixed

1. **Scalar/temperature MUSCL fell back to first order at the local array edge** (every rank seam, and the periodic edge at
   np=1). The far-upwind cell now comes from a padded copy (`Cpad`, upstream's fix, adapted so the exchange also works with
   device-resident data). Boussinesq `T` and passive `C` are identical across layouts.
2. **Wall rows of the seam ghost planes were computed from the rank's own wall-model coefficients**, so they disagreed with the
   owning rank. `apply_boundary_conditions` now ends with a pure-copy halo exchange and periodic wrap. Restart at np>1 reproduces
   an uninterrupted run to 1e-12 for every case (the loosened bound was removed), and anything sampling next to a wall and a seam
   is layout independent.
3. **IBM image-point stencils reached past the one ghost plane.** Ghost cells whose image lay beyond it were dropped or mislocated,
   and on x-split layouts the U-face neighbour test read one plane past the array (spurious ghost cells). The IBM now works on an
   extended halo of `ibm_E` extra x/z planes (`halo_pad.f90`; depth from the image distance, checked against the slab width):
   padded SDF, mask and axes for the ghost lists, and extended copies of U, V, W (and the scalar) for the image interpolation, the
   IBM wall model and the force diagnostics. The ghost-cell counts and the solution are identical for np=2, 3, 4, 1x4, 2x2, 4x1.
4. **Point particles.** The tracer RK3 and the IBM collision normal sample the fields at displaced positions that near a seam
   lie beyond the ghost cell; they now use padded copies. Seeding and reinjection use an ID-keyed hash instead of `Random_Number`
   (whose stream differs between compilers and per rank); reinjected IDs no longer collide with the seed IDs.
   Brownian motion and the SGS Langevin closure still use per-rank random streams and are only statistically layout independent.

Also fixed while merging the particle branch: a restart no longer re-applies the wall BCs before the first step (wall-model slip
lengths are not in the restart file), the UAV tilt model uses its own `uav_grav` (the global `grav` now defaults to 0), and the
`cfl_accel` monitor column is parsed by the parity script.

## Notes

- Multi-GPU needs `ibm_E` interior cells per rank in x and z (checked at setup, with a message).
- CPU-vs-GPU comparisons use a gfortran CPU build and an nvfortran GPU build; results agree to roundoff.
- Running several 4-rank GPU tests at once (`ctest -j`) can fail spuriously on 2 GPUs; run them sequentially.
