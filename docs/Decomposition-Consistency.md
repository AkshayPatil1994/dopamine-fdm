# Decomposition consistency (CPU / multi-CPU / GPU / multi-GPU)

Goal: the same case gives the same answer (to roundoff) at every rank count and layout, on CPU and GPU.
CPU-vs-GPU at equal layout is already at roundoff; the remaining inconsistencies are between **rank layouts**
and exist in the CPU build too.

## Tools

| Tool | Use |
|---|---|
| `DOPAMINE_TRACE_DIR=<dir>` (+ `DOPAMINE_TRACE_STEPS`, `DOPAMINE_TRACE_START`) | `src/debug_trace.f90`: dump U,V,W,P,nu_t,T,C (and the T RHS) at named points of each RK stage, assembled by global index; `seam.log` lists every ghost cell that disagrees with the owning rank. Off unless the variable is set. |
| `tests/regression/trace_case.sh <case> <exe> [layout_a] [layout_b] [steps] [tol]` | Run a case at two layouts (`2`, `2x2`, ...) and report the first trace point / field / cell that differs. |
| `tests/regression/trace_restart.sh <case> <exe> [layout]` | Same for hot-start restart vs an uninterrupted run. |
| `tests/regression/ibm_ghost_counts.py`, `ibm_ghost_diff.py` | Compare IBM ghost-cell bookkeeping between layouts (counts / cell by cell). |

`trace_compare.py` gates on the *interior* of the global array and scales by the field spread. The existing
`compare_monitor.py --field-tol` scales by the field maximum, which hides differences in a fluctuation on a large
offset (Boussinesq `T` around `T_ref`).

## GPU check

Repeated on the GPU build (2x RTX A6000, nvfortran, CUDA 12.2; run with the matching user-space libcuda on
`LD_LIBRARY_PATH`): the trace tools work unchanged, GPU np=1 vs np=2 shows the same Boussinesq divergence
(same stage, same cell, 1.45e-5), and the GPU np=2 trace equals the CPU np=2 trace to 1e-11 over 210 dumps
(bouss_small). So the layout dependence below is a property of the algorithm, not of the GPU port. Note: 4-rank tests
on 2 GPUs can fail spuriously when several GPU tests run concurrently (`ctest -j`); they pass when run alone.

## Findings (2-step traces, np=1 vs np=2 / 2x2 / 4x1, all regression cases)

| Case | Result |
|---|---|
| tgv, chan, chan_les, sem, rsb, uav, uavpath | interior identical to 1e-11 at every trace point |
| duct | np=2 identical; 2x2 not applicable (`p_col=1` is forced) |
| **bouss** (Boussinesq T) | `T` RHS differs at the two seam planes from RK stage 2 (rel. 1.5e-5); all inputs to the RHS are identical there |
| **sed** (passive scalar C) | identical at np=2, differs at 2x2 (1.7e-6 after 2 steps): same cause as bouss |
| **ibm** | np=2 identical; **2x2: 0.46 relative difference in U at step 1**, at the x seam |
| restart at np>1 | every wall-model case differs at step start in the y-ghost row of the plane next to the seam; DNS `chan_small` and all np=1 restarts are clean |

### 1. Scalar/temperature advection: first-order fallback at the local array edge (bouss, sed)

`compute_rhs_scalar_core` (`scalar_transport.f90`) uses a van Leer MUSCL reconstruction that needs the cell two
upwind of the face. Where that cell is outside the local array (`k == nzg-1`, `k == 2`, and the x equivalents) it
falls back to first order. With a 1-cell halo this happens at every rank seam, and also at the periodic domain edge
at np=1 (where the wrapped value exists but is not used). So the scheme order at the seam plane depends on the layout.

Evidence: the difference first appears in the stage-2 T RHS, only on the seam planes (141 cells on the last plane of
rank 0, 9 on the first plane of rank 1: the asymmetry follows the sign of W). Forcing the same fallback on the same two
global planes at np=1 removed the T, RHS and V differences (< 1e-11).

The fix must give the reconstruction its second upwind cell on all layouts (a second halo plane for scalars, or an
extra far-plane exchange), which also corrects the np=1 periodic edge. That changes np=1 results slightly.

### 2. IBM ghost-cell bookkeeping depends on the layout

`ibm_small`, global ghost cells (U, V, W) and dropped cells (image clips solid, U,V,W,CC):

| layout | kept U | kept V | kept W | dropped (U,V,W,CC) |
|---|---|---|---|---|
| 1 | 340 | 368 | 396 | 56 26 12 1 |
| 2 (z split) | 340 | 368 | 396 | 56 26 12 1 |
| 1x4 (z split) | 308 | 335 | 332 | 88 59 72 10 |
| 2x2, 4x1 (x split) | **470** | 368 | 396 | 56 26 12 1 |

- z split into 4 slabs (8 cells each) loses ghost cells: the image point is at least ~2 cells from the ghost cell,
  outside the 1-cell halo, so its stencil is rejected and that wall cell gets no boundary condition.
- x split creates 130 *spurious* U ghost cells (none in np=1), all on the last face of one rank next to the seam
  (`ibm_ghost_diff.py`). Not yet root-caused; this is the 0.46 difference at step 1.

### 3. Restart at np>1: y-ghost row of the seam plane

The trace shows the state read from the restart file differs from the uninterrupted state in `U(i, nyg, k)` (the
top-wall ghost row) on the plane next to the seam, only when the wall model is on. `seam.log` shows the same cells
disagree with the owning rank during normal runs: the y-ghost rows of the *seam ghost planes* are computed from
wall-model slip lengths whose ghost planes are a local copy (`alpha(:,:,1)=alpha(:,:,2)`, `wallmodel.f90`), not the
neighbour's values. The restart writer then lets the ghost-plane value overwrite the owner's on the shared plane.
The interior solution does not read these cells, which is why the interior agrees at np=1 vs np>1 but a restart does not.

### Not reproduced

UAV actuator disk (static across the seam, and moving) agrees to 1e-11 over 6 steps at np=2 and 4x1; the documented
UAV np-dependence did not show. Recycled precursor inflow (`inflow_type=2`) has no test case and was not traced.

## Plan status

- Phase 0 (diagnostics): done. The tools above, three localized causes, IBM ghost-count tests (`WILL_FAIL` until fixed).
- Next: fix (1) second upwind cell for scalar MUSCL; fix (3) writer takes owned cells only and the wall-model ghost-plane
  rows are made consistent; then the IBM stencil gather (2); then tighten the test bounds marked `KNOWN ISSUE`.
