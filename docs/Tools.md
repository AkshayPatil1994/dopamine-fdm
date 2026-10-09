# Pre- and Post-Processing Tools

[← Home](Home.md)

## GenSDF

`preProcessing/GenSDF/` is a standalone MPI-Fortran signed-distance-field (SDF)
generator used to build the cell-centre SDF the solver's IBM reads via `ibm_sdf_file`
(`&IBM` namelist — see [Input Parameters](Input-Parameters.md#ibm)). The solver
reads only this SDF, so `GenSDF` is the way to produce an IBM geometry.

Build:

```bash
cd preProcessing/GenSDF
make                                  # CPU build (mpif90)
make GPU=1 FC=/path/to/nvfortran      # optional GPU-accelerated build (OpenACC)
```

Configure via `parameters.in` (plain text, positional — comments above each line
describe the field): input geometry (a single OBJ/STL, or a `.list` manifest of
multiple solids for per-object BCs), grid dimensions/origin, whether to read a
non-uniform $z$ grid, whether to use the fast-sweep algorithm (narrow-band, faster) vs.
brute-force, narrow-band width, the solver's wall-normal axis (2=y-vertical,
3=z-vertical), and whether to emit face-staggered SDFs in addition to the cell-centre
one. See `examples/dns_ibm_wavyWall/geo/parameters.in` for a worked configuration
matching that example ([Examples § dns_ibm_wavyWall](Examples.md#dns_ibm_wavywall)).

**Stretched grids.** GenSDF reads the solver's `fields/geometry.out` and `fields/grid.out`
(copy them into `data/`) and treats the wall-normal axis as non-uniform when the
"non-uniform" flag is `.true.`. When the solver's spanwise grid is stretched
(`z_bc_type = 1` with `alpha_grid_z > 0`) it additionally writes `fields/grid_z.out` (same
five-column format as `grid.out`); copy that into `data/` as well and GenSDF picks it up
automatically (with `vertical_axis = 2`), using the per-cell spanwise spacing in the
narrow-band search and in the fast-sweep update. Without a `grid_z.out` the spanwise axis is
assumed uniform, as before. The solver removes a stale `grid_z.out` when z is uniform.

Run:

```bash
./gensdf   # reads parameters.in from the working directory
```

Output is written to `data/sdfp.bin` — the cell-centre SDF, big-endian float64, ready to
be pointed at by the solver's `ibm_sdf_file` (rename/symlink it to `SDF_in`, or set
`ibm_sdf_file` to its path directly) — and, when the geometry manifest defines more than
one solid, `data/sdfp_objid.bin` (per-solid object-ID field, consumed via
`ibm_objid_file` for per-object boundary conditions). When `compute_face_sdf = .true.`
in `parameters.in`, face-staggered SDFs (`sdfu`, `sdfv`, `sdfw`) are written alongside
`sdfp`.

> **Reference (fast-sweep):** Zhao, H., Osher, S. & Fedkiw, R. (2001/2005); Tsai, Y.-H.R.
> (2002) — see [Numerics § IBM](Numerics.md#6-immersed-boundary-method-ibm).
>
> **Reference (`GenSDF`):** Patil, A., Paranjothi, U.C.K. & García-Sánchez, C. (2025).
> *GenSDF: An MPI-Fortran based signed-distance-field generator for computational fluid
> dynamics applications*. SoftwareX 30, 102117.

## `postProcessing/` — the `dopamine_post` library

`postProcessing/` is a pip-installable package, `dopamine-fdm-post` (importable as
`dopamine_post`), rather than a set of independent scripts. Install it once:

```bash
pip install -e postProcessing/
```

then, from anywhere, either script against it directly:

```python
import dopamine_post as dp

pdat = dp.particles.ParticleData.from_case(".")
x, p = pdat.pdf("ax", bins=100, log=True)          # solver-computed acceleration PDF

series = dp.fields.FieldSeries(".")
snap = series.latest()
dp.fields.plot_profile(snap, save="profile.png")
```

or use the console script it installs, `dopamine-post`, which mirrors every module as a
subcommand group (`dopamine-post <module> <command> [args]`; `--help` at any level lists
what's available):

```bash
dopamine-post particles animate --last 40 --n 10 --view 3d --out tracks3d.gif
dopamine-post fields stats --start 60000 --end 120000 --interval 500
dopamine-post uav disk-animate uav_path_file.dat --radius 0.15
```

`postProcessing/examples/` has one short, argparse-free script per module
(`particles_example.py`, `fields_example.py`, ...) demonstrating the library calls
directly — read those instead of `--help` output to see idiomatic library usage.

### Modules

| Module | Covers | Key classes/functions |
|---|---|---|
| `dopamine_post.fields` | Binary field snapshots, time-averaged statistics, DNS comparison, XDMF export | `FieldSnapshot`, `FieldSeries`, `FieldStats`, `plot_profile`, `plot_stats`, `write_field_xmf` |
| `dopamine_post.particles` | Lagrangian point particles (`src/particles.f90`), including the solver-computed per-step acceleration (`ax`/`ay`/`az`, added alongside position/velocity/age) | `ParticleData`, `Geometry`, `plot_tracks`, `animate_tracks`, `animate_cloud`, `write_particles_xmf` |
| `dopamine_post.probes` | Line and slice probe output, XDMF export for slices | `LineProbe`, `SliceProbe` |
| `dopamine_post.ibm_surface` | Per-point IBM surface samples (pressure/viscous force) | `IBMSurface` |
| `dopamine_post.sdf` | Signed-distance-field reading and slice plotting | `SDF` |
| `dopamine_post.uav` | UAV actuator-disk path + scaled quadcopter body animation | `UAVPath` |
| `dopamine_post.rsb` | Reynolds-stress budget statistics | `RSBStats` |
| `dopamine_post.inflow` | SEM/ESEM inflow tooling: half-channel mirroring, donor-plane verification, inflow-optimization state | `mirror_half_channel`, `InflowDonor`/`check_donor`, `InflowOptState`/`read_inflow_opt` |
| `dopamine_post.runlog` | Solver diagnostics parsed from a run's stdout log | `RunLog` |
| `dopamine_post.vortex` | Q-criterion vortex identification and iso-surface movies | `q_criterion`, `q_volume_fraction`, `animate_q_isosurface` |

A shared internal module, `dopamine_post._core`, holds the low-level helpers every other
module builds on (Fortran-namelist parsing, `<prefix>.<step>` snapshot globbing, the
big-endian binary reader, ffmpeg movie encoding, percentile colour clipping, the cubic-
Hermite interpolation that matches `uav_actuator.f90`'s `uav_current_center`, and the
ParaView `.pvd` writer) — not part of the public API, but worth knowing about if you're
extending the library rather than just using it.

### `dopamine_post.particles`

Reads `fields/<fileout>_particles.<step>` (written by `src/particles.f90`, active
whenever `particles_active=1` in `&PARTICLES`). `ParticleData.from_case(case_dir)` loads
every snapshot into id-indexed `(Nt, Np, 3)` arrays (`pos`, `vel`, and `acc` — NaN where a
particle is absent that step). `acc` is the solver's own per-step acceleration when the
snapshot was written by a build with that feature (10-column binary layout: x,y,z,u,v,w,
ax,ay,az,age); for older 7-column snapshots (x,y,z,u,v,w,age, from before this feature
existed) it falls back to a finite difference of the stored velocity between snapshots —
only as good as the snapshot write cadence, so prefer a fresh run when trajectory
acceleration statistics matter.

`ParticleData.write_xmf()` (or the module function `write_particles_xmf`) writes
`paraview/particles.xmf`, pointing directly at the raw binary snapshots via byte-seek
HyperSlabs (no data duplication) — open it alongside `dopamine_post.fields.write_field_xmf`'s
output and use ParaView's shared time toolbar to scrub fields and particles together (both
use the solver step number, not physical time, as the XDMF `Time Value`). Exposes `id`,
`age`, `Velocity`, and (10-column snapshots only) `Acceleration` as point-cloud Attributes.

`Geometry` draws the solid a particle case is seeded over for track/cloud plots: an
analytic wavy wall (`Geometry.analytic`), nothing (`Geometry.flat`), or the `phi=0`
iso-surface of a solver SDF (`Geometry.from_sdf`, needs `scikit-image`).

### `dopamine_post.uav`

`UAVPath.read(uav_path_file)` then `.disk_animation(...)` (a moving-disk `.vtp`/`.pvd`
animation, using the same cubic-Hermite interpolation as the solver's own
`uav_current_center` in `uav_actuator.f90`, so the rendered disk sits exactly where the
actuator-disk force was applied) or `.drone_animation(...)` (a scaled quadcopter body,
`--stl-out`-style, with motion carried by a small `(time, position, rotation)` transform
table applied on the fly by a generated ParaView Programmable Source script — no per-frame
mesh duplication; needs `trimesh`, `manifold3d` optional for a proper boolean union). Both
accept `times_from=` a slice probe's `<base>_times.bin` (exact simulation write times) or
a line probe's plain-text time list, so the rendered path lands on the same frames as a
field/slice/particle series opened alongside it. `drone_animation(tilt=True)` additionally
replays `uav_disk_state`'s auto-tilt low-pass filter (only meaningful with
`uav_tilt_active=1`).

> **Reference (fast-sweep / cubic Hermite):** see [Numerics](Numerics.md) for the UAV
> actuator-disk model.

### `dopamine_post.fields`

`FieldSeries(case_dir).latest()` / `.read(step)` returns a `FieldSnapshot` (cell-centred
U,V,W,P and, when active, C/nu_t/T, ghost layers stripped) for wall-normal profile plots
(`plot_profile`) or arbitrary planes (`FieldSnapshot.plane`). `FieldSeries.time_average(
step_range, interval)` reproduces the former `compute_stats.py` (mean profiles, resolved
Reynolds stresses, resolved+SGS dissipation, centreline symmetry fold) as a `FieldStats`
you can `.save(stats_dir)` or hand to `plot_stats` (with an optional DNS/MKM overlay).
`FieldSeries.profile_vs_dns(dns_dir, x_stations)` reproduces the former
`analyse_channel.py` (multi-station spanwise+temporal-averaged profiles against
Moser–Kim–Mansour DNS and/or a measured wind-tunnel inflow — resolved stresses only; a
wall-modelled/coarse LES carries part of the stress in the SGS model, so a near-wall
deficit vs. DNS is expected). `FieldSeries.write_xmf()` writes `paraview/channel_test.xmf`
via the same zero-copy byte-seek approach as `dopamine_post.particles`.

### `dopamine_post.probes`

`LineProbe`/`SliceProbe` load 1-D line-probe and 2-D slice-probe output (`<base>.bin` +
`<base>_meta.txt`, written by the solver's probe module — see
[Input Parameters § STATISTICS](Input-Parameters.md#statistics));
`SliceProbe.write_xmf(...)` is the former `generate_slice_xmf.py`, reading exact write
times from `<base>_times.bin` when present (`--dt`/`--t0` remain a fallback for when it's
missing — line probes don't write one).

### `dopamine_post.rsb`

`RSBStats.read(case_dir).plot(last=None)` is the former `read_RSBstats.py`: production,
pressure-strain, viscous/turbulent/pressure diffusion, resolved/SGS dissipation, residual,
optionally averaging only the last `N` accumulated samples.

### `dopamine_post.sdf`

`SDF.read(case_dir)` plots the two orthogonal slices (`plot_xz`/`plot_xy`) of a
cell-centre SDF (`SDF_in`) for sanity-checking a `GenSDF` output before a run.

### `dopamine_post.inflow`

`mirror_half_channel` mirrors a half-channel SEM inflow profile into a full-channel one
(handling the sign flip of $V$ and any $v$-linear correlations under the reflection);
`InflowDonor`/`check_donor` verifies a SEM/ESEM inflow donor plane against a reference
profile; `InflowOptState`/`read_inflow_opt` reads/summarizes a SEM inflow-optimization
restart file.

### `dopamine_post.runlog`

`RunLog.read(path).plot(variables=...)` extracts and plots solver diagnostics (mean/max
velocity, divergence, convective/viscous CFL, `dt`) from a run's stdout log.

### `dopamine_post.vortex`

`q_criterion(snap)` returns $Q = \tfrac12(|\Omega|^2-|S|^2)$ on the cell centres of a
`FieldSnapshot` (central differences on the stretched grid); `q_volume_fraction(snap, qval)`
gives the volume fraction with $Q>$ `qval`. `animate_q_isosurface(series, "q.mp4", qval=5)`
(CLI: `dopamine-post vortex animate`) renders the $Q$ iso-surface coloured by $|u|$ beside
mid-span and cross-section $|u|$ slices. The movie needs VTK (`pip install dopamine-fdm-post[vortex]`),
ffmpeg, and an X display (wrap in `xvfb-run` on a headless node).

### `dopamine_post.ibm_surface`

`IBMSurface.read(path)` loads per-point IBM surface samples
(`ibm_surface/surface.<step>.bin`, see
[Input Parameters § IBM](Input-Parameters.md#ibm), `ibm_surface_nsampling`):
position, normal, pressure, and pressure/viscous force per point (summing these
reproduces the drag reported in `ibm_forces.csv`); `.to_vtp(out)` / `IBMSurface.to_pvd(
glob_pattern, out)` convert one file or a whole time series for ParaView.


## Two-phase VOF tools

* **`tests/regression/vof_check.py`** — runs a `tests/regression/vof_*` case at one or two rank counts and checks `vof_diag.dat`:
  volume drift, `C` in [0,1], the shape error of the reversal tests, agreement of selected columns between layouts
  (`--np-b`, `--p-grid`, `--same-cols`, `--same-tol`), final-row limits (`--limit COL:MIN:MAX`), the largest value over the run
  (`--rowmax`), and oscillation frequencies from zero crossings (`--freq COL:THEORY:TOL`). The ctest entries `vof_*` call it.
* **`scripts/vof_figures.py`** — `breaker RUN OUT.png t1 t2 ...` draws liquid fraction and interface from the `vof_snap_*.dat`
  slices (`vof_debug = 1`, `vof_snap_dt > 0`; per-rank files are merged).
* **`vof_diag.dat`** — the named 30-column diagnostic file of the two-fluid solver (columns in
  [Two-Phase VOF § Running](Two-Phase-VOF.md#2-running-and-output)); it is plain text and loads with `numpy.loadtxt(..., comments='#')`.

## dopamine-ESEM (precomputed SEM inflow)

`build/dopamine-ESEM` precomputes the SEM inflow plane of `&INFLOW inflow_type=1` once and writes it as a donor slice that a run with `inflow_type=2` reads, so the per-step eddy cost is paid only once.

```bash
mpirun -np N build/dopamine-ESEM --input=PATH [--out=BASENAME] [--samples=N | --dt=T] \
       [--duration=T] [--ensemble-samples=N] [--ensemble-periods=P]
```

`--input` (required) names the namelist file with the SEM settings (`sem_profile_format` 0 or 1). Pass `--dt` equal to the consuming run's time step (or an integer fraction of it): the default spacing is much finer than a CFL-limited step, and the consuming run would then sample near-uncorrelated points that the pressure projection removes (Reynolds stresses collapse just downstream of the inlet). Only z is decomposed, so the file works with any rank count in the consuming run.

## Showcase animations

The animations in the README are rendered from the validation-run snapshots by the scripts in
`postProcessing/animations/`. They read single planes straight out of the snapshot files
(memory-mapped, so a 4 GiB 512³ snapshot costs a few MiB of I/O per frame), render frames in
parallel with matplotlib, and encode an optimised GIF with `ffmpeg`.

| Script | Case | Content |
|--------|------|---------|
| `animate_tgv.py` | Taylor–Green vortex, Re = 1600, 512³ | vorticity magnitude on a `z = const` plane (central differences from three adjacent planes) |
| `animate_rb.py` | Rayleigh–Bénard convection (Boussinesq) | temperature: vertical slice and near-wall plan view |
| `animate_wavywall.py` | Wavy-wall DNS (ghost-cell IBM) | streamwise velocity, vertical slice with the solid masked from the SDF, and plan view |
| `animate_chan395.py` | Wall-modelled LES, Re_τ = 395 | streamwise velocity slice and near-wall streaks `u'` in plan view |
| `animate_vof_breaker.py` | Plunging breaker, two-fluid VOF (ak = 0.55, 1000:1, 128² cells, t ≤ 0.8 s) | liquid fraction and interface from the `vof_snap_*.dat` debug snapshots (`vof_debug = 1`, `vof_snap_dt`; set `nmonitor` so a snapshot is written every ~0.02 s, since snapshots are written at the monitor interval) |
| `animate_vof_zalesak.py` | Zalesak's slotted disk (`examples/vof_zalesak_disk`, 100² cells) | liquid fraction with the exact rotated shape as a dashed outline |

`anim_common.py` holds the shared code: the `Snapshot` plane reader, the step→time map
(interpolated from the solver's monitor lines in `run.log`, since snapshots carry no time stamp),
and the GIF encoder.

```bash
python postProcessing/animations/animate_tgv.py \
    --case-dir /path/to/tgv/512 --every 3 --out docs/animations/tgv.gif
```

Common options: `--case-dir`, `--prefix` (the `fileout` of the run), `--out`, `--every N`
(use every N-th snapshot), `--fps`. Case-specific: `--z-frac` (TGV plane), `--y-plan` (RB, wavy
wall plan-view height), `--yplus` and `--nstart` (channel). Requires `numpy`, `matplotlib`, and
`ffmpeg` on the `PATH`.

The two VOF scripts read the debug snapshots of a `vof_debug = 1` run (`--case-dir` = the run directory, `--tmax`, `--every`) and
draw the liquid fraction with bicubic interpolation at 256 GIF colours. The breaker is the `vof_break_small` input at N = 128 with
`vof_sigma = 0.072`, `vof_wave_amp = 0.0875`, `vof_snap_dt = 0.02`, `nmonitor = 4` and `sim_end_time = 1` (about 15 min on 4 ranks; the
animation stops at 0.8 s, after the splash); the Zalesak run takes about 30 s on one rank with `vof_snap_dt = 0.02`, `nmonitor = 20`.
