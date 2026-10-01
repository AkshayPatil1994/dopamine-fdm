# Pre- and Post-Processing Tools

← [[Home|Home]]

## GenSDF

`preProcessing/GenSDF/` is a standalone MPI-Fortran signed-distance-field (SDF)
generator used to build the cell-centre SDF the solver's IBM reads via `ibm_sdf_file`
(`&IBM` namelist — see [[Input Parameters|Input-Parameters#ibm-optional]]). The solver
no longer accepts a face-point mask input (`Umask_in`) — `GenSDF` is the only supported
way to produce an IBM geometry.

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
matching that example ([[Examples § dns_ibm_wavyWall|Examples#dns_ibm_wavywall]]).

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
> (2002) — see [[Numerics § IBM|Numerics#6-immersed-boundary-method-ibm]].

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

> **Reference (fast-sweep / cubic Hermite):** see [[Numerics|Numerics]] for the UAV
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
[[Input Parameters § STATISTICS|Input-Parameters#statistics-optional--omit-to-disable]]);
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
[[Input Parameters § IBM|Input-Parameters#ibm-optional]], `ibm_surface_nsampling`):
position, normal, pressure, and pressure/viscous force per point (summing these
reproduces the drag reported in `ibm_forces.csv`); `.to_vtp(out)` / `IBMSurface.to_pvd(
glob_pattern, out)` convert one file or a whole time series for ParaView.
