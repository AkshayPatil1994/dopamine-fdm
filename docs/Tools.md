# Pre- and Post-Processing Tools

← [Home](Home.md)

## GenSDF

`preProcessing/GenSDF/` is a standalone MPI-Fortran signed-distance-field (SDF)
generator used to build the cell-centre SDF the solver's IBM reads via `ibm_sdf_file`
(`&IBM` namelist — see [Input Parameters](Input-Parameters.md#ibm-optional)). The solver
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
matching that example ([Examples § dns_ibm_wavyWall](Examples.md#dns_ibm_wavywall)).

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

## `postProcessing/` scripts

All scripts are run with `python3 <script>.py [args]`, typically from the case directory
(the one containing `input_parameters` and `fields/`). Most that read binary snapshots
build on `snapshot_io.py`.

### `snapshot_io.py`

Reader **library** (not a standalone tool) for the solver's binary field-snapshot
format — parses `input_parameters`, lists available snapshots, and reads a snapshot into
NumPy arrays with ghost layers stripped. Imported by `plot_snapshot.py`, `plot_sdf.py`,
`analyse_channel.py`, etc.:

```python
from snapshot_io import parse_input_parameters, list_snapshots, read_snapshot
```

See the binary layout in [Installation § Output files](Installation.md#output-files).

### `generateXMF.py`

Writes XDMF metadata (`paraview/channel_test.xmf`) so ParaView can open the raw binary
snapshots directly via byte-seek HyperSlab selections — no field data is duplicated.
Auto-detects the snapshot filename prefix and available steps from `fields/`. Run it
**from the case directory** with no arguments:

```bash
python3 postProcessing/generateXMF.py
```

> The `--case`/`--nx`/`--ny`/`--nz` flags shown in the top-level README are stale — the
> current script takes no CLI arguments and auto-detects everything from `fields/`.

### `generate_slice_xmf.py`

XDMF time-series generator for 2-D slice-probe output (`<base>.bin` +
`<base>_meta.txt`, written by the solver's slice-probe module). Writes little-endian
coordinate arrays and a `.xmf` readable by ParaView ≥ 5 / VisIt ≥ 3:

```bash
python3 postProcessing/generate_slice_xmf.py <base>_meta.txt [<base2>_meta.txt ...] [--dt DT] [--t0 T0]
```

### `plot_snapshot.py`

Wall-normal profile plots (x-z averaged) from one or more field snapshots. Reads
`input_parameters` to determine active physics (`sgs_model`, `sediment_flag`):

```bash
python postProcessing/plot_snapshot.py             # latest snapshot
python postProcessing/plot_snapshot.py --step 500  # a specific step
python postProcessing/plot_snapshot.py --all       # time-averaged over all snapshots
python postProcessing/plot_snapshot.py --fields run2/fields --params run2/input_parameters
```

### `compute_stats.py`

Time- and plane-averaged channel-flow statistics (mean profiles, Reynolds stresses)
computed by averaging a window of field snapshots:

```bash
python postProcessing/compute_stats.py --avg_start START --avg_end END --interval N --prefix NAME \
    [--fields_dir fields] [--stats_dir stats] [--nu NU] [--no_sgs] [--out_fig plots/stats.png]
```

### `analyse_channel.py`

Wall-normal profiles from a channel run compared against Moser–Kim–Mansour DNS
reference data and/or a measured wind-tunnel inflow profile. Supports multiple
streamwise stations (spanwise+temporal averaging only, never streamwise), useful for
watching an SEM inflow develop with fetch:

```bash
python3 postProcessing/analyse_channel.py --x 2 4 6 8 \
    --dns-dir /path/to/validation/semChannel \
    --step-start 900 --re-tau 392.24 --output profiles.png
```

> Reported stresses are *resolved* only — a wall-modelled/coarse LES carries part of the
> stress in the SGS model, so a near-wall deficit vs. DNS is expected.

### `read_RSBstats.py`

Reads and plots the Reynolds stress budget (RSB) output described in
[Input Parameters § STATISTICS](Input-Parameters.md#statistics-optional--omit-to-disable)
(production, pressure-strain, viscous/turbulent/pressure diffusion, resolved/SGS
dissipation, residual). Run from the case directory:

```bash
python3 postProcessing/read_RSBstats.py [--last N]
```

`--last N` averages only the last `N` accumulated samples (default: all).

### `plot_log.py`

Extracts and plots solver diagnostics (mean/max velocity, divergence, convective/viscous
CFL, `dt`) from a run's stdout log:

```bash
python postProcessing/plot_log.py [run.log] [-v Umean,Umax,CFLc,CFLv,dt]
```

### `plot_sdf.py`

Loads and plots two orthogonal slices (an x-z plane and an x-y plane) of a cell-centre
SDF (`SDF_in`), useful for sanity-checking a `GenSDF` output before a run:

```bash
python postProcessing/plot_sdf.py [--sdf SDF_in] [--params input_parameters] [--y-slice Y] [--z-slice Z]
```

### `load_ibm_surface.py`

Reader/converter for the per-point IBM surface samples written to
`ibm_surface/surface.<step>.bin` (see
[Input Parameters § IBM](Input-Parameters.md#ibm-optional), `ibm_surface_nsampling`):
surface position, normal, pressure, and pressure/viscous force per point (summing these
reproduces the drag reported in `ibm_forces.csv`).

```python
from load_ibm_surface import read_ibm_surface
d = read_ibm_surface('ibm_surface/surface.00010000.bin')
```

```bash
python load_ibm_surface.py ibm_surface/surface.00010000.bin        # convert one file
python load_ibm_surface.py 'ibm_surface/surface.*.bin'             # convert a whole glob to ParaView .vtp
```

### `load_line_probes.py`

Loads 1-D line-probe output (`<base>.bin` + `<base>_meta.txt`) into a `(nsnaps, ncomp,
npts)` array:

```python
from load_line_probes import load_probe
probe = load_probe('mysim_line_meta.txt')
u = probe['data'][:, probe['comps'].index('U'), :]
```

```bash
python3 postProcessing/load_line_probes.py <meta.txt> [<meta2.txt> ...] [--snap S]
```

### `mirror_half_channel.py`

Mirrors a half-channel SEM inflow profile (one no-slip wall + symmetry plane at the
centreline, `y ∈ [0,h]`) into a full-channel profile (`y ∈ [0,2h]`, no-slip both walls),
for use as `inflow_profile_file` (see
[Input Parameters § INFLOW](Input-Parameters.md#inflow-optional--used-only-when-x_bc_type--1)).
Handles the sign flip of $V$ (and any $v$-linear correlation, `uv`/`vw`) under the
reflection:

```bash
python3 postProcessing/mirror_half_channel.py half_channel.csv full_channel.csv [--h H]
```
