# Pre-Processing: GenSDF

← [Home](Home.md)

`preProcessing/GenSDF/` is a standalone MPI-Fortran signed-distance-field (SDF) generator. It turns
a surface mesh (STL or OBJ) into the cell-centre SDF that the solver's IBM reads via `ibm_sdf_file`
(`&IBM` namelist — see [Input Parameters](Input-Parameters.md#ibm)), optionally with the
face-staggered SDFs and a per-solid object-ID field. It is the only supported way to produce an IBM
geometry. (Post-processing scripts live on the [Tools](Tools.md) page.)

## 1. Build

Requires an MPI Fortran compiler (`mpif90`, gfortran ≥ 9 or similar). No other libraries.

```bash
cd preProcessing/GenSDF
make            # CPU build -> ./gensdf
make clean      # remove objects/modules;  make allclean also removes ./gensdf
```

### GPU build (OpenACC, NVIDIA)

The GPU build offloads the narrow-band distance search and the fast-sweep fill; the mesh analysis and the
solid/fluid classification stay on the CPU. It needs the NVIDIA HPC SDK (`nvfortran`, tested with 23.3) and its
bundled OpenMPI. Put both on `PATH` — the system `mpif90` cannot find `nvfortran`, which is the usual cause of
"unable to find the specified compiler nvfortran":

```bash
NV=/opt/nvidia/hpc_sdk/Linux_x86_64/23.3          # adjust to your install
export PATH=$NV/compilers/bin:$NV/comm_libs/openmpi4/openmpi-4.0.5/bin:$PATH
cd preProcessing/GenSDF
make allclean                                     # never mix CPU and GPU objects
make GPU=1 FC=mpif90                              # -> ./gensdf (OpenACC)
```

Notes:

- `GPU=1` compiles with `-acc -gpu=cc86,rdc -DGPU_SDF`. **`cc86` is an Ampere card (RTX 30xx/A40).** For another
  GPU edit the two `-gpu=cc86,rdc` occurrences in the `Makefile` (e.g. `cc70` V100, `cc80` A100, `cc89` RTX 40xx,
  `cc90` H100). The binary must be built with the same SDK whose runtime you launch it with.
- Run it with the SDK's `mpirun` (same `PATH`), one MPI rank per GPU: `mpirun -np <ngpu> ./gensdf`. Ranks are
  bound to GPUs by node-local rank modulo the GPU count, and the run stops with an error if no GPU is visible.
- Use `make allclean` whenever you switch between CPU and GPU builds.
- A hash cell holding so many triangles that 27 cells exceed `max_cands_gpu` (32768, `spatial_hash.f90`) stops
  the GPU build with an error message; raise the constant or coarsen the mesh locally.
- The GPU and CPU builds give the same fields to round-off (~1e-14 in distance, identical signs and object IDs).
- Debug build (CPU): `make FFLAGS="-O0 -g -fcheck=all -fbacktrace -ffree-line-length-none -cpp"`.

## 2. Prepare a case directory

Work in a directory that contains `parameters.in` and a `data/` folder:

```
case/
├── parameters.in            # control file (section 3)
└── data/
    ├── geometry.out         # copy of the solver's fields/geometry.out (grid size and lengths)
    ├── grid.out             # only if the wall-normal grid is stretched (see "Stretched grids")
    ├── grid_z.out           # only if the spanwise grid is stretched (vertical_axis = 2)
    ├── objects.list         # or point parameters.in straight at one .stl/.obj
    └── body.stl, ...        # the geometry
```

`geometry.out` and `grid*.out` are written by the solver (`fields/`); generate the grid first (a short solver run
or the grid-only mode), copy them into `data/`, and make sure the geometry is in the **same coordinates and
units as the solver domain** (the grid origin is the `r0` line of `parameters.in`).

**Geometry input.** A single binary or ASCII `.stl`, a Wavefront `.obj` (any polygon is fan-triangulated; vertex
normals are optional), or a `.list` manifest with one path per line (`#` comments and blank lines skipped). The
order of the lines is the solid ID, used for per-object boundary conditions (at most 15 solids). Append ` ground`
to a line for an *open* surface — terrain, bathymetry, or a building shell without a floor: everything below its
topmost surface along the vertical axis is then solid, down to the domain floor.

```
data/building_a.stl
data/building_b.stl
data/terrain.stl ground
```

## 3. `parameters.in`

Plain text, **positional**: each non-comment line is one entry, lines starting with `#` and blank lines are skipped.

| # | Entry | Meaning |
|---|-------|---------|
| 1 | `'data/objects.list'` | geometry file or manifest (quoted) |
| 2 | `1d10  5  40` | far-field sentinel value, `buffer_points` (cells of margin around the geometry's bounding box), progress-bar width |
| 3 | `nx ny nz` | interior cell-centre counts, as in the solver's `geometry.out` |
| 4 | `0.0 0.0 0.0` | grid origin |
| 5 | `.false.` | read a non-uniform wall-normal grid from `data/grid.out` |
| 6 | `.true.` | fast sweep (narrow band + Eikonal fill, recommended) vs. brute force (distances only near the surface, no fill) |
| 7 | `5` | narrow-band width in cells (kept for compatibility; the search reach is 6 cells) |
| 8 | `3` | solver wall-normal axis: `2` = y vertical, `3` = z vertical |
| 9 | `.false.` | also write the face-staggered SDFs `sdfu`, `sdfv`, `sdfw` |
| 10 | `1` | *optional* `sign_method`: `0` three-axis parity, `1` signed rays with winding-number fallback (default), `2` winding number everywhere (slow check) |
| 11 | `1` | *optional* `cap_holes`: close simple planar holes with a lid (`1`, default) or leave the mesh as given (`0`) |

Lines 10–11 may be omitted; the defaults apply.

## 4. Run

```bash
cd case
mpirun -np 8 /path/to/preProcessing/GenSDF/gensdf        # CPU: any number of ranks >= 1
# GPU build: mpirun -np <number of GPUs> ./gensdf  (PATH as in section 1)
```

Ranks split the domain in x-slabs; the field, signs and object IDs do not depend on the rank count (distances
agree to ~1e-6 absolute at the very centre of a body, everything else to round-off).

Output in `data/` (big-endian float64, padded by one ghost layer in x and the wall-normal direction as the
solver's IBM reader expects):

| File | Content |
|------|---------|
| `sdfp.bin` | cell-centre SDF, negative inside the solid — point `ibm_sdf_file` at it (or rename to `SDF_in`) |
| `sdfp_objid.bin` | solid ID per cell near the surface (0 = fluid); for `ibm_objid_file` |
| `sdfu.bin`, `sdfv.bin`, `sdfw.bin` | face-staggered SDFs, only if entry 9 is `.true.` |
| `open_edges.obj` | boundary edges of the input surface (holes, cracks), only if the mesh has any |

Use `postProcessing` ([Tools](Tools.md)) to inspect `sdfp.bin` before a run.

## 5. Damaged and non-watertight geometry

Real surface meshes (CAD exports, scans, STL from a mesher) are rarely clean, so GenSDF does not assume a closed,
consistently oriented shell. At start-up it prints a **mesh check** per solid: vertices are welded, degenerate and
duplicated triangles are dropped, boundary edges (holes, cracks, missing faces) are grouped into loops with their
location, size and bounding box, non-manifold and inconsistently oriented edges are counted, and triangle
orientation is repaired. Read it: it tells you whether the surface you meant as a closed body is closed, and the
holes are in `open_edges.obj` (ParaView/MeshLab) to look at.

Solid/fluid is then decided in two tiers (`sign_method = 1`):

1. **Signed rays.** Six half-line rays (±x, ±y, ±z) from every grid point; the orientation-signed crossings along a
   ray are the point's winding number (0 outside, 1 inside, 2 where shells overlap or a patch is doubled) as long as
   the ray meets no defect, so each clean ray is an independent witness even if the mesh is open elsewhere. A ray
   that escapes an open shell without meeting a triangle abstains. A point is accepted when at least three rays vote
   and agree.
2. **Generalized winding number** (Jacobson et al. 2013; BVH with a dipole far field) for every point the rays cannot
   settle. It varies smoothly across a hole, so thresholding at ½ stays correct until a hole subtends more than half
   the sky.

`sign_method = 0` is the original three-axis parity vote: it is silently wrong behind any hole and cancels
completely for doubled or overlapping shells, so use it only for known-closed meshes. `sign_method = 2` evaluates the
winding number everywhere and reports how often it disagrees with a unanimous ray decision.

**Hole lids (`cap_holes = 1`).** A simple, nearly planar boundary loop — a missing floor or wall, a round opening, a
cut end — is closed with an ear-clipped lid before classification. Ray tests become exact and the distance field
stays consistent across the former opening (distance is measured to the lid, not to the rim), which the immersed
boundary needs. Loops that are not planar, not simple, or lie within two cells of another loop (a seam or crack,
where independent lids would cross) stay open for the winding number. The mesh check labels each loop `[lid]`,
`[seam]` or `[left open]`. A large, non-planar hole left open is genuinely ambiguous — the intended volume can only
be guessed — so close such openings in the mesh.

**Ground sheets.** `ground` solids are solid everywhere below their topmost surface, independent of holes.

## 6. Stretched grids and other details

**Stretched grids.** GenSDF reads the solver's `geometry.out` and `grid.out` (copy them into `data/`) and treats the
wall-normal axis as non-uniform when entry 5 is `.true.`. When the solver's spanwise grid is stretched
(`z_bc_type = 1` with `alpha_grid_z > 0`) it also writes `grid_z.out` (same five-column format); copy that in as well
and GenSDF uses it automatically with `vertical_axis = 2`. Without a `grid_z.out` the spanwise axis is assumed
uniform. The bounding-box buffer is `buffer_points` *local* cells.

**Resolution.** A solid narrower than ~2 cells cannot be represented; GenSDF warns per solid.

**Accuracy.** Distances are exact within 6 cells of the surface (the IBM needs ≥ 4); beyond that the fast sweep fills
the field to first-order accuracy (≤ ~0.5 cell far away). The sign is independent of the distance accuracy.

**Regression test.** After changing GenSDF run, from `preProcessing/GenSDF`:

```bash
make && python3 tests/test_gensdf.py            # needs numpy; add --np 1 4 for rank checks, --quick to skip slow cases
python3 tests/test_gensdf.py --exe ./gensdf     # GPU build: run with the SDK PATH set
```

It builds damaged test meshes (holes, missing faces, cracks, flipped/duplicated/overlapping triangles, deleted
triangles, terrain) and checks the masks, the OBJ/STL readers, `vertical_axis = 2`, object IDs and rank independence.

## 7. Worked examples

`preProcessing/GenSDF/examples/case05_sphere_DNS_Re200` and `case06_sphere_LES_Re3700` generate a sphere STL and a
matching `parameters.in`; run `bash run.sh [N_PROCS]` in either directory after building. For the wavy-wall case see
`examples/dns_ibm_wavyWall/geo/parameters.in` ([Examples § dns_ibm_wavyWall](Examples.md#dns_ibm_wavywall)).

> **Reference (fast-sweep):** Zhao, H., Osher, S. & Fedkiw, R. (2001/2005); Tsai, Y.-H.R. (2002) — see
> [Numerics § IBM](Numerics.md#6-immersed-boundary-method-ibm).
