# Installation & Running

← [[Home|Home]]

## Dependencies

| Library | Minimum version | Purpose |
|---------|----------------|---------|
| gfortran or ifort | GFortran ≥ 9 / Intel ≥ 2019 | Fortran compiler |
| MPI | any standard MPI-3 | Domain decomposition |
| FFTW3 | 3.3 (serial double-precision) | Local single-rank transforms/DCT, and the transform engine inside 2decomp&fft |
| [2decomp&fft](https://github.com/2decomp-fft/2decomp-fft) | `v2.1.0` | 2-D pencil domain decomposition and inter-rank transposes for the MPI-parallel pressure Poisson solve — see [[Numerics § MPI parallelism|Numerics#10-mpi-parallelism]] |
| LAPACK / BLAS | any | Linear algebra |
| CMake | 3.20 | Build system |
| git | any | Required at first configure to fetch 2decomp&fft (see below) |
| NVHPC SDK *(optional, GPU build only)* | 23.3+ | `nvfortran` + OpenACC + cuFFT + cuSPARSE |

On Debian/Ubuntu the CPU-build dependencies (excluding 2decomp&fft, fetched
automatically) can be installed with:

```bash
sudo apt install gfortran libopenmpi-dev libfftw3-dev liblapack-dev libblas-dev cmake git
```

## Building (CPU)

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

The executable `dopamine` is placed in `build/`.

> **2decomp&fft is fetched and built automatically** on first configure: CMake clones
> `2decomp-fft/2decomp-fft` (pinned to `v2.1.0`) into `build/_deps/2decomp_fft-src`,
> configures and builds it as an independent CMake invocation (`BUILD_TARGET=mpi`,
> `FFT_Choice=fftw`), and installs it under `build/_deps/2decomp_fft-install`. This needs
> `git` and network access the first time you configure the CMake build for a given build
> directory; once installed there, subsequent `cmake`/`cmake --build` calls reuse it
> without re-fetching. If you're on an offline/air-gapped build host, pre-populate
> `build/_deps/2decomp_fft-src` with a clone of the pinned tag before running `cmake -S`.

> **Note**: The solver uses `-fconvert=big-endian` globally so that field snapshots are
> big-endian. Reynolds stress budget output uses `CONVERT='little_endian'` on the file
> `OPEN` to produce little-endian float64 (consistent with the coordinate `.bin` files
> used by `generateXMF.py`).

## Building (GPU: OpenACC + cuFFT + cuSPARSE, single- or multi-GPU)

An optional `ENABLE_GPU` CMake target offloads the RHS/SGS/boundary-condition/projection
kernels and the pressure Poisson solve (cuFFT transforms + cuSPARSE batched tridiagonal
solve) to NVIDIA GPUs via OpenACC, one MPI rank per GPU. The CPU build above is untouched
and remains the default; this is a separate, opt-in build.

**Requirements**:
- [NVIDIA HPC SDK](https://developer.nvidia.com/hpc-sdk) 23.3 or later (provides
  `nvfortran`, its bundled **CUDA-aware** OpenMPI, and the `cufft`/`cusparse` device
  libraries) — e.g. installed under `/opt/nvidia/hpc_sdk`. Use its bundled `mpirun`/`mpif90`.
- NVIDIA GPUs whose driver supports the toolkit: the user-space `libcuda` must match the
  kernel driver (a stale `libcuda` crashes at start-up in `cuModuleLoad`). If the driver
  is older than the toolkit NVHPC defaults to, pass `-DGPU_CUDA_VERSION=<X.Y>` (nvfortran
  `-gpu=cudaX.Y`) with a matching toolkit installed.
- Compute capability defaults to `86` (Ampere, e.g. RTX A6000 / RTX 30-series); override
  with `-DCMAKE_CUDA_ARCHITECTURES=<cc>` (no dot) for other GPUs.

**How multi-GPU works**: the vendored 2decomp&fft is built in its GPU mode
(`-DDECOMP2D_GPU=ON`, default with `ENABLE_GPU`), so the pencil transposes of the
Poisson solve, and the velocity/scalar halo exchanges, pass *device* buffers straight to
MPI (GPU-aware MPI); each rank is bound to its own GPU by node-local rank
(`src/gpu_device.f90`). Launch with one rank per GPU.

**Scope — read before using**:
- **Multi-GPU (`nprocs>1`)**: supported for every GPU-supported BC combination, pencil-decomposed
  like the CPU build (`p_row`/`p_col`): periodic/wall `y`, periodic (FFT) or inflow/outflow
  (`x_bc_type=1`, DCT-IV) `x`, periodic `z`, and the 4-wall duct (`y_bc_type=1`, `z_bc_type=1`,
  `x_bc_type=0`; `p_col=1` is forced, `x` is split). Wall models, SGS, scalars etc. work as on
  one GPU, including the ghost-cell IBM and point particles (results match the CPU build and every rank layout; see
  [[Decomposition Consistency|Decomposition-Consistency]]). The IBM needs a few interior cells per rank in x and z (checked at start-up).
- **CPU-only**: a spanwise wall alone (`z_bc_type=1`, `y_bc_type=0`), or a duct with
  `x_bc_type=1`. A runtime guard `Stop`s immediately on unsupported combinations.

Put the NVHPC SDK's `nvfortran` and bundled MPI on your `PATH`/`LD_LIBRARY_PATH` first,
e.g.:

```bash
export PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/23.3/comm_libs/openmpi4/openmpi-4.0.5/bin:/opt/nvidia/hpc_sdk/Linux_x86_64/23.3/compilers/bin:$PATH
export LD_LIBRARY_PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/23.3/comm_libs/openmpi4/openmpi-4.0.5/lib:/opt/nvidia/hpc_sdk/Linux_x86_64/23.3/compilers/lib:$LD_LIBRARY_PATH
```

Then configure and build with `nvfortran` and `-DENABLE_GPU=ON`:

```bash
cmake -S . -B build_gpu -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_Fortran_COMPILER=nvfortran -DENABLE_GPU=ON \
      -DMPI_Fortran_COMPILER=$NVHPC/comm_libs/openmpi4/bin/mpif90
cmake --build build_gpu -j$(nproc)
```

The executable `dopamine` is placed in `build_gpu/`. `-Minfo=accel` output during the
build lists every OpenACC region the compiler offloaded — useful to confirm a kernel you
touched is still being generated for the GPU.

## Running

1. Edit `input_parameters` to set domain size, grid type, physics, and IC options — see
   the [[Input Parameters Reference|Input-Parameters]] for every field. For the GPU
   build, check the scope rules in "Building (GPU ...)" above (multi-GPU needs periodic
   x and z).
2. Create the required output directories (adjust paths to match your `fileout` and
   `rsb_fileout` settings):
   ```bash
   mkdir -p fields restart stats
   ```
3. Launch with MPI (CPU build, any rank count):
   ```bash
   mpirun -np <N> ./build/dopamine
   ```
   or, for the GPU build (one rank per GPU; use the NVHPC `mpirun`):
   ```bash
   mpirun -np 2 ./build_gpu/dopamine     # 2 GPUs
   ```
   The solver reads `input_parameters` from the working directory in both cases.

Every run prints a per-stage profiler summary (`src/profiler.f90`) at shutdown — wall
time and percent-of-tracked-time for SGS, wall model, RHS, boundary conditions, Poisson
FFT, Poisson tridiagonal solve, projection, and CFL check — useful for spotting where
time is going on either build.

## Testing

`ctest` (in the build directory) runs the unit drivers plus deterministic regression
cases in `tests/regression/` (TGV, LES/DNS channel, 4-wall duct, inflow/outflow, passive
scalar, Boussinesq temperature, UAV actuator disk, IBM sphere, Reynolds-stress budget) and point-particle cases
(tracers, inertial, Boussinesq-coupled, IBM collisions, inflow/outflow with reinjection; compared by particle ID).
Each case is run on 1 rank and on `TEST_PARITY_NPROCS` ranks (default 2), plus explicit
4-rank `2x2` and `4x1` pencil layouts (particles also on 3 ranks), and must agree within tolerance on the per-step
monitor diagnostics **and** field-by-field on the final snapshot (ghost layers excluded;
RSB output files are compared too). A hot-start restart must reproduce an uninterrupted
run. In a GPU build, pass a CPU build to also compare CPU vs GPU at the same rank
count and layout (agreement is to roundoff, ~1e-13):

```bash
cmake -S . -B build_gpu ... -DCPU_REFERENCE_EXE=$PWD/build/dopamine \
      -DCPU_REFERENCE_MPIRUN=$(which mpirun)   # optional; TEST_LD_LIBRARY_PATH adds a libcuda dir
ctest --test-dir build_gpu
```

**Layout independence**: results are independent of the rank count and pencil layout (np=1, 2, 3, 2x2, 4x1, ...) and agree
between the CPU and GPU builds, for LES, IBM, Boussinesq, scalars, UAV, inflow/outflow and point particles; restart reproduces an
uninterrupted run. Details, the debugging tools (`DOPAMINE_TRACE_DIR`, `tests/regression/trace_case.sh`) and what was fixed are in
[[Decomposition Consistency|Decomposition-Consistency]]. Known gaps: recycled precursor inflow (`inflow_type=2`) is not under test,
and Brownian motion / the SGS Langevin particle closure use per-rank random streams (statistically, not bitwise, layout independent).
Run 4-rank GPU tests on 2 GPUs one at a time (`ctest` without `-j`).

## Output files

| Location | Content | Format |
|----------|---------|--------|
| `fields/` | Velocity (U, V, W), pressure (P), and — when `sgs_model /= 0` — SGS turbulent viscosity (ν_t) snapshots, written every `nsave` steps (or every `tsave` time units if `nsave < 0`) | Big-endian float64 stream |
| `restart/` | Hot-restart fields | Big-endian float64 stream |
| `stats/` | Monitor statistics (text) and, if enabled, RSB budget files (little-endian float64) | See [[Input Parameters Reference § STATISTICS|Input-Parameters#statistics-optional--omit-to-disable]] |

Field snapshots are Fortran stream unformatted, big-endian float64 (no record markers).
Each field block is preceded by a 3-integer size header. The layout is:

| Block | Dimensions | Present |
|-------|-----------|---------|
| Grid arrays (x, y, z, xm, ym, zm) | 1-D, each preceded by 1 Int32 | always |
| U | `nx × nyg × nzg` | always |
| V | `nxg × ny × nzg` | always |
| W | `nxg × nyg × nz` | always |
| P | `nxg × nyg × nzg` | always |
| C | `nxg × nyg × nzg` | only when `sediment_flag >= 1` |
| ν_t | `nxg × nyg × nzg` | only when `sgs_model /= 0` |
| T | `nxg × nyg × nzg` | only when `boussinesq_flag >= 1` (written last) |

```python
import numpy as np
data = np.fromfile("fields/channel_test.1", dtype=">f8")
```

The reader library `postProcessing/snapshot_io.py` handles this layout (grid parsing,
ghost-cell stripping) automatically — see
[[Pre- and Post-Processing Tools|Tools#snapshot_iopy]].

### XDMF / ParaView post-processing

`postProcessing/generateXMF.py` generates XDMF metadata files that let ParaView open the
binary snapshots directly. It produces a single `Velocity` vector attribute using the
XDMF `JOIN` function rather than three separate scalars:

```bash
python postProcessing/generateXMF.py --case channel_test --nx 513 --ny 128 --nz 257
```

See [[Pre- and Post-Processing Tools|Tools]] for the full script reference.
