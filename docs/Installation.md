# Installation & Running

← [Home](Home.md)

## Dependencies

| Library | Minimum version | Purpose |
|---------|----------------|---------|
| gfortran or ifort | GFortran ≥ 9 / Intel ≥ 2019 | Fortran compiler |
| MPI | any standard MPI-3 | Domain decomposition |
| FFTW3 | 3.3 (serial double-precision) | Local single-rank transforms/DCT, and the transform engine inside 2decomp&fft |
| [2decomp&fft](https://github.com/2decomp-fft/2decomp-fft) | `v2.1.0` | 2-D pencil domain decomposition and inter-rank transposes for the MPI-parallel pressure Poisson solve — see [Numerics § MPI parallelism](Numerics.md#10-mpi-parallelism) |
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

## Building (single-GPU, OpenACC + cuFFT + cuSPARSE)

> Not tested extensively — watch for bugs.

An optional `ENABLE_GPU` CMake target offloads the RHS/SGS/boundary-condition/projection
kernels and the pressure Poisson solve (cuFFT transform + cuSPARSE batched tridiagonal
solve) to a single NVIDIA GPU via OpenACC. The CPU build above is untouched and remains
the default; this is a separate, opt-in build.

**Requirements**:
- [NVIDIA HPC SDK](https://developer.nvidia.com/hpc-sdk) 23.3 or later (provides
  `nvfortran`, its bundled OpenMPI, and the `cufft`/`cusparse` device libraries) — e.g.
  installed under `/opt/nvidia/hpc_sdk`.
- A CUDA-capable NVIDIA GPU. `CMakeLists.txt` targets compute capability `cc86` (Ampere,
  e.g. RTX 30-series); edit the `-gpu=cc86` flag in `CMakeLists.txt` to match a different
  GPU before building for another architecture.

**Scope — read before using**: the GPU Poisson solve currently only supports
**`nprocs=1`** (single GPU, no MPI domain decomposition), with either **`x_bc_type=0`**
(periodic, spectral FFT) or **`x_bc_type=1`** (inflow/outflow, DCT-IV) streamwise BC. A
runtime guard `Stop`s immediately on startup if `nprocs/=1` or `x_bc_type` is anything
other than 0/1 — every other solver feature (IBM, wall models, SGS, sediment transport,
RSB statistics, etc.) works normally in the GPU build.

Put the NVHPC SDK's `nvfortran` and bundled MPI on your `PATH`/`LD_LIBRARY_PATH` first,
e.g.:

```bash
export PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/23.3/comm_libs/openmpi4/openmpi-4.0.5/bin:/opt/nvidia/hpc_sdk/Linux_x86_64/23.3/compilers/bin:$PATH
export LD_LIBRARY_PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/23.3/comm_libs/openmpi4/openmpi-4.0.5/lib:/opt/nvidia/hpc_sdk/Linux_x86_64/23.3/compilers/lib:$LD_LIBRARY_PATH
```

Then configure and build with `nvfortran` and `-DENABLE_GPU=ON`:

```bash
cmake -S . -B build_gpu -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_Fortran_COMPILER=nvfortran -DENABLE_GPU=ON
cmake --build build_gpu -j$(nproc)
```

The executable `dopamine` is placed in `build_gpu/`. `-Minfo=accel` output during the
build lists every OpenACC region the compiler offloaded — useful to confirm a kernel you
touched is still being generated for the GPU.

## Running

1. Edit `input_parameters` to set domain size, grid type, physics, and IC options — see
   the [Input Parameters Reference](Input-Parameters.md) for every field. For the GPU
   build, also make sure `nprocs=1` (i.e. launch with `-np 1`) and `x_bc_type` is 0 or 1
   in `&BOUNDARY_CONDITIONS` (see "Building (single-GPU...)" above).
2. Create the required output directories (adjust paths to match your `fileout` and
   `rsb_fileout` settings):
   ```bash
   mkdir -p fields restart stats
   ```
3. Launch with MPI (CPU build, any rank count):
   ```bash
   mpirun -np <N> ./build/dopamine
   ```
   or, for the GPU build (single rank only):
   ```bash
   mpirun -np 1 ./build_gpu/dopamine
   ```
   The solver reads `input_parameters` from the working directory in both cases.

Every run prints a per-stage profiler summary (`src/profiler.f90`) at shutdown — wall
time and percent-of-tracked-time for SGS, wall model, RHS, boundary conditions, Poisson
FFT, Poisson tridiagonal solve, projection, and CFL check — useful for spotting where
time is going on either build.

## Output files

| Location | Content | Format |
|----------|---------|--------|
| `fields/` | Velocity (U, V, W), pressure (P), and — when `sgs_model /= 0` — SGS turbulent viscosity (ν_t) snapshots, written every `nsave` steps (or every `tsave` time units if `nsave < 0`) | Big-endian float64 stream |
| `restart/` | Hot-restart fields | Big-endian float64 stream |
| `stats/` | Monitor statistics (text) and, if enabled, RSB budget files (little-endian float64) | See [Input Parameters Reference § STATISTICS](Input-Parameters.md#statistics-optional--omit-to-disable) |

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

```python
import numpy as np
data = np.fromfile("fields/channel_test.1", dtype=">f8")
```

The reader library `postProcessing/snapshot_io.py` handles this layout (grid parsing,
ghost-cell stripping) automatically — see
[Pre- and Post-Processing Tools](Tools.md#snapshot_iopy).

### XDMF / ParaView post-processing

`postProcessing/generateXMF.py` generates XDMF metadata files that let ParaView open the
binary snapshots directly. It produces a single `Velocity` vector attribute using the
XDMF `JOIN` function rather than three separate scalars:

```bash
python postProcessing/generateXMF.py --case channel_test --nx 513 --ny 128 --nz 257
```

See [Pre- and Post-Processing Tools](Tools.md) for the full script reference.
