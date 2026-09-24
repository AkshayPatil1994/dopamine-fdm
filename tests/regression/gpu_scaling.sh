#!/usr/bin/env bash
# Strong-scaling timing of a dopamine GPU build: runs one case at several rank counts (one rank per GPU when
# enough GPUs are visible) and prints tracked profiler time per stage and the per-step wall time.
#
#   tests/regression/gpu_scaling.sh <case-dir-with-input_parameters> <exe> [nsteps=100] [np list = "1 2"]
#
# Needs the launcher of the build (NVHPC mpirun) first on PATH, e.g.
#   export PATH=$NVHPC/comm_libs/openmpi4/bin:$PATH
# The case's nsteps/nsave are overridden (fixed dt from the input is NOT changed; use a case with cfl_adaptive
# as needed). Output goes to a scratch directory that is removed afterwards.
set -euo pipefail
case_dir=${1:?case dir}; exe=${2:?dopamine executable}; nsteps=${3:-100}; nps=${4:-"1 2"}
exe=$(readlink -f "$exe"); case_dir=$(readlink -f "$case_dir")
scratch=$(mktemp -d)
trap 'rm -rf "$scratch"' EXIT
printf '%-4s %12s %10s %8s %8s %8s %8s %8s\n' np "tracked[s]" "s/step" RHS SGS PoissonFFT tridiag other
for np in $nps; do
  d=$scratch/np$np; mkdir -p "$d"/{fields,restart,stats}
  cp -r "$case_dir"/. "$d"/
  sed -i -E "s/nsteps *= *-?[0-9]+/nsteps = $nsteps/; s/nsave *= *-?[0-9]+/nsave = 100000000/" "$d/input_parameters"
  (cd "$d" && mpirun -np "$np" "$exe" > run.log 2>&1) || { echo "np=$np FAILED (see below)"; tail -20 "$d/run.log"; continue; }
  LC_ALL=C awk -v np="$np" -v n="$nsteps" '
    /Total tracked/ {tot=$(NF-1)}
    /RHS \(equations\)/ {rhs=$4} /SGS model +[0-9]/ {sgs=$3} /Poisson FFT/ {fft=$3} /Poisson tridiag/ {tri=$3}
    END {printf "%-4s %12.2f %10.3f %8.2f %8.2f %8.2f %8.2f %8.2f\n", np, tot, tot/n, rhs, sgs, fft, tri, tot-rhs-sgs-fft-tri}' "$d/run.log"
done
