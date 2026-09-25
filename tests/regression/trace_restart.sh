#!/usr/bin/env bash
# Restart bisection: run a case for 2N steps straight and as N steps + hot-start restart of N more (as restart_check.py),
# tracing the first step after the restart point in both, then report the first trace point/field/cell that differs.
# The case input must have nsteps = nsave = N. P and nu_t are not compared: both are recomputed before they are used (they are not restored from the file).
#
#   tests/regression/trace_restart.sh <case-dir> <exe> [layout=2] [tol=1e-11]
set -euo pipefail
case_dir=$(readlink -f "${1:?case dir}"); exe=$(readlink -f "${2:?dopamine executable}")
layout=${3:-2}; tol=${4:-1e-11}
here=$(cd "$(dirname "$0")" && pwd)
out=${TRACE_OUT:-$(mktemp -d)}; echo "trace output: $out"
if [[ $layout == *x* ]]; then row=${layout%x*}; col=${layout#*x}; np=$((row*col)); else np=$layout; row=; fi
n=$(grep -oE 'nsteps *= *[0-9]+' "$case_dir/input_parameters" | grep -oE '[0-9]+$')
fileout=$(grep -oE "fileout *= *'[^']+'" "$case_dir/input_parameters" | sed -E "s/.*'([^']+)'/\1/")
d=$out/run; mkdir -p "$d"/{fields,restart,stats} "$out/trace_a" "$out/trace_b"
cp -r "$case_dir"/. "$d"/
[[ -n $row ]] && sed -i -E "s/p_row *= *[0-9]+/p_row = $row/; s/p_col *= *[0-9]+/p_col = $col/" "$d/input_parameters"
cp "$d/input_parameters" "$out/input.orig"
setn() { sed -E "s/nsteps *= *[0-9]+/nsteps = $1/; s/nsave *= *[0-9]+/nsave = $1/" "$out/input.orig" > "$d/input_parameters"; }
go() { (cd "$d" && mpirun --oversubscribe -np "$np" "$exe" > "$1" 2>&1) || { echo "run failed"; tail -20 "$d/$1"; exit 1; }; }

# A: straight 2N steps, trace step N+1
setn $((2*n)); DOPAMINE_TRACE_DIR="$out/trace_a" DOPAMINE_TRACE_START=$((n+1)) DOPAMINE_TRACE_STEPS=1 go run_a.log
rm -f "$d"/fields/*.[0-9]*
# first half, write the restart file
setn "$n"; go run_half.log
cp "$d/fields/$fileout.$n" "$d/restart/$fileout.$n"
# B: restart, trace its first step
sed -E "s/nsteps *= *[0-9]+/nsteps = $n/; s/nsave *= *[0-9]+/nsave = $n/; s/restart *= *0/restart = 1/; s/nstep_init *= *[0-9]+/nstep_init = $n/; s#filein *= *'[^']*'#filein = 'restart/$fileout.$n'#" "$out/input.orig" > "$d/input_parameters"
DOPAMINE_TRACE_DIR="$out/trace_b" DOPAMINE_TRACE_START=1 DOPAMINE_TRACE_STEPS=1 go run_b.log
python3 "$here/trace_compare.py" "$out/trace_a" "$out/trace_b" --tol "$tol" --skip P --skip nut ${TRACE_ALL:+--all}
