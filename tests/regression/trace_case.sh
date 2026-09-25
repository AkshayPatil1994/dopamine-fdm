#!/usr/bin/env bash
# Decomposition bisection: run one regression case at two rank layouts with stage tracing on, then report the first
# trace point / field / cell at which the two runs differ (see src/debug_trace.f90, trace_compare.py).
#
#   tests/regression/trace_case.sh <case-dir> <exe> [layout_a=1] [layout_b=2] [steps=2] [tol=1e-10]
#
# A layout is a rank count ("2") or a pencil grid ("2x2" = p_row x p_col). Runs in scratch copies of the case;
# dumps are kept under $TRACE_OUT (default: a mktemp dir, path printed).
set -euo pipefail
case_dir=$(readlink -f "${1:?case dir}"); exe=$(readlink -f "${2:?dopamine executable}")
la=${3:-1}; lb=${4:-2}; steps=${5:-2}; tol=${6:-1e-10}
here=$(cd "$(dirname "$0")" && pwd)
out=${TRACE_OUT:-$(mktemp -d)}
echo "trace output: $out"

run_layout() {  # layout tag
  local layout=$1 tag=$2 np row col d
  if [[ $layout == *x* ]]; then row=${layout%x*}; col=${layout#*x}; np=$((row*col)); else np=$layout; row=; fi
  d=$out/run_$tag; mkdir -p "$d"/{fields,restart,stats} "$out/trace_$tag"
  cp -r "$case_dir"/. "$d"/
  sed -i -E "s/nsteps *= *-?[0-9]+/nsteps = $steps/; s/nsave *= *-?[0-9]+/nsave = 100000000/" "$d/input_parameters"
  if [[ -n $row ]]; then
    sed -i -E "s/p_row *= *[0-9]+/p_row = $row/; s/p_col *= *[0-9]+/p_col = $col/" "$d/input_parameters"
  fi
  (cd "$d" && DOPAMINE_TRACE_DIR="$out/trace_$tag" DOPAMINE_TRACE_STEPS=$steps \
     mpirun --oversubscribe -np "$np" "$exe" > run.log 2>&1) || { echo "layout $layout FAILED"; tail -20 "$d/run.log"; exit 1; }
}
run_layout "$la" a
run_layout "$lb" b
python3 "$here/trace_compare.py" "$out/trace_a" "$out/trace_b" --tol "$tol" ${TRACE_ALL:+--all}
