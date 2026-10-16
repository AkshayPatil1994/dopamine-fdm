#!/bin/bash
# usage: [FORCE=1] [MAKE_OPTS='--toe-tk ... --lz 0.6 ...'] ./setup.sh CASE [NP]    CASE = coarse | coarse_fit079 | medium | fine | test (resolution presets of make_case.py; extra make_case.py options via MAKE_OPTS); builds CASE/SDF_in for the solver's grid with GenSDF.
# 1) make_case.py writes input_parameters, bed.stl, gensdf_parameters.in
# 2) a solver start without SDF_in writes fields/geometry.out and fields/grid.out (and aborts at the missing SDF: expected)
# 3) GenSDF turns bed.stl into the cell-centre SDF (positive in the fluid) -> SDF_in
set -e
case=${1:?coarse|coarse_fit079|medium|fine|test}; np=${2:-8}
here=$(cd "$(dirname "$0")" && pwd)
repo=${DOPAMINE_REPO:-/home/Username/fdm-dopamine}
exe=${DOPAMINE_EXE:-$repo/build/dopamine}
gensdf=${GENSDF_EXE:-$repo/preProcessing/GenSDF/gensdf}
# never wipe a case that already holds results: setup.sh deletes CASE/fields and CASE/data (FORCE=1 overrides)
if [ -z "$FORCE" ] && { ls "$here/$case/fields"/tk94.* >/dev/null 2>&1 || [ -s "$here/$case/vof_gauges.dat" ]; }; then
  echo "$here/$case already contains snapshots / gauge records: refusing to overwrite (use another CASE copy or FORCE=1)"; exit 1
fi
cd "$here"; python3 make_case.py $case $MAKE_OPTS
cd $case; rm -rf fields data SDF_in; mkdir -p fields data restart probes
mpirun --oversubscribe -np $np $exe > prerun.log 2>&1 || true
test -s fields/geometry.out || { echo "pre-run did not write fields/geometry.out; see $case/prerun.log"; exit 1; }
cp fields/geometry.out fields/grid.out data/
cp bed.stl data/; cp gensdf_parameters.in parameters.in
sed -i "s#^'bed.stl'#'data/bed.stl'#" parameters.in
mpirun --oversubscribe -np $np $gensdf > gensdf.log 2>&1
mv data/sdfp.bin SDF_in
ls -l SDF_in; echo "$case ready: cd $here/$case && mpirun -np <N> $exe"
