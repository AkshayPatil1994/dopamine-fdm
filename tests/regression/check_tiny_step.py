#!/usr/bin/env python3
"""Run the ibm_stair_tiny_step case (ibm_pressure = 2, time-based saving with tsave = 0.01 + 1e-9: the save time lies 1e-9 s beyond the fifth step)
and require that every step runs and that the landing snapshot, the newest one, has a pressure of the order of the dynamic pressure. A landing
step of 1e-9 s would multiply the incremental pressure by dt_prev/dt (about 1e6 here, 1e7 in an urban case): the velocity stays sane in
that step and the next one blows up. The landing step is never shorter than dt/2 (monitor.f90 landing_dt)."""
import argparse, atexit, glob, os, re, sys, shutil, tempfile
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_monitor import parse_monitor, run_case, clear_snapshots, read_snapshot

p = argparse.ArgumentParser()
p.add_argument('--case-dir', required=True)
p.add_argument('--mpirun', default='mpirun')
p.add_argument('--exe', required=True)
p.add_argument('--np', type=int, default=1)
p.add_argument('--nsteps', type=int, default=8)
p.add_argument('--max-p', type=float, default=1e3)
a = p.parse_args()

scratch = tempfile.mkdtemp(prefix='case_')   # own work directory, so ctest -j cannot overlap runs of one case
shutil.copytree(a.case_dir, os.path.join(scratch, 'case'))
a.case_dir = os.path.join(scratch, 'case')
atexit.register(shutil.rmtree, scratch, ignore_errors=True)

clear_snapshots(a.case_dir)
rows = parse_monitor(run_case(a.mpirun, a.exe, a.np, a.case_dir, {}))
if len(rows) < a.nsteps:
    sys.exit(f'FAIL: only {len(rows)} of {a.nsteps} steps ran (blow-up)')
steps = sorted(int(m.group(1)) for f in glob.glob(os.path.join(a.case_dir, 'fields', '*.[0-9]*'))
               if (m := re.search(r'\.(\d+)$', f)))
if not steps or steps[-1] % 5 == 0:
    sys.exit(f'FAIL: no landing snapshot was written (snapshot steps {steps}); the test no longer exercises a short step')
umax = max(r[2] for r in rows.values())
pmax = max(np.abs(f).max() for f in read_snapshot(a.case_dir, interior_only=True)[3:4])
print(f'steps {len(rows)}, landing snapshot at step {steps[-1]}, max |U| {umax:.3f}, max |P| {pmax:.3e} (limit {a.max_p:.1e})')
sys.exit(0 if (umax < 100 and pmax < a.max_p) else 1)
