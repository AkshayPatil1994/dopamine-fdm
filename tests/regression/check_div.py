#!/usr/bin/env python3
"""Run one case and require the monitor's max |div| (all cells that the solver reports) to stay below --max-div at every step.
Used for ibm_pressure = 2, whose fluid cells must stay divergence-free to round-off even next to the body."""
import argparse, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_monitor import parse_monitor, run_case, clear_snapshots

p = argparse.ArgumentParser()
p.add_argument('--case-dir', required=True)
p.add_argument('--mpirun', default='mpirun')
p.add_argument('--exe', required=True)
p.add_argument('--np', type=int, default=1)
p.add_argument('--max-div', type=float, default=1e-9)
a = p.parse_args()

clear_snapshots(a.case_dir)
rows = parse_monitor(run_case(a.mpirun, a.exe, a.np, a.case_dir, {}))
if not rows:
    sys.exit('ERROR: no monitor lines parsed')
worst = max(r[3] for r in rows.values())   # columns: t, meanU, maxU, max_divergence, ...
print(f'max |div| over {len(rows)} steps: {worst:.3e} (limit {a.max_div:.1e})')
sys.exit(0 if worst <= a.max_div else 1)
