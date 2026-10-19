#!/usr/bin/env python3
"""Reversed flow through the outlet: with outflow_backflow_beta = 1 max |U| stays bounded; with 0 (Dirichlet p = 0) the same case blows up."""
import argparse, os, shutil, sys, tempfile
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from compare_monitor import parse_monitor, run_case

p = argparse.ArgumentParser()
p.add_argument('--mpirun', default='mpirun')
p.add_argument('--exe', required=True)
a = p.parse_args()
src = os.path.dirname(os.path.abspath(__file__))

def peak(beta):
    with tempfile.TemporaryDirectory() as d:
        for f in ('SDF_in', 'input_parameters'):
            shutil.copy(os.path.join(src, f), d)
        for sub in ('fields', 'restart'):
            os.makedirs(os.path.join(d, sub))
        inp = os.path.join(d, 'input_parameters')
        txt = open(inp).read().replace('outflow_backflow_beta = 1.0', f'outflow_backflow_beta = {beta}')
        open(inp, 'w').write(txt)
        rows = parse_monitor(run_case(a.mpirun, a.exe, 1, d, {}))
    return max(r[2] for r in rows.values())   # columns: t, meanU, maxU, ...

stable, unstable = peak(1.0), peak(0.0)
print(f'max |U|: beta=1 {stable:.3e}, beta=0 {unstable:.3e}')
sys.exit(0 if stable < 20.0 and unstable > 1e3 else 1)
