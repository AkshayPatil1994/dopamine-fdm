#!/usr/bin/env python3
"""IBM ghost-cell bookkeeping must not depend on the decomposition.

Runs a case's setup (1 step) at two layouts and compares the global IBM ghost-cell counts that setup_ibm prints:
the kept ghost cells (U,V,W) and the dropped ones (image clips solid / image outside the local array). A rank-count
dependence here means some wall cells get no (or a doubled) boundary condition at rank seams, which shows up as a
solution difference from step 1 (see trace_case.sh).

  ibm_ghost_counts.py --case-dir DIR --exe EXE --mpirun MPIRUN --layout-a 1 --layout-b 2x2
"""
import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile


def run_layout(case_dir, exe, mpirun, layout):
    if 'x' in layout:
        row, col = layout.split('x')
        np_ = int(row) * int(col)
    else:
        row = col = None
        np_ = int(layout)
    d = tempfile.mkdtemp(prefix='ibmcount_')
    try:
        for sub in ('fields', 'restart', 'stats'):
            os.makedirs(os.path.join(d, sub))
        for f in os.listdir(case_dir):
            fp = os.path.join(case_dir, f)
            if os.path.isfile(fp):
                shutil.copy(fp, d)
        inp = os.path.join(d, 'input_parameters')
        txt = open(inp).read()
        txt = re.sub(r'nsteps\s*=\s*-?\d+', 'nsteps = 1', txt)
        txt = re.sub(r'nsave\s*=\s*-?\d+', 'nsave = 100000000', txt)
        if row is not None:
            txt = re.sub(r'p_row\s*=\s*\d+', f'p_row = {row}', txt)
            txt = re.sub(r'p_col\s*=\s*\d+', f'p_col = {col}', txt)
        open(inp, 'w').write(txt)
        r = subprocess.run([mpirun, '--oversubscribe', '-np', str(np_), exe], cwd=d, stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, text=True, timeout=600)
        if r.returncode != 0:
            print(r.stdout)
            sys.exit(f'ERROR: layout {layout} exited {r.returncode}')
        out = r.stdout
    finally:
        shutil.rmtree(d, ignore_errors=True)
    kept = re.search(r'ghost cells \(GLOBAL\)\s+\S+\s+U:\s+(\d+)\s+V:\s+(\d+)\s+W:\s+(\d+)', out)
    solid = re.search(r'image clips solid:\s+(\d+) (\d+) (\d+) (\d+)', out)
    outside = re.search(r'image outside local array:\s+(\d+) (\d+) (\d+) (\d+)', out)
    if not (kept and solid and outside):
        sys.exit(f'ERROR: layout {layout}: IBM ghost-count lines not found in the output')
    return {'kept UVW': tuple(map(int, kept.groups())),
            'dropped (clips solid) UVWC': tuple(map(int, solid.groups())),
            'dropped (outside local array) UVWC': tuple(map(int, outside.groups()))}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case-dir', required=True)
    p.add_argument('--exe', required=True)
    p.add_argument('--mpirun', default='mpirun')
    p.add_argument('--layout-a', default='1')
    p.add_argument('--layout-b', default='2')
    a = p.parse_args()
    ca = run_layout(a.case_dir, a.exe, a.mpirun, a.layout_a)
    cb = run_layout(a.case_dir, a.exe, a.mpirun, a.layout_b)
    bad = False
    for k in ca:
        same = ca[k] == cb[k]
        bad |= not same
        print(f'  {k:38s} {a.layout_a}: {ca[k]}  {a.layout_b}: {cb[k]}  {"ok" if same else "DIFFER"}')
    print('RESULT: FAIL' if bad else 'RESULT: PASS')
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()
