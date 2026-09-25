#!/usr/bin/env python3
"""Particle parity: run a case with &PARTICLES in two configurations (rank count / pencil layout / executable, e.g. CPU vs
GPU) and compare the final particle state by particle ID.

The run writes `particles_restart` at the last step (nsave == nsteps in the regression cases): int32 count, int32 ids[count],
float64 (x, y, z, u, v, w, age) per particle. Particles are matched by ID (their order depends on the decomposition).
The two runs must hold the same set of IDs (particles removed or reinjected identically) and agree in every column to --tol.

Exit code 0 = agree, 1 = mismatch or a run failure.
"""
import argparse
import array
import os
import re
import shutil
import struct
import subprocess
import sys
import tempfile

COLS = ['x', 'y', 'z', 'u', 'v', 'w', 'age']


def read_particles(path):
    with open(path, 'rb') as f:
        raw = f.read()
    for order in ('<', '>'):
        (n,) = struct.unpack(order + 'i', raw[:4])
        if 0 <= n < 10**8 and 4 + 4 * n + 56 * n == len(raw):
            break
    else:
        sys.exit(f'ERROR: {path}: cannot parse particle file ({len(raw)} bytes)')
    ids = struct.unpack(order + f'{n}i', raw[4:4 + 4 * n])
    dat = array.array('d')
    dat.frombytes(raw[4 + 4 * n:])
    if (order == '>') == (sys.byteorder == 'little'):
        dat.byteswap()
    return {ids[i]: dat[7 * i:7 * i + 7].tolist() for i in range(n)}


def run(mpirun, exe, np_, case_dir, grid, env_extra, nsteps=None):
    d = tempfile.mkdtemp(prefix='part_')
    for sub in ('fields', 'restart', 'stats'):
        os.makedirs(os.path.join(d, sub))
    for f in os.listdir(case_dir):
        fp = os.path.join(case_dir, f)
        if os.path.isfile(fp):
            shutil.copy(fp, d)
    inp = os.path.join(d, 'input_parameters')
    txt = open(inp).read()
    if nsteps:
        txt = re.sub(r'nsteps\s*=\s*-?\d+', f'nsteps = {nsteps}', txt)
        txt = re.sub(r'nsave\s*=\s*-?\d+', f'nsave = {nsteps}', txt)
    if grid:
        txt = re.sub(r'p_row\s*=\s*\d+', f'p_row = {grid[0]}', txt)
        txt = re.sub(r'p_col\s*=\s*\d+', f'p_col = {grid[1]}', txt)
    open(inp, 'w').write(txt)
    env = os.environ.copy()
    env.update(env_extra)
    r = subprocess.run([mpirun, '--oversubscribe', '-np', str(np_), exe], cwd=d, env=env, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, text=True, timeout=900)
    if r.returncode != 0:
        print(r.stdout[-3000:])
        sys.exit(f'ERROR: run ({exe}, np={np_}) exited {r.returncode}')
    pf = os.path.join(d, 'particles_restart')
    if not os.path.exists(pf):
        sys.exit('ERROR: no particles_restart written (is particles_active=1 and nsave reached?)')
    parts = read_particles(pf)
    shutil.rmtree(d, ignore_errors=True)
    return parts


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case-dir', required=True)
    p.add_argument('--exe-a', required=True)
    p.add_argument('--exe-b', required=True)
    p.add_argument('--np-a', type=int, default=1)
    p.add_argument('--np-b', type=int, default=2)
    p.add_argument('--mpirun', default='mpirun')
    p.add_argument('--mpirun-b', default=None)
    p.add_argument('--p-grid', default=None, help='ROW,COL pencil grid for run B (and A with --p-grid-both)')
    p.add_argument('--p-grid-both', action='store_true')
    p.add_argument('--nsteps', type=int, default=None, help='override nsteps (and nsave) of the case')
    p.add_argument('--tol', type=float, default=1e-9, help='max absolute difference per column')
    p.add_argument('--env-b', action='append', default=[], help='KEY=VALUE for run B, repeatable')
    p.add_argument('--label-a', default='A')
    p.add_argument('--label-b', default='B')
    a = p.parse_args()

    grid = tuple(a.p_grid.split(',')) if a.p_grid else None
    envb = dict(kv.partition('=')[::2] for kv in a.env_b)
    pa = run(a.mpirun, a.exe_a, a.np_a, a.case_dir, grid if a.p_grid_both else None, {}, a.nsteps)
    pb = run(a.mpirun_b or a.mpirun, a.exe_b, a.np_b, a.case_dir, grid, envb, a.nsteps)

    ok = True
    print(f'{a.label_a}: {len(pa)} particles, {a.label_b}: {len(pb)} particles')
    if set(pa) != set(pb):
        only_a, only_b = sorted(set(pa) - set(pb)), sorted(set(pb) - set(pa))
        print(f'MISMATCH particle IDs: only in {a.label_a}: {only_a[:8]} ({len(only_a)}), '
              f'only in {a.label_b}: {only_b[:8]} ({len(only_b)})')
        ok = False
    common = sorted(set(pa) & set(pb))
    if not common:
        sys.exit('ERROR: no particles to compare')
    for c, name in enumerate(COLS):
        worst, wid = 0.0, None
        for i in common:
            d = abs(pa[i][c] - pb[i][c])
            if d > worst:
                worst, wid = d, i
        flag = '' if worst <= a.tol else f'  MISMATCH (> {a.tol:g}, particle {wid})'
        ok &= worst <= a.tol
        print(f'  {name:4s} max |diff| = {worst:.3e}{flag}')
    print('RESULT: PASS' if ok else 'RESULT: FAIL')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
