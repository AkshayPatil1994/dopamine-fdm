#!/usr/bin/env python3
"""Staircase IBM (ibm_method = 1): a wall one cell thick (SDF -h/2 next to +h/2, a tie in the face average) closes its faces like a thicker one.
Fences across the stream of 1, 2 and 3 cells stop it (mean U ~ 0) with round-off divergence and the pressure impulse of the fluid volume that stops;
1 and 4 ranks agree."""
import argparse, os, re, shutil, subprocess, sys, tempfile
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('--exe', required=True)
ap.add_argument('--mpirun', default='mpirun')
a = ap.parse_args()
src = os.path.dirname(os.path.abspath(__file__))
tmp = tempfile.mkdtemp(prefix='stairfence_')
NX, NY, NZ, I0, H = 18, 33, 18, 10, 1/32   # SDF file shape (NX+1, NY+1, NZ-1); first solid cell index (1-based)


def sdf(path, n):
    phi = np.full((NX + 1, NY + 1, NZ - 1), H/2)
    phi[I0-1:I0-1+n, :, :] = -H/2
    phi.astype('>f8').flatten(order='F').tofile(path)


def run(tag, np_, n):
    wd = os.path.join(tmp, tag)
    os.makedirs(os.path.join(wd, 'restart'))
    text = open(os.path.join(src, 'input_parameters')).read().replace('ibm_method = 1,', 'ibm_method = 1, nsampling = 1,')
    if np_ > 1:
        text = text.replace('p_row = 0, p_col = 0', 'p_row = 2, p_col = 2')
    open(os.path.join(wd, 'input_parameters'), 'w').write(text)
    sdf(os.path.join(wd, 'SDF_in'), n)
    r = subprocess.run([a.mpirun, '--oversubscribe', '-np', str(np_), a.exe], cwd=wd, capture_output=True, text=True, timeout=900)
    if r.returncode != 0:
        print(r.stdout[-2000:], r.stderr[-2000:]); sys.exit(1)
    rows = [l.split() for l in r.stdout.splitlines() if re.match(r'^\s*\d+\s+[-\d.E+]+\s+[-\d.E+]+\s+[-\d.E+]+\s', l)]
    force = np.loadtxt(os.path.join(wd, 'ibm_forces.csv'), delimiter=',', skiprows=1)
    return abs(float(rows[-1][2])), float(rows[-1][4]), force[:, 5].sum()*2e-3   # mean U, divergence, pressure impulse


try:
    res = {n: run('n%d' % n, 1, n) for n in (1, 2, 3)}
    u4, _, imp4 = run('n1_np4', 4, 1)
    ok = True
    for n, (u, div, imp) in res.items():
        print('fence %d cells: mean U %.2e  div %.1e  pressure impulse %.4f' % (n, u, div, imp))
        ok = ok and u < 1e-3 and div < 1e-9
    # the impulse is the momentum of the fluid between the fence copies: it scales with the fluid volume 16 - n cells
    ok = ok and all(abs(res[n][2]/res[1][2] - (16 - n)/15) < 0.02 for n in (2, 3))
    ok = ok and abs(imp4 - res[1][2]) < 1e-6*abs(res[1][2]) and abs(u4 - res[1][0]) < 1e-9
    print('np4 fence-1 impulse %.4f' % imp4)
    print('ok' if ok else 'FAIL')
    sys.exit(0 if ok else 1)
finally:
    shutil.rmtree(tmp, ignore_errors=True)
