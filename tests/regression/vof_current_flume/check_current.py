#!/usr/bin/env python3
"""Mean current + wave at the inlet of the two-fluid flume: the run completes, the liquid-volume ledger closes, and the largest velocity is that of
the log-law current at the still-water level, U ln(d/z0)/(ln(d/z0) - 1 + z0/d), plus the orbital velocity of the small wave, found near the free surface."""
import argparse, math, os, shutil, subprocess, sys, tempfile
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('--exe', required=True)
ap.add_argument('--mpirun', default='mpirun')
a = ap.parse_args()
src = os.path.dirname(os.path.abspath(__file__))
tmp = tempfile.mkdtemp(prefix='vofcur_')
try:
    os.makedirs(os.path.join(tmp, 'restart'))
    shutil.copy(os.path.join(src, 'input_parameters'), tmp)
    r = subprocess.run([a.mpirun, '--oversubscribe', '-np', '1', a.exe], cwd=tmp, capture_output=True, text=True, timeout=900)
    if r.returncode != 0:
        print(r.stdout[-2000:], r.stderr[-2000:]); sys.exit(1)
    d = np.loadtxt(os.path.join(tmp, 'vof_diag.dat'))
    resid = (d[-1, 3] - d[0, 3] - d[-1, 21] - d[-1, 22]) / d[0, 3]
    depth, z0, uc = 0.5, 1e-3, 0.2
    u_surf = uc*math.log(depth/z0)/(math.log(depth/z0) - 1 + z0/depth)
    umax = d[:, 15].max()
    ok = d[-1, 1] > 1.99 and abs(resid) < 1e-5 and u_surf < umax < u_surf + 0.1 and d[-1, 24] > 0.35
    print('t_end %.2f  ledger residual %.2e  max|U| %.4f  (surface current %.4f)  %s' % (d[-1, 1], resid, umax, u_surf, 'ok' if ok else 'FAIL'))
    sys.exit(0 if ok else 1)
finally:
    shutil.rmtree(tmp, ignore_errors=True)
