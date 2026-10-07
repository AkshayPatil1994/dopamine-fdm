#!/usr/bin/env python3
"""Run vof_tgv_mirror on one rank and require the last mid-plane snapshot to be exactly mirror-symmetric about y = Ly/2
(v antisymmetric, u symmetric): any asymmetric wall padding of the momentum stencils shows up at the 1e-4 level."""
import argparse, glob, os, subprocess, sys, tempfile, shutil
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('--exe', required=True)
ap.add_argument('--mpirun', default='mpirun')
ap.add_argument('--tol', type=float, default=1e-12)
a = ap.parse_args()
tmp = tempfile.mkdtemp(prefix='vofmirror_')
try:
    wd = os.path.join(tmp, 'a')
    shutil.copytree(os.path.dirname(os.path.abspath(__file__)), wd)
    os.makedirs(os.path.join(wd, 'restart'), exist_ok=True)
    r = subprocess.run([a.mpirun, '-np', '1', a.exe], cwd=wd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit('FAIL: run exited with %d\n%s' % (r.returncode, r.stdout[-1500:]))
    d = np.loadtxt(sorted(glob.glob(os.path.join(wd, 'vof_snap_*.dat')))[-1])
    nx, ny = len(np.unique(d[:, 0])), len(np.unique(d[:, 1]))
    u = d[:, 3].reshape(ny, nx);  v = d[:, 4].reshape(ny, nx)
    ev, eu = np.abs(v + v[::-1]).max(), np.abs(u - u[::-1]).max()
    print('mirror error: v %.3e  u %.3e  (max|v| %.3e)' % (ev, eu, np.abs(v).max()))
    if max(ev, eu) > a.tol:
        sys.exit('FAIL: solution is not mirror-symmetric')
    print('PASS')
finally:
    shutil.rmtree(tmp, ignore_errors=True)
