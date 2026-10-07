#!/usr/bin/env python3
"""Submerged sphere in still water: the IBM force output (ibm_forces.csv, Method 2: pressure + viscous) must be the buoyancy

rho_l g V pointing up (+y), and the horizontal loads must vanish. The sphere of the synthetic SDF has radius 0.8
(centre between two cell rows, so -min(sdf) = 0.7375 underestimates it); the staircase body of the estimator is within a few percent."""
import argparse, os, re, shutil, subprocess, sys, tempfile
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('--exe', required=True)
ap.add_argument('--mpirun', default='mpirun')
ap.add_argument('--tol', type=float, default=0.06)
ap.add_argument('--np', type=int, default=1)
a = ap.parse_args()
src = os.path.dirname(os.path.abspath(__file__))
text = re.sub(r'!.*', '', open(os.path.join(src, 'input_parameters')).read())
rho_l = float(re.search(r'vof_rho_l\s*=\s*([\d.eE+-]+)', text).group(1))
g = float(re.search(r'vof_grav\s*=\s*([\d.eE+-]+)', text).group(1))
R = 0.8
expected = rho_l*g*4.0/3.0*np.pi*R**3
tmp = tempfile.mkdtemp(prefix='vofbuoy_')
try:
    wd = os.path.join(tmp, 'a')
    shutil.copytree(src, wd)
    os.makedirs(os.path.join(wd, 'restart'), exist_ok=True)
    r = subprocess.run([a.mpirun, '--oversubscribe', '-np', str(a.np), a.exe], cwd=wd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit('FAIL: run exited with %d\n%s' % (r.returncode, r.stdout[-1500:]))
    d = np.loadtxt(os.path.join(wd, 'ibm_forces.csv'), delimiter=',', skiprows=1)
    d = d.reshape(-1, d.shape[-1])[-1]
    fp, fv = d[5:8], d[8:11]     # columns: step,t, Fx,Fy,Fz (impulse), pressure xyz, viscous xyz
    f = fp + fv
    print('R %.4f  expected buoyancy %.1f   Fy(pres) %.1f  Fy(visc) %.2f  Fx %.2f  Fz %.2f   ratio %.4f' %
          (R, expected, fp[1], fv[1], f[0], f[2], f[1]/expected))
    ok = abs(f[1]/expected - 1) < a.tol and abs(f[0]) < a.tol*expected and abs(f[2]) < a.tol*expected
    if not ok:
        sys.exit('FAIL: submerged-sphere load is not the buoyancy')
    print('PASS')
finally:
    shutil.rmtree(tmp, ignore_errors=True)
