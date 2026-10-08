#!/usr/bin/env python3
"""Mean current + wave at the inlet of the two-fluid flume: the run completes, the liquid-volume ledger closes, and the largest velocity is that of
the log-law current at the still-water level, U ln(d/z0)/(ln(d/z0) - 1 + z0/d), plus the orbital velocity of the small wave, found near the free surface.
Then, for current_type 1 (uniform), 2 (log law, z0) and 3 (power law, exponent n), with a negligible wave (H = 1e-4) and strong relaxation, the z-mean
U(y) of the water cells is compared with the analytic profile (depth mean 0.2) at the inlet face (2e-3) and in the absorption zone, whose target is the
current alone (1e-2: the zone relaxes U to its target within a fraction of the transit time)."""
import argparse, math, os, re, shutil, subprocess, sys, tempfile
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('--exe', required=True)
ap.add_argument('--mpirun', default='mpirun')
a = ap.parse_args()
src = os.path.dirname(os.path.abspath(__file__))
tmp = tempfile.mkdtemp(prefix='vofcur_')
depth, z0, uc, npow = 0.5, 1e-3, 0.2, 7.0


def analytic(typ, y):
    yc = np.minimum(y, depth)
    if typ == 1:
        return uc + 0*y
    if typ == 2:
        return uc*np.log(np.maximum(yc, z0)/z0)/(math.log(depth/z0) - 1 + z0/depth)
    return uc*(npow + 1)/npow*(yc/depth)**(1/npow)


def profiles(wd, typ):
    """z-mean U of the water cells at the inlet face and in the absorption zone, deviation from the analytic profile (snapshot fields/flume.*)"""
    with open(os.path.join(wd, 'fields', 'flume.%d' % nstep_snap(wd)), 'rb') as f:
        mesh = []
        for _ in range(6):
            n = int(np.fromfile(f, '>i4', 1)[0]); mesh.append(np.fromfile(f, '>f8', n))
        n = np.fromfile(f, '>i4', 3)
        U = np.fromfile(f, '>f8', int(np.prod(n))).reshape(tuple(n), order='F')
    x, ym = mesh[0], mesh[4]
    water = ym < depth
    ref = analytic(typ, ym)[water]
    um = U.mean(axis=2)[:, 1:len(ym) + 1][:, water]
    zone = x[:um.shape[0]] > 2.8
    return np.abs(um[0] - ref).max(), np.abs(um[zone] - ref).max()


def nstep_snap(wd):
    return max(int(f.split('.')[1]) for f in os.listdir(os.path.join(wd, 'fields')) if re.match(r'flume\.\d+$', f))


def run_profile(typ):
    wd = os.path.join(tmp, 'type%d' % typ)
    os.makedirs(os.path.join(wd, 'restart'))
    text = open(os.path.join(src, 'input_parameters')).read()
    for pat, rep in [(r'nsave = 100000', 'nsave = -1, tsave = 1.5'), (r'sim_end_time = 2.0', 'sim_end_time = 1.5'),
                     (r'wave_height = 0.01', 'wave_height = 1e-4'), (r'wave_ramp_time = 1.0', 'wave_ramp_time = 0.5'),
                     (r'wave_relax_rate = 20.0', 'wave_relax_rate = 100.0'),
                     (r'current_type = 2', 'current_type = %d, current_n = %g' % (typ, npow))]:
        text = re.sub(pat, rep, text, count=1)
    open(os.path.join(wd, 'input_parameters'), 'w').write(text)
    r = subprocess.run([a.mpirun, '--oversubscribe', '-np', '1', a.exe], cwd=wd, capture_output=True, text=True, timeout=900)
    if r.returncode != 0:
        print(r.stdout[-2000:], r.stderr[-2000:]); sys.exit(1)
    return profiles(wd, typ)


try:
    os.makedirs(os.path.join(tmp, 'restart'))
    shutil.copy(os.path.join(src, 'input_parameters'), tmp)
    r = subprocess.run([a.mpirun, '--oversubscribe', '-np', '1', a.exe], cwd=tmp, capture_output=True, text=True, timeout=900)
    if r.returncode != 0:
        print(r.stdout[-2000:], r.stderr[-2000:]); sys.exit(1)
    d = np.loadtxt(os.path.join(tmp, 'vof_diag.dat'))
    resid = (d[-1, 3] - d[0, 3] - d[-1, 21] - d[-1, 22]) / d[0, 3]
    u_surf = uc*math.log(depth/z0)/(math.log(depth/z0) - 1 + z0/depth)
    umax = d[:, 15].max()
    ok = d[-1, 1] > 1.99 and abs(resid) < 1e-5 and u_surf < umax < u_surf + 0.1 and d[-1, 24] > 0.35
    print('t_end %.2f  ledger residual %.2e  max|U| %.4f  (surface current %.4f)  %s' % (d[-1, 1], resid, umax, u_surf, 'ok' if ok else 'FAIL'))
    for typ, name in ((1, 'uniform'), (2, 'log law'), (3, 'power law')):
        d_in, d_zone = run_profile(typ)
        good = d_in < 2e-3 and d_zone < 1e-2
        ok = ok and good
        print('current_type %d (%s): max |U - analytic| inlet %.2e  absorption zone %.2e  %s' % (typ, name, d_in, d_zone, 'ok' if good else 'FAIL'))
    sys.exit(0 if ok else 1)
finally:
    shutil.rmtree(tmp, ignore_errors=True)
