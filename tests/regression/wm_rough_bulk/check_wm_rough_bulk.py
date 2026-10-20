#!/usr/bin/env python3
"""Column momentum balance of the rough EQWM on the Michel-Stadt precursor y grid (Ly = 240, z0 = 1.57, matching height j_m = 52, free-slip
lid, small horizontal extent). From a uniform stream U0 with no turbulence the wall is the only momentum sink, so the bulk speed obeys
d<U>/dt = dPdx - u_tau^2/Ly with u_tau = kappa u_m/ln(y_m/z0) at the matching cell. The change over 10 s must match the time-mean of that
(start and end u_m), the old single-cell stress gives ~+dPdx*T instead, the first cell must not be pinned, and 1 and 4 ranks must agree."""
import argparse, math, os, re, subprocess, sys, tempfile
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('--exe', required=True)
ap.add_argument('--mpirun', default='mpirun')
ap.add_argument('--tol', type=float, default=0.03, help='relative tolerance on the bulk speed change')
a = ap.parse_args()
src = os.path.dirname(os.path.abspath(__file__))
tmp = tempfile.mkdtemp(prefix='wmbulk_')
KAPPA, Z0, LY, DT, NSTEPS, U0, DPDX, JM = 0.41, 1.57, 240.0, 0.05, 200, 7.7, 1.450417e-03, 52


def run(tag, np_, edits):
    wd = os.path.join(tmp, tag)
    os.makedirs(os.path.join(wd, 'restart'))
    text = open(os.path.join(src, 'input_parameters')).read()
    for pat, rep in edits:
        text = re.sub(pat, rep, text)
    open(os.path.join(wd, 'input_parameters'), 'w').write(text)
    r = subprocess.run([a.mpirun, '--oversubscribe', '-np', str(np_), a.exe], cwd=wd, capture_output=True, text=True, timeout=900)
    if r.returncode != 0:
        print(r.stdout[-2000:], r.stderr[-2000:]); sys.exit(1)
    names = [f for f in os.listdir(os.path.join(wd, 'fields')) if re.match(r'x\.\d+$', f)]
    with open(os.path.join(wd, 'fields', max(names, key=lambda f: int(f.split('.')[1]))), 'rb') as f:
        hdr = []
        for _ in range(6):
            n = int(np.fromfile(f, '>i4', 1)[0]); hdr.append(np.fromfile(f, '>f8', n))
        n = np.fromfile(f, '>i4', 3)
        return hdr[1], np.fromfile(f, '>f8', int(np.prod(n))).reshape(tuple(n), order='F')


try:
    yf, U1 = run('np1', 1, [])
    _, U4 = run('np4', 4, [(r'p_row = 1, p_col = 1', 'p_row = 2, p_col = 2')])
    prof = U1[1:-1, :, 1:-1].mean(axis=(0, 2))           # index j-1 for Fortran j (ghost rows at both ends)
    dy = np.diff(yf)                                      # cell j spans yf[j-2] .. yf[j-1]
    bulk = (prof[1:-1]*dy).sum()/LY
    ym, ygm = yf[JM - 1], 0.5*(yf[JM - 2] + yf[JM - 1])
    ut2 = lambda u: (KAPPA*u/math.log(ygm/Z0))**2
    T = DT*NSTEPS
    pred = T*(DPDX - 0.5*(ut2(U0) + ut2(prof[JM - 1]))/LY)
    got = bulk - U0
    err = abs(got - pred)/abs(pred)
    print('bulk change over %.0f s: %.6f  predicted %.6f  rel.err %.3f  (no drag would give %+.6f)' % (T, got, pred, err, T*DPDX))
    first = prof[1]/prof[JM - 1]
    print('first cell / matching cell = %.3f' % first)
    d14 = np.abs(U1 - U4).max()
    print('np1 vs np4 max|dU| %.2e' % d14)
    if not (err < a.tol and first > 0.5 and d14 < 1e-9):
        print('FAILED'); sys.exit(1)
    print('PASSED')
finally:
    import shutil
    shutil.rmtree(tmp, ignore_errors=True)
