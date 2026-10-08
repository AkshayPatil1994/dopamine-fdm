#!/usr/bin/env python3
"""Rough-wall MOST heat-flux wall model (T_bc = 2, single phase): starting from uniform T = T_ref, the heat gained by the column (the
flux the discretisation applies at the wall face) must equal the iterated Businger-Dyer flux u_tau*theta_tau evaluated at the matching
height of the final state, for heated/cooled bottom and top walls, and be the same on 1 and 4 ranks."""
import argparse, math, os, re, subprocess, sys, tempfile
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('--exe', required=True)
ap.add_argument('--mpirun', default='mpirun')
ap.add_argument('--tol', type=float, default=0.01, help='relative tolerance on the applied wall flux')
a = ap.parse_args()
src = os.path.dirname(os.path.abspath(__file__))
tmp = tempfile.mkdtemp(prefix='wmthermal_')
KAPPA, G, TREF, Z0, Z0H, NY, LY, DT, NSTEPS = 0.41, 9.81, 300.0, 1e-3, 1e-4, 33, 1.0, 2e-3, 40


def psi(zeta):
    z = max(min(zeta, 2.0), -2.0)
    if z < 0:
        x = (1 - 16*z)**0.25
        return 2*math.log((1 + x)/2) + math.log((1 + x*x)/2) - 2*math.atan(x) + math.pi/2, 2*math.log((1 + x*x)/2)
    return -5*z, -5*z


def most(u, dtheta, y):
    L = 1e10
    for _ in range(200):
        pm, ph = psi(y/math.copysign(max(abs(L), 1e-3), L))
        ut = KAPPA*u/max(math.log(max(y, 2*Z0)/Z0) - pm, 1e-3)
        tt = KAPPA*dtheta/max(math.log(max(y, 2*Z0H)/Z0H) - ph, 1e-3)
        Ln = ut**2*TREF/(KAPPA*G*tt) if abs(tt) > 1e-12 else math.copysign(1e10, L)
        Ln = math.copysign(max(abs(Ln), 1e-3), Ln)
        done = abs(Ln - L) < 1e-9*abs(Ln)
        L = Ln
        if done:
            break
    return ut*tt


def read_fields(wd):
    names = [f for f in os.listdir(os.path.join(wd, 'fields')) if re.match(r'x\.\d+$', f)]
    with open(os.path.join(wd, 'fields', max(names, key=lambda f: int(f.split('.')[1]))), 'rb') as f:
        for _ in range(6):
            n = int(np.fromfile(f, '>i4', 1)[0]); np.fromfile(f, '>f8', n)
        blk = []
        for _ in range(6):   # U, V, W, P, nu_t, T
            n = np.fromfile(f, '>i4', 3)
            blk.append(np.fromfile(f, '>f8', int(np.prod(n))).reshape(tuple(n), order='F'))
    return blk[0], blk[5]


def run(tag, np_, edits):
    wd = os.path.join(tmp, tag)
    os.makedirs(os.path.join(wd, 'restart'))
    text = open(os.path.join(src, 'input_parameters')).read()
    for pat, rep in edits:
        text = re.sub(pat, rep, text, count=1)
    open(os.path.join(wd, 'input_parameters'), 'w').write(text)
    r = subprocess.run([a.mpirun, '--oversubscribe', '-np', str(np_), a.exe], cwd=wd, capture_output=True, text=True, timeout=900)
    if r.returncode != 0:
        print(r.stdout[-2000:], r.stderr[-2000:]); sys.exit(1)
    return read_fields(wd)


def analyse(U, T, top):
    nyg = NY + 1
    dy = LY/(NY - 1)
    yg = lambda j: (j - 1.5)*dy                       # cell centres, Fortran index j = 2 .. nyg-1
    jm = 2
    for j in range(2, nyg//2 + 1):                    # matching height: first cell with y/z0 >= 20 (lower half-channel)
        jm = j
        if yg(j)/Z0 >= 20:
            break
    ym = yg(jm)                                       # distance to the wall; the top wall mirrors the lower half-channel
    jm = nyg + 1 - jm if top else jm
    u = 0.5*(U[0:-1, jm - 1, :] + U[1:, jm - 1, :])[1:-1, 1:-1].mean()
    tm = T[1:-1, jm - 1, 1:-1].mean()
    heat = ((T[1:-1, 1:nyg - 1, 1:-1].mean(axis=(0, 2)) - TREF)*dy).sum()/(DT*NSTEPS)
    return heat, u, tm, ym


try:
    cases = [('heated', [], False), ('cooled', [(r'T_wall_bot = 301.0', 'T_wall_bot = 299.0')], False),
             ('strong', [(r'T_wall_bot = 301.0', 'T_wall_bot = 330.0')], False),
             ('top', [(r'T_bc_bot = 2, T_bc_top = 0', 'T_bc_bot = 0, T_bc_top = 2'), (r'T_wall_top = 300.0', 'T_wall_top = 299.0')], True)]
    bad = False
    for tag, edits, top in cases:
        U, T = run(tag, 1, edits)
        heat, u, tm, ym = analyse(U, T, top)
        twall = float(re.search(r'T_wall_%s = ([\d.]+)' % ('top' if top else 'bot'),
                                open(os.path.join(tmp, tag, 'input_parameters')).read()).group(1))
        qm = -most(u, tm - twall, ym)                 # heat gained by the fluid per unit area: -u_tau*theta_tau
        err = abs(heat - qm)/abs(qm)
        print('%-8s applied %.5e  MOST %.5e  rel.err %.3f' % (tag, heat, qm, err))
        bad |= not (err < a.tol)
    U1, T1 = run('np1', 1, [])
    U4, T4 = run('np4', 4, [(r'p_row = 1, p_col = 1', 'p_row = 2, p_col = 2')])
    d14 = np.abs(T1 - T4).max()
    print('np1 vs np4 max|dT| %.2e' % d14)
    bad |= not (d14 < 1e-9)
    if bad:
        print('FAILED'); sys.exit(1)
    print('PASSED')
finally:
    import shutil
    shutil.rmtree(tmp, ignore_errors=True)
