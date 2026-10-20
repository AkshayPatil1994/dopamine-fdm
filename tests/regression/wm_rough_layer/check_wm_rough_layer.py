#!/usr/bin/env python3
"""Rough-wall EQWM with the matching height above the first cell (flat_wall_model_flag = 2, free-slip lid): from a uniform stream with
no pressure gradient and no turbulence, each cell up to the matching height j_m obeys dU/dt = -u_tau^2/y_m with u_tau = kappa U/ln(yg_m/z0),
i.e. U(t) = U0/(1 + a U0 t), a = (kappa/ln(yg_m/z0))^2/y_m. The first cell must follow this too (not be pinned to zero), the cells above j_m
must stay at U0, and 1 and 4 ranks must agree. The mirrored case (rough wall on top, free-slip bottom) must give the mirrored profile."""
import argparse, math, os, re, subprocess, sys, tempfile
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('--exe', required=True)
ap.add_argument('--mpirun', default='mpirun')
ap.add_argument('--tol', type=float, default=5e-3)
ap.add_argument('--vof', action='store_true', help='run the two-fluid solver (all liquid) and also check independence of the liquid density')
a = ap.parse_args()
src = os.path.dirname(os.path.abspath(__file__))
tmp = tempfile.mkdtemp(prefix='wmrough_')
VOF = '''&VOF
  vof_active = 1, vof_flow = 1, vof_ic_type = 1, vof_level = 64.0,
  vof_rho_l = 1000.0, vof_rho_g = 1.2, vof_nu_l = 1.0e-5, vof_nu_g = 1.5e-5, vof_grav = 0.0, vof_pcg_tol = 1e-10, vof_pcg_iters = 400,
/
'''
KAPPA, Z0, NY, LY, DT, NSTEPS, U0 = 0.41, 0.1, 65, 64.0, 0.2, 500, 1.0


def run(tag, np_, edits):
    wd = os.path.join(tmp, tag)
    os.makedirs(os.path.join(wd, 'restart'))
    text = open(os.path.join(src, 'input_parameters')).read()
    if a.vof:
        text += VOF
    for pat, rep in edits:
        text = re.sub(pat, rep, text)
    open(os.path.join(wd, 'input_parameters'), 'w').write(text)
    r = subprocess.run([a.mpirun, '--oversubscribe', '-np', str(np_), a.exe], cwd=wd, capture_output=True, text=True, timeout=900)
    if r.returncode != 0:
        print(r.stdout[-2000:], r.stderr[-2000:]); sys.exit(1)
    names = [f for f in os.listdir(os.path.join(wd, 'fields')) if re.match(r'x\.\d+$', f)]
    with open(os.path.join(wd, 'fields', max(names, key=lambda f: int(f.split('.')[1]))), 'rb') as f:
        for _ in range(6):
            n = int(np.fromfile(f, '>i4', 1)[0]); np.fromfile(f, '>f8', n)
        n = np.fromfile(f, '>i4', 3)
        return np.fromfile(f, '>f8', int(np.prod(n))).reshape(tuple(n), order='F')


try:
    dy = LY/(NY - 1)
    jm = next(j for j in range(2, NY//2 + 2) if (j - 1.5)*dy/Z0 >= 20)   # first cell with y/z0 >= 20 (Fortran index)
    ym, ygm = (jm - 1)*dy, (jm - 1.5)*dy                                   # top face of cell jm, centre of cell jm
    aa = (KAPPA/math.log(ygm/Z0))**2/ym
    T = DT*NSTEPS
    exact = U0/(1 + aa*U0*T)
    U1 = run('np1', 1, [])
    U4 = run('np4', 4, [(r'p_row = 1, p_col = 1', 'p_row = 2, p_col = 2')])
    Ut = run('top', 1, [(r'bc_face_ylo = 1, bc_face_yhi = 2', 'bc_face_ylo = 2, bc_face_yhi = 1'), (r'z0_ylo = 0.1', 'z0_yhi = 0.1')])
    prof = U1[1:-1, :, 1:-1].mean(axis=(0, 2))     # U(j-1), Fortran j = 1..nyg
    bad = False
    if a.vof:
        Ur = run('rho', 1, [(r'vof_rho_l = 1000.0', 'vof_rho_l = 7.0')])
        drho = np.abs(Ur - U1).max()
        print('liquid density 1000 vs 7 max|dU| %.2e' % drho)
        bad |= not (drho < 1e-8)
    ptop = Ut[1:-1, :, 1:-1].mean(axis=(0, 2))
    mirror = np.abs(ptop[::-1][1:NY] - prof[1:NY]).max()   # cell j <-> nyg+1-j
    print('top-wall run vs mirrored bottom-wall run: max|dU| %.2e' % mirror)
    bad |= not (mirror < 1e-6)
    for j in range(2, jm + 1):
        err = abs(prof[j - 1] - exact)/exact
        print('cell j=%d  U %.6f  exact %.6f  rel.err %.2e' % (j, prof[j - 1], exact, err))
        bad |= not (err < a.tol)
    above = np.abs(prof[jm:NY] - U0).max()
    print('j_m = %d, exact U(T) = %.6f; max |U - U0| above j_m: %.2e' % (jm, exact, above))
    bad |= not (above < 1e-3)
    d14 = np.abs(U1 - U4).max()
    print('np1 vs np4 max|dU| %.2e' % d14)
    bad |= not (d14 < 1e-9)
    if bad:
        print('FAILED'); sys.exit(1)
    print('PASSED')
finally:
    import shutil
    shutil.rmtree(tmp, ignore_errors=True)
