#!/usr/bin/env python3
"""Wall-modelled spanwise (z) walls (single phase): a uniform stream along x decays by the log-law (Reichardt) wall stress u_tau^2 of the
first-cell velocity at the z walls (--duct: at the y walls too), identically on 1 and 4 ranks (z split; x split for the duct).
--reject-rough: the rough EQWM with z walls must be refused at input."""
import argparse, math, os, re, shutil, subprocess, sys, tempfile

ap = argparse.ArgumentParser()
ap.add_argument('--exe', required=True)
ap.add_argument('--mpirun', default='mpirun')
ap.add_argument('--duct', action='store_true', help='y walls as well (4-wall duct)')
ap.add_argument('--reject-rough', action='store_true')
ap.add_argument('--tol', type=float, default=0.05, help='relative tolerance on the velocity deficit')
a = ap.parse_args()
src = os.path.dirname(os.path.abspath(__file__))
tmp = tempfile.mkdtemp(prefix='wmstreamz_')


def run(tag, np_, edits, fail_ok=False):
    wd = os.path.join(tmp, tag)
    os.makedirs(os.path.join(wd, 'restart'))
    text = open(os.path.join(src, 'input_parameters')).read()
    for pat, rep in edits:
        text = re.sub(pat, rep, text, count=1)
    open(os.path.join(wd, 'input_parameters'), 'w').write(text)
    r = subprocess.run([a.mpirun, '--oversubscribe', '-np', str(np_), a.exe], cwd=wd, capture_output=True, text=True, timeout=900)
    if fail_ok:
        return r
    if r.returncode != 0:
        print(r.stdout[-2000:], r.stderr[-2000:]); sys.exit(1)
    rows = [l.split() for l in r.stdout.splitlines() if re.match(r'^\s*\d+\s+[-\d.E+]+\s+[-\d.E+]+\s+[-\d.E+]+\s', l)]
    return int(rows[-1][0]), float(rows[-1][1]), float(rows[-1][2])


def uplus(yp):
    k, c = 0.41, 5.2 - math.log(0.41)/0.41
    return math.log(1 + k*yp)/k + c*(1 - math.exp(-yp/11) - yp/11*math.exp(-0.33*yp))


def u_tau(u, y, nu):
    lo, hi = 1e-6*u, u
    for _ in range(200):
        m = 0.5*(lo + hi)
        lo, hi = (m, hi) if u/m - uplus(m*y/nu) > 0 else (lo, m)
    return 0.5*(lo + hi)


try:
    if a.reject_rough:
        r = run('rough', 1, [(r'flat_wall_model_flag = 1', 'flat_wall_model_flag = 2, z0_ylo = 1e-4, z0_yhi = 1e-4')], fail_ok=True)
        ok = r.returncode != 0 and 'not yet supported for z walls' in r.stdout + r.stderr
        print('rough EQWM with z walls %s' % ('rejected' if ok else 'NOT rejected: FAIL'))
        sys.exit(0 if ok else 1)
    if a.duct:
        e1 = [(r'y_bc_type = 0', 'y_bc_type = 1'), (r'ny = 18', 'ny = 34')]
        e4 = e1 + [(r'p_row = 1, p_col = 1', 'p_row = 4, p_col = 1')]
        nwall = 4
    else:
        e1 = []
        e4 = [(r'p_row = 1, p_col = 1', 'p_row = 1, p_col = 4')]
        nwall = 2
    step, t, u_np1 = run('np1', 1, e1)
    _, _, u_np4 = run('np4', 4, e4)
    l, y_ref, nu, dt, u = 1.0, 0.5*1.0/33, 1e-6, 2e-3, 1.0
    for _ in range(step):
        u -= dt*nwall*u_tau(u, y_ref, nu)**2/l
    ok = abs(u_np1 - u) <= a.tol*(1 - u) and abs(u_np4 - u_np1) < 1e-9
    print('step %d t=%.3f  mean U: run %.6f  log-law %.6f  np4 %.6f  %s' % (step, t, u_np1, u, u_np4, 'ok' if ok else 'FAIL'))
    sys.exit(0 if ok else 1)
finally:
    shutil.rmtree(tmp, ignore_errors=True)
