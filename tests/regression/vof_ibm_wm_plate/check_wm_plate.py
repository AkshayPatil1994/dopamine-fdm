#!/usr/bin/env python3
"""IBM wall model in the two-fluid solver: a uniform stream over a flat immersed plate (fluid layer 0.75 thick) decays by the wall stress
rho u_tau^2 of the log law (Reichardt) at the first fluid cell, whatever the liquid density, and identically on 1 and 4 ranks."""
import argparse, math, os, re, shutil, subprocess, sys, tempfile

ap = argparse.ArgumentParser()
ap.add_argument('--exe', required=True)
ap.add_argument('--mpirun', default='mpirun')
ap.add_argument('--rough', type=float, default=0.0, help='roughness length z0 of the bed (0: smooth Reichardt law)')
ap.add_argument('--tol', type=float, default=0.05, help='relative tolerance on the velocity deficit')
a = ap.parse_args()
src = os.path.dirname(os.path.abspath(__file__))
tmp = tempfile.mkdtemp(prefix='vofplate_')


def run(tag, np_, edits):
    wd = os.path.join(tmp, tag)
    os.makedirs(os.path.join(wd, 'restart'))
    text = open(os.path.join(src, 'input_parameters')).read()
    for pat, rep in edits:
        text = re.sub(pat, rep, text, count=1)
    open(os.path.join(wd, 'input_parameters'), 'w').write(text)
    shutil.copy(os.path.join(src, 'SDF_in'), wd)
    r = subprocess.run([a.mpirun, '--oversubscribe', '-np', str(np_), a.exe], cwd=wd, capture_output=True, text=True, timeout=900)
    if r.returncode != 0:
        print(r.stdout[-2000:], r.stderr[-2000:]); sys.exit(1)
    csv = os.path.join(wd, 'ibm_forces.csv')
    imp = sum(float(l.split(',')[8]) for l in open(csv).read().splitlines()[1:]) if os.path.exists(csv) else 0.0
    rows = [l.split() for l in r.stdout.splitlines() if re.match(r'^\s*\d+\s+[-\d.E+]+\s+[-\d.E+]+\s+[-\d.E+]+\s', l)]
    return int(rows[-1][0]), float(rows[-1][1]), float(rows[-1][2]), imp


def uplus(yp):
    k, c = 0.41, 5.2 - math.log(0.41)/0.41
    return math.log(1 + k*yp)/k + c*(1 - math.exp(-yp/11) - yp/11*math.exp(-0.33*yp))


def u_tau(u, y, nu):
    if a.rough > 0:
        return 0.41*u/math.log(y/a.rough)
    lo, hi = 1e-6*u, u
    for _ in range(200):
        m = 0.5*(lo + hi)
        lo, hi = (m, hi) if u/m - uplus(m*y/nu) > 0 else (lo, m)
    return 0.5*(lo + hi)


try:
    rough = [(r'ibm_input_mode = 1,', 'ibm_input_mode = 1, nsampling = 1,')]
    step, t, u_np1, imp = run('np1', 1, rough)
    _, _, u_rho, _ = run('rho', 1, rough + [(r'vof_rho_l = 1000.0', 'vof_rho_l = 7.0')])
    _, _, u_np4, _ = run('np4', 4, rough + [(r'p_row = 0, p_col = 0', 'p_row = 2, p_col = 2')])
    ly, y_ref, nu, dt, u = 0.75, 0.5*1.0/32, 1e-6, 2e-3, 1.0
    for _ in range(step):
        u -= dt*u_tau(u, y_ref, nu)**2/ly
    ok = abs((1 - u_np1) - (1 - u)) <= a.tol*(1 - u) and abs(u_rho - u_np1) < 1e-9 and abs(u_np4 - u_np1) < 1e-9
    # the viscous load on the plate times dt is the momentum the liquid lost; plate area (16/17)^2
    loss = 1000.0*(1 - u_np1)*ly*(16/17)**2
    ok = ok and abs(imp*dt - loss) <= 0.03*loss
    print('load impulse %.4e  momentum loss %.4e' % (imp*dt, loss))
    print('step %d t=%.3f  mean U: run %.6f  log-law %.6f  rho_l=7 %.6f  np4 %.6f  %s' % (step, t, u_np1, u, u_rho, u_np4, 'ok' if ok else 'FAIL'))
    sys.exit(0 if ok else 1)
finally:
    shutil.rmtree(tmp, ignore_errors=True)
