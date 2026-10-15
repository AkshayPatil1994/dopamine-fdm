#!/usr/bin/env python3
"""advection_scheme = 2 with adv_beta0 = adv_beta_ibm = 0 adds nothing: its monitor lines are identical to advection_scheme = 0."""
import argparse, os, re, shutil, subprocess, sys, tempfile

ap = argparse.ArgumentParser()
ap.add_argument('--exe', required=True)
ap.add_argument('--mpirun', default='mpirun')
a = ap.parse_args()
base = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'tgv_small', 'input_parameters')).read()
tmp = tempfile.mkdtemp(prefix='advbeta0_')


def mon(tag, extra):
    wd = os.path.join(tmp, tag)
    os.makedirs(os.path.join(wd, 'restart'))
    open(os.path.join(wd, 'input_parameters'), 'w').write(base.replace('flat_wall_model_flag = 0,', 'flat_wall_model_flag = 0, ' + extra))
    r = subprocess.run([a.mpirun, '-np', '1', a.exe], cwd=wd, capture_output=True, text=True, timeout=300)
    if r.returncode:
        print(r.stdout[-1500:], r.stderr[-500:]); sys.exit(1)
    return [' '.join(l.split()[:-1]) for l in r.stdout.splitlines() if re.match(r'^\s*\d+\s+[-\d.E+]+\s', l)]


try:
    m0 = mon('s0', 'advection_scheme = 0,')
    m2 = mon('s2', 'advection_scheme = 2, adv_order = 5, adv_beta0 = 0.0, adv_beta_ibm = 0.0,')
    ok = len(m0) > 0 and m0 == m2
    print('monitor lines identical: %s (%d steps)' % (ok, len(m0)))
    sys.exit(0 if ok else 1)
finally:
    shutil.rmtree(tmp, ignore_errors=True)
