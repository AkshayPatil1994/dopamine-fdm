#!/usr/bin/env python3
"""Single-phase input checks: each invalid setting must abort at input with its message, and the unmodified tgv_small must not."""
import argparse, os, re, shutil, subprocess, sys, tempfile

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'tgv_small')
CASES = [
    ('nu missing', r'nu\s*=\s*0.01,', '', 'nu \\(kinematic viscosity\\) must be set'),
    ('dt missing', r'dt\s*=\s*0.01,', '', 'dt must be > 0'),
    ('sgs_model range', r'sgs_model\s*=\s*0', 'sgs_model = 5', 'sgs_model must be'),
    ('advection_scheme range', r'sgs_model\s*=\s*0', 'sgs_model = 0, advection_scheme = 4', 'advection_scheme must be'),
    ('adv_order range', r'sgs_model\s*=\s*0', 'sgs_model = 0, advection_scheme = 2, adv_order = 4', 'adv_order must be'),
    ('adv_band range', r'sgs_model\s*=\s*0', 'sgs_model = 0, advection_scheme = 2, adv_band = 9', 'adv_band must be'),
    ('adv_beta0 negative', r'sgs_model\s*=\s*0', 'sgs_model = 0, advection_scheme = 2, adv_beta0 = -0.1', 'adv_beta0 and adv_beta_ibm'),
    ('flat_wall_model_flag range', r'flat_wall_model_flag\s*=\s*0', 'flat_wall_model_flag = 7', 'flat_wall_model_flag must be'),
    ('x_bc_type range', r'x_bc_type\s*=\s*0', 'x_bc_type = 3', 'x_bc_type must be'),
    ('ibm_pressure range', r'&PHYSICS', '&IBM\n  ibm_pressure = 5,\n/\n&PHYSICS', 'ibm_pressure must be'),
    ('ibm_pressure without body', r'&PHYSICS', '&IBM\n  ibm_pressure = 1,\n/\n&PHYSICS', 'needs an immersed body'),
    ('removed ibm_projection', r'&PHYSICS', '&IBM\n  ibm_projection = 1,\n/\n&PHYSICS', None),
]

ap = argparse.ArgumentParser()
ap.add_argument('--exe', required=True)
ap.add_argument('--mpirun', default='mpirun')
a = ap.parse_args()
tmp = tempfile.mkdtemp(prefix='inputabort_')
bad = 0
try:
    for name, pat, rep, msg in CASES:
        wd = os.path.join(tmp, re.sub(r'\W', '_', name))
        os.makedirs(os.path.join(wd, 'restart'))
        text = open(os.path.join(BASE, 'input_parameters')).read()
        new = re.sub(pat, lambda m: rep, text, count=1)
        assert new != text, name
        open(os.path.join(wd, 'input_parameters'), 'w').write(new)
        r = subprocess.run([a.mpirun, '-np', '1', a.exe], cwd=wd, capture_output=True, text=True, timeout=120)
        out = r.stdout + r.stderr
        ok = r.returncode != 0 and (msg is None or re.search(msg, out) is not None)
        print('%-30s %s' % (name, 'ok' if ok else 'FAIL'))
        bad += not ok
    wd = os.path.join(tmp, 'ref')
    os.makedirs(os.path.join(wd, 'restart'))
    text = re.sub(r'nsteps\s*=\s*-?\d+', 'nsteps = 1', open(os.path.join(BASE, 'input_parameters')).read())
    open(os.path.join(wd, 'input_parameters'), 'w').write(text)
    r = subprocess.run([a.mpirun, '-np', '1', a.exe], cwd=wd, capture_output=True, text=True, timeout=120)
    print('%-30s %s' % ('unmodified case runs', 'ok' if r.returncode == 0 else 'FAIL'))
    bad += r.returncode != 0
finally:
    shutil.rmtree(tmp, ignore_errors=True)
sys.exit(1 if bad else 0)
