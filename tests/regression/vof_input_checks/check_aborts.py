#!/usr/bin/env python3
"""Each unsupported feature combination of the VOF solver must abort at input with its message (and the unmodified case must not)."""
import argparse, os, re, shutil, subprocess, sys, tempfile

CASES = [
    ('vof_flow without vof_active', r'vof_active\s*=\s*1', 'vof_active = 0', 'requires vof_active=1'),
    ('UAV with VOF', r'&VOF', '&UAV\n  uav_active = 1,\n/\n&VOF', 'cannot be combined with the UAV'),
    ('wave inlet without wave', r'x_bc_type\s*=\s*0', 'x_bc_type = 1', None),
    ('wave_type without inflow', r'&VOF', '&WAVES\n  wave_type = 1, wave_height = 0.01, wave_period = 1.0,\n/\n&VOF', 'requires x_bc_type=1 and inflow_type=3'),
    ('staircase IBM with vof_flow', r'&VOF', '&IBM\n  ibm_method = 1,\n/\n&VOF', 'always uses the staircase wall stress'),
    ('y periodic', r'y_bc_type\s*=\s*1', 'y_bc_type = 0', 'wall-bounded y only'),
    ('IBM surface dump', r'&VOF', '&IBM\n  ibm_input_mode = 1, ibm_sdf_file = \'none\', ibm_surface_nsampling = 10,\n/\n&VOF', 'IBM surface field dump'),
]

ap = argparse.ArgumentParser()
ap.add_argument('--exe', required=True)
ap.add_argument('--mpirun', default='mpirun')
a = ap.parse_args()
src = os.path.dirname(os.path.abspath(__file__))
tmp = tempfile.mkdtemp(prefix='vofabort_')
bad = 0
try:
    for name, pat, rep, msg in CASES:
        if msg is None:
            continue
        wd = os.path.join(tmp, re.sub(r'\W', '_', name))
        shutil.copytree(src, wd)
        os.makedirs(os.path.join(wd, 'restart'), exist_ok=True)
        f = os.path.join(wd, 'input_parameters')
        text = open(f).read()
        new = re.sub(pat, lambda m: rep, text, count=1)
        assert new != text, name
        open(f, 'w').write(new)
        r = subprocess.run([a.mpirun, '-np', '1', a.exe], cwd=wd, capture_output=True, text=True, timeout=120)
        out = r.stdout + r.stderr
        ok = r.returncode != 0 and msg in out
        print('%-30s %s' % (name, 'ok' if ok else 'FAIL'))
        bad += not ok
finally:
    shutil.rmtree(tmp, ignore_errors=True)
sys.exit(1 if bad else 0)
