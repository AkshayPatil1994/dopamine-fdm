#!/usr/bin/env python3
"""Uniform T = T_ref under the two-fluid Boussinesq coupling must leave the run bit-identical to the same case without Boussinesq."""
import argparse, os, re, shutil, subprocess, sys, tempfile


def run(src, exe, mpirun, np, flag, work):
    shutil.copytree(src, work)
    os.makedirs(os.path.join(work, 'restart'), exist_ok=True)
    inp = os.path.join(work, 'input_parameters')
    text = re.sub(r'boussinesq_flag\s*=\s*\d', 'boussinesq_flag = %d' % flag, open(inp).read())
    open(inp, 'w').write(text)
    r = subprocess.run([mpirun, '--oversubscribe', '-np', str(np), exe], cwd=work, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit('FAIL: run exited with %d\n%s' % (r.returncode, r.stdout[-1500:]))
    return open(os.path.join(work, 'vof_diag.dat')).read()


ap = argparse.ArgumentParser()
ap.add_argument('--exe', required=True)
ap.add_argument('--mpirun', default='mpirun')
ap.add_argument('--np', type=int, default=1)
a = ap.parse_args()
tmp = tempfile.mkdtemp(prefix='vofbs_')
try:
    src = os.path.dirname(os.path.abspath(__file__))
    on = run(src, a.exe, a.mpirun, a.np, 1, os.path.join(tmp, 'on'))
    off = run(src, a.exe, a.mpirun, a.np, 0, os.path.join(tmp, 'off'))
    print('uniform T: vof_diag.dat %s' % ('identical' if on == off else 'DIFFERS'))
    sys.exit(0 if on == off else 'FAIL')
finally:
    shutil.rmtree(tmp, ignore_errors=True)
