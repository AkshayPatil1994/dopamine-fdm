#!/usr/bin/env python3
"""Start-up SGS advisory (src/monitor.f90 summary): one-step tiny runs must print exactly the advisory lines the flow case warrants
(fine/coarse streamwise resolution in wall units for the wall-modelled channel, stratified AMD, wall-resolved AMD) and nothing else."""
import argparse, os, re, subprocess, sys, tempfile

ap = argparse.ArgumentParser()
ap.add_argument('--exe', required=True)
ap.add_argument('--mpirun', default='mpirun')
a = ap.parse_args()
src = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'chan_les_amd_small', 'input_parameters')
text0 = open(src).read()
BOUSS = """&BOUSSINESQ
  boussinesq_flag = 1, beta_T = 3.0e-3, T_ref = 0.0, Pr = 0.7, Pr_t = 0.85,
  T_bc_bot = 1, T_bc_top = 1, T_wall_bot = 1.0, T_wall_top = 0.0, T_ic_type = 1, T_ic_grad = -0.5,
/
"""
KEYS = {'fine': 'AMD with fine streamwise resolution', 'coarse': 'coarse streamwise resolution', 'strat': 'stratified AMD',
        'notval': 'not validated for wall-resolved, IBM or two-phase'}
# tag: (edits, expected advisory keys)   dx+ = Lx/(nx-1)*u_tau/nu = 0.176*395 = 70 (Lx=3) or 279 (Lx=12) with u_tau = sqrt(dPdx*Ly/2) = 1
cases = {
    'fine_amd':    ([], {'fine'}),
    'fine_vreman': ([(r'sgs_model            = 2', 'sgs_model            = 1')], set()),
    'coarse_amd':  ([(r'Lx         = 3.0', 'Lx         = 12.0')], set()),
    'coarse_vre':  ([(r'Lx         = 3.0', 'Lx         = 12.0'), (r'sgs_model            = 2', 'sgs_model            = 1')], {'coarse'}),
    'resolved':    ([(r'flat_wall_model_flag = 1', 'flat_wall_model_flag = 0')], {'notval'}),
    'stratified':  ([(r'Lx         = 3.0', 'Lx         = 12.0'), ('APPEND', BOUSS)], {'strat'}),
}
tmp = tempfile.mkdtemp(prefix='sgsadv_')
bad = False
try:
    for tag, (edits, expect) in cases.items():
        wd = os.path.join(tmp, tag)
        os.makedirs(os.path.join(wd, 'restart'))
        text = text0
        for pat, rep in edits:
            text = text + '\n' + rep if pat == 'APPEND' else re.sub(pat, rep, text, count=1)
        text = re.sub(r'nsteps       = \d+', 'nsteps       = 1', text)
        open(os.path.join(wd, 'input_parameters'), 'w').write(text)
        r = subprocess.run([a.mpirun, '--oversubscribe', '-np', '1', a.exe], cwd=wd, capture_output=True, text=True, timeout=300)
        if r.returncode != 0:
            print(tag, 'run failed'); print(r.stdout[-1500:], r.stderr[-1500:]); sys.exit(1)
        lines = [l for l in r.stdout.splitlines() if 'ADVISORY:' in l]
        got = {k for k, s in KEYS.items() if any(s in l for l in lines)}
        ok = got == expect and len(lines) == len(expect)
        print('%-12s expected %-10s got %-10s %s' % (tag, sorted(expect), sorted(got), 'ok' if ok else 'MISMATCH'))
        bad |= not ok
finally:
    import shutil
    shutil.rmtree(tmp, ignore_errors=True)
if bad:
    print('FAILED'); sys.exit(1)
print('PASSED')
