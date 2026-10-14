#!/usr/bin/env python3
"""Passive scalar under the two-fluid solver: a uniform scalar (no settling) stays uniform in both phases while the interface sloshes,
and the liquid volume is conserved. Reads the last snapshot (blocks U,V,W,P,C,..., Cv last)."""
import argparse, os, shutil, subprocess, sys, tempfile
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from compare_monitor import read_snapshot


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--exe', required=True)
    ap.add_argument('--mpirun', default='mpirun')
    ap.add_argument('--np', type=int, default=1)
    a = ap.parse_args()
    tmp = tempfile.mkdtemp(prefix='vofsc_')
    try:
        w = os.path.join(tmp, 'run')
        shutil.copytree(os.path.dirname(os.path.abspath(__file__)), w)
        os.makedirs(os.path.join(w, 'restart'), exist_ok=True)
        r = subprocess.run([a.mpirun, '--oversubscribe', '-np', str(a.np), a.exe], cwd=w, capture_output=True, text=True)
        if r.returncode != 0:
            print(r.stdout[-1500:], r.stderr[-1500:])
            sys.exit('FAIL: run exited with %d' % r.returncode)
        blocks = read_snapshot(w, interior_only=True)
        c, cv = blocks[4], blocks[-1]
        dev = max(abs(x - 1.0) for x in c)
        print('uniform: max|C-1| = %.2e   liquid fraction range [%.3f, %.3f]' % (dev, min(cv), max(cv)))
        ok = dev < 1e-10
        # a slab of scalar below the still-water level: bounded, and its integral (uniform grid) conserved while the interface sloshes
        w2 = os.path.join(tmp, 'slab')
        shutil.copytree(w, w2, ignore=shutil.ignore_patterns('fields'))
        os.makedirs(os.path.join(w2, 'restart'), exist_ok=True)
        inp = os.path.join(w2, 'input_parameters')
        t = open(inp).read().replace('C_ic_type = 0', 'C_ic_type = 3, C_ic_height = 0.5').replace('tsave = 0.3', 'tsave = 0.15')
        open(inp, 'w').write(t)
        r = subprocess.run([a.mpirun, '--oversubscribe', '-np', str(a.np), a.exe], cwd=w2, capture_output=True, text=True)
        if r.returncode != 0:
            sys.exit('FAIL: slab run exited with %d' % r.returncode)
        import glob
        files = sorted(glob.glob(os.path.join(w2, 'fields', '*.[0-9]*')), key=lambda f: int(f.rsplit('.', 1)[1]))
        tot = []
        for f in (files[0], files[-1]):
            d = os.path.join(tmp, 'one'); shutil.rmtree(d, ignore_errors=True); os.makedirs(os.path.join(d, 'fields'))
            shutil.copy(f, os.path.join(d, 'fields'))
            cc = read_snapshot(d, interior_only=True)[4]
            tot.append((sum(cc), min(cc), max(cc)))
        drift = abs(tot[1][0] - tot[0][0]) / tot[0][0]
        print('slab: integral drift %.2e, range [%.4f, %.4f] -> [%.4f, %.4f]' % (drift, tot[0][1], tot[0][2], tot[1][1], tot[1][2]))
        ok = ok and drift < 1e-3 and tot[1][1] > -1e-3 and tot[1][2] < 1 + 1e-3
        sys.exit(0 if ok else 'FAIL')
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


main()
