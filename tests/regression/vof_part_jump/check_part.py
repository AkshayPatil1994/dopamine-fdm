#!/usr/bin/env python3
"""Heavy particles seeded in the gas of a two-fluid tank fall through the air and sink into the liquid."""
import argparse, os, shutil, subprocess, sys, tempfile
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from compare_particles import read_particles

ap = argparse.ArgumentParser()
ap.add_argument('--exe', required=True)
ap.add_argument('--mpirun', default='mpirun')
ap.add_argument('--np', type=int, default=1)
a = ap.parse_args()
tmp = tempfile.mkdtemp(prefix='vofpart_')
try:
    w = os.path.join(tmp, 'run')
    shutil.copytree(os.path.dirname(os.path.abspath(__file__)), w)
    for d in ('restart', 'stats'):
        os.makedirs(os.path.join(w, d), exist_ok=True)
    r = subprocess.run([a.mpirun, '--oversubscribe', '-np', str(a.np), a.exe], cwd=w, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-1500:], r.stderr[-1500:])
        sys.exit('FAIL: run exited with %d' % r.returncode)
    parts = read_particles(os.path.join(w, 'particles_restart'))
    for f in os.listdir(os.path.join(w, 'fields')):
        if '_particles' in f:
            os.remove(os.path.join(w, 'fields', f))
    ys = sorted(v[1] for v in parts.values())
    print('particle height range [%.4f, %.4f], median %.4f (surface at 0.5)' % (ys[0], ys[-1], ys[len(ys)//2]))
    sys.exit(0 if len(parts) == 200 and ys[-1] < 0.5 else 'FAIL: particles must have entered the liquid (y < 0.5)')
finally:
    shutil.rmtree(tmp, ignore_errors=True)
