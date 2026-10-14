#!/usr/bin/env python3
"""Tracers in a sloshing two-fluid tank: the count is conserved and every particle ends in (or next to) liquid."""
import argparse, os, shutil, subprocess, sys, tempfile
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from compare_monitor import read_snapshot
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
    cv = read_snapshot(w, interior_only=True)[-1]
    ncx, ncy, ncz = 33, 32, 5            # interior cells of nx=34, ny=33, nz=6 (Lx = 1.03125, Ly = 1)
    dx, dy = 1.03125 / 33, 1.0 / 32
    bad = 0
    for pid, v in parts.items():
        i = min(max(int(v[0] / dx), 0), ncx - 1)
        j = min(max(int(v[1] / dy), 0), ncy - 1)
        c = max(cv[i + ncx * (j + ncy * k)] for k in range(ncz))
        bad += c < 0.3
    print('particles %d, in gas: %d' % (len(parts), bad))
    sys.exit(0 if len(parts) == 400 and bad == 0 else 'FAIL')
finally:
    shutil.rmtree(tmp, ignore_errors=True)
