#!/usr/bin/env python3
"""The snapshot pressure of a 1-rank run and of a 2x2 run must agree everywhere, seam planes included (a stale ghost plane of the
pressure array used to overwrite the neighbour's last real plane in the file, leaving whole zero planes at the rank seams)."""
import argparse, glob, os, re, shutil, struct, subprocess, sys, tempfile
import numpy as np


def read_pressure(path):
    with open(path, 'rb') as fh:
        ri = lambda: struct.unpack('>i', fh.read(4))[0]
        for _ in range(3):
            fh.seek(8*ri(), 1)
        for _ in range(3):
            n = ri(); fh.seek(8*n, 1)
        for nm in ('U', 'V', 'W', 'P'):
            d = struct.unpack('>3i', fh.read(12))
            a = np.fromfile(fh, '>f8', d[0]*d[1]*d[2]).reshape(d, order='F')
    return a


def run(src, exe, mpirun, np_, pgrid, wd):
    shutil.copytree(src, wd)
    os.makedirs(os.path.join(wd, 'restart'), exist_ok=True)
    f = os.path.join(wd, 'input_parameters')
    t = open(f).read()
    if pgrid:
        t = re.sub(r'p_row\s*=\s*\d+', 'p_row = %s' % pgrid[0], t); t = re.sub(r'p_col\s*=\s*\d+', 'p_col = %s' % pgrid[1], t)
        open(f, 'w').write(t)
    cmd = [mpirun, '--oversubscribe', '-np', str(np_), exe]
    r = subprocess.run(cmd, cwd=wd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit('FAIL: run exited with %d (np=%d)\n%s' % (r.returncode, np_, r.stdout[-1500:]))
    return read_pressure(sorted(glob.glob(os.path.join(wd, 'fields', 'vof_p_seam.*')))[-1])


ap = argparse.ArgumentParser()
ap.add_argument('--exe', required=True)
ap.add_argument('--mpirun', default='mpirun')
ap.add_argument('--tol', type=float, default=1e-6)
a = ap.parse_args()
src = os.path.dirname(os.path.abspath(__file__))
tmp = tempfile.mkdtemp(prefix='vofpseam_')
try:
    p1 = run(src, a.exe, a.mpirun, 1, None, os.path.join(tmp, 'a'))
    p4 = run(src, a.exe, a.mpirun, 4, ('2', '2'), os.path.join(tmp, 'b'))
    inner = (slice(1, -1), slice(1, -1), slice(1, -1))
    zero_planes = [(ax, i) for ax in (0, 2) for i in range(1, p4.shape[ax]-1)
                   if not np.take(p4[inner], i-1, axis=ax).any() and np.take(p1[inner], i-1, axis=ax).any()]
    err = np.abs(p1[inner] - p4[inner]).max()/max(np.abs(p1[inner]).max(), 1e-300)
    print('max|P| %.3e, relative difference np1 vs 2x2 %.3e, zero planes only in 2x2: %s' % (np.abs(p1).max(), err, zero_planes))
    if zero_planes or err > a.tol:
        sys.exit('FAIL: seam error in the snapshot pressure')
    print('PASS')
finally:
    shutil.rmtree(tmp, ignore_errors=True)
