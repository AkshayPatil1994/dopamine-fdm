#!/usr/bin/env python3
"""Staircase IBM (ibm_method = 1): per-object roughness on the closed faces. A rough block (id 1) stands on a smooth plate (id 2); the same
problem mirrored in x with the stream reversed must give the same |U|: the roughness id of a closed face is the solid cell's, whichever of
the two cells the face joins is looked up."""
import argparse, os, re, shutil, subprocess, sys, tempfile
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('--exe', required=True)
ap.add_argument('--mpirun', default='mpirun')
a = ap.parse_args()
src = os.path.dirname(os.path.abspath(__file__))
tmp = tempfile.mkdtemp(prefix='stairids_')
NX, NY, NZ, H = 18, 33, 18, 1/32   # file shape (NX+1, NY+1, NZ-1)


def fields(wd, mirror):
    phi = np.full((NX + 1, NY + 1, NZ - 1), H/2)
    oid = np.zeros_like(phi)
    phi[:, :9, :] = -H/2;  oid[:, :9, :] = 2
    phi[9:13, 9:12, :] = -H/2;  oid[9:13, 9:12, :] = 1
    if mirror:
        phi, oid = phi[::-1].copy(), oid[::-1].copy()
        phi[-1], oid[-1] = phi[0], oid[0]   # periodic duplicate column
    for name, f in (('SDF_in', phi), ('OBJ_in', oid)):
        f.astype('>f8').flatten(order='F').tofile(os.path.join(wd, name))


def run(tag, mirror):
    wd = os.path.join(tmp, tag)
    os.makedirs(os.path.join(wd, 'restart'))
    text = open(os.path.join(src, 'input_parameters')).read()
    if mirror:
        text = text.replace('Utarget = 1.0', 'Utarget = -1.0')
    open(os.path.join(wd, 'input_parameters'), 'w').write(text)
    fields(wd, mirror)
    r = subprocess.run([a.mpirun, '--oversubscribe', '-np', '1', a.exe], cwd=wd, capture_output=True, text=True, timeout=900)
    if r.returncode != 0:
        print(r.stdout[-2000:], r.stderr[-2000:]); sys.exit(1)
    rows = [l.split() for l in r.stdout.splitlines() if re.match(r'^\s*\d+\s+[-\d.E+]+\s+[-\d.E+]+\s+[-\d.E+]+\s', l)]
    return abs(float(rows[-1][2])), float(rows[-1][3])   # mean |U|, max |U|


try:
    (u0, m0), (u1, m1) = run('fwd', False), run('mirror', True)
    ok = abs(u0 - u1) < 1e-4 and abs(m0 - m1) < 1e-3*m0 and u0 < 0.99
    print('mean U %.6f / %.6f  max |U| %.5f / %.5f  %s' % (u0, u1, m0, m1, 'ok' if ok else 'FAIL'))
    sys.exit(0 if ok else 1)
finally:
    shutil.rmtree(tmp, ignore_errors=True)
