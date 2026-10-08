#!/usr/bin/env python3
"""SDF cells at exactly phi = 0 are fluid everywhere (solid is phi < 0): a slab of phi = 0 cells in a noisy stream gives the same
flow as the same slab at phi = +1e-9 (no eddy viscosity zeroed in it)."""
import argparse, os, shutil, subprocess, sys, tempfile
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('--exe', required=True)
ap.add_argument('--mpirun', default='mpirun')
a = ap.parse_args()
src = os.path.dirname(os.path.abspath(__file__))
tmp = tempfile.mkdtemp(prefix='phizero_')
NX, NY, NZ, H = 18, 33, 18, 1/32   # SDF file shape (NX+1, NY+1, NZ-1)


def run(tag, np_, val):
    wd = os.path.join(tmp, tag)
    os.makedirs(os.path.join(wd, 'restart'))
    text = open(os.path.join(src, 'input_parameters')).read()
    open(os.path.join(wd, 'input_parameters'), 'w').write(text)
    phi = np.full((NX + 1, NY + 1, NZ - 1), H/2)
    phi[9:12, 10:20, :] = val
    phi.astype('>f8').flatten(order='F').tofile(os.path.join(wd, 'SDF_in'))
    r = subprocess.run([a.mpirun, '--oversubscribe', '-np', str(np_), a.exe], cwd=wd, capture_output=True, text=True, timeout=900)
    if r.returncode != 0:
        print(r.stdout[-2000:], r.stderr[-2000:]); sys.exit(1)
    b = open(os.path.join(wd, 'fields', 'x.20'), 'rb').read()
    return np.frombuffer(b[-8*(NX + 1)*(NY + 1)*(NZ - 1):], dtype='>f8')   # nu_t, the last block of the snapshot


try:
    z, p = run('zero', 1, 0.0), run('pos', 1, 1e-9)
    d = np.abs(z - p).max()
    print('max difference of nu_t between phi = 0 and phi = 1e-9 slabs: %.2e (max nu_t %.2e)' % (d, z.max()))
    ok = d < 1e-10 * z.max()
    print('ok' if ok else 'FAIL')
    sys.exit(0 if ok else 1)
finally:
    shutil.rmtree(tmp, ignore_errors=True)
