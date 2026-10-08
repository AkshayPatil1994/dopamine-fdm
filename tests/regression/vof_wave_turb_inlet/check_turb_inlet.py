#!/usr/bin/env python3
"""Recycled-precursor turbulence at the wave inlet (wave_turb = 1): a synthetic donor slice with spanwise-periodic fluctuations u', w' of amplitude 0.05
over the whole depth (the air included) is added to the wave + current inlet. The water receives it (max|W| close to the amplitude, none without wave_turb),
the liquid ledger still closes (1e-4)."""
import argparse, os, re, shutil, subprocess, sys, tempfile
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('--exe', required=True)
ap.add_argument('--mpirun', default='mpirun')
a = ap.parse_args()
src = os.path.dirname(os.path.abspath(__file__))
tmp = tempfile.mkdtemp(prefix='voftrb_')
n1, n1v, n2, nsn, amp = 32, 33, 5, 4, 0.05


def donor(wd):
    k = np.arange(n2) % (n2 - 1)   # the periodic last cell repeats the first
    u = np.tile(0.2 + amp*np.cos(2*np.pi*(k + 0.5)/(n2 - 1)), (n1, 1))
    w = np.tile(amp*np.sin(2*np.pi*k/(n2 - 1)), (n1, 1))
    v = np.zeros((n1v, n2))
    with open(os.path.join(wd, 'donor.bin'), 'wb') as f:
        for _ in range(nsn):
            for blk in (u, v, w):
                blk.T.copy().astype('>f8').tofile(f)   # big-endian (the build converts), j (first index) fastest
    np.linspace(0, 3, nsn).astype('>f8').tofile(os.path.join(wd, 'donor_times.bin'))
    open(os.path.join(wd, 'donor_meta.txt'), 'w').write('ncomp = 3\nn1 = %d\nn1_V = %d\nn2 = %d\ncomps = UVW\nnsnaps = %d\ndir = x\n' % (n1, n1v, n2, nsn))


def run(tag, edits):
    wd = os.path.join(tmp, tag)
    os.makedirs(os.path.join(wd, 'restart'))
    text = open(os.path.join(src, 'input_parameters')).read()
    for pat, rep in edits:
        text = re.sub(pat, rep, text, count=1)
    open(os.path.join(wd, 'input_parameters'), 'w').write(text)
    donor(wd)
    r = subprocess.run([a.mpirun, '--oversubscribe', '-np', '1', a.exe], cwd=wd, capture_output=True, text=True, timeout=900)
    if r.returncode != 0:
        print(r.stdout[-2000:], r.stderr[-2000:]); sys.exit(1)
    d = np.loadtxt(os.path.join(wd, 'vof_diag.dat'))
    return d[:, 17].max(), (d[-1, 3] - d[0, 3] - d[-1, 21] - d[-1, 22]) / d[0, 3], d[-1, 1]


try:
    w_on, res, t_end = run('on', [])
    w_off, res_off, _ = run('off', [(r'wave_turb = 1,', 'wave_turb = 0,')])
    ok = t_end > 1.49 and abs(res) < 1e-4 and 0.5*amp < w_on < 2.5*amp and w_off < 0.1*amp
    print('max|W| with %.4f without %.4f (amplitude %.2f)  ledger residual %.2e (%.2e without)  %s' % (w_on, w_off, amp, res, res_off, 'ok' if ok else 'FAIL'))
    sys.exit(0 if ok else 1)
finally:
    shutil.rmtree(tmp, ignore_errors=True)
