#!/usr/bin/env python3
"""Elevation-record wave inlet (wave_type = 4): the FFT decomposition reproduces the input record (fields/wave_series_check.out), the packet
reaches the gauges at x = 1.0 and 1.5 with the record amplitude (0.01) at the finite-depth group speed (+/-0.4 s), and the liquid ledger closes."""
import argparse, os, shutil, subprocess, sys, tempfile
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('--exe', required=True)
ap.add_argument('--mpirun', default='mpirun')
a = ap.parse_args()
src = os.path.dirname(os.path.abspath(__file__))
tmp = tempfile.mkdtemp(prefix='vofser_')
amp, depth, f0 = 0.01, 0.5, 1.0


def wave_number(w, d, g=9.81):
    k = w*w/g
    for _ in range(60):
        k -= (g*k*np.tanh(k*d) - w*w)/(g*np.tanh(k*d) + g*k*d/np.cosh(k*d)**2)
    return k


def group_speed(f, d, dw=1e-4):
    w = 2*np.pi*f
    return 2*dw/(wave_number(w + dw, d) - wave_number(w - dw, d))


cg = group_speed(f0, depth)
try:
    os.makedirs(os.path.join(tmp, 'restart'))
    shutil.copy(os.path.join(src, 'input_parameters'), tmp)
    t = np.arange(0, 6.0001, 0.05)
    np.savetxt(os.path.join(tmp, 'series.txt'), np.c_[t, amp*np.exp(-((t - 3.0)/1.0)**2)*np.cos(2*np.pi*(t - 3.0))], header='t eta')
    r = subprocess.run([a.mpirun, '--oversubscribe', '-np', '1', a.exe], cwd=tmp, capture_output=True, text=True, timeout=900)
    if r.returncode != 0:
        print(r.stdout[-2000:], r.stderr[-2000:]); sys.exit(1)
    hdr = open(os.path.join(tmp, 'fields', 'wave_series_check.out')).read().splitlines()[3].split()
    rms_err = float(hdr[6])
    d = np.loadtxt(os.path.join(tmp, 'vof_diag.dat'))
    resid = (d[-1, 3] - d[0, 3] - d[-1, 21] - d[-1, 22])/d[0, 3]
    g = np.loadtxt(os.path.join(tmp, 'vof_gauges.dat'))
    ok = rms_err < 0.02 and d[-1, 1] > 5.99 and abs(resid) < 1e-5
    print('t_end %.2f  ledger residual %.2e  reconstruction rms error %.2e' % (d[-1, 1], resid, rms_err))
    for col, x in ((1, 1.0), (2, 1.5)):
        e = g[:, col] - g[0, col]
        tpk = g[np.argmax(np.abs(e)), 0]
        t_exp = 3.0 + x/cg
        good = 0.7*amp < np.abs(e).max() < 1.3*amp and abs(tpk - t_exp) < 0.4
        ok = ok and good
        print('gauge x=%.1f: max |eta| %.4f (record %.3f), peak at t=%.2f (group-speed estimate %.2f)  %s'
              % (x, np.abs(e).max(), amp, tpk, t_exp, 'ok' if good else 'FAIL'))
    sys.exit(0 if ok else 1)
finally:
    shutil.rmtree(tmp, ignore_errors=True)
