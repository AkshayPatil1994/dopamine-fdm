#!/usr/bin/env python3
"""Figures of docs/Two-Phase-VOF.md from run directories of the two-fluid solver.

  vof_figures.py breaker RUN_DIR OUT.png [t1 t2 ...]   interface and liquid fraction at the snapshot times nearest to t1.. (vof_debug = 1,
                                                       vof_snap_dt > 0; one rank, or the per-rank vof_snap_NNNNN_rRRR.dat files)
"""
import glob, re, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def read_snapshot(run, n):
    files = sorted(glob.glob('%s/vof_snap_%05d.dat' % (run, n)) + glob.glob('%s/vof_snap_%05d_r*.dat' % (run, n)))
    t, rows = None, []
    for f in files:
        for line in open(f):
            if line.startswith('#'):
                t = float(line.split('=')[1])
            else:
                rows.append([float(v) for v in line.split()])
    a = np.array(rows)
    x, y = np.unique(np.round(a[:, 0], 9)), np.unique(np.round(a[:, 1], 9))
    c = np.full((len(y), len(x)), np.nan)
    ix, iy = np.searchsorted(x, np.round(a[:, 0], 9)), np.searchsorted(y, np.round(a[:, 1], 9))
    c[iy, ix] = a[:, 2]
    return t, x, y, c


def breaker(run, out, times):
    nsnap = len(set(re.findall(r'vof_snap_(\d+)', ' '.join(glob.glob(run + '/vof_snap_*.dat')))))
    snaps = [read_snapshot(run, n) for n in range(1, nsnap + 1)]
    ts = np.array([s[0] for s in snaps])
    pick = [int(np.argmin(abs(ts - t))) for t in times] if times else list(range(0, nsnap, max(1, nsnap // 6)))[:6]
    fig, ax = plt.subplots(2, (len(pick) + 1) // 2, figsize=(11, 6.2), sharex=True, sharey=True)
    for a, i in zip(ax.ravel(), pick):
        t, x, y, c = snaps[i]
        a.pcolormesh(x, y, c, cmap='Blues', vmin=0, vmax=1.3, shading='auto', rasterized=True)
        a.contour(x, y, c, levels=[0.5], colors='k', linewidths=0.8)
        a.set_title('t = %.2f s' % t, fontsize=10)
        a.set_aspect('equal'); a.set_ylim(0.2, 0.9)
    for a in ax[-1]:
        a.set_xlabel('x / $\\lambda$')
    for a in ax[:, 0]:
        a.set_ylabel('y / $\\lambda$')
    fig.tight_layout()
    fig.savefig(out, dpi=130)


breaker(sys.argv[2], sys.argv[3], [float(v) for v in sys.argv[4:]])
