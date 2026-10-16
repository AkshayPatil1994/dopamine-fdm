#!/usr/bin/env python3
"""Wave statistics of a TK94 run from vof_gauges.dat: per gauge the mean water level (set-up/set-down, relative to the still-water
level), crest and trough levels of the phase-averaged wave and the wave height, over whole periods after the spin-up
(time-weighted mean level, since the solver time step varies).

usage: postprocess.py CASE_DIR [--t0 T_START] [--out table.dat] [--exp exp/gauges.dat]
--exp: optional experimental table with columns x_TK  setup  crest  trough [m] (digitised TK94 data; not shipped) for comparison.
"""
import argparse, os, re
import numpy as np
from wavestats import namelist_values, resolve_restarts, phase_stats

ap = argparse.ArgumentParser()
ap.add_argument('case')
ap.add_argument('--t0', type=float, default=None, help='start of averaging window; trailing partial wave is excluded')
ap.add_argument('--out'); ap.add_argument('--exp')
ap.add_argument('--x0', type=float, default=None, help='solver x of x_TK = 0 (default: X0 from stations.dat)')
a = ap.parse_args()

X0 = a.x0                # solver x = x_TK + X0; vof_gauges.dat holds eta relative to the still-water level
if X0 is None:
    st = os.path.join(a.case, 'stations.dat')
    m = re.search(r'\bX0\s+([-\d.eE+]+)', open(st).read()) if os.path.exists(st) else None
    if not m:
        raise SystemExit('X0 not found: stations.dat missing or without a "# case ... X0" header; pass --x0')
    X0 = float(m.group(1))
inp = open(os.path.join(a.case, 'input_parameters')).read()
gx = namelist_values(inp, 'wave_gauge_x')
T_WAVE = (namelist_values(inp, 'wave_period') or [None])[0]
if not gx or not T_WAVE:
    raise SystemExit('wave_gauge_x / wave_period not found in input_parameters')
d = np.loadtxt(os.path.join(a.case, 'vof_gauges.dat'))
if d.shape[1] != 1 + len(gx):
    raise SystemExit(f'vof_gauges.dat has {d.shape[1] - 1} gauge columns, input_parameters lists {len(gx)} gauges')
kept = resolve_restarts(d[:, 0])
if len(kept) < len(d):
    print(f'WARNING: vof_gauges.dat contains restarts; dropped {len(d) - len(kept)} stale of {len(d)} rows')
t, eta = d[kept, 0], d[kept, 1:]
t1 = t[-1]
t0 = a.t0 if a.t0 is not None else 0.5*t1       # default: last half of the run
rows = []
for j, x in enumerate(gx):
    s = phase_stats(t, eta[:, j], t0, t1, T_WAVE)
    rows.append((x - X0, s['mean'], s['crest'], s['trough'], s['H']))
nper, t_end = s['nper'], s['t_end']
rows = np.array(rows)
hdr = 'x_TK[m]  mean_level[m]  crest[m]  trough[m]  H[m]   (%d periods, t = %.1f..%.1f s)' % (nper, t0, t_end)
print(hdr)
for r in rows:
    print('%8.2f %12.4f %10.4f %10.4f %8.4f' % tuple(r))
if a.out:
    np.savetxt(a.out, rows, header=hdr)
if a.exp and os.path.exists(a.exp):
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    ex = np.loadtxt(a.exp)
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    ax[0].plot(rows[:, 0], rows[:, 2], 'b-o', label='crest'); ax[0].plot(rows[:, 0], rows[:, 3], 'r-o', label='trough')
    ax[0].plot(ex[:, 0], ex[:, 2], 'bs', mfc='none'); ax[0].plot(ex[:, 0], ex[:, 3], 'rs', mfc='none')
    ax[1].plot(rows[:, 0], rows[:, 1], 'k-o', label='LES'); ax[1].plot(ex[:, 0], ex[:, 1], 'ks', mfc='none', label='TK94')
    for q, tl in zip(ax, ('crest / trough', 'mean level (set-up)')):
        q.set_xlabel('x_TK [m]'); q.set_ylabel('m'); q.set_title(tl); q.legend()
    fig.tight_layout(); fig.savefig(os.path.join(a.case, 'tk94_gauges.png'), dpi=130)
