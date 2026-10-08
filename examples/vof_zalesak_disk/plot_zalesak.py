#!/usr/bin/env python3
"""Zalesak slotted-disk validation: one figure comparing the PLIC-VOF solution with the exact solution.

  plot_zalesak.py [RUN_DIR] [OUT.png]      RUN_DIR holds vof_diag.dat and vof_snap_NNNNN.dat (vof_debug = 1, vof_snap_dt > 0)

The exact solution is the initial shape rigidly rotated (counter-clockwise about (0.5, 0.5), period T) and is drawn as the
reference at every snapshot time; the numerical interface is the C = 0.5 contour.
"""
import glob, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

XC, YC, R, W, YTOP, T = 0.5, 0.75, 0.15, 0.05, 0.85, 1.0   # disk centre, radius, slot width, slot top, period
OMEGA_C = (0.5, 0.5)


def read_snapshot(f):
    t, rows = None, []
    for line in open(f):
        if line.startswith('#'):
            t = float(line.split('=')[1])
        else:
            rows.append([float(v) for v in line.split()])
    a = np.array(rows)
    x, y = np.unique(np.round(a[:, 0], 9)), np.unique(np.round(a[:, 1], 9))
    c = np.full((len(y), len(x)), np.nan)
    c[np.searchsorted(y, np.round(a[:, 1], 9)), np.searchsorted(x, np.round(a[:, 0], 9))] = a[:, 2]
    return t, x, y, c


def exact_outline(t):
    """Closed polygon of the slotted disk rotated by 2 pi t/T about the domain centre."""
    th0 = np.arcsin(W/2/R)
    th = np.linspace(-0.5*np.pi + th0, 1.5*np.pi - th0, 400)
    px = np.r_[XC + R*np.cos(th), XC - W/2, XC + W/2, XC + R*np.cos(th[0])]
    py = np.r_[YC + R*np.sin(th), YTOP, YTOP, YC + R*np.sin(th[0])]
    a = 2*np.pi*t/T
    dx, dy = px - OMEGA_C[0], py - OMEGA_C[1]
    return OMEGA_C[0] + dx*np.cos(a) - dy*np.sin(a), OMEGA_C[1] + dx*np.sin(a) + dy*np.cos(a)


def exact_fraction(x, y, ns=16):
    """Cell-averaged exact initial shape on the cell-centre grids x, y (uniform spacing, ns x ns sub-samples)."""
    h = x[1] - x[0]
    o = (np.arange(ns) + 0.5)/ns - 0.5
    xs = (x[:, None] + h*o).ravel()
    ys = (y[:, None] + h*o).ravel()
    X, Y = np.meshgrid(xs, ys)
    inside = ((X - XC)**2 + (Y - YC)**2 <= R**2) & ~((np.abs(X - XC) <= W/2) & (Y <= YTOP))
    return inside.reshape(len(y), ns, len(x), ns).mean(axis=(1, 3))


def main():
    run = sys.argv[1] if len(sys.argv) > 1 else '.'
    out = sys.argv[2] if len(sys.argv) > 2 else 'zalesak_validation.png'
    snaps = [read_snapshot(f) for f in sorted(glob.glob(run + '/vof_snap_*.dat'))]
    snaps = [s for s in snaps if s[0] > 0.2]     # the first file is written at the first monitor step, not at t = 0
    diag = np.loadtxt(run + '/vof_diag.dat')

    t, x, y, c = snaps[-1]
    keep_x, keep_y = x < 1, y < 1                # drop the periodic duplicate cells outside the unit box
    x, y, c = x[keep_x], y[keep_y], c[np.ix_(keep_y, keep_x)]
    h = x[1] - x[0]
    c0 = exact_fraction(x, y)
    l1 = np.abs(c - c0).sum()*h*h
    area = c0.sum()*h*h
    drift = diag[-1, 4]

    col = ['#5b8def', '#e69f00', '#2a9d6f']      # t = 0.25, 0.5, 0.75 (needs vof_snap_dt = 0.25)
    fig = plt.figure(figsize=(14, 4.9), constrained_layout=True)
    ax = fig.subplots(1, 3, gridspec_kw={'width_ratios': [1, 1, 1.15]})

    # (a) the rotation: interface at several times against the exact solution
    ax[0].plot(*exact_outline(0), color='0.45', lw=1, ls=':', label='initial')
    for (ts, xs_, ys_, cs), k in zip(snaps[:-1], col):
        kx, ky = xs_ < 1, ys_ < 1
        ax[0].contour(xs_[kx], ys_[ky], cs[np.ix_(ky, kx)], [0.5], colors=k, linewidths=1.6)
        ax[0].plot(*exact_outline(ts), color='k', lw=0.9, ls='--')
        ax[0].plot([], [], color=k, lw=1.6, label='VOF  t/T = %.2f' % (ts/T))
    ax[0].contour(x, y, c, [0.5], colors='crimson', linewidths=1.6)
    ax[0].plot([], [], color='crimson', lw=1.6, label='VOF  t/T = %.2f' % (t/T))
    ax[0].plot([], [], color='k', lw=0.9, ls='--', label='exact')
    ax[0].set_xlim(0, 1);  ax[0].set_ylim(0, 1);  ax[0].set_aspect('equal')
    ax[0].set_title('(a) interface during one revolution')
    ax[0].legend(loc='lower left', fontsize=7.5, frameon=True)

    # (b) end of the revolution: liquid fraction and interface against the exact (initial) outline
    z = (x >= 0.3) & (x <= 0.7);  zy = (y >= 0.55) & (y <= 0.95)
    m = ax[1].pcolormesh(x[z] - h/2, y[zy] - h/2, c[np.ix_(zy, z)], cmap='Blues', vmin=0, vmax=1, shading='nearest')
    ax[1].contour(x[z], y[zy], c[np.ix_(zy, z)], [0.5], colors='crimson', linewidths=1.8)
    ax[1].plot(*exact_outline(0), color='k', lw=1.4, ls='--')
    ax[1].plot([], [], color='crimson', lw=1.8, label='VOF  C = 0.5, t = T')
    ax[1].plot([], [], color='k', lw=1.4, ls='--', label='exact (= initial shape)')
    ax[1].set_xlim(0.34, 0.66);  ax[1].set_ylim(0.58, 0.92);  ax[1].set_aspect('equal')
    ax[1].set_title('(b) after one revolution (zoom)')
    ax[1].legend(loc='lower right', fontsize=7.5, frameon=True)
    fig.colorbar(m, ax=ax[1], shrink=0.6, label='C (VOF)')

    # (c) liquid-fraction profile across the slot, y = 0.8, at t = T
    j = np.argmin(np.abs(y - 0.8))
    ax[2].step(x, c0[j], where='mid', color='k', ls='--', lw=1.4, label='exact (cell average)')
    ax[2].plot(x, c[j], 'o-', color='crimson', ms=3, lw=1.2, label='VOF')
    ax[2].set_xlim(0.33, 0.67)
    ax[2].set_xlabel('x');  ax[2].set_ylabel('C')
    ax[2].set_title('(c) profile at y = %.3f, t = T' % y[j])
    ax[2].legend(loc='center right', fontsize=8, bbox_to_anchor=(1, 0.45))
    ax[2].set_ylim(-0.05, 1.6)
    ax[2].text(0.03, 0.97, 'L1 error = %.2e (area units)\n= %.2f %% of the disk area\nvolume drift = %.1e\nC range [%.1f, %.1f]'
               % (l1, 100*l1/area, drift, np.nanmin(c), np.nanmax(c)), transform=ax[2].transAxes, va='top', fontsize=8.5,
               bbox=dict(boxstyle='round', fc='white', ec='0.7'))

    fig.suptitle("Zalesak's slotted disk, PLIC-VOF, %d x %d cells (R = 0.15, slot 0.05, one revolution)" % (len(x), len(y)))
    fig.savefig(out, dpi=170)
    print('L1 = %.3e  (%.2f %% of area)  volume drift = %.2e  -> %s' % (l1, 100*l1/area, drift, out))


if __name__ == '__main__':
    main()
