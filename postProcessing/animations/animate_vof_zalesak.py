#!/usr/bin/env python3
"""Zalesak's slotted disk (PLIC-VOF, one revolution of rigid rotation): liquid fraction against the exact (rotated) shape."""
import argparse
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

sys.path.insert(0, str(Path(__file__).parent))
from anim_common import read_vof_snapshot, render

XC, YC, R, W, YTOP, CEN = 0.5, 0.75, 0.15, 0.05, 0.85, (0.5, 0.5)

ap = argparse.ArgumentParser()
ap.add_argument("--case-dir", required=True, help="run directory with vof_snap_*.dat")
ap.add_argument("--out", default="vof_zalesak.gif")
ap.add_argument("--every", type=int, default=1)
ap.add_argument("--period", type=float, default=1.0, help="rotation period T")
ap.add_argument("--fps", type=float, default=12)
a = ap.parse_args()

run = a.case_dir
ids = sorted({int(m) for f in Path(run).glob("vof_snap_*.dat") for m in re.findall(r"vof_snap_(\d+)", f.name)})[::a.every]
cmap = LinearSegmentedColormap.from_list("liq", ["#0d1117", "#1f6fd0", "#8fd0ff"])


def exact_outline(t):
    th0 = np.arcsin(W / 2 / R)
    th = np.linspace(-0.5 * np.pi + th0, 1.5 * np.pi - th0, 400)
    px = np.r_[XC + R * np.cos(th), XC - W / 2, XC + W / 2, XC + R * np.cos(th[0])]
    py = np.r_[YC + R * np.sin(th), YTOP, YTOP, YC + R * np.sin(th[0])]
    ang = 2 * np.pi * t / a.period
    dx, dy = px - CEN[0], py - CEN[1]
    return CEN[0] + dx * np.cos(ang) - dy * np.sin(ang), CEN[1] + dx * np.sin(ang) + dy * np.cos(ang)


def frame(i):
    t, x, y, c = read_vof_snapshot(run, ids[i])
    kx, ky = x < 1, y < 1
    x, y, c = x[kx], y[ky], c[np.ix_(ky, kx)]
    fig = plt.figure(figsize=(4.6, 4.9))
    ax = fig.add_axes([0.03, 0.02, 0.94, 0.88])
    h = x[1] - x[0]
    ax.imshow(np.clip(c, 0, 1), cmap=cmap, vmin=0, vmax=1, origin="lower", interpolation="bicubic",
              extent=(x[0] - h / 2, x[-1] + h / 2, y[0] - h / 2, y[-1] + h / 2))
    ax.contour(x, y, c, levels=[0.5], colors="#e6edf3", linewidths=1.0)
    ax.plot(*exact_outline(t), color="#ffb454", lw=1.0, ls="--")
    ax.set_aspect("equal")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_axis_off()
    fig.suptitle(f"Zalesak's slotted disk, 100² cells (dashed: exact)     $t/T={t / a.period:4.2f}$", fontsize=10, y=0.97)
    return fig


render(frame, len(ids), a.out, fps=a.fps, width=560, dpi=220, colors=256, dither="sierra2_4a")
