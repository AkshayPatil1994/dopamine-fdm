#!/usr/bin/env python3
"""Plunging breaking wave (two-fluid PLIC-VOF, 1000:1): liquid fraction and interface from the debug snapshots of a vof_debug = 1 run."""
import argparse
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

sys.path.insert(0, str(Path(__file__).parent))
from anim_common import read_vof_snapshot, render

ap = argparse.ArgumentParser()
ap.add_argument("--case-dir", required=True, help="run directory with vof_snap_*.dat")
ap.add_argument("--out", default="vof_breaker.gif")
ap.add_argument("--every", type=int, default=1, help="use every N-th snapshot")
ap.add_argument("--tmax", type=float, default=1e9)
ap.add_argument("--fps", type=float, default=12)
a = ap.parse_args()

run = a.case_dir
n_all = sorted({int(m) for f in Path(run).glob("vof_snap_*.dat") for m in re.findall(r"vof_snap_(\d+)", f.name)})
ids = [n for n in n_all[::a.every] if read_vof_snapshot(run, n)[0] <= a.tmax]
cmap = LinearSegmentedColormap.from_list("liq", ["#0d1117", "#1f6fd0", "#8fd0ff"])


def frame(i):
    t, x, y, c = read_vof_snapshot(run, ids[i])
    fig = plt.figure(figsize=(6.4, 3.5))
    ax = fig.add_axes([0.02, 0.02, 0.96, 0.86])
    h = x[1] - x[0]
    ax.imshow(np.clip(c, 0, 1), cmap=cmap, vmin=0, vmax=1, origin="lower", interpolation="bicubic",
              extent=(x[0] - h / 2, x[-1] + h / 2, y[0] - h / 2, y[-1] + h / 2))
    ax.contour(x, y, c, levels=[0.5], colors="#e6edf3", linewidths=1.0)
    ax.set_aspect("equal")
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(0.35, 0.85)
    ax.set_axis_off()
    fig.suptitle(f"Plunging breaker, ak = 0.55, 1000:1, 128² cells — liquid fraction     $t={t:4.2f}$ s", fontsize=10.5, y=0.97)
    return fig


render(frame, len(ids), a.out, fps=a.fps, width=720, dpi=220, colors=256, dither="sierra2_4a")
