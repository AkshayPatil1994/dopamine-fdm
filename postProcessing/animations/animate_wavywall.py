#!/usr/bin/env python3
"""Turbulent flow over a wavy wall (IBM DNS) showcase: streamwise velocity, vertical slice + plan view."""
import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from anim_common import Snapshot, list_steps, step_time_map, render

ap = argparse.ArgumentParser()
ap.add_argument("--case-dir", default="/mnt/storage1/fdm-dopamine/validation/waveWall")
ap.add_argument("--prefix", default="channel")
ap.add_argument("--out", default="waveWall.gif")
ap.add_argument("--every", type=int, default=4)
ap.add_argument("--y-plan", type=float, default=0.2, help="height of the plan-view slice")
ap.add_argument("--fps", type=float, default=3.5)
a = ap.parse_args()

case = Path(a.case_dir)
steps = list_steps(case / "fields", a.prefix)[::a.every]
ss, tt = step_time_map([case / "run.log"])
x, y, z, xm, ym, zm = Snapshot(steps[0][1]).coords()
nxm, nym, nzm = len(xm), len(ym), len(zm)
phi = np.fromfile(case / "SDF_in", dtype=">f8").reshape((nxm + 2, nym + 2, nzm), order="F")[1:-1, 1:-1, :]
kz = nzm // 2
jy = int(np.argmin(abs(ym - a.y_plan)))
solid_side = phi[:, :, kz] < 0
solid_plan = phi[:, jy, :] < 0
Lx, Ly, Lz = x[-1], y[-1], z[-1]
cmap = plt.get_cmap("magma_r").copy()
cmap.set_bad("#3b4048")
VMIN, VMAX = -0.3, 1.4


def frame(i):
    s, p = steps[i]
    snap = Snapshot(p)
    side = np.ma.masked_where(solid_side, snap.plane("U", 2, kz)).T
    plan = np.ma.masked_where(solid_plan, snap.plane("U", 1, jy)).T
    fig = plt.figure(figsize=(7.4, 7.4 * 0.96 * (Ly + Lz) / Lx * 1.02 + 0.9))
    gs = fig.add_gridspec(2, 1, height_ratios=[Ly / Lx, Lz / Lx], hspace=0.06,
                          left=0.02, right=0.98, top=0.90, bottom=0.02)
    for k, (d, v, lab) in enumerate([(side, ym, "vertical slice, $z=L_z/2$"),
                                     (plan, zm, f"plan view, $y={ym[jy]:.2f}$")]):
        ax = fig.add_subplot(gs[k])
        ax.pcolormesh(xm, v, d, cmap=cmap, vmin=VMIN, vmax=VMAX, shading="gouraud", rasterized=True)
        if k == 0:
            ax.contour(xm, ym, phi[:, :, kz].T, levels=[0], colors="white", linewidths=0.8)
            ax.set_ylim(0, Ly)
        ax.set_aspect("equal")
        ax.set_axis_off()
        ax.text(0.005, 0.97, lab, transform=ax.transAxes, fontsize=9, va="top",
                bbox=dict(boxstyle="round,pad=0.25", fc="#0d1117", ec="none", alpha=0.7))
    t = np.interp(s, ss, tt)
    fig.suptitle(f"Turbulent flow over a wavy wall (IBM DNS) — streamwise velocity $u$     $t={t:6.1f}$",
                 fontsize=11.5, y=0.975)
    return fig


render(frame, len(steps), a.out, fps=a.fps, width=440)
