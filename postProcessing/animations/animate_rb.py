#!/usr/bin/env python3
"""Rayleigh-Benard convection showcase: temperature on a near-wall plan view (constant-y slice)."""
import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from anim_common import Snapshot, list_steps, step_time_map, render

ap = argparse.ArgumentParser()
ap.add_argument("--case-dir", default="/mnt/storage1/fdm-dopamine/validation/RBconvection")
ap.add_argument("--prefix", default="RBconv_Kunnen")
ap.add_argument("--out", default="RBconvection.gif")
ap.add_argument("--every", type=int, default=6, help="use every N-th snapshot")
ap.add_argument("--y-plan", type=float, default=0.08, help="height of the plan-view slice")
ap.add_argument("--fps", type=float, default=3.5)
a = ap.parse_args()

case = Path(a.case_dir)
steps = list_steps(case / "fields", a.prefix)[::a.every]
ss, tt = step_time_map([case / "run.log"])
first = Snapshot(steps[0][1], boussinesq_flag=1)
x, y, z, xm, ym, zm = first.coords()
jy = int(np.argmin(abs(ym - a.y_plan)))
Lx, Lz = x[-1], z[-1]


def frame(i):
    s, p = steps[i]
    snap = Snapshot(p, boussinesq_flag=1)
    plan = snap.plane("T", 1, jy).T
    fig = plt.figure(figsize=(5.6, 5.6 * Lz / Lx + 0.9))
    ax = fig.add_axes([0.02, 0.01, 0.96, 0.90])
    ax.pcolormesh(xm, zm, plan, cmap="RdYlBu_r", vmin=0, vmax=1, shading="gouraud", rasterized=True)
    ax.set_aspect("equal")
    ax.set_axis_off()
    ax.text(0.01, 0.98, f"plan view, $y={ym[jy]:.2f}$", transform=ax.transAxes, fontsize=9, va="top",
            bbox=dict(boxstyle="round,pad=0.25", fc="#0d1117", ec="none", alpha=0.7))
    t = np.interp(s, ss, tt)
    fig.suptitle(f"Rayleigh–Bénard convection (DNS) — temperature     $t={t:6.1f}$", fontsize=11.5, y=0.975)
    return fig


render(frame, len(steps), a.out, fps=a.fps, width=440)
