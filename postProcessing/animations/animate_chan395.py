#!/usr/bin/env python3
"""LES of channel flow at Re_tau=395 with wall model: near-wall streaks and a vertical slice."""
import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from anim_common import Snapshot, list_steps, step_time_map, render

ap = argparse.ArgumentParser()
ap.add_argument("--case-dir", default="/mnt/storage1/fdm-dopamine/validation/chan395")
ap.add_argument("--prefix", default="chan395")
ap.add_argument("--out", default="chan395.gif")
ap.add_argument("--every", type=int, default=1)
ap.add_argument("--yplus", type=float, default=40.0, help="wall distance (in wall units) of the plan view")
ap.add_argument("--fps", type=float, default=3.5)
ap.add_argument("--nstart", type=int, default=20000, help="skip snapshots before this step (spin-up)")
a = ap.parse_args()

case = Path(a.case_dir)
steps = [s for s in list_steps(case / "fields", a.prefix) if s[0] >= a.nstart][::a.every]
ss, tt = step_time_map([case / "run.log"])
x, y, z, xm, ym, zm = Snapshot(steps[0][1], sgs_model=1).coords()
Re_tau = 395.0
jy = int(np.argmin(abs(ym - a.yplus / Re_tau)))
kz = len(zm) // 2
Lx, Ly, Lz = x[-1], y[-1], z[-1]


def frame(i):
    s, p = steps[i]
    snap = Snapshot(p, sgs_model=1)
    plan = snap.plane("U", 1, jy)
    plan = (plan - plan.mean()).T
    side = snap.plane("U", 2, kz).T
    fig = plt.figure(figsize=(7.4, 7.4 * 0.96 * (Ly + Lz) / Lx * 1.02 + 0.9))
    gs = fig.add_gridspec(2, 1, height_ratios=[Ly / Lx, Lz / Lx], hspace=0.06,
                          left=0.02, right=0.98, top=0.90, bottom=0.02)
    panels = [(side, ym, "streamwise velocity $u$, vertical slice", "inferno", 2.0, 22.0),
              (plan, zm, f"$u'$, plan view, $y^+\\approx{ym[jy] * Re_tau:.0f}$", "RdBu_r", -4.0, 4.0)]
    for k, (d, v, lab, cm, lo, hi) in enumerate(panels):
        ax = fig.add_subplot(gs[k])
        ax.pcolormesh(xm, v, d, cmap=cm, vmin=lo, vmax=hi, shading="gouraud", rasterized=True)
        ax.set_aspect("equal")
        ax.set_axis_off()
        ax.text(0.005, 0.97, lab, transform=ax.transAxes, fontsize=9, va="top",
                bbox=dict(boxstyle="round,pad=0.25", fc="#0d1117", ec="none", alpha=0.7))
    t = np.interp(s, ss, tt)
    fig.suptitle(f"Wall-modelled LES, channel $Re_\\tau=395$     $t={t:6.1f}$", fontsize=11.5, y=0.975)
    return fig


render(frame, len(steps), a.out, fps=a.fps, width=520)
