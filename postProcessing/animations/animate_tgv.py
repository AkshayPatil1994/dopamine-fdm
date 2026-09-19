#!/usr/bin/env python3
"""Taylor-Green vortex (Re=1600) showcase: vorticity magnitude on the z=pi mid-plane."""
import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from anim_common import Snapshot, list_steps, step_time_map, render

ap = argparse.ArgumentParser()
ap.add_argument("--case-dir", default="/mnt/storage1/fdm-dopamine/validation/tgv/512")
ap.add_argument("--prefix", default="tgv")
ap.add_argument("--out", default="tgv.gif")
ap.add_argument("--every", type=int, default=3)
ap.add_argument("--z-frac", type=float, default=0.2, help="plane position as fraction of Lz")
ap.add_argument("--fps", type=float, default=3.5)
a = ap.parse_args()

case = Path(a.case_dir)
steps = list_steps(case / "fields", a.prefix)[::a.every]
ss, tt = step_time_map(sorted(case.glob("run*.log")))
x, y, z, xm, ym, zm = Snapshot(steps[0][1]).coords()
kz = int(round(a.z_frac * len(zm)))
dx, dy, dz = xm[1] - xm[0], ym[1] - ym[0], zm[1] - zm[0]


def ddx(f, ax, h):
    return (np.roll(f, -1, ax) - np.roll(f, 1, ax)) / (2 * h)


def vort_mag(snap):
    up, uc, um = (snap.plane("U", 2, kz + o) for o in (1, 0, -1))
    vp, vc, vm = (snap.plane("V", 2, kz + o) for o in (1, 0, -1))
    wp, wc, wm = (snap.plane("W", 2, kz + o) for o in (1, 0, -1))
    wx = ddx(wc, 1, dy) - (vp - vm) / (2 * dz)
    wy = (up - um) / (2 * dz) - ddx(wc, 0, dx)
    wz = ddx(vc, 0, dx) - ddx(uc, 1, dy)
    return np.sqrt(wx**2 + wy**2 + wz**2)


vmax = np.percentile(vort_mag(Snapshot(steps[int(0.45 * len(steps))][1])), 99.5)
L = x[-1]


def frame(i):
    s, p = steps[i]
    w = vort_mag(Snapshot(p))
    fig = plt.figure(figsize=(6.0, 6.5))
    ax = fig.add_axes([0.02, 0.02, 0.96, 0.88])
    ax.imshow(w.T, origin="lower", extent=[0, L, 0, L], cmap="inferno", vmin=0, vmax=vmax, interpolation="bicubic")
    ax.set_axis_off()
    t = np.interp(s, ss, tt)
    fig.suptitle(f"Taylor–Green vortex, Re = 1600, 512³ DNS — $|\\omega|$ at $z/L_z={a.z_frac:g}$     $t={t:5.2f}$",
                 fontsize=11, y=0.965)
    return fig


render(frame, len(steps), a.out, fps=a.fps, width=440)
