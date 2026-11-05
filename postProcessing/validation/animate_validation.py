#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
"""Light-theme (manuscript style) validation animations: tgv | chan395 | wavywall | rb | particles | urban | particles3d.

Usage: python3 animate_validation.py tgv --out docs/assets/validation/tgv/tgv_vorticity.gif
Frames reuse postProcessing/animations/anim_common.py (plane reader, parallel render, ffmpeg palette GIF).
"""
import argparse
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "animations"))
sys.path.insert(0, str(HERE.parent))
import plotstyle  # noqa: E402
from anim_common import Snapshot, list_steps, render, step_time_map  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

plotstyle.apply()
plt.rcParams.update({"figure.facecolor": "white", "savefig.facecolor": "white", "text.color": "black",
                     "axes.labelcolor": "black", "xtick.color": "black", "ytick.color": "black",
                     "axes.edgecolor": "black", "savefig.bbox": "standard", "savefig.dpi": 110})
V = "/mnt/storage1/fdm-dopamine/validation/"


def tgv(a):
    case = Path(a.case_dir or V + "tgv/256")
    steps = list_steps(case / "fields", "tgv")[::a.every]
    ss, tt = step_time_map(sorted(case.glob("run*.log")))
    x, y, z, xm, ym, zm = Snapshot(steps[0][1]).coords()
    kz = int(round(0.2 * len(zm)))
    dx, dy, dz = xm[1] - xm[0], ym[1] - ym[0], zm[1] - zm[0]
    ddx = lambda f, ax, h: (np.roll(f, -1, ax) - np.roll(f, 1, ax)) / (2 * h)

    def vort(snap):
        up, uc, um = (snap.plane("U", 2, kz + o) for o in (1, 0, -1))
        vp, vc, vm = (snap.plane("V", 2, kz + o) for o in (1, 0, -1))
        wp, wc, wm = (snap.plane("W", 2, kz + o) for o in (1, 0, -1))
        wx = ddx(wc, 1, dy) - (vp - vm) / (2 * dz)
        wy = (up - um) / (2 * dz) - ddx(wc, 0, dx)
        wz = ddx(vc, 0, dx) - ddx(uc, 1, dy)
        return np.sqrt(wx**2 + wy**2 + wz**2)

    vmax = np.percentile(vort(Snapshot(steps[int(0.45 * len(steps))][1])), 99.5)
    L = x[-1]

    def frame(i):
        s, p = steps[i]
        fig, ax = plt.subplots(figsize=(4.6, 4.9), constrained_layout=True)
        im = ax.imshow(vort(Snapshot(p)).T, origin="lower", extent=[0, L, 0, L], cmap="magma_r", vmin=0, vmax=vmax,
                       interpolation="bicubic")
        ax.set_xticks([]); ax.set_yticks([])
        cb = fig.colorbar(im, ax=ax, orientation="horizontal", fraction=0.045, pad=0.03)
        cb.set_label(r"$|\omega|\,L/U_0$")
        ax.set_title(rf"Taylor–Green vortex, $Re=1600$, $256^3$:  $t\,U_0/L={np.interp(s, ss, tt):5.2f}$", fontsize=10)
        return fig
    return steps, frame, 8, 460


def chan395(a):
    case = Path(a.case_dir or V + "chan395")
    steps = [s for s in list_steps(case / "fields", "chan395") if s[0] >= 20000][::a.every]
    ss, tt = step_time_map([case / "run.log"])
    x, y, z, xm, ym, zm = Snapshot(steps[0][1], sgs_model=1).coords()
    jy = int(np.argmin(abs(ym - 40.0 / 395.0)))
    Lx, Lz = x[-1], z[-1]

    def frame(i):
        s, p = steps[i]
        plan = Snapshot(p, sgs_model=1).plane("U", 1, jy)
        fig, ax = plt.subplots(figsize=(6.4, 6.4 * Lz / Lx + 1.3), constrained_layout=True)
        im = ax.imshow((plan - plan.mean()).T, origin="lower", extent=[0, Lx, 0, Lz], cmap="RdBu_r", vmin=-4, vmax=4,
                       interpolation="bicubic")
        ax.set_xlabel(r"$x/\delta$"); ax.set_ylabel(r"$z/\delta$")
        cb = fig.colorbar(im, ax=ax, orientation="horizontal", fraction=0.06, pad=0.05)
        cb.set_label(r"$u'/u_\tau$")
        ax.set_title(rf"Wall-modelled LES, $Re_\tau=395$, $y^+\approx40$:  $t\,u_\tau/\delta={np.interp(s, ss, tt):5.1f}$",
                     fontsize=10)
        return fig
    return steps, frame, 8, 520


def wavywall(a):
    case = Path(a.case_dir or V + "waveWall")
    steps = list_steps(case / "fields", "channel")[::a.every]
    ss, tt = step_time_map([case / "run.log"])
    x, y, z, xm, ym, zm = Snapshot(steps[0][1]).coords()
    nxm, nym, nzm = len(xm), len(ym), len(zm)
    phi = np.fromfile(case / "SDF_in", dtype=">f8").reshape((nxm + 2, nym + 2, nzm), order="F")[1:-1, 1:-1, :]
    kz = nzm // 2
    solid = phi[:, :, kz] < 0
    Lx, Ly = x[-1], y[-1]
    cmap = plt.get_cmap("magma_r").copy()
    cmap.set_bad("#9a9a9a")

    def frame(i):
        s, p = steps[i]
        d = np.ma.masked_where(solid, Snapshot(p).plane("U", 2, kz)).T
        fig, ax = plt.subplots(figsize=(7.0, 7.0 * Ly / Lx + 1.4), constrained_layout=True)
        im = ax.pcolormesh(xm, ym, d, cmap=cmap, vmin=-0.3, vmax=1.4, shading="gouraud", rasterized=True)
        ax.contour(xm, ym, phi[:, :, kz].T, levels=[0], colors="k", linewidths=0.8)
        ax.set_aspect("equal"); ax.set_ylim(0, Ly)
        ax.set_xlabel(r"$x/H$"); ax.set_ylabel(r"$y/H$")
        cb = fig.colorbar(im, ax=ax, orientation="horizontal", fraction=0.06, pad=0.05)
        cb.set_label(r"$u/U_b$")
        ax.set_title(rf"Flow over a wavy wall (IBM DNS), $Re_H=6760$:  $t={np.interp(s, ss, tt):6.1f}$", fontsize=10)
        return fig
    return steps, frame, 8, 560


def rb(a):
    case = Path(a.case_dir or V + "RBconvection")
    steps = list_steps(case / "fields", "RBconv_Kunnen")
    steps = steps[len(steps) // 2::a.every]
    ss, tt = step_time_map([case / "run.log"])
    x, y, z, xm, ym, zm = Snapshot(steps[0][1], boussinesq_flag=1).coords()
    jy = int(np.argmin(abs(ym - 0.08)))
    kz = len(zm) // 2
    Lx, Ly, Lz = x[-1], y[-1], z[-1]

    def frame(i):
        s, p = steps[i]
        snap = Snapshot(p, boussinesq_flag=1)
        fig = plt.figure(figsize=(6.0, 6.0 * (Ly + Lz) / Lx + 1.2), constrained_layout=True)
        gs = fig.add_gridspec(2, 1, height_ratios=[Ly / Lx, Lz / Lx])
        for k, (d, ext, lab, yl) in enumerate([
                (snap.plane("T", 2, kz).T, [0, Lx, 0, Ly], r"vertical slice, $z=L_z/2$", r"$y/H$"),
                (snap.plane("T", 1, jy).T, [0, Lx, 0, Lz], rf"plan view, $y/H={ym[jy]:.2f}$", r"$z/H$")]):
            ax = fig.add_subplot(gs[k])
            im = ax.imshow(d, origin="lower", extent=ext, cmap="RdYlBu_r", vmin=0, vmax=1, interpolation="bicubic")
            ax.set_ylabel(yl)
            ax.set_title(lab, fontsize=9, loc="left")
            if k == 1:
                ax.set_xlabel(r"$x/H$")
        cb = fig.colorbar(im, ax=fig.axes, orientation="horizontal", fraction=0.05, pad=0.03)
        cb.set_label(r"$T$")
        fig.suptitle(rf"Rayleigh–Bénard convection, $Ra=2.5\times10^6$:  $t={np.interp(s, ss, tt):6.1f}$", fontsize=10)
        return fig
    return steps, frame, 6, 480


def urban(a):
    run = Path(a.case_dir or "/mnt/storage3/voxel/MichelStadt/allruns/coarse/case_AMD")
    t = np.fromfile(run / "slices/LDA_2m_times.bin", dtype=">f8")
    mm = np.memmap(run / "slices/LDA_2m.bin", dtype=">f8", mode="r", shape=(t.size, 511, 799, 3))
    i0 = int(np.searchsorted(t, 1000.0))
    idx = np.arange(i0, t.size, a.every)[:a.nframes]
    uref, Lx, Lz = 6.1, 1997.5, 1277.5
    speed = lambda i: np.sqrt((np.asarray(mm[i], dtype=float) ** 2).sum(-1)) / uref
    phi = np.fromfile(run / "SDF_in", dtype=">f8").reshape((801, 63, 511), order="F")[1:-1, 1:-1, :]
    solid = (phi[:, 1, :] < 0).T          # cell-centre SDF of the staircase IBM at the slice height (y = 1.875 m)

    def frame(k):
        fig, ax = plt.subplots(figsize=(6.4, 6.4 * Lz / Lx + 1.2), constrained_layout=True)
        d = np.ma.masked_where(solid, speed(idx[k]))
        cm = plt.get_cmap("magma_r").copy()
        cm.set_bad("#9a9a9a")
        im = ax.imshow(d, origin="lower", extent=[0, Lx, 0, Lz], cmap=cm, vmin=0, vmax=1.3, interpolation="bilinear")
        ax.set_xlabel(r"$x$ [m]"); ax.set_ylabel(r"$z$ [m]")
        cb = fig.colorbar(im, ax=ax, orientation="horizontal", fraction=0.06, pad=0.04)
        cb.set_label(r"$|\mathbf{u}|/U_\mathrm{ref}$")
        fig.suptitle(rf"Michel-Stadt, AMD, $y=1.9$ m:  $t={t[idx[k]]:6.0f}$ s", fontsize=10)
        return fig
    return list(range(len(idx))), frame, 8, 440


def particles(a):
    from dopamine_post.particles import list_snapshots, read_snapshot
    root = Path(a.case_dir or V + "particleValidation")
    sets = []
    for st in (1, 5):
        snaps = list_snapshots(root / f"St{st}" / "fields", f"St{st}")[a.pstart:a.pstart + a.pframes * a.every:a.every]
        sets.append((st, snaps))
    n = min(len(sn) for _, sn in sets)
    re_tau, ywall = 150.0, 15.0 / 150.0
    x, y, z, xm, ym, zm = Snapshot(root / "St1" / "fields" / f"St1.{sets[0][1][0][0]}").coords()
    jy = int(np.argmin(abs(ym - 8.0 / re_tau)))
    Lx, Lz = x[-1], z[-1]
    steps = list(range(n))

    def frame(i):
        fig, axs = plt.subplots(2, 1, figsize=(6.2, 6.2 * Lz / Lx * 2 + 0.7), constrained_layout=True, sharex=True)
        tp = None
        for ax, (st, snaps) in zip(axs, sets):
            step, path = snaps[i]
            u = Snapshot(root / f"St{st}" / "fields" / f"St{st}.{step}").plane("U", 1, jy)
            im = ax.imshow((u - u.mean()).T, origin="lower", extent=[0, Lx, 0, Lz], cmap="RdBu_r", vmin=-3, vmax=3,
                           interpolation="bicubic", alpha=0.55)
            pid, d = read_snapshot(path)
            tp = float(np.median(d[:, 9])) * re_tau
            m = (d[:, 1] < ywall) & (pid % a.thin == 0)
            ax.scatter(d[m, 0], d[m, 2], s=0.7, c="k", alpha=0.7, linewidths=0, rasterized=True)
            ax.set_ylabel(r"$z/h$")
            ax.set_title(rf"$St^+={st}$, particles with $y^+<15$ ({m.sum()} shown)", fontsize=9, loc="left")
            ax.set_xlim(0, Lx); ax.set_ylim(0, Lz)
        axs[-1].set_xlabel(r"$x/h$")
        cb = fig.colorbar(im, ax=axs, orientation="vertical", fraction=0.04, pad=0.02)
        cb.set_label(r"$u'/u_\tau$ at $y^+\approx8$")
        fig.suptitle(rf"Inertial particles in channel flow, $Re_\tau=150$:  $t^+={tp:6.0f}$", fontsize=10)
        return fig
    return steps, frame, 8, 440


P3D_RUN = "/tmp/claude-1000/-home-akshay-CODE-fdm-dopamine-public/a8ddbc89-b15e-4490-8b28-15cf986d3e5a/scratchpad/p3d"


def particles3d(a):
    """Q-criterion vortices and 20 tracked St+=5 particles from a short fine-cadence restart (<case-dir>/fields/P3D.*)."""
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    from scipy.ndimage import gaussian_filter, map_coordinates
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    from skimage.measure import marching_cubes
    from dopamine_post.fields import FieldSnapshot
    from dopamine_post.particles import read_snapshot
    case = Path(a.case_dir or P3D_RUN)
    steps = list_steps(case / "fields", "P3D")[::a.every]
    nu, qplus, ymax, tail, vex = 1.0 / 150.0, a.qplus, 0.6, 14, 4.0
    first = FieldSnapshot.read(steps[0][1], 0, 0, 0, 0)
    xm, ym, zm = first.xm, first.ym, first.zm
    Lx, Lz = first.x[-1], first.z[-1]
    jm = int(np.searchsorted(ym, ymax))
    traj, tplus = [], []
    for s, _ in steps:
        ids, d = read_snapshot(case / "fields" / f"P3D_particles.{s}")
        o = np.argsort(ids)
        traj.append(d[o][:, :3])
        tplus.append(150.0 * float(np.median(d[:, 9])))          # particle age = t - 45, t+ = Re_tau * age
    traj = np.array(traj)                                   # (frames, N, 3)
    ok = traj[:, :, 1].max(0) < 0.3                         # stay within y+ < 45 over the whole window
    dxz = np.abs(np.diff(traj[:, :, [0, 2]], axis=0))
    ok &= ~((dxz[:, :, 0] > Lx / 2) | (dxz[:, :, 1] > Lz / 2)).any(0)   # no periodic wrap, so tracks stay whole
    ym_mean = traj[:, :, 1].mean(0)
    half = a.nparticles // 2
    groups = []
    for lo, hi in ((0.0, 0.1), (0.1, 0.3)):                 # y+ < 15 (buffer layer) and 15 < y+ < 45
        c = np.flatnonzero(ok & (ym_mean >= lo) & (ym_mean < hi))
        c = c[np.abs(traj[-1, c, 0] - traj[0, c, 0]) > np.percentile(np.abs(traj[-1, c, 0] - traj[0, c, 0]), 70)]   # moving ones, so tracks show
        c = c[np.argsort(traj[0, c, 2])]
        groups.append(c[np.linspace(0, len(c) - 1, half).astype(int)])
    tr = [traj[:, g, :] for g in groups]
    gcol = ("#D55E00", "#CC79A7")
    cm, norm = plt.get_cmap(a.cmap), plt.Normalize(a.umin, a.umax)
    el, az = np.radians(28), np.radians(-62)
    light = np.array([np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)])
    light = light + np.array([0.0, 0.0, 0.35])
    light /= np.linalg.norm(light)

    def frame(i):
        f = FieldSnapshot.read(steps[i][1], 0, 0, 0, 0)
        g = [[np.gradient(f[c], xyz, axis=ax) for ax, xyz in enumerate((xm, ym, zm))] for c in ("U", "V", "W")]
        Q = -0.5 * sum(g[p][q] * g[q][p] for p in range(3) for q in range(3))
        q = gaussian_filter(Q[:, :jm, :], a.smooth)
        v, fc, nv, _ = marching_cubes(q, qplus / nu**2, step_size=1)
        n_v = len(v)
        comp = connected_components(coo_matrix((np.ones(fc.size), (np.repeat(fc[:, 0], 3), fc.ravel())), shape=(n_v, n_v)),
                                    directed=False)[1]
        big = np.bincount(comp[fc[:, 0]]) >= a.min_faces
        fc = fc[big[comp[fc[:, 0]]]]
        P = np.c_[np.interp(v[:, 0], np.arange(len(xm)), xm), np.interp(v[:, 2], np.arange(len(zm)), zm),
                  np.interp(v[:, 1], np.arange(jm), ym[:jm])]         # plot axes: x, z, wall distance
        N = np.c_[nv[:, 0], nv[:, 2], nv[:, 1]]
        tri = P[fc]
        nrm = N[fc].mean(1)
        nrm /= np.linalg.norm(nrm, axis=1)[:, None] + 1e-30
        uv = map_coordinates(f["U"][:, :jm, :], v.T, order=1, mode="nearest")      # u/u_tau at the vertices (u_tau = 1)
        col = cm(norm(uv[fc].mean(1)))
        ndl = np.abs(nrm @ light)
        spec = 0.30 * ndl**24
        col[:, :3] = np.clip(col[:, :3] * (0.42 + 0.58 * ndl)[:, None] + spec[:, None], 0, 1)
        col[:, 3] = 0.9
        fig = plt.figure(figsize=(5.2, 3.2))
        ax = fig.add_axes([-0.13, -0.38, 1.26, 1.65], projection="3d", computed_zorder=False)
        ax.add_collection3d(Poly3DCollection([[(0, 0, -0.03), (Lx, 0, -0.03), (Lx, Lz, -0.03), (0, Lz, -0.03)]],
                                             facecolors="#ececec", edgecolors="#bbbbbb", linewidths=0.6, zorder=0))
        ax.add_collection3d(Poly3DCollection(tri, facecolors=col, linewidths=0, zorder=2))
        k0 = max(0, i - tail)
        for t_, c_ in zip(tr, gcol):
            for p in range(t_.shape[1]):
                for k in range(k0, i):
                    w = (k - k0 + 1) / (i - k0 + 1)
                    ax.plot(t_[k:k + 2, p, 0], t_[k:k + 2, p, 2], t_[k:k + 2, p, 1], color=c_, alpha=0.12 + 0.88 * w,
                            lw=0.5 + 2.3 * w, solid_capstyle="round", zorder=9)
            ax.scatter(t_[i, :, 0], t_[i, :, 2], t_[i, :, 1], s=75, c=c_, edgecolors="k", linewidths=0.8,
                       depthshade=False, zorder=10)
        ax.set_xlim(0, Lx); ax.set_ylim(0, Lz); ax.set_zlim(-0.03, ymax)
        ax.set_box_aspect((Lx, Lz, ymax * vex))
        ax.view_init(28, -62)
        ax.set_axis_off()
        cax = fig.add_axes([0.66, 0.82, 0.29, 0.03])
        cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cm), cax=cax, orientation="horizontal")
        cb.set_label(r"$u/u_\tau$", fontsize=8, labelpad=1)
        cb.ax.xaxis.set_label_position("top")
        cb.ax.tick_params(labelsize=7, length=2)
        fig.suptitle(rf"$St^+=5$, $Q^+={qplus:g}$ vortices, 20 particles:  $t^+={tplus[i]:6.0f}$",
                     fontsize=9, y=0.985)
        return fig
    return list(range(len(steps))), frame, 9, 600


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case", choices=["tgv", "chan395", "wavywall", "rb", "particles", "urban", "particles3d"])
    ap.add_argument("--out", required=True)
    ap.add_argument("--case-dir")
    ap.add_argument("--every", type=int, default=2)
    ap.add_argument("--nframes", type=int, default=48)
    ap.add_argument("--fps", type=float)
    ap.add_argument("--colors", type=int, default=64)
    ap.add_argument("--dither", default=None, help="ffmpeg paletteuse dither (default: bayer:bayer_scale=5, none for particles)")
    ap.add_argument("--pframes", type=int, default=20, help="particles: number of snapshots")
    ap.add_argument("--pstart", type=int, default=11, help="particles: first snapshot index (t+ = 300 (index + 1))")
    ap.add_argument("--qplus", type=float, default=0.016, help="particles3d: Q+ = Q nu^2/u_tau^4 level")
    ap.add_argument("--nparticles", type=int, default=20)
    ap.add_argument("--cmap", default="viridis")
    ap.add_argument("--umin", type=float, default=2.0)
    ap.add_argument("--umax", type=float, default=18.0)
    ap.add_argument("--dpi", type=int)
    ap.add_argument("--smooth", type=float, default=0.9)
    ap.add_argument("--min-faces", type=int, default=300)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--thin", type=int, default=3, help="particles: keep ids divisible by this")
    a = ap.parse_args()
    steps, frame, fps, width = {"tgv": tgv, "chan395": chan395, "wavywall": wavywall, "rb": rb, "particles": particles, "urban": urban,
                                   "particles3d": particles3d}[a.case](a)
    render(frame, len(steps), a.out, fps=a.fps or fps, width=width, colors=a.colors, workers=a.workers, dpi=a.dpi or (200 if a.case == "particles3d" else 110),
           dither=a.dither or ("none" if a.case.startswith("particles") else "bayer:bayer_scale=5"))
    print("wrote", a.out, len(steps), "frames")


if __name__ == "__main__":
    main()
