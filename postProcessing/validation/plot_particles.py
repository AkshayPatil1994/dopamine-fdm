#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
"""One-way-coupled heavy particles in channel flow (Re_tau = 150) against the Marchioli et al. (2008) benchmark.

Usage: python3 plot_particles.py --case /mnt/storage1/fdm-dopamine/validation/particleValidation \
           --out docs/assets/validation/particles
Reuses the slab binning and scoring of the case's own compare_particles.py / valcommon.py (imported, nothing is written there).
Grey bands: min-max envelope of the benchmark groups (UUD, TUE, TUD, HPU, ASU) digitised from the paper.
"""
import argparse
import json
import warnings
import os
import sys
from pathlib import Path

import numpy as np

import plotstyle

plotstyle.apply()
warnings.filterwarnings("ignore", category=RuntimeWarning)
import matplotlib.pyplot as plt  # noqa: E402

BLUE, ORANGE, GREY = plotstyle.BLUE, plotstyle.ORANGE, plotstyle.GREY
BAND = dict(color="0.65", alpha=0.5, lw=0)
WINDOWS = {1: (18000.0, 20250.0), 5: (12000.0, 14100.0)}     # t+ (time since seeding) pooled per Stokes number


def band(ax, curves, x, log_x=False):
    env = vc.envelope(curves, x, log_x)
    if env is not None:
        ax.fill_between(x, env[0], env[1], **BAND)


def main():
    global vc, cp
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-fluid", type=int, default=6, help="St1 field snapshots averaged for the fluid check")
    a = ap.parse_args()
    sys.dont_write_bytecode = True
    sys.path.insert(0, a.case)
    import valcommon as vc          # noqa: E402
    import compare_particles as cp  # noqa: E402
    os.makedirs(a.out, exist_ok=True)
    re_tau = vc.RE_TAU
    edges = vc.slab_edges(re_tau)
    zc, dz = 0.5 * (edges[1:] + edges[:-1]), np.diff(edges)
    out = {"fluid": {}, "particles": {}, "note_runs": "St1/St5 run 2026-10-03..05 with the pre-fix executable (x = z = 0 seeding); statistics pooled at t+ >= 12000, particles uniform in x and z by t+ ~ 3300"}

    # ---- fluid statistics seen by the particles (St1 run, last snapshots), both halves folded
    from dopamine_post.fields import FieldSnapshot
    fdir = Path(a.case) / "St1" / "fields"
    snaps = sorted((int(p.name.split(".")[1]), p) for p in fdir.iterdir()
                   if p.name.startswith("St1.") and p.name.split(".")[1].isdigit())[-a.n_fluid:]
    nu = 1.0 / re_tau
    acc = None
    for _, p in snaps:
        s = FieldSnapshot.read(p, sgs_model=0, sediment_flag=0, boussinesq_flag=0)
        U, V, W = (s[k][:-1, :, :-1] for k in ("U", "V", "W"))
        m = lambda f: f.mean(axis=(0, 2))
        cur = [m(U), m(V), m(W), m(U * U), m(V * V), m(W * W)]
        acc = cur if acc is None else [x + c for x, c in zip(acc, cur)]
    Um, Vm, Wm, UU, VV, WW = (x / len(snaps) for x in acc)
    n = len(s.ym)
    fold = lambda f: 0.5 * (f[:n // 2] + f[::-1][:n // 2])
    zp = s.ym[:n // 2] / nu
    Uf = fold(Um)
    rf = {"u": np.sqrt(fold(UU - Um ** 2)), "v": np.sqrt(fold(VV - Vm ** 2)), "w": np.sqrt(fold(WW - Wm ** 2))}
    mkm_m, mkm_r = np.loadtxt(vc.REF / "MKM_chan180_means.dat"), np.loadtxt(vc.REF / "MKM_chan180_reystress.dat")

    fig = plt.figure(figsize=(10.5, 7.6), constrained_layout=True)
    gs = fig.add_gridspec(3, 4)
    ax = fig.add_subplot(gs[0, 0:2])
    cv = vc.paper_curves("final-mean-streamwise-fluid-velocity-complete.csv")
    xx = np.logspace(0, np.log10(150), 200)
    band(ax, cv, xx, True)
    ax.semilogx(mkm_m[1:, 1], mkm_m[1:, 2], "k--", lw=1.1, label=r"MKM $Re_\tau=180$")
    ax.semilogx(zp, Uf, color=BLUE, lw=1.6, label="DOPAMINE")
    ax.fill_between([], [], [], label="benchmark groups", **BAND)
    ax.set(xlim=(0.5, 150), ylim=(0, 20), xlabel=r"$z^+$", ylabel=r"$\langle U\rangle^+$")
    ax.legend(loc="upper left")
    plotstyle.panel_label(ax, "(a)")
    out["fluid"]["U"] = vc.score(zp, Uf, cv, 1.0, 150.0, 0.05, True, label="fluid U+")
    ax = fig.add_subplot(gs[0, 2:4])
    for k, col, fn, mi, lab in (("u", "#D55E00", "final-streamwise-fluid-velocity-fluctuations.csv", 2, "u"),
                                ("v", "#009E73", "final-wallnormal-fluid-velocity-fluctuations.csv", 3, "v"),
                                ("w", "#CC79A7", "final-spanwise-fluid-velocity-fluctuations.csv", 4, "w")):
        cvk = vc.paper_curves(fn)
        xl = np.linspace(0.5, 150, 300)
        band(ax, cvk, xl)
        ax.plot(mkm_r[:, 1], np.sqrt(mkm_r[:, mi]), color=col, ls="--", lw=1.0)
        ax.plot(zp, rf[k], color=col, lw=1.6, label=rf"${lab}'_\mathrm{{rms}}$")
        out["fluid"][k + "rms"] = vc.score(zp, rf[k], cvk, 1.0, 150.0, 0.10, False, label=f"fluid {k} rms")
    ax.set(xlim=(0, 150), xlabel=r"$z^+$", ylabel=r"$\sigma_i^+$")
    ax.legend(loc="upper right")
    plotstyle.panel_label(ax, "(b)")

    # ---- particle statistics
    for row, st in enumerate((1, 5), start=1):
        case = Path(a.case) / f"St{st}"
        sn, Ly = cp.load_case(case, f"St{st}", re_tau)
        tp = np.array([t for _, t, _ in sn])
        lo, hi = WINDOWS[st]
        sel = [i for i, t in enumerate(tp) if lo <= t <= hi]
        ntot = len(sn[0][2])
        cnt = np.zeros(len(zc))
        S = {k: np.zeros(len(zc)) for k in ("u", "uu", "w", "ww", "v", "vv", "uv")}
        for i in sel:
            z_, u, w, vn = cp.fold(sn[i][2], Ly)
            z_ = z_ * re_tau          # fold() returns the wall distance in units of h (u_tau = h = 1)
            c, sm = cp.bin_moments(z_, u, w, vn, edges)
            cnt += c
            for k in S:
                S[k] += sm[k]
        with np.errstate(all="ignore"):
            mm = {k: S[k] / cnt for k in ("u", "w", "v")}
            rms = {k: np.sqrt(np.maximum(S[k + k] / cnt - mm[k] ** 2, 0)) for k in ("u", "w", "v")}
        conc = cnt / len(sel) / ntot * (re_tau / dz)
        ok = cnt >= 30
        panels = [
            ("conc", f"conc-t=20000-St={st}.csv", conc, 0.30, True, r"$C/C_0$", True),
            ("mean_u", f"final-mean-streamwise-particle-velocity-st={st}-complete.csv", mm["u"], 0.05, True,
             r"$\langle V_x\rangle^+$", False),
            ("rms_u", f"final-particle-streamwise-velocity-fluctuations-st={st}.csv", rms["u"], 0.10, False,
             r"$V_{x,\mathrm{rms}}^+$", False),
            ("rms_v", f"final-particle-wallnormal-velocity-fluctuations-st={st}.csv", rms["v"], 0.10, False,
             r"$V_{n,\mathrm{rms}}^+$", False)]
        res = {"t_window": [lo, hi], "n_snapshots": len(sel), "t_plus_last": float(tp[-1]), "n_particles": int(ntot)}
        for col, (key, fname, y, tol, logx, ylab, loglog) in enumerate(panels):
            ax = fig.add_subplot(gs[row, col])
            curves = vc.paper_curves(fname, kind="symbol" if key == "conc" else "line")
            lines = vc.paper_curves(fname)
            xs = np.logspace(-0.3, np.log10(150), 250) if logx else np.linspace(0.5, 150, 300)
            band(ax, lines if lines else curves, xs, logx)
            ax.plot(zc[ok], y[ok], color=BLUE if st == 1 else ORANGE, lw=1.6, marker="o", ms=2.6, mfc="white", mew=0.8)
            ax.set(xscale="log" if logx else "linear", xlim=(1, 150) if key == "conc" else ((0.5, 150) if logx else (0, 150)),
                   xlabel=r"$z^+$", ylabel=ylab)
            if loglog:
                ax.set_yscale("log")
            if col == 0:
                ax.text(0.04, 0.92, rf"$St^+={st}$", transform=ax.transAxes, va="top", fontsize=10)
            plotstyle.panel_label(ax, f"({'cd'[row - 1]}{col + 1})")
            res[key] = vc.score(zc[ok], y[ok], curves, 1.0, 150.0, tol, logx, label=f"St{st} {key}")
        out["particles"][f"St{st}"] = res
    fig.savefig(os.path.join(a.out, "particles_validation.png"), dpi=250)
    out["particles_profiles_zc"] = zc[::3].tolist()
    with open(os.path.join(a.out, "particles_validation.json"), "w") as fh:
        json.dump(out, fh, separators=(",", ":"), default=float)
    for g, d in (("fluid", out["fluid"]), *((k, v) for k, v in out["particles"].items())):
        for k, v in d.items():
            if isinstance(v, dict) and "inside" in v:
                print(f"{g:6s} {k:8s} n={v['n']:3d} inside={100 * v['inside']:5.1f}%  rms_rel={100 * v['rms_rel']:5.1f}%")


if __name__ == "__main__":
    main()
