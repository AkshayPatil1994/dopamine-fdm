#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
"""Michel-Stadt urban flow: mean and rms streamwise velocity on the tunnel centre plane (Y = 0) against the wind tunnel
(Leitl and Harms) and the SOD2D spectral-element LES of Teng et al. (2025, Case 2), for two grids and two SGS models.

Usage: python3 plot_michelstadt.py --root /mnt/storage3/voxel/MichelStadt --out docs/assets/validation/michelstadt
Runs (relative to --root/allruns): coarse/case (Vreman, dx=2.5 m), coarse/case_AMD (AMD, 2.5 m), fine/AMD (AMD, 1.25 m).
Reference data: post/paper_digitised (digitised from the paper's Figs. 7 to 9 and the tunnel points as plotted there).
"""
import argparse
import csv
import json
import os
import re
import sys
from pathlib import Path

import numpy as np

import plotstyle

plotstyle.apply()
import matplotlib.pyplot as plt  # noqa: E402

UREF, H = 6.1, 24.0
Z95 = 1.959964
EXCLUDE = ("b",)                       # station inside the first solid cell of a wall in the runs
RUNS = [  # label, directory under allruns, window [s], colour, linestyle
    ("Vreman, $\\Delta x=2.5$ m", "coarse/case", (1000.0, -1), plotstyle.ORANGE, "--"),
    ("AMD, $\\Delta x=2.5$ m", "coarse/case_AMD", (1000.0, -1), plotstyle.BLUE, "-"),
    ("AMD, $\\Delta x=1.25$ m", "fine/AMD", (300.0, -1), "#009E73", "-"),
]


def neff(x):
    n = x.shape[0]
    d = (x - x.mean(axis=0)).reshape(n, -1)
    f = np.fft.rfft(d, 2 * n, axis=0)
    ac = np.fft.irfft(f * np.conj(f), axis=0)[:n]
    var = ac[0]
    ok = var > 0
    rho = np.zeros_like(ac)
    rho[:, ok] = ac[:, ok] / var[ok]
    neg = rho[1:] <= 0
    first = np.where(neg.any(axis=0), neg.argmax(axis=0) + 1, n - 1)
    k = np.arange(1, n)[:, None]
    s = np.where(k < first[None, :], rho[1:], 0.0).sum(axis=0)
    return np.where(ok, n / np.maximum(1 + 2 * s, 1.0), 1.0).reshape(x.shape[1:])


def mean_ci(x):
    return x.mean(axis=0), Z95 * x.std(axis=0, ddof=1) / np.sqrt(neff(x))


def rms_ci(x):
    v = (x - x.mean(axis=0)) ** 2
    sig = np.sqrt(v.mean(axis=0))
    return sig, Z95 * v.std(axis=0, ddof=1) / np.sqrt(neff(v)) / np.maximum(2 * sig, 1e-300)


def profiles(run, post, window, stations):
    """mean u, v and u_rms (+95 % half widths) at the stations from the run's centre-plane slice XZ_y0"""
    blk = (run / "input_parameters").read_text().split("&DOMAIN")[1].split("/")[0]
    g = dict(re.findall(r"(\w+)\s*=\s*([0-9.eE+-]+)", blk))
    nx, Lx = int(g["nx"]), float(g["Lx"])
    xc = (np.arange(nx - 1) + 0.5) * Lx / (nx - 1)
    icell = np.array([int(np.argmin(abs(xc - float(s["x_solver_m"])))) for s in stations])
    gridf = run / "sdf" / "data" / "grid.out"
    if not gridf.is_file():
        gridf = post.parent / "coarse_runs" / "case_coarse_2p5m" / "sdf" / "data" / "grid.out"
    zH = np.loadtxt(gridf)[:, 1] / H
    meta = {k.strip(): v.strip() for k, v in (ln.split("=", 1) for ln in (run / "slices/XZ_y0_meta.txt").read_text().splitlines()
                                              if "=" in ln)}
    nc, n1, n2, comps = int(meta["ncomp"]), int(meta["n1"]), int(meta["n2"]), meta["comps"]
    t = np.fromfile(run / "slices/XZ_y0_times.bin", dtype=">f8")
    ns = min((run / "slices/XZ_y0.bin").stat().st_size // (8 * nc * n1 * n2), t.size)
    mm = np.memmap(run / "slices/XZ_y0.bin", dtype=">f8", mode="r", shape=(ns, n2, n1, nc))
    keep, tn = [], np.inf
    for i in range(ns - 1, -1, -1):
        if t[i] < tn:
            keep.append(i)
            tn = t[i]
    keep = np.array(sorted(keep))
    lo, hi = window
    sel = keep[(t[keep] >= lo) & ((t[keep] <= hi) if hi > 0 else (np.arange(len(keep)) < len(keep) - 1))]
    f = np.stack([np.asarray(mm[i][:, icell, :]) for i in sel])
    iu, iw = comps.index("U"), comps.index("W")
    U, V = f[..., iu] / UREF, -f[..., iw] / UREF
    out = {}
    for q, x in (("u", U), ("v", V)):
        m, h = mean_ci(x)
        out[q] = (m.T, h.T)
    s, h = rms_ci(U)
    out["urms"] = (s.T, h.T)
    return out, zH, (float(t[sel[0]]), float(t[sel[-1]]), len(sel))


def rms_diff(zH, P, zr, vr, lo, hi):
    """rms of (P - reference) on the run's z/H levels inside (lo, hi], over levels covered by the reference curve"""
    m = (zH > lo) & (zH <= hi) & (zH >= zr.min()) & (zH <= zr.max())
    return float(np.sqrt(np.mean((P[m] - np.interp(zH[m], zr, vr)) ** 2))), int(m.sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="/mnt/storage3/voxel/MichelStadt")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    root = Path(a.root)
    post = root / "post"
    os.makedirs(a.out, exist_ok=True)
    stations = list(csv.DictReader(open(post / "paper_digitised/stations.csv")))
    exp = list(csv.DictReader(open(post / "paper_digitised/exp_dots.csv")))
    cur = list(csv.DictReader(open(post / "paper_digitised/paper_fine_curve.csv")))
    res = {}
    for lab, d, win, col, ls in RUNS:
        res[lab] = profiles(root / "allruns" / d, post, win, stations)
        print(lab, "window %.0f..%.0f s, %d samples" % res[lab][2])
    shown = [(k, s["station"]) for k, s in enumerate(stations) if s["station"] not in EXCLUDE]

    def ref_curve(q, s):
        c = sorted((r for r in cur if r["quantity"] == q and r["station"] == s), key=lambda r: float(r["z_over_H"]))
        return np.array([float(r["z_over_H"]) for r in c]), np.array([float(r["value"]) for r in c])

    XL = {"u": (-0.42, 1.18), "v": (-0.21, 0.21), "urms": (0.04, 0.30)}
    LAB = {"u": r"$\overline{u}/U_\mathrm{ref}$", "v": r"$\overline{v}/U_\mathrm{ref}$", "urms": r"$u'_\mathrm{rms}/U_\mathrm{ref}$"}
    names = {"u": "mean", "urms": "rms"}
    for fig_q in ("u", "urms"):
        fig, ax = plt.subplots(2, 4, figsize=(7.2, 4.3), sharex=True, sharey=True, constrained_layout=True)
        for p, (k, s) in enumerate(shown):
            A = ax.flat[p]
            cz, cv = ref_curve(fig_q, s)
            gap = np.where(np.diff(cz) > 0.06)[0] + 1
            A.plot(np.insert(cv, gap, np.nan), np.insert(cz, gap, np.nan), color="k", lw=1.0, ls="-")
            for (lab, d, win, col, ls) in RUNS:
                P, zH, _ = res[lab]
                val, half = P[fig_q][0][k], P[fig_q][1][k]
                A.fill_betweenx(zH, val - half, val + half, color=col, alpha=0.25, lw=0)
                A.plot(val, zH, color=col, ls=ls, lw=1.0)
            e = [r for r in exp if r["quantity"] == fig_q and r["station"] == s]
            A.plot([float(r["value"]) for r in e], [float(r["z_over_H"]) for r in e], "o", mfc="none", mec="0.35", mew=0.6, ms=2.8)
            xh = float(stations[k]["x_over_H"])
            A.set_title(rf"({chr(97 + p)}) $x/H={0.0 if abs(xh) < 0.05 else xh:.1f}$", loc="left", fontsize=7.5, pad=2)
            A.axhline(2.0, color=plotstyle.GREY, lw=0.6, ls=":")
            A.set_ylim(0, 5)
            A.set_xlim(*XL[fig_q])
            if p >= 3:
                A.set_xlabel(LAB[fig_q])
            if p % 4 == 0:
                A.set_ylabel(r"$z/H$")
        leg = ax.flat[7]
        leg.axis("off")
        from matplotlib.lines import Line2D
        from matplotlib.patches import Patch
        h = [Line2D([], [], color="k", lw=1.0), Line2D([], [], marker="o", ls="", mfc="none", mec="0.35", mew=0.6, ms=2.8)]
        n = ["SOD2D (Teng et al.)", "Wind tunnel"]
        for (lab, d, win, col, ls) in RUNS:
            h.append((Patch(facecolor=col, alpha=0.25, lw=0), Line2D([], [], color=col, ls=ls, lw=1.0)))
            n.append(lab)
        from matplotlib.legend_handler import HandlerTuple
        leg.legend(h, n, loc="center", fontsize=6.8, handler_map={tuple: HandlerTuple(ndivide=None, pad=0)}, handlelength=2.2)
        fig.savefig(os.path.join(a.out, f"michelstadt_{names[fig_q]}_profiles.png"), dpi=240)
        plt.close(fig)

    metrics = {"windows": {lab: {"t0": res[lab][2][0], "t1": res[lab][2][1], "n": res[lab][2][2]} for lab, *_ in RUNS}, "rms_diff": {}}
    bands = {"z/H<=2": (0.0, 2.0), "z/H>2": (2.0, 5.0)}
    for q in ("u", "v", "urms"):
        for bname, (lo, hi) in bands.items():
            for (lab, *_r) in RUNS:
                P, zH, _ = res[lab]
                vals = []
                for k, s in shown:
                    zr, vr = ref_curve(q, s.__str__())
                    v, n = rms_diff(zH, P[q][0][k], zr, vr, lo, hi)
                    if n:
                        vals.append((v, n))
                tot = sum(n for _, n in vals)
                metrics["rms_diff"].setdefault(q, {}).setdefault(bname, {})[lab] = float(np.sqrt(sum(v * v * n for v, n in vals) / tot))
    # inter-resolution difference of the AMD runs (same SGS model)
    coarse, fine = res[RUNS[1][0]], res[RUNS[2][0]]
    metrics["amd_coarse_vs_fine_rms"] = {}
    for q in ("u", "v", "urms"):
        for bname, (lo, hi) in bands.items():
            acc = []
            for k, s in shown:
                m = (coarse[1] > lo) & (coarse[1] <= hi)
                acc.append(np.mean((coarse[0][q][0][k][m] - np.interp(coarse[1][m], fine[1], fine[0][q][0][k])) ** 2))
            metrics["amd_coarse_vs_fine_rms"].setdefault(q, {})[bname] = float(np.sqrt(np.mean(acc)))
    # approach-flow station (a, x/H = -18.3, upstream of the city): signed mean-u difference from SOD2D above z/H = 2
    metrics["approach_station_a_u_signed_diff_zH_gt2"] = {}
    ka = [k for k, st in shown if st == "a"][0]
    zr, vr = ref_curve("u", "a")
    for (lab, *_r) in RUNS:
        P, zH, _ = res[lab]
        m = (zH > 2.0) & (zH <= 4.5)
        metrics["approach_station_a_u_signed_diff_zH_gt2"][lab] = float(np.mean(P["u"][0][ka][m] - np.interp(zH[m], zr, vr)))
    metrics["signed_mean_u_diff_zH_gt2_all_stations"] = {}
    for (lab, *_r) in RUNS:
        P, zH, _ = res[lab]
        m = (zH > 2.0) & (zH <= 4.5)
        metrics["signed_mean_u_diff_zH_gt2_all_stations"][lab] = float(np.mean(
            [np.mean(P["u"][0][k][m] - np.interp(zH[m], *ref_curve("u", s))) for k, s in shown]))
    # convergence figure
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.5), constrained_layout=True, sharey=False)
    labs = [r[0] for r in RUNS]
    xs = np.arange(3)
    for A, q in zip(ax, ("u", "v", "urms")):
        for j, (bn, hatch) in enumerate((("z/H<=2", None), ("z/H>2", "//"))):
            vv = [metrics["rms_diff"][q][bn][l] for l in labs]
            A.bar(xs + (j - 0.5) * 0.38, vv, 0.36, color=[r[3] for r in RUNS], hatch=hatch, edgecolor="k", lw=0.5)
        A.set_xticks(xs)
        A.set_xticklabels(["Vreman\n2.5 m", "AMD\n2.5 m", "AMD\n1.25 m"], fontsize=7)
        A.set_ylabel(r"rms difference from SOD2D / $U_\mathrm{ref}$" if q == "u" else "")
        A.set_title(LAB[q], fontsize=8.5)
    from matplotlib.patches import Patch
    ax[2].legend([Patch(fc="0.8", ec="k", lw=0.5), Patch(fc="0.8", ec="k", lw=0.5, hatch="//")], [r"$z/H\leq2$", r"$z/H>2$"], fontsize=7)
    plotstyle.panel_label(ax[0], "(a)")
    plotstyle.panel_label(ax[1], "(b)")
    plotstyle.panel_label(ax[2], "(c)")
    plotstyle.save(fig, os.path.join(a.out, "michelstadt_convergence.png"))
    with open(os.path.join(a.out, "michelstadt_validation.json"), "w") as fh:
        json.dump(metrics, fh, separators=(",", ":"))
    print(json.dumps(metrics, indent=1))


if __name__ == "__main__":
    main()
