#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
"""Ting & Kirby (1994) spilling breaker: surface elevation, undertow and rms velocity against the digitised experiment.

Usage: python3 plot_tk94.py --case /mnt/storage1/fdm-dopamine/validation/voftest/tk94_spilling --run medium \
           --t0 40 --out docs/assets/validation/tk94
The statistics are evaluated by the case's own compare_simulation.py machinery (phase average over whole periods after --t0,
as in the paper); this script only re-plots them in the manuscript style. Nothing is written into the case directory.
"""
import argparse
import json
import os
import sys

import numpy as np

import plotstyle

plotstyle.apply()
import matplotlib.pyplot as plt  # noqa: E402

G = 9.81
STATIONS = (1, 3, 5, 7)      # 0-based indices of the LDA stations shown in (b), (c): x_TK = 5.945, 7.275, 8.495, 9.725 m
GREY_BAND = dict(color="0.6", alpha=0.45, lw=0, zorder=0)
EXP_KW = dict(marker="o", ls="none", color="k", mfc="none", ms=3.2, mew=0.7)


def load_machinery(case, run, t0):
    """Execute compare_simulation.py up to its plotting driver (the file calls main() at import time)."""
    sys.dont_write_bytecode = True
    sys.path.insert(0, case)
    src = open(os.path.join(case, "compare_simulation.py")).read()
    src = src[:src.index("\ndef main():")]
    argv, sys.argv = sys.argv, ["compare_simulation.py", os.path.join(case, run), "--t0", str(t0)]
    ns = {"__file__": os.path.join(case, "compare_simulation.py"), "__name__": "cs"}
    try:
        exec(compile(src, "compare_simulation.py", "exec"), ns)
    finally:
        sys.argv = argv
    return ns


def ci95(std, n):
    """Half-width of the 95 % confidence interval of a mean over n whole wave periods (Student t, n-1 dof)."""
    from scipy.stats import t as student
    return student.ppf(0.975, np.maximum(n - 1, 1)) * std / np.sqrt(n)


def profile_error(zz, v, ex, n):
    s = ex[ex[:, 0] == n + 1]
    ok = (s[:, 4] > zz.min()) & (s[:, 4] < zz.max())
    vi = np.interp(s[ok, 4], zz, v)
    return float(np.sqrt(np.mean((vi - s[ok, 5]) ** 2))), float(np.sqrt(np.mean(s[ok, 5] ** 2)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", required=True)
    ap.add_argument("--run", default="medium")
    ap.add_argument("--t0", type=float, default=40.0)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    ns = load_machinery(a.case, a.run, a.t0)
    run = os.path.join(a.case, a.run)
    d = ns["load_case"](run)
    prof = ns["stations_profiles"](run, d)
    gs = d["gstat"]                                  # x, mean, crest, trough, std..., 
    x, mean, crest, trough = gs[:, 0], gs[:, 1], gs[:, 2], gs[:, 3]
    e4 = np.loadtxt(os.path.join(a.case, "exp", "tk94_fig4_spilling.dat"))
    e5 = np.loadtxt(os.path.join(a.case, "exp", "tk94_fig5_undertow_spilling.dat"))
    e6 = np.loadtxt(os.path.join(a.case, "exp", "tk94_fig6_rms_spilling.dat"))

    fig = plt.figure(figsize=(7.2, 6.6), constrained_layout=True)
    gsp = fig.add_gridspec(3, len(STATIONS), height_ratios=[1.1, 1, 1])
    ax = fig.add_subplot(gsp[0, :])
    nper = float(prof[0]["nper"])
    for c0, sd in ((crest - mean, gs[:, 5]), (trough - mean, gs[:, 6]), (mean, gs[:, 4])):
        ax.fill_between(x, 100 * (c0 - ci95(sd, nper)), 100 * (c0 + ci95(sd, nper)), **GREY_BAND)
    ax.fill_between([], [], [], label="95% CI (DOPAMINE)", **{k: v for k, v in GREY_BAND.items() if k != "zorder"})
    ax.plot(x, 100 * (crest - mean), color=plotstyle.BLUE, lw=1.4, label="DOPAMINE")
    ax.plot(x, 100 * (trough - mean), color=plotstyle.BLUE, lw=1.4)
    ax.plot(x, 100 * mean, color=plotstyle.BLUE, lw=1.4)
    ax.plot(e4[:, 0], e4[:, 1], label="Ting & Kirby (1994)", **EXP_KW)
    ax.plot(e4[:, 0], e4[:, 2], **EXP_KW)
    ax.plot(e4[:, 0], e4[:, 3], **EXP_KW)
    ax.set_xlabel(r"$x_\mathrm{TK}$ [m]")
    ax.set_ylabel(r"surface elevation [cm]")
    ax.legend(loc="upper right")
    for txt, yy in ((r"$\zeta_\mathrm{max}-\bar\zeta$", 9.0), (r"$\bar\zeta$", 0.8), (r"$\zeta_\mathrm{min}-\bar\zeta$", -5.0)):
        ax.text(-1.9, yy, txt, ha="left", va="bottom", fontsize=8)
    ax.set_ylim(-5.6, 16)
    plotstyle.panel_label(ax, "(a)")

    # wave-height statistics against the experiment
    Hs = 100 * (crest - trough)
    xe = e4[:, 0]
    He = e4[:, 2] - e4[:, 3]
    inside = (x >= xe[0]) & (x <= xe[-1])
    Hi = np.interp(x[inside], xe, He)
    metrics = {
        "run": a.run, "t0": a.t0,
        "H_breaker_les_cm": float(Hs.max()), "x_breaker_les": float(x[np.argmax(Hs)]),
        "H_breaker_exp_cm": float(He.max()), "x_breaker_exp": float(xe[np.argmax(He)]),
        "rms_H_cm": float(np.sqrt(np.mean((Hs[inside] - Hi) ** 2))),
        "rms_mean_cm": float(np.sqrt(np.mean((100 * mean[inside] - np.interp(x[inside], xe, e4[:, 1])) ** 2))),
        "ci95_halfwidth_H_cm_median": float(np.median(100 * (ci95(gs[:, 5], nper) + ci95(gs[:, 6], nper)))),
        "n_periods": int(nper), "stations": {}, "x": x.tolist(), "mean_cm": (100 * mean).tolist(), "crest_cm": (100 * (crest - mean)).tolist(),
        "trough_cm": (100 * (trough - mean)).tolist(),
    }
    for row, (key, ex, xl, xlim) in enumerate(((("u", e5, r"$\bar u/\sqrt{gh}$", (-0.2, 0.2))),
                                               (("urms", e6, r"$\sqrt{\overline{u'^2}}/\sqrt{gh}$", (0, 0.12)))), start=1):
        for col, n in enumerate(STATIONS):
            q = fig.add_subplot(gsp[row, col])
            p = prof[n]
            zz = (p["y"] - p["depth"] - p["zeta"]) / p["h"]
            fl = (zz > -1.0) & (zz < 0.0) & (p["y"] > p["ybed"] + (p["y"][1] - p["y"][0]))
            v = p[key] / np.sqrt(G * p["h"])
            s = ex[ex[:, 0] == n + 1]
            sd = (p["u_std"] if key == "u" else p["urms_std"]) / np.sqrt(G * p["h"])
            hw = ci95(sd, p["nper"])
            q.fill_betweenx(zz[fl], (v - hw)[fl], (v + hw)[fl], **GREY_BAND)
            q.plot(v[fl], zz[fl], color=plotstyle.BLUE, lw=1.4)
            q.plot(s[:, 5], s[:, 4], **EXP_KW)
            q.set_ylim(-1, 0)
            q.set_xlim(*xlim)
            if key == "u":
                q.axvline(0, color=plotstyle.GREY, lw=0.7, ls=":")
                q.set_title(rf"$x_\mathrm{{TK}}={p['x']:.2f}$ m", fontsize=8.5, loc="right")
            q.set_xlabel(xl)
            if col == 0:
                q.set_ylabel(r"$(z-\bar\zeta)/h$")
            else:
                q.set_yticklabels([])
            plotstyle.panel_label(q, f"({'bc'[row - 1]}{col + 1})")
            err, rms_e = profile_error(zz[fl], v[fl], ex, n)
            metrics["stations"].setdefault(str(n + 1), {"x_TK": float(p["x"])})[f"{key}_rms_err"] = err
            metrics["stations"][str(n + 1)][f"{key}_exp_rms"] = rms_e
    plotstyle.save(fig, os.path.join(a.out, "tk94_validation.png"))
    with open(os.path.join(a.out, "tk94_validation.json"), "w") as fh:
        json.dump(metrics, fh, separators=(",", ":"))
    print("breaker: LES H=%.2f cm at x=%.2f; exp H=%.2f cm at x=%.2f; rms H %.2f cm, rms mean %.2f cm"
          % (metrics["H_breaker_les_cm"], metrics["x_breaker_les"], metrics["H_breaker_exp_cm"], metrics["x_breaker_exp"],
             metrics["rms_H_cm"], metrics["rms_mean_cm"]))
    for k, v in metrics["stations"].items():
        print(k, {kk: round(vv, 4) for kk, vv in v.items()})


if __name__ == "__main__":
    main()
