#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
"""Wall-modelled channel LES: Re_tau=395 against MKM DNS, and AMD vs Vreman at Re_tau=2000 and 5200 against Lee & Moser.

Usage: python3 plot_channel_wm.py --chan395 .../validation/chan395 --amdwl .../validation/AMDWallLayer \
           --out docs/assets/validation/channel
"""
import argparse
import json
import os
import sys

import numpy as np

import plotstyle

plotstyle.apply()
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from dopamine_post.fields import FieldSnapshot  # noqa: E402
from dopamine_post.tavg import read_tavg  # noqa: E402

BLUE, ORANGE, GREY = plotstyle.BLUE, plotstyle.ORANGE, plotstyle.GREY


def chan395(a, out):
    nu = 1.0 / 395.0
    z = np.load(os.path.join(a.chan395, "stats", "stats_chan395_80000_120000.npz"))
    m = z["ym"] <= 1.0
    yp = z["ym"][m] / nu
    mk = lambda f: np.loadtxt(os.path.join(a.chan395, f), comments="#")
    mean, rs = mk("MKM_MEANS.dat"), mk("MKM_REYSTRESS.dat")
    dyp, dU = mean[:, 1], mean[:, 2]
    kw = dict(color=BLUE, marker="o", ls="none", markersize=3.6)
    fig, ax = plt.subplots(1, 3, figsize=(11.5, 3.6), constrained_layout=True)
    ax[0].semilogx(dyp[1:], dU[1:], "k-", lw=1.6, label="DNS (MKM)")
    ax[0].semilogx(yp, z["U_mean"][m], label="DOPAMINE", **kw)
    yl = np.array([30.0, 395.0])
    ax[0].semilogx(yl, np.log(yl) / 0.41 + 5.2, color=GREY, ls="--", lw=1.1)
    ax[0].set_ylabel(r"$\langle U\rangle^+$")
    ax[0].legend(loc="upper left")
    for name, col, key, ref in (("u", "#D55E00", "urms", 2), ("v", "#009E73", "vrms", 3), ("w", "#CC79A7", "wrms", 4)):
        ax[1].plot(dyp, np.sqrt(np.maximum(rs[:, ref], 0)), color=col, lw=1.6)
        ax[1].plot(yp, z[key][m], color=col, marker="o", ls="none", markersize=3.6, label=rf"${name}'_\mathrm{{rms}}$")
    ax[1].set_ylabel(r"$\sigma_i^+$")
    ax[2].plot(dyp, -rs[:, 5], "k-", lw=1.6)
    ax[2].plot(yp, -z["uv"][m], **kw)
    ax[2].set_ylabel(r"$-\langle u'v'\rangle^+$")
    ax[1].legend(loc="upper right")
    for x, s in zip(ax, "abc"):
        x.set_xlabel(r"$y^+$")
        x.set_xlim(5, 395) if x is ax[0] else x.set_xlim(0, 395)
        plotstyle.panel_label(x, f"({s})")
    plotstyle.save(fig, os.path.join(out, "channel395_validation.png"))
    ym = np.interp(np.array([50.0, 100.0, 200.0, 395.0]), yp, z["U_mean"][m])
    rd = np.interp(np.array([50.0, 100.0, 200.0, 395.0]), dyp, dU)
    return {"yplus": [50, 100, 200, 395], "U_sim": ym.tolist(), "U_dns": rd.tolist(),
            "Ub_sim": float(np.trapz(z["U_mean"][m], z["ym"][m])),
            "Ub_dns": float(np.trapz(np.interp(z["ym"][m] / (1 / 392.24), dyp, dU), z["ym"][m])),
            "yp": yp[::2].tolist(), "U": z["U_mean"][m][::2].tolist()}


def wm(a, out):
    W = a.amdwl
    cases = [
        (2000, 5.0e-4, "LM_Channel_2000_mean_prof.dat",
         [("AMD", BLUE, "-", "s", W + "/cap/runs/re2000_amd_ctl", "b_tavg.250000", "b.250000")]),
        (5200, 1 / 5200.0, "LM_Channel_5200_mean_prof.dat",
         [("Vreman", ORANGE, "--", "o", W + "/re5200_vre", "b_tavg.100000", "b.100000"),
          ("AMD", BLUE, "-", "s", W + "/cap/runs/re5200_amd_ctl", "b_tavg.250000", "b.250000")]),
    ]
    fig, ax = plt.subplots(3, 2, figsize=(7.2, 8.0), constrained_layout=True)
    res = {}
    for j, (Re, nu, fn, runs) in enumerate(cases):
        d = np.loadtxt(os.path.join(W, "dns", fn), comments="%")
        yd, Ud = d[:, 0], d[:, 2]
        mm = yd <= 1
        Ubd = np.trapz(Ud[mm], yd[mm])
        ax[0, j].plot(yd[mm], Ud[mm], "k-", lw=1.8, label="DNS (Lee & Moser)")
        ax[1, j].axhline(0, color="k", lw=0.8)
        res[Re] = {"Ub_dns": float(Ubd), "runs": {}}
        for lab, c, ls, mk, run, tf, sf in runs:
            t, tv = read_tavg(os.path.join(run, "fields", tf))
            U = tv["U"].mean((0, 2))
            f = FieldSnapshot.read(os.path.join(run, "fields", sf))
            y, nt = f.ym, f.fields["nu_t"].mean((0, 2)) / nu
            kw = dict(color=c, ls=ls, marker=mk, markersize=3.4, markevery=3)
            err = 100 * (U.mean() / Ubd - 1)
            ax[0, j].plot(y, U, label=rf"{lab} ($U_b$ {err:+.1f}%)", **kw)
            ax[1, j].plot(y, U - np.interp(y, yd, Ud), **kw)
            ax[2, j].plot(y, nt, **kw)
            res[Re]["runs"][lab] = {"Ub": float(U.mean()), "Ub_err_pct": float(err), "t_avg": float(t),
                                    "first_cell_yplus": float(y[0] * Re)}
        ax[0, j].set_title(rf"$Re_\tau={Re}$, first cell $y^+={res[Re]['runs'][runs[0][0]]['first_cell_yplus']:.0f}$",
                           fontsize=10)
        for r in range(3):
            ax[r, j].set_xscale("log")
            ax[r, j].set_xlim(0.006, 1.0)
            plotstyle.panel_label(ax[r, j], f"({'abc'[r]}{j + 1})")
        ax[2, j].set_xlabel(r"$y/\delta$")
        ax[0, j].legend(loc="lower right")
        ax[1, j].set_ylim(-1.5, 3.0)
    ax[0, 0].set_ylabel(r"$\langle U\rangle^+$")
    ax[1, 0].set_ylabel(r"$\langle U\rangle^+-U^+_\mathrm{DNS}$")
    ax[2, 0].set_ylabel(r"$\langle\nu_t\rangle/\nu$")
    plotstyle.save(fig, os.path.join(out, "wallmodel_re2000_5200.png"))
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chan395", required=True)
    ap.add_argument("--amdwl", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    out = {"channel395": chan395(a, a.out), "wallmodel": wm(a, a.out)}
    with open(os.path.join(a.out, "channel_validation.json"), "w") as fh:
        json.dump(out, fh, separators=(",", ":"))
    c = out["channel395"]
    print("chan395 Ub sim %.3f dns %.3f (%+.1f%%)" % (c["Ub_sim"], c["Ub_dns"], 100 * (c["Ub_sim"] / c["Ub_dns"] - 1)))
    for Re, r in out["wallmodel"].items():
        print(Re, {k: round(v["Ub_err_pct"], 2) for k, v in r["runs"].items()})


if __name__ == "__main__":
    main()
