#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
"""Taylor-Green vortex Re=1600: kinetic energy and dissipation against the Incompact3d 512^3 DNS.

Usage: python3 plot_tgv.py --case /mnt/storage1/fdm-dopamine/validation/tgv --out docs/assets/validation/tgv
Reads <case>/{256,384,512}/tgv_stats_computed.dat and <case>/TGV_Re1600.dat.
"""
import argparse
import json
import os

import numpy as np

import plotstyle

plotstyle.apply()
import matplotlib.pyplot as plt  # noqa: E402

GRIDS = [(256, "#E69F00", "d"), (384, "#009E73", "s"), (512, "#0072B2", "^")]
TMAX = 20.0
EPS_SCALE = 1.0e3


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--sample", type=int, default=2, help="plot every n-th output as a marker")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    ref = np.loadtxt(os.path.join(a.case, "TGV_Re1600.dat"), comments="#")
    rt, rE, reps = ref[:, 0], ref[:, 1], ref[:, 2]
    fig, ax = plt.subplots(1, 3, figsize=(11.5, 3.5), constrained_layout=True)
    ax[0].plot(rt, rE, "k-", lw=1.6, label="DNS (Incompact3d, $512^3$)")
    eps = reps
    ax[1].plot(rt, EPS_SCALE * eps, "k-", lw=1.6)
    out = {"t_dns": rt[::5].tolist(), "eps_dns": (EPS_SCALE * eps[::5]).tolist(), "grids": {}}
    for n, col, mk in GRIDS:
        f = os.path.join(a.case, str(n), "tgv_stats_computed.dat")
        d = np.loadtxt(f, comments="#")
        t, E, e = d[:, 0], d[:, 1], d[:, 2]
        m = t <= TMAX
        t, E, e = t[m], E[m], e[m]
        kw = dict(color=col, marker=mk, ls="none", markersize=4.2, label=rf"DOPAMINE ${n}^3$")
        ax[0].plot(t[::a.sample], E[::a.sample], **kw)
        ax[1].plot(t[::a.sample], EPS_SCALE * e[::a.sample], **kw)
        eref = np.interp(t, rt, eps)
        Eref = np.interp(t, rt, rE)
        mm = t >= 0.5
        ax[2].plot(t[mm], 100 * (e[mm] - eref[mm]) / eps.max(), color=col, lw=1.3, label=rf"${n}^3$")
        ip = np.argmax(e)
        out["grids"][str(n)] = {
            "eps_peak": float(e[ip]), "t_peak": float(t[ip]),
            "eps_peak_err_pct": float(100 * (e[ip] / eps.max() - 1)),
            "t_peak_dns": float(rt[np.argmax(eps)]), "eps_peak_dns": float(eps.max()),
            "L2_eps_pct": float(100 * np.sqrt(np.mean((e[mm] - eref[mm]) ** 2)) / eps.max()),
            "L2_Ek_pct": float(100 * np.sqrt(np.mean((E - Eref) ** 2)) / rE[0]),
            "t_end": float(t[-1]),
        }
    ax[0].set_xlabel(r"$t\,U_0/L$"); ax[0].set_ylabel(r"$E_k$")
    ax[1].set_xlabel(r"$t\,U_0/L$"); ax[1].set_ylabel(r"$-\mathrm{d}E_k/\mathrm{d}t\times10^{3}$")
    ax[2].set_xlabel(r"$t\,U_0/L$"); ax[2].set_ylabel(r"$(\epsilon-\epsilon_\mathrm{DNS})/\epsilon_\mathrm{DNS,max}$ [%]")
    ax[2].axhline(0, color=plotstyle.GREY, lw=0.8, ls=":")
    ax[0].legend(loc="lower left"); ax[2].legend(loc="lower right")
    for x in ax:
        x.set_xlim(0, TMAX)
    ax[1].set_ylim(0, 14)
    for x, s in zip(ax, "abc"):
        plotstyle.panel_label(x, f"({s})")
    plotstyle.save(fig, os.path.join(a.out, "tgv_validation.png"))
    with open(os.path.join(a.out, "tgv_validation.json"), "w") as fh:
        json.dump(out, fh, separators=(",", ":"))
    for n, v in out["grids"].items():
        print(f"{n}^3: peak {v['eps_peak']:.5f} at t={v['t_peak']:.2f} ({v['eps_peak_err_pct']:+.1f}%), "
              f"L2(eps)={v['L2_eps_pct']:.2f}%, L2(Ek)={v['L2_Ek_pct']:.3f}%")


if __name__ == "__main__":
    main()
