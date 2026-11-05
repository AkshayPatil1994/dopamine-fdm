#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
"""Actuator-disk hover against the vortex-cylinder solution on three grids.

Usage: python3 plot_adm.py --case /mnt/storage1/fdm-dopamine/validation/adm --out docs/assets/validation/adm
The centreline probe of each grid is reduced with the case's compute_adm_profile.py (same physical-time window), in memory;
nothing is written into the case directory.
"""
import argparse
import json
import os
import sys

import numpy as np

import plotstyle

plotstyle.apply()
import matplotlib.pyplot as plt  # noqa: E402

GRIDS = [("32", r"$R/\Delta y=3.2$", "#E69F00", "d"), ("64", r"$R/\Delta y=6.4$", "#009E73", "s"),
         ("96", r"$R/\Delta y=9.6$", "#0072B2", "^")]


def reduce_case(cap, case_dir):
    p = cap.read_case_params(case_dir)
    vh = np.sqrt(p["Tkin"] / (2.0 * np.pi * p["R"] ** 2))
    steps, times = cap.read_step_time_table(os.path.join(case_dir, "run.log"))
    lo, hi = np.interp(cap.T_LO, times, steps), np.interp(cap.T_HI, times, steps)
    if cap.T_HI > times[-1]:
        hi = steps[-1]
    V = cap.load_line(os.path.join(case_dir, "lines", "centerline"))
    ns, _, npts = V.shape
    i_lo = int(np.clip(round(lo / p["line_freq"]), 0, ns - 1))
    i_hi = int(np.clip(round(hi / p["line_freq"]), i_lo + 1, ns))
    ym = (np.arange(npts) + 0.5) * p["Ly"] / npts
    s = (p["yc"] - ym) / p["R"]
    seg = V[i_lo:i_hi, 0, :]
    t_win = (float(np.interp(i_lo * p["line_freq"], steps, times)), float(np.interp(i_hi * p["line_freq"], steps, times)))
    return dict(s=s, v=seg.mean(0), sd=seg.std(0), vh=vh, R=p["R"], dy=p["Ly"] / npts, window=t_win,
                snaps=(i_lo, i_hi), t_end=float(times[-1]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    sys.dont_write_bytecode = True
    sys.path.insert(0, a.case)
    import compute_adm_profile as cap

    ref = np.loadtxt(os.path.join(a.case, "vortex_cylinder_reference.dat"), comments="#")
    fig, ax = plt.subplots(1, 2, figsize=(9.0, 3.8), constrained_layout=True, gridspec_kw=dict(width_ratios=[1.5, 1]))
    for x in ax:
        x.plot(ref[:, 0], ref[:, 1], "k-", lw=1.6, label="vortex-cylinder theory")
    out = {"grids": {}}
    for name, lab, col, mk in GRIDS:
        r = reduce_case(cap, os.path.join(a.case, name))
        y = -r["v"] / r["vh"]
        ri = np.interp(r["s"], ref[:, 0], ref[:, 1])
        e = y - ri
        m3, m05 = np.abs(r["s"]) <= 3.0, np.abs(r["s"]) <= 0.5
        for x in ax:
            x.plot(r["s"], y, color=col, marker=mk, ls="none", markersize=3.6, label=lab)
        out["grids"][name] = {
            "R_over_dy": r["R"] / r["dy"], "rmse_3": float(np.sqrt(np.mean(e[m3] ** 2))),
            "rmse_0p5": float(np.sqrt(np.mean(e[m05] ** 2))), "peak_V_over_vh": float(y[m05].max()),
            "far_wake_V_over_vh": float(y[(r["s"] >= 3) & (r["s"] <= 6)].mean()), "window_s": r["window"], "snapshots": r["snaps"],
            "t_end": r["t_end"]}
    ax[0].set_xlim(-10, 10)
    ax[1].set_xlim(-3, 3)
    ax[1].set_ylim(0, 2.2)
    ax[0].set_xlabel(r"$s/R$")
    ax[1].set_xlabel(r"$s/R$")
    ax[0].set_ylabel(r"$-V/v_h$")
    for x in ax:
        x.axvline(0, color=plotstyle.GREY, lw=0.8, ls=":")
    ax[0].legend(loc="upper left")
    plotstyle.panel_label(ax[0], "(a)")
    plotstyle.panel_label(ax[1], "(b)")
    plotstyle.save(fig, os.path.join(a.out, "adm_validation.png"))
    out["ref"] = {"s": ref[::20, 0].tolist(), "V": ref[::20, 1].tolist()}
    with open(os.path.join(a.out, "adm_validation.json"), "w") as fh:
        json.dump(out, fh, separators=(",", ":"))
    for k, v in out["grids"].items():
        print(f"{k}: R/dy={v['R_over_dy']:.1f} RMSE(|s/R|<=3)={v['rmse_3']:.4f} RMSE(<=0.5)={v['rmse_0p5']:.4f} "
              f"peak={v['peak_V_over_vh']:.3f} far={v['far_wake_V_over_vh']:.3f} window={v['window_s'][0]:.2f}-{v['window_s'][1]:.2f}s "
              f"run_end={v['t_end']:.2f}s")


if __name__ == "__main__":
    main()
