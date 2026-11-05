#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
"""Wavy-wall channel (Re_H=6760, 2a/lambda=0.10) against the Maass & Schumann DNS (ERCOFTAC Case 077).

Usage: python3 plot_wavywall.py --case /mnt/storage1/fdm-dopamine/validation/waveWall --out docs/assets/validation/wavywall
Station averaging (span, snapshots, the 4 periodic wavelengths) is reused from <case>/analyse_wavywall.py.
"""
import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np

import plotstyle

STATIONS = (0.195, 0.304, 0.492, 0.804)
# (label, DNS key, sim column in rey_profiles, sign, colour)
STRESSES = [(r"$\overline{u'u'}$", "uu", 0, 1.0, "#D55E00"),
            (r"$\overline{w'w'}$", "ww", 1, 1.0, "#0072B2"),
            (r"$\overline{v'v'}$", "vv", 2, 1.0, "#009E73"),
            (r"$-\overline{u'w'}$", "neg_uw", 3, -1.0, "#CC79A7")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--t-start", type=float, default=75.0, help="discard snapshots before this time (flow-through units)")
    ap.add_argument("--dpi", type=int, default=220)
    ap.add_argument("--stride", type=int, default=4, help="use every n-th snapshot")
    a = ap.parse_args()
    case = Path(a.case)
    sys.path.insert(0, str(case))
    import analyse_wavywall as aw  # noqa: E402  (sets usetex; plotstyle.apply() below resets it)
    import snapshot_io  # noqa: E402

    plotstyle.apply()
    import matplotlib.pyplot as plt  # noqa: E402

    params = snapshot_io.parse_input_parameters(case / "input_parameters")
    tsave = float(params.get("tsave", 0.5))
    sgs, sed = int(params.get("sgs_model", 0)), int(params.get("sediment_flag", 0))
    Lx, ub = float(params.get("lx", 4.0)), float(params.get("ub_target", 1.0))
    wavelength, amp, n_waves = aw.load_wave_geometry(case, argparse.Namespace(wavelength=None, amplitude=None,
                                                                              n_waves=None), Lx)
    fdir = case / "fields"
    snaps = snapshot_io.list_snapshots(str(fdir), "channel")
    steps = [s for s, _ in snaps][int(a.t_start / tsave)::a.stride]
    s0 = snapshot_io.read_snapshot(fdir / f"channel.{steps[0]}", sgs_model=sgs, sediment_flag=sed)
    xm, ym = s0["xm"], s0["ym"]
    stations = [(f, p) for f, p in ((x, case / "mass_schuman" / aw._local_name(x)) for x in STATIONS)]
    dns = [aw.parse_case077_profile(p) for _, p in stations]
    mean, rey, _, counts = aw.compute_station_profiles(str(fdir), "channel", steps, xm, list(STATIONS), wavelength,
                                                       n_waves, sgs, sed)

    fig, ax = plt.subplots(3, 4, figsize=(11.5, 8.2), constrained_layout=True)
    out = {"stations": list(STATIONS), "n_snapshots": len(steps), "avg_time": len(steps) * a.stride * tsave,
           "metrics": {}}
    ymax_stress = max(d["uu"].max() for d in dns)
    for j, (frac, m, r, d) in enumerate(zip(STATIONS, mean, rey, dns)):
        wall = amp * (1 + np.cos(2 * np.pi * frac))
        fluid = ym > wall
        zd = d["z_over_H"] + amp
        zs = d["z_over_H_stress"] + amp
        met = {}
        rows = [(m[:, 0] / ub, d["u"], zd, "U"), (m[:, 1] / ub, d["w"], zd, "V")]
        for i, (sv, dv, z, name) in enumerate(rows):
            ax[i, j].plot(dv, z, "o", color="k", markersize=2.8, markerfacecolor="none", markeredgewidth=0.7, zorder=2,
                          label="DNS (Maaß & Schumann)" if (i, j) == (0, 0) else None)
            ax[i, j].plot(sv[fluid], ym[fluid], "-", color=plotstyle.BLUE, lw=1.4, zorder=3,
                          label="DOPAMINE" if (i, j) == (0, 0) else None)
            ok = (z >= ym[fluid][0]) & (z <= ym[-1])
            e = np.interp(z[ok], ym[fluid], sv[fluid]) - dv[ok]
            met[f"L2_{name}"] = float(np.sqrt(np.mean(e ** 2)))
        for lab, key, col, sign, c in STRESSES:
            sv = sign * r[:, col] / ub ** 2
            ax[2, j].plot(d[key], zs, "o", color=c, markersize=2.8, markerfacecolor="none", markeredgewidth=0.7, zorder=2)
            ax[2, j].plot(sv[fluid], ym[fluid], "-", color=c, lw=1.4, zorder=3, label=lab if j == 0 else None)
            ok = (zs >= ym[fluid][0]) & (zs <= ym[-1])
            e = np.interp(zs[ok], ym[fluid], sv[fluid]) - d[key][ok]
            met[f"L2_{key}"] = float(np.sqrt(np.mean(e ** 2)))
            met[f"L2_{key}_rel"] = float(100 * met[f"L2_{key}"] / d[key].max())
        met["L2_U_pct_Ub"] = 100 * met["L2_U"]
        out["metrics"][f"{frac:.3f}"] = met
        ax[0, j].set_title(rf"$x_1/\lambda={frac:.3f}$", fontsize=10)
        ax[2, j].set_xlabel(r"$\overline{u_iu_j}/U_b^2$")
        for i in range(3):
            ax[i, j].set_ylim(0, 1.08)
            plotstyle.panel_label(ax[i, j], f"({'abc'[i]}{j + 1})")
            if j:
                ax[i, j].tick_params(labelleft=False)
        ax[0, j].set_xlabel(r"$U/U_b$")
        ax[1, j].set_xlabel(r"$V/U_b$")
        ax[2, j].set_xlim(left=-0.003)
    for i in range(3):
        ax[i, 0].set_ylabel(r"$x_3/H$")
    ax[0, 0].legend(loc="center left", fontsize=7.5)
    ax[2, 0].legend(loc="upper right", fontsize=7.5)
    out["max_uu_dns"] = float(ymax_stress)
    fig.savefig(os.path.join(a.out, "wavywall_validation.png"), dpi=a.dpi)
    with open(os.path.join(a.out, "wavywall_validation.json"), "w") as fh:
        json.dump(out, fh, separators=(",", ":"))
    for f, mt in out["metrics"].items():
        print(f"x/lambda={f}: L2(U)={mt['L2_U']:.4f} Ub, L2(V)={mt['L2_V']:.4f} Ub, "
              f"L2(uu)={mt['L2_uu_rel']:.1f}% of DNS max, L2(-uw)={mt['L2_neg_uw_rel']:.1f}%")


if __name__ == "__main__":
    main()
