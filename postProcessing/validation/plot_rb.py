#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
"""Rayleigh-Benard convection (Ra=2.5e6, Pr=1): mean temperature against Kunnen, Geurts & Clercx (2009).

Usage: python3 plot_rb.py --case /mnt/storage1/fdm-dopamine/validation/RBconvection --out docs/assets/validation/rb
The profile is the x-z and time average of the Boussinesq temperature over snapshots with t >= --t-start.
"""
import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np

import plotstyle

plotstyle.apply()
import matplotlib.pyplot as plt  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--t-start", type=float, default=75.0)
    a = ap.parse_args()
    case = Path(a.case)
    sys.path.insert(0, str(case))
    import snapshot_io as sio  # noqa: E402

    prm = sio.parse_input_parameters(case / "input_parameters")
    tsave = float(prm.get("tsave", 0.5))
    pre = str(prm.get("fileout", "RBconv_Kunnen"))
    snaps = sio.list_snapshots(str(case / "fields"), pre)[int(a.t_start / tsave):]
    acc, ym = 0.0, None
    for _, fp in snaps:
        s = sio.read_snapshot(fp, sgs_model=int(prm.get("sgs_model", 0)), sediment_flag=0, boussinesq_flag=1)
        ym = s["ym"]
        acc = acc + s["T"].mean(axis=(0, 2))
    T = acc / len(snaps)
    ref = np.loadtxt(case / "ref_data.dat", comments="#", delimiter=",")
    Tr, zr = ref[:, 0], ref[:, 1]
    m = (zr >= ym[0]) & (zr <= ym[-1])
    res = np.interp(zr[m], ym, T) - Tr[m]
    out = {"n_snapshots": len(snaps), "avg_time": len(snaps) * tsave, "n_ref": int(m.sum()),
           "L2": float(np.sqrt(np.mean(res ** 2))), "Linf": float(np.abs(res).max()),
           "z_Linf": float(zr[m][np.argmax(np.abs(res))]), "T_mid": float(np.interp(0.5, ym, T)),
           "T_mid_ref": float(np.interp(0.5, zr[::-1], Tr[::-1])),
           "z": ym[::2].tolist(), "T": T[::2].tolist()}

    fig, ax = plt.subplots(1, 2, figsize=(7.0, 3.8), gridspec_kw={"width_ratios": [1.6, 1]}, constrained_layout=True)
    ax[0].plot(Tr, zr, "k-", lw=1.6, label="Kunnen et al. (2009)")
    ax[0].plot(T, ym, marker="o", ls="none", color=plotstyle.BLUE, markersize=3.6, label="DOPAMINE")
    ax[0].set_xlabel(r"$\langle T\rangle$")
    ax[0].set_ylabel(r"$z/H$")
    ax[0].legend(loc="upper right")
    ax[1].plot(res, zr[m], color=plotstyle.BLUE, lw=1.2)
    ax[1].axvline(0, color=plotstyle.GREY, lw=0.8, ls=":")
    ax[1].set_xlabel(r"$\langle T\rangle-\langle T\rangle_\mathrm{ref}$")
    ax[1].tick_params(labelleft=False)
    for x, s in zip(ax, "ab"):
        x.set_ylim(0, 1)
        plotstyle.panel_label(x, f"({s})")
    os.makedirs(a.out, exist_ok=True)
    plotstyle.save(fig, os.path.join(a.out, "rb_validation.png"))
    with open(os.path.join(a.out, "rb_validation.json"), "w") as fh:
        json.dump(out, fh, separators=(",", ":"))
    print(f"{len(snaps)} snapshots (t_avg={out['avg_time']}), L2={out['L2']:.4f}, Linf={out['Linf']:.4f} at z={out['z_Linf']:.3f}, "
          f"T(0.5)={out['T_mid']:.4f}")


if __name__ == "__main__":
    main()
