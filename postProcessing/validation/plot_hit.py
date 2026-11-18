#!/usr/bin/env python3
"""plot_hit.py -- manuscript figures of the forced-HIT validation (finite-difference solver vs pseudo-spectral reference).

Cases (label=directory): the spectral reference (hit_spectral_reference.py output: hit_stats.dat, spectrum.*.dat) and any number of
finite-difference runs (hit_stats.dat + fields/hit.<step>). Window means over [--t0, --t1]. Output (png, plotstyle.py style):
  fig_hit_stationarity.png  TKE, eps/eps_in and Re_lambda against time
  fig_hit_spectra.png       (a) E(k) with Pope's model, (b) compensated spectrum, (c) Kolmogorov-scaled E(k eta) with the JHTDB
                            Re_lambda = 433 forced-isotropic spectrum (reference_data/jhtdb_spectrum.txt)
  fig_hit_convergence.png   (a) E_FD(k)/E_spectral(k), (b) window statistics against the resolution
  hit_summary.txt / hit_summary.json   window means (standard error) of the integral quantities, skewness and flatness
Usage: plot_hit.py --out DIR --nu NU --eps-in EPS --spectral spectral=DIR [--ref-dir DIR] [--t0 30 --t1 1e9] label=DIR ...
"""
import argparse
import glob
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
import plotstyle as ps                         # noqa: E402
from hit_plotutil import run_style, data_style, log_marks   # noqa: E402
import matplotlib.pyplot as plt                # noqa: E402
from dopamine_post import hit                  # noqa: E402

JHTDB_NU, JHTDB_EPS = 1.85e-4, 0.0928          # JHTDB forced isotropic 1024^3, README (t = 0 .. 2.048)


def fd_spectrum(case, nu, t0, t1, st):
    ks, E, mom, n = None, 0.0, [], 0
    for f in sorted(glob.glob(str(Path(case) / "fields" / "hit.*"))):
        step = int(f.rsplit(".", 1)[1])
        t = float(np.interp(step, st["step"], st["t"]))
        if not (t0 <= t <= t1):
            continue
        v = hit.read_velocity(f)
        k, e, _, _ = hit.spectrum(v, nu)
        ks, E, n = k, E + e, n + 1
        mom.append(hit.gradient_moments(v))
    return (None, None, None) if n == 0 else (ks, E / n, np.mean(mom, axis=0))


def spectral_moments(d):
    f = Path(d) / "fields.final.npz"
    if not f.exists():
        return None
    z = np.load(f)
    N = z["u"].shape[0]
    n = np.fft.fftfreq(N, 1.0 / N)
    kk = np.meshgrid(n, n, n, indexing="ij")
    g = [np.real(np.fft.ifftn(1j * kk[i] * np.fft.fftn(z[c]))) for i, c in enumerate("uvw")]
    return (np.mean([np.mean(x ** 3) / np.mean(x ** 2) ** 1.5 for x in g]),
            np.mean([np.mean(x ** 4) / np.mean(x ** 2) ** 2 for x in g]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cases", nargs="+", help="label=directory (finite-difference runs)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--nu", type=float, required=True)
    ap.add_argument("--eps-in", type=float, required=True)
    ap.add_argument("--spectral", default=None,
                    help="label=directory of the spectral reference; without it the first case is the reference (black line)")
    ap.add_argument("--ref-dir", default="/mnt/storage1/fdm-dopamine/validation/HIT/reference_data",
                    help="directory with jhtdb_spectrum.txt (downloaded by validation/HIT/fetch_reference_data.sh)")
    ap.add_argument("--t0", type=float, default=30.0)
    ap.add_argument("--t1", type=float, default=1e30)
    a = ap.parse_args()
    ps.apply()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)

    fds = [c.split("=", 1) for c in a.cases]
    rows = []
    for label, d in fds:
        st = hit.read_stats(Path(d) / "hit_stats.dat")
        k, E, mom = fd_spectrum(d, a.nu, a.t0, a.t1, st)
        if k is None:
            print(f"{label}: no snapshot in the window, skipped")
            continue
        rows.append((label, hit.time_average(st, a.t0, a.t1), mom, hit.scales(k, E, a.nu), k, E, st))
    if a.spectral:
        sp_label, sp_dir = a.spectral.split("=", 1)
        sp_stats, sp_spec = hit.read_spectral(sp_dir)
        sel = [s for s in sp_spec if a.t0 <= s[0] <= a.t1]
        k_sp, E_sp = sel[0][1], np.mean([s[2] for s in sel], axis=0)
        kcut = (2 * (len(k_sp) - 1)) // 3 - 1                    # shells fully inside the 2/3 truncation
        k_sp, E_sp = k_sp[: kcut + 1], E_sp[: kcut + 1]
        rows.insert(0, (sp_label, hit.time_average(sp_stats, a.t0, a.t1), spectral_moments(sp_dir),
                        hit.scales(k_sp, E_sp, a.nu), k_sp, E_sp, sp_stats))
    k_sp, E_sp = rows[0][4], rows[0][5]

    # ── stationarity
    fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.0))
    for i, (label, avg, mom, sc, k, E, st) in enumerate(rows):
        sty = data_style(i - 1, label, reference=(i == 0 and bool(a.spectral)))
        me = max(1, len(st["t"]) // 35)
        ax[0].plot(st["t"], st["tke"], markevery=me, label=label, **sty)
        ax[1].plot(st["t"], st["eps"] / a.eps_in, markevery=me, **sty)
        ax[2].plot(st["t"], st["rel"], markevery=me, **sty)
    ax[1].axhline(1.0, color=ps.GREY, ls=":", lw=1)
    for x, yl, pl in zip(ax, (r"$K$", r"$\epsilon/\epsilon_{in}$", r"$Re_\lambda$"), "abc"):
        x.set_xlabel(r"$t$"); x.set_ylabel(yl); x.axvspan(a.t0, x.get_xlim()[1], color="0.93", zorder=0); ps.panel_label(x, f"({pl})")
    ax[0].legend()
    fig.tight_layout(); ps.save(fig, out / "fig_hit_stationarity.png"); plt.close(fig)

    # ── spectra
    ref = rows[0][3]
    kk = np.logspace(-0.2, 2.2, 400)
    pm = hit.pope_spectrum(kk, ref["eps"], a.nu, ref["tke"])
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.4))
    for i, (label, avg, mom, sc, k, E, st) in enumerate(rows):
        sty = data_style(i - 1, label, reference=(i == 0 and bool(a.spectral)))
        me = log_marks(len(k) - 1)
        ax[0].loglog(k[1:], E[1:], markevery=me, label=label, **sty)
        ax[1].loglog(k[1:], k[1:] ** (5 / 3) * E[1:] / sc["eps"] ** (2 / 3), markevery=me, **sty)
        ax[2].loglog(k[1:] * sc["eta"], E[1:] / (a.nu ** 5 * sc["eps"]) ** 0.25, markevery=me, **sty)
    ax[0].loglog(kk, pm, color=ps.GREY, ls="--", lw=1.0, label="Pope model")
    ax[1].loglog(kk, kk ** (5 / 3) * pm / ref["eps"] ** (2 / 3), color=ps.GREY, ls="--", lw=1.0)
    ax[1].axhline(1.5, color=ps.GREY, ls=":", lw=0.8)
    etak = np.logspace(-3, 0.3, 300)
    pk = etak / ref["eta"]
    ax[2].loglog(etak, hit.pope_spectrum(pk, ref["eps"], a.nu, ref["tke"]) / (a.nu ** 5 * ref["eps"]) ** 0.25,
                 color=ps.GREY, ls="--", lw=1.0)
    jh = Path(a.ref_dir) / "jhtdb_spectrum.txt"
    if jh.exists():
        d = np.loadtxt(jh, skiprows=2)
        d = d[d[:, 0] <= 2 * 512 / 3 - 1]                       # inside the 2/3 truncation of the 1024^3 run
        eta_j = (JHTDB_NU ** 3 / JHTDB_EPS) ** 0.25
        ax[2].loglog(d[:, 0] * eta_j, d[:, 1] / (JHTDB_NU ** 5 * JHTDB_EPS) ** 0.25, color="k", ls=":", lw=1.0,
                     label=r"JHTDB $Re_\lambda=433$")
        ax[2].legend()
    ax[0].set_xlabel(r"$k$"); ax[0].set_ylabel(r"$E(k)$"); ax[0].legend()
    ax[1].set_xlabel(r"$k$"); ax[1].set_ylabel(r"$k^{5/3}E/\epsilon^{2/3}$"); ax[1].set_ylim(0.05, 4)
    ax[2].set_xlabel(r"$k\eta$"); ax[2].set_ylabel(r"$E(k)/(\epsilon\nu^5)^{1/4}$"); ax[2].set_xlim(2e-2, 2); ax[2].set_ylim(1e-6, 1e2)
    for x, pl in zip(ax, "abc"):
        ps.panel_label(x, f"({pl})")
    fig.tight_layout(); ps.save(fig, out / "fig_hit_spectra.png"); plt.close(fig)

    # ── convergence against the spectral reference
    fig, ax = plt.subplots(1, 2, figsize=(8.2, 3.2))
    names, vals = [], {q: [] for q in ("tke", "eps", "rel", "eta_dx")}
    for i, (label, avg, mom, sc, k, E, st) in enumerate(rows[1:]):
        n = min(len(k), len(k_sp))
        n = min(n, int(np.argmax(E_sp[1:] <= 1e-12)) + 1 if np.any(E_sp[1:] <= 1e-12) else n)   # spectral 2/3 truncation
        sty = data_style(i, label)
        ax[0].semilogy(k[1:n] * ref["eta"], E[1:n] / E_sp[1:n], markevery=max(1, n // 18), label=label, **sty)
        names.append(label)
        for q in vals:
            vals[q].append(avg[q][0])
    ax[0].axhline(1.0, color="k", lw=0.8)
    ax[0].set_xlabel(r"$k\eta$"); ax[0].set_ylabel(r"$E_{FD}(k)/E_{spectral}(k)$"); ax[0].set_ylim(0.3, 3); ax[0].legend()
    x = np.arange(len(names))
    for q, mk, lab in (("tke", "o", r"$K$"), ("eps", "s", r"$\epsilon$"), ("rel", "^", r"$Re_\lambda$")):
        ref_v = rows[0][1][q][0]
        ax[1].plot(x, np.array(vals[q]) / ref_v, mk, ls="none", label=lab)
    ax[1].axhline(1.0, color="k", lw=0.8); ax[1].set_xticks(x); ax[1].set_xticklabels(names, rotation=20)
    ax[1].set_ylabel("FD / spectral"); ax[1].legend()
    for x_, pl in zip(ax, "ab"):
        ps.panel_label(x_, f"({pl})")
    fig.tight_layout(); ps.save(fig, out / "fig_hit_convergence.png"); plt.close(fig)

    # ── table
    lines = [f"window t in [{a.t0}, {a.t1}], nu = {a.nu}, eps_in = {a.eps_in}; mean (standard error of the mean)",
             f"{'case':20s}{'TKE':>17s}{'eps/eps_in':>12s}{'u_rms':>9s}{'Re_lam':>9s}{'eta/dx':>9s}{'L11':>8s}{'T_e':>8s}"
             f"{'skew':>8s}{'flat':>7s}"]
    js = {}
    for label, avg, mom, sc, k, E, st in rows:
        sk, fl = (f"{mom[0]:.3f}", f"{mom[1]:.2f}") if mom is not None else ("-", "-")
        lines.append(f"{label:20s}{avg['tke'][0]:>10.4f}({avg['tke'][1]:.4f}){avg['eps'][0] / a.eps_in:>12.4f}"
                     f"{avg['urms'][0]:>9.4f}{avg['rel'][0]:>9.2f}{avg['eta_dx'][0]:>9.3f}{sc['L11']:>8.3f}"
                     f"{avg['Te'][0]:>8.3f}{sk:>8s}{fl:>7s}")
        js[label] = dict(tke=avg["tke"], eps=avg["eps"], urms=avg["urms"], rel=avg["rel"], eta_dx=avg["eta_dx"],
                         L11=sc["L11"], skew=None if mom is None else mom[0], flat=None if mom is None else mom[1])
    (out / "hit_summary.txt").write_text("\n".join(lines) + "\n")
    (out / "hit_summary.json").write_text(json.dumps(js, indent=1, default=float))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
