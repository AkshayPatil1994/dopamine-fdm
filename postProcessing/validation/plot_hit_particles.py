#!/usr/bin/env python3
"""plot_hit_particles.py -- manuscript figures of the point-particle validation in forced HIT.

Cases come from validation/HIT/particles/make_cases.py (tracer, heavy_St<S>_g<G>, light_rho<R>_St<S>_g<G>), one-way coupled, the same
flow realisation in every case. Snapshots fields/hit_particles.<step> every --dt-snap time units (periodic box of side 2*pi).
Output (png, plotstyle.py style):
  fig_particles_tracer.png    (a) Lagrangian velocity autocorrelation and its integral T_L, (b) single-particle dispersion
  fig_particles_variance.png  heavy particles without gravity: <v'^2>/<u'^2> against St with the Tchen-Hinze limit 1/(1+tau_p/T_L)
  fig_particles_settling.png  <v_y>/v_t against St for heavy (settling) and light (rise) particles; v_t = (1-beta) g tau_p is the
                              still-fluid Maxey-Riley terminal velocity, beta = 3 rho_f/(rho_f + 2 rho_p)
  fig_particles_rdf.png       radial distribution function g(r) (preferential concentration)
  particles_summary.txt / .json
Usage: plot_hit_particles.py --root CASES_DIR --out DIR [--nu 0.005 --eps 0.1 --t-spin 10]
"""
import argparse
import glob
import json
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
import plotstyle as ps                         # noqa: E402
from hit_plotutil import run_style             # noqa: E402
import matplotlib.pyplot as plt                          # noqa: E402
from scipy.optimize import brentq                        # noqa: E402
from scipy.spatial import cKDTree                        # noqa: E402
from dopamine_post.particles import read_snapshot        # noqa: E402

TWO_PI = 2 * np.pi


def load(case, dt_snap):
    files = sorted(glob.glob(str(Path(case) / "fields" / "hit_particles.*")), key=lambda f: int(f.rsplit(".", 1)[1]))
    pos, vel = [], []
    for f in files:
        i, d = read_snapshot(f)
        o = np.argsort(i)
        pos.append(d[o, 0:3]); vel.append(d[o, 3:6])
    if not pos:
        return None
    return dt_snap * np.arange(1, len(pos) + 1), np.array(pos), np.array(vel)


def case_params(case):
    txt = (Path(case) / "input_parameters").read_text()
    g = lambda k, d: float(m.group(1)) if (m := re.search(rf"{k}\s*=\s*([-+0-9.eEdD]+)", txt)) else d
    return g("particle_diam", np.nan), int(g("particle_stokes_drag", 0))


def terminal_velocity(g, rho, tau, d, nu, stokes):
    """signed still-fluid terminal velocity along y (heavy: negative): Stokes drag, or Schiller-Naumann as in the solver"""
    v0 = (1 - 3.0 / (1 + 2 * rho)) * g * tau
    if stokes or not np.isfinite(d) or v0 == 0:
        return -v0
    v = brentq(lambda x: x * (1 + 0.15 * (x * d / nu) ** 0.687) - abs(v0), 0, abs(v0))
    return -np.sign(v0) * v


def block_sem(x, n):
    """standard error of the mean of a time series from non-overlapping blocks of n samples"""
    nb = len(x) // n
    b = x[: nb * n].reshape(nb, n).mean(axis=1)
    return b.std(ddof=1) / np.sqrt(nb) if nb > 1 else np.nan


def unwrap(pos):
    inc = np.diff(pos, axis=0)
    inc -= TWO_PI * np.round(inc / TWO_PI)
    return np.concatenate([pos[:1], pos[:1] + np.cumsum(inc, axis=0)])


def rdf(pos, r):
    n = len(pos)
    tree = cKDTree(np.mod(pos, TWO_PI), boxsize=TWO_PI)
    cnt = tree.count_neighbors(tree, r, cumulative=True) - n
    return cnt / (n * (n - 1) * 4.0 / 3.0 * np.pi * r ** 3 / TWO_PI ** 3)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--nu", type=float, default=0.005)
    ap.add_argument("--eps", type=float, default=0.1)
    ap.add_argument("--t-spin", type=float, default=10.0, help="discard this long after the particle release")
    ap.add_argument("--dt-snap", type=float, default=0.25)
    ap.add_argument("--ref-dir", default="/mnt/storage1/fdm-dopamine/validation/HIT/reference_data",
                    help="directory with ireland2016_fig5_Rl90_digitised.csv")
    ap.add_argument("--rdf-ref", default="/mnt/storage1/fdm-dopamine/validation/HIT/reference_data/ireland2016_fig20_Rl88_rdf_digitised.csv")
    ap.add_argument("--nsub", type=int, default=20000, help="particles used for the pair counting")
    a = ap.parse_args()
    ps.apply()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    tau_eta = np.sqrt(a.nu / a.eps)
    eta = (a.nu ** 3 / a.eps) ** 0.25
    u_eta = (a.nu * a.eps) ** 0.25
    rr = eta * np.array([0.5, 1, 2, 4, 8, 16])

    cases = sorted(p.name for p in Path(a.root).iterdir() if (p / "input_parameters").exists())
    data = {c: load(Path(a.root) / c, a.dt_snap) for c in cases}
    summary, js = [], {}
    T_L, urms = np.nan, np.nan

    if data.get("tracer"):
        t, pos, vel = data["tracer"]
        sel = t >= a.t_spin
        v = vel[sel]
        u2 = (v ** 2).mean()
        urms = np.sqrt(u2)
        lags = np.arange(0, min(len(v) // 2, int(10 / a.dt_snap)))
        rho = np.array([(v[: len(v) - l] * v[l:]).mean() / u2 for l in lags])
        cut = int(np.argmax(rho < 0)) or len(rho)
        T_L = np.trapz(rho[:cut], dx=a.dt_snap)
        up = unwrap(pos[sel])
        msd = np.array([((up[l:] - up[: len(up) - l]) ** 2).sum(-1).mean() for l in lags])
        tt = lags * a.dt_snap
        fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.0))
        ax[0].plot(tt / T_L, rho, "o", color=ps.BLUE, markevery=max(1, len(tt) // 25), label="DNS tracers")
        ax[0].plot(tt / T_L, np.exp(-tt / T_L), color=ps.GREY, ls=":", label=r"$e^{-\tau/T_L}$")
        ax[0].set_xlabel(r"$\tau/T_L$"); ax[0].set_ylabel(r"$\rho_L(\tau)$"); ax[0].legend()
        ax[1].loglog(tt[1:] / T_L, msd[1:] / (3 * u2 * T_L ** 2), "o", color=ps.BLUE, markevery=max(1, len(tt) // 25), label="DNS tracers")
        ax[1].loglog(tt[1:] / T_L, (tt[1:] / T_L) ** 2, color=ps.GREY, ls="--", label="ballistic")
        ax[1].loglog(tt[1:] / T_L, 2 * tt[1:] / T_L, color=ps.GREY, ls=":", label=r"diffusive $2\tau/T_L$")
        ax[1].set_xlabel(r"$\tau/T_L$"); ax[1].set_ylabel(r"$\langle|\Delta x|^2\rangle/(3u'^2T_L^2)$"); ax[1].legend()
        for x, pl in zip(ax, "ab"):
            ps.panel_label(x, f"({pl})")
        fig.tight_layout(); ps.save(fig, out / "fig_particles_tracer.png"); plt.close(fig)
        C0 = 2 * u2 / (a.eps * T_L)
        summary.append(f"tracer: u' = {urms:.4f}, T_L = {T_L:.3f} ({T_L / tau_eta:.2f} tau_eta), C0 = 2u'^2/(eps T_L) = {C0:.2f}, "
                       f"D = u'^2 T_L = {u2 * T_L:.4f}")
        js["tracer"] = dict(urms=urms, T_L=T_L, C0=C0)

    res = []
    for c in cases:
        if c == "tracer" or data[c] is None:
            continue
        m = re.match(r"(heavy|light)_(?:rho([0-9.]+)_)?St([0-9.]+)_(?:Fr[0-9.]+_)?g([0-9.]+)", c)
        kind, rho, St, g = m.group(1), float(m.group(2) or 1000.0), float(m.group(3)), float(m.group(4))
        t, pos, vel = data[c]
        sel = t >= a.t_spin
        tau = St * tau_eta
        diam, stokes = case_params(Path(a.root) / c)
        vt = terminal_velocity(g, rho, tau, diam, a.nu, stokes)
        vy = vel[sel][..., 1].mean()
        nblk = max(1, int(np.ceil(2 * T_L / a.dt_snap)))
        sem = block_sem(vel[sel][..., 1].mean(axis=1), nblk)
        v2 = ((vel[sel] - vel[sel].mean(axis=(0, 1))) ** 2).mean() / urms ** 2
        step = max(1, int(sel.sum()) // 20)
        g_r = np.mean([rdf(pos[sel][k][: a.nsub], rr) for k in range(0, int(sel.sum()), step)], axis=0)
        th = 1.0 / (1.0 + tau / T_L)
        res.append(dict(case=c, fr="_Fr" in c, dv=(vt - vy) / u_eta, kind=kind, rho=rho, St=St, g=g, ratio=vy / vt if g > 0 else np.nan, v2=v2, th=th, grr=g_r, vy=vy, vt=vt,
                        sem=sem, d_eta=diam / eta, stokes=stokes))
        summary.append(f"{c:30s} <v_y> = {vy:+.4f} (+-{sem:.4f})  v_t = {vt:+.4f} ({'Stokes' if stokes else 'Schiller-Naumann'})  "
                       f"ratio = {(f'{vy / vt:.3f}' if g > 0 else '--')}  d/eta = {diam / eta:.2f}  "
                       f"<v'^2>/<u'^2> = {v2:.3f}  Tchen-Hinze = {th:.3f}  g(r/eta=0.5,1,2,4,8,16) = {np.round(g_r, 2)}")
        js[c] = dict(St=St, rho=rho, g=g, vy=vy, vy_sem=sem, vt=vt, stokes_drag=stokes, d_over_eta=diam / eta, ratio=res[-1]["ratio"], v2=v2, tchen_hinze=th, g_r=list(map(float, g_r)))
    (out / "particles_summary.txt").write_text("\n".join(summary) + "\n")
    (out / "particles_summary.json").write_text(json.dumps(js, indent=1, default=float))
    print("\n".join(summary))
    if not res:
        return

    # variance (heavy, g = 0) against Tchen-Hinze
    h0 = sorted((r for r in res if r["kind"] == "heavy" and r["g"] == 0), key=lambda r: r["St"])
    if h0:
        fig, ax = plt.subplots(figsize=(3.6, 3.0))
        Sts = np.logspace(-1.2, 0.7, 100)
        ax.semilogx(Sts, 1 / (1 + Sts * tau_eta / T_L), color=ps.GREY, ls="--", label="Tchen-Hinze")
        ax.semilogx([r["St"] for r in h0], [r["v2"] for r in h0], "o", color=ps.BLUE, label="DNS")
        ax.set_xlabel(r"$St_\eta$"); ax.set_ylabel(r"$\langle v'^2\rangle/\langle u'^2\rangle$"); ax.legend()
        fig.tight_layout(); ps.save(fig, out / "fig_particles_variance.png"); plt.close(fig)

    # same seeds and flow in every case: subtracting the g = 0 mean of the same species removes the sampling offset of <v_y>
    v0 = {(r["rho"], r["St"]): r["vy"] for r in res if r["g"] == 0 and r["kind"] == "heavy"}
    for r in res:
        key = (r["rho"], r["St"])
        r["ratio_c"] = (r["vy"] - v0[key]) / r["vt"] if r["g"] > 0 and not r["fr"] and key in v0 else r["ratio"]

    # settling / rise
    fr = sorted((r for r in res if r["fr"]), key=lambda r: r["St"])
    fig, ax = plt.subplots(1, 2 if fr else 1, figsize=(7.4 if fr else 3.8, 3.0), squeeze=False)
    ax = ax[0]
    groups = sorted({(r["kind"], r["rho"], r["g"]) for r in res if r["g"] > 0 and not r["fr"]})
    for i, (kind, rho, g) in enumerate(groups):
        pts = sorted((r["St"], r["ratio_c"]) for r in res if (r["kind"], r["rho"], r["g"]) == (kind, rho, g) and not r["fr"])
        lab = (f"heavy, $g={g:g}$" if kind == "heavy" else rf"light $\rho_p/\rho_f={rho:g}$, $g={g:g}$")
        ax[0].semilogx([p[0] for p in pts], [p[1] for p in pts], "o" if kind == "heavy" else "s",
                       color=ps.FALLBACK[i % len(ps.FALLBACK)], label=lab)
    ax[0].axhline(1.0, color=ps.GREY, ls=":")
    ax[0].set_xlabel(r"$St_\eta$"); ax[0].set_ylabel(r"$(\langle v_y\rangle-\langle v_y\rangle_{g=0})/v_t$"); ax[0].legend(fontsize=7)
    ps.panel_label(ax[0], "(a)")
    if fr:
        stk = all(r["stokes"] for r in fr)
        ax[1].errorbar([r["St"] for r in fr], [r["dv"] for r in fr], yerr=[r["sem"] / u_eta for r in fr], fmt="o", color=ps.BLUE,
                       capsize=2, label=rf"present, $Fr=0.052$, {'Stokes' if stk else 'Schiller-Naumann'} drag")
        ax[1].set_xscale("log")
        ref = Path(a.ref_dir) / "ireland2016_fig5_Rl90_digitised.csv"
        if ref.exists():
            d = np.array([[float(x) for x in ln.split(",")] for ln in ref.read_text().splitlines()
                          if ln.strip() and not ln.startswith("#") and not ln.startswith("St")])
            ax[1].semilogx(d[:, 0], d[:, 1], "s", color=ps.ORANGE, label=r"Ireland et al. (2016), $R_\lambda=90$ (digitised)")
        ax[1].axhline(0.0, color=ps.GREY, ls=":")
        ax[1].set_xlabel(r"$St_\eta$"); ax[1].set_ylabel(r"$(\langle v_y\rangle-v_t)/u_\eta$ (settling enhancement)"); ax[1].legend(fontsize=7)
        ps.panel_label(ax[1], "(b)")
    fig.tight_layout(); ps.save(fig, out / "fig_particles_settling.png"); plt.close(fig)

    # radial distribution functions
    fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.0))
    for i, r in enumerate(res):
        if r["g"] == 0 or r["kind"] == "light":
            sty = run_style(r["case"], i)
            ax[0 if r["kind"] == "heavy" else 1].loglog(rr / eta, r["grr"], ls="none", marker=sty["marker"], color=sty["color"], label=r["case"])
    ref = Path(a.rdf_ref)
    if ref.exists():
        d = np.array([[float(x) for x in ln.split(",")] for ln in ref.read_text().splitlines()
                      if ln.strip() and not ln.startswith("#") and not ln.startswith("St")])
        ax[0].loglog(d[:, 1], d[:, 2], "k*", ms=8, label=r"Ireland et al. (2016) $R_\lambda=88$, St=1, 3 (digitised)")
    for x, pl, t in zip(ax, "ab", ("heavy, no gravity", "light, gravity (finite size, qualitative)")):
        x.set_xlabel(r"$r/\eta$"); x.set_ylabel(r"$g(r)$"); x.legend(fontsize=6); ps.panel_label(x, f"({pl}) {t}")
    fig.tight_layout(); ps.save(fig, out / "fig_particles_rdf.png"); plt.close(fig)


if __name__ == "__main__":
    main()
