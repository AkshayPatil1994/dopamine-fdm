"""
rsb.py -- Reynolds-stress budget post-processing for fdm-dopamine.

Reads input_parameters (Fortran namelist, via `_core.parse_input_parameters`)
and stats/rsb.meta to obtain all simulation and output-grid parameters, loads
every RSB binary (mean velocity, Reynolds stresses, and every budget term),
and plots the temporal mean and standard deviation across available samples
against the MKM chan590 DNS reference (Moser, Kim & Mansour 1999,
Re_tau ~= 587), exploiting channel centre-plane symmetry to double the
effective sample count.

`rsb.meta`'s keys (out_nx, out_ny, out_nz, nsamples, ...) don't fit
`_core.read_meta`'s schema (written for probe_output.f90's line/slice-probe
meta files, which requires an npts/n1 key) -- this module keeps its own
small meta parser rather than force-fitting that helper.

Library use
-----------
    from dopamine_post.rsb import RSBStats

    stats = RSBStats.read(".")                 # reads stats/, input_parameters
    stats.plot()                                # all samples -> plots/*.png
    stats.plot(last=20)                         # last 20 snapshots only

Output (written under <case>/plots by default)
------------------------------------------------
    mean_velocity.png          U+ vs y+  (log-law overlay)
    reynolds_stresses.png      all Rij components vs y+
    budget_11.png  ...  budget_23.png   full stress budget per component
    budget_overview.png        2x3 panel overview of all budgets
"""
import io
import urllib.request
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from . import _core

# ── MKM DNS reference (Moser, Kim & Mansour 1999, Re_tau ~= 587) ─────────────
MKM_BASE = "https://turbulence.oden.utexas.edu/data/MKM/chan590/profiles"

# ── stress-component metadata ─────────────────────────────────────────────────
COMP_NAMES = ["11", "22", "33", "12", "13", "23"]
COMP_LABELS = [
    r"$\langle u'u' \rangle$",
    r"$\langle v'v' \rangle$",
    r"$\langle w'w' \rangle$",
    r"$\langle u'v' \rangle$",
    r"$\langle u'w' \rangle$",
    r"$\langle v'w' \rangle$",
]

# symmetry signs under the reflection  y -> Ly - y
# (+1 = symmetric, -1 = anti-symmetric)
SYM_UMEAN = np.array([+1., -1., +1.])          # U, V, W
SYM_RIJ = np.array([+1., +1., +1., -1., +1., -1.])   # R11 R22 R33 R12 R13 R23

# (file-stem, display label, matplotlib color)
BUDGET_TERMS = [
    ("Pij",      r"Production $P_{ij}$",                       "C0"),
    ("epsRes",   r"Resolved dissipation $-\varepsilon^r_{ij}$", "C3"),
    ("epsSGS",   r"SGS dissipation $-\varepsilon^{sgs}_{ij}$", "C4"),
    ("PiStrain", r"Pressure-strain $\Pi_{ij}$",                "C1"),
    ("DTij",     r"Turb. transport $D^T_{ij}$",                "C2"),
    ("Dnuij",    r"Viscous diffusion $D^\nu_{ij}$",            "C5"),
    ("PhiPij",   r"Pressure diffusion $D^p_{ij}$",             "C6"),
    ("Resid",    r"Residual",                                   "0.55"),
]


# ── rsb.meta parsing (doesn't fit _core.read_meta's probe-meta schema) ───────
def parse_meta(filepath):
    """Parse rsb.meta into a dict; integer and float values are converted."""
    meta = {}
    with open(filepath) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("#") or not line:
                continue
            key, _, val = line.partition("=")
            key = key.strip().lower()
            val = val.strip()
            for cast in (int, float):
                try:
                    meta[key] = cast(val)
                    break
                except ValueError:
                    pass
            else:
                meta[key] = val
    return meta


# ── MKM reference loader ─────────────────────────────────────────────────────
def load_mkm_reference(cache_dir):
    """Load MKM chan590 mean-velocity and Reynolds-stress profiles.

    On first call the two raw data files are downloaded from the JHTDB server
    and cached as CSVs under `cache_dir`; subsequent calls read the cached
    CSVs directly, avoiding repeated network access.

    Both files use wall-unit normalisation (U_tau, h), so the returned
    arrays are already in '+' units and can be plotted directly.

    Returns dict: yp_means, Umean, yp_stress, R_uu, R_vv, R_ww, R_uv, R_uw, R_vw.
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    def _load(name):
        cache_path = cache_dir / (name + ".csv")
        if cache_path.exists():
            print(f"  MKM cache hit: {cache_path}")
            return np.loadtxt(cache_path, delimiter=",")
        url = f"{MKM_BASE}/{name}"
        print(f"  Downloading MKM reference: {url}")
        with urllib.request.urlopen(url) as resp:
            text = resp.read().decode()
        data_lines = [ln for ln in text.splitlines()
                      if ln.strip() and not ln.strip().startswith("#")]
        arr = np.loadtxt(io.StringIO("\n".join(data_lines)))
        np.savetxt(cache_path, arr, delimiter=",")
        print(f"  Cached to {cache_path}")
        return arr

    means = _load("chan590.means")        # y  y+  Umean  dU/dy  Wmean  dW/dy  P
    stress = _load("chan590.reystress")   # y  y+  Ruu  Rvv  Rww  Ruv  Ruw  Rvw

    return {
        "yp_means":  means[:, 1],
        "Umean":     means[:, 2],
        "yp_stress": stress[:, 1],
        "R_uu":      stress[:, 2],
        "R_vv":      stress[:, 3],
        "R_ww":      stress[:, 4],
        "R_uv":      stress[:, 5],
        "R_uw":      stress[:, 6],
        "R_vw":      stress[:, 7],
    }


# ── grid helper ───────────────────────────────────────────────────────────────
def tanh_grid_faces(Ly, ny, alpha):
    """Face y-coordinates for grid_type=2 (tanh clustering at both walls).

    Mapping:  y_j = (Ly/2) * [1 + tanh(alpha*(j/(ny-1) - 1/2)) / tanh(alpha/2)]
    j = 0    ->  y = 0   (lower wall)
    j = ny-1 ->  y = Ly  (upper wall)
    """
    xi = np.arange(ny, dtype=float) / (ny - 1)
    return 0.5 * Ly * (1.0 + np.tanh(alpha * (xi - 0.5)) / np.tanh(0.5 * alpha))


# ── I/O helper ────────────────────────────────────────────────────────────────
def load_rsb(stats_dir, name, nc, out_nx, out_ny, out_nz, nsamples):
    """Load a single RSB binary.

    Disk layout: (nc, out_nx, out_ny, out_nz, nsamples), Fortran-contiguous.
    Returns array of shape (nc, out_ny, nsamples) after squeezing the
    homogenised x and z dimensions (this only supports rsb_hom_dir covering
    both x and z, i.e. out_nx == out_nz == 1).
    """
    path = Path(stats_dir) / f"rsb_{name}.bin"
    raw = np.fromfile(path, dtype="<f8")
    arr = raw.reshape((nc, out_nx, out_ny, out_nz, nsamples), order="F")
    if out_nx != 1 or out_nz != 1:
        raise ValueError(
            f"rsb_{name}.bin: out_nx={out_nx}, out_nz={out_nz} (expected 1,1); "
            "this only supports x and z both in rsb_hom_dir and only reads "
            "plane [0,:,0,:], which would silently discard the other planes.")
    return arr[:, 0, :, 0, :]          # shape (nc, out_ny, nsamples)


# ── stats helpers ─────────────────────────────────────────────────────────────
def _mean_std(arr):
    """Return (mean, std) over the last (sample) axis."""
    return arr.mean(axis=-1), arr.std(axis=-1)


def symmetry_fold(arr, sym_signs):
    """Double the effective sample count by exploiting channel centre-plane
    symmetry: the upper half of the wall-normal extent is reflected onto the
    lower half and appended as additional samples. The sign of each
    component under the reflection y -> Ly - y is given by `sym_signs`.

    arr (nc, out_ny, nsamples) -> (nc, n_half, 2*nsamples), n_half = out_ny//2
    """
    nc, out_ny, nsamples = arr.shape
    n_half = out_ny // 2
    lower = arr[:, :n_half, :]
    upper = arr[:, out_ny - n_half:, :][:, ::-1, :]
    upper = upper * sym_signs[:, np.newaxis, np.newaxis]
    return np.concatenate([lower, upper], axis=-1)


def _shade(ax, x, mean, std, color, alpha=0.18, **kw):
    """Plot a line with a shaded +-1 sigma band."""
    ax.plot(x, mean, color=color, **kw)
    ax.fill_between(x, mean - std, mean + std, color=color, alpha=alpha)


def _log_xaxis(ax, yp):
    ax.set_xscale("log")
    ax.set_xlim(yp.min() * 0.9, yp.max() * 1.05)
    ax.set_xlabel(r"$y^+$", fontsize=12)


class RSBStats:
    """A case's Reynolds-stress budget statistics: mean velocity, Reynolds
    stresses and every budget term, each (component, wall-normal, sample).

    `.plot(last=None)` writes the 9 figures listed in the module docstring;
    `last=N` restricts to the last N snapshots (default: all).
    """

    def __init__(self, case_dir, stats_dir, yp, Umean, Rij, budget_data,
                 u_tau, Re_tau, h, delta_nu, nu, nsamples):
        self.case_dir = Path(case_dir)
        self.stats_dir = Path(stats_dir)
        self.yp = yp                      # (out_ny,) full wall-normal range, both walls
        self.Umean = Umean                # (3, out_ny, nsamples)
        self.Rij = Rij                    # (6, out_ny, nsamples)
        self.budget_data = budget_data    # {stem: (6, out_ny, nsamples)}
        self.u_tau, self.Re_tau = u_tau, Re_tau
        self.h, self.delta_nu, self.nu = h, delta_nu, nu
        self.nsamples = nsamples

    @classmethod
    def read(cls, case_dir=".", stats_dir="stats", input_parameters="input_parameters"):
        """Read a case's stats/ directory (Umean, Rij, and every BUDGET_TERMS
        binary) plus input_parameters/rsb.meta for the grid and wall units."""
        case = Path(case_dir)
        stats = case / stats_dir
        p = _core.parse_input_parameters(case / input_parameters)
        m = parse_meta(stats / "rsb.meta")

        ny = p["ny"]           # face points in wall-normal direction
        Ly = p["ly"]           # full channel height
        alpha = p["alpha_grid"]  # tanh stretching parameter
        nu = p["nu"]            # kinematic viscosity
        dPdx = p["dpdx"]        # imposed streamwise pressure gradient

        out_nx, out_ny, out_nz = m["out_nx"], m["out_ny"], m["out_nz"]
        nsamples = m["nsamples"]

        h = 0.5 * Ly
        u_tau = np.sqrt(dPdx * h)
        delta_nu = nu / u_tau
        Re_tau = h / delta_nu

        y_face = tanh_grid_faces(Ly, ny, alpha)
        y_cc = 0.5 * (y_face[:-1] + y_face[1:])
        yp = y_cc / delta_nu

        Umean = load_rsb(stats, "Umean", 3, out_nx, out_ny, out_nz, nsamples)
        Rij = load_rsb(stats, "Rij", 6, out_nx, out_ny, out_nz, nsamples)
        budget_data = {stem: load_rsb(stats, stem, 6, out_nx, out_ny, out_nz, nsamples)
                        for stem, _, _ in BUDGET_TERMS}

        print("Parsed parameters:")
        print(f"  ny={ny}, Ly={Ly}, alpha_grid={alpha}, nu={nu:.6g}, dPdx={dPdx}")
        print(f"  out_ny={out_ny}, nsamples={nsamples}")
        print(f"  u_tau = {u_tau:.4f},  Re_tau = {Re_tau:.1f}")

        return cls(case, stats, yp, Umean, Rij, budget_data, u_tau, Re_tau, h, delta_nu,
                   nu, nsamples)

    def __repr__(self):
        return (f"RSBStats({self.nsamples} samples, Re_tau={self.Re_tau:.1f}, "
                f"case={self.case_dir})")

    def plot(self, last=None, out_dir="plots", mkm_cache_dir=None):
        """Average (optionally over just the last `last` snapshots), fold
        by channel-centre symmetry, and save the 9 figures described in the
        module docstring under `<case_dir>/<out_dir>` (default 'plots')."""
        out = self.case_dir / out_dir
        out.mkdir(exist_ok=True)
        mkm_cache_dir = Path(mkm_cache_dir) if mkm_cache_dir else self.stats_dir / "mkm_cache"

        Umean, Rij, budget_data = self.Umean, self.Rij, self.budget_data
        if last is not None:
            n_use = min(int(last), self.nsamples)
            if n_use < 1:
                raise ValueError(f"last must be >= 1, got {last}")
            Umean = Umean[:, :, -n_use:]
            Rij = Rij[:, :, -n_use:]
            budget_data = {k: v[:, :, -n_use:] for k, v in budget_data.items()}
            print(f"  Using last {n_use} of {self.nsamples} snapshots")

        out_ny = self.Umean.shape[1]
        n_half = out_ny // 2
        Umean = symmetry_fold(Umean, SYM_UMEAN)
        Rij = symmetry_fold(Rij, SYM_RIJ)
        budget_data = {k: symmetry_fold(v, SYM_RIJ) for k, v in budget_data.items()}
        yp = self.yp[:n_half]
        used = last if last is not None else self.nsamples
        print(f"  Symmetry folded: {out_ny} wall-normal points -> {n_half} half-channel, "
              f"nsamples {used} -> {2 * used} effective samples")

        u_tau, nu, Re_tau = self.u_tau, self.nu, self.Re_tau
        rij_scale = u_tau ** 2
        budget_scale = u_tau ** 4 / nu

        Um_mean, Um_std = _mean_std(Umean)
        Rij_mean, Rij_std = _mean_std(Rij)
        term_stats = {k: _mean_std(v) for k, v in budget_data.items()}
        term_color = {stem: col for stem, _, col in BUDGET_TERMS}
        term_label = {stem: lbl for stem, lbl, _ in BUDGET_TERMS}

        print("Fetching MKM chan590 DNS reference data ...")
        ref = load_mkm_reference(mkm_cache_dir)
        print("  Done.")

        # ── Figure 1: mean streamwise velocity ────────────────────────────
        fig, ax = plt.subplots(figsize=(6, 5))
        ax.plot(ref["yp_means"], ref["Umean"], color="C0", lw=1.5,
                label=r"MKM DNS $\langle U \rangle^+$")
        Up = Um_mean[0] / u_tau
        Up_err = Um_std[0] / u_tau
        ax.errorbar(yp, Up, yerr=Up_err,
                    fmt="o", ms=4, mfc="none", mec="k", mew=0.8,
                    ecolor="k", elinewidth=0.6, capsize=2, ls="none",
                    label=r"LES $\langle U \rangle^+$")
        kappa, B = 0.41, 5.2
        yp_log = np.array([30.0, yp.max()])
        ax.plot(yp_log, np.log(yp_log) / kappa + B,
                "k--", lw=1.2, label=rf"$(1/\kappa)\ln y^+ + {B}$")
        yp_vis = np.array([1.0, 12.0])
        ax.plot(yp_vis, yp_vis, "k:", lw=1.2, label=r"$U^+ = y^+$")
        _log_xaxis(ax, yp)
        ax.set_ylabel(r"$U^+$", fontsize=12)
        ax.set_title(rf"Mean velocity  ($Re_{{\tau}} \approx {Re_tau:.0f}$)", fontsize=12)
        ax.legend(fontsize=9)
        fig.tight_layout()
        fig.savefig(out / "mean_velocity.png", dpi=150)
        plt.close(fig)
        print("Saved  plots/mean_velocity.png")

        # ── Figure 2: Reynolds stress profiles ────────────────────────────
        fig, ax = plt.subplots(figsize=(7, 5))
        prop_cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]
        ref_stress = [ref["R_uu"], ref["R_vv"], ref["R_ww"],
                      ref["R_uv"], ref["R_uw"], ref["R_vw"]]
        for ic, (cname, clabel) in enumerate(zip(COMP_NAMES, COMP_LABELS)):
            col = prop_cycle[ic % len(prop_cycle)]
            ax.plot(ref["yp_stress"], ref_stress[ic], color=col, lw=1.5,
                    label=clabel + " (MKM DNS)")
            les_label = "LES (present)" if ic == 0 else "_"
            ax.errorbar(yp, Rij_mean[ic] / rij_scale,
                        yerr=Rij_std[ic] / rij_scale,
                        fmt="o", ms=4, mfc="none", mec="k", mew=0.8,
                        ecolor="k", elinewidth=0.6, capsize=2, ls="none",
                        label=les_label)
        ax.axhline(0, color="k", lw=0.6, ls="--")
        _log_xaxis(ax, yp)
        ax.set_ylabel(r"$\langle u_i'u_j' \rangle^+$", fontsize=12)
        ax.set_title(rf"Reynolds stresses  ($Re_{{\tau}} \approx {Re_tau:.0f}$)", fontsize=12)
        ax.legend(fontsize=9, ncol=2)
        fig.tight_layout()
        fig.savefig(out / "reynolds_stresses.png", dpi=150)
        plt.close(fig)
        print("Saved  plots/reynolds_stresses.png")

        # ── Figures 3-8: stress budgets, one per Rij component ────────────
        for ic, (cname, clabel) in enumerate(zip(COMP_NAMES, COMP_LABELS)):
            fig, ax = plt.subplots(figsize=(7, 5))
            for stem, _, _ in BUDGET_TERMS:
                mn, sd = term_stats[stem]
                mn_w, sd_w = mn[ic] / budget_scale, sd[ic] / budget_scale
                _shade(ax, yp, mn_w, sd_w, color=term_color[stem], lw=1.5,
                       label=term_label[stem] + "  (mean +- 1 sigma)")
            ax.axhline(0, color="k", lw=0.6, ls="--")
            _log_xaxis(ax, yp)
            ax.set_ylabel(r"Budget term $\times\;\nu/u_\tau^4$", fontsize=12)
            ax.set_title(rf"RSB {clabel}  ($Re_{{\tau}} \approx {Re_tau:.0f}$)", fontsize=12)
            ax.legend(fontsize=8, loc="best", ncol=1)
            fig.tight_layout()
            fname = out / f"budget_{cname}.png"
            fig.savefig(fname, dpi=150)
            plt.close(fig)
            print(f"Saved  plots/budget_{cname}.png")

        # ── Figure 9: budget overview panel (2x3) ──────────────────────────
        fig, axes = plt.subplots(2, 3, figsize=(14, 8), sharex=True)
        axes = axes.ravel()
        for ic, (cname, clabel) in enumerate(zip(COMP_NAMES, COMP_LABELS)):
            ax = axes[ic]
            for stem, _, _ in BUDGET_TERMS:
                mn, sd = term_stats[stem]
                mn_w, sd_w = mn[ic] / budget_scale, sd[ic] / budget_scale
                _shade(ax, yp, mn_w, sd_w, color=term_color[stem], lw=1.3,
                       label=term_label[stem])
            ax.axhline(0, color="k", lw=0.5, ls="--")
            ax.set_xscale("log")
            ax.set_xlim(yp.min() * 0.9, yp.max() * 1.05)
            ax.set_title(clabel, fontsize=11)
            if ic >= 3:
                ax.set_xlabel(r"$y^+$", fontsize=10)
            if ic % 3 == 0:
                ax.set_ylabel(r"Budget $\times\;\nu/u_\tau^4$", fontsize=9)
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="lower center", ncol=4,
                   fontsize=8.5, bbox_to_anchor=(0.5, -0.04))
        fig.suptitle(
            rf"Reynolds stress budgets  ($Re_{{\tau}} \approx {Re_tau:.0f}$)  --  mean +- 1 sigma shading",
            fontsize=13, y=1.01)
        fig.tight_layout()
        fig.savefig(out / "budget_overview.png", dpi=150, bbox_inches="tight")
        plt.close(fig)
        print("Saved  plots/budget_overview.png")

        return out
