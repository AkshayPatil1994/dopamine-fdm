"""
diagnostics.py -- the solver's plain-text time series, read into {column name: 1-D array} dicts.

    from dopamine_post import diagnostics as dg

    d = dg.read_vof_diag("vof_diag.dat")          # two-fluid runs: d["t"], d["Vliq"], d["KE_liquid"], ...
    d = dg.read_wall_diag("wall_diag.dat")        # flat wall model: d["u_tau_lo"], d["u_tau_hi"], ...
    d = dg.read_ibm_forces("ibm_forces.csv")      # IBM runs: d["Fx_pres"], d["Fx_visc"], ...
    d = dg.read_particle_counts("fields/particle_count.dat")

hit_stats.dat is read by `dopamine_post.hit.read_stats`. Each file is appended to by the solver and may hold rows of an earlier, restarted run; `last_wins=True` (the default)
keeps the last row of each repeated step/time so a restarted history is monotonic.
"""
import re
from pathlib import Path

import numpy as np

_VOF_DIAG = ("step", "t", "dt", "Vliq", "dVliq_rel", "Cmin", "Cmax", "n_interface", "clip_loss_cum", "Co_adv_max", "xc", "yc",
             "zc", "int_C1mC", "cos_moment", "max_U", "max_V", "max_W", "PCG_its_step", "PCG_res", "n_sub", "bnd_cum",
             "relax_cum", "x_of_max_U", "y_of_max_U", "KE_liquid", "KE_gas", "L1_C_minus_C0", "p_liquid_minus_p_gas",
             "HF_fallback_cells", "sum_p_dV", "sum_p2_dV")
_WALL_DIAG = ("step", "t", "u_tau_lo", "u_tau_hi", "qT_lo", "qT_hi", "invL_lo", "invL_hi")


def _dedupe(d, key, last_wins):
    """Rows of a restarted run repeat earlier keys: keep the last occurrence of each, in order."""
    if not last_wins or d[key].size == 0:
        return d
    keep = np.sort(len(d[key]) - 1 - np.unique(d[key][::-1], return_index=True)[1])
    return {k: v[keep] for k, v in d.items()}


def _table(path, names, key, last_wins):
    a = np.atleast_2d(np.loadtxt(path, comments="#"))
    if a.size == 0:
        a = np.empty((0, len(names)))
    if a.shape[1] != len(names):
        raise ValueError(f"{path}: {a.shape[1]} columns, expected {len(names)} ({', '.join(names)})")
    return _dedupe({n: a[:, i] for i, n in enumerate(names)}, key, last_wins)


def read_vof_diag(path="vof_diag.dat", last_wins=True):
    """vof_diag.dat (every two-fluid run): 32 columns, see _VOF_DIAG / the file's own header."""
    return _table(path, _VOF_DIAG, "step", last_wins)


def read_wall_diag(path="wall_diag.dat", last_wins=True):
    """wall_diag.dat (flat y-wall model): wall-averaged friction velocity, heat flux, inverse Obukhov length (low/high wall)."""
    return _table(path, _WALL_DIAG, "step", last_wins)


def read_vof_gauges(path="vof_gauges.dat", last_wins=True):
    """vof_gauges.dat (wave_gauge_x): {'t': ..., 'x': gauge positions, 'eta': (nt, ngauge) free-surface elevations}."""
    first = Path(path).read_text().splitlines()[0]
    x = np.array([float(v) for v in first.split("=", 1)[1].split()])
    a = np.atleast_2d(np.loadtxt(path, comments="#"))
    d = _dedupe({"t": a[:, 0], "i": np.arange(len(a))}, "t", last_wins)
    return {"t": d["t"], "x": x, "eta": a[d["i"], 1:]}


def read_ibm_forces(path="ibm_forces.csv", last_wins=True):
    """ibm_forces.csv: step, t, total (ibm), pressure and viscous force components (Fx_ibm ... Fz_visc)."""
    names = open(path).readline().strip().split(",")
    a = np.atleast_2d(np.loadtxt(path, delimiter=",", skiprows=1))
    return _dedupe({n: a[:, i] for i, n in enumerate(names)}, "step", last_wins)


def read_particle_counts(path="fields/particle_count.dat", last_wins=True):
    """particle_count.dat: istep, t and the counts exited / deposited / reinjected / active (since the previous write, active = now)."""
    names = open(path).readline().strip().split(",")
    a = np.atleast_2d(np.loadtxt(path, delimiter=",", skiprows=1))
    return _dedupe({n: a[:, i] for i, n in enumerate(names)}, "istep", last_wins)


def read_particle_deposit(path):
    """particle_deposit_file: {'istep', 't', 'count': (nt, nxg) deposits per x cell since the previous write}."""
    a = np.atleast_2d(np.loadtxt(path, delimiter=",", skiprows=1))
    return {"istep": a[:, 0], "t": a[:, 1], "count": a[:, 2:]}


def read_wave_series_check(path="fields/wave_series_check.out"):
    """wave_series_check.out (wave_type 4): the commented header as a dict of its numbers plus the columns
    t, eta_input_minus_mean, eta_reconstructed, error."""
    head = {}
    for ln in Path(path).read_text().splitlines():
        if ln.startswith("#"):
            for k, v in re.findall(r"([A-Za-z][\w ()/]*?)\s+(-?\d+\.?\d*(?:[Ee][-+]?\d+)?)", ln[1:]):
                head[k.strip()] = float(v)
    a = np.loadtxt(path, comments="#")
    cols = ("t", "eta_input_minus_mean", "eta_reconstructed", "error")
    return {"header": head, **{c: a[:, i] for i, c in enumerate(cols)}}


def plot_vof_diag(d, out=None):
    """Volume error, interface count, kinetic energies and PCG iterations against time."""
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 2, figsize=(10, 6), sharex=True)
    for a, (k, lab) in zip(ax.ravel(), (("dVliq_rel", "(Vliq-V0)/V0"), ("n_interface", "interface cells"),
                                        ("KE_liquid", "KE liquid"), ("PCG_its_step", "PCG its / step"))):
        a.plot(d["t"], d[k]); a.set_ylabel(lab)
        if k == "KE_liquid":
            a.plot(d["t"], d["KE_gas"], label="gas"); a.legend(["liquid", "gas"])
    for a in ax[1]:
        a.set_xlabel("t")
    fig.tight_layout()
    if out:
        fig.savefig(out, dpi=150)
    return fig
