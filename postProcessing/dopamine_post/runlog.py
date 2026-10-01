"""
runlog.py -- solver diagnostics from a dopamine run.log.

Each time-step in the log wraps across two physical lines; RSB sample
notifications and all header text are silently skipped. Every token that
can be parsed as a float is accumulated; once N_COLS values per step are
identified (N_COLS taken from the log's own column-header line, so this
adapts to whatever columns a given run actually wrote) the array is shaped
into (n_steps, n_cols).

This wrapped-line, header-driven parse is richer than `_core.times_from_log`
(which assumes one fixed-width monitor line per step and exists only to look
up physical time for a handful of solver steps elsewhere, e.g.
`dopamine_post.particles`) -- so this module keeps its own full parser
rather than force-fitting that narrower helper.

Library use
-----------
    from dopamine_post.runlog import RunLog

    log = RunLog.read("run.log")
    log.plot()                                   # Umean, Umax (default)
    log.plot(variables=["Umean", "Umax", "divergence"])

Available variable names
-------------------------
    Umean       mean streamwise velocity  <U>
    Umax        maximum velocity          |U|max
    divergence  max divergence            |div|
    CFLc        convective CFL            CFL_c
    CFLv        viscous CFL               CFL_v
    dt          time-step size            dt
"""
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── LaTeX labels for known column names ───────────────────────────────────────
LABEL_MAP = {
    "<U>":    r"$\langle U \rangle$",
    "|U|max": r"$|U|_{\mathrm{max}}$",
    "|div|":  r"$|\nabla \cdot \mathbf{u}|_{\mathrm{max}}$",
    "CFL_c":  r"$\mathrm{CFL}_c$",
    "CFL_v":  r"$\mathrm{CFL}_{\nu}$",
    "dt":     r"$\Delta t$",
    "t":      r"$t$",
    "step":   r"step",
}

# ── plot style ────────────────────────────────────────────────────────────────
RC = {
    "text.usetex":       False,
    "mathtext.fontset":  "stix",
    "font.family":       "serif",
    "font.size":         10,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "axes.linewidth":    0.8,
    "axes.labelsize":    11,
    "xtick.direction":   "out",
    "ytick.direction":   "out",
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.major.size":  4,
    "ytick.major.size":  4,
    "xtick.labelsize":   9,
    "ytick.labelsize":   9,
}

DATA_COLOR = "#2C5F8A"   # deep steel blue

# ── friendly name -> internal column name ─────────────────────────────────────
ALIAS_MAP = {
    "umean":      "<U>",
    "umax":       "|U|max",
    "divergence": "|div|",
    "cflc":       "CFL_c",
    "cflv":       "CFL_v",
    "dt":         "dt",
}

# Default set shown when no variables are given
DEFAULT_VARS = ("<U>", "|U|max")


def parse_log(path):
    """Parse a dopamine run.log file.

    Strategy
    --------
    1. Locate the column-header line(s) that begin with "step".
    2. Collect all header tokens (the header may be wrapped) up to the
       first separator line of dashes.
    3. From after the separator onward, accept only lines where *every*
       whitespace-separated token is a valid float. This transparently
       handles both the first and second physical halves of each wrapped
       step line while discarding RSB notifications, blank lines, and any
       other non-numeric text.
    4. Accumulate all accepted floats and reshape into (n_steps, n_cols).
    5. Also collect RSB sample step numbers for annotation.

    Returns (col_names, data, rsb_steps).
    """
    path = Path(path)
    lines = path.read_text().splitlines()

    header_idx = None
    for i, line in enumerate(lines):
        if re.match(r"\s+step\s+", line):
            header_idx = i
            break
    if header_idx is None:
        raise ValueError(f"Column header ('step ...') not found in {path}")

    header_tokens = []
    i = header_idx
    while i < len(lines) and not re.match(r"\s*-{3,}", lines[i]):
        header_tokens.extend(lines[i].split())
        i += 1
    col_names = header_tokens

    while i < len(lines) and re.match(r"\s*-{3,}", lines[i]):
        i += 1
    data_start = i

    values = []
    rsb_steps = []
    for line in lines[data_start:]:
        m_rsb = re.match(r"\s*RSB:\s+wrote sample\s+\d+\s+at step\s+(\d+)", line)
        if m_rsb:
            rsb_steps.append(int(m_rsb.group(1)))
            continue

        tokens = line.split()
        if not tokens:
            continue
        try:
            row_vals = [float(t) for t in tokens]
        except ValueError:
            continue
        values.extend(row_vals)

    n_cols = len(col_names)
    n_rows = len(values) // n_cols
    remainder = len(values) % n_cols
    if remainder:
        print(f"  Warning: {remainder} trailing value(s) discarded (incomplete last step).",
              file=sys.stderr)

    data = np.array(values[:n_rows * n_cols]).reshape(n_rows, n_cols)
    return col_names, data, rsb_steps


class RunLog:
    """A parsed run.log: `col_names` (as written by the solver, e.g. 'step',
    't', '<U>', '|U|max', ...), `data` (n_steps, n_cols), and `rsb_steps`
    (solver steps at which an RSB sample was written, for annotation).
    """

    def __init__(self, path, col_names, data, rsb_steps):
        self.path = Path(path)
        self.col_names = col_names
        self.data = data
        self.rsb_steps = rsb_steps
        self._idx = {name: j for j, name in enumerate(col_names)}

    @classmethod
    def read(cls, path="run.log"):
        col_names, data, rsb_steps = parse_log(path)
        return cls(path, col_names, data, rsb_steps)

    def __repr__(self):
        n_steps, n_cols = self.data.shape
        return f"RunLog({n_steps} steps x {n_cols} cols, path={self.path})"

    def get(self, name):
        """Column values (n_steps,) by raw column name (e.g. '<U>') or a
        friendly ALIAS_MAP alias (case-insensitive, e.g. 'Umean')."""
        col = ALIAS_MAP.get(name.lower(), name)
        if col not in self._idx:
            raise KeyError(f"unknown column '{name}' (have {self.col_names})")
        return self.data[:, self._idx[col]]

    @property
    def t(self):
        return self.get("t")

    def plot(self, variables=None, out=None):
        """Plot the requested variables (friendly ALIAS_MAP names, e.g.
        ["Umean", "Umax", "divergence"], or a comma-separated string of the
        same; default: Umean and Umax) vs. time, one subplot each in a
        2-column grid. Saves to `out` (default '<run.log's directory>/
        plots/run_log_diagnostics.png') and returns the Figure.
        """
        n_steps, n_cols = self.data.shape
        idx = self._idx
        t = self.data[:, idx["t"]]

        if variables is not None:
            requested = (variables.split(",") if isinstance(variables, str) else variables)
            requested = [tok.strip() for tok in requested if tok.strip()]
            plot_cols = []
            for alias in requested:
                col = ALIAS_MAP.get(alias.lower())
                if col is None:
                    print(f"  Warning: unknown variable '{alias}' -- ignored. "
                          f"Valid choices: {', '.join(ALIAS_MAP)}", file=sys.stderr)
                elif col not in self.col_names:
                    print(f"  Warning: column '{col}' (alias '{alias}') not found in log -- ignored.",
                          file=sys.stderr)
                else:
                    plot_cols.append(col)
            if not plot_cols:
                raise ValueError("no valid columns to plot")
        else:
            plot_cols = [c for c in DEFAULT_VARS if c in self.col_names]

        plt.rcParams.update(RC)

        n_plots = len(plot_cols)
        ncols_fig = 2
        nrows_fig = (n_plots + ncols_fig - 1) // ncols_fig

        fig, axes = plt.subplots(nrows_fig, ncols_fig, figsize=(11, 3.0 * nrows_fig), sharex=True)
        axes = np.asarray(axes).ravel()

        for ax, name in zip(axes, plot_cols):
            y = self.data[:, idx[name]]
            ax.yaxis.grid(True, color="0.88", linewidth=0.5, linestyle=":", zorder=0)
            ax.set_axisbelow(True)
            ax.plot(t, y, lw=1.0, color=DATA_COLOR)
            ax.set_ylabel(LABEL_MAP.get(name, name))
            if "div" in name.lower():
                ax.set_yscale("log")

        for ax in axes[n_plots - ncols_fig: n_plots]:
            ax.set_xlabel(r"$t$")
        for ax in axes[n_plots:]:
            ax.set_visible(False)

        fig.suptitle(f"Solver diagnostics – {self.path.name}  "
                     f"({n_steps} steps,  $t = {t[-1]:.2f}$)", fontsize=11)
        fig.tight_layout()

        if out is None:
            out_dir = self.path.parent / "plots"
            out_dir.mkdir(exist_ok=True)
            out = out_dir / "run_log_diagnostics.png"
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"Saved  {out}")
        return fig
