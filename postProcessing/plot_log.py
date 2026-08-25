#!/usr/bin/env python3
"""
plot_log.py  –  Extract and plot solver diagnostics from run.log.

Each time-step in the log wraps across two physical lines; RSB sample
notifications and all header text are silently skipped.  Every token that
can be parsed as a float is accumulated; once N_COLS values per step are
identified the array is shaped into (n_steps, n_cols) and plotted.

Usage
-----
    python plot_log.py                          # plots Umean and Umax
    python plot_log.py -v Umean,Umax,divergence # plots listed variables
    python plot_log.py path/to.log              # custom log path
    python plot_log.py -v CFLc,dt path/to.log   # both options together

Available variable names
------------------------
    Umean       mean streamwise velocity  <U>
    Umax        maximum velocity          |U|max
    divergence  max divergence            |div|
    CFLc        convective CFL            CFL_c
    CFLv        viscous CFL               CFL_v
    dt          time-step size            dt
"""

import argparse
import re
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

LOG_FILE = Path("run.log")
OUT_DIR  = Path("plots")

# ── LaTeX labels for known column names ───────────────────────────────────────
LABEL_MAP: dict[str, str] = {
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

# ── friendly name → internal column name ─────────────────────────────────────
ALIAS_MAP: dict[str, str] = {
    "umean":      "<U>",
    "umax":       "|U|max",
    "divergence": "|div|",
    "cflc":       "CFL_c",
    "cflv":       "CFL_v",
    "dt":         "dt",
}

# Default set shown when no --vars argument is given
DEFAULT_VARS = ("<U>", "|U|max")


# ── parser ─────────────────────────────────────────────────────────────────────

def parse_log(path: Path):
    """
    Parse a dopamine run.log file.

    Strategy
    --------
    1. Locate the column-header line(s) that begin with "step".
    2. Collect all header tokens (the header may be wrapped) up to the
       first separator line of dashes.
    3. From after the separator onward, accept only lines where *every*
       whitespace-separated token is a valid float.  This transparently
       handles both the first and second physical halves of each wrapped
       step line while discarding RSB notifications, blank lines, and any
       other non-numeric text.
    4. Accumulate all accepted floats and reshape into (n_steps, n_cols).
    5. Also collect RSB sample step numbers for annotation.

    Returns
    -------
    col_names : list[str]
    data      : ndarray, shape (n_steps, n_cols)
    rsb_steps : list[int]   solver steps at which RSB samples were written
    """
    lines = path.read_text().splitlines()

    # ── 1. find the column-header line ──────────────────────────────────────
    header_idx = None
    for i, line in enumerate(lines):
        if re.match(r"\s+step\s+", line):
            header_idx = i
            break
    if header_idx is None:
        raise ValueError(f"Column header ('step …') not found in {path}")

    # ── 2. collect the (possibly wrapped) header tokens ─────────────────────
    #   The header ends at the first line of dashes.
    header_tokens: list[str] = []
    i = header_idx
    while i < len(lines) and not re.match(r"\s*-{3,}", lines[i]):
        header_tokens.extend(lines[i].split())
        i += 1

    col_names = header_tokens          # e.g. ['step','t','<U>','|U|max', …]

    # ── 3. skip all separator lines (dashes, possibly wrapped) ──────────────
    while i < len(lines) and re.match(r"\s*-{3,}", lines[i]):
        i += 1
    data_start = i

    # ── 4. accumulate numeric values and RSB step numbers ───────────────────
    values:    list[float] = []
    rsb_steps: list[int]   = []

    for line in lines[data_start:]:
        # RSB notification: extract the step number for later annotation
        m_rsb = re.match(r"\s*RSB:\s+wrote sample\s+\d+\s+at step\s+(\d+)", line)
        if m_rsb:
            rsb_steps.append(int(m_rsb.group(1)))
            continue

        tokens = line.split()
        if not tokens:
            continue

        # Accept the line only if every token is a valid float; this rejects
        # all lines containing alphabetic text (header, notifications, etc.)
        try:
            row_vals = [float(t) for t in tokens]
        except ValueError:
            continue

        values.extend(row_vals)

    # ── 5. reshape ───────────────────────────────────────────────────────────
    n_cols    = len(col_names)
    n_rows    = len(values) // n_cols
    remainder = len(values) %  n_cols
    if remainder:
        print(
            f"  Warning: {remainder} trailing value(s) discarded "
            f"(incomplete last step).",
            file=sys.stderr,
        )

    data = np.array(values[: n_rows * n_cols]).reshape(n_rows, n_cols)
    return col_names, data, rsb_steps


# ── plotter ────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Plot solver diagnostics from a dopamine run.log."
    )
    parser.add_argument(
        "log",
        nargs="?",
        default=str(LOG_FILE),
        metavar="LOG_FILE",
        help="Path to the log file (default: run.log)",
    )
    parser.add_argument(
        "-v", "--vars",
        default=None,
        metavar="VAR1,VAR2,...",
        help=(
            "Comma-separated list of variables to plot. "
            "Choices: Umean, Umax, divergence, CFLc, CFLv, dt. "
            "If omitted, plots Umean and Umax."
        ),
    )
    args = parser.parse_args()
    log_path = Path(args.log)
    OUT_DIR.mkdir(exist_ok=True)

    print(f"Parsing {log_path} …")
    col_names, data, rsb_steps = parse_log(log_path)
    n_steps, n_cols = data.shape
    print(f"  {n_steps} steps  ×  {n_cols} columns: {col_names}")
    print(f"  {len(rsb_steps)} RSB samples found at steps: {rsb_steps[:5]}"
          + (" …" if len(rsb_steps) > 5 else ""))

    idx = {name: j for j, name in enumerate(col_names)}

    t = data[:, idx["t"]]

    # Columns to plot — filtered by --vars if provided, else default set
    if args.vars is not None:
        requested = [tok.strip() for tok in args.vars.split(",") if tok.strip()]
        plot_cols = []
        for alias in requested:
            col = ALIAS_MAP.get(alias.lower())
            if col is None:
                print(
                    f"  Warning: unknown variable '{alias}' — ignored. "
                    f"Valid choices: {', '.join(ALIAS_MAP)}",
                    file=sys.stderr,
                )
            elif col not in col_names:
                print(
                    f"  Warning: column '{col}' (alias '{alias}') not found in log — ignored.",
                    file=sys.stderr,
                )
            else:
                plot_cols.append(col)
        if not plot_cols:
            print("  Error: no valid columns to plot. Exiting.", file=sys.stderr)
            sys.exit(1)
    else:
        plot_cols = [c for c in DEFAULT_VARS if c in col_names]

    # ── global style ──────────────────────────────────────────────────────────
    plt.rcParams.update(RC)

    # ── layout ────────────────────────────────────────────────────────────────
    n_plots   = len(plot_cols)
    ncols_fig = 2
    nrows_fig = (n_plots + ncols_fig - 1) // ncols_fig

    fig, axes = plt.subplots(
        nrows_fig, ncols_fig,
        figsize=(11, 3.0 * nrows_fig),
        sharex=True,
    )
    axes = np.asarray(axes).ravel()

    for ax, name in zip(axes, plot_cols):
        y = data[:, idx[name]]

        # Subtle horizontal dotted grid, drawn behind everything
        ax.yaxis.grid(True, color="0.88", linewidth=0.5, linestyle=":", zorder=0)
        ax.set_axisbelow(True)

        # Data line
        ax.plot(t, y, lw=1.0, color=DATA_COLOR)

        # LaTeX y-axis label (fall back to raw name if not in the map)
        ax.set_ylabel(LABEL_MAP.get(name, name))

        # Log scale for divergence (spans many orders of magnitude)
        if "div" in name.lower():
            ax.set_yscale("log")

    # X-axis label only on the bottom row
    for ax in axes[n_plots - ncols_fig: n_plots]:
        ax.set_xlabel(r"$t$")

    # Hide the unused panel when there is an odd number of plots
    for ax in axes[n_plots:]:
        ax.set_visible(False)

    fig.suptitle(
        f"Solver diagnostics \u2013 {log_path.name}  "
        f"({n_steps} steps,  $t = {t[-1]:.2f}$)",
        fontsize=11,
    )
    fig.tight_layout()

    out_path = OUT_DIR / "run_log_diagnostics.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved  {out_path}")


if __name__ == "__main__":
    main()
