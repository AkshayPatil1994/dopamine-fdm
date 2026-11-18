"""Small helper of the HIT plotting scripts: per-run line styles on top of plotstyle.py (kept out of plotstyle.py itself)."""
import plotstyle as ps


def run_style(label, k):
    """dict(label, color, ls, marker) of the k-th run, colours from plotstyle.FALLBACK."""
    return dict(label=label, color=ps.FALLBACK[k % len(ps.FALLBACK)], ls="-", marker="o")


MARKERS = ["o", "s", "^", "D", "v", "P"]


def data_style(i, label, reference=False):
    """Numerical results are open symbols (no connecting line); the reference run is a solid black line."""
    if reference:
        return dict(color="k", ls="-", lw=1.4)
    return dict(color=ps.FALLBACK[i % len(ps.FALLBACK)], marker=MARKERS[i % len(MARKERS)], ls="none", ms=4.5)


def log_marks(n, m=22):
    """About m marker indices spread logarithmically over 1..n-1 (spectra live on log-k axes)."""
    import numpy as np
    return sorted({0} | set(int(round(x)) - 1 for x in np.logspace(0, np.log10(max(n, 2)), m))) if n > 3 else None
