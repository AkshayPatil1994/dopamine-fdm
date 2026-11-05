# SPDX-FileCopyrightText: 2026 DOPAMINE contributors
# SPDX-License-Identifier: AGPL-3.0-only
"""Shared manuscript figure style (serif text, Computer Modern maths, ticks on all sides pointing inwards).

Uses matplotlib's built-in mathtext, so no LaTeX installation is needed. Set DOPAMINE_USETEX=1 to
typeset with LaTeX instead.
"""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

BLUE = "#0072B2"
ORANGE = "#D55E00"
GREY = "#7f7f7f"

FALLBACK = ["#009E73", "#CC79A7", "#E69F00", "#56B4E9", BLUE, ORANGE]


def apply():
    plt.rcParams.update({
        "text.usetex": os.environ.get("DOPAMINE_USETEX") == "1",
        "font.family": "serif",
        "font.serif": ["STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "cm",
        "axes.formatter.use_mathtext": True,
        "font.size": 10,
        "axes.labelsize": 11,
        "axes.linewidth": 0.8,
        "legend.fontsize": 8.5,
        "legend.frameon": False,
        "legend.handlelength": 2.6,
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True,
        "xtick.minor.visible": True, "ytick.minor.visible": True,
        "xtick.major.size": 4, "ytick.major.size": 4,
        "xtick.minor.size": 2, "ytick.minor.size": 2,
        "xtick.major.width": 0.8, "ytick.major.width": 0.8,
        "xtick.minor.width": 0.6, "ytick.minor.width": 0.6,
        "lines.linewidth": 1.6,
        "lines.markersize": 4.5,
        "lines.markerfacecolor": "white",
        "lines.markeredgewidth": 1.1,
        "savefig.dpi": 400,
        "savefig.bbox": "tight",
        "pdf.fonttype": 42,
    })


def vorticity_cmap():
    """Diverging map for signed vorticity: cool blues for negative, the magma_r cream-to-purple ramp for positive."""
    from matplotlib.colors import LinearSegmentedColormap

    colors = ["#16224f", "#2c68a3", "#7fb3d5", "#d9e8ef", "#fbf3e4", "#fdd9a5", "#f68f5c", "#d1426d", "#7a1f7d",
              "#2a0a4d"]
    return LinearSegmentedColormap.from_list("vorticity", colors, N=256)


def panel_label(ax, text):
    ax.text(0.0, 1.015, text, transform=ax.transAxes, ha="left", va="bottom")


def save(fig, path):
    """Save to `path`; a PDF request also writes a PNG preview next to it."""
    fig.savefig(path)
    out = [path]
    if str(path).lower().endswith(".pdf"):
        png = str(path)[:-4] + ".png"
        fig.savefig(png)
        out.append(png)
    return out
