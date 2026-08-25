#!/usr/bin/env python3
"""
plot_snapshot.py — Wall-normal profile plots from fdm-dopamine snapshots.

Reads input_parameters from the working directory to determine active physics
(sgs_model, sediment_flag) and file layout, then loads snapshot(s) and plots
x-z averaged wall-normal profiles.  Figures are saved to postProcessing/.

Usage (run from the case directory):
    python postProcessing/plot_snapshot.py             # latest snapshot
    python postProcessing/plot_snapshot.py --step 500  # a specific step
    python postProcessing/plot_snapshot.py --all       # time-averaged over all snapshots
    python postProcessing/plot_snapshot.py --fields run2/fields --params run2/input_parameters
"""

import argparse
import sys
from pathlib import Path

import numpy as np

# ── locate snapshot_io next to this script ────────────────────────────────────
sys.path.insert(0, str(Path(__file__).parent))
from snapshot_io import parse_input_parameters, list_snapshots, read_snapshot

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except ImportError:
    sys.exit('matplotlib is required:  pip install matplotlib')


# ── CLI ───────────────────────────────────────────────────────────────────────

def _parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument('--params', default='input_parameters',
                   help='Fortran namelist file (default: input_parameters)')
    p.add_argument('--fields', default='fields',
                   help='directory containing binary snapshots (default: fields/)')
    p.add_argument('--step', type=int, default=None,
                   help='step number to load (default: latest available)')
    p.add_argument('--all', action='store_true',
                   help='average profiles over all available snapshots')
    p.add_argument('--out', default=None,
                   help='output figure path (default: postProcessing/profiles_<step>.png)')
    return p.parse_args()


# ── helpers ───────────────────────────────────────────────────────────────────

def _utau(dPdx, Ly):
    """Friction velocity from mean pressure gradient and channel half-height."""
    return np.sqrt(abs(dPdx) * 0.5 * Ly)


def _xz_mean(snap):
    """Return a dict of x-z averaged 1-D profiles for every 3-D field in snap."""
    return {k: v.mean(axis=(0, 2))
            for k, v in snap.items()
            if isinstance(v, np.ndarray) and v.ndim == 3}


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    args = _parse_args()

    # ── case parameters ───────────────────────────────────────────────────────
    prm          = parse_input_parameters(args.params)
    sgs_model    = int(prm.get('sgs_model',    0))
    sediment_flag = int(prm.get('sediment_flag', 0))
    fileout      = str(prm.get('fileout',       'channel_test'))
    nu           = float(prm.get('nu',    1e-6))
    dPdx         = float(prm.get('dpdx',  0.0))
    Ly           = float(prm.get('ly',    1.0))

    utau             = _utau(dPdx, Ly)
    use_wall_units   = utau > 0.0

    print(f'  sgs_model      = {sgs_model}')
    print(f'  sediment_flag  = {sediment_flag}')
    print(f'  fileout        = {fileout!r}')
    if use_wall_units:
        Retau = utau * 0.5 * Ly / nu
        print(f'  u_tau (est.)   = {utau:.4g}  (Re_tau ~ {Retau:.1f})')

    # ── discover snapshots ────────────────────────────────────────────────────
    snaps = list_snapshots(args.fields, fileout)
    if not snaps:
        sys.exit(f'No snapshots found in {args.fields!r} with prefix {fileout!r}')
    print(f'  found {len(snaps)} snapshot(s): '
          f'steps {snaps[0][0]} … {snaps[-1][0]}')

    read_kw = dict(sgs_model=sgs_model, sediment_flag=sediment_flag)

    # ── load snapshot(s) ──────────────────────────────────────────────────────
    if args.all:
        print('  loading all snapshots …')
        acc   = None
        count = 0
        ym    = None
        for step, fpath in snaps:
            s = read_snapshot(fpath, **read_kw)
            if ym is None:
                ym = s['ym']
            prof = _xz_mean(s)
            if acc is None:
                acc = {k: v.copy() for k, v in prof.items()}
            else:
                for k in acc:
                    acc[k] += prof[k]
            count += 1
        for k in acc:
            acc[k] /= count
        profile     = acc
        step_label  = f'time-average  ({count} snapshots)'
        step_tag    = 'avg'

    else:
        if args.step is not None:
            match = [fp for st, fp in snaps if st == args.step]
            if not match:
                sys.exit(f'Step {args.step} not found.')
            step, fpath = args.step, match[0]
        else:
            step, fpath = snaps[-1]
        print(f'  loading step {step}: {fpath}')
        s        = read_snapshot(fpath, **read_kw)
        ym       = s['ym']
        profile  = _xz_mean(s)
        step_label = f'step {step}'
        step_tag   = str(step)

    # ── figure layout ─────────────────────────────────────────────────────────
    panels = ['U']
    if sgs_model != 0:
        panels.append('nu_t')
    if sediment_flag >= 1:
        panels.append('C')

    ncols = len(panels)
    fig, axes = plt.subplots(1, ncols, figsize=(4.5 * ncols, 5.5), sharey=True)
    if ncols == 1:
        axes = [axes]

    fig.suptitle(f'{fileout}  —  {step_label}', fontsize=10, y=0.98)

    # y-axis: wall units or physical
    if use_wall_units:
        yplot  = ym * utau / nu
        ylabel = r'$y^+$'
    else:
        yplot  = ym
        ylabel = r'$y$  [m]'

    axes[0].set_ylabel(ylabel)

    # ── panel: mean streamwise velocity ───────────────────────────────────────
    ax = axes[panels.index('U')]
    U_mean = profile['U']
    if use_wall_units:
        ax.plot(U_mean / utau, yplot, color='C0')
        ax.set_xlabel(r'$\langle U \rangle^+$')
    else:
        ax.plot(U_mean, yplot, color='C0')
        ax.set_xlabel(r'$\langle U \rangle$  [m s$^{-1}$]')
    ax.set_title('Mean streamwise velocity')
    ax.grid(True, alpha=0.3)

    # ── panel: SGS turbulent viscosity ratio ──────────────────────────────────
    if sgs_model != 0:
        ax = axes[panels.index('nu_t')]
        ax.plot(profile['nu_t'] / nu, yplot, color='C1')
        ax.set_xlabel(r'$\langle \nu_t \rangle / \nu$')
        ax.set_title('SGS eddy-viscosity ratio')
        ax.grid(True, alpha=0.3)

    # ── panel: mean scalar concentration ──────────────────────────────────────
    if sediment_flag >= 1:
        ax = axes[panels.index('C')]
        ax.plot(profile['C'], yplot, color='C2')
        ax.set_xlabel(r'$\langle C \rangle$')
        ax.set_title('Mean scalar concentration')
        ax.grid(True, alpha=0.3)

    plt.tight_layout()

    # ── save ──────────────────────────────────────────────────────────────────
    out_path = args.out or str(Path(__file__).parent / f'profiles_{step_tag}.png')
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    print(f'\n  figure saved → {out_path}')


if __name__ == '__main__':
    main()
