#!/usr/bin/env python3
"""
plot_sdf.py — Load and plot two slices of the cell-centre SDF (SDF_in).

Slices produced:
  1. X-Z plane at y = Y_SLICE  (horizontal cut through the rough-wall layer)
  2. X-Y plane at z = Z_SLICE  (streamwise–wall-normal cross-section)

Binary layout of SDF_in (written by fdm-dopamine when ibm_input_mode=1):
  big-endian float64, Fortran column-major order
  shape  (nxg, nyg, nzm)  =  (nxm+2, nym+2, nzm)
  ghost layers present in x and y; no ghost layer in z
  → strip to (nxm, nym, nzm) by removing the first/last index in x and y

Usage (run from the case directory):
    python plot_sdf.py
    python plot_sdf.py --sdf SDF_in --params input_parameters
    python plot_sdf.py --y-slice 0.1 --z-slice 2.4
"""

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from snapshot_io import parse_input_parameters, coords_from_grid

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable
except ImportError:
    sys.exit('matplotlib is required: pip install matplotlib')

# ── default slice positions ───────────────────────────────────────────────────
Y_SLICE = 0.05   # X-Z plane: y position
Z_SLICE = 3.0    # X-Y plane: z position

# ── helpers ───────────────────────────────────────────────────────────────────

def _nearest_index(coords, value):
    """Return the index in *coords* closest to *value*."""
    return int(np.argmin(np.abs(coords - value)))


SENTINEL_FRAC = 0.5   # matches GenSDF's flood_fill_mod.f90 SENTINEL_FRAC

def _pcolor(ax, H, V, data, title, xlabel, ylabel, cmap='RdBu_r', symmetric=True):
    # GenSDF leaves cells outside its geometry-focused AABB at the raw background
    # sentinel (scalarvalue, e.g. 1e10 — "far fluid, never computed"); including
    # those in the color scale crushes the real near-wall SDF into invisibility.
    scalarvalue = np.abs(data).max()
    is_sentinel = np.abs(data) >= SENTINEL_FRAC * scalarvalue
    finite = data[~is_sentinel]
    masked = np.ma.masked_where(is_sentinel, data)

    if finite.size:
        vmax = np.abs(finite).max() or 1.0
    else:
        vmax = scalarvalue or 1.0
    vmin = -vmax if symmetric else (finite.min() if finite.size else data.min())

    ax.set_facecolor('lightgray')   # shows through masked (sentinel/uncomputed) cells
    pcm = ax.pcolormesh(H, V, masked.T, shading='auto', cmap=cmap,
                        vmin=vmin, vmax=vmax)
    ax.set_title(title, fontsize=10)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    div = make_axes_locatable(ax)
    cax = div.append_axes('right', size='3%', pad=0.06)
    ax.get_figure().colorbar(pcm, cax=cax, label='φ  [m]')
    return pcm


# ── load SDF ─────────────────────────────────────────────────────────────────

def load_sdf(sdf_path, params_path):
    """
    Return (phi, xm, ym, zm) where phi has shape (nxm, nym, nzm).

    File layout (Fortran column-major, big-endian float64):
      padded(nxm+2, nym+2, nzm)
      axis 0 = x-streamwise  (ghost at indices 0 and nxm+1)
      axis 1 = y-wall-normal (ghost at indices 0 and nym+1)
      axis 2 = z-spanwise    (no ghost)

    Grid coordinates are read from fields/grid.out + fields/geometry.out
    (written by genGridandIC.f90) so the same non-uniform y-grid used by
    the solver is reproduced exactly, without needing a solver snapshot.
    """
    p = parse_input_parameters(params_path)
    nxm = p['nx'] - 1
    nym = p['ny'] - 1
    nzm = p['nz'] - 1
    nxg = nxm + 2   # with ghost layers
    nyg = nym + 2

    # Raw array: (nxg, nyg, nzm) big-endian float64, Fortran column-major.
    # order='F' maps flat memory directly to (nxg, nyg, nzm) Fortran indices,
    # so phi_full[ix, iy, iz] = SDF at solver cell (ix, iy, iz).
    raw = np.fromfile(sdf_path, dtype='>f8')
    expected = nxg * nyg * nzm
    if raw.size != expected:
        sys.exit(
            f'SDF_in size mismatch: got {raw.size} elements, '
            f'expected {expected} ({nxg}×{nyg}×{nzm})'
        )
    phi_full = raw.reshape((nxg, nyg, nzm), order='F')

    # Strip ghost layers in x and y
    phi = phi_full[1:-1, 1:-1, :]   # (nxm, nym, nzm)

    grid_path = Path(params_path).parent / 'fields' / 'grid.out'
    if not grid_path.is_file():
        sys.exit(f'{grid_path} not found — needed for grid coordinates.')
    xm, ym, zm = coords_from_grid(grid_path)

    return phi, xm, ym, zm


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--sdf',     default='SDF_in',         help='SDF binary file')
    ap.add_argument('--params',  default='input_parameters', help='namelist file')
    ap.add_argument('--y-slice', type=float, default=Y_SLICE,
                    help=f'y position for the X-Z slice (default: {Y_SLICE})')
    ap.add_argument('--z-slice', type=float, default=Z_SLICE,
                    help=f'z position for the X-Y slice (default: {Z_SLICE})')
    ap.add_argument('--out-dir', default='plots', help='output directory')
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f'Loading SDF from {args.sdf} …')
    phi, xm, ym, zm = load_sdf(args.sdf, args.params)
    print(f'  phi  shape: {phi.shape}   range: [{phi.min():.4f}, {phi.max():.4f}]')

    # ── slice 1: X-Z at y ≈ y_slice ──────────────────────────────────────────
    iy = _nearest_index(ym, args.y_slice)
    y_actual = ym[iy]
    print(f'\nX-Z slice: requested y={args.y_slice}, using iy={iy}, y={y_actual:.5f}')

    fig1, ax1 = plt.subplots(figsize=(12, 5))
    _pcolor(ax1, xm, zm, phi[:, iy, :],
            title=f'SDF — X-Z plane  (y = {y_actual:.4f},  iy = {iy})',
            xlabel='x', ylabel='z')
    ax1.contour(xm, zm, phi[:, iy, :].T, levels=[0.0], colors='k', linewidths=0.8)
    fig1.tight_layout()
    path1 = out_dir / f'sdf_xz_y{y_actual:.4f}.png'
    fig1.savefig(path1, dpi=150)
    print(f'  Saved → {path1}')
    plt.close(fig1)

    # ── slice 2: X-Y at z ≈ z_slice ──────────────────────────────────────────
    iz = _nearest_index(zm, args.z_slice)
    z_actual = zm[iz]
    print(f'\nX-Y slice: requested z={args.z_slice}, using iz={iz}, z={z_actual:.5f}')

    fig2, ax2 = plt.subplots(figsize=(12, 4))
    _pcolor(ax2, xm, ym, phi[:, :, iz],
            title=f'SDF — X-Y plane  (z = {z_actual:.4f},  iz = {iz})',
            xlabel='x', ylabel='y  (wall-normal)')
    ax2.contour(xm, ym, phi[:, :, iz].T, levels=[0.0], colors='k', linewidths=0.8)
    fig2.tight_layout()
    path2 = out_dir / f'sdf_xy_z{z_actual:.4f}.png'
    fig2.savefig(path2, dpi=150)
    print(f'  Saved → {path2}')
    plt.close(fig2)

    print('\nDone.')


if __name__ == '__main__':
    main()
