#!/usr/bin/env python3
"""
load_line_probes.py  —  Load fdm-dopamine 1-D line probe data.

The probe_output.f90 module writes a  <base>.bin  file (big-endian float64,
Fortran stream, no record markers) and a companion  <base>_meta.txt.

Binary layout
-------------
Each snapshot is a Fortran array  out1d(ncomp, npts)  written in column-major
order, which in Python (C row-major) looks like  (npts, ncomp).
All nsnaps snapshots are appended consecutively, giving a flat array of
nsnaps × ncomp × npts  float64 values.

API
---
  probe = load_probe('mysim_line_meta.txt')   # or the .bin path
  # probe['data'] has shape (nsnaps, ncomp, npts)
  # probe['comps'] e.g. ['U', 'V', 'W']
  u = probe['data'][:, probe['comps'].index('U'), :]   # (nsnaps, npts)

CLI
---
  python3 load_line_probes.py  <meta.txt>  [<meta2.txt> ...]  [--snap S]

  --snap S   Print raw values for snapshot index S (0-based).
"""

import argparse
import sys
import numpy as np
from pathlib import Path


# ── meta parser ───────────────────────────────────────────────────────────────

def read_meta(path):
    """Parse a <base>_meta.txt line-probe metadata file.  Returns a dict."""
    meta = {}
    for line in Path(path).read_text().splitlines():
        if '=' not in line:
            continue
        k, _, v = line.partition('=')
        k, v = k.strip(), v.strip()
        if k in ('ncomp', 'npts', 'nsnaps'):
            meta[k] = int(v)
        else:
            meta[k] = v
    if 'npts' not in meta:
        raise ValueError(f'{path}: not a line-probe meta file (no npts key)')
    return meta


# ── component helper ──────────────────────────────────────────────────────────

def comp_names(s):
    """Return ordered component name list from the comps string."""
    upper = s.strip().upper()
    names = [c for c in ('U', 'V', 'W', 'P') if c in upper]
    return names or ['U', 'V', 'W']


# ── public API ────────────────────────────────────────────────────────────────

def load_probe(meta_or_bin_path):
    """
    Load a 1-D line probe written by probe_output.f90.

    Accepts either the ``_meta.txt`` path or the ``.bin`` path; the companion
    file is located automatically.

    Parameters
    ----------
    meta_or_bin_path : str or Path

    Returns
    -------
    dict
        'data'   : ndarray, shape (nsnaps, ncomp, npts), float64, host-endian
        'comps'  : list[str]  e.g. ['U', 'V', 'W']
        'dir'    : str        probe direction, e.g. 'y'
        'nsnaps' : int        number of snapshots in the file
        'npts'   : int        number of spatial points along the line
        'ncomp'  : int        number of field components

    Notes
    -----
    Individual components can be accessed conveniently as::

        u = probe['data'][:, probe['comps'].index('U'), :]   # (nsnaps, npts)
    """
    p = Path(meta_or_bin_path)

    if p.suffix == '.bin':
        bin_path  = p
        meta_path = p.parent / (p.stem + '_meta.txt')
    else:
        meta_path = p
        stem = p.stem
        base = stem[:-5] if stem.endswith('_meta') else stem
        bin_path = p.parent / (base + '.bin')

    meta   = read_meta(meta_path)
    nc     = meta['ncomp']
    npts   = meta['npts']
    nsnaps = meta['nsnaps']
    dirstr = meta.get('dir', '?').strip()
    comps  = comp_names(meta.get('comps', 'UVW'))

    # Verify nsnaps against actual file size
    if bin_path.is_file():
        file_size = bin_path.stat().st_size
        snap_bytes = nc * npts * 8
        actual = file_size // snap_bytes if snap_bytes else 0
        if actual != nsnaps:
            print(f'  NOTE: meta says {nsnaps} snaps, file has {actual}; using file count',
                  file=sys.stderr)
            nsnaps = actual
    else:
        raise FileNotFoundError(f'Binary file not found: {bin_path}')

    # Each snapshot: Fortran (nc, npts) col-major == C (npts, nc) row-major.
    # Read as big-endian float64, reshape to (nsnaps, npts, nc), transpose
    # to (nsnaps, nc, npts) for component-first indexing.
    raw  = np.fromfile(bin_path, dtype='>f8', count=nsnaps * nc * npts)
    data = raw.reshape(nsnaps, npts, nc).transpose(0, 2, 1).astype('f8', copy=False)

    return {
        'data':   data,
        'comps':  comps,
        'dir':    dirstr,
        'nsnaps': nsnaps,
        'npts':   npts,
        'ncomp':  nc,
    }


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(
        description='Load and inspect fdm-dopamine 1-D line probe files.')
    p.add_argument('files', nargs='+', metavar='FILE',
                   help='Path(s) to <base>_meta.txt or <base>.bin')
    p.add_argument('--snap', type=int, default=None, metavar='S',
                   help='Print all values for snapshot index S (0-based)')
    p.add_argument('--plot', action='store_true',
                   help='Show a quick matplotlib plot of the last snapshot')
    args = p.parse_args()

    for fpath in args.files:
        print(f'\n{fpath}')
        try:
            probe = load_probe(fpath)
        except Exception as e:
            print(f'  ERROR: {e}')
            continue

        print(f"  dir    : {probe['dir']}")
        print(f"  comps  : {probe['comps']}")
        print(f"  npts   : {probe['npts']}")
        print(f"  nsnaps : {probe['nsnaps']}")
        print(f"  shape  : {probe['data'].shape}  (nsnaps × ncomp × npts)")

        data  = probe['data']
        comps = probe['comps']

        if data.size > 0:
            print('  --- statistics (all snapshots) ---')
            for ci, c in enumerate(comps):
                arr = data[:, ci, :]
                print(f'  {c}: min={arr.min():.4g}  max={arr.max():.4g}'
                      f'  mean={arr.mean():.4g}  rms={float(np.sqrt(np.mean(arr**2))):.4g}')

        if args.snap is not None:
            s = args.snap
            if 0 <= s < probe['nsnaps']:
                print(f'  --- snapshot {s} ---')
                np.set_printoptions(precision=4, threshold=20, edgeitems=4)
                for ci, c in enumerate(comps):
                    vals = data[s, ci, :]
                    print(f'  {c}[{probe["npts"]}]: {vals}')
            else:
                print(f'  --snap {s}: index out of range (0..{probe["nsnaps"] - 1})')

        if args.plot:
            try:
                import matplotlib.pyplot as plt
            except ImportError:
                print('  --plot requires matplotlib; skipping')
                continue
            s_last = probe['nsnaps'] - 1
            pts    = np.arange(probe['npts'])
            fig, ax = plt.subplots(figsize=(8, 4))
            for ci, c in enumerate(comps):
                ax.plot(pts, data[s_last, ci, :], label=c)
            ax.set_xlabel(f'Point index along {probe["dir"]}')
            ax.set_ylabel('Value')
            ax.set_title(f'{Path(fpath).stem}  — snapshot {s_last}')
            ax.legend()
            plt.tight_layout()
            plt.show()


if __name__ == '__main__':
    main()
