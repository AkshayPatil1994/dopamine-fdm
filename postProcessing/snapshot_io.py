#!/usr/bin/env python3
"""
snapshot_io.py — Reader library for fdm-dopamine binary field snapshots.

Binary layout (big-endian float64, Fortran stream, no record markers):

  Grid header — 6 blocks of (Int32 count, float64[count] coordinates):
      nx  face points in x  →  x
      ny  face points in y  →  y
      nz  face points in z  →  z
      nxm cell centres in x →  xm
      nym cell centres in y →  ym
      nzm cell centres in z →  zm

  Field blocks — each block: (3 × Int32 shape), (float64[n1×n2×n3] data, Fortran/column-major order):
      U    (nx,  nyg, nzg)   x-face velocity         — always present
      V    (nxg, ny,  nzg)   y-face velocity         — always present
      W    (nxg, nyg, nz )   z-face velocity         — always present
      P    (nxg, nyg, nzg)   cell-centre pressure    — always present
      C    (nxg, nyg, nzg)   scalar concentration    — only if sediment_flag >= 1
      nu_t (nxg, nyg, nzg)   SGS turbulent viscosity — only if sgs_model != 0

  Ghost-cell convention: nxg = nxm+2, nyg = nym+2, nzg = nzm+2.
  ny and nz in the V and W headers are face-point counts (no ghost cells on
  the staggered direction), so ny = nym+1, nz = nzm+1.

  Returned velocity/pressure/scalar arrays are cell-centred with all ghost
  layers removed: shape (nxm, nym, nzm).
  Axis convention: 0 = x (streamwise), 1 = y (wall-normal), 2 = z (spanwise).
"""

import re
import struct
import numpy as np
from pathlib import Path


# ── namelist parser ───────────────────────────────────────────────────────────

def parse_input_parameters(fpath='input_parameters'):
    """
    Minimal Fortran-namelist parser.

    Returns a flat dict of lowercase key → Python value.
    Supports: integers, floats (including Fortran 'd'/'D' exponent notation),
    booleans (.TRUE./.FALSE.) and single-quoted strings.
    Comment lines (starting with !) are stripped before parsing.
    """
    text = Path(fpath).read_text()
    text = re.sub(r'!.*', '', text)      # strip ! comments
    text = text.replace('\n', ' ')       # flatten to one line

    params: dict = {}
    for m in re.finditer(r'(\w+)\s*=\s*([^,/]+)', text):
        key = m.group(1).strip().lower()
        raw = m.group(2).strip()
        if raw.upper() in ('.TRUE.', '.T.'):
            params[key] = True
        elif raw.upper() in ('.FALSE.', '.F.'):
            params[key] = False
        elif raw.startswith("'"):
            inner = re.match(r"'([^']*)'", raw)
            params[key] = inner.group(1).strip() if inner else raw
        else:
            norm = raw.replace('d', 'e').replace('D', 'e')
            try:
                params[key] = int(norm)
            except ValueError:
                try:
                    params[key] = float(norm)
                except ValueError:
                    params[key] = raw
    return params


def list_snapshots(fields_dir='fields', fileout_prefix='channel_test'):
    """
    Return a sorted list of (step, filepath) tuples for all snapshots whose
    filename matches '<fileout_prefix>.<integer>' inside fields_dir.
    """
    results = []
    for f in Path(fields_dir).iterdir():
        # filename pattern: <prefix>.<step>
        if f.name.startswith(fileout_prefix + '.'):
            suffix = f.name[len(fileout_prefix) + 1:]
            try:
                results.append((int(suffix), str(f)))
            except ValueError:
                pass
    return sorted(results)


# ── low-level binary reader ───────────────────────────────────────────────────

class _Reader:
    """Stateful reader for big-endian Fortran stream binary files."""

    def __init__(self, data: bytes):
        self._d = data
        self._p = 0

    def ri(self) -> int:
        v = struct.unpack_from('>i', self._d, self._p)[0]
        self._p += 4
        return v

    def ri3(self):
        v = struct.unpack_from('>3i', self._d, self._p)
        self._p += 12
        return v

    def rd(self, n: int) -> np.ndarray:
        arr = np.frombuffer(self._d, dtype='>f8', count=n, offset=self._p).copy()
        self._p += n * 8
        return arr

    def read_field(self):
        """
        Read one field block: 3-Int32 shape header followed by float64 data.
        Returns (shape_tuple, ndarray) with Fortran (column-major) axis order.
        """
        n1, n2, n3 = self.ri3()
        data = self.rd(n1 * n2 * n3).reshape((n1, n2, n3), order='F')
        return (n1, n2, n3), data


# ── ghost-layer stripping and face → cell-centre interpolation ────────────────

def _to_cell_centre(U_face, V_face, W_face, P_raw):
    """
    Interpolate staggered face values to cell centres and strip ghost layers.

        U_face : (nx,  nyg, nzg)   x-face velocity
        V_face : (nxg, ny,  nzg)   y-face velocity   (ny = nym+1 face points)
        W_face : (nxg, nyg, nz )   z-face velocity   (nz = nzm+1 face points)
        P_raw  : (nxg, nyg, nzg)   cell-centre pressure, ghost ring on all sides

    Returns U, V, W, P each of shape (nxm, nym, nzm).
    """
    # x-face → cell centre: average adjacent x-faces, strip y/z ghost layers
    U = 0.5 * (U_face[:-1, 1:-1, 1:-1] + U_face[1:,  1:-1, 1:-1])
    # y-face → cell centre: average adjacent y-faces, strip x/z ghost layers
    V = 0.5 * (V_face[1:-1, :-1,  1:-1] + V_face[1:-1,  1:, 1:-1])
    # z-face → cell centre: average adjacent z-faces, strip x/y ghost layers
    W = 0.5 * (W_face[1:-1, 1:-1, :-1 ] + W_face[1:-1, 1:-1,  1:])
    # pressure: strip one ghost layer on every side
    P = P_raw[1:-1, 1:-1, 1:-1]
    return U, V, W, P


# ── public API ────────────────────────────────────────────────────────────────

def read_snapshot(fpath, *, sgs_model=0, sediment_flag=0):
    """
    Read one fdm-dopamine binary snapshot and return cell-centred fields.

    Parameters
    ----------
    fpath : str or Path
        Path to the snapshot file.
    sgs_model : int
        0 = DNS (no nu_t block present).  Non-zero = LES (nu_t block appended
        after P, or after C if sediment is also active).
    sediment_flag : int
        0 = no scalar transport.  >= 1 = C block present after P.

    Returns
    -------
    dict with the following keys:

        x,  y,  z    — 1-D face-point coordinate arrays (float64)
        xm, ym, zm   — 1-D cell-centre coordinate arrays (float64)
        U, V, W      — cell-centred velocity components, shape (nxm, nym, nzm)
        P            — cell-centred pressure,            shape (nxm, nym, nzm)
        C            — scalar concentration,             shape (nxm, nym, nzm)
                       [present only when sediment_flag >= 1]
        nu_t         — SGS turbulent viscosity,          shape (nxm, nym, nzm)
                       [present only when sgs_model != 0]

    Axis convention: 0 = x (streamwise), 1 = y (wall-normal), 2 = z (spanwise).
    """
    r = _Reader(Path(fpath).read_bytes())

    # ── grid ─────────────────────────────────────────────────────────────────
    nx  = r.ri();  x  = r.rd(nx)
    ny  = r.ri();  y  = r.rd(ny)
    nz  = r.ri();  z  = r.rd(nz)
    nxm = r.ri();  xm = r.rd(nxm)
    nym = r.ri();  ym = r.rd(nym)
    nzm = r.ri();  zm = r.rd(nzm)

    # ── always-present fields ─────────────────────────────────────────────────
    _, U_face = r.read_field()    # (nx,  nyg, nzg)
    _, V_face = r.read_field()    # (nxg, ny,  nzg)
    _, W_face = r.read_field()    # (nxg, nyg, nz )
    _, P_raw  = r.read_field()    # (nxg, nyg, nzg)

    U, V, W, P = _to_cell_centre(U_face, V_face, W_face, P_raw)

    snap = dict(x=x, y=y, z=z, xm=xm, ym=ym, zm=zm,
                U=U, V=V, W=W, P=P)

    # ── optional: scalar concentration ───────────────────────────────────────
    if sediment_flag >= 1:
        _, C_raw = r.read_field()        # (nxg, nyg, nzg)
        snap['C'] = C_raw[1:-1, 1:-1, 1:-1]

    # ── optional: SGS turbulent viscosity ─────────────────────────────────────
    if sgs_model != 0:
        _, nut_raw = r.read_field()      # (nxg, nyg, nzg)
        snap['nu_t'] = nut_raw[1:-1, 1:-1, 1:-1]

    return snap


def coords_from_grid(grid_path):
    """
    Read (xm, ym, zm) cell-centre coordinates from fields/grid.out and
    fields/geometry.out — much faster than reading a full snapshot.

    grid.out format (written by genGridandIC.f90)::

        i  ym_i  y_face_{i+1}  dy_i  1/dy_i

    geometry.out format::

        nxm  nym  nzm
        Lx   Ly   Lz

    x and z grids are always uniform in fdm-dopamine.
    """
    grid_path = Path(grid_path)
    geo_path  = grid_path.parent / 'geometry.out'
    if not geo_path.is_file():
        raise FileNotFoundError(f'geometry.out not found alongside {grid_path}')

    tokens = geo_path.read_text().split()
    nxm, nym, nzm = int(tokens[0]), int(tokens[1]), int(tokens[2])
    Lx,  Ly,  Lz  = float(tokens[3]), float(tokens[4]), float(tokens[5])

    data = np.loadtxt(grid_path)
    ym   = data[:, 1]          # column 1 (0-based): cell-centre y

    xm = (np.arange(nxm) + 0.5) * (Lx / nxm)
    zm = (np.arange(nzm) + 0.5) * (Lz / nzm)
    return xm, ym, zm
