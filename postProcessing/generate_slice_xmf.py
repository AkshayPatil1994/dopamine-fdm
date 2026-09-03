#!/usr/bin/env python3
"""
NOTICE 

-- LLM GENERATED CODE BELOW --

generate_slice_xmf.py  —  XDMF time-series generator for fdm-dopamine 2-D slice probes.

Each 2-D slice probe writes a  <base>.bin  file (big-endian float64, no record
markers) with consecutive snapshots of Fortran shape (ncomp, n1, n2), and a
companion text file  <base>_meta.txt.

This script reads the meta file, writes 1-D coordinate arrays to
  <base>_ax1.bin  and  <base>_ax2.bin  (little-endian float64), and
generates  <base>.xmf  — an XDMF 2.0 temporal collection readable by
ParaView ≥ 5 and VisIt ≥ 3.

Binary layout
-------------
Fortran writes array  out2d(ncomp, n1, n2)  in column-major order.  In XDMF
the same bytes are treated as a C array of shape  (n2, n1, ncomp)  in row-
major order.  A HyperSlab selection of  Start="0 0 ci", Count="n2 n1 1"
extracts one scalar component, which maps onto the (n2 × n1) 2-D mesh.

Slice orientation
-----------------
  dir='x'  (x-normal):  axis-1 = y  (n1 = nym_global),  axis-2 = z  (n2 = nzm_global)
  dir='y'  (y-normal):  axis-1 = x  (n1 = nxm_global),  axis-2 = z  (n2 = nzm_global)
  dir='z'  (z-normal):  axis-1 = x  (n1 = nxm_global),  axis-2 = y  (n2 = nym_global)

Usage
-----
  python3 generate_slice_xmf.py <base>_meta.txt [more_meta.txt ...]
                                 [--snap SNAPSHOT_FILE]
                                 [--dt DT] [--t0 T0]

  --snap  Path to a solver snapshot from which to read exact physical
          grid coordinates (xm, ym, zm).  If omitted the script tries
          input_parameters in the case root for uniform-grid coordinates.

  --dt    Physical time between slice outputs = slice_freq × solver_dt  [default 1.0]
  --t0    Physical time of the very first snapshot  [default 0.0]

  --dt/--t0 are only a fallback: each meta file records a `times` key
  pointing at a companion <base>_times.bin (big-endian float64, one exact
  simulation time per snapshot -- written by the solver's slice-probe
  module). When that file is present, its times are used instead, so the
  resulting .xmf timeline matches the true (possibly adaptive-dt) write
  times. Point generate_UAVpath.py's --times-from at the same
  <base>_times.bin to get a UAV-disk animation whose frames land on
  exactly the same simulation times as this slice, so ParaView's shared
  time toolbar keeps the disk and the wake co-located frame-for-frame.
"""

import argparse
import os
import sys
import numpy as np
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))


# ── meta-file parser ──────────────────────────────────────────────────────────

def read_meta(path):
    """Parse a <base>_meta.txt file.  Returns a plain dict."""
    meta = {}
    for line in Path(path).read_text().splitlines():
        if '=' not in line:
            continue
        k, _, v = line.partition('=')
        k, v = k.strip(), v.strip()
        if k in ('ncomp', 'n1', 'n2', 'nsnaps'):
            meta[k] = int(v)
        elif k == 'pos':
            meta[k] = float(v)
        else:
            meta[k] = v
    if 'n1' not in meta:
        raise ValueError(f'{path}: does not look like a slice meta file (no n1 key)')
    return meta


# ── component-name helper ─────────────────────────────────────────────────────

def comp_names(s):
    """Return ordered list of component names present in the comps string."""
    upper = s.strip().upper()
    names = [c for c in ('U', 'V', 'W', 'P') if c in upper]
    return names or ['U', 'V', 'W']


# ── coordinate helpers ────────────────────────────────────────────────────────

def coords_from_snap(snap_path):
    """Read (xm, ym, zm) cell-centre coordinates from a solver snapshot or
    from fields/grid.out + fields/geometry.out (pass grid.out for speed)."""
    import snapshot_io
    if Path(snap_path).name == 'grid.out':
        return snapshot_io.coords_from_grid(snap_path)
    d = snapshot_io.read_snapshot(snap_path)
    return d['xm'], d['ym'], d['zm']


def _find_input_params(meta_path):
    """Search for input_parameters near the meta file."""
    for candidate in (
        Path(meta_path).parent / 'input_parameters',
        Path(meta_path).parent.parent / 'input_parameters',
        Path('input_parameters'),
    ):
        if candidate.is_file():
            return str(candidate)
    return None


def uniform_coords(meta_path='.'):
    """
    Approximate cell-centre coordinates from input_parameters.
    Returns (xm, ym, zm) or None if the file is not found.
    """
    ip_path = _find_input_params(meta_path)
    if ip_path is None:
        return None
    try:
        import snapshot_io
        ip = snapshot_io.parse_input_parameters(ip_path)
    except Exception:
        return None
    nx, ny, nz = int(ip.get('nx', 0)), int(ip.get('ny', 0)), int(ip.get('nz', 0))
    Lx = float(ip.get('lx', 1));  Ly = float(ip.get('ly', 1));  Lz = float(ip.get('lz', 1))
    if min(nx, ny, nz) <= 0:
        return None
    nxm, nym, nzm = nx - 1, ny - 1, nz - 1
    xm = (np.arange(nxm) + 0.5) * (Lx / nxm)
    ym = (np.arange(nym) + 0.5) * (Ly / nym)
    zm = (np.arange(nzm) + 0.5) * (Lz / nzm)
    if int(ip.get('grid_type', 1)) != 1:
        print('  NOTE: non-uniform y-grid detected; y-coordinates will be approximate.\n'
              '        Pass --snap SNAPSHOT for exact coordinates.')
    return xm, ym, zm


# ── XDMF writer ───────────────────────────────────────────────────────────────

_HEADER = """\
<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf Version="2.0">
<Domain>
  <Grid Name="{name}" GridType="Collection" CollectionType="Temporal">
"""

_FOOTER = """\
  </Grid>
</Domain>
</Xdmf>
"""


def _hyperslab_attr(cname, ci, nc, d0, d1, d2, seek, rel_bin):
    """Return XDMF lines for one scalar Attribute extracted via HyperSlab.

    Source layout (from Fortran stream write of array(nc, n1, n2) col-major)
    is treated as C array (d0, d1, d2, nc) row-major, where exactly one of
    d0, d1, d2 is 1 (the slice-normal dimension inserted for 3-D placement).
    HyperSlab Start="0 0 0 ci", Count="d0 d1 d2 1" selects component ci.
    """
    dims = f'{d0} {d1} {d2}'
    return (
        f'      <Attribute Name="{cname}" AttributeType="Scalar" Center="Node">\n'
        f'        <DataItem ItemType="HyperSlab" Dimensions="{dims}" Type="HyperSlab">\n'
        f'          <DataItem Dimensions="3 4" Format="XML">\n'
        f'            0 0 0 {ci}\n'
        f'            1 1 1 1\n'
        f'            {d0} {d1} {d2} 1\n'
        f'          </DataItem>\n'
        f'          <DataItem Dimensions="{dims} {nc}" Format="Binary"\n'
        f'                    DataType="Float" Precision="8" Endian="Big" Seek="{seek}">\n'
        f'            {rel_bin}\n'
        f'          </DataItem>\n'
        f'        </DataItem>\n'
        f'      </Attribute>\n'
    )


def make_xmf(meta_path, snap_path=None, dt=1.0, t0=0.0):
    """Generate <base>.xmf, <base>_ax1.bin, <base>_ax2.bin for one slice probe."""
    meta_path = Path(meta_path).resolve()
    stem = meta_path.stem
    # base = everything before '_meta'
    base     = meta_path.parent / (stem[:-5] if stem.endswith('_meta') else stem)
    bin_path = base.with_suffix('.bin')
    xmf_path = base.with_suffix('.xmf')
    ax1_path = base.parent / (base.name + '_ax1.bin')
    ax2_path = base.parent / (base.name + '_ax2.bin')

    meta   = read_meta(meta_path)
    nc     = meta['ncomp']
    n1     = meta['n1']
    n2     = meta['n2']
    nsnaps = meta['nsnaps']
    dirstr = meta.get('dir', 'z').strip().lower()[0]
    comps  = comp_names(meta.get('comps', 'UVW'))

    # Infer actual snapshot count from file size in case meta is stale
    if bin_path.is_file():
        snap_bytes = nc * n1 * n2 * 8
        actual = bin_path.stat().st_size // snap_bytes if snap_bytes else 0
        if actual != nsnaps:
            print(f'  NOTE: meta says {nsnaps} snaps, file has {actual}; using file count')
            nsnaps = actual

    if nsnaps == 0:
        print(f'  SKIP {meta_path.name}: no snapshots written yet')
        return

    # ── per-snapshot times ────────────────────────────────────────────────────
    # meta['times'] (written by the solver, see probe_output.f90) names a
    # companion <base>_times.bin: one big-endian float64 per snapshot, the
    # exact simulation time it was written at. Prefer that over the uniform
    # t0 + s*dt grid so this slice's animation timeline lines up exactly with
    # anything else synced to the same file (e.g. generate_UAVpath.py
    # --times-from <base>_times.bin) even when the solver's dt is adaptive.
    times = None
    times_key = meta.get('times')
    if times_key:
        # times_key is written by the solver relative to its run (case-root)
        # working directory, same as the slice_fileout base itself -- so it
        # is *not* relative to meta_path. Try, in order: alongside the
        # meta/bin files themselves (the common case, and robust to the
        # script being invoked from an unrelated cwd), then as given
        # (cwd- or absolute-relative), in case the case dir was moved and
        # only relative script invocation changed.
        candidates = [Path(times_key)] if Path(times_key).is_absolute() else [
            meta_path.parent / Path(times_key).name,
            Path(times_key),
        ]
        times_path = next((c for c in candidates if c.is_file()), candidates[0])
        if times_path.is_file():
            times = np.fromfile(times_path, dtype='>f8')
            if len(times) != nsnaps:
                print(f'  NOTE: {times_path.name} has {len(times)} times, expected {nsnaps}; '
                      'falling back to --dt/--t0')
                times = None
            else:
                print(f'  Times: from {times_path.name} (exact simulation time per snapshot)')
    if times is None:
        times = t0 + np.arange(nsnaps) * dt

    # ── physical coordinates ──────────────────────────────────────────────────
    xm = ym = zm = None

    if snap_path is not None:
        try:
            xm, ym, zm = coords_from_snap(snap_path)
            src = 'grid files' if Path(snap_path).name == 'grid.out' else f'snapshot {snap_path}'
            print(f'  Coordinates: from {src}')
        except Exception as e:
            print(f'  WARNING: snapshot read failed ({e}); trying input_parameters')

    if xm is None:
        result = uniform_coords(meta_path)
        if result is not None:
            xm, ym, zm = result
            print('  Coordinates: from input_parameters (uniform approximation)')
        else:
            print('  WARNING: no coordinates found; using integer indices')
            n_max = max(n1, n2) + 1
            xm = ym = zm = np.arange(n_max, dtype=np.float64)

    # select axis-1 (in-plane fast axis) and axis-2 (in-plane slow axis)
    if dirstr == 'x':
        ax1, ax2 = ym[:n1], zm[:n2]
        ax1_lbl, ax2_lbl = 'y', 'z'
    elif dirstr == 'y':
        ax1, ax2 = xm[:n1], zm[:n2]
        ax1_lbl, ax2_lbl = 'x', 'z'
    else:                           # 'z'
        ax1, ax2 = xm[:n1], ym[:n2]
        ax1_lbl, ax2_lbl = 'x', 'y'

    # nearest cell-centre coordinate on the slice-normal axis (for 3-D placement)
    pos = meta.get('pos', 0.0)
    if dirstr == 'x':
        ax_normal = np.array([xm[np.argmin(np.abs(xm - pos))]])
    elif dirstr == 'y':
        ax_normal = np.array([ym[np.argmin(np.abs(ym - pos))]])
    else:
        ax_normal = np.array([zm[np.argmin(np.abs(zm - pos))]])

    ax_n_path = base.parent / (base.name + '_axn.bin')

    # write coordinate binaries (little-endian, referenced with Endian="Little")
    np.asarray(ax1, dtype='<f8').tofile(ax1_path)
    np.asarray(ax2, dtype='<f8').tofile(ax2_path)
    np.asarray(ax_normal, dtype='<f8').tofile(ax_n_path)

    # relative paths (all files live in the same directory)
    xmf_dir = xmf_path.parent
    rel_bin  = os.path.relpath(bin_path, xmf_dir)
    rel_ax1  = os.path.relpath(ax1_path, xmf_dir)
    rel_ax2  = os.path.relpath(ax2_path, xmf_dir)
    rel_axn  = os.path.relpath(ax_n_path, xmf_dir)

    snap_bytes = nc * n1 * n2 * 8   # bytes per snapshot in the binary file

    # Embed slice as a 3DRectMesh with Dimensions="nz ny nx" (slowest→fastest),
    # one dimension = 1 for the slice-normal direction.  VxVyVz provides x, y, z
    # coordinate vectors; the slice is placed at its physical 3-D position.
    # Source binary (d0*d1*d2*nc floats) has the same byte layout as (n2,n1,nc)
    # since exactly one of d0,d1,d2 is 1.
    if dirstr == 'x':               # fixed x: spans y × z, embedded at x=ax_normal
        d0, d1, d2 = n2, n1, 1     # (nz, ny, 1)
        geo = (
            f'      <Geometry GeometryType="VxVyVz">\n'
            f'        <DataItem Dimensions="1" Format="Binary" DataType="Float"'
            f' Precision="8" Endian="Little">\n          {rel_axn}\n        </DataItem>\n'
            f'        <DataItem Dimensions="{n1}" Format="Binary" DataType="Float"'
            f' Precision="8" Endian="Little">\n          {rel_ax1}\n        </DataItem>\n'
            f'        <DataItem Dimensions="{n2}" Format="Binary" DataType="Float"'
            f' Precision="8" Endian="Little">\n          {rel_ax2}\n        </DataItem>\n'
            f'      </Geometry>\n'
        )
    elif dirstr == 'y':             # fixed y: spans x × z, embedded at y=ax_normal
        d0, d1, d2 = n2, 1, n1     # (nz, 1, nx)
        geo = (
            f'      <Geometry GeometryType="VxVyVz">\n'
            f'        <DataItem Dimensions="{n1}" Format="Binary" DataType="Float"'
            f' Precision="8" Endian="Little">\n          {rel_ax1}\n        </DataItem>\n'
            f'        <DataItem Dimensions="1" Format="Binary" DataType="Float"'
            f' Precision="8" Endian="Little">\n          {rel_axn}\n        </DataItem>\n'
            f'        <DataItem Dimensions="{n2}" Format="Binary" DataType="Float"'
            f' Precision="8" Endian="Little">\n          {rel_ax2}\n        </DataItem>\n'
            f'      </Geometry>\n'
        )
    else:                           # fixed z: spans x × y, embedded at z=ax_normal
        d0, d1, d2 = 1, n2, n1     # (1, ny, nx)
        geo = (
            f'      <Geometry GeometryType="VxVyVz">\n'
            f'        <DataItem Dimensions="{n1}" Format="Binary" DataType="Float"'
            f' Precision="8" Endian="Little">\n          {rel_ax1}\n        </DataItem>\n'
            f'        <DataItem Dimensions="{n2}" Format="Binary" DataType="Float"'
            f' Precision="8" Endian="Little">\n          {rel_ax2}\n        </DataItem>\n'
            f'        <DataItem Dimensions="1" Format="Binary" DataType="Float"'
            f' Precision="8" Endian="Little">\n          {rel_axn}\n        </DataItem>\n'
            f'      </Geometry>\n'
        )

    # ── build XDMF ────────────────────────────────────────────────────────────
    parts = [_HEADER.format(name=base.name)]

    for s in range(nsnaps):
        seek = s * snap_bytes
        parts.append(
            f'    <Grid Name="t{s:06d}">\n'
            f'      <Time Value="{times[s]:.6f}"/>\n'
            f'      <Topology TopologyType="3DRectMesh" Dimensions="{d0} {d1} {d2}"/>\n'
            + geo
        )
        for ci, cname in enumerate(comps):
            parts.append(_hyperslab_attr(cname, ci, nc, d0, d1, d2, seek, rel_bin))
        parts.append('    </Grid>\n')

    parts.append(_FOOTER)
    xmf_path.write_text(''.join(parts))

    print(f'  Written: {xmf_path.name}  ({nsnaps} snaps, '
          f'{ax1_lbl}[{n1}] × {ax2_lbl}[{n2}], comps={comps})')
    print(f'           {ax1_path.name}  {ax2_path.name}  {ax_n_path.name}')


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(
        description='Generate XDMF for fdm-dopamine 2-D slice probe output.')
    p.add_argument('meta_files', nargs='+', metavar='META',
                   help='Path(s) to <base>_meta.txt slice metadata files')
    p.add_argument('--snap', metavar='FILE',
                   help='fields/grid.out (fast) or a full solver snapshot — used for exact physical grid coordinates')
    p.add_argument('--dt', type=float, default=1.0,
                   help='Physical Δt between slice outputs = slice_freq × solver dt  [default: 1.0]')
    p.add_argument('--t0', type=float, default=0.0,
                   help='Physical time of first snapshot  [default: 0.0]')
    args = p.parse_args()

    for mf in args.meta_files:
        print(f'\nProcessing: {mf}')
        try:
            make_xmf(mf, snap_path=args.snap, dt=args.dt, t0=args.t0)
        except Exception as e:
            print(f'  ERROR: {e}')


if __name__ == '__main__':
    main()
