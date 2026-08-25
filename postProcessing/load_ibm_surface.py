#!/usr/bin/env python3
"""
load_ibm_surface.py — Reader/converter for fdm-dopamine IBM surface samples.

Each `ibm_surface/surface.<step>.bin` file is written by sample_ibm_surface
(src/ibm.f90).  It is a point cloud over the immersed body surface (phi=0),
one point per cell-centre interface cell.

Binary layout (big-endian, Fortran stream, no record markers):

    Int32    n            number of surface points
    float64  t            simulation time
    float64  x[n]         boundary-point coordinates
    float64  y[n]
    float64  z[n]
    float64  nx[n]        outward unit normal
    float64  ny[n]
    float64  nz[n]
    float64  p[n]         surface pressure
    float64  fp_x[n]      pressure force per point  (-p n dA)
    float64  fp_y[n]
    float64  fp_z[n]
    float64  fv_x[n]      viscous force per point   ((nu+nu_t) dU/dn dA)
    float64  fv_y[n]
    float64  fv_z[n]

Summing fp_* / fv_* over all points reproduces the Method-2 pressure / viscous
drag reported in ibm_forces.csv.

Usage
-----
    # Read into numpy:
    from load_ibm_surface import read_ibm_surface
    d = read_ibm_surface('ibm_surface/surface.00010000.bin')
    print(d['t'], d['xyz'].shape, d['p'].mean())

    # Convert one file (or a whole glob) to ParaView .vtp:
    python load_ibm_surface.py ibm_surface/surface.00010000.bin
    python load_ibm_surface.py 'ibm_surface/surface.*.bin'
"""

import sys
import glob
import struct
import numpy as np
from pathlib import Path

# big-endian: Int32 count, float64 time, then float64 blocks
_I4 = np.dtype('>i4')
_F8 = np.dtype('>f8')

_FIELDS = ('x', 'y', 'z', 'nx', 'ny', 'nz',
           'p', 'fp_x', 'fp_y', 'fp_z', 'fv_x', 'fv_y', 'fv_z')


def read_ibm_surface(fpath):
    """Read one surface .bin file. Returns a dict with keys:
       t, n, xyz (n,3), normal (n,3), p (n,),
       f_pres (n,3), f_visc (n,3), and the raw named columns.
    """
    raw = Path(fpath).read_bytes()
    off = 0
    n = struct.unpack('>i', raw[off:off + 4])[0]
    off += 4
    t = struct.unpack('>d', raw[off:off + 8])[0]
    off += 8

    cols = {}
    nbytes = n * 8
    for name in _FIELDS:
        cols[name] = np.frombuffer(raw[off:off + nbytes], dtype=_F8).astype(np.float64)
        off += nbytes

    xyz    = np.column_stack((cols['x'], cols['y'], cols['z']))
    normal = np.column_stack((cols['nx'], cols['ny'], cols['nz']))
    f_pres = np.column_stack((cols['fp_x'], cols['fp_y'], cols['fp_z']))
    f_visc = np.column_stack((cols['fv_x'], cols['fv_y'], cols['fv_z']))

    return dict(t=t, n=n, xyz=xyz, normal=normal,
                f_pres=f_pres, f_visc=f_visc, **cols)


def to_vtp(fpath, out=None):
    """Convert a surface .bin file to a VTK PolyData (.vtp) point cloud.

    Point arrays: pressure, normal, pressure_force, viscous_force,
                  total_force, pressure_force_mag, viscous_force_mag.
    Requires `pyvista` (pip install pyvista).
    """
    import pyvista as pv

    d = read_ibm_surface(fpath)
    cloud = pv.PolyData(d['xyz'])
    cloud['pressure']           = d['p']
    cloud['normal']             = d['normal']
    cloud['pressure_force']     = d['f_pres']
    cloud['viscous_force']      = d['f_visc']
    cloud['total_force']        = d['f_pres'] + d['f_visc']
    cloud['pressure_force_mag'] = np.linalg.norm(d['f_pres'], axis=1)
    cloud['viscous_force_mag']  = np.linalg.norm(d['f_visc'], axis=1)
    cloud.field_data['time'] = np.array([d['t']])

    if out is None:
        out = str(Path(fpath).with_suffix('.vtp'))
    cloud.save(out)
    return out


def to_pvd(pattern='ibm_surface/surface.*.bin', pvd_out='ibm_surface/surface.pvd'):
    """Convert all surface .bin files matching *pattern* to .vtp and write a
    ParaView Data (.pvd) collection file so the whole time series loads as an
    animation in ParaView.

    Usage
    -----
        python load_ibm_surface.py --pvd
        python load_ibm_surface.py --pvd 'run2/ibm_surface/surface.*.bin' run2/surface.pvd

    Steps in ParaView after running this:
        File > Open > surface.pvd  (group the .vtp files automatically)
        Apply, then use the time toolbar to animate.
    Requires `pyvista` (pip install pyvista).
    """
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(f'No files match: {pattern}')

    pvd_path = Path(pvd_out)
    pvd_path.parent.mkdir(parents=True, exist_ok=True)

    entries = []
    for f in files:
        vtp_out = to_vtp(f)
        d = read_ibm_surface(f)
        # path relative to the .pvd file location
        rel = Path(vtp_out).relative_to(pvd_path.parent)
        entries.append((d['t'], str(rel)))
        print(f'  {Path(f).name}  t={d["t"]:.6g}  -> {vtp_out}')

    lines = ['<?xml version="1.0"?>',
             '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">',
             '  <Collection>']
    for t, rel in entries:
        lines.append(f'    <DataSet timestep="{t:.10g}" file="{rel}"/>')
    lines += ['  </Collection>', '</VTKFile>', '']

    pvd_path.write_text('\n'.join(lines))
    print(f'\nWrote {pvd_out}  ({len(entries)} snapshots)')
    print('Open in ParaView: File > Open > surface.pvd')
    return pvd_out


def _summarise(fpath):
    d = read_ibm_surface(fpath)
    Fp = d['f_pres'].sum(axis=0)
    Fv = d['f_visc'].sum(axis=0)
    print(f'{Path(fpath).name}: t={d["t"]:.6g}  n={d["n"]}  '
          f'Fpres=({Fp[0]:.4e},{Fp[1]:.4e},{Fp[2]:.4e})  '
          f'Fvisc=({Fv[0]:.4e},{Fv[1]:.4e},{Fv[2]:.4e})')


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(0)

    # --pvd mode: convert full time series and write a .pvd collection
    if sys.argv[1] == '--pvd':
        pattern = sys.argv[2] if len(sys.argv) > 2 else 'ibm_surface/surface.*.bin'
        pvd_out = sys.argv[3] if len(sys.argv) > 3 else 'ibm_surface/surface.pvd'
        to_pvd(pattern, pvd_out)
        sys.exit(0)

    files = []
    for arg in sys.argv[1:]:
        files.extend(sorted(glob.glob(arg)))
    if not files:
        print('No matching files.')
        sys.exit(1)

    for f in files:
        _summarise(f)
        try:
            out = to_vtp(f)
            print(f'   wrote {out}')
        except ImportError:
            print('   (install pyvista to also write .vtp)')
