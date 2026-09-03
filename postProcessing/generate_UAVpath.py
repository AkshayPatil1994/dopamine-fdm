#!/usr/bin/env python3
"""
generate_UAVpath.py — Render a fdm-dopamine UAV actuator-disk path
(uav_path_file, see src/uav_actuator.f90 / docs/Input-Parameters.md's &UAV
group) as a moving-disk animation, so it can be loaded into ParaView
alongside the flow-field snapshots and visually overlaid on the wind data.

Interpolation exactly reproduces the solver's own uav_current_center /
hermite_eval (cubic Hermite, Catmull-Rom tangents, clamped outside the
waypoint file's time range) -- see uav_actuator.f90 -- so the rendered disk
sits exactly where the actuator-disk force was actually applied in the
simulation, not merely close to it.

Phase 1/2 scope note: the solver's disk stays horizontal (normal = +y) for
its whole path (no cruise-segment tilt yet, see docs/UAV_ActuatorDisk_
Design.md); this script matches that and always draws a flat, horizontal
disk. Update it alongside uav_actuator.f90 if/when tilt is added.

Output: one VTK PolyData (.vtp) disc per requested time, plus a ParaView
Data (.pvd) collection so the whole sequence loads as a single animated
object -- open uav_path.pvd (or your --pvd-out path) in ParaView alongside
the flow-field .pvd/.xmf (see generateXMF.py / generate_slice_xmf.py) and
use the shared time toolbar to animate both together.

Requires `pyvista` (pip install pyvista).

Usage
-----
    # Uniform time grid spanning the path file's own range, 50 frames:
    python generate_UAVpath.py uav_path_file.dat --radius 0.15

    # Sync frames exactly to an existing flow-field output's recorded
    # times (a slice probe's companion _times.bin -- line probes don't
    # write one -- or any plain-text file with one time per line):
    python generate_UAVpath.py uav_path_file.dat --radius 0.15 \\
        --times-from slices/y015_times.bin

    # Pull uav_disk_radius/uav_n_theta straight from the case's own
    # input_parameters instead of passing --radius by hand:
    python generate_UAVpath.py uav_path_file.dat --input-parameters input_parameters

Then in ParaView: File > Open > uav_path.pvd, Apply, and use the time
toolbar (shared with any other .pvd/.xmf you have open) to animate.
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np


def read_uav_path(fpath):
    """Read a uav_path_file: rows 't x y z', blank/'#' lines skipped."""
    rows = []
    with open(fpath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            rows.append([float(v) for v in line.split()[:4]])
    if len(rows) < 2:
        raise ValueError(f'{fpath}: fewer than 2 waypoints found')
    arr = np.array(rows)
    t, x, y, z = arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3]
    if np.any(np.diff(t) <= 0):
        raise ValueError(f'{fpath}: t column must be strictly increasing')
    return t, x, y, z


def hermite_eval(tt, vals_t, vals_p):
    """Cubic Hermite (Catmull-Rom tangent) interpolation, one component,
    at scalar or array time tt -- bit-for-bit the same construction as
    uav_actuator.f90's hermite_eval (clamped outside the endpoints,
    one-sided tangents at the first/last segment).
    """
    tt = np.atleast_1d(np.asarray(tt, dtype=float))
    n = len(vals_t)
    out = np.empty_like(tt)

    out[tt <= vals_t[0]] = vals_p[0]
    out[tt >= vals_t[-1]] = vals_p[-1]
    inside = (tt > vals_t[0]) & (tt < vals_t[-1])

    # bracket index for every interior point (vals_t[i] <= t < vals_t[i+1])
    idx = np.searchsorted(vals_t, tt[inside], side='right') - 1
    idx = np.clip(idx, 0, n - 2)

    t1, t2 = vals_t[idx], vals_t[idx + 1]
    p1, p2 = vals_p[idx], vals_p[idx + 1]
    h = t2 - t1
    s = (tt[inside] - t1) / h

    has_prev = idx > 0
    m1 = np.where(
        has_prev,
        (p2 - np.where(has_prev, vals_p[np.clip(idx - 1, 0, n - 1)], 0.0)) /
        np.where(has_prev, (t2 - np.where(has_prev, vals_t[np.clip(idx - 1, 0, n - 1)], 0.0)), 1.0),
        (p2 - p1) / h,
    )
    has_next = idx + 2 <= n - 1
    m2 = np.where(
        has_next,
        (np.where(has_next, vals_p[np.clip(idx + 2, 0, n - 1)], 0.0) - p1) /
        np.where(has_next, (np.where(has_next, vals_t[np.clip(idx + 2, 0, n - 1)], 0.0) - t1), 1.0),
        (p2 - p1) / h,
    )

    h00 = 2 * s**3 - 3 * s**2 + 1
    h10 = s**3 - 2 * s**2 + s
    h01 = -2 * s**3 + 3 * s**2
    h11 = s**3 - s**2

    out[inside] = h00 * p1 + h10 * h * m1 + h01 * p2 + h11 * h * m2
    return out


def disk_centre(times, path_t, path_x, path_y, path_z):
    xc = hermite_eval(times, path_t, path_x)
    yc = hermite_eval(times, path_t, path_y)
    zc = hermite_eval(times, path_t, path_z)
    return xc, yc, zc


def read_times_file(fpath):
    """Load output times: a numpy .bin (big-endian float64, e.g. a probe's
    own <fileout>_times.bin, see docs/Input-Parameters.md &STATISTICS) if
    the extension is .bin, else a plain text file with one time per line.
    """
    fpath = Path(fpath)
    if fpath.suffix == '.bin':
        return np.fromfile(fpath, dtype='>f8')
    return np.loadtxt(fpath).reshape(-1)


def read_uav_params(input_parameters):
    """Pull uav_disk_radius/uav_n_theta out of an input_parameters file's
    &UAV group with a small regex parser (good enough for the plain
    'key = value,' namelist rows this solver's input files use).
    """
    text = Path(input_parameters).read_text()
    m = re.search(r'&UAV\b(.*?)/', text, re.S | re.I)
    if not m:
        raise ValueError(f'{input_parameters}: no &UAV group found')
    block = m.group(1)

    def _get(name, default):
        mm = re.search(rf'\b{name}\s*=\s*([0-9.eE+-]+)', block)
        return float(mm.group(1)) if mm else default

    radius = _get('uav_disk_radius', 0.15)
    n_theta = int(_get('uav_n_theta', 48))
    return radius, n_theta


def write_pvd(entries, pvd_out):
    """entries: list of (time, relative_vtp_path)."""
    pvd_path = Path(pvd_out)
    pvd_path.parent.mkdir(parents=True, exist_ok=True)
    lines = ['<?xml version="1.0"?>',
             '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">',
             '  <Collection>']
    for t, rel in entries:
        lines.append(f'    <DataSet timestep="{t:.10g}" file="{rel}"/>')
    lines += ['  </Collection>', '</VTKFile>', '']
    pvd_path.write_text('\n'.join(lines))


def generate(path_file, times, radius, n_theta, out_dir, pvd_out):
    import pyvista as pv

    path_t, path_x, path_y, path_z = read_uav_path(path_file)
    xc, yc, zc = disk_centre(times, path_t, path_x, path_y, path_z)

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    pvd_path = Path(pvd_out)

    entries = []
    for i, (t, x, y, z) in enumerate(zip(times, xc, yc, zc)):
        # Phase 1/2: disk stays horizontal, normal=+y (see module docstring)
        disc = pv.Disc(center=(x, y, z), inner=0.0, outer=radius,
                        normal=(0.0, 1.0, 0.0), r_res=1, c_res=n_theta)
        disc.field_data['time'] = np.array([t])
        disc.field_data['centre'] = np.array([[x, y, z]])
        vtp_path = out_dir / f'uav_disk.{i:05d}.vtp'
        disc.save(str(vtp_path))
        rel = vtp_path.relative_to(pvd_path.parent)
        entries.append((t, str(rel)))

    write_pvd(entries, pvd_path)
    return pvd_path, len(entries)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('path_file', help='uav_path_file (rows: t x y z)')
    ap.add_argument('--radius', type=float, default=None,
                     help='disk radius [m] (default 0.15, or read from --input-parameters)')
    ap.add_argument('--n-theta', type=int, default=None,
                     help='circumferential resolution (default 48, or read from --input-parameters)')
    ap.add_argument('--input-parameters', default=None,
                     help="read uav_disk_radius/uav_n_theta from this case's &UAV namelist group "
                          '(overridden by --radius/--n-theta if also given)')
    ap.add_argument('--times-from', default=None,
                     help='sync frame times to this file: a probe _times.bin (big-endian float64) '
                          'or a plain text file with one time per line')
    ap.add_argument('--t0', type=float, default=None, help='first frame time (default: path start)')
    ap.add_argument('--t1', type=float, default=None, help='last frame time (default: path end)')
    ap.add_argument('--nframes', type=int, default=50, help='frame count for a uniform time grid (default 50)')
    ap.add_argument('--out-dir', default='uav_path', help='directory for the per-frame .vtp files')
    ap.add_argument('--pvd-out', default='uav_path.pvd', help='output .pvd collection path')
    args = ap.parse_args()

    radius, n_theta = 0.15, 48
    if args.input_parameters:
        radius, n_theta = read_uav_params(args.input_parameters)
    if args.radius is not None:
        radius = args.radius
    if args.n_theta is not None:
        n_theta = args.n_theta

    if args.times_from:
        times = read_times_file(args.times_from)
    else:
        path_t, _, _, _ = read_uav_path(args.path_file)
        t0 = args.t0 if args.t0 is not None else path_t[0]
        t1 = args.t1 if args.t1 is not None else path_t[-1]
        times = np.linspace(t0, t1, args.nframes)

    pvd_path, n = generate(args.path_file, times, radius, n_theta, args.out_dir, args.pvd_out)
    print(f'Wrote {n} frames -> {pvd_path}')
    print('Open in ParaView: File > Open > ' + str(pvd_path))
