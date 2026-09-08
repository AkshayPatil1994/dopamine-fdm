#!/usr/bin/env python3
"""
generate_UAVdrone.py -- Build a scaled quadcopter STL for a fdm-dopamine UAV
case, plus a per-snapshot rigid-body transform table, so the drone *body*
(not just its actuator disk, cf. generate_UAVpath.py) can be animated in
ParaView synced with the flow-field snapshots.

Only ONE copy of the drone geometry is ever written to disk. Motion is
carried entirely by a small per-time (position, 3x3 rotation) table
(--transforms-out, .npz + a human-readable .csv) -- there is no per-frame
mesh duplication. A ready-to-paste ParaView "Programmable Source" script
(--pv-source-out) reads the STL once and re-poses it in memory for whatever
the current animation time is.

Geometry
--------
Reuses a simple box/cylinder quadcopter model (central body, 4 diagonal
arms, motor housings, propeller discs, landing legs), boolean-unioned by
trimesh into one watertight mesh. The whole assembly is scaled uniformly so
its stock propeller radius (55 mm in the reference model) matches this
case's uav_disk_radius [m] (&UAV, see docs/Input-Parameters.md), then
rotated once into the solver's axis convention (x streamwise, y wall-
normal/"up", z spanwise) and written to --stl-out.

Position
--------
Exactly the solver's own uav_current_center: cubic-Hermite (Catmull-Rom
tangent) interpolation of uav_path_file, clamped outside the file's time
range -- see hermite_eval, reused from generate_UAVpath.py.

Orientation
-----------
--tilt off (default, or uav_tilt_active=0 in &UAV): body stays level,
y-up, identical to the solver's flat actuator-disk model.

--tilt on (only meaningful with uav_tilt_active=1): replays the *same*
recursive low-pass filter as uav_disk_state (uav_actuator.f90) over a
dense grid (this case's &NUMERICS dt) from the path's start time up to
each requested output time. The solver actually calls this filter once per
RK sub-stage (finer, and at slightly different instants) rather than once
per dt, so this is a close but not bit-exact reproduction of the in-run
tilt -- good enough to visualize, not to re-derive the applied force.

Yaw (heading about the vertical) has no counterpart in the solver's
actuator-disk model at all (a disk has no heading). By default this script
points the body's nose along the path's own horizontal velocity direction
purely for a plausible-looking animation; pass --no-heading to keep the
body's arms fixed in the X configuration instead.

Requires `trimesh` (pip install trimesh; `manifold3d` optional for a
proper boolean union, falls back to concatenation otherwise) and `numpy`.

Usage
-----
    # Geometry + a transform table synced to a slice probe's exact times:
    python3 postProcessing/generate_UAVdrone.py uav_path_file.dat \\
        --input-parameters input_parameters --times-from slices/y015_times.bin

    # Uniform 100-frame table spanning the path file's own time range:
    python3 postProcessing/generate_UAVdrone.py uav_path_file.dat \\
        --input-parameters input_parameters --nframes 100 --tilt

Then in ParaView: Sources > Programmable Source, paste in the generated
uav_drone_paraview_source.py, Apply -- open the flow-field .pvd/.xmf
alongside it and use the shared time toolbar; the Programmable Source
re-poses the single STL mesh from uav_transforms.npz on every time step.
"""

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from generate_UAVpath import read_uav_path, hermite_eval, read_times_file  # noqa: E402
from snapshot_io import parse_input_parameters  # noqa: E402


# ---------------------------------------------------------------
# Reference quadcopter dimensions [m] (the box/cylinder model given in the
# task prompt, converted from its native mm to m); scaled uniformly so
# PROP_RADIUS_REF maps onto this case's uav_disk_radius.
# ---------------------------------------------------------------
PROP_RADIUS_REF = 0.055
PROP_HEIGHT_REF = 0.002
LEG_RADIUS_REF = 0.004
LEG_LENGTH_REF = 0.030

# Rotates the reference model's local frame (arms in its xy-plane, z "up")
# into the solver's axis convention (x streamwise, y wall-normal/"up", z
# spanwise): local (x,y,z) -> solver (x, z, -y). Proper rotation (det=+1).
SOLVER_AXES = np.array([
    [1.0, 0.0, 0.0],
    [0.0, 0.0, 1.0],
    [0.0, -1.0, 0.0],
])


def build_drone_mesh(scale):
    """Build the watertight single-propeller+legs mesh at the given uniform scale
    (dimensionless multiplier on the *_REF dimensions above), already
    rotated into the solver's (x, y-up, z) axis convention. Origin is the
    body/disk centre, matching (uav_xc,uav_yc,uav_zc)/uav_path_file.
    """
    import trimesh

    s = scale
    prop_radius = PROP_RADIUS_REF * s
    prop_height = PROP_HEIGHT_REF * s
    leg_radius = LEG_RADIUS_REF * s
    leg_length = LEG_LENGTH_REF * s

    parts = []

    hub_radius = 0.18 * prop_radius
    hub_height = 4.0 * prop_height
    hub = trimesh.creation.cylinder(radius=hub_radius, height=hub_height, sections=24)
    parts.append(hub)

    blade_width = 0.16 * prop_radius
    blade = trimesh.creation.box(extents=(2.0 * prop_radius, blade_width, prop_height))
    parts.append(blade)

    ring_radius = 1.15 * prop_radius
    ring_width = 0.10 * prop_radius
    ring = trimesh.creation.annulus(
        r_min=ring_radius - ring_width, r_max=ring_radius, height=prop_height, sections=48)
    parts.append(ring)

    def strut_between(p0, p1, radius):
        p0, p1 = np.asarray(p0, dtype=float), np.asarray(p1, dtype=float)
        vec = p1 - p0
        length = np.linalg.norm(vec)
        strut = trimesh.creation.cylinder(radius=radius, height=length, sections=12)
        strut.apply_transform(trimesh.geometry.align_vectors([0, 0, 1], vec))
        strut.apply_translation((p0 + p1) / 2)
        return strut

    leg_offset = 0.7 * prop_radius
    leg_angles_deg = [45, 135, 225, 315]
    leg_top_z = -(hub_height / 2 - 0.002 * s)
    leg_z = -(hub_height / 2 + leg_length / 2 - 0.002 * s)
    for ang in leg_angles_deg:
        ang_rad = np.radians(ang)
        lx = leg_offset * np.cos(ang_rad)
        ly = leg_offset * np.sin(ang_rad)
        leg = trimesh.creation.cylinder(radius=leg_radius, height=leg_length, sections=16)
        leg.apply_translation((lx, ly, leg_z))
        parts.append(leg)

        rx, ry = ring_radius * np.cos(ang_rad), ring_radius * np.sin(ang_rad)
        parts.append(strut_between((rx, ry, 0.0), (lx, ly, leg_top_z), leg_radius))

    try:
        drone = trimesh.boolean.union(parts, engine="manifold")
    except Exception as e:
        print(f"Boolean union failed ({e}); falling back to concatenation.", file=sys.stderr)
        drone = trimesh.util.concatenate(parts)

    drone.remove_unreferenced_vertices()
    drone.merge_vertices()

    base_rot = np.eye(4)
    base_rot[:3, :3] = SOLVER_AXES
    drone.apply_transform(base_rot)

    return drone


# ---------------------------------------------------------------
# Orientation: reproduces uav_disk_state's tilt low-pass filter and
# uav_disk_frame (uav_actuator.f90) in Python.
# ---------------------------------------------------------------

def hermite_accel(tt, vals_t, vals_p):
    """Second derivative of the same cubic-Hermite segment as hermite_eval
    (generate_UAVpath.py) -- mirrors uav_actuator.f90's hermite_accel.
    Zero outside [vals_t[0], vals_t[-1]] (position is clamped there).
    """
    tt = np.atleast_1d(np.asarray(tt, dtype=float))
    n = len(vals_t)
    out = np.zeros_like(tt)
    inside = (tt > vals_t[0]) & (tt < vals_t[-1])
    if not np.any(inside):
        return out

    idx = np.searchsorted(vals_t, tt[inside], side='right') - 1
    idx = np.clip(idx, 0, n - 2)

    t1, t2 = vals_t[idx], vals_t[idx + 1]
    p1, p2 = vals_p[idx], vals_p[idx + 1]
    h = t2 - t1
    s = (tt[inside] - t1) / h

    has_prev = idx > 0
    idx_prev = np.clip(idx - 1, 0, n - 1)
    m1 = np.where(has_prev,
                  (p2 - vals_p[idx_prev]) / (t2 - vals_t[idx_prev]),
                  (p2 - p1) / h)

    has_next = idx + 2 <= n - 1
    idx_next = np.clip(idx + 2, 0, n - 1)
    m2 = np.where(has_next,
                  (vals_p[idx_next] - p1) / (vals_t[idx_next] - t1),
                  (p2 - p1) / h)

    out[inside] = ((12 * s - 6) * p1 + (6 * s - 4) * h * m1 +
                   (-12 * s + 6) * p2 + (6 * s - 2) * h * m2) / (h * h)
    return out


def disk_frame(nvec):
    """uav_disk_frame: orthonormal in-plane frame (e1,e2) for a unit disk
    normal nvec, degenerating only within a few degrees of +-xhat.
    """
    norm2 = np.hypot(nvec[1], nvec[2])
    if norm2 < 1e-8:
        return np.array([1.0, 0.0, 0.0]), np.array([0.0, 0.0, 1.0])
    e2 = np.array([0.0, -nvec[2] / norm2, nvec[1] / norm2])
    e1 = np.cross(nvec, e2)
    return e1, e2


def tilt_frames(sample_t, path_t, path_x, path_y, path_z, tau, grav):
    """Sequentially replay uav_disk_state's tilt low-pass filter over the
    (already sorted, strictly increasing) times in sample_t. Returns
    nvec/e1/e2 arrays, one row per sample_t entry.
    """
    n = len(sample_t)
    nvecs = np.empty((n, 3))
    e1s = np.empty((n, 3))
    e2s = np.empty((n, 3))

    ax = hermite_accel(sample_t, path_t, path_x)
    ay = hermite_accel(sample_t, path_t, path_y)
    az = hermite_accel(sample_t, path_t, path_z)

    smooth = np.array([ax[0], grav + ay[0], az[0]])
    norm = np.linalg.norm(smooth)
    if norm > 1e-30:
        smooth = smooth / norm
    e1, e2 = disk_frame(smooth)
    nvecs[0], e1s[0], e2s[0] = smooth, e1, e2

    for i in range(1, n):
        n_raw = np.array([ax[i], grav + ay[i], az[i]])
        norm_raw = np.linalg.norm(n_raw)
        if norm_raw > 1e-30:
            n_raw = n_raw / norm_raw
        dt_call = sample_t[i] - sample_t[i - 1]
        alpha = dt_call / (tau + dt_call)
        smooth = smooth + alpha * (n_raw - smooth)
        smooth = smooth / np.linalg.norm(smooth)
        e1, e2 = disk_frame(smooth)
        nvecs[i], e1s[i], e2s[i] = smooth, e1, e2

    return nvecs, e1s, e2s


def heading_yaw(times, path_t, path_x, path_z, eps):
    """Horizontal heading angle (about +y) from the path's own velocity
    direction, via a central-difference derivative of hermite_eval.
    Visualization-only -- the actuator-disk model has no heading.
    """
    vx = (hermite_eval(times + eps, path_t, path_x) -
          hermite_eval(times - eps, path_t, path_x)) / (2 * eps)
    vz = (hermite_eval(times + eps, path_t, path_z) -
          hermite_eval(times - eps, path_t, path_z)) / (2 * eps)
    speed = np.hypot(vx, vz)
    psi = np.where(speed > 1e-9, np.arctan2(vz, vx), 0.0)
    return psi


def compute_transforms(times, path_t, path_x, path_y, path_z,
                        tilt_active, tau, dt, grav, heading):
    """Per-time 4x4 rigid-body transforms (solver axis convention)."""
    xc = hermite_eval(times, path_t, path_x)
    yc = hermite_eval(times, path_t, path_y)
    zc = hermite_eval(times, path_t, path_z)

    n = len(times)
    if tilt_active:
        t_start = min(path_t[0], times.min())
        t_end = max(path_t[-1], times.max())
        n_dense = max(int(np.ceil((t_end - t_start) / dt)) + 1, 2)
        dense = t_start + dt * np.arange(n_dense)
        dense = dense[dense <= t_end]
        if dense[-1] < t_end:
            dense = np.append(dense, t_end)
        grid = np.union1d(dense, times)
        nvecs, e1s, e2s = tilt_frames(grid, path_t, path_x, path_y, path_z, tau, grav)
        sel = np.searchsorted(grid, times)
        nvecs, e1s, e2s = nvecs[sel], e1s[sel], e2s[sel]
    else:
        nvecs = np.tile([0.0, 1.0, 0.0], (n, 1))
        e1s = np.tile([1.0, 0.0, 0.0], (n, 1))
        e2s = np.tile([0.0, 0.0, 1.0], (n, 1))

    if heading:
        span = max(path_t[-1] - path_t[0], 1e-6)
        psi = heading_yaw(times, path_t, path_x, path_z, eps=1e-5 * span)
    else:
        psi = np.zeros(n)

    mats = np.zeros((n, 4, 4))
    mats[:, 3, 3] = 1.0
    mats[:, 0, 3], mats[:, 1, 3], mats[:, 2, 3] = xc, yc, zc

    cpsi, spsi = np.cos(psi), np.sin(psi)
    for i in range(n):
        r_yaw = np.array([[cpsi[i], 0.0, spsi[i]],
                           [0.0, 1.0, 0.0],
                           [-spsi[i], 0.0, cpsi[i]]])
        r_tilt = np.column_stack([e1s[i], nvecs[i], e2s[i]])
        mats[i, :3, :3] = r_tilt @ r_yaw

    return mats


def write_transforms(times, mats, npz_out, csv_out):
    np.savez(npz_out, times=times, matrices=mats)

    header = 'time,x,y,z,' + ','.join(f'r{i}{j}' for i in range(3) for j in range(3))
    rows = []
    for t, m in zip(times, mats):
        rot = m[:3, :3].flatten()
        rows.append([t, m[0, 3], m[1, 3], m[2, 3], *rot])
    np.savetxt(csv_out, rows, delimiter=',', header=header, comments='')


PV_SOURCE_TEMPLATE = '''\
# Auto-generated by generate_UAVdrone.py -- paste this whole file into a
# ParaView "Programmable Source" (Sources > Programmable Source, Output
# Type: vtkPolyData), then Apply. The single STL below is read once and
# re-posed in memory at every animation time step from the transform
# table -- no per-frame geometry is ever written to disk.

import vtk
import numpy as np

STL_PATH = r"{stl_path}"
TRANSFORMS_PATH = r"{npz_path}"

_data = np.load(TRANSFORMS_PATH)
_times = _data['times']
_mats = _data['matrices']

_reader = vtk.vtkSTLReader()
_reader.SetFileName(STL_PATH)
_reader.Update()
_base = _reader.GetOutput()

_executive = self.GetExecutive()
_outInfo = _executive.GetOutputInformation(0)
if _outInfo.Has(_executive.UPDATE_TIME_STEP()):
    t = _outInfo.Get(_executive.UPDATE_TIME_STEP())
else:
    t = float(_times[0])
i = int(np.argmin(np.abs(_times - t)))

transform = vtk.vtkTransform()
transform.SetMatrix(_mats[i].flatten().tolist())

tf = vtk.vtkTransformPolyDataFilter()
tf.SetInputData(_base)
tf.SetTransform(transform)
tf.Update()

output.ShallowCopy(tf.GetOutput())
'''

PV_SOURCE_INFO_TEMPLATE = '''\
# Paste into the same Programmable Source's "Script (RequestInformation)" box.

import vtk
import numpy as np

_data = np.load(r"{npz_path}")
_times = _data['times']

executive = self.GetExecutive()
outInfo = executive.GetOutputInformation(0)
outInfo.Remove(executive.TIME_STEPS())
for _t in _times:
    outInfo.Append(executive.TIME_STEPS(), float(_t))
outInfo.Remove(executive.TIME_RANGE())
outInfo.Append(executive.TIME_RANGE(), float(_times[0]))
outInfo.Append(executive.TIME_RANGE(), float(_times[-1]))
'''


def write_pv_source(stl_out, npz_out, pv_source_out):
    stl_abs = str(Path(stl_out).resolve())
    npz_abs = str(Path(npz_out).resolve())
    text = PV_SOURCE_TEMPLATE.format(stl_path=stl_abs, npz_path=npz_abs)
    text += '\n\n# --- Script (RequestInformation), paste separately ---\n\n'
    text += PV_SOURCE_INFO_TEMPLATE.format(npz_path=npz_abs)
    Path(pv_source_out).write_text(text)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('path_file', help='uav_path_file (rows: t x y z)')
    ap.add_argument('--input-parameters', default='input_parameters',
                     help="case's input_parameters, for uav_disk_radius/uav_tilt_* "
                          "and &NUMERICS dt (default: ./input_parameters)")
    ap.add_argument('--times-from', default=None,
                     help='sync frame times to this file: a probe _times.bin '
                          '(big-endian float64) or a plain text file with one time per line')
    ap.add_argument('--t0', type=float, default=None, help='first frame time (default: path start)')
    ap.add_argument('--t1', type=float, default=None, help='last frame time (default: path end)')
    ap.add_argument('--nframes', type=int, default=50, help='frame count for a uniform time grid (default 50)')
    ap.add_argument('--tilt', action='store_true',
                     help='enable auto-tilt orientation (requires uav_tilt_active=1 in &UAV to be meaningful)')
    ap.add_argument('--no-heading', dest='heading', action='store_false',
                     help='keep the body arms fixed in the X configuration instead of yawing to face travel direction')
    ap.add_argument('--stl-out', default='uav_drone.stl', help='output geometry path (written once)')
    ap.add_argument('--transforms-out', default='uav_transforms', help='output basename for .npz/.csv transform table')
    ap.add_argument('--pv-source-out', default='uav_drone_paraview_source.py',
                     help='output ParaView Programmable Source script')
    args = ap.parse_args()

    params = parse_input_parameters(args.input_parameters)
    disk_radius = params.get('uav_disk_radius', 0.15)
    tilt_active = bool(args.tilt and params.get('uav_tilt_active', 0))
    if args.tilt and not params.get('uav_tilt_active', 0):
        print(f"WARNING: --tilt given but uav_tilt_active is not set in "
              f"{args.input_parameters}; rendering flat (untilted).", file=sys.stderr)
    tau = float(params.get('uav_tilt_tau', 0.2))
    dt = float(params.get('dt', 5e-3))
    grav = float(params.get('grav', 9.81))

    scale = disk_radius / PROP_RADIUS_REF
    drone = build_drone_mesh(scale)
    drone.export(args.stl_out)
    print(f"Wrote geometry: {args.stl_out} "
          f"(scale={scale:.4g}, watertight={drone.is_watertight}, "
          f"vertices={len(drone.vertices)}, faces={len(drone.faces)})")

    path_t, path_x, path_y, path_z = read_uav_path(args.path_file)

    if args.times_from:
        times = read_times_file(args.times_from)
    else:
        t0 = args.t0 if args.t0 is not None else path_t[0]
        t1 = args.t1 if args.t1 is not None else path_t[-1]
        times = np.linspace(t0, t1, args.nframes)
    times = np.sort(np.asarray(times, dtype=float))

    mats = compute_transforms(times, path_t, path_x, path_y, path_z,
                               tilt_active, tau, dt, grav, args.heading)

    npz_out = Path(args.transforms_out).with_suffix('.npz')
    csv_out = Path(args.transforms_out).with_suffix('.csv')
    write_transforms(times, mats, npz_out, csv_out)
    print(f"Wrote {len(times)} transforms -> {npz_out}, {csv_out}")

    write_pv_source(args.stl_out, npz_out, args.pv_source_out)
    print(f"Wrote ParaView Programmable Source -> {args.pv_source_out}")
    print("In ParaView: Sources > Programmable Source (Output Type: vtkPolyData), "
          f"paste in {args.pv_source_out}'s two scripts (main + RequestInformation), Apply.")
