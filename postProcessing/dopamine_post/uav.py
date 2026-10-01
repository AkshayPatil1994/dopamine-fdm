"""
uav.py -- UAV actuator-disk path and drone geometry for fdm-dopamine cases.

Reads a `uav_path_file` (rows: t x y z, see src/uav_actuator.f90 / docs/
Input-Parameters.md's &UAV group) and renders it two ways for ParaView:

  UAVPath.disk_animation(...)   a moving flat disc (.vtp + .pvd), one per
                                 requested time -- what generate_UAVpath.py
                                 used to produce.
  UAVPath.drone_animation(...)  a single scaled quadcopter-body mesh (.stl,
                                 written once) plus a per-time rigid-body
                                 transform table (.npz/.csv) and a ParaView
                                 "Programmable Source" script that re-poses
                                 it in memory -- what generate_UAVdrone.py
                                 used to produce.

Both interpolate the path with `_core.hermite_eval` (cubic Hermite,
Catmull-Rom tangents, clamped outside the waypoint file's time range),
bit-for-bit the same construction as uav_actuator.f90's uav_current_center,
so the rendered disk/body sits exactly where the solver's own actuator-disk
force was applied -- not merely close to it. `drone_animation`'s optional
auto-tilt orientation additionally replays uav_disk_state's recursive
low-pass filter (uav_actuator.f90) using `_core.hermite_accel`'s second
derivative of the same Hermite segments; per that function's own module
docstring, the solver calls this filter once per RK sub-stage rather than
once per dt, so the replay is a close but not bit-exact reproduction of the
in-run tilt -- good enough to visualize, not to re-derive the applied force.

Library use
-----------
    from dopamine_post.uav import UAVPath

    path = UAVPath.read("uav_path_file.dat")
    path.disk_animation(input_parameters="input_parameters",
                         times_from="slices/y015_times.bin",
                         out_dir="uav_path", pvd_out="uav_path.pvd")
    path.drone_animation(input_parameters="input_parameters", tilt=True,
                          nframes=100)

Conventions
-----------
x streamwise, y wall-normal ("up"), z spanwise -- same as the rest of
dopamine_post. Phase 1/2 scope note (see docs/UAV_ActuatorDisk_Design.md):
the solver's disk stays horizontal (normal = +y) unless uav_tilt_active=1;
`disk_animation` always draws a flat, horizontal disk (it has no tilt
model of its own -- only `drone_animation` optionally tilts, for the body).

Requires: numpy; `pyvista` for disk_animation; `trimesh` (optional
`manifold3d` for a proper boolean union) for drone_animation.
"""
import sys
from pathlib import Path

import numpy as np

from . import _core

# ---------------------------------------------------------------
# Reference quadcopter dimensions [m] (box/cylinder model), converted from
# its native mm to m; scaled uniformly so PROP_RADIUS_REF maps onto a
# case's uav_disk_radius.
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


def read_times_file(fpath):
    """Load output times: a numpy .bin (big-endian float64, e.g. a probe's
    own <fileout>_times.bin, see docs/Input-Parameters.md &STATISTICS) if
    the extension is .bin, else a plain text file with one time per line.
    """
    fpath = Path(fpath)
    if fpath.suffix == '.bin':
        return np.fromfile(fpath, dtype='>f8')
    return np.loadtxt(fpath).reshape(-1)


# ════════════════════════════════════════════════════════════════════════════
# drone body: geometry
# ════════════════════════════════════════════════════════════════════════════
def build_drone_mesh(scale):
    """Build the watertight single-propeller+legs mesh at the given uniform
    scale (dimensionless multiplier on the *_REF dimensions above), already
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


# ════════════════════════════════════════════════════════════════════════════
# drone body: orientation -- reproduces uav_disk_state's tilt low-pass filter
# and uav_disk_frame (uav_actuator.f90) in Python.
# ════════════════════════════════════════════════════════════════════════════
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


def tilt_frames(sample_t, path_t, path_points, tau, grav):
    """Sequentially replay uav_disk_state's tilt low-pass filter over the
    (already sorted, strictly increasing) times in sample_t. Returns
    nvec/e1/e2 arrays, one row per sample_t entry.
    """
    n = len(sample_t)
    nvecs = np.empty((n, 3))
    e1s = np.empty((n, 3))
    e2s = np.empty((n, 3))

    accel = _core.hermite_accel(sample_t, path_t, path_points)
    ax, ay, az = accel[:, 0], accel[:, 1], accel[:, 2]

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


def heading_yaw(times, path_t, path_points, eps):
    """Horizontal heading angle (about +y) from the path's own velocity
    direction, via a central-difference derivative of _core.hermite_eval.
    Visualization-only -- the actuator-disk model has no heading.
    """
    fwd = _core.hermite_eval(times + eps, path_t, path_points)
    bwd = _core.hermite_eval(times - eps, path_t, path_points)
    vx, vz = (fwd[:, 0] - bwd[:, 0]) / (2 * eps), (fwd[:, 2] - bwd[:, 2]) / (2 * eps)
    speed = np.hypot(vx, vz)
    psi = np.where(speed > 1e-9, np.arctan2(vz, vx), 0.0)
    return psi


def compute_transforms(times, path_t, path_points, tilt_active, tau, dt, grav, heading):
    """Per-time 4x4 rigid-body transforms (solver axis convention)."""
    centers = _core.hermite_eval(times, path_t, path_points)
    xc, yc, zc = centers[:, 0], centers[:, 1], centers[:, 2]

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
        nvecs, e1s, e2s = tilt_frames(grid, path_t, path_points, tau, grav)
        sel = np.searchsorted(grid, times)
        nvecs, e1s, e2s = nvecs[sel], e1s[sel], e2s[sel]
    else:
        nvecs = np.tile([0.0, 1.0, 0.0], (n, 1))
        e1s = np.tile([1.0, 0.0, 0.0], (n, 1))
        e2s = np.tile([0.0, 0.0, 1.0], (n, 1))

    if heading:
        span = max(path_t[-1] - path_t[0], 1e-6)
        psi = heading_yaw(times, path_t, path_points, eps=1e-5 * span)
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
# Auto-generated by dopamine_post.uav.UAVPath.drone_animation -- paste this
# whole file into a ParaView "Programmable Source" (Sources > Programmable
# Source, Output Type: vtkPolyData), then Apply. The single STL below is
# read once and re-posed in memory at every animation time step from the
# transform table -- no per-frame geometry is ever written to disk.

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


# ════════════════════════════════════════════════════════════════════════════
# UAVPath
# ════════════════════════════════════════════════════════════════════════════
class UAVPath:
    """A uav_path_file: waypoint times `t` (N,) and centres `points` (N,3).

    `.center(times)` / `.accel(times)` evaluate the same cubic-Hermite
    construction the solver uses (see module docstring); `.disk_animation`
    and `.drone_animation` render the path for ParaView.
    """

    def __init__(self, t, points, path_file=None):
        self.t = np.asarray(t, float)
        self.points = np.asarray(points, float)
        self.path_file = Path(path_file) if path_file else None

    @classmethod
    def read(cls, path):
        """Read a uav_path_file: rows 't x y z', blank/'#' lines skipped."""
        rows = []
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                rows.append([float(v) for v in line.split()[:4]])
        if len(rows) < 2:
            raise ValueError(f'{path}: fewer than 2 waypoints found')
        arr = np.array(rows)
        t, points = arr[:, 0], arr[:, 1:4]
        if np.any(np.diff(t) <= 0):
            raise ValueError(f'{path}: t column must be strictly increasing')
        return cls(t, points, path_file=path)

    def __repr__(self):
        return (f"UAVPath({len(self.t)} waypoints, t=[{self.t[0]:.4g},{self.t[-1]:.4g}], "
                f"file={self.path_file})")

    def center(self, times):
        """Disk/body centre (M,3) at `times`: uav_actuator.f90's
        uav_current_center, exactly (see _core.hermite_eval)."""
        return _core.hermite_eval(times, self.t, self.points)

    def accel(self, times):
        """Second derivative (M,3) of the same Hermite path at `times`
        (see _core.hermite_accel); zero outside the waypoint time range."""
        return _core.hermite_accel(times, self.t, self.points)

    def _frame_times(self, times, times_from, t0, t1, nframes):
        if times is not None:
            return np.asarray(times, float)
        if times_from is not None:
            return read_times_file(times_from)
        t0 = t0 if t0 is not None else self.t[0]
        t1 = t1 if t1 is not None else self.t[-1]
        return np.linspace(t0, t1, nframes)

    # ════════════════════════════════════════════════════════════════════
    # moving-disk animation (generate_UAVpath.py)
    # ════════════════════════════════════════════════════════════════════
    def disk_animation(self, times=None, times_from=None, t0=None, t1=None, nframes=50,
                        radius=None, n_theta=None, input_parameters=None,
                        out_dir="uav_path", pvd_out="uav_path.pvd"):
        """Render the path as a moving flat disc: one .vtp per frame time
        plus a .pvd collection (`_core.write_pvd`) so the sequence loads in
        ParaView as a single animated object, synced to the shared time
        toolbar.

        times          explicit frame times (overrides times_from/t0/t1/nframes)
        times_from     sync frames exactly to an existing output's times: a
                       probe's <fileout>_times.bin (big-endian float64) or a
                       plain text file with one time per line
        t0, t1, nframes  uniform time grid spanning the path's own range by default
        radius, n_theta  disk geometry (default 0.15, 48; or read from
                       input_parameters's uav_disk_radius/uav_n_theta if given
                       and not overridden here)
        input_parameters  case namelist to pull uav_disk_radius/uav_n_theta from
        out_dir, pvd_out  where the per-frame .vtp files and the .pvd collection go

        Phase 1/2 scope: the disk stays horizontal (normal=+y) for the whole
        path -- matches the solver's flat actuator-disk model (see module
        docstring). Requires pyvista.
        """
        import pyvista as pv

        if input_parameters is not None and (radius is None or n_theta is None):
            params = _core.parse_input_parameters(input_parameters)
            if radius is None:
                radius = float(params.get('uav_disk_radius', 0.15))
            if n_theta is None:
                n_theta = int(params.get('uav_n_theta', 48))
        radius = 0.15 if radius is None else radius
        n_theta = 48 if n_theta is None else n_theta

        times = self._frame_times(times, times_from, t0, t1, nframes)
        centers = self.center(times)

        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        pvd_path = Path(pvd_out)

        entries = []
        for i, (t, c) in enumerate(zip(times, centers)):
            x, y, z = c
            # Phase 1/2: disk stays horizontal, normal=+y (see module docstring)
            disc = pv.Disc(center=(x, y, z), inner=0.0, outer=radius,
                            normal=(0.0, 1.0, 0.0), r_res=1, c_res=n_theta)
            disc.field_data['time'] = np.array([t])
            disc.field_data['centre'] = np.array([[x, y, z]])
            vtp_path = out_dir / f'uav_disk.{i:05d}.vtp'
            disc.save(str(vtp_path))
            entries.append((float(t), vtp_path.resolve()))

        _core.write_pvd(pvd_path, entries)
        print(f'{len(entries)} frames -> {pvd_path}')
        print('Open in ParaView: File > Open > ' + str(pvd_path))
        return pvd_path

    # ════════════════════════════════════════════════════════════════════
    # drone body animation (generate_UAVdrone.py)
    # ════════════════════════════════════════════════════════════════════
    def drone_animation(self, input_parameters="input_parameters", times=None, times_from=None,
                         t0=None, t1=None, nframes=50, tilt=False, heading=True,
                         stl_out="uav_drone.stl", transforms_out="uav_transforms",
                         pv_source_out="uav_drone_paraview_source.py"):
        """Build a scaled quadcopter body (single STL, written once) plus a
        per-time rigid-body transform table (.npz + a human-readable .csv),
        and a ready-to-paste ParaView "Programmable Source" script that reads
        the STL once and re-poses it in memory at the current animation time.

        input_parameters  case namelist: uav_disk_radius (sets the mesh scale),
                           uav_tilt_active/uav_tilt_tau, &NUMERICS dt, grav
        times / times_from / t0, t1, nframes  frame times, as in disk_animation
        tilt       replay the auto-tilt filter (only meaningful with
                   uav_tilt_active=1 in input_parameters; otherwise a warning
                   is printed and the body stays level)
        heading    yaw the body to face the path's own horizontal velocity
                   direction (visualization-only -- the actuator disk has no
                   heading); False keeps the arms fixed in the X configuration
        stl_out, transforms_out, pv_source_out  output paths

        Requires trimesh (optional manifold3d for a proper boolean union).
        """
        params = _core.parse_input_parameters(input_parameters)
        disk_radius = float(params.get('uav_disk_radius', 0.15))
        tilt_active = bool(tilt and params.get('uav_tilt_active', 0))
        if tilt and not params.get('uav_tilt_active', 0):
            print(f"WARNING: tilt=True given but uav_tilt_active is not set in "
                  f"{input_parameters}; rendering flat (untilted).", file=sys.stderr)
        tau = float(params.get('uav_tilt_tau', 0.2))
        dt = float(params.get('dt', 5e-3))
        grav = float(params.get('grav', 9.81))

        scale = disk_radius / PROP_RADIUS_REF
        drone = build_drone_mesh(scale)
        drone.export(stl_out)
        print(f"Wrote geometry: {stl_out} "
              f"(scale={scale:.4g}, watertight={drone.is_watertight}, "
              f"vertices={len(drone.vertices)}, faces={len(drone.faces)})")

        times = self._frame_times(times, times_from, t0, t1, nframes)
        times = np.sort(np.asarray(times, dtype=float))

        mats = compute_transforms(times, self.t, self.points, tilt_active, tau, dt, grav, heading)

        npz_out = Path(transforms_out).with_suffix('.npz')
        csv_out = Path(transforms_out).with_suffix('.csv')
        write_transforms(times, mats, npz_out, csv_out)
        print(f"Wrote {len(times)} transforms -> {npz_out}, {csv_out}")

        write_pv_source(stl_out, npz_out, pv_source_out)
        print(f"Wrote ParaView Programmable Source -> {pv_source_out}")
        print("In ParaView: Sources > Programmable Source (Output Type: vtkPolyData), "
              f"paste in {pv_source_out}'s two scripts (main + RequestInformation), Apply.")
        return Path(stl_out), npz_out, csv_out, Path(pv_source_out)
