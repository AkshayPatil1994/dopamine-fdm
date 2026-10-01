"""
_core.py -- shared low-level helpers for dopamine_post.

Every other module in this package imports from here rather than owning a private
copy of namelist parsing, `<prefix>.<step>` directory globbing, big-endian binary
reading, ffmpeg movie encoding, or colour-scale clipping. This is the one place
those concerns live.
"""
import os
import re
import shutil
import struct
import subprocess
import tempfile
import warnings
from multiprocessing import get_context
from pathlib import Path

import numpy as np


# ── namelist parsing ─────────────────────────────────────────────────────────

def parse_input_parameters(fpath='input_parameters'):
    """Minimal Fortran-namelist parser: flat dict of lowercase key -> Python value.

    Supports integers, floats (including Fortran 'd'/'D' exponent notation),
    booleans (.TRUE./.FALSE.) and single-quoted strings. Comment lines (starting
    with !) are stripped before parsing.
    """
    text = Path(fpath).read_text()
    text = re.sub(r'!.*', '', text)
    text = text.replace('\n', ' ')

    params = {}
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


def read_input(case, keys):
    """Scalar values from `<case>/input_parameters`, filtered to `keys`.

    A thin convenience wrapper around parse_input_parameters for callers that
    only want a handful of scalars (e.g. Lx/Ly/Lz) and don't care about the
    rest of the namelist.
    """
    out = {}
    try:
        params = parse_input_parameters(Path(case) / 'input_parameters')
    except OSError:
        return out
    for k in keys:
        if k.lower() in params:
            out[k] = params[k.lower()]
    return out


# ── `<prefix><sep><step>` directory globbing ─────────────────────────────────

def list_indexed(directory, regex, prefix=None):
    """Sorted [(step, Path)] of files under `directory` matching `regex`.

    `regex` must have exactly two capture groups: (prefix, step). When `prefix`
    is None, every distinct group-1 value found is a candidate; the first
    (alphabetically) is used and a warning is issued if more than one exists.
    """
    found = {}
    for p in Path(directory).iterdir():
        m = regex.match(p.name)
        if m and (prefix is None or m.group(1) == prefix):
            found.setdefault(m.group(1), []).append((int(m.group(2)), p))
    if not found:
        raise FileNotFoundError(f"no files matching {regex.pattern!r} in {directory}")
    if prefix is None and len(found) > 1:
        warnings.warn(f"multiple prefixes {sorted(found)}; using '{sorted(found)[0]}'")
    return sorted(found[prefix or sorted(found)[0]])


# ── run.log step -> physical time lookup ─────────────────────────────────────

_LOG_LINE = re.compile(r"^\s*(\d+)\s+([0-9.]+E[+-]\d+)\s+(?:\S+\s+){5}\S+\s+[0-9.]+\s*$")


def times_from_log(case, steps):
    """Physical time at each of `steps`, interpolated from a run.log's monitor
    lines (`<case>/run.log`); None if the log is missing or too short to use."""
    log = Path(case) / 'run.log'
    if not log.exists():
        return None
    rows = [(int(m.group(1)), float(m.group(2)))
            for m in map(_LOG_LINE.match, open(log, errors='ignore')) if m]
    if len(rows) < 2:
        return None
    r = np.array(sorted(rows))
    return np.interp(steps, r[:, 0], r[:, 1])


# ── low-level big-endian Fortran-stream binary reader ────────────────────────

class BinaryReader:
    """Stateful reader for big-endian Fortran stream binary files.

    Pass `memmap=path` to back the float reads with a read-only memmap instead
    of copying into RAM (useful for the large field snapshots in a showcase
    animation that only ever touches one plane per frame).
    """

    def __init__(self, data: bytes = None, memmap=None):
        if (data is None) == (memmap is None):
            raise ValueError("pass exactly one of data or memmap")
        if memmap is not None:
            self._mm = np.memmap(memmap, dtype=np.uint8, mode='r')
            self._d = self._mm
        else:
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
        arr = np.frombuffer(self._d, dtype='>f8', count=n, offset=self._p)
        self._p += n * 8
        return arr.copy() if not isinstance(self._d, np.memmap) else arr

    def read_field(self):
        """One field block: 3-Int32 shape header + float64 data, Fortran order."""
        n1, n2, n3 = self.ri3()
        data = self.rd(n1 * n2 * n3).reshape((n1, n2, n3), order='F')
        return (n1, n2, n3), data


# ── ffmpeg movie encoding ─────────────────────────────────────────────────────

_FRAME_FN = None      # frame builder handed to forked workers (inherited, not pickled)


def _job(args):
    i, tmp, dpi = args
    fig = _FRAME_FN(i)
    fig.savefig(f"{tmp}/frame_{i:04d}.png", dpi=dpi, facecolor="white")
    import matplotlib.pyplot as plt
    plt.close(fig)


def render_movie(frame_fn, n, out, fps=8, width=900, workers=8, dpi=100):
    """Render frame_fn(i)->Figure for i<n in parallel and encode .gif / .mp4 with ffmpeg."""
    global _FRAME_FN
    _FRAME_FN = frame_fn
    tmp = tempfile.mkdtemp(prefix="dopamine_post_frames_")
    try:
        with get_context("fork").Pool(max(1, workers)) as pool:
            pool.map(_job, [(i, tmp, dpi) for i in range(n)])
        base = ["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(fps), "-i", f"{tmp}/frame_%04d.png"]
        out = str(out)
        if out.endswith(".mp4"):
            subprocess.run(base + ["-pix_fmt", "yuv420p", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", out],
                           check=True)
        else:
            vf = f"fps={fps},scale={width}:-1:flags=lanczos"
            subprocess.run(base + ["-vf", vf + ",palettegen=max_colors=128", f"{tmp}/pal.png"], check=True)
            subprocess.run(base + ["-i", f"{tmp}/pal.png", "-lavfi",
                                   vf + "[x];[x][1:v]paletteuse=dither=bayer:bayer_scale=5",
                                   "-loop", "0", out], check=True)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print(f"Wrote {out}  ({n} frames)")


# ── colour-scale clipping ─────────────────────────────────────────────────────

def percentile_clim(values, lo=1, hi=99):
    """(vmin, vmax) percentile-clipped colour limits for a scatter/imshow colour
    array, ignoring NaNs; guarantees vmax > vmin."""
    vmin, vmax = np.nanpercentile(values, [lo, hi])
    if vmax <= vmin:
        vmax = vmin + 1.0
    return float(vmin), float(vmax)


# ── ParaView .pvd collection writer ───────────────────────────────────────────

def write_pvd(path, entries):
    """Write a ParaView .pvd time-series collection.

    `entries` is an iterable of (time, file_path) pairs; `file_path` is written
    relative to `path`'s own directory.
    """
    path = Path(path)
    out_dir = path.resolve().parent
    lines = ['<?xml version="1.0"?>',
             '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">',
             '  <Collection>']
    for t, f in entries:
        # Resolve both sides before computing the relative path: a caller-supplied
        # relative `f` is relative to the CURRENT WORKING DIRECTORY, not necessarily
        # to `path`'s own directory, so comparing them as given can double up a
        # shared subdirectory (e.g. out="ibm_surface/surface.pvd" with
        # f="ibm_surface/surface.100.vtp" -> wrongly "ibm_surface/surface.100.vtp"
        # instead of "surface.100.vtp", which ParaView resolves one level too deep).
        rel = os.path.relpath(Path(f).resolve(), out_dir)
        lines.append(f'    <DataSet timestep="{t:.10g}" group="" part="0" file="{rel}"/>')
    lines += ['  </Collection>', '</VTKFile>']
    path.write_text('\n'.join(lines) + '\n')
    print(f"Wrote {path}")


# ── cubic Hermite (Catmull-Rom tangent) interpolation ─────────────────────────

def _hermite_tangents(times, points):
    """Catmull-Rom tangents (dpoint/dt) at each knot: centred difference at
    interior knots, one-sided difference at the two end knots."""
    tang = np.empty_like(points)
    tang[1:-1] = (points[2:] - points[:-2]) / (times[2:] - times[:-2])[:, None]
    tang[0] = (points[1] - points[0]) / (times[1] - times[0])
    tang[-1] = (points[-1] - points[-2]) / (times[-1] - times[-2])
    return tang


def hermite_eval(t, times, points):
    """Cubic-Hermite (Catmull-Rom tangents) interpolation of `points` (N,3) at
    the knot times `times` (N,), evaluated at scalar or array `t`. Tangents at
    interior knots are the centred difference of neighbouring points; end
    knots use a one-sided difference. `t` outside [times[0], times[-1]] is
    clamped to the nearest endpoint value (flat extrapolation), not
    extrapolated by the cubic. Matches uav_actuator.f90's uav_current_center
    construction exactly, so a rendered path sits where the solver's own
    actuator-disk force was applied.
    """
    times = np.asarray(times, float)
    points = np.asarray(points, float)
    n = len(times)
    t = np.clip(np.atleast_1d(np.asarray(t, float)), times[0], times[-1])

    tang = _hermite_tangents(times, points)

    seg = np.clip(np.searchsorted(times, t, side='right') - 1, 0, n - 2)
    t0, t1 = times[seg], times[seg + 1]
    h = (t1 - t0)[:, None]
    s = ((t - t0) / (t1 - t0))[:, None]

    p0, p1 = points[seg], points[seg + 1]
    m0, m1 = tang[seg] * h, tang[seg + 1] * h

    s2, s3 = s * s, s * s * s
    h00 = 2 * s3 - 3 * s2 + 1
    h10 = s3 - 2 * s2 + s
    h01 = -2 * s3 + 3 * s2
    h11 = s3 - s2
    return h00 * p0 + h10 * m0 + h01 * p1 + h11 * m1


def hermite_accel(t, times, points):
    """Exact second derivative of the same cubic-Hermite segment used by
    hermite_eval (analytic, not a finite difference) -- mirrors
    uav_actuator.f90's hermite_accel. Zero outside [times[0], times[-1]]
    (position is clamped flat there, so acceleration is zero, not just the
    boundary value)."""
    times = np.asarray(times, float)
    points = np.asarray(points, float)
    n = len(times)
    t = np.atleast_1d(np.asarray(t, float))
    out = np.zeros((len(t),) + points.shape[1:])
    inside = (t > times[0]) & (t < times[-1])
    if not np.any(inside):
        return out

    tang = _hermite_tangents(times, points)

    seg = np.clip(np.searchsorted(times, t[inside], side='right') - 1, 0, n - 2)
    t0, t1 = times[seg], times[seg + 1]
    h = (t1 - t0)[:, None]
    s = ((t[inside] - t0) / (t1 - t0))[:, None]

    p0, p1 = points[seg], points[seg + 1]
    m0, m1 = tang[seg], tang[seg + 1]

    out[inside] = ((12 * s - 6) * p0 + (6 * s - 4) * h * m0
                   + (-12 * s + 6) * p1 + (6 * s - 2) * h * m1) / (h * h)
    return out


# ── probe_output.f90 meta.txt parsing ─────────────────────────────────────────

def read_meta(path):
    """Parse a probe `<base>_meta.txt` key=value file (line/slice probes,
    written by probe_output.f90) into a plain dict. Integer-valued keys
    (ncomp, npts, n1, n2, nsnaps) and the float `pos` key are coerced;
    everything else (dir, comps, times, ...) stays a stripped string.
    """
    meta = {}
    for line in Path(path).read_text().splitlines():
        if '=' not in line:
            continue
        k, _, v = line.partition('=')
        k, v = k.strip(), v.strip()
        if k in ('ncomp', 'npts', 'n1', 'n2', 'nsnaps'):
            meta[k] = int(v)
        elif k == 'pos':
            meta[k] = float(v)
        else:
            meta[k] = v
    if 'npts' not in meta and 'n1' not in meta:
        raise ValueError(f'{path}: not a recognised probe meta file (no npts/n1 key)')
    return meta


def comp_names(s):
    """Ordered component name list from a meta file's `comps` string.

    Order must match probe_output.f90's parse_comps/write loop, which always
    emits components in fixed U,V,W,P,T,C order (not the order the user typed
    them in the &... comps string) regardless of which subset is requested.
    """
    upper = s.strip().upper()
    names = [c for c in ('U', 'V', 'W', 'P', 'T', 'C') if c in upper]
    return names or ['U', 'V', 'W']
