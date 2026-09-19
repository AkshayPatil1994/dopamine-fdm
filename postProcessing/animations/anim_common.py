"""Shared helpers for the showcase animation scripts (fast plane reads, time map, GIF output)."""
import re
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.family": "sans-serif", "font.size": 11,
                     "figure.facecolor": "#0d1117", "savefig.facecolor": "#0d1117",
                     "text.color": "#e6edf3", "axes.labelcolor": "#e6edf3",
                     "xtick.color": "#e6edf3", "ytick.color": "#e6edf3",
                     "axes.edgecolor": "#30363d"})


class Snapshot:
    """Memory-mapped view of one fdm-dopamine snapshot: reads single planes without loading the file.

    Blocks after the header are U, V, W, P, then C (sediment), nu_t (sgs) and T (boussinesq) if enabled.
    """

    def __init__(self, path, sgs_model=0, sediment_flag=0, boussinesq_flag=0):
        self.path = str(path)
        pos = 0
        counts = []
        with open(self.path, "rb") as f:
            for _ in range(6):
                f.seek(pos)
                n = int(np.frombuffer(f.read(4), ">i4")[0])
                counts.append(n)
                pos += 4 + 8 * n
            self.nx, self.ny, self.nz, self.nxm, self.nym, self.nzm = counts
            names = ["U", "V", "W", "P"]
            if sediment_flag >= 1: names.append("C")
            if sgs_model != 0: names.append("nu_t")
            if boussinesq_flag >= 1: names.append("T")
            self.f = {}
            for nm in names:
                f.seek(pos)
                shp = tuple(int(v) for v in np.frombuffer(f.read(12), ">i4"))
                self.f[nm] = np.memmap(self.path, dtype=">f8", mode="r", offset=pos + 12,
                                       shape=shp, order="F")
                pos += 12 + 8 * shp[0] * shp[1] * shp[2]

    def coords(self):
        with open(self.path, "rb") as f:
            pos, out = 0, []
            for _ in range(6):
                f.seek(pos)
                n = int(np.frombuffer(f.read(4), ">i4")[0])
                out.append(np.frombuffer(f.read(8 * n), ">f8").astype(float))
                pos += 4 + 8 * n
        return out  # x, y, z, xm, ym, zm

    def plane(self, name, axis, i):
        """Cell-centred 2-D plane of `name` at cell index i (0-based, ghosts stripped) normal to `axis`."""
        stag = {"U": 0, "V": 1, "W": 2}.get(name)
        sl = []
        for d in range(3):
            if d == axis:
                sl.append(slice(i, i + 2) if stag == d else slice(i + 1, i + 2))
            else:
                sl.append(slice(None) if stag == d else slice(1, -1))
        v = np.asarray(self.f[name][tuple(sl)], dtype=float)
        for d in range(3):
            if stag == d:
                n = v.shape[d]
                v = 0.5 * (np.take(v, range(0, n - 1), axis=d) + np.take(v, range(1, n), axis=d))
        return v.squeeze(axis=axis)


def step_time_map(logs):
    """(step, t) rows from solver monitor lines in the given run logs, for interpolating snapshot times."""
    rows = []
    pat = re.compile(r"^\s*(\d+)\s+([0-9.]+E[+-]\d+)\s+(?:\S+\s+){5}\S+\s+[0-9.]+\s*$")
    for lg in logs:
        for line in open(lg, errors="ignore"):
            m = pat.match(line)
            if m:
                rows.append((int(m.group(1)), float(m.group(2))))
    rows.sort()
    a = np.array(rows)
    return a[:, 0], a[:, 1]


def list_steps(fields_dir, prefix):
    out = []
    for p in Path(fields_dir).iterdir():
        m = re.fullmatch(re.escape(prefix) + r"\.(\d+)", p.name)
        if m:
            out.append((int(m.group(1)), p))
    return sorted(out)


def save_gif(frames_dir, out, fps, width):
    """Encode PNG frames (frame_%04d.png) into an optimised GIF via ffmpeg palette."""
    out = str(out)
    pal = str(Path(frames_dir) / "pal.png")
    vf = f"fps={fps},scale={width}:-1:flags=lanczos"
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(fps), "-i",
                    f"{frames_dir}/frame_%04d.png", "-vf", vf + ",palettegen=max_colors=96:stats_mode=diff", pal], check=True)
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(fps), "-i",
                    f"{frames_dir}/frame_%04d.png", "-i", pal, "-lavfi",
                    vf + "[x];[x][1:v]paletteuse=dither=bayer:bayer_scale=5", "-loop", "0", out], check=True)


def render(frame_fn, n, out, fps=12, width=720, dpi=110, workers=8):
    """Call frame_fn(i) -> Figure for each frame in parallel, then encode to GIF/MP4 at `out`."""
    from concurrent.futures import ProcessPoolExecutor
    tmp = tempfile.mkdtemp(prefix="anim_")
    def job(i):
        fig = frame_fn(i)
        fig.savefig(f"{tmp}/frame_{i:04d}.png", dpi=dpi)
        plt.close(fig)
    from multiprocessing import get_context
    global _JOB
    _JOB = job
    with get_context("fork").Pool(workers) as p:
        p.map(_run, range(n))
    if str(out).endswith(".mp4"):
        subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(fps), "-i", f"{tmp}/frame_%04d.png",
                        "-pix_fmt", "yuv420p", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", str(out)], check=True)
    else:
        save_gif(tmp, out, fps, width)


def _run(i):
    _JOB(i)
