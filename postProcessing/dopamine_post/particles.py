"""
particles.py -- post-processing of the fdm-dopamine Lagrangian point particles.

One module, three layers (use the library from your own scripts, or the CLI):

  1. DATA      ParticleData   all snapshots '<prefix>_particles.<step>' as id-indexed
                              arrays (Nt, Np); NaN where a particle is absent.
  2. ANALYSIS  ParticleData   pdf / pdf2d / moments / msd / track / pick_ids / mask
  3. FIGURES   Geometry       dune wall (analytic) or sdf = 0 iso-surface (any solid)
               plot_pdfs, plot_tracks, animate_tracks, animate_cloud

Library use
-----------
    from dopamine_post.particles import ParticleData, Geometry, plot_tracks, animate_tracks

    pdat = ParticleData.from_case(".")                   # times=... if dt is adaptive
    x, p = pdat.pdf("y", frames=slice(-20, None), bins=80)            # location PDF
    x, p = pdat.pdf("ax", bins=100, log=True)                         # accel. PDF
    ids  = pdat.pick_ids(6, seed=1, where={"y": (0.1, 0.3)})
    geo  = Geometry.from_sdf(".", "SDF_in", grid_file="fields/grid.out")
    animate_tracks(pdat, ids, "tracks3d.gif", frames=slice(-40, None), geometry=geo, view="3d")

Command line  (python3 -m dopamine_post particles <command> -h for the options)
------------
    dopamine-post particles pdf y u v --last 20 --bins 80 --out pdf.png
    dopamine-post particles tracks   --last 20 --n 8 --out tracks.png
    dopamine-post particles animate  --last 40 --n 10 --view 3d --sdf-file SDF_in \\
                                  --grid-file ../waveWall/fields/grid.out --out tracks3d.gif
    dopamine-post particles cloud    --last 40 --color speed --out cloud.gif

Conventions
-----------
x streamwise, y wall-normal (drawn vertical in 3-D), z spanwise.  "last N" always
means the last N *snapshots* (not solver steps).  Quantities for pdf/get():
x y z  u v w  speed  ax ay az  amag  age.  ax/ay/az are the solver-computed
per-step acceleration (src/particles.f90) for snapshots written after that
feature was added (10-column binary layout); for older 7-column snapshots they
fall back to a finite difference of the stored velocity between snapshots,
only as good as the snapshot spacing (see ParticleData.acc).

Requires: numpy, matplotlib, ffmpeg (movies); scikit-image only for --sdf-file
(via dopamine_post.sdf.Geometry.from_sdf).
"""
import os
import re
import struct
import warnings
from pathlib import Path

import numpy as np

from . import _core
from ._core import list_indexed, read_input, render_movie

_PAT = re.compile(r"^([A-Za-z0-9_]+)_particles\.(\d+)$")
_AXES = {"x": 0, "y": 1, "z": 2}
_VEL = {"u": 0, "v": 1, "w": 2}
_ACC = {"ax": 0, "ay": 1, "az": 2}
SAND = (0.85, 0.74, 0.54)


# ════════════════════════════════════════════════════════════════════════════
# 1. DATA
# ════════════════════════════════════════════════════════════════════════════
def read_snapshot(path):
    """One snapshot -> (id[n], data[n,10]); columns x,y,z,u,v,w,ax,ay,az,age
    (big-endian). Detects legacy 7-column snapshots (x,y,z,u,v,w,age, written
    before src/particles.f90 stored per-step acceleration) from the record
    size and pads their ax/ay/az columns with NaN -- ParticleData.acc then
    falls back to a finite-difference reconstruction for those snapshots only.
    """
    with open(path, "rb") as f:
        n = struct.unpack(">i", f.read(4))[0]
        pid = np.frombuffer(f.read(4 * n), ">i4").astype(np.int64)
        rest = f.read()
    if n == 0:
        return pid, np.zeros((0, 10))
    ncol, remainder = divmod(len(rest), 8 * n)
    if remainder != 0 or ncol not in (7, 10):
        raise ValueError(f"{path}: unexpected particle record size ({len(rest)} bytes, n={n})")
    d = np.frombuffer(rest, ">f8", count=n * ncol).reshape(n, ncol).astype(float)
    if ncol == 7:
        out = np.full((n, 10), np.nan)
        out[:, 0:6] = d[:, 0:6]
        out[:, 9] = d[:, 6]
        d = out
    return pid, d


def list_snapshots(fields_dir, prefix=None):
    """Sorted [(step, path)] of one prefix (auto-detected if None)."""
    try:
        return list_indexed(fields_dir, _PAT, prefix)
    except FileNotFoundError:
        raise FileNotFoundError(f"no '<prefix>_particles.STEP' files in {fields_dir}")


def _times_from_log(case, steps):
    return _core.times_from_log(case, steps)


class ParticleData:
    """Id-indexed particle time series.

    steps (Nt,)  solver step of each snapshot     t (Nt,)  physical time
    ids   (Np,)  sorted particle ids              L  (Lx,Ly,Lz)   periodic (x,y,z flags)
    pos, vel (Nt,Np,3)  x,y,z and u,v,w           age (Nt,Np)
    """

    def __init__(self, steps, t, ids, pos, vel, age, L, periodic=(True, False, True),
                 acc_raw=None, case_dir=".", prefix=None):
        self.steps, self.t, self.ids = steps, t, ids
        self.pos, self.vel, self.age = pos, vel, age
        self.L, self.periodic = tuple(L), tuple(periodic)
        self.case_dir, self.prefix = case_dir, prefix
        self._acc_raw = acc_raw if acc_raw is not None else np.full_like(vel, np.nan)
        self._acc = None

    @classmethod
    def from_case(cls, case_dir=".", prefix=None, every=1, steps=None, times=None,
                  dt=None, periodic=(True, False, True), L=None):
        """Load all snapshots of a case.

        every  use every Nth snapshot           steps  (min, max) step window
        times  physical time per *selected* snapshot (else run.log, else step*dt)
        dt     fallback uniform dt (input_parameters dt if None)
        L      (Lx,Ly,Lz), read from input_parameters if None
        """
        case = Path(case_dir)
        snaps = list_snapshots(case / "fields", prefix)[::every]
        if steps is not None:
            snaps = [s for s in snaps if steps[0] <= s[0] <= steps[1]]
        st = np.array([s for s, _ in snaps])

        if L is None:
            v = read_input(case, ("Lx", "Ly", "Lz"))
            if len(v) < 3:
                raise ValueError("Lx/Ly/Lz not found in input_parameters; pass L=(Lx,Ly,Lz)")
            L = (v["Lx"], v["Ly"], v["Lz"])

        if times is not None:
            t = np.asarray(times, float)
        else:
            t = _times_from_log(case, st)
            if t is None:
                dt = dt or read_input(case, ("dt",)).get("dt")
                if dt is None:
                    raise ValueError("no run.log and no dt; pass times= or dt=")
                warnings.warn(f"no run.log: using t = step*dt with dt={dt}. Wrong if the run "
                              "used adaptive dt (cfl_adaptive=1); pass times= instead.")
                t = st * dt

        raw = [read_snapshot(p) for _, p in snaps]
        ids = np.unique(np.concatenate([r[0] for r in raw])) if raw else np.array([], dtype=np.int64)
        pos = np.full((len(raw), len(ids), 3), np.nan)
        vel = np.full_like(pos, np.nan)
        acc_raw = np.full_like(pos, np.nan)
        age = np.full((len(raw), len(ids)), np.nan)
        for k, (pid, d) in enumerate(raw):
            j = np.searchsorted(ids, pid)
            pos[k, j], vel[k, j], acc_raw[k, j], age[k, j] = d[:, 0:3], d[:, 3:6], d[:, 6:9], d[:, 9]
        prefix_used = prefix
        if prefix_used is None and snaps:
            m = _PAT.match(Path(snaps[0][1]).name)
            prefix_used = m.group(1) if m else _detect_prefix(case / "fields")
        return cls(st, t, ids, pos, vel, age, L, periodic, acc_raw, case_dir=str(case),
                   prefix=prefix_used)

    def __repr__(self):
        return (f"ParticleData({len(self.t)} snapshots, {len(self.ids)} ids, "
                f"t=[{self.t[0]:.4g},{self.t[-1]:.4g}], L={self.L})")

    # -- derived fields --
    @property
    def acc(self):
        """Acceleration (Nt,Np,3): the solver-stored value (src/particles.f90) when
        present, else a finite difference of the stored velocity on the snapshot
        times (only as good as the snapshot spacing) for legacy 7-column snapshots."""
        if self._acc is None:
            grad = np.gradient(self.vel, self.t, axis=0)
            self._acc = np.where(np.isnan(self._acc_raw), grad, self._acc_raw)
        return self._acc

    def unwrapped(self):
        """Positions with periodic jumps removed (continuous trajectories)."""
        p = self.pos.copy()
        d = np.diff(p, axis=0)
        for c in range(3):
            if self.periodic[c]:
                jump = -self.L[c] * np.round(d[..., c] / self.L[c])
                p[1:, :, c] += np.cumsum(np.nan_to_num(jump), axis=0)
        return p

    def get(self, name):
        """Quantity as an (Nt,Np) array."""
        if name in _AXES:
            return self.pos[..., _AXES[name]]
        if name in _VEL:
            return self.vel[..., _VEL[name]]
        if name in _ACC:
            return self.acc[..., _ACC[name]]
        if name == "speed":
            return np.linalg.norm(self.vel, axis=-1)
        if name == "amag":
            return np.linalg.norm(self.acc, axis=-1)
        if name == "age":
            return self.age
        raise KeyError(f"unknown quantity '{name}'")

    # ════════════════════════════════════════════════════════════════════════
    # 2. ANALYSIS
    # ════════════════════════════════════════════════════════════════════════
    def mask(self, frames=None, ids=None, where=None):
        """Boolean (Nt,Np) selector.  frames: slice/list of snapshots; ids: particle
        ids; where: {quantity: (lo, hi)}, e.g. {"y": (0.1, 0.2)}."""
        rows = np.zeros(len(self.t), bool)
        rows[np.arange(len(self.t))[frames] if frames is not None else slice(None)] = True
        cols = np.isin(self.ids, ids) if ids is not None else np.ones(len(self.ids), bool)
        m = rows[:, None] & cols[None, :]
        for q, (lo, hi) in (where or {}).items():
            v = self.get(q)
            m &= (v >= lo) & (v <= hi)
        return m

    def pick_ids(self, n, seed=0, frame=0, where=None):
        """n random ids present at snapshot `frame` (optionally satisfying `where`)."""
        ok = self.mask(frames=[frame], where=where)[frame] & ~np.isnan(self.pos[frame, :, 0])
        cand = self.ids[ok]
        if len(cand) < n:
            raise ValueError(f"only {len(cand)} candidate ids at frame {frame}")
        return np.sort(np.random.default_rng(seed).choice(cand, n, replace=False))

    def pdf(self, quantity, frames=None, ids=None, where=None, bins=64, range=None,
            density=True, log=False):
        """PDF (histogram if density=False) of a quantity pooled over the selected
        snapshots/particles.  Returns (bin_centres, values); log=True turns empty
        bins into NaN for semilog plots."""
        v = self.get(quantity)[self.mask(frames, ids, where)]
        h, e = np.histogram(v[np.isfinite(v)], bins=bins, range=range, density=density)
        h = h.astype(float)
        if log:
            h[h == 0] = np.nan
        return 0.5 * (e[1:] + e[:-1]), h

    def pdf2d(self, qx, qy, frames=None, ids=None, where=None, bins=64, range=None,
              density=True):
        """Joint PDF of two quantities.  Returns (xc, yc, H[nx,ny])."""
        m = self.mask(frames, ids, where)
        a, b = self.get(qx)[m], self.get(qy)[m]
        ok = np.isfinite(a) & np.isfinite(b)
        H, xe, ye = np.histogram2d(a[ok], b[ok], bins=bins, range=range, density=density)
        return 0.5 * (xe[1:] + xe[:-1]), 0.5 * (ye[1:] + ye[:-1]), H

    def moments(self, quantity, where=None):
        """Per-snapshot (mean, std) of a quantity, each (Nt,)."""
        v = np.where(self.mask(where=where), self.get(quantity), np.nan)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            return np.nanmean(v, axis=1), np.nanstd(v, axis=1)

    def msd(self, axes=(0, 2), ids=None):
        """Mean-square displacement from each particle's first valid (unwrapped)
        position.  Returns (t - t0, msd[Nt])."""
        p = self.unwrapped()
        if ids is not None:
            p = p[:, np.isin(self.ids, ids)]
        first = p[np.argmax(np.isfinite(p[..., 0]), axis=0), np.arange(p.shape[1])]
        d2 = np.nansum((p - first)[..., list(axes)] ** 2, axis=-1)
        d2[np.isnan(p[..., 0])] = np.nan
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            return self.t - self.t[0], np.nanmean(d2, axis=1)

    def track(self, ids, unwrap=True):
        """Time series of selected particles: dict t, steps, ids, pos/vel/acc (Nt,n,3), age."""
        ids = np.atleast_1d(ids)
        j = np.searchsorted(self.ids, ids)
        if np.any(self.ids[np.clip(j, 0, len(self.ids) - 1)] != ids):
            raise KeyError(f"ids not in data: {ids[~np.isin(ids, self.ids)]}")
        pos = self.unwrapped() if unwrap else self.pos
        return dict(t=self.t, steps=self.steps, ids=ids, pos=pos[:, j], vel=self.vel[:, j],
                    acc=self.acc[:, j], age=self.age[:, j])

    # ════════════════════════════════════════════════════════════════════════
    # XDMF export (zero-copy byte-seek into the raw snapshot files)
    # ════════════════════════════════════════════════════════════════════════
    def write_xmf(self, out=None):
        """Write a ParaView XDMF time series for this case's particle snapshots.
        See write_particles_xmf for details; uses this instance's case_dir/prefix."""
        write_particles_xmf(self.case_dir, self.prefix, out)


def _detect_prefix(fields_dir):
    for p in Path(fields_dir).iterdir():
        m = _PAT.match(p.name)
        if m:
            return m.group(1)
    return None


def _probe_layout(fpath):
    """(total, ncol) for one snapshot file, from its size alone (no full read)."""
    size = Path(fpath).stat().st_size
    with open(fpath, "rb") as f:
        total = struct.unpack(">i", f.read(4))[0]
    if total == 0:
        return 0, 10
    ncol, remainder = divmod(size - 4 - 4 * total, 8 * total)
    if remainder != 0 or ncol not in (7, 10):
        raise ValueError(f"{fpath}: unexpected particle record size ({size} bytes, n={total})")
    return total, ncol


def write_particles_xmf(case_dir=".", prefix=None, out=None):
    """Write a ParaView XDMF time series (<out>, default 'paraview/particles.xmf')
    for '<prefix>_particles.<step>' snapshots, pointing directly at the original
    binary files via byte-seek HyperSlabs -- no particle data is copied.

    Point count varies from one snapshot to the next (particles exit/deposit/
    reinject), so every timestep declares its own Polyvertex Topology/Geometry
    sized to that snapshot's own particle count. Time Value is the step number
    (matching dopamine_post.fields.write_field_xmf's own convention), so
    ParaView's shared time toolbar scrubs fields and particles together.
    Handles both the legacy 7-column (x,y,z,u,v,w,age) and current 10-column
    (x,y,z,u,v,w,ax,ay,az,age) snapshot layouts, file by file.
    """
    case = Path(case_dir)
    fields_dir = case / "fields"
    out_dir = case / "paraview"
    out_dir.mkdir(exist_ok=True)
    out_path = Path(out) if out else out_dir / "particles.xmf"

    snaps = list_snapshots(fields_dir, prefix)
    lines = [
        '<?xml version="1.0" ?>',
        '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>',
        '<Xdmf Version="2.0">',
        '  <Domain>',
        '    <Grid Name="ParticleTimeSeries" GridType="Collection" CollectionType="Temporal">',
    ]
    n_written = 0
    for step, fpath in snaps:
        total, ncol = _probe_layout(fpath)
        if total == 0:
            continue
        frel = os.path.relpath(fpath, out_path.parent)
        age_col = 6 if ncol == 7 else 9
        dat_offset = 4 + 4 * total
        lines += [
            '',
            f'      <Grid Name="t{float(step):.6g}" GridType="Uniform">',
            f'        <Time Value="{float(step):.6g}"/>',
            f'        <Topology TopologyType="Polyvertex" NumberOfElements="{total}"/>',
            '        <Geometry GeometryType="XYZ">',
            f'          <DataItem ItemType="HyperSlab" Dimensions="{total} 3" Type="HyperSlab">',
            '            <DataItem Dimensions="3 2" Format="XML">',
            '              0 0',
            '              1 1',
            f'              {total} 3',
            '            </DataItem>',
            f'            <DataItem Dimensions="{total} {ncol}" Format="Binary"',
            f'                     DataType="Float" Precision="8" Endian="Big" Seek="{dat_offset}">',
            f'              {frel}',
            '            </DataItem>',
            '          </DataItem>',
            '        </Geometry>',
            '',
            '        <Attribute Name="Velocity" Center="Node" AttributeType="Vector">',
            f'          <DataItem ItemType="HyperSlab" Dimensions="{total} 3" Type="HyperSlab">',
            '            <DataItem Dimensions="3 2" Format="XML">',
            '              0 3',
            '              1 1',
            f'              {total} 3',
            '            </DataItem>',
            f'            <DataItem Dimensions="{total} {ncol}" Format="Binary"',
            f'                     DataType="Float" Precision="8" Endian="Big" Seek="{dat_offset}">',
            f'              {frel}',
            '            </DataItem>',
            '          </DataItem>',
            '        </Attribute>',
            '',
            '        <Attribute Name="age" Center="Node" AttributeType="Scalar">',
            f'          <DataItem ItemType="HyperSlab" Dimensions="{total} 1" Type="HyperSlab">',
            '            <DataItem Dimensions="3 2" Format="XML">',
            f'              0 {age_col}',
            '              1 1',
            f'              {total} 1',
            '            </DataItem>',
            f'            <DataItem Dimensions="{total} {ncol}" Format="Binary"',
            f'                     DataType="Float" Precision="8" Endian="Big" Seek="{dat_offset}">',
            f'              {frel}',
            '            </DataItem>',
            '          </DataItem>',
            '        </Attribute>',
            '',
            '        <Attribute Name="id" Center="Node" AttributeType="Scalar">',
            f'          <DataItem Dimensions="{total}" Format="Binary"',
            '                   DataType="Int" Precision="4" Endian="Big" Seek="4">',
            f'            {frel}',
            '          </DataItem>',
            '        </Attribute>',
            '      </Grid>',
        ]
        n_written += 1
    lines += ['    </Grid>', '  </Domain>', '</Xdmf>']

    if n_written == 0:
        raise SystemExit("No non-empty particle snapshots found -- nothing to write.")
    out_path.write_text('\n'.join(lines) + '\n')
    print(f"Wrote {out_path}  ({n_written} snapshots)")


# ════════════════════════════════════════════════════════════════════════════
# 3a. GEOMETRY  (solid wall / immersed body)
# ════════════════════════════════════════════════════════════════════════════
class Geometry:
    """The solid, for drawing.  Build with one of:

        Geometry.analytic(L, amp=0.05, wavelength=1.0, phase=0.0)   y_w = a + a cos(2 pi x/lambda + phase)
        Geometry.flat(L)                                            nothing drawn
        Geometry.from_sdf(case, sdf_file, grid_file=None, stride=2) sdf = 0 iso-surface of the solver SDF

    draw_side(ax)  -> solid in the x-y plane (analytic profile / mid-z SDF slice)
    draw_3d(ax)    -> opaque flat-shaded solid on 3-D axes (X=x, Y=z, Z=y)
    """

    def __init__(self, L, kind="flat", wall=None, sdf=None, stride=2):
        self.L, self.kind, self._wall, self._sdf, self.stride = tuple(L), kind, wall, sdf, stride
        self._mesh = None

    # -- constructors --
    @classmethod
    def flat(cls, L):
        return cls(L, "flat")

    @classmethod
    def analytic(cls, L, amp=0.05, wavelength=1.0, phase=0.0):
        return cls(L, "analytic",
                   wall=lambda x: amp + amp * np.cos(2 * np.pi * np.asarray(x) / wavelength + phase))

    @classmethod
    def from_sdf(cls, case=".", sdf_file="SDF_in", grid_file=None, stride=2):
        """Read the solver SDF (big-endian f8, Fortran order, shape (nxm+2, nym+2, nzm),
        phi<0 solid).  The y grid is stretched, so a grid.out is required (default
        <case>/fields/grid.out); the dimensions come from input_parameters."""
        case = Path(case)
        v = read_input(case, ("nx", "ny", "nz", "Lx", "Ly", "Lz"))
        if len(v) < 6:
            raise ValueError(f"nx,ny,nz,Lx,Ly,Lz not all found in {case}/input_parameters")
        n = (int(v["nx"]) - 1, int(v["ny"]) - 1, int(v["nz"]) - 1)
        L = (v["Lx"], v["Ly"], v["Lz"])
        gf = Path(grid_file) if grid_file else case / "fields" / "grid.out"
        if not gf.is_file():
            raise FileNotFoundError(f"{gf} not found: the y-grid is stretched, so pass grid_file= "
                                    "(--grid-file), a grid.out of a run with the same grid settings")
        ym = np.loadtxt(gf)[:, 1]
        if len(ym) != n[1]:
            raise ValueError(f"{gf} has {len(ym)} y cells, input_parameters implies {n[1]}")
        xm, zm = (np.arange(n[0]) + 0.5) * L[0] / n[0], (np.arange(n[2]) + 0.5) * L[2] / n[2]
        sdf = Path(sdf_file)
        sdf = sdf if sdf.is_absolute() or sdf.exists() else case / sdf
        raw = np.fromfile(sdf, dtype=">f8")
        ng = (n[0] + 2, n[1] + 2, n[2])
        if raw.size != np.prod(ng):
            raise ValueError(f"{sdf}: {raw.size} values, expected {ng} = {int(np.prod(ng))}")
        phi = raw.reshape(ng, order="F")[1:-1, 1:-1, :]
        return cls(L, "sdf", sdf=(np.clip(phi, -L[1], L[1]), xm, ym, zm), stride=stride)

    # -- drawing --
    def draw_side(self, ax, color=SAND):
        Lx = self.L[0]
        if self.kind == "analytic":
            xs = np.linspace(0, Lx, 600)
            ax.fill_between(xs, 0, self._wall(xs), color=color, zorder=0)
        elif self.kind == "sdf":
            phi, xm, ym, zm = self._sdf
            ax.contourf(xm, ym, phi[:, :, len(zm) // 2].T, levels=[-1e9, 0], colors=[color], zorder=0)

    def draw_3d(self, ax, zorder=1):
        if self.kind == "flat":
            return
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
        verts, faces = self.mesh()
        tri = verts[faces]
        nrm = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
        nrm /= np.linalg.norm(nrm, axis=1, keepdims=True) + 1e-30
        light = np.array([-0.4, 0.8, 0.45]); light /= np.linalg.norm(light)
        shade = 0.45 + 0.55 * np.abs(nrm @ light)                     # two-sided
        fc = np.clip(np.asarray(SAND)[None, :] * shade[:, None], 0, 1)
        ax.add_collection3d(Poly3DCollection(tri[:, :, [0, 2, 1]], facecolors=fc, edgecolors="none",
                                             linewidths=0, antialiaseds=False, zorder=zorder))

    def mesh(self):
        """Triangle mesh (verts[M,3] as x,y,z ; faces[K,3]); built once."""
        if self._mesh is None:
            self._mesh = self._wall_mesh() if self.kind == "analytic" else self._sdf_mesh()
        return self._mesh

    def _wall_mesh(self):
        Lx, _, Lz = self.L
        xs = np.linspace(0, Lx, 500)
        X, Z = np.meshgrid(xs, [0.0, Lz], indexing="ij")
        verts = np.stack([X.ravel(), self._wall(X).ravel(), Z.ravel()], axis=1)
        i = np.arange(len(xs) - 1) * 2
        faces = np.concatenate([np.stack([i, i + 1, i + 2], 1), np.stack([i + 1, i + 3, i + 2], 1)])
        return verts, faces

    def _sdf_mesh(self):
        try:
            from skimage.measure import marching_cubes
        except ImportError:
            raise SystemExit("--sdf-file needs scikit-image: pip install scikit-image")
        phi, xs, ym, zs = self._sdf
        s, L = self.stride, self.L
        if s > 1:
            phi, xs, zs = phi[::s, :, ::s], xs[::s], zs[::s]
        # One layer of linear extrapolation to the domain faces, so a surface that reaches
        # a face (e.g. dune troughs at y=0) is still found instead of leaving a gap.
        for ax_, c, hi_ext in ((0, xs, True), (1, ym, False), (2, zs, True)):
            a0, a1 = phi.take([0], ax_), phi.take([1], ax_)
            lo = a0 + (a0 - a1) * (c[0] / (c[1] - c[0]))
            b0, b1 = phi.take([-1], ax_), phi.take([-2], ax_)
            hi = b0 + (b0 - b1) * ((L[ax_] - c[-1]) / (c[-1] - c[-2])) if hi_ext else b0
            phi = np.concatenate([lo, phi, hi], axis=ax_)
        coords = [np.concatenate([[0.0], c, [L[a]]]) for a, c in enumerate((xs, ym, zs))]
        v, f, _, _ = marching_cubes(phi, level=0.0)                   # index-space vertices
        verts = np.stack([np.interp(v[:, a], np.arange(len(c)), c) for a, c in enumerate(coords)], axis=1)
        return verts, f


# ════════════════════════════════════════════════════════════════════════════
# 3b. FIGURES / MOVIES
# ════════════════════════════════════════════════════════════════════════════
def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def _broken(xy, L):
    """Insert NaN rows where a periodic coordinate jumps by more than L/2 (no wrap streaks)."""
    jump = np.abs(np.diff(xy[:, 0])) > 0.5 * L[0]
    if L[1] is not None:
        jump |= np.abs(np.diff(xy[:, 1])) > 0.5 * L[1]
    return np.insert(xy, np.where(jump)[0] + 1, np.nan, axis=0)


def _axes3d(fig, L, ytop, yexag, zoom, elev, azim):
    ax = fig.add_subplot(projection="3d", facecolor="white")
    ax.computed_zorder = False                      # honour zorder: solid first, particles on top
    ax.set_proj_type("persp", focal_length=0.35)
    ax.view_init(elev=elev, azim=azim)
    ax.set_box_aspect((L[0], L[2], ytop * yexag), zoom=zoom)
    ax.set_xlim(0, L[0]); ax.set_ylim(0, L[2]); ax.set_zlim(0, ytop)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.set_pane_color((1, 1, 1, 0))
        axis._axinfo["grid"]["color"] = (0.85, 0.85, 0.85, 1)
        axis.line.set_color("#888888")
    for lab, f in (("$x$", ax.set_xlabel), ("$z$", ax.set_ylabel), ("$y$", ax.set_zlabel)):
        f(lab, labelpad=-8, fontsize=13)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
    for zz in (0.0, ytop):
        ax.plot([0, L[0], L[0], 0, 0], [0, 0, L[2], L[2], 0], [zz] * 5, color="#999999", lw=0.7, zorder=0)
    for cx, cy in ((0, 0), (L[0], 0), (L[0], L[2]), (0, L[2])):
        ax.plot([cx, cx], [cy, cy], [0, ytop], color="#999999", lw=0.7, zorder=0)
    return ax


def _window(pdat, frames):
    return np.arange(len(pdat.t))[slice(-40, None) if frames is None else frames]


def _view_opts(o):
    d = dict(tail=0, bg=0, ymax=None, yexag=1.0, zoom=1.0, elev=22.0, azim=-62.0, rotate=0.0,
             size=30.0, seed=1, fps=8, width=900, workers=8)
    bad = set(o) - set(d)
    if bad:
        raise TypeError(f"unknown options {sorted(bad)}")
    d.update(o)
    return d


# -- tracks ------------------------------------------------------------------
def _tracks_frame(pdat, ids, idx, geo, o, view, i, start_marks=False):
    """Figure showing the tracks of `ids` up to window frame i."""
    plt = _plt()
    L = pdat.L
    ytop = o["ymax"] or L[1]
    pos = pdat.track(ids, unwrap=False)["pos"][idx]                 # (nw, n, 3)
    pos = np.where((pos[..., 1] > ytop)[..., None], np.nan, pos)    # 3-D axes do not clip
    cmap = plt.get_cmap("tab10")
    lo = max(0, i - o["tail"]) if o["tail"] > 0 else 0
    title = (f"$t = {pdat.t[idx[i]]:.2f}$   (step {pdat.steps[idx[i]]})   "
             f"{len(ids)} tracked particles")

    if view == "2d":
        fig, (a1, a2) = plt.subplots(2, 1, figsize=(10, 7), sharex=True,
                                     gridspec_kw=dict(height_ratios=[1, 1.4]), facecolor="white")
        geo.draw_side(a1)
        for k, pid in enumerate(ids):
            c, p = cmap(k % 10), pos[lo:i + 1, k]
            for ax, col, Lc in ((a1, 1, None), (a2, 2, L[2])):
                seg = _broken(p[:, [0, col]], (L[0], Lc))
                ax.plot(seg[:, 0], seg[:, 1], color=c, lw=1.2, label=f"id {pid}" if ax is a1 else None)
                v = np.where(np.isfinite(p[:, 0]))[0]
                if len(v):
                    if start_marks:
                        ax.plot(p[v[0], 0], p[v[0], col], "o", mfc="none", color=c, ms=6)
                    ax.plot(p[v[-1], 0], p[v[-1], col], "o", color=c, ms=5)
        a1.set_xlim(0, L[0]); a1.set_ylim(0, ytop); a1.set_ylabel("$y$")
        a2.set_ylim(0, L[2]); a2.set_ylabel("$z$"); a2.set_xlabel("$x$")
        a1.legend(ncol=4, fontsize=8, loc="upper right")
        a1.set_title(title + ("   (○ start, ● end)" if start_marks else ""), fontsize=10)
        fig.tight_layout()
        return fig

    fig = plt.figure(figsize=(9, 5.6), dpi=100, facecolor="white")
    ax = _axes3d(fig, L, ytop, o["yexag"], o["zoom"], o["elev"],
                 o["azim"] + o["rotate"] * i / max(1, len(idx) - 1))
    geo.draw_3d(ax, zorder=1)
    if o["bg"] > 0:
        rest = np.setdiff1d(pdat.ids, ids)
        sel = np.sort(np.random.default_rng(o["seed"] + 1).choice(rest, min(o["bg"], len(rest)), replace=False))
        b = pdat.pos[idx[i]][np.searchsorted(pdat.ids, sel)]
        ok = np.isfinite(b[:, 0]) & (b[:, 1] <= ytop)
        ax.scatter(b[ok, 0], b[ok, 2], b[ok, 1], s=2, c="#999999", alpha=0.35, linewidths=0,
                   depthshade=True, zorder=3)
    segs, cols = [], []
    for k in range(len(ids)):
        p = pos[lo:i + 1, k]
        for j in range(len(p) - 1):
            s = p[j:j + 2]
            if np.isnan(s).any() or abs(s[1, 0] - s[0, 0]) > 0.5 * L[0] or abs(s[1, 2] - s[0, 2]) > 0.5 * L[2]:
                continue                                              # gap or periodic wrap
            segs.append(s[:, [0, 2, 1]])
            cols.append(cmap(k % 10)[:3] + (0.25 + 0.75 * (j + 1) / len(p),))
    if segs:
        from mpl_toolkits.mplot3d.art3d import Line3DCollection
        ax.add_collection3d(Line3DCollection(segs, colors=cols, linewidths=1.8, zorder=4))
    cur = pos[i]
    ok = np.isfinite(cur[:, 0])
    ax.scatter(cur[ok, 0], cur[ok, 2], cur[ok, 1], s=o["size"], c=[cmap(k % 10) for k in np.where(ok)[0]],
               edgecolors="k", linewidths=0.5, depthshade=False, zorder=6)
    fig.suptitle(title, fontsize=10, color="#444444")
    fig.tight_layout(pad=0.4)
    return fig


def plot_tracks(pdat, ids, save=None, frames=slice(-20, None), geometry=None, view="2d", **opts):
    """Static picture of the tracks over `frames` (default last 20 snapshots).
    view '2d' (side + top views, ○ start ● end) or '3d'.  Returns the Figure."""
    o = _view_opts(opts)
    idx = _window(pdat, frames)
    geo = geometry or Geometry.flat(pdat.L)
    fig = _tracks_frame(pdat, np.atleast_1d(ids), idx, geo, o, view, len(idx) - 1, start_marks=True)
    if save:
        fig.savefig(save, dpi=150, facecolor="white")
        print(f"Wrote {save}")
    return fig


def animate_tracks(pdat, ids, out, frames=slice(-40, None), geometry=None, view="3d", **opts):
    """Movie of the tracks over `frames`: growing trails + current-position markers.

    view 3d|2d;  options: tail (trail length in snapshots, 0 = all), bg (faint background
    particles, 3d), ymax, yexag, zoom, elev, azim, rotate, size, seed, fps, width, workers.
    Points above ymax are hidden."""
    o = _view_opts(opts)
    idx = _window(pdat, frames)
    geo = geometry or Geometry.flat(pdat.L)
    ids = np.atleast_1d(ids)
    if view == "3d":
        geo.mesh() if geo.kind != "flat" else None                    # build before forking
    render_movie(lambda i: _tracks_frame(pdat, ids, idx, geo, o, view, i), len(idx), out,
                 o["fps"], o["width"], o["workers"])


# -- particle cloud ----------------------------------------------------------
def animate_cloud(pdat, out, frames=slice(-40, None), geometry=None, color="speed", max_particles=0,
                  s=3.0, **opts):
    """3-D movie of all particles coloured by 'speed' | 'height' | 'id' (global colour range)."""
    plt = _plt()
    o = _view_opts(opts)
    idx = _window(pdat, frames)
    geo = geometry or Geometry.flat(pdat.L)
    L = pdat.L
    ytop = o["ymax"] or L[1]
    val = {"speed": lambda: pdat.get("speed"), "height": lambda: pdat.get("y"),
           "id": lambda: np.broadcast_to(pdat.ids.astype(float), pdat.get("y").shape)}[color]()
    lab = {"speed": r"particle speed $|\mathbf{u}_p|$", "height": "particle height $y$", "id": "particle id"}[color]
    v = val[idx]
    vmin, vmax = _core.percentile_clim(v, 1, 99)
    cmap = plt.get_cmap({"speed": "plasma", "height": "viridis", "id": "turbo"}[color])
    keep = (np.sort(np.random.default_rng(0).choice(len(pdat.ids), max_particles, replace=False))
            if 0 < max_particles < len(pdat.ids) else slice(None))
    if geo.kind != "flat":
        geo.mesh()

    def frame(i):
        fig = plt.figure(figsize=(9, 5.6), dpi=100, facecolor="white")
        ax = _axes3d(fig, L, ytop, o["yexag"], o["zoom"], o["elev"],
                     o["azim"] + o["rotate"] * i / max(1, len(idx) - 1))
        geo.draw_3d(ax, zorder=1)
        p, c = pdat.pos[idx[i]][keep], val[idx[i]][keep]
        ok = np.isfinite(p[:, 0]) & (p[:, 1] <= ytop)
        p, c = p[ok], c[ok]
        order = np.argsort(c)                                         # high values drawn last
        sc = ax.scatter(p[order, 0], p[order, 2], p[order, 1], c=c[order], cmap=cmap, vmin=vmin, vmax=vmax,
                        s=s, alpha=0.9, linewidths=0, depthshade=True, zorder=5)
        cb = fig.colorbar(sc, ax=ax, shrink=0.5, pad=0.02)
        cb.set_label(lab, fontsize=9)
        cb.ax.tick_params(labelsize=8, colors="#444444")
        fig.suptitle(f"$t = {pdat.t[idx[i]]:.2f}$   (step {pdat.steps[idx[i]]})   N = {ok.sum():,}",
                     fontsize=10, color="#444444")
        fig.tight_layout(pad=0.4)
        return fig

    render_movie(frame, len(idx), out, o["fps"], o["width"], o["workers"])


# -- PDFs --------------------------------------------------------------------
def plot_pdfs(pdat, quantities, save=None, frames=slice(-20, None), log=False, **pdf_kw):
    """Overlay the PDFs of several quantities (pooled over `frames`) on one axes."""
    plt = _plt()
    fig, ax = plt.subplots(figsize=(5.5, 3.8))
    for q in quantities:
        x, p = pdat.pdf(q, frames=frames, log=log, **pdf_kw)
        ax.plot(x, p, label=q)
    ax.set_ylabel("PDF"); ax.set_xlabel(" / ".join(quantities))
    if log:
        ax.set_yscale("log")
    ax.legend()
    fig.tight_layout()
    if save:
        fig.savefig(save, dpi=150)
        print(f"Wrote {save}")
    return fig
