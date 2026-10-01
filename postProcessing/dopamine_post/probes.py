"""
probes.py -- fdm-dopamine probe_output.f90 line and slice probes.

Two probe kinds, one module:

  LineProbe   a 1-D line probe: `<base>.bin` (big-endian float64, Fortran
              stream, no record markers) + companion `<base>_meta.txt`.
              Each snapshot is Fortran array out1d(ncomp, npts).

  SliceProbe  a 2-D plane probe: `<base>.bin` + `<base>_meta.txt`, each
              snapshot Fortran array out2d(ncomp, n1, n2).  `.write_xmf(...)`
              ports generate_slice_xmf.py's byte-seek XDMF writer (no data
              copy -- ParaView reads straight out of the original .bin).

Library use
-----------
    from dopamine_post.probes import LineProbe, SliceProbe

    line = LineProbe.read("case/line_meta.txt")
    u = line.get("U")                       # (nsnaps, npts)

    slc = SliceProbe.read("case/slice_meta.txt")
    slc.write_xmf(snap_path="case/fields/grid.out")

Both share `_core.read_meta`/`_core.comp_names` for the identical meta.txt
key=value parsing that generate_slice_xmf.py and load_line_probes.py used to
each carry a private (verbatim-duplicate) copy of.
"""
import os
import warnings
from pathlib import Path

import numpy as np

from . import _core
from ._core import comp_names, read_meta


# ════════════════════════════════════════════════════════════════════════════
# line probes
# ════════════════════════════════════════════════════════════════════════════
class LineProbe:
    """A 1-D line probe: `data` has shape (nsnaps, ncomp, npts).

    u = probe.get("U")                          # (nsnaps, npts)
    u = probe.data[:, probe.comps.index("U"), :] # equivalent, by hand
    """

    def __init__(self, data, comps, direction, meta, meta_path, bin_path):
        self.data = data
        self.comps = comps
        self.dir = direction
        self.meta = meta
        self.meta_path, self.bin_path = Path(meta_path), Path(bin_path)
        self.nsnaps, self.ncomp, self.npts = data.shape

    def __repr__(self):
        return (f"LineProbe(dir={self.dir!r}, comps={self.comps}, "
                f"nsnaps={self.nsnaps}, npts={self.npts})")

    @classmethod
    def read(cls, path):
        """Load a line probe. `path` is either the `_meta.txt` or the `.bin`
        file; the companion file is located automatically alongside it."""
        meta_path, bin_path = _companion_paths(path, "_meta.txt")
        meta = read_meta(meta_path)
        nc, npts, nsnaps = meta["ncomp"], meta["npts"], meta["nsnaps"]
        direction = meta.get("dir", "?").strip()
        comps = comp_names(meta.get("comps", "UVW"))

        if not bin_path.is_file():
            raise FileNotFoundError(f"binary file not found: {bin_path}")
        snap_bytes = nc * npts * 8
        actual = bin_path.stat().st_size // snap_bytes if snap_bytes else 0
        if actual != nsnaps:
            warnings.warn(f"{meta_path}: meta says {nsnaps} snaps, file has {actual}; "
                          "using file count")
            nsnaps = actual

        # Each snapshot: Fortran (nc, npts) col-major == C (npts, nc) row-major.
        raw = np.fromfile(bin_path, dtype=">f8", count=nsnaps * nc * npts)
        data = raw.reshape(nsnaps, npts, nc).transpose(0, 2, 1).astype("f8", copy=False)
        return cls(data, comps, direction, meta, meta_path, bin_path)

    def get(self, name):
        """One component's (nsnaps, npts) array, by name (e.g. 'U')."""
        return self.data[:, self.comps.index(name), :]


# ════════════════════════════════════════════════════════════════════════════
# slice probes
# ════════════════════════════════════════════════════════════════════════════
class SliceProbe:
    """A 2-D plane probe: `data` has shape (nsnaps, ncomp, n1, n2).

    dir='x' (x-normal): axis-1 = y (n1 = nym_global), axis-2 = z (n2 = nzm_global)
    dir='y' (y-normal): axis-1 = x (n1 = nxm_global), axis-2 = z (n2 = nzm_global)
    dir='z' (z-normal): axis-1 = x (n1 = nxm_global), axis-2 = y (n2 = nym_global)
    """

    def __init__(self, data, comps, direction, meta, meta_path, bin_path, times):
        self.data = data
        self.comps = comps
        self.dir = direction
        self.meta = meta
        self.meta_path, self.bin_path = Path(meta_path), Path(bin_path)
        self.pos = meta.get("pos", 0.0)
        self.times = times
        self.nsnaps, self.ncomp, self.n1, self.n2 = data.shape

    def __repr__(self):
        return (f"SliceProbe(dir={self.dir!r}, comps={self.comps}, "
                f"nsnaps={self.nsnaps}, shape=({self.n1},{self.n2}))")

    @classmethod
    def read(cls, path):
        """Load a slice probe. `path` is either the `_meta.txt` or the `.bin`
        file; the companion file is located automatically alongside it."""
        meta_path, bin_path = _companion_paths(path, "_meta.txt")
        meta = read_meta(meta_path)
        nc, n1, n2, nsnaps = meta["ncomp"], meta["n1"], meta["n2"], meta["nsnaps"]
        direction = meta.get("dir", "z").strip().lower()[0]
        comps = comp_names(meta.get("comps", "UVW"))

        if not bin_path.is_file():
            raise FileNotFoundError(f"binary file not found: {bin_path}")
        snap_bytes = nc * n1 * n2 * 8
        actual = bin_path.stat().st_size // snap_bytes if snap_bytes else 0
        if actual != nsnaps:
            warnings.warn(f"{meta_path}: meta says {nsnaps} snaps, file has {actual}; "
                          "using file count")
            nsnaps = actual

        times = _times_for(meta, meta_path, nsnaps)

        # Each snapshot: Fortran (nc, n1, n2) col-major == C (n2, n1, nc) row-major.
        raw = np.fromfile(bin_path, dtype=">f8", count=nsnaps * nc * n1 * n2)
        data = (raw.reshape(nsnaps, n2, n1, nc).transpose(0, 3, 2, 1)
                   .astype("f8", copy=False))
        return cls(data, comps, direction, meta, meta_path, bin_path, times)

    def get(self, name):
        """One component's (nsnaps, n1, n2) array, by name (e.g. 'U')."""
        return self.data[:, self.comps.index(name), :, :]

    # ── XDMF export (zero-copy byte-seek into the raw .bin, like ParticleData.write_xmf) ──
    def write_xmf(self, out=None, snap_path=None, dt=1.0, t0=0.0):
        """Write `<base>.xmf` plus `<base>_ax1.bin`/`_ax2.bin`/`_axn.bin` coordinate
        files for this slice, readable by ParaView >= 5 / VisIt >= 3.

        snap_path  fields/grid.out (fast) or a full solver snapshot, used for
                   exact physical grid coordinates; falls back to input_parameters
                   (uniform-grid approximation) and finally integer indices.
        dt, t0     fallback uniform time axis (t = t0 + snapshot_index * dt),
                   used only when this meta file's `times` key (pointing at a
                   companion `<base>_times.bin` of exact simulation times,
                   written by probe_output.f90) is absent or mismatched.
        """
        meta, meta_path, bin_path = self.meta, self.meta_path, self.bin_path
        nc, n1, n2 = meta["ncomp"], self.n1, self.n2
        nsnaps = self.nsnaps
        dirstr = self.dir
        comps = self.comps

        stem = meta_path.stem
        base = meta_path.parent / (stem[:-5] if stem.endswith("_meta") else stem)
        xmf_path = Path(out) if out else base.with_suffix(".xmf")
        ax1_path = base.parent / (base.name + "_ax1.bin")
        ax2_path = base.parent / (base.name + "_ax2.bin")
        ax_n_path = base.parent / (base.name + "_axn.bin")

        if nsnaps == 0:
            raise SystemExit(f"{meta_path.name}: no snapshots written yet -- nothing to write")

        times = self.times
        if times is None:
            times = t0 + np.arange(nsnaps) * dt

        # ── physical coordinates ──
        xm = ym = zm = None
        if snap_path is not None:
            xm, ym, zm = _coords_from_snap(snap_path)
        if xm is None:
            result = _uniform_coords(meta_path)
            if result is not None:
                xm, ym, zm = result
            else:
                n_max = max(n1, n2) + 1
                xm = ym = zm = np.arange(n_max, dtype=np.float64)

        if dirstr == "x":
            ax1, ax2 = ym[:n1], zm[:n2]
        elif dirstr == "y":
            ax1, ax2 = xm[:n1], zm[:n2]
        else:
            ax1, ax2 = xm[:n1], ym[:n2]

        pos = self.pos
        if dirstr == "x":
            ax_normal = np.array([xm[np.argmin(np.abs(xm - pos))]])
        elif dirstr == "y":
            ax_normal = np.array([ym[np.argmin(np.abs(ym - pos))]])
        else:
            ax_normal = np.array([zm[np.argmin(np.abs(zm - pos))]])

        np.asarray(ax1, dtype="<f8").tofile(ax1_path)
        np.asarray(ax2, dtype="<f8").tofile(ax2_path)
        np.asarray(ax_normal, dtype="<f8").tofile(ax_n_path)

        xmf_dir = xmf_path.parent
        rel_bin = os.path.relpath(bin_path, xmf_dir)
        rel_ax1 = os.path.relpath(ax1_path, xmf_dir)
        rel_ax2 = os.path.relpath(ax2_path, xmf_dir)
        rel_axn = os.path.relpath(ax_n_path, xmf_dir)

        snap_bytes = nc * n1 * n2 * 8

        # Embed slice as a 3DRectMesh Dimensions="d0 d1 d2" (slowest->fastest),
        # one dimension = 1 for the slice-normal direction.
        if dirstr == "x":
            d0, d1, d2 = n2, n1, 1
            geo = (
                '      <Geometry GeometryType="VxVyVz">\n'
                f'        <DataItem Dimensions="1" Format="Binary" DataType="Float"'
                f' Precision="8" Endian="Little">\n          {rel_axn}\n        </DataItem>\n'
                f'        <DataItem Dimensions="{n1}" Format="Binary" DataType="Float"'
                f' Precision="8" Endian="Little">\n          {rel_ax1}\n        </DataItem>\n'
                f'        <DataItem Dimensions="{n2}" Format="Binary" DataType="Float"'
                f' Precision="8" Endian="Little">\n          {rel_ax2}\n        </DataItem>\n'
                '      </Geometry>\n'
            )
        elif dirstr == "y":
            d0, d1, d2 = n2, 1, n1
            geo = (
                '      <Geometry GeometryType="VxVyVz">\n'
                f'        <DataItem Dimensions="{n1}" Format="Binary" DataType="Float"'
                f' Precision="8" Endian="Little">\n          {rel_ax1}\n        </DataItem>\n'
                f'        <DataItem Dimensions="1" Format="Binary" DataType="Float"'
                f' Precision="8" Endian="Little">\n          {rel_axn}\n        </DataItem>\n'
                f'        <DataItem Dimensions="{n2}" Format="Binary" DataType="Float"'
                f' Precision="8" Endian="Little">\n          {rel_ax2}\n        </DataItem>\n'
                '      </Geometry>\n'
            )
        else:
            d0, d1, d2 = 1, n2, n1
            geo = (
                '      <Geometry GeometryType="VxVyVz">\n'
                f'        <DataItem Dimensions="{n1}" Format="Binary" DataType="Float"'
                f' Precision="8" Endian="Little">\n          {rel_ax1}\n        </DataItem>\n'
                f'        <DataItem Dimensions="{n2}" Format="Binary" DataType="Float"'
                f' Precision="8" Endian="Little">\n          {rel_ax2}\n        </DataItem>\n'
                f'        <DataItem Dimensions="1" Format="Binary" DataType="Float"'
                f' Precision="8" Endian="Little">\n          {rel_axn}\n        </DataItem>\n'
                '      </Geometry>\n'
            )

        parts = [f'<?xml version="1.0" ?>\n<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n'
                 f'<Xdmf Version="2.0">\n<Domain>\n'
                 f'  <Grid Name="{base.name}" GridType="Collection" CollectionType="Temporal">\n']
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
        parts.append('  </Grid>\n</Domain>\n</Xdmf>\n')

        xmf_path.write_text(''.join(parts))
        print(f"Wrote {xmf_path}  ({nsnaps} snaps, comps={comps})")
        return xmf_path


def _hyperslab_attr(cname, ci, nc, d0, d1, d2, seek, rel_bin):
    """XDMF lines for one scalar Attribute extracted via HyperSlab.

    Source layout (Fortran stream write of array(nc, n1, n2) col-major) is
    treated as C array (d0, d1, d2, nc) row-major, where exactly one of
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


def _times_for(meta, meta_path, nsnaps):
    """Exact simulation times from the `<base>_times.bin` a meta's `times` key
    points at (big-endian float64, one per snapshot), or None if absent/stale."""
    times_key = meta.get("times")
    if not times_key:
        return None
    candidates = ([Path(times_key)] if Path(times_key).is_absolute() else
                  [meta_path.parent / Path(times_key).name, Path(times_key)])
    times_path = next((c for c in candidates if c.is_file()), candidates[0])
    if not times_path.is_file():
        return None
    times = np.fromfile(times_path, dtype=">f8")
    if len(times) != nsnaps:
        warnings.warn(f"{times_path.name} has {len(times)} times, expected {nsnaps}; "
                      "falling back to --dt/--t0 at write_xmf time")
        return None
    return times


def _coords_from_snap(snap_path):
    """(xm, ym, zm) cell-centre coordinates, from fields/grid.out (fast path,
    needs a sibling geometry.out) or the grid header of a full solver snapshot."""
    snap_path = Path(snap_path)
    if snap_path.name == "grid.out":
        geo_path = snap_path.parent / "geometry.out"
        if not geo_path.is_file():
            raise FileNotFoundError(f"geometry.out not found alongside {snap_path}")
        tokens = geo_path.read_text().split()
        nxm, nym, nzm = int(tokens[0]), int(tokens[1]), int(tokens[2])
        Lx, Lz = float(tokens[3]), float(tokens[5])
        ym = np.loadtxt(snap_path)[:, 1]
        xm = (np.arange(nxm) + 0.5) * (Lx / nxm)
        zm = (np.arange(nzm) + 0.5) * (Lz / nzm)
        return xm, ym, zm
    # Full snapshot: only the grid header is needed, so skip past it without
    # reading the (much larger) field blocks that follow.
    r = _core.BinaryReader(snap_path.read_bytes())
    for _ in range(3):
        r.rd(r.ri())          # x, y, z face-point coords (discarded)
    xm = r.rd(r.ri())
    ym = r.rd(r.ri())
    zm = r.rd(r.ri())
    return xm, ym, zm


def _uniform_coords(meta_path):
    """(xm, ym, zm) approximated from a nearby input_parameters, or None."""
    for candidate in (meta_path.parent / "input_parameters",
                       meta_path.parent.parent / "input_parameters",
                       Path("input_parameters")):
        if candidate.is_file():
            ip_path = candidate
            break
    else:
        return None
    try:
        ip = _core.parse_input_parameters(ip_path)
    except Exception:
        return None
    nx, ny, nz = int(ip.get("nx", 0)), int(ip.get("ny", 0)), int(ip.get("nz", 0))
    Lx, Ly, Lz = float(ip.get("lx", 1)), float(ip.get("ly", 1)), float(ip.get("lz", 1))
    if min(nx, ny, nz) <= 0:
        return None
    nxm, nym, nzm = nx - 1, ny - 1, nz - 1
    xm = (np.arange(nxm) + 0.5) * (Lx / nxm)
    ym = (np.arange(nym) + 0.5) * (Ly / nym)
    zm = (np.arange(nzm) + 0.5) * (Lz / nzm)
    return xm, ym, zm


def _companion_paths(path, meta_suffix):
    """(meta_path, bin_path) whichever of the two `path` names."""
    p = Path(path)
    if p.suffix == ".bin":
        return p.parent / (p.stem + meta_suffix), p
    stem = p.stem
    base = stem[:-5] if stem.endswith("_meta") else stem
    return p, p.parent / (base + ".bin")
