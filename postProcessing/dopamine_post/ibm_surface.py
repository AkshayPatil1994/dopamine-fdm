"""
ibm_surface.py -- fdm-dopamine immersed-boundary surface samples.

Each `ibm_surface/surface.<step>.bin` file is written by sample_ibm_surface
(src/ibm.f90): a point cloud over the immersed body surface (phi=0), one
point per cell-centre interface cell.

Binary layout (big-endian, Fortran stream, no record markers)
---------------------------------------------------------------
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

Summing fp_* / fv_* over all points reproduces the Method-2 pressure/viscous
drag reported in ibm_forces.csv.

Library use
-----------
    from dopamine_post.ibm_surface import IBMSurface

    surf = IBMSurface.read("ibm_surface/surface.00010000.bin")
    print(surf.t, surf.xyz.shape, surf.p.mean())

    surf.to_vtp("surface.00010000.vtp")                       # one snapshot
    IBMSurface.to_pvd("ibm_surface/surface.*.bin", "ibm_surface/surface.pvd")  # whole series

Requires `pyvista` (pip install pyvista) for `.to_vtp`/`.to_pvd`.
"""
import glob
import struct
from pathlib import Path

import numpy as np

from ._core import write_pvd

# big-endian: Int32 count, float64 time, then 13 float64 blocks of n values each
_F8 = np.dtype(">f8")
_FIELDS = ("x", "y", "z", "nx", "ny", "nz",
           "p", "fp_x", "fp_y", "fp_z", "fv_x", "fv_y", "fv_z")


class IBMSurface:
    """One IBM surface sample: `t` (time), `n` (point count), and per-point
    arrays `xyz`, `normal`, `p`, `f_pres`, `f_visc` (each (n,3) except `p`)."""

    def __init__(self, t, xyz, normal, p, f_pres, f_visc, cols, path=None):
        self.t = t
        self.n = xyz.shape[0]
        self.xyz, self.normal, self.p = xyz, normal, p
        self.f_pres, self.f_visc = f_pres, f_visc
        self.cols = cols            # raw named columns (x,y,z,nx,...,fv_z), 1-D each
        self.path = Path(path) if path else None

    def __repr__(self):
        return f"IBMSurface(t={self.t:.6g}, n={self.n}, path={self.path})"

    @classmethod
    def read(cls, path):
        """Read one `surface.<step>.bin` file."""
        raw = Path(path).read_bytes()
        off = 0
        n = struct.unpack(">i", raw[off:off + 4])[0]
        off += 4
        t = struct.unpack(">d", raw[off:off + 8])[0]
        off += 8

        cols = {}
        nbytes = n * 8
        for name in _FIELDS:
            cols[name] = np.frombuffer(raw[off:off + nbytes], dtype=_F8).astype(np.float64)
            off += nbytes

        xyz = np.column_stack((cols["x"], cols["y"], cols["z"]))
        normal = np.column_stack((cols["nx"], cols["ny"], cols["nz"]))
        f_pres = np.column_stack((cols["fp_x"], cols["fp_y"], cols["fp_z"]))
        f_visc = np.column_stack((cols["fv_x"], cols["fv_y"], cols["fv_z"]))
        return cls(t, xyz, normal, cols["p"], f_pres, f_visc, cols, path=path)

    def forces(self):
        """(F_pressure, F_viscous), each a (3,) array: the per-point forces
        summed over the whole surface (reproduces ibm_forces.csv's Method-2)."""
        return self.f_pres.sum(axis=0), self.f_visc.sum(axis=0)

    def to_vtp(self, out=None):
        """Write this sample as a VTK PolyData (.vtp) point cloud.

        Point arrays: pressure, normal, pressure_force, viscous_force,
        total_force, pressure_force_mag, viscous_force_mag. Needs pyvista.
        """
        import pyvista as pv

        cloud = pv.PolyData(self.xyz)
        cloud["pressure"] = self.p
        cloud["normal"] = self.normal
        cloud["pressure_force"] = self.f_pres
        cloud["viscous_force"] = self.f_visc
        cloud["total_force"] = self.f_pres + self.f_visc
        cloud["pressure_force_mag"] = np.linalg.norm(self.f_pres, axis=1)
        cloud["viscous_force_mag"] = np.linalg.norm(self.f_visc, axis=1)
        cloud.field_data["time"] = np.array([self.t])

        if out is None:
            if self.path is None:
                raise ValueError("pass out= (no source path to derive it from)")
            out = str(Path(self.path).with_suffix(".vtp"))
        cloud.save(out)
        return out

    @staticmethod
    def to_pvd(pattern="ibm_surface/surface.*.bin", out="ibm_surface/surface.pvd"):
        """Convert every `surface.<step>.bin` matching `pattern` to .vtp and
        write a ParaView `.pvd` time-series collection (via `_core.write_pvd`)
        so the whole run loads as one animation in ParaView. Needs pyvista.
        """
        files = sorted(glob.glob(pattern))
        if not files:
            raise FileNotFoundError(f"no files match: {pattern}")

        pvd_path = Path(out)
        pvd_path.parent.mkdir(parents=True, exist_ok=True)

        entries = []
        for f in files:
            surf = IBMSurface.read(f)
            vtp_out = surf.to_vtp()
            entries.append((surf.t, vtp_out))
            print(f"  {Path(f).name}  t={surf.t:.6g}  -> {vtp_out}")

        write_pvd(pvd_path, entries)
        return pvd_path
