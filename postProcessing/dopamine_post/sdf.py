"""
sdf.py -- read and plot the cell-centre signed-distance field (SDF_in) fdm-dopamine
reads when ibm_input_mode=1.

    from dopamine_post.sdf import SDF

    s = SDF.read(".")                               # reads SDF_in + input_parameters
    fig = s.plot_xz(y=0.05)                          # X-Z plane at y ~= 0.05
    fig = s.plot_xy(z=3.0)                           # X-Y plane at z ~= 3.0

Binary layout of SDF_in (written by GenSDF when ibm_input_mode=1):
    big-endian float64, Fortran column-major order
    shape  (nxg, nyg, nzm)  =  (nxm+2, nym+2, nzm)
    ghost layers present in x and y; no ghost layer in z
    -> strip to (nxm, nym, nzm) by removing the first/last index in x and y

Grid coordinates are read from fields/grid.out + fields/geometry.out (written by
genGridandIC.f90) so the same non-uniform y-grid used by the solver is reproduced
exactly, without needing a solver snapshot.

This module owns its own SDF reader; it is deliberately independent from
`dopamine_post.particles.Geometry.from_sdf` (which reads SDF_in for drawing a solid
under particle tracks) to avoid a cross-module dependency during parallel development.
"""
from pathlib import Path

import numpy as np

from ._core import parse_input_parameters

# GenSDF's flood_fill_mod.f90 SENTINEL_FRAC: cells outside its geometry-focused AABB
# are left at the raw background sentinel value (e.g. 1e10, "far fluid, never
# computed"); a cell is treated as sentinel once |phi| exceeds this fraction of the
# array's own max |phi|.
SENTINEL_FRAC = 0.5


def _nearest_index(coords, value):
    return int(np.argmin(np.abs(coords - value)))


def coords_from_grid(grid_path):
    """(xm, ym, zm) cell-centre coordinates from fields/grid.out + fields/geometry.out.

    grid.out format (written by genGridandIC.f90)::

        i  ym_i  y_face_{i+1}  dy_i  1/dy_i

    geometry.out format::

        nxm  nym  nzm
        Lx   Ly   Lz

    x and z grids are always uniform in fdm-dopamine.
    """
    grid_path = Path(grid_path)
    geo_path = grid_path.parent / "geometry.out"
    if not geo_path.is_file():
        raise FileNotFoundError(f"geometry.out not found alongside {grid_path}")

    tokens = geo_path.read_text().split()
    nxm, nym, nzm = int(tokens[0]), int(tokens[1]), int(tokens[2])
    Lx, Ly, Lz = float(tokens[3]), float(tokens[4]), float(tokens[5])

    data = np.loadtxt(grid_path)
    ym = data[:, 1]  # column 1 (0-based): cell-centre y

    xm = (np.arange(nxm) + 0.5) * (Lx / nxm)
    zm = (np.arange(nzm) + 0.5) * (Lz / nzm)
    return xm, ym, zm


class SDF:
    """A cell-centre signed-distance field: phi[nxm,nym,nzm], phi<0 solid, plus
    the (xm, ym, zm) coordinates it was sampled on."""

    def __init__(self, phi, xm, ym, zm):
        self.phi, self.xm, self.ym, self.zm = phi, xm, ym, zm

    @classmethod
    def read(cls, case=".", sdf_file="SDF_in", params="input_parameters", grid_file=None):
        """Read SDF_in from a case directory (or an explicit sdf_file/params pair).

        `case` is the case directory; `sdf_file`/`params`/`grid_file` are relative
        to it unless given as absolute paths. grid_file defaults to
        <case>/fields/grid.out.
        """
        case = Path(case)
        sdf_path = Path(sdf_file)
        sdf_path = sdf_path if sdf_path.is_absolute() else case / sdf_path
        params_path = Path(params)
        params_path = params_path if params_path.is_absolute() else case / params_path

        p = parse_input_parameters(params_path)
        nxm = p["nx"] - 1
        nym = p["ny"] - 1
        nzm = p["nz"] - 1
        nxg = nxm + 2  # with ghost layers
        nyg = nym + 2

        # Raw array: (nxg, nyg, nzm) big-endian float64, Fortran column-major.
        # order='F' maps flat memory directly to (nxg, nyg, nzm) Fortran indices,
        # so phi_full[ix, iy, iz] = SDF at solver cell (ix, iy, iz).
        raw = np.fromfile(sdf_path, dtype=">f8")
        expected = nxg * nyg * nzm
        if raw.size != expected:
            raise ValueError(
                f"{sdf_path}: SDF size mismatch: got {raw.size} elements, "
                f"expected {expected} ({nxg}x{nyg}x{nzm})"
            )
        phi_full = raw.reshape((nxg, nyg, nzm), order="F")

        # Strip ghost layers in x and y.
        phi = phi_full[1:-1, 1:-1, :]  # (nxm, nym, nzm)

        gf = Path(grid_file) if grid_file else case / "fields" / "grid.out"
        if not gf.is_file():
            raise FileNotFoundError(f"{gf} not found -- needed for grid coordinates.")
        xm, ym, zm = coords_from_grid(gf)

        return cls(phi, xm, ym, zm)

    # -- slicing --
    def slice_xz(self, y):
        """(x_actual, z, phi[nxm,nzm]) at the y index nearest `y`."""
        iy = _nearest_index(self.ym, y)
        return self.ym[iy], self.phi[:, iy, :]

    def slice_xy(self, z):
        """(z_actual, phi[nxm,nym]) at the z index nearest `z`."""
        iz = _nearest_index(self.zm, z)
        return self.zm[iz], self.phi[:, :, iz]

    # -- plotting --
    def _pcolor(self, ax, H, V, data, title, xlabel, ylabel, cmap="RdBu_r", symmetric=True):
        # GenSDF leaves cells outside its geometry-focused AABB at the raw background
        # sentinel (e.g. 1e10 -- "far fluid, never computed"); including those in the
        # colour scale crushes the real near-wall SDF into invisibility, so they are
        # masked out (and shown as a flat grey) rather than colour-scaled.
        from mpl_toolkits.axes_grid1 import make_axes_locatable

        scalarvalue = np.abs(data).max()
        is_sentinel = np.abs(data) >= SENTINEL_FRAC * scalarvalue
        finite = data[~is_sentinel]
        masked = np.ma.masked_where(is_sentinel, data)

        if finite.size:
            vmax = np.abs(finite).max() or 1.0
        else:
            vmax = scalarvalue or 1.0
        vmin = -vmax if symmetric else (finite.min() if finite.size else data.min())

        ax.set_facecolor("lightgray")  # shows through masked (sentinel/uncomputed) cells
        pcm = ax.pcolormesh(H, V, masked.T, shading="auto", cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        div = make_axes_locatable(ax)
        cax = div.append_axes("right", size="3%", pad=0.06)
        ax.get_figure().colorbar(pcm, cax=cax, label="phi  [m]")
        return pcm

    def plot_xz(self, y=0.05, contour=True):
        """X-Z plane at y nearest `y` (horizontal cut through the rough-wall layer).
        Returns the Figure."""
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        y_actual, data = self.slice_xz(y)
        iy = _nearest_index(self.ym, y)
        fig, ax = plt.subplots(figsize=(12, 5))
        self._pcolor(ax, self.xm, self.zm, data,
                     title=f"SDF -- X-Z plane  (y = {y_actual:.4f},  iy = {iy})",
                     xlabel="x", ylabel="z")
        if contour:
            ax.contour(self.xm, self.zm, data.T, levels=[0.0], colors="k", linewidths=0.8)
        fig.tight_layout()
        return fig

    def plot_xy(self, z=3.0, contour=True):
        """X-Y plane at z nearest `z` (streamwise-wall-normal cross-section).
        Returns the Figure."""
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        z_actual, data = self.slice_xy(z)
        iz = _nearest_index(self.zm, z)
        fig, ax = plt.subplots(figsize=(12, 4))
        self._pcolor(ax, self.xm, self.ym, data,
                     title=f"SDF -- X-Y plane  (z = {z_actual:.4f},  iz = {iz})",
                     xlabel="x", ylabel="y  (wall-normal)")
        if contour:
            ax.contour(self.xm, self.ym, data.T, levels=[0.0], colors="k", linewidths=0.8)
        fig.tight_layout()
        return fig
