"""
inflow.py -- SEM/ESEM inflow tooling: three related but distinct tools operating on
three different file formats, kept as separate focused functions/classes in one
module rather than merged into one class.

    from dopamine_post.inflow import mirror_half_channel, InflowDonor, InflowOptState

    # 1. half-channel -> full-channel SEM inflow_profile_file mirroring
    mirror_half_channel("half_channel.csv", "full_channel.csv")

    # 2. donor-plane sanity check (dopamine-ESEM / probe_output.f90 slice)
    donor = InflowDonor.read("inflow_data/inflow_planes")
    mean, rey = donor.time_span_stats()
    fig = donor.plot(reference="full_inflow.csv", grid="fields/grid.out")

    # 3. SEM inflow-optimization restart file
    state = InflowOptState.read("fields/inflow_opt_data.dat")
    print(state.phase_name)
"""
import re
import struct
import sys
import warnings
from pathlib import Path

import numpy as np


# ════════════════════════════════════════════════════════════════════════════
# 1. mirror_half_channel -- half-channel -> full-channel SEM profile mirroring
# ════════════════════════════════════════════════════════════════════════════
#
# Input/output format matches src/sem.f90's inflow_profile_file reader:
# free-form text, '#'-comment/blank lines skipped, data rows
#
#     y  U  V  W  uu  vv  ww  uv  uw  vw
#
# The solver only reads columns 1,2,5,6,7,8 (y,U,uu,vv,ww,uv); the other columns
# (V,W,uw,vw) are carried through/mirrored too for a complete, physically
# consistent file.
#
# Under the reflection y -> 2h-y about the centreline, V (and any correlation
# linear in v: uv, vw) flips sign since v itself changes sign under the
# reflection while u,w don't; U,W,uu,vv,ww,uw stay the same (EVEN).

def read_profile(path):
    """Parse a half/full-channel SEM inflow_profile_file. Returns (header_lines,
    rows) where rows is a list of [y,U,V,W,uu,vv,ww,uv,uw,vw], sorted by y."""
    rows = []
    header_lines = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                header_lines.append(line.rstrip("\n"))
                continue
            parts = s.split()
            if len(parts) < 8:
                raise ValueError(f"line has fewer than 8 columns: {line!r}")
            # pad missing uw,vw with 0 if only 8 columns given
            while len(parts) < 10:
                parts.append("0.0")
            rows.append([float(p) for p in parts[:10]])
    if len(rows) < 2:
        raise ValueError("fewer than 2 data rows found")
    rows.sort(key=lambda r: r[0])
    return header_lines, rows


def _mirror_rows(rows, h):
    """Build the full-channel row list: originals (y<=h) plus the mirror image
    (y -> 2h-y) of every row with y < h (the y=h row, if present, maps to itself
    and must not be duplicated)."""
    out = []
    for y, U, V, W, uu, vv, ww, uv, uw, vw in rows:
        if y > h + 1e-10:
            raise ValueError(
                f"input row at y={y} exceeds the given half-height h={h} -- pass "
                "h explicitly if the profile's own max(y) isn't h"
            )
        out.append([y, U, V, W, uu, vv, ww, uv, uw, vw])
        if y < h - 1e-10:
            neg = lambda x: -x if x != 0.0 else 0.0  # avoid printing "-0.000000e+00"
            out.append([2 * h - y, U, neg(V), W, uu, vv, ww, neg(uv), uw, neg(vw)])
    out.sort(key=lambda r: r[0])
    return out


def mirror_half_channel(input_path, output_path, h=None):
    """Mirror a half-channel SEM inflow profile (y in [0,h], one no-slip wall +
    symmetry/free-slip plane at the centreline) into a full-channel profile
    (y in [0,2h], no-slip at both walls). `h` defaults to max(y) in the input
    file. Writes `output_path` and returns the full-channel row list."""
    header_lines, rows = read_profile(input_path)
    h = h if h is not None else rows[-1][0]

    full_rows = _mirror_rows(rows, h)

    with open(output_path, "w") as f:
        f.write(f"# Full-channel profile mirrored from {input_path} about y=h={h:.6g}\n")
        f.write("# U,W,uu,vv,ww,uw mirrored even; V,uv,vw mirrored odd (negated)\n")
        for hl in header_lines:
            f.write(hl + "\n")
        f.write("# y            U             V             W             "
                "uu            vv            ww            uv            uw            vw\n")
        for y, U, V, W, uu, vv, ww, uv, uw, vw in full_rows:
            f.write(f"{y:.6e}  {U:.6e}  {V:.6e}  {W:.6e}  "
                    f"{uu:.6e}  {vv:.6e}  {ww:.6e}  {uv:.6e}  {uw:.6e}  {vw:.6e}\n")

    print(f"wrote {len(full_rows)} rows spanning y=[0, {2*h:.6g}] to {output_path}")
    print(f"  (from {len(rows)} input rows spanning y=[0, {h:.6g}])")
    return full_rows


# ════════════════════════════════════════════════════════════════════════════
# shared: cell-centre coordinates from fields/grid.out + geometry.out
# (ported from snapshot_io.coords_from_grid, which check_inflow_donor.py imported)
# ════════════════════════════════════════════════════════════════════════════
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


# ════════════════════════════════════════════════════════════════════════════
# 2. InflowDonor -- SEM/ESEM inflow donor-plane check
# ════════════════════════════════════════════════════════════════════════════
#
# Reads the donor slice written either by dopamine-ESEM (the standalone precursor
# generator, src/dopamine_esem_main.f90) or by the main solver's probe_output.f90
# (same binary format: a <base>.bin data stream, <base>_times.bin sample times,
# and a <base>_meta.txt header) -- i.e. exactly the file inflow_type=2 later reads
# back via inflow_recycle_file. Averages the injected mean flow and Reynolds
# stresses over time and the spanwise direction (dir=x donor: axes are (y,z), so
# this is a wall-normal profile) and overlays them on a reference target profile,
# so a bad donor (wrong normalisation, insufficient sem_ensemble_samples, wrong
# profile file, ...) is caught *before* spending a full CFD run on it.
#
# Binary layout (big-endian float64, Fortran stream, no record markers) -- matches
# probe_output.f90's write_slice_n exactly:
#
#     <base>.bin        : nsnaps blocks of (ncomp, n1, n2) Fortran/column-major
#                          float64, one block per sample
#     <base>_times.bin   : nsnaps float64 sample times
#     <base>_meta.txt     : text key = value header (ncomp, n1, n2, dir, pos,
#                          comps, nsnaps, times)

def read_donor_meta(meta_path):
    """Parse a dopamine-ESEM-style '<key>  = <value>' text header."""
    meta = {}
    for line in Path(meta_path).read_text().splitlines():
        if "=" not in line:
            continue
        key, val = line.split("=", 1)
        meta[key.strip()] = val.strip()
    meta["ncomp"] = int(meta["ncomp"])
    meta["n1"] = int(meta["n1"])
    meta["n2"] = int(meta["n2"])
    if "n1_V" not in meta:
        raise ValueError(f"{meta_path}: missing n1_V -- this looks like a donor from an older "
                          f"(pre native-staggered-grid) build; regenerate it with the current dopamine-ESEM")
    meta["n1_V"] = int(meta["n1_V"])
    meta["nsnaps"] = int(meta["nsnaps"])
    meta["pos"] = float(meta["pos"])
    return meta


class InflowDonor:
    """A SEM/ESEM inflow donor slice: per-snapshot U/W on an (n1,n2) cell-centre
    grid and V on its own (n1_V,n2) y-face grid, plus the sample times and the
    parsed meta header.

    Build with `InflowDonor.read(base)`; use `.time_span_stats()` for the
    time-/spanwise-averaged mean and Reynolds stresses, and `.plot()` for the
    six-panel comparison-vs-reference figure that check_inflow_donor.py produced.
    """

    def __init__(self, t, U, V, W, meta):
        self.t, self.U, self.V, self.W, self.meta = t, U, V, W, meta

    @classmethod
    def read(cls, base):
        """Read a donor slice given its basename (no .bin/_meta.txt suffix), as
        written by dopamine-ESEM's per-component native-staggered-grid writer:
        each snapshot is a sequence of separate per-component blocks (fixed
        U,V,W[,T][,C] order), U/W/T/C sharing one (n1,n2) cell-centre grid and V
        on its own (n1_V,n2) y-face grid (n1_V = n1+1) -- no shared grid, no
        interpolation on write.
        """
        base = Path(base)
        meta = read_donor_meta(base.parent / f"{base.name}_meta.txt")

        times_path = base.parent / meta["times"]
        if not times_path.is_file():
            times_path = Path(meta["times"])  # meta stores it relative to cwd at write time
        t = np.fromfile(times_path, dtype=">f8")

        n1, n1v, n2, nsnaps = meta["n1"], meta["n1_V"], meta["n2"], meta["nsnaps"]
        raw = np.fromfile(f"{base}.bin", dtype=">f8")
        expect = nsnaps * (n1 * n2 + n1v * n2 + n1 * n2)  # U + V + W blocks per snapshot
        if raw.size != expect:
            raise ValueError(f"{base}.bin holds {raw.size} floats, expected {expect} "
                              f"(n1={n1}, n1_V={n1v}, n2={n2}, nsnaps={nsnaps})")
        if t.size != nsnaps:
            print(f"WARNING: {times_path} holds {t.size} samples, meta says nsnaps={nsnaps}", file=sys.stderr)

        frame_floats = n1 * n2 + n1v * n2 + n1 * n2
        raw = raw.reshape(nsnaps, frame_floats)
        off = 0
        U = raw[:, off:off + n1 * n2].reshape(nsnaps, n2, n1).transpose(0, 2, 1); off += n1 * n2
        V = raw[:, off:off + n1v * n2].reshape(nsnaps, n2, n1v).transpose(0, 2, 1); off += n1v * n2
        W = raw[:, off:off + n1 * n2].reshape(nsnaps, n2, n1).transpose(0, 2, 1)

        return cls(t, U, V, W, meta)

    def time_span_stats(self, step_start=None, step_end=None):
        """Time- and spanwise- (n2 axis) averaged mean and resolved Reynolds
        stresses, as a function of the n1 axis (wall-normal for a dir='x' donor).

        Returns (mean, rey): mean has columns [U,V,W], rey has columns
        [uu,vv,ww,uv,uw,vw] (uw,vw are zero -- this donor format has no cross
        term between the spanwise-separated component and V's own grid).
        """
        t = self.t

        sel = np.ones(t.size, dtype=bool)
        if step_start is not None:
            sel &= (np.arange(t.size) >= step_start)
        if step_end is not None:
            sel &= (np.arange(t.size) <= step_end)

        U = self.U[sel]  # (nsel, n1,   n2)
        W = self.W[sel]  # (nsel, n1,   n2)
        Vf = self.V[sel]  # (nsel, n1_V, n2) -- on its own y-face grid
        # V interpolated onto U/W's cell-centre grid (post-hoc, for this diagnostic
        # plot only -- the solver itself never does this, it injects V on its own
        # native y-face grid exactly, see sem.f90's recycle_value)
        V = 0.5 * (Vf[:, :-1, :] + Vf[:, 1:, :])  # (nsel, n1, n2)

        meanU = U.mean(axis=(0, 2)); meanV = V.mean(axis=(0, 2)); meanW = W.mean(axis=(0, 2))
        fu = U - meanU[None, :, None]
        fv = V - meanV[None, :, None]
        fw = W - meanW[None, :, None]

        rey = np.stack([
            (fu*fu).mean(axis=(0, 2)), (fv*fv).mean(axis=(0, 2)), (fw*fw).mean(axis=(0, 2)),
            (fu*fv).mean(axis=(0, 2)), (fu*fw).mean(axis=(0, 2)), (fv*fw).mean(axis=(0, 2)),
        ], axis=1)  # (n1,6): uu,vv,ww,uv,uw,vw

        mean_out = np.stack([meanU, meanV, meanW], axis=1)  # (n1,3)

        return mean_out, rey

    def y_coords(self, grid="fields/grid.out"):
        """Wall-normal cell-centre coordinates for the n1 axis, from grid_path;
        falls back to a raw index axis (with a warning) if unavailable or of the
        wrong size."""
        n1 = self.meta["n1"]
        grid_path = Path(grid)
        if not grid_path.is_file():
            warnings.warn(f"{grid_path} not found -- using a raw index axis, not physical y")
            return np.arange(n1)
        try:
            _, y, _ = coords_from_grid(grid_path)
        except Exception as exc:
            warnings.warn(f"could not read {grid_path} ({exc}) -- using a raw index axis")
            return np.arange(n1)
        if y.size != n1:
            warnings.warn(f"{grid_path} has {y.size} cell centres, donor n1={n1} -- "
                           "using a raw index axis")
            return np.arange(n1)
        return y

    def plot(self, reference=None, grid="fields/grid.out", step_start=None, step_end=None,
              output="inflow_donor_check.png", csv=None):
        """Six-panel figure of mean U, -uv, TKE, uu, vv, ww vs. wall-normal
        position, optionally overlaid with a reference profile (a y U V W uu vv ww
        uv [uw vw] text file, the format read/written by sem.f90's
        read_mean_profile). Saves to `output`, prints an error summary vs.
        `reference` if given, and optionally writes the extracted profile to
        `csv`. Returns the Figure."""
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        meta = self.meta
        mean, rey = self.time_span_stats(step_start, step_end)
        n_used = meta["nsnaps"] if step_end is None else min(meta["nsnaps"], step_end + 1)
        n_used -= (step_start or 0)
        print(f"averaged {n_used} of {meta['nsnaps']} samples over the n2={meta['n2']} spanwise points "
              f"({n_used * meta['n2']} samples/level)")

        if meta["dir"] != "x":
            print(f"WARNING: donor dir='{meta['dir']}' -- n1/n2 axis meaning may not be (y,z); "
                  f"check probe_output.f90's convention before trusting the y axis below", file=sys.stderr)

        y = self.y_coords(grid)
        ref = load_reference(reference) if reference else None
        if reference and ref is None:
            print(f"WARNING: could not load reference profile {reference}", file=sys.stderr)

        report_donor_errors(y, mean, rey, ref)

        fig, ax = plt.subplots(2, 3, figsize=(15, 8.5))

        ax[0, 0].plot(y, mean[:, 0], color="C0", label="donor")
        ax[0, 1].plot(y, -rey[:, 3], color="C0", label="donor")
        tke = 0.5 * (rey[:, 0] + rey[:, 1] + rey[:, 2])
        ax[0, 2].plot(y, tke, color="C0", label="donor")

        ax[1, 0].plot(y, rey[:, 0], color="C0", label="donor")
        ax[1, 1].plot(y, rey[:, 1], color="C0", label="donor")
        ax[1, 2].plot(y, rey[:, 2], color="C0", label="donor")

        if ref is not None:
            ax[0, 0].plot(ref["y"], ref["U"], "k--", lw=1.6, label="reference")
            ax[0, 1].plot(ref["y"], -ref["uv"], "k--", lw=1.6, label="reference")
            tke_ref = 0.5 * (ref["uu"] + ref["vv"] + ref["ww"])
            ax[0, 2].plot(ref["y"], tke_ref, "k--", lw=1.6, label="reference")
            ax[1, 0].plot(ref["y"], ref["uu"], "k--", lw=1.6, label="reference")
            ax[1, 1].plot(ref["y"], ref["vv"], "k--", lw=1.6, label="reference")
            ax[1, 2].plot(ref["y"], ref["ww"], "k--", lw=1.6, label="reference")

        titles = [(r"mean $U$", r"$U$"),
                  (r"shear stress $-\overline{u'v'}$", r"$-\overline{u'v'}$"),
                  (r"TKE $\frac{1}{2}\overline{u_i'u_i'}$", r"TKE"),
                  (r"$\overline{u'u'}$ (streamwise)", r"$\overline{u'u'}$"),
                  (r"$\overline{v'v'}$ (wall-normal)", r"$\overline{v'v'}$"),
                  (r"$\overline{w'w'}$ (spanwise)", r"$\overline{w'w'}$")]
        for a, (t, yl) in zip(ax.ravel(), titles):
            a.set_title(t)
            a.set_xlabel(r"$y$")
            a.set_ylabel(yl)
            a.grid(alpha=0.25)
            a.legend(fontsize=8, frameon=False)

        fig.suptitle(f"inflow donor check  (x = {meta['pos']:.4f}, "
                     f"samples {step_start if step_start is not None else 'first'}"
                     f"..{step_end if step_end is not None else 'last'})", y=0.99)
        fig.tight_layout()
        fig.savefig(output, dpi=150)
        print(f"wrote {output}")

        if csv:
            header = "y,U,V,W,uu,vv,ww,uv,uw,vw"
            out = np.column_stack([y, mean, rey])
            np.savetxt(csv, out, delimiter=",", header=header, comments="")
            print(f"wrote {csv}")

        return fig


def load_reference(path):
    """y U V W uu vv ww uv [uw vw], whitespace-delimited, '#' comments -- the
    format read/written by sem.f90's read_mean_profile (e.g. full_inflow.csv)."""
    prof = np.loadtxt(path, comments="#")
    if prof.shape[1] < 8:
        raise ValueError(f"{path}: expected >=8 columns (y U V W uu vv ww uv), got {prof.shape[1]}")
    out = dict(y=prof[:, 0], U=prof[:, 1], V=prof[:, 2], W=prof[:, 3],
               uu=prof[:, 4], vv=prof[:, 5], ww=prof[:, 6], uv=prof[:, 7])
    if prof.shape[1] >= 10:
        out["uw"] = prof[:, 8]
        out["vw"] = prof[:, 9]
    return out


def report_donor_errors(y, mean, rey, ref):
    """Print a quick relative-error summary against the reference profile,
    interpolated onto y. No-op if ref is None."""
    if ref is None:
        return
    fields = {"U": mean[:, 0], "uu": rey[:, 0], "vv": rey[:, 1], "ww": rey[:, 2], "-uv": -rey[:, 3]}
    ref_fields = {"U": ref["U"], "uu": ref["uu"], "vv": ref["vv"], "ww": ref["ww"], "-uv": -ref["uv"]}

    print(f"\n{'field':>6}  {'peak(donor)':>12}  {'peak(ref)':>12}  {'ratio':>8}  {'L2 rel.err':>10}")
    for name, sim in fields.items():
        r = np.interp(y, ref["y"], ref_fields[name])
        peak_i = np.argmax(np.abs(r))
        ratio = sim[peak_i] / r[peak_i] if r[peak_i] != 0 else np.nan
        denom = np.sqrt(np.mean(r ** 2))
        l2 = np.sqrt(np.mean((sim - r) ** 2)) / denom if denom > 0 else np.nan
        print(f"{name:>6}  {sim[peak_i]:12.4f}  {r[peak_i]:12.4f}  {ratio:8.3f}  {l2:10.3f}")


def check_donor(donor, reference=None, grid="fields/grid.out", step_start=None, step_end=None,
                 output="inflow_donor_check.png", csv=None):
    """Functional convenience wrapper: `InflowDonor.read(donor)` then `.plot(...)`."""
    return InflowDonor.read(donor).plot(reference=reference, grid=grid, step_start=step_start,
                                         step_end=step_end, output=output, csv=csv)


# ════════════════════════════════════════════════════════════════════════════
# 3. InflowOptState -- SEM inflow-optimization restart file
# ════════════════════════════════════════════════════════════════════════════
#
# Reads fields/inflow_opt_data.dat, the SEM inflow-optimization restart file
# written by write_inflow_opt_restart / read_inflow_opt_restart in src/sem.f90.
#
# Binary layout (Fortran unformatted stream, no record markers, big-endian --
# this codebase's CMakeLists.txt builds with -fconvert=big-endian / convert
# big_endian / -Mbyteswapio), all Real values 8-byte (Real(Int64) in the source is
# double precision), all Integer values 4-byte:
#
#     n_bezier                                   int32
#     inflow_opt_phase, inflow_opt_step_count,   4 x int32
#       inflow_opt_iter, inflow_opt_no_improve
#     x_cp_R22(n_bezier), x_cp_R33(n_bezier)     2 x n_bezier float64   (frozen/current control points)
#     x_prev_R22(n_bezier), x_prev_R33(n_bezier) 2 x n_bezier float64   (control points before the last-applied correction)
#     stats_step0(3,n_bezier)                    3 x n_bezier float64  (u'^2, v'^2, w'^2 at baseline)
#     stats_step1(3,n_bezier)                    3 x n_bezier float64  (measured after doubling v'^2 and w'^2 together)
#     stats_prev(3,n_bezier)                     3 x n_bezier float64  (measured stats before the last-applied correction)
#     slope_v(n_bezier), slope_w(n_bezier)       2 x n_bezier float64  (per-control-point scalar secant slopes)
#     best_resid                                 float64               (worst-case relative residual at the best iterate)
#     best_x_cp_R22(n_bezier), best_x_cp_R33(n_bezier)  2 x n_bezier float64 (control points at the best iterate)
#     prof_R22(n_profile), prof_R33(n_profile)   2 x n_profile float64 (frozen full-resolution profile)
#
# n_profile is not stored explicitly; it is recovered from the remaining file size.

PHASE_NAME = {
    1: "1 (accumulating step0 / baseline)",
    2: "2 (accumulating step1 / v'^2 and w'^2 doubled together)",
    3: "3 (verify -- the paper's single correction, or the experimental iter>1 extension)",
    4: "4 (frozen -- best iterate applied)",
}


class InflowOptState:
    """Parsed SEM inflow-optimization restart file (fields/inflow_opt_data.dat).

    Build with `InflowOptState.read(path)`. Fields mirror the Fortran restart
    layout 1:1 (see module docstring); `.phase_name` is the human-readable
    PHASE_NAME string for `.phase`.
    """

    def __init__(self, **kw):
        for k, v in kw.items():
            setattr(self, k, v)

    @property
    def phase_name(self):
        return PHASE_NAME.get(self.phase, self.phase)

    @classmethod
    def read(cls, path):
        data = Path(path).read_bytes()
        off = 0

        def take(fmt, n=1):
            nonlocal off
            size = struct.calcsize(fmt) * n
            # CMakeLists.txt builds with -fconvert=big-endian (gfortran) / convert
            # big_endian (ifort) / -Mbyteswapio (nvfortran), so all unformatted
            # stream I/O in this codebase is big-endian regardless of host arch
            vals = struct.unpack_from(f">{n}{fmt}", data, off)
            off += size
            return vals if n > 1 else vals[0]

        n_bezier = take("i")
        phase, step_count, opt_iter, no_improve = take("i", 4)
        x_cp_R22 = take("d", n_bezier)
        x_cp_R33 = take("d", n_bezier)
        x_prev_R22 = take("d", n_bezier)
        x_prev_R33 = take("d", n_bezier)
        stats_step0 = take("d", 3 * n_bezier)
        stats_step1 = take("d", 3 * n_bezier)
        stats_prev = take("d", 3 * n_bezier)
        slope_v = take("d", n_bezier)
        slope_w = take("d", n_bezier)
        best_resid = take("d")
        best_x_cp_R22 = take("d", n_bezier)
        best_x_cp_R33 = take("d", n_bezier)

        remaining = len(data) - off
        if remaining % (2 * 8) != 0:
            raise ValueError(
                f"unexpected trailing byte count ({remaining}); file may not match "
                "the current write_inflow_opt_restart layout"
            )
        n_profile = remaining // 16
        prof_R22 = take("d", n_profile)
        prof_R33 = take("d", n_profile)

        def reshape3(flat):
            # Fortran column-major (3,n_bezier): fastest-varying index is the row (u/v/w)
            return [flat[i::3] for i in range(3)]

        return cls(
            n_bezier=n_bezier, phase=phase, step_count=step_count, opt_iter=opt_iter,
            no_improve=no_improve, x_cp_R22=x_cp_R22, x_cp_R33=x_cp_R33,
            x_prev_R22=x_prev_R22, x_prev_R33=x_prev_R33,
            stats_step0=reshape3(stats_step0), stats_step1=reshape3(stats_step1),
            stats_prev=reshape3(stats_prev), slope_v=slope_v, slope_w=slope_w,
            best_resid=best_resid, best_x_cp_R22=best_x_cp_R22, best_x_cp_R33=best_x_cp_R33,
            n_profile=n_profile, prof_R22=prof_R22, prof_R33=prof_R33,
        )

    def summary(self):
        """Print the same human-readable summary read_inflow_opt.py's CLI did."""
        print(f"n_bezier:      {self.n_bezier}")
        print(f"n_profile:     {self.n_profile}")
        print(f"phase:         {self.phase_name}")
        print(f"step_count:    {self.step_count} (steps accumulated in the current phase so far)")
        print(f"opt_iter:      {self.opt_iter} (correction steps applied so far)")
        print(f"no_improve:    {self.no_improve} (consecutive non-improving corrections)")
        print(f"best_resid:    {self.best_resid:.5f} (worst-case relative residual at the frozen/best iterate)")

    def control_point_heights(self, profile_path):
        """Bezier control-point y positions spanning the y range of `profile_path`
        (a y U V W uu vv ww uv text profile), matching sem.f90's control-point
        placement (quadratic clustering toward the wall)."""
        prof_y, _, _ = read_target_profile(profile_path, fmt=0)
        return bezier_heights(prof_y[0], prof_y[-1], self.n_bezier)


def namelist_value(text, group, key):
    """A scalar value from one Fortran namelist group in `text`, e.g.
    namelist_value(text, "INFLOW", "sem_profile_format")."""
    # namelist terminator '/' stands alone on its own line; a bare '.*?/' would
    # instead stop at the first '/' inside a quoted path value (e.g. a filename)
    m = re.search(rf"&{group}\b(.*?)^\s*/\s*$", text, re.S | re.I | re.M)
    if not m:
        return None
    m2 = re.search(rf"\b{key}\s*=\s*'?([^,'\n]+?)'?\s*(?:,|$)", m.group(1), re.I | re.M)
    return m2.group(1).strip() if m2 else None


def read_target_profile(path, fmt):
    """The target profile referenced by inflow_profile_file/sem_profile_format,
    reduced to (y, R22, R33). fmt=1: wind-tunnel TI format (z U Iu Iv Iw [...]);
    fmt=0: Reynolds-stress format (y U V W uu vv ww uv)."""
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            rows.append([float(x) for x in line.split()])

    y = [r[0] for r in rows]
    U = [r[1] for r in rows]
    if fmt == 1:
        R22 = [(r[3] * u) ** 2 for r, u in zip(rows, U)]
        R33 = [(r[4] * u) ** 2 for r, u in zip(rows, U)]
    else:
        R22 = [r[5] for r in rows]
        R33 = [r[6] for r in rows]
    return y, R22, R33


def linterp(xs, ys, xq):
    if xq <= xs[0]:
        return ys[0]
    if xq >= xs[-1]:
        return ys[-1]
    for i in range(1, len(xs)):
        if xq <= xs[i]:
            t = (xq - xs[i - 1]) / (xs[i] - xs[i - 1])
            return ys[i - 1] + t * (ys[i] - ys[i - 1])
    return ys[-1]


def bezier_heights(y0, y1, n_bezier):
    return [y0 + (y1 - y0) * (i / (n_bezier - 1)) ** 2 for i in range(n_bezier)]


def read_inflow_opt(case_dir, restart="fields/inflow_opt_data.dat", profile=None, fmt=None):
    """Read an inflow-optimization restart file plus (optionally auto-detected
    from `<case_dir>/input_parameters`) its target profile. Returns
    (state, y_cp, target_v, target_w) where the latter three are None if no
    profile could be resolved. This is `InflowOptState.read` plus the
    profile/control-point resolution read_inflow_opt.py's CLI did before
    printing/plotting."""
    case_dir = Path(case_dir)
    restart_path = Path(restart)
    if not restart_path.is_absolute():
        restart_path = case_dir / restart_path

    state = InflowOptState.read(restart_path)

    profile_path = profile
    input_parameters = case_dir / "input_parameters"
    if (profile_path is None or fmt is None) and input_parameters.exists():
        text = input_parameters.read_text()
        if profile_path is None:
            profile_path = namelist_value(text, "INFLOW", "inflow_profile_file")
            if profile_path:
                profile_path = case_dir / profile_path
        if fmt is None:
            fmt_str = namelist_value(text, "INFLOW", "sem_profile_format")
            fmt = int(float(fmt_str)) if fmt_str else 0

    y_cp = target_v = target_w = None
    if profile_path and Path(profile_path).exists():
        prof_y, prof_R22_t, prof_R33_t = read_target_profile(profile_path, fmt or 0)
        y_cp = bezier_heights(prof_y[0], prof_y[-1], state.n_bezier)
        target_v = [linterp(prof_y, prof_R22_t, y) for y in y_cp]
        target_w = [linterp(prof_y, prof_R33_t, y) for y in y_cp]

    return state, y_cp, target_v, target_w


def plot_inflow_opt(state, y_cp=None, profile=None, fmt=None, output=None):
    """The v'^2/w'^2 vs. y comparison figure read_inflow_opt.py's --plot produced:
    target profile (if `profile` given) vs. the frozen full-resolution injected
    profile. Returns the Figure; saves to `output` if given."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), sharey=True)
    if profile and Path(profile).exists():
        prof_y, prof_R22_t, prof_R33_t = read_target_profile(profile, fmt or 0)
        axes[0].plot(prof_R22_t, prof_y, "k--", label="target (wind tunnel)")
        axes[1].plot(prof_R33_t, prof_y, "k--", label="target (wind tunnel)")
        axes[0].plot(state.prof_R22, prof_y, "C0-o", label="frozen inflow")
        axes[1].plot(state.prof_R33, prof_y, "C1-o", label="frozen inflow")
    else:
        idx = list(range(state.n_profile))
        axes[0].plot(state.prof_R22, idx, "C0-o", label="frozen inflow")
        axes[1].plot(state.prof_R33, idx, "C1-o", label="frozen inflow")
    axes[0].set_xlabel("v'^2"); axes[0].set_ylabel("y")
    axes[1].set_xlabel("w'^2")
    for ax in axes:
        ax.legend(); ax.grid(True, alpha=0.3)
    fig.suptitle(f"Inflow optimization -- phase {state.phase}")
    fig.tight_layout()
    if output:
        fig.savefig(output, dpi=150)
        print(f"saved plot to {output}")
    return fig
