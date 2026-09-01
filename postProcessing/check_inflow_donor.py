#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 DOPAMINE contributors
# SPDX-License-Identifier: AGPL-3.0-only

"""
check_inflow_donor.py — Verify a SEM/ESEM inflow donor plane before running anything.

Reads the donor slice written either by dopamine-ESEM (the standalone precursor
generator, src/dopamine_esem_main.f90) or by the main solver's probe_output.f90
(same binary format: a <base>.bin data stream, <base>_times.bin sample times,
and a <base>_meta.txt header) -- i.e. exactly the file inflow_type=2 later reads
back via inflow_recycle_file. Averages the injected mean flow and Reynolds
stresses over time and the spanwise direction (dir=x donor: axes are (y,z), so
this is a wall-normal profile) and overlays them on a reference target profile,
so a bad donor (wrong normalisation, insufficient sem_ensemble_samples, wrong
profile file, ...) is caught *before* spending a full CFD run on it.

Binary layout (big-endian float64, Fortran stream, no record markers) --
matches probe_output.f90's write_slice_n exactly:

    <base>.bin        : nsnaps blocks of (ncomp, n1, n2) Fortran/column-major
                         float64, one block per sample
    <base>_times.bin   : nsnaps float64 sample times
    <base>_meta.txt     : text key = value header (ncomp, n1, n2, dir, pos,
                         comps, nsnaps, times)

Reference profile: a whitespace-delimited text file with columns
    y  U  V  W  uu  vv  ww  uv  [uw  vw]
(the format read/written by mean_profile/read_mean_profile in sem.f90, e.g.
this case's full_inflow.csv/inflow.csv) -- already in solver units, no rescaling
applied.

Usage
-----
    python3 check_inflow_donor.py --donor inflow_data/inflow_planes \\
        --reference full_inflow.csv --grid fields/grid.out \\
        --output inflow_donor_check.png
"""

import argparse
import sys
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
import snapshot_io


# ----------------------------------------------------------------------------
# donor slice I/O
# ----------------------------------------------------------------------------
def read_meta(meta_path):
    """Parse a dopamine-ESEM-style '<key>  = <value>' text header."""
    meta = {}
    for line in Path(meta_path).read_text().splitlines():
        if "=" not in line:
            continue
        key, val = line.split("=", 1)
        meta[key.strip()] = val.strip()
    meta["ncomp"]  = int(meta["ncomp"])
    meta["n1"]     = int(meta["n1"])
    meta["n2"]     = int(meta["n2"])
    if "n1_V" not in meta:
        raise ValueError(f"{meta_path}: missing n1_V -- this looks like a donor from an older "
                          f"(pre native-staggered-grid) build; regenerate it with the current dopamine-ESEM")
    meta["n1_V"]   = int(meta["n1_V"])
    meta["nsnaps"] = int(meta["nsnaps"])
    meta["pos"]    = float(meta["pos"])
    return meta


def read_donor(base):
    """
    Read a donor slice given its basename (no .bin/_meta.txt suffix), as written
    by dopamine-ESEM's per-component native-staggered-grid writer: each snapshot
    is a sequence of separate per-component blocks (fixed U,V,W[,T][,C] order),
    U/W/T/C sharing one (n1,n2) cell-centre grid and V on its own (n1_V,n2)
    y-face grid (n1_V = n1+1) -- no shared grid, no interpolation on write.

    Returns dict with:
        t     : (nsnaps,) sample times
        U, W  : (nsnaps, n1,   n2) float64
        V     : (nsnaps, n1_V, n2) float64
        meta  : the parsed header dict
    """
    base = Path(base)
    meta = read_meta(base.parent / f"{base.name}_meta.txt")

    times_path = base.parent / meta["times"]
    if not times_path.is_file():
        times_path = Path(meta["times"])  # meta stores it relative to cwd at write time
    t = np.fromfile(times_path, dtype=">f8")

    n1, n1v, n2, nsnaps = meta["n1"], meta["n1_V"], meta["n2"], meta["nsnaps"]
    raw = np.fromfile(f"{base}.bin", dtype=">f8")
    expect = nsnaps * (n1*n2 + n1v*n2 + n1*n2)   # U + V + W blocks per snapshot
    if raw.size != expect:
        raise ValueError(f"{base}.bin holds {raw.size} floats, expected {expect} "
                          f"(n1={n1}, n1_V={n1v}, n2={n2}, nsnaps={nsnaps})")
    if t.size != nsnaps:
        print(f"WARNING: {times_path} holds {t.size} samples, meta says nsnaps={nsnaps}", file=sys.stderr)

    frame_floats = n1*n2 + n1v*n2 + n1*n2
    raw = raw.reshape(nsnaps, frame_floats)
    off = 0
    U = raw[:, off:off+n1*n2].reshape(nsnaps, n2, n1).transpose(0, 2, 1);  off += n1*n2
    V = raw[:, off:off+n1v*n2].reshape(nsnaps, n2, n1v).transpose(0, 2, 1); off += n1v*n2
    W = raw[:, off:off+n1*n2].reshape(nsnaps, n2, n1).transpose(0, 2, 1)

    return dict(t=t, U=U, V=V, W=W, meta=meta)


# ----------------------------------------------------------------------------
# statistics
# ----------------------------------------------------------------------------
def time_span_stats(donor, step_start=None, step_end=None):
    """
    Time- and spanwise- (n2 axis) averaged mean and resolved Reynolds
    stresses, as a function of the n1 axis (wall-normal for a dir='x' donor).

    Returns (mean, rey): mean has columns [U,V,W] (or fewer if ncomp<3),
    rey has columns [uu,vv,ww,uv,uw,vw] (zero where a component is absent).
    """
    t = donor["t"]

    sel = np.ones(t.size, dtype=bool)
    if step_start is not None:
        sel &= (np.arange(t.size) >= step_start)
    if step_end is not None:
        sel &= (np.arange(t.size) <= step_end)

    U = donor["U"][sel]  # (nsel, n1,   n2)
    W = donor["W"][sel]  # (nsel, n1,   n2)
    Vf = donor["V"][sel]  # (nsel, n1_V, n2) -- on its own y-face grid
    # V interpolated onto U/W's cell-centre grid (post-hoc, for this diagnostic
    # plot only -- the solver itself never does this, it injects V on its own
    # native y-face grid exactly, see sem.f90's recycle_value)
    V = 0.5 * (Vf[:, :-1, :] + Vf[:, 1:, :])   # (nsel, n1, n2)

    n1 = U.shape[1]
    meanU = U.mean(axis=(0, 2));  meanV = V.mean(axis=(0, 2));  meanW = W.mean(axis=(0, 2))
    fu = U - meanU[None, :, None]
    fv = V - meanV[None, :, None]
    fw = W - meanW[None, :, None]

    rey = np.stack([
        (fu*fu).mean(axis=(0, 2)), (fv*fv).mean(axis=(0, 2)), (fw*fw).mean(axis=(0, 2)),
        (fu*fv).mean(axis=(0, 2)), (fu*fw).mean(axis=(0, 2)), (fv*fw).mean(axis=(0, 2)),
    ], axis=1)  # (n1,6): uu,vv,ww,uv,uw,vw

    mean_out = np.stack([meanU, meanV, meanW], axis=1)  # (n1,3)

    return mean_out, rey


# ----------------------------------------------------------------------------
# reference profile
# ----------------------------------------------------------------------------
def load_reference(path):
    """
    y U V W uu vv ww uv [uw vw], whitespace-delimited, '#' comments -- the
    format read/written by sem.f90's read_mean_profile (e.g. full_inflow.csv).
    """
    prof = np.loadtxt(path, comments="#")
    if prof.shape[1] < 8:
        raise ValueError(f"{path}: expected >=8 columns (y U V W uu vv ww uv), got {prof.shape[1]}")
    out = dict(y=prof[:, 0], U=prof[:, 1], V=prof[:, 2], W=prof[:, 3],
               uu=prof[:, 4], vv=prof[:, 5], ww=prof[:, 6], uv=prof[:, 7])
    if prof.shape[1] >= 10:
        out["uw"] = prof[:, 8]
        out["vw"] = prof[:, 9]
    return out


# ----------------------------------------------------------------------------
# plotting
# ----------------------------------------------------------------------------
def make_figure(y, mean, rey, ref, args):
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

    fig.suptitle(f"inflow donor check: {args.donor}  (x = {args.pos:.4f}, "
                 f"samples {args.step_start if args.step_start is not None else 'first'}"
                 f"..{args.step_end if args.step_end is not None else 'last'})", y=0.99)
    fig.tight_layout()
    fig.savefig(args.output, dpi=150)
    print(f"wrote {args.output}")


def report_errors(y, mean, rey, ref):
    """Print a quick relative-error summary against the reference profile, interpolated onto y."""
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


# ----------------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Time- and spanwise-averaged statistics of a SEM/ESEM inflow donor plane, "
                    "checked against a reference profile before running any CFD.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--donor", required=True,
                     help="donor basename, e.g. inflow_data/inflow_planes "
                          "(reads <donor>.bin, <donor>_meta.txt, and the times file it names)")
    ap.add_argument("--reference", default=None,
                     help="y U V W uu vv ww uv [uw vw] reference profile, e.g. full_inflow.csv")
    ap.add_argument("--grid", default="fields/grid.out",
                     help="grid.out (with geometry.out alongside) giving the wall-normal cell-centre "
                          "coordinates for a dir='x' donor's n1 axis; if absent, plots use raw index")
    ap.add_argument("--step-start", type=int, default=None, help="discard donor samples before this index")
    ap.add_argument("--step-end", type=int, default=None, help="discard donor samples after this index")
    ap.add_argument("--output", default="inflow_donor_check.png")
    ap.add_argument("--csv", default=None, help="also write the extracted donor profile here")
    args = ap.parse_args()

    donor = read_donor(args.donor)
    meta = donor["meta"]
    print(f"donor: dir={meta['dir']}  pos={meta['pos']:.4f}  n1={meta['n1']}  n2={meta['n2']}  "
          f"nsnaps={meta['nsnaps']}  comps={meta['comps']}")
    args.pos = meta["pos"]

    mean, rey = time_span_stats(donor, args.step_start, args.step_end)
    n_used = meta["nsnaps"] if args.step_end is None else min(meta["nsnaps"], args.step_end + 1)
    n_used -= (args.step_start or 0)
    print(f"averaged {n_used} of {meta['nsnaps']} samples over the n2={meta['n2']} spanwise points "
          f"({n_used * meta['n2']} samples/level)")

    if meta["dir"] != "x":
        print(f"WARNING: donor dir='{meta['dir']}' -- n1/n2 axis meaning may not be (y,z); "
              f"check probe_output.f90's convention before trusting the y axis below", file=sys.stderr)

    grid_path = Path(args.grid)
    if grid_path.is_file():
        try:
            _, y, _ = snapshot_io.coords_from_grid(grid_path)
            if y.size != meta["n1"]:
                print(f"WARNING: {grid_path} has {y.size} cell centres, donor n1={meta['n1']} -- "
                      f"falling back to a raw index axis", file=sys.stderr)
                y = np.arange(meta["n1"])
        except Exception as exc:
            print(f"WARNING: could not read {grid_path} ({exc}) -- falling back to a raw index axis",
                  file=sys.stderr)
            y = np.arange(meta["n1"])
    else:
        print(f"NOTE: {grid_path} not found -- plotting against donor grid index, not physical y",
              file=sys.stderr)
        y = np.arange(meta["n1"])

    ref = load_reference(args.reference) if args.reference else None
    if args.reference and ref is None:
        print(f"WARNING: could not load reference profile {args.reference}", file=sys.stderr)

    report_errors(y, mean, rey, ref)
    make_figure(y, mean, rey, ref, args)

    if args.csv:
        header = "y,U,V,W,uu,vv,ww,uv,uw,vw"
        out = np.column_stack([y, mean, rey])
        np.savetxt(args.csv, out, delimiter=",", header=header, comments="")
        print(f"wrote {args.csv}")


if __name__ == "__main__":
    main()
