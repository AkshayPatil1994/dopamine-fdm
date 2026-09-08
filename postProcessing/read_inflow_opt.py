#!/usr/bin/env python3
"""Read/plot fields/inflow_opt_data.dat, the SEM inflow-optimization restart file
written by write_inflow_opt_restart / read_inflow_opt_restart in src/sem.f90.

Binary layout (Fortran unformatted stream, no record markers, big-endian --
this codebase's CMakeLists.txt builds with -fconvert=big-endian / convert
big_endian / -Mbyteswapio), all Real values 8-byte (Real(Int64) in the source
is double precision), all Integer values 4-byte:

    n_bezier                                   int32
    inflow_opt_phase, inflow_opt_step_count,   4 x int32
      inflow_opt_iter, inflow_opt_no_improve
    x_cp_R22(n_bezier), x_cp_R33(n_bezier)     2 x n_bezier float64   (frozen/current control points)
    x_prev_R22(n_bezier), x_prev_R33(n_bezier) 2 x n_bezier float64   (control points before the last-applied correction)
    stats_step0(3,n_bezier)                    3 x n_bezier float64  (u'^2, v'^2, w'^2 at baseline)
    stats_step1(3,n_bezier)                    3 x n_bezier float64  (measured after doubling v'^2 and w'^2 together)
    stats_prev(3,n_bezier)                     3 x n_bezier float64  (measured stats before the last-applied correction)
    slope_v(n_bezier), slope_w(n_bezier)       2 x n_bezier float64  (per-control-point scalar secant slopes)
    best_resid                                 float64               (worst-case relative residual at the best iterate)
    best_x_cp_R22(n_bezier), best_x_cp_R33(n_bezier)  2 x n_bezier float64 (control points at the best iterate)
    prof_R22(n_profile), prof_R33(n_profile)   2 x n_profile float64 (frozen full-resolution profile)

n_profile is not stored explicitly; it is recovered from the remaining file size.

Usage:
    python3 read_inflow_opt.py /mnt/storage1/fdm-dopamine/validation/aij/caseC
    python3 read_inflow_opt.py /path/to/case --restart fields/inflow_opt_data.dat \
        --profile referenceData/windtunnel_inflow.csv --format 1
"""

import argparse
import re
import struct
import sys
from pathlib import Path

PHASE_NAME = {
    1: "1 (accumulating step0 / baseline)",
    2: "2 (accumulating step1 / v'^2 and w'^2 doubled together)",
    3: "3 (verify -- the paper's single correction, or the experimental iter>1 extension)",
    4: "4 (frozen -- best iterate applied)",
}


def read_restart(path):
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

    return {
        "n_bezier": n_bezier,
        "phase": phase,
        "step_count": step_count,
        "opt_iter": opt_iter,
        "no_improve": no_improve,
        "x_cp_R22": x_cp_R22,
        "x_cp_R33": x_cp_R33,
        "x_prev_R22": x_prev_R22,
        "x_prev_R33": x_prev_R33,
        "stats_step0": reshape3(stats_step0),
        "stats_step1": reshape3(stats_step1),
        "stats_prev": reshape3(stats_prev),
        "slope_v": slope_v,
        "slope_w": slope_w,
        "best_resid": best_resid,
        "best_x_cp_R22": best_x_cp_R22,
        "best_x_cp_R33": best_x_cp_R33,
        "n_profile": n_profile,
        "prof_R22": prof_R22,
        "prof_R33": prof_R33,
    }


def namelist_value(text, group, key):
    # namelist terminator '/' stands alone on its own line; a bare '.*?/' would
    # instead stop at the first '/' inside a quoted path value (e.g. a filename)
    m = re.search(rf"&{group}\b(.*?)^\s*/\s*$", text, re.S | re.I | re.M)
    if not m:
        return None
    m2 = re.search(rf"\b{key}\s*=\s*'?([^,'\n]+?)'?\s*(?:,|$)", m.group(1), re.I | re.M)
    return m2.group(1).strip() if m2 else None


def read_profile(path, fmt):
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
        # wind-tunnel TI format: z U Iu Iv Iw [...]
        R22 = [(r[3] * u) ** 2 for r, u in zip(rows, U)]
        R33 = [(r[4] * u) ** 2 for r, u in zip(rows, U)]
    else:
        # Reynolds-stress format: y U V W uu vv ww uv
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


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case_dir", help="case directory containing fields/inflow_opt_data.dat")
    ap.add_argument("--restart", default="fields/inflow_opt_data.dat",
                     help="restart file path, relative to case_dir unless absolute")
    ap.add_argument("--profile", default=None,
                     help="inflow_profile_file to overlay targets (default: read from input_parameters)")
    ap.add_argument("--format", type=int, default=None, choices=[0, 1],
                     help="sem_profile_format (default: read from input_parameters)")
    ap.add_argument("--plot", default=None, help="save a comparison PNG to this path (requires matplotlib)")
    args = ap.parse_args()

    case_dir = Path(args.case_dir)
    restart_path = Path(args.restart)
    if not restart_path.is_absolute():
        restart_path = case_dir / restart_path

    d = read_restart(restart_path)

    print(f"file:          {restart_path}")
    print(f"n_bezier:      {d['n_bezier']}")
    print(f"n_profile:     {d['n_profile']}")
    print(f"phase:         {PHASE_NAME.get(d['phase'], d['phase'])}")
    print(f"step_count:    {d['step_count']} (steps accumulated in the current phase so far)")
    print(f"opt_iter:      {d['opt_iter']} (correction steps applied so far)")
    print(f"no_improve:    {d['no_improve']} (consecutive non-improving corrections)")
    print(f"best_resid:    {d['best_resid']:.5f} (worst-case relative residual at the frozen/best iterate)")
    print()

    profile_path = args.profile
    fmt = args.format
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
        prof_y, prof_R22_t, prof_R33_t = read_profile(profile_path, fmt or 0)
        y_cp = bezier_heights(prof_y[0], prof_y[-1], d["n_bezier"])
        target_v = [linterp(prof_y, prof_R22_t, y) for y in y_cp]
        target_w = [linterp(prof_y, prof_R33_t, y) for y in y_cp]
    else:
        print("(no profile file found -- pass --profile/--format to see target-vs-corrected control points)\n")

    header = f"{'idx':>3} {'y_cp':>10} {'target_v2':>12} {'corrected_v2':>13} {'target_w2':>12} {'corrected_w2':>13}"
    if y_cp is not None:
        print("Control-point correction (frozen inflow vs. target):")
        print(header)
        for i in range(d["n_bezier"]):
            print(f"{i+1:>3} {y_cp[i]:>10.4f} {target_v[i]:>12.5f} {d['x_cp_R22'][i]:>13.5f} "
                  f"{target_w[i]:>12.5f} {d['x_cp_R33'][i]:>13.5f}")
        print()

    print("Measured downstream stats at each control point, by phase (rows: u'^2, v'^2, w'^2):")
    for name in ("stats_step0", "stats_step1"):
        print(f"  {name}:")
        for row_name, row in zip(("u'^2", "v'^2", "w'^2"), d[name]):
            print(f"    {row_name}: " + " ".join(f"{v:8.5f}" for v in row))
    print()

    print(f"Frozen full-resolution injected profile (n_profile={d['n_profile']}):")
    print("  prof_R22: " + " ".join(f"{v:8.5f}" for v in d["prof_R22"]))
    print("  prof_R33: " + " ".join(f"{v:8.5f}" for v in d["prof_R33"]))

    if args.plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print("\nmatplotlib not available; skipping --plot", file=sys.stderr)
            return

        fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), sharey=True)
        if profile_path and Path(profile_path).exists():
            prof_y, prof_R22_t, prof_R33_t = read_profile(profile_path, fmt or 0)
            axes[0].plot(prof_R22_t, prof_y, "k--", label="target (wind tunnel)")
            axes[1].plot(prof_R33_t, prof_y, "k--", label="target (wind tunnel)")
            axes[0].plot(d["prof_R22"], prof_y, "C0-o", label="frozen inflow")
            axes[1].plot(d["prof_R33"], prof_y, "C1-o", label="frozen inflow")
        else:
            idx = list(range(d["n_profile"]))
            axes[0].plot(d["prof_R22"], idx, "C0-o", label="frozen inflow")
            axes[1].plot(d["prof_R33"], idx, "C1-o", label="frozen inflow")
        axes[0].set_xlabel("v'^2");  axes[0].set_ylabel("y")
        axes[1].set_xlabel("w'^2")
        for ax in axes:
            ax.legend();  ax.grid(True, alpha=0.3)
        fig.suptitle(f"Inflow optimization -- phase {d['phase']}")
        fig.tight_layout()
        fig.savefig(args.plot, dpi=150)
        print(f"\nsaved plot to {args.plot}")


if __name__ == "__main__":
    main()
