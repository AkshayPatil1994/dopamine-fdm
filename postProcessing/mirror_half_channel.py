#!/usr/bin/env python3
"""
NOTICE 

-- LLM GENERATED CODE BELOW --

Mirror a half-channel SEM inflow profile (y in [0,h], one no-slip wall +
symmetry/free-slip plane at the centreline) into a full-channel profile
(y in [0,2h], no-slip at both walls).

Input/output format matches src/sem.f90's inflow_profile_file reader:
free-form text, '#'-comment/blank lines skipped, data rows

    y  U  V  W  uu  vv  ww  uv  uw  vw

The solver only reads columns 1,2,5,6,7,8 (y,U,uu,vv,ww,uv); the other
columns (V,W,uw,vw) are carried through/mirrored too for a complete,
physically consistent file.

Under the reflection y -> 2h-y about the centreline, V (and any correlation
linear in v: uv, vw) flips sign since v itself changes sign under the
reflection while u,w don't; U,W,uu,vv,ww,uw stay the same (EVEN).

Usage:
    python3 mirror_half_channel.py half_channel.csv full_channel.csv
    python3 mirror_half_channel.py half_channel.csv full_channel.csv --h 1.0
"""

import argparse
import sys


def read_profile(path):
    rows = []
    header_lines = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith('#'):
                header_lines.append(line.rstrip('\n'))
                continue
            parts = s.split()
            if len(parts) < 8:
                sys.exit(f"ERROR: line has fewer than 8 columns: {line!r}")
            # pad missing uw,vw with 0 if only 8 columns given
            while len(parts) < 10:
                parts.append('0.0')
            rows.append([float(p) for p in parts[:10]])
    if len(rows) < 2:
        sys.exit("ERROR: fewer than 2 data rows found")
    rows.sort(key=lambda r: r[0])
    return header_lines, rows


def mirror(rows, h):
    """Build the full-channel row list: originals (y<=h) plus the mirror
    image (y -> 2h-y) of every row with y < h (the y=h row, if present,
    maps to itself and must not be duplicated)."""
    out = []
    for y, U, V, W, uu, vv, ww, uv, uw, vw in rows:
        if y > h + 1e-10:
            sys.exit(f"ERROR: input row at y={y} exceeds the given half-height h={h} "
                      "-- pass --h explicitly if the profile's own max(y) isn't h")
        out.append([y, U, V, W, uu, vv, ww, uv, uw, vw])
        if y < h - 1e-10:
            neg = lambda x: -x if x != 0.0 else 0.0   # avoid printing "-0.000000e+00"
            out.append([2*h - y, U, neg(V), W, uu, vv, ww, neg(uv), uw, neg(vw)])
    out.sort(key=lambda r: r[0])
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help="half-channel profile (y in [0,h])")
    ap.add_argument("output", help="full-channel profile to write (y in [0,2h])")
    ap.add_argument("--h", type=float, default=None,
                     help="channel half-height (default: max(y) in the input file)")
    args = ap.parse_args()

    header_lines, rows = read_profile(args.input)
    h = args.h if args.h is not None else rows[-1][0]

    full_rows = mirror(rows, h)

    with open(args.output, "w") as f:
        f.write(f"# Full-channel profile mirrored from {args.input} about y=h={h:.6g}\n")
        f.write("# U,W,uu,vv,ww,uw mirrored even; V,uv,vw mirrored odd (negated)\n")
        for hl in header_lines:
            f.write(hl + "\n")
        f.write("# y            U             V             W             "
                "uu            vv            ww            uv            uw            vw\n")
        for y, U, V, W, uu, vv, ww, uv, uw, vw in full_rows:
            f.write(f"{y:.6e}  {U:.6e}  {V:.6e}  {W:.6e}  "
                     f"{uu:.6e}  {vv:.6e}  {ww:.6e}  {uv:.6e}  {uw:.6e}  {vw:.6e}\n")

    print(f"wrote {len(full_rows)} rows spanning y=[0, {2*h:.6g}] to {args.output}")
    print(f"  (from {len(rows)} input rows spanning y=[0, {h:.6g}])")


if __name__ == "__main__":
    main()
