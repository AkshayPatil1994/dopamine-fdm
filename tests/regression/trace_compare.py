#!/usr/bin/env python3
"""Compare two DOPAMINE_TRACE_DIR dumps (src/debug_trace.f90) and report where the decomposition first changes the result.

  trace_compare.py <dir_a> <dir_b> [--tol 1e-10] [--all]

Files are <seq>_<tag>.<field>: 3 int32 (n1,n2,n3) + float64 global array (ghost cells included). Points are visited
in sequence order; the first field whose max |a-b| exceeds tol*(max(a)-min(a)) is reported with its global (i,j,k), split
into 'interior' (away from the outer edge of the global array) and 'edge' (physical-boundary ghost layer).
Also prints the ghost/owner mismatches recorded in seam.log (stale halos seen by the run itself).
Exit 0 = no interior difference above tol, 1 = difference found or files missing.
'interior' = cells away from the outer edge of the global array in x and z (ghost layers, incl. seam ghost planes'
y-ghost rows and physical-boundary ghosts, are reported as 'ghost-layer-only' and do not fail the comparison).
"""
import argparse
import array
import os
import struct
import sys


def read_block(path):
    with open(path, 'rb') as f:
        hdr = f.read(12)
        # the build may convert unformatted I/O to big-endian; pick the byte order that gives a sane header
        for order in ('<', '>'):
            n1, n2, n3 = struct.unpack(order + '3i', hdr)
            if 0 < n1 < 100000 and 0 < n2 < 100000 and 0 < n3 < 100000:
                break
        a = array.array('d')
        a.fromfile(f, n1 * n2 * n3)
        if (order == '>') == (sys.byteorder == 'little'):
            a.byteswap()
    return (n1, n2, n3), a


def compare_field(pa, pb):
    (n1, n2, n3), a = read_block(pa)
    shape_b, b = read_block(pb)
    if shape_b != (n1, n2, n3):
        return None, ('shape', (n1, n2, n3), shape_b)
    # scale by the field's spread, not its magnitude: a fluctuation field on a large offset (T ~ T_ref) would otherwise hide real differences
    # ... with a floor of 1e-8*max|a| so a nearly uniform field (spread ~ roundoff) does not turn roundoff into a 'difference'
    scale = max(max(a) - min(a), 1e-8 * max(abs(max(a)), abs(min(a))), 1e-30)
    worst, loc = 0.0, None
    worst_int, loc_int = 0.0, None
    for idx in range(len(a)):
        d = abs(a[idx] - b[idx])
        if d > worst:
            worst = d
            loc = idx
        if d > worst_int:
            i = idx % n1
            k = idx // (n1 * n2)
            if 0 < i < n1 - 1 and 0 < k < n3 - 1:
                worst_int, loc_int = d, idx
    def ijk(idx):
        return (idx % n1 + 1, (idx // n1) % n2 + 1, idx // (n1 * n2) + 1) if idx is not None else None
    return (worst / scale, ijk(loc), worst_int / scale, ijk(loc_int)), None


def read_seam(d):
    out = []
    p = os.path.join(d, 'seam.log')
    if os.path.exists(p):
        for line in open(p):
            if 'n_bad=0' not in line:
                out.append(line.rstrip())
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('dir_a')
    ap.add_argument('dir_b')
    ap.add_argument('--tol', type=float, default=1e-10)
    ap.add_argument('--skip', action='append', default=[], help='field name to ignore (e.g. P), repeatable')
    ap.add_argument('--all', action='store_true', help='list every differing point, not only the first')
    args = ap.parse_args()

    names_a = sorted(f for f in os.listdir(args.dir_a) if f[:4].isdigit())
    names_b = set(f for f in os.listdir(args.dir_b) if f[:4].isdigit())
    if not names_a or set(names_a) != names_b:
        print('ERROR: dumps missing or different point/field sets', file=sys.stderr)
        return 1

    found = 0
    ghost_only = []
    for n in names_a:
        if n.rsplit('.', 1)[1] in args.skip:
            continue
        res, err = compare_field(os.path.join(args.dir_a, n), os.path.join(args.dir_b, n))
        if err:
            print(f'{n}: {err}')
            found += 1
            continue
        rel, loc, rel_int, loc_int = res
        # gate on the interior: ghost layers (seam ghost-plane y-rows, physical-edge ghosts) are reported separately
        if rel_int > args.tol:
            found += 1
            print(f'{n:32s} INTERIOR rel diff {rel_int:9.2e} at (i,j,k)={loc_int}')
            if not args.all:
                break
        elif rel > args.tol:
            ghost_only.append((n, rel, loc))
    for tag, d in (('A', args.dir_a), ('B', args.dir_b)):
        bad = read_seam(d)
        if bad:
            print(f'--- ghost/owner mismatches recorded by run {tag} ({len(bad)} lines, first 8):')
            for line in bad[:8]:
                print('   ', line)
    if ghost_only:
        n, rel, loc = ghost_only[0]
        print(f'ghost-layer-only differences in {len(ghost_only)} dumps (first: {n} {rel:.2e} at {loc}); not gated')
    if not found:
        print(f'no interior difference above {args.tol:g} in {len(names_a)} dumps')
    return 1 if found else 0


if __name__ == '__main__':
    sys.exit(main())
