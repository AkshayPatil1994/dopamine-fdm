#!/usr/bin/env python3
"""Run the tgv_small regression case at two rank/device counts and compare
output_monitor's per-step diagnostics (src/monitor.f90:61) within tolerance.

Used both as the CPU-vs-CPU MPI rank-count parity check (Phase 1) and,
unchanged, as the CPU-vs-GPU and single-GPU-vs-multi-GPU parity check later
(Phase 3+) -- just point --exe-a/--exe-b/--np-a/--np-b/--env-a/--env-b at the
GPU build and a CUDA_VISIBLE_DEVICES-restricted rank.

Exit code 0 = all monitored columns agree within tolerance at every step
printed by both runs; 1 = a mismatch or a run failure.
"""
import argparse
import array
import glob
import os
import re
import shutil
import struct
import subprocess
import sys
import tempfile

# src/monitor.f90:63 -- istep, t, meanU, maxU, max_divergence, cfl_conv, cfl_visc, cfl_accel, dt, wall_dt_s
MONITOR_LINE = re.compile(
    r'^\s*(\d+)\s+([-\d.E+]+)\s+([-\d.E+]+)\s+([-\d.E+]+)\s+([-\d.E+]+)\s+'
    r'([-\d.E+]+)\s+([-\d.E+]+)\s+([-\d.E+]+)\s+([-\d.E+]+)\s+[-\d.]+\s*$'
)
COLUMNS = ['t', 'meanU', 'maxU', 'max_divergence', 'cfl_conv', 'cfl_visc', 'cfl_accel', 'dt']


def parse_monitor(log_text):
    rows = {}
    for line in log_text.splitlines():
        m = MONITOR_LINE.match(line)
        if m:
            istep = int(m.group(1))
            rows[istep] = [float(g) for g in m.groups()[1:]]
    return rows


def clear_snapshots(case_dir):
    for f in glob.glob(os.path.join(case_dir, 'fields', '*.[0-9]*')):
        os.remove(f)


def read_snapshot(case_dir, interior_only=False):
    """Parse the newest big-endian stream snapshot written by output_data (src/input_output.f90): six mesh
    blocks (int32 count + float64 array: x,y,z,xm,ym,zm), then one (int32 nx,ny,nz + float64 nx*ny*nz) block
    per field (U,V,W,P, then C/nu_t/T when active). Returns a list of float64 arrays, one per field block.

    interior_only drops the ghost layers of cell-centred (ghosted) dimensions. Those layers hold boundary/seam
    ghost values (e.g. the wall-ghost row at a rank seam) that the interior update never reads and that depend
    on the order of the BC and halo passes, so they are not a meaningful parity quantity.
    """
    files = sorted(glob.glob(os.path.join(case_dir, 'fields', '*.[0-9]*')),
                   key=lambda f: int(f.rsplit('.', 1)[1]))
    if not files:
        print('ERROR: --fields given but no snapshot was written (set nsave in the case input)', file=sys.stderr)
        sys.exit(1)
    b = open(files[-1], 'rb').read()
    p = 0
    mesh = []
    for _ in range(6):
        n, = struct.unpack_from('>i', b, p)
        mesh.append(n)
        p += 4 + 8 * n
    ncell = mesh[3:6]  # xm, ym, zm lengths = interior cell counts; a ghosted dimension has ncell + 2 points
    blocks = []
    while p < len(b):
        nx, ny, nz = struct.unpack_from('>iii', b, p)
        p += 12
        n = nx * ny * nz
        a = array.array('d')
        a.frombytes(b[p:p + 8 * n])
        a.byteswap()
        p += 8 * n
        if interior_only:
            dims = (nx, ny, nz)
            lo = [1 if dims[d] == ncell[d] + 2 else 0 for d in range(3)]
            hi = [dims[d] - lo[d] for d in range(3)]
            a = array.array('d', (a[i + nx * (j + ny * k)]
                                  for k in range(lo[2], hi[2])
                                  for j in range(lo[1], hi[1])
                                  for i in range(lo[0], hi[0])))
        blocks.append(a)
    return blocks


def clear_files(case_dir, patterns):
    for pat in patterns:
        for f in glob.glob(os.path.join(case_dir, pat)):
            os.remove(f)


def read_raw_files(case_dir, patterns):
    """Little-endian float64 raw files (e.g. RSB/slice output, no record markers), keyed by relative path."""
    out = {}
    for pat in patterns:
        for f in sorted(glob.glob(os.path.join(case_dir, pat))):
            a = array.array('d')
            b = open(f, 'rb').read()
            a.frombytes(b[:len(b) // 8 * 8])
            out[os.path.relpath(f, case_dir)] = a
    return out


def run_case(mpirun, exe, np, case_dir, extra_env):
    env = os.environ.copy()
    env.update(extra_env)
    cmd = [mpirun, '-np', str(np), exe]
    result = subprocess.run(cmd, cwd=case_dir, env=env,
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                             text=True, timeout=300)
    if result.returncode != 0:
        print(result.stdout)
        print(f'ERROR: {cmd} exited {result.returncode}', file=sys.stderr)
        sys.exit(1)
    return result.stdout


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case-dir', required=True)
    p.add_argument('--mpirun', default='mpirun')
    p.add_argument('--mpirun-b', default=None, help='mpirun for run B (default: --mpirun); needed when A and B use different MPI stacks, e.g. CPU vs NVHPC GPU build')
    p.add_argument('--exe-a', required=True)
    p.add_argument('--exe-b', required=True)
    p.add_argument('--np-a', type=int, default=1)
    p.add_argument('--np-b', type=int, default=2)
    p.add_argument('--env-a', action='append', default=[], help='KEY=VALUE, repeatable')
    p.add_argument('--env-b', action='append', default=[], help='KEY=VALUE, repeatable')
    p.add_argument('--rtol', type=float, default=1e-2,
                    help='relative tolerance (default 1e-2: covers MPI-reduction-order '
                         'and GPU-vs-CPU FP-associativity differences on a chaotic-ish '
                         'nonlinear flow, not just roundoff)')
    p.add_argument('--atol', type=float, default=1e-6,
                    help='absolute tolerance floor (default 1e-6), combined with --rtol '
                         'as |a-b| <= atol + rtol*max(|a|,|b|) -- needed because meanU '
                         'oscillates through zero on this decaying-TGV case, where a '
                         'pure relative tolerance is meaningless')
    p.add_argument('--fields', action='store_true',
                    help='also compare the final field snapshot block by block, ghost layers excluded (needs nsave in the case input); '
                         'catches errors in P, nu_t, scalars, ghost cells that the monitor columns cannot see')
    p.add_argument('--field-tol', type=float, default=1e-7,
                    help='per-block bound on max|a-b| relative to max|a| (default 1e-7)')
    p.add_argument('--files', action='append', default=[],
                    help='glob (relative to the case dir) of raw little-endian float64 output files to also compare, '
                         'e.g. stats/rsb_*.bin; repeatable')
    p.add_argument('--p-grid', default=None, metavar='ROW,COL',
                    help='force the 2decomp p_row,p_col pencil grid (e.g. 2,2) in a temporary copy of the case, for '
                         'BOTH runs (a 2D rank layout needs np-a/np-b == ROW*COL wherever it is used; the reference '
                         'run keeps the case default when its np-a differs -- see --p-grid-b-only)')
    p.add_argument('--p-grid-b-only', action='store_true',
                    help='apply --p-grid to run B only; run A keeps the case default (auto) decomposition')
    p.add_argument('--label-a', default='A')
    p.add_argument('--label-b', default='B')
    args = p.parse_args()

    def parse_env(pairs):
        out = {}
        for kv in pairs:
            k, _, v = kv.partition('=')
            out[k] = v
        return out

    case_a = case_b = args.case_dir
    tmp_dirs = []
    if args.p_grid:
        row, col = args.p_grid.split(',')

        def with_grid(src):
            d = tempfile.mkdtemp(prefix='pgrid_')
            tmp_dirs.append(d)
            shutil.copytree(src, os.path.join(d, 'case'))
            f = os.path.join(d, 'case', 'input_parameters')
            txt = open(f).read()
            txt = re.sub(r'p_row\s*=\s*\d+', f'p_row = {row}', txt)
            txt = re.sub(r'p_col\s*=\s*\d+', f'p_col = {col}', txt)
            open(f, 'w').write(txt)
            return os.path.join(d, 'case')

        case_b = with_grid(args.case_dir)
        case_a = args.case_dir if args.p_grid_b_only else with_grid(args.case_dir)
    if args.fields:
        clear_snapshots(case_a)
    clear_files(case_a, args.files)
    out_a = run_case(args.mpirun, args.exe_a, args.np_a, case_a, parse_env(args.env_a))
    snap_a = read_snapshot(case_a, interior_only=True) if args.fields else None
    raw_a = read_raw_files(case_a, args.files)
    if args.fields:
        clear_snapshots(case_b)
    clear_files(case_b, args.files)
    out_b = run_case(args.mpirun_b or args.mpirun, args.exe_b, args.np_b, case_b, parse_env(args.env_b))
    snap_b = read_snapshot(case_b, interior_only=True) if args.fields else None
    raw_b = read_raw_files(case_b, args.files)
    clear_files(case_b, args.files)
    for d in tmp_dirs:
        shutil.rmtree(d, ignore_errors=True)

    rows_a = parse_monitor(out_a)
    rows_b = parse_monitor(out_b)

    if not rows_a or not rows_b:
        print('ERROR: no monitor lines parsed from one or both runs', file=sys.stderr)
        sys.exit(1)

    common_steps = sorted(set(rows_a) & set(rows_b))
    if not common_steps:
        print('ERROR: no common monitored steps between the two runs', file=sys.stderr)
        sys.exit(1)

    worst = {}
    failed = False
    for istep in common_steps:
        for col, va, vb in zip(COLUMNS, rows_a[istep], rows_b[istep]):
            diff = abs(va - vb)
            bound = args.atol + args.rtol * max(abs(va), abs(vb))
            rel = diff / max(abs(va), abs(vb), 1e-12)  # reported for visibility only
            if rel > worst.get(col, 0.0):
                worst[col] = rel
            if diff > bound:
                failed = True
                print(f'MISMATCH step {istep} col {col}: '
                      f'{args.label_a}={va:.6E} {args.label_b}={vb:.6E} '
                      f'|diff|={diff:.3E} > bound={bound:.3E} (atol={args.atol:.1E}, rtol={args.rtol:.1E})')

    print(f'Compared {len(common_steps)} steps, {args.label_a} (np={args.np_a}) vs '
          f'{args.label_b} (np={args.np_b}), rtol={args.rtol:.1E}')
    for col in COLUMNS:
        print(f'  worst-case relative diff [{col}]: {worst.get(col, 0.0):.3E}')

    if args.fields:
        if len(snap_a) != len(snap_b):
            failed = True
            print(f'MISMATCH: snapshot has {len(snap_a)} field blocks vs {len(snap_b)}')
        else:
            for ib, (fa, fb) in enumerate(zip(snap_a, snap_b)):
                scale = max(max(abs(v) for v in fa), 1e-30)
                dmax = max(abs(x - y) for x, y in zip(fa, fb))
                print(f'  field block {ib}: max|a-b|/max|a| = {dmax / scale:.3E}')
                if dmax / scale > args.field_tol:
                    failed = True
                    print(f'MISMATCH field block {ib}: {dmax / scale:.3E} > {args.field_tol:.1E}')

    if args.files:
        if not raw_a or set(raw_a) != set(raw_b):
            failed = True
            print(f'MISMATCH: raw output files differ: {sorted(raw_a)} vs {sorted(raw_b)}')
        for name in sorted(set(raw_a) & set(raw_b)):
            fa, fb = raw_a[name], raw_b[name]
            if len(fa) != len(fb):
                failed = True
                print(f'MISMATCH {name}: {len(fa)} vs {len(fb)} values')
                continue
            scale = max(max((abs(v) for v in fa), default=0.0), 1e-30)
            dmax = max((abs(x - y) for x, y in zip(fa, fb)), default=0.0)
            print(f'  file {name}: max|a-b|/max|a| = {dmax / scale:.3E}')
            if dmax / scale > args.field_tol:
                failed = True
                print(f'MISMATCH file {name}: {dmax / scale:.3E} > {args.field_tol:.1E}')

    if failed:
        print('RESULT: FAIL')
        sys.exit(1)
    print('RESULT: PASS')


if __name__ == '__main__':
    main()
