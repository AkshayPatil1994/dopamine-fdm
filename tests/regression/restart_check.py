#!/usr/bin/env python3
"""Restart parity: run a case for 2N steps straight, then as N steps + a hot-start restart of N more from the
snapshot written at step N, and require the two final snapshots to agree field-by-field (all blocks: U,V,W,P,
scalars, nu_t). Catches state that is not saved/restored (or not re-uploaded to the device) on restart.

The case input must have nsteps = nsave = N (as the tgv/chan/... parity cases do) and a deterministic IC.
"""
import argparse
import array
import os
import re
import shutil
import struct
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_monitor import read_snapshot, clear_snapshots  # noqa: E402
from compare_particles import read_particles  # noqa: E402


def run(mpirun, exe, np, cwd, env_extra):
    env = os.environ.copy()
    env.update(env_extra)
    r = subprocess.run([mpirun, '-np', str(np), exe], cwd=cwd, env=env, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, text=True, timeout=600)
    if r.returncode != 0:
        print(r.stdout)
        sys.exit(f'ERROR: run exited {r.returncode}')


def read_tavg_file(path):
    """(averaging time, [per-field list of values]) of a fields/<fileout>_tavg.<step> file (src/time_stats.f90)."""
    raw = open(path, 'rb').read()
    nx, ny, nz, nf = struct.unpack('<4i', raw[:16])
    t_avg = struct.unpack('<d', raw[16:24])[0]
    vals = array.array('d')
    off = 24 + 8 * nf
    vals.frombytes(raw[off:off + 8 * nx * ny * nz * nf])
    n = nx * ny * nz
    return t_avg, [vals[i * n:(i + 1) * n] for i in range(nf)]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case-dir', required=True)
    p.add_argument('--mpirun', default='mpirun')
    p.add_argument('--exe', required=True)
    p.add_argument('--np', type=int, default=1)
    p.add_argument('--tol', type=float, default=1e-12, help='per-block max|a-b|/max|a| bound (default 1e-12)')
    p.add_argument('--tavg-tol', type=float, default=1e-10,
                   help='bound on the time-averaged statistics (central moments lose digits through the restart)')
    p.add_argument('--particles', action='store_true', help='also compare the particles_restart state by particle ID')
    p.add_argument('--env', action='append', default=[], help='KEY=VALUE, repeatable')
    args = p.parse_args()
    env = dict(kv.partition('=')[::2] for kv in args.env)

    src = open(os.path.join(args.case_dir, 'input_parameters')).read()
    n = int(re.search(r'nsteps\s*=\s*(\d+)', src).group(1))
    fileout = re.search(r"fileout\s*=\s*'([^']+)'", src).group(1)

    def make_input(nsteps, restart):
        s = re.sub(r'nsteps\s*=\s*\d+', f'nsteps = {nsteps}', src)
        s = re.sub(r'nsave\s*=\s*\d+', f'nsave = {nsteps}', s)
        if restart:
            s = re.sub(r'restart\s*=\s*0', 'restart = 1', s)
            s = re.sub(r'nstep_init\s*=\s*\d+', f'nstep_init = {n}', s)
            s = re.sub(r"filein\s*=\s*'[^']*'", f"filein = 'restart/{fileout}.{n}'", s)
        return s

    tmp = tempfile.mkdtemp(prefix='restart_check_')
    try:
        for d in ('fields', 'restart', 'stats'):
            os.makedirs(os.path.join(tmp, d))
        write = lambda s: open(os.path.join(tmp, 'input_parameters'), 'w').write(s)
        # auxiliary input files referenced by the case (UAV paths, profiles, ...)
        for f in os.listdir(args.case_dir):
            fp = os.path.join(args.case_dir, f)
            if os.path.isfile(fp) and f != 'input_parameters':
                shutil.copy(fp, tmp)

        write(make_input(2 * n, False))
        run(args.mpirun, args.exe, args.np, tmp, env)
        straight = read_snapshot(tmp, interior_only=True)
        straight_p = read_particles(os.path.join(tmp, 'particles_restart')) if args.particles else None
        tavg_path = lambda step: os.path.join(tmp, 'fields', f'{fileout}_tavg.{step}')
        straight_tavg = read_tavg_file(tavg_path(2 * n)) if os.path.exists(tavg_path(2 * n)) else None

        clear_snapshots(tmp)
        write(make_input(n, False))
        run(args.mpirun, args.exe, args.np, tmp, env)
        # clear_snapshots below would delete the average the restart must reload
        tavg_n = open(tavg_path(n), 'rb').read() if straight_tavg else None
        shutil.copy(os.path.join(tmp, 'fields', f'{fileout}.{n}'), os.path.join(tmp, 'restart', f'{fileout}.{n}'))
        wmf = os.path.join(tmp, 'fields', f'{fileout}.{n}.wmf')   # wall-model filter state next to the snapshot
        if os.path.exists(wmf):
            shutil.copy(wmf, os.path.join(tmp, 'restart', f'{fileout}.{n}.wmf'))

        clear_snapshots(tmp)
        if tavg_n:
            open(tavg_path(n), 'wb').write(tavg_n)
        write(make_input(n, True))
        run(args.mpirun, args.exe, args.np, tmp, env)
        restarted = read_snapshot(tmp, interior_only=True)
        restarted_tavg = read_tavg_file(tavg_path(2 * n)) if straight_tavg else None
        restarted_p = read_particles(os.path.join(tmp, 'particles_restart')) if args.particles else None
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    if len(straight) != len(restarted):
        sys.exit(f'FAIL: {len(straight)} vs {len(restarted)} field blocks')
    failed = False
    for i, (a, b) in enumerate(zip(straight, restarted)):
        scale = max(max(abs(v) for v in a), 1e-30)
        d = max(abs(x - y) for x, y in zip(a, b)) / scale
        print(f'  field block {i}: restart vs straight max|a-b|/max|a| = {d:.3E}')
        failed |= d > args.tol
    if straight_tavg:
        dt_avg = abs(straight_tavg[0] - restarted_tavg[0]) / straight_tavg[0]
        print(f'  tavg averaging time: restart vs straight relative difference = {dt_avg:.3E}')
        failed |= dt_avg > 1e-12
        for i, (a, b) in enumerate(zip(straight_tavg[1], restarted_tavg[1])):
            scale = max(max(abs(v) for v in a), 1e-30)
            d = max(abs(x - y) for x, y in zip(a, b)) / scale
            print(f'  tavg field {i}: restart vs straight max|a-b|/max|a| = {d:.3E}')
            failed |= d > args.tavg_tol
    if args.particles:
        if set(straight_p) != set(restarted_p):
            sys.exit('FAIL: the particle ID sets differ')
        dp = max(abs(x - y) for i in straight_p for x, y in zip(straight_p[i], restarted_p[i]))
        print(f'  particles ({len(straight_p)}): restart vs straight max|a-b| = {dp:.3E}')
        failed |= dp > args.tol
    print('RESULT: FAIL' if failed else 'RESULT: PASS')
    sys.exit(1 if failed else 0)


if __name__ == '__main__':
    main()
