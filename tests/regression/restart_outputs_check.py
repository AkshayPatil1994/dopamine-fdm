#!/usr/bin/env python3
"""Restart of the time-series outputs: a run that is stopped after the restart snapshot and resumed from it must give the same
probe files, vof_diag.dat, vof_gauges.dat (and ibm_forces.csv) as an uninterrupted run, with no duplicated or lost records.

The case input (a VOF case with fixed dt, nsteps = N) gets two line probes, a slice probe and wave gauges added. Runs, all in one
directory per scenario:  straight  2N steps;  stopped  N+K steps (snapshot at N), then restart from the snapshot for N steps
(K = --extra, not a multiple of the probe intervals; N is a multiple of them, since a restarted run counts its probe
interval from its own first step);  plain  N steps, then restart (K = 0).
"""
import argparse, glob, os, re, shutil, subprocess, sys, tempfile
import numpy as np

STAT = """
&STATISTICS
  n_lines = 2, line_freq = 8,
  line_dir(1) = 'y', line_pos1(1) = 1.0, line_pos2(1) = 1.0, line_comps(1) = 'UVW', line_fileout(1) = 'probes/l1',
  line_dir(2) = 'x', line_pos1(2) = 1.0, line_pos2(2) = 1.0, line_comps(2) = 'UP', line_fileout(2) = 'probes/l2',
  n_slices = 1, slice_freq = 10, slice_dir(1) = 'z', slice_pos(1) = 1.5, slice_comps(1) = 'UVW', slice_fileout(1) = 'probes/s1',
/
&WAVES
  wave_gauge_x = 1.0, 3.0,
/
"""


def run(mpirun, exe, nproc, cwd):
    cmd = [mpirun, '-np', str(nproc), exe]
    if 'oversubscribe' in subprocess.run([mpirun, '--help'], capture_output=True, text=True).stdout:
        cmd.insert(1, '--oversubscribe')
    r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=900)
    open(os.path.join(cwd, 'run.log'), 'a').write(r.stdout + r.stderr)
    if r.returncode != 0:
        print(r.stdout[-3000:], r.stderr[-2000:])
        sys.exit(f'ERROR: run exited {r.returncode} in {cwd}')


def text_rows(path):
    heads, rows = 0, []
    for l in open(path):
        if l.startswith('#') or l.startswith('step'):
            heads += 1
        elif l.strip():
            rows.append([float(x) for x in l.replace(',', ' ').split()])
    return heads, np.array(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--case-dir', required=True)
    ap.add_argument('--exe', required=True)
    ap.add_argument('--mpirun', default='mpirun')
    ap.add_argument('--np', type=int, default=2)
    ap.add_argument('--extra', type=int, default=13, help='steps run after the snapshot before the "stop" (K)')
    ap.add_argument('--tol', type=float, default=1e-9)
    a = ap.parse_args()

    src = open(os.path.join(a.case_dir, 'input_parameters')).read()
    n = 40
    fileout = re.search(r"fileout\s*=\s*'([^']+)'", src).group(1)

    def make(nsteps, restart):
        s = re.sub(r'nsteps\s*=\s*\d+', f'nsteps = {nsteps}', src)
        s = re.sub(r'nsave\s*=\s*\d+', f'nsave = {n}', s)
        if restart:
            s = re.sub(r'restart\s*=\s*0', 'restart = 1', s)
            s = re.sub(r'(restart = 1,?)', rf'\1 nstep_init = {n},', s, count=1)
            s = re.sub(r"filein\s*=\s*'[^']*'", f"filein = 'restart/{fileout}.{n}'", s)
        if 'nsampling' not in s:
            s = s.replace('ibm_input_mode = 1,', 'ibm_input_mode = 1, nsampling = 5,')
        return s + STAT

    def newdir(tmp, name):
        d = os.path.join(tmp, name)
        shutil.copytree(a.case_dir, d)
        for sub in ('fields', 'restart', 'stats', 'probes'):
            os.makedirs(os.path.join(d, sub), exist_ok=True)
        return d

    def go(d, text):
        open(os.path.join(d, 'input_parameters'), 'w').write(text)
        run(a.mpirun, a.exe, a.np, d)

    tmp = tempfile.mkdtemp(prefix='restart_outputs_')
    failed = False
    try:
        ref = newdir(tmp, 'straight')
        go(ref, make(2 * n, False))
        scen = {}
        for name, extra in (('stopped', a.extra), ('plain', 0)):
            d = newdir(tmp, name)
            go(d, make(n + extra, False))
            shutil.copy(os.path.join(d, 'fields', f'{fileout}.{n}'), os.path.join(d, 'restart', f'{fileout}.{n}'))
            go(d, make(n, True))
            scen[name] = d

        names = sorted(os.path.relpath(p, ref) for p in glob.glob(os.path.join(ref, 'probes', '*')))
        names += [f for f in ('vof_diag.dat', 'vof_gauges.dat', 'ibm_forces.csv') if os.path.exists(os.path.join(ref, f))]
        print('compared files:', ', '.join(names))
        for name, d in scen.items():
            for f in names:
                pa, pb = os.path.join(ref, f), os.path.join(d, f)
                if not os.path.exists(pb):
                    print(f'  [{name}] {f}: MISSING'); failed = True; continue
                if f.endswith('.bin'):
                    A, B = np.fromfile(pa, '>f8'), np.fromfile(pb, '>f8')
                    ok = A.shape == B.shape and np.max(np.abs(A - B)) <= a.tol * max(np.max(np.abs(A)), 1e-30)
                    info = f'{A.size} vs {B.size} values'
                elif f.endswith('_meta.txt'):
                    ok = open(pa).read() == open(pb).read()
                    info = ''
                else:
                    ha, A = text_rows(pa); hb, B = text_rows(pb)
                    if f != 'vof_gauges.dat':   # the step counter restarts at 0 in a restarted run
                        A, B = A[:, 1:], B[:, 1:]
                    ok = hb == ha and A.shape == B.shape and np.allclose(A, B, rtol=a.tol, atol=a.tol * max(np.nanmax(np.abs(A)), 1e-30), equal_nan=True)
                    info = f'{A.shape[0]} vs {B.shape[0]} rows, headers {ha} vs {hb}'
                print(f'  [{name}] {f}: {"ok" if ok else "FAIL"} {info}')
                failed |= not ok
    finally:
        if not failed or os.environ.get('KEEP') is None:
            shutil.rmtree(tmp, ignore_errors=True)
        else:
            print('kept', tmp)
    print('RESULT: FAIL' if failed else 'RESULT: PASS')
    sys.exit(1 if failed else 0)


if __name__ == '__main__':
    main()
