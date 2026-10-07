#!/usr/bin/env python3
"""Run a VOF regression case (tests/regression/vof_*) at one or two rank counts and check the vof_diag.dat summary.

Checks on the last row of vof_diag.dat of every run: |relative liquid-volume drift| <= --max-dvol, C within [0,1] (+1e-12),
and, when given, the L1 distance of C to its initial shape (column 'L1(C-C_init)') <= --max-l1. With --np-b the second run
(optionally on --p-grid) must reproduce the first one: every column named in --same-cols agrees to --same-tol (relative).
Column numbers are the 1-based numbers of the header line of vof_diag.dat.
Further checks on the first run: --limit COL:MIN:MAX (last row), --rowmax COL:MAX (largest value over all rows),
--freq COL:THEORY:TOL (mean half period from the zero crossings of the column against pi/THEORY, relative tolerance TOL).
"""
import argparse, os, re, shutil, subprocess, sys, tempfile


def run(case_dir, exe, mpirun, np, pgrid, workdir):
    shutil.copytree(case_dir, workdir)
    os.makedirs(os.path.join(workdir, 'restart'), exist_ok=True)
    inp = os.path.join(workdir, 'input_parameters')
    text = open(inp).read()
    if pgrid:
        text = re.sub(r'p_row\s*=\s*\d+', 'p_row = %s' % pgrid.split(',')[0], text)
        text = re.sub(r'p_col\s*=\s*\d+', 'p_col = %s' % pgrid.split(',')[1], text)
        open(inp, 'w').write(text)
    cmd = [mpirun, '--oversubscribe', '-np', str(np), exe] if 'oversubscribe' in subprocess.run(
        [mpirun, '--help'], capture_output=True, text=True).stdout else [mpirun, '-np', str(np), exe]
    r = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True)
    open(os.path.join(workdir, 'run.log'), 'w').write(r.stdout + r.stderr)
    if r.returncode != 0:
        print(r.stdout[-2000:], r.stderr[-2000:])
        sys.exit('FAIL: run exited with %d (np=%d)' % (r.returncode, np))
    rows = [[float(x) for x in l.split()] for l in open(os.path.join(workdir, 'vof_diag.dat')) if not l.startswith('#')]
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--case-dir', required=True)
    ap.add_argument('--exe', required=True)
    ap.add_argument('--mpirun', default='mpirun')
    ap.add_argument('--np-a', type=int, default=1)
    ap.add_argument('--np-b', type=int, default=0)
    ap.add_argument('--p-grid', default='')
    ap.add_argument('--max-dvol', type=float, default=1e-12)
    ap.add_argument('--max-l1', type=float, default=-1)
    ap.add_argument('--same-cols', default='4,5,6,7,14,15,28')
    ap.add_argument('--same-tol', type=float, default=1e-9)
    ap.add_argument('--limit', action='append', default=[])
    ap.add_argument('--rowmax', action='append', default=[])
    ap.add_argument('--freq', default='')
    a = ap.parse_args()
    tmp = tempfile.mkdtemp(prefix='vofchk_')
    try:
        runs = [run(a.case_dir, a.exe, a.mpirun, a.np_a, '', os.path.join(tmp, 'a'))]
        if a.np_b:
            runs.append(run(a.case_dir, a.exe, a.mpirun, a.np_b, a.p_grid, os.path.join(tmp, 'b')))
        allrows = runs
        runs = [r[-1] for r in allrows]
        for k, row in enumerate(runs):
            dvol, cmin, cmax = row[4], row[5], row[6]
            l1 = row[27] if len(row) > 27 else 0.0
            print('run %d: t=%.4g dVol=%.3e C=[%.3g,%.6g] L1=%.6e' % (k, row[1], dvol, cmin, cmax, l1))
            if abs(dvol) > a.max_dvol:
                sys.exit('FAIL: volume drift %.3e > %.3e' % (abs(dvol), a.max_dvol))
            if cmin < -1e-12 or cmax > 1 + 1e-12:
                sys.exit('FAIL: C outside [0,1]')
            if a.max_l1 > 0 and l1 > a.max_l1:
                sys.exit('FAIL: L1(C-C_init) %.3e > %.3e' % (l1, a.max_l1))
        for lim in a.limit:
            c, lo, hi = lim.split(':')
            v = runs[0][int(c)-1]
            if not float(lo) <= v <= float(hi):
                sys.exit('FAIL: column %s = %.6e outside [%s, %s]' % (c, v, lo, hi))
        for lim in a.rowmax:
            c, hi = lim.split(':')
            v = max(abs(r[int(c)-1]) for r in allrows[0])
            if v > float(hi):
                sys.exit('FAIL: max |column %s| = %.6e > %s' % (c, v, hi))
        if a.freq:
            c, w, tol = a.freq.split(':')
            t = [r[1] for r in allrows[0]]
            f = [r[int(c)-1] for r in allrows[0]]
            zc = [t[i] - f[i]*(t[i+1] - t[i])/(f[i+1] - f[i]) for i in range(len(t)-1) if f[i]*f[i+1] < 0]
            if len(zc) < 2:
                sys.exit('FAIL: fewer than two zero crossings of column %s' % c)
            hp = (zc[-1] - zc[0])/(len(zc) - 1)
            err = 3.141592653589793/hp/float(w) - 1
            print('frequency error %.3f %%' % (100*err))
            if abs(err) > float(tol):
                sys.exit('FAIL: frequency error %.3f %% > %.3f %%' % (100*err, 100*float(tol)))
        if len(runs) == 2:
            for c in [int(x) for x in a.same_cols.split(',')]:
                u, v = runs[0][c-1], runs[1][c-1]
                if abs(u - v) > a.same_tol*max(abs(u), abs(v), 1e-30) and abs(u - v) > 1e-14:
                    sys.exit('FAIL: column %d differs between np=%d (%.12e) and np=%d (%.12e)' % (c, a.np_a, u, a.np_b, v))
        print('PASS')
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


main()
