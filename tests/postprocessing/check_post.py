#!/usr/bin/env python3
"""dopamine_post reads what the solver writes: a tiny AMD (default sgs_model) Boussinesq case with box snapshots and time averages is run, then
every reader is checked (flag autodetect, wrong flags raise, box snapshot == full snapshot, tavg, prefix detection) and importing the package must
leave the matplotlib backend alone."""
import argparse, os, re, shutil, subprocess, sys, tempfile
import numpy as np
try:
    import matplotlib  # noqa: F401
except ImportError:
    print('SKIP: matplotlib not installed')
    sys.exit(77)

p = argparse.ArgumentParser()
p.add_argument('--mpirun', default='mpirun')
p.add_argument('--exe', required=True)
p.add_argument('--case-dir', required=True)
p.add_argument('--post-dir', required=True)
a = p.parse_args()
sys.path.insert(0, a.post_dir)

# the package must not switch the backend: set a file backend that is not Agg, import, compare
import matplotlib
matplotlib.use('pdf')
import dopamine_post  # noqa: E402
from dopamine_post.fields import FieldSnapshot, FieldSeries, _detect_prefix  # noqa: E402
from dopamine_post.tavg import read_tavg  # noqa: E402
assert matplotlib.get_backend().lower() == 'pdf', f'import dopamine_post changed the matplotlib backend to {matplotlib.get_backend()}'

tmp = tempfile.mkdtemp(prefix='post_check_')
try:
    for sub in ('fields', 'restart', 'stats'):
        os.makedirs(os.path.join(tmp, sub))
    txt = open(os.path.join(a.case_dir, 'input_parameters')).read()
    txt = re.sub(r'sgs_model\s*=\s*\d+,', '', txt)                  # missing key: the solver default (AMD) applies
    txt = re.sub(r'wm_filter_mode\s*=\s*-?\d+,', '', txt)
    open(os.path.join(tmp, 'input_parameters'), 'w').write(txt)
    shutil.copy(os.path.join(a.case_dir, 'input_parameters'), os.path.join(tmp, 'input_parameters.orig'))
    r = subprocess.run([a.mpirun, '-np', '1', a.exe], cwd=tmp, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if r.returncode != 0:
        print(r.stdout); sys.exit('FAIL: solver run failed')
    fdir = os.path.join(tmp, 'fields')
    names = sorted(os.listdir(fdir))
    prefix = _detect_prefix(fdir)
    assert prefix == 'chunk_small', f'prefix {prefix!r} from {names}'
    series = FieldSeries(tmp)
    step = series.steps()[-1]
    snap = series.read(step)
    assert {'U', 'V', 'W', 'P', 'nu_t', 'T'} <= set(snap.fields), f'default flags: {sorted(snap.fields)} (nu_t must be present with the AMD default)'
    assert np.all(snap['nu_t'] >= 0), 'nu_t negative'
    full = os.path.join(fdir, f'{prefix}.{step}')
    try:
        FieldSnapshot.read(full, sgs_model=0, sediment_flag=0, boussinesq_flag=1, vof_active=0)   # lenient explicit flags: T read from the nu_t block position
    except ValueError:
        pass
    try:
        FieldSnapshot.read(full, sgs_model=2, sediment_flag=1, boussinesq_flag=1, vof_active=0)   # more blocks than the file has
        sys.exit('FAIL: a flag set that needs more blocks than the file has must raise')
    except ValueError:
        pass
    box = FieldSnapshot.from_box(os.path.join(fdir, f'{prefix}Block.{step}'))
    f1, f2 = box.box['f1'], box.box['f2']
    sl = tuple(slice(f1[i] - 1, f2[i] - 1) for i in range(3))     # box cells inside the full interior grid (1-based faces -> 0-based cells)
    for name in ('U', 'V', 'W', 'P', 'T', 'nu_t'):
        ref = snap[name][sl]
        d = np.abs(box[name] - ref).max() / max(np.abs(ref).max(), 1e-30)
        assert d < 1e-12, f'box vs full snapshot {name}: {d:.2e}'
    t_avg, tf = read_tavg(os.path.join(fdir, f'{prefix}_tavg.{step}'))
    assert t_avg > 0 and "u'u'" in tf and 'T' in tf, sorted(tf)
    from dopamine_post.fields import plot_profile
    fig = plot_profile(snap)
    assert fig is not None and matplotlib.get_backend().lower() == 'pdf', 'plot_profile must return a live figure and keep the backend'
    print('PASS: dopamine_post readers (flags autodetect incl. AMD default, wrong flags raise, box == full, tavg, backend untouched)')
finally:
    shutil.rmtree(tmp, ignore_errors=True)
