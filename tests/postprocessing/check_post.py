#!/usr/bin/env python3
"""dopamine_post reads what the solver writes: a tiny AMD (default sgs_model) Boussinesq case with box snapshots and time averages is run, then
every reader is checked (flag autodetect, wrong flags raise, box snapshot == full snapshot, tavg, prefix detection, XDMF for single-phase/Boussinesq/box/VOF/particle cases) and importing the package must
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

    # XDMF writers: every Seek offset must lie inside the binary file and the attributes must hold what the solver wrote
    import xml.etree.ElementTree as ET
    from dopamine_post.fields import write_field_xmf, write_box_xmf
    from dopamine_post.particles import write_particles_xmf

    def xmf_attrs(xmf):
        out = {}
        for at in ET.parse(xmf).iter('Attribute'):
            for di in at.iter('DataItem'):
                if di.get('Format') == 'Binary' and di.get('Seek'):
                    n = int(np.prod([int(v) for v in di.get('Dimensions').split()]))
                    fn = os.path.join(os.path.dirname(str(xmf)), di.text.strip())
                    assert int(di.get('Seek')) + 8 * n <= os.path.getsize(fn), f'{xmf}: {at.get("Name")} runs past the end of {fn}'
                    out[at.get('Name')] = np.fromfile(fn, '>f8', n, offset=int(di.get('Seek')))
        return out

    attrs = xmf_attrs(write_field_xmf(tmp))
    assert {'Velocity', 'P', 'nu_t', 'T'} <= set(attrs), f'single-phase Boussinesq XMF attributes {sorted(attrs)}'
    assert {'Velocity', 'P'} <= set(xmf_attrs(write_box_xmf(tmp)))

    def tiny_case(name, edits):
        src = os.path.join(os.path.dirname(a.case_dir), name)
        if not os.path.exists(os.path.join(src, 'input_parameters')):
            return None
        wd = os.path.join(tmp, name)
        for sub in ('fields', 'restart', 'stats'):
            os.makedirs(os.path.join(wd, sub))
        text = open(os.path.join(src, 'input_parameters')).read()
        for pat, rep in edits:
            text = re.sub(pat, rep, text, count=1)
        open(os.path.join(wd, 'input_parameters'), 'w').write(text)
        r = subprocess.run([a.mpirun, '-np', '1', a.exe], cwd=wd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        if r.returncode != 0:
            print(r.stdout); sys.exit(f'FAIL: {name} run failed')
        return wd

    wd = tiny_case('vof_drop_small', [(r'nsteps\s*=\s*-1', 'nsteps = 4'), (r'nsave\s*=\s*\d+', 'nsave = 4'),
                                      (r'nmonitor\s*=\s*\d+', 'nmonitor = 1'), (r'&VOF', "&STATISTICS\n  n_slices = 1, slice_freq = 2, slice_dir(1) = 'z', slice_pos(1) = 0.5, "
                                       "slice_comps(1) = 'UVWA', slice_fileout(1) = 'slices/zA',\n/\n&VOF")])
    if wd:
        from dopamine_post import diagnostics as dg
        from dopamine_post.probes import SliceProbe
        vd = dg.read_vof_diag(os.path.join(wd, 'vof_diag.dat'))
        assert len(vd) == 32 and vd['step'].size == 4 and abs(vd['dVliq_rel']).max() < 1e-2, 'vof_diag.dat columns/rows'
        sp = SliceProbe.read(os.path.join(wd, 'slices', 'zA.bin'))
        assert sp.comps == ['U', 'V', 'W', 'A'] and 0 <= sp.get('A').min() and sp.get('A').max() <= 1, 'slice probe component A'
        vattrs = xmf_attrs(write_field_xmf(wd))
        assert 'vof' in vattrs, f'VOF case: the XMF has no liquid-fraction attribute ({sorted(vattrs)})'
        vs = FieldSeries(wd).latest()['vof']
        assert 0 <= vattrs['vof'].min() and vattrs['vof'].max() <= 1 and abs(vattrs['vof'].max() - vs.max()) < 1e-12, 'vof attribute does not match the snapshot'
    wd = tiny_case('part_tracer_small', [])
    if wd:
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter('error')                      # fields next to <prefix>_particles.<step> must not trigger the multiple-prefix warning
            write_field_xmf(wd)
        write_particles_xmf(wd)
        topo = ET.parse(os.path.join(wd, 'paraview', 'particles.xmf')).getroot().findall('.//Topology')
        assert len(topo) == 1 and int(topo[0].get('NumberOfElements')) == 300, 'particle XMF: one saved step with n_particles_init points expected'
        assert 'Velocity' in xmf_attrs(os.path.join(wd, 'paraview', 'particles.xmf'))
    # stretched-z grid.out/grid_z.out, legacy probe meta (requested letters incl. A, A not written), restarted-history dedupe
    import warnings
    from dopamine_post._core import coords_from_grid, comp_names
    gd = os.path.join(tmp, 'gridtest'); os.makedirs(gd)
    open(os.path.join(gd, 'geometry.out'), 'w').write('4 3 2\n 1.0 3.0 2.0\n')
    np.savetxt(os.path.join(gd, 'grid.out'), np.c_[np.arange(1, 4), [.5, 1.5, 2.5], [1, 2, 3], np.ones(3), np.ones(3)])
    assert np.allclose(coords_from_grid(os.path.join(gd, 'grid.out'))[2], [.5, 1.5])
    np.savetxt(os.path.join(gd, 'grid_z.out'), np.c_[[1, 2], [.2, 1.2], [.4, 2.0], [.4, 1.6], [2.5, .625]])
    assert np.allclose(coords_from_grid(os.path.join(gd, 'grid.out'))[2], [.2, 1.2]), 'grid_z.out ignored'
    with warnings.catch_warnings(record=True):
        warnings.simplefilter('always')
        assert comp_names('UVWA', 3) == ['U', 'V', 'W']
    try:
        comp_names('UVWPT', 3); sys.exit('FAIL: comps/ncomp mismatch must raise')
    except ValueError:
        pass
    open(os.path.join(gd, 'ibm_forces.csv'), 'w').write('step,t,Fx_ibm\n1,0.1,1\n2,0.2,2\n2,0.2,3\n')
    assert dg.read_ibm_forces(os.path.join(gd, 'ibm_forces.csv'))['Fx_ibm'].tolist() == [1, 3], 'restart rows must keep the last'
    print('PASS: dopamine_post readers (flags autodetect incl. AMD default, wrong flags raise, box == full, tavg, backend untouched)')
finally:
    shutil.rmtree(tmp, ignore_errors=True)
