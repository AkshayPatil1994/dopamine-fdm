#!/usr/bin/env python3
"""Each unsupported feature combination of the VOF solver must abort at input with its message (and the unmodified case must not)."""
import argparse, os, re, shutil, subprocess, sys, tempfile

CASES = [
    ('vof_flow without vof_active', r'vof_active\s*=\s*1', 'vof_active = 0', 'requires vof_active=1'),
    ('UAV with VOF', r'&VOF', '&UAV\n  uav_active = 1,\n/\n&VOF', 'cannot be combined with the UAV'),
    ('wave inlet without wave', r'x_bc_type\s*=\s*0', 'x_bc_type = 1', None),
    ('wave_type without inflow', r'&VOF', '&WAVES\n  wave_type = 1, wave_height = 0.01, wave_period = 1.0,\n/\n&VOF', 'requires x_bc_type=1 and inflow_type=3'),
    ('current without wave', r'&VOF', '&WAVES\n  current_type = 1, current_U = 0.1,\n/\n&VOF', 'current_type>0 requires wave_type>0'),
    ('turbulence without donor', r'&VOF', '&WAVES\n  wave_turb = 1,\n/\n&VOF', 'wave_turb=1 needs'),
    ('staircase IBM with vof_flow', r'&VOF', '&IBM\n  ibm_method = 1,\n/\n&VOF', 'always uses the staircase wall stress'),
    ('Boussinesq grav mismatch', r'&VOF', '&BOUSSINESQ\n  boussinesq_flag = 1, grav = 3.0,\n/\n&VOF', 'needs grav'),
    ('Boussinesq with vof_hsplit', r'&VOF', '&BOUSSINESQ\n  boussinesq_flag = 1, grav = 0.0,\n/\n&VOF\n  vof_hsplit = 1,', 'vof_hsplit = 0'),
    ('particles rho_f mismatch', r'&VOF', '&PARTICLES\n  particles_active = 1, particle_rho_f = 500.0,\n/\n&VOF', 'particle_rho_f must equal vof_rho_l'),
    ('particles seeded in gas', r'&VOF', '&PARTICLES\n  particles_active = 1, n_particles_init = 50, particle_rho_f = 1000.0,\n/\n&VOF', 'seeded in the gas|not ported to the GPU'),
    ('inflow_type out of range', r'&VOF', '&INFLOW\n  inflow_type = 7,\n/\n&VOF', 'inflow_type must be 0..3'),
    ('y periodic', r'y_bc_type\s*=\s*1', 'y_bc_type = 0', 'wall-bounded y only'),
    ('IBM surface dump', r'&VOF', '&IBM\n  ibm_input_mode = 1, ibm_sdf_file = \'none\', ibm_surface_nsampling = 10,\n/\n&VOF', 'IBM surface field dump'),
]

# accepted inputs that carry ignored settings must say so (the run goes on for one step)
WARNS = [
    ('waves without wave_type', r'&VOF', '&WAVES\n  wave_height = 0.02,\n/\n&VOF', 'wave_type = 0: no waves'),
    ('particles disabled', r'&VOF', '&PARTICLES\n  n_particles_init = 10,\n/\n&VOF', 'particles_active = 0'),
    ('nu unused with vof_flow', r'&VOF', '&VOF', 'nu +not used'),
]

ap = argparse.ArgumentParser()
ap.add_argument('--exe', required=True)
ap.add_argument('--mpirun', default='mpirun')
a = ap.parse_args()
src = os.path.dirname(os.path.abspath(__file__))
tmp = tempfile.mkdtemp(prefix='vofabort_')
bad = 0
try:
    for name, pat, rep, msg in CASES:
        if msg is None:
            continue
        wd = os.path.join(tmp, re.sub(r'\W', '_', name))
        shutil.copytree(src, wd)
        os.makedirs(os.path.join(wd, 'restart'), exist_ok=True)
        f = os.path.join(wd, 'input_parameters')
        text = open(f).read()
        new = re.sub(pat, lambda m: rep, text, count=1)
        assert new != text, name
        open(f, 'w').write(new)
        r = subprocess.run([a.mpirun, '-np', '1', a.exe], cwd=wd, capture_output=True, text=True, timeout=120)
        out = r.stdout + r.stderr
        ok = r.returncode != 0 and re.search(msg, out) is not None
        print('%-30s %s' % (name, 'ok' if ok else 'FAIL'))
        bad += not ok
    for name, pat, rep, msg in WARNS:
        wd = os.path.join(tmp, re.sub(r'\W', '_', name))
        shutil.copytree(src, wd)
        os.makedirs(os.path.join(wd, 'restart'), exist_ok=True)
        f = os.path.join(wd, 'input_parameters')
        text = re.sub(r'nsteps\s*=\s*-?\d+', 'nsteps = 1', open(f).read())
        open(f, 'w').write(re.sub(pat, lambda m: rep, text, count=1))
        r = subprocess.run([a.mpirun, '-np', '1', a.exe], cwd=wd, capture_output=True, text=True, timeout=120)
        ok = r.returncode == 0 and re.search(msg, r.stdout + r.stderr) is not None
        print('%-30s %s' % (name, 'ok' if ok else 'FAIL'))
        bad += not ok
finally:
    shutil.rmtree(tmp, ignore_errors=True)
sys.exit(1 if bad else 0)
