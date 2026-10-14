#!/usr/bin/env python3
"""Feature-combination matrix of the VOF solver: {passive vof_flow=0, two-fluid vof_flow=1} x {free-surface tank, IBM body} x
{sediment, Boussinesq temperature, particles, all three} must run a few steps without error and without NaN, and with particles the
count stays at the seeded number. Cases that must be rejected are covered by vof_input_checks."""
import argparse, itertools, os, re, shutil, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
REG = os.path.dirname(HERE)

SED = '&SEDIMENT\n  sediment_flag = 1, d_s = 1.0e-4, rho_s = 2650.0, rho_f = %(rf)s, grav = %(g)s, Sc = 1.0, Sc_t = 1.0, C_ref = 1.0, sed_bc_bot = 0, C_ic_type = 0,\n/\n'
BOU = '&BOUSSINESQ\n  boussinesq_flag = 1, beta_T = 2.0e-4, T_ref = 20.0, grav = %(g)s, Pr = 7.0, Pr_t = 0.85, T_ic_type = 1, T_ic_grad = 1.0,\n/\n'
PAR = '&PARTICLES\n  particles_active = 1, n_particles_init = 60, particle_seed_seed = 5, particle_seed_ymin = 0.05, particle_seed_ymax = 0.4,\n  bc_particle_x = 2, bc_particle_y = 2, bc_particle_z = 2, particle_mode = 0, particle_rho_f = %(rf)s,\n/\n'

BASES = {  # name -> (case dir, density of the liquid, gravity)
    'tank': ('vof_slosh_small', '1000.0', '9.81'),
    'ibm': ('vof_ibm_small', '1.0', '0.0'),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--exe', required=True)
    ap.add_argument('--mpirun', default='mpirun')
    ap.add_argument('--np', type=int, default=1)
    a = ap.parse_args()
    tmp = tempfile.mkdtemp(prefix='vofcombo_')
    bad = 0
    try:
        for flow, base, feats in itertools.product((0, 1), BASES, ('sed', 'bou', 'par', 'all')):
            src, rf, g = BASES[base]
            if flow == 0 and base == 'ibm':
                continue          # the passive IBM case needs its own prescribed flow; the single-phase IBM + features are covered elsewhere
            wd = os.path.join(tmp, '%s_f%d_%s' % (base, flow, feats))
            shutil.copytree(os.path.join(REG, src), wd)
            for d in ('restart', 'stats', 'fields'):
                os.makedirs(os.path.join(wd, d), exist_ok=True)
            inp = os.path.join(wd, 'input_parameters')
            t = open(inp).read()
            t = re.sub(r'nsteps\s*=\s*-?\d+', 'nsteps = 8', t)
            t = re.sub(r'nsave\s*=\s*-?\d+', 'nsave = 8', t)
            t = re.sub(r'vof_flow\s*=\s*\d', 'vof_flow = %d' % flow, t)
            if base == 'ibm':
                t = re.sub(r'vof_grav\s*=\s*[\d.]+', 'vof_grav = 0.0', t)
            sub = {'rf': rf, 'g': g}
            t += (SED % sub if feats in ('sed', 'all') else '') + (BOU % sub if feats in ('bou', 'all') else '') \
                + (PAR % sub if feats in ('par', 'all') else '')
            if base == 'ibm' and feats in ('par', 'all'):
                t = t.replace('particle_seed_ymax = 0.4', 'particle_seed_ymax = 1.9')
            open(inp, 'w').write(t)
            r = subprocess.run([a.mpirun, '--oversubscribe', '-np', str(a.np), a.exe], cwd=wd, capture_output=True, text=True)
            out = r.stdout + r.stderr
            ok = r.returncode == 0 and 'NaN' not in out and 'Infinity' not in out
            msg = ''
            if ok and feats in ('par', 'all'):
                pf = os.path.join(wd, 'particles_restart')
                if os.path.exists(pf):
                    sys.path.insert(0, REG)
                    from compare_particles import read_particles
                    n = len(read_particles(pf))
                    if n != 60:
                        ok = False
                        msg = ' particle count %d' % n
            print('%-6s vof_flow=%d %-4s %s%s' % (base, flow, feats, 'ok' if ok else 'FAIL', msg))
            if not ok:
                print(out[-800:])
            bad += not ok
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    sys.exit(1 if bad else 0)


main()
