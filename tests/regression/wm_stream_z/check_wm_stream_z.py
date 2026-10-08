#!/usr/bin/env python3
"""Wall-modelled spanwise (z) walls (single phase): a uniform stream along x decays by the log-law (Reichardt) wall stress u_tau^2 of the
first-cell velocity at the z walls (--duct: at the y walls too), identically on 1 and 4 ranks (z split; x split for the duct).
--mode v|diag|prof|prof_stretch: the stream is imposed through a restart file built from a one-step seed run, on 1 and 4 ranks (z split),
and the result must agree between the two and with the log-law integration: v/diag a uniform stream along y / the diagonal (mean decay,
exercises tau_zv), prof/prof_stretch an unequal, non-uniform (U,V) profile whose first cell at each wall decays by its own log-law stress
(exercises the wall-owning rank and, for prof_stretch, the stretched z grid: first-cell height and width).
--reject-rough: the rough EQWM with z walls must be refused at input."""
import argparse, math, os, re, shutil, struct, subprocess, sys, tempfile
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('--exe', required=True)
ap.add_argument('--mpirun', default='mpirun')
ap.add_argument('--duct', action='store_true', help='y walls as well (4-wall duct)')
ap.add_argument('--reject-rough', action='store_true')
ap.add_argument('--mode', choices=['v', 'diag', 'prof', 'prof_stretch'])
ap.add_argument('--tol', type=float, default=0.05, help='relative tolerance on the velocity deficit')
a = ap.parse_args()
src = os.path.dirname(os.path.abspath(__file__))
tmp = tempfile.mkdtemp(prefix='wmstreamz_')


def run(tag, np_, edits, fail_ok=False, restart=None, snap=False):
    wd = os.path.join(tmp, tag)
    os.makedirs(os.path.join(wd, 'restart'))
    if restart:
        write_snapshot(os.path.join(wd, 'restart', 'x'), *restart)
    text = open(os.path.join(src, 'input_parameters')).read()
    for pat, rep in edits:
        text = re.sub(pat, rep, text, count=1)
    open(os.path.join(wd, 'input_parameters'), 'w').write(text)
    r = subprocess.run([a.mpirun, '--oversubscribe', '-np', str(np_), a.exe], cwd=wd, capture_output=True, text=True, timeout=900)
    if fail_ok:
        return r
    if r.returncode != 0:
        print(r.stdout[-2000:], r.stderr[-2000:]); sys.exit(1)
    if snap:
        return read_snapshot(os.path.join(wd, 'fields'))
    rows = [l.split() for l in r.stdout.splitlines() if re.match(r'^\s*\d+\s+[-\d.E+]+\s+[-\d.E+]+\s+[-\d.E+]+\s', l)]
    return int(rows[-1][0]), float(rows[-1][1]), float(rows[-1][2])


def uplus(yp):
    k, c = 0.41, 5.2 - math.log(0.41)/0.41
    return math.log(1 + k*yp)/k + c*(1 - math.exp(-yp/11) - yp/11*math.exp(-0.33*yp))


def u_tau(u, y, nu):
    lo, hi = 1e-6*u, u
    for _ in range(200):
        m = 0.5*(lo + hi)
        lo, hi = (m, hi) if u/m - uplus(m*y/nu) > 0 else (lo, m)
    return 0.5*(lo + hi)


def read_snapshot(d):
    """Newest big-endian snapshot (== restart file): six mesh arrays, then U, V, W, ... as Fortran-ordered (nx,ny,nz) blocks."""
    f = max((x for x in os.listdir(d) if x.startswith('x.')), key=lambda x: int(x.split('.')[1]))
    b = open(os.path.join(d, f), 'rb').read()
    p, mesh, fields = 0, [], []
    for _ in range(6):
        n, = struct.unpack_from('>i', b, p)
        mesh.append(np.frombuffer(b, '>f8', n, p + 4).astype(float)); p += 4 + 8*n
    while p < len(b):
        sh = struct.unpack_from('>3i', b, p)
        fields.append(np.frombuffer(b, '>f8', sh[0]*sh[1]*sh[2], p + 12).astype(float).reshape(sh, order='F')); p += 12 + 8*sh[0]*sh[1]*sh[2]
    return mesh, fields


def write_snapshot(path, mesh, fields):
    with open(path, 'wb') as fh:
        for m in mesh:
            fh.write(struct.pack('>i', len(m)) + m.astype('>f8').tobytes())
        for f in fields:
            fh.write(struct.pack('>3i', *f.shape) + f.astype('>f8').tobytes(order='F'))


def interior(f):
    """Interior points: ghost layers of the cell-centred dimensions dropped (the seam ghosts depend on the rank layout)."""
    return f[1:-1, 1:-1, 1:-1]


def euler(vel, y, dz, nsteps, dt, nu):
    """Explicit log-law decay of the velocity vector (u,v) of one cell of height y and width dz."""
    u, v = vel
    for _ in range(nsteps):
        s = math.hypot(u, v)
        ut2 = u_tau(s, y, nu)**2
        u, v = u - dt*ut2*u/s/dz, v - dt*ut2*v/s/dz
    return np.array([u, v])


def stream_modes(mode):
    base = [(r'nsteps = 150, nsave = 100000', 'nsteps = 1, nsave = 1,')]
    if mode == 'prof_stretch':
        base.append((r'alpha_grid = 1.0,', 'alpha_grid = 1.0, alpha_grid_z = 2.0,'))
    mesh, f = run('seed', 1, base, snap=True)
    z, zc = mesh[2], mesh[5]
    zg = np.concatenate([zc[:1], zc, zc[-1:]])[None, None, :]
    lz, nsteps, dt, nu = z[-1], 150, 2e-3, 1e-6
    if mode == 'v':
        U, V = 0*zg, 1 + 0*zg
    elif mode == 'diag':
        U, V = (1/math.sqrt(2)) + 0*zg, (1/math.sqrt(2)) + 0*zg
    else:
        U, V = 0.5 + zg/lz, 1.0 - 0.6*zg/lz
    fields = [np.broadcast_to(U, f[0].shape).copy(), np.broadcast_to(V, f[1].shape).copy()] + f[2:]
    ed = [(r'nsteps = 150, nsave = 100000', 'nsteps = 150, nsave = 150,'), (r'restart = 0', 'restart = 1, nstep_init = 0')]
    ed += base[1:]
    e4 = ed + [(r'p_row = 1, p_col = 1', 'p_row = 1, p_col = 4')]
    out = {}
    for tag, np_, e in (('np1', 1, ed), ('np4', 4, e4)):
        out[tag] = run(tag, np_, e, restart=(mesh, fields), snap=True)[1]
    d = max(np.abs(interior(p) - interior(q)).max() for p, q in zip(out['np1'][:3], out['np4'][:3]))
    Uo, Vo = (interior(out['np1'][i]).mean(axis=(0, 1)) for i in (0, 1))
    wmax = np.abs(interior(out['np1'][2])).max()
    msg = 'np1-np4 max diff %.2e  max|W| %.2e' % (d, wmax)
    ok = d < 1e-9 and wmax < 1e-8
    if mode in ('v', 'diag'):
        ref0 = np.array([U.flat[0], V.flat[0]])
        vel = euler(ref0, zc[0], lz/2, nsteps, dt, nu)
        um, vm = (np.sum(p*np.diff(z))/lz for p in (Uo, Vo))
        ok = ok and max(abs(um - vel[0]), abs(vm - vel[1])) <= a.tol*np.hypot(*(ref0 - vel))
        if mode == 'v':
            ok = ok and abs(um) < 1e-8
        msg += '  mean (U,V) run (%.6f, %.6f)  log-law (%.6f, %.6f)' % (um, vm, vel[0], vel[1])
    else:
        for name, k, y, dz in (('lo', 0, zc[0], z[1] - z[0]), ('hi', -1, lz - zc[-1], z[-1] - z[-2])):
            ref0 = np.array([U[0, 0, k], V[0, 0, k]])
            vel = euler(ref0, y, dz, nsteps, dt, nu)
            got = np.array([Uo[k], Vo[k]])
            good = np.hypot(*(got - vel)) <= a.tol*np.hypot(*(ref0 - vel))
            ok = ok and good
            msg += '\n  %s wall cell (U,V) run (%.6f, %.6f)  log-law (%.6f, %.6f)  %s' % (name, got[0], got[1], vel[0], vel[1], 'ok' if good else 'FAIL')
    print('%s: %s  %s' % (mode, msg, 'ok' if ok else 'FAIL'))
    sys.exit(0 if ok else 1)


try:
    if a.mode:
        stream_modes(a.mode)
    if a.reject_rough:
        r = run('rough', 1, [(r'flat_wall_model_flag = 1', 'flat_wall_model_flag = 2, z0_ylo = 1e-4, z0_yhi = 1e-4')], fail_ok=True)
        ok = r.returncode != 0 and 'not yet supported for z walls' in r.stdout + r.stderr
        print('rough EQWM with z walls %s' % ('rejected' if ok else 'NOT rejected: FAIL'))
        sys.exit(0 if ok else 1)
    if a.duct:
        e1 = [(r'y_bc_type = 0', 'y_bc_type = 1'), (r'ny = 18', 'ny = 34')]
        e4 = e1 + [(r'p_row = 1, p_col = 1', 'p_row = 4, p_col = 1')]
        nwall = 4
    else:
        e1 = []
        e4 = [(r'p_row = 1, p_col = 1', 'p_row = 1, p_col = 4')]
        nwall = 2
    step, t, u_np1 = run('np1', 1, e1)
    _, _, u_np4 = run('np4', 4, e4)
    l, y_ref, nu, dt, u = 1.0, 0.5*1.0/33, 1e-6, 2e-3, 1.0
    for _ in range(step):
        u -= dt*nwall*u_tau(u, y_ref, nu)**2/l
    ok = abs(u_np1 - u) <= a.tol*(1 - u) and abs(u_np4 - u_np1) < 1e-9
    print('step %d t=%.3f  mean U: run %.6f  log-law %.6f  np4 %.6f  %s' % (step, t, u_np1, u, u_np4, 'ok' if ok else 'FAIL'))
    sys.exit(0 if ok else 1)
finally:
    shutil.rmtree(tmp, ignore_errors=True)
