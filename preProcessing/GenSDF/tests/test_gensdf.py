#!/usr/bin/env python3
"""Regression tests for GenSDF's solid/fluid classification on damaged geometry (needs numpy and mpirun).

Builds small analytic shapes (sphere, box, torus, hollow body, terrain) in Python, damages them the way real meshes are
damaged (holes, missing faces, cracks, flipped or duplicated or overlapping triangles, scattered deleted triangles),
runs ../gensdf on a 64^3 grid and counts the cells whose sign disagrees with the intended solid (cells within ~1.5 cells of
the surface are not scored).  Also checks rank-count independence of the whole output, the OBJ/ASCII-STL readers,
vertical_axis = 2, multi-solid object IDs and the `ground` keyword.

    cd preProcessing/GenSDF && make && python3 tests/test_gensdf.py [--exe ./gensdf] [--np 1 4] [--keep]
"""
import argparse, os, shutil, struct, subprocess, sys, tempfile
import numpy as np

N = (64, 64, 64)
DX = 1.0 / 64
C = (0.5, 0.5, 0.5)
R = 0.3

# ------------------------------------------------------------------------------------------------ geometry builders
def icosphere(c, r, sub=3):
    t = (1 + 5 ** 0.5) / 2
    v = np.array([[-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0], [0, -1, t], [0, 1, t], [0, -1, -t], [0, 1, -t],
                  [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1]], float)
    f = np.array([[0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11], [1, 5, 9], [5, 11, 4], [11, 10, 2],
                  [10, 7, 6], [7, 1, 8], [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9], [4, 9, 5], [2, 4, 11],
                  [6, 2, 10], [8, 6, 7], [9, 8, 1]])
    v = v / np.linalg.norm(v, axis=1)[:, None]
    for _ in range(sub):
        cache, nv, nf = {}, list(v), []
        def mid(a, b):
            k = (min(a, b), max(a, b))
            if k not in cache:
                m = v[a] + v[b]; m /= np.linalg.norm(m); cache[k] = len(nv); nv.append(m)
            return cache[k]
        for a, b, cc in f:
            ab, bc, ca = mid(a, b), mid(b, cc), mid(cc, a)
            nf += [[a, ab, ca], [b, bc, ab], [cc, ca, bc], [ab, bc, ca]]
        v, f = np.array(nv), np.array(nf)
    return v[f] * r + np.asarray(c)

def box(lo, hi, skip=()):
    X, Y, Z = [lo[0], hi[0]], [lo[1], hi[1]], [lo[2], hi[2]]
    P = lambda i, j, k: np.array([X[i], Y[j], Z[k]])
    quads = {'-x': [P(0,0,0), P(0,0,1), P(0,1,1), P(0,1,0)], '+x': [P(1,0,0), P(1,1,0), P(1,1,1), P(1,0,1)],
             '-y': [P(0,0,0), P(1,0,0), P(1,0,1), P(0,0,1)], '+y': [P(0,1,0), P(0,1,1), P(1,1,1), P(1,1,0)],
             '-z': [P(0,0,0), P(0,1,0), P(1,1,0), P(1,0,0)], '+z': [P(0,0,1), P(1,0,1), P(1,1,1), P(0,1,1)]}
    tri = []
    for k, q in quads.items():
        if k not in skip:
            tri += [[q[0], q[1], q[2]], [q[0], q[2], q[3]]]
    return np.array(tri)

def torus(Rm=0.25, r=0.1, nu=48, nv=24):
    u = np.linspace(0, 2*np.pi, nu, endpoint=False); v = np.linspace(0, 2*np.pi, nv, endpoint=False)
    P = lambda a, b: np.array([C[0]+(Rm+r*np.cos(b))*np.cos(a), C[1]+(Rm+r*np.cos(b))*np.sin(a), C[2]+r*np.sin(b)])
    tri = []
    for i in range(nu):
        for j in range(nv):
            p00, p10 = P(u[i], v[j]), P(u[(i+1) % nu], v[j])
            p11, p01 = P(u[(i+1) % nu], v[(j+1) % nv]), P(u[i], v[(j+1) % nv])
            tri += [[p00, p10, p11], [p00, p11, p01]]
    return np.array(tri)

def terrain(n=40, amp=0.08, base=0.3):
    xs = np.linspace(0, 1, n + 1)
    h = lambda x, y: base + amp*np.sin(2*np.pi*x)*np.cos(2*np.pi*y)
    tri = []
    for i in range(n):
        for j in range(n):
            p = lambda a, b: np.array([xs[a], xs[b], h(xs[a], xs[b])])
            tri += [[p(i, j), p(i+1, j), p(i+1, j+1)], [p(i, j), p(i+1, j+1), p(i, j+1)]]
    return np.array(tri), h

def cut_cone(tri, c, d, half_angle):
    d = np.asarray(d, float); d /= np.linalg.norm(d)
    cen = tri.mean(axis=1) - np.asarray(c)
    return tri[(cen @ d) / np.linalg.norm(cen, axis=1) < np.cos(np.radians(half_angle))]

def flip_random(tri, frac, seed=1):
    m = np.random.default_rng(seed).random(len(tri)) < frac
    t = tri.copy(); t[m] = t[m][:, [0, 2, 1], :]; return t

def flip_all(tri): return tri[:, [0, 2, 1], :]
def delete_random(tri, frac, seed=3): return tri[np.random.default_rng(seed).random(len(tri)) >= frac]

def write_stl(path, tri):
    n = np.cross(tri[:, 1]-tri[:, 0], tri[:, 2]-tri[:, 0]); l = np.linalg.norm(n, axis=1); l[l == 0] = 1; n /= l[:, None]
    with open(path, 'wb') as f:
        f.write(b'\0'*80); f.write(struct.pack('<I', len(tri)))
        for i in range(len(tri)):
            f.write(struct.pack('<12fH', *n[i], *tri[i].ravel(), 0))

def write_ascii_stl(path, tri):
    n = np.cross(tri[:, 1]-tri[:, 0], tri[:, 2]-tri[:, 0]); n /= np.linalg.norm(n, axis=1)[:, None]
    with open(path, 'w') as f:
        f.write('solid s\n')
        for t, nn in zip(tri, n):
            f.write('facet normal %e %e %e\n outer loop\n' % tuple(nn))
            for p in t: f.write('  vertex %e %e %e\n' % tuple(p))
            f.write(' endloop\nendfacet\n')
        f.write('endsolid s\n')

def write_quad_obj_cube(path, lo=(0.3,)*3, hi=(0.7,)*3):          # quads, no vertex normals
    with open(path, 'w') as f:
        for z in (lo[2], hi[2]):
            for y in (lo[1], hi[1]):
                for x in (lo[0], hi[0]): f.write('v %f %f %f\n' % (x, y, z))
        for q in [(1,3,4,2), (5,6,8,7), (1,2,6,5), (3,7,8,4), (1,5,7,3), (2,4,8,6)]: f.write('f %d %d %d %d\n' % q)

# ------------------------------------------------------------------------------------------------ truth and harness
def grid(n=N, L=(1, 1, 1)):
    return np.meshgrid(*[(np.arange(m) + 0.5) * l / m for m, l in zip(n, L)], indexing='ij')

def sphere_truth(holes=(), c=C, r=R):
    def f(X, Y, Z):
        rr = np.sqrt((X-c[0])**2 + (Y-c[1])**2 + (Z-c[2])**2)
        ins, band = rr < r, np.abs(rr - r) < 1.5*DX
        for d, a in holes:                     # a planar lid across every hole: what lies beyond the rim plane is outside
            d = np.asarray(d, float); d /= np.linalg.norm(d)
            h = (X-c[0])*d[0] + (Y-c[1])*d[1] + (Z-c[2])*d[2]; cut = r*np.cos(np.radians(a))
            band = band | ((np.abs(h - cut) < 2*DX) & (rr < r)); ins = ins & (h <= cut)
        return ins, band
    return f

def box_truth(lo=(0.25, 0.3, 0.2), hi=(0.7, 0.65, 0.6)):
    def f(X, Y, Z):
        ins = (X > lo[0]) & (X < hi[0]) & (Y > lo[1]) & (Y < hi[1]) & (Z > lo[2]) & (Z < hi[2])
        dd = np.minimum.reduce([np.abs(X-lo[0]), np.abs(X-hi[0]), np.abs(Y-lo[1]), np.abs(Y-hi[1]), np.abs(Z-lo[2]), np.abs(Z-hi[2])])
        return ins, dd < 0.6*DX
    return f

def hollow_truth(rin):
    def f(X, Y, Z):
        rr = np.sqrt((X-C[0])**2 + (Y-C[1])**2 + (Z-C[2])**2)
        return (rr < R) & (rr > rin), (np.abs(rr-R) < 1.5*DX) | (np.abs(rr-rin) < 1.5*DX)
    return f

def torus_truth(X, Y, Z, Rm=0.25, r=0.1):
    d = np.sqrt((np.sqrt((X-C[0])**2 + (Y-C[1])**2) - Rm)**2 + (Z-C[2])**2)
    return d < r, np.abs(d - r) < 2.0*DX

def union_truth(spheres):
    def f(X, Y, Z):
        ins = np.zeros_like(X, bool); band = np.zeros_like(X, bool)
        for c, r in spheres:
            rr = np.sqrt((X-c[0])**2 + (Y-c[1])**2 + (Z-c[2])**2); ins |= rr < r; band |= np.abs(rr - r) < 1.5*DX
        return ins, band
    return f

def setup(cd, solids, n=N, L=(1, 1, 1), vaxis=3, faces=False, extra=()):
    d = os.path.join(cd, 'data'); os.makedirs(d, exist_ok=True)
    lines = []
    for s in solids:
        name, geom = s[0], s[1]
        if callable(geom): geom(os.path.join(d, name))
        else: write_stl(os.path.join(d, name), geom)
        lines.append('data/%s%s' % (name, ' ground' if len(s) > 2 and s[2] == 'ground' else ''))
    open(os.path.join(d, 'objects.list'), 'w').write('\n'.join(lines) + '\n')
    open(os.path.join(d, 'geometry.out'), 'w').write('# n\n%d %d %d\n# l\n%g %g %g\n' % (*n, *L))
    p = ["'data/objects.list'", "1d10 5 40", "%d %d %d" % n, "0 0 0", ".false.", ".true.", "5", str(vaxis),
         ".true." if faces else ".false."] + list(extra)
    open(os.path.join(cd, 'parameters.in'), 'w').write('\n'.join(p) + '\n')

def run(exe, cd, np_):
    r = subprocess.run(['mpirun', '--oversubscribe', '-np', str(np_), os.path.abspath(exe)], cwd=cd, capture_output=True, text=True)
    open(os.path.join(cd, 'log.txt'), 'w').write(r.stdout + r.stderr)
    if r.returncode != 0: raise RuntimeError('gensdf failed (np=%d):\n%s' % (np_, (r.stdout + r.stderr)[-800:]))

def read(cd, name='sdfp.bin', n=N):                                      # -> internal (x, y, z=vertical) order
    a = np.fromfile(os.path.join(cd, 'data', name), '>f8').reshape((n[0]+2, n[2]+2, n[1]), order='F')
    return a[1:-1, 1:-1, :].transpose(0, 2, 1)

def wrong(s, truth):
    X, Y, Z = grid(); ins, band = truth(X, Y, Z)
    return int((((s < 0) != ins) & ~band).sum())

# ------------------------------------------------------------------------------------------------ the cases
BOX = lambda **k: box((0.25, 0.3, 0.2), (0.7, 0.65, 0.6), **k)
SPH = lambda: icosphere(C, R, 3)
def holes_sphere(specs):
    t = SPH()
    for d, a in specs: t = cut_cone(t, C, d, a)
    return t
def hemispheres_gap(gap):
    t = SPH(); cen = t.mean(axis=1); up, dn = t[cen[:, 2] >= C[2]].copy(), t[cen[:, 2] < C[2]].copy()
    up[:, :, 2] += gap/2; dn[:, :, 2] -= gap/2; return np.concatenate([up, dn])
def hollow(rin=0.15): return np.concatenate([SPH(), flip_all(icosphere(C, rin, 3))])
def ring_crack(t, w): return t[np.abs(t.mean(axis=1)[:, 2] - C[2]) > w*R]

# name -> (triangles, truth, max wrong cells allowed)
CASES = {
 'sphere closed':               (SPH, sphere_truth(), 0),
 'sphere, mixed orientation':   (lambda: flip_random(SPH(), 0.3), sphere_truth(), 0),
 'sphere, all inverted':        (lambda: flip_all(SPH()), sphere_truth(), 0),
 'sphere, one 25deg hole':      (lambda: holes_sphere([((0, 0, 1), 25)]), sphere_truth([((0, 0, 1), 25)]), 2),
 'sphere, one 60deg hole':      (lambda: holes_sphere([((0, 0, 1), 60)]), sphere_truth([((0, 0, 1), 60)]), 2),
 'sphere, six 20deg holes':     (lambda: holes_sphere([((1,0,0),20), ((-1,0,0),20), ((0,1,0),20), ((0,-1,0),20), ((0,0,1),20), ((0,0,-1),20)]),
                                 sphere_truth([((1,0,0),20), ((-1,0,0),20), ((0,1,0),20), ((0,-1,0),20), ((0,0,1),20), ((0,0,-1),20)]), 2),
 'sphere, holes + flips':       (lambda: flip_random(holes_sphere([((0, 0, 1), 25)]), 0.3), sphere_truth([((0, 0, 1), 25)]), 2),
 'sphere, duplicated shell':    (lambda: np.concatenate([SPH(), SPH()]), sphere_truth(), 0),
 'sphere, ring crack':          (lambda: ring_crack(SPH(), 0.03), sphere_truth(), 2),
 'sphere, 0.4-cell crack':      (lambda: hemispheres_gap(0.4*DX), sphere_truth(), 2),
 'sphere, 5% triangles lost':   (lambda: delete_random(SPH(), 0.05), sphere_truth(), 2),
 'sphere, 15% triangles lost':  (lambda: delete_random(SPH(), 0.15), sphere_truth(), 40),
 'box closed':                  (BOX, box_truth(), 0),
 'box, no floor':               (lambda: BOX(skip=('-z',)), box_truth(), 0),
 'box, no top and no floor':    (lambda: BOX(skip=('+z', '-z')), box_truth(), 0),
 'box, no top and one side':    (lambda: BOX(skip=('+z', '-x')), box_truth(), 0),
 'box, three faces missing':    (lambda: BOX(skip=('-z', '-x', '-y')), box_truth(), 0),
 'two overlapping spheres':     (lambda: np.concatenate([icosphere((0.42, 0.5, 0.5), 0.22), icosphere((0.6, 0.5, 0.5), 0.22)]),
                                 union_truth([((0.42, 0.5, 0.5), 0.22), ((0.6, 0.5, 0.5), 0.22)]), 0),
 'hollow body, scrambled':      (lambda: flip_random(hollow(), 0.3), hollow_truth(0.15), 0),
 'torus':                       (torus, torus_truth, 0),
 'torus, 8% triangles lost':    (lambda: delete_random(torus(), 0.08), torus_truth, 10),
}

def run_cases(exe, np_, work, names):
    bad = 0
    for name in names:
        mk, truth, allowed = CASES[name]
        cd = os.path.join(work, 'np%d_%s' % (np_, name.replace(' ', '_').replace(',', '').replace('%', 'pct')))
        setup(cd, [('obj.stl', mk())]); run(exe, cd, np_)
        w = wrong(read(cd), truth); ok = w <= allowed; bad += (not ok)
        print('  %-30s np=%d  wrong cells %5d (allowed %d)  %s' % (name, np_, w, allowed, 'ok' if ok else 'FAIL'))
    return bad

def test_rank_independence(exe, work, nps):
    """the written fields must not depend on the number of ranks (classification and fast sweep are both slab-parallel).
    Signs, object ids and every distance outside the converging-front region at the very centre of a body agree to round-off;
    inside that region a multi-rank run takes one more confirming sweep round than a serial one, hence the 1e-5 tolerance."""
    ref = None; bad = 0
    for np_ in nps:
        cd = os.path.join(work, 'par%d' % np_)
        setup(cd, [('obj.stl', flip_random(holes_sphere([((0, 0, 1), 25)]), 0.2))], faces=True); run(exe, cd, np_)
        cur = [np.fromfile(os.path.join(cd, 'data', f), '>f8') for f in ('sdfp.bin', 'sdfu.bin', 'sdfv.bin', 'sdfw.bin', 'sdfp_objid.bin')]
        if ref is None: ref = cur
        d = max(np.abs(a - b).max() for a, b in zip(ref, cur)); ok = d < 1e-5; bad += (not ok)
        print('  all five fields at np=%d vs np=%d: max difference %.1e  %s' % (np_, nps[0], d, 'ok' if ok else 'FAIL'))
    return bad

def test_ground_and_multisolid(exe, work, np_):
    bad = 0
    tri, h = terrain()
    cd = os.path.join(work, 'terrain'); setup(cd, [('t.stl', tri, 'ground')]); run(exe, cd, np_)
    X, Y, Z = grid(); s = read(cd); w = int((((s < 0) != (Z < h(X, Y))) & ~(np.abs(Z - h(X, Y)) < 1.5*DX)).sum())
    print('  ground sheet (everything below is solid, incl. below its lowest vertex): wrong %d  %s' % (w, 'ok' if w == 0 else 'FAIL')); bad += (w != 0)
    s1 = icosphere((0.3, 0.55, 0.6), 0.15); s2 = cut_cone(icosphere((0.7, 0.55, 0.6), 0.15), (0.7, 0.55, 0.6), (0, 0, 1), 30)
    tri, h = terrain(base=0.2, amp=0.05)
    cd = os.path.join(work, 'multi'); setup(cd, [('a.stl', s1), ('b.stl', s2), ('t.stl', tri, 'ground')]); run(exe, cd, np_)
    s, oid = read(cd), read(cd, 'sdfp_objid.bin')
    r1 = np.sqrt((X-.3)**2 + (Y-.55)**2 + (Z-.6)**2); r2 = np.sqrt((X-.7)**2 + (Y-.55)**2 + (Z-.6)**2)
    in2 = (r2 < .15) & ((Z-.6) <= .15*np.cos(np.radians(30))); ins = (r1 < .15) | in2 | (Z < h(X, Y))
    band = (np.abs(r1-.15) < 1.5*DX) | (np.abs(r2-.15) < 1.5*DX) | (np.abs(Z-h(X, Y)) < 1.5*DX) | (np.abs((Z-.6) - .15*np.cos(np.radians(30))) < 2*DX)
    w = int((((s < 0) != ins) & ~band).sum())
    exp = np.where(r1 < .15, 1, np.where(in2, 2, np.where(Z < h(X, Y), 3, 0)))
    wid = int(((np.rint(oid) != exp) & (s < 0) & ins & ~band & (s > -4*DX)).sum())
    print('  3 solids (closed sphere, holed sphere, ground): wrong cells %d, wrong object ids near surface %d  %s' % (w, wid, 'ok' if w == wid == 0 else 'FAIL'))
    return bad + (w != 0 or wid != 0)

def test_formats_and_axis(exe, work, np_):
    bad = 0
    cd = os.path.join(work, 'objquad'); setup(cd, [('c.obj', write_quad_obj_cube)]); run(exe, cd, np_)
    X, Y, Z = grid(); s = read(cd); ins = (X > .3) & (X < .7) & (Y > .3) & (Y < .7) & (Z > .3) & (Z < .7)
    band = np.minimum.reduce([np.abs(X-.3), np.abs(X-.7), np.abs(Y-.3), np.abs(Y-.7), np.abs(Z-.3), np.abs(Z-.7)]) < .6*DX
    w = int((((s < 0) != ins) & ~band).sum()); print('  OBJ with quad faces and no vertex normals: wrong %d  %s' % (w, 'ok' if w == 0 else 'FAIL')); bad += (w != 0)
    cd = os.path.join(work, 'asciistl'); setup(cd, [('s.stl', lambda p: write_ascii_stl(p, SPH()))]); run(exe, cd, np_)
    w = wrong(read(cd), sphere_truth()); print('  ASCII STL: wrong %d  %s' % (w, 'ok' if w == 0 else 'FAIL')); bad += (w != 0)
    # vertical_axis = 2: sphere at solver (0.5, 0.3, 0.25) in a 1 x 0.75 x 0.5 box, 64 x 48 x 32 cells, all four staggerings
    nx, ny, nz = 64, 48, 32; c = (0.5, 0.3, 0.25); rr = 0.2
    cd = os.path.join(work, 'vaxis2'); setup(cd, [('s.stl', icosphere(c, rr))], n=(nx, ny, nz), L=(1.0, 0.75, 0.5), vaxis=2, faces=True); run(exe, cd, np_)
    dxs, dys, dzs = 1/nx, 0.75/ny, 0.5/nz; res = {}
    for f, (a, b, d) in {'sdfp.bin': (0, 0, 0), 'sdfu.bin': (.5, 0, 0), 'sdfv.bin': (0, .5, 0), 'sdfw.bin': (0, 0, .5)}.items():
        arr = np.fromfile(os.path.join(cd, 'data', f), '>f8').reshape((nx+2, ny+2, nz), order='F')[1:-1, 1:-1, :]
        XX, YY, ZZ = np.meshgrid((np.arange(nx)+.5+a)*dxs, (np.arange(ny)+.5+b)*dys, (np.arange(nz)+.5+d)*dzs, indexing='ij')
        q = np.sqrt((XX-c[0])**2 + (YY-c[1])**2 + (ZZ-c[2])**2)
        res[f] = int((((arr < 0) != (q < rr)) & ~(np.abs(q-rr) < 1.5*max(dxs, dys, dzs))).sum())
    print('  vertical_axis = 2, all four staggerings: wrong %s  %s' % (res, 'ok' if not any(res.values()) else 'FAIL')); bad += any(res.values())
    return bad

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--exe', default=os.path.join(os.path.dirname(__file__), '..', 'gensdf'))
    ap.add_argument('--np', type=int, nargs='+', default=[1, 4], help='rank counts; the first one runs every case')
    ap.add_argument('--keep', action='store_true', help='keep the work directory')
    ap.add_argument('--quick', action='store_true', help='skip the 15%% deletion / torus cases')
    a = ap.parse_args()
    if not os.path.exists(a.exe): sys.exit('gensdf not built: run make in preProcessing/GenSDF first (or pass --exe)')
    work = tempfile.mkdtemp(prefix='gensdf_test_'); bad = 0
    try:
        names = [n for n in CASES if not (a.quick and ('15%' in n or 'torus' in n))]
        print('classification of damaged meshes:'); bad += run_cases(a.exe, a.np[0], work, names)
        for np_ in a.np[1:]:
            print('same, subset on more ranks:'); bad += run_cases(a.exe, np_, work, ['sphere, six 20deg holes', 'box, three faces missing', 'two overlapping spheres'])
        print('rank independence:'); bad += test_rank_independence(a.exe, work, a.np if len(a.np) > 1 else [1, 3])
        print('ground / multi-solid:'); bad += test_ground_and_multisolid(a.exe, work, a.np[0])
        print('readers / vertical axis:'); bad += test_formats_and_axis(a.exe, work, a.np[0])
    finally:
        if a.keep: print('work directory kept:', work)
        else: shutil.rmtree(work, ignore_errors=True)
    print('\n%s' % ('ALL PASSED' if bad == 0 else '%d FAILED' % bad)); sys.exit(1 if bad else 0)

if __name__ == '__main__':
    main()
