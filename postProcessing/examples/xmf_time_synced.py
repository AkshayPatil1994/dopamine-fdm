#!/usr/bin/env python3
"""Fields and particles in ParaView on one shared physical-time axis (single-phase, Boussinesq, sediment or two-fluid VOF; particles optional).

Run:   python xmf_time_synced.py [case_dir] [--stride N] [--first STEP] [--last STEP] [--sdf FILE|none]       (default case_dir = .)
Writes in <case>/paraview/:
  xm.bin, ym.bin, zm.bin      cell-centre coordinates (the only data written; the solver files are not copied)
  <prefix>.xmf                fields  Velocity, P, [C], [nu_t], [T], [alpha], [sdf]   from fields/<prefix>.<step>
  <prefix>_particles.xmf      particles (point cloud: Velocity, Acceleration, age, id) from fields/<prefix>_particles.<step>, if present
Both files list the same snapshots with the same Time values, so the ParaView time toolbar plays them together. Open both with the XDMF
Reader; for VOF contour alpha = 0.5 (free surface) and sdf = 0 (IBM bed), and show the particles with Glyph or Point Gaussian.

Time value: the monitor line of <case>/run.log just before each "Writing snapshot: fields/<prefix>.<step>" line, which stays right after a
restart (the monitor step counter restarts at 0, the file names do not). Without a log: index*tsave when nsave < 0, else the step number.
Which optional blocks a snapshot holds (C, nu_t, T, alpha) is read from <case>/input_parameters; a missing sgs_model means the default
AMD, so nu_t is present. Unlike dopamine_post.fields.write_field_xmf this script is standalone (numpy only) and uses physical time.
Velocity is the lower-face value of the staggered components, ghost layers stripped, no averaging.
"""
import argparse, os, re, struct
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('case', nargs='?', default=os.getcwd())
ap.add_argument('--stride', type=int, default=1)
ap.add_argument('--first', type=int, default=None)
ap.add_argument('--last', type=int, default=None)
ap.add_argument('--sdf', default=None, help='signed-distance file (default: ibm_sdf_file of input_parameters); "none" skips it')
args = ap.parse_args()
CASE = os.path.abspath(args.case)
FIELDS, OUT = os.path.join(CASE, 'fields'), os.path.join(CASE, 'paraview')
os.makedirs(OUT, exist_ok=True)

text = re.sub(r'!.*', '', open(os.path.join(CASE, 'input_parameters')).read())


def get(name, default=0, cast=int):
    m = re.search(rf'\b{name}\s*=\s*([^,/\s]+)', text, re.IGNORECASE)
    return cast(m.group(1).replace('d', 'e').replace('D', 'e').strip("'\"")) if m else default


SEDIMENT, SGS, BOUSS = get('sediment_flag') >= 1, get('sgs_model', 2) != 0, get('boussinesq_flag') >= 1
VOF = get('vof_active') >= 1
NSAVE, TSAVE = get('nsave', 0), get('tsave', 0.0, float)
SDF_FILE = re.search(r"ibm_sdf_file\s*=\s*'([^']+)'", text, re.IGNORECASE)
SDF_FILE = os.path.join(CASE, SDF_FILE.group(1)) if SDF_FILE and get('ibm_input_mode') >= 1 else None
if args.sdf:
    SDF_FILE = None if args.sdf.lower() == 'none' else os.path.join(CASE, args.sdf)

pat = re.compile(r'^([A-Za-z0-9_]+)\.(\d+)$')
found = {}
for f in os.listdir(FIELDS):
    m = pat.match(f)
    if m and not f.startswith(('grid', 'geometry')) and not m.group(1).endswith(('_particles', '_tavg', '_tavgBlock')):
        found.setdefault(m.group(1), []).append(int(m.group(2)))
if not found:
    raise SystemExit(f'no "prefix.STEP" snapshots in {FIELDS}')
prefix = max(found, key=lambda k: len(found[k]))
steps = sorted(found[prefix])
steps = [s for s in steps if (args.first is None or s >= args.first) and (args.last is None or s <= args.last)][::args.stride]
print(f'prefix "{prefix}", {len(steps)} snapshots, steps {steps[0]}..{steps[-1]}')


def log_times():
    """snapshot step -> physical time: the monitor line before each 'Writing snapshot: fields/<name>.<step>' line of run.log (the
    monitor step counter restarts at 0 on a restart, the file name does not); {} if there is no log."""
    tl, t = {}, None
    try:
        for ln in open(os.path.join(CASE, 'run.log'), errors='ignore'):
            m = re.match(r'\s*\d+\s+([-+0-9.]+E[-+]\d+)\s', ln)
            if m:
                t = float(m.group(1))
                continue
            m = re.search(r'Writing snapshot:\s*\S*\.(\d+)\b', ln)
            if m and t is not None:
                tl[int(m.group(1))] = t
    except OSError:
        pass
    return tl


TLOG = log_times()


def probe(path):
    """Grid and byte offsets of every block from the header and the 3-int block headers (seeks over the data)."""
    with open(path, 'rb') as fh:
        def ri(): return struct.unpack('>i', fh.read(4))[0]
        def rd(n): return np.fromfile(fh, dtype='>f8', count=n)
        for _ in range(3):
            n = ri(); fh.seek(8*n, 1)
        xm = rd(ri()); ym = rd(ri()); zm = rd(ri())
        blocks = []
        names = ['U', 'V', 'W', 'P'] + (['C'] if SEDIMENT else []) + (['nu_t'] if SGS else []) + (['T'] if BOUSS else []) \
            + (['alpha'] if VOF else [])
        for nm in names:
            dims = struct.unpack('>3i', fh.read(12))
            blocks.append((nm, dims, fh.tell()))
            fh.seek(8*dims[0]*dims[1]*dims[2], 1)
    return xm, ym, zm, {b[0]: b[1:] for b in blocks}


xm, ym, zm, blk = probe(os.path.join(FIELDS, f'{prefix}.{steps[0]}'))
nxm, nym, nzm = len(xm), len(ym), len(zm)
for n, a in (('xm', xm), ('ym', ym), ('zm', zm)):
    a.astype('<f8').tofile(os.path.join(OUT, n + '.bin'))
print(f'grid {nxm} x {nym} x {nzm} cell centres; blocks: {", ".join(blk)}')
count = f'{nzm} {nym} {nxm}'
START = {'U': '1 1 0', 'V': '1 0 1', 'W': '0 1 1'}      # staggered components: lower face, ghost layers stripped


def slab(dims, start, seek, path, ind='        ', ds='Float'):
    return [f'{ind}<DataItem ItemType="HyperSlab" Dimensions="{count}" Type="HyperSlab">',
            f'{ind}  <DataItem Dimensions="3 3" Format="XML">{start}  1 1 1  {count}</DataItem>',
            f'{ind}  <DataItem Dimensions="{dims[2]} {dims[1]} {dims[0]}" Format="Binary" DataType="{ds}" Precision="8" '
            f'Endian="Big" Seek="{seek}">{path}</DataItem>',
            f'{ind}</DataItem>']


def scalar(name, dims, start, seek, path):
    return ([f'        <Attribute Name="{name}" Center="Node" AttributeType="Scalar">'] + slab(dims, start, seek, path, '          ')
            + ['        </Attribute>'])


rel = lambda p: os.path.relpath(p, OUT)
L = ['<?xml version="1.0" ?>', '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>', '<Xdmf Version="2.0">', '  <Domain>',
     f'    <Topology name="topo" TopologyType="3DRECTMesh" Dimensions="{count}"/>', '    <Geometry name="geo" Type="VxVyVz">']
for ax, n in (('x', nxm), ('y', nym), ('z', nzm)):
    L.append(f'      <DataItem Dimensions="{n}" Format="Binary" DataType="Float" Precision="8" Endian="Little">{ax}m.bin</DataItem>')
L += ['    </Geometry>', '    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">']

if SDF_FILE:   # (nxg, nyg, nzm): ghost cells in x and y, interior planes only in z
    sdf_dims = (nxm + 2, nym + 2, nzm)
    if os.path.getsize(SDF_FILE) != 8*sdf_dims[0]*sdf_dims[1]*sdf_dims[2]:
        print('warning: SDF file size does not match the grid; sdf skipped')
        SDF_FILE = None

def pcol(name, c0, nc, n, f):
    """Particle attribute: columns c0..c0+nc-1 of the (n, 10) block x,y,z,u,v,w,ax,ay,az,age (big endian, after int32 n and n int32 ids)."""
    ind = '          '
    return [f'{ind}<Attribute Name="{name}" Center="Node" AttributeType="{"Vector" if nc == 3 else "Scalar"}">',
            f'{ind}  <DataItem ItemType="HyperSlab" Dimensions="{n} {nc}" Type="HyperSlab">',
            f'{ind}    <DataItem Dimensions="3 2" Format="XML">0 {c0}  1 1  {n} {nc}</DataItem>',
            f'{ind}    <DataItem Dimensions="{n} 10" Format="Binary" DataType="Float" Precision="8" Endian="Big" '
            f'Seek="{4 + 4*n}">{f}</DataItem>', f'{ind}  </DataItem>', f'{ind}</Attribute>']


def particle_grid(n, f, t):
    g = ['', f'      <Grid Name="t{t:.5g}" GridType="Uniform">', f'        <Time Value="{t:.6g}"/>',
         f'        <Topology TopologyType="Polyvertex" NumberOfElements="{n}"/>', '        <Geometry GeometryType="XYZ">']
    g += pcol('x', 0, 3, n, f)[1:-1] + ['        </Geometry>']
    g += pcol('Velocity', 3, 3, n, f) + pcol('Acceleration', 6, 3, n, f) + pcol('age', 9, 1, n, f)
    g += ['          <Attribute Name="id" Center="Node" AttributeType="Scalar">',
          f'            <DataItem Dimensions="{n}" Format="Binary" DataType="Int" Precision="4" Endian="Big" Seek="4">{f}</DataItem>',
          '          </Attribute>', '      </Grid>']
    return g


TIMES, P = {}, []
for s in steps:
    path = rel(os.path.join(FIELDS, f'{prefix}.{s}'))
    t = TLOG.get(s, (sorted(found[prefix]).index(s) + 1)*TSAVE if (NSAVE < 0 and TSAVE > 0) else float(s))
    pfile = os.path.join(FIELDS, f'{prefix}_particles.{s}')
    npart = struct.unpack('>i', open(pfile, 'rb').read(4))[0] if os.path.exists(pfile) else 0
    TIMES[s] = t
    L += ['', f'      <Grid Name="t{t:.5g}">', f'        <Time Value="{t:.6g}"/>',
          '        <Topology Reference="/Xdmf/Domain/Topology[1]"/>', '        <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>',
          '        <Attribute Name="Velocity" Center="Node" AttributeType="Vector">',
          f'          <DataItem ItemType="Function" Dimensions="{count} 3" Function="JOIN($0, $1, $2)">']
    for nm in ('U', 'V', 'W'):
        d, off = blk[nm]
        L += slab(d, START[nm], off, path, '            ')
    L += ['          </DataItem>', '        </Attribute>']
    for nm, an in (('P', 'P'), ('C', 'C'), ('nu_t', 'nu_t'), ('T', 'T'), ('alpha', 'alpha')):
        if nm in blk:
            d, off = blk[nm]
            L += scalar(an, d, '1 1 1', off, path)
    if SDF_FILE:
        L += scalar('sdf', sdf_dims, '0 1 1', 0, rel(SDF_FILE))
    L.append('      </Grid>')
    if npart:
        P += particle_grid(npart, rel(pfile), t)
L += ['    </Grid>', '  </Domain>', '</Xdmf>']
out = os.path.join(OUT, f'{prefix}.xmf')
open(out, 'w').write('\n'.join(L) + '\n')
print(f'wrote {out}')
if P:
    pout = os.path.join(OUT, f'{prefix}_particles.xmf')
    head = ['<?xml version="1.0" ?>', '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>', '<Xdmf Version="2.0">', '  <Domain>',
            '    <Grid Name="ParticleTimeSeries" GridType="Collection" CollectionType="Temporal">']
    open(pout, 'w').write('\n'.join(head + P + ['    </Grid>', '  </Domain>', '</Xdmf>']) + '\n')
    print(f'wrote {pout}  ({sum("<Time" in l for l in P)} time steps)')
