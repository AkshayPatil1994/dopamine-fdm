#!/usr/bin/env python3
"""
Generate XDMF for fdm-dopamine binary snapshots — NO DATA DUPLICATION.

The generated channel_test.xmf points directly to the original Fortran
binary field files via byte-seek offsets and XDMF HyperSlab selections.
Only the small 1-D coordinate arrays (xm.bin, ym.bin, zm.bin) are written
to the paraview/ subdirectory; the field data is never copied.

Binary format: big-endian, Fortran stream-access (no record markers).
Fields are on a staggered grid; ghost-layer stripping and the selection of
the lower face value per cell are handled entirely within the XDMF HyperSlab
— no averaging and no temporary field files are created.

Usage:
    python3 generateXMF.py

Output (in ./paraview/ subdirectory):
    xm.bin, ym.bin, zm.bin   – cell-centre coordinate arrays (float64 LE)
    channel_test.xmf         – XDMF time-series; open this in ParaView
"""

import struct
import numpy as np
import os

# ── paths ──────────────────────────────────────────────────────────────────────
CASE_DIR   = os.path.dirname(os.path.abspath(__file__))
FIELDS_DIR = os.path.join(CASE_DIR, 'fields')
OUT_DIR    = os.path.join(CASE_DIR, 'paraview')
os.makedirs(OUT_DIR, exist_ok=True)

# ── auto-detect snapshot prefix and steps ─────────────────────────────────
def detect_field_prefix_and_steps():
    """Scan FIELDS_DIR for '<prefix>.<step>' files; return (prefix, sorted step array)."""
    if not os.path.isdir(FIELDS_DIR):
        raise SystemExit(f'Fields directory not found: {FIELDS_DIR}')

    files = os.listdir(FIELDS_DIR)

    import re
    pattern = re.compile(r'^([a-zA-Z0-9_]+)\.(\d+)$')

    detected_prefix = None
    steps_list = []

    for fname in files:
        match = pattern.match(fname)
        if match:
            prefix, step_str = match.groups()
            if detected_prefix is None:
                detected_prefix = prefix
            elif detected_prefix != prefix:
                print(f'Warning: Multiple prefixes found ({detected_prefix}, {prefix}). '
                      f'Using {detected_prefix}.')
            if prefix == detected_prefix:
                steps_list.append(int(step_str))
    
    if not steps_list:
        raise SystemExit(f'No field files matching pattern "prefix.STEP" found in {FIELDS_DIR}')
    
    steps = np.array(sorted(steps_list))
    print(f'Auto-detected prefix: "{detected_prefix}", found {len(steps)} files')
    print(f'  Steps: {steps[0]} to {steps[-1]} (count: {len(steps)})')
    
    return detected_prefix, steps


# ── snapshot list ─────────────────────────────────────────────────────────
field_name_prefix, STEPS = detect_field_prefix_and_steps()
DT = 1

# ── optional fields ────────────────────────────────────────────────────────────
# Must match the run that produced the files: blocks are appended C then nu_t, after P.
HAS_SEDIMENT = False   # scalar C block after P (sediment_flag >= 1)
HAS_SGS      = True    # nu_t block after P/C (sgs_model != 0)

# ── header probe ──────────────────────────────────────────────────────────────
def probe_snapshot(fpath):
    """
    Read only the binary header of one snapshot file.

    Returns a dict with:
      nxm, nym, nzm        – number of cell-centre points per axis
      xm, ym, zm           – 1-D coordinate arrays
      un, vn, wn, pn       – (nx, ny, nz) raw shapes of U, V, W, P
      off_U, off_V,
      off_W, off_P         – byte offsets to the start of each field's
                             float64 data within the file
    """
    with open(fpath, 'rb') as fh:
        data = fh.read()
    pos = [0]

    def ri():
        v = struct.unpack_from('>i', data, pos[0])[0]
        pos[0] += 4
        return v

    def skip8(n):
        pos[0] += n * 8

    def rd(n):
        arr = np.frombuffer(data, dtype='>f8', count=n, offset=pos[0]).copy()
        pos[0] += n * 8
        return arr

    def ri3():
        v = struct.unpack_from('>3i', data, pos[0])
        pos[0] += 12
        return v

    # ── grid header ───────────────────────────────────────────────────────────
    nx  = ri();  skip8(nx)       # x face points (discard)
    ny  = ri();  skip8(ny)       # y face points (discard)
    nz  = ri();  skip8(nz)       # z face points (discard)
    nxm = ri();  xm = rd(nxm)   # x cell centres
    nym = ri();  ym = rd(nym)   # y cell centres
    nzm = ri();  zm = rd(nzm)   # z cell centres

    # ── field shapes and byte offsets ─────────────────────────────────────────
    un = ri3();  off_U = pos[0];  skip8(un[0] * un[1] * un[2])
    vn = ri3();  off_V = pos[0];  skip8(vn[0] * vn[1] * vn[2])
    wn = ri3();  off_W = pos[0];  skip8(wn[0] * wn[1] * wn[2])
    pn = ri3();  off_P = pos[0];  skip8(pn[0] * pn[1] * pn[2])

    off_nut = None
    nutn    = None
    cn      = None
    off_C   = None
    if HAS_SEDIMENT:
        cn = ri3();  off_C = pos[0];  skip8(cn[0] * cn[1] * cn[2])
    if HAS_SGS:
        nutn = ri3();  off_nut = pos[0]

    return dict(nxm=nxm, nym=nym, nzm=nzm,
                xm=xm, ym=ym, zm=zm,
                un=un, vn=vn, wn=wn, pn=pn,
                off_U=off_U, off_V=off_V, off_W=off_W, off_P=off_P,
                cn=cn, off_C=off_C,
                nutn=nutn, off_nut=off_nut)


# ── coordinate writer ─────────────────────────────────────────────────────────
def write_1d(arr, fpath):
    """Write a 1-D coordinate array as little-endian float64."""
    arr.astype('<f8').tofile(fpath)


# ── main ──────────────────────────────────────────────────────────────────────
# Probe the first available snapshot to learn grid dimensions and byte offsets.
# No field data is read beyond the header.
print('Probing snapshots …')

info = None
for step in STEPS:
    fpath = os.path.join(FIELDS_DIR, f'{field_name_prefix}.{step}')
    if not os.path.exists(fpath):
        print(f'  step {step}: file not found, skipping')
        continue
    if info is None:
        info = probe_snapshot(fpath)
        nxm, nym, nzm = info['nxm'], info['nym'], info['nzm']
        write_1d(info['xm'], os.path.join(OUT_DIR, 'xm.bin'))
        write_1d(info['ym'], os.path.join(OUT_DIR, 'ym.bin'))
        write_1d(info['zm'], os.path.join(OUT_DIR, 'zm.bin'))
        print(f'  grid: {nxm} x {nym} x {nzm} cell-centre nodes')
    print(f'  registered step {step}')

if info is None:
    raise SystemExit('No snapshot files found — check FIELDS_DIR and field_name_prefix.')

un   = info['un'];   vn   = info['vn']
wn   = info['wn'];   pn   = info['pn']
offU = info['off_U']; offV = info['off_V']
offW = info['off_W']; offP = info['off_P']
cn     = info['cn']
offC   = info['off_C']
nutn   = info['nutn']
offNut = info['off_nut']

# ── write XDMF ────────────────────────────────────────────────────────────────
# Binary files are big-endian Fortran column-major (x varies fastest).
# XDMF Dimensions are slowest→fastest, so a Fortran (nx,ny,nz) array is
# declared as  Dimensions="nz ny nx"  to preserve the correct axis mapping.
#
# HyperSlab selects the cell-centre interior, stripping ghost layers and
# picking the *lower* face for each staggered velocity component.
# No averaging is performed — acceptable for visualisation.
#
# Ghost-layer layout (0-based indices in each axis):
#   U (staggered x): shape (nx,  nyg, nzg)
#     → XDMF full Dims="nzg nyg nx",   Start="1 1 0", Count="nzm nym nxm"
#   V (staggered y): shape (nxg, ny,  nzg)
#     → XDMF full Dims="nzg ny  nxg",  Start="1 0 1", Count="nzm nym nxm"
#   W (staggered z): shape (nxg, nyg, nz)
#     → XDMF full Dims="nz  nyg nxg",  Start="0 1 1", Count="nzm nym nxm"
#   P (cell-centre): shape (nxg, nyg, nzg)
#     → XDMF full Dims="nzg nyg nxg",  Start="1 1 1", Count="nzm nym nxm"

def rel(path):
    """Relative path from OUT_DIR to the given absolute path."""
    return os.path.relpath(path, OUT_DIR)


def xdmf_hyperslab(attr_name, full_dims_str, start_str, seek, fpath_rel):
    """Return XDMF lines for one scalar Attribute using a HyperSlab."""
    count = f'{nzm} {nym} {nxm}'
    return [
        f'        <Attribute Name="{attr_name}" Center="Node" AttributeType="Scalar">',
        f'          <DataItem ItemType="HyperSlab" Dimensions="{count}" Type="HyperSlab">',
        f'            <DataItem Dimensions="3 3" Format="XML">',
        f'              {start_str}',
        f'              1 1 1',
        f'              {count}',
        f'            </DataItem>',
        f'            <DataItem Dimensions="{full_dims_str}" Format="Binary"',
        f'                     DataType="Float" Precision="8" Endian="Big" Seek="{seek}">',
        f'              {fpath_rel}',
        f'            </DataItem>',
        f'          </DataItem>',
        f'        </Attribute>',
    ]


def _hyperslab_item(full_dims_str, start_str, seek, fpath_rel):
    """Bare HyperSlab DataItem (no Attribute wrapper) — used inside JOIN."""
    count = f'{nzm} {nym} {nxm}'
    return [
        f'            <DataItem ItemType="HyperSlab" Dimensions="{count}" Type="HyperSlab">',
        f'              <DataItem Dimensions="3 3" Format="XML">',
        f'                {start_str}',
        f'                1 1 1',
        f'                {count}',
        f'              </DataItem>',
        f'              <DataItem Dimensions="{full_dims_str}" Format="Binary"',
        f'                       DataType="Float" Precision="8" Endian="Big" Seek="{seek}">',
        f'                {fpath_rel}',
        f'              </DataItem>',
        f'            </DataItem>',
    ]


def xdmf_velocity(un, vn, wn, offU, offV, offW, fpath_rel):
    """Return XDMF lines for a 3-component Velocity vector via XDMF JOIN."""
    count = f'{nzm} {nym} {nxm}'
    ux, uy, uz = un
    vx, vy, vz = vn
    wx, wy, wz = wn
    lines = [
        f'        <Attribute Name="Velocity" Center="Node" AttributeType="Vector">',
        f'          <DataItem ItemType="Function" Dimensions="{count} 3"',
        f'                   Function="JOIN($0, $1, $2)">',
    ]
    lines += _hyperslab_item(f'{uz} {uy} {ux}', '1 1 0', offU, fpath_rel)
    lines += _hyperslab_item(f'{vz} {vy} {vx}', '1 0 1', offV, fpath_rel)
    lines += _hyperslab_item(f'{wz} {wy} {wx}', '0 1 1', offW, fpath_rel)
    lines += [
        f'          </DataItem>',
        f'        </Attribute>',
    ]
    return lines


ux, uy, uz = un   # e.g. 130, 66, 67  (U: no ghost in x, ghost in y and z)
vx, vy, vz = vn   # e.g. 131, 65, 67  (V: ghost in x and z, staggered in y)
wx, wy, wz = wn   # e.g. 131, 66, 66  (W: ghost in x and y, staggered in z)
px, py, pz = pn   # e.g. 131, 66, 67  (P: ghost in x, y and z)

lines = ['<?xml version="1.0" ?>',
         '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>',
         '<Xdmf Version="2.0">',
         '  <Domain>',
         '',
         f'    <!-- Rectilinear mesh: {nxm} x {nym} x {nzm} cell-centre nodes -->',
         f'    <Topology name="topo" TopologyType="3DRECTMesh"',
         f'              Dimensions="{nzm} {nym} {nxm}"/>',
         '',
         '    <Geometry name="geo" Type="VxVyVz">',
         f'      <!-- Vx → x-axis ({nxm} values), Vy → y-axis, Vz → z-axis -->',
         f'      <DataItem Name="Vx" Dimensions="{nxm}" Format="Binary"',
         '               DataType="Float" Precision="8" Endian="Little">',
         '        xm.bin',
         '      </DataItem>',
         f'      <DataItem Name="Vy" Dimensions="{nym}" Format="Binary"',
         '               DataType="Float" Precision="8" Endian="Little">',
         '        ym.bin',
         '      </DataItem>',
         f'      <DataItem Name="Vz" Dimensions="{nzm}" Format="Binary"',
         '               DataType="Float" Precision="8" Endian="Little">',
         '        zm.bin',
         '      </DataItem>',
         '    </Geometry>',
         '',
         '    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">',
         ]

for step in STEPS:
    fpath = os.path.join(FIELDS_DIR, f'{field_name_prefix}.{step}')
    if not os.path.exists(fpath):
        continue
    frel = rel(fpath)
    t    = step * DT
    lines += [
        '',
        f'      <Grid Name="t{t:.4g}">',
        f'        <Time Value="{t:.6g}"/>',
        '        <Topology Reference="/Xdmf/Domain/Topology[1]"/>',
        '        <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>',
    ]
    lines += xdmf_velocity(un, vn, wn, offU, offV, offW, frel)
    lines += xdmf_hyperslab('P', f'{pz} {py} {px}', '1 1 1', offP, frel)
    if HAS_SEDIMENT and offC is not None:
        cx, cy, cz = cn   # same shape as P: (nxg, nyg, nzg)
        lines += xdmf_hyperslab('C', f'{cz} {cy} {cx}', '1 1 1', offC, frel)
    if HAS_SGS and offNut is not None:
        ntx, nty, ntz = nutn   # same shape as P: (nxg, nyg, nzg)
        lines += xdmf_hyperslab('nu_t', f'{ntz} {nty} {ntx}', '1 1 1', offNut, frel)
    lines.append('      </Grid>')

lines += [
    '    </Grid>',
    '',
    '  </Domain>',
    '</Xdmf>',
]

xmf_path = os.path.join(OUT_DIR, 'channel_test.xmf')
with open(xmf_path, 'w') as f:
    f.write('\n'.join(lines) + '\n')

print(f'\nXDMF written → {xmf_path}')
print('Open channel_test.xmf in ParaView (File → Open).')
