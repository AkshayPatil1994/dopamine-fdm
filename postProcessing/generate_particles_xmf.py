#!/usr/bin/env python3
"""
Generate XDMF for fdm-dopamine particle snapshots (src/particles.f90).

Binary format written by particles.f90's write_particle_snapshot (also used,
same layout, by the hot-restart file from write_particle_restart):

    Int32              total            -- particle count this snapshot
    Int32[total]       id               -- global particle ID
    Float64[total, 7]  x,y,z,u,v,w,age  -- row-major (C order): one row per
                                           particle, matching the on-disk
                                           layout exactly (no reshape needed)

All integers/floats are big-endian (this build's global -fconvert=big-endian/
-convert big_endian/-Mbyteswapio Fortran flag, same convention generateXMF.py
already relies on for the field snapshots).

Particle count varies from one snapshot to the next (particles exit/deposit/
reinject). By default this script PADS every snapshot up to the run's own
maximum particle count and writes one small auxiliary binary file per
timestep into ./paraview/ (id=-1, active=0 for padding rows) -- unlike
generateXMF.py's fields, which are large enough that zero-copy byte-seeking
into the original files matters, particle snapshots are tiny, and giving
every timestep the SAME Topology/Geometry size is required for vtkXdmfReader
to treat the series as one homogeneous time-varying dataset rather than
promoting it to a vtkMultiBlockDataSet (a real ParaView reader limitation --
without this, ParaView prints "Data type generated (vtkMultiBlockDataSet)
does not match data type expected (vtkUnstructuredGrid)" and animation can
misbehave). Pass --no-pad to fall back to the original zero-duplication,
variable-size-per-timestep behaviour if you don't hit that warning.

Usage:
    python3 generate_particles_xmf.py            # padded (default, fixes the ParaView warning)
    python3 generate_particles_xmf.py --no-pad    # original byte-seek-only, no auxiliary files

Output (in ./paraview/ subdirectory, alongside generateXMF.py's own output):
    particles.xmf                        – XDMF time series; open in ParaView
    <prefix>_particles_padded.<step>.bin – padded per-timestep data (--pad only)

Open particles.xmf alongside channel_test.xmf (generateXMF.py's output) and use
the shared time toolbar to animate fields and particles together -- Time Value
is the step number, matching generateXMF.py's own DT=1 convention exactly, not
physical time, so the two time axes agree.

Point-cloud Attributes: "id" (track individual particles, -1 on padding rows
with --pad), "age", "Velocity", and (--pad only) "active" (1=real particle,
0=padding -- Threshold on active>0.5 to hide the padding rows). A Glyph filter
(Sphere, small radius) on the point cloud is usually clearer than the raw
Points representation in ParaView.
"""

import argparse
import os
import re
import struct

import numpy as np

CASE_DIR = os.path.dirname(os.path.abspath(__file__))
FIELDS_DIR = os.path.join(CASE_DIR, 'fields')
OUT_DIR = os.path.join(CASE_DIR, 'paraview')
os.makedirs(OUT_DIR, exist_ok=True)


def detect_particle_prefix_and_steps():
    """Scan FIELDS_DIR for '<prefix>_particles.STEP' files; return (prefix, sorted steps)."""
    if not os.path.isdir(FIELDS_DIR):
        raise SystemExit(f'Fields directory not found: {FIELDS_DIR}')

    pattern = re.compile(r'^([a-zA-Z0-9_]+)_particles\.(\d+)$')
    detected_prefix = None
    steps_list = []

    for fname in os.listdir(FIELDS_DIR):
        match = pattern.match(fname)
        if match:
            prefix, step_str = match.groups()
            if detected_prefix is None:
                detected_prefix = prefix
            steps_list.append(int(step_str))

    if not steps_list:
        raise SystemExit(
            f'No particle snapshot files matching "*_particles.STEP" found in {FIELDS_DIR} '
            '-- confirm particles_active=1 in &PARTICLES and that the run has saved at '
            'least one snapshot (nsave in &NUMERICS).')

    steps = np.array(sorted(steps_list))
    print(f'Auto-detected prefix: "{detected_prefix}", found {len(steps)} particle snapshots')
    print(f'  Steps: {steps[0]} to {steps[-1]} (count: {len(steps)})')
    return detected_prefix, steps


def read_snapshot(fpath):
    """Read one full snapshot: (total, id[total] int32, dat[total,7] float64)."""
    with open(fpath, 'rb') as fh:
        data = fh.read()
    total = struct.unpack_from('>i', data, 0)[0]
    if total == 0:
        return 0, np.zeros(0, dtype='>i4'), np.zeros((0, 7), dtype='>f8')
    ids = np.frombuffer(data, dtype='>i4', count=total, offset=4)
    dat = np.frombuffer(data, dtype='>f8', count=total * 7, offset=4 + 4 * total).reshape(total, 7)
    return total, ids, dat


def rel(path):
    return os.path.relpath(path, OUT_DIR)


def write_padded(fpath, ids, dat, max_total):
    """Write one padded auxiliary file: Int32[max_total] id, Int32[max_total] active,
    Float64[max_total,7] dat -- padding rows get id=-1, active=0, dat=0."""
    total = len(ids)
    id_pad = np.full(max_total, -1, dtype='>i4')
    id_pad[:total] = ids
    active_pad = np.zeros(max_total, dtype='>i4')
    active_pad[:total] = 1
    dat_pad = np.zeros((max_total, 7), dtype='>f8')
    dat_pad[:total, :] = dat

    with open(fpath, 'wb') as fh:
        fh.write(id_pad.tobytes())
        fh.write(active_pad.tobytes())
        fh.write(dat_pad.tobytes())


def xdmf_grid_for_step(n, frel, t, padded):
    """One <Grid> (Polyvertex point cloud + attributes) for a single (already-sized) timestep.
    n is the Topology/Geometry element count: the padded ceiling if padded=True, else this
    snapshot's own (possibly different from every other step's) real particle count."""
    if padded:
        id_off = 0
        active_off = 4 * n
        dat_off = 8 * n   # past id[n] (Int32) and active[n] (Int32)
    else:
        id_off = 4          # past the Int32 total header
        dat_off = 4 + 4 * n   # past the Int32 total header and the Int32 id[] array

    lines = [
        '',
        f'      <Grid Name="t{t:.6g}" GridType="Uniform">',
        f'        <Time Value="{t:.6g}"/>',
        f'        <Topology TopologyType="Polyvertex" NumberOfElements="{n}"/>',
        '        <Geometry GeometryType="XYZ">',
        f'          <DataItem ItemType="HyperSlab" Dimensions="{n} 3" Type="HyperSlab">',
        '            <DataItem Dimensions="3 2" Format="XML">',
        '              0 0',
        '              1 1',
        f'              {n} 3',
        '            </DataItem>',
        f'            <DataItem Dimensions="{n} 7" Format="Binary"',
        f'                     DataType="Float" Precision="8" Endian="Big" Seek="{dat_off}">',
        f'              {frel}',
        '            </DataItem>',
        '          </DataItem>',
        '        </Geometry>',
        '',
        '        <Attribute Name="Velocity" Center="Node" AttributeType="Vector">',
        f'          <DataItem ItemType="HyperSlab" Dimensions="{n} 3" Type="HyperSlab">',
        '            <DataItem Dimensions="3 2" Format="XML">',
        '              0 3',
        '              1 1',
        f'              {n} 3',
        '            </DataItem>',
        f'            <DataItem Dimensions="{n} 7" Format="Binary"',
        f'                     DataType="Float" Precision="8" Endian="Big" Seek="{dat_off}">',
        f'              {frel}',
        '            </DataItem>',
        '          </DataItem>',
        '        </Attribute>',
        '',
        '        <Attribute Name="age" Center="Node" AttributeType="Scalar">',
        f'          <DataItem ItemType="HyperSlab" Dimensions="{n} 1" Type="HyperSlab">',
        '            <DataItem Dimensions="3 2" Format="XML">',
        '              0 6',
        '              1 1',
        f'              {n} 1',
        '            </DataItem>',
        f'            <DataItem Dimensions="{n} 7" Format="Binary"',
        f'                     DataType="Float" Precision="8" Endian="Big" Seek="{dat_off}">',
        f'              {frel}',
        '            </DataItem>',
        '          </DataItem>',
        '        </Attribute>',
        '',
        '        <Attribute Name="id" Center="Node" AttributeType="Scalar">',
        f'          <DataItem Dimensions="{n}" Format="Binary"',
        f'                   DataType="Int" Precision="4" Endian="Big" Seek="{id_off}">',
        f'            {frel}',
        '          </DataItem>',
        '        </Attribute>',
    ]
    if padded:
        lines += [
            '',
            '        <Attribute Name="active" Center="Node" AttributeType="Scalar">',
            f'          <DataItem Dimensions="{n}" Format="Binary"',
            f'                   DataType="Int" Precision="4" Endian="Big" Seek="{active_off}">',
            f'            {frel}',
            '          </DataItem>',
            '        </Attribute>',
        ]
    lines.append('      </Grid>')
    return lines


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--no-pad', dest='pad', action='store_false',
                     help='disable padding: original zero-duplication, variable-size-per-timestep '
                          'behaviour (may trigger a vtkXdmfReader vtkMultiBlockDataSet warning in ParaView)')
    args = ap.parse_args()

    prefix, steps = detect_particle_prefix_and_steps()

    max_total = 0
    if args.pad:
        for step in steps:
            fpath = os.path.join(FIELDS_DIR, f'{prefix}_particles.{step}')
            if os.path.exists(fpath):
                max_total = max(max_total, struct.unpack('>i', open(fpath, 'rb').read(4))[0])
        print(f'Padding every snapshot up to {max_total} particles (this run\'s maximum)')

    lines = [
        '<?xml version="1.0" ?>',
        '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>',
        '<Xdmf Version="2.0">',
        '  <Domain>',
        '    <Grid Name="ParticleTimeSeries" GridType="Collection" CollectionType="Temporal">',
    ]

    n_written = 0
    for step in steps:
        fpath = os.path.join(FIELDS_DIR, f'{prefix}_particles.{step}')
        if not os.path.exists(fpath):
            continue
        t = float(step)

        if args.pad:
            total, ids, dat = read_snapshot(fpath)
            if max_total == 0:
                print(f'  step {step}: 0 particles in every snapshot, skipping')
                continue
            padded_path = os.path.join(OUT_DIR, f'{prefix}_particles_padded.{step}.bin')
            write_padded(padded_path, ids, dat, max_total)
            lines += xdmf_grid_for_step(max_total, rel(padded_path), t, padded=True)
            print(f'  registered step {step} ({total}/{max_total} active particles)')
        else:
            total = probe_count_only(fpath)
            if total == 0:
                print(f'  step {step}: 0 particles, skipping')
                continue
            lines += xdmf_grid_for_step(total, rel(fpath), t, padded=False)
            print(f'  registered step {step} ({total} particles)')
        n_written += 1

    lines += [
        '    </Grid>',
        '  </Domain>',
        '</Xdmf>',
    ]

    if n_written == 0:
        raise SystemExit('No non-empty particle snapshots found -- nothing to write.')

    xmf_path = os.path.join(OUT_DIR, 'particles.xmf')
    with open(xmf_path, 'w') as f:
        f.write('\n'.join(lines) + '\n')

    print(f'\nXDMF written -> {xmf_path}')
    print('Open particles.xmf in ParaView (File -> Open), or open it alongside '
          "channel_test.xmf (generateXMF.py's output) and use the shared time "
          'toolbar to animate fields and particles together.')
    if args.pad:
        print('Padding rows have active=0 -- apply Threshold (active > 0.5) to hide them.')


def probe_count_only(fpath):
    with open(fpath, 'rb') as fh:
        return struct.unpack('>i', fh.read(4))[0]


if __name__ == '__main__':
    main()
