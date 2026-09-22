#!/usr/bin/env python3
"""
Generate XDMF for fdm-dopamine particle snapshots (src/particles.f90) — NO DATA
DUPLICATION, same philosophy as generateXMF.py for the flow fields: the
generated .xmf points directly at the original Fortran binary snapshot files
via byte-seek HyperSlabs; no particle data is ever copied.

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

Point count varies from one snapshot to the next (particles exit/deposit/
reinject), so — unlike the fixed field-grid XDMF, which shares one Topology/
Geometry across the whole time series — every timestep here declares its own
Polyvertex Topology and XYZ Geometry, sized to that snapshot's own particle
count.

Usage:
    python3 generate_particles_xmf.py

Output (in ./paraview/ subdirectory, alongside generateXMF.py's own output):
    particles.xmf   – XDMF time series; open this in ParaView (or combine
                       with channel_test.xmf's time toolbar to scrub fields
                       and particles together -- Time Value is the step
                       number, matching generateXMF.py's own DT=1 convention
                       exactly, not physical time, so the two time axes agree)

Requires the id/age fields shown as point-cloud Attributes: colour by
"id" to track individual particles, or by "age"/"Velocity" (magnitude or
component) for dispersion/turbophoresis-style plots. A Glyph filter
(Sphere, small radius) on the point cloud is usually clearer than the raw
Points representation in ParaView.
"""

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


def probe_count(fpath):
    """Read just the leading Int32 particle count from one snapshot file."""
    with open(fpath, 'rb') as fh:
        return struct.unpack('>i', fh.read(4))[0]


def rel(path):
    return os.path.relpath(path, OUT_DIR)


def xdmf_grid_for_step(total, frel, t):
    """One <Grid> (Polyvertex point cloud + attributes) for a single snapshot."""
    id_bytes = 4          # Int32
    dat_offset = 4 + 4 * total   # past the Int32 total header and the Int32 id[] array

    lines = [
        '',
        f'      <Grid Name="t{t:.6g}" GridType="Uniform">',
        f'        <Time Value="{t:.6g}"/>',
        f'        <Topology TopologyType="Polyvertex" NumberOfElements="{total}"/>',
        '        <Geometry GeometryType="XYZ">',
        f'          <DataItem ItemType="HyperSlab" Dimensions="{total} 3" Type="HyperSlab">',
        '            <DataItem Dimensions="3 2" Format="XML">',
        '              0 0',
        '              1 1',
        f'              {total} 3',
        '            </DataItem>',
        f'            <DataItem Dimensions="{total} 7" Format="Binary"',
        f'                     DataType="Float" Precision="8" Endian="Big" Seek="{dat_offset}">',
        f'              {frel}',
        '            </DataItem>',
        '          </DataItem>',
        '        </Geometry>',
        '',
        '        <Attribute Name="Velocity" Center="Node" AttributeType="Vector">',
        f'          <DataItem ItemType="HyperSlab" Dimensions="{total} 3" Type="HyperSlab">',
        '            <DataItem Dimensions="3 2" Format="XML">',
        '              0 3',
        '              1 1',
        f'              {total} 3',
        '            </DataItem>',
        f'            <DataItem Dimensions="{total} 7" Format="Binary"',
        f'                     DataType="Float" Precision="8" Endian="Big" Seek="{dat_offset}">',
        f'              {frel}',
        '            </DataItem>',
        '          </DataItem>',
        '        </Attribute>',
        '',
        '        <Attribute Name="age" Center="Node" AttributeType="Scalar">',
        f'          <DataItem ItemType="HyperSlab" Dimensions="{total} 1" Type="HyperSlab">',
        '            <DataItem Dimensions="3 2" Format="XML">',
        '              0 6',
        '              1 1',
        f'              {total} 1',
        '            </DataItem>',
        f'            <DataItem Dimensions="{total} 7" Format="Binary"',
        f'                     DataType="Float" Precision="8" Endian="Big" Seek="{dat_offset}">',
        f'              {frel}',
        '            </DataItem>',
        '          </DataItem>',
        '        </Attribute>',
        '',
        '        <Attribute Name="id" Center="Node" AttributeType="Scalar">',
        f'          <DataItem Dimensions="{total}" Format="Binary"',
        f'                   DataType="Int" Precision="4" Endian="Big" Seek="{id_bytes}">',
        f'            {frel}',
        '          </DataItem>',
        '        </Attribute>',
        '      </Grid>',
    ]
    return lines


def main():
    prefix, steps = detect_particle_prefix_and_steps()
    # Time Value = step (not physical time): matches generateXMF.py's own DT=1
    # convention exactly, which is what lets ParaView's shared time toolbar line
    # fields and particles up frame-for-frame when both .xmf files are open
    # together. (generateXMF.py doesn't read dt from input_parameters either --
    # see its own DT=1 comment -- so scaling by dt here would only make the two
    # time axes disagree.)

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
        total = probe_count(fpath)
        if total == 0:
            # An empty Polyvertex grid confuses some ParaView versions; skip
            # steps with no active particles rather than emitting one.
            print(f'  step {step}: 0 particles, skipping')
            continue
        t = float(step)
        lines += xdmf_grid_for_step(total, rel(fpath), t)
        n_written += 1
        print(f'  registered step {step} ({total} particles)')

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


if __name__ == '__main__':
    main()
