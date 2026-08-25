#!/usr/bin/env python3
"""
NOTICE
    -- LLM GENERATED CODE --

generate_wavy_wall.py — create a sinusoidal wavy-bottom-wall channel solid
as a closed, watertight binary STL.

Replicates the ERCOFTAC Classic Collection Case 076 "Flow over a wavy wall"
geometry (Buckles, Hanratty & Adrian, JFM 1984; Re_H = U_b*H/nu = 6760,
LDV measurements over a sinusoidal Plexiglas wavy wall, flat opposite wall):

    (y - y0) / a = cos(2*pi*x/lambda + phase)     a = amplitude, y0 = mean wave height
    2a/lambda = 0.10                              (height-to-wavelength ratio, case076 spec)

Only the bottom wall is wavy; the top wall stays the solver's flat no-slip
BC (bc_face_yhi=1) — this script only builds the IBM solid for the bottom.

Coordinate convention (matching the dopamine solver): x streamwise,
y wall-normal (vertical), z spanwise. The wave surface is placed so its
TROUGH sits exactly at y=0 (the domain floor) and its CREST at y=2a, i.e.
    y0 = a
    wave_y(x) = a + a*cos(2*pi*x/lambda + phase)
so the entire wavy surface lies within the domain (no wave excursion
below y=0) and the solid fills everything below it, down to a safety
margin below y=0 so the solid has nonzero thickness even at the trough.

The streamwise domain length Lx MUST equal an integer number of
wavelengths (N_WAVES) for the wave to close exactly under the solver's
periodic x boundary condition (x_bc_type=0) — wave_y(0) == wave_y(Lx)
by construction here, so the geometry is periodicity-safe regardless of
mesh resolution.

Output: wavy_wall.stl (this directory)

Usage (run from this directory):
    python3 generate_wavy_wall.py

The generated STL is consumed by GenSDF (gensdf_mpi, see
preProcessing/GenSDF) to compute the signed-distance field on the
solver's Cartesian grid; the resulting cell-centre SDF must be renamed
'SDF_in' in the run directory with ibm_input_mode=1 in input_parameters.
"""

import os
import struct
import numpy as np

# ---------------------------------------------------------------------------
# Domain — must match this case's input_parameters &DOMAIN block
# ---------------------------------------------------------------------------
LX, LY, LZ = 4.0, 1.05, 2.0   # m

# ---------------------------------------------------------------------------
# Wave parameters (ERCOFTAC case076 / Buckles-Hanratty-Adrian 1984)
# ---------------------------------------------------------------------------
N_WAVES              = 4      # integer number of full wavelengths across Lx (periodicity)
HEIGHT_TO_WAVELENGTH = 0.10   # 2*a/lambda, case076 spec
PHASE                = 0.0    # rad, phase offset (coordinate-origin dependent; case076's
                               # published profile uses 0.144 rad relative to their own
                               # origin — irrelevant here since our origin is arbitrary)

WAVELENGTH = LX / N_WAVES
AMPLITUDE  = 0.5 * HEIGHT_TO_WAVELENGTH * WAVELENGTH   # a = 0.05*lambda
Y_WAVE_MEAN = AMPLITUDE                                 # trough at y=0, crest at y=2a

# ---------------------------------------------------------------------------
# Solid-block / mesh parameters
# ---------------------------------------------------------------------------
Y_BOTTOM     = -0.25 * LY   # flat base of the solid block, safely below y=0
Z_MARGIN     = 0.05 * LZ    # extend end caps beyond [0,Lz] to avoid SDF edge ambiguity
N_X          = 800          # streamwise mesh resolution (points along the wave)
OUT_FILE     = 'wavy_wall.stl'
# ---------------------------------------------------------------------------


def wave_y(x: np.ndarray) -> np.ndarray:
    """Wavy-wall surface height y(x)."""
    return Y_WAVE_MEAN + AMPLITUDE * np.cos(2.0 * np.pi * x / WAVELENGTH + PHASE)


def _tri(v0, v1, v2, ref):
    """Triangle (v0,v1,v2) with its normal flipped to align with reference outward direction."""
    n = np.cross(v1 - v0, v2 - v0)
    if np.dot(n, ref) < 0.0:
        v1, v2 = v2, v1
        n = -n
    norm = np.linalg.norm(n)
    return (n / max(norm, 1.0e-30), v0, v1, v2)


def build_triangles():
    """Build the closed wavy-bottom-wall solid: top (wave), bottom, left, right, and two end caps."""
    x  = np.linspace(0.0, LX, N_X + 1)
    yw = wave_y(x)
    z0, z1 = -Z_MARGIN, LZ + Z_MARGIN

    tris = []

    # --- Top surface (the wavy wall itself): ruled surface between z0 and z1 ---
    for i in range(N_X):
        A = np.array([x[i],   yw[i],   z0])
        B = np.array([x[i+1], yw[i+1], z0])
        C = np.array([x[i+1], yw[i+1], z1])
        D = np.array([x[i],   yw[i],   z1])
        ref = np.array([0.0, 1.0, 0.0])
        tris.append(_tri(A, B, C, ref))
        tris.append(_tri(A, C, D, ref))

    # --- Bottom surface (flat base, y=Y_BOTTOM) ---
    for i in range(N_X):
        A = np.array([x[i],   Y_BOTTOM, z0])
        B = np.array([x[i+1], Y_BOTTOM, z0])
        C = np.array([x[i+1], Y_BOTTOM, z1])
        D = np.array([x[i],   Y_BOTTOM, z1])
        ref = np.array([0.0, -1.0, 0.0])
        tris.append(_tri(A, B, C, ref))
        tris.append(_tri(A, C, D, ref))

    # --- Left side (x=0 plane), single planar quad ---
    A = np.array([0.0, Y_BOTTOM, z0]); B = np.array([0.0, Y_BOTTOM, z1])
    C = np.array([0.0, yw[0],    z1]); D = np.array([0.0, yw[0],    z0])
    ref = np.array([-1.0, 0.0, 0.0])
    tris.append(_tri(A, B, C, ref))
    tris.append(_tri(A, C, D, ref))

    # --- Right side (x=Lx plane), single planar quad ---
    A = np.array([LX, Y_BOTTOM, z0]); B = np.array([LX, Y_BOTTOM, z1])
    C = np.array([LX, yw[-1],   z1]); D = np.array([LX, yw[-1],   z0])
    ref = np.array([1.0, 0.0, 0.0])
    tris.append(_tri(A, B, C, ref))
    tris.append(_tri(A, C, D, ref))

    # --- End caps (z=z0 and z=z1): vertical strip triangulation ---
    for zc, ref in ((z0, np.array([0.0, 0.0, -1.0])), (z1, np.array([0.0, 0.0, 1.0]))):
        for i in range(N_X):
            A = np.array([x[i],   Y_BOTTOM, zc])
            B = np.array([x[i+1], Y_BOTTOM, zc])
            C = np.array([x[i+1], yw[i+1],  zc])
            D = np.array([x[i],   yw[i],    zc])
            tris.append(_tri(A, B, C, ref))
            tris.append(_tri(A, C, D, ref))

    return tris


def write_stl_binary(tris: list, path: str) -> None:
    """Write triangles to a binary STL file (little-endian)."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)

    hdr_str = b'Binary STL wavy wall - generate_wavy_wall.py (dopamine IBM, ERCOFTAC case076)'
    hdr_str = hdr_str[:80].ljust(80, b'\x00')

    with open(path, 'wb') as f:
        f.write(hdr_str)
        f.write(struct.pack('<I', len(tris)))
        for n, v0, v1, v2 in tris:
            data = np.array([*n, *v0, *v1, *v2], dtype=np.float32)
            f.write(data.tobytes())
            f.write(struct.pack('<H', 0))   # attribute byte count = 0


if __name__ == '__main__':
    print("Generating wavy-wall solid STL (ERCOFTAC case076-style):")
    print(f"  Domain (Lx,Ly,Lz)     = ({LX}, {LY}, {LZ}) m")
    print(f"  N_WAVES               = {N_WAVES}")
    print(f"  Wavelength lambda     = {WAVELENGTH:.6f} m")
    print(f"  Amplitude a           = {AMPLITUDE:.6f} m   (2a = {2*AMPLITUDE:.6f} m)")
    print(f"  2a/lambda             = {2*AMPLITUDE/WAVELENGTH:.4f}  (target 0.10)")
    print(f"  Wave y-range          = [0, {2*AMPLITUDE:.6f}] m  (trough at y=0, crest at y=2a)")
    print(f"  Solid base y          = {Y_BOTTOM:.6f} m")
    print(f"  z end-cap margin      = {Z_MARGIN:.6f} m  (caps at z={-Z_MARGIN:.4f}, {LZ+Z_MARGIN:.4f})")

    tris = build_triangles()
    write_stl_binary(tris, OUT_FILE)

    expected = 8 * N_X + 4  # top+bottom+2 caps (2*N_X triangles each) + left/right (2 each)
    print(f"  Triangles written: {len(tris)}  (expected {expected})")
    print(f"  Output: {OUT_FILE}")
    print()
    print("Next step: mpirun -np <N> ../../../preProcessing/GenSDF/gensdf_mpi   (from this directory,")
    print("with a parameters.in pointing at wavy_wall.stl and this case's grid nx/ny/nz)")
