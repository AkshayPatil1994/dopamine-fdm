#!/usr/bin/env python3
"""Ting & Kirby (1994) spilling-breaker case for the two-fluid (VOF) solver: input files and bed geometry.

Experiment (Ting & Kirby 1994, Coastal Eng. 24, 51-80; Sect. 2, Fig. 1, Table 1): 40 m x 0.6 m x 1.0 m wave
tank, piston (bulkhead) wavemaker driven with a cnoidal signal (Goring 1978); still-water depth 0.40 m over
the flat part, then a plywood false bottom with slope 1:35. Spilling case: T = 2.0 s, H = 0.125 m in the
horizontal region (H0 = 0.127 m in deep water, Table 1), breaking at x_b = 6.40 m where d_b = 0.196 m.

Coordinates: x_TK is measured shoreward from the origin drawn in Fig. 1, where the still depth is 0.38 m
(not at the 0.40 m toe: the slope starts 0.02 m x 35 = 0.7 m upstream of the origin). The still depths d of
Table 1 (7 values for x_TK >= 5.9) fit a 1:35 plane with the slope starting at x_TK = -0.79 m (max.
deviation 0.3 cm; a free fit gives slope 1:35.1, d = 0 at x_TK = 13.2 m). The default X_TOE_TK is 0 (the
slope starts at the origin, as the paper's text says: x measured from the toe); this bed reproduces the
measured breaking point and wave-height evolution, whereas the -0.79 bed breaks ~1 m too early (--toe-tk
-0.79 or the coarse_exact_slope preset select it). The still-water shoreline is then at x_TK = 14 m.

Numerical flume (solver coordinate x = x_TK + X0, x = 0 at the wavemaker side of the domain, y vertical from the flat bed, z spanwise periodic):

    x = 0 ...... gen zone (wave_gen_len) ....... toe (x = X0 + X_TOE_TK) ........ beach, IBM, slope 1:35 ........ x = Lx (bed at SWL)
                                                                                       absorption zone (wave_abs_len) at the downstream end

The wave is the stream-function solution (order 16, wave_type = 2, zero mean volume flux = closed flume, wave_current_mode = 1) of H = 0.125 m, T = 2 s,
h = 0.4 m, which represents the cnoidal wave of the experiment (Ursell number ~27). The generation zone is one wavelength; the absorption zone ends the
domain at the shoreline and starts shoreward of the last gauge (x_TK = 11.5).

Known differences from the experiment: periodic span Lz = 0.2 m instead of the 0.6 m tank (no side walls; --lz), no piston wavemaker (generation zone
instead), no swash zone (the domain ends where the bed reaches the still-water level).

Outputs for the comparison with Figs. 4-6 of the paper (spilling breaker, post-processed by compare_simulation.py):
  * surface-elevation gauges (vof_gauges.dat, up to 32): the 8 LDA stations of Table 1 plus the stations of Fig. 4 (x_TK = -1.5 ... 11.5 m, GAUGES_FIG4_TK;
    the nominal positions are the digitised ones rounded to 0.5 m / 0.01 m)
  * 8 vertical line probes (U, V, W of the cell centres, mid-span) at the LDA stations: probes/line_N.bin -> phase-averaged undertow and rms
    velocity profiles, normalised with sqrt(g h) as in the paper (the gauge records give the local mean level / trough for masking).
  --gauges X1 X2 ... (x_TK) replaces the gauge list.

Usage:   python3 make_case.py CASE [--dx ... --dy ... --dz ... --gauges ...]      CASE = coarse | coarse_exact_slope | medium | fine | test (PRESETS)
Writes in ./CASE/: input_parameters, bed.stl, gensdf_parameters.in, stations.dat. setup.sh then builds SDF_in with GenSDF.
"""
import argparse, math, os, struct
import numpy as np

# ---- experiment ------------------------------------------------------------------------------------------------------------
H_WAVE, T_WAVE, DEPTH, SLOPE = 0.125, 2.0, 0.40, 1.0/35.0
X_TOE_TK = 0.0                                 # x_TK where the 1:35 slope starts (0 = toe at the origin; -0.79 = fit of the Table 1 depths, see docstring); overridden by --toe-tk
X_SHORE_TK = X_TOE_TK + DEPTH/SLOPE            # still-water shoreline (x_TK), 13.21 m
# measuring stations of the spilling case, Table 1 of the paper: x_TK, still depth d, mean depth h (with set-up)
STATIONS = [(-1.265, 0.400, 0.400), (5.945, 0.208, 0.212), (6.665, 0.185, 0.192), (7.275, 0.169, 0.175),
            (7.885, 0.152, 0.161), (8.495, 0.137, 0.148), (9.110, 0.119, 0.133), (9.725, 0.097, 0.112)]
X_BREAK_TK, D_BREAK = 6.4, 0.196
# wave gauges of Fig. 4 (digitised x of the circles, rounded to the nominal station); those within 0.03 m of an LDA station are dropped
GAUGES_FIG4_TK = [-1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.0, 6.4, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5]
MAX_GAUGES = 32                                 # MAX_WAVEGAUGES of the solver (src/global.f90)

# ---- numerical flume --------------------------------------------------------------------------------------------------------
X0 = 6.0                                        # solver x of x_TK = 0: gives 0.7 m of flat bed between the generation zone and the first gauge (x_TK = -1.5)
LY, LZ = 0.80, 0.20                             # height (0.4 m of air above SWL), periodic spanwise width (tank: 0.6 m; --lz)
GEN_LEN, ABS_LEN = 3.8, 1.5                     # generation zone = one wavelength (L = 3.7 m at h = 0.4, T = 2 s); absorption zone (x_TK > 11.7, behind the last gauge)
GAUGES_TK = sorted([s[0] for s in STATIONS] + [g for g in GAUGES_FIG4_TK if min(abs(g - s[0]) for s in STATIONS) > 0.03])

PRESETS = {   # cell sizes [m], end time [s], snapshot interval [s], line-probe interval [steps]; optional p_row, p_col (default 0 = automatic)
    # three resolutions for the grid study (dy = dx/2, dz = dx, except fine: dz = 0.01; cells across the 0.4 m depth: 40 / 80 / 160)
    'coarse': dict(dx=0.02,   dy=0.01,   dz=0.015,  t_end=80.0, tsave=0.2, line_freq=5),    # coarse mesh
    'medium': dict(dx=0.01,   dy=0.005,  dz=0.01,  t_end=80.0, tsave=0.2, line_freq=10),    # medium mesh
    'fine':   dict(dx=0.008,  dy=0.004, dz=0.008,  t_end=80.0, tsave=0.2, line_freq=20),    # fine mesh
    # same grid as 'coarse' with the slope at the Table 1 fit (x_TK = -0.79), for the bed A/B test (README); the default toe is x_TK = 0
    'coarse_exact_slope': dict(dx=0.02, dy=0.01, dz=0.015, t_end=80.0, tsave=0.2, line_freq=5, toe_tk=-0.79),
    'test':   dict(dx=0.04,   dy=0.02,   dz=0.04,  t_end=30.0, tsave=0.2, line_freq=2),     # set-up check only
}
# 80 s = 40 T; take the statistics (compare_simulation.py --t0) after >= 12 T of spin-up. All runs use the same physics and wave generation.


def write_bed_stl(path, lx, dz_pad=0.1):
    """Closed watertight STL of the solid below the bed: flat at y = 0 up to x = X0 + X_TOE_TK, then slope 1:35 to x = lx; extruded in z over
    [-dz_pad, LZ + dz_pad]; the polygon is extended 0.1 m below the floor and beyond both ends."""
    xt = X0 + X_TOE_TK
    yb_end = (lx - xt)*SLOPE
    poly = [(-0.1, -0.1), (lx + 0.1, -0.1), (lx + 0.1, yb_end + 0.1*SLOPE), (lx, yb_end), (xt, 0.0), (-0.1, 0.0)]
    z0, z1 = -dz_pad, LZ + dz_pad
    tris = []
    n = len(poly)
    c = np.mean(poly, axis=0)
    for i in range(n):                              # side walls
        a, b = poly[i], poly[(i + 1) % n]
        tris += [((a[0], a[1], z0), (b[0], b[1], z0), (b[0], b[1], z1)), ((a[0], a[1], z0), (b[0], b[1], z1), (a[0], a[1], z1))]
    for i in range(1, n - 1):                       # caps (fan from vertex 0; polygon is convex up to the toe kink, still a valid fan)
        a, b, d = poly[0], poly[i], poly[i + 1]
        tris += [((a[0], a[1], z0), (d[0], d[1], z0), (b[0], b[1], z0)), ((a[0], a[1], z1), (b[0], b[1], z1), (d[0], d[1], z1))]
    with open(path, 'wb') as f:
        f.write(b'tk94 bed'.ljust(80, b' '))
        f.write(struct.pack('<I', len(tris)))
        for t in tris:
            p = np.array(t, float)
            nv = np.cross(p[1] - p[0], p[2] - p[0])
            nv = nv/(np.linalg.norm(nv) + 1e-300)
            f.write(struct.pack('<12fH', *nv, *p.ravel(), 0))
    return len(tris)


def main():
    global X_TOE_TK, X_SHORE_TK, LZ
    ap = argparse.ArgumentParser()
    ap.add_argument('case', choices=list(PRESETS))
    for k in ('dx', 'dy', 'dz', 't_end', 'tsave'):
        ap.add_argument('--' + k, type=float)
    ap.add_argument('--line-freq', type=int, help='line-probe interval in time steps')
    ap.add_argument('--gauges', type=float, nargs='+', help=f'x_TK of up to {MAX_GAUGES} surface gauges (default: the 8 LDA stations of Table 1 + the Fig. 4 stations)')
    ap.add_argument('--lz', type=float, help='spanwise width [m] (periodic; the tank is 0.6 m wide)')
    ap.add_argument('--p-row', type=int, help='pencil decomposition rows (0 = automatic)')
    ap.add_argument('--p-col', type=int, help='pencil decomposition columns (0 = automatic)')
    ap.add_argument('--restart-step', type=int, help='resume from fields/tk94.<STEP> (needs --restart-time); probe samples still buffered at a crash are lost, '
                    'compare_simulation.py corrects the probe times after each restart (nstep_init in the logs)')
    ap.add_argument('--restart-time', type=float, help='simulation time of the restart file [s]')
    ap.add_argument('--out')
    ap.add_argument('--toe-tk', type=float, default=None, help='x_TK of the start of the slope (default 0 = origin at the toe; -0.79 = fit of the Table 1 depths, used by the coarse_exact_slope preset)')
    a = ap.parse_args()
    if a.lz: LZ = a.lz
    p = dict(PRESETS[a.case]); p.setdefault('p_row', 0); p.setdefault('p_col', 0)
    X_TOE_TK = a.toe_tk if a.toe_tk is not None else p.pop('toe_tk', X_TOE_TK); p.pop('toe_tk', None); X_SHORE_TK = X_TOE_TK + DEPTH/SLOPE
    p.update({k: getattr(a, k) for k in p if getattr(a, k, None) is not None})
    assert (a.restart_step is None) == (a.restart_time is None), '--restart-step and --restart-time go together'
    if a.restart_step is None:
        restart_nml, filein = 'restart = 0,', 'restart/tk94'
    else:
        restart_nml, filein = f'restart = 1,\n  t_start = {a.restart_time}, nstep_init = {a.restart_step},', f'fields/tk94.{a.restart_step}'
    gx_tk = a.gauges or GAUGES_TK
    assert len(gx_tk) <= MAX_GAUGES, f'the solver supports at most {MAX_GAUGES} gauges'
    out = a.out or os.path.join(os.path.dirname(os.path.abspath(__file__)), a.case)
    os.makedirs(os.path.join(out, 'restart'), exist_ok=True)
    os.makedirs(os.path.join(out, 'probes'), exist_ok=True)

    lx = X0 + X_SHORE_TK
    nx, ny, nz = round(lx/p['dx']), round(LY/p['dy']) + 1, round(LZ/p['dz'])
    if nx % 2: nx += 1                          # the solver needs an even nx too (Lx grows by one cell)
    if nz % 2:                                  # the solver needs an even nz: take the next even count and keep Lz (dz adjusts)
        nz += 1
        p['dz'] = LZ/(nz - 1)
    lx, ly, lz = (nx - 1)*p['dx'], (ny - 1)*p['dy'], (nz - 1)*p['dz']
    ntri = write_bed_stl(os.path.join(out, 'bed.stl'), lx)
    gauges = ', '.join('%.4f' % (g + X0) for g in gx_tk)
    zmid = 0.5*lz
    lines = ''.join(f"  line_dir({n+1}) = 'y', line_pos1({n+1}) = {s[0] + X0:.4f}, line_pos2({n+1}) = {zmid:.4f}, line_comps({n+1}) = 'UVW',\n"
                    f"  line_fileout({n+1}) = 'probes/line_{n+1}',\n" for n, s in enumerate(STATIONS))
    with open(os.path.join(out, 'stations.dat'), 'w') as f:
        f.write('# n  x_TK  x_solver  d_still(paper)  h_mean(paper)   [m]; gauges: ' + ' '.join('%.4f' % g for g in gx_tk) + f'\n'
                f'# case {a.case}  X0 {X0}  X_TOE_TK {X_TOE_TK}  SLOPE {SLOPE}  DEPTH {DEPTH}  T {T_WAVE}  H {H_WAVE}  line_freq {p["line_freq"]}  Lz {lz}\n')
        for n, s in enumerate(STATIONS):
            f.write('%d %.4f %.4f %.4f %.4f\n' % (n + 1, s[0], s[0] + X0, s[1], s[2]))

    open(os.path.join(out, 'input_parameters'), 'w').write(f"""! Ting & Kirby (1994) spilling breaker, T = {T_WAVE} s, H = {H_WAVE} m, h = {DEPTH} m, beach 1:35, LES.
! Generated by make_case.py ({a.case}): solver x = x_TK + {X0}; slope starts at x = {X0 + X_TOE_TK:.2f} (x_TK = {X_TOE_TK}); still-water shoreline at x = {X0 + X_SHORE_TK:.4f}; Lx = {lx:.4f}.
&DOMAIN
  nx = {nx}, ny = {ny}, nz = {nz},
  Lx = {lx:.10f}, Ly = {ly:.10f}, Lz = {lz:.10f},
  grid_type = 1, alpha_grid = 1.0, alpha_grid_z = 0.0, p_row = {p['p_row']}, p_col = {p['p_col']},
/
&PHYSICS
  nu = 1.0e-6, sgs_model = 1, flat_wall_model_flag = 1,
/
&NUMERICS
  dt = 1.0e-3, nsteps = -1, nsave = -1, tsave = {p['tsave']}, nmonitor = 1,
  cfl_adaptive = 1, cfl_target = 0.2, cfl_safety = 0.9, dt_max = 1.0, sim_end_time = {p['t_end']},
/
&BOUNDARY_CONDITIONS
  bc_face_ylo = 1, bc_face_yhi = 2, x_bc_type = 1, y_bc_type = 1, z_bc_type = 0,
/
&INFLOW
  inflow_type = 3,
/
&INITIAL_CONDITIONS
  Utarget = 0.0, ic_type = 3, noise_percent = 0.0, {restart_nml}
/
&STATISTICS
  n_lines = {len(STATIONS)}, line_freq = {p['line_freq']},
{lines}/
&IO
  filein = '{filein}', fileout = 'tk94',
/
&IBM
  ibm_input_mode = 1, ibm_wall_model_flag = 1, ibm_sdf_file = 'SDF_in',
/
&VOF
  vof_active = 1, vof_flow = 1, vof_ic_type = 1, vof_level = {DEPTH}, vof_inertia = 1,
  vof_pcg_iters = 800, vof_adv_iters = 100, vof_div_tol = 1.0e-8,
  vof_rho_l = 1000.0, vof_rho_g = 1.2, vof_nu_l = 1.0e-6, vof_nu_g = 1.5e-5, vof_grav = 9.81, vof_sigma = 0.072,
  vof_rk_mom = 2, vof_tr_floor = 1e-12,
/
&WAVES
  wave_type = 2, wave_height = {H_WAVE}, wave_period = {T_WAVE}, wave_sf_n = 16, wave_current_mode = 1,
  wave_gen_len = {GEN_LEN}, wave_abs_len = {ABS_LEN}, wave_relax_rate = 20.0, wave_ramp_time = {2*T_WAVE},
  wave_gauge_x = {gauges},
/
""")
    open(os.path.join(out, 'gensdf_parameters.in'), 'w').write(f"""# GenSDF configuration for the TK94 bed (written by make_case.py)
'bed.stl'
1d10  50  40
{nx - 1} {ny - 1} {nz - 1}
0.0 0.0 0.0
.true.
.true.
10
2
.false.
""")
    dmis = max(abs(s[1] - (DEPTH - max(0.0, (s[0] - X_TOE_TK)*SLOPE))) for s in STATIONS)
    print(f'toe at x_TK = {X_TOE_TK} (solver x = {X0 + X_TOE_TK}); max |d_paper(Table 1) - d_geometry| over the stations = {dmis*100:.1f} cm')
    print(f"{a.case}: {nx}x{ny}x{nz} cells ({(nx-1)*(ny-1)*(nz-1)/1e6:.2f} M), L = {lx:.3f} x {ly:.3f} x {lz:.3f}, "
          f"dx,dy,dz = {p['dx']},{p['dy']},{p['dz']}, bed STL {ntri} triangles -> {out}")


main()
