#!/usr/bin/env python3
"""Figs. 4-8 of Ting & Kirby (1994), spilling breaker, from one or several runs (e.g. coarse / medium / fine):

  Fig. 4  wave crest, trough and mean water level along the flume   (vof_gauges.dat; optionally the free surface of the snapshots)
  Fig. 5  time-mean horizontal velocity (undertow) u/sqrt(g h) vs (z - zeta_mean)/h at the 8 LDA stations
  Fig. 6  time-mean rms horizontal turbulent velocity sqrt(u'^2)/sqrt(g h) at the same stations
  Fig. 7  turbulent kinetic energy (k/gh)^(1/2), k = (1.33/2)(u'^2 + w'^2) as estimated in the paper, at stations 4-8 (needs the W component of the line probes)
  Fig. 8  phase-averaged surface elevation (zeta - zeta_mean)/h vs t/T at station 8 from the wave gauge there (phase aligned on the bore front)

Data: stations.dat, vof_gauges.dat and probes/line_N.bin written by a case of make_case.py.  Definitions as in the paper (sect. 2):
the phase average is taken at fixed phase t/T over whole periods after the spin-up (--t0), u' = u - phase average, the velocity is set
to zero when the probe level is above the free surface (gauge record at the same x), the time mean is the mean over the phase bins.
Only the resolved velocity is used (LES); it is sampled at one span position (mid-span line), so ~N periods x 1 point per phase bin.

usage:  python3 compare_simulation.py CASE [CASE ...] [--labels coarse medium fine] [--t0 24] [--out tk94] [--snap]
output: OUT_fig4.png OUT_fig5.png OUT_fig6.png OUT_fig7.png OUT_fig8.png and OUT_tables.txt
"""
import argparse, os, re, struct, sys
import numpy as np
from wavestats import phase_stats

G, T_WAVE = 9.81, 2.0
K_TK = 1.33/2          # TK94 estimate of k from the two measured components
ap = argparse.ArgumentParser()
ap.add_argument('cases', nargs='+')
ap.add_argument('--labels', nargs='+')
ap.add_argument('--t0', type=float, default=None, help='start of averaging window [s]; exclude trailing partial waves (default: last half of the run)')
ap.add_argument('--probe-step0', type=int, default=None, help='solver steps before the first probe sample (default: detected from nstep_init in the case logs)')
ap.add_argument('--no-restart-align', action='store_true', help='do not correct the probe sample times after restarts (see align_probe_restarts)')
ap.add_argument('--nphase', type=int, default=40, help='phase bins per period')
ap.add_argument('--nphase-wave', type=int, default=100, help='phase bins per period for the phase-averaged surface elevation (Fig. 8)')
ap.add_argument('--no-align-phase', action='store_true', help='Fig. 8: do not shift the simulated phase so that the zero up-crossing of the bore front coincides with the experiment (t/T = 0 is arbitrary in the runs)')
ap.add_argument('--snap', action='store_true', help='Fig. 4 also from the free surface of the snapshots (alpha columns; all x, crest/trough limited by tsave)')
ap.add_argument('--shift-x', type=float, default=0.0, help='diagnostic: translate the simulated Fig. 4 data by this x_TK offset [m] (Figs. 5-6 are per station and unaffected); the legend and tables are labelled')
ap.add_argument('--out', default='tk94')
args = ap.parse_args()
labels = args.labels or [os.path.basename(os.path.abspath(c)) for c in args.cases]
assert len(labels) == len(args.cases)

EXP = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'exp')   # digitised TK94 data (exp/digitise)

# reference numbers quoted in the text of the paper (spilling breaker)
REF_H = [(6.4, 0.1625, 0.0029), (11.5, 0.0338, 0.0187), (9.675, 0.064, 0.0134)]       # x_TK, H, mean level  (breaking point, bore, station 8)
REF_U8 = [(-0.298, -0.017), (-0.834, -0.10)]                   # station 8: u/sqrt(g h) at (z - zeta)/h  ("0.017 C" and "0.1 C" undertow)
REF_UR8 = [(-0.298, 0.078), (-0.834, 0.052)]                   # station 8: rms u'/sqrt(g h)


def read_stations(case):
    meta, rows = {}, []
    for ln in open(os.path.join(case, 'stations.dat')):
        if ln.startswith('# n'):
            meta['gauges'] = [float(v) for v in ln.split('gauges:')[1].split()]
        elif ln.startswith('# case'):
            w = ln.split()[2:]
            meta.update({w[i]: float(w[i + 1]) for i in range(1, len(w) - 1, 2)})
        elif not ln.startswith('#'):
            rows.append([float(v) for v in ln.split()])
    return meta, np.array(rows)


def bed_height(xtk, m):                               # y of the bed above the flat floor at x_TK
    return np.clip((xtk - m['X_TOE_TK'])*m['SLOPE'], 0, None)


def mean_ignoring_nan(values, axis):
    valid = np.isfinite(values)
    count = valid.sum(axis=axis)
    total = np.where(valid, values, 0.0).sum(axis=axis)
    return np.divide(total, count, out=np.full(total.shape, np.nan), where=count > 0)


def std_ignoring_nan(values, axis):
    mean = mean_ignoring_nan(values, axis)
    valid = np.isfinite(values)
    count = valid.sum(axis=axis)
    squared = np.where(valid, (values - np.expand_dims(mean, axis))**2, 0.0).sum(axis=axis)
    return np.sqrt(np.divide(squared, count - 1, out=np.zeros(squared.shape), where=count > 1))


def read_grid(case):
    """y cell centres from the header of a snapshot."""
    fdir = os.path.join(case, 'fields')
    f = sorted(f for f in os.listdir(fdir) if re.fullmatch(r'[A-Za-z0-9_]+\.\d+', f) and not f.startswith(('grid', 'geometry')))[0]
    with open(os.path.join(fdir, f), 'rb') as fh:
        ri = lambda: struct.unpack('>i', fh.read(4))[0]
        for _ in range(3):
            fh.seek(8*ri(), 1)
        xm = np.fromfile(fh, '>f8', ri()); ym = np.fromfile(fh, '>f8', ri())
    return xm, ym


def monotonic(t):
    """mask of the samples whose time is not below any earlier sample."""
    prev = np.concatenate(([-np.inf], np.maximum.accumulate(t)[:-1]))
    return t >= prev


def restart_rows(t):
    """raw row indices of vof_gauges.dat after resolving restarts: when the time steps back (a run resumed from an earlier
    step, appending to the file) the stale rows after that time are discarded, so row i of the result is solver step i + 1."""
    kept = []
    for i, ti in enumerate(t):
        if kept and ti < t[kept[-1]]:
            kept = kept[:np.searchsorted(t[kept], ti, side='left')]
        kept.append(i)
    return np.array(kept)


def probe_step_offset(case, nrows, last_step):
    """solver steps preceding the first probe sample. A resumed run counts its steps (and writes the probe files) from 1 again,
    so the probe step must be shifted by nstep_init of the resume; candidates are 0 and nstep_init from the case logs, and the one
    whose last probe sample lands closest to (not after) the end of the gauge record is taken."""
    if args.probe_step0 is not None:
        return args.probe_step0
    cands = {0}
    for fn in ('resume.log', 'run.log'):
        f = os.path.join(case, fn)
        if os.path.exists(f):
            cands.update(int(v) for v in re.findall(r'nstep_init\s*=\s*(\d+)', open(f, errors='replace').read()))
    fits = [(nrows - (c + last_step), c) for c in cands if 0 <= nrows - (c + last_step)]
    return min(fits)[1] if fits else 0


def restart_times(case, t):
    """times of the restart snapshots (nstep_init > 0 in the case logs) on the gauge time axis"""
    steps = set()
    for fn in ('resume.log', 'run.log'):
        f = os.path.join(case, fn)
        if os.path.exists(f):
            steps.update(int(v) for v in re.findall(r'nstep_init\s*=\s*(\d+)', open(f, errors='replace').read()))
    return sorted(t[n - 1] for n in steps if 0 < n <= len(t))


def fundamental_phase(tt, x):
    w = 2*np.pi*tt/T_WAVE
    c = np.linalg.lstsq(np.c_[np.cos(w), np.sin(w), np.ones_like(w)], x, rcond=None)[0]
    return np.arctan2(c[1], c[0])


def restart_gap(tp, ucol, t, eta, tr):
    """Probe files carry no time stamp. A crashed run loses the probe samples still buffered, so after a restart the n-th sample is
    from a LATER step than n*line_freq: the times assigned from the step count are too early by dt, while the gauge file is truncated
    exactly at the restart snapshot (time tr) and is not affected. For a periodic wave the phase of u (below the trough) relative to the
    surface elevation is constant, so dt follows from its change across the restart (5 periods before / after). Returns dt [s] or None."""
    def dphi(a, b):
        q = (tp >= a) & (tp < b)
        if q.sum() < 100:
            return None
        return fundamental_phase(tp[q], ucol[q]) - fundamental_phase(tp[q], np.interp(tp[q], t, eta))
    pre, post = dphi(tr - 6*T_WAVE, tr - T_WAVE), dphi(tr, tr + 5*T_WAVE)
    if pre is None or post is None:
        return None
    return ((pre - post + np.pi) % (2*np.pi) - np.pi)*T_WAVE/(2*np.pi)


def load_probe(case, n, m, t, comp=0):
    """raw line probe n (1-based): (U[sample, y, component comp], solver step of every sample), None if absent (comp 0 = u, 2 = w for line_comps = 'UVW')"""
    f = os.path.join(case, 'probes', f'line_{n}')
    if not os.path.exists(f + '.bin'):
        return None
    meta = {k.strip(): v.strip() for k, v in (l.split('=', 1) for l in open(f + '_meta.txt').read().splitlines() if '=' in l)}
    npts, nc = int(meta['npts']), int(meta['ncomp'])
    raw = np.fromfile(f + '.bin', '>f8')
    if comp >= nc:
        return None
    ns = raw.size//(npts*nc)
    U = raw[:ns*npts*nc].reshape(ns, npts, nc)[:, :, comp]
    step = (np.arange(ns) + 1)*int(m['line_freq'])
    step = step + probe_step_offset(case, len(t), step[-1])
    return U, step


def restart_step_shifts(case, d, ym):
    """{restart step: steps lost from the probe files}. Estimated at every station (restart_gap), converted to steps and combined
    by the median (the loss is the same for all lines); rounded to whole probe samples."""
    m, t, eta = d['meta'], d['t'], d['eta']
    gx = np.array(m['gauges']); lf = int(m['line_freq'])
    shifts = {}
    if args.no_restart_align:
        return shifts
    for tr in restart_times(case, t):
        nr = int(np.searchsorted(t, tr)) + 1
        est = []
        for n, row in enumerate(d['st']):
            pr = load_probe(case, n + 1, m, t)
            if pr is None:
                continue
            U, step = pr; ok = step <= len(t)
            j = int(np.argmin(abs(gx - row[1])))
            kref = int(np.argmin(abs(ym - 0.5*(bed_height(row[1], m) + m['DEPTH'] + eta[:, j].min() - (ym[1] - ym[0])))))   # always wet: below the trough
            dt = restart_gap(t[step[ok] - 1], U[ok, kref], t, eta[:, j], tr)
            if dt is not None and dt > 0.005:
                i0, i1 = np.searchsorted(t, [tr, tr + 5*T_WAVE])
                est.append((n + 1, dt, dt*(i1 - i0)/(5*T_WAVE)))      # steps lost = time lost x mean step rate after the restart
        if not est:
            continue
        S = int(np.median([e[2] for e in est]))
        S = lf*int(round(S/lf))
        shifts[nr] = S
        print(f'NOTE {case}: probe samples after the restart at t = {tr:.2f} s (step {nr}) are {S} steps later than their count gives '
              f'(per station, s: {", ".join(f"{e[0]}:{e[1]:.3f}" for e in est)}); step-based correction, spread of the per-station estimates {min(e[2] for e in est):.0f}..{max(e[2] for e in est):.0f} steps', file=sys.stderr)
    return shifts


def gauge_stats(t, eta, t0, t_last=None):
    t_last = t[-1] if t_last is None else min(t[-1], t_last)
    s = phase_stats(t, eta, t0, t_last, T_WAVE)
    std = lambda values: np.std(values, ddof=1) if len(values) > 1 else 0.0
    pm, pc, pt = s['period_mean'], s['period_crest'], s['period_trough']
    return (s['mean'], s['crest'], s['trough'], s['nper'], t0,
            std(pm), std(pc - pm), std(pt - pm), std(pc - pt), s['t_end'])


def load_case(case):
    m, st = read_stations(case)
    g = np.loadtxt(os.path.join(case, 'vof_gauges.dat'))
    rows = restart_rows(g[:, 0])
    if len(rows) < len(g):
        print(f'WARNING {case}: vof_gauges.dat contains restarts; dropped {len(g) - len(rows)} stale of {len(g)} rows', file=sys.stderr)
    if g.shape[1] - 1 != len(m['gauges']):
        raise SystemExit(f'{case}: vof_gauges.dat has {g.shape[1] - 1} gauge columns but stations.dat lists {len(m["gauges"])} gauges')
    t, eta = g[rows, 0], g[rows, 1:]            # row index = solver step - 1 (used for the probe step -> time lookup)
    t0 = args.t0 if args.t0 is not None else 0.5*t[-1]
    stats = gauge_stats(t, eta[:, 0], t0)
    out = dict(meta=m, st=st, t=t, eta=eta, t0=t0, t_end=stats[-1])
    # Fig. 4 from the gauges
    gx = m['gauges']
    gstat = []
    for j, x in enumerate(gx):
        stats = gauge_stats(t, eta[:, j], t0)
        gstat.append((x, *stats[:3], *stats[5:]))
    out['gstat'] = np.asarray(gstat)
    return out


def snapshot_surface(case, d):
    """free-surface elevation (above SWL) of every snapshot at all x: eta = y_bed + sum(alpha dy over fluid cells) - 0.4 (z-mean)."""
    m = d['meta']
    fdir = os.path.join(case, 'fields')
    pre = re.compile(r'([A-Za-z0-9_]+)\.(\d+)$')
    files = sorted((int(pre.match(f).group(2)), f) for f in os.listdir(fdir) if pre.match(f) and not f.startswith(('grid', 'geometry')))
    tsave = 0.2
    X0 = m['X0']; level = m['DEPTH']
    res = []
    for i, (step, f) in enumerate(files):
        tt = (i + 1)*tsave
        if tt < d['t0'] - 1e-9 or tt >= d['t_end'] - 1e-9:
            continue
        with open(os.path.join(fdir, f), 'rb') as fh:
            size = os.path.getsize(os.path.join(fdir, f))
            ri = lambda: struct.unpack('>i', fh.read(4))[0]
            for _ in range(3):
                fh.seek(8*ri(), 1)
            xm = np.fromfile(fh, '>f8', ri()); ym = np.fromfile(fh, '>f8', ri()); zm = np.fromfile(fh, '>f8', ri())
            blocks = []
            while fh.tell() < size:
                dd = struct.unpack('>3i', fh.read(12))
                blocks.append((dd, fh.tell())); fh.seek(8*dd[0]*dd[1]*dd[2], 1)
            dd, pos = blocks[-1]
            fh.seek(pos)
            a = np.fromfile(fh, '>f8', dd[0]*dd[1]*dd[2]).reshape(dd[2], dd[1], dd[0])
        nz, ny, nx = len(zm), len(ym), len(xm)
        a = a[1:1 + nz, 1:1 + ny, 1:1 + nx].mean(axis=0)          # (ny, nx), z-mean
        dy = np.gradient(ym)
        yb = bed_height(xm - X0, m)
        fl = ym[:, None] > yb[None, :]
        res.append((tt, xm - X0, yb + (a*dy[:, None]*fl).sum(axis=0) - level))
    return res


def stations_profiles(case, d):
    """phase-averaged undertow and rms u' at the stations. Returns list of dicts."""
    m, st = d['meta'], d['st']
    t, eta = d['t'], d['eta']
    gx = np.array(m['gauges'])
    xm, ym = read_grid(case)
    shifts = restart_step_shifts(case, d, ym)
    out = []
    for n, row in enumerate(st):
        x = row[1]
        pr = load_probe(case, n + 1, m, t)
        if pr is None:
            out.append(None); continue
        U, step = pr
        pw = load_probe(case, n + 1, m, t, comp=2)
        W = pw[0] if pw is not None else None
        npts = U.shape[1]
        for nr, S in shifts.items():
            step = np.where(step >= nr - S, step + S, step)               # samples written after the restart
        ok = step <= len(t)
        U, step = U[ok], step[ok]
        W = W[ok] if W is not None else None
        tp = t[step - 1]
        j = int(np.argmin(abs(gx - x)))
        if abs(gx[j] - x) > 1e-3:
            raise SystemExit(f'no gauge at station {x}: the masking of the dry probe points needs one (run with the default --gauges)')
        ok = monotonic(tp)                                               # drop probe samples written across a restart overlap
        if not ok.all():
            print(f'WARNING station {n + 1}: dropped {(~ok).sum()} of {len(ok)} probe samples with non-monotonic time', file=sys.stderr)
        U, tp, step = U[ok], tp[ok], step[ok]
        W = W[ok] if W is not None else None
        et = eta[step - 1, j]
        # the window ends with the probe record (may be shorter than the gauge record); gauge stats use the same window
        t0_eff = d['t0']
        if tp[0] > t0_eff + 1e-6:                                        # probe record starts later: begin at the next whole period
            t0_eff += np.ceil((tp[0] - t0_eff)/T_WAVE - 1e-9)*T_WAVE
        stats = gauge_stats(t, eta[:, j], t0_eff, t_last=tp[-1])
        if t0_eff > d['t0'] + 1e-6 or tp[-1] < d['t'][-1] - T_WAVE:
            print(f'WARNING station {n + 1}: probes cover t = {tp[0]:.2f}..{tp[-1]:.2f} s, gauges {d["t0"]:.2f}..{d["t"][-1]:.2f} s; averaging over {stats[3]} periods from t0 = {t0_eff:.2f} s', file=sys.stderr)
        mean_lvl, cr, tr, nper, t0 = stats[:5]
        t_end = stats[-1]
        s = (tp >= t0 - 1e-9) & (tp < t_end - 1e-9)
        U, tp, et = U[s], tp[s], et[s]
        W = W[s] if W is not None else None
        if npts != len(ym):
            raise SystemExit(f'station {n + 1}: line probe has {npts} points but the grid has {len(ym)} y cells')
        yc = ym[:npts]
        wet = (yc[None, :] - m['DEPTH']) < et[:, None]
        ph = np.minimum((((tp - t0)/T_WAVE) % 1.0*args.nphase).astype(int), args.nphase - 1)
        period = np.floor((tp - t0)/T_WAVE + 1e-9).clip(0, nper - 1).astype(int)
        uph = np.full((args.nphase, npts), np.nan); rph = np.full_like(uph, np.nan)      # empty bins stay NaN and are skipped in the mean
        uph_period = np.full((nper, args.nphase, npts), np.nan)
        rph_period = np.full_like(uph_period, np.nan)
        kph = np.full_like(uph, np.nan); kph_period = np.full_like(uph_period, np.nan)      # phase variance of w (k needs u'^2 + w'^2)
        for b in range(args.nphase):
            q = ph == b
            if not q.any():
                continue
            Um = U[q]*wet[q]                                          # velocity set to zero above the free surface
            uph[b] = Um.mean(axis=0)
            nb = q.sum()
            if nb > 1:                                                # unbiased (ddof = 1) variance of the masked velocity about its phase mean
                rph[b] = np.sqrt(((Um - uph[b])**2).sum(axis=0)/(nb - 1))
            pq = period[q]
            for pidx in range(nper):
                qp = pq == pidx
                if qp.any():
                    uph_period[pidx, b] = Um[qp].mean(axis=0)
                    rph_period[pidx, b] = np.sqrt(((Um[qp] - uph[b])**2).mean(axis=0))   # per-period spread only (std bars), not bias-corrected
            if W is not None:
                Wm = W[q]*wet[q]; wm = Wm.mean(axis=0)
                if nb > 1:
                    kph[b] = ((Wm - wm)**2).sum(axis=0)/(nb - 1)                         # phase variance of w
                for pidx in range(nper):
                    qp = pq == pidx
                    if qp.any():
                        kph_period[pidx, b] = ((Wm[qp] - wm)**2).mean(axis=0)
        cycle_u = mean_ignoring_nan(uph_period, axis=1)
        cycle_urms = mean_ignoring_nan(rph_period, axis=1)
        std_u = std_ignoring_nan(cycle_u, axis=0)
        std_urms = std_ignoring_nan(cycle_urms, axis=0)
        if W is not None:      # k = (1.33/2)(u'^2 + w'^2) as estimated in TK94 (no transverse component); time mean over the phase bins
            k_mean = K_TK*(mean_ignoring_nan(rph**2, axis=0) + mean_ignoring_nan(kph, axis=0))
            k_cycle = K_TK*(mean_ignoring_nan(rph_period**2, axis=1) + mean_ignoring_nan(kph_period, axis=1))
            k_std = std_ignoring_nan(k_cycle, axis=0)
        else:
            k_mean = k_std = np.full(npts, np.nan)
        dstill = m['DEPTH'] - bed_height(x, m)                         # x = x_TK (stations.dat col 2): geometry of the run
        out.append(dict(x=row[1], d=dstill, h=dstill + mean_lvl, zeta=mean_lvl, trough=tr, crest=cr, nper=nper, ns=len(tp),
                        y=yc, ybed=bed_height(row[1], m), u=mean_ignoring_nan(uph, axis=0), urms=mean_ignoring_nan(rph, axis=0), depth=m['DEPTH'],
                        u_std=std_u, urms_std=std_urms, k=k_mean, k_std=k_std))
    return out


def main():
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt, matplotlib.colors as mc, colorsys
    plt.rcParams.update({
        'font.size': 9, 'axes.labelsize': 10, 'axes.titlesize': 10,
        'xtick.labelsize': 8, 'ytick.labelsize': 8, 'legend.fontsize': 8,
    })
    # Okabe-Ito colourblind-safe palette (black is kept for the experiment, yellow is dropped: poor contrast on white);
    # the marker also changes with the case so that runs stay distinguishable in greyscale and with many resolutions
    palette = ['#D55E00', '#0072B2', '#009E73', '#CC79A7', '#E69F00', '#56B4E9']      # vermillion, blue, bluish green, reddish purple, orange, sky blue
    markers = ['D', 's', '^', 'v', 'P', 'X']                                          # 'o' is the experiment
    style = lambda i: (palette[i % len(palette)], markers[i % len(markers)] if len(args.cases) > 1 else 'D')
    def darker(c, f=0.6):
        r, g, b = mc.to_rgb(c); hh, ll, ss = colorsys.rgb_to_hls(r, g, b)
        return colorsys.hls_to_rgb(hh, ll*f, ss)
    LSTY = ['-', '--', '-.', ':']

    def shade(ax, i, x, y, err, vertical=True, label=None, snap=False, every=1):
        """mean line (+ markers) with a +-err band: translucent fill of the case colour, bordered by a darker line of the same hue
        (solid for the first case, dashed for the second, ...). vertical: err is in y (Fig. 4, y(x)); else in x (profiles, z(x))."""
        col, mk = style(i); edge = darker(col); ls = LSTY[i % len(LSTY)] if len(args.cases) > 1 else '-'
        x, y, err = (np.asarray(v, float) for v in (x, y, np.nan_to_num(err)))
        if vertical:
            ax.fill_between(x, y - err, y + err, color=col, alpha=0.12 if snap else 0.28, lw=0)
            for sgn in (-1, 1): ax.plot(x, y + sgn*err, ls=ls, color=edge, lw=0.8, alpha=0.7 if snap else 1)
        else:
            ax.fill_betweenx(y, x - err, x + err, color=col, alpha=0.28, lw=0)
            for sgn in (-1, 1): ax.plot(x + sgn*err, y, ls=ls, color=edge, lw=0.8)
        if snap:
            ax.plot(x, y, ls=':', color=edge, lw=1.0)
        else:
            ax.plot(x, y, ls=ls, color=edge, lw=1.3, marker=mk, ms=3.5, markevery=every, markerfacecolor=col, markeredgecolor=edge, label=label)
    data = [load_case(c) for c in args.cases]
    for d in data:
        d['gstat'] = d['gstat'].copy(); d['gstat'][:, 0] += args.shift_x
    prof = [stations_profiles(c, d) for c, d in zip(args.cases, data)]
    tab = open(args.out + '_tables.txt', 'w')

    # ---- Fig. 4
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
    sh = f' (x shifted {args.shift_x:+.2f} m)' if args.shift_x else ''
    for i, (c, d) in enumerate(zip(args.cases, data)):
        gs = d['gstat']
        col, mk = style(i)
        shade(ax[0], i, gs[:, 0], 100*(gs[:, 2] - gs[:, 1]), 100*gs[:, 5], label=labels[i] + sh)
        shade(ax[0], i, gs[:, 0], 100*(gs[:, 3] - gs[:, 1]), 100*gs[:, 6])
        shade(ax[0], i, gs[:, 0], 100*gs[:, 1], 100*gs[:, 4])
        shade(ax[1], i, gs[:, 0], 100*(gs[:, 2] - gs[:, 3]), 100*gs[:, 7], label=labels[i] + sh)
        tab.write(f'# Fig. 4, {labels[i]}{sh}: x_TK  mean  crest  trough  H [m] (relative to SWL; {int((d["t_end"]-d["t0"])/T_WAVE)} whole periods from t0 = {d["t0"]:.1f} to {d["t_end"]:.1f} s)\n')
        for r in gs:
            tab.write('%8.3f %9.4f %9.4f %9.4f %9.4f\n' % (r[0], r[1], r[2], r[3], r[2] - r[3]))
        if args.snap:
            sn = snapshot_surface(c, d)
            if sn:
                x = sn[0][1] + args.shift_x
                by_period = {}
                for tt, _, surface in sn:
                    pidx = int((tt - d['t0'])//T_WAVE)
                    if pidx >= 0:
                        by_period.setdefault(pidx, []).append(surface)
                period_surfaces = [np.array(by_period[p]) for p in sorted(by_period)]
                if period_surfaces:
                    means = np.array([e.mean(axis=0) for e in period_surfaces])
                    crests = np.array([e.max(axis=0) for e in period_surfaces])
                    troughs = np.array([e.min(axis=0) for e in period_surfaces])
                    snapshot_stats = (crests - means, troughs - means, means, crests - troughs)
                    for axis, values, error in (
                        (ax[0], snapshot_stats[0], std_ignoring_nan(snapshot_stats[0], axis=0)),
                        (ax[0], snapshot_stats[1], std_ignoring_nan(snapshot_stats[1], axis=0)),
                        (ax[0], snapshot_stats[2], std_ignoring_nan(snapshot_stats[2], axis=0)),
                        (ax[1], snapshot_stats[3], std_ignoring_nan(snapshot_stats[3], axis=0)),
                    ):
                        shade(axis, i, x, 100*values.mean(axis=0), 100*error, snap=True)
    e4 = os.path.join(EXP, 'tk94_fig4_spilling.dat')
    if os.path.exists(e4):
        e = np.loadtxt(e4)
        ax[0].plot(e[:, 0], e[:, 1], 'ko', mfc='none', ms=4, label='TK94 (digitised)')
        ax[0].plot(e[:, 0], e[:, 2], 'ko', mfc='none', ms=4); ax[0].plot(e[:, 0], e[:, 3], 'ko', mfc='none', ms=4)
        ax[1].plot(e[:, 0], e[:, 2] - e[:, 3], 'ko', mfc='none', ms=4, label='TK94 (digitised)')
        for i, d in enumerate(data):                   # error at the gauges (mean level, crest, trough) against the experiment
            gs = d['gstat']
            for r in gs:                               # experimental curves interpolated linearly to the gauge position
                if e[0, 0] <= r[0] <= e[-1, 0]:
                    em, ec, et = (np.interp(r[0], e[:, 0], e[:, k]) for k in (1, 2, 3))
                    tab.write(f'# Fig. 4 vs TK94, {labels[i]}: x_TK={r[0]:.3f}  mean {100*r[1]:.2f}/{em:.2f}  '
                              f'crest-mean {100*(r[2]-r[1]):.2f}/{ec:.2f}  trough-mean {100*(r[3]-r[1]):.2f}/{et:.2f}  [cm, LES/exp]\n')
    else:
        ax[0].plot([r[0] for r in REF_H], [100*r[2] for r in REF_H], 'ks', mfc='none', label='paper (text values)')
        ax[1].plot([r[0] for r in REF_H], [100*r[1] for r in REF_H], 'ks', mfc='none', label='paper (text values)')
    ax[0].set_ylabel(r'$\zeta_{max}-\bar\zeta,\ \zeta_{min}-\bar\zeta,\ \bar\zeta$  [cm]'); ax[1].set_ylabel('H [cm]')
    for q in ax:
        q.set_xlabel(r'$x_{TK}$ [m]'); q.axvline(6.4, c='0.7', ls=':'); q.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(args.out + '_fig4.png', dpi=300, bbox_inches='tight')

    # ---- Figs. 5, 6
    EXPD = {k: np.loadtxt(os.path.join(EXP, fn)) for k, fn in ((5, 'tk94_fig5_undertow_spilling.dat'), (6, 'tk94_fig6_rms_spilling.dat')) if os.path.exists(os.path.join(EXP, fn))}
    for fig_no, key, xl, ref, lim in ((5, 'u', r'$\bar u/\sqrt{gh}$', REF_U8, None), (6, 'urms', r"$\sqrt{\widetilde{u'^2}}/\sqrt{gh}$", REF_UR8, (0, None))):
        fig, axs = plt.subplots(2, 4, figsize=(13, 7), sharey=True)
        for n, q in enumerate(axs.ravel()):
            for i, pr in enumerate(prof):
                p = pr[n]
                if p is None:
                    continue
                zz = (p['y'] - p['depth'] - p['zeta'])/p['h']
                fl = (zz > -1.0) & (zz < 0.0)
                fl &= p['y'] > p['ybed'] + (p['y'][1] - p['y'][0])
                v = p[key]/np.sqrt(G*p['h'])
                shade(q, i, v[fl], zz[fl], p[f'{key}_std'][fl]/np.sqrt(G*p['h']), vertical=False, label=labels[i], every=3)
                if fig_no == 5:
                    dyp = p['y'][1] - p['y'][0]; wetc = p['y'] > p['ybed'] + dyp
                    tab.write(f'# station {n + 1} {labels[i]}: net volume flux Q = {np.nansum(p["u"][wetc])*dyp*1e3:.2f} L/s/m (zero in a closed flume; wave crest transport ~ g H^2/8c ~ 10 L/s/m)\n')
                    tab.write(f'# station {n + 1} {labels[i]}: x_TK={p["x"]} d={p["d"]:.3f} zeta={p["zeta"]:.4f} h={p["h"]:.3f} trough={p["trough"]:.4f}'
                              f' ({p["nper"]} periods, {p["ns"]} samples)\n')
            ex = EXPD.get(fig_no)
            if ex is not None:
                s = ex[ex[:, 0] == n + 1]
                q.plot(s[:, 5], s[:, 4], 'ko', mfc='none', ms=4, label='TK94 (digitised)')
                for i, pr in enumerate(prof):
                    p = pr[n]
                    if p is None: continue
                    zz = (p['y'] - p['depth'] - p['zeta'])/p['h']; v = p[key]/np.sqrt(G*p['h'])
                    fl = (zz > -1) & (zz < 0) & (p['y'] > p['ybed'] + (p['y'][1] - p['y'][0]))
                    zs = s[:, 4]
                    ok = (zs > zz[fl].min()) & (zs < zz[fl].max())
                    if ok.any():
                        vi = np.interp(zs[ok], zz[fl], v[fl])
                        err = np.sqrt(np.mean((vi - s[ok, 5])**2))
                        tab.write(f'# Fig. {fig_no} vs TK94, station {n + 1} {labels[i]}: rms error {err:.4f} (exp rms {np.sqrt(np.mean(s[ok, 5]**2)):.4f}) in u/sqrt(gh), {ok.sum()} points between the bed and the mean level\n')
            elif n == 7:
                q.plot([r[1] for r in ref], [r[0] for r in ref], 'ks', mfc='none', label='paper (text)')
            pp = next((pr[n] for pr in prof if pr[n] is not None), None)       # depth used for the normalisation (geometry of the run, not Table 1)
            q.set_title(f'x_TK = {data[0]["st"][n][1]:.3f} m' + (f', d = {pp["d"]:.3f} m' if pp else ''), fontsize=9)
            q.set_xlabel(xl); q.set_ylim(-1, 0)
            if fig_no == 5: q.set_xlim(-0.2, 0.2)
            else: q.set_xlim(left=0)                                  # rms >= 0; upper limit autoscaled per panel
            if fig_no == 5: q.axvline(0, c='0.8', lw=0.8)
            if n % 4 == 0: q.set_ylabel(r'$(z-\bar\zeta)/h$')
        axs[0, 0].legend(fontsize=8)
        fig.tight_layout(); fig.savefig(f'{args.out}_fig{fig_no}.png', dpi=300, bbox_inches='tight')

    # ---- Fig. 7: turbulent kinetic energy k = (1.33/2)(u'^2 + w'^2), (k/gh)^(1/2), at the 5 stations of Fig. 7a (stations 4-8)
    f7 = os.path.join(EXP, 'tk94_fig7_tke_spilling.dat')
    ex7 = np.loadtxt(f7) if os.path.exists(f7) else None
    exp_mk = {4: ('o', 'none'), 5: ('o', 'k'), 6: ('P', 'k'), 7: ('*', 'k'), 8: ('X', 'k')}       # as in Fig. 7a of the paper (marker = station)
    st7 = [3, 4, 5, 6, 7]
    fig, axs = plt.subplots(2, 3, figsize=(11, 7), sharey=True)
    for q, n in zip(axs.ravel(), st7 + [None]):
        stl = st7 if n is None else [n]
        for i, pr in enumerate(prof if n is not None else []):      # the overlay panel shows the experiment only (Fig. 7a of the paper)
            for nn in stl:
                p = pr[nn]
                if p is None or not np.isfinite(p['k']).any():
                    continue
                zz = (p['y'] - p['depth'] - p['zeta'])/p['h']
                fl = (zz > -1.0) & (zz < 0.0) & (p['y'] > p['ybed'] + (p['y'][1] - p['y'][0]))
                v = np.sqrt(p['k']/(G*p['h']))
                with np.errstate(divide='ignore', invalid='ignore'):
                    verr = 0.5*v*p['k_std']/p['k']                     # d sqrt(k) = sqrt(k) dk/(2k)
                shade(q, i, v[fl], zz[fl], verr[fl], vertical=False, label=labels[i], every=3)
        if ex7 is not None:
            for nn in stl:
                s7 = ex7[ex7[:, 0] == nn + 1]
                mk, fc = exp_mk[nn + 1]
                q.plot(s7[:, 5], s7[:, 4], mk, color='k', mfc=fc, ms=4 if mk == 'o' else 5, ls='none', label='TK94 (digitised)' if n is not None else None)
                if n is not None:
                    for i, pr in enumerate(prof):
                        p = pr[nn]
                        if p is None or not np.isfinite(p['k']).any(): continue
                        zz = (p['y'] - p['depth'] - p['zeta'])/p['h']; v = np.sqrt(p['k']/(G*p['h']))
                        fl = (zz > -1) & (zz < 0) & (p['y'] > p['ybed'] + (p['y'][1] - p['y'][0]))
                        zs = s7[:, 4]; ok = (zs > zz[fl].min()) & (zs < zz[fl].max())
                        if ok.any():
                            err = np.sqrt(np.mean((np.interp(zs[ok], zz[fl], v[fl]) - s7[ok, 5])**2))
                            tab.write(f'# Fig. 7 vs TK94, station {nn + 1} {labels[i]}: rms error {err:.4f} (exp rms {np.sqrt(np.mean(s7[ok, 5]**2)):.4f}) in (k/gh)^1/2, {ok.sum()} points\n')
        if n is None:
            q.set_title('all stations (as Fig. 7a of the paper)', fontsize=9)
        else:
            pp = next((pr[n] for pr in prof if pr[n] is not None), None)
            q.set_title(f'x_TK = {data[0]["st"][n][1]:.3f} m' + (f', d = {pp["d"]:.3f} m' if pp else ''), fontsize=9)
        q.set_xlabel(r'$(\bar k/gh)^{1/2}$'); q.set_xlim(0, 0.1); q.set_ylim(-1, 0)
        if q in axs[:, 0]: q.set_ylabel(r'$(z-\bar\zeta)/h$')
    axs[0, 0].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(f'{args.out}_fig7.png', dpi=300, bbox_inches='tight')

    # ---- Fig. 8: phase-averaged surface elevation (zeta - zeta_mean)/h at station 8 ((x - xb)/hb = 16.709, h/hb = 0.563) from the wave gauge there
    f8 = os.path.join(EXP, 'tk94_fig8_phase_avg_spilling.dat')
    e8 = np.loadtxt(f8) if os.path.exists(f8) else None
    fig, ax8 = plt.subplots(figsize=(7, 4.5))
    nb = args.nphase_wave

    def front_phase(ph, y):
        """phase [0,1) of the steepest zero up-crossing (the front of the bore) of the periodic curve y(ph)"""
        y1 = np.roll(y, -1); up = np.where((y < 0) & (y1 >= 0))[0]
        if len(up) == 0: return np.nan
        k = up[np.argmax((y1 - y)[up])]; f = -y[k]/(y1[k] - y[k])
        return ((ph[k] + f*(ph[1] - ph[0])) % 1.0)

    if e8 is not None:
        te = np.argsort(e8[:, 0], kind='stable'); ge = np.linspace(0, 1, 2001)[:-1]
        ye = np.interp(ge, e8[te, 0], e8[te, 1]); phi_e = front_phase(ge, ye)
        ax8.plot(e8[:, 0], e8[:, 1], 'k-', lw=1.4, label='TK94 (digitised)', zorder=5)
    for i, d in enumerate(data):
        m = d['meta']; t, eta = d['t'], d['eta']
        xs = d['st'][7][1]; j = int(np.argmin(abs(np.array(m['gauges']) - xs)))
        nper = int((d['t_end'] - d['t0'])/T_WAVE + 1e-9)
        if nper < 1:
            print(f'WARNING {labels[i]}: no whole period for Fig. 8', file=sys.stderr); continue
        sel = (t >= d['t0']) & (t < d['t0'] + nper*T_WAVE)
        tt, ee = t[sel] - d['t0'], eta[sel, j]
        per = np.minimum((tt//T_WAVE).astype(int), nper - 1); ph = np.minimum(((tt/T_WAVE) % 1.0*nb).astype(int), nb - 1)
        cnt = np.bincount(per*nb + ph, minlength=nper*nb).astype(float); tot = np.bincount(per*nb + ph, weights=ee, minlength=nper*nb)
        cyc = np.divide(tot, cnt, out=np.full(nper*nb, np.nan), where=cnt > 0).reshape(nper, nb)          # (period, phase bin)
        zmean = np.nanmean(cyc); h = m['DEPTH'] - bed_height(xs, m) + zmean
        mean_c, std_c = np.nanmean(cyc, axis=0)/h - zmean/h, (np.nanstd(cyc, axis=0, ddof=1)/h if nper > 1 else np.zeros(nb))
        phc = (np.arange(nb) + 0.5)/nb
        shift = 0.0
        if e8 is not None and not args.no_align_phase:
            shift = (phi_e - front_phase(phc, mean_c) + 0.5) % 1.0 - 0.5
        pp = (phc + shift) % 1.0; o = np.argsort(pp)
        shade(ax8, i, pp[o], mean_c[o], std_c[o], label=f'{labels[i]} ({nper} periods)', every=max(1, nb//20))
        msg = f'# Fig. 8, {labels[i]}: gauge x_TK = {m["gauges"][j]:.3f} m, h = {h:.4f} m, zeta_mean = {zmean:.4f} m, {nper} whole periods; phase shift applied {shift:+.4f} T'
        msg += f'; crest {mean_c.max():.3f}, trough {mean_c.min():.3f} (x h)'
        if e8 is not None:
            es = np.interp(pp[o], ge, ye, period=1.0)
            msg += f'; vs TK94: rms difference {np.sqrt(np.mean((mean_c[o] - es)**2)):.4f}, crest {mean_c.max():.3f}/{ye.max():.3f}, trough {mean_c.min():.3f}/{ye.min():.3f} [(zeta-zeta_mean)/h, LES/exp]'
        tab.write(msg + '\n')
    ax8.axhline(0, c='0.6', lw=0.8); ax8.set_xlim(0, 1); ax8.set_ylim(-0.3, 0.4)
    ax8.set_xlabel('t/T'); ax8.set_ylabel(r'$(\zeta-\bar\zeta)/h$'); ax8.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(f'{args.out}_fig8.png', dpi=300, bbox_inches='tight')
    tab.close()
    print('wrote', ', '.join(f'{args.out}_{s}' for s in ('fig4.png', 'fig5.png', 'fig6.png', 'fig7.png', 'fig8.png', 'tables.txt')))


main()
