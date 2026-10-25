"""Shared wave-gauge statistics for the TK94 post-processing (postprocess.py, compare_simulation.py)."""
import re
import numpy as np


def namelist_values(text, key):
    """numbers following 'key =' in a Fortran namelist (comments after '!' stripped; value list ends at the next 'name =' or '/')"""
    text = re.sub(r'!.*', '', text)
    m = re.search(r'\b%s\s*=\s*' % re.escape(key), text)
    if not m:
        return None
    rest = text[m.end():]
    end = re.search(r'[A-Za-z_]\w*\s*=|/', rest)
    return [float(v) for v in (rest[:end.start()] if end else rest).replace(',', ' ').split()]


def resolve_restarts(t):
    """indices of the rows kept after a restart: rows after a step back in time are stale and dropped (equal times are kept)"""
    kept = []
    for i, ti in enumerate(t):
        if kept and ti < t[kept[-1]]:
            kept = kept[:np.searchsorted(t[kept], ti, side='left')]
        kept.append(i)
    return np.array(kept, dtype=int)


def phase_stats(t, eta, t0, t_last, T, nphase=200):
    """Statistics over the whole periods in [t0, t0 + nper*T], nper = floor((t_last - t0)/T). t may have a variable time step.
    * mean level: time-weighted (trapezoid) mean over the window, split exactly at the period boundaries
    * phase average: eta linearly resampled at nphase equidistant phases of every period and averaged over the periods;
      crest / trough = max / min of that phase-averaged wave, H = crest - trough
    * per-period values (same resampling) give the cycle-to-cycle scatter.
    returns dict(nper, t_end, mean, crest, trough, H, phase_avg, period_mean, period_crest, period_trough)"""
    nper = int(np.floor((t_last - t0)/T + 1e-9))
    if nper < 1:
        raise SystemExit(f'no complete {T:g} s wave period between t0={t0:g} s and t_end={t_last:g} s')
    t_end = t0 + nper*T
    if t[0] > t0 + 1e-9 or t[-1] < t_end - 1e-9:
        raise SystemExit(f'gauge record {t[0]:g}..{t[-1]:g} s does not cover the window {t0:g}..{t_end:g} s')
    pmean, resampled = [], []
    for p in range(nper):
        a, b = t0 + p*T, t0 + (p + 1)*T
        inside = (t > a) & (t < b)
        tt = np.concatenate(([a], t[inside], [b]))
        ee = np.interp(tt, t, eta)
        pmean.append(np.sum(0.5*(ee[1:] + ee[:-1])*np.diff(tt))/T)
        resampled.append(np.interp(a + T*np.arange(nphase)/nphase, t, eta))
    resampled = np.array(resampled)
    pa = resampled.mean(axis=0)
    return dict(nper=nper, t_end=t_end, mean=float(np.mean(pmean)), crest=pa.max(), trough=pa.min(), H=pa.max() - pa.min(),
                phase_avg=pa, period_mean=np.array(pmean), period_crest=resampled.max(axis=1), period_trough=resampled.min(axis=1))
