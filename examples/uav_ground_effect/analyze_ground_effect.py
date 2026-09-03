#!/usr/bin/env python3
"""Compare the IGE (near_ground) and OGE (far_ground) vertical-velocity
line probes recorded through the disk centre, after both cases have been
run (see each case's input_parameters for setup). Run from this directory:

    python3 analyze_ground_effect.py

Expects near_ground/lines/centerline.bin(+_meta.txt) and the far_ground
equivalent to already exist (i.e. both cases have been run).
"""
import numpy as np
import os

def load_line(base):
    meta = {}
    with open(base + '_meta.txt') as f:
        for line in f:
            if '=' in line:
                k, v = line.split('=', 1)
                meta[k.strip()] = v.strip()
    npts   = int(meta['npts'])
    nsnaps = int(meta['nsnaps'])
    # big-endian float64, no record markers (see docs/Input-Parameters.md &STATISTICS)
    d = np.fromfile(base + '.bin', dtype='>f8')
    return d.reshape(nsnaps, 1, npts)[:, 0, :]

def main():
    here = os.path.dirname(os.path.abspath(__file__))
    near = load_line(os.path.join(here, 'near_ground/lines/centerline'))
    far  = load_line(os.path.join(here, 'far_ground/lines/centerline'))

    Ly, npts = 2.0, near.shape[1]
    ym = (np.arange(npts) + 0.5) * Ly / npts   # cell centres, uniform y grid
    h_near, h_far = 0.3, 1.2                    # must match each case's uav_yc
    offset = 0.03                                # sample this far below the disk

    def just_below(prof_series, h):
        j = np.argmin(np.abs(ym - (h - offset)))
        return prof_series[:, j]

    vb_near, vb_far = just_below(near, h_near), just_below(far, h_far)
    n_avg = min(10, near.shape[0], far.shape[0])

    print(f"V (downward, m/s) at {offset} m below the disk, last {n_avg} snapshots:")
    print(f"  near_ground (h/R=2): {vb_near[-n_avg:].mean():.3f} +/- {vb_near[-n_avg:].std():.3f}")
    print(f"  far_ground  (h/R=8): {vb_far[-n_avg:].mean():.3f} +/- {vb_far[-n_avg:].std():.3f}")
    ratio = vb_near[-n_avg:].mean() / vb_far[-n_avg:].mean()
    print(f"  ratio (near/far)   : {ratio:.3f}  (ground effect at constant thrust predicts <1)")

if __name__ == '__main__':
    main()
