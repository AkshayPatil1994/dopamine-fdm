"""
tavg.py -- on-the-fly time-averaged statistics written by src/time_stats.f90.

    from dopamine_post.tavg import read_tavg

    t_avg, f = read_tavg("fields/chan_tavg.20000")
    f["U"], f["u'u'"], f["u'w'"]        # (nxm, nym, nzm) cell-centre arrays
    urms = f["u'u'"] ** 0.5

Fields: U V W, the Reynolds stresses u'u' v'v' w'w' u'v' u'w' v'w', and per enabled scalar S (T or C)
S, S'S', u'S', v'S', w'S'. `t_avg` is the averaging time of the file.
"""
import numpy as np


def read_tavg(path):
    """Return (averaging_time, {name: array(nxm, nym, nzm)}) from a fields/<fileout>_tavg.<step> file."""
    with open(path, "rb") as fh:
        nx, ny, nz, nf = np.fromfile(fh, "<i4", 4)
        t_avg = float(np.fromfile(fh, "<f8", 1)[0])
        names = [fh.read(8).decode().strip() for _ in range(nf)]
        data = np.fromfile(fh, "<f8", nf * nx * ny * nz).reshape((nf, nz, ny, nx))
    return t_avg, {n: data[i].transpose(2, 1, 0) for i, n in enumerate(names)}
