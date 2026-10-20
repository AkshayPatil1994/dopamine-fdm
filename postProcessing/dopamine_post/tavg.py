"""
tavg.py -- on-the-fly time-averaged statistics written by src/time_stats.f90.

    from dopamine_post.tavg import read_tavg

    t_avg, f = read_tavg("fields/chan_tavg.20000")
    f["U"], f["u'u'"], f["u'w'"]        # (nxm, nym, nzm) cell-centre arrays
    urms = f["u'u'"] ** 0.5

Fields: U V W, the Reynolds stresses u'u' v'v' w'w' u'v' u'w' v'w', and per enabled scalar S (T or C)
S, S'S', u'S', v'S', w'S'. `t_avg` is the averaging time of the file.

A chunk-box file (`<fileout>_tavgBlock.<step>`, written with `chunked_snapshot`) has the same layout for the box and a
trailer with its cell-centre coordinates; `read_tavg(path, return_coords=True)` returns them as a third value
`{"f1", "f2", "xm", "ym", "zm"}` (None for a full-domain file).
"""
import numpy as np


def read_tavg(path, return_coords=False):
    """Return (averaging_time, {name: array(nx, ny, nz)}) from a fields/<fileout>_tavg[Block].<step> file."""
    with open(path, "rb") as fh:
        nx, ny, nz, nf = np.fromfile(fh, "<i4", 4)
        t_avg = float(np.fromfile(fh, "<f8", 1)[0])
        names = [fh.read(8).decode().strip() for _ in range(nf)]
        data = np.fromfile(fh, "<f8", nf * nx * ny * nz).reshape((nf, nz, ny, nx))
        coords = None
        if fh.read(1):
            fh.seek(-1, 1)
            f1, f2 = np.fromfile(fh, "<i4", 3), np.fromfile(fh, "<i4", 3)
            xm, ym, zm = (np.fromfile(fh, "<f8", n) for n in (nx, ny, nz))
            coords = {"f1": f1, "f2": f2, "xm": xm, "ym": ym, "zm": zm}
    fields = {n: data[i].transpose(2, 1, 0) for i, n in enumerate(names)}
    return (t_avg, fields, coords) if return_coords else (t_avg, fields)
