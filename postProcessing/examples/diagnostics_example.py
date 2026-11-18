"""Example: the solver's plain-text time series (vof_diag.dat, wall_diag.dat, ibm_forces.csv, particle_count.dat, vof_gauges.dat).

Usage: python diagnostics_example.py path/to/case
"""
import sys
from pathlib import Path

from dopamine_post import diagnostics as dg

case = Path(sys.argv[1])

if (case / "vof_diag.dat").is_file():
    d = dg.read_vof_diag(case / "vof_diag.dat")
    print("liquid volume drift:", d["dVliq_rel"][-1], " interface cells:", d["n_interface"][-1])
    dg.plot_vof_diag(d, "vof_diag.png")
if (case / "wall_diag.dat").is_file():
    d = dg.read_wall_diag(case / "wall_diag.dat")
    print("u_tau at the low wall (last row):", d["u_tau_lo"][-1])
if (case / "ibm_forces.csv").is_file():
    d = dg.read_ibm_forces(case / "ibm_forces.csv")
    print("pressure drag Fx_pres (last row):", d["Fx_pres"][-1])
if (case / "fields" / "particle_count.dat").is_file():
    d = dg.read_particle_counts(case / "fields" / "particle_count.dat")
    print("active particles:", d["active"][-1])
if (case / "vof_gauges.dat").is_file():
    d = dg.read_vof_gauges(case / "vof_gauges.dat")
    print("gauges at x =", d["x"], " eta shape", d["eta"].shape)
