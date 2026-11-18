#!/usr/bin/env python3
"""Maxey-Riley check (tests/regression/part_mr_neutral): in a forced HIT box a neutrally buoyant particle (rho_p = rho_f, beta = 1)
released at the fluid velocity follows a fluid tracer, dv/dt = Du/Dt + (u-v)/tau. Three runs of the same flow, after 40 steps:
tracer (particle_mode = 0), Maxey-Riley (particle_added_mass = 2) and drag-only inertial (particle_added_mass = 0). Asserts
  1. the Maxey-Riley displacement from the tracer is below --abs-tol and at most --ratio of the drag-only one (the Du/Dt term works),
  2. halving dt reduces the Maxey-Riley displacement (the slip is a time-integration error, not a model error)."""
import argparse, glob, os, re, shutil, subprocess, sys, tempfile

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--case-dir", required=True)
ap.add_argument("--exe", required=True)
ap.add_argument("--mpirun", default="mpirun")
ap.add_argument("--np", type=int, default=2)
ap.add_argument("--abs-tol", type=float, default=3e-3)
ap.add_argument("--ratio", type=float, default=0.6)
a = ap.parse_args()
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../postProcessing"))
from dopamine_post.particles import read_snapshot   # noqa: E402

base = open(os.path.join(a.case_dir, "input_parameters")).read()
PART = {"tracer": "particle_mode = 0,",
        "mr": "particle_mode = 1, particle_diam = 0.05, particle_rho = 1.0, particle_rho_f = 1.0, particle_added_mass = 2, particle_init_fluid = 1,",
        "drag": "particle_mode = 1, particle_diam = 0.05, particle_rho = 1.0, particle_rho_f = 1.0, particle_added_mass = 0, particle_init_fluid = 1,"}


def run(kind, dt, nsteps):
    w = tempfile.mkdtemp(prefix=f"mr_{kind}_")
    s = re.sub(r"particle_mode\s*=\s*0,", PART[kind], base)
    s = re.sub(r"dt\s*=\s*[0-9.eE+-]+,", f"dt = {dt},", s, count=1)
    s = re.sub(r"nsteps\s*=\s*\d+", f"nsteps = {nsteps}", s)
    s = re.sub(r"nsave\s*=\s*\d+", f"nsave = {nsteps}", s)
    open(os.path.join(w, "input_parameters"), "w").write(s)
    r = subprocess.run([a.mpirun, "-np", str(a.np), a.exe], cwd=w, capture_output=True, text=True)
    if r.returncode:
        print(r.stdout[-1500:], r.stderr[-800:]); sys.exit(f"{kind}: solver failed")
    f = glob.glob(os.path.join(w, "fields", "*_particles.*"))[0]
    i, d = read_snapshot(f)
    shutil.rmtree(w)
    return d[np.argsort(i), 0:3]


L = None
m = re.search(r"Lx\s*=\s*([0-9.eE+-]+)", base); nx = int(re.search(r"nx\s*=\s*(\d+)", base).group(1))
P = float(m.group(1)) * (nx - 2) / (nx - 1)


def dist(x, y):
    d = x - y
    d -= P * np.round(d / P)
    return np.sqrt((d ** 2).sum(1).mean())


dt0 = float(re.search(r"dt\s*=\s*([0-9.eE+-]+),", base).group(1))
tr, mr, dr = run("tracer", dt0, 40), run("mr", dt0, 40), run("drag", dt0, 40)
tr2, mr2 = run("tracer", dt0 / 2, 80), run("mr", dt0 / 2, 80)
e_mr, e_dr, e_mr2 = dist(tr, mr), dist(tr, dr), dist(tr2, mr2)
print(f"rms displacement from the tracer: Maxey-Riley {e_mr:.2e}, drag-only {e_dr:.2e}, Maxey-Riley at dt/2 {e_mr2:.2e}")
ok = e_mr < a.abs_tol and e_mr < a.ratio * e_dr and e_mr2 < e_mr
print("PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
