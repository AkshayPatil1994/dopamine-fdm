#!/usr/bin/env python3
"""Physics check of the HIT forcing and IC (hit_stats.dat of tests/regression/hit_budget):
  1. the IC has the requested TKE 1.5*Utarget^2 within --ic-tol (sampling of a few modes),
  2. the TKE budget closes: dK/dt = eps_in - eps (central differences of the monitor rows) within --budget-tol,
     which also fixes the periodic length of all three directions (a wrong period breaks the energy conservation of
     the skew-symmetric convection),
  3. the discrete divergence stays at round-off (max_divergence column of the monitor)."""
import argparse, os, re, shutil, subprocess, sys, tempfile

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--case-dir", required=True)
ap.add_argument("--exe", required=True)
ap.add_argument("--mpirun", default="mpirun")
ap.add_argument("--np", type=int, default=4)
ap.add_argument("--ic-tol", type=float, default=0.10)
ap.add_argument("--budget-tol", type=float, default=0.02)
a = ap.parse_args()

work = tempfile.mkdtemp(prefix="hit_check_")
shutil.copy(os.path.join(a.case_dir, "input_parameters"), work)
txt = open(os.path.join(work, "input_parameters")).read()
eps_in = float(re.search(r"hit_eps\s*=\s*([0-9.eE+-]+)", txt).group(1))
u0 = float(re.search(r"Utarget\s*=\s*([0-9.eE+-]+)", txt).group(1))
r = subprocess.run([a.mpirun, "-np", str(a.np), a.exe], cwd=work, capture_output=True, text=True)
if r.returncode != 0:
    print(r.stdout[-2000:], r.stderr[-2000:]); sys.exit("solver failed")
d = np.loadtxt(os.path.join(work, "hit_stats.dat"))
t, tke, eps = d[:, 1], d[:, 2], d[:, 3]
ok = True
e0 = abs(tke[0]/(1.5*u0**2) - 1)
print(f"IC TKE / (1.5 u0^2) - 1 = {e0:.3f}  (tol {a.ic_tol})"); ok &= e0 < a.ic_tol
dk = (tke[2:] - tke[:-2])/(t[2:] - t[:-2])
res = np.abs(dk - (eps_in - eps[1:-1]))/eps_in
print(f"budget |dK/dt - (eps_in - eps)|/eps_in: max {res.max():.3f}  (tol {a.budget_tol})"); ok &= res.max() < a.budget_tol
divs = [float(l.split()[4]) for l in r.stdout.splitlines() if re.match(r"^\s*\d+\s+[0-9.]+E[+-]\d+\s", l)]
print(f"max divergence {max(divs):.2e}"); ok &= max(divs) < 1e-9
shutil.rmtree(work)
sys.exit(0 if ok else 1)
