#!/usr/bin/env python3
"""Donor-slice column check (tests/regression/sem_recycle_p): an inflow_type=2 donor whose comps are UVWPT (P=99, T=7, interleaved
as probe_output.f90 writes them) must feed T=7 at the inlet, not the P column."""
import argparse, glob, os, shutil, subprocess, sys, tempfile

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--case-dir", required=True)
ap.add_argument("--exe", required=True)
ap.add_argument("--mpirun", default="mpirun")
a = ap.parse_args()
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../postProcessing"))
from dopamine_post.fields import FieldSnapshot   # noqa: E402

n1, n2 = 16, 17    # nym_global, nzm_global of this case
w = tempfile.mkdtemp(prefix="recycle_p_")
shutil.copy(os.path.join(a.case_dir, "input_parameters"), w)
shutil.copy(os.path.join(a.case_dir, "inflow_profile.dat"), w)
frame = np.zeros((n1 * n2, 5)); frame[:, 0] = 5.0; frame[:, 1:3] = 0.01; frame[:, 3] = 99.0; frame[:, 4] = 7.0
np.concatenate([frame, frame]).astype(">f8").tofile(os.path.join(w, "donor.bin"))   # the solver reads big-endian
np.array([0.0, 1000.0]).astype(">f8").tofile(os.path.join(w, "donor_times.bin"))
open(os.path.join(w, "donor_meta.txt"), "w").write(
    f"ncomp  = 5\nn1     = {n1}\nn2     = {n2}\ndir    = x\npos    = 1.0\ncomps  = UVWPT\nnsnaps = 2\ntimes  = donor_times.bin\n")
r = subprocess.run([a.mpirun, "-np", "1", a.exe], cwd=w, capture_output=True, text=True)
if r.returncode:
    print(r.stdout[-1500:], r.stderr[-800:]); sys.exit("solver failed")
snap = FieldSnapshot.read(glob.glob(os.path.join(w, "fields", "*.20"))[0], sgs_model=1, boussinesq_flag=1,
                          sediment_flag=0, vof_active=0)
T = snap["T"][0]
shutil.rmtree(w)
print(f"inlet-plane T in [{T.min():.3f}, {T.max():.3f}] (donor T=7, donor P=99)")
sys.exit(0 if abs(T.mean() - 7.0) < 0.5 else "FAIL: inlet T is not the donor T column")
