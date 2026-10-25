"""Read ONE full snapshot (fields/<fileout>.<step>) of any case type and look at it.

Run:   python fields_single_snapshot.py [case_dir]        (default: examples/les_chan395)

Which optional blocks the file holds (C, nu_t, T, vof, in this order after U,V,W,P) comes from <case>/input_parameters:
sediment_flag, sgs_model (missing = solver default AMD, so nu_t exists), boussinesq_flag, vof_active.
"""
from _case import case_dir
from dopamine_post.fields import FieldSeries, FieldSnapshot

case = case_dir("les_chan395")

# 1. the flags are read from <case>/input_parameters: single-phase, Boussinesq, sediment and VOF cases all use this call
series = FieldSeries(case)
snap = series.latest()                      # or series.read(step)
print(snap)                                 # step, grid and the fields present
print("fields:", sorted(snap.fields))       # e.g. U V W P nu_t [C] [T] [vof]; all cell-centred, ghost layers removed, shape (nx, ny, nz)

ym, U = snap.profile("U")                   # x-z averaged wall-normal profile
plane = snap.plane("U", axis="z", index=snap.fields["U"].shape[2] // 2)   # an x-y plane

# 2. no input_parameters next to the file (copied snapshot) or a deliberate override: give the flags of the run that wrote it.
#    Wrong flags raise ValueError (the number of blocks does not match) instead of silently shifting the fields.
snap2 = FieldSnapshot.read(snap.path, sgs_model=2, sediment_flag=0, boussinesq_flag=0, vof_active=0)

# 3. a chunked_snapshot=.true. case writes <fileout>Block.<step>: the box fields are averaged to the cell centres
#    box = FieldSnapshot.from_box("fields/<fileout>Block.1000")       # same profile()/plane() interface
