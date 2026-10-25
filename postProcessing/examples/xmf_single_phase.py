"""ParaView time series of the full snapshots of a single-phase case (also Boussinesq T, sediment C, SGS nu_t).

Run:   python xmf_single_phase.py [case_dir]        (default: examples/les_chan395)
Writes <case>/paraview/channel_test.xmf (+ xm/ym/zm.bin); open the .xmf in ParaView ("Xdmf3 Reader" or "XDMF Reader").
The .xmf points at the solver's own binary files (no data is copied), so keep fields/ where it is.
CLI equivalent:  dopamine-post fields xmf --case-dir <case>
"""
from _case import case_dir
from dopamine_post.fields import FieldSeries, write_field_xmf

case = case_dir("les_chan395")
write_field_xmf(case)                      # the attributes (Velocity, P, nu_t, T, C, vof) follow input_parameters
# FieldSeries(case).write_xmf()            # same, as a method
