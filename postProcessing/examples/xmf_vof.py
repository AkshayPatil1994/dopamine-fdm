"""ParaView time series of a two-fluid VOF case: velocity, pressure and the liquid fraction `vof` (0 = gas, 1 = liquid).

Run any VOF case with nsave > 0 first (input_parameters with vof_active = 1), then:
       python xmf_vof.py <case_dir>
Writes <case>/paraview/channel_test.xmf; in ParaView colour by `vof` and add a Contour at 0.5 for the interface.
The block `vof` is the last one in the snapshot (after C, nu_t, T), which write_field_xmf takes from input_parameters.
CLI equivalent:  dopamine-post fields xmf --case-dir <case>
The vof_debug = 1 snapshots (vof_snap_*.dat) are plain text and are plotted by the scripts in examples/vof_zalesak_disk and
postProcessing/animations, not through XMF.
"""
import sys
from dopamine_post.fields import FieldSeries, write_field_xmf

if len(sys.argv) < 2:
    sys.exit(__doc__)
case = sys.argv[1]
write_field_xmf(case)
print("vof range of the latest snapshot:", FieldSeries(case).latest()["vof"].min(), FieldSeries(case).latest()["vof"].max())
