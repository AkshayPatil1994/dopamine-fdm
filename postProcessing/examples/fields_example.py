"""Overview of the field readers; each use case has its own script (see README.md in this directory):
fields_single_snapshot.py, fields_multiple_snapshots.py, xmf_single_phase.py, xmf_vof.py, xmf_particles.py, xmf_box.py.
This one compares the LES channel profiles with the DNS reference and needs examples/les_chan395/dns_ref."""
from _case import case_dir
from dopamine_post.fields import FieldSeries, plot_profile_vs_dns

case = case_dir("les_chan395")
result = FieldSeries(case).profile_vs_dns(dns_dir=f"{case}/dns_ref", x_stations=[2, 4, 6, 8])
plot_profile_vs_dns(result, save="channel_profiles.png")
