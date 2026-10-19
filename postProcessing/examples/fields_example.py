from dopamine_post.fields import FieldSeries, plot_profile, plot_stats, plot_profile_vs_dns

case = "../../examples/les_chan395"

series = FieldSeries(case)

snap = series.latest()
ym, U = snap.profile("U")
fig = plot_profile(snap, save="profile_latest.png")

stats = series.time_average((100000, 150000), 500)
stats.save("stats")
plot_stats(stats, case_dir=case, save="stats.png")

result = series.profile_vs_dns(dns_dir=f"{case}/dns_ref", x_stations=[2, 4, 6, 8])
plot_profile_vs_dns(result, save="channel_profiles.png")

series.write_xmf()

# chunked_snapshot=.true. output (<fileout>Block.<step>): read one box and write its XDMF time series
from dopamine_post.fields import read_box_snapshot, write_box_xmf
box = read_box_snapshot(f"{case}/fields/channel_testBlock.1000")
print(box["f1"], box["f2"], {k: v[0].shape for k, v in box["fields"].items()})
write_box_xmf(case)
