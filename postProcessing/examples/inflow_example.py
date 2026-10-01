from dopamine_post.inflow import mirror_half_channel, InflowDonor, InflowOptState

case = "../../examples/precursor_successor/precursor"

mirror_half_channel(f"{case}/half_channel.csv", "full_channel.csv")

donor = InflowDonor.read(f"{case}/inflow_data/inflow_planes")
mean, rey = donor.time_span_stats()
donor.plot(reference=f"{case}/full_inflow.csv", grid=f"{case}/fields/grid.out",
           output="inflow_donor_check.png")

state = InflowOptState.read(f"{case}/fields/inflow_opt_data.dat")
state.summary()
