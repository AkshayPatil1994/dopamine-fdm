from dopamine_post.sdf import SDF

case = "../../examples/dns_ibm_wavyWall"

s = SDF.read(case, sdf_file="geo/SDF_in", params="input_parameters")

fig1 = s.plot_xz(y=0.05)
fig1.savefig("sdf_xz.png", dpi=150)

fig2 = s.plot_xy(z=s.zm[len(s.zm) // 2])
fig2.savefig("sdf_xy.png", dpi=150)

y_actual, data = s.slice_xz(y=0.1)
print(f"X-Z slice at y={y_actual:.4f}: phi range [{data.min():.4f}, {data.max():.4f}]")
