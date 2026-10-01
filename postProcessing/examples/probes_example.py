from dopamine_post.probes import LineProbe, SliceProbe

case = "../../examples/les_chan395"

line = LineProbe.read(f"{case}/line_meta.txt")
u = line.get("U")
print(f"line probe: dir={line.dir}  comps={line.comps}  shape={u.shape}")

slc = SliceProbe.read(f"{case}/slice_meta.txt")
w = slc.get("W")
print(f"slice probe: dir={slc.dir}  comps={slc.comps}  shape={w.shape}")

slc.write_xmf(snap_path=f"{case}/fields/grid.out")
