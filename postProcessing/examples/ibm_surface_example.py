from dopamine_post.ibm_surface import IBMSurface

case = "../../examples/dns_ibm_wavyWall"

surf = IBMSurface.read(f"{case}/ibm_surface/surface.00010000.bin")
Fp, Fv = surf.forces()
print(f"t={surf.t:.6g}  n={surf.n}  Fpres={Fp}  Fvisc={Fv}")

surf.to_vtp("surface_00010000.vtp")

IBMSurface.to_pvd(f"{case}/ibm_surface/surface.*.bin", "surface.pvd")
