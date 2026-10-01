"""Example: render a uav_hover_disk case's UAV path for ParaView.

Target case: ../../examples/uav_hover_disk (a solver run produces
uav_path_file and input_parameters there; this script has nothing to run
against until that case has actually been simulated in this checkout).
"""
from pathlib import Path

from dopamine_post.uav import UAVPath

case = Path(__file__).resolve().parents[2] / "examples" / "uav_hover_disk"

path = UAVPath.read(case / "uav_path_file.dat")
print(path)

centre_mid = path.center(0.5 * (path.t[0] + path.t[-1]))
print("centre at path midpoint time:", centre_mid)

path.disk_animation(input_parameters=case / "input_parameters", nframes=60,
                     out_dir=case / "uav_path", pvd_out=case / "uav_path.pvd")

path.drone_animation(input_parameters=case / "input_parameters", nframes=60, tilt=True,
                      stl_out=case / "uav_drone.stl",
                      transforms_out=case / "uav_transforms",
                      pv_source_out=case / "uav_drone_paraview_source.py")
