"""Particle snapshots (<prefix>_particles.<step>): read them as arrays and write the point-cloud XMF.
Run:   python particles_example.py [case_dir]        (default: examples/particle_channel_turbophoresis)
For particles together with the flow field in ParaView see xmf_particles.py."""
from _case import case_dir
from dopamine_post.particles import ParticleData, write_particles_xmf

case = case_dir("particle_channel_turbophoresis")

pdat = ParticleData.from_case(case)
print(pdat)

# ParaView: one point cloud per saved step (id, age, Velocity, Acceleration), time = solver step number
write_particles_xmf(case)

# same thing as a method, writing to a custom path
pdat.write_xmf(out=f"{case}/paraview/particles_copy.xmf")
