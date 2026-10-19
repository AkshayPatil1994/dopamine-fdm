from dopamine_post.particles import ParticleData, write_particles_xmf

case = "../../examples/particle_channel_turbophoresis"

pdat = ParticleData.from_case(case)
print(pdat)

# ParaView: one point cloud per saved step (id, age, Velocity, Acceleration), time = solver step number
write_particles_xmf(case)

# same thing as a method, writing to a custom path
pdat.write_xmf(out=f"{case}/paraview/particles_copy.xmf")
