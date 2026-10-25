"""Particles together with the flow field in ParaView (single-phase or VOF carrier flow with particles_active = 1).

Run:   python xmf_particles.py [case_dir]        (default: examples/particle_channel_turbophoresis)
Writes in <case>/paraview/: channel_test.xmf (fields, from fields/<fileout>.<step>) and particles.xmf (point clouds, from
fields/<fileout>_particles.<step>). Open both files; the time value of each frame is the solver step in both, so they
play in sync. Particle attributes: id, age, Velocity (and Acceleration for the 10-column files).
Keep the solver's stdout as <case>/run.log (`mpirun ... dopamine > run.log`) so ParticleData gets the physical time of each step;
without it t = step*dt is used and a warning says so.
CLI equivalent:  dopamine-post fields xmf --case-dir <case>;  dopamine-post particles xmf --case-dir <case>
"""
from _case import case_dir
from dopamine_post.fields import write_field_xmf
from dopamine_post.particles import ParticleData, write_particles_xmf

case = case_dir("particle_channel_turbophoresis")
write_field_xmf(case)
write_particles_xmf(case)

pdat = ParticleData.from_case(case)          # the same data as arrays: pos / vel / acc of shape (Nt, Np, 3), times pdat.t
print(pdat)
