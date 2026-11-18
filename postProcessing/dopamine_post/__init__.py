"""
dopamine_post -- post-processing library for the dopamine-fdm CFD solver.

Import once, call methods -- no more hunting through independent scripts:

    import dopamine_post as dp

    pdat = dp.particles.ParticleData.from_case(".")
    snap = dp.fields.FieldSnapshot.read("fields/channel_test.100000")

Submodules (also reachable directly, e.g. `from dopamine_post.particles import ParticleData`):
    fields        binary field snapshots + channel statistics
    particles     Lagrangian point particles
    probes        line + slice probes
    ibm_surface   IBM surface samples
    sdf           signed-distance-field reading/plotting/meshing
    uav           UAV actuator-disk path + drone geometry
    vortex        Q-criterion + iso-surface movies
    rsb           Reynolds-stress budget
    tavg          on-the-fly time-averaged statistics
    inflow        SEM/ESEM inflow tooling
    runlog        run.log diagnostics
    diagnostics   plain-text series (vof_diag, wall_diag, ibm_forces, particle_count, gauges, wave check)
    hit           forced homogeneous isotropic turbulence statistics and spectra
"""
from . import (
    diagnostics,
    fields,
    hit,
    ibm_surface,
    inflow,
    particles,
    probes,
    rsb,
    runlog,
    sdf,
    tavg,
    uav,
    vortex,
)

__all__ = [
    "fields", "particles", "probes", "ibm_surface", "sdf", "uav", "rsb",
    "inflow", "runlog", "vortex", "tavg", "diagnostics", "hit",
]
