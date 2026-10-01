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
    rsb           Reynolds-stress budget
    inflow        SEM/ESEM inflow tooling
    runlog        run.log diagnostics
"""
from . import (
    fields,
    ibm_surface,
    inflow,
    particles,
    probes,
    rsb,
    runlog,
    sdf,
    uav,
)

__all__ = [
    "fields", "particles", "probes", "ibm_surface", "sdf", "uav", "rsb",
    "inflow", "runlog",
]
