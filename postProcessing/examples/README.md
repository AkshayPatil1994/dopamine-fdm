# dopamine_post examples

Every script takes the case directory as its first argument (`python xmf_vof.py path/to/case`). The case must have been run
first: the snapshots are written by the solver (`nsave > 0`, `fileout = '<prefix>'` in `input_parameters`). Install once with
`pip install -e postProcessing` (or set `PYTHONPATH=postProcessing`). For particle runs keep the solver output as
`<case>/run.log` so step numbers are converted to physical time.

| What you ran | Read it | Open in ParaView |
|---|---|---|
| single-phase LES/DNS, `fields/<prefix>.<step>` | `fields_single_snapshot.py`, `fields_multiple_snapshots.py` | `xmf_single_phase.py` |
| + `boussinesq_flag`, `sediment_flag`, `sgs_model` | same scripts; the extra fields `T`, `C`, `nu_t` appear automatically | same script |
| two-fluid VOF (`vof_active = 1`) | same scripts; field `vof` | `xmf_vof.py` |
| particles (`particles_active = 1`, `<prefix>_particles.<step>`) | `particles_example.py` | `xmf_particles.py` (fields + particles) |
| `chunked_snapshot = .true.` (`<prefix>Block.<step>`) | `FieldSnapshot.from_box`, see `fields_single_snapshot.py` | `xmf_box.py` |
| slice / line probes | `probes_example.py` | `SliceProbe.write_xmf` (same script) |
| IBM surface, SDF, UAV path, RSB, SEM inflow, run log, Q-criterion | `ibm_surface_example.py`, `sdf_example.py`, `uav_example.py`, `rsb_example.py`, `inflow_example.py`, `runlog_example.py`, `vortex_example.py` | `ibm_surface_example.py` writes `.vtp/.pvd` |

Which optional blocks a full snapshot contains (after U, V, W, P: `C`, `nu_t`, `T`, `vof`) is taken from `<case>/input_parameters`
(`sediment_flag`, `sgs_model`, `boussinesq_flag`, `vof_active`; a missing `sgs_model` means the default AMD, so `nu_t` is present).
If the file holds a different number of blocks `FieldSnapshot.read` raises a `ValueError`; pass the flags explicitly then.
The `.xmf` files only point at the solver's binary files, so do not move or delete `fields/`.
