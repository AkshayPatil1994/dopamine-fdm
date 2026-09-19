# Showcase Animations

The animations on the README are rendered from the validation-run snapshots by the scripts in
`postProcessing/animations/`. They read single planes straight out of the snapshot files
(memory-mapped, so a 4 GiB 512³ snapshot costs a few MiB of I/O per frame), render frames in
parallel with matplotlib, and encode an optimised GIF with `ffmpeg`.

| Script | Case | Content |
|--------|------|---------|
| `animate_tgv.py` | Taylor–Green vortex, Re = 1600, 512³ | vorticity magnitude on a `z = const` plane (central differences from three adjacent planes) |
| `animate_rb.py` | Rayleigh–Bénard convection (Boussinesq) | temperature: vertical slice and near-wall plan view |
| `animate_wavywall.py` | Wavy-wall DNS (ghost-cell IBM) | streamwise velocity, vertical slice with the solid masked from the SDF, and plan view |
| `animate_chan395.py` | Wall-modelled LES, Re_τ = 395 | streamwise velocity slice and near-wall streaks `u'` in plan view |

`anim_common.py` holds the shared code: the `Snapshot` plane reader, the step→time map
(interpolated from the solver's monitor lines in `run.log`, since snapshots carry no time stamp),
and the GIF encoder.

```bash
python postProcessing/animations/animate_tgv.py \
    --case-dir /path/to/tgv/512 --every 3 --out docs/animations/tgv.gif
```

Common options: `--case-dir`, `--prefix` (the `fileout` of the run), `--out`, `--every N`
(use every N-th snapshot), `--fps`. Case-specific: `--z-frac` (TGV plane), `--y-plan` (RB, wavy
wall plan-view height), `--yplus` and `--nstart` (channel). Requires `numpy`, `matplotlib`, and
`ffmpeg` on the `PATH`.
