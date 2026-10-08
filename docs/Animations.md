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
| `animate_vof_breaker.py` | Plunging breaker, two-fluid VOF (ak = 0.55, 1000:1, 128² cells, t ≤ 0.8 s) | liquid fraction and interface from the `vof_snap_*.dat` debug snapshots (`vof_debug = 1`, `vof_snap_dt`; set `nmonitor` so a snapshot is written every ~0.02 s, since snapshots are written at the monitor interval) |
| `animate_vof_zalesak.py` | Zalesak's slotted disk (`examples/vof_zalesak_disk`, 100² cells) | liquid fraction with the exact rotated shape as a dashed outline |

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

The two VOF scripts read the debug snapshots of a `vof_debug = 1` run (`--case-dir` = the run directory, `--tmax`, `--every`) and
draw the liquid fraction with bicubic interpolation at 256 GIF colours. The breaker is the `vof_break_small` input at N = 128 with
`vof_sigma = 0.072`, `vof_wave_amp = 0.0875`, `vof_snap_dt = 0.02`, `nmonitor = 4` and `sim_end_time = 1` (about 15 min on 4 ranks; the
animation stops at 0.8 s, after the splash); the Zalesak run takes about 30 s on one rank with `vof_snap_dt = 0.02`, `nmonitor = 20`.
