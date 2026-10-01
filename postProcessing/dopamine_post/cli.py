"""
cli.py -- the `dopamine-post` command-line entry point.

Thin argument parsing only: every subcommand here just unpacks flags and calls
straight into the corresponding dopamine_post module. The library is fully
usable without this file -- see postProcessing/examples/ for the same
operations driven by `import dopamine_post as dp` instead.
"""
import argparse
import sys

from . import fields, ibm_surface, inflow, particles, probes, rsb, runlog, sdf, uav, vortex


def _slice_arg(last):
    return slice(-last, None) if last else slice(None)


# ── particles ──────────────────────────────────────────────────────────────

def _particles_pdf(a):
    pdat = particles.ParticleData.from_case(a.case_dir, a.prefix, every=a.every)
    particles.plot_pdfs(pdat, a.quantities, a.out or "pdf.png", _slice_arg(a.last),
                        log=a.log, bins=a.bins)


def _particles_tracks(a):
    pdat = particles.ParticleData.from_case(a.case_dir, a.prefix, every=a.every)
    frames = _slice_arg(a.last)
    first = range(len(pdat.t))[frames][0]
    ids = (particles.np.array(a.ids) if a.ids else
           pdat.pick_ids(a.n, a.seed, first))
    geo = particles.Geometry.flat(pdat.L)
    particles.plot_tracks(pdat, ids, a.out or "tracks.png", frames, geo, a.view)


def _particles_animate(a):
    pdat = particles.ParticleData.from_case(a.case_dir, a.prefix, every=a.every)
    frames = _slice_arg(a.last)
    first = range(len(pdat.t))[frames][0]
    ids = (particles.np.array(a.ids) if a.ids else
           pdat.pick_ids(a.n, a.seed, first))
    geo = particles.Geometry.flat(pdat.L)
    particles.animate_tracks(pdat, ids, a.out or "tracks.gif", frames, geo, a.view,
                             tail=a.tail, fps=a.fps)


def _particles_cloud(a):
    pdat = particles.ParticleData.from_case(a.case_dir, a.prefix, every=a.every)
    geo = particles.Geometry.flat(pdat.L)
    particles.animate_cloud(pdat, a.out or "cloud.gif", _slice_arg(a.last), geo,
                            color=a.color, fps=a.fps)


def _particles_xmf(a):
    particles.write_particles_xmf(a.case_dir, a.prefix, a.out)


def _add_particles(sub):
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--case-dir", default=".")
    common.add_argument("--prefix", default=None)
    common.add_argument("--every", type=int, default=1)
    common.add_argument("--last", type=int, default=20)
    common.add_argument("--out", default=None)

    sel = argparse.ArgumentParser(add_help=False)
    sel.add_argument("--ids", type=int, nargs="+", default=None)
    sel.add_argument("--n", type=int, default=8)
    sel.add_argument("--seed", type=int, default=1)

    p = sub.add_parser("pdf", parents=[common], help="PDFs of quantities pooled over the last N snapshots")
    p.add_argument("quantities", nargs="+")
    p.add_argument("--bins", type=int, default=80)
    p.add_argument("--log", action="store_true")
    p.set_defaults(func=_particles_pdf)

    p = sub.add_parser("tracks", parents=[common, sel], help="static track figure")
    p.add_argument("--view", choices=["2d", "3d"], default="2d")
    p.set_defaults(func=_particles_tracks)

    p = sub.add_parser("animate", parents=[common, sel], help="track movie")
    p.add_argument("--view", choices=["2d", "3d"], default="3d")
    p.add_argument("--tail", type=int, default=0)
    p.add_argument("--fps", type=float, default=8)
    p.set_defaults(func=_particles_animate)

    p = sub.add_parser("cloud", parents=[common], help="all particles coloured by a quantity")
    p.add_argument("--color", choices=["speed", "height", "id"], default="speed")
    p.add_argument("--fps", type=float, default=8)
    p.set_defaults(func=_particles_cloud)

    p = sub.add_parser("xmf", parents=[common], help="write a ParaView XDMF time series")
    p.set_defaults(func=_particles_xmf)


# ── fields ─────────────────────────────────────────────────────────────────

def _fields_profile(a):
    series = fields.FieldSeries(a.case_dir)
    snap = series.read(a.step) if a.step is not None else series.latest()
    fields.plot_profile(snap, save=a.out or "profile.png")


def _fields_stats(a):
    series = fields.FieldSeries(a.case_dir)
    stats = series.time_average((a.start, a.end), a.interval)
    stats.save(a.stats_dir)
    fields.plot_stats(stats, a.case_dir, dns_dir=a.dns_dir, save=a.out or "plots/stats.png")


def _fields_xmf(a):
    fields.write_field_xmf(a.case_dir, a.prefix, a.out)


def _add_fields(sub):
    p = sub.add_parser("profile", help="wall-normal profile plot from one snapshot")
    p.add_argument("--case-dir", default=".")
    p.add_argument("--step", type=int, default=None, help="default: latest snapshot")
    p.add_argument("--out", default=None)
    p.set_defaults(func=_fields_profile)

    p = sub.add_parser("stats", help="time/plane-averaged channel statistics")
    p.add_argument("--case-dir", default=".")
    p.add_argument("--start", type=int, required=True)
    p.add_argument("--end", type=int, required=True)
    p.add_argument("--interval", type=int, default=1)
    p.add_argument("--stats-dir", default="stats")
    p.add_argument("--dns-dir", default=None)
    p.add_argument("--out", default=None)
    p.set_defaults(func=_fields_stats)

    p = sub.add_parser("xmf", help="write a ParaView XDMF time series for field snapshots")
    p.add_argument("--case-dir", default=".")
    p.add_argument("--prefix", default=None)
    p.add_argument("--out", default=None)
    p.set_defaults(func=_fields_xmf)


# ── vortex ─────────────────────────────────────────────────────────────────

def _vortex_animate(a):
    vortex.animate_q_isosurface(
        fields.FieldSeries(a.case_dir), a.out, start=a.start, end=a.end, stride=a.stride,
        max_frames=a.max_frames, qval=a.qval, vmax=a.vmax, cmap=a.cmap, fps=a.fps,
        width=a.width, height=a.height, plot=a.plot, view=a.view, cam_pos=a.cam_pos,
        azimuth=a.azimuth, elevation=a.elevation, roll=a.roll, zoom=a.zoom,
        zoom_slice=a.zoom_slice, slice_x=a.slice_x, parallel=a.parallel)


def _add_vortex(sub):
    p = sub.add_parser("animate", help="Q-criterion iso-surface movie (needs VTK + an X display)")
    p.add_argument("--case-dir", default=".")
    p.add_argument("--out", default="q.mp4", help=".mp4 or .gif")
    p.add_argument("--start", type=int, default=None, help="first step to include")
    p.add_argument("--end", type=int, default=None, help="last step to include")
    p.add_argument("--stride", type=int, default=1, help="use every Nth snapshot")
    p.add_argument("--max-frames", type=int, default=None)
    p.add_argument("--qval", type=float, default=5.0, help="Q iso-value")
    p.add_argument("--vmax", type=float, default=None, help="top of the |u| colour range (default 1.7*Ub_target)")
    p.add_argument("--cmap", default="RdYlBu_r")
    p.add_argument("--fps", type=int, default=10)
    p.add_argument("--width", type=int, default=None)
    p.add_argument("--height", type=int, default=600)
    p.add_argument("--plot", choices=["both", "q", "slice"], default="both")
    p.add_argument("--view", choices=list(vortex.VIEWS), default="iso")
    p.add_argument("--cam-pos", type=float, nargs=3, metavar=("DX", "DY", "DZ"))
    p.add_argument("--azimuth", type=float, default=0.0)
    p.add_argument("--elevation", type=float, default=0.0)
    p.add_argument("--roll", type=float, default=0.0)
    p.add_argument("--zoom", type=float, default=1.3, help="zoom of the Q panel")
    p.add_argument("--zoom-slice", type=float, default=1.0)
    p.add_argument("--slice-x", type=float, default=0.5, help="x-normal slice location as a fraction of Lx")
    p.add_argument("--parallel", action="store_true", help="orthographic projection")
    p.set_defaults(func=_vortex_animate)


# ── probes ─────────────────────────────────────────────────────────────────

def _probes_slice_xmf(a):
    slc = probes.SliceProbe.read(a.meta)
    slc.write_xmf(a.out, dt=a.dt, t0=a.t0)


def _add_probes(sub):
    p = sub.add_parser("slice-xmf", help="write a ParaView XDMF time series for a slice probe")
    p.add_argument("meta", help="<base>_meta.txt")
    p.add_argument("--out", default=None)
    p.add_argument("--dt", type=float, default=1.0)
    p.add_argument("--t0", type=float, default=0.0)
    p.set_defaults(func=_probes_slice_xmf)


# ── ibm_surface ────────────────────────────────────────────────────────────

def _ibm_surface_convert(a):
    if a.pvd:
        ibm_surface.IBMSurface.to_pvd(a.path, a.pvd)
    else:
        surf = ibm_surface.IBMSurface.read(a.path)
        surf.to_vtp(a.out)


def _add_ibm_surface(sub):
    p = sub.add_parser("convert", help="convert IBM surface sample(s) to .vtp/.pvd")
    p.add_argument("path", help="one surface.<step>.bin file, or a glob for --pvd")
    p.add_argument("--out", default=None, help=".vtp output for a single file")
    p.add_argument("--pvd", default=None, help=".pvd output; treats `path` as a glob pattern")
    p.set_defaults(func=_ibm_surface_convert)


# ── sdf ────────────────────────────────────────────────────────────────────

def _sdf_plot(a):
    s = sdf.SDF.read(a.case_dir, a.sdf)
    if a.y_slice is not None:
        s.plot_xz(a.y_slice).savefig(a.out or "sdf_xz.png", dpi=150)
    if a.z_slice is not None:
        s.plot_xy(a.z_slice).savefig(a.out or "sdf_xy.png", dpi=150)


def _add_sdf(sub):
    p = sub.add_parser("plot", help="plot x-z and x-y SDF slices")
    p.add_argument("--case-dir", default=".")
    p.add_argument("--sdf", default="SDF_in")
    p.add_argument("--y-slice", type=float, default=None)
    p.add_argument("--z-slice", type=float, default=None)
    p.add_argument("--out", default=None)
    p.set_defaults(func=_sdf_plot)


# ── uav ────────────────────────────────────────────────────────────────────

def _uav_disk_animate(a):
    path = uav.UAVPath.read(a.path_file)
    path.disk_animation(radius=a.radius, times_from=a.times_from, pvd_out=a.out)


def _uav_drone_animate(a):
    path = uav.UAVPath.read(a.path_file)
    path.drone_animation(times_from=a.times_from, tilt=a.tilt, stl_out=a.out)


def _add_uav(sub):
    p = sub.add_parser("disk-animate", help="render a uav_path_file as a moving-disk animation")
    p.add_argument("path_file")
    p.add_argument("--radius", type=float, default=0.15)
    p.add_argument("--times-from", default=None)
    p.add_argument("--out", default="uav_path.pvd")
    p.set_defaults(func=_uav_disk_animate)

    p = sub.add_parser("drone-animate", help="render a uav_path_file as a quadcopter body animation")
    p.add_argument("path_file")
    p.add_argument("--times-from", default=None)
    p.add_argument("--tilt", action="store_true")
    p.add_argument("--out", default="uav_drone.stl")
    p.set_defaults(func=_uav_drone_animate)


# ── rsb ────────────────────────────────────────────────────────────────────

def _rsb_plot(a):
    stats = rsb.RSBStats.read(a.case_dir)
    stats.plot(last=a.last, out_dir=a.out_dir)


def _add_rsb(sub):
    p = sub.add_parser("plot", help="plot the Reynolds-stress budget")
    p.add_argument("--case-dir", default=".")
    p.add_argument("--last", type=int, default=None)
    p.add_argument("--out-dir", default="plots")
    p.set_defaults(func=_rsb_plot)


# ── inflow ─────────────────────────────────────────────────────────────────

def _inflow_mirror(a):
    inflow.mirror_half_channel(a.input, a.output, h=a.h)


def _inflow_check_donor(a):
    donor = inflow.InflowDonor.read(a.base)
    inflow.check_donor(donor, reference=a.reference, grid=a.grid, output=a.out)


def _inflow_opt(a):
    state = inflow.read_inflow_opt(a.case_dir)
    print(state.summary())


def _add_inflow(sub):
    p = sub.add_parser("mirror", help="mirror a half-channel SEM inflow profile to a full channel")
    p.add_argument("input")
    p.add_argument("output")
    p.add_argument("--h", type=float, default=None)
    p.set_defaults(func=_inflow_mirror)

    p = sub.add_parser("check-donor", help="verify a SEM/ESEM inflow donor plane vs a reference profile")
    p.add_argument("base", help="<base>_meta.txt (without the suffix)")
    p.add_argument("--reference", default=None)
    p.add_argument("--grid", default="fields/grid.out")
    p.add_argument("--out", default="inflow_donor_check.png")
    p.set_defaults(func=_inflow_check_donor)

    p = sub.add_parser("opt-state", help="summarize a SEM inflow-optimization restart file")
    p.add_argument("--case-dir", default=".")
    p.set_defaults(func=_inflow_opt)


# ── runlog ─────────────────────────────────────────────────────────────────

def _runlog_plot(a):
    log = runlog.RunLog.read(a.path)
    variables = a.variables.split(",") if a.variables else None
    log.plot(variables=variables, out=a.out or "run_log.png")


def _add_runlog(sub):
    p = sub.add_parser("plot", help="plot solver diagnostics from a run.log")
    p.add_argument("path", nargs="?", default="run.log")
    p.add_argument("-v", "--variables", default=None,
                   help="comma-separated variable names (default: Umean,Umax,CFLc,CFLv,dt)")
    p.add_argument("--out", default=None)
    p.set_defaults(func=_runlog_plot)


# ── top-level dispatch ────────────────────────────────────────────────────

_GROUPS = {
    "particles": _add_particles,
    "fields": _add_fields,
    "probes": _add_probes,
    "ibm-surface": _add_ibm_surface,
    "sdf": _add_sdf,
    "uav": _add_uav,
    "rsb": _add_rsb,
    "inflow": _add_inflow,
    "runlog": _add_runlog,
    "vortex": _add_vortex,
}


def main(argv=None):
    ap = argparse.ArgumentParser(prog="dopamine-post",
                                 description="Post-processing for dopamine-fdm cases.")
    group_sub = ap.add_subparsers(dest="group", required=True)
    for name, add_fn in _GROUPS.items():
        gp = group_sub.add_parser(name)
        cmd_sub = gp.add_subparsers(dest="command", required=True)
        add_fn(cmd_sub)

    args = ap.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    sys.exit(main())
