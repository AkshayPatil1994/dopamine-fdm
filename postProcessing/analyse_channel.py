#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 DOPAMINE contributors
# SPDX-License-Identifier: AGPL-3.0-only

"""
NOTICE 

-- LLM GENERATED CODE BELOW --

Wall-normal profiles from a channel run, compared against channel-DNS reference data.

Reads binary field snapshots from fields/ (see snapshot_io.py), extracts one or more
streamwise (x) stations, averages over z and over a window of snapshots, and forms the
resolved Reynolds stresses. Averaging is spanwise+temporal only (never streamwise), so
profiles at successive x show how an inflow develops with fetch (useful with SEM).
Profiles are overlaid on Moser-Kim-Mansour DNS and on a measured wind-tunnel inflow
profile (windtunnel_inflow.csv: y, U, Iu, Iv, Iw); both already use the solver's axes
(x streamwise, y wall-normal, z spanwise), so no remapping is needed.

Note: stresses reported here are *resolved* only — a wall-modelled/coarse LES carries
part of the stress in the subgrid model, so a near-wall deficit vs DNS is expected.

Usage
-----
    python3 analyse_channel.py --x 2 4 6 8 \\
        --dns-dir /mnt/storage1/dopamine/validation/semChannel \\
        --step-start 900 --re-tau 392.24 --output profiles.png
"""

import argparse
import sys
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

sys.path.insert(0, str(Path(__file__).parent))
import snapshot_io


# ----------------------------------------------------------------------------
# statistics
# ----------------------------------------------------------------------------
def profile_at_station(fields_dir, prefix, steps, x_target, xm, sgs_model, sediment_flag):
    """Spanwise- and snapshot-averaged mean and resolved stresses at the x-plane nearest x_target."""
    ix = int(np.argmin(np.abs(xm - x_target)))
    x_sel = float(xm[ix])

    s1 = s2 = None
    n = 0
    for step in steps:
        fpath = Path(fields_dir) / f"{prefix}.{step}"
        snap = snapshot_io.read_snapshot(fpath, sgs_model=sgs_model, sediment_flag=sediment_flag)
        U, V, W = snap["U"][ix], snap["V"][ix], snap["W"][ix]   # each (nym, nzm)

        if s1 is None:
            nym = U.shape[0]
            s1 = np.zeros((nym, 3))
            s2 = np.zeros((nym, 6))

        s1[:, 0] += U.sum(axis=1)
        s1[:, 1] += V.sum(axis=1)
        s1[:, 2] += W.sum(axis=1)
        s2[:, 0] += (U * U).sum(axis=1)
        s2[:, 1] += (V * V).sum(axis=1)
        s2[:, 2] += (W * W).sum(axis=1)
        s2[:, 3] += (U * V).sum(axis=1)
        s2[:, 4] += (U * W).sum(axis=1)
        s2[:, 5] += (V * W).sum(axis=1)
        
        n += U.shape[1]

    if n == 0:
        sys.exit(f"ERROR: no snapshots found near x = {x_target}")

    mean = s1 / n
    prod_mean = s2 / n
    rey = prod_mean - np.stack(
        [mean[:, 0] * mean[:, 0], mean[:, 1] * mean[:, 1], mean[:, 2] * mean[:, 2],
         mean[:, 0] * mean[:, 1], mean[:, 0] * mean[:, 2], mean[:, 1] * mean[:, 2]], axis=1)
    return x_sel, mean, rey, n


# Sign each component picks up under the wall-to-wall mirror y -> Ly - y
# (V flips sign, U and W don't; a stress component's sign is the product of
# its two velocity signs).
_MEAN_MIRROR_SIGN = np.array([1.0, -1.0, 1.0])                    # U, V, W
_REY_MIRROR_SIGN = np.array([1.0, 1.0, 1.0, -1.0, 1.0, -1.0])      # uu,vv,ww,uv,uw,vw


def fold_profile(mean, rey):
    """Average a profile with its mirror image about the centreline y = Ly/2.

    Only valid for a channel with matching (no-slip/no-slip) walls and a
    wall-normal grid that is itself symmetric about the centreline -- both
    checked by the caller before this is used.
    """
    mean_f = 0.5 * (mean + _MEAN_MIRROR_SIGN * mean[::-1])
    rey_f = 0.5 * (rey + _REY_MIRROR_SIGN * rey[::-1])
    return mean_f, rey_f


# ----------------------------------------------------------------------------
# immersed-boundary SDF
# ----------------------------------------------------------------------------
def load_sdf(sdf_path, nxm, nym, nzm):
    """Read the cell-centre SDF written by fdm-dopamine (ibm_input_mode=1).

    Binary layout: big-endian float64, Fortran column-major order,
    shape (nxm+2, nym+2, nzm) with ghost layers in x and y only. Returns
    phi with shape (nxm, nym, nzm); phi < 0 marks solid cells.
    """
    nxg, nyg = nxm + 2, nym + 2
    raw = np.fromfile(sdf_path, dtype=">f8")
    expected = nxg * nyg * nzm
    if raw.size != expected:
        sys.exit(f"ERROR: {sdf_path} size mismatch: got {raw.size} elements, "
                  f"expected {expected} ({nxg}x{nyg}x{nzm})")
    phi_full = raw.reshape((nxg, nyg, nzm), order="F")
    return phi_full[1:-1, 1:-1, :]


# Contrasting marker colors for the tracked x-locations, chosen to stay legible
# against the magma_r colorplot background (ColorBrewer Set1, first four).
_LINE_COLORS = ["#e41a1c", "#377eb8", "#4daf4a", "#ff7f00"]


def select_equidistant_x(x_list, max_n=4):
    """Pick up to max_n locations, equally spaced between min(x_list) and max(x_list).

    With max_n or fewer stations requested, all of them are used as-is.
    Otherwise max_n equidistant targets spanning [min, max] are returned
    (the nearest cell-centre plane to each is used by the caller).
    """
    xs = sorted(x_list)
    if len(xs) <= max_n:
        return xs
    return list(np.linspace(xs[0], xs[-1], max_n))


# ----------------------------------------------------------------------------
# z-normal velocity-magnitude animation
# ----------------------------------------------------------------------------
def animate_velocity_gif(fields_dir, prefix, steps, xm, ym, zm, sgs_model, sediment_flag,
                          z_target, output, fps, phi=None, line_x=None, inflow=None):
    """Animate the (x, y) velocity-magnitude slice at the z-plane nearest z_target.

    When phi (the cell-centre SDF, shape (nxm, nym, nzm)) is given, solid
    cells (phi < 0) at that z-plane are masked out and rendered grey.

    line_x, when given, marks each x-location with a dashed vertical line on
    the colorplot and adds smaller subplots below it showing the wall-normal
    profiles of spanwise-averaged streamwise velocity and turbulence
    intensities (Iu, Iv, Iw, from spanwise fluctuations) at that x-plane,
    redrawn every frame. Each location and its profile lines share a
    contrasting color. When inflow (from load_inflow) is given, the
    wind-tunnel reference profile is overlaid (static) on each subplot.
    """
    iz = int(np.argmin(np.abs(zm - z_target)))
    z_sel = float(zm[iz])

    solid = phi[:, :, iz] < 0 if phi is not None else None

    line_x = line_x or []
    ix_lines = [int(np.argmin(np.abs(xm - xq))) for xq in line_x]
    x_line_sel = [float(xm[ix]) for ix in ix_lines]
    colors = [_LINE_COLORS[i % len(_LINE_COLORS)] for i in range(len(ix_lines))]

    frames = []
    # profiles[quantity][location] -> list of per-frame (nym,) arrays
    profiles_U = [[] for _ in ix_lines]
    profiles_Iu = [[] for _ in ix_lines]
    profiles_Iv = [[] for _ in ix_lines]
    profiles_Iw = [[] for _ in ix_lines]
    for step in steps:
        snap = snapshot_io.read_snapshot(Path(fields_dir) / f"{prefix}.{step}",
                                          sgs_model=sgs_model, sediment_flag=sediment_flag)
        U, V, W = snap["U"][:, :, iz], snap["V"][:, :, iz], snap["W"][:, :, iz]
        speed = np.sqrt(U * U + V * V + W * W)
        if solid is not None:
            speed = np.ma.masked_where(solid, speed)
        frames.append(speed)

        for k, ix in enumerate(ix_lines):
            Ul, Vl, Wl = snap["U"][ix], snap["V"][ix], snap["W"][ix]   # (nym, nzm)
            Ubar_y = Ul.mean(axis=1)
            up = Ul - Ubar_y[:, None]
            vp = Vl - Vl.mean(axis=1)[:, None]
            wp = Wl - Wl.mean(axis=1)[:, None]
            with np.errstate(divide="ignore", invalid="ignore"):
                safe_U = np.where(np.abs(Ubar_y) > 1e-8, Ubar_y, np.nan)
                Iu = np.sqrt(np.mean(up * up, axis=1)) / safe_U
                Iv = np.sqrt(np.mean(vp * vp, axis=1)) / safe_U
                Iw = np.sqrt(np.mean(wp * wp, axis=1)) / safe_U
            profiles_U[k].append(Ubar_y)
            profiles_Iu[k].append(Iu)
            profiles_Iv[k].append(Iv)
            profiles_Iw[k].append(Iw)

    vmin, vmax = 0.0, 20.0  # fixed colour scale so successive frames stay comparable

    cmap = matplotlib.colormaps["magma_r"].copy()
    cmap.set_bad(color="0.5")

    # Size the figure so the (x, y) colorplot -- drawn with aspect="equal" --
    # fills its row edge-to-edge instead of letterboxing inside a row height
    # picked without regard to the data's actual aspect ratio.
    fig_w = 12.0
    img_w_frac = 0.90  # fraction of fig_w left for the image after the colorbar
    data_aspect = (ym[-1] - ym[0]) / (xm[-1] - xm[0])
    img_h = fig_w * img_w_frac * data_aspect

    if ix_lines:
        prof_h = 2.4  # inches, fixed height for the profile row
        fig = plt.figure(figsize=(fig_w, img_h + prof_h + 1.1))
        gs = fig.add_gridspec(2, 4, height_ratios=[img_h, prof_h], hspace=0.08, wspace=0.35,
                               left=0.06, right=0.98, top=0.95, bottom=0.10)
        ax = fig.add_subplot(gs[0, :])
        ax_prof = [fig.add_subplot(gs[1, j]) for j in range(4)]
    else:
        fig, ax = plt.subplots(figsize=(fig_w, img_h + 0.9))
        ax_prof = []

    im = ax.imshow(frames[0].T, origin="lower", aspect="equal",
                    extent=[xm[0], xm[-1], ym[0], ym[-1]],
                    vmin=vmin, vmax=vmax, cmap=cmap)
    fig.colorbar(im, ax=ax, label=r"$|U|$", shrink=0.5)
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")
    title = ax.set_title(f"z = {z_sel:.3g}, step {steps[0]}")

    for x_sel, c in zip(x_line_sel, colors):
        ax.axvline(x_sel, color=c, ls="--", lw=1.6)

    prof_titles = [r"$U$", r"$I_u$", r"$I_v$", r"$I_w$"]
    prof_data = [profiles_U, profiles_Iu, profiles_Iv, profiles_Iw]
    inflow_keys = [None, "Iu", "Iv", "Iw"]
    prof_lines = [[] for _ in ax_prof]
    for j, (a, ttl, data) in enumerate(zip(ax_prof, prof_titles, prof_data)):
        for k, c in enumerate(colors):
            line, = a.plot(data[k][0], ym, color=c, lw=1.3,
                            label=f"x = {x_line_sel[k]:.2f}")
            prof_lines[j].append(line)
        if inflow is not None:
            ikey = "U" if j == 0 else inflow_keys[j]
            a.plot(inflow[ikey], inflow["y"], "ko", markerfacecolor="none",
                   markeredgecolor="black", markevery=2, lw=1.0, ms=4,
                   label="wind tunnel", zorder=5)
        a.set_xlabel(ttl, fontsize=9)
        if j == 0:
            a.set_ylabel(r"$y$", fontsize=9)
        a.tick_params(labelsize=7)
        a.grid(alpha=0.25)
    if ax_prof:
        ax_prof[0].legend(fontsize=6, frameon=False)
        # Not tight_layout: it would re-space the rows using their nominal
        # height_ratios, undoing the aspect-matched sizing set above.
    else:
        fig.tight_layout()

    def update(i):
        im.set_data(frames[i].T)
        title.set_text(f"z = {z_sel:.3g}, step {steps[i]}")
        for j, data in enumerate(prof_data):
            for k in range(len(ix_lines)):
                prof_lines[j][k].set_xdata(data[k][i])
        return (im, title, *[ln for group in prof_lines for ln in group])

    ani = FuncAnimation(fig, update, frames=len(frames), blit=False)
    ani.save(output, writer=PillowWriter(fps=fps))
    plt.close(fig)
    print(f"wrote {output} ({len(frames)} frames, z = {z_sel:.4g})")


# ----------------------------------------------------------------------------
# DNS reference
# ----------------------------------------------------------------------------
def read_dat(path, ncol):
    rows = []
    for line in Path(path).read_text().splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        p = s.split()
        if len(p) < ncol:
            continue
        try:
            rows.append([float(v) for v in p[:ncol]])
        except ValueError:
            continue
    return np.array(rows)


def load_dns(dns_dir, u_tau, h):
    """Return y, U, and stresses in solver axes, or None when the files are absent."""
    d = Path(dns_dir)
    means, rey = d / "MKM_MEANS.dat", d / "MKM_REYSTRESS.dat"
    if not (means.exists() and rey.exists()):
        return None
    m, r = read_dat(means, 7), read_dat(rey, 8)
    if len(m) != len(r):
        print(f"WARNING: DNS row counts differ ({len(m)} vs {len(r)}), skipping reference")
        return None
    return {
        "y": m[:, 0] * h,
        "U": m[:, 2] * u_tau,
        # DNS axes (y wall-normal, z spanwise) already match the solver's
        # binary field snapshots, so no remap is needed here.
        "uu": r[:, 2] * u_tau**2,
        "vv": r[:, 3] * u_tau**2,
        "ww": r[:, 4] * u_tau**2,
        "uv": r[:, 5] * u_tau**2,
    }


def load_inflow(path):
    """Read the wind-tunnel inflow profile: y, U, Iu, Iv, Iw (turbulence intensities)."""
    rows = read_dat(path, 5)
    if rows.size == 0:
        return None
    y, U, Iu, Iv, Iw = rows.T
    return {"y": y, "U": U, "Iu": Iu, "Iv": Iv, "Iw": Iw}


# ----------------------------------------------------------------------------
# plotting
# ----------------------------------------------------------------------------
def make_figure(stations, dns, inflow, args):
    fig, ax = plt.subplots(1, 4, figsize=(18, 5.5))
    colors = plt.cm.viridis(np.linspace(0.05, 0.85, len(stations)))

    for (x_sel, y, mean, rey, _), c in zip(stations, colors):
        lbl = f"x = {x_sel:.2f}"
        U = mean[:, 0]
        ax[0].plot(U, y, color=c, label=lbl)
        ax[1].plot(np.sqrt(np.clip(rey[:, 0], 0, None)) / U, y, color=c, label=lbl)
        ax[2].plot(np.sqrt(np.clip(rey[:, 1], 0, None)) / U, y, color=c, label=lbl)
        ax[3].plot(np.sqrt(np.clip(rey[:, 2], 0, None)) / U, y, color=c, label=lbl)

    keys = ["U", "Iu", "Iv", "Iw"]
    for a, key in zip(ax, keys):
        if dns is not None and key == "U":
            a.plot(dns["U"], dns["y"], "k-.", lw=1.4, label="DNS", zorder=5)
        elif dns is not None:
            dns_I = {"Iu": "uu", "Iv": "vv", "Iw": "ww"}[key]
            a.plot(np.sqrt(np.clip(dns[dns_I], 0, None)) / dns["U"], dns["y"],
                   "k-.", lw=1.4, label="DNS", zorder=5)
        if inflow is not None:
            a.plot(inflow[key], inflow["y"], "ko", lw=1.6,
                   markerfacecolor="none", markeredgecolor="black", markevery=2,
                   label="wind tunnel", zorder=5)

    titles = [r"$U$", r"$I_u$", r"$I_v$", r"$I_w$"]
    for a, xl in zip(ax, titles):
        a.set_xlabel(xl)
        a.set_ylabel(r"$y$")
        a.grid(alpha=0.25)
        a.legend(fontsize=8, frameon=False)
    step_range = f"{args.step_start if args.step_start is not None else 'first'}" \
                 f"..{args.step_end if args.step_end is not None else 'last'}"
    fig.suptitle(f"{args.prefix}: spanwise- and snapshot-averaged profiles "
                 f"(steps {step_range})", y=0.99)
    fig.tight_layout()
    fig.savefig(args.output, dpi=150)
    print(f"wrote {args.output}")


def _str2bool(v):
    if v.lower() in ("true", "1", "yes"):
        return True
    if v.lower() in ("false", "0", "no"):
        return False
    raise argparse.ArgumentTypeError(f"expected true/false, got {v!r}")


# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Compare simulated channel profiles (binary fields/ snapshots) against DNS reference data.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fields-dir", default="fields", help="directory holding the binary snapshots")
    ap.add_argument("--prefix", default="channel", help="snapshot filename prefix (fields/<prefix>.<step>)")
    ap.add_argument("--x", type=float, nargs="+", required=True,
                    help="streamwise station(s); the nearest cell-centre plane is used")
    ap.add_argument("--dns-dir", default=None,
                    help="directory holding MKM_MEANS.dat and MKM_REYSTRESS.dat")
    ap.add_argument("--inflow", default=None,
                    help="wind-tunnel inflow profile CSV, y U Iu Iv Iw "
                         "(default: windtunnel_inflow.csv beside this script)")
    ap.add_argument("--step-start", type=int, default=None, help="discard snapshot steps below this")
    ap.add_argument("--step-end", type=int, default=None, help="discard snapshot steps above this")
    ap.add_argument("--u-tau", type=float, default=1.0, help="friction velocity of the case")
    ap.add_argument("--h", type=float, default=1.0, help="channel half-height")
    ap.add_argument("--re-tau", type=float, default=None,
                    help="enables the law-of-the-wall panel and y+ scaling")
    ap.add_argument("--sgs-model", type=int, default=None,
                    help="override sgs_model from input_parameters (0=DNS, nonzero=LES with nu_t block)")
    ap.add_argument("--sediment-flag", type=int, default=None,
                    help="override sediment_flag from input_parameters")
    ap.add_argument("--no-fold", action="store_true",
                    help="disable the automatic no-slip/no-slip symmetry fold")
    ap.add_argument("--output", default="channel_profiles.png")
    ap.add_argument("--csv", default=None, help="also write the extracted profiles here")
    ap.add_argument("--gif-z", type=float, default=None,
                     help="z-station for the velocity-magnitude animation (default: mid-span)")
    ap.add_argument("--gif-output", default="velocity.gif",
                     help="output path for the animated z-normal velocity-magnitude slice")
    ap.add_argument("--gif-fps", type=float, default=1.0, help="frame rate for --gif-output")
    ap.add_argument("--no-gif", action="store_true", help="skip writing the velocity-magnitude animation")
    ap.add_argument("--ibm", type=_str2bool, nargs="?", const=True, default=False,
                     help="mask solid (immersed-boundary) cells grey in the animation using SDF_in")
    ap.add_argument("--sdf-in", default="SDF_in",
                     help="cell-centre SDF file used to mask solid cells when --ibm true (default: SDF_in)")
    args = ap.parse_args()

    case_dir = Path(__file__).parent
    params = {}
    ip_path = case_dir / "input_parameters"
    if ip_path.exists():
        params = snapshot_io.parse_input_parameters(ip_path)
    sgs_model = args.sgs_model if args.sgs_model is not None else int(params.get("sgs_model", 0))
    sediment_flag = args.sediment_flag if args.sediment_flag is not None else int(params.get("sediment_flag", 0))

    # &BOUNDARY_CONDITIONS defaults to both walls no-slip when the group is omitted.
    bc_face_ylo = int(params.get("bc_face_ylo", 1))
    bc_face_yhi = int(params.get("bc_face_yhi", 1))
    fold_walls = bc_face_ylo == 1 and bc_face_yhi == 1 and not args.no_fold

    all_steps = [s for s, _ in snapshot_io.list_snapshots(args.fields_dir, args.prefix)]
    if not all_steps:
        sys.exit(f"ERROR: no snapshots matching '{args.prefix}.<step>' found in {args.fields_dir}")
    steps = [s for s in all_steps
             if (args.step_start is None or s >= args.step_start)
             and (args.step_end is None or s <= args.step_end)]
    if not steps:
        sys.exit(f"ERROR: no snapshots in steps = [{args.step_start}, {args.step_end}]; "
                 f"available range is {all_steps[0]}..{all_steps[-1]}")
    print(f"averaging {len(steps)} of {len(all_steps)} snapshots, steps = {steps[0]}..{steps[-1]}")

    s0 = snapshot_io.read_snapshot(Path(args.fields_dir) / f"{args.prefix}.{steps[0]}",
                                    sgs_model=sgs_model, sediment_flag=sediment_flag)
    xm, ym = s0["xm"], s0["ym"]
    print(f"mesh: nxm={xm.size}, nym={ym.size}, nzm={s0['zm'].size}, "
          f"y = {ym[0]:.4g}..{ym[-1]:.4g}")

    if fold_walls:
        Ly = float(s0["y"][-1])
        grid_asym = np.max(np.abs(ym + ym[::-1] - Ly))
        if grid_asym > 1e-6 * Ly:
            print(f"WARNING: wall-normal grid is not mirror-symmetric about the centreline "
                  f"(max deviation {grid_asym:.3g}); skipping the no-slip/no-slip fold")
            fold_walls = False
        else:
            print("both walls no-slip (bc_face_ylo=bc_face_yhi=1): "
                  "folding profiles about the centreline")

    y_out = ym[:len(ym) // 2] if fold_walls else ym

    stations = []
    for xq in args.x:
        x_sel, mean, rey, n = profile_at_station(args.fields_dir, args.prefix, steps, xq, xm,
                                                   sgs_model, sediment_flag)
        if fold_walls:
            mean, rey = fold_profile(mean, rey)
            half = len(ym) // 2
            mean, rey = mean[:half], rey[:half]
        stations.append((x_sel, y_out, mean, rey, n))
        umax = np.nanmax(mean[:, 0])
        print(f"  x = {xq:8.3f} -> plane {x_sel:8.3f}, {n} samples/level, max U = {umax:.4f}")

    dns = load_dns(args.dns_dir, args.u_tau, args.h) if args.dns_dir else None
    if args.dns_dir and dns is None:
        print(f"WARNING: no DNS files found in {args.dns_dir}; plotting simulation only")
    inflow_path = Path(args.inflow) if args.inflow else case_dir / "windtunnel_inflow.csv"
    inflow = load_inflow(inflow_path) if inflow_path.exists() else None
    if inflow is None:
        print(f"WARNING: no readable wind-tunnel inflow profile found at {inflow_path}; plotting simulation only")

    make_figure(stations, dns, inflow, args)

    if not args.no_gif:
        phi = None
        if args.ibm:
            sdf_path = Path(args.sdf_in)
            if not sdf_path.exists():
                sys.exit(f"ERROR: --ibm true but SDF file not found: {sdf_path}")
            zm = s0["zm"]
            phi = load_sdf(sdf_path, xm.size, ym.size, zm.size)
            print(f"loaded SDF from {sdf_path}: shape {phi.shape}, "
                  f"{np.count_nonzero(phi < 0)} solid cells")
        zm = s0["zm"]
        z_target = args.gif_z if args.gif_z is not None else float(zm[len(zm) // 2])
        line_x = select_equidistant_x(args.x, max_n=4)
        animate_velocity_gif(args.fields_dir, args.prefix, steps, xm, ym, zm,
                              sgs_model, sediment_flag, z_target, args.gif_output, args.gif_fps,
                              phi=phi, line_x=line_x, inflow=inflow)

    if args.csv:
        with open(args.csv, "w") as fh:
            fh.write("# x, y, U, V, W, uu, vv, ww, uv, uw, vw, k  (resolved stresses)\n")
            for x_sel, y, mean, rey, _ in stations:
                for k in range(len(y)):
                    fh.write(f"{x_sel:.6e},{y[k]:.6e}," +
                             ",".join(f"{v:.6e}" for v in mean[k]) + "," +
                             ",".join(f"{v:.6e}" for v in rey[k]) + "," +
                             f"{0.5 * rey[k, :3].sum():.6e}\n")
        print(f"wrote {args.csv}")

    if dns is not None:
        print("\nNote: stresses above are resolved only; a wall-modelled LES keeps part")
        print("of the near-wall stress in the subgrid model, so a deficit there is expected.")


if __name__ == "__main__":
    main()
