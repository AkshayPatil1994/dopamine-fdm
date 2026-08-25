#!/usr/bin/env python3
"""
compute_stats.py — Time- and plane-averaged statistics for fdm-dopamine channel flow.

Reads snapshots from fields/, averages over time and the two homogeneous
directions (x, z), then applies channel symmetry about y = Ly/2 to improve
statistical convergence.

Computed quantities (all as wall-normal y-profiles):
  U_mean, V_mean, W_mean  — time- and plane-averaged velocity components
  urms, vrms, wrms        — RMS velocity fluctuations
  uv                      — Reynolds shear stress  <u'v'>
  TKE                     — turbulent kinetic energy  0.5*(urms²+vrms²+wrms²)
  prod                    — TKE production            -<u'v'> dU/dy
  eps_res                 — resolved pseudo-dissipation
  eps_sgs                 — SGS dissipation  (requires nu_t)

Outputs:
  stats/<prefix>_<avg_start>_<avg_end>.npz   — all profiles in one NumPy archive
  stats/<name>.bin                           — individual little-endian float64 files
  plots/<prefix>_<avg_start>_<avg_end>.png   — wall-unit profile plots (automatic)

Usage:
  python3 compute_stats.py \\
      --avg_start 100000 --interval 500 --avg_end 150000 --prefix leschan_550
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
import snapshot_io


# ── argument parsing ──────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description='Time- and plane-averaged channel flow statistics.'
    )
    p.add_argument('--avg_start', type=int, required=True,
                   help='First snapshot step included in the average.')
    p.add_argument('--interval',  type=int, required=True,
                   help='Step interval between consecutive snapshots (nsave).')
    p.add_argument('--avg_end',   type=int, required=True,
                   help='Last snapshot step included in the average.')
    p.add_argument('--prefix',    type=str, required=True,
                   help='Snapshot filename prefix (e.g. leschan_550).')
    p.add_argument('--fields_dir', type=str, default='fields',
                   help='Directory that contains the snapshot files (default: fields).')
    p.add_argument('--stats_dir',  type=str, default='stats',
                   help='Directory for output files (default: stats).')
    p.add_argument('--nu',   type=float, default=None,
                   help='Kinematic viscosity ν. Defaults to value in input_parameters.')
    p.add_argument('--no_sgs', action='store_true',
                   help='Do not read nu_t even when sgs_model != 0.')
    p.add_argument('--out_fig', type=str, default='plots/stats.png')
    return p.parse_args()


# ── helpers ───────────────────────────────────────────────────────────────────

def xz_mean(arr):
    """Average a (nxm, nym, nzm) array over axes 0 and 2 → (nym,)."""
    return arr.mean(axis=(0, 2))


def grad_all(U, V, W, ym, dx, dz):
    """
    Return all nine velocity-gradient components at cell centres.

    Uses numpy.gradient (central differences in the interior, first-order
    one-sided at boundaries).  x and z are uniform with spacings dx and dz;
    y is non-uniform and requires the full ym coordinate array.

    Returns a dict  {'dU_dx', 'dU_dy', 'dU_dz',
                     'dV_dx', 'dV_dy', 'dV_dz',
                     'dW_dx', 'dW_dy', 'dW_dz'}.
    """
    return {
        'dU_dx': np.gradient(U, dx, axis=0),
        'dU_dy': np.gradient(U, ym, axis=1),
        'dU_dz': np.gradient(U, dz, axis=2),
        'dV_dx': np.gradient(V, dx, axis=0),
        'dV_dy': np.gradient(V, ym, axis=1),
        'dV_dz': np.gradient(V, dz, axis=2),
        'dW_dx': np.gradient(W, dx, axis=0),
        'dW_dy': np.gradient(W, ym, axis=1),
        'dW_dz': np.gradient(W, dz, axis=2),
    }


# ── symmetry folding ──────────────────────────────────────────────────────────

def sym_fold(a):
    """Symmetric (even) fold about the channel centreline."""
    return 0.5 * (a + a[::-1])


def antisym_fold(a):
    """Anti-symmetric (odd) fold about the channel centreline.

    For <u'v'>: negative in the lower half, positive in the upper half.
    After folding the result carries the sign of the lower half (y < Ly/2).
    """
    return 0.5 * (a - a[::-1])


# ── plotting ─────────────────────────────────────────────────────────────────

_LES_STYLE = dict(
    color='k',
    marker='o',
    linestyle='none',
    markerfacecolor='none',
    markeredgecolor='k',
    markersize=4.5,
    markeredgewidth=0.9,
)


def half_channel(data, u_tau, delta_nu, nu):
    """Return wall-unit profiles for the lower half-channel."""
    nym = data['ym'].size
    half = nym // 2

    yp = data['ym'][:half] / delta_nu
    Up = data['U_mean'][:half] / u_tau
    urp = data['urms'][:half] / u_tau
    vrp = data['vrms'][:half] / u_tau
    wrp = data['wrms'][:half] / u_tau
    uvp = -data['uv'][:half] / u_tau**2
    TKEp = data['TKE'][:half] / u_tau**2

    fac = nu / u_tau**4
    prod_p = data['prod'][:half] * fac
    epsr_p = data['eps_res'][:half] * fac
    epss_p = data['eps_sgs'][:half] * fac
    eps_tot = epsr_p + epss_p

    return dict(
        yp=yp,
        Up=Up,
        urp=urp,
        vrp=vrp,
        wrp=wrp,
        uvp=uvp,
        TKEp=TKEp,
        prod_p=prod_p,
        epsr_p=epsr_p,
        epss_p=epss_p,
        eps_tot=eps_tot,
    )


def load_dns_reference(dns_dir):
    mean = np.loadtxt(dns_dir / 'mean_prof.txt', comments='%')
    vel = np.loadtxt(dns_dir / 'vel_fluc_prof.txt', comments='%')
    rste_k = np.loadtxt(dns_dir / 'RSTE_k_prof.txt', comments='%')

    uv_budget = rste_k[:, 7]
    if np.nanmean(uv_budget) > 0.0:
        uv_budget = -uv_budget

    return {
        'yp': mean[:, 1],
        'Up': mean[:, 2],
        'urp': np.sqrt(np.clip(vel[:, 2], 0.0, None)),
        'vrp': np.sqrt(np.clip(vel[:, 3], 0.0, None)),
        'wrp': np.sqrt(np.clip(vel[:, 4], 0.0, None)),
        'uvp': -vel[:, 5],
        'TKEp': vel[:, 8],
        'prod_p': rste_k[:, 2],
        'eps_budget_p': uv_budget,
        'eps_p': -uv_budget,
    }


def load_mkm_reference(case_dir):
    means_path = case_dir / 'MKM_MEANS.dat'
    reystress_path = case_dir / 'MKM_REYSTRESS.dat'
    if not means_path.exists() or not reystress_path.exists():
        return None

    means = np.loadtxt(means_path, comments='#')
    reystress = np.loadtxt(reystress_path, comments='#')

    return {
        'yp': means[:, 1],
        'Up': means[:, 2],
        'urp': np.sqrt(np.clip(reystress[:, 2], 0.0, None)),
        'vrp': np.sqrt(np.clip(reystress[:, 3], 0.0, None)),
        'wrp': np.sqrt(np.clip(reystress[:, 4], 0.0, None)),
        'uvp': -reystress[:, 5],
        'TKEp': 0.5 * (reystress[:, 2] + reystress[:, 3] + reystress[:, 4]),
    }


def make_figure():
    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(2, 3, hspace=0.38, wspace=0.32)
    ax_U = fig.add_subplot(gs[0, 0])
    ax_rms = fig.add_subplot(gs[0, 1])
    ax_uv = fig.add_subplot(gs[0, 2])
    ax_TKE = fig.add_subplot(gs[1, 0])
    ax_bud = fig.add_subplot(gs[1, 1])
    ax_eps = fig.add_subplot(gs[1, 2])
    return fig, ax_U, ax_rms, ax_uv, ax_TKE, ax_bud, ax_eps


def plot_dns(ax_U, ax_rms, ax_uv, ax_TKE, ax_bud, ax_eps, dns):
    ax_U.semilogx(dns['yp'], dns['Up'], color='tab:blue', lw=1.6,
                  label='DNS')

    ax_rms.plot(dns['urp'], dns['yp'], color='tab:red', lw=1.6,
                label=r"DNS $u'_{rms}$")
    ax_rms.plot(dns['vrp'], dns['yp'], color='tab:green', lw=1.6,
                label=r"DNS $v'_{rms}$")
    ax_rms.plot(dns['wrp'], dns['yp'], color='tab:purple', lw=1.6,
                label=r"DNS $w'_{rms}$")

    ax_uv.plot(dns['uvp'], dns['yp'], color='tab:blue', lw=1.6,
               label='DNS')
    ax_TKE.plot(dns['TKEp'], dns['yp'], color='tab:orange', lw=1.6,
                label='DNS')

    ax_bud.plot(dns['prod_p'], dns['yp'], color='tab:red', lw=1.6,
                label=r'DNS $P^+$')
    ax_bud.plot(dns['eps_budget_p'], dns['yp'], color='tab:blue', lw=1.6,
                label=r'DNS $-\varepsilon^+$')

    ax_eps.plot(dns['eps_p'], dns['yp'], color='tab:blue', lw=1.6,
                label=r'DNS $\varepsilon^+$')


def plot_mkm(ax_U, ax_rms, ax_uv, ax_TKE, mkm):
    ax_U.semilogx(mkm['yp'], mkm['Up'], color='tab:brown', lw=1.3,
                  ls='--', label='MKM DNS')

    ax_rms.plot(mkm['urp'], mkm['yp'], color='tab:brown', lw=1.3,
                ls='--', label=r"MKM $u'_{rms}$")
    ax_rms.plot(mkm['vrp'], mkm['yp'], color='tab:brown', lw=1.3,
                ls='--', label=r"MKM $v'_{rms}$")
    ax_rms.plot(mkm['wrp'], mkm['yp'], color='tab:brown', lw=1.3,
                ls='--', label=r"MKM $w'_{rms}$")

    ax_uv.plot(mkm['uvp'], mkm['yp'], color='tab:brown', lw=1.3,
               ls='--', label='MKM DNS')
    ax_TKE.plot(mkm['TKEp'], mkm['yp'], color='tab:brown', lw=1.3,
                ls='--', label='MKM DNS')


def plot_les_case(ax_U, ax_rms, ax_uv, ax_TKE, ax_bud, ax_eps, h, label):
    ax_U.semilogx(h['yp'], h['Up'], label=label, **_LES_STYLE)

    ax_rms.plot(h['urp'], h['yp'], label=rf"{label} $u'_{{rms}}$", **_LES_STYLE)
    ax_rms.plot(h['vrp'], h['yp'], label=rf"{label} $v'_{{rms}}$", **_LES_STYLE)
    ax_rms.plot(h['wrp'], h['yp'], label=rf"{label} $w'_{{rms}}$", **_LES_STYLE)

    ax_uv.plot(h['uvp'], h['yp'], label=label, **_LES_STYLE)
    ax_TKE.plot(h['TKEp'], h['yp'], label=label, **_LES_STYLE)

    ax_bud.plot(h['prod_p'], h['yp'], label=rf'{label} $P^+$', **_LES_STYLE)
    ax_bud.plot(-h['eps_tot'], h['yp'],
                label=rf'{label} $-\varepsilon^+_{{tot}}$', **_LES_STYLE)

    ax_eps.plot(h['epsr_p'], h['yp'],
                label=rf'{label} $\varepsilon^+_{{res}}$', **_LES_STYLE)
    ax_eps.plot(h['epss_p'], h['yp'],
                label=rf'{label} $\varepsilon^+_{{SGS}}$', **_LES_STYLE)


def plot_profiles(npz_path, case_dir, nu, label, save_path, dpi=150):
    """Plot wall-unit profiles for one stats file, with DNS/MKM references."""
    save_path = Path(save_path)
    dns = None
    dns_dir = case_dir / 'dns_ref'
    dns_files = ['mean_prof.txt', 'vel_fluc_prof.txt', 'RSTE_k_prof.txt']
    if all((dns_dir / name).exists() for name in dns_files):
        dns = load_dns_reference(dns_dir)

    mkm = load_mkm_reference(case_dir)

    fig, ax_U, ax_rms, ax_uv, ax_TKE, ax_bud, ax_eps = make_figure()

    if dns is not None:
        plot_dns(ax_U, ax_rms, ax_uv, ax_TKE, ax_bud, ax_eps, dns)
    if mkm is not None:
        plot_mkm(ax_U, ax_rms, ax_uv, ax_TKE, mkm)

    data = np.load(npz_path)
    u_tau = 1.0
    delta_nu = nu / u_tau
    Re_tau = u_tau * (data['ym'].max() / 2.0) / nu
    print(f'{label}:  u_τ = {u_tau:.5f},  Re_τ = {Re_tau:.1f}')

    h = half_channel(data, u_tau, delta_nu, nu)
    plot_les_case(ax_U, ax_rms, ax_uv, ax_TKE, ax_bud, ax_eps, h, label)

    ax_U.set_xlabel(r'$y^+$')
    ax_U.set_ylabel(r'$U^+$')
    ax_U.set_title('Mean streamwise velocity')
    ax_U.legend(fontsize=7)

    ax_rms.set_xlabel(r'Fluctuation amplitude in wall units')
    ax_rms.set_ylabel(r'$y^+$')
    ax_rms.set_title('Velocity RMS')

    ax_uv.set_xlabel(r"$-\langle u'v' \rangle^+$")
    ax_uv.set_ylabel(r'$y^+$')
    ax_uv.set_title('Reynolds shear stress')

    ax_TKE.set_xlabel(r'$k^+$')
    ax_TKE.set_ylabel(r'$y^+$')
    ax_TKE.set_title('Turbulent kinetic energy')

    ax_bud.set_xlabel(r'$P^+,\;-\varepsilon^+$')
    ax_bud.set_ylabel(r'$y^+$')
    ax_bud.axvline(0, color='k', lw=0.6, ls=':')
    ax_bud.set_title('TKE production & dissipation')

    ax_eps.set_xlabel(r'$\varepsilon^+$')
    ax_eps.set_ylabel(r'$y^+$')
    ax_eps.set_title('Dissipation breakdown (res. / SGS)')

    ax_U.grid(True, which='both', ls=':', lw=0.4, alpha=0.6)
    ax_U.set_xlim(left=1)

    for ax in (ax_rms, ax_uv, ax_TKE, ax_bud, ax_eps):
        ax.grid(True, ls=':', lw=0.4, alpha=0.6)
        ax.set_ylim(bottom=1)

    fig.suptitle('Channel flow statistics — wall units', fontsize=12, y=1.01)

    save_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(save_path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f'Plot saved: {save_path}')


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    case_dir   = Path(__file__).parent
    fields_dir = case_dir / args.fields_dir
    stats_dir  = case_dir / args.stats_dir
    stats_dir.mkdir(exist_ok=True)

    # ── physical parameters ───────────────────────────────────────────────────
    params    = snapshot_io.parse_input_parameters(case_dir / 'input_parameters')
    nu        = args.nu if args.nu is not None else float(params['nu'])
    sgs_model = 0 if args.no_sgs else int(params.get('sgs_model', 0))
    print(f'ν = {nu:.6e},  sgs_model = {sgs_model}')

    # ── build snapshot list ───────────────────────────────────────────────────
    steps     = range(args.avg_start, args.avg_end + 1, args.interval)
    snapfiles = []
    for step in steps:
        fpath = fields_dir / f'{args.prefix}.{step}'
        if fpath.is_file():
            snapfiles.append((step, fpath))
        else:
            print(f'  WARNING: {fpath.name} not found, skipping.')

    if not snapfiles:
        sys.exit('ERROR: No snapshot files found in the specified range.')

    print(f'Found {len(snapfiles)} snapshots: '
          f'steps {snapfiles[0][0]} – {snapfiles[-1][0]}')

    # ── probe grid from first snapshot ───────────────────────────────────────
    print(f'Reading grid from {snapfiles[0][1].name} …')
    s0  = snapshot_io.read_snapshot(snapfiles[0][1], sgs_model=sgs_model)
    xm  = s0['xm'];  ym  = s0['ym'];  zm  = s0['zm']
    nxm = xm.size;   nym = ym.size;   nzm = zm.size
    dx  = xm[1] - xm[0]   # uniform x
    dz  = zm[1] - zm[0]   # uniform z
    Ly  = float(s0['y'][-1])
    print(f'Grid : nxm={nxm}, nym={nym}, nzm={nzm}')
    print(f'Domain: Ly={Ly:.4f},  dx={dx:.6f},  dz={dz:.6f}')

    # ── accumulation buffers (y-profiles) ────────────────────────────────────
    acc_U   = np.zeros(nym)
    acc_V   = np.zeros(nym)
    acc_W   = np.zeros(nym)
    acc_UU  = np.zeros(nym)
    acc_VV  = np.zeros(nym)
    acc_WW  = np.zeros(nym)
    acc_UV  = np.zeros(nym)

    # Resolved pseudo-dissipation: accumulate Σᵢⱼ <(∂uᵢ/∂xⱼ)²>
    # The mean-gradient correction is applied after all snapshots are read.
    acc_grad_sq = np.zeros(nym)

    # SGS dissipation: accumulate <νₜ · 2SᵢⱼSᵢⱼ>
    acc_eps_sgs = np.zeros(nym)
    has_sgs     = False

    N = 0  # snapshot counter

    # ── main averaging loop ───────────────────────────────────────────────────
    for step, fpath in snapfiles:
        print(f'  [{N+1:>4d}/{len(snapfiles)}]  step {step}', end='\r', flush=True)

        snap = snapshot_io.read_snapshot(fpath, sgs_model=sgs_model)
        U    = snap['U']   # (nxm, nym, nzm)
        V    = snap['V']
        W    = snap['W']

        # First moments
        acc_U  += xz_mean(U)
        acc_V  += xz_mean(V)
        acc_W  += xz_mean(W)

        # Second moments (for variances and Reynolds stress)
        acc_UU += xz_mean(U * U)
        acc_VV += xz_mean(V * V)
        acc_WW += xz_mean(W * W)
        acc_UV += xz_mean(U * V)

        # Velocity gradients
        g = grad_all(U, V, W, ym, dx, dz)

        # Sum of squared velocity gradients (resolved pseudo-dissipation)
        grad_sq = (g['dU_dx']**2 + g['dU_dy']**2 + g['dU_dz']**2 +
                   g['dV_dx']**2 + g['dV_dy']**2 + g['dV_dz']**2 +
                   g['dW_dx']**2 + g['dW_dy']**2 + g['dW_dz']**2)
        acc_grad_sq += xz_mean(grad_sq)

        # SGS dissipation
        if sgs_model != 0 and 'nu_t' in snap:
            has_sgs = True
            nu_t = snap['nu_t']
            # Strain-rate tensor components
            S_12 = 0.5 * (g['dU_dy'] + g['dV_dx'])
            S_13 = 0.5 * (g['dU_dz'] + g['dW_dx'])
            S_23 = 0.5 * (g['dV_dz'] + g['dW_dy'])
            two_SijSij = (2.0 * (g['dU_dx']**2 + g['dV_dy']**2 + g['dW_dz']**2) +
                          4.0 * (S_12**2 + S_13**2 + S_23**2))
            acc_eps_sgs += xz_mean(nu_t * two_SijSij)

        N += 1

    print(f'\nProcessed {N} snapshot(s).')

    # ── time averages ─────────────────────────────────────────────────────────
    U_mean = acc_U  / N
    V_mean = acc_V  / N
    W_mean = acc_W  / N

    UU_mean = acc_UU / N
    VV_mean = acc_VV / N
    WW_mean = acc_WW / N
    UV_mean = acc_UV / N

    # Reynolds stress components  <u_i'u_j'> = <u_i u_j> - <u_i><u_j>
    uu = np.maximum(UU_mean - U_mean**2, 0.0)
    vv = np.maximum(VV_mean - V_mean**2, 0.0)
    ww = np.maximum(WW_mean - W_mean**2, 0.0)
    uv = UV_mean - U_mean * V_mean

    urms = np.sqrt(uu)
    vrms = np.sqrt(vv)
    wrms = np.sqrt(ww)
    TKE  = 0.5 * (uu + vv + ww)

    # Resolved pseudo-dissipation with mean-shear correction:
    # (only y-derivatives of the mean are non-zero for a fully-developed channel)
    dU_dy_mean = np.gradient(U_mean, ym)
    dV_dy_mean = np.gradient(V_mean, ym)
    dW_dy_mean = np.gradient(W_mean, ym)

    eps_res = nu * (acc_grad_sq / N
                    - dU_dy_mean**2
                    - dV_dy_mean**2
                    - dW_dy_mean**2)

    eps_sgs = acc_eps_sgs / N if has_sgs else np.zeros(nym)

    # TKE production  P = -<u'v'> · dU/dy
    prod = -uv * dU_dy_mean

    # ── channel symmetry folding about y = Ly/2 ───────────────────────────────
    # Symmetric quantities (even about centreline):
    U_mean_s  = sym_fold(U_mean)
    V_mean_s  = sym_fold(V_mean)
    W_mean_s  = sym_fold(W_mean)
    uu_s      = sym_fold(uu)
    vv_s      = sym_fold(vv)
    ww_s      = sym_fold(ww)
    TKE_s     = sym_fold(TKE)
    eps_res_s = sym_fold(eps_res)
    eps_sgs_s = sym_fold(eps_sgs)

    # Anti-symmetric quantity (odd about centreline):
    uv_s = antisym_fold(uv)

    urms_s = np.sqrt(np.maximum(uu_s, 0.0))
    vrms_s = np.sqrt(np.maximum(vv_s, 0.0))
    wrms_s = np.sqrt(np.maximum(ww_s, 0.0))

    # Recompute production using symmetrised profiles
    dU_dy_s = np.gradient(U_mean_s, ym)
    prod_s  = -uv_s * dU_dy_s

    # ── wall-unit summary ─────────────────────────────────────────────────────
    u_tau   = np.sqrt(nu * np.abs(dU_dy_s[0]))
    Re_tau  = u_tau * (Ly / 2.0) / nu
    print(f'\nEstimated wall units:')
    print(f'  u_τ   = {u_tau:.6f}')
    print(f'  Re_τ  = {Re_tau:.2f}')

    # ── save results ──────────────────────────────────────────────────────────
    results = {
        'ym'      : ym,
        'U_mean'  : U_mean_s,
        'V_mean'  : V_mean_s,
        'W_mean'  : W_mean_s,
        'urms'    : urms_s,
        'vrms'    : vrms_s,
        'wrms'    : wrms_s,
        'uv'      : uv_s,
        'TKE'     : TKE_s,
        'eps_res' : eps_res_s,
        'eps_sgs' : eps_sgs_s,
        'prod'    : prod_s,
    }

    # Compressed NumPy archive (load with np.load)
    npz_path = stats_dir / f'stats_{args.prefix}_{args.avg_start}_{args.avg_end}.npz'
    np.savez_compressed(npz_path, **results)
    print(f'\nSaved: {npz_path}')

    # Individual little-endian float64 binary files
    print('Saved binary profiles:')
    for name, arr in results.items():
        bin_path = stats_dir / f'{name}.bin'
        arr.astype('<f8').tofile(bin_path)
        print(f'  {bin_path}')

    # ── automatic plotting ────────────────────────────────────────────────────
    plot_path = Path(args.out_fig)
    if not plot_path.is_absolute():
        plot_path = case_dir / plot_path
    plot_profiles(npz_path, case_dir, nu, args.prefix, plot_path)


if __name__ == '__main__':
    main()
