"""hit.py -- statistics of forced homogeneous isotropic turbulence (src/hit.f90, triply periodic box).

read_stats        hit_stats.dat (monitor-cadence TKE, eps, Re_lambda, ...) and its time average over a window
read_velocity     U,V,W of a full snapshot on the N^3 lattices of distinct points, in the solver's own staggering
spectrum          shell-summed E(k) (unit shells of the fundamental wavenumber) of the staggered components, and the
                  modified-wavenumber dissipation spectrum the finite-difference solver actually applies
scales            u', eps, eta, lambda, Re_lambda, L11 and the velocity-derivative skewness/flatness of one field
pope_spectrum     Pope (2000, eq. 6.246) model spectrum, the universal-shape reference
read_spectral     the output directory of postProcessing/validation/hit_spectral_reference.py
Periodic geometry: N = n-2 distinct points per direction (n faces), period P = N*dx (see src/hit.f90). Requires numpy.
"""
import glob
import re
from pathlib import Path

import numpy as np

from . import _core


# ════════════════════════════════════════════════════════════════════════════
# 1. DATA
# ════════════════════════════════════════════════════════════════════════════
def read_stats(path):
    """hit_stats.dat -> dict of columns: step t tke eps urms lam rel eta eta_dx L Te."""
    d = np.loadtxt(path)
    names = "step t tke eps urms lam rel eta eta_dx L Te".split()
    return {n: d[:, i] for i, n in enumerate(names)}


def time_average(stats, t0, t1=np.inf):
    """Mean of every column over t0 <= t <= t1 and the standard error of the mean from 10 blocks."""
    m = (stats["t"] >= t0) & (stats["t"] <= t1)
    out = {}
    for k, v in stats.items():
        if k in ("step", "t"):
            continue
        x = v[m]
        blocks = np.array_split(x, 10) if len(x) >= 20 else [x]
        out[k] = (float(x.mean()), float(np.std([b.mean() for b in blocks]) / np.sqrt(len(blocks))))
    out["n"] = int(m.sum())
    return out


def read_velocity(path):
    """Full snapshot -> dict t-less: U,V,W (N,N,N) on their own staggered lattices of distinct points, dx, N."""
    r = _core.BinaryReader(data=Path(path).read_bytes())
    nx = r.ri(); x = r.rd(nx)
    ny = r.ri(); r.rd(ny)
    nz = r.ri(); r.rd(nz)
    for _ in range(3):
        n = r.ri(); r.rd(n)
    _, U = r.read_field()
    _, V = r.read_field()
    _, W = r.read_field()
    Un = U[1:-1, 1:-2, 1:-2]
    Vn = V[1:-2, 1:-1, 1:-2]
    Wn = W[1:-2, 1:-2, 1:-1]
    N = Un.shape[0]
    if not (Un.shape == Vn.shape == Wn.shape == (N, N, N)):
        raise ValueError(f"{path}: not a cubic periodic box ({Un.shape}, {Vn.shape}, {Wn.shape})")
    return dict(U=Un, V=Vn, W=Wn, dx=float(x[1] - x[0]), N=N)


def read_spectral(case_dir):
    """Spectral-reference directory -> (stats dict, list of (t, k, E, D))."""
    case_dir = Path(case_dir)
    st = read_stats(case_dir / "hit_stats.dat")
    spec = []
    for f in sorted(glob.glob(str(case_dir / "spectrum.*.dat"))):
        t = float(re.search(r"t = ([0-9.eE+-]+)", open(f).readline()).group(1))
        a = np.loadtxt(f)
        spec.append((t, a[:, 0], a[:, 1], a[:, 2]))
    return st, spec


# ════════════════════════════════════════════════════════════════════════════
# 2. ANALYSIS
# ════════════════════════════════════════════════════════════════════════════
def spectrum(vel, nu):
    """E(k), the physical dissipation 2 nu k^2 E and the discrete one 2 nu k_mod^2 E on unit shells, k = 1..N/2.

    E(k) dk is the energy of the modes in the shell (|n| rounded to k) with k in units of 2*pi/P; the returned wavenumber
    is the physical k = 2*pi*n/P. k_mod = (2/dx) sin(k dx/2) is what the staggered second-order operators see.
    """
    N, dx = vel["N"], vel["dx"]
    P = N * dx
    n = np.fft.fftfreq(N, 1.0 / N)
    nr = np.fft.rfftfreq(N, 1.0 / N)
    nx_, ny_, nz_ = np.meshgrid(n, n, nr, indexing="ij")
    nmag = np.sqrt(nx_ ** 2 + ny_ ** 2 + nz_ ** 2)
    wgt = np.where((nz_ == 0) | (nz_ == N // 2), 1.0, 2.0)
    e = np.zeros_like(nmag)
    for c in "UVW":
        e += 0.5 * wgt * np.abs(np.fft.rfftn(vel[c])) ** 2 / N ** 6
    shell = np.rint(nmag).astype(int).ravel()
    ks = np.arange(N // 2 + 1)
    kphys = 2 * np.pi * ks / P
    E = np.bincount(shell, weights=e.ravel(), minlength=N)[: N // 2 + 1] * (P / (2 * np.pi))   # per unit physical k
    # exact dissipation sums with the modified wavenumber of each mode
    kmod2 = sum((2.0 / dx * np.sin(np.pi * m * dx / P)) ** 2 for m in (nx_, ny_, nz_))
    Dd = np.bincount(shell, weights=(2 * nu * kmod2 * e).ravel(), minlength=N)[: N // 2 + 1]
    return kphys, E, 2 * nu * kphys ** 2 * E, Dd / (2 * np.pi / P)


def scales(k, E, nu):
    """Integral scales from a spectrum E(k) (per unit k, k[0] = 0 shell ignored)."""
    dk = k[1] - k[0]
    K = np.sum(E[1:]) * dk
    eps = 2 * nu * np.sum(k[1:] ** 2 * E[1:]) * dk
    urms = np.sqrt(2 * K / 3)
    lam = np.sqrt(15 * nu * urms ** 2 / eps)
    L11 = np.pi / (2 * urms ** 2) * np.sum(E[1:] / k[1:]) * dk
    return dict(tke=K, eps=eps, urms=urms, lam=lam, rel=urms * lam / nu, eta=(nu ** 3 / eps) ** 0.25, L11=L11,
                L=K ** 1.5 / eps)


def gradient_moments(vel):
    """Skewness and flatness of the longitudinal velocity gradients (mean of the three), from the staggered differences."""
    dx = vel["dx"]
    g = [(np.roll(vel[c], -1, axis=a) - vel[c]) / dx for a, c in enumerate("UVW")]
    # forward difference of the staggered component along its own axis lives at the cell centre: that is du_i/dx_i
    sk = np.mean([np.mean(x ** 3) / np.mean(x ** 2) ** 1.5 for x in g])
    fl = np.mean([np.mean(x ** 4) / np.mean(x ** 2) ** 2 for x in g])
    return float(sk), float(fl)


def pope_spectrum(k, eps, nu, tke, c=1.5, beta=5.2, c_eta=0.40, c_l=6.78, p0=2.0):
    """Pope (2000) eq. (6.246) model spectrum E(k) = C eps^(2/3) k^(-5/3) f_L(kL) f_eta(k eta), L = tke^(3/2)/eps."""
    L = tke ** 1.5 / eps
    eta = (nu ** 3 / eps) ** 0.25
    fl = (k * L / np.sqrt((k * L) ** 2 + c_l)) ** (5.0 / 3.0 + p0)
    fe = np.exp(-beta * (((k * eta) ** 4 + c_eta ** 4) ** 0.25 - c_eta))
    return c * eps ** (2.0 / 3.0) * k ** (-5.0 / 3.0) * fl * fe
