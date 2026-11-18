#!/usr/bin/env python3
"""hit_spectral_reference.py -- independent pseudo-spectral DNS of forced homogeneous isotropic turbulence, the reference for
the finite-difference solver's HIT validation (src/hit.f90).

Same problem as the solver: periodic box of side 2*pi with N^3 points, viscosity nu, band-limited deterministic forcing
f_hat = eps_in/(2 E_band) u_hat on the solenoidal modes 0 < |n| <= kf (Lamorgese, Caughey & Pope 2005), random-phase
solenoidal Passot-Pouquet initial field with rms u0 and peak kp. Numerics: rotational form of the convection, 2/3 truncation,
integrating-factor RK4, scipy.fft with worker threads.

Output (in --out): hit_stats.dat (columns as the solver's: step t TKE eps u_rms lambda Re_lambda eta eta/dx L T_e),
spectrum.<k>.dat snapshots (k, E(k) over unit shells, D(k)) every --tspec, and fields.<n>.npz of the final state (--save).
Usage: hit_spectral_reference.py --N 128 --nu 0.005 --eps 0.1 --T 80 --out DIR [--workers 4]
"""
import argparse
import os
import time

import numpy as np
import scipy.fft as sfft


def build_grid(N):
    n = np.fft.fftfreq(N, 1.0 / N)
    nr = np.fft.rfftfreq(N, 1.0 / N)
    kx, ky, kz = np.meshgrid(n, n, nr, indexing="ij")
    return kx, ky, kz


def shell_spectrum(uh, vh, wh, kmag, N, wgt):
    """E(k) on unit shells (k = 1..N/2) from half-complex arrays (wgt = 2 for the doubled kz>0 planes, 1 for kz=0, Nyquist)."""
    e = 0.5 * wgt * (np.abs(uh) ** 2 + np.abs(vh) ** 2 + np.abs(wh) ** 2) / N ** 6
    idx = np.rint(kmag).astype(int).ravel()
    E = np.bincount(idx, weights=e.ravel(), minlength=N)
    return E[: N // 2 + 1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=128)
    ap.add_argument("--nu", type=float, default=0.005)
    ap.add_argument("--eps", type=float, default=0.1, help="forcing power eps_in")
    ap.add_argument("--kf", type=float, default=2.5)
    ap.add_argument("--u0", type=float, default=0.8, help="initial rms velocity per component")
    ap.add_argument("--kp", type=float, default=3.0)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--T", type=float, default=60.0)
    ap.add_argument("--dt", type=float, default=0.0, help="0: CFL-adaptive (cfl 0.5)")
    ap.add_argument("--cfl", type=float, default=0.5)
    ap.add_argument("--tstat", type=float, default=0.5, help="statistics output interval")
    ap.add_argument("--tspec", type=float, default=5.0, help="spectrum output interval")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--save", action="store_true", help="save the final velocity field")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    os.makedirs(a.out, exist_ok=True)
    N, nu = a.N, a.nu
    W = a.workers
    kx, ky, kz = build_grid(N)
    k2 = kx ** 2 + ky ** 2 + kz ** 2
    kmag = np.sqrt(k2)
    k2s = np.where(k2 == 0, 1.0, k2)
    dealias = (kmag <= N / 3.0)
    wgt = np.where((kz == 0) | (kz == N // 2), 1.0, 2.0)
    band = (k2 > 0) & (k2 <= a.kf ** 2)

    def ifft(fh):
        return sfft.irfftn(fh, s=(N, N, N), workers=W)

    def fft(f):
        return sfft.rfftn(f, workers=W)

    def project(uh, vh, wh):
        d = (kx * uh + ky * vh + kz * wh) / k2s
        return uh - kx * d, vh - ky * d, wh - kz * d

    # random solenoidal IC (Passot-Pouquet), Hermitian by construction via a real random field in physical space
    rng = np.random.default_rng(a.seed)
    uh = [fft(rng.standard_normal((N, N, N))) for _ in range(3)]
    uh = list(project(*uh))
    E_target = np.where(kmag > 0, kmag ** 4 * np.exp(-2.0 * k2 / a.kp ** 2) / (4 * np.pi * np.maximum(k2, 1e-30)), 0.0)
    for c in range(3):
        uh[c] = uh[c] * np.sqrt(E_target) * dealias
    uh = list(project(*uh))
    ke = (0.5 * wgt * sum(np.abs(c) ** 2 for c in uh) / N ** 6).sum()
    s = np.sqrt(1.5 * a.u0 ** 2 / ke)
    uh = [c * s for c in uh]

    def rhs(uh):
        """-P[u x omega] + forcing (physical-space product, 2/3 truncation)."""
        u = [ifft(c) for c in uh]
        om = (ifft(1j * (ky * uh[2] - kz * uh[1])), ifft(1j * (kz * uh[0] - kx * uh[2])),
              ifft(1j * (kx * uh[1] - ky * uh[0])))
        nl = (fft(u[1] * om[2] - u[2] * om[1]), fft(u[2] * om[0] - u[0] * om[2]), fft(u[0] * om[1] - u[1] * om[0]))
        nl = project(*[c * dealias for c in nl])
        Eb = (0.5 * wgt * band * sum(np.abs(c) ** 2 for c in uh) / N ** 6).sum()
        sc = a.eps / (2.0 * Eb) if Eb > 1e-300 else 0.0
        return [nl[i] + sc * band * uh[i] for i in range(3)], u

    tag = open(os.path.join(a.out, "hit_stats.dat"), "w")
    tag.write("# step t TKE eps u_rms lambda Re_lambda eta eta/dx L=u^3/eps T_e=L/u_rms\n")
    t, step, next_stat, next_spec, wall0 = 0.0, 0, 0.0, 0.0, time.time()
    dx = 2 * np.pi / N
    nspec = 0
    umax = a.u0 * 3
    while t < a.T + 1e-12:
        if t >= next_stat - 1e-12:
            e = 0.5 * wgt * sum(np.abs(c) ** 2 for c in uh) / N ** 6
            tke = e.sum()
            eps = 2 * nu * (k2 * e).sum()
            ur = np.sqrt(2 * tke / 3)
            lam = np.sqrt(15 * nu * ur ** 2 / eps)
            eta = (nu ** 3 / eps) ** 0.25
            L = ur ** 3 / eps
            tag.write(f"{step:10d}{t:15.7E}{tke:15.7E}{eps:15.7E}{ur:15.7E}{lam:15.7E}{ur * lam / nu:15.7E}{eta:15.7E}"
                      f"{eta / dx:15.7E}{L:15.7E}{L / ur:15.7E}\n")
            tag.flush()
            next_stat += a.tstat
            print(f"t={t:8.3f} TKE={tke:.4f} eps={eps:.4f} Re_l={ur * lam / nu:.1f} kmax*eta={N / 3 * eta:.2f} "
                  f"umax={umax:.2f} wall={time.time() - wall0:.0f}s", flush=True)
        if t >= next_spec - 1e-12:
            E = shell_spectrum(*uh, kmag, N, wgt)
            k = np.arange(len(E))
            np.savetxt(os.path.join(a.out, f"spectrum.{nspec:04d}.dat"), np.c_[k, E, 2 * nu * k ** 2 * E],
                       header=f"t = {t:.6f}  k E(k) D(k)")
            nspec += 1
            next_spec += a.tspec
        dt = a.dt if a.dt > 0 else min(a.cfl * dx / umax, 0.02)
        # integrating-factor RK4: the viscous term is exact, E = exp(-nu k^2 dt/2)
        E1 = np.exp(-nu * k2 * dt / 2.0)
        E2 = E1 * E1
        k1, u = rhs(uh)
        umax = max(np.abs(c).max() for c in u) * np.sqrt(3.0)
        k2_, _ = rhs([E1 * (uh[i] + 0.5 * dt * k1[i]) for i in range(3)])
        k3_, _ = rhs([E1 * uh[i] + 0.5 * dt * k2_[i] for i in range(3)])
        k4_, _ = rhs([E2 * uh[i] + dt * E1 * k3_[i] for i in range(3)])
        uh = [E2 * uh[i] + dt / 6.0 * (E2 * k1[i] + 2.0 * E1 * (k2_[i] + k3_[i]) + k4_[i]) for i in range(3)]
        t += dt
        step += 1
    tag.close()
    if a.save:
        u = [ifft(c) for c in uh]
        np.savez(os.path.join(a.out, "fields.final.npz"), u=u[0], v=u[1], w=u[2], t=t)


if __name__ == "__main__":
    main()
