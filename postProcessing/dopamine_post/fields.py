"""
fields.py -- post-processing of fdm-dopamine binary field (U,V,W,P,...) snapshots.

One module, three layers (use the library from your own scripts, or the future CLI):

  1. DATA      FieldSnapshot   one snapshot's U,V,W,P,(C,nu_t,T) + grid, cell-centred,
                                ghost layers stripped.
               FieldSeries     a case's fields/ directory: list + iterate snapshots.
  2. ANALYSIS  FieldSeries.time_average     time- and (x,z)-plane averaged channel
                                            statistics with centreline symmetry folding
                                            (compute_stats.py).
               FieldSeries.profile_vs_dns   spanwise+snapshot-averaged profiles at one
                                            or more x-stations, with an optional DNS
                                            overlay (analyse_channel.py).
  3. FIGURES   plot_profile, plot_stats, plot_profile_vs_dns, animate_velocity_slice
               write_field_xmf / FieldSeries.write_xmf   (generateXMF.py)

Binary layout (big-endian float64, Fortran stream, no record markers) -- see
_core.BinaryReader.read_field for the low-level block reader:

  Grid header -- 6 blocks of (Int32 count, float64[count] coordinates):
      nx  face points in x  ->  x        nxm cell centres in x -> xm
      ny  face points in y  ->  y        nym cell centres in y -> ym
      nz  face points in z  ->  z        nzm cell centres in z -> zm

  Field blocks -- each: (3 x Int32 shape), (float64[n1*n2*n3] data, Fortran order):
      U    (nx,  nyg, nzg)   x-face velocity         -- always present
      V    (nxg, ny,  nzg)   y-face velocity         -- always present
      W    (nxg, nyg, nz )   z-face velocity         -- always present
      P    (nxg, nyg, nzg)   cell-centre pressure    -- always present
      C    (nxg, nyg, nzg)   scalar concentration    -- only if sediment_flag >= 1
      nu_t (nxg, nyg, nzg)   SGS turbulent viscosity -- only if sgs_model != 0
      T    (nxg, nyg, nzg)   Boussinesq temperature  -- only if boussinesq_flag >= 1
  (block order C, nu_t, T after P matches input_output.f90's snapshot writer exactly.)

  Ghost-cell convention: nxg = nxm+2, nyg = nym+2, nzg = nzm+2.  ny and nz in the V
  and W headers are face-point counts (no ghost cells on the staggered direction),
  so ny = nym+1, nz = nzm+1.  Returned arrays are cell-centred with every ghost layer
  removed: shape (nxm, nym, nzm).  Axis convention: 0 = x (streamwise),
  1 = y (wall-normal), 2 = z (spanwise).

Library use
-----------
    from dopamine_post.fields import FieldSnapshot, FieldSeries, plot_profile

    series = FieldSeries(".")
    snap = series.latest()
    ym, U = snap.profile("U")                       # x-z averaged mean profile
    stats = series.time_average((100000, 150000), 500)
    stats.save("stats")
    result = series.profile_vs_dns("dns_ref", x_stations=[2, 4, 6, 8])
    series.write_xmf()

Requires: numpy, matplotlib (figures/movies); ffmpeg for animate_velocity_slice.
"""
import os
import re
import struct
import warnings
from pathlib import Path

import numpy as np

from . import _core
from ._core import list_indexed, render_movie

_PAT = re.compile(r"^([A-Za-z0-9_]+)\.(\d+)$")
_AXIS = {"x": 0, "y": 1, "z": 2}


# ════════════════════════════════════════════════════════════════════════════
# 1. DATA
# ════════════════════════════════════════════════════════════════════════════
def _to_cell_centre(U_face, V_face, W_face, P_raw):
    """Interpolate staggered face values to cell centres and strip ghost layers.

        U_face : (nx,  nyg, nzg)   x-face velocity
        V_face : (nxg, ny,  nzg)   y-face velocity   (ny = nym+1 face points)
        W_face : (nxg, nyg, nz )   z-face velocity   (nz = nzm+1 face points)
        P_raw  : (nxg, nyg, nzg)   cell-centre pressure, ghost ring on all sides

    Returns U, V, W, P each of shape (nxm, nym, nzm).
    """
    U = 0.5 * (U_face[:-1, 1:-1, 1:-1] + U_face[1:, 1:-1, 1:-1])
    V = 0.5 * (V_face[1:-1, :-1, 1:-1] + V_face[1:-1, 1:, 1:-1])
    W = 0.5 * (W_face[1:-1, 1:-1, :-1] + W_face[1:-1, 1:-1, 1:])
    P = P_raw[1:-1, 1:-1, 1:-1]
    return U, V, W, P


def _step_from_name(name):
    m = _PAT.match(name)
    return int(m.group(2)) if m else None


def _detect_prefix(fields_dir):
    for p in Path(fields_dir).iterdir():
        m = _PAT.match(p.name)
        if m:
            return m.group(1)
    return None


class FieldSnapshot:
    """One snapshot: U,V,W,P,(C,nu_t,T) cell-centred fields + grid, ghost layers
    stripped (snapshot_io.py's binary layout knowledge)."""

    def __init__(self, x, y, z, xm, ym, zm, fields, path=None, step=None,
                 sgs_model=0, sediment_flag=0, boussinesq_flag=0):
        self.x, self.y, self.z = x, y, z
        self.xm, self.ym, self.zm = xm, ym, zm
        self.fields = fields
        self.path = Path(path) if path is not None else None
        self.step = step
        self.sgs_model, self.sediment_flag, self.boussinesq_flag = (
            sgs_model, sediment_flag, boussinesq_flag)

    @classmethod
    def read(cls, path, sgs_model=None, sediment_flag=None, boussinesq_flag=None):
        """Read one binary snapshot.

        sgs_model/sediment_flag/boussinesq_flag select which optional field blocks
        (nu_t, C, T) are present; any left as None are looked up from
        `<case>/input_parameters` (case = the fields/ directory's parent), matching
        the run that produced the file.  Pass them explicitly to override, or when
        no input_parameters is available next to the snapshot.
        """
        path = Path(path)
        if sgs_model is None or sediment_flag is None or boussinesq_flag is None:
            prm = _core.read_input(path.parent.parent,
                                    ("sgs_model", "sediment_flag", "boussinesq_flag"))
            if sgs_model is None:
                sgs_model = int(prm.get("sgs_model", 0))
            if sediment_flag is None:
                sediment_flag = int(prm.get("sediment_flag", 0))
            if boussinesq_flag is None:
                boussinesq_flag = int(prm.get("boussinesq_flag", 0))

        r = _core.BinaryReader(data=path.read_bytes())

        nx = r.ri(); x = r.rd(nx)
        ny = r.ri(); y = r.rd(ny)
        nz = r.ri(); z = r.rd(nz)
        nxm = r.ri(); xm = r.rd(nxm)
        nym = r.ri(); ym = r.rd(nym)
        nzm = r.ri(); zm = r.rd(nzm)

        _, U_face = r.read_field()
        _, V_face = r.read_field()
        _, W_face = r.read_field()
        _, P_raw = r.read_field()
        U, V, W, P = _to_cell_centre(U_face, V_face, W_face, P_raw)
        fields = dict(U=U, V=V, W=W, P=P)

        if sediment_flag >= 1:
            _, C_raw = r.read_field()
            fields["C"] = C_raw[1:-1, 1:-1, 1:-1]
        if sgs_model != 0:
            _, nut_raw = r.read_field()
            fields["nu_t"] = nut_raw[1:-1, 1:-1, 1:-1]
        if boussinesq_flag >= 1:
            _, T_raw = r.read_field()
            fields["T"] = T_raw[1:-1, 1:-1, 1:-1]

        return cls(x, y, z, xm, ym, zm, fields, path=path, step=_step_from_name(path.name),
                   sgs_model=sgs_model, sediment_flag=sediment_flag,
                   boussinesq_flag=boussinesq_flag)

    def __repr__(self):
        return (f"FieldSnapshot(step={self.step}, vars={sorted(self.fields)}, "
                f"shape={self.fields['U'].shape})")

    def __getitem__(self, var):
        return self.fields[var]

    def __contains__(self, var):
        return var in self.fields

    def profile(self, var, average=("x", "z")):
        """x-z (or other pair of axes) averaged 1-D profile of `var`.  Returns
        (coord, values) along the one axis not averaged over (plot_snapshot.py logic)."""
        axes = tuple(sorted(_AXIS[a] for a in average))
        remaining = [a for a in (0, 1, 2) if a not in axes]
        if len(remaining) != 1:
            raise ValueError("profile() must average over exactly 2 of the 3 axes")
        coord = (self.xm, self.ym, self.zm)[remaining[0]]
        return coord, self.fields[var].mean(axis=axes)

    def plane(self, var, axis, index):
        """2-D slice of `var` at cell-centre index `index` along `axis` ('x'|'y'|'z')."""
        return np.take(self.fields[var], index, axis=_AXIS[axis])


class FieldSeries:
    """A case's fields/ directory: list + iterate binary snapshots."""

    def __init__(self, case_dir=".", prefix=None):
        self.case_dir = Path(case_dir)
        self.fields_dir = self.case_dir / "fields"
        self.prefix = prefix
        self._snaps = None

    @property
    def snapshots(self):
        """Sorted [(step, Path)], auto-detecting the prefix if not given at construction."""
        if self._snaps is None:
            self._snaps = list_indexed(self.fields_dir, _PAT, self.prefix)
        return self._snaps

    def steps(self):
        return [s for s, _ in self.snapshots]

    def __len__(self):
        return len(self.snapshots)

    def __iter__(self):
        """Yield (step, FieldSnapshot) for every snapshot, in step order."""
        for step, path in self.snapshots:
            yield step, FieldSnapshot.read(path)

    def read(self, step):
        """FieldSnapshot at an exact step."""
        match = [p for s, p in self.snapshots if s == step]
        if not match:
            raise FileNotFoundError(f"step {step} not found in {self.fields_dir}")
        return FieldSnapshot.read(match[0])

    def latest(self):
        """FieldSnapshot of the most recent step."""
        return FieldSnapshot.read(self.snapshots[-1][1])

    # ════════════════════════════════════════════════════════════════════════
    # 2. ANALYSIS
    # ════════════════════════════════════════════════════════════════════════
    def time_average(self, step_range, interval, nu=None, no_sgs=False, no_sediment=False,
                      no_boussinesq=False):
        """Time- and (x,z)-plane averaged channel statistics over
        `step_range = (start, end)` at the given snapshot `interval` (compute_stats.py's
        averaging loop, unchanged numerics), with channel-symmetry folding about
        y = Ly/2 applied to the momentum statistics.  Returns a FieldStats."""
        start, end = step_range
        ip = self.case_dir / "input_parameters"
        prm = _core.parse_input_parameters(ip) if ip.exists() else {}
        nu = nu if nu is not None else float(prm["nu"])
        sgs_model = 0 if no_sgs else int(prm.get("sgs_model", 0))
        sediment_flag = 0 if no_sediment else int(prm.get("sediment_flag", 0))
        boussinesq_flag = 0 if no_boussinesq else int(prm.get("boussinesq_flag", 0))

        prefix = self.prefix or _detect_prefix(self.fields_dir)
        steps = range(start, end + 1, interval)
        snapfiles = [(s, self.fields_dir / f"{prefix}.{s}") for s in steps]
        snapfiles = [(s, p) for s, p in snapfiles if p.is_file()]
        if not snapfiles:
            raise FileNotFoundError(
                f"no '{prefix}.<step>' snapshots in steps {start}..{end} "
                f"(interval {interval}) under {self.fields_dir}")

        read_kw = dict(sgs_model=sgs_model, sediment_flag=sediment_flag,
                       boussinesq_flag=boussinesq_flag)

        s0 = FieldSnapshot.read(snapfiles[0][1], **read_kw)
        xm, ym, zm = s0.xm, s0.ym, s0.zm
        nym = ym.size
        dx, dz = xm[1] - xm[0], zm[1] - zm[0]      # uniform x, z
        Ly = float(s0.y[-1])

        acc_U = np.zeros(nym); acc_V = np.zeros(nym); acc_W = np.zeros(nym)
        acc_UU = np.zeros(nym); acc_VV = np.zeros(nym); acc_WW = np.zeros(nym)
        acc_UV = np.zeros(nym)
        acc_grad_sq = np.zeros(nym)
        acc_eps_sgs = np.zeros(nym)
        has_sgs = False
        acc_T = np.zeros(nym); acc_TT = np.zeros(nym)
        has_T = False
        N = 0

        for step, fpath in snapfiles:
            snap = s0 if fpath == snapfiles[0][1] else FieldSnapshot.read(fpath, **read_kw)
            U, V, W = snap["U"], snap["V"], snap["W"]

            acc_U += _xz_mean(U); acc_V += _xz_mean(V); acc_W += _xz_mean(W)
            acc_UU += _xz_mean(U * U); acc_VV += _xz_mean(V * V); acc_WW += _xz_mean(W * W)
            acc_UV += _xz_mean(U * V)

            g = _grad_all(U, V, W, ym, dx, dz)
            grad_sq = (g["dU_dx"] ** 2 + g["dU_dy"] ** 2 + g["dU_dz"] ** 2 +
                      g["dV_dx"] ** 2 + g["dV_dy"] ** 2 + g["dV_dz"] ** 2 +
                      g["dW_dx"] ** 2 + g["dW_dy"] ** 2 + g["dW_dz"] ** 2)
            acc_grad_sq += _xz_mean(grad_sq)

            if sgs_model != 0 and "nu_t" in snap:
                has_sgs = True
                nu_t = snap["nu_t"]
                S_12 = 0.5 * (g["dU_dy"] + g["dV_dx"])
                S_13 = 0.5 * (g["dU_dz"] + g["dW_dx"])
                S_23 = 0.5 * (g["dV_dz"] + g["dW_dy"])
                two_SijSij = (2.0 * (g["dU_dx"] ** 2 + g["dV_dy"] ** 2 + g["dW_dz"] ** 2) +
                             4.0 * (S_12 ** 2 + S_13 ** 2 + S_23 ** 2))
                acc_eps_sgs += _xz_mean(nu_t * two_SijSij)

            if boussinesq_flag >= 1 and "T" in snap:
                has_T = True
                T = snap["T"]
                acc_T += _xz_mean(T); acc_TT += _xz_mean(T * T)

            N += 1

        U_mean, V_mean, W_mean = acc_U / N, acc_V / N, acc_W / N
        UU_mean, VV_mean, WW_mean, UV_mean = acc_UU / N, acc_VV / N, acc_WW / N, acc_UV / N

        uu = np.maximum(UU_mean - U_mean ** 2, 0.0)
        vv = np.maximum(VV_mean - V_mean ** 2, 0.0)
        ww = np.maximum(WW_mean - W_mean ** 2, 0.0)
        uv = UV_mean - U_mean * V_mean

        dU_dy_mean = np.gradient(U_mean, ym)
        dV_dy_mean = np.gradient(V_mean, ym)
        dW_dy_mean = np.gradient(W_mean, ym)

        eps_res = nu * (acc_grad_sq / N - dU_dy_mean ** 2 - dV_dy_mean ** 2 - dW_dy_mean ** 2)
        eps_sgs = acc_eps_sgs / N if has_sgs else np.zeros(nym)
        T_mean = acc_T / N if has_T else None
        Trms = np.sqrt(np.maximum(acc_TT / N - T_mean ** 2, 0.0)) if has_T else None
        prod = -uv * dU_dy_mean

        U_mean_s, V_mean_s, W_mean_s = _sym_fold(U_mean), _sym_fold(V_mean), _sym_fold(W_mean)
        uu_s, vv_s, ww_s = _sym_fold(uu), _sym_fold(vv), _sym_fold(ww)
        TKE_s = _sym_fold(0.5 * (uu + vv + ww))
        eps_res_s, eps_sgs_s = _sym_fold(eps_res), _sym_fold(eps_sgs)
        uv_s = _antisym_fold(uv)
        urms_s = np.sqrt(np.maximum(uu_s, 0.0))
        vrms_s = np.sqrt(np.maximum(vv_s, 0.0))
        wrms_s = np.sqrt(np.maximum(ww_s, 0.0))
        dU_dy_s = np.gradient(U_mean_s, ym)
        prod_s = -uv_s * dU_dy_s

        u_tau = np.sqrt(nu * np.abs(dU_dy_s[0]))
        Re_tau = u_tau * (Ly / 2.0) / nu

        profiles = dict(U_mean=U_mean_s, V_mean=V_mean_s, W_mean=W_mean_s,
                        urms=urms_s, vrms=vrms_s, wrms=wrms_s, uv=uv_s, TKE=TKE_s,
                        eps_res=eps_res_s, eps_sgs=eps_sgs_s, prod=prod_s)
        # T is deliberately not symmetrised: Boussinesq flows are generally NOT
        # symmetric about the centreline (e.g. bottom-heated/top-cooled), unlike the
        # momentum statistics above for an unstratified channel.
        if has_T:
            profiles["T_mean"] = T_mean
            profiles["Trms"] = Trms

        return FieldStats(ym, profiles, nu=nu, u_tau=u_tau, Re_tau=Re_tau, prefix=prefix,
                          step_range=(start, end), interval=interval, n=N)

    def profile_vs_dns(self, dns_dir=None, x_stations=(), step_start=None, step_end=None,
                       fold=None, u_tau=1.0, h=1.0, sgs_model=None, sediment_flag=None):
        """Spanwise- and snapshot-averaged mean/resolved-stress profiles at one or more
        x-stations (analyse_channel.py's profile_at_station), with an optional
        Moser-Kim-Mansour DNS overlay loaded from `dns_dir`.  Averaging is
        spanwise+temporal only (never streamwise), so profiles at successive x show how
        an inflow develops with fetch.  `fold=None` auto-enables the no-slip/no-slip
        centreline fold when both walls are no-slip (bc_face_ylo=bc_face_yhi=1) and the
        wall-normal grid is itself mirror-symmetric; pass True/False to override.

        Returns a dict: stations=[(x_sel, y, mean[:,3], rey[:,6], n), ...], dns (or
        None), fold_walls, prefix, steps.
        """
        ip = self.case_dir / "input_parameters"
        prm = _core.parse_input_parameters(ip) if ip.exists() else {}
        sgs_model = sgs_model if sgs_model is not None else int(prm.get("sgs_model", 0))
        sediment_flag = sediment_flag if sediment_flag is not None else int(prm.get("sediment_flag", 0))

        all_steps = self.steps()
        if not all_steps:
            raise FileNotFoundError(f"no snapshots in {self.fields_dir}")
        steps = [s for s in all_steps
                if (step_start is None or s >= step_start) and (step_end is None or s <= step_end)]
        if not steps:
            raise ValueError(f"no snapshots in steps=[{step_start},{step_end}]; "
                             f"available {all_steps[0]}..{all_steps[-1]}")

        prefix = self.prefix or _detect_prefix(self.fields_dir)
        s0 = FieldSnapshot.read(self.fields_dir / f"{prefix}.{steps[0]}", sgs_model=sgs_model,
                                sediment_flag=sediment_flag, boussinesq_flag=0)
        xm, ym = s0.xm, s0.ym

        bc_face_ylo = int(prm.get("bc_face_ylo", 1))
        bc_face_yhi = int(prm.get("bc_face_yhi", 1))
        fold_walls = (bc_face_ylo == 1 and bc_face_yhi == 1) if fold is None else fold
        if fold_walls:
            Ly = float(s0.y[-1])
            grid_asym = np.max(np.abs(ym + ym[::-1] - Ly))
            if grid_asym > 1e-6 * Ly:
                warnings.warn("wall-normal grid is not mirror-symmetric about the centreline "
                              f"(max deviation {grid_asym:.3g}); skipping the no-slip/no-slip fold")
                fold_walls = False

        y_out = ym[:len(ym) // 2] if fold_walls else ym
        stations = []
        for xq in x_stations:
            x_sel, mean, rey, n = _profile_at_station(self.fields_dir, prefix, steps, xq, xm,
                                                       sgs_model, sediment_flag)
            if fold_walls:
                mean, rey = _fold_profile(mean, rey)
                half = len(ym) // 2
                mean, rey = mean[:half], rey[:half]
            stations.append((x_sel, y_out, mean, rey, n))

        dns = load_dns_mkm(dns_dir, u_tau, h) if dns_dir else None
        return dict(stations=stations, dns=dns, fold_walls=fold_walls, prefix=prefix, steps=steps)

    # ════════════════════════════════════════════════════════════════════════
    # XDMF export
    # ════════════════════════════════════════════════════════════════════════
    def write_xmf(self, out=None):
        """Write a ParaView XDMF time series for this case's field snapshots.  See
        write_field_xmf for details; uses this instance's case_dir/prefix."""
        return write_field_xmf(self.case_dir, self.prefix, out)


class FieldStats:
    """Time- and (x,z)-plane averaged channel statistics (compute_stats.py's output):
    U_mean, V_mean, W_mean, urms, vrms, wrms, uv, TKE, prod, eps_res, eps_sgs (all (nym,),
    symmetry-folded about the centreline), plus T_mean/Trms when boussinesq_flag >= 1
    (not folded -- see FieldSeries.time_average)."""

    def __init__(self, ym, profiles, nu, u_tau, Re_tau, prefix, step_range, interval, n):
        self.ym = ym
        self.nu, self.u_tau, self.Re_tau = nu, u_tau, Re_tau
        self.prefix, self.step_range, self.interval, self.n = prefix, step_range, interval, n
        self._profiles = dict(profiles)
        for k, v in self._profiles.items():
            setattr(self, k, v)

    def __repr__(self):
        return (f"FieldStats(prefix={self.prefix!r}, steps={self.step_range}, n={self.n}, "
                f"u_tau={self.u_tau:.5f}, Re_tau={self.Re_tau:.1f})")

    def save(self, out_dir="stats"):
        """Write stats_<prefix>_<start>_<end>.npz (load with np.load) and individual
        little-endian float64 .bin profiles under out_dir (compute_stats.py's outputs)."""
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        results = dict(ym=self.ym, **self._profiles)
        start, end = self.step_range
        npz_path = out_dir / f"stats_{self.prefix}_{start}_{end}.npz"
        np.savez_compressed(npz_path, **results)
        for name, arr in results.items():
            (out_dir / f"{name}.bin").write_bytes(np.asarray(arr).astype("<f8").tobytes())
        print(f"Wrote {npz_path}  (+ {len(results)} .bin profiles)")
        return npz_path


# ── small numeric helpers used by time_average ────────────────────────────────

def _xz_mean(arr):
    """Average a (nxm, nym, nzm) array over axes 0 and 2 -> (nym,)."""
    return arr.mean(axis=(0, 2))


def _grad_all(U, V, W, ym, dx, dz):
    """All nine velocity-gradient components at cell centres (central differences in
    the interior, one-sided at boundaries; x/z uniform, y from the stretched ym)."""
    return {
        "dU_dx": np.gradient(U, dx, axis=0), "dU_dy": np.gradient(U, ym, axis=1),
        "dU_dz": np.gradient(U, dz, axis=2),
        "dV_dx": np.gradient(V, dx, axis=0), "dV_dy": np.gradient(V, ym, axis=1),
        "dV_dz": np.gradient(V, dz, axis=2),
        "dW_dx": np.gradient(W, dx, axis=0), "dW_dy": np.gradient(W, ym, axis=1),
        "dW_dz": np.gradient(W, dz, axis=2),
    }


def _sym_fold(a):
    """Symmetric (even) fold about the channel centreline."""
    return 0.5 * (a + a[::-1])


def _antisym_fold(a):
    """Anti-symmetric (odd) fold about the channel centreline (e.g. <u'v'>): carries
    the sign of the lower half (y < Ly/2) after folding."""
    return 0.5 * (a - a[::-1])


# ── profile_vs_dns helpers ─────────────────────────────────────────────────────

# Sign each mean/Reynolds-stress component picks up under the wall-to-wall mirror
# y -> Ly - y (V flips sign, U and W don't; a stress component's sign is the product
# of its two velocity signs).
_MEAN_MIRROR_SIGN = np.array([1.0, -1.0, 1.0])                     # U, V, W
_REY_MIRROR_SIGN = np.array([1.0, 1.0, 1.0, -1.0, 1.0, -1.0])       # uu,vv,ww,uv,uw,vw

# Contrasting marker colours for tracked x-locations in animate_velocity_slice,
# chosen to stay legible against the magma_r colorplot background (ColorBrewer
# Set1, first four).
_LINE_COLORS = ["#e41a1c", "#377eb8", "#4daf4a", "#ff7f00"]


def _profile_at_station(fields_dir, prefix, steps, x_target, xm, sgs_model, sediment_flag):
    """Spanwise- and snapshot-averaged mean and resolved stresses at the x-plane
    nearest x_target.  Returns (x_sel, mean[:,3] (U,V,W), rey[:,6]
    (uu,vv,ww,uv,uw,vw), n)."""
    ix = int(np.argmin(np.abs(xm - x_target)))
    x_sel = float(xm[ix])

    s1 = s2 = None
    n = 0
    for step in steps:
        fpath = Path(fields_dir) / f"{prefix}.{step}"
        snap = FieldSnapshot.read(fpath, sgs_model=sgs_model, sediment_flag=sediment_flag,
                                  boussinesq_flag=0)
        U, V, W = snap["U"][ix], snap["V"][ix], snap["W"][ix]       # each (nym, nzm)

        if s1 is None:
            nym = U.shape[0]
            s1 = np.zeros((nym, 3)); s2 = np.zeros((nym, 6))

        s1[:, 0] += U.sum(axis=1); s1[:, 1] += V.sum(axis=1); s1[:, 2] += W.sum(axis=1)
        s2[:, 0] += (U * U).sum(axis=1); s2[:, 1] += (V * V).sum(axis=1)
        s2[:, 2] += (W * W).sum(axis=1); s2[:, 3] += (U * V).sum(axis=1)
        s2[:, 4] += (U * W).sum(axis=1); s2[:, 5] += (V * W).sum(axis=1)
        n += U.shape[1]

    if n == 0:
        raise ValueError(f"no snapshots found near x = {x_target}")

    mean = s1 / n
    prod_mean = s2 / n
    rey = prod_mean - np.stack(
        [mean[:, 0] * mean[:, 0], mean[:, 1] * mean[:, 1], mean[:, 2] * mean[:, 2],
         mean[:, 0] * mean[:, 1], mean[:, 0] * mean[:, 2], mean[:, 1] * mean[:, 2]], axis=1)
    return x_sel, mean, rey, n


def _fold_profile(mean, rey):
    """Average a profile with its mirror image about the centreline y = Ly/2.  Only
    valid for a channel with matching (no-slip/no-slip) walls and a wall-normal grid
    that is itself symmetric about the centreline -- both checked by the caller."""
    mean_f = 0.5 * (mean + _MEAN_MIRROR_SIGN * mean[::-1])
    rey_f = 0.5 * (rey + _REY_MIRROR_SIGN * rey[::-1])
    return mean_f, rey_f


def read_dat(path, ncol):
    """Whitespace-delimited numeric table, skipping blank/'#'-comment lines."""
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


def load_dns_mkm(dns_dir, u_tau, h):
    """Moser-Kim-Mansour DNS reference (MKM_MEANS.dat, MKM_REYSTRESS.dat) rescaled into
    solver units by u_tau, h; returns None when the files are absent."""
    d = Path(dns_dir)
    means, rey = d / "MKM_MEANS.dat", d / "MKM_REYSTRESS.dat"
    if not (means.exists() and rey.exists()):
        return None
    m, r = read_dat(means, 7), read_dat(rey, 8)
    if len(m) != len(r):
        warnings.warn(f"DNS row counts differ ({len(m)} vs {len(r)}), skipping reference")
        return None
    return {
        "y": m[:, 0] * h, "U": m[:, 2] * u_tau,
        "uu": r[:, 2] * u_tau ** 2, "vv": r[:, 3] * u_tau ** 2,
        "ww": r[:, 4] * u_tau ** 2, "uv": r[:, 5] * u_tau ** 2,
    }


def load_wind_tunnel_inflow(path):
    """Measured wind-tunnel inflow profile CSV: y, U, Iu, Iv, Iw (turbulence
    intensities); returns None if the file has no readable rows."""
    rows = read_dat(path, 5)
    if rows.size == 0:
        return None
    y, U, Iu, Iv, Iw = rows.T
    return {"y": y, "U": U, "Iu": Iu, "Iv": Iv, "Iw": Iw}


def select_equidistant_x(x_list, max_n=4):
    """Pick up to max_n locations, equally spaced between min(x_list) and max(x_list)
    (used to choose which stations get a profile subplot in animate_velocity_slice)."""
    xs = sorted(x_list)
    if len(xs) <= max_n:
        return xs
    return list(np.linspace(xs[0], xs[-1], max_n))


def load_sdf(sdf_path, nxm, nym, nzm):
    """Read the cell-centre SDF written by fdm-dopamine (ibm_input_mode=1).

    Binary layout: big-endian float64, Fortran column-major order, shape
    (nxm+2, nym+2, nzm) with ghost layers in x and y only.  Returns phi with
    shape (nxm, nym, nzm); phi < 0 marks solid cells."""
    nxg, nyg = nxm + 2, nym + 2
    raw = np.fromfile(sdf_path, dtype=">f8")
    expected = nxg * nyg * nzm
    if raw.size != expected:
        raise ValueError(f"{sdf_path}: size mismatch: got {raw.size} elements, "
                         f"expected {expected} ({nxg}x{nyg}x{nzm})")
    phi_full = raw.reshape((nxg, nyg, nzm), order="F")
    return phi_full[1:-1, 1:-1, :]


# ════════════════════════════════════════════════════════════════════════════
# 3. FIGURES / MOVIES
# ════════════════════════════════════════════════════════════════════════════
def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def plot_profile(snap, nu=None, dPdx=None, Ly=None, label=None, save=None):
    """Wall-normal profile figure for one FieldSnapshot: mean streamwise velocity
    (in wall units when a friction velocity can be estimated from dPdx/Ly), plus
    nu_t/C/T panels for whichever of those fields are present (plot_snapshot.py logic).
    nu/dPdx/Ly/label default to `<case>/input_parameters` (case = snap.path's
    grandparent) when snap was produced by FieldSnapshot.read/FieldSeries."""
    plt = _plt()
    case_dir = snap.path.parent.parent if snap.path is not None else None
    ip = case_dir / "input_parameters" if case_dir is not None else None
    prm = _core.parse_input_parameters(ip) if ip is not None and ip.exists() else {}
    nu = nu if nu is not None else float(prm.get("nu", 1e-6))
    dPdx = dPdx if dPdx is not None else float(prm.get("dpdx", 0.0))
    Ly = Ly if Ly is not None else float(prm.get("ly", float(snap.y[-1])))
    label = label if label is not None else str(prm.get("fileout", "case"))

    utau = np.sqrt(abs(dPdx) * 0.5 * Ly)
    use_wall_units = utau > 0.0

    panels = ["U"] + [v for v in ("nu_t", "C", "T") if v in snap]
    ncols = len(panels)
    fig, axes = plt.subplots(1, ncols, figsize=(4.5 * ncols, 5.5), sharey=True)
    axes = np.atleast_1d(axes)

    ym, U_mean = snap.profile("U")
    step_label = f"step {snap.step}" if snap.step is not None else "snapshot"
    fig.suptitle(f"{label}  --  {step_label}", fontsize=10, y=0.98)

    yplot, ylabel = (ym * utau / nu, r"$y^+$") if use_wall_units else (ym, r"$y$  [m]")
    axes[0].set_ylabel(ylabel)

    ax = axes[panels.index("U")]
    if use_wall_units:
        ax.plot(U_mean / utau, yplot, color="C0")
        ax.set_xlabel(r"$\langle U \rangle^+$")
    else:
        ax.plot(U_mean, yplot, color="C0")
        ax.set_xlabel(r"$\langle U \rangle$  [m s$^{-1}$]")
    ax.set_title("Mean streamwise velocity")
    ax.grid(True, alpha=0.3)

    if "nu_t" in panels:
        _, nut_mean = snap.profile("nu_t")
        ax = axes[panels.index("nu_t")]
        ax.plot(nut_mean / nu, yplot, color="C1")
        ax.set_xlabel(r"$\langle \nu_t \rangle / \nu$")
        ax.set_title("SGS eddy-viscosity ratio")
        ax.grid(True, alpha=0.3)

    if "C" in panels:
        _, C_mean = snap.profile("C")
        ax = axes[panels.index("C")]
        ax.plot(C_mean, yplot, color="C2")
        ax.set_xlabel(r"$\langle C \rangle$")
        ax.set_title("Mean scalar concentration")
        ax.grid(True, alpha=0.3)

    if "T" in panels:
        _, T_mean = snap.profile("T")
        ax = axes[panels.index("T")]
        ax.plot(T_mean, yplot, color="C3")
        ax.set_xlabel(r"$\langle T \rangle$")
        ax.set_title("Mean temperature")
        ax.grid(True, alpha=0.3)

    fig.tight_layout()
    if save:
        fig.savefig(save, dpi=150, bbox_inches="tight")
        print(f"Wrote {save}")
    return fig


_LES_STYLE = dict(color="k", marker="o", linestyle="none", markerfacecolor="none",
                  markeredgecolor="k", markersize=4.5, markeredgewidth=0.9)


def _half_channel(stats, u_tau=1.0):
    """Wall-unit profiles for the lower half-channel, from a FieldStats.

    u_tau defaults to 1.0, matching compute_stats.py's plot_profiles (which
    normalizes by a fixed u_tau=1.0 for the *plot*, not by the u_tau derived from
    the mean wall shear during time_average -- that one is reported separately as
    stats.u_tau/stats.Re_tau). Pass u_tau=stats.u_tau explicitly to normalize the
    figure by the measured friction velocity instead.
    """
    nu = stats.nu
    delta_nu = nu / u_tau
    nym = stats.ym.size
    half = nym // 2

    yp = stats.ym[:half] / delta_nu
    fac = nu / u_tau ** 4
    epsr_p = stats.eps_res[:half] * fac
    epss_p = stats.eps_sgs[:half] * fac
    return dict(
        yp=yp, Up=stats.U_mean[:half] / u_tau,
        urp=stats.urms[:half] / u_tau, vrp=stats.vrms[:half] / u_tau, wrp=stats.wrms[:half] / u_tau,
        uvp=-stats.uv[:half] / u_tau ** 2, TKEp=stats.TKE[:half] / u_tau ** 2,
        prod_p=stats.prod[:half] * fac, epsr_p=epsr_p, epss_p=epss_p, eps_tot=epsr_p + epss_p,
    )


def _make_stats_figure(plt):
    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(2, 3, hspace=0.38, wspace=0.32)
    return (fig, fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[0, 2]),
           fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1]), fig.add_subplot(gs[1, 2]))


def _plot_dns_stats(ax_U, ax_rms, ax_uv, ax_TKE, ax_bud, ax_eps, dns):
    ax_U.semilogx(dns["yp"], dns["Up"], color="tab:blue", lw=1.6, label="DNS")
    ax_rms.plot(dns["urp"], dns["yp"], color="tab:red", lw=1.6, label=r"DNS $u'_{rms}$")
    ax_rms.plot(dns["vrp"], dns["yp"], color="tab:green", lw=1.6, label=r"DNS $v'_{rms}$")
    ax_rms.plot(dns["wrp"], dns["yp"], color="tab:purple", lw=1.6, label=r"DNS $w'_{rms}$")
    ax_uv.plot(dns["uvp"], dns["yp"], color="tab:blue", lw=1.6, label="DNS")
    ax_TKE.plot(dns["TKEp"], dns["yp"], color="tab:orange", lw=1.6, label="DNS")
    ax_bud.plot(dns["prod_p"], dns["yp"], color="tab:red", lw=1.6, label=r"DNS $P^+$")
    ax_bud.plot(dns["eps_budget_p"], dns["yp"], color="tab:blue", lw=1.6, label=r"DNS $-\varepsilon^+$")
    ax_eps.plot(dns["eps_p"], dns["yp"], color="tab:blue", lw=1.6, label=r"DNS $\varepsilon^+$")


def load_mkm_stats_reference(case_dir):
    """MKM_MEANS.dat/MKM_REYSTRESS.dat reference used by plot_stats, in *wall units
    already* (unlike load_dns_mkm's MKM_*.dat loader, which rescales by u_tau/h for
    profile_vs_dns's physical-unit overlay -- these are two distinct column
    conventions from the two original scripts, kept separate on purpose).  Returns
    None when the files are absent."""
    case_dir = Path(case_dir)
    means_path = case_dir / "MKM_MEANS.dat"
    reystress_path = case_dir / "MKM_REYSTRESS.dat"
    if not means_path.exists() or not reystress_path.exists():
        return None
    means = np.loadtxt(means_path, comments="#")
    reystress = np.loadtxt(reystress_path, comments="#")
    return {
        "yp": means[:, 1], "Up": means[:, 2],
        "urp": np.sqrt(np.clip(reystress[:, 2], 0.0, None)),
        "vrp": np.sqrt(np.clip(reystress[:, 3], 0.0, None)),
        "wrp": np.sqrt(np.clip(reystress[:, 4], 0.0, None)),
        "uvp": -reystress[:, 5],
        "TKEp": 0.5 * (reystress[:, 2] + reystress[:, 3] + reystress[:, 4]),
    }


def _plot_mkm_stats(ax_U, ax_rms, ax_uv, ax_TKE, mkm):
    ax_U.semilogx(mkm["yp"], mkm["Up"], color="tab:brown", lw=1.3, ls="--", label="MKM DNS")
    ax_rms.plot(mkm["urp"], mkm["yp"], color="tab:brown", lw=1.3, ls="--", label=r"MKM $u'_{rms}$")
    ax_rms.plot(mkm["vrp"], mkm["yp"], color="tab:brown", lw=1.3, ls="--", label=r"MKM $v'_{rms}$")
    ax_rms.plot(mkm["wrp"], mkm["yp"], color="tab:brown", lw=1.3, ls="--", label=r"MKM $w'_{rms}$")
    ax_uv.plot(mkm["uvp"], mkm["yp"], color="tab:brown", lw=1.3, ls="--", label="MKM DNS")
    ax_TKE.plot(mkm["TKEp"], mkm["yp"], color="tab:brown", lw=1.3, ls="--", label="MKM DNS")


def load_channel_dns_budget(dns_dir):
    """DNS mean/fluctuation/budget reference used by plot_stats (mean_prof.txt,
    vel_fluc_prof.txt, RSTE_k_prof.txt -- a different reference set than
    load_dns_mkm's MKM_*.dat, used by profile_vs_dns)."""
    dns_dir = Path(dns_dir)
    files = ["mean_prof.txt", "vel_fluc_prof.txt", "RSTE_k_prof.txt"]
    if not all((dns_dir / name).exists() for name in files):
        return None
    mean = np.loadtxt(dns_dir / "mean_prof.txt", comments="%")
    vel = np.loadtxt(dns_dir / "vel_fluc_prof.txt", comments="%")
    rste_k = np.loadtxt(dns_dir / "RSTE_k_prof.txt", comments="%")
    uv_budget = rste_k[:, 7]
    if np.nanmean(uv_budget) > 0.0:
        uv_budget = -uv_budget
    return {
        "yp": mean[:, 1], "Up": mean[:, 2],
        "urp": np.sqrt(np.clip(vel[:, 2], 0.0, None)), "vrp": np.sqrt(np.clip(vel[:, 3], 0.0, None)),
        "wrp": np.sqrt(np.clip(vel[:, 4], 0.0, None)), "uvp": -vel[:, 5], "TKEp": vel[:, 8],
        "prod_p": rste_k[:, 2], "eps_budget_p": uv_budget, "eps_p": -uv_budget,
    }


def plot_stats(stats, case_dir=".", dns_dir=None, u_tau=1.0, save="plots/stats.png", dpi=150):
    """Wall-unit profile figure for a FieldStats: mean U, velocity RMS, Reynolds shear
    stress, TKE, production/dissipation budget, and the resolved/SGS dissipation
    breakdown, each vs y+; overlaid with a DNS reference (mean_prof.txt/
    vel_fluc_prof.txt/RSTE_k_prof.txt under `dns_dir`, default `<case_dir>/dns_ref`)
    and/or an MKM reference (MKM_MEANS.dat/MKM_REYSTRESS.dat directly under
    `case_dir`) when present -- compute_stats.py's plot_profiles.  See
    _half_channel's docstring for why u_tau defaults to 1.0 here rather than
    stats.u_tau."""
    plt = _plt()
    case_dir = Path(case_dir)
    save = Path(save)
    dns_dir = Path(dns_dir) if dns_dir is not None else case_dir / "dns_ref"
    dns = load_channel_dns_budget(dns_dir) if dns_dir.is_dir() else None
    mkm = load_mkm_stats_reference(case_dir)

    fig, ax_U, ax_rms, ax_uv, ax_TKE, ax_bud, ax_eps = _make_stats_figure(plt)
    if dns is not None:
        _plot_dns_stats(ax_U, ax_rms, ax_uv, ax_TKE, ax_bud, ax_eps, dns)
    if mkm is not None:
        _plot_mkm_stats(ax_U, ax_rms, ax_uv, ax_TKE, mkm)

    Re_tau = u_tau * (stats.ym.max() / 2.0) / stats.nu
    print(f"{stats.prefix}:  u_tau = {u_tau:.5f},  Re_tau = {Re_tau:.1f}")

    h = _half_channel(stats, u_tau=u_tau)
    label = stats.prefix
    ax_U.semilogx(h["yp"], h["Up"], label=label, **_LES_STYLE)
    ax_rms.plot(h["urp"], h["yp"], label=rf"{label} $u'_{{rms}}$", **_LES_STYLE)
    ax_rms.plot(h["vrp"], h["yp"], label=rf"{label} $v'_{{rms}}$", **_LES_STYLE)
    ax_rms.plot(h["wrp"], h["yp"], label=rf"{label} $w'_{{rms}}$", **_LES_STYLE)
    ax_uv.plot(h["uvp"], h["yp"], label=label, **_LES_STYLE)
    ax_TKE.plot(h["TKEp"], h["yp"], label=label, **_LES_STYLE)
    ax_bud.plot(h["prod_p"], h["yp"], label=rf"{label} $P^+$", **_LES_STYLE)
    ax_bud.plot(-h["eps_tot"], h["yp"], label=rf"{label} $-\varepsilon^+_{{tot}}$", **_LES_STYLE)
    ax_eps.plot(h["epsr_p"], h["yp"], label=rf"{label} $\varepsilon^+_{{res}}$", **_LES_STYLE)
    ax_eps.plot(h["epss_p"], h["yp"], label=rf"{label} $\varepsilon^+_{{SGS}}$", **_LES_STYLE)

    ax_U.set_xlabel(r"$y^+$"); ax_U.set_ylabel(r"$U^+$")
    ax_U.set_title("Mean streamwise velocity"); ax_U.legend(fontsize=7)
    ax_rms.set_xlabel("Fluctuation amplitude in wall units"); ax_rms.set_ylabel(r"$y^+$")
    ax_rms.set_title("Velocity RMS")
    ax_uv.set_xlabel(r"$-\langle u'v' \rangle^+$"); ax_uv.set_ylabel(r"$y^+$")
    ax_uv.set_title("Reynolds shear stress")
    ax_TKE.set_xlabel(r"$k^+$"); ax_TKE.set_ylabel(r"$y^+$"); ax_TKE.set_title("Turbulent kinetic energy")
    ax_bud.set_xlabel(r"$P^+,\;-\varepsilon^+$"); ax_bud.set_ylabel(r"$y^+$")
    ax_bud.axvline(0, color="k", lw=0.6, ls=":"); ax_bud.set_title("TKE production & dissipation")
    ax_eps.set_xlabel(r"$\varepsilon^+$"); ax_eps.set_ylabel(r"$y^+$")
    ax_eps.set_title("Dissipation breakdown (res. / SGS)")

    ax_U.grid(True, which="both", ls=":", lw=0.4, alpha=0.6); ax_U.set_xlim(left=1)
    for ax in (ax_rms, ax_uv, ax_TKE, ax_bud, ax_eps):
        ax.grid(True, ls=":", lw=0.4, alpha=0.6)
        ax.set_ylim(bottom=1)

    fig.suptitle("Channel flow statistics -- wall units", fontsize=12, y=1.01)
    save.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(save, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {save}")
    return fig


def plot_profile_vs_dns(result, inflow=None, title=None, save="channel_profiles.png"):
    """U, Iu, Iv, Iw panels for every station in `result` (from
    FieldSeries.profile_vs_dns), optionally overlaid with a DNS reference
    (result['dns']) and/or a measured wind-tunnel inflow profile (analyse_channel.py's
    make_figure)."""
    plt = _plt()
    stations, dns = result["stations"], result["dns"]
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
            a.plot(inflow[key], inflow["y"], "ko", lw=1.6, markerfacecolor="none",
                  markeredgecolor="black", markevery=2, label="wind tunnel", zorder=5)

    for a, xl in zip(ax, [r"$U$", r"$I_u$", r"$I_v$", r"$I_w$"]):
        a.set_xlabel(xl); a.set_ylabel(r"$y$"); a.grid(alpha=0.25); a.legend(fontsize=8, frameon=False)

    if title is None:
        title = f"{result['prefix']}: spanwise- and snapshot-averaged profiles"
    fig.suptitle(title, y=0.99)
    fig.tight_layout()
    fig.savefig(save, dpi=150)
    print(f"Wrote {save}")
    return fig


def animate_velocity_slice(series, steps, z_target=None, out="velocity.gif", fps=1.0,
                           phi=None, line_x=(), inflow=None, sgs_model=None, sediment_flag=None,
                           workers=8, vmin=0.0, vmax=20.0):
    """Animate the (x, y) velocity-magnitude slice at the z-plane nearest z_target, over
    `steps` (analyse_channel.py's animate_velocity_gif).

    When `phi` (the cell-centre SDF, shape (nxm,nym,nzm), see load_sdf) is given, solid
    cells at that z-plane are masked grey.  `line_x`, when given, marks each
    x-location with a dashed vertical line and adds subplots below showing the
    wall-normal profiles of spanwise-averaged streamwise velocity and turbulence
    intensities (Iu, Iv, Iw) at that x-plane, redrawn every frame; `inflow` (from
    load_wind_tunnel_inflow) overlays a static reference profile on each subplot.

    Uses _core.render_movie's ffmpeg two-pass palette encoder (parallel per-frame
    Figure rendering) in place of the original matplotlib FuncAnimation/PillowWriter,
    for consistency with the rest of the package; this changes the *rendering*
    mechanism but not the visual content (same frames, same fixed [vmin, vmax] colour
    scale so successive frames stay comparable).
    """
    import matplotlib

    prm = {}
    ip = series.case_dir / "input_parameters"
    if ip.exists():
        prm = _core.parse_input_parameters(ip)
    sgs_model = sgs_model if sgs_model is not None else int(prm.get("sgs_model", 0))
    sediment_flag = sediment_flag if sediment_flag is not None else int(prm.get("sediment_flag", 0))
    prefix = series.prefix or _detect_prefix(series.fields_dir)
    read_kw = dict(sgs_model=sgs_model, sediment_flag=sediment_flag, boussinesq_flag=0)

    s0 = FieldSnapshot.read(series.fields_dir / f"{prefix}.{steps[0]}", **read_kw)
    xm, ym, zm = s0.xm, s0.ym, s0.zm
    iz = int(np.argmin(np.abs(zm - (z_target if z_target is not None else zm[len(zm) // 2]))))
    z_sel = float(zm[iz])
    solid = phi[:, :, iz] < 0 if phi is not None else None

    ix_lines = [int(np.argmin(np.abs(xm - xq))) for xq in line_x]
    x_line_sel = [float(xm[ix]) for ix in ix_lines]
    colors = [_LINE_COLORS[i % len(_LINE_COLORS)] for i in range(len(ix_lines))]

    frames = []
    profiles_U = [[] for _ in ix_lines]
    profiles_Iu = [[] for _ in ix_lines]
    profiles_Iv = [[] for _ in ix_lines]
    profiles_Iw = [[] for _ in ix_lines]
    for step in steps:
        snap = FieldSnapshot.read(series.fields_dir / f"{prefix}.{step}", **read_kw)
        U, V, W = snap["U"][:, :, iz], snap["V"][:, :, iz], snap["W"][:, :, iz]
        speed = np.sqrt(U * U + V * V + W * W)
        if solid is not None:
            speed = np.ma.masked_where(solid, speed)
        frames.append(speed)

        for k, ix in enumerate(ix_lines):
            Ul, Vl, Wl = snap["U"][ix], snap["V"][ix], snap["W"][ix]     # (nym, nzm)
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

    cmap = matplotlib.colormaps["magma_r"].copy()
    cmap.set_bad(color="0.5")

    fig_w = 12.0
    img_w_frac = 0.90
    data_aspect = (ym[-1] - ym[0]) / (xm[-1] - xm[0])
    img_h = fig_w * img_w_frac * data_aspect

    prof_titles = [r"$U$", r"$I_u$", r"$I_v$", r"$I_w$"]
    prof_data = [profiles_U, profiles_Iu, profiles_Iv, profiles_Iw]
    inflow_keys = [None, "Iu", "Iv", "Iw"]

    def frame(i):
        plt = _plt()
        if ix_lines:
            prof_h = 2.4
            fig = plt.figure(figsize=(fig_w, img_h + prof_h + 1.1))
            gs = fig.add_gridspec(2, 4, height_ratios=[img_h, prof_h], hspace=0.08, wspace=0.35,
                                  left=0.06, right=0.98, top=0.95, bottom=0.10)
            ax = fig.add_subplot(gs[0, :])
            ax_prof = [fig.add_subplot(gs[1, j]) for j in range(4)]
        else:
            fig, ax = plt.subplots(figsize=(fig_w, img_h + 0.9))
            ax_prof = []

        im = ax.imshow(frames[i].T, origin="lower", aspect="equal",
                       extent=[xm[0], xm[-1], ym[0], ym[-1]], vmin=vmin, vmax=vmax, cmap=cmap)
        fig.colorbar(im, ax=ax, label=r"$|U|$", shrink=0.5)
        ax.set_xlabel(r"$x$"); ax.set_ylabel(r"$y$")
        ax.set_title(f"z = {z_sel:.3g}, step {steps[i]}")
        for x_sel, c in zip(x_line_sel, colors):
            ax.axvline(x_sel, color=c, ls="--", lw=1.6)

        for j, (a, ttl, data) in enumerate(zip(ax_prof, prof_titles, prof_data)):
            for k, c in enumerate(colors):
                a.plot(data[k][i], ym, color=c, lw=1.3, label=f"x = {x_line_sel[k]:.2f}")
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
        return fig

    render_movie(frame, len(frames), out, fps=fps, workers=workers)
    print(f"Wrote {out}  ({len(frames)} frames, z = {z_sel:.4g})")


# ════════════════════════════════════════════════════════════════════════════
# XDMF export (zero-copy byte-seek into the raw snapshot files)
# ════════════════════════════════════════════════════════════════════════════
def _probe_field_header(fpath, has_sediment, has_sgs, has_boussinesq):
    """Header-only probe: grid + field shapes/byte offsets for the XDMF writer, without
    materializing the (possibly large) field blocks into arrays.

    Implemented directly with struct rather than _core.BinaryReader.rd (which always
    copies into an ndarray): this must SKIP multi-megabyte field blocks cheaply, which
    is the one place in this module where that distinction matters.
    """
    data = Path(fpath).read_bytes()
    pos = [0]

    def ri():
        v = struct.unpack_from(">i", data, pos[0])[0]
        pos[0] += 4
        return v

    def ri3():
        v = struct.unpack_from(">3i", data, pos[0])
        pos[0] += 12
        return v

    def skip8(n):
        pos[0] += n * 8

    def rd(n):
        arr = np.frombuffer(data, dtype=">f8", count=n, offset=pos[0]).copy()
        pos[0] += n * 8
        return arr

    nx = ri(); skip8(nx)
    ny = ri(); skip8(ny)
    nz = ri(); skip8(nz)
    nxm = ri(); xm = rd(nxm)
    nym = ri(); ym = rd(nym)
    nzm = ri(); zm = rd(nzm)

    un = ri3(); off_U = pos[0]; skip8(un[0] * un[1] * un[2])
    vn = ri3(); off_V = pos[0]; skip8(vn[0] * vn[1] * vn[2])
    wn = ri3(); off_W = pos[0]; skip8(wn[0] * wn[1] * wn[2])
    pn = ri3(); off_P = pos[0]; skip8(pn[0] * pn[1] * pn[2])

    cn = off_C = nutn = off_nut = Tn = off_T = None
    if has_sediment:
        cn = ri3(); off_C = pos[0]; skip8(cn[0] * cn[1] * cn[2])
    if has_sgs:
        nutn = ri3(); off_nut = pos[0]; skip8(nutn[0] * nutn[1] * nutn[2])
    if has_boussinesq:
        Tn = ri3(); off_T = pos[0]

    return dict(nxm=nxm, nym=nym, nzm=nzm, xm=xm, ym=ym, zm=zm, un=un, vn=vn, wn=wn, pn=pn,
               off_U=off_U, off_V=off_V, off_W=off_W, off_P=off_P,
               cn=cn, off_C=off_C, nutn=nutn, off_nut=off_nut, Tn=Tn, off_T=off_T)


def _xdmf_hyperslab(attr_name, full_dims_str, start_str, seek, fpath_rel, nxm, nym, nzm):
    count = f"{nzm} {nym} {nxm}"
    return [
        f'        <Attribute Name="{attr_name}" Center="Node" AttributeType="Scalar">',
        f'          <DataItem ItemType="HyperSlab" Dimensions="{count}" Type="HyperSlab">',
        '            <DataItem Dimensions="3 3" Format="XML">',
        f'              {start_str}',
        '              1 1 1',
        f'              {count}',
        '            </DataItem>',
        f'            <DataItem Dimensions="{full_dims_str}" Format="Binary"',
        f'                     DataType="Float" Precision="8" Endian="Big" Seek="{seek}">',
        f'              {fpath_rel}',
        '            </DataItem>',
        '          </DataItem>',
        '        </Attribute>',
    ]


def _xdmf_hyperslab_item(full_dims_str, start_str, seek, fpath_rel, nxm, nym, nzm):
    count = f"{nzm} {nym} {nxm}"
    return [
        f'            <DataItem ItemType="HyperSlab" Dimensions="{count}" Type="HyperSlab">',
        '              <DataItem Dimensions="3 3" Format="XML">',
        f'                {start_str}',
        '                1 1 1',
        f'                {count}',
        '              </DataItem>',
        f'              <DataItem Dimensions="{full_dims_str}" Format="Binary"',
        f'                       DataType="Float" Precision="8" Endian="Big" Seek="{seek}">',
        f'                {fpath_rel}',
        '              </DataItem>',
        '            </DataItem>',
    ]


def _xdmf_velocity(un, vn, wn, offU, offV, offW, fpath_rel, nxm, nym, nzm):
    count = f"{nzm} {nym} {nxm}"
    ux, uy, uz = un
    vx, vy, vz = vn
    wx, wy, wz = wn
    lines = [
        '        <Attribute Name="Velocity" Center="Node" AttributeType="Vector">',
        f'          <DataItem ItemType="Function" Dimensions="{count} 3"',
        '                   Function="JOIN($0, $1, $2)">',
    ]
    lines += _xdmf_hyperslab_item(f"{uz} {uy} {ux}", "1 1 0", offU, fpath_rel, nxm, nym, nzm)
    lines += _xdmf_hyperslab_item(f"{vz} {vy} {vx}", "1 0 1", offV, fpath_rel, nxm, nym, nzm)
    lines += _xdmf_hyperslab_item(f"{wz} {wy} {wx}", "0 1 1", offW, fpath_rel, nxm, nym, nzm)
    lines += ['          </DataItem>', '        </Attribute>']
    return lines


def write_field_xmf(case_dir=".", prefix=None, out=None):
    """Write a ParaView XDMF time series for `<case_dir>/fields/<prefix>.<step>`
    binary field snapshots (generateXMF.py).  Points directly at the original binary
    files via byte-seek HyperSlab selections -- no field data is duplicated.  Writes
    `paraview/xm.bin`, `ym.bin`, `zm.bin` (little-endian cell-centre coordinate
    arrays) alongside `<out>` (default `paraview/channel_test.xmf`).

    Fields are on a staggered grid; ghost-layer stripping and the selection of the
    lower face value per cell (no averaging) are handled entirely within the XDMF
    HyperSlab, exactly as in generateXMF.py.  Which optional blocks (C, nu_t, T) are
    present is auto-detected from `<case_dir>/input_parameters` (sediment_flag,
    sgs_model, boussinesq_flag) when available.
    """
    case = Path(case_dir)
    fields_dir = case / "fields"
    out_dir = case / "paraview"
    out_dir.mkdir(exist_ok=True)
    out_path = Path(out) if out else out_dir / "channel_test.xmf"

    snaps = list_indexed(fields_dir, _PAT, prefix)

    ip_path = case / "input_parameters"
    has_sediment = has_sgs = has_boussinesq = False
    if ip_path.exists():
        prm = _core.parse_input_parameters(ip_path)
        has_sediment = int(prm.get("sediment_flag", 0)) >= 1
        has_sgs = int(prm.get("sgs_model", 0)) != 0
        has_boussinesq = int(prm.get("boussinesq_flag", 0)) >= 1

    info = None
    body_lines = []
    nxm = nym = nzm = None
    for step, fpath in snaps:
        h = _probe_field_header(fpath, has_sediment, has_sgs, has_boussinesq)
        if info is None:
            info = h
            nxm, nym, nzm = info["nxm"], info["nym"], info["nzm"]
            (out_dir / "xm.bin").write_bytes(info["xm"].astype("<f8").tobytes())
            (out_dir / "ym.bin").write_bytes(info["ym"].astype("<f8").tobytes())
            (out_dir / "zm.bin").write_bytes(info["zm"].astype("<f8").tobytes())

        frel = os.path.relpath(fpath, out_path.parent)
        body_lines += [
            '',
            f'      <Grid Name="t{float(step):.4g}">',
            f'        <Time Value="{float(step):.6g}"/>',
            '        <Topology Reference="/Xdmf/Domain/Topology[1]"/>',
            '        <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>',
        ]
        body_lines += _xdmf_velocity(h["un"], h["vn"], h["wn"], h["off_U"], h["off_V"], h["off_W"],
                                     frel, nxm, nym, nzm)
        px, py, pz = h["pn"]
        body_lines += _xdmf_hyperslab("P", f"{pz} {py} {px}", "1 1 1", h["off_P"], frel, nxm, nym, nzm)
        if has_sediment and h["off_C"] is not None:
            cx, cy, cz = h["cn"]
            body_lines += _xdmf_hyperslab("C", f"{cz} {cy} {cx}", "1 1 1", h["off_C"], frel, nxm, nym, nzm)
        if has_sgs and h["off_nut"] is not None:
            ntx, nty, ntz = h["nutn"]
            body_lines += _xdmf_hyperslab("nu_t", f"{ntz} {nty} {ntx}", "1 1 1", h["off_nut"], frel, nxm, nym, nzm)
        if has_boussinesq and h["off_T"] is not None:
            Tx, Ty, Tz = h["Tn"]
            body_lines += _xdmf_hyperslab("T", f"{Tz} {Ty} {Tx}", "1 1 1", h["off_T"], frel, nxm, nym, nzm)
        body_lines.append('      </Grid>')

    if info is None:
        raise FileNotFoundError(f"no field snapshots found in {fields_dir}")

    lines = [
        '<?xml version="1.0" ?>',
        '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>',
        '<Xdmf Version="2.0">',
        '  <Domain>',
        '',
        f'    <!-- Rectilinear mesh: {nxm} x {nym} x {nzm} cell-centre nodes -->',
        '    <Topology name="topo" TopologyType="3DRECTMesh"',
        f'              Dimensions="{nzm} {nym} {nxm}"/>',
        '',
        '    <Geometry name="geo" Type="VxVyVz">',
        f'      <DataItem Name="Vx" Dimensions="{nxm}" Format="Binary"',
        '               DataType="Float" Precision="8" Endian="Little">',
        '        xm.bin',
        '      </DataItem>',
        f'      <DataItem Name="Vy" Dimensions="{nym}" Format="Binary"',
        '               DataType="Float" Precision="8" Endian="Little">',
        '        ym.bin',
        '      </DataItem>',
        f'      <DataItem Name="Vz" Dimensions="{nzm}" Format="Binary"',
        '               DataType="Float" Precision="8" Endian="Little">',
        '        zm.bin',
        '      </DataItem>',
        '    </Geometry>',
        '',
        '    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">',
    ] + body_lines + [
        '    </Grid>',
        '',
        '  </Domain>',
        '</Xdmf>',
    ]

    out_path.write_text('\n'.join(lines) + '\n')
    print(f"Wrote {out_path}  ({len(snaps)} snapshots)")
    return out_path
