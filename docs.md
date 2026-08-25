# fdm-dopamine — Documentation

Notice: This file was generated using an LLM and subsequently edited for correctness. If you find a mistake please create an issue on github. Thank you!

## Overview

`fdm-dopamine` is a parallel, finite-difference solver for the three-dimensional incompressible Navier–Stokes equations, targeting turbulent channel and open-channel flows.  It runs on distributed-memory systems via MPI (1-D z-slab decomposition) and supports large-eddy simulation (LES) as well as direct numerical simulation (DNS).

---

## 1. Governing equations

The solver integrates the incompressible Navier–Stokes equations in non-dimensionalised or dimensional form:

$$\frac{\partial u_i}{\partial t} + \frac{\partial (u_i u_j)}{\partial x_j} = -\frac{1}{\rho}\frac{\partial p}{\partial x_i} + \frac{\partial}{\partial x_j}\left[(\nu + \nu_t)\left(\frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i}\right)\right] + f_i$$

$$\frac{\partial u_i}{\partial x_i} = 0$$

where $\nu$ is the kinematic viscosity, $\nu_t$ is the sub-grid scale (SGS) eddy viscosity (zero for DNS), and $f_i$ is an imposed body force (typically a mean pressure gradient).

Streamwise ($x$) and spanwise ($z$) forcing:

$$f_x(t) = \frac{\mathrm{d}P}{\mathrm{d}x}\bigg|_0 + U_{b,x}\,\omega_x\cos(\omega_x t + \varphi_x), \qquad f_z(t) = \frac{\mathrm{d}P}{\mathrm{d}z}\bigg|_0 + U_{b,z}\,\omega_z\cos(\omega_z t + \varphi_z)$$

where $\omega_x = 2\pi / T_{\mathrm{wave},x}$ and $\omega_z = 2\pi / T_{\mathrm{wave},z}$ are derived from the user-supplied wave periods, and $\varphi_x$, $\varphi_z$ are optional phase offsets (radians, default 0).  Setting $T_{\mathrm{wave},x} = 0$ (or $T_{\mathrm{wave},z} = 0$) disables oscillation in that direction, recovering steady forcing.

---

## 2. Spatial discretisation

### 2.1 Staggered MAC grid

Velocities are stored on a **marker-and-cell (MAC)** staggered grid (Harlow & Welch 1965):

| Variable | Location |
|----------|----------|
| $U$ | $x$-face centres — array `(nx, nyg, nzg)` |
| $V$ | $y$-face centres — array `(nxg, ny, nzg)` |
| $W$ | $z$-face centres — array `(nxg, nyg, nz)` |
| $P$, scalars, $\nu_t$ | Cell centres — array `(nxg, nyg, nzg)` |

The grid is **uniform** in the homogeneous $x$ and $z$ directions (spacing $\Delta x$, $\Delta z$) and **non-uniform** in the wall-normal $y$ direction.  Seven vertical-grid options are provided (uniform, symmetric/single-sided hyperbolic tangent stretching, and roughness-sublayer variants), all configured via `grid_type` and `alpha_grid` in `input_parameters`.

### 2.2 Finite differences

All spatial derivatives are discretised with **second-order central differences**.  Cross-derivative terms (e.g. $\partial u/\partial y$ on a $y$-face) require interpolation between the two staggered positions; this is performed by arithmetic averaging (linear interpolation in uniform-grid directions, weight-averaged in $y$).

---

## 3. Time integration

### 3.1 Low-storage explicit Runge–Kutta 3

Time advancement uses the **three-stage, low-storage Runge–Kutta** scheme of Wray (1990) — the same scheme used in the Kim, Moin & Moser (1987) channel flow code and many subsequent LES solvers:

$$u^{(s)} = u^{(s-1)} + \alpha_s \Delta t\, R^{(s-1)}, \quad s = 1,2,3$$

Coefficients: $\alpha_1 = 8/15$, $\alpha_2 = 5/12$, $\alpha_3 = 3/4$.

> **Reference:** Wray, A.A. (1990). *Minimal storage time advancement schemes for spectral methods*. NASA Ames Research Center Report.

### 3.2 CFL condition and adaptive time-stepping

The convective CFL $C_\text{conv} = |u|\Delta t / \Delta x$ and viscous CFL $C_\text{visc} = (\nu+\nu_t)\Delta t / \Delta x^2$ are evaluated at the start of each step.  With `cfl_adaptive = 1` the time-step is adjusted as:

$$\Delta t_\text{new} = \Delta t_\text{old} \cdot \frac{C_\text{target}}{C_\text{current}} \cdot f_\text{safety}$$

### 3.3 Stopping and save-interval criteria

By default the run stops after a fixed number of steps (`nsteps`) and writes a field snapshot every fixed number of steps (`nsave`). Two physical-time-based alternatives are available:

- **`nsteps < 0`** — the step count is ignored; the run instead stops once the physical time reaches `sim_end_time` (checked once per step, after `t` is advanced).
- **`nsave < 0`** — the step count is ignored; a snapshot is instead written every `tsave` physical time units. Because `dt` can vary (`cfl_adaptive = 1`, or a fixed `dt` that does not divide `tsave` evenly), the step immediately before a save is shrunk — `dt = tsave_next - t` — so that `t` lands exactly on multiples of `tsave` and the physical-time spacing between snapshots stays perfectly uniform.

These two switches are independent: `nsteps` and `nsave` can each be positive or negative in any combination.

---

## 4. Pressure projection (fractional step method)

The incompressibility constraint is enforced via a **fractional step / projection method** (Chorin 1968; Kim & Moin 1985):

1. Compute an intermediate velocity $u^*$ from the RK3 momentum update (without enforcing $\nabla\cdot u^* = 0$).
2. Solve the pseudo-pressure Poisson equation:

$$\nabla^2 \phi = \frac{1}{\Delta t}\nabla\cdot u^*$$

3. Project to a divergence-free field:

$$u^{n+1} = u^* - \Delta t\,\nabla\phi$$

### 4.1 Spectral Poisson solver

Taking discrete Fourier transforms in the homogeneous $x$ and $z$ directions (via **FFTW3-MPI**, Frigo & Johnson 2005) converts the Poisson equation into a set of independent 1-D tridiagonal systems in $y$, one per Fourier mode pair $(k_x, k_z)$:

$$\left(-k_x^2 - k_z^2 + \frac{\partial^2}{\partial y^2}\right)\hat{\phi}(k_x, y, k_z) = \hat{f}(k_x, y, k_z)$$

Each tridiagonal system is solved by the LAPACK complex tridiagonal routine `Zgtsv`.  The implementation uses FFTW's **transposed-output** layout so that the FFT-distributed data aligns with the tridiagonal solve without extra communication.

> **Reference (FFTW3):** Frigo, M. & Johnson, S.G. (2005). *The design and implementation of FFTW3*. Proc. IEEE 93(2), 216–231. DOI: [10.1109/JPROC.2004.840301](https://doi.org/10.1109/JPROC.2004.840301)

> **Reference (fractional step):** Kim, J. & Moin, P. (1985). *Application of a fractional-step method to incompressible Navier–Stokes equations*. J. Comput. Phys. 59, 308–323. DOI: [10.1016/0021-9991(85)90148-2](https://doi.org/10.1016/0021-9991(85)90148-2)

---

## 5. Sub-grid scale model

The SGS eddy viscosity $\nu_t$ is provided by the **Vreman (2004) model**, controlled by `sgs_model = 1` in `input_parameters` (setting `sgs_model = 0` gives DNS, $\nu_t = 0$):

$$\nu_t = c_V \sqrt{\frac{B_\beta}{\alpha_{ij}\alpha_{ij}}}$$

where $\alpha_{ij} = \partial u_j / \partial x_i$, $\beta_{mn} = \sum_l \Delta_l^2\,\alpha_{lm}\alpha_{ln}$, $B_\beta = \beta_{11}\beta_{22} - \beta_{12}^2 + \beta_{11}\beta_{33} - \beta_{13}^2 + \beta_{22}\beta_{33} - \beta_{23}^2$, and $c_V = 2.5 C_s^2$.

All velocity gradients are evaluated at cell centres by second-order central differences; off-diagonal gradients are averaged from the two adjacent face-point values.

> **Reference:** Vreman, A.W. (2004). *An eddy-viscosity subgrid-scale model for turbulent shear flow: Algebraic theory and applications*. Phys. Fluids 16(10), 3670–3681. DOI: [10.1063/1.1785131](https://doi.org/10.1063/1.1785131)

---

## 6. Immersed boundary method (IBM)

Complex solid geometries (e.g. rough surfaces) are represented via a **ghost-cell immersed boundary method** following Tseng & Ferziger (2003).

### 6.1 Signed-distance function

The interface geometry is described by a signed-distance field $\phi$ at cell centres ($\phi < 0$ inside solid, $\phi > 0$ in fluid).  Two input modes are supported:

- **Mode 1** (`ibm_input_mode = 1`): read a binary face-point mask (`Umask_in`); $\phi$ is computed automatically by the fast-sweep algorithm of Zhao et al. (2005).
- **Mode 2** (`ibm_input_mode = 2`): read a precomputed cell-centre SDF (`SDF_in`) directly.

The `GenSDF` preprocessing tool in `preProcessing/GenSDF/` generates SDFs for user-defined geometries.

### 6.2 Ghost-cell interpolation

For each **ghost cell** $G$ (solid cell immediately adjacent to the fluid–solid interface):

1. Locate the boundary point $B$ at the zero crossing of $\phi$ along the wall-normal direction: $|GB| = |\phi_G|$.
2. Locate the image point $I$, the mirror of $G$ about $B$: $|BI| = |GB|$.
3. Interpolate the fluid velocity at $I$ using trilinear interpolation.
4. Set the ghost-cell velocity so that linear interpolation between $G$ and $I$ satisfies the boundary condition at $B$:

$$U(G) = 2 U_\text{wall} - U(I)$$

No-slip ($U_\text{wall} = 0$) is the default; moving-wall BCs can be set by changing `U_wall`, `V_wall`, `W_wall` in `ibm.f90`.

> **Reference:** Tseng, Y.-H. & Ferziger, J.H. (2003). *A ghost-cell immersed boundary method for flow in complex geometry*. J. Comput. Phys. 192(2), 593–623. DOI: [10.1016/j.jcp.2003.07.023](https://doi.org/10.1016/j.jcp.2003.07.023)

> **Reference (fast-sweep SDF):** Zhao, H., Osher, S. & Fedkiw, R. (2001 / 2005). *Fast surface reconstruction using the level set method*. Proc. IEEE ICCV; see also Tsai, Y.-H.R. (2002). *Rapid and accurate computation of the distance function using grids*. J. Comput. Phys. 178, 175–195. DOI: [10.1006/jcph.2002.7028](https://doi.org/10.1006/jcph.2002.7028)

---

## 7. Wall models

Both flat-wall and IBM-surface wall models use an **equilibrium wall model (EQWM)** based on the classical log-law.  The model is applied when `flat_wall_model_flag = 1` (flat walls) or `ibm_wall_model_flag = 1` (IBM surfaces).

### 7.1 Robin ghost-cell BC

The wall-normal velocity gradient at the boundary is related to the local wall-shear stress through a Robin (slip-length) condition parameterised by $\alpha$:

$$U_\text{ghost} = \alpha_y\,U_\text{interior} + (1-\alpha_y)\,U_\text{wall}$$

For no-slip ($\alpha_y = 0$) and free-slip ($\alpha_y = 1$) this reduces to the respective Dirichlet/Neumann limits.

### 7.2 Log-law Newton iteration

The friction velocity $u_\tau$ is obtained by a Newton iteration on the log-law / viscous-sublayer composite:

$$u^+ = \begin{cases} y^+ & y^+ < 5 \\ \frac{1}{\kappa}\ln(y^+) + B & y^+ \geq 5 \end{cases}$$

with $\kappa = 0.41$ and $B = 5.2$ (Prandtl–Kármán constants).  The iteration uses up to 20 Newton steps; $\alpha$ is then back-computed from $u_\tau$.

For IBM surfaces the reference velocity is the wall-tangential component of the velocity interpolated at the image-point location, projected onto the surface tangent plane.

---

## 8. Scalar transport (suspended sediment)

A passive scalar $C$ (suspended-sediment concentration) is advected and diffused on cell centres alongside the velocity field:

$$\frac{\partial C}{\partial t} + \frac{\partial (u_j C)}{\partial x_j} + \frac{\partial (w_s C)}{\partial y} = \frac{\partial}{\partial x_j}\left[\kappa_\text{eff}\frac{\partial C}{\partial x_j}\right]$$

where $\kappa_\text{eff} = \nu/Sc + \nu_t/Sc_t$ is the effective diffusivity and $w_s$ is the gravitational settling velocity (positive downward).

The equation is discretised in finite-volume form on cell centres.  Both the advective and diffusive flux divergences are normalised by the **local cell width** (the face-to-face spacing $x_{i}-x_{i-1}$, $y_{j}-y_{j-1}$, $z_{k}-z_{k-1}$), which guarantees discrete conservation on the stretched wall-normal grid and is consistent with the momentum solver.

### 8.1 Advection — van Leer slope-limited MUSCL

Advective face values are obtained from a MUSCL reconstruction of the upwind cell,

$$C_\text{face} = C_\text{up} + \sigma_\text{up}\,\bigl(x_\text{face} - x_{g,\text{up}}\bigr),$$

where the limited cell-centred slope $\sigma_\text{up}$ uses the **van Leer (1974) harmonic-mean limiter** applied to the forward and backward gradients formed with the *actual* cell-centre distances:

$$g_f = \frac{C_{m+1}-C_m}{x_{g,m+1}-x_{g,m}}, \quad g_b = \frac{C_m-C_{m-1}}{x_{g,m}-x_{g,m-1}}, \quad \sigma = \begin{cases}\dfrac{2\,g_f g_b}{g_f+g_b}, & g_f g_b > 0,\\[2mm] 0, & \text{otherwise.}\end{cases}$$

Because real geometric distances enter both the gradients and the projection to the face, the scheme remains **second-order accurate and TVD on non-uniform (stretched) grids** — unlike ratio-based limiters that assume uniform spacing.  The reconstruction is face-consistent (the value used as a cell's high face equals the neighbour's low face), so the discretisation is locally conservative.  At single-ghost boundaries (periodic seams, MPI slab interfaces, walls) where the far-upwind cell is unavailable, the stencil falls back to first-order upwind.

> **References:** van Leer, B. (1974). *Towards the ultimate conservative difference scheme. II. Monotonicity and conservation combined in a second-order scheme*. J. Comput. Phys. **14**(4), 361–370. DOI: [10.1016/0021-9991(74)90019-9](https://doi.org/10.1016/0021-9991(74)90019-9).  See also Toro, E.F. (2009). *Riemann Solvers and Numerical Methods for Fluid Dynamics*, 3rd ed., Springer, ch. 13–14, and LeVeque, R.J. (2002). *Finite Volume Methods for Hyperbolic Problems*, Cambridge University Press, ch. 6.

### 8.2 Settling velocity — Soulsby (1997)

The particle settling velocity is computed from the Soulsby (1997) formula:

$$D_* = d_s\left[\frac{(s-1)g}{\nu^2}\right]^{1/3}, \quad w_s = \frac{\nu}{d_s}\left[\sqrt{10.36^2 + 1.049\,D_*^3} - 10.36\right]$$

where $s = \rho_s/\rho_f$ is the specific gravity of sediment and $d_s$ is the particle diameter.

> **Reference:** Soulsby, R.L. (1997). *Dynamics of Marine Sands*. Thomas Telford, London. ISBN: 978-0-7277-2584-5

---

## 9. Reynolds stress budget

The `reynolds_stress_budget` module computes all terms in the **Reynolds stress transport equation** following Pope (2000) §7.4, Eq. 7.176:

$$\frac{DR_{ij}}{Dt} = P_{ij} + \Pi_{ij} + D^\nu_{ij} + D^T_{ij} + \Phi^P_{ij} - \varepsilon^\text{res}_{ij} - \varepsilon^\text{sgs}_{ij}$$

| Term | Symbol | Description |
|------|--------|-------------|
| Production | $P_{ij}$ | $-R_{ik}\partial\langle U_j\rangle/\partial x_k - R_{jk}\partial\langle U_i\rangle/\partial x_k$ |
| Pressure–strain | $\Pi_{ij}$ | $\langle p'(\partial u_i'/\partial x_j + \partial u_j'/\partial x_i)\rangle$ |
| Viscous diffusion | $D^\nu_{ij}$ | $\nu\,\nabla^2 R_{ij}$ |
| Turbulent diffusion | $D^T_{ij}$ | $-\partial\langle u_i' u_j' u_k'\rangle/\partial x_k$ |
| Pressure diffusion | $\Phi^P_{ij}$ | $-\partial(\langle p' u_i'\rangle\delta_{jk} + \langle p' u_j'\rangle\delta_{ik})/\partial x_k$ |
| Resolved dissipation | $\varepsilon^\text{res}_{ij}$ | $2\nu\langle(\partial u_i'/\partial x_k)(\partial u_j'/\partial x_k)\rangle$ |
| SGS dissipation | $\varepsilon^\text{sgs}_{ij}$ | $2\langle\nu_t s_{ij}'\rangle$ |

Statistics are accumulated as raw moments and converted to central moments at write time.  Homogeneous spatial averaging (configurable over $x$, $z$, or both) reduces the output to 1-D profiles for canonical channel flows.  Output files are little-endian float64, appended each window for restart continuity.

> **Reference:** Pope, S.B. (2000). *Turbulent Flows*. Cambridge University Press. ISBN: 978-0-521-59886-6. DOI: [10.1017/CBO9780511840531](https://doi.org/10.1017/CBO9780511840531)

---

## 10. MPI parallelism

The domain is decomposed into **z-slabs**: each MPI rank owns the full $(x, y)$ extent and a contiguous range of $z$-planes.  Ghost-cell exchanges in $z$ are performed with `MPI_Sendrecv` at the start of each RK3 sub-step.  The FFT-based Poisson solver uses FFTW's MPI transposed-output layout, which redistributes Fourier modes across ranks without additional all-to-all communication.

---

## 11. Streamwise inflow/outflow BC (`x_bc_type = 1`)

An alternative to the default periodic streamwise BC (`x_bc_type = 0`): a prescribed velocity is imposed at the inlet ($x=0$) and a convective condition advects the flow out at $x=L_x$. The pressure Poisson solve switches from a periodic FFT to a DCT-IV transform in $x$ (Section 4.1), consistent with the non-periodic velocity BC.

### 11.1 Inlet — constant or synthetic eddy method

`inflow_type = 0` imposes a uniform constant velocity. `inflow_type = 1` drives the inlet with the **Ensemble Synthetic Eddy Method (ESEM)** of Schau, Johnson, Muller & Oefelein (2022), an extension of the original SEM of Jarrin (2006): a population of eddies convects through a virtual box upstream of the inlet plane, and each eddy contributes a localised velocity-fluctuation kernel so that the superposed signal reproduces a target mean profile $\langle U(y)\rangle$ and Reynolds-stress profile $R_{ij}(y)$ read from `inflow_profile_file`. Relative to classical SEM, ESEM normalises the raw eddy-kernel sum empirically (sampled over an ensemble window) rather than analytically, giving exact one-point statistics regardless of eddy placement/length-scale assumptions. Optional extensions (`sem_sigma_file`, `sem_eddy_placement`, `sem_divergence_free`) add inhomogeneous eddy length scales, PDF-weighted placement, and a divergence-free (curl-of-vector-potential) construction after Poletto, Craft & Revell (2013). 

> **References:** Jarrin, N. (2006). *Synthetic inflow boundary conditions for the numerical simulation of turbulence*. PhD thesis, University of Manchester. Schau, K.A., Johnson, R.F., Muller, S. & Oefelein, J.C. (2022). *An ensemble Synthetic Eddy Method for accurate treatment of inhomogeneous turbulence*. Computers & Fluids 248, 105671. DOI: [10.1016/j.compfluid.2022.105671](https://doi.org/10.1016/j.compfluid.2022.105671). Poletto, R., Craft, T. & Revell, A. (2013). *A new divergence free synthetic eddy method for the reproduction of inlet flow conditions for LES*. Flow Turb. Combust. 91(3), 519–539. DOI: [10.1007/s10494-013-9488-2](https://doi.org/10.1007/s10494-013-9488-2)

### 11.2 Outlet — convective (advective) BC

At the outlet face, each velocity component is relaxed toward a one-sided upwind extrapolation from the interior at a rate set by the local Courant number:

$$F(n_\text{last}) \leftarrow F(n_\text{last}) - C\,\bigl[F(n_\text{last}) - F(n_\text{last}-1)\bigr], \qquad C = \min\!\bigl(\max(U_c,0),\,1\bigr)\frac{\Delta t}{\Delta x}$$

which discretises $\partial F/\partial t + U_c\,\partial F/\partial x = 0$ by first-order upwinding; clamping $C$ to $[0,1]$ keeps the update stable and avoids pulling from outside the domain under local backflow.

### 11.3 In-situ TI-profile rescaling (`ti_rescale_active = 1`)

Because the ESEM inlet Reynolds-stress profile decays/adjusts as it convects downstream before reaching equilibrium, an optional feedback loop samples the resolved turbulence intensity at a station `ti_rescale_x` downstream and nudges the *injected* target profile so the resolved profile matches the intended one at that station. Every `ti_rescale_freq` steps, raw velocity moments accumulated since the last window are reduced across ranks, converted to variances, and used to scale $R_{11}, R_{22}, R_{33}$ multiplicatively (damped by exponent `ti_rescale_relax`, clipped per-window to `[1/ti_rescale_clip, ti_rescale_clip]`), with an additional anti-windup clamp (`ti_rescale_abs_clip`) bounding the cumulative drift of the injected profile relative to its original target.

---

## 12. Initial conditions

| `ic_type` | Profile |
|-----------|---------|
| 1 | Log-law profile + white noise |
| 2 | Linear ramp (tent) + white noise |
| 3 | Zero mean + white noise |
| 4 | Reichardt (1951) composite law-of-the-wall + structured perturbation |
| 5 | Inverse-linear (anti-tent, promotes rapid transition) + white noise |

The Reichardt (1951) profile for IC type 4:

$$U^+(y^+) = \frac{1}{\kappa}\ln(1 + \kappa y^+) + 7.8\left[1 - e^{-y^+/11} - \frac{y^+}{11}e^{-0.33 y^+}\right]$$

> **Reference:** Reichardt, H. (1951). *Vollständige Darstellung der turbulenten Geschwindigkeitsverteilung in glatten Leitungen*. Z. Angew. Math. Mech. 31(7), 208–219. DOI: [10.1002/zamm.19510310704](https://doi.org/10.1002/zamm.19510310704)

---

## Acknowledgments

The original form of the channel flow solver was shared by Adrian Lozano-Duran during my PhD and I am grateful for his generous help and support with the code. 

