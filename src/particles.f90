!> Lagrangian point-particle tracking (Phase 1: tracer, one-way, full boundary handling).
!  Particle ownership follows the pencil (a rank owns a particle whose (x,z) falls in its
!  own cell-centre interior range, extended to +/-infinity at the true global x/z edges --
!  see owns_particle); y is never decomposed, so no y-ownership test is needed. Advances
!  once per full RK3 step against the frozen, already-consistent post-step U/V/W (see
!  advance_particles's caller in time_integration.f90), using a self-contained classical
!  RK3 for the tracer ODE dx/dt = u(x). Migration is a two-pass (x-then-z) dimension-split
!  exchange reusing decomp.f90's existing halo-neighbour/periodic-partner topology.
Module particles

  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use mpi
  Use decomp, Only : x_halo_neighbors, z_halo_neighbors, x_periodic_partner, z_periodic_partner
  Use synthetic_eddy_method, Only : random_seed_from

  Implicit None

  ! Structure-of-Arrays particle storage (per-rank; capacity grows by doubling). Particle IDs
  ! are Int32: seed particles get 0..n_particles_init-1 (identical candidate list generated on
  ! every rank, see seed_particles_fresh, so IDs are globally unique without communication);
  ! inflow-reinjected particles get myid*10_000_000 + a per-rank running counter, unique as
  ! long as n_particles_init and per-rank reinjection counts both stay under 10 million.
  Integer(Int32) :: n_particles_local = 0
  Integer(Int32) :: particle_capacity = 0
  Integer(Int32), Allocatable :: p_id(:)
  Real   (Int64), Allocatable :: p_x(:), p_y(:), p_z(:)
  Real   (Int64), Allocatable :: p_u(:), p_v(:), p_w(:)
  Real   (Int64), Allocatable :: p_age(:)
  ! Phase 5 (sgs_particle_model==1 only): persistent SGS velocity fluctuation state, evolved
  ! by the simplified Langevin model (advance_sgs_velocity) and added to the resolved fluid
  ! velocity used for advection
  Real   (Int64), Allocatable :: p_sgs_u(:), p_sgs_v(:), p_sgs_w(:)

  ! Resolved per-direction BC codes (0=periodic,1=exit,2=reflect,3=absorb) and this rank's
  ! position in the periodic topology, both set once in setup_particles
  Integer(Int32) :: bcx_resolved, bcy_resolved, bcz_resolved
  Logical        :: is_first_x_m, is_last_x_m, is_first_z_m, is_last_z_m
  Integer(Int32) :: partner_x_m, partner_z_m

  ! Per-rank event counters since the last report_particle_counts call
  Integer(Int64) :: n_exited_local = 0, n_deposited_local = 0, n_reinjected_local = 0
  Integer(Int32) :: n_exit_outflow_local = 0   ! outflow (is_last_x, x exit) exits this step, feeds reinject_at_inflow
  Integer(Int32) :: n_reinject_counter  = 0    ! per-rank running counter for reinjected-particle IDs

  ! Phase 2 (particle_mode==1 only): Stokes response time, computed once in setup_particles
  Real(Int64) :: tau_p_stokes

  ! Phase 4: per-rank streamwise deposition-rate accumulator, sized nxg_global (indexed by
  ! GLOBAL cell-centre index via bracket_index against xg_global -- see apply_ibm_collision).
  ! Every rank's array only ever gets incremented at indices inside its own owned x-range, so
  ! summing all ranks' copies (report_deposit_profile) gives the correct combined profile with
  ! no double-counting. Reset to 0 after each write, so the written value is a RATE per
  ! particle_deposit_freq*nmonitor steps, not a running total.
  Real(Int64), Allocatable :: particle_deposit_x(:)

Contains

  !> Seed (or restart-load) particles; no-op unless particles_active>=1. Called once from
  !  initialize() (initialization.f90), after the grid and x/z periodic-partner topology exist.
  Subroutine setup_particles

    Logical :: loaded

    If ( particles_active < 1 ) Return

    Call x_periodic_partner(is_first_x_m, is_last_x_m, partner_x_m)
    Call z_periodic_partner(is_first_z_m, is_last_z_m, partner_z_m)

    bcx_resolved = bc_particle_x
    bcy_resolved = bc_particle_y
    bcz_resolved = bc_particle_z

    ! Stokes response time tau_p = rho_p*d_p^2/(18*mu_f), mu_f = rho_f*nu (particle_mode==1 only,
    ! but harmless to precompute unconditionally)
    tau_p_stokes = particle_rho * particle_diam**2 / (18d0 * particle_rho_f * nu)

    n_particles_local = 0
    particle_capacity = 0
    n_reinject_counter = 0

    If ( ibm_input_mode >= 1 ) Then
       Allocate ( particle_deposit_x(nxg_global) )
       particle_deposit_x = 0d0
    End If

    loaded = .False.
    If ( restart == 1 ) Call read_particle_restart(loaded)
    If ( .Not. loaded ) Call seed_particles_fresh

    If ( myid == 0 ) Write(*,'(A,I8,A)') ' particles: ', n_particles_init, ' seed particles configured (tracer, one-way)'

  End Subroutine setup_particles

  !> Uniform-random seed within the particle_seed_[xyz]min/max box: every rank draws the
  !  SAME n_particles_init candidates (identical RNG seed/sequence), then keeps only the
  !  ones it owns -- avoids any seed-time communication while keeping candidate index i-1
  !  (globally unique, no rank collision) as each kept particle's ID.
  Subroutine seed_particles_fresh

    Integer(Int32) :: i
    Real   (Int64) :: xp, yp, zp, rx, ry, rz

    Call random_seed_from(particle_seed_seed)

    Do i = 1, n_particles_init
       Call random_number(rx)
       Call random_number(ry)
       Call random_number(rz)
       xp = particle_seed_xmin + rx*(particle_seed_xmax - particle_seed_xmin)
       yp = particle_seed_ymin + ry*(particle_seed_ymax - particle_seed_ymin)
       zp = particle_seed_zmin + rz*(particle_seed_zmax - particle_seed_zmin)
       If ( owns_particle(xp, zp) ) Then
          Call ensure_capacity(n_particles_local + 1)
          n_particles_local = n_particles_local + 1
          p_id (n_particles_local) = i - 1
          p_x  (n_particles_local) = xp
          p_y  (n_particles_local) = yp
          p_z  (n_particles_local) = zp
          p_u  (n_particles_local) = 0d0
          p_v  (n_particles_local) = 0d0
          p_w  (n_particles_local) = 0d0
          p_age(n_particles_local) = 0d0
          If ( sgs_particle_model == 1 ) Then
             p_sgs_u(n_particles_local) = 0d0
             p_sgs_v(n_particles_local) = 0d0
             p_sgs_w(n_particles_local) = 0d0
          End If
       End If
    End Do

  End Subroutine seed_particles_fresh

  !> True if this rank owns (xp,zp): the half-open range [xg(2),xg(nxg))/[zg(2),zg(nzg)),
  !  extended to +/-infinity at the true global x/z edges (is_first/is_last) so seeding/
  !  restart never drops an out-of-domain point. The upper bound is the GHOST coordinate
  !  (nxg/nzg), not the last true-interior one (nxg-1/nzg-1): local xg(nxg) on rank A is
  !  exactly the same physical point as local xg(2) on A's right neighbour (both slice
  !  xg_global via kg1/kg2_global, decomp.f90 -- kg2_global(A)=kg1_global(A_right)+1, so
  !  A's ghost coordinate and the neighbour's first interior coordinate coincide). Using
  !  the last-true-interior point as an exclusive bound instead leaves a one-grid-cell-wide
  !  band at every internal seam that neither rank claims -- caught empirically: a 4-rank
  !  z-split smoke test lost ~11% of seeded particles (~1 cell out of ~8 per rank at each
  !  of 3 seams) before this fix.
  Logical Function owns_particle(xp, zp) Result(owns)

    Real(Int64), Intent(In) :: xp, zp
    Logical :: owns_x, owns_z

    owns_x = ( xp >= xg(2) .Or. is_first_x_m ) .And. ( xp < xg(nxg) .Or. is_last_x_m )
    owns_z = ( zp >= zg(2) .Or. is_first_z_m ) .And. ( zp < zg(nzg) .Or. is_last_z_m )
    owns = owns_x .And. owns_z

  End Function owns_particle

  !> Binary search for the bracket lo such that arr(lo)<=pos<arr(lo+1) (clamped at the ends),
  !  for a monotonically increasing arr. Same structure as sem.f90's invert_cdf and
  !  uav_actuator.f90's nearest_zg_index; reused here for every interpolation axis (x,y,z,
  !  whichever grid -- face or cell-centre -- a given field component is staggered on).
  Integer(Int32) Function bracket_index(arr, n, pos) Result(lo)

    Integer(Int32), Intent(In) :: n
    Real   (Int64), Intent(In) :: arr(n)
    Real   (Int64), Intent(In) :: pos
    Integer(Int32) :: hi, m

    If ( pos <= arr(1) ) Then
       lo = 1
       Return
    Else If ( pos >= arr(n) ) Then
       lo = n - 1
       Return
    End If

    lo = 1
    hi = n - 1
    Do While ( lo < hi )
       m = (lo + hi) / 2
       If ( arr(m+1) < pos ) Then
          lo = m + 1
       Else
          hi = m
       End If
    End Do

  End Function bracket_index

  !> Trilinear interpolation of F at (xp,yp,zp), where F is staggered on (ax,ay,az) --
  !  whichever face/cell-centre grid arrays the caller passes for F's own staggering.
  Real(Int64) Function interp3(F, ax, na, ay, nb, az, nc, xp, yp, zp) Result(val)

    Integer(Int32), Intent(In) :: na, nb, nc
    Real   (Int64), Intent(In) :: F(na,nb,nc)
    Real   (Int64), Intent(In) :: ax(na), ay(nb), az(nc)
    Real   (Int64), Intent(In) :: xp, yp, zp

    Integer(Int32) :: i0, j0, k0
    Real   (Int64) :: fx, fy, fz

    i0 = bracket_index(ax, na, xp)
    fx = Min( Max( (xp-ax(i0)) / (ax(i0+1)-ax(i0)), 0d0 ), 1d0 )
    j0 = bracket_index(ay, nb, yp)
    fy = Min( Max( (yp-ay(j0)) / (ay(j0+1)-ay(j0)), 0d0 ), 1d0 )
    k0 = bracket_index(az, nc, zp)
    fz = Min( Max( (zp-az(k0)) / (az(k0+1)-az(k0)), 0d0 ), 1d0 )

    val = (1d0-fx)*(1d0-fy)*(1d0-fz)*F(i0  ,j0  ,k0  ) + &
          (   fx)*(1d0-fy)*(1d0-fz)*F(i0+1,j0  ,k0  ) + &
          (1d0-fx)*(   fy)*(1d0-fz)*F(i0  ,j0+1,k0  ) + &
          (1d0-fx)*(1d0-fy)*(   fz)*F(i0  ,j0  ,k0+1) + &
          (   fx)*(   fy)*(1d0-fz)*F(i0+1,j0+1,k0  ) + &
          (   fx)*(1d0-fy)*(   fz)*F(i0+1,j0  ,k0+1) + &
          (1d0-fx)*(   fy)*(   fz)*F(i0  ,j0+1,k0+1) + &
          (   fx)*(   fy)*(   fz)*F(i0+1,j0+1,k0+1)

  End Function interp3

  !> Fluid velocity at (xp,yp,zp), one trilinear read per staggered component (U on x-face/
  !  yg/zg centres, V on xg/y-face/zg, W on xg/yg/z-face -- this solver's usual staggering).
  Subroutine interpolate_velocity(xp, yp, zp, up, vp, wp)

    Real(Int64), Intent(In)  :: xp, yp, zp
    Real(Int64), Intent(Out) :: up, vp, wp

    up = interp3(U, x,  nx,  yg, nyg, zg, nzg, xp, yp, zp)
    vp = interp3(V, xg, nxg, y,  ny,  zg, nzg, xp, yp, zp)
    wp = interp3(W, xg, nxg, yg, nyg, z,  nz,  xp, yp, zp)

  End Subroutine interpolate_velocity

  !> Advance local particle i one step under the inertial (Maxey-Riley-reduced) model:
  !  dv/dt = (u_f-v)/tau_eff + a_other, with u_f, tau_eff (nonlinear Schiller-Naumann drag
  !  correction) and a_other (gravity + optional added-mass) all frozen at their start-of-
  !  step values -- this makes the ODE linear over the step, so it has an EXACT exponential
  !  solution with no sub-cycling needed even for tau_eff << dt (small-Stokes-number
  !  particles). Optional Brownian motion (particle_brownian==1) is added as an independent
  !  random walk on top of this deterministic update. Saffman-Mei lift is deferred to
  !  Phase 3 (needs the IBM wall-distance field to gate it near-wall).
  Subroutine advance_particle_inertial(i)

    Integer(Int32), Intent(In) :: i

    Real(Int64) :: uf, vf, wf, uf_o, vf_o, wf_o
    Real(Int64) :: urel, Re_p, f_drag, tau_eff
    Real(Int64) :: ax, ay, az, rho_f_local
    Real(Int64) :: ex, ey, ez, dv, dsig

    Call interpolate_velocity(p_x(i), p_y(i), p_z(i), uf, vf, wf)
    If ( sgs_particle_model == 1 ) Then
       uf = uf + p_sgs_u(i);  vf = vf + p_sgs_v(i);  wf = wf + p_sgs_w(i)
    End If

    urel = Sqrt( (uf-p_u(i))**2 + (vf-p_v(i))**2 + (wf-p_w(i))**2 )
    Re_p = urel * particle_diam / nu
    f_drag = 1d0 + 0.15d0 * Re_p**0.687d0            ! Schiller-Naumann, valid Re_p ~< 1000
    tau_eff = tau_p_stokes / f_drag

    rho_f_local = particle_rho_f
    If ( particle_boussinesq_coupling == 1 .And. boussinesq_flag >= 1 ) Then
       ! local Boussinesq fluid density at the particle, buoyancy term only (see global.f90)
       rho_f_local = particle_rho_f * ( 1d0 - beta_T * &
            (interp3(Tscal, xg, nxg, yg, nyg, zg, nzg, p_x(i), p_y(i), p_z(i)) - T_ref) )
    End If
    ax = 0d0;  ay = -grav*(1d0 - rho_f_local/particle_rho);  az = 0d0   ! gravity acts along -y (Boussinesq convention)

    If ( particle_added_mass == 1 ) Then
       ! Simplified added-mass: local Eulerian dU/dt at the particle's own position (Uo/Vo/Wo
       ! are the previous RK-step's start-of-step field, global.f90) -- omits the convective
       ! (u.grad)u term of the true material derivative, a deliberate simplification to avoid
       ! an expensive on-the-fly velocity-gradient evaluation at an arbitrary particle position.
       Call interpolate_velocity_old(p_x(i), p_y(i), p_z(i), uf_o, vf_o, wf_o)
       ax = ax + 0.5d0*(particle_rho_f/particle_rho) * (uf-uf_o)/Max(dt,1d-14)
       ay = ay + 0.5d0*(particle_rho_f/particle_rho) * (vf-vf_o)/Max(dt,1d-14)
       az = az + 0.5d0*(particle_rho_f/particle_rho) * (wf-wf_o)/Max(dt,1d-14)
    End If

    Call exp_integrate(p_x(i), p_u(i), uf, ax, tau_eff, dt)
    Call exp_integrate(p_y(i), p_v(i), vf, ay, tau_eff, dt)
    Call exp_integrate(p_z(i), p_w(i), wf, az, tau_eff, dt)

    If ( particle_brownian == 1 ) Then
       ! Isotropic Stokes-Einstein Brownian kick: D = kB*T/(3*pi*mu_f*d_p), displacement
       ! std dev = sqrt(2*D*dt) per axis (independent Box-Muller draws)
       dsig = Sqrt( 2d0 * (1.380649d-23*particle_temp_abs) / &
            (3d0*pi*particle_rho_f*nu*particle_diam) * dt )
       Call gaussian_pair(ex, ey);  Call gaussian_pair(ez, dv)
       p_x(i) = p_x(i) + dsig*ex
       p_y(i) = p_y(i) + dsig*ey
       p_z(i) = p_z(i) + dsig*ez
    End If

  End Subroutine advance_particle_inertial

  !> Exact solution of dv/dt=(u-v)/tau+a (u,a frozen) over [0,dt], advancing both x and v in place.
  Subroutine exp_integrate(x, v, u, a, tau, dt_in)

    Real(Int64), Intent(InOut) :: x, v
    Real(Int64), Intent(In)    :: u, a, tau, dt_in

    Real(Int64) :: vterm, dv0, edt

    vterm = u + a*tau
    dv0   = v - vterm
    edt   = Exp(-dt_in/Max(tau,1d-14))

    x = x + vterm*dt_in + tau*dv0*(1d0 - edt)
    v = vterm + dv0*edt

  End Subroutine exp_integrate

  !> Fluid velocity at (xp,yp,zp) from the PREVIOUS step's start-of-step field (Uo,Vo,Wo,
  !  global.f90) -- same staggering/interpolation as interpolate_velocity, used only by
  !  the simplified added-mass term above.
  Subroutine interpolate_velocity_old(xp, yp, zp, up, vp, wp)

    Real(Int64), Intent(In)  :: xp, yp, zp
    Real(Int64), Intent(Out) :: up, vp, wp

    up = interp3(Uo, x,  nx,  yg, nyg, zg, nzg, xp, yp, zp)
    vp = interp3(Vo, xg, nxg, y,  ny,  zg, nzg, xp, yp, zp)
    wp = interp3(Wo, xg, nxg, yg, nyg, z,  nz,  xp, yp, zp)

  End Subroutine interpolate_velocity_old

  !> Two independent standard-normal draws via Box-Muller (Fortran's Random_Number is uniform only)
  Subroutine gaussian_pair(g1, g2)

    Real(Int64), Intent(Out) :: g1, g2
    Real(Int64) :: u1, u2, r

    Call random_number(u1);  Call random_number(u2)
    u1 = Max(u1, 1d-300)   ! avoid Log(0)
    r  = Sqrt(-2d0*Log(u1))
    g1 = r*Cos(2d0*pi*u2)
    g2 = r*Sin(2d0*pi*u2)

  End Subroutine gaussian_pair

  !> Diagnose SGS turbulent kinetic energy and dissipation rate at (xp,yp,zp) from the eddy-
  !  viscosity SGS model's nu_t (sgs_model.f90/global.f90), via a Deardorff-style mixing-length
  !  closure nu_t=C_k*sqrt(k_sgs)*Delta => k_sgs=(nu_t/(C_k*Delta))^2, and eps_sgs=C_eps*
  !  k_sgs^1.5/Delta. C_k=0.1, C_eps=1.0 are standard textbook LES values, not case-tuned.
  Subroutine compute_sgs_stats(xp, yp, zp, k_sgs, eps_sgs)

    Real(Int64), Intent(In)  :: xp, yp, zp
    Real(Int64), Intent(Out) :: k_sgs, eps_sgs

    Real(Int64) :: nut_local
    Real(Int64), Parameter :: C_k = 0.1d0, C_eps = 1.0d0

    nut_local = interp3(nu_t, xg, nxg, yg, nyg, zg, nzg, xp, yp, zp)
    k_sgs   = ( Max(nut_local,0d0) / (C_k*Max(Delta,1d-14)) )**2
    eps_sgs = C_eps * k_sgs**1.5d0 / Max(Delta,1d-14)

  End Subroutine compute_sgs_stats

  !> Explicit-Euler step of the simplified (isotropic) Langevin SGS model (simplified Thomson/
  !  Weil-Sullivan-Moeng well-mixed model): du_sgs = -u_sgs/T_L*dt + sqrt(C0*eps_sgs)*dW, with
  !  T_L=2*k_sgs/(C0*eps_sgs). Assumes dt<<T_L (typical for LES time steps); a documented
  !  simplification rather than the exact-exponential treatment used for drag (exp_integrate).
  Subroutine advance_sgs_velocity(i)

    Integer(Int32), Intent(In) :: i

    Real(Int64) :: k_sgs, eps_sgs, T_L, b, g1, g2, g3, g4

    Call compute_sgs_stats(p_x(i), p_y(i), p_z(i), k_sgs, eps_sgs)
    If ( k_sgs < 1d-20 .Or. eps_sgs < 1d-20 ) Then
       p_sgs_u(i) = 0d0;  p_sgs_v(i) = 0d0;  p_sgs_w(i) = 0d0
       Return
    End If

    T_L = 2d0*k_sgs / (particle_langevin_C0*eps_sgs)
    b   = Sqrt(particle_langevin_C0*eps_sgs)

    Call gaussian_pair(g1, g2)
    Call gaussian_pair(g3, g4)
    p_sgs_u(i) = p_sgs_u(i)*(1d0 - dt/T_L) + b*Sqrt(dt)*g1
    p_sgs_v(i) = p_sgs_v(i)*(1d0 - dt/T_L) + b*Sqrt(dt)*g2
    p_sgs_w(i) = p_sgs_w(i)*(1d0 - dt/T_L) + b*Sqrt(dt)*g3

  End Subroutine advance_sgs_velocity

  !> Grow every per-particle array to at least n_needed (doubling), preserving existing data.
  Subroutine ensure_capacity(n_needed)

    Integer(Int32), Intent(In) :: n_needed
    Integer(Int32) :: new_cap

    If ( n_needed <= particle_capacity ) Return
    new_cap = Max(n_needed, Max(2*particle_capacity, 16))

    Call grow_int (p_id,  new_cap)
    Call grow_real(p_x,   new_cap)
    Call grow_real(p_y,   new_cap)
    Call grow_real(p_z,   new_cap)
    Call grow_real(p_u,   new_cap)
    Call grow_real(p_v,   new_cap)
    Call grow_real(p_w,   new_cap)
    Call grow_real(p_age, new_cap)
    If ( sgs_particle_model == 1 ) Then
       Call grow_real(p_sgs_u, new_cap)
       Call grow_real(p_sgs_v, new_cap)
       Call grow_real(p_sgs_w, new_cap)
    End If
    particle_capacity = new_cap

  End Subroutine ensure_capacity

  Subroutine grow_real(arr, new_cap)
    Real(Int64), Allocatable, Intent(InOut) :: arr(:)
    Integer(Int32), Intent(In) :: new_cap
    Real(Int64), Allocatable :: tmp(:)
    Allocate(tmp(new_cap))
    tmp = 0d0
    If ( Allocated(arr) ) tmp(1:Size(arr)) = arr
    Call Move_Alloc(tmp, arr)
  End Subroutine grow_real

  Subroutine grow_int(arr, new_cap)
    Integer(Int32), Allocatable, Intent(InOut) :: arr(:)
    Integer(Int32), Intent(In) :: new_cap
    Integer(Int32), Allocatable :: tmp(:)
    Allocate(tmp(new_cap))
    tmp = 0
    If ( Allocated(arr) ) tmp(1:Size(arr)) = arr
    Call Move_Alloc(tmp, arr)
  End Subroutine grow_int

  !> Remove local particle i by swapping the last particle into its slot (order doesn't
  !  matter physically) and shrinking the count; O(1), no compaction pass needed.
  Subroutine remove_particle(i)

    Integer(Int32), Intent(In) :: i
    Integer(Int32) :: n

    n = n_particles_local
    If ( i /= n ) Then
       p_id (i) = p_id (n)
       p_x  (i) = p_x  (n)
       p_y  (i) = p_y  (n)
       p_z  (i) = p_z  (n)
       p_u  (i) = p_u  (n)
       p_v  (i) = p_v  (n)
       p_w  (i) = p_w  (n)
       p_age(i) = p_age(n)
       If ( sgs_particle_model == 1 ) Then
          p_sgs_u(i) = p_sgs_u(n)
          p_sgs_v(i) = p_sgs_v(n)
          p_sgs_w(i) = p_sgs_w(n)
       End If
    End If
    n_particles_local = n - 1

  End Subroutine remove_particle

  !> Per-particle boundary handling on the just-updated position, applied per direction
  !  (x,y,z independently); returns .True. if the particle should be removed (exit/absorb/
  !  aged-out). Reflect mirrors the position and flips the matching velocity component;
  !  periodic (case 0) does nothing here -- the true-global-edge crossing is detected and
  !  handled (with wrap) by exchange_particles_1d, which alone knows the periodic partner.
  Logical Function apply_particle_bc(i) Result(do_remove)

    Integer(Int32), Intent(In) :: i

    do_remove = .False.

    If ( p_age(i) > particle_max_age ) Then
       do_remove = .True.
       n_exited_local = n_exited_local + 1
       Return
    End If

    ! x
    Select Case (bcx_resolved)
    Case (1)
       If ( p_x(i) < x_global(1) .Or. p_x(i) >= x_global(nx_global) ) Then
          do_remove = .True.
          n_exited_local = n_exited_local + 1
          If ( is_last_x_m .And. p_x(i) >= x_global(nx_global) ) n_exit_outflow_local = n_exit_outflow_local + 1
          Return
       End If
    Case (2)
       If ( p_x(i) < x_global(1) ) Then
          p_x(i) = 2d0*x_global(1) - p_x(i);  p_u(i) = -p_u(i)
       Else If ( p_x(i) >= x_global(nx_global) ) Then
          p_x(i) = 2d0*x_global(nx_global) - p_x(i);  p_u(i) = -p_u(i)
       End If
    Case (3)
       If ( p_x(i) < x_global(1) .Or. p_x(i) >= x_global(nx_global) ) Then
          do_remove = .True.
          n_deposited_local = n_deposited_local + 1
          Return
       End If
    End Select

    ! y
    Select Case (bcy_resolved)
    Case (1)
       If ( p_y(i) < y_global(1) .Or. p_y(i) >= y_global(ny_global) ) Then
          do_remove = .True.
          n_exited_local = n_exited_local + 1
          Return
       End If
    Case (2)
       If ( p_y(i) < y_global(1) ) Then
          p_y(i) = 2d0*y_global(1) - p_y(i);  p_v(i) = -p_v(i)
       Else If ( p_y(i) >= y_global(ny_global) ) Then
          p_y(i) = 2d0*y_global(ny_global) - p_y(i);  p_v(i) = -p_v(i)
       End If
    Case (3)
       If ( p_y(i) < y_global(1) .Or. p_y(i) >= y_global(ny_global) ) Then
          do_remove = .True.
          n_deposited_local = n_deposited_local + 1
          Return
       End If
    End Select

    ! z
    Select Case (bcz_resolved)
    Case (1)
       If ( p_z(i) < z_global(1) .Or. p_z(i) >= z_global(nz_global) ) Then
          do_remove = .True.
          n_exited_local = n_exited_local + 1
          Return
       End If
    Case (2)
       If ( p_z(i) < z_global(1) ) Then
          p_z(i) = 2d0*z_global(1) - p_z(i);  p_w(i) = -p_w(i)
       Else If ( p_z(i) >= z_global(nz_global) ) Then
          p_z(i) = 2d0*z_global(nz_global) - p_z(i);  p_w(i) = -p_w(i)
       End If
    Case (3)
       If ( p_z(i) < z_global(1) .Or. p_z(i) >= z_global(nz_global) ) Then
          do_remove = .True.
          n_deposited_local = n_deposited_local + 1
          Return
       End If
    End Select

  End Function apply_particle_bc

  !> IBM/SDF collision handling (Phase 3, ibm_input_mode>=1 only): phi/Umask_cc/ibm_obj_id
  !  (global.f90, populated by ibm.f90's setup_ibm) are read directly -- no `Use ibm` needed,
  !  which keeps particles.f90 usable by dopamine-ESEM (which never compiles ibm.f90).
  !  phi<=0 means the particle has penetrated a solid; the surface normal is the (finite-
  !  difference) gradient of the interpolated phi field, and the particle is pushed back out
  !  along it by approximately its own penetration depth (phi is a signed distance field).
  !  Per-object BC (particle_ibm_bc, global.f90) selects absorb/reflect/deposit_resuspend;
  !  reflect uses a simple Stokes-number-dependent restitution heuristic (see global.f90's
  !  comment on particle_ibm_tau_crit -- a documented simplification, not a validated
  !  closed-form model). Returns .True. if the particle should be removed (absorbed/deposited).
  Logical Function apply_ibm_collision(i) Result(do_remove)

    Integer(Int32), Intent(In) :: i

    Real(Int64) :: phi_c, eps, gx, gy, gz, gmag, vn, e_rest, urel, obj_r
    Integer(Int32) :: obj_id, bc_code

    do_remove = .False.

    phi_c = interp3(phi, xg, nxg, yg, nyg, zg, nzg, p_x(i), p_y(i), p_z(i))
    If ( phi_c > 0d0 ) Return   ! still on the fluid side, nothing to do

    eps = 0.5d0*dxmin
    gx = interp3(phi, xg, nxg, yg, nyg, zg, nzg, p_x(i)+eps, p_y(i), p_z(i)) - &
         interp3(phi, xg, nxg, yg, nyg, zg, nzg, p_x(i)-eps, p_y(i), p_z(i))
    gy = interp3(phi, xg, nxg, yg, nyg, zg, nzg, p_x(i), p_y(i)+eps, p_z(i)) - &
         interp3(phi, xg, nxg, yg, nyg, zg, nzg, p_x(i), p_y(i)-eps, p_z(i))
    gz = interp3(phi, xg, nxg, yg, nyg, zg, nzg, p_x(i), p_y(i), p_z(i)+eps) - &
         interp3(phi, xg, nxg, yg, nyg, zg, nzg, p_x(i), p_y(i), p_z(i)-eps)
    gmag = Sqrt(gx*gx + gy*gy + gz*gz)
    If ( gmag < 1d-14 ) Return   ! degenerate gradient (shouldn't happen inside a real solid); leave the particle as-is rather than divide by ~0
    gx = gx/gmag;  gy = gy/gmag;  gz = gz/gmag   ! unit normal, points from solid toward fluid

    obj_r = interp3(ibm_obj_id, xg, nxg, yg, nyg, zg, nzg, p_x(i), p_y(i), p_z(i))
    obj_id = Max(0, Min(max_ibm_objects, Nint(obj_r)))
    bc_code = particle_ibm_bc(obj_id)

    If ( bc_code == 3 ) Then
       ! deposit_resuspend: resuspend (treat as reflect) only if the local relative speed
       ! exceeds particle_resuspend_ucrit, otherwise deposit -- simplified Shields/van Rijn proxy
       urel = Sqrt( p_u(i)**2 + p_v(i)**2 + p_w(i)**2 )
       bc_code = Merge(2, 1, urel > particle_resuspend_ucrit)
    End If

    If ( bc_code == 1 ) Then
       do_remove = .True.
       n_deposited_local = n_deposited_local + 1
       particle_deposit_x(bracket_index(xg_global, nxg_global, p_x(i))) = &
            particle_deposit_x(bracket_index(xg_global, nxg_global, p_x(i))) + 1d0
       Return
    End If

    ! reflect (bc_code==2, the default): push back out along the normal by the penetration
    ! depth (+ a small margin so it doesn't land exactly on phi=0), and mirror the velocity's
    ! normal component with the Stokes-dependent restitution coefficient
    p_x(i) = p_x(i) + (-phi_c + eps)*gx
    p_y(i) = p_y(i) + (-phi_c + eps)*gy
    p_z(i) = p_z(i) + (-phi_c + eps)*gz

    e_rest = Min(1d0, tau_p_stokes/particle_ibm_tau_crit)
    vn = p_u(i)*gx + p_v(i)*gy + p_w(i)*gz
    p_u(i) = p_u(i) - (1d0+e_rest)*vn*gx
    p_v(i) = p_v(i) - (1d0+e_rest)*vn*gy
    p_w(i) = p_w(i) - (1d0+e_rest)*vn*gz

  End Function apply_ibm_collision

  !> Advance every local particle one full step against the already-consistent post-
  !  projection U/V/W (frozen for the step), apply boundary handling, migrate, and (at
  !  the monitor cadence) report counters. Called once per step from time_integration.f90's
  !  compute_time_step_RK3, after the final host U/V/W sync. particle_mode==0 (tracer):
  !  classical 3-stage RK on dx/dt=u_fluid(x), Phase 1 behaviour. particle_mode==1
  !  (inertial): see advance_particle_inertial.
  Subroutine advance_particles

    Integer(Int32) :: i
    Real(Int64) :: x0, y0, z0
    Real(Int64) :: k1u, k1v, k1w, k2u, k2v, k2w, k3u, k3v, k3w
    Logical :: do_remove

    If ( particles_active < 1 ) Return

    n_exit_outflow_local = 0

    i = 1
    Do While ( i <= n_particles_local )

       If ( sgs_particle_model == 1 ) Call advance_sgs_velocity(i)

       If ( particle_mode == 1 ) Then
          Call advance_particle_inertial(i)
       Else
          x0 = p_x(i);  y0 = p_y(i);  z0 = p_z(i)

          Call interpolate_velocity(x0, y0, z0, k1u, k1v, k1w)
          Call interpolate_velocity(x0 + 0.5d0*dt*k1u, y0 + 0.5d0*dt*k1v, z0 + 0.5d0*dt*k1w, k2u, k2v, k2w)
          Call interpolate_velocity(x0 - dt*k1u + 2d0*dt*k2u, y0 - dt*k1v + 2d0*dt*k2v, z0 - dt*k1w + 2d0*dt*k2w, k3u, k3v, k3w)

          ! SGS fluctuation held frozen across the 3 RK evaluations (see advance_sgs_velocity)
          If ( sgs_particle_model == 1 ) Then
             k1u = k1u+p_sgs_u(i);  k1v = k1v+p_sgs_v(i);  k1w = k1w+p_sgs_w(i)
             k2u = k2u+p_sgs_u(i);  k2v = k2v+p_sgs_v(i);  k2w = k2w+p_sgs_w(i)
             k3u = k3u+p_sgs_u(i);  k3v = k3v+p_sgs_v(i);  k3w = k3w+p_sgs_w(i)
          End If

          p_x(i) = x0 + dt/6d0*(k1u + 4d0*k2u + k3u)
          p_y(i) = y0 + dt/6d0*(k1v + 4d0*k2v + k3v)
          p_z(i) = z0 + dt/6d0*(k1w + 4d0*k2w + k3w)
          p_u(i) = k1u;  p_v(i) = k1v;  p_w(i) = k1w
       End If
       p_age(i) = p_age(i) + dt

       ! IBM collision first (Phase 3): may reflect/push the position back to the fluid side
       ! before the domain-boundary check below sees it
       do_remove = .False.
       If ( ibm_input_mode >= 1 ) do_remove = apply_ibm_collision(i)
       If ( .Not. do_remove .And. apply_particle_bc(i) ) do_remove = .True.

       If ( do_remove ) Then
          Call remove_particle(i)   ! swaps the last particle into slot i; re-check the same slot
       Else
          i = i + 1
       End If
    End Do

    Call migrate_particles

    If ( particle_reinit_on_exit == 1 .And. bcx_resolved == 1 ) Call reinject_at_inflow

    If ( Mod(istep, Max(nmonitor,1)) == 0 ) Call report_particle_counts

    If ( ibm_input_mode >= 1 .And. Mod(istep, Max(nmonitor*particle_deposit_freq,1)) == 0 ) &
         Call write_deposit_profile

  End Subroutine advance_particles

  Subroutine migrate_particles

    Call exchange_particles_1d(1)   ! x
    Call exchange_particles_1d(2)   ! z

  End Subroutine migrate_particles

  !> Two-directional (up/down) non-blocking-free (Mpi_sendrecv) particle exchange along one
  !  axis. On the true global edge (is_first/is_last), the "missing" neighbour side (which
  !  x_halo_neighbors/z_halo_neighbors report as MPI_PROC_NULL) is only used at all when that
  !  direction is periodic, in which case it IS the periodic partner and the crossing
  !  particle's coordinate is unwrapped (+/- domain length) at hand-off -- this makes the true
  !  global edge a same-rank-count "6th neighbour" for exactly the ranks that touch it, with
  !  no relay/multi-hop needed. Non-periodic true-edge crossings never reach here: they were
  !  already resolved (exit/absorb/reflect) by apply_particle_bc.
  Subroutine exchange_particles_1d(dir)

    Integer(Int32), Intent(In) :: dir

    Integer(Int32) :: up, down
    Integer(Int32) :: i, n_up, n_down, n_recv_up, n_recv_down
    Integer(Int32), Allocatable :: id_up(:), id_down(:), id_ru(:), id_rd(:)
    Real   (Int64), Allocatable :: dat_up(:,:), dat_down(:,:), dat_ru(:,:), dat_rd(:,:)
    Logical :: goes_up, goes_down, is_first_m, is_last_m, periodic
    Integer(Int32) :: partner
    Real   (Int64) :: pos, domain_len, xp_send

    If ( dir == 1 ) Then
       Call x_halo_neighbors(up, down)
       is_first_m = is_first_x_m;  is_last_m = is_last_x_m;  partner = partner_x_m
       periodic = ( bcx_resolved == 0 )
       domain_len = x_global(nx_global) - x_global(1)
    Else
       Call z_halo_neighbors(up, down)
       is_first_m = is_first_z_m;  is_last_m = is_last_z_m;  partner = partner_z_m
       periodic = ( bcz_resolved == 0 )
       domain_len = z_global(nz_global) - z_global(1)
    End If
    If ( is_last_m  .And. periodic ) up   = partner
    If ( is_first_m .And. periodic ) down = partner

    Allocate ( id_up(Max(n_particles_local,1)), id_down(Max(n_particles_local,1)) )
    Allocate ( dat_up(7,Max(n_particles_local,1)), dat_down(7,Max(n_particles_local,1)) )
    n_up = 0;  n_down = 0

    i = 1
    Do While ( i <= n_particles_local )
       If ( dir == 1 ) Then
          pos = p_x(i)
       Else
          pos = p_z(i)
       End If

       If ( is_last_m ) Then
          goes_up = ( periodic .And. pos >= x_edge_hi(dir) )
       Else
          goes_up = ( pos >= xg_edge_hi(dir) )
       End If
       If ( is_first_m ) Then
          goes_down = ( periodic .And. pos < x_edge_lo(dir) )
       Else
          goes_down = ( pos < xg_edge_lo(dir) )
       End If

       If ( goes_up ) Then
          n_up = n_up + 1
          xp_send = pos
          If ( is_last_m ) xp_send = xp_send - domain_len   ! unwrap on hand-off to the periodic partner
          Call pack_particle(i, dir, xp_send, id_up(n_up), dat_up(:,n_up))
          Call remove_particle(i)
       Else If ( goes_down ) Then
          n_down = n_down + 1
          xp_send = pos
          If ( is_first_m ) xp_send = xp_send + domain_len   ! unwrap on hand-off to the periodic partner
          Call pack_particle(i, dir, xp_send, id_down(n_down), dat_down(:,n_down))
          Call remove_particle(i)
       Else
          i = i + 1
       End If
    End Do

    n_recv_down = 0
    Call Mpi_sendrecv(n_up,   1, MPI_INTEGER, up,   200, n_recv_down, 1, MPI_INTEGER, down, 200, MPI_COMM_WORLD, istat, ierr)
    Allocate ( id_rd(Max(n_recv_down,1)), dat_rd(7,Max(n_recv_down,1)) )
    Call Mpi_sendrecv(id_up,  n_up,   MPI_INTEGER, up, 201, id_rd,  n_recv_down,   MPI_INTEGER, down, 201, &
         MPI_COMM_WORLD, istat, ierr)
    Call Mpi_sendrecv(dat_up, 7*n_up, Mpi_real8,   up, 202, dat_rd, 7*n_recv_down, Mpi_real8,   down, 202, &
         MPI_COMM_WORLD, istat, ierr)

    n_recv_up = 0
    Call Mpi_sendrecv(n_down,   1, MPI_INTEGER, down, 203, n_recv_up, 1, MPI_INTEGER, up, 203, MPI_COMM_WORLD, istat, ierr)
    Allocate ( id_ru(Max(n_recv_up,1)), dat_ru(7,Max(n_recv_up,1)) )
    Call Mpi_sendrecv(id_down,  n_down,   MPI_INTEGER, down, 204, id_ru,  n_recv_up,   MPI_INTEGER, up, 204, &
         MPI_COMM_WORLD, istat, ierr)
    Call Mpi_sendrecv(dat_down, 7*n_down, Mpi_real8,   down, 205, dat_ru, 7*n_recv_up, Mpi_real8,   up, 205, &
         MPI_COMM_WORLD, istat, ierr)

    Do i = 1, n_recv_down
       Call unpack_particle(id_rd(i), dat_rd(:,i))
    End Do
    Do i = 1, n_recv_up
       Call unpack_particle(id_ru(i), dat_ru(:,i))
    End Do

    Deallocate ( id_up, id_down, dat_up, dat_down, id_rd, id_ru, dat_rd, dat_ru )

  End Subroutine exchange_particles_1d

  !> This rank's ordinary-interior-neighbour ownership bound on axis dir (1=x,2=z), from
  !  the cell-centre ghost arrays (xg/zg) -- see owns_particle for why xg/zg, not x/z.
  Real(Int64) Function xg_edge_lo(dir) Result(v)
    Integer(Int32), Intent(In) :: dir
    If ( dir == 1 ) Then;  v = xg(2);  Else;  v = zg(2);  End If
  End Function xg_edge_lo

  Real(Int64) Function xg_edge_hi(dir) Result(v)
    Integer(Int32), Intent(In) :: dir
    If ( dir == 1 ) Then;  v = xg(nxg);  Else;  v = zg(nzg);  End If
  End Function xg_edge_hi

  !> The true global domain edge on axis dir (1=x,2=z) -- only consulted on is_first/is_last.
  Real(Int64) Function x_edge_lo(dir) Result(v)
    Integer(Int32), Intent(In) :: dir
    If ( dir == 1 ) Then;  v = x_global(1);  Else;  v = z_global(1);  End If
  End Function x_edge_lo

  Real(Int64) Function x_edge_hi(dir) Result(v)
    Integer(Int32), Intent(In) :: dir
    If ( dir == 1 ) Then;  v = x_global(nx_global);  Else;  v = z_global(nz_global);  End If
  End Function x_edge_hi

  Subroutine pack_particle(i, dir, pos_send, id_out, dat_out)
    Integer(Int32), Intent(In)  :: i, dir
    Real   (Int64), Intent(In)  :: pos_send
    Integer(Int32), Intent(Out) :: id_out
    Real   (Int64), Intent(Out) :: dat_out(7)
    id_out = p_id(i)
    dat_out(1) = p_x(i);  dat_out(2) = p_y(i);  dat_out(3) = p_z(i)
    dat_out(4) = p_u(i);  dat_out(5) = p_v(i);  dat_out(6) = p_w(i)
    dat_out(7) = p_age(i)
    If ( dir == 1 ) Then
       dat_out(1) = pos_send
    Else
       dat_out(3) = pos_send
    End If
  End Subroutine pack_particle

  Subroutine unpack_particle(id_in, dat_in)
    Integer(Int32), Intent(In) :: id_in
    Real   (Int64), Intent(In) :: dat_in(7)
    Call ensure_capacity(n_particles_local + 1)
    n_particles_local = n_particles_local + 1
    p_id (n_particles_local) = id_in
    p_x  (n_particles_local) = dat_in(1)
    p_y  (n_particles_local) = dat_in(2)
    p_z  (n_particles_local) = dat_in(3)
    p_u  (n_particles_local) = dat_in(4)
    p_v  (n_particles_local) = dat_in(5)
    p_w  (n_particles_local) = dat_in(6)
    p_age(n_particles_local) = dat_in(7)
    ! SGS fluctuation state is not carried across migration (documented simplification --
    ! the well-mixed model's timescale T_L is normally short compared to a full-domain
    ! transit, so restarting it at 0 on migration has limited practical impact)
    If ( sgs_particle_model == 1 ) Then
       p_sgs_u(n_particles_local) = 0d0
       p_sgs_v(n_particles_local) = 0d0
       p_sgs_w(n_particles_local) = 0d0
    End If
  End Subroutine unpack_particle

  !> Replace particles that exited through the outflow face (n_exit_outflow_local, counted on
  !  the is_last-x-row rank that detected each exit) with new ones at the inflow plane, on the
  !  is_first-x-row rank in the SAME z-column (rank = Mod(myid,p_col) when row-major
  !  rank=row*p_col+col) -- a direct, communication-cheap pairing that needs no search, since
  !  x- and z-decomposition are independent. New particles are seeded at a small offset past
  !  the inflow face, with y/z drawn uniformly from this receiving rank's own owned z-range
  !  (not the full inflow plane) and the seed box's y-range, so no further ownership check or
  !  hand-off is needed after creation. A no-op when p_row==1 (is_first==is_last==this rank).
  Subroutine reinject_at_inflow

    Integer(Int32) :: target_rank, source_rank, n_new, i, new_id
    Real   (Int64) :: eps, yp, zp, ry, rz

    target_rank = Mod(myid, p_col)          ! row 0, same column -- only meaningful if this rank is on the outflow row
    source_rank = (p_row-1)*p_col + Mod(myid, p_col)   ! the outflow-row rank in this rank's own column

    If ( is_first_x_m .And. is_last_x_m ) Then
       ! single x-row: reinject locally using this rank's own exit count, no communication
       n_new = n_exit_outflow_local
    Else If ( is_first_x_m ) Then
       Call Mpi_recv(n_new, 1, MPI_INTEGER, source_rank, 300, MPI_COMM_WORLD, istat, ierr)
    Else If ( is_last_x_m ) Then
       Call Mpi_send(n_exit_outflow_local, 1, MPI_INTEGER, target_rank, 300, MPI_COMM_WORLD, ierr)
       Return
    Else
       Return
    End If

    eps = 1d-6 * (x_global(2) - x_global(1))
    Do i = 1, n_new
       Call random_number(ry);  Call random_number(rz)
       yp = particle_seed_ymin + ry*(particle_seed_ymax - particle_seed_ymin)
       zp = zg(2) + rz*(zg(nzg-1) - zg(2))   ! this rank's own owned z-range only, see header note
       n_reinject_counter = n_reinject_counter + 1
       new_id = myid*10000000 + n_reinject_counter
       Call ensure_capacity(n_particles_local + 1)
       n_particles_local = n_particles_local + 1
       p_id (n_particles_local) = new_id
       p_x  (n_particles_local) = x_global(1) + eps
       p_y  (n_particles_local) = yp
       p_z  (n_particles_local) = zp
       p_u  (n_particles_local) = 0d0
       p_v  (n_particles_local) = 0d0
       p_w  (n_particles_local) = 0d0
       p_age(n_particles_local) = 0d0
       If ( sgs_particle_model == 1 ) Then
          p_sgs_u(n_particles_local) = 0d0
          p_sgs_v(n_particles_local) = 0d0
          p_sgs_w(n_particles_local) = 0d0
       End If
       n_reinjected_local = n_reinjected_local + 1
    End Do

  End Subroutine reinject_at_inflow

  !> Reduce and print per-rank event counters (global sums), then reset the running totals.
  Subroutine report_particle_counts

    Integer(Int64) :: local_buf(4), global_buf(4)
    Integer(Int32) :: n_active_global

    local_buf = (/ n_exited_local, n_deposited_local, n_reinjected_local, Int(n_particles_local, Int64) /)
    Call MPI_Reduce(local_buf, global_buf, 4, MPI_INTEGER8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    If ( myid == 0 ) Then
       Write(*,'(A,I10,A,I10,A,I10,A,I12)') ' particles: exited=', global_buf(1), &
            ' deposited=', global_buf(2), ' reinjected=', global_buf(3), ' active=', global_buf(4)
    End If

    n_exited_local = 0;  n_deposited_local = 0;  n_reinjected_local = 0

  End Subroutine report_particle_counts

  !> Sum every rank's particle_deposit_x (each rank only ever populated indices in its own
  !  owned x-range, so the sum has no double-counting) and append one row to
  !  particle_deposit_file: istep, t, then nxg_global deposit counts since the last write.
  Subroutine write_deposit_profile

    Real(Int64), Allocatable :: global_x(:)
    Integer(Int32) :: funit, k
    Logical :: file_exists

    Allocate ( global_x(nxg_global) )
    Call MPI_Reduce(particle_deposit_x, global_x, nxg_global, Mpi_real8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    If ( myid == 0 ) Then
       Inquire(file=Trim(particle_deposit_file), exist=file_exists)
       Open(newunit=funit, file=Trim(particle_deposit_file), status='unknown', &
            position='append', action='write')
       If ( .Not. file_exists ) Then
          Write(funit,'(A)') 'istep,t,deposit_count_per_x_cell...'
       End If
       Write(funit,'(I10,A,ES14.6,999999(A,ES14.6))') istep, ',', t, (',', global_x(k), k=1,nxg_global)
       Close(funit)
    End If

    particle_deposit_x = 0d0
    Deallocate ( global_x )

  End Subroutine write_deposit_profile

  !> Rank-0 gather of every particle's (id) and (x,y,z,u,v,w,age) into flat arrays,
  !  shared by write_particle_restart and write_particle_snapshot. On return, id_all/
  !  dat_all are allocated and populated on EVERY rank the same way write_particle_restart
  !  always did (harmless on non-root ranks: Mpi_gatherv only reads counts/displs at the
  !  root) -- total is 0 (and the arrays size-1 placeholders) on non-root ranks.
  Subroutine gather_all_particles(id_all, dat_all, total)

    Integer(Int32), Allocatable, Intent(Out) :: id_all(:)
    Real   (Int64), Allocatable, Intent(Out) :: dat_all(:)
    Integer(Int32), Intent(Out) :: total

    Integer(Int32) :: counts(0:nprocs-1), displs(0:nprocs-1), rcounts(0:nprocs-1), rdispls(0:nprocs-1)
    Integer(Int32) :: r, i
    Real(Int64), Allocatable :: dat_local(:)

    Call Mpi_gather(n_particles_local, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    total = 0
    If ( myid == 0 ) Then
       displs(0) = 0
       Do r = 1, nprocs-1
          displs(r) = displs(r-1) + counts(r-1)
       End Do
       total = Sum(counts)
       rcounts = counts * 7
       rdispls = displs * 7
    End If
    Allocate ( id_all(Max(total,1)), dat_all(Max(total*7,1)) )

    Allocate ( dat_local(7*Max(n_particles_local,1)) )
    Do i = 1, n_particles_local
       dat_local(7*(i-1)+1) = p_x(i);   dat_local(7*(i-1)+2) = p_y(i);   dat_local(7*(i-1)+3) = p_z(i)
       dat_local(7*(i-1)+4) = p_u(i);   dat_local(7*(i-1)+5) = p_v(i);   dat_local(7*(i-1)+6) = p_w(i)
       dat_local(7*(i-1)+7) = p_age(i)
    End Do

    Call Mpi_gatherv(p_id,      n_particles_local,   MPI_INTEGER, id_all,  counts,  displs,  MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    Call Mpi_gatherv(dat_local, 7*n_particles_local, Mpi_real8,   dat_all, rcounts, rdispls, Mpi_real8,   0, MPI_COMM_WORLD, ierr)

    Deallocate ( dat_local )

  End Subroutine gather_all_particles

  !> Write a flat particle-restart file (rank-0 gather, per the repo's own rank-0-centric
  !  restart I/O convention, input_output.f90): a global count header, then that many
  !  (id, x,y,z,u,v,w,age) records. Overwrites particle_restart_file every call (a single
  !  latest-state snapshot for hot-restart, NOT a time series -- see write_particle_snapshot
  !  for the visualization time series). Called from output_data at the field-snapshot cadence.
  Subroutine write_particle_restart

    Integer(Int32) :: total, funit
    Integer(Int32), Allocatable :: id_all(:)
    Real   (Int64), Allocatable :: dat_all(:)

    If ( particles_active < 1 ) Return

    Call gather_all_particles(id_all, dat_all, total)

    If ( myid == 0 ) Then
       Open(newunit=funit, file=Trim(particle_restart_file), access='stream', form='unformatted', status='replace')
       Write(funit) total
       If ( total > 0 ) Then
          Write(funit) id_all(1:total)
          Write(funit) dat_all(1:7*total)
       End If
       Close(funit)
    End If

    Deallocate ( id_all, dat_all )

  End Subroutine write_particle_restart

  !> Write a timestamped particle snapshot for visualization (same (id,x,y,z,u,v,w,age)
  !  binary layout as write_particle_restart, but to fname -- normally 'fields/<fileout>
  !  _particles.<istep>', one file per saved step, called from output_data alongside the
  !  main field snapshot so postProcessing/generate_particles_xmf.py can build a ParaView
  !  time series that scrubs together with the field snapshots' own XDMF timeline).
  Subroutine write_particle_snapshot(fname)

    Character(*), Intent(In) :: fname

    Integer(Int32) :: total, funit
    Integer(Int32), Allocatable :: id_all(:)
    Real   (Int64), Allocatable :: dat_all(:)

    If ( particles_active < 1 ) Return

    Call gather_all_particles(id_all, dat_all, total)

    If ( myid == 0 ) Then
       Open(newunit=funit, file=Trim(fname), access='stream', form='unformatted', status='replace')
       Write(funit) total
       If ( total > 0 ) Then
          Write(funit) id_all(1:total)
          Write(funit) dat_all(1:7*total)
       End If
       Close(funit)
    End If

    Deallocate ( id_all, dat_all )

  End Subroutine write_particle_snapshot

  !> Read the particle-restart file (if present) and redistribute every record to whichever
  !  rank owns it now (owns_particle), independent of p_row/p_col matching the writing run.
  !  loaded=.False. (no file found) tells setup_particles to seed fresh particles instead.
  Subroutine read_particle_restart(loaded)

    Logical, Intent(Out) :: loaded

    Logical :: file_exists
    Integer(Int32) :: total, funit, i
    Integer(Int32), Allocatable :: id_all(:)
    Real   (Int64), Allocatable :: dat_all(:)
    Real   (Int64) :: xp, zp

    file_exists = .False.
    If ( myid == 0 ) Inquire(file=Trim(particle_restart_file), exist=file_exists)
    Call Mpi_bcast(file_exists, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

    loaded = .False.
    If ( .Not. file_exists ) Then
       If ( myid == 0 ) Write(*,'(A)') ' particles: no restart file found, seeding fresh particles instead'
       Return
    End If

    total = 0
    If ( myid == 0 ) Then
       Open(newunit=funit, file=Trim(particle_restart_file), access='stream', form='unformatted', status='old', action='read')
       Read(funit) total
    End If
    Call Mpi_bcast(total, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    Allocate ( id_all(Max(total,1)), dat_all(Max(total*7,1)) )
    If ( myid == 0 .And. total > 0 ) Then
       Read(funit) id_all(1:total)
       Read(funit) dat_all(1:7*total)
    End If
    If ( myid == 0 ) Close(funit)

    If ( total > 0 ) Then
       Call Mpi_bcast(id_all,  total,   MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
       Call Mpi_bcast(dat_all, 7*total, Mpi_real8,   0, MPI_COMM_WORLD, ierr)
    End If

    Do i = 1, total
       xp = dat_all(7*(i-1)+1);  zp = dat_all(7*(i-1)+3)
       If ( owns_particle(xp, zp) ) Then
          Call ensure_capacity(n_particles_local + 1)
          n_particles_local = n_particles_local + 1
          p_id (n_particles_local) = id_all(i)
          p_x  (n_particles_local) = dat_all(7*(i-1)+1)
          p_y  (n_particles_local) = dat_all(7*(i-1)+2)
          p_z  (n_particles_local) = dat_all(7*(i-1)+3)
          p_u  (n_particles_local) = dat_all(7*(i-1)+4)
          p_v  (n_particles_local) = dat_all(7*(i-1)+5)
          p_w  (n_particles_local) = dat_all(7*(i-1)+6)
          p_age(n_particles_local) = dat_all(7*(i-1)+7)
          If ( sgs_particle_model == 1 ) Then
             p_sgs_u(n_particles_local) = 0d0
             p_sgs_v(n_particles_local) = 0d0
             p_sgs_w(n_particles_local) = 0d0
          End If
       End If
    End Do

    Deallocate ( id_all, dat_all )
    loaded = .True.
    If ( myid == 0 ) Write(*,'(A,I8,A)') ' particles: loaded ', total, ' particles from restart file'

  End Subroutine read_particle_restart

End Module particles
