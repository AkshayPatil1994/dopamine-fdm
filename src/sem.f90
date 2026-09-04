!> Streamwise inflow BC (x_bc_type==1): constant uniform flow, Ensemble SEM, or recycled precursor slice
Module synthetic_eddy_method

  ! Modules
  Use iso_fortran_env, Only : Int32, Int64
  Use global
  Use mpi

  ! prevent implicit typing
  Implicit None

  ! ---- mean/Reynolds-stress reference profile (inflow_type==1) ----------
  Integer(Int32) :: n_profile = 0
  Real   (Int64), Allocatable :: prof_y(:), prof_U(:)
  Real   (Int64), Allocatable :: prof_R11(:), prof_R22(:), prof_R33(:), prof_R12(:)
  Logical :: inflow_profile_loaded = .False.   ! guards init_inflow_profile against a second (redundant) file read

  ! optional SEM inflow mean-temperature profile (inflow_temperature_file); n_profile_T==0 -> mean_profile_T falls back to T_ref
  Integer(Int32) :: n_profile_T = 0
  Real   (Int64), Allocatable :: prof_y_T(:), prof_T(:)
  !$acc declare create(n_profile_T, prof_y_T, prof_T)

  ! inhomogeneous eddy length-scale profile sigma_ij(y) (i=u,v,w; j=x,y,z); n_sigma==0 (sem_sigma_file unset) => homogeneous mode
  Integer(Int32) :: n_sigma = 0
  Real   (Int64), Allocatable :: sig_y(:)
  Real   (Int64), Allocatable :: sig_ux(:), sig_uy(:), sig_uz(:)
  Real   (Int64), Allocatable :: sig_vx(:), sig_vy(:), sig_vz(:)
  Real   (Int64), Allocatable :: sig_wx(:), sig_wy(:), sig_wz(:)

  ! placement CDF table (Section 4.2), built once if sem_eddy_placement==1
  Integer(Int32) :: n_cdf = 0
  Real   (Int64), Allocatable :: cdf_y(:), cdf_val(:)

  ! eddy population (inflow_type==1): fixed size class/x-phase, but (y,z) and sign re-drawn every recycle -- eddy-recycling/phase-locking rationale
  Real   (Int64), Allocatable :: eddy_x0(:)
  Real   (Int64), Allocatable :: eddy_y0(:)        ! (sem_n_eddies): placement height this eddy's size class (eddy_sig) was calibrated at; anchors eddy_realize's per-cycle y_r redraw (inhomogeneous mode) so a wall-sized eddy can't wander into the outer layer and vice versa
  Real   (Int64), Allocatable :: eddy_sig(:,:,:)   ! (3,3,sem_n_eddies): (i,j,k)=sigma_ij of eddy k
  Real   (Int64), Allocatable :: eddy_smax(:,:)    ! (3,sem_n_eddies): (j,k)=max_i(sigma_ij) of eddy k
  Real   (Int64), Allocatable :: eddy_Tperiod(:)   ! (sem_n_eddies): 2*eddy_smax(1,k)/Uconv_sem
  Real   (Int64) :: Uconv_sem = 0d0                ! eddy-box convection velocity

  Real   (Int64), Parameter :: sem_placement_band_factor = 3d0  ! half-width (in units of eddy_smax(2,k)) of the y-band eddy_realize confines eddy k's re-realized position to, inhomogeneous mode only

  ! placement box bounds, needed at runtime (not just at init) by
  ! eddy_realize's per-cycle re-randomization
  Real   (Int64) :: box_y_lo, box_y_hi, box_z_lo, box_z_hi

  ! ESEM empirical normalisation (Schau et al. 2022 sec. 3): per-grid-point mean/std of the raw eddy-kernel sum over an ensemble window, computed once at start-up; each staggered component (U,V,W) needs its own table(s)
  Real   (Int64), Allocatable :: ens_mean_uU(:,:), ens_std_uU(:,:)   ! U's grid (nyg,nzg): u* stats
  Real   (Int64), Allocatable :: ens_mean_uV(:,:), ens_std_uV(:,:)   ! V's grid (ny,nzg):  u* stats
  Real   (Int64), Allocatable :: ens_mean_vV(:,:), ens_std_vV(:,:)   ! V's grid (ny,nzg):  v* stats
  Real   (Int64), Allocatable :: ens_mean_wW(:,:), ens_std_wW(:,:)   ! W's grid (nyg,nz):  w* stats

  ! classical (oSEM, Jarrin 2006) analytical normalisation, used only when sem_use_esem==0, for comparison against ESEM (homogeneous mode only, see init_inflow)
  Real   (Int64) :: sem_norm = 1d0   ! sqrt(V_B); tent_kernel already carries sigma^-3/2

  Real   (Int64), Parameter :: sem_clip_sigma = 5d0  ! safety-net clip, local target std-devs

  ! near-wall no-slip enforcement for the fluctuating field (Section sem_fluctuation): precomputed once (host-side, init_inflow) from y_global/y_bc_type/bc_face_ylo/yhi so sem_fluctuation's !$acc routine seq body only needs these device-resident scalars, not y_global itself
  Integer(Int32) :: wall_active_lo = 0, wall_active_hi = 0
  Real   (Int64) :: wall_y_lo = 0d0, wall_y_hi = 0d0, wall_Ltaper_lo = 1d0, wall_Ltaper_hi = 1d0

  ! device residency for the state sem_fluctuation's per-step (!$acc routine seq) call chain reads; populated once by init_inflow (see its !$acc update device calls), read-only thereafter
  !$acc declare create(n_profile, prof_y, prof_U, prof_R11, prof_R22, prof_R33, prof_R12)
  !$acc declare create(n_sigma, sig_y, sig_ux, sig_uy, sig_uz, sig_vx, sig_vy, sig_vz, sig_wx, sig_wy, sig_wz)
  !$acc declare create(n_cdf, cdf_y, cdf_val)
  !$acc declare create(eddy_x0, eddy_y0, eddy_sig, eddy_smax, eddy_Tperiod, Uconv_sem)
  !$acc declare create(box_y_lo, box_y_hi, box_z_lo, box_z_hi)
  !$acc declare create(ens_mean_uU, ens_std_uU, ens_mean_uV, ens_std_uV, ens_mean_vV, ens_std_vV, ens_mean_wW, ens_std_wW)
  !$acc declare create(sem_norm)
  !$acc declare create(wall_active_lo, wall_active_hi, wall_y_lo, wall_y_hi, wall_Ltaper_lo, wall_Ltaper_hi)

  ! in-situ TI-profile rescaling (ti_rescale_active==1): host-only bookkeeping
  ! prof_R11_target/R22_target/R33_target: immutable copy of the originally-read profile; prof_R11/R22/R33 are the mutable, currently-injected profile that gets nudged
  Real   (Int64), Allocatable :: prof_R11_target(:), prof_R22_target(:), prof_R33_target(:)
  Integer(Int32) :: ti_i0 = 0                ! resolved global cc x-index of the sampling station
  Real   (Int64), Allocatable :: ti_yg_cc(:)   ! cell-centre y-coords (nym_global), for interpolation
  ! raw-moment accumulators, one entry per y cell-centre, summed over local z and time
  Real   (Int64), Allocatable :: acc_ti_U(:), acc_ti_V(:), acc_ti_W(:)
  Real   (Int64), Allocatable :: acc_ti_UU(:), acc_ti_VV(:), acc_ti_WW(:), acc_ti_n(:)
  ! EMA-filtered measured variance (on prof_y) and window counter, carried across apply_ti_rescale calls
  Real   (Int64), Allocatable :: ti_R11_filt(:), ti_R22_filt(:), ti_R33_filt(:)
  ! gain-decay is gated per-y-point on that point's own convergence, not a single profile-wide flag: a point whose target is physically unreachable (e.g. permanently anti-windup-clamped) must not hold the rest of the profile at undecayed gain forever
  Integer(Int32), Allocatable :: ti_R_decay_k(:)   ! per-m windows since the variance loop's gain decay started at that y
  Logical,        Allocatable :: ti_R_decaying(:)  ! per-m: has the variance loop's gain decay started at that y
  Integer(Int32), Allocatable :: ti_U_decay_k(:)   ! per-m windows since the mean-U loop's gain decay started at that y
  Logical,        Allocatable :: ti_U_decaying(:)  ! per-m: has the mean-U loop's gain decay started at that y
  ! true once the EMA filters have been seeded with a real measurement; kept separate from the per-m decay counters so resetting one of them doesn't re-trigger the cold-start (unfiltered) branch and discard the accumulated filter history
  Logical :: ti_filt_seeded = .False.
  ! mean-profile (U) companion state: immutable target snapshot and EMA-filtered measured mean, mirrors prof_R11_target/ti_R11_filt above
  Real   (Int64), Allocatable :: prof_U_target(:)
  Real   (Int64), Allocatable :: ti_U_filt(:)

  ! persisted controller state for restart (see write/read_ti_rescale_restart)
  Character(*), Parameter :: ti_rescale_restart_file = 'fields/ti_rescale_data.dat'
  ! leading sentinel + version tag, so read_ti_rescale_restart can tell a current-format file (which starts with this magic) from a pre-ti_gain_decay_started file (which starts directly with n_profile, a small positive count that can never equal this magic)
  Integer(Int32), Parameter :: ti_rescale_restart_magic = -987654321
  Integer(Int32), Parameter :: ti_rescale_restart_version = 4

  ! ---- recycled precursor inflow (inflow_type==2) ------------------------
  ! rec_active is only .True. on the rank owning the x=1 face (row==0 in the
  ! p_row/p_col grid); every other rank leaves it .False. and every routine
  ! below is a no-op there. Can't tell this apart from rec_unit's sign alone:
  ! Open(newunit=...) always hands back a negative unit (the standard's
  ! reserved range), so "rec_unit<0" is true both before AND after a
  ! successful open -- that previously made update_inflow_recycle return
  ! before ever reading a donor frame, on every rank, every call.
  Logical        :: rec_active  = .False.
  Integer(Int32) :: rec_unit    = -1
  Integer(Int32) :: rec_ncomp   = 0
  Integer(Int32) :: rec_n1      = 0   ! U/W/T/C's native n1 = nym_global (cell-centre y)
  Integer(Int32) :: rec_n1_v    = 0   ! V's own native n1 = ny_global (y-face, no half-cell shift needed) when rec_v_native; else == rec_n1 (legacy shared-grid donor, see rec_v_native)
  Logical        :: rec_v_native = .True.  ! .True.: donor meta had 'n1_V' (dopamine-ESEM / per-component-native-grid writer) -- V read at its own y-face index, no shift. .False.: older/generic probe_output.f90 x-normal slice donor, which cell-centre-interpolates V onto the same (nym_global) grid as U/W -- V must then be indexed with the same j-1 cell-centre shift as U/W (see recycle_value)
  Integer(Int32) :: rec_n2_global = 0 ! all components' donor n2 = nzm_global (U/V cell-centre z; W left-face-of-cell z)
  Integer(Int32) :: rec_col_U   = 0, rec_col_V = 0, rec_col_W = 0   ! presence flags (1 if the component is in the donor's comps, 0 if absent) -- U,V,W always required
  Integer(Int32) :: rec_col_T   = 0   ! 0 if the donor slice didn't record T (fine unless boussinesq_flag>=1, checked in init_inflow_recycle)
  Integer(Int32) :: rec_col_C   = 0   ! 0 if the donor slice didn't record C (fine unless sediment_flag>=1, checked in init_inflow_recycle)
  Integer(Int32) :: rec_nsnaps  = 0
  Real   (Int64), Allocatable :: rec_times(:)   ! full donor snapshot-time array; tiny (one Real64 per snapshot), kept resident on the host only
  Integer(Int32) :: rec_k1 = 0, rec_k2 = 0       ! this rank's local donor z range (global 1-based interior cc/face indices, same convention as kg1_global/kg2_global-2)
  Integer(Int64) :: rec_frame_bytes = 0_Int64    ! bytes per full per-snapshot frame (sum of the selected components' own (n1_c,n2) blocks, fixed U,V,W,T,C order)
  Integer(Int64) :: rec_off_U = 0_Int64, rec_off_V = 0_Int64, rec_off_W = 0_Int64   ! each component's byte offset of its block's start within a frame
  Integer(Int64) :: rec_off_T = 0_Int64, rec_off_C = 0_Int64
  Integer(Int32) :: rec_idx_lo = -1, rec_idx_hi = -1   ! cached bracketing donor frame indices, so a step that doesn't cross a donor sample re-reads nothing
  Integer(Int32) :: rec_shift_z   = 0                   ! spanwise shift (donor z-cells) for the current pass through the donor timeline; re-drawn every donor-loop wrap (inflow_recycle_shift_z==1)
  Integer(Int32) :: rec_loop_idx  = -999999              ! which pass through the donor timeline rec_shift_z was drawn for
  Integer(Int32) :: rec_idx_shift = -999999              ! shift value currently baked into rec_lo_*/rec_hi_* (vs. rec_shift_z, the shift that should be active now)
  ! bracketing frames per component -- each on its OWN native (n1_c, local nz) grid (no
  ! half-cell shift needed on read: U/W index directly via rec_n1, V via rec_n1_v)
  Real   (Int64), Allocatable :: rec_lo_U(:,:), rec_hi_U(:,:)
  Real   (Int64), Allocatable :: rec_lo_V(:,:), rec_hi_V(:,:)
  Real   (Int64), Allocatable :: rec_lo_W(:,:), rec_hi_W(:,:)
  Real   (Int64), Allocatable :: rec_lo_T(:,:), rec_hi_T(:,:)
  Real   (Int64), Allocatable :: rec_lo_C(:,:), rec_hi_C(:,:)
  Real   (Int64) :: rec_frac = 0d0   ! time-interpolation weight between rec_lo_* and rec_hi_*

  !$acc declare create(rec_lo_U, rec_hi_U, rec_lo_V, rec_hi_V, rec_lo_W, rec_hi_W, rec_lo_T, rec_hi_T, rec_lo_C, rec_hi_C, rec_frac, rec_col_U, rec_col_V, rec_col_W, rec_col_T, rec_col_C, rec_n1, rec_n1_v, rec_v_native)

Contains

  !> Minimum grid spacing over the whole domain (x,y,z faces); resolvability floor for the SEM/TI-rescale auto-tuning heuristics
  Real(Int64) Function min_grid_spacing() Result(dmin)

    dmin = Min( Minval(x_global(2:nx_global)-x_global(1:nx_global-1)), &
                Minval(y_global(2:ny_global)-y_global(1:ny_global-1)), &
                Minval(z_global(2:nz_global)-z_global(1:nz_global-1)) )

  End Function min_grid_spacing

  !> Per-direction maximum grid spacing (x,y,z), the building block behind both max_grid_spacing()
  !> (isotropic sem_length_scale's floor) and enforce_sigma_floor (sem_sigma_file's per-direction
  !> sigma_ij floor, which matches each component to the spacing in ITS OWN direction rather than one
  !> isotropic value -- a wall-clustered y-grid's coarse mid-channel spacing shouldn't be diluted by a
  !> fine spanwise z-grid, or vice versa).
  Subroutine max_grid_spacing_xyz(dx_max, dy_max, dz_max)

    Real(Int64), Intent(Out) :: dx_max, dy_max, dz_max

    dx_max = Maxval(x_global(2:nx_global)-x_global(1:nx_global-1))
    dy_max = Maxval(y_global(2:ny_global)-y_global(1:ny_global-1))
    dz_max = Maxval(z_global(2:nz_global)-z_global(1:nz_global-1))

  End Subroutine max_grid_spacing_xyz

  !> Maximum grid spacing over the whole domain (x,y,z faces): the resolvability floor for
  !> sem_length_scale needs THIS, not min_grid_spacing() -- an eddy's footprint is ~2*sem_length_scale
  !> in every direction (homogeneous mode: eddy_sig is isotropic, see place_eddies), so what matters is
  !> whether the COARSEST relevant spacing (typically the streamwise spacing, or a wall-clustered grid's
  !> mid-channel y-spacing) can resolve it; min_grid_spacing() is dominated by the tiny near-wall cell,
  !> which an eddy doesn't need to fit inside and so is far too permissive a floor on its own.
  Real(Int64) Function max_grid_spacing() Result(dmax)

    Real(Int64) :: dx_max, dy_max, dz_max

    Call max_grid_spacing_xyz(dx_max, dy_max, dz_max)
    dmax = Max(dx_max, dy_max, dz_max)

  End Function max_grid_spacing

  !> Mean streamwise (x) grid spacing; used (paired with Uconv_sem) to estimate dt for the ti_rescale_freq/nstart auto-tuning -- unlike min_grid_spacing(), this deliberately excludes the wall-clustered y-spacing, since a near-wall cell's tiny spacing pairs with a near-zero local velocity there and is not representative of the convective step size that limits dt
  Real(Int64) Function streamwise_grid_spacing() Result(dx)

    dx = ( x_global(nx_global) - x_global(1) ) / Real(nx_global-1,8)

  End Function streamwise_grid_spacing

  !> Friction-velocity estimate from the imposed pressure gradient, or (dPdx~0, e.g. open-channel/free-stream) from the wall-normal derivative of the mean inflow profile; shared by init_ti_rescale's nstart auto-tuning and place_eddies' wall damping
  Real(Int64) Function estimate_u_tau() Result(u_tau)

    Real(Int64) :: dUdy_wall

    If ( Abs(dPdx) > 1d-12 ) Then
       u_tau = Sqrt( Abs(dPdx)*Ly/2d0 )
    Else
       dUdy_wall = ( prof_U(2)-prof_U(1) ) / ( prof_y(2)-prof_y(1) )
       u_tau = Sqrt( nu*Abs(dUdy_wall) )
    End If

  End Function estimate_u_tau

  !> Provisional sem_length_scale from domain/grid geometry alone, used while parsing the inflow profile when sem_length_scale<=0 requests auto-tuning; refined by refine_sem_length_scale once the profile is available
  Real(Int64) Function geometric_length_scale_estimate() Result(L0)

    L0 = 0.3d0 * Min(Ly, Lz)
    L0 = Max( L0, 5d0*min_grid_spacing() )

  End Function geometric_length_scale_estimate

  !> Refines the provisional sem_length_scale using the now-parsed inflow profile: average of the measured/derived streamwise length-scale profile (sig_ux) if available, else a direct mixing-length estimate sqrt(R11)/|dU/dy| from prof_R11/prof_U (mirrors derive_mixing_length); clamped to [5*min grid spacing, 0.3*min(Ly,Lz)]
  Subroutine refine_sem_length_scale

    Real(Int64) :: Lfloor, Lcap, dUdy, Lsum
    Integer(Int32) :: i

    Lfloor = 5d0*min_grid_spacing()
    Lcap   = 0.3d0*Min(Ly,Lz)

    If ( n_sigma > 0 ) Then
       sem_length_scale = Sum(sig_ux) / Real(n_sigma,8)
    Else
       Lsum = 0d0
       Do i = 1, n_profile
          If ( i == 1 ) Then
             dUdy = ( prof_U(2) - prof_U(1) ) / ( prof_y(2) - prof_y(1) )
          Else If ( i == n_profile ) Then
             dUdy = ( prof_U(n_profile) - prof_U(n_profile-1) ) / ( prof_y(n_profile) - prof_y(n_profile-1) )
          Else
             dUdy = ( prof_U(i+1) - prof_U(i-1) ) / ( prof_y(i+1) - prof_y(i-1) )
          End If
          Lsum = Lsum + Sqrt(Max(prof_R11(i),0d0)) / Max(Abs(dUdy),1d-12)
       End Do
       sem_length_scale = Lsum / Real(n_profile,8)
    End If

    sem_length_scale = Min( Max(sem_length_scale, Lfloor), Lcap )

    If ( myid == 0 ) Write(*,'(A,E12.4)') &
         ' INFO: auto-tuned sem_length_scale (sem_length_scale<=0 requested) = ', sem_length_scale

  End Subroutine refine_sem_length_scale

  !> Unconditional resolvability floor on sem_length_scale, applied regardless of whether it was
  !> auto-tuned or set explicitly in the namelist: an eddy's footprint is ~2*sem_length_scale in every
  !> direction (homogeneous mode: eddy_sig is isotropic, see place_eddies), so a value narrower than a
  !> few grid cells cannot be represented by the mesh and is destroyed by advection/pressure-projection
  !> within the first cell or two of the domain -- observed in practice as injected Reynolds stresses
  !> collapsing right at the inflow even though the target statistics (verified offline, e.g. via
  !> postProcessing/check_inflow_donor.py) were exact. The auto-tune path (geometric_length_scale_estimate/
  !> refine_sem_length_scale) already has its own internal floor, but that one is keyed off
  !> min_grid_spacing() -- dominated by the smallest cell anywhere in the domain, almost always a
  !> wall-clustered near-wall cell an eddy never needs to fit inside, so it under-floors whenever the
  !> coarsest *relevant* direction (typically streamwise) is far coarser than the finest one. This uses
  !> max_grid_spacing() instead, which tracks the actual bottleneck. A user-supplied value had no floor
  !> at all before this -- only the auto-tune path was ever checked.
  Subroutine enforce_sem_length_scale_floor

    Real(Int64) :: Lfloor

    Lfloor = 5d0*max_grid_spacing()
    If ( sem_length_scale < Lfloor ) Then
       If ( myid == 0 ) Write(*,'(A,E12.4,A,E12.4,A)') &
            ' WARNING: sem_length_scale (', sem_length_scale, ') is below the grid''s resolvability floor ' // &
            '(5*max_grid_spacing = ', Lfloor, ') -- eddies this small cannot be represented by the mesh ' // &
            'and would be destroyed by advection/pressure-projection within the first cell or two of ' // &
            'the domain; raising it to the floor.'
       sem_length_scale = Lfloor
    End If

  End Subroutine enforce_sem_length_scale_floor

  !> Per-direction resolvability floor on the inhomogeneous eddy length-scale profile (sem_sigma_file,
  !> or sem_profile_format=1's derived sig_* -- see read_sigma_profile/derive_mixing_length/
  !> fill_length_scale), mirroring enforce_sem_length_scale_floor above for the homogeneous case: no-op
  !> when n_sigma==0 (homogeneous mode, already covered by enforce_sem_length_scale_floor). sigma_ij sets
  !> eddy component i's footprint in direction j (see place_eddies/unified_kernel), so unlike the isotropic
  !> sem_length_scale floor, each of the 9 profiles is floored against the spacing in ITS OWN direction --
  !> sig_*x vs the streamwise spacing, sig_*y vs the (possibly wall-clustered) y-spacing, sig_*z vs the
  !> spanwise spacing -- rather than one shared value, so a fine z-grid can't be diluted by a coarse y-grid
  !> or vice versa. Clamps in place; warns once per profile (not once per point) when at least one node
  !> needed raising, reporting that profile's pre-clamp minimum against the floor it was raised to.
  Subroutine enforce_sigma_floor

    Real(Int64) :: dx_max, dy_max, dz_max, Lfloor_x, Lfloor_y, Lfloor_z

    If ( n_sigma == 0 ) Return

    Call max_grid_spacing_xyz(dx_max, dy_max, dz_max)
    Lfloor_x = 5d0*dx_max;  Lfloor_y = 5d0*dy_max;  Lfloor_z = 5d0*dz_max

    Call floor_sigma_component('sig_ux', sig_ux, Lfloor_x)
    Call floor_sigma_component('sig_uy', sig_uy, Lfloor_y)
    Call floor_sigma_component('sig_uz', sig_uz, Lfloor_z)
    Call floor_sigma_component('sig_vx', sig_vx, Lfloor_x)
    Call floor_sigma_component('sig_vy', sig_vy, Lfloor_y)
    Call floor_sigma_component('sig_vz', sig_vz, Lfloor_z)
    Call floor_sigma_component('sig_wx', sig_wx, Lfloor_x)
    Call floor_sigma_component('sig_wy', sig_wy, Lfloor_y)
    Call floor_sigma_component('sig_wz', sig_wz, Lfloor_z)

  Contains

    Subroutine floor_sigma_component(name, sig, Lfloor)

      Character(*),   Intent(In)    :: name
      Real   (Int64), Intent(InOut) :: sig(:)
      Real   (Int64), Intent(In)    :: Lfloor

      Real(Int64) :: Lmin_before

      Lmin_before = Minval(sig)
      If ( Lmin_before < Lfloor ) Then
         If ( myid == 0 ) Write(*,'(A,A,A,E12.4,A,E12.4,A)') &
              ' WARNING: ', name, ' has a minimum (', Lmin_before, ') below the grid''s resolvability ' // &
              'floor in its direction (', Lfloor, ') -- eddies this small cannot be represented by the ' // &
              'mesh; clamping up.'
         sig = Max(sig, Lfloor)
      End If

    End Subroutine floor_sigma_component

  End Subroutine enforce_sigma_floor

  !> Auto-sets sem_n_eddies (sem_n_eddies<=0 requests this) from Jarrin's box-fill density guideline: ~1 eddy per sem_length_scale^3 cell of the padded inflow box; mirrors build_analytical_normalisation's recommended_N diagnostic
  Subroutine auto_set_n_eddies

    Real(Int64) :: y_lo, y_hi, z_lo, z_hi, VB, recommended_N
    Integer(Int32), Parameter :: n_eddies_cap = 200000

    y_lo = y_global(1)         - sem_length_scale
    y_hi = y_global(ny_global) + sem_length_scale
    z_lo = z_global(1)         - sem_length_scale
    z_hi = z_global(nz_global) + sem_length_scale

    VB            = (2d0*sem_length_scale) * (y_hi-y_lo) * (z_hi-z_lo)
    recommended_N = VB / sem_length_scale**3

    sem_n_eddies = Max( 1, Ceiling(recommended_N) )
    If ( sem_n_eddies > n_eddies_cap ) Then
       If ( myid == 0 ) Write(*,'(A,I8,A)') &
            ' WARNING: auto-tuned sem_n_eddies capped at ', n_eddies_cap, ' -- increase sem_length_scale to relax this'
       sem_n_eddies = n_eddies_cap
    End If

    If ( myid == 0 ) Write(*,'(A,I8)') &
         ' INFO: auto-tuned sem_n_eddies (sem_n_eddies<=0 requested) = ', sem_n_eddies

  End Subroutine auto_set_n_eddies

  ! One-time inflow setup. No-op unless x_bc_type==1.
  !> Read just the mean-profile file (prof_y/prof_U/prof_R**) so an IC generator can seed the interior from the same mean profile the inflow BC will impose; safe to call before the rest of init_inflow's eddy-population setup, and idempotent (init_inflow skips re-reading via inflow_profile_loaded)
  Subroutine init_inflow_profile

    If ( x_bc_type /= 1 ) Return
    If ( inflow_profile_loaded ) Return

    If ( inflow_type == 1 ) Then
       If ( sem_profile_format == 0 ) Then
          Call read_mean_profile
       Else If ( sem_profile_format == 1 ) Then
          Call read_mean_profile_TI
       Else
          Stop 'ERROR: sem_profile_format must be 0 (Reynolds-stress) or 1 (wind-tunnel TI)'
       End If
       inflow_profile_loaded = .True.
    Else If ( inflow_type == 2 ) Then
       ! mean profile only (z-average of the donor's first frame), so genGridandIC can
       ! seed the interior mean without a step-1 inflow/IC divergence spike, same as inflow_type==1
       Call read_recycle_mean_profile
       inflow_profile_loaded = .True.
    End If

  End Subroutine init_inflow_profile

  Subroutine init_inflow

    Logical :: auto_length_scale

    If ( x_bc_type /= 1 ) Return

    If ( inflow_type == 2 ) Then
       Call init_inflow_profile
       Call init_inflow_recycle
       Return
    End If

    If ( inflow_type /= 1 ) Return

    auto_length_scale = ( sem_length_scale <= 0d0 )
    If ( auto_length_scale ) sem_length_scale = geometric_length_scale_estimate()

    Call init_inflow_profile
    Call read_sigma_profile   ! sem_sigma_file, if set, overrides any length scales read_mean_profile_TI derived from inflow_profile_file
    Call enforce_sigma_floor  ! unconditional: catches a sig_* profile (explicit file or mixing-length-derived) too small for the grid to resolve, same rationale as enforce_sem_length_scale_floor below

    If ( boussinesq_flag >= 1 .And. Len_trim(inflow_temperature_file) > 0 ) Call read_mean_profile_T

    If ( sem_use_esem == 0 .And. ( n_sigma > 0 .Or. sem_wall_damping == 1 ) ) Stop &
         'ERROR: sem_use_esem=0 (classical SEM) requires sem_sigma_file to be unset and sem_wall_damping=0 -- ' // &
         'classical SEM has no well-defined normalisation for inhomogeneous length scales'

    If ( sem_divergence_free == 1 .And. sem_use_esem == 0 ) Stop &
         'ERROR: sem_divergence_free=1 requires sem_use_esem=1 -- classical (analytical) ' // &
         'normalisation has no closed form for the curl construction'

    Uconv_sem = Sum(prof_U) / Real(n_profile,8)
    If ( Uconv_sem <= 0d0 ) Stop 'ERROR: SEM inflow profile has non-positive bulk velocity'

    ! near-wall no-slip enforcement setup (see wall_taper/sem_fluctuation): only active at faces that are literal Dirichlet (no-slip) walls, not periodic or free-slip
    wall_active_lo = Merge( 1, 0, y_bc_type == 1 .And. bc_face_ylo == 1 )
    wall_active_hi = Merge( 1, 0, y_bc_type == 1 .And. bc_face_yhi == 1 )
    wall_y_lo      = y_global(1)
    wall_y_hi      = y_global(ny_global)
    wall_Ltaper_lo = 2d0*( y_global(2)         - y_global(1) )
    wall_Ltaper_hi = 2d0*( y_global(ny_global) - y_global(ny_global-1) )

    If ( auto_length_scale ) Call refine_sem_length_scale
    Call enforce_sem_length_scale_floor   ! unconditional: also catches an explicitly-set sem_length_scale too small for the grid to resolve
    If ( sem_n_eddies <= 0 ) Call auto_set_n_eddies

    Call place_eddies

    If ( sem_use_esem == 1 ) Then
       Call build_ensemble_normalisation
    Else
       Call build_analytical_normalisation
    End If

    ! one-time push of everything sem_fluctuation's per-step !$acc routine seq chain reads
    !$acc update device(n_profile, prof_y, prof_U, prof_R11, prof_R22, prof_R33, prof_R12)
    !$acc update device(wall_active_lo, wall_active_hi, wall_y_lo, wall_y_hi, wall_Ltaper_lo, wall_Ltaper_hi)
    If ( n_sigma > 0 ) Then
       !$acc update device(n_sigma, sig_y, sig_ux, sig_uy, sig_uz, sig_vx, sig_vy, sig_vz, sig_wx, sig_wy, sig_wz)
    Else
       !$acc update device(n_sigma)
    End If
    !$acc update device(Uconv_sem, box_y_lo, box_y_hi, box_z_lo, box_z_hi)
    !$acc update device(eddy_x0, eddy_y0, eddy_sig, eddy_smax, eddy_Tperiod)
    If ( sem_eddy_placement == 1 .And. n_sigma > 0 ) Then
       !$acc update device(n_cdf, cdf_y, cdf_val)
    End If
    !$acc update device(sem_norm)
    If ( sem_use_esem == 1 ) Then
       !$acc update device(ens_mean_uU, ens_std_uU, ens_mean_uV, ens_std_uV, ens_mean_vV, ens_std_vV, ens_mean_wW, ens_std_wW)
    End If
    If ( n_profile_T > 0 ) Then
       !$acc update device(n_profile_T, prof_y_T, prof_T)
    Else
       !$acc update device(n_profile_T)
    End If

  End Subroutine init_inflow

  !> Parse a recycled-inflow donor's <inflow_recycle_file>_meta.txt. Two donor formats are
  !> accepted: (1) the per-component-native-grid format written by dopamine-ESEM, which
  !> carries an explicit 'n1_V' key (V's own native y-face count, ny_global) alongside the
  !> shared U/W/T/C native n1 (cell-centre y, nym_global); (2) the older/generic format
  !> written by probe_output.f90's x-normal slice output (e.g. examples/precursor_successor),
  !> which cell-centre-interpolates V onto the same nym_global grid as U/W and has no 'n1_V'
  !> key -- for that format n1_v is defaulted to n1 and v_native returned .False., so the
  !> caller (init_inflow_recycle) and recycle_value read/index V using the same cell-centre
  !> convention as U/W instead of V's true (unavailable) face grid.
  Subroutine read_recycle_meta(ncomp, n1, n1_v, n2, nsnaps, colU, colV, colW, colT, colC, v_native)

    Integer(Int32), Intent(Out) :: ncomp, n1, n1_v, n2, nsnaps, colU, colV, colW, colT, colC
    Logical,        Intent(Out) :: v_native

    Integer(Int32) :: u_meta, ios, eqpos, ic
    Character(300) :: fname, line
    Character(20)  :: key
    Character(8)   :: comps_str

    If ( Len_trim(inflow_recycle_file) == 0 ) &
         Stop 'ERROR: inflow_type=2 (recycled inflow) requires inflow_recycle_file to be set'

    Write(fname,'(A,A)') Trim(inflow_recycle_file), '_meta.txt'
    Open(newunit=u_meta, file=Trim(fname), form='formatted', status='old', action='read', iostat=ios)
    If ( ios /= 0 ) Stop 'ERROR: cannot open inflow_recycle_file meta (expected <inflow_recycle_file>_meta.txt)'

    ncomp = 0; n1 = 0; n1_v = 0; n2 = 0; nsnaps = 0; comps_str = ''
    v_native = .False.
    Do
       Read(u_meta,'(A)', iostat=ios) line
       If ( ios /= 0 ) Exit
       eqpos = Index(line,'=')
       If ( eqpos < 1 ) Cycle
       key = Adjustl(line(1:eqpos-1))
       Select Case ( Trim(key) )
       Case ('ncomp');  Read(line(eqpos+1:),*) ncomp
       Case ('n1');     Read(line(eqpos+1:),*) n1
       Case ('n1_V');   Read(line(eqpos+1:),*) n1_v;  v_native = .True.
       Case ('n2');     Read(line(eqpos+1:),*) n2
       Case ('comps');  comps_str = Adjustl(line(eqpos+1:))
       Case ('nsnaps'); Read(line(eqpos+1:),*) nsnaps
       Case ('dir')
          line = Adjustl(line(eqpos+1:))
          If ( line(1:1) /= 'x' .And. line(1:1) /= 'X' ) &
               Stop 'ERROR: inflow_recycle_file meta is not an x-normal slice (dir/=x)'
       End Select
    End Do
    Close(u_meta)

    If ( ncomp < 1 .Or. n1 < 1 .Or. n2 < 1 .Or. nsnaps < 1 ) &
         Stop 'ERROR: inflow_recycle_file meta is incomplete or malformed (expected ncomp/n1/n2/comps/nsnaps, ' // &
         'as written by probe_output.f90''s slice output or dopamine-ESEM)'

    If ( .Not. v_native ) Then
       n1_v = n1   ! legacy/generic donor: V shares U/W's cell-centre grid on disk
       If ( myid == 0 ) Write(*,'(A)') ' NOTE: inflow_recycle_file has no n1_V -- treating as a generic ' // &
            '(non-ESEM) donor: V is read cell-centre-interpolated like U/W, not on its own native face grid'
    End If

    Do ic = 1, Len_trim(comps_str)
       If ( comps_str(ic:ic) >= 'a' .And. comps_str(ic:ic) <= 'z' ) &
            comps_str(ic:ic) = Achar(Iachar(comps_str(ic:ic)) - 32)
    End Do

    ! presence flags (1/0), fixed U,V,W,T,C block order on disk (see read_recycle_frame_comp)
    colU = Merge(1,0, Index(comps_str,'U') > 0)
    colV = Merge(1,0, Index(comps_str,'V') > 0)
    colW = Merge(1,0, Index(comps_str,'W') > 0)
    colT = Merge(1,0, Index(comps_str,'T') > 0)
    colC = Merge(1,0, Index(comps_str,'C') > 0)
    If ( colU == 0 .Or. colV == 0 .Or. colW == 0 ) Stop &
         'ERROR: inflow_recycle_file must contain U,V,W (set slice_comps="UVW..." on the donor run)'

  End Subroutine read_recycle_meta

  !> Mean-profile seed for inflow_type==2 (recycled precursor inflow): every rank independently reads the donor's first frame and z-averages U(y), mirroring read_mean_profile's per-rank-independent read of a small shared file (here the read is one frame, not the whole donor file); R11/R22/R33/R12 are left at zero since sem_fluctuation is never called under inflow_type==2 -- these only feed mean_profile_U, used by genGridandIC to seed the IC mean
  Subroutine read_recycle_mean_profile

    Integer(Int32) :: ncomp, n1, n1_v, n2, nsnaps, colU, colV, colW, colT, colC, unit_in, ios, jy
    Logical :: v_native
    Real(Int64), Allocatable :: frame_U(:,:)
    Character(300) :: fname

    Call read_recycle_meta(ncomp, n1, n1_v, n2, nsnaps, colU, colV, colW, colT, colC, v_native)

    If ( n1 /= nym_global ) Stop 'ERROR: inflow_recycle_file ny (n1) does not match this run''s nym_global'
    If ( n2 /= nzm_global ) Stop 'ERROR: inflow_recycle_file nz (n2) does not match this run''s nzm_global'

    n_profile = n1
    Allocate( prof_y(n_profile), prof_U(n_profile) )
    Allocate( prof_R11(n_profile), prof_R22(n_profile), prof_R33(n_profile), prof_R12(n_profile) )
    prof_y = ym_global
    prof_R11 = 0d0;  prof_R22 = 0d0;  prof_R33 = 0d0;  prof_R12 = 0d0

    ! only rank 0 pays the O(n1*n2) transient allocation+read for the donor's first frame's
    ! U block (the first block on disk, see the per-component block order written by
    ! dopamine-ESEM); the result (prof_U, n1 reals) is tiny, so broadcast it instead
    If ( myid == 0 ) Then
       Allocate( frame_U(n1,n2) )
       Write(fname,'(A,A)') Trim(inflow_recycle_file), '.bin'
       Open(newunit=unit_in, file=Trim(fname), access='stream', form='unformatted', &
            status='old', action='read', iostat=ios)
       If ( ios /= 0 ) Stop 'ERROR: cannot open inflow_recycle_file data (expected <inflow_recycle_file>.bin)'
       Read(unit_in) frame_U   ! U is always the first block of the first frame
       Close(unit_in)
       Do jy = 1, n1
          prof_U(jy) = Sum(frame_U(jy,:)) / Real(n2,8)
       End Do
       Deallocate(frame_U)
    End If

    Call Mpi_bcast( prof_U, n_profile, MPI_real8, 0, MPI_COMM_WORLD, ierr )

  End Subroutine read_recycle_mean_profile

  !> One-time setup for inflow_type==2's runtime BC buffers: no-op on every rank except the one owning the x=1 face (row==0), which opens a persistent stream handle on the donor .bin and primes the bracketing time-interpolation frames (rec_lo/rec_hi) for this rank's own local z-slab only
  Subroutine init_inflow_recycle

    Integer(Int32) :: ios, u_times
    Character(300) :: fname

    rec_unit = -1
    rec_active = .False.
    ! this rank owns the x=1 face iff it's in row 0 of the p_row/p_col grid (mirrors decomp's x_periodic_partner, avoiding a new decomp.f90 dependency here since sem.f90 is also compiled standalone by the tests/verify_* targets)
    If ( myid/p_col /= 0 ) Return

    Call read_recycle_meta(rec_ncomp, rec_n1, rec_n1_v, rec_n2_global, rec_nsnaps, &
         rec_col_U, rec_col_V, rec_col_W, rec_col_T, rec_col_C, rec_v_native)

    If ( rec_n1 /= nym_global ) Stop 'ERROR: inflow_recycle_file ny (n1) does not match this run''s nym_global'
    If ( rec_v_native ) Then
       If ( rec_n1_v /= ny_global ) &
            Stop 'ERROR: inflow_recycle_file n1_V (V''s own y-face count) does not match this run''s ny_global'
    End If
    If ( rec_n2_global /= nzm_global ) Stop 'ERROR: inflow_recycle_file nz (n2) does not match this run''s nzm_global'
    If ( rec_nsnaps < 2 ) Stop 'ERROR: inflow_recycle_file must contain at least 2 snapshots to interpolate in time'
    If ( boussinesq_flag >= 1 .And. rec_col_T == 0 ) Stop 'ERROR: boussinesq_flag>=1 with ' // &
         'inflow_type=2 requires T in the donor slice (donor slice_comps must include ''T'')'
    If ( sediment_flag >= 1 .And. rec_col_C == 0 ) Stop 'ERROR: sediment_flag>=1 with ' // &
         'inflow_type=2 requires C in the donor slice (donor slice_comps must include ''C'')'

    Allocate( rec_times(rec_nsnaps) )
    Write(fname,'(A,A)') Trim(inflow_recycle_file), '_times.bin'
    Open(newunit=u_times, file=Trim(fname), access='stream', form='unformatted', &
         status='old', action='read', iostat=ios)
    If ( ios /= 0 ) Stop 'ERROR: cannot open inflow_recycle_file times (expected <inflow_recycle_file>_times.bin)'
    Read(u_times) rec_times
    Close(u_times)

    ! this rank's local donor z-slab: same global 1-based interior cc/face convention the writer used (kg_g = kg1_global(myid)+ka-2)
    rec_k1 = kg1_global(myid)
    rec_k2 = kg2_global(myid) - 2

    ! each component's own block, fixed U,V,W,T,C disk order (only present ones actually occupy space)
    rec_off_U = 0_Int64
    rec_off_V = rec_off_U + 8_Int64*Int(rec_n1,  Int64)*Int(rec_n2_global,Int64)
    rec_off_W = rec_off_V + 8_Int64*Int(rec_n1_v,Int64)*Int(rec_n2_global,Int64)
    rec_off_T = rec_off_W + 8_Int64*Int(rec_n1,  Int64)*Int(rec_n2_global,Int64)
    rec_off_C = rec_off_T + Merge(8_Int64*Int(rec_n1,Int64)*Int(rec_n2_global,Int64), 0_Int64, rec_col_T==1)
    rec_frame_bytes = rec_off_C + Merge(8_Int64*Int(rec_n1,Int64)*Int(rec_n2_global,Int64), 0_Int64, rec_col_C==1)

    Write(fname,'(A,A)') Trim(inflow_recycle_file), '.bin'
    Open(newunit=rec_unit, file=Trim(fname), access='stream', form='unformatted', &
         status='old', action='read', iostat=ios)
    If ( ios /= 0 ) Stop 'ERROR: cannot open inflow_recycle_file data (expected <inflow_recycle_file>.bin)'
    rec_active = .True.

    Allocate( rec_lo_U(rec_n1,   rec_k2-rec_k1+1) );  Allocate( rec_hi_U(rec_n1,   rec_k2-rec_k1+1) )
    Allocate( rec_lo_V(rec_n1_v, rec_k2-rec_k1+1) );  Allocate( rec_hi_V(rec_n1_v, rec_k2-rec_k1+1) )
    Allocate( rec_lo_W(rec_n1,   rec_k2-rec_k1+1) );  Allocate( rec_hi_W(rec_n1,   rec_k2-rec_k1+1) )
    Allocate( rec_lo_T(Merge(rec_n1,0,rec_col_T==1), Merge(rec_k2-rec_k1+1,0,rec_col_T==1)) )
    Allocate( rec_hi_T(Merge(rec_n1,0,rec_col_T==1), Merge(rec_k2-rec_k1+1,0,rec_col_T==1)) )
    Allocate( rec_lo_C(Merge(rec_n1,0,rec_col_C==1), Merge(rec_k2-rec_k1+1,0,rec_col_C==1)) )
    Allocate( rec_hi_C(Merge(rec_n1,0,rec_col_C==1), Merge(rec_k2-rec_k1+1,0,rec_col_C==1)) )
    rec_idx_lo = -1;  rec_idx_hi = -1

    !$acc update device(rec_col_U, rec_col_V, rec_col_W, rec_col_T, rec_col_C, rec_n1, rec_n1_v, rec_v_native)

    Call update_inflow_recycle(t)   ! prime rec_lo_*/rec_hi_*/rec_frac before the first BC application

    Write(*,'(A,I0,A,A,A,I0,A)') ' Rank ', myid, ': recycled precursor inflow from ', &
         Trim(inflow_recycle_file), ' (', rec_nsnaps, ' donor snapshots)'

  End Subroutine init_inflow_recycle

  !> Read one contiguous donor-z span (z0..z0+count-1, 1-based global donor z index, no shift/wrap applied by this routine) of one component's block (comp_off = that component's byte offset within a frame, comp_n1 = its own native n1) into buf(:,zoff:zoff+count-1); each component's own (n1_c,n2) block is plain column-major (n1_c fastest), so any single z-span is a contiguous disk read
  Subroutine read_recycle_span(frame_idx, comp_off, comp_n1, z0, count, buf, zoff)

    Integer(Int32), Intent(In)    :: frame_idx, comp_n1, z0, count, zoff
    Integer(Int64), Intent(In)    :: comp_off
    Real   (Int64), Intent(InOut) :: buf(:,:)

    Integer(Int64) :: byte_pos

    If ( count < 1 ) Return
    byte_pos = Int(frame_idx-1,Int64)*rec_frame_bytes + comp_off + &
               Int(z0-1,Int64)*8_Int64*Int(comp_n1,Int64) + 1_Int64
    Read(rec_unit, pos=byte_pos) buf(:,zoff:zoff+count-1)

  End Subroutine read_recycle_span

  !> Read one component's block from one donor frame's local z-slab, circularly shifted by shift_z donor z-cells (spanwise shift, see update_inflow_recycle); this rank's unshifted span is [rec_k1,rec_k2], so the shifted span is that window rotated by shift_z within the donor's periodic z range [1,rec_n2_global] -- at most one wrap, split into two contiguous reads when it crosses the donor's z boundary
  Subroutine read_recycle_frame_comp(frame_idx, shift_z, comp_off, comp_n1, buf)

    Integer(Int32), Intent(In)    :: frame_idx, shift_z, comp_n1
    Integer(Int64), Intent(In)    :: comp_off
    Real   (Int64), Intent(InOut) :: buf(:,:)

    Integer(Int32) :: nzloc, z_start, nA

    nzloc   = rec_k2 - rec_k1 + 1
    z_start = Modulo( rec_k1 - 1 + shift_z, rec_n2_global ) + 1

    If ( z_start + nzloc - 1 <= rec_n2_global ) Then
       Call read_recycle_span(frame_idx, comp_off, comp_n1, z_start, nzloc, buf, 1)
    Else
       nA = rec_n2_global - z_start + 1
       Call read_recycle_span(frame_idx, comp_off, comp_n1, z_start, nA,       buf, 1)
       Call read_recycle_span(frame_idx, comp_off, comp_n1, 1,       nzloc-nA, buf, nA+1)
    End If

  End Subroutine read_recycle_frame_comp

  !> Read all active components' blocks for one donor frame (see read_recycle_frame_comp); dst_lo=.True. writes into the rec_lo_* buffers, .False. into rec_hi_*
  Subroutine read_recycle_frame(frame_idx, shift_z, dst_lo)

    Integer(Int32), Intent(In) :: frame_idx, shift_z
    Logical,        Intent(In) :: dst_lo

    If ( dst_lo ) Then
       Call read_recycle_frame_comp(frame_idx, shift_z, rec_off_U, rec_n1,   rec_lo_U)
       Call read_recycle_frame_comp(frame_idx, shift_z, rec_off_V, rec_n1_v, rec_lo_V)
       Call read_recycle_frame_comp(frame_idx, shift_z, rec_off_W, rec_n1,   rec_lo_W)
       If ( rec_col_T == 1 ) Call read_recycle_frame_comp(frame_idx, shift_z, rec_off_T, rec_n1, rec_lo_T)
       If ( rec_col_C == 1 ) Call read_recycle_frame_comp(frame_idx, shift_z, rec_off_C, rec_n1, rec_lo_C)
    Else
       Call read_recycle_frame_comp(frame_idx, shift_z, rec_off_U, rec_n1,   rec_hi_U)
       Call read_recycle_frame_comp(frame_idx, shift_z, rec_off_V, rec_n1_v, rec_hi_V)
       Call read_recycle_frame_comp(frame_idx, shift_z, rec_off_W, rec_n1,   rec_hi_W)
       If ( rec_col_T == 1 ) Call read_recycle_frame_comp(frame_idx, shift_z, rec_off_T, rec_n1, rec_hi_T)
       If ( rec_col_C == 1 ) Call read_recycle_frame_comp(frame_idx, shift_z, rec_off_C, rec_n1, rec_hi_C)
    End If

  End Subroutine read_recycle_frame

  !> Advance the recycled-inflow bracket (rec_lo/rec_hi/rec_frac) to time t_now; no-op except on the rank owning the x=1 face. Re-reads donor frames from disk only when the bracket or the spanwise shift actually changes, and only ever reads the frame(s) needed rather than the whole donor plane -- when the bracket simply slides forward (same loop pass, same shift), the old "hi" frame becomes the new "lo" and only one fresh frame is read
  Subroutine update_inflow_recycle(t_now)

    Real(Int64), Intent(In) :: t_now

    Real(Int64) :: t0, t1, duration, t_target, u
    Integer(Int32) :: ilo, ihi, lo, hi, mid, loop_idx, shift_z
    Logical :: changed

    If ( .Not. rec_active ) Return

    t0 = rec_times(1);  t1 = rec_times(rec_nsnaps)
    duration = t1 - t0

    t_target = t_now + inflow_recycle_t_offset
    loop_idx = 0
    If ( inflow_recycle_loop == 1 .And. duration > 0d0 ) Then
       loop_idx = Floor( (t_target - t0) / duration )
       t_target = t0 + Modulo(t_target - t0, duration)
    Else
       t_target = Min( Max(t_target, t0), t1 )
    End If

    ! new spanwise shift each time we start a fresh pass through the donor timeline (deterministic hash of (seed,loop_idx), same tool sem.f90 already uses to re-randomize eddies each recycle -- restart-safe, no RNG state to persist)
    shift_z = 0
    If ( inflow_recycle_loop == 1 .And. inflow_recycle_shift_z == 1 .And. rec_n2_global > 1 ) Then
       If ( loop_idx /= rec_loop_idx ) Then
          u = hash_uniform(inflow_recycle_seed, loop_idx, 0, 97)
          rec_shift_z  = Int( u * Real(rec_n2_global,8) )
          rec_shift_z  = Min( Max(rec_shift_z,0), rec_n2_global-1 )
          rec_loop_idx = loop_idx
       End If
       shift_z = rec_shift_z
    End If

    ! binary search rec_times(1:rec_nsnaps) (monotone by construction) for the bracketing pair
    lo = 1;  hi = rec_nsnaps
    Do While ( hi - lo > 1 )
       mid = (lo+hi)/2
       If ( rec_times(mid) <= t_target ) Then
          lo = mid
       Else
          hi = mid
       End If
    End Do
    ilo = lo;  ihi = hi
    If ( ilo == ihi ) ihi = Min(ilo+1, rec_nsnaps)

    changed = .False.
    If ( ilo /= rec_idx_lo .Or. ihi /= rec_idx_hi .Or. shift_z /= rec_idx_shift ) Then
       If ( ilo == rec_idx_hi .And. shift_z == rec_idx_shift ) Then
          ! bracket slid forward, same shift: reuse the old "hi" as the new "lo", don't re-read it
          rec_lo_U = rec_hi_U;  rec_lo_V = rec_hi_V;  rec_lo_W = rec_hi_W
          If ( rec_col_T == 1 ) rec_lo_T = rec_hi_T
          If ( rec_col_C == 1 ) rec_lo_C = rec_hi_C
          Call read_recycle_frame(ihi, shift_z, .False.)
       Else
          Call read_recycle_frame(ilo, shift_z, .True.)    ! first call, a jump, or a fresh shift at a loop wrap
          Call read_recycle_frame(ihi, shift_z, .False.)
       End If
       rec_idx_lo = ilo;  rec_idx_hi = ihi;  rec_idx_shift = shift_z
       changed = .True.
    End If

    If ( rec_times(ihi) > rec_times(ilo) ) Then
       rec_frac = ( t_target - rec_times(ilo) ) / ( rec_times(ihi) - rec_times(ilo) )
    Else
       rec_frac = 0d0
    End If
    rec_frac = Min( Max(rec_frac,0d0), 1d0 )

    If ( changed ) Then
       !$acc update device(rec_lo_U, rec_hi_U, rec_lo_V, rec_hi_V, rec_lo_W, rec_hi_W, rec_lo_T, rec_hi_T, rec_lo_C, rec_hi_C)
    End If
    !$acc update device(rec_frac)

  End Subroutine update_inflow_recycle

  !> Time-interpolated recycled-inflow value for component comp (1=U,2=V,3=W,4=T,5=C) at this rank's local array position (j,k). Each component now reads its OWN native-grid donor block (see dopamine-ESEM's per-component staggered-grid writer), so this is an EXACT lookup, no half-cell interpolation/index-shift approximation: U/W use apply_inflow_bc_x's (yg,zg)-ghost-inclusive index space (j-1/k-1, matching rec_n1/rec_k1..rec_k2's cell-centre convention exactly like before); V uses its OWN (y,zg) index space -- j directly (V's local array has no y-ghost, so j already IS the y-face index, no -1 shift) and k-1 in z (V is still z ghost-cell-centred); W keeps U's y-index convention (j-1) but uses k directly (its local array has no z-ghost, so k already IS the z-face index, matching this donor's z(:) left-face-of-cell storage)
  Real(Int64) Function recycle_value(comp, j, k) Result(val)
    !$acc routine seq

    Integer(Int32), Intent(In) :: comp, j, k

    Integer(Int32) :: jy, kz, nzloc

    nzloc = rec_k2 - rec_k1 + 1

    Select Case (comp)
    Case (1)   ! U: (yg,zg) ghost cell-centre
       jy = Max( 1, Min(rec_n1, j-1) );  kz = Max( 1, Min(nzloc, k-1) )
       val = rec_lo_U(jy,kz) + rec_frac*( rec_hi_U(jy,kz) - rec_lo_U(jy,kz) )
    Case (2)   ! V: native-grid donor -- (y,zg) y-face (no shift); legacy/generic donor -- V was
               ! cell-centre-interpolated onto the same grid as U/W, so index it the same way (j-1)
       If ( rec_v_native ) Then
          jy = Max( 1, Min(rec_n1_v, j) )
       Else
          jy = Max( 1, Min(rec_n1, j-1) )
       End If
       kz = Max( 1, Min(nzloc, k-1) )
       val = rec_lo_V(jy,kz) + rec_frac*( rec_hi_V(jy,kz) - rec_lo_V(jy,kz) )
    Case (3)   ! W: (yg,z) y ghost cell-centre, z-face (no shift)
       jy = Max( 1, Min(rec_n1, j-1) );  kz = Max( 1, Min(nzloc, k) )
       val = rec_lo_W(jy,kz) + rec_frac*( rec_hi_W(jy,kz) - rec_lo_W(jy,kz) )
    Case (4)   ! T: cell-centre, same convention as U
       jy = Max( 1, Min(rec_n1, j-1) );  kz = Max( 1, Min(nzloc, k-1) )
       val = rec_lo_T(jy,kz) + rec_frac*( rec_hi_T(jy,kz) - rec_lo_T(jy,kz) )
    Case Default   ! C: cell-centre, same convention as U
       jy = Max( 1, Min(rec_n1, j-1) );  kz = Max( 1, Min(nzloc, k-1) )
       val = rec_lo_C(jy,kz) + rec_frac*( rec_hi_C(jy,kz) - rec_lo_C(jy,kz) )
    End Select

  End Function recycle_value

  !> Read the SEM inflow mean-temperature profile (two columns: y T); optional companion to read_mean_profile
  Subroutine read_mean_profile_T

    Integer(Int32) :: unit_in, ios, n
    Real   (Int64) :: col(2)
    Character(300) :: line

    If ( myid == 0 ) Write(*,'(A,A)') ' Reading SEM inflow temperature profile from ', Trim(inflow_temperature_file)

    Open(newunit=unit_in, file=Trim(inflow_temperature_file), status='old', action='read', iostat=ios)
    If ( ios /= 0 ) Stop 'ERROR: cannot open inflow_temperature_file'

    n_profile_T = 0
    Do
       Read(unit_in,'(A)',iostat=ios) line
       If ( ios /= 0 ) Exit
       line = Adjustl(line)
       If ( Len_trim(line) == 0 .Or. line(1:1) == '#' ) Cycle
       n_profile_T = n_profile_T + 1
    End Do
    If ( n_profile_T < 2 ) Stop 'ERROR: inflow_temperature_file: fewer than 2 data rows found'

    Allocate( prof_y_T(n_profile_T), prof_T(n_profile_T) )

    Rewind(unit_in)
    n = 0
    Do
       Read(unit_in,'(A)',iostat=ios) line
       If ( ios /= 0 ) Exit
       line = Adjustl(line)
       If ( Len_trim(line) == 0 .Or. line(1:1) == '#' ) Cycle
       n = n + 1
       Read(line,*,iostat=ios) col(1:2)
       If ( ios /= 0 ) Stop 'ERROR: inflow_temperature_file: failed to parse a data row (need 2 columns: y T)'
       prof_y_T(n) = col(1)
       prof_T(n)   = col(2)
    End Do
    Close(unit_in)

  End Subroutine read_mean_profile_T

  ! Read the mean-velocity / Reynolds-stress reference profile (unchanged
  ! format/semantics from the original SEM implementation).
  Subroutine read_mean_profile

    Integer(Int32) :: unit_in, ios, n
    Real   (Int64) :: col(8)
    Character(300) :: line

    If ( myid == 0 ) Write(*,'(A,A)') ' Reading SEM inflow profile from ', Trim(inflow_profile_file)

    ! every rank reads the (small, shared) profile file independently; free-form text, skip blank/'#' lines, each row: y U V W uu vv ww uv [uw vw] (extra cols read but unused, see sem_fluctuation's R13=R23=0 assumption)
    Open(newunit=unit_in, file=Trim(inflow_profile_file), status='old', action='read', iostat=ios)
    If ( ios /= 0 ) Stop 'ERROR: cannot open inflow_profile_file'

    n_profile = 0
    Do
       Read(unit_in,'(A)',iostat=ios) line
       If ( ios /= 0 ) Exit
       line = Adjustl(line)
       If ( Len_trim(line) == 0 .Or. line(1:1) == '#' ) Cycle
       n_profile = n_profile + 1
    End Do
    If ( n_profile < 2 ) Stop 'ERROR: inflow_profile_file: fewer than 2 data rows found'

    Allocate( prof_y(n_profile), prof_U(n_profile) )
    Allocate( prof_R11(n_profile), prof_R22(n_profile), prof_R33(n_profile), prof_R12(n_profile) )

    Rewind(unit_in)
    n = 0
    Do
       Read(unit_in,'(A)',iostat=ios) line
       If ( ios /= 0 ) Exit
       line = Adjustl(line)
       If ( Len_trim(line) == 0 .Or. line(1:1) == '#' ) Cycle
       n = n + 1
       Read(line,*,iostat=ios) col(1:8)
       If ( ios /= 0 ) Stop 'ERROR: inflow_profile_file: failed to parse a data row (need >=8 columns)'
       prof_y(n)   = col(1)   ! y
       prof_U(n)   = col(2)   ! U
       prof_R11(n) = col(5)   ! uu
       prof_R22(n) = col(6)   ! vv
       prof_R33(n) = col(7)   ! ww
       prof_R12(n) = col(8)   ! uv

       ! PSD check (R13=R23=0 assumed): R11,R33>=0 and the 2x2 leading minor R11*R22-R12^2>=0
       If ( prof_R11(n) < 0d0 .Or. prof_R33(n) < 0d0 .Or. &
            prof_R11(n)*prof_R22(n) - prof_R12(n)**2 < -1d-10*Max(prof_R11(n)*prof_R22(n),1d-30) ) &
            Stop 'ERROR: inflow_profile_file: Reynolds-stress tensor is not positive-semidefinite at some y'
    End Do
    Close(unit_in)

    If ( prof_y(1) > y_global(1)+1d-10 .Or. prof_y(n_profile) < y_global(ny_global)-1d-10 ) Then
       If ( myid == 0 ) Write(*,'(A)') &
            ' WARNING: inflow_profile_file does not cover the full domain height -- clamping at the endpoints'
    End If

  End Subroutine read_mean_profile

  !> Read a wind-tunnel style profile (sem_profile_format==1): header "# z U Iu Iv Iw [Lux Luy Luz Lvx Lvy Lvz Lwx Lwy Lwz]", intensities converted to Reynolds stresses (uv=0, no shear data), length scales optional/partial with the tiered fallback documented in sem.md
  Subroutine read_mean_profile_TI

    Integer(Int32), Parameter :: n_slot = 14
    Character(4) :: slot_name(n_slot) = (/ &
         'z   ', 'U   ', 'Iu  ', 'Iv  ', 'Iw  ', &
         'Lux ', 'Luy ', 'Luz ', 'Lvx ', 'Lvy ', 'Lvz ', 'Lwx ', 'Lwy ', 'Lwz ' /)
    Integer(Int32) :: col_of_slot(n_slot), ncols
    Integer(Int32) :: unit_in, ios, n
    Character(2000) :: line
    Real   (Int64), Allocatable :: row(:)
    Real   (Int64) :: Iu_v, Iv_v, Iw_v
    Logical :: any_Lscale

    If ( myid == 0 ) Write(*,'(A,A)') ' Reading wind-tunnel TI SEM inflow profile from ', Trim(inflow_profile_file)

    Open(newunit=unit_in, file=Trim(inflow_profile_file), status='old', action='read', iostat=ios)
    If ( ios /= 0 ) Stop 'ERROR: cannot open inflow_profile_file'

    ! ---- header: first non-blank line, must start with '#' and list column names ----
    Do
       Read(unit_in,'(A)',iostat=ios) line
       If ( ios /= 0 ) Stop 'ERROR: inflow_profile_file: missing "# z U Iu Iv Iw ..." header'
       line = Adjustl(line)
       If ( Len_trim(line) == 0 ) Cycle
       If ( line(1:1) /= '#' ) Stop 'ERROR: inflow_profile_file: first line must be "# ..." header'
       Call parse_TI_header(line(2:), slot_name, n_slot, col_of_slot, ncols)
       Exit
    End Do

    If ( Any( col_of_slot(1:5) == 0 ) ) &
         Stop 'ERROR: inflow_profile_file: header must include z (or y), U, Iu, Iv, Iw'
    any_Lscale = Any( col_of_slot(6:14) > 0 )

    ! ---- count data rows ----
    n_profile = 0
    Do
       Read(unit_in,'(A)',iostat=ios) line
       If ( ios /= 0 ) Exit
       line = Adjustl(line)
       If ( Len_trim(line) == 0 .Or. line(1:1) == '#' ) Cycle
       n_profile = n_profile + 1
    End Do
    If ( n_profile < 2 ) Stop 'ERROR: inflow_profile_file: fewer than 2 data rows found'

    Allocate( prof_y(n_profile), prof_U(n_profile) )
    Allocate( prof_R11(n_profile), prof_R22(n_profile), prof_R33(n_profile), prof_R12(n_profile) )

    ! sig_* is populated whenever measured/ratio-derived (any_Lscale) or mixing-length-derived data will be used; the mixing-length fallback needs sem_use_esem==1 (classical SEM has no normalisation for inhomogeneous length scales, see init_inflow)
    If ( any_Lscale .Or. sem_use_esem == 1 ) Then
       n_sigma = n_profile
       Allocate( sig_y(n_sigma) )
       Allocate( sig_ux(n_sigma), sig_uy(n_sigma), sig_uz(n_sigma) )
       Allocate( sig_vx(n_sigma), sig_vy(n_sigma), sig_vz(n_sigma) )
       Allocate( sig_wx(n_sigma), sig_wy(n_sigma), sig_wz(n_sigma) )
    End If

    Allocate( row(ncols) )

    Rewind(unit_in)
    Do   ! re-skip the header
       Read(unit_in,'(A)') line
       line = Adjustl(line)
       If ( Len_trim(line) > 0 ) Exit
    End Do

    n = 0
    Do
       Read(unit_in,'(A)',iostat=ios) line
       If ( ios /= 0 ) Exit
       line = Adjustl(line)
       If ( Len_trim(line) == 0 .Or. line(1:1) == '#' ) Cycle
       n = n + 1
       Read(line,*,iostat=ios) row(1:ncols)
       If ( ios /= 0 ) Stop 'ERROR: inflow_profile_file: failed to parse a data row'

       prof_y(n) = row(col_of_slot(1))
       prof_U(n) = row(col_of_slot(2))
       Iu_v = row(col_of_slot(3)); Iv_v = row(col_of_slot(4)); Iw_v = row(col_of_slot(5))
       If ( Iu_v < 0d0 .Or. Iv_v < 0d0 .Or. Iw_v < 0d0 ) &
            Stop 'ERROR: inflow_profile_file: turbulence intensities must be non-negative'
       prof_R11(n) = (Iu_v*prof_U(n))**2
       prof_R22(n) = (Iv_v*prof_U(n))**2
       prof_R33(n) = (Iw_v*prof_U(n))**2
       prof_R12(n) = 0d0   ! shear stress not available from typical wind-tunnel TI data

       If ( n_sigma > 0 ) sig_y(n) = prof_y(n)
       If ( any_Lscale ) Then
          Call fill_length_scale( col_of_slot(6),  col_of_slot(7),  col_of_slot(8),  row, ncols, sig_ux(n), sig_uy(n), sig_uz(n) )
          Call fill_length_scale( col_of_slot(9),  col_of_slot(10), col_of_slot(11), row, ncols, sig_vx(n), sig_vy(n), sig_vz(n) )
          Call fill_length_scale( col_of_slot(12), col_of_slot(13), col_of_slot(14), row, ncols, sig_wx(n), sig_wy(n), sig_wz(n) )
       End If
    End Do
    Close(unit_in)
    Deallocate(row)

    ! no L** column at all: derive an inhomogeneous length-scale profile from mixing-length theory instead of the old isotropic sem_length_scale fallback, see derive_mixing_length
    If ( .Not. any_Lscale .And. n_sigma > 0 ) Call derive_mixing_length

    If ( prof_y(1) > y_global(1)+1d-10 .Or. prof_y(n_profile) < y_global(ny_global)-1d-10 ) Then
       If ( myid == 0 ) Write(*,'(A)') &
            ' WARNING: inflow_profile_file does not cover the full domain height -- clamping at the endpoints'
    End If

  End Subroutine read_mean_profile_TI

  !> Prandtl mixing-length fallback for sem_profile_format=1 when no L** column is measured: L(y) = sigma_u(y)/|dU/dy|, sigma_u = sqrt(R11); same L used for all three components' x-scale (no independent shear estimate for v,w), y/z-scales via the usual sem_Lscale_ratio_y/z; clamped to [sem_length_scale, outer_frac*profile height] since |dU/dy|->0 blows the estimate up near a free-stream/outer edge
  Subroutine derive_mixing_length

    Real(Int64), Parameter :: outer_frac = 0.3d0
    Real(Int64) :: dUdy, Lcap, Lmix
    Integer(Int32) :: i

    Lcap = outer_frac * ( prof_y(n_profile) - prof_y(1) )

    Do i = 1, n_profile
       If ( i == 1 ) Then
          dUdy = ( prof_U(2) - prof_U(1) ) / ( prof_y(2) - prof_y(1) )
       Else If ( i == n_profile ) Then
          dUdy = ( prof_U(n_profile) - prof_U(n_profile-1) ) / ( prof_y(n_profile) - prof_y(n_profile-1) )
       Else
          dUdy = ( prof_U(i+1) - prof_U(i-1) ) / ( prof_y(i+1) - prof_y(i-1) )
       End If

       Lmix = Sqrt(prof_R11(i)) / Max( Abs(dUdy), 1d-12 )
       Lmix = Min( Max( Lmix, sem_length_scale ), Lcap )

       sig_ux(i) = Lmix;  sig_uy(i) = sem_Lscale_ratio_y*Lmix;  sig_uz(i) = sem_Lscale_ratio_z*Lmix
       sig_vx(i) = Lmix;  sig_vy(i) = sem_Lscale_ratio_y*Lmix;  sig_vz(i) = sem_Lscale_ratio_z*Lmix
       sig_wx(i) = Lmix;  sig_wy(i) = sem_Lscale_ratio_y*Lmix;  sig_wz(i) = sem_Lscale_ratio_z*Lmix
    End Do

  End Subroutine derive_mixing_length

  !> Tokenize a TI-format header line into slot->column-index map; unknown/duplicate column names are fatal
  Subroutine parse_TI_header(hdr, slot_name, n_slot, col_of_slot, ncols)

    Character(*),   Intent(In)  :: hdr
    Character(*),   Intent(In)  :: slot_name(:)
    Integer(Int32), Intent(In)  :: n_slot
    Integer(Int32), Intent(Out) :: col_of_slot(n_slot)
    Integer(Int32), Intent(Out) :: ncols

    Character(2000) :: buf
    Character(20)   :: tok
    Character(1), Parameter :: tab = Achar(9)
    Integer(Int32)  :: p, islot

    col_of_slot = 0
    ncols = 0
    buf = hdr

    Do
       Do While ( Len_trim(buf) > 0 )   ! strip leading spaces/tabs (Adjustl only moves spaces, not tabs)
          If ( buf(1:1) /= ' ' .And. buf(1:1) /= tab ) Exit
          buf = buf(2:)
       End Do
       If ( Len_trim(buf) == 0 ) Exit
       p = Scan(buf, ' '//tab)            ! first whitespace char, i.e. end of this token
       If ( p == 0 ) p = Len_trim(buf) + 1
       tok = buf(1:p-1)
       buf = buf(p:)
       ncols = ncols + 1

       islot = token_slot(tok, slot_name, n_slot)
       If ( islot < 0 ) Stop 'ERROR: inflow_profile_file: unrecognized header column name "'//Trim(tok)//'"'
       If ( col_of_slot(islot) /= 0 ) Stop 'ERROR: inflow_profile_file: duplicate header column name'
       col_of_slot(islot) = ncols
    End Do

  End Subroutine parse_TI_header

  !> Map a header token to its slot id (1..n_slot), 'y' accepted as an alias for slot 1 ('z'); -1 if unrecognized
  Integer(Int32) Function token_slot(tok, slot_name, n_slot) Result(islot)

    Character(*),   Intent(In) :: tok
    Character(*),   Intent(In) :: slot_name(:)
    Integer(Int32), Intent(In) :: n_slot
    Integer(Int32) :: m

    If ( Trim(Adjustl(tok)) == 'y' ) Then
       islot = 1
       Return
    End If
    Do m = 1, n_slot
       If ( Trim(Adjustl(tok)) == Trim(slot_name(m)) ) Then
          islot = m
          Return
       End If
    End Do
    islot = -1

  End Function token_slot

  !> Fill one component's (Lox,Loy,Loz) length scales with the tiered fallback (measured > ratio*Lox > sem_length_scale), see sem.md
  Subroutine fill_length_scale(cx, cy, cz, row, ncols, Lx_out, Ly_out, Lz_out)

    Integer(Int32), Intent(In)  :: cx, cy, cz, ncols
    Real   (Int64), Intent(In)  :: row(ncols)
    Real   (Int64), Intent(Out) :: Lx_out, Ly_out, Lz_out

    If ( cx > 0 ) Then
       Lx_out = row(cx)
    Else
       Lx_out = sem_length_scale
    End If

    If ( cy > 0 ) Then
       Ly_out = row(cy)
    Else If ( cx > 0 ) Then
       Ly_out = sem_Lscale_ratio_y * row(cx)
    Else
       Ly_out = sem_length_scale
    End If

    If ( cz > 0 ) Then
       Lz_out = row(cz)
    Else If ( cx > 0 ) Then
       Lz_out = sem_Lscale_ratio_z * row(cx)
    Else
       Lz_out = sem_length_scale
    End If

    If ( Lx_out <= 0d0 .Or. Ly_out <= 0d0 .Or. Lz_out <= 0d0 ) &
         Stop 'ERROR: inflow_profile_file: length scales must be strictly positive'

  End Subroutine fill_length_scale

  !> Read the optional inhomogeneous eddy length-scale profile (ESEM Section 4.1); empty sem_sigma_file => homogeneous mode (n_sigma=0)
  Subroutine read_sigma_profile

    Integer(Int32) :: unit_in, ios, n
    Real   (Int64) :: col(10)
    Character(300) :: line

    If ( Len_trim(sem_sigma_file) == 0 ) Return   ! keep whatever read_mean_profile_TI already derived (n_sigma stays 0 if none)

    If ( myid == 0 ) Write(*,'(A,A)') ' Reading inhomogeneous SEM length-scale profile from ', Trim(sem_sigma_file)

    If ( Allocated(sig_y) ) &
         Deallocate( sig_y, sig_ux, sig_uy, sig_uz, sig_vx, sig_vy, sig_vz, sig_wx, sig_wy, sig_wz )

    Open(newunit=unit_in, file=Trim(sem_sigma_file), status='old', action='read', iostat=ios)
    If ( ios /= 0 ) Stop 'ERROR: cannot open sem_sigma_file'

    n_sigma = 0
    Do
       Read(unit_in,'(A)',iostat=ios) line
       If ( ios /= 0 ) Exit
       line = Adjustl(line)
       If ( Len_trim(line) == 0 .Or. line(1:1) == '#' ) Cycle
       n_sigma = n_sigma + 1
    End Do
    If ( n_sigma < 2 ) Stop 'ERROR: sem_sigma_file: fewer than 2 data rows found'

    Allocate( sig_y(n_sigma) )
    Allocate( sig_ux(n_sigma), sig_uy(n_sigma), sig_uz(n_sigma) )
    Allocate( sig_vx(n_sigma), sig_vy(n_sigma), sig_vz(n_sigma) )
    Allocate( sig_wx(n_sigma), sig_wy(n_sigma), sig_wz(n_sigma) )

    Rewind(unit_in)
    n = 0
    Do
       Read(unit_in,'(A)',iostat=ios) line
       If ( ios /= 0 ) Exit
       line = Adjustl(line)
       If ( Len_trim(line) == 0 .Or. line(1:1) == '#' ) Cycle
       n = n + 1
       Read(line,*,iostat=ios) col(1:10)
       If ( ios /= 0 ) Stop 'ERROR: sem_sigma_file: failed to parse a data row (need 10 columns)'
       If ( Minval(col(2:10)) <= 0d0 ) Stop 'ERROR: sem_sigma_file: length scales must be strictly positive'
       sig_y(n)  = col(1)
       sig_ux(n) = col(2); sig_uy(n) = col(3); sig_uz(n) = col(4)
       sig_vx(n) = col(5); sig_vy(n) = col(6); sig_vz(n) = col(7)
       sig_wx(n) = col(8); sig_wy(n) = col(9); sig_wz(n) = col(10)
    End Do
    Close(unit_in)

    If ( sig_y(1) > y_global(1)+1d-10 .Or. sig_y(n_sigma) < y_global(ny_global)-1d-10 ) Then
       If ( myid == 0 ) Write(*,'(A)') &
            ' WARNING: sem_sigma_file does not cover the full domain height -- clamping at the endpoints'
    End If

    If ( sem_eddy_placement == 1 .And. myid == 0 ) &
         Write(*,'(A)') ' Using PDF-weighted eddy placement (ESEM Section 4.2)'

  End Subroutine read_sigma_profile

  ! Generic piecewise-linear interpolation, clamped at the endpoints.
  Real(Int64) Function linterp(xarr, yarr, n, xq) Result(val)
    !$acc routine seq

    Real(Int64),    Intent(In) :: xarr(:), yarr(:), xq
    Integer(Int32), Intent(In) :: n
    Integer(Int32) :: m

    If ( xq <= xarr(1) ) Then
       val = yarr(1)
       Return
    Else If ( xq >= xarr(n) ) Then
       val = yarr(n)
       Return
    End If
    Do m = 1, n-1
       If ( xq >= xarr(m) .And. xq <= xarr(m+1) ) Then
          val = yarr(m) + ( yarr(m+1)-yarr(m) ) * ( xq-xarr(m) ) / ( xarr(m+1)-xarr(m) )
          Return
       End If
    End Do
    val = yarr(n)   ! unreachable, keeps the function well-defined

  End Function linterp

  Real(Int64) Function interp_profile(col, yc) Result(val)
    !$acc routine seq
    Real(Int64), Intent(In) :: col(:), yc
    val = linterp(prof_y, col, n_profile, yc)
  End Function interp_profile

  Real(Int64) Function interp_sigma(col, yc) Result(val)
    !$acc routine seq
    Real(Int64), Intent(In) :: col(:), yc
    val = linterp(sig_y, col, n_sigma, yc)
  End Function interp_sigma

  !> Build each eddy's fixed size class (9-component length scales, Section 4.1) and x-phase/period; (y,z) and sign are re-drawn every cycle by eddy_realize
  Subroutine place_eddies

    Integer(Int32) :: kk
    Real   (Int64) :: r1, r3, y0
    Real   (Int64) :: smaxy_glob, smaxz_glob
    Real   (Int64) :: u_tau, wall_dist, yplus, damp
    Real   (Int64), Parameter :: sem_wall_damping_floor = 0.05d0   ! numerical floor, keeps eddy_sig away from zero

    If ( n_sigma == 0 ) Then
       smaxy_glob = sem_length_scale
       smaxz_glob = sem_length_scale
    Else
       smaxy_glob = Max( Maxval(sig_uy), Maxval(sig_vy), Maxval(sig_wy) )
       smaxz_glob = Max( Maxval(sig_uz), Maxval(sig_vz), Maxval(sig_wz) )
    End If

    box_y_lo = y_global(1)         - smaxy_glob
    box_y_hi = y_global(ny_global) + smaxy_glob
    box_z_lo = z_global(1)         - smaxz_glob
    box_z_hi = z_global(nz_global) + smaxz_glob

    Allocate( eddy_x0(sem_n_eddies), eddy_y0(sem_n_eddies) )
    Allocate( eddy_sig(3,3,sem_n_eddies), eddy_smax(3,sem_n_eddies), eddy_Tperiod(sem_n_eddies) )

    ! identical seed/sequence on every rank: eddy population is a single shared virtual field upstream of the (MPI-decomposed-in-z) inflow plane
    Call random_seed_from(sem_seed)

    If ( sem_eddy_placement == 1 .And. n_sigma > 0 ) Call build_placement_cdf(box_y_lo,box_y_hi)

    If ( sem_wall_damping == 1 ) u_tau = estimate_u_tau()

    Do kk = 1, sem_n_eddies

       ! ---- initial placement height, used only to fix this eddy's size
       ! class (uniform or PDF-weighted, Section 4.2) ----
       Call random_number(r1)
       If ( sem_eddy_placement == 1 .And. n_sigma > 0 ) Then
          y0 = invert_cdf(r1)
       Else
          y0 = box_y_lo + (box_y_hi-box_y_lo)*r1
       End If
       eddy_y0(kk) = y0

       ! ---- this eddy's own 9-component length scales (Section 4.1);
       ! homogeneous mode: all nine = sem_length_scale
       If ( n_sigma == 0 ) Then
          eddy_sig(:,:,kk) = sem_length_scale
       Else
          eddy_sig(1,1,kk) = interp_sigma(sig_ux,y0)
          eddy_sig(1,2,kk) = interp_sigma(sig_uy,y0)
          eddy_sig(1,3,kk) = interp_sigma(sig_uz,y0)
          eddy_sig(2,1,kk) = interp_sigma(sig_vx,y0)
          eddy_sig(2,2,kk) = interp_sigma(sig_vy,y0)
          eddy_sig(2,3,kk) = interp_sigma(sig_vz,y0)
          eddy_sig(3,1,kk) = interp_sigma(sig_wx,y0)
          eddy_sig(3,2,kk) = interp_sigma(sig_wy,y0)
          eddy_sig(3,3,kk) = interp_sigma(sig_wz,y0)
       End If

       ! Van Driest-style near-wall damping: shrinks this eddy's whole size class (preserving its shape/anisotropy) based on its placement height's distance to the nearest wall, so injected turbulence intensity (set independently by prof_R11/22/33 in sem_fluctuation) stays on target while oversized near-wall eddies stop biasing the mean flow
       If ( sem_wall_damping == 1 ) Then
          wall_dist = Max( Min( y0-y_global(1), y_global(ny_global)-y0 ), 0d0 )
          yplus     = wall_dist * u_tau / nu
          damp      = Max( 1d0 - Exp(-yplus/sem_wall_damping_Aplus), sem_wall_damping_floor )
          eddy_sig(:,:,kk) = eddy_sig(:,:,kk) * damp
       End If

       ! bounding scale per direction (Section 4.3): the largest of the
       ! three components' length scales in that direction
       eddy_smax(1,kk) = Max( eddy_sig(1,1,kk), eddy_sig(2,1,kk), eddy_sig(3,1,kk) )  ! x
       eddy_smax(2,kk) = Max( eddy_sig(1,2,kk), eddy_sig(2,2,kk), eddy_sig(3,2,kk) )  ! y
       eddy_smax(3,kk) = Max( eddy_sig(1,3,kk), eddy_sig(2,3,kk), eddy_sig(3,3,kk) )  ! z

       Call random_number(r3)
       eddy_x0(kk) = -eddy_smax(1,kk) + 2d0*eddy_smax(1,kk)*r3

       eddy_Tperiod(kk) = 2d0*eddy_smax(1,kk) / Uconv_sem

    End Do

  End Subroutine place_eddies

  !> Deterministic hash-based uniform(0,1), keyed by (seed,k,n,salt) via MurmurHash3's fmix64, so raw_eddy_sum stays a pure, repeatably-callable function with no mutable per-eddy state
  Real(Int64) Function hash_uniform(seed_in, k, n, salt) Result(u)
    !$acc routine seq

    Integer(Int32), Intent(In) :: seed_in, k, n, salt
    Integer(Int64) :: h

    h = Int(seed_in,8)
    h = Ieor( h, Ishft(Int(k,8),    20) )
    h = Ieor( h, Ishft(Int(n,8),    40) )
    h = Ieor( h, Ishft(Int(salt,8),  5) )
    h = Ieor( h, Ishft(h,-33) ); h = h * (-49064778989728563_8)
    h = Ieor( h, Ishft(h,-33) ); h = h * (-4265267296055464877_8)
    h = Ieor( h, Ishft(h,-33) )

    u = Real( Iand(h, Huge(1_8)), 8 ) / Real( Huge(1_8), 8 )   ! clear sign bit -> [0,1]

  End Function hash_uniform

  !> Fresh transverse position and fluctuation sign for eddy k's n-th cycle (uniform or PDF-weighted, see place_eddies/build_placement_cdf)
  Subroutine eddy_realize(k, n, y_r, z_r, eps_r)
    !$acc routine seq

    Integer(Int32), Intent(In)  :: k, n
    Real   (Int64), Intent(Out) :: y_r, z_r, eps_r(3)

    Real(Int64) :: u1, u2, u3, u4, u5, band

    u1 = hash_uniform(sem_seed, k, n, 1)
    If ( sem_eddy_placement == 1 .And. n_sigma > 0 ) Then
       y_r = invert_cdf(u1)
    Else
       y_r = box_y_lo + (box_y_hi-box_y_lo)*u1
    End If

    ! Inhomogeneous mode: keep eddy k statistically anchored to the
    ! wall-normal region its size class (eddy_sig, fixed at place_eddies
    ! from placement height eddy_y0(k)) was calibrated for. Without this,
    ! the draw above is independent of k's own size, so an outer-layer-sized
    ! eddy can be re-realized at the wall (injecting an oversized fluctuation
    ! there and biasing the mean profile) while a wall-sized eddy is wasted
    ! out in the free stream (under-populating outer-layer variance). Also needed under sem_wall_damping alone (n_sigma==0), since place_eddies gives each eddy a per-eddy-damped eddy_sig/eddy_smax there too.
    If ( n_sigma > 0 .Or. sem_wall_damping == 1 ) Then
       band = sem_placement_band_factor * eddy_smax(2,k)
       y_r  = Min( Max( y_r, eddy_y0(k)-band ), eddy_y0(k)+band )
       y_r  = Min( Max( y_r, box_y_lo ), box_y_hi )
    End If

    u2 = hash_uniform(sem_seed, k, n, 2); z_r = box_z_lo + (box_z_hi-box_z_lo)*u2

    u3 = hash_uniform(sem_seed, k, n, 3); eps_r(1) = Merge(1d0,-1d0, u3 >= 0.5d0)
    u4 = hash_uniform(sem_seed, k, n, 4); eps_r(2) = Merge(1d0,-1d0, u4 >= 0.5d0)
    u5 = hash_uniform(sem_seed, k, n, 5); eps_r(3) = Merge(1d0,-1d0, u5 >= 0.5d0)

  End Subroutine eddy_realize

  !> PDF-weighted eddy placement (ESEM Section 4.2, Eq. 9): trapezoidal CDF over sigma_file y-nodes, weighted inversely by eddy bounding volume so smaller eddies get proportionally more placements
  Subroutine build_placement_cdf(y_lo, y_hi)

    Real(Int64), Intent(In) :: y_lo, y_hi

    Real(Int64), Allocatable :: Vvol(:), w(:)
    Real(Int64) :: Vmax, Vmin, total
    Integer(Int32) :: m

    n_cdf = n_sigma + 2
    Allocate( cdf_y(n_cdf), cdf_val(n_cdf), Vvol(n_cdf), w(n_cdf) )

    cdf_y(1)          = y_lo
    cdf_y(2:n_cdf-1)  = sig_y(1:n_sigma)
    cdf_y(n_cdf)      = y_hi

    Do m = 2, n_cdf-1
       Vvol(m) = Max( sig_ux(m-1), sig_vx(m-1), sig_wx(m-1) ) &
               * Max( sig_uy(m-1), sig_vy(m-1), sig_wy(m-1) ) &
               * Max( sig_uz(m-1), sig_vz(m-1), sig_wz(m-1) )
    End Do
    Vvol(1)     = Vvol(2)          ! flat extrapolation into the padding
    Vvol(n_cdf) = Vvol(n_cdf-1)

    Vmax = Maxval(Vvol); Vmin = Minval(Vvol)
    If ( Vmin <= 0d0 ) Stop 'ERROR: sem_sigma_file implies a non-positive eddy volume somewhere'
    w = ( -Vvol + Vmax )/Vmin + 1d0   ! Eq. (9), Schau et al. 2022

    cdf_val(1) = 0d0
    Do m = 2, n_cdf
       cdf_val(m) = cdf_val(m-1) + 0.5d0*( w(m)+w(m-1) )*( cdf_y(m)-cdf_y(m-1) )
    End Do
    total = cdf_val(n_cdf)
    cdf_val = cdf_val / total

    Deallocate(Vvol,w)

  End Subroutine build_placement_cdf

  Real(Int64) Function invert_cdf(u) Result(yq)
    !$acc routine seq

    Real(Int64), Intent(In) :: u
    Integer(Int32) :: m

    If ( u <= cdf_val(1) ) Then
       yq = cdf_y(1)
       Return
    Else If ( u >= cdf_val(n_cdf) ) Then
       yq = cdf_y(n_cdf)
       Return
    End If
    Do m = 1, n_cdf-1
       If ( u >= cdf_val(m) .And. u <= cdf_val(m+1) ) Then
          yq = cdf_y(m) + ( cdf_y(m+1)-cdf_y(m) ) * ( u-cdf_val(m) ) / ( cdf_val(m+1)-cdf_val(m) )
          Return
       End If
    End Do
    yq = cdf_y(n_cdf)

  End Function invert_cdf

  ! Deterministic RNG seed, identical on every MPI rank
  Subroutine random_seed_from(seed_in)

    Integer(Int32), Intent(In) :: seed_in
    Integer, Allocatable :: seed(:)
    Integer :: n

    Call random_seed(size=n)
    Allocate(seed(n))
    seed = seed_in
    Call random_seed(put=seed)
    Deallocate(seed)

  End Subroutine random_seed_from

  !> ESEM core (Section 3): sample the raw eddy-kernel sum over an ensemble window and store its empirical mean/std, replacing analytical normalisation; signal is periodic in timing but not value (eddy_realize redraws each cycle), so sem_ensemble_periods cycles of the slowest eddy are sampled
  Subroutine build_ensemble_normalisation

    Integer(Int32) :: j, k, s
    Real   (Int64) :: T_E, t_s, ustar, vstar, wstar
    Real   (Int64) :: su, su2, sv, sv2, sw, sw2, suv, corr, max_abs_corr
    Real   (Int64) :: Ns_r

    T_E   = Real(sem_ensemble_periods,8) * Maxval(eddy_Tperiod)
    Ns_r  = Real(sem_ensemble_samples,8)

    If ( myid == 0 ) Write(*,'(A,E12.4,A,I6,A)') &
         ' Building ESEM ensemble normalisation (T_E=', T_E, ', ', sem_ensemble_samples, ' samples/point) ...'

    ! ---- U's own grid (yg,zg): needs u* stats only ----
    Allocate( ens_mean_uU(nyg,nzg), ens_std_uU(nyg,nzg) )
    Do k = 1, nzg
       Do j = 1, nyg
          su = 0d0; su2 = 0d0
          Do s = 0, sem_ensemble_samples-1
             t_s = T_E * Real(s,8) / Ns_r
             Call raw_eddy_sum( yg(j), zg(k), t_s, ustar, vstar, wstar )
             su = su + ustar; su2 = su2 + ustar*ustar
          End Do
          ens_mean_uU(j,k) = su / Ns_r
          ens_std_uU(j,k)  = Sqrt( Max( su2/Ns_r - ens_mean_uU(j,k)**2, 0d0 ) )
       End Do
    End Do

    ! ---- V's own grid (y,zg): needs u* and v* stats ----
    Allocate( ens_mean_uV(ny,nzg), ens_std_uV(ny,nzg) )
    Allocate( ens_mean_vV(ny,nzg), ens_std_vV(ny,nzg) )
    max_abs_corr = 0d0
    Do k = 1, nzg
       Do j = 1, ny
          su = 0d0; su2 = 0d0; sv = 0d0; sv2 = 0d0; suv = 0d0
          Do s = 0, sem_ensemble_samples-1
             t_s = T_E * Real(s,8) / Ns_r
             Call raw_eddy_sum( y(j), zg(k), t_s, ustar, vstar, wstar )
             su = su + ustar; su2 = su2 + ustar*ustar
             sv = sv + vstar; sv2 = sv2 + vstar*vstar
             suv = suv + ustar*vstar
          End Do
          ens_mean_uV(j,k) = su / Ns_r
          ens_std_uV(j,k)  = Sqrt( Max( su2/Ns_r - ens_mean_uV(j,k)**2, 0d0 ) )
          ens_mean_vV(j,k) = sv / Ns_r
          ens_std_vV(j,k)  = Sqrt( Max( sv2/Ns_r - ens_mean_vV(j,k)**2, 0d0 ) )

          ! diagnostic only: sem_fluctuation's Cholesky mixing assumes <u*v*> ~ 0
          If ( ens_std_uV(j,k) > 1d-300 .And. ens_std_vV(j,k) > 1d-300 ) Then
             corr = ( suv/Ns_r - ens_mean_uV(j,k)*ens_mean_vV(j,k) ) / ( ens_std_uV(j,k)*ens_std_vV(j,k) )
             max_abs_corr = Max( max_abs_corr, Abs(corr) )
          End If
       End Do
    End Do
    If ( myid == 0 .And. max_abs_corr > 0.3d0 ) Write(*,'(A,F6.3)') &
         ' WARNING: SEM raw <u*v*> cross-covariance up to ', max_abs_corr

    ! ---- W's own grid (yg,z): needs w* stats only ----
    Allocate( ens_mean_wW(nyg,nz), ens_std_wW(nyg,nz) )
    Do k = 1, nz
       Do j = 1, nyg
          sw = 0d0; sw2 = 0d0
          Do s = 0, sem_ensemble_samples-1
             t_s = T_E * Real(s,8) / Ns_r
             Call raw_eddy_sum( yg(j), z(k), t_s, ustar, vstar, wstar )
             sw = sw + wstar; sw2 = sw2 + wstar*wstar
          End Do
          ens_mean_wW(j,k) = sw / Ns_r
          ens_std_wW(j,k)  = Sqrt( Max( sw2/Ns_r - ens_mean_wW(j,k)**2, 0d0 ) )
       End Do
    End Do

    If ( myid == 0 ) Write(*,'(A)') ' ESEM ensemble normalisation ready.'

  End Subroutine build_ensemble_normalisation

  !> Classical SEM (Jarrin 2006) analytical normalisation
  Subroutine build_analytical_normalisation

    Real(Int64) :: y_lo, y_hi, z_lo, z_hi, VB, recommended_N

    y_lo = y_global(1)         - sem_length_scale
    y_hi = y_global(ny_global) + sem_length_scale
    z_lo = z_global(1)         - sem_length_scale
    z_hi = z_global(nz_global) + sem_length_scale

    VB       = (2d0*sem_length_scale) * (y_hi-y_lo) * (z_hi-z_lo)
    sem_norm = Sqrt( VB )   ! tent_kernel already carries its own sigma^-3/2 factor

    ! Diagnostic: Jarrin's density guideline wants ~1 eddy per sigma^3 cell of the box
    recommended_N = VB / sem_length_scale**3
    If ( Real(sem_n_eddies,8) < 0.05d0*recommended_N ) Then
       If ( myid == 0 ) Then
          Write(*,'(A)')       ' WARNING: classical SEM (sem_use_esem=0): sem_n_eddies is far below the'
          Write(*,'(A)')       '          eddy density this sem_length_scale needs for well-behaved statistics.'
          Write(*,'(A,E12.4)') '          recommended sem_n_eddies (>=)  : ', recommended_N
          Write(*,'(A,I8)')    '          current   sem_n_eddies         : ', sem_n_eddies
          Write(*,'(A)')       '          Increase sem_length_scale/sem_n_eddies, or switch to sem_use_esem=1.'
       End If
    End If

  End Subroutine build_analytical_normalisation

  !> Raw eddy-kernel sum at (yc,zc,t_now), Eq. (1) of Jarrin (2006)/Schau et al. (2022); homogeneous mode uses the tent kernel, inhomogeneous mode uses the unified (cosine) shape function (Section 4.3)
  Subroutine raw_eddy_sum(yc, zc, t_now, ustar, vstar, wstar)
    !$acc routine seq

    Real(Int64), Intent(In)  :: yc, zc, t_now
    Real(Int64), Intent(Out) :: ustar, vstar, wstar

    Real(Int64) :: xk, fx, fy, fz, fprod, phase, halfwidth
    Real(Int64) :: y_r, z_r, eps_r(3)
    Integer(Int32) :: kk, n_cycle

    ustar = 0d0; vstar = 0d0; wstar = 0d0

    If ( sem_divergence_free == 1 ) Then
       Call raw_eddy_curl_sum(yc, zc, t_now, ustar, vstar, wstar)
       Return
    End If

    If ( n_sigma == 0 .And. sem_wall_damping == 0 ) Then

       ! ---- homogeneous: original Jarrin tent-kernel formula, with a
       ! freshly re-randomized (y,z,sign) every cycle (see eddy_realize) ----
       ! (sem_wall_damping alone routes through the inhomogeneous branch below instead, since this branch ignores the per-eddy damped eddy_sig/eddy_smax and would otherwise make sem_wall_damping a no-op)
       halfwidth = sem_length_scale
       Do kk = 1, sem_n_eddies
          phase   = (eddy_x0(kk)+halfwidth) + Uconv_sem*t_now
          n_cycle = Int( Floor( phase / (2d0*halfwidth) ), Int32 )
          xk = -halfwidth + Modulo( phase, 2d0*halfwidth )
          fx = tent_kernel( xk, halfwidth )
          If ( fx == 0d0 ) Cycle

          Call eddy_realize(kk, n_cycle, y_r, z_r, eps_r)

          fy = tent_kernel( yc-y_r, halfwidth )
          If ( fy == 0d0 ) Cycle
          fz = tent_kernel( zc-z_r, halfwidth )
          If ( fz == 0d0 ) Cycle
          fprod = fx*fy*fz
          ustar = ustar + eps_r(1)*fprod
          vstar = vstar + eps_r(2)*fprod
          wstar = wstar + eps_r(3)*fprod
       End Do

    Else

       ! inhomogeneous: unified shape function, per-component/per-direction length scales (ESEM Sections 4.1/4.3), fresh (y,z,sign) every cycle
       Do kk = 1, sem_n_eddies
          phase   = (eddy_x0(kk)+eddy_smax(1,kk)) + Uconv_sem*t_now
          n_cycle = Int( Floor( phase / (2d0*eddy_smax(1,kk)) ), Int32 )
          xk = -eddy_smax(1,kk) + Modulo( phase, 2d0*eddy_smax(1,kk) )
          If ( Abs(xk) >= eddy_smax(1,kk) ) Cycle

          Call eddy_realize(kk, n_cycle, y_r, z_r, eps_r)

          If ( Abs(yc-y_r) >= eddy_smax(2,kk) ) Cycle
          If ( Abs(zc-z_r) >= eddy_smax(3,kk) ) Cycle

          fx = unified_kernel( xk,     eddy_smax(1,kk), eddy_sig(1,1,kk) )
          fy = unified_kernel( yc-y_r, eddy_smax(2,kk), eddy_sig(1,2,kk) )
          fz = unified_kernel( zc-z_r, eddy_smax(3,kk), eddy_sig(1,3,kk) )
          ustar = ustar + eps_r(1)*fx*fy*fz

          fx = unified_kernel( xk,     eddy_smax(1,kk), eddy_sig(2,1,kk) )
          fy = unified_kernel( yc-y_r, eddy_smax(2,kk), eddy_sig(2,2,kk) )
          fz = unified_kernel( zc-z_r, eddy_smax(3,kk), eddy_sig(2,3,kk) )
          vstar = vstar + eps_r(2)*fx*fy*fz

          fx = unified_kernel( xk,     eddy_smax(1,kk), eddy_sig(3,1,kk) )
          fy = unified_kernel( yc-y_r, eddy_smax(2,kk), eddy_sig(3,2,kk) )
          fz = unified_kernel( zc-z_r, eddy_smax(3,kk), eddy_sig(3,3,kk) )
          wstar = wstar + eps_r(3)*fx*fy*fz
       End Do

    End If

  End Subroutine raw_eddy_sum

  !> Divergence-free (Poletto, Craft & Revell 2013) raw eddy-kernel sum: curl of a separable vector potential built from the same eddy_sig/eddy_realize machinery as raw_eddy_sum
  Subroutine raw_eddy_curl_sum(yc, zc, t_now, ustar, vstar, wstar)
    !$acc routine seq

    Real(Int64), Intent(In)  :: yc, zc, t_now
    Real(Int64), Intent(Out) :: ustar, vstar, wstar

    Real(Int64) :: xk, phase, halfwidth, y_r, z_r, eps_r(3)
    Real(Int64) :: fx1,fy1,fz1,dfx1,dfy1,dfz1
    Real(Int64) :: fx2,fy2,fz2,dfx2,dfy2,dfz2
    Real(Int64) :: fx3,fy3,fz3,dfx3,dfy3,dfz3
    Integer(Int32) :: kk, n_cycle

    ustar = 0d0; vstar = 0d0; wstar = 0d0

    ! xk represents x_eddy_k(t)-x (it grows with t as the eddy convects downstream), so d(xk)/dx=-1 a.e. and dPsi/dx = -eps*f'(xk)*f*f (verified: d(xk)/dx=+1 would break div-free)
    If ( n_sigma == 0 ) Then

       ! homogeneous: single sigma per eddy, shared by all three Psi components (tent kernel)
       halfwidth = sem_length_scale
       Do kk = 1, sem_n_eddies
          phase   = (eddy_x0(kk)+halfwidth) + Uconv_sem*t_now
          n_cycle = Int( Floor( phase / (2d0*halfwidth) ), Int32 )
          xk = -halfwidth + Modulo( phase, 2d0*halfwidth )
          If ( Abs(xk) >= halfwidth ) Cycle

          Call eddy_realize(kk, n_cycle, y_r, z_r, eps_r)
          If ( Abs(yc-y_r) >= halfwidth ) Cycle
          If ( Abs(zc-z_r) >= halfwidth ) Cycle

          fx1 = tent_kernel(xk,halfwidth);      dfx1 = tent_kernel_deriv(xk,halfwidth)
          fy1 = tent_kernel(yc-y_r,halfwidth);   dfy1 = tent_kernel_deriv(yc-y_r,halfwidth)
          fz1 = tent_kernel(zc-z_r,halfwidth);   dfz1 = tent_kernel_deriv(zc-z_r,halfwidth)
          fx2=fx1; fy2=fy1; fz2=fz1; dfx2=dfx1; dfy2=dfy1; dfz2=dfz1
          fx3=fx1; fy3=fy1; fz3=fz1; dfx3=dfx1; dfy3=dfy1; dfz3=dfz1

          ustar = ustar + ( eps_r(3)*fx3*dfy3*fz3 - eps_r(2)*fx2*fy2*dfz2 )
          vstar = vstar + ( eps_r(1)*fx1*fy1*dfz1 + eps_r(3)*dfx3*fy3*fz3 )
          wstar = wstar + (-eps_r(2)*dfx2*fy2*fz2 - eps_r(1)*fx1*dfy1*fz1 )
       End Do

    Else

       ! inhomogeneous: unified kernel, per-component-per-direction sigma_ij
       Do kk = 1, sem_n_eddies
          phase   = (eddy_x0(kk)+eddy_smax(1,kk)) + Uconv_sem*t_now
          n_cycle = Int( Floor( phase / (2d0*eddy_smax(1,kk)) ), Int32 )
          xk = -eddy_smax(1,kk) + Modulo( phase, 2d0*eddy_smax(1,kk) )
          If ( Abs(xk) >= eddy_smax(1,kk) ) Cycle

          Call eddy_realize(kk, n_cycle, y_r, z_r, eps_r)
          If ( Abs(yc-y_r) >= eddy_smax(2,kk) ) Cycle
          If ( Abs(zc-z_r) >= eddy_smax(3,kk) ) Cycle

          fx1  = unified_kernel(xk,eddy_smax(1,kk),eddy_sig(1,1,kk))
          dfx1 = unified_kernel_deriv(xk,eddy_smax(1,kk),eddy_sig(1,1,kk))
          fy1  = unified_kernel(yc-y_r,eddy_smax(2,kk),eddy_sig(1,2,kk))
          dfy1 = unified_kernel_deriv(yc-y_r,eddy_smax(2,kk),eddy_sig(1,2,kk))
          fz1  = unified_kernel(zc-z_r,eddy_smax(3,kk),eddy_sig(1,3,kk))
          dfz1 = unified_kernel_deriv(zc-z_r,eddy_smax(3,kk),eddy_sig(1,3,kk))

          fx2  = unified_kernel(xk,eddy_smax(1,kk),eddy_sig(2,1,kk))
          dfx2 = unified_kernel_deriv(xk,eddy_smax(1,kk),eddy_sig(2,1,kk))
          fy2  = unified_kernel(yc-y_r,eddy_smax(2,kk),eddy_sig(2,2,kk))
          dfy2 = unified_kernel_deriv(yc-y_r,eddy_smax(2,kk),eddy_sig(2,2,kk))
          fz2  = unified_kernel(zc-z_r,eddy_smax(3,kk),eddy_sig(2,3,kk))
          dfz2 = unified_kernel_deriv(zc-z_r,eddy_smax(3,kk),eddy_sig(2,3,kk))

          fx3  = unified_kernel(xk,eddy_smax(1,kk),eddy_sig(3,1,kk))
          dfx3 = unified_kernel_deriv(xk,eddy_smax(1,kk),eddy_sig(3,1,kk))
          fy3  = unified_kernel(yc-y_r,eddy_smax(2,kk),eddy_sig(3,2,kk))
          dfy3 = unified_kernel_deriv(yc-y_r,eddy_smax(2,kk),eddy_sig(3,2,kk))
          fz3  = unified_kernel(zc-z_r,eddy_smax(3,kk),eddy_sig(3,3,kk))
          dfz3 = unified_kernel_deriv(zc-z_r,eddy_smax(3,kk),eddy_sig(3,3,kk))

          ustar = ustar + ( eps_r(3)*fx3*dfy3*fz3 - eps_r(2)*fx2*fy2*dfz2 )
          vstar = vstar + ( eps_r(1)*fx1*fy1*dfz1 + eps_r(3)*dfx3*fy3*fz3 )
          wstar = wstar + (-eps_r(2)*dfx2*fy2*fz2 - eps_r(1)*fx1*dfy1*fz1 )
       End Do

    End If

  End Subroutine raw_eddy_curl_sum

  !> Target mean streamwise inflow profile: inflow_type==0 constant (inflow_Uconst), inflow_type==1 interpolated reference profile U(y)
  Real(Int64) Function mean_profile_U(yc) Result(Uc)
    !$acc routine seq

    Real(Int64), Intent(In) :: yc

    If ( inflow_type == 0 ) Then
       Uc = inflow_Uconst
    Else
       Uc = interp_profile(prof_U, yc)
    End If

  End Function mean_profile_U

  !> Target mean inflow temperature profile: n_profile_T==0 falls back to T_ref (uniform), else interpolated prof_T(y)
  Real(Int64) Function mean_profile_T(yc) Result(Tc)
    !$acc routine seq

    Real(Int64), Intent(In) :: yc

    If ( n_profile_T > 0 ) Then
       Tc = linterp(prof_y_T, prof_T, n_profile_T, yc)
    Else
       Tc = T_ref
    End If

  End Function mean_profile_T

  !> Inflow velocity fluctuation at (yc,zc) for component comp (1=U,2=V,3=W): inflow_type==0 none, inflow_type==1 ESEM raw eddy sum normalised (Section 3) then Cholesky-correlated to the target Reynolds-stress tensor (R13=R23=0 assumed)
  Subroutine sem_fluctuation(comp, j, k, yc, zc, t_now, up, vp, wp)
    !$acc routine seq

    Integer(Int32), Intent(In)  :: comp, j, k
    Real   (Int64), Intent(In)  :: yc, zc, t_now
    Real   (Int64), Intent(Out) :: up, vp, wp

    Real(Int64) :: ustar, vstar, wstar, ui, vi, wi
    Real(Int64) :: R11, R22, R33, R12, a11, a21, a22, a33, taper

    If ( inflow_type == 0 ) Then
       up = 0d0; vp = 0d0; wp = 0d0
       Return
    End If

    Call raw_eddy_sum(yc, zc, t_now, ustar, vstar, wstar)

    ui = 0d0; vi = 0d0; wi = 0d0
    If ( sem_use_esem == 1 ) Then
       ! ESEM: exact empirical per-grid-point normalisation
       If ( comp == 1 ) Then
          If ( ens_std_uU(j,k) > 1d-300 ) ui = (ustar - ens_mean_uU(j,k)) / ens_std_uU(j,k)
       Else If ( comp == 2 ) Then
          If ( ens_std_uV(j,k) > 1d-300 ) ui = (ustar - ens_mean_uV(j,k)) / ens_std_uV(j,k)
          If ( ens_std_vV(j,k) > 1d-300 ) vi = (vstar - ens_mean_vV(j,k)) / ens_std_vV(j,k)
       Else
          If ( ens_std_wW(j,k) > 1d-300 ) wi = (wstar - ens_mean_wW(j,k)) / ens_std_wW(j,k)
       End If
    Else
       ! classical SEM: single analytical normalisation factor, applied
       ! uniformly (no per-grid-point lookup, no ensemble precompute)
       ui = ustar * sem_norm / Sqrt( Real(sem_n_eddies,8) )
       vi = vstar * sem_norm / Sqrt( Real(sem_n_eddies,8) )
       wi = wstar * sem_norm / Sqrt( Real(sem_n_eddies,8) )
    End If

    ! correlate the (now exactly zero-mean, unit-variance) fluctuation
    ! with the local target Reynolds-stress tensor via Cholesky
    R11 = interp_profile(prof_R11,yc); R22 = interp_profile(prof_R22,yc)
    R33 = interp_profile(prof_R33,yc); R12 = interp_profile(prof_R12,yc)

    ! Force no-slip on the injected fluctuation at literal solid walls,
    ! independent of whether the user-supplied Reynolds-stress profile
    ! itself decays to zero there: linterp clamps to the profile's
    ! endpoint value outside its data range, so a profile that doesn't
    ! extend flush to the wall would otherwise hold a non-decaying,
    ! nonzero R_ii right up to the wall face.
    taper = 1d0
    If ( wall_active_lo == 1 ) taper = taper * wall_taper( yc-wall_y_lo, wall_Ltaper_lo )
    If ( wall_active_hi == 1 ) taper = taper * wall_taper( wall_y_hi-yc, wall_Ltaper_hi )
    R11 = R11*taper; R22 = R22*taper; R33 = R33*taper; R12 = R12*taper

    a11 = Sqrt( Max(R11,0d0) )
    a21 = Merge( R12/a11, 0d0, a11 > 0d0 )
    a22 = Sqrt( Max(R22 - a21*a21, 0d0) )
    a33 = Sqrt( Max(R33,0d0) )

    up = a11*ui
    vp = a21*ui + a22*vi
    wp = a33*wi

    ! safety-net clip at a few std-devs of the local target stress; rarely triggers with exact ESEM normalisation, kept as defense-in-depth against an under-resolved ensemble window
    up = Sign( Min(Abs(up), sem_clip_sigma*a11), up )
    vp = Sign( Min(Abs(vp), sem_clip_sigma*Sqrt(Max(R22,0d0))), vp )
    wp = Sign( Min(Abs(wp), sem_clip_sigma*a33), wp )

  End Subroutine sem_fluctuation

  !> Smoothstep ramp from 0 (at the wall, d=0) to 1 (by d>=Ltaper); used to enforce no-slip on injected SEM fluctuations at literal solid walls regardless of the input Reynolds-stress profile's near-wall coverage
  Real(Int64) Function wall_taper(d, Ltaper) Result(f)
    !$acc routine seq

    Real(Int64), Intent(In) :: d, Ltaper
    Real(Int64) :: tt

    tt = Min( Max( d/Ltaper, 0d0 ), 1d0 )
    f  = tt*tt*(3d0 - 2d0*tt)

  End Function wall_taper

  !> 1-D tent (triangular) SEM shape function (Jarrin 2006), normalised so integral of f^2 over its support equals 1; homogeneous mode only
  Real(Int64) Function tent_kernel(r, sigma) Result(f)
    !$acc routine seq

    Real(Int64), Intent(In) :: r, sigma

    If ( Abs(r) >= sigma ) Then
       f = 0d0
    Else
       f = Sqrt(1.5d0/sigma) * ( 1d0 - Abs(r)/sigma )
    End If

  End Function tent_kernel

  !> d(tent_kernel)/dr; sign(0)=+1 convention at the r=0 kink (measure-zero)
  Real(Int64) Function tent_kernel_deriv(r, sigma) Result(df)
    !$acc routine seq

    Real(Int64), Intent(In) :: r, sigma

    If ( Abs(r) >= sigma ) Then
       df = 0d0
    Else
       df = -Sqrt(1.5d0/sigma) / sigma * Sign(1d0, r)
    End If

  End Function tent_kernel_deriv

  !> "Unified" inhomogeneous shape function (ESEM Section 4.3, Eq. 11): bounded by sigma_max in this direction, with an extra cosine modulation at sigma_this if smaller, giving every eddy one bounding volume per direction instead of nine
  Real(Int64) Function unified_kernel(r, sigma_max, sigma_this) Result(f)
    !$acc routine seq

    Real(Int64), Intent(In) :: r, sigma_max, sigma_this

    If ( Abs(r) >= sigma_max ) Then
       f = 0d0
    Else If ( sigma_this >= sigma_max - 1d-14 ) Then
       f = Cos( 0.5d0*pi*r/sigma_max )
    Else
       f = Cos( 0.5d0*pi*r/sigma_max ) * Cos( 0.5d0*pi*r/sigma_this )
    End If

  End Function unified_kernel

  !> d(unified_kernel)/dr; same zero-outside-sigma_max convention as unified_kernel at the support boundary (measure-zero)
  Real(Int64) Function unified_kernel_deriv(r, sigma_max, sigma_this) Result(df)
    !$acc routine seq

    Real(Int64), Intent(In) :: r, sigma_max, sigma_this
    Real(Int64) :: a, b

    If ( Abs(r) >= sigma_max ) Then
       df = 0d0
    Else If ( sigma_this >= sigma_max - 1d-14 ) Then
       a  = 0.5d0*pi/sigma_max
       df = -a * Sin(a*r)
    Else
       a  = 0.5d0*pi/sigma_max; b = 0.5d0*pi/sigma_this
       df = -a*Sin(a*r)*Cos(b*r) - b*Cos(a*r)*Sin(b*r)
    End If

  End Function unified_kernel_deriv

  !> One-time TI-rescale setup: resolve the sampling x-station, snapshot the immutable target profile, allocate accumulators, auto-tune ti_rescale_freq/ti_rescale_nstart if requested. No-op unless ti_rescale_active==1.
  Subroutine init_ti_rescale

    Integer(Int32) :: i, best
    Real(Int64) :: d, d_best
    Logical :: auto_freq, auto_nstart
    Real(Int64) :: dt_est, T_eddy, u_tau, t_adv, t_nstart_phys
    Real(Int64), Parameter :: N_eff_target = 100d0

    If ( ti_rescale_active /= 1 ) Return

    If ( x_bc_type /= 1 .Or. inflow_type /= 1 ) Stop &
         'ERROR: ti_rescale_active=1 requires x_bc_type=1 and inflow_type=1 (SEM inflow)'

    ! nearest global interior cc x-index to ti_rescale_x; duplicates probe_output.f90's nearest_cc_global rather than Use it (compiles after sem.f90 / unlinked in standalone sem test targets)
    best = 1
    d_best = Abs(xg_global(2) - ti_rescale_x)
    Do i = 2, nxg_global-1
       d = Abs(xg_global(i) - ti_rescale_x)
       If ( d < d_best ) Then;  d_best = d;  best = i - 1;  End If
    End Do
    ti_i0 = Max(1, Min(nxm_global, best))

    ! auto-tune ti_rescale_freq (<=0) / ti_rescale_nstart (<0) from Uconv_sem, sem_length_scale, and the profile's wall shear
    auto_freq   = ( ti_rescale_freq   <= 0 )
    auto_nstart = ( ti_rescale_nstart < 0 )
    If ( auto_freq .Or. auto_nstart ) Then
       dt_est = Min( Max( cfl_target*streamwise_grid_spacing()/Uconv_sem, dt_min ), dt_max )
       T_eddy = 2d0*sem_length_scale / Uconv_sem
    End If

    If ( auto_freq ) Then
       ti_rescale_freq = Max( 1, Nint( N_eff_target*T_eddy/dt_est ) )
       If ( myid == 0 ) Write(*,'(A,I10)') &
            ' INFO: auto-tuned ti_rescale_freq (ti_rescale_freq<=0 requested) = ', ti_rescale_freq
    End If

    If ( auto_nstart ) Then
       u_tau = estimate_u_tau()
       t_adv              = ti_rescale_x / Uconv_sem
       t_nstart_phys      = Max( 3d0*t_adv, 15d0*(Ly/2d0)/Max(u_tau,1d-12) )
       ti_rescale_nstart  = Max( 0, Nint( t_nstart_phys/dt_est ) )
       If ( myid == 0 ) Write(*,'(A,I10)') &
            ' INFO: auto-tuned ti_rescale_nstart (ti_rescale_nstart<0 requested) = ', ti_rescale_nstart
    End If

    ! heuristic placement checks: no adaptation-length model exists (SEM literature calibrates this empirically)
    If ( ti_rescale_x < 10d0*sem_length_scale ) Then
       If ( myid == 0 ) Write(*,'(A,F10.4,A)') &
            ' WARNING: ti_rescale_x is within 10*sem_length_scale of the inflow (', ti_rescale_x, &
            ' m) -- SEM turbulence may not have adapted to the target statistics yet'
    End If
    If ( (Lx - ti_rescale_x) < 2d0*sem_length_scale ) Then
       If ( myid == 0 ) Write(*,'(A,F10.4,A)') &
            ' WARNING: ti_rescale_x is within 2*sem_length_scale of the outflow (', ti_rescale_x, &
            ' m) -- convective outflow BC may contaminate the sampled statistics'
    End If

    Allocate( prof_R11_target(n_profile), prof_R22_target(n_profile), prof_R33_target(n_profile) )
    prof_R11_target = prof_R11;  prof_R22_target = prof_R22;  prof_R33_target = prof_R33

    Allocate( prof_U_target(n_profile) )
    prof_U_target = prof_U
    Allocate( ti_U_filt(n_profile) )

    Allocate( ti_yg_cc(nym_global) )
    ti_yg_cc = yg_global(2:nyg_global-1)

    Allocate( acc_ti_U(nym_global), acc_ti_V(nym_global), acc_ti_W(nym_global) )
    Allocate( acc_ti_UU(nym_global), acc_ti_VV(nym_global), acc_ti_WW(nym_global), acc_ti_n(nym_global) )
    acc_ti_U = 0d0;  acc_ti_V = 0d0;  acc_ti_W = 0d0
    acc_ti_UU = 0d0; acc_ti_VV = 0d0; acc_ti_WW = 0d0; acc_ti_n = 0d0

    Allocate( ti_R11_filt(n_profile), ti_R22_filt(n_profile), ti_R33_filt(n_profile) )

    Allocate( ti_R_decay_k(n_profile), ti_R_decaying(n_profile) )
    ti_R_decay_k = 0;  ti_R_decaying = .False.
    Allocate( ti_U_decay_k(n_profile), ti_U_decaying(n_profile) )
    ti_U_decay_k = 0;  ti_U_decaying = .False.

    ! resume controller state (live profile, EMA filters, gain-schedule counter) from a prior run; no-op if the file is missing
    If ( restart == 1 ) Call read_ti_rescale_restart

    If ( myid == 0 ) Write(*,'(A,F10.4,A,I6)') &
         ' TI-rescale: sampling station x = ', xg_global(ti_i0+1), ', global cc index = ', ti_i0

  End Subroutine init_ti_rescale

  !> Per-step raw-moment accumulation at the TI-rescale sampling station, summed over the local z range (no MPI; cheap host loop, mirrors accumulate_rsb's style). No-op unless ti_rescale_active==1, istep>=ti_rescale_nstart, and this rank's x-row owns the sampling station.
  Subroutine accumulate_ti_rescale

    Integer(Int32) :: j, k, ia, jg
    Real(Int64) :: uc, vc, wc

    If ( ti_rescale_active /= 1 ) Return
    If ( istep < ti_rescale_nstart ) Return
    ! ti_i0 is a global interior cc x-index; skip on ranks whose x-row doesn't own it
    If ( ti_i0 < ig1_global(myid) .Or. ti_i0 > ig2_global(myid)-2 ) Return

    ia = ti_i0 - ig1_global(myid) + 2   ! global cc index -> local ghost-array index

    Do k = 2, nzg-1
       Do j = 2, nyg-1
          If ( ibm_input_mode >= 1 ) Then
             If ( phi(ia,j,k) < 0d0 ) Cycle
          End If
          jg = j - 1
          uc = 0.5d0*( U(ia-1,j,k) + U(ia,j,k) )
          vc = 0.5d0*( V(ia,j-1,k) + V(ia,j,k) )
          wc = 0.5d0*( W(ia,j,k-1) + W(ia,j,k) )
          acc_ti_U(jg)  = acc_ti_U(jg)  + uc
          acc_ti_V(jg)  = acc_ti_V(jg)  + vc
          acc_ti_W(jg)  = acc_ti_W(jg)  + wc
          acc_ti_UU(jg) = acc_ti_UU(jg) + uc*uc
          acc_ti_VV(jg) = acc_ti_VV(jg) + vc*vc
          acc_ti_WW(jg) = acc_ti_WW(jg) + wc*wc
          acc_ti_n(jg)  = acc_ti_n(jg)  + 1d0
       End Do
    End Do

  End Subroutine accumulate_ti_rescale

  !> Periodic (every ti_rescale_freq steps) reduce + EMA-filtered, deadbanded, anti-windup-clamped multiplicative nudge of the injected SEM target profile toward the wind-tunnel target (per-y-point full gain until that point first lands inside the deadband, then per-point Robbins-Monro gain decay), plus an additive companion nudge of the mean profile prof_U when ti_rescale_u_active==1. No-op unless ti_rescale_active==1.
  Subroutine apply_ti_rescale

    Real(Int64) :: g_U(nym_global), g_V(nym_global), g_W(nym_global)
    Real(Int64) :: g_UU(nym_global), g_VV(nym_global), g_WW(nym_global), g_n(nym_global)
    Real(Int64) :: mean_U(nym_global), mean_V(nym_global), mean_W(nym_global)
    Real(Int64) :: var_UU(nym_global), var_VV(nym_global), var_WW(nym_global)
    Real(Int64) :: R11m, R22m, R33m, ratio11, ratio22, ratio33, dev11, dev22, dev33
    Real(Int64) :: ratio11_raw, ratio22_raw, ratio33_raw
    Real(Int64) :: clip_lo, clip_hi, max_dev, relax_k, max_dev_gain
    Real(Int64) :: Um, bias, bias_raw, bias_clip, relax_k_u, max_dev_u, max_dev_u_gain
    Character(3) :: max_dev_comp
    Integer(Int32) :: max_dev_m, max_dev_u_m
    Logical :: max_dev_sat, max_dev_u_sat
    Integer(Int32) :: j, m

    If ( ti_rescale_active /= 1 ) Return

    Call MPI_Allreduce(acc_ti_U,  g_U,  nym_global, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    Call MPI_Allreduce(acc_ti_V,  g_V,  nym_global, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    Call MPI_Allreduce(acc_ti_W,  g_W,  nym_global, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    Call MPI_Allreduce(acc_ti_UU, g_UU, nym_global, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    Call MPI_Allreduce(acc_ti_VV, g_VV, nym_global, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    Call MPI_Allreduce(acc_ti_WW, g_WW, nym_global, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    Call MPI_Allreduce(acc_ti_n,  g_n,  nym_global, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

    Do j = 1, nym_global
       If ( g_n(j) > 0d0 ) Then
          mean_U(j) = g_U(j)/g_n(j);  mean_V(j) = g_V(j)/g_n(j);  mean_W(j) = g_W(j)/g_n(j)
          var_UU(j) = Max( g_UU(j)/g_n(j) - mean_U(j)**2, 0d0 )
          var_VV(j) = Max( g_VV(j)/g_n(j) - mean_V(j)**2, 0d0 )
          var_WW(j) = Max( g_WW(j)/g_n(j) - mean_W(j)**2, 0d0 )
       Else
          var_UU(j) = 0d0;  var_VV(j) = 0d0;  var_WW(j) = 0d0
       End If
    End Do

    clip_lo = 1d0/ti_rescale_clip;  clip_hi = ti_rescale_clip
    max_dev = 0d0;  max_dev_comp = '---';  max_dev_m = 0;  max_dev_sat = .False.;  max_dev_gain = 0d0
    Do m = 1, n_profile
       ! Below the first resolved LES cell centre, linterp just clamps to that cell's value -- comparing it against a target sampled at (or inside) the wall-modelled sublayer isn't a like-for-like measurement, so skip it rather than let it dominate the diagnostics/correction with an unresolvable mismatch
       If ( prof_y(m) < ti_yg_cc(1) ) Cycle

       ! full undecayed gain at this y until its own deviation first falls inside the deadband (real, deterministic mismatch during the initial transient, not sampling noise); gated per-point rather than profile-wide, since one persistently-unreachable y (e.g. anti-windup-clamped) must not hold every other point at undecayed gain forever
       If ( ti_R_decaying(m) ) Then
          relax_k = Max( ti_rescale_relax / Sqrt(Real(ti_R_decay_k(m)+1,8)), ti_rescale_relax_min )
       Else
          relax_k = ti_rescale_relax
       End If

       R11m = linterp( ti_yg_cc, var_UU, nym_global, prof_y(m) )
       R22m = linterp( ti_yg_cc, var_VV, nym_global, prof_y(m) )
       R33m = linterp( ti_yg_cc, var_WW, nym_global, prof_y(m) )

       ! EMA low-pass filter of the window's measured variance, to reject sampling noise instead of reacting to each noisy window
       If ( .Not. ti_filt_seeded ) Then
          ti_R11_filt(m) = R11m;  ti_R22_filt(m) = R22m;  ti_R33_filt(m) = R33m
       Else
          ti_R11_filt(m) = ti_rescale_filter_alpha*R11m + (1d0-ti_rescale_filter_alpha)*ti_R11_filt(m)
          ti_R22_filt(m) = ti_rescale_filter_alpha*R22m + (1d0-ti_rescale_filter_alpha)*ti_R22_filt(m)
          ti_R33_filt(m) = ti_rescale_filter_alpha*R33m + (1d0-ti_rescale_filter_alpha)*ti_R33_filt(m)
       End If

       ! unclipped ratios, kept only for the true-deviation diagnostic below -- the clipped ratio11/22/33 still drive the actual profile update
       ratio11_raw = prof_R11_target(m) / Max(ti_R11_filt(m),1d-12)
       ratio22_raw = prof_R22_target(m) / Max(ti_R22_filt(m),1d-12)
       ratio33_raw = prof_R33_target(m) / Max(ti_R33_filt(m),1d-12)

       ratio11 = Min( Max( ratio11_raw, clip_lo ), clip_hi )
       ratio22 = Min( Max( ratio22_raw, clip_lo ), clip_hi )
       ratio33 = Min( Max( ratio33_raw, clip_lo ), clip_hi )

       ! deadband: don't perturb the profile chasing a ratio that's within noise of unity
       If ( Abs(ratio11-1d0) < ti_rescale_deadband ) ratio11 = 1d0
       If ( Abs(ratio22-1d0) < ti_rescale_deadband ) ratio22 = 1d0
       If ( Abs(ratio33-1d0) < ti_rescale_deadband ) ratio33 = 1d0

       prof_R11(m) = prof_R11(m) * ratio11**relax_k
       prof_R22(m) = prof_R22(m) * ratio22**relax_k
       prof_R33(m) = prof_R33(m) * ratio33**relax_k

       ! anti-windup: clamp cumulative drift against prof_R**_target directly, since ti_rescale_clip alone only bounds one window's step
       prof_R11(m) = Min( Max( prof_R11(m), prof_R11_target(m)/ti_rescale_abs_clip ), prof_R11_target(m)*ti_rescale_abs_clip )
       prof_R22(m) = Min( Max( prof_R22(m), prof_R22_target(m)/ti_rescale_abs_clip ), prof_R22_target(m)*ti_rescale_abs_clip )
       prof_R33(m) = Min( Max( prof_R33(m), prof_R33_target(m)/ti_rescale_abs_clip ), prof_R33_target(m)*ti_rescale_abs_clip )

       ! track the worst-offending (y, component) pair and whether it's already pinned at the anti-windup bound, for the saturation diagnostic printed below; deviation uses the unclipped ratio so a saturated window still reports its true error instead of the clip bound
       dev11 = Abs(ratio11_raw-1d0);  dev22 = Abs(ratio22_raw-1d0);  dev33 = Abs(ratio33_raw-1d0)
       If ( dev11 >= max_dev ) Then
          max_dev = dev11;  max_dev_comp = 'R11';  max_dev_m = m;  max_dev_gain = relax_k
          max_dev_sat = ( prof_R11(m) <= (1d0+1d-6)*prof_R11_target(m)/ti_rescale_abs_clip .Or. &
                           prof_R11(m) >= (1d0-1d-6)*prof_R11_target(m)*ti_rescale_abs_clip )
       End If
       If ( dev22 >= max_dev ) Then
          max_dev = dev22;  max_dev_comp = 'R22';  max_dev_m = m;  max_dev_gain = relax_k
          max_dev_sat = ( prof_R22(m) <= (1d0+1d-6)*prof_R22_target(m)/ti_rescale_abs_clip .Or. &
                           prof_R22(m) >= (1d0-1d-6)*prof_R22_target(m)*ti_rescale_abs_clip )
       End If
       If ( dev33 >= max_dev ) Then
          max_dev = dev33;  max_dev_comp = 'R33';  max_dev_m = m;  max_dev_gain = relax_k
          max_dev_sat = ( prof_R33(m) <= (1d0+1d-6)*prof_R33_target(m)/ti_rescale_abs_clip .Or. &
                           prof_R33(m) >= (1d0-1d-6)*prof_R33_target(m)*ti_rescale_abs_clip )
       End If

       ! per-point convergence gate: this y's own deviation (not the profile-wide worst case) decides whether its gain starts decaying
       If ( .Not. ti_R_decaying(m) ) Then
          If ( Max(dev11,dev22,dev33) < ti_rescale_deadband ) Then
             ti_R_decaying(m) = .True.
             ti_R_decay_k(m) = -1
          End If
       End If
       ti_R_decay_k(m) = ti_R_decay_k(m) + 1
    End Do

    !$acc update device(prof_R11, prof_R22, prof_R33)

    ! mean-profile companion nudge: additive (not multiplicative -- U crosses/approaches zero at the wall, where a ratio-based law is ill-conditioned), same EMA filter and per-point Robbins-Monro gain schedule as the variance loop above
    max_dev_u = 0d0;  max_dev_u_m = 0;  max_dev_u_sat = .False.;  max_dev_u_gain = 0d0
    If ( ti_rescale_u_active == 1 ) Then
       bias_clip = ti_rescale_u_clip * Uconv_sem

       Do m = 1, n_profile
          ! same rationale as the R_ii loop above: skip points below the first resolved LES cell centre
          If ( prof_y(m) < ti_yg_cc(1) ) Cycle

          ! full undecayed gain at this y until its own bias first falls inside the deadband, gated per-point for the same reason as the variance loop above
          If ( ti_U_decaying(m) ) Then
             relax_k_u = Max( ti_rescale_u_relax / Sqrt(Real(ti_U_decay_k(m)+1,8)), ti_rescale_u_relax_min )
          Else
             relax_k_u = ti_rescale_u_relax
          End If

          Um = linterp( ti_yg_cc, mean_U, nym_global, prof_y(m) )

          If ( .Not. ti_filt_seeded ) Then
             ti_U_filt(m) = Um
          Else
             ti_U_filt(m) = ti_rescale_filter_alpha*Um + (1d0-ti_rescale_filter_alpha)*ti_U_filt(m)
          End If

          bias_raw = prof_U_target(m) - ti_U_filt(m)
          bias = Min( Max(bias_raw, -bias_clip), bias_clip )

          ! deadband: don't perturb the profile chasing a bias that's within noise of zero
          If ( Abs(bias)/Max(Uconv_sem,1d-12) < ti_rescale_u_deadband ) bias = 0d0

          prof_U(m) = prof_U(m) + relax_k_u*bias

          ! anti-windup: clamp cumulative drift against prof_U_target directly, mirrors ti_rescale_abs_clip above
          prof_U(m) = Min( Max( prof_U(m), prof_U_target(m)-ti_rescale_u_abs_clip*Uconv_sem ), &
                                            prof_U_target(m)+ti_rescale_u_abs_clip*Uconv_sem )

          ! diagnostic uses the unclipped bias so a saturated window still reports its true error instead of the clip bound
          If ( Abs(bias_raw)/Max(Uconv_sem,1d-12) >= max_dev_u ) Then
             max_dev_u = Abs(bias_raw)/Max(Uconv_sem,1d-12);  max_dev_u_m = m;  max_dev_u_gain = relax_k_u
             max_dev_u_sat = ( prof_U(m) <= prof_U_target(m)-(1d0-1d-6)*ti_rescale_u_abs_clip*Uconv_sem .Or. &
                                prof_U(m) >= prof_U_target(m)+(1d0-1d-6)*ti_rescale_u_abs_clip*Uconv_sem )
          End If

          ! per-point convergence gate, mirrors the variance loop above
          If ( .Not. ti_U_decaying(m) ) Then
             If ( Abs(bias_raw)/Max(Uconv_sem,1d-12) < ti_rescale_u_deadband ) Then
                ti_U_decaying(m) = .True.
                ti_U_decay_k(m) = -1
             End If
          End If
          ti_U_decay_k(m) = ti_U_decay_k(m) + 1
       End Do

       !$acc update device(prof_U)
    End If

    acc_ti_U = 0d0;  acc_ti_V = 0d0;  acc_ti_W = 0d0
    acc_ti_UU = 0d0; acc_ti_VV = 0d0; acc_ti_WW = 0d0; acc_ti_n = 0d0

    ti_filt_seeded = .True.

    ! persist controller state after every update so a restart can resume from here
    Call write_ti_rescale_restart

    If ( myid == 0 ) Then
       Write(*,'(A,I10,A,F8.4,A,F6.3,A,F8.4)') &
            ' TI-rescale: nudged inflow profile at istep = ', istep, ', max |ratio-1| = ', max_dev, &
            ', gain at worst y = ', max_dev_gain, ', max |U-bias|/Uconv = ', max_dev_u
       Write(*,'(A,A,A,F10.4,A,L1)') &
            '   -> R_ii driver: ', max_dev_comp, ' at y = ', prof_y(max_dev_m), ', pinned at ti_rescale_abs_clip = ', max_dev_sat
       If ( ti_rescale_u_active == 1 ) Write(*,'(A,F10.4,A,F6.3,A,L1)') &
            '   -> U driver: y = ', prof_y(max_dev_u_m), ', gain at worst y = ', max_dev_u_gain, &
            ', pinned at ti_rescale_u_abs_clip = ', max_dev_u_sat
    End If

  End Subroutine apply_ti_rescale

  !> Writes the live TI-rescale controller state (nudged profile, EMA filters, gain-schedule counter) to fields/ti_rescale_data.dat; called after every apply_ti_rescale update so the file always reflects the latest state. Host-only; overwrites the previous snapshot.
  Subroutine write_ti_rescale_restart

    Integer(Int32) :: unit_out

    If ( myid /= 0 ) Return

    Call system('mkdir -p fields/')

    Open( newunit=unit_out, file=Trim(ti_rescale_restart_file), access='stream', &
          form='unformatted', status='replace', action='write' )
    Write(unit_out) ti_rescale_restart_magic
    Write(unit_out) ti_rescale_restart_version
    Write(unit_out) n_profile
    Write(unit_out) ti_filt_seeded
    Write(unit_out) ti_R_decay_k, ti_R_decaying
    Write(unit_out) ti_U_decay_k, ti_U_decaying
    Write(unit_out) prof_R11, prof_R22, prof_R33, prof_U
    Write(unit_out) ti_R11_filt, ti_R22_filt, ti_R33_filt, ti_U_filt
    Close(unit_out)

  End Subroutine write_ti_rescale_restart

  !> Reads controller state written by write_ti_rescale_restart, so a restart resumes the TI-rescale adaptation (nudged profile, EMA filters, per-y gain-decay counters) instead of restarting it from scratch; no-op (fresh start) if the file is missing, and falls back to restarting the per-y gain decay fresh (while keeping the restored profile/filters) for files predating the version-4 per-point layout. Called from init_ti_rescale when restart==1.
  Subroutine read_ti_rescale_restart

    Integer(Int32) :: unit_in, n_profile_f, tag, version
    Logical :: file_exists

    If ( myid == 0 ) Inquire( file=Trim(ti_rescale_restart_file), exist=file_exists )
    Call MPI_Bcast( file_exists, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr )

    If ( .Not. file_exists ) Then
       If ( myid == 0 ) Write(*,'(A)') &
            ' TI-rescale restart: no '//Trim(ti_rescale_restart_file)//' found; starting fresh.'
       Return
    End If

    If ( myid == 0 ) Then
       Open( newunit=unit_in, file=Trim(ti_rescale_restart_file), access='stream', &
             form='unformatted', status='old', action='read' )
       Read(unit_in) tag
       If ( tag == ti_rescale_restart_magic ) Then
          Read(unit_in) version
          Read(unit_in) n_profile_f
       Else
          ! pre-magic file: the leading word we already consumed as "tag" is actually n_profile
          n_profile_f = tag
          version = 1
          Write(*,'(A)') ' TI-rescale restart: pre-versioning file detected, reading as version 1 (no gain-decay-started flag).'
       End If
       If ( n_profile_f /= n_profile ) Stop &
            'ERROR: fields/ti_rescale_data.dat n_profile mismatch with current inflow profile'

       If ( version >= 4 ) Then
          Read(unit_in) ti_filt_seeded
          Read(unit_in) ti_R_decay_k, ti_R_decaying
          Read(unit_in) ti_U_decay_k, ti_U_decaying
       Else
          ! pre-version-4 files carried a single profile-wide gain-decay flag/counter that doesn't map onto the per-y state below -- consume and discard those legacy fields (present unconditionally in v1-3), restart every point's gain decay fresh instead
          Block
            Integer(Int32) :: legacy_k
            Logical :: legacy_gain_decay_started
            Read(unit_in) legacy_k
            If ( version >= 2 ) Read(unit_in) legacy_gain_decay_started
          End Block
          If ( version >= 3 ) Then
             Read(unit_in) ti_filt_seeded
          Else
             ti_filt_seeded = .True.   ! filters restored below already hold real history from the prior run
          End If
          ti_R_decay_k = 0;  ti_R_decaying = .False.
          ti_U_decay_k = 0;  ti_U_decaying = .False.
       End If

       Read(unit_in) prof_R11, prof_R22, prof_R33, prof_U
       Read(unit_in) ti_R11_filt, ti_R22_filt, ti_R33_filt, ti_U_filt
       Close(unit_in)
       Write(*,'(A,I8,A)') ' TI-rescale restart: resumed controller state (format version ', &
            version, ').'
    End If

    Call MPI_Bcast( ti_filt_seeded, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr )
    Call MPI_Bcast( ti_R_decay_k,  n_profile, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
    Call MPI_Bcast( ti_R_decaying, n_profile, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr )
    Call MPI_Bcast( ti_U_decay_k,  n_profile, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
    Call MPI_Bcast( ti_U_decaying, n_profile, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr )
    Call MPI_Bcast( prof_R11,     n_profile, MPI_REAL8,   0, MPI_COMM_WORLD, ierr )
    Call MPI_Bcast( prof_R22,     n_profile, MPI_REAL8,   0, MPI_COMM_WORLD, ierr )
    Call MPI_Bcast( prof_R33,     n_profile, MPI_REAL8,   0, MPI_COMM_WORLD, ierr )
    Call MPI_Bcast( prof_U,       n_profile, MPI_REAL8,   0, MPI_COMM_WORLD, ierr )
    Call MPI_Bcast( ti_R11_filt,  n_profile, MPI_REAL8,   0, MPI_COMM_WORLD, ierr )
    Call MPI_Bcast( ti_R22_filt,  n_profile, MPI_REAL8,   0, MPI_COMM_WORLD, ierr )
    Call MPI_Bcast( ti_R33_filt,  n_profile, MPI_REAL8,   0, MPI_COMM_WORLD, ierr )
    Call MPI_Bcast( ti_U_filt,    n_profile, MPI_REAL8,   0, MPI_COMM_WORLD, ierr )

    !$acc update device(prof_R11, prof_R22, prof_R33, prof_U)

    ! flow field is already turbulent/adapted on restart -- skip the SEM-eddy advection warm-up wait so accumulation resumes immediately
    ti_rescale_nstart = 0

  End Subroutine read_ti_rescale_restart

End Module synthetic_eddy_method
