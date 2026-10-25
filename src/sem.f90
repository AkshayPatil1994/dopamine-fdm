!> Streamwise inflow BC (x_bc_type==1): constant uniform flow, Ensemble SEM, or recycled precursor slice
Module synthetic_eddy_method

  Use iso_fortran_env, Only : Int32, Int64
  Use global
  Use mpi

  Implicit None

  ! ---- mean/Reynolds-stress reference profile (inflow_type==1) ----------
  Integer(Int32) :: n_profile = 0
  Real   (Int64), Allocatable :: prof_y(:), prof_U(:)
  ! donor first-frame z-mean over unique cells (the periodic last cell repeats the first); the wave inlet
  ! subtracts it
  Real   (Int64), Allocatable :: prof_U_fl(:)
  Real   (Int64), Allocatable :: prof_R11(:), prof_R22(:), prof_R33(:), prof_R12(:)
  Logical :: inflow_profile_loaded = .False.   ! guards init_inflow_profile against a second (redundant) file read

  ! optional mean-temperature profile (inflow_temperature_file); n_profile_T==0 -> mean_profile_T falls back
  ! to T_ref
  Integer(Int32) :: n_profile_T = 0
  Real   (Int64), Allocatable :: prof_y_T(:), prof_T(:)
  !$acc declare create(n_profile_T, prof_y_T, prof_T)

  ! eddy length-scale profile sigma_ij(y) (i=u,v,w; j=x,y,z); n_sigma==0 (sem_sigma_file unset) =>
  ! homogeneous mode
  Integer(Int32) :: n_sigma = 0
  Real   (Int64), Allocatable :: sig_y(:)
  Real   (Int64), Allocatable :: sig_ux(:), sig_uy(:), sig_uz(:)
  Real   (Int64), Allocatable :: sig_vx(:), sig_vy(:), sig_vz(:)
  Real   (Int64), Allocatable :: sig_wx(:), sig_wy(:), sig_wz(:)

  ! placement CDF table (Section 4.2), built once if sem_eddy_placement==1
  Integer(Int32) :: n_cdf = 0
  Real   (Int64), Allocatable :: cdf_y(:), cdf_val(:)

  ! eddy population (inflow_type==1): size class and x-phase fixed; (y,z) and sign re-drawn every recycle
  ! (eddy_realize) so the signal does not repeat with the eddy transit period
  Real   (Int64), Allocatable :: eddy_x0(:)
  ! placement height the eddy's size class was calibrated at; anchors eddy_realize's y redraw so a wall-sized
  ! eddy cannot wander into the outer layer and vice versa (inhomogeneous mode)
  Real   (Int64), Allocatable :: eddy_y0(:)
  Real   (Int64), Allocatable :: eddy_sig(:,:,:)   ! (3,3,sem_n_eddies): (i,j,k)=sigma_ij of eddy k
  Real   (Int64), Allocatable :: eddy_smax(:,:)    ! (3,sem_n_eddies): (j,k)=max_i(sigma_ij) of eddy k
  Real   (Int64), Allocatable :: eddy_Tperiod(:)   ! (sem_n_eddies): 2*eddy_smax(1,k)/Uconv_sem
  Real   (Int64) :: Uconv_sem = 0d0                ! eddy-box convection velocity

  ! per-eddy realisation at time sem_cache_t (x position, active flag, re-drawn y/z/sign), shared by all
  ! inlet points of a BC call
  Integer(Int32), Allocatable :: sc_act(:)
  Real   (Int64), Allocatable :: sc_xk(:), sc_y(:), sc_z(:), sc_eps(:,:)
  Real   (Int64) :: sem_cache_t = -Huge(1d0)
  !$acc declare create(sc_act, sc_xk, sc_y, sc_z, sc_eps)

  Real   (Int64), Parameter :: sem_placement_band_factor = 3d0   ! y-band half-width, units of eddy_smax(2,k)

  ! placement box bounds, needed at runtime by eddy_realize
  Real   (Int64) :: box_y_lo, box_y_hi, box_z_lo, box_z_hi

  ! ESEM normalisation (Schau et al. 2022 sec. 3): per-point mean/std of the raw eddy sum, per staggered
  ! component
  Real   (Int64), Allocatable :: ens_mean_uU(:,:), ens_std_uU(:,:)   ! U's grid (nyg,nzg): u* stats
  Real   (Int64), Allocatable :: ens_mean_uV(:,:), ens_std_uV(:,:)   ! V's grid (ny,nzg):  u* stats
  Real   (Int64), Allocatable :: ens_mean_vV(:,:), ens_std_vV(:,:)   ! V's grid (ny,nzg):  v* stats
  Real   (Int64), Allocatable :: ens_mean_wW(:,:), ens_std_wW(:,:)   ! W's grid (nyg,nz):  w* stats

  ! classical (oSEM, Jarrin 2006) analytical normalisation, only for sem_use_esem==0 (homogeneous mode, see
  ! init_inflow)
  Real   (Int64) :: sem_norm = 1d0   ! sqrt(V_B); tent_kernel already carries sigma^-3/2

  Real   (Int64), Parameter :: sem_clip_sigma = 5d0  ! safety-net clip, local target std-devs

  ! near-wall no-slip enforcement (sem_fluctuation): set once in init_inflow so the device routine needs no
  ! y_global; only at Dirichlet walls
  Integer(Int32) :: wall_active_lo = 0, wall_active_hi = 0
  Real   (Int64) :: wall_y_lo = 0d0, wall_y_hi = 0d0, wall_Ltaper_lo = 1d0, wall_Ltaper_hi = 1d0

  ! device-resident state of sem_fluctuation's call chain; set once by init_inflow, read-only afterwards
  !$acc declare create(n_profile, prof_y, prof_U, prof_U_fl, prof_R11, prof_R22, prof_R33, prof_R12)
  !$acc declare create(n_sigma, sig_y, sig_ux, sig_uy, sig_uz, sig_vx, sig_vy, sig_vz, sig_wx, sig_wy, sig_wz)
  !$acc declare create(n_cdf, cdf_y, cdf_val)
  !$acc declare create(eddy_x0, eddy_y0, eddy_sig, eddy_smax, eddy_Tperiod, Uconv_sem)
  !$acc declare create(box_y_lo, box_y_hi, box_z_lo, box_z_hi)
  !$acc declare create(ens_mean_uU, ens_std_uU, ens_mean_uV, ens_std_uV, ens_mean_vV, ens_std_vV, ens_mean_wW, ens_std_wW)
  !$acc declare create(sem_norm)
  !$acc declare create(wall_active_lo, wall_active_hi, wall_y_lo, wall_y_hi, wall_Ltaper_lo, wall_Ltaper_hi)

  ! Recycled precursor inflow (inflow_type==2). rec_active: only the rank owning the x=1 face (row==0); other
  ! ranks no-op. Not rec_unit<0: Open(newunit=) always returns a negative unit.
  Logical        :: rec_active  = .False.
  Integer(Int32) :: rec_unit    = -1
  Integer(Int32) :: rec_ncomp   = 0
  Integer(Int32) :: rec_n1      = 0   ! U/W/T/C's native n1 = nym_global (cell-centre y)
  Integer(Int32) :: rec_n1_v    = 0   ! V's native n1 = ny_global when rec_v_native, else rec_n1
  ! donor has 'n1_V': V on its own y-face grid, no shift; else V is cell-centre-interpolated like U/W and
  ! needs the j-1 shift (see recycle_value)
  Logical        :: rec_v_native = .True.
  Integer(Int32) :: rec_n2_global = 0 ! all components' donor n2 = nzm_global (U/V cell-centre z; W left-face-of-cell z)
  Integer(Int32) :: rec_col_U   = 0, rec_col_V = 0, rec_col_W = 0   ! component present in the donor (U,V,W required)
  Integer(Int32) :: rec_col_T   = 0   ! 0: donor has no T (needed if boussinesq_flag>=1)
  Integer(Int32) :: rec_col_C   = 0   ! 0: donor has no C (needed if sediment_flag>=1)
  Integer(Int32) :: rec_nsnaps  = 0
  Real   (Int64), Allocatable :: rec_times(:)
  ! this rank's donor z range, global 1-based interior cell/face indices (as kg1_global/kg2_global-2)
  Integer(Int32) :: rec_k1 = 0, rec_k2 = 0
  Integer(Int64) :: rec_frame_bytes = 0_Int64
  Integer(Int64) :: rec_off_U = 0_Int64, rec_off_V = 0_Int64, rec_off_W = 0_Int64
  Integer(Int64) :: rec_off_T = 0_Int64, rec_off_C = 0_Int64
  Integer(Int32) :: rec_idx_lo = -1, rec_idx_hi = -1
  Integer(Int32) :: rec_shift_z   = 0   ! donor z-cells; re-drawn at every loop wrap
  Integer(Int32) :: rec_loop_idx  = -999999
  Integer(Int32) :: rec_idx_shift = -999999   ! shift baked into rec_lo_*/rec_hi_*
  ! bracketing frames per component on its own native (n1_c, local nz) grid
  Real   (Int64), Allocatable :: rec_lo_U(:,:), rec_hi_U(:,:)
  Real   (Int64), Allocatable :: rec_lo_V(:,:), rec_hi_V(:,:)
  Real   (Int64), Allocatable :: rec_lo_W(:,:), rec_hi_W(:,:)
  Real   (Int64), Allocatable :: rec_lo_T(:,:), rec_hi_T(:,:)
  Real   (Int64), Allocatable :: rec_lo_C(:,:), rec_hi_C(:,:)
  Real   (Int64) :: rec_frac = 0d0   ! time-interpolation weight between rec_lo_* and rec_hi_*

  !$acc declare create(rec_lo_U, rec_hi_U, rec_lo_V, rec_hi_V, rec_lo_W, rec_hi_W, rec_lo_T, rec_hi_T, rec_lo_C, rec_hi_C, rec_frac, rec_col_U, rec_col_V, rec_col_W, rec_col_T, rec_col_C, rec_n1, rec_n1_v, rec_v_native, rec_k1, rec_k2)

Contains

  !> Abort every rank on a fatal error (a bare Stop exits with status 0 and, on a rank subset, hangs the
  !  others)
  Subroutine sem_abort(msg)
    Character(*), Intent(In) :: msg
    Write(*,'(A)') Trim(msg)
    Call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  End Subroutine sem_abort

  !> Minimum grid spacing (x,y,z faces); floor for the SEM/TI-rescale auto-tuning
  Real(Int64) Function min_grid_spacing() Result(dmin)

    dmin = Min( Minval(x_global(2:nx_global)-x_global(1:nx_global-1)), &
                Minval(y_global(2:ny_global)-y_global(1:ny_global-1)), &
                Minval(z_global(2:nz_global)-z_global(1:nz_global-1)) )

  End Function min_grid_spacing

  !> Per-direction maximum grid spacing, for max_grid_spacing() and enforce_sigma_floor (each sigma_ij is
  !  floored against the spacing of its own direction)
  Subroutine max_grid_spacing_xyz(dx_max, dy_max, dz_max)

    Real(Int64), Intent(Out) :: dx_max, dy_max, dz_max

    dx_max = Maxval(x_global(2:nx_global)-x_global(1:nx_global-1))
    dy_max = Maxval(y_global(2:ny_global)-y_global(1:ny_global-1))
    dz_max = Maxval(z_global(2:nz_global)-z_global(1:nz_global-1))

  End Subroutine max_grid_spacing_xyz

  !> Maximum grid spacing (x,y,z faces): an eddy spans ~2*sem_length_scale, so its floor needs the coarsest
  !  spacing, not min_grid_spacing()
  Real(Int64) Function max_grid_spacing() Result(dmax)

    Real(Int64) :: dx_max, dy_max, dz_max

    Call max_grid_spacing_xyz(dx_max, dy_max, dz_max)
    dmax = Max(dx_max, dy_max, dz_max)

  End Function max_grid_spacing

  !> Friction velocity from dPdx or, if dPdx~0, the mean-profile wall gradient (place_eddies wall damping)
  Real(Int64) Function estimate_u_tau() Result(u_tau)

    Real(Int64) :: dUdy_wall

    If ( Abs(dPdx) > 1d-12 ) Then
       u_tau = Sqrt( Abs(dPdx)*Ly/2d0 )
    Else
       dUdy_wall = ( prof_U(2)-prof_U(1) ) / ( prof_y(2)-prof_y(1) )
       u_tau = Sqrt( nu*Abs(dUdy_wall) )
    End If

  End Function estimate_u_tau

  !> Provisional sem_length_scale from geometry while the profile is parsed (<=0); refined by
  !  refine_sem_length_scale
  Real(Int64) Function geometric_length_scale_estimate() Result(L0)

    L0 = 0.3d0 * Min(Ly, Lz)
    L0 = Max( L0, 5d0*min_grid_spacing() )

  End Function geometric_length_scale_estimate

  !> Refines sem_length_scale from the mean sig_ux, else sqrt(R11)/|dU/dy| (as derive_mixing_length); clamped
  !  to [5*min spacing, 0.3*min(Ly,Lz)]
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

  !> Floor on sem_length_scale: an eddy spans ~2*sem_length_scale, and one narrower than a few cells is
  !  destroyed within the first cells, collapsing the injected stresses. Uses max_grid_spacing() (a
  !  wall-clustered min spacing would be too permissive).
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

  !> Per-direction floor on the sigma_ij profiles (sem_sigma_file or derived), each against the spacing of its
  !  own direction; clamps in place and warns once per profile. No-op when n_sigma==0.
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

  !> Auto-sets sem_n_eddies (<=0): ~1 eddy per sem_length_scale^3 of the padded inflow box (Jarrin); mirrors
  !  build_analytical_normalisation's recommended_N
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

  ! One-time inflow setup (no-op unless x_bc_type==1). Reading the mean profile first lets an IC generator
  ! seed the interior from it; idempotent via inflow_profile_loaded
  Subroutine init_inflow_profile

    If ( x_bc_type /= 1 ) Return
    If ( inflow_profile_loaded ) Return

    If ( inflow_type == 1 ) Then
       If ( sem_profile_format == 0 ) Then
          Call read_mean_profile
       Else If ( sem_profile_format == 1 ) Then
          Call read_mean_profile_TI
       Else
          Call sem_abort( 'ERROR: sem_profile_format must be 0 (Reynolds-stress) or 1 (wind-tunnel TI)' )
       End If
       inflow_profile_loaded = .True.
    Else If ( inflow_type == 2 ) Then
       ! mean profile only (z-average of the donor's first frame), so genGridandIC seeds the interior without
       ! a step-1 divergence spike
       Call read_recycle_mean_profile
       inflow_profile_loaded = .True.
    End If

  End Subroutine init_inflow_profile

  Subroutine init_inflow

    Logical :: auto_length_scale

    If ( x_bc_type /= 1 ) Return

    If ( inflow_type == 2 .Or. ( inflow_type == 3 .And. wave_turb == 1 ) ) Then
       ! the wave inlet takes only the donor fluctuation, so its mean (z-mean of the first frame) is needed
       If ( inflow_type == 3 ) Call read_recycle_mean_profile
       If ( inflow_type == 2 ) Call init_inflow_profile
       Call init_inflow_recycle
       Return
    End If

    If ( inflow_type /= 1 ) Return

    auto_length_scale = ( sem_length_scale <= 0d0 )
    If ( auto_length_scale ) sem_length_scale = geometric_length_scale_estimate()

    Call init_inflow_profile
    Call read_sigma_profile   ! sem_sigma_file overrides derived length scales
    Call enforce_sigma_floor   ! clamp sig_* the grid cannot resolve

    If ( boussinesq_flag >= 1 .And. Len_trim(inflow_temperature_file) > 0 ) Call read_mean_profile_T

    If ( sem_use_esem == 0 .And. ( n_sigma > 0 .Or. sem_wall_damping == 1 ) ) Call sem_abort( &
         'ERROR: sem_use_esem=0 (classical SEM) requires sem_sigma_file to be unset and sem_wall_damping=0 -- ' // &
         'classical SEM has no well-defined normalisation for inhomogeneous length scales' )

    If ( sem_divergence_free == 1 .And. sem_use_esem == 0 ) Call sem_abort( &
         'ERROR: sem_divergence_free=1 requires sem_use_esem=1 -- classical (analytical) ' // &
         'normalisation has no closed form for the curl construction' )

    Uconv_sem = Sum(prof_U) / Real(n_profile,8)
    If ( Uconv_sem <= 0d0 ) Call sem_abort( 'ERROR: SEM inflow profile has non-positive bulk velocity' )

    ! near-wall no-slip setup (wall_taper/sem_fluctuation): only at literal Dirichlet walls, not periodic or
    ! free-slip
    wall_active_lo = Merge( 1, 0, y_bc_type == 1 .And. bc_face_ylo == 1 )
    wall_active_hi = Merge( 1, 0, y_bc_type == 1 .And. bc_face_yhi == 1 )
    wall_y_lo      = y_global(1)
    wall_y_hi      = y_global(ny_global)
    wall_Ltaper_lo = 2d0*( y_global(2)         - y_global(1) )
    wall_Ltaper_hi = 2d0*( y_global(ny_global) - y_global(ny_global-1) )

    If ( auto_length_scale ) Call refine_sem_length_scale
    Call enforce_sem_length_scale_floor   ! also clamps an explicit sem_length_scale
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

  !> Parse a recycled-inflow donor's _meta.txt: (1) dopamine-ESEM per-component grids with an 'n1_V' key (V's
  !  y-face count); (2) generic probe_output slice, V cell-centre interpolated onto the U/W grid, no 'n1_V'
  !  (n1_v = n1, v_native = .False.)
  Subroutine read_recycle_meta(ncomp, n1, n1_v, n2, nsnaps, colU, colV, colW, colT, colC, v_native)

    Integer(Int32), Intent(Out) :: ncomp, n1, n1_v, n2, nsnaps, colU, colV, colW, colT, colC
    Logical,        Intent(Out) :: v_native

    Integer(Int32) :: u_meta, ios, eqpos, ic
    Character(300) :: fname, line
    Character(20)  :: key
    Character(8)   :: comps_str

    If ( Len_trim(inflow_recycle_file) == 0 ) &
         Call sem_abort('ERROR: inflow_type=2 (recycled inflow) requires inflow_recycle_file to be set')

    Write(fname,'(A,A)') Trim(inflow_recycle_file), '_meta.txt'
    Open(newunit=u_meta, file=Trim(fname), form='formatted', status='old', action='read', iostat=ios)
    If ( ios /= 0 ) Call sem_abort('ERROR: cannot open inflow_recycle_file meta (expected <inflow_recycle_file>_meta.txt)')

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
               Call sem_abort('ERROR: inflow_recycle_file meta is not an x-normal slice (dir/=x)')
       End Select
    End Do
    Close(u_meta)

    If ( ncomp < 1 .Or. n1 < 1 .Or. n2 < 1 .Or. nsnaps < 1 ) &
         Call sem_abort('ERROR: inflow_recycle_file meta is incomplete or malformed (expected ncomp/n1/n2/comps/nsnaps, ' // &
         'as written by probe_output.f90''s slice output or dopamine-ESEM)')

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
    If ( colU == 0 .Or. colV == 0 .Or. colW == 0 ) Call sem_abort(&
         'ERROR: inflow_recycle_file must contain U,V,W (set slice_comps="UVW..." on the donor run)')

  End Subroutine read_recycle_meta

  !> Mean-profile seed for inflow_type==2: every rank z-averages U(y) of the donor's first frame; R11..R12
  !  stay zero since sem_fluctuation is not called, the profile only seeds the IC mean in genGridandIC
  Subroutine read_recycle_mean_profile

    Integer(Int32) :: ncomp, n1, n1_v, n2, nsnaps, colU, colV, colW, colT, colC, unit_in, ios, jy, n2_unique
    Logical :: v_native
    Real(Int64), Allocatable :: frame_U(:,:), frame_il(:,:,:)
    Character(300) :: fname

    Call read_recycle_meta(ncomp, n1, n1_v, n2, nsnaps, colU, colV, colW, colT, colC, v_native)

    If ( n1 /= nym_global ) Call sem_abort('ERROR: inflow_recycle_file ny (n1) does not match this run''s nym_global')
    If ( n2 /= nzm_global ) Call sem_abort('ERROR: inflow_recycle_file nz (n2) does not match this run''s nzm_global')

    n_profile = n1
    Allocate( prof_y(n_profile), prof_U(n_profile), prof_U_fl(n_profile) )
    Allocate( prof_R11(n_profile), prof_R22(n_profile), prof_R33(n_profile), prof_R12(n_profile) )
    prof_y = ym_global
    prof_R11 = 0d0;  prof_R22 = 0d0;  prof_R33 = 0d0;  prof_R12 = 0d0

    ! only rank 0 reads the donor's first U block (first on disk); the result is tiny, so it is broadcast
    If ( myid == 0 ) Then
       n2_unique = n2 - Merge(1, 0, z_bc_type == 0)
       Allocate( frame_U(n1,n2) )
       Write(fname,'(A,A)') Trim(inflow_recycle_file), '.bin'
       Open(newunit=unit_in, file=Trim(fname), access='stream', form='unformatted', &
            status='old', action='read', iostat=ios)
       If ( ios /= 0 ) Call sem_abort('ERROR: cannot open inflow_recycle_file data (expected <inflow_recycle_file>.bin)')
       If ( v_native ) Then
          Read(unit_in) frame_U   ! U is always the first block of the first frame
       Else
          Allocate( frame_il(ncomp,n1,n2) )
          Read(unit_in) frame_il
          frame_U = frame_il(1,:,:)
          Deallocate(frame_il)
       End If
       Close(unit_in)
       Do jy = 1, n1
          prof_U(jy) = Sum(frame_U(jy,:)) / Real(n2,8)
          prof_U_fl(jy) = Sum(frame_U(jy,1:n2_unique)) / Real(n2_unique,8)
       End Do
       Deallocate(frame_U)
    End If

    Call Mpi_bcast( prof_U, n_profile, MPI_real8, 0, MPI_COMM_WORLD, ierr )
    Call Mpi_bcast( prof_U_fl, n_profile, MPI_real8, 0, MPI_COMM_WORLD, ierr )
    !$acc update device(prof_U_fl)

  End Subroutine read_recycle_mean_profile

  !> One-time setup of the inflow_type==2 buffers; only the rank owning the x=1 face opens the donor .bin and
  !  primes the rec_lo/rec_hi frames for its z-slab
  Subroutine init_inflow_recycle

    Integer(Int32) :: ios, u_times
    Character(300) :: fname
    ! nvfortran 23.3 parser bug: an inline rec_* integer in Merge(...,rec_col_T==1) gives a bogus 'Syntax
    ! error at or near REC'; precompute the masks into plain Logicals
    Logical :: has_T, has_C

    rec_unit = -1
    rec_active = .False.
    ! this rank owns the x=1 face iff it is in row 0 of the p_row/p_col grid (as decomp's x_periodic_partner;
    ! sem.f90 is also built standalone by tests)
    If ( myid/p_col /= 0 ) Return

    Call read_recycle_meta(rec_ncomp, rec_n1, rec_n1_v, rec_n2_global, rec_nsnaps, &
         rec_col_U, rec_col_V, rec_col_W, rec_col_T, rec_col_C, rec_v_native)

    If ( rec_n1 /= nym_global ) Call sem_abort('ERROR: inflow_recycle_file ny (n1) does not match this run''s nym_global')
    If ( rec_v_native ) Then
       If ( rec_n1_v /= ny_global ) &
            Call sem_abort('ERROR: inflow_recycle_file n1_V (V''s own y-face count) does not match this run''s ny_global')
    End If
    If ( rec_n2_global /= nzm_global ) Call sem_abort('ERROR: inflow_recycle_file nz (n2) does not match this run''s nzm_global')
    If ( rec_nsnaps < 2 ) Call sem_abort('ERROR: inflow_recycle_file must contain at least 2 snapshots to interpolate in time')
    If ( boussinesq_flag >= 1 .And. rec_col_T == 0 ) Call sem_abort('ERROR: boussinesq_flag>=1 with ' // &
         'inflow_type=2 requires T in the donor slice (donor slice_comps must include ''T'')')
    If ( sediment_flag >= 1 .And. rec_col_C == 0 ) Call sem_abort('ERROR: sediment_flag>=1 with ' // &
         'inflow_type=2 requires C in the donor slice (donor slice_comps must include ''C'')')

    Allocate( rec_times(rec_nsnaps) )
    Write(fname,'(A,A)') Trim(inflow_recycle_file), '_times.bin'
    Open(newunit=u_times, file=Trim(fname), access='stream', form='unformatted', &
         status='old', action='read', iostat=ios)
    If ( ios /= 0 ) Call sem_abort('ERROR: cannot open inflow_recycle_file times (expected <inflow_recycle_file>_times.bin)')
    Read(u_times) rec_times
    Close(u_times)

    ! this rank's donor z-slab, same convention as the writer (kg_g = kg1_global(myid)+ka-2)
    rec_k1 = kg1_global(myid)
    rec_k2 = kg2_global(myid) - 2
    !$acc update device(rec_k1, rec_k2)

    has_T = ( rec_col_T == 1 );  has_C = ( rec_col_C == 1 )

    ! blocks in fixed U,V,W,T,C disk order (only present ones occupy space)
    rec_off_U = 0_Int64
    rec_off_V = rec_off_U + 8_Int64*Int(rec_n1,  Int64)*Int(rec_n2_global,Int64)
    rec_off_W = rec_off_V + 8_Int64*Int(rec_n1_v,Int64)*Int(rec_n2_global,Int64)
    rec_off_T = rec_off_W + 8_Int64*Int(rec_n1,  Int64)*Int(rec_n2_global,Int64)
    rec_off_C = rec_off_T + Merge(8_Int64*Int(rec_n1,Int64)*Int(rec_n2_global,Int64), 0_Int64, has_T)
    rec_frame_bytes = rec_off_C + Merge(8_Int64*Int(rec_n1,Int64)*Int(rec_n2_global,Int64), 0_Int64, has_C)
    If ( .Not. rec_v_native ) Then
       ! probe_output.f90 slices are interleaved (ncomp,n1,n2): rec_off_* hold the 0-based component column,
       ! not a byte offset
       rec_off_U = 0_Int64;  rec_off_V = 1_Int64;  rec_off_W = 2_Int64
       rec_off_T = 3_Int64;  rec_off_C = Merge(4_Int64, 3_Int64, has_T)
       rec_frame_bytes = 8_Int64*Int(rec_ncomp,Int64)*Int(rec_n1,Int64)*Int(rec_n2_global,Int64)
    End If

    Write(fname,'(A,A)') Trim(inflow_recycle_file), '.bin'
    Open(newunit=rec_unit, file=Trim(fname), access='stream', form='unformatted', &
         status='old', action='read', iostat=ios)
    If ( ios /= 0 ) Call sem_abort('ERROR: cannot open inflow_recycle_file data (expected <inflow_recycle_file>.bin)')
    rec_active = .True.

    Allocate( rec_lo_U(rec_n1,   rec_k2-rec_k1+1) );  Allocate( rec_hi_U(rec_n1,   rec_k2-rec_k1+1) )
    Allocate( rec_lo_V(rec_n1_v, rec_k2-rec_k1+1) );  Allocate( rec_hi_V(rec_n1_v, rec_k2-rec_k1+1) )
    Allocate( rec_lo_W(rec_n1,   rec_k2-rec_k1+1) );  Allocate( rec_hi_W(rec_n1,   rec_k2-rec_k1+1) )
    Allocate( rec_lo_T(Merge(rec_n1,0,has_T), Merge(rec_k2-rec_k1+1,0,has_T)) )
    Allocate( rec_hi_T(Merge(rec_n1,0,has_T), Merge(rec_k2-rec_k1+1,0,has_T)) )
    Allocate( rec_lo_C(Merge(rec_n1,0,has_C), Merge(rec_k2-rec_k1+1,0,has_C)) )
    Allocate( rec_hi_C(Merge(rec_n1,0,has_C), Merge(rec_k2-rec_k1+1,0,has_C)) )
    rec_idx_lo = -1;  rec_idx_hi = -1

    !$acc update device(rec_col_U, rec_col_V, rec_col_W, rec_col_T, rec_col_C, rec_n1, rec_n1_v, rec_v_native)

    Call update_inflow_recycle(t)   ! prime rec_lo_*/rec_hi_*/rec_frac before the first BC application

    Write(*,'(A,I0,A,A,A,I0,A)') ' Rank ', myid, ': recycled precursor inflow from ', &
         Trim(inflow_recycle_file), ' (', rec_nsnaps, ' donor snapshots)'

  End Subroutine init_inflow_recycle

  !> Read one contiguous donor-z span (z0..z0+count-1, global 1-based, no shift/wrap) of one component's block
  !  into buf(:,zoff:); each block is column-major so a z-span is one contiguous read
  Subroutine read_recycle_span(frame_idx, comp_off, comp_n1, z0, count, buf, zoff)

    Integer(Int32), Intent(In)    :: frame_idx, comp_n1, z0, count, zoff
    Integer(Int64), Intent(In)    :: comp_off
    Real   (Int64), Intent(InOut) :: buf(:,:)

    Integer(Int64) :: byte_pos
    Real   (Int64), Allocatable :: il(:,:,:)

    If ( count < 1 ) Return
    If ( .Not. rec_v_native ) Then
       Allocate( il(rec_ncomp, comp_n1, count) )
       byte_pos = Int(frame_idx-1,Int64)*rec_frame_bytes + &
                  Int(z0-1,Int64)*8_Int64*Int(rec_ncomp,Int64)*Int(comp_n1,Int64) + 1_Int64
       Read(rec_unit, pos=byte_pos) il
       buf(:,zoff:zoff+count-1) = il(comp_off+1,:,:)
       Return
    End If
    byte_pos = Int(frame_idx-1,Int64)*rec_frame_bytes + comp_off + &
               Int(z0-1,Int64)*8_Int64*Int(comp_n1,Int64) + 1_Int64
    Read(rec_unit, pos=byte_pos) buf(:,zoff:zoff+count-1)

  End Subroutine read_recycle_span

  !> Read one component's block of the local z-slab, circularly shifted by shift_z donor cells (periodic z
  !  range [1,rec_n2_global]); split into two reads when it wraps
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

  !> Read all active components' blocks of one donor frame into rec_lo_* (dst_lo) or rec_hi_*
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

  !> Advance the donor bracket (rec_lo/rec_hi/rec_frac) to t_now (rank owning x=1 only). Frames are re-read
  !  only when the bracket or shift changes; when it slides forward the old hi becomes the new lo
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

    ! new spanwise shift per donor pass: hash of (seed,loop_idx), restart-safe, no RNG state to persist
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

    ! binary search of rec_times for the bracketing pair
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
          ! bracket slid forward, same shift: reuse the old hi as the new lo
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

  !> Time-interpolated recycled-inflow value of component comp (1=U,2=V,3=W,4=T,5=C) at local (j,k). Each
  !  component reads its own donor grid, so there is no half-cell shift: U/W use j-1 (y ghost), V uses j
  !  (y-face), U/V use k-1, W uses k (z-face, donor stores left faces).
  Real(Int64) Function recycle_value(comp, j, k) Result(val)
    !$acc routine seq

    Integer(Int32), Intent(In) :: comp, j, k

    Integer(Int32) :: jy, kz, nzloc

    nzloc = rec_k2 - rec_k1 + 1

    Select Case (comp)
    Case (1)   ! U: (yg,zg) ghost cell-centre
       jy = Max( 1, Min(rec_n1, j-1) );  kz = Max( 1, Min(nzloc, k-1) )
       val = rec_lo_U(jy,kz) + rec_frac*( rec_hi_U(jy,kz) - rec_lo_U(jy,kz) )
    Case (2)   ! V: native donor -- y-face (no shift); generic donor -- cell-centre like U/W (j-1)
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
    If ( ios /= 0 ) Call sem_abort( 'ERROR: cannot open inflow_temperature_file' )

    n_profile_T = 0
    Do
       Read(unit_in,'(A)',iostat=ios) line
       If ( ios /= 0 ) Exit
       line = Adjustl(line)
       If ( Len_trim(line) == 0 .Or. line(1:1) == '#' ) Cycle
       n_profile_T = n_profile_T + 1
    End Do
    If ( n_profile_T < 2 ) Call sem_abort( 'ERROR: inflow_temperature_file: fewer than 2 data rows found' )

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
       If ( ios /= 0 ) Call sem_abort( 'ERROR: inflow_temperature_file: failed to parse a data row (need 2 columns: y T)' )
       prof_y_T(n) = col(1)
       prof_T(n)   = col(2)
    End Do
    Close(unit_in)

  End Subroutine read_mean_profile_T

  ! Read the mean-velocity / Reynolds-stress reference profile
  Subroutine read_mean_profile

    Integer(Int32) :: unit_in, ios, n
    Real   (Int64) :: col(8)
    Character(300) :: line

    If ( myid == 0 ) Write(*,'(A,A)') ' Reading SEM inflow profile from ', Trim(inflow_profile_file)

    ! every rank reads the (small, shared) profile file independently; free-form text, skip blank/'#' lines,
    ! each row: y U V W uu vv ww uv [uw vw] (extra cols read but unused, see sem_fluctuation's R13=R23=0
    ! assumption)
    Open(newunit=unit_in, file=Trim(inflow_profile_file), status='old', action='read', iostat=ios)
    If ( ios /= 0 ) Call sem_abort( 'ERROR: cannot open inflow_profile_file' )

    n_profile = 0
    Do
       Read(unit_in,'(A)',iostat=ios) line
       If ( ios /= 0 ) Exit
       line = Adjustl(line)
       If ( Len_trim(line) == 0 .Or. line(1:1) == '#' ) Cycle
       n_profile = n_profile + 1
    End Do
    If ( n_profile < 2 ) Call sem_abort( 'ERROR: inflow_profile_file: fewer than 2 data rows found' )

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
       If ( ios /= 0 ) Call sem_abort( 'ERROR: inflow_profile_file: failed to parse a data row (need >=8 columns)' )
       prof_y(n)   = col(1)   ! y
       prof_U(n)   = col(2)   ! U
       prof_R11(n) = col(5)   ! uu
       prof_R22(n) = col(6)   ! vv
       prof_R33(n) = col(7)   ! ww
       prof_R12(n) = col(8)   ! uv

       ! PSD check (R13=R23=0): R11,R33>=0 and R11*R22-R12^2>=0
       If ( prof_R11(n) < 0d0 .Or. prof_R33(n) < 0d0 .Or. &
            prof_R11(n)*prof_R22(n) - prof_R12(n)**2 < -1d-10*Max(prof_R11(n)*prof_R22(n),1d-30) ) &
            Call sem_abort( 'ERROR: inflow_profile_file: Reynolds-stress tensor is not positive-semidefinite at some y' )
    End Do
    Close(unit_in)

    If ( prof_y(1) > y_global(1)+1d-10 .Or. prof_y(n_profile) < y_global(ny_global)-1d-10 ) Then
       If ( myid == 0 ) Write(*,'(A)') &
            ' WARNING: inflow_profile_file does not cover the full domain height -- clamping at the endpoints'
    End If

  End Subroutine read_mean_profile

  !> Read a wind-tunnel profile (sem_profile_format==1): header "# z U Iu Iv Iw [Lux ... Lwz]"; intensities
  !  become Reynolds stresses (uv=0), length scales optional with the fallback tiers in sem.md
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
    If ( ios /= 0 ) Call sem_abort( 'ERROR: cannot open inflow_profile_file' )

    ! ---- header: first non-blank line, must start with '#' and list column names ----
    Do
       Read(unit_in,'(A)',iostat=ios) line
       If ( ios /= 0 ) Call sem_abort( 'ERROR: inflow_profile_file: missing "# z U Iu Iv Iw ..." header' )
       line = Adjustl(line)
       If ( Len_trim(line) == 0 ) Cycle
       If ( line(1:1) /= '#' ) Call sem_abort( 'ERROR: inflow_profile_file: first line must be "# ..." header' )
       Call parse_TI_header(line(2:), slot_name, n_slot, col_of_slot, ncols)
       Exit
    End Do

    If ( Any( col_of_slot(1:5) == 0 ) ) &
         Call sem_abort( 'ERROR: inflow_profile_file: header must include z (or y), U, Iu, Iv, Iw' )
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
    If ( n_profile < 2 ) Call sem_abort( 'ERROR: inflow_profile_file: fewer than 2 data rows found' )

    Allocate( prof_y(n_profile), prof_U(n_profile) )
    Allocate( prof_R11(n_profile), prof_R22(n_profile), prof_R33(n_profile), prof_R12(n_profile) )

    ! sig_* is filled when length scales are measured/ratio-derived or mixing-length-derived; the latter
    ! needs sem_use_esem==1 (classical SEM cannot normalise inhomogeneous scales)
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
       If ( ios /= 0 ) Call sem_abort( 'ERROR: inflow_profile_file: failed to parse a data row' )

       prof_y(n) = row(col_of_slot(1))
       prof_U(n) = row(col_of_slot(2))
       Iu_v = row(col_of_slot(3)); Iv_v = row(col_of_slot(4)); Iw_v = row(col_of_slot(5))
       If ( Iu_v < 0d0 .Or. Iv_v < 0d0 .Or. Iw_v < 0d0 ) &
            Call sem_abort( 'ERROR: inflow_profile_file: turbulence intensities must be non-negative' )
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

    ! no L** column: derive an inhomogeneous length-scale profile from mixing-length theory
    ! (derive_mixing_length)
    If ( .Not. any_Lscale .And. n_sigma > 0 ) Call derive_mixing_length

    If ( prof_y(1) > y_global(1)+1d-10 .Or. prof_y(n_profile) < y_global(ny_global)-1d-10 ) Then
       If ( myid == 0 ) Write(*,'(A)') &
            ' WARNING: inflow_profile_file does not cover the full domain height -- clamping at the endpoints'
    End If

  End Subroutine read_mean_profile_TI

  !> Mixing-length fallback when no L** column is given: L(y) = sqrt(R11)/|dU/dy| for all components' x-scale,
  !  y/z via sem_Lscale_ratio_y/z; clamped to [sem_length_scale, outer_frac*profile height] because |dU/dy|->0
  !  near the outer edge
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
       If ( islot < 0 ) Call sem_abort( 'ERROR: inflow_profile_file: unrecognized header column name "'//Trim(tok)//'"' )
       If ( col_of_slot(islot) /= 0 ) Call sem_abort( 'ERROR: inflow_profile_file: duplicate header column name' )
       col_of_slot(islot) = ncols
    End Do

  End Subroutine parse_TI_header

  !> Map a header token to its slot id (1..n_slot), 'y' accepted as an alias for slot 1 ('z'); -1 if
  !  unrecognized
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

  !> Fill one component's (Lox,Loy,Loz): measured > ratio*Lox > sem_length_scale (sem.md)
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
         Call sem_abort( 'ERROR: inflow_profile_file: length scales must be strictly positive' )

  End Subroutine fill_length_scale

  !> Read the optional eddy length-scale profile (ESEM 4.1); empty sem_sigma_file => homogeneous (n_sigma=0)
  Subroutine read_sigma_profile

    Integer(Int32) :: unit_in, ios, n
    Real   (Int64) :: col(10)
    Character(300) :: line

    If ( Len_trim(sem_sigma_file) == 0 ) Return   ! keep what read_mean_profile_TI derived

    If ( myid == 0 ) Write(*,'(A,A)') ' Reading inhomogeneous SEM length-scale profile from ', Trim(sem_sigma_file)

    If ( Allocated(sig_y) ) &
         Deallocate( sig_y, sig_ux, sig_uy, sig_uz, sig_vx, sig_vy, sig_vz, sig_wx, sig_wy, sig_wz )

    Open(newunit=unit_in, file=Trim(sem_sigma_file), status='old', action='read', iostat=ios)
    If ( ios /= 0 ) Call sem_abort( 'ERROR: cannot open sem_sigma_file' )

    n_sigma = 0
    Do
       Read(unit_in,'(A)',iostat=ios) line
       If ( ios /= 0 ) Exit
       line = Adjustl(line)
       If ( Len_trim(line) == 0 .Or. line(1:1) == '#' ) Cycle
       n_sigma = n_sigma + 1
    End Do
    If ( n_sigma < 2 ) Call sem_abort( 'ERROR: sem_sigma_file: fewer than 2 data rows found' )

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
       If ( ios /= 0 ) Call sem_abort( 'ERROR: sem_sigma_file: failed to parse a data row (need 10 columns)' )
       If ( Minval(col(2:10)) <= 0d0 ) Call sem_abort( 'ERROR: sem_sigma_file: length scales must be strictly positive' )
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

  ! Piecewise-linear interpolation clamped at the endpoints; xarr non-decreasing, bracket by binary search
  ! (called at every inflow point and RK sub-stage)
  Real(Int64) Function linterp(xarr, yarr, n, xq) Result(val)
    !$acc routine seq

    Real(Int64),    Intent(In) :: xarr(:), yarr(:), xq
    Integer(Int32), Intent(In) :: n
    Integer(Int32) :: m, lo, hi

    If ( xq <= xarr(1) ) Then
       val = yarr(1)
       Return
    Else If ( xq >= xarr(n) ) Then
       val = yarr(n)
       Return
    End If
    lo = 1
    hi = n-1
    Do While ( lo < hi )
       m = (lo+hi)/2
       If ( xarr(m+1) < xq ) Then
          lo = m+1
       Else
          hi = m
       End If
    End Do
    m = lo
    val = yarr(m) + ( yarr(m+1)-yarr(m) ) * ( xq-xarr(m) ) / ( xarr(m+1)-xarr(m) )

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

  !> Build each eddy's fixed size class (ESEM 4.1) and x-phase/period; (y,z) and sign are redrawn per cycle by
  !  eddy_realize
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
    Allocate( sc_act(sem_n_eddies), sc_xk(sem_n_eddies), sc_y(sem_n_eddies), sc_z(sem_n_eddies), sc_eps(3,sem_n_eddies) )
    sem_cache_t = -Huge(1d0)

    ! identical seed on every rank: one shared virtual eddy field upstream of the z-decomposed inflow plane
    Call random_seed_from(sem_seed)

    If ( sem_eddy_placement == 1 .And. n_sigma > 0 ) Call build_placement_cdf(box_y_lo,box_y_hi)

    If ( sem_wall_damping == 1 ) u_tau = estimate_u_tau()

    Do kk = 1, sem_n_eddies

       ! initial placement height, used only to fix this eddy's size class (uniform or PDF-weighted)
       Call random_number(r1)
       If ( sem_eddy_placement == 1 .And. n_sigma > 0 ) Then
          y0 = invert_cdf(r1)
       Else
          y0 = box_y_lo + (box_y_hi-box_y_lo)*r1
       End If
       eddy_y0(kk) = y0

       ! this eddy's 9 length scales (Section 4.1); homogeneous: all = sem_length_scale
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

       ! Van Driest-style damping: shrinks the whole size class (shape kept) by the wall distance of its
       ! placement height, so oversized near-wall eddies stop biasing the mean while the intensity stays set
       ! by prof_R11/22/33
       If ( sem_wall_damping == 1 ) Then
          wall_dist = Max( Min( y0-y_global(1), y_global(ny_global)-y0 ), 0d0 )
          yplus     = wall_dist * u_tau / nu
          damp      = Max( 1d0 - Exp(-yplus/sem_wall_damping_Aplus), sem_wall_damping_floor )
          eddy_sig(:,:,kk) = eddy_sig(:,:,kk) * damp
       End If

       ! bounding scale per direction: largest of the three components' scales (Section 4.3)
       eddy_smax(1,kk) = Max( eddy_sig(1,1,kk), eddy_sig(2,1,kk), eddy_sig(3,1,kk) )  ! x
       eddy_smax(2,kk) = Max( eddy_sig(1,2,kk), eddy_sig(2,2,kk), eddy_sig(3,2,kk) )  ! y
       eddy_smax(3,kk) = Max( eddy_sig(1,3,kk), eddy_sig(2,3,kk), eddy_sig(3,3,kk) )  ! z

       Call random_number(r3)
       eddy_x0(kk) = -eddy_smax(1,kk) + 2d0*eddy_smax(1,kk)*r3

       eddy_Tperiod(kk) = 2d0*eddy_smax(1,kk) / Uconv_sem

    End Do

  End Subroutine place_eddies

  !> Deterministic hash uniform(0,1) of (seed,k,n,salt) via MurmurHash3 fmix64, so raw_eddy_sum stays pure
  !  with no per-eddy state
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

  !> Fresh transverse position and sign for eddy k's n-th cycle (uniform or PDF-weighted,
  !  place_eddies/build_placement_cdf)
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

    ! Inhomogeneous mode (and sem_wall_damping alone): keep eddy k near the wall-normal region its size class
    ! was calibrated for, else an outer-layer eddy lands at the wall (biased mean) and a wall-sized one is
    ! wasted outside
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

  !> PDF-weighted placement (ESEM 4.2, Eq. 9): trapezoidal CDF over sigma_file y-nodes, weighted inversely by
  !  eddy volume so small eddies are placed more often
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
    If ( Vmin <= 0d0 ) Call sem_abort( 'ERROR: sem_sigma_file implies a non-positive eddy volume somewhere' )
    w = ( -Vvol + Vmax )/Vmin + 1d0   ! Eq. (9), Schau et al. 2022

    cdf_val(1) = 0d0
    Do m = 2, n_cdf
       cdf_val(m) = cdf_val(m-1) + 0.5d0*( w(m)+w(m-1) )*( cdf_y(m)-cdf_y(m-1) )
    End Do
    total = cdf_val(n_cdf)
    cdf_val = cdf_val / total

    Deallocate(Vvol,w)

  End Subroutine build_placement_cdf

  ! cdf_val is non-decreasing (trapezoid integral of a positive weight); binary search per eddy per placement
  ! draw
  Real(Int64) Function invert_cdf(u) Result(yq)
    !$acc routine seq

    Real(Int64), Intent(In) :: u
    Integer(Int32) :: m, lo, hi

    If ( u <= cdf_val(1) ) Then
       yq = cdf_y(1)
       Return
    Else If ( u >= cdf_val(n_cdf) ) Then
       yq = cdf_y(n_cdf)
       Return
    End If
    lo = 1
    hi = n_cdf-1
    Do While ( lo < hi )
       m = (lo+hi)/2
       If ( cdf_val(m+1) < u ) Then
          lo = m+1
       Else
          hi = m
       End If
    End Do
    m = lo
    yq = cdf_y(m) + ( cdf_y(m+1)-cdf_y(m) ) * ( u-cdf_val(m) ) / ( cdf_val(m+1)-cdf_val(m) )

  End Function invert_cdf

  ! identical RNG seed on every rank
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

  !> ESEM core (Sec. 3): empirical mean/std of the raw eddy-kernel sum over sem_ensemble_periods cycles of the
  !  slowest eddy (periodic in timing, not value, as eddy_realize redraws each cycle)
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

          ! diagnostic only: Cholesky mixing assumes <u*v*> ~ 0
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

    ! diagnostic: ~1 eddy per sigma^3 of the box (Jarrin)
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

  !> Realise every eddy at t_now once (x position, y/z/sign) so each inlet point only evaluates kernels
  Subroutine sem_cache_refresh(t_now)

    Real(Int64), Intent(In) :: t_now

    Real(Int64) :: halfwidth, phase, xk, fx
    Integer(Int32) :: kk, n_cycle
    Logical :: homog

    If ( t_now == sem_cache_t ) Return
    sem_cache_t = t_now
    homog = ( n_sigma == 0 .And. sem_wall_damping == 0 )
    halfwidth = sem_length_scale
    !$acc parallel loop default(present) private(phase,n_cycle,xk,fx)
    Do kk = 1, sem_n_eddies
       If ( homog ) Then
          phase   = (eddy_x0(kk)+halfwidth) + Uconv_sem*t_now
          n_cycle = Int( Floor( phase / (2d0*halfwidth) ), Int32 )
          xk = -halfwidth + Modulo( phase, 2d0*halfwidth )
          fx = tent_kernel( xk, halfwidth )
          sc_act(kk) = Merge(0, 1, fx == 0d0)
       Else
          phase   = (eddy_x0(kk)+eddy_smax(1,kk)) + Uconv_sem*t_now
          n_cycle = Int( Floor( phase / (2d0*eddy_smax(1,kk)) ), Int32 )
          xk = -eddy_smax(1,kk) + Modulo( phase, 2d0*eddy_smax(1,kk) )
          sc_act(kk) = Merge(0, 1, Abs(xk) >= eddy_smax(1,kk))
       End If
       sc_xk(kk) = xk
       If ( sc_act(kk) == 1 ) Call eddy_realize(kk, n_cycle, sc_y(kk), sc_z(kk), sc_eps(:,kk))
    End Do

  End Subroutine sem_cache_refresh

  !> raw_eddy_sum from the cached realisations, summing only the terms comp needs (V mixes u* and v*); same
  !  order as raw_eddy_sum
  Subroutine raw_eddy_sum_cached(yc, zc, comp, ustar, vstar, wstar)
    !$acc routine seq

    Real(Int64), Intent(In)  :: yc, zc
    Integer(Int32), Intent(In) :: comp
    Real(Int64), Intent(Out) :: ustar, vstar, wstar

    Real(Int64) :: fx, fy, fz, fprod, halfwidth
    Integer(Int32) :: kk
    Logical :: need_u, need_v, need_w

    ustar = 0d0; vstar = 0d0; wstar = 0d0
    need_u = ( comp <= 2 );  need_v = ( comp == 2 );  need_w = ( comp == 3 )

    If ( n_sigma == 0 .And. sem_wall_damping == 0 ) Then
       halfwidth = sem_length_scale
       Do kk = 1, sem_n_eddies
          If ( sc_act(kk) == 0 ) Cycle
          fx = tent_kernel( sc_xk(kk), halfwidth )
          fy = tent_kernel( yc-sc_y(kk), halfwidth )
          If ( fy == 0d0 ) Cycle
          fz = tent_kernel( zc-sc_z(kk), halfwidth )
          If ( fz == 0d0 ) Cycle
          fprod = fx*fy*fz
          If ( need_u ) ustar = ustar + sc_eps(1,kk)*fprod
          If ( need_v ) vstar = vstar + sc_eps(2,kk)*fprod
          If ( need_w ) wstar = wstar + sc_eps(3,kk)*fprod
       End Do
    Else
       Do kk = 1, sem_n_eddies
          If ( sc_act(kk) == 0 ) Cycle
          If ( Abs(yc-sc_y(kk)) >= eddy_smax(2,kk) ) Cycle
          If ( Abs(zc-sc_z(kk)) >= eddy_smax(3,kk) ) Cycle
          If ( need_u ) Then
             fx = unified_kernel( sc_xk(kk),    eddy_smax(1,kk), eddy_sig(1,1,kk) )
             fy = unified_kernel( yc-sc_y(kk), eddy_smax(2,kk), eddy_sig(1,2,kk) )
             fz = unified_kernel( zc-sc_z(kk), eddy_smax(3,kk), eddy_sig(1,3,kk) )
             ustar = ustar + sc_eps(1,kk)*fx*fy*fz
          End If
          If ( need_v ) Then
             fx = unified_kernel( sc_xk(kk),    eddy_smax(1,kk), eddy_sig(2,1,kk) )
             fy = unified_kernel( yc-sc_y(kk), eddy_smax(2,kk), eddy_sig(2,2,kk) )
             fz = unified_kernel( zc-sc_z(kk), eddy_smax(3,kk), eddy_sig(2,3,kk) )
             vstar = vstar + sc_eps(2,kk)*fx*fy*fz
          End If
          If ( need_w ) Then
             fx = unified_kernel( sc_xk(kk),    eddy_smax(1,kk), eddy_sig(3,1,kk) )
             fy = unified_kernel( yc-sc_y(kk), eddy_smax(2,kk), eddy_sig(3,2,kk) )
             fz = unified_kernel( zc-sc_z(kk), eddy_smax(3,kk), eddy_sig(3,3,kk) )
             wstar = wstar + sc_eps(3,kk)*fx*fy*fz
          End If
       End Do
    End If

  End Subroutine raw_eddy_sum_cached

  !> Raw eddy-kernel sum at (yc,zc,t_now), Eq. (1) of Jarrin (2006)/Schau et al. (2022): tent kernel
  !  (homogeneous) or unified cosine shape (inhomogeneous, Sec. 4.3)
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

       ! homogeneous: Jarrin tent kernel, (y,z,sign) redrawn per cycle (eddy_realize). sem_wall_damping alone
       ! uses the inhomogeneous branch below, which honours the per-eddy damped eddy_sig/eddy_smax
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

       ! inhomogeneous: unified shape function, per-component/per-direction length scales (ESEM Sections
       ! 4.1/4.3), fresh (y,z,sign) every cycle
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

  !> Divergence-free (Poletto, Craft & Revell 2013) eddy sum: curl of a separable vector potential from the
  !  same eddy_sig/eddy_realize machinery
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

    ! xk = x_eddy_k(t)-x grows with t, so d(xk)/dx=-1 and dPsi/dx = -eps*f'(xk)*f*f (d(xk)/dx=+1 would break
    ! div-free)
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

  !> Target mean streamwise inflow profile: inflow_type==0 constant (inflow_Uconst), inflow_type==1
  !  interpolated reference profile U(y)
  Real(Int64) Function mean_profile_U(yc) Result(Uc)
    !$acc routine seq

    Real(Int64), Intent(In) :: yc

    If ( inflow_type == 0 ) Then
       Uc = inflow_Uconst
    Else If ( inflow_type == 3 ) Then
       Uc = 0d0
    Else
       Uc = interp_profile(prof_U, yc)
    End If

  End Function mean_profile_U

  !> Target mean inflow temperature profile: n_profile_T==0 falls back to T_ref (uniform), else interpolated
  !  prof_T(y)
  Real(Int64) Function mean_profile_T(yc) Result(Tc)
    !$acc routine seq

    Real(Int64), Intent(In) :: yc

    If ( n_profile_T > 0 ) Then
       Tc = linterp(prof_y_T, prof_T, n_profile_T, yc)
    Else
       Tc = T_ref
    End If

  End Function mean_profile_T

  !> Inflow fluctuation at (yc,zc), component comp: none (inflow_type==0) or ESEM-normalised eddy sum
  !  Cholesky-correlated to the target stresses (R13=R23=0)
  Subroutine sem_fluctuation(comp, j, k, yc, zc, t_now, up, vp, wp, use_cache)
    !$acc routine seq

    Integer(Int32), Intent(In)  :: comp, j, k
    Real   (Int64), Intent(In)  :: yc, zc, t_now
    Real   (Int64), Intent(Out) :: up, vp, wp
    Logical, Intent(In)  :: use_cache   ! use the realisations cached by sem_cache_refresh (not div-free)

    Real(Int64) :: ustar, vstar, wstar, ui, vi, wi
    Real(Int64) :: R11, R22, R33, R12, a11, a21, a22, a33, taper

    If ( inflow_type == 0 ) Then
       up = 0d0; vp = 0d0; wp = 0d0
       Return
    End If

    If ( use_cache ) Then
       Call raw_eddy_sum_cached(yc, zc, comp, ustar, vstar, wstar)
    Else
       Call raw_eddy_sum(yc, zc, t_now, ustar, vstar, wstar)
    End If

    ui = 0d0; vi = 0d0; wi = 0d0
    If ( sem_use_esem == 1 ) Then
       ! ESEM: empirical per-point normalisation
       If ( comp == 1 ) Then
          If ( ens_std_uU(j,k) > 1d-300 ) ui = (ustar - ens_mean_uU(j,k)) / ens_std_uU(j,k)
       Else If ( comp == 2 ) Then
          If ( ens_std_uV(j,k) > 1d-300 ) ui = (ustar - ens_mean_uV(j,k)) / ens_std_uV(j,k)
          If ( ens_std_vV(j,k) > 1d-300 ) vi = (vstar - ens_mean_vV(j,k)) / ens_std_vV(j,k)
       Else
          If ( ens_std_wW(j,k) > 1d-300 ) wi = (wstar - ens_mean_wW(j,k)) / ens_std_wW(j,k)
       End If
    Else
       ! classical SEM: one analytical factor, no per-point lookup
       ui = ustar * sem_norm / Sqrt( Real(sem_n_eddies,8) )
       vi = vstar * sem_norm / Sqrt( Real(sem_n_eddies,8) )
       wi = wstar * sem_norm / Sqrt( Real(sem_n_eddies,8) )
    End If

    ! correlate the zero-mean, unit-variance fluctuation with the target stresses via Cholesky
    R11 = interp_profile(prof_R11,yc); R22 = interp_profile(prof_R22,yc)
    R33 = interp_profile(prof_R33,yc); R12 = interp_profile(prof_R12,yc)

    ! force no-slip at solid walls: linterp clamps to the profile end, which could leave a nonzero R_ii at
    ! the wall face
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

    ! safety-net clip at a few std-devs of the target stress; rarely active with exact ESEM normalisation
    up = Sign( Min(Abs(up), sem_clip_sigma*a11), up )
    vp = Sign( Min(Abs(vp), sem_clip_sigma*Sqrt(Max(R22,0d0))), vp )
    wp = Sign( Min(Abs(wp), sem_clip_sigma*a33), wp )

  End Subroutine sem_fluctuation

  !> Smoothstep ramp 0 (wall) to 1 (d>=Ltaper): no-slip on the injected fluctuation at solid walls whatever
  !  the profile's near-wall coverage
  Real(Int64) Function wall_taper(d, Ltaper) Result(f)
    !$acc routine seq

    Real(Int64), Intent(In) :: d, Ltaper
    Real(Int64) :: tt

    tt = Min( Max( d/Ltaper, 0d0 ), 1d0 )
    f  = tt*tt*(3d0 - 2d0*tt)

  End Function wall_taper

  !> 1-D tent (triangular) SEM shape function (Jarrin 2006), normalised so integral of f^2 over its support
  !  equals 1; homogeneous mode only
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

  !> Unified inhomogeneous shape function (ESEM 4.3, Eq. 11): bounded by sigma_max in the direction, cosine
  !  modulation at sigma_this if smaller; one bounding volume per direction
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

  !> d(unified_kernel)/dr; same zero-outside-sigma_max convention as unified_kernel at the support boundary
  !  (measure-zero)
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
End Module synthetic_eddy_method
