!> Standalone verification driver for src/sem.f90's Bezier-parametrized inflow-optimization
!> helpers (bezier_point, bezier_eval, init_bezier_control_points, build_initial_slopes,
!> scalar_correction_step, secant_update_slopes, inflow_opt_wall_excluded,
!> inflow_opt_max_rel_residual); built via the CMake target verify_inflow_opt (see
!> CMakeLists.txt), not part of the dopamine executable
Program verify_inflow_opt

  Use iso_fortran_env, Only : Int32, Int64
  Use mpi
  Use global
  Use synthetic_eddy_method

  Implicit None

  Real(Int64), Parameter :: tol = 1d-6
  Logical :: all_pass

  Call Mpi_Init(ierr)
  Call Mpi_Comm_rank(MPI_COMM_WORLD, myid, ierr)
  Call Mpi_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

  all_pass = .True.

  Call case_bezier_endpoints
  Call case_bezier_symmetric_hump
  Call case_init_bezier_control_points
  Call case_slopes_and_correction
  Call case_correction_scale
  Call case_secant_update_slopes
  Call case_wall_excluded
  Call case_max_rel_residual

  If ( myid == 0 ) Then
     If ( all_pass ) Then
        Write(*,'(A)') 'verify_inflow_opt: all cases PASSED'
     Else
        Write(*,'(A)') 'verify_inflow_opt: FAILED'
        Call Mpi_Finalize(ierr)
        Stop 1
     End If
  End If

  Call Mpi_Finalize(ierr)

Contains

  Subroutine expect(cond, label)
    Logical,      Intent(In) :: cond
    Character(*), Intent(In) :: label
    If ( myid /= 0 ) Return
    If ( cond ) Then
       Write(*,'(A,A)') '  PASS: ', label
    Else
       Write(*,'(A,A)') '  FAIL: ', label
       all_pass = .False.
    End If
  End Subroutine expect

  !> bezier_eval must reproduce the first/last control-point value exactly at (and beyond) the endpoints
  Subroutine case_bezier_endpoints

    Real(Int64) :: cp_y(4), cp_x(4)

    If ( myid == 0 ) Write(*,'(A)') '=== bezier_eval endpoint clamping ==='
    cp_y = (/ 0d0, 1d0, 2d0, 3d0 /)
    cp_x = (/ 0.2d0, 0.9d0, 0.7d0, 0.1d0 /)

    Call expect( Abs(bezier_eval(cp_y,cp_x,4,0d0)  - 0.2d0) < tol, 'bezier_eval at y=y_cp(1) == cp_x(1)' )
    Call expect( Abs(bezier_eval(cp_y,cp_x,4,3d0)  - 0.1d0) < tol, 'bezier_eval at y=y_cp(n) == cp_x(n)' )
    Call expect( Abs(bezier_eval(cp_y,cp_x,4,-1d0) - 0.2d0) < tol, 'bezier_eval below y_cp(1) clamps to cp_x(1)' )
    Call expect( Abs(bezier_eval(cp_y,cp_x,4,5d0)  - 0.1d0) < tol, 'bezier_eval above y_cp(n) clamps to cp_x(n)' )

  End Subroutine case_bezier_endpoints

  !> With uniformly spaced control-point heights, the curve's y(t) is exactly linear (a Bezier
  !> curve degree-elevates a linear control polygon to itself), so the parameter t at a query
  !> height is known analytically; cross-check bezier_eval against the closed-form cubic
  !> Bernstein-basis evaluation at that t for a symmetric hump control net.
  Subroutine case_bezier_symmetric_hump

    Real(Int64) :: cp_y(4), cp_x(4), t, expected

    If ( myid == 0 ) Write(*,'(A)') '=== bezier_eval matches closed-form cubic Bernstein value ==='
    cp_y = (/ 0d0, 1d0, 2d0, 3d0 /)
    cp_x = (/ 0d0, 1d0, 1d0, 0d0 /)

    t = 0.5d0   ! y(t) = 3t is exactly linear here, so yq=1.5 <-> t=0.5
    expected = (1d0-t)**3*cp_x(1) + 3d0*(1d0-t)**2*t*cp_x(2) &
             + 3d0*(1d0-t)*t**2*cp_x(3) + t**3*cp_x(4)
    Call expect( Abs(bezier_eval(cp_y,cp_x,4,1.5d0) - expected) < tol, &
         'bezier_eval(yq=1.5) == closed-form cubic Bernstein value (0.75)' )
    Call expect( Abs(expected-0.75d0) < tol, 'sanity: closed-form value is 0.75' )

  End Subroutine case_bezier_symmetric_hump

  !> init_bezier_control_points must: pin control-point heights to the profile's own bounds, place
  !> them monotonically increasing and clustered toward the wall, and sample (not fit) the
  !> inflow profile at each height for the control-point target values.
  Subroutine case_init_bezier_control_points

    Integer(Int32) :: i

    If ( myid == 0 ) Write(*,'(A)') '=== init_bezier_control_points control-point placement ==='

    n_profile = 5
    Allocate( prof_y(n_profile) )
    prof_y = (/ 0d0, 0.25d0, 0.5d0, 1d0, 2d0 /)
    Allocate( prof_R11(n_profile), prof_R22(n_profile), prof_R33(n_profile) )
    prof_R11 = (/ 0.9d0, 0.7d0, 0.6d0, 0.4d0, 0.3d0 /)
    prof_R22 = (/ 0.5d0, 0.4d0, 0.3d0, 0.2d0, 0.1d0 /)
    prof_R33 = (/ 0.8d0, 0.6d0, 0.5d0, 0.3d0, 0.2d0 /)

    n_bezier = 5
    Call init_bezier_control_points

    Call expect( Allocated(y_cp) .And. Size(y_cp) == n_bezier, 'y_cp allocated with n_bezier entries' )
    Call expect( Abs(y_cp(1)-prof_y(1)) < tol, 'y_cp(1) == prof_y(1)' )
    Call expect( Abs(y_cp(n_bezier)-prof_y(n_profile)) < tol, 'y_cp(n_bezier) == prof_y(n_profile)' )

    Block
      Logical :: monotonic, clustered_near_wall
      monotonic = .True.
      Do i = 2, n_bezier
         If ( y_cp(i) <= y_cp(i-1) ) monotonic = .False.
      End Do
      Call expect( monotonic, 'y_cp strictly increasing' )
      ! quadratic clustering: the first interior spacing must be smaller than the last
      clustered_near_wall = ( (y_cp(2)-y_cp(1)) < (y_cp(n_bezier)-y_cp(n_bezier-1)) )
      Call expect( clustered_near_wall, 'y_cp clustered toward the wall (first spacing < last spacing)' )
    End Block

    Call expect( Abs(x_cp_R22_target(1)-prof_R22(1)) < tol, &
         'x_cp_R22_target(1) samples prof_R22 at y_cp(1)' )
    Call expect( Abs(x_cp_R22_target(n_bezier)-prof_R22(n_profile)) < tol, &
         'x_cp_R22_target(n_bezier) samples prof_R22 at y_cp(n_bezier)' )
    Call expect( All(x_cp_R22 == x_cp_R22_target), 'x_cp_R22 initialized to the target (no correction yet)' )
    Call expect( All(x_cp_R33 == x_cp_R33_target), 'x_cp_R33 initialized to the target (no correction yet)' )

    Deallocate( prof_y, prof_R11, prof_R22, prof_R33 )
    Deallocate( y_cp, x_cp_R22_target, x_cp_R33_target, u_cp_target, x_cp_R22, x_cp_R33 )

  End Subroutine case_init_bezier_control_points

  !> build_initial_slopes + scalar_correction_step(scale=1) must reproduce an exact Newton step
  !> at an interior control point: doubling v'^2 (resp. w'^2) at the inflow together with the
  !> other (step1) changes only the downstream v'^2 (resp. w'^2) with unit slope, so the corrected
  !> control point should exactly cancel the step0 residual. Endpoints must stay pinned to target.
  Subroutine case_slopes_and_correction

    Real(Int64) :: u_t, v_t, w_t, u_b0, v_b0, w_b0

    If ( myid == 0 ) Write(*,'(A)') '=== scalar_correction_step matches a closed-form secant Newton step ==='

    n_bezier = 4
    Allocate( y_cp(n_bezier) )
    y_cp = (/ 0d0, 0.5d0, 1.5d0, 3d0 /)

    u_t = 0.6d0;  v_t = 0.5d0;  w_t = 0.3d0
    u_b0 = 0.55d0;  v_b0 = 0.4d0;  w_b0 = 0.35d0   ! step0 measurement: v'^2 low, w'^2 high vs. target

    Allocate( u_cp_target(n_bezier), x_cp_R22_target(n_bezier), x_cp_R33_target(n_bezier) )
    Allocate( x_cp_R22(n_bezier), x_cp_R33(n_bezier) )
    u_cp_target = u_t;  x_cp_R22_target = v_t;  x_cp_R33_target = w_t
    x_cp_R22 = v_t;  x_cp_R33 = w_t

    Allocate( stats_step0(3,n_bezier), stats_step1(3,n_bezier) )
    stats_step0(1,:) = u_b0;  stats_step0(2,:) = v_b0;  stats_step0(3,:) = w_b0
    ! step1: v'^2 and w'^2 doubled together; each downstream component responds only to its own
    ! input with unit slope (u'^2 stays at its step0 value -- it has no slope in this scheme)
    stats_step1 = stats_step0
    stats_step1(2,:) = stats_step0(2,:) + v_t
    stats_step1(3,:) = stats_step0(3,:) + w_t

    Allocate( slope_v(n_bezier), slope_w(n_bezier) )
    inflow_opt_iter = 1
    Call build_initial_slopes
    Call scalar_correction_step( stats_step0, 1d0 )

    Call expect( Abs(x_cp_R22(2) - (v_t - (v_b0-v_t))) < tol, &
         'interior control point 2: v''^2 correction cancels the step0 residual' )
    Call expect( Abs(x_cp_R33(2) - (w_t - (w_b0-w_t))) < tol, &
         'interior control point 2: w''^2 correction cancels the step0 residual' )
    Call expect( Abs(x_cp_R22(3) - (v_t - (v_b0-v_t))) < tol, &
         'interior control point 3: v''^2 correction cancels the step0 residual' )
    Call expect( Abs(x_cp_R33(3) - (w_t - (w_b0-w_t))) < tol, &
         'interior control point 3: w''^2 correction cancels the step0 residual' )

    Call expect( Abs(x_cp_R22(1)-v_t) < tol .And. Abs(x_cp_R33(1)-w_t) < tol, &
         'endpoint control point 1 stays pinned to the target' )
    Call expect( Abs(x_cp_R22(n_bezier)-v_t) < tol .And. Abs(x_cp_R33(n_bezier)-w_t) < tol, &
         'endpoint control point n_bezier stays pinned to the target' )

    Deallocate( y_cp, u_cp_target, x_cp_R22_target, x_cp_R33_target, x_cp_R22, x_cp_R33 )
    Deallocate( stats_step0, stats_step1, slope_v, slope_w )

  End Subroutine case_slopes_and_correction

  !> scalar_correction_step's scale argument must linearly scale the applied Newton step -- the
  !> mechanism the experimental iter>1 extension uses for its Robbins-Monro-style decaying step
  !> size (inflow_opt_relax/iter) instead of a fixed trust-region magnitude cap.
  Subroutine case_correction_scale

    Real(Int64) :: v_t, w_t, v_b0, w_b0, dx_full

    If ( myid == 0 ) Write(*,'(A)') '=== scalar_correction_step scale argument linearly scales the step ==='

    n_bezier = 4
    Allocate( y_cp(n_bezier) )
    y_cp = (/ 0d0, 0.5d0, 1.5d0, 3d0 /)

    v_t = 0.5d0;  w_t = 0.3d0
    v_b0 = 0.1d0;  w_b0 = 0.2d0

    Allocate( u_cp_target(n_bezier), x_cp_R22_target(n_bezier), x_cp_R33_target(n_bezier) )
    Allocate( x_cp_R22(n_bezier), x_cp_R33(n_bezier) )
    u_cp_target = 0.6d0;  x_cp_R22_target = v_t;  x_cp_R33_target = w_t
    x_cp_R22 = v_t;  x_cp_R33 = w_t

    Allocate( stats_step0(3,n_bezier), stats_step1(3,n_bezier) )
    stats_step0(1,:) = 0.6d0;  stats_step0(2,:) = v_b0;  stats_step0(3,:) = w_b0
    stats_step1 = stats_step0
    stats_step1(2,:) = stats_step0(2,:) + v_t
    stats_step1(3,:) = stats_step0(3,:) + w_t

    Allocate( slope_v(n_bezier), slope_w(n_bezier) )
    inflow_opt_iter = 1
    Call build_initial_slopes
    Call scalar_correction_step( stats_step0, 1d0 )
    dx_full = x_cp_R22(2) - v_t   ! full (unscaled) step just applied

    x_cp_R22 = v_t;  x_cp_R33 = w_t   ! reset before re-applying at half scale
    Call scalar_correction_step( stats_step0, 0.5d0 )

    Call expect( Abs( (x_cp_R22(2)-v_t) - 0.5d0*dx_full ) < tol, &
         'scale=0.5 applies exactly half the scale=1 step' )

    Deallocate( y_cp, u_cp_target, x_cp_R22_target, x_cp_R33_target, x_cp_R22, x_cp_R33 )
    Deallocate( stats_step0, stats_step1, slope_v, slope_w )

  End Subroutine case_correction_scale

  !> secant_update_slopes must satisfy the secant condition exactly at the control point whose
  !> step was updated: slope_new * Delta_x == Delta_y, independently for v'^2 and w'^2, for
  !> arbitrary (wrong) starting slopes and arbitrary Delta_x/Delta_y. Used only by the experimental
  !> iter>1 extension, but the update rule itself is unconditional.
  Subroutine case_secant_update_slopes

    Real(Int64) :: dxv, dy, stats_now(3,3)

    If ( myid == 0 ) Write(*,'(A)') '=== secant_update_slopes satisfies the secant condition ==='

    n_bezier = 3   ! single interior control point (m=2), minimal setup
    Allocate( x_cp_R22(n_bezier), x_cp_R33(n_bezier) )
    Allocate( x_prev_R22(n_bezier), x_prev_R33(n_bezier) )
    Allocate( slope_v(n_bezier), slope_w(n_bezier) )
    Allocate( stats_prev(3,n_bezier) )

    ! deliberately "wrong" starting slope guesses -- the update must correct them
    slope_v(2) = 0.2d0;  slope_w(2) = 0.5d0

    x_prev_R22(2) = 0.5d0;  x_cp_R22(2) = 0.7d0   ! Delta_x_v = 0.2
    x_prev_R33(2) = 0.6d0;  x_cp_R33(2) = 0.9d0   ! Delta_x_w = 0.3
    dxv = 0.2d0

    stats_prev(:,2) = (/ 1.0d0, 1.0d0, 1.0d0 /)
    dy = 0.05d0
    stats_now(:,2) = stats_prev(:,2)
    stats_now(2,2) = stats_prev(2,2) + dy    ! v'^2 response
    stats_now(3,2) = stats_prev(3,2) - 0.02d0  ! w'^2 response

    Call secant_update_slopes( stats_now )

    Call expect( Abs( slope_v(2)*dxv - dy ) < tol, &
         'updated v''^2 slope reproduces the observed response exactly along the step direction' )
    Call expect( Abs( slope_w(2)*0.3d0 - (-0.02d0) ) < tol, &
         'updated w''^2 slope reproduces the observed response exactly along the step direction' )

    Deallocate( x_cp_R22, x_cp_R33, x_prev_R22, x_prev_R33, slope_v, slope_w, stats_prev )

  End Subroutine case_secant_update_slopes

  !> inflow_opt_wall_excluded must flag a control point within inflow_opt_wall_exclude*wall_Ltaper
  !> of an active no-slip wall, and only when that wall is actually active (y_bc_type==1-style
  !> no-slip, not a free-slip/symmetry boundary) -- mirrors sem_fluctuation's own taper zone.
  Subroutine case_wall_excluded

    If ( myid == 0 ) Write(*,'(A)') '=== inflow_opt_wall_excluded matches the no-slip taper zone ==='

    n_bezier = 3
    Allocate( y_cp(n_bezier) )
    y_cp = (/ 0.01d0, 0.05d0, 0.20d0 /)

    inflow_opt_wall_exclude = 1d0
    wall_active_lo = 1;  wall_y_lo = 0d0;  wall_Ltaper_lo = 0.1d0   ! taper zone: y in [0, 0.1]
    wall_active_hi = 0

    Call expect( inflow_opt_wall_excluded(1), 'y_cp=0.01 (inside taper zone) is excluded' )
    Call expect( inflow_opt_wall_excluded(2), 'y_cp=0.05 (inside taper zone) is excluded' )
    Call expect( .Not. inflow_opt_wall_excluded(3), 'y_cp=0.20 (outside taper zone) is not excluded' )

    wall_active_lo = 0   ! inactive wall (e.g. free-slip/symmetry) excludes nothing regardless of height
    Call expect( .Not. inflow_opt_wall_excluded(1), 'no exclusion when the wall is not active (wall_active_lo=0)' )

    wall_active_lo = 0;  wall_active_hi = 0   ! reset module state for subsequent test cases
    Deallocate( y_cp )

  End Subroutine case_wall_excluded

  !> inflow_opt_max_rel_residual must return the largest |measured-target|/target over
  !> u'^2/v'^2/w'^2 and all interior control points, ignoring the endpoints.
  Subroutine case_max_rel_residual

    Real(Int64) :: stats_ref(3,4), r

    If ( myid == 0 ) Write(*,'(A)') '=== inflow_opt_max_rel_residual picks the worst-case relative error ==='

    n_bezier = 4
    Allocate( u_cp_target(n_bezier), x_cp_R22_target(n_bezier), x_cp_R33_target(n_bezier) )
    u_cp_target = 1.0d0;  x_cp_R22_target = 2.0d0;  x_cp_R33_target = 4.0d0

    ! endpoints (1,4) deliberately way off target -- must be ignored (loop only covers m=2,3)
    stats_ref(:,1) = (/ 100d0, 100d0, 100d0 /)
    stats_ref(:,4) = (/ 100d0, 100d0, 100d0 /)
    ! interior point 2: u off by 10%, v off by 25% (the worst), w exact
    stats_ref(:,2) = (/ 1.1d0, 2.5d0, 4.0d0 /)
    ! interior point 3: all within 5%
    stats_ref(:,3) = (/ 1.05d0, 2.05d0, 4.05d0 /)

    r = inflow_opt_max_rel_residual( stats_ref )

    Call expect( Abs(r - 0.25d0) < tol, 'worst-case relative residual is the 25% v''^2 error at control point 2' )

    Deallocate( u_cp_target, x_cp_R22_target, x_cp_R33_target )

  End Subroutine case_max_rel_residual

End Program verify_inflow_opt
