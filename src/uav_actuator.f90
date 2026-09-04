!> UAV actuator-disk body force: loading can be uniform or parabolic
!  tip-tapered, radially; the reaction force can be purely vertical (flat
!  disk) or a full 3-component vector when the disk auto-tilts to follow
!  its own prescribed path kinematics, with an optional in-plane swirl
!  (rotor-torque reaction) component. The disk centre can follow a
!  prescribed path, and the total thrust can follow its own prescribed
!  schedule instead of staying fixed.
!
!  The disk is represented by a ring of Lagrangian markers carrying a share
!  of the total kinematic thrust, spread onto the U/V/W momentum grids with
!  a normalized Gaussian kernel so that, for every marker, the sum of the
!  force density times the local cell volume over its kernel support
!  exactly reproduces that marker's thrust share -- conservative regardless
!  of grid non-uniformity or how the kernel support is split across MPI
!  ranks: each rank deposits force only into the cells it actually owns, but
!  normalizes against the marker's FULL kernel-box sum (computed arithmetically
!  from the uniform x/z spacing, not by reading neighbour-rank cells) so a
!  marker split across two ranks still contributes exactly its own thrust
!  share in total, not double.
!
!  Tilt (uav_tilt_active=1): the disk normal is derived each step from the
!  path's own kinematic acceleration (a differentially-flat point-mass
!  argument: the thrust vector must equal m*(a + g*yhat), so its direction
!  is n = normalize(ax, grav+ay, az), independent of mass). The Catmull-Rom
!  path interpolation is only C1 (velocity-continuous, not acceleration-
!  continuous), so the raw per-step acceleration has knot-to-knot jump
!  discontinuities; n is low-pass filtered (time constant uav_tilt_tau)
!  before use to keep the disk orientation (and hence the force it applies)
!  from jumping at path waypoints. With uav_tilt_active=0 (default) n is
!  fixed at (0,1,0) and every formula below reduces exactly to the original
!  flat, vertical-only-force model.
!
!  Known simplifications:
!   - swirl (uav_swirl_frac) direction is an arbitrary modelling choice, not
!     derived from any tracked rotor RPM or rotation sense
!   - kernel support is not wrapped across periodic x/z boundaries; keep the
!     disk at least uav_kernel_ncell*max(dx,dz) away from a periodic edge
!     for its whole path
!   - no OpenACC offload internal to this module; apply_uav_forcing_u/v/w run
!     on the host and the caller (equations.f90) syncs rhs_u/rhs_v/rhs_w with
!     the device around each call on GPU_POISSON builds
Module uav_actuator

  Use iso_fortran_env, Only : Int32, Int64
  Use global
  Use mpi, Only : myid

  Implicit None

  Integer(Int32) :: n_uav_markers = 0
  ! marker offsets relative to the (possibly moving/tilted) disk centre, in
  ! the disk-local frame (e1,e2 below); uav_mk_theta is each marker's
  ! azimuthal angle in that same frame (stored, not recomputed, since it is
  ! needed every force-application call for the swirl tangential direction).
  ! uav_mk_frac: each marker's FRACTIONAL share of the disk's total thrust
  ! (dimensionless, sums to 1 over all markers) -- the total is evaluated
  ! separately each call (fixed uav_hover_thrust, or uav_current_thrust(t)
  ! when uav_thrust_active=1) so a time-varying schedule doesn't require
  ! rebuilding the marker ring.
  Real   (Int64), Allocatable :: uav_mk_dx(:), uav_mk_dz(:), uav_mk_frac(:), uav_mk_theta(:)

  ! Phase 2 path (uav_path_active==1 only): waypoints "t x y z", sorted by
  ! increasing t; the disk centre is cubic-Hermite (Catmull-Rom tangent)
  ! interpolated between them and clamped outside [path_t(1),path_t(n_path)]
  Integer(Int32) :: n_path = 0
  Real   (Int64), Allocatable :: path_t(:), path_x(:), path_y(:), path_z(:)

  ! Thrust schedule (uav_thrust_active==1 only): waypoints "t T" [s, m^4/s^2
  ! kinematic thrust], same interpolation/clamping convention as the path
  Integer(Int32) :: n_thrust = 0
  Real   (Int64), Allocatable :: thrust_t(:), thrust_val(:)

  ! Persisted low-pass state for the auto-tilt disk normal (uav_tilt_active
  ! only); uav_tilt_last_t starts at a large-negative sentinel so the first
  ! call's effective filter gain saturates to 1 (i.e. the filter output is
  ! just the first raw sample, no special-case branch needed).
  Real   (Int64) :: uav_tilt_smooth(3) = (/ 0d0, 1d0, 0d0 /)
  Real   (Int64) :: uav_tilt_last_t    = -1d300

Contains

  !> Build the marker ring (uav_n_r radial bands x uav_n_theta azimuthal
  !  sectors), each carrying a thrust share proportional to its annulus-
  !  sector area (uniform disk loading) times a radial weight (uav_load_profile),
  !  and read the path/thrust files if active.
  Subroutine setup_uav

    Integer(Int32) :: ir, ith, n
    Real(Int64) :: dr, dtheta, r_mid, theta, area_disk, area_k, w_r, frac_sum

    n_uav_markers = uav_n_r * uav_n_theta
    Allocate( uav_mk_dx(n_uav_markers), uav_mk_dz(n_uav_markers), &
              uav_mk_frac(n_uav_markers), uav_mk_theta(n_uav_markers) )

    dr        = uav_disk_radius / Real(uav_n_r, Int64)
    dtheta    = 2d0*pi / Real(uav_n_theta, Int64)
    area_disk = pi * uav_disk_radius**2

    n = 0
    Do ir = 1, uav_n_r
       r_mid  = ( Real(ir,Int64) - 0.5d0 ) * dr
       ! annulus-sector area for this radial band, split evenly over its uav_n_theta sectors
       area_k = pi * ( (Real(ir,Int64)*dr)**2 - (Real(ir-1,Int64)*dr)**2 ) / Real(uav_n_theta, Int64)
       If ( uav_load_profile == 1 ) Then
          w_r = Max( 0d0, 1d0 - (r_mid/uav_disk_radius)**2 )   ! parabolic tip taper, zero at the tip
       Else
          w_r = 1d0                                             ! uniform loading (default)
       End If
       Do ith = 1, uav_n_theta
          theta = ( Real(ith,Int64) - 0.5d0 ) * dtheta
          n = n + 1
          uav_mk_dx(n)    = r_mid*Cos(theta)
          uav_mk_dz(n)    = r_mid*Sin(theta)
          uav_mk_theta(n) = theta
          uav_mk_frac(n)  = area_k*w_r/area_disk   ! renormalized below (exact only for w_r==1)
       End Do
    End Do

    frac_sum = Sum(uav_mk_frac)
    uav_mk_frac = uav_mk_frac / frac_sum   ! exact no-op when uniform; restores sum=1 under tip-tapering

    If ( uav_path_active >= 1 ) Then
       If ( Len_Trim(uav_path_file) == 0 ) Stop 'ERROR: uav_path_active=1 but uav_path_file is empty'
       Call read_uav_path
    End If

    If ( uav_thrust_active >= 1 ) Then
       If ( Len_Trim(uav_thrust_file) == 0 ) Stop 'ERROR: uav_thrust_active=1 but uav_thrust_file is empty'
       Call read_uav_thrust
    End If

    If ( myid==0 ) Then
       Write(*,*) 'UAV: actuator disk active -- markers:', n_uav_markers, ' radius:', uav_disk_radius
       If ( uav_path_active >= 1 ) Then
          Write(*,*) 'UAV: following path from ', Trim(uav_path_file), ' (', n_path, ' waypoints)'
       Else
          Write(*,*) 'UAV: static centre:', uav_xc, uav_yc, uav_zc
       End If
       If ( uav_thrust_active >= 1 ) Then
          Write(*,*) 'UAV: following thrust schedule from ', Trim(uav_thrust_file), ' (', n_thrust, ' points)'
       Else
          Write(*,*) 'UAV: fixed kinematic thrust:', uav_hover_thrust
       End If
       If ( uav_load_profile == 1 ) Write(*,*) 'UAV: parabolic tip-tapered radial loading'
       If ( uav_tilt_active >= 1 ) Write(*,*) 'UAV: auto-tilt active, filter tau =', uav_tilt_tau
       If ( uav_swirl_frac /= 0d0 ) Write(*,*) 'UAV: swirl fraction =', uav_swirl_frac
    End If

  End Subroutine setup_uav

  !> Read the waypoint file: free-form text, blank/'#' lines skipped, each
  !  data row "t x y z", sorted by strictly increasing t. Every rank reads
  !  the (small, shared) file independently, mirroring sem.f90's
  !  read_mean_profile convention for inflow_profile_file.
  Subroutine read_uav_path

    Integer(Int32) :: unit_in, ios, n
    Real   (Int64) :: col(4)
    Character(300) :: line

    If ( myid==0 ) Write(*,'(A,A)') ' Reading UAV path from ', Trim(uav_path_file)

    Open(newunit=unit_in, file=Trim(uav_path_file), status='old', action='read', iostat=ios)
    If ( ios /= 0 ) Stop 'ERROR: cannot open uav_path_file'

    n_path = 0
    Do
       Read(unit_in,'(A)',iostat=ios) line
       If ( ios /= 0 ) Exit
       line = Adjustl(line)
       If ( Len_Trim(line) == 0 .Or. line(1:1) == '#' ) Cycle
       n_path = n_path + 1
    End Do
    If ( n_path < 2 ) Stop 'ERROR: uav_path_file: fewer than 2 waypoints found'

    Allocate( path_t(n_path), path_x(n_path), path_y(n_path), path_z(n_path) )

    Rewind(unit_in)
    n = 0
    Do
       Read(unit_in,'(A)',iostat=ios) line
       If ( ios /= 0 ) Exit
       line = Adjustl(line)
       If ( Len_Trim(line) == 0 .Or. line(1:1) == '#' ) Cycle
       n = n + 1
       Read(line,*,iostat=ios) col(1:4)
       If ( ios /= 0 ) Stop 'ERROR: uav_path_file: failed to parse a data row (need 4 columns: t x y z)'
       path_t(n) = col(1);  path_x(n) = col(2);  path_y(n) = col(3);  path_z(n) = col(4)
       If ( n > 1 ) Then
          If ( path_t(n) <= path_t(n-1) ) Stop 'ERROR: uav_path_file: t column must be strictly increasing'
       End If
    End Do
    Close(unit_in)

  End Subroutine read_uav_path

  !> Read the thrust schedule: free-form text, blank/'#' lines skipped,
  !  each data row "t T" (T = kinematic thrust [m^4/s^2]), sorted by
  !  strictly increasing t. Same per-rank-independent-read convention as
  !  read_uav_path.
  Subroutine read_uav_thrust

    Integer(Int32) :: unit_in, ios, n
    Real   (Int64) :: col(2)
    Character(300) :: line

    If ( myid==0 ) Write(*,'(A,A)') ' Reading UAV thrust schedule from ', Trim(uav_thrust_file)

    Open(newunit=unit_in, file=Trim(uav_thrust_file), status='old', action='read', iostat=ios)
    If ( ios /= 0 ) Stop 'ERROR: cannot open uav_thrust_file'

    n_thrust = 0
    Do
       Read(unit_in,'(A)',iostat=ios) line
       If ( ios /= 0 ) Exit
       line = Adjustl(line)
       If ( Len_Trim(line) == 0 .Or. line(1:1) == '#' ) Cycle
       n_thrust = n_thrust + 1
    End Do
    If ( n_thrust < 2 ) Stop 'ERROR: uav_thrust_file: fewer than 2 points found'

    Allocate( thrust_t(n_thrust), thrust_val(n_thrust) )

    Rewind(unit_in)
    n = 0
    Do
       Read(unit_in,'(A)',iostat=ios) line
       If ( ios /= 0 ) Exit
       line = Adjustl(line)
       If ( Len_Trim(line) == 0 .Or. line(1:1) == '#' ) Cycle
       n = n + 1
       Read(line,*,iostat=ios) col(1:2)
       If ( ios /= 0 ) Stop 'ERROR: uav_thrust_file: failed to parse a data row (need 2 columns: t T)'
       thrust_t(n) = col(1);  thrust_val(n) = col(2)
       If ( n > 1 ) Then
          If ( thrust_t(n) <= thrust_t(n-1) ) Stop 'ERROR: uav_thrust_file: t column must be strictly increasing'
       End If
    End Do
    Close(unit_in)

  End Subroutine read_uav_thrust

  !> Cubic Hermite (Catmull-Rom tangent) interpolation of one path
  !  component at time tt, clamped outside [vals_t(1),vals_t(n)].
  !  Endpoint tangents use the one-sided (forward/backward) difference in
  !  place of the missing neighbour.
  Pure Function hermite_eval(tt, vals_t, vals_p, n) Result(val)

    Integer(Int32), Intent(In) :: n
    Real   (Int64), Intent(In) :: tt, vals_t(n), vals_p(n)
    Real   (Int64) :: val

    Integer(Int32) :: i
    Real(Int64) :: t0,t1,t2,t3, p0,p1,p2,p3, h, s, m1, m2
    Real(Int64) :: h00, h10, h01, h11

    If ( tt <= vals_t(1) ) Then
       val = vals_p(1);  Return
    Else If ( tt >= vals_t(n) ) Then
       val = vals_p(n);  Return
    End If

    ! locate bracket i such that vals_t(i) <= tt <= vals_t(i+1) (small n_path expected: linear scan)
    i = 1
    Do While ( i < n-1 .And. vals_t(i+1) < tt )
       i = i + 1
    End Do

    t1 = vals_t(i);    t2 = vals_t(i+1)
    p1 = vals_p(i);    p2 = vals_p(i+1)
    h  = t2 - t1
    s  = (tt - t1) / h

    If ( i > 1 ) Then
       t0 = vals_t(i-1);  p0 = vals_p(i-1)
       m1 = (p2 - p0) / (t2 - t0)
    Else
       m1 = (p2 - p1) / h   ! forward difference at the first segment
    End If

    If ( i+2 <= n ) Then
       t3 = vals_t(i+2);  p3 = vals_p(i+2)
       m2 = (p3 - p1) / (t3 - t1)
    Else
       m2 = (p2 - p1) / h   ! backward difference at the last segment
    End If

    h00 =  2d0*s**3 - 3d0*s**2 + 1d0
    h10 =       s**3 - 2d0*s**2 + s
    h01 = -2d0*s**3 + 3d0*s**2
    h11 =       s**3 -      s**2

    val = h00*p1 + h10*h*m1 + h01*p2 + h11*h*m2

  End Function hermite_eval

  !> Second time-derivative of the same cubic-Hermite segment used by
  !  hermite_eval, for the auto-tilt acceleration term. A cubic's second
  !  derivative is linear in s, so this is exact on each segment, but (like
  !  the underlying Catmull-Rom tangent itself) is only C0 continuous across
  !  segment knots -- the low-pass filter in uav_disk_state exists precisely
  !  to tame that. Zero outside [vals_t(1),vals_t(n)] (position is clamped
  !  constant there, so acceleration is identically zero).
  Pure Function hermite_accel(tt, vals_t, vals_p, n) Result(acc)

    Integer(Int32), Intent(In) :: n
    Real   (Int64), Intent(In) :: tt, vals_t(n), vals_p(n)
    Real   (Int64) :: acc

    Integer(Int32) :: i
    Real(Int64) :: t0,t1,t2,t3, p0,p1,p2,p3, h, s, m1, m2

    If ( tt <= vals_t(1) .Or. tt >= vals_t(n) ) Then
       acc = 0d0;  Return
    End If

    i = 1
    Do While ( i < n-1 .And. vals_t(i+1) < tt )
       i = i + 1
    End Do

    t1 = vals_t(i);    t2 = vals_t(i+1)
    p1 = vals_p(i);    p2 = vals_p(i+1)
    h  = t2 - t1
    s  = (tt - t1) / h

    If ( i > 1 ) Then
       t0 = vals_t(i-1);  p0 = vals_p(i-1)
       m1 = (p2 - p0) / (t2 - t0)
    Else
       m1 = (p2 - p1) / h
    End If

    If ( i+2 <= n ) Then
       t3 = vals_t(i+2);  p3 = vals_p(i+2)
       m2 = (p3 - p1) / (t3 - t1)
    Else
       m2 = (p2 - p1) / h
    End If

    ! d2/dt2 = (1/h^2) * [ h00''(s)*p1 + h10''(s)*h*m1 + h01''(s)*p2 + h11''(s)*h*m2 ]
    acc = ( (12d0*s-6d0)*p1 + (6d0*s-4d0)*h*m1 + (-12d0*s+6d0)*p2 + (6d0*s-2d0)*h*m2 ) / (h*h)

  End Function hermite_accel

  !> Current disk centre: the path-interpolated position (uav_path_active=1)
  !  or the fixed (uav_xc,uav_yc,uav_zc) otherwise.
  Subroutine uav_current_center(tt, xc, yc, zc)

    Real(Int64), Intent(In)  :: tt
    Real(Int64), Intent(Out) :: xc, yc, zc

    If ( uav_path_active >= 1 ) Then
       xc = hermite_eval(tt, path_t, path_x, n_path)
       yc = hermite_eval(tt, path_t, path_y, n_path)
       zc = hermite_eval(tt, path_t, path_z, n_path)
    Else
       xc = uav_xc;  yc = uav_yc;  zc = uav_zc
    End If

  End Subroutine uav_current_center

  !> Current total kinematic thrust: the schedule interpolated from
  !  uav_thrust_file (uav_thrust_active=1) or the fixed uav_hover_thrust
  !  otherwise.
  Pure Function uav_current_thrust(tt) Result(Qtot)

    Real(Int64), Intent(In) :: tt
    Real(Int64) :: Qtot

    If ( uav_thrust_active >= 1 ) Then
       Qtot = hermite_eval(tt, thrust_t, thrust_val, n_thrust)
    Else
       Qtot = uav_hover_thrust
    End If

  End Function uav_current_thrust

  !> Orthonormal in-disk-plane frame (e1,e2) for a given unit disk normal n,
  !  chosen so that at n=(0,1,0) (no tilt) e1=(1,0,0), e2=(0,0,1) exactly --
  !  i.e. identical to the original flat-disk marker layout (dx along x,
  !  dz along z). Degenerates only if n is within a few degrees of +-xhat,
  !  which a physically modest rotor tilt never approaches.
  Pure Subroutine uav_disk_frame(nvec, e1, e2)

    Real(Int64), Intent(In)  :: nvec(3)
    Real(Int64), Intent(Out) :: e1(3), e2(3)

    Real(Int64) :: norm2

    norm2 = Sqrt( nvec(2)**2 + nvec(3)**2 )

    e2(1) = 0d0
    e2(2) = -nvec(3)/norm2
    e2(3) =  nvec(2)/norm2

    ! e1 = n x e2 (completes a right-handed frame; reduces to (1,0,0) at n=(0,1,0))
    e1(1) = nvec(2)*e2(3) - nvec(3)*e2(2)
    e1(2) = nvec(3)*e2(1) - nvec(1)*e2(3)
    e1(3) = nvec(1)*e2(2) - nvec(2)*e2(1)

  End Subroutine uav_disk_frame

  !> Evaluate the disk's current centre, total thrust, normal and in-plane
  !  frame at time tt. With uav_tilt_active=0 the normal is fixed at
  !  (0,1,0) (flat disk, exactly the original model). With uav_tilt_active=1
  !  the normal is derived from the path's own kinematic acceleration
  !  (differentially-flat point-mass tilt) and low-pass filtered in time
  !  (persisted module state uav_tilt_smooth/uav_tilt_last_t) so repeated
  !  calls at the same tt (once per rhs_u/rhs_v/rhs_w per RK sub-stage)
  !  reuse the same filtered value instead of re-blending it three times.
  Subroutine uav_disk_state(tt, xc, yc, zc, Qtot, nvec, e1, e2)

    Real(Int64), Intent(In)  :: tt
    Real(Int64), Intent(Out) :: xc, yc, zc, Qtot, nvec(3), e1(3), e2(3)

    Real(Int64) :: ax, ay, az, n_raw(3), norm_raw, dt_call, alpha

    Call uav_current_center(tt, xc, yc, zc)
    Qtot = uav_current_thrust(tt)

    If ( uav_tilt_active >= 1 .And. uav_path_active >= 1 ) Then

       ax = hermite_accel(tt, path_t, path_x, n_path)
       ay = hermite_accel(tt, path_t, path_y, n_path)
       az = hermite_accel(tt, path_t, path_z, n_path)

       n_raw(1) = ax;  n_raw(2) = grav + ay;  n_raw(3) = az
       norm_raw = Sqrt( n_raw(1)**2 + n_raw(2)**2 + n_raw(3)**2 )
       If ( norm_raw > 1d-30 ) n_raw = n_raw / norm_raw

       If ( tt /= uav_tilt_last_t ) Then
          dt_call = tt - uav_tilt_last_t
          alpha   = dt_call / (uav_tilt_tau + dt_call)   ! in (0,1]; ~1 on the first (huge dt_call) call
          uav_tilt_smooth = uav_tilt_smooth + alpha*(n_raw - uav_tilt_smooth)
          uav_tilt_smooth = uav_tilt_smooth / Sqrt( Sum(uav_tilt_smooth**2) )
          uav_tilt_last_t = tt
       End If

       nvec = uav_tilt_smooth
       Call uav_disk_frame(nvec, e1, e2)

    Else

       nvec = (/ 0d0, 1d0, 0d0 /)
       e1   = (/ 1d0, 0d0, 0d0 /)
       e2   = (/ 0d0, 0d0, 1d0 /)

    End If

  End Subroutine uav_disk_state

  !> Add the disk's reaction force into rhs_u (x-component). Called once
  !  per RK3 sub-stage from compute_rhs_u, guarded by uav_active. See
  !  apply_uav_forcing_v for the detailed normalization/deposit convention;
  !  this mirrors it exactly but on the x-face/yg-centre/zg-centre grid u
  !  lives on (only the x index is a face; y,z are cell centres, same
  !  convention as rhs_v's x,z indices).
  Subroutine apply_uav_forcing_u(rhs_u)

    Real(Int64), Dimension(2:nx-1,2:nyg-1,2:nzg-1), Intent(InOut) :: rhs_u

    Integer(Int32) :: n, i, j, k, i0, j0, k0
    Integer(Int32) :: ilo, ihi, jlo, jhi, klo, khi, ifull_lo, ifull_hi, kfull_lo, kfull_hi
    Real(Int64) :: xc, yc, zc, Qtot, nvec(3), e1(3), e2(3)
    Real(Int64) :: xp, yp, zp, Qn, e_theta(3), Fx
    Real(Int64) :: sigx, sigz, sigy, dxp, dyp, dzp, wgt, norm, dVj, best

    Call uav_disk_state(t, xc, yc, zc, Qtot, nvec, e1, e2)

    Do n = 1, n_uav_markers

       Qn = Qtot * uav_mk_frac(n)
       xp = xc + uav_mk_dx(n)*e1(1) + uav_mk_dz(n)*e2(1)
       yp = yc + uav_mk_dx(n)*e1(2) + uav_mk_dz(n)*e2(2)
       zp = zc + uav_mk_dx(n)*e1(3) + uav_mk_dz(n)*e2(3)

       e_theta = -Sin(uav_mk_theta(n))*e1 + Cos(uav_mk_theta(n))*e2
       Fx = -Qn*nvec(1) + uav_swirl_frac*Qn*e_theta(1)

       ! nearest local x-face index (uniform spacing, face offset by half a cell from the xg centre grid)
       i0 = Nint( (xp - xg(2))/dx + 1.5d0 )
       k0 = Nint( (zp - zg(2))/dz )        + 2

       ! nearest yg cell-centre index (non-uniform grid: linear search, yg not MPI-decomposed)
       j0 = 2
       best = Abs( yg(2) - yp )
       Do j = 3, nyg-1
          If ( Abs(yg(j)-yp) < best ) Then
             best = Abs(yg(j)-yp)
             j0   = j
          End If
       End Do

       ilo = Max(2,   i0-uav_kernel_ncell);  ihi = Min(nx-1,  i0+uav_kernel_ncell)
       jlo = Max(2,   j0-uav_kernel_ncell);  jhi = Min(nyg-1, j0+uav_kernel_ncell)
       klo = Max(2,   k0-uav_kernel_ncell);  khi = Min(nzg-1, k0+uav_kernel_ncell)

       If ( ilo > ihi .Or. jlo > jhi .Or. klo > khi ) Cycle

       ifull_lo = i0 - uav_kernel_ncell;  ifull_hi = i0 + uav_kernel_ncell
       kfull_lo = k0 - uav_kernel_ncell;  kfull_hi = k0 + uav_kernel_ncell

       sigx = Max(uav_kernel_ncell,1) * dx                * 0.5d0
       sigz = Max(uav_kernel_ncell,1) * dz                * 0.5d0
       sigy = Max(uav_kernel_ncell,1) * (y(j0)-y(j0-1))   * 0.5d0

       norm = 0d0
       Do k = kfull_lo, kfull_hi
          dzp = zg(2) + Real(k-2,Int64)*dz - zp
          Do j = jlo, jhi
             dVj = dx * dz * ( y(j) - y(j-1) )
             dyp = yg(j) - yp
             Do i = ifull_lo, ifull_hi
                dxp = xg(2) + ( Real(i,Int64) - 1.5d0 )*dx - xp
                wgt = Exp( -0.5d0*( (dxp/sigx)**2 + (dyp/sigy)**2 + (dzp/sigz)**2 ) )
                norm = norm + wgt*dVj
             End Do
          End Do
       End Do

       If ( norm < 1d-30 ) Cycle

       Do k = klo, khi
          dzp = zg(k) - zp
          Do j = jlo, jhi
             dyp = yg(j) - yp
             Do i = ilo, ihi
                dxp = x(i) - xp
                wgt = Exp( -0.5d0*( (dxp/sigx)**2 + (dyp/sigy)**2 + (dzp/sigz)**2 ) )
                rhs_u(i,j,k) = rhs_u(i,j,k) + Fx * wgt / norm
             End Do
          End Do
       End Do

    End Do

  End Subroutine apply_uav_forcing_u

  !> Add the disk's reaction force into rhs_v (y-component, the only
  !  component the original flat/untilted model had). Called once per RK3
  !  sub-stage from compute_rhs_v, guarded by uav_active.
  Subroutine apply_uav_forcing_v(rhs_v)

    Real(Int64), Dimension(2:nxg-1,2:ny-1,2:nzg-1), Intent(InOut) :: rhs_v

    Integer(Int32) :: n, i, j, k, i0, j0, k0
    Integer(Int32) :: ilo, ihi, jlo, jhi, klo, khi
    Integer(Int32) :: ifull_lo, ifull_hi, kfull_lo, kfull_hi
    Real(Int64) :: xc, yc, zc, Qtot, nvec(3), e1(3), e2(3)
    Real(Int64) :: xp, yp, zp, Qn, e_theta(3), Fy
    Real(Int64) :: sigx, sigz, sigy, dxp, dyp, dzp, wgt, norm, dVj
    Real(Int64) :: best

    Call uav_disk_state(t, xc, yc, zc, Qtot, nvec, e1, e2)

    Do n = 1, n_uav_markers

       Qn = Qtot * uav_mk_frac(n)
       xp = xc + uav_mk_dx(n)*e1(1) + uav_mk_dz(n)*e2(1)
       yp = yc + uav_mk_dx(n)*e1(2) + uav_mk_dz(n)*e2(2)
       zp = zc + uav_mk_dx(n)*e1(3) + uav_mk_dz(n)*e2(3)

       e_theta = -Sin(uav_mk_theta(n))*e1 + Cos(uav_mk_theta(n))*e2
       Fy = -Qn*nvec(2) + uav_swirl_frac*Qn*e_theta(2)

       ! locate the nearest local cell-centre index in x/z (uniform spacing)
       i0 = Nint( (xp - xg(2))/dx ) + 2
       k0 = Nint( (zp - zg(2))/dz ) + 2

       ! locate the nearest local y-face index (non-uniform grid: linear search over this rank's small local y array)
       j0 = 2
       best = Abs( y(2) - yp )
       Do j = 3, ny-1
          If ( Abs(y(j)-yp) < best ) Then
             best = Abs(y(j)-yp)
             j0   = j
          End If
       End Do

       ilo = Max(2,    i0-uav_kernel_ncell);  ihi = Min(nxg-1, i0+uav_kernel_ncell)
       jlo = Max(2,    j0-uav_kernel_ncell);  jhi = Min(ny-1,  j0+uav_kernel_ncell)
       klo = Max(2,    k0-uav_kernel_ncell);  khi = Min(nzg-1, k0+uav_kernel_ncell)

       If ( ilo > ihi .Or. jlo > jhi .Or. klo > khi ) Cycle   ! marker's support does not overlap this rank

       ! full (un-clipped by this rank's x/z sub-domain) kernel box, for the normalization below
       ifull_lo = i0 - uav_kernel_ncell;  ifull_hi = i0 + uav_kernel_ncell
       kfull_lo = k0 - uav_kernel_ncell;  kfull_hi = k0 + uav_kernel_ncell

       sigx = Max(uav_kernel_ncell,1) * dx                    * 0.5d0
       sigz = Max(uav_kernel_ncell,1) * dz                    * 0.5d0
       sigy = Max(uav_kernel_ncell,1) * (yg(j0+1)-yg(j0))     * 0.5d0

       ! Pass 1: normalization so that sum(wgt*dV) over the support equals the marker's thrust share exactly.
       ! Summed over the marker's FULL kernel box (ifull_lo:ifull_hi, kfull_lo:kfull_hi), not just the
       ! ilo:ihi/klo:khi cells this rank owns: x/z are uniform and dx/dz are the same on every rank, so
       ! xg(i)-xp for any integer i is exactly xg(2)+(i-2)*dx-xp regardless of whether i is a valid local
       ! index here -- computed arithmetically below without touching xg/zg out of bounds. Using only the
       ! locally-owned range here would make each rank normalize its own partial sum back up to the full
       ! thrust share, double-counting (or worse) whenever a marker's support straddles two ranks in x/z.
       ! y is not MPI-decomposed, so jlo:jhi (already clipped at the physical wall, not a rank seam) is fine.
       norm = 0d0
       Do k = kfull_lo, kfull_hi
          dzp = zg(2) + Real(k-2,Int64)*dz - zp
          Do j = jlo, jhi
             dVj = dx * dz * ( yg(j+1) - yg(j) )
             dyp = y(j) - yp
             Do i = ifull_lo, ifull_hi
                dxp = xg(2) + Real(i-2,Int64)*dx - xp
                wgt = Exp( -0.5d0*( (dxp/sigx)**2 + (dyp/sigy)**2 + (dzp/sigz)**2 ) )
                norm = norm + wgt*dVj
             End Do
          End Do
       End Do

       If ( norm < 1d-30 ) Cycle   ! kernel support present but entirely between/outside sampled points (degenerate); nothing to add

       ! Pass 2: spread the marker's y-force share as an acceleration density
       Do k = klo, khi
          dzp = zg(k) - zp
          Do j = jlo, jhi
             dyp = y(j) - yp
             Do i = ilo, ihi
                dxp = xg(i) - xp
                wgt = Exp( -0.5d0*( (dxp/sigx)**2 + (dyp/sigy)**2 + (dzp/sigz)**2 ) )
                rhs_v(i,j,k) = rhs_v(i,j,k) + Fy * wgt / norm
             End Do
          End Do
       End Do

    End Do

  End Subroutine apply_uav_forcing_v

  !> Add the disk's reaction force into rhs_w (z-component). Called once
  !  per RK3 sub-stage from compute_rhs_w, guarded by uav_active. Mirrors
  !  apply_uav_forcing_v on the xg-centre/yg-centre/z-face grid w lives on
  !  (only the z index is a face; x,y are cell centres).
  Subroutine apply_uav_forcing_w(rhs_w)

    Real(Int64), Dimension(2:nxg-1,2:nyg-1,2:nz-1), Intent(InOut) :: rhs_w

    Integer(Int32) :: n, i, j, k, i0, j0, k0
    Integer(Int32) :: ilo, ihi, jlo, jhi, klo, khi, ifull_lo, ifull_hi, kfull_lo, kfull_hi
    Real(Int64) :: xc, yc, zc, Qtot, nvec(3), e1(3), e2(3)
    Real(Int64) :: xp, yp, zp, Qn, e_theta(3), Fz
    Real(Int64) :: sigx, sigz, sigy, dxp, dyp, dzp, wgt, norm, dVj, best

    Call uav_disk_state(t, xc, yc, zc, Qtot, nvec, e1, e2)

    Do n = 1, n_uav_markers

       Qn = Qtot * uav_mk_frac(n)
       xp = xc + uav_mk_dx(n)*e1(1) + uav_mk_dz(n)*e2(1)
       yp = yc + uav_mk_dx(n)*e1(2) + uav_mk_dz(n)*e2(2)
       zp = zc + uav_mk_dx(n)*e1(3) + uav_mk_dz(n)*e2(3)

       e_theta = -Sin(uav_mk_theta(n))*e1 + Cos(uav_mk_theta(n))*e2
       Fz = -Qn*nvec(3) + uav_swirl_frac*Qn*e_theta(3)

       i0 = Nint( (xp - xg(2))/dx )        + 2
       k0 = Nint( (zp - zg(2))/dz + 1.5d0 )

       j0 = 2
       best = Abs( yg(2) - yp )
       Do j = 3, nyg-1
          If ( Abs(yg(j)-yp) < best ) Then
             best = Abs(yg(j)-yp)
             j0   = j
          End If
       End Do

       ilo = Max(2,   i0-uav_kernel_ncell);  ihi = Min(nxg-1, i0+uav_kernel_ncell)
       jlo = Max(2,   j0-uav_kernel_ncell);  jhi = Min(nyg-1, j0+uav_kernel_ncell)
       klo = Max(2,   k0-uav_kernel_ncell);  khi = Min(nz-1,  k0+uav_kernel_ncell)

       If ( ilo > ihi .Or. jlo > jhi .Or. klo > khi ) Cycle

       ifull_lo = i0 - uav_kernel_ncell;  ifull_hi = i0 + uav_kernel_ncell
       kfull_lo = k0 - uav_kernel_ncell;  kfull_hi = k0 + uav_kernel_ncell

       sigx = Max(uav_kernel_ncell,1) * dx                * 0.5d0
       sigz = Max(uav_kernel_ncell,1) * dz                * 0.5d0
       sigy = Max(uav_kernel_ncell,1) * (y(j0)-y(j0-1))   * 0.5d0

       norm = 0d0
       Do k = kfull_lo, kfull_hi
          dzp = zg(2) + ( Real(k,Int64) - 1.5d0 )*dz - zp
          Do j = jlo, jhi
             dVj = dx * dz * ( y(j) - y(j-1) )
             dyp = yg(j) - yp
             Do i = ifull_lo, ifull_hi
                dxp = xg(2) + Real(i-2,Int64)*dx - xp
                wgt = Exp( -0.5d0*( (dxp/sigx)**2 + (dyp/sigy)**2 + (dzp/sigz)**2 ) )
                norm = norm + wgt*dVj
             End Do
          End Do
       End Do

       If ( norm < 1d-30 ) Cycle

       Do k = klo, khi
          dzp = z(k) - zp
          Do j = jlo, jhi
             dyp = yg(j) - yp
             Do i = ilo, ihi
                dxp = xg(i) - xp
                wgt = Exp( -0.5d0*( (dxp/sigx)**2 + (dyp/sigy)**2 + (dzp/sigz)**2 ) )
                rhs_w(i,j,k) = rhs_w(i,j,k) + Fz * wgt / norm
             End Do
          End Do
       End Do

    End Do

  End Subroutine apply_uav_forcing_w

End Module uav_actuator
