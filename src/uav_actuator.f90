!> UAV actuator-disk body force (Phase 1: uniform loading, vertical-only
!  reaction force; Phase 2: the disk centre can follow a prescribed path,
!  and the total thrust can follow its own prescribed schedule instead of
!  staying fixed -- see docs/UAV_ActuatorDisk_Design.md).
!
!  The disk is represented by a ring of Lagrangian markers carrying a share
!  of the total kinematic thrust, spread onto the V (wall-normal momentum)
!  grid with a normalized Gaussian kernel so that, for every marker, the sum
!  of the force density times the local cell volume over its kernel support
!  exactly reproduces that marker's thrust share -- conservative regardless
!  of grid non-uniformity or how the kernel support is split across MPI
!  ranks (a marker whose support straddles two ranks contributes only the
!  cells each rank actually owns; no explicit rank-ownership bookkeeping is
!  needed because xg/zg/y are already each rank's own local slice of the
!  global grid).
!
!  Known simplifications (see design doc for later phases):
!   - the disk stays horizontal (normal = +y) even while translating along
!     a path -- correct for vertical takeoff/climb/hover/descent, not yet
!     for a tilted cruise segment
!   - only the vertical (rhs_v) reaction force is applied
!   - kernel support is not wrapped across periodic x/z boundaries; keep the
!     disk at least uav_kernel_ncell*max(dx,dz) away from a periodic edge
!     for its whole path
!   - no OpenACC offload yet (CPU builds only; GPU_POISSON builds still run
!     this module on the host)
Module uav_actuator

  Use iso_fortran_env, Only : Int32, Int64
  Use global
  Use mpi, Only : myid

  Implicit None

  Integer(Int32) :: n_uav_markers = 0
  ! marker offsets relative to the (possibly moving) disk centre, and each
  ! marker's kinematic thrust share [m^4/s^2]
  ! uav_mk_frac: each marker's FRACTIONAL share of the disk's total thrust
  ! (dimensionless, sums to 1 over all markers) -- the total is evaluated
  ! separately each call (fixed uav_hover_thrust, or uav_current_thrust(t)
  ! when uav_thrust_active=1) so a time-varying schedule doesn't require
  ! rebuilding the marker ring.
  Real   (Int64), Allocatable :: uav_mk_dx(:), uav_mk_dz(:), uav_mk_frac(:)

  ! Phase 2 path (uav_path_active==1 only): waypoints "t x y z", sorted by
  ! increasing t; the disk centre is cubic-Hermite (Catmull-Rom tangent)
  ! interpolated between them and clamped outside [path_t(1),path_t(n_path)]
  Integer(Int32) :: n_path = 0
  Real   (Int64), Allocatable :: path_t(:), path_x(:), path_y(:), path_z(:)

  ! Thrust schedule (uav_thrust_active==1 only): waypoints "t T" [s, m^4/s^2
  ! kinematic thrust], same interpolation/clamping convention as the path
  Integer(Int32) :: n_thrust = 0
  Real   (Int64), Allocatable :: thrust_t(:), thrust_val(:)

Contains

  !> Build the marker ring (uav_n_r radial bands x uav_n_theta azimuthal
  !  sectors), each carrying a thrust share proportional to its annulus-
  !  sector area (uniform disk loading), and read the path file if active.
  Subroutine setup_uav

    Integer(Int32) :: ir, ith, n
    Real(Int64) :: dr, dtheta, r_mid, theta, area_disk, area_k

    n_uav_markers = uav_n_r * uav_n_theta
    Allocate( uav_mk_dx(n_uav_markers), uav_mk_dz(n_uav_markers), uav_mk_frac(n_uav_markers) )

    dr        = uav_disk_radius / Real(uav_n_r, Int64)
    dtheta    = 2d0*pi / Real(uav_n_theta, Int64)
    area_disk = pi * uav_disk_radius**2

    n = 0
    Do ir = 1, uav_n_r
       r_mid  = ( Real(ir,Int64) - 0.5d0 ) * dr
       ! annulus-sector area for this radial band, split evenly over its uav_n_theta sectors
       area_k = pi * ( (Real(ir,Int64)*dr)**2 - (Real(ir-1,Int64)*dr)**2 ) / Real(uav_n_theta, Int64)
       Do ith = 1, uav_n_theta
          theta = ( Real(ith,Int64) - 0.5d0 ) * dtheta
          n = n + 1
          uav_mk_dx(n)   = r_mid*Cos(theta)
          uav_mk_dz(n)   = r_mid*Sin(theta)
          uav_mk_frac(n) = area_k/area_disk
       End Do
    End Do

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

  !> Add the disk's vertical reaction force (downward on the fluid, opposing
  !  the upward thrust that supports the vehicle) into rhs_v, at the disk's
  !  current centre and total thrust (global `t`, from equations.f90's use
  !  of global). Called once per RK3 sub-stage from compute_rhs_v, guarded
  !  by uav_active.
  Subroutine apply_uav_forcing(rhs_v)

    Real(Int64), Dimension(2:nxg-1,2:ny-1,2:nzg-1), Intent(InOut) :: rhs_v

    Integer(Int32) :: n, i, j, k, i0, j0, k0
    Integer(Int32) :: ilo, ihi, jlo, jhi, klo, khi
    Real(Int64) :: xc, yc, zc, xp, yp, zp, Qtot, Qn
    Real(Int64) :: sigx, sigz, sigy, dxp, dyp, dzp, wgt, norm, dVj
    Real(Int64) :: best

    Call uav_current_center(t, xc, yc, zc)
    Qtot = uav_current_thrust(t)

    Do n = 1, n_uav_markers

       Qn = Qtot * uav_mk_frac(n)
       xp = xc + uav_mk_dx(n);  yp = yc;  zp = zc + uav_mk_dz(n)

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

       sigx = Max(uav_kernel_ncell,1) * dx                    * 0.5d0
       sigz = Max(uav_kernel_ncell,1) * dz                    * 0.5d0
       sigy = Max(uav_kernel_ncell,1) * (yg(j0+1)-yg(j0))     * 0.5d0

       ! Pass 1: normalization so that sum(wgt*dV) over the support equals the marker's thrust share exactly
       norm = 0d0
       Do k = klo, khi
          dzp = zg(k) - zp
          Do j = jlo, jhi
             dVj = dx * dz * ( yg(j+1) - yg(j) )
             dyp = y(j) - yp
             Do i = ilo, ihi
                dxp = xg(i) - xp
                wgt = Exp( -0.5d0*( (dxp/sigx)**2 + (dyp/sigy)**2 + (dzp/sigz)**2 ) )
                norm = norm + wgt*dVj
             End Do
          End Do
       End Do

       If ( norm < 1d-30 ) Cycle   ! kernel support present but entirely between/outside sampled points (degenerate); nothing to add

       ! Pass 2: spread the marker's thrust share as a downward (-y) acceleration density
       Do k = klo, khi
          dzp = zg(k) - zp
          Do j = jlo, jhi
             dyp = y(j) - yp
             Do i = ilo, ihi
                dxp = xg(i) - xp
                wgt = Exp( -0.5d0*( (dxp/sigx)**2 + (dyp/sigy)**2 + (dzp/sigz)**2 ) )
                rhs_v(i,j,k) = rhs_v(i,j,k) - Qn * wgt / norm
             End Do
          End Do
       End Do

    End Do

  End Subroutine apply_uav_forcing

End Module uav_actuator
