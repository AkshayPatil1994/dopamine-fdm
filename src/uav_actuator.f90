!> UAV actuator-disk body force (Phase 1: static disk, uniform loading,
!  vertical-only reaction force -- see docs/UAV_ActuatorDisk_Design.md).
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
!  Known Phase 1 simplifications (see design doc for later phases):
!   - the disk does not move (uav_xc/yc/zc are fixed for the whole run)
!   - only the vertical (rhs_v) reaction force is applied -- correct for a
!     level hovering/climbing/descending disk, not yet for a tilted one
!   - kernel support is not wrapped across periodic x/z boundaries; keep the
!     disk at least uav_kernel_ncell*max(dx,dz) away from a periodic edge
!   - no OpenACC offload yet (CPU builds only; GPU_POISON builds still run
!     this module on the host)
Module uav_actuator

  Use iso_fortran_env, Only : Int32, Int64
  Use global
  Use mpi, Only : myid

  Implicit None

  Integer(Int32) :: n_uav_markers = 0
  Real   (Int64), Allocatable :: uav_mk_x(:), uav_mk_y(:), uav_mk_z(:)
  Real   (Int64), Allocatable :: uav_mk_Q(:)   ! kinematic thrust share per marker [m^4/s^2]

Contains

  !> Build the marker ring (uav_n_r radial bands x uav_n_theta azimuthal
  !  sectors), each carrying a thrust share proportional to its annulus-
  !  sector area (uniform disk loading). Must be called after uav_active
  !  is known (post namelist read/broadcast); no grid dependency.
  Subroutine setup_uav

    Integer(Int32) :: ir, ith, n
    Real(Int64) :: dr, dtheta, r_mid, theta, area_disk, area_k

    n_uav_markers = uav_n_r * uav_n_theta
    Allocate( uav_mk_x(n_uav_markers), uav_mk_y(n_uav_markers), &
              uav_mk_z(n_uav_markers), uav_mk_Q(n_uav_markers) )

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
          uav_mk_x(n) = uav_xc + r_mid*Cos(theta)
          uav_mk_y(n) = uav_yc
          uav_mk_z(n) = uav_zc + r_mid*Sin(theta)
          uav_mk_Q(n) = uav_hover_thrust * (area_k/area_disk)
       End Do
    End Do

    If ( myid==0 ) Then
       Write(*,*) 'UAV: actuator disk active -- markers:', n_uav_markers, &
                  ' radius:', uav_disk_radius, ' centre:', uav_xc, uav_yc, uav_zc, &
                  ' kinematic thrust:', uav_hover_thrust
    End If

  End Subroutine setup_uav

  !> Add the disk's vertical reaction force (downward on the fluid, opposing
  !  the upward thrust that supports the vehicle) into rhs_v. Called once per
  !  RK3 sub-stage from compute_rhs_v (equations.f90), guarded by uav_active.
  Subroutine apply_uav_forcing(rhs_v)

    Real(Int64), Dimension(2:nxg-1,2:ny-1,2:nzg-1), Intent(InOut) :: rhs_v

    Integer(Int32) :: n, i, j, k, i0, j0, k0
    Integer(Int32) :: ilo, ihi, jlo, jhi, klo, khi
    Real(Int64) :: xp, yp, zp
    Real(Int64) :: sigx, sigz, sigy, dxp, dyp, dzp, w, norm, dVj
    Real(Int64) :: best

    Do n = 1, n_uav_markers

       xp = uav_mk_x(n);  yp = uav_mk_y(n);  zp = uav_mk_z(n)

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

       ! Pass 1: normalization so that sum(w*dV) over the support equals the marker's thrust share exactly
       norm = 0d0
       Do k = klo, khi
          dzp = zg(k) - zp
          Do j = jlo, jhi
             dVj = dx * dz * ( yg(j+1) - yg(j) )
             dyp = y(j) - yp
             Do i = ilo, ihi
                dxp = xg(i) - xp
                w = Exp( -0.5d0*( (dxp/sigx)**2 + (dyp/sigy)**2 + (dzp/sigz)**2 ) )
                norm = norm + w*dVj
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
                w = Exp( -0.5d0*( (dxp/sigx)**2 + (dyp/sigy)**2 + (dzp/sigz)**2 ) )
                rhs_v(i,j,k) = rhs_v(i,j,k) - uav_mk_Q(n) * w / norm
             End Do
          End Do
       End Do

    End Do

  End Subroutine apply_uav_forcing

End Module uav_actuator
