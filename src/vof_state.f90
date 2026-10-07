!> Solver-side state and coupling of the geometric PLIC VOF liquid volume fraction C: the padded field, its halo exchange on the
!  decomposed grid, initial interface shapes, the RK3 sub-step advance and the global diagnostics. The kernels (vof_plic,
!  vof_normals, vof_advect) know nothing of MPI or the solver globals; this module adapts them. Host-resident for now: the
!  halo exchange reuses the scalar-transport host path, so a GPU build needs the field moved on-device first.
!
!  Time stepping: C is advanced in three sub-steps (8/15, 2/15, 1/3 dt) that end at the RK stage times, each with the
!  stage-consistent linearly extrapolated midpoint velocity (tests/vof/spike_c_timeint.f90, observed second order):
!    sub-step 1  u^n and the previous step's stage-2 velocity,  midpoint t_n + 4/15 dt
!    sub-step 2  u^n and u^(1),                                 midpoint t_n + 3/5 dt
!    sub-step 3  u^(1) and u^(2),                               midpoint t_n + 5/6 dt
!  Sub-step s must precede the stage-s projection so that C^(s) exists when that projection needs the density.
Module vof_state

  Use iso_fortran_env, Only : Int32, Int64
  Use global
  Use mpi
  Use decomp, Only : x_periodic_partner, z_periodic_partner
  Use halo_pad, Only : pad_field
  Use boundary_conditions, Only : update_ghost_interior_planes_x
  Use scalar_transport, Only : finish_scalar_halos
  Use vof_plic
  Use vof_normals
  Use vof_advect
  Use vof_pressure, Only : vp_init, vp_set_density, vp_selftest, vp_rho, vp_w, vp_iters_last, vp_res_last, vp_its_total
  Use input_output, Only : read_vof_restart
  Use waves, Only : wave_init, wave_eta

  Implicit None

  ! padded C: cell indices match the solver's (1:nxg,1:nyg,1:nzg), the extra layer 0 / n+1 is the second ghost layer
  Real(Int64), Allocatable, Dimension(:,:,:) :: Cv, Cw, Pw, Cinit
  Real(Int64), Allocatable, Dimension(:)     :: hx, hy, hz
  ! transporting velocities: stage 1 and 2 of this step, stage 2 of the previous step, and the extrapolated midpoint field
  Real(Int64), Allocatable, Dimension(:,:,:) :: Us1, Vs1, Ws1, Us2, Vs2, Ws2, Ue, Ve, We

  Integer(Int32) :: vof_unit = 0, gauge_unit = 0
  Logical :: vof_unit_open = .False., gauge_unit_open = .False.   ! newunit numbers are negative, so they cannot flag 'not opened'
  Real(Int64), Allocatable, Dimension(:) :: vof_rho_ref   ! still-water row densities of the initial state (hydrostatic reference)
  Logical        :: vof_have_prev = .False.
  Real(Int64)    :: vof_bnd_cum = 0d0, vof_relax_cum = 0d0   ! cumulative liquid volume through the x boundaries and from relaxation
  Character(*), Parameter :: diag_fmt = '(I9,2ES18.9E3,ES23.14E3,ES13.4E3,2ES12.3E3,I8,2ES12.3E3,3ES18.9E3,ES15.6E3,' // &
                                       'ES18.9E3,3ES14.5E3,I6,ES12.3E3,I4,*(ES16.7E3))'
  Integer(Int32) :: vof_nsub_last = 1
  Integer(Int64) :: vof_hf_fail = 0   ! interface cells without a valid height function in the last curvature evaluation (this rank)
  Integer(Int64) :: vof_its_prev = 0
  Real(Int64)    :: vof_dt_prev = 0d0, vof_vol0 = 0d0, vof_clip_total = 0d0, vof_co_max = 0d0
  Integer(Int32) :: vof_nadv = 0

Contains

  !> Allocate C and the work arrays, set the initial interface and make the halos consistent
  Subroutine vof_init

    Integer(Int32) :: i, j, k

    Real(Int64) :: vliq, cmin, cmax, mom(7)
    Integer(Int64) :: nint

#ifdef GPU_POISSON
    If ( myid == 0 ) Write(*,'(A)') ' ERROR: the VOF / two-fluid solver is host-only, not available in GPU builds'
    Call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
#endif

    Allocate( Cv(0:nxg+1,0:nyg+1,0:nzg+1), Cw(nxg,nyg,nzg), Pw(0:nxg+1,nyg,0:nzg+1) )
    Allocate( hx(nxg), hy(nyg), hz(nzg) )
    Allocate( Us1(nx,nyg,nzg), Vs1(nxg,ny,nzg), Ws1(nxg,nyg,nz) )
    Allocate( Us2(nx,nyg,nzg), Vs2(nxg,ny,nzg), Ws2(nxg,nyg,nz) )
    Allocate( Ue (nx,nyg,nzg), Ve (nxg,ny,nzg), We (nxg,nyg,nz) )

    hx = dx
    Do j = 2, nyg-1
       hy(j) = y(j) - y(j-1)
    End Do
    hy(1) = hy(2);  hy(nyg) = hy(nyg-1)
    Do k = 2, nzg-1
       hz(k) = z(k) - z(k-1)
    End Do
    hz(1) = hz(2);  hz(nzg) = hz(nzg-1)

    Cv = 0d0
    Do k = 2, nzg-1
       Do j = 2, nyg-1
          Do i = 2, nxg-1
             Cv(i,j,k) = initial_fraction( x(i-1), x(i), y(j-1), y(j), z(k-1), z(k) )
          End Do
       End Do
    End Do
    Call vof_fill_pad(Cv, nxg, nyg, nzg)
    Call vof_advect_init(nxg, nyg, nzg)
    vof_flux_scheme = vof_method - 1
    vof_thinc_beta = vof_beta

    If ( z_bc_type == 1 ) Then
       If ( myid == 0 ) Write(*,'(A)') ' ERROR: the VOF field supports periodic z only (z_bc_type=0)'
       Call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    End If
    If ( wave_type > 0 ) Call wave_init
    Allocate( Cvof_io(nxg,nyg,nzg), vof_rho_ref(nyg) )
    If ( vof_flow >= 1 ) Then
       Call vp_init
       Call vof_fill_pad(Cv, nxg, nyg, nzg)
       Call vp_set_density(Cv)
       Call row_reference_density
       If ( vof_selftest == 1 ) Call vp_selftest
    End If

    If ( restart == 1 ) Then
       Call read_vof_restart
       Cv(1:nxg,1:nyg,1:nzg) = Cvof_io
       Call vof_fill_pad(Cv, nxg, nyg, nzg)
    End If
    Cvof_io = Cv(1:nxg,1:nyg,1:nzg)
    If ( vof_prescribed > 0 ) Then
       Allocate( Cinit(nxg,nyg,nzg) )
       Cinit = Cv(1:nxg,1:nyg,1:nzg)
    End If

    Call vof_diagnostics(vliq, cmin, cmax, nint, mom)
    vof_vol0 = vliq
    If ( myid == 0 ) Write(*,'(A,ES14.6,A,ES10.3,A,ES10.3)') &
         '   VOF initial liquid volume = ', vliq, '  C range ', cmin, ' ..', cmax

  End Subroutine vof_init


  !> x-z mean density of every row of the initial state (volume weighted, periodic duplicate cells excluded)
  Subroutine row_reference_density

    Integer(Int32) :: i, j, k
    Real(Int64), Allocatable :: s1(:), s2(:), g1(:), g2(:)

    Allocate( s1(nyg), s2(nyg), g1(nyg), g2(nyg) )
    s1 = 0d0;  s2 = 0d0
    Do k = 2, nzg-1
       Do j = 2, nyg-1
          Do i = 2, nxg-1
             s1(j) = s1(j) + vp_w(i,j,k)*vp_rho(i,j,k)
             s2(j) = s2(j) + vp_w(i,j,k)
          End Do
       End Do
    End Do
    Call MPI_Allreduce(s1, g1, nyg, MPI_real8, MPI_SUM, MPI_COMM_WORLD, ierr)
    Call MPI_Allreduce(s2, g2, nyg, MPI_real8, MPI_SUM, MPI_COMM_WORLD, ierr)
    vof_rho_ref = vof_rho_g
    Do j = 2, nyg-1
       If ( g2(j) > 0d0 ) vof_rho_ref(j) = g1(j)/g2(j)
    End Do
    vof_rho_ref(1) = vof_rho_ref(2);  vof_rho_ref(nyg) = vof_rho_ref(nyg-1)
    Deallocate( s1, s2, g1, g2 )

  End Subroutine row_reference_density


  !> Liquid fraction of the box [x0,x1]x[y0,y1]x[z0,z1] for the configured initial shape
  Function initial_fraction(x0, x1, y0, y1, z0, z1) Result(c)

    Real(Int64), Intent(In) :: x0, x1, y0, y1, z0, z1
    Real(Int64) :: c

    Integer(Int32), Parameter :: ns = 10
    Integer(Int32) :: a, b, d, cnt
    Real(Int64) :: px, py, pz, r2, dmin2, dmax2, qx, qy, qz
    Real(Int64), Parameter :: pi_ = 3.14159265358979323846d0

    If ( vof_ic_type == 1 ) Then
       c = Min(1d0, Max(0d0, (vof_level - y0)/(y1 - y0)))
       Return
    End If
    If ( vof_ic_type == 4 ) Then
       ! standing/sloshing wave: liquid below y = level + amp cos(2 pi x/lambda), column average over 32 sub-columns
       If ( vof_smooth_w > 0d0 ) Then
          ! smoothed Heaviside of the signed vertical distance (interface normal distance ~ dy/sqrt(1+eta_x^2)), width vof_smooth_w cells
          px = 0.5d0*(x0 + x1);  py = 0.5d0*(y0 + y1);  pz = y1 - y0
          qx = vof_level + vof_wave_amp*Cos(2d0*pi_*px/vof_wave_lambda)
          qy = -vof_wave_amp*(2d0*pi_/vof_wave_lambda)*Sin(2d0*pi_*px/vof_wave_lambda)
          qz = (qx - py)/Sqrt(1d0 + qy*qy)/(vof_smooth_w*pz)
          If ( qz <= -1d0 ) Then
             c = 0d0
          Else If ( qz >= 1d0 ) Then
             c = 1d0
          Else
             c = 0.5d0*( 1d0 + qz + Sin(pi_*qz)/pi_ )
          End If
          Return
       End If
       c = 0d0
       Do a = 1, 32
          px = x0 + (x1 - x0)*(Real(a,Int64) - 0.5d0)/32d0
          c = c + Min(1d0, Max(0d0, (vof_level + wave_surface(px) - y0)/(y1 - y0)))
       End Do
       c = c/32d0
       Return
    End If

    If ( vof_ic_type == 5 ) Then
       ! liquid box x < vof_center(1), y < vof_center(2) (dam-break column): exact overlap fraction
       c = Min(1d0, Max(0d0, (vof_center(1) - x0)/(x1 - x0)))
       ! vof_radius > 0: mirror image against the periodic wrap (a column of twice the width centred on the x = 0 plane)
       If ( vof_radius > 0d0 ) c = Min(1d0, c + Min(1d0, Max(0d0, (x1 - (Lx - dx - vof_center(1)))/(x1 - x0))))
       c = c*Min(1d0, Max(0d0, (vof_center(2) - y0)/(y1 - y0)))
       Return
    End If

    If ( vof_ic_type == 6 ) Then
       ! disk in the x-y plane (uniform in z): vof_center(1:2), vof_radius; sampled 20 x 20 per cell
       cnt = 0
       Do b = 1, 20
          py = y0 + (y1 - y0)*(Real(b,Int64) - 0.5d0)/20d0
          Do a = 1, 20
             px = x0 + (x1 - x0)*(Real(a,Int64) - 0.5d0)/20d0
             If ( (px - vof_center(1))**2 + (py - vof_center(2))**2 <= vof_radius**2 ) cnt = cnt + 1
          End Do
       End Do
       c = Real(cnt,Int64)/400d0
       Return
    End If

    ! sphere: classify the box by its nearest and farthest point first, sub-sample only the cut boxes
    r2 = vof_radius*vof_radius
    qx = Min(Max(vof_center(1), x0), x1) - vof_center(1)
    qy = Min(Max(vof_center(2), y0), y1) - vof_center(2)
    qz = Min(Max(vof_center(3), z0), z1) - vof_center(3)
    dmin2 = qx*qx + qy*qy + qz*qz
    qx = Max(Abs(x0 - vof_center(1)), Abs(x1 - vof_center(1)))
    qy = Max(Abs(y0 - vof_center(2)), Abs(y1 - vof_center(2)))
    qz = Max(Abs(z0 - vof_center(3)), Abs(z1 - vof_center(3)))
    dmax2 = qx*qx + qy*qy + qz*qz
    If ( dmax2 <= r2 ) Then
       c = 1d0
    Else If ( dmin2 >= r2 ) Then
       c = 0d0
    Else
       cnt = 0
       Do d = 1, ns
          pz = z0 + (z1 - z0)*(Real(d,Int64) - 0.5d0)/ns
          Do b = 1, ns
             py = y0 + (y1 - y0)*(Real(b,Int64) - 0.5d0)/ns
             Do a = 1, ns
                px = x0 + (x1 - x0)*(Real(a,Int64) - 0.5d0)/ns
                If ( (px - vof_center(1))**2 + (py - vof_center(2))**2 + (pz - vof_center(3))**2 <= r2 ) cnt = cnt + 1
             End Do
          End Do
       End Do
       c = Real(cnt,Int64)/Real(ns*ns*ns,Int64)
    End If
    If ( vof_ic_type == 3 ) c = 1d0 - c

  End Function initial_fraction


  !> Initial surface elevation of the vof_ic_type=4 wave: cosine, plus the second-order Stokes harmonic if vof_wave_stokes = 1
  Pure Function wave_surface(xx) Result(eta)

    Real(Int64), Intent(In) :: xx
    Real(Int64) :: eta, kk, pi_

    pi_ = 4d0*Atan(1d0)
    kk = 2d0*pi_/vof_wave_lambda
    eta = vof_wave_amp*Cos(kk*xx)
    If ( vof_wave_stokes == 1 ) eta = eta + 0.5d0*kk*vof_wave_amp**2*Cos(2d0*kk*xx)

  End Function wave_surface


  !> Refresh both ghost layers of the padded C: rank seams and periodic wraps through the scalar-transport host path, zero-gradient
  !  at x inflow/outflow and y walls (periodic in y when y_bc_type=0). Matches vof_advect's fill_pad_iface.
  Subroutine vof_fill_pad(Cp, n1, n2, n3)

    Integer(Int32), Intent(In)    :: n1, n2, n3
    Real(Int64),    Intent(InOut) :: Cp(0:n1+1,0:n2+1,0:n3+1)

    Logical :: is_first, is_last
    Integer(Int32) :: partner

    Cw = Cp(1:n1,1:n2,1:n3)
    Call update_ghost_interior_planes_x(Cw, 4)
    If ( x_bc_type == 1 ) Then
       Call x_periodic_partner(is_first, is_last, partner)
       If ( is_first ) Then
          Cw(1,:,:) = Cw(2,:,:)
          If ( inflow_type == 3 .And. wave_type > 0 ) Call wave_inlet_fraction
       End If
       If ( is_last  ) Cw(n1,:,:) = Cw(n1-1,:,:)
    End If
    Call finish_scalar_halos(Cw)
    If ( y_bc_type == 0 ) Then
       Cw(:,1,:)    = Cw(:,n2-2,:)
       Cw(:,n2-1,:) = Cw(:,2,:)
       Cw(:,n2,:)   = Cw(:,3,:)
    Else
       Cw(:,1,:)  = Cw(:,2,:)
       Cw(:,n2,:) = Cw(:,n2-1,:)
    End If

    Call pad_field(Cw, n1, n2, n3, .False., .False., 1, Pw)
    Cp(0:n1+1,1:n2,0:n3+1) = Pw
    If ( y_bc_type == 0 ) Then
       Cp(:,0,:)    = Cp(:,n2-3,:)
       Cp(:,n2+1,:) = Cp(:,4,:)
    Else
       Cp(:,0,:)    = Cp(:,3,:)
       Cp(:,n2+1,:) = Cp(:,n2-2,:)
    End If

  End Subroutine vof_fill_pad


  !> Inlet ghost cell of C: liquid below the target wave surface at the ghost-cell centre
  Subroutine wave_inlet_fraction

    Integer(Int32) :: j
    Real(Int64) :: top

    top = vof_level + wave_eta(xg(1), t)
    Do j = 2, nyg-1
       Cw(1,j,:) = Min(1d0, Max(0d0, (top - y(j-1))/(y(j) - y(j-1))))
    End Do

  End Subroutine wave_inlet_fraction


  !> Store the projected stage velocity (stage 1 or 2) used by the later sub-steps' midpoint extrapolation
  Subroutine vof_save_stage(s)

    Integer(Int32), Intent(In) :: s

    !$acc update host(U,V,W)
    If ( s == 1 ) Then
       Us1 = U;  Vs1 = V;  Ws1 = W
    Else
       Us2 = U;  Vs2 = V;  Ws2 = W
    End If

  End Subroutine vof_save_stage


  !> Sub-step s (1..3) of the C advance, ending at the stage time rk_t(s); call before the stage-s projection
  Subroutine vof_advance_substep(s)

    Integer(Int32), Intent(In) :: s

    Real(Int64) :: dts, a, cl, co

    !$acc update host(U,V,W,Uo,Vo,Wo)
    Select Case(s)
    Case(1)
       dts = dt*rk_t(1)
       If ( vof_have_prev ) Then
          a = Min(0.8d0*dt/vof_dt_prev, 2d0)
          Ue = Uo + a*(Uo - Us2);  Ve = Vo + a*(Vo - Vs2);  We = Wo + a*(Wo - Ws2)
       Else
          Ue = Uo;  Ve = Vo;  We = Wo
       End If
    Case(2)
       dts = dt*(rk_t(2) - rk_t(1))
       a = (0.6d0 - rk_t(1))/rk_t(1)
       Ue = Us1 + a*(Us1 - Uo);  Ve = Vs1 + a*(Vs1 - Vo);  We = Ws1 + a*(Ws1 - Wo)
    Case Default
       dts = dt*(rk_t(3) - rk_t(2))
       a = (5d0/6d0 - rk_t(2))/(rk_t(2) - rk_t(1))
       Ue = Us2 + a*(Us2 - Us1);  Ve = Vs2 + a*(Vs2 - Vs1);  We = Ws2 + a*(Ws2 - Ws1)
    End Select

    If ( vof_prescribed > 0 ) Call prescribed_velocity(t - 0.5d0*dts, Ue, Ve, We)
    Call vof_advect_step(Cv, nxg, nyg, nzg, Ue, Ve, We, hx, hy, hz, dts, vof_nadv, vof_normal_scheme, vof_fill_pad, co, cl)
    vof_nadv = vof_nadv + 1
    vof_co_max = Max(vof_co_max, co)
    vof_clip_total = vof_clip_total + cl

  End Subroutine vof_advance_substep


  !> Analytic divergence-free test velocity at time tm on the staggered faces, evaluated as the discrete curl of a vector potential
  !  at the cell edges so that the discrete divergence is zero to round-off: 1 LeVeque vortex reversal, 2 Enright deformation
  Subroutine prescribed_velocity(tm, Uf, Vf, Wf)

    Real(Int64), Intent(In)  :: tm
    Real(Int64), Intent(Out) :: Uf(nx,nyg,nzg), Vf(nxg,ny,nzg), Wf(nxg,nyg,nz)
    Integer(Int32) :: i, j, k
    Real(Int64) :: g

    g = Cos(4d0*Atan(1d0)*tm/vof_presc_T)
    Uf = 0d0;  Vf = 0d0;  Wf = 0d0
    Do k = 2, nzg-1
       Do j = 2, nyg-1
          Do i = 1, nx
             Uf(i,j,k) = g*( ( pot3(x(i), y(j), zg(k)) - pot3(x(i), y(j-1), zg(k)) )/hy(j) &
                           - ( pot2(x(i), yg(j), z(k)) - pot2(x(i), yg(j), z(k-1)) )/hz(k) )
          End Do
       End Do
    End Do
    Do k = 2, nzg-1
       Do j = 1, ny
          Do i = 2, nxg-1
             Vf(i,j,k) = -g*( pot3(x(i), y(j), zg(k)) - pot3(x(i-1), y(j), zg(k)) )/dx
          End Do
       End Do
    End Do
    If ( vof_prescribed == 2 ) Then
       Do k = 1, nz
          Do j = 2, nyg-1
             Do i = 2, nxg-1
                Wf(i,j,k) = g*( pot2(x(i), yg(j), z(k)) - pot2(x(i-1), yg(j), z(k)) )/dx
             End Do
          End Do
       End Do
    End If

  End Subroutine prescribed_velocity


  !> z-component of the vector potential: vortex-reversal stream function / Enright field
  Pure Function pot3(xx, yy, zz) Result(a)
    Real(Int64), Intent(In) :: xx, yy, zz
    Real(Int64) :: a, pi
    pi = 4d0*Atan(1d0)
    If ( vof_prescribed == 1 ) Then
       a = Sin(pi*xx)**2*Sin(pi*yy)**2/pi
    Else
       a = -Cos(2d0*pi*xx)/(2d0*pi)*Sin(pi*yy)**2*Sin(2d0*pi*zz)
    End If
  End Function pot3


  !> y-component of the vector potential (Enright field only)
  Pure Function pot2(xx, yy, zz) Result(a)
    Real(Int64), Intent(In) :: xx, yy, zz
    Real(Int64) :: a, pi
    pi = 4d0*Atan(1d0)
    a = 0d0
    If ( vof_prescribed == 2 ) a = Cos(2d0*pi*xx)/(2d0*pi)*Sin(2d0*pi*yy)*Sin(pi*zz)**2 + Sin(2d0*pi*yy)*Cos(2d0*pi*zz)/(2d0*pi)
  End Function pot2


  !> End of step: remember stage 2 as the previous-step velocity for the next step's first sub-step
  Subroutine vof_end_step

    Cvof_io = Cv(1:nxg,1:nyg,1:nzg)
    vof_dt_prev = dt
    vof_have_prev = .True.

  End Subroutine vof_end_step


  !> Append one row to vof_diag.dat (rank 0): time, liquid volume, relative drift against the initial volume, C range, interface
  !  cells, cumulative clip loss, largest sub-step Courant number seen since the last row
  Subroutine vof_output_monitor

    Real(Int64) :: vliq, cmin, cmax, mom(7), vv, vel_loc(3), vel_glb(3), xloc, yloc, serr, pjump
    Integer(Int64) :: nfail
    Integer(Int64) :: nint

    Call vof_diagnostics(vliq, cmin, cmax, nint, mom)
    serr = shape_error()
    pjump = pressure_jump()
    Call MPI_Allreduce(vof_hf_fail, nfail, 1, MPI_integer8, MPI_SUM, MPI_COMM_WORLD, ierr)
    vv = Max(vliq, 1d-300)
    vel_loc = (/ MaxVal(Abs(U(2:nx-1,2:nyg-1,2:nzg-1))), MaxVal(Abs(V(2:nxg-1,2:ny-1,2:nzg-1))), &
                 MaxVal(Abs(W(2:nxg-1,2:nyg-1,2:nz-1))) /)
    Call MPI_Allreduce(vel_loc, vel_glb, 3, MPI_real8, MPI_MAX, MPI_COMM_WORLD, ierr)
    ! where the fastest |U| is (x, y of the face), reduced with MAXLOC
    Block
      Integer(Int32) :: ia(3)
      Real(Int64) :: pin(2), pout(2)
      ia = MaxLoc(Abs(U(2:nx-1,2:nyg-1,2:nzg-1)))
      pin(1) = MaxVal(Abs(U(2:nx-1,2:nyg-1,2:nzg-1)))
      pin(2) = Real(myid,Int64)
      Call MPI_Allreduce(pin, pout, 1, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, MPI_COMM_WORLD, ierr)
      xloc = 0d0;  yloc = 0d0
      If ( Real(myid,Int64) == pout(2) ) Then
         xloc = x(ia(1)+1);  yloc = yg(ia(2)+1)
      End If
      pin = (/ xloc, yloc /)
      Call MPI_Allreduce(pin, pout, 2, MPI_real8, MPI_SUM, MPI_COMM_WORLD, ierr)
      xloc = pout(1);  yloc = pout(2)
    End Block
    If ( myid == 0 ) Then
       If ( .Not. vof_unit_open ) Then
          vof_unit_open = .True.
          Open(newunit=vof_unit, file='vof_diag.dat', status='unknown', position='append', action='write')
          Write(vof_unit,'(A)') '# 1 step  2 t  3 dt  4 Vliq  5 (Vliq-V0)/V0  6 Cmin  7 Cmax  8 n_interface  9 clip_loss_cum  ' // &
               '10 Co_adv_max  11 xc  12 yc  13 zc  14 int C(1-C)  15 cos-moment  16 max|U|  17 max|V|  18 max|W|  ' // &
               '19 PCG_its_step  20 last_PCG_res  21 n_sub  22 bnd_cum  23 relax_cum  24 x_of_max|U|  25 y_of_max|U|  ' // &
               '26 KE_liquid  27 KE_gas  28 L1(C-C_init)  29 p_liquid-p_gas  30 HF_fallback_cells'
       End If
       Write(vof_unit, diag_fmt) &
            istep, t, dt, vliq, (vliq - vof_vol0)/Max(vof_vol0, 1d-300), cmin, cmax, Int(nint), vof_clip_total, vof_co_max, &
            mom(1)/vv, mom(2)/vv, mom(3)/vv, mom(4), mom(5), vel_glb, Int(vp_its_total - vof_its_prev), vp_res_last, &
            vof_nsub_last, vof_bnd_cum, vof_relax_cum, xloc, yloc, mom(6), mom(7), serr, pjump, Real(nfail, Int64)
       Flush(vof_unit)
    End If
    vof_its_prev = vp_its_total
    If ( vof_debug == 1 .And. vof_snap_dt > 0d0 .And. nprocs == 1 ) Call write_snapshot
    vof_co_max = 0d0
    ! diagnostic (vof_debug = 1): deviation of the velocity from the uniform start-up stream vof_u0
    If ( vof_debug == 1 ) Then
    Block
      Real(Int64) :: dl(3), dg(3)
      Integer, Save :: udv = 0
      Logical, Save :: udv_open = .False.
      Integer(Int32) :: ia(3) = 0
      dl = (/ MaxVal(Abs(U(2:nx-1,2:nyg-1,2:nzg-1) - vof_u0)), MaxVal(Abs(V(2:nxg-1,2:ny-1,2:nzg-1))), &
              MaxVal(Abs(W(2:nxg-1,2:nyg-1,2:nz-1))) /)
      Call MPI_Allreduce(dl, dg, 3, MPI_real8, MPI_MAX, MPI_COMM_WORLD, ierr)
      If ( nprocs == 1 ) Then
         ia = MaxLoc(Abs(U(2:nx-1,2:nyg-1,2:nzg-1) - vof_u0))
      End If
      If ( myid == 0 ) Then
         If ( .Not. udv_open ) Then
            Open(newunit=udv, file='vof_dev.dat', status='unknown', action='write')
            udv_open = .True.
         End If
         Write(udv,'(I8,4ES16.8,3I6)') istep, t, dg, ia(1)+1, ia(2)+1, ia(3)+1
         Flush(udv)
      End If
    End Block
    End If
    ! diagnostic (vof_debug = 1, single rank): where the fastest |V| is and the phase fractions above and below that face
    If ( vof_debug == 1 .And. nprocs == 1 ) Then
       Block
          Integer(Int32) :: ia(3)
          Integer, Save :: uvm = 0
          Logical, Save :: uvm_open = .False.
          ia = MaxLoc(Abs(V(2:nxg-1,2:ny-1,2:nzg-1))) + 1
          If ( .Not. uvm_open ) Then
             Open(newunit=uvm, file='vof_vmax.dat', status='unknown', action='write')
             uvm_open = .True.
             Write(uvm,'(A)') '# step t V(max|V|) i j k  x y  C_below C_above  V_neighbours(j-1,j+1) U_at_cell'
          End If
          Write(uvm,'(I8,ES14.6,ES14.6,3I6,2F9.4,2ES11.3,3ES13.5)') istep, t, V(ia(1),ia(2),ia(3)), ia, xg(ia(1)), y(ia(2)), &
               Cv(ia(1),ia(2),ia(3)), Cv(ia(1),ia(2)+1,ia(3)), V(ia(1),ia(2)-1,ia(3)), V(ia(1),ia(2)+1,ia(3)), &
               0.5d0*(U(ia(1),ia(2),ia(3)) + U(ia(1)-1,ia(2),ia(3)))
          Flush(uvm)
       End Block
    End If
    ! diagnostic (vof_debug = 1): leading edge of the liquid in the bottom cell row inside the first half period (dam break)
    If ( vof_debug == 1 ) Then
       Block
          Real(Int64) :: xf, xfg
          Integer, Save :: ufr = 0
          Logical, Save :: ufr_open = .False.
          Integer(Int32) :: i, k
          xf = 0d0
          Do k = 2, nzg-1
             Do i = 2, nxg-1
                If ( Cv(i,2,k) > 0.5d0 .And. xg(i) < 0.5d0*(Lx - dx) ) xf = Max(xf, xg(i) + 0.5d0*hx(i))
             End Do
          End Do
          Call MPI_Allreduce(xf, xfg, 1, MPI_real8, MPI_MAX, MPI_COMM_WORLD, ierr)
          If ( myid == 0 ) Then
             If ( .Not. ufr_open ) Then
                Open(newunit=ufr, file='vof_front.dat', status='unknown', action='write')
                ufr_open = .True.
                Write(ufr,'(A)') '# t  x_front (liquid edge in the bottom row, x < half period)'
             End If
             Write(ufr,'(2ES16.8)') t, xfg
             Flush(ufr)
          End If
       End Block
    End If
    ! diagnostic (vof_debug = 1, single rank): u(y), C(y), v(y) through the column at x ~ Lx/4 every 5 steps
    If ( vof_debug == 1 .And. nprocs == 1 .And. Mod(istep, 5) == 0 ) Then
       Block
          Integer(Int32) :: ic, jj
          Integer, Save :: upr = 0
          Logical, Save :: upr_open = .False.
          ic = 2
          Do jj = 2, nx-1
             If ( Abs(x(jj) - 0.25d0) < Abs(x(ic) - 0.25d0) ) ic = jj
          End Do
          If ( .Not. upr_open ) Then
             Open(newunit=upr, file='vof_prof.dat', status='unknown', action='write')
             upr_open = .True.
          End If
          Write(upr,'(A,ES14.6,A,I6,A,F8.4)') '# t=', t, ' step=', istep, ' x=', x(ic)
          Do jj = 2, nyg-1
             Write(upr,'(F10.5,2ES16.8,ES16.8)') yg(jj), U(ic,jj,3), 0.5d0*(Cv(ic,jj,3)+Cv(ic+1,jj,3)), V(ic,Min(jj,ny),3)
          End Do
          Flush(upr)
       End Block
    End If

  End Subroutine vof_output_monitor


  !> Debug snapshot (one rank): x, y, C and the cell-centred u, v of the middle z plane
  Subroutine write_snapshot

    Integer(Int32) :: i, j, k, iu
    Integer, Save :: nsnap = 0
    Real(Int64), Save :: tnext = 0d0
    Character(len=32) :: fn

    If ( t + 1d-12 < tnext ) Return
    nsnap = nsnap + 1
    tnext = tnext + vof_snap_dt
    k = nzg/2
    Write(fn,'(A,I5.5,A)') 'vof_snap_', nsnap, '.dat'
    Open(newunit=iu, file=Trim(fn), status='replace', action='write')
    Write(iu,'(A,ES14.6)') '# t = ', t
    Do j = 2, nyg-1
       Do i = 2, nxg-1
          Write(iu,'(5ES14.6)') xg(i), yg(j), Cv(i,j,k), 0.5d0*(U(i,j,k) + U(i-1,j,k)), 0.5d0*(V(i,j,k) + V(i,j-1,k))
       End Do
    End Do
    Close(iu)

  End Subroutine write_snapshot


  !> Volume integral of |C - C_init| over the distinct cells (vof_prescribed runs: the reversal test returns to the initial shape)
  Function shape_error() Result(e)

    Real(Int64) :: e, el
    Integer(Int32) :: i, j, k, ihi, khi
    Logical :: is_first, is_last
    Integer(Int32) :: partner

    e = 0d0
    If ( .Not. Allocated(Cinit) ) Return
    ihi = nxg-1;  khi = nzg-1
    If ( x_bc_type == 0 ) Then
       Call x_periodic_partner(is_first, is_last, partner)
       If ( is_last ) ihi = nxg-2
    End If
    If ( z_bc_type == 0 ) Then
       Call z_periodic_partner(is_first, is_last, partner)
       If ( is_last ) khi = nzg-2
    End If
    el = 0d0
    Do k = 2, khi
       Do j = 2, nyg-1
          Do i = 2, ihi
             el = el + Abs(Cv(i,j,k) - Cinit(i,j,k))*hx(i)*hy(j)*hz(k)
          End Do
       End Do
    End Do
    Call MPI_Allreduce(el, e, 1, MPI_real8, MPI_SUM, MPI_COMM_WORLD, ierr)

  End Function shape_error


  !> Volume-mean pressure of the liquid cells (C > 0.99) minus that of the gas cells (C < 0.01): the Laplace jump of a static drop
  Function pressure_jump() Result(pj)

    Real(Int64) :: pj, loc(4), glb(4), vc
    Integer(Int32) :: i, j, k, ihi, khi
    Logical :: is_first, is_last
    Integer(Int32) :: partner

    ihi = nxg-1;  khi = nzg-1
    If ( x_bc_type == 0 ) Then
       Call x_periodic_partner(is_first, is_last, partner)
       If ( is_last ) ihi = nxg-2
    End If
    If ( z_bc_type == 0 ) Then
       Call z_periodic_partner(is_first, is_last, partner)
       If ( is_last ) khi = nzg-2
    End If
    loc = 0d0
    Do k = 2, khi
       Do j = 2, nyg-1
          Do i = 2, ihi
             vc = hx(i)*hy(j)*hz(k)
             If ( Cv(i,j,k) > 0.99d0 ) Then
                loc(1) = loc(1) + vc*P(i,j,k);  loc(2) = loc(2) + vc
             Else If ( Cv(i,j,k) < 0.01d0 ) Then
                loc(3) = loc(3) + vc*P(i,j,k);  loc(4) = loc(4) + vc
             End If
          End Do
       End Do
    End Do
    Call MPI_Allreduce(loc, glb, 4, MPI_real8, MPI_SUM, MPI_COMM_WORLD, ierr)
    pj = 0d0
    If ( glb(2) > 0d0 .And. glb(4) > 0d0 ) pj = glb(1)/glb(2) - glb(3)/glb(4)

  End Function pressure_jump


  !> Surface elevation at the wave gauges (column liquid height, mean over z), appended to vof_gauges.dat
  Subroutine vof_write_gauges

    Integer(Int32) :: ig, i, j, k, ng, ihi, khi
    Real(Int64) :: loc(8), glb(8), cnt_loc(8), cnt(8), col
    Logical :: is_first, is_last
    Integer(Int32) :: partner

    ng = Count(wave_gauge_x >= 0d0)
    If ( ng == 0 ) Return
    ihi = nxg-1;  khi = nzg-1
    If ( z_bc_type == 0 ) Then
       Call z_periodic_partner(is_first, is_last, partner)
       If ( is_last ) khi = nzg-2
    End If
    loc = 0d0;  cnt_loc = 0d0
    Do ig = 1, ng
       Do i = 2, ihi
          ! owner of the gauge: the cell whose faces bracket it (cells are [x(i-1), x(i)])
          If ( wave_gauge_x(ig) >= x(i-1) .And. wave_gauge_x(ig) < x(i) ) Then
             Do k = 2, khi
                col = 0d0
                Do j = 2, nyg-1
                   col = col + Cv(i,j,k)*hy(j)
                End Do
                loc(ig) = loc(ig) + col - vof_level
                cnt_loc(ig) = cnt_loc(ig) + 1d0
             End Do
          End If
       End Do
    End Do
    Call MPI_Allreduce(loc, glb, 8, MPI_real8, MPI_SUM, MPI_COMM_WORLD, ierr)
    Call MPI_Allreduce(cnt_loc, cnt, 8, MPI_real8, MPI_SUM, MPI_COMM_WORLD, ierr)
    If ( myid == 0 ) Then
       If ( .Not. gauge_unit_open ) Then
          gauge_unit_open = .True.
          Open(newunit=gauge_unit, file='vof_gauges.dat', status='unknown', position='append', action='write')
          Write(gauge_unit,'(A,8ES12.4)') '# t eta(x_g) at x =', wave_gauge_x(1:ng)
       End If
       Write(gauge_unit,'(ES16.8,8ES14.6)') t, (glb(ig)/Max(cnt(ig), 1d0), ig = 1, ng)
       Flush(gauge_unit)
    End If

  End Subroutine vof_write_gauges


  !> Global liquid volume, extremes of C and interface-cell count over the owned cells (periodic duplicate cells left out), plus
  !  mom = (first moments of the liquid in x, y, z, integral of C(1-C) as a smearing measure)
  Subroutine vof_diagnostics(vliq, cmin, cmax, nint, mom)

    Real(Int64),    Intent(Out) :: vliq, cmin, cmax, mom(7)
    Integer(Int64), Intent(Out) :: nint

    Integer(Int32) :: i, j, k, ihi, jhi, khi
    Logical :: is_first, is_last
    Integer(Int32) :: partner
    Real(Int64) :: c, vc, buf(8), gbuf(8), mn, mx, kd, ke
    Integer(Int64) :: ni, gni

    ihi = nxg-1;  jhi = nyg-1;  khi = nzg-1
    If ( x_bc_type == 0 ) Then
       Call x_periodic_partner(is_first, is_last, partner)
       If ( is_last ) ihi = nxg-2
    End If
    If ( z_bc_type == 0 ) Then
       Call z_periodic_partner(is_first, is_last, partner)
       If ( is_last ) khi = nzg-2
    End If
    If ( y_bc_type == 0 ) jhi = nyg-2

    buf = 0d0;  mn = 1d300;  mx = -1d300;  ni = 0
    kd = 0d0
    If ( vof_wave_lambda > 0d0 ) kd = 8d0*Atan(1d0)/vof_wave_lambda
    Do k = 2, khi
       Do j = 2, jhi
          Do i = 2, ihi
             c = Cv(i,j,k)
             vc = hx(i)*hy(j)*hz(k)
             buf(1) = buf(1) + c*vc
             buf(2) = buf(2) + c*vc*xg(i)
             buf(3) = buf(3) + c*vc*yg(j)
             buf(4) = buf(4) + c*vc*zg(k)
             buf(5) = buf(5) + c*(1d0 - c)*vc
             buf(6) = buf(6) + c*vc*Cos(kd*xg(i))
             ke = 0.125d0*( (U(i,j,k) + U(i-1,j,k))**2 + (V(i,j,k) + V(i,j-1,k))**2 + (W(i,j,k) + W(i,j,k-1))**2 )*vc
             buf(7) = buf(7) + c*vof_rho_l*ke
             buf(8) = buf(8) + (1d0 - c)*vof_rho_g*ke
             mn = Min(mn, c);  mx = Max(mx, c)
             If ( c > vof_eps .And. c < 1d0 - vof_eps ) ni = ni + 1
          End Do
       End Do
    End Do

    Call MPI_Allreduce(buf, gbuf, 8, MPI_real8, MPI_SUM, MPI_COMM_WORLD, ierr)
    Call MPI_Allreduce(mn, cmin, 1, MPI_real8, MPI_MIN, MPI_COMM_WORLD, ierr)
    Call MPI_Allreduce(mx, cmax, 1, MPI_real8, MPI_MAX, MPI_COMM_WORLD, ierr)
    Call MPI_Allreduce(ni, gni, 1, MPI_integer8, MPI_SUM, MPI_COMM_WORLD, ierr)
    vliq = gbuf(1);  mom = gbuf(2:8);  nint = gni

  End Subroutine vof_diagnostics

End Module vof_state
