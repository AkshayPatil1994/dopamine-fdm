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
  Use vof_pressure, Only : vp_init, vp_set_density, vp_selftest, vp_rho, vp_w, vp_iters_last, vp_res_last
  Use input_output, Only : read_vof_restart
  Use waves, Only : wave_init, wave_eta

  Implicit None

  ! padded C: cell indices match the solver's (1:nxg,1:nyg,1:nzg), the extra layer 0 / n+1 is the second ghost layer
  Real(Int64), Allocatable, Dimension(:,:,:) :: Cv, Cw, Pw
  Real(Int64), Allocatable, Dimension(:)     :: hx, hy, hz
  ! transporting velocities: stage 1 and 2 of this step, stage 2 of the previous step, and the extrapolated midpoint field
  Real(Int64), Allocatable, Dimension(:,:,:) :: Us1, Vs1, Ws1, Us2, Vs2, Ws2, Ue, Ve, We

  Integer(Int32) :: vof_unit = 0, gauge_unit = 0
  Logical :: vof_unit_open = .False., gauge_unit_open = .False.   ! newunit numbers are negative, so they cannot flag 'not opened'
  Real(Int64), Allocatable, Dimension(:) :: vof_rho_ref   ! still-water row densities of the initial state (hydrostatic reference)
  Logical        :: vof_have_prev = .False.
  Real(Int64)    :: vof_bnd_cum = 0d0, vof_relax_cum = 0d0   ! cumulative liquid volume through the x boundaries and from relaxation
  Real(Int64)    :: vof_dt_prev = 0d0, vof_vol0 = 0d0, vof_clip_total = 0d0, vof_co_max = 0d0
  Integer(Int32) :: vof_nadv = 0

Contains

  !> Allocate C and the work arrays, set the initial interface and make the halos consistent
  Subroutine vof_init

    Integer(Int32) :: i, j, k
    Real(Int64) :: vliq, cmin, cmax, mom(5)
    Integer(Int64) :: nint

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

    If ( z_bc_type == 1 ) Then
       If ( myid == 0 ) Write(*,'(A)') ' ERROR: the VOF field supports periodic z only (z_bc_type=0)'
       Call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    End If
    If ( wave_type > 0 ) Call wave_init
    Allocate( Cvof_io(nxg,nyg,nzg), vof_rho_ref(nyg) )
    If ( vof_flow >= 1 ) Then
       Call vp_init
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
       c = 0d0
       Do a = 1, 32
          px = x0 + (x1 - x0)*(Real(a,Int64) - 0.5d0)/32d0
          c = c + Min(1d0, Max(0d0, (vof_level + vof_wave_amp*Cos(2d0*pi_*px/vof_wave_lambda) - y0)/(y1 - y0)))
       End Do
       c = c/32d0
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

    Call vof_advect_step(Cv, nxg, nyg, nzg, Ue, Ve, We, hx, hy, hz, dts, vof_nadv, vof_normal_scheme, vof_fill_pad, co, cl)
    vof_nadv = vof_nadv + 1
    vof_co_max = Max(vof_co_max, co)
    vof_clip_total = vof_clip_total + cl

  End Subroutine vof_advance_substep


  !> End of step: remember stage 2 as the previous-step velocity for the next step's first sub-step
  Subroutine vof_end_step

    Cvof_io = Cv(1:nxg,1:nyg,1:nzg)
    vof_dt_prev = dt
    vof_have_prev = .True.

  End Subroutine vof_end_step


  !> Append one row to vof_diag.dat (rank 0): time, liquid volume, relative drift against the initial volume, C range, interface
  !  cells, cumulative clip loss, largest sub-step Courant number seen since the last row
  Subroutine vof_output_monitor

    Real(Int64) :: vliq, cmin, cmax, mom(5), vv, vel_loc(3), vel_glb(3), xloc, yloc
    Integer(Int64) :: nint

    Call vof_diagnostics(vliq, cmin, cmax, nint, mom)
    vv = Max(vliq, 1d-300)
    vel_loc = (/ -MinVal(U(2:nx-1,2:nyg-1,2:nzg-1)), MaxVal(Abs(V(2:nxg-1,2:ny-1,2:nzg-1))), &
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
          Write(vof_unit,'(A)') '# step t Vliq (Vliq-V0)/V0 Cmin Cmax n_interface clip_loss Co_max xc yc zc ' // &
               ' int C(1-C) cos-moment -Umin |V|max |W|max last-PCG-its last-PCG-res bnd_cum relax_cum x_Umax y_Umax'
       End If
       Write(vof_unit,'(I10,ES18.10,ES22.14,3ES14.5,I10,2ES14.5,3ES20.12,ES16.8,ES20.12,3ES14.5,I6,ES11.3,2ES16.8,2ES11.3)') &
            istep, t, vliq, &
            (vliq - vof_vol0)/Max(vof_vol0, 1d-300), cmin, cmax, Int(nint), vof_clip_total, vof_co_max, &
            mom(1)/vv, mom(2)/vv, mom(3)/vv, mom(4), mom(5), vel_glb, vp_iters_last, vp_res_last, &
            vof_bnd_cum, vof_relax_cum, xloc, yloc
       Flush(vof_unit)
    End If
    vof_co_max = 0d0

  End Subroutine vof_output_monitor


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

    Real(Int64),    Intent(Out) :: vliq, cmin, cmax, mom(5)
    Integer(Int64), Intent(Out) :: nint

    Integer(Int32) :: i, j, k, ihi, jhi, khi
    Logical :: is_first, is_last
    Integer(Int32) :: partner
    Real(Int64) :: c, vc, buf(6), gbuf(6), mn, mx, kd
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
             mn = Min(mn, c);  mx = Max(mx, c)
             If ( c > vof_eps .And. c < 1d0 - vof_eps ) ni = ni + 1
          End Do
       End Do
    End Do

    Call MPI_Allreduce(buf, gbuf, 6, MPI_real8, MPI_SUM, MPI_COMM_WORLD, ierr)
    Call MPI_Allreduce(mn, cmin, 1, MPI_real8, MPI_MIN, MPI_COMM_WORLD, ierr)
    Call MPI_Allreduce(mx, cmax, 1, MPI_real8, MPI_MAX, MPI_COMM_WORLD, ierr)
    Call MPI_Allreduce(ni, gni, 1, MPI_integer8, MPI_SUM, MPI_COMM_WORLD, ierr)
    vliq = gbuf(1);  mom = gbuf(2:6);  nint = gni

  End Subroutine vof_diagnostics

End Module vof_state
