!> Two-fluid time step on the staggered grid (VOF liquid fraction C, one-fluid momentum):
!>
!>   Strang step:  A(dt/2) - project - RK3(forces, frozen density) - A(dt/2) - project
!>
!>   A(tau)  advances C with the Weymouth-Yue sweeps and the momentum rho*u with the same mass fluxes, so mass and momentum are
!>           transported consistently (uniform velocity is reproduced exactly at any density ratio); the c-tilde*div(u) term of
!>           the C update enters the momentum update with rho-tilde. Central face interpolation of the transported velocity.
!>   RK3     Wray stages with the solver's coefficients: viscous stress (variable viscosity), gravity with the still-water
!>           hydrostatic part removed, and the pressure as an explicit predictor (previous stage) plus an increment from a few
!>           PCG iterations preconditioned by the fast solver (vof_pressure).
!>
!> Supported: DNS/LES without wall model, periodic or wall boundaries, ghost-cell IBM (velocity condition), wave inlet and
!> relaxation zones. Not supported (checked at input): scalars, particles, UAV, wall models.
Module vof_twofluid

  Use iso_fortran_env, Only : Int32, Int64
  Use global
  Use mpi
  Use decomp, Only : x_periodic_partner, z_periodic_partner
  Use boundary_conditions, Only : apply_boundary_conditions, outflow_relax_on, update_ghost_interior_planes, &
                                  update_ghost_interior_planes_x, apply_periodic_bc_x, apply_periodic_bc_z
  Use sgs_models, Only : compute_sgs_model
  Use wallmodel, Only : compute_wall_model, solve_u_tau_wall
  Use ibm_stress
  Use waves, Only : wave_eta, wave_profile_at, current_profile
  Use monitor, Only : compute_cfl, select_dt, write_force_csv
  Use vof_plic
  Use vof_normals
  Use vof_advect
  Use vof_state
  Use vof_pressure
  Use profiler
  Use vof_curv, Only : vof_curvature
  Use mom_recon
  Use halo_pad, Only : pad_field
  Use ibm, Only : apply_ghost_cell_ibm, apply_ghost_cell_ibm_scalar, apply_ghost_cell_ibm_scalar_noflux
  Use scalar_transport, Only : compute_rhs_scalar_core, apply_scalar_bc, kf_on, kfield
  Use thermal_transport, Only : apply_temperature_bc
  Use particles, Only : advance_particles, particles_load_surface, particles_check_gas_seeds
  Use vof_ibm, Only : vof_compute_ibm_forces
  Use projection, Only : compute_pseudo_pressure_rhs, solve_poisson_equation, project_velocity

  Implicit None

  ! transported momentum density and its face velocity, mass flux of the current sweep, c-tilde density, volume-flux divergence
  Real(Int64), Allocatable, Dimension(:,:,:) :: qu, qv, qw, uqu, uqv, uqw, Md, rt, Dd
  ! pseudo-time RK3 of the momentum update: start-of-sweep momentum and the stage tendency
  Real(Int64), Allocatable, Dimension(:,:,:) :: q0u, q0v, q0w, dqu, dqv, dqw, Fl
  ! padded copies (E=3 planes each side) of the transported velocity for the wide momentum-reconstruction stencils
  Integer(Int32), Parameter :: EP = 3
  Real(Int64), Allocatable, Dimension(:,:,:) :: PadU, PadV, PadW
  ! staggered densities before the sweep and the smaller of before/after (mass of a control volume that empties or fills in the sweep)
  Real(Int64), Allocatable, Dimension(:,:,:) :: rou, rov, row, rmu, rmv, rmw
  ! transporting velocity (exactly divergence-free) and a scratch copy of the velocity
  Real(Int64), Allocatable, Dimension(:,:,:) :: Ut, Vt, Wt, Utmp, Vtmp, Wtmp
  ! RK3: start-of-step velocity is Uo/Vo/Wo (global); stage force combinations H = F_explicit - beta grad p of stages 1 and 2
  Real(Int64), Allocatable, Dimension(:,:,:) :: Hu1, Hv1, Hw1, Hu2, Hv2, Hw2, Fu_, Fv_, Fw_, Gu_, Gv_, Gw_
  Real(Int64), Allocatable, Dimension(:,:,:) :: ppre, dphi, fdiv, mu_c, mr_c, kap_c
  ! stage force workspace (predictor pressure gradient), padded-transport staging, transport-projection potential
  Real(Int64), Allocatable, Dimension(:,:,:) :: Bu_, Bv_, Bw_, tmpU, tmpV, tmpW, phT
  Real(Int64), Allocatable, Dimension(:)     :: rhos, rsf, cw0, cw1
  ! flat-wall EQWM stress fluxes on the u and w faces of the bottom (3rd index 1) and top (2) wall, and where they replace the molecular flux
  Real(Int64), Allocatable, Dimension(:,:,:) :: wt_u, wt_w
  ! scalar (sediment, temperature) work: start-of-step copy, stage RHS, effective y velocity with the settling folded in
  Real(Int64), Allocatable, Dimension(:,:,:) :: Fsc, Vsc, Csc0
  Logical,     Allocatable, Dimension(:,:,:) :: wt_on_u, wt_on_w
  ! IBM wall-model stress (ibm_stress): on when the IBM wall model is active; kinematic viscosity of the phase at the cells (0 = off)
  Logical :: ibm_wm_on = .False.
  ! relaxation-zone targets per x station (vof_relax)
  Real(Int64), Allocatable, Dimension(:)   :: rl_fr_c, rl_top, rl_fr_u
  Real(Int64), Allocatable, Dimension(:,:) :: rl_vf, rl_uc
  Real(Int64), Allocatable, Dimension(:,:,:) :: nu_cw
  Integer(Int32) :: vf_proj_its_last = 0, vf_proj_its_sum = 0
  ! IBM load output: whether this step is a sampling step
  Logical :: ibm_sampling_now = .False.
  Integer(Int32) :: vf_tr_warn = 0   ! transporting-velocity projections that did not converge (warnings are limited to the first few)
  Real(Int64) :: vf_bnd_loc = 0d0, vf_relax_loc = 0d0   ! rank-local liquid volume through the x boundaries / added by relaxation this step

Contains

  Subroutine vof_flow_init

    Integer(Int32) :: i, j, k
    Real(Int64) :: kt(3), om

    If ( uav_active >= 1 .Or. rotation_active >= 1 ) Then
       If ( myid == 0 ) Write(*,'(A)') ' ERROR: vof_flow=1 does not yet support ' // &
            'UAV or rotation'
       If ( myid == 0 ) Write(*,'(A,7I3)') ' ibm_wall_model, sediment, boussinesq, particles, uav, flat_wall, rotation: ', &
            ibm_wall_model_flag, sediment_flag, boussinesq_flag, particles_active, uav_active, flat_wall_model_flag, rotation_active
       Call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    End If
    If ( dPdx /= 0d0 .Or. dPdz /= 0d0 .Or. flow_forcing_mode /= 0 ) Then
       If ( myid == 0 ) Write(*,'(A)') ' ERROR: vof_flow=1 does not apply dPdx, dPdz or mass-flux forcing (set them to zero)'
       Call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    End If

    Allocate( wt_u(nx,nzg,2), wt_w(nxg,nz,2), wt_on_u(nx,nzg,2), wt_on_w(nxg,nz,2) )
    wt_u = 0d0;  wt_w = 0d0;  wt_on_u = .False.;  wt_on_w = .False.
    ibm_wm_on = ( ibm_input_mode >= 1 .And. ibm_wall_model_flag == 1 )
    If ( ibm_wm_on ) Then
       Allocate( nu_cw(nxg,nyg,nzg) )
       Call ibm_stress_init(.True.)
    End If
    Allocate( qu(nx,nyg,nzg), qv(nxg,ny,nzg), qw(nxg,nyg,nz), uqu(nx,nyg,nzg), uqv(nxg,ny,nzg), uqw(nxg,nyg,nz) )
    Allocate( Md(nxg,nyg,nzg), rt(nxg,nyg,nzg), Dd(nxg,nyg,nzg) )
    Allocate( Ut(nx,nyg,nzg), Vt(nxg,ny,nzg), Wt(nxg,nyg,nz), Utmp(nx,nyg,nzg), Vtmp(nxg,ny,nzg), Wtmp(nxg,nyg,nz) )
    Allocate( Hu1(nx,nyg,nzg), Hv1(nxg,ny,nzg), Hw1(nxg,nyg,nz), Hu2(nx,nyg,nzg), Hv2(nxg,ny,nzg), Hw2(nxg,nyg,nz) )
    Allocate( Fu_(nx,nyg,nzg), Fv_(nxg,ny,nzg), Fw_(nxg,nyg,nz), Gu_(nx,nyg,nzg), Gv_(nxg,ny,nzg), Gw_(nxg,nyg,nz) )
    Allocate( ppre(nxg,nyg,nzg), dphi(nxg,nyg,nzg), fdiv(nxg,nyg,nzg), mu_c(nxg,nyg,nzg), mr_c(nxg,nyg,nzg), kap_c(nxg,nyg,nzg) )
    kap_c = 0d0
    Allocate( Bu_(nx,nyg,nzg), Bv_(nxg,ny,nzg), Bw_(nxg,nyg,nz), phT(nxg,nyg,nzg) )
    Allocate( tmpU(1-EP:nx+EP,nyg,1-EP:nzg+EP), tmpV(1-EP:nxg+EP,ny,1-EP:nzg+EP), tmpW(1-EP:nxg+EP,nyg,1-EP:nz+EP) )
    Allocate( rhos(nyg), rsf(ny), cw0(nyg), cw1(nyg) )
    Allocate( PadU(1-EP:nx+EP,1-EP:nyg+EP,1-EP:nzg+EP), PadV(1-EP:nxg+EP,1-EP:ny+EP,1-EP:nzg+EP), &
              PadW(1-EP:nxg+EP,1-EP:nyg+EP,1-EP:nz+EP) )
    PadU = 0d0;  PadV = 0d0;  PadW = 0d0
    Allocate( rou(nx,nyg,nzg), rov(nxg,ny,nzg), row(nxg,nyg,nz), rmu(nx,nyg,nzg), rmv(nxg,ny,nzg), rmw(nxg,nyg,nz) )
    rou = 1d0;  rov = 1d0;  row = 1d0;  rmu = 1d0;  rmv = 1d0;  rmw = 1d0
    Allocate( q0u(nx,nyg,nzg), q0v(nxg,ny,nzg), q0w(nxg,nyg,nz), dqu(nx,nyg,nzg), dqv(nxg,ny,nzg), dqw(nxg,nyg,nz) )
    qu = 0d0;  qv = 0d0;  qw = 0d0;  Md = 0d0;  rt = 0d0;  Dd = 0d0
    Allocate( Fl(0:nxg+1,0:nyg+1,0:nzg+1) )
    q0u = 0d0;  q0v = 0d0;  q0w = 0d0;  dqu = 0d0;  dqv = 0d0;  dqw = 0d0;  Fl = 0d0
    ppre = 0d0;  dphi = 0d0

    ! cell-centre weights of the y-face quantities (v) of the stretched grid
    cw0 = 0.5d0
    Do j = 2, nyg-1
       cw0(j) = ( y(j) - yg(j) )/( y(j) - y(j-1) )
    End Do
    cw1 = 1d0 - cw0

    ! still-water reference density per row of the initial state (hydrostatic part removed from gravity)
    rhos = vof_rho_ref
    Do j = 1, ny
       rsf(j) = ( rhos(j)*vp_hy(j) + rhos(j+1)*vp_hy(j+1) )/( vp_hy(j) + vp_hy(j+1) )
    End Do
    If ( restart == 0 .And. vof_tgv == 1 ) Then
       kt(1) = 8d0*Atan(1d0)/(Lx - dx);  kt(2) = 8d0*Atan(1d0)/Ly;  kt(3) = 8d0*Atan(1d0)/(Lz - hz(2))
       Do k = 1, nzg
          Do j = 1, nyg
             Do i = 1, nx
                U(i,j,k) = vof_u0*Sin(kt(1)*x(i))*Cos(kt(2)*yg(j))*Cos(kt(3)*zg(k))
             End Do
          End Do
       End Do
       Do k = 1, nzg
          Do j = 1, ny
             Do i = 1, nxg
                V(i,j,k) = -vof_u0*Cos(kt(1)*xg(i))*Sin(kt(2)*y(j))*Cos(kt(3)*zg(k))
             End Do
          End Do
       End Do
       W = 0d0
    Else If ( restart == 0 .And. vof_wave_stokes == 1 ) Then
       kt(1) = 8d0*Atan(1d0)/vof_wave_lambda
       om = Sqrt(vof_grav*kt(1))
       V = 0d0;  W = 0d0
       Do k = 1, nzg
          Do j = 1, nyg
             Do i = 1, nx
                U(i,j,k) = Merge(om*vof_wave_amp*Exp(kt(1)*(yg(j) - vof_level))*Cos(kt(1)*x(i)), 0d0, &
                                 yg(j) < vof_level + wave_surface(x(i)))
             End Do
          End Do
       End Do
       Do k = 1, nzg
          Do j = 1, ny
             Do i = 1, nxg
                V(i,j,k) = Merge(om*vof_wave_amp*Exp(kt(1)*(y(j) - vof_level))*Sin(kt(1)*xg(i)), 0d0, &
                                 y(j) < vof_level + wave_surface(xg(i)))
             End Do
          End Do
       End Do
    Else If ( restart == 0 .And. vof_u0 /= 0d0 ) Then
       U = vof_u0;  V = 0d0;  W = 0d0
    Else If ( restart == 0 .And. current_type > 0 ) Then
       Block
          Real(Int64) :: ucur(nyg)
          Call current_profile(t, ucur)
          V = 0d0;  W = 0d0
          Do j = 1, nyg
             U(:,j,:) = ucur(j)
          End Do
       End Block
    End If

    ! the initial fields above were set on the host
    !$acc update device(U,V,W)

    ! pseudo-time RK3 of the momentum needs less than half of a control volume's mass to leave in one sweep, which fails for
    ! density ratios above a few thousand (drop in a uniform stream, 1e4 .. 1e6): forward Euler with small sub-steps there
    If ( vof_rk_nth <= 0 ) vof_rk_nth = 1
    If ( vof_rho_l/vof_rho_g >= 2d3 .And. vof_rk_mom >= 1 ) Then
       vof_rk_mom = 0;  vof_co_sub = Min(vof_co_sub, 3d-2)
    End If

    ! the viscous time-step limit sees the largest kinematic viscosity of the two fluids
    nu = Max(vof_nu_l, vof_nu_g)
    !$acc update device(nu)

    !$acc enter data copyin(qu,qv,qw,uqu,uqv,uqw,Md,rt,Dd,Ut,Vt,Wt,Utmp,Vtmp,Wtmp,Hu1,Hv1,Hw1,Hu2,Hv2,Hw2,Fu_,Fv_,Fw_,Gu_,Gv_,Gw_) &
    !$acc& copyin(ppre,dphi,fdiv,mu_c,mr_c,kap_c,Bu_,Bv_,Bw_,phT,tmpU,tmpV,tmpW,rsf,cw0,cw1,wt_u,wt_w,wt_on_u,wt_on_w) &
    !$acc& copyin(PadU,PadV,PadW,rou,rov,row,rmu,rmv,rmw,q0u,q0v,q0w,dqu,dqv,dqw,Fl)
    If ( ibm_wm_on ) Then
       !$acc enter data copyin(nu_cw)
    End If
    If ( restart == 1 ) Then
       ppre = P
       !$acc update device(ppre)
       Call vp_halo(ppre, .True.)
    Else
       Call vof_initial_pressure
    End If
    If ( particles_active >= 1 ) Then
       !$acc update host(Cv)
       Call particles_load_surface(Cv(1:nxg,1:nyg,1:nzg))
       Call particles_check_gas_seeds
    End If

  End Subroutine vof_flow_init


  !> Exact dynamic pressure of the initial state (u = 0): div(beta grad p) = div(gravity deficit), the first predictor
  Subroutine vof_initial_pressure

    Real(Int64), Allocatable :: gz(:,:,:), zu(:,:,:), zw(:,:,:)
    Integer(Int32) :: i, j, k

    Call vof_fill_pad(Cv, nxg, nyg, nzg)
    Call vp_set_density(Cv)
    Allocate( gz(nxg,ny,nzg), zu(nx,nyg,nzg), zw(nxg,nyg,nz) )
    zu = 0d0;  zw = 0d0
    !$acc data copyin(zu,zw) create(gz)
    !$acc parallel loop collapse(3) present(gz,vp_rfv,vp_riv,rsf)
    Do k = 1, nzg
       Do j = 1, ny
          Do i = 1, nxg
             gz(i,j,k) = -vof_grav*( vp_rfv(i,j,k) - rsf(j) )/vp_riv(i,j,k)
          End Do
       End Do
    End Do
    If ( vof_hsplit /= 0 ) Then
       !$acc kernels present(gz)
       gz = 0d0
       !$acc end kernels
    End If
    If ( vof_sigma > 0d0 ) Then
       Call vof_curvature(kap_c)
       Call add_surface_tension(zu, gz, zw)
    End If
    Call vp_div(zu, gz, zw, fdiv)
    Call vp_pcg(fdiv, 400, 1d-12, ppre)
    Call vp_halo(ppre, .True.)
    !$acc end data
    If ( myid == 0 ) Write(*,'(A,I4,A,ES10.2)') '   VOF initial pressure: PCG iterations', vp_iters_last, '  residual', vp_res_last
    Deallocate( gz, zu, zw )

  End Subroutine vof_initial_pressure


  !> Halo exchange of face fields shaped like U, V, W (rank seams and periodic wraps); wall rows are left to the caller
  Subroutine face_halo(Fu, Fv, Fw)

    Real(Int64), Intent(InOut) :: Fu(nx,nyg,nzg), Fv(nxg,ny,nzg), Fw(nxg,nyg,nz)

    Call update_ghost_interior_planes(Fu, 1)
    Call update_ghost_interior_planes(Fv, 2)
    Call update_ghost_interior_planes(Fw, 3)
    Call update_ghost_interior_planes_x(Fu, 1)
    Call update_ghost_interior_planes_x(Fv, 2)
    Call update_ghost_interior_planes_x(Fw, 2)
    If ( x_bc_type == 0 ) Then
       Call apply_periodic_bc_x(Fu, 1)
       Call apply_periodic_bc_x(Fv, 2)
       Call apply_periodic_bc_x(Fw, 2)
    End If
    If ( z_bc_type == 0 ) Then
       Call apply_periodic_bc_z(Fu, 1)
       Call apply_periodic_bc_z(Fv, 2)
       Call apply_periodic_bc_z(Fw, 3)
    End If

  End Subroutine face_halo


  !> Ghost-cell IBM on U, V, W followed by a halo refresh: the image-point corrections change interior values next to the rank
  !  seams, whose neighbours' ghost planes must follow before the next stencil reads them
  Subroutine enforce_ibm

    Call apply_ghost_cell_ibm(U, V, W)
    Call face_halo(U, V, W)

  End Subroutine enforce_ibm


  !> Advance C and the momentum rho*u by tau with the current velocity U, V, W as the (divergence-free) transporting field
  Subroutine vof_advect_half(tau, keep_transport)

    Real(Int64), Intent(In) :: tau
    Logical, Intent(In), Optional :: keep_transport

    Integer(Int32) :: isw, d, order(3), i, j, k
    Real(Int64) :: co, cl

    If ( .Not. Present(keep_transport) ) Then
       Call profiler_start(PROF_VOF_TRANSPORT)
       Call make_transport_velocity
       Call profiler_stop(PROF_VOF_TRANSPORT)
    Else If ( .Not. keep_transport ) Then
       Call profiler_start(PROF_VOF_TRANSPORT)
       Call make_transport_velocity
       Call profiler_stop(PROF_VOF_TRANSPORT)
    End If
    Call profiler_start(PROF_VOF_ADVECT)
    Call profiler_start(PROF_ADV_DENSITY)
    Call vof_fill_pad(Cv, nxg, nyg, nzg)
    Call vp_set_density(Cv, sweep_only=.True.)
    !$acc kernels present(qu,qv,qw,vp_rau,vp_rav,vp_raw,U,V,W,vof_cc,Cv,rt)
    qu = vp_rau*U;  qv = vp_rav*V;  qw = vp_raw*W
    vof_cc = Merge(1d0, 0d0, Cv(1:nxg,1:nyg,1:nzg) > 0.5d0)
    rt = vof_rho_g + (vof_rho_l - vof_rho_g)*vof_cc
    !$acc end kernels
    Call profiler_stop(PROF_ADV_DENSITY)
    co = 0d0;  cl = 0d0

    If ( Mod(vof_nadv, 2) == 0 ) Then
       order = (/ 1, 2, 3 /)
    Else
       order = (/ 3, 2, 1 /)
    End If

    Do isw = 1, 3
       d = order(isw)
       ! Cv, its pad and the density are already current here (set above for the first sweep, by momentum_update after the
       ! previous one); the geometric density has also reconstructed the planes
       Call profiler_start(PROF_ADV_RECON)
       If ( vof_geo_density < 1 ) Call vof_reconstruct(Cv, nxg, nyg, nzg, vof_normal_scheme)
       Call profiler_stop(PROF_ADV_RECON)
       !$acc kernels present(uqu,uqv,uqw,qu,qv,qw,vp_rau,vp_rav,vp_raw,rou,rov,row)
       uqu = qu/vp_rau;  uqv = qv/vp_rav;  uqw = qw/vp_raw
       !$acc end kernels
       Call dev_copy(rou, vp_rau)
       Call dev_copy(rov, vp_rav)
       Call dev_copy(row, vp_raw)
       Call profiler_start(PROF_ADV_SWEEP)
       If ( d == 1 ) Then
          Call vof_sweep(1, Cv, nxg, nyg, nzg, Ut, hx, hy, hz, tau, co, cl)
       Else If ( d == 2 ) Then
          Call vof_sweep(2, Cv, nxg, nyg, nzg, Vt, hx, hy, hz, tau, co, cl)
       Else
          Call vof_sweep(3, Cv, nxg, nyg, nzg, Wt, hx, hy, hz, tau, co, cl)
       End If
       Call profiler_stop(PROF_ADV_SWEEP)
       Call profiler_start(PROF_ADV_MASSFLUX)
       If ( d == 1 .And. x_bc_type == 1 ) Call tally_x_boundary_flux
       Call sweep_mass_flux(d, tau)
       Call vp_halo(Md, .False.)
       Call vp_halo(Dd, .False.)
       !$acc kernels present(Md)
       Md(:,1,:) = 0d0;  Md(:,nyg,:) = 0d0
       !$acc end kernels
       Call profiler_stop(PROF_ADV_MASSFLUX)
       Call profiler_start(PROF_ADV_MOMENTUM)
       Call momentum_update(d)
       Call face_halo(qu, qv, qw)
       Call profiler_stop(PROF_ADV_MOMENTUM)
    End Do
    vof_nadv = vof_nadv + 1
    vof_co_max = Max(vof_co_max, co)
    vof_clip_total = vof_clip_total + cl

    !$acc parallel loop collapse(3) present(U,qu,vp_rau)
    Do k = 2, nzg-1
       Do j = 2, nyg-1
          Do i = 2, nx-1
             U(i,j,k) = qu(i,j,k)/vp_rau(i,j,k)
          End Do
       End Do
    End Do
    !$acc parallel loop collapse(3) present(V,qv,vp_rav)
    Do k = 2, nzg-1
       Do j = 2, ny-1
          Do i = 2, nxg-1
             V(i,j,k) = qv(i,j,k)/vp_rav(i,j,k)
          End Do
       End Do
    End Do
    !$acc parallel loop collapse(3) present(W,qw,vp_raw)
    Do k = 2, nz-1
       Do j = 2, nyg-1
          Do i = 2, nxg-1
             W(i,j,k) = qw(i,j,k)/vp_raw(i,j,k)
          End Do
       End Do
    End Do

    ! the interior was rewritten: refresh the seam and periodic ghost faces before the next sub-step rebuilds q = rho*u from them
    Call face_halo(U, V, W)
    Call profiler_stop(PROF_VOF_ADVECT)

  End Subroutine vof_advect_half


  !> Transporting velocity: the current velocity made exactly divergence-free by one constant-coefficient fast solve, so the
  !  volume of liquid is conserved to round-off whatever residual the variable-density PCG leaves
  Subroutine make_transport_velocity

    If ( vp_masked ) Then
       Call make_transport_velocity_masked
       Return
    End If
    Call dev_copy(Utmp, U)
    Call dev_copy(Vtmp, V)
    Call dev_copy(Wtmp, W)
    Call compute_pseudo_pressure_rhs
    Call solve_poisson_equation(skip_p_save=.True.)
    Call project_velocity
    Call dev_copy(Ut, U)
    Call dev_copy(Vt, V)
    Call dev_copy(Wt, W)
    Call face_halo(Ut, Vt, Wt)
    Call dev_copy(U, Utmp)
    Call dev_copy(V, Vtmp)
    Call dev_copy(W, Wtmp)

  End Subroutine make_transport_velocity


  !> Transporting velocity with an immersed body: the fast solver cannot honour the closed faces, so the constant-coefficient
  !  projection is done by PCG on the masked operator (Neumann at the body), to round-off. Closed faces carry no flux, hence no
  !  liquid or momentum crosses the body, and the divergence-free field is exact in every fluid cell.
  Subroutine make_transport_velocity_masked

    Real(Int64) :: beta0_save, umax

    beta0_save = vp_beta0
    Call dev_copy(vp_bu, vp_mu)
    Call dev_copy(vp_bv, vp_mv)
    Call dev_copy(vp_bw, vp_mw)
    Call dev_copy(Utmp, U)
    Call dev_copy(Vtmp, V)
    Call dev_copy(Wtmp, W)
    !$acc kernels present(vp_mu,vp_mv,vp_mw,U,V,W)
    U = U*vp_mu;  V = V*vp_mv;  W = W*vp_mw
    !$acc end kernels
    vp_beta0 = 1d0
    vp_mg_stale = ( pcg_precond >= 1 )

    Call vp_div(U, V, W, fdiv)
    ! the iteration stops at the round-off level of the velocity (a relative tolerance alone is unreachable once the field is
    ! already nearly divergence free, and PCG then amplifies noise)
    umax = Max( maxabs3(U, nx, nyg, nzg), maxabs3(V, nxg, ny, nzg), maxabs3(W, nxg, nyg, nz) )
    Call MPI_Allreduce(MPI_IN_PLACE, umax, 1, MPI_real8, MPI_MAX, MPI_COMM_WORLD, ierr)
    Call profiler_start(PROF_PCG)
    Call profiler_start(PROF_PCG_TRANSPORT)
    Call vp_pcg(fdiv, 800, 1d-13, phT, rfloor=vof_tr_floor*umax*Sqrt(vp_wsum)/Min(dx, dymin, dzmin))
    Call profiler_stop(PROF_PCG_TRANSPORT)
    Call profiler_stop(PROF_PCG)
    ! Ut must be solenoidal for the liquid volume to be conserved: report a solve that hit the cap or broke down at the start
    If ( ( vp_iters_last >= 800 .And. vp_res_last > 1d-6 ) .Or. &
         ( vp_iters_last == 0 .And. maxabs3(fdiv, nxg, nyg, nzg) > 0d0 ) ) Then
       vf_tr_warn = vf_tr_warn + 1
       If ( myid == 0 .And. vf_tr_warn <= 5 ) Then
          Write(*,'(A,I0,A,I0,A,ES9.2)') ' WARNING: transporting-velocity projection ended after ', vp_iters_last, &
               ' iterations (cap 800) at step ', istep, ', relative residual ', vp_res_last
       End If
    End If
    Call vp_halo(phT, .True.)
    Call vp_grad(phT, Gu_, Gv_, Gw_)
    Call apply_face_gradient(Gu_, Gv_, Gw_, 1d0)
    Call dev_copy(Ut, U)
    Call dev_copy(Vt, V)
    Call dev_copy(Wt, W)
    Call face_halo(Ut, Vt, Wt)
    Call dev_copy(U, Utmp)
    Call dev_copy(V, Vtmp)
    Call dev_copy(W, Wtmp)

    ! vp_bu.. keep the mask coefficients: every solve with the density coefficients starts with a full vp_set_density
    vp_beta0 = beta0_save
    vp_mg_stale = ( pcg_precond >= 1 )

  End Subroutine make_transport_velocity_masked


  !> max|a| over a whole face-shaped array (ghost planes included), this rank only
  Function maxabs3(a, n1, n2, n3) Result(s)

    Integer(Int32), Intent(In) :: n1, n2, n3
    Real(Int64), Intent(In) :: a(n1,n2,n3)
    Real(Int64) :: s
    Integer(Int32) :: i, j, k

    s = 0d0
    !$acc parallel loop collapse(3) reduction(max:s) present(a)
    Do k = 1, n3
       Do j = 1, n2
          Do i = 1, n1
             s = Max(s, Abs(a(i,j,k)))
          End Do
       End Do
    End Do

  End Function maxabs3


  !> Liquid volume that crossed the inlet (positive in) and outlet (positive out) faces in the x-sweep just done
  Subroutine tally_x_boundary_flux

    Logical :: is_first_x, is_last_x, is_first_z, is_last_z
    Integer(Int32) :: partner, khi, j, k
    Real(Int64) :: s

    khi = nzg-1
    If ( z_bc_type == 0 ) Then
       Call z_periodic_partner(is_first_z, is_last_z, partner)
       If ( is_last_z ) khi = nzg-2
    End If
    Call x_periodic_partner(is_first_x, is_last_x, partner)
    If ( is_first_x ) Then
       s = 0d0
       !$acc parallel loop collapse(2) reduction(+:s) present(vof_flux)
       Do k = 2, khi
          Do j = 2, nyg-1
             s = s + vof_flux(1,j,k)
          End Do
       End Do
       vf_bnd_loc = vf_bnd_loc + s
    End If
    If ( is_last_x ) Then
       s = 0d0
       !$acc parallel loop collapse(2) reduction(+:s) present(vof_flux)
       Do k = 2, khi
          Do j = 2, nyg-1
             s = s + vof_flux(nx,j,k)
          End Do
       End Do
       vf_bnd_loc = vf_bnd_loc - s
    End If

  End Subroutine tally_x_boundary_flux


  !> Mass crossing each face of the sweep direction d in this sweep, and the volume-flux divergence Dd of the sweep
  Subroutine sweep_mass_flux(d, tau)

    Integer(Int32), Intent(In) :: d
    Real(Int64),    Intent(In) :: tau
    Integer(Int32) :: i, j, k
    Real(Int64) :: drho

    drho = vof_rho_l - vof_rho_g
    !$acc kernels present(Md,Dd)
    Md = 0d0
    Dd = 0d0
    !$acc end kernels
    Select Case(d)
    Case(1)
       !$acc parallel loop collapse(3) present(Md,Ut,hy,hz,vof_flux)
       Do k = 2, nzg-1
          Do j = 2, nyg-1
             Do i = 1, nx
                Md(i,j,k) = vof_rho_g*Ut(i,j,k)*tau*hy(j)*hz(k) + drho*vof_flux(i,j,k)
             End Do
          End Do
       End Do
       !$acc parallel loop collapse(3) present(Dd,Ut,hy,hz)
       Do k = 2, nzg-1
          Do j = 2, nyg-1
             Do i = 2, nx
                Dd(i,j,k) = tau*hy(j)*hz(k)*( Ut(i,j,k) - Ut(i-1,j,k) )
             End Do
          End Do
       End Do
    Case(2)
       !$acc parallel loop collapse(3) present(Md,Vt,hz,vof_flux)
       Do k = 2, nzg-1
          Do j = 1, ny
             Do i = 2, nxg-1
                Md(i,j,k) = vof_rho_g*Vt(i,j,k)*tau*dx*hz(k) + drho*vof_flux(i,j,k)
             End Do
          End Do
       End Do
       !$acc parallel loop collapse(3) present(Dd,Vt,hz)
       Do k = 2, nzg-1
          Do j = 2, ny
             Do i = 2, nxg-1
                Dd(i,j,k) = tau*dx*hz(k)*( Vt(i,j,k) - Vt(i,j-1,k) )
             End Do
          End Do
       End Do
    Case Default
       !$acc parallel loop collapse(3) present(Md,Wt,hy,vof_flux)
       Do k = 1, Min(nz, nzg-1)
          Do j = 2, nyg-1
             Do i = 2, nxg-1
                Md(i,j,k) = vof_rho_g*Wt(i,j,k)*tau*dx*hy(j) + drho*vof_flux(i,j,k)
             End Do
          End Do
       End Do
       !$acc parallel loop collapse(3) present(Dd,Wt,hy)
       Do k = 2, nz
          Do j = 2, nyg-1
             Do i = 2, nxg-1
                Dd(i,j,k) = tau*dx*hy(j)*( Wt(i,j,k) - Wt(i,j,k-1) )
             End Do
          End Do
       End Do
    End Select

  End Subroutine sweep_mass_flux


  !> Momentum update of the sweep along d with the exact mass fluxes Md of the geometric sweep. The staggered density moves
  !  linearly from rou (before) to vp_rau (after) along the sweep pseudo-time theta, so the momentum ODE
  !  dq/dtheta = tendency(q/rho(theta)) is integrated by forward Euler (vof_rk_mom = 0) SSP-RK3 (1; theta = 0, 1, 1/2) or SSP-RK2 (2; theta = 0, 1)
  !  with the same fluxes at every stage; vof_rk_nth > 1 splits the pseudo-time into equal segments (a stage may step
  !  past the end of the density path, which needs less than half of a control volume's mass to leave in one segment)
  Subroutine momentum_update(d)

    Integer(Int32), Intent(In) :: d
    Integer(Int32) :: stage, nst, seg
    Real(Int64) :: a0, a1, th

    Call profiler_start(PROF_MOM_DENSITY)
    Call vof_fill_pad(Cv, nxg, nyg, nzg)
    Call vp_set_density(Cv, sweep_only=.True.)
    Call profiler_stop(PROF_MOM_DENSITY)
    !$acc kernels present(rmu,rmv,rmw,rou,rov,row,vp_rau,vp_rav,vp_raw)
    rmu = Min(rou, vp_rau);  rmv = Min(rov, vp_rav);  rmw = Min(row, vp_raw)
    !$acc end kernels
    nst = Merge(3, Merge(2, 1, vof_rk_mom == 2), vof_rk_mom == 1)
    Do seg = 1, vof_rk_nth
       Call profiler_start(PROF_MOM_QPREP)
       Call dev_copy(q0u, qu)
       Call dev_copy(q0v, qv)
       Call dev_copy(q0w, qw)
       Call profiler_stop(PROF_MOM_QPREP)
       Do stage = 1, nst
          If ( stage > 1 .Or. seg > 1 ) Then
             Call profiler_start(PROF_MOM_HALO)
             Call face_halo(qu, qv, qw)
             Call profiler_stop(PROF_MOM_HALO)
             Call profiler_start(PROF_MOM_QPREP)
             th = ( Real(seg - 1, Int64) + Merge(0d0, Merge(1d0, 0.5d0, stage == 2), stage == 1) )/vof_rk_nth
             !$acc kernels present(uqu,uqv,uqw,qu,qv,qw,rou,rov,row,vp_rau,vp_rav,vp_raw)
             uqu = qu/( rou + th*(vp_rau - rou) );  uqv = qv/( rov + th*(vp_rav - rov) );  uqw = qw/( row + th*(vp_raw - row) )
             !$acc end kernels
             Call profiler_stop(PROF_MOM_QPREP)
          End If
          Call profiler_start(PROF_MOM_PAD)
          Call pad_transported
          Call profiler_stop(PROF_MOM_PAD)
          Call profiler_start(PROF_MOM_TENDENCY)
          Call momentum_tendency(d)
          Call profiler_stop(PROF_MOM_TENDENCY)
          Call profiler_start(PROF_MOM_UPDATE)
          a0 = Merge(0d0, Merge(Merge(0.5d0, 0.75d0, vof_rk_mom == 2), 1d0/3d0, stage == 2), stage == 1);  a1 = 1d0 - a0
          !$acc kernels present(qu,qv,qw,q0u,q0v,q0w,dqu,dqv,dqw)
          qu = a0*q0u + a1*( qu + dqu/vof_rk_nth );  qv = a0*q0v + a1*( qv + dqv/vof_rk_nth )
          qw = a0*q0w + a1*( qw + dqw/vof_rk_nth )
          !$acc end kernels
          Call profiler_stop(PROF_MOM_UPDATE)
       End Do
    End Do

  End Subroutine momentum_update


  !> Flux of q through the faces of each staggered control volume (mass flux averaged onto the CV face, face value of the
  !  transported velocity uq) plus the rho-tilde * volume-flux-divergence correction, divided by the control volume
  Subroutine momentum_tendency(d)

    Integer(Int32), Intent(In) :: d

    Call tend_comp(1, d, dqu, uqu, 2, nx-1,  2, nyg-1, 2, nzg-1)
    Call tend_comp(2, d, dqv, uqv, 2, nxg-1, 2, ny-1,  2, nzg-1)
    Call tend_comp(3, d, dqw, uqw, 2, nxg-1, 2, nyg-1, 2, nz-1)

  End Subroutine momentum_tendency


  !> Tendency of component c (1 u, 2 v, 3 w) for the sweep along d: each control-volume face flux is evaluated once and used
  !  by the two control volumes that share it
  Subroutine tend_comp(c, d, dq, uc, i0, i1, j0, j1, k0, k1)

    Integer(Int32), Intent(In)    :: c, d, i0, i1, j0, j1, k0, k1
    Real(Int64),    Intent(InOut) :: dq(:,:,:)
    Real(Int64),    Intent(In)    :: uc(:,:,:)
    Integer(Int32) :: i, j, k, ec1, ec2, ec3, ed1, ed2, ed3
    Real(Int64) :: corr, vcv

    ec1 = Merge(1,0,c==1);  ec2 = Merge(1,0,c==2);  ec3 = Merge(1,0,c==3)
    ed1 = Merge(1,0,d==1);  ed2 = Merge(1,0,d==2);  ed3 = Merge(1,0,d==3)
    Select Case(c)
    Case(1)
       Call face_fluxes(c, d, i0, i1, j0, j1, k0, k1, PadU, rmu, weight_y_0, weight_y_1, 0)
    Case(2)
       Call face_fluxes(c, d, i0, i1, j0, j1, k0, k1, PadV, rmv, cw0, cw1, 1)
    Case Default
       Call face_fluxes(c, d, i0, i1, j0, j1, k0, k1, PadW, rmw, weight_y_0, weight_y_1, 0)
    End Select
    !$acc parallel loop collapse(3) private(vcv,corr) present(dq,uc,rt,Dd,Fl,hy,hz,yg,zg)
    Do k = k0, k1
       Do j = j0, j1
          Do i = i0, i1
             vcv = dx*Merge(yg(j+1) - yg(j), hy(j), c==2)*Merge(zg(k+1) - zg(k), hz(k), c==3)
             corr = uc(i,j,k)*0.5d0*( rt(i,j,k)*Dd(i,j,k) + rt(i+ec1,j+ec2,k+ec3)*Dd(i+ec1,j+ec2,k+ec3) )
             dq(i,j,k) = ( corr - ( Fl(i,j,k) - Fl(i-ed1,j-ed2,k-ed3) ) )/vcv
          End Do
       End Do
    End Do

  End Subroutine tend_comp


  !> Fl = momentum flux of component c through the face of its control volume on the high side in direction d (the flux through the
  !  low face of the next control volume): centre-plane mass flux of the two cells the CV spans (Md), face value of the padded
  !  transported velocity pad blended to first-order upwind where the control volume is refilled within a step (rm: the smaller
  !  of the CV densities before and after the sweep). vof_mom_scheme 0 is the weighted central value (stretched-grid weights
  !  wy0/wy1 at j+wsh along y), others come from mom_recon (index-space formulas); WENO5-Z uses the branch-free form of mom_recon.
  Subroutine face_fluxes(c, d, i0, i1, j0, j1, k0, k1, pad, rm, wy0, wy1, wsh)

    Integer(Int32), Intent(In) :: c, d, i0, i1, j0, j1, k0, k1, wsh
    Real(Int64),    Intent(In) :: pad(1-EP:,1-EP:,1-EP:), rm(:,:,:), wy0(:), wy1(:)
    Integer(Int32) :: i, j, k, m, ec1, ec2, ec3, ed1, ed2, ed3, ii, jj, kk
    Real(Int64) :: ivcv, vcv, mfc, wa, wb, r1, r2, cm, up, t, rv, s(-2:3), icm, cm0, cm1

    ec1 = Merge(1,0,c==1);  ec2 = Merge(1,0,c==2);  ec3 = Merge(1,0,c==3)
    ed1 = Merge(1,0,d==1);  ed2 = Merge(1,0,d==2);  ed3 = Merge(1,0,d==3)
    cm0 = vof_mom_cm0;  cm1 = vof_mom_cm1;  icm = 1d0/(cm1 - cm0)
    If ( vof_mom_scheme == MOM_WENO5 ) Then
       !$acc parallel loop collapse(3) private(ivcv) present(pad,rm,Md,Fl,hy,hz,yg,zg)
       Do k = k0-ed3, k1
          Do j = j0-ed2, j1
             Do i = i0-ed1, i1
                ivcv = 1d0/( dx*Merge(yg(j+1) - yg(j), hy(j), c==2)*Merge(zg(k+1) - zg(k), hz(k), c==3) )
                Fl(i,j,k) = mom_flux_weno5z_pt( pad(i-2*ed1, j-2*ed2, k-2*ed3), pad(i-ed1, j-ed2, k-ed3), pad(i,j,k), &
                     pad(i+ed1, j+ed2, k+ed3), pad(i+2*ed1, j+2*ed2, k+2*ed3), pad(i+3*ed1, j+3*ed2, k+3*ed3), &
                     Md(i,j,k), Md(i+ec1,j+ec2,k+ec3), rm(i,j,k), rm(i+ed1,j+ed2,k+ed3), ivcv, cm0, icm )
             End Do
          End Do
       End Do
    Else
       !$acc parallel loop collapse(3) private(s,mfc,wa,wb,r1,r2,cm,up,t,rv,vcv) &
       !$acc& present(pad,rm,Md,Fl,hy,hz,yg,zg,wy0,wy1)
       Do k = k0-ed3, k1
          Do j = j0-ed2, j1
             Do i = i0-ed1, i1
                ii = i + ed1;  jj = j + ed2;  kk = k + ed3
                mfc = 0.5d0*( Md(i,j,k) + Md(i+ec1,j+ec2,k+ec3) )
                Do m = -2, 3
                   s(m) = pad(i + m*ed1, j + m*ed2, k + m*ed3)
                End Do
                r1 = rm(i,j,k);  r2 = rm(ii,jj,kk)
                vcv = dx*Merge(yg(j+1) - yg(j), hy(j), c==2)*Merge(zg(k+1) - zg(k), hz(k), c==3)
                wa = 0.5d0;  wb = 0.5d0
                If ( d == 2 ) Then;  wa = wy0(j+wsh);  wb = wy1(j+wsh);  End If
                cm = Abs(mfc)*Abs(1d0/r1 - 1d0/r2)/vcv
                up = Merge(s(0), s(1), mfc >= 0d0)
                If ( vof_mom_scheme == 0 ) Then
                   rv = wa*s(0) + wb*s(1)
                Else
                   rv = mom_face(s, mfc, vof_mom_scheme)
                End If
                t = Min(1d0, Max(0d0, (cm - cm0)/(cm1 - cm0)))
                rv = up + (1d0 - t*t*(3d0 - 2d0*t))*(rv - up)
                Fl(i,j,k) = mfc*rv
             End Do
          End Do
       End Do
    End If

  End Subroutine face_fluxes


  !> Padded copies of uqu, uqv, uqw: EP planes beyond the ghost plane in x and z (neighbour rank / periodic partner via pad_field),
  !  EP rows beyond the wall ghost row in y by reflection (tangential: odd for no-slip, even for free-slip; normal: odd about the wall face)
  Subroutine pad_transported

    Real(Int64) :: s_lo, s_hi
    Integer(Int32) :: m

    s_lo = Merge(-1d0, 1d0, bc_face_ylo == 1);  s_hi = Merge(-1d0, 1d0, bc_face_yhi == 1)
    Call pad_field_dev(uqu, tmpU, nx, nyg, nzg, .True., .False., EP)
    Call pad_field_dev(uqv, tmpV, nxg, ny, nzg, .False., .False., EP)
    Call pad_field_dev(uqw, tmpW, nxg, nyg, nz, .False., .True., EP)
    !$acc kernels present(PadU,PadV,PadW,tmpU,tmpV,tmpW)
    PadU(:,1:nyg,:) = tmpU
    PadV(:,1:ny,:) = tmpV
    PadW(:,1:nyg,:) = tmpW
    ! the wall ghost rows: the stored ones are the boundary-condition values of the start of the step, not the reflection of the
    ! rows just advected (alpha = 0 / large for no-slip / free-slip gives exactly -1 / +1)
    PadU(:,1,:) = s_lo*PadU(:,2,:);  PadU(:,nyg,:) = s_hi*PadU(:,nyg-1,:)
    PadW(:,1,:) = s_lo*PadW(:,2,:);  PadW(:,nyg,:) = s_hi*PadW(:,nyg-1,:)
    !$acc end kernels
    Do m = 1, EP
       !$acc kernels present(PadU,PadV,PadW)
       PadU(:,1-m,:)     = s_lo*PadU(:,2+m,:);        PadU(:,nyg+m,:) = s_hi*PadU(:,nyg-1-m,:)
       PadW(:,1-m,:)     = s_lo*PadW(:,2+m,:);        PadW(:,nyg+m,:) = s_hi*PadW(:,nyg-1-m,:)
       PadV(:,1-m,:)     = -PadV(:,1+m,:);             PadV(:,ny+m,:)  = -PadV(:,ny-m,:)
       !$acc end kernels
    End Do

  End Subroutine pad_transported


  !> Project U, V, W onto the divergence-free space of the current density: solve div(beta grad phi) = div u by PCG (at most
  !  maxit iterations, relative residual tol) and subtract beta grad phi from the interior faces, as project_velocity does
  Subroutine vof_project(tol, maxit)

    Real(Int64),    Intent(In) :: tol
    Integer(Int32), Intent(In) :: maxit

    Call vof_fill_pad(Cv, nxg, nyg, nzg)
    Call vp_set_density(Cv)
    Call vp_div(U, V, W, fdiv)
    Call profiler_start(PROF_PCG)
    Call profiler_start(PROF_PCG_PROJECT)
    Call vp_pcg(fdiv, maxit, tol, dphi, dt)
    Call profiler_stop(PROF_PCG_PROJECT)
    Call profiler_stop(PROF_PCG)
    vf_proj_its_last = vp_iters_last
    vf_proj_its_sum = vf_proj_its_sum + vp_iters_last
    Call vp_halo(dphi, .True.)
    Call vp_grad(dphi, Gu_, Gv_, Gw_)
    Call apply_face_gradient(Gu_, Gv_, Gw_, 1d0)

  End Subroutine vof_project




  !> Dynamic viscosity (molecular mixture plus the SGS contribution scaled by the local density) at the cell centres
  Subroutine set_viscosity

    Integer(Int32) :: i, j, k
    Real(Int64) :: mul, mug

    mul = vof_rho_l*vof_nu_l
    mug = vof_rho_g*vof_nu_g
    !$acc parallel loop collapse(3) present(mu_c,mr_c,Cv,vp_rho,nu_t)
    Do k = 1, nzg
       Do j = 1, nyg
          Do i = 1, nxg
             mu_c(i,j,k) = mug + (mul - mug)*Cv(i,j,k) + vp_rho(i,j,k)*nu_t(i,j,k)
             mr_c(i,j,k) = 1d0/mu_c(i,j,k)
          End Do
       End Do
    End Do
    ! the SGS model leaves nu_t at a free-slip wall's ghost row unset: use the nearest interior value there (no-slip keeps zero)
    If ( bc_face_ylo == 2 ) Then
       !$acc parallel loop collapse(2) present(mu_c,mr_c,Cv,vp_rho,nu_t)
       Do k = 1, nzg
          Do i = 1, nxg
             mu_c(i,1,k) = mug + (mul - mug)*Cv(i,1,k) + vp_rho(i,1,k)*nu_t(i,2,k);  mr_c(i,1,k) = 1d0/mu_c(i,1,k)
          End Do
       End Do
    End If
    If ( bc_face_yhi == 2 ) Then
       !$acc parallel loop collapse(2) present(mu_c,mr_c,Cv,vp_rho,nu_t)
       Do k = 1, nzg
          Do i = 1, nxg
             mu_c(i,nyg,k) = mug + (mul - mug)*Cv(i,nyg,k) + vp_rho(i,nyg,k)*nu_t(i,nyg-1,k)
             mr_c(i,nyg,k) = 1d0/mu_c(i,nyg,k)
          End Do
       End Do
    End If

  End Subroutine set_viscosity

  !> Kinematic viscosity of the phase at each cell for the IBM wall model: that of the fluid when the cell and its 26 neighbours hold
  !  one fluid (C within 0.05 of 0 or 1), else 0 (the wall model is off and the stencil flux stays), as for the flat wall in set_wall_stress
  Subroutine set_cell_phase_nu

    Real(Int64), Parameter :: cgate = 0.05d0
    Integer(Int32) :: i, j, k
    Real(Int64), Allocatable :: mnx(:,:,:), mxx(:,:,:), mny(:,:,:), mxy(:,:,:)
    Real(Int64) :: cmn, cmx

    ! min and max over the clipped 3x3x3 window, taken one direction at a time
    Allocate( mnx(nxg,nyg,nzg), mxx(nxg,nyg,nzg), mny(nxg,nyg,nzg), mxy(nxg,nyg,nzg) )
    !$acc data create(mnx,mxx,mny,mxy)
    !$acc parallel loop collapse(3) present(Cv,mnx,mxx)
    Do k = 1, nzg
       Do j = 1, nyg
          Do i = 1, nxg
             mnx(i,j,k) = Min(Cv(i,j,k), Cv(Max(i-1,1),j,k), Cv(Min(i+1,nxg),j,k))
             mxx(i,j,k) = Max(Cv(i,j,k), Cv(Max(i-1,1),j,k), Cv(Min(i+1,nxg),j,k))
          End Do
       End Do
    End Do
    !$acc parallel loop collapse(3) present(mnx,mxx,mny,mxy)
    Do k = 1, nzg
       Do j = 1, nyg
          Do i = 1, nxg
             mny(i,j,k) = Min(mnx(i,j,k), mnx(i,Max(j-1,1),k), mnx(i,Min(j+1,nyg),k))
             mxy(i,j,k) = Max(mxx(i,j,k), mxx(i,Max(j-1,1),k), mxx(i,Min(j+1,nyg),k))
          End Do
       End Do
    End Do
    !$acc parallel loop collapse(3) private(cmn,cmx) present(mny,mxy,nu_cw)
    Do k = 1, nzg
       Do j = 1, nyg
          Do i = 1, nxg
             cmn = Min(mny(i,j,k), mny(i,j,Max(k-1,1)), mny(i,j,Min(k+1,nzg)))
             cmx = Max(mxy(i,j,k), mxy(i,j,Max(k-1,1)), mxy(i,j,Min(k+1,nzg)))
             nu_cw(i,j,k) = 0d0
             If ( cmn >= 1d0 - cgate ) nu_cw(i,j,k) = vof_nu_l
             If ( cmx <= cgate ) nu_cw(i,j,k) = vof_nu_g
          End Do
       End Do
    End Do
    !$acc end data
    Deallocate( mnx, mxx, mny, mxy )

  End Subroutine set_cell_phase_nu


  !> Flat-wall EQWM stress at the y walls: wt_u / wt_w are the tangential stress fluxes (rho u_tau**2 along the local tangential velocity,
  !  signed as the momentum flux through the wall-edge of the first interior row) that viscous_accel uses in place of the molecular
  !  edge flux. u_tau is solved from the first interior row with the viscosity of the phase there. The log law holds in a single
  !  fluid only: where the interface lies within the wall cell, the one above it or one cell sideways the face stays no-slip.
  Subroutine set_wall_stress

    Real(Int64), Parameter :: cgate = 0.05d0
    Integer(Int32) :: i, k, iw, jw, jn, jm, ia, ib, ka, kb, ii, kk, ph(nxg,nzg), pa, pb, ph_keep
    Real(Int64) :: ym, z0w, ut, ux, wz, ur, um, nul, rhof, sgn

    !$acc data create(ph) present(Cv,U,V,W,vp_rfu,vp_rfw,wt_u,wt_w,wt_on_u,wt_on_w)
    !$acc kernels
    wt_on_u = .False.;  wt_on_w = .False.
    !$acc end kernels
    Do iw = 1, 2
       If ( Merge(bc_face_ylo, bc_face_yhi, iw == 1) == 2 ) Cycle
       jw = Merge(2, nyg-1, iw == 1);  jn = Merge(3, nyg-2, iw == 1)
       jm = Merge(j_match_ylo, j_match_yhi, iw == 1);  z0w = Merge(z0_ylo, z0_yhi, iw == 1)
       ym = Merge(yg(jm), Ly - yg(jm), iw == 1);  sgn = Merge(1d0, -1d0, iw == 1)
       !$acc parallel loop collapse(2) default(present)
       Do k = 1, nzg
          Do i = 1, nxg
             ph(i,k) = 0
             If ( Min(Cv(i,jw,k), Cv(i,jn,k), Cv(i,jm,k)) >= 1d0 - cgate ) ph(i,k) = 1
             If ( Max(Cv(i,jw,k), Cv(i,jn,k), Cv(i,jm,k)) <= cgate ) ph(i,k) = -1
          End Do
       End Do
       ! the in-place sweep below depends on the visiting order, so it stays sequential
       !$acc serial default(present)
       Do k = 1, nzg
          ka = Max(k-1, 1);  kb = Min(k+1, nzg)
          Do i = 1, nxg
             ia = Max(i-1, 1);  ib = Min(i+1, nxg)
             ph_keep = ph(i,k)
             Do kk = ka, kb
                Do ii = ia, ib
                   If ( ph(ii,kk) /= ph_keep ) ph(i,k) = 0
                End Do
             End Do
          End Do
       End Do
       !$acc end serial
       !$acc parallel loop collapse(2) default(present) private(pa,pb,ux,wz,ur,nul,um,ut,rhof)
       Do k = 2, nzg-1
          Do i = 2, nx-1
             pa = ph(i,k);  pb = ph(i+1,k)
             If ( pa == 0 .Or. pa /= pb ) Cycle
             ux = U(i,jw,k);  wz = 0.5d0*( W(i,jw,k-1) + W(i,jw,k) );  ur = Sqrt(ux*ux + wz*wz)
             If ( ur == 0d0 ) Cycle
             nul = Merge(vof_nu_l, vof_nu_g, pa == 1)
             um = Sqrt(U(i,jm,k)**2 + (0.5d0*( W(i,jm,k-1) + W(i,jm,k) ))**2)
             Call solve_u_tau_wall(um, ym, z0w, nul, ut)
             rhof = vp_rfu(i,jw,k)
             wt_u(i,k,iw) = sgn*rhof*ut*ut*ux/ur;  wt_on_u(i,k,iw) = .True.
          End Do
       End Do
       !$acc parallel loop collapse(2) default(present) private(pa,pb,ux,wz,ur,nul,um,ut,rhof)
       Do k = 2, nz-1
          Do i = 2, nxg-1
             pa = ph(i,k);  pb = ph(i,k+1)
             If ( pa == 0 .Or. pa /= pb ) Cycle
             wz = W(i,jw,k);  ux = 0.5d0*( U(i-1,jw,k) + U(i,jw,k) );  ur = Sqrt(ux*ux + wz*wz)
             If ( ur == 0d0 ) Cycle
             nul = Merge(vof_nu_l, vof_nu_g, pa == 1)
             um = Sqrt(W(i,jm,k)**2 + (0.5d0*( U(i-1,jm,k) + U(i,jm,k) ))**2)
             Call solve_u_tau_wall(um, ym, z0w, nul, ut)
             rhof = vp_rfw(i,jw,k)
             wt_w(i,k,iw) = sgn*rhof*ut*ut*wz/ur;  wt_on_w(i,k,iw) = .True.
          End Do
       End Do
    End Do
    !$acc end data

  End Subroutine set_wall_stress



  !> Acceleration from the viscous stress, (1/rho_face) div( mu (grad u + grad u^T) ), at the interior faces
  Subroutine viscous_accel(Fu, Fv, Fw)
  ! shear stresses at the edges use the weighted harmonic mean of mu (series layers carry a continuous shear stress), the
  ! arithmetic mean over-weights the dense fluid next to an interface by the viscosity ratio

    Real(Int64), Intent(Out) :: Fu(nx,nyg,nzg), Fv(nxg,ny,nzg), Fw(nxg,nyg,nz)

    Integer(Int32) :: i, j, k
    Real(Int64) :: inv_dx, inv_dx2, mx1, mx2, my1, my2, mz1, mz2, hyj, hzk, ihy, ihz
    Real(Int64) :: dy1, dy2, dyc, dz1, dz2, dzc, fx1, fx2, fy1, fy2, fz1, fz2

    inv_dx = 1d0/dx
    inv_dx2 = inv_dx*inv_dx
    !$acc kernels present(Fu,Fv,Fw)
    Fu = 0d0;  Fv = 0d0;  Fw = 0d0
    !$acc end kernels

    ! ---- u faces
    !$acc parallel loop collapse(3) private(hzk,ihz,hyj,ihy,mx1,mx2,my1,my2,mz1,mz2,fy1,fy2,fz1,fz2) &
    !$acc& present(Fu,U,V,W,mu_c,mr_c,weight_y_0,weight_y_1,weight_z_0,weight_z_1,y,yg,z,zg,wt_on_u,wt_u,vp_rfu,vp_riu,ovr_u)
    Do k = 2, nzg-1
       Do j = 2, nyg-1
          Do i = 2, nx-1
             hzk = z(k) - z(k-1);  ihz = 1d0/hzk
             hyj = y(j) - y(j-1);  ihy = 1d0/hyj
             mx1 = mu_c(i,j,k);  mx2 = mu_c(i+1,j,k)
             my1 = 1d0/( 0.5d0*( weight_y_0(j-1)*( mr_c(i,j-1,k) + mr_c(i+1,j-1,k) ) + weight_y_1(j-1)*( mr_c(i,j,k) &
                   + mr_c(i+1,j,k) ) ) )
             my2 = 1d0/( 0.5d0*( weight_y_0(j)*( mr_c(i,j,k) + mr_c(i+1,j,k) ) + weight_y_1(j)*( mr_c(i,j+1,k) &
                   + mr_c(i+1,j+1,k) ) ) )
             mz1 = 1d0/( 0.5d0*( weight_z_0(k-1)*( mr_c(i,j,k-1) + mr_c(i+1,j,k-1) ) + weight_z_1(k-1)*( mr_c(i,j,k) &
                   + mr_c(i+1,j,k) ) ) )
             mz2 = 1d0/( 0.5d0*( weight_z_0(k)*( mr_c(i,j,k) + mr_c(i+1,j,k) ) + weight_z_1(k)*( mr_c(i,j,k+1) &
                   + mr_c(i+1,j,k+1) ) ) )
             fy2 = my2*( ( U(i,j+1,k) - U(i,j,k) )/( yg(j+1) - yg(j) ) + ( V(i+1,j,k) - V(i,j,k) )*inv_dx )
             fy1 = my1*( ( U(i,j,k) - U(i,j-1,k) )/( yg(j) - yg(j-1) ) + ( V(i+1,j-1,k) - V(i,j-1,k) )*inv_dx )
             If ( j == nyg-1 ) Then;  If ( wt_on_u(i,k,2) ) fy2 = wt_u(i,k,2);  End If
             If ( j == 2 ) Then;  If ( wt_on_u(i,k,1) ) fy1 = wt_u(i,k,1);  End If
             fz2 = mz2*( ( U(i,j,k+1) - U(i,j,k) )/( zg(k+1) - zg(k) ) + ( W(i+1,j,k) - W(i,j,k) )*inv_dx )
             fz1 = mz1*( ( U(i,j,k) - U(i,j,k-1) )/( zg(k) - zg(k-1) ) + ( W(i+1,j,k-1) - W(i,j,k-1) )*inv_dx )
             If ( ibm_wm_on ) Then
                If ( ovr_u(i,j,k,1) < ovr_none ) fy1 = ovr_u(i,j,k,1)*vp_rfu(i,j,k)
                If ( ovr_u(i,j,k,2) < ovr_none ) fy2 = ovr_u(i,j,k,2)*vp_rfu(i,j,k)
                If ( ovr_u(i,j,k,3) < ovr_none ) fz1 = ovr_u(i,j,k,3)*vp_rfu(i,j,k)
                If ( ovr_u(i,j,k,4) < ovr_none ) fz2 = ovr_u(i,j,k,4)*vp_rfu(i,j,k)
             End If
             Fu(i,j,k) = ( 2d0*inv_dx2*( mx2*( U(i+1,j,k) - U(i,j,k) ) - mx1*( U(i,j,k) - U(i-1,j,k) ) )              &
                  + ihy*( fy2 - fy1 ) + ihz*( fz2 - fz1 ) )/vp_riu(i,j,k)
          End Do
       End Do
    End Do

    ! ---- v faces
    !$acc parallel loop collapse(3) private(hzk,ihz,dy1,dy2,dyc,mx1,mx2,my1,my2,mz1,mz2,fx1,fx2,fz1,fz2) &
    !$acc& present(Fv,U,V,W,mu_c,mr_c,weight_y_0,weight_y_1,weight_z_0,weight_z_1,y,yg,z,zg,vp_rfv,vp_riv,ovr_v)
    Do k = 2, nzg-1
       Do j = 2, ny-1
          Do i = 2, nxg-1
             hzk = z(k) - z(k-1);  ihz = 1d0/hzk
             dy1 = y(j) - y(j-1);  dy2 = y(j+1) - y(j);  dyc = yg(j+1) - yg(j)
             mx1 = 1d0/( 0.5d0*( weight_y_0(j)*( mr_c(i-1,j,k) + mr_c(i,j,k) ) + weight_y_1(j)*( mr_c(i-1,j+1,k) &
                   + mr_c(i,j+1,k) ) ) )
             mx2 = 1d0/( 0.5d0*( weight_y_0(j)*( mr_c(i,j,k) + mr_c(i+1,j,k) ) + weight_y_1(j)*( mr_c(i,j+1,k) &
                   + mr_c(i+1,j+1,k) ) ) )
             my1 = mu_c(i,j,k);  my2 = mu_c(i,j+1,k)
             mz1 = 1d0/( weight_z_0(k-1)*( weight_y_0(j)*mr_c(i,j,k-1) + weight_y_1(j)*mr_c(i,j+1,k-1) ) &
                   + weight_z_1(k-1)*( weight_y_0(j)*mr_c(i,j,k) + weight_y_1(j)*mr_c(i,j+1,k) ) )
             mz2 = 1d0/( weight_z_0(k)*( weight_y_0(j)*mr_c(i,j,k) + weight_y_1(j)*mr_c(i,j+1,k) ) &
                   + weight_z_1(k)*( weight_y_0(j)*mr_c(i,j,k+1) + weight_y_1(j)*mr_c(i,j+1,k+1) ) )
             fx2 = mx2*( ( V(i+1,j,k) - V(i,j,k) )*inv_dx + ( U(i,j+1,k) - U(i,j,k) )/dyc )
             fx1 = mx1*( ( V(i,j,k) - V(i-1,j,k) )*inv_dx + ( U(i-1,j+1,k) - U(i-1,j,k) )/dyc )
             fz2 = mz2*( ( V(i,j,k+1) - V(i,j,k) )/( zg(k+1) - zg(k) ) + ( W(i,j+1,k) - W(i,j,k) )/dyc )
             fz1 = mz1*( ( V(i,j,k) - V(i,j,k-1) )/( zg(k) - zg(k-1) ) + ( W(i,j+1,k-1) - W(i,j,k-1) )/dyc )
             If ( ibm_wm_on ) Then
                If ( ovr_v(i,j,k,1) < ovr_none ) fx1 = ovr_v(i,j,k,1)*vp_rfv(i,j,k)
                If ( ovr_v(i,j,k,2) < ovr_none ) fx2 = ovr_v(i,j,k,2)*vp_rfv(i,j,k)
                If ( ovr_v(i,j,k,3) < ovr_none ) fz1 = ovr_v(i,j,k,3)*vp_rfv(i,j,k)
                If ( ovr_v(i,j,k,4) < ovr_none ) fz2 = ovr_v(i,j,k,4)*vp_rfv(i,j,k)
             End If
             Fv(i,j,k) = ( inv_dx*( fx2 - fx1 )                                                                       &
                  + 2d0/dyc*( my2*( V(i,j+1,k) - V(i,j,k) )/dy2 - my1*( V(i,j,k) - V(i,j-1,k) )/dy1 )                 &
                  + ihz*( fz2 - fz1 ) )/vp_riv(i,j,k)
          End Do
       End Do
    End Do

    ! ---- w faces
    !$acc parallel loop collapse(3) private(dz1,dz2,dzc,hyj,ihy,mx1,mx2,my1,my2,mz1,mz2,fx1,fx2,fy1,fy2) &
    !$acc& present(Fw,U,V,W,mu_c,mr_c,weight_y_0,weight_y_1,weight_z_0,weight_z_1,y,yg,z,zg,wt_on_w,wt_w,vp_rfw,vp_riw,ovr_w)
    Do k = 2, nz-1
       Do j = 2, nyg-1
          Do i = 2, nxg-1
             dz1 = z(k) - z(k-1);  dz2 = z(k+1) - z(k);  dzc = zg(k+1) - zg(k)
             hyj = y(j) - y(j-1);  ihy = 1d0/hyj
             mx1 = 1d0/( 0.5d0*( weight_z_0(k)*( mr_c(i-1,j,k) + mr_c(i,j,k) ) + weight_z_1(k)*( mr_c(i-1,j,k+1) &
                   + mr_c(i,j,k+1) ) ) )
             mx2 = 1d0/( 0.5d0*( weight_z_0(k)*( mr_c(i,j,k) + mr_c(i+1,j,k) ) + weight_z_1(k)*( mr_c(i,j,k+1) &
                   + mr_c(i+1,j,k+1) ) ) )
             my1 = 1d0/( weight_z_0(k)*( weight_y_0(j-1)*mr_c(i,j-1,k) + weight_y_1(j-1)*mr_c(i,j,k) ) &
                   + weight_z_1(k)*( weight_y_0(j-1)*mr_c(i,j-1,k+1) + weight_y_1(j-1)*mr_c(i,j,k+1) ) )
             my2 = 1d0/( weight_z_0(k)*( weight_y_0(j)*mr_c(i,j,k) + weight_y_1(j)*mr_c(i,j+1,k) ) &
                   + weight_z_1(k)*( weight_y_0(j)*mr_c(i,j,k+1) + weight_y_1(j)*mr_c(i,j+1,k+1) ) )
             mz1 = mu_c(i,j,k);  mz2 = mu_c(i,j,k+1)
             fy2 = my2*( ( W(i,j+1,k) - W(i,j,k) )/( yg(j+1) - yg(j) ) + ( V(i,j,k+1) - V(i,j,k) )/dzc )
             fy1 = my1*( ( W(i,j,k) - W(i,j-1,k) )/( yg(j) - yg(j-1) ) + ( V(i,j-1,k+1) - V(i,j-1,k) )/dzc )
             If ( j == nyg-1 ) Then;  If ( wt_on_w(i,k,2) ) fy2 = wt_w(i,k,2);  End If
             If ( j == 2 ) Then;  If ( wt_on_w(i,k,1) ) fy1 = wt_w(i,k,1);  End If
             fx2 = mx2*( ( W(i+1,j,k) - W(i,j,k) )*inv_dx + ( U(i,j,k+1) - U(i,j,k) )/dzc )
             fx1 = mx1*( ( W(i,j,k) - W(i-1,j,k) )*inv_dx + ( U(i-1,j,k+1) - U(i-1,j,k) )/dzc )
             If ( ibm_wm_on ) Then
                If ( ovr_w(i,j,k,1) < ovr_none ) fx1 = ovr_w(i,j,k,1)*vp_rfw(i,j,k)
                If ( ovr_w(i,j,k,2) < ovr_none ) fx2 = ovr_w(i,j,k,2)*vp_rfw(i,j,k)
                If ( ovr_w(i,j,k,3) < ovr_none ) fy1 = ovr_w(i,j,k,3)*vp_rfw(i,j,k)
                If ( ovr_w(i,j,k,4) < ovr_none ) fy2 = ovr_w(i,j,k,4)*vp_rfw(i,j,k)
             End If
             Fw(i,j,k) = ( inv_dx*( fx2 - fx1 )                                                                       &
                  + ihy*( fy2 - fy1 )  &
                  + 2d0/dzc*( mz2*( W(i,j,k+1) - W(i,j,k) )/dz2 - mz1*( W(i,j,k) - W(i,j,k-1) )/dz1 )                 &
                  )/vp_riw(i,j,k)
          End Do
       End Do
    End Do

  End Subroutine viscous_accel


  !> Balanced-force surface tension acceleration sigma kappa grad(C)/rho on the faces (the face gradient of C, the mean of the
  !  curvatures of the interface cells next to the face, the same face density as the pressure gradient)
  Subroutine add_surface_tension(Fx, Fy, Fz)

    Real(Int64), Intent(InOut) :: Fx(nx,nyg,nzg), Fy(nxg,ny,nzg), Fz(nxg,nyg,nz)
    Integer(Int32) :: i, j, k

    !$acc parallel loop collapse(3) present(Fx,Cv,kap_c,vp_riu)
    Do k = 2, nzg-1
       Do j = 2, nyg-1
          Do i = 2, nx-1
             Fx(i,j,k) = Fx(i,j,k) + vof_sigma*face_kappa(Cv(i,j,k), Cv(i+1,j,k), kap_c(i,j,k), kap_c(i+1,j,k)) &
                         *( Cv(i+1,j,k) - Cv(i,j,k) )/( dx*vp_riu(i,j,k) )
          End Do
       End Do
    End Do
    !$acc parallel loop collapse(3) present(Fy,Cv,kap_c,vp_riv,yg)
    Do k = 2, nzg-1
       Do j = 2, ny-1
          Do i = 2, nxg-1
             Fy(i,j,k) = Fy(i,j,k) + vof_sigma*face_kappa(Cv(i,j,k), Cv(i,j+1,k), kap_c(i,j,k), kap_c(i,j+1,k)) &
                         *( Cv(i,j+1,k) - Cv(i,j,k) )/( ( yg(j+1) - yg(j) )*vp_riv(i,j,k) )
          End Do
       End Do
    End Do
    !$acc parallel loop collapse(3) present(Fz,Cv,kap_c,vp_riw,zg)
    Do k = 2, nz-1
       Do j = 2, nyg-1
          Do i = 2, nxg-1
             Fz(i,j,k) = Fz(i,j,k) + vof_sigma*face_kappa(Cv(i,j,k), Cv(i,j,k+1), kap_c(i,j,k), kap_c(i,j,k+1)) &
                         *( Cv(i,j,k+1) - Cv(i,j,k) )/( ( zg(k+1) - zg(k) )*vp_riw(i,j,k) )
          End Do
       End Do
    End Do

  End Subroutine add_surface_tension


  Function face_kappa(c1, c2, kap1, kap2) Result(kf)
    !$acc routine seq

    Real(Int64), Intent(In) :: c1, c2, kap1, kap2
    Real(Int64) :: kf, a1, a2

    a1 = Merge(1d0, 0d0, c1 > vof_eps .And. c1 < 1d0 - vof_eps)
    a2 = Merge(1d0, 0d0, c2 > vof_eps .And. c2 < 1d0 - vof_eps)
    kf = 0d0
    If ( a1 + a2 > 0d0 ) kf = ( a1*kap1 + a2*kap2 )/( a1 + a2 )

  End Function face_kappa


  !> Wray RK3 stages with frozen density: explicit viscous stress and hydrostatic-reduced gravity, pressure predictor + PCG increment
  Subroutine vof_forces_rk3(to)

    Real(Int64), Intent(In) :: to
    Integer(Int32) :: s, i, j, k
    Real(Int64) :: ass, a1, a2, dta

    Call vof_fill_pad(Cv, nxg, nyg, nzg)
    Call vp_set_density(Cv)
    Call profiler_start(PROF_F_BODY)
    If ( vof_sigma > 0d0 ) Call vof_curvature(kap_c)
    Call profiler_stop(PROF_F_BODY)
    ! Cv does not change in the force stages
    Call profiler_start(PROF_F_WALL)
    If ( ibm_wm_on ) Call set_cell_phase_nu
    Call profiler_stop(PROF_F_WALL)

    Do s = 1, 3
       rk_step = s
       Call profiler_start(PROF_F_SGS)
       Call compute_sgs_model(U, V, W, nu_t)
       Call profiler_stop(PROF_F_SGS)
       Call profiler_start(PROF_F_WALL)
       If ( ibm_wm_on ) Then
          Call ibm_stress_update(U, V, W, vp_mu, vp_mv, vp_mw, nu_cw)
       End If
       Call compute_wall_model(U, V, W)
       Call set_viscosity
       If ( flat_wall_model_flag > 0 ) Call set_wall_stress
       Call profiler_stop(PROF_F_WALL)
       Call profiler_start(PROF_F_VISC)
       Call viscous_accel(Fu_, Fv_, Fw_)
       Call profiler_stop(PROF_F_VISC)
       Call profiler_start(PROF_F_BODY)
       If ( vof_hsplit == 0 ) Then
          !$acc parallel loop collapse(3) present(Fv_,vp_rfv,vp_riv,rsf)
          Do k = 2, nzg-1
             Do j = 2, ny-1
                Do i = 2, nxg-1
                   Fv_(i,j,k) = Fv_(i,j,k) - vof_grav*( vp_rfv(i,j,k) - rsf(j) )/vp_riv(i,j,k)
                End Do
             End Do
          End Do
       Else
          !$acc parallel loop collapse(3) present(Fu_,vp_rho,vp_riu,yg)
          Do k = 2, nzg-1
             Do j = 2, nyg-1
                Do i = 2, nx-1
                   Fu_(i,j,k) = Fu_(i,j,k) + vof_grav*( yg(j) - vof_level )*( vp_rho(i+1,j,k) - vp_rho(i,j,k) )/( dx*vp_riu(i,j,k) )
                End Do
             End Do
          End Do
          !$acc parallel loop collapse(3) present(Fv_,vp_rho,vp_riv,y,yg)
          Do k = 2, nzg-1
             Do j = 2, ny-1
                Do i = 2, nxg-1
                   Fv_(i,j,k) = Fv_(i,j,k) + vof_grav*( y(j) - vof_level )*( vp_rho(i,j+1,k) - vp_rho(i,j,k) ) &
                                /( ( yg(j+1) - yg(j) )*vp_riv(i,j,k) )
                End Do
             End Do
          End Do
          !$acc parallel loop collapse(3) present(Fw_,vp_rho,vp_riw,yg,zg)
          Do k = 2, nz-1
             Do j = 2, nyg-1
                Do i = 2, nxg-1
                   Fw_(i,j,k) = Fw_(i,j,k) + vof_grav*( yg(j) - vof_level )*( vp_rho(i,j,k+1) - vp_rho(i,j,k) ) &
                                /( ( zg(k+1) - zg(k) )*vp_riw(i,j,k) )
                End Do
             End Do
          End Do
       End If

       ! Boussinesq buoyancy of the liquid only; the sweeps and the pressure operator keep the unperturbed phase densities
       If ( boussinesq_flag >= 1 ) Then
          !$acc parallel loop collapse(3) present(Fv_,Tscal,Cv,vp_riv,weight_y_0,weight_y_1)
          Do k = 2, nzg-1
             Do j = 2, ny-1
                Do i = 2, nxg-1
                   Fv_(i,j,k) = Fv_(i,j,k) + vof_grav*beta_T*vof_rho_l &
                        *( weight_y_0(j)*(Tscal(i,j,k) - T_ref) + weight_y_1(j)*(Tscal(i,j+1,k) - T_ref) ) &
                        *0.5d0*( Min(1d0, Max(0d0, Cv(i,j,k))) + Min(1d0, Max(0d0, Cv(i,j+1,k))) )/vp_riv(i,j,k)
                End Do
             End Do
          End Do
       End If

       If ( vof_sigma > 0d0 ) Call add_surface_tension(Fu_, Fv_, Fw_)
       Call profiler_stop(PROF_F_BODY)

       Call profiler_start(PROF_F_PRED)
       ass = rk_coef(s,s)
       dta = dt*ass
       ! stage predictor u* = u_o + dt sum_{j<s} a_sj H_j + dt a_ss F_s  (interior faces; seams and walls come from the BCs)
       a1 = rk_coef(s,1);  a2 = rk_coef(s,2)
       !$acc parallel loop collapse(3) present(U,Uo,Fu_,Hu1,Hu2)
       Do k = 2, nzg-1
          Do j = 2, nyg-1
             Do i = 2, nx-1
                U(i,j,k) = Uo(i,j,k) + dta*Fu_(i,j,k)
                If ( s >= 2 ) U(i,j,k) = U(i,j,k) + dt*a1*Hu1(i,j,k)
                If ( s == 3 ) U(i,j,k) = U(i,j,k) + dt*a2*Hu2(i,j,k)
             End Do
          End Do
       End Do
       !$acc parallel loop collapse(3) present(V,Vo,Fv_,Hv1,Hv2)
       Do k = 2, nzg-1
          Do j = 2, ny-1
             Do i = 2, nxg-1
                V(i,j,k) = Vo(i,j,k) + dta*Fv_(i,j,k)
                If ( s >= 2 ) V(i,j,k) = V(i,j,k) + dt*a1*Hv1(i,j,k)
                If ( s == 3 ) V(i,j,k) = V(i,j,k) + dt*a2*Hv2(i,j,k)
             End Do
          End Do
       End Do
       !$acc parallel loop collapse(3) present(W,Wo,Fw_,Hw1,Hw2)
       Do k = 2, nz-1
          Do j = 2, nyg-1
             Do i = 2, nxg-1
                W(i,j,k) = Wo(i,j,k) + dta*Fw_(i,j,k)
                If ( s >= 2 ) W(i,j,k) = W(i,j,k) + dt*a1*Hw1(i,j,k)
                If ( s == 3 ) W(i,j,k) = W(i,j,k) + dt*a2*Hw2(i,j,k)
             End Do
          End Do
       End Do
       t = to + rk_t(s)*dt

       outflow_relax_on = .True.
       Call apply_boundary_conditions
       If ( ibm_input_mode >= 1 ) Call enforce_ibm
       ! explicit predictor force, then the increment
       Call vp_halo(ppre, .True.)
       Call vp_grad(ppre, Bu_, Bv_, Bw_)
       Call apply_face_gradient(Bu_, Bv_, Bw_, dta)
       Call face_halo(U, V, W)
       Call vp_div(U, V, W, fdiv)
       !$acc kernels present(fdiv)
       fdiv = fdiv/dta
       !$acc end kernels
       Call profiler_stop(PROF_F_PRED)
       Call profiler_start(PROF_PCG)
       Call profiler_start(PROF_PCG_FORCES)
       Call vp_pcg(fdiv, Max(vof_pcg_iters, 1), vof_pcg_tol, dphi, dt*dta)
       Call profiler_stop(PROF_PCG_FORCES)
       Call profiler_stop(PROF_PCG)
       Call profiler_start(PROF_F_POST)
       Call vp_halo(dphi, .True.)
       Call vp_grad(dphi, Gu_, Gv_, Gw_)
       Call apply_face_gradient(Gu_, Gv_, Gw_, dta)
       ! H_s = F_s - beta grad(p_s), with beta grad(p_s) = b + g
       !$acc kernels present(ppre,dphi,Fu_,Fv_,Fw_,Bu_,Bv_,Bw_,Gu_,Gv_,Gw_,Hu1,Hv1,Hw1,Hu2,Hv2,Hw2)
       ppre = ppre + dphi
       If ( s == 1 ) Then
          Hu1 = Fu_ - (Bu_ + Gu_);  Hv1 = Fv_ - (Bv_ + Gv_);  Hw1 = Fw_ - (Bw_ + Gw_)
       Else If ( s == 2 ) Then
          Hu2 = Fu_ - (Bu_ + Gu_);  Hv2 = Fv_ - (Bv_ + Gv_);  Hw2 = Fw_ - (Bw_ + Gw_)
       End If
       !$acc end kernels
       Call apply_boundary_conditions(after_projection=.True.)
       If ( ibm_input_mode >= 1 ) Call enforce_ibm
       Call profiler_stop(PROF_F_POST)
    End Do
    outflow_relax_on = .False.

  End Subroutine vof_forces_rk3



  !> waves2Foam-style strength of a relaxation zone, r = 1 at the zone's outer boundary, 0 at its inner edge
  Function relax_strength(r) Result(f)

    Real(Int64), Intent(In) :: r
    Real(Int64) :: f

    f = ( Exp(Min(Max(r, 0d0), 1d0)**3.5d0) - 1d0 )/( Exp(1d0) - 1d0 )

  End Function relax_strength


  !> Zone strength at the streamwise position xs and whether it is the generation zone (target: the wave) or the absorption zone
  Subroutine zone_at(xs, strength, generation)

    Real(Int64), Intent(In)  :: xs
    Real(Int64), Intent(Out) :: strength
    Logical,     Intent(Out) :: generation

    strength = 0d0
    generation = .False.
    If ( wave_gen_len > 0d0 .And. xs < wave_gen_len ) Then
       strength = relax_strength(1d0 - xs/wave_gen_len)
       generation = .True.
    Else If ( wave_abs_len > 0d0 .And. xs > Lx_i - wave_abs_len ) Then
       strength = relax_strength((xs - (Lx_i - wave_abs_len))/wave_abs_len)
    End If

  End Subroutine zone_at


  !> Relax C and the velocity toward the target (wave in the generation zone, still water at rest in the absorption zone) over dts;
  !  the liquid volume added to C is tallied for the ledger
  Subroutine vof_relax(dts)

    Real(Int64), Intent(In) :: dts

    Integer(Int32) :: i, j, k
    Real(Int64) :: strength, fr, top, ct, eta, vol
    Real(Int64), Allocatable :: uc(:), vf(:)
    Logical :: gen
    Logical :: is_first_z, is_last_z
    Integer(Int32) :: partner, khi

    If ( wave_gen_len <= 0d0 .And. wave_abs_len <= 0d0 ) Return
    ! the wave targets are evaluated on the host per x station and handed to the kernels as tables
    If ( .Not. Allocated(rl_fr_c) ) Then
       Allocate( rl_fr_c(nxg), rl_top(nxg), rl_vf(ny,nxg), rl_fr_u(nx), rl_uc(nyg,nx) )
       !$acc enter data create(rl_fr_c,rl_top,rl_vf,rl_fr_u,rl_uc)
    End If
    Allocate( uc(nyg), vf(ny) )
    khi = nzg-1
    If ( z_bc_type == 0 ) Then
       Call z_periodic_partner(is_first_z, is_last_z, partner)
       If ( is_last_z ) khi = nzg-2
    End If

    rl_fr_c = 0d0;  rl_top = vof_level;  rl_vf = 0d0;  rl_fr_u = 0d0;  rl_uc = 0d0
    Do i = 2, nxg-1
       Call zone_at(xg(i), strength, gen)
       If ( strength <= 0d0 ) Cycle
       rl_fr_c(i) = 1d0 - Exp(-wave_relax_rate*strength*dts)
       If ( gen .And. wave_type > 0 ) Then
          rl_top(i) = vof_level + wave_eta(xg(i), t)
          Call wave_profile_at(xg(i), t, eta, uc, vf)
          rl_vf(:,i) = vf
       End If
    End Do
    Do i = 2, nx-1
       Call zone_at(x(i), strength, gen)
       If ( strength <= 0d0 ) Cycle
       rl_fr_u(i) = 1d0 - Exp(-wave_relax_rate*strength*dts)
       If ( gen .And. wave_type > 0 ) Then
          Call wave_profile_at(x(i), t, eta, uc, vf)
       Else If ( current_type > 0 ) Then
          Call current_profile(t, uc)
       Else
          uc = 0d0
       End If
       rl_uc(:,i) = uc
    End Do
    !$acc update device(rl_fr_c,rl_top,rl_vf,rl_fr_u,rl_uc)

    If ( vp_masked ) Then
       !$acc parallel loop collapse(3) default(present) private(fr,ct,vol) reduction(+:vf_relax_loc)
       Do i = 2, nxg-1
          Do j = 2, nyg-1
             Do k = 2, nzg-1
                fr = rl_fr_c(i)
                ct = Min(1d0, Max(0d0, (rl_top(i) - y(j-1))/(y(j) - y(j-1))))
                If ( fr > 0d0 .And. vp_act(i,j,k) >= 0.5d0 ) Then   ! IBM solid cell: its liquid fraction is a continuation, not a state
                   vol = Cv(i,j,k)
                   Cv(i,j,k) = Cv(i,j,k) + fr*(ct - Cv(i,j,k))
                   If ( k <= khi ) vf_relax_loc = vf_relax_loc + (Cv(i,j,k) - vol)*dx*hy(j)*hz(k)
                End If
             End Do
          End Do
       End Do
    Else
       !$acc parallel loop collapse(3) default(present) private(fr,ct,vol) reduction(+:vf_relax_loc)
       Do i = 2, nxg-1
          Do j = 2, nyg-1
             Do k = 2, nzg-1
                fr = rl_fr_c(i)
                ct = Min(1d0, Max(0d0, (rl_top(i) - y(j-1))/(y(j) - y(j-1))))
                If ( fr > 0d0 ) Then
                   vol = Cv(i,j,k)
                   Cv(i,j,k) = Cv(i,j,k) + fr*(ct - Cv(i,j,k))
                   If ( k <= khi ) vf_relax_loc = vf_relax_loc + (Cv(i,j,k) - vol)*dx*hy(j)*hz(k)
                End If
             End Do
          End Do
       End Do
    End If
    !$acc parallel loop collapse(3) default(present)
    Do i = 2, nxg-1
       Do k = 2, nzg-1
          Do j = 2, ny-1
             V(i,j,k) = V(i,j,k) + rl_fr_c(i)*(rl_vf(j,i) - V(i,j,k))
          End Do
       End Do
    End Do
    !$acc parallel loop collapse(3) default(present)
    Do i = 2, nxg-1
       Do k = 2, nz-1
          Do j = 2, nyg-1
             W(i,j,k) = W(i,j,k)*(1d0 - rl_fr_c(i))
          End Do
       End Do
    End Do
    !$acc parallel loop collapse(3) default(present)
    Do i = 2, nx-1
       Do k = 2, nzg-1
          Do j = 2, nyg-1
             U(i,j,k) = U(i,j,k) + rl_fr_u(i)*(rl_uc(j,i) - U(i,j,k))
          End Do
       End Do
    End Do
    Deallocate( uc, vf )

  End Subroutine vof_relax


  !> Fold this step's rank-local boundary-flux and relaxation tallies into the global cumulative ledger
  Subroutine tally_ledger

    Real(Int64) :: loc(2), glb(2)

    loc = (/ vf_bnd_loc, vf_relax_loc /)
    Call MPI_Allreduce(loc, glb, 2, MPI_real8, MPI_SUM, MPI_COMM_WORLD, ierr)
    vof_bnd_cum = vof_bnd_cum + glb(1)
    vof_relax_cum = vof_relax_cum + glb(2)
    vf_bnd_loc = 0d0;  vf_relax_loc = 0d0

  End Subroutine tally_ledger


  !> Explicit SSP-RK3 step of dts for one cell-centred scalar (kind 1 sediment with settling, kind 2 temperature) with the end-of-step
  !  velocity. The flux form -div(u C) + C div(u) keeps a uniform field exactly uniform whatever divergence the velocity has; the
  !  molecular diffusivity is nu/Sc (nu/Pr) of the phase at each cell blended by the liquid fraction, harmonic on the faces.
  Subroutine vof_scalar_rk3(Cs, kind, dts)

    Real(Int64), Intent(InOut) :: Cs(nxg,nyg,nzg)
    Integer(Int32), Intent(In) :: kind
    Real(Int64), Intent(In) :: dts

    Integer(Int32) :: stage, i, j, k
    Real(Int64) :: a0, divu, kl, kg, fr, cf1, cf2, wset

    If ( .Not. Allocated(Fsc) ) Then
       Allocate( Fsc(2:nxg-1,2:nyg-1,2:nzg-1), Vsc(nxg,ny,nzg), Csc0(nxg,nyg,nzg) )
       If ( Allocated(kfield) ) Then
          !$acc exit data delete(kfield)
          Deallocate( kfield )
       End If
       Allocate( kfield(nxg,nyg,nzg) )
       !$acc enter data create(kfield,Fsc,Vsc,Csc0)
    End If
    kl = vof_nu_l/Merge(Sc, Pr, kind == 1);  kg = vof_nu_g/Merge(Sc, Pr, kind == 1)
    wset = Merge(ws, 0d0, kind == 1)
    !$acc parallel loop collapse(3) present(kfield,Cv)
    Do k = 1, nzg
       Do j = 1, nyg
          Do i = 1, nxg
             kfield(i,j,k) = kg + (kl - kg)*Min(1d0, Max(0d0, Cv(i,j,k)))
          End Do
       End Do
    End Do
    kf_on = .True.
    ! settling acts in the liquid only
    !$acc parallel loop collapse(3) present(Vsc,V,Cv) private(cf1,cf2,fr)
    Do k = 1, nzg
       Do j = 1, ny
          Do i = 1, nxg
             cf1 = Min(1d0, Max(0d0, Cv(i,j,k)));  cf2 = Min(1d0, Max(0d0, Cv(i,j+1,k)))
             fr = 0.5d0*( cf1 + cf2 )
             Vsc(i,j,k) = V(i,j,k) - wset*fr
          End Do
       End Do
    End Do
    !$acc kernels present(Csc0,Cs)
    Csc0 = Cs
    !$acc end kernels
    Do stage = 1, 3
       Call compute_rhs_scalar_core(Cs, U, Vsc, W, 0d0, 0d0, Merge(1d0/Sc_t, 1d0/Pr_t, kind == 1), kind == 2, Fsc)
       a0 = Merge(0.75d0, 1d0/3d0, stage == 2)
       ! written as (u + dt F) + a0 (u0 - (u + dt F)) so that a uniform field with F = 0 is reproduced exactly
       !$acc parallel loop collapse(3) present(Cs,Csc0,Fsc,U,V,W,x,y,z) private(divu)
       Do k = 2, nzg-1
          Do j = 2, nyg-1
             Do i = 2, nxg-1
                divu = ( U(i,j,k) - U(i-1,j,k) )/( x(i) - x(i-1) ) + ( V(i,j,k) - V(i,j-1,k) )/( y(j) - y(j-1) ) &
                     + ( W(i,j,k) - W(i,j,k-1) )/( z(k) - z(k-1) )
                Cs(i,j,k) = Cs(i,j,k) + dts*( Fsc(i,j,k) + Cs(i,j,k)*divu )
                If ( stage > 1 ) Cs(i,j,k) = Cs(i,j,k) + a0*( Csc0(i,j,k) - Cs(i,j,k) )
             End Do
          End Do
       End Do
       If ( kind == 1 ) Then
          Call apply_scalar_bc(Cs)
       Else
          Call apply_temperature_bc(Cs)
       End If
       If ( ibm_input_mode >= 1 ) Then
          If ( kind == 1 ) Then
             Call apply_ghost_cell_ibm_scalar_noflux(Cs)
          Else
             Call apply_ghost_cell_ibm_scalar(Cs)
          End If
       End If
    End Do
    kf_on = .False.

  End Subroutine vof_scalar_rk3


  !> One two-fluid time step: dt from the CFL (the Courant limit sees the gas currents too), then the Strang sequence
  Subroutine compute_time_step_vof

    Real(Int64) :: to, cfl_conv, cfl_visc, cfl_accel, dt_presnap, dt_in, hmin_w, om2
    Integer(Int32) :: isub, nsub

    ibm_sampling_now = ( ibm_input_mode >= 1 .And. nsampling > 0 )
    If ( ibm_sampling_now ) ibm_sampling_now = ( Mod(istep, nsampling) == 0 )
    If ( ibm_input_mode >= 1 ) Call enforce_ibm
    Call compute_sgs_model(U, V, W, nu_t)
    Call compute_wall_model(U, V, W)
    Call compute_cfl(cfl_conv, cfl_visc, cfl_accel)
    ! the temperature diffuses faster than the momentum when Pr < 1
    If ( boussinesq_flag >= 1 .And. Pr < 1d0 ) cfl_visc = cfl_visc/Pr
    cfl_conv_last  = cfl_conv
    cfl_visc_last  = cfl_visc
    cfl_accel_last = cfl_accel
    ! shortest resolvable gravity-capillary wave (k = pi/h, omega^2 = g k + sigma k^3/(rho_l + rho_g)): dt stays below (pi/2)/omega,
    ! the form of Brackbill's capillary limit sqrt((rho_l + rho_g) h^3/(4 pi sigma)) that this reduces to for g = 0. Without it the
    ! step grows with the (small) flow velocities of a nearly still free surface until the explicit interface coupling goes unstable.
    If ( vof_sigma > 0d0 .Or. vof_grav > 0d0 ) Then
       hmin_w = Min(dx, dymin, dzmin)
       om2 = vof_grav*4d0*Atan(1d0)/hmin_w + vof_sigma*(4d0*Atan(1d0))**3/((vof_rho_l + vof_rho_g)*hmin_w**3)
       cfl_accel = Max(cfl_accel, dt*Sqrt(om2)/(2d0*Atan(1d0)))
       cfl_accel_last = cfl_accel
    End If
    cfl_current    = Max(cfl_conv, cfl_visc, cfl_accel)
    dt_in = dt
    Call select_dt(cfl_current, dt_presnap)
    to = t
    ! the convective outflow relaxation runs once per RK stage inside vof_forces_rk3 (stage shares of dt sum to dt)
    outflow_relax_on = .False.

    nsub = vof_nsub
    If ( nsub <= 0 ) nsub = Max(1, Ceiling( 0.5d0*cfl_conv*dt/dt_in/vof_co_sub ))
    vof_nsub_last = nsub
    Do isub = 1, nsub
       Call vof_advect_half(0.5d0*dt/nsub, keep_transport=( isub > 1 .And. vof_freeze_ut == 1 ))
    End Do
    Call apply_boundary_conditions
    Call profiler_start(PROF_VOF_PROJECT)
    Call vof_project(vof_adv_tol, vof_adv_iters)
    Call profiler_stop(PROF_VOF_PROJECT)
    Call apply_boundary_conditions(after_projection=.True.)
    If ( ibm_input_mode >= 1 ) Call enforce_ibm

    Call dev_copy(Uo, U)
    Call dev_copy(Vo, V)
    Call dev_copy(Wo, W)
    Call profiler_start(PROF_VOF_FORCES)
    Call vof_forces_rk3(to)
    Call profiler_stop(PROF_VOF_FORCES)

    nsub = vof_nsub
    If ( nsub <= 0 ) nsub = Max(1, Ceiling( 0.5d0*cfl_conv*dt/dt_in/vof_co_sub ))
    Do isub = 1, nsub
       Call vof_advect_half(0.5d0*dt/nsub, keep_transport=( isub > 1 .And. vof_freeze_ut == 1 ))
    End Do
    t = to + dt
    Call vof_relax(dt)
    Call apply_boundary_conditions
    Call profiler_start(PROF_VOF_PROJECT)
    Call vof_project(vof_adv_tol, vof_adv_iters)
    Call profiler_stop(PROF_VOF_PROJECT)
    Call apply_boundary_conditions(after_projection=.True.)
    If ( ibm_input_mode >= 1 ) Call enforce_ibm
    Call profiler_start(PROF_VOF_END)
    Call vof_fill_pad(Cv, nxg, nyg, nzg)
    If ( particles_active >= 1 ) Then
       ! the particles are host code: they read the host velocity and liquid fraction
       !$acc update host(U,V,W,Cv)
       Call particles_load_surface(Cv(1:nxg,1:nyg,1:nzg))
       Call advance_particles
    End If
    If ( sediment_flag >= 1 ) Call vof_scalar_rk3(Cscal, 1, dt)
    If ( boussinesq_flag >= 1 ) Call vof_scalar_rk3(Tscal, 2, dt)
    Call tally_ledger
    ! the whole array, ghost planes included: snapshots copy every rank's full local block, so a stale (zero) seam ghost plane would
    ! overwrite the neighbour's last real plane in the file, and a restart reads that file back as the pressure predictor
    Call dev_copy(P, ppre)
    !$acc kernels present(P)
    P(:,1,:) = P(:,2,:);  P(:,nyg,:) = P(:,nyg-1,:)
    !$acc end kernels
    Call vof_end_step
    If ( Mod(istep, nmonitor) == 0 ) Call vof_output_monitor
    Call vof_write_gauges
    outflow_relax_on = .True.
    If ( ibm_sampling_now ) Then
       Block
          Real(Int64) :: fi(3), fp(3), fv(3)
          Call vof_compute_ibm_forces(fi(1), fi(2), fi(3), fp(1), fp(2), fp(3), fv(1), fv(2), fv(3))
          Call write_force_csv(fi(1), fi(2), fi(3), fp(1), fp(2), fp(3), fv(1), fv(2), fv(3))
       End Block
       ibm_sampling_now = .False.
    End If
    Call profiler_stop(PROF_VOF_END)
    dt_step = dt
    dt = dt_presnap

  End Subroutine compute_time_step_vof

End Module vof_twofluid
