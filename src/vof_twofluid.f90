!> Two-fluid time step on the staggered grid (VOF liquid fraction C, one-fluid momentum). Validated design from tests/vof/spike_c_momentum.f90:
!>
!>   Strang step:  A(dt/2) - project - RK3(forces, frozen density) - A(dt/2) - project
!>
!>   A(tau)  advances C with the Weymouth-Yue sweeps and the momentum rho*u with the same mass fluxes, so mass and momentum are
!>           transported consistently (uniform velocity is reproduced exactly at any density ratio); the c-tilde*div(u) term of
!>           the C update enters the momentum update with rho-tilde. Central face interpolation of the transported velocity.
!>   RK3     Wray stages with the solver's coefficients: viscous stress (variable viscosity), gravity with the still-water
!>           hydrostatic part removed, and the pressure as an explicit predictor (previous stage) plus an increment from a few PCG
!>           iterations preconditioned by the constant-coefficient fast solver (vof_pressure).
!>
!>  Supported: DNS/LES without wall model, periodic or wall boundaries, ghost-cell IBM (velocity condition), wave inlet and relaxation
!>  zones; not scalars, particles, UAV or the wall models (checked at input).
Module vof_twofluid

  Use iso_fortran_env, Only : Int32, Int64
  Use global
  Use mpi
  Use decomp, Only : x_periodic_partner, z_periodic_partner
  Use boundary_conditions, Only : apply_boundary_conditions, outflow_relax_on, update_ghost_interior_planes, &
                                  update_ghost_interior_planes_x, apply_periodic_bc_x, apply_periodic_bc_z
  Use sgs_models, Only : compute_sgs_model
  Use wallmodel, Only : compute_wall_model
  Use waves, Only : wave_eta, wave_profile_at
  Use monitor, Only : compute_cfl, write_force_csv
  Use vof_plic
  Use vof_normals
  Use vof_advect
  Use vof_state
  Use vof_pressure
  Use vof_curv, Only : vof_curvature
  Use mom_recon
  Use halo_pad, Only : pad_field
  Use ibm, Only : apply_ghost_cell_ibm
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
  Real(Int64), Allocatable, Dimension(:)     :: rhos, rsf, cw0, cw1
  Integer(Int32) :: vf_proj_its_last = 0, vf_proj_its_sum = 0
  ! IBM load output: whether this step is a sampling step
  Logical :: ibm_sampling_now = .False.
  Real(Int64) :: vf_bnd_loc = 0d0, vf_relax_loc = 0d0   ! rank-local liquid volume through the x boundaries / added by relaxation this step

Contains

  Subroutine vof_flow_init

    Integer(Int32) :: i, j, k
    Real(Int64) :: kt(3), om

    If ( ibm_wall_model_flag /= 0 .Or. sediment_flag >= 1 .Or. boussinesq_flag >= 1 .Or. particles_active >= 1 &
         .Or. uav_active >= 1 .Or. flat_wall_model_flag /= 0 .Or. rotation_active >= 1 ) Then
       If ( myid == 0 ) Write(*,'(A)') ' ERROR: vof_flow=1 does not yet support the IBM wall model, sediment, ' // &
            'Boussinesq, particles, UAV, flat wall models or rotation'
       If ( myid == 0 ) Write(*,'(A,7I3)') ' ibm_wall_model, sediment, boussinesq, particles, uav, flat_wall, rotation: ', &
            ibm_wall_model_flag, sediment_flag, boussinesq_flag, particles_active, uav_active, flat_wall_model_flag, rotation_active
       Call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    End If
    If ( dPdx /= 0d0 .Or. dPdz /= 0d0 .Or. flow_forcing_mode /= 0 ) Then
       If ( myid == 0 ) Write(*,'(A)') ' ERROR: vof_flow=1 does not apply dPdx, dPdz or mass-flux forcing (set them to zero)'
       Call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    End If

    Allocate( qu(nx,nyg,nzg), qv(nxg,ny,nzg), qw(nxg,nyg,nz), uqu(nx,nyg,nzg), uqv(nxg,ny,nzg), uqw(nxg,nyg,nz) )
    Allocate( Md(nxg,nyg,nzg), rt(nxg,nyg,nzg), Dd(nxg,nyg,nzg) )
    Allocate( Ut(nx,nyg,nzg), Vt(nxg,ny,nzg), Wt(nxg,nyg,nz), Utmp(nx,nyg,nzg), Vtmp(nxg,ny,nzg), Wtmp(nxg,nyg,nz) )
    Allocate( Hu1(nx,nyg,nzg), Hv1(nxg,ny,nzg), Hw1(nxg,nyg,nz), Hu2(nx,nyg,nzg), Hv2(nxg,ny,nzg), Hw2(nxg,nyg,nz) )
    Allocate( Fu_(nx,nyg,nzg), Fv_(nxg,ny,nzg), Fw_(nxg,nyg,nz), Gu_(nx,nyg,nzg), Gv_(nxg,ny,nzg), Gw_(nxg,nyg,nz) )
    Allocate( ppre(nxg,nyg,nzg), dphi(nxg,nyg,nzg), fdiv(nxg,nyg,nzg), mu_c(nxg,nyg,nzg), mr_c(nxg,nyg,nzg), kap_c(nxg,nyg,nzg) )
    kap_c = 0d0
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
    Else If ( restart == 0 .And. vof_shear == 1 ) Then
       V = 0d0;  W = 0d0
       Do j = 1, nyg
          U(:,j,:) = Merge(vof_u0, 0d0, yg(j) > vof_level)
       End Do
    Else If ( restart == 0 .And. vof_u0 /= 0d0 ) Then
       U = vof_u0;  V = 0d0;  W = 0d0
    End If

    ! pseudo-time RK3 of the momentum needs less than half of a control volume's mass to leave in one sweep, which fails for
    ! density ratios above a few thousand (drop in a uniform stream, 1e4 .. 1e6): forward Euler with small sub-steps there
    If ( vof_rk_nth <= 0 ) vof_rk_nth = 1
    If ( vof_rho_l/vof_rho_g >= 2d3 .And. vof_rk_mom == 1 ) Then
       vof_rk_mom = 0;  vof_co_sub = Min(vof_co_sub, 3d-2)
    End If

    ! the viscous time-step limit sees the largest kinematic viscosity of the two fluids
    nu = Max(vof_nu_l, vof_nu_g)

    If ( restart == 1 ) Then
       ppre = P
       Call vp_halo(ppre, .True.)
    Else
       Call vof_initial_pressure
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
    Do k = 1, nzg
       Do j = 1, ny
          Do i = 1, nxg
             gz(i,j,k) = -vof_grav*( vp_rfv(i,j,k) - rsf(j) )/vp_rfv(i,j,k)
          End Do
       End Do
    End Do
    If ( vof_hsplit /= 0 ) gz = 0d0
    If ( vof_sigma > 0d0 ) Then
       Call vof_curvature(kap_c)
       Call add_surface_tension(zu, gz, zw)
    End If
    Call vp_div(zu, gz, zw, fdiv)
    Call vp_pcg(fdiv, 400, 1d-12, ppre)
    Call vp_halo(ppre, .True.)
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
       Call make_transport_velocity
    Else If ( .Not. keep_transport ) Then
       Call make_transport_velocity
    End If
    Call vof_fill_pad(Cv, nxg, nyg, nzg)
    Call vp_set_density(Cv)
    qu = vp_rau*U;  qv = vp_rav*V;  qw = vp_raw*W
    vof_cc = Merge(1d0, 0d0, Cv(1:nxg,1:nyg,1:nzg) > 0.5d0)
    rt = vof_rho_g + (vof_rho_l - vof_rho_g)*vof_cc
    co = 0d0;  cl = 0d0

    If ( Mod(vof_nadv, 2) == 0 ) Then
       order = (/ 1, 2, 3 /)
    Else
       order = (/ 3, 2, 1 /)
    End If

    Do isw = 1, 3
       d = order(isw)
       Call vof_fill_pad(Cv, nxg, nyg, nzg)
       Call vof_reconstruct(Cv, nxg, nyg, nzg, vof_normal_scheme)
       Call vp_set_density(Cv)
       uqu = qu/vp_rau;  uqv = qv/vp_rav;  uqw = qw/vp_raw
       rou = vp_rau;  rov = vp_rav;  row = vp_raw
       If ( d == 1 ) Then
          Call vof_sweep(1, Cv, nxg, nyg, nzg, Ut, hx, hy, hz, tau, co, cl)
       Else If ( d == 2 ) Then
          Call vof_sweep(2, Cv, nxg, nyg, nzg, Vt, hx, hy, hz, tau, co, cl)
       Else
          Call vof_sweep(3, Cv, nxg, nyg, nzg, Wt, hx, hy, hz, tau, co, cl)
       End If
       If ( d == 1 .And. x_bc_type == 1 ) Call tally_x_boundary_flux
       Call sweep_mass_flux(d, tau)
       Call vp_halo(Md, .False.)
       Call vp_halo(Dd, .False.)
       Md(:,1,:) = 0d0;  Md(:,nyg,:) = 0d0
       Call momentum_update(d)
       Call face_halo(qu, qv, qw)
    End Do
    vof_nadv = vof_nadv + 1
    vof_co_max = Max(vof_co_max, co)
    vof_clip_total = vof_clip_total + cl

    Call vof_fill_pad(Cv, nxg, nyg, nzg)
    Call vp_set_density(Cv)
    Do k = 2, nzg-1
       Do j = 2, nyg-1
          Do i = 2, nx-1
             U(i,j,k) = qu(i,j,k)/vp_rau(i,j,k)
          End Do
       End Do
    End Do
    Do k = 2, nzg-1
       Do j = 2, ny-1
          Do i = 2, nxg-1
             V(i,j,k) = qv(i,j,k)/vp_rav(i,j,k)
          End Do
       End Do
    End Do
    Do k = 2, nz-1
       Do j = 2, nyg-1
          Do i = 2, nxg-1
             W(i,j,k) = qw(i,j,k)/vp_raw(i,j,k)
          End Do
       End Do
    End Do

  End Subroutine vof_advect_half


  !> Transporting velocity: the current velocity made exactly divergence-free by one constant-coefficient fast solve, so the
  !  volume of liquid is conserved to round-off whatever residual the variable-density PCG leaves
  Subroutine make_transport_velocity

    If ( vp_masked ) Then
       Call make_transport_velocity_masked
       Return
    End If
    Utmp = U;  Vtmp = V;  Wtmp = W
    Call compute_pseudo_pressure_rhs
    Call solve_poisson_equation(skip_p_save=.True.)
    Call project_velocity
    Ut = U;  Vt = V;  Wt = W
    Call face_halo(Ut, Vt, Wt)
    U = Utmp;  V = Vtmp;  W = Wtmp

  End Subroutine make_transport_velocity


  !> Transporting velocity with an immersed body: the fast solver cannot honour the closed faces, so the constant-coefficient
  !  projection is done by PCG on the masked operator (Neumann at the body), to round-off. Closed faces carry no flux, hence no
  !  liquid or momentum crosses the body, and the divergence-free field is exact in every fluid cell.
  Subroutine make_transport_velocity_masked

    Real(Int64), Allocatable :: sbu(:,:,:), sbv(:,:,:), sbw(:,:,:), gu(:,:,:), gv(:,:,:), gw(:,:,:), ph(:,:,:)
    Real(Int64) :: beta0_save, umax

    Allocate( sbu(nx,nyg,nzg), sbv(nxg,ny,nzg), sbw(nxg,nyg,nz), gu(nx,nyg,nzg), gv(nxg,ny,nzg), gw(nxg,nyg,nz), ph(nxg,nyg,nzg) )
    sbu = vp_bu;  sbv = vp_bv;  sbw = vp_bw;  beta0_save = vp_beta0
    vp_bu = vp_mu;  vp_bv = vp_mv;  vp_bw = vp_mw;  vp_beta0 = 1d0;  vp_use_layered = .False.

    Utmp = U;  Vtmp = V;  Wtmp = W
    U = U*vp_mu;  V = V*vp_mv;  W = W*vp_mw
    Call vp_div(U, V, W, fdiv)
    ! the iteration stops at the round-off level of the velocity (a relative tolerance alone is unreachable once the field is
    ! already nearly divergence free, and PCG then amplifies noise)
    umax = Max( MaxVal(Abs(U)), MaxVal(Abs(V)), MaxVal(Abs(W)) )
    Call MPI_Allreduce(MPI_IN_PLACE, umax, 1, MPI_real8, MPI_MAX, MPI_COMM_WORLD, ierr)
    Call vp_pcg(fdiv, 800, 1d-13, ph, rfloor=1d-14*umax*Sqrt(vp_wsum)/Min(dx, dymin, dzmin))
    Call vp_halo(ph, .True.)
    Call vp_grad(ph, gu, gv, gw)
    Call apply_face_gradient(gu, gv, gw, 1d0)
    Ut = U;  Vt = V;  Wt = W
    Call face_halo(Ut, Vt, Wt)
    U = Utmp;  V = Vtmp;  W = Wtmp

    vp_bu = sbu;  vp_bv = sbv;  vp_bw = sbw;  vp_beta0 = beta0_save;  vp_use_layered = .True.
    Deallocate( sbu, sbv, sbw, gu, gv, gw, ph )

  End Subroutine make_transport_velocity_masked


  !> Liquid volume that crossed the inlet (positive in) and outlet (positive out) faces in the x-sweep just done
  Subroutine tally_x_boundary_flux

    Logical :: is_first_x, is_last_x, is_first_z, is_last_z
    Integer(Int32) :: partner, khi

    khi = nzg-1
    If ( z_bc_type == 0 ) Then
       Call z_periodic_partner(is_first_z, is_last_z, partner)
       If ( is_last_z ) khi = nzg-2
    End If
    Call x_periodic_partner(is_first_x, is_last_x, partner)
    If ( is_first_x ) vf_bnd_loc = vf_bnd_loc + Sum( vof_flux(1,2:nyg-1,2:khi) )
    If ( is_last_x )  vf_bnd_loc = vf_bnd_loc - Sum( vof_flux(nx,2:nyg-1,2:khi) )

  End Subroutine tally_x_boundary_flux


  !> Mass crossing each face of the sweep direction d in this sweep, and the volume-flux divergence Dd of the sweep
  Subroutine sweep_mass_flux(d, tau)

    Integer(Int32), Intent(In) :: d
    Real(Int64),    Intent(In) :: tau
    Integer(Int32) :: i, j, k
    Real(Int64) :: drho

    drho = vof_rho_l - vof_rho_g
    Md = 0d0
    Dd = 0d0
    Select Case(d)
    Case(1)
       Do k = 2, nzg-1
          Do j = 2, nyg-1
             Do i = 1, nx
                Md(i,j,k) = vof_rho_g*Ut(i,j,k)*tau*hy(j)*hz(k) + drho*vof_flux(i,j,k)
             End Do
             Do i = 2, nx
                Dd(i,j,k) = tau*hy(j)*hz(k)*( Ut(i,j,k) - Ut(i-1,j,k) )
             End Do
          End Do
       End Do
    Case(2)
       Do k = 2, nzg-1
          Do j = 1, ny
             Do i = 2, nxg-1
                Md(i,j,k) = vof_rho_g*Vt(i,j,k)*tau*dx*hz(k) + drho*vof_flux(i,j,k)
             End Do
          End Do
          Do j = 2, ny
             Do i = 2, nxg-1
                Dd(i,j,k) = tau*dx*hz(k)*( Vt(i,j,k) - Vt(i,j-1,k) )
             End Do
          End Do
       End Do
    Case Default
       Do k = 1, Min(nz, nzg-1)
          Do j = 2, nyg-1
             Do i = 2, nxg-1
                Md(i,j,k) = vof_rho_g*Wt(i,j,k)*tau*dx*hy(j) + drho*vof_flux(i,j,k)
             End Do
          End Do
       End Do
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
  !  dq/dtheta = tendency(q/rho(theta)) is integrated by forward Euler (vof_rk_mom = 0) or SSP-RK3 (1; theta = 0, 1, 1/2)
  !  with the same fluxes at every stage; vof_rk_nth > 1 splits the pseudo-time into equal segments (a stage may step
  !  past the end of the density path, which needs less than half of a control volume's mass to leave in one segment)
  Subroutine momentum_update(d)

    Integer(Int32), Intent(In) :: d
    Integer(Int32) :: stage, nst, seg
    Real(Int64) :: a0, a1, th

    Call vof_fill_pad(Cv, nxg, nyg, nzg)
    Call vp_set_density(Cv)
    rmu = Min(rou, vp_rau);  rmv = Min(rov, vp_rav);  rmw = Min(row, vp_raw)
    nst = Merge(3, 1, vof_rk_mom == 1)
    Do seg = 1, vof_rk_nth
       q0u = qu;  q0v = qv;  q0w = qw
       Do stage = 1, nst
          If ( stage > 1 .Or. seg > 1 ) Then
             Call face_halo(qu, qv, qw)
             th = ( Real(seg - 1, Int64) + Merge(0d0, Merge(1d0, 0.5d0, stage == 2), stage == 1) )/vof_rk_nth
             uqu = qu/( rou + th*(vp_rau - rou) );  uqv = qv/( rov + th*(vp_rav - rov) );  uqw = qw/( row + th*(vp_raw - row) )
          End If
          Call pad_transported
          Call momentum_tendency(d)
          a0 = Merge(0d0, Merge(0.75d0, 1d0/3d0, stage == 2), stage == 1);  a1 = 1d0 - a0
          qu = a0*q0u + a1*( qu + dqu/vof_rk_nth );  qv = a0*q0v + a1*( qv + dqv/vof_rk_nth )
          qw = a0*q0w + a1*( qw + dqw/vof_rk_nth )
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
    Integer(Int32) :: i, j, k, ec(3), ed(3)
    Real(Int64) :: corr, vcv

    ec = (/ Merge(1,0,c==1), Merge(1,0,c==2), Merge(1,0,c==3) /)
    ed = (/ Merge(1,0,d==1), Merge(1,0,d==2), Merge(1,0,d==3) /)
    Do k = k0-ed(3), k1
       Do j = j0-ed(2), j1
          Do i = i0-ed(1), i1
             Fl(i,j,k) = face_flux(c, d, i, j, k)
          End Do
       End Do
    End Do
    Do k = k0, k1
       Do j = j0, j1
          Do i = i0, i1
             vcv = dx*Merge(yg(j+1) - yg(j), hy(j), c==2)*Merge(zg(k+1) - zg(k), hz(k), c==3)
             corr = uc(i,j,k)*0.5d0*( rt(i,j,k)*Dd(i,j,k) + rt(i+ec(1),j+ec(2),k+ec(3))*Dd(i+ec(1),j+ec(2),k+ec(3)) )
             dq(i,j,k) = ( corr - ( Fl(i,j,k) - Fl(i-ed(1),j-ed(2),k-ed(3)) ) )/vcv
          End Do
       End Do
    End Do

  End Subroutine tend_comp


  !> Momentum flux of component c through the face of its control volume on the high side in direction d (the flux through the
  !  low face of the next control volume): centre-plane mass flux of the two cells the CV spans, blended face value
  Function face_flux(c, d, i, j, k) Result(f)

    Integer(Int32), Intent(In) :: c, d, i, j, k
    Real(Int64) :: f, mfc, wa, wb, vcv, r1, r2, s(-2:3)
    Integer(Int32) :: ii, jj, kk

    ii = i + Merge(1,0,d==1);  jj = j + Merge(1,0,d==2);  kk = k + Merge(1,0,d==3)
    wa = 0.5d0;  wb = 0.5d0
    Select Case(c)
    Case(1)
       mfc = 0.5d0*( Md(i,j,k) + Md(i+1,j,k) )
       s = gatu(i,j,k, d);  r1 = rmu(i,j,k);  r2 = rmu(ii,jj,kk)
       vcv = dx*hy(j)*hz(k)
       If ( d == 2 ) Then;  wa = weight_y_0(j);  wb = weight_y_1(j);  End If
    Case(2)
       mfc = 0.5d0*( Md(i,j,k) + Md(i,j+1,k) )
       s = gatv(i,j,k, d);  r1 = rmv(i,j,k);  r2 = rmv(ii,jj,kk)
       vcv = dx*( yg(j+1) - yg(j) )*hz(k)
       If ( d == 2 ) Then;  wa = cw0(j+1);  wb = cw1(j+1);  End If
    Case Default
       mfc = 0.5d0*( Md(i,j,k) + Md(i,j,k+1) )
       s = gatw(i,j,k, d);  r1 = rmw(i,j,k);  r2 = rmw(ii,jj,kk)
       vcv = dx*hy(j)*( zg(k+1) - zg(k) )
       If ( d == 2 ) Then;  wa = weight_y_0(j);  wb = weight_y_1(j);  End If
    End Select
    f = mfc*fmom( mfc, s, wa, wb, cmcv(mfc, r1, r2, vcv) )

  End Function face_flux


  !> Face value of the transported velocity for the momentum flux from the six-point stencil s(-2:3) around the face:
  !  vof_mom_scheme 0 weighted central (stretched-grid weights wa, wb), others from mom_recon (index-space formulas). cm is the
  !  face mass flux over the volume times the jump of 1/rho between the adjacent control volumes: the Courant number of the
  !  single-fluid scheme is excluded, so only control volumes that are refilled within a step (update weights must stay
  !  positive there) get the face value blended to first-order upwind
  Pure Function fmom(mf, s, wa, wb, cm) Result(r)

    Real(Int64), Intent(In) :: mf, s(-2:3), wa, wb, cm
    Real(Int64) :: r, up, t

    up = Merge(s(0), s(1), mf >= 0d0)
    If ( vof_mom_scheme == 0 ) Then
       r = wa*s(0) + wb*s(1)
    Else
       r = mom_face(s, mf, vof_mom_scheme)
    End If
    t = Min(1d0, Max(0d0, (cm - vof_mom_cm0)/(vof_mom_cm1 - vof_mom_cm0)))
    r = up + (1d0 - t*t*(3d0 - 2d0*t))*(r - up)

  End Function fmom


  Pure Function cmcv(mf, r1, r2, vcv) Result(c)

    Real(Int64), Intent(In) :: mf, r1, r2, vcv
    Real(Int64) :: c

    c = Abs(mf)*Abs(1d0/r1 - 1d0/r2)/vcv

  End Function cmcv


  Pure Function gatu(i, j, k, dir) Result(s)
    Integer(Int32), Intent(In) :: i, j, k, dir
    Real(Int64) :: s(-2:3)
    Integer(Int32) :: m
    Do m = -2, 3
       s(m) = PadU( i + m*Merge(1,0,dir==1), j + m*Merge(1,0,dir==2), k + m*Merge(1,0,dir==3) )
    End Do
  End Function gatu

  Pure Function gatv(i, j, k, dir) Result(s)
    Integer(Int32), Intent(In) :: i, j, k, dir
    Real(Int64) :: s(-2:3)
    Integer(Int32) :: m
    Do m = -2, 3
       s(m) = PadV( i + m*Merge(1,0,dir==1), j + m*Merge(1,0,dir==2), k + m*Merge(1,0,dir==3) )
    End Do
  End Function gatv

  Pure Function gatw(i, j, k, dir) Result(s)
    Integer(Int32), Intent(In) :: i, j, k, dir
    Real(Int64) :: s(-2:3)
    Integer(Int32) :: m
    Do m = -2, 3
       s(m) = PadW( i + m*Merge(1,0,dir==1), j + m*Merge(1,0,dir==2), k + m*Merge(1,0,dir==3) )
    End Do
  End Function gatw


  !> Padded copies of uqu, uqv, uqw: EP planes beyond the ghost plane in x and z (neighbour rank / periodic partner via pad_field),
  !  EP rows beyond the wall ghost row in y by reflection (tangential: odd for no-slip, even for free-slip; normal: odd about the wall face)
  Subroutine pad_transported
    Real(Int64), Allocatable :: tmp(:,:,:)
    Real(Int64) :: s_lo, s_hi
    Integer(Int32) :: m

    s_lo = Merge(-1d0, 1d0, bc_face_ylo == 1);  s_hi = Merge(-1d0, 1d0, bc_face_yhi == 1)
    Allocate( tmp(1-EP:nx+EP,nyg,1-EP:nzg+EP) )
    Call pad_field(uqu, nx, nyg, nzg, .True., .False., EP, tmp)
    PadU(:,1:nyg,:) = tmp
    Deallocate( tmp )
    Allocate( tmp(1-EP:nxg+EP,ny,1-EP:nzg+EP) )
    Call pad_field(uqv, nxg, ny, nzg, .False., .False., EP, tmp)
    PadV(:,1:ny,:) = tmp
    Deallocate( tmp )
    Allocate( tmp(1-EP:nxg+EP,nyg,1-EP:nz+EP) )
    Call pad_field(uqw, nxg, nyg, nz, .False., .True., EP, tmp)
    PadW(:,1:nyg,:) = tmp
    Deallocate( tmp )
    Do m = 1, EP
       PadU(:,1-m,:)     = s_lo*PadU(:,2+m,:);        PadU(:,nyg+m,:) = s_hi*PadU(:,nyg-1-m,:)
       PadW(:,1-m,:)     = s_lo*PadW(:,2+m,:);        PadW(:,nyg+m,:) = s_hi*PadW(:,nyg-1-m,:)
       PadV(:,1-m,:)     = -PadV(:,1+m,:);             PadV(:,ny+m,:)  = -PadV(:,ny-m,:)
    End Do

  End Subroutine pad_transported


  !> Project U, V, W onto the divergence-free space of the current density: solve div(beta grad phi) = div u by PCG (at most
  !  maxit iterations, relative residual tol) and subtract beta grad phi from the interior faces, as project_velocity does
  Subroutine vof_project(tol, maxit)

    Real(Int64),    Intent(In) :: tol
    Integer(Int32), Intent(In) :: maxit

    Real(Int64), Allocatable :: gu(:,:,:), gv(:,:,:), gw(:,:,:)
    Logical :: is_first_x, is_last_x
    Integer(Int32) :: partner_x

    Allocate( gu(nx,nyg,nzg), gv(nxg,ny,nzg), gw(nxg,nyg,nz) )
    Call vof_fill_pad(Cv, nxg, nyg, nzg)
    Call vp_set_density(Cv)
    Call vp_div(U, V, W, fdiv)
    Call vp_pcg(fdiv, maxit, tol, dphi, dt)
    vf_proj_its_last = vp_iters_last
    vf_proj_its_sum = vf_proj_its_sum + vp_iters_last
    Call vp_halo(dphi, .True.)
    Call vp_grad(dphi, gu, gv, gw)
    Call apply_face_gradient(gu, gv, gw, 1d0)
    Deallocate( gu, gv, gw )

  End Subroutine vof_project


  !> U -= s*gu etc. on the faces the solver's projection updates (interior faces, the row-seam face, the outflow face)
  Subroutine apply_face_gradient(gu, gv, gw, s)

    Real(Int64), Intent(In) :: gu(nx,nyg,nzg), gv(nxg,ny,nzg), gw(nxg,nyg,nz), s
    Logical :: is_first_x, is_last_x
    Integer(Int32) :: partner_x

    Call x_periodic_partner(is_first_x, is_last_x, partner_x)
    U(2:nx-1,2:nyg-1,2:nzg-1) = U(2:nx-1,2:nyg-1,2:nzg-1) - s*gu(2:nx-1,2:nyg-1,2:nzg-1)
    If ( .Not. is_last_x .Or. x_bc_type == 1 ) Then
       U(nx,2:nyg-1,2:nzg-1) = U(nx,2:nyg-1,2:nzg-1) - s*gu(nx,2:nyg-1,2:nzg-1)
    End If
    V(2:nxg-1,2:ny-1,2:nzg-1) = V(2:nxg-1,2:ny-1,2:nzg-1) - s*gv(2:nxg-1,2:ny-1,2:nzg-1)
    W(2:nxg-1,2:nyg-1,2:nz-1) = W(2:nxg-1,2:nyg-1,2:nz-1) - s*gw(2:nxg-1,2:nyg-1,2:nz-1)

  End Subroutine apply_face_gradient


  !> Dynamic viscosity (molecular mixture plus the SGS contribution scaled by the local density) at the cell centres
  Subroutine set_viscosity

    Integer(Int32) :: i, j, k
    Real(Int64) :: mul, mug

    mul = vof_rho_l*vof_nu_l
    mug = vof_rho_g*vof_nu_g
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
       Do k = 1, nzg
          Do i = 1, nxg
             mu_c(i,1,k) = mug + (mul - mug)*Cv(i,1,k) + vp_rho(i,1,k)*nu_t(i,2,k);  mr_c(i,1,k) = 1d0/mu_c(i,1,k)
          End Do
       End Do
    End If
    If ( bc_face_yhi == 2 ) Then
       Do k = 1, nzg
          Do i = 1, nxg
             mu_c(i,nyg,k) = mug + (mul - mug)*Cv(i,nyg,k) + vp_rho(i,nyg,k)*nu_t(i,nyg-1,k)
             mr_c(i,nyg,k) = 1d0/mu_c(i,nyg,k)
          End Do
       End Do
    End If

  End Subroutine set_viscosity


  !> Acceleration from the viscous stress, (1/rho_face) div( mu (grad u + grad u^T) ), at the interior faces
  Subroutine viscous_accel(Fu, Fv, Fw)
  ! shear stresses at the edges use the weighted harmonic mean of mu (series layers carry a continuous shear stress), the
  ! arithmetic mean over-weights the dense fluid next to an interface by the viscosity ratio

    Real(Int64), Intent(Out) :: Fu(nx,nyg,nzg), Fv(nxg,ny,nzg), Fw(nxg,nyg,nz)

    Integer(Int32) :: i, j, k
    Real(Int64) :: inv_dx, inv_dx2, mx1, mx2, my1, my2, mz1, mz2, hyj, hzk, ihy, ihz
    Real(Int64) :: dy1, dy2, dyc, dz1, dz2, dzc

    inv_dx = 1d0/dx
    inv_dx2 = inv_dx*inv_dx
    Fu = 0d0;  Fv = 0d0;  Fw = 0d0

    ! ---- u faces
    Do k = 2, nzg-1
       hzk = z(k) - z(k-1);  ihz = 1d0/hzk
       Do j = 2, nyg-1
          hyj = y(j) - y(j-1);  ihy = 1d0/hyj
          Do i = 2, nx-1
             mx1 = mu_c(i,j,k);  mx2 = mu_c(i+1,j,k)
             my1 = 1d0/( 0.5d0*( weight_y_0(j-1)*( mr_c(i,j-1,k) + mr_c(i+1,j-1,k) ) + weight_y_1(j-1)*( mr_c(i,j,k) &
                   + mr_c(i+1,j,k) ) ) )
             my2 = 1d0/( 0.5d0*( weight_y_0(j)*( mr_c(i,j,k) + mr_c(i+1,j,k) ) + weight_y_1(j)*( mr_c(i,j+1,k) &
                   + mr_c(i+1,j+1,k) ) ) )
             mz1 = 1d0/( 0.5d0*( weight_z_0(k-1)*( mr_c(i,j,k-1) + mr_c(i+1,j,k-1) ) + weight_z_1(k-1)*( mr_c(i,j,k) &
                   + mr_c(i+1,j,k) ) ) )
             mz2 = 1d0/( 0.5d0*( weight_z_0(k)*( mr_c(i,j,k) + mr_c(i+1,j,k) ) + weight_z_1(k)*( mr_c(i,j,k+1) &
                   + mr_c(i+1,j,k+1) ) ) )
             Fu(i,j,k) = ( 2d0*inv_dx2*( mx2*( U(i+1,j,k) - U(i,j,k) ) - mx1*( U(i,j,k) - U(i-1,j,k) ) )              &
                  + ihy*( my2*( ( U(i,j+1,k) - U(i,j,k) )/( yg(j+1) - yg(j) ) + ( V(i+1,j,k) - V(i,j,k) )*inv_dx )   &
                        - my1*( ( U(i,j,k) - U(i,j-1,k) )/( yg(j) - yg(j-1) ) + ( V(i+1,j-1,k) - V(i,j-1,k) )*inv_dx ) ) &
                  + ihz*( mz2*( ( U(i,j,k+1) - U(i,j,k) )/( zg(k+1) - zg(k) ) + ( W(i+1,j,k) - W(i,j,k) )*inv_dx )   &
                        - mz1*( ( U(i,j,k) - U(i,j,k-1) )/( zg(k) - zg(k-1) ) + ( W(i+1,j,k-1) - W(i,j,k-1) )*inv_dx ) ) &
                  )/vp_rfu(i,j,k)
          End Do
       End Do
    End Do

    ! ---- v faces
    Do k = 2, nzg-1
       hzk = z(k) - z(k-1);  ihz = 1d0/hzk
       Do j = 2, ny-1
          dy1 = y(j) - y(j-1);  dy2 = y(j+1) - y(j);  dyc = yg(j+1) - yg(j)
          Do i = 2, nxg-1
             mx1 = 1d0/( 0.5d0*( weight_y_0(j)*( mr_c(i-1,j,k) + mr_c(i,j,k) ) + weight_y_1(j)*( mr_c(i-1,j+1,k) &
                   + mr_c(i,j+1,k) ) ) )
             mx2 = 1d0/( 0.5d0*( weight_y_0(j)*( mr_c(i,j,k) + mr_c(i+1,j,k) ) + weight_y_1(j)*( mr_c(i,j+1,k) &
                   + mr_c(i+1,j+1,k) ) ) )
             my1 = mu_c(i,j,k);  my2 = mu_c(i,j+1,k)
             mz1 = 1d0/( weight_z_0(k-1)*( weight_y_0(j)*mr_c(i,j,k-1) + weight_y_1(j)*mr_c(i,j+1,k-1) ) &
                   + weight_z_1(k-1)*( weight_y_0(j)*mr_c(i,j,k) + weight_y_1(j)*mr_c(i,j+1,k) ) )
             mz2 = 1d0/( weight_z_0(k)*( weight_y_0(j)*mr_c(i,j,k) + weight_y_1(j)*mr_c(i,j+1,k) ) &
                   + weight_z_1(k)*( weight_y_0(j)*mr_c(i,j,k+1) + weight_y_1(j)*mr_c(i,j+1,k+1) ) )
             Fv(i,j,k) = ( inv_dx*( mx2*( ( V(i+1,j,k) - V(i,j,k) )*inv_dx + ( U(i,j+1,k) - U(i,j,k) )/dyc )          &
                                  - mx1*( ( V(i,j,k) - V(i-1,j,k) )*inv_dx + ( U(i-1,j+1,k) - U(i-1,j,k) )/dyc ) )      &
                  + 2d0/dyc*( my2*( V(i,j+1,k) - V(i,j,k) )/dy2 - my1*( V(i,j,k) - V(i,j-1,k) )/dy1 )                 &
                  + ihz*( mz2*( ( V(i,j,k+1) - V(i,j,k) )/( zg(k+1) - zg(k) ) + ( W(i,j+1,k) - W(i,j,k) )/dyc )       &
                        - mz1*( ( V(i,j,k) - V(i,j,k-1) )/( zg(k) - zg(k-1) ) + ( W(i,j+1,k-1) - W(i,j,k-1) )/dyc ) ) &
                  )/vp_rfv(i,j,k)
          End Do
       End Do
    End Do

    ! ---- w faces
    Do k = 2, nz-1
       dz1 = z(k) - z(k-1);  dz2 = z(k+1) - z(k);  dzc = zg(k+1) - zg(k)
       Do j = 2, nyg-1
          hyj = y(j) - y(j-1);  ihy = 1d0/hyj
          Do i = 2, nxg-1
             mx1 = 1d0/( 0.5d0*( weight_z_0(k)*( mr_c(i-1,j,k) + mr_c(i,j,k) ) + weight_z_1(k)*( mr_c(i-1,j,k+1) &
                   + mr_c(i,j,k+1) ) ) )
             mx2 = 1d0/( 0.5d0*( weight_z_0(k)*( mr_c(i,j,k) + mr_c(i+1,j,k) ) + weight_z_1(k)*( mr_c(i,j,k+1) &
                   + mr_c(i+1,j,k+1) ) ) )
             my1 = 1d0/( weight_z_0(k)*( weight_y_0(j-1)*mr_c(i,j-1,k) + weight_y_1(j-1)*mr_c(i,j,k) ) &
                   + weight_z_1(k)*( weight_y_0(j-1)*mr_c(i,j-1,k+1) + weight_y_1(j-1)*mr_c(i,j,k+1) ) )
             my2 = 1d0/( weight_z_0(k)*( weight_y_0(j)*mr_c(i,j,k) + weight_y_1(j)*mr_c(i,j+1,k) ) &
                   + weight_z_1(k)*( weight_y_0(j)*mr_c(i,j,k+1) + weight_y_1(j)*mr_c(i,j+1,k+1) ) )
             mz1 = mu_c(i,j,k);  mz2 = mu_c(i,j,k+1)
             Fw(i,j,k) = ( inv_dx*( mx2*( ( W(i+1,j,k) - W(i,j,k) )*inv_dx + ( U(i,j,k+1) - U(i,j,k) )/dzc )          &
                                  - mx1*( ( W(i,j,k) - W(i-1,j,k) )*inv_dx + ( U(i-1,j,k+1) - U(i-1,j,k) )/dzc ) )      &
                  + ihy*( my2*( ( W(i,j+1,k) - W(i,j,k) )/( yg(j+1) - yg(j) ) + ( V(i,j,k+1) - V(i,j,k) )/dzc )        &
                        - my1*( ( W(i,j,k) - W(i,j-1,k) )/( yg(j) - yg(j-1) ) + ( V(i,j-1,k+1) - V(i,j-1,k) )/dzc ) )  &
                  + 2d0/dzc*( mz2*( W(i,j,k+1) - W(i,j,k) )/dz2 - mz1*( W(i,j,k) - W(i,j,k-1) )/dz1 )                 &
                  )/vp_rfw(i,j,k)
          End Do
       End Do
    End Do

  End Subroutine viscous_accel


  !> Balanced-force surface tension acceleration sigma kappa grad(C)/rho on the faces (the face gradient of C, the mean of the
  !  curvatures of the interface cells next to the face, the same face density as the pressure gradient)
  Subroutine add_surface_tension(Fx, Fy, Fz)

    Real(Int64), Intent(InOut) :: Fx(nx,nyg,nzg), Fy(nxg,ny,nzg), Fz(nxg,nyg,nz)
    Integer(Int32) :: i, j, k

    Do k = 2, nzg-1
       Do j = 2, nyg-1
          Do i = 2, nx-1
             Fx(i,j,k) = Fx(i,j,k) + vof_sigma*face_kappa(i,j,k, i+1,j,k)*( Cv(i+1,j,k) - Cv(i,j,k) )/( dx*vp_rfu(i,j,k) )
          End Do
       End Do
    End Do
    Do k = 2, nzg-1
       Do j = 2, ny-1
          Do i = 2, nxg-1
             Fy(i,j,k) = Fy(i,j,k) + vof_sigma*face_kappa(i,j,k, i,j+1,k)*( Cv(i,j+1,k) - Cv(i,j,k) ) &
                         /( ( yg(j+1) - yg(j) )*vp_rfv(i,j,k) )
          End Do
       End Do
    End Do
    Do k = 2, nz-1
       Do j = 2, nyg-1
          Do i = 2, nxg-1
             Fz(i,j,k) = Fz(i,j,k) + vof_sigma*face_kappa(i,j,k, i,j,k+1)*( Cv(i,j,k+1) - Cv(i,j,k) ) &
                         /( ( zg(k+1) - zg(k) )*vp_rfw(i,j,k) )
          End Do
       End Do
    End Do

  End Subroutine add_surface_tension


  Function face_kappa(i1, j1, k1, i2, j2, k2) Result(kf)

    Integer(Int32), Intent(In) :: i1, j1, k1, i2, j2, k2
    Real(Int64) :: kf, a1, a2

    a1 = Merge(1d0, 0d0, Cv(i1,j1,k1) > vof_eps .And. Cv(i1,j1,k1) < 1d0 - vof_eps)
    a2 = Merge(1d0, 0d0, Cv(i2,j2,k2) > vof_eps .And. Cv(i2,j2,k2) < 1d0 - vof_eps)
    kf = 0d0
    If ( a1 + a2 > 0d0 ) kf = ( a1*kap_c(i1,j1,k1) + a2*kap_c(i2,j2,k2) )/( a1 + a2 )

  End Function face_kappa


  !> Wray RK3 stages with frozen density: explicit viscous stress and hydrostatic-reduced gravity, pressure predictor + PCG increment
  Subroutine vof_forces_rk3(to)

    Real(Int64), Intent(In) :: to
    Integer(Int32) :: s, i, j, k
    Real(Int64) :: ass, a1, a2, dta
    Real(Int64), Allocatable :: bu(:,:,:), bv(:,:,:), bw(:,:,:), gu(:,:,:), gv(:,:,:), gw(:,:,:)

    Allocate( bu(nx,nyg,nzg), bv(nxg,ny,nzg), bw(nxg,nyg,nz), gu(nx,nyg,nzg), gv(nxg,ny,nzg), gw(nxg,nyg,nz) )
    Call vof_fill_pad(Cv, nxg, nyg, nzg)
    Call vp_set_density(Cv)
    If ( vof_sigma > 0d0 ) Call vof_curvature(kap_c)

    Do s = 1, 3
       rk_step = s
       Call compute_sgs_model(U, V, W, nu_t)
       Call compute_wall_model(U, V, W, nu_t)
       Call set_viscosity
       Call viscous_accel(Fu_, Fv_, Fw_)
       If ( vof_hsplit == 0 ) Then
          Do k = 2, nzg-1
             Do j = 2, ny-1
                Do i = 2, nxg-1
                   Fv_(i,j,k) = Fv_(i,j,k) - vof_grav*( vp_rfv(i,j,k) - rsf(j) )/vp_rfv(i,j,k)
                End Do
             End Do
          End Do
       Else
          Do k = 2, nzg-1
             Do j = 2, nyg-1
                Do i = 2, nx-1
                   Fu_(i,j,k) = Fu_(i,j,k) + vof_grav*( yg(j) - vof_level )*( vp_rho(i+1,j,k) - vp_rho(i,j,k) )/( dx*vp_rfu(i,j,k) )
                End Do
             End Do
          End Do
          Do k = 2, nzg-1
             Do j = 2, ny-1
                Do i = 2, nxg-1
                   Fv_(i,j,k) = Fv_(i,j,k) + vof_grav*( y(j) - vof_level )*( vp_rho(i,j+1,k) - vp_rho(i,j,k) ) &
                                /( ( yg(j+1) - yg(j) )*vp_rfv(i,j,k) )
                End Do
             End Do
          End Do
          Do k = 2, nz-1
             Do j = 2, nyg-1
                Do i = 2, nxg-1
                   Fw_(i,j,k) = Fw_(i,j,k) + vof_grav*( yg(j) - vof_level )*( vp_rho(i,j,k+1) - vp_rho(i,j,k) ) &
                                /( ( zg(k+1) - zg(k) )*vp_rfw(i,j,k) )
                End Do
             End Do
          End Do
       End If

       If ( vof_sigma > 0d0 ) Call add_surface_tension(Fu_, Fv_, Fw_)

       ass = rk_coef(s,s)
       dta = dt*ass
       ! stage predictor u* = u_o + dt sum_{j<s} a_sj H_j + dt a_ss F_s  (interior faces; seams and walls come from the BCs)
       a1 = rk_coef(s,1);  a2 = rk_coef(s,2)
       Do k = 2, nzg-1
          Do j = 2, nyg-1
             Do i = 2, nx-1
                U(i,j,k) = Uo(i,j,k) + dta*Fu_(i,j,k)
                If ( s >= 2 ) U(i,j,k) = U(i,j,k) + dt*a1*Hu1(i,j,k)
                If ( s == 3 ) U(i,j,k) = U(i,j,k) + dt*a2*Hu2(i,j,k)
             End Do
          End Do
       End Do
       Do k = 2, nzg-1
          Do j = 2, ny-1
             Do i = 2, nxg-1
                V(i,j,k) = Vo(i,j,k) + dta*Fv_(i,j,k)
                If ( s >= 2 ) V(i,j,k) = V(i,j,k) + dt*a1*Hv1(i,j,k)
                If ( s == 3 ) V(i,j,k) = V(i,j,k) + dt*a2*Hv2(i,j,k)
             End Do
          End Do
       End Do
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
       Call vp_grad(ppre, bu, bv, bw)
       Call apply_face_gradient(bu, bv, bw, dta)
       Call face_halo(U, V, W)
       Call vp_div(U, V, W, fdiv)
       fdiv = fdiv/dta
       Call vp_pcg(fdiv, Max(vof_pcg_iters, 1), vof_pcg_tol, dphi, dt*dta)
       Call vp_halo(dphi, .True.)
       Call vp_grad(dphi, gu, gv, gw)
       Call apply_face_gradient(gu, gv, gw, dta)
       ppre = ppre + dphi
       ! H_s = F_s - beta grad(p_s), with beta grad(p_s) = b + g
       If ( s == 1 ) Then
          Hu1 = Fu_ - (bu + gu);  Hv1 = Fv_ - (bv + gv);  Hw1 = Fw_ - (bw + gw)
       Else If ( s == 2 ) Then
          Hu2 = Fu_ - (bu + gu);  Hv2 = Fv_ - (bv + gv);  Hw2 = Fw_ - (bw + gw)
       End If
       Call apply_boundary_conditions(after_projection=.True.)
       If ( ibm_input_mode >= 1 ) Call enforce_ibm
    End Do
    outflow_relax_on = .False.
    Deallocate( bu, bv, bw, gu, gv, gw )

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
    Allocate( uc(nyg), vf(ny) )
    khi = nzg-1
    If ( z_bc_type == 0 ) Then
       Call z_periodic_partner(is_first_z, is_last_z, partner)
       If ( is_last_z ) khi = nzg-2
    End If

    Do i = 2, nxg-1
       Call zone_at(xg(i), strength, gen)
       If ( strength <= 0d0 ) Cycle
       fr = 1d0 - Exp(-wave_relax_rate*strength*dts)
       If ( gen .And. wave_type > 0 ) Then
          top = vof_level + wave_eta(xg(i), t)
       Else
          top = vof_level
       End If
       Do j = 2, nyg-1
          ct = Min(1d0, Max(0d0, (top - y(j-1))/(y(j) - y(j-1))))
          Do k = 2, nzg-1
             vol = Cv(i,j,k)
             Cv(i,j,k) = Cv(i,j,k) + fr*(ct - Cv(i,j,k))
             If ( k <= khi ) vf_relax_loc = vf_relax_loc + (Cv(i,j,k) - vol)*dx*hy(j)*hz(k)
          End Do
       End Do
       If ( gen .And. wave_type > 0 ) Then
          Call wave_profile_at(xg(i), t, eta, uc, vf)
       Else
          vf = 0d0
       End If
       Do k = 2, nzg-1
          Do j = 2, ny-1
             V(i,j,k) = V(i,j,k) + fr*(vf(j) - V(i,j,k))
          End Do
       End Do
       Do k = 2, nz-1
          Do j = 2, nyg-1
             W(i,j,k) = W(i,j,k)*(1d0 - fr)
          End Do
       End Do
    End Do

    Do i = 2, nx-1
       Call zone_at(x(i), strength, gen)
       If ( strength <= 0d0 ) Cycle
       fr = 1d0 - Exp(-wave_relax_rate*strength*dts)
       If ( gen .And. wave_type > 0 ) Then
          Call wave_profile_at(x(i), t, eta, uc, vf)
       Else
          uc = 0d0
       End If
       Do k = 2, nzg-1
          Do j = 2, nyg-1
             U(i,j,k) = U(i,j,k) + fr*(uc(j) - U(i,j,k))
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


  !> One two-fluid time step: dt from the CFL (the Courant limit sees the gas currents too), then the Strang sequence
  Subroutine compute_time_step_vof

    Real(Int64) :: to, cfl_conv, cfl_visc, cfl_accel, dt_new, dt_presnap, dt_in
    Integer(Int32) :: isub, nsub

    ibm_sampling_now = ( ibm_input_mode >= 1 .And. nsampling > 0 )
    If ( ibm_sampling_now ) ibm_sampling_now = ( Mod(istep, nsampling) == 0 )
    If ( ibm_input_mode >= 1 ) Call enforce_ibm
    Call compute_sgs_model(U, V, W, nu_t)
    Call compute_wall_model(U, V, W, nu_t)
    Call compute_cfl(cfl_conv, cfl_visc, cfl_accel)
    cfl_conv_last  = cfl_conv
    cfl_visc_last  = cfl_visc
    cfl_accel_last = cfl_accel
    If ( vof_sigma > 0d0 ) Then
       ! capillary wave limit of the explicit surface tension (Brackbill): dt < sqrt((rho_l + rho_g) h^3/(4 pi sigma))
       cfl_accel = Max(cfl_accel, dt*Sqrt(16d0*Atan(1d0)*vof_sigma/((vof_rho_l + vof_rho_g)*Min(dx, dymin, dzmin)**3)))
       cfl_accel_last = cfl_accel
    End If
    cfl_current    = Max(cfl_conv, cfl_visc, cfl_accel)
    dt_in = dt
    If ( cfl_adaptive == 1 .And. cfl_current > 0d0 ) Then
       dt_new = dt * cfl_target / cfl_current * cfl_safety
       dt = Max(dt_min, Min(dt_max, dt_new))
    End If
    dt_presnap = dt
    If ( nsave < 0 .And. t + dt > tsave_next ) dt = tsave_next - t
    If ( nsteps < 0 .And. t + dt > sim_end_time ) dt = sim_end_time - t
    to = t
    ! the convective outflow relaxation runs once per RK stage inside vof_forces_rk3 (stage shares of dt sum to dt)
    outflow_relax_on = .False.

    If ( vof_frozen == 0 ) Then
       nsub = vof_nsub
       If ( nsub <= 0 ) nsub = Max(1, Ceiling( 0.5d0*cfl_conv*dt/dt_in/vof_co_sub ))
       vof_nsub_last = nsub
       Do isub = 1, nsub
          Call vof_advect_half(0.5d0*dt/nsub, keep_transport=( isub > 1 .And. vof_freeze_ut == 1 ))
       End Do
    End If
    Call apply_boundary_conditions
    Call vof_project(vof_adv_tol, vof_adv_iters)
    Call apply_boundary_conditions(after_projection=.True.)
    If ( ibm_input_mode >= 1 ) Call enforce_ibm

    Uo = U;  Vo = V;  Wo = W
    Call vof_forces_rk3(to)

    If ( vof_frozen == 0 ) Then
       nsub = vof_nsub
       If ( nsub <= 0 ) nsub = Max(1, Ceiling( 0.5d0*cfl_conv*dt/dt_in/vof_co_sub ))
       Do isub = 1, nsub
          Call vof_advect_half(0.5d0*dt/nsub, keep_transport=( isub > 1 .And. vof_freeze_ut == 1 ))
       End Do
    End If
    t = to + dt
    Call vof_relax(dt)
    Call apply_boundary_conditions
    Call vof_project(vof_adv_tol, vof_adv_iters)
    Call apply_boundary_conditions(after_projection=.True.)
    If ( ibm_input_mode >= 1 ) Call enforce_ibm
    Call vof_fill_pad(Cv, nxg, nyg, nzg)
    Call tally_ledger
    ! the whole array, ghost planes included: snapshots copy every rank's full local block, so a stale (zero) seam ghost plane would
    ! overwrite the neighbour's last real plane in the file, and a restart reads that file back as the pressure predictor
    P = ppre
    P(:,1,:) = P(:,2,:);  P(:,nyg,:) = P(:,nyg-1,:)
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
    dt_step = dt
    dt = dt_presnap

  End Subroutine compute_time_step_vof

End Module vof_twofluid
