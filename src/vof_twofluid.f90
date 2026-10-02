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
!>  Supported so far: DNS/LES without wall model, periodic or wall boundaries; no IBM, scalars, particles or UAV (checked at start).
Module vof_twofluid

  Use iso_fortran_env, Only : Int32, Int64
  Use global
  Use mpi
  Use decomp, Only : x_periodic_partner, z_periodic_partner
  Use boundary_conditions, Only : apply_boundary_conditions, update_ghost_interior_planes, update_ghost_interior_planes_x, &
                                  apply_periodic_bc_x, apply_periodic_bc_z
  Use sgs_models, Only : compute_sgs_model
  Use wallmodel, Only : compute_wall_model
  Use waves, Only : wave_eta, wave_profile_at
  Use monitor, Only : compute_cfl
  Use vof_plic
  Use vof_normals
  Use vof_advect
  Use vof_state
  Use vof_pressure
  Use projection, Only : compute_pseudo_pressure_rhs, solve_poisson_equation, project_velocity

  Implicit None

  ! transported momentum density and its face velocity, mass flux of the current sweep, c-tilde density, volume-flux divergence
  Real(Int64), Allocatable, Dimension(:,:,:) :: qu, qv, qw, uqu, uqv, uqw, Md, rt, Dd
  ! transporting velocity (exactly divergence-free) and a scratch copy of the velocity
  Real(Int64), Allocatable, Dimension(:,:,:) :: Ut, Vt, Wt, Utmp, Vtmp, Wtmp
  ! RK3: start-of-step velocity is Uo/Vo/Wo (global); stage force combinations H = F_explicit - beta grad p of stages 1 and 2
  Real(Int64), Allocatable, Dimension(:,:,:) :: Hu1, Hv1, Hw1, Hu2, Hv2, Hw2, Fu_, Fv_, Fw_, Gu_, Gv_, Gw_
  Real(Int64), Allocatable, Dimension(:,:,:) :: ppre, dphi, fdiv, mu_c
  Real(Int64), Allocatable, Dimension(:)     :: rhos, rsf, cw0, cw1
  Integer(Int32) :: vf_proj_its_last = 0, vf_proj_its_sum = 0
  Real(Int64) :: vf_bnd_loc = 0d0, vf_relax_loc = 0d0   ! rank-local liquid volume through the x boundaries / added by relaxation this step

Contains

  Subroutine vof_flow_init

    Integer(Int32) :: i, j, k

    If ( ibm_input_mode >= 1 .Or. sediment_flag >= 1 .Or. boussinesq_flag >= 1 .Or. particles_active >= 1 &
         .Or. uav_active >= 1 .Or. flat_wall_model_flag /= 0 .Or. rotation_active >= 1 ) Then
       If ( myid == 0 ) Write(*,'(A)') ' ERROR: vof_flow=1 does not yet support IBM, sediment, Boussinesq, particles, UAV, ' // &
            'wall models or rotation'
       Call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    End If

    Allocate( qu(nx,nyg,nzg), qv(nxg,ny,nzg), qw(nxg,nyg,nz), uqu(nx,nyg,nzg), uqv(nxg,ny,nzg), uqw(nxg,nyg,nz) )
    Allocate( Md(nxg,nyg,nzg), rt(nxg,nyg,nzg), Dd(nxg,nyg,nzg) )
    Allocate( Ut(nx,nyg,nzg), Vt(nxg,ny,nzg), Wt(nxg,nyg,nz), Utmp(nx,nyg,nzg), Vtmp(nxg,ny,nzg), Wtmp(nxg,nyg,nz) )
    Allocate( Hu1(nx,nyg,nzg), Hv1(nxg,ny,nzg), Hw1(nxg,nyg,nz), Hu2(nx,nyg,nzg), Hv2(nxg,ny,nzg), Hw2(nxg,nyg,nz) )
    Allocate( Fu_(nx,nyg,nzg), Fv_(nxg,ny,nzg), Fw_(nxg,nyg,nz), Gu_(nx,nyg,nzg), Gv_(nxg,ny,nzg), Gw_(nxg,nyg,nz) )
    Allocate( ppre(nxg,nyg,nzg), dphi(nxg,nyg,nzg), fdiv(nxg,nyg,nzg), mu_c(nxg,nyg,nzg) )
    Allocate( rhos(nyg), rsf(ny), cw0(nyg), cw1(nyg) )
    qu = 0d0;  qv = 0d0;  qw = 0d0;  Md = 0d0;  rt = 0d0;  Dd = 0d0
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
    If ( restart == 0 .And. vof_u0 /= 0d0 ) Then
       U = vof_u0;  V = 0d0;  W = 0d0
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


  !> Advance C and the momentum rho*u by tau with the current velocity U, V, W as the (divergence-free) transporting field
  Subroutine vof_advect_half(tau)

    Real(Int64), Intent(In) :: tau

    Integer(Int32) :: isw, d, order(3), i, j, k
    Real(Int64) :: co, cl

    Call make_transport_velocity
    Call vof_fill_pad(Cv, nxg, nyg, nzg)
    Call vp_set_density(Cv)
    qu = vp_rfu*U;  qv = vp_rfv*V;  qw = vp_rfw*W
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
       uqu = qu/vp_rfu;  uqv = qv/vp_rfv;  uqw = qw/vp_rfw
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
             U(i,j,k) = qu(i,j,k)/vp_rfu(i,j,k)
          End Do
       End Do
    End Do
    Do k = 2, nzg-1
       Do j = 2, ny-1
          Do i = 2, nxg-1
             V(i,j,k) = qv(i,j,k)/vp_rfv(i,j,k)
          End Do
       End Do
    End Do
    Do k = 2, nz-1
       Do j = 2, nyg-1
          Do i = 2, nxg-1
             W(i,j,k) = qw(i,j,k)/vp_rfw(i,j,k)
          End Do
       End Do
    End Do

  End Subroutine vof_advect_half


  !> Transporting velocity: the current velocity made exactly divergence-free by one constant-coefficient fast solve, so the
  !  volume of liquid is conserved to round-off whatever residual the variable-density PCG leaves
  Subroutine make_transport_velocity

    Utmp = U;  Vtmp = V;  Wtmp = W
    Call compute_pseudo_pressure_rhs
    Call solve_poisson_equation(skip_p_save=.True.)
    Call project_velocity
    Ut = U;  Vt = V;  Wt = W
    Call face_halo(Ut, Vt, Wt)
    U = Utmp;  V = Vtmp;  W = Wtmp

  End Subroutine make_transport_velocity


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


  !> Conservative update of the three momentum components for the sweep along d: flux of q through the faces of each staggered
  !  control volume (mass flux averaged onto the CV face, central velocity) plus the rho-tilde * volume-flux-divergence correction
  Subroutine momentum_update(d)

    Integer(Int32), Intent(In) :: d
    Integer(Int32) :: i, j, k
    Real(Int64) :: hi, lo, corr, vcv

    ! ---- u faces (i,j,k): i = 2..nx-1
    Do k = 2, nzg-1
       Do j = 2, nyg-1
          Do i = 2, nx-1
             vcv = dx*hy(j)*hz(k)
             Select Case(d)
             Case(1)
                hi = 0.5d0*( Md(i,j,k) + Md(i+1,j,k) )*0.5d0*( uqu(i,j,k) + uqu(i+1,j,k) )
                lo = 0.5d0*( Md(i-1,j,k) + Md(i,j,k) )*0.5d0*( uqu(i-1,j,k) + uqu(i,j,k) )
                corr = uqu(i,j,k)*0.5d0*( rt(i,j,k)*Dd(i,j,k) + rt(i+1,j,k)*Dd(i+1,j,k) )
             Case(2)
                hi = 0.5d0*( Md(i,j,k) + Md(i+1,j,k) )*( weight_y_0(j)*uqu(i,j,k) + weight_y_1(j)*uqu(i,j+1,k) )
                lo = 0.5d0*( Md(i,j-1,k) + Md(i+1,j-1,k) )*( weight_y_0(j-1)*uqu(i,j-1,k) + weight_y_1(j-1)*uqu(i,j,k) )
                corr = uqu(i,j,k)*0.5d0*( rt(i,j,k)*Dd(i,j,k) + rt(i+1,j,k)*Dd(i+1,j,k) )
             Case Default
                hi = 0.5d0*( Md(i,j,k) + Md(i+1,j,k) )*0.5d0*( uqu(i,j,k) + uqu(i,j,k+1) )
                lo = 0.5d0*( Md(i,j,k-1) + Md(i+1,j,k-1) )*0.5d0*( uqu(i,j,k-1) + uqu(i,j,k) )
                corr = uqu(i,j,k)*0.5d0*( rt(i,j,k)*Dd(i,j,k) + rt(i+1,j,k)*Dd(i+1,j,k) )
             End Select
             qu(i,j,k) = qu(i,j,k) + ( corr - (hi - lo) )/vcv
          End Do
       End Do
    End Do

    ! ---- v faces (i,j,k): j = 2..ny-1
    Do k = 2, nzg-1
       Do j = 2, ny-1
          Do i = 2, nxg-1
             vcv = dx*( yg(j+1) - yg(j) )*hz(k)
             Select Case(d)
             Case(1)
                hi = 0.5d0*( Md(i,j,k) + Md(i,j+1,k) )*0.5d0*( uqv(i,j,k) + uqv(i+1,j,k) )
                lo = 0.5d0*( Md(i-1,j,k) + Md(i-1,j+1,k) )*0.5d0*( uqv(i-1,j,k) + uqv(i,j,k) )
                corr = uqv(i,j,k)*0.5d0*( rt(i,j,k)*Dd(i,j,k) + rt(i,j+1,k)*Dd(i,j+1,k) )
             Case(2)
                ! Md(j) is the face flux above cell j; the centre-plane flux of cell m is the mean of its two face fluxes
                hi = 0.5d0*( Md(i,j,k) + Md(i,j+1,k) )*( cw0(j+1)*uqv(i,j,k) + cw1(j+1)*uqv(i,j+1,k) )
                lo = 0.5d0*( Md(i,j-1,k) + Md(i,j,k) )*( cw0(j)*uqv(i,j-1,k) + cw1(j)*uqv(i,j,k) )
                corr = uqv(i,j,k)*0.5d0*( rt(i,j,k)*Dd(i,j,k) + rt(i,j+1,k)*Dd(i,j+1,k) )
             Case Default
                hi = 0.5d0*( Md(i,j,k) + Md(i,j+1,k) )*0.5d0*( uqv(i,j,k) + uqv(i,j,k+1) )
                lo = 0.5d0*( Md(i,j,k-1) + Md(i,j+1,k-1) )*0.5d0*( uqv(i,j,k-1) + uqv(i,j,k) )
                corr = uqv(i,j,k)*0.5d0*( rt(i,j,k)*Dd(i,j,k) + rt(i,j+1,k)*Dd(i,j+1,k) )
             End Select
             qv(i,j,k) = qv(i,j,k) + ( corr - (hi - lo) )/vcv
          End Do
       End Do
    End Do

    ! ---- w faces (i,j,k): k = 2..nz-1
    Do k = 2, nz-1
       Do j = 2, nyg-1
          Do i = 2, nxg-1
             vcv = dx*hy(j)*( zg(k+1) - zg(k) )
             Select Case(d)
             Case(1)
                hi = 0.5d0*( Md(i,j,k) + Md(i,j,k+1) )*0.5d0*( uqw(i,j,k) + uqw(i+1,j,k) )
                lo = 0.5d0*( Md(i-1,j,k) + Md(i-1,j,k+1) )*0.5d0*( uqw(i-1,j,k) + uqw(i,j,k) )
                corr = uqw(i,j,k)*0.5d0*( rt(i,j,k)*Dd(i,j,k) + rt(i,j,k+1)*Dd(i,j,k+1) )
             Case(2)
                hi = 0.5d0*( Md(i,j,k) + Md(i,j,k+1) )*( weight_y_0(j)*uqw(i,j,k) + weight_y_1(j)*uqw(i,j+1,k) )
                lo = 0.5d0*( Md(i,j-1,k) + Md(i,j-1,k+1) )*( weight_y_0(j-1)*uqw(i,j-1,k) + weight_y_1(j-1)*uqw(i,j,k) )
                corr = uqw(i,j,k)*0.5d0*( rt(i,j,k)*Dd(i,j,k) + rt(i,j,k+1)*Dd(i,j,k+1) )
             Case Default
                hi = 0.5d0*( Md(i,j,k) + Md(i,j,k+1) )*0.5d0*( uqw(i,j,k) + uqw(i,j,k+1) )
                lo = 0.5d0*( Md(i,j,k-1) + Md(i,j,k) )*0.5d0*( uqw(i,j,k-1) + uqw(i,j,k) )
                corr = uqw(i,j,k)*0.5d0*( rt(i,j,k)*Dd(i,j,k) + rt(i,j,k+1)*Dd(i,j,k+1) )
             End Select
             qw(i,j,k) = qw(i,j,k) + ( corr - (hi - lo) )/vcv
          End Do
       End Do
    End Do

  End Subroutine momentum_update


  !> Project U, V, W onto the divergence-free space of the current density: solve div(beta grad phi) = div u by PCG (at most
  !  maxit iterations, relative residual tol) and subtract beta grad phi from the interior faces, as project_velocity does
  Subroutine vof_project(tol, maxit)

    Real(Int64),    Intent(In) :: tol
    Integer(Int32), Intent(In) :: maxit

    Real(Int64), Allocatable :: gu(:,:,:), gv(:,:,:), gw(:,:,:)
    Logical :: is_first_x, is_last_x
    Integer(Int32) :: partner_x

    Allocate( gu(nx,nyg,nzg), gv(nxg,ny,nzg), gw(nxg,nyg,nz) )
    Call vp_set_density(Cv)
    Call vp_div(U, V, W, fdiv)
    Call vp_pcg(fdiv, maxit, tol, dphi)
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
          End Do
       End Do
    End Do

  End Subroutine set_viscosity


  !> Acceleration from the viscous stress, (1/rho_face) div( mu (grad u + grad u^T) ), at the interior faces
  Subroutine viscous_accel(Fu, Fv, Fw)

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
             my1 = 0.5d0*( weight_y_0(j-1)*( mu_c(i,j-1,k) + mu_c(i+1,j-1,k) ) + weight_y_1(j-1)*( mu_c(i,j,k) + mu_c(i+1,j,k) ) )
             my2 = 0.5d0*( weight_y_0(j)*( mu_c(i,j,k) + mu_c(i+1,j,k) ) + weight_y_1(j)*( mu_c(i,j+1,k) + mu_c(i+1,j+1,k) ) )
             mz1 = 0.5d0*( weight_z_0(k-1)*( mu_c(i,j,k-1) + mu_c(i+1,j,k-1) ) + weight_z_1(k-1)*( mu_c(i,j,k) + mu_c(i+1,j,k) ) )
             mz2 = 0.5d0*( weight_z_0(k)*( mu_c(i,j,k) + mu_c(i+1,j,k) ) + weight_z_1(k)*( mu_c(i,j,k+1) + mu_c(i+1,j,k+1) ) )
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
             mx1 = 0.5d0*( weight_y_0(j)*( mu_c(i-1,j,k) + mu_c(i,j,k) ) + weight_y_1(j)*( mu_c(i-1,j+1,k) + mu_c(i,j+1,k) ) )
             mx2 = 0.5d0*( weight_y_0(j)*( mu_c(i,j,k) + mu_c(i+1,j,k) ) + weight_y_1(j)*( mu_c(i,j+1,k) + mu_c(i+1,j+1,k) ) )
             my1 = mu_c(i,j,k);  my2 = mu_c(i,j+1,k)
             mz1 = weight_z_0(k-1)*( weight_y_0(j)*mu_c(i,j,k-1) + weight_y_1(j)*mu_c(i,j+1,k-1) ) &
                 + weight_z_1(k-1)*( weight_y_0(j)*mu_c(i,j,k) + weight_y_1(j)*mu_c(i,j+1,k) )
             mz2 = weight_z_0(k)*( weight_y_0(j)*mu_c(i,j,k) + weight_y_1(j)*mu_c(i,j+1,k) ) &
                 + weight_z_1(k)*( weight_y_0(j)*mu_c(i,j,k+1) + weight_y_1(j)*mu_c(i,j+1,k+1) )
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
             mx1 = 0.5d0*( weight_z_0(k)*( mu_c(i-1,j,k) + mu_c(i,j,k) ) + weight_z_1(k)*( mu_c(i-1,j,k+1) + mu_c(i,j,k+1) ) )
             mx2 = 0.5d0*( weight_z_0(k)*( mu_c(i,j,k) + mu_c(i+1,j,k) ) + weight_z_1(k)*( mu_c(i,j,k+1) + mu_c(i+1,j,k+1) ) )
             my1 = weight_z_0(k)*( weight_y_0(j-1)*mu_c(i,j-1,k) + weight_y_1(j-1)*mu_c(i,j,k) ) &
                 + weight_z_1(k)*( weight_y_0(j-1)*mu_c(i,j-1,k+1) + weight_y_1(j-1)*mu_c(i,j,k+1) )
             my2 = weight_z_0(k)*( weight_y_0(j)*mu_c(i,j,k) + weight_y_1(j)*mu_c(i,j+1,k) ) &
                 + weight_z_1(k)*( weight_y_0(j)*mu_c(i,j,k+1) + weight_y_1(j)*mu_c(i,j+1,k+1) )
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


  !> Wray RK3 stages with frozen density: explicit viscous stress and hydrostatic-reduced gravity, pressure predictor + PCG increment
  Subroutine vof_forces_rk3(to)

    Real(Int64), Intent(In) :: to
    Integer(Int32) :: s, i, j, k
    Real(Int64) :: ass, a1, a2, dta
    Real(Int64), Allocatable :: bu(:,:,:), bv(:,:,:), bw(:,:,:), gu(:,:,:), gv(:,:,:), gw(:,:,:)

    Allocate( bu(nx,nyg,nzg), bv(nxg,ny,nzg), bw(nxg,nyg,nz), gu(nx,nyg,nzg), gv(nxg,ny,nzg), gw(nxg,nyg,nz) )
    Call vp_set_density(Cv)

    Do s = 1, 3
       rk_step = s
       Call compute_sgs_model(U, V, W, nu_t)
       Call compute_wall_model(U, V, W, nu_t)
       Call set_viscosity
       Call viscous_accel(Fu_, Fv_, Fw_)
       Do k = 2, nzg-1
          Do j = 2, ny-1
             Do i = 2, nxg-1
                Fv_(i,j,k) = Fv_(i,j,k) - vof_grav*( vp_rfv(i,j,k) - rsf(j) )/vp_rfv(i,j,k)
             End Do
          End Do
       End Do

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

       Call apply_boundary_conditions
       ! explicit predictor force, then the increment
       Call vp_halo(ppre, .True.)
       Call vp_grad(ppre, bu, bv, bw)
       Call apply_face_gradient(bu, bv, bw, dta)
       Call vp_div(U, V, W, fdiv)
       fdiv = fdiv/dta
       Call vp_pcg(fdiv, Max(vof_pcg_iters, 1), 0d0, dphi)
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
    End Do
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

    Real(Int64) :: to, cfl_conv, cfl_visc, cfl_accel, dt_new, dt_presnap

    Call compute_sgs_model(U, V, W, nu_t)
    Call compute_wall_model(U, V, W, nu_t)
    Call compute_cfl(cfl_conv, cfl_visc, cfl_accel)
    cfl_conv_last  = cfl_conv
    cfl_visc_last  = cfl_visc
    cfl_accel_last = cfl_accel
    cfl_current    = Max(cfl_conv, cfl_visc, cfl_accel)
    If ( cfl_adaptive == 1 .And. cfl_current > 0d0 ) Then
       dt_new = dt * cfl_target / cfl_current * cfl_safety
       dt = Max(dt_min, Min(dt_max, dt_new))
    End If
    dt_presnap = dt
    If ( nsave < 0 .And. t + dt > tsave_next ) dt = tsave_next - t
    If ( nsteps < 0 .And. t + dt > sim_end_time ) dt = sim_end_time - t
    to = t

    Call vof_advect_half(0.5d0*dt)
    Call apply_boundary_conditions
    Call vof_project(vof_adv_tol, vof_adv_iters)
    Call apply_boundary_conditions(after_projection=.True.)

    Uo = U;  Vo = V;  Wo = W
    Call vof_forces_rk3(to)

    Call vof_advect_half(0.5d0*dt)
    t = to + dt
    Call vof_relax(dt)
    Call apply_boundary_conditions
    Call vof_project(vof_adv_tol, vof_adv_iters)
    Call apply_boundary_conditions(after_projection=.True.)
    Call vof_fill_pad(Cv, nxg, nyg, nzg)
    Call tally_ledger
    P(2:nxg-1,2:nyg-1,2:nzg-1) = ppre(2:nxg-1,2:nyg-1,2:nzg-1)
    P(:,1,:) = P(:,2,:);  P(:,nyg,:) = P(:,nyg-1,:)
    Call vof_end_step
    If ( Mod(istep, nmonitor) == 0 ) Call vof_output_monitor
    Call vof_write_gauges
    dt_step = dt
    dt = dt_presnap

  End Subroutine compute_time_step_vof

End Module vof_twofluid
