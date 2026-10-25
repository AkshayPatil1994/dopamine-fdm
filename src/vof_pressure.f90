!> Variable-density pressure operator of the two-fluid solver: face coefficients beta = 1/rho_face, the
!  operator A phi = div(beta grad phi) (the fast Poisson solver's Laplacian for constant beta) and a PCG
!  preconditioned by that solver with the lightest-fluid coefficient (about 17/28/39 iterations to 1e-6 at
!  density ratios 10/100/1000). Cell vectors have shape (nxg,nyg,nzg); bu, bv, bw sit on the U, V, W faces.
!  The periodic duplicate cell and the ghost cells have zero weight in the inner products, so the CG sees only
!  the distinct unknowns.
Module vof_pressure

  Use iso_fortran_env, Only : Int32, Int64
  Use global
  Use mpi
  Use profiler
  Use decomp, Only : x_periodic_partner, z_periodic_partner
  Use boundary_conditions, Only : update_ghost_interior_planes_x, update_ghost_interior_planes, &
                                  apply_periodic_bc_x, apply_periodic_bc_z
  Use scalar_transport, Only : finish_scalar_halos
  Use projection, Only : compute_projection_step, solve_poisson_equation
  Use vof_plic
  Use vof_advect, Only : vof_reconstruct, vof_mx, vof_my, vof_mz, vof_al
  Use halo_pad, Only : pad_field
  Use vof_mg, Only : mg_init, mg_set_coef, mg_precond

  Implicit None

  Real(Int64), Allocatable, Dimension(:,:,:) :: vp_rho, vp_bu, vp_bv, vp_bw, vp_w, vp_rfu, vp_rfv, vp_rfw, vp_rau, vp_rav, vp_raw
  Real(Int64), Allocatable, Dimension(:,:,:) :: vp_riu, vp_riv, vp_riw
  Real(Int64), Allocatable, Dimension(:) :: vp_hy
  Real(Int64), Allocatable, Dimension(:,:,:) :: vp_r, vp_z, vp_d, vp_ap
  Real(Int64) :: vp_beta0 = 1d0, vp_wsum = 1d0
  Integer(Int32) :: vp_iters_last = 0
  Integer(Int64) :: vp_its_total = 0   ! PCG iterations since the start of the run
  Real(Int64) :: vp_res_last = 0d0
  ! Ghost-cell IBM (ibm_input_mode>=1): faces the IBM treats as solid (face-averaged signed distance < 0, as
  ! in apply_ghost_cell_ibm) carry no flux in the pressure operator, the transporting velocity and the
  ! interface transport. vp_mu/mv/mw are the 1/0 face masks (shaped like U, V, W), vp_act the 1/0 cell mask
  ! (0 = every face closed: solid cell), padded by one plane like C.
  Logical :: vp_masked = .False.
  ! masked single-phase path on the GPU: the PCG vectors and coefficients are device-resident
  Logical :: vp_dev = .False.
  Real(Int64), Allocatable, Dimension(:,:,:) :: vp_mu, vp_mv, vp_mw, vp_act
  ! vp_bu/bv/bw changed since the multigrid conductances were built (rebuilt by the next preconditioner call)
  Logical :: vp_mg_stale = .False.

Contains

  !> Whole-array copy as a parallel loop: the compiler runs a multi-array `a = b` kernels region in scalar
  !  mode (one GPU thread)
  Subroutine dev_copy(dst, src)

    Real(Int64), Intent(InOut) :: dst(:,:,:)
    Real(Int64), Intent(In) :: src(:,:,:)
    Integer(Int32) :: i, j, k

    !$acc parallel loop collapse(3) present(dst,src)
    Do k = 1, Size(dst,3)
       Do j = 1, Size(dst,2)
          Do i = 1, Size(dst,1)
             dst(i,j,k) = src(i,j,k)
          End Do
       End Do
    End Do

  End Subroutine dev_copy


  Subroutine vp_init

    Integer(Int32) :: i, j, k, ihi, jhi, khi
    Logical :: is_first, is_last
    Integer(Int32) :: partner
    Real(Int64) :: wl

    Allocate( vp_rho(nxg,nyg,nzg), vp_bu(nx,nyg,nzg), vp_bv(nxg,ny,nzg), vp_bw(nxg,nyg,nz), vp_w(nxg,nyg,nzg) )
    Allocate( vp_r(nxg,nyg,nzg), vp_z(nxg,nyg,nzg), vp_d(nxg,nyg,nzg), vp_ap(nxg,nyg,nzg) )
    Allocate( vp_rfu(nx,nyg,nzg), vp_rfv(nxg,ny,nzg), vp_rfw(nxg,nyg,nz), vp_hy(nyg) )
    Allocate( vp_rau(nx,nyg,nzg), vp_rav(nxg,ny,nzg), vp_raw(nxg,nyg,nz) )
    Allocate( vp_riu(nx,nyg,nzg), vp_riv(nxg,ny,nzg), vp_riw(nxg,nyg,nz) )
    vp_riu = vof_rho_g;  vp_riv = vof_rho_g;  vp_riw = vof_rho_g
    vp_rau = vof_rho_g;  vp_rav = vof_rho_g;  vp_raw = vof_rho_g
    vp_rfu = vof_rho_g;  vp_rfv = vof_rho_g;  vp_rfw = vof_rho_g
    Do j = 2, nyg-1
       vp_hy(j) = y(j) - y(j-1)
    End Do
    vp_hy(1) = vp_hy(2);  vp_hy(nyg) = vp_hy(nyg-1)
    vp_rho = vof_rho_g;  vp_bu = 1d0/vof_rho_g;  vp_bv = 1d0/vof_rho_g;  vp_bw = 1d0/vof_rho_g
    vp_r = 0d0;  vp_z = 0d0;  vp_d = 0d0;  vp_ap = 0d0
    vp_beta0 = 1d0/Min(vof_rho_l, vof_rho_g)

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

    Allocate( vp_act(0:nxg+1,0:nyg+1,0:nzg+1) )
    vp_act = 1d0
    vp_masked = ( ibm_input_mode >= 1 )
    If ( vp_masked ) Call vp_build_ibm_masks

    vp_w = 0d0
    wl = 0d0
    Do k = 2, khi
       Do j = 2, jhi
          Do i = 2, ihi
             vp_w(i,j,k) = dx*(y(j) - y(j-1))*(z(k) - z(k-1))*vp_act(i,j,k)
             wl = wl + vp_w(i,j,k)
          End Do
       End Do
    End Do
    Call MPI_Allreduce(wl, vp_wsum, 1, MPI_real8, MPI_SUM, MPI_COMM_WORLD, ierr)
    If ( pcg_precond >= 1 ) Call mg_init(ihi, jhi, khi)
#ifdef GPU_POISSON
    If ( vof_flow >= 1 ) Then   ! the single-phase staircase path copies its masks in vp_init_masked
       vp_dev = .True.
       !$acc enter data copyin(vp_rho,vp_bu,vp_bv,vp_bw,vp_w,vp_rfu,vp_rfv,vp_rfw,vp_rau,vp_rav,vp_raw,vp_riu,vp_riv,vp_riw,vp_hy,vp_act)
       !$acc enter data create(vp_r,vp_z,vp_d,vp_ap)
       If ( vp_masked ) Then
          !$acc enter data copyin(vp_mu,vp_mv,vp_mw)
       End If
    End If
#endif

  End Subroutine vp_init


  !> U -= s*gu etc. on the faces the solver's projection updates (interior faces, the row-seam face, the
  !  outflow face)
  Subroutine apply_face_gradient(gu, gv, gw, s)

    Real(Int64), Intent(In) :: gu(nx,nyg,nzg), gv(nxg,ny,nzg), gw(nxg,nyg,nz), s
    Logical :: is_first_x, is_last_x
    Integer(Int32) :: partner_x

    Call x_periodic_partner(is_first_x, is_last_x, partner_x)
    !$acc kernels present(U,V,W,gu,gv,gw) if(vp_dev)
    U(2:nx-1,2:nyg-1,2:nzg-1) = U(2:nx-1,2:nyg-1,2:nzg-1) - s*gu(2:nx-1,2:nyg-1,2:nzg-1)
    If ( .Not. is_last_x .Or. x_bc_type == 1 ) Then
       U(nx,2:nyg-1,2:nzg-1) = U(nx,2:nyg-1,2:nzg-1) - s*gu(nx,2:nyg-1,2:nzg-1)
    End If
    V(2:nxg-1,2:ny-1,2:nzg-1) = V(2:nxg-1,2:ny-1,2:nzg-1) - s*gv(2:nxg-1,2:ny-1,2:nzg-1)
    W(2:nxg-1,2:nyg-1,2:nz-1) = W(2:nxg-1,2:nyg-1,2:nz-1) - s*gw(2:nxg-1,2:nyg-1,2:nz-1)
    !$acc end kernels

  End Subroutine apply_face_gradient


  !> Single-phase setup for ibm_method = 1: the face masks (and the device copies); the pressure is the fast
  !  whole-box solve
  Subroutine vp_init_masked

    Allocate( vp_act(0:nxg+1,0:nyg+1,0:nzg+1) )
    vp_masked = .True.
    Call vp_build_ibm_masks
    vp_dev = .True.
    !$acc enter data copyin(vp_mu,vp_mv,vp_mw)

  End Subroutine vp_init_masked


  !> Projection of the single-phase velocity with a staircase body (ibm_method = 1): the fast whole-box solve,
  !  the closed faces zeroed before it (its divergence is that of the masked velocity) and, for ibm_pressure <
  !  2, again after it
  Subroutine vp_project_staircase

    !$acc kernels present(U,V,W,vp_mu,vp_mv,vp_mw) if(vp_dev)
    U = U*vp_mu;  V = V*vp_mv;  W = W*vp_mw
    !$acc end kernels
    Call compute_projection_step
    ! ibm_pressure = 2: the body faces keep the flux of the solve, so the fluid cells stay divergence-free;
    ! they are zeroed again at the start of the next stage, before its predictor is built (zeroing here would
    ! discard that flux, as ibm_pressure = 0 does)
    If ( ibm_pressure >= 2 ) Return
    !$acc kernels present(U,V,W,vp_mu,vp_mv,vp_mw) if(vp_dev)
    U = U*vp_mu;  V = V*vp_mv;  W = W*vp_mw
    !$acc end kernels

  End Subroutine vp_project_staircase


  !> IBM face and cell masks from the signed distance: a face is open only if both cells it joins have phi >=
  !  0 (cell voxelisation, so a wall one cell thick closes its faces); a cell is active if any of its faces is
  !  open.
  Subroutine vp_build_ibm_masks

    Real(Int64), Allocatable :: pp(:,:,:)
    Integer(Int32) :: i, j, k
    Logical :: ou1, ou2, ov1, ov2, ow1, ow2

    If ( .Not. Allocated(phi) ) Then
       If ( myid == 0 ) Write(*,'(A)') ' ERROR: IBM active but the signed distance is not available for the VOF masks'
       Call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    End If
    ! phi padded by two planes in x and z: neighbour ranks / periodic partners supply the planes beyond the
    ! ghost plane
    Allocate( pp(-1:nxg+2,nyg,-1:nzg+2) )
    Call pad_field(phi, nxg, nyg, nzg, .False., .False., 2, pp)

    Allocate( vp_mu(nx,nyg,nzg), vp_mv(nxg,ny,nzg), vp_mw(nxg,nyg,nz) )
    vp_mu = 1d0;  vp_mv = 1d0;  vp_mw = 1d0
    Do k = 1, nzg
       Do j = 1, nyg
          Do i = 1, nx
             If ( Min( pp(i,j,k), pp(i+1,j,k) ) < 0d0 ) vp_mu(i,j,k) = 0d0
          End Do
       End Do
    End Do
    Do k = 1, nzg
       Do j = 2, ny-1
          Do i = 1, nxg
             If ( Min( pp(i,j,k), pp(i,j+1,k) ) < 0d0 ) vp_mv(i,j,k) = 0d0
          End Do
       End Do
    End Do
    Do k = 1, nz
       Do j = 1, nyg
          Do i = 1, nxg
             If ( Min( pp(i,j,k), pp(i,j,k+1) ) < 0d0 ) vp_mw(i,j,k) = 0d0
          End Do
       End Do
    End Do

    ! cell activity on the padded range (0:n+1) so the ghost layers agree with the owning ranks; the wall
    ! rows (j = 1, nyg) stay active
    vp_act = 1d0
    Do k = 0, nzg+1
       Do j = 2, nyg-1
          Do i = 0, nxg+1
             ou1 = ( Min( pp(i-1,j,k), pp(i,j,k) ) >= 0d0 );  ou2 = ( Min( pp(i,j,k), pp(i+1,j,k) ) >= 0d0 )
             ! the wall faces (j = 1 and ny) carry no flux: they do not make a cell active
             ov1 = .False.;  ov2 = .False.
             If ( j-1 >= 2 ) ov1 = ( Min( pp(i,j-1,k), pp(i,j,k) ) >= 0d0 )
             If ( j <= ny-1 ) ov2 = ( Min( pp(i,j,k), pp(i,j+1,k) ) >= 0d0 )
             ow1 = ( Min( pp(i,j,k-1), pp(i,j,k) ) >= 0d0 );  ow2 = ( Min( pp(i,j,k), pp(i,j,k+1) ) >= 0d0 )
             If ( .Not. ( ou1 .Or. ou2 .Or. ov1 .Or. ov2 .Or. ow1 .Or. ow2 ) ) vp_act(i,j,k) = 0d0
          End Do
       End Do
    End Do
    Deallocate( pp )

  End Subroutine vp_build_ibm_masks


  !> Cell densities (ghosts included), face densities (arithmetic in x and z, volume-weighted across the
  !  stretched y cells so that rho_face * u is the mass-consistent momentum of the staggered control volume)
  !  and beta = 1/rho_face, from the padded C Density of the cells and the faces from the liquid fraction.
  !  With sweep_only the face densities and pressure coefficients of the forces and the projections are left
  !  alone: the geometric sweeps need only the arithmetic densities and the planes
  Subroutine vp_set_density(Cp, sweep_only)

    Real(Int64), Intent(In) :: Cp(0:nxg+1,0:nyg+1,0:nzg+1)
    Logical, Intent(In), Optional :: sweep_only
    Integer(Int32) :: i, j, k
    Real(Int64) :: dr, hl, hu

    Call profiler_start(PROF_DENS_BASE)
    !$acc parallel loop collapse(3) present(Cp,vp_rho)
    Do k = 1, nzg
       Do j = 1, nyg
          Do i = 1, nxg
             vp_rho(i,j,k) = vof_rho_g + (vof_rho_l - vof_rho_g)*Cp(i,j,k)
          End Do
       End Do
    End Do
    !$acc parallel loop collapse(3) present(vp_rho,vp_rau)
    Do k = 1, nzg
       Do j = 1, nyg
          Do i = 1, Min(nx, nxg-1)
             vp_rau(i,j,k) = 0.5d0*( vp_rho(i,j,k) + vp_rho(i+1,j,k) )
          End Do
       End Do
    End Do
    ! on a rank with a +x neighbour the last face is the duplicate of the neighbour's first one: its far cell
    ! is the pad layer
    If ( nx == nxg ) Then
       !$acc kernels present(Cp,vp_rho,vp_rau)
       vp_rau(nx,:,:) = 0.5d0*( vp_rho(nx,:,:) + vof_rho_g + (vof_rho_l - vof_rho_g)*Cp(nxg+1,1:nyg,1:nzg) )
       !$acc end kernels
    End If
    !$acc parallel loop collapse(3) present(vp_rho,vp_rav,vp_hy)
    Do k = 1, nzg
       Do j = 1, ny
          Do i = 1, nxg
             vp_rav(i,j,k) = ( vp_rho(i,j,k)*vp_hy(j) + vp_rho(i,j+1,k)*vp_hy(j+1) )/( vp_hy(j) + vp_hy(j+1) )
          End Do
       End Do
    End Do
    !$acc parallel loop collapse(3) present(vp_rho,vp_raw)
    Do k = 1, nz
       Do j = 1, nyg
          Do i = 1, nxg
             vp_raw(i,j,k) = 0.5d0*( vp_rho(i,j,k) + vp_rho(i,j,Min(k+1,nzg)) )
          End Do
       End Do
    End Do
    Call profiler_stop(PROF_DENS_BASE)
    If ( Present(sweep_only) ) Then
       If ( sweep_only ) Then
          If ( vof_geo_density >= 1 ) Call vof_reconstruct(Cp, nxg, nyg, nzg, vof_normal_scheme)
          Return
       End If
    End If
    If ( vof_geo_density >= 1 ) Then
       ! geometric staggered-cell density (Fuster/Arrufat): the liquid in the half cells that make up the
       ! control volume of the face is taken from each cell's reconstructed plane, so that a flat interface
       ! is hydrostatically exact whatever C is
       Call profiler_start(PROF_DENS_RECON)
       Call vof_reconstruct(Cp, nxg, nyg, nzg, vof_normal_scheme)
       Call profiler_stop(PROF_DENS_RECON)
       Call profiler_start(PROF_DENS_GEO)
       dr = vof_rho_l - vof_rho_g
       !$acc parallel loop collapse(3) present(Cp,vp_rfu,vof_mx,vof_my,vof_mz,vof_al)
       Do k = 1, nzg
          Do j = 1, nyg
             Do i = 1, Min(nx, nxg-1)
                vp_rfu(i,j,k) = vof_rho_g + dr*0.5d0*( half_frac(Cp(i,j,k), vof_mx(i,j,k), vof_my(i,j,k), vof_mz(i,j,k), &
                                                                 vof_al(i,j,k), 1, .True.) &
                                                     + half_frac(Cp(i+1,j,k), vof_mx(i+1,j,k), vof_my(i+1,j,k), vof_mz(i+1,j,k), &
                                                                 vof_al(i+1,j,k), 1, .False.) )
             End Do
          End Do
       End Do
       If ( nx == nxg ) Then
          !$acc kernels present(vp_rfu,vp_rau)
          vp_rfu(nx,:,:) = vp_rau(nx,:,:)
          !$acc end kernels
       End If
       !$acc parallel loop collapse(3) private(hl,hu) present(Cp,vp_rfv,vp_hy,vof_mx,vof_my,vof_mz,vof_al)
       Do k = 1, nzg
          Do j = 1, ny
             Do i = 1, nxg
                hl = vp_hy(j);  hu = vp_hy(j+1)
                vp_rfv(i,j,k) = vof_rho_g + dr*( hl*half_frac(Cp(i,j,k), vof_mx(i,j,k), vof_my(i,j,k), vof_mz(i,j,k), &
                                                              vof_al(i,j,k), 2, .True.) &
                                               + hu*half_frac(Cp(i,j+1,k), vof_mx(i,j+1,k), vof_my(i,j+1,k), vof_mz(i,j+1,k), &
                                                              vof_al(i,j+1,k), 2, .False.) )/( hl + hu )
             End Do
          End Do
       End Do
       !$acc parallel loop collapse(3) present(Cp,vp_rfw,vof_mx,vof_my,vof_mz,vof_al)
       Do k = 1, nz
          Do j = 1, nyg
             Do i = 1, nxg
                vp_rfw(i,j,k) = vof_rho_g + dr*0.5d0*( half_frac(Cp(i,j,k), vof_mx(i,j,k), vof_my(i,j,k), vof_mz(i,j,k), &
                                                                 vof_al(i,j,k), 3, .True.) &
                     + half_frac(Cp(i,j,Min(k+1,nzg)), vof_mx(i,j,Min(k+1,nzg)), vof_my(i,j,Min(k+1,nzg)), &
                                 vof_mz(i,j,Min(k+1,nzg)), vof_al(i,j,Min(k+1,nzg)), 3, .False.) )
             End Do
          End Do
       End Do
       Call profiler_stop(PROF_DENS_GEO)
    Else
       Call dev_copy(vp_rfu, vp_rau)
       Call dev_copy(vp_rfv, vp_rav)
       Call dev_copy(vp_rfw, vp_raw)
    End If
    Call profiler_start(PROF_DENS_COEF)
    If ( vof_inertia >= 1 ) Then
       Call dev_copy(vp_riu, vp_rau)
       Call dev_copy(vp_riv, vp_rav)
       Call dev_copy(vp_riw, vp_raw)
    Else
       Call dev_copy(vp_riu, vp_rfu)
       Call dev_copy(vp_riv, vp_rfv)
       Call dev_copy(vp_riw, vp_rfw)
    End If
    !$acc kernels present(vp_riu,vp_riv,vp_riw,vp_bu,vp_bv,vp_bw)
    vp_bu = 1d0/vp_riu;  vp_bv = 1d0/vp_riv;  vp_bw = 1d0/vp_riw
    !$acc end kernels
    If ( vp_masked ) Then
       !$acc kernels present(vp_bu,vp_bv,vp_bw,vp_mu,vp_mv,vp_mw)
       vp_bu = vp_bu*vp_mu;  vp_bv = vp_bv*vp_mv;  vp_bw = vp_bw*vp_mw
       !$acc end kernels
    End If
    vp_mg_stale = ( pcg_precond >= 1 )
    Call profiler_stop(PROF_DENS_COEF)

  End Subroutine vp_set_density


  !> Liquid fraction of the upper (upper=.True.) or lower half of the cell (i,j,k) along direction dir, from
  !  its plane
  Function half_frac(c, mx, my, mz, al, dir, upper) Result(f)
    !$acc routine seq

    Real(Int64),    Intent(In) :: c, mx, my, mz, al
    Integer(Int32), Intent(In) :: dir
    Logical,        Intent(In) :: upper
    Real(Int64) :: f, xlo(3), w(3)

    If ( c <= vof_eps ) Then
       f = 0d0
    Else If ( c >= 1d0 - vof_eps ) Then
       f = 1d0
    Else
       xlo = -0.5d0;  w = 1d0
       w(dir) = 0.5d0
       If ( upper ) xlo(dir) = 0d0
       f = plic_subbox_fraction( mx, my, mz, al, xlo(1), xlo(2), xlo(3), w(1), w(2), w(3) )
    End If

  End Function half_frac


  !> Ghost layers of a cell-centred field: rank seams and periodic wraps, zero gradient at y walls (periodic
  !  wrap for y_bc_type=0), zero gradient at the x inlet and, at the x outlet, mirror (antisym=.False.) or the
  !  Dirichlet-zero reflection of the pressure problem (antisym=.True.)
  Subroutine vp_halo(F, antisym)

    Real(Int64), Intent(InOut) :: F(nxg,nyg,nzg)
    Logical,     Intent(In)    :: antisym

    Logical :: is_first, is_last
    Integer(Int32) :: partner

    If ( vp_dev .And. nprocs == 1 ) Then
       ! one rank: the wraps of the host path below as device copies (z always wraps, as in
       ! finish_scalar_halos)
       !$acc kernels present(F)
       If ( x_bc_type == 0 ) Then
          F(1,:,:) = F(nxg-2,:,:);  F(nxg-1,:,:) = F(2,:,:);  F(nxg,:,:) = F(3,:,:)
       Else
          F(1,:,:) = F(2,:,:)
          If ( antisym ) Then
             F(nxg,:,:) = -F(nxg-1,:,:)
          Else
             F(nxg,:,:) = F(nxg-1,:,:)
          End If
       End If
       F(:,:,1) = F(:,:,nzg-2);  F(:,:,nzg-1) = F(:,:,2);  F(:,:,nzg) = F(:,:,3)
       F(:,1,:) = F(:,2,:);  F(:,nyg,:) = F(:,nyg-1,:)
       !$acc end kernels
       Return
    End If
    If ( vp_dev ) Then
       ! several ranks, all on the device: the order of vp_halo_host (x seam, x walls, z seam, periodic
       ! wraps, y)
       Call update_ghost_interior_planes_x(F, 4)
       If ( x_bc_type == 1 ) Then
          Call x_periodic_partner(is_first, is_last, partner)
          !$acc kernels present(F)
          If ( is_first ) F(1,:,:) = F(2,:,:)
          If ( is_last .And. antisym ) F(nxg,:,:) = -F(nxg-1,:,:)
          If ( is_last .And. .Not. antisym ) F(nxg,:,:) = F(nxg-1,:,:)
          !$acc end kernels
       End If
       Call update_ghost_interior_planes(F, 4)
       If ( x_bc_type == 0 ) Call apply_periodic_bc_x(F, 2)
       Call apply_periodic_bc_z(F, 4)
       !$acc kernels present(F)
       If ( y_bc_type == 0 ) Then
          F(:,1,:) = F(:,nyg-2,:);  F(:,nyg-1,:) = F(:,2,:);  F(:,nyg,:) = F(:,3,:)
       Else
          F(:,1,:) = F(:,2,:);  F(:,nyg,:) = F(:,nyg-1,:)
       End If
       !$acc end kernels
       Return
    End If
    Call vp_halo_host(F, antisym)

  End Subroutine vp_halo


  Subroutine vp_halo_host(F, antisym)

    Real(Int64), Intent(InOut) :: F(nxg,nyg,nzg)
    Logical,     Intent(In)    :: antisym

    Logical :: is_first, is_last
    Integer(Int32) :: partner

    Call update_ghost_interior_planes_x(F, 4)
    If ( x_bc_type == 1 ) Then
       Call x_periodic_partner(is_first, is_last, partner)
       If ( is_first ) F(1,:,:) = F(2,:,:)
       If ( is_last ) Then
          If ( antisym ) Then
             F(nxg,:,:) = -F(nxg-1,:,:)
          Else
             F(nxg,:,:) = F(nxg-1,:,:)
          End If
       End If
    End If
    Call finish_scalar_halos(F)
    If ( y_bc_type == 0 ) Then
       F(:,1,:)       = F(:,nyg-2,:)
       F(:,nyg-1,:)   = F(:,2,:)
       F(:,nyg,:)     = F(:,3,:)
    Else
       F(:,1,:)   = F(:,2,:)
       F(:,nyg,:) = F(:,nyg-1,:)
    End If

  End Subroutine vp_halo_host


  !> out = div(beta grad phi) at the interior cells; phi needs valid ghost layers (vp_halo with
  !  antisym=.True.)
  Subroutine vp_apply_A(phi, out)

    Real(Int64), Intent(In)  :: phi(nxg,nyg,nzg)
    Real(Int64), Intent(Out) :: out(nxg,nyg,nzg)

    Integer(Int32) :: i, j, k
    Real(Int64) :: inv_dx2, hy, hz, ihy_hi, ihy_lo, ihz_hi, ihz_lo, fxh, fxl, fyh, fyl, fzh, fzl

    inv_dx2 = 1d0/(dx*dx)
    ! the interior cells are overwritten below, only the outer shell needs to be zero
    !$acc kernels present(out) if(vp_dev)
    out(1,:,:) = 0d0;  out(nxg,:,:) = 0d0;  out(:,1,:) = 0d0;  out(:,nyg,:) = 0d0;  out(:,:,1) = 0d0;  out(:,:,nzg) = 0d0
    !$acc end kernels
    !$acc parallel loop collapse(3) present(phi,out,vp_bu,vp_bv,vp_bw,y,yg,z,zg) if(vp_dev)
    Do k = 2, nzg-1
       Do j = 2, nyg-1
          Do i = 2, nxg-1
             hz = z(k) - z(k-1)
             ihz_hi = 1d0/(zg(k+1) - zg(k));  ihz_lo = 1d0/(zg(k) - zg(k-1))
             hy = y(j) - y(j-1)
             ihy_hi = 1d0/(yg(j+1) - yg(j));  ihy_lo = 1d0/(yg(j) - yg(j-1))
             fxh = vp_bu(i,j,k)*( phi(i+1,j,k) - phi(i,j,k) )
             fxl = vp_bu(i-1,j,k)*( phi(i,j,k) - phi(i-1,j,k) )
             fyh = vp_bv(i,j,k)*( phi(i,j+1,k) - phi(i,j,k) )*ihy_hi
             fyl = vp_bv(i,j-1,k)*( phi(i,j,k) - phi(i,j-1,k) )*ihy_lo
             fzh = vp_bw(i,j,k)*( phi(i,j,k+1) - phi(i,j,k) )*ihz_hi
             fzl = vp_bw(i,j,k-1)*( phi(i,j,k) - phi(i,j,k-1) )*ihz_lo
             out(i,j,k) = (fxh - fxl)*inv_dx2 + (fyh - fyl)/hy + (fzh - fzl)/hz
          End Do
       End Do
    End Do

  End Subroutine vp_apply_A


  !> Face fields beta grad(phi): gu(1:nx), gv(1:ny), gw(1:nz) (ghost layers of phi must be valid)
  Subroutine vp_grad(phi, gu, gv, gw)

    Real(Int64), Intent(In)  :: phi(nxg,nyg,nzg)
    Real(Int64), Intent(Out) :: gu(nx,nyg,nzg), gv(nxg,ny,nzg), gw(nxg,nyg,nz)

    Integer(Int32) :: i, j, k
    Real(Int64) :: inv_dx

    inv_dx = 1d0/dx
    ! the loops below write gu(1:Min(nx,nxg-1)), gv(:,1:ny), gw(:,:,1:Min(nz,nzg-1)) completely: only a last
    ! face that is not one of them (nx == nxg, nz == nzg) is left to zero
    !$acc kernels present(gu,gv,gw) if(vp_dev)
    If ( nx == nxg ) gu(nx,:,:) = 0d0
    If ( nz == nzg ) gw(:,:,nz) = 0d0
    !$acc end kernels
    !$acc parallel loop collapse(3) present(phi,gu,vp_bu) if(vp_dev)
    Do k = 1, nzg
       Do j = 1, nyg
          Do i = 1, Min(nx, nxg-1)
             gu(i,j,k) = vp_bu(i,j,k)*( phi(i+1,j,k) - phi(i,j,k) )*inv_dx
          End Do
       End Do
    End Do
    !$acc parallel loop collapse(3) present(phi,gv,vp_bv,yg) if(vp_dev)
    Do k = 1, nzg
       Do j = 1, ny
          Do i = 1, nxg
             gv(i,j,k) = vp_bv(i,j,k)*( phi(i,j+1,k) - phi(i,j,k) )/( yg(j+1) - yg(j) )
          End Do
       End Do
    End Do
    !$acc parallel loop collapse(3) present(phi,gw,vp_bw,zg) if(vp_dev)
    Do k = 1, Min(nz, nzg-1)
       Do j = 1, nyg
          Do i = 1, nxg
             gw(i,j,k) = vp_bw(i,j,k)*( phi(i,j,k+1) - phi(i,j,k) )/( zg(k+1) - zg(k) )
          End Do
       End Do
    End Do

  End Subroutine vp_grad


  !> d = div(F) at the interior cells for face fields with the shapes of U, V, W
  Subroutine vp_div(Fu, Fv, Fw, d)

    Real(Int64), Intent(In)  :: Fu(nx,nyg,nzg), Fv(nxg,ny,nzg), Fw(nxg,nyg,nz)
    Real(Int64), Intent(Out) :: d(nxg,nyg,nzg)

    Integer(Int32) :: i, j, k
    Real(Int64) :: inv_dx, inv_hy, inv_hz

    inv_dx = 1d0/dx
    !$acc kernels present(d) if(vp_dev)
    d = 0d0
    !$acc end kernels
    !$acc parallel loop collapse(3) present(Fu,Fv,Fw,d,y,z) if(vp_dev)
    Do k = 2, nzg-1
       Do j = 2, nyg-1
          Do i = 2, nxg-1
             inv_hz = 1d0/(z(k) - z(k-1))
             inv_hy = 1d0/(y(j) - y(j-1))
             d(i,j,k) = ( Fu(i,j,k) - Fu(i-1,j,k) )*inv_dx + ( Fv(i,j,k) - Fv(i,j-1,k) )*inv_hy &
                      + ( Fw(i,j,k) - Fw(i,j,k-1) )*inv_hz
          End Do
       End Do
    End Do
    If ( vp_masked ) Then
       ! blocked (IBM solid) faces carry no flux: every connected fluid region keeps a compatible right-hand
       ! side
       !$acc parallel loop collapse(3) present(Fu,Fv,Fw,d,vp_mu,vp_mv,vp_mw,y,z) if(vp_dev)
       Do k = 2, nzg-1
          Do j = 2, nyg-1
             Do i = 2, nxg-1
                inv_hz = 1d0/(z(k) - z(k-1))
                inv_hy = 1d0/(y(j) - y(j-1))
                d(i,j,k) = ( Fu(i,j,k)*vp_mu(i,j,k) - Fu(i-1,j,k)*vp_mu(i-1,j,k) )*inv_dx &
                         + ( Fv(i,j,k)*vp_mv(i,j,k) - Fv(i,j-1,k)*vp_mv(i,j-1,k) )*inv_hy &
                         + ( Fw(i,j,k)*vp_mw(i,j,k) - Fw(i,j,k-1)*vp_mw(i,j,k-1) )*inv_hz
             End Do
          End Do
       End Do
    End If

  End Subroutine vp_div


  Function vp_dot(a, b) Result(s)

    Real(Int64), Intent(In) :: a(nxg,nyg,nzg), b(nxg,nyg,nzg)
    Real(Int64) :: s, sl
    Integer(Int32) :: i, j, k

    sl = 0d0
    !$acc parallel loop collapse(3) reduction(+:sl) present(vp_w,a,b) if(vp_dev)
    Do k = 1, nzg
       Do j = 1, nyg
          Do i = 1, nxg
             sl = sl + vp_w(i,j,k)*a(i,j,k)*b(i,j,k)
          End Do
       End Do
    End Do
    Call MPI_Allreduce(sl, s, 1, MPI_real8, MPI_SUM, MPI_COMM_WORLD, ierr)

  End Function vp_dot


  !> max|a| over the distinct unknowns (zero-weight ghost and duplicate cells excluded)
  Function vp_maxabs(a) Result(s)

    Real(Int64), Intent(In) :: a(nxg,nyg,nzg)
    Real(Int64) :: s, sl
    Integer(Int32) :: i, j, k

    sl = 0d0
    !$acc parallel loop collapse(3) reduction(max:sl) present(vp_w,a) if(vp_dev)
    Do k = 1, nzg
       Do j = 1, nyg
          Do i = 1, nxg
             If ( vp_w(i,j,k) > 0d0 ) sl = Max(sl, Abs(a(i,j,k)))
          End Do
       End Do
    End Do
    Call MPI_Allreduce(sl, s, 1, MPI_real8, MPI_MAX, MPI_COMM_WORLD, ierr)

  End Function vp_maxabs


  Subroutine vp_remove_mean(a)

    Real(Int64), Intent(InOut) :: a(nxg,nyg,nzg)
    Real(Int64) :: sl, s
    Integer(Int32) :: i, j, k

    ! the pressure is pinned by the Dirichlet outlet when x is not periodic, so the operator is not singular:
    ! a mean removal would perturb the system and stall the PCG at a finite residual
    If ( x_bc_type == 1 ) Return
    sl = 0d0
    !$acc parallel loop collapse(3) reduction(+:sl) present(vp_w,a) if(vp_dev)
    Do k = 1, nzg
       Do j = 1, nyg
          Do i = 1, nxg
             sl = sl + vp_w(i,j,k)*a(i,j,k)
          End Do
       End Do
    End Do
    Call MPI_Allreduce(sl, s, 1, MPI_real8, MPI_SUM, MPI_COMM_WORLD, ierr)
    !$acc kernels present(a,vp_act) if(vp_dev)
    a = a - s/vp_wsum
    ! the solid cells stay exactly zero (the weighted sum is unchanged)
    If ( vp_masked ) a = a*vp_act(1:nxg,1:nyg,1:nzg)
    !$acc end kernels

  End Subroutine vp_remove_mean


  !> z = M^-1 r: constant-coefficient fast solve with beta0, zero weighted mean (keeps the preconditioner
  !  symmetric)
  Subroutine vp_precond(r, z)

    Real(Int64), Intent(In)  :: r(nxg,nyg,nzg)
    Real(Int64), Intent(Out) :: z(nxg,nyg,nzg)

    If ( pcg_precond >= 1 ) Then
       If ( vp_mg_stale ) Then
          Call mg_set_coef(vp_bu, vp_bv, vp_bw)
          vp_mg_stale = .False.
       End If
       Call mg_precond(r, z)
    Else
       z = 0d0
       rhs_p(2:nxg,2:nyg-1,2:nzg) = r(2:nxg,2:nyg-1,2:nzg)/vp_beta0
       Call solve_poisson_equation(skip_p_save=.True.)
       z(2:nxg,2:nyg-1,2:nzg) = rhs_p(2:nxg,2:nyg-1,2:nzg)
    End If
    ! the fast solver knows nothing of the body: keep the correction out of the solid cells (rows of the
    ! operator are empty there)
    !$acc kernels present(z,vp_act) if(vp_dev)
    If ( vp_masked ) z = z*vp_act(1:nxg,1:nyg,1:nzg)
    !$acc end kernels
    Call vp_remove_mean(z)

  End Subroutine vp_precond


  !> Up to nit PCG iterations on A x = f from x = 0 (mean of f removed), stopping at residual < tol * initial
  !  (tol = 0: always nit). With vof_div_tol > 0 and rscale, stops at rscale*max|r| < vof_div_tol (r =
  !  divergence left after the correction). rfloor (optional): absolute round-off residual norm below which
  !  the solve stops. Breakdown (rz or curvature no longer negative) ends the loop. Sets vp_iters_last and the
  !  final relative residual.
  Subroutine vp_pcg(f, nit, tol, x, rscale, rfloor)

    Real(Int64),    Intent(In)  :: f(nxg,nyg,nzg), tol
    Integer(Int32), Intent(In)  :: nit
    Real(Int64),    Intent(Out) :: x(nxg,nyg,nzg)
    Real(Int64),    Intent(In), Optional :: rscale, rfloor

    Integer(Int32) :: it
    Real(Int64) :: rz, rzn, alpha, r0, rn, den
    Logical :: abs_stop

    abs_stop = ( vof_div_tol > 0d0 .And. Present(rscale) )
    !$acc kernels present(x,vp_r,f) if(vp_dev)
    x = 0d0
    vp_r = f
    !$acc end kernels
    vp_iters_last = 0
    vp_res_last = 0d0
    Call vp_remove_mean(vp_r)
    r0 = Sqrt(vp_dot(vp_r, vp_r))
    If ( r0 == 0d0 ) Return
    If ( abs_stop ) Then
       If ( rscale*vp_maxabs(vp_r) < vof_div_tol ) Return
    End If
    Call profiler_start(PROF_PCG_PRECOND)
    Call vp_precond(vp_r, vp_z)
    Call profiler_stop(PROF_PCG_PRECOND)
    !$acc kernels present(vp_d,vp_z) if(vp_dev)
    vp_d = vp_z
    !$acc end kernels
    rz = vp_dot(vp_r, vp_z)
    If ( .Not. ( rz < 0d0 ) ) Return   ! the operator and the preconditioner are negative definite (Laplacians)
    Do it = 1, nit
       Call profiler_start(PROF_PCG_OPER)
       Call vp_halo(vp_d, .True.)
       Call vp_apply_A(vp_d, vp_ap)
       Call profiler_stop(PROF_PCG_OPER)
       den = vp_dot(vp_d, vp_ap)
       If ( .Not. ( den < 0d0 ) ) Exit
       alpha = rz/den
       !$acc kernels present(x,vp_d,vp_r,vp_ap) if(vp_dev)
       x = x + alpha*vp_d
       vp_r = vp_r - alpha*vp_ap
       !$acc end kernels
       Call vp_remove_mean(vp_r)
       vp_iters_last = it
       vp_its_total = vp_its_total + 1
       rn = Sqrt(vp_dot(vp_r, vp_r))
       vp_res_last = rn/r0
       If ( it == nit ) Exit
       If ( Present(rfloor) ) Then
          If ( rn <= rfloor ) Exit
       End If
       If ( abs_stop ) Then
          If ( rscale*vp_maxabs(vp_r) < vof_div_tol ) Exit
       Else If ( rn < tol*r0 ) Then
          Exit
       End If
       Call profiler_start(PROF_PCG_PRECOND)
       Call vp_precond(vp_r, vp_z)
       Call profiler_stop(PROF_PCG_PRECOND)
       rzn = vp_dot(vp_r, vp_z)
       If ( .Not. ( rzn < 0d0 ) ) Exit
       !$acc kernels present(vp_d,vp_z) if(vp_dev)
       vp_d = vp_z + (rzn/rz)*vp_d
       !$acc end kernels
       rz = rzn
    End Do

  End Subroutine vp_pcg


  !> Development check, run from vof_init when vof_selftest=1: PCG iteration counts to 1e-6 and 1e-10 for a
  !  rough zero-mean right-hand side on the current density field, and the symmetry of A (a.Ab - b.Aa)
  Subroutine vp_selftest

    Real(Int64), Allocatable :: f(:,:,:), x(:,:,:), a(:,:,:), b(:,:,:), Aa(:,:,:), Ab(:,:,:)
    Integer(Int32) :: i, j, k
    Real(Int64) :: sab, sba, tol_list(2)
    Integer(Int32) :: it

    Allocate( f(nxg,nyg,nzg), x(nxg,nyg,nzg), a(nxg,nyg,nzg), b(nxg,nyg,nzg), Aa(nxg,nyg,nzg), Ab(nxg,nyg,nzg) )
    f = 0d0
    Do k = 2, nzg-1
       Do j = 2, nyg-1
          Do i = 2, nxg-1
             f(i,j,k) = Sin(6.283d0*xg(i) + 3.14d0*yg(j) + 4d0*zg(k)) + Cos(2d0*xg(i) - 9.42d0*yg(j) + 6d0*zg(k))
             a(i,j,k) = Sin(3d0*xg(i) + 2d0*yg(j)) + Cos(5d0*zg(k))
             b(i,j,k) = Cos(2d0*xg(i) - 3d0*yg(j) + zg(k))
          End Do
       End Do
    End Do
    f = f*vp_w/Max(vp_w, 1d-300)
    a = a*vp_w/Max(vp_w, 1d-300)
    b = b*vp_w/Max(vp_w, 1d-300)
    Call vp_halo(a, .True.);  Call vp_halo(b, .True.)
    Call vp_apply_A(a, Aa);  Call vp_apply_A(b, Ab)
    sab = vp_dot(a, Ab);  sba = vp_dot(b, Aa)
    If ( myid == 0 ) Write(*,'(A,3ES14.5)') '   vp_selftest symmetry a.Ab, b.Aa, rel diff = ', sab, sba, &
         Abs(sab - sba)/Max(Abs(sab), 1d-300)

    tol_list = (/ 1d-6, 1d-10 /)
    Do it = 1, 2
       Call vp_pcg(f, 200, tol_list(it), x)
       If ( myid == 0 ) Write(*,'(A,ES9.1,A,I4,A,ES10.2)') '   vp_selftest PCG tol', tol_list(it), ' iterations', &
            vp_iters_last, ' final residual', vp_res_last
       Call vp_remove_mean(x)
       sab = vp_dot(x, x);  sba = vp_dot(f, x)
       If ( myid == 0 ) Write(*,'(A,2ES20.12)') '   vp_selftest solution norm^2, f.x = ', sab, sba
    End Do
    Deallocate( f, x, a, b, Aa, Ab )

  End Subroutine vp_selftest

End Module vof_pressure
