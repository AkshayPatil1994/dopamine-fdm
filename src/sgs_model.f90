!> SGS eddy-viscosity models for LES (sgs_model: 0=DNS, 1=Vreman)
Module sgs_models

  Use iso_fortran_env, Only : Int32, Int64
  Use global
  Use mpi
  Use wallmodel, Only : kappa_wm, solve_u_tau_wall
  Use decomp, Only : z_periodic_partner, x_periodic_partner
  Use boundary_conditions, Only : apply_periodic_bc_z, update_ghost_interior_planes_x, update_ghost_interior_planes
#ifdef GPU_POISSON
  Use boundary_conditions, Only : gpu_periodic_wrap
#endif

  Implicit None

  ! cells where the IBM pass can differ from the standard one (solid neighbour or filter clamp); built on
  ! first use
  Logical, Allocatable, Dimension(:,:,:) :: sgs_ibm_band

  ! horizontal (x,z) mean of the temperature per y plane, the reference of the AMD buoyancy term
  Real(Int64), Allocatable, Dimension(:) :: Tm_amd

Contains

  !> Dispatch to the SGS model (0 DNS no-op, 1 Vreman, 2 AMD; they share compute_vreman's gradient/IBM/ghost
  !  machinery)
  Subroutine compute_sgs_model(U_, V_, W_, nu_t_)

    Real(Int64), Dimension(nx,  nyg, nzg), Intent(In)    :: U_
    Real(Int64), Dimension(nxg, ny,  nzg), Intent(In)    :: V_
    Real(Int64), Dimension(nxg, nyg, nz ), Intent(In)    :: W_
    Real(Int64), Dimension(nxg, nyg, nzg), Intent(InOut) :: nu_t_

    Select Case (sgs_model)
    Case (0)
       ! DNS: nu_t stays zero (set during initialisation, never modified here)
       Return
    Case (1, 2)
       Call compute_vreman(U_, V_, W_, nu_t_)
    Case Default
       If (myid == 0) Write(*,*) 'WARNING: unknown sgs_model =', sgs_model, ' — using DNS (nu_t=0).'
       Return
    End Select

  End Subroutine compute_sgs_model


  !> Cells whose Pass 2 value can differ from Pass 1: a solid cell among the six neighbours, or body distance
  !  below half the largest cell width (filter clamp)
  Subroutine build_sgs_ibm_band

    Integer(Int32) :: i, j, k
    Real(Int64) :: wmax

    Allocate( sgs_ibm_band(nxg,nyg,nzg) )
    sgs_ibm_band = .True.
    Do k = 2, nzg-1
       Do j = 2, nyg-1
          Do i = 2, nxg-1
             wmax = Max(dx, y(j)-y(j-1), z(k)-z(k-1))
             sgs_ibm_band(i,j,k) = ( 2d0*phi(i,j,k) < wmax ) .Or. &
                  Any( (/ Umask_cc(i-1,j,k), Umask_cc(i+1,j,k), Umask_cc(i,j-1,k), Umask_cc(i,j+1,k), &
                          Umask_cc(i,j,k-1), Umask_cc(i,j,k+1) /) < 0.5d0 )
          End Do
       End Do
    End Do
    !$acc enter data copyin(sgs_ibm_band)

  End Subroutine build_sgs_ibm_band


  !> Vreman (2004) SGS model: nu_t=c_V*sqrt(B_beta/(alpha_ij alpha_ij))
  Subroutine compute_vreman(U_, V_, W_, nu_t_)

    Real(Int64), Dimension(nx,  nyg, nzg), Intent(In)    :: U_
    Real(Int64), Dimension(nxg, ny,  nzg), Intent(In)    :: V_
    Real(Int64), Dimension(nxg, nyg, nz ), Intent(In)    :: W_
    Real(Int64), Dimension(nxg, nyg, nzg), Intent(InOut) :: nu_t_

    ! local scalars
    Integer(Int32) :: i, j, k
    Real   (Int64) :: dx_c, dy_c, dz_c          ! local filter widths (IBM pass only)
    Real   (Int64) :: dx2, dz2                  ! dx^2 (uniform) and dz_c^2 (per k, z may be stretched)
    Real   (Int64) :: inv_dx                    ! 1/dx
    Real   (Int64) :: inv_zg_km, inv_zg_kp      ! per-k z-grid reciprocals
    Real   (Int64) :: inv_dz_k                  ! 1/(z(k)-z(k-1)) — per k
    Real   (Int64) :: inv_yg_jm, inv_yg_jp      ! per-j y-grid reciprocals
    Real   (Int64) :: inv_dy_j                  ! 1/(y(j)-y(j-1)) — per j
    Real   (Int64) :: dy2_c                     ! dy_c^2 — per j (standard path)
    Real   (Int64) :: a11, a12, a13             ! alpha_1j = du_j/dx_1 (at cell ctr)
    Real   (Int64) :: a21, a22, a23             ! alpha_2j = du_j/dx_2
    Real   (Int64) :: a31, a32, a33             ! alpha_3j = du_j/dx_3
    Real   (Int64) :: b11, b12, b13, b22, b23, b33  ! beta_mn
    Real   (Int64) :: B_beta, alpha_sq, c_V, nu_t_loc, amd_num
    Logical        :: use_amd                  ! sgs_model == 2 (AMD) instead of Vreman
    Logical        :: use_buoy                 ! AMD with the Boussinesq buoyancy term
    Real   (Int64) :: buoy_c, th1, th2, th3
    Logical        :: ibm_active               ! Guard Umask_cc accesses
    Logical        :: wm_lo, wm_hi             ! flat-wall EQWM on the low / high y wall
    ! first-cell tangential velocity, friction velocity, wall distance
    Real   (Int64) :: u_c, w_c, u_mag, u_tau, fac_w, y_w

    c_V = 2.5d0 * Cs_vreman**2
    use_amd = sgs_model == 2
    use_buoy = use_amd .And. amd_buoyancy == 1 .And. boussinesq_flag >= 1 .And. vof_flow < 1
    buoy_c = beta_T * grav
    If ( use_buoy ) Call amd_temperature_mean
    ibm_active = ibm_input_mode >= 1 .And. Allocated(Umask_cc)
    wm_lo = vof_flow < 1 .And. flat_wall_model_flag > 0 .And. y_bc_type == 1 .And. bc_face_ylo == 1
    wm_hi = vof_flow < 1 .And. flat_wall_model_flag > 0 .And. y_bc_type == 1 .And. bc_face_yhi == 1

    inv_dx = 1d0 / dx
    dx2    = dx * dx

    ! Pass 1: standard stencil, unconditional so it stays GPU-offloadable without phi; Pass 2 / solid-cell
    ! zeroing discards it where IBM applies
    !$acc parallel loop collapse(2) present(U_,V_,W_,nu_t_,y,yg,z,zg,xg)
    Do k = 2, nzg-1
       Do j = 2, nyg-1

          ! Per-(k,j) constants: y and z spacings (non-uniform y; z if z_bc_type==1 with alpha_grid_z>0)
          dz_c      = z(k) - z(k-1)
          dz2       = dz_c * dz_c
          inv_dz_k  = 1d0 / dz_c
          inv_zg_km = 1d0 / ( zg(k)   - zg(k-1) )
          inv_zg_kp = 1d0 / ( zg(k+1) - zg(k  ) )
          dy_c      = y(j) - y(j-1)
          dy2_c     = dy_c * dy_c
          inv_yg_jm = 1d0 / ( yg(j)   - yg(j-1) )
          inv_yg_jp = 1d0 / ( yg(j+1) - yg(j  ) )
          inv_dy_j  = 1d0 / dy_c

          Do i = 2, nxg-1

             ! ------ velocity gradient tensor (standard stencils) -------
             a11 = ( U_(i,j,k)   - U_(i-1,j,k)   ) * inv_dx
             ! Non-uniform y: use per-j precomputed inverses
             a21 = 0.5d0*( ( U_(i,j,k)   - U_(i,j-1,k)  ) * inv_yg_jm + &
                           ( U_(i,j+1,k) - U_(i,j,  k)  ) * inv_yg_jp )
             a31 = 0.5d0*( ( U_(i,j,k)   - U_(i,j,k-1)  ) * inv_zg_km + &
                           ( U_(i,j,k+1) - U_(i,j,k  )  ) * inv_zg_kp )
             a12 = 0.5d0*( ( V_(i,j,k)   - V_(i-1,j,k)  ) + &
                           ( V_(i+1,j,k) - V_(i,  j,k)  ) ) * inv_dx
             ! Non-uniform y: V on y-faces
             a22 = ( V_(i,j,k)   - V_(i,j-1,k)   ) * inv_dy_j
             a32 = 0.5d0*( ( V_(i,j,k)   - V_(i,j,k-1)  ) * inv_zg_km + &
                           ( V_(i,j,k+1) - V_(i,j,k  )  ) * inv_zg_kp )
             a13 = 0.5d0*( ( W_(i,j,k)   - W_(i-1,j,k)  ) + &
                           ( W_(i+1,j,k) - W_(i,  j,k)  ) ) * inv_dx
             ! Non-uniform y
             a23 = 0.5d0*( ( W_(i,j,k)   - W_(i,j-1,k)  ) * inv_yg_jm + &
                           ( W_(i,j+1,k) - W_(i,j,  k)  ) * inv_yg_jp )
             ! W on z-faces
             a33 = ( W_(i,j,k)   - W_(i,j,k-1)   ) * inv_dz_k

             ! wall-modelled wall cell: the no-slip mirror ghost overstates du/dy several-fold at y+ ~
             ! 50-100, so use the log-law gradient u_tau/(kappa y) (u_tau re-solved from the first-cell
             ! velocity)
             If ( (wm_lo .And. j == 2) .Or. (wm_hi .And. j == nyg-1) ) Then
                y_w   = Merge(yg(2), Ly - yg(nyg-1), j == 2)
                u_c   = 0.5d0*( U_(i-1,j,k) + U_(i,j,k) )
                w_c   = 0.5d0*( W_(i,j,k-1) + W_(i,j,k) )
                u_mag = Sqrt(u_c**2 + w_c**2)
                If ( u_mag > 1d-14 ) Then
                   Call solve_u_tau_wall(u_mag, y_w, Merge(z0_ylo, z0_yhi, j == 2), nu, u_tau)
                   fac_w = Merge(1d0, -1d0, j == 2) * u_tau / (kappa_wm*y_w*u_mag)
                   a21 = fac_w*u_c
                   a23 = fac_w*w_c
                End If
             End If

             ! ------ beta_mn (standard: dx_c=dx, dz_c=z(k)-z(k-1), dy_c=y(j)-y(j-1)) --
             b11 = dx2*(a11*a11) + dy2_c*(a21*a21) + dz2*(a31*a31)
             b22 = dx2*(a12*a12) + dy2_c*(a22*a22) + dz2*(a32*a32)
             b33 = dx2*(a13*a13) + dy2_c*(a23*a23) + dz2*(a33*a33)
             b12 = dx2*(a11*a12) + dy2_c*(a21*a22) + dz2*(a31*a32)
             b13 = dx2*(a11*a13) + dy2_c*(a21*a23) + dz2*(a31*a33)
             b23 = dx2*(a12*a13) + dy2_c*(a22*a23) + dz2*(a32*a33)

             ! ------ B_beta: second invariant of beta -------------------
             B_beta = b11*b22 - b12**2 &
                    + b11*b33 - b13**2 &
                    + b22*b33 - b23**2

             ! ------ alpha_ij alpha_ij = sum of squared gradients -------
             alpha_sq = a11**2 + a12**2 + a13**2 &
                      + a21**2 + a22**2 + a23**2 &
                      + a31**2 + a32**2 + a33**2

             ! ------ eddy viscosity -------------------------------------
             If ( use_amd ) Then
                ! AMD (Rozema et al. 2015): beta_mn is the delta-scaled gradient product, numerator -beta_mn
                ! S_mn
                amd_num = -( b11*a11 + b22*a22 + b33*a33 + b12*(a12 + a21) + b13*(a13 + a31) + b23*(a23 + a32) )
                If ( use_buoy ) Then
                   ! AMD buoyancy production (Abkar & Moin 2017): + beta g (d_k v)(d_k theta'), theta' = T -
                   ! <T>_xz
                   th1 = ( Tscal(i+1,j,k) - Tscal(i-1,j,k) ) / ( xg(i+1) - xg(i-1) )
                   th2 = ( Tscal(i,j+1,k) - Tscal(i,j-1,k) - Tm_amd(j+1) + Tm_amd(j-1) ) / ( yg(j+1) - yg(j-1) )
                   th3 = ( Tscal(i,j,k+1) - Tscal(i,j,k-1) ) / ( zg(k+1) - zg(k-1) )
                   amd_num = amd_num + buoy_c*( dx2*a12*th1 + dy2_c*a22*th2 + dz2*a32*th3 )
                End If
                ! the numerator cancels to round-off in near-irrotational flow: drop values below 1e-12 of the
                ! Smagorinsky-type scale
                If ( alpha_sq > 1d-20 .And. amd_num > 1d-12*(dx2 + dy2_c + dz2)*alpha_sq*Sqrt(alpha_sq) ) Then
                   nu_t_loc = C_amd * amd_num / alpha_sq
                Else
                   nu_t_loc = 0d0
                End If
             Else If ( alpha_sq > 1d-20 .And. B_beta > 0d0 ) Then
                nu_t_loc = c_V * Sqrt( B_beta / alpha_sq )
             Else
                nu_t_loc = 0d0
             End If

             nu_t_(i,j,k) = nu_t_loc

          End Do
       End Do
    End Do
    !$acc end parallel loop

    ! Pass 2: IBM corrections (ibm_active); re-visits fluid cells next to the body
    If ( ibm_active ) Then
       If ( .Not. Allocated(sgs_ibm_band) ) Call build_sgs_ibm_band
       !$acc parallel loop collapse(2) present(U_,V_,W_,nu_t_,y,yg,z,zg,xg,phi,Umask_cc,sgs_ibm_band)
       Do k = 2, nzg-1
          Do j = 2, nyg-1
             dz_c      = z(k) - z(k-1)
             inv_dz_k  = 1d0 / dz_c
             inv_zg_km = 1d0 / ( zg(k)   - zg(k-1) )
             inv_zg_kp = 1d0 / ( zg(k+1) - zg(k  ) )
             dy_c      = y(j) - y(j-1)
             inv_yg_jm = 1d0 / ( yg(j)   - yg(j-1) )
             inv_yg_jp = 1d0 / ( yg(j+1) - yg(j  ) )
             inv_dy_j  = 1d0 / dy_c

             Do i = 2, nxg-1
                If ( phi(i,j,k) < 0d0 ) Cycle   ! solid cell: zeroed in post-processing
                If ( .Not. sgs_ibm_band(i,j,k) ) Cycle

                ! Re-evaluate gradients; then apply one-sided corrections where needed
                a11 = ( U_(i,j,k)   - U_(i-1,j,k)   ) * inv_dx
                a21 = 0.5d0*( (U_(i,j,k)-U_(i,j-1,k))*inv_yg_jm + (U_(i,j+1,k)-U_(i,j,k))*inv_yg_jp )
                a31 = 0.5d0*( (U_(i,j,k)-U_(i,j,k-1))*inv_zg_km + (U_(i,j,k+1)-U_(i,j,k))*inv_zg_kp )
                a12 = 0.5d0*( (V_(i,j,k)-V_(i-1,j,k)) + (V_(i+1,j,k)-V_(i,j,k)) ) * inv_dx
                a22 = ( V_(i,j,k)   - V_(i,j-1,k)   ) * inv_dy_j
                a32 = 0.5d0*( (V_(i,j,k)-V_(i,j,k-1))*inv_zg_km + (V_(i,j,k+1)-V_(i,j,k))*inv_zg_kp )
                a13 = 0.5d0*( (W_(i,j,k)-W_(i-1,j,k)) + (W_(i+1,j,k)-W_(i,j,k)) ) * inv_dx
                a23 = 0.5d0*( (W_(i,j,k)-W_(i,j-1,k))*inv_yg_jm + (W_(i,j+1,k)-W_(i,j,k))*inv_yg_jp )
                a33 = ( W_(i,j,k)   - W_(i,j,k-1)   ) * inv_dz_k

                ! One-sided stencil corrections near IBM solid faces
                If ( Umask_cc(i,j-1,k) < 0.5d0 .And. Umask_cc(i,j+1,k) < 0.5d0 ) Then
                   a21 = 0d0;  a23 = 0d0
                Else If ( Umask_cc(i,j-1,k) < 0.5d0 ) Then
                   a21 = ( U_(i,j+1,k) - U_(i,j,k) ) * inv_yg_jp
                   a23 = ( W_(i,j+1,k) - W_(i,j,k) ) * inv_yg_jp
                Else If ( Umask_cc(i,j+1,k) < 0.5d0 ) Then
                   a21 = ( U_(i,j,k)   - U_(i,j-1,k) ) * inv_yg_jm
                   a23 = ( W_(i,j,k)   - W_(i,j-1,k) ) * inv_yg_jm
                End If

                If ( Umask_cc(i,j,k-1) < 0.5d0 .And. Umask_cc(i,j,k+1) < 0.5d0 ) Then
                   a31 = 0d0;  a32 = 0d0
                Else If ( Umask_cc(i,j,k-1) < 0.5d0 ) Then
                   a31 = ( U_(i,j,k+1) - U_(i,j,k) ) * inv_zg_kp
                   a32 = ( V_(i,j,k+1) - V_(i,j,k) ) * inv_zg_kp
                Else If ( Umask_cc(i,j,k+1) < 0.5d0 ) Then
                   a31 = ( U_(i,j,k)   - U_(i,j,k-1) ) * inv_zg_km
                   a32 = ( V_(i,j,k)   - V_(i,j,k-1) ) * inv_zg_km
                End If

                If ( Umask_cc(i-1,j,k) < 0.5d0 .And. Umask_cc(i+1,j,k) < 0.5d0 ) Then
                   a12 = 0d0;  a13 = 0d0
                Else If ( Umask_cc(i-1,j,k) < 0.5d0 ) Then
                   a12 = ( V_(i+1,j,k) - V_(i,j,k) ) * inv_dx
                   a13 = ( W_(i+1,j,k) - W_(i,j,k) ) * inv_dx
                Else If ( Umask_cc(i+1,j,k) < 0.5d0 ) Then
                   a12 = ( V_(i,j,k)   - V_(i-1,j,k) ) * inv_dx
                   a13 = ( W_(i,j,k)   - W_(i-1,j,k) ) * inv_dx
                End If

                ! IBM filter-width clamping: collapse filter to zero at the surface
                dx_c = Min(dx, 2d0*phi(i,j,k))
                dy_c = Min(y(j)-y(j-1), 2d0*phi(i,j,k))
                dz_c = Min(z(k)-z(k-1), 2d0*phi(i,j,k))

                b11 = dx_c*dx_c*(a11*a11) + dy_c*dy_c*(a21*a21) + dz_c*dz_c*(a31*a31)
                b22 = dx_c*dx_c*(a12*a12) + dy_c*dy_c*(a22*a22) + dz_c*dz_c*(a32*a32)
                b33 = dx_c*dx_c*(a13*a13) + dy_c*dy_c*(a23*a23) + dz_c*dz_c*(a33*a33)
                b12 = dx_c*dx_c*(a11*a12) + dy_c*dy_c*(a21*a22) + dz_c*dz_c*(a31*a32)
                b13 = dx_c*dx_c*(a11*a13) + dy_c*dy_c*(a21*a23) + dz_c*dz_c*(a31*a33)
                b23 = dx_c*dx_c*(a12*a13) + dy_c*dy_c*(a22*a23) + dz_c*dz_c*(a32*a33)

                B_beta = b11*b22 - b12**2 &
                       + b11*b33 - b13**2 &
                       + b22*b33 - b23**2

                alpha_sq = a11**2 + a12**2 + a13**2 &
                         + a21**2 + a22**2 + a23**2 &
                         + a31**2 + a32**2 + a33**2

                If ( use_amd ) Then
                   amd_num = -( b11*a11 + b22*a22 + b33*a33 + b12*(a12 + a21) + b13*(a13 + a31) + b23*(a23 + a32) )
                   If ( use_buoy ) Then
                      th1 = ( Tscal(i+1,j,k) - Tscal(i-1,j,k) ) / ( xg(i+1) - xg(i-1) )
                      th2 = ( Tscal(i,j+1,k) - Tscal(i,j-1,k) - Tm_amd(j+1) + Tm_amd(j-1) ) / ( yg(j+1) - yg(j-1) )
                      th3 = ( Tscal(i,j,k+1) - Tscal(i,j,k-1) ) / ( zg(k+1) - zg(k-1) )
                      amd_num = amd_num + buoy_c*( dx_c*dx_c*a12*th1 + dy_c*dy_c*a22*th2 + dz_c*dz_c*a32*th3 )
                   End If
                   If ( alpha_sq > 1d-20 .And. amd_num > 1d-12*(dx_c*dx_c + dy_c*dy_c + dz_c*dz_c)*alpha_sq*Sqrt(alpha_sq) ) Then
                      nu_t_(i,j,k) = C_amd * amd_num / alpha_sq
                   Else
                      nu_t_(i,j,k) = 0d0
                   End If
                Else If ( alpha_sq > 1d-20 .And. B_beta > 0d0 ) Then
                   nu_t_(i,j,k) = c_V * Sqrt( B_beta / alpha_sq )
                Else
                   nu_t_(i,j,k) = 0d0
                End If

             End Do
          End Do
       End Do
    End If

    Call sgs_fill_ghosts(nu_t_, ibm_active)
#ifdef GPU_POISSON
    ! host copy only needed by host consumers (snapshot write syncs it itself, see input_output.f90)
    If ( rsb_active == 1 ) Then
       !$acc update host(nu_t_)
    End If
#endif

  End Subroutine compute_vreman

  !> Horizontal (x,z) mean of T per y plane over all ranks (Tm_amd), reference of the AMD buoyancy term
  Subroutine amd_temperature_mean

    Integer(Int32) :: i, j, k
    Real   (Int64) :: acc, sums(nyg), cnt(1), cnt_g(1), sums_g(nyg)

    If ( .Not. Allocated(Tm_amd) ) Then
       Allocate( Tm_amd(nyg) )
       !$acc enter data create(Tm_amd)
    End If
    !$acc parallel loop present(Tscal) private(acc) copyout(sums)
    Do j = 1, nyg
       acc = 0d0
       !$acc loop collapse(2) reduction(+:acc)
       Do k = 2, nzg-1
          Do i = 2, nxg-1
             acc = acc + Tscal(i,j,k)
          End Do
       End Do
       sums(j) = acc
    End Do
    !$acc end parallel loop
    cnt(1) = Real( (nxg-2)*(nzg-2), Int64 )
    Call Mpi_allreduce(sums, sums_g, nyg, Mpi_real8, Mpi_sum, MPI_COMM_WORLD, ierr)
    Call Mpi_allreduce(cnt, cnt_g, 1, Mpi_real8, Mpi_sum, MPI_COMM_WORLD, ierr)
    Tm_amd = sums_g / cnt_g(1)
    !$acc update device(Tm_amd)

  End Subroutine amd_temperature_mean


  !> Turbulent scalar diffusivity kappa_t for compute_rhs_scalar_core: nu_t/Pr_t (1/Sc_t) for Vreman and DNS;
  !  AMD scalar closure (Abkar, Bae & Moin 2016) for sgs_model == 2: kappa_t = C_amd max(0, -(delta_k d_k
  !  u_i)(delta_k d_k S) d_i S)/(d_l S d_l S)
  Subroutine compute_scalar_diffusivity(U_, V_, W_, S_, kap_inv)

    Real(Int64), Dimension(nx,  nyg, nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg, ny,  nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg, nyg, nz ), Intent(In) :: W_
    Real(Int64), Dimension(nxg, nyg, nzg), Intent(In) :: S_
    Real(Int64),                           Intent(In) :: kap_inv

    Integer(Int32) :: i, j, k
    Real   (Int64) :: dx2, dy2_c, dz2, inv_dx, inv_yg_jm, inv_yg_jp, inv_dy_j, inv_zg_km, inv_zg_kp, inv_dz_k
    Real   (Int64) :: a11, a12, a13, a21, a22, a23, a31, a32, a33, s1, s2, s3, num, den
    Logical        :: ibm_active

    If ( sgs_model /= 2 ) Then
       !$acc kernels present(kappa_t,nu_t)
       kappa_t = nu_t * kap_inv
       !$acc end kernels
       Return
    End If

    ibm_active = ibm_input_mode >= 1 .And. Allocated(Umask_cc)
    If ( ibm_active .And. .Not. Allocated(sgs_ibm_band) ) Call build_sgs_ibm_band
    inv_dx = 1d0 / dx
    dx2    = dx * dx

    !$acc parallel loop collapse(2) present(U_,V_,W_,S_,kappa_t,nu_t,y,yg,z,zg,xg)
    Do k = 2, nzg-1
       Do j = 2, nyg-1
          dz2       = ( z(k) - z(k-1) )**2
          inv_dz_k  = 1d0 / ( z(k) - z(k-1) )
          inv_zg_km = 1d0 / ( zg(k)   - zg(k-1) )
          inv_zg_kp = 1d0 / ( zg(k+1) - zg(k  ) )
          dy2_c     = ( y(j) - y(j-1) )**2
          inv_dy_j  = 1d0 / ( y(j) - y(j-1) )
          inv_yg_jm = 1d0 / ( yg(j)   - yg(j-1) )
          inv_yg_jp = 1d0 / ( yg(j+1) - yg(j  ) )

          Do i = 2, nxg-1
             a11 = ( U_(i,j,k)   - U_(i-1,j,k) ) * inv_dx
             a21 = 0.5d0*( ( U_(i,j,k) - U_(i,j-1,k) ) * inv_yg_jm + ( U_(i,j+1,k) - U_(i,j,k) ) * inv_yg_jp )
             a31 = 0.5d0*( ( U_(i,j,k) - U_(i,j,k-1) ) * inv_zg_km + ( U_(i,j,k+1) - U_(i,j,k) ) * inv_zg_kp )
             a12 = 0.5d0*( ( V_(i,j,k) - V_(i-1,j,k) ) + ( V_(i+1,j,k) - V_(i,j,k) ) ) * inv_dx
             a22 = ( V_(i,j,k)   - V_(i,j-1,k) ) * inv_dy_j
             a32 = 0.5d0*( ( V_(i,j,k) - V_(i,j,k-1) ) * inv_zg_km + ( V_(i,j,k+1) - V_(i,j,k) ) * inv_zg_kp )
             a13 = 0.5d0*( ( W_(i,j,k) - W_(i-1,j,k) ) + ( W_(i+1,j,k) - W_(i,j,k) ) ) * inv_dx
             a23 = 0.5d0*( ( W_(i,j,k) - W_(i,j-1,k) ) * inv_yg_jm + ( W_(i,j+1,k) - W_(i,j,k) ) * inv_yg_jp )
             a33 = ( W_(i,j,k)   - W_(i,j,k-1) ) * inv_dz_k

             s1 = ( S_(i+1,j,k) - S_(i-1,j,k) ) / ( xg(i+1) - xg(i-1) )
             s2 = ( S_(i,j+1,k) - S_(i,j-1,k) ) / ( yg(j+1) - yg(j-1) )
             s3 = ( S_(i,j,k+1) - S_(i,j,k-1) ) / ( zg(k+1) - zg(k-1) )

             num = -( dx2*s1*(a11*s1 + a12*s2 + a13*s3) + dy2_c*s2*(a21*s1 + a22*s2 + a23*s3) + dz2*s3*(a31*s1 + a32*s2 + a33*s3) )
             den = s1*s1 + s2*s2 + s3*s3
             If ( den > 1d-30 .And. num > 0d0 ) Then
                kappa_t(i,j,k) = C_amd * num / den
             Else
                kappa_t(i,j,k) = 0d0
             End If
          End Do
       End Do
    End Do
    !$acc end parallel loop

    ! next to a body the momentum model's one-sided pass has no scalar counterpart: keep nu_t/Pr_t
    If ( ibm_active ) Then
       !$acc parallel loop collapse(3) present(kappa_t,nu_t,sgs_ibm_band)
       Do k = 2, nzg-1
          Do j = 2, nyg-1
             Do i = 2, nxg-1
                If ( sgs_ibm_band(i,j,k) ) kappa_t(i,j,k) = nu_t(i,j,k) * kap_inv
             End Do
          End Do
       End Do
       !$acc end parallel loop
    End If

    Call sgs_fill_ghosts(kappa_t, ibm_active)

  End Subroutine compute_scalar_diffusivity


  !> Wall, periodic, IBM-solid and rank-seam ghost values of a cell-centred SGS diffusivity (nu_t or AMD
  !  kappa_t)
  Subroutine sgs_fill_ghosts(f_, zero_solid)

    Real(Int64), Dimension(nxg, nyg, nzg), Intent(InOut) :: f_
    Logical,                               Intent(In)    :: zero_solid
    Logical        :: is_first, is_last
    Integer(Int32) :: partner
    Real   (Int64) :: xbuf2(nyg,nzg), xbuf3(nyg,nzg,2)

#ifdef GPU_POISSON
    Call vreman_ghosts_device(f_, zero_solid)
#else
    ! No-slip flat walls only (j=1, nyg): zero the field so it does not pollute the Robin BC; free-slip
    ! (bc_face_y*==2) is left alone
    If ( y_bc_type == 1 .And. bc_face_ylo == 1 ) f_(:,  1,:) = 0d0
    If ( y_bc_type == 1 .And. bc_face_yhi == 1 ) f_(:,nyg,:) = 0d0

    ! y-periodicity (y_bc_type==0): fill ghost planes j=1,nyg (never written by
    ! Pass 1/2 above, but read by compute_rhs_v/w at the y boundaries)
    If ( y_bc_type == 0 ) Then
       f_(:,    1,:) = f_(:,nyg-2,:)
       f_(:,nyg-1,:) = f_(:,    2,:)
       f_(:,nyg  ,:) = f_(:,    3,:)
    End If

    ! IBM solid cells: suppress the SGS stress (after the y fills, before the halo copy)
    If ( zero_solid ) Then
       Where ( Umask_cc(2:nxg-1, 2:nyg-1, 2:nzg-1) < 0.5d0 )
          f_(2:nxg-1, 2:nyg-1, 2:nzg-1) = 0d0
       End Where
    End If

    ! x-periodicity: fill ghost planes i=1,nxg (never written above, but read by compute_rhs_v/w at the x
    ! boundaries) (x-split: seam planes come from the x-neighbour, the wrap from the partner rank at the
    ! opposite domain edge)
    Call update_ghost_interior_planes_x(f_, 2)
    Call x_periodic_partner(is_first, is_last, partner)
    If ( x_bc_type == 0 ) Then
       ! periodic: ghost 1 <- cell nxg-2, cells nxg-1,nxg <- 2,3 (nxg-1 duplicates 2), as apply_periodic_bc_x
       If ( is_first .And. is_last ) Then
          f_(1,    :,:) = f_(nxg-2,:,:)
          f_(nxg-1,:,:) = f_(2,    :,:)
          f_(nxg,  :,:) = f_(3,    :,:)
       Else If ( is_first ) Then
          xbuf3(:,:,1) = f_(2,:,:);  xbuf3(:,:,2) = f_(3,:,:)
          Call Mpi_sendrecv(xbuf3, 2*nyg*nzg, Mpi_real8, partner, 13, xbuf2, nyg*nzg, Mpi_real8, partner, 14, &
                            MPI_COMM_WORLD, istat, ierr)
          f_(1,:,:) = xbuf2
       Else If ( is_last ) Then
          xbuf2 = f_(nxg-2,:,:)
          Call Mpi_sendrecv(xbuf2, nyg*nzg, Mpi_real8, partner, 14, xbuf3, 2*nyg*nzg, Mpi_real8, partner, 13, &
                            MPI_COMM_WORLD, istat, ierr)
          f_(nxg-1,:,:) = xbuf3(:,:,1);  f_(nxg,:,:) = xbuf3(:,:,2)
       End If
    Else
       ! inflow/outflow: zero-gradient at the domain edges (interior-seam planes were filled by the exchange
       ! above)
       If ( is_first ) f_(1,  :,:) = f_(2,    :,:)
       If ( is_last  ) f_(nxg,:,:) = f_(nxg-1,:,:)
    End If
    ! Ring exchange for intermediate ranks (host-only); rank-0/rank-(nprocs-1) wrap handled below.
    Call update_ghost_interior_planes(f_, 4)

    If ( z_bc_type == 0 ) Then
       Call apply_periodic_bc_z(f_, 4)
    Else
       ! No-slip/wall-model z walls: zero the field at the wall faces (like the y walls) instead of wrapping;
       ! only the ranks owning z=0/z=Lz act
       Call z_periodic_partner(is_first, is_last, partner)
       If ( is_first ) f_(:,:,1)   = 0d0
       If ( is_last  ) f_(:,:,nzg) = 0d0
    End If
#endif

  End Subroutine sgs_fill_ghosts



#ifdef GPU_POISSON
  !> Device-resident equivalent of compute_vreman's post-Pass-1 host block (non-IBM): y ghosts, x seam + wrap,
  !  z seam + BC on nu_t_, MPI on device buffers
  Subroutine vreman_ghosts_device(nu_t_, zero_solid)

    Real(Int64), Dimension(nxg, nyg, nzg), Intent(InOut) :: nu_t_
    Logical,                               Intent(In)    :: zero_solid
    Logical        :: is_first, is_last
    Integer(Int32) :: partner, i, j, k

    ! No-slip flat walls only (j=1, nyg): zero nu_t; leave it alone at a free-slip boundary
    If ( y_bc_type == 1 .And. bc_face_ylo == 1 ) Then
       !$acc kernels present(nu_t_)
       nu_t_(:,  1,:) = 0d0
       !$acc end kernels
    End If
    If ( y_bc_type == 1 .And. bc_face_yhi == 1 ) Then
       !$acc kernels present(nu_t_)
       nu_t_(:,nyg,:) = 0d0
       !$acc end kernels
    End If
    If ( y_bc_type == 0 ) Then
       !$acc kernels present(nu_t_)
       nu_t_(:,    1,:) = nu_t_(:,nyg-2,:)
       nu_t_(:,nyg-1,:) = nu_t_(:,    2,:)
       nu_t_(:,nyg  ,:) = nu_t_(:,    3,:)
       !$acc end kernels
    End If

    ! IBM solid cells: suppress SGS stress to avoid polluting adjacent fluid (after the y fills, before the
    ! halos copy the values)
    If ( zero_solid ) Then
       !$acc parallel loop collapse(3) present(nu_t_,Umask_cc)
       Do k = 2, nzg-1
          Do j = 2, nyg-1
             Do i = 2, nxg-1
                If ( Umask_cc(i,j,k) < 0.5d0 ) nu_t_(i,j,k) = 0d0
             End Do
          End Do
       End Do
    End If

    ! x seam planes from the x-neighbour, then the periodic wrap (single-plane, cell-centred)
    Call update_ghost_interior_planes_x(nu_t_, 2)
    Call x_periodic_partner(is_first, is_last, partner)
    If ( x_bc_type == 0 ) Then
       If ( is_first .And. is_last ) Then
          !$acc kernels present(nu_t_)
          nu_t_(1,    :,:) = nu_t_(nxg-2,:,:)
          nu_t_(nxg-1,:,:) = nu_t_(2,    :,:)
          nu_t_(nxg,  :,:) = nu_t_(3,    :,:)
          !$acc end kernels
       Else
          Call gpu_periodic_wrap(nu_t_, .True., 2, nxg, is_first, is_last, partner, 13)
       End If
    Else
       !$acc kernels present(nu_t_)
       If ( is_first ) nu_t_(1,  :,:) = nu_t_(2,    :,:)
       If ( is_last  ) nu_t_(nxg,:,:) = nu_t_(nxg-1,:,:)
       !$acc end kernels
    End If

    ! z seam planes, then z periodic wrap or wall zeroing
    Call update_ghost_interior_planes(nu_t_, 2)
    If ( z_bc_type == 0 ) Then
       Call apply_periodic_bc_z(nu_t_, 4)
    Else
       Call z_periodic_partner(is_first, is_last, partner)
       !$acc kernels present(nu_t_)
       If ( is_first ) nu_t_(:,:,1)   = 0d0
       If ( is_last  ) nu_t_(:,:,nzg) = 0d0
       !$acc end kernels
    End If

  End Subroutine vreman_ghosts_device
#endif

End Module sgs_models
