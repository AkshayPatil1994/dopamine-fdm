!     Module for LES wall-models
Module wallmodel

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use mpi
  Use interpolation
  Use boundary_conditions

  ! prevent implicit typing
  Implicit None

  ! Log-law constants
  Real(Int64), Parameter :: kappa_wm  = 0.41d0   ! von Karman constant
  Real(Int64), Parameter :: B_wm      = 5.2d0    ! log-law intercept
  Integer(Int32), Parameter :: n_iter_wm = 20    ! Newton iterations

  ! Reichardt (1951) composite law-of-the-wall constants
  Real(Int64), Parameter :: REICH_A1  = 11.0d0
  Real(Int64), Parameter :: REICH_A2  = 0.33d0
  Real(Int64), Parameter :: c_buf_wm  = B_wm - Log(kappa_wm)/kappa_wm

  ! Minimum y/z0 ratio for a grid point to be trusted as the rough-EQWM matching
  ! height (see j_match_ylo/yhi in global.f90)
  Real(Int64), Parameter :: MATCH_RATIO_MIN = 20d0

Contains

  !> Reichardt (1951) smooth u+(y+) law of the wall
  Pure Subroutine reichardt_uplus(yplus, uplus, duplus)
    !$acc routine seq
    Real(Int64), Intent(In)  :: yplus
    Real(Int64), Intent(Out) :: uplus, duplus
    Real(Int64) :: e1, e2, a1i

    a1i = 1d0 / REICH_A1
    e1  = Exp(-yplus * a1i)
    e2  = Exp(-REICH_A2 * yplus)

    uplus  = Log(1d0 + kappa_wm*yplus) / kappa_wm &
           + c_buf_wm * (1d0 - e1 - yplus*a1i*e2)
    duplus = 1d0 / (1d0 + kappa_wm*yplus) &
           + c_buf_wm * (a1i*e1 - a1i*e2 + yplus*a1i*REICH_A2*e2)
  End Subroutine reichardt_uplus

  !> Newton solve for u_tau from (u_ref, y_ref) via the Reichardt profile
  Subroutine solve_u_tau_reichardt(u_ref, y_ref, u_tau)
    !$acc routine seq
    Real(Int64), Intent(In)  :: u_ref, y_ref
    Real(Int64), Intent(Out) :: u_tau

    Real(Int64) :: yplus, uplus, duplus, f, fp, u_tau_new
    Integer(Int32) :: iter

    u_tau = Max(1d-3*u_ref, Sqrt(nu*u_ref/Max(y_ref,1d-14)))

    Do iter = 1, n_iter_wm
       yplus = u_tau * y_ref / nu
       Call reichardt_uplus(yplus, uplus, duplus)

       ! F(u_tau) = u_ref/u_tau - u+(y+) = 0; dF/du_tau adds the y+ chain term.
       f  =  u_ref / Max(u_tau, 1d-20) - uplus
       fp = -u_ref / Max(u_tau**2, 1d-20) - duplus * y_ref / nu

       u_tau_new = u_tau - f / ( fp + Sign(1d-14, fp) )
       u_tau_new = Min(Max(u_tau_new, 0.1d0*u_tau), 10d0*u_tau)

       If ( Abs(u_tau_new - u_tau) < 1d-6*u_tau_new ) Then
          u_tau = u_tau_new
          Exit
       End If
       u_tau = u_tau_new
    End Do

  End Subroutine solve_u_tau_reichardt

  !> Explicit fully-rough log-law u_tau from (u_ref, y_ref, z0); no viscous
  !  sublayer term -- valid once y_ref clears the roughness sublayer (y_ref >> z0)
  Pure Subroutine solve_u_tau_rough(u_ref, y_ref, z0, u_tau)
    !$acc routine seq
    Real(Int64), Intent(In)  :: u_ref, y_ref, z0
    Real(Int64), Intent(Out) :: u_tau

    Real(Int64) :: log_ratio

    log_ratio = Log( Max(y_ref, 2d0*z0) / Max(z0, 1d-8) )
    u_tau     = kappa_wm * u_ref / Max(log_ratio, 1d-3)

  End Subroutine solve_u_tau_rough

  !> Dispatch the friction-velocity solve: smooth Reichardt EQWM (mode 1) or
  !  rough z0 log law (mode 2); z0 is unused (but still passed) in mode 1
  Subroutine solve_u_tau_wall(u_ref, y_ref, z0, u_tau)
    !$acc routine seq
    Real(Int64), Intent(In)  :: u_ref, y_ref, z0
    Real(Int64), Intent(Out) :: u_tau

    If ( flat_wall_model_flag == 2 ) Then
       Call solve_u_tau_rough(u_ref, y_ref, z0, u_tau)
    Else
       Call solve_u_tau_reichardt(u_ref, y_ref, u_tau)
    End If

  End Subroutine solve_u_tau_wall

  !> Explicit fully-rough log law for the temperature scale theta_tau, from the
  !  wall-normal potential-temperature difference (T_ref - T_wall) and the thermal
  !  roughness z0h. Neutral limit only (psi_h=0) -- decoupled from u_tau, since in
  !  this limit the momentum and thermal log laws don't interact.
  Pure Subroutine solve_theta_tau_rough(dtheta, y_ref, z0h, theta_tau)
    !$acc routine seq
    Real(Int64), Intent(In)  :: dtheta, y_ref, z0h
    Real(Int64), Intent(Out) :: theta_tau

    Real(Int64) :: log_ratio

    log_ratio = Log( Max(y_ref, 2d0*z0h) / Max(z0h, 1d-8) )
    theta_tau = kappa_wm * dtheta / Max(log_ratio, 1d-3)

  End Subroutine solve_theta_tau_rough

  !> Businger-Dyer MOST stability functions psi_m/psi_h(zeta), zeta=y/L.
  !  Unstable (Paulson 1970 / Dyer 1974 integrated form), stable (linear).
  !  zeta is clipped to +-2 -- these empirical forms aren't trusted, and can
  !  blow up, outside roughly that range.
  Pure Subroutine most_stability_functions(zeta, psi_m, psi_h)
    !$acc routine seq
    Real(Int64), Intent(In)  :: zeta
    Real(Int64), Intent(Out) :: psi_m, psi_h

    Real(Int64) :: x, zeta_c

    zeta_c = Max(Min(zeta, 2d0), -2d0)

    If ( zeta_c < 0d0 ) Then
       x     = (1d0 - 16d0*zeta_c)**0.25d0
       psi_m = 2d0*Log((1d0+x)/2d0) + Log((1d0+x**2)/2d0) - 2d0*Atan(x) + pi/2d0
       psi_h = 2d0*Log((1d0+x**2)/2d0)
    Else
       psi_m = -5d0*zeta_c
       psi_h = -5d0*zeta_c
    End If

  End Subroutine most_stability_functions

  !> Coupled MOST solve for (u_tau, theta_tau) and the Obukhov length L, via
  !  fixed-point iteration of the Businger-Dyer-corrected log laws:
  !    u_tau     = kappa*u_ref  / [ln(y_ref/z0 ) - psi_m(y_ref/L)]
  !    theta_tau = kappa*dtheta / [ln(y_ref/z0h) - psi_h(y_ref/L)]
  !    L         = u_tau^2 * T_ref / (kappa*grav*theta_tau)
  !  L is seeded from the caller's persisted value (near-converged already once
  !  the flow is quasi-steady, so this typically only needs 1-2 passes) and
  !  updated in place. Reduces exactly to solve_u_tau_rough/solve_theta_tau_rough
  !  in the neutral limit (psi_m=psi_h=0).
  Subroutine solve_most(u_ref, dtheta, y_ref, z0, z0h, u_tau, theta_tau, L)
    !$acc routine seq
    Real(Int64), Intent(In)    :: u_ref, dtheta, y_ref, z0, z0h
    Real(Int64), Intent(Out)   :: u_tau, theta_tau
    Real(Int64), Intent(InOut) :: L

    Integer(Int32), Parameter :: MOST_MAX_ITER = 10
    Integer(Int32) :: iter
    Real(Int64) :: zeta, psi_m, psi_h, log_m, log_h, L_new

    Do iter = 1, MOST_MAX_ITER
       zeta = y_ref / Sign(Max(Abs(L), 1d-3), L)
       Call most_stability_functions(zeta, psi_m, psi_h)

       log_m = Log( Max(y_ref, 2d0*z0 ) / Max(z0,  1d-8) ) - psi_m
       log_h = Log( Max(y_ref, 2d0*z0h) / Max(z0h, 1d-8) ) - psi_h

       u_tau     = kappa_wm * u_ref  / Max(log_m, 1d-3)
       theta_tau = kappa_wm * dtheta / Max(log_h, 1d-3)

       If ( Abs(theta_tau) > 1d-12 ) Then
          L_new = u_tau**2 * T_ref / (kappa_wm * grav * theta_tau)
       Else
          L_new = Sign(1d10, L)   ! no flux -> neutral, huge |L|
       End If
       L_new = Sign( Max(Abs(L_new), 1d-3), L_new )

       If ( Abs(L_new - L) < 1d-3*Abs(L_new) ) Then
          L = L_new
          Exit
       End If
       L = L_new
    End Do

  End Subroutine solve_most

  !              Select wall model
  Subroutine compute_wall_model(U_,V_,W_,nu_t_)

    Real(Int64), Dimension(nx,  nyg, nzg), Intent(InOut) :: U_
    Real(Int64), Dimension(nxg, ny,  nzg), Intent(InOut) :: V_
    Real(Int64), Dimension(nxg, nyg, nz ), Intent(InOut) :: W_
    Real(Int64), Dimension(nxg, nyg, nzg), Intent(In)    :: nu_t_

    If ( y_bc_type == 0 ) Return   ! no wall model meaning for periodic y

    If ( flat_wall_model_flag == 1 .Or. flat_wall_model_flag == 2 ) Then
       ! Flat-wall log-law EQWM (smooth Reichardt or rough z0): compute alpha from local u_tau
       Call compute_flat_wall_eqwm(U_, W_)
    Else
       ! Constant Robin alpha (no-slip or free-slip depending on bc_face_y*)
       Call compute_constant_alpha
    End If
    ! alpha_y = 0 enforced inside compute_constant_alpha and compute_flat_wall_eqwm

    ! Rough-wall thermal coupling (neutral EQWM): only meaningful alongside the rough
    ! momentum EQWM, and only guarded here since Tscal is unallocated when Boussinesq is off
    If ( boussinesq_flag >= 1 .And. flat_wall_model_flag == 2 .And. &
         ( T_bc_bot == 2 .Or. T_bc_top == 2 ) ) Then
       Call compute_flat_wall_thermal_eqwm(U_, W_, Tscal)
    End If

    ! compute_pseudo_pressure_bc_for_robin_bc is host-only and reads alpha_y just written on-device
    !$acc update host(alpha_y)

    ! Compute boundary conditions for pseudo-pressure (flat walls)
    Call compute_pseudo_pressure_bc_for_robin_bc

    ! IBM surface wall model (ghost-cell EQWM)
    If ( ibm_input_mode >= 1 .And. ibm_wall_model_flag == 1 ) Then
       Call compute_ibm_wall_model(U_, V_, W_, nu_t_)
    End If

  End Subroutine compute_wall_model

  !> Equilibrium wall model for IBM ghost cells
  Subroutine compute_ibm_wall_model(U_, V_, W_, nu_t_)

    Use ibm, Only : ghost_u_idx, ghost_u_img, ghost_u_wgt, ghost_u_ref, n_ghost_u, &
                    ghost_u_nrm, ghost_u_yref,                                      &
                    ghost_v_idx, ghost_v_img, ghost_v_wgt, ghost_v_ref, n_ghost_v, &
                    ghost_v_nrm, ghost_v_yref,                                      &
                    ghost_w_idx, ghost_w_img, ghost_w_wgt, ghost_w_ref, n_ghost_w, &
                    ghost_w_nrm, ghost_w_yref,                                      &
                    trilinear_interp_u, trilinear_interp_v, trilinear_interp_w, &
                    U_wall, V_wall, W_wall

    Real(Int64), Dimension(nx,  nyg, nzg), Intent(InOut) :: U_
    Real(Int64), Dimension(nxg, ny,  nzg), Intent(InOut) :: V_
    Real(Int64), Dimension(nxg, nyg, nz ), Intent(InOut) :: W_
    Real(Int64), Dimension(nxg, nyg, nzg), Intent(In)    :: nu_t_

    Integer(Int32) :: n, i, j, k, ir, jr, kr
    Real   (Int64) :: u_ref, u_tau, y_ref
    Real   (Int64) :: uI_x, uI_y, uI_z   ! image-point velocity (all 3 components)
    Real   (Int64) :: uI_n                ! wall-normal projection of image velocity
    Real   (Int64) :: uI_tx, uI_ty, uI_tz ! tangential part of image velocity
    Real   (Int64) :: u_tan               ! tangential magnitude at image point
    Real   (Int64) :: u_I_eqwm
    Real   (Int64) :: yplus, uplus_img, duplus_img
    Real   (Int64) :: nx_, ny_, nz_       ! unit inward wall normal

    !--- EQWM for U ghost cells ---
    !$acc parallel loop present(U_,V_,W_,ghost_u_idx,ghost_u_ref,ghost_u_nrm,ghost_u_yref,ghost_u_dGB)
    Do n = 1, n_ghost_u
       i  = ghost_u_idx(1,n);  j  = ghost_u_idx(2,n);  k  = ghost_u_idx(3,n)
       ir = ghost_u_ref(1,n);  jr = ghost_u_ref(2,n);  kr = ghost_u_ref(3,n)
       nx_ = ghost_u_nrm(1,n);  ny_ = ghost_u_nrm(2,n);  nz_ = ghost_u_nrm(3,n)

       ! Reference-cell velocity averaging
       uI_x = 0.5d0 * (U_(ir-1, jr, kr) + U_(ir, jr, kr))
       uI_y = 0.5d0 * (V_(ir, Max(jr-1,2), kr) + V_(ir, jr, kr))
       uI_z = 0.5d0 * (W_(ir, jr, Max(kr-1,2)) + W_(ir, jr, kr))

       ! Project out wall-normal component: u_tan = u_I - (u_I · n) n
       uI_n  = uI_x*nx_ + uI_y*ny_ + uI_z*nz_
       uI_tx = uI_x - uI_n*nx_
       uI_ty = uI_y - uI_n*ny_
       uI_tz = uI_z - uI_n*nz_
       u_tan = Sqrt(uI_tx**2 + uI_ty**2 + uI_tz**2)

       ! Wall-normal distance to reference cell (n_image_layers from surface)
       y_ref = ghost_u_yref(n)
       u_ref = u_tan

       ! Newton solve using the Reichardt profile and nu only (excluding
       ! nu_t avoids SGS contamination of the wall-law).
       Call solve_u_tau_reichardt(u_ref, y_ref, u_tau)

       ! EQWM image-point capping at y_ref
       yplus = Min(ghost_u_dGB(n), y_ref) * u_tau / nu
       Call reichardt_uplus(yplus, uplus_img, duplus_img)
       u_I_eqwm = u_tau * uplus_img

       ! Mirror formula: U_G = 2*U_wall - U_I_eqwm (r = 0.5 always).
       If ( u_tan > 1d-14 ) Then
          U_(i, j, k) = 2d0*U_wall - uI_tx * u_I_eqwm / u_tan
       Else
          U_(i, j, k) = U_wall
       End If
    End Do
    !$acc end parallel loop

    !--- EQWM for V ghost cells ---
    !$acc parallel loop present(U_,V_,W_,ghost_v_idx,ghost_v_ref,ghost_v_nrm,ghost_v_yref,ghost_v_dGB)
    Do n = 1, n_ghost_v
       i  = ghost_v_idx(1,n);  j  = ghost_v_idx(2,n);  k  = ghost_v_idx(3,n)
       ir = ghost_v_ref(1,n);  jr = ghost_v_ref(2,n);  kr = ghost_v_ref(3,n)
       nx_ = ghost_v_nrm(1,n);  ny_ = ghost_v_nrm(2,n);  nz_ = ghost_v_nrm(3,n)

       uI_x = 0.5d0 * (U_(ir-1, jr, kr) + U_(ir, jr, kr))
       uI_y = 0.5d0 * (V_(ir, Max(jr-1,2), kr) + V_(ir, jr, kr))
       uI_z = 0.5d0 * (W_(ir, jr, Max(kr-1,2)) + W_(ir, jr, kr))

       uI_n  = uI_x*nx_ + uI_y*ny_ + uI_z*nz_
       uI_tx = uI_x - uI_n*nx_
       uI_ty = uI_y - uI_n*ny_
       uI_tz = uI_z - uI_n*nz_
       u_tan = Sqrt(uI_tx**2 + uI_ty**2 + uI_tz**2)

       y_ref = ghost_v_yref(n)
       u_ref = u_tan

       Call solve_u_tau_reichardt(u_ref, y_ref, u_tau)

       yplus = Min(ghost_v_dGB(n), y_ref) * u_tau / nu
       Call reichardt_uplus(yplus, uplus_img, duplus_img)
       u_I_eqwm = u_tau * uplus_img

       If ( u_tan > 1d-14 ) Then
          V_(i, j, k) = 2d0*V_wall - uI_ty * u_I_eqwm / u_tan
       Else
          V_(i, j, k) = V_wall
       End If
    End Do
    !$acc end parallel loop

    !--- EQWM for W ghost cells ---
    !$acc parallel loop present(U_,V_,W_,ghost_w_idx,ghost_w_ref,ghost_w_nrm,ghost_w_yref,ghost_w_dGB)
    Do n = 1, n_ghost_w
       i  = ghost_w_idx(1,n);  j  = ghost_w_idx(2,n);  k  = ghost_w_idx(3,n)
       ir = ghost_w_ref(1,n);  jr = ghost_w_ref(2,n);  kr = ghost_w_ref(3,n)
       nx_ = ghost_w_nrm(1,n);  ny_ = ghost_w_nrm(2,n);  nz_ = ghost_w_nrm(3,n)

       uI_x = 0.5d0 * (U_(ir-1, jr, kr) + U_(ir, jr, kr))
       uI_y = 0.5d0 * (V_(ir, Max(jr-1,2), kr) + V_(ir, jr, kr))
       uI_z = 0.5d0 * (W_(ir, jr, Max(kr-1,2)) + W_(ir, jr, kr))

       uI_n  = uI_x*nx_ + uI_y*ny_ + uI_z*nz_
       uI_tx = uI_x - uI_n*nx_
       uI_ty = uI_y - uI_n*ny_
       uI_tz = uI_z - uI_n*nz_
       u_tan = Sqrt(uI_tx**2 + uI_ty**2 + uI_tz**2)

       y_ref = ghost_w_yref(n)
       u_ref = u_tan

       Call solve_u_tau_reichardt(u_ref, y_ref, u_tau)

       yplus = Min(ghost_w_dGB(n), y_ref) * u_tau / nu
       Call reichardt_uplus(yplus, uplus_img, duplus_img)
       u_I_eqwm = u_tau * uplus_img

       If ( u_tan > 1d-14 ) Then
          W_(i, j, k) = 2d0*W_wall - uI_tz * u_I_eqwm / u_tan
       Else
          W_(i, j, k) = W_wall
       End If
    End Do
    !$acc end parallel loop

  End Subroutine compute_ibm_wall_model

  !  Flat-wall equilibrium wall model — local per-point EQWM
  Subroutine compute_flat_wall_eqwm(U_, W_)

    ! Local per-point EQWM for flat channel walls
    ! u_tau is solved from nu only -- excluding nu_t avoids SGS contamination of the wall-law, mirrors compute_ibm_wall_model

    Real(Int64), Dimension(nx,  nyg, nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg, nyg, nz ), Intent(In) :: W_

    Integer(Int32) :: i, k
    Real   (Int64) :: u_ref, u_match
    Real   (Int64) :: y_ref_lo, y_ref_hi, y_match_lo, y_match_hi
    Real   (Int64) :: u_tau
    Real   (Int64) :: Delta_yg_lo, Delta_yg_hi
    Real   (Int64) :: alpha_lo, alpha_hi
    Real   (Int64) :: W_at_pt, W_match

    ! Wall-normal reference distances (same for the entire wall plane). y_ref_*
    ! is the actual first-interior-cell height the Robin BC extrapolates from
    ! (j=2/nyg-1, always); y_match_* is where u_tau/theta_tau are sampled from --
    ! the same point when flat_wall_model_flag/=2 (j_match_*=2/nyg-1, unshifted),
    ! or a point further from the wall for the rough EQWM (see j_match_ylo/yhi).
    y_ref_lo    = yg(2)
    Delta_yg_lo = yg(2) - yg(1)
    y_ref_hi    = Ly - yg(nyg-1)
    Delta_yg_hi = yg(nyg) - yg(nyg-1)
    y_ref_lo    = Max(y_ref_lo, 1d-14)
    y_ref_hi    = Max(y_ref_hi, 1d-14)
    y_match_lo  = Max(yg(j_match_ylo), 1d-14)
    y_match_hi  = Max(Ly - yg(j_match_yhi), 1d-14)

    !  alpha_x : Robin slip-length for U (x-faces, nx × nyg × nzg)
    !$acc parallel loop collapse(2) present(U_,W_,alpha_x)
    Do k = 2, nzg-1
       Do i = 2, nx-1

          ! ---- bottom wall: tangential speed at first interior cell (j=2) ----
          If ( bc_face_ylo /= 2 ) Then
             ! W approximated by averaging the two bracketing z-faces.
             W_at_pt  = 0.5d0*(W_(i, 2, k-1) + W_(i, 2, k))
             u_ref    = Sqrt(U_(i, 2, k)**2 + W_at_pt**2)

             ! Matching-height sample for the u_tau solve (== u_ref/j=2 unless
             ! flat_wall_model_flag=2 shifted it further from the wall)
             W_match  = 0.5d0*(W_(i, j_match_ylo, k-1) + W_(i, j_match_ylo, k))
             u_match  = Sqrt(U_(i, j_match_ylo, k)**2 + W_match**2)

             ! Newton solve (smooth) or explicit log law (rough) for u_tau using nu only
             Call solve_u_tau_wall(u_match, y_match_lo, z0_ylo, u_tau)

             ! Robin alpha derivation -- always referenced to the actual first
             ! interior cell (u_ref/y_ref_lo), since that's what apply_Robin_bc_y
             ! extrapolates from; only the u_tau solve above uses the matching height
             alpha_lo = nu * u_ref / Max(u_tau**2, 1d-20) &
                      - Delta_yg_lo*0.5d0
             alpha_x(i, 1, k) = Max(alpha_lo, 0d0)
          Else
             alpha_x(i, 1, k) = 1.0e10_8   ! Neumann free-slip
          End If

          ! ---- top wall -------------------------------------------
          If ( bc_face_yhi /= 2 ) Then
             W_at_pt  = 0.5d0*(W_(i, nyg-1, k-1) + W_(i, nyg-1, k))
             u_ref    = Sqrt(U_(i, nyg-1, k)**2 + W_at_pt**2)

             W_match  = 0.5d0*(W_(i, j_match_yhi, k-1) + W_(i, j_match_yhi, k))
             u_match  = Sqrt(U_(i, j_match_yhi, k)**2 + W_match**2)

             Call solve_u_tau_wall(u_match, y_match_hi, z0_yhi, u_tau)

             alpha_hi = nu * u_ref / Max(u_tau**2, 1d-20) &
                      - Delta_yg_hi*0.5d0
             alpha_x(i, 2, k) = Max(alpha_hi, 0d0)
          Else
             alpha_x(i, 2, k) = 1.0e10_8   ! Neumann free-slip
          End If

       End Do
    End Do
    !$acc end parallel loop

    ! Fill x-halo planes (periodic in x)
    !$acc kernels present(alpha_x)
    alpha_x( 1, :, :) = alpha_x(   2, :, :)
    alpha_x(nx, :, :) = alpha_x(nx-1, :, :)
    ! Fill z-halo planes (copy from nearest interior; avoids extra MPI exchange)
    alpha_x(:, :,     1) = alpha_x(:, :,       2)
    alpha_x(:, :, nzg  ) = alpha_x(:, :, nzg-1  )
    !$acc end kernels

    !  alpha_z : Robin slip-length for W (z-faces, nxg × nyg × nz)
    !$acc parallel loop collapse(2) present(U_,W_,alpha_z)
    Do k = 2, nz-1
       Do i = 2, nxg-1

          ! ---- bottom wall ----------------------------------------
          If ( bc_face_ylo /= 2 ) Then
             ! U at this x-centre: average the two bracketing x-faces.
             u_ref    = Sqrt((0.5d0*(U_(i-1, 2, k) + U_(i, 2, k)))**2 &
                           + W_(i, 2, k)**2)

             u_match  = Sqrt((0.5d0*(U_(i-1, j_match_ylo, k) + U_(i, j_match_ylo, k)))**2 &
                           + W_(i, j_match_ylo, k)**2)

             Call solve_u_tau_wall(u_match, y_match_lo, z0_ylo, u_tau)

             alpha_lo = nu * u_ref / Max(u_tau**2, 1d-20) &
                      - Delta_yg_lo*0.5d0
             alpha_z(i, 1, k) = Max(alpha_lo, 0d0)
          Else
             alpha_z(i, 1, k) = 1.0e10_8   ! Neumann free-slip
          End If

          ! ---- top wall -------------------------------------------
          If ( bc_face_yhi /= 2 ) Then
             u_ref    = Sqrt((0.5d0*(U_(i-1, nyg-1, k) + U_(i, nyg-1, k)))**2 &
                           + W_(i, nyg-1, k)**2)

             u_match  = Sqrt((0.5d0*(U_(i-1, j_match_yhi, k) + U_(i, j_match_yhi, k)))**2 &
                           + W_(i, j_match_yhi, k)**2)

             Call solve_u_tau_wall(u_match, y_match_hi, z0_yhi, u_tau)

             alpha_hi = nu * u_ref / Max(u_tau**2, 1d-20) &
                      - Delta_yg_hi*0.5d0
             alpha_z(i, 2, k) = Max(alpha_hi, 0d0)
          Else
             alpha_z(i, 2, k) = 1.0e10_8   ! Neumann free-slip
          End If

       End Do
    End Do
    !$acc end parallel loop

    ! Fill x-halo planes (periodic in x)
    !$acc kernels present(alpha_z,alpha_y)
    alpha_z(   1, :, :) = alpha_z(     2, :, :)
    alpha_z( nxg, :, :) = alpha_z( nxg-1, :, :)
    ! Fill z-halo planes (copy from nearest interior)
    alpha_z(:, :,   1) = alpha_z(:, :,     2)
    alpha_z(:, :,  nz) = alpha_z(:, :, nz-1)

    alpha_y = 0d0
    !$acc end kernels

  End Subroutine compute_flat_wall_eqwm

  !> Flat-wall rough-EQWM thermal coupling (neutral limit: psi_h=0).
  !  Computes the Robin slip-length alpha_T for Tscal's y-ghost cells from the
  !  z0h log law, using the same numerical device as the momentum alpha_x/alpha_z:
  !  alpha_T is chosen so that the discrete molecular-only flux (nu/Pr)*(T_ref-T_ghost)/dy
  !  reproduces the target kinematic heat flux Q = u_tau*theta_tau. Only active on
  !  walls where T_bc_bot/top==2; other walls leave alpha_T untouched (unused there).
  Subroutine compute_flat_wall_thermal_eqwm(U_, W_, T_)

    Real(Int64), Dimension(nx,  nyg, nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg, nyg, nz ), Intent(In) :: W_
    Real(Int64), Dimension(nxg, nyg, nzg), Intent(In) :: T_

    Integer(Int32) :: i, k
    Real   (Int64) :: u_match, u_tau, theta_tau, q_target
    Real   (Int64) :: y_ref_lo, y_ref_hi, y_match_lo, y_match_hi, Delta_yg_lo, Delta_yg_hi
    Real   (Int64) :: U_at_pt, W_at_pt, T_here, T_match, alpha_lo, alpha_hi

    y_ref_lo    = Max(yg(2), 1d-14)
    Delta_yg_lo = yg(2) - yg(1)
    y_ref_hi    = Max(Ly - yg(nyg-1), 1d-14)
    Delta_yg_hi = yg(nyg) - yg(nyg-1)
    y_match_lo  = Max(yg(j_match_ylo), 1d-14)
    y_match_hi  = Max(Ly - yg(j_match_yhi), 1d-14)

    ! u_tau/theta_tau are solved at the matching height (j_match_*, further from
    ! the wall than j=2/nyg-1 when the near-wall grid is fine relative to z0/z0h),
    ! via the coupled Businger-Dyer MOST iteration (reduces to the neutral log
    ! laws when the surface buoyancy flux is ~0). The alpha_T formula still
    ! references the actual first interior cell (T_here, y_ref_*), since that's
    ! what apply_Robin_bc_y_scalar extrapolates from.
    !
    ! Known scope limitation: this stability correction feeds alpha_T only.
    ! alpha_x/alpha_z (momentum) stay on the neutral rough z0 EQWM from
    ! compute_flat_wall_eqwm -- they sample at different index spaces (x-faces,
    ! z-faces) than this cell-centred pass, so consistently stability-correcting
    ! them needs their own persisted L state and is deferred.
    !$acc parallel loop collapse(2) present(U_,W_,T_,alpha_T,L_obukhov_ylo,L_obukhov_yhi,yg)
    Do k = 2, nzg-1
       Do i = 2, nxg-1

          ! ---- bottom wall ----
          If ( T_bc_bot == 2 ) Then
             U_at_pt = 0.5d0*(U_(i-1, j_match_ylo, k) + U_(i, j_match_ylo, k))
             W_at_pt = 0.5d0*(W_(i, j_match_ylo, Max(k-1,2)) + W_(i, j_match_ylo, k))
             u_match = Sqrt(U_at_pt**2 + W_at_pt**2)
             T_here  = T_(i, 2, k)
             T_match = T_(i, j_match_ylo, k)

             Call solve_most(u_match, T_match - T_wall_bot, y_match_lo, z0_ylo, z0h_ylo, &
                              u_tau, theta_tau, L_obukhov_ylo(i,k))

             q_target = u_tau * theta_tau
             If ( Abs(q_target) > 1d-12 ) Then
                alpha_lo = (nu/Pr) * (T_here - T_wall_bot) / q_target - Delta_yg_lo*0.5d0
             Else
                alpha_lo = 1.0e10_8   ! no resolved flux -> effectively adiabatic ghost
             End If
             alpha_T(i, 1, k) = Max(alpha_lo, 0d0)
          End If

          ! ---- top wall ----
          If ( T_bc_top == 2 ) Then
             U_at_pt = 0.5d0*(U_(i-1, j_match_yhi, k) + U_(i, j_match_yhi, k))
             W_at_pt = 0.5d0*(W_(i, j_match_yhi, Max(k-1,2)) + W_(i, j_match_yhi, k))
             u_match = Sqrt(U_at_pt**2 + W_at_pt**2)
             T_here  = T_(i, nyg-1, k)
             T_match = T_(i, j_match_yhi, k)

             Call solve_most(u_match, T_match - T_wall_top, y_match_hi, z0_yhi, z0h_yhi, &
                              u_tau, theta_tau, L_obukhov_yhi(i,k))

             q_target = u_tau * theta_tau
             If ( Abs(q_target) > 1d-12 ) Then
                alpha_hi = (nu/Pr) * (T_here - T_wall_top) / q_target - Delta_yg_hi*0.5d0
             Else
                alpha_hi = 1.0e10_8
             End If
             alpha_T(i, 2, k) = Max(alpha_hi, 0d0)
          End If

       End Do
    End Do
    !$acc end parallel loop

    ! Fill x-halo planes (periodic in x)
    !$acc kernels present(alpha_T)
    alpha_T(  1, :, :) = alpha_T(    2, :, :)
    alpha_T(nxg, :, :) = alpha_T(nxg-1, :, :)
    ! Fill z-halo planes (copy from nearest interior)
    alpha_T(:, :,   1) = alpha_T(:, :,     2)
    alpha_T(:, :, nzg) = alpha_T(:, :, nzg-1)
    !$acc end kernels

  End Subroutine compute_flat_wall_thermal_eqwm

  !> Set constant Robin BC alpha for y-walls (bc_face=1 -> Dirichlet no-slip alpha=0; bc_face=2 -> Neumann free-slip alpha=1e10)
  Subroutine compute_constant_alpha

    Real(Int64) :: alpha_lo, alpha_hi

    alpha_lo = Merge(1.0e10_8, 0d0, bc_face_ylo == 2)
    alpha_hi = Merge(1.0e10_8, 0d0, bc_face_yhi == 2)

    !$acc kernels present(alpha_x,alpha_y,alpha_z)
    alpha_x(1:nx,  1, 1:nzg) = alpha_lo
    alpha_x(1:nx,  2, 1:nzg) = alpha_hi
    alpha_y                  = 0d0
    alpha_z(1:nxg, 1, 1:nz ) = alpha_lo
    alpha_z(1:nxg, 2, 1:nz ) = alpha_hi
    !$acc end kernels

  End Subroutine compute_constant_alpha

  !            Compute Neumann boundary conditions
  !    for pseudo-pressure when slip-wall model is active
  Subroutine compute_pseudo_pressure_bc_for_robin_bc

    ! local variables
    Real   (Int64) :: a, b, c
    Real   (Int64) :: beta, Delta_r, alphad
    Real   (Int64) :: p_bc2, p_bc3, p_bcn, p_bcn1
    Integer(Int64) :: j    

    ! bottom wall
    j        = 2 
    a        = 1d0/( y(j)-y(j-1) )/( yg(j+1) - yg(j) )
    b        = 1d0/( y(j)-y(j-1) )*( -1d0/( yg(j+1) - yg(j) ) -1d0/( yg(j) - yg(j-1) ) )
    c        = 1d0/( y(j)-y(j-1) )/( yg(j) - yg(j-1) ) 
    Delta_r  = ( yg(2)-yg(1) )/( yg(3)-yg(2) )
    alphad   = alpha_y(2,1,2)
    beta     = alphad/(alphad + y(2)-y(1) )
    p_bc2    = 1d0 + beta*Delta_r
    p_bc3    =     - beta*Delta_r

    Dyy(2,2) = b + c*p_bc2
    Dyy(2,3) = a + c*p_bc3
    
    ! top wall
    j        = nyg-1
    a        = 1d0/( y(j)-y(j-1) )/( yg(j+1) - yg(j) )
    b        = 1d0/( y(j)-y(j-1) )*( -1d0/( yg(j+1) - yg(j) ) -1d0/( yg(j) - yg(j-1) ) )
    c        = 1d0/( y(j)-y(j-1) )/( yg(j) - yg(j-1) )     
    alphad   = -alpha_y(2,2,2)
    Delta_r  = ( yg(nyg) - yg(nyg-1) )/( yg(nyg-1) - yg(nyg-2) )
    beta     = alphad/( alphad - (y(ny)-y(ny-1)) )
    p_bcn    = 1d0 + beta*Delta_r
    p_bcn1   =     - beta*Delta_r
    
    Dyy(nyg-1,nyg-1) = b + a*p_bcn
    Dyy(nyg-1,nyg-2) = c + a*p_bcn1
    
  End Subroutine compute_pseudo_pressure_bc_for_robin_bc

End Module wallmodel

