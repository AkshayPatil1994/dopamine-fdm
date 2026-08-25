!> Standalone verification driver for src/sem.f90's raw_eddy_curl_sum (sem_divergence_free=1); built via the CMake target verify_sem_divergence_free (see CMakeLists.txt), not part of the dopamine executable
Program verify_sem_divergence_free

  Use iso_fortran_env, Only : Int32, Int64
  Use mpi
  Use global
  Use synthetic_eddy_method

  Implicit None

  Real(Int64) :: dy_, dz_

  Call Mpi_Init(ierr)
  Call Mpi_Comm_rank(MPI_COMM_WORLD, myid, ierr)
  Call Mpi_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

  Call setup_grid_and_profile

  ! ---- homogeneous mode ----
  Call reset_sem_state
  sem_divergence_free = 0
  Call run_case('homogeneous, sem_divergence_free=0 (regression baseline)')

  Call reset_sem_state
  sem_divergence_free = 1
  Call run_case('homogeneous, sem_divergence_free=1')

  ! ---- inhomogeneous mode ----
  Call setup_inhomogeneous_sigma

  Call reset_sem_state
  sem_divergence_free = 0
  Call run_case('inhomogeneous, sem_divergence_free=0 (regression baseline)')

  Call reset_sem_state
  sem_divergence_free = 1
  Call run_case('inhomogeneous, sem_divergence_free=1')

  If ( myid == 0 ) Write(*,'(A)') 'verify_sem_divergence_free: all cases completed'

  Call Mpi_Finalize(ierr)

Contains

  !> Small synthetic grid + a shear-flow-like reference profile (not representative of the production staggered layout, only self-consistent enough to exercise sem.f90)
  Subroutine setup_grid_and_profile

    Integer(Int32) :: i

    ny_global = 20; nz_global = 16
    Allocate( y_global(ny_global), z_global(nz_global) )
    dy_ = 1d0 / Real(ny_global-1,8)
    Do i = 1, ny_global
       y_global(i) = Real(i-1,8)*dy_
    End Do
    dz_ = 1d0 / Real(nz_global-1,8)
    Do i = 1, nz_global
       z_global(i) = Real(i-1,8)*dz_
    End Do

    ny = ny_global; nz = nz_global; nyg = ny_global; nzg = nz_global
    Allocate( y(ny), z(nz), yg(nyg), zg(nzg) )
    y = y_global; z = z_global; yg = y_global; zg = z_global

    n_profile = 2
    Allocate( prof_y(2), prof_U(2), prof_R11(2), prof_R22(2), prof_R33(2), prof_R12(2) )
    prof_y   = (/ 0d0,   1d0   /)
    prof_U   = (/ 1d0,   1d0   /)
    prof_R11 = (/ 0.05d0, 0.05d0 /)
    prof_R22 = (/ 0.03d0, 0.03d0 /)
    prof_R33 = (/ 0.04d0, 0.04d0 /)
    prof_R12 = (/ -0.01d0, -0.01d0 /)   ! R11*R22-R12^2 = 0.0014 > 0, PSD

    inflow_type = 1
    Uconv_sem   = Sum(prof_U) / Real(n_profile,8)   ! matches init_inflow's own formula

    sem_n_eddies         = 150
    sem_length_scale      = 0.12d0
    sem_seed              = 12345
    sem_ensemble_samples  = 30
    sem_ensemble_periods  = 4
    sem_use_esem           = 1
    sem_eddy_placement     = 0
    n_sigma                = 0   ! homogeneous by default; setup_inhomogeneous_sigma overrides

  End Subroutine setup_grid_and_profile

  !> Switches to a 9-component sigma_ij(y) profile (ESEM Section 4.1) for the inhomogeneous-mode cases
  Subroutine setup_inhomogeneous_sigma

    n_sigma = 3
    Allocate( sig_y(3) )
    Allocate( sig_ux(3), sig_uy(3), sig_uz(3) )
    Allocate( sig_vx(3), sig_vy(3), sig_vz(3) )
    Allocate( sig_wx(3), sig_wy(3), sig_wz(3) )

    sig_y  = (/ 0d0, 0.5d0, 1d0 /)
    sig_ux = (/ 0.06d0, 0.10d0, 0.06d0 /)
    sig_uy = (/ 0.05d0, 0.12d0, 0.05d0 /)
    sig_uz = (/ 0.08d0, 0.08d0, 0.08d0 /)
    sig_vx = (/ 0.07d0, 0.09d0, 0.07d0 /)
    sig_vy = (/ 0.06d0, 0.10d0, 0.06d0 /)
    sig_vz = (/ 0.08d0, 0.08d0, 0.08d0 /)
    sig_wx = (/ 0.07d0, 0.09d0, 0.07d0 /)
    sig_wy = (/ 0.08d0, 0.08d0, 0.08d0 /)
    sig_wz = (/ 0.06d0, 0.11d0, 0.06d0 /)

  End Subroutine setup_inhomogeneous_sigma

  !> Deallocate per-eddy/ensemble state so place_eddies/build_ensemble_normalisation can be re-run for the next case
  Subroutine reset_sem_state

    If ( Allocated(eddy_x0)     ) Deallocate(eddy_x0)
    If ( Allocated(eddy_sig)    ) Deallocate(eddy_sig)
    If ( Allocated(eddy_smax)   ) Deallocate(eddy_smax)
    If ( Allocated(eddy_Tperiod)) Deallocate(eddy_Tperiod)
    If ( Allocated(ens_mean_uU) ) Deallocate(ens_mean_uU, ens_std_uU)
    If ( Allocated(ens_mean_uV) ) Deallocate(ens_mean_uV, ens_std_uV)
    If ( Allocated(ens_mean_vV) ) Deallocate(ens_mean_vV, ens_std_vV)
    If ( Allocated(ens_mean_wW) ) Deallocate(ens_mean_wW, ens_std_wW)
    If ( Allocated(cdf_y)       ) Deallocate(cdf_y, cdf_val)
    n_cdf = 0

  End Subroutine reset_sem_state

  Subroutine run_case(label)

    Character(*), Intent(In) :: label

    If ( myid == 0 ) Write(*,'(A)') '=== ' // label // ' ==='

    Call place_eddies
    Call build_ensemble_normalisation

    Call check_divergence
    Call check_reynolds_stress

  End Subroutine run_case

  !> FD divergence check: dv/dy, dw/dz by central FD in space; du/dx via Taylor's-hypothesis FD in time (du/dx = -(1/Uconv_sem)*du/dt), exact for this xk-based construction regardless of internal sign conventions
  Subroutine check_divergence

    Real(Int64), Parameter :: hy = 1d-6, hz = 1d-6, ht = 1d-6
    Real(Int64) :: yc, zc, t0
    Real(Int64) :: uP,vP,wP, uM,vM,wM, dvdy, dwdz, dudx, divg, umag
    Real(Int64) :: max_divg, max_umag
    Integer(Int32) :: isamp

    max_divg = 0d0; max_umag = 0d0

    Do isamp = 1, 5
       yc = 0.15d0 + 0.15d0*Real(isamp-1,8)
       zc = 0.20d0 + 0.10d0*Real(isamp-1,8)
       t0 = 0.37d0 + 0.05d0*Real(isamp-1,8)

       Call raw_eddy_sum(yc, zc, t0, uP, vP, wP)
       umag = Sqrt(uP*uP+vP*vP+wP*wP)
       max_umag = Max(max_umag, umag)

       Call raw_eddy_sum(yc+hy, zc, t0, uP, vP, wP)
       Call raw_eddy_sum(yc-hy, zc, t0, uM, vM, wM)
       dvdy = (vP-vM)/(2d0*hy)

       Call raw_eddy_sum(yc, zc+hz, t0, uP, vP, wP)
       Call raw_eddy_sum(yc, zc-hz, t0, uM, vM, wM)
       dwdz = (wP-wM)/(2d0*hz)

       Call raw_eddy_sum(yc, zc, t0+ht, uP, vP, wP)
       Call raw_eddy_sum(yc, zc, t0-ht, uM, vM, wM)
       dudx = -(uP-uM)/(2d0*ht) / Uconv_sem

       divg = Abs(dudx+dvdy+dwdz)
       max_divg = Max(max_divg, divg)
    End Do

    If ( myid == 0 ) Write(*,'(A,E12.4,A,E12.4)') &
         '  max |div(u'')| over 5 sample points = ', max_divg, ', max |u''| = ', max_umag

  End Subroutine check_divergence

  !> Statistical smoke test: sample sem_fluctuation over time at one (y,z) point, compare realized second moments to the target profile
  Subroutine check_reynolds_stress

    Integer(Int32), Parameter :: nsamp = 4000
    Integer(Int32) :: jc, kc
    Real(Int64) :: yc, zc, t0, dt_s
    Real(Int64) :: up, vp, wp
    Real(Int64) :: su, sv, sw, suu, svv, sww, suv
    Integer(Int32) :: s

    ! (j,k) must index the same grid point as (yc,zc) -- ens_mean/ens_std lookups are per-grid-point (see build_ensemble_normalisation); this synthetic grid has y=yg=z=zg so one (jc,kc) works for all three components
    jc = nyg/2; kc = nzg/2
    yc = yg(jc); zc = zg(kc)
    dt_s = Maxval(eddy_Tperiod) * Real(sem_ensemble_periods,8) / Real(nsamp,8)

    su=0d0; sv=0d0; sw=0d0; suu=0d0; svv=0d0; sww=0d0; suv=0d0
    Do s = 0, nsamp-1
       t0 = dt_s * Real(s,8)
       Call sem_fluctuation(1, jc, kc, yc, zc, t0, up, vp, wp)
       su = su + up; suu = suu + up*up
       Call sem_fluctuation(2, jc, kc, yc, zc, t0, up, vp, wp)
       sv = sv + vp; svv = svv + vp*vp
       Call sem_fluctuation(3, jc, kc, yc, zc, t0, up, vp, wp)
       sw = sw + wp; sww = sww + wp*wp
    End Do

    If ( myid == 0 ) Write(*,'(A,3F9.5,A,3F9.5)') &
         '  realized (uu,vv,ww) = ', suu/nsamp-(su/nsamp)**2, svv/nsamp-(sv/nsamp)**2, sww/nsamp-(sw/nsamp)**2, &
         '   target (uu,vv,ww) = ', prof_R11(1), prof_R22(1), prof_R33(1)

  End Subroutine check_reynolds_stress

End Program verify_sem_divergence_free
