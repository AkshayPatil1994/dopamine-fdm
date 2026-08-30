! Module to monitor status of the simulation
Module monitor

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use mpi
  Use projection, Only : check_divergence
  Use ibm

  ! prevent implicit typing
  Implicit None

Contains

  ! Output some key values during simulation
  Subroutine output_monitor

    Real(Int64)    :: local_max
    Real(Int64)    :: meanU, maxU
    Real(Int64)    :: max_divergence
    Integer(Int32) :: local_nan, global_nan
    Real(Int64)    :: local_sc, total_sc

    If ( Mod(istep, nmonitor) == 0 ) Then

       ! ── NaN / Inf guard: x /= x detects NaN; Abs(x) > 1d8 catches overflow ──
       local_nan = 0
       If ( Any(U(2:nx-1, 2:nyg-1, 2:nzg-1) /= U(2:nx-1, 2:nyg-1, 2:nzg-1)) .Or. &
            Any(V(2:nxg-1,2:ny-1,  2:nzg-1) /= V(2:nxg-1,2:ny-1,  2:nzg-1)) .Or. &
            Any(W(2:nxg-1,2:nyg-1, 2:nz-1 ) /= W(2:nxg-1,2:nyg-1, 2:nz-1 )) .Or. &
            Any(Abs(U(2:nx-1, 2:nyg-1, 2:nzg-1)) > 1d8) .Or. &
            Any(Abs(V(2:nxg-1,2:ny-1,  2:nzg-1)) > 1d8) .Or. &
            Any(Abs(W(2:nxg-1,2:nyg-1, 2:nz-1 )) > 1d8) ) Then
          local_nan = 1
       End If
       Call MPI_Allreduce(local_nan, global_nan, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
       If (global_nan /= 0) Then
          If (myid == 0) Then
             Write(*,'(A)') ' '
             Write(*,'(A,I0,A,E13.6)') ' FATAL: NaN or |U| > 1e8 at step ', istep, &
                  '  t = ', t
             Write(*,'(A)') ' Aborting simulation.'
             Write(*,'(A)') ' '
          End If
          Call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
       End If

       ! ── mean (volume-weighted, IBM-solid-excluded) and max of streamwise velocity ──
       Call compute_bulk_velocity(meanU)

       local_max = MaxVal( Abs(U(2:nx-1, 2:nyg-1, 2:nzg-1)) )
       Call MPI_Reduce(local_max, maxU, 1, MPI_real8, MPI_MAX, 0, MPI_COMM_WORLD, ierr)

       ! ── divergence ───────────────────────────────────────────────────
       Call check_divergence(max_divergence)

       time2 = MPI_WTIME()

       If ( myid == 0 ) Then
          Write(*,'(I8,3X,ES12.5,3X,ES11.4,3X,ES11.4,3X,ES10.3,3X,ES9.2,3X,ES9.2,3X,ES12.5,3X,F7.2)') &
               istep, t,                                                     &
               meanU,                                                        &
               maxU, max_divergence,                                         &
               cfl_conv_last, cfl_visc_last, dt, time2 - time1
       End If

       time1 = MPI_WTIME()

       ! scalar-mass reporting disabled to avoid log clutter

    End If

  End Subroutine output_monitor

  !   Output summary of initial parameters
  Subroutine summary

    Real(Int64)   :: u_tau_est, Re_tau_est, h_ref, visc_cfl, t_end_s, dy_max
    Character(40) :: str_sgs, str_ibm, str_ibmwm, str_wm, str_bc_top

    If ( myid==0 ) Then

       ! ── derived quantities ──────────────────────────────────────────
       If ( bc_face_yhi == 1 ) Then
          h_ref = 0.5d0 * Ly   ! half-height for closed channel
       Else
          h_ref = Ly            ! full height for open channel
       End If
       u_tau_est  = Sqrt( Max(Abs(dPdx) * h_ref, 0d0) )
       If ( nu > 0d0 .And. u_tau_est > 0d0 ) Then
          Re_tau_est = u_tau_est * h_ref / nu
          visc_cfl   = u_tau_est * dt / dymin
       Else
          Re_tau_est = 0d0
          visc_cfl   = 0d0
       End If
       t_end_s = Real(nsteps, Int64) * dt
       dy_max  = Maxval( yg_global(2:nyg_global) - yg_global(1:nyg_global-1) )

       ! ── label strings ───────────────────────────────────────────────
       Select Case (sgs_model)
       Case (0);        str_sgs   = '0  (DNS — no eddy viscosity)'
       Case (1);        str_sgs   = '1  (Vreman)'
       Case Default; Write(str_sgs,'(I2,A)') sgs_model, '  (unknown)'
       End Select

       Select Case (ibm_input_mode)
       Case (0);        str_ibm   = '0  (no IBM body — smooth walls)'
       Case (1);        str_ibm   = '1  (SDF-based ghost-cell IBM)'
       Case Default; Write(str_ibm,'(I2,A)') ibm_input_mode, '  (unknown)'
       End Select

       If (ibm_wall_model_flag == 0) Then
          str_ibmwm = '0  (no-slip, DNS)'
       Else
          str_ibmwm = '1  (log-law EQWM)'
       End If

       If (flat_wall_model_flag == 0) Then
          str_wm = '0  (no-slip, DNS)'
       Else
          str_wm = '1  (log-law EQWM)'
       End If

       If ( y_bc_type == 1 ) Then
          If (bc_face_yhi == 1) Then
             str_bc_top = '1  (no-slip — closed channel)'
          Else
             str_bc_top = '2  (free-slip — open channel)'
          End If
       Else
          str_bc_top = '0  (periodic, y_bc_type=0)'
       End If

       ! ── print ───────────────────────────────────────────────────────
       Write(*,'(A)') ' '
       Write(*,'(A)') ' +====================================================================+'
       Write(*,'(A)') '  SIMULATION PARAMETERS SUMMARY'
       Write(*,'(A)') ' +====================================================================+'
       Write(*,'(A,I0,A)') '  MPI ranks          :  ', nprocs, ' ranks'
       Write(*,'(A)') '  --------------------------------------------------------------------'
       Write(*,'(A)') '  DOMAIN'
       Write(*,'(A,3F10.4)') '    (Lx, Ly, Lz)     : ', Lx, Ly, Lz
       Write(*,'(A,3I6)')    '    faces  (nx,ny,nz) : ', nx_global, ny_global, nz_global
       Write(*,'(A,3I6)')    '    cells (nxm,nym,nzm): ', nxm_global, nym_global, nzm_global
       Write(*,'(A,I2)')     '    grid type          :  ', grid_type
       Write(*,'(A,E12.4)')  '    dx                 :  ', dxmin
       Write(*,'(A,E12.4,A,E12.4)') '    dy  min / max      :  ', dymin, '  /  ', dy_max
       Write(*,'(A,E12.4)')  '    dz                 :  ', dzmin
       Write(*,'(A)') '  --------------------------------------------------------------------'
       Write(*,'(A)') '  PHYSICS'
       Write(*,'(A,E12.4)')  '    nu                 :  ', nu
       Write(*,'(A,E12.4)')  '    dPdx               :  ', dPdx
       Write(*,'(A,E12.4)')  '    dPdz               :  ', dPdz
       If ( waveOmega_x > 0d0 ) Then
          Write(*,'(A,E12.4)') '    Ub_x               :  ', Ub_x
          Write(*,'(A,E12.4)') '    wave_omega_x       :  ', waveOmega_x
          Write(*,'(A,E12.4)') '    phi_wave_x         :  ', phi_wave_x
       End If
       If ( waveOmega_z > 0d0 ) Then
          Write(*,'(A,E12.4)') '    Ub_z               :  ', Ub_z
          Write(*,'(A,E12.4)') '    wave_omega_z       :  ', waveOmega_z
          Write(*,'(A,E12.4)') '    phi_wave_z         :  ', phi_wave_z
       End If
       If (Re_tau_est > 0d0) Then
          Write(*,'(A,F9.2)') '    Re_tau  (est.)     :  ', Re_tau_est
          Write(*,'(A,E12.4)') '    u_tau   (est.)     :  ', u_tau_est
       End If
       Write(*,'(A)') '  --------------------------------------------------------------------'
       Write(*,'(A)') '  TIME STEPPING'
       Write(*,'(A,E12.4)')  '    dt                 :  ', dt
       If (cfl_adaptive == 1) Then
          Write(*,'(A)') '    adaptive dt        :  ON'
          Write(*,'(A,F8.4)')  '    cfl_target         :  ', cfl_target
          Write(*,'(A,F8.4)')  '    cfl_safety         :  ', cfl_safety
          Write(*,'(A,E12.4)') '    dt_min             :  ', dt_min
          Write(*,'(A,E12.4)') '    dt_max             :  ', dt_max
       Else
          Write(*,'(A)') '    adaptive dt        :  OFF (fixed dt)'
       End If
       If (visc_cfl > 0d0) &
          Write(*,'(A,E12.4,A)') '    viscous CFL (est.) :  ', visc_cfl, '  (u_tau*dt/dy_min)'
       If ( sediment_flag >= 1 .And. ws > 0d0 ) Then
          Write(*,'(A,E12.4,A)') '    settling CFL (est.):  ', ws*dt/dymin, '  (ws*dt/dy_min)'
       End If
       Write(*,'(A,I8)')     '    nsteps             :  ', nsteps
       Write(*,'(A,F12.4)')  '    t_end              :  ', t_end_s
       Write(*,'(A,I8)')     '    nsave              :  ', nsave
       Write(*,'(A,I8)')     '    nmonitor           :  ', nmonitor
       Write(*,'(A)') '  --------------------------------------------------------------------'
       Write(*,'(A)') '  MODELS'
       Write(*,'(A,A)') '    SGS model          :  ', Trim(str_sgs)
       If (sgs_model == 1) &
          Write(*,'(A,E12.4)') '    Cs_vreman          :  ', Cs_vreman
       Write(*,'(A,A)') '    Flat-wall model    :  ', Trim(str_wm)
       Write(*,'(A,A)') '    IBM mode           :  ', Trim(str_ibm)
       If (ibm_input_mode >= 1) &
          Write(*,'(A,A)') '    IBM wall model     :  ', Trim(str_ibmwm)
       If (ibm_input_mode >= 1 .And. nsampling > 0) &
          Write(*,'(A,I8)') '    IBM force sampling :  every ', nsampling, ' steps'
       If ( y_bc_type == 1 ) Then
          Write(*,'(A,A)') '    Top wall BC        :  ', Trim(str_bc_top)
       Else
          Write(*,'(A,A)') '    y-direction BC     :  ', Trim(str_bc_top)
       End If
       If ( sediment_flag >= 1 ) Then
          Write(*,'(A)') '  --------------------------------------------------------------------'
          Write(*,'(A)') '  SEDIMENT'
          Write(*,'(A,E12.4)')  '    d_s                :  ', d_s
          Write(*,'(A,E12.4)')  '    rho_s              :  ', rho_s
          Write(*,'(A,E12.4)')  '    ws (Soulsby 1997)  :  ', ws
          Write(*,'(A,E12.4)')  '    Sc                 :  ', Sc
          Write(*,'(A,E12.4)')  '    Sc_t               :  ', Sc_t
          Write(*,'(A,E12.4)')  '    C_ref              :  ', C_ref
          Write(*,'(A,I2)')     '    sed_bc_bot         :  ', sed_bc_bot
       End If
       Write(*,'(A)') '  --------------------------------------------------------------------'
       Write(*,'(A)') '  FILES'
       If (restart == 0) Then
          Write(*,'(A,I2)')     '    ic_type            :  ', ic_type
          Write(*,'(A,F6.1,A)') '    noise              :  ', noise_percent, ' %'
       Else
          Write(*,'(A,A)')      '    restart from       :  ', Trim(filein)
       End If
       Write(*,'(A,A)')         '    output prefix      :  ', Trim(fileout)
       Write(*,'(A)') ' +====================================================================+'
       Write(*,'(A)') ' '
       Write(*,'(A,I0,A,F0.4)') '  Starting time loop   nsteps=', nsteps, '  t_end=', t_end_s
       Write(*,'(A)') ' '
       Write(*,'(A)') '    step             t           <U>        |U|max        |div|       CFL_c       CFL_v' // &
                      '             dt   wall(s)'
       Write(*,'(A)') '--------   ------------   -----------   -----------   ----------   ---------   ---------' // &
                      '   ------------   -------'

    End If

  End Subroutine summary

  !> Compute convective/viscous CFL (global max over MPI ranks)
  Subroutine compute_cfl(cfl_conv_out, cfl_visc_out)

    Real(Int64), Intent(Out) :: cfl_conv_out, cfl_visc_out

    Real(Int64)    :: local_conv, local_visc, conv_ijk, visc_ijk
    Real(Int64)    :: inv_dx, inv_dy, inv_dz, inv_dx2, inv_dy2, inv_dz2
    Real(Int64)    :: local_buf(2), global_buf(2)
    Integer(Int32) :: i, j, k

    local_conv = 0d0
    local_visc = 0d0

    ! x and z grids are uniform: hoist their inverses outside all loops.
    inv_dx  = 1d0 / dx
    inv_dz  = 1d0 / dz
    inv_dx2 = inv_dx * inv_dx
    inv_dz2 = inv_dz * inv_dz

    If ( ibm_input_mode >= 1 ) Then
       ! Kept as a separate branch (not merged behind a runtime .And. guard) since Fortran doesn't guarantee short-circuit evaluation and phi mustn't be touched when not IBM-active
       !$acc parallel loop collapse(2) present(U,V,W,nu_t,y,phi) reduction(max:local_conv,local_visc)
       Do k = 2, nzg-1
          Do j = 2, nyg-1
             ! y-grid is non-uniform: recompute inv_dy only once per j.
             inv_dy  = 1d0 / Max(y(j)-y(j-1), 1d-14)
             inv_dy2 = inv_dy * inv_dy
             Do i = 2, nxg-1
                ! Skip solid cells
                If ( phi(i,j,k) < 0d0 ) Cycle

                ! Convective: use centre-cell averaged velocities
                conv_ijk = Abs(0.5d0*(U(i,j,k)+U(i-1,j,k)))*inv_dx + &
                           Abs(0.5d0*(V(i,j,k)+V(i,j-1,k)))*inv_dy + &
                           Abs(0.5d0*(W(i,j,k)+W(i,j,k-1)))*inv_dz
                local_conv = Max(local_conv, conv_ijk)

                ! Viscous
                visc_ijk = (nu + nu_t(i,j,k)) * 2d0 * (inv_dx2 + inv_dy2 + inv_dz2)
                local_visc = Max(local_visc, visc_ijk)
             End Do
          End Do
       End Do
       !$acc end parallel loop
    Else
       !$acc parallel loop collapse(2) present(U,V,W,nu_t,y) reduction(max:local_conv,local_visc)
       Do k = 2, nzg-1
          Do j = 2, nyg-1
             inv_dy  = 1d0 / Max(y(j)-y(j-1), 1d-14)
             inv_dy2 = inv_dy * inv_dy
             Do i = 2, nxg-1
                conv_ijk = Abs(0.5d0*(U(i,j,k)+U(i-1,j,k)))*inv_dx + &
                           Abs(0.5d0*(V(i,j,k)+V(i,j-1,k)))*inv_dy + &
                           Abs(0.5d0*(W(i,j,k)+W(i,j,k-1)))*inv_dz
                local_conv = Max(local_conv, conv_ijk)

                visc_ijk = (nu + nu_t(i,j,k)) * 2d0 * (inv_dx2 + inv_dy2 + inv_dz2)
                local_visc = Max(local_visc, visc_ijk)
             End Do
          End Do
       End Do
       !$acc end parallel loop
    End If

    local_conv = local_conv * dt
    local_visc = local_visc * dt

    ! Single Allreduce for both quantities to halve MPI collective latency.
    local_buf(1) = local_conv
    local_buf(2) = local_visc
    Call MPI_Allreduce(local_buf, global_buf, 2, MPI_real8, MPI_MAX, MPI_COMM_WORLD, ierr)
    cfl_conv_out = global_buf(1)
    cfl_visc_out = global_buf(2)

  End Subroutine compute_cfl

  !> Volume-weighted (non-uniform y, IBM-solid-excluded) domain bulk streamwise velocity, for constant-mass-flux forcing
  Subroutine compute_bulk_velocity(Ub_out)

    Real(Int64), Intent(Out) :: Ub_out

    Real(Int64)    :: local_sum, local_wgt, dy_j
    Real(Int64)    :: local_buf(2), global_buf(2)
    Integer(Int32) :: i, j, k

    local_sum = 0d0
    local_wgt = 0d0

    If ( ibm_input_mode >= 1 ) Then
       !$acc parallel loop collapse(2) present(U,y,phi) reduction(+:local_sum,local_wgt)
       Do k = 2, nzg-1
          Do j = 2, nyg-1
             dy_j = y(j) - y(j-1)
             Do i = 2, nxg-1
                If ( phi(i,j,k) < 0d0 ) Cycle
                local_sum = local_sum + 0.5d0*(U(i,j,k)+U(i-1,j,k)) * dy_j
                local_wgt = local_wgt + dy_j
             End Do
          End Do
       End Do
       !$acc end parallel loop
    Else
       !$acc parallel loop collapse(2) present(U,y) reduction(+:local_sum,local_wgt)
       Do k = 2, nzg-1
          Do j = 2, nyg-1
             dy_j = y(j) - y(j-1)
             Do i = 2, nxg-1
                local_sum = local_sum + 0.5d0*(U(i,j,k)+U(i-1,j,k)) * dy_j
                local_wgt = local_wgt + dy_j
             End Do
          End Do
       End Do
       !$acc end parallel loop
    End If

    local_buf(1) = local_sum
    local_buf(2) = local_wgt
    Call MPI_Allreduce(local_buf, global_buf, 2, MPI_real8, MPI_SUM, MPI_COMM_WORLD, ierr)
    Ub_out = global_buf(1) / Max(global_buf(2), 1d-14)

  End Subroutine compute_bulk_velocity

  !  Open ibm_forces.csv on rank 0 and write header row.
  !  Call once after setup_ibm completes.
  Subroutine open_force_csv

    Logical :: file_exists
    If ( myid /= 0 ) Return
    Inquire(file='ibm_forces.csv', exist=file_exists)
    If ( file_exists ) Then
       ! Restart: append to existing file (no new header)
       Open(newunit=ibm_force_unit, file='ibm_forces.csv', status='old', &
            position='append', action='write')
    Else
       ! Fresh start: create file and write header
       Open(newunit=ibm_force_unit, file='ibm_forces.csv', status='new', action='write')
       Write(ibm_force_unit, '(A)') &
            'step,t,'  // &
            'Fx_ibm,Fy_ibm,Fz_ibm,'  // &
            'Fx_pres,Fy_pres,Fz_pres,'  // &
            'Fx_visc,Fy_visc,Fz_visc'
    End If

  End Subroutine open_force_csv

  !  Append one row to ibm_forces.csv (rank 0 only).
  Subroutine write_force_csv(Fx_ibm, Fy_ibm, Fz_ibm, &
                              Fx_pres, Fy_pres, Fz_pres, &
                              Fx_visc, Fy_visc, Fz_visc)

    Real(Int64), Intent(In) :: Fx_ibm,  Fy_ibm,  Fz_ibm
    Real(Int64), Intent(In) :: Fx_pres, Fy_pres, Fz_pres
    Real(Int64), Intent(In) :: Fx_visc, Fy_visc, Fz_visc

    If ( myid /= 0 ) Return
    Write(ibm_force_unit, '(I10,10(",",ES18.10))') &
         istep, t, &
         Fx_ibm,  Fy_ibm,  Fz_ibm, &
         Fx_pres, Fy_pres, Fz_pres, &
         Fx_visc, Fy_visc, Fz_visc

  End Subroutine write_force_csv

  !  Close ibm_forces.csv (rank 0 only).
  Subroutine close_force_csv

    If ( myid /= 0 ) Return
    If ( ibm_force_unit > 0 ) Close(ibm_force_unit)

  End Subroutine close_force_csv

End Module monitor
