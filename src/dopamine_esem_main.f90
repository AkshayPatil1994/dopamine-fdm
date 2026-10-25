!> Standalone ESEM inflow precursor generator (dopamine-ESEM): precomputes an x-normal inflow-plane time
!  series (mean + Cholesky-correlated fluctuation, as apply_inflow_bc_x injects under inflow_type==1) in
!  probe_output.f90's donor slice format, read unchanged by the inflow_type==2 reader of sem.f90, so
!  raw_eddy_sum's cost is paid once. Usage: mpirun -np N dopamine-ESEM --input=PATH [--out=BASENAME]
!  [--samples=N | --dt=T] [--duration=T] [--ensemble-samples=N] [--ensemble-periods=P]. --input is required
!  (the namelist of an inflow_type=1/sem_use_esem=1 run) so a leftover ./input_parameters is never used. Only
!  z is decomposed: the donor file is independent of the rank layout. --dt should equal the consuming run's dt
!  (or an integer fraction): the default spacing (~1/8 of the fastest eddy's transit period) is much finer
!  than a CFL-limited dt, so the consuming run would sample near-uncorrelated points and the projection
!  removes that incoherent forcing (Reynolds stresses collapse near the inlet); matching it keeps the time
!  correlation.
Program dopamine_esem

  Use iso_fortran_env, Only : Int32, Int64, output_unit
  Use mpi
  Use global
  Use decomp
  Use genGridAndIC
  Use input_output
  Use synthetic_eddy_method

  Implicit None

  Character(300) :: out_base, arg, input_path
  Integer(Int32) :: n_samples, iarg, nargs, eqpos
  Real   (Int64) :: duration, dt_sample, dt_sample_in, t_s
  Logical :: duration_given, samples_given, input_given, dt_given

  Integer(Int32) :: nz_local, j, k, isamp, u_bin, u_times, u_meta
  Integer(Int32) :: ja, ka, jg, jf, kg_g
  Real   (Int64), Allocatable :: lbuf_U(:,:), gbuf_U(:,:)
  Real   (Int64), Allocatable :: lbuf_V(:,:), gbuf_V(:,:)
  Real   (Int64), Allocatable :: lbuf_W(:,:), gbuf_W(:,:)
  Real   (Int64) :: up, vp, wp
  Character(300) :: fname

  Call Mpi_Init(ierr)
  Call Mpi_Comm_rank(MPI_COMM_WORLD, myid, ierr)
  Call Mpi_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

  ! --key=value only; --input is required (no fallback, see the header)
  out_base       = 'esem_precursor'
  duration_given = .False.
  samples_given  = .False.
  input_given    = .False.
  dt_given       = .False.
  duration       = 0d0
  n_samples      = 0
  dt_sample_in   = 0d0

  nargs = Command_Argument_Count()
  Do iarg = 1, nargs
     Call Get_Command_Argument(iarg, arg)
     If ( arg(1:2) /= '--' ) Then
        If ( myid == 0 ) Write(*,'(A,A)') ' ERROR: unrecognised argument (expected --key=value): ', Trim(arg)
        Call Mpi_Abort(MPI_COMM_WORLD, 1, ierr)
     End If
     eqpos = Index(arg,'=')
     If ( eqpos < 1 ) Then
        If ( myid == 0 ) Write(*,'(A,A)') ' ERROR: malformed option (expected --key=value): ', Trim(arg)
        Call Mpi_Abort(MPI_COMM_WORLD, 1, ierr)
     End If
     Select Case ( arg(3:eqpos-1) )
     Case ('input');            input_path = Trim(arg(eqpos+1:));         input_given    = .True.
     Case ('out');              out_base = Trim(arg(eqpos+1:))
     Case ('samples');          Read(arg(eqpos+1:),*) n_samples;          samples_given  = .True.
     Case ('duration');         Read(arg(eqpos+1:),*) duration;           duration_given = .True.
     Case ('dt');               Read(arg(eqpos+1:),*) dt_sample_in;       dt_given       = .True.
     Case ('ensemble-samples'); Read(arg(eqpos+1:),*) sem_ensemble_samples
     Case ('ensemble-periods'); Read(arg(eqpos+1:),*) sem_ensemble_periods
     Case Default
        If ( myid == 0 ) Write(*,'(A,A)') ' ERROR: unrecognised option: ', Trim(arg)
        Call Mpi_Abort(MPI_COMM_WORLD, 1, ierr)
     End Select
  End Do

  If ( .Not. input_given ) Then
     If ( myid == 0 ) Write(*,'(A)') &
          ' ERROR: no input_parameters file specified -- pass --input=PATH ' // &
          '(no default is read, so this never silently picks up whatever ./input_parameters happens to hold)'
     Call Mpi_Abort(MPI_COMM_WORLD, 1, ierr)
  End If

  If ( dt_given .And. samples_given ) Then
     If ( myid == 0 ) Write(*,'(A)') &
          ' ERROR: --dt and --samples both set the donor sample spacing -- pass only one'
     Call Mpi_Abort(MPI_COMM_WORLD, 1, ierr)
  End If
  If ( dt_given .And. dt_sample_in <= 0d0 ) Call abort_input( 'ERROR: --dt must be > 0' )

  If ( myid == 0 ) Then
     Write(*,'(A)') ' dopamine-ESEM: standalone ESEM inflow precursor generator'
     Write(*,'(A,A)') '   input_parameters = ', Trim(input_path)
     Write(*,'(A,A)') '   donor basename   = ', Trim(out_base)
  End If

  Call read_input_parameters(Trim(input_path))

  If ( x_bc_type /= 1 ) Call abort_input( 'ERROR: dopamine-ESEM requires x_bc_type=1 in the namelist' )
  If ( inflow_type /= 1 ) Call abort_input( 'ERROR: dopamine-ESEM requires inflow_type=1 (SEM) in the namelist -- ' // &
       'this generates the donor file a *separate* inflow_type=2 successor run will consume, ' // &
       'it does not itself read inflow_type=2 settings' )
  If ( sem_use_esem /= 1 ) Call abort_input( 'ERROR: dopamine-ESEM requires sem_use_esem=1 (classical SEM has no ' // &
       'well-defined per-grid-point normalisation to precompute from, see init_inflow)' )

  Call decomp_init_pencil
  Call decomp_build_xz_ranges
  If ( Any( (kg2_global-kg1_global-1) < 2 ) ) &
       Call abort_input( 'ERROR: each rank needs at least 2 interior z-cells -- use fewer ranks or a larger nz' )

  ! grid: face points then cell centers (global)
  Allocate( x_global(nx_global), y_global(ny_global), z_global(nz_global) )

  ! generateGrid's grid.out/geometry.out writer loops over nxm_global/nym_global/nzm_global, so these are set
  ! before the call (as initialize() does) or it silently writes an empty grid.out.
  nxm_global = nx_global - 1
  nym_global = ny_global - 1
  nzm_global = nz_global - 1
  Call generateGrid

  Allocate( xm_global(nxm_global), ym_global(nym_global), zm_global(nzm_global) )
  Do j = 1, nym_global
     ym_global(j) = 0.5d0*( y_global(j) + y_global(j+1) )
  End Do
  Do k = 1, nzm_global
     zm_global(k) = 0.5d0*( z_global(k) + z_global(k+1) )
  End Do

  ! this rank's local z-slab of the donor grid (interior cell/face count); global index kg_g =
  ! kg1_global(myid)+ka-2 as in write_slice_n
  nz_local = kg2_global(myid) - kg1_global(myid) - 1

  ! Per-component staggered grids (same ghost extrapolation as initialization.f90) so ens_mean_*/ens_std_*
  ! (sized per component) are evaluated where apply_inflow_bc_x's live SEM evaluates them; avoids a half-cell
  ! shift in sem.f90's recycle_value (inflow_type=2).
  Allocate( yg_global(nym_global+2), zg_global(nzm_global+2) )
  yg_global(2:nym_global+1) = ym_global
  yg_global(1)              = ym_global(1)          - 2d0*(ym_global(1)-y_global(1))
  yg_global(nym_global+2)   = ym_global(nym_global) + 2d0*(y_global(ny_global)-ym_global(nym_global))
  zg_global(2:nzm_global+1) = zm_global
  zg_global(1)              = zm_global(1)          - 2d0*(zm_global(1)-z_global(1))
  zg_global(nzm_global+2)   = zm_global(nzm_global) + 2d0*(z_global(nz_global)-zm_global(nzm_global))

  ny  = ny_global                                 ! V's native (y-face, no ghost) count
  nyg = nym_global + 2                            ! U/W's native (y ghost cell-centre) count
  nzg = kg2_global(myid) - kg1_global(myid) + 1   ! U/V's native LOCAL (z ghost cell-centre) count
  ! W's z count equals nzg (both use interior ka=2..nzg-1); differs slightly from the true local face count on
  ! the last z-rank
  nz  = nzg

  Allocate( y(ny) );   y  = y_global
  Allocate( yg(nyg) ); yg = yg_global
  Allocate( zg(nzg) ); zg = zg_global( kg1_global(myid) : kg2_global(myid) )
  Allocate( z(nz) )
  Do k = 1, nz
     ! z(ka) = left z-face of the cell centred at zg(ka) (global index kg1_global(myid)+ka-2), W's staggered
     ! location; ghost slots ka=1/nz are never queried.
     z(k) = z_global( Max(1, Min(nzm_global, kg1_global(myid)+k-2)) )
  End Do

  Call init_inflow_profile   ! dispatches read_mean_profile / read_mean_profile_TI on sem_profile_format
  ! sigma file, wall taper, Uconv_sem, length-scale/n_eddies auto-tune, place_eddies, ensemble normalisation
  ! at this rank's slab
  Call init_inflow

  If ( .Not. duration_given ) duration = Real(sem_ensemble_periods,8) * Maxval(eddy_Tperiod)
  If ( samples_given .And. n_samples < 2 ) Then
     If ( myid == 0 ) Write(*,'(A,I0,A)') ' ERROR: --samples=', n_samples, ' must be >= 2'
     Call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  End If
  If ( dt_given ) Then
     ! pin the sample spacing to (a submultiple of) the consuming run's dt instead of the fastest-eddy default
     ! below, so donor frames are the instants recycle_value will query (see --dt in the header)
     n_samples = Max( 2, Ceiling(duration/dt_sample_in) + 1 )
  Else If ( .Not. samples_given ) Then
     ! ~8 samples per the fastest eddy's transit period, capped so a very fine near-wall length scale cannot
     ! blow up the table
     n_samples = Min( 200000, Max( 50, Ceiling( 8d0 * duration / Minval(eddy_Tperiod) ) ) )
  End If
  dt_sample = duration / Real(n_samples-1,8)

  If ( myid == 0 ) Then
     Write(*,'(A,E12.4,A,I8,A,E12.4)') ' Donor timeline: duration = ', duration, &
          ',  samples = ', n_samples, ',  dt_sample = ', dt_sample
     Write(*,'(A,I0,A,I0,A)') ' Donor grid: ', nym_global, ' x ', nzm_global, ' (ny x nz cell centers)'
  End If

  ! rank 0 only; mirrors probe_output's write_slice_n/init_probes format exactly (ncomp=3, dir=x, comps=UVW)
  If ( myid == 0 ) Then
     Write(fname,'(A,A)') Trim(out_base), '.bin'
     Open(newunit=u_bin, file=Trim(fname), access='stream', form='unformatted', status='replace', action='write')
     Write(fname,'(A,A)') Trim(out_base), '_times.bin'
     Open(newunit=u_times, file=Trim(fname), access='stream', form='unformatted', status='replace', action='write')
  End If

  ! U,W: native (yg,zg) cell-centre grid (nym_global x nzm_global). V: native (y,zg) grid (ny_global x
  ! nzm_global) in its own block written after U's, so the consumer needs no half-cell shift.
  Allocate( lbuf_U(nym_global,nzm_global), gbuf_U(nym_global,nzm_global) )
  Allocate( lbuf_V(ny_global, nzm_global), gbuf_V(ny_global, nzm_global) )
  Allocate( lbuf_W(nym_global,nzm_global), gbuf_W(nym_global,nzm_global) )

  Do isamp = 1, n_samples
     t_s = Real(isamp-1,8) * dt_sample

     lbuf_U = 0d0;  lbuf_V = 0d0;  lbuf_W = 0d0

     ! Ranks sharing col = Mod(myid,p_col) own the same z-slab; only row==0 ranks compute, otherwise the
     ! Mpi_Reduce(MPI_SUM) below would inflate the field by p_row (the rest leave lbuf_* at zero).
     If ( myid/p_col == 0 ) Then

        ! U: cell-centred in y and z
        Do ka = 2, nzg-1
           kg_g = kg1_global(myid) + ka - 2
           Do ja = 2, nyg-1
              jg = ja - 1
              Call sem_fluctuation( 1, ja, ka, yg(ja), zg(ka), t_s, up, vp, wp, .False. )
              lbuf_U(jg, kg_g) = mean_profile_U( yg(ja) ) + up
           End Do
        End Do

        ! V: y-face (all ny_global points, no ghost), z cell-centred
        Do ka = 2, nzg-1
           kg_g = kg1_global(myid) + ka - 2
           Do jf = 1, ny
              Call sem_fluctuation( 2, jf, ka, y(jf), zg(ka), t_s, up, vp, wp, .False. )
              lbuf_V(jf, kg_g) = vp
           End Do
        End Do

        ! W: y cell-centred, z-face (left face of each cell, see z(:) above)
        Do ka = 2, nz-1
           kg_g = kg1_global(myid) + ka - 2
           Do ja = 2, nyg-1
              jg = ja - 1
              Call sem_fluctuation( 3, ja, ka, yg(ja), z(ka), t_s, up, vp, wp, .False. )
              lbuf_W(jg, kg_g) = wp
           End Do
        End Do

     End If

     Call Mpi_Reduce( lbuf_U, gbuf_U, nym_global*nzm_global, MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
     Call Mpi_Reduce( lbuf_V, gbuf_V, ny_global*nzm_global,  MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
     Call Mpi_Reduce( lbuf_W, gbuf_W, nym_global*nzm_global, MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, ierr )

     If ( myid == 0 ) Then
        Write(u_bin) gbuf_U
        Write(u_bin) gbuf_V
        Write(u_bin) gbuf_W
        Write(u_times) t_s
        Flush(u_bin); Flush(u_times)
        Call print_progress(isamp, n_samples)
     End If
  End Do

  If ( myid == 0 ) Then
     Close(u_bin); Close(u_times)

     Write(fname,'(A,A)') Trim(out_base), '_meta.txt'
     Open(newunit=u_meta, file=Trim(fname), form='formatted', status='replace')
     Write(u_meta,'(A,I0)')    'ncomp  = ', 3
     Write(u_meta,'(A,I0)')    'n1     = ', nym_global
     Write(u_meta,'(A,I0)')    'n2     = ', nzm_global
     Write(u_meta,'(A,I0)')    'n1_V   = ', ny_global
     Write(u_meta,'(A,A)')     'dir    = ', 'x'
     Write(u_meta,'(A,F16.8)') 'pos   = ', x_global(1)
     Write(u_meta,'(A,A)')     'comps  = ', 'UVW'
     Write(u_meta,'(A,I0)')    'nsnaps = ', n_samples
     Write(u_meta,'(A,A)')     'times  = ', Trim(out_base)//'_times.bin'
     Close(u_meta)

     Write(*,'(A)') ' '
     Write(*,'(A,A,A)') ' Done. To consume this donor, set inflow_type=2, inflow_recycle_file=''', &
          Trim(out_base), ''' in the successor run''s &INFLOW (add inflow_recycle_loop=1, ' // &
          'inflow_recycle_shift_z=1 to avoid spurious periodicity when the donor timeline wraps).'
  End If

  Call Mpi_Finalize(ierr)

Contains

  !> Rank-0-only progress bar: carriage-return refresh, no trailing newline until 100%
  Subroutine print_progress(i, n)

    Integer(Int32), Intent(In) :: i, n
    Integer(Int32), Parameter :: barw = 30
    Integer(Int32) :: filled
    Character(barw) :: bar

    filled = Int( Real(barw,8) * Real(i,8) / Real(n,8) )
    bar = Repeat('#', filled) // Repeat('-', barw-filled)

    Write(output_unit,'(A,A,A,A,I4,A,I8,A,I8,A)',advance='no') Achar(13), &
         ' [', bar, '] ', Int(100d0*Real(i,8)/Real(n,8)), '%  (', i, '/', n, ' samples)'
    Flush(output_unit)
    If ( i == n ) Write(output_unit,'(A)') ' '

  End Subroutine print_progress

End Program dopamine_esem
