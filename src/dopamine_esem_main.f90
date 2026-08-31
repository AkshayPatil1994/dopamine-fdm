!> Standalone ESEM inflow precursor generator: precomputes an x-normal inflow-plane time series (mean + Cholesky-correlated fluctuation, as apply_inflow_bc_x would inject under inflow_type==1) and writes it in probe_output.f90's donor slice format, consumed unmodified by the existing inflow_type==2 (recycled precursor inflow) reader in sem.f90 -- amortizes raw_eddy_sum's per-eddy cost to a one-time build instead of every timestep.
!> Usage: mpirun -np N dopamine-ESEM [--out=BASENAME] [--samples=N] [--duration=T] [--ensemble-samples=N] [--ensemble-periods=P]; reads ./input_parameters (fixed filename, see read_input_parameters), the same namelist a live inflow_type=1/sem_use_esem=1 run would use (sem_profile_format=0 or 1 both supported). Only z is decomposed (y never is, matching the live solver and probe_output's own writer), so the donor file is decomposition-independent -- the consuming run can use a different nprocs/p_row/p_col.
Program dopamine_esem

  Use iso_fortran_env, Only : Int32, Int64, output_unit
  Use mpi
  Use global
  Use decomp
  Use genGridAndIC
  Use input_output
  Use synthetic_eddy_method

  Implicit None

  Character(300) :: out_base, arg
  Integer(Int32) :: n_samples, iarg, nargs, eqpos
  Real   (Int64) :: duration, dt_sample, t_s
  Logical :: duration_given, samples_given

  Integer(Int32) :: nz_local, j, k, isamp, u_bin, u_times, u_meta
  Real   (Int64), Allocatable :: lbuf(:,:,:), gbuf(:,:,:)
  Real   (Int64) :: up, vp, wp
  Character(300) :: fname

  Call Mpi_Init(ierr)
  Call Mpi_Comm_rank(MPI_COMM_WORLD, myid, ierr)
  Call Mpi_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

  ! --key=value only; the namelist path is fixed (read_input_parameters always opens ./input_parameters)
  out_base       = 'esem_precursor'
  duration_given = .False.
  samples_given  = .False.
  duration       = 0d0
  n_samples      = 0

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
     Case ('out');              out_base = Trim(arg(eqpos+1:))
     Case ('samples');          Read(arg(eqpos+1:),*) n_samples;          samples_given  = .True.
     Case ('duration');         Read(arg(eqpos+1:),*) duration;           duration_given = .True.
     Case ('ensemble-samples'); Read(arg(eqpos+1:),*) sem_ensemble_samples
     Case ('ensemble-periods'); Read(arg(eqpos+1:),*) sem_ensemble_periods
     Case Default
        If ( myid == 0 ) Write(*,'(A,A)') ' ERROR: unrecognised option: ', Trim(arg)
        Call Mpi_Abort(MPI_COMM_WORLD, 1, ierr)
     End Select
  End Do

  If ( myid == 0 ) Then
     Write(*,'(A)') ' dopamine-ESEM: standalone ESEM inflow precursor generator'
     Write(*,'(A,A)') '   donor basename = ', Trim(out_base)
  End If

  Call read_input_parameters

  If ( x_bc_type /= 1 ) Stop 'ERROR: dopamine-ESEM requires x_bc_type=1 in the namelist'
  If ( inflow_type /= 1 ) Stop 'ERROR: dopamine-ESEM requires inflow_type=1 (SEM) in the namelist -- ' // &
       'this generates the donor file a *separate* inflow_type=2 successor run will consume, ' // &
       'it does not itself read inflow_type=2 settings'
  If ( sem_use_esem /= 1 ) Stop 'ERROR: dopamine-ESEM requires sem_use_esem=1 (classical SEM has no ' // &
       'well-defined per-grid-point normalisation to precompute from, see init_inflow)'

  Call decomp_init_pencil
  Call decomp_build_xz_ranges
  If ( Any( (kg2_global-kg1_global-1) < 2 ) ) &
       Stop 'ERROR: each rank needs at least 2 interior z-cells -- use fewer ranks or a larger nz'

  ! grid: face points then cell centers (global)
  Allocate( x_global(nx_global), y_global(ny_global), z_global(nz_global) )
  Call generateGrid

  nxm_global = nx_global - 1
  nym_global = ny_global - 1
  nzm_global = nz_global - 1
  Allocate( xm_global(nxm_global), ym_global(nym_global), zm_global(nzm_global) )
  Do j = 1, nym_global
     ym_global(j) = 0.5d0*( y_global(j) + y_global(j+1) )
  End Do
  Do k = 1, nzm_global
     zm_global(k) = 0.5d0*( z_global(k) + z_global(k+1) )
  End Do

  ! this rank's local z-slab of the donor grid, same global-index convention write_slice_n uses (kg_g = kg1_global(myid)+ka-2)
  nz_local = kg2_global(myid) - kg1_global(myid) - 1

  ! substitute the module's (y,z,yg,zg) with the donor's shared cell-center grid before calling init_inflow*, so all three velocity components land on the ONE shared (n1,n2) grid the donor format needs (matching a real precursor's cc_val-interpolated probe slice) instead of each component's own staggered location; init_inflow's build_ensemble_normalisation call then builds ens_mean/ens_std directly at these points, keeping sem_fluctuation's (j,k) lookups self-consistent
  ny = nym_global;  nz = nz_local
  nyg = ny;          nzg = nz
  Allocate( y(ny), yg(ny) );  y = ym_global;  yg = ym_global
  Allocate( z(nz), zg(nz) );  z = zm_global( kg1_global(myid) : kg2_global(myid)-2 )
  zg = z

  Call init_inflow_profile   ! dispatches read_mean_profile / read_mean_profile_TI on sem_profile_format
  Call init_inflow           ! sigma file, wall taper, Uconv_sem, length-scale/n_eddies auto-tune, place_eddies, build_ensemble_normalisation (now at this rank's cell-center slab)

  If ( .Not. duration_given ) duration = Real(sem_ensemble_periods,8) * Maxval(eddy_Tperiod)
  ! ~8 samples per the fastest eddy's transit period, capped so a very fine near-wall length scale can't blow the table up unboundedly
  If ( .Not. samples_given ) n_samples = Min( 200000, Max( 50, Ceiling( 8d0 * duration / Minval(eddy_Tperiod) ) ) )
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

  Allocate( lbuf(3,nym_global,nzm_global), gbuf(3,nym_global,nzm_global) )

  Do isamp = 1, n_samples
     t_s = Real(isamp-1,8) * dt_sample

     lbuf = 0d0
     Do k = 1, nz_local
        Do j = 1, nym_global
           Call sem_fluctuation( 1, j, k, y(j), z(k), t_s, up, vp, wp )
           lbuf(1, j, kg1_global(myid)-1+k) = mean_profile_U( y(j) ) + up
           Call sem_fluctuation( 2, j, k, y(j), z(k), t_s, up, vp, wp )
           lbuf(2, j, kg1_global(myid)-1+k) = vp
           Call sem_fluctuation( 3, j, k, y(j), z(k), t_s, up, vp, wp )
           lbuf(3, j, kg1_global(myid)-1+k) = wp
        End Do
     End Do

     Call Mpi_Reduce( lbuf, gbuf, 3*nym_global*nzm_global, MPI_real8, MPI_SUM, 0, MPI_COMM_WORLD, ierr )

     If ( myid == 0 ) Then
        Write(u_bin) gbuf
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
