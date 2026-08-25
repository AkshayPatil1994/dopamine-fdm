!> Per-stage MPI_Wtime() profiler
Module profiler

  Use iso_fortran_env, Only : Int32, Int64
  Use mpi

  Implicit None

  Integer(Int32), Parameter :: PROF_SGS            = 1
  Integer(Int32), Parameter :: PROF_WALLMODEL      = 2
  Integer(Int32), Parameter :: PROF_RHS            = 3
  Integer(Int32), Parameter :: PROF_IBM            = 4
  Integer(Int32), Parameter :: PROF_BC             = 5
  Integer(Int32), Parameter :: PROF_POISSON_FFT    = 6
  Integer(Int32), Parameter :: PROF_POISSON_TRIDIAG = 7
  Integer(Int32), Parameter :: PROF_PROJECTION     = 8
  Integer(Int32), Parameter :: PROF_SCALAR         = 9
  Integer(Int32), Parameter :: PROF_CFL            = 10
  Integer(Int32), Parameter :: PROF_RK_UPDATE      = 11
  Integer(Int32), Parameter :: PROF_CMFR           = 12
  Integer(Int32), Parameter :: PROF_N              = 12

  Character(16), Dimension(PROF_N), Parameter :: prof_names = [ Character(16) :: &
       'SGS model', 'Wall model', 'RHS (equations)', 'IBM', &
       'Boundary cond.', 'Poisson FFT', 'Poisson tridiag', 'Projection', &
       'Scalar transport', 'CFL check', 'RK update+sync', 'Const. flow rate' ]

  Real(Int64) :: prof_elapsed(PROF_N) = 0d0
  Real(Int64) :: prof_t0(PROF_N)      = 0d0

Contains

  !> Start the timer for stage `tag` (one of the PROF_* constants above)
  Subroutine profiler_start(tag)
    Integer(Int32), Intent(In) :: tag
    prof_t0(tag) = MPI_WTIME()
  End Subroutine profiler_start

  !> Accumulate elapsed time for stage `tag` since the matching profiler_start
  Subroutine profiler_stop(tag)
    Integer(Int32), Intent(In) :: tag
    prof_elapsed(tag) = prof_elapsed(tag) + ( MPI_WTIME() - prof_t0(tag) )
  End Subroutine profiler_stop

  !> Reduce per-rank timings to rank 0 and print a summary table (call once, at shutdown)
  Subroutine profiler_print_summary

    Real(Int64) :: prof_sum(PROF_N), prof_max(PROF_N), total_sum
    Integer(Int32) :: i

    Call MPI_Reduce(prof_elapsed, prof_sum, PROF_N, MPI_real8, MPI_sum, 0, MPI_COMM_WORLD, ierr)
    Call MPI_Reduce(prof_elapsed, prof_max, PROF_N, MPI_real8, MPI_max, 0, MPI_COMM_WORLD, ierr)

    If ( myid == 0 ) Then
       total_sum = Sum(prof_sum)
       Write(*,'(A)') ' '
       Write(*,'(A)') ' +====================================================================+'
       Write(*,'(A)') '  PROFILER SUMMARY (per-rank average, and slowest-rank max)'
       Write(*,'(A)') ' +====================================================================+'
       Write(*,'(A)') '  Stage               avg/rank [s]   max-rank [s]   % of tracked time'
       Write(*,'(A)') '  --------------------------------------------------------------------'
       Do i = 1, PROF_N
          If ( prof_sum(i) > 0d0 ) Write(*,'(2X,A16,3X,F10.4,4X,F10.4,7X,F6.2)') &
               prof_names(i), prof_sum(i)/Real(nprocs,8), prof_max(i), &
               100d0*prof_sum(i)/Max(total_sum,1d-14)
       End Do
       Write(*,'(A)') '  --------------------------------------------------------------------'
       Write(*,'(A,F10.4,A)') '  Total tracked (avg/rank)  :  ', total_sum/Real(nprocs,8), ' s'
       Write(*,'(A)') ' +====================================================================+'
       Write(*,'(A)') ' '
    End If

  End Subroutine profiler_print_summary

End Module profiler
