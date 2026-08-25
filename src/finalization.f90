! Module to finalize FFT and MPI
Module finalization
  
  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use mpi
  Use monitor,   Only : close_force_csv
  Use profiler,  Only : profiler_print_summary
  
  ! prevent implicit typing
  Implicit None
  
Contains

  ! Finalize FFTW plans and MPI
  Subroutine finalize
  
    Real(Int64) :: total_wall_time

    ! finalize FFTW
    Call dfftw_destroy_plan(plan_d)
    Call dfftw_destroy_plan(plan_i)

    Call Mpi_barrier(MPI_COMM_WORLD, ierr)
    total_wall_time = MPI_WTIME() - time_wall_start
    
    If ( myid==0 ) Then
       Write(*,'(A)') ' '
       Write(*,'(A)') ' +====================================================================+'
       Write(*,'(A)') '  RUN COMPLETE'
       Write(*,'(A)') ' +====================================================================+'
       Write(*,'(A,I0)')  '  Steps completed    :  ', istep
       Write(*,'(A,F12.6)') '  Physical time      :  ', t
       Write(*,'(A,F10.2,A)') '  Total wall time    :  ', total_wall_time, ' s'
       Write(*,'(A,F10.4,A)') '  Avg wall/step      :  ', &
            total_wall_time / Max(Real(istep, 8), 1d0), ' s'
       Write(*,'(A)') ' +====================================================================+'
       Write(*,'(A)') ' '
    End If

    Call close_force_csv

    Call profiler_print_summary

    ! finalized MPI
    Call Mpi_finalize(ierr)
  End Subroutine finalize
  
End Module finalization
