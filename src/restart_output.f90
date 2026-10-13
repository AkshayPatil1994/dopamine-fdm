!> Restart handling of the time-series output files: on restart=1 the run continues from the snapshot, which can be older than
!  the point where the previous run stopped, so the records written after the snapshot are dropped before new ones are appended.
Module restart_output

  Use iso_fortran_env, Only : Int32, Int64
  Use global
  Use mpi

  Implicit None
  Private
  Public :: trim_restart_text_files, trim_restart_file, truncate_file, file_exists_nonempty

Contains

  Logical Function file_exists_nonempty(fname)

    Character(*), Intent(In) :: fname
    Logical :: ex
    Integer(Int64) :: sz

    Inquire(file=fname, exist=ex, size=sz)
    file_exists_nonempty = ex .And. sz > 0

  End Function file_exists_nonempty


  !> Cut a file to nbytes (rank 0)
  Subroutine truncate_file(fname, nbytes)

    Character(*), Intent(In) :: fname
    Integer(Int64), Intent(In) :: nbytes
    Character(40) :: sz
    Integer(Int32) :: stat

    Write(sz,'(I0)') nbytes
    Call execute_command_line('truncate -s ' // Trim(sz) // ' ' // Trim(fname), wait=.True., exitstat=stat)
    If ( stat /= 0 ) Then
       Write(*,'(A,A)') ' ERROR: cannot truncate ', Trim(fname)
       Call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    End If

  End Subroutine truncate_file


  !> Drop the rows of a text file whose time (field tcol, blanks or commas as separators) is later than tmax; comment lines (#) and
  !  the header line of a csv (first line, not numeric) stay
  Subroutine trim_rows_after(fname, tcol, tmax)

    Character(*), Intent(In) :: fname
    Integer(Int32), Intent(In) :: tcol
    Real(Int64), Intent(In) :: tmax

    Character(2000) :: line, work
    Character(len=Len_Trim(fname)+4) :: tmp
    Integer(Int32) :: fi, fo, ios, ios2, ic, nkept, ndrop
    Real(Int64) :: tr
    Logical :: keep

    tmp = Trim(fname) // '.tmp'
    Open(newunit=fi, file=fname, status='old', action='read')
    Open(newunit=fo, file=tmp, status='replace', action='write')
    nkept = 0;  ndrop = 0
    Do
       Read(fi,'(A)', iostat=ios) line
       If ( ios /= 0 ) Exit
       keep = .True.
       If ( line(1:1) /= '#' ) Then
          work = line
          Do ic = 1, Len(work)
             If ( work(ic:ic) == ',' ) work(ic:ic) = ' '
          End Do
          Read(work, *, iostat=ios2) ( tr, ic = 1, tcol )   ! the last value read is the time column
          If ( ios2 == 0 ) keep = ( tr <= tmax )
       End If
       If ( keep ) Then
          Write(fo,'(A)') Trim(line)
          nkept = nkept + 1
       Else
          ndrop = ndrop + 1
       End If
    End Do
    Close(fi)
    Close(fo)
    Call execute_command_line('mv ' // Trim(tmp) // ' ' // Trim(fname), wait=.True.)
    If ( ndrop > 0 ) Write(*,'(A,A,A,I0,A)') ' Restart: ', Trim(fname), ' dropped ', ndrop, &
         ' rows written after the restart snapshot'

  End Subroutine trim_rows_after


  !> Restart (rank 0): drop the rows after the restart time from a time-series file, which must not be open yet
  Subroutine trim_restart_file(fname, tcol)

    Character(*), Intent(In) :: fname
    Integer(Int32), Intent(In) :: tcol

    If ( restart /= 1 .Or. myid /= 0 ) Return
    If ( file_exists_nonempty(fname) ) Call trim_rows_after(fname, tcol, t + 1d-6*Max(1d0, Abs(t)))

  End Subroutine trim_restart_file


  Subroutine trim_restart_text_files

    Call trim_restart_file('vof_diag.dat', 2)
    Call trim_restart_file('vof_gauges.dat', 1)

  End Subroutine trim_restart_text_files

End Module restart_output
