Module mpi

  ! General Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64

  ! prevent implicit typing
  Implicit None

  Include 'mpif.h'

  ! declarations
  Integer(Int32) :: ierr, myid, nprocs
  Integer        :: istat ( MPI_STATUS_SIZE )

  ! global faces index range for each processor (z-split, via p_col)
  Integer(Int32), Dimension(:), Allocatable ::  k1_global,  k2_global

  ! global centers index range for each processor (z-split, via p_col)
  Integer(Int32), Dimension(:), Allocatable :: kg1_global, kg2_global

  ! global faces index range for each processor (x-split, via p_row)
  Integer(Int32), Dimension(:), Allocatable ::  i1_global,  i2_global

  ! global centers index range for each processor (x-split, via p_row)
  Integer(Int32), Dimension(:), Allocatable :: ig1_global, ig2_global

End Module mpi
