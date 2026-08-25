Module mpi

  ! General Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64

  ! prevent implicit typing
  Implicit None

  Include 'mpif.h'

  ! declarations
  Integer(Int32) :: ierr, myid, nprocs, nslices_z
  Integer        :: istat ( MPI_STATUS_SIZE )

  ! global faces index range for each processor
  Integer(Int32), Dimension(:), Allocatable ::  k1_global,  k2_global

  ! global centers index range for each processor
  Integer(Int32), Dimension(:), Allocatable :: kg1_global, kg2_global

End Module mpi
