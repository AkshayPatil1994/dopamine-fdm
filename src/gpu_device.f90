!> Rank-to-GPU binding for multi-GPU runs: each MPI rank on a node picks a
!> distinct physical GPU by node-local rank, modulo the GPU count on that
!> node, before any device memory is allocated or any cuFFT/cuSPARSE plan is
!> created. Compiled only in ENABLE_GPU builds (CMakeLists.txt sets
!> -DGPU_POISSON); see initialization.f90's call right after MPI_Init.
Module gpu_device

  Use iso_fortran_env, Only : Int32
  Use openacc
  Use mpi

  Implicit None

Contains

  !> Splits MPI_COMM_WORLD by shared-memory locality (one sub-communicator per
  !> node) and assigns device (local_rank mod ndevices_on_node) to this rank.
  !> acc_set_device_num binds both the OpenACC runtime and the underlying CUDA
  !> context that cuFFT/cuSPARSE calls pick up (Use cufft/cusparse in
  !> poisson_gpu.f90) -- no separate cudaSetDevice call is needed.
  Subroutine assign_gpu_device

    Integer :: local_comm, local_rank, local_nprocs, ndevices, device_id, ierr_local

    Call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, &
                              MPI_INFO_NULL, local_comm, ierr_local)
    Call MPI_Comm_rank(local_comm, local_rank,  ierr_local)
    Call MPI_Comm_size(local_comm, local_nprocs, ierr_local)

    ndevices = acc_get_num_devices(acc_device_nvidia)
    If ( ndevices < 1 ) Stop 'ERROR: ENABLE_GPU build found no NVIDIA GPU visible to this rank ' // &
                             '(check CUDA_VISIBLE_DEVICES and that the node has a working driver)'

    device_id = Mod(local_rank, ndevices)
    Call acc_set_device_num(device_id, acc_device_nvidia)

    If ( local_nprocs > ndevices ) Then
       Write(*,'(A,I0,A,I0,A,I0,A)') ' WARNING: rank ', myid, ': ', local_nprocs, &
            ' MPI ranks share this node but only ', ndevices, &
            ' GPU(s) are visible -- multiple ranks will oversubscribe the same device'
    End If

    Write(*,'(A,I0,A,I0,A,I0,A,I0,A)') ' GPU bind: global rank ', myid, ' (node-local rank ', local_rank, &
         ' of ', local_nprocs, ') -> device ', device_id

    Call MPI_Comm_free(local_comm, ierr_local)

  End Subroutine assign_gpu_device

End Module gpu_device
