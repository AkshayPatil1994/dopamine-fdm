!> Rank-to-GPU binding for multi-GPU runs: each MPI rank on a node picks a
!> distinct physical GPU by node-local rank, modulo the GPU count on that
!> node, before any device memory is allocated or any cuFFT/cuSPARSE plan is
!> created. The binding is made BEFORE MPI_Init when the launcher exports the
!> node-local rank (Open MPI/MVAPICH2/Slurm/PMI): a CUDA-aware MPI otherwise
!> initialises its CUDA context on device 0 for every rank first, wasting
!> memory there and breaking peer/IPC transfers between ranks. Compiled only in ENABLE_GPU builds (CMakeLists.txt sets
!> -DGPU_POISSON); see initialization.f90's call right after MPI_Init.
Module gpu_device

  Use iso_fortran_env, Only : Int32
  Use openacc
  Use mpi

  Implicit None

  Logical :: bound_before_mpi = .False.

Contains

  !> Bind this rank to a GPU from the launcher's node-local-rank environment variable, before MPI_Init; a no-op (bound_before_mpi stays false) if none is set, in which case assign_gpu_device binds after MPI_Init as a fallback
  Subroutine assign_gpu_device_pre_mpi

    Character(32) :: val
    Character(24), Parameter :: names(4) = [ Character(24) :: 'OMPI_COMM_WORLD_LOCAL_RANK', &
         'MV2_COMM_WORLD_LOCAL_RANK', 'MPI_LOCALRANKID', 'SLURM_LOCALID' ]
    Integer :: i, stat, length, local_rank, ndevices

    Do i = 1, Size(names)
       Call Get_Environment_Variable( Trim(names(i)), val, length, stat )
       If ( stat == 0 .And. length > 0 ) Then
          Read(val,*,iostat=stat) local_rank
          If ( stat /= 0 ) Cycle
          ndevices = acc_get_num_devices(acc_device_nvidia)
          If ( ndevices < 1 ) Return
          Call acc_set_device_num(Mod(local_rank, ndevices), acc_device_nvidia)
          bound_before_mpi = .True.
          Return
       End If
    End Do

  End Subroutine assign_gpu_device_pre_mpi

  !> Splits MPI_COMM_WORLD by shared-memory locality (one sub-communicator per
  !> node) and assigns device (local_rank mod ndevices_on_node) to this rank.
  !> acc_set_device_num binds both the OpenACC runtime and the underlying CUDA
  !> context that cuFFT/cuSPARSE calls pick up (Use cufft/cusparse in
  !> poisson_gpu.f90) -- no separate cudaSetDevice call is needed.
  Subroutine assign_gpu_device

    Integer :: local_comm, local_rank, local_nprocs, ndevices, device_id, ierr_local, nprocs_w
    Character(16) :: cuda_env
    Integer :: env_len, env_stat

    Call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, &
                              MPI_INFO_NULL, local_comm, ierr_local)
    Call MPI_Comm_rank(local_comm, local_rank,  ierr_local)
    Call MPI_Comm_size(local_comm, local_nprocs, ierr_local)

    ndevices = acc_get_num_devices(acc_device_nvidia)
    If ( ndevices < 1 ) Then
       Write(*,'(A)') 'ERROR: ENABLE_GPU build found no NVIDIA GPU visible to this rank (check CUDA_VISIBLE_DEVICES and driver)'
       Call MPI_Abort(MPI_COMM_WORLD, 1, ierr_local)
    End If

    device_id = Mod(local_rank, ndevices)
    ! (re)binding after MPI_Init is only needed when the launcher gave no local-rank variable; when
    ! bound_before_mpi is set, the same device is already active and MPI's context lives on it
    If ( .Not. bound_before_mpi ) Call acc_set_device_num(device_id, acc_device_nvidia)

    If ( local_nprocs > ndevices ) Then
       Write(*,'(A,I0,A,I0,A,I0,A)') ' WARNING: rank ', myid, ': ', local_nprocs, &
            ' MPI ranks share this node but only ', ndevices, &
            ' GPU(s) are visible -- multiple ranks will oversubscribe the same device'
    End If

    Write(*,'(A,I0,A,I0,A,I0,A,I0,A)') ' GPU bind: global rank ', myid, ' (node-local rank ', local_rank, &
         ' of ', local_nprocs, ') -> device ', device_id

    Call MPI_Comm_free(local_comm, ierr_local)

    ! device pointers are handed to MPI (host_data use_device in the halo/transposes): without CUDA-aware MPI this crashes or
    ! silently corrupts. No portable Fortran query exists, so warn when no launcher/MCA hint of CUDA support is visible.
    Call MPI_Comm_size(MPI_COMM_WORLD, nprocs_w, ierr_local)
    If ( nprocs_w > 1 .And. myid == 0 ) Then
       cuda_env = ''
       Call Get_environment_variable('OMPI_MCA_opal_cuda_support', cuda_env, env_len, env_stat)
       If ( env_stat /= 0 ) Call Get_environment_variable('MV2_USE_CUDA', cuda_env, env_len, env_stat)
       If ( env_stat /= 0 ) Call Get_environment_variable('MPICH_GPU_SUPPORT_ENABLED', cuda_env, env_len, env_stat)
       If ( env_stat /= 0 ) Call Get_environment_variable('UCX_TLS', cuda_env, env_len, env_stat)
       If ( env_stat /= 0 ) Write(*,'(A)') ' NOTE: multi-rank GPU run; make sure MPI is CUDA-aware (e.g. Open MPI built with ' // &
            'CUDA, or OMPI_MCA_opal_cuda_support=true), since device buffers are passed to MPI directly'
    End If

  End Subroutine assign_gpu_device

End Module gpu_device
