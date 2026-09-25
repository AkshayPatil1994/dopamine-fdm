!> Debug-only decomposition tracing: dumps every field at named points of the time step, assembled by global
!  index on rank 0, and checks each rank's ghost cells against the owning rank's interior value.
!
!  Enabled by the environment variable DOPAMINE_TRACE_DIR (unset = one cheap logical test per call, nothing else).
!  DOPAMINE_TRACE_STEPS (default 2) limits tracing to the first steps.
!
!  Outputs, per trace point and field:
!    <dir>/<seq>_<tag>.<field>   3 int32 (n1,n2,n3) + float64 global array, ghost cells included
!    <dir>/seam.log              one line per point/field: largest |ghost - owner interior| and where
!  Run the same case at two rank counts and compare the dump directories with tests/regression/trace_compare.py
!  to find the first point/field/cell at which the decomposition changes the result.
Module debug_trace

  Use iso_fortran_env, Only : Int32, Int64
  Use global
  Use mpi
  Use input_output, Only : distributed_block_range

  Implicit None

  Logical,             Save :: trace_on = .False., trace_ready = .False.
  Character(Len=1024), Save :: trace_dir = ''
  Integer(Int32),      Save :: trace_nsteps = 2, trace_seq = 0

Contains

  !> Trace point: no-op unless DOPAMINE_TRACE_DIR is set and istep is within DOPAMINE_TRACE_STEPS
  Subroutine trace_stage(tag)

    Character(Len=*), Intent(In) :: tag

    If ( .Not. trace_ready ) Call trace_setup
    If ( .Not. trace_on ) Return
    If ( istep > trace_nsteps ) Return

    trace_seq = trace_seq + 1

    !$acc update host(U,V,W,P)
    Call trace_field(tag, 'U', U, .True.,  .False.)
    Call trace_field(tag, 'V', V, .False., .False.)
    Call trace_field(tag, 'W', W, .False., .True. )
    Call trace_field(tag, 'P', P, .False., .False.)
    If ( Allocated(nu_t) .And. sgs_model /= 0 ) Then
       !$acc update host(nu_t)
       Call trace_field(tag, 'nut', nu_t, .False., .False.)
    End If
    If ( boussinesq_flag >= 1 ) Then
       !$acc update host(Tscal)
       Call trace_field(tag, 'T', Tscal, .False., .False.)
    End If
    If ( sediment_flag >= 1 ) Then
       !$acc update host(Cscal)
       Call trace_field(tag, 'C', Cscal, .False., .False.)
    End If

  End Subroutine trace_stage

  !> True while tracing is enabled and within the traced steps (guard for callers that must first sync device data to the host)
  Logical Function trace_active()
    If ( .Not. trace_ready ) Call trace_setup
    trace_active = trace_on .And. ( istep <= trace_nsteps )
  End Function trace_active

  !> Dump an interior-only array (indices 2..n-1 in each dim, as the RHS arrays) under the current sequence number.
  !  kind: 1 = cell-centred (nxg,nyg,nzg)
  Subroutine trace_interior(tag, name, A, kind)

    Character(Len=*), Intent(In) :: tag, name
    Real(Int64), Dimension(2:,2:,2:), Intent(In) :: A
    Integer(Int32), Intent(In) :: kind

    Real(Int64), Allocatable :: full(:,:,:)

    If ( .Not. trace_on ) Return
    If ( istep > trace_nsteps ) Return
    If ( kind /= 1 ) Return
    Allocate( full(nxg,nyg,nzg) )
    full = 0d0
    full(2:nxg-1,2:nyg-1,2:nzg-1) = A
    Call trace_field(tag, name, full, .False., .False.)
    Deallocate(full)

  End Subroutine trace_interior

  Subroutine trace_setup

    Character(Len=1024) :: val
    Integer(Int32)      :: length, stat, flag

    trace_ready = .True.
    flag = 0
    If ( myid == 0 ) Then
       Call Get_Environment_Variable('DOPAMINE_TRACE_DIR', val, length, stat)
       If ( stat == 0 .And. length > 0 ) Then
          trace_dir = val(1:length)
          flag = 1
       End If
       Call Get_Environment_Variable('DOPAMINE_TRACE_STEPS', val, length, stat)
       If ( stat == 0 .And. length > 0 ) Read(val(1:length),*) trace_nsteps
    End If
    Call Mpi_bcast(flag,         1, MPI_integer, 0, MPI_COMM_WORLD, ierr)
    Call Mpi_bcast(trace_nsteps, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)
    trace_on = ( flag == 1 )

  End Subroutine trace_setup

  !> Assemble one field by global index on rank 0. Interior cells (local 2..n-1 in x and z) are written first; every ghost
  !  cell is then compared with the owning rank's interior value (a mismatch is a stale or wrong halo), or, where no rank
  !  owns the cell (physical boundary ghost), stored as is.
  Subroutine trace_field(tag, name, F, is_x_face, is_z_face)

    Character(Len=*), Intent(In) :: tag, name
    Real(Int64), Dimension(:,:,:), Intent(In) :: F
    Logical, Intent(In) :: is_x_face, is_z_face

    Type :: blk_t
       Real(Int64), Allocatable :: d(:,:,:)
    End Type blk_t
    Type(blk_t), Allocatable :: blk(:)
    Real(Int64), Allocatable :: G(:,:,:)
    Logical,     Allocatable :: own(:,:), set(:,:)
    Integer(Int32) :: iproc, ix1, ix2, iz1, iz2, nxl, nzl, n1, n2, n3, i, j, k, gi, gk, nbad, worst(3), u
    Real(Int64)    :: dmax, dd, tol
    Character(Len=1024) :: fname

    n2 = Size(F,2)

    If ( myid /= 0 ) Then
       Call Mpi_send(F, Size(F), Mpi_real8, 0, 777, MPI_COMM_WORLD, ierr)
       Return
    End If

    n1 = Merge(nx_global, nxg_global, is_x_face)
    n3 = Merge(nz_global, nzg_global, is_z_face)
    Allocate( blk(0:nprocs-1), G(n1,n2,n3), own(n1,n3), set(n1,n3) )
    G = 0d0;  own = .False.;  set = .False.

    Do iproc = 0, nprocs-1
       Call distributed_block_range(iproc, is_x_face, is_z_face, ix1, ix2, iz1, iz2)
       nxl = ix2 - ix1 + 1;  nzl = iz2 - iz1 + 1
       Allocate( blk(iproc)%d(nxl,n2,nzl) )
       If ( iproc == 0 ) Then
          If ( Size(F,1) /= nxl .Or. Size(F,3) /= nzl ) Then
             Write(*,'(A,A,A,6(I0,1X))') ' trace: local shape mismatch for ', name, ': ', Size(F,1), nxl, Size(F,3), nzl
             Call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
          End If
          blk(0)%d = F
       Else
          Call Mpi_recv(blk(iproc)%d, nxl*n2*nzl, Mpi_real8, iproc, 777, MPI_COMM_WORLD, istat, ierr)
       End If
    End Do

    nbad = 0;  dmax = 0d0;  worst = 0

    ! Pass A: interiors (two ranks owning the same cell, e.g. a shared face, must also agree)
    Do iproc = 0, nprocs-1
       Call distributed_block_range(iproc, is_x_face, is_z_face, ix1, ix2, iz1, iz2)
       nxl = ix2 - ix1 + 1;  nzl = iz2 - iz1 + 1
       Do k = 2, nzl-1
          Do i = 2, nxl-1
             gi = ix1 + i - 1;  gk = iz1 + k - 1
             If ( own(gi,gk) ) Then
                Call note_mismatch(G(gi,:,gk), blk(iproc)%d(i,:,k), gi, gk)
             Else
                G(gi,:,gk) = blk(iproc)%d(i,:,k);  own(gi,gk) = .True.
             End If
          End Do
       End Do
    End Do

    ! Pass B: ghost cells
    Do iproc = 0, nprocs-1
       Call distributed_block_range(iproc, is_x_face, is_z_face, ix1, ix2, iz1, iz2)
       nxl = ix2 - ix1 + 1;  nzl = iz2 - iz1 + 1
       Do k = 1, nzl
          Do i = 1, nxl
             If ( i > 1 .And. i < nxl .And. k > 1 .And. k < nzl ) Cycle
             gi = ix1 + i - 1;  gk = iz1 + k - 1
             If ( own(gi,gk) .Or. set(gi,gk) ) Then
                Call note_mismatch(G(gi,:,gk), blk(iproc)%d(i,:,k), gi, gk)
             Else
                G(gi,:,gk) = blk(iproc)%d(i,:,k);  set(gi,gk) = .True.
             End If
          End Do
       End Do
    End Do

    Write(fname,'(A,A,I4.4,A,A,A,A)') Trim(trace_dir), '/', trace_seq, '_', Trim(tag), '.', Trim(name)
    Open(newunit=u, file=Trim(fname), access='stream', form='unformatted', status='replace')
    Write(u) n1, n2, n3
    Write(u) G
    Close(u)

    Open(newunit=u, file=Trim(trace_dir)//'/seam.log', position='append', action='write')
    Write(u,'(I4.4,1X,A,1X,A,1X,A,ES11.3,A,3(I0,1X),A,I0)') trace_seq, Trim(tag), Trim(name), &
         'seam_maxdiff=', dmax, ' at(i,j,k)=', worst, ' n_bad=', nbad
    Close(u)
    If ( nbad > 0 ) Write(*,'(A,I4.4,1X,A,1X,A,A,ES10.3,A,3(I0,1X),A,I0)') ' [trace] ', trace_seq, Trim(tag), Trim(name), &
         ' ghost/owner mismatch: max ', dmax, ' at ', worst, ' n_bad=', nbad

    Deallocate( blk, G, own, set )

  Contains

    Subroutine note_mismatch(a, b, gi_, gk_)
      Real(Int64),    Intent(In) :: a(:), b(:)
      Integer(Int32), Intent(In) :: gi_, gk_
      Integer(Int32) :: jj
      Do jj = 1, Size(a)
         dd  = Abs(a(jj) - b(jj))
         tol = 1d-13*(1d0 + Abs(a(jj)))
         If ( dd > tol ) Then
            nbad = nbad + 1
            If ( dd > dmax ) Then
               dmax = dd;  worst = [gi_, jj, gk_]
            End If
         End If
      End Do
    End Subroutine note_mismatch

  End Subroutine trace_field

End Module debug_trace
