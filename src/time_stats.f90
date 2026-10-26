!> On-the-fly time-averaged velocity statistics (mean, Reynolds stresses) and optional scalar statistics
!  (mean, variance, flux), accumulated on the device each step and written, little-endian, next to every full field
!  snapshot as fields/<fileout>_tavg.<step>. File layout: Int32 nxm,nym,nzm,nf; Real64 t_acc; nf x
!  Character(8) names; then nf blocks of Real64 (nxm,nym,nzm), Fortran order. With chunked_snapshot the same
!  layout holds for the chunk box (sizes of the box) plus a trailer: Int32 f1(3), f2(3); Real64 xm, ym, zm of
!  the box cells.
Module time_stats

  Use iso_fortran_env, Only : Int32, Int64
  Use global
  Use mpi
  Use reynolds_stress_budget, Only : reduce_to_rank0_4d

  Implicit None
  Private
  Public :: init_tavg, accumulate_tavg, write_tavg, finalize_tavg

  Integer(Int32), Parameter :: NAME_LEN = 8
  Integer(Int32) :: nf = 0, iT0 = 0, iC0 = 0   ! field count; first field of the T and C blocks (0: block absent)
  Real(Int64)    :: t_acc = 0d0   ! time accumulated (the sums are dt-weighted, so adaptive dt does not bias them)

  ! Sums of dt*(x): U,V,W, UU,VV,WW,UV,UW,VW, then per scalar S: S, SS, uS, vS, wS
  Real(Int64), Allocatable :: acc(:,:,:,:)
  ! output field c = <raw c> - <raw mp(c)><raw mq(c)>; mp = 0 for plain means
  Integer(Int32), Allocatable :: mp(:), mq(:)
  Character(NAME_LEN), Allocatable :: tname(:)

Contains

  !> Allocate the sums, set the field table and, on a cumulative restart, reload the sums from the
  !  snapshot-step file
  Subroutine init_tavg

    Integer(Int32) :: c, a, b, ab(2,6), base
    Character(300) :: path
    Character(1)   :: sname(2)
    Logical        :: have
    Character(NAME_LEN) :: vname(9)

    If ( tavg_active /= 1 ) Return
    If ( tavg_T == 1 .And. boussinesq_flag < 1 ) Call tavg_abort('tavg_T = 1 needs boussinesq_flag >= 1')
    If ( tavg_C == 1 .And. sediment_flag < 1 )   Call tavg_abort('tavg_C = 1 needs sediment_flag >= 1')

    nf = 9
    If ( tavg_T == 1 ) Then
       iT0 = nf + 1;  nf = nf + 5
    End If
    If ( tavg_C == 1 ) Then
       iC0 = nf + 1;  nf = nf + 5
    End If
    Allocate( acc(nxm, nym_global, nzm, nf), mp(nf), mq(nf), tname(nf) )
    acc = 0d0;  mp = 0;  mq = 0

    vname = [ 'U       ', 'V       ', 'W       ', "u'u'    ", "v'v'    ", "w'w'    ", "u'v'    ", "u'w'    ", "v'w'    " ]
    tname(1:9) = vname
    ab = Reshape( [ 1,1, 2,2, 3,3, 1,2, 1,3, 2,3 ], [ 2, 6 ] )
    Do c = 1, 6
       mp(3+c) = ab(1,c);  mq(3+c) = ab(2,c)
    End Do

    sname = [ 'T', 'C' ]
    Do b = 1, 2
       base = Merge(iT0, iC0, b == 1)
       If ( base == 0 ) Cycle
       tname(base)   = sname(b)
       tname(base+1) = sname(b) // "'" // sname(b) // "'"
       mp(base+1) = base;  mq(base+1) = base
       Do a = 1, 3
          tname(base+1+a) = Achar(Iachar('u')+a-1) // "'" // sname(b) // "'"
          mp(base+1+a) = a;  mq(base+1+a) = base
       End Do
    End Do

    If ( restart == 1 .And. tavg_reset == 0 ) Then
       path = tavg_path(nstep_init)
       have = .False.
       If ( myid == 0 ) Inquire(file=Trim(path), exist=have)
       Call MPI_Bcast(have, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
       If ( have ) Then
          Call read_tavg(path)
       Else If ( tavg_nstart <= nstep_init ) Then
          Call tavg_abort('restart file ' // Trim(path) // ' is missing and the earlier average would be lost; ' &
               // 'set tavg_reset = 1 or tavg_nstart > nstep_init to start fresh')
       End If
    End If

    !$acc enter data copyin(acc)

  End Subroutine init_tavg


  Subroutine finalize_tavg

    If ( tavg_active /= 1 ) Return
    !$acc exit data delete(acc)
    Deallocate( acc, mp, mq, tname )

  End Subroutine finalize_tavg


  !> Add the step just completed (call after the step, with the final U,V,W and scalars)
  Subroutine accumulate_tavg

    Integer(Int32) :: i, j, k
    Real(Int64)    :: uc, vc, wc, dts

    If ( tavg_active /= 1 ) Return
    If ( istep + nstep_init < tavg_nstart ) Return

    dts = dt_step
    !$acc parallel loop collapse(3) present(acc,U,V,W)
    Do k = 2, nzg-1
       Do j = 2, nyg-1
          Do i = 2, nxg-1
             uc = 0.5d0 * ( U(i-1,j,k) + U(i,j,k) )
             vc = 0.5d0 * ( V(i,j-1,k) + V(i,j,k) )
             wc = 0.5d0 * ( W(i,j,k-1) + W(i,j,k) )
             acc(i-1,j-1,k-1,1) = acc(i-1,j-1,k-1,1) + dts*uc
             acc(i-1,j-1,k-1,2) = acc(i-1,j-1,k-1,2) + dts*vc
             acc(i-1,j-1,k-1,3) = acc(i-1,j-1,k-1,3) + dts*wc
             acc(i-1,j-1,k-1,4) = acc(i-1,j-1,k-1,4) + dts*uc*uc
             acc(i-1,j-1,k-1,5) = acc(i-1,j-1,k-1,5) + dts*vc*vc
             acc(i-1,j-1,k-1,6) = acc(i-1,j-1,k-1,6) + dts*wc*wc
             acc(i-1,j-1,k-1,7) = acc(i-1,j-1,k-1,7) + dts*uc*vc
             acc(i-1,j-1,k-1,8) = acc(i-1,j-1,k-1,8) + dts*uc*wc
             acc(i-1,j-1,k-1,9) = acc(i-1,j-1,k-1,9) + dts*vc*wc
          End Do
       End Do
    End Do
    !$acc end parallel loop

    If ( iT0 > 0 ) Call accumulate_scalar(Tscal, iT0)
    If ( iC0 > 0 ) Call accumulate_scalar(Cscal, iC0)
    t_acc = t_acc + dt_step

  End Subroutine accumulate_tavg


  Subroutine accumulate_scalar(S, base)

    Real(Int64),    Intent(In) :: S(:,:,:)
    Integer(Int32), Intent(In) :: base

    Integer(Int32) :: i, j, k
    Real(Int64)    :: uc, vc, wc, sc, dts

    dts = dt_step
    !$acc parallel loop collapse(3) present(acc,S,U,V,W)
    Do k = 2, nzg-1
       Do j = 2, nyg-1
          Do i = 2, nxg-1
             uc = 0.5d0 * ( U(i-1,j,k) + U(i,j,k) )
             vc = 0.5d0 * ( V(i,j-1,k) + V(i,j,k) )
             wc = 0.5d0 * ( W(i,j,k-1) + W(i,j,k) )
             sc = S(i,j,k)
             acc(i-1,j-1,k-1,base  ) = acc(i-1,j-1,k-1,base  ) + dts*sc
             acc(i-1,j-1,k-1,base+1) = acc(i-1,j-1,k-1,base+1) + dts*sc*sc
             acc(i-1,j-1,k-1,base+2) = acc(i-1,j-1,k-1,base+2) + dts*uc*sc
             acc(i-1,j-1,k-1,base+3) = acc(i-1,j-1,k-1,base+3) + dts*vc*sc
             acc(i-1,j-1,k-1,base+4) = acc(i-1,j-1,k-1,base+4) + dts*wc*sc
          End Do
       End Do
    End Do
    !$acc end parallel loop

  End Subroutine accumulate_scalar


  !> Collective: write the averages as of now: the full domain to fields/<fileout>_tavg.<step> when do_full and,
  !  with chunked_snapshot, the chunk box to fields/<fileout>_tavgBlock.<step> when do_chunk (plain _tavg when
  !  full_snapshot is off). With tavg_reset = 1 the sums restart after the primary file (the full one, else the
  !  chunk one), so a chunk-only write cannot drop the average the next restart file needs
  Subroutine write_tavg(step, do_full, do_chunk)

    Integer(Int32), Intent(In) :: step
    Logical,        Intent(In) :: do_full, do_chunk

    Logical :: reset

    If ( tavg_active /= 1 ) Return
    If ( t_acc <= 0d0 ) Return

    If ( do_full .And. ( full_snapshot .Or. .Not. chunked_snapshot ) ) Then
       !$acc update host(acc)
       Call write_tavg_file(tavg_path(step), .False.)
    End If
    If ( do_chunk ) Then
       !$acc update host(acc)
       If ( full_snapshot ) Then
          Call write_tavg_file(tavg_path(step, 'Block'), .True.)
       Else
          Call write_tavg_file(tavg_path(step), .True.)
       End If
    End If

    reset = ( do_full .And. full_snapshot ) .Or. ( do_chunk .And. .Not. full_snapshot )
    If ( tavg_reset == 1 .And. reset ) Then
       t_acc = 0d0
       !$acc kernels present(acc)
       acc = 0d0
       !$acc end kernels
    End If

  End Subroutine write_tavg


  !> Collective: write all fields of the whole domain, or of the chunk box followed by the trailer f1, f2, xm,
  !  ym, zm of its cells
  Subroutine write_tavg_file(path, box)

    Character(*), Intent(In) :: path
    Logical,      Intent(In) :: box

    Real(Int64), Allocatable :: buf(:,:,:,:), gbuf(:,:,:,:), bbuf(:,:,:)
    Real(Int64)    :: dn
    Integer(Int32) :: c, unit_no, nb(3)

    dn = 1d0 / t_acc
    nb = [ nxm_global, nym_global, nzm_global ]
    If ( box ) nb = chunk_f2 - chunk_f1

    Allocate( buf(1, nxm, nym_global, nzm) )
    If ( myid == 0 ) Then
       If ( box ) Then
          Allocate( bbuf(nb(1), nb(2), nb(3)) )
       Else
          Allocate( gbuf(1, nb(1), nb(2), nb(3)) )
       End If
       Open(newunit=unit_no, file=Trim(path), access='stream', form='unformatted', action='write', status='replace', &
            convert='little_endian')
       Write(unit_no) nb, nf
       Write(unit_no) t_acc
       Write(unit_no) tname
    Else
       Allocate( gbuf(1,1,1,1), bbuf(1,1,1) )
    End If

    Do c = 1, nf
       buf(1,:,:,:) = acc(:,:,:,c) * dn
       If ( mp(c) > 0 ) buf(1,:,:,:) = buf(1,:,:,:) - acc(:,:,:,mp(c)) * acc(:,:,:,mq(c)) * dn * dn
       If ( ibm_input_mode >= 1 ) Then
          Where ( phi(2:nxg-1,2:nyg-1,2:nzg-1) < 0d0 ) buf(1,:,:,:) = 0d0
       End If
       If ( box ) Then
          Call gather_box(buf(1,:,:,:), bbuf)
          If ( myid == 0 ) Write(unit_no) bbuf
       Else
          Call reduce_to_rank0_4d(buf, gbuf, 1, nxm, nym_global, nzm)
          If ( myid == 0 ) Write(unit_no) gbuf
       End If
    End Do

    If ( myid == 0 ) Then
       If ( box ) Write(unit_no) chunk_f1, chunk_f2, xm_global(chunk_f1(1):chunk_f2(1)-1), &
            ym_global(chunk_f1(2):chunk_f2(2)-1), zm_global(chunk_f1(3):chunk_f2(3)-1)
       Close(unit_no)
       Write(*,'(A,A,A,ES10.3)') '  Wrote ', Trim(path), '  (averaging time ', t_acc
    End If

  End Subroutine write_tavg_file


  !> Global cell range lo:hi of the chunk box owned by rank r (empty when any hi < lo)
  Subroutine box_cells(r, lo, hi)

    Integer(Int32), Intent(In)  :: r
    Integer(Int32), Intent(Out) :: lo(3), hi(3)

    lo = [ Max(chunk_f1(1), ig1_global(r)), chunk_f1(2), Max(chunk_f1(3), kg1_global(r)) ]
    hi = [ Min(chunk_f2(1)-1, ig2_global(r)-2), chunk_f2(2)-1, Min(chunk_f2(3)-1, kg2_global(r)-2) ]

  End Subroutine box_cells


  !> Collective: assemble the chunk-box cells of the local field buf(nxm,nym_global,nzm) in gbuf on rank 0
  Subroutine gather_box(buf, gbuf)

    Real(Int64), Intent(In)    :: buf(:,:,:)
    Real(Int64), Intent(InOut) :: gbuf(:,:,:)

    Real(Int64), Allocatable :: tmp(:,:,:)
    Integer(Int32) :: r, lo(3), hi(3), i0, k0

    Do r = 0, nprocs-1
       If ( myid /= 0 .And. myid /= r ) Cycle
       Call box_cells(r, lo, hi)
       If ( Any(hi < lo) ) Cycle
       i0 = ig1_global(r) - 1
       k0 = kg1_global(r) - 1
       Allocate( tmp(hi(1)-lo(1)+1, hi(2)-lo(2)+1, hi(3)-lo(3)+1) )
       If ( myid == r ) tmp = buf(lo(1)-i0:hi(1)-i0, lo(2):hi(2), lo(3)-k0:hi(3)-k0)
       If ( myid == r .And. r /= 0 ) Then
          Call MPI_Send(tmp, Size(tmp), MPI_REAL8, 0, r, MPI_COMM_WORLD, ierr)
       Else If ( myid == 0 .And. r /= 0 ) Then
          Call MPI_Recv(tmp, Size(tmp), MPI_REAL8, r, r, MPI_COMM_WORLD, istat, ierr)
       End If
       If ( myid == 0 ) gbuf(lo(1)-chunk_f1(1)+1:hi(1)-chunk_f1(1)+1, lo(2)-chunk_f1(2)+1:hi(2)-chunk_f1(2)+1, &
                             lo(3)-chunk_f1(3)+1:hi(3)-chunk_f1(3)+1) = tmp
       Deallocate( tmp )
    End Do

  End Subroutine gather_box


  !> Reload the sums from a file written by write_tavg: rank 0 reads, each rank receives its block
  Subroutine read_tavg(path)

    Character(*), Intent(In) :: path

    Real(Int64), Allocatable :: buf(:,:,:,:), gbuf(:,:,:,:)
    Character(NAME_LEN), Allocatable :: fnames(:)
    Integer(Int32) :: c, unit_no, hdr(4)

    Allocate( buf(1, nxm, nym_global, nzm) )
    If ( myid == 0 ) Then
       Open(newunit=unit_no, file=Trim(path), access='stream', form='unformatted', action='read', convert='little_endian')
       Read(unit_no) hdr
       If ( hdr(1) /= nxm_global .Or. hdr(2) /= nym_global .Or. hdr(3) /= nzm_global .Or. hdr(4) /= nf ) &
            Call tavg_abort('restart file ' // Trim(path) // ' does not match the grid or the tavg_T/tavg_C setting')
       Read(unit_no) t_acc
       Allocate( fnames(nf) )
       Read(unit_no) fnames
       If ( Any( fnames /= tname ) ) Call tavg_abort('restart file ' // Trim(path) // ' has different field names')
       Allocate( gbuf(1, nxm_global, nym_global, nzm_global) )
    Else
       Allocate( gbuf(1,1,1,1) )
    End If
    Call MPI_Bcast(t_acc, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

    Do c = 1, nf
       If ( myid == 0 ) Read(unit_no) gbuf
       Call scatter_from_rank0(gbuf, buf)
       acc(:,:,:,c) = buf(1,:,:,:) * t_acc
       ! the file holds central moments; the sums need the raw ones (the means are fields 1-3 and the scalar
       ! means, read first)
       If ( mp(c) > 0 ) acc(:,:,:,c) = acc(:,:,:,c) + acc(:,:,:,mp(c)) * acc(:,:,:,mq(c)) / t_acc
    End Do

    If ( myid == 0 ) Then
       Close(unit_no)
       Write(*,'(A,A,A,ES10.3)') ' tavg restart: read ', Trim(path), ', averaging time ', t_acc
    End If
    Deallocate( buf, gbuf )

  End Subroutine read_tavg


  !> Inverse of reduce_to_rank0_4d for one component: rank 0 sends each rank its interior block
  Subroutine scatter_from_rank0(gbuf, buf)

    Real(Int64), Intent(In)  :: gbuf(:,:,:,:)
    Real(Int64), Intent(Out) :: buf(:,:,:,:)

    Real(Int64), Allocatable :: tmp(:,:,:,:)
    Integer(Int32) :: iproc, rnx, rnz, i0, k0

    If ( myid /= 0 ) Then
       Call MPI_Recv(buf, Size(buf), MPI_REAL8, 0, myid, MPI_COMM_WORLD, istat, ierr)
       Return
    End If

    Do iproc = 0, nprocs-1
       rnx = ig2_global(iproc) - ig1_global(iproc) - 1
       rnz = kg2_global(iproc) - kg1_global(iproc) - 1
       i0  = ig1_global(iproc)
       k0  = kg1_global(iproc)
       If ( iproc == 0 ) Then
          buf = gbuf(:, i0:i0+rnx-1, :, k0:k0+rnz-1)
       Else
          Allocate( tmp(1, rnx, nym_global, rnz) )
          tmp = gbuf(:, i0:i0+rnx-1, :, k0:k0+rnz-1)
          Call MPI_Send(tmp, Size(tmp), MPI_REAL8, iproc, iproc, MPI_COMM_WORLD, ierr)
          Deallocate( tmp )
       End If
    End Do

  End Subroutine scatter_from_rank0


  Function tavg_path(step, tag) Result(path)

    Integer(Int32), Intent(In) :: step
    Character(*), Intent(In), Optional :: tag
    Character(300) :: path
    Character(12)  :: ext

    Write(ext,'(I0)') step
    path = 'fields/' // Trim(Adjustl(fileout)) // '_tavg'
    If ( Present(tag) ) path = Trim(path) // tag
    path = Trim(path) // '.' // Trim(ext)

  End Function tavg_path


  Subroutine tavg_abort(msg)

    Character(*), Intent(In) :: msg

    If ( myid == 0 ) Write(*,'(A,A)') ' ERROR: time_stats: ', msg
    Call MPI_Abort(MPI_COMM_WORLD, 1, ierr)

  End Subroutine tavg_abort

End Module time_stats
