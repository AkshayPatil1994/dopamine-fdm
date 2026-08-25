!> Slice/line probe output to binary files; format and layout
Module probe_output

  Use iso_fortran_env, Only : Int32, Int64
  Use global
  Use mpi

  Implicit None
  Private
  Public :: init_probes, write_probes, finalize_probes

  ! ── per-slice resolved state (set once in init_probes) ─────────────────
  Integer(Int32) :: sli_ax   (MAX_PROBES)    ! 1=x-normal 2=y-normal 3=z-normal
  Integer(Int32) :: sli_idx  (MAX_PROBES)    ! global 1-based interior cc index along normal
  Integer(Int32) :: sli_n1   (MAX_PROBES)    ! first  output dimension
  Integer(Int32) :: sli_n2   (MAX_PROBES)    ! second output dimension
  Logical        :: sli_cmask(4,MAX_PROBES)  ! component mask: U,V,W,P
  Integer(Int32) :: sli_ncomp(MAX_PROBES)    ! number of selected components
  Integer(Int32) :: sli_funit(MAX_PROBES)    ! Fortran file unit (rank 0)
  Integer(Int32) :: sli_nwrit(MAX_PROBES)    ! snapshots written so far

  ! ── per-line resolved state ─────────────────────────────────────────────
  Integer(Int32) :: lin_ax   (MAX_PROBES)    ! 1=along-x 2=along-y 3=along-z
  Integer(Int32) :: lin_t1   (MAX_PROBES)    ! first  transverse global cc index
  Integer(Int32) :: lin_t2   (MAX_PROBES)    ! second transverse global cc index
  Integer(Int32) :: lin_ks   (MAX_PROBES)    ! start index along line (global 1-based cc)
  Integer(Int32) :: lin_ke   (MAX_PROBES)    ! end   index along line (global 1-based cc)
  Integer(Int32) :: lin_npts (MAX_PROBES)    ! = ke - ks + 1
  Logical        :: lin_cmask(4,MAX_PROBES)  ! component mask
  Integer(Int32) :: lin_ncomp(MAX_PROBES)
  Integer(Int32) :: lin_funit(MAX_PROBES)
  Integer(Int32) :: lin_nwrit(MAX_PROBES)

Contains

  !> Resolve physical positions to grid indices, open output files, and write initial metadata (call once after grid setup)
  Subroutine init_probes

    Integer(Int32) :: n, meta_unit
    Character(300) :: fname
    Character(4)   :: dstr
    Character(8)   :: cstr

    sli_funit = 0;  sli_nwrit = 0
    lin_funit = 0;  lin_nwrit = 0

    ! ── 2-D slices ─────────────────────────────────────────────────────────
    Do n = 1, n_slices

      dstr = slice_dir(n)(1:4)
      cstr = slice_comps(n)(1:8)

      Select Case (dstr(1:1))
      Case ('x','X');  sli_ax(n) = 1
      Case ('y','Y');  sli_ax(n) = 2
      Case Default;    sli_ax(n) = 3   ! 'z' or unrecognised
      End Select

      sli_idx(n) = nearest_cc_global(sli_ax(n), slice_pos(n))

      Select Case (sli_ax(n))
      Case (1);  sli_n1(n) = nym_global;  sli_n2(n) = nzm_global
      Case (2);  sli_n1(n) = nxm_global;  sli_n2(n) = nzm_global
      Case (3);  sli_n1(n) = nxm_global;  sli_n2(n) = nym_global
      End Select

      Call parse_comps(cstr, sli_cmask(:,n), sli_ncomp(n))

      If (myid == 0) Then
        Write(fname,'(A,A)') Trim(slice_fileout(n)), '.bin'
        Call ensure_dir(Trim(fname))
        Open(newunit=sli_funit(n), file=Trim(fname), access='stream', &
             form='unformatted', status='replace', action='write')

        Write(fname,'(A,A)') Trim(slice_fileout(n)), '_meta.txt'
        Open(newunit=meta_unit, file=Trim(fname), form='formatted', status='replace')
        Write(meta_unit,'(A,I0)')   'ncomp  = ', sli_ncomp(n)
        Write(meta_unit,'(A,I0)')   'n1     = ', sli_n1(n)
        Write(meta_unit,'(A,I0)')   'n2     = ', sli_n2(n)
        Write(meta_unit,'(A,A)')    'dir    = ', Trim(dstr)
        Write(meta_unit,'(A,F16.8)') 'pos   = ', slice_pos(n)
        Write(meta_unit,'(A,A)')    'comps  = ', Trim(cstr)
        Write(meta_unit,'(A,I0)')   'nsnaps = 0'
        Close(meta_unit)

        Write(*,'(A,I2,A,A,A,I6,A,I6,A,I2)') &
          ' Probe slice ', n, ': ', Trim(dstr), '-normal  ', &
          sli_n1(n), ' x ', sli_n2(n), ',  ncomp = ', sli_ncomp(n)
      End If

    End Do

    ! ── 1-D lines ──────────────────────────────────────────────────────────
    Do n = 1, n_lines

      dstr = line_dir(n)(1:4)
      cstr = line_comps(n)(1:8)

      Select Case (dstr(1:1))
      Case ('x','X');  lin_ax(n) = 1
      Case ('y','Y');  lin_ax(n) = 2
      Case Default;    lin_ax(n) = 3
      End Select

      ! transverse positions → nearest cc indices
      Select Case (lin_ax(n))
      Case (1)   ! along x: pos1=y, pos2=z
        lin_t1(n) = nearest_cc_global(2, line_pos1(n))
        lin_t2(n) = nearest_cc_global(3, line_pos2(n))
      Case (2)   ! along y: pos1=x, pos2=z
        lin_t1(n) = nearest_cc_global(1, line_pos1(n))
        lin_t2(n) = nearest_cc_global(3, line_pos2(n))
      Case (3)   ! along z: pos1=x, pos2=y
        lin_t1(n) = nearest_cc_global(1, line_pos1(n))
        lin_t2(n) = nearest_cc_global(2, line_pos2(n))
      End Select

      Call line_extent(lin_ax(n), line_start(n), line_end(n), lin_ks(n), lin_ke(n))
      lin_npts(n) = lin_ke(n) - lin_ks(n) + 1

      Call parse_comps(cstr, lin_cmask(:,n), lin_ncomp(n))

      If (myid == 0) Then
        Write(fname,'(A,A)') Trim(line_fileout(n)), '.bin'
        Call ensure_dir(Trim(fname))
        Open(newunit=lin_funit(n), file=Trim(fname), access='stream', &
             form='unformatted', status='replace', action='write')

        Write(fname,'(A,A)') Trim(line_fileout(n)), '_meta.txt'
        Open(newunit=meta_unit, file=Trim(fname), form='formatted', status='replace')
        Write(meta_unit,'(A,I0)')  'ncomp  = ', lin_ncomp(n)
        Write(meta_unit,'(A,I0)')  'npts   = ', lin_npts(n)
        Write(meta_unit,'(A,A)')   'dir    = ', Trim(dstr)
        Write(meta_unit,'(A,A)')   'comps  = ', Trim(cstr)
        Write(meta_unit,'(A,I0)')  'nsnaps = 0'
        Close(meta_unit)

        Write(*,'(A,I2,A,A,A,I6,A,I2)') &
          ' Probe line  ', n, ': along-', Trim(dstr), ',  npts = ', &
          lin_npts(n), ',  ncomp = ', lin_ncomp(n)
      End If

    End Do

  End Subroutine init_probes


  !  write_probes
  !  Called every step.  Output is gated on slice_freq / line_freq.
  Subroutine write_probes

    Integer(Int32) :: n

    If ( n_slices > 0 .And. slice_freq > 0 ) Then
      If ( Mod(istep, slice_freq) == 0 ) Then
        Do n = 1, n_slices
          Call write_slice_n(n)
        End Do
      End If
    End If

    If ( n_lines > 0 .And. line_freq > 0 ) Then
      If ( Mod(istep, line_freq) == 0 ) Then
        Do n = 1, n_lines
          Call write_line_n(n)
        End Do
      End If
    End If

  End Subroutine write_probes


  !  finalize_probes
  !  Close all open file units (rank 0 only).
  Subroutine finalize_probes

    Integer(Int32) :: n

    If ( myid /= 0 ) Return

    Do n = 1, n_slices
      If (sli_funit(n) > 0) Close(sli_funit(n))
    End Do

    Do n = 1, n_lines
      If (lin_funit(n) > 0) Close(lin_funit(n))
    End Do

  End Subroutine finalize_probes


  !  Private: write_slice_n
  !  Extract, MPI-reduce, and append one 2-D slice snapshot.
  Subroutine write_slice_n(n)

    Integer(Int32), Intent(In) :: n

    Real(Int64), Allocatable :: lbuf(:,:), gbuf(:,:), out2d(:,:,:)
    Integer(Int32) :: comp, c, ia, ja, ka, ig, jg, kg_g
    Integer(Int32) :: n1, n2, nc, fix_ia, fix_ja, fix_ka

    n1 = sli_n1(n);  n2 = sli_n2(n);  nc = sli_ncomp(n)

    Allocate( lbuf(n1,n2) )
    Allocate( gbuf(n1,n2) )
    If (myid == 0) Allocate( out2d(nc,n1,n2) )

    c = 0
    Do comp = 1, 4
      If ( .Not. sli_cmask(comp,n) ) Cycle
      c = c + 1
      lbuf = 0d0

      Select Case (sli_ax(n))

      !── x-normal: fixed x-index, output(jg, kg_g) ─────────────────────!
      Case (1)
        fix_ia = sli_idx(n) + 1   ! global cc 1-based → ghost array index
        Do ka = 2, nzg-1
          kg_g = kg1_global(myid) + ka - 2   ! global 1-based interior cc z
          If (kg_g < 1 .Or. kg_g > nzm_global) Cycle
          Do ja = 2, nyg-1
            jg = ja - 1
            lbuf(jg, kg_g) = cc_val(comp, fix_ia, ja, ka)
          End Do
        End Do

      !── y-normal: fixed y-index, output(ig, kg_g) ─────────────────────!
      Case (2)
        fix_ja = sli_idx(n) + 1
        Do ka = 2, nzg-1
          kg_g = kg1_global(myid) + ka - 2
          If (kg_g < 1 .Or. kg_g > nzm_global) Cycle
          Do ia = 2, nxg-1
            ig = ia - 1
            lbuf(ig, kg_g) = cc_val(comp, ia, fix_ja, ka)
          End Do
        End Do

      !── z-normal: fixed z-index, output(ig, jg) ───────────────────────!
      Case (3)
        ! global cc index sli_idx(n) (1-based) → local array index
        fix_ka = sli_idx(n) - kg1_global(myid) + 2
        If (fix_ka >= 2 .And. fix_ka <= nzg-1) Then
          Do ja = 2, nyg-1
            jg = ja - 1
            Do ia = 2, nxg-1
              ig = ia - 1
              lbuf(ig, jg) = cc_val(comp, ia, ja, fix_ka)
            End Do
          End Do
        End If

      End Select

      Call MPI_Reduce(lbuf(1,1), gbuf(1,1), n1*n2, MPI_REAL8, MPI_SUM, &
                      0, MPI_COMM_WORLD, ierr)

      If (myid == 0) out2d(c,:,:) = gbuf

    End Do

    If (myid == 0) Then
      Write(sli_funit(n)) out2d
      sli_nwrit(n) = sli_nwrit(n) + 1
      Call update_slice_meta(n)
      Deallocate(out2d)
    End If

    Deallocate(lbuf, gbuf)

  End Subroutine write_slice_n


  !  Private: write_line_n
  !  Extract, MPI-reduce, and append one 1-D line snapshot.
  Subroutine write_line_n(n)

    Integer(Int32), Intent(In) :: n

    Real(Int64), Allocatable :: lbuf(:), gbuf(:), out1d(:,:)
    Integer(Int32) :: comp, c, ia, ja, ka, ig_g, jg_g, kg_g, lp
    Integer(Int32) :: npts, nc, fix_ka, fix_ia, fix_ja

    npts = lin_npts(n);  nc = lin_ncomp(n)

    Allocate( lbuf(npts) )
    Allocate( gbuf(npts) )
    If (myid == 0) Allocate( out1d(nc,npts) )

    c = 0
    Do comp = 1, 4
      If ( .Not. lin_cmask(comp,n) ) Cycle
      c = c + 1
      lbuf = 0d0

      Select Case (lin_ax(n))

      !── line along x: fixed y (t1) and z (t2) ─────────────────────────!
      Case (1)
        fix_ka = lin_t2(n) - kg1_global(myid) + 2   ! local z-array index
        If (fix_ka >= 2 .And. fix_ka <= nzg-1) Then
          fix_ja = lin_t1(n) + 1   ! y ghost-array index
          Do lp = 1, npts
            ig_g = lin_ks(n) + lp - 1   ! global 1-based cc x index
            ia   = ig_g + 1
            lbuf(lp) = cc_val(comp, ia, fix_ja, fix_ka)
          End Do
        End If

      !── line along y: fixed x (t1) and z (t2) ─────────────────────────!
      Case (2)
        fix_ka = lin_t2(n) - kg1_global(myid) + 2
        If (fix_ka >= 2 .And. fix_ka <= nzg-1) Then
          fix_ia = lin_t1(n) + 1
          Do lp = 1, npts
            jg_g = lin_ks(n) + lp - 1
            ja   = jg_g + 1
            lbuf(lp) = cc_val(comp, fix_ia, ja, fix_ka)
          End Do
        End If

      !── line along z: fixed x (t1) and y (t2), each rank fills slab ───!
      Case (3)
        fix_ia = lin_t1(n) + 1
        fix_ja = lin_t2(n) + 1
        Do ka = 2, nzg-1
          kg_g = kg1_global(myid) + ka - 2   ! global 1-based cc z
          If (kg_g < lin_ks(n) .Or. kg_g > lin_ke(n)) Cycle
          lp = kg_g - lin_ks(n) + 1
          lbuf(lp) = cc_val(comp, fix_ia, fix_ja, ka)
        End Do

      End Select

      Call MPI_Reduce(lbuf(1), gbuf(1), npts, MPI_REAL8, MPI_SUM, &
                      0, MPI_COMM_WORLD, ierr)

      If (myid == 0) out1d(c,:) = gbuf

    End Do

    If (myid == 0) Then
      Write(lin_funit(n)) out1d
      lin_nwrit(n) = lin_nwrit(n) + 1
      Call update_line_meta(n)
      Deallocate(out1d)
    End If

    Deallocate(lbuf, gbuf)

  End Subroutine write_line_n


  !  cc_val — cell-centre interpolated value for component comp
  !  at ghost-array position (ia, ja, ka).
  Real(Int64) Function cc_val(comp, ia, ja, ka)

    Integer(Int32), Intent(In) :: comp, ia, ja, ka

    Select Case (comp)
    Case (1);  cc_val = 0.5d0 * (U(ia-1,ja,ka) + U(ia,ja,ka))
    Case (2);  cc_val = 0.5d0 * (V(ia,ja-1,ka) + V(ia,ja,ka))
    Case (3);  cc_val = 0.5d0 * (W(ia,ja,ka-1) + W(ia,ja,ka))
    Case Default;  cc_val = P(ia,ja,ka)
    End Select

  End Function cc_val


  !  nearest_cc_global — find the nearest global 1-based interior cc index
  !  for physical coordinate pos along axis (1=x, 2=y, 3=z).
  Integer(Int32) Function nearest_cc_global(axis, pos)

    Integer(Int32), Intent(In) :: axis
    Real(Int64),    Intent(In) :: pos

    Integer(Int32) :: i, best, n_pts
    Real(Int64)    :: d, d_best

    Select Case (axis)
    Case (1)   ! x cell centres: xg_global(2..nxg_global-1)
      n_pts  = nxm_global
      best   = 1
      d_best = Abs(xg_global(2) - pos)
      Do i = 2, nxg_global-1
        d = Abs(xg_global(i) - pos)
        If (d < d_best) Then;  d_best = d;  best = i - 1;  End If
      End Do
    Case (2)   ! y cell centres: yg_global(2..nyg_global-1)
      n_pts  = nym_global
      best   = 1
      d_best = Abs(yg_global(2) - pos)
      Do i = 2, nyg_global-1
        d = Abs(yg_global(i) - pos)
        If (d < d_best) Then;  d_best = d;  best = i - 1;  End If
      End Do
    Case Default   ! z cell centres: zg_global(2..nzg_global-1)
      n_pts  = nzm_global
      best   = 1
      d_best = Abs(zg_global(2) - pos)
      Do i = 2, nzg_global-1
        d = Abs(zg_global(i) - pos)
        If (d < d_best) Then;  d_best = d;  best = i - 1;  End If
      End Do
    End Select

    nearest_cc_global = Max(1, Min(n_pts, best))

  End Function nearest_cc_global


  !  line_extent — convert physical start/end to 1-based interior cc index
  !  range along the given axis.
  Subroutine line_extent(axis, pos_s, pos_e, ks, ke)

    Integer(Int32), Intent(In)  :: axis
    Real(Int64),    Intent(In)  :: pos_s, pos_e
    Integer(Int32), Intent(Out) :: ks, ke

    Integer(Int32) :: n_pts

    Select Case (axis)
    Case (1);  n_pts = nxm_global
    Case (2);  n_pts = nym_global
    Case Default;  n_pts = nzm_global
    End Select

    If (pos_s <= 0d0 .And. pos_e >= 1d29) Then
      ks = 1;  ke = n_pts   ! default: full extent
    Else
      ks = nearest_cc_global(axis, pos_s)
      ke = nearest_cc_global(axis, pos_e)
      If (ks > ke) Then;  ks = ke;  ke = n_pts;  End If
    End If

  End Subroutine line_extent


  !> Build a 4-element logical mask (U,V,W,P) from str ('U'/'V'/'W'/'P', case-insensitive; empty = all three velocities)
  Subroutine parse_comps(str, mask, ncomp)

    Character(*),   Intent(In)  :: str
    Logical,        Intent(Out) :: mask(4)
    Integer(Int32), Intent(Out) :: ncomp

    Character(Len(str)) :: ustr
    Integer(Int32)      :: i

    ustr = str
    Do i = 1, Len_Trim(ustr)
      If (ustr(i:i) >= 'a' .And. ustr(i:i) <= 'z') &
          ustr(i:i) = Achar(Iachar(ustr(i:i)) - 32)
    End Do

    mask(1) = Index(Trim(ustr), 'U') > 0
    mask(2) = Index(Trim(ustr), 'V') > 0
    mask(3) = Index(Trim(ustr), 'W') > 0
    mask(4) = Index(Trim(ustr), 'P') > 0

    If ( .Not. (mask(1) .Or. mask(2) .Or. mask(3) .Or. mask(4)) ) Then
      mask = [.True., .True., .True., .False.]   ! default: U,V,W
    End If

    ncomp = Count(mask)

  End Subroutine parse_comps


  !  ensure_dir — create directory component of file path if missing.
  !  Called on rank 0 only.
  Subroutine ensure_dir(fname)

    Character(*), Intent(In) :: fname

    Integer(Int32) :: slash
    Logical        :: exists

    slash = Index(Trim(fname), '/', Back=.True.)
    If (slash < 1) Return
    Inquire(file=fname(1:slash-1)//'/.', exist=exists)
    If (.Not. exists) &
      Call execute_command_line('mkdir -p ' // fname(1:slash-1), wait=.True.)

  End Subroutine ensure_dir


  !  update_slice_meta — rewrite slice n's meta file with current nsnaps.
  Subroutine update_slice_meta(n)

    Integer(Int32), Intent(In) :: n

    Integer(Int32) :: u
    Character(300) :: fname

    Write(fname,'(A,A)') Trim(slice_fileout(n)), '_meta.txt'
    Open(newunit=u, file=Trim(fname), form='formatted', status='replace')
    Write(u,'(A,I0)')    'ncomp  = ', sli_ncomp(n)
    Write(u,'(A,I0)')    'n1     = ', sli_n1(n)
    Write(u,'(A,I0)')    'n2     = ', sli_n2(n)
    Write(u,'(A,A)')     'dir    = ', Trim(slice_dir(n))
    Write(u,'(A,F16.8)') 'pos    = ', slice_pos(n)
    Write(u,'(A,A)')     'comps  = ', Trim(slice_comps(n))
    Write(u,'(A,I0)')    'nsnaps = ', sli_nwrit(n)
    Close(u)

  End Subroutine update_slice_meta


  !  update_line_meta — rewrite line n's meta file with current nsnaps.
  Subroutine update_line_meta(n)

    Integer(Int32), Intent(In) :: n

    Integer(Int32) :: u
    Character(300) :: fname

    Write(fname,'(A,A)') Trim(line_fileout(n)), '_meta.txt'
    Open(newunit=u, file=Trim(fname), form='formatted', status='replace')
    Write(u,'(A,I0)') 'ncomp  = ', lin_ncomp(n)
    Write(u,'(A,I0)') 'npts   = ', lin_npts(n)
    Write(u,'(A,A)')  'dir    = ', Trim(line_dir(n))
    Write(u,'(A,A)')  'comps  = ', Trim(line_comps(n))
    Write(u,'(A,I0)') 'nsnaps = ', lin_nwrit(n)
    Close(u)

  End Subroutine update_line_meta

End Module probe_output
