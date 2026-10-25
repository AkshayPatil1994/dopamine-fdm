!> Host-side padded copies of a decomposed field: E extra planes beyond the one ghost plane, in x and z.
!  Needed by anything that samples at displaced points (particle RK stages, IBM image points, SDF
!  differences), which would clamp near a rank seam and make the result depend on the rank layout. Plane rules
!  follow apply_periodic_bc_x/z. Across a seam the high side takes the neighbour's planes 3..2+E, the low side
!  its planes n-1-E..n-2 (n = its plane count). Across the periodic wrap the first rank supplies planes 4..3+E
!  (cell-centred) or 3..2+E (face-located), the last rank n-2-E..n-3 or n-1-E..n-2. At a non-periodic edge the
!  edge plane is repeated.
Module halo_pad

  Use iso_fortran_env, Only : Int32, Int64
  Use global
  Use mpi
  Use decomp, Only : x_halo_neighbors, z_halo_neighbors, x_periodic_partner, z_periodic_partner

  Implicit None

  ! staging planes of pad_field_dev (device copies kept in step)
  Real(Int64), Allocatable, Dimension(:) :: pad_bs, pad_br

Contains

  !> Make the staging buffers hold at least n values (device copies included)
  Subroutine pad_buffers(n)

    Integer(Int32), Intent(In) :: n

    If ( Allocated(pad_bs) ) Then
       If ( Size(pad_bs) >= n ) Return
       !$acc exit data delete(pad_bs,pad_br)
       Deallocate( pad_bs, pad_br )
    End If
    Allocate( pad_bs(n), pad_br(n) )
    pad_bs = 0d0;  pad_br = 0d0
    !$acc enter data create(pad_bs,pad_br)

  End Subroutine pad_buffers


  !> Padded axis apad(1-E:n+E) from the local axis a(1:n): the neighbouring points from the global axis where
  !  they exist, otherwise extrapolated with the edge spacing (a periodic edge; the axes involved there are
  !  uniform)
  Subroutine pad_axis(a, n, aglob, ng, g1, E, apad)

    Integer(Int32), Intent(In)  :: n, ng, g1, E
    Real   (Int64), Intent(In)  :: a(n), aglob(ng)
    Real   (Int64), Intent(Out) :: apad(1-E:n+E)

    Integer(Int32) :: m

    apad(1:n) = a
    Do m = 1, E
       If ( g1 - m >= 1 ) Then
          apad(1-m) = aglob(g1-m)
       Else
          apad(1-m) = a(1) - Real(m,Int64)*( a(2) - a(1) )
       End If
       If ( g1 + n - 1 + m <= ng ) Then
          apad(n+m) = aglob(g1+n-1+m)
       Else
          apad(n+m) = a(n) + Real(m,Int64)*( a(n) - a(n-1) )
       End If
    End Do

  End Subroutine pad_axis

  !> P(1-E:n1+E, :, 1-E:n3+E) = F, with the extra planes filled from the neighbouring rank / periodic partner
  !  / edge (see above). xface/zface: F is located at faces (not cell centres) in that direction. Needs at
  !  least E interior planes on every rank.
  Subroutine pad_field(F, n1, n2, n3, xface, zface, E, P)

    Integer(Int32), Intent(In)  :: n1, n2, n3, E
    Real   (Int64), Intent(In)  :: F(n1,n2,n3)
    Logical,        Intent(In)  :: xface, zface
    Real   (Int64), Intent(Out) :: P(1-E:n1+E,n2,1-E:n3+E)

    Logical        :: is_first, is_last, per
    Integer(Int32) :: up, down, partner, dst_dn, dst_up, src_dn, src_up, sdn, sup, m
    Real   (Int64), Allocatable :: sb(:,:,:), rb(:,:,:)

    P(1:n1,:,1:n3) = F

    !-- x -----------------------------------------------------------------
    per = ( x_bc_type == 0 )
    Call x_halo_neighbors(up, down)
    Call x_periodic_partner(is_first, is_last, partner)
    If ( is_first .And. is_last ) Then
       If ( per ) Then
          sdn = Merge(3, 4, xface);  sup = Merge(n1-2, n1-3, xface)
          Do m = 1, E
             P(n1+m,:,1:n3) = F(sdn+m-1,:,:)
             P(1-m, :,1:n3) = F(sup-m+1,:,:)
          End Do
       Else
          Do m = 1, E
             P(n1+m,:,1:n3) = F(n1,:,:)
             P(1-m, :,1:n3) = F(1, :,:)
          End Do
       End If
    Else
       dst_dn = down;  src_up = up;  dst_up = up;  src_dn = down
       sdn = 3;  sup = n1-2
       If ( per .And. is_first ) Then
          dst_dn = partner;  src_dn = partner;  sdn = Merge(3, 4, xface)
       End If
       If ( per .And. is_last ) Then
          dst_up = partner;  src_up = partner;  sup = Merge(n1-2, n1-3, xface)
       End If
       Allocate( sb(E,n2,n3), rb(E,n2,n3) )
       sb = F(sdn:sdn+E-1,:,:)
       Call Mpi_sendrecv( sb, E*n2*n3, Mpi_real8, dst_dn, 61, rb, E*n2*n3, Mpi_real8, src_up, 61, MPI_COMM_WORLD, istat, ierr )
       If ( src_up /= MPI_PROC_NULL ) Then
          P(n1+1:n1+E,:,1:n3) = rb
       Else
          Do m = 1, E
             P(n1+m,:,1:n3) = F(n1,:,:)
          End Do
       End If
       sb = F(sup-E+1:sup,:,:)
       Call Mpi_sendrecv( sb, E*n2*n3, Mpi_real8, dst_up, 62, rb, E*n2*n3, Mpi_real8, src_dn, 62, MPI_COMM_WORLD, istat, ierr )
       If ( src_dn /= MPI_PROC_NULL ) Then
          P(1-E:0,:,1:n3) = rb
       Else
          Do m = 1, E
             P(1-m,:,1:n3) = F(1,:,:)
          End Do
       End If
       Deallocate( sb, rb )
    End If

    !-- z (over the x-extended extent, so the corners are consistent) ---------
    per = ( z_bc_type == 0 )
    Call z_halo_neighbors(up, down)
    Call z_periodic_partner(is_first, is_last, partner)
    If ( is_first .And. is_last ) Then
       If ( per ) Then
          sdn = Merge(3, 4, zface);  sup = Merge(n3-2, n3-3, zface)
          Do m = 1, E
             P(:,:,n3+m) = P(:,:,sdn+m-1)
             P(:,:,1-m ) = P(:,:,sup-m+1)
          End Do
       Else
          Do m = 1, E
             P(:,:,n3+m) = P(:,:,n3)
             P(:,:,1-m ) = P(:,:,1)
          End Do
       End If
    Else
       dst_dn = down;  src_up = up;  dst_up = up;  src_dn = down
       sdn = 3;  sup = n3-2
       If ( per .And. is_first ) Then
          dst_dn = partner;  src_dn = partner;  sdn = Merge(3, 4, zface)
       End If
       If ( per .And. is_last ) Then
          dst_up = partner;  src_up = partner;  sup = Merge(n3-2, n3-3, zface)
       End If
       Allocate( sb(n1+2*E,n2,E), rb(n1+2*E,n2,E) )
       sb = P(:,:,sdn:sdn+E-1)
       Call Mpi_sendrecv( sb, (n1+2*E)*n2*E, Mpi_real8, dst_dn, 63, rb, (n1+2*E)*n2*E, Mpi_real8, src_up, 63, &
                          MPI_COMM_WORLD, istat, ierr )
       If ( src_up /= MPI_PROC_NULL ) Then
          P(:,:,n3+1:n3+E) = rb
       Else
          Do m = 1, E
             P(:,:,n3+m) = P(:,:,n3)
          End Do
       End If
       sb = P(:,:,sup-E+1:sup)
       Call Mpi_sendrecv( sb, (n1+2*E)*n2*E, Mpi_real8, dst_up, 64, rb, (n1+2*E)*n2*E, Mpi_real8, src_dn, 64, &
                          MPI_COMM_WORLD, istat, ierr )
       If ( src_dn /= MPI_PROC_NULL ) Then
          P(:,:,1-E:0) = rb
       Else
          Do m = 1, E
             P(:,:,1-m) = P(:,:,1)
          End Do
       End If
       Deallocate( sb, rb )
    End If

  End Subroutine pad_field

  !> Fe(1-E:n1+E, :, 1-E:n3+E) = F plus E extra planes in x and z from the neighbour rank / periodic partner /
  !  edge (the rules of halo_pad.pad_field). Device-aware (F and Fe present on the device): the copy and the
  !  pack/unpack run on the device in GPU builds, only the small exchanged planes go through the host for MPI.
  Subroutine pad_field_dev(F, Fe, n1, n2, n3, xface, zface, E)

    Integer(Int32), Intent(In)    :: n1, n2, n3, E
    Real   (Int64), Intent(In)    :: F(n1,n2,n3)
    Real   (Int64), Intent(InOut) :: Fe(1-E:n1+E, n2, 1-E:n3+E)
    Logical,        Intent(In)    :: xface, zface

    Call pad_field_dev_y(F, Fe, n1, n2, n3, xface, zface, E, 0)

  End Subroutine pad_field_dev


  !> pad_field_dev into an Fe that also has ye rows beyond y = 1 and y = n2 (left untouched)
  Subroutine pad_field_dev_y(F, Fe, n1, n2, n3, xface, zface, E, ye)

    Integer(Int32), Intent(In)    :: n1, n2, n3, E, ye
    Real   (Int64), Intent(In)    :: F(n1,n2,n3)
    Real   (Int64), Intent(InOut) :: Fe(1-E:n1+E, 1-ye:n2+ye, 1-E:n3+E)
    Logical,        Intent(In)    :: xface, zface

    Logical        :: is_first, is_last, per
    Integer(Int32) :: up, down, partner, dst_dn, dst_up, src_dn, src_up, sdn, sup, m, i, j, k, cnt

    cnt = Max((n1+2*E)*n2*E, E*n2*n3)
    Call pad_buffers(cnt)

    !$acc parallel loop collapse(3) present(F,Fe)
    Do k = 1, n3
       Do j = 1, n2
          Do i = 1, n1
             Fe(i,j,k) = F(i,j,k)
          End Do
       End Do
    End Do
    !$acc end parallel loop

    !-- x --
    per = ( x_bc_type == 0 )
    Call x_halo_neighbors(up, down)
    Call x_periodic_partner(is_first, is_last, partner)
    If ( is_first .And. is_last ) Then
       If ( per ) Then
          sdn = Merge(3, 4, xface);  sup = Merge(n1-2, n1-3, xface)
          !$acc parallel loop collapse(3) present(Fe)
          Do k = 1, n3
             Do j = 1, n2
                Do m = 1, E
                   Fe(n1+m,j,k) = Fe(sdn+m-1,j,k)
                   Fe(1-m, j,k) = Fe(sup-m+1,j,k)
                End Do
             End Do
          End Do
          !$acc end parallel loop
       Else
          !$acc parallel loop collapse(3) present(Fe)
          Do k = 1, n3
             Do j = 1, n2
                Do m = 1, E
                   Fe(n1+m,j,k) = Fe(n1,j,k)
                   Fe(1-m, j,k) = Fe(1, j,k)
                End Do
             End Do
          End Do
          !$acc end parallel loop
       End If
    Else
       dst_dn = down;  src_up = up;  dst_up = up;  src_dn = down
       sdn = 3;  sup = n1-2
       If ( per .And. is_first ) Then
          dst_dn = partner;  src_dn = partner;  sdn = Merge(3, 4, xface)
       End If
       If ( per .And. is_last ) Then
          dst_up = partner;  src_up = partner;  sup = Merge(n1-2, n1-3, xface)
       End If
       cnt = E*n2*n3
       !$acc parallel loop collapse(3) present(Fe,pad_bs)
       Do k = 1, n3
          Do j = 1, n2
             Do m = 1, E
                pad_bs(m + E*((j-1) + n2*(k-1))) = Fe(sdn+m-1,j,k)
             End Do
          End Do
       End Do
       !$acc end parallel loop
       !$acc update host(pad_bs(1:cnt))
       Call MPI_Sendrecv( pad_bs, cnt, MPI_real8, dst_dn, 71, pad_br, cnt, MPI_real8, src_up, 71, MPI_COMM_WORLD, istat, ierr )
       !$acc update device(pad_br(1:cnt))
       !$acc parallel loop collapse(3) present(Fe,pad_br)
       Do k = 1, n3
          Do j = 1, n2
             Do m = 1, E
                If ( src_up /= MPI_PROC_NULL ) Then
                   Fe(n1+m,j,k) = pad_br(m + E*((j-1) + n2*(k-1)))
                Else
                   Fe(n1+m,j,k) = Fe(n1,j,k)
                End If
             End Do
          End Do
       End Do
       !$acc end parallel loop
       !$acc parallel loop collapse(3) present(Fe,pad_bs)
       Do k = 1, n3
          Do j = 1, n2
             Do m = 1, E
                pad_bs(m + E*((j-1) + n2*(k-1))) = Fe(sup-E+m,j,k)
             End Do
          End Do
       End Do
       !$acc end parallel loop
       !$acc update host(pad_bs(1:cnt))
       Call MPI_Sendrecv( pad_bs, cnt, MPI_real8, dst_up, 72, pad_br, cnt, MPI_real8, src_dn, 72, MPI_COMM_WORLD, istat, ierr )
       !$acc update device(pad_br(1:cnt))
       !$acc parallel loop collapse(3) present(Fe,pad_br)
       Do k = 1, n3
          Do j = 1, n2
             Do m = 1, E
                If ( src_dn /= MPI_PROC_NULL ) Then
                   Fe(m-E,j,k) = pad_br(m + E*((j-1) + n2*(k-1)))
                Else
                   Fe(m-E,j,k) = Fe(1,j,k)
                End If
             End Do
          End Do
       End Do
       !$acc end parallel loop
    End If

    !-- z (over the x-extended extent) --
    per = ( z_bc_type == 0 )
    Call z_halo_neighbors(up, down)
    Call z_periodic_partner(is_first, is_last, partner)
    If ( is_first .And. is_last ) Then
       If ( per ) Then
          sdn = Merge(3, 4, zface);  sup = Merge(n3-2, n3-3, zface)
          !$acc parallel loop collapse(3) present(Fe)
          Do m = 1, E
             Do j = 1, n2
                Do i = 1-E, n1+E
                   Fe(i,j,n3+m) = Fe(i,j,sdn+m-1)
                   Fe(i,j,1-m ) = Fe(i,j,sup-m+1)
                End Do
             End Do
          End Do
          !$acc end parallel loop
       Else
          !$acc parallel loop collapse(3) present(Fe)
          Do m = 1, E
             Do j = 1, n2
                Do i = 1-E, n1+E
                   Fe(i,j,n3+m) = Fe(i,j,n3)
                   Fe(i,j,1-m ) = Fe(i,j,1)
                End Do
             End Do
          End Do
          !$acc end parallel loop
       End If
    Else
       dst_dn = down;  src_up = up;  dst_up = up;  src_dn = down
       sdn = 3;  sup = n3-2
       If ( per .And. is_first ) Then
          dst_dn = partner;  src_dn = partner;  sdn = Merge(3, 4, zface)
       End If
       If ( per .And. is_last ) Then
          dst_up = partner;  src_up = partner;  sup = Merge(n3-2, n3-3, zface)
       End If
       cnt = (n1+2*E)*n2*E
       !$acc parallel loop collapse(3) present(Fe,pad_bs)
       Do m = 1, E
          Do j = 1, n2
             Do i = 1, n1+2*E
                pad_bs(i + (n1+2*E)*((j-1) + n2*(m-1))) = Fe(i-E,j,sdn+m-1)
             End Do
          End Do
       End Do
       !$acc end parallel loop
       !$acc update host(pad_bs(1:cnt))
       Call MPI_Sendrecv( pad_bs, cnt, MPI_real8, dst_dn, 73, pad_br, cnt, MPI_real8, src_up, 73, MPI_COMM_WORLD, istat, ierr )
       !$acc update device(pad_br(1:cnt))
       !$acc parallel loop collapse(3) present(Fe,pad_br)
       Do m = 1, E
          Do j = 1, n2
             Do i = 1, n1+2*E
                If ( src_up /= MPI_PROC_NULL ) Then
                   Fe(i-E,j,n3+m) = pad_br(i + (n1+2*E)*((j-1) + n2*(m-1)))
                Else
                   Fe(i-E,j,n3+m) = Fe(i-E,j,n3)
                End If
             End Do
          End Do
       End Do
       !$acc end parallel loop
       !$acc parallel loop collapse(3) present(Fe,pad_bs)
       Do m = 1, E
          Do j = 1, n2
             Do i = 1, n1+2*E
                pad_bs(i + (n1+2*E)*((j-1) + n2*(m-1))) = Fe(i-E,j,sup-E+m)
             End Do
          End Do
       End Do
       !$acc end parallel loop
       !$acc update host(pad_bs(1:cnt))
       Call MPI_Sendrecv( pad_bs, cnt, MPI_real8, dst_up, 74, pad_br, cnt, MPI_real8, src_dn, 74, MPI_COMM_WORLD, istat, ierr )
       !$acc update device(pad_br(1:cnt))
       !$acc parallel loop collapse(3) present(Fe,pad_br)
       Do m = 1, E
          Do j = 1, n2
             Do i = 1, n1+2*E
                If ( src_dn /= MPI_PROC_NULL ) Then
                   Fe(i-E,j,m-E) = pad_br(i + (n1+2*E)*((j-1) + n2*(m-1)))
                Else
                   Fe(i-E,j,m-E) = Fe(i-E,j,1)
                End If
             End Do
          End Do
       End Do
       !$acc end parallel loop
    End If

  End Subroutine pad_field_dev_y

End Module halo_pad
