!> Host-side padded copies of a decomposed field: E extra planes beyond the one ghost plane the solver arrays carry, in x and z.
!  The single ghost plane is enough for the stencils of the solver itself, but not for anything that samples the field at
!  displaced points (a particle's RK stages, the image points of ghost-cell IBM, finite differences of the SDF); with only the
!  ghost plane such a sample near a rank seam clamps, and the result depends on the rank layout.
!
!  Plane rules follow apply_periodic_bc_x/z. Across a rank seam the extra planes come from the neighbour: the high side from its
!  planes 3..2+E, the low side from its planes n-1-E..n-2 (n = the neighbour's own plane count). Across the periodic wrap the
!  first rank supplies planes 4..3+E (cell-centred) or 3..2+E (face-located) and the last rank planes n-2-E..n-3 (cell-centred)
!  or n-1-E..n-2 (face-located). At a non-periodic domain edge the extra planes repeat the edge plane.
Module halo_pad

  Use iso_fortran_env, Only : Int32, Int64
  Use global
  Use mpi
  Use decomp, Only : x_halo_neighbors, z_halo_neighbors, x_periodic_partner, z_periodic_partner

  Implicit None

Contains

  !> Padded axis apad(1-E:n+E) from the local axis a(1:n): the neighbouring points from the global axis where they exist,
  !  otherwise extrapolated with the edge spacing (a periodic edge; the axes involved there are uniform)
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

  !> P(1-E:n1+E, :, 1-E:n3+E) = F, with the extra planes filled from the neighbouring rank / periodic partner / edge (see above).
  !  xface/zface: F is located at faces (not cell centres) in that direction. Needs at least E interior planes on every rank.
  Subroutine pad_field(F, n1, n2, n3, xface, zface, E, P)

    Integer(Int32), Intent(In)  :: n1, n2, n3, E
    Real   (Int64), Intent(In)  :: F(n1,n2,n3)
    Logical,        Intent(In)  :: xface, zface
    Real   (Int64), Intent(Out) :: P(1-E:n1+E,n2,1-E:n3+E)

    Logical        :: is_first, is_last, per
    Integer(Int32) :: up, down, partner, dst_dn, dst_up, src_dn, src_up, sdn, sup, m, tagbase
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

End Module halo_pad
