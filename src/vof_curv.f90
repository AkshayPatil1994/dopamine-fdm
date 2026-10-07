!> Interface curvature of the VOF field C for surface tension: height functions ((2 NH + 1)-cell columns, 3x3 columns, Popinet 2009) in
!  the coordinate direction closest to the interface normal; cells whose columns do not terminate in a full and an empty cell
!  (thin ligaments, drops of a few cells) take the mean of the valid height-function curvatures of their neighbours, else zero.
!  Convention: kappa = div(n_out) with n_out pointing out of the liquid, i.e. 2/R for a spherical drop; the force on the fluid
!  is sigma kappa grad(C).
Module vof_curv

  Use iso_fortran_env, Only : Int32, Int64
  Use global
  Use halo_pad, Only : pad_field
  Use vof_plic, Only : vof_eps
  Use vof_state, Only : Cv, hy, hz, vof_hf_fail
  Use vof_pressure, Only : vp_halo

  Implicit None

  Integer(Int32), Parameter :: NH = 4        ! half length of the height-function columns (cells)
  Integer(Int32), Parameter :: EH = NH - 1   ! planes beyond the ghost plane the columns reach from an interior cell

Contains

  !> kap(1:nxg,1:nyg,1:nzg) for the interface cells (0 < C < 1) of the padded field Cv, halo filled
  Subroutine vof_curvature(kap)

    Real(Int64), Intent(Out) :: kap(nxg,nyg,nzg)

    Real(Int64), Allocatable :: P(:,:,:), tmp(:,:,:), cw(:,:,:), kap0(:,:,:)
    Logical, Allocatable :: valid(:,:,:)
    Integer(Int32) :: i, j, k, m
    Real(Int64) :: c, kv, acc
    Integer(Int32) :: ii, jj, kk, cnt
    Logical :: ok

    Allocate( P(1-EH:nxg+EH,1-EH:nyg+EH,1-EH:nzg+EH), tmp(1-EH:nxg+EH,nyg,1-EH:nzg+EH), cw(nxg,nyg,nzg) )
    Allocate( valid(nxg,nyg,nzg), kap0(nxg,nyg,nzg) )
    cw = Cv(1:nxg,1:nyg,1:nzg)
    Call pad_field(cw, nxg, nyg, nzg, .False., .False., EH, tmp)
    P = 0d0
    P(:,1:nyg,:) = tmp
    Do m = 1, EH
       P(:,1-m,:)   = P(:,2+m,:)
       P(:,nyg+m,:) = P(:,nyg-1-m,:)
    End Do
    Deallocate( tmp, cw )

    kap = 0d0;  valid = .False.;  vof_hf_fail = 0
    Do k = 2, nzg-1
       Do j = 2, nyg-1
          Do i = 2, nxg-1
             c = P(i,j,k)
             If ( c <= vof_eps .Or. c >= 1d0 - vof_eps ) Cycle
             Call height_curvature(P, i, j, k, kv, ok)
             If ( ok ) Then
                kap(i,j,k) = kv;  valid(i,j,k) = .True.
             Else
                vof_hf_fail = vof_hf_fail + 1
             End If
          End Do
       End Do
    End Do
    Call vp_halo(kap, .False.)
    kap0 = kap

    ! interface cells without a height function: mean of the valid neighbours (the halo carries the seam cells)
    Do k = 2, nzg-1
       Do j = 2, nyg-1
          Do i = 2, nxg-1
             c = P(i,j,k)
             If ( c <= vof_eps .Or. c >= 1d0 - vof_eps .Or. valid(i,j,k) ) Cycle
             acc = 0d0;  cnt = 0
             Do kk = k-1, k+1
                Do jj = j-1, j+1
                   Do ii = i-1, i+1
                      If ( P(ii,jj,kk) > vof_eps .And. P(ii,jj,kk) < 1d0 - vof_eps .And. kap0(ii,jj,kk) /= 0d0 ) Then
                         acc = acc + kap0(ii,jj,kk);  cnt = cnt + 1
                      End If
                   End Do
                End Do
             End Do
             If ( cnt > 0 ) kap(i,j,k) = acc/cnt
          End Do
       End Do
    End Do
    Call vp_halo(kap, .False.)
    Deallocate( P, valid, kap0 )

  End Subroutine vof_curvature


  !> Height-function curvature of interface cell (i,j,k): direction d of the largest central difference of C, heights of the
  !  3x3 neighbouring columns, valid only if every column starts in a full (or empty) and ends in an empty (full) cell
  Subroutine height_curvature(P, i, j, k, kv, ok)

    Real(Int64),    Intent(In)  :: P(1-EH:,1-EH:,1-EH:)
    Integer(Int32), Intent(In)  :: i, j, k
    Real(Int64),    Intent(Out) :: kv
    Logical,        Intent(Out) :: ok

    Integer(Int32) :: d, t1, t2, a, b, m, q(3), e(3,3)
    Real(Int64) :: g(3), H(-1:1,-1:1), h1, h2, lo, hi, Hx, Hz, Hxx, Hzz, Hxz, den
    Logical :: below

    e = 0
    e(1,1) = 1;  e(2,2) = 1;  e(3,3) = 1
    g(1) = P(i+1,j,k) - P(i-1,j,k)
    g(2) = P(i,j+1,k) - P(i,j-1,k)
    g(3) = P(i,j,k+1) - P(i,j,k-1)
    kv = 0d0;  ok = .False.
    d = MaxLoc(Abs(g), 1)
    If ( Abs(g(d)) < 1d-12 ) Return
    below = g(d) < 0d0               ! liquid on the low side along d
    t1 = Mod(d, 3) + 1;  t2 = Mod(d+1, 3) + 1

    Do b = -1, 1
       Do a = -1, 1
          H(a,b) = 0d0
          Do m = -NH, NH
             q = (/ i, j, k /) + a*e(:,t1) + b*e(:,t2) + m*e(:,d)
             H(a,b) = H(a,b) + P(q(1),q(2),q(3))*cell_size(d, q)
          End Do
          q = (/ i, j, k /) + a*e(:,t1) + b*e(:,t2) - NH*e(:,d)
          lo = P(q(1),q(2),q(3))
          q = (/ i, j, k /) + a*e(:,t1) + b*e(:,t2) + NH*e(:,d)
          hi = P(q(1),q(2),q(3))
          If ( below ) Then
             If ( lo < 1d0 - 1d-9 .Or. hi > 1d-9 ) Return
          Else
             If ( lo > 1d-9 .Or. hi < 1d0 - 1d-9 ) Return
          End If
       End Do
    End Do

    q = (/ i, j, k /)
    h1 = cell_size(t1, q);  h2 = cell_size(t2, q)
    Hx  = ( H(1,0) - H(-1,0) )/(2d0*h1)
    Hz  = ( H(0,1) - H(0,-1) )/(2d0*h2)
    Hxx = ( H(1,0) - 2d0*H(0,0) + H(-1,0) )/h1**2
    Hzz = ( H(0,1) - 2d0*H(0,0) + H(0,-1) )/h2**2
    Hxz = ( H(1,1) - H(1,-1) - H(-1,1) + H(-1,-1) )/(4d0*h1*h2)
    den = (1d0 + Hx**2 + Hz**2)**1.5d0
    kv = -( Hxx*(1d0 + Hz**2) + Hzz*(1d0 + Hx**2) - 2d0*Hxz*Hx*Hz )/den
    ok = .True.

  End Subroutine height_curvature


  !> Size of the cell q(1:3) along direction d (the padded planes repeat the edge sizes)
  Pure Function cell_size(d, q) Result(h)

    Integer(Int32), Intent(In) :: d, q(3)
    Real(Int64) :: h

    If ( d == 1 ) Then
       h = dx
    Else If ( d == 2 ) Then
       h = hy(Min(Max(q(2), 1), nyg))
    Else
       h = hz(Min(Max(q(3), 1), nzg))
    End If

  End Function cell_size

End Module vof_curv
