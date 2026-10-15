!> Upwind-biased dissipation for the momentum convection (advection_scheme = 2).
!  The upwind-biased face value is phi_central + s*D, s the sign of the convecting velocity, D = third difference/12 (adv_order = 3)
!  or -fifth difference/60 (adv_order = 5). The skew-symmetric second-order central convection of compute_rhs_* is kept and only
!  beta*|uc|*D is added as a conservative flux difference, so beta = 0 is exactly advection_scheme = 0 and beta = 1 is the
!  dissipation of UP3/UP5 on top of the second-order central part (the formal order stays 2; the added error is O(dx^3), O(dx^5)).
!  For constant u it is -beta*|u|*dx^3/12 d4phi/dx4 (order 5: -beta*|u|*dx^5/60 d6phi/dx6): symbol (4/3) beta |u|/dx sin^4(theta/2)
!  resp. (16/15) beta |u|/dx sin^6(theta/2), so only modes near the grid scale feel it.
!  beta = max(adv_beta0, adv_beta_ibm in the adv_band cells around a solid cell): the damping sits where the sharp velocity jump
!  of a body excites grid-scale ringing, and is small elsewhere. Rai & Moin (1991) use this construction.
Module adv_hybrid

  Use iso_fortran_env, Only : Int32, Int64
  Use global
  Use mpi
  Use halo_pad, Only : pad_field, pad_field_dev

  Implicit None

  Integer(Int32), Parameter :: EA = 3   ! halo planes in x and z (5th-order stencil reaches 3 cells beyond the face)
  Real(Int64), Allocatable, Dimension(:,:,:) :: Up, Vp, Wp          ! velocities with EA extra planes in x and z
  Real(Int64), Allocatable, Dimension(:,:,:) :: adv_bu, adv_bv, adv_bw   ! beta at the U, V, W points

Contains

  !> beta fields: adv_beta0 everywhere, adv_beta_ibm within adv_band cells of a solid cell (phi < 0)
  Subroutine adv_init

    Real(Int64), Allocatable :: sd(:,:,:), sp(:,:,:), sn(:,:,:)
    Integer(Int32) :: i, j, k, a, b, c, E

    Allocate( Up(1-EA:nx+EA,nyg,1-EA:nzg+EA), Vp(1-EA:nxg+EA,ny,1-EA:nzg+EA), Wp(1-EA:nxg+EA,nyg,1-EA:nz+EA) )
    Allocate( adv_bu(nx,nyg,nzg), adv_bv(nxg,ny,nzg), adv_bw(nxg,nyg,nz) )
    Up = 0d0;  Vp = 0d0;  Wp = 0d0
    ! solid indicator on the cells, padded by adv_band+1 planes in x and z from the neighbouring ranks: the face betas read one cell
    ! beyond the ghost plane (U has nx = nxg points on all but the last rank), and the box maximum needs adv_band more
    E = adv_band + 1
    Allocate( sp(1-E:nxg+E,nyg,1-E:nzg+E), sd(1:nxg+1,nyg,1:nzg+1) )
    sd = 0d0
    If ( ibm_input_mode >= 1 .And. Allocated(phi) ) Then
       Allocate( sn(nxg,nyg,nzg) )
       sn = 0d0
       Where ( phi < 0d0 ) sn = 1d0
       Call pad_field(sn, nxg, nyg, nzg, .False., .False., E, sp)
       Deallocate( sn )
       Do k = 1, nzg+1
          Do j = 1, nyg
             Do i = 1, nxg+1
                Do c = -adv_band, adv_band
                   Do b = -adv_band, adv_band
                      If ( j+b < 1 .Or. j+b > nyg ) Cycle
                      Do a = -adv_band, adv_band
                         sd(i,j,k) = Max( sd(i,j,k), sp(i+a,j+b,k+c) )
                      End Do
                   End Do
                End Do
             End Do
          End Do
       End Do
    End If
    adv_bu = adv_beta0;  adv_bv = adv_beta0;  adv_bw = adv_beta0
    Do k = 1, nzg
       Do j = 1, nyg
          Do i = 1, nx
             adv_bu(i,j,k) = Max( adv_beta0, adv_beta_ibm*Max( sd(i,j,k), sd(i+1,j,k) ) )
          End Do
       End Do
    End Do
    Do k = 1, nzg
       Do j = 1, ny
          Do i = 1, nxg
             adv_bv(i,j,k) = Max( adv_beta0, adv_beta_ibm*Max( sd(i,j,k), sd(i,j+1,k) ) )
          End Do
       End Do
    End Do
    Do k = 1, nz
       Do j = 1, nyg
          Do i = 1, nxg
             adv_bw(i,j,k) = Max( adv_beta0, adv_beta_ibm*Max( sd(i,j,k), sd(i,j,k+1) ) )
          End Do
       End Do
    End Do
    Deallocate( sd, sp )
    !$acc enter data copyin(adv_bu,adv_bv,adv_bw) create(Up,Vp,Wp)

  End Subroutine adv_init


  !> Padded copies of the velocity, once per RK stage before the compute_rhs_* calls
  Subroutine adv_pad(U_, V_, W_)

    Real(Int64), Intent(In) :: U_(nx,nyg,nzg), V_(nxg,ny,nzg), W_(nxg,nyg,nz)

    Call pad_field_dev(U_, Up, nx, nyg, nzg, .True., .False., EA)
    Call pad_field_dev(V_, Vp, nxg, ny, nzg, .False., .False., EA)
    Call pad_field_dev(W_, Wp, nxg, nyg, nz, .False., .True., EA)

  End Subroutine adv_pad


  !> Stencil weights of the dissipative flux for the order o (3, 5; any other value: none), for s(-2:3) = a..f
  Subroutine adv_wsel(o, w)
    !$acc routine seq

    Integer(Int32), Intent(In)  :: o
    Real(Int64),    Intent(Out) :: w(6)

    If ( o == 5 ) Then
       w = (/ 1d0, -5d0, 10d0, -10d0, 5d0, -1d0 /)/60d0
    Else If ( o == 3 ) Then
       w = (/ 0d0, -1d0, 3d0, -3d0, 1d0, 0d0 /)/12d0
    Else
       w = 0d0
    End If

  End Subroutine adv_wsel


  !> Dissipative flux through a face between c = s(0) and d = s(1) of the advected variable: beta*|uc|*sum(w*s). Order 3:
  !  (s(2) - 3 s(1) + 3 s(0) - s(-1))/12 (UP3 minus C4);
  !  order 5: -(s(3) - 5 s(2) + 10 s(1) - 10 s(0) + 5 s(-1) - s(-2))/60 (UP5 minus C6)
  Pure Function adv_dflux(a, b, c, d, e, f, uc, beta, w) Result(r)
    !$acc routine seq

    Real(Int64), Intent(In) :: a, b, c, d, e, f, uc, beta, w(6)
    Real(Int64) :: r

    r = beta*Abs(uc)*( w(1)*a + w(2)*b + w(3)*c + w(4)*d + w(5)*e + w(6)*f )

  End Function adv_dflux


  !> Order usable for the face between rows m and m+1 of 1..n (no halo in y): reduced next to the walls
  Pure Function adv_ord_y(m, n, ord) Result(r)
    !$acc routine seq

    Integer(Int32), Intent(In) :: m, n, ord
    Integer(Int32) :: r

    r = 0
    If ( ord == 5 .And. m-2 >= 1 .And. m+3 <= n ) Then
       r = 5
    Else If ( m-1 >= 1 .And. m+2 <= n ) Then
       r = 3
    End If

  End Function adv_ord_y


  !> rhs_u += -d/dx_j of the dissipative fluxes
  Subroutine adv_add_u(rhs_u)

    Real(Int64), Intent(InOut) :: rhs_u(2:nx-1,2:nyg-1,2:nzg-1)
    Integer(Int32) :: i, j, k, m, ii, jj, kk, oy, n, p, q
    Real(Int64) :: r, uc, bt, sg, wx(6), wy(6), dyj, dzk, inv_dx

    n = adv_order
    inv_dx = 1d0/dx
    Call adv_wsel(n, wx)
    !$acc parallel loop collapse(3) present(rhs_u,Up,Vp,Wp,adv_bu,y,z) &
    !$acc& private(r,m,ii,jj,kk,uc,bt,sg,oy,p,q,wy,dyj,dzk) firstprivate(wx,n,inv_dx)
    Do k = 2, nzg-1
       Do j = 2, nyg-1
          Do i = 2, nx-1
             dyj = 1d0/( y(j) - y(j-1) );  dzk = 1d0/( z(k) - z(k-1) )
             r = 0d0
             Do m = 0, 1
                sg = 2d0*m - 1d0
                ii = i - 1 + m
                uc = 0.5d0*( Up(ii,j,k) + Up(ii+1,j,k) )
                bt = Max( adv_bu(ii,j,k), adv_bu(ii+1,j,k) )
                r = r - sg*adv_dflux( Up(ii-2,j,k), Up(ii-1,j,k), Up(ii,j,k), Up(ii+1,j,k), Up(ii+2,j,k), &
                                      Up(ii+3,j,k), uc, bt, wx )*inv_dx
                jj = j - 1 + m
                uc = 0.5d0*( Vp(i,jj,k) + Vp(i+1,jj,k) )
                bt = Max( adv_bu(i,jj,k), adv_bu(i,jj+1,k) )
                oy = adv_ord_y(jj, nyg, n)
                Call adv_wsel(oy, wy)
                r = r - sg*adv_dflux( Up(i,Max(jj-2,1),k), Up(i,Max(jj-1,1),k), Up(i,jj,k), Up(i,jj+1,k), &
                                      Up(i,Min(jj+2,nyg),k), Up(i,Min(jj+3,nyg),k), uc, bt, wy )*dyj
                kk = k - 1 + m
                uc = 0.5d0*( Wp(i,j,kk) + Wp(i+1,j,kk) )
                bt = Max( adv_bu(i,j,kk), adv_bu(i,j,kk+1) )
                r = r - sg*adv_dflux( Up(i,j,kk-2), Up(i,j,kk-1), Up(i,j,kk), Up(i,j,kk+1), Up(i,j,kk+2), &
                                      Up(i,j,kk+3), uc, bt, wx )*dzk
             End Do
             rhs_u(i,j,k) = rhs_u(i,j,k) + r
          End Do
       End Do
    End Do
    !$acc end parallel loop

  End Subroutine adv_add_u

  !> rhs_v += -d/dx_j of the dissipative fluxes
  Subroutine adv_add_v(rhs_v)

    Real(Int64), Intent(InOut) :: rhs_v(2:nxg-1,2:ny-1,2:nzg-1)
    Integer(Int32) :: i, j, k, m, ii, jj, kk, oy, n, p, q
    Real(Int64) :: r, uc, bt, sg, wx(6), wy(6), dyj, dzk, inv_dx

    n = adv_order
    inv_dx = 1d0/dx
    Call adv_wsel(n, wx)
    !$acc parallel loop collapse(3) present(rhs_v,Up,Vp,Wp,adv_bv,yg,z) &
    !$acc& private(r,m,ii,jj,kk,uc,bt,sg,oy,p,q,wy,dyj,dzk) firstprivate(wx,n,inv_dx)
    Do k = 2, nzg-1
       Do j = 2, ny-1
          Do i = 2, nxg-1
             dyj = 1d0/( yg(j+1) - yg(j) );  dzk = 1d0/( z(k) - z(k-1) )
             r = 0d0
             Do m = 0, 1
                sg = 2d0*m - 1d0
                ii = i - 1 + m
                uc = 0.5d0*( Up(ii,j,k) + Up(ii,j+1,k) )
                bt = Max( adv_bv(ii,j,k), adv_bv(ii+1,j,k) )
                r = r - sg*adv_dflux( Vp(ii-2,j,k), Vp(ii-1,j,k), Vp(ii,j,k), Vp(ii+1,j,k), Vp(ii+2,j,k), &
                                      Vp(ii+3,j,k), uc, bt, wx )*inv_dx
                jj = j - 1 + m
                uc = 0.5d0*( Vp(i,jj,k) + Vp(i,jj+1,k) )
                bt = Max( adv_bv(i,jj,k), adv_bv(i,jj+1,k) )
                oy = adv_ord_y(jj, ny, n)
                Call adv_wsel(oy, wy)
                r = r - sg*adv_dflux( Vp(i,Max(jj-2,1),k), Vp(i,Max(jj-1,1),k), Vp(i,jj,k), Vp(i,jj+1,k), &
                                      Vp(i,Min(jj+2,ny),k), Vp(i,Min(jj+3,ny),k), uc, bt, wy )*dyj
                kk = k - 1 + m
                uc = 0.5d0*( Wp(i,j,kk) + Wp(i,j+1,kk) )
                bt = Max( adv_bv(i,j,kk), adv_bv(i,j,kk+1) )
                r = r - sg*adv_dflux( Vp(i,j,kk-2), Vp(i,j,kk-1), Vp(i,j,kk), Vp(i,j,kk+1), Vp(i,j,kk+2), &
                                      Vp(i,j,kk+3), uc, bt, wx )*dzk
             End Do
             rhs_v(i,j,k) = rhs_v(i,j,k) + r
          End Do
       End Do
    End Do
    !$acc end parallel loop

  End Subroutine adv_add_v

  !> rhs_w += -d/dx_j of the dissipative fluxes
  Subroutine adv_add_w(rhs_w)

    Real(Int64), Intent(InOut) :: rhs_w(2:nxg-1,2:nyg-1,2:nz-1)
    Integer(Int32) :: i, j, k, m, ii, jj, kk, oy, n, p, q
    Real(Int64) :: r, uc, bt, sg, wx(6), wy(6), dyj, dzk, inv_dx

    n = adv_order
    inv_dx = 1d0/dx
    Call adv_wsel(n, wx)
    !$acc parallel loop collapse(3) present(rhs_w,Up,Vp,Wp,adv_bw,y,zg) &
    !$acc& private(r,m,ii,jj,kk,uc,bt,sg,oy,p,q,wy,dyj,dzk) firstprivate(wx,n,inv_dx)
    Do k = 2, nz-1
       Do j = 2, nyg-1
          Do i = 2, nxg-1
             dyj = 1d0/( y(j) - y(j-1) );  dzk = 1d0/( zg(k+1) - zg(k) )
             r = 0d0
             Do m = 0, 1
                sg = 2d0*m - 1d0
                ii = i - 1 + m
                uc = 0.5d0*( Up(ii,j,k) + Up(ii,j,k+1) )
                bt = Max( adv_bw(ii,j,k), adv_bw(ii+1,j,k) )
                r = r - sg*adv_dflux( Wp(ii-2,j,k), Wp(ii-1,j,k), Wp(ii,j,k), Wp(ii+1,j,k), Wp(ii+2,j,k), &
                                      Wp(ii+3,j,k), uc, bt, wx )*inv_dx
                jj = j - 1 + m
                uc = 0.5d0*( Vp(i,jj,k) + Vp(i,jj,k+1) )
                bt = Max( adv_bw(i,jj,k), adv_bw(i,jj+1,k) )
                oy = adv_ord_y(jj, nyg, n)
                Call adv_wsel(oy, wy)
                r = r - sg*adv_dflux( Wp(i,Max(jj-2,1),k), Wp(i,Max(jj-1,1),k), Wp(i,jj,k), Wp(i,jj+1,k), &
                                      Wp(i,Min(jj+2,nyg),k), Wp(i,Min(jj+3,nyg),k), uc, bt, wy )*dyj
                kk = k - 1 + m
                uc = 0.5d0*( Wp(i,j,kk) + Wp(i,j,kk+1) )
                bt = Max( adv_bw(i,j,kk), adv_bw(i,j,kk+1) )
                r = r - sg*adv_dflux( Wp(i,j,kk-2), Wp(i,j,kk-1), Wp(i,j,kk), Wp(i,j,kk+1), Wp(i,j,kk+2), &
                                      Wp(i,j,kk+3), uc, bt, wx )*dzk
             End Do
             rhs_w(i,j,k) = rhs_w(i,j,k) + r
          End Do
       End Do
    End Do
    !$acc end parallel loop

  End Subroutine adv_add_w

End Module adv_hybrid
