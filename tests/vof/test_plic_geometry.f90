!> Unit test of src/vof_plic.f90: alpha(C) round trip, box volume vs brute-force quadrature, sub-box volumes
Program test_plic_geometry

  Use iso_fortran_env, Only : Int32, Int64
  Use vof_plic

  Implicit None

  Integer(Int32), Parameter :: ntest = 40, nq = 160
  Integer(Int32) :: it, i, j, k, ndeg
  Real(Int64) :: m(3), c, alpha, err, errmax_rt, errmax_q, errmax_sb, f, fq, xlo(3), w(3), xx, yy, zz, cnt, r(8)
  Logical :: ok

  ok = .True.
  errmax_rt = 0d0;  errmax_q = 0d0;  errmax_sb = 0d0
  Call Random_seed()

  ! round trip: volume of the reconstructed plane reproduces C, including 2-D/1-D (zero-component) normals
  Do ndeg = 0, 2
     Do it = 1, 20000
        Call Random_number(r)
        m = 2d0*r(1:3) - 1d0
        If ( ndeg >= 1 ) m(3) = 0d0
        If ( ndeg >= 2 ) m(2) = 0d0
        If ( Abs(m(1)) + Abs(m(2)) + Abs(m(3)) < 1d-3 ) m(1) = 1d0
        c = Max(1d-9, Min(1d0 - 1d-9, r(4)))
        alpha = plic_alpha(c, m(1), m(2), m(3))
        f = plic_subbox_fraction(m(1)/Sum(Abs(m)), m(2)/Sum(Abs(m)), m(3)/Sum(Abs(m)), alpha, &
                                 -0.5d0, -0.5d0, -0.5d0, 1d0, 1d0, 1d0)
        errmax_rt = Max(errmax_rt, Abs(f - c))
     End Do
  End Do

  ! quadrature check of full-cell and sub-box fractions
  Do it = 1, ntest
     Call Random_number(r)
     m = 2d0*r(1:3) - 1d0
     m = m/Sum(Abs(m))
     alpha = (r(4) - 0.5d0)*1.2d0
     xlo = -0.5d0 + 0.5d0*r(5:7)
     w = Min(1d0 - (xlo + 0.5d0), 0.3d0 + 0.7d0*r(8))
     w = Max(w, 0.05d0)
     w = Min(w, 0.5d0 - xlo)

     cnt = 0d0
     Do k = 1, nq
        zz = xlo(3) + w(3)*(Real(k,Int64) - 0.5d0)/nq
        Do j = 1, nq
           yy = xlo(2) + w(2)*(Real(j,Int64) - 0.5d0)/nq
           Do i = 1, nq
              xx = xlo(1) + w(1)*(Real(i,Int64) - 0.5d0)/nq
              If ( m(1)*xx + m(2)*yy + m(3)*zz <= alpha ) cnt = cnt + 1d0
           End Do
        End Do
     End Do
     fq = cnt/Real(nq,Int64)**3
     f = plic_subbox_fraction(m(1), m(2), m(3), alpha, xlo(1), xlo(2), xlo(3), w(1), w(2), w(3))
     err = Abs(f - fq)
     errmax_sb = Max(errmax_sb, err)
  End Do

  Write(*,'(a,es10.3)') 'alpha round-trip max |V(alpha(C))-C|   : ', errmax_rt
  Write(*,'(a,es10.3)') 'sub-box fraction vs quadrature max err   : ', errmax_sb
  If ( errmax_rt > 1d-11 ) ok = .False.
  If ( errmax_sb > 2d-3 ) ok = .False.
  If ( ok ) Then
     Write(*,'(a)') 'PASS'
  Else
     Write(*,'(a)') 'FAIL'
     Stop 1
  End If

End Program test_plic_geometry
