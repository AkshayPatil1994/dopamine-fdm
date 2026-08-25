!> Standalone verification driver for src/sem.f90's read_mean_profile_TI (sem_profile_format=1); built via the CMake target verify_ti_profile (see CMakeLists.txt), not part of the dopamine executable
Program verify_ti_profile

  Use iso_fortran_env, Only : Int32, Int64
  Use mpi
  Use global
  Use synthetic_eddy_method

  Implicit None

  Real(Int64), Parameter :: tol = 1d-10
  Logical :: all_pass

  Call Mpi_Init(ierr)
  Call Mpi_Comm_rank(MPI_COMM_WORLD, myid, ierr)
  Call Mpi_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

  ny_global = 3; nz_global = 3
  Allocate( y_global(ny_global), z_global(nz_global) )
  y_global = (/ 0d0, 0.5d0, 1d0 /)
  z_global = (/ 0d0, 0.5d0, 1d0 /)

  sem_length_scale  = 0.01d0
  sem_Lscale_ratio_y = 0.3d0
  sem_Lscale_ratio_z = 0.2d0

  all_pass = .True.

  Call case_full_length_scales
  Call case_partial_length_scales
  Call case_no_length_scales
  Call case_tab_separated_header

  If ( myid == 0 ) Then
     If ( all_pass ) Then
        Write(*,'(A)') 'verify_ti_profile: all cases PASSED'
     Else
        Write(*,'(A)') 'verify_ti_profile: FAILED'
        Call Mpi_Finalize(ierr)
        Stop 1
     End If
  End If

  Call Mpi_Finalize(ierr)

Contains

  Subroutine reset_profile_state
    If ( Allocated(prof_y)  ) Deallocate( prof_y, prof_U, prof_R11, prof_R22, prof_R33, prof_R12 )
    If ( Allocated(sig_y)   ) Deallocate( sig_y, sig_ux, sig_uy, sig_uz, &
                                           sig_vx, sig_vy, sig_vz, sig_wx, sig_wy, sig_wz )
    n_sigma = 0
  End Subroutine reset_profile_state

  Subroutine expect(cond, label)
    Logical,      Intent(In) :: cond
    Character(*), Intent(In) :: label
    If ( myid /= 0 ) Return
    If ( cond ) Then
       Write(*,'(A,A)') '  PASS: ', label
    Else
       Write(*,'(A,A)') '  FAIL: ', label
       all_pass = .False.
    End If
  End Subroutine expect

  !> Header lists all 14 columns; every sigma slot should be the measured value (tier 1), no fallback
  Subroutine case_full_length_scales

    Integer(Int32) :: unit_out

    If ( myid == 0 ) Write(*,'(A)') '=== full length-scale columns (all measured, tier 1) ==='
    Call reset_profile_state
    inflow_profile_file = 'ti_full.dat'
    sem_profile_format  = 1

    If ( myid == 0 ) Then
       Open(newunit=unit_out, file=inflow_profile_file, status='replace', action='write')
       Write(unit_out,'(A)') '# z U Iu Iv Iw Lux Luy Luz Lvx Lvy Lvz Lwx Lwy Lwz'
       Write(unit_out,'(A)') '0.0  5.0  0.20  0.15  0.10  0.05 0.02 0.01 0.04 0.02 0.01 0.03 0.02 0.01'
       Write(unit_out,'(A)') '1.0 10.0  0.10  0.08  0.06  0.07 0.03 0.02 0.06 0.03 0.02 0.05 0.03 0.02'
       Close(unit_out)
    End If
    Call Mpi_Barrier(MPI_COMM_WORLD, ierr)

    Call read_mean_profile_TI

    Call expect( n_profile == 2, 'n_profile == 2' )
    Call expect( Abs(prof_U(1)-5d0)  < tol, 'prof_U(1) == 5.0' )
    Call expect( Abs(prof_R11(1) - (0.20d0*5d0)**2) < tol, 'prof_R11(1) == (Iu*U)^2' )
    Call expect( Abs(prof_R22(1) - (0.15d0*5d0)**2) < tol, 'prof_R22(1) == (Iv*U)^2' )
    Call expect( Abs(prof_R33(1) - (0.10d0*5d0)**2) < tol, 'prof_R33(1) == (Iw*U)^2' )
    Call expect( prof_R12(1) == 0d0, 'prof_R12 == 0 (no shear data)' )
    Call expect( n_sigma == 2, 'n_sigma == n_profile' )
    Call expect( Abs(sig_uy(1)-0.02d0) < tol, 'sig_uy(1) == measured Luy (tier 1, not ratio-derived)' )
    Call expect( Abs(sig_wz(2)-0.02d0) < tol, 'sig_wz(2) == measured Lwz (tier 1)' )

  End Subroutine case_full_length_scales

  !> Header lists only Lux/Lvx/Lwx; Loy/Loz must fall back to sem_Lscale_ratio_{y,z}*Lox (tier 2)
  Subroutine case_partial_length_scales

    Integer(Int32) :: unit_out

    If ( myid == 0 ) Write(*,'(A)') '=== partial length-scale columns (x-only, tier 2 ratio fallback) ==='
    Call reset_profile_state
    inflow_profile_file = 'ti_partial.dat'
    sem_profile_format  = 1

    If ( myid == 0 ) Then
       Open(newunit=unit_out, file=inflow_profile_file, status='replace', action='write')
       Write(unit_out,'(A)') '# z U Iu Iv Iw Lux'
       Write(unit_out,'(A)') '0.0  5.0  0.20  0.15  0.10  0.05'
       Write(unit_out,'(A)') '1.0 10.0  0.10  0.08  0.06  0.07'
       Close(unit_out)
    End If
    Call Mpi_Barrier(MPI_COMM_WORLD, ierr)

    Call read_mean_profile_TI

    Call expect( n_sigma == 2, 'n_sigma == n_profile' )
    Call expect( Abs(sig_ux(1)-0.05d0) < tol, 'sig_ux(1) == measured Lux (tier 1)' )
    Call expect( Abs(sig_uy(1)-sem_Lscale_ratio_y*0.05d0) < tol, 'sig_uy(1) == ratio_y*Lux (tier 2)' )
    Call expect( Abs(sig_uz(1)-sem_Lscale_ratio_z*0.05d0) < tol, 'sig_uz(1) == ratio_z*Lux (tier 2)' )
    Call expect( Abs(sig_vx(1)-sem_length_scale) < tol, 'sig_vx(1) == sem_length_scale (tier 3, v never measured)' )
    Call expect( Abs(sig_vy(1)-sem_length_scale) < tol, 'sig_vy(1) == sem_length_scale (tier 3)' )

  End Subroutine case_partial_length_scales

  !> No L** columns at all: sig_* must be filled from Prandtl mixing-length theory, L = sigma_u/|dU/dy| (see derive_mixing_length in sem.f90)
  Subroutine case_no_length_scales

    Integer(Int32) :: unit_out
    Real(Int64) :: Lmix_expected

    If ( myid == 0 ) Write(*,'(A)') '=== no length-scale columns (mixing-length fallback) ==='
    Call reset_profile_state
    inflow_profile_file = 'ti_none.dat'
    sem_profile_format  = 1
    sem_use_esem        = 1

    If ( myid == 0 ) Then
       Open(newunit=unit_out, file=inflow_profile_file, status='replace', action='write')
       Write(unit_out,'(A)') '# z U Iu Iv Iw'
       Write(unit_out,'(A)') '0.0  5.0  0.20  0.15  0.10'
       Write(unit_out,'(A)') '1.0 10.0  0.10  0.08  0.06'
       Close(unit_out)
    End If
    Call Mpi_Barrier(MPI_COMM_WORLD, ierr)

    Call read_mean_profile_TI

    ! sigma_u = Iu*U = 1.0 at both points, dU/dy = (10-5)/(1-0) = 5 -> Lmix = 1.0/5 = 0.2, within [sem_length_scale, 0.3*(1-0)]
    Lmix_expected = 0.2d0

    Call expect( n_sigma == 2, 'n_sigma == n_profile (mixing-length profile always filled)' )
    Call expect( Abs(sig_ux(1)-Lmix_expected) < tol, 'sig_ux(1) == sigma_u/|dU/dy| (mixing length)' )
    Call expect( Abs(sig_vx(1)-Lmix_expected) < tol, 'sig_vx(1) == same mixing length (no independent shear estimate for v)' )
    Call expect( Abs(sig_uy(1)-sem_Lscale_ratio_y*Lmix_expected) < tol, 'sig_uy(1) == ratio_y*Lmix' )
    Call expect( Abs(sig_uz(2)-sem_Lscale_ratio_z*Lmix_expected) < tol, 'sig_uz(2) == ratio_z*Lmix' )

  End Subroutine case_no_length_scales

  !> Header/data delimited by tabs, as e.g. Excel's "Save as tab-delimited" produces; regression for a tokenizer bug where only spaces were treated as column separators
  Subroutine case_tab_separated_header

    Integer(Int32) :: unit_out
    Character(1), Parameter :: tab = Achar(9)

    If ( myid == 0 ) Write(*,'(A)') '=== tab-separated header/data ==='
    Call reset_profile_state
    inflow_profile_file = 'ti_tabs.dat'
    sem_profile_format  = 1

    If ( myid == 0 ) Then
       Open(newunit=unit_out, file=inflow_profile_file, status='replace', action='write')
       Write(unit_out,'(A)') '#'//tab//'z'//tab//'U'//tab//'Iu'//tab//'Iv'//tab//'Iw'
       Write(unit_out,'(A)') '0.0'//tab//'5.0'//tab//'0.20'//tab//'0.15'//tab//'0.10'
       Write(unit_out,'(A)') '1.0'//tab//'10.0'//tab//'0.10'//tab//'0.08'//tab//'0.06'
       Close(unit_out)
    End If
    Call Mpi_Barrier(MPI_COMM_WORLD, ierr)

    Call read_mean_profile_TI

    Call expect( n_profile == 2, 'n_profile == 2 (tab-separated header parsed)' )
    Call expect( Abs(prof_U(1)-5d0) < tol, 'prof_U(1) == 5.0 (tab-separated data parsed)' )
    Call expect( Abs(prof_R11(2) - (0.10d0*10d0)**2) < tol, 'prof_R11(2) == (Iu*U)^2' )

  End Subroutine case_tab_separated_header

End Program verify_ti_profile
