! Module for generating grid and initial conditions
Module genGridAndIC

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use mpi
  Use synthetic_eddy_method, Only : mean_profile_U, init_inflow_profile
  !Use ifport   ! Only for intel compiler

  ! prevent implicit typing
  Implicit None

Contains
	!Subroutine to generate the grid!
	Subroutine generateGrid

		! Local variables
		Integer(Int32) :: i, j, n_smooth, io_unit, N_sp, iter_geom
		Real   (Int64) :: eta, dy_k, Ly_smooth, y_start, yp_c, dy_c, &
		                  r_geom, r_lo, r_hi, r_mid, f_lo, f_mid

		If (myid==0) Write(*,'(A,I2)') ' Generating the grid, grid_type =', grid_type

		! ---- x grid: uniform face locations --------------------------
		Do i = 1, nx_global
			x_global(i) = Lx_i * Real(i-1,8) / Real(nx_global-1,8)
		End Do

		! ---- z grid: uniform face locations --------------------------
		Do i = 1, nz_global
			z_global(i) = Lz_i * Real(i-1,8) / Real(nz_global-1,8)
		End Do

		! ---- y grid: V-face locations (tanh stretching; options 5-7 add a
		!  uniform roughness sublayer [0,ks] before the stretched region) ----

		Select Case (grid_type)

		Case (1)
			! Uniform spacing across the full channel height
			Do i = 1, ny_global
				y_global(i) = Ly_i * Real(i-1,8) / Real(ny_global-1,8)
			End Do

		Case (2)
			! Symmetric tanh: fine at both walls, coarse at centre
			Do i = 1, ny_global
				eta = Real(i-1,8) / Real(ny_global-1,8)
				y_global(i) = Ly_i * 0.5d0 * &
					( 1d0 + dtanh(alphaGrid*(eta - 0.5d0)) &
					       / dtanh(0.5d0*alphaGrid) )
			End Do

		Case (3)
			! Single-sided tanh: fine at bottom wall, coarse at top
			Do i = 1, ny_global
				eta = Real(i-1,8) / Real(ny_global-1,8)
				y_global(i) = Ly_i * ( 1d0 - dtanh(alphaGrid*(1d0-eta)) &
				                              / dtanh(alphaGrid) )
			End Do

		Case (4)
			! Single-sided tanh: fine at top wall, coarse at bottom
			Do i = 1, ny_global
				eta = Real(i-1,8) / Real(ny_global-1,8)
				y_global(i) = Ly_i * dtanh(alphaGrid*eta) / dtanh(alphaGrid)
			End Do

		Case (5, 6, 7)
			! ----- Roughness sublayer + stretched smooth region: validate inputs -----
			If ( nks_global < 1 ) Then
				If (myid==0) Write(*,*) 'ERROR: grid_type', grid_type, &
					'requires nks >= 1 in &IBM'
				Call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
			End If

			dy_k      = ks / Real(nks_global, 8)   ! uniform spacing in sublayer
			y_start   = ks + dy_k                   ! first face of smooth region
			Ly_smooth = Ly_i - y_start              ! height of smooth region
			n_smooth  = ny_global - nks_global - 1  ! face-point count in smooth region

			If ( n_smooth < 2 ) Then
				If (myid==0) Write(*,*) 'ERROR: not enough y-points above roughness'
				Call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
			End If

			! Roughness sublayer — uniform
			Do i = 1, nks_global + 1
				y_global(i) = Real(i-1,8) * dy_k
			End Do

			! Smooth region above roughness — stretched
			Select Case (grid_type)

			Case (5)
				! Symmetric tanh: fine at rough interface and at top wall
				Do i = nks_global+2, ny_global
					j   = i - (nks_global+2)
					eta = Real(j,8) / Real(n_smooth-1,8)
					y_global(i) = y_start + Ly_smooth * 0.5d0 * &
						( 1d0 + dtanh(alphaGrid*(eta - 0.5d0)) &
						       / dtanh(0.5d0*alphaGrid) )
				End Do

			Case (6)
				! Geometric stretching, expanding (fine at interface, coarse at top); root-finding details
				N_sp = n_smooth - 1

				If ( Ly_smooth <= 0d0 ) Then
					If (myid==0) Write(*,*) 'ERROR: Ly_smooth <= 0 for grid_type=6'
					Call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
				End If

				If ( N_sp <= 1 .Or. &
				     Abs(Real(N_sp,8)*dy_k - Ly_smooth) < 1d-10*Ly_smooth ) Then
					r_geom = 1.0d0
				Else
					If ( Real(N_sp,8)*dy_k >= Ly_smooth ) Then
						If (myid==0) Then
							Write(*,*) 'ERROR: grid_type=6 requires Ly_smooth > N_sp*dy_k'
							Write(*,*) '  N_sp*dy_k =', Real(N_sp,8)*dy_k, &
							           '  Ly_smooth =', Ly_smooth
						End If
						Call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
					End If
					r_lo = 1.0d0 + 1d-10
					r_hi = Max(2.0d0, &
					       2.0d0*(Ly_smooth/dy_k)**(1.0d0/Real(Max(N_sp-1,1),8)))
					Do iter_geom = 1, 60   ! grow r_hi until f(r_hi) >= 0
						If ( dy_k*(r_hi**N_sp - 1.0d0)/(r_hi - 1.0d0) &
						     >= Ly_smooth ) Exit
						r_hi = 2.0d0 * r_hi
					End Do
					Do iter_geom = 1, 200
						r_mid = 0.5d0*(r_lo + r_hi)
						f_mid = dy_k*(r_mid**N_sp - 1.0d0)/(r_mid - 1.0d0) - Ly_smooth
						If ( Abs(f_mid) < 1d-12*Ly_smooth ) Exit
						f_lo  = dy_k*(r_lo **N_sp - 1.0d0)/(r_lo  - 1.0d0) - Ly_smooth
						If ( f_lo*f_mid <= 0.0d0 ) Then
							r_hi = r_mid
						Else
							r_lo = r_mid
						End If
					End Do
					r_geom = r_mid
				End If

				Do i = nks_global+2, ny_global
					j = i - (nks_global+2)
					If ( Abs(r_geom - 1.0d0) < 1d-10 ) Then
						y_global(i) = y_start + Real(j,8)*dy_k
					Else
						y_global(i) = y_start + &
							dy_k*(r_geom**j - 1.0d0)/(r_geom - 1.0d0)
					End If
				End Do

				If (myid==0) Write(*,'(A,ES12.5,A,ES12.5,A,ES12.5)') &
					'   [grid_type=6] geom ratio r=', r_geom, &
					'  dy_interface=', dy_k, &
					'  dy_top=', dy_k*r_geom**Max(N_sp-1,1)

			Case (7)
				! Geometric stretching, contracting (coarse at interface, fine at top); root-finding details
				N_sp = n_smooth - 1

				If ( Ly_smooth <= 0d0 ) Then
					If (myid==0) Write(*,*) 'ERROR: Ly_smooth <= 0 for grid_type=7'
					Call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
				End If

				If ( N_sp <= 1 .Or. &
				     Abs(Real(N_sp,8)*dy_k - Ly_smooth) < 1d-10*Ly_smooth ) Then
					r_geom = 1.0d0
				Else
					If ( Real(N_sp,8)*dy_k <= Ly_smooth ) Then
						If (myid==0) Then
							Write(*,*) 'ERROR: grid_type=7 requires Ly_smooth < N_sp*dy_k'
							Write(*,*) '  N_sp*dy_k =', Real(N_sp,8)*dy_k, &
							           '  Ly_smooth =', Ly_smooth
						End If
						Call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
					End If
					r_lo = 1d-8
					r_hi = 1.0d0 - 1d-10
					Do iter_geom = 1, 200
						r_mid = 0.5d0*(r_lo + r_hi)
						f_mid = dy_k*(r_mid**N_sp - 1.0d0)/(r_mid - 1.0d0) - Ly_smooth
						If ( Abs(f_mid) < 1d-12*Ly_smooth ) Exit
						f_lo  = dy_k*(r_lo **N_sp - 1.0d0)/(r_lo  - 1.0d0) - Ly_smooth
						If ( f_lo*f_mid <= 0.0d0 ) Then
							r_hi = r_mid
						Else
							r_lo = r_mid
						End If
					End Do
					r_geom = r_mid
				End If

				Do i = nks_global+2, ny_global
					j = i - (nks_global+2)
					If ( Abs(r_geom - 1.0d0) < 1d-10 ) Then
						y_global(i) = y_start + Real(j,8)*dy_k
					Else
						y_global(i) = y_start + &
							dy_k*(r_geom**j - 1.0d0)/(r_geom - 1.0d0)
					End If
				End Do

				If (myid==0) Write(*,'(A,ES12.5,A,ES12.5,A,ES12.5)') &
					'   [grid_type=7] geom ratio r=', r_geom, &
					'  dy_interface=', dy_k, &
					'  dy_top=', dy_k*r_geom**Max(N_sp-1,1)

			End Select

		Case Default
			If (myid==0) Write(*,*) 'ERROR: unknown grid_type =', grid_type, &
				'(valid range: 1-7)'
			Call MPI_Abort(MPI_COMM_WORLD, 1, ierr)

		End Select

		! --- Write geometry.out (cell counts + domain lengths) and
		!     grid.out (y-grid table) for GenSDF ---
		If (myid == 0) Then
			Call system('mkdir -p fields/')
			Open(newunit=io_unit, file='fields/geometry.out', status='replace')
			Write(io_unit, '(3(I0,1X))') nxm_global, nym_global, nzm_global
			Write(io_unit, '(3(ES24.16,1X))') Lx_i, Ly_i, Lz_i
			Close(io_unit)
			Open(newunit=io_unit, file='fields/grid.out', status='replace')
			Do i = 1, nym_global
				yp_c = 0.5d0 * (y_global(i) + y_global(i+1))
				dy_c = y_global(i+1) - y_global(i)
				Write(io_unit, '(I6,4(1X,ES24.16))') i, yp_c, y_global(i+1), dy_c, 1.0d0/dy_c
			End Do
			Close(io_unit)
			Write(*,'(A)') '   Wrote fields/geometry.out and fields/grid.out for GenSDF'
		End If

	End Subroutine generateGrid


	! Generate the initial condition
	Subroutine generateIC

		! Local variables
		Integer(Int32) :: ii, jj, kk
		Real   (Int64) :: y_c, y_wall, y_plus, U_base
		Real   (Int64) :: u_tau, noise_frac
		Real   (Int64) :: rnd   ! white-noise scratch (Random_number, not GNU rand(), for compiler portability)

		! Log-law constants (must match wallmodel.f90)
		Real(Int64), Parameter :: kappa_ic = 0.41d0
		Real(Int64), Parameter :: B_ic     = 5.2d0

		! channel_perturb IC (ic_type==4)
		Integer(Int32) :: kk_glob, ii_glob
		Real   (Int64) :: max_u_g, max_v_g, max_w_g, sum_u_g
		Integer(Int32) :: n_int_ch, ii_int_ch
		Real   (Int64) :: Re_tau_ch, Re_b_ch, C_ch, k_ch
		Real   (Int64) :: h_ch, u_tau_ch
		Real   (Int64) :: dyp_ch, yp_int, up_int, sum_ubulk
		Real   (Int64) :: y_norm, yp_ch
		Real   (Int64) :: llx_ch, llz_ch
		Real   (Int64) :: alpha_ls, beta_ls, alpha_ss, beta_ss
		Real   (Int64) :: x_c, z_c, z_f, ran_ch

		! Taylor-Green Vortex IC (ic_type==6)
		Real   (Int64) :: kx_tgv, ky_tgv, kz_tgv, x_f

		noise_frac = noise_percent * 0.01d0

		! ---- Taylor-Green Vortex IC (ic_type==6): fully analytic, deterministic,
		! exactly divergence-free on any box aspect ratio -- no noise added.
		! Requires x_bc_type==0 and y_bc_type==0 (validated in initialization.f90).
		If ( ic_type == 6 ) Then

			kx_tgv = 2d0*pi / Lx_i
			ky_tgv = 2d0*pi / Ly_i
			kz_tgv = 2d0*pi / Lz_i

			If ( myid == 0 ) Write(*,'(A)') '   IC = Taylor-Green Vortex (ic_type=6), Utarget used as amplitude U0'

			! -- U: x-face, y-center, z-center --
			U = 0d0
			Do kk = 2, nzg-1
				kk_glob = k1_global(myid) + kk - 1
				z_c = 0.5d0 * ( z_global(kk_glob-1) + z_global(kk_glob) )
				Do jj = 2, nyg_global-1
					y_c = 0.5d0 * ( y_global(jj-1) + y_global(jj) )
					Do ii = 1, nx
						ii_glob = i1_global(myid) + ii - 1
						x_f = x_global(ii_glob)
						U(ii,jj,kk) = Utarget * dsin(kx_tgv*x_f) * dcos(ky_tgv*y_c) * dcos(kz_tgv*z_c)
					End Do
				End Do
			End Do

			! -- V: x-center, y-face, z-center -- faces 1 and ny_global are the
			! same physical periodic point (y=0==y=Ly), so the plain analytic
			! formula is exact there with no wraparound special-casing needed
			V = 0d0
			Do kk = 2, nzg-1
				kk_glob = k1_global(myid) + kk - 1
				z_c = 0.5d0 * ( z_global(kk_glob-1) + z_global(kk_glob) )
				Do jj = 2, ny_global-1
					Do ii = 2, nxg-1
						ii_glob = ig1_global(myid) + ii - 2
						x_c = 0.5d0 * ( x_global(ii_glob-1) + x_global(ii_glob) )
						V(ii,jj,kk) = -Utarget*(kx_tgv/ky_tgv) * dcos(kx_tgv*x_c) * dsin(ky_tgv*y_global(jj)) * dcos(kz_tgv*z_c)
					End Do
				End Do
			End Do

			! -- W = 0 everywhere (TGV has no spanwise velocity) --
			W = 0d0

			Call Mpi_allreduce( MaxVal(Abs(U)), max_u_g, 1, MPI_real8, MPI_MAX, MPI_COMM_WORLD, ierr )
			Call Mpi_allreduce( MaxVal(Abs(V)), max_v_g, 1, MPI_real8, MPI_MAX, MPI_COMM_WORLD, ierr )
			If ( myid==0 ) Then
				Write(*,'(A,E12.4)') '   IC Max |U| = ', max_u_g
				Write(*,'(A,E12.4)') '   IC Max |V| = ', max_v_g
				Write(*,'(A,E12.4)') '   IC Max |W| = ', 0d0
			End If

			Return  ! skip log-law/Reichardt setup, Select Case, and noise loops entirely

		End If

		! x_bc_type==1: load the inflow mean profile now so it can seed U_base below, avoiding a step-1 IC/BC divergence spike
		If ( x_bc_type == 1 ) Call init_inflow_profile

		! ---- Reichardt profile parameters for channel_perturb IC --------
		C_ch   = 5.17d0
		k_ch   = 0.41d0
		llx_ch = Lx_i
		llz_ch = Lz_i
		alpha_ls = 3.0d0  * 2.0d0*pi / llx_ch
		beta_ls  = 4.0d0  * 2.0d0*pi / llz_ch
		alpha_ss = 17.0d0 * 2.0d0*pi / llx_ch
		beta_ss  = 13.0d0 * 2.0d0*pi / llz_ch

		! Re_tau derived from dPdx, nu, and channel reference height
		! (same physics as the u_tau estimate for ic_type==1)
		If ( bc_face_yhi == 2 ) Then
			h_ch = Ly_i           ! open channel: full height
		Else
			h_ch = 0.5d0 * Ly_i   ! closed channel: half-height
		End If
		If ( Abs(dPdx) > 0d0 .And. nu > 0d0 ) Then
			u_tau_ch  = Sqrt( Abs(dPdx) * h_ch )
			Re_tau_ch = u_tau_ch * h_ch / nu
		Else
			! Fallback when dPdx is not set: rough estimate from Utarget
			u_tau_ch  = 0.05d0 * Max( Abs(Utarget), 1d-10 )
			Re_tau_ch = u_tau_ch * h_ch / Max( nu, 1d-15 )
		End If

		! Re_b from trapezoidal integration of the Reichardt profile
		! (exact formula; consistent across all Re_tau).
		n_int_ch  = 10000
		dyp_ch    = Re_tau_ch / Real(n_int_ch, Int64)
		sum_ubulk = 0d0
		Do ii_int_ch = 0, n_int_ch
			yp_int = Max( Real(ii_int_ch, Int64) * dyp_ch, 1d-10 )  ! avoid log(0)
			up_int = (1d0/k_ch) * log(1d0 + k_ch*yp_int) + &
			         (C_ch - (1d0/k_ch)*log(k_ch)) * &
			         (1d0 - exp(-yp_int/11d0) - yp_int/11d0*exp(-yp_int/3d0))
			If ( ii_int_ch == 0 .Or. ii_int_ch == n_int_ch ) Then
				sum_ubulk = sum_ubulk + 0.5d0 * up_int  ! trapezoidal endpoints
			Else
				sum_ubulk = sum_ubulk + up_int
			End If
		End Do
		Re_b_ch = sum_ubulk * dyp_ch
		If ( myid==0 ) Write(*,'(A,F9.1,A,F9.1,A,E12.4)') &
			'   IC channel_perturb: Re_tau =', Re_tau_ch, &
			'  Re_b =', Re_b_ch, '  u_tau =', u_tau_ch

		! ---- Estimate friction velocity for log-law IC ---------------
		u_tau = 0d0
		If ( ic_type == 1 ) Then
			If ( Abs(dPdx) > 0d0 ) Then
				If ( bc_face_yhi == 2 ) Then
					! Open channel / free-surface: tau_w = |dPdx| * Ly
					u_tau = Sqrt( Abs(dPdx) * Ly_i )
				Else
					! Full channel (both walls no-slip): tau_w = |dPdx| * Ly/2
					u_tau = Sqrt( Abs(dPdx) * 0.5d0 * Ly_i )
				End If
			Else
				! No imposed pressure gradient — rough estimate from Utarget
				u_tau = 0.05d0 * Abs(Utarget)
			End If
			If ( myid==0 ) Write(*,'(A,E12.4)') '   IC log-law u_tau = ', u_tau
		End If

		! ---- Build mean U profile; y_c approximated as midpoint of V-face points ----

		U = 0d0
		Do jj = 2, nyg_global-1

			y_c = 0.5d0 * ( y_global(jj-1) + y_global(jj) )

			Select Case (ic_type)

			Case (1)  ! Log-law profile
				If ( bc_face_yhi == 1 ) Then
					! Both walls no-slip: mirror about channel centreline
					y_wall = Min( y_c, Ly_i - y_c )
				Else
					! Free-slip top: log-law from bottom wall only
					y_wall = y_c
				End If
				y_wall = Max( y_wall, 0d0 )
				If ( y_wall > 0d0 .And. u_tau > 0d0 .And. nu > 0d0 ) Then
					y_plus = Max( y_wall * u_tau / nu, 1d0 )
					U_base = u_tau * ( Log(y_plus) / kappa_ic + B_ic )
				Else
					U_base = 0d0
				End If

			Case (2)  ! Linear profile
				If ( bc_face_yhi == 1 ) Then
					! Both walls no-slip: tent function peaking at Ly/2
					U_base = Utarget * 2d0 * Min( y_c, Ly_i - y_c ) / Ly_i
				Else
					! Free-slip top: ramp from 0 at bottom to Utarget at top
					U_base = Utarget * y_c / Ly_i
				End If

			Case (3)  ! Zero mean — noise only
				U_base = 0d0

			Case (4)  ! Reichardt turbulent mean profile (channel_perturb)
				y_norm = 2.0d0 * y_c / Ly_i - 1.0d0
				If (y_norm .ge. 0.0d0) Then
					yp_ch = (1.0d0 - y_norm) * Re_tau_ch
				Else
					yp_ch = (1.0d0 + y_norm) * Re_tau_ch
				End If
				U_base = 1.0d0/k_ch * log(1.0d0 + k_ch*yp_ch) + &
				         (C_ch - (1.0d0/k_ch)*log(k_ch)) * &
				         (1.0d0 - exp(-yp_ch/11.0d0) - yp_ch/11.0d0*exp(-yp_ch/3.0d0))
			! Scale by u_tau (profile is U+; centerline ≈ U_cl+ * u_tau).
			U_base = U_base * u_tau_ch

			Case (5)  ! Inverse-linear profile: large near walls, zero at centre
				! Closed channel: anti-tent (U=Utarget at walls, 0 at centre); open channel: inverted ramp (U=Utarget at bottom, 0 at top); strong shear promotes transition
				If ( bc_face_yhi == 1 ) Then
					U_base = Utarget * ( 1d0 - 2d0*Min(y_c, Ly_i - y_c)/Ly_i )
				Else
					U_base = Utarget * ( 1d0 - y_c/Ly_i )
				End If

			End Select

			! x_bc_type==1: seed the interior mean with the inflow BC's own profile instead of ic_type's, so step 1 doesn't see a large inflow-face divergence
			If ( x_bc_type == 1 ) U_base = mean_profile_U( y_c )

			U(:, jj, :) = U_base

		End Do

		! ---- Ghost-cell BCs for U ------------------------------------
		! (y_bc_type==0: left unfilled here -- apply_periodic_bc_y fills them
		! on the first RK-stage call, before any RHS reads a ghost cell)
		If ( y_bc_type == 1 ) Then
			U(:, 1, :) = -U(:, 2, :)                        ! bottom: always no-slip
			If ( bc_face_yhi == 1 ) Then
				U(:, nyg_global, :) = -U(:, nyg_global-1, :)  ! top: no-slip
			Else
				U(:, nyg_global, :) =  U(:, nyg_global-1, :)  ! top: free-slip
			End If
		End If

		! ---- channel_perturb IC: 3-D structured + random perturbations ----
		! Replaces the white-noise loops entirely for ic_type == 4.
		If ( ic_type == 4 ) Then

			! -- U: add structured perturbations (x-face, y-center, z-center) --
			Do kk = 2, nzg-1
				kk_glob = k1_global(myid) + kk - 1
				z_c = 0.5d0 * ( z_global(kk_glob-1) + z_global(kk_glob) )
				Do jj = 2, nyg_global-1
					Do ii = 1, nx
						ii_glob = i1_global(myid) + ii - 1
						x_c = x_global(ii_glob)
						U(ii,jj,kk) = U(ii,jj,kk) &
							+ noise_frac      *beta_ls  * sin(alpha_ls*x_c)*cos(beta_ls*z_c) &
							+ 0.1d0*noise_frac*beta_ss * sin(alpha_ss*x_c)*cos(beta_ss*z_c)
					End Do
				End Do
			End Do

			! Re-enforce ghost-cell BCs after perturbations
			If ( y_bc_type == 1 ) Then
				U(:, 1, :) = -U(:, 2, :)
				If ( bc_face_yhi == 1 ) Then
					U(:, nyg_global, :) = -U(:, nyg_global-1, :)
				Else
					U(:, nyg_global, :) =  U(:, nyg_global-1, :)
				End If
			End If

			! -- V: structured + random perturbations (x-center, y-face, z-center) --
			V = 0.0d0
			Do kk = 2, nzg-1
				kk_glob = k1_global(myid) + kk - 1
				z_c = 0.5d0 * ( z_global(kk_glob-1) + z_global(kk_glob) )
				Do jj = 2, ny_global-1
					y_c    = y_global(jj)
					y_norm = 2.0d0 * y_c / Ly_i - 1.0d0
					Do ii = 2, nxg-1
						ii_glob = ig1_global(myid) + ii - 2
						x_c = 0.5d0 * ( x_global(ii_glob-1) + x_global(ii_glob) )
						! Structured perturbations (large- and small-scale)
						V(ii,jj,kk) = noise_frac       * sin(alpha_ls*x_c)*sin(beta_ls*z_c) &
						            + 0.1d0*noise_frac  * sin(alpha_ss*x_c)*sin(beta_ss*z_c)
						! Random perturbation
						ran_ch = sin( -20.0d0*x_c*z_c + y_norm**3*tan(x_c*z_c**2) + &
						              100.0d0*z_c*y_norm - 20.0d0*sin(x_c*y_norm*z_c)**5 )
						V(ii,jj,kk) = V(ii,jj,kk) + 0.02d0*noise_frac * ran_ch
					End Do
				End Do
			End Do

			! -- W: structured perturbations (x-center, y-center, z-face) --
			W = 0.0d0
			Do kk = 1, nz
				z_f = z_global( k1_global(myid) + kk - 1 )
				Do jj = 2, nyg_global-1
					Do ii = 2, nxg-1
						ii_glob = ig1_global(myid) + ii - 2
						x_c = 0.5d0 * ( x_global(ii_glob-1) + x_global(ii_glob) )
						W(ii,jj,kk) = -noise_frac      *alpha_ls * cos(alpha_ls*x_c)*sin(beta_ls*z_f) &
						              -0.1d0*noise_frac*alpha_ss * cos(alpha_ss*x_c)*sin(beta_ss*z_f)
					End Do
				End Do
			End Do

			! global-sum diagnostics: MPI_Allreduce over each rank's local extent
			Call Mpi_allreduce( MaxVal(Abs(U)), max_u_g, 1, MPI_real8, MPI_MAX, MPI_COMM_WORLD, ierr )
			Call Mpi_allreduce( MaxVal(Abs(V)), max_v_g, 1, MPI_real8, MPI_MAX, MPI_COMM_WORLD, ierr )
			Call Mpi_allreduce( MaxVal(Abs(W)), max_w_g, 1, MPI_real8, MPI_MAX, MPI_COMM_WORLD, ierr )
			Call Mpi_allreduce( Sum(U(1:nx, 2:nyg_global-1, 2:nzg-1)), sum_u_g, 1, MPI_real8, MPI_SUM, MPI_COMM_WORLD, ierr )
			If ( myid==0 ) Then
				Write(*,'(A,E12.4)') '   IC Max |U| = ', max_u_g
				Write(*,'(A,E12.4)') '   IC Max |V| = ', max_v_g
				Write(*,'(A,E12.4)') '   IC Max |W| = ', max_w_g
				Write(*,'(A,E12.4)') '   IC Mean U  = ', &
					sum_u_g / Real(nx_global*(nyg_global-2)*(nzm_global), Int64)
			End If

			Return  ! skip white-noise loops

		End If

		! ---- Add white noise to interior U planes --------------------
		Do kk = 1, nzg
			Do jj = 2, nyg_global-1
				Do ii = 1, nx
					Call Random_number(rnd)
					U(ii,jj,kk) = U(ii,jj,kk) + &
						noise_frac * Utarget * 2d0*(rnd-0.5d0)
				End Do
			End Do
		End Do

		! Re-impose ghost-cell BCs after noise (keeps BCs clean)
		If ( y_bc_type == 1 ) Then
			U(:, 1, :) = -U(:, 2, :)
			If ( bc_face_yhi == 1 ) Then
				U(:, nyg_global, :) = -U(:, nyg_global-1, :)
			Else
				U(:, nyg_global, :) =  U(:, nyg_global-1, :)
			End If
		End If

		! ---- V: zero mean + noise (y-face velocity component) --------
		V = 0d0
		Do kk = 1, nzg
			Do jj = 1, ny_global
				Do ii = 1, nxg
					Call Random_number(rnd)
					V(ii,jj,kk) = noise_frac * Utarget * 2d0*(rnd-0.5d0)
				End Do
			End Do
		End Do

		! ---- W: zero mean + noise (z-face velocity component) --------
		W = 0d0
		Do kk = 1, nz
			Do jj = 1, nyg_global
				Do ii = 1, nxg
					Call Random_number(rnd)
					W(ii,jj,kk) = noise_frac * Utarget * 2d0*(rnd-0.5d0)
				End Do
			End Do
		End Do

		! ---- Diagnostics: MPI_Allreduce over each rank's local extent, print on rank 0 ----
		Call Mpi_allreduce( MaxVal(Abs(U)), max_u_g, 1, MPI_real8, MPI_MAX, MPI_COMM_WORLD, ierr )
		Call Mpi_allreduce( MaxVal(Abs(V)), max_v_g, 1, MPI_real8, MPI_MAX, MPI_COMM_WORLD, ierr )
		Call Mpi_allreduce( MaxVal(Abs(W)), max_w_g, 1, MPI_real8, MPI_MAX, MPI_COMM_WORLD, ierr )
		Call Mpi_allreduce( Sum(U(1:nx, 2:nyg_global-1, 2:nzg-1)), sum_u_g, 1, MPI_real8, MPI_SUM, MPI_COMM_WORLD, ierr )
		If ( myid==0 ) Then
			Write(*,'(A,E12.4)') '   IC Max |U| = ', max_u_g
			Write(*,'(A,E12.4)') '   IC Max |V| = ', max_v_g
			Write(*,'(A,E12.4)') '   IC Max |W| = ', max_w_g
			Write(*,'(A,E12.4)') '   IC Mean U  = ', &
				sum_u_g / Real(nx_global*(nyg_global-2)*(nzm_global), Int64)
		End If

	End Subroutine generateIC

End Module genGridAndIC
