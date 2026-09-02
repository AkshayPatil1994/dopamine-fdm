! Module with initialization of global variables
Module initialization

  ! Modules
  Use, Intrinsic :: iso_c_binding
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use mpi
  Use decomp, Only : decomp_init_pencil, decomp_build_xz_ranges, decomp_init_poisson_pencil, decomp_poisson, z_periodic_partner
  Use input_output
  Use ibmSetup
  Use scalar_transport, Only : compute_settling_velocity
  Use synthetic_eddy_method, Only : init_inflow, init_ti_rescale

  ! prevent implicit typing
  Implicit None

  ! declarations
Contains

  ! Initialize everything
  Subroutine initialize(input_file)

    Character(*), Intent(In), Optional :: input_file
    Integer(Int32) :: i, j, k, kk, nzpe, pos, ipos
    Real   (Int64) :: dy1, dy2, det, a, b, c, r, Qflow_ref
    Real   (Int64) :: dy   ! uniform y spacing, periodic-y path only (y_bc_type==0)
    Integer(Int32), Dimension(:,:), Allocatable :: A_kmodes, A_kmodes_local
    Logical        :: is_first_z, is_last_z
    Integer(Int32) :: z_partner

    ! first initialize MPI
    call Mpi_init(ierr)
    call Mpi_comm_size(MPI_COMM_WORLD, nprocs, ierr)
    call Mpi_comm_rank(MPI_COMM_WORLD,   myid, ierr)

    time_wall_start = MPI_WTIME()

    If (myid==0) Then
       Write(*,'(A)') ' '
       Write(*,'(A)') ' +--------------------------------------------------------------------+'
       Write(*,'(A)') ' |      ___                                                           |'
       Write(*,'(A)') ' |     (  _`\                                  _                      |'
       Write(*,'(A)') ' |     | | ) |   _    _ _      _ _   ___ ___  (_)  ___     __         |'
       Write(*,'(A)') " |     | | | ) /'_`\ ( '_`\  /'_` )/' _ ` _ `\| |/' _ `\ /'__`\       |"
       Write(*,'(A)') ' |     | |_) |( (_) )| (_) )( (_| || ( ) ( ) || || ( ) |(  ___/       |'
       Write(*,'(A)') " |     (____/'`\___/'| ,__/'`\__,_)(_) (_) (_)(_)(_) (_)`\____)       |"
       Write(*,'(A)') ' |                   | |                                              |'
       Write(*,'(A)') ' |                   (_)   DNS / LES / IBM  Navier-Stokes Solver      |'
       Write(*,'(A)') ' +--------------------------------------------------------------------+'
       Write(*,'(A)') ' '
       Write(*,'(A,I0,A)') ' Running on ', nprocs, ' MPI ranks'
       Write(*,'(A)') ' '
    End If

    ! read parameters from standard input
    If ( myid==0 ) Then
       If ( Present(input_file) ) Then
          Write(*,'(A)') ' [1/4] Reading input parameters from '//Trim(input_file)//' ...'
       Else
          Write(*,'(A)') ' [1/4] Reading input parameters ...'
       End If
    End If
    If ( Present(input_file) ) Then
       Call read_input_parameters(input_file)
    Else
       Call read_input_parameters
    End If

#ifdef GPU_POISSON
    ! single-GPU cuFFT Poisson solve (periodic or DCT-IV) only supports nprocs==1
    If ( nprocs /= 1 ) Stop 'ERROR: GPU_POISSON build only supports nprocs=1 (single-GPU); rebuild without ENABLE_GPU for multi-rank runs'
    If ( x_bc_type /= 0 .And. x_bc_type /= 1 ) Stop 'ERROR: GPU_POISSON build only supports x_bc_type=0 or 1'
    If ( y_bc_type == 0 .And. x_bc_type /= 0 ) Stop 'ERROR: GPU_POISSON with y_bc_type=0 (periodic y) requires x_bc_type=0 too'
#endif

    ! time: on restart advance physical time to match nstep_init * dt
    If ( restart == 1 ) Then
       t = Real(nstep_init, Int64) * dt
    Else
       t = 0d0
    End If

    ! time-based save cadence (only meaningful when nsave < 0): first
    ! snapshot is due one tsave interval after the run's starting time
    If ( nsave < 0 ) tsave_next = t + tsave
    
    ! restrictions
    If ( Mod( nx_global, 2 )/=0 ) Stop 'Error: nx must be even'
    If ( Mod( nz_global, 2 )/=0 ) Stop 'Error: nz must be even'
    If ( y_bc_type == 0 .And. grid_type /= 1 ) Stop 'ERROR: y_bc_type=0 (periodic y) requires grid_type=1 (uniform y grid)'
    If ( ic_type == 6 .And. ( x_bc_type /= 0 .Or. y_bc_type /= 0 ) ) &
       Stop 'ERROR: ic_type=6 (Taylor-Green Vortex) requires x_bc_type=0 and y_bc_type=0 (fully periodic box)'
    If ( flat_wall_model_flag == 2 ) Then
       If ( bc_face_ylo /= 2 .And. z0_ylo <= 0d0 ) &
          Stop 'ERROR: flat_wall_model_flag=2 (rough z0 EQWM) requires z0_ylo > 0 for a no-slip bottom wall'
       If ( bc_face_yhi /= 2 .And. z0_yhi <= 0d0 ) &
          Stop 'ERROR: flat_wall_model_flag=2 (rough z0 EQWM) requires z0_yhi > 0 for a no-slip top wall'
    End If
    If ( T_bc_bot == 2 .Or. T_bc_top == 2 ) Then
       If ( boussinesq_flag < 1 ) &
          Stop 'ERROR: T_bc_bot/top=2 (rough EQWM flux BC) requires boussinesq_flag >= 1'
       If ( flat_wall_model_flag /= 2 ) &
          Stop 'ERROR: T_bc_bot/top=2 (rough EQWM flux BC) requires flat_wall_model_flag=2'
       If ( T_bc_bot == 2 .And. z0h_ylo <= 0d0 ) &
          Stop 'ERROR: T_bc_bot=2 requires z0h_ylo > 0'
       If ( T_bc_top == 2 .And. z0h_yhi <= 0d0 ) &
          Stop 'ERROR: T_bc_top=2 requires z0h_yhi > 0'
    End If

    ! domain decomposition: resolves p_row/p_col (auto-factorized when both are 0, see decomp_auto_factorize), then builds this rank's x/z ownership ranges via 2decomp&fft's own distribute() algorithm, replicated locally -- no restriction that nx_global/nz_global divide evenly by p_row/p_col
    Call decomp_init_pencil
    If ( myid==0 ) Write(*,'(A,I0,A,I0)') '   p_row, p_col (resolved)     =  ', p_row, '  ', p_col
    Call decomp_build_xz_ranges

    ! restriction for MPI boundaries (ghost-cell stencils need >=2 interior cells per rank)
    If ( Any( (kg2_global-kg1_global-1) < 2 ) ) Stop 'Error: each rank needs at least 2 interior z-cells'
    If ( Any( (ig2_global-ig1_global-1) < 2 ) ) Stop 'Error: each rank needs at least 2 interior x-cells'

    ! face points
    nx = i2_global(myid) - i1_global(myid) + 1
    ny = ny_global
    nz = k2_global(myid) - k1_global(myid) + 1

    ! middle points
    nxm_global = nx_global - 1
    nym_global = ny_global - 1
    nzm_global = nz_global - 1

    nxm = ig2_global(myid) - ig1_global(myid) + 1 - 2
    nym = ny - 1
    nzm = kg2_global(myid) - kg1_global(myid) + 1 - 2

    ! middle points + ghost cells
    nxg_global = nxm_global + 2
    nyg_global = nym_global + 2
    nzg_global = nzm_global + 2

    nxg = nxm + 2
    nyg = nym + 2
    nzg = kg2_global(myid) - kg1_global(myid) + 1

    ! Allocate main arrays
    If ( myid==0 ) Write(*,'(A)') ' [2/4] Allocating arrays ...'
    Allocate ( x_global (  nx_global),  y_global (  ny_global),  z_global (  nz_global)  )
    Allocate ( xm_global( nxm_global),  ym_global( nym_global),  zm_global( nzm_global)  )
    Allocate ( xg_global(nxm_global+2), yg_global(nym_global+2), zg_global(nzm_global+2) )

    Allocate (  x (  nx),  y (  ny),  z (  nz) )
    Allocate (  xm( nxm),  ym( nym),  zm( nzm) )
    Allocate ( xg(nxm+2), yg(nym+2), zg(nzm+2) )

    Allocate ( yg_m (nyg-1) )
    Allocate ( yg_mm(nyg-2) )
    


    ! global interior + boundary + ghost points
    Allocate (U (    nx, nym+2, nzm+2) )
    Allocate (V ( nxm+2,    ny, nzm+2) )
    Allocate (W ( nxm+2, nym+2,    nz) )
    Allocate (P ( nxm+2, nym+2, nzm+2) )
    P = 0d0

    Allocate (Uo  (    nx,  nym+2, nzm+2) )
    Allocate (Vo  ( nxm+2,     ny, nzm+2) )
    Allocate (Wo  ( nxm+2,  nym+2,    nz) )
    Allocate (Po  ( nxm+2,  nym+2, nzm+2) )

    ! Auxiliary arrays
    Allocate ( term   ( nxg, nyg, nzm+2 ) ) 
    Allocate ( term_1 ( nxg, nyg, nzm+2 ) ) 
    Allocate ( term_2 ( nxg, nyg, nzm+2 ) ) 

    ! RHS: interior points only
    Allocate ( rhs_uo ( 2:nx-1,  2:nyg-1, 2:nzg-1 ) ) 
    Allocate ( rhs_vo ( 2:nxg-1, 2:ny-1,  2:nzg-1 ) )
    Allocate ( rhs_wo ( 2:nxg-1, 2:nyg-1, 2:nz-1  ) )
    Allocate ( rhs_p  ( 2:nxg,   2:nyg-1, 2:nzg   ) ) ! ONE EXTRA PLANE IN X AND Z FOR GHOST CELLS
    rhs_p = 0d0
    
    ! read data
    If ( myid==0 ) Then
       If ( restart == 1 ) Then
          Write(*,'(A,A,A)') ' [3/4] Restarting from ', Trim(Adjustl(filein)), ' ...'
       Else
          Write(*,'(A)') ' [3/4] Generating initial condition ...'
       End If
    End If
    If ( restart == 1 ) Then
       Call read_input_data
    Else
       Call init_flow
    End If

    ! definie global grids from x_global, y_global and z_global (face to centers)
    ! local faces
    x = x_global( i1_global(myid):i2_global(myid) )
    y = y_global
    z = z_global( k1_global(myid):k2_global(myid) )

    ! global interior centers
    Do i = 1, nxm_global
       xm_global(i) = 0.5d0*( x_global(i) + x_global(i+1) )
    End Do
    Do j=1,nym_global
       ym_global(j) = 0.5d0*( y_global(j) + y_global(j+1) )
    End Do
    Do k=1,nzm_global
       zm_global(k) = 0.5d0*( z_global(k) + z_global(k+1) )
    End Do

    ! local interior centers
    xm = xm_global( ig1_global(myid):ig2_global(myid)-2 )
    ym = ym_global
    zm = zm_global( kg1_global(myid):kg2_global(myid)-2 )

    ! global 
    xg_global(2:nxm_global+1) = xm_global
    xg_global(1)              = xm_global(1)          - 2d0*(xm_global(1)-x_global(1))
    xg_global(nxm_global+2)   = xm_global(nxm_global) + 2d0*(x_global(nx_global)-xm_global(nxm_global))
    
    yg_global(2:nym_global+1) = ym_global
    yg_global(1)              = ym_global(1)          - 2d0*(ym_global(1)-y_global(1))
    yg_global(nym_global+2)   = ym_global(nym_global) + 2d0*(y_global(ny_global)-ym_global(nym_global))
    
    zg_global(2:nzm_global+1) = zm_global
    zg_global(1)              = zm_global(1)          - 2d0*(zm_global(1)-z_global(1))
    zg_global(nzm_global+2)   = zm_global(nzm_global) + 2d0*(z_global(nz_global)-zm_global(nzm_global))    

    xg = xg_global( ig1_global(myid):ig2_global(myid) )
    yg = yg_global
    zg = zg_global( kg1_global(myid):kg2_global(myid) )

    ! middle points for yg (.not. equal to y in general)
    yg_m = 0.5d0*( yg(2:nyg) + yg(1:nyg-1) )

    ! middle points for yg_m (.not. equal to ym in general)
    yg_mm = 0.5d0*( yg_m(2:nyg-1) + yg_m(1:nyg-2) )

    ! local minimum grid size for CFL
    dxmin = Minval ( xg_global(2:nxg_global) - xg_global(1:nxg_global-1) )
    dymin = Minval ( yg_global(2:nyg_global) - yg_global(1:nyg_global-1) )
    dzmin = Minval ( zg_global(2:nzg_global) - zg_global(1:nzg_global-1) )

    ! total domain size
    Lx = x_global(nx_global) - x_global(1)
    Ly = y_global(ny_global) - y_global(1)
    Lz = z_global(nz_global) - z_global(1)

    ! SEM inflow profile + eddy population (no-op unless x_bc_type==1 and
    ! inflow_type==1): must come after y_global/z_global/Ly/Lz are set
    Call init_inflow

    ! in-situ TI-profile rescaling setup (no-op unless ti_rescale_active==1); must come after init_inflow (needs prof_R11/22/33 already read)
    Call init_ti_rescale

    ! IBM setup: must come after xg/yg/zg/dxmin/dymin/dzmin are set
    ! (compute_normal_at_face_* divides by grid spacings).
    If ( ibm_input_mode >= 1 ) Then
       Call readSDF
    End If

    ! Boundary conditions
    ! local velocity, initial z-planes
    Allocate ( buffer_ui(nx,nyg,2:3), buffer_vi(nxg,ny,2:3), buffer_wi(nxg,nyg), buffer_ci(nxg,nyg,2:3) )
    ! local velocity, ending  z-planes
    Allocate ( buffer_ue(nx,nyg),     buffer_ve(nxg,ny),     buffer_we(nxg,nyg), buffer_ce(nxg,nyg) )
    ! local pressure z-plane
    Allocate ( buffer_p(2:nxg-1,2:nyg-1) ) 

    ! Interior communications; 3rd dim: 1=+z exchange, 2=-z exchange, issued concurrently
    Allocate ( buffer_us(nx ,nyg,2), buffer_ur(nx ,nyg,2) )
    Allocate ( buffer_vs(nxg, ny,2), buffer_vr(nxg, ny,2) )
    Allocate ( buffer_ws(nxg,nyg,2), buffer_wr(nxg,nyg,2) )
    Allocate ( buffer_ps(2:nxg-1,2:nyg-1), buffer_pr(2:nxg-1,2:nyg-1) )

    ! Fourier transform
    If ( myid==0 ) Write(*,'(A)') ' [4/4] Initializing FFT (2decomp&fft pencil transposes) ...'

    ! Fourier constant grid spacing
    dx = dxmin
    dz = dzmin

    ! Length for periodic domain; periodic-vs-DCT-IV point-count rationale
    If ( x_bc_type == 0 ) Then
       Lxp        = Lx - dx
       nxp_global = nxm_global - 1
    Else
       Lxp        = Lx
       nxp_global = nxm_global
    End If
    Lzp = Lz - dz
    nzp_global = nzm_global - 1

    ! global indices for fourier modes starting from 0
    mx_global = nxp_global - 1
    mz_global = nzp_global - 1

    ! local z-slab count for physical-space rhs_p ghost fill (apply_periodic_xz_pressure), independent of the pencil-transpose chain below
    Call z_periodic_partner(is_first_z, is_last_z, z_partner)
    nzp = nzm
    If ( is_last_z ) nzp = nzm - 1

    ! Poisson pencil grid: sized to the Fourier-transform grid (nxp_global x nym_global x nzp_global), not the face-point grid decomp_init_pencil used
    Call decomp_init_poisson_pencil( Int(nxp_global,Int32), nym_global, Int(nzp_global,Int32) )

    ! this rank's local mode-count range in the y-pencil (post transpose chain)
    mx = decomp_poisson%ysz(1) - 1
    mz = decomp_poisson%ysz(3) - 1

    ! y-first layout: rhs_p_hat(2:nyg-1, 0:mx, 0:mz) keeps each y-pencil
    ! contiguous so Zgtsv can operate in-place without a temporary copy.
    Allocate ( rhs_p_hat ( 2:nyg-1, 0:mx, 0:mz ) )

    ! pencil work arrays for the transpose chain
    Allocate ( poisson_y_r ( decomp_poisson%ysz(1), decomp_poisson%ysz(2), decomp_poisson%ysz(3) ) )
    Allocate ( poisson_x_r ( decomp_poisson%xsz(1), decomp_poisson%xsz(2), decomp_poisson%xsz(3) ) )
    Allocate ( poisson_x_c ( decomp_poisson%xsz(1), decomp_poisson%xsz(2), decomp_poisson%xsz(3) ) )
    Allocate ( poisson_y_c ( decomp_poisson%ysz(1), decomp_poisson%ysz(2), decomp_poisson%ysz(3) ) )
    Allocate ( poisson_z_c ( decomp_poisson%zsz(1), decomp_poisson%zsz(2), decomp_poisson%zsz(3) ) )

    ! local (non-MPI) batched complex 1-D FFT in z, dimension 3 (fully local in the z-pencil)
    plan_fz_fwd = fftw_plan_many_dft( 1_C_INT, [Int(nzp_global,C_INT)],                          &
             Int(decomp_poisson%zsz(1)*decomp_poisson%zsz(2),C_INT),                             &
             poisson_z_c, [Int(nzp_global,C_INT)], Int(decomp_poisson%zsz(1)*decomp_poisson%zsz(2),C_INT), 1_C_INT, &
             poisson_z_c, [Int(nzp_global,C_INT)], Int(decomp_poisson%zsz(1)*decomp_poisson%zsz(2),C_INT), 1_C_INT, &
             FFTW_FORWARD, FFTW_MEASURE )
    plan_fz_inv = fftw_plan_many_dft( 1_C_INT, [Int(nzp_global,C_INT)],                          &
             Int(decomp_poisson%zsz(1)*decomp_poisson%zsz(2),C_INT),                             &
             poisson_z_c, [Int(nzp_global,C_INT)], Int(decomp_poisson%zsz(1)*decomp_poisson%zsz(2),C_INT), 1_C_INT, &
             poisson_z_c, [Int(nzp_global,C_INT)], Int(decomp_poisson%zsz(1)*decomp_poisson%zsz(2),C_INT), 1_C_INT, &
             FFTW_BACKWARD, FFTW_MEASURE )

    If ( x_bc_type == 0 ) Then

       ! ---- periodic x: local batched complex 1-D FFT in x, dimension 1 (fully local in the x-pencil) ----------
       plan_fx_fwd = fftw_plan_many_dft( 1_C_INT, [Int(nxp_global,C_INT)],                    &
                Int(decomp_poisson%xsz(2)*decomp_poisson%xsz(3),C_INT),                       &
                poisson_x_c, [Int(nxp_global,C_INT)], 1_C_INT, Int(nxp_global,C_INT),         &
                poisson_x_c, [Int(nxp_global,C_INT)], 1_C_INT, Int(nxp_global,C_INT),         &
                FFTW_FORWARD, FFTW_MEASURE )
       plan_fx_inv = fftw_plan_many_dft( 1_C_INT, [Int(nxp_global,C_INT)],                    &
                Int(decomp_poisson%xsz(2)*decomp_poisson%xsz(3),C_INT),                       &
                poisson_x_c, [Int(nxp_global,C_INT)], 1_C_INT, Int(nxp_global,C_INT),         &
                poisson_x_c, [Int(nxp_global,C_INT)], 1_C_INT, Int(nxp_global,C_INT),         &
                FFTW_BACKWARD, FFTW_MEASURE )

       ! global Fourier coeficients with modified wave-number for the second derivative
       Allocate ( kxx(0:mx_global), kzz(0:mz_global) )
       kxx = 0d0
       kzz = 0d0
       Do i = 0, Ceiling( Real(nxp_global)/2d0 )
          kxx(i) = 2d0*( dcos(2d0*pi*Real(i,8)/Real(nxp_global,8)) - 1d0 )/dx**2d0
       End do
       Do i = Ceiling( Real(nxp_global)/2d0 )+1, mx_global
          kxx(i) = 2d0*( dcos(2d0*pi*Real(-nxp_global+i,8)/Real(nxp_global,8)) - 1d0 )/dx**2d0
       End do

       Do k = 0, Ceiling( Real(nzp_global)/2d0 )
          kzz(k) = 2d0*( dcos(2d0*pi*Real(k,8)/Real(nzp_global,8)) - 1d0 )/dz**2d0
       End do
       Do k = Ceiling( Real(nzp_global)/2d0 )+1, mz_global
          kzz(k) = 2d0*( dcos(2d0*pi*Real(-nzp_global+k,8)/Real(nzp_global,8)) - 1d0 )/dz**2d0
       End do

    Else

       ! Inflow/outflow x: local (non-MPI) batched DCT-IV in x, dimension 1 (fully local in the x-pencil); DCT-IV (FFTW_REDFT11) is self-inverse (up to 1/(2N) FFTW normalisation), so reused forward/backward
       plan_dct = fftw_plan_many_r2r( 1_C_INT, [Int(nxp_global,C_INT)],                       &
                  Int(decomp_poisson%xsz(2)*decomp_poisson%xsz(3),C_INT),                     &
                  poisson_x_r, [Int(nxp_global,C_INT)], 1_C_INT, Int(nxp_global,C_INT),       &
                  poisson_x_r, [Int(nxp_global,C_INT)], 1_C_INT, Int(nxp_global,C_INT),       &
                  [FFTW_REDFT11], FFTW_MEASURE )

       ! DCT-IV (quarter-wave) modified wavenumber; BC/null-space rationale
       Allocate ( kxx(0:mx_global), kzz(0:mz_global) )
       Do i = 0, mx_global
          kxx(i) = 2d0*( dcos( pi*Real(2*i+1,8)/(2d0*Real(nxp_global,8)) ) - 1d0 )/dx**2d0
       End Do
       kzz = 0d0
       Do k = 0, Ceiling( Real(nzp_global)/2d0 )
          kzz(k) = 2d0*( dcos(2d0*pi*Real(k,8)/Real(nzp_global,8)) - 1d0 )/dz**2d0
       End do
       Do k = Ceiling( Real(nzp_global)/2d0 )+1, mz_global
          kzz(k) = 2d0*( dcos(2d0*pi*Real(-nzp_global+k,8)/Real(nzp_global,8)) - 1d0 )/dz**2d0
       End do

    End If

    ! Wall-normal (y) periodicity (y_bc_type==0): local batched complex 1-D FFT
    ! in y, dimension 2 (fully local in the y-pencil -- y is never MPI-split).
    ! Same "last interior cell is a redundant duplicate of the first" convention
    ! as the periodic x/z directions (nxp_global=nxm_global-1): the periodic FFT
    ! length is nyp = decomp_poisson%ysz(2)-1 = nym_global-1, NOT the full
    ! nym_global -- the excluded last cell (index nyg-1) is reconstructed by a
    ! plain copy from index 2 after the solve (see solve_poisson_equation),
    ! exactly mirroring apply_periodic_xz_pressure's rhs_p(nxg-1,:,:)=rhs_p(2,:,:).
    ! Batched over the stride-1 x-index only (dim1); the z-index (dim3) is
    ! looped over at execute time in solve_poisson_equation, since
    ! fftw_plan_many_dft can only batch over one additional stride/dist pair.
    If ( y_bc_type == 0 ) Then

       plan_fy_fwd = fftw_plan_many_dft( 1_C_INT, [Int(decomp_poisson%ysz(2)-1,C_INT)],                      &
                Int(decomp_poisson%ysz(1),C_INT),                                                            &
                poisson_y_c, [Int(decomp_poisson%ysz(2)-1,C_INT)], Int(decomp_poisson%ysz(1),C_INT), 1_C_INT, &
                poisson_y_c, [Int(decomp_poisson%ysz(2)-1,C_INT)], Int(decomp_poisson%ysz(1),C_INT), 1_C_INT, &
                FFTW_FORWARD, FFTW_MEASURE )
       plan_fy_inv = fftw_plan_many_dft( 1_C_INT, [Int(decomp_poisson%ysz(2)-1,C_INT)],                      &
                Int(decomp_poisson%ysz(1),C_INT),                                                            &
                poisson_y_c, [Int(decomp_poisson%ysz(2)-1,C_INT)], Int(decomp_poisson%ysz(1),C_INT), 1_C_INT, &
                poisson_y_c, [Int(decomp_poisson%ysz(2)-1,C_INT)], Int(decomp_poisson%ysz(1),C_INT), 1_C_INT, &
                FFTW_BACKWARD, FFTW_MEASURE )

       ! global Fourier coefficients with modified wave-number for the second derivative; uniform y grid (grid_type=1 enforced above)
       Allocate ( kyy(0:decomp_poisson%ysz(2)-2) )
       kyy = 0d0
       dy  = y_global(2) - y_global(1)
       Do j = 0, Ceiling( Real(decomp_poisson%ysz(2)-1)/2d0 )
          kyy(j) = 2d0*( dcos(2d0*pi*Real(j,8)/Real(decomp_poisson%ysz(2)-1,8)) - 1d0 )/dy**2d0
       End Do
       Do j = Ceiling( Real(decomp_poisson%ysz(2)-1)/2d0 )+1, decomp_poisson%ysz(2)-2
          kyy(j) = 2d0*( dcos(2d0*pi*Real(-(decomp_poisson%ysz(2)-1)+j,8)/Real(decomp_poisson%ysz(2)-1,8)) - 1d0 )/dy**2d0
       End Do

    End If

    ! Tridiagonal linear solver
    If ( myid==0 ) Write(*,*) 'initializing pressure solver...'
    Allocate ( pivot(nyg) )
    Allocate ( Dyy(2:nyg-1,2:nyg-1) )
    Allocate ( D(2:nyg-1), DL(2:nyg-2), DU(2:nyg-2) )
     
    ! second derivative matrix for pressure (full data in y assumed)
    Dyy = 0d0
    Do j=3,nyg-2

       a = 1d0/( y(j)-y(j-1) )/( yg(j+1) - yg(j) )
       b = 1d0/( y(j)-y(j-1) )*( -1d0/( yg(j+1) - yg(j) ) -1d0/( yg(j) - yg(j-1) ) )
       c = 1d0/( y(j)-y(j-1) )/( yg(j) - yg(j-1) ) 

       Dyy(j,j+1) = a
       Dyy(j,j-1) = c 
       Dyy(j,j  ) = b

    End Do

    ! Boundary conditions for pressure (full data in y assumed)
    j = 2
    a = 1d0/( y(j)-y(j-1) )/( yg(j+1) - yg(j) )
    b = 1d0/( y(j)-y(j-1) )*( -1d0/( yg(j+1) - yg(j) ) -1d0/( yg(j) - yg(j-1) ) )
    c = 1d0/( y(j)-y(j-1) )/( yg(j) - yg(j-1) ) 
    ! Dirichlet in V: p(1)==p(2) 
    Dyy(2,2)   = b + c 
    Dyy(2,3)   = a
    coef_bc_1  = c

    j = nyg-1
    a = 1d0/( y(j)-y(j-1) )/( yg(j+1) - yg(j) )
    b = 1d0/( y(j)-y(j-1) )*( -1d0/( yg(j+1) - yg(j) ) -1d0/( yg(j) - yg(j-1) ) )
    c = 1d0/( y(j)-y(j-1) )/( yg(j) - yg(j-1) )     
    ! Dirichlet in V: p(nyg)==p(nyg-1) 
    Dyy(nyg-1,nyg-1) = a + b
    Dyy(nyg-1,nyg-2) = c
    coef_bc_2        = a

    Allocate ( bc_1(2:nxg-1,2:nzg-1), bc_2(2:nxg-1,2:nzg-1) )
    Allocate ( bc_1_hat(0:mx,0:mz),   bc_2_hat(0:mx,0:mz)   )

    ! some parameters for linear solver
    nr   = nym
    nrhs = 1

    ! interpolation weights
    in1 = 1
    in2 = 2 
    If ( in2==1 ) Write(*,*) 'Conservative interpolations'
    Allocate ( weight_y_0(ny), weight_y_1(ny) )
    weight_y_0 = ( yg(2:nyg) - y(1:ny) ) / ( yg(2:nyg) - yg(1:nyg-1)  )
    weight_y_1 = 1d0 - weight_y_0

    ! Runge-Kutta 3
    If ( myid==0 ) Write(*,*) 'initializing time integration...'
    Allocate( rk_coef(3,3), rk_t(3) )

    rk_t(1)      =  8d0/15d0
    rk_t(2)      =  2d0/3d0
    rk_t(3)      =  1d0

    rk_coef      =  0d0
    rk_coef(1,1) =  8d0/15d0
    rk_coef(2,1) =  1d0/4d0
    rk_coef(2,2) =  5d0/12d0
    rk_coef(3,1) =  1d0/4d0
    rk_coef(3,2) =  0d0
    rk_coef(3,3) =  3d0/4d0

    Allocate ( Fu1 ( 2:nx-1,  2:nyg-1, 2:nzg-1 ) )
    Allocate ( Fu2 ( 2:nx-1,  2:nyg-1, 2:nzg-1 ) )
    Allocate ( Fu3 ( 2:nx-1,  2:nyg-1, 2:nzg-1 ) )

    Allocate ( Fv1 ( 2:nxg-1,  2:ny-1, 2:nzg-1 ) )
    Allocate ( Fv2 ( 2:nxg-1,  2:ny-1, 2:nzg-1 ) )
    Allocate ( Fv3 ( 2:nxg-1,  2:ny-1, 2:nzg-1 ) )

    Allocate ( Fw1 ( 2:nxg-1,  2:nyg-1, 2:nz-1 ) )
    Allocate ( Fw2 ( 2:nxg-1,  2:nyg-1, 2:nz-1 ) )
    Allocate ( Fw3 ( 2:nxg-1,  2:nyg-1, 2:nz-1 ) )

    dPdy      = 0d0

    ! Eddy viscosity
    Allocate( nu_t     (1:nxg,1:nyg,1:nzg) )

    !	Define eddy viscosity to 0
    nu_t = 0d0

    ! Suspended sediment
    If ( sediment_flag >= 1 ) Then
       Allocate( Cscal  (1:nxg, 1:nyg, 1:nzg) )
       Allocate( Cscal_o(1:nxg, 1:nyg, 1:nzg) )
       Allocate( Fcs1(2:nxg-1, 2:nyg-1, 2:nzg-1) )
       Allocate( Fcs2(2:nxg-1, 2:nyg-1, 2:nzg-1) )
       Allocate( Fcs3(2:nxg-1, 2:nyg-1, 2:nzg-1) )
       If ( restart == 1 .And. scalar_restart == 1 ) Then
          Call read_scalar_restart
       Else
          Call init_scalar_profile
       End If
       Call compute_settling_velocity
       If (myid == 0) Write(*,'(A,E12.4,A)') '   ws (Soulsby 1997) = ', ws, ' m/s'
    End If

    ! Boussinesq temperature
    If ( boussinesq_flag >= 1 ) Then
       Allocate( Tscal  (1:nxg, 1:nyg, 1:nzg) )
       Allocate( Tscal_o(1:nxg, 1:nyg, 1:nzg) )
       Allocate( Ft1(2:nxg-1, 2:nyg-1, 2:nzg-1) )
       Allocate( Ft2(2:nxg-1, 2:nyg-1, 2:nzg-1) )
       Allocate( Ft3(2:nxg-1, 2:nyg-1, 2:nzg-1) )
       If ( restart == 1 .And. scalar_restart == 1 ) Then
          Call read_temperature_restart
       Else
          Call init_temperature_profile
       End If
    End If

    ! Wall-models

    Delta = ( (x(2)-x(1))*(y(2)-y(1))*(z(2)-z(1)) )**(1d0/3d0)

    ! ui = alpha_i dui/dy
    Allocate( alpha_x(1:nx ,1:2,1:nzg) )
    Allocate( alpha_y(1:nxg,1:2,1:nzg) )
    Allocate( alpha_z(1:nxg,1:2,1:nz ) )    

    alpha_x = 0d0
    alpha_y = 0d0
    alpha_z = 0d0

    If ( boussinesq_flag >= 1 ) Then
       Allocate( alpha_T(1:nxg,1:2,1:nzg) )
       alpha_T = 0d0
    End If

    ! Rough-EQWM matching height: smallest interior j (bottom) / largest interior j
    ! (top) whose distance to the wall clears 20*z0, so the log-law u_tau/theta_tau
    ! solves in wallmodel.f90 aren't sampled from inside the roughness sublayer just
    ! because the near-wall grid happens to be fine (same threshold as wallmodel.f90's
    ! MATCH_RATIO_MIN -- kept as a literal here to avoid a new module dependency).
    If ( flat_wall_model_flag == 2 ) Then
       j_match_ylo = 2
       If ( bc_face_ylo /= 2 .And. z0_ylo > 0d0 ) Then
          Do j = 2, nyg/2
             j_match_ylo = j
             If ( yg(j)/z0_ylo >= 20d0 ) Exit
          End Do
          If ( myid==0 .And. yg(j_match_ylo)/z0_ylo < 20d0 ) Write(*,'(A)') &
             ' WARNING: no grid point within the lower half-channel clears y/z0_ylo >= 20 for the rough EQWM matching height'
       End If

       j_match_yhi = nyg-1
       If ( bc_face_yhi /= 2 .And. z0_yhi > 0d0 ) Then
          Do j = nyg-1, nyg/2, -1
             j_match_yhi = j
             If ( (Ly - yg(j))/z0_yhi >= 20d0 ) Exit
          End Do
          If ( myid==0 .And. (Ly-yg(j_match_yhi))/z0_yhi < 20d0 ) Write(*,'(A)') &
             ' WARNING: no grid point within the upper half-channel clears y/z0_yhi >= 20 for the rough EQWM matching height'
       End If

       If ( myid == 0 ) Then
          Write(*,'(A,I0,A,E12.4)') '   rough EQWM matching height: j_match_ylo = ', j_match_ylo, &
             ', y_match_lo = ', yg(j_match_ylo)
          Write(*,'(A,I0,A,E12.4)') '   rough EQWM matching height: j_match_yhi = ', j_match_yhi, &
             ', y_match_hi = ', Ly - yg(j_match_yhi)
       End If
    End If

    ! Push host values once for scalars `!$acc declare create`d in global.f90 (needed by !$acc routine seq procedures)
    !$acc update device(nx,ny,nz,nxg,nyg,nzg,nu,pi,inflow_type,inflow_Uconst,sem_n_eddies,sem_length_scale,sem_seed,sem_eddy_placement,sem_use_esem,sem_divergence_free)
    !$acc update device(flat_wall_model_flag,z0_ylo,z0_yhi,z0h_ylo,z0h_yhi,j_match_ylo,j_match_yhi)
    !$acc update device(boussinesq_flag,beta_T,T_ref,Pr,Pr_t,grav,ibm_T_bc_type,ibm_T_wall)
    !$acc update device(T_bc_bot,T_bc_top,T_wall_bot,T_wall_top)
    ! OpenACC: static grid arrays copied once, scratch arrays created once, never re-transferred
    !$acc enter data copyin(y,yg,z,zg,weight_y_0,weight_y_1)
    ! x,xg needed on device by compute_rhs_scalar_core (non-uniform-grid-style stencil, even though x itself is uniform)
    !$acc enter data copyin(x,xg)
    !$acc enter data create(term,term_1,term_2)
    ! Evolving per-substage fields: allocated here (create, not copyin) so time_integration.f90's per-RK-substage update device/host calls have a target
    !$acc enter data create(U,V,W,nu_t,Fu1,Fv1,Fw1,Fu2,Fv2,Fw2,Fu3,Fv3,Fw3)
    ! One-time initial sync: enter data create only allocates, so push host U,V,W or step 1 runs on uninitialized device data
    !$acc update device(U,V,W)
    ! rhs_p, RK3 base-state snapshots (Uo/Vo/Wo), and Robin-BC slip-length coefficients (host-written, GPU-read)
    !$acc enter data create(rhs_p,Uo,Vo,Wo,alpha_x,alpha_y,alpha_z)
    ! Boussinesq temperature: stays device-resident end-to-end so compute_rhs_v can read it directly
    If ( boussinesq_flag >= 1 ) Then
       !$acc enter data create(Tscal,Tscal_o,Ft1,Ft2,Ft3)
       !$acc update device(Tscal,Tscal_o)
       !$acc enter data create(alpha_T)
       !$acc update device(alpha_T)
    End If
    ! Sediment scalar: device residency required by the shared compute_rhs_scalar_core (see scalar_transport.f90)
    If ( sediment_flag >= 1 ) Then
       !$acc enter data create(Cscal,Cscal_o,Fcs1,Fcs2,Fcs3)
       !$acc update device(Cscal,Cscal_o)
    End If
#ifdef GPU_POISSON
    ! rhs_p_hat stays device-resident end-to-end through the whole Poisson stage
    !$acc enter data create(rhs_p_hat)
#endif

    ! Done
    Call Mpi_barrier(MPI_COMM_WORLD,ierr)

    ! Measure time
    time1 = MPI_WTIME()
    
  End Subroutine initialize

  ! Set the initial scalar concentration profile; C_ic_type 0-3 profile formulas
  Subroutine init_scalar_profile

    Integer(Int32) :: jj
    Real   (Int64) :: y_c, y0, a_ref, Rouse_Z, u_tau_ic
    Real   (Int64), Parameter :: kappa_sc = 0.41d0

    ! Bottom wall y-position (ghost cell boundary)
    y0 = yg(1)

    Select Case (C_ic_type)

    Case (0)   ! Uniform
       Cscal = C_ref

    Case (1)   ! Rouse equilibrium profile
       ! Estimate u_tau the same way generateIC does
       If ( Abs(dPdx) > 0d0 ) Then
          If ( bc_face_yhi == 2 ) Then
             u_tau_ic = Sqrt( Abs(dPdx) * Ly_i )
          Else
             u_tau_ic = Sqrt( Abs(dPdx) * 0.5d0 * Ly_i )
          End If
       Else
          u_tau_ic = 0.05d0 * Abs(Utarget)
       End If
       ! Reference height a: first interior cell centre above bottom wall
       a_ref = 0.5d0 * ( yg(2) - yg(1) )   ! half first cell height
       a_ref = Max( a_ref, 1d-10 )
       If ( u_tau_ic > 0d0 .And. ws > 0d0 ) Then
          Rouse_Z = ws / ( kappa_sc * u_tau_ic )
       Else
          Rouse_Z = 0d0
       End If
       If (myid == 0) Write(*,'(A,E12.4,A,E12.4)') &
          '   Rouse Z = ', Rouse_Z, '  a_ref = ', a_ref
       Do jj = 2, nyg-1
          y_c = yg(jj) - y0   ! distance above bottom wall
          y_c = Max( y_c, a_ref )
          If ( Rouse_Z > 0d0 ) Then
             Cscal(:,jj,:) = C_ref * ( (a_ref/y_c) * (Ly_i - y_c)/(Ly_i - a_ref) )**Rouse_Z
          Else
             Cscal(:,jj,:) = C_ref
          End If
       End Do
       ! Ghost cells: copy nearest interior value
       Cscal(:,1,:)   = Cscal(:,2,:)
       Cscal(:,nyg,:) = Cscal(:,nyg-1,:)

    Case (2)   ! Linear ramp C_ref (bed) -> 0 (top)
       Do jj = 2, nyg-1
          y_c = yg(jj) - y0
          Cscal(:,jj,:) = C_ref * Max( 0d0, 1d0 - y_c/Ly_i )
       End Do
       Cscal(:,1,:)   = Cscal(:,2,:)
       Cscal(:,nyg,:) = 0d0

    Case (3)   ! Slab: C_ref below C_ic_height, 0 above
       Do jj = 2, nyg-1
          y_c = yg(jj) - y0
          If ( y_c <= C_ic_height ) Then
             Cscal(:,jj,:) = C_ref
          Else
             Cscal(:,jj,:) = 0d0
          End If
       End Do
       Cscal(:,1,:)   = Cscal(:,2,:)
       Cscal(:,nyg,:) = 0d0

    Case Default
       Cscal = C_ref

    End Select

    Cscal_o = Cscal

    ! Zero scalar concentration inside IBM solid cells so that the IC value
    ! never persists inside solid and pollutes adjacent fluid via diffusion.
    If ( ibm_input_mode >= 1 .And. Allocated(phi) ) Then
       Where ( phi(2:nxg-1, 2:nyg-1, 2:nzg-1) <= 0d0 )
          Cscal(2:nxg-1, 2:nyg-1, 2:nzg-1) = 0d0
       End Where
       Cscal_o = Cscal
    End If

    If (myid == 0) Write(*,'(A,I1,A,E12.4)') &
       '   C_ic_type = ', C_ic_type, '  Max Cscal = ', MaxVal(Cscal)

  End Subroutine init_scalar_profile

  ! Set the initial temperature profile; T_ic_type 0=uniform (T_ref), 1=linear gradient in y from the bottom wall
  Subroutine init_temperature_profile

    Integer(Int32) :: jj
    Real   (Int64) :: y0

    y0 = yg(1)

    Select Case (T_ic_type)

    Case (1)   ! Linear gradient: T = T_ref + T_ic_grad*(y-y0)
       Do jj = 1, nyg
          Tscal(:,jj,:) = T_ref + T_ic_grad * ( yg(jj) - y0 )
       End Do

    Case Default   ! Uniform
       Tscal = T_ref

    End Select

    Tscal_o = Tscal

    ! Zero temperature perturbation inside IBM solid cells so the IC never persists inside solid and pollutes adjacent fluid via diffusion
    If ( ibm_input_mode >= 1 .And. Allocated(phi) ) Then
       Where ( phi(2:nxg-1, 2:nyg-1, 2:nzg-1) <= 0d0 )
          Tscal(2:nxg-1, 2:nyg-1, 2:nzg-1) = T_ref
       End Where
       Tscal_o = Tscal
    End If

    If (myid == 0) Write(*,'(A,I1,A,E12.4)') &
       '   T_ic_type = ', T_ic_type, '  Max Tscal = ', MaxVal(Tscal)

  End Subroutine init_temperature_profile

End Module initialization
