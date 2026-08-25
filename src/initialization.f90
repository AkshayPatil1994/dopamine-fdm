! Module with initialization of global variables
Module initialization

  ! Modules
  Use, Intrinsic :: iso_c_binding
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use mpi
  Use input_output
  Use ibmSetup
  Use scalar_transport, Only : compute_settling_velocity
  Use synthetic_eddy_method, Only : init_inflow, init_ti_rescale

  ! prevent implicit typing
  Implicit None

  ! declarations
Contains

  ! Initialize everything
  Subroutine initialize
    
    Integer(Int32) :: i, j, k, kk, nzpe, pos, ipos, nze, nzme
    Real   (Int64) :: dy1, dy2, det, a, b, c, r, Qflow_ref
    Integer(Int32), Dimension(:,:), Allocatable :: A_kmodes, A_kmodes_local

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
    If ( myid==0 ) Write(*,'(A)') ' [1/4] Reading input parameters ...'
    Call read_input_parameters

#ifdef GPU_POISSON
    ! single-GPU cuFFT Poisson solve (periodic or DCT-IV) only supports nprocs==1
    If ( nprocs /= 1 ) Stop 'ERROR: GPU_POISSON build only supports nprocs=1 (single-GPU); rebuild without ENABLE_GPU for multi-rank runs'
    If ( x_bc_type /= 0 .And. x_bc_type /= 1 ) Stop 'ERROR: GPU_POISSON build only supports x_bc_type=0 or 1'
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
    
    ! grid definitions
    Allocate (  k1_global(0:nprocs-1),  k2_global(0:nprocs-1) )
    Allocate ( kg1_global(0:nprocs-1), kg2_global(0:nprocs-1) )

    ! restrictions for FFTW mapping
    If ( Mod( nx_global   , 2     )/=0 ) Stop 'Error: nx must be even'
    If ( Mod( nz_global   , 2     )/=0 ) Stop 'Error: nz must be even'
    If ( Mod( nz_global-2 , nprocs)/=0 ) Stop 'nz-2 should be divisible by number of processors'

    ! number of interior z-planes per processor based on fftw decomposition
    nslices_z = Nint( Real((nz_global-2))/Real(nprocs) ) 
    
    ! restriction for MPI boundaries
    If ( nslices_z<2 ) Stop 'Error: nslices_z must be at least 2' 

    ! domain decomposition. Must be consistent with fftw
    Do i = 0, nprocs-1
       ! range index for faces in each processor
       k1_global(i)  = i*nslices_z  + 1
       k2_global(i)  = k1_global(i) + nslices_z + 1
       ! range index for centers in each processor
       kg1_global(i) = i*nslices_z   + 1
       kg2_global(i) = kg1_global(i) + nslices_z + 1
    End Do    

    ! remaining planes in last processor
    k2_global (nprocs-1) = nz_global 
    kg2_global(nprocs-1) = nz_global + 1

    ! face points
    nx = nx_global
    ny = ny_global
    nz = k2_global(myid) - k1_global(myid) + 1 

    ! middle points
    nxm_global = nx_global - 1
    nym_global = ny_global - 1
    nzm_global = nz_global - 1

    nxm = nx - 1
    nym = ny - 1
    nzm = kg2_global(myid) - kg1_global(myid) + 1 - 2  
    
    ! middle points + ghost cells
    nxg_global = nxm_global + 2
    nyg_global = nym_global + 2
    nzg_global = nzm_global + 2

    nxg = nxm + 2
    nyg = nym + 2
    nzg = kg2_global(myid) - kg1_global(myid) + 1 

    ! size for last proccesor nz and nzm -> nze and nzme
    nze  = nz
    nzme = nzm
    Call Mpi_bcast (  nze,1,MPI_integer,nprocs-1,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( nzme,1,MPI_integer,nprocs-1,MPI_COMM_WORLD,ierr )
   
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

    If (myid == 0) Then
       Allocate (Uoo (    nx,  nym+2, nzme+2) ) ! z-planes modified for I/O
       Allocate (Voo ( nxm+2,     ny, nzme+2) )
       Allocate (Woo ( nxm+2,  nym+2,    nze) )
       Allocate (Poo ( nxm+2,  nym+2, nzme+2) )
    End If

    ! Auxiliary arrays
    Allocate ( term   ( nxg, nyg, nzm+2 ) ) 
    Allocate ( term_1 ( nxg, nyg, nzm+2 ) ) 
    Allocate ( term_2 ( nxg, nyg, nzm+2 ) ) 

    ! RHS: interior points only
    Allocate ( rhs_uo ( 2:nx-1,  2:nyg-1, 2:nzg-1 ) ) 
    Allocate ( rhs_vo ( 2:nxg-1, 2:ny-1,  2:nzg-1 ) )
    Allocate ( rhs_wo ( 2:nxg-1, 2:nyg-1, 2:nz-1  ) )
    Allocate ( rhs_p  ( 2:nxg-1, 2:nyg-1, 2:nzg   ) ) ! ONE EXTRA PLANE IN Z FOR GHOST CELL
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
    x = x_global
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
    xm = xm_global
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

    xg = xg_global
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

    ! Interior communications
    Allocate ( buffer_us(nx ,nyg), buffer_ur(nx ,nyg) )
    Allocate ( buffer_vs(nxg, ny), buffer_vr(nxg, ny) )
    Allocate ( buffer_ws(nxg,nyg), buffer_wr(nxg,nyg) )
    Allocate ( buffer_ps(2:nxg-1,2:nyg-1), buffer_pr(2:nxg-1,2:nyg-1) ) 

    ! Fourier transform
    If ( myid==0 ) Write(*,'(A)') ' [4/4] Initializing FFT (FFTW-MPI) ...'
    ! initialize MPI FFTW
    Call fftw_mpi_init()

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

    ! Get local sizes: local x-direction size (x is never domain-decomposed, in either path)
    nxp = nxp_global
    mx  =  mx_global

    If ( x_bc_type == 0 ) Then

       ! ---- periodic x: unchanged joint 2-D complex FFT path ----------
       ! local data size in z direction (note dimension reversal)
       alloc_local = fftw_mpi_local_size_2d(nzp_global, nxp_global, MPI_COMM_WORLD, nzp, local_k_offset)
       mz  = nzp - 1

       ! sanity check and restrictions in fftw
       If ( (nzp/=nzm .And. myid/=nprocs-1) .Or. (nzp/=nzm-1 .And. myid==nprocs-1) ) Then
          Write(*,*) nzp,nzm
          Stop 'Error: something wrong in FFTW size'
       End If

       ! allocate variables
       cplane_fft = fftw_alloc_complex(alloc_local)
       Call c_f_pointer(cplane_fft,plane,[nxp,nzp])
       plane_hat(0:,0:) => plane
       ! y-first layout: rhs_p_hat(2:nyg-1, 0:mx, 0:mz) keeps each y-pencil
       ! contiguous so Zgtsv can operate in-place without a temporary copy.
       Allocate ( rhs_p_hat ( 2:nyg-1, 0:mx, 0:mz ) )

       ! create MPI plan for forward DFT (note dimension reversal and transposed_out/in)
       plan_d = fftw_mpi_plan_dft_2d( nzp_global, nxp_global, plane, plane_hat,           &
                MPI_COMM_WORLD,  FFTW_FORWARD, ior(FFTW_MEASURE, FFTW_MPI_TRANSPOSED_OUT) )
       plan_i = fftw_mpi_plan_dft_2d( nzp_global, nxp_global, plane_hat, plane,           &
                MPI_COMM_WORLD, FFTW_BACKWARD, ior(FFTW_MEASURE, FFTW_MPI_TRANSPOSED_IN)  )

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

       ! Inflow/outflow x: DCT-IV (quarter-wave cosine transform); joint-vs-batched transform rationale
       alloc_local = fftw_mpi_local_size_many( 1_C_INT, [nzp_global], nxp_global,        &
                     FFTW_MPI_DEFAULT_BLOCK, MPI_COMM_WORLD, nzp, local_k_offset )
       mz  = nzp - 1

       If ( (nzp/=nzm .And. myid/=nprocs-1) .Or. (nzp/=nzm-1 .And. myid==nprocs-1) ) Then
          Write(*,*) nzp,nzm
          Stop 'Error: something wrong in FFTW size'
       End If

       cplane_fft = fftw_alloc_complex(alloc_local)
       Call c_f_pointer(cplane_fft,plane,[nxp,nzp])
       plane_hat(0:,0:) => plane
       Allocate ( rhs_p_hat ( 2:nyg-1, 0:mx, 0:mz ) )

       ! Batched (howmany=nxp_global) distributed 1-D complex FFT in z; no-transpose-flags rationale
       plan_d = fftw_mpi_plan_many_dft( 1_C_INT, [nzp_global], nxp_global,               &
                FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, plane, plane_hat,         &
                MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE )
       plan_i = fftw_mpi_plan_many_dft( 1_C_INT, [nzp_global], nxp_global,               &
                FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, plane_hat, plane,         &
                MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE )

       ! Local (non-MPI) DCT-IV plan for x, batched over local z-lines; DCT-IV (FFTW_REDFT11) is self-inverse (up to 1/(2N) FFTW normalisation), so reused forward/backward
       rplane_fft = fftw_alloc_real( int(nxp_global,C_SIZE_T) * int(nzp,C_SIZE_T) )
       Call c_f_pointer( rplane_fft, xplane, [nxp_global,nzp] )
       plan_dct = fftw_plan_many_r2r( 1_C_INT, [Int(nxp_global,C_INT)], Int(nzp,C_INT),  &
                  xplane, [Int(nxp_global,C_INT)], 1_C_INT, Int(nxp_global,C_INT),       &
                  xplane, [Int(nxp_global,C_INT)], 1_C_INT, Int(nxp_global,C_INT),       &
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

    ! FFTW+MPI local-to-global mode mapping (transposed_out/in); periodic-only applicability rationale
    If ( x_bc_type == 0 ) Then

       ! this needs (mz_global+1)*(mx_global+1)/nprocs to be an integer
       If ( Mod((mz_global+1)*(mx_global+1),nprocs)/=0 ) Stop 'Error: (mz_global+1)*(mx_global+1)/nprocs should be an integer'
       Allocate ( imode_map_fft(0:mx_global,0:mz) )
       Allocate ( kmode_map_fft(0:mx_global,0:mz) )
       Do i = 0, mx_global
          Do k = 0, mz
             pos = i + (mx_global+1)*k + (mz_global+1)*(mx_global+1)/nprocs*myid
             imode_map_fft(i,k) = Floor( Real(pos/(mz_global+1)) )
             kmode_map_fft(i,k) = Mod  ( pos, mz_global+1 )
          end Do
       End Do

       ! Sanity check for FFTW mapping
       Allocate(A_kmodes      (0:mx_global,0:mz_global))
       Allocate(A_kmodes_local(0:mx_global,0:mz_global))
       A_kmodes       = 0
       A_kmodes_local = 0
       Do i = 0, mx_global
          Do k = 0, mz
             A_kmodes_local( imode_map_fft(i,k), kmode_map_fft(i,k) ) =  A_kmodes_local(imode_map_fft(i,k), kmode_map_fft(i,k) ) + 1
          end Do
       End Do
       Call MPI_AllReduce(A_kmodes_local,A_kmodes,Int((mx_global+1)*(mz_global+1),Int32),MPI_integer,MPI_sum,MPI_COMM_WORLD,ierr)
       If ( Any(A_kmodes>1) .Or. Any(A_kmodes==0) ) Stop 'Error: wrong combination of nx, nz and processors'
       Deallocate(A_kmodes)
       Deallocate(A_kmodes_local)

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
       If (myid == 0) Then
          Allocate( Cscal_oo(1:nxg, 1:nyg, 1:kg2_global(nprocs-1)-kg1_global(nprocs-1)+1) )
       End If
       If ( restart == 1 .And. scalar_restart == 1 ) Then
          Call read_scalar_restart
       Else
          Call init_scalar_profile
       End If
       Call compute_settling_velocity
       If (myid == 0) Write(*,'(A,E12.4,A)') '   ws (Soulsby 1997) = ', ws, ' m/s'
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

    ! Push host values once for scalars `!$acc declare create`d in global.f90 (needed by !$acc routine seq procedures)
    !$acc update device(nx,ny,nz,nxg,nyg,nzg,nu,pi,inflow_type,inflow_Uconst,sem_n_eddies,sem_length_scale,sem_seed,sem_eddy_placement,sem_use_esem,sem_divergence_free)
    ! OpenACC: static grid arrays copied once, scratch arrays created once, never re-transferred
    !$acc enter data copyin(y,yg,z,zg,weight_y_0,weight_y_1)
    !$acc enter data create(term,term_1,term_2)
    ! Evolving per-substage fields: allocated here (create, not copyin) so time_integration.f90's per-RK-substage update device/host calls have a target
    !$acc enter data create(U,V,W,nu_t,Fu1,Fv1,Fw1,Fu2,Fv2,Fw2,Fu3,Fv3,Fw3)
    ! One-time initial sync: enter data create only allocates, so push host U,V,W or step 1 runs on uninitialized device data
    !$acc update device(U,V,W)
    ! rhs_p, RK3 base-state snapshots (Uo/Vo/Wo), and Robin-BC slip-length coefficients (host-written, GPU-read)
    !$acc enter data create(rhs_p,Uo,Vo,Wo,alpha_x,alpha_y,alpha_z)
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

End Module initialization
