!> Exact Reynolds stress transport budget for all six R_ij = <u_i' u_j'> (Pope 2000 Section 7.4, Eq. 7.176)
Module reynolds_stress_budget

  Use iso_fortran_env, Only : Int32, Int64
  Use global
  Use mpi

  Implicit None
  Private

  ! ── public interface ─────────────────────────────────────────────────────
  Public :: init_rsb, accumulate_rsb, output_rsb, finalize_rsb
  Public :: n_accum    ! exposed so main.f90 can gate output_rsb calls

  ! component index helpers: 6 independent symmetric-tensor components in order 1=11,2=22,3=33,4=12,5=13,6=23
  Integer(Int32), Parameter :: N_COMP = 6

  ! ── homogeneous direction flags ──────────────────────────────────────────
  Logical :: hom_x = .False., hom_z = .False.

  ! ── output grid sizes after homogeneous averaging ───────────────────────
  !  These are set in init_rsb and never change.
  Integer(Int32) :: out_nx, out_ny, out_nz   ! = nxm_global, nym_global, nzm_global
                                              !   or 1 along homogeneous dirs

  ! ── accumulation counter ─────────────────────────────────────────────────
  Integer(Int32) :: n_accum = 0        ! steps accumulated in current window
  Integer(Int32) :: nsamples = 0       ! windows written so far (across restarts)

  ! accumulator arrays (cell-centre, local z-slab): shape (nxm, nym_global, nzm) -- x and y always full since MPI decomposition is along z only
  !  1st moments
  Real(Int64), Allocatable :: acc_U(:,:,:), acc_V(:,:,:), acc_W(:,:,:)
  Real(Int64), Allocatable :: acc_P(:,:,:)

  !  2nd moments  (UiUj raw: UU, VV, WW, UV, UW, VW)
  Real(Int64), Allocatable :: acc_UiUj(:,:,:,:)    ! (6, nxm, nym_global, nzm)

  !  Velocity-gradient products for resolved dissipation and production
  !  (dUi/dxk)(dUj/dxk)  →  6 components × 3 grad directions
  Real(Int64), Allocatable :: acc_GijRes(:,:,:,:)   ! (6, nxm, nym_global, nzm)

  !  SGS: 6 components of  2 nu_t s_ij'   where s_ij' = 0.5(du_i'/dx_j + du_j'/dx_i)
  Real(Int64), Allocatable :: acc_GijSGS(:,:,:,:)   ! (6, nxm, nym_global, nzm)

  !  Pressure-strain raw: p * (dui/dxj + duj/dxi) for each of 6 components
  Real(Int64), Allocatable :: acc_PiStrain(:,:,:,:) ! (6, nxm, nym_global, nzm)

  !  Pressure-velocity for pressure diffusion: PU, PV, PW
  Real(Int64), Allocatable :: acc_PVel(:,:,:,:)     ! (3, nxm, nym_global, nzm)

  !  Triple correlations for turbulent diffusion: d<u_i'u_j'u_k'>/dx_k needs UiUjUk products; flux divergences stored per direction (y always; UiUj*V' for 6 pairs)
  Real(Int64), Allocatable :: acc_TijY(:,:,:,:)    ! (6, nxm, nym_global, nzm)
  Real(Int64), Allocatable :: acc_TijX(:,:,:,:)    ! (6, nxm, nym_global, nzm)  ! only if !hom_x
  Real(Int64), Allocatable :: acc_TijZ(:,:,:,:)    ! (6, nxm, nym_global, nzm)  ! only if !hom_z

Contains

  !> init_rsb: called once after grid allocation; parses rsb_hom_dir, allocates accumulator arrays, opens/appends output files depending on restart
  Subroutine init_rsb

    Character(200) :: hd
    Integer(Int32) :: ios, meta_unit, existing_samples
    Character(300) :: meta_path
    Logical        :: meta_exists

    If ( rsb_active /= 1 ) Return

    ! ── parse homogeneous direction string ───────────────────────────────
    hd = Trim(Adjustl(rsb_hom_dir))
    hom_x = ( Index(hd,'x') > 0 .Or. Index(hd,'X') > 0 )
    hom_z = ( Index(hd,'z') > 0 .Or. Index(hd,'Z') > 0 )

    out_nx = Merge(1, nxm_global, hom_x)
    out_ny = nym_global
    out_nz = Merge(1, nzm_global, hom_z)

    If ( myid == 0 ) Then
       Write(*,'(A)')       ' RSB: Reynolds stress budget enabled.'
       Write(*,'(A,A)')     '   hom_dir   = ', Trim(hd)
       Write(*,'(A,3I5)')   '   out shape = ', out_nx, out_ny, out_nz
       Write(*,'(A,I8)')    '   rsb_freq  = ', rsb_freq
       Write(*,'(A,I8)')    '   rsb_nstart= ', rsb_nstart

       ! ── ensure output directory exists ──────────────────────────────
       Call ensure_rsb_dir
    End If

    ! ── restart: read existing sample count from meta file ───────────────
    nsamples = 0
    If ( restart == 1 .And. myid == 0 ) Then
       meta_path = Trim(rsb_fileout) // '.meta'
       Inquire(file=Trim(meta_path), exist=meta_exists)
       If ( meta_exists ) Then
          Open(newunit=meta_unit, file=Trim(meta_path), status='old', &
               action='read', iostat=ios)
          If ( ios == 0 ) Then
             Call parse_meta_nsamples(meta_unit, existing_samples)
             Close(meta_unit)
             nsamples = existing_samples
             Write(*,'(A,I8,A)') '   RSB restart: found ', nsamples, &
                  ' existing samples; will append.'
          End If
       Else
          Write(*,'(A)') '   RSB restart: no meta file found; starting fresh.'
       End If
    End If
    ! Broadcast nsamples so all ranks are consistent
    Call MPI_Bcast(nsamples, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! ── allocate accumulator arrays (local z-slab) ────────────────────────
    Allocate( acc_U      (  nxm, nym_global, nzm) )
    Allocate( acc_V      (  nxm, nym_global, nzm) )
    Allocate( acc_W      (  nxm, nym_global, nzm) )
    Allocate( acc_P      (  nxm, nym_global, nzm) )
    Allocate( acc_UiUj   (6,nxm, nym_global, nzm) )
    Allocate( acc_GijRes (6,nxm, nym_global, nzm) )
    Allocate( acc_GijSGS (6,nxm, nym_global, nzm) )
    Allocate( acc_PiStrain(6,nxm,nym_global, nzm) )
    Allocate( acc_PVel   (3,nxm, nym_global, nzm) )
    Allocate( acc_TijY   (6,nxm, nym_global, nzm) )
    If ( .Not. hom_x ) Allocate( acc_TijX(6,nxm,nym_global,nzm) )
    If ( .Not. hom_z ) Allocate( acc_TijZ(6,nxm,nym_global,nzm) )

    Call zero_accumulators

    ! fresh run: truncate any leftover files
    If ( restart == 0 .And. myid == 0 ) Then
       Call truncate_rsb_files
    End If

  End Subroutine init_rsb


  !> accumulate_rsb: called every rsb_freq steps when istep >= rsb_nstart; interpolates staggered U/V/W to cell centres, computes gradients, accumulates raw moments
  Subroutine accumulate_rsb

    Integer(Int32) :: i, j, k
    Real(Int64)    :: uc, vc, wc, pc
    Real(Int64)    :: nut_c
    Real(Int64)    :: dudx, dudy, dudz
    Real(Int64)    :: dvdx, dvdy, dvdz
    Real(Int64)    :: dwdx, dwdy, dwdz
    Real(Int64)    :: sij_11, sij_22, sij_33
    Real(Int64)    :: sij_12, sij_13, sij_23
    Real(Int64)    :: inv_dx, inv_dz             ! uniform-grid reciprocals
    Real(Int64)    :: inv_yg_jm, inv_yg_jp       ! per-j y-grid reciprocals (cell-centre spacing)
    Real(Int64)    :: inv_dy_j                    ! per-j face-to-face spacing: 1/(y(j)-y(j-1))

    If ( rsb_active /= 1 ) Return
    If ( istep < rsb_nstart ) Return

    inv_dx = 1d0 / dx
    inv_dz = 1d0 / dz

    ! loop over interior cell centres i,j,k in [2..nxg-1]x[2..nyg-1]x[2..nzg-1]; cell-centre (i,j,k) maps to accumulator index (i-1,j-1,k-1)
    Do k = 2, nzg-1
       Do j = 2, nyg-1
          ! per-j y-grid reciprocals
          inv_yg_jm = 1d0 / ( yg(j)   - yg(j-1) )
          inv_yg_jp = 1d0 / ( yg(j+1) - yg(j  ) )
          inv_dy_j  = 1d0 / ( y (j)   - y (j-1) )  ! face-to-face: correct for V at y-faces
          Do i = 2, nxg-1

             ! IBM: skip solid cells
             If ( ibm_input_mode >= 1 ) Then
                If ( phi(i,j,k) < 0d0 ) Cycle
             End If

             ! ── cell-centre velocities and pressure ──────────────────────
             uc = 0.5d0 * ( U(i-1, j,   k  ) + U(i,   j,   k  ) )
             vc = 0.5d0 * ( V(i,   j-1, k  ) + V(i,   j,   k  ) )
             wc = 0.5d0 * ( W(i,   j,   k-1) + W(i,   j,   k  ) )
             pc = P(i, j, k)

             ! ── SGS eddy viscosity at cell centre ─────────────────────
             nut_c = 0d0
             If ( sgs_model > 0 ) nut_c = nu_t(i, j, k)

             ! ── accumulate 1st moments ───────────────────────────────
             acc_U(i-1,j-1,k-1) = acc_U(i-1,j-1,k-1) + uc
             acc_V(i-1,j-1,k-1) = acc_V(i-1,j-1,k-1) + vc
             acc_W(i-1,j-1,k-1) = acc_W(i-1,j-1,k-1) + wc
             acc_P(i-1,j-1,k-1) = acc_P(i-1,j-1,k-1) + pc

             ! ── accumulate raw 2nd velocity moments (order: 11,22,33,12,13,23) ──
             acc_UiUj(1,i-1,j-1,k-1) = acc_UiUj(1,i-1,j-1,k-1) + uc*uc
             acc_UiUj(2,i-1,j-1,k-1) = acc_UiUj(2,i-1,j-1,k-1) + vc*vc
             acc_UiUj(3,i-1,j-1,k-1) = acc_UiUj(3,i-1,j-1,k-1) + wc*wc
             acc_UiUj(4,i-1,j-1,k-1) = acc_UiUj(4,i-1,j-1,k-1) + uc*vc
             acc_UiUj(5,i-1,j-1,k-1) = acc_UiUj(5,i-1,j-1,k-1) + uc*wc
             acc_UiUj(6,i-1,j-1,k-1) = acc_UiUj(6,i-1,j-1,k-1) + vc*wc

             ! ── pressure-velocity for pressure diffusion ─────────────
             acc_PVel(1,i-1,j-1,k-1) = acc_PVel(1,i-1,j-1,k-1) + pc*uc
             acc_PVel(2,i-1,j-1,k-1) = acc_PVel(2,i-1,j-1,k-1) + pc*vc
             acc_PVel(3,i-1,j-1,k-1) = acc_PVel(3,i-1,j-1,k-1) + pc*wc

             ! ── triple correlations for turbulent diffusion ───────────
             !  y-direction flux (always stored): u_i u_j v
             acc_TijY(1,i-1,j-1,k-1) = acc_TijY(1,i-1,j-1,k-1) + uc*uc*vc
             acc_TijY(2,i-1,j-1,k-1) = acc_TijY(2,i-1,j-1,k-1) + vc*vc*vc
             acc_TijY(3,i-1,j-1,k-1) = acc_TijY(3,i-1,j-1,k-1) + wc*wc*vc
             acc_TijY(4,i-1,j-1,k-1) = acc_TijY(4,i-1,j-1,k-1) + uc*vc*vc
             acc_TijY(5,i-1,j-1,k-1) = acc_TijY(5,i-1,j-1,k-1) + uc*wc*vc
             acc_TijY(6,i-1,j-1,k-1) = acc_TijY(6,i-1,j-1,k-1) + vc*wc*vc

             ! IBM: skip gradient stencils crossing the solid-fluid interface
             If ( ibm_input_mode >= 1 ) Then
                If ( phi(i-1,j,k) < 0d0 .Or. phi(i+1,j,k) < 0d0 .Or. &
                     phi(i,j-1,k) < 0d0 .Or. phi(i,j+1,k) < 0d0 .Or. &
                     phi(i,j,k-1) < 0d0 .Or. phi(i,j,k+1) < 0d0 ) Cycle
             End If

             ! ── velocity gradients ─────────────────────────────────────────────
             dudx = ( U(i,j,k) - U(i-1,j,k) ) * inv_dx
             dudy = 0.5d0*( (U(i,j,k)-U(i,j-1,k))*inv_yg_jm + &
                            (U(i,j+1,k)-U(i,j,k))*inv_yg_jp )
             dudz = 0.5d0*( (U(i,j,k)-U(i,j,k-1)) + &
                            (U(i,j,k+1)-U(i,j,k)) ) * inv_dz

             dvdx = 0.5d0*( (V(i,j,k)-V(i-1,j,k)) + &
                            (V(i+1,j,k)-V(i,j,k)) ) * inv_dx
             dvdy = ( V(i,j,k) - V(i,j-1,k) ) * inv_dy_j  ! V on y-faces: use face spacing
             dvdz = 0.5d0*( (V(i,j,k)-V(i,j,k-1)) + &
                            (V(i,j,k+1)-V(i,j,k)) ) * inv_dz

             dwdx = 0.5d0*( (W(i,j,k)-W(i-1,j,k)) + &
                            (W(i+1,j,k)-W(i,j,k)) ) * inv_dx
             dwdy = 0.5d0*( (W(i,j,k)-W(i,j-1,k))*inv_yg_jm + &
                            (W(i,j+1,k)-W(i,j,k))*inv_yg_jp )
             dwdz = ( W(i,j,k) - W(i,j,k-1) ) * inv_dz

             ! ── strain-rate tensor (symmetric part of grad u) ──────────
             sij_11 = dudx
             sij_22 = dvdy
             sij_33 = dwdz
             sij_12 = 0.5d0 * (dudy + dvdx)
             sij_13 = 0.5d0 * (dudz + dwdx)
             sij_23 = 0.5d0 * (dvdz + dwdy)

             ! ── resolved dissipation products: (dui/dxk)(duj/dxk) ─────
             ! Component 11: sum_k (du/dxk)^2
             acc_GijRes(1,i-1,j-1,k-1) = acc_GijRes(1,i-1,j-1,k-1) + &
                  dudx*dudx + dudy*dudy + dudz*dudz
             ! Component 22: sum_k (dv/dxk)^2
             acc_GijRes(2,i-1,j-1,k-1) = acc_GijRes(2,i-1,j-1,k-1) + &
                  dvdx*dvdx + dvdy*dvdy + dvdz*dvdz
             ! Component 33: sum_k (dw/dxk)^2
             acc_GijRes(3,i-1,j-1,k-1) = acc_GijRes(3,i-1,j-1,k-1) + &
                  dwdx*dwdx + dwdy*dwdy + dwdz*dwdz
             ! Component 12: sum_k (du/dxk)(dv/dxk)
             acc_GijRes(4,i-1,j-1,k-1) = acc_GijRes(4,i-1,j-1,k-1) + &
                  dudx*dvdx + dudy*dvdy + dudz*dvdz
             ! Component 13: sum_k (du/dxk)(dw/dxk)
             acc_GijRes(5,i-1,j-1,k-1) = acc_GijRes(5,i-1,j-1,k-1) + &
                  dudx*dwdx + dudy*dwdy + dudz*dwdz
             ! Component 23: sum_k (dv/dxk)(dw/dxk)
             acc_GijRes(6,i-1,j-1,k-1) = acc_GijRes(6,i-1,j-1,k-1) + &
                  dvdx*dwdx + dvdy*dwdy + dvdz*dwdz

             ! ── SGS dissipation: 2*nu_t*s_ij ──────────────────────────────
             acc_GijSGS(1,i-1,j-1,k-1) = acc_GijSGS(1,i-1,j-1,k-1) + 2d0*nut_c*sij_11
             acc_GijSGS(2,i-1,j-1,k-1) = acc_GijSGS(2,i-1,j-1,k-1) + 2d0*nut_c*sij_22
             acc_GijSGS(3,i-1,j-1,k-1) = acc_GijSGS(3,i-1,j-1,k-1) + 2d0*nut_c*sij_33
             acc_GijSGS(4,i-1,j-1,k-1) = acc_GijSGS(4,i-1,j-1,k-1) + 2d0*nut_c*sij_12
             acc_GijSGS(5,i-1,j-1,k-1) = acc_GijSGS(5,i-1,j-1,k-1) + 2d0*nut_c*sij_13
             acc_GijSGS(6,i-1,j-1,k-1) = acc_GijSGS(6,i-1,j-1,k-1) + 2d0*nut_c*sij_23

             ! ── pressure-strain: p * (du_i/dx_j + du_j/dx_i) ────────
             acc_PiStrain(1,i-1,j-1,k-1) = acc_PiStrain(1,i-1,j-1,k-1) + 2d0*pc*dudx
             acc_PiStrain(2,i-1,j-1,k-1) = acc_PiStrain(2,i-1,j-1,k-1) + 2d0*pc*dvdy
             acc_PiStrain(3,i-1,j-1,k-1) = acc_PiStrain(3,i-1,j-1,k-1) + 2d0*pc*dwdz
             acc_PiStrain(4,i-1,j-1,k-1) = acc_PiStrain(4,i-1,j-1,k-1) + pc*(dudy+dvdx)
             acc_PiStrain(5,i-1,j-1,k-1) = acc_PiStrain(5,i-1,j-1,k-1) + pc*(dudz+dwdx)
             acc_PiStrain(6,i-1,j-1,k-1) = acc_PiStrain(6,i-1,j-1,k-1) + pc*(dvdz+dwdy)

          End Do
       End Do
    End Do

    ! ── acc_TijX: x-direction triple flux (only if x not homogeneous) ────
    If ( .Not. hom_x ) Then
       Do k = 2, nzg-1
          Do j = 2, nyg-1
             Do i = 2, nxg-1
                If ( ibm_input_mode >= 1 ) Then
                   If ( phi(i,j,k) < 0d0 ) Cycle
                End If
                uc = 0.5d0 * ( U(i-1,j,k) + U(i,j,k) )
                vc = 0.5d0 * ( V(i,j-1,k) + V(i,j,k) )
                wc = 0.5d0 * ( W(i,j,k-1) + W(i,j,k) )
                acc_TijX(1,i-1,j-1,k-1) = acc_TijX(1,i-1,j-1,k-1) + uc*uc*uc
                acc_TijX(2,i-1,j-1,k-1) = acc_TijX(2,i-1,j-1,k-1) + vc*vc*uc
                acc_TijX(3,i-1,j-1,k-1) = acc_TijX(3,i-1,j-1,k-1) + wc*wc*uc
                acc_TijX(4,i-1,j-1,k-1) = acc_TijX(4,i-1,j-1,k-1) + uc*vc*uc
                acc_TijX(5,i-1,j-1,k-1) = acc_TijX(5,i-1,j-1,k-1) + uc*wc*uc
                acc_TijX(6,i-1,j-1,k-1) = acc_TijX(6,i-1,j-1,k-1) + vc*wc*uc
             End Do
          End Do
       End Do
    End If

    ! ── acc_TijZ: z-direction triple flux (only if z not homogeneous) ────
    If ( .Not. hom_z ) Then
       Do k = 2, nzg-1
          Do j = 2, nyg-1
             Do i = 2, nxg-1
                If ( ibm_input_mode >= 1 ) Then
                   If ( phi(i,j,k) < 0d0 ) Cycle
                End If
                uc = 0.5d0 * ( U(i-1,j,k) + U(i,j,k) )
                vc = 0.5d0 * ( V(i,j-1,k) + V(i,j,k) )
                wc = 0.5d0 * ( W(i,j,k-1) + W(i,j,k) )
                acc_TijZ(1,i-1,j-1,k-1) = acc_TijZ(1,i-1,j-1,k-1) + uc*uc*wc
                acc_TijZ(2,i-1,j-1,k-1) = acc_TijZ(2,i-1,j-1,k-1) + vc*vc*wc
                acc_TijZ(3,i-1,j-1,k-1) = acc_TijZ(3,i-1,j-1,k-1) + wc*wc*wc
                acc_TijZ(4,i-1,j-1,k-1) = acc_TijZ(4,i-1,j-1,k-1) + uc*vc*wc
                acc_TijZ(5,i-1,j-1,k-1) = acc_TijZ(5,i-1,j-1,k-1) + uc*wc*wc
                acc_TijZ(6,i-1,j-1,k-1) = acc_TijZ(6,i-1,j-1,k-1) + vc*wc*wc
             End Do
          End Do
       End Do
    End If

    n_accum = n_accum + 1

  End Subroutine accumulate_rsb


  !> output_rsb: called every rsb_freq steps; MPI-reduces accumulators, converts to central moments, computes all budget terms, appends one slice per .bin file, updates meta, resets accumulators
  Subroutine output_rsb

    ! ── global (rank-0) arrays after MPI-reduce and hom. averaging ───────
    Real(Int64), Allocatable :: g_U(:,:,:),  g_V(:,:,:),  g_W(:,:,:)
    Real(Int64), Allocatable :: g_P(:,:,:)
    Real(Int64), Allocatable :: g_UiUj(:,:,:,:)
    Real(Int64), Allocatable :: g_GijRes(:,:,:,:), g_GijSGS(:,:,:,:)
    Real(Int64), Allocatable :: g_PiStrain(:,:,:,:), g_PVel(:,:,:,:)
    Real(Int64), Allocatable :: g_TijY(:,:,:,:)
    Real(Int64), Allocatable :: g_TijX(:,:,:,:), g_TijZ(:,:,:,:)

    ! ── budget term arrays (output grid: out_nx × out_ny × out_nz) ───────
    Real(Int64), Allocatable :: Rij(:,:,:,:)     ! Reynolds stresses (central moments)
    Real(Int64), Allocatable :: Pij(:,:,:,:)     ! production
    Real(Int64), Allocatable :: epsRes(:,:,:,:)  ! resolved dissipation
    Real(Int64), Allocatable :: epsSGS(:,:,:,:)  ! SGS dissipation
    Real(Int64), Allocatable :: PiStrain(:,:,:,:)! pressure-strain
    Real(Int64), Allocatable :: DTij(:,:,:,:)    ! turbulent diffusion
    Real(Int64), Allocatable :: Dnuij(:,:,:,:)   ! viscous diffusion
    Real(Int64), Allocatable :: PhiPij(:,:,:,:)  ! pressure diffusion
    Real(Int64), Allocatable :: Umean(:,:,:,:)   ! mean velocities (3 components)
    Real(Int64), Allocatable :: Resid(:,:,:,:)   ! budget residual

    Integer(Int32) :: c
    Real(Int64)    :: dn   ! = 1 / n_accum (for averaging)

    If ( rsb_active /= 1 ) Return
    If ( n_accum == 0 ) Return  ! nothing accumulated yet

    dn = 1d0 / Real(n_accum, Int64)

    ! ── MPI-reduce accumulators to rank 0 ─────────────────────────────
    Allocate( g_U      (  nxm_global, nym_global, nzm_global) )
    Allocate( g_V      (  nxm_global, nym_global, nzm_global) )
    Allocate( g_W      (  nxm_global, nym_global, nzm_global) )
    Allocate( g_P      (  nxm_global, nym_global, nzm_global) )
    Allocate( g_UiUj   (6,nxm_global, nym_global, nzm_global) )
    Allocate( g_GijRes (6,nxm_global, nym_global, nzm_global) )
    Allocate( g_GijSGS (6,nxm_global, nym_global, nzm_global) )
    Allocate( g_PiStrain(6,nxm_global,nym_global, nzm_global) )
    Allocate( g_PVel   (3,nxm_global, nym_global, nzm_global) )
    Allocate( g_TijY   (6,nxm_global, nym_global, nzm_global) )
    Allocate( g_TijX(6,nxm_global,nym_global,nzm_global) ); g_TijX = 0d0
    Allocate( g_TijZ(6,nxm_global,nym_global,nzm_global) ); g_TijZ = 0d0

    Call reduce_to_rank0(acc_U,      g_U,       nxm*nym_global*nzm)
    Call reduce_to_rank0(acc_V,      g_V,       nxm*nym_global*nzm)
    Call reduce_to_rank0(acc_W,      g_W,       nxm*nym_global*nzm)
    Call reduce_to_rank0(acc_P,      g_P,       nxm*nym_global*nzm)
    Call reduce_to_rank0_4d(acc_UiUj,    g_UiUj,    6,nxm,nym_global,nzm)
    Call reduce_to_rank0_4d(acc_GijRes,  g_GijRes,  6,nxm,nym_global,nzm)
    Call reduce_to_rank0_4d(acc_GijSGS,  g_GijSGS,  6,nxm,nym_global,nzm)
    Call reduce_to_rank0_4d(acc_PiStrain,g_PiStrain,6,nxm,nym_global,nzm)
    Call reduce_to_rank0_4d(acc_PVel,    g_PVel,    3,nxm,nym_global,nzm)
    Call reduce_to_rank0_4d(acc_TijY,    g_TijY,    6,nxm,nym_global,nzm)
    If ( .Not. hom_x ) Call reduce_to_rank0_4d(acc_TijX,g_TijX,6,nxm,nym_global,nzm)
    If ( .Not. hom_z ) Call reduce_to_rank0_4d(acc_TijZ,g_TijZ,6,nxm,nym_global,nzm)

    ! ── rank 0: compute budget and write ─────────────────────────────────
    If ( myid == 0 ) Then

       ! scale by 1/n_accum
       g_U       = g_U       * dn;  g_V  = g_V  * dn;  g_W = g_W * dn
       g_P       = g_P       * dn
       g_UiUj    = g_UiUj    * dn
       g_GijRes  = g_GijRes  * dn
       g_GijSGS  = g_GijSGS  * dn
       g_PiStrain= g_PiStrain* dn
       g_PVel    = g_PVel    * dn
       g_TijY    = g_TijY    * dn
       If ( .Not. hom_x ) g_TijX = g_TijX * dn
       If ( .Not. hom_z ) g_TijZ = g_TijZ * dn

       ! allocate output budget arrays on the reduced (out_nx x out_ny x out_nz) grid
       Allocate( Rij     (N_COMP, out_nx, out_ny, out_nz) )
       Allocate( Pij     (N_COMP, out_nx, out_ny, out_nz) )
       Allocate( epsRes  (N_COMP, out_nx, out_ny, out_nz) )
       Allocate( epsSGS  (N_COMP, out_nx, out_ny, out_nz) )
       Allocate( PiStrain(N_COMP, out_nx, out_ny, out_nz) )
       Allocate( DTij    (N_COMP, out_nx, out_ny, out_nz) )
       Allocate( Dnuij   (N_COMP, out_nx, out_ny, out_nz) )
       Allocate( PhiPij  (N_COMP, out_nx, out_ny, out_nz) )
       Allocate( Umean   (3,      out_nx, out_ny, out_nz) )
       Allocate( Resid   (N_COMP, out_nx, out_ny, out_nz) )

       ! compute central moments and all budget terms on the output grid
       Call compute_budget_terms( &
            g_U, g_V, g_W, g_P, g_UiUj, g_GijRes, g_GijSGS, &
            g_PiStrain, g_PVel, g_TijY, g_TijX, g_TijZ,      &
            Rij, Pij, epsRes, epsSGS, PiStrain,               &
            DTij, Dnuij, PhiPij, Umean, Resid)

       nsamples = nsamples + 1

       ! write one slice of each budget group to its .bin file
       Call write_rsb_slice('Umean',    Umean,    3)
       Call write_rsb_slice('Rij',      Rij,      N_COMP)
       Call write_rsb_slice('Pij',      Pij,      N_COMP)
       Call write_rsb_slice('epsRes',   epsRes,   N_COMP)
       Call write_rsb_slice('epsSGS',   epsSGS,   N_COMP)
       Call write_rsb_slice('PiStrain', PiStrain, N_COMP)
       Call write_rsb_slice('DTij',     DTij,     N_COMP)
       Call write_rsb_slice('Dnuij',    Dnuij,    N_COMP)
       Call write_rsb_slice('PhiPij',   PhiPij,   N_COMP)
       Call write_rsb_slice('Resid',    Resid,    N_COMP)

       ! update meta file
       Call write_meta

       Write(*,'(A,I8,A,I8,A)') ' RSB: wrote sample ', nsamples, &
            ' at step ', istep, '.'

       Deallocate(Rij, Pij, epsRes, epsSGS, PiStrain)
       Deallocate(DTij, Dnuij, PhiPij, Umean, Resid)
    End If

    ! ── deallocate global buffers ─────────────────────────────────────────
    Deallocate(g_U, g_V, g_W, g_P, g_UiUj, g_GijRes, g_GijSGS)
    Deallocate(g_PiStrain, g_PVel, g_TijY)
    Deallocate(g_TijX, g_TijZ)

    ! ── reset accumulators for next window ───────────────────────────────
    Call zero_accumulators
    n_accum = 0

  End Subroutine output_rsb


  !  finalize_rsb
  Subroutine finalize_rsb

    If ( rsb_active /= 1 ) Return

    Deallocate(acc_U, acc_V, acc_W, acc_P)
    Deallocate(acc_UiUj, acc_GijRes, acc_GijSGS)
    Deallocate(acc_PiStrain, acc_PVel, acc_TijY)
    If ( Allocated(acc_TijX) ) Deallocate(acc_TijX)
    If ( Allocated(acc_TijZ) ) Deallocate(acc_TijZ)

  End Subroutine finalize_rsb


  !  Private helper: zero_accumulators
  Subroutine zero_accumulators
    acc_U        = 0d0;  acc_V = 0d0;  acc_W = 0d0;  acc_P = 0d0
    acc_UiUj     = 0d0
    acc_GijRes   = 0d0
    acc_GijSGS   = 0d0
    acc_PiStrain = 0d0
    acc_PVel     = 0d0
    acc_TijY     = 0d0
    If ( Allocated(acc_TijX) ) acc_TijX = 0d0
    If ( Allocated(acc_TijZ) ) acc_TijZ = 0d0
  End Subroutine zero_accumulators


  !> ensure_rsb_dir (rank 0 only): extracts the directory component of rsb_fileout and creates it (with parents) if missing
  Subroutine ensure_rsb_dir

    Character(200) :: dir
    Integer(Int32) :: slash_pos
    Logical        :: dir_exists

    slash_pos = Index(Trim(rsb_fileout), '/', Back=.True.)
    If ( slash_pos < 1 ) Return   ! no directory component; output goes to cwd

    dir = rsb_fileout(1:slash_pos-1)
    Inquire(file=Trim(dir)//'/.', exist=dir_exists)
    If ( .Not. dir_exists ) Then
       Call execute_command_line('mkdir -p ' // Trim(dir), wait=.True.)
       Write(*,'(A,A)') '   RSB: created output directory ', Trim(dir)
    End If

  End Subroutine ensure_rsb_dir


  !  Gather a local z-slab array into the full global array on rank 0
  !  using MPI_Gatherv along the z decomposition.
  Subroutine reduce_to_rank0(local_arr, global_arr, local_sz)

    Integer(Int32), Intent(In)  :: local_sz
    Real(Int64),    Intent(In)  :: local_arr(local_sz)
    Real(Int64),    Intent(Out) :: global_arr(*)

    Integer(Int32) :: sendcount
    Integer(Int32), Allocatable :: recvcounts(:), displs(:)
    Integer(Int32) :: iproc, off

    sendcount = local_sz

    ! rank 0 builds recvcounts and displacements for Gatherv
    If ( myid == 0 ) Then
       Allocate(recvcounts(nprocs), displs(nprocs))
       off = 0
       Do iproc = 0, nprocs-1
          ! each rank owns (kg2_global(iproc)-kg1_global(iproc)+1-2) z-centre cells
          ! (subtract 2 ghost planes; accumulator arrays have size nzm, not nzg)
          recvcounts(iproc+1) = nxm_global * nym_global * &
               (kg2_global(iproc) - kg1_global(iproc) + 1 - 2)
          displs(iproc+1) = off
          off = off + recvcounts(iproc+1)
       End Do
       Call MPI_Gatherv(local_arr, sendcount, MPI_REAL8, &
                        global_arr(1), recvcounts, displs, MPI_REAL8, &
                        0, MPI_COMM_WORLD, ierr)
       Deallocate(recvcounts, displs)
    Else
       Call MPI_Gatherv(local_arr, sendcount, MPI_REAL8, &
                        global_arr(1), sendcount, 0, MPI_REAL8, &
                        0, MPI_COMM_WORLD, ierr)
    End If

  End Subroutine reduce_to_rank0


  !  reduce_to_rank0_4d — wraps reduce_to_rank0 for 4D (nc,nx,ny,nz) arrays
  Subroutine reduce_to_rank0_4d(local_arr, global_arr, nc, lnx, lny, lnz)

    Integer(Int32), Intent(In)  :: nc, lnx, lny, lnz
    Real(Int64),    Intent(In)  :: local_arr (nc, lnx, lny, lnz)
    Real(Int64),    Intent(Out) :: global_arr(nc, nxm_global, nym_global, nzm_global)

    Integer(Int32) :: c

    Do c = 1, nc
       Call reduce_to_rank0(local_arr(c,:,:,:), global_arr(c,:,:,:), lnx*lny*lnz)
    End Do

  End Subroutine reduce_to_rank0_4d


  !> compute_budget_terms: all arithmetic on rank-0 only over the full global grid; budget terms are homogeneously averaged after central moments are computed
  Subroutine compute_budget_terms( &
       gU, gV, gW, gP, gUiUj, gGijRes, gGijSGS,     &
       gPiStrain, gPVel, gTijY, gTijX, gTijZ,         &
       Rij, Pij, epsResOut, epsSGSOut, PiStrainOut,    &
       DTijOut, DnuijOut, PhiPijOut, UmeanOut, ResidOut)

    Real(Int64), Intent(In) :: gU(nxm_global,nym_global,nzm_global)
    Real(Int64), Intent(In) :: gV(nxm_global,nym_global,nzm_global)
    Real(Int64), Intent(In) :: gW(nxm_global,nym_global,nzm_global)
    Real(Int64), Intent(In) :: gP(nxm_global,nym_global,nzm_global)
    Real(Int64), Intent(In) :: gUiUj   (N_COMP,nxm_global,nym_global,nzm_global)
    Real(Int64), Intent(In) :: gGijRes (N_COMP,nxm_global,nym_global,nzm_global)
    Real(Int64), Intent(In) :: gGijSGS (N_COMP,nxm_global,nym_global,nzm_global)
    Real(Int64), Intent(In) :: gPiStrain(N_COMP,nxm_global,nym_global,nzm_global)
    Real(Int64), Intent(In) :: gPVel   (3,nxm_global,nym_global,nzm_global)
    Real(Int64), Intent(In) :: gTijY   (N_COMP,nxm_global,nym_global,nzm_global)
    Real(Int64), Intent(In) :: gTijX   (N_COMP,nxm_global,nym_global,nzm_global)
    Real(Int64), Intent(In) :: gTijZ   (N_COMP,nxm_global,nym_global,nzm_global)
    Real(Int64), Intent(Out):: Rij     (N_COMP,out_nx,out_ny,out_nz)
    Real(Int64), Intent(Out):: Pij     (N_COMP,out_nx,out_ny,out_nz)
    Real(Int64), Intent(Out):: epsResOut(N_COMP,out_nx,out_ny,out_nz)
    Real(Int64), Intent(Out):: epsSGSOut(N_COMP,out_nx,out_ny,out_nz)
    Real(Int64), Intent(Out):: PiStrainOut(N_COMP,out_nx,out_ny,out_nz)
    Real(Int64), Intent(Out):: DTijOut (N_COMP,out_nx,out_ny,out_nz)
    Real(Int64), Intent(Out):: DnuijOut(N_COMP,out_nx,out_ny,out_nz)
    Real(Int64), Intent(Out):: PhiPijOut(N_COMP,out_nx,out_ny,out_nz)
    Real(Int64), Intent(Out):: UmeanOut(3,out_nx,out_ny,out_nz)
    Real(Int64), Intent(Out):: ResidOut(N_COMP,out_nx,out_ny,out_nz)

    ! ── working arrays on the full global grid (computed before averaging) ──
    Real(Int64), Allocatable :: Rij_g(:,:,:,:)
    Real(Int64), Allocatable :: Pij_g(:,:,:,:)
    Real(Int64), Allocatable :: epsRes_g(:,:,:,:)
    Real(Int64), Allocatable :: epsSGS_g(:,:,:,:)
    Real(Int64), Allocatable :: PiStrain_g(:,:,:,:)
    Real(Int64), Allocatable :: DTij_g(:,:,:,:)
    Real(Int64), Allocatable :: Dnuij_g(:,:,:,:)
    Real(Int64), Allocatable :: PhiP_g(:,:,:,:)
    Real(Int64), Allocatable :: Umean_g(:,:,:,:)
    ! Central-moment triple-flux arrays (corrected from raw accumulator moments)
    Real(Int64), Allocatable :: cTijY(:,:,:,:), cTijX(:,:,:,:), cTijZ(:,:,:,:)
    ! Central-moment pressure-velocity fields: <p'u_i'> = <p*u_i> - <p><u_i>
    Real(Int64), Allocatable :: pu_c(:,:,:), pv_c(:,:,:), pw_c(:,:,:)

    Integer(Int32) :: i, j, k, c
    Real(Int64)    :: Ub, Vb, Wb          ! local mean velocities
    Real(Int64)    :: Rij_val(N_COMP)     ! local central-moment vector
    Real(Int64)    :: dRij_dy, dRij_dx, dRij_dz   ! spatial derivatives for Dnu
    Real(Int64)    :: dTy_dy, dTx_dx, dTz_dz       ! triple-flux divergences
    Real(Int64)    :: dPux_dx, dPuy_dy, dPuz_dz    ! pressure-velocity flux divs
    Real(Int64)    :: dPvx_dx, dPvy_dy, dPvz_dz
    Real(Int64)    :: dPwx_dx, dPwy_dy, dPwz_dz
    Real(Int64)    :: dUdx, dUdy, dUdz              ! mean-velocity gradients
    Real(Int64)    :: dVdx, dVdy, dVdz
    Real(Int64)    :: dWdx, dWdy, dWdz

    Allocate( Rij_g     (N_COMP,nxm_global,nym_global,nzm_global) )
    Allocate( Pij_g     (N_COMP,nxm_global,nym_global,nzm_global) )
    Allocate( epsRes_g  (N_COMP,nxm_global,nym_global,nzm_global) )
    Allocate( epsSGS_g  (N_COMP,nxm_global,nym_global,nzm_global) )
    Allocate( PiStrain_g(N_COMP,nxm_global,nym_global,nzm_global) )
    Allocate( DTij_g    (N_COMP,nxm_global,nym_global,nzm_global) )
    Allocate( Dnuij_g   (N_COMP,nxm_global,nym_global,nzm_global) )
    Allocate( PhiP_g    (N_COMP,nxm_global,nym_global,nzm_global) )
    Allocate( Umean_g   (3,     nxm_global,nym_global,nzm_global) )
    Allocate( cTijY     (N_COMP,nxm_global,nym_global,nzm_global) )
    Allocate( cTijX     (N_COMP,nxm_global,nym_global,nzm_global) )
    Allocate( cTijZ     (N_COMP,nxm_global,nym_global,nzm_global) )
    Allocate( pu_c      (nxm_global,nym_global,nzm_global) )
    Allocate( pv_c      (nxm_global,nym_global,nzm_global) )
    Allocate( pw_c      (nxm_global,nym_global,nzm_global) )

    ! strides on the global cell-centre grid

    Umean_g(1,:,:,:) = gU
    Umean_g(2,:,:,:) = gV
    Umean_g(3,:,:,:) = gW

    Do k = 1, nzm_global
       Do j = 1, nym_global
          Do i = 1, nxm_global

             Ub = gU(i,j,k);  Vb = gV(i,j,k);  Wb = gW(i,j,k)

             ! ── Reynolds stress: R_ij = <U_iU_j> - <U_i><U_j> ────────
             Rij_g(1,i,j,k) = gUiUj(1,i,j,k) - Ub*Ub
             Rij_g(2,i,j,k) = gUiUj(2,i,j,k) - Vb*Vb
             Rij_g(3,i,j,k) = gUiUj(3,i,j,k) - Wb*Wb
             Rij_g(4,i,j,k) = gUiUj(4,i,j,k) - Ub*Vb
             Rij_g(5,i,j,k) = gUiUj(5,i,j,k) - Ub*Wb
             Rij_g(6,i,j,k) = gUiUj(6,i,j,k) - Vb*Wb

             ! resolved dissipation eps_ij = 2*nu*<(dui/dxk)(duj/dxk)>; mean-gradient term (d<ui>/dxk)(d<uj>/dxk) subtracted below via centred differences
             epsRes_g(:,i,j,k) = 2d0 * nu * gGijRes(:,i,j,k)
             ! (mean-gradient correction is done in the production loop below
             !  using the same mean-gradient estimates)

             ! pressure-strain (raw <p(dui/dxj+duj/dxi)>); mean part <p><dui/dxj>+<p><duj/dxi> ~=0 for periodic, so raw value equals the fluctuation part to leading order
             PiStrain_g(:,i,j,k) = gPiStrain(:,i,j,k)

             ! ── SGS dissipation ────────────────────────────────────────
             epsSGS_g(:,i,j,k) = gGijSGS(:,i,j,k)

          End Do
       End Do
    End Do

    ! ── Pre-compute central-moment pressure-velocity fields for PhiP ──────
    pu_c = gPVel(1,:,:,:) - gP * gU
    pv_c = gPVel(2,:,:,:) - gP * gV
    pw_c = gPVel(3,:,:,:) - gP * gW

    ! convert raw triple-flux correlations to central moments
    cTijX = 0d0
    cTijZ = 0d0
    Do k = 1, nzm_global
       Do j = 1, nym_global
          Do i = 1, nxm_global
             Ub = gU(i,j,k); Vb = gV(i,j,k); Wb = gW(i,j,k)
             ! ── y-flux: <u_i u_j v> → <u_i' u_j' v'> ─────────────────
             ! c=1 (u,u): -R11*<v> - 2R12*<u> - <u>^2<v>
             cTijY(1,i,j,k) = gTijY(1,i,j,k) - Rij_g(1,i,j,k)*Vb - 2d0*Rij_g(4,i,j,k)*Ub - Ub**2*Vb
             ! c=2 (v,v): -3R22*<v> - <v>^3
             cTijY(2,i,j,k) = gTijY(2,i,j,k) - 3d0*Rij_g(2,i,j,k)*Vb - Vb**3
             ! c=3 (w,w): -R33*<v> - 2R23*<w> - <w>^2<v>
             cTijY(3,i,j,k) = gTijY(3,i,j,k) - Rij_g(3,i,j,k)*Vb - 2d0*Rij_g(6,i,j,k)*Wb - Wb**2*Vb
             ! c=4 (u,v): -2R12*<v> - R22*<u> - <u><v>^2
             cTijY(4,i,j,k) = gTijY(4,i,j,k) - 2d0*Rij_g(4,i,j,k)*Vb - Rij_g(2,i,j,k)*Ub - Ub*Vb**2
             ! c=5 (u,w): -R13*<v> - R12*<w> - R23*<u> - <u><w><v>
             cTijY(5,i,j,k) = gTijY(5,i,j,k) - Rij_g(5,i,j,k)*Vb - Rij_g(4,i,j,k)*Wb &
                              - Rij_g(6,i,j,k)*Ub - Ub*Wb*Vb
             ! c=6 (v,w): -2R23*<v> - R22*<w> - <v>^2<w>
             cTijY(6,i,j,k) = gTijY(6,i,j,k) - 2d0*Rij_g(6,i,j,k)*Vb - Rij_g(2,i,j,k)*Wb - Vb**2*Wb

             If ( .Not. hom_x ) Then
                ! ── x-flux: <u_i u_j u> → <u_i' u_j' u'> ────────────────
                ! c=1 (u,u): -3R11*<u> - <u>^3
                cTijX(1,i,j,k) = gTijX(1,i,j,k) - 3d0*Rij_g(1,i,j,k)*Ub - Ub**3
                ! c=2 (v,v): -R22*<u> - 2R12*<v> - <v>^2<u>
                cTijX(2,i,j,k) = gTijX(2,i,j,k) - Rij_g(2,i,j,k)*Ub - 2d0*Rij_g(4,i,j,k)*Vb - Vb**2*Ub
                ! c=3 (w,w): -R33*<u> - 2R13*<w> - <w>^2<u>
                cTijX(3,i,j,k) = gTijX(3,i,j,k) - Rij_g(3,i,j,k)*Ub - 2d0*Rij_g(5,i,j,k)*Wb - Wb**2*Ub
                ! c=4 (u,v): -2R12*<u> - R11*<v> - <u>^2<v>
                cTijX(4,i,j,k) = gTijX(4,i,j,k) - 2d0*Rij_g(4,i,j,k)*Ub - Rij_g(1,i,j,k)*Vb - Ub**2*Vb
                ! c=5 (u,w): -2R13*<u> - R11*<w> - <u>^2<w>
                cTijX(5,i,j,k) = gTijX(5,i,j,k) - 2d0*Rij_g(5,i,j,k)*Ub - Rij_g(1,i,j,k)*Wb - Ub**2*Wb
                ! c=6 (v,w): -R23*<u> - R12*<w> - R13*<v> - <u><v><w>
                cTijX(6,i,j,k) = gTijX(6,i,j,k) - Rij_g(6,i,j,k)*Ub - Rij_g(4,i,j,k)*Wb &
                                 - Rij_g(5,i,j,k)*Vb - Ub*Vb*Wb
             End If

             If ( .Not. hom_z ) Then
                ! ── z-flux: <u_i u_j w> → <u_i' u_j' w'> ────────────────
                ! c=1 (u,u): -R11*<w> - 2R13*<u> - <u>^2<w>
                cTijZ(1,i,j,k) = gTijZ(1,i,j,k) - Rij_g(1,i,j,k)*Wb - 2d0*Rij_g(5,i,j,k)*Ub - Ub**2*Wb
                ! c=2 (v,v): -R22*<w> - 2R23*<v> - <v>^2<w>
                cTijZ(2,i,j,k) = gTijZ(2,i,j,k) - Rij_g(2,i,j,k)*Wb - 2d0*Rij_g(6,i,j,k)*Vb - Vb**2*Wb
                ! c=3 (w,w): -3R33*<w> - <w>^3
                cTijZ(3,i,j,k) = gTijZ(3,i,j,k) - 3d0*Rij_g(3,i,j,k)*Wb - Wb**3
                ! c=4 (u,v): -R12*<w> - R13*<v> - R23*<u> - <u><v><w>
                cTijZ(4,i,j,k) = gTijZ(4,i,j,k) - Rij_g(4,i,j,k)*Wb - Rij_g(5,i,j,k)*Vb &
                                 - Rij_g(6,i,j,k)*Ub - Ub*Vb*Wb
                ! c=5 (u,w): -2R13*<w> - R33*<u> - <u><w>^2
                cTijZ(5,i,j,k) = gTijZ(5,i,j,k) - 2d0*Rij_g(5,i,j,k)*Wb - Rij_g(3,i,j,k)*Ub - Ub*Wb**2
                ! c=6 (v,w): -2R23*<w> - R33*<v> - <v><w>^2
                cTijZ(6,i,j,k) = gTijZ(6,i,j,k) - 2d0*Rij_g(6,i,j,k)*Wb - Rij_g(3,i,j,k)*Vb - Vb*Wb**2
             End If
          End Do
       End Do
    End Do

    ! mean-velocity gradients + production P_ij=-R_ik dUj/dxk-R_jk dUi/dxk (Pope 7.176) and viscous diffusion D^nu_ij=nu d^2 R_ij/dxk dxk
    Pij_g   = 0d0
    Dnuij_g = 0d0

    Do k = 1, nzm_global
       Do j = 1, nym_global
          Do i = 1, nxm_global

             ! -- mean velocity gradients (centred where possible, 1-sided at BCs) --
             If ( i > 1 .And. i < nxm_global ) Then
                dUdx = (gU(i+1,j,k) - gU(i-1,j,k)) / (xm_global(i+1) - xm_global(i-1))
                dVdx = (gV(i+1,j,k) - gV(i-1,j,k)) / (xm_global(i+1) - xm_global(i-1))
                dWdx = (gW(i+1,j,k) - gW(i-1,j,k)) / (xm_global(i+1) - xm_global(i-1))
             Else If ( i == 1 ) Then
                dUdx = (gU(2,j,k) - gU(1,j,k)) / (xm_global(2) - xm_global(1))
                dVdx = (gV(2,j,k) - gV(1,j,k)) / (xm_global(2) - xm_global(1))
                dWdx = (gW(2,j,k) - gW(1,j,k)) / (xm_global(2) - xm_global(1))
             Else
                dUdx = (gU(nxm_global,j,k)-gU(nxm_global-1,j,k)) / &
                       (xm_global(nxm_global)-xm_global(nxm_global-1))
                dVdx = (gV(nxm_global,j,k)-gV(nxm_global-1,j,k)) / &
                       (xm_global(nxm_global)-xm_global(nxm_global-1))
                dWdx = (gW(nxm_global,j,k)-gW(nxm_global-1,j,k)) / &
                       (xm_global(nxm_global)-xm_global(nxm_global-1))
             End If

             If ( j > 1 .And. j < nym_global ) Then
                dUdy = (gU(i,j+1,k) - gU(i,j-1,k)) / (ym_global(j+1) - ym_global(j-1))
                dVdy = (gV(i,j+1,k) - gV(i,j-1,k)) / (ym_global(j+1) - ym_global(j-1))
                dWdy = (gW(i,j+1,k) - gW(i,j-1,k)) / (ym_global(j+1) - ym_global(j-1))
             Else If ( j == 1 ) Then
                dUdy = (gU(i,2,k) - gU(i,1,k)) / (ym_global(2) - ym_global(1))
                dVdy = (gV(i,2,k) - gV(i,1,k)) / (ym_global(2) - ym_global(1))
                dWdy = (gW(i,2,k) - gW(i,1,k)) / (ym_global(2) - ym_global(1))
             Else
                dUdy = (gU(i,nym_global,k)-gU(i,nym_global-1,k)) / &
                       (ym_global(nym_global)-ym_global(nym_global-1))
                dVdy = (gV(i,nym_global,k)-gV(i,nym_global-1,k)) / &
                       (ym_global(nym_global)-ym_global(nym_global-1))
                dWdy = (gW(i,nym_global,k)-gW(i,nym_global-1,k)) / &
                       (ym_global(nym_global)-ym_global(nym_global-1))
             End If

             If ( k > 1 .And. k < nzm_global ) Then
                dUdz = (gU(i,j,k+1) - gU(i,j,k-1)) / (zm_global(k+1) - zm_global(k-1))
                dVdz = (gV(i,j,k+1) - gV(i,j,k-1)) / (zm_global(k+1) - zm_global(k-1))
                dWdz = (gW(i,j,k+1) - gW(i,j,k-1)) / (zm_global(k+1) - zm_global(k-1))
             Else If ( k == 1 ) Then
                dUdz = (gU(i,j,2) - gU(i,j,1)) / (zm_global(2) - zm_global(1))
                dVdz = (gV(i,j,2) - gV(i,j,1)) / (zm_global(2) - zm_global(1))
                dWdz = (gW(i,j,2) - gW(i,j,1)) / (zm_global(2) - zm_global(1))
             Else
                dUdz = (gU(i,j,nzm_global)-gU(i,j,nzm_global-1)) / &
                       (zm_global(nzm_global)-zm_global(nzm_global-1))
                dVdz = (gV(i,j,nzm_global)-gV(i,j,nzm_global-1)) / &
                       (zm_global(nzm_global)-zm_global(nzm_global-1))
                dWdz = (gW(i,j,nzm_global)-gW(i,j,nzm_global-1)) / &
                       (zm_global(nzm_global)-zm_global(nzm_global-1))
             End If

             ! production P_ij=-R_ik dUj/dxk-R_jk dUi/dxk; component order 1=11,2=22,3=33,4=12,5=13,6=23 (R21=R12, R31=R13, R32=R23 by symmetry)
             Rij_val = Rij_g(:,i,j,k)

             !  P_11 = -2(R11 dU/dx + R12 dU/dy + R13 dU/dz)
             Pij_g(1,i,j,k) = -2d0*( Rij_val(1)*dUdx + Rij_val(4)*dUdy + Rij_val(5)*dUdz )
             !  P_22 = -2(R21 dV/dx + R22 dV/dy + R23 dV/dz)
             Pij_g(2,i,j,k) = -2d0*( Rij_val(4)*dVdx + Rij_val(2)*dVdy + Rij_val(6)*dVdz )
             !  P_33 = -2(R31 dW/dx + R32 dW/dy + R33 dW/dz)
             Pij_g(3,i,j,k) = -2d0*( Rij_val(5)*dWdx + Rij_val(6)*dWdy + Rij_val(3)*dWdz )
             !  P_12 = -(R11 dV/dx + R12 dV/dy + R13 dV/dz)
             !        -(R21 dU/dx + R22 dU/dy + R23 dU/dz)
             Pij_g(4,i,j,k) = -( Rij_val(1)*dVdx + Rij_val(4)*dVdy + Rij_val(5)*dVdz &
                                + Rij_val(4)*dUdx + Rij_val(2)*dUdy + Rij_val(6)*dUdz )
             !  P_13 = -(R11 dW/dx + R12 dW/dy + R13 dW/dz)
             !        -(R31 dU/dx + R32 dU/dy + R33 dU/dz)
             Pij_g(5,i,j,k) = -( Rij_val(1)*dWdx + Rij_val(4)*dWdy + Rij_val(5)*dWdz &
                                + Rij_val(5)*dUdx + Rij_val(6)*dUdy + Rij_val(3)*dUdz )
             !  P_23 = -(R21 dW/dx + R22 dW/dy + R23 dW/dz)
             !        -(R31 dV/dx + R32 dV/dy + R33 dV/dz)
             Pij_g(6,i,j,k) = -( Rij_val(4)*dWdx + Rij_val(2)*dWdy + Rij_val(6)*dWdz &
                                + Rij_val(5)*dVdx + Rij_val(6)*dVdy + Rij_val(3)*dVdz )

             ! viscous diffusion nu*Laplacian(R_ij): three-point second derivative using cell-centre spacing, one-sided first differences at boundary nodes
             Do c = 1, N_COMP
                ! d^2 R_ij / dy^2  (always computed)
                If ( j > 1 .And. j < nym_global ) Then
                   dRij_dy = ( (Rij_g(c,i,j+1,k)-Rij_g(c,i,j,k))/(ym_global(j+1)-ym_global(j)) &
                              -(Rij_g(c,i,j,k)-Rij_g(c,i,j-1,k))/(ym_global(j)-ym_global(j-1)) ) &
                             / ( 0.5d0*(ym_global(j+1)-ym_global(j-1)) )
                Else
                   dRij_dy = 0d0   ! wall: Rij=0 by no-slip; d2/dy2=0 at boundary node
                End If

                ! d^2 R_ij / dx^2  (only if x is not homogeneous)
                dRij_dx = 0d0
                If ( .Not. hom_x .And. i > 1 .And. i < nxm_global ) Then
                   dRij_dx = ( (Rij_g(c,i+1,j,k)-Rij_g(c,i,j,k))/(xm_global(i+1)-xm_global(i)) &
                              -(Rij_g(c,i,j,k)-Rij_g(c,i-1,j,k))/(xm_global(i)-xm_global(i-1)) ) &
                             / ( 0.5d0*(xm_global(i+1)-xm_global(i-1)) )
                End If

                ! d^2 R_ij / dz^2  (only if z is not homogeneous)
                dRij_dz = 0d0
                If ( .Not. hom_z .And. k > 1 .And. k < nzm_global ) Then
                   dRij_dz = ( (Rij_g(c,i,j,k+1)-Rij_g(c,i,j,k))/(zm_global(k+1)-zm_global(k)) &
                              -(Rij_g(c,i,j,k)-Rij_g(c,i,j,k-1))/(zm_global(k)-zm_global(k-1)) ) &
                             / ( 0.5d0*(zm_global(k+1)-zm_global(k-1)) )
                End If

                Dnuij_g(c,i,j,k) = nu * (dRij_dx + dRij_dy + dRij_dz)
             End Do

             ! correct epsRes: subtract mean-gradient contribution, eps_ij(fluct) = 2nu*<(dui/dxk)(duj/dxk)> - 2nu*(d<ui>/dxk)(d<uj>/dxk)
             epsRes_g(1,i,j,k) = epsRes_g(1,i,j,k) - 2d0*nu*(dUdx**2 + dUdy**2 + dUdz**2)
             epsRes_g(2,i,j,k) = epsRes_g(2,i,j,k) - 2d0*nu*(dVdx**2 + dVdy**2 + dVdz**2)
             epsRes_g(3,i,j,k) = epsRes_g(3,i,j,k) - 2d0*nu*(dWdx**2 + dWdy**2 + dWdz**2)
             epsRes_g(4,i,j,k) = epsRes_g(4,i,j,k) - 2d0*nu*(dUdx*dVdx + dUdy*dVdy + dUdz*dVdz)
             epsRes_g(5,i,j,k) = epsRes_g(5,i,j,k) - 2d0*nu*(dUdx*dWdx + dUdy*dWdy + dUdz*dWdz)
             epsRes_g(6,i,j,k) = epsRes_g(6,i,j,k) - 2d0*nu*(dVdx*dWdx + dVdy*dWdy + dVdz*dWdz)

             ! correct PiStrain: subtract mean-pressure x mean-strain-rate, Pi_ij(fluct) = <p(dui/dxj+duj/dxi)> - <p>*(d<ui>/dxj + d<uj>/dxi)
             PiStrain_g(1,i,j,k) = PiStrain_g(1,i,j,k) - 2d0*gP(i,j,k)*dUdx
             PiStrain_g(2,i,j,k) = PiStrain_g(2,i,j,k) - 2d0*gP(i,j,k)*dVdy
             PiStrain_g(3,i,j,k) = PiStrain_g(3,i,j,k) - 2d0*gP(i,j,k)*dWdz
             PiStrain_g(4,i,j,k) = PiStrain_g(4,i,j,k) - gP(i,j,k)*(dUdy + dVdx)
             PiStrain_g(5,i,j,k) = PiStrain_g(5,i,j,k) - gP(i,j,k)*(dUdz + dWdx)
             PiStrain_g(6,i,j,k) = PiStrain_g(6,i,j,k) - gP(i,j,k)*(dVdz + dWdy)

          End Do
       End Do
    End Do

    ! ── Turbulent diffusion: D^T_ij = -d<u_i'u_j'u_k'>/dx_k ─────────────
    DTij_g = 0d0
    Do k = 1, nzm_global
       Do j = 1, nym_global
          Do i = 1, nxm_global
             Do c = 1, N_COMP
                ! y-direction (always)
                If ( j > 1 .And. j < nym_global ) Then
                   dTy_dy = (cTijY(c,i,j+1,k)-cTijY(c,i,j-1,k)) / &
                            (ym_global(j+1)-ym_global(j-1))
                Else
                   dTy_dy = 0d0
                End If

                dTx_dx = 0d0
                If ( .Not. hom_x .And. i > 1 .And. i < nxm_global ) Then
                   dTx_dx = (cTijX(c,i+1,j,k)-cTijX(c,i-1,j,k)) / &
                            (xm_global(i+1)-xm_global(i-1))
                End If

                dTz_dz = 0d0
                If ( .Not. hom_z .And. k > 1 .And. k < nzm_global ) Then
                   dTz_dz = (cTijZ(c,i,j,k+1)-cTijZ(c,i,j,k-1)) / &
                            (zm_global(k+1)-zm_global(k-1))
                End If

                DTij_g(c,i,j,k) = -(dTy_dy + dTx_dx + dTz_dz)
             End Do
          End Do
       End Do
    End Do

    ! pressure diffusion Phi^P_ij=-(d<p'u_i'>/dx_j+d<p'u_j'>/dx_i), using central-moment arrays pu_c/pv_c/pw_c = <pu_k>-<p><u_k>, all 6 components via centred differences
    PhiP_g = 0d0
    Do k = 1, nzm_global
       Do j = 1, nym_global
          Do i = 1, nxm_global

             ! -- pre-compute all 9 directional derivatives at this node -------
             dPux_dx = 0d0
             If ( .Not. hom_x .And. i > 1 .And. i < nxm_global ) &
                dPux_dx = (pu_c(i+1,j,k)-pu_c(i-1,j,k))/(xm_global(i+1)-xm_global(i-1))
             dPuy_dy = 0d0
             If ( j > 1 .And. j < nym_global ) &
                dPuy_dy = (pu_c(i,j+1,k)-pu_c(i,j-1,k))/(ym_global(j+1)-ym_global(j-1))
             dPuz_dz = 0d0
             If ( .Not. hom_z .And. k > 1 .And. k < nzm_global ) &
                dPuz_dz = (pu_c(i,j,k+1)-pu_c(i,j,k-1))/(zm_global(k+1)-zm_global(k-1))

             dPvx_dx = 0d0
             If ( .Not. hom_x .And. i > 1 .And. i < nxm_global ) &
                dPvx_dx = (pv_c(i+1,j,k)-pv_c(i-1,j,k))/(xm_global(i+1)-xm_global(i-1))
             dPvy_dy = 0d0
             If ( j > 1 .And. j < nym_global ) &
                dPvy_dy = (pv_c(i,j+1,k)-pv_c(i,j-1,k))/(ym_global(j+1)-ym_global(j-1))
             dPvz_dz = 0d0
             If ( .Not. hom_z .And. k > 1 .And. k < nzm_global ) &
                dPvz_dz = (pv_c(i,j,k+1)-pv_c(i,j,k-1))/(zm_global(k+1)-zm_global(k-1))

             dPwx_dx = 0d0
             If ( .Not. hom_x .And. i > 1 .And. i < nxm_global ) &
                dPwx_dx = (pw_c(i+1,j,k)-pw_c(i-1,j,k))/(xm_global(i+1)-xm_global(i-1))
             dPwy_dy = 0d0
             If ( j > 1 .And. j < nym_global ) &
                dPwy_dy = (pw_c(i,j+1,k)-pw_c(i,j-1,k))/(ym_global(j+1)-ym_global(j-1))
             dPwz_dz = 0d0
             If ( .Not. hom_z .And. k > 1 .And. k < nzm_global ) &
                dPwz_dz = (pw_c(i,j,k+1)-pw_c(i,j,k-1))/(zm_global(k+1)-zm_global(k-1))

             PhiP_g(1,i,j,k) = -2d0 * dPux_dx              ! -2 d<p'u'>/dx
             PhiP_g(2,i,j,k) = -2d0 * dPvy_dy              ! -2 d<p'v'>/dy
             PhiP_g(3,i,j,k) = -2d0 * dPwz_dz              ! -2 d<p'w'>/dz
             PhiP_g(4,i,j,k) = -(dPvx_dx + dPuy_dy)        ! -(d<p'v'>/dx+d<p'u'>/dy)
             PhiP_g(5,i,j,k) = -(dPwx_dx + dPuz_dz)        ! -(d<p'w'>/dx+d<p'u'>/dz)
             PhiP_g(6,i,j,k) = -(dPwy_dy + dPvz_dz)        ! -(d<p'w'>/dy+d<p'v'>/dz)
          End Do
       End Do
    End Do

    ! ── Homogeneous averaging and assign output ───────────────────────────
    Call hom_avg_4d(Umean_g,    3,      UmeanOut)
    Call hom_avg_4d(Rij_g,      N_COMP, Rij)
    Call hom_avg_4d(Pij_g,      N_COMP, Pij)
    Call hom_avg_4d(epsRes_g,   N_COMP, epsResOut)
    Call hom_avg_4d(epsSGS_g,   N_COMP, epsSGSOut)
    Call hom_avg_4d(PiStrain_g, N_COMP, PiStrainOut)
    Call hom_avg_4d(DTij_g,     N_COMP, DTijOut)
    Call hom_avg_4d(Dnuij_g,    N_COMP, DnuijOut)
    Call hom_avg_4d(PhiP_g,     N_COMP, PhiPijOut)

    ! ── Budget residual: Pij + Pi_ij + D^nu + D^T + PhiP - eps_res - eps_sgs
    !    (advection and unsteady terms are not tracked instantaneously)
    ResidOut = Pij + PiStrainOut + DnuijOut + DTijOut + PhiPijOut &
             - epsResOut - epsSGSOut

    Deallocate(Rij_g, Pij_g, epsRes_g, epsSGS_g, PiStrain_g)
    Deallocate(DTij_g, Dnuij_g, PhiP_g, Umean_g)
    Deallocate(cTijY, cTijX, cTijZ)
    Deallocate(pu_c, pv_c, pw_c)

  End Subroutine compute_budget_terms


  !> hom_avg_4d: average a (nc, nxm_global, nym_global, nzm_global) array along the homogeneous directions, returning (nc, out_nx, out_ny, out_nz)
  Subroutine hom_avg_4d(src, nc, dst)

    Integer(Int32), Intent(In)  :: nc
    Real(Int64),    Intent(In)  :: src(nc, nxm_global, nym_global, nzm_global)
    Real(Int64),    Intent(Out) :: dst(nc, out_nx, out_ny, out_nz)

    Integer(Int32) :: i, j, k, c

    If ( hom_x .And. hom_z ) Then
       ! Average over both x and z → output shape (nc, 1, nym_global, 1)
       Do j = 1, nym_global
          Do c = 1, nc
             dst(c, 1, j, 1) = Sum( src(c, :, j, :) ) / Real(nxm_global * nzm_global, Int64)
          End Do
       End Do

    Else If ( hom_x ) Then
       ! Average over x only → output shape (nc, 1, nym_global, nzm_global)
       Do k = 1, nzm_global
          Do j = 1, nym_global
             Do c = 1, nc
                dst(c, 1, j, k) = Sum( src(c, :, j, k) ) / Real(nxm_global, Int64)
             End Do
          End Do
       End Do

    Else If ( hom_z ) Then
       ! Average over z only → output shape (nc, nxm_global, nym_global, 1)
       Do i = 1, nxm_global
          Do j = 1, nym_global
             Do c = 1, nc
                dst(c, i, j, 1) = Sum( src(c, i, j, :) ) / Real(nzm_global, Int64)
             End Do
          End Do
       End Do

    Else
       ! No averaging → copy straight through
       dst = src

    End If

  End Subroutine hom_avg_4d


  !> write_rsb_slice: append one (nc, out_nx, out_ny, out_nz) slice to rsb_<base>_<tag>.bin in little-endian float64 (direct-access stream)
  Subroutine write_rsb_slice(tag, arr, nc)

    Character(*),   Intent(In) :: tag
    Integer(Int32), Intent(In) :: nc
    Real(Int64),    Intent(In) :: arr(nc, out_nx, out_ny, out_nz)

    Character(300)  :: fpath
    Integer(Int32)  :: funit, ios
    Integer(Int64)  :: rec_len
    Logical         :: fexists

    fpath = Trim(rsb_fileout) // '_' // Trim(tag) // '.bin'

    ! convert to little-endian before writing (overrides the solver-wide -fconvert=big-endian flag)

    Open(newunit=funit, file=Trim(fpath), form='unformatted', access='stream', &
         status='unknown', position='append', convert='little_endian', iostat=ios)
    If ( ios /= 0 ) Then
       Write(*,'(A,A)') ' RSB WARNING: cannot open file ', Trim(fpath)
       Return
    End If
    Write(funit) arr
    Close(funit)

  End Subroutine write_rsb_slice


  !> write_meta: (re)write the .meta companion file with current grid dimensions and sample count, called after every successful output_rsb
  Subroutine write_meta

    Character(300)  :: fpath
    Integer(Int32)  :: funit, ios

    fpath = Trim(rsb_fileout) // '.meta'
    Open(newunit=funit, file=Trim(fpath), form='formatted', &
         status='replace', action='write', iostat=ios)
    If ( ios /= 0 ) Then
       Write(*,'(A)') ' RSB WARNING: cannot write meta file.'
       Return
    End If
    Write(funit,'(A)')      '# fdm-dopamine Reynolds stress budget metadata'
    Write(funit,'(A,A)')    'hom_dir   = ', Trim(rsb_hom_dir)
    Write(funit,'(A,3I8)')  'out_shape = ', out_nx, out_ny, out_nz
    Write(funit,'(A,I8)')   'nxm       = ', nxm_global
    Write(funit,'(A,I8)')   'nym       = ', nym_global
    Write(funit,'(A,I8)')   'nzm       = ', nzm_global
    Write(funit,'(A,I8)')   'out_nx    = ', out_nx
    Write(funit,'(A,I8)')   'out_ny    = ', out_ny
    Write(funit,'(A,I8)')   'out_nz    = ', out_nz
    Write(funit,'(A,I8)')   'nsamples  = ', nsamples
    Write(funit,'(A,I8)')   'rsb_freq  = ', rsb_freq
    Write(funit,'(A,I8)')   'rsb_nstart= ', rsb_nstart
    Write(funit,'(A,I8)')   'last_step = ', istep
    Write(funit,'(A)')      'endian    = little'
    Write(funit,'(A)')      'dtype     = float64'
    Write(funit,'(A)')      '# component order: 11 22 33 12 13 23'
    Write(funit,'(A)')      '# files: Umean(3) Rij Pij epsRes epsSGS PiStrain DTij Dnuij PhiPij Resid'
    Write(funit,'(A)')      '# array layout in each file: (nc, out_nx, out_ny, out_nz, nsamples)'
    Write(funit,'(A)')      '#   i.e. one slice of shape (nc, out_nx, out_ny, out_nz) per sample'
    Close(funit)

  End Subroutine write_meta


  !  Private: parse_meta_nsamples
  !  Read the 'nsamples' line from an open meta file unit.
  Subroutine parse_meta_nsamples(funit, n)

    Integer(Int32), Intent(In)  :: funit
    Integer(Int32), Intent(Out) :: n

    Character(300) :: line
    Integer(Int32) :: ios, eq_pos

    n = 0
    Do
       Read(funit, '(A)', iostat=ios) line
       If ( ios /= 0 ) Exit
       If ( Index(line, 'nsamples') > 0 ) Then
          eq_pos = Index(line, '=')
          If ( eq_pos > 0 ) Read(line(eq_pos+1:), *, iostat=ios) n
          Exit
       End If
    End Do

  End Subroutine parse_meta_nsamples


  !  Private: truncate_rsb_files
  !  Delete (truncate) any pre-existing RSB output files for a fresh run.
  Subroutine truncate_rsb_files

    Character(300)  :: fpath
    Integer(Int32)  :: funit
    Character(10), Dimension(10), Parameter :: tags = [ &
         'Umean     ', 'Rij       ', 'Pij       ', 'epsRes    ', 'epsSGS    ', &
         'PiStrain  ', 'DTij      ', 'Dnuij     ', 'PhiPij    ', 'Resid     ' ]
    Integer(Int32) :: t

    Do t = 1, 10
       fpath = Trim(rsb_fileout) // '_' // Trim(Adjustl(tags(t))) // '.bin'
       Open(newunit=funit, file=Trim(fpath), form='unformatted', &
            access='stream', status='replace', iostat=funit)
       Close(funit)
    End Do
    ! meta file is written/replaced at the first output call
    Write(*,'(A)') ' RSB: truncated pre-existing output files for fresh run.'

  End Subroutine truncate_rsb_files

End Module reynolds_stress_budget
