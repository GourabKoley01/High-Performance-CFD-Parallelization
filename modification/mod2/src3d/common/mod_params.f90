module mod_params
  implicit none


  !!INPUT PARAMETERS!!
  real(kind =8), parameter :: &
    & GAM = 1.4, &
    & GAM1 = 0.4, &
    & GAM1I = 2.5, &
    & PRANDTL = 0.71, &
    & SMUD = 110.4, & 
    & Rgas = 287.053, &
    & XPI = 2.*ATAN2(1.0,0.0), &
    & onethird = 1.0 / 3.0, &
    & eps = 1.0e-16, &
    & zero = 0.0d0, &
    & hf = 0.50d0, &
    & fth = 0.250d0, &
    & recon_ratio = 1.0d0/3.0d0
  real(kind =8), parameter, dimension(4) :: RK4 = (/0.25, 1.0/3.0, 0.5, 0.5/)
  real(kind =8), parameter, dimension(3) :: RK3A = (/0.0, 3.0/4.0, 1.0/3.0/)
  real(kind =8), parameter, dimension(3) :: RK3B = (/1.0, 1.0/4.0, 2.0/3.0/)
  !Compact Schemes
  !C3-AC4-C4-AC4-C3
  !alpha,a,b,c
  real(kind =8), parameter, dimension(3) :: C2B = (/1.0, -2.0, 2.0/)
  real(kind =8), parameter, dimension(4) :: C3B = (/2.0, -2.5, 2.0, 0.5/)
  real(kind =8), parameter, dimension(4) :: AC4B  = (/0.25, -0.75, 0.0, 0.75/)
  real(kind =8), parameter, dimension(2) :: C4 = (/0.25, 1.5/)
  real(kind =8), parameter, dimension(3) :: C6 = (/1.0/3.0, 14.0/9.0, 1.0/9.0/)
  !Compact Filters
  integer, parameter :: nvars = 5, ndim = 3
  integer, parameter :: wunit = 11, wunitp = 22, wunitm =99
  character(len=100) :: gridFile
  !WENO
  integer, parameter :: wenoexp = 2
  real (kind=8) :: RE, XMACH, GXM, FACTOR, REI, SMU, SMUP1
  real (kind=8) :: Tinf, rhoinf, uinf, vinf, winf, pinf
  integer :: iSolver
  integer :: imetric_scheme, itime_scheme, iflux_scheme, ivis_scheme
  integer :: nitr, itrst
  integer :: movie_freq, iprint_freq, movie2d_freq
  integer :: ifprobe
  integer :: iper, jper, kper
  integer :: kClip, kPlane
  real (kind=8) :: dt
  !!
  integer :: iflux_order
  integer :: ifilt_order 
  integer :: iroe_order 
  integer :: imean, iDebug 
  real(kind=8) :: alpha, a, b, c
  real(kind=8) :: alphaf(6)
  real(kind=8), dimension(6,6) :: afilt
  !!afilt(4,5) = F8(=4*2), a0(1-1), a1, a2, a3, a4, 
  !integer, allocatable, dimension(:) :: mfilt
  !!mfilt decides the order of filter for grid point
  integer :: igfstart !!interior grid point on which filter will operate  
  !!GRID PARAMETERS!!
  integer :: nx, ny, nz, nmax
  !integer, allocatable, dimension(:,:) :: iflag
  real (kind =8), allocatable, dimension(:,:,:) :: x,y, z
  real (kind =8), allocatable, dimension(:,:,:) :: xix,etax,ztax
  real (kind =8), allocatable, dimension(:,:,:) :: xiy,etay,ztay
  real (kind =8), allocatable, dimension(:,:,:) :: xiz,etaz,ztaz
  real (kind =8), allocatable, dimension(:,:,:) :: xjac
 ! real (kind =8), allocatable, dimension(:,:,:) :: xjacmi, xjacmj
 ! real (kind =8), allocatable, dimension(:,:,:) :: xixmi, xiymi, etaxmi, etaymi
 ! real (kind =8), allocatable, dimension(:,:,:) :: xixmj, xiymj, etaxmj, etaymj

  !!FLOW PARAMETERS!!
  real (kind =8), allocatable, dimension(:,:,:) :: T,xmu
  real (kind =8), allocatable, dimension(:,:,:,:) :: rhs
  real (kind =8), allocatable, dimension(:,:,:) :: deltaq
  real (kind =8), allocatable, dimension(:,:,:,:) :: q, q0, prim
  real (kind =8), allocatable, dimension(:,:,:,:) :: primmean

  !!ROE FLUX
  real (kind =8), allocatable, dimension(:,:,:,:) :: primL, primR
  real (kind =8), allocatable, dimension(:,:,:,:) :: xiflux, etflux, ztflux

  !!TDG SYSTEM!!
  real (kind =8), allocatable, dimension(:,:) :: atdg,btdg,ctdg
  real (kind =8), allocatable, dimension(:,:,:) :: rhstdg
  real (kind =8), allocatable, dimension(:,:,:) :: tdg

  !!TDG_field SYSTEM!!
  real (kind =8), allocatable, dimension(:,:,:) :: atdg_field,btdg_field,ctdg_field,gamtdg_field
  real (kind =8), allocatable, dimension(:,:,:,:) :: rhstdg_field
  
  !!MULTIPLE BLOCKS
  integer, allocatable, dimension(:,:,:) :: ifblock, iShock
  integer :: nblnksg
  integer, allocatable, dimension(:) :: isblkg, ieblkg, jsblkg, jeblkg, ksblkg,&
  keblkg
  integer :: nblnks
  integer, allocatable, dimension(:) :: isblk, ieblk, jsblk, jeblk, ksblk, keblk

  !!TEMP VARIABLES!!
  real (kind =8), allocatable, dimension(:,:,:,:) :: temp
  real (kind =8), allocatable, dimension(:,:,:,:) :: temp2
  !real (kind =8), allocatable, dimension(:,:) :: rhsfilt
  integer :: ierr1, ierr2, ierr3
  !!TIME
   real(kind=8) :: &
      start_wtime ! Starting wall-clock time of the code
  !!character, dimension(10) :: chrA
  character(len=30), dimension(10) :: pfile
  integer :: nprobesg, nprobes
  integer, allocatable, dimension(:,:) :: probeg,probe
  real (kind =8), allocatable, dimension(:,:,:) :: probedatag


  !!$acc declare create(T(:,:,:), xmu(:,:,:))
  !!$acc declare create(temp)
  
  !!$acc declare create(primL(:,:,:,:), primR(:,:,:,:), xiflux(:,:,:,:), etflux(:,:,:,:), ztflux(:,:,:,:))
  
  !!$acc declare create(atdg_field(:,:,:), btdg_field(:,:,:), ctdg_field(:,:,:), gamtdg_field(:,:,:), rhstdg_field(:,:,:,:))
  !!$acc declare copyin(C6(:), AC4B(:), C3B(:))

end module
