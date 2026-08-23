!####################################################################################################
! CommunicateAll: packed, non-blocking halo exchange for a 4D array (nx,ny,nz,nv)
! Replaces the old pattern of calling Communicate() once per variable with blocking
! MPI_Send/MPI_Recv. Fixes (per LOG discussion):
!   1) one packed message per face-direction instead of one message per variable
!   2) MPI_Isend/MPI_Irecv + MPI_Waitall instead of blocking MPI_Send/MPI_Recv
!   3) no full-array (nx*ny*nz) copy -- only the noverlap-thick halo slabs are
!      packed/unpacked, and packing is done directly from/into the actual array
!      (no arrayin/arrayout aliasing games)
! Still host-staged (no CUDA-aware MPI) by design -- the !$acc update self/device
! calls remain around the call site in iterate.f90.
!####################################################################################################


subroutine CommunicateAll(arr, nv)
  use mod_params
  use mod_spmd
  implicit none
 
  integer, intent(in) :: nv
  real(kind=8), intent(inout) :: arr(nx,ny,nz,nv)
 
  integer, parameter :: TAG_XPLUS=201, TAG_XMINUS=202, &
                        TAG_YPLUS=203, TAG_YMINUS=204, &
                        TAG_ZPLUS=205, TAG_ZMINUS=206
 
  integer :: ierr, nreq, i,j,k,v,ictr
  integer :: reqs(12)
  integer :: statuses(mpi_status_size,12)
 
  integer :: nxsz, nysz, nzsz   ! packed-buffer sizes for x-,y-,z-face exchanges
 
  real(kind=8), allocatable, save :: sbuf_xl(:), sbuf_xr(:), rbuf_xl(:), rbuf_xr(:)
  real(kind=8), allocatable, save :: sbuf_yl(:), sbuf_yr(:), rbuf_yl(:), rbuf_yr(:)
  real(kind=8), allocatable, save :: sbuf_zl(:), sbuf_zr(:), rbuf_zl(:), rbuf_zr(:)
  logical, save :: firstcall = .true.
  integer, save :: nv_alloc = -1
 
  nxsz = noverlap*ny*nz*nv
  nysz = nx*noverlap*nz*nv
  nzsz = nx*ny*noverlap*nv
 
  !! (Re)allocate persistent buffers the first time, or if nv changes
  !! (nv is always nvars=5 in practice, but this keeps it safe/general)
  if (firstcall .or. nv .ne. nv_alloc) then
    if (allocated(sbuf_xl)) deallocate(sbuf_xl,sbuf_xr,rbuf_xl,rbuf_xr, &
                                        sbuf_yl,sbuf_yr,rbuf_yl,rbuf_yr, &
                                        sbuf_zl,sbuf_zr,rbuf_zl,rbuf_zr)
    allocate(sbuf_xl(nxsz), sbuf_xr(nxsz), rbuf_xl(nxsz), rbuf_xr(nxsz))
    allocate(sbuf_yl(nysz), sbuf_yr(nysz), rbuf_yl(nysz), rbuf_yr(nysz))
    allocate(sbuf_zl(nzsz), sbuf_zr(nzsz), rbuf_zl(nzsz), rbuf_zr(nzsz))
    firstcall = .false.
    nv_alloc = nv
  endif
 
  nreq = 0
 
  !=========================================================
  ! Post all receives FIRST (reduces the chance of unexpected
  ! message buffering / early-arrival overhead on some MPI stacks)
  !=========================================================
  if (pleft  .ge. 0) then
    nreq = nreq+1
    call MPI_IRECV(rbuf_xl, nxsz, mpir8, pleft,  TAG_XPLUS,  mpicom, reqs(nreq), ierr)
  endif
  if (pright .ge. 0) then
    nreq = nreq+1
    call MPI_IRECV(rbuf_xr, nxsz, mpir8, pright, TAG_XMINUS, mpicom, reqs(nreq), ierr)
  endif
  if (pbelow .ge. 0) then
    nreq = nreq+1
    call MPI_IRECV(rbuf_yl, nysz, mpir8, pbelow, TAG_YPLUS,  mpicom, reqs(nreq), ierr)
  endif
  if (pabove .ge. 0) then
    nreq = nreq+1
    call MPI_IRECV(rbuf_yr, nysz, mpir8, pabove, TAG_YMINUS, mpicom, reqs(nreq), ierr)
  endif
  if (pdown  .ge. 0) then
    nreq = nreq+1
    call MPI_IRECV(rbuf_zl, nzsz, mpir8, pdown,  TAG_ZPLUS,  mpicom, reqs(nreq), ierr)
  endif
  if (pup    .ge. 0) then
    nreq = nreq+1
    call MPI_IRECV(rbuf_zr, nzsz, mpir8, pup,    TAG_ZMINUS, mpicom, reqs(nreq), ierr)
  endif
 
  !=========================================================
  ! Pack local boundary slabs (ALL variables in one buffer) and
  ! fire off the matching non-blocking sends
  !=========================================================
  ! ---- +x face: send to pright ----
  if (pright .ge. 0) then
    ictr = 0
    do v=1,nv
      do k=1,nz
        do j=1,ny
          do i=iegrid-noverlap+1,iegrid
            ictr = ictr+1
            sbuf_xr(ictr) = arr(i,j,k,v)
          enddo
        enddo
      enddo
    enddo
    nreq = nreq+1
    call MPI_ISEND(sbuf_xr, nxsz, mpir8, pright, TAG_XPLUS, mpicom, reqs(nreq), ierr)
  endif
 
  ! ---- -x face: send to pleft ----
  if (pleft .ge. 0) then
    ictr = 0
    do v=1,nv
      do k=1,nz
        do j=1,ny
          do i=isgrid,isgrid+noverlap-1
            ictr = ictr+1
            sbuf_xl(ictr) = arr(i,j,k,v)
          enddo
        enddo
      enddo
    enddo
    nreq = nreq+1
    call MPI_ISEND(sbuf_xl, nxsz, mpir8, pleft, TAG_XMINUS, mpicom, reqs(nreq), ierr)
  endif
 
  ! ---- +y face: send to pabove ----
  if (pabove .ge. 0) then
    ictr = 0
    do v=1,nv
      do k=1,nz
        do j=jegrid-noverlap+1,jegrid
          do i=1,nx
            ictr = ictr+1
            sbuf_yr(ictr) = arr(i,j,k,v)
          enddo
        enddo
      enddo
    enddo
    nreq = nreq+1
    call MPI_ISEND(sbuf_yr, nysz, mpir8, pabove, TAG_YPLUS, mpicom, reqs(nreq), ierr)
  endif
 
  ! ---- -y face: send to pbelow ----
  if (pbelow .ge. 0) then
    ictr = 0
    do v=1,nv
      do k=1,nz
        do j=jsgrid,jsgrid+noverlap-1
          do i=1,nx
            ictr = ictr+1
            sbuf_yl(ictr) = arr(i,j,k,v)
          enddo
        enddo
      enddo
    enddo
    nreq = nreq+1
    call MPI_ISEND(sbuf_yl, nysz, mpir8, pbelow, TAG_YMINUS, mpicom, reqs(nreq), ierr)
  endif
 
  ! ---- +z face: send to pup ----
  if (pup .ge. 0) then
    ictr = 0
    do v=1,nv
      do k=kegrid-noverlap+1,kegrid
        do j=1,ny
          do i=1,nx
            ictr = ictr+1
            sbuf_zr(ictr) = arr(i,j,k,v)
          enddo
        enddo
      enddo
    enddo
    nreq = nreq+1
    call MPI_ISEND(sbuf_zr, nzsz, mpir8, pup, TAG_ZPLUS, mpicom, reqs(nreq), ierr)
  endif
 
  ! ---- -z face: send to pdown ----
  if (pdown .ge. 0) then
    ictr = 0
    do v=1,nv
      do k=ksgrid,ksgrid+noverlap-1
        do j=1,ny
          do i=1,nx
            ictr = ictr+1
            sbuf_zl(ictr) = arr(i,j,k,v)
          enddo
        enddo
      enddo
    enddo
    nreq = nreq+1
    call MPI_ISEND(sbuf_zl, nzsz, mpir8, pdown, TAG_ZMINUS, mpicom, reqs(nreq), ierr)
  endif
 
  !=========================================================
  ! Wait for everything (sends+recvs) to complete
  !=========================================================
  call MPI_WAITALL(nreq, reqs, statuses, ierr)
 
  !=========================================================
  ! Unpack received buffers into the ghost regions
  !=========================================================
  if (pleft .ge. 0) then
    ictr = 0
    do v=1,nv
      do k=1,nz
        do j=1,ny
          do i=1,isgrid-1
            ictr = ictr+1
            arr(i,j,k,v) = rbuf_xl(ictr)
          enddo
        enddo
      enddo
    enddo
  endif
 
  if (pright .ge. 0) then
    ictr = 0
    do v=1,nv
      do k=1,nz
        do j=1,ny
          do i=iegrid+1,nx
            ictr = ictr+1
            arr(i,j,k,v) = rbuf_xr(ictr)
          enddo
        enddo
      enddo
    enddo
  endif
 
  if (pbelow .ge. 0) then
    ictr = 0
    do v=1,nv
      do k=1,nz
        do j=1,jsgrid-1
          do i=1,nx
            ictr = ictr+1
            arr(i,j,k,v) = rbuf_yl(ictr)
          enddo
        enddo
      enddo
    enddo
  endif
 
  if (pabove .ge. 0) then
    ictr = 0
    do v=1,nv
      do k=1,nz
        do j=jegrid+1,ny
          do i=1,nx
            ictr = ictr+1
            arr(i,j,k,v) = rbuf_yr(ictr)
          enddo
        enddo
      enddo
    enddo
  endif
 
  if (pdown .ge. 0) then
    ictr = 0
    do v=1,nv
      do k=1,ksgrid-1
        do j=1,ny
          do i=1,nx
            ictr = ictr+1
            arr(i,j,k,v) = rbuf_zl(ictr)
          enddo
        enddo
      enddo
    enddo
  endif
 
  if (pup .ge. 0) then
    ictr = 0
    do v=1,nv
      do k=kegrid+1,nz
        do j=1,ny
          do i=1,nx
            ictr = ictr+1
            arr(i,j,k,v) = rbuf_zr(ictr)
          enddo
        enddo
      enddo
    enddo
  endif
 
end subroutine CommunicateAll
 
