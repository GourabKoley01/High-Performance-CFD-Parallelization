!####################################################################################################
! CommunicateAll: packed, non-blocking halo exchange for a 4D array (nx,ny,nz,nv)
! Replaces the old pattern of calling Communicate() once per variable with blocking
! MPI_Send/MPI_Recv. Fixes (per LOG discussion):
!   1) one packed message per face-direction instead of one message per variable
!   2) MPI_Isend/MPI_Irecv + MPI_Waitall instead of blocking MPI_Send/MPI_Recv
!   3) no full-array (nx*ny*nz) copy -- only the noverlap-thick halo slabs are
!      packed/unpacked, and packing is done directly from/into the actual array
!      (no arrayin/arrayout aliasing games)
! Still host-staged (no CUDA-aware MPI) BY DESIGN -- this must run correctly on
! machines without GPUDirect/CUDA-aware MPI, so all !$acc update self/device
! calls needed for the halo exchange are now done INSIDE this subroutine, using
! EXACT ranges that match what is packed/unpacked (isgrid/iegrid-based), instead
! of the generic 1:noverlap+1 / nx-noverlap:nx ranges that used to live in
! iterate.f90. Those generic ranges only overlapped the real send/recv slabs by
! a single point and were the likely source of the NaN/Inf blowup -- the old
! blocking, once-per-variable Communicate() call happened to hide this because
! its extra serialization gave outstanding async GPU work time to finish before
! the (wrong-range) host copy was read; the faster consolidated path removed
! that incidental delay and exposed the stale-data race.
!
! An explicit "!$acc wait" is issued before every update self, so we do not
! depend on incidental timing/serialization elsewhere in the code to guarantee
! that the device values we are about to copy to host are actually finished.
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

  !=========================================================
  ! Make sure the device has actually finished writing 'arr'
  ! before we copy any of it to host. Without this, "update self"
  ! is only guaranteed to fetch whatever is currently in device
  ! memory -- if a prior kernel touching arr is still in flight on
  ! an async queue, we can grab stale/half-written values.
  !=========================================================
  !$acc wait

  !=========================================================
  ! Pull EXACTLY the slabs we are about to pack into the send
  ! buffers from device to host. Ranges match the pack loops
  ! below one-for-one (isgrid/iegrid-based, not the old generic
  ! 1:noverlap+1 / nx-noverlap:nx approximation).
  !=========================================================
  if (pright .ge. 0) then
    !$acc update self(arr(iegrid-noverlap+1:iegrid,1:ny,1:nz,1:nv))
  endif
  if (pleft .ge. 0) then
    !$acc update self(arr(isgrid:isgrid+noverlap-1,1:ny,1:nz,1:nv))
  endif
  if (pabove .ge. 0) then
    !$acc update self(arr(1:nx,jegrid-noverlap+1:jegrid,1:nz,1:nv))
  endif
  if (pbelow .ge. 0) then
    !$acc update self(arr(1:nx,jsgrid:jsgrid+noverlap-1,1:nz,1:nv))
  endif
  if (pup .ge. 0) then
    !$acc update self(arr(1:nx,1:ny,kegrid-noverlap+1:kegrid,1:nv))
  endif
  if (pdown .ge. 0) then
    !$acc update self(arr(1:nx,1:ny,ksgrid:ksgrid+noverlap-1,1:nv))
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
  ! fire off the matching non-blocking sends. These loops now
  ! read host data that was just explicitly refreshed above, so
  ! they are guaranteed current.
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
  ! Unpack received buffers into the ghost regions (host side)
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

  !=========================================================
  ! Push the freshly-unpacked ghost regions back to the device.
  ! Ranges match the unpack loops above exactly.
  !=========================================================
  if (pleft .ge. 0) then
    !$acc update device(arr(1:isgrid-1,1:ny,1:nz,1:nv))
  endif
  if (pright .ge. 0) then
    !$acc update device(arr(iegrid+1:nx,1:ny,1:nz,1:nv))
  endif
  if (pbelow .ge. 0) then
    !$acc update device(arr(1:nx,1:jsgrid-1,1:nz,1:nv))
  endif
  if (pabove .ge. 0) then
    !$acc update device(arr(1:nx,jegrid+1:ny,1:nz,1:nv))
  endif
  if (pdown .ge. 0) then
    !$acc update device(arr(1:nx,1:ny,1:ksgrid-1,1:nv))
  endif
  if (pup .ge. 0) then
    !$acc update device(arr(1:nx,1:ny,kegrid+1:nz,1:nv))
  endif

end subroutine CommunicateAll
