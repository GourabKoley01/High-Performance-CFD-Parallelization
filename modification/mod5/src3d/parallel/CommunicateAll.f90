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
  integer :: nxsz, nysz, nzsz

  real(kind=8), allocatable, save :: sbuf_xl(:), sbuf_xr(:), rbuf_xl(:), rbuf_xr(:)
  real(kind=8), allocatable, save :: sbuf_yl(:), sbuf_yr(:), rbuf_yl(:), rbuf_yr(:)
  real(kind=8), allocatable, save :: sbuf_zl(:), sbuf_zr(:), rbuf_zl(:), rbuf_zr(:)
  logical, save :: firstcall = .true.
  integer, save :: nv_alloc = -1

  nxsz = noverlap*ny*nz*nv
  nysz = nx*noverlap*nz*nv
  nzsz = nx*ny*noverlap*nv

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

  ! Ensure async compute queues are finished before accessing host memory
  !$acc wait

  ! =========================================================
  ! PHASE 1: X-DIRECTION EXCHANGE
  ! =========================================================
  if (pright .ge. 0) then
    !$acc update self(arr(iegrid-noverlap+1:iegrid,1:ny,1:nz,1:nv))
  endif
  if (pleft .ge. 0) then
    !$acc update self(arr(isgrid:isgrid+noverlap-1,1:ny,1:nz,1:nv))
  endif

  nreq = 0
  if (pleft .ge. 0) then
    nreq = nreq+1
    call MPI_IRECV(rbuf_xl, nxsz, mpir8, pleft, TAG_XPLUS, mpicom, reqs(nreq), ierr)
  endif
  if (pright .ge. 0) then
    nreq = nreq+1
    call MPI_IRECV(rbuf_xr, nxsz, mpir8, pright, TAG_XMINUS, mpicom, reqs(nreq), ierr)
  endif

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

  call MPI_WAITALL(nreq, reqs, statuses, ierr)

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
    !$acc update device(arr(1:isgrid-1,1:ny,1:nz,1:nv))
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
    !$acc update device(arr(iegrid+1:nx,1:ny,1:nz,1:nv))
  endif

  ! =========================================================
  ! PHASE 2: Y-DIRECTION EXCHANGE (Now contains fresh X corners)
  ! =========================================================
  if (pabove .ge. 0) then
    !$acc update self(arr(1:nx,jegrid-noverlap+1:jegrid,1:nz,1:nv))
  endif
  if (pbelow .ge. 0) then
    !$acc update self(arr(1:nx,jsgrid:jsgrid+noverlap-1,1:nz,1:nv))
  endif

  nreq = 0
  if (pbelow .ge. 0) then
    nreq = nreq+1
    call MPI_IRECV(rbuf_yl, nysz, mpir8, pbelow, TAG_YPLUS, mpicom, reqs(nreq), ierr)
  endif
  if (pabove .ge. 0) then
    nreq = nreq+1
    call MPI_IRECV(rbuf_yr, nysz, mpir8, pabove, TAG_YMINUS, mpicom, reqs(nreq), ierr)
  endif

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

  call MPI_WAITALL(nreq, reqs, statuses, ierr)

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
    !$acc update device(arr(1:nx,1:jsgrid-1,1:nz,1:nv))
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
    !$acc update device(arr(1:nx,jegrid+1:ny,1:nz,1:nv))
  endif

  ! =========================================================
  ! PHASE 3: Z-DIRECTION EXCHANGE (Now contains fresh X and Y corners)
  ! =========================================================
  if (pup .ge. 0) then
    !$acc update self(arr(1:nx,1:ny,kegrid-noverlap+1:kegrid,1:nv))
  endif
  if (pdown .ge. 0) then
    !$acc update self(arr(1:nx,1:ny,ksgrid:ksgrid+noverlap-1,1:nv))
  endif

  nreq = 0
  if (pdown .ge. 0) then
    nreq = nreq+1
    call MPI_IRECV(rbuf_zl, nzsz, mpir8, pdown, TAG_ZPLUS, mpicom, reqs(nreq), ierr)
  endif
  if (pup .ge. 0) then
    nreq = nreq+1
    call MPI_IRECV(rbuf_zr, nzsz, mpir8, pup, TAG_ZMINUS, mpicom, reqs(nreq), ierr)
  endif

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

  call MPI_WAITALL(nreq, reqs, statuses, ierr)

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
    !$acc update device(arr(1:nx,1:ny,1:ksgrid-1,1:nv))
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
    !$acc update device(arr(1:nx,1:ny,kegrid+1:nz,1:nv))
  endif

end subroutine CommunicateAll
