subroutine CommunicateAll(arr, nv)
  use mod_params
  use mod_spmd
  implicit none

  integer, intent(in) :: nv
  real(kind=8), intent(inout) :: arr(nx,ny,nz,nv)

  integer, parameter :: TAG_XPLUS=201, TAG_XMINUS=202, &
                        TAG_YPLUS=203, TAG_YMINUS=204, &
                        TAG_ZPLUS=205, TAG_ZMINUS=206

  integer :: ierr, nreq, i,j,k,v,ii,jj,kk,ictr
  integer :: reqs(12)
  integer :: statuses(mpi_status_size,12)
  integer :: nxsz, nysz, nzsz
  integer :: ibase, jbase, kbase

  real(kind=8), allocatable, save :: sbuf_xl(:), sbuf_xr(:), rbuf_xl(:), rbuf_xr(:)
  real(kind=8), allocatable, save :: sbuf_yl(:), sbuf_yr(:), rbuf_yl(:), rbuf_yr(:)
  real(kind=8), allocatable, save :: sbuf_zl(:), sbuf_zr(:), rbuf_zl(:), rbuf_zr(:)
  logical, save :: firstcall = .true.
  integer, save :: nv_alloc = -1

  nxsz = noverlap*ny*nz*nv
  nysz = nx*noverlap*nz*nv
  nzsz = nx*ny*noverlap*nv

  if (firstcall .or. nv .ne. nv_alloc) then
    if (allocated(sbuf_xl)) then
      !$acc exit data delete(sbuf_xl,sbuf_xr,rbuf_xl,rbuf_xr, &
      !$acc                  sbuf_yl,sbuf_yr,rbuf_yl,rbuf_yr, &
      !$acc                  sbuf_zl,sbuf_zr,rbuf_zl,rbuf_zr)
      deallocate(sbuf_xl,sbuf_xr,rbuf_xl,rbuf_xr, &
                 sbuf_yl,sbuf_yr,rbuf_yl,rbuf_yr, &
                 sbuf_zl,sbuf_zr,rbuf_zl,rbuf_zr)
    endif
    allocate(sbuf_xl(nxsz), sbuf_xr(nxsz), rbuf_xl(nxsz), rbuf_xr(nxsz))
    allocate(sbuf_yl(nysz), sbuf_yr(nysz), rbuf_yl(nysz), rbuf_yr(nysz))
    allocate(sbuf_zl(nzsz), sbuf_zr(nzsz), rbuf_zl(nzsz), rbuf_zr(nzsz))

    ! Create persistent device mirrors for the flat, CONTIGUOUS staging
    ! buffers only. These are small (payload-sized), unlike the old
    ! approach of update self/device-ing irregular slabs of the full
    ! (nx,ny,nz,nv) array.
    !$acc enter data create(sbuf_xl,sbuf_xr,rbuf_xl,rbuf_xr)
    !$acc enter data create(sbuf_yl,sbuf_yr,rbuf_yl,rbuf_yr)
    !$acc enter data create(sbuf_zl,sbuf_zr,rbuf_zl,rbuf_zr)

    firstcall = .false.
    nv_alloc = nv
  endif

  ! Ensure async compute queues are finished before touching arr
  !$acc wait

  ! =========================================================
  ! PHASE 1: X-DIRECTION EXCHANGE
  ! =========================================================
  ! ---- pack (GPU, parallel, writes into flat contiguous buffers) ----
  if (pright .ge. 0) then
    ibase = iegrid-noverlap+1
    !$acc parallel loop collapse(4) present(arr, sbuf_xr)
    do v=1,nv
      do k=1,nz
        do j=1,ny
          do ii=1,noverlap
            ictr = ii + noverlap*((j-1) + ny*((k-1) + nz*(v-1)))
            sbuf_xr(ictr) = arr(ibase+ii-1,j,k,v)
          enddo
        enddo
      enddo
    enddo
    !$acc update self(sbuf_xr)   ! small, contiguous -> fast bulk transfer
  endif

  if (pleft .ge. 0) then
    ibase = isgrid
    !$acc parallel loop collapse(4) present(arr, sbuf_xl)
    do v=1,nv
      do k=1,nz
        do j=1,ny
          do ii=1,noverlap
            ictr = ii + noverlap*((j-1) + ny*((k-1) + nz*(v-1)))
            sbuf_xl(ictr) = arr(ibase+ii-1,j,k,v)
          enddo
        enddo
      enddo
    enddo
    !$acc update self(sbuf_xl)
  endif

  ! ---- MPI (unchanged) ----
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
    nreq = nreq+1
    call MPI_ISEND(sbuf_xr, nxsz, mpir8, pright, TAG_XPLUS, mpicom, reqs(nreq), ierr)
  endif
  if (pleft .ge. 0) then
    nreq = nreq+1
    call MPI_ISEND(sbuf_xl, nxsz, mpir8, pleft, TAG_XMINUS, mpicom, reqs(nreq), ierr)
  endif
  call MPI_WAITALL(nreq, reqs, statuses, ierr)

  ! ---- unpack (GPU, parallel, reads flat contiguous buffers) ----
  if (pleft .ge. 0) then
    !$acc update device(rbuf_xl)
    !$acc parallel loop collapse(4) present(arr, rbuf_xl)
    do v=1,nv
      do k=1,nz
        do j=1,ny
          do ii=1,noverlap
            ictr = ii + noverlap*((j-1) + ny*((k-1) + nz*(v-1)))
            arr(ii,j,k,v) = rbuf_xl(ictr)
          enddo
        enddo
      enddo
    enddo
  endif

  if (pright .ge. 0) then
    !$acc update device(rbuf_xr)
    !$acc parallel loop collapse(4) present(arr, rbuf_xr)
    do v=1,nv
      do k=1,nz
        do j=1,ny
          do ii=1,noverlap
            ictr = ii + noverlap*((j-1) + ny*((k-1) + nz*(v-1)))
            arr(iegrid+ii,j,k,v) = rbuf_xr(ictr)
          enddo
        enddo
      enddo
    enddo
  endif

  ! =========================================================
  ! PHASE 2: Y-DIRECTION EXCHANGE (Now contains fresh X corners)
  ! =========================================================
  if (pabove .ge. 0) then
    jbase = jegrid-noverlap+1
    !$acc parallel loop collapse(4) present(arr, sbuf_yr)
    do v=1,nv
      do k=1,nz
        do jj=1,noverlap
          do i=1,nx
            ictr = i + nx*((jj-1) + noverlap*((k-1) + nz*(v-1)))
            sbuf_yr(ictr) = arr(i,jbase+jj-1,k,v)
          enddo
        enddo
      enddo
    enddo
    !$acc update self(sbuf_yr)
  endif

  if (pbelow .ge. 0) then
    jbase = jsgrid
    !$acc parallel loop collapse(4) present(arr, sbuf_yl)
    do v=1,nv
      do k=1,nz
        do jj=1,noverlap
          do i=1,nx
            ictr = i + nx*((jj-1) + noverlap*((k-1) + nz*(v-1)))
            sbuf_yl(ictr) = arr(i,jbase+jj-1,k,v)
          enddo
        enddo
      enddo
    enddo
    !$acc update self(sbuf_yl)
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
    nreq = nreq+1
    call MPI_ISEND(sbuf_yr, nysz, mpir8, pabove, TAG_YPLUS, mpicom, reqs(nreq), ierr)
  endif
  if (pbelow .ge. 0) then
    nreq = nreq+1
    call MPI_ISEND(sbuf_yl, nysz, mpir8, pbelow, TAG_YMINUS, mpicom, reqs(nreq), ierr)
  endif
  call MPI_WAITALL(nreq, reqs, statuses, ierr)

  if (pbelow .ge. 0) then
    !$acc update device(rbuf_yl)
    !$acc parallel loop collapse(4) present(arr, rbuf_yl)
    do v=1,nv
      do k=1,nz
        do jj=1,noverlap
          do i=1,nx
            ictr = i + nx*((jj-1) + noverlap*((k-1) + nz*(v-1)))
            arr(i,jj,k,v) = rbuf_yl(ictr)
          enddo
        enddo
      enddo
    enddo
  endif

  if (pabove .ge. 0) then
    !$acc update device(rbuf_yr)
    !$acc parallel loop collapse(4) present(arr, rbuf_yr)
    do v=1,nv
      do k=1,nz
        do jj=1,noverlap
          do i=1,nx
            ictr = i + nx*((jj-1) + noverlap*((k-1) + nz*(v-1)))
            arr(i,jegrid+jj,k,v) = rbuf_yr(ictr)
          enddo
        enddo
      enddo
    enddo
  endif

  ! =========================================================
  ! PHASE 3: Z-DIRECTION EXCHANGE (Now contains fresh X and Y corners)
  ! =========================================================
  if (pup .ge. 0) then
    kbase = kegrid-noverlap+1
    !$acc parallel loop collapse(4) present(arr, sbuf_zr)
    do v=1,nv
      do kk=1,noverlap
        do j=1,ny
          do i=1,nx
            ictr = i + nx*((j-1) + ny*((kk-1) + noverlap*(v-1)))
            sbuf_zr(ictr) = arr(i,j,kbase+kk-1,v)
          enddo
        enddo
      enddo
    enddo
    !$acc update self(sbuf_zr)
  endif

  if (pdown .ge. 0) then
    kbase = ksgrid
    !$acc parallel loop collapse(4) present(arr, sbuf_zl)
    do v=1,nv
      do kk=1,noverlap
        do j=1,ny
          do i=1,nx
            ictr = i + nx*((j-1) + ny*((kk-1) + noverlap*(v-1)))
            sbuf_zl(ictr) = arr(i,j,kbase+kk-1,v)
          enddo
        enddo
      enddo
    enddo
    !$acc update self(sbuf_zl)
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
    nreq = nreq+1
    call MPI_ISEND(sbuf_zr, nzsz, mpir8, pup, TAG_ZPLUS, mpicom, reqs(nreq), ierr)
  endif
  if (pdown .ge. 0) then
    nreq = nreq+1
    call MPI_ISEND(sbuf_zl, nzsz, mpir8, pdown, TAG_ZMINUS, mpicom, reqs(nreq), ierr)
  endif
  call MPI_WAITALL(nreq, reqs, statuses, ierr)

  if (pdown .ge. 0) then
    !$acc update device(rbuf_zl)
    !$acc parallel loop collapse(4) present(arr, rbuf_zl)
    do v=1,nv
      do kk=1,noverlap
        do j=1,ny
          do i=1,nx
            ictr = i + nx*((j-1) + ny*((kk-1) + noverlap*(v-1)))
            arr(i,j,kk,v) = rbuf_zl(ictr)
          enddo
        enddo
      enddo
    enddo
  endif

  if (pup .ge. 0) then
    !$acc update device(rbuf_zr)
    !$acc parallel loop collapse(4) present(arr, rbuf_zr)
    do v=1,nv
      do kk=1,noverlap
        do j=1,ny
          do i=1,nx
            ictr = i + nx*((j-1) + ny*((kk-1) + noverlap*(v-1)))
            arr(i,j,kegrid+kk,v) = rbuf_zr(ictr)
          enddo
        enddo
      enddo
    enddo
  endif

end subroutine CommunicateAll
