! working on this
SUBROUTINE OutputMean3D(itr)
!SUBROUTINE OutputMean3D
  use mod_params
  use mod_spmd
  implicit none
  integer, intent(in) :: itr
  integer:: i,j,k,ivar,ictr,ierr,nsize, iproc
  integer :: nvarst, length
  character(len=60) :: filename
  integer:: mpistat(mpi_status_size)
  real(kind=8), allocatable, dimension(:,:,:) :: uu, vv, ww, uv, vw, uw !rhou, rhov

  nvarst = nvars+6
  if(.not. masterproc)then
    nsize = nvarst*(iegrid-isgrid+1)*(jegrid-jsgrid+1)*(kegrid-ksgrid+1)
    allocate(send_array(nsize))
    ictr = 0
    do ivar =1,nvarst
      do k=ksgrid,kegrid
        do j=jsgrid,jegrid
          do i=isgrid,iegrid
            ictr = ictr+1
            send_array(ictr) = primmean(i,j,k,ivar)
          enddo
        enddo
      enddo
    enddo
    call MPI_SEND(send_array,nsize,mpir8,0,15,mpicom,ierr)
    deallocate(send_array)
  endif

  if(masterproc)then
          allocate(uu(nxg,nyg,nzg))
          allocate(vv(nxg,nyg,nzg))
          allocate(ww(nxg,nyg,nzg))
          allocate(uv(nxg,nyg,nzg))
          allocate(vw(nxg,nyg,nzg))
          allocate(uw(nxg,nyg,nzg))
    do iproc=2,nprocs
      nsize = nvarst*(ipende(iproc)-ipstarte(iproc)+1)*(jpende(iproc)-jpstarte(iproc)+1)*(kpende(iproc)-kpstarte(iproc)+1) 
      allocate(recv_array(nsize)) 
      call MPI_RECV(recv_array,nsize,mpir8,iproc-1,15,mpicom,mpistat,ierr)
      ictr = 0
      do k=kpstarte(iproc),kpende(iproc)
        do j=jpstarte(iproc),jpende(iproc)
          do i=ipstarte(iproc),ipende(iproc)
            ictr = ictr+1
            rhog(i,j,k) = recv_array(ictr)
          enddo
        enddo
      enddo
      do k=kpstarte(iproc),kpende(iproc)
        do j=jpstarte(iproc),jpende(iproc)
          do i=ipstarte(iproc),ipende(iproc)
            ictr = ictr+1
            ug(i,j,k) = recv_array(ictr) 
          enddo
        enddo
      enddo
      do k=kpstarte(iproc),kpende(iproc)
        do j=jpstarte(iproc),jpende(iproc)
          do i=ipstarte(iproc),ipende(iproc)
            ictr = ictr+1
            vg(i,j,k) = recv_array(ictr) 
          enddo
        enddo
      enddo
      do k=kpstarte(iproc),kpende(iproc)
        do j=jpstarte(iproc),jpende(iproc)
          do i=ipstarte(iproc),ipende(iproc)
            ictr = ictr+1
            wg(i,j,k) = recv_array(ictr) 
          enddo
        enddo
      enddo
      do k=kpstarte(iproc),kpende(iproc)
        do j=jpstarte(iproc),jpende(iproc)
          do i=ipstarte(iproc),ipende(iproc)
            ictr = ictr+1
            Pg(i,j,k) = recv_array(ictr) 
          enddo
        enddo
      enddo
      do k=kpstarte(iproc),kpende(iproc)
        do j=jpstarte(iproc),jpende(iproc)
          do i=ipstarte(iproc),ipende(iproc)
            ictr = ictr+1
            uu(i,j,k) = recv_array(ictr) 
          enddo
        enddo
      enddo
      do k=kpstarte(iproc),kpende(iproc)
        do j=jpstarte(iproc),jpende(iproc)
          do i=ipstarte(iproc),ipende(iproc)
            ictr = ictr+1
            vv(i,j,k) = recv_array(ictr) 
          enddo
        enddo
      enddo
      do k=kpstarte(iproc),kpende(iproc)
        do j=jpstarte(iproc),jpende(iproc)
          do i=ipstarte(iproc),ipende(iproc)
            ictr = ictr+1
            ww(i,j,k) = recv_array(ictr) 
          enddo
        enddo
      enddo
      do k=kpstarte(iproc),kpende(iproc)
        do j=jpstarte(iproc),jpende(iproc)
          do i=ipstarte(iproc),ipende(iproc)
            ictr = ictr+1
            uv(i,j,k) = recv_array(ictr) 
          enddo
        enddo
      enddo
      do k=kpstarte(iproc),kpende(iproc)
        do j=jpstarte(iproc),jpende(iproc)
          do i=ipstarte(iproc),ipende(iproc)
            ictr = ictr+1
            vw(i,j,k) = recv_array(ictr) 
          enddo
        enddo
      enddo
      do k=kpstarte(iproc),kpende(iproc)
        do j=jpstarte(iproc),jpende(iproc)
          do i=ipstarte(iproc),ipende(iproc)
            ictr = ictr+1
            uw(i,j,k) = recv_array(ictr) 
          enddo
        enddo
      enddo
      deallocate(recv_array)
    enddo

 ! allocate(tmp(nxg,nyg,nzg))
    rhog(ipstarte(1):ipende(1),jpstarte(1):jpende(1),kpstarte(1):kpend(1))=&
            primmean(isgrid:iegrid,jsgrid:jegrid,ksgrid:kegrid,1)
    ug(ipstarte(1):ipende(1),jpstarte(1):jpende(1),kpstarte(1):kpend(1))=&
            primmean(isgrid:iegrid,jsgrid:jegrid,ksgrid:kegrid,2)
    vg(ipstarte(1):ipende(1),jpstarte(1):jpende(1),kpstarte(1):kpend(1))=&
            primmean(isgrid:iegrid,jsgrid:jegrid,ksgrid:kegrid,3)
    wg(ipstarte(1):ipende(1),jpstarte(1):jpende(1),kpstarte(1):kpend(1))=&
            primmean(isgrid:iegrid,jsgrid:jegrid,ksgrid:kegrid,4)
    Pg(ipstarte(1):ipende(1),jpstarte(1):jpende(1),kpstarte(1):kpend(1))=&
            primmean(isgrid:iegrid,jsgrid:jegrid,ksgrid:kegrid,5)

    filename='meanflow'
    length=len_trim(filename) 
    write(filename(length+1:length+1),'(A1)') '_'
    length=len_trim(filename) 
    write(filename(length+1:length+7),'(I7.7)') itrst
    length=len_trim(filename) 
    write(filename(length+1:length+1),'(A1)') '_'
    length=len_trim(filename) 
    !write(filename(length+1:length+7),'(I7.7)') itrst+nitr
    write(filename(length+1:length+7),'(I7.7)') itr
    length=len_trim(filename) 
    write(filename(length+1:length+4),'(A3)') '.fn'
    open(16,file=trim(filename),form='unformatted',convert='big_endian')
    !write(16) nblivars
    write(16) nxg,nyg,nzg,nvars
    write(16) ug,vg,wg,Pg,rhog
    close(16)
    print*,filename, 'written'
    uu(ipstarte(1):ipende(1),jpstarte(1):jpende(1),kpstarte(1):kpend(1))=&
            primmean(isgrid:iegrid,jsgrid:jegrid,ksgrid:kegrid,6)
    vv(ipstarte(1):ipende(1),jpstarte(1):jpende(1),kpstarte(1):kpend(1))=&
            primmean(isgrid:iegrid,jsgrid:jegrid,ksgrid:kegrid,7)
    ww(ipstarte(1):ipende(1),jpstarte(1):jpende(1),kpstarte(1):kpend(1))=&
            primmean(isgrid:iegrid,jsgrid:jegrid,ksgrid:kegrid,8)
    uv(ipstarte(1):ipende(1),jpstarte(1):jpende(1),kpstarte(1):kpend(1))=&
            primmean(isgrid:iegrid,jsgrid:jegrid,ksgrid:kegrid,9)
    vw(ipstarte(1):ipende(1),jpstarte(1):jpende(1),kpstarte(1):kpend(1))=&
            primmean(isgrid:iegrid,jsgrid:jegrid,ksgrid:kegrid,10)
    uw(ipstarte(1):ipende(1),jpstarte(1):jpende(1),kpstarte(1):kpend(1))=&
            primmean(isgrid:iegrid,jsgrid:jegrid,ksgrid:kegrid,11)
    filename='meancorrel'
    length=len_trim(filename) 
    write(filename(length+1:length+1),'(A1)') '_'
    length=len_trim(filename) 
    write(filename(length+1:length+7),'(I7.7)') itrst
    length=len_trim(filename) 
    write(filename(length+1:length+1),'(A1)') '_'
    length=len_trim(filename) 
    !write(filename(length+1:length+7),'(I7.7)') itrst+nitr
    write(filename(length+1:length+7),'(I7.7)') itr
    length=len_trim(filename) 
    write(filename(length+1:length+4),'(A3)') '.fn'
    open(16,file=trim(filename),form='unformatted',convert='big_endian')
    !write(16) nblivars
    write(16) nxg,nyg,nzg,nvarst-nvars
    write(16) uu,vv,ww,uv,vw,uw
    close(16)
  endif

  !write(16) ((u,i=1,nx),j=1,ny),&
  !  & ((v,i=1,nx),j=1,ny),&
  !  & ((P,i=1,nx),j=1,ny),&
  !  & ((rho,i=1,nx),j=1,ny)

END SUBROUTINE OutputMean3D
