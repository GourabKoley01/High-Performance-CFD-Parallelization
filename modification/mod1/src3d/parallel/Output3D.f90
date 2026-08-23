SUBROUTINE outputmovie3d(itr)
  use mod_params
  use mod_spmd
  implicit none
  integer, intent(in) :: itr
  integer:: i,j,k,ivar,ictr,ierr,nsize, iproc
  integer :: length
  character(len=60) :: filename
  character(len=20) :: path='./snapshots/'
  integer:: mpistat(mpi_status_size)
  logical :: fileexists
!  real (kind =4), allocatable, dimension(:,:,:) :: tmp

  !write in single precision
  if(.not. masterproc)then
    nsize = nvars*(iegrid-isgrid+1)*(jegrid-jsgrid+1)*(kegrid-ksgrid+1)
    allocate(send_array(nsize))
    ictr = 0
    do ivar =1,nvars
      do k=ksgrid,kegrid
        do j=jsgrid,jegrid
          do i=isgrid,iegrid
            ictr = ictr+1
            send_array(ictr) = prim(i,j,k,ivar)
          enddo
        enddo
      enddo
    enddo
    call MPI_SEND(send_array,nsize,mpir8,0,15,mpicom,ierr)
    deallocate(send_array)
  endif

  if(masterproc)then
    do iproc=2,nprocs
      nsize = nvars*(ipende(iproc)-ipstarte(iproc)+1)*(jpende(iproc)-jpstarte(iproc)+1)*(kpende(iproc)-kpstarte(iproc)+1) 
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
      deallocate(recv_array)
    enddo
    rhog(ipstarte(1):ipende(1),jpstarte(1):jpende(1),kpstarte(1):kpend(1))=prim(:,:,:,1)
    ug(ipstarte(1):ipende(1),jpstarte(1):jpende(1),kpstarte(1):kpend(1))=prim(:,:,:,2)
    vg(ipstarte(1):ipende(1),jpstarte(1):jpende(1),kpstarte(1):kpend(1))=prim(:,:,:,3)
    wg(ipstarte(1):ipende(1),jpstarte(1):jpende(1),kpstarte(1):kpend(1))=prim(:,:,:,4)
    Pg(ipstarte(1):ipende(1),jpstarte(1):jpende(1),kpstarte(1):kpend(1))=prim(:,:,:,5)

 ! allocate(tmp(nxg,nyg,nzg))
    filename='movsp'
    length=len_trim(filename) 
    write(filename(length+1:length+1),'(A1)') '_'
    length=len_trim(filename) 
    write(filename(length+1:length+7),'(I7.7)') itr
    length=len_trim(filename) 
    write(filename(length+1:length+4),'(A3)') '.fn'

    !open(16,file=trim(filename),form='unformatted',convert='big_endian')
    open(16, file=trim(path)//trim(filename),form='UNFORMATTED',convert='BIG_ENDIAN')
    !write(16) nblks
    write(16) nxg,nyg,nzg,nvars
    write(16) real(ug),real(vg),real(wg),real(Pg),real(rhog)
   ! tmp = ug
   ! write(16) tmp
   ! tmp = vg
   ! write(16) tmp
   ! tmp = wg
   ! write(16) tmp
   ! tmp = Pg
   ! write(16) tmp
   ! tmp = rhog
   ! write(16) tmp
  !  write(16) (((real(ug(i,j,k)),i=1,nxg),j=1,nyg),k=1,nzg)
  !  write(16) (((real(vg(i,j,k)),i=1,nxg),j=1,nyg),k=1,nzg)
  !  write(16) (((real(wg(i,j,k)),i=1,nxg),j=1,nyg),k=1,nzg)
  !  write(16) (((real(Pg(i,j,k)),i=1,nxg),j=1,nyg),k=1,nzg)
  !  write(16) (((real(rhog(i,j,k)),i=1,nxg),j=1,nyg),k=1,nzg)
  !  16 format(f6.3)
    close(16)
    print*,filename, 'written'
 ! deallocate(tmp, stat =  ierr1)
  endif

  !write(16) ((u,i=1,nx),j=1,ny),&
  !  & ((v,i=1,nx),j=1,ny),&
  !  & ((P,i=1,nx),j=1,ny),&
  !  & ((rho,i=1,nx),j=1,ny)

END SUBROUTINE outputmovie3d

SUBROUTINE Output3d(itr)
  use mod_params
  use mod_spmd
  implicit none
  integer, intent(in) :: itr
  integer:: i,j,k,ivar,ictr,ierr,nsize, iproc
  integer :: length
  character(len=60) :: filename
  integer:: mpistat(mpi_status_size)
  character(len=20) :: path='./snapshots/'
!  real (kind =4), allocatable, dimension(:,:,:) :: tmp

  !write in single precision
  if(.not. masterproc)then
    nsize = nvars*(iegrid-isgrid+1)*(jegrid-jsgrid+1)*(kegrid-ksgrid+1)
    allocate(send_array(nsize))
    ictr = 0
    do ivar =1,nvars
      do k=ksgrid,kegrid
        do j=jsgrid,jegrid
          do i=isgrid,iegrid
            ictr = ictr+1
            send_array(ictr) = prim(i,j,k,ivar)
          enddo
        enddo
      enddo
    enddo
    call MPI_SEND(send_array,nsize,mpir8,0,15,mpicom,ierr)
    deallocate(send_array)
  endif

  if(masterproc)then
    do iproc=2,nprocs
      nsize = nvars*(ipende(iproc)-ipstarte(iproc)+1)*(jpende(iproc)-jpstarte(iproc)+1)*(kpende(iproc)-kpstarte(iproc)+1) 
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
      deallocate(recv_array)
    enddo
    rhog(ipstarte(1):ipende(1),jpstarte(1):jpende(1),kpstarte(1):kpend(1))=prim(:,:,:,1)
    ug(ipstarte(1):ipende(1),jpstarte(1):jpende(1),kpstarte(1):kpend(1))=prim(:,:,:,2)
    vg(ipstarte(1):ipende(1),jpstarte(1):jpende(1),kpstarte(1):kpend(1))=prim(:,:,:,3)
    wg(ipstarte(1):ipende(1),jpstarte(1):jpende(1),kpstarte(1):kpend(1))=prim(:,:,:,4)
    Pg(ipstarte(1):ipende(1),jpstarte(1):jpende(1),kpstarte(1):kpend(1))=prim(:,:,:,5)

 ! allocate(tmp(nxg,nyg,nzg))
    filename='mov'
    length=len_trim(filename) 
    write(filename(length+1:length+1),'(A1)') '_'
    length=len_trim(filename) 
    write(filename(length+1:length+7),'(I7.7)') itr
    length=len_trim(filename) 
    write(filename(length+1:length+4),'(A3)') '.fn'

    open(16, file=trim(path)//trim(filename),form='UNFORMATTED',convert='BIG_ENDIAN')
    !open(16,file=trim(filename),form='unformatted',convert='big_endian')
    !write(16) nblks
    write(16) nxg,nyg,nzg,nvars
    write(16) ug,vg,wg,Pg,rhog
    close(16)
    print*,filename, 'written'
 ! deallocate(tmp, stat =  ierr1)
  endif
END SUBROUTINE Output3d

SUBROUTINE outputplane2d(itr,iplane,kp)
  use mod_params
  use mod_spmd
  implicit none
  integer, intent(in) :: itr,kp,iplane
  integer:: i,j,k,ivar,ictr,ierr,nsize, iproc
  integer :: length
  character(len=60) :: filename
  integer:: mpistat(mpi_status_size)
  character(len=20) :: path='./snapshots2D/'
  logical :: fileexists
!  real (kind =4), allocatable, dimension(:,:,:) :: tmp

  !write in single precision
  if(.not. masterproc)then
    nsize = nvars*(iegrid-isgrid+1)*(jegrid-jsgrid+1)*(kegrid-ksgrid+1)
    allocate(send_array(nsize))
    ictr = 0
    do ivar =1,nvars
      do k=ksgrid,kegrid
        do j=jsgrid,jegrid
          do i=isgrid,iegrid
            ictr = ictr+1
            send_array(ictr) = prim(i,j,k,ivar)
          enddo
        enddo
      enddo
    enddo
    call MPI_SEND(send_array,nsize,mpir8,0,15,mpicom,ierr)
    deallocate(send_array)
  endif

  if(masterproc)then
    do iproc=2,nprocs
      nsize = nvars*(ipende(iproc)-ipstarte(iproc)+1)*(jpende(iproc)-jpstarte(iproc)+1)*(kpende(iproc)-kpstarte(iproc)+1) 
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
      deallocate(recv_array)
    enddo
    rhog(ipstarte(1):ipende(1),jpstarte(1):jpende(1),kpstarte(1):kpend(1))=prim(:,:,:,1)
    ug(ipstarte(1):ipende(1),jpstarte(1):jpende(1),kpstarte(1):kpend(1))=prim(:,:,:,2)
    vg(ipstarte(1):ipende(1),jpstarte(1):jpende(1),kpstarte(1):kpend(1))=prim(:,:,:,3)
    wg(ipstarte(1):ipende(1),jpstarte(1):jpende(1),kpstarte(1):kpend(1))=prim(:,:,:,4)
    Pg(ipstarte(1):ipende(1),jpstarte(1):jpende(1),kpstarte(1):kpend(1))=prim(:,:,:,5)

 ! allocate(tmp(nxg,nyg,nzg))
 !if(iplane .eq. 1)then
 !        plane = 'i'
! elseif(idir .eq. 2)then
!         plane = 'j'
! else
!         plane = 'k'
! endif

    filename='kmov'
    length=len_trim(filename) 
    write(filename(length+1:length+1),'(A1)') '_'
    length=len_trim(filename) 
    write(filename(length+1:length+7),'(I7.7)') itr
    length=len_trim(filename) 
    write(filename(length+1:length+4),'(A3)') '.fn'

    open(16,file=trim(path)//trim(filename),form='unformatted',convert='big_endian')
    !write(16) nblks
    write(16) nxg,nyg,nvars
    write(16) real(ug(:,:,kp)),real(vg(:,:,kp)),real(wg(:,:,kp)),real(Pg(:,:,kp)),real(rhog(:,:,kp))
    close(16)
    print*,filename, 'written'

      INQUIRE( FILE='kgrid2D.xyz', EXIST=fileexists )
  IF ( fileexists )then
  else
    filename='kgrid2D.xyz'
    open(17,file=trim(path)//trim(filename),form='unformatted',convert='big_endian')
    !write(16) nblks
    write(17) nxg,nyg,nvars
    write(17) real(xg(:,:,kp)),real(yg(:,:,kp))
    close(17)
    endif
   ! tmp = ug
   ! write(16) tmp
   ! tmp = vg
   ! write(16) tmp
   ! tmp = wg
   ! write(16) tmp
   ! tmp = Pg
   ! write(16) tmp
   ! tmp = rhog
   ! write(16) tmp
  !  write(16) (((real(ug(i,j,k)),i=1,nxg),j=1,nyg),k=1,nzg)
  !  write(16) (((real(vg(i,j,k)),i=1,nxg),j=1,nyg),k=1,nzg)
  !  write(16) (((real(wg(i,j,k)),i=1,nxg),j=1,nyg),k=1,nzg)
  !  write(16) (((real(Pg(i,j,k)),i=1,nxg),j=1,nyg),k=1,nzg)
  !  write(16) (((real(rhog(i,j,k)),i=1,nxg),j=1,nyg),k=1,nzg)
  !  16 format(f6.3)
 ! deallocate(tmp, stat =  ierr1)
  endif

  !write(16) ((u,i=1,nx),j=1,ny),&
  !  & ((v,i=1,nx),j=1,ny),&
  !  & ((P,i=1,nx),j=1,ny),&
  !  & ((rho,i=1,nx),j=1,ny)

END SUBROUTINE outputplane2d
