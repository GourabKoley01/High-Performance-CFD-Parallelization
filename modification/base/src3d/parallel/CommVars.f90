subroutine CommVars
  use mod_params
  use mod_spmd
  implicit none
  integer:: i,j,k,ivar
  integer :: iproc, ictr, nsize, ierr
  integer :: length
  integer:: nxeach, nyeach, nzeach, neighbor
  integer:: mpistat(mpi_status_size)


  if(masterproc)then
    x(isgrid:iegrid,jsgrid:jegrid,ksgrid:kegrid)=& 
      xg(ipstarte(1):ipende(1),jpstarte(1):jpende(1),kpstarte(1):kpende(1))
    y(isgrid:iegrid,jsgrid:jegrid,ksgrid:kegrid)=& 
      yg(ipstarte(1):ipende(1),jpstarte(1):jpende(1),kpstarte(1):kpende(1))
    z(isgrid:iegrid,jsgrid:jegrid,ksgrid:kegrid)=& 
      zg(ipstarte(1):ipende(1),jpstarte(1):jpende(1),kpstarte(1):kpende(1))
    prim(isgrid:iegrid,jsgrid:jegrid,ksgrid:kegrid,1)=& 
      rhog(ipstarte(1):ipende(1),jpstarte(1):jpende(1),kpstarte(1):kpende(1))
    prim(isgrid:iegrid,jsgrid:jegrid,ksgrid:kegrid,2)=& 
      ug(ipstarte(1):ipende(1),jpstarte(1):jpende(1),kpstarte(1):kpende(1))
    prim(isgrid:iegrid,jsgrid:jegrid,ksgrid:kegrid,3)=& 
      vg(ipstarte(1):ipende(1),jpstarte(1):jpende(1),kpstarte(1):kpende(1))
    prim(isgrid:iegrid,jsgrid:jegrid,ksgrid:kegrid,4)=& 
      wg(ipstarte(1):ipende(1),jpstarte(1):jpende(1),kpstarte(1):kpende(1))
    prim(isgrid:iegrid,jsgrid:jegrid,ksgrid:kegrid,5)=& 
      Pg(ipstarte(1):ipende(1),jpstarte(1):jpende(1),kpstarte(1):kpende(1))

    do iproc=2,nprocs
      nsize = (nvars+ndim)*(ipende(iproc)-ipstarte(iproc)+1)*&
        (jpende(iproc)-jpstarte(iproc)+1)*& 
        (kpende(iproc)-kpstarte(iproc)+1) 
      allocate(send_array(nsize)) 
      ictr = 0
      do k=kpstarte(iproc),kpende(iproc)
        do j=jpstarte(iproc),jpende(iproc)
          do i=ipstarte(iproc),ipende(iproc)
            ictr = ictr+1
            send_array(ictr) = xg(i,j,k)
          enddo
        enddo
      enddo
      do k=kpstarte(iproc),kpende(iproc)
        do j=jpstarte(iproc),jpende(iproc)
          do i=ipstarte(iproc),ipende(iproc)
            ictr = ictr+1
            send_array(ictr) = yg(i,j,k)
          enddo
        enddo
      enddo
      do k=kpstarte(iproc),kpende(iproc)
        do j=jpstarte(iproc),jpende(iproc)
          do i=ipstarte(iproc),ipende(iproc)
            ictr = ictr+1
            send_array(ictr) = zg(i,j,k)
          enddo
        enddo
      enddo
      do k=kpstarte(iproc),kpende(iproc)
        do j=jpstarte(iproc),jpende(iproc)
          do i=ipstarte(iproc),ipende(iproc)
            ictr = ictr+1
            send_array(ictr) = rhog(i,j,k)
          enddo
        enddo
      enddo
      do k=kpstarte(iproc),kpende(iproc)
        do j=jpstarte(iproc),jpende(iproc)
          do i=ipstarte(iproc),ipende(iproc)
            ictr = ictr+1
            send_array(ictr) = ug(i,j,k)
          enddo
        enddo
      enddo
      do k=kpstarte(iproc),kpende(iproc)
        do j=jpstarte(iproc),jpende(iproc)
          do i=ipstarte(iproc),ipende(iproc)
            ictr = ictr+1
            send_array(ictr) = vg(i,j,k)
          enddo
        enddo
      enddo
      do k=kpstarte(iproc),kpende(iproc)
        do j=jpstarte(iproc),jpende(iproc)
          do i=ipstarte(iproc),ipende(iproc)
            ictr = ictr+1
            send_array(ictr) = wg(i,j,k)
          enddo
        enddo
      enddo
      do k=kpstarte(iproc),kpende(iproc)
        do j=jpstarte(iproc),jpende(iproc)
          do i=ipstarte(iproc),ipende(iproc)
            ictr = ictr+1
            send_array(ictr) = Pg(i,j,k)
          enddo
        enddo
      enddo
      call MPI_SEND(send_array,nsize,mpir8,iproc-1,15,mpicom,ierr)
      deallocate(send_array)
    enddo
  endif

  if(.not. masterproc)then
    !nsize = nvars*nx*ny
    nsize = (nvars+ndim)*(iegrid-isgrid+1)*&
      (jegrid-jsgrid+1)*(kegrid-ksgrid+1)
    allocate(recv_array(nsize))
    call MPI_RECV(recv_array,nsize,mpir8,0,15,mpicom,mpistat,ierr)
    ictr = 0
    !do k=1,nvars
    do k=ksgrid,kegrid
      do j=jsgrid,jegrid
        do i=isgrid,iegrid
          ictr = ictr+1
          x(i,j,k) = recv_array(ictr)
        enddo
      enddo
    enddo
    do k=ksgrid,kegrid
      do j=jsgrid,jegrid
        do i=isgrid,iegrid
          ictr = ictr+1
          y(i,j,k) = recv_array(ictr)
        enddo
      enddo
    enddo
    do k=ksgrid,kegrid
      do j=jsgrid,jegrid
        do i=isgrid,iegrid
          ictr = ictr+1
          z(i,j,k) = recv_array(ictr)
        enddo
      enddo
    enddo
    do ivar=1,nvars
      do k=ksgrid,kegrid
        do j=jsgrid,jegrid
          do i=isgrid,iegrid
            ictr = ictr+1
            prim(i,j,k,ivar) = recv_array(ictr)
          enddo
        enddo
      enddo
    enddo
    !enddo
    deallocate(recv_array)
  endif

  !!write back and check

  !print*,'rank, umax, rhomax,Pmax=',myrank, maxval(prim(:,:,2)),maxval(prim(:,:,1)),maxval(prim(:,:,4))
  !print*,'rank, umin, rhomin,Pmin=',myrank, minval(prim(:,:,2)),minval(prim(:,:,1)),minval(prim(:,:,4))

  !xmu = 0.0
  !rhs = 0.0
  call MPI_BARRIER(mpicom, ierr)

  call Communicate(x,x)
  call Communicate(y,y)
  call Communicate(z,z)
  do ivar=1,nvars
    call Communicate(prim(:,:,:,ivar),prim(:,:,:,ivar))
  enddo

   call MPI_BARRIER(mpicom, ierr)
  !only for debug
  !  print*,'myrank,nx,ny,nz',myrank,nx,ny,nz
  !  print*,ipstart(1),ipend(1),isgrid, iegrid
  !  print*,jpstart(1),jpend(1),jsgrid, jegrid
  !  print*,kpstart(1),kpend(1),ksgrid, kegrid
  !  print*,maxval(x),maxval(y),maxval(z)
  !  print*,maxval(xg),maxval(yg),maxval(zg)
end subroutine CommVars
