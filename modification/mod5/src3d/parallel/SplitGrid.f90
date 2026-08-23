! TO DO: Add multiblockGrid.f90 
subroutine SplitGrid
  use mod_params
  use mod_spmd
  implicit none
  integer:: i,j,k,ierr
  integer:: nxeach, nyeach, nzeach, neighbor
  integer:: istart,iend    
  integer:: jstart,jend
  integer:: kstart,kend
  integer:: mpistat(mpi_status_size)
  integer :: ijkctr

  call mpi_inputs

  nprocs = np_x*np_y*np_z
!  print*,'nprocs= ', nprocs
  if(masterproc) then
    print*, np_x, np_y,np_z, nprocs
    allocate(ipstart(nprocs))
    allocate(ipend(nprocs))
    allocate(jpstart(nprocs))
    allocate(jpend(nprocs))
    allocate(kpstart(nprocs))
    allocate(kpend(nprocs))

    allocate(ipstarte(nprocs))
    allocate(ipende(nprocs))
    allocate(jpstarte(nprocs))
    allocate(jpende(nprocs))
    allocate(kpstarte(nprocs))
    allocate(kpende(nprocs))

    nxeach = nxg/np_x
    nyeach = nyg/np_y
    nzeach = nzg/np_z

    !print*, nxeach, nyeach, nzeach

    ijkctr = 0
    istart = 1; iend = 1
    jstart = 1; jend = 1
    kstart = 1; kend = 1
    neighbor = 0
    do k=1,np_z
      do j=1,np_y
        do i=1,np_x
          ijkctr = ijkctr + (iend -istart+1)*(jend-jstart+1)*(kend-kstart+1)
          istart = (i-1)*nxeach +1
          iend  = istart+nxeach-1
          jstart = (j-1)*nyeach +1
          jend  = jstart+nyeach-1
          kstart = (k-1)*nzeach +1
          kend  = kstart+nzeach-1

          if(i.eq.np_x) iend=nxg
          if(j.eq.np_y) jend=nyg
          if(k.eq.np_z) kend=nzg

          isubgrid(1)=i
          isubgrid(2)=istart
          isubgrid(3)=iend
          isubgrid(4)=j
          isubgrid(5)=jstart
          isubgrid(6)=jend
          isubgrid(7)=k
          isubgrid(8)=kstart
          isubgrid(9)=kend
          isubgrid(10)=ijkctr

          !For receiving, no overlap points..Exclusive points
          ipstarte(neighbor+1) = istart 
          ipende(neighbor+1) = iend
          jpstarte(neighbor+1) = jstart 
          jpende(neighbor+1) = jend 
          kpstarte(neighbor+1) = kstart 
          kpende(neighbor+1) = kend 

          !For communication to all, include overlap points
          ipstart(neighbor+1) = istart -noverlap 
          ipend(neighbor+1) = iend + noverlap
          jpstart(neighbor+1) = jstart - noverlap
          jpend(neighbor+1) = jend + noverlap
          kpstart(neighbor+1) = kstart - noverlap
          kpend(neighbor+1) = kend + noverlap


          !if(i.eq.0) ipstart(neighbor+1)=istart
          !if(i.eq.np_x) ipend(neighbor+1)=iend
          !if(j.eq.0) jpstart(neighbor+1)=jstart
          !if(j.eq.np_y) jpend(neighbor+1)=jend

          if(istart.eq.1) ipstart(neighbor+1)=istart
          if(iend.eq.nxg) ipend(neighbor+1)=iend
          if(jstart.eq.1) jpstart(neighbor+1)=jstart
          if(jend.eq.nyg) jpend(neighbor+1)=jend
          if(kstart.eq.1) kpstart(neighbor+1)=kstart
          if(kend.eq.nzg) kpend(neighbor+1)=kend

          if(neighbor .ne. 0) then
            call MPI_SEND(isubgrid,10,mpiint,neighbor,15,mpicom,ierr)
          endif
          neighbor = neighbor +1
        enddo
      enddo
    enddo

    isubgrid(1)=1
    isubgrid(2)=1
    isubgrid(3)=nxeach
    isubgrid(4)=1
    isubgrid(5)=1
    isubgrid(6)=nyeach
    isubgrid(7)=1
    isubgrid(8)=1
    isubgrid(9)=nzeach
    isubgrid(10)=1

    ipstarte(1) = 1
    ipende(1) = nxeach 
    jpstarte(1) = 1
    jpende(1) = nyeach 
    kpstarte(1) = 1
    kpende(1) = nzeach 

    ipstart(1) = 1
    ipend(1) = nxeach + noverlap
    jpstart(1) = 1
    jpend(1) = nyeach + noverlap
    kpstart(1) = 1
    kpend(1) = nzeach + noverlap
  endif

  if(.not. masterproc) then
    call MPI_RECV(isubgrid,10,mpiint,0,15,mpicom,mpistat,ierr)
  endif

  !call MPI_BCAST(nxg,1,mpiint,0,mpicom,ierr)
  !call MPI_BCAST(nyg,1,mpiint,0,mpicom,ierr)
  !print*,'nxg,nyg', myrank, nxg,nyg

  pi = isubgrid(1) !iproc
  isglobal = isubgrid(2) !istart
  ieglobal = isubgrid(3)  !iend 
  pj = isubgrid(4) !jproc
  jsglobal = isubgrid(5) !jstart
  jeglobal = isubgrid(6)  !jend 
  pk = isubgrid(7) !jproc
  ksglobal = isubgrid(8) !jstart
  keglobal = isubgrid(9)  !jend 
  ijkdisp = isubgrid(10) !disp for mpi io

  !print*,'i,nx,ny',myrank,(ieglobal-isglobal+1),(jeglobal-jsglobal+1)

  pleft = myrank-1
  pright = myrank+1
  pabove = myrank+np_x
  pbelow = myrank-np_x
  pup = myrank+np_x*np_y
  pdown = myrank-np_x*np_y

  !!below think it through. In UNNI's code its i,j instead of proc
  if(pi .eq. 1)pleft=-1
  if(pi .eq. np_x)pright=-1
  if(pj .eq. 1)pbelow=-1
  if(pj .eq. np_y)pabove=-1
  if(pk .eq. 1)pdown=-1
  if(pk .eq. np_z)pup=-1
  !if(isglobal .eq. 1)pleft=-1
  !if(ieglobal .eq. nxg)pright=-1
  !if(jsglobal .eq. 1)pbelow=-1
  !if(jeglobal .eq. nyg)pabove=-1

  !!Now include overlap
  nx = (ieglobal - isglobal)+1+2*noverlap
  ny = (jeglobal - jsglobal)+1+2*noverlap
  nz = (keglobal - ksglobal)+1+2*noverlap
  isgrid = noverlap +1
  iegrid = nx - noverlap 
  jsgrid = noverlap +1
  jegrid = ny - noverlap 
  ksgrid = noverlap +1
  kegrid = nz - noverlap 

  !!Boundaries, one sided-overlap
!!Boundaries, one sided-overlap if not periodic/axisymmetric
if(iper .eq. 0)then
  if(pleft .eq. -1) then
    nx = (ieglobal - isglobal)+1+noverlap
    isgrid = 1
    iegrid = nx-noverlap
  endif
  if(pright .eq. -1) then
    nx = (ieglobal - isglobal)+1+noverlap
    isgrid = noverlap+1
    iegrid = nx
  endif

endif
if(jper .eq. 0)then
  if(pabove .eq. -1) then
    ny = (jeglobal - jsglobal)+1+noverlap
    jsgrid = noverlap+1
    jegrid = ny
  endif
  if(pbelow .eq. -1) then
    ny = (jeglobal - jsglobal)+1+noverlap
    jsgrid = 1
    jegrid = ny - noverlap
  endif
endif

if(kper .eq. 0)then
  if(pup .eq. -1) then
    nz = (keglobal - ksglobal)+1+noverlap
    ksgrid = noverlap+1
    kegrid = nz
  endif
  if(pdown .eq. -1) then
    nz = (keglobal - ksglobal)+1+noverlap
    ksgrid = 1
    kegrid = nz - noverlap
  endif
endif
   
  !NEW BLOCK - LOVESH
  !for pencil decompositon
 if(np_y .eq. 1)then
        ny = (jeglobal - jsglobal) + 1
        jsgrid = 1
        jegrid = ny
 endif

 if(np_z .eq. 1)then
        nz= (keglobal - ksglobal) + 1
        ksgrid = 1
        kegrid = nz
 endif 

  !add code to change nx ny nz limits to excude ovelap on both sides
  
 
  allocate(x(nx,ny,nz))
  allocate(y(nx,ny,nz))
  allocate(z(nx,ny,nz))

  !#############
  allocate(ifblock(nx,ny,nz))
  ifblock = 0

  call mpi_barrier(mpicom,ierr)
  if(nblnksg .gt. 0)call MultiblockGrid 
  call mpi_barrier(mpicom,ierr)
  !#############
  !print*,'i,nx,ny',myrank,(ieglobal-isglobal+1),(jeglobal-jsglobal+1)
  !print*,myrank,nx,ny


contains
  subroutine mpi_inputs
    use mod_params
    implicit none

    open(10, file='mpi_input.in',status='OLD', action='READ')
    read(10,*) np_x, np_y, np_z
    !read(10,*) noverlap
    close(10)

  noverlap = 11
  if(iflux_scheme .eq. 2) noverlap=5
    if(mod(np_x,2) .eq. 0)then
      jcomodd =1; jcomeven=0
    endif
    if(mod(np_x,2) .eq. 1)then
      jcomodd =0; jcomeven=1
    endif

    if(mod(np_x*np_y,2) .eq. 0)then
      kcomodd =1; kcomeven=0
    endif
    if(mod(np_x*np_y,2) .eq. 1)then
      kcomodd =0; kcomeven=1
    endif
  end subroutine mpi_inputs
end subroutine SplitGrid
