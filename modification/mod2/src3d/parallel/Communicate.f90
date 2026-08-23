!####################################################################################################
! periodic/axisymmetric to be updated like 2D
! Temporay array. Need to implement own like in ANUROOP
! Implemented for 3D by RR on 15Feb 2020
subroutine Communicate(arrayin,arrayout)
use mod_params 
use mod_spmd 

implicit none
!include 'mpif.h'

integer :: ierr, status(MPI_STATUS_SIZE)
real(kind=8),intent(in)::arrayin(nx,ny,nz)
real(kind=8),intent(out)::arrayout(nx,ny,nz)

arrayout=arrayin
!##########################################################################################
!##########################################################################################
!##########################################################################################


!if (np_x .ge. 1)then   !NEW LINE --- LOVESH
!Send to east and receieve from west
if (mod(myrank,2) .eq. 1) then
   if (pright .ge. 0) then
      call MPI_send(arrayout(iegrid-noverlap+1:iegrid,1:ny,1:nz), ny*nz*noverlap, mpir8, pright, 15, mpicom, ierr)
      !print*,myrank,'sending to',pright,nyp,arrayout(nxp-1,5)
   endif
endif

if (mod(myrank,2) .eq. 0) then
   if (pleft .ge. 0) then
      call MPI_recv(arrayout(1:isgrid-1,1:ny,1:nz), ny*nz*noverlap, mpir8, pleft, 15, mpicom,status, ierr )
      !print*,myrank,'receiving from',pleft,nyp,pass_dummyout2c(1,5)  
   endif
endif

if (mod(myrank,2) .eq. 0) then
   if (pright .ge. 0) then
      call MPI_send(arrayout(iegrid-noverlap+1:iegrid,1:ny,1:nz), ny*nz*noverlap, mpir8, pright, 15, mpicom, ierr)
      !print*,myrank,'sending to',pright,nyp,arrayout(nxp-1,5)
   endif
endif

if (mod(myrank,2) .eq. 1) then
   if (pleft .ge. 0) then
      call MPI_recv(arrayout(1:isgrid-1,1:ny,1:nz), ny*nz*noverlap, mpir8, pleft, 15, mpicom,status, ierr )   
      !print*,myrank,'receiving from',pleft,nyp,pass_dummyout2c(1,5)
   endif
endif


!Send to west and receieve from east

if (mod(myrank,2) .eq. 1) then
   if (pleft .ge. 0) then
      call MPI_send(arrayout(isgrid:isgrid+noverlap-1,1:ny,1:nz), ny*nz*noverlap, mpir8, pleft, 15, mpicom, ierr)
      !print*,myrank,'sending to',pleft,nyp,arrayout(2,5)
   endif
endif

if (mod(myrank,2) .eq. 0) then
   if (pright .ge. 0) then
      call MPI_recv(arrayout(iegrid+1:nx,1:ny,1:nz), ny*nz*noverlap, mpir8, pright, 15, mpicom,status, ierr )
      !print*,myrank,'receiving from',pright,nyp,pass_dummyout2c(nxp,5)
   endif
endif

if (mod(myrank,2) .eq. 0) then
   if (pleft .ge. 0) then
      call MPI_send(arrayout(isgrid:isgrid+noverlap-1,1:ny,1:nz), ny*nz*noverlap, mpir8, pleft, 15, mpicom, ierr)
      !print*,myrank,'sending to',pleft,nyp,arrayout(2,5)
   endif
endif

if (mod(myrank,2) .eq. 1) then
   if (pright .ge. 0) then
      call MPI_recv(arrayout(iegrid+1:nx,1:ny,1:nz), ny*nz*noverlap, mpir8, pright, 15, mpicom,status, ierr )
      !print*,myrank,'receiving from',pright,nyp,pass_dummyout2c(nxp,5)
   endif
endif

!endif !NEW BLOCK -LOVESH

!call MPI_BARRIER(mpicom, ierr)

!############################################################################################################

!Send to north and receieve from south
!if(np_y .eq. 1)then
if (mod(myrank,2) .eq. 1) then
   if (pabove .ge. 0) then
      call MPI_send(arrayout(1:nx,jegrid-noverlap+1:jegrid,1:nz), nx*nz*noverlap, mpir8, pabove, 15, mpicom, ierr)
      !print*,myrank,'sending to',pabove,nxp,arrayout(5,nyp-1)
   endif
endif

if (mod(myrank,2) .eq. jcomodd) then
   if (pbelow .ge. 0) then
      call MPI_recv(arrayout(1:nx,1:jsgrid-1,1:nz), nx*nz*noverlap, mpir8, pbelow, 15, mpicom,status, ierr )
      !print*,myrank,'receiving from',pbelow,nxp,pass_dummyout2c(5,1)
   endif
endif

if (mod(myrank,2) .eq. 0) then
   if (pabove .ge. 0) then
      call MPI_send(arrayout(1:nx,jegrid-noverlap+1:jegrid,1:nz), nx*nz*noverlap, mpir8, pabove, 15, mpicom, ierr)
      !print*,myrank,'sending to',pabove,nxp,arrayout(5,nyp-1)
   endif
endif

if (mod(myrank,2) .eq. jcomeven) then
   if (pbelow .ge. 0) then
      call MPI_recv(arrayout(1:nx,1:jsgrid-1,1:nz), nx*nz*noverlap, mpir8, pbelow, 15, mpicom,status, ierr )
      !print*,myrank,'receiving from',pbelow,nxp,pass_dummyout2c(5,1)   
   endif
endif


!Send to south and receieve from north
if (mod(myrank,2) .eq. 1) then
   if (pbelow .ge. 0) then
      call MPI_send(arrayout(1:nx,jsgrid:jsgrid+noverlap-1,1:nz), nx*nz*noverlap, mpir8, pbelow, 15, mpicom, ierr)
      !print*,myrank,'sending to',pbelow,pass_dummyout2c(1,2,5)
   endif
endif
if (mod(myrank,2) .eq. jcomodd) then
   if (pabove .ge. 0) then
      !print*,myrank,'receiving from',pabove
      call MPI_recv(arrayout(1:nx,jegrid+1:ny,1:nz), nx*nz*noverlap, mpir8, pabove, 15, mpicom,status, ierr )
   endif
endif

if (mod(myrank,2) .eq. 0) then
   if (pbelow .ge. 0) then
       call MPI_send(arrayout(1:nx,jsgrid:jsgrid+noverlap-1,1:nz), nx*nz*noverlap, mpir8, pbelow, 15, mpicom, ierr)
      !print*,myrank,'sending to',pbelow,nxp,pass_dummyout(5,2)
   endif
endif

if (mod(myrank,2) .eq. jcomeven) then
   if (pabove .ge. 0) then      
      call MPI_recv(arrayout(1:nx,jegrid+1:ny,1:nz), nx*nz*noverlap, mpir8, pabove, 15, mpicom,status, ierr )   
      !print*,myrank,'receiving from',pabove,nxp,pass_dummyout2(5,nyp)
   endif
endif

!endif
!call MPI_BARRIER(mpicom, ierr)
!############################################################################################################
!Send to top and receieve from bottom
!if( np_z .eq. 1)then
if (mod(myrank,2) .eq. 1) then
   if (pup .ge. 0) then
      call MPI_send(arrayout(1:nx,1:ny,kegrid-noverlap+1:kegrid), nx*ny*noverlap, mpir8, pup, 15, mpicom, ierr)
      !print*,myrank,'sending to',pup,nxp,arrayout(5,nyp-1)
   endif
endif

if (mod(myrank,2) .eq. kcomodd) then
   if (pdown .ge. 0) then
      call MPI_recv(arrayout(1:nx,1:ny,1:ksgrid-1), nx*ny*noverlap, mpir8, pdown, 15, mpicom,status, ierr )
      !print*,myrank,'receiving from',pdown,nxp,pass_dummyout2c(5,1)
   endif
endif

if (mod(myrank,2) .eq. 0) then
   if (pup .ge. 0) then
      call MPI_send(arrayout(1:nx,1:ny,kegrid-noverlap+1:kegrid), nx*ny*noverlap, mpir8, pup, 15, mpicom, ierr)
      !print*,myrank,'sending to',pup,nxp,arrayout(5,nyp-1)
   endif
endif

if (mod(myrank,2) .eq. kcomeven) then
   if (pdown .ge. 0) then
      call MPI_recv(arrayout(1:nx,1:ny,1:ksgrid-1), nx*ny*noverlap, mpir8, pdown, 15, mpicom,status, ierr )
      !print*,myrank,'receiving from',pdown,nxp,pass_dummyout2c(5,1)   
   endif
endif


!Send to down and receieve from up
if (mod(myrank,2) .eq. 1) then
   if (pdown .ge. 0) then
      call MPI_send(arrayout(1:nx,1:ny,ksgrid:ksgrid+noverlap-1), nx*ny*noverlap, mpir8, pdown, 15, mpicom, ierr)
      !print*,myrank,'sending to',pdown,pass_dummyout2c(1,2,5)
   endif
endif
if (mod(myrank,2) .eq. kcomodd) then
   if (pup .ge. 0) then
      !print*,myrank,'receiving from',pup
      call MPI_recv(arrayout(1:nx,1:ny,kegrid+1:nz), nx*ny*noverlap, mpir8, pup, 15, mpicom,status, ierr )
   endif
endif

if (mod(myrank,2) .eq. 0) then
   if (pdown .ge. 0) then
       call MPI_send(arrayout(1:nx,1:ny,ksgrid:ksgrid+noverlap-1), nx*ny*noverlap, mpir8, pdown, 15, mpicom, ierr)
      !print*,myrank,'sending to',pdown,nxp,pass_dummyout(5,2)
   endif
endif

if (mod(myrank,2) .eq. kcomeven) then
   if (pup .ge. 0) then      
      call MPI_recv(arrayout(1:nx,1:ny,kegrid+1:nz), nx*ny*noverlap, mpir8, pup, 15, mpicom,status, ierr )   
      !print*,myrank,'receiving from',pup,nxp,pass_dummyout2(5,nyp)
   endif
endif
!endif

!call MPI_BARRIER(mpicom, ierr)


!##########################################################################################
!##########################################################################################
!##########################################################################################

end subroutine Communicate
