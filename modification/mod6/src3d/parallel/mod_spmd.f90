module mod_spmd
  use mod_params
   implicit none
  include 'mpif.h'

   logical :: &
      masterproc ! Flag for myrank = 0
   integer :: &
      np_x, & ! Number of processors along x direction
      np_y, & ! Number of processors along y direction
      np_z, & ! Number of processors along z direction
      nprocs, & ! Total number of processors
      myrank, & ! myrank of calling process in communication
      pi, pj, pk, & ! Co-ordinates defining a processor
      mpicom, & ! MPI communicator
      mpicom_pi0, & ! New MPI communicator for myranks with pi = 0
                    ! Credits: Patricia Balle (CRAY Inc.)
      mpinull, & ! For mpi_info_null
      mpiundef, & ! For mpi_undefined
      mpilog, & ! For mpi_logical
      mpiint, & ! For mpi_integer
      mpir8, & ! For 8 byte real: mpi_double/mpi_real8
      mpisum, & ! For mpi_sum
      mpimax ! For mpi_max

    integer :: &
      nxg, nyg, nzg
      
    real(kind=8), allocatable, dimension(:,:,:) :: xg,yg,zg, ug, vg,wg, Pg, rhog
    real(kind=8), allocatable, dimension(:) :: recv_array, send_array

    integer, allocatable, dimension(:) :: &
      ipstart,ipend,jpstart,jpend,kpstart, kpend
    integer, allocatable, dimension(:) :: &
      ipstarte,ipende,jpstarte,jpende, kpstarte, kpende

    integer :: &
      pleft, pright, pabove, pbelow,pup,pdown, & 
      isglobal, ieglobal,ijkdisp, &
      jsglobal, jeglobal, noverlap, &
      ksglobal, keglobal,  &
      isgrid, iegrid, jsgrid, jegrid, ksgrid, kegrid

    integer :: jcomodd, jcomeven
    integer :: kcomodd, kcomeven
    integer, dimension(10) :: isubgrid

    public :: &
      spmd_init, & ! Initialize SPMD
      spmd_final, & ! Finalize SPMD
  !    calc_pipjpk, & ! Calculate co-ordinates of each processor
!      spmd_sendrecv, & ! MPI send-receive routine
      mpi_shorthand, & ! Shorthand stuff for MPI
!      spmd_comm_split_pi0, & ! New MPI communicator for myranks with pi = 0
      spmd_wtime ! Return wall-clock kime
!      spmd_file_write ! Write data on single file

!//////////////////////////////////////////////////


!------------------------------
contains
!------------------------------
!++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine spmd_init

!----------------------------------------
! DESCRIPTION:
!  Initialize SPMD
!----------------------------------------

!   use mpi
!   use shr_vars, only: myrank, nprocs, masterproc, wunitp

   implicit none
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   integer :: &
      ierr
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!
! Initialize MPI
!
   call mpi_init( ierr )
   if( ierr /= mpi_success ) then
      write(wunitp,*)"mpi_init failed, ierr = ",ierr," aborting"
      call abort
   end if

! 
! Get my id   
!
   call mpi_comm_rank( mpi_comm_world,myrank,ierr )
   if( ierr /= mpi_success ) then
      write(wunitp,*)"mpi_comm_myrank failed, ierr = ",ierr," aborting"
      call abort
   end if
   if( myrank == 0 ) then
      masterproc = .true.
   else
      masterproc = .false.
   end if

! 
! Get number of processors 
! 
   call mpi_comm_size( mpi_comm_world,nprocs,ierr )
   if( ierr /= mpi_success ) then
      write(wunitp,*)"mpi_comm_size failed, ierr = ",ierr," aborting"
      call abort
   end if

!
! Set some shorthand parameters
!
   call mpi_shorthand

   return

end subroutine spmd_init
!++++++++++++++++++++++++++++++++++++++++++++++++++


!++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine spmd_final

!----------------------------------------
! DESCRIPTION:
!  Finalize SPMD
!----------------------------------------

!   use mpi
!   use shr_vars, only: wunitp

   implicit none
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   integer ::&
      ierr
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   call mpi_finalize( ierr )

   if( ierr /= mpi_success ) then
      write(wunitp,*)"mpi_finalize failed, ierr = ",ierr," aborting"
      call abort
   end if

   return

end subroutine spmd_final
!++++++++++++++++++++++++++++++++++++++++++++++++++


!++++++++++++++++++++++++++++++++++++++++++++++++++
!subroutine calc_pipjpk
!shifted to read_grid.f90 where partition is happening

!!----------------------------------------
!! DESCRIPTION:
!!  Calculate co-ordinates (pi, pj, and pk) of each processor based on
!! processor myrank. This should be called only after namelist_cloud has
!!! been read, since np_x and np_y are read from namelist_cloud
!!----------------------------------------

!!   use shr_vars, only: myrank, np_x, np_y, pi, pj, pk

!   implicit none

!   pk = myrank / (np_x*np_y) ! z co-ordinates
!   pj = ( myrank - pk*(np_x*np_y) ) / np_x ! y co-ordinates
!   pi = myrank - pj*np_x - pk*(np_x*np_y) ! x co-ordinates

!   return

!end subroutine calc_pipjpk
!!++++++++++++++++++++++++++++++++++++++++++++++++++


!++++++++++++++++++++++++++++++++++++++++++++++++++
!subroutine spmd_sendrecv( axis,direction,dim1b,dim1e,dim2b,dim2e,fac, &
!                          arr2d_in,arr2d_out )

!----------------------------------------
! DESCRIPTION:
!  Subroutine to send and receive messages using mpi_send and
! mpi_receive
!----------------------------------------
! NOTES:
! Limited to 2-dimensional in and out arrays
! tag_* variables denote source processors
!----------------------------------------

!   use mpi
!   use shr_kind_mod, only: kind=8 => shr_kind_kind=8
!   use shr_vars, only: np_x, np_y, np_z, pi, pj, pk, mpicom, mpir8, zero

!   implicit none
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   integer, intent(in) :: &
!      dim1b, & ! Lower bound of 1st dimension of input array
!      dim1e, & ! Upper bound of 1st dimension of input array
!      dim2b, & ! Lower bound of 2nd dimension of input array
!      dim2e, & ! Upper bound of 2nd dimension of input array
!      fac ! Factor to multiply with tag
!   integer :: &
!      mpistat(mpi_status_size)
!   real(kind=8), intent(in) :: &
!      arr2d_in(dim1b:dim1e,dim2b:dim2e) ! Array to be sent
!   real(kind=8), intent(out) :: &
!      arr2d_out(dim1b:dim1e,dim2b:dim2e) ! Array to be received
!   character(len=1), intent(in) :: &
!      axis, & ! Axis of information flow
!      direction ! Direction of information flow 
!   integer :: &
!      len1d, & ! Length of 1st dimension of input array
!      len2d, & ! Length of 2nd dimension of input array
!      rp, & ! myrank of processor for send-receive communications
!      tag, & ! Tag associated with send-receive communications
!      ierr
!   real(kind=8), allocatable, dimension(:) :: &
!      arr1d_send, arr1d_recv ! Local send-receive arrays
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!!
!! Generate a random number
!!
!!   call random_number( rand )
!!   rand = rand * itr * 1e6
!!   fac = anint( rand )
!
!!
!! Allocate memory
!!
!   len1d = size( arr2d_in,1 ); len2d = size( arr2d_in,2 )
!   allocate( arr1d_send(len1d*len2d),arr1d_recv(len1d*len2d),stat=ierr )
!   if( ierr /= 0) &
!      call mpiabort("'spmd_sendrecv': unable to allocate memory to arr1d_send,arr1d_recv; aborting")
!
!!
!! Initialize to zero, *** do not change ***
!!
!   arr1d_send(:) = zero; arr1d_recv(:) = zero
!
!!
!! Feed 2d input array to 1d array for sending
!!
!   arr1d_send(:) = pack( arr2d_in(:,:),.true. )
!!   do j = dim2b,dim2e
!!      do i = dim1b,dim1e
!!         jj = (i-dim1b+1) + (j-dim2b)*len1d
!!         arr1d_send(jj) = arr2d_in(i,j)
!!      end do
!!   end do
!
!!
!! Calculate source and destination processors for each direction.
!! Do a blocking send and blocking receive simultaneously
!!
!   if( axis == 'x' ) then
!
!      if( pi /= np_x-1 ) then ! Use 0 to np_x-2 processors
!         rp = (pi+1) + pj*np_x + pk*np_x*np_y ! pi = 1 to np_x-1
!         if( direction == '+' ) then
!            tag = pi + pj*np_x + pk*np_x*np_y ! Source processors
!            call mpi_send( arr1d_send,len1d*len2d,mpir8,rp,tag*fac, &
!                           mpicom,ierr )
!         else if( direction == '-' ) then
!            tag = (pi+1) + pj*np_x + pk*np_x*np_y ! Source processors
!            call mpi_recv( arr1d_recv,len1d*len2d,mpir8,rp,tag*fac, &
!                           mpicom,mpistat,ierr )
!         end if
!      end if
!      if( pi /= 0 ) then ! Use 1 to np_x-1 processors
!         rp = (pi-1) + pj*np_x + pk*np_x*np_y ! pi = 0 to np_x-2
!         if( direction == '+' ) then
!            tag = (pi-1) + pj*np_x + pk*np_x*np_y ! Source processors
!            call mpi_recv( arr1d_recv,len1d*len2d,mpir8,rp,tag*fac, &
!                           mpicom,mpistat,ierr )
!         else if( direction == '-' ) then
!            tag = pi + pj*np_x + pk*np_x*np_y ! Source processors
!            call mpi_send( arr1d_send,len1d*len2d,mpir8,rp,tag*fac, &
!                           mpicom,ierr )
!         end if
!      end if
!
!   else if( axis == 'y' ) then
!
!      if( pj /= np_y-1 ) then ! Use 0 to np_y-2 processors
!         rp = pi + (pj+1)*np_x + pk*np_x*np_y ! pj = 1 to np_y-1
!         if( direction == '+' ) then
!            tag = pi + pj*np_x + pk*np_x*np_y ! Source processors
!            call mpi_send( arr1d_send,len1d*len2d,mpir8,rp,tag*fac, &
!                           mpicom,ierr )
!         else if( direction == '-' ) then
!            tag = pi + (pj+1)*np_x + pk*np_x*np_y ! Source processors
!            call mpi_recv( arr1d_recv,len1d*len2d,mpir8,rp,tag*fac, &
!                           mpicom,mpistat,ierr )
!         end if
!      end if
!      if( pj /= 0 ) then ! Use 1 to np_y-1 processors
!         rp = pi + (pj-1)*np_x + pk*np_x*np_y ! pj = 0 to np_y-2
!         if( direction == '+' ) then
!            tag = pi + (pj-1)*np_x + pk*np_x*np_y ! Source processors
!            call mpi_recv( arr1d_recv,len1d*len2d,mpir8,rp,tag*fac, &
!                           mpicom,mpistat,ierr )
!         else if( direction == '-' ) then
!            tag = pi + pj*np_x + pk*np_x*np_y ! Source processors
!            call mpi_send( arr1d_send,len1d*len2d,mpir8,rp,tag*fac, &
!                           mpicom,ierr )
!         end if
!      end if
!
!   else if( axis == 'z' ) then
!
!      if( pk /= np_z-1 ) then ! Use 0 to np_z-2 processors
!         rp = pi + pj*np_x + (pk+1)*np_x*np_y ! pk = 1 to np_z-1
!         if( direction == '+' ) then
!            tag = pi + pj*np_x + pk*np_x*np_y ! Source processors
!            call mpi_send( arr1d_send,len1d*len2d,mpir8,rp,tag*fac, &
!                           mpicom,ierr )
!         else if( direction == '-' ) then
!            tag = pi + pj*np_x + (pk+1)*np_x*np_y ! Source processors
!            call mpi_recv( arr1d_recv,len1d*len2d,mpir8,rp,tag*fac, &
!                           mpicom,mpistat,ierr )
!         end if
!      end if
!      if( pk /= 0 ) then ! Use 1 to np_z-1 processors
!         rp = pi + pj*np_x + (pk-1)*np_x*np_y ! pk = 0 to np_z-2
!         if( direction == '+' ) then
!            tag = pi + pj*np_x + (pk-1)*np_x*np_y ! Source processors
!            call mpi_recv( arr1d_recv,len1d*len2d,mpir8,rp,tag*fac, &
!                           mpicom,mpistat,ierr )
!         else if( direction == '-' ) then
!            tag = pi + pj*np_x + pk*np_x*np_y ! Source processors
!            call mpi_send( arr1d_send,len1d*len2d,mpir8,rp,tag*fac, &
!                           mpicom,ierr )
!         end if
!!      end if
!
!!   end if
!
!!
!! Feed received 1d array to 2d array to be sent out
!!
!!   arr2d_out(:,:) = reshape( arr1d_recv,[len1d,len2d] )
!!   do j = dim2b,dim2e
!!      do i = dim1b,dim1e
!!         jj = (i-dim1b+1) + (j-dim2b)*len1d
!!         arr2d_out(i,j) = arr1d_recv(jj)
!!      end do
!!   end do
!
!!   deallocate( arr1d_send,arr1d_recv )
!
!!   return
!
!!end subroutine spmd_sendrecv
!!++++++++++++++++++++++++++++++++++++++++++++++++++


!++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine mpi_shorthand

!----------------------------------------
! DESCRIPTION:
!  Define shorthand variables for MPI. Idea taken from CAM-3.1
!----------------------------------------

!   use mpi
!   use shr_vars, only: mpicom, mpinull, mpilog, mpiint, mpir8, &
!                       mpiundef, mpisum, mpimax

   implicit none
   save
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   integer :: &
      ierr
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!
! Need to set shorthands as variables rather than parameters since
! some MPI implementations set values for MPI tags at runtime
!
   call mpi_comm_dup( mpi_comm_world,mpicom,ierr )
   mpinull = mpi_info_null ! For mpi_info_null
   mpilog = mpi_logical ! For mpi_logical
   mpiint = mpi_integer ! For integer data type
   mpir8 = mpi_double ! For double/real*8 data type
   mpiundef = mpi_undefined ! For mpi_ndefined
   mpisum = mpi_sum ! For mpi_sum
   mpimax = mpi_max ! For mpi_max

   return

end subroutine mpi_shorthand
!++++++++++++++++++++++++++++++++++++++++++++++++++


!++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine spmd_comm_split_pi0

!----------------------------------------
! DESCRIPTION:
!  Define an MOI communicator for myranks with pi = 0
! Credits: Patricia Balle (CRAY Inc.)
!----------------------------------------

!   use mpi
!   use shr_vars, only: mpicom, mpicom_pi0, myrank, pi, mpiundef

   implicit none
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   integer :: &
      colour, & ! Colour for subgroups defined with and without myrank
                ! with pi = 0
      ierr
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   if( pi == 0 ) then
      colour = 1
   else
      colour = mpiundef
   end if

   call mpi_comm_split( mpicom,colour,myrank,mpicom_pi0,ierr )

   return

end subroutine spmd_comm_split_pi0
!++++++++++++++++++++++++++++++++++++++++++++++++++


!==================================================
real(kind=8) function spmd_wtime( w_time )

!----------------------------------------
! DESCRIPTION:
!  Return wall clock time for MPI jobs
!----------------------------------------   

   implicit none
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!#include<mpif.h>
   real(kind=8), intent(in) :: &
      w_time
   integer :: &
      ierr
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!   call mpi_barrier( ierr )
   spmd_wtime = mpi_wtime( ierr ) - w_time

   return

end function spmd_wtime
!==================================================


!++++++++++++++++++++++++++++++++++++++++++++++++++
!subroutine spmd_file_write( fld_name )

!----------------------------------------
! DESCRIPTION:
!  Write data on single file using MPI routines
!----------------------------------------

!   implicit none

!   character(len=*), intent(in) :: &
!      fld_name ! Name of field to be written
!   integer :: &
!      cmode, & ! file create mode
!      disp, & ! Displacement for each processor
!      ierr
!   character(cs) :: &
!      fl_name, & ! Name of data file
!      fmt_str
!//////////////////////////////////////////////////

!   call endrun("'spmd_file_write' not implemented; aborting")

!   fmt_str = '(a,a,a1,a,a4,i0)'
!   write(fl_name,fmt_str)'output/',trim(caseid),'_',fld_name,'_itr',itr
!   fl_name = trim(fl_name)//'.bin'

!
! Open file
!
!   cmode = mpi_mode_wronly + mpi_mode_create
!   call mpi_file_open( mpicom,trim(fl_name),cmode,mpinull,mpistat,ierr )
!
!   disp = (1 + pi*(nx-1)) * (1 + pj*(ny-1)) * (1 + pk*(nz-1))
!   call mpi_file_set_view( mpistat,disp,mpir8,mpiint,'native', & 
!                          mpinull,ierr)
!   call mpi_file_write( mpistat, buf, BUFSIZE, MPI_INTEGER, & 
!                        MPI_STATUS_IGNORE, ierr)

!
! Close the file
!
!   call mpi_file_close( mpistat,ierr )

!   return

!end subroutine spmd_file_write
!++++++++++++++++++++++++++++++++++++++++++++++++++


end module mod_spmd
