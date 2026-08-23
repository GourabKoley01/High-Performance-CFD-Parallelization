subroutine printscreen(itr) 
  use mod_spmd
  use mod_params
  implicit none
  integer, intent(in):: itr
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  integer :: i,j,k, ierr
  real(kind=8)::rhomin, umin, pmin, resnorm
  real(kind=8)::uming, pming, resnormg
  !  integer :: i,j,k,itr, ierr, ntot

  umin=minval(prim(isgrid:iegrid,jsgrid:jegrid,ksgrid:kegrid,2))
  Pmin=minval(prim(isgrid:iegrid,jsgrid:jegrid,ksgrid:kegrid,5))
  !print*,'Pmin=',myrank,Pmin
  resnorm = 0.0d0
  do k=ksgrid,kegrid
  do j=jsgrid,jegrid
    do i=isgrid,iegrid
      resnorm = resnorm + rhs(i,j,k,1)*rhs(i,j,k,1)
    enddo
  enddo
enddo
  !call MPI_REDUCE(rhomin,rhoming,1,mpir8,MPI_MIN,0,mpicom,ierr) 
  call MPI_REDUCE(umin,uming,1,mpir8,MPI_MIN,0,mpicom,ierr) 
  call MPI_REDUCE(Pmin,Pming,1,mpir8,MPI_MIN,0,mpicom,ierr) 
  call MPI_REDUCE(resnorm,resnormg,1,mpir8,MPI_SUM,0,mpicom,ierr) 
  resnormg = sqrt(resnormg)
  if(masterproc)then
    print 100, itr,itr*dt, uming, Pming, log10(resnormg+eps)
    write(wunitp,100)itr,itr*dt, uming, Pming, log10(resnormg+eps)

    100 format(I7.1,4F14.7)
    !100 format(I7.1,2F11.4,2F11.4)
    !100 format(I5.1,2F11.4,2F11.4)
  endif

  rhomin=minval(prim(isgrid:iegrid,jsgrid:jegrid,ksgrid:kegrid,1))
  if(rhomin .lt. 0)then
    print*,'Density going negative!!', myrank, rhomin
    stop
  endif
  !
end subroutine
