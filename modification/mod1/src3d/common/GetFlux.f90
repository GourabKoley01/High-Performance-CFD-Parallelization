subroutine GetFlux
  use mod_params
  use nvtx
  implicit none
  integer:: i,j,ivar, iblk
  integer:: ibs,ibe,jbs, jbe

  !$acc kernels !!async(1)
  rhs=0.0
  !$acc end kernels

  !tdg=0.0
  
  !call ApplyBC
  ! Above Apply BC call for Multiblock --  
  !print*,'1',maxval(rhs)
  !print*,'bc applied'
  !call compact_flux
  !print*,'2',maxval(rhs)
  !print*,'compact flux done'
  if (iflux_scheme .eq. 10) then
    call nvtxStartRange("CompFlux")
    call CompactFlux
    call nvtxEndRange
  else if (iflux_scheme .eq. 2) then
    call nvtxStartRange("RoeFlux")
    call RoeFlux
    call nvtxEndRange
  endif

  if (ivis_scheme .eq. 10) then
    call CompactVflux
  else if (ivis_scheme .eq. 2) then
    call nvtxStartRange("SecondVflux")
    call SecondVflux
    call nvtxEndRange
  else if (ivis_scheme .eq. 1) then
    call FirstVflux
  endif


  !rhs = -dt*rhs !this operation is done in RK3 update

  !if (nblnks .ne. 0)then
  !  do iblk=1,nblnks
  !    ibs = isblk(iblk)
  !    ibe = ieblk(iblk)
  !    jbs = jsblk(iblk)
  !    jbe = jeblk(iblk)
  !    do ivar =1,nvars
  !      do i=ibs,ibe
  !        do j=jbs,jbe
  !          rhs(i,j,ivar) = 0.0
  !        enddo
  !      enddo
  !    enddo
  !  enddo
  !endif
  !if(iflux_scheme .eq. 10)call filter(itr)

  !print*,'3',maxval(rhs)
end subroutine GetFlux
