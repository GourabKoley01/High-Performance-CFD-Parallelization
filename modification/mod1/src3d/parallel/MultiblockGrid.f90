subroutine MultiblockGrid
  use mod_params
  use mod_spmd
  implicit none

  integer::  iblk, icount,i,j,k 
  integer:: ibs,ibe,jbs,jbe, kbs, kbe
  integer:: ibs1,ibe1,jbs1,jbe1, kbs1, kbe1 
  !  nblnksg = nblnks

  if(myrank .eq. 0)then
    print*, 'Its multiblock problem with nblnks= ', nblnksg  
    !print*, 'Its multiblock problem with nblnks= ,' nblnksg  
    do iblk=1, nblnksg
      print*,'block = ', iblk
      print*, isblkg(iblk), ieblkg(iblk)
      print*, jsblkg(iblk), jeblkg(iblk)
      print*, ksblkg(iblk), keblkg(iblk)
    enddo
  endif

  nblnks = 0

  do iblk=1,nblnksg
    ibs = isblkg(iblk)
    ibe = ieblkg(iblk)
    jbs = jsblkg(iblk)
    jbe = jeblkg(iblk)
    kbs = ksblkg(iblk)
    kbe = keblkg(iblk)

    if(((ibs .ge. isglobal) .and. (ibs .le. ieglobal)) &
      & .or.((ibe .ge. isglobal) .and. (ibe .le. ieglobal))&
      & .or.((ibs .le. isglobal) .and. (ibe .ge. ieglobal)))then

    if(((jbs .ge. jsglobal) .and. (jbs .le. jeglobal)) &
      & .or.((jbe .ge. jsglobal) .and. (jbe .le. jeglobal))&
      & .or.((jbs .le. jsglobal) .and. (jbe .ge. jeglobal)))then

    if(((kbs .ge. ksglobal) .and. (kbs .le. keglobal)) &
      & .or.((kbe .ge. ksglobal) .and. (kbe .le. keglobal))&
      & .or.((kbs .le. ksglobal) .and. (kbe .ge. keglobal)))then
    nblnks = nblnks+1
  endif
endif
endif
enddo

if(nblnks .gt. 0)then
  allocate(isblk(nblnks))
  allocate(ieblk(nblnks))
  allocate(jsblk(nblnks))
  allocate(jeblk(nblnks))
  allocate(ksblk(nblnks))
  allocate(keblk(nblnks))
endif

print*,'------------------------------------------------------------'
print*, 'myrank = ', myrank,'nx,ny,nz = ', nx, ny,nz, 'nblnks = ', nblnks
print*, 'global=', isglobal, ieglobal, jsglobal, jeglobal, ksglobal, keglobal
print*, 'local =', isgrid, iegrid, jsgrid, jegrid, ksgrid, kegrid
if(nblnks .gt. 0)then
  !print*, 'local ij = 1 ', nx, ' 1 ', ny
  icount=0;
  do iblk=1,nblnksg
    ibs = isblkg(iblk)
    ibe = ieblkg(iblk)
    jbs = jsblkg(iblk)
    jbe = jeblkg(iblk)
    kbs = ksblkg(iblk)
    kbe = keblkg(iblk)

    if(((ibs .ge. isglobal) .and. (ibs .le. ieglobal)) &
      & .or.((ibe .ge. isglobal) .and. (ibe .le. ieglobal))&
      & .or.((ibs .le. isglobal) .and. (ibe .ge. ieglobal)))then

    if(((jbs .ge. jsglobal) .and. (jbs .le. jeglobal)) &
      & .or.((jbe .ge. jsglobal) .and. (jbe .le. jeglobal))&
      & .or.((jbs .le. jsglobal) .and. (jbe .ge. jeglobal)))then

    if(((kbs .ge. ksglobal) .and. (kbs .le. keglobal)) &
      & .or.((kbe .ge. ksglobal) .and. (kbe .le. keglobal))&
      & .or.((kbs .le. ksglobal) .and. (kbe .ge. keglobal)))then

    print*,'Global Block # ',iblk
    !print*,'global= ', isblkg(iblk), ieblkg(iblk), jsblkg(iblk), jeblkg(iblk)
    print*,'global= ', max(isblkg(iblk),isglobal), min(ieblkg(iblk),ieglobal),&
      & max(jsblkg(iblk),jsglobal), min(jeblkg(iblk),jeglobal)
    icount = icount+1
    !ibs1 = (ibs - isglobal)+1 
    !ibe1 = (ibe - isglobal)+1
    !jbs1 = (jbs - jsglobal)+1
    !jbe1 = (jbe - jsglobal)+1
    ibs1 = (ibs - isglobal)+isgrid 
    ibe1 = (ibe - isglobal)+isgrid
    jbs1 = (jbs - jsglobal)+jsgrid
    jbe1 = (jbe - jsglobal)+jsgrid
    kbs1 = (kbs - ksglobal)+ksgrid
    kbe1 = (kbe - ksglobal)+ksgrid

    isblk(icount) = max(ibs1,isgrid)
    ieblk(icount) = min(ibe1,iegrid)
    jsblk(icount) = max(jbs1,jsgrid)
    jeblk(icount) = min(jbe1,jegrid)
    ksblk(icount) = max(kbs1,ksgrid)
    keblk(icount) = min(kbe1,kegrid)

    !take care of the overlap region
    if(isblk(icount) .eq. isgrid)isblk(icount) = 1
    if(ieblk(icount) .eq. iegrid)ieblk(icount) = nx
    if(jsblk(icount) .eq. jsgrid)jsblk(icount) = 1
    if(jeblk(icount) .eq. jegrid)jeblk(icount) = ny
    if(ksblk(icount) .eq. ksgrid)ksblk(icount) = 1
    if(keblk(icount) .eq. kegrid)keblk(icount) = nz
    !isblk(icount) = max(ibs1,1)
    !ieblk(icount) = min(ibe1,nx)
    !jsblk(icount) = max(jbs1,1)
    !jeblk(icount) = min(jbe1,ny)

    !isblk(nblnks) = ibs1 - noverlap
    !ieblk(nblnks) = ibe1 - noverlap
    !jsblk(nblnks) = jbs1 - noverlap
    !jeblk(nblnks) = jbe1 - noverlap
  endif
endif
endif
  enddo

  do iblk=1,nblnks
    !print*,'block= ', iblk
    print*,'Local Block # ',iblk
    print*,'local= ', isblk(iblk), ieblk(iblk), jsblk(iblk), jeblk(iblk),&
      ksblk(iblk), keblk(iblk)
  enddo

  do iblk=1,nblnks
    ibs = isblk(iblk)
    ibe = ieblk(iblk)
    jbs = jsblk(iblk)
    jbe = jeblk(iblk)
    kbs = ksblk(iblk)
    kbe = keblk(iblk)
    do i=ibs,ibe
      do j=jbs,jbe
        do k=kbs,kbe
          ifblock(i,j,k) = 1
        enddo
      enddo
    enddo
  enddo
endif
print*,'------------------------------------------------------------'


!if(nblnks .ne. 0)then
!do iblk=1,nblnksg
!  print*,'block= ',iblk, isblkg(nblnks), ieblkg(nblnks), jsblkg(nblnks), jeblkg(nblnks)
!enddo
!note here nblnks and nblnksg may not be the same
!endif

!print*, myrank, isglobal, ieglobal, nblnks
!print*, 'myrank, isglobal, ieglobal, nblanks'
!print*, myrank, isglobal, ieglobal, nblanks

! & .or. (ibe .le. ieglobal) &
!   & .or. (jbs .ge. jsglobal) .or. 
!  do i = isglobal,ieglobal
!    if((ibs .ge. isglobal)   
end subroutine MultiblockGrid
