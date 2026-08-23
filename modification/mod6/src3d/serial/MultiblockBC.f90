subroutine MultiblockBC
  use mod_params
  implicit none
  integer :: i,j,k
  integer :: iblk,ibs,ibe,jbs,jbe,kbs,kbe, ivar

  !block is the blanked region and block boundaries are walls
  !make sure block boundary does not fall on processor boundary
  
  do iblk=1,nblnksg
    ibs = isblkg(iblk)
    ibe = ieblkg(iblk)
    jbs = jsblkg(iblk)
    jbe = jeblkg(iblk)
    ! Addition for 3D BC ----> By Mritunjay [May-24-2024]
    kbs = ksblkg(iblk)
    kbe = keblkg(iblk)
    !$acc parallel present(prim)
    if(ibs .ne. 1)then
      !$acc loop gang collapse(2) 
      do j=jbs,jbe
        do k = kbs, kbe
        prim(ibs,j,k,1) = prim(ibs-1,j,k,1)
        prim(ibs,j,k,2) = 0.0
        prim(ibs,j,k,3) = 0.0
        prim(ibs,j,k,4) = 0.0
        prim(ibs,j,k,5) = prim(ibs-1,j,k,5)
        end do
      enddo
    endif
    if(ibe .ne. nx)then
      !$acc loop gang collapse(2) 
      do j=jbs,jbe
        do k=kbs,kbe
        prim(ibe,j,k,1) = prim(ibe+1,j,k,1)
        prim(ibe,j,k,2) = 0.0
        prim(ibe,j,k,3) = 0.0
        prim(ibe,j,k,4) = 0.0
        prim(ibe,j,k,5) = prim(ibe+1,j,k,5)
        end do
      enddo
    endif

    if(jbs .ne. 1)then
      !$acc loop gang collapse(2) 
      do i=ibs,ibe
        do k=kbs, kbe
        prim(i,jbs,k,1) = prim(i,jbs-1,k,1)
        prim(i,jbs,k,2) = 0.0
        prim(i,jbs,k,3) = 0.0
        prim(i,jbs,k,4) = 0.0
        prim(i,jbs,k,5) = prim(i,jbs-1,k,5)
        end do
      enddo
    endif
    if(jbe .ne. ny)then
      !$acc loop gang collapse(2) 
      do i=ibs,ibe
        do k=kbs, kbe
        prim(i,jbe,k,1) = prim(i,jbe+1,k,1)
        prim(i,jbe,k,2) = 0.0
        prim(i,jbe,k,3) = 0.0
        prim(i,jbe,k,4) = 0.0
        prim(i,jbe,k,5) = prim(i,jbe+1,k,5)
        end do
      enddo
    endif

    ! Addition of wall in z-direction block
    if(kbs .ne. 1)then
      !$acc loop gang collapse(2) 
      do i=ibs,ibe
        do j=jbs,jbe
           prim(i,j,kbs,1) = prim(i,j,kbs-1,1)
           prim(i,j,kbs,2) = 0.0
           prim(i,j,kbs,3) = 0.0
           prim(i,j,kbs,4) = 0.0
           prim(i,j,kbs,5) = prim(i,j,kbs-1,5)
        enddo
      enddo
    endif
    if(kbe .ne. nz)then
      !$acc loop gang collapse(2) 
      do i=ibs,ibe
        do j=jbs,jbe
           prim(i,j,kbe,1) = prim(i,j,kbe+1,1)
           prim(i,j,kbe,2) = 0.0
           prim(i,j,kbe,3) = 0.0
           prim(i,j,kbe,4) = 0.0
           prim(i,j,kbe,5) = prim(i,j,kbe+1,5)
        enddo
      enddo
    endif
    !$acc end parallel
  enddo


  !put something inside !internal of the block
  do iblk=1,nblnksg
    ibs = isblkg(iblk)
    ibe = ieblkg(iblk)
    jbs = jsblkg(iblk)
    jbe = jeblkg(iblk)
    ! Addition of 3D BC ---> [May,24-2024]
    kbs = ksblkg(iblk)
    kbe = keblkg(iblk)
    !$acc parallel loop gang collapse(3) present(prim)
    do i=ibs+1,ibe-1
      do j=jbs+1,jbe-1
        do k=kbs+1,kbe-1
        prim(i,j,k,1) = 1.0
        prim(i,j,k,2) = 0.0
        prim(i,j,k,3) = 0.0
        prim(i,j,k,4) = 0.0
        prim(i,j,k,5) = pinf
        enddo
      enddo
    enddo
  enddo


  !do we need to make gradients zero next to the boundary, may be needed for Roe
  do iblk=1,nblnksg
    ibs = isblkg(iblk)
    ibe = ieblkg(iblk)
    jbs = jsblkg(iblk)
    jbe = jeblkg(iblk)
    ! Addition of 3D BC --> [May-24-2024]
    kbs = ksblkg(iblk)
    kbe = keblkg(iblk)
    !$acc parallel present(prim)

    ! BC in x-direction
    if(ibs .ne. 1)then
      !i = min(ibs+1,nx)
      i = ibs+1
      !$acc loop gang collapse(3)
      do j=jbs,jbe
        do k=kbs,kbe            ! <-----
          do ivar=1,5
            prim(i,j,k,ivar) = prim(i-1,j,k,ivar)
          enddo
        enddo
      enddo
    else
      i=1
      !$acc loop gang collapse(2)
      do j=jbs,jbe
        do k=kbs,kbe             ! <-------
          prim(i,j,k,1) = 1.0
          prim(i,j,k,2) = 0.0
          prim(i,j,k,3) = 0.0
          prim(i,j,k,4) = 0.0
          prim(i,j,k,5) = pinf
        enddo
      enddo
    endif

    if(ibe .ne. nx)then
      !i = max(ibe-1,1)
      i = ibe-1
      !$acc loop gang collapse(3)
      do j=jbs,jbe
        do k=kbs,kbe
          do ivar=1,5
            prim(i,j,k,ivar) = prim(i+1,j,k,ivar)
          enddo
        enddo
      enddo
    else
      i=ibe
      !$acc loop gang collapse(2)
      do j=jbs,jbe
        do k=kbs,kbe       ! <------
          prim(i,j,k,1) = 1.0
          prim(i,j,k,2) = 0.0
          prim(i,j,k,3) = 0.0
          prim(i,j,k,4) = 0.0
          prim(i,j,k,5) = pinf
        enddo
      enddo
    endif

    ! BC for y-direction
    if(jbs .ne. 1)then
      !j = min(jbs+1,ny)
      j = jbs+1
      !$acc loop gang collapse(3)
      do i=ibs,ibe
        do k=kbs,kbe             ! <-----
          do ivar=1,5
            prim(i,j,k,ivar) = prim(i,j-1,k,ivar)
          enddo
        enddo
      enddo
    else
      j=1
      !$acc loop gang collapse(2)
      do i=ibs,ibe
        do k=kbs,kbe              ! <----
          prim(i,j,k,1) = 1.0
          prim(i,j,k,2) = 0.0
          prim(i,j,k,3) = 0.0
          prim(i,j,k,4) = 0.0
          prim(i,j,k,5) = pinf
        enddo
      enddo
    endif

    if(jbe .ne. ny)then
      !j = max(jbe-1,1)
      j=jbe-1
      !$acc loop gang collapse(3)
      do i=ibs,ibe
        do k=kbs,kbe             ! <----
          do ivar=1,5
            prim(i,j,k,ivar) = prim(i,j+1,k,ivar)
          enddo
        enddo
      enddo
    else
      j=jbe
      !$acc loop gang collapse(2)
      do i=ibs,ibe
        do k=kbs,kbe 
          prim(i,j,k,1) = 1.0
          prim(i,j,k,2) = 0.0
          prim(i,j,k,3) = 0.0
          prim(i,j,k,4) = 0.0
          prim(i,j,k,5) = pinf
        enddo
      enddo 
    endif

    ! Addition of BC for z-direction
    if(kbs .ne. 1)then
      !k = min(kbs+1,nz)
      k=kbs+1
      !$acc loop gang collapse(3)
      do i=ibs,ibe
        do j=jbs,jbe             ! <-----
          do ivar=1,5
            prim(i,j,k,ivar) = prim(i,j,k-1,ivar)
          enddo
        enddo
      enddo
    else
      k=1
      !$acc loop gang collapse(2)
      do i=ibs,ibe
        do j=jbs,jbe              ! <----
          prim(i,j,k,1) = 1.0
          prim(i,j,k,2) = 0.0
          prim(i,j,k,3) = 0.0
          prim(i,j,k,4) = 0.0
          prim(i,j,k,5) = pinf
        enddo
      enddo
    endif

    if(kbe .ne. nz)then
      !k = max(kbe-1,1)
      k=kbe-1
      !$acc loop gang collapse(3)
      do i=ibs,ibe
        do j=jbs,jbe             ! <----
          do ivar=1,5
            prim(i,j,k,ivar) = prim(i,j,k+1,ivar)
          enddo
        enddo
      enddo
    else
      k=kbe
      !$acc loop gang collapse(2)
      do i=ibs,ibe
        do j=jbs,jbe 
          prim(i,j,k,1) = 1.0
          prim(i,j,k,2) = 0.0
          prim(i,j,k,3) = 0.0
          prim(i,j,k,4) = 0.0
          prim(i,j,k,5) = pinf
        enddo
      enddo 
    endif
    !$acc end parallel
  enddo

  !Redo inside
  do iblk=1,nblnksg
    ibs = isblkg(iblk)
    ibe = ieblkg(iblk)
    jbs = jsblkg(iblk)
    jbe = jeblkg(iblk)
    ! Addition of 3D BC ----> [May-24-2024]
    kbs = ksblkg(iblk)
    kbe = keblkg(iblk)
    !$acc parallel loop gang collapse(3) present(prim)
    do i=ibs+1,ibe-1
      do j=jbs+1,jbe-1
        do k=kbs+1,kbe-1
          prim(i,j,k,1) = 1.0
          prim(i,j,k,2) = 0.0
          prim(i,j,k,3) = 0.0
          prim(i,j,k,4) = 0.0
          prim(i,j,k,5) = pinf
        enddo
      enddo
    enddo
  enddo

  !prim to cons
  do iblk=1,nblnksg
    ibs = isblkg(iblk)
    ibe = ieblkg(iblk)
    jbs = jsblkg(iblk)
    jbe = jeblkg(iblk)
    ! Addition of 3D BC
    kbs = ksblkg(iblk)
    kbe = keblkg(iblk)
    !$acc parallel loop gang collapse(3) present(prim, q)
    do i=ibs,ibe
      do j=jbs,jbe
        do k=kbs,kbe
          q(i,j,k,1) = prim(i,j,k,1)
          q(i,j,k,2) = prim(i,j,k,2)*prim(i,j,k,1) 
          q(i,j,k,3) = prim(i,j,k,3)*prim(i,j,k,1)
          q(i,j,k,4) = prim(i,j,k,4)*prim(i,j,k,1)
          q(i,j,k,5) = GAM1I * prim(i,j,k,5)  + 0.50d0*prim(i,j,k,1)*&
            & (prim(i,j,k,2)**2+prim(i,j,k,3)**2+prim(i,j,k,4)**2)
        enddo
      enddo
    enddo
  enddo

end subroutine MultiblockBC
