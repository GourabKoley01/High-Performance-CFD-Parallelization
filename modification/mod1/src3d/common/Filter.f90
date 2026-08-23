subroutine filter(itr)
  use mod_params
  use nvtx
  implicit none
  integer, intent(in):: itr

  !print*, itr, maxval(q(:,:,1)), minval(q(:,:,1)), q(81,81,1)
  !print*, itr, maxval(q(4:nx-4,4:ny-4,1)), minval(q(4:nx-4,4:ny-4,1)), q(81,81,1)
  
  call nvtxStartRange("xifilter")
  call xifilter
  call nvtxEndRange
  call nvtxStartRange("etfilter")
  call etfilter
  call nvtxEndRange
  call nvtxStartRange("ztfilter")
  call ztfilter
  call nvtxEndRange
  
  !if(mod(itr,2) .eq. 0)then
  !  call ifilter
  !  call jfilter
  !else
  !  call jfilter
  !  call ifilter
  !endif
  !print*, itr, maxval(q(:,:,1)), minval(q(:,:,1)), q(81,81,1)
  !print*, itr, maxval(q(2:nx-1,2:ny-1,1)), minval(q(2:nx-1,2:ny-1,1)), q(81,81,1)
  !print*, itr, maxval(q(4:nx-4,4:ny-4,1)), minval(q(4:nx-4,4:ny-4,1)), q(81,81,1)
  !print*, itr, maxval(q(2:nx-1,2:ny-1,1)), minval(q(2:nx-1,2:ny-1,1)), q(81,81,1)

  ! print*, itr, maxval(q(:,:,1)), minval(q(:,:,1)), q(81,81,1)

  !allocate(temp(nx,ny,nvars))
  !temp = q
  !q=temp
  !do ivar =1,nvars
  !  do j = 2,ny-1
  !    do i = 2,nx-1
  !      q(i,j,ivar) = temp(i,j,ivar)
  !    enddo
  !  enddo
  !enddo
  !deallocate(temp)

contains
  subroutine xifilter
    integer :: i,j,k,ifil,ivar,ifilm,nfilt 
   
    !$acc parallel loop gang collapse(3)
    do k = 2,nz-1
      do j = 2,ny-1
        do i = 1,nx
          if(i .lt. 5)then
            nfilt = i
          elseif(i .gt. nx-4)then
            nfilt = nx - i +1
          else
            nfilt = ifilt_order/2 +1
          endif

          atdg_field(i,j,k) = alphaf(nfilt)
          btdg_field(i,j,k) = 1.0d0
          ctdg_field(i,j,k) = alphaf(nfilt)

          do ivar =1,nvars
            rhstdg_field(i,j,k,ivar) = 0.0d0
            do ifil = 1,nfilt
              ifilm = ifil-1
              rhstdg_field(i,j,k,ivar) = rhstdg_field(i,j,k,ivar) +&
                & afilt(nfilt,ifil)*(q(i+ifilm,j,k,ivar)+q(i-ifilm,j,k,ivar))
            enddo
          enddo
        enddo
      enddo
    enddo

    call ThomasAlg(nvars,'xi')

    !$acc parallel loop gang collapse(4)
    do ivar =1,nvars
      do k = 2,nz-1
        do j = 2,ny-1
          do i = 1,nx
            q(i,j,k,ivar) = temp(i,j,k,ivar)
          enddo
        enddo
      enddo
    enddo
  end subroutine xifilter
  
  subroutine xifilter_old
    integer :: i,j,k,ifil,ivar,ifilm,nfilt 
    do j = 2,ny-1
      do k = 2,nz-1
        do i = 1,nx
          if(i .lt. 5)then
            nfilt = i
          elseif(i .gt. nx-4)then
            nfilt = nx - i +1
          else
            nfilt = ifilt_order/2 +1
          endif

          atdg(k,i) = alphaf(nfilt)
          btdg(k,i) = 1.0d0
          ctdg(k,i) = alphaf(nfilt)

          do ivar =1,nvars
            rhstdg(k,i,ivar) = 0.0d0
            do ifil = 1,nfilt
              ifilm = ifil-1
              rhstdg(k,i,ivar) = rhstdg(k,i,ivar) +&
                & afilt(nfilt,ifil)*(q(i+ifilm,j,k,ivar)+q(i-ifilm,j,k,ivar))
            enddo
          enddo
        enddo
      enddo

      call ThomasAlg(nvars,nx, 2, nz-1)

      do k = 2,nz-1
        do i = 1,nx
          do ivar =1,nvars
            q(i,j,k,ivar) = tdg(k,i,ivar)
          enddo
        enddo
      enddo
    enddo
  end subroutine xifilter_old

  subroutine etfilter
    integer :: i,j,k,ifil,ivar,ifilm,nfilt,joff

    !$acc parallel loop gang collapse(3)
    do k = 2,nz-1
      do j = 1,ny 
        do i = 2,nx-1
          if(j .lt. 5)then
            nfilt = j
          elseif(j .gt. ny-4)then
            nfilt = ny - j +1
          else
            nfilt = ifilt_order/2 +1
          endif
          atdg_field(i,j,k) = alphaf(nfilt)
          btdg_field(i,j,k) = 1.0d0
          ctdg_field(i,j,k) = alphaf(nfilt)
          !if((j .eq. 1) .or. (j .eq. ny))then
          !  atdg(j) = 0.0d0
          !  ctdg(j) = 0.0d0
          !endif
          do ivar =1,nvars
            rhstdg_field(i,j,k,ivar) = 0.0d0
            do ifil = 1,nfilt
              ifilm = ifil-1
              rhstdg_field(i,j,k,ivar) = rhstdg_field(i,j,k,ivar) +&
                & afilt(nfilt,ifil)*(q(i,j+ifilm,k,ivar)+q(i,j-ifilm,k,ivar))
            enddo
          enddo
        enddo
      enddo
    enddo

    call ThomasAlg(nvars,'eta')

    !$acc parallel loop gang collapse(4)
    do ivar =1,nvars
      do k = 2,nz-1
        do j = 1,ny
          do i = 2,nx-1
            q(i,j,k,ivar) = temp(i,j,k,ivar)
          enddo
        enddo
      enddo
    enddo
  end subroutine etfilter

  subroutine ztfilter
    integer :: i,j,k,ifil,ivar,ifilm,nfilt,joff

    !$acc parallel loop gang collapse(3)
    do k = 1,nz
      do j = 2,ny-1
        do i = 2,nx-1
          if(k .lt. 5)then
            nfilt = k
          elseif(k .gt. nz-4)then
            nfilt = nz - k +1
          else
            nfilt = ifilt_order/2 +1
          endif
          atdg_field(i,j,k) = alphaf(nfilt)
          btdg_field(i,j,k) = 1.0d0
          ctdg_field(i,j,k) = alphaf(nfilt)
          !if((j .eq. 1) .or. (j .eq. ny))then
          !  atdg(j) = 0.0d0
          !  ctdg(j) = 0.0d0
          !endif
          do ivar =1,nvars
            rhstdg_field(i,j,k,ivar) = 0.0d0
            do ifil = 1,nfilt
              ifilm = ifil-1
              rhstdg_field(i,j,k,ivar) = rhstdg_field(i,j,k,ivar) +&
                & afilt(nfilt,ifil)*(q(i,j,k+ifilm,ivar)+q(i,j,k-ifilm,ivar))
            enddo
          enddo
        enddo
      enddo
    enddo

    call ThomasAlg(nvars,'zeta')

    !$acc parallel loop gang collapse(4)
    do ivar =1,nvars
      do k = 1,nz
        do j = 2,ny-1
          do i = 2,nx-1
            q(i,j,k,ivar) = temp(i,j,k,ivar)
          enddo
        enddo
      enddo
    enddo 
  end subroutine ztfilter

end subroutine filter
