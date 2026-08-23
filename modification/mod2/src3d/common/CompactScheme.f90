subroutine CompactDeriv(nvartdg,dir)
    use mod_params
    use nvtx

    integer, intent(in) :: nvartdg
    character(len=*), intent(in) :: dir 
    integer :: i,j,k,ivar
    call nvtxStartRange("CompactDeriv")

    !$acc parallel loop gang collapse(3) 
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          atdg_field(i,j,k) = C6(1)
          btdg_field(i,j,k) = 1.0d0
          ctdg_field(i,j,k) = C6(1)
        enddo
      enddo             
    enddo

    if (dir .eq. "xi") then
      !$acc parallel loop gang collapse(4) 
      do ivar = 1,nvartdg
        do k = 2,nz-1
          do j = 2,ny-1
            do i = 3,nx-2
              rhstdg_field(i,j,k,ivar) = 0.5*C6(2)*( temp(i+1,j,k,ivar)-temp(i-1,j,k,ivar) ) + &
                             & 0.25*C6(3)*( temp(i+2,j,k,ivar)-temp(i-2,j,k,ivar) )
            enddo
          enddo
        enddo             
      enddo
    
    elseif (dir .eq. "eta") then
      !$acc parallel loop gang collapse(4) 
      do ivar = 1,nvartdg
        do k = 2,nz-1
          do j = 3,ny-2
            do i = 2,nx-1
              rhstdg_field(i,j,k,ivar) = 0.5*C6(2)*( temp(i,j+1,k,ivar)-temp(i,j-1,k,ivar) ) + &
                             & 0.25*C6(3)*( temp(i,j+2,k,ivar)-temp(i,j-2,k,ivar) )
            enddo
          enddo             
        enddo
      enddo
    
    elseif (dir .eq. "zeta") then
      !$acc parallel loop gang collapse(4) 
      do ivar = 1,nvartdg
        do k = 3,nz-2
          do j = 2,ny-1
            do i = 2,nx-1
              rhstdg_field(i,j,k,ivar) = 0.5*C6(2)*( temp(i,j,k+1,ivar)-temp(i,j,k-1,ivar) ) + &
                             & 0.25*C6(3)*( temp(i,j,k+2,ivar)-temp(i,j,k-2,ivar) )
            enddo
          enddo             
        enddo
      enddo

    endif

    if (iper .eq. 0) then
      call DerivBound(nvartdg,dir)
    endif

    call nvtxEndRange
  end subroutine CompactDeriv


  subroutine DerivBound(nvartdg,dir)
    use mod_params

    integer, intent(in) :: nvartdg
    character(len=*), intent(in) :: dir
    integer :: i,j,k,ivar

    if (dir .eq. "xi") then
      !$acc parallel loop gang collapse(2)
      do k = 2,nz-1
        do j = 2,ny-1
          atdg_field(1,j,k) = 0.0d0; btdg_field(1,j,k) = 1.0d0; ctdg_field(1,j,k) = C3B(1)
          atdg_field(nx,j,k) = C3B(1); btdg_field(nx,j,k) = 1.0d0; ctdg_field(nx,j,k) = 0.0d0

          atdg_field(2,j,k) = AC4B(1); btdg_field(2,j,k) = 1.0d0; ctdg_field(2,j,k) = AC4B(1)
          atdg_field(nx-1,j,k) = AC4B(1); btdg_field(nx-1,j,k) = 1.0d0; ctdg_field(nx-1,j,k) = AC4B(1)
        enddo
      enddo
      !$acc parallel loop gang collapse(3)
      do ivar = 1,nvartdg
        do k = 2,nz-1
          do j = 2,ny-1
            rhstdg_field(1,j,k,ivar) = C3B(2)*temp(1,j,k,ivar) + C3B(3)*temp(2,j,k,ivar) + &
                        & C3B(4)*temp(3,j,k,ivar)
            rhstdg_field(nx,j,k,ivar) = -C3B(2)*temp(nx,j,k,ivar) - C3B(3)*temp(nx-1,j,k,ivar) - &
                          & C3B(4)*temp(nx-2,j,k,ivar)
            !rhstdg_field(2,j,k,ivar) = 0.50d0*AC4B(2)*(temp(3,j,k,ivar) - temp(1,j,k,ivar))                
            !rhstdg_field(nx-1,j,k,ivar) = 0.50d0*AC4B(2)*(temp(nx,j,k,ivar) - temp(nx-2,j,k,ivar))
            rhstdg_field(2,j,k,ivar) = AC4B(2)*temp(1,j,k,ivar) + AC4B(4)*temp(3,j,k,ivar)
            rhstdg_field(nx-1,j,k,ivar) = -AC4B(2)*temp(nx,j,k,ivar) - AC4B(4)*temp(nx-2,j,k,ivar)
          enddo
        enddo
      enddo

    elseif (dir .eq. "eta") then
      !$acc parallel loop gang collapse(2)
      do k = 2,nz-1
        do i = 2,nx-1
          atdg_field(i,1,k) = 0.0d0; btdg_field(i,1,k) = 1.0d0; ctdg_field(i,1,k) = C3B(1)
          atdg_field(i,ny,k) = C3B(1); btdg_field(i,ny,k) = 1.0d0; ctdg_field(i,ny,k) = 0.0d0
    
          atdg_field(i,2,k) = AC4B(1); btdg_field(i,2,k) = 1.0d0; ctdg_field(i,2,k) = AC4B(1)
          atdg_field(i,ny-1,k) = AC4B(1); btdg_field(i,ny-1,k) = 1.0d0; ctdg_field(i,ny-1,k) = AC4B(1)
        enddo
      enddo
      !$acc parallel loop gang collapse(3)
      do ivar = 1,nvartdg
        do k = 2,nz-1
          do i = 2,nx-1
            rhstdg_field(i,1,k,ivar) = C3B(2)*temp(i,1,k,ivar) + C3B(3)*temp(i,2,k,ivar) + &
                        & C3B(4)*temp(i,3,k,ivar)
            rhstdg_field(i,ny,k,ivar) = -C3B(2)*temp(i,ny,k,ivar) - C3B(3)*temp(i,ny-1,k,ivar) - &
                          & C3B(4)*temp(i,ny-2,k,ivar)
            !rhstdg_field(i,2,k,ivar) = 0.50d0*AC4B(2)*(temp(i,3,k,ivar) - temp(i,1,k,ivar))
            !rhstdg_field(i,ny-1,k,ivar) = 0.50d0*AC4B(2)*(temp(i,ny,k,ivar) - temp(i,ny-2,k,ivar))
            rhstdg_field(i,2,k,ivar) = AC4B(2)*temp(i,1,k,ivar) + AC4B(4)*temp(i,3,k,ivar)
            rhstdg_field(i,ny-1,k,ivar) = -AC4B(2)*temp(i,ny,k,ivar) - AC4B(4)*temp(i,ny-2,k,ivar)
          enddo
        enddo
      enddo

    elseif (dir .eq. "zeta") then
      !$acc parallel loop gang collapse(2)
      do j = 2,ny-1
        do i = 2,nx-1
          atdg_field(i,j,1) = 0.0d0; btdg_field(i,j,1) = 1.0d0; ctdg_field(i,j,1) = C3B(1)
          atdg_field(i,j,nz) = C3B(1); btdg_field(i,j,nz) = 1.0d0; ctdg_field(i,j,nz) = 0.0d0
    
          atdg_field(i,j,2) = AC4B(1); btdg_field(i,j,2) = 1.0d0; ctdg_field(i,j,2) = AC4B(1)
          atdg_field(i,j,nz-1) = AC4B(1); btdg_field(i,j,nz-1) = 1.0d0; ctdg_field(i,j,nz-1) = AC4B(1)
        enddo
      enddo
      !$acc parallel loop gang collapse(3)
      do ivar = 1,nvartdg
        do j = 2,ny-1
          do i = 2,nx-1
            rhstdg_field(i,j,1,ivar) = C3B(2)*temp(i,j,1,ivar) + C3B(3)*temp(i,j,2,ivar) + &
                        & C3B(4)*temp(i,j,3,ivar)
            rhstdg_field(i,j,nz,ivar) = -C3B(2)*temp(i,j,nz,ivar) - C3B(3)*temp(i,j,nz-1,ivar) - &
                          & C3B(4)*temp(i,j,nz-2,ivar)
            !rhstdg_field(i,j,2,ivar) = 0.50d0*AC4B(2)*(temp(i,j,3,ivar) - temp(i,j,1,ivar))
            !rhstdg_field(i,j,nz-1,ivar) = 0.50d0*AC4B(2)*(temp(i,j,nz,ivar) - temp(i,j,nz-2,ivar))
            rhstdg_field(i,j,2,ivar) = AC4B(2)*temp(i,j,1,ivar) + AC4B(4)*temp(i,j,3,ivar)
            rhstdg_field(i,j,nz-1,ivar) = -AC4B(2)*temp(i,j,nz,ivar) - AC4B(4)*temp(i,j,nz-2,ivar)
          enddo
        enddo
      enddo

    endif

  end subroutine DerivBound



  subroutine ThomasAlg(nvartdg,dir)
    use mod_params
    use nvtx

    integer, intent(in)::nvartdg
    character(len=*), intent(in) :: dir
    integer :: i,j,k, ivar,l
    real (kind =8) :: betinv
    call nvtxStartRange("ThomasAlg")
    
    if (dir .eq. "xi") then
      !$acc parallel loop gang vector_length(32)
      do k = 1,nz 
      !$acc loop worker vector 
      do j = 1,ny
        betinv = 1.0d0/btdg_field(1,j,k)                   !! initial pivot
        !$acc loop seq
        do ivar = 1,nvartdg
          temp(1,j,k,ivar) = rhstdg_field(1,j,k,ivar)*betinv
        enddo

        !$acc loop seq
        do i = 2,nx                                      !!Forward Substitution
          gamtdg_field(i,j,k) = ctdg_field(i-1,j,k)*betinv
          betinv = 1.0d0/(btdg_field(i,j,k)-atdg_field(i,j,k)*gamtdg_field(i,j,k))
          !$acc loop seq
          do ivar = 1, nvartdg
            temp(i,j,k,ivar) = (rhstdg_field(i,j,k,ivar) - atdg_field(i,j,k)*temp(i-1,j,k,ivar))*betinv
          enddo
        enddo

        !$acc loop seq
        do i=nx-1,1,-1                                   !!Inverse Substitution
          !$acc loop seq
          do ivar=1,nvartdg
            temp(i,j,k,ivar) = temp(i,j,k,ivar) - gamtdg_field(i+1,j,k)*temp(i+1,j,k,ivar)
          enddo
        enddo
      enddo
      enddo

    elseif (dir .eq. "xi_old") then
      !$acc parallel loop gang collapse(2) vector_length(256)
      do k = 1,nz
      do j = 1,ny
        betinv = 1.0d0/btdg_field(1,j,k)                   !! initial pivot
        !$acc loop seq
        do ivar = 1,nvartdg
          temp(1,j,k,ivar) = rhstdg_field(1,j,k,ivar)*betinv
        enddo

        !$acc loop seq
        do i = 2,nx                                      !!Forward Substitution
          gamtdg_field(i,j,k) = ctdg_field(i-1,j,k)*betinv
          betinv = 1.0d0/(btdg_field(i,j,k)-atdg_field(i,j,k)*gamtdg_field(i,j,k))
          !$acc loop seq
          do ivar = 1, nvartdg
            temp(i,j,k,ivar) = (rhstdg_field(i,j,k,ivar) - atdg_field(i,j,k)*temp(i-1,j,k,ivar))*betinv
          enddo
        enddo

        !$acc loop seq
        do i=nx-1,1,-1                                   !!Inverse Substitution
          !$acc loop seq
          do ivar=1,nvartdg
            temp(i,j,k,ivar) = temp(i,j,k,ivar) - gamtdg_field(i+1,j,k)*temp(i+1,j,k,ivar)
          enddo
        enddo
      enddo
      enddo

    elseif (dir .eq. "eta") then
      !$acc parallel loop gang collapse(2)
      do k = 1,nz
      do i = 1,nx
        betinv = 1.0d0/btdg_field(i,1,k)                       !!initial pivot
        !$acc loop seq
        do l = 1,nvartdg
          temp(i,1,k,l) = rhstdg_field(i,1,k,l)*betinv
        enddo
    
        !$acc loop seq
        do j = 2,ny                                          !!Forward Substitution
          gamtdg_field(i,j,k) = ctdg_field(i,j-1,k)*betinv
          betinv = 1.0d0/(btdg_field(i,j,k)-atdg_field(i,j,k)*gamtdg_field(i,j,k))
          !$acc loop seq
          do l = 1, nvartdg
            temp(i,j,k,l) = (rhstdg_field(i,j,k,l) - atdg_field(i,j,k)*temp(i,j-1,k,l))*betinv
          enddo
        enddo
    
        !$acc loop seq
        do j=ny-1,1,-1                                        !!Inverse Substitution
          !$acc loop seq
          do l=1,nvartdg
            temp(i,j,k,l) = temp(i,j,k,l) - gamtdg_field(i,j+1,k)*temp(i,j+1,k,l)
          enddo
        enddo
      enddo
      enddo

    elseif (dir .eq. "zeta") then
      !$acc parallel loop gang collapse(2)
      do j = 1,ny
      do i = 1,nx
        betinv = 1.0d0/btdg_field(i,j,1)                       !!initial pivot
        !$acc loop seq
        do l = 1,nvartdg
          temp(i,j,1,l) = rhstdg_field(i,j,1,l)*betinv
        enddo
    
        !$acc loop seq
        do k = 2,nz                                          !!Forward Substitution
          gamtdg_field(i,j,k) = ctdg_field(i,j,k-1)*betinv
          betinv = 1.0d0/(btdg_field(i,j,k)-atdg_field(i,j,k)*gamtdg_field(i,j,k))
          !$acc loop seq
          do l = 1, nvartdg
            temp(i,j,k,l) = (rhstdg_field(i,j,k,l) - atdg_field(i,j,k)*temp(i,j,k-1,l))*betinv
          enddo
        enddo
    
        !$acc loop seq
        do k=nz-1,1,-1                                        !!Inverse Substitution
          !$acc loop seq
          do l=1,nvartdg
            temp(i,j,k,l) = temp(i,j,k,l) - gamtdg_field(i,j,k+1)*temp(i,j,k+1,l)
          enddo
        enddo
      enddo
      enddo

    endif
    call nvtxEndRange
    return
  end subroutine ThomasAlg
