subroutine CompactFlux
  use mod_params
  use nvtx
  implicit none

  !allocate(temp(nmax,nmax,nmax,nvars))
  call nvtxStartRange("xicomp")
  call xicomp
  call nvtxEndRange
  
  call nvtxStartRange("etcomp")
  call etcomp
  call nvtxEndRange
  
  call nvtxStartRange("ztcomp")
  call ztcomp
  call nvtxEndRange
  !deallocate(temp, stat=ierr1)

contains

  subroutine xicomp
    integer :: i,j,k,ivar 
    real (kind =8) :: unor, rhou, rhoe

    !$acc parallel loop gang collapse(3) present(prim,xix,xiy,xiz)
    do k = 2,nz-1
      do j = 2,ny-1
        do i = 1,nx
          unor = xix(i,j,k)*prim(i,j,k,2) + xiy(i,j,k)*prim(i,j,k,3) + xiz(i,j,k)*prim(i,j,k,4)
          rhou = prim(i,j,k,1)*unor
          rhoe = (1+GAM1I)*prim(i,j,k,5) + 0.5*prim(i,j,k,1)*(prim(i,j,k,2)**2 + prim(i,j,k,3)**2 + prim(i,j,k,4)**2)
          temp(i,j,k,1) = rhou
          temp(i,j,k,2) = rhou*prim(i,j,k,2) + xix(i,j,k)*prim(i,j,k,5)
          temp(i,j,k,3) = rhou*prim(i,j,k,3) + xiy(i,j,k)*prim(i,j,k,5)
          temp(i,j,k,4) = rhou*prim(i,j,k,4) + xiz(i,j,k)*prim(i,j,k,5)
          temp(i,j,k,5) = rhoe*unor
        enddo
      enddo
    enddo

    call CompactDeriv(nvars,'xi')
    call ThomasAlg(nvars,'xi') 

    !$acc parallel loop gang collapse(4) present(rhs)
    do ivar = 1,nvars
      do k = 2,nz-1
        do j = 2,ny-1
          do i = 1,nx
            rhs(i,j,k,ivar) = rhs(i,j,k,ivar) + temp(i,j,k,ivar)
          enddo
        enddo
      enddo
    enddo
    
    return
  end subroutine xicomp
  
  subroutine etcomp
    integer :: i,j,k,ivar 
    real (kind =8) :: vnor, rhov, rhoe
    
    !$acc parallel loop gang collapse(3) present(prim,etax,etay,etaz)
    do k = 2,nz-1
      do j = 1,ny
        do i = 2,nx-1
          vnor = etax(i,j,k)*prim(i,j,k,2) + etay(i,j,k)*prim(i,j,k,3) + etaz(i,j,k)*prim(i,j,k,4)
          rhov = prim(i,j,k,1)*vnor
          rhoe = (1+GAM1I)*prim(i,j,k,5) + 0.5*prim(i,j,k,1)*(prim(i,j,k,2)**2 + prim(i,j,k,3)**2 + prim(i,j,k,4)**2)
          temp(i,j,k,1) = rhov
          temp(i,j,k,2) = rhov*prim(i,j,k,2) + etax(i,j,k)*prim(i,j,k,5)
          temp(i,j,k,3) = rhov*prim(i,j,k,3) + etay(i,j,k)*prim(i,j,k,5)
          temp(i,j,k,4) = rhov*prim(i,j,k,4) + etaz(i,j,k)*prim(i,j,k,5)
          temp(i,j,k,5) = rhoe*vnor
        enddo
      enddo
    enddo

    call CompactDeriv(nvars,'eta')
    call ThomasAlg(nvars,'eta') 

    !$acc parallel loop gang collapse(4) present(rhs)
    do ivar = 1,nvars
      do k = 2,nz-1
        do j = 1,ny
          do i = 2,nx-1
            rhs(i,j,k,ivar) = rhs(i,j,k,ivar) + temp(i,j,k,ivar)
          enddo
        enddo
      enddo
    enddo
    
    return
  end subroutine etcomp
  
  subroutine ztcomp
    integer :: i,j,k,ivar 
    real (kind =8) :: wnor, rhow, rhoe
    
    !$acc parallel loop gang collapse(3) present(prim,ztax,ztay,ztaz)
    do k = 1,nz
      do j = 2,ny-1
        do i = 2,nx-1
          wnor = ztax(i,j,k)*prim(i,j,k,2) + ztay(i,j,k)*prim(i,j,k,3) + ztaz(i,j,k)*prim(i,j,k,4)
          rhow = prim(i,j,k,1)*wnor
          rhoe = (1+GAM1I)*prim(i,j,k,5) + 0.5*prim(i,j,k,1)*(prim(i,j,k,2)**2 + prim(i,j,k,3)**2 + prim(i,j,k,4)**2)
          temp(i,j,k,1) = rhow
          temp(i,j,k,2) = rhow*prim(i,j,k,2) + ztax(i,j,k)*prim(i,j,k,5)
          temp(i,j,k,3) = rhow*prim(i,j,k,3) + ztay(i,j,k)*prim(i,j,k,5)
          temp(i,j,k,4) = rhow*prim(i,j,k,4) + ztaz(i,j,k)*prim(i,j,k,5)
          temp(i,j,k,5) = rhoe*wnor
        enddo
      enddo
    enddo

    call CompactDeriv(nvars,'zeta')
    call ThomasAlg(nvars,'zeta') 

    !$acc parallel loop gang collapse(4) present(rhs)
    do ivar = 1,nvars
      do k = 1,nz
        do j = 2,ny-1
          do i = 2,nx-1
            rhs(i,j,k,ivar) = rhs(i,j,k,ivar) + temp(i,j,k,ivar)
          enddo
        enddo
      enddo
    enddo
    
    return
  end subroutine ztcomp
  
end subroutine CompactFlux
