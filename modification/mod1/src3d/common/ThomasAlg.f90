subroutine ThomasAlg_(nvartdg,ntdg,istart,iend)
  use mod_params
  implicit none

  integer, intent(in)::nvartdg,ntdg
  integer, intent(in)::istart,iend
  integer :: i,j,k, ivar
  real (kind =8), allocatable, dimension(:,:) :: gamtdg, betinv  
  real (kind =8), dimension(nmax) :: betinv_n  

  allocate(gamtdg(nmax,nmax)) 
  allocate(betinv(nmax,nmax)) 
  !!Solve TDG system with Thomas Algorithm
  
  do k=istart,iend
    betinv_n(1) = 1.0/btdg(k,1)                         ! initial pivot
    do ivar = 1,nvartdg
      tdg(k,1,ivar) = rhstdg(k,1,ivar)*betinv_n(1)
    enddo

    do i = 2,ntdg                                       ! Forward Elimination
      gamtdg(k,i) = ctdg(k,i-1)*betinv_n(i-1) 
      betinv_n(i) = 1.0/(btdg(k,i)-atdg(k,i)*gamtdg(k,i)) 
      do ivar = 1, nvartdg
        tdg(k,i,ivar) = (rhstdg(k,i,ivar) - atdg(k,i)*tdg(k,i-1,ivar))*betinv_n(i)
      enddo
    enddo
    
    do ivar=1,nvartdg                                   ! Backward Substitution
      do i=ntdg-1,1,-1
        tdg(k,i,ivar) = tdg(k,i,ivar) - gamtdg(k,i+1)*tdg(k,i+1,ivar)
      enddo
    enddo
  enddo

  deallocate(gamtdg) 
  deallocate(betinv) 
  return
end subroutine ThomasAlg_
