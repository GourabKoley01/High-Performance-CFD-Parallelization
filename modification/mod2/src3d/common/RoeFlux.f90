subroutine RoeFlux
  use mod_params
  use nvtx
  implicit none

  integer :: i,j,k,l
  real(kind=8):: fluxi, fluxj, fluxk

!  call Recon
  !print*,'6recon - roe',maxval(rhs)

  
  !$acc parallel loop collapse(4) !!async(2)
  do l=1,nvars
    do k=1,nz
      do j=1,ny
        do i=1,nx
          xiflux(i,j,k,l) = 0.0d0
          etflux(i,j,k,l) = 0.0d0
          ztflux(i,j,k,l) = 0.0d0
        enddo
      enddo
    enddo
  enddo
  !!$acc wait(2)
  
  !xiflux = 0.0d0
  !etflux = 0.0d0
  !ztflux = 0.0d0

!  do k=1,nz-1
!    do j=1,ny-1
!      do i=1,nx-1
!        call Roe3D(1,xiflux(i,j,k,:))
!        call Roe3D(2,etflux(i,j,k,:))
!        call Roe3D(3,ztflux(i,j,k,:))
        !        call Roe3D(primLi(i,j,k,:),primRi(i,j,k,:),xix(i,j,k),xiy(i,j,k),xiz(i,j,k),&
        !          &xiflux(i,j,k,:))
        !        call Roe3D(primLj(i,j,k,:),primRj(i,j,k,:),etax(i,j,k),etay(i,j,k),etaz(i,j,k),&
        !          &etflux(i,j,k,:))
        !        call Roe3D(primLk(i,j,k,:),primRk(i,j,k,:),ztax(i,j,k),ztay(i,j,k),ztaz(i,j,k),&
        !          &ztflux(i,j,k,:))
!      enddo
!    enddo
!  enddo
  
  call nvtxStartRange("xiroe")
  call xiroe
  call nvtxEndRange
 
  call nvtxStartRange("etroe")
  call etroe
  call nvtxEndRange

  call nvtxStartRange("ztroe")
  call ztroe
  call nvtxEndRange

  !$acc parallel loop collapse(4) &
  !$acc& present(xjac) &
  !$acc& private(fluxk,fluxj,fluxi)
  do l=1,nvars
    do k=2,nz-1
      do j=2,ny-1
        do i=2,nx-1
          fluxi = (xiflux(i,j,k,l)-xiflux(i-1,j,k,l))
          fluxj = (etflux(i,j,k,l)-etflux(i,j-1,k,l))
          fluxk = (ztflux(i,j,k,l)-ztflux(i,j,k-1,l))
          rhs(i,j,k,l) = rhs(i,j,k,l) + (fluxi+fluxj+fluxk)/xjac(i,j,k)
        enddo
      enddo
    enddo
  enddo

end subroutine RoeFlux
