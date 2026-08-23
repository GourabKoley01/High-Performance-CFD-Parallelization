!subroutine ConsToPrim
!  use mod_params
!  implicit none
!  real(kind=8),intent(in)::qt(5)
!  real(kind=8),intent(out)::primt(5)
!  real(kind = 8):: rhoinv
!  ! print*,'ConsToPrim',maxval(prim(:,:,:,1)), minval(prim(:,:,:,1))
!  primt(1) = qt(1)
!  rhoinv = 1.0d0/primt(1)
!  primt(2) = qt(2)*rhoinv 
!  primt(3) = qt(3)*rhoinv
!  primt(4) = qt(4)*rhoinv
!  primt(5) = GAM1 * (qt(5) -0.50d0*primt(1)*&
!    & (primt(2)**2+primt(3)**2+primt(4)**2))
!  return
!end subroutine ConsToPrim

subroutine ConsToPrim
  use mod_params
  implicit none
  real(kind = 8):: rhoinv
  integer:: i,j,k
  
  !$acc parallel loop gang collapse(3) present(prim,q) private(rhoinv) 
  do k = 2,nz-1 
    do j = 2,ny-1
      do i = 2,nx-1
        prim(i,j,k,1) = q(i,j,k,1)
        rhoinv = 1.0d0/prim(i,j,k,1)
        prim(i,j,k,2) = q(i,j,k,2)*rhoinv
        prim(i,j,k,3) = q(i,j,k,3)*rhoinv
        prim(i,j,k,4) = q(i,j,k,4)*rhoinv
        prim(i,j,k,5) = GAM1 * (q(i,j,k,5) -0.50d0*prim(i,j,k,1)*&
                & (prim(i,j,k,2)**2 + prim(i,j,k,3)**2 + prim(i,j,k,4)**2))
      enddo
    enddo
  enddo

  return
end subroutine ConsToPrim
