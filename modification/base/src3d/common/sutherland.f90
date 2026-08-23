subroutine sutherland
  use mod_params
  implicit none

  integer :: i,j,k
  !xmu = 0.0
  !real (kind=8):: T
  ! Call Sutherland's law here

  !$acc parallel loop collapse(3) present(prim)
  do k = 1,nz
    do j = 1,ny
      do i = 1,nx
        T(i,j,k) = GXM * prim(i,j,k,5)/prim(i,j,k,1)
        xmu(i,j,k) = T(i,j,k)**1.5d0*SMUP1/(T(i,j,k)+SMU)
        !xmu(i,j) = T(i,j)**0.70d0
      enddo
    enddo
  enddo
end subroutine sutherland
