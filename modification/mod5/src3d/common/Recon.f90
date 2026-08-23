subroutine Recon(idir)
  use mod_params
  implicit none

  integer, intent(in) :: idir
  integer :: i,j,k, ivar
  real(kind =8):: reconm, reconp
  real(kind =8):: dprimm, dprimp, rm, rp, vl

  reconp = 1.0d0 + recon_ratio
  reconm = 1.0d0 - recon_ratio

  !!$acc kernels present(prim)
  !primL = prim
  !primR = prim
  !!$acc end kernels

  if(idir .eq. 1)then
    if(iroe_order .eq. 1)then
      !$acc parallel loop collapse(4) &
      !$acc& present(prim)
      do ivar=1,nvars
        do k =1, nz-1
          do j =1, ny-1
            do i =1, nx-1
              primL(i,j,k,ivar) = prim(i,j,k,ivar) 
              primR(i,j,k,ivar) = prim(i+1,j,k,ivar) 
            enddo
          enddo
        enddo
      enddo

    else if(iroe_order .eq. 3)then
      !$acc parallel loop collapse(4) &
      !$acc& present(prim, ifblock) &
      !$acc& private(dprimp,dprimm,rp,rm)
      do ivar = 1, nvars
        do k =1, nz
          do j =1, ny
            do i =2, nx-1
              dprimm = prim(i,j,k,ivar) - prim(i-1,j,k,ivar)
              dprimp = prim(i+1,j,k,ivar) - prim(i,j,k,ivar)
              rm = dprimm / (dprimp+eps)
              rp = dprimp / (dprimm+eps)
              dprimm = dprimm*vl(rp)
              dprimp = dprimp*vl(rm)
              if (ifblock(i,j,k) .eq. 1)then
                primL(i,j,k,ivar) = prim(i,j,k,ivar) 
                primR(i-1,j,k,ivar) = prim(i,j,k,ivar)
              else
                primL(i,j,k,ivar) = prim(i,j,k,ivar) +0.250d0*(reconm*dprimm +&
                  & reconp*dprimp) 
                primR(i-1,j,k,ivar) = prim(i,j,k,ivar) -0.250d0*(reconp*dprimm +&
                  & reconm*dprimp) 
              endif
            enddo
          enddo
        enddo
      enddo

      !$acc parallel loop collapse(3) &
      !$acc& present(prim)
      do ivar = 1, nvars
        do k =1, nz
          do j =1, ny
            primL(1,j,k,ivar) = prim(1,j,k,ivar) 
            !        primR(1,j,ivar) = prim(2,j,ivar) 
            !        primL(nx-1,j,ivar) = prim(nx-1,j,ivar) 
            primR(nx-1,j,k,ivar) = prim(nx,j,k,ivar) 
          enddo
        enddo
      enddo
    endif
  endif

  if(idir .eq. 2)then
    if(iroe_order .eq. 1)then
      !$acc parallel loop collapse(4) &
      !$acc& present(prim)
      do ivar=1,nvars
        do k =1, nz-1
          do j =1, ny-1
            do i =1, nx-1
              primL(i,j,k,ivar) = prim(i,j,k,ivar) 
              primR(i,j,k,ivar) = prim(i,j+1,k,ivar) 
            enddo
          enddo
        enddo
      enddo

    else if(iroe_order .eq. 3)then
      !$acc parallel loop collapse(4) &
      !$acc& present(prim, ifblock) &
      !$acc& private(dprimp,dprimm,rp,rm)
      do ivar=1,nvars
        do k =1, nz
          do j =2, ny-1
            do i =1, nx
              dprimm = prim(i,j,k,ivar) - prim(i,j-1,k,ivar)
              dprimp = prim(i,j+1,k,ivar) - prim(i,j,k,ivar)
              rm = dprimm / (dprimp+eps)
              rp = dprimp / (dprimm+eps)
              dprimm = dprimm*vl(rp)
              dprimp = dprimp*vl(rm)
              if (ifblock(i,j,k) .eq. 1)then
                primL(i,j,k,ivar) = prim(i,j,k,ivar) 
                primR(i,j-1,k,ivar) = prim(i,j,k,ivar) 
              else
                primL(i,j,k,ivar) = prim(i,j,k,ivar) +0.250d0*(reconm*dprimm +&
                  & reconp*dprimp) 
                primR(i,j-1,k,ivar) = prim(i,j,k,ivar) -0.250d0*(reconp*dprimm +&
                  & reconm*dprimp) 
              endif
            enddo
          enddo
        enddo
      enddo
      
      !$acc parallel loop collapse(3) &
      !$acc& present(prim)
      do ivar=1,nvars 
        do k =1, nz
          do i =1, nx
            primL(i,1,k,ivar) = prim(i,1,k,ivar) 
            !        primR(i,1,ivar) = prim(i,2,ivar) 
            !        primL(i,ny-1,ivar) = prim(i,ny-1,ivar) 
            primR(i,ny-1,k,ivar) = prim(i,ny,k,ivar) 
          enddo
        enddo
      enddo
      !print*, reconm, reconp, rm, rp, vl(rm), vl(rp)
    endif
  endif

  if(idir .eq. 3)then
    if(iroe_order .eq. 1)then
      !$acc parallel loop collapse(4) &
      !$acc& present(prim)
      do ivar=1,nvars
        do k =1, nz-1
          do j =1, ny-1
            do i =1, nx-1
              primL(i,j,k,ivar) = prim(i,j,k,ivar) 
              primR(i,j,k,ivar) = prim(i,j,k+1,ivar) 
            enddo
          enddo
        enddo
      enddo

    else if(iroe_order .eq. 3)then
      !$acc parallel loop collapse(4) &
      !$acc& present(prim, ifblock) &
      !$acc& private(dprimp,dprimm,rp,rm)
      do ivar=1,nvars
        do k =2, nz-1
          do j =1, ny
            do i =1, nx
              dprimm = prim(i,j,k,ivar) - prim(i,j,k-1,ivar)
              dprimp = prim(i,j,k+1,ivar) - prim(i,j,k,ivar)
              rm = dprimm / (dprimp+eps)
              rp = dprimp / (dprimm+eps)
              dprimm = dprimm*vl(rp)
              dprimp = dprimp*vl(rm)
              if (ifblock(i,j,k) .eq. 1)then
                primL(i,j,k,ivar) = prim(i,j,k,ivar) 
                primR(i,j,k-1,ivar) = prim(i,j,k,ivar) 
              else
                primL(i,j,k,ivar) = prim(i,j,k,ivar) +0.250d0*(reconm*dprimm +&
                  & reconp*dprimp) 
                primR(i,j,k-1,ivar) = prim(i,j,k,ivar) -0.250d0*(reconp*dprimm +&
                  & reconm*dprimp) 
              endif
            enddo
          enddo
        enddo
      enddo

      !$acc parallel loop collapse(3) &
      !$acc& present(prim)
      do ivar=1,nvars
        do j =1, ny
          do i =1, nx
            primL(i,j,1,ivar) = prim(i,j,1,ivar) 
            !        primR(i,1,ivar) = prim(i,2,ivar) 
            !        primL(i,ny-1,ivar) = prim(i,ny-1,ivar) 
            primR(i,j,nz-1,ivar) = prim(i,j,nz,ivar) 
          enddo
        enddo
      enddo
      !print*, reconm, reconp, rm, rp, vl(rm), vl(rp)
    endif
  endif
end subroutine Recon


pure function vl(a)
  implicit none
  !$acc routine
  real(kind=8) :: vl
  real(kind=8),intent(in) :: a
  vl = (a+dabs(a))/(1.0d0+dabs(a))
end function vl

pure function va(a)
  implicit none
  !$acc routine
  real(kind=8) :: va
  real(kind=8),intent(in) :: a
  if(a .gt. 0)then 
    va = (a+a*a)/(1.0d0+a*a)
  else
    va = 0.0d0
  endif
end function va
