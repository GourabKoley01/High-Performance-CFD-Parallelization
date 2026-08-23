subroutine ApplyBC
  use mod_params
  use mod_spmd
  implicit none
  integer:: i,j,k
  real(kind=8) :: aval, velnor, xnor, ynor, mag
  !real(kind=8) :: uprof(133)
  !character(len=80)::filename

  !!!!!YET TO PUT MULTIBLOCK 3D
  !if(nblnks .gt. 0) call MultiblockBC


!FOR LDC all walls as of now

!$acc parallel present(prim, q)
!LEFT
  if(pleft .lt. 0)then
    !$acc loop collapse(3)      
    do k=1,nz
      do j = 1,ny
        do i = 1,1
          prim(i,j,k,1) = prim(i+1,j,k,1)
          prim(i,j,k,2) = 0.0
          prim(i,j,k,3) = 0.0
          prim(i,j,k,4) = 0.0
          prim(i,j,k,5) = prim(i+1,j,k,5)
        enddo
      enddo
    enddo
    !$acc loop collapse(3)      
    do k=1,nz
      do j = 1,ny
        do i = 1,1
          q(i,j,k,1) = prim(i,j,k,1)
          q(i,j,k,2) = prim(i,j,k,2)*prim(i,j,k,1) 
          q(i,j,k,3) = prim(i,j,k,3)*prim(i,j,k,1)
          q(i,j,k,4) = prim(i,j,k,4)*prim(i,j,k,1)
          q(i,j,k,5) = GAM1I * prim(i,j,k,5)  + 0.50d0*prim(i,j,k,1)*&
            & (prim(i,j,k,2)**2+prim(i,j,k,3)**2+prim(i,j,k,4)**2)
        enddo
      enddo
    enddo
  endif

!RIGHT
  if(pright .lt. 0)then
    !$acc loop collapse(3)      
    do k=1,nz
      do j = 1,ny
        do i = nx,nx
          prim(nx,j,k,1) = prim(nx-1,j,k,1)
          prim(nx,j,k,2) = 0.0
          prim(nx,j,k,3) = 0.0
          prim(nx,j,k,4) = 0.0
          prim(nx,j,k,5) = prim(nx-1,j,k,5)
        enddo
      enddo
    enddo
    !$acc loop collapse(3)      
    do k=1,nz
      do j = 1,ny
        do i = nx,nx
          q(i,j,k,1) = prim(i,j,k,1)
          q(i,j,k,2) = prim(i,j,k,2)*prim(i,j,k,1) 
          q(i,j,k,3) = prim(i,j,k,3)*prim(i,j,k,1)
          q(i,j,k,4) = prim(i,j,k,4)*prim(i,j,k,1)
          q(i,j,k,5) = GAM1I * prim(i,j,k,5)  + 0.50d0*prim(i,j,k,1)*&
            & (prim(i,j,k,2)**2+prim(i,j,k,3)**2+prim(i,j,k,4)**2)
        enddo
      enddo
    enddo
  endif

!NORTH
  if(pabove .lt. 0)then
    !prim(:,ny,:,2) = 0.0d0
    !do k=1,nz
    !  do i=1,nx
        !prim(i,ny,k,2) = 16*x(i,ny,k)**2*(1-x(i,ny,k))**2
      !enddo
    !enddo
    !prim(:,ny,:,3) = 0.0d0
    
    !$acc loop collapse(3)      
    do k=1,nz
      do j = ny,ny
        do i = 1,nx
          prim(i,ny,k,1) = prim(i,ny-1,k,1)
          prim(i,ny,k,2) = 0.0d0
          prim(i,ny,k,3) = 0.0d0
          prim(i,ny,k,4) = 0.0d0
          prim(i,ny,k,5) = prim(i,ny-1,k,5)
        enddo
      enddo
    enddo
    !$acc loop collapse(3)      
    do k=1,nz
      do j = ny,ny
        do i = 1,nx
          q(i,j,k,1) = prim(i,j,k,1)
          q(i,j,k,2) = prim(i,j,k,2)*prim(i,j,k,1) 
          q(i,j,k,3) = prim(i,j,k,3)*prim(i,j,k,1)
          q(i,j,k,4) = prim(i,j,k,4)*prim(i,j,k,1)
          q(i,j,k,5) = GAM1I * prim(i,j,k,5)  + 0.50d0*prim(i,j,k,1)*&
            & (prim(i,j,k,2)**2+prim(i,j,k,3)**2+prim(i,j,k,4)**2)
        enddo
      enddo
    enddo
  endif

!SOUTH
  if(pbelow .lt. 0)then
    
    !$acc loop collapse(3)      
    do k=1,nz
      do j = 1,1
        do i = 1,nx
          prim(i,1,k,1) = prim(i,2,k,1)
          prim(i,1,k,2) = 0.0d0
          prim(i,1,k,3) = 0.0d0
          prim(i,1,k,4) = 0.0d0
          prim(i,1,k,5) = prim(i,2,k,5)
        enddo
      enddo
    enddo
    !$acc loop collapse(3)      
    do k=1,nz
      do j = 1,1
        do i = 1,nx
          q(i,j,k,1) = prim(i,j,k,1)
          q(i,j,k,2) = prim(i,j,k,2)*prim(i,j,k,1) 
          q(i,j,k,3) = prim(i,j,k,3)*prim(i,j,k,1)
          q(i,j,k,4) = prim(i,j,k,4)*prim(i,j,k,1)
          q(i,j,k,5) = GAM1I * prim(i,j,k,5)  + 0.50d0*prim(i,j,k,1)*&
            & (prim(i,j,k,2)**2+prim(i,j,k,3)**2+prim(i,j,k,4)**2)
        enddo
      enddo
    enddo
  endif

!UP slip
  if(pup .lt. 0)then
    
    !$acc loop collapse(3)      
    do k=nz,nz
      do j = 1,ny
        do i = 1,nx
          prim(i,j,nz,1) = prim(i,j,nz-1,1)
          prim(i,j,nz,2) = 1.0d0
          prim(i,j,nz,3) = 0.0d0
          prim(i,j,nz,4) = 0.0d0
          prim(i,j,nz,5) = prim(i,j,nz-1,5)
        enddo
      enddo
    enddo
    !$acc loop collapse(3)      
    do k=nz,nz
      do j = 1,ny
        do i = 1,nx
          q(i,j,k,1) = prim(i,j,k,1)
          q(i,j,k,2) = prim(i,j,k,2)*prim(i,j,k,1) 
          q(i,j,k,3) = prim(i,j,k,3)*prim(i,j,k,1)
          q(i,j,k,4) = prim(i,j,k,4)*prim(i,j,k,1)
          q(i,j,k,5) = GAM1I * prim(i,j,k,5)  + 0.50d0*prim(i,j,k,1)*&
            & (prim(i,j,k,2)**2+prim(i,j,k,3)**2+prim(i,j,k,4)**2)
        enddo
      enddo
    enddo
  endif

  !DOWN
  if(pdown .lt. 0)then
    
    !$acc loop collapse(3)      
    do k=1,1
      do j = 1,ny
        do i = 1,nx
          prim(i,j,k,1) = prim(i,j,k+1,1)
          prim(i,j,k,2) = 0.0d0
          prim(i,j,k,3) = 0.0d0
          prim(i,j,k,4) = 0.0d0
          prim(i,j,k,5) = prim(i,j,k+1,5)
        enddo
      enddo
    enddo
    !$acc loop collapse(3)      
    do k=1,1
      do j = 1,ny
        do i = 1,nx
          q(i,j,k,1) = prim(i,j,k,1)
          q(i,j,k,2) = prim(i,j,k,2)*prim(i,j,k,1) 
          q(i,j,k,3) = prim(i,j,k,3)*prim(i,j,k,1)
          q(i,j,k,4) = prim(i,j,k,4)*prim(i,j,k,1)
          q(i,j,k,5) = GAM1I * prim(i,j,k,5)  + 0.50d0*prim(i,j,k,1)*&
            & (prim(i,j,k,2)**2+prim(i,j,k,3)**2+prim(i,j,k,4)**2)
        enddo
      enddo
    enddo
  endif

 !$acc end parallel
end subroutine ApplyBC
