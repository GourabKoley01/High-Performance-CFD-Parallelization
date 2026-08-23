subroutine metric3d
  use mod_params
  implicit none

  real(kind = 8):: xxi, xeta, xzta
  real(kind = 8):: yxi, yeta, yzta
  real(kind = 8):: zxi, zeta, zzta
  real(kind = 8):: xixj, xiyj, xizj
  real(kind = 8):: xjacmin, xjacinv
  integer:: i,j,k

  if (imetric_scheme .eq. 10) then
    ! xi CompactDerivatives
    do j =1,ny
      do k =1,nz
        do i=1,nx
          tdg(k,i,1) = x(i,j,k)
          tdg(k,i,2) = y(i,j,k)
          tdg(k,i,3) = z(i,j,k)
        enddo
      enddo

      call CompactDeriv(3,nx,iper,1,nz)

      do k =1,nz
        do i=1,nx
          xix(i,j,k) = tdg(k,i,1)
          xiy(i,j,k) = tdg(k,i,2)
          xiz(i,j,k) = tdg(k,i,3)
        enddo
      enddo
    enddo


    ! eta CompactDerivatives
    do k =1,nz
      do i =1,nx
        do j=1,ny
          tdg(i,j,1) = x(i,j,k)
          tdg(i,j,2) = y(i,j,k)
          tdg(i,j,3) = z(i,j,k)
        enddo
      enddo

      call CompactDeriv(3,ny,jper,1,nx)

      do i =1,nx
        do j=1,ny
          etax(i,j,k) = tdg(i,j,1)
          etay(i,j,k) = tdg(i,j,2)
          etaz(i,j,k) = tdg(i,j,3)
        enddo
      enddo
    enddo

    ! zeta CompactDerivatives
    do j=1,ny
      do k =1,nz
        do i =1,nx
          tdg(i,k,1) = x(i,j,k)
          tdg(i,k,2) = y(i,j,k)
          tdg(i,k,3) = z(i,j,k)
        enddo
      enddo

      call CompactDeriv(3,nz,kper,1,nx)

      do k =1,nz
        do i=1,nx
          ztax(i,j,k) = tdg(i,k,1)
          ztay(i,j,k) = tdg(i,k,2)
          ztaz(i,j,k) = tdg(i,k,3)
        enddo
      enddo
    enddo

  do k =1,nz
    do j =1,ny
      do i=1,nx
        xxi = xix(i,j,k)
        yxi = xiy(i,j,k)
        zxi = xiz(i,j,k)

        xeta = etax(i,j,k)
        yeta = etay(i,j,k)
        zeta = etaz(i,j,k)

        xzta = ztax(i,j,k)
        yzta = ztay(i,j,k)
        zzta = ztaz(i,j,k)

        xixj = (yeta*zzta - yzta*zeta)
        xiyj = (xzta*zeta - xeta*zzta)
        !xizj = (xeta*yzta - zzta*yeta)
        xizj = (xeta*yzta - xzta*yeta)

        xjac(i,j,k) = (xxi*xixj+yxi*xiyj+zxi*xizj)

        xjacinv = 1.0d0/xjac(i,j,k)

        xix(i,j,k) = xjacinv*(yeta*zzta - yzta*zeta)
        xiy(i,j,k) = xjacinv*(xzta*zeta - xeta*zzta)
        xiz(i,j,k) = xjacinv*(xeta*yzta - xzta*yeta)
        etax(i,j,k) = xjacinv*(yzta*zxi - yxi*zzta)
        etay(i,j,k) = xjacinv*(xxi*zzta - xzta*zxi)
        etaz(i,j,k) = xjacinv*(xzta*yxi - xxi*yzta)
        ztax(i,j,k) = xjacinv*(yxi*zeta - yeta*zxi)
        ztay(i,j,k) = xjacinv*(xeta*zxi - xxi*zeta)
        ztaz(i,j,k) = xjacinv*(xxi*yeta - xeta*yxi)

      enddo
    enddo
  enddo

  else if (imetric_scheme .eq. 2) then
    !call met2nd
    call MetricSecond
!    do k=1,nz
!      do j=1,ny
!        do i=1,nx
!          xjacinv = 1.0d0/xjac(i,j,k)
!          xix(i,j,k) = xjacinv*xix(i,j,k) 
!          xiy(i,j,k) = xjacinv*xiy(i,j,k) 
!          xiz(i,j,k) = xjacinv*xiz(i,j,k) 
!          etax(i,j,k) = xjacinv*etax(i,j,k) 
!          etay(i,j,k) = xjacinv*etay(i,j,k) 
!          etaz(i,j,k) = xjacinv*etaz(i,j,k) 
!          ztax(i,j,k) = xjacinv*ztax(i,j,k) 
!          ztay(i,j,k) = xjacinv*ztay(i,j,k) 
!          ztaz(i,j,k) = xjacinv*ztaz(i,j,k) 
!        enddo
 !     enddo
  !  enddo
  !  do j = 2,ny-1 
  !    jp = j+1
  !    jm = j-1
  !    do i = 2,nx-1 
  !      ip = i+1
  !      im = i-1
  !      xxi(i,j) = x(ip,j) - x(im,j)
  !      yxi(i,j) = y(ip,j) - y(im,j)
  !      xeta(i,j) = x(i,jp) - x(i,jm)
  !      yeta(i,j) = y(i,jp) - y(i,jm)
  !    enddo
  !  enddo

    !!Edges one side differentation
  !  do i = 2,nx-1 
  !    ip = i+1
  !    im = i-1
  !    xxi(i,1) = x(ip,1) - x(im,1)
  !    yxi(i,1) = y(ip,1) - y(im,1)
  !    xeta(i,1) = -x(i,3) + 4.*x(i,2) -3.*x(i,1)
  !    yeta(i,1) = -y(i,3) + 4.*y(i,2) -3.*y(i,1)

  !    xxi(i,ny) = x(ip,ny) - x(im,ny)
  !    yxi(i,ny) = y(ip,ny) - y(im,ny)
  !    xeta(i,ny) = -1*(-x(i,ny-2) + 4.*x(i,ny-1) -3.*x(i,ny))
  !    yeta(i,ny) = -1*(-y(i,ny-2) + 4.*y(i,ny-1) -3.*y(i,ny))
  !  enddo

  !  do j = 2,ny-1 
  !    jp = j+1
  !    jm = j-1
  !    xxi(1,j) = -x(3,j) + 4.*x(2,j) - 3.*x(1,j)
  !    yxi(1,j) = -y(3,j) + 4.*y(2,j) - 3.*y(1,j)
  !    xeta(1,j) = x(1,jp) - x(1,jm)
  !    yeta(1,j) = y(1,jp) - y(1,jm)

  !    xxi(nx,j) = -1*(-x(nx-2,j) + 4.*x(nx-1,j) - 3.*x(nx,j))
  !    yxi(nx,j) = -1*(-y(nx-2,j) + 4.*y(nx-1,j) - 3.*y(nx,j))
  !    xeta(nx,j) = x(nx,jp) - x(nx,jm)
  !    yeta(nx,j) = y(nx,jp) - y(nx,jm)
  !  enddo

  !  xxi(1,1) = -x(3,1) + 4.*x(2,1) - 3.*x(1,1)
  !  yxi(1,1) = -y(3,1) + 4.*y(2,1) - 3.*y(1,1)
  !  xeta(1,1) = -x(1,3) + 4.*x(1,2) - 3.*x(1,1)
  !  yeta(1,1) = -y(1,3) + 4.*y(1,2) - 3.*y(1,1)

  !  xxi(1,ny) = -x(3,ny) + 4.*x(2,ny) - 3.*x(1,ny)
  !  yxi(1,ny) = -y(3,ny) + 4.*y(2,ny) - 3.*y(1,ny)
  !  xeta(1,ny) = -1*(-x(1,ny-2) + 4.*x(1,ny-1) - 3.*x(1,ny))
  !  yeta(1,ny) = -1*(-y(1,ny-2) + 4.*y(1,ny-1) - 3.*y(1,ny))

  !  xxi(nx,1) = -1*(-x(nx-2,1) + 4.*x(nx-1,1) - 3.*x(nx,1))
  !  yxi(nx,1) = -1*(-y(nx-2,1) + 4.*y(nx-1,1) - 3.*y(nx,1))
  !  xeta(nx,1) = -x(nx,3) + 4.*x(nx,2) - 3.*x(nx,1)
  !  yeta(nx,1) = -y(nx,3) + 4.*y(nx,2) - 3.*y(nx,1)

  !  xxi(nx,ny) = -1*(-x(nx-2,ny) + 4.*x(nx-1,ny) - 3.*x(nx,ny))
  !  yxi(nx,ny) = -1*(-y(nx-2,ny) + 4.*y(nx-1,ny) - 3.*y(nx,ny))
  !  xeta(nx,ny) = -1*(-x(nx,ny-2) + 4.*x(nx,ny-1) - 3.*x(nx,ny))
  !  yeta(nx,ny) = -1*(-y(nx,ny-2) + 4.*y(nx,ny-1) - 3.*y(nx,ny))

  !  xxi = xxi*0.5
  !  yxi = yxi*0.5
  !  xeta = xeta*0.5
  !  yeta = yeta*0.5
  !print*,'This scheme is not implemented yet'
  endif


  if(idebug .eq. 1)then
    print*,'max(xix), min(xix)',maxval(xix), minval(xix)
    print*,'max(xiy), min(xiy)',maxval(xiy), minval(xiy)
    print*,'max(xiz), min(xiz)',maxval(xiz), minval(xiz)
    print*,'max(etax),min(etax)',maxval(etax), minval(etax)
    print*,'max(etay),min(etay)',maxval(etay), minval(etay)
    print*,'max(etaz),min(etaz)',maxval(etaz), minval(etaz)
    print*,'max(ztax),min(ztax)',maxval(ztax), minval(ztax)
    print*,'max(ztay),min(ztay)',maxval(ztay), minval(ztay)
    print*,'max(ztaz),min(ztaz)',maxval(ztaz), minval(ztaz)
    print*,'Jmax, Jmin', maxval(xjac), minval(xjac)
  endif


  xjacmin = minval(xjac)

  if(xjacmin .le. 0.0) then
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          if (xjac(i,j,k) .le. 0.0) then
            print*,'Negative Jacobian Encountered in i=',i, ' j=',j, ' k=',k
            print*,'Jacobian =' , xjac(i,j,k)
            print*,'PROGRAM STOPPED'
          endif
        enddo
      enddo
    enddo
    stop
  endif



  !do i = 1,nx
  !  do j = 1,ny
  !    print*, xjac(i,j)
  !  enddo
  !enddo
  !do k = 50,55
  !  do j = 50,55
  !    do i = 50,55
  !      print*, xix(i,j,k),etax(i,j,k),ztax(i,j,k),xjac(i,j,k)
  !    enddo
  !  enddo
  !enddo
end subroutine metric3d

subroutine MetricSecond
  use mod_params
  implicit none
  integer:: i,j,k
  real (kind=8),allocatable,dimension(:,:,:) :: xxi, yxi, zxi
  real (kind=8),allocatable,dimension(:,:,:) :: xeta, yeta, zeta
  real (kind=8),allocatable,dimension(:,:,:) :: xzta, yzta, zzta
  !real (kind=8) :: xixj, xiyj, xizj
  real (kind=8) :: xjacinv

   integer :: ip,im,jp,jm,kp,km,iFlg,jFlg,kFlg
   !Subroutine calculates the grid metrics for the entire grid.
   !Uses 2 pt central difference to calculate partials for interior points.
   !Use 2 pt forward or backward differences to calculate partials at boundary,
   !unless that boundary is defined as axisymmetric, in which case the central
   !difference is again used. 

   allocate(xxi(nx,ny,nz))
   allocate(yxi(nx,ny,nz))
   allocate(zxi(nx,ny,nz))
   allocate(xeta(nx,ny,nz))
   allocate(yeta(nx,ny,nz))
   allocate(zeta(nx,ny,nz))
   allocate(xzta(nx,ny,nz))
   allocate(yzta(nx,ny,nz))
   allocate(zzta(nx,ny,nz))

   do i=1,nx
      ip=i+1; im=i-1; iFlg=2.0d0
      if (i .eq. nx) then
         if (iper .eq. 1) then
            ip=2; iFlg=2.0d0
         else
            ip=i; iFlg=1.0d0
         endif
      endif
      if (i .eq. 1) then
         if (iper .eq. 1) then
            im=nx-1; iFlg=2.0d0
         else
            im=i; iFlg=1.0d0
         endif
      endif
      do j=1,ny
         jp=j+1; jm=j-1; jFlg=2.0d0
         if (j .eq. ny) then
            if (jper .eq. 1) then
               jp=2; jFlg=2.0d0
            else
               jp=j; jFlg=1.0d0
            endif
         endif
         if (j .eq. 1) then
            if (jper .eq. 1) then
               jm=ny-1; jFlg=2.0d0
            else
               jm=j; jFlg=1.0d0
            endif
         endif
         do k=1,nz
            kp=k+1; km=k-1; kFlg=2.0d0
            if (k .eq. nz) then
               if (kper .eq. 1) then
                  kp=2; kFlg=2.0d0
               else
                  kp=k; kFlg=1.0d0 
               endif
            endif
            if (k .eq. 1) then
               if (kper .eq. 1) then
                  km=nz-1; kFlg=2.0d0
               else
                  km=k; kFlg=1.0d0
               endif
            endif

            !Calculate partial derivatives
            xxi(i,j,k)=(x(ip,j,k)-x(im,j,k))/iFlg
            xeta(i,j,k)=(x(i,jp,k)-x(i,jm,k))/jFlg
            xzta(i,j,k)=(x(i,j,kp)-x(i,j,km))/kFlg
            yxi(i,j,k)=(y(ip,j,k)-y(im,j,k))/iFlg     
            yeta(i,j,k)=(y(i,jp,k)-y(i,jm,k))/jFlg                  
            yzta(i,j,k)=(y(i,j,kp)-y(i,j,km))/kFlg                  
            zxi(i,j,k)=(z(ip,j,k)-z(im,j,k))/iFlg
            zeta(i,j,k)=(z(i,jp,k)-z(i,jm,k))/jFlg
            zzta(i,j,k)=(z(i,j,kp)-z(i,j,km))/kFlg
         enddo
      enddo
   enddo
   !Calculate determinate of the jacobian matrix and use Cramer's rule to find the
   !inverse matrix with the grid metrics
   do i=1,nx
      do j=1,ny
         do k=1,nz
            xjac(i,j,k)=&
               &xxi(i,j,k)*(yeta(i,j,k)*zzta(i,j,k)-yzta(i,j,k)*zeta(i,j,k))+&
               &xeta(i,j,k)*(yzta(i,j,k)*zxi(i,j,k)-yxi(i,j,k)*zzta(i,j,k))+&
               &xzta(i,j,k)*(yxi(i,j,k)*zeta(i,j,k)-yeta(i,j,k)*zxi(i,j,k))
            xjacinv =1.0d0/xjac(i,j,k)
            xix(i,j,k)=(yeta(i,j,k)*zzta(i,j,k)-yzta(i,j,k)*zeta(i,j,k))*xjacinv
            xiy(i,j,k)=-(xeta(i,j,k)*zzta(i,j,k)-xzta(i,j,k)*zeta(i,j,k))*xjacinv
            xiz(i,j,k)=(xeta(i,j,k)*yzta(i,j,k)-xzta(i,j,k)*yeta(i,j,k))*xjacinv
            etax(i,j,k)=-(yxi(i,j,k)*zzta(i,j,k)-yzta(i,j,k)*zxi(i,j,k))*xjacinv           
            etay(i,j,k)=(xxi(i,j,k)*zzta(i,j,k)-xzta(i,j,k)*zxi(i,j,k))*xjacinv            
            etaz(i,j,k)=-(xxi(i,j,k)*yzta(i,j,k)-xzta(i,j,k)*yxi(i,j,k))*xjacinv           
            ztax(i,j,k)=(yxi(i,j,k)*zeta(i,j,k)-yeta(i,j,k)*zxi(i,j,k))*xjacinv
            ztay(i,j,k)=-(xxi(i,j,k)*zeta(i,j,k)-xeta(i,j,k)*zxi(i,j,k))*xjacinv
            ztaz(i,j,k)=(xxi(i,j,k)*yeta(i,j,k)-xeta(i,j,k)*yxi(i,j,k))*xjacinv 
         enddo
      enddo
   enddo
end subroutine MetricSecond
