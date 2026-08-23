subroutine CompactVflux
  use mod_params
  implicit none

  integer :: i,j,k,ivar
  real (kind=8):: uder, vder, wder, Tder
  real (kind=8):: derx2, dery2, derz2
  real (kind=8):: derxy, deryz, derzx, dertot
  real (kind=8):: xmuxj
  real (kind=8):: twothird

  !allocate(temp(nmax,nmax,nvars))

  twothird = 2.0d0*onethird
  call sutherland
  !print*,'twothird = ',twothird, 'max(mu) =', maxval(xmu), minval(xmu)
  call xivcomp
  call etvcomp
  call ztvcomp

  !call crossvcomp

  deallocate(temp, stat =  ierr1)

contains
  subroutine xivcomp
    ! In the energy equation, k = mu/pr, taken care by FACTOR. Check again
    do j = 2,ny-1
      do k = 2,nz-1
        do i = 1,nx
          tdg(k,i,1) = prim(i,j,k,2)
          tdg(k,i,2) = prim(i,j,k,3)
          tdg(k,i,3) = prim(i,j,k,4)
          tdg(k,i,4) = T(i,j,k)
        enddo
      enddo

      call CompactDeriv(4,nx,iper,2,nz-1)

      do k = 2,nz-1
        do i = 1,nx
          xmuxj = xmu(i,j,k)
          derx2 = xix(i,j,k)**2
          dery2 = xiy(i,j,k)**2
          derz2 = xiz(i,j,k)**2
          dertot = derx2 + dery2 +derz2
          derxy = xix(i,j,k)*xiy(i,j,k)
          deryz = xiz(i,j,k)*xiy(i,j,k)
          derzx = xix(i,j,k)*xiz(i,j,k)

          uder = tdg(k,i,1)
          vder = tdg(k,i,2)
          wder = tdg(k,i,3)
          Tder = tdg(k,i,4)

          tdg(k,i,1) = xmuxj*((dertot + onethird*derx2)*uder &
            & + onethird*(derxy*vder+derzx*wder))
          tdg(k,i,2) = xmuxj*((dertot + onethird*dery2)*vder &
            & + onethird*(derxy*uder+deryz*wder))
          tdg(k,i,3) = xmuxj*((dertot + onethird*derz2)*wder &
            & + onethird*(derzx*uder+deryz*vder))
          tdg(k,i,4) = prim(i,j,k,2)*tdg(k,i,1) + prim(i,j,k,3)*tdg(k,i,2)  &
            & + prim(i,j,k,4)*tdg(k,i,3) + xmuxj*FACTOR*dertot*Tder
        enddo
      enddo

      call CompactDeriv(4,nx,iper,2,nz-1)

      do ivar =2,nvars
        do k = 2,nz-1
          do i =2,nx-1
            rhs(i,j,k,ivar) = rhs(i,j,k,ivar) - REI*tdg(k,i,ivar-1)
          enddo
        enddo
        !print*,'after', j, maxval(tdg),minval(tdg), maxval(rhs), maxval(xmu)
      enddo
    enddo
  end subroutine xivcomp

  subroutine etvcomp
    ! In the energy equation, k = mu/pr, taken care by FACTOR. Check again
    do k = 2,nz-1
      do i = 2,nx-1
        do j = 1,ny
          tdg(i,j,1) = prim(i,j,k,2)
          tdg(i,j,2) = prim(i,j,k,3)
          tdg(i,j,3) = prim(i,j,k,4)
          tdg(i,j,4) = T(i,j,k)
        enddo
      enddo

      call CompactDeriv(4,ny,jper,2,nx-1)

      do i = 2,nx-1
        do j = 1,ny
          xmuxj = xmu(i,j,k)
          derx2 = etax(i,j,k)**2
          dery2 = etay(i,j,k)**2
          derz2 = etaz(i,j,k)**2
          dertot = derx2 + dery2 +derz2
          derxy = etax(i,j,k)*etay(i,j,k)
          deryz = etaz(i,j,k)*etay(i,j,k)
          derzx = etax(i,j,k)*etaz(i,j,k)

          uder = tdg(i,j,1)
          vder = tdg(i,j,2)
          wder = tdg(i,j,3)
          Tder = tdg(i,j,4)

          tdg(i,j,1) = xmuxj*((dertot + onethird*derx2)*uder &
            & + onethird*(derxy*vder+derzx*wder))
          tdg(i,j,2) = xmuxj*((dertot + onethird*dery2)*vder &
            & + onethird*(derxy*uder+deryz*wder))
          tdg(i,j,3) = xmuxj*((dertot + onethird*derz2)*wder &
            & + onethird*(derzx*uder+deryz*vder))
          tdg(i,j,4) = prim(i,j,k,2)*tdg(i,j,1) + prim(i,j,k,3)*tdg(i,j,2)  &
            & + prim(i,j,k,4)*tdg(i,j,3) + xmuxj*FACTOR*dertot*Tder
        enddo
      enddo

      call CompactDeriv(4,ny,jper,2,nx-1)

      do ivar =2,nvars
        do i =2,nx-1
          do j = 1,ny
            rhs(i,j,k,ivar) = rhs(i,j,k,ivar) - REI*tdg(i,j,ivar-1)
          enddo
        enddo
        !print*,'after', j, maxval(tdg),minval(tdg), maxval(rhs), maxval(xmu)
      enddo
    enddo
  end subroutine etvcomp

  subroutine ztvcomp
    ! In the energy equation, k = mu/pr, taken care by FACTOR. Check again
    do j = 2,ny-1
      do i = 2,nx-1
        do k = 1,nz
          tdg(i,k,1) = prim(i,j,k,2)
          tdg(i,k,2) = prim(i,j,k,3)
          tdg(i,k,3) = prim(i,j,k,4)
          tdg(i,k,4) = T(i,j,k)
        enddo
      enddo

      call CompactDeriv(4,nz,kper,2,nx-1)

      do i = 2,nx-1
        do k = 1,nz
          xmuxj = xmu(i,j,k)
          derx2 = ztax(i,j,k)**2
          dery2 = ztay(i,j,k)**2
          derz2 = ztaz(i,j,k)**2
          dertot = derx2 + dery2 +derz2
          derxy = ztax(i,j,k)*ztay(i,j,k)
          deryz = ztaz(i,j,k)*ztay(i,j,k)
          derzx = ztax(i,j,k)*ztaz(i,j,k)

          uder = tdg(i,k,1)
          vder = tdg(i,k,2)
          wder = tdg(i,k,3)
          Tder = tdg(i,k,4)

          tdg(i,k,1) = xmuxj*((dertot + onethird*derx2)*uder &
            & + onethird*(derxy*vder+derzx*wder))
          tdg(i,k,2) = xmuxj*((dertot + onethird*dery2)*vder &
            & + onethird*(derxy*uder+deryz*wder))
          tdg(i,k,3) = xmuxj*((dertot + onethird*derz2)*wder &
            & + onethird*(derzx*uder+deryz*vder))
          tdg(i,k,4) = prim(i,j,k,2)*tdg(i,k,1) + prim(i,j,k,3)*tdg(i,k,2)  &
            & + prim(i,j,k,4)*tdg(i,k,3) + xmuxj*FACTOR*dertot*Tder
        enddo
      enddo

      call CompactDeriv(4,nz,kper,2,nx-1)

      do ivar =2,nvars
        do i =2,nx-1
          do k = 1,nz
            rhs(i,j,k,ivar) = rhs(i,j,k,ivar) - REI*tdg(i,k,ivar-1)
          enddo
        enddo
        !print*,'after', j, maxval(tdg),minval(tdg), maxval(rhs), maxval(xmu)
      enddo
    enddo
  end subroutine ztvcomp
end subroutine CompactVflux
