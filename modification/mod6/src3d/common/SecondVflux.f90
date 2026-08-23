subroutine SecondVflux
  use mod_params
  use nvtx
  implicit none

  integer :: i,j,k,ivar
  integer :: jp, jm, kp, km
  integer :: viscross
  real (kind=8):: uder, vder,wder, Tder
  real (kind=8):: derx2, dery2, derz2, derxy, deryz, derzx, dertot
  real (kind=8):: xmuxj
  real (kind=8):: twothird
  real (kind=8), dimension(4):: tmp1, tmp2, tmp3

  viscross = 0
  !allocate(temp(nmax,nmax,nmax,nvars))
  twothird = 2.0d0*onethird
  
  call nvtxStartRange("sutherland")
  call sutherland
  call nvtxEndRange
  !print*,'twothird = ',twothird, 'max(mu) =', maxval(xmu), minval(xmu)
  
  call nvtxStartRange("xivsecond")
  call xivsecond
  call nvtxEndRange
  
  call nvtxStartRange("etvsecond")
  call etvsecond
  call nvtxEndRange
  
  call nvtxStartRange("ztvsecond")
  call ztvsecond
  call nvtxEndRange
  
  if(viscross .eq. 1)then
    allocate(temp2(nmax,nmax,nmax,nvars))
    call nvtxStartRange("viscrossz")
    call viscrossz !takes case of cross-terms in et
    call nvtxEndRange
  
    call nvtxStartRange("viscrossy")
    call viscrossy !takes case of cross-terms in xi and zt
    call nvtxEndRange
    deallocate(temp2, stat =  ierr1)
  endif

  !deallocate(temp, stat =  ierr1)

  contains
    subroutine xivsecond
      !! dF_v/dxi
      ! In the energy equation, k = mu/pr, taken care by FACTOR. Check again
      ! j=constant plane
      
        !!!$acc& private(xmuxj,derx2,dery2,derz2,dertot,derxy,deryz,derzx,uder,vder,wder,Tder)
      
      !$acc parallel loop collapse(3) &
      !$acc& present(prim, rhs, ifblock, xix, xiy, xiz) &
      !$acc& private(xmuxj,derx2,dery2,derz2,dertot,derxy,deryz,derzx,uder,vder,wder,Tder)
      do k = 2,nz-1
        do j = 2,ny-1
          do i = 2,nx
            !xmuxj = xmu(i,j,k)/xjac(i,j,k)
            xmuxj = xmu(i,j,k)
            derx2 = xix(i,j,k)**2
            dery2 = xiy(i,j,k)**2
            derz2 = xiz(i,j,k)**2
            dertot = derx2 + dery2 +derz2
            derxy = xix(i,j,k)*xiy(i,j,k)
            deryz = xiz(i,j,k)*xiy(i,j,k)
            derzx = xix(i,j,k)*xiz(i,j,k)
            uder = prim(i,j,k,2) - prim(i-1,j,k,2)
            vder = prim(i,j,k,3) - prim(i-1,j,k,3)
            wder = prim(i,j,k,4) - prim(i-1,j,k,4)
            Tder = GXM*(prim(i,j,k,5)/prim(i,j,k,1) - prim(i-1,j,k,5)/prim(i-1,j,k,1))

            temp(i,j,k,1) = xmuxj*((dertot + onethird*derx2)*uder &
              & + onethird*(derxy*vder+derzx*wder))
            temp(i,j,k,2) = xmuxj*((dertot + onethird*dery2)*vder &
              & + onethird*(derxy*uder+deryz*wder))
            temp(i,j,k,3) = xmuxj*((dertot + onethird*derz2)*wder &
              & + onethird*(derzx*uder+deryz*vder))
            temp(i,j,k,4) = prim(i,j,k,2)*temp(i,j,k,1) + prim(i,j,k,3)*temp(i,j,k,2)  &
              & + prim(i,j,k,4)*temp(i,j,k,3) + xmuxj*FACTOR*dertot*Tder
          enddo
        enddo
      enddo

      !$acc parallel loop collapse(4) &
      !$acc& present(rhs, ifblock)
      do ivar =2,nvars
        do k = 2,nz-1
          do j = 2,ny-1
            do i =2,nx-1
              if (ifblock(i,j,k) .eq. 0)then
                rhs(i,j,k,ivar) = rhs(i,j,k,ivar) - REI*(temp(i+1,j,k,ivar-1)-temp(i,j,k,ivar-1))
              endif
            enddo
          enddo
        enddo
      enddo

      !cross terms are taken care in viscrossy and viscrossz
    end subroutine xivsecond

    subroutine etvsecond
      !! dG_v/deta
      ! In the energy equation, ivar = mu/pr, taken care by FACTOR. Checivar again
      !k=constant plane
      
      !$acc parallel loop collapse(3) &
      !$acc& present(prim, rhs, ifblock, etax, etay, etaz) &
      !$acc& private(xmuxj,derx2,dery2,derz2,dertot,derxy,deryz,derzx,uder,vder,wder,Tder)
      do k = 2,nz-1
        do j = 2,ny
          do i = 2,nx-1
            !xmuxj = xmu(i,j,k)/xjac(i,j,k)
            xmuxj = xmu(i,j,k)
            derx2 = etax(i,j,k)**2
            dery2 = etay(i,j,k)**2
            derz2 = etaz(i,j,k)**2
            dertot = derx2 + dery2 +derz2
            derxy = etax(i,j,k)*etay(i,j,k)
            deryz = etay(i,j,k)*etaz(i,j,k)
            derzx = etaz(i,j,k)*etax(i,j,k)
            uder = prim(i,j,k,2) - prim(i,j-1,k,2)
            vder = prim(i,j,k,3) - prim(i,j-1,k,3)
            wder = prim(i,j,k,4) - prim(i,j-1,k,4)
            Tder = GXM*(prim(i,j,k,5)/prim(i,j,k,1) - prim(i,j-1,k,5)/prim(i,j-1,k,1))

            temp(i,j,k,1) = xmuxj*((dertot + onethird*derx2)*uder &
              & + onethird*(derxy*vder+derzx*wder))
            temp(i,j,k,2) = xmuxj*((dertot + onethird*dery2)*vder &
              & + onethird*(derxy*uder+deryz*wder))
            temp(i,j,k,3) = xmuxj*((dertot + onethird*derz2)*wder &
              & + onethird*(derzx*uder+deryz*vder))
            temp(i,j,k,4) = prim(i,j,k,2)*temp(i,j,k,1) + prim(i,j,k,3)*temp(i,j,k,2)  &
              & + prim(i,j,k,4)*temp(i,j,k,3) + xmuxj*FACTOR*dertot*Tder
          enddo
        enddo
      enddo

      !$acc parallel loop collapse(4) &
      !$acc& present(rhs, ifblock)
      do ivar =2,nvars
        do k = 2,nz-1
          do j =2,ny-1
            do i = 2,nx-1
              if (ifblock(i,j,k) .eq. 0)then
                rhs(i,j,k,ivar) = rhs(i,j,k,ivar) - REI*(temp(i,j+1,k,ivar-1)-temp(i,j,k,ivar-1))
              endif
            enddo
          enddo
        enddo
      enddo
      
      !    call viscrossz(k)
    end subroutine etvsecond

    subroutine ztvsecond
      !! dG_v/deta
      ! In the energy equation, ivar = mu/pr, taivaren care by FACTOR. Checivar again
      ! j=constant plane
      
      !$acc parallel loop collapse(3) &
      !$acc& present(prim, rhs, ifblock, ztax, ztay, ztaz) &
      !$acc& private(xmuxj,derx2,dery2,derz2,dertot,derxy,deryz,derzx,uder,vder,wder,Tder)
      do k = 2,nz
        do j = 2,ny-1
          do i = 2,nx-1
            !xmuxj = xmu(i,j,k)/xjac(i,j,k)
            xmuxj = xmu(i,j,k)
            derx2 = ztax(i,j,k)**2
            dery2 = ztay(i,j,k)**2
            derz2 = ztaz(i,j,k)**2
            dertot = derx2 + dery2 +derz2
            derxy = ztax(i,j,k)*ztay(i,j,k)
            deryz = ztay(i,j,k)*ztaz(i,j,k)
            derzx = ztaz(i,j,k)*ztax(i,j,k)
            uder = prim(i,j,k,2) - prim(i,j,k-1,2)
            vder = prim(i,j,k,3) - prim(i,j,k-1,3)
            wder = prim(i,j,k,4) - prim(i,j,k-1,4)
            Tder = GXM*(prim(i,j,k,5)/prim(i,j,k,1) - prim(i,j,k,5)/prim(i,j,k-1,1))

            temp(i,j,k,1) = xmuxj*((dertot + onethird*derx2)*uder &
              & + onethird*(derxy*vder+derzx*wder))
            temp(i,j,k,2) = xmuxj*((dertot + onethird*dery2)*vder &
              & + onethird*(derxy*uder+deryz*wder))
            temp(i,j,k,3) = xmuxj*((dertot + onethird*derz2)*wder &
              & + onethird*(derzx*uder+deryz*vder))
            temp(i,j,k,4) = prim(i,j,k,2)*temp(i,j,k,1) + prim(i,j,k,3)*temp(i,j,k,2)  &
              & + prim(i,j,k,4)*temp(i,j,k,3) + xmuxj*FACTOR*dertot*Tder
          enddo
        enddo
      enddo

      !$acc parallel loop collapse(4) &
      !$acc& present(rhs, ifblock)
      do ivar =2,nvars
        do k =2,nz-1
          do j = 2,ny-1
            do i = 2,nx-1
              if (ifblock(i,j,k) .eq. 0)then
                rhs(i,j,k,ivar) = rhs(i,j,k,ivar) - REI*(temp(i,j,k+1,ivar-1)-temp(i,j,k,ivar-1))
              endif
            enddo
          enddo
        enddo
      enddo
      
      !    call viscrossy(j)
    end subroutine ztvsecond


    subroutine viscrossz
      !compute cross derivative on k constant plane
      do k = 2,nz-1
        do i = 2,nx-1
          do j = 2,ny
            !xmuxj = xmu(i,j,k)/xjac(i,j,k)
            kp = k+1
            km = k-1
            xmuxj = xmu(i,j,k)

            !compute d/dzeta
            uder = prim(i,j,kp,2) - prim(i,j,km,2)
            vder = prim(i,j,kp,3) - prim(i,j,km,3)
            wder = prim(i,j,kp,4) - prim(i,j,km,4)
            Tder = GXM*(prim(i,j,kp,5)/prim(i,j,kp,1) - prim(i,j,km,5)/prim(i,j,km,1))
            dertot = uder*ztax(i,j,k) + vder*ztay(i,j,k) + wder*ztaz(i,j,k) 


            !---- FV(zta) ---------- 
            tmp1(1) = xmuxj*(2.0d0*uder*ztax(i,j,k) - twothird*dertot) 
            tmp1(2) = xmuxj*(uder*ztay(i,j,k) + vder*ztax(i,j,k)) 
            tmp1(3) = xmuxj*(uder*ztaz(i,j,k) + wder*ztax(i,j,k)) 
            tmp1(4) = prim(i,j,k,2)*tmp1(1) + prim(i,j,k,3)*tmp1(2)  &
              & + prim(i,j,k,4)*tmp1(3) + xmuxj*FACTOR*Tder*ztax(i,j,k)

            !---- GV(zta) ---------- 
            tmp2(1) = tmp1(2) 
            tmp2(2) = xmuxj*(2.0d0*vder*ztay(i,j,k) - twothird*dertot) 
            tmp2(3) = xmuxj*(vder*ztaz(i,j,k) + wder*ztay(i,j,k)) 
            tmp2(4) = prim(i,j,k,2)*tmp2(1) + prim(i,j,k,3)*tmp2(2)  &
              & + prim(i,j,k,4)*tmp2(3) + xmuxj*FACTOR*Tder*ztay(i,j,k)

            !---- HV(zta) ---------- 
            tmp3(1) = tmp1(3) 
            tmp3(2) = tmp2(3) 
            tmp3(3) = xmuxj*(2.0d0*wder*ztaz(i,j,k) - twothird*dertot) 
            tmp3(4) = prim(i,j,k,2)*tmp3(1) + prim(i,j,k,3)*tmp3(2)  &
              & + prim(i,j,k,4)*tmp3(3) + xmuxj*FACTOR*Tder*ztaz(i,j,k)

            do ivar = 1,nvars-1
              temp(i,j,k,ivar) = xix(i,j,k)*tmp1(ivar) + xiy(i,j,k)*tmp2(ivar) + &
                xiz(i,j,k)*tmp3(ivar) 
              temp2(i,j,k,ivar) = etax(i,j,k)*tmp1(ivar) + etay(i,j,k)*tmp2(ivar) + &
                etaz(i,j,k)*tmp3(ivar) 
            enddo
          enddo
        enddo


        do ivar =2,nvars
          do i = 2,nx-1
            do j =2,ny-1
              if (ifblock(i,j,k) .eq. 0)then
                rhs(i,j,k,ivar) = rhs(i,j,k,ivar) - REI*((temp(i+1,j,k,ivar-1)-temp(i-1,j,k,ivar-1))+&
                  temp2(i,j+1,k,ivar-1)-temp2(i,j-1,k,ivar-1))
              endif
            enddo
          enddo
        enddo
      enddo !end k loop
    end subroutine viscrossz


    subroutine viscrossy
      !compute cross derivative on j constant plane
      do j = 2,ny-1
        do k = 2,nz-1
          do i = 2,nx-1
            !xmuxj = xmu(i,j,k)/xjac(i,j,k)
            jp = j+1
            jm = j-1
            xmuxj = xmu(i,j,k)

            !compute d/dzeta
            uder = prim(i,jp,k,2) - prim(i,jm,k,2)
            vder = prim(i,jp,k,3) - prim(i,jm,k,3)
            wder = prim(i,jp,k,4) - prim(i,jm,k,4)
            Tder = GXM*(prim(i,jp,k,5)/prim(i,jp,k,1) - prim(i,jm,k,5)/prim(i,jm,k,1))
            dertot = uder*ztax(i,j,k) + vder*ztay(i,j,k) + wder*ztaz(i,j,k) 


            !---- FV(zta) ---------- 
            tmp1(1) = xmuxj*(2.0d0*uder*etax(i,j,k) - twothird*dertot) 
            tmp1(2) = xmuxj*(uder*etay(i,j,k) + vder*etax(i,j,k)) 
            tmp1(3) = xmuxj*(uder*etaz(i,j,k) + wder*etax(i,j,k)) 
            tmp1(4) = prim(i,j,k,2)*tmp1(1) + prim(i,j,k,3)*tmp1(2)  &
              & + prim(i,j,k,4)*tmp1(3) + xmuxj*FACTOR*Tder*etax(i,j,k)

            !---- GV(zta) ---------- 
            tmp2(1) = tmp1(2) 
            tmp2(2) = xmuxj*(2.0d0*vder*etay(i,j,k) - twothird*dertot) 
            tmp2(3) = xmuxj*(vder*etaz(i,j,k) + wder*etay(i,j,k)) 
            tmp2(4) = prim(i,j,k,2)*tmp2(1) + prim(i,j,k,3)*tmp2(2)  &
              & + prim(i,j,k,4)*tmp2(3) + xmuxj*FACTOR*Tder*etay(i,j,k)

            !---- HV(zta) ---------- 
            tmp3(1) = tmp1(3) 
            tmp3(2) = tmp2(3) 
            tmp3(3) = xmuxj*(2.0d0*wder*etaz(i,j,k) - twothird*dertot) 
            tmp3(4) = prim(i,j,k,2)*tmp3(1) + prim(i,j,k,3)*tmp3(2)  &
              & + prim(i,j,k,4)*tmp3(3) + xmuxj*FACTOR*Tder*etaz(i,j,k)

            do ivar = 1,nvars-1
              temp(i,j,k,ivar) = xix(i,j,k)*tmp1(ivar) + xiy(i,j,k)*tmp2(ivar) + &
                xiz(i,j,k)*tmp3(ivar) 
              temp2(i,j,k,ivar) = ztax(i,j,k)*tmp1(ivar) + ztay(i,j,k)*tmp2(ivar) + &
                ztaz(i,j,k)*tmp3(ivar) 
            enddo
          enddo
        enddo
        do ivar =2,nvars
          do k =2,nz-1
            do i = 2,nx-1
              if (ifblock(i,j,k) .eq. 0)then
                rhs(i,j,k,ivar) = rhs(i,j,k,ivar) - REI*((temp(i+1,j,k,ivar-1)-temp(i-1,j,k,ivar-1))+&
                  temp2(i,j,k+1,ivar-1)-temp2(i,j,k-1,ivar-1))
              endif
            enddo
          enddo
        enddo
      enddo !end j loop
    end subroutine viscrossy
  end subroutine SecondVflux
