subroutine FirstVflux
  use mod_params
  implicit none
  integer :: i,j,k,ivar
  real(kind=8):: fluxi, fluxj, fluxk
  real(kind=8):: tauxx, tauyy, tauzz
  real(kind=8):: tauxy, tauyz, tauzx
  real(kind=8):: uderi, vderi, wderi, Tderi 
  real(kind=8):: uderj, vderj, wderj, Tderj 
  real(kind=8):: uderk, vderk, wderk, Tderk 
  real(kind=8):: uhf, vhf, whf,xmuhf
  real(kind=8):: dudx, dudy, dvdx, dvdy, dwdx, dwdy
  real(kind=8):: dudz, dvdz, dwdz, divU
  real(kind=8):: dTdx, dTdy, dTdz
  real(kind=8):: xixm, xiym, xizm
  real(kind=8):: etaxm, etaym, etazm
  real(kind=8):: ztaxm, ztaym, ztazm
  real(kind=8):: xmuxj

  !integer:: ivscs,jvscs,kvscs,ip1,im1,jp1,jm1
  !double precision::flgi,flgj,tv1,tv2

  xiflux = 0.0d0
  etflux = 0.0d0
  ztflux = 0.0d0

  call sutherland
  call xiv
  call etv
  call ztv

  !Take flux difference
  if (ifblock(i,j,k) .eq. 0)then
    do ivar=1,nvars
      do k=2,nz-1
        do j=2,ny-1
          do i=2,nx-1
            fluxi=(xiflux(i,j,k,ivar)-xiflux(i-1,j,k,ivar))
            fluxj=(etflux(i,j,k,ivar)-etflux(i,j-1,k,ivar))
            fluxk=(ztflux(i,j,k,ivar)-ztflux(i,j,k-1,ivar))
            rhs(i,j,k,ivar) = rhs(i,j,k,ivar) - REI*(fluxi+fluxj+fluxk)/xjac(i,j,k)
          enddo  
        enddo
      enddo
    enddo
  endif

  contains
  !! First xiflux
subroutine xiv
  do k=2,nz-1
    do j=2,ny-1
      do i=1,nx-1
        !xjacmi = hf*(xjac(i,j,k)+xjac(i+1,j,k))
        !xjaci = 1.0d0/xjacmi

        xixm =hf*(xix(i,j,k)*xjac(i,j,k)+xix(i+1,j,k)*xjac(i+1,j,k)) 
        xiym =hf*(xiy(i,j,k)*xjac(i,j,k)+xiy(i+1,j,k)*xjac(i+1,j,k)) 
        xizm = hf*(xiz(i,j,k)*xjac(i,j,k)+xiz(i+1,j,k)*xjac(i+1,j,k)) 
        etaxm =hf*(etax(i,j,k)*xjac(i,j,k)+etax(i+1,j,k)*xjac(i+1,j,k)) 
        etaym =hf*(etay(i,j,k)*xjac(i,j,k)+etay(i+1,j,k)*xjac(i+1,j,k)) 
        etazm =hf*(etaz(i,j,k)*xjac(i,j,k)+etaz(i+1,j,k)*xjac(i+1,j,k)) 
        ztaxm = hf*(ztax(i,j,k)*xjac(i,j,k)+ztax(i+1,j,k)*xjac(i+1,j,k)) 
        ztaym =hf*(ztay(i,j,k)*xjac(i,j,k)+ztay(i+1,j,k)*xjac(i+1,j,k)) 
        ztazm = hf*(ztaz(i,j,k)*xjac(i,j,k)+ztaz(i+1,j,k)*xjac(i+1,j,k)) 

        uhf=(prim(i+1,j,k,2)+prim(i,j,k,2))*hf
        vhf=(prim(i+1,j,k,3)+prim(i,j,k,3))*hf
        whf=(prim(i+1,j,k,4)+prim(i,j,k,4))*hf
        xmuhf=(xmu(i+1,j,k)+xmu(i,j,k))*hf
        !xmuxj = xmuhf*xjacmi
        xmuxj = xmuhf

        uderi=prim(i+1,j,k,2)-prim(i,j,k,2)
        vderi=prim(i+1,j,k,3)-prim(i,j,k,3)
        wderi=prim(i+1,j,k,4)-prim(i,j,k,4)
        Tderi = T(i+1,j,k) - T(i,j,k)

        uderj=(prim(i,j+1,k,2)+prim(i+1,j+1,k,2)&
          &-prim(i+1,j-1,k,2)-prim(i,j-1,k,2))*fth
        vderj=(prim(i,j+1,k,3)+prim(i+1,j+1,k,3)&
          &-prim(i+1,j-1,k,3)-prim(i,j-1,k,3))*fth
        wderj=(prim(i,j+1,k,4)+prim(i+1,j+1,k,4)&
          &-prim(i+1,j-1,k,4)-prim(i,j-1,k,4))*fth
        Tderj=(T(i,j+1,k)+T(i+1,j+1,k)&
          &-T(i+1,j-1,k)-T(i,j-1,k))*fth

        uderk=(prim(i,j,k+1,2)+prim(i+1,j,k+1,2)&
          &-prim(i+1,j,k-1,2)-prim(i,j,k-1,2))*fth
        vderk=(prim(i,j,k+1,3)+prim(i+1,j,k+1,3)&
          &-prim(i+1,j,k-1,3)-prim(i,j,k-1,3))*fth
        wderk=(prim(i,j,k+1,4)+prim(i+1,j,k+1,4)&
          &-prim(i+1,j,k-1,4)-prim(i,j,k-1,4))*fth
        Tderk=(T(i,j,k+1)+T(i+1,j,k+1)&
          &-T(i+1,j,k-1)-T(i,j,k-1))*fth

        dudx=xixm*uderi+etaxm*uderj+ztaxm*uderk
        dvdx=xixm*vderi+etaxm*vderj+ztaxm*vderk
        dwdx=xixm*wderi+etaxm*wderj+ztaxm*wderk
        dTdx=xixm*Tderi+etaxm*Tderj+ztaxm*Tderk

        dudy=xiym*uderi+etaym*uderj+ztaym*uderk
        dvdy=xiym*vderi+etaym*vderj+ztaym*vderk
        dwdy=xiym*wderi+etaym*wderj+ztaym*wderk
        dTdy=xiym*Tderi+etaym*Tderj+ztaym*Tderk

        dudz=xizm*uderi+etazm*uderj+ztazm*uderk
        dvdz=xizm*vderi+etazm*vderj+ztazm*vderk
        dwdz=xizm*wderi+etazm*wderj+ztazm*wderk
        dTdz=xizm*Tderi+etazm*Tderj+ztazm*Tderk

        divU = dudx + dvdy + dwdz
        tauxx=(2.0d0*dudx-(2.0d0/3.0d0)*divU)
        tauyy=(2.0d0*dvdy-(2.0d0/3.0d0)*divU)
        tauzz=(2.0d0*dwdz-(2.0d0/3.0d0)*divU)

        tauxy=(dudy+dvdx)
        tauyz=(dvdz+dwdy)
        tauzx=(dwdx+dudz)

        xiflux(i,j,k,1)=0.0d0
        xiflux(i,j,k,2)=(xixm*tauxx+xiym*tauxy+xizm*tauzx)*xmuxj
        xiflux(i,j,k,3)=(xixm*tauxy+xiym*tauyy+xizm*tauyz)*xmuxj
        xiflux(i,j,k,4)=(xixm*tauzx+xiym*tauyz+xizm*tauzz)*xmuxj
        xiflux(i,j,k,5)=(uhf*xiflux(i,j,k,2)+vhf*xiflux(i,j,k,3)+whf*xiflux(i,j,k,4))&
          &+FACTOR*xmuxj*(xixm*dTdx+xiym*dTdy+xizm*dTdz)
      enddo
    enddo
  enddo
end subroutine xiv


subroutine etv
  !! etflux
  do k=2,nz-1
    do i=2,nx-1
      do j=1,ny-1
        xixm = hf*(xix(i,j,k)*xjac(i,j,k)+xix(i,j+1,k)*xjac(i,j+1,k)) 
        xiym = hf*(xiy(i,j,k)*xjac(i,j,k)+xiy(i,j+1,k)*xjac(i,j+1,k)) 
        xizm = hf*(xiz(i,j,k)*xjac(i,j,k)+xiz(i,j+1,k)*xjac(i,j+1,k)) 
        etaxm = hf*(etax(i,j,k)*xjac(i,j,k)+etax(i,j+1,k)*xjac(i,j+1,k)) 
        etaym = hf*(etay(i,j,k)*xjac(i,j,k)+etay(i,j+1,k)*xjac(i,j+1,k)) 
        etazm = hf*(etaz(i,j,k)*xjac(i,j,k)+etaz(i,j+1,k)*xjac(i,j+1,k)) 
        ztaxm = hf*(ztax(i,j,k)*xjac(i,j,k)+ztax(i,j+1,k)*xjac(i,j+1,k)) 
        ztaym = hf*(ztay(i,j,k)*xjac(i,j,k)+ztay(i,j+1,k)*xjac(i,j+1,k)) 
        ztazm = hf*(ztaz(i,j,k)*xjac(i,j,k)+ztaz(i,j+1,k)*xjac(i,j+1,k)) 
        uhf=(prim(i,j+1,k,2)+prim(i,j,k,2))*hf
        vhf=(prim(i,j+1,k,3)+prim(i,j,k,3))*hf
        whf=(prim(i,j+1,k,4)+prim(i,j,k,4))*hf
        xmuhf=(xmu(i,j+1,k)+xmu(i,j,k))*hf
        xmuxj = xmuhf

        uderi=(prim(i+1,j,k,2)+prim(i+1,j+1,k,2)&
          &-prim(i-1,j,k,2)-prim(i-1,j+1,k,2))*fth
        vderi=(prim(i+1,j,k,3)+prim(i+1,j+1,k,3)-&
          &prim(i-1,j,k,3)-prim(i-1,j+1,k,3))*fth
        wderi=(prim(i+1,j,k,4)+prim(i+1,j+1,k,4)-&
          &prim(i-1,j,k,4)-prim(i-1,j+1,k,4))*fth
        Tderi=(T(i+1,j,k)+T(i+1,j+1,k)&
          &-T(i-1,j,k)-T(i-1,j+1,k))*fth

        uderj=prim(i,j+1,k,2)-prim(i,j,k,2)
        vderj=prim(i,j+1,k,3)-prim(i,j,k,3)
        wderj=prim(i,j+1,k,4)-prim(i,j,k,4)
        Tderj=T(i,j+1,k)-T(i,j,k)

        uderk=(prim(i,j,k+1,2)+prim(i,j+1,k+1,2)&
          &-prim(i,j,k-1,2)-prim(i,j+1,k-1,2))*fth
        vderk=(prim(i,j,k+1,3)+prim(i,j+1,k+1,3)&
          &-prim(i,j,k-1,3)-prim(i,j+1,k-1,3))*fth
        wderk=(prim(i,j,k+1,4)+prim(i,j+1,k+1,4)&
          &-prim(i,j,k-1,4)-prim(i,j+1,k-1,4))*fth
        Tderk=(T(i,j,k+1)+T(i,j+1,k+1)&
          &-T(i,j,k-1)-T(i,j+1,k-1))*fth

        dudx=xixm*uderi+etaxm*uderj+ztaxm*uderk
        dvdx=xixm*vderi+etaxm*vderj+ztaxm*vderk
        dwdx=xixm*wderi+etaxm*wderj+ztaxm*wderk
        dTdx=xixm*Tderi+etaxm*Tderj+ztaxm*Tderk

        dudy=xiym*uderi+etaym*uderj+ztaym*uderk
        dvdy=xiym*vderi+etaym*vderj+ztaym*vderk
        dwdy=xiym*wderi+etaym*wderj+ztaym*wderk
        dTdy=xiym*Tderi+etaym*Tderj+ztaym*Tderk

        dudz=xizm*uderi+etazm*uderj+ztazm*uderk
        dvdz=xizm*vderi+etazm*vderj+ztazm*vderk
        dwdz=xizm*wderi+etazm*wderj+ztazm*wderk
        dTdz=xizm*Tderi+etazm*Tderj+ztazm*Tderk


        divU = dudx + dvdy + dwdz
        tauxx=(2.0d0*dudx-(2.0d0/3.0d0)*divU)
        tauyy=(2.0d0*dvdy-(2.0d0/3.0d0)*divU)
        tauzz=(2.0d0*dwdz-(2.0d0/3.0d0)*divU)

        tauxy=(dudy+dvdx)
        tauyz=(dvdz+dwdy)
        tauzx=(dwdx+dudz)

        etflux(i,j,k,1)=0.0d0
        etflux(i,j,k,2)=(etaxm*tauxx+etaym*tauxy+etazm*tauzx)*xmuxj
        etflux(i,j,k,3)=(etaxm*tauxy+etaym*tauyy+etazm*tauyz)*xmuxj
        etflux(i,j,k,4)=(etaxm*tauzx+etaym*tauyz+etazm*tauzz)*xmuxj
        etflux(i,j,k,5)=(uhf*etflux(i,j,k,2)+vhf*etflux(i,j,k,3)+whf*etflux(i,j,k,4))&
          &+FACTOR*xmuxj*(etaxm*dTdx+etaym*dTdy+etazm*dTdz)
      enddo
    enddo
  enddo
end subroutine etv

subroutine ztv
  !! ztflux
  do j=2,ny-1
    do i=2,nx-1
      do k=1,nz-1
        xixm = hf*(xix(i,j,k)*xjac(i,j,k)+xix(i,j,k+1)*xjac(i,j,k+1)) 
        xiym = hf*(xiy(i,j,k)*xjac(i,j,k)+xiy(i,j,k+1)*xjac(i,j,k+1)) 
        xizm = hf*(xiz(i,j,k)*xjac(i,j,k)+xiz(i,j,k+1)*xjac(i,j,k+1)) 
        etaxm = hf*(etax(i,j,k)*xjac(i,j,k)+etax(i,j,k+1)*xjac(i,j,k+1)) 
        etaym = hf*(etay(i,j,k)*xjac(i,j,k)+etay(i,j,k+1)*xjac(i,j,k+1)) 
        etazm = hf*(etaz(i,j,k)*xjac(i,j,k)+etaz(i,j,k+1)*xjac(i,j,k+1)) 
        ztaxm = hf*(ztax(i,j,k)*xjac(i,j,k)+ztax(i,j,k+1)*xjac(i,j,k+1)) 
        ztaym = hf*(ztay(i,j,k)*xjac(i,j,k)+ztay(i,j,k+1)*xjac(i,j,k+1)) 
        ztazm = hf*(ztaz(i,j,k)*xjac(i,j,k)+ztaz(i,j,k+1)*xjac(i,j,k+1)) 
        uhf=(prim(i,j,k+1,2)+prim(i,j,k,2))*hf
        vhf=(prim(i,j,k+1,3)+prim(i,j,k,3))*hf
        whf=(prim(i,j,k+1,4)+prim(i,j,k,4))*hf
        xmuhf=(xmu(i,j,k+1)+xmu(i,j,k))*hf
        xmuxj = xmuhf

        uderi=(prim(i+1,j,k,2)+prim(i+1,j,k+1,2)&
          &-prim(i-1,j,k,2)-prim(i-1,j,k+1,2))*fth
        vderi=(prim(i+1,j,k,3)+prim(i+1,j,k+1,3)-&
          &prim(i-1,j,k,3)-prim(i-1,j,k+1,3))*fth
        wderi=(prim(i+1,j,k,4)+prim(i+1,j,k+1,4)-&
          &prim(i-1,j,k,4)-prim(i-1,j,k+1,4))*fth
        Tderi=(T(i+1,j,k)+T(i+1,j,k+1)&
          &-T(i-1,j,k)-T(i-1,j,k+1))*fth

        uderj=(prim(i,j+1,k,2)+prim(i,j+1,k+1,2)&
          &-prim(i,j-1,k,2)-prim(i,j-1,k+1,2))*fth
        vderj=(prim(i,j+1,k,3)+prim(i,j+1,k+1,3)-&
          &prim(i,j-1,k,3)-prim(i,j-1,k+1,3))*fth
        wderj=(prim(i,j+1,k,4)+prim(i,j+1,k+1,4)-&
          &prim(i,j-1,k,4)-prim(i,j-1,k+1,4))*fth
        Tderj=(T(i,j+1,k)+T(i,j+1,k+1)&
          &-T(i,j-1,k)-T(i,j-1,k+1))*fth

        uderk=prim(i,j,k+1,2)-prim(i,j,k,2)
        vderk=prim(i,j,k+1,3)-prim(i,j,k,3)
        wderk=prim(i,j,k+1,4)-prim(i,j,k,4)
        Tderk=T(i,j,k+1)-T(i,j,k)

        dudx=xixm*uderi+etaxm*uderj+ztaxm*uderk
        dvdx=xixm*vderi+etaxm*vderj+ztaxm*vderk
        dwdx=xixm*wderi+etaxm*wderj+ztaxm*wderk
        dTdx=xixm*Tderi+etaxm*Tderj+ztaxm*Tderk

        dudy=xiym*uderi+etaym*uderj+ztaym*uderk
        dvdy=xiym*vderi+etaym*vderj+ztaym*vderk
        dwdy=xiym*wderi+etaym*wderj+ztaym*wderk
        dTdy=xiym*Tderi+etaym*Tderj+ztaym*Tderk

        dudz=xizm*uderi+etazm*uderj+ztazm*uderk
        dvdz=xizm*vderi+etazm*vderj+ztazm*vderk
        dwdz=xizm*wderi+etazm*wderj+ztazm*wderk
        dTdz=xizm*Tderi+etazm*Tderj+ztazm*Tderk


        divU = dudx + dvdy + dwdz
        tauxx=(2.0d0*dudx-(2.0d0/3.0d0)*divU)
        tauyy=(2.0d0*dvdy-(2.0d0/3.0d0)*divU)
        tauzz=(2.0d0*dwdz-(2.0d0/3.0d0)*divU)

        tauxy=(dudy+dvdx)
        tauyz=(dvdz+dwdy)
        tauzx=(dwdx+dudz)

        ztflux(i,j,k,1)=0.0d0
        ztflux(i,j,k,2)=(ztaxm*tauxx+ztaym*tauxy+ztazm*tauzx)*xmuxj
        ztflux(i,j,k,3)=(ztaxm*tauxy+ztaym*tauyy+ztazm*tauyz)*xmuxj
        ztflux(i,j,k,4)=(ztaxm*tauzx+ztaym*tauyz+ztazm*tauzz)*xmuxj
        ztflux(i,j,k,5)=(uhf*ztflux(i,j,k,2)+vhf*ztflux(i,j,k,3)+whf*ztflux(i,j,k,4))&
          &+FACTOR*xmuxj*(ztaxm*dTdx+ztaym*dTdy+ztazm*dTdz)
      enddo
    enddo
  enddo

end subroutine ztv

end subroutine FirstVflux
