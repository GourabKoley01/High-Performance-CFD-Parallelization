subroutine xiroe
  use mod_params
  use nvtx
  implicit none
  !  integer, intent(in) :: idir
  ! real (kind=8),intent(out) :: Roeflux
  real (kind=8) :: &
    &kx,ky,kz,magk,xjacm,xjaci,xjacip,kxn,kyn,kzn,&
    &al,rhol,ul,vl,wl,pl,rhoel,hl,uconl,&
    &ar,rhor,ur,vr,wr,pr,rhoer,hr,uconr,&
    &fql(nvars),fqr(nvars),&
    &sc1,sc2,sc3,&
    &rhoh,uh,vh,wh,hh,ch,&
    &alfa,beta,theta,phi2,& 
    &lamh(nvars,nvars),tk(nvars,nvars),tkInv(nvars,nvars),&
    &ql(nvars),qr(nvars),&
    &dq(nvars),tm1(nvars),tm2(nvars),tm3(nvars)
  real (kind=8) :: epHarten = 0.0010d0
  integer:: i,j,k,ivar,invars
  !This subroutine implements the Roe scheme to determine the flux at an
  !interface given the left and right states. This subroutine is generalized so
  !that it can calculate the flux in any direction given the proper grid
  !metrics as inputs.
  
  call nvtxStartRange("recon(1)")
  call recon(1)
  call nvtxEndRange
  !$acc parallel loop collapse(3) &
  !$acc& present(xix,xiy,xiz,xjac) &
  !$acc& private(fql,fqr,ql,qr,dq,tm1,tm2,tm3,lamh,tk,tkInv) &
  !$acc& private(magk,sc3,alfa,beta,ch,vr,ur,uconr,rhoer,sc2,sc1,rhor,uconl,hl,rhoel,theta,hh,pr,pl,wl,vl,ul,rhol,phi2,wh,vh,uh,rhoh,kz,ky,kx,kzn,kyn,kxn,hr,xjaci,wr,xjacm,xjacip)
  do k=1,nz-1
    do j=1,ny-1
      do i=1,nx-1
        xjacm = hf*(xjac(i,j,k)+xjac(i+1,j,k))
        xjaci = xjac(i,j,k)/xjacm
        xjacip = xjac(i+1,j,k)/xjacm
        kx =hf*(xix(i,j,k)*xjaci+xix(i+1,j,k)*xjacip) 
        ky =hf*(xiy(i,j,k)*xjaci+xiy(i+1,j,k)*xjacip) 
        kz =hf*(xiz(i,j,k)*xjaci+xiz(i+1,j,k)*xjacip) 

        magk = sqrt(kx**2 + ky**2 +kz**2) 
        kxn = kx/magk
        kyn = ky/magk
        kzn = kz/magk

        !!  Left state
        rhol = primL(i,j,k,1) 
        ul = primL(i,j,k,2) 
        vl = primL(i,j,k,3) 
        wl = primL(i,j,k,4) 
        !uconl = ul*kx + vl*ky + wl*kz 
        pl = primL(i,j,k,5) 
        !al = sqrt(GAM*pl/rhol) 
        !rhoel = pl/GAM1 + hf*(ul*ul+vl*vl+wl*wl) 

        !!  Right state
        rhor = primR(i,j,k,1) 
        ur = primR(i,j,k,2) 
        vr = primR(i,j,k,3) 
        wr = primR(i,j,k,4) 
        !uconr = ur*kx + vr*ky + wr*kz 
        pr = primR(i,j,k,5) 
        !ar = sqrt(GAM*pr/rhor) 
        !rhoer = pr/GAM1 + hf*(ur*ur+vr*vr+wr*wr) 

        !!Extract primative variables for left state from ql
        !rhol=ql(1)
        !ul=ql(2)
        !vl=ql(3)
        !wl=ql(4)
        !pl=ql(5)
        !!Calculate rho*e, enthalpy, and contravariant velocity at left state
        rhoel=(pl/GAM1)+0.50d0*rhol*(ul**2+vl**2+wl**2)
        hl=(rhoel+pl)/rhol
        uconl=kx*ul+ky*vl+kz*wl

        !!Extract primative variables for right state from qr
        !rhor=qr(1)
        !ur=qr(2)
        !vr=qr(3)
        !wr=qr(4)
        !pr=qr(5)
        !!Calculate rho*e, enthalpy, and contravariant velocity at right state
        rhoer=(pr/GAM1)+0.50d0*rhor*(ur**2+vr**2+wl**2)
        hr=(rhoer+pr)/rhor
        uconr=kx*ur+ky*vr+kz*wr

        !Calculate the flux based on the left state
        fql(1)=rhol*uconl
        fql(2)=rhol*ul*uconl+pl*kx
        fql(3)=rhol*vl*uconl+pl*ky
        fql(4)=rhol*wl*uconl+pl*kz
        fql(5)=(rhoel+pl)*uconl

        !Calculate the flux based on the right state
        fqr(1)=rhor*uconr
        fqr(2)=rhor*ur*uconr+pr*kx
        fqr(3)=rhor*vr*uconr+pr*ky
        fqr(4)=rhor*wr*uconr+pr*kz
        fqr(5)=(rhoer+pr)*uconr

        !Calculate Roe averaged variables
        sc1=SQRT(rhol); sc2=SQRT(rhor); sc3=sc1+sc2
        rhoh=sc1*sc2
        uh=(ul*sc1+ur*sc2)/sc3
        vh=(vl*sc1+vr*sc2)/sc3
        wh=(wl*sc1+wr*sc2)/sc3
        hh=(hl*sc1+hr*sc2)/sc3
        ch=SQRT(GAM1*(hh-0.50d0*(uh**2+vh**2+wh**2)))

        !Eigenvalue matrix
        lamh(1,1)=kx*uh+ky*vh+kz*wh
        lamh(2,2)=lamh(1,1)
        lamh(3,3)=lamh(1,1)
        lamh(4,4)=lamh(1,1)+ch*magk
        lamh(5,5)=lamh(1,1)-ch*magk
        lamh=abs(lamh)

        !Apply entropy fix according to formula by Harten
        do invars=1,nvars
          if (lamh(invars,invars) .lt. 2.0d0*epHarten*ch) then
            lamh(invars,invars)=(lamh(invars,invars)**2)/(4.0d0*epHarten*ch)+epHarten*ch
          endif
        enddo

        !Convenient variables for assembling matrices
        alfa=rhoh/(ch*SQRT(2.0d0))
        beta=1.0d0/(rhoh*ch*SQRT(2.0d0))
        theta=(kx*uh+ky*vh+kz*wh)/magk
        phi2=0.50d0*GAM1*(uh**2+vh**2+wh**2)

        !Assemble [T_k] matrix (See writeup for details)
        tk(1,1)=kxn
        tk(1,2)=kyn
        tk(1,3)=kzn
        tk(1,4)=alfa
        tk(1,5)=alfa

        tk(2,1)=kxn*uh
        tk(2,2)=kyn*uh-kzn*rhoh
        tk(2,3)=kzn*uh+kyn*rhoh
        tk(2,4)=alfa*(uh+kxn*ch)
        tk(2,5)=alfa*(uh-kxn*ch)

        tk(3,1)=kxn*vh+kzn*rhoh
        tk(3,2)=kyn*vh
        tk(3,3)=kzn*vh-kxn*rhoh
        tk(3,4)=alfa*(vh+kyn*ch)
        tk(3,5)=alfa*(vh-kyn*ch)

        tk(4,1)=kxn*wh-kyn*rhoh
        tk(4,2)=kyn*wh+kxn*rhoh
        tk(4,3)=kzn*wh
        tk(4,4)=alfa*(wh+kzn*ch)
        tk(4,5)=alfa*(wh-kzn*ch)

        tk(5,1)=kxn*phi2/GAM1+rhoh*(kzn*vh-kyn*wh)
        tk(5,2)=kyn*phi2/GAM1+rhoh*(kxn*wh-kzn*uh)
        tk(5,3)=kzn*phi2/GAM1+rhoh*(kyn*uh-kxn*vh)
        tk(5,4)=alfa*((phi2+ch**2)/GAM1+theta*ch)
        tk(5,5)=alfa*((phi2+ch**2)/GAM1-theta*ch)

        !Assemble [T_k]^-1 matrix (See writeup for details)
        tkInv(1,1)=kxn*(1.0d0-phi2/ch**2)-(kzn*vh-kyn*wh)/rhoh
        tkInv(1,2)=kxn*GAM1*uh/ch**2
        tkInv(1,3)=kxn*GAM1*vh/ch**2+kzn/rhoh
        tkInv(1,4)=kxn*GAM1*wh/ch**2-kyn/rhoh
        tkInv(1,5)=-kxn*GAM1/ch**2

        tkInv(2,1)=kyn*(1.0d0-phi2/ch**2)-(kxn*wh-kzn*uh)/rhoh
        tkInv(2,2)=kyn*GAM1*uh/ch**2-kzn/rhoh
        tkInv(2,3)=kyn*GAM1*vh/ch**2
        tkInv(2,4)=kyn*GAM1*wh/ch**2+kxn/rhoh
        tkInv(2,5)=-kyn*GAM1/ch**2

        tkInv(3,1)=kzn*(1.0d0-phi2/ch**2)-(kyn*uh-kxn*vh)/rhoh
        tkInv(3,2)=kzn*GAM1*uh/ch**2+kyn/rhoh
        tkInv(3,3)=kzn*GAM1*vh/ch**2-kxn/rhoh
        tkInv(3,4)=kzn*GAM1*wh/ch**2
        tkInv(3,5)=-kzn*GAM1/ch**2

        tkInv(4,1)=beta*(phi2-theta*ch)
        tkInv(4,2)=-beta*(GAM1*uh-kxn*ch)
        tkInv(4,3)=-beta*(GAM1*vh-kyn*ch)
        tkInv(4,4)=-beta*(GAM1*wh-kzn*ch)
        tkInv(4,5)=beta*GAM1

        tkInv(5,1)=beta*(phi2+theta*ch)
        tkInv(5,2)=-beta*(GAM1*uh+kxn*ch)
        tkInv(5,3)=-beta*(GAM1*vh+kyn*ch)
        tkInv(5,4)=-beta*(GAM1*wh+kzn*ch)
        tkInv(5,5)=beta*GAM1

        !Calculate left state in conservative variables
        ql(1)=rhol
        ql(2)=rhol*ul
        ql(3)=rhol*vl
        ql(4)=rhol*wl
        ql(5)=rhoel

        !Calculate right state in conservative variables
        qr(1)=rhor
        qr(2)=rhor*ur
        qr(3)=rhor*vr
        qr(4)=rhor*wr
        qr(5)=rhoer

        !Calculate difference between left and right states 
        dq=qr-ql

        !Calculate flux at interface: F=0.5*(fql+fqr-[tk][lamh][tkInv](dq))
        tm1(1)=tkInv(1,1)*dq(1)+tkInv(1,2)*dq(2)+tkInv(1,3)*dq(3)+tkInv(1,4)*dq(4)+tkInv(1,5)*dq(5)
        tm1(2)=tkInv(2,1)*dq(1)+tkInv(2,2)*dq(2)+tkInv(2,3)*dq(3)+tkInv(2,4)*dq(4)+tkInv(2,5)*dq(5)
        tm1(3)=tkInv(3,1)*dq(1)+tkInv(3,2)*dq(2)+tkInv(3,3)*dq(3)+tkInv(3,4)*dq(4)+tkInv(3,5)*dq(5)
        tm1(4)=tkInv(4,1)*dq(1)+tkInv(4,2)*dq(2)+tkInv(4,3)*dq(3)+tkInv(4,4)*dq(4)+tkInv(4,5)*dq(5)
        tm1(5)=tkInv(5,1)*dq(1)+tkInv(5,2)*dq(2)+tkInv(5,3)*dq(3)+tkInv(5,4)*dq(4)+tkInv(5,5)*dq(5)

        tm2(1)=lamh(1,1)*tm1(1)
        tm2(2)=lamh(2,2)*tm1(2)
        tm2(3)=lamh(3,3)*tm1(3)
        tm2(4)=lamh(4,4)*tm1(4)
        tm2(5)=lamh(5,5)*tm1(5)

        tm3(1)=tk(1,1)*tm2(1)+tk(1,2)*tm2(2)+tk(1,3)*tm2(3)+tk(1,4)*tm2(4)+tk(1,5)*tm2(5)
        tm3(2)=tk(2,1)*tm2(1)+tk(2,2)*tm2(2)+tk(2,3)*tm2(3)+tk(2,4)*tm2(4)+tk(2,5)*tm2(5)
        tm3(3)=tk(3,1)*tm2(1)+tk(3,2)*tm2(2)+tk(3,3)*tm2(3)+tk(3,4)*tm2(4)+tk(3,5)*tm2(5)
        tm3(4)=tk(4,1)*tm2(1)+tk(4,2)*tm2(2)+tk(4,3)*tm2(3)+tk(4,4)*tm2(4)+tk(4,5)*tm2(5)
        tm3(5)=tk(5,1)*tm2(1)+tk(5,2)*tm2(2)+tk(5,3)*tm2(3)+tk(5,4)*tm2(4)+tk(5,5)*tm2(5)

        !roeFlux=0.50d0*(fql+fqr-tm3)        
        !roeFlux=0.50d0*(fql+fqr-tm3)*xjacm*mag        
        !roeFlux=0.50d0*(fql+fqr-tm3)*xjacm        
        do ivar=1,nvars
          xiflux(i,j,k,ivar) = (hf*(fql(ivar) + fqr(ivar) -tm3(ivar)))*xjacm
        enddo
      enddo
    enddo
  enddo
end subroutine xiroe

subroutine etroe
  use mod_params
  use nvtx
  implicit none
  real (kind=8) :: &
    &kx,ky,kz,magk,xjacm,xjaci,xjacip,kxn,kyn,kzn,&
    &al,rhol,ul,vl,wl,pl,rhoel,hl,uconl,&
    &ar,rhor,ur,vr,wr,pr,rhoer,hr,uconr,&
    &fql(nvars),fqr(nvars),&
    &sc1,sc2,sc3,&
    &rhoh,uh,vh,wh,hh,ch,&
    &alfa,beta,theta,phi2,& 
    &lamh(nvars,nvars),tk(nvars,nvars),tkInv(nvars,nvars),&
    &ql(nvars),qr(nvars),&
    &dq(nvars),tm1(nvars),tm2(nvars),tm3(nvars)
  real (kind=8) :: epHarten = 0.0010d0
  integer:: i,j,k,ivar,invars
  !This subroutine implements the Roe scheme to determine the flux at an
  !interface given the left and right states. This subroutine is generalized so
  !that it can calculate the flux in any direction given the proper grid
  !metrics as inputs.

  call nvtxStartRange("recon(2)")
  call recon(2)
  call nvtxEndRange
  !$acc parallel loop collapse(3) &
  !$acc& present(etax,etay,etaz,xjac) &
  !$acc& private(fql,fqr,ql,qr,dq,tm1,tm2,tm3,lamh,tk,tkInv) &
  !$acc& private(magk,sc3,alfa,beta,ch,vr,ur,uconr,rhoer,sc2,sc1,rhor,uconl,hl,rhoel,theta,hh,pr,pl,wl,vl,ul,rhol,phi2,wh,vh,uh,rhoh,kz,ky,kx,kzn,kyn,kxn,hr,xjaci,wr,xjacm,xjacip)
  do k=1,nz-1
    do j=1,ny-1
      do i=1,nx-1
        xjacm = hf*(xjac(i,j,k)+xjac(i,j+1,k))
        xjaci = xjac(i,j,k)/xjacm
        xjacip = xjac(i,j+1,k)/xjacm
        kx = hf*(etax(i,j,k)*xjaci+etax(i,j+1,k)*xjacip) 
        ky = hf*(etay(i,j,k)*xjaci+etay(i,j+1,k)*xjacip) 
        kz = hf*(etaz(i,j,k)*xjaci+etaz(i,j+1,k)*xjacip) 

        magk = sqrt(kx**2 + ky**2 +kz**2) 
        kxn = kx/magk
        kyn = ky/magk
        kzn = kz/magk

        !!  Left state
        rhol = primL(i,j,k,1) 
        ul = primL(i,j,k,2) 
        vl = primL(i,j,k,3) 
        wl = primL(i,j,k,4) 
        pl = primL(i,j,k,5) 

        !!  Right state
        rhor = primR(i,j,k,1) 
        ur = primR(i,j,k,2) 
        vr = primR(i,j,k,3) 
        wr = primR(i,j,k,4) 
        pr = primR(i,j,k,5) 

        !!Extract primative variables for left state from ql
        !rhol=ql(1)
        !ul=ql(2)
        !vl=ql(3)
        !wl=ql(4)
        !pl=ql(5)
        !!Calculate rho*e, enthalpy, and contravariant velocity at left state
        rhoel=(pl/GAM1)+0.50d0*rhol*(ul**2+vl**2+wl**2)
        hl=(rhoel+pl)/rhol
        uconl=kx*ul+ky*vl+kz*wl

        !!Extract primative variables for right state from qr
        !rhor=qr(1)
        !ur=qr(2)
        !vr=qr(3)
        !wr=qr(4)
        !pr=qr(5)
        !!Calculate rho*e, enthalpy, and contravariant velocity at right state
        rhoer=(pr/GAM1)+0.50d0*rhor*(ur**2+vr**2+wl**2)
        hr=(rhoer+pr)/rhor
        uconr=kx*ur+ky*vr+kz*wr

        !Calculate the flux based on the left state
        fql(1)=rhol*uconl
        fql(2)=rhol*ul*uconl+pl*kx
        fql(3)=rhol*vl*uconl+pl*ky
        fql(4)=rhol*wl*uconl+pl*kz
        fql(5)=(rhoel+pl)*uconl

        !Calculate the flux based on the right state
        fqr(1)=rhor*uconr
        fqr(2)=rhor*ur*uconr+pr*kx
        fqr(3)=rhor*vr*uconr+pr*ky
        fqr(4)=rhor*wr*uconr+pr*kz
        fqr(5)=(rhoer+pr)*uconr

        !Calculate Roe averaged variables
        sc1=SQRT(rhol); sc2=SQRT(rhor); sc3=sc1+sc2
        rhoh=sc1*sc2
        uh=(ul*sc1+ur*sc2)/sc3
        vh=(vl*sc1+vr*sc2)/sc3
        wh=(wl*sc1+wr*sc2)/sc3
        hh=(hl*sc1+hr*sc2)/sc3
        ch=SQRT(GAM1*(hh-0.50d0*(uh**2+vh**2+wh**2)))

        !Eigenvalue matrix
        lamh(1,1)=kx*uh+ky*vh+kz*wh
        lamh(2,2)=lamh(1,1)
        lamh(3,3)=lamh(1,1)
        lamh(4,4)=lamh(1,1)+ch*magk
        lamh(5,5)=lamh(1,1)-ch*magk
        lamh=abs(lamh)

        !Apply entropy fix according to formula by Harten
        do invars=1,nvars
          if (lamh(invars,invars) .lt. 2.0d0*epHarten*ch) then
            lamh(invars,invars)=(lamh(invars,invars)**2)/(4.0d0*epHarten*ch)+epHarten*ch
          endif
        enddo

        !Convenient variables for assembling matrices
        alfa=rhoh/(ch*SQRT(2.0d0))
        beta=1.0d0/(rhoh*ch*SQRT(2.0d0))
        theta=(kx*uh+ky*vh+kz*wh)/magk
        phi2=0.50d0*GAM1*(uh**2+vh**2+wh**2)

        !Assemble [T_k] matrix (See writeup for details)
        tk(1,1)=kxn
        tk(1,2)=kyn
        tk(1,3)=kzn
        tk(1,4)=alfa
        tk(1,5)=alfa

        tk(2,1)=kxn*uh
        tk(2,2)=kyn*uh-kzn*rhoh
        tk(2,3)=kzn*uh+kyn*rhoh
        tk(2,4)=alfa*(uh+kxn*ch)
        tk(2,5)=alfa*(uh-kxn*ch)

        tk(3,1)=kxn*vh+kzn*rhoh
        tk(3,2)=kyn*vh
        tk(3,3)=kzn*vh-kxn*rhoh
        tk(3,4)=alfa*(vh+kyn*ch)
        tk(3,5)=alfa*(vh-kyn*ch)

        tk(4,1)=kxn*wh-kyn*rhoh
        tk(4,2)=kyn*wh+kxn*rhoh
        tk(4,3)=kzn*wh
        tk(4,4)=alfa*(wh+kzn*ch)
        tk(4,5)=alfa*(wh-kzn*ch)

        tk(5,1)=kxn*phi2/GAM1+rhoh*(kzn*vh-kyn*wh)
        tk(5,2)=kyn*phi2/GAM1+rhoh*(kxn*wh-kzn*uh)
        tk(5,3)=kzn*phi2/GAM1+rhoh*(kyn*uh-kxn*vh)
        tk(5,4)=alfa*((phi2+ch**2)/GAM1+theta*ch)
        tk(5,5)=alfa*((phi2+ch**2)/GAM1-theta*ch)

        !Assemble [T_k]^-1 matrix (See writeup for details)
        tkInv(1,1)=kxn*(1.0d0-phi2/ch**2)-(kzn*vh-kyn*wh)/rhoh
        tkInv(1,2)=kxn*GAM1*uh/ch**2
        tkInv(1,3)=kxn*GAM1*vh/ch**2+kzn/rhoh
        tkInv(1,4)=kxn*GAM1*wh/ch**2-kyn/rhoh
        tkInv(1,5)=-kxn*GAM1/ch**2

        tkInv(2,1)=kyn*(1.0d0-phi2/ch**2)-(kxn*wh-kzn*uh)/rhoh
        tkInv(2,2)=kyn*GAM1*uh/ch**2-kzn/rhoh
        tkInv(2,3)=kyn*GAM1*vh/ch**2
        tkInv(2,4)=kyn*GAM1*wh/ch**2+kxn/rhoh
        tkInv(2,5)=-kyn*GAM1/ch**2

        tkInv(3,1)=kzn*(1.0d0-phi2/ch**2)-(kyn*uh-kxn*vh)/rhoh
        tkInv(3,2)=kzn*GAM1*uh/ch**2+kyn/rhoh
        tkInv(3,3)=kzn*GAM1*vh/ch**2-kxn/rhoh
        tkInv(3,4)=kzn*GAM1*wh/ch**2
        tkInv(3,5)=-kzn*GAM1/ch**2

        tkInv(4,1)=beta*(phi2-theta*ch)
        tkInv(4,2)=-beta*(GAM1*uh-kxn*ch)
        tkInv(4,3)=-beta*(GAM1*vh-kyn*ch)
        tkInv(4,4)=-beta*(GAM1*wh-kzn*ch)
        tkInv(4,5)=beta*GAM1

        tkInv(5,1)=beta*(phi2+theta*ch)
        tkInv(5,2)=-beta*(GAM1*uh+kxn*ch)
        tkInv(5,3)=-beta*(GAM1*vh+kyn*ch)
        tkInv(5,4)=-beta*(GAM1*wh+kzn*ch)
        tkInv(5,5)=beta*GAM1

        !Calculate left state in conservative variables
        ql(1)=rhol
        ql(2)=rhol*ul
        ql(3)=rhol*vl
        ql(4)=rhol*wl
        ql(5)=rhoel

        !Calculate right state in conservative variables
        qr(1)=rhor
        qr(2)=rhor*ur
        qr(3)=rhor*vr
        qr(4)=rhor*wr
        qr(5)=rhoer

        !Calculate difference between left and right states 
        dq=qr-ql

        !Calculate flux at interface: F=0.5*(fql+fqr-[tk][lamh][tkInv](dq))
        tm1(1)=tkInv(1,1)*dq(1)+tkInv(1,2)*dq(2)+tkInv(1,3)*dq(3)+tkInv(1,4)*dq(4)+tkInv(1,5)*dq(5)
        tm1(2)=tkInv(2,1)*dq(1)+tkInv(2,2)*dq(2)+tkInv(2,3)*dq(3)+tkInv(2,4)*dq(4)+tkInv(2,5)*dq(5)
        tm1(3)=tkInv(3,1)*dq(1)+tkInv(3,2)*dq(2)+tkInv(3,3)*dq(3)+tkInv(3,4)*dq(4)+tkInv(3,5)*dq(5)
        tm1(4)=tkInv(4,1)*dq(1)+tkInv(4,2)*dq(2)+tkInv(4,3)*dq(3)+tkInv(4,4)*dq(4)+tkInv(4,5)*dq(5)
        tm1(5)=tkInv(5,1)*dq(1)+tkInv(5,2)*dq(2)+tkInv(5,3)*dq(3)+tkInv(5,4)*dq(4)+tkInv(5,5)*dq(5)

        tm2(1)=lamh(1,1)*tm1(1)
        tm2(2)=lamh(2,2)*tm1(2)
        tm2(3)=lamh(3,3)*tm1(3)
        tm2(4)=lamh(4,4)*tm1(4)
        tm2(5)=lamh(5,5)*tm1(5)

        tm3(1)=tk(1,1)*tm2(1)+tk(1,2)*tm2(2)+tk(1,3)*tm2(3)+tk(1,4)*tm2(4)+tk(1,5)*tm2(5)
        tm3(2)=tk(2,1)*tm2(1)+tk(2,2)*tm2(2)+tk(2,3)*tm2(3)+tk(2,4)*tm2(4)+tk(2,5)*tm2(5)
        tm3(3)=tk(3,1)*tm2(1)+tk(3,2)*tm2(2)+tk(3,3)*tm2(3)+tk(3,4)*tm2(4)+tk(3,5)*tm2(5)
        tm3(4)=tk(4,1)*tm2(1)+tk(4,2)*tm2(2)+tk(4,3)*tm2(3)+tk(4,4)*tm2(4)+tk(4,5)*tm2(5)
        tm3(5)=tk(5,1)*tm2(1)+tk(5,2)*tm2(2)+tk(5,3)*tm2(3)+tk(5,4)*tm2(4)+tk(5,5)*tm2(5)

        !roeFlux=0.50d0*(fql+fqr-tm3)        
        !roeFlux=0.50d0*(fql+fqr-tm3)*xjacm*mag        
        !roeFlux=0.50d0*(fql+fqr-tm3)*xjacm        
        do ivar=1,nvars
          etflux(i,j,k,ivar) = (hf*(fql(ivar) + fqr(ivar) -tm3(ivar)))*xjacm
        enddo
      enddo
    enddo
  enddo
end subroutine etroe

subroutine ztroe
  use mod_params
  use nvtx
  implicit none
  real (kind=8) :: &
    &kx,ky,kz,magk,xjacm,xjaci,xjacip,kxn,kyn,kzn,&
    &al,rhol,ul,vl,wl,pl,rhoel,hl,uconl,&
    &ar,rhor,ur,vr,wr,pr,rhoer,hr,uconr,&
    &fql(nvars),fqr(nvars),&
    &sc1,sc2,sc3,&
    &rhoh,uh,vh,wh,hh,ch,&
    &alfa,beta,theta,phi2,& 
    &lamh(nvars,nvars),tk(nvars,nvars),tkInv(nvars,nvars),&
    &ql(nvars),qr(nvars),&
    &dq(nvars),tm1(nvars),tm2(nvars),tm3(nvars)
  real (kind=8) :: epHarten = 0.0010d0
  integer:: i,j,k,ivar,invars
  !This subroutine implements the Roe scheme to determine the flux at an
  !interface given the left and right states. This subroutine is generalized so
  !that it can calculate the flux in any direction given the proper grid
  !metrics as inputs.

  call nvtxStartRange("recon(3)")
  call recon(3)
  call nvtxEndRange
  !$acc parallel loop collapse(3) &
  !$acc& present(ztax,ztay,ztaz,xjac) &
  !$acc& private(fql,fqr,ql,qr,dq,tm1,tm2,tm3,lamh,tk,tkInv) &
  !$acc& private(magk,sc3,alfa,beta,ch,vr,ur,uconr,rhoer,sc2,sc1,rhor,uconl,hl,rhoel,theta,hh,pr,pl,wl,vl,ul,rhol,phi2,wh,vh,uh,rhoh,kz,ky,kx,kzn,kyn,kxn,hr,xjaci,wr,xjacm,xjacip)
  do k=1,nz-1
    do j=1,ny-1
      do i=1,nx-1
        xjacm = hf*(xjac(i,j,k)+xjac(i,j,k+1))
        xjaci = xjac(i,j,k)/xjacm
        xjacip = xjac(i,j,k+1)/xjacm
        kx = hf*(ztax(i,j,k)*xjaci+ztax(i,j,k+1)*xjacip) 
        ky = hf*(ztay(i,j,k)*xjaci+ztay(i,j,k+1)*xjacip) 
        kz = hf*(ztaz(i,j,k)*xjaci+ztaz(i,j,k+1)*xjacip) 

        magk = sqrt(kx**2 + ky**2 +kz**2) 
        kxn = kx/magk
        kyn = ky/magk
        kzn = kz/magk

        !!  Left state
        rhol = primL(i,j,k,1) 
        ul = primL(i,j,k,2) 
        vl = primL(i,j,k,3) 
        wl = primL(i,j,k,4) 
        pl = primL(i,j,k,5) 

        !!  Right state
        rhor = primR(i,j,k,1) 
        ur = primR(i,j,k,2) 
        vr = primR(i,j,k,3) 
        wr = primR(i,j,k,4) 
        pr = primR(i,j,k,5) 

        !!Extract primative variables for left state from ql
        !rhol=ql(1)
        !ul=ql(2)
        !vl=ql(3)
        !wl=ql(4)
        !pl=ql(5)
        !!Calculate rho*e, enthalpy, and contravariant velocity at left state
        rhoel=(pl/GAM1)+0.50d0*rhol*(ul**2+vl**2+wl**2)
        hl=(rhoel+pl)/rhol
        uconl=kx*ul+ky*vl+kz*wl

        !!Extract primative variables for right state from qr
        !rhor=qr(1)
        !ur=qr(2)
        !vr=qr(3)
        !wr=qr(4)
        !pr=qr(5)
        !!Calculate rho*e, enthalpy, and contravariant velocity at right state
        rhoer=(pr/GAM1)+0.50d0*rhor*(ur**2+vr**2+wl**2)
        hr=(rhoer+pr)/rhor
        uconr=kx*ur+ky*vr+kz*wr

        !Calculate the flux based on the left state
        fql(1)=rhol*uconl
        fql(2)=rhol*ul*uconl+pl*kx
        fql(3)=rhol*vl*uconl+pl*ky
        fql(4)=rhol*wl*uconl+pl*kz
        fql(5)=(rhoel+pl)*uconl

        !Calculate the flux based on the right state
        fqr(1)=rhor*uconr
        fqr(2)=rhor*ur*uconr+pr*kx
        fqr(3)=rhor*vr*uconr+pr*ky
        fqr(4)=rhor*wr*uconr+pr*kz
        fqr(5)=(rhoer+pr)*uconr

        !Calculate Roe averaged variables
        sc1=SQRT(rhol); sc2=SQRT(rhor); sc3=sc1+sc2
        rhoh=sc1*sc2
        uh=(ul*sc1+ur*sc2)/sc3
        vh=(vl*sc1+vr*sc2)/sc3
        wh=(wl*sc1+wr*sc2)/sc3
        hh=(hl*sc1+hr*sc2)/sc3
        ch=SQRT(GAM1*(hh-0.50d0*(uh**2+vh**2+wh**2)))

        !Eigenvalue matrix
        lamh(1,1)=kx*uh+ky*vh+kz*wh
        lamh(2,2)=lamh(1,1)
        lamh(3,3)=lamh(1,1)
        lamh(4,4)=lamh(1,1)+ch*magk
        lamh(5,5)=lamh(1,1)-ch*magk
        lamh=abs(lamh)

        !Apply entropy fix according to formula by Harten
        do invars=1,nvars
          if (lamh(invars,invars) .lt. 2.0d0*epHarten*ch) then
            lamh(invars,invars)=(lamh(invars,invars)**2)/(4.0d0*epHarten*ch)+epHarten*ch
          endif
        enddo

        !Convenient variables for assembling matrices
        alfa=rhoh/(ch*SQRT(2.0d0))
        beta=1.0d0/(rhoh*ch*SQRT(2.0d0))
        theta=(kx*uh+ky*vh+kz*wh)/magk
        phi2=0.50d0*GAM1*(uh**2+vh**2+wh**2)

        !Assemble [T_k] matrix (See writeup for details)
        tk(1,1)=kxn
        tk(1,2)=kyn
        tk(1,3)=kzn
        tk(1,4)=alfa
        tk(1,5)=alfa

        tk(2,1)=kxn*uh
        tk(2,2)=kyn*uh-kzn*rhoh
        tk(2,3)=kzn*uh+kyn*rhoh
        tk(2,4)=alfa*(uh+kxn*ch)
        tk(2,5)=alfa*(uh-kxn*ch)

        tk(3,1)=kxn*vh+kzn*rhoh
        tk(3,2)=kyn*vh
        tk(3,3)=kzn*vh-kxn*rhoh
        tk(3,4)=alfa*(vh+kyn*ch)
        tk(3,5)=alfa*(vh-kyn*ch)

        tk(4,1)=kxn*wh-kyn*rhoh
        tk(4,2)=kyn*wh+kxn*rhoh
        tk(4,3)=kzn*wh
        tk(4,4)=alfa*(wh+kzn*ch)
        tk(4,5)=alfa*(wh-kzn*ch)

        tk(5,1)=kxn*phi2/GAM1+rhoh*(kzn*vh-kyn*wh)
        tk(5,2)=kyn*phi2/GAM1+rhoh*(kxn*wh-kzn*uh)
        tk(5,3)=kzn*phi2/GAM1+rhoh*(kyn*uh-kxn*vh)
        tk(5,4)=alfa*((phi2+ch**2)/GAM1+theta*ch)
        tk(5,5)=alfa*((phi2+ch**2)/GAM1-theta*ch)

        !Assemble [T_k]^-1 matrix (See writeup for details)
        tkInv(1,1)=kxn*(1.0d0-phi2/ch**2)-(kzn*vh-kyn*wh)/rhoh
        tkInv(1,2)=kxn*GAM1*uh/ch**2
        tkInv(1,3)=kxn*GAM1*vh/ch**2+kzn/rhoh
        tkInv(1,4)=kxn*GAM1*wh/ch**2-kyn/rhoh
        tkInv(1,5)=-kxn*GAM1/ch**2

        tkInv(2,1)=kyn*(1.0d0-phi2/ch**2)-(kxn*wh-kzn*uh)/rhoh
        tkInv(2,2)=kyn*GAM1*uh/ch**2-kzn/rhoh
        tkInv(2,3)=kyn*GAM1*vh/ch**2
        tkInv(2,4)=kyn*GAM1*wh/ch**2+kxn/rhoh
        tkInv(2,5)=-kyn*GAM1/ch**2

        tkInv(3,1)=kzn*(1.0d0-phi2/ch**2)-(kyn*uh-kxn*vh)/rhoh
        tkInv(3,2)=kzn*GAM1*uh/ch**2+kyn/rhoh
        tkInv(3,3)=kzn*GAM1*vh/ch**2-kxn/rhoh
        tkInv(3,4)=kzn*GAM1*wh/ch**2
        tkInv(3,5)=-kzn*GAM1/ch**2

        tkInv(4,1)=beta*(phi2-theta*ch)
        tkInv(4,2)=-beta*(GAM1*uh-kxn*ch)
        tkInv(4,3)=-beta*(GAM1*vh-kyn*ch)
        tkInv(4,4)=-beta*(GAM1*wh-kzn*ch)
        tkInv(4,5)=beta*GAM1

        tkInv(5,1)=beta*(phi2+theta*ch)
        tkInv(5,2)=-beta*(GAM1*uh+kxn*ch)
        tkInv(5,3)=-beta*(GAM1*vh+kyn*ch)
        tkInv(5,4)=-beta*(GAM1*wh+kzn*ch)
        tkInv(5,5)=beta*GAM1

        !Calculate left state in conservative variables
        ql(1)=rhol
        ql(2)=rhol*ul
        ql(3)=rhol*vl
        ql(4)=rhol*wl
        ql(5)=rhoel

        !Calculate right state in conservative variables
        qr(1)=rhor
        qr(2)=rhor*ur
        qr(3)=rhor*vr
        qr(4)=rhor*wr
        qr(5)=rhoer

        !Calculate difference between left and right states 
        dq=qr-ql

        !Calculate flux at interface: F=0.5*(fql+fqr-[tk][lamh][tkInv](dq))
        tm1(1)=tkInv(1,1)*dq(1)+tkInv(1,2)*dq(2)+tkInv(1,3)*dq(3)+tkInv(1,4)*dq(4)+tkInv(1,5)*dq(5)
        tm1(2)=tkInv(2,1)*dq(1)+tkInv(2,2)*dq(2)+tkInv(2,3)*dq(3)+tkInv(2,4)*dq(4)+tkInv(2,5)*dq(5)
        tm1(3)=tkInv(3,1)*dq(1)+tkInv(3,2)*dq(2)+tkInv(3,3)*dq(3)+tkInv(3,4)*dq(4)+tkInv(3,5)*dq(5)
        tm1(4)=tkInv(4,1)*dq(1)+tkInv(4,2)*dq(2)+tkInv(4,3)*dq(3)+tkInv(4,4)*dq(4)+tkInv(4,5)*dq(5)
        tm1(5)=tkInv(5,1)*dq(1)+tkInv(5,2)*dq(2)+tkInv(5,3)*dq(3)+tkInv(5,4)*dq(4)+tkInv(5,5)*dq(5)

        tm2(1)=lamh(1,1)*tm1(1)
        tm2(2)=lamh(2,2)*tm1(2)
        tm2(3)=lamh(3,3)*tm1(3)
        tm2(4)=lamh(4,4)*tm1(4)
        tm2(5)=lamh(5,5)*tm1(5)

        tm3(1)=tk(1,1)*tm2(1)+tk(1,2)*tm2(2)+tk(1,3)*tm2(3)+tk(1,4)*tm2(4)+tk(1,5)*tm2(5)
        tm3(2)=tk(2,1)*tm2(1)+tk(2,2)*tm2(2)+tk(2,3)*tm2(3)+tk(2,4)*tm2(4)+tk(2,5)*tm2(5)
        tm3(3)=tk(3,1)*tm2(1)+tk(3,2)*tm2(2)+tk(3,3)*tm2(3)+tk(3,4)*tm2(4)+tk(3,5)*tm2(5)
        tm3(4)=tk(4,1)*tm2(1)+tk(4,2)*tm2(2)+tk(4,3)*tm2(3)+tk(4,4)*tm2(4)+tk(4,5)*tm2(5)
        tm3(5)=tk(5,1)*tm2(1)+tk(5,2)*tm2(2)+tk(5,3)*tm2(3)+tk(5,4)*tm2(4)+tk(5,5)*tm2(5)

        !roeFlux=0.50d0*(fql+fqr-tm3)        
        !roeFlux=0.50d0*(fql+fqr-tm3)*xjacm*mag        
        !roeFlux=0.50d0*(fql+fqr-tm3)*xjacm        
        do ivar=1,nvars
          ztflux(i,j,k,ivar) = (hf*(fql(ivar) + fqr(ivar) -tm3(ivar)))*xjacm
        enddo
      enddo
    enddo
  enddo
end subroutine ztroe
