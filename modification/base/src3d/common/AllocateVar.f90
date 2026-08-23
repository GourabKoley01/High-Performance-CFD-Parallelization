subroutine allocatevar
  use mod_params
  implicit none

  nmax = max(nx,ny,nz)
  allocate(prim(nx,ny,nz,nvars))
  if(imean .eq. 1)allocate(primmean(nx,ny,nz,nvars+6))
  allocate(T(nx,ny,nz))
  allocate(xmu(nx,ny,nz))
  allocate(rhs(nx,ny,nz,nvars))
  allocate(q(nx,ny,nz,nvars))
  allocate(q0(nx,ny,nz,nvars))
  allocate(xiflux(nx,ny,nz,nvars))
  allocate(etflux(nx,ny,nz,nvars))
  allocate(ztflux(nx,ny,nz,nvars))

  if(iflux_scheme .eq. 2)then
    allocate(primL(nx,ny,nz,nvars))
    allocate(primR(nx,ny,nz,nvars))
  endif
  !   allocate(primLj(nx,ny,nz,nvars))
  !   allocate(primRj(nx,ny,nz,nvars))
  !   allocate(xjacmi(nx,ny,nz))
  !   allocate(xjacmj(nx,ny,nz))

  ! allocate(xxi(nx,ny,nz))
  ! allocate(yxi(nx,ny,nz))
  ! allocate(xeta(nx,ny,nz))
  ! allocate(yeta(nx,ny,nz))

  ! allocate(xixmj(nx,ny,nz))
  ! allocate(etaxmj(nx,ny,nz))
  ! allocate(xiymj(nx,ny,nz))
  ! allocate(etaymj(nx,ny,nz))
  ! allocate(xixmi(nx,ny,nz))
  ! allocate(etaxmi(nx,ny,nz))
  ! allocate(xiymi(nx,ny,nz))
  ! allocate(etaymi(nx,ny,nz))

  allocate(xix(nx,ny,nz))
  allocate(etax(nx,ny,nz))
  allocate(ztax(nx,ny,nz))

  allocate(xiy(nx,ny,nz))
  allocate(etay(nx,ny,nz))
  allocate(ztay(nx,ny,nz))

  allocate(xiz(nx,ny,nz))
  allocate(etaz(nx,ny,nz))
  allocate(ztaz(nx,ny,nz))


  allocate(xjac(nx,ny,nz))

  xix = 0.
  xiy = 0.
  xiz = 0.

  etax = 0.
  etay = 0.
  etaz = 0.

  ztax = 0.
  ztay = 0.
  ztaz = 0.

  xjac = 0.
  !!TDG system
  allocate(atdg(nmax,nmax)) 
  allocate(btdg(nmax,nmax)) 
  allocate(ctdg(nmax,nmax)) 
  allocate(tdg(nmax,nmax,nvars)) 
  allocate(rhstdg(nmax,nmax,nvars)) 

  !!TDG_field system
  allocate(atdg_field(nx,ny,nz))
  allocate(btdg_field(nx,ny,nz))
  allocate(ctdg_field(nx,ny,nz))
  allocate(gamtdg_field(nx,ny,nz))
  allocate(rhstdg_field(nx,ny,nz,nvars))
  
  allocate(temp(nmax,nmax,nmax,nvars))

  !Inlet BC
  !allocate(uprof_in(ny))
  !allocate(uprof_out(ny))
end subroutine allocatevar
