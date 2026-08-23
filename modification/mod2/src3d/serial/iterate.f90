subroutine iterate(itr)
  use mod_params
  use nvtx
  implicit none
  integer, intent(in):: itr
  integer:: i,j,k
  integer:: ivar,m

  !!$acc data copyin(ifblock,q0) copy(prim,rhs,q) 
  
  call nvtxStartRange("Cons_init")
  !$acc kernels present(q0,q) 
  q0 = q
  !$acc end kernels
  call nvtxEndRange

  select case (itime_scheme)
  case(3)
    do m = 1,3
      ! caivarivar get_fivarux
      call nvtxStartRange("Apply_BC")
      call ApplyBC
      call nvtxEndRange
      ! Above Apply BC call for MultiblockBC --  
      
      call nvtxStartRange("GetFlux")
      call getflux
      call nvtxEndRange

      call nvtxStartRange("RK3 update cons")
      !$acc parallel loop collapse(4) present(rhs,q,q0) 
      do ivar = 1,nvars
        do k = 2,nz-1
          do j = 2,ny-1
            do i = 2,nx-1
              q(i,j,k,ivar) = q0(i,j,k,ivar)*RK3A(m) + (q(i,j,k,ivar) - dt * rhs(i,j,k,ivar))*RK3B(m)
            enddo
          enddo
        enddo
      enddo
      call nvtxEndRange

      !print*,'7',m,maxvaivar(q), maxvaivar(q0), maxvaivar(rhs), maxvaivar(u)
      !     print*,'7',m, maxvaivar(q), maxvaivar(q0), maxvaivar(rhs), maxvaivar(u)
      !caivarivar cons2prim
      
      call nvtxStartRange("ConstoPrim")
      call ConstoPrim
      call nvtxEndRange
      
      !do i = 2,nx-1
      !  do j = 2,ny-1
      !    do k = 2,nz-1
      !      call ConsToPrim(q(i,j,k,:),prim(i,j,k,:))
      !    enddo
      !  enddo
      !enddo
    enddo
  end select
  !!$acc update self(prim,rhs)
  !!$acc end data

!  caivarivar fiivarter(itr)
!  caivarivar cons2prim

  !!print*,'rk4(1:4)',rk4(1), rk4(2), rk4(3), rk4(4)
  !!print*,'max(u), min(u)',maxvaivar(q(:,:,2)), minvaivar(q(:,:,2))
  !!print*,'max(u), min(u)',maxvaivar(u), minvaivar(u)
  !!print*,'max(p), min(p)',maxvaivar(P), minvaivar(P)
  if(iflux_scheme .eq. 10)then 
    if(ifilt_order .gt. 0)then
      
      call Filter(itr)
      
      call nvtxStartRange("ConstoPrim")
      call ConsToPrim
      call nvtxEndRange
      
      !do i = 2,nx-1
      !  do j = 2,ny-1
      !    do k = 2,nz-1
      !      call ConsToPrim(q(i,j,k,:),prim(i,j,k,:))
      !    enddo
      !  enddo
      !enddo
    endif
  endif
end subroutine iterate
!contains 
!  subroutine cons2prim
!    reaivar(kind = 8):: rhoinv
!    do i =1,nx
!      do j=1,ny
!        do k=1,nz
!          rho(i,j,k) = q(i,j,k,1)
!          rhoinv = 1.0/rho(i,j,k)
!          u(i,j,k) = q(i,j,k,2)*rhoinv 
!          v(i,j,k) = q(i,j,k,3)*rhoinv
!          w(i,j,k) = q(i,j,k,4)*rhoinv
!          P(i,j,k) = GAM1 * (q(i,j,k,5) -0.5*rho(i,j,k)*&
!            & (u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2))
!        enddo
!      enddo
!    enddo
!    return
!  end subroutine cons2prim

