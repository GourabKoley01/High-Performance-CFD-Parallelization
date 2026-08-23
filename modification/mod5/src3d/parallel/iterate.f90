subroutine iterate(itr)
  use mod_params
  use mod_spmd
  use nvtx
  implicit none
  integer, intent(in):: itr
  integer:: i,j,k
  integer:: ivar,m

  call nvtxStartRange("Cons_init")
  !$acc kernels present(q0, q)
  q0 = q
  !$acc end kernels
  call nvtxEndRange

  select case (itime_scheme)

  case(3)
    do m = 1,3
      
      call nvtxStartRange("Comms prim")
      ! NOTE: the !$acc update self/device calls that used to live here (with
      ! generic 1:noverlap+1 / nx-noverlap:nx ranges) have been moved INSIDE
      ! CommunicateAll, using exact isgrid/iegrid-based ranges that match what
      ! is actually packed/unpacked, plus an explicit "!$acc wait" before the
      ! host copy. That generic-range mismatch was the likely source of the
      ! NaN/Inf blowup. Do not re-add update self/device calls here.
      call CommunicateAll(prim, nvars)
      call nvtxEndRange
      
      call nvtxStartRange("Apply_BC")
      call ApplyBC
      call nvtxEndRange
      
      !Above Apply BC call for MultiblockBC --  
      
      call nvtxStartRange("GetFlux")
      call getflux
      call nvtxEndRange
      
      call nvtxStartRange("RK3 update cons")
      !$acc parallel loop collapse(4) present(rhs, q, q0)
      do ivar = 1,nvars
        do k = ksgrid,kegrid
          do j = jsgrid,jegrid
            do i = isgrid,iegrid
              q(i,j,k,ivar) = q0(i,j,k,ivar)*RK3A(m) + (q(i,j,k,ivar) - dt * rhs(i,j,k,ivar))*RK3B(m)
            enddo
          enddo
        enddo
      enddo
      call nvtxEndRange
      
      call nvtxStartRange("ConstoPrim")
      call ConsToPrim
      call nvtxEndRange
    enddo

 ! case(4)
 !   do m = 1,4
 !     call get_flux
 !     !call filter(itr)
 !     !print*,'7 itr,dt,max(q), min(q0), max(rhs)',itr, dt,maxval(q),&
 !!     !&maxval(q0),maxval(rhs)
 !     do l = 1,nvars
 !       do i = 2,nx-1
 !         do j = 2,ny-1
 !           do k = 2,nz-1
 !             q(i,j,k,l) = q0(i,j,k,l) - dt * rhs(i,j,k,l)*RK4(m)
 !           enddo
 !         enddo
 !       enddo
 !     enddo
 !     call cons2prim
 !   enddo

  end select

  if(iflux_scheme .eq. 10)then 
    if(ifilt_order .gt. 0)then
      call nvtxStartRange("Filter")
      call Filter(itr)
      call nvtxEndRange
      
      call nvtxStartRange("Comms Cons")
      !$acc update self(q(1:noverlap+1,:,:,:), q(nx-noverlap:nx,:,:,:))
      !$acc update self(q(:,1:noverlap+1,:,:), q(:,ny-noverlap:ny,:,:))
      !$acc update self(q(:,:,1:noverlap+1,:), q(:,:,nz-noverlap:nz,:))
      do ivar=1,nvars
        call Communicate(q(:,:,:,ivar),q(:,:,:,ivar))
      enddo
      !$acc update device(q(1:noverlap+1,:,:,:), q(nx-noverlap:nx,:,:,:))
      !$acc update device(q(:,1:noverlap+1,:,:), q(:,ny-noverlap:ny,:,:))
      !$acc update device(q(:,:,1:noverlap+1,:), q(:,:,nz-noverlap:nz,:))
      call nvtxEndRange
      
      call nvtxStartRange("ConstoPrim")
      call ConsToPrim
      call nvtxEndRange
    endif
  endif

  !call ApplyBC
  !do ivar=1,nvars
  !  call Communicate(q(:,:,:,ivar),q(:,:,:,ivar))
  !enddo
  !call ConsToPrim
end subroutine iterate
