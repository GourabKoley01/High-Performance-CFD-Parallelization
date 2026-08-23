subroutine Initialize
  use mod_params
  use mod_spmd
  implicit none
  integer :: i,j,k
  integer :: iproc, ictr, nsize, ierr
  integer :: length
  character(len=60) :: filename
  character(len=20) :: path='./snapshots/'
  integer:: mpistat(mpi_status_size)
  integer :: nblk1,nx1, ny1,nz1, nvar1
  logical there
  real(kind=8) :: Pmin_global

  print*, itrst
  !GXM = 1.0/Rgas
  !REI = 1.0
  !FACTOR=1.0
  rhoinf = 1.0d0
  uinf = 0.0d0
  vinf = 0.0d0
  winf = 0.0d0
  pinf = 1.0d0/GXM
  !pinf = 1.0d0
  
  if(masterproc)then
    allocate(ug(nxg,nyg,nzg))
    allocate(vg(nxg,nyg,nzg))
    allocate(wg(nxg,nyg,nzg))
    allocate(Pg(nxg,nyg,nzg))
    allocate(rhog(nxg,nyg,nzg))

    !allocate(uprof_in(nyg))
    !  filename='Inlet_profile_i251xj173.txt'
    !  open(10, file=trim(filename))
    !  read(10,*) uprof_in
    !  close(10)
    ! Note boundary profile is taken from initial file here
    if (itrst .eq. 0)then 
        ug=uinf
        vg=vinf
        wg=winf
        Pg=pinf
        rhog=rhoinf
    else
      print*,'ITS A RESTART FROM ITR=',itrst
      filename='movsp'
      length=len_trim(filename) 
      write(filename(length+1:length+1),'(A1)') '_'
      length=len_trim(filename) 
      write(filename(length+1:length+7),'(I7.7)') itrst
      length=len_trim(filename) 
      write(filename(length+1:length+3),'(A3)') '.fn'
      open(10, file=trim(path)//trim(filename),form='UNFORMATTED',status='OLD',convert='BIG_ENDIAN')
      !open(10, file=trim(filename),form='UNFORMATTED',status='OLD',convert='BIG_ENDIAN')
      read(10) nx1,ny1,nz1,nvar1
      read(10) ug,vg,wg,Pg,rhog
      close(10)
      print*, 'Pmax =', maxval(Pg(:,:,:))
      !print*, 'Pmin =', minval(Pg(:,:,:))
      !print*, 'Umax =', maxval(ug(:,:,:))
      !print*, 'Umin =', minval(ug(:,:,:))
      !print*, 'Vmax =', maxval(vg(:,:,:))
      !print*, 'Vmin =', minval(vg(:,:,:))
      !print*, 'Wmax =', maxval(wg(:,:,:))
      !print*, 'Wmin =', minval(wg(:,:,:))
      !print*, 'Rhomax =', maxval(rhog(:,:,:))
      !print*, 'Rhomin =', minval(rhog(:,:,:))
      !call restart(itrst)
    endif
  endif

  call CommVars
  ! right after "call CommVars" in Initialize.f90, on EVERY rank (no masterproc guard)
  !print*, 'myrank=', myrank, 'local Pmax,Pmin=', maxval(prim(:,:,:,5)), minval(prim(:,:,:,5))

  Pmin_global = minval(Pg(:,:,:))
  if (Pmin_global < 0.0d0) then
    print *, 'WARNING: negative pressure detected, clipping'
    where (prim(:,:,:,5) < 1.0d-6) prim(:,:,:,5) = 1.0d-6
  end if
  
  !print *, "After CommVars"

  !print *, "Pmax =", maxval(Pg)
  !print *, "Pmin =", minval(Pg)
  !print *, "Any NaN P =", any(Pg /= Pg)

  !print *, "Any NaN U =", any(ug /= ug)

  !print *, "Any NaN V =", any(vg /= vg)

  !print *, "Any NaN W =", any(wg /= wg)

  !print *, "Any NaN rho =", any(rhog /= rhog)
  !do i =1,nx
  !  do j=1,ny
  !    do k=1,nz
  !      call PrimToCons(prim(i,j,k,:),q(i,j,k,:))
        !q0(i,j,k,1) = prim(i,j,k,1)
        !q0(i,j,k,2) = prim(i,j,k,1)*prim(i,j,k,2) 
        !q0(i,j,k,3) = prim(i,j,k,1)*prim(i,j,k,3)
        !q0(i,j,k,4) = prim(i,j,k,1)*prim(i,j,k,4)
        !q0(i,j,k,5) = prim(i,j,k,5)*GAM1I + 0.5*prim(i,j,k,1)*&
       ! & (prim(i,j,k,2)**2+prim(i,j,k,3)**2+prim(i,j,k,4)**2)
   !   enddo
   ! enddo
  !enddo

  q0 = q

  !print *, "q max =", maxval(q)
  !print *, "q min =", minval(q)
  !print *, "Any NaN q =", any(q /= q)
  !uprof_in = prim(1,:,2)
  !uprof_out = prim(nx,:,2)
  !above will be used for BC
  ! Convert Primitive to Conservative variables before passing to GPU
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        q(i,j,k,1) = prim(i,j,k,1)
        q(i,j,k,2) = prim(i,j,k,1) * prim(i,j,k,2)
        q(i,j,k,3) = prim(i,j,k,1) * prim(i,j,k,3)
        q(i,j,k,4) = prim(i,j,k,1) * prim(i,j,k,4)
        q(i,j,k,5) = prim(i,j,k,5) * GAM1I + 0.5d0 * prim(i,j,k,1) * &
                     (prim(i,j,k,2)**2 + prim(i,j,k,3)**2 + prim(i,j,k,4)**2)
      enddo
    enddo
  enddo

  q0 = q
end subroutine Initialize
