PROGRAM make_grid
  implicit none
  ! program to generate profiles of disturbances at inlet

  !  real, allocatable, dimension(:,:) :: x2d,y2d

  real(kind=8), allocatable, dimension(:,:,:) :: x,y,z
  real (kind=8):: PIVL, XM1, XM12, GAM, GXM 
  real (kind=8):: xmin, xmax, ymin, ymax, zmin, zmax 
  real (kind=8):: dx, dy, dz 
  integer :: nx, ny, nz, ibias, jbias, kbias
  integer :: i, j, k
  integer :: nblks, nfcns
  real (kind=8):: u, v, w, rho, P 

  PIVL=2.*ATAN2(1.0,0.0)
  PRINT*,'PIVL=',PIVL

  !  open(10,file='2D_grid.grd',form='UNFORMATTED',STATUS='OLD')
  !  read(10) nx,ny
  !  print*,'nx,ny=',nx,ny
  !  allocate (x2d(nx,ny))
  !  allocate (y2d(nx,ny))
  !  read(10) x2d,y2d
  !  close(10)
  !  print*,'minmax ',minval(x2d),maxval(x2d),minval(y2d),maxval(y2d)

  open(10, file='mesh_input.in',status='OLD', action='READ')
  read(10,*) xmin, xmax, nx 
  read(10,*) ymin, ymax, ny 
  read(10,*) zmin, zmax, nz 
  read(10,*) ibias, jbias, kbias 
  close(10)

  print*, xmin, xmax, nx 
  print*, ymin, ymax, ny 
  print*, zmin, zmax, nz 
  print*, ibias, jbias, kbias 
  !print*,iflux_scheme, itime_scheme, nitr
  !print*,'enter nx --- suggest 100'
  !read*,nx
  !print*,'enter (xmin,xmax)/PI   suggest 0,4'
  !read*,xmin,xmax
  !print*,'enter ny --- suggest 100'
  !read*,ny
  !print*,'enter ymin,ymax.  suggest -1,1'
  !read*,ymin,ymax
  !print*,'enter nz --- suggest 100 '
  !read*,nz
  !print*,'enter 3(zmin,zmax)/PI.  suggest 0,4'
  !read*,zmin,zmax

  !#if channel
  !print*, 'x and z values wnxl be multiplied by PI' 
  !xmax = xmax*PIVL;
  !zmax = zmax*PIVL/3.;
  !#endif

  if (ibias .eq. 0) then
    dx=(xmax-xmin)/dble(nx-1)
  endif
  if (jbias .eq. 0) then
    dy=(ymax-ymin)/dble(ny-1)
  endif
  if (kbias .eq. 0) then
    dz=(zmax-zmin)/dble(nz-1)
  endif

  allocate (x(nx,ny,nz))
  allocate (y(nx,ny,nz))
  allocate (z(nx,ny,nz))

  do i=1,nx
    do j=1,ny
      do k=1,nz
        !           x(i,j,k)=x2d(i,j)
        !           y(i,j,k)=y2d(i,j)
        x(i,j,k)=xmin+dble(i-1)*dx
        y(i,j,k)=ymin+dble(j-1)*dy
        !           y(i,j,k)=tanh(4.2*(real(j-1)*1.0/real(ny-1)-0.5))/tanh(4.2*0.5)
        z(i,j,k)=zmin+dble(k-1)*dz
      enddo
    enddo
  enddo

  print*,'minmax ',minval(x),maxval(x),minval(y),maxval(y),minval(z),maxval(z)
  open(15,file='SOD3d.xyz',form='UNFORMATTED',convert='BIG_ENDIAN')
  nblks=1
!  write(15) nblks
  write(15) nx,ny,nz
  write(15) x,y,z
  close(15)
  open(16,file='grid3d.fcn',form='UNFORMATTED',convert='BIG_ENDIAN')
  nfcns=5
 ! write(16) nblks
  write(16) nx,ny,nz,nfcns

  XM1=0.1
  XM12=XM1**2
  GAM=1.4
  GXM=GAM*XM12
  u=0.0
  v=0.0
  w=0.0
  P=1.0/GXM
  rho=1.0

  write(16) (((u,i=1,nx),j=1,ny),k=1,nz),&
    & (((v,i=1,nx),j=1,ny),k=1,nz),&
    & (((w,i=1,nx),j=1,ny),k=1,nz),&
    & (((p,i=1,nx),j=1,ny),k=1,nz),&
    & (((rho,i=1,nx),j=1,ny),k=1,nz)
  close(16)

  STOP
END PROGRAM make_grid
