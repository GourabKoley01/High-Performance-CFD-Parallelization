PROGRAM make_grid
  implicit none
  real(kind=8), allocatable, dimension(:,:,:) :: x,y,z
  !real (kind=8):: xmin, xmax, ymin, ymax, zmin, zmax 
  integer :: nx, ny, nz
  integer :: i, j, k

!  print*,'nprocs= ', nprocs
    open(10, file='grid3d.grd',form='UNFORMATTED',status='OLD',convert='BIG_ENDIAN')
    read(10) nx,ny,nz
    print*,'nx,ny,nz=',nx,ny,nz
    allocate(x(nx,ny,nz))
    allocate(y(nx,ny,nz))
    allocate(z(nx,ny,nz))
    read(10) x,y,z
    close(10)
    print*,'xmax,ymax,zmax=',maxval(x),maxval(y),maxval(z)
    print*,'xmin,ymin,zmin=',minval(x),minval(y),minval(z)

  open(15,file='grid3d_sp.xyz',form='UNFORMATTED',convert='BIG_ENDIAN')
  write(15) nx,ny,nz
  write(15) real(x),real(y),real(z)
  close(15)
  end
