subroutine gridInput
  use mod_params
  implicit none
  character(len=20) :: filename
  character(len=20) :: path='./snapshots/'

!  print*,'nprocs= ', nprocs
    !open(10, file='grid3d.grd',form='UNFORMATTED',status='OLD',convert='BIG_ENDIAN')
    open(10, file=trim(gridFile),form='UNFORMATTED',status='OLD',convert='BIG_ENDIAN')
!    read(10) nblks
    read(10) nx,ny,nz
    print*,'nx,ny,nz=',nx,ny,nz
    allocate(x(nx,ny,nz))
    allocate(y(nx,ny,nz))
    allocate(z(nx,ny,nz))
    read(10) x,y,z
    close(10)
    !print*,'npx, npy, npz=', np_x, np_y, np_z
    print*,'xmax,ymax,zmax=',maxval(x),maxval(y),maxval(z)
    print*,'xmin,ymin,zmin=',minval(x),minval(y),minval(z)
 ! endif
    filename='gridsp.xyz'
    open(16, file=trim(path)//trim(filename),form='UNFORMATTED',convert='BIG_ENDIAN')
    write(16) nx,ny,nz
    write(16) real(x),real(y),real(z)
    close(16)
    allocate(ifblock(nx,ny,nz))
    ifblock = 0
end subroutine gridInput
