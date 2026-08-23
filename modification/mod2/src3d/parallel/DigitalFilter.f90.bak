subroutine gridInput
  use mod_params
  use mod_spmd
  implicit none
  character(len=20) :: filename
  character(len=20) :: path='./snapshots/'

!  print*,'nprocs= ', nprocs
  if(masterproc) then
    !open(10, file='grid3d.grd',form='UNFORMATTED',status='OLD',convert='BIG_ENDIAN')
    open(10, file=trim(gridFile),form='UNFORMATTED',status='OLD',convert='BIG_ENDIAN')
!    read(10) nblks
    read(10) nxg,nyg,nzg
    print*,'nx,ny,nz=',nxg,nyg,nzg
    allocate(xg(nxg,nyg,nzg))
    allocate(yg(nxg,nyg,nzg))
    allocate(zg(nxg,nyg,nzg))
    read(10) xg,yg,zg
    close(10)
    !print*,'npx, npy, npz=', np_x, np_y, np_z
    print*,'xmax,ymax,zmax=',maxval(xg),maxval(yg),maxval(zg)
    print*,'xmin,ymin,zmin=',minval(xg),minval(yg),minval(zg)
 ! endif
    filename='gridsp.xyz'
    open(16, file=trim(path)//trim(filename),form='UNFORMATTED',convert='BIG_ENDIAN')
    write(16) nxg,nyg,nzg
    write(16) real(xg),real(yg),real(zg)
    close(16)
  endif
end subroutine gridInput
