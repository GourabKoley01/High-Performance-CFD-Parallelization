SUBROUTINE Output3d(itr)
  use mod_params
  implicit none
  integer, intent(in) :: itr
  integer :: length
  character(len=60) :: filename
  character(len=20) :: path='./snapshots/'

    filename='mov'
    length=len_trim(filename) 
    write(filename(length+1:length+1),'(A1)') '_'
    length=len_trim(filename) 
    write(filename(length+1:length+7),'(I7.7)') itr
    length=len_trim(filename) 
    write(filename(length+1:length+4),'(A4)') '.fcn'

    open(16,file=trim(path)//trim(filename),form='unformatted',convert='big_endian')
    write(16) nx,ny,nz,nvars
    write(16) prim(:,:,:,2),prim(:,:,:,3),prim(:,:,:,4),prim(:,:,:,5),prim(:,:,:,1)
    close(16)
    print*,filename, 'written'

END SUBROUTINE Output3d

SUBROUTINE OutputPlane2D(itr, iplane, kp)
  use mod_params
  implicit none
  integer, intent(in) :: itr,kp,iplane
  integer :: length
  character(len=60) :: filename
  character(len=20) :: path='./snapshots2D/'
  logical :: fileexists
    
    filename='kmov'
    length=len_trim(filename) 
    write(filename(length+1:length+1),'(A1)') '_'
    length=len_trim(filename) 
    write(filename(length+1:length+7),'(I7.7)') itr
    length=len_trim(filename) 
    write(filename(length+1:length+4),'(A4)') '.fcn'

    open(16,file=trim(path)//trim(filename),form='unformatted',convert='big_endian')
    write(16) nx,ny,nvars
    write(16) real(prim(:,:,kp,2)),real(prim(:,:,kp,3)),real(prim(:,:,kp,4)),real(prim(:,:,kp,5)),real(prim(:,:,kp,1))
    close(16)
    print*,filename, 'written'

    INQUIRE( FILE='kgrid2D.xyz', EXIST=fileexists )
    
    if (fileexists)then
    else
      filename='kgrid2D.xyz'
      open(17,file=trim(path)//trim(filename),form='unformatted',convert='big_endian')
      !write(16) nblks
      write(17) nx,ny
      write(17) real(x(:,:,kp)),real(y(:,:,kp))
      close(17)
    endif


END SUBROUTINE OutputPlane2D

