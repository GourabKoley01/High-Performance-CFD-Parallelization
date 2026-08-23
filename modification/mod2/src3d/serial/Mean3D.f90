SUBROUTINE CalcMean3d
  use mod_params
  implicit none
  integer :: i,j,k

  !$acc parallel loop collapse(3) present(prim, primmean)
  do k = 1,nz
    do j = 1,ny
      do i = 1,nx 
        primmean(i,j,k,1) = primmean(i,j,k,1) + prim(i,j,k,1)  !rho
        primmean(i,j,k,2) = primmean(i,j,k,2) + prim(i,j,k,2)  !u
        primmean(i,j,k,3) = primmean(i,j,k,3) + prim(i,j,k,3)  !v
        primmean(i,j,k,4) = primmean(i,j,k,4) + prim(i,j,k,4)  !w
        primmean(i,j,k,5) = primmean(i,j,k,5) + prim(i,j,k,5)  !P
        primmean(i,j,k,6) = primmean(i,j,k,6) + prim(i,j,k,2)*prim(i,j,k,2) !uu
        primmean(i,j,k,7) = primmean(i,j,k,7) + prim(i,j,k,3)*prim(i,j,k,3) !vv
        primmean(i,j,k,8) = primmean(i,j,k,8) + prim(i,j,k,4)*prim(i,j,k,4) !ww
        primmean(i,j,k,9) = primmean(i,j,k,9) + prim(i,j,k,2)*prim(i,j,k,3) !uv
        primmean(i,j,k,10) = primmean(i,j,k,10) + prim(i,j,k,3)*prim(i,j,k,4) !vw
        primmean(i,j,k,11) = primmean(i,j,k,11) + prim(i,j,k,2)*prim(i,j,k,4) !uw
      enddo
    enddo
  enddo

END SUBROUTINE CalcMean3d


SUBROUTINE OutputMean3d(itr)
  use mod_params
  implicit none
  integer, intent(in) :: itr
  integer :: length, nvarst
  character(len=60) :: filename
  nvarst = nvars + 6
          
! Mean Flow
    filename='meanflow'
    length=len_trim(filename) 
    write(filename(length+1:length+1),'(A1)') '_'
    length=len_trim(filename) 
    write(filename(length+1:length+7),'(I7.7)') itrst
    length=len_trim(filename) 
    write(filename(length+1:length+1),'(A1)') '_'
    length=len_trim(filename) 
    write(filename(length+1:length+7),'(I7.7)') itr
    length=len_trim(filename) 
    write(filename(length+1:length+4),'(A4)') '.fcn'

    open(16,file=trim(filename),form='unformatted',convert='big_endian')
    write(16) nx,ny,nz,nvars
    write(16) primmean(:,:,:,2),primmean(:,:,:,3),primmean(:,:,:,4),primmean(:,:,:,5),primmean(:,:,:,1)
    close(16)
    print*,filename, 'written'


! Mean correlation
    filename='meancorrel'
    length=len_trim(filename) 
    write(filename(length+1:length+1),'(A1)') '_'
    length=len_trim(filename) 
    write(filename(length+1:length+7),'(I7.7)') itrst
    length=len_trim(filename) 
    write(filename(length+1:length+1),'(A1)') '_'
    length=len_trim(filename) 
    write(filename(length+1:length+7),'(I7.7)') itr
    length=len_trim(filename) 
    write(filename(length+1:length+4),'(A4)') '.fcn'

    open(16,file=trim(filename),form='unformatted',convert='big_endian')
    write(16) nx,ny,nz,nvars
    write(16) primmean(:,:,:,6),primmean(:,:,:,7),primmean(:,:,:,8),primmean(:,:,:,9),primmean(:,:,:,10), primmean(:,:,:,11)
    close(16)
    print*,filename, 'written'
END SUBROUTINE OutputMean3d

