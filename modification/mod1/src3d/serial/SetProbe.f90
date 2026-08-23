SUBROUTINE SetProbe
  use mod_params
  implicit none
  integer:: ipr
  integer:: ib,jb,kb, length
  logical :: fileexists

  open(10, file='probe_input.in',status='OLD', action='READ')
  read(10,*) nprobesg
  allocate(probeg(nprobesg,3))
  allocate(probedatag(nprobesg,iprint_freq,nvars))
    do ipr =1, nprobesg
      read(10,*) probeg(ipr,1), probeg(ipr,2), probeg(ipr,3)
    enddo
  close(10)

  ipr=0;
  do ipr=1,nprobesg
    ib = probeg(ipr,1)
    jb = probeg(ipr,2)
    kb = probeg(ipr,3)

    pfile(ipr)='probe'
    length=len_trim(pfile(ipr)) 
    write(pfile(ipr)(length+1:length+1),'(A1)') '_'
    length=len_trim(pfile(ipr)) 
    write(pfile(ipr)(length+1:length+7),'(I2.2)') ipr
    length=len_trim(pfile(ipr)) 
    write(pfile(ipr)(length+1:length+5),'(A4)') '.txt'

  INQUIRE( FILE=trim(pfile(ipr)), EXIST=fileexists ) 
  IF (.not. fileexists )then 
    open(200+ipr,file=trim(pfile(ipr)))
    write(200+ipr,"('# I, J, K', 3(I4,2X))") ib,jb,kb 
    write(200+ipr,"('# X, Y', 3(F8.5,2X))") x(probeg(ipr,1),probeg(ipr,2),probeg(ipr,3)),&
            y(probeg(ipr,1),probeg(ipr,2),probeg(ipr,3)),&
            z(probeg(ipr,1),probeg(ipr,2),probeg(ipr,3))
    write(200+ipr,"('#    ITR', 5X, 'TIME', 8X, 'U', 10X, 'V', 10X, 'W', 10X, 'P', 10X, 'Rho')") 
    !write(200+ipr,"(a,F4.2,2X,a,F4.2)") 'x=',x(probe(ipr,1),probe(ipr,2)),"y=",y(probe(ipr,1),probe(ipr,2)) 
    close(200+ipr)
   !ELSE
   ! open(200+ipr,file=trim(pfile(ipr)),position='append')
   ! close(200+ipr)
   ENDIF
  enddo

END SUBROUTINE SetProbe


SUBROUTINE WriteProbe(itr)
  use mod_params
  implicit none
  integer :: ipr, itr, itr0, wt
  
  111 format(I7,3X,6(F8.5,3X))
  itr0 = itr-iprint_freq

  do ipr = 1,nprobesg
    open(200+ipr,file=trim(pfile(ipr)),position='append')
    do wt = 1,iprint_freq
      write(200+ipr,111) itr0+wt, (itr0+wt)*dt, probedatag(ipr,wt,2), probedatag(ipr,wt,3), &
              probedatag(ipr,wt,4), probedatag(ipr,wt,5), probedatag(ipr,wt,1)
    enddo
    close(200+ipr)
  enddo
        
  !    write(200+ipr,111) itr, itr*dt, &
  !          prim(probeg(ipr,1),probeg(ipr,2),probeg(ipr,3),2), &
  !          prim(probeg(ipr,1),probeg(ipr,2),probeg(ipr,3),3), &
  !          prim(probeg(ipr,1),probeg(ipr,2),probeg(ipr,3),4), &
  !          prim(probeg(ipr,1),probeg(ipr,2),probeg(ipr,3),5), &
  !          prim(probeg(ipr,1),probeg(ipr,2),probeg(ipr,3),1)

END SUBROUTINE WriteProbe
