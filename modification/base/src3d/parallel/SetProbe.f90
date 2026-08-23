SUBROUTINE SetProbe
  use mod_params
  use mod_spmd
  implicit none
  integer:: ipr,ict
  integer:: ib,jb,kb,length
  logical :: fileexists

  open(10, file='probe_input.in',status='OLD', action='READ')
  read(10,*) nprobesg
  allocate(probeg(nprobesg,3))
    do ipr =1, nprobesg
      read(10,*) probeg(ipr,1), probeg(ipr,2), probeg(ipr,3)
    enddo
  close(10)

  do ipr=1,nprobesg
    ib = probeg(ipr,1)
    jb = probeg(ipr,2)
    kb = probeg(ipr,3)

    if(((ib .ge. isglobal) .and. (ib .le. ieglobal)))then 
    if(((jb .ge. jsglobal) .and. (jb .le. jeglobal)))then 
    if(((kb .ge. ksglobal) .and. (kb .le. keglobal)))then 
    nprobes = nprobes+1
  endif
endif
endif
enddo

if(nprobes .gt. 0)then
  allocate(probe(nprobes,3))
endif

print*,'------------------------------------------------------------'
print*, 'myrank = ', myrank,'nx,ny = ', nx, ny, 'nprobes = ', nprobes
print*, 'global=', isglobal, ieglobal, jsglobal, jeglobal
print*, 'local =', isgrid, iegrid, jsgrid, jegrid
if(nprobes .gt. 0)then
  !print*, 'local ij = 1 ', nx, ' 1 ', ny
  ict=0;
  do ipr=1,nprobesg
    ib = probeg(ipr,1)
    jb = probeg(ipr,2)
    kb = probeg(ipr,3)
    if(((ib .ge. isglobal) .and. (ib .le. ieglobal)))then 
    if(((jb .ge. jsglobal) .and. (jb .le. jeglobal)))then 
    if(((kb .ge. ksglobal) .and. (kb .le. keglobal)))then 
    !print*,'Global Block # ',ipr
    !print*,'global= ', isblkg(ipr), ieblkg(ipr), jsblkg(ipr), jeblkg(ipr)
    !print*,'global= ', max(isblkg(ipr),isglobal), min(ieblkg(ipr),ieglobal),&
    !  & max(jsblkg(ipr),jsglobal), min(jeblkg(ipr),jeglobal)
    ict = ict+1
    !ibs1 = (ibs - isglobal)+1 
    !ibe1 = (ibe - isglobal)+1
    !jbs1 = (jbs - jsglobal)+1
    !jbe1 = (jbe - jsglobal)+1
    !ib1 = (ib - isglobal)+isgrid 
    !jb1 = (jb - jsglobal)+jsgrid

    !probe(ict,1) = max(ib1,isgrid)
    !probe(ict,2) = max(jb1,jsgrid)

    probe(ict,1) = (ib - isglobal)+isgrid 
    probe(ict,2) = (jb - jsglobal)+jsgrid
    probe(ict,3) = (kb - ksglobal)+ksgrid
    !take care of the overlap region
!    if(isblk(ict) .eq. isgrid)isblk(ict) = 1
!    if(ieblk(ict) .eq. iegrid)ieblk(ict) = nx
!    if(jsblk(ict) .eq. jsgrid)jsblk(ict) = 1
!    if(jeblk(ict) .eq. jegrid)jeblk(ict) = ny
    !isblk(ict) = max(ibs1,1)
    !ieblk(ict) = min(ibe1,nx)
    !jsblk(ict) = max(jbs1,1)
    !jeblk(ict) = min(jbe1,ny)

    !isblk(nblnks) = ibs1 - noverlap
    !ieblk(nblnks) = ibe1 - noverlap
    !jsblk(nblnks) = jbs1 - noverlap
    !jeblk(nblnks) = jbe1 - noverlap
    pfile(ict)='probe'
    length=len_trim(pfile(ict)) 
    write(pfile(ict)(length+1:length+1),'(A1)') '_'
    length=len_trim(pfile(ict)) 
    write(pfile(ict)(length+1:length+7),'(I2.2)') ipr
    length=len_trim(pfile(ict)) 
    write(pfile(ict)(length+1:length+5),'(A4)') '.txt'

  INQUIRE( FILE=trim(pfile(ict)), EXIST=fileexists ) 
  IF (.not. fileexists )then 
    open(200+ict,file=trim(pfile(ict)))
    !write(16,"(a,I4,2X,a,I4)") 'I=',probe(ict,1),"J=",probe(ict,2) 
    !write(200+ict,"('#', a,I4,2X,a,I4)") 'I=',ib,"J=",jb 
    !write(200+ict,"(a,F4.2,2X,a,F4.2)") 'x=',x(probe(ict,1),probe(ict,2)),"y=",y(probe(ict,1),probe(ict,2)) 
 !   write(200+ict,10) ib,jb 
    write(200+ict,"('# I, J, K', 3(I4,2X))") ib,jb, kb 
!10	format('# I, J', 2(I4,2X)) 
    write(200+ict,"('# X, Y, Z', 3(F8.5,2X))") x(probe(ict,1),probe(ict,2),probe(ict,3)), &
            y(probe(ict,1),probe(ict,2),probe(ict,3)), z(probe(ict,1),probe(ict,2),probe(ict,3))
    write(200+ict,"('# ITR, TIME, U, P')") 
    !write(200+ict,"(a,F4.2,2X,a,F4.2)") 'x=',x(probe(ict,1),probe(ict,2)),"y=",y(probe(ict,1),probe(ict,2)) 
   ELSE
    open(200+ict,file=trim(pfile(ict)),position='append')
    !close(200+ict)
   ENDIF
  endif
endif
endif
  enddo
  !do ipr=1,nblnks
  !  !print*,'block= ', ipr
  !  print*,'Local Block # ',ipr
  !  print*,'local= ', isblk(ipr), ieblk(ipr), jsblk(ipr), jeblk(ipr)
  !enddo

  !do ipr=1,nblnks
  !  ibs = isblk(ipr)
  !  ibe = ieblk(ipr)
  !  jbs = jsblk(ipr)
  !  jbe = jeblk(ipr)
  !  do i=ibs,ibe
  !    do j=jbs,jbe
  !      ifblock(i,j) = 1
  !    enddo
  !  enddo
  !enddo
endif

END SUBROUTINE SetProbe
