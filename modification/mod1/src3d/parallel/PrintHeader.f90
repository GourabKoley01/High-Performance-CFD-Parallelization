subroutine printheader(ictr) 
  use mod_spmd
  use mod_params
  implicit none
  integer, intent(in) :: ictr 
  integer :: i
  integer :: &
    vals(8) ! Date and time
  real(kind=8) :: &
    w_time ! Wall-clock time
  if(ictr .eq. 1)then
    if( masterproc ) then
      write(wunitp,*)"+~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+"
      write(wunitp,*)"| Program: ANUROOP FD 3D MPI                     |"
      write(wunitp,*)"| Rajesh Ranjan, IIT Kanpur                     |"
      write(wunitp,*)"+~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+"
      write(wunitp,*)'--------------------------------------------------'
      call date_and_time( values=vals )
      write(wunitp,1001)vals(3),vals(2),vals(1)
      1001  format(1x,"Starting date: ",i0,".",i0,".",i0," (DD.MM.CCYY)")
      write(wunitp,1002)vals(5),vals(6),vals(7)
      1002  format(1x,"Starting time: ",i0,":",i0,":",i0," (HH:MM:SS)")
      write(wunitp,*)'--------------------------------------------------'
      write(wunitp,*)'!!!!!! SIMULATION INPUTS !!!!!!'
      write(wunitp,*)'--------------------------------------------------'
      write(wunitp,*)'Re, Mach, Pr =', RE, XMACH, PRANDTL
      if (imetric_scheme .eq. 10) then
        write(wunitp,*)'metric_scheme = COMPACT'
      else if (imetric_scheme .eq. 2) then
        write(wunitp,*)'metric_scheme = SECOND'
      else
        print*,'NO SUITABLE JACOBIAN SCHEME FOUND. STOPPING!!!'
        STOP
      endif

      if (iflux_scheme .eq. 10) then
        write(wunitp,*)'flux_scheme = COMPACT WITH ORDER =', iflux_order
        write(wunitp,*)'and FILTER WITH ORDER =', ifilt_order, ' and alphaf =', alphaf
      else if (iflux_scheme .eq. 2) then
        write(wunitp,*)'flux_scheme = ROE WITH ORDER =', iroe_order
      else
        print*,'NO SUITABLE FLUX SCHEME FOUND. STOPPING!!!'
        STOP
      endif

      if (ivis_scheme .eq. 10) then
        write(wunitp,*)'viscous_scheme = COMPACT'
      else if (ivis_scheme .eq. 2) then
        write(wunitp,*)'viscous_scheme = SECOND'
      else if (ivis_scheme .eq. 1) then
        write(wunitp,*)'viscous_scheme = FIRST ORDER'
      else if (ivis_scheme .eq. 0) then
        write(wunitp,*)'WARNING: NO VISCOUS TERM COMPUTED'
      else
        print*,'NO SUITABLE VISCOUS SCHEME FOUND. STOPPING!!!'
        STOP
      endif

      if (itime_scheme .eq. 2) then
        write(wunitp,*)'time_scheme = BW2'
      elseif (itime_scheme .eq. 3) then
        write(wunitp,*)'time_scheme = RK3'
      elseif (itime_scheme .eq. 4) then
        write(wunitp,*)'time_scheme = RK4'
      else
        print*,'NO SUITABLE TIME SCHEME FOUND. STOPPING!!!'
        STOP
      endif

      if (idebug .eq. 1) then
        write(wunitp,*)'DEBUG MODE ON!!!'
      endif

      write(wunitp,*)'Nitr, dt = ',nitr, dt
      write(wunitp,*)'itrst, movie_freq, iprint_freq =',itrst, movie_freq, iprint_freq
      write(wunitp,*)'iper, jper, kper =', iper, jper,kper
      write(wunitp,*)'--------------------------------------------------'
      write(wunitp,*)'!!!!!! GRID INPUTS !!!!!!'
      write(wunitp,*)'--------------------------------------------------'
      write(wunitp,*)'Nx =', nxg, 'Ny =', nyg, 'Nz =', nzg
      write(wunitp,*)'xmax,ymax,zmax=',maxval(x),maxval(y),maxval(z)
      write(wunitp,*)'xmin,ymin,zmin=',minval(x),minval(y),minval(z)
      write(wunitp,*)'Npx =', np_x, 'Npy =', np_y, 'Npz =', np_z
      write(wunitp,*)'Noverlap =', noverlap
      write(wunitp,*)'proc nx ny nz nxt nyt nzt'
      do i=1,nprocs
        write(wunitp,*) i, ipende(i)-ipstarte(i)+1, jpende(i)-jpstarte(i)+1,kpende(i)-kpstarte(i)+1, & 
          ipend(i)-ipstart(i)+1, jpend(i)-jpstart(i)+1,kpend(i)-kpstart(i)+1 
          !write(wunitp,*) 'proc= ',i, 'nxt(nx)= ',ipend(i)-ipstart(i)+1, ipende(i)-ipstarte(i)+1,&
        !  & 'nyt(ny)= ', jpend(i)-jpstart(i)+1,jpende(i)-jpstarte(i)+1,&
        !'nzt(nz)= ', kpend(i)-kpstart(i)+1,kpende(i)-kpstarte(i)+1 
      enddo
      write(wunitp,*)'--------------------------------------------------'
      write(wunitp,*)'!!!!!! ITERATION BEGINS HERE !!!!!!'
      write(wunitp,*)"--------------------------------------------------"
      write(wunitp,*)'itr     time     min(u)    min(p)  norm(rhs)'
    endif

  else
    if( masterproc ) then
      w_time = spmd_wtime( start_wtime )
      call date_and_time( values=vals )
      write(wunitp,*)'--------------------------------------------------'
      write(wunitp,1004)vals(3),vals(2),vals(1)
      1004  format(1x,"Completion date: ",i0,".",i0,".",i0," (DD.MM.CCYY)")
      write(wunitp,1005)vals(5),vals(6),vals(7)
      1005  format(1x,"Completion time: ",i0,":",i0,":",i0," (HH:MM:SS)")
      write(wunitp,*)"--------------------------------------------------"
      write(wunitp,*)"Code execution completed"
      write(wunitp,1003)w_time
      1003  format(1x,"Code execution time: ",f0.3," seconds")
      write(wunitp,*)"--------------------------------------------------"
    end if
  endif
end subroutine
