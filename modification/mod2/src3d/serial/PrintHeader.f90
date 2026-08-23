subroutine printheader(ictr) 
  use mod_params
  implicit none
  integer, intent(in) :: ictr 
  integer :: i
  integer :: &
    vals(8) ! Date and time
  real(kind=8) :: &
    w_time ! Wall-clock time
      print*,"+~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+"
      print*,"| Program: ANUROOP FD 3D                         |"
      print*,"| Rajesh Ranjan, IIT Kanpur                     |"
      print*,"| rajeshranjanarya@gmail.com                     |"
      print*,"+~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+"
  if(ictr .eq. 1)then
      write(wunit,*)"+~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+"
      write(wunit,*)"| Program: ANUROOP FD 3D                         |"
      write(wunit,*)"| Rajesh Ranjan, IIT Kanpur                     |"
      write(wunit,*)"| rajeshranjanarya@gmail.com                     |"
      write(wunit,*)"+~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+"
      write(wunit,*)'--------------------------------------------------'
      call date_and_time( values=vals )
      write(wunit,1001)vals(3),vals(2),vals(1)
      1001  format(1x,"Starting date: ",i0,".",i0,".",i0," (DD.MM.CCYY)")
      write(wunit,1002)vals(5),vals(6),vals(7)
      1002  format(1x,"Starting time: ",i0,":",i0,":",i0," (HH:MM:SS)")
      write(wunit,*)'--------------------------------------------------'
      write(wunit,*)'!!!!!! SIMULATION INPUTS !!!!!!'
      write(wunit,*)'--------------------------------------------------'
      write(wunit,*)'Re, Mach, Pr =', RE, XMACH, PRANDTL
      if (imetric_scheme .eq. 10) then
        write(wunit,*)'metric_scheme = COMPACT'
      else if (imetric_scheme .eq. 2) then
        write(wunit,*)'metric_scheme = SECOND'
      else
        print*,'NO SUITABLE JACOBIAN SCHEME FOUND. STOPPING!!!'
        STOP
      endif

      if (iflux_scheme .eq. 10) then
        write(wunit,*)'flux_scheme = COMPACT WITH ORDER =', iflux_order
        write(wunit,*)'and FILTER WITH ORDER =', ifilt_order, ' and alphaf =', alphaf
      else if (iflux_scheme .eq. 2) then
        write(wunit,*)'flux_scheme = ROE WITH ORDER =', iroe_order
      else
        print*,'NO SUITABLE FLUX SCHEME FOUND. STOPPING!!!'
        STOP
      endif

      if (ivis_scheme .eq. 10) then
        write(wunit,*)'viscous_scheme = COMPACT'
      else if (ivis_scheme .eq. 2) then
        write(wunit,*)'viscous_scheme = SECOND'
      else if (ivis_scheme .eq. 1) then
        write(wunit,*)'viscous_scheme = FIRST ORDER'
      else if (ivis_scheme .eq. 0) then
        write(wunit,*)'WARNING: NO VISCOUS TERM COMPUTED'
      else
        print*,'NO SUITABLE VISCOUS SCHEME FOUND. STOPPING!!!'
        STOP
      endif

      if (itime_scheme .eq. 2) then
        write(wunit,*)'time_scheme = BW2'
      elseif (itime_scheme .eq. 3) then
        write(wunit,*)'time_scheme = RK3'
      elseif (itime_scheme .eq. 4) then
        write(wunit,*)'time_scheme = RK4'
      else
        print*,'NO SUITABLE TIME SCHEME FOUND. STOPPING!!!'
        STOP
      endif

      if (idebug .eq. 1) then
        write(wunit,*)'DEBUG MODE ON!!!'
      endif

      write(wunit,*)'Nitr, dt = ',nitr, dt
      write(wunit,*)'itrst, movie_freq =',itrst, movie_freq
      write(wunit,*)'iper, jper =', iper, jper
      write(wunit,*)'--------------------------------------------------'
      write(wunit,*)'!!!!!! GRID INPUTS !!!!!!'
      write(wunit,*)'--------------------------------------------------'
  write(wunit,*)'Nx =', nx, 'Ny =', ny
  write(wunit,*)'xmax,ymax=',maxval(x),maxval(y)
  write(wunit,*)'xmin,ymin=',minval(x),minval(y)
      write(wunit,*)'--------------------------------------------------'
      write(wunit,*)'!!!!!! ITERATION BEGINS HERE !!!!!!'
      write(wunit,*)"--------------------------------------------------"
      write(wunit,*)'itr     time     min(u)    min(p)  norm(rhs)'

  else
  call cpu_time(w_time)
 !     w_time = spmd_wtime( start_wtime )
      !w_time = spmd_wtime( wtime0)
      call date_and_time( values=vals )
      write(wunit,*)'--------------------------------------------------'
      write(wunit,1004)vals(3),vals(2),vals(1)
      1004  format(1x,"Completion date: ",i0,".",i0,".",i0," (DD.MM.CCYY)")
      write(wunit,1005)vals(5),vals(6),vals(7)
      1005  format(1x,"Completion time: ",i0,":",i0,":",i0," (HH:MM:SS)")
      write(wunit,*)"--------------------------------------------------"
      write(wunit,*)"Code execution completed"
      write(wunit,1003)w_time-start_wtime
      1003  format(1x,"Code execution time: ",f0.3," seconds")
      write(wunit,*)"--------------------------------------------------"
  endif
end subroutine
