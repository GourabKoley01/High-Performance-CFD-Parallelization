program main3d
        use mod_spmd
        use mod_params
        use openacc
        use nvtx
        implicit none
        integer :: ngpus, igpu
        integer :: ip,itr, ivar, imeanctr
        logical :: fileexists
        real(kind=8) :: Pmin
        integer:: i,j,k
        real(kind=8) :: ttnan
        real(kind=4) :: startw_time, endw_time

        call spmd_init
        start_wtime = spmd_wtime( zero )
        !  print*,myrank, ': I am in'
        !  print*,'!!!!!! SIMULATION INPUTS !!!!!!'
        !call cpu_time(start_wtime)
        
        if(masterproc)call system("mkdir -p snapshots snapshots2D")
        
        call UserInputs
        call GridInput
        !if(masterproc)call GridInput
        call SplitGrid
        call AllocateVar
        call initialize
        call metric3d
        !call read_grid
        !call allocateVar
        call PrintHeader(1)

        !!$acc data copy(prim,q)
        !call applybc
        !!$acc end data
        
        nprobes = 0
        INQUIRE( FILE='probe_input.in', EXIST=fileexists )
        IF ( fileexists ) CALL SetProbe
        
        do ivar=1,nvars
          call Communicate(prim(:,:,:,ivar),prim(:,:,:,ivar))
        enddo
        
        ngpus =acc_get_num_devices(acc_device_nvidia)
        igpu =  mod(myrank,ngpus)
        !print*, "myrank = ", myrank, "igpu = ", igpu, "ngpus = ", ngpus
        call acc_set_device_num(igpu, acc_device_nvidia)

        !if (itrst .eq. 0) call Output3D(itrst)
        if (itrst .eq. 0) call OutputMovie3D(itrst)
        !if (itrst .eq. 0) call Output3D(itrst)
        ! Ouput is single precision
        !print*,'in comm2', myrank, minval(prim(:,:,1)),maxval(prim(:,:,1))
        
        if(masterproc)then
                print*,'                               '
                print*,nitr                               
                print*,'!!!!!! ITERATION BEGINS HERE !!!!!!'
                print*,'                               '
                print*,'itr     time     min(u)    min(p)  norm(rhs)'
        endif
        
        call printscreen(itrst)
        if(imean .eq. 1) then
                primmean = 0.0d0
                imeanctr = 0
        endif

        !$acc data create(q0) copyin(prim, rhs, q, ifblock) copy(primmean)
        !$acc enter data copyin(xix,etax,ztax,xiy,etay,ztay,xiz,etaz,ztaz,xjac)
        !$acc enter data create(T, xmu, temp)
        if(iflux_scheme .eq. 2) then
                !$acc enter data create(primL, primR, xiflux, etflux, ztflux)
        elseif(iflux_scheme .eq. 10) then
                !$acc enter data create(atdg_field, btdg_field, ctdg_field, gamtdg_field, rhstdg_field)
                !$acc enter data copyin(afilt,alphaf)
                !!$acc declare copyin(C6, AC4B, C3B)
        endif
        
        if(masterproc)call cpu_time(startw_time) 
        !Iteration begins here
        do itr = itrst+1,itrst+nitr
          !do i=1,nx
          !do j=1,ny
          !!do k=1,nz
          !ttnan =  rhs(i,j,1,1)
          !if (isnan(ttnan)) stop '"rhs" is a NaN'
          !! enddo
          !enddo
          !enddo

          call nvtxStartRange("iterate")
          call iterate(itr)
          call nvtxEndRange

          if(mod(itr,iprint_freq) .eq. 0) then
                !$acc update self(prim,rhs)
                call printscreen(itr)
          endif
          if(mod(itr,movie_freq) .eq. 0) then
                call OutputMovie3D(itr)
                !call Output3D(itr)
          endif
          if(mod(itr,movie2D_freq) .eq. 0) then
                call OutputPlane2D(itr,kClip,kPlane)
          endif
          !call OutputMovie2D(itr)
          !call OutputMovie3D(itr)
          if (nprobes .ne. 0)then
                !only writing u and P in probe
                do ip=1,nprobes
                  write(200+ip,"(I7,3X,F11.7,3X,2(F11.7,3X))") itr, itr*dt,&
                        prim(probe(ip,1),probe(ip,2),probe(ip,3),2), &
                        !prim(probe(ip,1),probe(ip,2),3), &
                        !prim(probe(ip,1),probe(ip,2),4), &
                        prim(probe(ip,1),probe(ip,2),probe(ip,3),5)
                enddo
          endif

          if(imean .eq. 1) then
                imeanctr = imeanctr +1
                primmean(:,:,:,1:5) = primmean(:,:,:,1:5) + prim
                primmean(:,:,:,6) = primmean(:,:,:,6) + prim(:,:,:,2)*prim(:,:,:,2) !uu
                primmean(:,:,:,7) = primmean(:,:,:,7) + prim(:,:,:,3)*prim(:,:,:,3) !vv
                primmean(:,:,:,8) = primmean(:,:,:,8) + prim(:,:,:,4)*prim(:,:,:,4) !ww
                primmean(:,:,:,9) = primmean(:,:,:,9) + prim(:,:,:,2)*prim(:,:,:,3) !uv
                primmean(:,:,:,10) = primmean(:,:,:,10) + prim(:,:,:,3)*prim(:,:,:,4) !vw
                primmean(:,:,:,11) = primmean(:,:,:,11) + prim(:,:,:,2)*prim(:,:,:,4) !uw
          endif
        
          INQUIRE( FILE='stoprun', EXIST=fileexists )
          if ( fileexists ) then
            if(masterproc) print*,'STOP RUN REQUESTED'
            call Output3D(itr)
            exit
          endif
        enddo
        !Iteration ends here
        if(masterproc)call cpu_time(endw_time) 

        if(masterproc) then
              print*, " "
              print*, "Iteration time = ", endw_time-startw_time, "seconds" 
        endif
        
        !$acc end data
        !$acc exit data delete(xix,etax,ztax,xiy,etay,ztay,xiz,etaz,ztaz,xjac)

        if(imean .eq. 1) then
                primmean = primmean / imeanctr
                call OutputMean3D(itr)
                print*, 'mean collected from ',itrst, 'to ', itrst+imeanctr, 'iterations'
        endif
        
        !if(mod(itr,isave_freq) .eq. 0) then
        !  call Output3D(itr-1)
        !endif

        !    if(imean .eq. 1) then
        !      imeanctr = imeanctr +1
        !      primmean = primmean + prim
        !    endif
        !  Pmin=minval(prim(isgrid:iegrid,jsgrid:jegrid,ksgrid:kegrid,5))
        !Pmin=minval(prim(:,:,:,5))
        !    print*,myrank,Pmin
        !    print 100, itr,minval(rho),minval(u), minval(P), minval(rhs)

        !    if(mod(itr,movie_freq) .eq. 0) then
        !    call output(itr)
        !  endif
        !enddo
        !call Output3D(itr)
        !if(imean .eq. 1) then
        !  primmean = primmean / imeanctr
        !  call OutputMean2D
        !  print*, 'mean collected from ',itrst, 'to ', itrst+imeanctr, 'iterations'   
        !endif
        !call PrintScreen(itrst) 
        !  print*,'max(u), min(u)',maxval(u), minval(u)
        !  print*,'max(p), min(p)',maxval(P), minval(P)
        !  print*,'                               '
        !  print*,'!!!!!! JACOBIAN VALUES, USING',imetric_scheme, '-ND ORDER SCHEME!!!!!!'
        !  call metric3d
        !  print*,'!!!!!! ITERATION BEGINS HERE1 !!!!!!'
        call PrintHeader(2)
        call spmd_final
end program main3d
