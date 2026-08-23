subroutine UserInputs
  use mod_params
  implicit none
  real(kind=8) :: alpha_f
  real (kind=8) :: Linf, ainf, xmuinf
  integer :: i

  open(10, file='solver_input.in',status='OLD', action='READ')
  read(10,*) iSolver 
  read(10,*) gridFile 
  read(10,*) RE 
  read(10,*) XMACH
  read(10,*) TINF
  read(10,*) dt
  read(10,*) itrst
  read(10,*) nitr
  read(10,*) movie_freq
  read(10,*) iprint_freq
  read(10,*) iper, jper, kper
  read(10,*) imetric_scheme 
  read(10,*) iflux_scheme 
  read(10,*) iflux_order
  read(10,*) ifilt_order
  read(10,*) alpha_f
  read(10,*) iroe_order
  read(10,*) ivis_scheme 
  read(10,*) itime_scheme
  read(10,*) imean
  read(10,*) iDebug
  read(10,*) kClip, kPlane, movie2d_freq
  read(10,*) nblnksg
  if(nblnksg .gt. 0)then
    allocate(isblkg(nblnksg))
    allocate(ieblkg(nblnksg))
    allocate(jsblkg(nblnksg))
    allocate(jeblkg(nblnksg))
    allocate(ksblkg(nblnksg))
    allocate(keblkg(nblnksg))
    do i =1, nblnksg
      read(10,*) isblkg(i), ieblkg(i), jsblkg(i), jeblkg(i), ksblkg(i), keblkg(i)
    enddo
  endif
  close(10)

  if(iSolver .eq. 0)then
          SMU = SMUD
          SMUP1 = 1.458e-6
          GXM = 1.0/Rgas
          REI=1.0d0
          FACTOR=1.0d0
          xmuinf = Tinf**1.5d0*SMUP1/(Tinf+SMU)
          Linf = 1.0d0
          ainf = sqrt(GAM*Rgas*Tinf)
          uinf = XMACH*ainf
          vinf = 0.0d0
          winf = 0.0d0
          rhoinf = RE*xmuinf/(uinf*Linf)
          Pinf = rhoinf*Tinf/GXM
  else
          SMU = SMUD/Tinf
          SMUP1=SMU+1.0d0
          GXM = GAM*XMACH**2 
          REI = 1.0d0 / RE 
          !Factor for energy equation
          !G1PRM2I = GAM1I/(PR*XMACH**2) 
          FACTOR = GAM1I/(PRANDTL*XMACH**2) 
          rhoinf = 1.0d0
          uinf = 1.0d0
          vinf = 0.0d0
          winf = 0.0d0
          pinf = 1.0d0/GXM
  endif


  if(iflux_scheme .eq. 10)call filter_coeff(alpha_f)
end subroutine UserInputs

subroutine filter_coeff(alpha_f)
  use mod_params
  implicit none
  real(kind=8), intent(in) :: alpha_f
!  integer :: i,j
!  integer :: ifilt_array
  !!Initialize Filter Parameters
  afilt = 0.0d0
  !10th Order
  alphaf(6) = alpha_f
  afilt(6,1) = (193.0d0 + 126.0d0*alphaf(6)) / 256.0d0
  afilt(6,2) = (105.0d0 + 302.0d0*alphaf(6)) / 256.0d0
  afilt(6,3) = (-15.0d0 + 30.0d0*alphaf(6)) / 64.0d0
  afilt(6,4) = (45.0d0 - 90.0d0*alphaf(6))/ 512.0d0 
  afilt(6,5) = (-5.0d0 + 10.0d0*alphaf(6))/ 256.0d0 
  afilt(6,6) = (1.0d0 - 2.0d0*alphaf(6))/ 512.0d0 

  !8th Order
  alphaf(5) = alpha_f
  afilt(5,1) = (93.0d0 + 70.0d0*alphaf(5)) / 128.0d0
  afilt(5,2) = (7.0d0 + 18.0d0*alphaf(5)) / 16.0d0
  afilt(5,3) = (-7.0d0 + 14.0d0*alphaf(5)) / 32.0d0
  afilt(5,4) = 1.0d0 / 16.0d0 - alphaf(5)/ 8.0d0 
  afilt(5,5) = -1.0d0 / 128.0d0 + alphaf(5)/ 64.0d0 

  !6th Order
  alphaf(4) = alpha_f
  afilt(4,1) = 11.0d0 / 16.0d0 + 5.0d0*alphaf(4)/ 8.0d0 
  afilt(4,2) = 15.0d0 / 32.0d0 + 17.0d0*alphaf(4)/ 16.0d0 
  afilt(4,3) = -3.0d0 / 16.0d0 + 3.0d0*alphaf(4)/ 8.0d0 
  afilt(4,4) = 1.0d0 / 32.0d0 - alphaf(4)/ 16.0d0 

  !4th Order
  alphaf(3) = alpha_f
  afilt(3,1) = 5.0d0 / 8.0d0 + 3.0d0*alphaf(3)/ 4.0d0 
  afilt(3,2) = 0.50d0 + alphaf(3) 
  afilt(3,3) = -1.0d0 / 8.0d0 + alphaf(3)/ 4.0d0 

  !2nd Order
  alphaf(2) = alpha_f
  afilt(2,1) = 0.50d0 + alphaf(2) 
  afilt(2,2) = 0.50d0 + alphaf(2) 

  !2nd Order
  alphaf(1) = 0.0d0
  afilt(1,1) = 1.0d0 
  !multiplication by half done already here (eq. 2.2 of Gaitonde & Visbal FDL
  !report)

!  do i=1,6
!    do j=1,6
!      print*,i,j,afilt(i,j)
!    enddo
!  enddo
  afilt = 0.50d0*afilt
end subroutine filter_coeff
