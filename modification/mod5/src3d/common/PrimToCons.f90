subroutine PrimToCons(primt,qt)
  use mod_params
  implicit none
  real(kind=8),intent(in)::primt(5)
  real(kind=8),intent(out)::qt(5)
  qt(1) = primt(1)
  qt(2) = primt(2)*primt(1)
  qt(3) = primt(3)*primt(1)
  qt(4) = primt(4)*primt(4)
  qt(5) = GAM1I * primt(5)  + 0.50d0*primt(1)*&
    & (primt(2)**2+primt(3)**2+primt(4)**2)
  return
end subroutine PrimToCons
