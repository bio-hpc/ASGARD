! <compile=optimized>
#include "../include/dprec.fh"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ bond energy place holder
subroutine bond(eb)

   implicit none

   ! passed variables
   _REAL_ eb
   
   eb = 0.d0
   

end subroutine bond 
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ angle energy place holder
subroutine angl(eba)

   implicit none

   ! passed variables
   _REAL_ eba

   eba = 0.d0


end subroutine angl
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ dihedral energy place holder
subroutine ephi(ep,enbp,eelp)

   implicit none

   ! passed variables
   _REAL_   ep,enbp,eelp

   enbp = 0.d0
   eelp = 0.d0
   ep   = 0.d0


end subroutine ephi
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ capwat restraining force place holder
subroutine capwat()

   implicit none


end subroutine capwat
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ position restraining force place holder
subroutine xconst(econ)

   implicit none

   ! passed variables
   _REAL_ econ

   econ = 0.d0


end subroutine xconst
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ belly restraining force place holder
subroutine bellyf()

   implicit none


end subroutine bellyf
