! <compile=optimized>
#include "../include/dprec.fh"

module genborn

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ generalized Born nonbonded routine place holder
subroutine egb(epol, eelt, evdw, esurf)

   implicit none
  
   ! passed variables
   _REAL_ epol, eelt, evdw, esurf
 
   ! local variables   
   
   epol = 0.0d0
   eelt = 0.0d0
   evdw = 0.0d0
   esurf = 0.0d0
 
 
end subroutine egb 

end module genborn
