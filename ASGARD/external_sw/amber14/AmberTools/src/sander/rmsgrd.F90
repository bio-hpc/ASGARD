! <compile=optimized>
#include "copyright.h"
#include "../include/assert.fh"
#include "../include/dprec.fh"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Calculate the RMS gradient
!-----------------------------------------------------------------------
!     --- RMSGRD ---
!-----------------------------------------------------------------------

subroutine rmsgrd(forces, grms)
   
   use constants, only : ZERO
   use poisson_boltzmann, only : outwat, oution
   
   implicit none
   
#include "../include/md.h"
#include "box.h"
#include "../include/memory.h"

   _REAL_, intent(in) :: forces(*)
   _REAL_, intent(out) :: grms
   
   ! the ddot function is external
   _REAL_ ddot
   
   _REAL_ :: dotprod
   
   integer :: numcomponents !, shakecomponents
   
   if (ibelly > 0) then
      numcomponents = natbel * 3
   else
      numcomponents = nrp * 3
   end if
   
   ! Ben Roberts: Took out this block on 14 April 2011.
   ! Need to find out whether SHAKEn atoms should be included in
   ! RMS gradient calculations or not.
   !shakecomponents = 0
   !if (ntc == 2) then
   !   shakecomponents = nbonh
   !else if (ntc == 3) then
   !   shakecomponents = nbonh + nbona
   !end if
   !numcomponents = numcomponents - shakecomponents
   
   if (ifcap == 2 .or. ifcap == 5) then
      numcomponents = numcomponents - 3*(outwat + oution)
   end if
      
   ! Initialise rmsgrad so it at least has a valid value
   grms = ZERO
   dotprod = ddot(3*nrp, forces, 1, forces, 1)
   if (numcomponents /= 0) grms = sqrt(dotprod / numcomponents)
   
end subroutine rmsgrd


