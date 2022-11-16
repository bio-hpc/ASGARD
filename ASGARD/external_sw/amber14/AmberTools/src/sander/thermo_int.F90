! <compile=optimized>
#include "copyright.h"
#include "../include/dprec.fh"
#include "../include/assert.fh"

#ifdef MPI
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Mix the two states with a lambda-based mixing
subroutine mix_frcti(f,ener,fcopy,ecopy,nr3,clambda,klambda)

   use state
   use softcore, only : ifsc, mix_frc_sc,log_dvdl

   implicit none
   integer nr3,klambda
   type(state_rec)  :: ener, ecopy
   _REAL_ f(*),fcopy(*)
   _REAL_ clambda,dvdl,w0,w1

   !     ---we are in final state:   

   w0 = (1.d0 - clambda)**klambda
   w1 = 1.d0 - w0
   if( klambda == 1 )then
     dvdl = ener%pot%tot - ecopy%pot%tot
   else
     dvdl = klambda*(ener%pot%tot - ecopy%pot%tot)*(1.d0 - clambda)**(klambda-1)
   end if
      
   ener = (ener*w1) + (ecopy*w0)

   if (ifsc == 1) then
      ! This subroutine also adds the softcore contribution to dvdl
      call mix_frc_sc(dvdl,w0,w1,f,fcopy)
   else
      f(1:nr3) = w1*f(1:nr3) + w0*fcopy(1:nr3)
   end if

   !write(6,*) "thermo_int:: assigning ener%dvdl with value of :", dvdl
   ener%pot%dvdl = dvdl

   ecopy  = ener

   if (ifsc == 0) then
      fcopy(1:nr3) = f(1:nr3)
   end if

   call log_dvdl(dvdl)      

   return
end subroutine mix_frcti 
#else

! This function should never be called in serial
subroutine mix_frcti

   REQUIRE( .false. )

end subroutine
#endif /* MPI */
