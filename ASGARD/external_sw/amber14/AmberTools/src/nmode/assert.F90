#include "assert.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Assertion failure reporter

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine aass here]
subroutine aass( condition, file, line )

   implicit none
   character(*) condition
   character(*) file
   integer      line

   write( 6,'(5A,I6,A)') 'ASSERTion ''', condition, ''' failed in ', &
         file, ' at line ', line, '.'
   call mexit(6,1)
   return
end subroutine aass
