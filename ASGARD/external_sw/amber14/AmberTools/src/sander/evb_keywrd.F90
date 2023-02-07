! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Check EVB keywords                                                     |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#include "../include/dprec.fh"

   subroutine evb_keywrd ( dia_type, xch_type, evb_dyn ) 

   use evb_allowed

   implicit none

#  include "parallel.h"

   character(512), intent(in) :: dia_type, xch_type, evb_dyn

   !  ..........................................................................

   integer :: n, lfound 

!  +---------------------------------------------------------------------------+
!  |  Check diabatic state type inputs                                         |
!  +---------------------------------------------------------------------------+

   lfound = 0

   do n = 1, ndia_type
      if( trim( adjustl( dia_type_key(n) ) ) == &
          trim( adjustl( dia_type        ) )      ) lfound = 1
   enddo

   if( lfound == 0 ) then
      write(6,'(A)') 'ERROR:: EVB keyword ' // trim( adjustl( dia_type ) ) &
                  // ' has no meaning'
      write(6,'(A)') 'Permitted keywords for dia_type are :: '
      do n = 1, ndia_type
         write(6,'(A)') '   ' // trim( adjustl( dia_type_key(n) ) )
      enddo
      write(6,'(A)')
      call mexit(6,1)
   endif
 
!  +---------------------------------------------------------------------------+
!  |  Check exchange type inputs                                               |
!  +---------------------------------------------------------------------------+

   lfound = 0

   do n = 1, nxch_type
      if( trim( adjustl( xch_type_key(n) ) ) == &
          trim( adjustl( xch_type        ) )      ) lfound = 1
   enddo

   if( lfound == 0 ) then 
      write(6,'(A)') 'ERROR:: EVB keyword ' // trim( adjustl( xch_type ) ) &
                  // ' has no meaning'
      write(6,'(A)') 'Permitted keywords for xch_type are :: '
      do n = 1, nxch_type
         write(6,'(A)') '   ' // trim( adjustl( xch_type_key(n) ) )
      enddo 
      write(6,'(A)') 
      call mexit(6,1)

   endif 

!  +---------------------------------------------------------------------------+
!  |  Check EVB dynamics inputs                                                |
!  +---------------------------------------------------------------------------+

   lfound = 0

   do n = 1, nevb_dyn
      if( trim( adjustl( evb_dyn_key(n) ) ) == &
          trim( adjustl( evb_dyn        ) )      ) lfound = 1
   enddo

   if( lfound == 0 ) then
      write(6,'(A)') 'ERROR:: EVB keyword ' // trim( adjustl( evb_dyn ) ) &
                  // ' has no meaning'
      write(6,'(A)') 'Permitted keywords for evb_dyn are :: '
      do n = 1, nevb_dyn
         write(6,'(A)') '   ' // trim( adjustl( evb_dyn_key(n) ) )
      enddo
      write(6,'(A)')
      call mexit(6,1)
   endif

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine evb_keywrd

