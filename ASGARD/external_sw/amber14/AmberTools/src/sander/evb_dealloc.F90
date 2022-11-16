! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  EVB_DEALLOC                                                            |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#include "../include/dprec.fh"

   subroutine evb_dealloc

   use evb_parm,  only: input_ready 
   use evb_amber, only: evb_amber_ready
   use evb_xchff, only: xchff_warshel_ready, xchff_gauss_ready
   use schlegel,  only: init_ready
   use evb_data,  only: evb_Hmat, evb_frc_ready
   use evb_check, only: full_evb_debug
   use file_io_dat

   implicit none

#  include "parallel.h"

   if( full_evb_debug ) write(6,'(A)') '|'

!  +---------------------------------------------------------------------------+
!  |  Deallocate memory for data arrays in evb_input                           |
!  |  :: Allocation done in evb_input ::                                       |  
!  +---------------------------------------------------------------------------+

   if( input_ready ) then
      call evb_input_dealloc
      input_ready = .false. 
      if( full_evb_debug ) write(6,'(A)') '| DONE deallocation for evb_input'
   endif 

!  +---------------------------------------------------------------------------+
!  |  Deallocate memory for xdat_dia and xdat_xch in evb_init                  |
!  |  :: Allocation done in evb_init ::                                        |
!  +---------------------------------------------------------------------------+

   if( init_ready ) then
      call evb_init_dealloc
      init_ready = .false.
      if( full_evb_debug ) write(6,'(A)') '| DONE deallocating for evb_init'
   endif

!  +---------------------------------------------------------------------------+
!  |  Deallocate memory for xch_warshel in exchange_warshel                    |
!  +---------------------------------------------------------------------------+

   if( xchff_warshel_ready ) then
      call xch_warshel_dealloc
      xchff_warshel_ready = .false.
      if( full_evb_debug ) &
         write(6,'(A)') '| DONE deallocating for exchange_warshel'
   endif

!  +---------------------------------------------------------------------------+
!  |  Deallocate memory for xch_gauss in exchange_gauss                        |
!  +---------------------------------------------------------------------------+

   if( xchff_gauss_ready ) then
      call xch_gauss_dealloc
      xchff_gauss_ready = .false.
      if( full_evb_debug ) &
         write(6,'(A)') '| DONE deallocating for exchange_gauss'
   endif

!  +---------------------------------------------------------------------------+
!  |  Deallocate memory for evb_Hmat in evb_matrix                             |
!  +---------------------------------------------------------------------------+

   if( evb_Hmat%ready ) then
      call evb_mat_type_dealloc
      evb_Hmat%ready = .false.
      if( full_evb_debug ) write(6,'(A)') '| DONE deallocating for evb_matrix'
   endif 

!  +---------------------------------------------------------------------------+
!  |  Deallocate memory for evb_frc in evb_force                               |
!  +---------------------------------------------------------------------------+

   if( evb_frc_ready ) then
      call evb_frc_type_dealloc
      evb_frc_ready = .false.
      if( full_evb_debug ) write(6,'(A)') '| DONE deallocating for evb_force'
   endif

!  +---------------------------------------------------------------------------+
!  |  Deallocate memory for xnrg, xf, and xq in force                          |
!  +---------------------------------------------------------------------------+

   if( evb_amber_ready ) then
      call evb_amber_dealloc
      evb_amber_ready = .false.
      if( full_evb_debug ) write(6,'(A)') '| DONE deallocating for evb_amber'
   endif

!  +---------------------------------------------------------------------------+
!  |  Close EVB output file unit                                               |
!  +---------------------------------------------------------------------------+

   if( worldrank == 0 ) then
      write(6,'(A)') '| Closing evb_unit associated with file ' &
                   // trim( adjustl( evbout ) )
      close( evb_unit )
   endif

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine evb_dealloc



