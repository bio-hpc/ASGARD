#include "../include/dprec.fh"

! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  KAPPA_KECK: reactive flux via the Keck prescription            |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   subroutine kappa_keck ( p_ndx, RC_ndx_all )

   use wigner, only: rflux

   implicit none

   integer, intent(in) :: p_ndx, RC_ndx_all

   !  +---------------------------------------------------------------+

   integer :: n, nRC
   _REAL_  :: RC_shifted(RC_ndx_all)
   _REAL_  :: RC_reactant, RC_product, recross, alpha, diff
  
!  +---------------------------------------------------------------+
!  |  Shift all RC value with respect to dividing surface:         |
!  |    RC = 0    @ dividing surface                               |    
!  |    RC < 0    reactant                                         |
!  |    RC > 0    product                                          |
!  +---------------------------------------------------------------+

   nRC = 0

   do n = 1, RC_ndx_all

      diff = rflux%RC_sampled(n) - rflux%RC_s0

      if( diff /= 0.0d0 ) then

         nRC = nRC + 1
         RC_shifted(nRC) = diff

      endif

   enddo 

   write(555,*)
   write(555,*) 'p_ndx = ', p_ndx

   do n = 1, nRC
      write(555,*) n, RC_shifted(n)
   enddo 

!  +---------------------------------------------------------------+
!  |  xi(p,q;n) = 1/alpha  if (p,q;n) has alpha forward crossing   |
!  |                       and (1-alpha) reverse crossing          |
!  |                                                               |
!  |            = 0        otherwise                               |
!  +---------------------------------------------------------------+

   RC_reactant = RC_shifted(1)
   RC_product  = RC_shifted(nRC)

   rflux%xi(p_ndx) = 0.0d0

   write(666,*)
   write(666,*) 'p_ndx = ', p_ndx
   write(666,*) RC_reactant, RC_product

   if( RC_reactant < 0.0d0 .and. RC_product > 0.0d0 ) then

      alpha = 0.0d0

      do n = 2, nRC

         recross = RC_shifted(n-1) * RC_shifted(n)

         if( recross < 0.0d0 ) alpha = alpha + 1.0d0

          write(666,*) n, RC_shifted(n-1), RC_shifted(n), recross, alpha

      enddo

      alpha = 0.50 * ( alpha + 1.0d0 )

      rflux%xi(p_ndx) = 1.0d0 / alpha

   endif

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine kappa_keck


