! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Output data read from the Gaussian formatted checkpoint file           |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#include "../include/assert.fh"
#include "../include/dprec.fh"

   subroutine bmat_debug ( bmat, bmat_bond, bmat_angle, bmat_dihed &
                         , binv, dbmat, ioe )

   use schlegel, only: ncoord, natm, nbond, nangle, ndihed

   implicit none

#  include "parallel.h"

   integer, intent(in) :: ioe

   _REAL_, intent(in) :: bmat(ncoord,natm*3), bmat_bond(nbond,natm*3) &
                       , bmat_angle(nangle,natm*3), bmat_dihed(ndihed,natm*3) &
                       , binv(natm*3,ncoord), dbmat(ncoord,natm*3,natm*3)

   !  ..........................................................................

   integer :: m, n 

   open( ioe, file = 'bmat.debug' )

!  +---------------------------------------------------------------------------+
!  |  Write B matrix contributions from bonds                                  |
!  +---------------------------------------------------------------------------+

   write( ioe, '( A ) ') '[B matrix elements from bonds]'
   write( ioe, 1000 ) ( bmat_bond(:,n), n = 1,  natm*3 )
   write( ioe, * )

!  +---------------------------------------------------------------------------+
!  |  Write B matrix contributions from angles                                 |
!  +---------------------------------------------------------------------------+

   write( ioe, '( A ) ') '[B matrix elements from angles]'
   write( ioe, 1000 ) ( bmat_angle(:,n), n = 1,  natm*3 )
   write( ioe, * )

!  +---------------------------------------------------------------------------+
!  |  Write B matrix contributions from dihedrals                              |
!  +---------------------------------------------------------------------------+

   write( ioe, '( A ) ') '[B matrix elements from dihedrals]'
   write( ioe, 1000 ) ( bmat_dihed(:,n), n = 1,  natm*3 ) 
   write( ioe, * )

!  +---------------------------------------------------------------------------+
!  |  Write B matrix contributions from dihedrals                              |
!  +---------------------------------------------------------------------------+

   write( ioe, '( A ) ') '[B matrix]'
   write( ioe, 1000 ) ( bmat(:,n), n = 1,  natm*3 )
   write( ioe, * )

!  +---------------------------------------------------------------------------+
!  |  Write inverse of B matrix                                                |
!  +---------------------------------------------------------------------------+

   write( ioe, '( A ) ') '[B^-1 matrix]'
   write( ioe, 1000 ) ( binv(:,n), n = 1,  ncoord )
   write( ioe, * )

!  +---------------------------------------------------------------------------+
!  |  Write gradient of B matrix                                               |
!  +---------------------------------------------------------------------------+

   write( ioe, '( A ) ') '[d/dx B matrix]'
   do m = 1, natm*3
   write( ioe, 1000 ) ( dbmat(n,m,:), n = 1, ncoord )
   enddo
   write( ioe, * )

   close( ioe )

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

 1000 format( ( 5( 1PE16.8 ) ) )

   end subroutine bmat_debug

