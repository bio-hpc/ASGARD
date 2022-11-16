! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Output data read from the Gaussian formatted checkpoint file           |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#include "../include/assert.fh"
#include "../include/dprec.fh"

   subroutine evb_io_debug2 ( qcart, grad_cart, hess_cart, mass, xdat, ioe )

   use schlegel, only: berny_type, ncoord, natm &
                     , ibond, iangle, idihed, nbond, nangle, ndihed

   implicit none

#  include "parallel.h"

   type( berny_type ) :: xdat

   integer, intent(in) :: ioe

   _REAL_, intent(in) :: qcart(natm*3), grad_cart(natm*3) &
                       , hess_cart(natm*3,natm*3), mass(natm)

   !  ..........................................................................

   integer :: m, n, nn 
   _REAL_  :: bohr2ang
   _REAL_  :: th(ncoord*(ncoord+1)/2), th_cart(natm*3*(natm*3+1)/2)

   bohr2ang = 0.529177249
   open( ioe, file = trim( adjustl( xdat%filename ) )//'.debug2' )

!  +---------------------------------------------------------------------------+
!  |  Write external EVB redundant internal data dimension                     |
!  +---------------------------------------------------------------------------+

   write( ioe, '( A )' ) '[external evb redundant internal data dimension]'
   write( ioe, '(I12)' ) ncoord 
   write( ioe, * )

!  +---------------------------------------------------------------------------+
!  |  Write external EVB Cartesian data dimension                              |
!  +---------------------------------------------------------------------------+

   write( ioe, '( A )' ) '[external evb cartesian data dimension]'
   write( ioe, '(I12)' ) natm * 3
   write( ioe, * )

!  +---------------------------------------------------------------------------+
!  |  Write internal coordinate specifications                                 |
!  +---------------------------------------------------------------------------+

   write( ioe, '( A )' ) '[redundant internal coordinates]'
   do n = 1, nbond
      write( ioe, '(4I9)' ) ibond(n,1), ibond(n,2), 0, 0
   enddo
   do n = 1, nangle
      write( ioe, '(4I9)' ) iangle(n,1), iangle(n,3), iangle(n,2), 0
   enddo
   do n = 1, ndihed
      write( ioe, '(4I9)' ) idihed(n,1), idihed(n,2), idihed(n,3), idihed(n,4)
   enddo
   write( ioe, * )

!  +---------------------------------------------------------------------------+
!  |  Write redundant internal coordinates                                     |
!  +---------------------------------------------------------------------------+

   write( ioe, '( A )' ) '[redundant internal coordinates]'
   write( ioe, 1000 ) ( xdat%q(n), n = 1, ncoord )
   write( ioe, * )

!  +---------------------------------------------------------------------------+
!  |  Write Cartesian coordinate                                               |
!  +---------------------------------------------------------------------------+

   write( ioe, '( A )' ) '[cartesian coordinates]' 
   write( ioe, 1000 ) ( qcart(n) * bohr2ang, n = 1, natm * 3 )
   write( ioe, * )

!  +---------------------------------------------------------------------------+
!  |  Write energy                                                             |
!  +---------------------------------------------------------------------------+

   write( ioe, '( A )' ) '[electronic energy]'
   write( ioe, '( 1PE22.15 )' ) xdat%v
   write( ioe, * )

!  +---------------------------------------------------------------------------+
!  |  Write redundant internal gradient                                        |
!  +---------------------------------------------------------------------------+

   write( ioe, '( A ) ') '[redundant internal gradient]' 
   write( ioe, 1000 ) ( xdat%d(n), n = 1, ncoord ) 
   write( ioe, * )

!  +---------------------------------------------------------------------------+
!  |  Write Cartesian gradient                                                 |
!  +---------------------------------------------------------------------------+

   write( ioe, '( A ) ') '[cartesian gradient]'
   write( ioe, 1000 ) ( grad_cart(n), n = 1, natm*3 )
   write( ioe, * )

!  +---------------------------------------------------------------------------+
!  |  Write internal hessian                                                   |
!  +---------------------------------------------------------------------------+

   nn = 0
   do m = 1, ncoord
      do n = 1, m
         nn = nn + 1
         th(nn) = xdat%k(m,n)
      enddo
   enddo
   write( ioe, '( A ) ') '[redundant internal hessian]'
   write( ioe, 9999 ) ( th(n), n = 1,  ncoord*(ncoord+1)/2 )
   write( ioe, * )

!  write( ioe, '( A ) ') '[redundant internal hessian]'
!  write( ioe, 1000 ) ( xdat%k(:,n), n = 1,  ncoord ) 
!  write( ioe, * )

!  +---------------------------------------------------------------------------+
!  |  Write Cartesian hessian                                                  |
!  +---------------------------------------------------------------------------+

   write( ioe, '( A ) ') '[cartesian hessian]'
   nn = 0
   do m = 1, natm*3
      do n = 1, m
         nn = nn + 1
         th_cart(nn) = hess_cart(m,n)
      enddo
   enddo
   write( ioe, '( A ) ') '[cartesian hessian]'
   write( ioe, 9999 ) ( th_cart(n), n = 1, natm*3*(natm*3+1)/2 )
   write( ioe, * )

!  write( ioe, '( A ) ') '[cartesian hessian]'
!  write( ioe, 1000 ) ( hess_cart(:,n), n = 1, natm*3 )
!  write( ioe, * )

!  +---------------------------------------------------------------------------+
!  |  Write atomic mass                                                        |
!  +---------------------------------------------------------------------------+

   write( ioe, '( A ) ' ) '[mass]'
   write( ioe, 1000 ) ( mass(n), n = 1, natm )
   write( ioe, * )

   close( ioe )

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

 1000 format( ( 5( F16.8 ) ) )
 9999 format( 5( 1PE16.8 ) )

   end subroutine evb_io_debug2

