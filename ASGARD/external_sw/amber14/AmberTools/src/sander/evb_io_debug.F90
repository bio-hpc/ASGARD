! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Output data read from the Gaussian formatted checkpoint file           |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#include "../include/assert.fh"
#include "../include/dprec.fh"

   subroutine evb_io_debug ( qcart, grad_cart, hess_cart, mass, bmat, binv &
                           , dbmat, xdat, ioe )

   use schlegel, only: berny_type, ncoord, natm &
                     , ibond, iangle, idihed, nbond, nangle, ndihed

   implicit none

#  include "parallel.h"

   type( berny_type ) :: xdat

   integer, intent(in) :: ioe

   _REAL_, intent(in) :: qcart(natm*3), grad_cart(natm*3) &
                       , hess_cart(natm*3,natm*3), mass(natm) &
                       , bmat(ncoord,natm*3), binv(natm*3,ncoord) &
                       , dbmat(ncoord,natm*3*natm*3)

   !  ..........................................................................

   integer ::  m, n, nn 
   _REAL_  :: th(ncoord*(ncoord+1)/2), th_cart(natm*3*(natm*3+1)/2)

   open( ioe, file = trim( adjustl( xdat%filename ) )//'.debug' )

!  +---------------------------------------------------------------------------+
!  |  Write coordinate type                                                    |
!  +---------------------------------------------------------------------------+
   write( ioe, '( A )' ) '[coordinate type]'
   write( ioe, '( A )' ) '  use_internals'
   write( ioe, * )

!  +---------------------------------------------------------------------------+
!  |  Write external EVB redundant internal data dimension                     |
!  +---------------------------------------------------------------------------+
   write( ioe, '( A )' ) '[external evb data dimension]'
   write( ioe, '(5I12)' ) ncoord, natm, nbond, nangle, ndihed 
   write( ioe, * )

!  +---------------------------------------------------------------------------+
!  |  Write internal coordinate specifications                                 |
!  +---------------------------------------------------------------------------+

   write( ioe, '( A )' ) '[redundant internal indices]'
   do n = 1, nbond
      write( ioe, '(4I12)' ) ibond(n,1), ibond(n,2), 0, 0
   enddo
   do n = 1, nangle
      write( ioe, '(4I12)' ) iangle(n,1), iangle(n,3), iangle(n,2), 0
   enddo
   do n = 1, ndihed
      write( ioe, '(4I12)' ) idihed(n,1), idihed(n,2), idihed(n,3), idihed(n,4)
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
   write( ioe, 1000 ) ( qcart(n), n = 1, natm * 3 )
   write( ioe, * )

!  +---------------------------------------------------------------------------+
!  |  Write energy                                                             |
!  +---------------------------------------------------------------------------+
   write( ioe, '( A )' ) '[electronic energy]'
   write( ioe, '( X, 1PE22.15 )' ) xdat%v
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
!  |  Write internal hessian (lower triangle)                                  |
!  +---------------------------------------------------------------------------+
   nn = 0
   do m = 1, ncoord
      do n = 1, m
         nn = nn + 1
         th(nn) = xdat%k(m,n)
      enddo
   enddo
   write( ioe, '( A ) ') '[redundant internal hessian]'
   write( ioe, 1000 ) ( th(n), n = 1,  ncoord*(ncoord+1)/2 ) 
   write( ioe, * )

!  +---------------------------------------------------------------------------+
!  |  Write Cartesian hessian (lower triangle)                                 |
!  +---------------------------------------------------------------------------+
   nn = 0
   do m = 1, natm*3 
      do n = 1, m
         nn = nn + 1
         th_cart(nn) = hess_cart(m,n)
      enddo
   enddo
   write( ioe, '( A ) ') '[cartesian hessian]'
   write( ioe, 1000 ) ( th_cart(n), n = 1, natm*3*(natm*3+1)/2 )
   write( ioe, * )

!  +---------------------------------------------------------------------------+
!  |  Write atomic mass                                                        |
!  +---------------------------------------------------------------------------+
   write( ioe, '( A ) ' ) '[mass]'
   write( ioe, 1000 ) ( mass(n), n = 1, natm )
   write( ioe, * )

  goto 888 

!  +---------------------------------------------------------------------------+
!  |  Write B matrix                                                           |
!  +---------------------------------------------------------------------------+
   write( ioe, '( A ) ') '[B matrix]'
   write( ioe, 1000 ) ( bmat(:,n), n = 1, natm*3 )
   write( ioe, * )

!  +---------------------------------------------------------------------------+
!  |  Write B^-1 matrix                                                        |
!  +---------------------------------------------------------------------------+
   write( ioe, '( A ) ') '[B^-1 matrix]'
   write( ioe, 1000 ) ( binv(:,n), n = 1, ncoord )
   write( ioe, * )

!  +---------------------------------------------------------------------------+
!  |  Write gradient of B matrix                                               |
!  +---------------------------------------------------------------------------+
   write( ioe, '( A ) ') '[gradient of B matrix]'
   write( ioe, 1000 ) ( dbmat(:,n), n = 1, natm*3*natm*3 )
   write( ioe, * )

 888 continue

   close( ioe )

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

 1000 format( ( 5( 1PE16.8 ) ) )

   end subroutine evb_io_debug

