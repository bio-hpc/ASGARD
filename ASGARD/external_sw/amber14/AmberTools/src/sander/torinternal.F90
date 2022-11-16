! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Transform from Cartesian coordinates to redundant internals            |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#include "../include/dprec.fh"

   subroutine torinternal ( qcart, q, bmat, ibond, iangle, idihed &
                          , natm, ncoord, nbond, nangle, ndihed )

   implicit none

   integer, intent(in ) :: natm, ncoord, nbond, nangle, ndihed
   integer, intent(in ) :: ibond(nbond,2), iangle(nangle,3), idihed(ndihed,4)
   _REAL_ , intent(in ) :: qcart(3,natm)
   _REAL_ , intent(out) :: q(ncoord), bmat(ncoord,natm*3)

   !  ..........................................................................

   integer tdx

!  +---------------------------------------------------------------------------+
!  |  Form B matrix                                                            |
!  +---------------------------------------------------------------------------+

   tdx = 0
   bmat(:,:) = 0.0d0 

   call bnd_strtch ( qcart, q, bmat, ibond, natm, ncoord, nbond, tdx )
   call ang_bend ( qcart, q, bmat, iangle, natm, ncoord, nangle, tdx)
   call dihed_tors ( qcart, q, bmat, idihed, natm, ncoord, ndihed, tdx )

!  +---------------------------------------------------------------------------+
!  |  Form B^-1 matrix                                                         |
!  +---------------------------------------------------------------------------+

!     nrint = natm3 - 6

!     call bmat_inv ( bmat, binv, umat, mass, natm, natm3, nintc &
!    &              , nrint, ldebug )


!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine torinternal

