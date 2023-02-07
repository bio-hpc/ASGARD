! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Convert Cartesians to redundant internal coordinates                   |#
! #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#+#+#+#+#

#include "../include/dprec.fh"

   subroutine cart2internal ( qcart, q )

   use schlegel, only: ncoord, natm, ibond, iangle, idihed, nbond, nangle, ndihed
   use evb_math, only: dist, angle, dihedral

   implicit none

   _REAL_ , intent( in) :: qcart(3,natm)
   _REAL_ , intent(out) :: q(ncoord)

   !  ..........................................................................

   integer :: j, k, l, m, n, nn 


   nn = 0

!  +---------------------------------------------------------------------------+
!  |  Bond stretch                                                             |
!  +---------------------------------------------------------------------------+

   do n = 1, nbond

      j = ibond(n,1)
      k = ibond(n,2)

      nn = nn + 1
      q(nn) = dist(qcart,j,k)

   enddo

!  +---------------------------------------------------------------------------+
!  |  Angle bend                                                               |
!  +---------------------------------------------------------------------------+

   do n = 1, nangle

      j = iangle(n,1)
      k = iangle(n,2)
      l = iangle(n,3)

      nn = nn + 1
      q(nn) = angle(qcart,j,l,k)

   enddo

!  +---------------------------------------------------------------------------+
!  |  Proper dihedral torsion                                                  |
!  +---------------------------------------------------------------------------+

   do n = 1, ndihed

      j = idihed(n,1)
      k = idihed(n,2)
      l = idihed(n,3)
      m = idihed(n,4)

      nn = nn + 1
      q(nn) = dihedral(qcart,j,k,l,m)

   enddo

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine cart2internal

