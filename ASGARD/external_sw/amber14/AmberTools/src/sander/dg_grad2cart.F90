! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Convert distributed gaussians gradient from redundant internal |#
! #|  coordinate to Cartesian                                        |#
! #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

#include "../include/assert.fh"
#include "../include/dprec.fh"

   subroutine dg_grad2cart ( dV_EVB, bmat, dV_Cart, natm, ncoord, ibond &
                           , iangle, idihed, nbond, nangle, ndihed )

   implicit none 

   integer, intent(in ) :: natm, ncoord, nbond, nangle, ndihed
   integer, intent(in ) :: ibond(nbond,2), iangle(nangle,3), idihed(ndihed,4)
   _REAL_ , intent(in ) :: dV_EVB(ncoord), bmat(ncoord,natm*3)
   _REAL_ , intent(out) :: dV_Cart(3,natm)

   !  +---------------------------------------------------------------+

   integer :: i, t, jb, kb, lb, mb, jndx, kndx, lndx, mndx, tdx

!  +---------------------------------------------------------------+
!  |  Bond contribution                                            |
!  +---------------------------------------------------------------+

   dV_Cart(:,:) = 0.0d0

   tdx = 0

   do t = 1, nbond

      tdx = tdx + 1

      jb = ibond(t,1) 
      kb = ibond(t,2) 

      jndx = ( jb - 1 ) * 3
      kndx = ( kb - 1 ) * 3

      do i = 1, 3 
         dV_Cart(i,jb) = dV_Cart(i,jb) + dV_EVB(tdx) * bmat(tdx,jndx+i)
         dV_Cart(i,kb) = dV_Cart(i,kb) + dV_EVB(tdx) * bmat(tdx,kndx+i)
      enddo 

   enddo 

!  +---------------------------------------------------------------+
!  |  Angle contribution                                           |
!  +---------------------------------------------------------------+

   do t = 1, nangle

      tdx = tdx + 1

      jb = iangle(t,1)
      kb = iangle(t,2)
      lb = iangle(t,3)

      jndx = ( jb - 1 ) * 3
      kndx = ( kb - 1 ) * 3
      lndx = ( lb - 1 ) * 3

      do i = 1, 3

         dV_Cart(i,jb) = dV_Cart(i,jb) + dV_EVB(tdx) * bmat(tdx,jndx+i)
         dV_Cart(i,kb) = dV_Cart(i,kb) + dV_EVB(tdx) * bmat(tdx,kndx+i)
         dV_Cart(i,lb) = dV_Cart(i,lb) + dV_EVB(tdx) * bmat(tdx,lndx+i)

      enddo

   enddo

!  +---------------------------------------------------------------+
!  |  Dihedral contribution                                        |
!  +---------------------------------------------------------------+

   do t = 1, ndihed

      tdx = tdx + 1

      jb = idihed(t,1)
      kb = idihed(t,2)
      lb = idihed(t,3)
      mb = idihed(t,4)

      jndx = ( jb - 1 ) * 3
      kndx = ( kb - 1 ) * 3
      lndx = ( lb - 1 ) * 3
      mndx = ( mb - 1 ) * 3

      do i = 1, 3

         dV_Cart(i,jb) = dV_Cart(i,jb) + dV_EVB(tdx) * bmat(tdx,jndx+i)
         dV_Cart(i,kb) = dV_Cart(i,kb) + dV_EVB(tdx) * bmat(tdx,kndx+i)
         dV_Cart(i,lb) = dV_Cart(i,lb) + dV_EVB(tdx) * bmat(tdx,lndx+i)
         dV_Cart(i,mb) = dV_Cart(i,mb) + dV_EVB(tdx) * bmat(tdx,mndx+i)

      enddo

   enddo

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine dg_grad2cart

