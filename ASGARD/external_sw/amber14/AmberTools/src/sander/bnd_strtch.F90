! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Bond stretch elements of the Wilson B matrix                           |#
! #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#+#+#+#+#

#include "../include/dprec.fh"

   subroutine bnd_strtch ( qcart, q, bmat, ibond, natm, ncoord, nbond, tdx )

   use evb_math, only: dist, uvec

   implicit none 

   integer, intent(in   ) :: natm, ncoord, nbond
   integer, intent(inout) :: tdx
   integer, intent(in   ) :: ibond(nbond,2)
   _REAL_ , intent(in   ) :: qcart(3,natm)
   _REAL_ , intent(inout) :: q(ncoord), bmat(ncoord,natm*3)

   !  ..........................................................................

   integer :: i, t, jb, kb, jndx, kndx
   _REAL_  :: e_jk(3)

!  +---------------------------------------------------------------------------+
!  |  Compute the B matrix elements for a bond stretch                         |
!  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  |  Wilson, Decius, and Cross, "Molecular Vibrations" p. 55-56               |
!  |                                                                           |
!  |     j --> k                        r_jk = r_k - r_j                       |
!  |                                                                           |
!  |  s_tj = e_kj = - e_jk           s_tk = e_jk = - e_kj                      |
!  +---------------------------------------------------------------------------+

   do t = 1, nbond

      tdx = tdx + 1

      jb = ibond(t,1) 
      kb = ibond(t,2) 

      e_jk(:) = uvec( qcart(:,:), kb, jb )
      q(tdx)  = dist( qcart(:,:), kb, jb )

      jndx = ( jb - 1 ) * 3
      kndx = ( kb - 1 ) * 3

      do i = 1, 3 
         bmat(tdx,jndx+i) =  - e_jk(i) 
         bmat(tdx,kndx+i) =    e_jk(i) 
      enddo 

   enddo 

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine bnd_strtch

