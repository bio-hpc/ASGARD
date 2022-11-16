! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Angle bend elements of the Wilson B matrix                             |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#include "../include/dprec.fh"

   subroutine ang_bend ( qcart, q, bmat, iangle, natm, ncoord, nangle, tdx )

   implicit none

   integer, intent(in   ) :: natm, ncoord, nangle
   integer, intent(inout) :: tdx
   integer, intent(in   ) :: iangle(nangle,3)
   _REAL_ , intent(in   ) :: qcart(3,natm)
   _REAL_ , intent(inout) :: q(ncoord), bmat(ncoord,natm*3)

   !  ..........................................................................

   integer :: i, t, jb, kb, lb, jndx, kndx, lndx
   _REAL_  :: e_lj(3), e_lk(3), rlj, rlk, rlj_inv, rlk_inv &
            , cos_phi, sin_phi_inv, s_tj, s_tk, s_tl  
   _REAL_  , parameter :: eps = 1.0d-8
   intrinsic :: dot_product, max, sqrt, acos

!  +---------------------------------------------------------------------------+
!  |  Compute the B matrix elements for an angle bend                          |
!  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  |  Wilson, Decius, and Cross, "Molecular Vibrations" p. 56-58               |
!  |                                                                           |
!  |       l              s_tj = ( cos(phi) e_lj - e_lk )                      | 
!  |     /   \                 / ( r_lj sin(phi) )                             |
!  |   j       k                                                               | 
!  |                      s_tk = ( cos(phi) e_lk - e_lj )                      | 
!  |                           / ( r_lk sin(phi) )                             | 
!  |                                                                           |
!  |                      s_tl = - s_tj - s_tk                                 |
!  |                                                                           | 
!  |   cos(phi) = e_lj * e_lk   sin(phi) = SQRT( 1 - cos^2(phi) )              |
!  +---------------------------------------------------------------------------+

!  +---------------------------------------------------------------------------+
!  |  Obtain unit vectors and bend angle                                       |
!  +---------------------------------------------------------------------------+

   do t = 1, nangle

      jb = iangle(t,1) 
      kb = iangle(t,2) 
      lb = iangle(t,3)

      call unitv ( qcart(1,lb), qcart(1,jb), e_lj, rlj, rlj_inv )
      call unitv ( qcart(1,lb), qcart(1,kb), e_lk, rlk, rlk_inv )

      cos_phi = dot_product( e_lj, e_lk )
      sin_phi_inv = 1.0d0 / max( sqrt( 1.0d0 - cos_phi**2 ), eps ) 

!  +---------------------------------------------------------------------------+
!  |  Compute elements and populate B matrix                                   |
!  +---------------------------------------------------------------------------+

      tdx = tdx + 1

      q(tdx) = acos( cos_phi ) 

      jndx = ( jb - 1 ) * 3
      kndx = ( kb - 1 ) * 3
      lndx = ( lb - 1 ) * 3

      do i = 1, 3 

         s_tj = ( cos_phi * e_lj(i) - e_lk(i) ) * rlj_inv * sin_phi_inv
         s_tk = ( cos_phi * e_lk(i) - e_lj(i) ) * rlk_inv * sin_phi_inv
         s_tl = - s_tj - s_tk 

         bmat(tdx,jndx+i) = s_tj
         bmat(tdx,kndx+i) = s_tk
         bmat(tdx,lndx+i) = s_tl

      enddo          

   enddo 

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine ang_bend


