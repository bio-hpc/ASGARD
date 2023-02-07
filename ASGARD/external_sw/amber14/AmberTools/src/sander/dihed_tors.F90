! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Proper dihedral elements of the Wilson B matrix                          |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#include "../include/dprec.fh"

   subroutine dihed_tors ( qcart, q, bmat, idihed, natm, ncoord, ndihed, tdx )

   implicit none

   integer, intent(in   ) :: natm, ncoord, ndihed
   integer, intent(inout) :: tdx
   integer, intent(in   ) :: idihed(ndihed,4)
   _REAL_ , intent(in   ) :: qcart(3,natm)
   _REAL_ , intent(inout) :: q(ncoord), bmat(ncoord,natm*3) 

   !  ..........................................................................

   integer :: i, t, jb, kb, lb, mb, jndx, kndx, lndx, mndx

   _REAL_ :: e_jk(3), e_kl(3), e_lm(3), e_kj(3), e_lk(3), e_ml(3)
   _REAL_ :: ejk_x_ekl(3), eml_x_elk(3), ekl_x_elm(3)
   _REAL_ :: rjk, rkl, rlm, rjk_inv, rkl_inv, rlm_inv &
           , cos_phi2, cos_phi3, sinsq_phi2_inv, sinsq_phi3_inv &
           , s_tj, s_tk, s_tl, s_tm, cos_tau, pi
   intrinsic :: dot_product, sqrt, acos, abs

!  +---------------------------------------------------------------------------+
!  |  Compute the B matrix elements for a proper dihedral                      |
!  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  |  Wilson, Decius, and Cross, "Molecular Vibrations" p. 60-61               |
!  |                                                                           |
!  |                 s_tj = - e_jk X e_kl / r_jk / sin^2(phi2)                 |
!  |  j  (phi2)                                                                |
!  |    \            s_tk = ( r_kl - r_jk cos(phi2) )                          |
!  |      k -- l          * e_jk X e_kl / r_kl / r_jk                          |
!  |             \        / sin^2(phi2) + cos(phi3) e_ml x e_lk                |
!  |       (phi3)  m      / r_kl / sin^2(phi3)                                 |
!  |                                                                           |
!  |                 s_tl = [(jm)(kl)]_s_tk                                    | 
!  |                                                                           | 
!  |                 s_tm = [(jm)(kl)]_s_tj                                    |
!  |                                                                           |
!  |     cos(tau) = (e_jk X e_kl) * (e_kl X e_lm)                              |
!  |              / sin(phi2) / sin(phi3)                                      |
!  |                                                                           |
!  |   s_tl & s_tm are obtained via permutation of j with m and                |
!  |   k with l and correponding phi2 and phi3 angles                          |
!  +---------------------------------------------------------------------------+

   pi = acos( -1.0d0 ) 

!  +---------------------------------------------------------------------------+
!  |  Obtain unit vectors, cross products, and angles                          |
!  +---------------------------------------------------------------------------+

   do t = 1, ndihed

      jb = idihed(t,1) 
      kb = idihed(t,2) 
      lb = idihed(t,3)
      mb = idihed(t,4)

      call unitv ( qcart(1,jb), qcart(1,kb), e_jk, rjk, rjk_inv ) 
      call unitv ( qcart(1,kb), qcart(1,lb), e_kl, rkl, rkl_inv )   
      call unitv ( qcart(1,lb), qcart(1,mb), e_lm, rlm, rlm_inv )   

      do i = 1, 3 
         e_kj(i) = - e_jk(i) 
         e_lk(i) = - e_kl(i) 
         e_ml(i) = - e_lm(i) 
      enddo 

      call crossv ( e_jk, e_kl, ejk_x_ekl )
      call crossv ( e_ml, e_lk, eml_x_elk )
      ekl_x_elm = - eml_x_elk

      cos_phi2 = dot_product( e_kj, e_kl )
      cos_phi3 = dot_product( e_lk, e_lm )

      sinsq_phi2_inv = 1.0d0 / ( 1.0d0 - cos_phi2**2 )
      sinsq_phi3_inv = 1.0d0 / ( 1.0d0 - cos_phi3**2 )

!  +---------------------------------------------------------------------------+
!  |  Compute elements and populate B matrix                                   |
!  +---------------------------------------------------------------------------+

      tdx = tdx + 1 

      cos_tau = dot_product( ejk_x_ekl, ekl_x_elm ) * sqrt( sinsq_phi2_inv ) &
                                                    * sqrt( sinsq_phi3_inv )
  
      if( abs(cos_tau) > 1.0d0 ) then
         if( cos_tau > 0.0d0 ) then
            q(tdx) = 0.0d0
         else
            q(tdx) = pi
         endif 
      else
         q(tdx) = acos( cos_tau )
      endif  

      jndx = ( jb - 1 ) * 3
      kndx = ( kb - 1 ) * 3
      lndx = ( lb - 1 ) * 3
      mndx = ( mb - 1 ) * 3

      do i = 1, 3 

         s_tj = - ejk_x_ekl(i) * rjk_inv * sinsq_phi2_inv 

         s_tk = ( rkl - rjk * cos_phi2 ) *  ejk_x_ekl(i) * rkl_inv &
              * rjk_inv * sinsq_phi2_inv &
              + cos_phi3 * eml_x_elk(i) * rkl_inv * sinsq_phi3_inv

         s_tl = ( rkl - rlm * cos_phi3 ) * eml_x_elk(i) * rkl_inv &
              * rlm_inv * sinsq_phi3_inv &
              + cos_phi2 * ejk_x_ekl(i) * rkl_inv * sinsq_phi2_inv

         s_tm = - eml_x_elk(i) * rlm_inv * sinsq_phi3_inv 

         bmat(tdx,jndx+i) = s_tj
         bmat(tdx,kndx+i) = s_tk
         bmat(tdx,lndx+i) = s_tl
         bmat(tdx,mndx+i) = s_tm

      enddo 

   enddo 

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine dihed_tors

