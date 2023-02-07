! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Compute the B matrix                                                   |#
! #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#+#+#+#+#

#include "../include/dprec.fh"

   subroutine wdc_bmat ( qcart, bmat )

   use schlegel, only: ncoord, natm, ibond, iangle, idihed, nbond, nangle &
                     , ndihed
   use evb_math, only: uvec, orthovec, vrot, dist

   implicit none

   _REAL_ , intent( in) :: qcart(3,natm)
   _REAL_ , intent(out) :: bmat(ncoord,natm*3)

   !  ..........................................................................

   integer :: i, j, k, l, m, n, jj, kk, ll, mm, nn
   _REAL_  :: s_tj(3), s_tk(3), s_tl(3), s_tm(3)

   bmat(:,:) = 0.0d0
   nn = 0

!  +---------------------------------------------------------------------------+
!  |  Compute the B matrix elements for a bond stretch                         |
!  :...........................................................................:
!  |  Wilson, Decius, and Cross, "Molecular Vibrations" p. 55-56               |
!  |                                                                           |
!  |     j --> k                        r_jk = r_k - r_j                       |
!  |                                                                           |
!  |  s_tj = e_kj = - e_jk           s_tk = e_jk = - e_kj                      |
!  +---------------------------------------------------------------------------+

   do n = 1, nbond

      j = ibond(n,1)
      k = ibond(n,2)

      s_tj(:) = uvec(qcart,k,j)

      jj = ( j - 1 ) * 3
      kk = ( k - 1 ) * 3
      nn = nn + 1

      do i = 1, 3
         bmat(nn,jj+i) =  - s_tj(i)
         bmat(nn,kk+i) =    s_tj(i)
      enddo

   enddo

!  +---------------------------------------------------------------------------+
!  |  Compute the B matrix elements for an angle bend                          |
!  :...........................................................................:
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

   do n = 1, nangle

      j = iangle(n,1)
      k = iangle(n,2)
      l = iangle(n,3)

      jj = ( j - 1 ) * 3
      kk = ( k - 1 ) * 3
      ll = ( l - 1 ) * 3
      nn = nn + 1

      s_tj(:) = - orthovec(qcart,k,l,j) / dist(qcart,j,l)
      s_tk(:) = - orthovec(qcart,j,l,k) / dist(qcart,k,l)
      s_tl(:) = - s_tj(:) - s_tk(:)

      do i = 1, 3 
         bmat(nn,jj+i) = s_tj(i)
         bmat(nn,kk+i) = s_tk(i)
         bmat(nn,ll+i) = s_tl(i)
      enddo

   enddo

!  +---------------------------------------------------------------------------+
!  |  Compute the B matrix elements for a proper dihedral                      |
!  :...........................................................................:
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

   do n = 1, ndihed

      j = idihed(n,1)
      k = idihed(n,2)
      l = idihed(n,3)
      m = idihed(n,4)

      jj = ( j - 1 ) * 3
      kk = ( k - 1 ) * 3
      ll = ( l - 1 ) * 3
      mm = ( m - 1 ) * 3
      nn = nn + 1

      s_tj(:) = vrot(qcart,j,k,l)

      s_tk(:) = - vrot(qcart,j,k,l) * ( 1.0d0 - dot_product( uvec(qcart,j,k) &
              , uvec(qcart,l,k) ) * dist(qcart,j,k) / dist(qcart,k,l) ) &
                - vrot(qcart,m,l,k) * dot_product( uvec(qcart,k,l) &
              , uvec(qcart,m,l) ) * dist(qcart,l,m) / dist(qcart,k,l) 

      s_tl(:) = - vrot(qcart,m,l,k) * ( 1.0d0 - dot_product( uvec(qcart,k,l) &
              , uvec(qcart,m,l) ) * dist(qcart,l,m) / dist(qcart,k,l) ) &
                - vrot(qcart,j,k,l) * dot_product( uvec(qcart,j,k) &
              , uvec(qcart,l,k) ) * dist(qcart,j,k) / dist(qcart,k,l)

      s_tm(:) = vrot(qcart,m,l,k)

      do i = 1, 3
         bmat(nn,jj+i) = s_tj(i)
         bmat(nn,kk+i) = s_tk(i)
         bmat(nn,ll+i) = s_tl(i)
         bmat(nn,mm+i) = s_tm(i)
      enddo

   enddo

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine wdc_bmat

