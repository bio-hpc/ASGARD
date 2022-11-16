
!************************************************************************
!                              AMBER                                   **
!                                                                      **
!                  Copyright (c) 1986, 1991, 1995                      **
!             Regents of the University of California                  **
!                       All Rights Reserved.                           **
!                                                                      **
!  This software provided pursuant to a license agreement containing   **
!  restrictions on its disclosure, duplication, and use. This software **
!  contains confidential and proprietary information, and may not be   **
!  extracted or distributed, in whole or in part, for any purpose      **
!  whatsoever, without the express written permission of the authors.  **
!  This notice, and the associated author list, must be attached to    **
!  all copies, or extracts, of this software. Any additional           **
!  restrictions set forth in the license agreement also apply to this  **
!  software.                                                           **
!************************************************************************


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine difang here]
subroutine difang(cst,s,dc,ddc,imat)
   implicit double precision (a-h,o-z)
   
   !     ----- sets up second derivatives for angle-like potentials -----
   
   !           cst ... cosine of theta
   !           s   ... bond vectors ij and kj
   !           dc  ... first der. of costheta w/respct to the cartesian
   !                   differences in s
   !           ddc ... output of second derivative of costheta w/respect
   !                   to the cartesian differences
   
   !           coded by d.a.case
   
   dimension s(6),dc(6),ddc(6,6), imat(*)
   
   !      ----- first set up distances needed -----
   
   bij = s(1)*s(1)+s(2)*s(2)+s(3)*s(3)
   boi2 = 1.0e0/bij
   bij = sqrt(bij)
   bkj = s(4)*s(4)+s(5)*s(5)+s(6)*s(6)
   boj2 = 1.0e0/bkj
   bkj = sqrt(bkj)
   br = boi2
   
   !     ----- set up ddc first as second derivative of costheta w/respect
   !           to cartesian differences -----
   
   do i=1,6
      if(i == 4) br = boj2
      ddc(i,i) = -br*(2.0e0*s(i)*dc(i)+cst*(1.0e0-s(i)*s(i)*br))
   end do
   
   do i=1,3
      ddc(i,i+3) = -boi2*s(i)*dc(i+3)-boj2*s(i+3)*dc(i) &
            -boi2*boj2*(s(i)*s(i+3)*cst-bij*bkj)
   end do
   
   do i=1,2
      i1 = i+1
      do j=i1,3
         ddc(i,j) = -boi2*(s(i)*dc(j)+s(j)*dc(i)-boi2*cst*s(i)*s(j))
      end do
   end do
   
   do i=4,5
      i1 = i+1
      do j=i1,6
         ddc(i,j) = -boj2*(s(i)*dc(j)+s(j)*dc(i)-boj2*cst*s(i)*s(j))
      end do
   end do
   
   boij3 = boj2/(bij*bkj)
   boji3 = boi2/(bij*bkj)
   bt = boi2*boj2*cst
   do i=1,3
      do j=4,6
         if(j /= i+3) &
               ddc(i,j) = -boij3*s(i+3)*s(j)-boji3*s(i)*s(j-3)+bt*s(i)*s(j)
      end do
   end do
   
   !     ----- symmetrize the ddc array -----
   
   do i=1,6
      do j=i,6
         ddc(j,i) = ddc(i,j)
      end do
   end do
   return
end subroutine difang 
