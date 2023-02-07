
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
!+ [Enter a one-line description of subroutine level here]
subroutine level(h,x,amass,n,natom)
   
   !------adjust Hessian matrix to put translations and rotations at
   !       a relatively high frequency
   
   implicit double precision (a-h,o-z)
   dimension h(*),x(*),amass(*),d(6,9000)
   
   xx(i)=x(3*(i-1)+1)
   yy(i)=x(3*(i-1)+2)
   zz(i)=x(3*(i-1)+3)
   
   !---normalize transaltion and rotation degrees of freedom
   
   xsum = 0.0d0
   ysum = 0.0d0
   zsum = 0.0d0
   tmas = 0.0d0
   do 10 i = 1,natom
      xsum = xsum + xx(i)*xx(i)*amass(i)
      ysum = ysum + yy(i)*yy(i)*amass(i)
      zsum = zsum + zz(i)*zz(i)*amass(i)
      tmas = tmas + amass(i)
   10 continue
   vnor1 = dsqrt(xsum+ysum)
   vnor2 = dsqrt(ysum+zsum)
   vnor3 = dsqrt(xsum+zsum)
   vnor4 = dsqrt(tmas)
   
   ! --- set up the trsl. and rot.  vectors
   !     see C. Eckart , Phys. Rev. 47, 552 (1935)
   
   do 20 i=1,n
      do 15 j=1,6
         d(j,i) = 0.0
      15 continue
   20 continue
   iat = 1
   do 30 i=1,n,3
      d(1,i) = dsqrt(amass(iat))/vnor4
      d(2,i+1) = dsqrt(amass(iat))/vnor4
      d(3,i+2) = dsqrt(amass(iat))/vnor4
      d(4,i) = -dsqrt(amass(iat))*yy(iat)/vnor1
      d(4,i+1) = dsqrt(amass(iat))*xx(iat)/vnor1
      d(5,i+1) = -dsqrt(amass(iat))*zz(iat)/vnor2
      d(5,i+2) = dsqrt(amass(iat))*yy(iat)/vnor2
      d(6,i) = dsqrt(amass(iat))*zz(iat)/vnor3
      d(6,i+2) = -dsqrt(amass(iat))*xx(iat)/vnor3
      iat = iat + 1
   30 continue
   
   ! --- here do level shift
   
   k = 0
   do 40 i=1,n
      do 35 j=1,i
         k = k + 1
         do 33 ic=1,6
            h(k) = h(k) + 3500.00d0*d(ic,i)*d(ic,j)
         33 continue
      35 continue
   40 continue
   return
end subroutine level 
