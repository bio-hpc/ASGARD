
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
!+ [Enter a one-line description of subroutine remrot here]
subroutine remrot (natom, nreal, d, x, z, wr, wi, ind, dnorm, &
      amass, dbl, iclass)
   
   !     ----- remove translational and rotational parts from real
   !     ----- Langevin modes
   
   implicit double precision (a-h,o-z)
   dimension d(3*natom,6), x(3*natom), z(6*natom,6*natom), a(6), &
         ind(6*natom), amass(natom), dbl(6*natom), &
         wr(6*natom), wi(6*natom), dnorm(6*natom), &
         iclass(6*natom)
   
   xx(i) = x(3*(i-1)+1)
   yy(i) = x(3*(i-1)+2)
   zz(i) = x(3*(i-1)+3)
   
   n = 3 * natom
   
   !     ---- normalize transaltion and rotation degrees of freedom
   
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
   if (vnor1 == 0.0) vnor1 = 1.0
   if (vnor2 == 0.0) vnor2 = 1.0
   if (vnor3 == 0.0) vnor3 = 1.0
   
   ! --- set up the trsl. and rot.  vectors
   !     see C. Eckart , Phys. Rev. 47, 552 (1935)
   
   do 20 i=1,n
      do 15 j=1,6
         d(i,j) = 0.0
      15 continue
   20 continue
   iat = 1
   do 30 i=1,n,3
      d(i,1) = dsqrt(amass(iat))/vnor4
      d(i+1,2) = dsqrt(amass(iat))/vnor4
      d(i+2,3) = dsqrt(amass(iat))/vnor4
      d(i,4) = -dsqrt(amass(iat))*yy(iat)/vnor1
      d(i+1,4) = dsqrt(amass(iat))*xx(iat)/vnor1
      d(i+1,5) = -dsqrt(amass(iat))*zz(iat)/vnor2
      d(i+2,5) = dsqrt(amass(iat))*yy(iat)/vnor2
      d(i,6) = dsqrt(amass(iat))*zz(iat)/vnor3
      d(i+2,6) = -dsqrt(amass(iat))*xx(iat)/vnor3
      iat = iat + 1
   30 continue
   
   !     ----- subtract translational and rotational contributions
   !     ----- from each eigenvector
   
   do 90 k = 7, nreal
      
      ik = ind(k)
      
      do 50 j = 1, 6
         a(j) = 0.0
         do 40 i = 1, n
            a(j) = a(j) + z(i,ik)*d(i,j)
         40 continue
      50 continue
      
      do 70 i = 1, n
         sum = 0.0
         do 60 j = 1, 6
            sum = sum + a(j) * d(i,j)
         60 continue
         diff = z(i,ik) - sum
         dnorm(ik) = dnorm(ik) + diff*diff
      70 continue
      sum = dnorm(ik)
      do 80 j = 1, 6
         sum = sum + a(j)*a(j)
      80 continue
      if (sum /= 0.0) dnorm(ik) = dnorm(ik) / sum
      
   90 continue
   
   !     ----- see if anything remains after that subtraction
   !     ----- of course a tiny bit may remain always
   !     ----- so remove 12 modes with smallest remainders
   
   do 100 i = 7, nreal
      dbl(i) = dnorm(ind(i))
   100 continue
   
   do 120 i = 7, nreal-1
      ip = i
      do 110 j = i+1, nreal
         if (dbl(j) < dbl(ip)) ip = j
      110 continue
      dt = dbl(i)
      dbl(i) = dbl(ip)
      dbl(ip) = dt
      it = ind(i)
      ind(i) = ind(ip)
      ind(ip) = it
      it = iclass(i)
      iclass(i) = iclass(ip)
      iclass(ip) = it
   120 continue
   
   return
end subroutine remrot 
