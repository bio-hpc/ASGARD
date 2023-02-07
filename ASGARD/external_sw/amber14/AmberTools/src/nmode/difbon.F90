
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
!+ [Enter a one-line description of subroutine difbon here]
subroutine difbon(dd,df,ddf,xij,b,i3p,j3p,nbel)
   implicit double precision (a-h,o-z)
   
   !     ----- calculate second derivatives of bond-like terms -----
   
   !           dd  ... second derivative matrix
   !           df  ... first der. of v with respect to bond length
   !           ddf ... corresponding second der.
   !           xij ... bond vector x(i)-x(j) etc
   !           b   ... bond length
   !           i3p ... atom coordinate index for atom i (3*i-3)
   !           j3p ... atom coordinate for atom j
   
   !           coded by d.a.case
   
   dimension dr(6),ddr(6,6), nbel(*), istore(6,6)
   dimension dd(2),xij(2)
   ibel(m3) = nbel(m3/3+1)
   

   i3 = ibel(min(i3p,j3p))
   j3 = ibel(max(i3p,j3p))
   
   !     ----- first set up array dr, which holds first derivative of bond
   !           length with respect to cartesian coordinates -----
   
   do k = 1,3
      dr(k) = xij(k)/b
      dr(k+3) = -dr(k)
   end do
   
   !     ----- now set up array ddr, which contains the second derivatives
   !           of bond length with respect to cartesians -----
   
   a = (1.0e0-dr(1)*dr(1))/b
   ddr(1,1) = a
   ddr(4,4) = a
   ddr(1,4) = -a
   a = (1.0e0-dr(2)*dr(2))/b
   ddr(2,2) = a
   ddr(5,5) = a
   ddr(2,5) = -a
   a = (1.0e0-dr(3)*dr(3))/b
   ddr(3,3) = a
   ddr(6,6) = a
   ddr(3,6) = -a
   a = -dr(1)*dr(2)/b
   ddr(1,2) = a
   ddr(4,5) = a
   ddr(1,5) = -a
   ddr(2,4) = -a
   a = -dr(1)*dr(3)/b
   ddr(1,3) = a
   ddr(4,6) = a
   ddr(1,6) = -a
   ddr(3,4) = -a
   a = -dr(2)*dr(3)/b
   ddr(2,3) = a
   ddr(5,6) = a
   ddr(2,6) = -a
   ddr(3,5) = -a
   
   !     ----- FORM THE INDEX ARRAY AND MERGE DF AND DDF -----
   
   istore(1,1) = i3*(i3+1)/2+i3+1
   istore(1,2) = istore(1,1)+i3+1
   istore(2,2) = istore(1,2)+1
   istore(1,3) = istore(1,2)+i3+2
   istore(2,3) = istore(1,3)+1
   istore(3,3) = istore(2,3)+1
   istore(1,4) = j3*(j3+1)/2+i3+1
   istore(2,4) = istore(1,4)+1
   istore(3,4) = istore(2,4)+1
   istore(1,5) = istore(1,4)+j3+1
   istore(2,5) = istore(1,5)+1
   istore(3,5) = istore(2,5)+1
   istore(1,6) = istore(1,5)+j3+2
   istore(2,6) = istore(1,6)+1
   istore(3,6) = istore(2,6)+1
   istore(4,4) = istore(1,4)-i3+j3
   istore(4,5) = istore(4,4)+j3+1
   istore(5,5) = istore(4,5)+1
   istore(4,6) = istore(4,5)+j3+2
   istore(5,6) = istore(4,6)+1
   istore(6,6) = istore(5,6)+1
   
   !     ----- NOW SUM UP SECOND DERIVATIVES -----
   
   is = 1
   il = 6
   if(i3 == -1) is = 4
   if(j3 == -1) il = 3
   do i = is,il
      do j = i,il
         k = istore(i,j)
         dd(k) = dd(k)+df*ddr(i,j)+ddf*dr(i)*dr(j)
      end do
   end do
   
   return
end subroutine difbon 
