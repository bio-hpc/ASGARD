
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
!+ [Enter a one-line description of subroutine daxpy here]
subroutine daxpy(n,da,dx,incx,dy,incy)
   
   !     ----- overwrite double precision dy with double precision
   !           da*dx + dy -----
   
   double precision dx(1),dy(1),da
   if(n <= 0.or.da == 0.d0) return
   if(incx == incy) if(incx-1) 5,20,60
   5 continue
   
   !     ----- code for nonequal or nonpositive increments -----
   
   ix = 1
   iy = 1
   if(incx < 0)ix = (-n+1)*incx + 1
   if(incy < 0)iy = (-n+1)*incy + 1
   do 10 i = 1,n
      dy(iy) = dy(iy) + da*dx(ix)
      ix = ix + incx
      iy = iy + incy
   10 continue
   return
   
   !     ----- code for both increments equal to 1 -----
   
   !     ----- clean-up loop so remaining vector length is a
   !           multiple of 4 -----
   
   20 m = n - (n/4)*4
   if( m == 0 ) goto 40
   do 30 i = 1,m
      dy(i) = dy(i) + da*dx(i)
   30 continue
   if( n < 4 ) return
   40 mp1 = m + 1
   do 50 i = mp1,n,4
      dy(i) = dy(i) + da*dx(i)
      dy(i + 1) = dy(i + 1) + da*dx(i + 1)
      dy(i + 2) = dy(i + 2) + da*dx(i + 2)
      dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
   return
   
   !     ----- code for equal, positive, nonunit increments -----
   
   60 continue
   ns = n*incx
   do 70 i=1,ns,incx
      dy(i) = da*dx(i) + dy(i)
   70 continue
   return
end subroutine daxpy 
