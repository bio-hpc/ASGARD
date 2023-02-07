
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

double precision function ddot(n,dx,incx,dy,incy)

!     ----- RETURNS THE DOT PRODUCT OF DOUBLE PRECISION DX AND DY -----

double precision dx(1),dy(1)

ddot = 0.d0
if(n <= 0)return
if(incx == incy) if(incx-1) 5,20,60
5 continue

!     ----- CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS -----

ix = 1
iy = 1
if(incx < 0)ix = (-n+1)*incx + 1
if(incy < 0)iy = (-n+1)*incy + 1
do 10 i = 1,n
   ddot = ddot + dx(ix)*dy(iy)
   ix = ix + incx
   iy = iy + incy
10 continue
return

!     ----- CODE FOR BOTH INCREMENTS EQUAL TO 1 -----

!     ----- CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A
!           MULTIPLE OF 5 -----

20 m = n - (n/5)*5
if( m == 0 ) goto 40
do 30 i = 1,m
   ddot = ddot + dx(i)*dy(i)
30 continue
if( n < 5 ) return
40 mp1 = m + 1
do 50 i = mp1,n,5
   ddot = ddot + dx(i)*dy(i) + dx(i+1)*dy(i+1) + &
         dx(i+2)*dy(i+2) + dx(i+3)*dy(i+3) + dx(i+4)*dy(i+4)
50 continue
return

!     ----- CODE FOR POSITIVE EQUAL INCREMENTS .NE.1 -----

60 continue
ns = n*incx
do 70 i=1,ns,incx
   ddot = ddot + dx(i)*dy(i)
70 continue
return
end subroutine daxpy !FIXME:               ddot.f, line   70: Indentation error: end subroutine daxpy 
