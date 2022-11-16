
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
!+ [Enter a one-line description of subroutine dirch here]
subroutine dirch(vect,x,n3,idir,xinit,xdir,alph1)
   
   ! -- this routine chooses idir for downhill direction
   !    to go to a NEW minima.
   
   implicit double precision (a-h,o-z)
   dimension vect(n3,1),x(n3),xinit(n3),xdir(n3)
   do 10 i=1,n3
      xinit(i) = x(i)-xinit(i)
   10 continue
   proj = ddot(n3,xinit,1,vect,1)
   if(proj <= 0.00d0) then
      idir = -1
   else
      idir = 1
   end if
#ifdef debug
   write(6,101) idir
   101 format('Choosing idir =',i2,' for downhill walk')
#endif
   
   do 20 i=1,n3
      xdir(i) = idir*alph1*vect(i,1)
   20 continue
   return
end subroutine dirch 
