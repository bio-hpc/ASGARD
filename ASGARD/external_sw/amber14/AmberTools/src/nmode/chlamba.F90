
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
!+ [Enter a one-line description of subroutine chlamba here]
subroutine chlamba(tlamba,vect,roots,scale,fg, &
      xdir,n3,b1,b2)
   
   !----this subroutine calculates lambda and does coordinates transformation
   !          if necessary, returns correction to xdir later needed for
   !          scaling purposes
   
   implicit double precision (a-h,o-z)
#  include "inpdat.h"
   dimension roots(*),vect(n3,*),fg(*),xdir(*)
   do 15 i =1,n3
      xdir(i) = 0.0d0
   15 continue
   
   tlamba =  b1 + 0.50d0*(0.5*b2-b1)
   if(b1 <= buphl) tlamba = 0.00d0
   write(6,10) tlamba
   10 format('   lambda = ',f10.4)
   
   !----here do coordinate transformation
   
   if(scale == 1.00d0) return
   jdot =isdir
   xsum = ddot(n3,fg,1,vect,jdot)
   dem1 = tlamba - scale*scale*roots(isdir)
   dem2 = tlamba - roots(isdir)
   beta = (scale*scale*xsum/dem1) - (xsum/dem2)
   do 20 ians = 1,n3
      xdir(ians) = beta*vect(ians,isdir)
   20 continue
   return
end subroutine chlamba 
