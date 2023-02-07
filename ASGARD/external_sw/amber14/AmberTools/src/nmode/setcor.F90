
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
!+ [Enter a one-line description of subroutine setcor here]
subroutine setcor(nat,igrp,xp,x)
   double precision xp,x
   dimension xp(1),x(1),igrp(1)
   j3 = 0
   do i = 1,nat
      if(igrp(i) > 0) then
         i3 = 3*i-3
         x(j3+1) = xp(i3+1)
         x(j3+2) = xp(i3+2)
         x(j3+3) = xp(i3+3)
         j3 = j3+3
      end if
   end do
   return
end subroutine setcor 
