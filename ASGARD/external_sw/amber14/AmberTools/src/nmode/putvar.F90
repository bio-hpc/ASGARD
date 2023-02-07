
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
!+ [Enter a one-line description of subroutine putvar here]
subroutine putvar(nat,igrp,xbel,x)
   real*8 xbel,x
   
   !     ----- ROUTINE TO STUFF THE VARYING (belly) COORDINATES INTO THE
   !           FULL COORDINATE ARRAY -----
   
   dimension igrp(*),xbel(*),x(*)
   
   j3 = 0
   do 10 i = 1,nat
      if(igrp(i) > 0) then
         i3 = 3*i-3
         x(i3+1) = xbel(j3+1)
         x(i3+2) = xbel(j3+2)
         x(i3+3) = xbel(j3+3)
         j3 = j3+3
      end if
   10 continue
   return
end subroutine putvar 
