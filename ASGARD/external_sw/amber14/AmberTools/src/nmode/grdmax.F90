
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
!+ [Enter a one-line description of subroutine grdmax here]
function grdmax(g,gn,n)
   implicit double precision (a-h,o-z)
   dimension g(2)
   dum = 0.0e+00
   dumg = dum
   do i = 1,n
      gi = dabs(g(i))
      if(gi > dum) dum = gi
      dumg = dumg+gi*gi
   end do
   gn = dsqrt(dumg/dble(n))
   grdmax = dum
   return
end function grdmax 
