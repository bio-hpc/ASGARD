
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
!+ [Enter a one-line description of subroutine cdiv here]
subroutine cdiv(ar,ai,br,bi,cr,ci)
   double precision ar,ai,br,bi,cr,ci
   
   !     complex division, (cr,ci) = (ar,ai)/(br,bi)
   
   double precision s,ars,ais,brs,bis
   s = dabs(br) + dabs(bi)
   ars = ar/s
   ais = ai/s
   brs = br/s
   bis = bi/s
   s = brs**2 + bis**2
   cr = (ars*brs + ais*bis)/s
   ci = (ais*brs - ars*bis)/s
   return
end subroutine cdiv 
