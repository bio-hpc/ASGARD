
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
!+ [Enter a one-line description of subroutine zxcgr here]
subroutine zxcgr()
   
   write(6,*) 'Calling subroutine zxcgr (ntrun=3)....'
   write(6,*) '  You need to link the IMSL routine in place'
   write(6,*) '  of the dummy zxcgr.f routine supplied in '
   write(6,*) '  AMBER distribution.  Or, use the Newton-Raphson'
   write(6,*) '  minimizer (ntrun=4) or the sander program'
   write(6,*) '  to carry out your minimizations.'
   call mexit(6, 1)
end subroutine zxcgr 
