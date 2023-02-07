
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
!+ [Enter a one-line description of subroutine debug here]
subroutine debug (dbl, n1, int, n2)
   double precision dbl(n1)
   integer int(n2)
   logical first
   data first /.true./
   
   if (first) then
      first = .false.
      call amopen(15, 'debug.dat', 'N', 'F', 'W')
   end if
   write (15,*) ' debug:'
   if (n1 /= 0) write (15,1) (dbl(i),i=1,n1)
   if (n2 /= 0) write (15,*) (int(i),i=1,n2)
   1 format (3(2x,e11.4))
   
end subroutine debug 
