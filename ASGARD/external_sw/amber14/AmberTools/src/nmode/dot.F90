
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
!+ [Enter a one-line description of subroutine dot here]
function dot(v1,v2,n,i,j)
   implicit none
   double precision :: v1, v2, dot, sum
   integer :: isum,i,j,n

   dimension v1(n,*),v2(n,*)
   
   !----this function performs the dot product of
   !      v1(k,i).v2(k,j)
   
   sum = 0.0d0
   do isum=1,n
      sum =sum +v1(isum,i)*v2(isum,j)
   end do
   dot =sum
   return
end function dot 
