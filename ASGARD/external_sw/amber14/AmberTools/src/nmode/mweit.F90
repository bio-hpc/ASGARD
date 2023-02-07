
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
!+ [Enter a one-line description of subroutine mweit here]
subroutine mweit (n, neig, ind, z, amass, winv)
   implicit double precision (a-h,o-z)
   
   dimension ind(n), z(n,n), amass(n/6), winv(n/2)
   
   k = 1
   do 10 i=1,n/6
      winv(k) = 1.0/sqrt(amass(i))
      winv(k+1) = winv(k)
      winv(k+2) = winv(k)
      k = k + 3
   10 continue
   
   do 30 kj = 1, neig
      j = ind (kj)
      do 20 i = 1, n/2
         z(i,j) = z(i,j) * winv(i)
      20 continue
   30 continue
   
   return
end subroutine mweit 
