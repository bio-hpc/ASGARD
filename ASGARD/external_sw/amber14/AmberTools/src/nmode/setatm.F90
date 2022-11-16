
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
!+ [Enter a one-line description of subroutine setatm here]
subroutine setatm(nat,nres,natb,ipres,igrp,igres)
   logical active
   dimension ipres(*),igrp(*),igres(*)
   
   idum = 0
   do 10 i = 1,nat
      if(igrp(i) > 0) idum = idum+1
   10 continue
   natb = idum
   
   do 20 i = 1,nres
      i1 = ipres(i)
      i2 = ipres(i+1)-1
      igres(i) = 0
      active = .false.
      do 15 j = i1,i2
         if(igrp(j) > 0) active = .true.
      15 continue
      if(active) igres(i) = 1
   20 continue
   return
end subroutine setatm 
