
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
!+ [Enter a one-line description of subroutine setbon here]
subroutine setbon(nb,mb,ib,jb,icb,igrp)
   dimension ib(2),jb(2),icb(2),igrp(2)
   
   nba = 0
   nca = 0
   do 10 i = 1,nb
      iat = ib(i)/3+1
      jat = jb(i)/3+1
      if(igrp(iat) > 0 .or. igrp(jat) > 0) then
         if(i > mb) nca = nca+1
         nba = nba+1
         ib(nba) = ib(i)
         jb(nba) = jb(i)
         icb(nba) = icb(i)
      end if
   10 continue
   nb = nba
   mb = nba-nca
   return
end subroutine setbon 
