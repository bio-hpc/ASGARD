
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
!+ [Enter a one-line description of subroutine setdih here]
subroutine setdih(np,mp,ip,jp,kp,lp,icp,igrp)
   dimension ip(*),jp(*),kp(*),lp(*),icp(*),igrp(*)
   
   npa = 0
   mpa = 0
   do 10 i = 1,np
      iat = ip(i)/3+1
      jat = jp(i)/3+1
      kat = iabs(kp(i))/3+1
      lat = iabs(lp(i))/3+1
      if(igrp(iat) > 0 .or. igrp(jat) > 0 .or. igrp(kat) > 0 &
            .or. igrp(lat) > 0) then
         if(i > mp) mpa = mpa+1
         npa = npa+1
         ip(npa) = ip(i)
         jp(npa) = jp(i)
         kp(npa) = kp(i)
         lp(npa) = lp(i)
         icp(npa) = icp(i)
      end if
   10 continue
   np = npa
   mp = npa-mpa
   return
end subroutine setdih 
