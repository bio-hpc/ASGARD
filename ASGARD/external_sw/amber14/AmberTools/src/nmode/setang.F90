
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
!+ [Enter a one-line description of subroutine setang here]
subroutine setang(nt,mt,it,jt,kt,ict,igrp)
   dimension it(2),jt(2),kt(2),ict(2),igrp(2)
   
   nta = 0
   mta = 0
   do 10 i = 1,nt
      iat = it(i)/3+1
      jat = jt(i)/3+1
      kat = kt(i)/3+1
      if(igrp(iat) > 0 .or. igrp(jat) > 0 .or. igrp(kat) > 0) &
            then
         if(i > mt) mta = mta+1
         nta = nta+1
         it(nta) = it(i)
         jt(nta) = jt(i)
         kt(nta) = kt(i)
         ict(nta) = ict(i)
      end if
   10 continue
   nt = nta
   mt = nta-mta
   return
end subroutine setang 
