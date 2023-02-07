
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
!+ [Enter a one-line description of subroutine percon here]
subroutine percon(rij2,xij)
   implicit double precision (a-h,o-z)
   common/setbox/box(3),beta,boxh(3),boxhm,boxhm2,boxoh,boxoq, &
         cosb,cosb2,ntm,ntb
   dimension xij(3)
   
   !     ----- periodic boundary condition through minimum image
   !           convention -----
   
   if(rij2 < boxhm2) return
   if(ntm /= 0) call traco (1,0,xij,beta,1)
   do 10 m = 1,3
      if (xij(m) >= boxh(m)) then
         xij(m) = xij(m)-box(m)
      else if (xij(m) < -boxh(m)) then
         xij(m) = xij(m)+box(m)
      end if
   10 continue
   if(ntm /= 0) call traco (1,0,xij,beta,-1)
   if(ntb <= 0) then
      if(dabs(xij(1))+ dabs(xij(2))+ dabs(xij(3)) >= boxoq) then
         do 20 m = 1,3
            xij(m) = xij(m)- sign(boxoh,xij(m))
         20 continue
      end if
   end if
   rij2 = 0.e0
   do 30 m = 1,3
      rij2 = rij2+xij(m)**2
   30 continue
   return
end subroutine percon 
