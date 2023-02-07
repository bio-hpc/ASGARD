
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
!+ [Enter a one-line description of subroutine chck here]
subroutine chck(x,ix,ih,n3,xdir,ipair,xbel)
   !     call chck      (x,ix,n3,xdir,4,x(mxbel))

   implicit double precision (a-h,o-z)
   dimension x(*),xdir(*),ix(*),xbel(*)
   character(len=4) ih(*)
#  include "pointer.h"
#  include "epot.h"
#  include "inpdat.h"
   
   ! --- make the move then check for step length
   
   do 10 i=1,n3
      xbel(i) = xbel(i) + xdir(i)
   10 continue
   
   energy = etot
   slen2 = ddot(n3,xdir,1,xdir,1)
   rmslen = dsqrt(slen2/dble(n3))
   write(6,11) rmslen
   
   !   --- scale down step length to rmax
   
   if(rmslen > smx) then
      alph1 = smx/rmslen
      do 20 i=1,n3
         xbel(i) = xbel(i) - (1.00d0 - alph1)*xdir(i)
         xdir(i) = alph1*xdir(i)
      20 continue
   end if
   
   !  --- get the energy at the new position, and repeatedly scale
   !      in order to get the energy change to be less than edmax
   
   do 30 loop=1,20
      ndrv = 0
      call forces(x,ix,ih,ipair,ndrv)
      write(6,45)energy,etot
      
      !        edtrue ---- change in true energy after this step
      
      edtrue = etot-energy
      if(edtrue < emx) return
      do 25 i=1,n3
         delta = (1.00d0 - alpha)*xdir(i)
         xbel(i) = xbel(i) - delta
         xdir(i) = alpha*xdir(i)
      25 continue
   30 continue
   
   !   --- scaling was unsuccessful if normal exit from do loop
   
   write(6,61)
   call mexit(6, 1)
   11 format(10x,'rms of step length = ',f12.10)
   45 format(10x,'ene @step k-1 = ',f15.6,'  ene @k = ',f15.6)
   61 format(10x,'... scaling was unsuccessful ... stop')
end subroutine chck 
