
!************************************************************************
!                              AMBER                                   **
!                                                                      **
!                  Copyright (c) 1986, 1991, 1995, 1997                **
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
!+ [Enter a one-line description of subroutine tstatf here]
subroutine tstatf(x,ix,ih,fg,h,dd1,vect)
   implicit double precision (a-h,o-z)
   dimension x(*) ,ix(*)
   character(len=4) ih(*)
   dimension fg(*),h(*),dd1(*),vect(*)
#  include "pointer.h"
#  include "inpdat.h"
   
   !    ----- this is a driver program for the transition state procedure
   
   
   ith = 1
   ntt = (ns3*(ns3+1))/2
   if(ibelly == 0) call movecm(x(mx),x(mamass),natsys)
   if (iflag == 0) call savei(x(mx),x(mxinit),ns3)
   
   !----major loop for search routine begins here
   
   do 100 ith=1,maxcyc
      
      ndrv = 2
      ipair = 4
      call forces(x,ix,ih,ipair,ndrv)
      gmax = grdmax(fg,rms,ns3)
      if (mod(ith,nprint) == 0) then
         call eprnt(ith,natom,x(mf))
         call dihed(x(mx),ih(mgraph),idih)
      end if
      do i = 1,ns3
         fg(i) = - fg(i)
      end do
      
      !  ---if we reach a zero gradient, or if we run out of steps,
      !        save the coordinates
      
      if((ith /= 1.and.rms < drms) .or. ith == maxcyc) then
         call savet(natom,x(mx))
         
         !  ---if we haven't run out steps, presumably must be at a transition
         !       state:
         
         if (ith < maxcyc) then
            write (6,5)
            5 format (/ /'Found transition state:'/ /)
            call eprnt(ith,natom,x(mf))
            call dihed(x(mx),ih(mgraph),idih)
            if (iflag == -1) then
               call mexit(6, 1)
            end if
            
            !*********************************************************************
            
            !    Following section is mainly for convenience:  starting from
            !      this transition state, go downhill to the next minimum.
            !      The "dirch" routine attempts to get things started so that
            !      the program will not return to the previous minimum, but
            !      this is not foolproof.  (Still, usually works.)
            
            !----       Now take the first step starting from a
            !           transition state and set up variables for the nrch routine
            
            call mweight(h,fg,x(mxdir),0,ns3,x(mamass), &
                  x(mvect),nvect)
            if (ibelly == 0) call level(h,x(mx),x(mamass),ns3,natsys)
            
            !  --- find the lowest eigenvector to take a step
            
            call findl(ns3,h,bdwnhl,tlamba,dd1,x(mb),x(mroots), &
                  x(mvect))
            call mweight(h,fg,x(mxdir),1,ns3,x(mamass), &
                  x(mvect),nvect)
            alph1 = hnot*dsqrt(dble(natsys))
            
            ! --- now chose idir for downhill direction
            
            call dirch(x(mvect),x(mx),ns3,idir,x(mxinit),x(mxdir), &
                  alph1)
            
            !---- check for goodness of step, then call nrch to complete the
            !       minimization
            
            call chck(x,ix,ns3,x(mxdir),4,x(mxbel))
            call nrch(x,ix,h,fg,x(mxdir))
            write (6,6)
            6 format(/ /'Found minimum:'/ /)
            call eprnt(ith,natom,x(mf))
            call dihed(x(mx),ih(mgraph),idih)
            call savec(natom,x(mx),4)
            
            !*********************************************************************
            
         else
            
            !  ---we have run out of steps:
            
            write(6,*) ' Did not find transition state....quitting.'
         end if  ! (ith < maxcyc)
         return
      else
         
         !  ---we have not yet reached a transition state: call "tsearch"
         !       to take a single step:
         
         call mweight(h,fg,x(mxdir),0,ns3,x(mamass), &
               x(mvect),nvect)
         if (ibelly == 0) call level(h,x(mx),x(mamass),ns3,natsys)
         do 30 i=1,ntt
            dd1(i) = h(i)
         30 continue
         call tsearch(dd1,x(mvect),x(mroots),x(mb), &
               h,ns3,natsys,x(mxdir),fg,x,ix)
      end if
   100 continue
   
end subroutine tstatf 
