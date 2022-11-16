
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
!+ [Enter a one-line description of subroutine nrch here]
subroutine nrch(x,ix,ih,h,f,xdir)
   
   implicit double precision (a-h,o-z)
   dimension ix(*)
   character(len=4) ih(*)
#  include "pointer.h"
#  include "inpdat.h"
   character(len=1) uplo
   data uplo/'U'/
   dimension x(*),xdir(*),h(*),f(*)
   data eps /1.d-4/
   
   !  ---  loop over maxcyc steps of Newton-Raphson
   
   ndrv = 2
   ipair = 4
   if (ibelly == 0) call movecm(x(mx),x(mamass),natom)
   do 50 ifn=1,maxcyc
      call mweight(h,f,x(1),0,ns3,x(mamass),dummy,dummy)
      if(ibelly == 0) call level(x(mh),x(mxbel),x(mamass),ns3,natom)
      
      !  --- add a constant into the diagonal elements of Hessian
      !          to make it positive defininte.
      
      if (ndiag /= 0) then
         call findl(ns3,x(mh),bdwnhl,tlamba,x(mdd),x(mb), &
               x(mroots),x(mvect))
      else
         tlamba = -bdwnhl
      end if
      jadd = 0
      do 10 iadd=1,ns3
         jadd= jadd + iadd
         h(jadd) = h(jadd) - tlamba
      10 continue
      
      call dppsv(uplo,ns3,1,h,f,ns3,info)
      if (info /= 0) then
         write(6,*) 'dppsv returns error code: ',info
         call mexit(6, 1)
      end if
      
      !  --- check for goodness of step
      
      call mweight(h,x(1),f,1,ns3,x(mamass),dummy,dummy)
      do 20 i=1,ns3
         xdir(i) = f(i)
      20 continue
      call chck(x,ix,ih,ns3,xdir,ipair,x(mxbel))
      
      !   ---  get, print the new energy, and save coord.
      
      if(ifn == maxcyc) ndrv = 1
      call forces(x,ix,ih,ipair,ndrv)
      call eprnt(ifn,natom,x(mf))
      call dihed(x(mx),ih(mgraph),idih)
      gmax = grdmax(x(mf),gnorm,ns3)
      if(mod(ifn,nsave) == 0) call savec(natom,x(mx),3)
      if(gnorm <= drms) return
   50 continue
   101 format(' ..... error in dppfa... info = ',i5)
   ! 102 format(' ier = ',i5,' irank = ',i5)
   return
end subroutine nrch 
