
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
!+ [Enter a one-line description of subroutine chdirec here]
subroutine chdirec(vect,roots,xdir,n3,scale,b1,b2)
   
   !---- calculate which eigenvector to follow and scale eigenvalues
   !       if necessary
   
   implicit double precision (a-h,o-z)
#  include "inpdat.h"
   dimension roots(*),xdir(*),vect(n3,*)
   
   ! ---  get projection of the eigenvectors on the previous step's
   !        move vector
   
   dotmax = 0.0d0
   do 10 jdot = 1,ivect
      dotref = dabs(ddot(n3,xdir,1,vect,jdot))
      if(dotref > dotmax) then
         dotmax = dotref
         itemp = jdot
      end if
   10 continue
   
   if(itemp <= isdir) isdir = itemp
   if(roots(1) < 0.00d0) isdir = 1
   !     if(ith.eq.isw) isdir = 1
   if(isdir == 1) isec = 2
   if(isdir /= 1) isec = 1
   
   write(6,11) isdir,isec
   
   !--- scale eigenvalue spectrum if b(isdir) > 0.30*b(isec)
   
   if(roots(isdir) < (0.30d0*roots(isec))) then
      b1 = roots(isdir)
      b2 = roots(isec)
      scale = 1.00d0
   else
      scale = 0.5d0*dsqrt(abs(roots(isec)/roots(isdir)))
      write(6,12)scale
      b1 = scale*scale*roots(isdir)
      b2 = roots(isec)
   end if
   write(6,40) b1,b2
   return
   
   11 format('     isdir = ',i3,'  isec = ',i3)
   12 format(5x,'scale coordinates by ',f10.5,'  to find lambda')
   40 format('     b1 = ',f10.5,' b2 = ',f10.5)
end subroutine chdirec 
