
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
!+ [Enter a one-line description of subroutine mweight here]
subroutine mweight(dd,fg,xdir,ichoice,n3,amass,vecs,nvect)
   
   !-----mass weighting of coordinate system
   
   implicit double precision (a-h,o-z)
#  include "sizes2.h"
   logical first
   common/rm/rwinv(3*maxatom)
   dimension dd(*),fg(*),xdir(*),amass(*),vecs(n3,*)
   data first/.true./
   
   !---- ichoice = 0  convert into mass-weighted coordinate system
   !               1  unweight the coordinate
   !               2  unweight the normal mode vectors
   
   if (first) then
      if (n3 > 3*maxatom) then
         write(6,*) '3*MAXATOM is too small in subroutine mweight!'
         write(6,*) 'n3, MAXATOM =', n3, maxatom
         call mexit(6, 1)
      end if
      natom = n3/3
      k = 1
      do i=1,natom
         rwinv(k) = 1./sqrt(amass(i))
         rwinv(k+1) = rwinv(k)
         rwinv(k+2) = rwinv(k)
         k = k + 3
      end do
      first = .false.
   end if
   
   if(ichoice == 0) then
      k = 0
      do i = 1,n3
         fg(i) = fg(i)*rwinv(i)
         do j = 1,i
            k = k+1
            dd(k) = rwinv(i)*dd(k)*rwinv(j)
         end do
      end do

   else if (ichoice == 1) then
      do i=1,n3
         xdir(i) = xdir(i)*rwinv(i)
      end do

   else if (ichoice == 2) then
      do iv=1,nvect
         do i=1,n3
            vecs(i,iv) = vecs(i,iv)*rwinv(i)
         end do
      end do

   else
      write(6,*) 'Error in call to mweight'
      call mexit(6, 1)
   end if
   return
end subroutine mweight 
