
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
!+ [Enter a one-line description of subroutine xconst here]
subroutine xconst(e,x,f,dd,ndrv,xref,wref,natom)
   implicit double precision (a-h,o-z)
   
   !     ----- routine for constraints of form (1/2)*wtcons*(x-xcons)**2
   
   dimension x(*),f(*),dd(*),xref(*),wref(*)
   
   e = 0.0d+00
   
   nr = 3*natom
   do i=1,nr
      wtcons = wref(((i-1)/3)+1)
      del = x(i) - xref(i)
      f(i) = f(i) - wtcons*del
      e = e + 0.5*wtcons*del*del
   end do
   if (ndrv == 2) then
      k = 0
      do i=1,nr
         k = k + i
         dd(k) = dd(k) + wref(((i-1)/3)+1)
      end do
   end if
   return
end subroutine xconst 
