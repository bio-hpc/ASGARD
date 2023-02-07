
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
!+ [Enter a one-line description of subroutine traco here]
subroutine traco(nr,npa,x,beta,it)
   implicit double precision (a-h,o-z)
   logical new
   
   !     ----- routine to make transformation from cartesian to oblique
   !           coordinates and vice versa -----
   
   dimension x(2)
   data new/.false./
   
   !     ----- initialize -----
   
   if (.not.new) then
      new = .true.
      one = 1.e0
      pye = 4.e0* atan(one)
      conv = 1.8e2/pye
      betar = beta/conv
      cosb = cos(betar)
      sinb = sin(betar)
   end if
   
   !     ----- perform the transformation -----
   
   if (it > 0) then
      
      i3 = 3*npa+1
      do 10 j = 1,nr
         xh = x(i3+2)/sinb
         x(i3+2) = xh
         x(i3) = x(i3)-cosb*xh
         i3 = i3+3
      10 continue
   else if (it < 0) then
      i3 = 3*npa+1
      do 20 j = 1,nr
         x(i3) = x(i3)+cosb*x(i3+2)
         x(i3+2) = x(i3+2)*sinb
         i3 = i3+3
      20 continue
   end if
   return
end subroutine traco 
