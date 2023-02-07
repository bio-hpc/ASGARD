
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
!+ [Enter a one-line description of subroutine langout here]
subroutine langout (neig, ind, mn, n, wr, wi, z, iclass, dnorm)
   
   !     ----- to output Langevin modes to the file LMODE
   
   implicit double precision (a-h,o-z)
   character(len=10) type
   dimension wr(n), wi(n), z(mn,n), ind(n), iclass(neig), dnorm(n)
   
   !     ----- write eigenvalues
   
   do 10 i = 1, n
      if (dnorm(ind(i)) /= 0.0) then
         write (11,'(1x, i4, 3x, f11.4, 2x, f11.4, 2x, f10.4)') &
               i, wr(ind(i))*108.59, wi(ind(i))*108.59, dnorm(ind(i))
      else
         write (11,'(1x, i4, 3x, f11.4, 2x, f11.4)') &
               i, wr(ind(i))*108.59, wi(ind(i))*108.59
      end if
   10 continue
   
   !     ----- write eigenvectors
   
   kb = 1
   do 30 k = kb, neig
      j = ind(k)
      if (iclass(k) == 1) then
         type = ' complex '
      else if (iclass(k) == 2) then
         type = ' real'
      else if (iclass(k) == 3) then
         type = ' imaginary'
      else if (iclass(k) == 4) then
         type = ' real'
      else
         type = ' zero'
      end if
      write (11,*)
      write (11,'(i4, a1, a10)') k, ':', type
      do 20 i = 1, n/2, 6
         write (11,'(6(2x,e11.4))') (z(i+l,j),l=0,5)
      20 continue
      write (11, '(2x, i10)') iclass(k)
   30 continue
   
   return
end subroutine langout 
