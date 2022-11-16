
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
!+ [Enter a one-line description of subroutine resdis here]
subroutine resdis(ia,ib,ja,jb,x,val)
   !     -- subroutine to calculate the smallest distance
   !        between the atoms of two residues
   
   implicit double precision (a-h,o-z)
   dimension x(3,*)
   
   
   val = 1.0e38
   do 100 i = ia,ib
      do 100 j = ja,jb
         scrach = (x(1,j)-x(1,i))**2 + (x(2,j)-x(2,i))**2 + &
               (x(3,j)-x(3,i))**2
         val = min (val,scrach)
         !          ii = ii+1
   100 continue
   
   return
end subroutine resdis 
