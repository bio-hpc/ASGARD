
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
!+ [Enter a one-line description of subroutine getcmn here]
subroutine getcmn(nr,x,ntx,box,beta)
   implicit double precision (a-h,o-z)
   
   !     ----- ROUTINE TO READ THE COORDINATES -----
   
   common/runhed/ihead(20),ihead1(20)
   
   dimension x(*),box(*)
#  include "files.h"
   
   nr3 = 3*nr
   
   write(6,91)
   call amopen(50,vecs,'O','F','R')
   read(50,40) ihead1
   read(50,94) nr3
   natom = nr3/3
   if(natom /= nr) goto 200
   read(50,92) (x(i),i=1,nr3)
   close(50)
   return
   200 continue
   write(6,1000)
   call mexit(iout, 1)
   40 format(20a4)
   91 format(/ /,'   3.  A T O M I C   C O O R D I N A T E S',/ /)
   92 format(7f11.5)
   94 format(i5,f10.5)
   1000 format(/10x,'****** ERROR IN INPUT STOPPED IN GETCMN ******')
end subroutine getcmn 
