
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
!+ [Enter a one-line description of subroutine prtma here]
subroutine prtma(d,ia,ib,na,nb,iw)
   
   !     ----- PRINT OUT A MATRIX -----
   
   implicit double precision (a-h,o-z)
   dimension d(50,100),ia(*),ib(*)
   max = 4
   imax = ib(1)-1
   100 imin = imax+1
   imax = imax+max
   if (imax > ib(nb)) then
      idif = imax - ib(nb)
      max = max - idif
      imax = ib(nb)
   end if
   write (iw,9028) (i,i = imin,imax)
   do 160 i = imin,imax,max
      do 150 k=1,na
         icnt = ia(k)
         iend = i+(max-1)
         write(iw,9048) icnt,(d(lp,icnt),lp=i,iend)
      150 continue
   160 continue
   if (imax < ib(nb)) goto 100
   return
   9028 format(9x,10(7x,i3,7x))
   9048 format(i5,1x,10f17.7)
end subroutine prtma 
