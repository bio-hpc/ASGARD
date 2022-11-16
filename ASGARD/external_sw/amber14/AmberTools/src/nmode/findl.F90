
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
!+ [Enter a one-line description of subroutine findl here]
subroutine findl(n,h,bdwnhl,tlamba,d,b,roots,vect)
   
   ! --- routine to find lambda for the down hill search
   
   implicit double precision (a-h,o-z)
   character(len=1) uplo,jobz
   data uplo,jobz / 'U','N' /
   dimension h(*),d(*),b(*),roots(*),vect(n,*)
   
   ntt = n*(n+1)/2
   do 10 i = 1,ntt
      d(i) = h(i)
   10 continue
   
   call dspev(jobz,uplo,n,d,roots,vect,n,b,ier)
   if(ier /= 0) then
      write(6,11) ier
      call mexit(6, 1)
   else
      bdif = roots(1) - bdwnhl
      tlamba = min(0.00d0,bdif)
      write(6,12) roots(1),roots(2),tlamba
   end if
   return
   11 format('... stop in routine findl... ier = ',i5)
   12 format(10x,'b1 = ',f10.5,' b2 =',f10.5,' tlamba = ',f10.5)
end subroutine findl 
