
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
!+ [Enter a one-line description of subroutine dotref here]
subroutine dotref(vect,rvec,nvect,nrcar,ncart, &
      ibeg,iend,ncary)
   
   !     ----- routine to take dot product of eigenvector with
   !            reference eigenvector---
   
   implicit double precision (a-h,o-z)
#  include "sizes2.h"
   common/ntrun5/nrvec,nrat,iat,jat,imov
   common/coords/c(3*maxatom),c1(3*maxatom)
   common/mass/atmas(maxatom)
   dimension vect(ncary,nvect),rvec(nrcar,nrvec), &
         dp(50,100),ia(500),ib(500)
   if(nrcar /= ncart) then
      write(6,999) nrcar,ncart
      call mexit(6, 1)
   end if
   if(imov == 1) then
      call reori(vect,nvect,ncart,c,ncary,nrat)
      call reori(rvec,nrvec,nrcar,c1,nrcar,nrat)
   end if
   ie = 1
   do 500 l=1,nrvec
      ib(ie)=l
      ie = ie+1
      id=1
      do 400 k=ibeg,iend
         ia(id) = k
         id=id+1
         dp(l,k) = 0.0
         do 300 m=1,ncart,3
            im = ((m-1)/3)+1
            dp(l,k) = dp(l,k) + (vect(m,k)*rvec(m,l) + &
                  vect(m+1,k)*rvec(m+1,l) + &
                  vect(m+2,k)*rvec(m+2,l))*atmas(im)
         300 continue
         !         if (abs(dp(l,k)).lt.0.2) dp(l,k) = 0.0
      400 continue
   500 continue
   nb=nrvec
   na=(iend-ibeg)+1
   write(6,997)
   call prtma(dp,ia,ib,na,nb,6)
   return
   999 format(1x,' number of atoms do not match in dotref', &
         2i9)
   998 format(/,1x,' dot product of ref vector',i4,' with vector', &
         i4,' equals',e20.5)
   997 format(/,10x,'DOT PRODUCT ARRAY',/,1x,'SYST',/,1x,'VECT', &
         10x,'REFERENCE VECTOR',/)
end subroutine dotref 
