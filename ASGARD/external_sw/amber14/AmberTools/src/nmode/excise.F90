
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
!+ [Enter a one-line description of subroutine excise here]
subroutine excise(vect,ncart,nvect)
   
   !     -----routine to excise part of the eigenvector to be compared
   !          to the reference vector
   
   implicit double precision (a-h,o-z)
#  include "sizes2.h"
#  include "infoa.h"
   common/ntrun5/nrvec,nrat,iat,jat,imov
   common/mass/atmas(maxatom)
   common/coords/c(3*maxatom),c1(3*maxatom)
   dimension vect(ncart,nvect),sum(900)
   
   if(iat == 1.and.jat == natom) then
      write(6,97)
      do 100 l=1,natom
      100 atmas(l) = amass(l)
      return
   end if
   nat = (jat-iat)+1
   if(nrat /= nat) then
      write(6,96) nrat,nat
      call mexit(6, 1)
   end if
   ncar = nat*3
   do 500 l=1,nvect
      i3 = ((iat-1)*3)+1
      ia = iat
      do 400 k=1,ncar,3
         vect(k,l) = vect(i3,l)
         vect(k+1,l) = vect(i3+1,l)
         vect(k+2,l) = vect(i3+2,l)
         c(k) = c(i3)
         c(k+1) = c(i3+1)
         c(k+2) = c(i3+2)
         ik = ((k-1)/3)+1
         atmas(ik) = amass(ia)
         ia = ia+1
      400 i3=i3+3
   500 continue
   
   !     ----- normalize this new vector ---
   
   do 700 l=1,nvect
      sum(l) = 0.0
      do 600 k=1,ncar,3
         ik = ((k-1)/3)+1
         sum(l) = sum(l) + (vect(k,l)*vect(k,l) + vect(k+1,l)* &
               vect(k+1,l) + vect(k+2,l)*vect(k+2,l))*atmas(ik)
      600 continue
   700 continue
   do 900 l=1,nvect
      do 800 k=1,ncar
         vect(k,l) = vect(k,l)/(dsqrt(sum(l)))
      800 continue
   900 continue
   ncart=ncar
   return
   96 format(1x,'  nrat does not match nat in excise:',2i8)
   97 format(1x,'  the entire system vector will be used in dotref')
end subroutine excise 
