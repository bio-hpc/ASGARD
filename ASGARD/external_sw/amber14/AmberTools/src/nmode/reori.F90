
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
!+ [Enter a one-line description of subroutine reori here]
subroutine reori(ve,nv,nc,crd,ndim,nrat)
   
   !     -----subroutine to reorient the eigenvector + atomic coordinates
   !      along the principle axes of the molecule - this way the
   !      system and reference vibrational motions can be directly
   !      compared, irrespective of rigid body rotation/translation
   !      of the system (relative to some other part of the total
   !      system, ie, use only when you have used excise with iat
   !      and jat other than 1 and natom, respectively.) Confused?
   !      me too!-----
   
   implicit double precision (a-h,o-z)
#  include "sizes2.h"
   common/mass/atmas(maxatom)
   dimension ve(ndim,nv),crd(3*maxatom),sum(500)
   
   !     -----first add the eigenvectors to the original coordinates----
   
   write(6,9008)
   do 500 l=1,nv
      do 400 m=1,nc
      400 ve(m,l) = ve(m,l) + crd(m)
      
      !       -----now reorient molecule-----
      
      call movecm2(ve(1,l),nrat,ndim,crd)
      
      !       -----subtract out the reoriented molecular coords
      
      do 425 m=1,nc
      425 ve(m,l) = ve(m,l) - crd(m)
      
      !       -----now normalize these new "vectors"-----
      
      do 450 li=1,nv
         sum(li)=0.0
         do 430 mi=1,nc,3
            ii = ((mi-1)/3)+1
            sum(li) = sum(li) + (ve(mi,li)*ve(mi,li) + &
                  ve(mi+1,li)*ve(mi+1,li) + ve(mi+2,li)*ve(mi+2,li)) &
                  *atmas(ii)
         430 continue
      450 continue
      do 480 li=1,nv
         do 470 mi=1,nc,3
            ve(mi,li) = ve(mi,li)/(dsqrt(sum(li)))
         470 continue
      480 continue
   500 continue
   return
   9008 format(5x, 'MOLECULE IS REORIENTED ALONG THE PRINCIPAL AXES')
end subroutine reori 
