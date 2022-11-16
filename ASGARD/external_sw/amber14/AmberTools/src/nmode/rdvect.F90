
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
!+ [Enter a one-line description of subroutine rdvect here]
subroutine rdvect(vect,freq,nvect,n3,nf,xj,ieff,ivform,langevin)
   
   !--- this routine will read in eigenvectors
   
   implicit double precision (a-h,o-z)
#  include "sizes2.h"
   character(len=40) title
   character(len=4) star
   real freq4,xw(3*maxatom)
   dimension vect(n3,nvect),xj(*),freq(*)
   
   if (ivform == 0) then
      read(nf,err=990) n
   else
      read(nf,1) title
      read(nf,2) n
   end if
   if(n /= n3) then
      write(6,5) n,n3
      call mexit(6, 1)
   end if
   if (ivform == 0) then
      read(nf,err=990) (xw(j),j=1,n3)
      do 40 k=1,n3
         xj(k) = xw(k)
      40 continue
      write(6,*) 'start of xj: ',xj(1),xj(2),xj(3),xj(4)
   else
      read(nf,3,err=990) (xj(j),j=1,n3)
   end if
   
   if( langevin == 1 ) return
   do 30 l = 1,nvect
      if (ivform == 0) then
         read(nf,end=90,err=990) ivec,freq4
         freq(l) = freq4
         write(6,*) l,freq(l)
         read(nf,err=990) (xw(k),k=1,n3)
         do 50 k=1,n3
            vect(k,l) = xw(k)
         50 continue
      else
         read(nf,4,end=90) star
         read(nf,2) ivec,freq(ivec),eff
         if (ieff == 1 .and. eff /= 0.0) freq(ivec) = eff
         read(nf,3) (vect(k,l),k=1,n3)
      end if
   30 continue
   
   1 format(a40)
   2 format(i5,2f12.5)
   3 format(7f11.5)
   4 format(a4)
   5 format (1x,' n and n3 do not match ',2i10)
   close(unit=nf)
   return
   90 write(6,6) l-1
   6 format(' Found only',i4,' eigenvectors on vecs file')
   nvect = l-1
   close(unit=nf)
   return
   990 write(6,*) 'Error in reading vecs file'
   call mexit(6,1)
end subroutine rdvect 
