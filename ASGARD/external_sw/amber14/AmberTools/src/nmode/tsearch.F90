
!************************************************************************
!                              AMBER                                   **
!                                                                      **
!                  Copyright (c) 1986, 1991, 1995, 1997                **
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
!+ [Enter a one-line description of subroutine tsearch here]
subroutine tsearch(dd,vect,       roots,     b, h,n3, &
      natbel,xdir,    fg,x,ix)
   !           call tsearch(dd1,x(mvect),x(mroots),x(mb),h,ns3,
   !                        natsys,x(mxdir),fg,x,ix)
   
   !----this subroutine performs a single  step in the search
   !       based on the Cerjan-Miller proc.
   !          See Simons, Taylor, Jorgensen and Ozment,
   !          J. Phys. Chem. 87,2745-2753(1983)
   
   implicit double precision(a-h,o-z)
   dimension x(*)
   integer ix(*)
   character(len=1) uplo,jobz
   data uplo,jobz / 'U','V' /

#  include "pointer.h"
#  include "inpdat.h"

   dimension dd(*),vect(n3,*),roots(*),b(n3,*)
   dimension h(*),xdir(*),fg(*)
   
   !     xdir ---- displacement vector
   !     n3   ---- 3*number of atoms in moving part of molecule
   !     fg   ---- gradient
   !     slen2---- step length of move
   !     rmslen---  "     "        "   expressed as rms atomic movement
   
   
   !  --- find the lowest ivect eigenvalues and vectors
   
   call dspev(jobz,uplo,n3,dd,roots,vect,n3,b,ier)
   if(ier /= 0) then
      write(6,1)
      call mexit(6, 1)
   end if
   
   !  --- choose lambda and the search direction
   
   sqrnat = dsqrt(dble(natbel))
   if (istart == 0) then
      alph1 = hnot*sqrnat
      do 10 i=1,n3
         xdir(i)= idir*alph1*vect(i,isdir)
      10 continue
      istart = 1
   else
      call chdirec(vect,roots,xdir,n3,scale,b1,b2)
      call chlamba(tlamba,vect,roots,scale,fg,xdir, &
            n3,b1,b2)
      j = 0
      do 20 i = 1,n3
         j = j + i
         h(j) = h(j) - tlamba
      20 continue
      call dspsv(uplo,n3,1,h,ix(mkpvt),fg,n3,info)
      do 30 i=1,n3
         xdir(i) = xdir(i) - fg(i)
      30 continue
   end if
   call mweight(h,fg,xdir,1,ns3,x(mamass),x(mvect),nvect)
   
   !  --- now check for goodness of step
   
   call chck(x,ix,n3,xdir,4,x(mxbel))
   return
   
   1 format(10x,'... error in giveis.... stop.....')
end subroutine tsearch 
