
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
!+ [Enter a one-line description of subroutine grpprj here]
subroutine grpprj(b,vect,freq,nint,ncart,nvect)
   implicit double precision (a-h,o-z)
#  include "anal.h"
   dimension b(ncart,nint),vect(ncart,nvect),freq(*)
   dimension prjct(2000)
   
   if(nint > 2000) then
      write(6,100) nint
      call mexit(6, 1)
   end if
   
   ! --- here start big loop over all eigenvectors
   
   do 10  ivec= ibeg,iend
      
      ! --- zero out all arrays
      
      do 15 izero = 1,nint
         prjct(izero)  = 0.00d0
      15 continue
      
      ! --- calculate the projection
      
      do 25 l=1,nint
         
         do 20 k = 1,ncart
         20 prjct(l) = prjct(l) + vect(k,ivec)*b(k,l)
         
      25 continue
      
      ! --- here do sort then print
      
      
      write(6,110) ivec,freq(ivec)
      write(6,111)
      write(6,112)
      
      ip = 1
      do 200 i=1,nrgrp
         r = prjct(ip)*prjct(ip) + prjct(ip+1)*prjct(ip+1) &
               + prjct(ip+2)*prjct(ip+2)
         r = dsqrt(r)
         write(6,120) i,(prjct(lprt),lprt=ip,ip+2),r &
               , (prjct(jprt),jprt=ip+3,ip+5)
         ip = ip+6
      200 continue
      
   10 continue
   return
   
   100 format('  nint > 2000 ... stop... nint =',i10)
   110 format(/ /' NORMAL MODE ',i3,'  WITH FREQUENCY =',f12.4/)
   111 format(27x,'TRANSLATION',27x,'ROTATION')
   112 format(17x,'X',9x,'Y',9x,'Z',9x,'R',9x,'X',9x,'Y',9x,'Z')
   120 format(' group',i5,7(2x,f8.3))
end subroutine grpprj 
