
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
!+ [Enter a one-line description of subroutine group here]
subroutine group(b,nrgrp,nint,ncart)
   implicit double precision (a-h,o-z)
#  include "sizes2.h"
#  include "infoa.h"
   common/grp/jbelly,jatbel,jgrp(400),jgrp1(400),jgres(400)
   common/belly/ibelly,natbel,igroup(maxatom),isymbl(maxatom), &
         itree(maxatom),igres(600)
   dimension b(ncart,nint),xbel1(100)
   dimension weit(2),xc(2),d(100,6),d1(400,6)
   
   !     -- define some macros in fortran's subtle way
   !        (xbel1 will be initialized in setcor):
   
   xx(i) = xbel1(3*(i-1)+1)
   yy(i) = xbel1(3*(i-1)+2)
   zz(i) = xbel1(3*(i-1)+3)
   
   !  --- zero out b matrix
   
   do 1 i=1,nint
      do 2 j=1,ncart
      2 b(j,i) = 0.00d0
   1 continue
   
   !  --- read in a rigid group at a time in full atom space
   
   int = 0
   do 10 i = 1,nrgrp
      call rgroup2(natom,jatbel,nres,ngrp,ipres,lbres,igraph, &
            isymbl,itree,jgrp,weit,xc,.false.,.true.,4)
      call setatm(natom,nres,jatbel,ipres,jgrp,jgres)
      
      ! -- convert into belly space
      
      call stigrp(natbel,igroup,jgrp,jgrp1)
      call setcor(natbel,jgrp1,x,xbel1)
      n = jatbel*3
      
      ! --- now calculate the dmatrix (C. Eckart - Phys. Rev. 47,552(1935)
      !     Here assume that all atoms have an atomic mass of 1.00d0
      
      do 20  j = 1, n
         do 20 k = 1, 6
      20 d(j,k) = 0.00d0
      
      iat=1
      do 30 l=1,n,3
         d(l,1)     =  1.000d0/jatbel
         d(l+1,2)   =  1.000d0/jatbel
         d(l+2,3)   =  1.000d0/jatbel
         d(l,4)     = -yy(iat)/jatbel
         d(l+1,4)   =  xx(iat)/jatbel
         d(l+1,5)   = -zz(iat)/jatbel
         d(l+2,5)   =  yy(iat)/jatbel
         d(l,6)     =  zz(iat)/jatbel
         d(l+2,6)   = -xx(iat)/jatbel
         iat = iat + 1
      30 continue
      
      ! -- unpack d  and put in b for projection
      
      do 40 k = 1,6
         k1 = int + k
         call putvar(natbel,jgrp1,d(1,k),b(1,k1))
      40 continue
      int = int + 6
   10 continue
   return
end subroutine group 
