
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
!+ [Enter a one-line description of subroutine savec here]
subroutine savec(natom,x,ipoint)
   
   ! ---- this little subroutine saves coordinates of the entire molecule
   
   implicit double precision (a-h,o-z)
   dimension x(*)
#  include "files.h"
#  include "inpdat.h"
   
   write(6,'(t2,''in savec point = '',i3)')ipoint
   if (ntxo == 0) then
      open(unit=23,file=restrt,form='unformatted',status='unknown')
      write(23) ititl
      nr3 = 3*natom
      write(23) natom
      write(23) (x(i),i=1,nr3)
   else
      open(unit=23,file=restrt,form='FORMATTED',status='unknown')
      write(23,1) ititl
      nr3 = 3*natom
      write(23,2) natom
      write(23,3) (x(i),i=1,nr3)
      1 format(20a4)
      2 format(i5)
      3 format(6f12.7)
   end if
   close(23)
   return
end subroutine savec 
