
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
!+ [Enter a one-line description of subroutine pairw here]
subroutine pairw (natom,npair,kpair,nhb,iar1,iar2)
   dimension iar1(*),iar2(*)
#  include "files.h"
   
   if (kpair == 1) then
      call amopen(56, 'PRLIST', owrite, 'F', 'W')
      write(56,2) natom,npair,nhb
      write(56,2) (iar1(i),i=1,natom)
      write(56,2) (iar2(i),i=1,npair)
      close (56)
   else
      call amopen(56, 'PRLIST', owrite, 'U', 'W')
      write(56) natom,npair,nhb
      write(56) (iar1(i),i=1,natom)
      write(56) (iar2(i),i=1,npair)
      close (56)
   end if
   
   write(6,1) npair
   1 format('0WRITING ',i6,' PAIRS TO FILE PRLIST (unit 56)')
   
   return
   2 format (12(1x,i5))
end subroutine pairw 
