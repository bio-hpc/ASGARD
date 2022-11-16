
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
!+ [Enter a one-line description of subroutine pairr here]
subroutine pairr (natom,npair,kpair,nhb,iar1,iar2)

   dimension iar1(*),iar2(*)
   
   if (kpair == 1) then
      call amopen(56, 'PRLIST', 'O', 'F', 'R')
      read(56,3) latom,npair,nhb
      if (latom /= natom) then
         write(6,2) latom,npair,nhb
         call mexit(6, 1)
      end if
      read(56,3) (iar1(i),i=1,natom)
      read(56,3) (iar2(i),i=1,npair)
   else
      call amopen(56, 'PRLIST', 'O', 'U', 'R')
      read(56) latom,npair,nhb
      if (latom /= natom) then
         write(6,2) latom,npair,nhb
         call mexit(6, 1)
      end if
      read(56) (iar1(i),i=1,natom)
      read(56) (iar2(i),i=1,npair)
   end if
   write(6,1) npair
   1 format(' Reading ',i6,' pairs from file PRLIST')
   close (56)
   
   return
   2 format('0ERROR IN READING NONBON LIST:',3i6)
   3 format (12(1x,i5))
end subroutine pairr 
