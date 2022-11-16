
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
!+ [Enter a one-line description of subroutine errstp here]
subroutine errstp (mode,info)
   
   !     ----- for error handling
   
   character(len=6) v
   
   if (mode == 1) then
      write (6,1) info
   else if (mode == 2) then
      write (6,2) info
   else if (mode == 3) then
      write (6,3) info
   else if (mode == 4) then
      if (info > 0) then
         v = 'value '
      else
         info = - info
         v = 'vector'
      end if
      write (6,4) info, v
   else if (mode == 5) then
      write (6,5)
   else if (mode == 6) then
      write (6,6)
   else if (mode == 7) then
      write (6,7) info
   else  if (mode == 8) then
      if (info > 0) then
         v = 'value '
      else
         info = - info
         v = 'vector'
      end if
      write (6,8) info, v
   else if (mode == 9) then
      write (6,9)
   else
      write (6,98)
   end if  ! (mode == 1)
   
   write (6,99)
   call mexit(6, 1)
   
   !     ----- error messages
   
   1 format (/,5x, ' Forget it.  Insufficient memory even to',/, &
         5x, ' hold the small arrays. ',/, &
         5x, ' Memory needed so far = ', i15)
   2 format (5x, ' Insufficient memory to hold ',/, &
         5x, ' nonbonded pair list',/, &
         5x, ' Memory available for list = ', i15)
   3 format (5x, ' Memory overflow ',/, &
         5x, ' Total memory required = ', i15, &
         2x,'plus space for non-bon list')
   4 format (5x, ' Iteration for ', i6, ' -th eigen', a6, ' failed' &
         ,/,  5x, 'during diagonalization of chi')
   5 format (5x, ' Disallowed input value: ntrun = ', i5 )
   6 format (5x, ' programming error : dissallowed value for',/, &
         5x, 'an internal flag')
   7 format (5x, ' memory insufficient for running dynamics', /, &
         5x, 'memory required = ', i15)
   8 format (5x, ' Iteration for ', i6, ' -th eigen', a6, ' failed' &
         ,/,  5x, 'during diagonalization of beta')
   9 format (5x, 'Incorrect format flag for one of the files', &
         /,  5x, 'set 1 for FORMATTED and 0 for unformatted')
   98 format (5x, 'Something wrong somewhere')
   99 format (/,5x,  '     ***** ERROR STOP *****')
   
   !     7 and 9 are probably not used
   
end subroutine errstp 
