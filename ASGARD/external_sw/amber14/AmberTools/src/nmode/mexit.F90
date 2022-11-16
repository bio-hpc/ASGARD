
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine mexit here]
subroutine mexit(output_unit, status)
   
   !  mexit() - machine-dependent exit() procedure, designed to return an
   !            appropriate (success/failure) value to the operating system.
   
   !************************************************************************
   !                              AMBER                                   **
   !                                                                      **
   !               Copyright (c) 1986, 1991, 1995, 1997, 1999             **
   !                Regents of the University of California               **
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
   
   implicit none
   integer output_unit  ! close this unit if greater than zero
   integer status       ! exit status; error if non-zero


#ifdef IBM_VM_CMS
   if (output_unit /= 0 .and. output_unit /= 6) then
      close(unit=output_unit)
   end if
#else
   if (output_unit > 0) then
      close(unit=output_unit)
   end if
#endif

#if XLF90 || IBM3090 || F2C
   if (status /= 0) then
      stop 1
   else
      stop 0
   end if
#else
#ifdef UXP
   call setrcd(status)
#else
   call exit(status)
#endif
#endif
end subroutine mexit 



