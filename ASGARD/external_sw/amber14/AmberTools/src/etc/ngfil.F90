subroutine ngfil
   
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
   
   ! Author: George Seibel
   ! gets unix command line input for program NUCGEN.
   
   use nucgen_files_module
   
   implicit none
   
#ifndef NOGETARG
   
   ! NUCGEN files as described in amber v 3.0 Rev A documentation
   !
   ! INTERNAL:
   
   ! ... temp for each of the whitespace delimited command line words
   character(len=80) :: arg
   
   ! arg pointer, final number of arguments
   integer :: iarg, narg, index
   
#endif
   
   ! --- initialize file names ---
   ngin   = 'ngin'
   ngout  = 'ngout'
   ngdat  = 'ngdat'
   pdbout = 'pdbout'
   
   ! --- default output file status: 'N'ew
   owrite = 'N'
   
#ifndef NOGETARG
!
!     --- get com line arguments ---
!
# ifdef HITACHI_GETARG
   iarg = 1
# else
   iarg = 0
# endif
   index = iargc()
   
   10 continue
   iarg = iarg + 1
   call getarg(iarg,arg)
   if (arg .eq. '-O') then
      owrite = 'U'
   else if (arg .eq. '-i') then
      iarg = iarg + 1
      call getarg(iarg,ngin)
   else if (arg .eq. '-o') then
      iarg = iarg + 1
      call getarg(iarg,ngout)
   else if (arg .eq. '-d') then
      iarg = iarg + 1
      call getarg(iarg,ngdat)
   else if (arg .eq. '-p') then
      iarg = iarg + 1
      call getarg(iarg,pdbout)
   else
      if (arg .eq. ' ') go to 20
      write(6,'(/,5x,a,a)') 'unknown flag: ',arg
      write(6,9000)
      call mexit(6, 1)
   end if
   
   if (iarg .lt. index) go to 10
   
   20 continue
   
   narg = iarg - 1
   if (narg .lt. 2) then
      write(6,9000)
      call mexit(6, 1)
   end if
   
#endif
   
   return
   
   9000 format(/5x, 'usage: nucgen [-O] -i ngin -o ngout -d ngdat -p pdbout')
   
end subroutine ngfil
