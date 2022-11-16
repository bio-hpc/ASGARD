
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

!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine nmdfil here]
subroutine nmdfil
   
   !     OUTPUT: (to common)
   
#  include "files.h"
   
   !        MD files as described in amber v 3.0 Rev A documentation
   
   !     INTERNAL:
   
   character(len=80) arg
   !        ... temp for each of the whitespace delimited command line words
   integer iarg, narg
   !        ... arg pointer, final number of arguments
   
   !     --- initialize file names ---
   
   nmdin   = 'nmdin'
   nmdout  = 'screen'
   inpcrd = 'inpcrd'
   parm   = 'prmtop'
   restrt = 'restrt'
   refc   = 'refc'
   vecs   = 'vecs'
   lmod  = 'lmode'
   tstat = 'tstate'
   expfil = 'expfile'
   
   !     --- default output files: 'N'ew
   
   owrite = 'N'
#ifndef NOGETARG
   
   !     --- get com line arguments ---
   
# ifdef HITACHI_GETARG
   iarg = 1
# else
   iarg = 0
# endif
   indx = iargc()
   10 continue
   iarg = iarg + 1
   call getarg(iarg,arg)
   if (arg == '-O') then
      owrite = 'U'
   else if (arg == '-i') then
      iarg = iarg + 1
      call getarg(iarg,nmdin)
   else if (arg == '-o') then
      iarg = iarg + 1
      call getarg(iarg,nmdout)
   else if (arg == '-p') then
      iarg = iarg + 1
      call getarg(iarg,parm)
   else if (arg == '-c') then
      iarg = iarg + 1
      call getarg(iarg,inpcrd)
   else if (arg == '-r') then
      iarg = iarg + 1
      call getarg(iarg,restrt)
   else if (arg == '-ref' .or. arg == '-z') then
      iarg = iarg + 1
      call getarg(iarg,refc)
   else if (arg == '-v') then
      iarg = iarg + 1
      call getarg(iarg,vecs)
   else if (arg == '-l') then
      iarg = iarg + 1
      call getarg(iarg,lmod)
   else if (arg == '-t') then
      iarg = iarg + 1
      call getarg(iarg,tstat)
   else if (arg == '-e') then
      iarg = iarg + 1
      call getarg(iarg,expfil)
   else
      if (arg == ' ') goto 20
      write(6,'(/,5x,a,a)') 'unknown flag: ',arg
      write(6,9000)
      call mexit(6, 1)
   end if  ! (arg == '-O')
   if (iarg < indx) goto 10
   
   20 continue
   narg = iarg - 1
   if (narg < 1) then
      write(6,9000)
      call mexit(6, 1)
   end if
   
#endif
   return
   9000 format(/,5x,'usage: nmode [-O] -i nmdin -o nmdout -p prmtop', &
         ' -c inpcrd -r restrt', &
         /,18x,'-ref refc -v vecs -t tstate -l lmode -e expfile')
end subroutine nmdfil 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine nmdfil2 here]
subroutine nmdfil2
   
#  include "files.h"
   character(len=20) arg
   
   !     --- initialize file names ---
   
   nmdin   = 'nmdin'
   nmdout  = 'nmdout'
   inpcrd = 'no_inpcrd'
   parm   = 'prmtop'
   rvecs   = 'rvecs'
   vecs   = 'vecs'
   prjfil = 'proj'
   masfil = 'mass'
   
   !     --- default for output files: 'N'ew
   
   owrite = 'N'
#ifndef NOGETARG
   
   !     --- get com line arguments ---
   
# ifdef HITACHI_GETARG
   iarg = 1
# else
   iarg = 0
# endif
   indx = iargc()
   10 continue
   iarg = iarg + 1
   call getarg(iarg,arg)
   if (arg == '-O') then
      owrite = 'U'
   else if (arg == '-i') then
      iarg = iarg + 1
      call getarg(iarg,nmdin)
   else if (arg == '-o') then
      iarg = iarg + 1
      call getarg(iarg,nmdout)
   else if (arg == '-p') then
      iarg = iarg + 1
      call getarg(iarg,parm)
   else if (arg == '-r') then
      iarg = iarg + 1
      call getarg(iarg,rvecs)
   else if (arg == '-v') then
      iarg = iarg + 1
      call getarg(iarg,vecs)
   else if (arg == '-m') then
      iarg = iarg + 1
      call getarg(iarg,masfil)
   else if (arg == '-proj') then
      iarg = iarg + 1
      call getarg(iarg,prjfil)
   else
      if (arg == ' ') goto 20
      write(6,'(/,5x,a,a)') 'unknown flag: ',arg
      write(6,9000)
      call mexit(6, 1)
   end if
   if (iarg < indx) goto 10
   
   20 continue
   narg = iarg - 1
   if (narg < 2) then
      write(6,9000)
      call mexit(6, 1)
   end if
   
#endif
   return
   9000 format(/,5x, &
         'usage: nmanal [-O] -i nmdin -o nmdout -p prmtop -r rvecs', &
         ' -v vecs -m mass -proj proj')
end subroutine nmdfil2 
