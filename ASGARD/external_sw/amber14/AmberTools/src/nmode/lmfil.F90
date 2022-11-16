!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine lmfil here]
subroutine lmfil
   !     Author: George Seibel
   !     gets unix command line input for program nmode
   !     implicit none
   
   !     OUTPUT: (to common)
   
#  include "files.h"
   
   !        MD files as described in amber v 3.0 Rev A documentation
   
   !     INTERNAL:
   
   character(len=80) arg
   !        ... temp for each of the whitespace delimited command line words
   integer iarg, narg
   !        ... arg pointer, final number of arguments
   
   !     --- initialize file names ---
   
   lmdin   = 'lmdin'
   lmdout  = 'lmdout'
   inpcrd = 'inpcrd'
   lmod  = 'lmode'
   plot  = 'plot'
   
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
      call getarg(iarg,lmdin)
   else if (arg == '-o') then
      iarg = iarg + 1
      call getarg(iarg,lmdout)
   else if (arg == '-c') then
      iarg = iarg + 1
      call getarg(iarg,inpcrd)
   else if (arg == '-l') then
      iarg = iarg + 1
      call getarg(iarg,lmod)
   else if (arg == '-p') then
      iarg = iarg + 1
      call getarg(iarg,plot)
   else
      if (arg == ' ') goto 20
      write(6,'(/,5x,a,a)') 'unknown flag: ',arg
      write(0,9000)
      call mexit(0, 1)
   end if
   if (iarg < indx) goto 10
   
   20 continue
   narg = iarg - 1
   if (narg < 2) then
      write(0,9000)
      call mexit(0, 1)
   end if
   
#endif
   return
   9000 format(/,5x,'usage: lmanal [-O] -i lmdin -o lmdout -c inpcrd ', &
         '-l lmode -p plot')
end subroutine lmfil 
