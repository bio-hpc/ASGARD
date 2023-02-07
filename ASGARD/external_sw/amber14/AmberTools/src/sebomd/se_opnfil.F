      subroutine se_opnfil(ifile,ierror)

!     routine to open file number ifile.  see block data routine for
!     unit numbers, file names, etc.  error messages are written
!     to the main output file and to standard output.  this assumes
!     that the main output file is opened before any other file, and
!     that no problems are encountered opening it.

      implicit none
      integer ierror, ifile
#include "sebomd.dim"
#include "se_funit.h"
#include "se_fopen.h"
#include "se_parproc.h"
c
      ierror = 0
      if(ifile.lt.1.or.ifile.gt.nfiles)then
        ierror = 1
        write(iout,10) ifile
  10    format(/' ERROR IN OPNFIL -- FILE NUMBER ',i5,' OUT OF RANGE'/)
        go to 1000
      endif
      open(iunit(ifile),file=fname(ifile),status=fstat(ifile),err=100)
      go to 1000
 100  ierror = 1
      write(iout,110) iunit(ifile),fname(ifile)
  110 FORMAT(/' ERROR IN OPNFIL -- CANNOT OPEN UNIT ',i3,' FILE ',a20/)
 1000 return
      end
