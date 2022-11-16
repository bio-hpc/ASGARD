      subroutine se_rdall(ierror)

!     driver subroutine for reading and processing all input.

      implicit none
#include "sebomd.dim"
#include "se_fulsys.h"
#include "se_logic.h"
#include "se_parproc.h"
#include "se_resdef.h"
      integer ierror
      integer i

      if(nres.gt.0)then
         if(screen) then
            write(iout,'(//" NUMBER OF RESIDUES = ",i6,
     .           /" RESIDUE POINTERS:")') nres
            write(iout,'(12i6)') (irpnt(i),i=1,nres)
            write(iout,*)
         endif
      endif

!     read in trailing lines of input file (subsetting etc)

      call se_rdtail(ierror)
      if (ierror.ne.0) return
!     if ((.not.recipintr).and.pme) then
!        write(iout,'(/" ERROR: NO GROUPS ARE SPECIFIED")')
!        ierror = 1
!        return
!     endif

      natoms3 = 3*natoms
      natoms4 = 4*natoms
      nres3 = 3*nres

      return
      end
