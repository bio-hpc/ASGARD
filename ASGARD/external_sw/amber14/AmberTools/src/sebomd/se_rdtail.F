      subroutine se_rdtail(ierror)

      implicit none
#include "sebomd.dim"
#include "se_logic.h"
#include "se_parproc.h"
#include "se_resdef.h"
#include "se_subpar.h"

      integer ierror, icoreln

! locals:
      logical readsome, nreadcl
      integer i

      nreadcl = .true.

      if (clust.or.resclust) then
        readsome = .false.
      else
        readsome = .true.
      endif

      i = 0
 100  i = i + 1
      if(clust.or.resclust.and..not.readsome)then
        if (ierror.ne.0) return
        icoreln = nres
        icorel1(1) = 1
        ncore(1) = nncore
        ncores = 1
        do i=1,nres
          icorel(i) = i
        enddo
        icorel1(ncores+1) = icoreln + 1
        readsome = .true.
        nreadcl = .false.
        goto 100
      endif
      
      if (readsome) then
        if ((clust.or.resclust).and.nreadcl) then
          write(iout,'(/" ERROR: NO CLUSTER PARAMETERS FOUND")')
          ierror = 1
        endif
        if (.not.stand) then
          if(nreadcl)then
            ierror=1
            write(iout,'(/" CHOOSE EITHER FOR A GRID BASED ",
     &              " SUBSETTING SCHEME A CLUSTER",/,
     &              " BASED SUBSETTING SCHEME, A",
     &              " COMBINED SUBSETTING SCHEME OR USE",
     &              " THE SUB .. END_SUB PARAMETERS")')
          endif
        endif
        return
      else
        ierror = 1
        write(iOUT,'(/" UNEXPECTED END-OF-FILE ENCOUNTERED")')
        return
      endif

      end
