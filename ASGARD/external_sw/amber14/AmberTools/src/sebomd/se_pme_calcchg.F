      subroutine se_pme_calcchg

      implicit none
#include "sebomd.dim"
#include "se_fulsys.h"
#include "se_elemnt.h"
#include "se_xyzcrd.h"
#include "se_glbmat1.h"

      integer i, ii, iai, norbsi, iiorb
      double precision zi, psumii

!     get the Mulliken charges
      ii = 0
      do i=1,natoms
        iai = iatnum(i)
!       skip dummy atom.
        if(iai.ne.0) then
          norbsi = natorb(iai)
          zi = zchg(iai)
          psumii = 0.0d0
          if(norbsi.gt.0)then
            do iiorb=1,norbsi
              ii = ii + iiorb
              psumii = psumii + pdiag(ii)
            end do
          endif
          atchg(i) = zi - psumii
          atchg2(i) = atchg(i)
          atchg3(i) = atchg(i)
        endif
      end do
      end
