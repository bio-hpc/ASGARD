      subroutine se_gensub(ierror)

      implicit none
#include "sebomd.dim"
#include "se_elemnt.h"
#include "se_fulsys.h"
#include "se_logic.h"
#include "se_subsys.h"
      integer ierror

      integer i, iai
      integer natot
      integer ipoint
      integer nbig
      integer ibig

      nsub = 0
      natot = 0
      ipoint = 1
      nbig = 0
      ibig = 0

      if (stand) then
        nsub = 1
        iatom1(1) = 1
        do i=1,natoms
!         skip dummies and sparkles.
          iai = iatnum(i)
          if(natorb(iai).ne.0) then
            iatoms(ipoint) = i
            iabuff(ipoint) = 0
            ipoint = ipoint + 1
          endif
        end do
        iatom1(2) = ipoint
      elseif (clust) then
         call se_clustsub(ipoint, natot, nbig, ibig, ierror)
      elseif (resclust) then
         call se_resclustsub(ipoint, natot, nbig, ibig, ierror)
      endif

      end
