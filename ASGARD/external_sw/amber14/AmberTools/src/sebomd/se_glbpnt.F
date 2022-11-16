      subroutine se_glbpnt(ierror)

!     sets up pointers to access diatomic and diagonal blocks of
!     global matrices.  should be called after subsystems are defined
!     and the bonded atom pairlist is created.
!
!     error flag is set to 1 if storage limits for global matrices
!     would be exceeded.

      implicit none
#include "sebomd.dim"
#include "se_bpairs.h"
#include "se_elemnt.h"
#include "se_logic.h"
#include "se_fulsys.h"
#include "se_glbmat2.h"
#include "se_parproc.h"

      integer ierror

      integer iim, norbsi, iatm, iai
      integer ijm, norbsj, jatm, iaj
      integer j, j1, j2
      integer i, ipmax

      logical first
      data first /.true./
      save first

      ierror = 0

!     first set up pointers for entire system:


!     iimat(k) --> start of diagonal block in hdiag, pdiag, and fdiag
!                  for atom k (only needs to be determined once).

!     ijmat(i) --> start of diatomic blocks in hdiat, pdiat, and
!                  fdiat for bonded atom pairlist entry i.

!     ijrep(i) --> start of diatomic repulsions in eerep for full
!                  pairlist entry i (done only once, and only if
!                  it's not a direct calculation).


      if(first)then

!       assign iimat and ijrep.

        first = .false.
        iim = 1
        do 10 i=1,natoms
          iimat(i) = iim
          norbsi = natorb(iatnum(i))
          iim = iim + (norbsi*(norbsi+1))/2
 10     continue
        iimat(natoms+1) = iim
        if(iim-1.gt.mxdiag)then
          ierror = 1
          write(iout,'(/" MAXIMUM STORAGE EXCEEDED FOR DIAGONAL BLOCKS",
     .                  " OF GLOBAL H, F, P MATRICES"/" -- INCREASE",
     .                  " MXDIAG PARAMETER IN sebomd.dim TO AT LEAST ",
     .                  I6)') iim
        endif
      endif

!     assign ijmat.  this gets done whenever susbystems are updated
!     and is based on the bonded atom pairlist ipair.

      if(natoms.eq.1) return
      ijm = 1
      do 200 iatm=2,natoms
        j1 = ip1(iatm)
        j2 = ip1(iatm+1)-1

!       skip iatm if it nas no pairs.

        if(j2.lt.j1) go to 200
        iai = iatnum(iatm)
        norbsi = natorb(iai)
        do 180 j=j1,j2
          jatm = ipair(j)
          iaj = iatnum(jatm)
          norbsj = natorb(iaj)
          ijmat(j) = ijm
          ijm = ijm + norbsi*norbsj
 180    continue
 200  continue
      if(ijm-1.gt.mxdiat)then
        ierror = 1
        write(iout,'(/" MAXIMUM STORAGE EXCEEDED FOR DIATOMIC BLOCKS",
     .                " OF GLOBAL H, F, P MATRICES"/" -- INCREASE",
     .                " MXDIAT PARAMETER IN sebomd.dim TO AT LEAST ",
     .                  I7)') ijm
      endif
      ipmax = ip1(natoms+1)
      ijmat(ipmax) = ijm
      return
      end
