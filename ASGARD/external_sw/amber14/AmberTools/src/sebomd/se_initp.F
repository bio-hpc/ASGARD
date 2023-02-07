      subroutine se_initp(init,ierror)
!
!     initializes global density matrix.  this routine should be called
!     after each pairlist update and global pointer update.
!
!     if init=0 then the density matrix will be fully initialized to
!     diagonal form.  otherwise, the old density matrix will be used
!     to build the new one.
!
!     note that this routine uses old pairlist information stored in
!     common /pairij/ to extract diatomic blocks from the old density
!     matrix.  common /pairij/ reverts back to its normal purpose with
!     the next call to ijmake.
!
      implicit none
#include "sebomd.dim"
#include "se_bpairs.h"
#include "se_elemnt.h"
#include "se_glbmat1.h"
#include "se_glbmat2.h"
#include "se_logic.h"
#include "se_fulsys.h"
#include "se_pairij.h"
#include "se_convrg.h"

      integer init, ierror

      integer i, iai, norbsi, ii, iatm
      integer jj, jnew, j1new, j2new, jatm
      integer jmatch, jold, j1old, j2old
      integer ij, ipmax, maxold, maxnew
      integer ij1new, ij2new, ijshft
      double precision pii

!     create a diagonal density matrix if init=0.  otherwise use information
!     from the old density matrix.
!
      logical first
      data first /.true./
      save first

      if(init.eq.0)then
        if (first.and.guess) then
!         call rddmx(ierror)
          if (ierror.ne.0) return
        else
!         diagonal blocks:
          do i=1,natoms
            iai = iatnum(i)
            norbsi = natorb(iai)
            if(norbsi.ne.0) then
!             divide the number of valence electrons evenly among the diagonal
!             atomic orbital elements.
              pii = dfloat(nvelec(iai))/dfloat(norbsi)
              ij = iimat(i)
              do ii=1,norbsi
                do jj=1,ii
                  pdiag(ij) = 0.0d0
                  if(jj.eq.ii) pdiag(ij) = pii
                  ij = ij+1
                enddo
              enddo
            endif
          enddo
!         off-diagonal blocks:
          ipmax = ip1(natoms+1)
          do i=1,ijmat(ipmax)-1
            pdiat(i) = 0.0d0
          enddo
        endif
      else
!       use what we can from the old density matrix.
        maxold = ijold(ip1old(natoms+1))-1
        do ij=1,maxold
          pijold(ij) = pdiat(ij)
        enddo
        maxnew = ijmat(ip1(natoms+1))-1
        do ij=1,maxnew
          pdiat(ij) = 0.0d0
        enddo
!       copy old blocks to new blocks for atom pairs that the old and
!       new pairlists have in common.
        do iatm=2,natoms
!         range of elements to search in new pairlist:
          j1new = ip1(iatm)
          j2new = ip1(iatm+1)-1
          if(j2new.ge.j1new) then
!           range of elments to search in old pairlist:
            j1old = ip1old(iatm)
            j2old = ip1old(iatm+1)-1
            if(j2old.ge.j1old) then
              do jnew=j1new,j2new
!               go to the next value of iatm if we've looked at all of
!               the iatm pairs in the old pairlist.
!               if(j1old.gt.j2old) go to 400
                jatm = ipair(jnew)
!               if jatm is a higher atom number than the last atom in this
!               block of the old pairlist, then don't bother searching for
!               any more matches with iatm,jatm pairs for the current value
!               of iatm.
                if(jatm.gt.ipold(j2old)) go to 400
                jmatch = 0
                do jold=j1old,j2old
                  if(ipold(jold).eq.jatm)then
!                   element jold of the old pairlist matches element jnew of
!                   the new pairlist.  copy the corresponding block of the
!                   density matrix.
                    jmatch = jold
                    ij1new = ijmat(jnew)
                    ij2new = ijmat(jnew+1)-1
                    ijshft = ijold(jold)-ij1new
                    do ij=ij1new,ij2new
                      pdiat(ij) = pijold(ij+ijshft)
                    enddo
                    go to 310
                  endif
                enddo
!               update pointer for search in old pairlist if we found a match.
 310            if(jmatch.ne.0) j1old = jmatch + 1
              enddo
            endif
          endif
 400      continue
        enddo
      endif
      first = .false.
      return
      end
