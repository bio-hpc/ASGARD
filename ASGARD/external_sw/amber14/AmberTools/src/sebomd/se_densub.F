      subroutine se_densub(k,norbsk,evec1)
!
!     adds in density matrix contribution from subsystem k.  it is
!     assumed that it will be faster to subtract electrons from a
!     saturated density matrix, which is just 2.0 times the identity
!     matrix.
!
      implicit none
#include "sebomd.dim"
#include "se_bal_mpi.h"
#include "se_bpairs.h"
#include "se_elemnt.h"
#include "se_electr.h"
#include "se_fulsys.h"
#include "se_subsys.h"
#include "se_glbmat1.h"
#include "se_glbmat2.h"
#include "se_parproc.h"
#ifdef MPI
#include "mpif.h"
#endif
      integer k, norbsk
      double precision evec1(norbsk,*)
!
!     local:
!
      double precision hole(msorb)
      logical dodens
      integer ipmax, ijmax, nlower, natms, npairs, nsubsi
      integer i, norbsi, ii, ij, i1, i2, ip, ijp, i0, i0i, iatm0, iatm
      integer j, norbsj, jj, jp, i0j, ijsum, jatm, ijaddr, j0
      integer l, l0, ievec1, jevec1, ijpmat
      double precision dij, pij, dii

!
!     initialize a saturated density matrix if this is the first
!     subsystem.
!
      dodens = .true.
#ifdef MPI 
!     jv first time through each pe must init density matrix
      if (k .eq. my_subs(1) ) then
#else
      if (k.eq.1) then
#endif
        ipmax = ip1(natoms+1)
        ijmax = ijmat(ipmax)-1
        do i=1,natoms
          norbsi = natorb(iatnum(i))
          if(norbsi.ne.0) then
            ij = iimat(i)
            do ii=1,norbsi
              do jj=1,ii
                pdiag(ij) = 0.0d0
!               jv only master inits diag so that we don't sum from all pe's
#ifdef MPI 
                if((jj.eq.ii) .and. (myid .eq. 0))
#else      
                if((jj.eq.ii))
#endif     
     $          pdiag(ij) = 2.0d0
                ij = ij+1
              enddo
            enddo
          endif
        enddo
        do ij=1,ijmax
          pdiat(ij) = 0.0d0
        enddo
      endif

      if (dodens) then
!       calculate the density matrix,
!       determine lower mo limit for electron subtraction and assign
!       hole occupation numbers from fermi occupation numbers.
!
        nlower = 0
        do i=iorbpt(k),iorbpt(k+1)-1
          nlower = nlower + 1
          if(fermi(i).lt.0.999999d0) go to 100
        enddo
 100    l0 = iorbpt(k)-1
        do l=nlower,norbsk
          hole(l) = fermi(l0+l) - 1.0d0
        enddo

!       diagonal blocks first.
        i1 = iatom1(k)
        i2 = iatom1(k+1)-1
        do i=i1,i2
!         skip density contribution if this is a buffer atom.
          if(iabuff(i).eq.0) then
            iatm = iatoms(i)
            norbsi = natorb(iatnum(iatm))
            nsubsi = isub1(iatm+1) - isub1(iatm)
            dii = 2.0d0/nsubsi
            ijp = iimat(iatm)
            i0 = iorb1(i)-1
            j0 = i0
            do ip=1,norbsi
              ievec1 = i0+ip
              do jp=1,ip
                jevec1 = j0+jp
                pij = 0.0d0
                do l=nlower,norbsk
                  pij = pij + evec1(ievec1,l)*evec1(jevec1,l)*hole(l)
                enddo
                pdiag(ijp) = pdiag(ijp) + pij*dii
                ijp = ijp + 1
              enddo
            enddo
          endif
        enddo

!       now off-diagonal blocks.
        natms = iatom1(k+1) - iatom1(k)
        if(natms.eq.1) return
        npairs = ip1(natoms+1)-1
        iatm0 = iatom1(k)-1
        do i=2,natms
          i0i = iatm0+i
!         no bonding with i0i if it's an outer layer buffer atom.
          if(iabuff(i0i).ne.2) then
            iatm = iatoms(i0i)
            norbsi = natorb(iatnum(iatm))
            i0 = iorb1(i0i)-1
            do j=1,i-1
              i0j = iatm0+j
!             no bonding with i0j if it's an outer layer buffer atom or
!             if both i0j and i0i have any sort of buffer status.
              ijsum = iabuff(i0i) + iabuff(i0j)
              if(ijsum.lt.2) then
                jatm = iatoms(i0j)
!               find postion in pairlist of (iatm,jatm) pair.
                call se_ijfind(npairs,iatm,jatm,ijaddr)
!               skip this pair if it hasn't been stored as a bonded pair.
                if(ijaddr.ne.0) then
                  norbsj = natorb(iatnum(jatm))
                  j0 = iorb1(i0j)-1
                  ijpmat = ijmat(ijaddr)
                  dij = 2.0d0/nshare(1,ijaddr)
                  do ip=1,norbsi
                    ievec1 = i0+ip
                    do jp=1,norbsj
                      pij = 0.0d0
                      jevec1 = j0+jp
                      do l=nlower,norbsk
                        pij = pij + evec1(ievec1,l)
     $                             *evec1(jevec1,l)*hole(l)
                      enddo
                      pdiat(ijpmat) = pdiat(ijpmat) + pij*dij
                      ijpmat = ijpmat + 1
                    enddo
                  enddo
                endif
              endif
            enddo
          endif
        enddo
      endif
      return
      end
