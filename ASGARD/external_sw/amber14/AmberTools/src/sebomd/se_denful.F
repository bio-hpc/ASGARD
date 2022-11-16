      subroutine se_denful(norbs,evec1)
!
!     computes the global density matrix (pdiag,pdiat) for the full system
!     in a calculation that does not involve separate subsystems.
!
      implicit none
#include "sebomd.dim"
#include "se_bpairs.h"
#include "se_electr1.h"
#include "se_elemnt.h"
#include "se_fulsys.h"
#include "se_subsys.h"
#include "se_glbmat1.h"
#include "se_glbmat2.h"
      integer norbs
      double precision evec1(norbs,*)

      integer nocc, nocc1, nocc2
      integer nhalf, natms, npairs

      integer iimax, ipmax, ijmax
      integer ijaddr, ijpmat, idiag, ijp, iatm1

      integer i, norbsi, ii, iatm, i0, ip, ievec1
      integer j, norbsj, ij, jatm, j0, jp, jevec1
      integer l

      double precision psign
      double precision pij

      nocc = nelec/2
      nhalf = norbs/2
      if(nocc.gt.nhalf)then
!
!       density matrix will be formed by subtracting electrons from
!       a saturated system.  this is usually the case.
!
        nocc1 = nocc + 1
        nocc2 = norbs
        psign = -2.0d0
      else
!
!       density matrix will be formed by adding electrons to an empty system.
!
        nocc1 = 1
        nocc2 = nocc
        psign = 2.0
      endif
      iimax = iimat(natoms+1)-1
      ipmax = ip1(natoms+1)
      ijmax = ijmat(ipmax)-1
      do i=1,iimax
        pdiag(i) = 0.0d0
      end do
      do ij=1,ijmax
        pdiat(ij) = 0.0d0
      end do
      if(psign.lt.0.0d0)then
!
!       assign a 'saturated' density matrix with each diagonal=2.0.
!
        ii = 0
        do i=1,natoms
          norbsi = natorb(iatnum(i))
!         skip if dummy or sparkle.
          if(norbsi.ne.0) then
            do idiag=1,norbsi
              ii = ii + idiag
              pdiag(ii) = 2.0d0
            end do
          endif
        end do
      endif
      natms = iatom1(2)-1
!
!     diagonal blocks first.
!
      do i=1,natms
        iatm = iatoms(i)
        norbsi = natorb(iatnum(iatm))
        ijp = iimat(iatm)
        i0 = iorb1(i)-1
        j0 = i0
        do ip=1,norbsi
          ievec1 = i0+ip
          do jp=1,ip
            pij = 0.0d0
            jevec1 = j0+jp
            do l=nocc1,nocc2
              pij = pij + evec1(ievec1,l)*evec1(jevec1,l)
            end do
            pdiag(ijp) = pdiag(ijp) + psign*pij
            ijp = ijp + 1
          end do
        end do
      end do
!
!     now off-diagonal blocks.
!
      if(natms.eq.1) return
      npairs = ip1(natoms+1)-1
      iatm1 = iatoms(iatom1(1))
      do i=2,natms
        iatm = iatoms(i)
        norbsi = natorb(iatnum(iatm))
        i0 = iorb1(i)-1
        do j=1,i-1
          jatm = iatoms(j)
          norbsj = natorb(iatnum(jatm))
!
!         find postion in pairlist of (iatm,jatm) pair.
!
          call se_ijfind(npairs,iatm,jatm,ijaddr)
          if(ijaddr.ne.0) then
            ijpmat = ijmat(ijaddr)
            j0 = iorb1(j)-1
            do ip=1,norbsi
              ievec1 = i0+ip
              do jp=1,norbsj
                pij = 0.0d0
                jevec1 = j0+jp
                do l=nocc1,nocc2
                  pij = pij + evec1(ievec1,l)*evec1(jevec1,l)
                end do
                pdiat(ijpmat) = pdiat(ijpmat) + psign*pij
                ijpmat = ijpmat + 1
              end do
            end do
          endif
        end do
      end do
      return
      end
