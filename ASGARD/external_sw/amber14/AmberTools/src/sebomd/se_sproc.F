!------------------------------------------------------------------------
!     subsystem processing routines.
!
!     sproc1 should be called whenever a new subsystem atom list is
!     created.  it should be called before attempting to create or
!     update the pairlist.  sproc1 requires the subsystem information
!     nsub, iatom1, iatoms, iabuff.
!
!     sproc2 should be called after each pairlist update.  it requires
!     pairlist arrays ip1 and ipair and the subsystem arrays isub1
!     and isubs.
!
!
!******************************************************************************
!
!
      subroutine se_sproc1(ierror)
!
!     does a variety of subsystem setup work prior to creating
!     pairlist.
!
!     bug fix: 6/29/98 -- s. dixon.  when dummy atoms appear after the
!     first atom, the is1atm pointer wasn't being updated properly.
!     this occurred in the "do 500 iatm=1,natoms loop".  now fixed.
!    7/26/99 fix also applies to sparkles
!
      implicit none
#include "sebomd.dim"
#include "se_fulsys.h"
#include "se_subsys.h"
#include "se_elemnt.h"
#include "se_parproc.h"
#include "se_sublst.h"
      integer ierror

      integer iorbk, nsorb, k, i1, i2, norbs, iorbi, i, iatm
      integer norbsi, ipt, ksub, j1, j2, j, jatm, nsubs
!
!
!     assign atomic orbital pointers for each subsystem:
!
!     iorbpt(k) = global pointer to first atomic orbital in subsystem k.
!     iorb1(i) = local pointer in a given subsystem for the first
!                atomic orbital on atom iatom(i).
!
!     set error flag if global or local storage limits are exceeded.
!
      iorbk = 1
      nsorb = 0
      do k=1,nsub
        iorbpt(k) = iorbk
        i1 = iatom1(k)
        i2 = iatom1(k+1)-1
        norbs = 0
        iorbi = 1
        do i=i1,i2
          iatm = iatoms(i)
          norbsi = natorb(iatnum(iatm))
          norbs = norbs + norbsi
          iorb1(i) = iorbi
          iorbi = iorbi + norbsi
        enddo
        nsorb = max(nsorb,norbs)
        if(norbs.gt.msorb)then
          write(iout,'(" SUBSYSTEM ",I4," HAS ",I4," ATOMIC ORBITALS,",
     .                 " WHICH EXCEEDS MAXIMUM ALLOWED")') k,norbs
        endif
        iorbk = iorbk + norbs
      enddo
      iorbpt(nsub+1) = iorbk
      if(nsorb.gt.msorb)then
        ierror = 1
        write(iout,'(/" LOCAL STORAGE FOR SUBSYSTEM ATOMIC ORBITALS",
     .                " EXCEEDED --"/" INCREASE MSORB PARAMETER IN",
     .                " sebomd.dim TO AT LEAST ",I4)') nsorb
      endif
      if(iorbk.gt.msval)then
        ierror = 1
        write(iout,'(/" GLOBAL STORAGE FOR SUBSYSTEM EIGENVALUES",
     .                " EXCEEDED --"/" INCREASE MSVAL PARAMETER IN",
     .                " sebomd.dim TO AT LEAST ",I7)') iorbk
      endif
      if(ierror.ne.0) return
!
!     create lists of subsystems isubs in which each atom appears as
!     a non-buffer atom.
!
      ipt = 1
      do iatm=1,natoms
!
!       isub1(iatm) is a pointer to the position in isubs that
!       contains the first subsystem in which iatm appears as a
!       non-buffer atom.
!
        isub1(iatm) = ipt
!
!       don't bother to look for subsystems if iatm is not a real atom.
!
        if(natorb(iatnum(iatm)).ne.0) then
          do ksub=1,nsub
            j1 = iatom1(ksub)
            j2 = iatom1(ksub+1)-1
!        
!           see if iatm coincides with any non-buffer atoms in ksub.
!        
            do j=j1,j2
              jatm = iatoms(j)
!        
!             if jatm is a higher atom number than iatm, we don't need to
!             continue the search in ksub because the atom list is sorted.
!        
              if(jatm.gt.iatm) exit ! exit loop
!        
!             skip current atom if it's a buffer in ksub.
!        
              if(iabuff(j).eq.0) then
!        
                if(jatm.eq.iatm)then
!              
!                 store ksub in the list for iatm and check the next subsystem.
!              
                  isubs(ipt) = ksub
                  ipt = ipt + 1
                  exit ! exit loop
                endif
              endif
            enddo
          enddo
        endif
      enddo
!
!     store overflow pointer.
!
      isub1(natoms+1) = ipt
!
!     now check to make sure that each real atom appears at least once
!     as a non-buffer atom.
!
      do i=1,natoms
        if(natorb(iatnum(i)).ne.0) then
          nsubs = isub1(i+1)-isub1(i)
          if(nsubs.le.0)then
            ierror = 1
            write(iout,'(" ATOM ",I5," DOES NOT APPEAR AS A NON-BUFFER",
     .                   " ATOM IN ANY SUBSYSTEM")') i
          endif
        endif
      enddo
!
!     for each atom, generate the list of subsystems iasubs in which
!     the atom appears.  iasubs differs from isubs in that iasubs
!     does not distinguish whether or not the atom has buffer status.
!
      ipt = 1
      do iatm=1,natoms
!
!       pointer for first subsystem belonging to iatm:
!
        is1atm(iatm) = ipt
!
!       skip iatm if it's not a real atom.  this is the bug fix.  the go to
!       statement appeared before the pointer assignment in the previous
!       version, and it caused an error with dummy atoms  -- s. dixon.
!
        if(natorb(iatnum(iatm)).ne.0) then
!
          do ksub=1,nsub
            j1 = iatom1(ksub)
            j2 = iatom1(ksub+1)-1
!        
!           see if iatm matches any atoms in ksub.
!        
            do j=j1,j2
              jatm = iatoms(j)
!        
!             because the subsystem atom lists are sorted, we can discontinue
!             search in ksub if we've 'passed up' iatm.
!        
              if(jatm.gt.iatm) exit ! exit current loop
!        
              if(jatm.eq.iatm)then
!        
!               iatm is in ksub.  store ksub in its list.
!        
                iasubs(ipt) = ksub
                ipt = ipt + 1
                exit ! exit loop
              endif
            enddo
          enddo
        endif
      enddo
!
!     overflow pointer:
!
      is1atm(natoms+1) = ipt
      return
      end
!
!
!
      subroutine se_sproc2
!
!     for each pair of atoms (iatm,jatm) in the pairlist, this routine
!     assigns the number of subsystems nshare in which both atoms appear:
!
!     k = pairlist entry for (iatm,jatm).
!     nshare(1,k) = number of times (iatm,jatm) appear with at least
!                   one atom having non-buffer status, and the other
!                   atom having non-buffer status or inner buffer layer
!                   status.
!     nshare(2,k) = number of times (iatm,jatm) appear with any status.
!
!
!     nshare(1,k) keeps track of how many times a pair of atoms contribute
!     to the density matrix, while nshare(2,k) keeps track of contributions
!     to the 1-electron and fock matrices.
!
!     should be called after each pairlist update, but only after
!     calling se_ijmake.
!
      implicit none
#include "sebomd.dim"
#include "se_fulsys.h"
#include "se_subsys.h"
#include "se_elemnt.h"
#include "se_parproc.h"
#include "se_sublst.h"
#include "se_bpairs.h"

      integer npairs, i, k, natms, iatm0, i0i, iatm, ibuffi, i0j, jatm
      integer ijaddr, ibuffj, ijsum, j
!
      if(natoms.eq.1) return
!
      npairs = ip1(natoms+1)-1
      do i=1,npairs
        nshare(1,i) = 0
        nshare(2,i) = 0
      enddo
!
!     loop over each pair of atoms in each subsystem and accumulate
!     nshare based on buffer status.
!
      do k=1,nsub
        natms = iatom1(k+1) - iatom1(k)
        if(natms.ne.1) then
          iatm0 = iatom1(k)-1
          do i=2,natms
            i0i = iatm0+i
            iatm = iatoms(i0i)
            ibuffi = iabuff(i0i)
            do j=1,i-1
              i0j = iatm0+j
              jatm = iatoms(i0j)
!             skip the iatm,jatm pair if it hasn't been stored as a
!             bonded pair.
              call se_ijfind(npairs,iatm,jatm,ijaddr)
              if(ijaddr.ne.0) then
                ibuffj = iabuff(i0j)
                ijsum = ibuffi + ibuffj
                if(ijsum.le.1) nshare(1,ijaddr) = nshare(1,ijaddr) + 1
                nshare(2,ijaddr) = nshare(2,ijaddr) + 1
              endif
            enddo
          enddo
        endif
      enddo
      return
      end
