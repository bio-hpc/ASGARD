      subroutine se_bpair(ierror)

!     bonded atom pairlist generation routine.  two atoms are stored
!     in the bonded atom pairlist if they ever appear together in the
!     same subsystem and they are separated by less than the user
!     parameter cutbond.
!
!     for i=2,3,...,natoms, the atoms that "bond" with i and which have
!     lower atom numbers than i are:
!
!         ipair(k), for k=ip1(i),ip1(i)+1,...,ip1(i+1)-1
!
!     if the maximum allowed number of pairs is exceeded, ierror is
!     returned with the atom number i at which the maximum was reached.
!     this routine should be called each time a new set of subsystems
!     is defined.
!
!     note that the common block /pairij/ is used temporarily here to
!     store the old pairlist and pointers to diatomic blocks in the
!     old density matrix.  this information will be used in initp to
!     construct a new density matrix from the old one.
!
      implicit none
#include "sebomd.dim"
#include "se_bpairs.h"
#include "se_fulsys.h"
#include "se_glbmat2.h"
#include "se_pairij.h"
#include "se_sublst.h"

      integer ierror
     
!     local:
      integer iatm
      integer ijmax, npairs
      integer i, i1, i2, ilist, isub
      integer j, j1, j2, jlist, jsub
      logical first
      data first /.true./
      save first

      if(first)then
        first = .false.
      else
!       temporarily store old pairlist and diatomic block pointers.
        do iatm=2,natoms
          j1 = ip1(iatm)
          j2 = ip1(iatm+1)-1
          ip1old(iatm) = j1
          if(j2.ge.j1) then
            do j=j1,j2
              ipold(j) = ipair(j)
              ijold(j) = ijmat(j)
            enddo
          endif
        enddo
        ip1old(natoms+1) = ip1(natoms+1)
        ijmax = ijmat(ip1(natoms+1))-1
        ijold(ijmax) = ijmat(ijmax)
      endif

      ierror = 0
      npairs = 0
      do i=2,natoms
        ip1(i) = npairs + 1
        if(iatnum(i).ne.0) then
!         pointers for subsystems belonging to atom i:
          i1 = is1atm(i)
          i2 = is1atm(i+1)-1
          if(i2.ge.i1) then
!           locate atoms < i that share a subsystem with atom i.
            do j=1,i-1
              if(iatnum(j).ne.0) then
!               pointers for subsystems belonging to atom j:
                j1 = is1atm(j)
                j2 = is1atm(j+1)-1
!               skip i,j pair if subsystem lists have no chance of overlap.
                if(j2.ge.j1) then
                  if(iasubs(j2).ge.iasubs(i1).and.
     .               iasubs(i2).ge.iasubs(j1)) then
!                   loop over i-j subsystems and look for one in common.
                    do ilist=i1,i2
                      isub = iasubs(ilist)
                      do jlist=j1,j2
                        jsub = iasubs(jlist)
                        if(isub.eq.jsub)then
!                         we have a match -- store the i,j pair.
                          npairs = npairs + 1
                          if(npairs.gt.mbpair)then
                            ierror = i
                            return
                          endif
                          ipair(npairs) = j
                          go to 90
                        endif
                      enddo
                    enddo
                  endif
                endif
              endif
 90           continue
            enddo
          endif
        endif
      enddo
!     store extra pointer for convenience in looping.
      ip1(natoms+1) = npairs + 1
      return
      end
