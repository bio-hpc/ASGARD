      subroutine se_resclustsub(ipoint, natot, nbig, ibig, ierror)
!
!     routine to do cluster subsetting.
!
!     March 31 1998, Arjan van der Vaart: added multiple core feature.
!     the core will consist of ncore(1) residues selected from the
!     group of residues icorel(1)..icorel(icorel1(2)-1), ncore(2) residues
!     selected form the group of residues 
!     icorel(icorel1(2))..icorel(icorel1(3)-1) etc. Buffers will always
!     be build from the entire set (1..nres) of residues.
!
      implicit none
#include "sebomd.dim"
#include "se_logic.h"
#include "se_xyzcrd.h"
#include "se_pbcbox.h"
#include "se_fulsys.h"
#include "se_clst.h"
#include "se_subsys.h"
#include "se_endecomp.h"
#include "se_parproc.h"
#include "se_subpar.h"
#include "se_elemnt.h"
#include "se_resdef.h"
!
!
!     local:
!
      integer ipoint, natot, nbig, ibig, ierror
      integer ircore(maxres),ir1(maxres),
     .          ilook(maxres),icatms(maxatm),ibatms(2,maxatm),
     &     ibuflist(maxres), itmp(maxres)
      logical incore(maxres),insub(maxres),ingsub(maxres),harsh, 
     &     addineighb
      integer ineighb, nlarge, isubmx, ncsub, ncatot, nresn, nacore
      integer istore, ir, ires, ia1, ia2, iatm
      integer jr, jres, ja1, ja2, jatm
      integer i, i1, icr, i0, iadd, is, ipass, icore, ileft, iright
      integer j, j1, j2
      integer nbuff1, nbuff2, nval, nused, iseed, iseed0, iseed1
      integer nadded, k, nat, jrmin, isub, natavg, ntmp, j0, istart
      double precision rcut, rbuff1, rbuff2, rmin, rsqr, rij
      double precision xi, yi, zi
      double precision xj, yj, zj
!


!c---------------------------------------------------------------------cc
      addineighb = .false.

      neighbor(1) = 0
      ineighb = 1


      if(prtsub) write(iout,'(" subsystem definitions:")')
!
!     rcut = square of actual cutoff in angstroms.
      rcut = 49.0d0

      rbuff1 = dcbuff1**2
      rbuff2 = (dcbuff1 + dcbuff2)**2
      nlarge = 0
      isubmx = 0
      ncsub = 0
      ncatot = 0
      isubend(0) = 0

!     make cores for all different values of ncore
      do 5000 icr=1,ncores

         if (addineighb) then
            ineighb = ineighb + 1
            addineighb = .false.
         endif
!     .  i0 is the offset for the residue list and nresn is the
!     .  number of residues with ncore=ncore(icr)
         i0 = icorel1(icr)-1
         nresn = icorel1(icr+1)-icorel1(icr)


         if(ncore(icr).ne.1.and.ncore(icr).ne.nresn)then
!
!     .     construct a residue-based pairlist and store the mininum
!     .     inter-residue interatomic distance in each case.  store
!     .     pairs that are within 5.0 angstroms.  if there are isolated
!     .     residues that are located a great distance from all other
!     .     residues, then it may be necessary to increase this cutoff
!     .     in order to get successful subsetting.
!     
            istore = 1
            do 80 ir=1,nresn
               ires = icorel(ir+i0)
               ir1(ir) = istore
               ia1 = irpnt(ires)
               ia2 = irpnt(ires+1)-1
               do 70 jr=1,nresn
                  jres = icorel(jr+i0)
                  if(ires.eq.jres) go to 70
                  ja1 = irpnt(jres)
                  ja2 = irpnt(jres+1)-1
                  rmin = 1.0d10
                  do 60 iatm=ia1,ia2
                     if(natorb(iatnum(iatm)).eq.0) go to 60
                     do 50 jatm=ja1,ja2
                        if(natorb(iatnum(jatm)).eq.0) go to 50
                        xj = gcres(1,jatm)
                        yj = gcres(2,jatm)
                        zj = gcres(3,jatm)
                        rsqr = (gcres(1,iatm)-xj)**2
     .                       + (gcres(2,iatm)-yj)**2
     .                       + (gcres(3,iatm)-zj)**2
                        rmin = min(rmin,rsqr)
 50                  continue
 60               continue
                  if(rmin.lt.rcut)then
                     irpair(istore) = jr
                     rdist(istore) = rmin
                     istore = istore + 1
                     if(istore.gt.mxpres)then
                        ierror = 1
                        write(iout,'(/" maximum storage for",
     .                       " residue-based",
     .                       " pairlist exceeded",/" -- increase",
     .                       " mxpres parameter in sebomd.dim and",
     .                       " recompile")')
                        return
                     endif
                  endif
 70            continue
!     
!     .        sort pairlist entries for ires on the basis of increasing
!     .        distance.
!     
               i1 = ir1(ir)
               nval = istore - i1
               call se_rsort(nval,rdist(i1),irpair(i1))
 80         continue
            ir1(nresn+1) = istore
         endif
!
!     .  build each subsystem by generating a core collection of residues,
!     .  then extending out using the buffer thickness.
!     
!     .  pbcs will be used (if requested) to generate buffer regions.
!     
!     
         nused = 0
         iseed1 = 1
         do 90 ires=1,nresn
            incore(ires) = .false.
 90      continue

!     
!     .  nused = count of residues that appear in any core.
!     
!     .  iseed1 = starting point for search in residue list to find
!     .  seed residue to start next core.
!     
!     .  incore = flag for residues that appear in any core.
!     
!     .  here all arrays and iseed is used locally i.e. filled in from
!     .  1 till nresn. this means that e.g. incore(i) means that residue
!     .  icorel(i+i0) is in the core.
!     .  ingsub(i) is the only global array i.e. the indices go from 1 
!     .  till nres. this means that ingsub(i) means that residue i is 
!     .  included in the current subsystem.
!     
         do 500 is=1,nresn
            isub = icorel(is+i0)
            nsub = nsub + 1
            ncsub = ncsub + 1
            harsh = .true.
!     
!     .     initially, a harsh criterion will be used for addition of
!     .     residues to core once the seed residue has been selected.
!     .     the harsh criterion is that a residue will not be added
!     .     to the core of isub if it appears in the core of any other
!     .     subsystem.  if residues become "isolated" or if ncore does
!     .     not divide evenly into nres, then it may become necessary to
!     .     relax criterion and allow previously used residues to be added.
!     
 100        nadded = 0
            iseed0 = iseed1
!     
!     .     nadded = the number of residues in the core of isub.
!     
!     .     iseed0 will be used to reset iseed1 if the search has to
!     .     be repeated because the harsh criterion was too harsh.
!     
            do 110 ires=1,nresn
               insub(ires) = .false.
 110        continue
            do ires=1,nres
               ingsub(ires) = .false.
            enddo
!     
!     .     insub flags residues that have been added to the core of isub.
!     
            if(ncore(icr).eq.1)then
!     
!     .     only one residue per core -- isub will just correspond to
!     .     the core residue.
!     
               nused = nused + 1
               nadded = nadded + 1
               ircore(1) = is
               insub(is) = .true.
               ingsub(isub) = .true.
            elseif(ncore(icr).eq.nresn)then
!     
!     .        only one subsystem -- all residues.
!     
               nused = nresn
               nadded = nresn
               do 120 ires=1,nresn
                  ircore(ires) = ires
                  ingsub(icorel(ires+i0)) = .true.
 120           continue
            else
!     
!     .        search residues that do not already appear in any core list to
!     .        find a seed as a starting point for the next core.
!     
               do 140 iseed=iseed1,nresn
                  if(.not.incore(iseed))then
                     incore(iseed) = .true.
                     insub(iseed) = .true.
                     ingsub(icorel(iseed+i0)) = .true.
                     ilook(iseed) = ir1(iseed)
!     
!     .              ilook(ires) = starting point for search in ires pairlist
!     .              to locate the next "close" residue.
!     
                     nused = nused + 1
                     nadded = 1
                     ircore(1) = iseed
!     
!     .              ircore(i) = ith core residue in isub.
!     
                     go to 150
                  endif
 140           continue
 150           iseed1 = iseed + 1
!     
!     .        now make ncore-1 passes through the pairlist of the core
!     .        residues, bringing in the closest residue jrmin at each pass.
!     
               do 200 ipass=2,ncore(icr)
                  rmin = 1.0d10
                  jrmin = 0
                  do 180 icore=1,ipass-1
                     ir = ircore(icore)
                     ires = icorel(ir+i0)
!     
!     .              skip this core residue if we've reached the end of
!     .              its pairlist.
!     
 160                 j1 = ilook(ir)
                     j2 = ir1(ir+1)
                     if(j1.eq.j2) go to 180
                     ilook(ir) = j1+1
!     
                     jres = irpair(j1)
!     
!     .              if jres has already been used (whichever criterion),
!     .              then go back and check the next entry in the ires
!     .              pairlist.
!     
                     if(harsh)then
                        if(incore(jres)) go to 160
                     else
                        if(insub(jres)) go to 160
                     endif
!     
                     rij = rdist(j1)
                     if(rij.lt.rmin)then
                        rmin = rij
                        jrmin = jres
                     endif
 180              continue
                  if(jrmin.eq.0.and..not.harsh)then
!     
!     .              return with an error if the relaxed criterion is
!     .              in effect and we still can't find a close residue
!     .              in the core pairlist.  it is probably a result of
!     .              the current cores being far away from everything else.
!     .              it should be possible to resolve this problem by
!     .              increasing rcut.
!     
                     ierror = 1
                     write(iout,'(/" cannot build core for subsystem ",
     .                    i4,/" -- try increasing cutoff value in",
     .                    " clustsub")') isub
                     return
                  elseif(jrmin.eq.0)then
!     
!     .              couldn't find a close residue using harsh criterion.
!     .              switch to relaxed criterion and try again.
!     
                     harsh = .false.
                     do 190 i=1,ipass-1
                        ires = ircore(i)
                        incore(ires) = .false.
 190                 continue
                     nused = nused - nadded
                     iseed1 = iseed0
                     go to 100
                  endif
!     
!     .           add closest residue to core.
!     
                  ircore(ipass) = jrmin
                  ilook(jrmin) = ir1(jrmin)
                  insub(jrmin) = .true.
                  ingsub(icorel(jrmin+i0)) = .true.
                  if(harsh.and.incore(jrmin)) go to 200
                  if(.not.incore(jrmin)) nused = nused + 1
                  nadded = nadded + 1
                  incore(jrmin) = .true.
 200           continue
            endif
!     
!     .     core is complete.  store core atoms in icatms.
!     
            nacore = 0
            do 220 iadd=1,nadded
               ir = ircore(iadd)
               ires = icorel(ir+i0)
               ia1 = irpnt(ires)
               ia2 = irpnt(ires+1)-1
               do 210 iatm=ia1,ia2
                  if(natorb(iatnum(iatm)).eq.0) go to 210
                  nacore = nacore + 1
                  icatms(nacore) = iatm
 210           continue
 220        continue
!     
!     .     create two sets of buffer atoms for this subsystem.  the
!     .     first set will contain all atoms that are within sqrt(rbuff1)
!     .     of the core.  the second set will contain atoms that are
!     .     beyond sqrt(rbuff1) but within sqrt(rbuff2) of core.
!
!     .     if noovrlp is .false. then
!     .     note that atoms for the buffers are extracted from the global
!     .     set of residues, i.e. from residues 1..nres (and not only
!     .     from residues icorel(1+i0)..icorel(nresn+i0))
!
!     .     if noovrlp is .true. then 
!     .     buffers are created from "own" local list
!     
            nbuff1 = 0
            nbuff2 = 0
            if(ncore(icr).eq.nres) go to 400

            if (.not.noovrlp) then
               do 300 ires=1,nres

                  if(ingsub(ires)) go to 300
                  
                  ia1 = irpnt(ires)
                  ia2 = irpnt(ires+1)-1
                  do 280 iatm=ia1,ia2
!     
!     .              iatm is the candidate atom for addition to buffer lists.
!     
                     if(natorb(iatnum(iatm)).eq.0) go to 280
                     xi = gcres(1,iatm)
                     yi = gcres(2,iatm)
                     zi = gcres(3,iatm)
!     
!     .              get minimum distance between iatm and core atoms.
!     
                     rmin = 1.0d10
                     do 250 j=1,nacore
                        jatm = icatms(j)
                        if(pbc)then
                           call se_pbcgcres(iatm,jatm,xj,yj,zj)
                        else
                           xj = gcres(1,jatm)
                           yj = gcres(2,jatm)
                           zj = gcres(3,jatm)
                        endif
                        rij = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
                        rmin = min(rmin,rij)
 250                 continue
                     if(rmin.le.rbuff1)then
                        nbuff1 = nbuff1 + 1
                        ibatms(1,nbuff1) = iatm
                     elseif(rmin.le.rbuff2)then
                        nbuff2 = nbuff2 + 1
                        ibatms(2,nbuff2) = iatm
                     endif
 280              continue
 300           continue
            else
!              if (neighbor(ineighb).eq.icr) then
!                 ir = ircore(1)
!                 ires = icorel(i0+ir)
!                 ileft = ires - neighborn(ineighb)
!                 ileft = max(1,ileft)
!                 ntmp = 0
!                 do j=ileft,ires-1
!                    ntmp = ntmp + 1
!                    itmp(ntmp) = j
!                 enddo
!                 iright = ires + neighborn(ineighb)
!                 iright = min(nres,iright)
!                 do j=ires+1,iright
!                    ntmp = ntmp + 1
!                    itmp(ntmp) = j
!                 enddo
!                 j0 = 0
!                 if (ntmp.gt.0) then
!                    k = 1
!                    do j=1,nresn
!301                    if (icorel(j+i0).gt.itmp(k)) then
!                          k = k + 1
!                          if (k.gt.ntmp) goto 302
!                          goto 301
!                       elseif (icorel(j+i0).eq.itmp(k)) then
!                          j0 = j0 + 1
!                          ibuflist(j0) = itmp(k)
!                          k = k + 1
!                          if (k.gt.ntmp) goto 302
!                       endif
!                    enddo
!                 endif
!302              addineighb = .true.
!              else
                  j0 = nresn
                  do j=1,nresn
                     ibuflist(j) = icorel(j+i0)
                  enddo
!              endif

!     the following commented code is without neighbors
!               do 390 ir=1,nresn
!                 ires = icorel(ir+i0)

               do 390 ir=1,j0
                  ires = ibuflist(ir)

                  if(ingsub(ires)) go to 390
                  
                  ia1 = irpnt(ires)
                  ia2 = irpnt(ires+1)-1
                  do 380 iatm=ia1,ia2
!     
!     .              iatm is the candidate atom for addition to buffer lists.
!     
                     if(natorb(iatnum(iatm)).eq.0) go to 380
                     xi = gcres(1,iatm)
                     yi = gcres(2,iatm)
                     zi = gcres(3,iatm)
!     
!     .              get minimum distance between iatm and core atoms.
!     
                     rmin = 1.0d10
                     do 350 j=1,nacore
                        jatm = icatms(j)
                        if(pbc)then
                           call se_pbcgcres(iatm,jatm,xj,yj,zj)
                        else
                           xj = gcres(1,jatm)
                           yj = gcres(2,jatm)
                           zj = gcres(3,jatm)
                        endif
                        rij = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
                        rmin = min(rmin,rij)
 350                 continue
                     if(rmin.le.rbuff1)then
                        nbuff1 = nbuff1 + 1
                        ibatms(1,nbuff1) = iatm
                     elseif(rmin.le.rbuff2)then
                        nbuff2 = nbuff2 + 1
                        ibatms(2,nbuff2) = iatm
                     endif
 380              continue
 390           continue
               
            endif

!     
!     .     write out latest subsystem information if requested.
!     
 400        if(prtsub)then
               write(iout,'(/" subsystem ",i4,
     .              /" --------------",
     .              /" core atoms:")') nsub
               write(iout,'(16i5)') (icatms(i),i=1,nacore)
               if(nbuff1.gt.0)then
                  write(iout,'(" inner buffer layer:")')
                  write(iout,'(16i5)') (ibatms(1,i),i=1,nbuff1)
               endif
               if(nbuff2.gt.0)then
                  write(iout,'(" outer buffer layer:")')
                  write(iout,'(16i5)') (ibatms(2,i),i=1,nbuff2)
               endif
               write(iout,*)
            endif
!     
            nat = nacore + nbuff1 + nbuff2
!     
!     .     make sure we have enough space to store the nat atoms in
!     .     this subsystem.
!     
            if((ipoint+nat).ge.mslist)then
               ierror = 1
               write(iout,'(/" atom list storage limit reached ",
     .              "at subsystem",i4,
     .              /" -- increase mslist parameter in dicon.dim",
     .              " and recompile")') nsub
               return
            endif
!     
!     .     store core and buffer atoms in subsystem atom list.  assign
!     .     buffer status iabuff as we go
!     
            iatom1(nsub) = ipoint
            do 420 i=1,nacore
               iatoms(ipoint) = icatms(i)
               iabuff(ipoint) = 0
               ipoint = ipoint + 1
 420        continue
            if(nbuff1.gt.0)then
               do 440 i=1,nbuff1
                  iatoms(ipoint) = ibatms(1,i)
                  iabuff(ipoint) = 1
                  ipoint = ipoint + 1
 440           continue
            endif
            if(nbuff2.gt.0)then
               do 460 i=1,nbuff2
                  iatoms(ipoint) = ibatms(2,i)
                  iabuff(ipoint) = 2
                  ipoint = ipoint + 1
 460           continue
            endif
!     
!     .     sort atom list and carry along buffer status.
!     
            istart = iatom1(nsub)
            call se_bsort(nat,iatoms(istart),iabuff(istart))
!     
            natot = natot + nat
            ncatot = ncatot + nat
            if(nat.gt.nlarge)then
               nlarge = nat
               isubmx = nsub
               if (nat.gt.nbig) then
                  nbig = nat
                  ibig = nsub
               endif
            endif
!     
            if(nused.eq.nresn)then
!     
!     .        all residues have been used -- we are done subsetting.
!     
               goto 4900
            endif

 500     continue

 4900    isubend(icr) = nsub

 5000 continue

!     overflow pointer:
!     
      iatom1(nsub+1) = ipoint
!
      natavg = ncatot/ncsub
!-rdc       write(iout,'(//" summary of cluster based subsetting:",
!-rdc      .     /" -----------------------------------",
!-rdc      .     //" number of subsystems   =  ",i4,
!-rdc      .     /" average subsystem size =  ",i4,
!-rdc      .     /" maximum subsystem size =  ",i4,
!-rdc      .     /" largest subsystem      = #",i4)')
!-rdc      .     ncsub,natavg,nlarge,isubmx
      
      return

      end
