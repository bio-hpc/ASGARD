      subroutine addspace(ipick,curlesid)
      implicit none
c
#include "SIZE.h"
#include "TOP.h"
#include "TTOP.h"
#include "MISC.h"
      integer i1,numcb,numca,numcd,totalex,izero,inon, tempint, chunk
      real*8 velfac,randnum
      logical debug
      logical size_error
c
c variables for sorting the list of excluded atoms- needs to be
c sorted list for each atom, so we will sort the sub-list.
c
       integer j4,incr,temp,bgn,end
c
       debug=.false.
c
c zero the new list of atoms in molecules
c
      do 5 i=1,nspm
              tnsp(i)=0
5     continue

c
c make the actual copies of the currently defined subspace
c ipick(i)= number of copies to make (1= no change)
c
      tnatom=0

      do 100 i=1,natom
c
c set revpoint, which points to the FIRST atom that is a copy
c of the particle. This first particle may already be a LES copy
c from a previous call to this routine- it is only the first
c atom of the current set of copies. the index is the atom
c number before this set of copies was made.
c
c we'll use it when making bond/angle etc lists- it will show us
c for a given non-copied atom, where to start in the new atom
c list to find the copies of the atom.
c
              revpoint(i)=tnatom+1

              if( ipick(i).ne.1.and.mod( ipick(i), npack ) .ne. 0 ) then
                  write(*,*) 'warning:copy number of atom ', i,
     &                       'is not 1 or times of npack'
                  stop
              else
                  chunk = ipick(i) / npack
              end if

              do 110 j=1,ipick(i)
                 tnatom=tnatom+1
c
c set pointer to original atom
c
                 poiatom(tnatom)=i
c
c copy coords
c
                 if (rcrdok.or.rcvdok.or.rcvbok.or.rcbdok) then
                     if( npack .gt. 1 ) then
                         tempint = ( ( j - 1 ) / chunk ) * natom
                         tx(tnatom)=x(tempint+i)
                         ty(tnatom)=y(tempint+i)
                         tz(tnatom)=z(tempint+i)
                     else
                         tx(tnatom)=x(i)
                         ty(tnatom)=y(i)
                         tz(tnatom)=z(i)
                     endif
                 endif
c
                 if (rcvdok.or.rcvbok) then
c
c velocities were read
c
                     velfac=1.0d0
c
c to keep the KE right, we have to scale each of these UP by sqrt(number of 
c copies)- mv^2 is smaller for these by a factor of N for N copies
c due to the mass scaling. this may make them too fast...
c
                     if (.not.allmas.and.ipick(i).gt.1) 
     &               velfac=velfac* sqrt(float(ipick(i)))

c
c now we have a problem- we need the velocities to be different since the 
c coords are the same- or else the copies will not spread apart. we have to
c introduce some randomness, can't be avoided. the smaller we perturb the
c velocities, the longer it will take for the copies to become uncorrelated.
c but larger perts will change the equilibration by more....
c and we also want to preserve the temperature. well,  for big sets
c of copies they will need to re-equilibrate. we'll perturb the TOTAL VELOCITY
c MAGNITUDE of atoms by say 1% up or down... don't mess with the direction of
c the velocity vector.
c
c do this unless user requested NOMODV option
c
                     if (.not.nomodv.and.ipick(i).gt.1) then
c
c get a random number between 0.9 and 1.1
c
                            call amrand(randnum)
                            velfac=velfac*(randnum*0.2+0.9)

                write (6,'(a,f9.5)') 
     &              "Modifying velocities, factor is",velfac
c
c endif for checking nomodv
c
                      endif
c
                      tvx(tnatom)=vx(i)*velfac
                      tvy(tnatom)=vy(i)*velfac
                      tvz(tnatom)=vz(i)*velfac
c
c endif for rcvbok
c
               endif
c
c set original particle number for subsequent picks
c
               torigpt(tnatom)=origpt(i)
c
c track total # copies made
c
               if (j.gt.1) totcop(origpt(i))=totcop(origpt(i))+1
c
c flag the first copy for ease of making les-noles pert topologies
c only flag first copy of a first copy (for multiple level les subspaces)
c
                treal(tnatom)=(lreal(i).and.j.eq.1)
c
c set # les levels, lesid's
c
                tnlev(tnatom)=nlev(i)
                do 125 k=1,nlev(i)
                    tlesid(tnatom,k)=lesid(i,k)
125             continue
c
c add a level if appropriate, set the copy number
c
                if (ipick(i).gt.1) then
                    tnlev(tnatom)=nlev(i)+1
                    tlesid(tnatom,tnlev(tnatom))=curlesid
                endif
c
c set parameters to that of orig. atom, divided by # copies
c
                tchrg(tnatom)=chrg(i)/dble(ipick(i))
c pert stuff
                if (ifpert.gt.0) then
                    tcgper(tnatom)=cgper(i)/dble(ipick(i))
                    tiaper(tnatom)=iaper(i)
                    tiacper(tnatom)=iacper(i)
                    tiacpfac(tnatom)=iacpfac(i)*ipick(i)
                    talmper(tnatom)=almper(i)
                    tismper(tnatom)=ismper(i)
                    tigrper(tnatom)=igrper(i)
                endif
                if (allmas) then
                   tamass(tnatom)=amass(i)
                   if (itimass) ttimass(tnatom)=timass(i)
                else
                   tamass(tnatom)=amass(i)/dble(ipick(i))
                endif
c
c born radii and screening factor
c
                trborn(tnatom)=rborn(i)
                tfs(tnatom)=fs(i)
c
c set names and other things not scaled
c
                tigraph(tnatom)=igraph(i)
                tisymbl(tnatom)=isymbl(i)
                titree(tnatom)=itree(i)

c L-J pointer, leave to same atom type but add scale factor.
c then after all LES is done we'l create new types as needed using
c iac and iacfac  (in main program)
c we also need to save the original type for perturbation
c
                tiac(tnatom)=iac(i)
                tiacfac(tnatom)=iacfac(i)*ipick(i)
c box
                if (ifbox.gt.0) then
c
c now change the NSP() for this molecule to reflect that more particles
c are in the molecule
c
c find what molecule it is
c                
                     k=0
                     do 800 i1=1,nspm
c
c l is the molecule #, k is the last atom in the l molecule #
c this pointer is not available directly, have to add them up
c
                           k2=k
                           k=k+nsp(i1)
                           if (i.gt.k2.and.i.le.k) then
c
c match, add one for every atom we create
c
                                tnsp(i1)=tnsp(i1)+1
                           endif
800                   continue
c
c endif for ifbox check
c
                endif
c
c next j
c
110         continue
c
c next i
c
100     continue
c
        write (6,'(a,i6,a,a,i6,a)') 'there were ',natom,' particles; ',
     .            'currently ',tnatom,' particles'
c
c re-set residue pointers for new atom numbering
c
        do 130 i=1,nres
             j=ipres(i)
             do 140 k=1,tnatom
                  if (j.eq.poiatom(k)) then
c 
c this is 1st atom that corresponds to the original first atom of that residue
c
                         ipres(i)=k
c
c skip other matches
c
                         go to 130
                   endif
140          continue
130     continue
c
c bonds! start with zero bonds- this way all bonds will be 
c in the same order as before, just like the particles.
c we are NOT adding bonds as the end of the list, which
c would mess up the requirement in amber that constraint bonds/angles/etc
c must be at the end of the list
c
      tnbonh=0
      do 300 i=1,nbonh
c
c find how many copies of this bond will be made
c
          ncop=max(ipick(ibh(i)),ipick(jbh(i)))
c
          do 310 j=1,ncop
              tnbonh=tnbonh+1
c
c only increment particle number for LES copies,
c non-LES particles will get multiple bonds to same particle
c (or single bond if no LES at all for this bond)
c
              if (ipick(ibh(i)).eq.1) then
                  tibh(tnbonh)=revpoint(ibh(i))
              else
                  tibh(tnbonh)=revpoint(ibh(i))+j-1
              endif
              if (ipick(jbh(i)).eq.1) then
                  tjbh(tnbonh)=revpoint(jbh(i))
              else
                  tjbh(tnbonh)=revpoint(jbh(i))+j-1
              endif
c
c leave pointer for now, will fix all bond pointers at once.
c save value corresponding to the scale factor that should be used
c allow the scale factor to grow if the bond is modified over
c several LES spaces
c
              ticbh(tnbonh)=icbh(i)
              tbfach(tnbonh)=bfach(i)*ncop
310            continue
300        continue
c
      write (6,'(a,i6,a,a,i6)') 'there were ',nbonh,' nbonh bonds, ',
     .            'now there are ',tnbonh
c
c now the nbona bonds, the same way as above
c tnumc will count the number of constraints, allow us to 
c properly set nbona and nbona from tnbona
c note that there are no constraint bonds in nbonh
c
      numcb=0
c
      tnbona=0
      do 350 i=1,nbona
          ncop=max(ipick(ib(i)),ipick(jb(i)))
          do 360 j=1,ncop
              tnbona=tnbona+1
              if (ipick(ib(i)).eq.1) then
                  tib(tnbona)=revpoint(ib(i))
              else
                  tib(tnbona)=revpoint(ib(i))+j-1
              endif
              if (ipick(jb(i)).eq.1) then
                  tjb(tnbona)=revpoint(jb(i))
              else
                  tjb(tnbona)=revpoint(jb(i))+j-1
              endif
              ticb(tnbona)=icb(i)
              tbfac(tnbona)=bfac(i)*ncop
c
c increment the # constraints if needed
c        
              if (i.gt.mbona) then
                      numcb=numcb+1
              endif
360            continue
350        continue
      write (6,'(a,i6,a,a,i6)') 'there were ',nbona,' nbona bonds, ',
     .            'now there are ',tnbona
c
c now pert bonds
c
c numbb - number of boundary bonds
c
      numbb=0
      tnbper=0
      do 2350 i=1,nbper
          ncop=max(ipick(ibper(i)),ipick(jbper(i)))
          do 2360 j=1,ncop
              tnbper=tnbper+1
              if (ipick(ibper(i)).eq.1) then
                  tibper(tnbper)=revpoint(ibper(i))
              else
                  tibper(tnbper)=revpoint(ibper(i))+j-1
              endif
              if (ipick(jbper(i)).eq.1) then
                  tjbper(tnbper)=revpoint(jbper(i))
              else
                  tjbper(tnbper)=revpoint(jbper(i))+j-1
              endif
              ticbper(tnbper)=icbper(i)
c
c 2 pointers for type, for two endpoints of pert 
c (but only need one factor bfacp) but since we don't yet know how many 
c bonds, angles, etc there will be, can't store at end. have to use
c second storage area.
c
              t2icbper(tnbper)=icbper(i+nbper)
              tbfacp(tnbper)=bfacp(i)*ncop
c
c increment the # constraints if needed
c        
              if (i.gt.mbper) then
                      numbb=numbb+1
              endif
2360            continue
2350        continue
c
      write (6,'(a,i6,a,a,i6)') 'there were ',nbper,' nbper bonds, ',
     .            'now there are ',tnbper
c
c done with the bonds, move to angles
c
      tntheth=0
      do 400 i=1,ntheth
c
c find how many copies of this angle will be made
c
          ncop=max(ipick(ith(i)),ipick(jth(i)))
          ncop=max(ipick(kth(i)),ncop)
c
          do 410 j=1,ncop
              tntheth=tntheth+1
c
c only increment the particle numbers that are LES, the others
c get multiple angles created to the same particle
c
              if (ipick(ith(i)).eq.1) then
                  tith(tntheth)=revpoint(ith(i))
              else
                  tith(tntheth)=revpoint(ith(i))+j-1
              endif
              if (ipick(jth(i)).eq.1) then
                  tjth(tntheth)=revpoint(jth(i))
              else
                  tjth(tntheth)=revpoint(jth(i))+j-1
              endif
              if (ipick(kth(i)).eq.1) then
                  tkth(tntheth)=revpoint(kth(i))
              else
                  tkth(tntheth)=revpoint(kth(i))+j-1
              endif
              ticth(tntheth)=icth(i)
              tafach(tntheth)=afach(i)*ncop
410            continue
400        continue
      write (6,'(a,i6,a,a,i6)') 'there were ',ntheth,' ntheth angles, ',
     .            'now there are ',tntheth
c
c same as above for ntheta angles
c
      numca=0
      tntheta=0
      do 420 i=1,ntheta
          ncop=max(ipick(it(i)),ipick(jt(i)))
          ncop=max(ipick(kt(i)),ncop)
          do 430 j=1,ncop
              tntheta=tntheta+1
              if (ipick(it(i)).eq.1) then
                  tit(tntheta)=revpoint(it(i))
              else
                  tit(tntheta)=revpoint(it(i))+j-1
              endif
              if (ipick(jt(i)).eq.1) then
                  tjt(tntheta)=revpoint(jt(i))
              else
                  tjt(tntheta)=revpoint(jt(i))+j-1
              endif
              if (ipick(kt(i)).eq.1) then
                  tkt(tntheta)=revpoint(kt(i))
              else
                  tkt(tntheta)=revpoint(kt(i))+j-1
              endif
              tict(tntheta)=ict(i)
              tafac(tntheta)=afac(i)*ncop
c
c increment the # constraints if needed
c        
              if (i.gt.mtheta) then
                      numca=numca+1
              endif
430            continue
420        continue
      write (6,'(a,i6,a,a,i6)') 'there were ',ntheta,' ntheta angles, ',
     .            'now there are ',tntheta
c
c pert angles
c
      numba=0
      tngper=0
      do 2420 i=1,ngper
          ncop=max(ipick(itper(i)),ipick(jtper(i)))
          ncop=max(ipick(ktper(i)),ncop)
          do 2430 j=1,ncop
              tngper=tngper+1
              if (ipick(itper(i)).eq.1) then
                  titper(tngper)=revpoint(itper(i))
              else
                  titper(tngper)=revpoint(itper(i))+j-1
              endif
              if (ipick(jtper(i)).eq.1) then
                  tjtper(tngper)=revpoint(jtper(i))
              else
                  tjtper(tngper)=revpoint(jtper(i))+j-1
              endif
              if (ipick(ktper(i)).eq.1) then
                  tktper(tngper)=revpoint(ktper(i))
              else
                  tktper(tngper)=revpoint(ktper(i))+j-1
              endif
              tictper(tngper)=ictper(i)
              t2ictper(tngper)=ictper(i+ngper)
              tafacp(tngper)=afacp(i)*ncop
              if (i.gt.mgper) then
                      numba=numba+1
              endif
2430            continue
2420        continue
      write (6,'(a,i6,a,a,i6)') 'there were ',ngper,' ngper angles, ',
     .            'now there are ',tngper
c
c done with angles
c on to torsions
c
      tnphih=0
      do 500 i=1,nphih
          ncop=max(ipick(iph(i)),ipick(jph(i)))
          ncop=max(ncop,ipick(kph(i)))
          ncop=max(ncop,ipick(lph(i)))
          do 510 j=1,ncop
              tnphih=tnphih+1
              if (ipick(iph(i)).eq.1) then
                  tiph(tnphih)=revpoint(iph(i))
              else
                  tiph(tnphih)=revpoint(iph(i))+j-1
              endif
              if (ipick(jph(i)).eq.1) then
                  tjph(tnphih)=revpoint(jph(i))
              else
                  tjph(tnphih)=revpoint(jph(i))+j-1
              endif
              if (ipick(kph(i)).eq.1) then
                  tkph(tnphih)=revpoint(kph(i))
              else
                  tkph(tnphih)=revpoint(kph(i))+j-1
              endif
              if (ipick(lph(i)).eq.1) then
                  tlph(tnphih)=revpoint(lph(i))
              else
                  tlph(tnphih)=revpoint(lph(i))+j-1
              endif
c
c set the signs for impropers and 1-4
c
              tinegh(tnphih)=inegh(i)
              tjnegh(tnphih)=jnegh(i)
              tknegh(tnphih)=knegh(i)
              tlnegh(tnphih)=lnegh(i)
c
c IMPORTANT!
c check for multiple copies made of torsion when end atoms do not
c have copies- this would result in multiple 1-4 interactions
c between the non-copied and atoms. change this by making lnegh=-1
c this is because amber uses the torsion list for the 1-4 list
c unless the 3rd atom is negative
c
              if ((ipick(iph(i)).eq.1.and.ipick(lph(i)).eq.1) 
     &                        .and.(j.gt.1)) then
                      tknegh(tnphih)=-1
              endif
c
              ticph(tnphih)=icph(i)
              ttfach(tnphih)=tfach(i)*ncop
510            continue
500        continue
      write (6,'(a,i6,a,a,i6)') 'there were ',nphih,' nphih torsions, ',
     .            'now there are ',tnphih
c
c now nphia torsions
c
      numcd=0
      tnphia=0
      do 520 i=1,nphia
          ncop=max(ipick(ip(i)),ipick(jp(i)))
          ncop=max(ncop,ipick(kp(i)))
          ncop=max(ncop,ipick(lp(i)))
          do 530 j=1,ncop
              tnphia=tnphia+1
              if (ipick(ip(i)).eq.1) then
                  tip(tnphia)=revpoint(ip(i))
              else
                  tip(tnphia)=revpoint(ip(i))+j-1
              endif
              if (ipick(jp(i)).eq.1) then
                  tjp(tnphia)=revpoint(jp(i))
              else
                  tjp(tnphia)=revpoint(jp(i))+j-1
              endif
              if (ipick(kp(i)).eq.1) then
                  tkp(tnphia)=revpoint(kp(i))
              else
                  tkp(tnphia)=revpoint(kp(i))+j-1
              endif
              if (ipick(lp(i)).eq.1) then
                  tlp(tnphia)=revpoint(lp(i))
              else
                  tlp(tnphia)=revpoint(lp(i))+j-1
              endif
              tineg(tnphia)=ineg(i)
              tjneg(tnphia)=jneg(i)
              tkneg(tnphia)=kneg(i)
              tlneg(tnphia)=lneg(i)
              if ((ipick(ip(i)).eq.1.and.ipick(lp(i)).eq.1) 
     &                        .and.(j.gt.1)) then
                      tkneg(tnphia)=-1
              endif
              ticp(tnphia)=icp(i)
              ttfac(tnphia)=tfac(i)*ncop
c
c increment the # constraints if needed
c        
              if (i.gt.mphia) then
                      numcd=numcd+1
              endif
530            continue
520        continue
      write (6,'(a,i6,a,a,i6)') 'there were ',nphia,' nphia torsions, ',
     .            'now there are ',tnphia
c
c pert torsions
c
      numbd=0
      tndper=0
      do 2520 i=1,ndper
          ncop=max(ipick(ipper(i)),ipick(jpper(i)))
          ncop=max(ncop,ipick(kpper(i)))
          ncop=max(ncop,ipick(lpper(i)))
          do 2530 j=1,ncop
              tndper=tndper+1
              if (ipick(ipper(i)).eq.1) then
                  tipper(tndper)=revpoint(ipper(i))
              else
                  tipper(tndper)=revpoint(ipper(i))+j-1
              endif
              if (ipick(jpper(i)).eq.1) then
                  tjpper(tndper)=revpoint(jpper(i))
              else
                  tjpper(tndper)=revpoint(jpper(i))+j-1
              endif
              if (ipick(kpper(i)).eq.1) then
                  tkpper(tndper)=revpoint(kpper(i))
              else
                  tkpper(tndper)=revpoint(kpper(i))+j-1
              endif
              if (ipick(lpper(i)).eq.1) then
                  tlpper(tndper)=revpoint(lpper(i))
              else
                  tlpper(tndper)=revpoint(lpper(i))+j-1
              endif
              tinegp(tndper)=inegp(i)
              tjnegp(tndper)=jnegp(i)
              tknegp(tndper)=knegp(i)
              tlnegp(tndper)=lnegp(i)
              if ((ipick(ipper(i)).eq.1.and.ipick(lpper(i)).eq.1) 
     &                        .and.(j.gt.1)) then
                      tkneg(tndper)=-1
              endif
              ticpper(tndper)=icpper(i)
              t2icpper(tndper)=icpper(i+ndper)
              ttfacp(tndper)=tfacp(i)*ncop
              if (i.gt.mdper) then
                      numbd=numbd+1
              endif
2530            continue
2520        continue
      write (6,'(a,i6,a,a,i6)') 'there were ',ndper,' ndper torsions, ',
     .            'now there are ',tndper
c
c add to the exclusion list - very complicated!
c
      nexcl=0
      tnexcl=0
      totalex=0
      write (6,'(a)') 'processing exclusion list ' 
c
      do 600 i=1,natom
c
         if (debug) 
     &   write (6,*) "processing next atom: orig atom ",i
c
c first convert the particle numbers in the existing list,
c make sure all copies of a particle get the list. also exclude
c all copies of a particle that was excluded. only then do exclusion
c of 'other copies of same subspace' for LES
c
c loop over (possible) copies of this particle
c
         do 610 j=1,ipick(i)

           if (debug)
     &     write (6,*) "processing copy ",j," of ",ipick(i)
c
c j1 is the new particle number(s) corresponding to the original particle #
c
           j1=revpoint(i)+j-1

           if (debug) write (6,*) "it is NEW ATOM # ",j1

           tnumex(j1)=0
c
c loop over the original exclusions (nexcl is number of excl already
c assigned to other particles, since numex() is not a direct pointer
c in to the natex list) 
c
           if (debug) then
             write (6,*) "looping over original exclusions of ",i
             write (6,*) "there are ",nexcl+numex(i)-nexcl
           endif

           do 620 k=nexcl+1,nexcl+numex(i)
c
c oddly, if no exclusions, 1 dummy exclusion of zero is given. 
c support it... but later we have to remove these unless a particle
c STILL only has one exclusion (especially after LES exclusions).
c or is it that a particle is ALWAYS excluded from itself?
c
              if (natex(k).eq.0)  then
                 if (debug) write (6,*) "adding a dummy exclusion"
                 tnexcl=tnexcl+1
                 totalex=totalex+1
                 tnatex(tnexcl)=0
                 tnumex(j1)=tnumex(j1)+1
c
              else
c
c natex(k) is the excluded particle (orig. #, new # is revpoint(natex(k)) )
c loop over its (possible) copies
c
                if (debug) then
                  write (6,*) "atom ",natex(k)," was excluded from ",i
                endif

                if(pimd) then
                   tnexcl=tnexcl+1
                   totalex=totalex+1
                   tnatex(tnexcl)=revpoint(natex(k))+j-1
                   tnumex(j1)=tnumex(j1)+1
                else
                   if (debug)
     &                write (6,*) "adding ",ipick(natex(k)),  
     &                   " copies of ",natex(k)
                
                   do 630 k2=1,ipick(natex(k))
                      tnexcl=tnexcl+1
                      totalex=totalex+1
                      tnatex(tnexcl)=revpoint(natex(k))+k2-1

                     if (debug)
     &               write (6,*) "tnatex(",tnexcl,") is ",tnatex(tnexcl)
c
c add one to number of excluded particles for this particle
c
                      tnumex(j1)=tnumex(j1)+1
c
                      if (debug)
     &                   write (6,*) j1," now has ",tnumex(j1),
     &                      " exclusions"

630                continue
                endif
c
c endif for checking if no exclusions
c
             endif
c
c next original excluded atom 
c
620        continue
      
           if (debug)
     &      write (6,*) "done with copying original exclusions"
c
c now we have the original exclusions. add the new LES exclusions.
c loop over particles with ipick(i).ne.1 (LES in THIS SUBSPACE)
c then exclude all other particles with ipick(i).ne.1 that
c have different copy number (which is simply j in this part of the code)
c don't duplicate exclusions already in list for this particle....
c
           if (debug)
     &       write (6,*) "now excluding other copies of atom (if need)"

           if (.not.pimd.and.ipick(i).ne.1) then
c
c this is a LES particle, add the exclusions
c
c ONLY NEED TO EXCLUDE ATOMS WITH _HIGHER_ NUMBER FOR SANDER!
c
            do 640 i2=1,natom
c
c look at all other atoms: i is original atom, i2 is partner
c
             if (ipick(i2).ne.1) then
c 
c this particle also had copies made, exclude those that are not same copy #
c
               if (debug) then
                  write (6,*) "looking at orig atom",i2,
     &               " with ",ipick(i2)," copies"
                  write (6,*) "exclude copies with copy# diff from ",j
               endif
c
               do 650 i3=1,ipick(i2)
c 
                 if (debug) write (6,*) "considering copy ",i3
c
                 if (i3.ne.j) then
c
c different copy number.   see if it's already in list. if not, add it.
c
                    inlist=.false.
c
c k5 is the particle we are looking at (to possibly exclude)
c
                    k5=revpoint(i2)+i3-1

                    if (debug) write (6,*) "should exclude this atom, #"
     &               ,k5," from ",j1
c
c skip it if the atom number is <j1
c
                    if (k5.lt.j1) then
c
                      if (debug) write (6,*) 
     &                 "number is less than ",j1,", skipping"
c
                      go to 650
                    endif

                    if (debug) then
                      write (6,*) "number is larger than ",j1
                      write (6,*) 
     &                 "checking to see if it is already excluded"
                    endif
c
c loop over current exclusions for j1
c
                    do 660 i4=tnexcl-tnumex(j1)+1,tnexcl
                       if (tnatex(i4).eq.k5) then
                          inlist=.true.
                          if (debug) write (6,*) 
     &                     "it was already there, not adding it"
                       endif
660                 continue
c
                    if (.not.inlist) then
c
c add it
c
                       tnexcl=tnexcl+1
                       totalex=totalex+1
                       tnatex(tnexcl)=k5
                       tnumex(j1)=tnumex(j1)+1
c
                       if (debug) write (6,*) 
     &                  "it was not there, adding it as #",tnexcl
c
                    endif
c
c endif for i3.ne.j (different copy #)

                 endif
c
c next i3 (next copy of atom i2)
c
650            continue
c
c endif for i2 being a LES atom
c
             endif
c
c next i2
c
640        continue
c
c endif for i being a LES atom
c
          endif
c
c j1 is new number for old atom i, it is the atom for which we are curently
c adding exclusions 
c
          if (debug) then
            write (6,*) "final # of exclusions for new atom ",
     &       j1," is ",tnumex(j1)
            write (6,*) "Now sorting the exclusion list"
          endif
c
c now sort the list of exclusions for j1 in ascending order
c
c        write (6,*) "list before sorting"
c        write (6,607) (tnatex(j4),j4=tnexcl-tnumex(j1)+1,tnexcl)
607     format (12(1x,i5))
c        write (6,*) 

           bgn=tnexcl-tnumex(j1)+1
           end=tnexcl
c
           incr=1
601        incr=incr*3+1
           if (incr.le.end) go to 601
602        continue
             incr=incr/3
             do i4=incr+bgn,end
               temp=tnatex(i4)
               j4=i4
603            if (tnatex(j4-incr).gt.temp) then
                   tnatex(j4)=tnatex(j4-incr)
                   j4=j4-incr
                   if ((j4-bgn+1).le.incr) go to 604 
                   go to 603
               endif
604            tnatex(j4)=temp
             enddo
          if (incr.gt.1) goto 602
c
c        write (6,*) "list after sorting"
c        write (6,607) (tnatex(j4),j4=tnexcl-tnumex(j1)+1,tnexcl)
c        write (6,*) 
c
c done with sorting new exclusion list for atom j1
c
          if (debug) write (6,*) "Finished sorting list"
c
c get next copy of atom i (next j1)
c
610      continue
c
c get next i (next real particle)
c increment pointer into original arrays (for reading)
c
         nexcl=nexcl+numex(i)
c
600    continue
c
c now remove those annoying extra zeros from list
c
      tnexcl=0
      do 670 i=1,tnatom
c
c izero is for whether this atom has a zero for an exclusion (among others)
c inon is for non-zero exclusions
c
         izero=0
         inon=0
         k=0
         do 680 j=tnexcl+1, tnexcl+tnumex(i)
c
c see if there are any non-zero for this particle
c
                if (tnatex(j).eq.0) then
                              izero=izero+1
                else
                              inon=inon+1
                endif
680      continue
c
c keep it if only one and it is zero, otherwise remove extras
c
         if (inon.eq.0.and.izero.gt.1) then
c
c multiple zero exclusions and no real ones,
c save 1 of the zeros (delete numex(i)-1 others)
c
                  k=tnumex(i)-1
                  tnumex(i)=1
                  do 685 j=tnexcl+2,totalex-k
                              tnatex(j)=tnatex(j+k)
685               continue
                  totalex=totalex-k
c
         elseif (izero.eq.0) then
c
c do nothing, no zero exclusions are present
c
         elseif (izero.ne.0.and.inon.ne.0) then
c
c this means that there are non-zero elements as well as zeroes
c go through again, and remove ALL zeroes
c
                  i1=0
                  do 690 j=tnexcl+1, tnexcl+tnumex(i)
                        if (tnatex(j).ne.0) then
                           i1=i1+1
c move it forward
                           tnatex(tnexcl+i1)=tnatex(j)
                        endif
690               continue
                  tnumex(i)=inon
c
c move the rest of the list forward
c
                  do 695 j=tnexcl+inon+1,totalex-izero
                           tnatex(j)=tnatex(j+izero)
695               continue
                  totalex=totalex-izero
c
c endif for checking for zeroes
c
         endif
c
c next i (next atom)
c
        tnexcl=tnexcl+tnumex(i)
670     continue
c
        if (debug) write (6,*) "Finished processing exclusion list"
c
c done with exclusions!
c
c wrap it all up by reverting to normal arrays - 
c copy new prmtop info to old
c
      nbonh=tnbonh
      nbona=tnbona
      nbper=tnbper
      mbona=nbona-numcb
      mbper=nbper-numbb
      ntheth=tntheth
      ntheta=tntheta
      ngper=tngper
      mtheta=ntheta-numca
      mgper=ngper-numba
      nphih=tnphih
      nphia=tnphia
      ndper=tndper
      mphia=nphia-numcd
      mdper=ndper-numbd
      do 700 i=1,nbonh
          ibh(i)=tibh(i)
          jbh(i)=tjbh(i)
          icbh(i)=ticbh(i)
          bfach(i)=tbfach(i)
700        continue
      do 710 i=1,nbona
          ib(i)=tib(i)
          jb(i)=tjb(i)
          icb(i)=ticb(i)
          bfac(i)=tbfac(i)
710        continue
      do 715 i=1,nbper
          ibper(i)=tibper(i)
          jbper(i)=tjbper(i)
          icbper(i)=ticbper(i)
          icbper(i+nbper)=t2icbper(i)
c
c set ticbper properly for pert type now that we know nbper
c
          ticbper(i+nbper)=t2icbper(i)

          bfacp(i)=tbfacp(i)
715        continue
      do 720 i=1,ntheth
          ith(i)=tith(i)
          jth(i)=tjth(i)
          kth(i)=tkth(i)
          icth(i)=ticth(i)
          afach(i)=tafach(i)
720        continue
      do 730 i=1,ntheta
          it(i)=tit(i)
          jt(i)=tjt(i)
          kt(i)=tkt(i)
          ict(i)=tict(i)
          afac(i)=tafac(i)
730        continue
      do 735 i=1,ngper 
          itper(i)=titper(i)
          jtper(i)=tjtper(i)
          ktper(i)=tktper(i)
          ictper(i)=tictper(i)
          ictper(i+ngper)=t2ictper(i)
          tictper(i+ngper)=t2ictper(i)
          afacp(i)=tafacp(i)
735        continue
      do 740 i=1,nphih
          iph(i)=tiph(i)
          jph(i)=tjph(i)
          kph(i)=tkph(i)
          lph(i)=tlph(i)
          icph(i)=ticph(i)
          tfach(i)=ttfach(i)
          inegh(i)=tinegh(i)
          jnegh(i)=tjnegh(i)
          knegh(i)=tknegh(i)
          lnegh(i)=tlnegh(i)
740        continue
      do 750 i=1,nphia
          ip(i)=tip(i)
          jp(i)=tjp(i)
          kp(i)=tkp(i)
          lp(i)=tlp(i)
          icp(i)=ticp(i)
          tfac(i)=ttfac(i)
          ineg(i)=tineg(i)
          jneg(i)=tjneg(i)
          kneg(i)=tkneg(i)
          lneg(i)=tlneg(i)
750        continue
      do 755 i=1,ndper
          ipper(i)=tipper(i)
          jpper(i)=tjpper(i)
          kpper(i)=tkpper(i)
          lpper(i)=tlpper(i)
          icpper(i)=ticpper(i)
          icpper(i+ndper)=t2icpper(i)
          ticpper(i+ndper)=t2icpper(i)
          tfacp(i)=ttfacp(i)
          inegp(i)=tinegp(i)
          jnegp(i)=tjnegp(i)
          knegp(i)=tknegp(i)
          lnegp(i)=tlnegp(i)
755        continue
c
      natom=tnatom
      do 760 i=1,natom
          x(i)=tx(i)
          y(i)=ty(i)
          z(i)=tz(i)
          vx(i)=tvx(i)
          vy(i)=tvy(i)
          vz(i)=tvz(i)
          nlev(i)=tnlev(i)
          origpt(i)=torigpt(i)
             chrg(i)=tchrg(i)
          lreal(i)=treal(i)
             if (ifpert.gt.0) then
              cgper(i)=tcgper(i)
              igrper(i)=tigrper(i)
              ismper(i)=tismper(i)
              almper(i)=talmper(i)
              iaper(i)=tiaper(i)
              iacper(i)=tiacper(i)
                  iacpfac(i)=tiacpfac(i)
          endif
          amass(i)=tamass(i)
          if (itimass) timass(i)=ttimass(i)
          igraph(i)=tigraph(i)
          rborn(i)=trborn(i)
          fs(i)=tfs(i)
          isymbl(i)=tisymbl(i)
          iac(i)=tiac(i)
          iacfac(i)=tiacfac(i)
          itree(i)=titree(i)
          numex(i)=tnumex(i)
          do 770 j=1,nlev(i)
              lesid(i,j)=tlesid(i,j)
770            continue
760        continue
c
      next=tnexcl
      do 780 i=1,next
              natex(i)=tnatex(i)
780        continue
c
c molecule numbers for box
c
      if (ifbox.gt.0) then
              do 830 i=1,nspm
                      nsp(i)=tnsp(i)
830                continue
      endif

      write (6,'(a)') 'finished creating LES subspace '
c
c check array sizes, though this may not always catch overflow
c
      call checksz( size_error )
      if (size_error) then
         stop
      endif

      return
      end
