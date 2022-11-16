c
c OPTIONS:
c
c        file: open a file, also use one of :
c                rcrd :         read coords from this file
c                rcvd:        read coords + velo from file
c                rcvb:        read coords, velo and box from file
c                wcrd:        write coords (and more if rcvd, rcvb) to file
c                wprm:        write new topology file
c                les0:         write gibbs topology for single to multiple 
c                        corresponding to original lambda=0
c                les1:         write gibbs topology for single to multiple 
c                        corresponding to original lambda=1
c        
c        action:        start run, all of the following options must come AFTER action
c
c        spac:         add a new subspace definition, using a pick command
c
c        defp:         define which atoms should be in the 'perturbed'
c                 group during the single -> multiple perts.
c                uses a pick command to select the atoms.
c                IN ADDITION TO LES ATOMS WHICH MUST BE INCLUDED!
c        
c        omas:        leave all masses at original values (otherwise scale 1/N)
c        pimd:        write topology file for pimd(do not exclude atoms from other copy)
c        scbx:        scale box size during run (no longer used)
c        bigm:        allows user to pick a set of atoms and assign a new mass
c        pert:        must be specified to use perturbation (gibbs) topologies
c
c here are a few that control how the scaling of the LES extra copies
c are handled for gibbs topologies: note that one must be very careful 
c when changing the relative weight of different terms in potential
c functions, and also when creating a LES system where the endpoints do not
c exactly match the real system, introducing error into to thermodynamic
c cycle being modeled.
c
c        sdih:        dihedral terms will be scaled to 1/N during the multiple
c                copy to multiple copy stage.  during the single to multiple
c                stages, leave 'extra' copies at 1/N while 'real' copy goes to 1
c
c        sdi2:        dihedral terms will be scaled to 1/N during the multiple
c                copy to multiple copy stage.  during the single to multiple
c                stages, all copies will go to scale factor of 1 (not just 
c                the 1 'real' copy)
c
c        sdi3:   dihedral terms will be scaled to 1/N during the multiple
c             copy to multiple copy stage.  during the single to multiple
c             stages, scale extra copies to 0 while 'real' copy goes to 1
c             THIS MATCHES THE ENDPOINT OF THE REAL SYSTEM!
c
c        ndih:        dihedral terms will be left at 1 (not scaled) during the 
c                multiple copy to multiple copy stage, and also left at 1 during
c                single to multiple stages.
c
c        ndi2:        dihedral terms will be left at 1/N (scaled) during the 
c                multiple copy to multiple copy stage, and also left at 1/N 
c                during single to multiple stages.
c
c         EITHER SDIH, SDI2 NDIH or NDI2 MUST BE USED (for perturbations)!
c
c        sang:        angle terms will be scaled to 1/N during the multiple
c                copy to multiple copy stage.  during the single to multiple
c                stages, leave 'extra' copies at 1/N while 'real' copy goes to 1
c
c        san2:        angle terms will be scaled to 1/N during the multiple
c                copy to multiple copy stage.  during the single to multiple
c                stages, all copies will go to scale factor of 1 (not just 
c                the 1 'real' copy)
c
c        nang:        angle terms will be left at 1 (not scaled) during the 
c                multiple copy to multiple copy stage, and also left at 1 during
c                single to multiple stages.
c
c        nan2:        angle terms will be left at 1/N (scaled) during the 
c                multiple copy to multiple copy stage, and also left at 1/N during
c                single to multiple stages.
c
c         EITHER SANG, SAN2 NANG or NAN2 MUST BE USED (for perturbations)!
c
c        sbon:        bond terms will be scaled to 1/N during the multiple
c                copy to multiple copy stage.  during the single to multiple
c                stages, leave 'extra' copies at 1/N while 'real' copy goes to 1
c
c        sbo2:        bond terms will be scaled to 1/N during the multiple
c                copy to multiple copy stage.  during the single to multiple
c                stages, all copies will go to scale factor of 1 (not just 
c                the 1 'real' copy)
c
c        nbon:        bond terms will be left at 1 (not scaled) during the 
c                multiple copy to multiple copy stage, and also left at 1 during
c                single to multiple stages.
c
c        nbo2:        bond terms will be left at 1/N (scaled) during the 
c                multiple copy to multiple copy stage, and also left at 1/N during
c                single to multiple stages.
c
c         EITHER SBON, SBO2 NBON or NBO2 MUST BE USED (for perturbations)!
c
c        nonbond force constants and charges are always scaled 1/N in the LES
c        state, and scaled to either 1 for 'real' copy or 0 for 'extra' copies
c        during the single to multiple stages. other options would not be
c        appropriate.        
c
c main program for addles
c
#include "SIZE.h"
#include "TOP.h"
#include "LINE.h"
#include "TTOP.h"
#include "MISC.h"
#include "DEBUG.h"
      integer hidx,imatch,l,i1
      character*1 chgbox
      logical wprm2ok,wprm3ok,defspc
      
      
      call amrset(178753)

      allmas=.false.
      defspc=.false.
      bigmas=.false.
      lespert=.false.
      nomodv=.false.

      rcrdok=.false.
      rcvdok=.false.
      rcvbok=.false.
      wcrdok=.false.
      wnmrok=.false.
      rprmok=.false.
      wprmok=.false.
      wprm2ok=.false.
      wprm3ok=.false.
      pimd = .false.
      atm1st=.false.
      npack=1
c
c get all the input info, filenames, etc.
c
      write(6,'(t2,''Amber12 Module: addles'')')
      write(6,'(t2,''set up Locally Enhanced Sampling topology'')')

      stdi = 5
        stderr = 0
      name='add_les'
      namel=7
      dumiac=0
      dumdih=0
c
c potential scaling options for perturbations
c defaults: if value=1, scale bonds,angles,dihedrals by 1/N
c and while removing copies, IF the copies were scaled by 1/N,
c leave 'extra' copies at 1/N while 'real' copy force constants go to 1
c
c for now FORCE user to choose behavior for each flag. valid values
c will be 1, 2 or 3.
c
      dscale=0
      ascale=0
      bscale=0
c
c
c open junk file for rline
c
        jnkf = 25
        open (unit=jnkf,status='scratch')
c
c get the input 
c
1       continue
        call rline(name,namel,stdi)
        if (find('debu')) then
           debug = .true.
        else if (find('file')) then
           if (find('rprm')) then
c
c get topology 
c
              urprm = of()
              call readprm(urprm)
c
c check errors
c
              rprmok=.true.
          endif
          if (find('rcrd')) then
c
c read coords
c
              urcrd=of()
              print *, 'urcrd:', urcrd
c
c check errors
c
              rcrdok=.true.
              write (6,'(a,i4)') 'Coordinates only from unit',
     &                urcrd
              call readcrd(urcrd)
          endif
c
c read coords and velocities (restart file format)
c
          if (find('rcvd')) then
              urcrd=of()
              rcvdok=.true.
              write (6,'(a,i4)') 'Coordinates and velocities from unit',
     &                urcrd
              call readcrd(urcrd)
          endif
c
c read coords, velocities and box size (restart file format)
c
          if (find('rcvb')) then
              urcrd=of()
              rcvbok=.true.
              write (6,'(a,i4)') 'Coords, velocities and box from unit',
     &                urcrd
              call readcrd(urcrd)
          endif

          if (find('rcbd')) then
              urcrd=of()
              rcbdok=.true.
              write(6,'(a,i4)') 'Coords, box from unit', urcrd
              call readcrd(urcrd)
          endif            

          if (find('wprm')) then
              uwprm=of()
              wprmok=.true.
          endif        
          if (find('les0')) then
              uw2prm=of()
              wprm2ok=.true.
          endif        
          if (find('les1')) then
              uw3prm=of()
              wprm3ok=.true.
          endif        
          if (find('wcrd')) then
              uwcrd=of()
              wcrdok=.true.
          endif        
          if (find('wnmr')) then
              uwnmr=of()
              wnmrok=.true.
          endif        
      endif
      if (find('acti')) go to 2
c
c loop through next line
c
        go to 1
c

2       continue
c
c check for proper file existence
c
      if (.not.rprmok) then
              write (6,*) 'File for old topology not specified'
              stop
      endif
      if (.not.wprmok) then
              write (6,*) 'File for new topology not specified'
              stop
      endif
      if (.not.rcrdok.and..not.rcvdok.and..not.rcvbok) then
              write (6,*) 'Warning- old coordinate file not specified'
              write (6,*) 'New coordinate file  will not be created'
      endif
      if (.not.wcrdok) then
              write (6,*) 'Warning- new coordinate file not specified'
              write (6,*) 'New coordinate file  will not be created'
      endif
c
c set default LES params for the particles
c then set the initial scale factors for bond, angle and torsion types
c
      do 100 i=1,natom
          nlev(i)=0
c
c set the first element in the lesid array for each atom to be 
c zero. this way we can write lesid(i,1) for all atoms i
c in the topology file, it will be used for NMR restraints.
c if lesid(i,1)=0 for atom i, then LES was not used for that atom.
c if it is not zero, it corresponds to the subspace number for the
c FIRST subspace the atom was copied in. Currently the NMR stuff will
c assume 1 subspace per atom, although more are possible.
c
          lesid(i,1)=0
c
c set atom type scale factor, used to create new atom type list
c
          iacfac(i)=1
          iacpfac(i)=1
c
c set pointer to orig particle number, this will be saved and used
c for subsequent picks that may specify particle number for the
c original topology
c
          origpt(i)=i
c
c set orig charges,mass
c
          ochrg(i)=chrg(i)
          ocgper(i)=cgper(i)
          omass(i)=amass(i)
c            if (ifpert.gt.0) ocgper(i)=cgper(i)
          totcop(i)=1
          lreal(i)=.true.
c
100        continue
c
      do 105 i=1,nbonh
          bfach(i)=1
105   continue
      do 110 i=1,nbona
          bfac(i)=1
110   continue
      do 112 i=1,nbper
          bfacp(i)=1
112   continue
      do 115 i=1,ntheth
          afach(i)=1
115   continue
      do 120 i=1,ntheta
          afac(i)=1
120   continue
      do 122 i=1,ngper
          afacp(i)=1
122   continue
      do 125 i=1,nphih
          tfach(i)=1
125   continue
      do 130 i=1,nphia
          tfac(i)=1
130   continue
      do 132 i=1,ndper
          tfacp(i)=1
132   continue
c 
c write molecule info
c
c get ready to make copies
c
      curlesid=0
c
c see if a les space definition is present (or other commands)
c
20        call rline(name,namel,stdi)
c
        if (find('pert')) then
           if (defspc) then
               write (6,*) 'no commands allowed after spac commands!'
               stop
           endif
           lespert=.true.
           write (6,*) '********************************************'
           write (6,*) 'ADDLES does not yet support perturbations'
           write (6,*) 'Check the AMBER web page for updates'
           write (6,*) '********************************************'
           stop
c
c           write (6,*) '********************************************'
c           write (6,*) '3 perturbation topology files will be created'
c           write (6,*) '********************************************'
c           write (6,*) 'Use of the omas option STRONGLY recommended!'
c           write (6,*) '********************************************'
c           if (.not.wprm2ok) then
c                   write (6,*) 'wpr1 filename not specified'
c                   stop
c           endif
c           if (.not.wprm3ok) then
c                   write (6,*) 'wpr2 filename not specified'
c                   stop
c           endif
c           go to 20
c
        else if (find('pimd')) then
           pimd = .true.
           write(*,*) 'writing pimd prmtop'
           go to 20
c
        else if (find('a1st')) then
           atm1st = .true.
           write(*,*) 'writing copy atoms in atom 1st order'
           go to 20
c  
        else if (find('omas')) then
           if (defspc) then
               write (6,*) 'no commands allowed after spac commands!'
               stop
           endif
           allmas=.true.
           write (6,*) 'All masses will be left at initial values' 
           go to 20
c
        else if (find('sdih')) then
          if (defspc) then
               write (6,*) 'no commands allowed after spac commands!'
               stop
           endif
           write (6,*) 'Dihedral force constants will be scaled 1/N'
           write (6,*) 'Dummy copies will be left at 1/N when removed'
           write (6,*) 'Real copy will go to 1 when dummies removed'
           dscale=1
           go to 20
c
        else if (find('sdi2')) then
           if (defspc) then
               write (6,*) 'no commands allowed after spac commands!'
               stop
           endif
           write (6,*) 'Dihedral force constants will be scaled 1/N'
           write (6,*) 
     &          'Dummy copies will scale to full value when removed'
           write (6,*) 'Real copy will go to 1 when dummies removed'
           dscale=2
           go to 20
c
        else if (find('sdi3')) then
              if (defspc) then
                  write (6,*) 'no commands allowed after spac commands!'
                  stop
              endif
              write (6,*) 'Dihedral force constants will be scaled 1/N'
              write (6,*)
     &                   'Dummy copies will scale to zero when removed'
              write (6,*) 'Real copy will go to 1 when dummies removed'
              dscale=5
              go to 20
c
        else if (find('ndih')) then
          if (defspc) then
               write (6,*) 'no commands allowed after spac commands!'
               stop
           endif
           write (6,*) 'Dihedral force constants will NEVER be scaled'
           dscale=3
           go to 20
c
        else if (find('ndi2')) then
          if (defspc) then
               write (6,*) 'no commands allowed after spac commands!'
               stop
           endif
           write (6,*) 'Dihedral force constants will ALWAYS be scaled'
           dscale=4
           go to 20
c
        else if (find('sang')) then
          if (defspc) then
               write (6,*) 'no commands allowed after spac commands!'
               stop
           endif
           write (6,*) 'Angle force constants will be scaled 1/N'
           write (6,*) 'Dummy copies will be left at 1/N when removed'
           write (6,*) 'Real copy will go to 1 when dummies removed'
           ascale=1
           go to 20
c
        else if (find('san2')) then
          if (defspc) then
               write (6,*) 'no commands allowed after spac commands!'
               stop
           endif
           write (6,*) 'Angle force constants will be scaled 1/N'
           write (6,*) 
     &          'Dummy copies will scale to full value when removed'
           write (6,*) 'Real copy will go to 1 when dummies removed'
           ascale=2
           go to 20
c
        else if (find('nang')) then
          if (defspc) then
               write (6,*) 'no commands allowed after spac commands!'
               stop
           endif
           write (6,*) 'Angle force constants will NEVER be scaled'
           ascale=3
           go to 20
c
        else if (find('nan2')) then
          if (defspc) then
               write (6,*) 'no commands allowed after spac commands!'
               stop
           endif
           write (6,*) 'Angle force constants will ALWAYS be scaled'
           ascale=4
           go to 20
c
        else if (find('sbon')) then
          if (defspc) then
               write (6,*) 'no commands allowed after spac commands!'
               stop
           endif
           write (6,*) 'Bond force constants will be scaled 1/N'
           write (6,*) 'Dummy copies will be left at 1/N when removed'
           write (6,*) 'Real copy will go to 1 when dummies removed'
           bscale=1
           go to 20
c
        else if (find('sbo2')) then
           if (defspc) then
               write (6,*) 'no commands allowed after spac commands!'
               stop
           endif
           write (6,*) 'Bond force constants will be scaled 1/N'
           write (6,*) 
     &          'Dummy copies will scale to full value when removed'
           write (6,*) 'Real copy will go to 1 when dummies removed'
           bscale=2
           go to 20
c
        else if (find('nbon')) then
           if (defspc) then
               write (6,*) 'no commands allowed after spac commands!'
               stop
           endif
           write (6,*) 'Bond force constants will NEVER be scaled'
           bscale=3
           go to 20
c
        else if (find('nbo2')) then
           if (defspc) then
               write (6,*) 'no commands allowed after spac commands!'
               stop
           endif
           write (6,*) 'Bond force constants will ALWAYS be scaled'
           bscale=4
           go to 20
c
        else if (find('defp')) then
           if (defspc) then
               write (6,*) 'no commands allowed after spac commands!'
               stop
           endif
           call pick(ipert,1,origpt)
           write (6,*) 'User-defined perturbed group'
           j=0
           do 1650 i=1,natom
               if (ipert(i).gt.0) j=j+1
1650       continue
           write (*,1651) j,natom
1651       format ('Picked ',i5,' out of ',i5,' atoms for pert group')
           go to 20
c
        else if (find('nomodv')) then
          write (6,*) 'copy velocities will NOT be slightly modfied'
          nomodv=.true.
c
        else if (find('bigm')) then
           if (defspc) then
               write (6,*) 'no commands allowed after spac commands!'
               stop
           endif
           bigmas=.true.
           newm=getd('newm',10.d0)
           call pick(imass,1,origpt)
           write (6,*) 'Masses of picked atoms will be ',newm
           go to 20
c
      else if (find('spac')) then
c
c don't allow any other commands anymore- atom numbering, etc changes.
c more space definitions are ok, but nothing else.
c we want the pick arrays from options like bigm to correspond to the
c original particle numbering scheme. then we can reference them by
c imass(origpt(i)), for example.
c
          defspc=.true.
          if (lespert.and.(dscale.eq.0.or.ascale.eq.0
     .         .or.bscale.eq.0)) then
              write (6,*) 
     &        'YOU MUST SET bond, angle and dihedral scaling options'
              stop
          endif
c
          numcop=geti('numc',1)
          if (numcop.lt.1) then
              write (6,*) 'Illegal number of copies: ',numcop
              stop
          elseif (numcop.gt.maxcopy) then
              write (6,*) 'exceeded maxcopy',maxcopy
              stop
          elseif (numcop .eq. 1) then
              write (6,*) 'partial REM topology will be created'
          endif
c
c get a list of particles to be copied in this subspace
c
          call pick(ipick,numcop,origpt)
c
c define partial REM atoms
c
         if(numcop .eq. 1) then
            do i = 1, natom
                nrem(i) = ipick(i)
            enddo
         endif
c
c write some info on the pick
c
          write (6,*) 'picking from ',natom,' particles'
          npick=0
          do 22 i=1,natom
              if (ipick(i).eq.numcop) then
                   npick=npick+1
              elseif (ipick(i).eq.0) then
                   ipick(i)=1
              else
                   write (6,*) 'picking error ',i,ipick(i)
              endif
22        continue
67        format ('atom, #, orig atom #, copies to make: ',a,
     .        3(1x,i5),1x)
c
          write (6,*) 'Picked ',npick,' particles'
          write (6,*) 'Making ',numcop,' copies'
c
          curlesid=curlesid+1
          ncopies(curlesid)=numcop
          if( atm1st ) then
              call addspace_atm1st(ipick,curlesid)
          else
              call addspace(ipick,curlesid)
          end if
          go to 20
c
      elseif (.not.find('*EOD')) then
          write (6,*) 'Problem in line read...'
          stop
      endif
c
      if (curlesid.eq.0) then
              write (6,*) 'No subspaces defined!'
              stop
      endif
      write (6,*) 'Finished reading subspace definitions. '
      write (6,*) 'Looking for unique atom and covalent types'
c
c must be done with LES. clean up, set params, write topology and coords.
c
c set the bond type list- make extra types as necessary and
c assign the types to the bond pointers
c didn't do this earlier since there might be duplication of
c types- best to wait until after all new bonds were made
c in the addspace part of the program, the bond pointers were left
c as is but a scale factor was added as a separate vector (bfach and bfac)
c
      tnumbnd=0
      do 200 i=1,nbonh
c
c look at the original bond pointer, and the scaling factor bfach()
c loop through the new list of bonds and see if it is already here
c if not add the force constant and equ length trk() and treq(), getting
c them from the 'real' rk() and req(). of course treq() is NOT scaled
c only trk is scaled. use tbon1() to tell which of the orig. bond types
c it came from, and tbon2() to tell what the scale factor was.
c that's how we can tell if a new type has been added yet or not.
c
c an alternate method would be to just look for unique req and rk, but
c that is more likely to be influenced by roundoff, etc. which might
c make a larger list with bond types that are nearly identical
c
          imatch=0
          do 210 j=1,tnumbnd
c
              if (icbh(i).eq.tbon1(j)) then
c
c these came from same original bond type
c
                  if (bfach(i).eq.tbon2(j).or.
     &                        (lespert.and.bscale.eq.3)) then
c
c they had the same scale factor too (ignore if we aren't scaling bonds)
c
                      imatch=j
c
                  endif
              endif
210       continue
          if (imatch.eq.0) then
c
c add this as a new bond type to tnumbnd
c
              tnumbnd=tnumbnd+1
              treq(tnumbnd)=req(icbh(i))
c
c if doing LES pert, scale it based on user choice of bscale
c bscale=1: scale it
c bscale=2: scale it
c bscale=3: don't scale
c
              if (lespert.and.bscale.eq.3) then
                      trk(tnumbnd)=rk(icbh(i))
c
c tbon2 should change based on bscale behavior - whether this bond WAS scaled
c
                      tbon2(tnumbnd)=1
              else
                      trk(tnumbnd)=rk(icbh(i))/dble(bfach(i))
                      tbon2(tnumbnd)=bfach(i)
              endif

              print *, icbh(i), tnumbnd, treq(tnumbnd),trk(tnumbnd)
c
c set the pointers to tell us where this bond came from
c (in the orig bond type list)
c
              tbon1(tnumbnd)=icbh(i)
c
c set the bond pointer to this new type
c
              icbh(i)=tnumbnd
          else
              icbh(i)=imatch
          endif
c
c now, if (ifpert) and want to write the les-noles pert topologies
c and unscaled bonds are ever needed (bscale = 1 or 2)
c we should LEAVE THE ORIGINAL TYPE IN THE LIST! will make list
c bigger but then we only need one type list for the whole set
c of perts needed for LES calculation of free energies.
c
c          if (.not.lespert.or.bscale.eq.3) go to 200
c
          imatch=0
          do 1210 j=1,tnumbnd
c
c using ticbh since it was left as icbh last addspace and icbh() has been 
c changed in the lines just above...
c so it still contains the original pointer. just set scale factor to 1
c
              if (ticbh(i).eq.tbon1(j)) then
                  if (1.eq.tbon2(j)) then
                      imatch=j
                  endif
              endif
1210      continue
          if (imatch.eq.0) then
              tnumbnd=tnumbnd+1
              treq(tnumbnd)=req(ticbh(i))
              trk(tnumbnd)=rk(ticbh(i))
              tbon1(tnumbnd)=ticbh(i)
              tbon2(tnumbnd)=1
              ticbh(i)=tnumbnd
              print *, tnumbnd, treq(tnumbnd), trk(tnumbnd)
          else
              ticbh(i)=imatch
          endif
200   continue

      print *, 'bond with hydrogen:', tnumbnd
c
c now do the nbona bonds, same way
c
      do 220 i=1,nbona
          imatch=0
          do 230 j=1,tnumbnd
              if (icb(i).eq.tbon1(j)) then
                  if (bfac(i).eq.tbon2(j).or.
     &                        (lespert.and.bscale.eq.3)) then
                      imatch=j
                  endif
              endif
230       continue
          if (imatch.eq.0) then
              tnumbnd=tnumbnd+1
              treq(tnumbnd)=req(icb(i))
              if (lespert.and.bscale.eq.3) then
                      trk(tnumbnd)=rk(icb(i))
                      tbon2(tnumbnd)=1
              else
                      trk(tnumbnd)=rk(icb(i))/dble(bfac(i))
                      tbon2(tnumbnd)=bfac(i)
              endif
              tbon1(tnumbnd)=icb(i)
              icb(i)=tnumbnd
          else
                  icb(i)=imatch
          endif
c
c          if (.not.lespert.or.bscale.eq.3) go to 220
          imatch=0
          do 1230 j=1,tnumbnd
              if (ticb(i).eq.tbon1(j)) then
                  if (1.eq.tbon2(j)) then
                      imatch=j
                  endif
              endif
1230      continue
          if (imatch.eq.0) then
              tnumbnd=tnumbnd+1
              treq(tnumbnd)=req(ticb(i))
              trk(tnumbnd)=rk(ticb(i))
              tbon1(tnumbnd)=ticb(i)
              tbon2(tnumbnd)=1
              ticb(i)=tnumbnd
          else
              ticb(i)=imatch
          endif
c
220   continue
c
c now do the nbper bonds, same way except two types per bond!
c
      do 2220 i=1,nbper
          imatch=0
          do 2230 j=1,tnumbnd
              if (icbper(i).eq.tbon1(j)) then
                  if (bfacp(i).eq.tbon2(j).or.
     &                        (lespert.and.bscale.eq.3)) then
                      imatch=j
                  endif
              endif
2230      continue
          if (imatch.eq.0) then
              tnumbnd=tnumbnd+1
              treq(tnumbnd)=req(icbper(i))
              if (lespert.and.bscale.eq.3) then
                      trk(tnumbnd)=rk(icbper(i))
                      tbon2(tnumbnd)=1
              else
                      trk(tnumbnd)=rk(icbper(i))/dble(bfacp(i))
                      tbon2(tnumbnd)=bfacp(i)
              endif
              tbon1(tnumbnd)=icbper(i)
              icbper(i)=tnumbnd
          else
                  icbper(i)=imatch
          endif
c now lambda=1
          imatch=0
          do 2235 j=1,tnumbnd
              if (icbper(i+nbper).eq.tbon1(j)) then
                  if (bfacp(i).eq.tbon2(j).or.
     &                        (lespert.and.bscale.eq.3)) then
                      imatch=j
                  endif
              endif
2235      continue
          if (imatch.eq.0) then
              tnumbnd=tnumbnd+1
              treq(tnumbnd)=req(icbper(i+nbper))
              if (lespert.and.bscale.eq.3) then
                      trk(tnumbnd)=rk(icbper(i+nbper))
                      tbon2(tnumbnd)=1
              else
                      trk(tnumbnd)=rk(icbper(i+nbper))/dble(bfacp(i))
                      tbon2(tnumbnd)=bfacp(i)
              endif
              tbon1(tnumbnd)=icbper(i+nbper)
              icbper(i+nbper)=tnumbnd
          else
                  icbper(i+nbper)=imatch
          endif
c
c          if (.not.lespert.or.bscale.eq.3) go to 2220
          imatch=0
          do 1240 j=1,tnumbnd
              if (ticbper(i).eq.tbon1(j)) then
                  if (1.eq.tbon2(j)) then
                      imatch=j
                  endif
              endif
1240      continue
          if (imatch.eq.0) then
              tnumbnd=tnumbnd+1
              treq(tnumbnd)=req(ticbper(i))
              trk(tnumbnd)=rk(ticbper(i))
              tbon1(tnumbnd)=ticbper(i)
              tbon2(tnumbnd)=1
              ticbper(i)=tnumbnd
          else
              ticbper(i)=imatch
          endif
          imatch=0
          do 1250 j=1,tnumbnd
              if (ticbper(i+nbper).eq.tbon1(j)) then
                  if (1.eq.tbon2(j)) then
                      imatch=j
                  endif
              endif
1250      continue
          if (imatch.eq.0) then
              tnumbnd=tnumbnd+1
              treq(tnumbnd)=req(ticbper(i+nbper))
              trk(tnumbnd)=rk(ticbper(i+nbper))
              tbon1(tnumbnd)=ticbper(i+nbper)
              tbon2(tnumbnd)=1
              ticbper(i+nbper)=tnumbnd
          else
              ticbper(i+nbper)=imatch
          endif
c
2220  continue
c
      write (6,'(a,i5,a,a,i5)') 'There were ',numbnd,' bond types, ',
     .            'now there are ',tnumbnd
c
c move these to the real arrays.
c bfac(), bfach(), tbon1() and tbon2() not needed anymore
c
      numbnd=tnumbnd
      do 240 i=1,numbnd
          rk(i)=trk(i)
          req(i)=treq(i)
240   continue
c
c done with bonds, do angles essentially the same, so no explanations
c still check to see if we need to keep original types in list
c
      tnumang=0
      do 300 i=1,ntheth
          imatch=0
          do 310 j=1,tnumang
              if (icth(i).eq.tang1(j)) then
                  if (afach(i).eq.tang2(j).or.
     &                        (lespert.and.ascale.eq.3)) then
                      imatch=j
                  endif
              endif
310       continue
          if (imatch.eq.0) then
              tnumang=tnumang+1
              tteq(tnumang)=teq(icth(i))
              if (lespert.and.ascale.eq.3) then
c
c icth should point to the scaled angle, and ticth will be the full angle
c
                      ttk(tnumang)=tk(icth(i))
                      tang2(tnumang)=1
              else
                      ttk(tnumang)=tk(icth(i))/dble(afach(i))
                      tang2(tnumang)=afach(i)
              endif
              tang1(tnumang)=icth(i)
              icth(i)=tnumang
          else
              icth(i)=imatch
          endif
c          if (.not.lespert.or.ascale.eq.3) go to 300
c
c now add in the full force constant angle type for single -> mult perts
c use ticth since icth has changed
c
          imatch=0
          do 1255 j=1,tnumang
              if (ticth(i).eq.tang1(j)) then
                  if (1.eq.tang2(j)) then
                      imatch=j
                  endif
              endif
1255      continue
          if (imatch.eq.0) then
              tnumang=tnumang+1
              tteq(tnumang)=teq(ticth(i))
              ttk(tnumang)=tk(ticth(i))
              tang1(tnumang)=ticth(i)
              tang2(tnumang)=1
              ticth(i)=tnumang
c
c now ticth(i) is the original, full force constant angle for i
c
          else
              ticth(i)=imatch
          endif
300   continue
      do 320 i=1,ntheta
          imatch=0
          do 330 j=1,tnumang
              if (ict(i).eq.tang1(j)) then
                  if (afac(i).eq.tang2(j).or.
     &                        (lespert.and.ascale.eq.3)) then
                      imatch=j
                  endif
              endif
330       continue
          if (imatch.eq.0) then
              tnumang=tnumang+1
              tteq(tnumang)=teq(ict(i))
              if (lespert.and.ascale.eq.3) then
                      ttk(tnumang)=tk(ict(i))
                      tang2(tnumang)=1
              else
                      ttk(tnumang)=tk(ict(i))/dble(afac(i))
                      tang2(tnumang)=afac(i)
              endif
              tang1(tnumang)=ict(i)
              ict(i)=tnumang
          else
              ict(i)=imatch
          endif
c          if (.not.lespert.or.ascale.eq.3) go to 320
          imatch=0
          do 1260 j=1,tnumang
              if (tict(i).eq.tang1(j)) then
                  if (1.eq.tang2(j)) then
                      imatch=j
                  endif
              endif
1260      continue
          if (imatch.eq.0) then
              tnumang=tnumang+1
              tteq(tnumang)=teq(tict(i))
              ttk(tnumang)=tk(tict(i))
              tang1(tnumang)=tict(i)
              tang2(tnumang)=1
              tict(i)=tnumang
          else
              tict(i)=imatch
          endif
320   continue
      do 2320 i=1,ngper
          imatch=0
          do 2330 j=1,tnumang
              if (ictper(i).eq.tang1(j)) then
                  if (afacp(i).eq.tang2(j).or.
     &                        (lespert.and.ascale.eq.3)) then
                      imatch=j
                  endif
              endif
2330      continue
          if (imatch.eq.0) then
              tnumang=tnumang+1
              tteq(tnumang)=teq(ictper(i))
              if (lespert.and.ascale.eq.3) then
                      ttk(tnumang)=tk(ictper(i))
                      tang2(tnumang)=1
              else
                      ttk(tnumang)=tk(ictper(i))/dble(afacp(i))
                      tang2(tnumang)=afacp(i)
              endif
              tang1(tnumang)=ictper(i)
              ictper(i)=tnumang
          else
              ictper(i)=imatch
          endif
c
c now other lambda
c
          imatch=0
          do 2332 j=1,tnumang
              if (ictper(i+ngper).eq.tang1(j)) then
                  if (afacp(i).eq.tang2(j).or.
     &                        (lespert.and.ascale.eq.3)) then
                      imatch=j
                  endif
              endif
2332      continue
          if (imatch.eq.0) then
              tnumang=tnumang+1
              tteq(tnumang)=teq(ictper(i+ngper))
              if (lespert.and.ascale.eq.3) then
                      ttk(tnumang)=tk(ictper(i+ngper))
                      tang2(tnumang)=1
              else
                      ttk(tnumang)=tk(ictper(i+ngper))/dble(afacp(i))
                      tang2(tnumang)=afacp(i)
              endif
              tang1(tnumang)=ictper(i+ngper)
              ictper(i+ngper)=tnumang
          else
              ictper(i+ngper)=imatch
          endif
c          if (.not.lespert.or.ascale.eq.3) go to 2320
c
c now get the unscaled versions for single -> mult perts
c use tictper since ictper was changed above
c
          imatch=0
          do 1265 j=1,tnumang
              if (tictper(i).eq.tang1(j)) then
                  if (1.eq.tang2(j)) then
                      imatch=j
                  endif
              endif
1265      continue
          if (imatch.eq.0) then
              tnumang=tnumang+1
              tteq(tnumang)=teq(tictper(i))
              ttk(tnumang)=tk(tictper(i))
              tang1(tnumang)=tictper(i)
              tang2(tnumang)=1
              tictper(i)=tnumang
          else
              tictper(i)=imatch
          endif
          imatch=0
          do 1270 j=1,tnumang
              if (tictper(i+ngper).eq.tang1(j)) then
                  if (1.eq.tang2(j)) then
                      imatch=j
                  endif
              endif
1270      continue
          if (imatch.eq.0) then
              tnumang=tnumang+1
              tteq(tnumang)=teq(tictper(i+ngper))
              ttk(tnumang)=tk(tictper(i+ngper))
              tang1(tnumang)=tictper(i+ngper)
              tang2(tnumang)=1
              tictper(i+ngper)=tnumang
          else
              tictper(i+ngper)=imatch
          endif
c
c        
2320  continue
c
c
      write (6,'(a,i5,a,a,i5)') 'There were ',numang,' angle types, ',
     .            'now there are ',tnumang
c
c move these to the real arrays.
c
      numang=tnumang
      do 340 i=1,numang
          tk(i)=ttk(i)
          teq(i)=tteq(i)
340   continue
c
c now the torsions
c
      tnptra=0
      do 400 i=1,nphih
          imatch=0
          do 410 j=1,tnptra
c
c tphi1 is the comparison for type icph(), tphi2 is for scale factor tfach()
c
                if (icph(i).eq.tphi1(j)) then
                    if (tfach(i).eq.tphi2(j).or.
     &                        lespert.and.dscale.eq.3) then
c 
c same type and scale factor
c
                        imatch=j
                    endif
                endif
410         continue
            if (imatch.eq.0) then
c
c k is an offset in case there are multiple terms
c
                k=0
c
c save icph(i) for next term (if multi term) since we'll need it
c to find the next term's phase, etc. and icph(i) will get changed
c with first term
c
                l=icph(i)
415             tnptra=tnptra+1
                if (k.eq.0) then
c
c only set icph the first term (also don't use tphi1 and tphi2 for extra terms
c since we never want to match an second (or 3rd, etc) term)
c
                      tphi1(tnptra)=icph(i)
                      if (lespert.and.dscale.eq.3) then
                                tphi2(tnptra)=1
                      else
                                tphi2(tnptra)=tfach(i)
                      endif
                else
                      write (6,*) 'multi-term torsion',i,tnptra
                        tphi1(tnptra)=-1
                        tphi2(tnptra)=-1
                endif
                if (dscale.eq.3.and.lespert) then
                        tpk(tnptra)=pk(l+k)
                else
                        tpk(tnptra)=pk(l+k)/dble(tfach(i))
                endif
                tpn(tnptra)=pn(l+k)
                tphase(tnptra)=phase(l+k)
c
c set this torsion to this new type
c
                if (k.eq.0) then
                      icph(i)=tnptra
                endif
c
                if (tpn(tnptra).lt.0) then
c
c more than one term to the series, get the others
c
                        k=k+1
                        go to 415
                endif
            else
                icph(i)=imatch
            endif
C
c          if (.not.lespert.or.dscale.eq.3) go to 400
          imatch=0
          do 1275 j=1,tnptra
                if (ticph(i).eq.tphi1(j)) then
                    if (1.eq.tphi2(j)) then
                        imatch=j
                    endif
                endif
1275      continue
            if (imatch.eq.0) then
                k=0
                l=ticph(i)
1280            tnptra=tnptra+1
                if (k.eq.0) then
                        tphi1(tnptra)=ticph(i)
                        tphi2(tnptra)=1
                else
                        tphi1(tnptra)=-1
                        tphi2(tnptra)=-1
                endif
                tpn(tnptra)=pn(l+k)
                tpk(tnptra)=pk(l+k)
                tphase(tnptra)=phase(l+k)
                if (k.eq.0) then
                      ticph(i)=tnptra
                endif
                if (tpn(tnptra).lt.0) then
                        k=k+1
                        go to 1280
                endif
            else
                ticph(i)=imatch
            endif
c
400   continue
      do 420 i=1,nphia
          imatch=0
          do 430 j=1,tnptra
              if (icp(i).eq.tphi1(j)) then
                  if (tfac(i).eq.tphi2(j).or.
     &                        lespert.and.dscale.eq.3) then
                      imatch=j
                  endif
              endif
430       continue
          if (imatch.eq.0) then
                k=0
                l=icp(i)
435             tnptra=tnptra+1
                if (k.eq.0) then
                      tphi1(tnptra)=icp(i)
                    if (lespert.and.dscale.eq.3) then
                              tphi2(tnptra)=1
                    else
                              tphi2(tnptra)=tfac(i)
                    endif
                else  
                    write (6,*) 'multi-term torsion',i,tnptra
                      tphi1(tnptra)=-1
                      tphi2(tnptra)=-1
                endif
                if (lespert.and.dscale.eq.3) then
                          tpk(tnptra)=pk(l+k)
                else
                          tpk(tnptra)=pk(l+k)/dble(tfac(i))
                endif
                tpn(tnptra)=pn(l+k)
                tphase(tnptra)=phase(l+k)
                if (k.eq.0) icp(i)=tnptra
                if (tpn(tnptra).lt.0) then
                        k=k+1
                        go to 435
                endif
          else
                icp(i)=imatch
          endif
c
c          if (.not.lespert.or.dscale.eq.3) go to 420
          imatch=0
          do 1290 j=1,tnptra
                if (ticp(i).eq.tphi1(j)) then
                    if (1.eq.tphi2(j)) then
                        imatch=j
                    endif
                endif
1290        continue
            if (imatch.eq.0) then
                k=0
                l=ticp(i)
1295            tnptra=tnptra+1
                if (k.eq.0) then
                        tphi1(tnptra)=ticp(i)
                        tphi2(tnptra)=1
                else
                        tphi1(tnptra)=-1
                        tphi2(tnptra)=-1
                endif
                tpn(tnptra)=pn(l+k)
                tpk(tnptra)=pk(l+k)
                tphase(tnptra)=phase(l+k)
                if (k.eq.0) then
                      ticp(i)=tnptra
                endif
                if (tpn(tnptra).lt.0) then
                        k=k+1
                        go to 1295
                endif
            else
                ticp(i)=imatch
            endif
c
420   continue
c
      do 2420 i=1,ndper
          imatch=0
          do 2430 j=1,tnptra
              if (icpper(i).eq.tphi1(j)) then
                  if (tfacp(i).eq.tphi2(j).or.
     &                        lespert.and.dscale.eq.3) then
                      imatch=j
                  endif
              endif
2430      continue
          if (imatch.eq.0) then
                k=0
                l=icpper(i)
2435            tnptra=tnptra+1
                if (k.eq.0) then
                      tphi1(tnptra)=icpper(i)
                      if (lespert.and.dscale.eq.3) then
                                tphi2(tnptra)=1
                      else
                                tphi2(tnptra)=tfacp(i)
                      endif
                else  
                      tphi1(tnptra)=-1
                      tphi2(tnptra)=-1
                endif
                if (lespert.and.dscale.eq.3) then
                          tpk(tnptra)=pk(l+k)
                else
                          tpk(tnptra)=pk(l+k)/dble(tfacp(i))
                endif
                tpn(tnptra)=pn(l+k)
                tphase(tnptra)=phase(l+k)
                if (k.eq.0) icpper(i)=tnptra
                if (tpn(tnptra).lt.0) then
                        k=k+1
                        go to 2435
                endif
          else
                icpper(i)=imatch
          endif
          imatch=0
          do 2437 j=1,tnptra
              if (icpper(i+ndper).eq.tphi1(j)) then
                  if (tfacp(i).eq.tphi2(j).or.
     &                        lespert.and.dscale.eq.3) then
                      imatch=j
                  endif
              endif
2437      continue
          if (imatch.eq.0) then
                k=0
                l=icpper(i+ndper)
2440            tnptra=tnptra+1
                if (k.eq.0) then
                      tphi1(tnptra)=icpper(i+ndper)
                      if (lespert.and.dscale.eq.3) then
                                tphi2(tnptra)=1
                      else
                                tphi2(tnptra)=tfacp(i)
                      endif
                else  
                      tphi1(tnptra)=-1
                      tphi2(tnptra)=-1
                endif
                if (lespert.and.dscale.eq.3) then
                        tpk(tnptra)=pk(l+k)
                else
                        tpk(tnptra)=pk(l+k)/dble(tfacp(i))
                endif
                tpn(tnptra)=pn(l+k)
                tphase(tnptra)=phase(l+k)
                if (k.eq.0) icpper(i+ndper)=tnptra
                if (tpn(tnptra).lt.0) then
                          k=k+1
                          go to 2440
                endif
          else
                icpper(i+ndper)=imatch
          endif
c
c          if (.not.lespert) go to 2420
          imatch=0
          do 1300 j=1,tnptra
                if (ticpper(i).eq.tphi1(j)) then
                    if (1.eq.tphi2(j)) then
                        imatch=j
                    endif
                endif
1300      continue
          if (imatch.eq.0) then
                k=0
                l=ticpper(i)
1305            tnptra=tnptra+1
                if (k.eq.0) then
                        tphi1(tnptra)=ticpper(i)
                        tphi2(tnptra)=1
                else
                        tphi1(tnptra)=-1
                        tphi2(tnptra)=-1
                endif
                tpn(tnptra)=pn(l+k)
                tpk(tnptra)=pk(l+k)
                tphase(tnptra)=phase(l+k)
                if (k.eq.0) then
                      ticpper(i)=tnptra
                endif
                if (tpn(tnptra).lt.0) then
                        k=k+1
                        go to 1305
                endif
          else
                ticpper(i)=imatch
          endif
c
          imatch=0
          do 1310 j=1,tnptra
                if (ticpper(i+ndper).eq.tphi1(j)) then
                    if (1.eq.tphi2(j)) then
                        imatch=j
                    endif
                endif
1310      continue
          if (imatch.eq.0) then
                k=0
                l=ticpper(i+ndper)
1320            tnptra=tnptra+1
                if (k.eq.0) then
                        tphi1(tnptra)=ticpper(i+ndper)
                        tphi2(tnptra)=1
                else
                        tphi1(tnptra)=-1
                        tphi2(tnptra)=-1
                endif
                tpn(tnptra)=pn(l+k)
                tpk(tnptra)=pk(l+k)
                tphase(tnptra)=phase(l+k)
                if (k.eq.0) then
                      ticpper(i+ndper)=tnptra
                endif
                if (tpn(tnptra).lt.0) then
                        k=k+1
                        go to 1320
                endif
          else
                ticpper(i+ndper)=imatch
          endif
            
c
2420  continue
c
c now add dummy dihedrals- no barriers, but all else same. this is for 
c free rotation about torsions for LES copies when they 'do not exist' 
c (dummy copies) only do it if doing a LES perturbation. make a FULL
c second set of torsions although some might not be needed. that could
c be cleaned up but for now it is just fine.
c
c THIS IS NOT SUPPORTED UNTIL FURTHER TESTS ON CONVERGENCE 
c        if (lespert) then
c                realtors=tnptra
c                do 2600 i=1,realtors
c                   tnptra=tnptra+1
c                   tpk(tnptra)=0.d0
c                   tphase(tnptra)=tphase(i)
c                   tpn(tnptra)=tpn(i)
c2600                continue
c        endif
c
      write (6,'(a,i5,a,a,i5)') 'There were ',nptra,' dihedral types, ',
     .            'now there are ',tnptra
c
c move these to the real arrays.
c
      nptra=tnptra
      do 440 i=1,nptra
          pk(i)=tpk(i)
          pn(i)=tpn(i)
          phase(i)=tphase(i)
440   continue
c
      tntypes=0
      do 500 i=1,natom
c
c assume unique type
c
          imatch=0
          do 510 j=1,tntypes
c
c look for the type in list, both atom type and scale factor
c
              if (iac(i).eq.newiac(j)) then
                  if (iacfac(i).eq.newiacfac(j)) then
c
c set the pointer to this type
c
                      imatch=j
                  endif
              endif
510       continue
          if (imatch.eq.0) then
c
c add this as a new atom type
c
                tntypes = tntypes + 1
c
c set the pointers to tell us the orig type and scale factor
c to use for comparison of next atoms as well as creating new
c L-J terms
c
c example: atom 1 was type 3. now, it has 5 copies and therefore
c it needs a new type. let's say it is type 7.
c originally, iac(1)=3, iacfac(1)=5. now we change that.
c after changes, we woud have iac(1)=7 , newiac(7)=3, newiacfac(7)=5
c
                newiac(tntypes)=iac(i)
                newiacfac(tntypes)=iacfac(i)
                imatch=tntypes
          endif
c
c set the atom type pointer to this new type
c this means that the new atom type is for the new atom - sounds obvious but
c it is the one that gets used, not some pointer to the old atom type. 
c we do need to save the old one, updated for the new atom type list, for
c gibbs. do that down a bit...
c
          iac(i)=imatch
c
c
c look for original types if we are saving them to write les pert 
c topology files. just like bonds, etc above.
c
          if (.not.lespert) go to 1400
c
          imatch=0
          do 1410 j=1,tntypes
              if (tiac(i).eq.newiac(j)) then
                  if (1.eq.newiacfac(j)) then
                      imatch=j
                  endif
              endif
1410      continue
          if (imatch.eq.0) then
              tntypes = tntypes + 1
              newiac(tntypes)=tiac(i)
              newiacfac(tntypes)=1
              imatch=tntypes
          endif
          tiac(i)=imatch
c
c so tiac is what the atom was originally. need it for perturbing 
c 1 to many copies
c
c now look for the perturbed types (if they exist)
c
1400      if (ifpert.eq.0) go to 500
c
          imatch=0
          do 2510 j=1,tntypes
              if (iacper(i).eq.newiac(j)) then
                  if (iacpfac(i).eq.newiacfac(j)) then
                      imatch=j
                  endif
              endif
2510      continue
          if (imatch.eq.0) then
              tntypes = tntypes + 1
              newiac(tntypes)=iacper(i)
              newiacfac(tntypes)=iacpfac(i)
              imatch=tntypes
          endif
          iacper(i)=imatch
c
c now the original perturbed type, if we are saving it
c
1420      if (.not.lespert) go to 500
c
          imatch=0
          do 1430 j=1,tntypes
              if (tiacper(i).eq.newiac(j)) then
                  if (1.eq.newiacfac(j)) then
                      imatch=j
                  endif
              endif
1430      continue
          if (imatch.eq.0) then
              tntypes = tntypes + 1
              newiac(tntypes)=tiacper(i)
              newiacfac(tntypes)=1
              imatch=tntypes
          endif
          tiacper(i)=imatch
c
c next atom
c
500   continue
c
c now add a dummy type for disappearing atoms in LES- extra copies when
c not used
c
      if (lespert) then
          tntypes=tntypes+1
          dumiac=tntypes
      endif
c
      write (6,'(a,i5,a,a,i5)') 'There were ',ntypes,' atom types',
     .            'Now there are ',tntypes
c
c now the CN1 and CN2 lists. just use the index to pull the originals,
c scale them if needed, and store them. we'll make the index last
c
      hidx=0
      do 520 i=1,tntypes
          do 530 j=i,tntypes
c
c tntypes*(tntypes+1)/2 of these
c
c
c pull the original type for each atom, use it to get the CN1 (or CN2) index
c
c AMBER USES TWO WAYS TO GET CN1 AND CN2: ONE IS BY USING 
c IDX=ICO(NTYPES*(IAC(I)-1)+IAC(J)) AND THE OTHER IS TO IGNORE ICO AND
c SAY IDX=IAC(I)*(IAC(I)-1)/2+IAC(J) WHERE IAC(I)>IAC(J)
c HAVE TO SUPPORT BOTH!!
c
c if this involves trhe dummy type for LES perturbations, just set
c the interactions to zero.
c
            if (i.eq.dumiac.or.j.eq.dumiac) then
                    idx=j*(j-1)/2+i
                    tcn1(idx)=0
                    tcn2(idx)=0
                    tico(tntypes*(i-1)+j)=idx
                    tico(tntypes*(j-1)+i)=idx
                    go to 530
            endif
c normal case                
            k=ico(ntypes*(newiac(i)-1)+newiac(j))
            save=(k.gt.0)
            if (k.gt.0) then
c
c set idx to follow rule for new type numbers (j >= i here)
c
                idx=j*(j-1)/2+i
                tcn1(idx)=cn1(k)/dble(newiacfac(i))/dble(newiacfac(j))
                tcn2(idx)=cn2(k)/dble(newiacfac(i))/dble(newiacfac(j))
                tico(tntypes*(i-1)+j)=idx
                tico(tntypes*(j-1)+i)=idx
            else
c
c this is a 10-12 term not a 6-12, no funky rule here.
c
              hidx=hidx+1
              tasol(hidx)=asol(-k)/dble(newiacfac(i))/dble(newiacfac(j))
              tbsol(hidx)=bsol(-k)/dble(newiacfac(i))/dble(newiacfac(j))
              thbcut(hidx)=hbcut(-k)
              tico(tntypes*(i-1)+j)=-hidx
              tico(tntypes*(j-1)+i)=-hidx
c
c NOW APPARENTLY WE STILL NEED A 6-12 PARAMETER SET FOR THE 1-4 BETWEEN
c THIS PAIR, WHICH DOES NOT USE ICO() -SO THEREFORE THERE ARE NO 1-4
c HYDROGEN BONDS, AND VDW IS STILL CALCULATED. SO WE HAVE TO USE THE 'RULE'
c TO GET THE VDW PARAMS FOR THIS PAIR, CAN'T USE ICO, AND THEN USE
c THE RULE AGAIN TO PLACE THE PARAMS FOR THE 1-4'S IN THE RIGHT SPOT
c
              idx=j*(j-1)/2+i
c
c need new k since it was gotten from ico and it refers to hbond param here
c
              if (newiac(i).gt.newiac(j)) then        
                      k=newiac(i)*(newiac(i)-1)/2+newiac(j)
              else
                      k=newiac(j)*(newiac(j)-1)/2+newiac(i)
              endif
              tcn1(idx)=cn1(k)/dble(newiacfac(i))/dble(newiacfac(j))
              tcn2(idx)=cn2(k)/dble(newiacfac(i))/dble(newiacfac(j))
c
c leave tico as it was above, set from hidx!
c
            endif
530       continue
520   continue
c
c copy back to real arrays
c
      ntypes=tntypes
      do 540 i=1,ntypes*ntypes
              ico(i)=tico(i)
540   continue
      do 560 i=1,idx
          cn1(i)=tcn1(i)
          cn2(i)=tcn2(i)
560   continue
      nphb=hidx
      do 565 i=1,nphb
          asol(i)=tasol(i)
          bsol(i)=tbsol(i)
          hbcut(i)=thbcut(i)
565   continue
c
c done with L-J
c        do 38 i=1,natom
c       do 38 j=1,natom
c                k=ico(ntypes*(iac(i)-1)+iac(j))
c                if (k.gt.0) then
c                  write (11,36) i,j,' 6-12 ',cn1(k),cn2(k)
c                else
c                  write (11,36) i,j,'10-12 ',asol(-k),bsol(-k)
c36      format (2(1x,i4),a6,2(1x,f12.4))
c                endif
c
c                write (*,39) i,j,cn1(k),cn2(k)
c39      format (2(1x,i3),2(1x,f11.4))
c38      continue
c
c                
88        continue
c
c
c now on to the lestypes and mult matrix for nonbonded.
c
      nlestyp=0
      do 600 i=1,natom
c
c loop over current les types defined, look for match. assume match
c exists
c
          match=.false.
          do 610 j=1,nlestyp
c
c see if the particle has the same number les levels (number lesid's)
c as the les type, and if so, make sure all id's match. copy
c number irrelevant here.
c
c match2 tells whether THIS les type is appropriate
c match was for whether ANY type is appropriate
c
              match2=.true.
              if (numlev(j).eq.nlev(i)) then
                  do 620 k=1,nlev(i)
                      if (lesid(i,k).ne.typid(j,k)) then
                          match2=.false.
                      endif
620               continue
              else
                  match2=.false.
              endif
              if (match2) then
                  match=.true.
c
c give this particle this les type
c
                  lestyp(i)=j
c
c could issue a 'goto' here to skip the rest of the loop
c
                  go to 600
              endif
610       continue
c
c see if a match was found
c
          if (.not.match) then
c
c make a new les type
c
              nlestyp=nlestyp+1
              numlev(nlestyp)=nlev(i)
              do 630 k=1,nlev(i)
                  typid(nlestyp,k)=lesid(i,k)
630           continue
              lestyp(i)=nlestyp
          endif
600   continue
c
c so now all les types are assigned. need the mult factor matrix
c for the potential calculation. this is not trivially the
c diagonal elements only (scale up the same les type) since
c you need to scale up for any pair in the same les ID. this means
c some of the off-diagonal elements are non-zero, since different
c les types can share some les id's. this factor will be used to
c multiply ALL interactions when LES is used. not just the same
c les type! and one of the les types corresponds to 'no LES' (nlev=0).
c
      do 640 i=1,nlestyp
          do 650 j=1,nlestyp
             fac=1
             do 660 k=1,numlev(i)
               do 670 k2=1,numlev(j)
c
c see if this pair has ANY les id's in common. if so, multiply
c the scale factor
c
                  if (typid(i,k).eq.typid(j,k2)) then
                      fac=fac*dble(ncopies(typid(i,k)))
                  endif
670            continue
660          continue
             lesfac(nlestyp*(i-1)+j)=fac
650       continue
640   continue
c
c add a check to ensure that solvent molecules are not copied(for les
c only, we allow solvent molecules to be copied for PIMD)
c this is not supported in the current AMBER6 release, so we should
c not allow users to create the prmtop here.
c
c loop over solvent residues and make sure all atoms have
c nlev=0 (no LES levels defined)
c
      if( .not. pimd ) then 
         if (iptres.gt.0) then
            do i=iptres+1,nres
c
               if (i.eq.nres) then
                  j=natom
               else
                  j=ipres(i+1)-1
               endif
c
               do k=ipres(i),j
                  if (nlev(k).gt.0) then
                     write (6,*) "LES was applied to solvent, atom ",k
                     write (6,*) "This is not allowed."
                     stop
                  endif
               enddo
            enddo
         endif
      endif
c
c end checking for LES and solvent residues
c
c
c now calc the MAX number of atoms in a single residue- amber needs to 
c know this....
c
      nmxrs=0
      do 900 i=1,nres-1
              if (ipres(i+1)-ipres(i).gt.nmxrs) 
     &                        nmxrs=ipres(i+1)-ipres(i)
900   continue
c
c check last residue
c
      if (natom-(ipres(nres)-1).gt.nmxrs)
     &                 nmxrs=natom-(ipres(nres)-1)
        write (6,'(a,i7)') 'MAX # ATOMS IN A SINGLE RESIDUE = ',nmxrs
c
      write(6,*) ' ATOM  origpt'
!      do 1020 i=1,natom
!              write (6,1021) i,origpt(i)
!1021          format ('ATOM ',i6,' original atom # ',i6)
!1020  continue
c
c write extra topology files if requested for perturbations
c do this before stuff gets changed by writprm()
c
      if (lespert) then
              call lesprm(uw3prm)
              call les2prm(uw2prm)
      endif
c
c write new topology file
c
      call writprm(uwprm,uwnmr,wnmrok)
c
c
      if ((rcrdok.or.rcvdok.or.rcvbok.or.rcbdok).and.wcrdok) 
     &       call writcrd(uwcrd)
c
c        write (6,*) 'There are ',nspm,' molecules, nsp(i):'
c        write (6,800) (nsp(i),i=1,nspm)
c        write (6,*) 'mbper,mgper,mdper ',mbper,mgper,mgper
      write (6,*) 'Successful completion of ADDLES'
      end
