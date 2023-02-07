      subroutine writprm(uwprm,uwnmr,wnmrok)

#include "SIZE.h"
#include "TOP.h"
#include "MISC.h"
      integer cnum(MAXNATOM),count
      integer invcnt(MAXNATOM), invcop(MAXNATOM,MAXCOPY)

10    format (20a4)
15    format (a)
20    format (10i8)
30    format (5e16.8)

      write(uwprm,'(a)') '%VERSION  VERSION_STAMP = V0001.000  '
      write(uwprm,'(a)') '%FLAG TITLE'
      write( uwprm,'(a)') '%FORMAT(20a4)'
      write (uwprm,'(a)',err=1000) ititl
      nparm = 1
      write( uwprm,'(a)') '%FLAG POINTERS'
      write( uwprm,'(a)') '%FORMAT(10I8)'
      write (uwprm,20,err=1000) natom,ntypes,nbonh,mbona,ntheth,mtheta,
     &    nphih,mphia,nhparm,nparm,next,nres,nbona,ntheta,nphia,numbnd,
     &    numang,nptra,natyp,nphb,ifpert,nbper,ngper,ndper,mbper,mgper,
     &    mdper,ifbox,nmxrs,ifcap,numextra,numcop
c write atom names
      write( uwprm,'(a)') '%FLAG ATOM_NAME'
      write( uwprm,'(a)') '%FORMAT(20a4)'
      write (uwprm,10,err=1000) (igraph(i),i=1,natom)
c charges
      write( uwprm,'(a)') '%FLAG CHARGE'
      write( uwprm,'(a)') '%FORMAT(5E16.8)'
      write (uwprm,30,err=1000) (chrg(i),i=1,natom)
c masses
      if (bigmas) then
        do 200 i=1,natom
          if (imass(origpt(i)).ne.0) amass(i)=newm
200      continue
      endif    
      write( uwprm,'(a)') '%FLAG MASS'
      write( uwprm,'(a)') '%FORMAT(5E16.8)'
      write (uwprm,30,err=1000) (amass(i),i=1,natom)
c perturbed masses (for TI, optional)
      if (itimass) then
         write( uwprm,'(a)') '%FLAG TI_MASS'
         write( uwprm,'(a)') '%FORMAT(5E16.8)'
         write (uwprm,30,err=1000) (timass(i),i=1,natom)
      endif
c L-J index
      write( uwprm,'(a)') '%FLAG ATOM_TYPE_INDEX'
      write( uwprm,'(a)') '%FORMAT(10I8)'
      write (uwprm,20,err=1000) (iac(i),i=1,natom)
c exclusions
      write( uwprm,'(a)') '%FLAG NUMBER_EXCLUDED_ATOMS'
      write( uwprm,'(a)') '%FORMAT(10I8)'
      write (uwprm,20,err=1000) (numex(i),i=1,natom)
c another L-J index?
      write( uwprm,'(a)') '%FLAG NONBONDED_PARM_INDEX'
      write( uwprm,'(a)') '%FORMAT(10I8)'
      write (uwprm,20,err=1000) (ico(i),i=1,ntypes*ntypes)
c residue labels
      write( uwprm,'(a)') '%FLAG RESIDUE_LABEL'
      write( uwprm,'(a)') '%FORMAT(20a4)'
      write (uwprm,10,err=1000) (labres(i),i=1,nres)
c first atom in residue i
      write( uwprm,'(a)') '%FLAG RESIDUE_POINTER'
      write( uwprm,'(a)') '%FORMAT(10I8)'
      write (uwprm,20,err=1000) (ipres(i),i=1,nres)
c bond force constants
      write( uwprm,'(a)') '%FLAG BOND_FORCE_CONSTANT'
      write( uwprm,'(a)') '%FORMAT(5E16.8)'
      write (uwprm,30,err=1000) (rk(i),i=1,numbnd)
c bond equ length
      write( uwprm,'(a)') '%FLAG BOND_EQUIL_VALUE'
      write( uwprm,'(a)') '%FORMAT(5E16.8)'
      write (uwprm,30,err=1000) (req(i),i=1,numbnd)
c angle force constants
      write( uwprm,'(a)') '%FLAG ANGLE_FORCE_CONSTANT'
      write( uwprm,'(a)') '%FORMAT(5E16.8)'
      write (uwprm,30,err=1000) (tk(i),i=1,numang)
c angle equ value
      write( uwprm,'(a)') '%FLAG ANGLE_EQUIL_VALUE'
      write( uwprm,'(a)') '%FORMAT(5E16.8)'
      write (uwprm,30,err=1000) (teq(i),i=1,numang)
c dihedral force constant
      write( uwprm,'(a)') '%FLAG DIHEDRAL_FORCE_CONSTANT'
      write( uwprm,'(a)') '%FORMAT(5E16.8)'
      write (uwprm,30,err=1000) (pk(i),i=1,nptra)
c dihedral periodicity
      write( uwprm,'(a)') '%FLAG DIHEDRAL_PERIODICITY'
      write( uwprm,'(a)') '%FORMAT(5E16.8)'
      write (uwprm,30,err=1000) (pn(i),i=1,nptra)
c dihedral phase
      write( uwprm,'(a)') '%FLAG DIHEDRAL_PHASE'
      write( uwprm,'(a)') '%FORMAT(5E16.8)'
      write (uwprm,30,err=1000) (phase(i),i=1,nptra)
c solty unused
      write( uwprm,'(a)') '%FLAG SOLTY'
      write( uwprm,'(a)') '%FORMAT(5E16.8)'
      write (uwprm,30,err=1000) (solty(i),i=1,natyp)
c L-J r**12 terms
      write( uwprm,'(a)') '%FLAG LENNARD_JONES_ACOEF'
      write( uwprm,'(a)') '%FORMAT(5E16.8)'
      write (uwprm,30,err=1000) (cn1(i),i=1,ntypes*(ntypes+1)/2)
c L-J r**6 terms
      write( uwprm,'(a)') '%FLAG LENNARD_JONES_BCOEF'
      write( uwprm,'(a)') '%FORMAT(5E16.8)'
      write (uwprm,30,err=1000) (cn2(i),i=1,ntypes*(ntypes+1)/2)
c atoms numbers and index for bonds to hydrogen
      do 100 i=1,nbonh
          ibh(i)=(ibh(i)-1)*3
          jbh(i)=(jbh(i)-1)*3
100    continue
      write( uwprm,'(a)') '%FLAG BONDS_INC_HYDROGEN'
      write( uwprm,'(a)') '%FORMAT(10I8)'
      write (uwprm,20,err=1000) (ibh(i),jbh(i),icbh(i),i=1,nbonh)
c atom numbers and index for bonds not to hydrogen
      do 110 i=1,nbona
          ib(i)=(ib(i)-1)*3
          jb(i)=(jb(i)-1)*3
110    continue
      write( uwprm,'(a)') '%FLAG BONDS_WITHOUT_HYDROGEN'
      write( uwprm,'(a)') '%FORMAT(10I8)'
      write (uwprm,20,err=1000) (ib(i),jb(i),icb(i),i=1,nbona)
c angle info for hydrogens
      do 120 i=1,ntheth
          ith(i)=(ith(i)-1)*3
          jth(i)=(jth(i)-1)*3
          kth(i)=(kth(i)-1)*3
120    continue
      write( uwprm,'(a)') '%FLAG ANGLES_INC_HYDROGEN'
      write( uwprm,'(a)') '%FORMAT(10I8)'
      write (uwprm,20,err=1000) (ith(i),jth(i),kth(i),icth(i),
     .        i=1,ntheth)
c angles
      do 130 i=1,ntheta
          it(i)=(it(i)-1)*3
          jt(i)=(jt(i)-1)*3
          kt(i)=(kt(i)-1)*3
130   continue
      write( uwprm,'(a)') '%FLAG ANGLES_WITHOUT_HYDROGEN'
      write( uwprm,'(a)') '%FORMAT(10I8)'
      write (uwprm,20,err=1000) (it(i),jt(i),kt(i),ict(i),i=1,ntheta)
c dihedrals with hydrogens
      do 140 i=1,nphih
          iph(i)=(iph(i)-1)*3
          jph(i)=(jph(i)-1)*3
          kph(i)=(kph(i)-1)*3
          lph(i)=(lph(i)-1)*3
          if (inegh(i).lt.0) iph(i)=-(iph(i))
          if (jnegh(i).lt.0) jph(i)=-(jph(i))
          if (knegh(i).lt.0) kph(i)=-(kph(i))
          if (lnegh(i).lt.0) lph(i)=-(lph(i))
140     continue
      write( uwprm,'(a)') '%FLAG DIHEDRALS_INC_HYDROGEN'
      write( uwprm,'(a)') '%FORMAT(10I8)'
      write (uwprm,20,err=1000) (iph(i),jph(i),kph(i),lph(i),icph(i)
     &    ,i=1,nphih)
c dihedrals
      do 150 i=1,nphia
          ip(i)=(ip(i)-1)*3
          jp(i)=(jp(i)-1)*3
          kp(i)=(kp(i)-1)*3
          lp(i)=(lp(i)-1)*3
          if (ineg(i).lt.0) ip(i)=-(ip(i))
          if (jneg(i).lt.0) jp(i)=-(jp(i))
          if (kneg(i).lt.0) kp(i)=-(kp(i))
          if (lneg(i).lt.0) lp(i)=-(lp(i))
150     continue
      write( uwprm,'(a)') '%FLAG DIHEDRALS_WITHOUT_HYDROGEN'
      write( uwprm,'(a)') '%FORMAT(10I8)'
      write (uwprm,20,err=1000)(ip(i),jp(i),kp(i),lp(i),icp(i),
     &       i=1,nphia)
c excluded atoms
      write( uwprm,'(a)') '%FLAG EXCLUDED_ATOMS_LIST'
      write( uwprm,'(a)') '%FORMAT(10I8)'
      write (uwprm,20,err=1000) (natex(i),i=1,next)
c r**12 hbond term
      write( uwprm,'(a)') '%FLAG HBOND_ACOEF'
      write( uwprm,'(a)') '%FORMAT(5E16.8)'
      write (uwprm,30,err=1000) (asol(i),i=1,nphb)
c r**10 term
      write( uwprm,'(a)') '%FLAG HBOND_BCOEF'
      write( uwprm,'(a)') '%FORMAT(5E16.8)'
      write (uwprm,30,err=1000) (bsol(i),i=1,nphb)
c hbcut, not used
      write( uwprm,'(a)') '%FLAG HBCUT'
      write( uwprm,'(a)') '%FORMAT(5E16.8)'
      write (uwprm,30,err=1000) (hbcut(i),i=1,nphb)
c amber atom types
      write( uwprm,'(a)') '%FLAG AMBER_ATOM_TYPE'
      write( uwprm,'(a)') '%FORMAT(20a4)'
      write (uwprm,10,err=1000) (isymbl(i),i=1,natom)
c tree info types
      write( uwprm,'(a)') '%FLAG TREE_CHAIN_CLASSIFICATION'
      write( uwprm,'(a)') '%FORMAT(20a4)'
      write (uwprm,10,err=1000) (itree(i),i=1,natom)
c tree joining info
      write( uwprm,'(a)') '%FLAG JOIN_ARRAY'
      write( uwprm,'(a)') '%FORMAT(10I8)'
      write (uwprm,20,err=1000) (join(i),i=1,natom)
c last atom to move (not used)
      write( uwprm,'(a)') '%FLAG IROTAT'
      write( uwprm,'(a)') '%FORMAT(10I8)'
      write (uwprm,20,err=1000) (irotat(i),i=1,natom)
c born radii
      write( uwprm,'(a)') '%FLAG RADII'
      write( uwprm,'(a)') '%FORMAT(5E16.8)'
      write (uwprm,30,err=1000) (rborn(i),i=1,natom)
c screening factor
      write( uwprm,'(a)') '%FLAG SCREEN'
      write( uwprm,'(a)') '%FORMAT(5E16.8)'
      write (uwprm,30,err=1000) (fs(i),i=1,natom)
c only if ibox.gt.0)
      if (ifbox.gt.0) then
        write( uwprm,'(a)') '%FLAG SOLVENT_POINTERS'
        write( uwprm,'(a)') '%FORMAT(10I8)'
        write (uwprm,20,err=1000) iptres, nspm,nspsol
        write( uwprm,'(a)') '%FLAG ATOMS_PER_MOLECULE'
        write( uwprm,'(a)') '%FORMAT(10I8)'
        write (uwprm,20,err=1000) (nsp(i),i=1,nspm)
        write( uwprm,'(a)') '%FLAG BOX_DIMENSIONS'
        write( uwprm,'(a)') '%FORMAT(5E16.8)'
        write (uwprm,30,err=1000) beta,box(1),box(2),box(3)
      endif
c cap info
      if (ifcap.gt.0) then
        write( uwprm,'(a)') '%FLAG CAP_INFO'
        write( uwprm,'(a)') '%FORMAT(10I8)'
        write (uwprm,20,err=1000) natcap
        write( uwprm,'(a)') '%FLAG CAP_INFO2'
        write( uwprm,'(a)') '%FORMAT(5E16.8)'
        write (uwprm,30,err=1000) cutcap, xcap, ycap, zcap
      endif
c pert info
      if (ifpert.gt.0) then
          do 2100 i=1,nbper
            ibper(i)=(ibper(i)-1)*3
            jbper(i)=(jbper(i)-1)*3
2100      continue
          write( uwprm,'(a)') '%FLAG PERT_BOND_ATOMS'
          write( uwprm,'(a)') '%FORMAT(10I8)'
          write (uwprm,20,err=1000) (ibper(i),jbper(i),i=1,nbper)
          write( uwprm,'(a)') '%FLAG PERT_BOND_PARAMS'
          write (uwprm,20,err=1000) (icbper(i),i=1,2*nbper)
          do 2110 i=1,ngper
            itper(i)=(itper(i)-1)*3
            jtper(i)=(jtper(i)-1)*3
            ktper(i)=(ktper(i)-1)*3
2110      continue
          write( uwprm,'(a)') '%FLAG PERT_ANGLE_ATOMS'
          write( uwprm,'(a)') '%FORMAT(10I8)'
          write (uwprm,20,err=1000) (itper(i),jtper(i),ktper(i),
     &        i=1,ngper)
          write( uwprm,'(a)') '%FLAG PERT_BOND_PARAMS'
          write( uwprm,'(a)') '%FORMAT(10I8)'
          write (uwprm,20,err=1000) (ictper(i),i=1,2*ngper)
          do 2120 i=1,ndper
                ipper(i)=(ipper(i)-1)*3
                jpper(i)=(jpper(i)-1)*3
                kpper(i)=(kpper(i)-1)*3
                lpper(i)=(lpper(i)-1)*3
                if (inegp(i).lt.0) then
                        ipper(i)=-(ipper(i))
                endif
                if (jnegp(i).lt.0) then
                        jpper(i)=-(jpper(i))
                endif
                if (knegp(i).lt.0) then
                        kpper(i)=-(kpper(i))
                endif
                if (lnegp(i).lt.0) then
                        lpper(i)=-(lpper(i))
                endif
2120      continue
          write( uwprm,'(a)') '%FLAG PERT_DIHEDRAL_ATOMS'
          write( uwprm,'(a)') '%FORMAT(10I8)'
          write (uwprm,20,err=1000) (ipper(i),jpper(i),kpper(i),
     &        lpper(i),i=1,ndper)
          write( uwprm,'(a)') '%FLAG PERT_DIHEDRAL_PARAMS'
          write (uwprm,20,err=1000) (icpper(i),i=1,2*ndper)
          write( uwprm,'(a)') '%FLAG PERT_RESIDUE_NAME'
          write( uwprm,'(a)') '%FORMAT(20A4)'
          write (uwprm,10,err=1000) (labres(i),i=1,nres)
          write( uwprm,'(a)') '%FLAG PERT_ATOM_NAME'
          write( uwprm,'(a)') '%FORMAT(20A4)'
          write (uwprm,10,err=1000) (igrper(i),i=1,natom)
          write( uwprm,'(a)') '%FLAG PERT_ATOM_SYMBOL'
          write( uwprm,'(a)') '%FORMAT(20A4)'
          write (uwprm,10,err=1000) (ismper(i),i=1,natom)
          write( uwprm,'(a)') '%FLAG ALMPER'
          write( uwprm,'(a)') '%FORMAT(5E16.8)'
          write (uwprm,30,err=1000) (almper(i),i=1,natom)
          write( uwprm,'(a)') '%FLAG IAPER'
          write( uwprm,'(a)') '%FORMAT(10I8)'
          write (uwprm,20,err=1000) (iaper(i),i=1,natom)
          write( uwprm,'(a)') '%FLAG PERT_ATOM_TYPE_INDEX'
          write (uwprm,20,err=1000) (iacper(i),i=1,natom)
          write( uwprm,'(a)') '%FLAG PERT_CHARGE'
          write( uwprm,'(a)') '%FORMAT(5E16.8)'
          write (uwprm,30,err=1000) (cgper(i),i=1,natom)
      endif
      if (ipol.eq.1) then
          stop
c        write (uwprm,30,err=1000) (atpol(i),i=1,natom)
      endif
      if (ipol.eq.1.and.ifpert.eq.1) then
          stop
c        write (uwprm,30,err=1000) (atpol1(i),i=1,natom)
      endif
      write( uwprm,'(a)') '%FLAG LES_NTYP'
      write( uwprm,'(a)') '%FORMAT(10I8)'
      write (uwprm,20,err=1000) nlestyp
      write( uwprm,'(a)') '%FLAG LES_TYPE'
      write( uwprm,'(a)') '%FORMAT(10I8)'
      write (uwprm,20,err=1000) (lestyp(i),i=1,natom)
      write( uwprm,'(a)') '%FLAG LES_FAC'
      write( uwprm,'(a)') '%FORMAT(5E16.8)'
      write (uwprm,30,err=1000) (lesfac(i),i=1,nlestyp*nlestyp)
c
c now write a copy number for each atom
c
        count=1
        do 1990 i=1,natom
          invcnt(i) = 0
1990    continue

        if(numcop .eq. 1) then
           do i=1,natom
             cnum(i) = nrem(i)
           enddo
        else
           do 2000 i=1,natom
             if (totcop(origpt(i)).eq.1) then
               cnum(i)=0
             else
               cnum(i)=count
               count=count+1
             endif
             if (count.gt.totcop(origpt(i))) count=1
             if( .not.atm1st ) then
                invcnt(origpt(i)) = invcnt(origpt(i)) + 1
                invcop(origpt(i),invcnt(origpt(i))) = i
             end if
2000       continue
        endif
c
      write( uwprm,'(a)') '%FLAG LES_CNUM'
      write( uwprm,'(a)') '%FORMAT(10I8)'
      write (uwprm,20,err=1000) (cnum(i),i=1,natom)
c
c write out the subspace to which each particle belongs, warn
c if particle is in more than 1 subspace (nlev>1) 
c nlev=0 for non-LES atoms
c for particles where nlev(i)=0, lesid(i,1) has been set to zero
c by default.
c
c check to make sure each atom is in only a single subspace for now
c
       do 2014 i=1,natom
           if (nlev(i).gt.1) then
               write (6,*) 'WARNING: atom ',i
               write (6,*) 'is in more than 1 LES subspace!'
               write (6,*) 'NMR restraints may not be supported.'
           endif
2014   continue
c
      write( uwprm,'(a)') '%FLAG LES_ID'
      write( uwprm,'(a)') '%FORMAT(10I8)'
      write (uwprm,20,err=1000) (lesid(i,1),i=1,natom)
c
c  write out material to process NMR restraints:
c
        if( wnmrok) then
c
c         get an "inverse" mapping for subsequent maniupalation of
c           NMR constraints, etc.
c           invcnt  contains the number of copies, indexed by original atom #
c           invcop  contains the new atom numbers
c
          write (uwnmr,20,err=1000) natomo
          do 2010 i=1,natomo
            write( uwprm,'(a)') '%FLAG LES_INVERSE'
            write( uwprm,'(a)') '%FORMAT(10I8)'
            write (uwnmr,20,err=1000) invcnt(i),
     .                  (invcop(i,j),j=1,invcnt(i))
2010      continue
        end if
      return
c
1000  write (6,*) 'problem writing parameter file'
      stop
      end
