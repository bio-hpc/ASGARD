      subroutine readprm(urprm)
      
      implicit none

#include "SIZE.h"
#include "TOP.h"
#include "MISC.h"

      CHARACTER*80 FMT
      CHARACTER*80 FMTIN,IFMT,AFMT,RFMT,TYPE
      integer nf
      logical size_error
      integer IOK

      IFMT = '(12I6)'
      AFMT = '(20A4)'
      RFMT = '(5E16.8)'
c
      nf = urprm
      FMTIN = AFMT
      TYPE = 'TITLE'
      CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
      READ(NF,FMT) ITITL
c
c require that ADDLES uses new style prmtop files only
c this is a good idea but also needed because
c addles was converting old style prmtop to new ones
c and old prmtops do not have RADII values. If it writes
c a new style without radii, sander will have trouble
c because it assumes that the radii are correct
c if it is a new style prmtop
c
      if( iok == -1 ) then
         write(0,*) 'ADDLES now requires a new-style prmtop file'
         write(6,*) 'ADDLES now requires a new-style prmtop file'
         stop
      end if
c
      FMTIN = IFMT
      TYPE = 'POINTERS'
      CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
      read (urprm,fmt,err=1000) natom,ntypes,nbonh,mbona,ntheth,mtheta,
     &  nphih,mphia,nhparm,nparm,next,nres,nbona,ntheta,nphia,numbnd,
     &  numang,nptra,natyp,nphb,ifpert,nbper,ngper,ndper,mbper,mgper,
     &  mdper,ifbox,nmxrs,ifcap,numextra
c
c save orignal natom value
      natomo = natom
c
c check to make sure that nparm is zero- otherwise
c either LES was already used for this topology or else
c nparm is used- which it should not be.
c
      if (nparm.ne.0) then
              write (6,*) 'ERROR! NPARM IS NOT 0 ',nparm
              stop
      endif
c
c check size before actually reading the values
c
      call checksz( size_error )
      if (size_error) then
         close (urprm)
         stop
      endif

c
c read atom names
      FMTIN = AFMT
      TYPE = 'ATOM_NAME'
      CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
      read (urprm,fmt,err=1000) (igraph(i),i=1,natom)
c charges
      FMTIN = RFMT
      TYPE = 'CHARGE'
      CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
      read (urprm,fmt,err=1000) (chrg(i),i=1,natom)
c masses
      FMTIN = RFMT
      TYPE = 'MASS'
      CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
      read (urprm,fmt,err=1000) (amass(i),i=1,natom)
c perturbed masses (for TI, optional)
      FMTIN = RFMT
      TYPE = 'TI_MASS'
c JVAN: IONERR=1 (3rd argument to NXTSEC below) guarantees 
c       that TI_MASS section is optional.
      CALL NXTSEC(NF,  6,  1,FMTIN,  TYPE,  FMT,  IOK)
      if (IOK.EQ.0) then
         itimass=.TRUE.
         read (urprm,fmt,err=1000) (timass(i),i=1,natom)
      else 
         itimass=.FALSE.
      endif
c L-J index
      FMTIN = IFMT
      TYPE = 'ATOM_TYPE_INDEX'
      CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
      read (urprm,fmt,err=1000) (iac(i),i=1,natom)
c exclusions
      FMTIN = IFMT
      TYPE = 'NUMBER_EXCLUDED_ATOMS'
      CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
      read (urprm,fmt,err=1000) (numex(i),i=1,natom)
c intermediate L-J index
      FMTIN = IFMT
      TYPE = 'NONBONDED_PARM_INDEX'
      CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
      read (urprm,fmt,err=1000) (ico(i),i=1,ntypes*ntypes)
c residue labels
      FMTIN = AFMT
      TYPE = 'RESIDUE_LABEL'
      CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
      read (urprm,fmt,err=1000) (labres(i),i=1,nres)
c first atom in residue i
      FMTIN = IFMT
      TYPE = 'RESIDUE_POINTER'
      CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
      read (urprm,fmt,err=1000) (ipres(i),i=1,nres)
c bond force constants
      FMTIN = RFMT
      TYPE = 'BOND_FORCE_CONSTANT'
      CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
      read (urprm,fmt,err=1000) (rk(i),i=1,numbnd)
c bond equ length
      FMTIN = RFMT
      TYPE = 'BOND_EQUIL_VALUE'
      CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
      read (urprm,fmt,err=1000) (req(i),i=1,numbnd)
c angle force constants
      FMTIN = RFMT
      TYPE = 'ANGLE_FORCE_CONSTANT'
      CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
      read (urprm,fmt,err=1000) (tk(i),i=1,numang)
c angle equ value
      FMTIN = RFMT
      TYPE = 'ANGLE_EQUIL_VALUE'
      CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
      read (urprm,fmt,err=1000) (teq(i),i=1,numang)
c dihedral force constant
      FMTIN = RFMT
      TYPE = 'DIHEDRAL_FORCE_CONSTANT'
      CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
      read (urprm,fmt,err=1000) (pk(i),i=1,nptra)
c dihedral periodicity
      FMTIN = RFMT
      TYPE = 'DIHEDRAL_PERIODICITY'
      CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
      read (urprm,fmt,err=1000) (pn(i),i=1,nptra)
c dihedral phase
      FMTIN = RFMT
      TYPE = 'DIHEDRAL_PHASE'
      CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
      read (urprm,fmt,err=1000) (phase(i),i=1,nptra)
c solty unused
      FMTIN = RFMT
      TYPE = 'SOLTY'
      CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
      read (urprm,fmt,err=1000) (solty(i),i=1,natyp)
c L-J r**12 terms
      FMTIN = RFMT
      TYPE = 'LENNARD_JONES_ACOEF'
      CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
      read (urprm,fmt,err=1000) (cn1(i),i=1,ntypes*(ntypes+1)/2)
c L-J r**6 terms
      FMTIN = RFMT
      TYPE = 'LENNARD_JONES_BCOEF'
      CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
      read (urprm,fmt,err=1000) (cn2(i),i=1,ntypes*(ntypes+1)/2)
c atoms numers and index for bonds to hydrogen
      FMTIN = IFMT
      TYPE = 'BONDS_INC_HYDROGEN'
      CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
      read (urprm,fmt,err=1000) (ibh(i),jbh(i),icbh(i),i=1,nbonh)
      do 100 i=1,nbonh
        ibh(i)=ibh(i)/3+1
        jbh(i)=jbh(i)/3+1
 100  continue
c atom numbers and index for bonds not to hydrogen
      FMTIN = IFMT
      TYPE = 'BONDS_WITHOUT_HYDROGEN'
      CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
      read (urprm,fmt,err=1000) (ib(i),jb(i),icb(i),i=1,nbona)
      do 110 i=1,nbona
        ib(i)=ib(i)/3+1
        jb(i)=jb(i)/3+1
 110  continue
c angle info for hydrogens
      FMTIN = IFMT
      TYPE = 'ANGLES_INC_HYDROGEN'
      CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
      read (urprm,fmt,err=1000)(ith(i),jth(i),kth(i),icth(i),i=1,ntheth)
      do 120 i=1,ntheth
        ith(i)=ith(i)/3+1
        jth(i)=jth(i)/3+1
        kth(i)=kth(i)/3+1
 120  continue
c angles
      FMTIN = IFMT
      TYPE = 'ANGLES_WITHOUT_HYDROGEN'
      CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
      read (urprm,fmt,err=1000) (it(i),jt(i),kt(i),ict(i),i=1,ntheta)
      do 130 i=1,ntheta
        it(i)=it(i)/3+1
        jt(i)=jt(i)/3+1
        kt(i)=kt(i)/3+1
 130  continue
c
c dihedrals with hydrogens
c
      FMTIN = IFMT
      TYPE = 'DIHEDRALS_INC_HYDROGEN'
      CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
      read (urprm,fmt,err=1000) (iph(i),jph(i),kph(i),lph(i),icph(i)
     &  ,i=1,nphih)

      do 140 i=1,nphih
c
c save the negative info, it is for 1-4 and improper torsion flagging
c but make the atom numbers positive for now, then make negative again 
c when finished
c
        if (iph(i).lt.0) then
          iph(i)=abs(iph(i))
          inegh(i)=-1
        endif
        if (jph(i).lt.0) then
          jph(i)=abs(jph(i))
          jnegh(i)=-1
        endif
        if (kph(i).lt.0) then
          kph(i)=abs(kph(i))
          knegh(i)=-1
        endif
        if (lph(i).lt.0) then
          lph(i)=abs(lph(i))
          lnegh(i)=-1
        endif
        iph(i)=iph(i)/3+1
        jph(i)=jph(i)/3+1
        kph(i)=kph(i)/3+1
        lph(i)=lph(i)/3+1
 140  continue
c
c dihedrals
c
      FMTIN = IFMT
      TYPE = 'DIHEDRALS_WITHOUT_HYDROGEN'
      CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
      read (urprm,fmt,err=1000)(ip(i),jp(i),kp(i),lp(i),icp(i),
     &     i=1,nphia)
      do 150 i=1,nphia
        if (ip(i).lt.0) then
          ip(i)=abs(ip(i))
          ineg(i)=-1
        endif
        if (jp(i).lt.0) then
          jp(i)=abs(jp(i))
          jneg(i)=-1
        endif
        if (kp(i).lt.0) then
          kp(i)=abs(kp(i))
          kneg(i)=-1
        endif
        if (lp(i).lt.0) then
          lp(i)=abs(lp(i))
          lneg(i)=-1
        endif
        ip(i)=ip(i)/3+1
        jp(i)=jp(i)/3+1
        kp(i)=kp(i)/3+1
        lp(i)=lp(i)/3+1
 150  continue
c
c excluded atoms
c
      FMTIN = IFMT
      TYPE = 'EXCLUDED_ATOMS_LIST'
      CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
      read (urprm,fmt,err=1000) (natex(i),i=1,next)
c
c r**12 hbond term
c
      FMTIN = RFMT
      TYPE = 'HBOND_ACOEF'
      CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
      read (urprm,fmt,err=1000) (asol(i),i=1,nphb)
c
c r**10 term
c
      FMTIN = RFMT
      TYPE = 'HBOND_BCOEF'
      CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
      read (urprm,fmt,err=1000) (bsol(i),i=1,nphb)
c
c hbcut, not used
c
      FMTIN = RFMT
      TYPE = 'HBCUT'
      CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
      read (urprm,fmt,err=1000) (hbcut(i),i=1,nphb)
c
c amber atom types
c
      FMTIN = AFMT
      TYPE = 'AMBER_ATOM_TYPE'
      CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
      read (urprm,fmt,err=1000) (isymbl(i),i=1,natom)
c
c tree info types
c
      FMTIN = AFMT
      TYPE = 'TREE_CHAIN_CLASSIFICATION'
      CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
      read (urprm,fmt,err=1000) (itree(i),i=1,natom)
c
c tree joining info
c
      FMTIN = IFMT
      TYPE = 'JOIN_ARRAY'
      CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
      read (urprm,fmt,err=1000) (join(i),i=1,natom)
c
c last atom to move (not used)
c
      FMTIN = IFMT
      TYPE = 'IROTAT'
      CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
      read (urprm,fmt,err=1000) (irotat(i),i=1,natom)
c
c born radii and screening factor
c
c SRB: these require a new-style prmtop file;
c      if an old-style prmtop or { a new-style prmtop and
c      the flag was not found }, ie, iok not equal to 0,
c      then do Not attempt the read since it will fail.
 
      FMTIN = RFMT
      TYPE = 'RADII'
      CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
      if( iok .eq. 0 ) then
        read (urprm,fmt,err=1000) (rborn(i),i=1,natom)
      end if

      FMTIN = RFMT
      TYPE = 'SCREEN'
      CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
      if( iok .eq. 0 ) then
        read (urprm,fmt,err=1000) (fs(i),i=1,natom)
      end if
c
c only if ifbox.gt.0)
c
      if (ifbox.gt.0) then
        FMTIN = IFMT
        TYPE = 'SOLVENT_POINTERS'
        CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
        read (urprm,fmt,err=1000) iptres, nspm,nspsol
c
c check array bound
c
        if (nspm.gt.maxnspm) then
           write (6,*) 'Exceeded MAXNSPM'
           write (6,*) 'NSPM=',nspm
           write (6,*) 'MAXNSPM=',maxnspm
           stop
        endif
        FMTIN = IFMT
        TYPE = 'ATOMS_PER_MOLECULE'
        CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
        read (urprm,fmt,err=1000) (nsp(i),i=1,nspm)
        FMTIN = RFMT
        TYPE = 'BOX_DIMENSIONS'
        CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
        read (urprm,fmt,err=1000) beta,box(1),box(2),box(3)
      endif
c
c cap info
c
      if (ifcap.gt.0) then
        FMTIN = '(I6)'
        TYPE = 'CAP_INFO'
        CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
        read (urprm,fmt,err=1000) natcap
        FMTIN = '(4E16.8)'
        TYPE = 'CAP_INFO2'
        CALL NXTSEC(NF,  6,  0,FMTIN,  TYPE,  FMT,  IOK)
        read (urprm,fmt,err=1000) cutcap, xcap, ycap, zcap
      endif
c
c pert info
c
      if (ifpert.gt.0) then
        FMTIN = IFMT
        TYPE = 'PERT_BOND_ATOMS'
        CALL NXTSEC(NF,  6,  1,FMTIN,  TYPE,  FMT,  IOK)
        read (urprm,fmt,err=1000) (ibper(i),jbper(i),i=1,nbper)
        do 2100 i=1,nbper
          ibper(i)=ibper(i)/3+1
          jbper(i)=jbper(i)/3+1
2100  continue
        FMTIN = IFMT
        TYPE = 'PERT_BOND_PARAMS'
        CALL NXTSEC(NF,  6,  1,FMTIN,  TYPE,  FMT,  IOK)
        read (urprm,fmt,err=1000) (icbper(i),i=1,2*nbper)
C
        FMTIN = IFMT
        TYPE = 'PERT_ANGLE_ATOMS'
        CALL NXTSEC(NF,  6,  1,FMTIN,  TYPE,  FMT,  IOK)
        read (urprm,fmt,err=1000) (itper(i),jtper(i),ktper(i),
     &      i=1,ngper)
        do 2115 i=1,ngper
          itper(i)=itper(i)/3+1
          jtper(i)=jtper(i)/3+1
          ktper(i)=ktper(i)/3+1
2115  continue
        FMTIN = IFMT
        TYPE = 'PERT_ANGLE_PARAMS'
        CALL NXTSEC(NF,  6,  1,FMTIN,  TYPE,  FMT,  IOK)
        read (urprm,fmt,err=1000) (ictper(i),i=1,2*ngper)
C
        FMTIN = IFMT
        TYPE = 'PERT_DIHEDRAL_ATOMS'
        CALL NXTSEC(NF,  6,  1,FMTIN,  TYPE,  FMT,  IOK)
        read (urprm,fmt,err=1000) (ipper(i),jpper(i),kpper(i),
     &      lpper(i),i=1,ndper)
        do 2120 i=1,ndper
          if (ipper(i).lt.0) then
          ipper(i)=abs(ipper(i))
          inegp(i)=-1
          endif
          if (jpper(i).lt.0) then
          jpper(i)=abs(jpper(i))
          jnegp(i)=-1
          endif
          if (kpper(i).lt.0) then
          kpper(i)=abs(kpper(i))
          knegp(i)=-1
          endif
          if (lpper(i).lt.0) then
          lpper(i)=abs(lpper(i))
          lnegp(i)=-1
          endif
          ipper(i)=ipper(i)/3+1
          jpper(i)=jpper(i)/3+1
          kpper(i)=kpper(i)/3+1
          lpper(i)=lpper(i)/3+1
2120  continue
        FMTIN = IFMT
        TYPE = 'PERT_DIHEDRAL_PARAMS'
        CALL NXTSEC(NF,  6,  1,FMTIN,  TYPE,  FMT,  IOK)
        read (urprm,fmt,err=1000) (icpper(i),i=1,2*ndper)
C
        FMTIN = AFMT
        TYPE = 'PERT_RESIDUE_NAME'
        CALL NXTSEC(NF,  6,  1,FMTIN,  TYPE,  FMT,  IOK)
        read (urprm,fmt,err=1000) (labres(i),i=1,nres)
        FMTIN = AFMT
        TYPE = 'PERT_ATOM_NAME'
        CALL NXTSEC(NF,  6,  1,FMTIN,  TYPE,  FMT,  IOK)
        read (urprm,fmt,err=1000) (igrper(i),i=1,natom)
        FMTIN = AFMT
        TYPE = 'PERT_ATOM_SYMBOL'
        CALL NXTSEC(NF,  6,  1,FMTIN,  TYPE,  FMT,  IOK)
        read (urprm,fmt,err=1000) (ismper(i),i=1,natom)
        FMTIN = RFMT
        TYPE = 'ALMPER'
        CALL NXTSEC(NF,  6,  1,FMTIN,  TYPE,  FMT,  IOK)
        read (urprm,fmt,err=1000) (almper(i),i=1,natom)
        FMTIN = IFMT
        TYPE = 'IAPER'
        CALL NXTSEC(NF,  6,  1,FMTIN,  TYPE,  FMT,  IOK)
        read (urprm,fmt,err=1000) (iaper(i),i=1,natom)
        FMTIN = IFMT
        TYPE = 'PERT_ATOM_TYPE_INDEX'
        CALL NXTSEC(NF,  6,  1,FMTIN,  TYPE,  FMT,  IOK)
        read (urprm,fmt,err=1000) (iacper(i),i=1,natom)
        FMTIN = RFMT
        TYPE = 'PERT_CHARGE'
        CALL NXTSEC(NF,  6,  1,FMTIN,  TYPE,  FMT,  IOK)
        read (urprm,fmt,err=1000) (cgper(i),i=1,natom)
      endif
      if (ipol.eq.1) then
        write (6,*) 'ipol=1 not supported'
        stop
c      read (urprm,30,err=1000) (atpol(i),i=1,natom)
      endif
      if (ipol.eq.1.and.ifpert.eq.1) then
        write (6,*) 'ipol=1 not supported'
        stop
c      read (urprm,30,err=1000) (atpol1(i),i=1,natom)
      endif
c
c check array sizes in this initial topology file
c
      call checksz( size_error )
      if (size_error) then
         close (urprm)
         stop
      endif

      return
c
1000  write (6,*) 'problem reading parameter file'
      stop
      end
