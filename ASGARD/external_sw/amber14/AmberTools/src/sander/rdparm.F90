#include "copyright.h"
#include "../include/dprec.fh"
#include "../include/assert.fh"

#ifdef API
logical function need_to_bail(i, title, ierr)

   implicit none

   integer, intent(in) :: i
   character(len=*), intent(in) :: title
   integer, intent(out) :: ierr

   ierr = 0
   need_to_bail = .false.

   if (i == -2) then
      write(6, '(3a)') 'ERROR: Flag "', trim(title), '" not found in PARM file'
      need_to_bail = .true.
      ierr = 1
   end if

end function need_to_bail

#  define CHECK_REQUIRED(val, ttl) if (need_to_bail(val, ttl, ierr)) return
#else /* NOT API */

subroutine need_to_bail(i, title)

   implicit none

   integer, intent(in) :: i
   character(len=*), intent(in) :: title

   if (i == -2) then
      write(6, '(3a)') 'ERROR: Flag "', trim(title), '" not found in PARM file'
      call mexit(6, 1)
   end if
end subroutine need_to_bail

#  define CHECK_REQUIRED(val, ttl) call need_to_bail(val, ttl)
#endif /* API */

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Read the pointers (ONLY) from the topology file
#ifdef API
subroutine rdparm1(nf, ierr)
#else
subroutine rdparm1(nf)
#endif

   use parms, only : numbnd, numang, nptra, nphb, nttyp
   use charmm_mod, only : check_for_charmm
   use ff11_mod, only : check_cmap
   use file_io_dat
   use constants, only : RETIRED_INPUT_OPTION
   use amoeba_mdin,only:iamoeba

   implicit none
   
#  include "../lib/nxtsec.h"      
#  include "../include/md.h"
#  include "../include/memory.h"
#  include "box.h"
#  include "nmr.h"
#  include "extra_pts.h"
#  include "ew_cntrl.h"
   integer nf
   integer i,nspsol,iok
   integer nhparm
   integer mbper,mgper,mdper,mbona,mtheta,mphia ! read but ignored
   integer nlines, iscratch
   character(len=80) fmt,line
   character(len=80) fmtin,ifmt,afmt,rfmt,type
   character(len=80), allocatable, dimension(:) :: ffdesc
#ifdef API
   logical :: need_to_bail
   integer, intent(out) :: ierr
#else
   integer :: ierr
#endif
   ierr = 0

   initprmtop = .true.

   ifmt = '(12I6)'
   afmt = '(20A4)'
   rfmt = '(5E16.8)'
   
   !     ----- READ THE MOLECULAR TOPOLOGY -----
   
   nspsol = 0
   
   !     ----- FORMATTED INPUT -----

   ! Support both TITLE and CTITLE - CTITLE is used for a chamber 
   ! prmtop file. Essentially to prevent earlier versions of the code
   ! from loading such files.

   ! Search for CTITLE first
   fmtin = '(A80)'
   type = 'CTITLE'
   call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)

   if ( iok /= 0) then
     fmtin = '(A80)'
     type = 'TITLE'
     call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
     CHECK_REQUIRED(iok, 'TITLE')
   end if
   read(nf,fmtin) title
   
   fmtin = ifmt
   type = 'POINTERS'
   call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
   CHECK_REQUIRED(iok, 'POINTERS')
   read(nf,fmt) natom,ntypes,nbonh,mbona,ntheth,mtheta,nphih,mphia, &
         nhparm,nparm,nnb,nres,nbona,ntheta,nphia, &
         numbnd,numang,nptra,natyp,nphb,ifpert,nbper,ngper,ndper, &
         mbper,mgper,mdper,ifbox,nmxrs,ifcap,numextra,ncopy
   ! Back up our original numbnd value, since QM/MM simulations can increase
   ! this, and we need to know how much it increased by
#ifdef LES
   if (nparm /= 1) then
      write(6,*) ' *** THIS VERSION ONLY ACCEPTS TOPOLOGY FILES'
      write(6,*) '     THAT WERE CREATED BY ADDLES, WITH NPARM=1'
      write(6,*) '     USE A VERSION COMPILED WITHOUT -DLES '
#  ifdef API
      ierr = 1
      return
#  else
      call mexit(6,1)
#  endif /* API */
   end if
#else
   if (nparm == 1) then
      write(6,*) ' *** THIS VERSION WILL NOT ACCEPT TOPOLOGY FILES'
      write(6,*) '     THAT WERE CREATED BY ADDLES, WITH NPARM=1'
      write(6,*) '     USE A VERSION COMPILED WITH -DLES '
#  ifdef API
      ierr = 1
      return
#  else
      call mexit(6,1)
#  endif /* API */
   end if
#endif /* LES */
#ifndef API
   write(6,8118) natom,ntypes,nbonh,mbona,ntheth,mtheta,nphih,mphia, &
         nhparm,nparm,nnb,nres,nbona,ntheta,nphia,numbnd, &
         numang,nptra,natyp,nphb,ifbox,nmxrs,ifcap,numextra,ncopy
   8118 format(t2, &
         'NATOM  = ',i7,' NTYPES = ',i7,' NBONH = ',i7,' MBONA  = ',i7, &
         /' NTHETH = ',i7,' MTHETA = ',i7,' NPHIH = ',i7,' MPHIA  = ',i7, &
         /' NHPARM = ',i7,' NPARM  = ',i7,' NNB   = ',i7,' NRES   = ',i7, &
         /' NBONA  = ',i7,' NTHETA = ',i7,' NPHIA = ',i7,' NUMBND = ',i7, &
         /' NUMANG = ',i7,' NPTRA  = ',i7,' NATYP = ',i7,' NPHB   = ',i7, &
         /' IFBOX  = ',i7,' NMXRS  = ',i7,' IFCAP = ',i7,' NEXTRA = ',i7 &
         ,/' NCOPY  = ',i7/)
#endif /* API */
   
   nttyp = ntypes * (ntypes + 1) / 2 ! number of LJ type pairs
   
   ! Check that we don't have geometric constraints in the prmtop anymore

   if(nbona /= mbona .or. ntheta /= mtheta .or. nphia /= mphia) then
      write(6,*) 'Sander no longer allows constraints in prmtop'
      write(6,*) '...must have nbona=mbona, ntheta=mtheta, nphi=mphi'
      call mexit(6,1)
   end if

!! IPOL is now entered in prmtop. YD
!! If IPOL is set in prmtop, it replaces the one set in mdin, regardless.
   fmtin = ifmt
   type = 'IPOL'
   call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
   if (iok == 0) then
     read(nf,fmt) i
     if (ipol /= RETIRED_INPUT_OPTION .and. i /= ipol) &
         write(6,'(/,a,/,a,/,a)') 'Warning: ipol has been retired.', &
                                  '  ipol is replaced by the value in prmtop.'
     ipol = i
   end if
!! If not set in prmtop, it's set to 0 or the one entered from mdin
   if (ipol == RETIRED_INPUT_OPTION) then
     ipol = 0
   end if
! These two are actually used in the code.
   mpoltype = ipol
   induced  = ipol
   if( iamoeba > 0 ) induced = 1

#ifndef API
   ! Write implicit solvent radius and screening info to mdout
   if (( (igb /= 0 .or. ipb /= 0) .and. (ifcap == 0 .or. ifcap == 5)) &
                                  .or.hybridgb>0.or.icnstph.gt.1) then
      fmtin = afmt
      type = 'RADIUS_SET'
      call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
      if (iok == 0) then         ! Allow failure to support pre-AMBER9 prmtop
         read(nf,fmt) type      ! Reuse type var to avoid declaring a new one
         write(6,'(A,A)') ' Implicit solvent radii are ',type
      end if
      if ( igb == 7 ) then
         write(6,'(A)') ' Replacing prmtop screening parameters with GBn (igb=7) values'
      end if
      !Hai Nguyen: add igb = 8 here
      if ( igb == 8 ) then
         write(6,'(A)') ' Replacing prmtop screening parameters with GBn2 (igb=8) values'
      end if
   end if
#endif /* API */

   !    --- Read the force field information from the prmtop if available
   fmtin = '(i,a)'
   type = 'FORCE_FIELD_TYPE'
   call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
   if(iok == 0) then
     ! We found a force field description. Should be 1 or more lines
     ! of text which will be echoed back to the user. In some case
     ! we also take some decisions based on this.
#ifndef API
     write (6,'(a)') '| Force field information read from topology file: '
#endif
     ! format should be (nlines, text) where we read the first lines nlines
     ! to know how many more lines to read - we ignore nlines on subsequent
     ! lines but they still need to be there.
     ! e.g.
     ! 2 PARM99
     ! 2 GAFF
     read(nf,fmt) nlines,line
     allocate (ffdesc(nlines), stat = ierr)
     REQUIRE(ierr==0)
     ffdesc(1) = line
#ifndef API
     write(6,'(a,a)') '| ',ffdesc(1)
#endif
     do i = 2, nlines
       read(nf,fmt) iscratch,ffdesc(i)
#ifndef API
       write(6,'(a,a)') '| ',ffdesc(i)
#endif
     end do

     !Test to see if any special force fields are in use.
     
     !1) Check for CHARMM
     !Sets charmm_active = .true. if the charmm force field is in use.
     call check_for_charmm(nlines,ffdesc)


#ifndef API
     !End Test
     deallocate (ffdesc, stat = ierr)
     REQUIRE(ierr==0)
#endif /* API */

   end if 

   fmtin = '(i,a)'
   type = 'CMAP_COUNT'
   call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
! This is not necessary....
   if(iok == 0) then
#ifndef API
     write (6,'(a)') '| CMAP information read from topology file: '
#endif /* API */
     ! Check for CMAP only, sets cmap_active = .true. if cmap is used.
     call check_cmap !I should pass nf to validate the format.
   endif

   return
end subroutine rdparm1 

!=======================================================================

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Read in the rest of the prmtop file.
#ifdef API
subroutine rdparm2(x,ix,ih,nf,ierr)
#else
subroutine rdparm2(x,ix,ih,nf)
#endif
   use parms
   use molecule, only : mol_info
   use charmm_mod, only : charmm_active, read_charmm_params
   use ff11_mod, only : cmap_active, read_cmap_params
   use nblist, only: a,b,c
   use pimd_vars, only : dmdlm, itimass
#ifdef LES
   use les_data, only : lestyp, lestmp, lesfac, lfac, nlesty, maxlestyp, &
                        cnum, maxles, ileslst, subsp, nlesadj, maxlesadj, &
                        jleslst
   use pimd_vars, only : ipimd
   use amoeba_mdin,only: iamoeba
#endif
   use file_io_dat
#ifdef DSSP
   use dssp, only: ipepc, npepc
#endif
   implicit none
   integer nf
   _REAL_ x(*)
   integer ix(*)
   character(len=4) ih(*)
   
#  include "../lib/nxtsec.h"
#  include "../include/md.h"
#  include "../include/memory.h"
#  include "box.h"
#  include "nmr.h"
#ifdef LES
   integer iexcl,numex,k1,j1
#endif
#  include "ew_mpole.h"
   integer ntype,i,iok
   integer j,jj,k,l,n,nn
   integer iptres,nspsol,natsm
   _REAL_ oldbeta,duma,dumb,dumc,dumd
   integer mol_atm_cnt
   _REAL_ massdiff   ! = mass[perturbed] - mass[original] for TI w.r.t. mass.
   integer irotat( natom )  ! dummy: disappears when leaving subroutine


   integer :: res_start, res_end

   integer ier ! Allocation status.

   character(len=80) fmt
   character(len=80) fmtin,ifmt,afmt,rfmt,type
#ifdef API
   logical :: need_to_bail
   integer, intent(out) :: ierr
   ierr = 0
#endif

   ifmt = '(12I6)'
   afmt = '(20A4)'
   rfmt = '(5E16.8)'
   ntype = ntypes*ntypes
   
   !     ----- READ THE SYMBOLS AND THE CHARGES AND THE MASSES -----
   
   fmtin = afmt
   type = 'ATOM_NAME'
   call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
   CHECK_REQUIRED(iok, 'ATOM_NAME')
   read(nf,fmt) (ih(m04+i-1),i = 1,natom)
   
   fmtin = rfmt
   type = 'CHARGE'
   call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
   CHECK_REQUIRED(iok, 'CHARGE')
   read(nf,fmt) (x(l15+i-1),i = 1,natom)
   
   fmtin = ifmt
   type = 'ATOMIC_NUMBER'
   call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
   if (iok == 0) then
      read(nf,fmt) (ix(i100+i),i = 1,natom)
      ix(i100)=1
   else
      ix(i100)=0
#ifndef API
      if (igb .eq. 7 .or. igb .eq. 8 .or. gbsa .eq. 1) then
         write(6,'(a)') &
           '| Warning: ATOMIC_NUMBER section not found'
         write(6,'(a)') &
           '|          Guessing atomic numbers from masses for GB parameters'
         write(6,'(a)') &
           '|          Remake topology file with AmberTools 12 or add ATOMIC_NUMBERS'
         write(6,'(a)') &
           '|          with ParmEd to remove this warning'
      end if
#endif /* API */
   end if
   
   fmtin = rfmt
   type = 'MASS'
   call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
   CHECK_REQUIRED(iok, 'MASS')
   read(nf,fmt) (x(lwinv+i-1),i = 1,natom)
   
   fmtin = ifmt
   type = 'ATOM_TYPE_INDEX'
   call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
   CHECK_REQUIRED(iok, 'ATOM_TYPE_INDEX')
   read(nf,fmt) (ix(i04+i-1),i = 1,natom)
   
   fmtin = ifmt
   type = 'NUMBER_EXCLUDED_ATOMS'
   call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
   CHECK_REQUIRED(iok, 'NUMBER_EXCLUDED_ATOMS')
   read(nf,fmt) (ix(i+i08-1),i = 1,natom)
   
   fmtin = ifmt
   type = 'NONBONDED_PARM_INDEX'
   call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
   CHECK_REQUIRED(iok, 'NONBONDED_PARM_INDEX')
   read(nf,fmt) (ix(i+i06-1),i = 1,ntype)
   
   fmtin = afmt
   type = 'RESIDUE_LABEL'
   call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
   CHECK_REQUIRED(iok, 'RESIDUE_LABEL')
   read(nf,fmt) (ih(i+m02-1),i=1,nres)
   
   fmtin = ifmt
   type = 'RESIDUE_POINTER'
   call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
   CHECK_REQUIRED(iok, 'RESIDUE_POINTER')
   read(nf,fmt) (ix(i+i02-1),i=1,nres)
   ix(i02+nres) = natom+1
 
   !     ----- READ THE PARAMETERS -----
   
   fmtin = rfmt
   type = 'BOND_FORCE_CONSTANT'
   call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
   CHECK_REQUIRED(iok, 'BOND_FORCE_CONSTANT')
   read(nf,fmt) (rk(i),    i = 1,numbnd)
   
   fmtin = rfmt
   type = 'BOND_EQUIL_VALUE'
   call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
   CHECK_REQUIRED(iok, 'BOND_EQUIL_VALUE')
   read(nf,fmt) (req(i),   i = 1,numbnd)
   
   fmtin = rfmt
   type = 'ANGLE_FORCE_CONSTANT'
   call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
   CHECK_REQUIRED(iok, 'ANGLE_FORCE_CONSTANT')
   read(nf,fmt) (tk(i),    i = 1,numang)
   
   fmtin = rfmt
   type = 'ANGLE_EQUIL_VALUE'
   call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
   CHECK_REQUIRED(iok, 'ANGLE_EQUIL_VALUE')
   read(nf,fmt) (teq(i),   i = 1,numang)
   
   fmtin = rfmt
   type = 'DIHEDRAL_FORCE_CONSTANT'
   call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
   CHECK_REQUIRED(iok, 'DIHEDRAL_FORCE_CONSTANT')
   read(nf,fmt) (pk(i),    i = 1,nptra)
   
   fmtin = rfmt
   type = 'DIHEDRAL_PERIODICITY'
   call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
   CHECK_REQUIRED(iok, 'DIHEDRAL_PERIODICITY')
   read(nf,fmt) (pn(i),    i = 1,nptra)
   
   fmtin = rfmt
   type = 'DIHEDRAL_PHASE'
   call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
   CHECK_REQUIRED(iok, 'DIHEDRAL_PHASE')
   read(nf,fmt) (phase(i), i = 1,nptra)

   !   ----- READ VARIABLE SCEE AND SCNB VALUES IF THEY EXIST -----
   ! --- SCEE ---
   fmtin = rfmt
   type = 'SCEE_SCALE_FACTOR'
   call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
   if (iok == 0) then
     !we found the SCEE scale factor data so read it in, there should be
     !a value for each dihedral type. Note while there is one for each
     !dihedral type not all are actually used in the calculation since
     !1-4's are only done over unique dihedral terms. If multiple dihedrals
     !for a given set of atom types exist, with different pn's for example
     !then the last value of SCEE/SCNB should be used.
#ifndef API
     write (6,'(a)') ''
     write (6,'(a)') '| Note: 1-4 EEL scale factors are being read from the topology file.'
#endif

     read(nf,fmt) (one_scee(i), i = 1,nptra)
     ! We used to use vdinv here but this could in principal cause an issue since
     ! nptra includes improper dihedrals which have 0.0d0 in the prmtop file.
     ! call vdinv(nptra,one_scee,one_scee)
     do i = 1, nptra
       if (one_scee(i) /= 0.0d0) then
         one_scee(i) = 1.0d0/one_scee(i)
       !else it gets left at 0.0d0.
       end if
     end do
   else
     !We will use default scee of 1.2
#ifndef API
     write (6,'(a)') ''
     write (6,'(a)') '| Note: 1-4 EEL scale factors were NOT found in the topology file.'
     write (6,'(a)') '|       Using default value of 1.2.'
#endif
     do i=1,nptra
       one_scee(i)=1.0d0/1.2d0
     end do 
   end if !if iok==0

   ! --- SCNB ---
   fmtin = rfmt
   type = 'SCNB_SCALE_FACTOR'
   call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
   if (iok == 0) then
#ifndef API
     write (6,'(a)') ''
     write (6,'(a)') '| Note: 1-4 VDW scale factors are being read from the topology file.'
#endif

     read(nf,fmt) (one_scnb(i), i = 1,nptra)
     ! We used to use vdinv here but this could in principal cause an issue
     ! since
     ! nptra includes improper dihedrals which have 0.0d0 in the prmtop file.
     ! call vdinv(nptra,one_scnb,one_scnb)
     do i = 1, nptra
       if (one_scnb(i) /= 0.0d0) then
         one_scnb(i) = 1.0d0/one_scnb(i)
       !else it gets left at 0.0d0.
       end if
     end do
   else
     !We will use default scnb of 2.0d0
#ifndef API
     write (6,'(a)') ''
     write (6,'(a)') '| Note: 1-4 VDW scale factors were NOT found in the topology file.'
     write (6,'(a)') '|       Using default value of 2.0.'
#endif
     do i=1,nptra
       one_scnb(i)=1.0d0/2.0d0
     end do 
   end if !if iok==0
   
   !   ----- END READ VARIABLE SCEE AND SCNB VALUES IF THEY EXIST -----

   fmtin = rfmt
   type = 'SOLTY'
   call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
   CHECK_REQUIRED(iok, 'SOLTY')
   read(nf,fmt) (solty(i), i = 1,natyp)
   
   fmtin = rfmt
   type = 'LENNARD_JONES_ACOEF'
   call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
   CHECK_REQUIRED(iok, 'LENNARD_JONES_ACOEF')
   read(nf,fmt) (cn1(i),   i = 1,nttyp)
   
   fmtin = rfmt
   type = 'LENNARD_JONES_BCOEF'
   call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
   CHECK_REQUIRED(iok, 'LENNARD_JONES_BCOEF')
   read(nf,fmt) (cn2(i),   i = 1,nttyp)

   if (lj1264 /= 0) then
      fmtin=rfmt
      type = 'LENNARD_JONES_CCOEF'
      call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
      CHECK_REQUIRED(iok, 'LENNARD_JONES_CCOEF')

      if (iok /= 0) then
         write(6, '(a)') 'Could not find LENNARD_JONES_CCOEF in topology &
                      &file! Recreate topology for lj1264.'
         call mexit(6,1)
      endif

      read(nf, fmt) (cn6(i), i= 1,nttyp)

   end if

   if ( vdwmodel == 1 ) then
      fmtin = rfmt
      type = 'EXPVDWMODEL_BETA'
      call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
      CHECK_REQUIRED(iok, 'EXPVDWMODEL_BETA')
      read(nf,fmt) (cn3(i),   i = 1,nttyp)
      fmtin = rfmt
      type = 'EXPVDWMODEL_A'
      call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
      CHECK_REQUIRED(iok, 'EXPVDWMODEL_A')
      read(nf,fmt) (cn4(i),   i = 1,nttyp)
      fmtin = rfmt
      type = 'EXPVDWMODEL_B'
      call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
      CHECK_REQUIRED(iok, 'EXPVDWMODEL_B')
      read(nf,fmt) (cn5(i),   i = 1,nttyp)
   endif
   
   !     ----- READ THE BONDING INFORMATIONS -----
   
   fmtin = ifmt
   type = 'BONDS_INC_HYDROGEN'
   call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
   CHECK_REQUIRED(iok, 'BONDS_INC_HYDROGEN')
   read(nf,fmt) (ix(i+iibh-1),ix(i+ijbh-1),ix(i+iicbh-1), &
         i = 1,nbonh)
   
   fmtin = ifmt
   type = 'BONDS_WITHOUT_HYDROGEN'
   call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
   CHECK_REQUIRED(iok, 'BONDS_WITHOUT_HYDROGEN')
   read(nf,fmt)(ix(i+iiba-1),ix(i+ijba-1),ix(i+iicba-1),i = 1,nbona)
   
   fmtin = ifmt
   type = 'ANGLES_INC_HYDROGEN'
   call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
   CHECK_REQUIRED(iok, 'ANGLES_INC_HYDROGEN')
   read(nf,fmt) (ix(i+i24-1),ix(i+i26-1),ix(i+i28-1),ix(i+i30-1), &
         i = 1,ntheth)
   
   fmtin = ifmt
   type = 'ANGLES_WITHOUT_HYDROGEN'
   call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
   CHECK_REQUIRED(iok, 'ANGLES_WITHOUT_HYDROGEN')
   read(nf,fmt) (ix(i+i32-1),ix(i+i34-1),ix(i+i36-1),ix(i+i38-1), &
         i = 1,ntheta)
   
   fmtin = ifmt
   type = 'DIHEDRALS_INC_HYDROGEN'
   call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
   CHECK_REQUIRED(iok, 'DIHEDRALS_INC_HYDROGEN')
   read(nf,fmt) (ix(i+i40-1),ix(i+i42-1),ix(i+i44-1),ix(i+i46-1), &
         ix(i+i48-1),i = 1,nphih)
   
   fmtin = ifmt
   type = 'DIHEDRALS_WITHOUT_HYDROGEN'
   call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
   CHECK_REQUIRED(iok, 'DIHEDRALS_WITHOUT_HYDROGEN')
   read(nf,fmt) (ix(i+i50-1),ix(i+i52-1),ix(i+i54-1),ix(i+i56-1), &
         ix(i+i58-1),i = 1,nphia)
   
   fmtin = ifmt
   type = 'EXCLUDED_ATOMS_LIST'
   call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
   CHECK_REQUIRED(iok, 'EXCLUDED_ATOMS_LIST')
   read(nf,fmt) (ix(i+i10-1),i=1,nnb)
   
   !     ----- READ THE H-BOND PARAMETERS -----
   
   fmtin = rfmt
   type = 'HBOND_ACOEF'
   call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
   CHECK_REQUIRED(iok, 'HBOND_ACOEF')
   read(nf,fmt) (asol(i),i=1,nphb)
   
#ifndef HAS_10_12
   do i=1,nphb
      if( asol(i) /= 0.d0 ) then
         write(6,*) 'Found a non-zero 10-12 coefficient, but source', &
               ' was not compiled with -DHAS_10_12.'
         write(6,*) 'If you are using a pre-1994 force field, you', &
               ' will need to re-compile with this flag.'
         call mexit(6,1)
      end if
   end do
#endif
   
   fmtin = rfmt
   type = 'HBOND_BCOEF'
   call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
   CHECK_REQUIRED(iok, 'HBOND_BCOEF')
   read(nf,fmt) (bsol(i),i=1,nphb)
   
   fmtin = rfmt
   type = 'HBCUT'
   call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
   CHECK_REQUIRED(iok, 'HBCUT')
   read(nf,fmt) (hbcut(i),i=1,nphb)
   
   !     ----- READ ISYMBL,ITREE,JOIN AND IROTAT ARRAYS -----
   
   fmtin = afmt
   type = 'AMBER_ATOM_TYPE'
   call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
   CHECK_REQUIRED(iok, 'AMBER_ATOM_TYPE')
   read(nf,fmt) (ih(i+m06-1),i=1,natom)
   
   fmtin = afmt
   type = 'TREE_CHAIN_CLASSIFICATION'
   call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
   CHECK_REQUIRED(iok, 'TREE_CHAIN_CLASSIFICATION')
   read(nf,fmt) (ih(i+m08-1),i=1,natom)
   
   fmtin = ifmt
   type = 'JOIN_ARRAY'
   call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
   CHECK_REQUIRED(iok, 'JOIN_ARRAY')
   read(nf,fmt) (ix(i+i64-1),i=1,natom)
   
   fmtin = ifmt
   type = 'IROTAT'
   call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
   CHECK_REQUIRED(iok, 'IROTAT')
   read(nf,fmt) (irotat(i), i=1,natom)
   
   !     ----- READ THE BOUNDARY CONDITION STUFF -----
   
   nspm = 1
   ix(i70) = natom
   if (ifbox > 0) then
      
      fmtin = ifmt
      type = 'SOLVENT_POINTERS'
      call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
      CHECK_REQUIRED(iok, 'SOLVENT_POINTERS')
      read(nf,fmt) iptres,nspm,nspsol
      
      fmtin = ifmt
      type = 'ATOMS_PER_MOLECULE'
      call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
      CHECK_REQUIRED(iok, 'ATOMS_PER_MOLECULE')
      read(nf,fmt) (ix(i+i70-1),i=1,nspm)
      
      ! Box dimensions in the prmtop file are redundant. For new format prmtop files
      ! we simply ignore them. For old format files we have to read them for compatibility
      ! but they are not used.
      fmtin = rfmt
      type = 'BOX_DIMENSIONS'
      call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
      !-1 = amber 6 prmtop file.
      if (iok==-1) read(nf,fmt) oldbeta,duma,dumb,dumc
      
      if( igb /= 0  .or. ipb /= 0 .or.  ntb == 0 )then
         box(1)=0.0d0
         box(2)=0.0d0
         box(3)=0.0d0
      else
         box(1)=a
         box(2)=b
         box(3)=c
      end if
      
   end if  ! (ifbox > 0)
   
   !     ----- LOAD THE CAP INFORMATION IF NEEDED -----
   
   if(ifcap == 1) then
      fmtin = '(I6)'
      type = 'CAP_INFO'
      call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
      CHECK_REQUIRED(iok, 'CAP_INFO')
      read(nf,fmt) natcap
      
      fmtin = '(4E16.8)'
      type = 'CAP_INFO2'
      call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
      CHECK_REQUIRED(iok, 'CAP_INFO2')
      read(nf,fmt) cutcap,xcap,ycap,zcap
   end if
   
   if( (igb /= 0 .or. ipb /= 0) .and. ifcap == 0 .and. iok == -1 ) then
      write(0,*) 'GB/PB calculations now require a new-style prmtop file'
      write(6,*) 'GB/PB calculations now require a new-style prmtop file'
      call mexit(6,1)
   end if
   
   if ((igb /= 0 .or. ipb /= 0) .and. ifcap /= 0 .and. ifcap /= 5) then
      write(0,*) 'GB/PB calculations are incompatible with spherical solvent caps'
      write(6,*) 'GB/PB calculations are incompatible with spherical solvent caps'
      call mexit(6,1)
   end if

   if (( (igb /= 0 .or. ipb /= 0) .and. (ifcap == 0 .or. ifcap == 5)) &
                                  .or.hybridgb>0.or.icnstph.gt.1) then
      fmtin = rfmt
      type = 'RADII'
      call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
      CHECK_REQUIRED(iok, 'RADII')
      read(nf,fmt) (x(l97+i-1),i=1,natom)
      type = 'SCREEN'
      call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
      CHECK_REQUIRED(iok, 'SCREEN')
      read(nf,fmt) (x(l96+i-1),i=1,natom)
   end if
   
!
! IPOL is now specified in prmtop YD
!
   if (ipol > 0) then
      if (igb /= 0 .or. ipb /= 0) then
         write(0,*) 'GB/PB calculations are incompatible with polarizable force fields'
         write(6,*) 'GB/PB calculations are incompatible with polarizable force fields'
         call mexit(6,1)
      end if
      fmtin = rfmt
      type = 'POLARIZABILITY'
      call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
      CHECK_REQUIRED(iok, 'POLARIZABILITY')
      read(nf,fmt) (x(lpol+i-1),i=1,natom)
   end if

! Modified by WJM
   if (ipol > 1) then
      fmtin = rfmt
      type = 'DIPOLE_DAMP_FACTOR'
      call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
      if ( iok .eq. 0 ) then
         read(nf,fmt) (x(ldf+i-1),i=1,natom)
      end if
   end if
#ifdef API
   call load_ewald_pol_info(ipol, x(lpol), x(lpol2), x(ldf), natom)
#else
   call load_ewald_pol_info(ipol, x(lpol), x(lpol2), x(ldf), natom, iok)
#endif

   ! Check that every atom is assigned to a molecule for NTP simulations. If
   ! not, segfaults or chaos will ensue.

   if (ntp .gt. 0) then

      mol_atm_cnt = 0

      do i = 1, nspm
         mol_atm_cnt = mol_atm_cnt + ix(i70+i-1)
      end do

      if (mol_atm_cnt .ne. natom) then
         write(6,'(a)') 'Error: Bad topology file. Sum of ATOMS_PER_MOLECULE &
                        &does not equal NATOM.'
         call mexit(6, 1)
      end if

   end if

   !     ----- READ THE PERTURBED MASSES IF NEEDED  -----
   if (itimass > 0) then
      ! JVAN: dmdlm must be allocated here, not in pimd_init.
      allocate( dmdlm(1:natom),stat=ier )
      REQUIRE( ier == 0 )
      fmtin = rfmt
      type = 'TI_MASS'
      call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
      CHECK_REQUIRED(iok, 'TI_MASS')
      read(nf,fmt) (dmdlm(i) ,i = 1,natom)
      do i=1,natom
         massdiff = (dmdlm(i) - x(lwinv+i-1))
         x(lwinv+i-1) = x(lwinv+i-1) + clambda * massdiff
         dmdlm(i) = massdiff/x(lwinv+i-1)
      end do
   end if
   
#ifdef DSSP
   !   ----- construct an array containing the atom numbers of the carbon atoms of all peptide
   !         groups
   allocate( ipepc(1:natom),stat=ier )
   REQUIRE( ier == 0 )
   k = 0
   do i=1,natom
      if( ih(m04+i-1) == 'C   ' ) then
         k = k+1
         ipepc(k) = i
      end if
   end do
   npepc = k
#endif
#ifdef LES

   if (nparm == 1.and.iamoeba.eq.0) then
      fmtin = ifmt
      type = 'LES_NTYP'
      call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
      CHECK_REQUIRED(iok, 'LES_NTYP')
      read (nf,fmt) nlesty
      lestmp=nlesty*nlesty
      
      ! check the array sizes to make sure we do not overflow
      
      ! LES types
      
      if (nlesty > maxlestyp) then
         write (6,*) 'Exceeded MAXLESTYP',nlesty
         stop
      end if
      
      ! LES atoms
      
      if (natom > maxles) then
         write (6,*) 'Exceeded MAXLES',natom
         stop
      end if

      fmtin = ifmt
      type = 'LES_TYPE'
      call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
      CHECK_REQUIRED(iok, 'LES_TYPE')
      read (nf,fmt) (lestyp(i),i=1,natom)
      fmtin = rfmt
      type = 'LES_FAC'
      call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
      CHECK_REQUIRED(iok, 'LES_FAC')
      read (nf,fmt) (lesfac(i),i=1,lestmp)
      fmtin = ifmt
      type = 'LES_CNUM'
      call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
      CHECK_REQUIRED(iok, 'LES_CNUM')
      read (nf,fmt) (cnum(i), i=1,natom)
      fmtin = ifmt
      type = 'LES_ID'
      call nxtsec(nf,  6,  1,fmtin,  type,  fmt,  iok)
      CHECK_REQUIRED(iok, 'LES_ID')
      read (nf,fmt) (subsp(i), i=1,natom)
#ifndef API
      write (6, '(a)') ' LES parameters were found'
#endif
      
      ! now create the list of atoms that have non-unitary scaling factors.
      ! this will be used in the Ewald calculation to correct for the
      ! lack of use of the intra-copy scaling factor in the charge grid.
      ! all of these pairs will need correction. The list will not change and
      ! is therefore only calculated once (here).
      
      nlesadj=0
      
      iexcl=0
      
      ! pairs are listed in two arrays, for i and j, rather than using
      ! a set of pointers like the nonbond and exclusion lists. This is 
      ! since many atoms will not have any correction partners (since 
      ! they are not in LES).

      if( ipimd.eq.0 ) then
         ! pimd are not going to use nb_adjust_les, so we do not need to generate
         ! les adjust list
         !
         do k=1,natom

            lestmp=nlesty*(lestyp(k)-1)
            !         write (6,*) 'atom1 : ',k,lestmp
            
            ! need to sum all f the number of exclusions even if non-LES atoms.
            ! see below.

            numex=ix(k+i08-1)

            DO_J: do j=k+1,natom
               lfac=lesfac(lestmp+lestyp(j))

               ! check for non-zero scaling factor (meaning a correction will 
               ! be required)

               if (abs(lfac-1.0d0) > 0.01) then

                  ! check to make sure these aren't excluded atoms (since then 
                  ! no correction is wanted)

                  !  FORMAT(12I6)  (NATEX(i), i=1,NEXT)
                  !  the excluded atom list.  To get the excluded list for atom
                  !  "i" you need to traverse the NUMEX list, adding up all
                  !  the previous NUMEX values, since NUMEX(i) holds the number
                  !  of excluded atoms for atom "i", not the index into the
                  !  NATEX list.  Let IEXCL = SUM(NUMEX(j), j=1,i-1), then
                  !  excluded atoms are NATEX(IEXCL)+1 to NATEX(IEXCL+NUMEX(i)).


                  do k1=iexcl+1,iexcl+numex

                     ! get exclusion

                     j1=ix(k1+i10-1)

                     ! check against atom j

                     if (j1 == j) then
                        ! excluded, get next atom j
                        ! write (6,*) 'Exclusion list match'
                        cycle DO_J
                     end if

                     ! check next entry

                  end do

                  ! if we arrived here, the atom was not in the exclusion list
                  ! so this pair will need correction
                  ! (should add boundary checking for variables here)

                  if (nlesadj == maxlesadj) then
                     write (6,*) 'EXCEEDED MAXLESADJ!'
                     stop
                  end if

                  nlesadj=nlesadj+1
                  ileslst(nlesadj)=k
                  jleslst(nlesadj)=j
               end if

               ! next j

            end do DO_J

            ! increment the exclusion list pointer for atom i
            iexcl=iexcl+numex

            ! next i

         end do

      end if !(ipimd == 0 )

#ifndef API
      write (6,6520) nlesadj
      6520 format (1x,i7,' LES atom pairs require adjustment')
#endif
      
      ! end creation of LES adjustment list and reading LES info

   end if  ! (nparm==1)

#endif /* LES */
   
   !     ----- CALCULATE INVERSE, TOTAL MASSES -----
   !       -- save the masses for removal of mass weighted velocity,
   !          leaving the inverse masses in the legacy, Lwinv area
   
   tmass = 0.0d0
   !     -- index over molecules
   j = l75-1
   jj = i70-1
   !     -- index over mass->invmass
   k = lwinv-1
   !     -- index over saved mass
   l = lmass-1
   do n = 1,nspm
      j = j + 1
      jj = jj + 1
      x(j) = 0.0d0
      natsm = ix(jj)
      do nn = 1,natsm
         k = k+1
         l = l+1
         
         ! -- sum molecule
         x(j) = x(j) + x(k)
         
         ! -- save mass in "new" Lmass area
         x(l) = x(k)
         
         ! -- make inverse in "old" Lwinv area
         if( x(k) /= 0.d0 ) x(k) = 1.0d0 / x(k)
      end do
      tmass = tmass + x(j)
   end do
   tmassinv = 1.0d0 / tmass

 
   !     ----- SCALE THE CHARGES IF DIELC.GT.1.0E0 -----
   
   if (dielc /= 1.0e0 .and. igb == 0 .and. ipb == 0) then
      dumd = sqrt(dielc)
      do i = 1,natom
         x(i+l15-1) = x(i+l15-1)/dumd
      end do
   end if
   
   !     ----- INVERT THE HBCUT ARRAY -----
   
   do i = 1,nphb
      if(hbcut(i) <= 0.001e0) hbcut(i) = 1.0d-10
      hbcut(i) = 1.0e0/hbcut(i)
   end do
   
   !     ----- duplicate dihedral pointers for vector ephi -----
   
   call dihdup(nphih,ix(i40),ix(i42),ix(i44),ix(i46),ix(i48),pn)
   call dihdup(nphia,ix(i50),ix(i52),ix(i54),ix(i56),ix(i58),pn)
   
   !     --- pre-calculate some parameters for vector ephi ---
   
   call dihpar(nptra,pk,pn,phase,gamc,gams,ipn)
   
   if (charmm_active) then
     call read_charmm_params(nf)
   end if
   if (cmap_active) then
     call read_cmap_params(nf)
   end if
   
   ! -----------------------------
   ! Fill module 'molecule' arrays
   ! -----------------------------
   !
   ! Atomic masses
   do i=1,natom
     mol_info%atom_mass(i) = x(lwinv+i-1)
   end do
   
   ! atom_to_resid_map provides residue number 
   ! given a specific atom number.
   do i = 1, nres
     ! natom_res
     mol_info%natom_res(i) = ix(i+i70-1)

     ! atom_to_resid_map
     res_start = ix(i+i02-1)
     res_end = ix(i+i02) - 1
     do j = res_start, res_end
       mol_info%atom_to_resid_map(j) = i
     end do
   end do
   
   return
end subroutine rdparm2

#undef CHECK_REQUIRED
