!OPTIMIZEME??
#include "copyright.h"
#include "../include/dprec.fh"
#include "../include/assert.fh"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Read the pointers (ONLY) from the topology file
subroutine cpparm1(parmdata, ierr)

   use amoeba_mdin, only : iamoeba
   use charmm_mod, only : charmm_active
   use constants, only : RETIRED_INPUT_OPTION
   use file_io_dat
   use ff11_mod, only : check_cmap
   use prmtop_type, only : prmtop_struct
   use parms, only : numbnd, numang, nptra, nphb, nttyp

   implicit none
   
#  include "../lib/nxtsec.h"      
#  include "../include/md.h"
#  include "../include/memory.h"
#  include "box.h"
#  include "nmr.h"
#  include "extra_pts.h"
#  include "ew_cntrl.h"
   type(prmtop_struct), intent(in) :: parmdata
   integer nspsol, nhparm
   integer mbper,mgper,mdper,mbona,mtheta,mphia ! read but ignored
   character(len=80) ifmt,afmt,rfmt
   integer, intent(out) :: ierr
   ierr = 0

   initprmtop = .true.

   ifmt = '(12I6)'
   afmt = '(20A4)'
   rfmt = '(5E16.8)'
   
   !     ----- READ THE MOLECULAR TOPOLOGY -----
   
   nspsol = 0
   
   !     ----- FORMATTED INPUT -----

   title = parmdata%title
   natom = parmdata%natom
   ntypes = parmdata%ntypes
   nbonh = parmdata%nbonh
   mbona = parmdata%mbona
   ntheth = parmdata%ntheth
   mtheta = parmdata%mtheta
   nphih = parmdata%nphih
   mphia = parmdata%mphia
   nhparm = parmdata%nhparm
   nparm = parmdata%nparm
   nnb = parmdata%nnb
   nres = parmdata%nres
   nbona = parmdata%nbona
   ntheta = parmdata%ntheta
   nphia = parmdata%nphia
   numbnd = parmdata%numbnd
   numang = parmdata%numang
   nptra = parmdata%nptra
   natyp = parmdata%natyp
   nphb = parmdata%nphb
   ifpert = parmdata%ifpert
   nbper = parmdata%nbper
   ngper = parmdata%ngper
   ndper = parmdata%ndper
   mbper = parmdata%mbper
   mgper = parmdata%mgper
   mdper = parmdata%mdper
   ifbox = parmdata%ifbox
   nmxrs = parmdata%nmxrs
   ifcap = parmdata%ifcap
   numextra = parmdata%numextra
   ncopy = parmdata%ncopy

   ! Back up our original numbnd value, since QM/MM simulations can increase
   ! this, and we need to know how much it increased by
#ifdef LES
   if (nparm /= 1) then
      write(6,*) ' *** THIS VERSION ONLY ACCEPTS TOPOLOGY FILES'
      write(6,*) '     THAT WERE CREATED BY ADDLES, WITH NPARM=1'
      write(6,*) '     USE A VERSION COMPILED WITHOUT -DLES '
      ierr = 1
      return
   end if
#else
   if (nparm == 1) then
      write(6,*) ' *** THIS VERSION WILL NOT ACCEPT TOPOLOGY FILES'
      write(6,*) '     THAT WERE CREATED BY ADDLES, WITH NPARM=1'
      write(6,*) '     USE A VERSION COMPILED WITH -DLES '
      ierr = 1
      return
   end if
#endif
   
   nttyp = ntypes * (ntypes + 1) / 2 ! number of LJ type pairs
   
   ! Check that we don't have geometric constraints in the prmtop anymore

   if(nbona /= mbona .or. ntheta /= mtheta .or. nphia /= mphia) then
      write(6,*) 'Sander no longer allows constraints in prmtop'
      write(6,*) '...must have nbona=mbona, ntheta=mtheta, nphi=mphi'
      ierr = 1
      return
   end if

   ipol = parmdata%ipol

! These two are actually used in the code.
   mpoltype = ipol
   induced  = ipol
   if( iamoeba > 0 ) induced = 1

   charmm_active = parmdata%is_chamber /= 0

   return
end subroutine cpparm1 

!=======================================================================

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Read in the rest of the prmtop file.
subroutine cpparm2(x,ix,ih,parmdata,ierr)
!  use charmm_mod, only : charmm_active, read_charmm_params2
#ifdef DSSP
   use dssp, only: ipepc, npepc
#endif
!  use ff11_mod, only : cmap_active
   use molecule, only : mol_info
   use nblist, only : a,b,c
   use parms
   use pimd_vars, only : dmdlm, itimass
   use prmtop_type, only : prmtop_struct
#ifdef LES
   use amoeba_mdin, only : iamoeba
   use pimd_vars, only : ipimd
   use les_data, only : lestyp, lesfac, cnum, subsp, nlesty, nlesadj, &
                        maxlestyp, maxles, lfac, lesfac, lestmp, ileslst, &
                        jleslst, maxlesadj
#endif
   implicit none
   type(prmtop_struct), intent(in) :: parmdata
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
   integer ntype,i
   integer i4,j,jj,k,l,n,nn
   integer iptres,nspsol,natsm
   _REAL_ dumd
   integer mol_atm_cnt
   _REAL_ massdiff   ! = mass[perturbed] - mass[original] for TI w.r.t. mass.

   integer :: res_start, res_end

   integer ier ! Allocation status.
   integer, intent(out) :: ierr
   ierr = 0

   ntype = ntypes*ntypes
   
#ifndef NO_ALLOCATABLES_IN_TYPE
   ! Copy the sections that _must_ be present
   
   x(l15:l15+natom-1) = parmdata%charge(1:natom) 
   x(lwinv:lwinv+natom-1) = parmdata%mass(1:natom)
   ix(i04:i04+natom-1) = parmdata%atom_type_index(1:natom)
   ix(i08:i08+natom-1) = parmdata%number_excluded_atoms(1:natom)
   ix(i06:i06+ntype-1) = parmdata%nonbonded_parm_index(1:ntype)
   ix(i02:i02+nres-1) = parmdata%residue_pointer(1:nres)
   ix(i02+nres) = natom + 1
   rk(1:numbnd) = parmdata%bond_force_constant(1:numbnd)
   req(1:numbnd) = parmdata%bond_equil_value(1:numbnd)
   tk(1:numang) = parmdata%angle_force_constant(1:numang)
   teq(1:numang) = parmdata%angle_equil_value(1:numang)
   pk(1:nptra) = parmdata%dihedral_force_constant(1:nptra)
   pn(1:nptra) = parmdata%dihedral_periodicity(1:nptra)
   phase(1:nptra) = parmdata%dihedral_phase(1:nptra)
   cn1(1:nttyp) = parmdata%lennard_jones_acoef(1:nttyp)
   cn2(1:nttyp) = parmdata%lennard_jones_bcoef(1:nttyp)
   ix(i10:i10+nnb-1) = parmdata%excluded_atoms_list(1:nnb)
   if (allocated(parmdata%hbond_acoef)) asol(1:nphb) = parmdata%hbond_acoef(1:nphb)
   if (allocated(parmdata%hbond_bcoef)) bsol(1:nphb) = parmdata%hbond_bcoef(1:nphb)
   if (allocated(parmdata%hbcut)) hbcut(1:nphb) = parmdata%hbcut(1:nphb)

   ! Copy the sections that need not be present

   if (allocated(parmdata%atomic_number)) then
      ix(i100+1:i100+natom) = parmdata%atomic_number(1:natom)
      ix(i100)=1
   else
      ix(i100)=0
   end if
   if (allocated(parmdata%atom_name)) then
      do i = 0, natom - 1
         i4 = i * 4 + 1
         write(ih(m04+i), '(4a1)') (parmdata%atom_name(i4+j), j=0, 3)
      end do
   else
      ih(m04:m04+natom-1) = '    '
   end if
   if (allocated(parmdata%residue_label)) then
      do i = 0, nres-1
         i4 = i * 4 + 1
         write(ih(m02+i), '(4a1)') (parmdata%residue_label(i4+j), j=0, 3)
      end do
   else
      ih(m02:m02+nres-1) = '    '
   end if
   if (allocated(parmdata%scee_scale_factor)) then
      one_scee(1:nptra) = 1.d0 / parmdata%scee_scale_factor(1:nptra)
   else
      one_scee(1:nptra) = 1.d0 / 1.2d0
   end if
   if (allocated(parmdata%scnb_scale_factor)) then
      one_scnb(1:nptra) = 1.d0 / parmdata%scnb_scale_factor(1:nptra)
   else
      one_scnb(1:nptra) = 1.d0 / 1.2d0
   end if
!  if (allocated(parmdata%solty)) then
!     solty(1:natyp) = parmdata%solty(1:natyp)
!  else
      solty(1:natyp) = 0.d0
!  end if
   if (lj1264 == 1) then
      cn6(1:nttyp) = parmdata%lennard_jones_ccoef(1:nttyp)
   else
      cn6(1:nttyp) = 0.d0
   end if
   if (vdwmodel == 1) then
      cn3(1:nttyp) = parmdata%expvdwmodel_beta(1:nttyp)
      cn4(1:nttyp) = parmdata%expvdwmodel_a(1:nttyp)
      cn5(1:nttyp) = parmdata%expvdwmodel_b(1:nttyp)
   end if

   ! Copy the parameter arrays

   do i = 1, nbonh
      j = i * 3 - 2
      ix(iibh+i-1)  = parmdata%bonds_inc_hydrogen(j  )
      ix(ijbh+i-1)  = parmdata%bonds_inc_hydrogen(j+1)
      ix(iicbh+i-1) = parmdata%bonds_inc_hydrogen(j+2)
   end do
   do i = 1, nbona
      j = i * 3 - 2
      ix(iiba+i-1)  = parmdata%bonds_without_hydrogen(j  )
      ix(ijba+i-1)  = parmdata%bonds_without_hydrogen(j+1)
      ix(iicba+i-1) = parmdata%bonds_without_hydrogen(j+2)
   end do
   do i = 1, ntheth
      j = i * 4 - 3
      ix(i24+i-1) = parmdata%angles_inc_hydrogen(j  )
      ix(i26+i-1) = parmdata%angles_inc_hydrogen(j+1)
      ix(i28+i-1) = parmdata%angles_inc_hydrogen(j+2)
      ix(i30+i-1) = parmdata%angles_inc_hydrogen(j+3)
   end do
   do i = 1, ntheta
      j = i * 4 - 3
      ix(i32+i-1) = parmdata%angles_without_hydrogen(j  )
      ix(i34+i-1) = parmdata%angles_without_hydrogen(j+1)
      ix(i36+i-1) = parmdata%angles_without_hydrogen(j+2)
      ix(i38+i-1) = parmdata%angles_without_hydrogen(j+3)
   end do
   do i = 1, nphih
      j = i * 5 - 4
      ix(i40+i-1) = parmdata%dihedrals_inc_hydrogen(j  )
      ix(i42+i-1) = parmdata%dihedrals_inc_hydrogen(j+1)
      ix(i44+i-1) = parmdata%dihedrals_inc_hydrogen(j+2)
      ix(i46+i-1) = parmdata%dihedrals_inc_hydrogen(j+3)
      ix(i48+i-1) = parmdata%dihedrals_inc_hydrogen(j+4)
   end do
   do i = 1, nphia
      j = i * 5 - 4
      ix(i50+i-1) = parmdata%dihedrals_without_hydrogen(j  )
      ix(i52+i-1) = parmdata%dihedrals_without_hydrogen(j+1)
      ix(i54+i-1) = parmdata%dihedrals_without_hydrogen(j+2)
      ix(i56+i-1) = parmdata%dihedrals_without_hydrogen(j+3)
      ix(i58+i-1) = parmdata%dihedrals_without_hydrogen(j+4)
   end do

   ! Copy other largely unecessary data if present

   if (allocated(parmdata%amber_atom_type)) then
      do i = 0, natom - 1
         i4 = i * 4 + 1
         write(ih(m06+i), '(4a1)') (parmdata%amber_atom_type(i4+j), j=0,3)
      end do
   else
      ih(m06:m06+natom-1) = '    '
   end if
!  if (allocated(parmdata%tree_chain_classification)) then
!     do i = 0, natom - 1
!        i4 = i * 4 + 1
!        write(ih(m08+i), '(4a1)') &
!              (parmdata%tree_chain_classification(i4+j), j=0,3)
!     end do
!  else
      ih(m08:m08+natom-1) = 'BLA '
!  end if
!  if (allocated(parmdata%join_array)) then
!     ix(i64:i64+natom-1) = parmdata%join_array(1:natom)
!  else
      ix(i64:i64+natom-1) = 0
!  end if
!  if (allocated(parmdata%irotat)) then
!     ix(i64:i64+natom-1) = parmdata%irotat(1:natom)
!  else
      ix(i64:i64+natom-1) = 0
!  end if
   
   ! PBC stuff
   
   nspm = 1
   ix(i70) = natom
   if (ifbox > 0) then
      iptres = parmdata%iptres
      nspm = parmdata%nspm
      nspsol = parmdata%nspsol
      ix(i70:i70+nspm-1) = parmdata%atoms_per_molecule(1:nspm)
      if( igb /= 0  .or. ipb /= 0 .or.  ntb == 0 )then
         box(1)=0.0d0
         box(2)=0.0d0
         box(3)=0.0d0
      else
         box(1)=a
         box(2)=b
         box(3)=c
      end if
   end if
   
   ! CAP stuff
   
   if(ifcap == 1) then
      natcap = parmdata%natcap
      cutcap = parmdata%cutcap
      xcap   = parmdata%xcap
      ycap   = parmdata%ycap
      zcap   = parmdata%zcap
   end if
   
   if ((igb /= 0 .or. ipb /= 0) .and. ifcap /= 0 .and. ifcap /= 5) then
      write(0,*) 'GB/PB calculations are incompatible with spherical solvent caps'
      ierr = 1
      return
   end if

   if (( (igb /= 0 .or. ipb /= 0) .and. (ifcap == 0 .or. ifcap == 5)) &
                                  .or.hybridgb>0.or.icnstph.gt.1) then
      x(l97:l97+natom-1) = parmdata%radii(1:natom)
      x(l96:l96+natom-1) = parmdata%screen(1:natom)
   end if
   
!
! IPOL is now specified in prmtop YD
!
   if (ipol > 0) then
      if (igb /= 0 .or. ipb /= 0) then
         write(0,*) 'GB/PB calculations are incompatible with polarizable force fields'
         ierr = 1
         return
      end if
      x(lpol:lpol+natom-1) = parmdata%polarizability(1:natom)
   end if

! Modified by WJM
   if (ipol > 1) then
      if (allocated(parmdata%dipole_damp_factor)) &
         x(ldf:ldf+natom-1) = parmdata%dipole_damp_factor(1:natom)
   end if

   call load_ewald_pol_info(ipol, x(lpol), x(lpol2), x(ldf), natom)

   ! Check that every atom is assigned to a molecule for NTP simulations. If
   ! not, segfaults or chaos may ensue.

   if (ntp .gt. 0) then

      mol_atm_cnt = 0

      do i = 1, nspm
         mol_atm_cnt = mol_atm_cnt + ix(i70+i-1)
      end do

      if (mol_atm_cnt .ne. natom) then
         write(6,'(a)') 'Error: Bad topology file. Sum of ATOMS_PER_MOLECULE &
                        &does not equal NATOM.'
         ierr = 1
         return
      end if

   end if

   !     ----- READ THE PERTURBED MASSES IF NEEDED  -----
   if (itimass > 0) then
      ! JVAN: dmdlm must be allocated here, not in pimd_init.
      allocate( dmdlm(1:natom),stat=ier )
      REQUIRE( ier == 0 )
      dmdlm(1:natom) = parmdata%ti_mass(1:natom)
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
      nlesty = parmdata%nlesty
      lestmp=nlesty*nlesty
      
      ! check the array sizes to make sure we do not overflow
      
      ! LES types
      
      if (nlesty > maxlestyp) then
         write (6,*) 'Exceeded MAXLESTYP',nlesty
         ierr = 1
         return
      end if
      
      ! LES atoms
      
      if (natom > maxles) then
         write (6,*) 'Exceeded MAXLES',natom
         ierr = 1
         return
      end if

      lestyp(1:natom) = parmdata%les_type(1:natom)
      lesfac(1:lestmp) = parmdata%les_fac(1:lestmp)
      cnum(1:natom) = parmdata%les_cnum(1:natom)
      subsp(1:natom) = parmdata%les_id(1:natom)
      
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
   
!  if (charmm_active) then
!    call read_charmm_params2(parmdata)
!  end if
!  if (cmap_active) then
!    call read_cmap_params2(parmdata)
!  end if
   
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
   
#endif /* NO_ALLOCATABLES_IN_TYPE */

   return
end subroutine cpparm2 
