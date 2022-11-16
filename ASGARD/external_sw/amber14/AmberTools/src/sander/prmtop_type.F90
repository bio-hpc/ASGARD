#define SAFE_READ read(LUN, fmt, err=666, end=666)
#define SAFE_DEALLOC(thing) if(allocated(parm%thing)) deallocate(parm%thing)

module prmtop_type

   type :: prmtop_struct

      sequence

      ! Integer pointers

      integer :: natom, ntypes, nbonh,  mbona, ntheth, mtheta, nphih
      integer :: mphia, nhparm, nparm,  nnb  , nres  , nbona , ntheta
      integer :: nphia, numbnd, numang, nptra, natyp , nphb  , ifpert
      integer :: nbper, ngper , ndper , mbper, mgper , mdper , ifbox
      integer :: nmxrs, ifcap
      integer :: numextra, ncopy

      integer :: iptres, nspm, nspsol ! solvent pointers

      integer :: natcap, nlesty

      ! Some bool flags

      integer :: is_chamber
      integer :: ipol

      ! Other single parameters

      double precision :: cutcap, xcap, ycap, zcap

      ! CHARMM-specific pointers

      integer :: cmap_term_count, cmap_type_count, charmm_nub, charmm_nubtypes
      integer :: charmm_nimphi, charmm_nimprtyp

      ! Data arrays

      character(len=80)       :: title
#ifndef NO_ALLOCATABLES_IN_TYPE
      character, allocatable, dimension(:) :: atom_name
      character, allocatable, dimension(:) :: residue_label
      character, allocatable, dimension(:) :: amber_atom_type
      character, allocatable, dimension(:) :: tree_chain_classification

      integer, allocatable, dimension(:)   :: atomic_number
      integer, allocatable, dimension(:)   :: atom_type_index
      integer, allocatable, dimension(:)   :: number_excluded_atoms
      integer, allocatable, dimension(:)   :: nonbonded_parm_index
      integer, allocatable, dimension(:)   :: residue_pointer
      integer, allocatable, dimension(:)   :: bonds_inc_hydrogen
      integer, allocatable, dimension(:)   :: bonds_without_hydrogen
      integer, allocatable, dimension(:)   :: angles_inc_hydrogen
      integer, allocatable, dimension(:)   :: angles_without_hydrogen
      integer, allocatable, dimension(:)   :: dihedrals_inc_hydrogen
      integer, allocatable, dimension(:)   :: dihedrals_without_hydrogen
      integer, allocatable, dimension(:)   :: excluded_atoms_list
      integer, allocatable, dimension(:)   :: join_array
      integer, allocatable, dimension(:)   :: irotat
      integer, allocatable, dimension(:)   :: atoms_per_molecule
      integer, allocatable, dimension(:)   :: les_type
      integer, allocatable, dimension(:)   :: les_fac
      integer, allocatable, dimension(:)   :: les_cnum
      integer, allocatable, dimension(:)   :: les_id
      integer, allocatable, dimension(:)   :: charmm_cmap_resolution
      integer, allocatable, dimension(:)   :: charmm_cmap_index
      integer, allocatable, dimension(:)   :: charmm_impropers
      integer, allocatable, dimension(:)   :: charmm_urey_bradley
      
      double precision, allocatable, dimension(:) :: charge
      double precision, allocatable, dimension(:) :: mass
      double precision, allocatable, dimension(:) :: bond_force_constant
      double precision, allocatable, dimension(:) :: bond_equil_value
      double precision, allocatable, dimension(:) :: angle_force_constant
      double precision, allocatable, dimension(:) :: angle_equil_value
      double precision, allocatable, dimension(:) :: dihedral_force_constant
      double precision, allocatable, dimension(:) :: dihedral_periodicity
      double precision, allocatable, dimension(:) :: dihedral_phase
      double precision, allocatable, dimension(:) :: scee_scale_factor
      double precision, allocatable, dimension(:) :: scnb_scale_factor
      double precision, allocatable, dimension(:) :: solty
      double precision, allocatable, dimension(:) :: lennard_jones_acoef
      double precision, allocatable, dimension(:) :: lennard_jones_bcoef
      double precision, allocatable, dimension(:) :: lennard_jones_ccoef
      double precision, allocatable, dimension(:) :: expvdwmodel_beta
      double precision, allocatable, dimension(:) :: expvdwmodel_a
      double precision, allocatable, dimension(:) :: expvdwmodel_b
      double precision, allocatable, dimension(:) :: hbond_acoef
      double precision, allocatable, dimension(:) :: hbond_bcoef
      double precision, allocatable, dimension(:) :: hbcut
      double precision, allocatable, dimension(:) :: box_dimensions
      double precision, allocatable, dimension(:) :: radii
      double precision, allocatable, dimension(:) :: screen
      double precision, allocatable, dimension(:) :: polarizability
      double precision, allocatable, dimension(:) :: dipole_damp_factor
      double precision, allocatable, dimension(:) :: ti_mass
      double precision, allocatable, dimension(:) :: charmm_cmap_parameter
      double precision, allocatable, dimension(:) :: charmm_improper_force_constant
      double precision, allocatable, dimension(:) :: charmm_improper_phase
      double precision, allocatable, dimension(:) :: lennard_jones_14_acoef
      double precision, allocatable, dimension(:) :: lennard_jones_14_bcoef
      double precision, allocatable, dimension(:) :: charmm_urey_bradley_equil_value
      double precision, allocatable, dimension(:) :: charmm_urey_bradley_force_constant

#endif /* NO_ALLOCATABLES_IN_TYPE */

   end type prmtop_struct

   contains

#ifndef NO_ALLOCATABLES_IN_TYPE
   ! This subroutine will parse a prmtop file and populate a prmtop struct with
   ! the data contained inside. The parm struct should be completely deallocated
   ! (or segfaults will likely occur). The flag ierr is set to 0 upon success or
   ! 1 upon failure
   subroutine read_prmtop_file(filename, parm, ierr)

      implicit none

      ! Passed parameters

      type(prmtop_struct), intent(in out) :: parm
      character(len=*), intent(in)        :: filename
      integer, intent(out)                :: ierr

      ! Local variables

      integer            :: i, j, i4
      integer            :: iok
      integer            :: nttyp, ntypes2
      integer            :: cmap_ptr
      integer, parameter :: LUN = 8
      character(len=80)  :: fmt, ifmt, afmt, rfmt
      character(len=80)  :: prmtop_flag, word
      character(len=4), allocatable :: charholder(:)
      double precision   :: dum1, dum2, dum3, dum4

      ierr = 0
      ! old-style prmtop file formats
      ifmt = '(12I6)'
      afmt = '(20A4)'
      rfmt = '(5E16.8)'

      ! Prevent nxtsec from printing out any information

      call nxtsec_silence
      call nxtsec_crd_silence

      ! Make sure all arrays are deallocated and all pointers zeroed

!     call destroy_prmtop_struct(parm)

      open(unit=LUN, file=filename, status='OLD', iostat=ierr)

      if (ierr /= 0) then
         ierr = 1
         return
      end if

      ! Start parsing the sections (in the same order as sander does)

      call nxtsec(LUN, 6, 1, '(A80)', 'CTITLE', fmt, iok)

      if (iok /= 0) then
         call nxtsec(LUN, 6, 1, '(A80)', 'TITLE', fmt, iok)
         if (iok /= 0) goto 666
         parm%is_chamber = 0
      else
         parm%is_chamber = 1
      end if
      SAFE_READ parm%title

      call nxtsec(LUN, 6, 1, ifmt, 'POINTERS', fmt, iok)
      if (iok == -2) goto 666
      SAFE_READ &
         parm%natom,  parm%ntypes, parm%nbonh, parm%mbona,  parm%ntheth, &
         parm%mtheta, parm%nphih,  parm%mphia, parm%nhparm, parm%nparm,  &
         parm%nnb,    parm%nres,   parm%nbona, parm%ntheta, parm%nphia,  &
         parm%numbnd, parm%numang, parm%nptra, parm%natyp,  parm%nphb,   &
         parm%ifpert, parm%nbper,  parm%ngper, parm%ndper,  parm%mbper,  &
         parm%mgper,  parm%mdper,  parm%ifbox, parm%nmxrs,  parm%ifcap,  &
         parm%numextra, parm%ncopy

      nttyp = parm%ntypes * (parm%ntypes + 1) / 2 ! number of LJ type pairs
      ntypes2 = parm%ntypes * parm%ntypes

      call nxtsec(LUN, 6, 1, ifmt, 'IPOL', fmt, iok)
      if (iok == 0) SAFE_READ parm%ipol

      allocate(charholder(parm%natom))
      call nxtsec(LUN, 6, 1, afmt, 'ATOM_NAME', fmt, iok)
      if (iok == -2) goto 666
      SAFE_READ (charholder(i), i=1, parm%natom)
      allocate(parm%atom_name(4*parm%natom))
      do i = 0, parm%natom - 1
         i4 = i * 4 + 1
         read(charholder(i+1), '(4a1)') (parm%atom_name(i4+j), j=0, 3)
      end do
      deallocate(charholder)

      call nxtsec(LUN, 6, 1, rfmt, 'CHARGE', fmt, iok)
      if (iok == -2) goto 666
      allocate(parm%charge(parm%natom))
      SAFE_READ (parm%charge(i), i=1, parm%natom)

      call nxtsec(LUN, 6, 1, ifmt, 'ATOMIC_NUMBER', fmt, iok)
      if (iok == 0) then
         allocate(parm%atomic_number(parm%natom))
         SAFE_READ (parm%atomic_number(i),i=1,parm%natom)
      end if

      call nxtsec(LUN, 6, 1, rfmt, 'MASS', fmt, iok)
      if (iok == -2) goto 666
      allocate(parm%mass(parm%natom))
      SAFE_READ (parm%mass(i), i=1, parm%natom)

      call nxtsec(LUN, 6, 1, ifmt, 'ATOM_TYPE_INDEX', fmt, iok)
      if (iok == -2) goto 666
      allocate(parm%atom_type_index(parm%natom))
      SAFE_READ (parm%atom_type_index(i),i=1,parm%natom)

      call nxtsec(LUN, 6, 1, ifmt, 'NUMBER_EXCLUDED_ATOMS', fmt, iok)
      if (iok == -2) goto 666
      allocate(parm%number_excluded_atoms(parm%natom))
      SAFE_READ (parm%number_excluded_atoms(i), i=1, parm%natom)

      call nxtsec(LUN, 6, 1, ifmt, 'NONBONDED_PARM_INDEX', fmt, iok)
      if (iok == -2) goto 666
      allocate(parm%nonbonded_parm_index(ntypes2))
      SAFE_READ (parm%nonbonded_parm_index(i), i=1, ntypes2)

      call nxtsec(LUN, 6, 1, afmt, 'RESIDUE_LABEL', fmt, iok)
      if (iok == -2) goto 666
      allocate(charholder(parm%nres))
      SAFE_READ (charholder(i), i=1, parm%nres)
      allocate(parm%residue_label(4*parm%nres))
      do i = 0, parm%nres - 1
         i4 = i*4 + 1
         read(charholder(i+1), '(4a1)') (parm%residue_label(i4+j), j=0, 3)
      end do
      deallocate(charholder)

      call nxtsec(LUN, 6, 1, ifmt, 'RESIDUE_POINTER', fmt, iok)
      if (iok == -2) goto 666
      allocate(parm%residue_pointer(parm%nres))
      SAFE_READ (parm%residue_pointer(i), i=1, parm%nres)

      call nxtsec(LUN, 6, 1, rfmt, 'BOND_FORCE_CONSTANT', fmt, iok)
      if (iok == -2) goto 666
      allocate(parm%bond_force_constant(parm%numbnd))
      SAFE_READ (parm%bond_force_constant(i), i=1, parm%numbnd)

      call nxtsec(LUN, 6, 1, rfmt, 'BOND_EQUIL_VALUE', fmt, iok)
      if (iok == -2) goto 666
      allocate(parm%bond_equil_value(parm%numbnd))
      SAFE_READ (parm%bond_equil_value(i), i=1, parm%numbnd)

      call nxtsec(LUN, 6, 1, rfmt, 'ANGLE_FORCE_CONSTANT', fmt, iok)
      if (iok == -2) goto 666
      allocate(parm%angle_force_constant(parm%numang))
      SAFE_READ (parm%angle_force_constant(i), i=1, parm%numang)

      call nxtsec(LUN, 6, 1, rfmt, 'ANGLE_EQUIL_VALUE', fmt, iok)
      if (iok == -2) goto 666
      allocate(parm%angle_equil_value(parm%numang))
      SAFE_READ (parm%angle_equil_value(i), i=1, parm%numang)

      call nxtsec(LUN, 6, 1, rfmt, 'DIHEDRAL_FORCE_CONSTANT', fmt, iok)
      if (iok == -2) goto 666
      allocate(parm%dihedral_force_constant(parm%nptra))
      SAFE_READ (parm%dihedral_force_constant(i), i=1, parm%nptra)

      call nxtsec(LUN, 6, 1, rfmt, 'DIHEDRAL_PERIODICITY', fmt, iok)
      if (iok == -2) goto 666
      allocate(parm%dihedral_periodicity(parm%nptra))
      SAFE_READ (parm%dihedral_periodicity(i), i=1, parm%nptra)

      call nxtsec(LUN, 6, 1, rfmt, 'DIHEDRAL_PHASE', fmt, iok)
      if (iok == -2) goto 666
      allocate(parm%dihedral_phase(parm%nptra))
      SAFE_READ (parm%dihedral_phase(i), i=1, parm%nptra)

      call nxtsec(LUN, 6, 1, rfmt, 'SCEE_SCALE_FACTOR', fmt, iok)
      allocate(parm%scee_scale_factor(parm%nptra))
      if (iok == 0) then
         SAFE_READ (parm%scee_scale_factor(i), i=1, parm%nptra)
      else
         parm%scee_scale_factor(:) = 1.2d0
      end if

      call nxtsec(LUN, 6, 1, rfmt, 'SCNB_SCALE_FACTOR', fmt, iok)
      allocate(parm%scnb_scale_factor(parm%nptra))
      if (iok == 0) then
         SAFE_READ (parm%scnb_scale_factor(i), i=1, parm%nptra)
      else
         parm%scnb_scale_factor(:) = 2.d0
      end if

      call nxtsec(LUN, 6, 1, rfmt, 'SOLTY', fmt, iok)
      if (iok == -2) goto 666
      allocate(parm%solty(parm%natyp))
      SAFE_READ (parm%solty(i), i=1, parm%natyp)

      call nxtsec(LUN, 6, 1, rfmt, 'LENNARD_JONES_ACOEF', fmt, iok)
      if (iok == -2) goto 666
      allocate(parm%lennard_jones_acoef(nttyp))
      SAFE_READ (parm%lennard_jones_acoef(i), i=1, nttyp)

      call nxtsec(LUN, 6, 1, rfmt, 'LENNARD_JONES_BCOEF', fmt, iok)
      if (iok == -2) goto 666
      allocate(parm%lennard_jones_bcoef(nttyp))
      SAFE_READ (parm%lennard_jones_bcoef(i), i=1, nttyp)

      call nxtsec(LUN, 6, 1, rfmt, 'LENNARD_JONES_CCOEF', fmt, iok)
      if (iok == 0) then
         allocate(parm%lennard_jones_ccoef(nttyp))
         SAFE_READ (parm%lennard_jones_ccoef(i), i=1, nttyp)
      end if

      call nxtsec(LUN, 6, 1, rfmt, 'EXPVDWMODEL_BETA', fmt, iok)
      if (iok == 0) then
         allocate(parm%expvdwmodel_beta(nttyp))
         SAFE_READ (parm%expvdwmodel_beta(i), i=1, nttyp)
      end if

      call nxtsec(LUN, 6, 1, rfmt, 'EXPVDWMODEL_A', fmt, iok)
      if (iok == 0) then
         allocate(parm%expvdwmodel_a(nttyp))
         SAFE_READ (parm%expvdwmodel_a(i), i=1, nttyp)
      end if

      call nxtsec(LUN, 6, 1, rfmt, 'EXPVDWMODEL_B', fmt, iok)
      if (iok == 0) then
         allocate(parm%expvdwmodel_b(nttyp))
         SAFE_READ (parm%expvdwmodel_b(i), i=1, nttyp)
      end if

      call nxtsec(LUN, 6, 1, ifmt, 'BONDS_INC_HYDROGEN', fmt, iok)
      if (iok == -2) goto 666
      allocate(parm%bonds_inc_hydrogen(parm%nbonh*3))
      SAFE_READ (parm%bonds_inc_hydrogen(i), i=1, parm%nbonh*3)

      call nxtsec(LUN, 6, 1, ifmt, 'BONDS_WITHOUT_HYDROGEN', fmt, iok)
      if (iok == -2) goto 666
      allocate(parm%bonds_without_hydrogen(parm%nbona*3))
      SAFE_READ (parm%bonds_without_hydrogen(i), i=1, parm%nbona*3)

      call nxtsec(LUN, 6, 1, ifmt, 'ANGLES_INC_HYDROGEN', fmt, iok)
      if (iok == -2) goto 666
      allocate(parm%angles_inc_hydrogen(parm%ntheth*4))
      SAFE_READ (parm%angles_inc_hydrogen(i), i=1, parm%ntheth*4)

      call nxtsec(LUN, 6, 1, ifmt, 'ANGLES_WITHOUT_HYDROGEN', fmt, iok)
      if (iok == -2) goto 666
      allocate(parm%angles_without_hydrogen(parm%ntheta*4))
      SAFE_READ (parm%angles_without_hydrogen(i), i=1, parm%ntheta*4)

      call nxtsec(LUN, 6, 1, ifmt, 'DIHEDRALS_INC_HYDROGEN', fmt, iok)
      if (iok == -2) goto 666
      allocate(parm%dihedrals_inc_hydrogen(parm%nphih*5))
      SAFE_READ (parm%dihedrals_inc_hydrogen(i), i=1, parm%nphih*5)

      call nxtsec(LUN, 6, 1, ifmt, 'DIHEDRALS_WITHOUT_HYDROGEN', fmt, iok)
      if (iok == -2) goto 666
      allocate(parm%dihedrals_without_hydrogen(parm%nphia*5))
      SAFE_READ (parm%dihedrals_without_hydrogen(i), i=1, parm%nphia*5)

      call nxtsec(LUN, 6, 1, ifmt, 'EXCLUDED_ATOMS_LIST', fmt, iok)
      if (iok == -2) goto 666
      allocate(parm%excluded_atoms_list(parm%nnb))
      SAFE_READ (parm%excluded_atoms_list(i), i=1, parm%nnb)

      call nxtsec(LUN, 6, 1, rfmt, 'HBOND_ACOEF', fmt, iok)
      if (iok == -2) goto 666
      allocate(parm%hbond_acoef(parm%nphb))
      SAFE_READ (parm%hbond_acoef(i), i=1, parm%nphb)

      call nxtsec(LUN, 6, 1, rfmt, 'HBOND_BCOEF', fmt, iok)
      if (iok == -2) goto 666
      allocate(parm%hbond_bcoef(parm%nphb))
      SAFE_READ (parm%hbond_bcoef(i), i=1, parm%nphb)

      call nxtsec(LUN, 6, 1, rfmt, 'HBCUT', fmt, iok)
      if (iok == -2) goto 666
      allocate(parm%hbcut(parm%nphb))
      SAFE_READ (parm%hbcut(i), i=1, parm%nphb)

      call nxtsec(LUN, 6, 1, afmt, 'AMBER_ATOM_TYPE', fmt, iok)
      if (iok == -2) goto 666
      allocate(charholder(parm%natom))
      SAFE_READ (charholder(i), i=1, parm%natom)
      allocate(parm%amber_atom_type(parm%natom*4))
      do i = 0, parm%natom - 1
         i4 = i * 4 + 1
         read(charholder(i+1), '(4a1)') (parm%amber_atom_type(i4+j), j=0, 3)
      end do

      call nxtsec(LUN, 6, 1, afmt, 'TREE_CHAIN_CLASSIFICATION', fmt, iok)
      if (iok == -2) goto 666
      SAFE_READ (charholder(i), i=1, parm%natom)
      allocate(parm%tree_chain_classification(parm%natom*4))
      do i = 0, parm%natom - 1
         i4 = i * 4 + 1
         read(charholder(i+1), '(4a1)') &
               (parm%tree_chain_classification(i4+j), j=0, 3)
      end do
      deallocate(charholder)

      call nxtsec(LUN, 6, 1, ifmt, 'JOIN_ARRAY', fmt, iok)
      if (iok == -2) goto 666
      allocate(parm%join_array(parm%natom))
      SAFE_READ (parm%join_array(i), i=1, parm%natom)

      call nxtsec(LUN, 6, 1, ifmt, 'IROTAT', fmt, iok)
      if (iok == -2) goto 666
      allocate(parm%irotat(parm%natom))
      SAFE_READ (parm%irotat(i), i=1, parm%natom)

      if (parm%ifbox > 0) then
         call nxtsec(LUN, 6, 1, ifmt, 'SOLVENT_POINTERS', fmt, iok)
         if (iok == -2) goto 666
         SAFE_READ parm%iptres, parm%nspm, parm%nspsol

         call nxtsec(LUN, 6, 1, ifmt, 'ATOMS_PER_MOLECULE', fmt, iok)
         if (iok == -2) goto 666
         allocate(parm%atoms_per_molecule(parm%nspm))
         SAFE_READ (parm%atoms_per_molecule(i), i=1, parm%nspm)

         call nxtsec(LUN, 6, 1, rfmt, 'BOX_DIMENSIONS', fmt, iok)
         if (iok == -1) &
            SAFE_READ dum1, dum2, dum3, dum4
      end if

      if (parm%ifcap == 1) then
         call nxtsec(LUN, 6, 1, '(I6)', 'CAP_INFO', fmt, iok)
         if (iok == -2) goto 666
         SAFE_READ parm%natcap

         call nxtsec(LUN, 6, 1, '(4E16.8)', 'CAP_INFO2', fmt, iok)
         if (iok == -2) goto 666
         SAFE_READ parm%cutcap, parm%xcap, parm%ycap, parm%zcap
      end if

      call nxtsec(LUN, 6, 1, rfmt, 'RADII', fmt, iok)
      if (iok == 0) then
         allocate(parm%radii(parm%natom))
         SAFE_READ (parm%radii(i), i=1, parm%natom)
      end if

      call nxtsec(LUN, 6, 1, rfmt, 'SCREEN', fmt, iok)
      if (iok == 0) then
         allocate(parm%screen(parm%natom))
         SAFE_READ (parm%screen(i), i=1, parm%natom)
      end if

      if (parm%ipol > 0) then
         call nxtsec(LUN, 6, 1, rfmt, 'POLARIZABILITY', fmt, iok)
         if (iok == 0) then
            allocate(parm%polarizability(parm%natom))
            SAFE_READ (parm%polarizability(i), i=1, parm%natom)
         end if

         if (parm%ipol > 1) then
            call nxtsec(LUN, 6, 1, rfmt, 'DIPOLE_DAMP_FACTOR', fmt, iok)
            if (iok == 0) then
               allocate(parm%dipole_damp_factor(parm%natom))
               SAFE_READ (parm%dipole_damp_factor(i), i=1, parm%natom)
            end if
         end if
      end if

      call nxtsec(LUN, 6, 1, rfmt, 'TI_MASS', fmt, iok)
      if (iok == 0) then
         allocate(parm%ti_mass(parm%natom))
         SAFE_READ (parm%ti_mass(i), i=1, parm%natom)
      end if

      call nxtsec(LUN, 6, 1, ifmt, 'LES_NTYP', fmt, iok)
      if (iok == 0) SAFE_READ parm%nlesty

      call nxtsec(LUN, 6, 1, ifmt, 'LES_TYPE', fmt, iok)
      if (iok == 0) then
         allocate(parm%les_type(parm%natom))
         SAFE_READ (parm%les_type(i), i=1, parm%natom)
      end if

      call nxtsec(LUN, 6, 1, rfmt, 'LES_FAC', fmt, iok)
      if (iok == 0) then
         allocate(parm%les_fac(parm%nlesty*parm%nlesty))
         SAFE_READ (parm%les_fac(i), i=1, parm%nlesty*parm%nlesty)
      end if

      call nxtsec(LUN, 6, 1, ifmt, 'LES_CNUM', fmt, iok)
      if (iok == 0) then
         allocate(parm%les_cnum(parm%natom))
         SAFE_READ (parm%les_cnum(i), i=1, parm%natom)
      end if

      call nxtsec(LUN, 6, 1, ifmt, 'LES_ID', fmt, iok)
      if (iok == 0) then
         allocate(parm%les_id(parm%natom))
         SAFE_READ (parm%les_id(i), i=1, parm%natom)
      end if

      if (parm%is_chamber > 0) then
         call nxtsec(LUN, 6, 1, ifmt, 'CHARMM_UREY_BRADLEY_COUNT', fmt, iok)
         if (iok /= 0) goto 666
         SAFE_READ parm%charmm_nub, parm%charmm_nubtypes

         call nxtsec(LUN, 6, 1, ifmt, 'CHARMM_NUM_IMPROPERS', fmt, iok)
         if (iok /= 0) goto 666
         SAFE_READ parm%charmm_nimphi

         call nxtsec(LUN, 6, 1, ifmt, 'CHARMM_NUM_IMPR_TYPES', fmt, iok)
         if (iok /= 0) goto 666
         SAFE_READ parm%charmm_nimprtyp

         prmtop_flag = 'CHARMM_IMPROPER_FORCE_CONSTANT'
         call nxtsec(LUN, 6, 1, rfmt, prmtop_flag, fmt, iok)
         if (iok /= 0) goto 666
         allocate(parm%charmm_improper_force_constant(parm%charmm_nimprtyp))
         SAFE_READ (parm%charmm_improper_force_constant(i), &
                    i=1, parm%charmm_nimprtyp)

         prmtop_flag = 'CHARMM_IMPROPER_PHASE'
         call nxtsec(LUN, 6, 1, rfmt, prmtop_flag, fmt, iok)
         if (iok /= 0) goto 666
         allocate(parm%charmm_improper_phase(parm%charmm_nimprtyp))
         SAFE_READ (parm%charmm_improper_phase(i), i=1, parm%charmm_nimprtyp)

         prmtop_flag = 'CHARMM_IMPROPERS'
         call nxtsec(LUN, 6, 1, ifmt, prmtop_flag, fmt, iok)
         if (iok /= 0) goto 666
         allocate(parm%charmm_impropers(parm%charmm_nimphi*5))
         SAFE_READ (parm%charmm_impropers(i), i=1, parm%charmm_nimphi*5)

         prmtop_flag = 'LENNARD_JONES_14_ACOEF'
         call nxtsec(LUN, 6, 1, rfmt, prmtop_flag, fmt, iok)
         if (iok /= 0) goto 666
         allocate(parm%lennard_jones_14_acoef(nttyp))
         SAFE_READ (parm%lennard_jones_14_acoef(i), i=1, nttyp)

         prmtop_flag = 'LENNARD_JONES_14_BCOEF'
         call nxtsec(LUN, 6, 1, rfmt, prmtop_flag, fmt, iok)
         if (iok /= 0) goto 666
         allocate(parm%lennard_jones_14_bcoef(nttyp))
         SAFE_READ (parm%lennard_jones_14_bcoef(i), i=1, nttyp)

         prmtop_flag = 'CHARMM_UREY_BRADLEY'
         call nxtsec(LUN, 6, 1, ifmt, prmtop_flag, fmt, iok)
         if (iok /= 0) goto 666
         allocate(parm%charmm_urey_bradley(parm%charmm_nub*3))
         SAFE_READ (parm%charmm_urey_bradley(i), i=1, parm%charmm_nub*3)

         prmtop_flag = 'CHARMM_UREY_BRADLEY_EQUIL_VALUE'
         call nxtsec(LUN, 6, 1, rfmt, prmtop_flag, fmt, iok)
         if (iok /= 0) goto 666
         allocate(parm%charmm_urey_bradley_equil_value(parm%charmm_nubtypes))
         SAFE_READ (parm%charmm_urey_bradley_equil_value(i), &
                    i=1, parm%charmm_nubtypes)

         prmtop_flag = 'CHARMM_UREY_BRADLEY_FORCE_CONSTANT'
         call nxtsec(LUN, 6, 1, rfmt, prmtop_flag, fmt, iok)
         if (iok /= 0) goto 666
         allocate(parm%charmm_urey_bradley_force_constant(parm%charmm_nubtypes))
         SAFE_READ (parm%charmm_urey_bradley_force_constant(i), &
                    i=1, parm%charmm_nubtypes)

         prmtop_flag = 'CHARMM_CMAP_COUNT'
         call nxtsec(LUN, 6, 1, ifmt, prmtop_flag, fmt, iok)
         if (iok == 0) then
            SAFE_READ parm%cmap_term_count, parm%cmap_type_count

            prmtop_flag = 'CHARMM_CMAP_RESOLUTION'
            call nxtsec(LUN, 6, 1, ifmt, prmtop_flag, fmt, iok)
            allocate(parm%charmm_cmap_resolution(parm%cmap_type_count))
            SAFE_READ (parm%charmm_cmap_resolution(i), &
                       i=1, parm%cmap_type_count)

            ! Figure out how many numbers are in each CMAP based on their
            ! resolution and allocate enough space for all cmaps
            cmap_ptr = 0
            do j = 1, parm%cmap_type_count
               i4 = parm%charmm_cmap_resolution(j)
               cmap_ptr = cmap_ptr + i4*i4
            end do
            allocate(parm%charmm_cmap_parameter(cmap_ptr))
            cmap_ptr = 0

            ! Now parse all of the sections
            do j = 1, parm%cmap_type_count
               write(word, '(i2.2)') j
               prmtop_flag = 'CHARMM_CMAP_PARAMETER_' // trim(word)
               call nxtsec(LUN, 6, 1, rfmt, prmtop_flag, fmt, iok)
               if (iok /= 0) goto 666
               i4 = parm%charmm_cmap_resolution(j)
               i4 = i4 * i4
               SAFE_READ (parm%charmm_cmap_parameter(cmap_ptr+i), i=1, i4)
               cmap_ptr = cmap_ptr + i4
            end do

            prmtop_flag = 'CHARMM_CMAP_INDEX'
            call nxtsec(LUN, 6, 1, ifmt, prmtop_flag, fmt, iok)
            if (iok /= 0) goto 666
            allocate(parm%charmm_cmap_index(parm%cmap_term_count*6))
            SAFE_READ (parm%charmm_cmap_index(i), i=1, parm%cmap_term_count*6)
         end if
      end if

      ! Now make sure we're not doing AMOEBA

      call nxtsec(LUN, 6, 1, ifmt, 'AMOEBA_FORCEFIELD', fmt, iok)
      if (iok == 0) goto 666

      return

666   continue
      ! Error condition
      close(LUN)
      if (allocated(charholder)) deallocate(charholder)
      call destroy_prmtop_struct(parm)
      ierr = -1
      return

   end subroutine read_prmtop_file


   ! Destroys a parm instance by deallocating all of its memory and resetting
   ! its flags
   subroutine destroy_prmtop_struct(parm)

      implicit none

      type(prmtop_struct), intent(in out) :: parm

      parm%natom = 0
      parm%ntypes = 0
      parm%nbonh = 0
      parm%mbona = 0
      parm%ntheth = 0
      parm%mtheta = 0
      parm%nphih = 0
      parm%mphia = 0
      parm%nhparm = 0
      parm%nparm = 0
      parm%nnb = 0
      parm%nres = 0
      parm%nbona = 0
      parm%ntheta = 0
      parm%nphia = 0
      parm%numbnd = 0
      parm%numang = 0
      parm%nptra = 0
      parm%natyp = 0
      parm%nphb = 0
      parm%ifpert = 0
      parm%nbper = 0
      parm%ngper = 0
      parm%ndper = 0
      parm%mbper = 0
      parm%mgper = 0
      parm%mdper = 0
      parm%ifbox = 0
      parm%nmxrs = 0
      parm%ifcap = 0
      parm%numextra = 0
      parm%ncopy = 0

      parm%iptres = 0
      parm%nspm = 0
      parm%nspsol = 0
      parm%natcap = 0
      parm%nlesty = 0

      parm%is_chamber = 0
      parm%ipol = 0

      parm%cutcap = 0.d0
      parm%xcap = 0.d0
      parm%ycap = 0.d0
      parm%zcap = 0.d0

      parm%cmap_term_count = 0
      parm%cmap_type_count = 0
      parm%charmm_nub = 0
      parm%charmm_nubtypes = 0
      parm%charmm_nimphi = 0
      parm%charmm_nimprtyp = 0

      parm%title = ' '

      SAFE_DEALLOC(atom_name)
      SAFE_DEALLOC(residue_label)
      SAFE_DEALLOC(amber_atom_type)
      SAFE_DEALLOC(tree_chain_classification)

      SAFE_DEALLOC(atomic_number)
      SAFE_DEALLOC(atom_type_index)
      SAFE_DEALLOC(number_excluded_atoms)
      SAFE_DEALLOC(nonbonded_parm_index)
      SAFE_DEALLOC(residue_pointer)
      SAFE_DEALLOC(bonds_inc_hydrogen)
      SAFE_DEALLOC(bonds_without_hydrogen)
      SAFE_DEALLOC(angles_inc_hydrogen)
      SAFE_DEALLOC(angles_without_hydrogen)
      SAFE_DEALLOC(dihedrals_inc_hydrogen)
      SAFE_DEALLOC(dihedrals_without_hydrogen)
      SAFE_DEALLOC(excluded_atoms_list)
      SAFE_DEALLOC(join_array)
      SAFE_DEALLOC(irotat)
      SAFE_DEALLOC(atoms_per_molecule)
      SAFE_DEALLOC(les_type)
      SAFE_DEALLOC(les_fac)
      SAFE_DEALLOC(les_cnum)
      SAFE_DEALLOC(les_id)
      SAFE_DEALLOC(charmm_cmap_resolution)
      SAFE_DEALLOC(charmm_cmap_index)
      SAFE_DEALLOC(charmm_impropers)
      SAFE_DEALLOC(charmm_urey_bradley)
      
      SAFE_DEALLOC(charge)
      SAFE_DEALLOC(mass)
      SAFE_DEALLOC(bond_force_constant)
      SAFE_DEALLOC(bond_equil_value)
      SAFE_DEALLOC(angle_force_constant)
      SAFE_DEALLOC(angle_equil_value)
      SAFE_DEALLOC(dihedral_force_constant)
      SAFE_DEALLOC(dihedral_periodicity)
      SAFE_DEALLOC(dihedral_phase)
      SAFE_DEALLOC(scee_scale_factor)
      SAFE_DEALLOC(scnb_scale_factor)
      SAFE_DEALLOC(solty)
      SAFE_DEALLOC(lennard_jones_acoef)
      SAFE_DEALLOC(lennard_jones_bcoef)
      SAFE_DEALLOC(lennard_jones_ccoef)
      SAFE_DEALLOC(expvdwmodel_beta)
      SAFE_DEALLOC(expvdwmodel_a)
      SAFE_DEALLOC(expvdwmodel_b)
      SAFE_DEALLOC(hbond_acoef)
      SAFE_DEALLOC(hbond_bcoef)
      SAFE_DEALLOC(hbcut)
      SAFE_DEALLOC(box_dimensions)
      SAFE_DEALLOC(radii)
      SAFE_DEALLOC(screen)
      SAFE_DEALLOC(polarizability)
      SAFE_DEALLOC(dipole_damp_factor)
      SAFE_DEALLOC(ti_mass)
      SAFE_DEALLOC(charmm_cmap_parameter)
      SAFE_DEALLOC(charmm_improper_force_constant)
      SAFE_DEALLOC(charmm_improper_phase)
      SAFE_DEALLOC(lennard_jones_14_acoef)
      SAFE_DEALLOC(lennard_jones_14_bcoef)
      SAFE_DEALLOC(charmm_urey_bradley_equil_value)
      SAFE_DEALLOC(charmm_urey_bradley_force_constant)

   end subroutine destroy_prmtop_struct

#else

   subroutine read_prmtop_file(filename, parm, ierr)

      implicit none

      ! Passed parameters

      type(prmtop_struct), intent(in out) :: parm
      character(len=*), intent(in)        :: filename
      integer, intent(out)                :: ierr

      write(0,*) 'Compiler does not support allocatables in types'
      stop 1
   end subroutine read_prmtop_file

   ! Destroys a parm instance by deallocating all of its memory and resetting
   ! its flags
   subroutine destroy_prmtop_struct(parm)

      implicit none

      type(prmtop_struct), intent(in out) :: parm

   end subroutine destroy_prmtop_struct

#endif /* NO_ALLOCATABLES_IN_TYPE */

end module prmtop_type

! remove the helper macros
#undef SAFE_READ
#undef SAFE_DEALLOC
