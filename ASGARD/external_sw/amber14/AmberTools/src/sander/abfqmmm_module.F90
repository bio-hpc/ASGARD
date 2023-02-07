! <compile=optimized>
#include "copyright.h"
#include "../include/assert.fh"
#include "../include/dprec.fh"
module abfqmmm_module
!----------------------------------------------------------------------
! Adaptive buffered-force QM/MM
! Author: Letif Mones
!         <lam81@cam.ac.uk>
! Date: December 2011
!
!----------------------------------------------------------------------

  use charmm_mod, only: charmm_active
  use charmm_mod, only: charmm_nub, chm_ang_ub_struct, charmm_ang_ub
  use charmm_mod, only: charmm_nimphi, chm_imp_struct, charmm_imp
  use charmm_mod, only: cmap_term_count, cmap_index

  implicit none

  private

  public :: abfqmmm_parameters
  public :: abfqmmm_param
  public :: abfqmmm_init_param
  public :: abfqmmm_connection_list
  public :: abfqmmm_cut_bond_list
  public :: abfqmmm_oxidation_number_list
  public :: abfqmmm_setup
  public :: abfqmmm_residue_labelling
  public :: abfqmmm_recursive_labelling
  public :: abfqmmm_recursive_find_label_heavy_atom
  public :: abfqmmm_update_qmatoms
  public :: abfqmmm_atom_group_min_distance
  public :: abfqmmm_residue_group_min_distance
  public :: abfqmmm_atom_labelling_by_cutoff
  public :: abfqmmm_res_labelling_by_cutoff
  public :: abfqmmm_atom_atom_dist
  public :: abfqmmm_point_point_dist
  public :: abfqmmm_point_atom_dist
  public :: abfqmmm_point_res_dist
  public :: abfqmmm_center_and_radius_of_atom_list
  public :: abfqmmm_select_system_qmatoms
  public :: abfqmmm_allocate_arrays_of_parameters
  public :: abfqmmm_store_parameters
  public :: abfqmmm_set_parameters
  public :: abfqmmm_vel_verlet1
  public :: abfqmmm_combine_forces
  public :: abfqmmm_vel_verlet2
  public :: abfqmmm_read_idrst
  public :: abfqmmm_write_idrst
  public :: abfqmmm_write_pdb
! public :: abfqmmm_write_outputs
  public :: abfqmmm_next_step
  public :: abfqmmm_diffusion_restraint
! public :: abfqmmm_diffusion_potential
  public :: abfqmmm_diffusion_force
  public :: abfqmmm_calc_diff_coeff

  type abfqmmm_parameters

    integer :: abfqmmm                                       ! is calculation ABF QM/MM? 0 - no (default), 1 - yes
    integer :: hot_spot                                      ! is calculation hot spot? 0 - no (default), 1 - yes
    integer :: natom                                         ! number of atoms
    integer :: nres                                          ! number of residues
    integer, dimension(:), pointer :: id => null()           ! id of each atom: 1 - user defined (and automatically corrected) core atoms
                                                             !                  2 - core atoms around user specified core atoms (r_core_in/r_core_out)
                                                             !                  3 - user defined (and automatically corrected) qm atoms
                                                             !                  4 - qm atoms around user specified qm atoms (r_qm_in/r_qm_out)
                                                             !                  5 - user defined (and automatically corrected) buffer atoms 
                                                             !                  6 - buffer atoms around qm atoms (r_buffer_in/r_buffer_out)
                                                             !                  7 - mm atoms
    integer, dimension(:), pointer :: id_orig => null()      ! original id of each atom stored
    integer, dimension(:), pointer :: isqm => null()         ! is given atom QM in a given simulation? 0 - no, 1 - yes
    integer :: n_user_core                                   ! number of user defined core atoms (with possible correction)
    integer, dimension(:), pointer :: user_core => null()    ! array of user defined core atoms (with possible correction)     
    integer :: n_user_qm                                     ! number of user defined qm + user defind core atoms (with possible correction)
    integer, dimension(:), pointer :: user_qm => null()      ! array of user defined qm + user defined core atoms (with possible correction)
    integer :: n_user_buffer                                 ! number of user defined buffer atoms (with possible correction)
    integer, dimension(:), pointer :: user_buffer => null()  ! array of user defined buffer atoms (with possible correction)
    integer :: n_core                                        ! number of user defined + extended core atoms
    integer, dimension(:), pointer :: core => null()         ! array of user defined + extended core atoms
    integer :: n_qm                                          ! number of user defined + extended core + user defined + extended qm atoms
    integer, dimension(:), pointer :: qm => null()           ! array of user defined + extended core + user defined + extended qm atoms
    integer :: n_buffer                                      ! number of user defined + extended buffer atoms
    integer, dimension(:), pointer :: buffer => null()       ! array of user defined + extended buffer atoms

    integer :: n_subset_core                                   ! number of core subset atoms / subset residues (Fixed atom list): natoms / nres (fixed atom list) if not specified by the user
    integer, dimension(:), pointer :: subset_core => null()    ! array of user defined core subset atoms / subset residues (fixed atom list)
    integer :: n_subset_qm                                     ! number of qm subset atoms / subset residues (fixed atom list): natoms / nres (fixed atom list) if not specified by the user
    integer, dimension(:), pointer :: subset_qm => null()      ! array of user defined qm subset atoms / subset residues (fixed atom list)
    integer :: n_subset_buffer                                 ! number of buffer subset atoms / subset residues (fixed atom list): natoms /nres (fixed atom list) if not specified by the user
    integer, dimension(:), pointer :: subset_buffer => null()  ! array of user defined buffer subset atoms / subset residues (fixed atom list)

    integer :: n_center                                        ! number of atoms for calculating com/geometric center for atom selection (which is n_user_center > n_user_core > n_user_qm)
    integer, dimension(:), pointer :: center => null()         ! array of atoms for calculating com/geometric center for atom selection (which is user_center > user_core > user_qm)

    integer :: selection_type                                ! atom selection type: 1 - atom-atom distance (default)
                                                             !                      2 - flexible center selection
                                                             !                      3 - fixed center selection

    integer :: center_type                                   ! type of center selection: 1 - center of mass (default)
                                                             !                           2 - geometric center

    integer :: initial_selection_type                        ! type of initial selection based on the chosen atom selection type: -1 - inner sphere selection
                                                             !                                                                     0 - middle sphere selection (default)
                                                             !                                                                     1 - outer sphere selection

    integer :: system                                        ! current system to be calculated: 1 - extended (with buffer), 2 - core (without buffer) 
    integer :: qmstep                                        ! current number of md step: 0 for not ABF QM/MM calculation
    integer :: maxqmstep                                     ! maximum number of ABF QM/MM steps = nstlim when ABF QM/MM is active

    integer :: class_charge                                  ! total classical charge of system
    integer :: abfcharge                                     ! total charge of current qm system
    integer :: corecharge                                    ! total charge of core atoms
    integer :: qmcharge                                      ! total charge of qm atoms
    integer :: buffercharge                                  ! total charge of buffer atoms

    character(len=256) :: cut_bond_list_file                 ! list file of breakable bonds, if not specified then default is used
    character(len=256) :: oxidation_number_list_file         ! list file of oxidation numbers
    integer :: max_bonds_per_atom                            ! maximum coordination number (default is 4)
    integer :: n_max_recursive                               ! maximum number of atoms can be recursively labelled (default is 10000)
    integer :: n_recursive                                   ! actual number of recursively labelled atom

    integer :: mom_cons_type                                 ! type of force correction for momentum conservation: 0 - no force correction is applied
                                                             !                                                     1 - equal acceleration on each atom (default)
                                                             !                                                     2 - equal force on each atom
                                                             !                                                    -1 - proportional acceleration each atom
                                                             !                                                    -2 - proportional force on each atom

    integer :: mom_cons_region                               ! region where force correction is applied for momentum conservation: 0 - only core region 
                                                             !                                                                     1 - only core+qm region (default)
                                                             !                                                                     2 - only core+qm+buffer region
                                                             !                                                                     3 - all atom

    integer :: fix_atom_list                                 ! atom list defined in mask is: 0 - not fixed (default), 1 - fixed
    integer :: solvent_atom_number                           ! number of atoms in each solvent molecule (default is 3, water)

    character(len=256) :: read_idrst_file                    ! id restart file name to start simulation containing atom assignment
    character(len=256) :: write_idrst_file                   ! id restart file name for printing atom assignment
    integer :: ntwidrst                                      ! print id restart frequency (in the last md step automatically printed)

    character(len=256) :: pdb_file                           ! output pdb file name with atom assignment
    integer :: ntwpdb                                        ! print pdb frequency

    _REAL_, dimension(:), pointer :: x => null()             ! positions
    _REAL_, dimension(:), pointer :: v => null()             ! velocities
    _REAL_, dimension(:), pointer :: f => null()             ! final, combined forces
    _REAL_, dimension(:), pointer :: f1 => null()            ! forces of system 1 (extended)
    _REAL_, dimension(:), pointer :: f2 => null()            ! forces of system 2 (core)

    _REAL_  :: r_core_in                                     ! inner core radius around user specified (corrected) core atoms
    _REAL_  :: r_core_out                                    ! outer core radius around user specified (corrected) core atoms
    _REAL_  :: r_qm_in                                       ! inner extended radius around user specified (corrected) qm atoms
    _REAL_  :: r_qm_out                                      ! outer extended radius around user specified (corrected) qm atoms
    _REAL_  :: r_buffer_in                                   ! inner buffer radius around qm atoms (user defined + not user defined)
    _REAL_  :: r_buffer_out                                  ! outer buffer radius around qm atoms (user defined + not user defined)

    _REAL_  :: r_hot_spot_in                                 ! inner radius of hot spot layer
    _REAL_  :: r_hot_spot_out                                ! outer radius of hot spot layer
    _REAL_ , dimension(:), pointer :: r_hot_spot => null()   ! distance of residues from center of hot spot

    integer, dimension(:), pointer :: res_pointers => null()                ! first atom number of each residue
    integer, dimension(:), pointer :: res_atom_number => null()             ! number of atoms of each residue
    character(len=4), dimension(:), pointer :: atom_name => null()          ! array of atom names
    character(len=4), dimension(:), pointer :: res_name_of_atom => null()   ! array of res names
    integer, dimension(:), pointer :: res_id_of_atom => null()              ! array of res ids
    integer, dimension(:), pointer :: numbond => null()                     ! number of unbreakable bonds 
    integer, dimension(:,:), pointer :: listbond => null()                  ! list of unbreakable bonds
    integer, dimension(:), pointer :: cutnumbond => null()                  ! number of breakable bonds 
    integer, dimension(:,:), pointer :: cutlistbond => null()               ! list of breakable bonds
    integer, dimension(:), pointer :: oxidation_number => null()            ! oxidation number of atoms 

    _REAL_ :: min_heavy_mass                                          ! minimum mass of heavy atom (default is 4.0)
    _REAL_, dimension(:), pointer :: mass => null()                   ! atom masses

    _REAL_ :: gamma_ln_qm                                    ! friction coefficient for QM region

    ! diffusion restraint
    _REAL_ :: r_diff_in                                      ! inner radius for diffusion restraint
    _REAL_ :: r_diff_out                                     ! outer radius for diffusion restraint
    _REAL_ :: diff_k                                         ! diffusion restraint force constant
    _REAL_ :: diff_potential                                 ! diffusion restraint potential

    integer :: water_o_n
    _REAL_, dimension(:), pointer :: x0 => null()
    _REAL_, dimension(:), pointer :: dx => null()
    integer, dimension(:), pointer :: water_o_id => null()
    integer, dimension(:), pointer :: water_o_is_qm => null()
    integer, dimension(:), pointer :: water_o_time => null()
    integer, dimension(:), pointer :: ndiffmm => null()
    integer, dimension(:), pointer :: ndiffqm => null()
    _REAL_, dimension(:), pointer :: diffmm => null()
    _REAL_, dimension(:), pointer :: diffqm => null()

     ! general bond parameters
     integer :: numbnd
     _REAL_, dimension(:), pointer :: rk => null()
     _REAL_, dimension(:), pointer :: req => null()

     ! bonds with hydrogens
     integer :: nbonh
     integer, dimension(:), pointer :: iibh => null()
     integer, dimension(:), pointer :: ijbh => null()
     integer, dimension(:), pointer :: icbh => null()

     ! bonds without hydrogens
     integer :: nbona
     integer, dimension(:), pointer :: iiba => null()
     integer, dimension(:), pointer :: ijba => null()
     integer, dimension(:), pointer :: icba => null()

     ! angles with hydrogens
     integer :: ntheth
     integer, dimension(:), pointer :: iith => null() ! i24
     integer, dimension(:), pointer :: ijth => null() ! i26
     integer, dimension(:), pointer :: ikth => null() ! i28
     integer, dimension(:), pointer :: icth => null() ! i30

     ! angles without hydrogens
     integer :: ntheta
     integer, dimension(:), pointer :: iita => null() ! i32
     integer, dimension(:), pointer :: ijta => null() ! i34
     integer, dimension(:), pointer :: ikta => null() ! i36
     integer, dimension(:), pointer :: icta => null() ! i38

     ! dihedrals with hydrogens
     integer :: nphih
     integer, dimension(:), pointer :: iiph => null() ! i40
     integer, dimension(:), pointer :: ijph => null() ! i42
     integer, dimension(:), pointer :: ikph => null() ! i44
     integer, dimension(:), pointer :: ilph => null() ! i46
     integer, dimension(:), pointer :: icph => null() ! i48

     ! dihedrals without hydrogens
     integer :: nphia
     integer, dimension(:), pointer :: iipa => null() ! i50
     integer, dimension(:), pointer :: ijpa => null() ! i52
     integer, dimension(:), pointer :: ikpa => null() ! i54
     integer, dimension(:), pointer :: ilpa => null() ! i56
     integer, dimension(:), pointer :: icpa => null() ! i58

     ! mm charges
     _REAL_, dimension(:), pointer :: charge => null()

     ! charmm Urey-Bradley
     integer :: charmm_nub
     type(chm_ang_ub_struct), dimension(:), pointer :: charmm_ang_ub => null()

     ! charmm improper
     integer :: charmm_nimphi
     type(chm_imp_struct), dimension(:), pointer :: charmm_imp => null()

     ! charmm CMAP
     integer :: cmap_term_count
     integer, dimension(:,:), pointer :: cmap_index => null()

  end type abfqmmm_parameters

  type(abfqmmm_parameters), save :: abfqmmm_param

contains

subroutine abfqmmm_init_param()

  implicit none

   abfqmmm_param%abfqmmm=0
   abfqmmm_param%maxqmstep=1
   abfqmmm_param%qmstep=1
   abfqmmm_param%system=1 

end subroutine abfqmmm_init_param


subroutine abfqmmm_connection_list(natom, nbonh, nbona, iibh, ijbh, iiba, ijba)

  implicit none

  integer, intent(in) :: natom
  integer, intent(in) :: nbonh
  integer, intent(in) :: nbona 
  integer, intent(in) :: iibh(*)
  integer, intent(in) :: ijbh(*)
  integer, intent(in) :: iiba(*)
  integer, intent(in) :: ijba(*)

  integer :: ier=0
  integer :: i, j
  integer :: atom1, atom2

  ! first step: allocate arrays
  allocate(abfqmmm_param%numbond(natom), stat=ier)
  REQUIRE(ier==0)

  allocate(abfqmmm_param%listbond(natom, abfqmmm_param%max_bonds_per_atom), stat=ier)
  REQUIRE(ier==0)

  ! allocate arrays also for breakable bonds  
  allocate(abfqmmm_param%cutnumbond(natom), stat=ier)
  REQUIRE(ier==0)

  allocate(abfqmmm_param%cutlistbond(natom, abfqmmm_param%max_bonds_per_atom), stat=ier)
  REQUIRE(ier==0)

  ! second step: clean the table
  do i=1,natom
   abfqmmm_param%numbond(i)=0
   abfqmmm_param%cutnumbond(i)=0
   do j=1,abfqmmm_param%max_bonds_per_atom
    abfqmmm_param%listbond(i,j)=0
    abfqmmm_param%cutlistbond(i,j)=0
   enddo
  enddo

  ! third step: build the table from heavy atom - heavy atom bond list
  do i=1,nbona
   atom1=iiba(i)/3+1
   atom2=ijba(i)/3+1
   abfqmmm_param%numbond(atom1)=abfqmmm_param%numbond(atom1)+1
   abfqmmm_param%numbond(atom2)=abfqmmm_param%numbond(atom2)+1
   abfqmmm_param%listbond(atom1,abfqmmm_param%numbond(atom1))=atom2
   abfqmmm_param%listbond(atom2,abfqmmm_param%numbond(atom2))=atom1
  enddo

  ! fourth step: build the table from heavy atom - hydrogen atom bond list
  do i=1,nbonh
   atom1=iibh(i)/3+1
   atom2=ijbh(i)/3+1
   abfqmmm_param%numbond(atom1)=abfqmmm_param%numbond(atom1)+1
   abfqmmm_param%numbond(atom2)=abfqmmm_param%numbond(atom2)+1
   abfqmmm_param%listbond(atom1,abfqmmm_param%numbond(atom1))=atom2
   abfqmmm_param%listbond(atom2,abfqmmm_param%numbond(atom2))=atom1
  enddo
  
end subroutine abfqmmm_connection_list


subroutine abfqmmm_cut_bond_list

  implicit none

  type cut_bonds_type 
   character(len=5) :: spec
   character(len=4) :: type1
   character(len=4) :: type2
   integer :: index1
   integer :: index2
  end type cut_bonds_type

  type(cut_bonds_type), dimension(:), allocatable :: cut_bonds

  integer :: i, j, k, l, m, ios, ier

  integer :: n_cut_bond 

  character(len=256) :: line
  character(len=1) :: first_char
  character(len=4) :: type1, type2
  integer :: index1, index2
  character(len=3) :: direction
  character(len=6) :: test1, test2 

  logical :: found

  if(abfqmmm_param%cut_bond_list_file == '') then
   write(6,'(a45)') 'WARNING: cut_bond_list_file is not specified!'
   write(6,*)
   return
  endif

  open(unit=1983,file=abfqmmm_param%cut_bond_list_file,status='old',iostat=ios)

  if(ios /=0) then
   write(6,*) 'ERROR: cannot open cut_bond_list_file: ', trim(abfqmmm_param%cut_bond_list_file)
   stop
  end if

  n_cut_bond = 0
  do while (ios==0)
   read(unit=1983,fmt='(a256)',iostat=ios) line
   if(ios==0) then
    line=adjustl(line)
    first_char=line(1:1)
    if(first_char /= '!' .and. first_char /= '#') then
     ! we do not know yet if it is an index or type specification, anyway, now we are looking for directions
     read(line,*,iostat=ios) test1, direction, test2
     if(direction == '<=>') then
      n_cut_bond = n_cut_bond + 2
     else if(direction == '=>' .or. direction == '<=') then
      n_cut_bond = n_cut_bond + 1
     endif
    endif
   endif
  enddo

  allocate(cut_bonds(n_cut_bond), stat=ier)
  REQUIRE(ier==0)

  rewind(unit=1983)

  ios=0
  i = 0
  do while (ios==0)
   read(unit=1983,fmt='(a256)',iostat=ios) line
   if(ios==0) then
    line=adjustl(line)
    first_char=line(1:1)
    if(first_char /= '!' .and. first_char /= '#') then
     ! first try if the specification is index based
     read(line,*,iostat=ios) index1, direction, index2
     if(ios==0) then
      if(direction == '<=>') then
       i = i + 1
       cut_bonds(i)%spec='index'
       cut_bonds(i)%index1=index1
       cut_bonds(i)%index2=index2
       i = i + 1
       cut_bonds(i)%spec='index'
       cut_bonds(i)%index1=index2
       cut_bonds(i)%index2=index1
      else if(direction == '=>') then
       i = i + 1
       cut_bonds(i)%spec='index'
       cut_bonds(i)%index1=index1
       cut_bonds(i)%index2=index2
      else if(direction == '<=') then
       i = i + 1
       cut_bonds(i)%spec='index'
       cut_bonds(i)%index1=index2
       cut_bonds(i)%index2=index1
      end if
     else
      ! second try if the specification is atom type based
      ios=0
      read(line,*,iostat=ios) type1, direction, type2
      if(ios==0) then
       if(direction == '<=>') then
        i = i + 1
        cut_bonds(i)%spec='type'
        cut_bonds(i)%type1=type1
        cut_bonds(i)%type2=type2
        i = i + 1
        cut_bonds(i)%spec='type'
        cut_bonds(i)%type1=type2
        cut_bonds(i)%type2=type1
       else if(direction == '=>') then
        i = i + 1
        cut_bonds(i)%spec='type'
        cut_bonds(i)%type1=type1
        cut_bonds(i)%type2=type2
       else if(direction == '<=') then
        i = i + 1
        cut_bonds(i)%spec='type'
        cut_bonds(i)%type1=type2
        cut_bonds(i)%type2=type1
       end if
      else
       write(6,*) 'ERROR: wrong specification of cut bond: ', line
       stop
      endif
     endif
    endif
   endif
  enddo
     
  write(6,*)
  write(6,'(a)') 'Breakable bonds'
  write(6,'(a)') 'QM atom type/index  =>  MM atom type/index'
  write(6,'(a)') '------------------      ------------------'
  do i=1, n_cut_bond
   if(cut_bonds(i)%spec=='type') then
    write(6,'(6x,a4,20x,a4)') adjustr(cut_bonds(i)%type1), adjustr(cut_bonds(i)%type2)
   else if(cut_bonds(i)%spec=='index') then
    write(6,'(4x,i6,18x,i6)') cut_bonds(i)%index1, cut_bonds(i)%index2
   end if
  end do

  ! now delete breakable QM-MM bonds from abfqmmm_param%listbond 
  do i=1, n_cut_bond
   if(cut_bonds(i)%spec=='type') then
    do j=1, abfqmmm_param%natom
     if(abfqmmm_param%atom_name(j) == cut_bonds(i)%type1) then
      k=1
      do while(k<=abfqmmm_param%numbond(j))
       if(abfqmmm_param%atom_name(abfqmmm_param%listbond(j,k)) == cut_bonds(i)%type2) then
        m = abfqmmm_param%listbond(j,k)
        do l=k, abfqmmm_param%numbond(j)-1
         abfqmmm_param%listbond(j,l) = abfqmmm_param%listbond(j,l+1)
        enddo
        abfqmmm_param%listbond(j,abfqmmm_param%numbond(j)) = 0
        abfqmmm_param%numbond(j) = abfqmmm_param%numbond(j)-1
        ! store breakable bonds for printing out link atom
        found=.false.
        do l=1, abfqmmm_param%cutnumbond(j)
         if(abfqmmm_param%cutlistbond(j,l) == m) then
          found=.true.
          exit
         end if
        end do
        if(.not. found) then
         abfqmmm_param%cutnumbond(j) = abfqmmm_param%cutnumbond(j) + 1
         abfqmmm_param%cutlistbond(j,abfqmmm_param%cutnumbond(j)) = m
         abfqmmm_param%cutnumbond(m) = abfqmmm_param%cutnumbond(m) + 1
         abfqmmm_param%cutlistbond(m,abfqmmm_param%cutnumbond(m)) = j
        end if
       else
        k=k+1
       end if
      enddo
     endif
    enddo
   else if(cut_bonds(i)%spec=='index') then
    j=cut_bonds(i)%index1
    k=1
    do while(k<=abfqmmm_param%numbond(j))
     if(abfqmmm_param%listbond(j,k) == cut_bonds(i)%index2) then
      m = abfqmmm_param%listbond(j,k)
      do l=k, abfqmmm_param%numbond(j)-1
       abfqmmm_param%listbond(j,l) = abfqmmm_param%listbond(j,l+1)
      enddo
      abfqmmm_param%listbond(j,abfqmmm_param%numbond(j)) = 0
      abfqmmm_param%numbond(j) = abfqmmm_param%numbond(j)-1
      ! store breakable bonds for printing out link atom
      found=.false.
      do l=1, abfqmmm_param%cutnumbond(j)
       if(abfqmmm_param%cutlistbond(j,l) == m) then
        found=.true.
        exit
       end if
      end do
      if(.not. found) then
       abfqmmm_param%cutnumbond(j) = abfqmmm_param%cutnumbond(j) + 1
       abfqmmm_param%cutlistbond(j,abfqmmm_param%cutnumbond(j)) = m
       abfqmmm_param%cutnumbond(m) = abfqmmm_param%cutnumbond(m) + 1
       abfqmmm_param%cutlistbond(m,abfqmmm_param%cutnumbond(m)) = j
      end if
      exit
     else
      k=k+1
     end if
    enddo
   endif
  enddo

  close(1983)

end subroutine abfqmmm_cut_bond_list


subroutine abfqmmm_oxidation_number_list

  implicit none

  type oxidation_number_type
   character(len=4) :: atomtype
   character(len=4) :: resname
   integer          :: resid
   integer          :: atomid
   integer          :: oxnum
  end type oxidation_number_type

  type(oxidation_number_type), dimension(:), allocatable :: oxidation_number

  integer :: i, j, k, ier, ios
  integer :: n_oxidation_number, total_charge

  character(len=256) :: line
  character(len=1) :: first_char
  character(len=4) :: atomtype
  character(len=4) :: resname
  integer          :: resid
  integer          :: atomid
  integer :: oxnum
  integer :: star_pos

  allocate(abfqmmm_param%oxidation_number(abfqmmm_param%natom), stat=ier)
  REQUIRE(ier==0)
  
  if(abfqmmm_param%oxidation_number_list_file == '') then
   write(6,'(a53)') 'WARNING: oxidation_number_list_file is not specified!'
   write(6,'(a55)') 'WARNING: all atom oxidation numbers are assumed to be 0'
   write(6,'(a55)') '         and only region specified charges are applied!'
   write(6,*)

   do i=1,abfqmmm_param%natom
    abfqmmm_param%oxidation_number(i)=0
   end do 
   return
  endif

  open(unit=1985,file=abfqmmm_param%oxidation_number_list_file,status='old',iostat=ios)

  if(ios /=0 ) then
   write(6,*) 'ERROR: cannot open oxidation_number_list_file: ', trim(abfqmmm_param%oxidation_number_list_file)
   stop
  end if

  n_oxidation_number = 0
  do while (ios==0)
   read(unit=1985,fmt='(a256)',iostat=ios) line
   if(ios==0) then
    line=adjustl(line)
    first_char=line(1:1)
    if(first_char /= '!' .and. first_char /= '#') then
     n_oxidation_number = n_oxidation_number + 1
    endif
   endif
  enddo

  allocate(oxidation_number(n_oxidation_number), stat=ier)
  REQUIRE(ier==0)

  rewind(unit=1985)

  ios=0
  i = 0
  do while (ios==0)
   read(unit=1985,fmt='(a256)',iostat=ios) line
   if(ios==0) then
    line=adjustl(line)
    first_char=line(1:1)
    if(first_char /= '!' .and. first_char /= '#') then
     read(line,*,iostat=ios) resid, atomtype, oxnum
     if(ios == 0 .and. resid /= 0) then
      if(resid <= 0 .or. resid > abfqmmm_param%nres) then
       write(6,*) 'ERROR: oxidation number cannot be identified from this line:'
       write(6,*) line
       stop
      end if
      do j=1,i
       if(resid == oxidation_number(j)%resid .and. atomtype == oxidation_number(j)%atomtype) then
        write(6,*) 'ERROR: oxidation number of atomtype ', atomtype , ' in residue ', resid, 'is defined more than once!'
        stop
       end if
      end do
      i=i+1
      oxidation_number(i)%resid=resid
      oxidation_number(i)%resname=''
      oxidation_number(i)%atomtype=atomtype
      oxidation_number(i)%atomid=0
      oxidation_number(i)%oxnum=oxnum
     else
      ios=0
      read(line,*,iostat=ios) resname, atomid, oxnum
      if(ios == 0 .and. atomid /= 0) then
       if(resname /= 'atom' .or. atomid <= 0 .or. atomid > abfqmmm_param%natom) then
        write(6,*) 'ERROR: oxidation number cannot be identified from this line:'
        write(6,*) line
        stop
       end if
       do j=1,i
        if(resname == oxidation_number(j)%resname .and. atomid == oxidation_number(j)%atomid) then
         write(6,*) 'ERROR: oxidation number of atom ', atomid, 'is defined more than once!'
         stop
        end if
       end do
       i=i+1
       oxidation_number(i)%resname=resname
       oxidation_number(i)%resid=0
       oxidation_number(i)%atomtype=''
       oxidation_number(i)%atomid=atomid
       oxidation_number(i)%oxnum=oxnum
      else
       ios=0
       read(line,*,iostat=ios) resname, atomtype, oxnum
       if(ios /= 0) then
        write(6,*) 'ERROR: oxidation number cannot be identified from this line:'
        write(6,*) line
        stop
       end if
       do j=1,i
       if(resname == oxidation_number(j)%resname .and. atomtype == oxidation_number(j)%atomtype) then
         write(6,*) 'ERROR: oxidation number of atomtype ', atomtype , ' in residue type ', resname, 'is defined more than once!'
         stop
        end if
       end do
       i=i+1
       oxidation_number(i)%resname=resname
       oxidation_number(i)%resid=0
       oxidation_number(i)%atomtype=atomtype
       oxidation_number(i)%atomid=0
       oxidation_number(i)%oxnum=oxnum
      end if
     end if
    end if
   end if
  end do

  write(6,*)
  write(6,'(a)') 'Oxidation numbers'
  write(6,'(a)') 'Residue type/id      Atom type/id       Ox. number '
  write(6,'(a)') '---------------      ------------      ------------'
  do i=1, n_oxidation_number
   if(oxidation_number(i)%resid==0) then
    if(oxidation_number(i)%atomid==0) then
     write(6,'(4x,a6,15x,a4,14x,i4)') oxidation_number(i)%resname, oxidation_number(i)%atomtype, oxidation_number(i)%oxnum
    else
     write(6,'(4x,a6,11x,i8,14x,i4)') oxidation_number(i)%resname, oxidation_number(i)%atomid, oxidation_number(i)%oxnum
    end if
   else
    write(6,'(4x,i6,15x,a4,14x,i4)') oxidation_number(i)%resid, oxidation_number(i)%atomtype, oxidation_number(i)%oxnum
   end if
  end do

  ! now set up abfqmmm_param%oxidation_number
  do i=1, abfqmmm_param%natom
   abfqmmm_param%oxidation_number(i) = 0
   ! first: go over the 'all' residues and 'X*'->'XY*'->'XYZ*' atom types
   do k=2,4
    do j=1, n_oxidation_number
     star_pos=scan(oxidation_number(j)%atomtype,'*')
     if(oxidation_number(j)%resname == 'all' .and. star_pos == k) then
      if(abfqmmm_param%atom_name(i)(1:(star_pos-1)) == oxidation_number(j)%atomtype(1:(star_pos-1))) then
       abfqmmm_param%oxidation_number(i) = oxidation_number(j)%oxnum
       exit
      end if
     end if
    end do
   end do
   ! second: go over the 'all' residues and specific atom type 
   do j=1, n_oxidation_number
    star_pos=scan(oxidation_number(j)%atomtype,'*')
    if(oxidation_number(j)%resname == 'all' .and. star_pos == 0) then
     if(abfqmmm_param%atom_name(i) == oxidation_number(j)%atomtype) then
      abfqmmm_param%oxidation_number(i) = oxidation_number(j)%oxnum
      exit
     end if
    end if
   end do
   ! third: go over specific residue name and 'X*'->'XY*'->'XYZ*' atom types
   do k=2,4
    do j=1, n_oxidation_number
     star_pos=scan(oxidation_number(j)%atomtype,'*')
     if(oxidation_number(j)%resname /= 'all' .and. oxidation_number(j)%resid == 0 &
        .and. oxidation_number(j)%atomid == 0 .and. star_pos == k) then
      if(abfqmmm_param%atom_name(i)(1:(star_pos-1)) == oxidation_number(j)%atomtype(1:(star_pos-1)) .and. &
         abfqmmm_param%res_name_of_atom(i) == oxidation_number(j)%resname) then
        abfqmmm_param%oxidation_number(i) = oxidation_number(j)%oxnum
        exit
      end if
     end if
    end do
   end do
   ! fourth: go over specific residue name and specific atom type
   do j=1, n_oxidation_number
    star_pos=scan(oxidation_number(j)%atomtype,'*')
    if(oxidation_number(j)%resname /= 'all' .and. oxidation_number(j)%resid == 0 &
       .and. oxidation_number(j)%atomid == 0 .and. star_pos == 0) then
     if(abfqmmm_param%atom_name(i) == oxidation_number(j)%atomtype .and. &
        abfqmmm_param%res_name_of_atom(i) == oxidation_number(j)%resname) then
      abfqmmm_param%oxidation_number(i) = oxidation_number(j)%oxnum
      exit
     end if
    end if
   end do
   ! fifth: given residue index and 'X*'->'XY*'->'XYZ*' atom types
   do k=2,4
    do j=1, n_oxidation_number
     star_pos=scan(oxidation_number(j)%atomtype,'*')
     if(oxidation_number(j)%resid /= 0 .and. star_pos == k) then
      if(abfqmmm_param%atom_name(i)(1:(star_pos-1)) == oxidation_number(j)%atomtype(1:(star_pos-1)) .and. &
         abfqmmm_param%res_id_of_atom(i) == oxidation_number(j)%resid) then
       abfqmmm_param%oxidation_number(i) = oxidation_number(j)%oxnum
       exit
      end if
     end if
    end do
   end do
   ! sixth: given residue index and specific atom type
   do j=1, n_oxidation_number
    star_pos=scan(oxidation_number(j)%atomtype,'*')
    if(oxidation_number(j)%resid /= 0 .and. star_pos == 0) then
     if(abfqmmm_param%atom_name(i) == oxidation_number(j)%atomtype .and. &
        abfqmmm_param%res_id_of_atom(i) == oxidation_number(j)%resid) then
      abfqmmm_param%oxidation_number(i) = oxidation_number(j)%oxnum
      exit
     end if
    end if
   end do
   ! seventh: given atom index
   do j=1, n_oxidation_number
    if(oxidation_number(j)%atomid /= 0) then
     if(i == oxidation_number(j)%atomid) then
      abfqmmm_param%oxidation_number(i) = oxidation_number(j)%oxnum
      exit
     end if
    end if
   end do
  end do 

  total_charge = 0
  do i=1, abfqmmm_param%natom
   total_charge = total_charge + abfqmmm_param%oxidation_number(i)
  end do

  write(6,'(a55,i4)') 'Total charge of system according to oxidation numbers: ', total_charge

  close(1985)

end subroutine abfqmmm_oxidation_number_list


subroutine abfqmmm_setup(natom,nres,res_pointers,atom_name,res_name,mass,nbonh,nbona,iibh,ijbh,iiba,ijba)

  use qmmm_module, only : qmmm_nml, qmmm_struct

  implicit none

  integer, intent(in) :: natom
  integer, intent(in) :: nres
  integer, intent(in) :: res_pointers(nres+1)
  character(len=4), intent(in) :: atom_name(natom)
  character(len=4), intent(in) :: res_name(nres)
  _REAL_, intent(in) :: mass(natom)
  integer, intent(in) :: nbonh
  integer, intent(in) :: nbona
  integer, intent(in) :: iibh(*)
  integer, intent(in) :: ijbh(*)
  integer, intent(in) :: iiba(*)
  integer, intent(in) :: ijba(*)
  integer :: i,j,k,l
  integer :: ier=0
  character(10) :: abftag

  abfqmmm_param%abfqmmm = qmmm_struct%abfqmmm
  abfqmmm_param%hot_spot = qmmm_struct%hot_spot
  abfqmmm_param%natom = natom
  abfqmmm_param%nres = nres
  abfqmmm_param%r_core_in = qmmm_struct%r_core_in
  abfqmmm_param%r_core_out = qmmm_struct%r_core_out
  abfqmmm_param%r_qm_in = qmmm_struct%r_qm_in
  abfqmmm_param%r_qm_out = qmmm_struct%r_qm_out
  abfqmmm_param%r_buffer_in = qmmm_struct%r_buffer_in
  abfqmmm_param%r_buffer_out = qmmm_struct%r_buffer_out

  abfqmmm_param%corecharge = qmmm_nml%corecharge
  abfqmmm_param%qmcharge = qmmm_nml%qmcharge
  abfqmmm_param%buffercharge = qmmm_nml%buffercharge

  abfqmmm_param%cut_bond_list_file = qmmm_struct%cut_bond_list_file
  abfqmmm_param%oxidation_number_list_file = qmmm_struct%oxidation_number_list_file

  abfqmmm_param%mom_cons_type = qmmm_struct%mom_cons_type
  abfqmmm_param%mom_cons_region = qmmm_struct%mom_cons_region

  abfqmmm_param%selection_type = qmmm_struct%selection_type
  abfqmmm_param%center_type = qmmm_struct%center_type
  abfqmmm_param%initial_selection_type = qmmm_struct%initial_selection_type

  abfqmmm_param%gamma_ln_qm = qmmm_struct%gamma_ln_qm
! abfqmmm_param%r_diff_in = qmmm_struct%r_diff_in
! abfqmmm_param%r_diff_out = qmmm_struct%r_diff_out
! abfqmmm_param%diff_k = qmmm_struct%diff_k

  abfqmmm_param%fix_atom_list = qmmm_struct%fix_atom_list
  abfqmmm_param%solvent_atom_number = qmmm_struct%solvent_atom_number

  abfqmmm_param%max_bonds_per_atom = qmmm_struct%max_bonds_per_atom
  abfqmmm_param%n_max_recursive = qmmm_struct%n_max_recursive

  abfqmmm_param%min_heavy_mass = qmmm_struct%min_heavy_mass

  abfqmmm_param%read_idrst_file = qmmm_struct%read_idrst_file
  abfqmmm_param%write_idrst_file = qmmm_struct%write_idrst_file
  abfqmmm_param%ntwidrst = qmmm_struct%ntwidrst

  abfqmmm_param%ntwpdb = qmmm_struct%ntwpdb
  abfqmmm_param%pdb_file = qmmm_struct%pdb_file

  allocate(abfqmmm_param%x(3*natom), stat=ier)
  REQUIRE(ier==0)

  allocate(abfqmmm_param%id(natom), stat=ier)
  REQUIRE(ier==0)

  allocate(abfqmmm_param%id_orig(natom), stat=ier)
  REQUIRE(ier==0)

  allocate(abfqmmm_param%isqm(natom), stat=ier)
  REQUIRE(ier==0)

  allocate(abfqmmm_param%atom_name(natom), stat=ier)
  REQUIRE(ier==0)

  do i=1, natom
   abfqmmm_param%atom_name(i)=atom_name(i)
  end do

  allocate(abfqmmm_param%res_pointers(nres+1), stat=ier)
  REQUIRE(ier==0)

  do i=1, nres+1
   abfqmmm_param%res_pointers(i) = res_pointers(i)
  end do

  allocate(abfqmmm_param%res_atom_number(nres), stat=ier)
  REQUIRE(ier==0)

  do i=1, nres
   abfqmmm_param%res_atom_number(i) = abfqmmm_param%res_pointers(i+1) - abfqmmm_param%res_pointers(i)
  end do

  allocate(abfqmmm_param%res_name_of_atom(natom), stat=ier)
  REQUIRE(ier==0)

  allocate(abfqmmm_param%res_id_of_atom(natom), stat=ier)
  REQUIRE(ier==0)

  do i=1, natom
   do j=1, nres
    if(i >= res_pointers(j)) then
     abfqmmm_param%res_name_of_atom(i) = res_name(j)
     abfqmmm_param%res_id_of_atom(i) = j
    end if
   end do
  end do

  allocate(abfqmmm_param%mass(natom), stat=ier)
  REQUIRE(ier==0)

  do i=1, natom
   abfqmmm_param%mass(i)=mass(i)
  end do

  ! subset specifications

  ! => breakable bond based atom selection 
  if(abfqmmm_param%fix_atom_list == 0) then
   ! core atoms subset if user specified
   if(qmmm_struct%core_nsubset > 0) then
    abfqmmm_param%n_subset_core = qmmm_struct%core_nsubset
    allocate(abfqmmm_param%subset_core(abfqmmm_param%n_subset_core), stat=ier)
    REQUIRE(ier==0)
    do i=1, abfqmmm_param%n_subset_core
     abfqmmm_param%subset_core(i) = qmmm_struct%core_subsetatoms(i)
    end do
   ! by default all atoms are in core subset
   else
    abfqmmm_param%n_subset_core = abfqmmm_param%natom
    allocate(abfqmmm_param%subset_core(abfqmmm_param%n_subset_core), stat=ier)
    REQUIRE(ier==0)
    do i=1, abfqmmm_param%n_subset_core
     abfqmmm_param%subset_core(i) = i
    end do
   end if
   ! qm atoms subset if user specified
   if(qmmm_struct%qm_nsubset > 0) then
    abfqmmm_param%n_subset_qm = qmmm_struct%qm_nsubset
    allocate(abfqmmm_param%subset_qm(abfqmmm_param%n_subset_qm), stat=ier)
    REQUIRE(ier==0)
    do i=1, abfqmmm_param%n_subset_qm
     abfqmmm_param%subset_qm(i) = qmmm_struct%qm_subsetatoms(i)
    end do
   ! by default all atoms are in qm subset
   else
    abfqmmm_param%n_subset_qm = abfqmmm_param%natom
    allocate(abfqmmm_param%subset_qm(abfqmmm_param%n_subset_qm), stat=ier)
    REQUIRE(ier==0)
    do i=1, abfqmmm_param%n_subset_qm
     abfqmmm_param%subset_qm(i) = i
    end do
   end if
   ! buffer atoms subset if user specified
   if(qmmm_struct%buffer_nsubset > 0) then
    abfqmmm_param%n_subset_buffer = qmmm_struct%buffer_nsubset
    allocate(abfqmmm_param%subset_buffer(abfqmmm_param%n_subset_buffer), stat=ier)
    REQUIRE(ier==0)
    do i=1, abfqmmm_param%n_subset_buffer
     abfqmmm_param%subset_buffer(i) = qmmm_struct%buffer_subsetatoms(i)
    end do
   ! by default all atoms are in buffer subset
   else
    abfqmmm_param%n_subset_buffer = abfqmmm_param%natom
    allocate(abfqmmm_param%subset_buffer(abfqmmm_param%n_subset_buffer), stat=ier)
    REQUIRE(ier==0)
    do i=1, abfqmmm_param%n_subset_buffer
     abfqmmm_param%subset_buffer(i) = i
    end do
   end if
  ! => fixed atom list based selection
  else
   ! core residues subset if user specified
   if(qmmm_struct%core_nsubset > 0) then
    do i=1, qmmm_struct%core_nsubset 
     qmmm_struct%core_subsetatoms(i) = abfqmmm_param%res_id_of_atom(i) 
    end do
    abfqmmm_param%n_subset_core = 1
    j=qmmm_struct%core_subsetatoms(1)
    do i=2, qmmm_struct%core_nsubset
     if(qmmm_struct%core_subsetatoms(i) /= j) then
      j=qmmm_struct%core_subsetatoms(i)
      abfqmmm_param%n_subset_core = abfqmmm_param%n_subset_core + 1
     end if
    end do
    allocate(abfqmmm_param%subset_core(abfqmmm_param%n_subset_core), stat=ier)
    abfqmmm_param%subset_core(1) = qmmm_struct%core_subsetatoms(1)
    j=1
    do i=2, qmmm_struct%core_nsubset
     if(qmmm_struct%core_subsetatoms(i) /= abfqmmm_param%subset_core(j)) then
      j=j+1
      abfqmmm_param%subset_core(j) = qmmm_struct%core_subsetatoms(i)
     end if
    end do
   ! by default all residues are in core subset
   else
    abfqmmm_param%n_subset_core = abfqmmm_param%nres
    allocate(abfqmmm_param%subset_core(abfqmmm_param%n_subset_core), stat=ier)
    do i=1, abfqmmm_param%n_subset_core
     abfqmmm_param%subset_core(i) = i
    end do
   end if
   ! qm residues subset if user specified
   if(qmmm_struct%qm_nsubset > 0) then
    do i=1, qmmm_struct%qm_nsubset
     qmmm_struct%qm_subsetatoms(i) = abfqmmm_param%res_id_of_atom(i)
    end do
    abfqmmm_param%n_subset_qm = 1
    j=qmmm_struct%qm_subsetatoms(1)
    do i=2, qmmm_struct%qm_nsubset
     if(qmmm_struct%qm_subsetatoms(i) /= j) then
      j=qmmm_struct%qm_subsetatoms(i)
      abfqmmm_param%n_subset_qm = abfqmmm_param%n_subset_qm + 1
     end if
    end do
    allocate(abfqmmm_param%subset_qm(abfqmmm_param%n_subset_qm), stat=ier)
    abfqmmm_param%subset_qm(1) = qmmm_struct%qm_subsetatoms(1)
    j=1
    do i=2, qmmm_struct%qm_nsubset
     if(qmmm_struct%qm_subsetatoms(i) /= abfqmmm_param%subset_qm(j)) then
      j=j+1
      abfqmmm_param%subset_qm(j) = qmmm_struct%qm_subsetatoms(i)
     end if
    end do
   ! by default all residues are in qm subset
   else
    abfqmmm_param%n_subset_qm = abfqmmm_param%nres
    allocate(abfqmmm_param%subset_qm(abfqmmm_param%n_subset_qm), stat=ier)
    do i=1, abfqmmm_param%n_subset_qm
     abfqmmm_param%subset_qm(i) = i
    end do
   end if
   ! buffer residues subset if user specified
   if(qmmm_struct%buffer_nsubset > 0) then
    do i=1, qmmm_struct%buffer_nsubset
     qmmm_struct%buffer_subsetatoms(i) = abfqmmm_param%res_id_of_atom(i)
    end do
    abfqmmm_param%n_subset_buffer = 1
    j=qmmm_struct%buffer_subsetatoms(1)
    do i=2, qmmm_struct%buffer_nsubset
     if(qmmm_struct%buffer_subsetatoms(i) /= j) then
      j=qmmm_struct%buffer_subsetatoms(i)
      abfqmmm_param%n_subset_buffer = abfqmmm_param%n_subset_buffer + 1
     end if
    end do
    allocate(abfqmmm_param%subset_buffer(abfqmmm_param%n_subset_buffer), stat=ier)
    abfqmmm_param%subset_buffer(1) = qmmm_struct%buffer_subsetatoms(1)
    j=1
    do i=2, qmmm_struct%buffer_nsubset
     if(qmmm_struct%buffer_subsetatoms(i) /= abfqmmm_param%subset_buffer(j)) then
      j=j+1
      abfqmmm_param%subset_buffer(j) = qmmm_struct%buffer_subsetatoms(i)
     end if
    end do
   ! by default all residues are in buffer subset
   else
    abfqmmm_param%n_subset_buffer = abfqmmm_param%nres
    allocate(abfqmmm_param%subset_buffer(abfqmmm_param%n_subset_buffer), stat=ier)
    do i=1, abfqmmm_param%n_subset_buffer
     abfqmmm_param%subset_buffer(i) = i
    end do
   end if
  end if

  ! first of all call abfqmmm_connection_list if fix_atom_list == 0

  if(abfqmmm_param%fix_atom_list == 0) call abfqmmm_connection_list(natom, nbonh, nbona, iibh, ijbh, iiba, ijba)

  ! by default in the beginning all atoms are MM 

  do i=1, natom
   abfqmmm_param%id(i) = 7
  end do

  ! now let's use the user specifications

  ! first we set up the user defined buffer atoms
  do i=1, qmmm_struct%buffer_nquant
   abfqmmm_param%id(qmmm_struct%buffer_iqmatoms(i)) = 5
  enddo

  ! second we set up the user defined qm atoms (possible buffer id will be overwritten)
  do i=1, qmmm_struct%nquant
   abfqmmm_param%id(qmmm_struct%iqmatoms(i)) = 3
  enddo

  ! third we set up the user defined core atoms (possible buffer/qm id will be overwritten)
  do i=1, qmmm_struct%core_nquant
   abfqmmm_param%id(qmmm_struct%core_iqmatoms(i)) = 1
  enddo

  call abfqmmm_read_idrst()

  write(6,*)
  write(6,'(a)') '------------------------------'
  if (abfqmmm_param%read_idrst_file == '') then
   write(6,'(a)') 'A. User defined atoms         '
  else
   write(6,'(a)') 'A. Atom ids from idrst file   '
  end if
  write(6,'(a)') '------------------------------'
  write(6,*)
  write(6,'(a,2x,a)') ' MM_NO. ', '   TYPE   '
  write(6,'(a,2x,a)') '--------', '----------'
  do i=1, natom 
   if(abfqmmm_param%id(i) == 1) then
    abftag='   core   '
   else if(abfqmmm_param%id(i) == 2) then
    abftag=' ext-core '
   else if(abfqmmm_param%id(i) == 3) then
    abftag='    qm    '
   else if(abfqmmm_param%id(i) == 4) then
    abftag='  ext-qm  '
   else if(abfqmmm_param%id(i) == 5) then
    abftag='  buffer  '
   else if(abfqmmm_param%id(i) == 6) then
    abftag='ext-buffer'
   end if
   if(abfqmmm_param%id(i) /= 7) write(6,'(i8,2x,10a)') i, abftag
  end do
  write(6,'(a)') '-----------------------------'
  write(6,*)
  write(6,'(a22,f5.2,a4)') 'Core inner radius:    ', abfqmmm_param%r_core_in, ' [A]'
  write(6,'(a22,f5.2,a4)') 'Core outer radius:    ', abfqmmm_param%r_core_out, ' [A]'
  write(6,'(a22,f5.2,a4)') 'Quantum inner radius: ', abfqmmm_param%r_qm_in, ' [A]'
  write(6,'(a22,f5.2,a4)') 'Quantum outer radius: ', abfqmmm_param%r_qm_out, ' [A]'
  write(6,'(a22,f5.2,a4)') 'Buffer inner radius:  ', abfqmmm_param%r_buffer_in, ' [A]'
  write(6,'(a22,f5.2,a4)') 'Buffer outer radius:  ', abfqmmm_param%r_buffer_out, ' [A]'
  write(6,*)

  if(abfqmmm_param%hot_spot == 1) then
   if (abfqmmm_param%r_buffer_in .ne. abfqmmm_param%r_buffer_out) then
    write(6,'(a)') 'ERROR: Hot spot is active: r_buffer_in must be equal to r_buffer_out!'
    stop
   else
    write(6,'(a43,f5.2,a4)') 'Hot spot width (r_buffer_out = r_buffer_in): ', &
                             abfqmmm_param%r_buffer_out, ' [A]'
    write(6,*)
   end if
  end if

  ! selection type
  select case(abfqmmm_param%selection_type)
   case(1)
    write(6,'(a)') 'Selection type is atom-atom distance selection'
   case(2)
    write(6,'(a)') 'Selection type is flexible-sphere center selection'
   case(3)
    write(6,'(a)') 'Selection type is fixed-sphere center selection'
  end select
  write(6,*)

  ! center type
  select case(abfqmmm_param%center_type)
   case(1)
    write(6,'(a)') 'Center type is center of mass'
   case(2)
    write(6,'(a)') 'Center type is geometric center'
  end select
  write(6,*)

  ! initial selection type
  if (abfqmmm_param%read_idrst_file == '') then
   select case(abfqmmm_param%initial_selection_type)
    case(-1)
     write(6,'(a)') 'Initial selection type is inner sphere selection'
    case(0)
     write(6,'(a)') 'Initial selection type is middle sphere selection'
    case(1)
     write(6,'(a)') 'Initial selection type is outer sphere selection'
   end select
  end if
  write(6,*)

  if(abfqmmm_param%gamma_ln_qm == 0.0d0) then
   write(6,'(a)') 'Fricition coefficient for QM region was not specified: gamma_ln_qm = gamma_ln'
  else if(abfqmmm_param%gamma_ln_qm > 0.0d0) then
   write(6,'(a40,f6.2,a8)') 'Fricition coefficient for QM region is: ', abfqmmm_param%gamma_ln_qm, ' ps^(-1)'
  else
   write(6,'(a)') 'Fricition coefficient for QM region must be positive!'
   stop
  end if
  write(6,*)
  ! diffusion restraint
! write(6,'(a30,f8.3,a15)') 'Diffusion restraint constant: ', abfqmmm_param%diff_k, ' [kcal/mol]'
! if(abfqmmm_param%r_diff_in == 0.0d0) then
!  write(6,'(a)') 'Diffusion restraint inner radius is calculated according to the actual qm size'
! else
!  write(6,'(a34,f4.1,a4)') 'Diffusion restraint inner radius: ', abfqmmm_param%r_diff_in, ' [A]'
! end if
! if(abfqmmm_param%r_diff_out == 0.0d0) then
!  write(6,'(a)') 'Diffusion restraint outer radius is calculated according to the actual buffer size'
! else
!  write(6,'(a34,f4.1,a4)') 'Diffusion restraint outer radius: ', abfqmmm_param%r_diff_out, ' [A]'
! end if
! write(6,*)

  ! force correction type for momentum consercation
  select case(abfqmmm_param%mom_cons_type)
   case(0)
    write(6,'(a)') 'WARNING: no force correction is applied! Momentum conservation is not held!'
   case(1)
    write(6,'(a)') 'Equal acceleration is applied for momentum conservation'
   case(2)
    write(6,'(a)') 'Equal force is applied for momentum conservation'
   case(-1)
    write(6,'(a)') 'Proportional acceleration is applied for momentum conservation'
   case(-2)
    write(6,'(a)') 'Proportional force is applied for momentum conservation'
   case default
    write(6,'(a)') 'ERROR: mom_cons_type must be 0, 1, 2, -1 or -2!'
    stop
  end select

  if(abfqmmm_param%mom_cons_type /= 0) then 
   ! force correction region for momentum consercation
   select case(abfqmmm_param%mom_cons_region)
    case(0)
     write(6,'(a)') 'Force correction is distributed on CORE atoms'
    case(1)
     write(6,'(a)') 'Force correction is distributed on CORE+QM atoms'
    case(2)
     write(6,'(a)') 'Force correction is distributed on CORE+QM+BUFFER atoms'
    case(3)
     write(6,'(a)') 'Force correction is distributed only on all atoms'
    case default
     write(6,'(a)') 'ERROR: mom_cons_region must be  0, 1, 2 or 3!'
     stop
   end select
  end if
  write(6,*)

  ! which method are we going to use

  if(abfqmmm_param%fix_atom_list == 0) then
   write(6,'(a)') 'Breakable bond based atom selection is applied:'
   call abfqmmm_cut_bond_list
  else
   write(6,'(a)') 'Fix atom list based selection is applied:'
   write(6,'(a)') '1. fixed atom list for solutes defined in coremask, qmmask and buffermask'
   write(6,'(a)') '2. residue based selection for solvent molecules'
   write(6,*)
   write(6,'(a40,i2)') 'Number of atoms in solvent molecule is: ', abfqmmm_param%solvent_atom_number
  end if

  call abfqmmm_oxidation_number_list

  ! correction of user specifications

  if (abfqmmm_param%read_idrst_file == '') then

   ! start again from assuming all atoms as MM

   do i=1, natom
    abfqmmm_param%id(i) = 7
   end do

   ! first correct the buffer
   if(qmmm_struct%buffer_nquant > 0) then
    do i=1, qmmm_struct%buffer_nquant
     abfqmmm_param%id(qmmm_struct%buffer_iqmatoms(i)) = 5 
    end do

    if(abfqmmm_param%fix_atom_list == 0) then
     abfqmmm_param%n_recursive = 0
     do i=1, qmmm_struct%buffer_nquant 
      call abfqmmm_recursive_labelling(qmmm_struct%buffer_iqmatoms(i),5)
     end do
    else
     do i=1, nres
      call abfqmmm_residue_labelling(i,5)
     end do
    end if
   end if

   ! second correct the qm atoms (possible buffer id will be overwritten)
   do i=1, qmmm_struct%nquant
    abfqmmm_param%id(qmmm_struct%iqmatoms(i)) = 3
   end do 

   if(abfqmmm_param%fix_atom_list == 0) then
    abfqmmm_param%n_recursive = 0
    do i=1, qmmm_struct%nquant
     call abfqmmm_recursive_labelling(qmmm_struct%iqmatoms(i),3)
    end do
   else
    do i=1, nres
     call abfqmmm_residue_labelling(i,3)
    end do 
   end if

   ! third correct the core atoms (possible buffer/qm id will be overwritten)
   if(qmmm_struct%core_nquant > 0) then
    do i=1, qmmm_struct%core_nquant
     abfqmmm_param%id(qmmm_struct%core_iqmatoms(i)) = 1
    end do

    if(abfqmmm_param%fix_atom_list == 0) then
      abfqmmm_param%n_recursive = 0
     do i=1, qmmm_struct%core_nquant
      call abfqmmm_recursive_labelling(qmmm_struct%core_iqmatoms(i),1)
     end do
    else
     do i=1, nres
      call abfqmmm_residue_labelling(i,1)
     end do
    end if
   end if

  end if

  abfqmmm_param%n_user_buffer=0
  abfqmmm_param%n_user_qm=0
  abfqmmm_param%n_user_core=0
  do i=1, natom
   if(abfqmmm_param%id(i) == 5) then
    abfqmmm_param%n_user_buffer=abfqmmm_param%n_user_buffer+1
   else if(abfqmmm_param%id(i) == 3) then
    abfqmmm_param%n_user_qm=abfqmmm_param%n_user_qm+1
   else if(abfqmmm_param%id(i) == 1) then
    abfqmmm_param%n_user_core=abfqmmm_param%n_user_core+1
    abfqmmm_param%n_user_qm=abfqmmm_param%n_user_qm+1
   end if
  end do

  if(abfqmmm_param%n_user_buffer /= 0) then
   allocate(abfqmmm_param%user_buffer(abfqmmm_param%n_user_buffer), stat=ier)
   REQUIRE(ier==0)
  end if

  if(abfqmmm_param%n_user_qm /= 0) then
   allocate(abfqmmm_param%user_qm(abfqmmm_param%n_user_qm), stat=ier)
   REQUIRE(ier==0)
  end if

  if(abfqmmm_param%n_user_core /= 0) then
   allocate(abfqmmm_param%user_core(abfqmmm_param%n_user_core), stat=ier)
   REQUIRE(ier==0)
  end if

  j=1
  k=1
  l=1
  do i=1, natom
   if(abfqmmm_param%id(i) == 5) then
    abfqmmm_param%user_buffer(j) = i
    j=j+1
   else if(abfqmmm_param%id(i) == 3) then
    abfqmmm_param%user_qm(k) = i
    k=k+1
   else if(abfqmmm_param%id(i) == 1) then
    abfqmmm_param%user_core(l) = i
    l=l+1
    abfqmmm_param%user_qm(k) = i
    k=k+1
   end if
  end do

  if( (abfqmmm_param%n_user_core == abfqmmm_param%n_user_qm) .and. (abfqmmm_param%n_user_buffer == 0) .and. &
      (abfqmmm_param%r_qm_out == 0.0d0) .and. (abfqmmm_param%r_buffer_out == 0.0d0) ) then
      write(6,'(a)') 'Adaptive (unbuffered) QM/MM is active: reduced calculation is a full MM'
      write(6,'(a)') 'but it is not taken into account in this type of simulation'
      write(6,*)
  end if

  if (abfqmmm_param%read_idrst_file == '') then

   write(6,'(a)') '------------------------------'
   write(6,'(a)') 'B. Automatic completion       '
   write(6,'(a)') '------------------------------'
   write(6,'(a37,i8)') 'Number of user defined core atoms:   ', abfqmmm_param%n_user_core
   write(6,'(a37,i8)') 'Number of user defined qm atoms:     ', abfqmmm_param%n_user_qm-abfqmmm_param%n_user_core
   write(6,'(a37,i8)') 'Number of user defined buffer atoms: ', abfqmmm_param%n_user_buffer 
   write(6,*)
   write(6,'(a,2x,a)') ' MM_NO. ', '   TYPE   '
   write(6,'(a,2x,a)') '--------', '----------'
   do i=1, natom
    if(abfqmmm_param%id(i) == 1) then
     abftag='   core   '
    else if(abfqmmm_param%id(i) == 3) then
     abftag='    qm    '
    else if(abfqmmm_param%id(i) == 5) then
     abftag='  buffer  '
    end if
    if(abfqmmm_param%id(i) /= 7) write(6,'(i8,2x,10a)') i, abftag
   end do
   write(6,'(a)') '-----------------------------'

   do i=1, natom
     abfqmmm_param%id_orig(i) = abfqmmm_param%id(i)
   end do

  end if

  ! center list for com and geometric center based selections
  ! if defined by user then use the specified set as center list
  if(qmmm_struct%center_nsubset > 0) then
   abfqmmm_param%n_center = qmmm_struct%center_nsubset
   allocate(abfqmmm_param%center(abfqmmm_param%n_center), stat=ier)
   REQUIRE(ier==0)
   do i=1, abfqmmm_param%n_center
    abfqmmm_param%center(i) = qmmm_struct%center_subsetatoms(i)
   end do
  ! if not defined use the user defined core list as center list
  else if(abfqmmm_param%n_user_core > 0) then
   abfqmmm_param%n_center = abfqmmm_param%n_user_core
   allocate(abfqmmm_param%center(abfqmmm_param%n_center), stat=ier)
   REQUIRE(ier==0)
   do i=1, abfqmmm_param%n_center
    abfqmmm_param%center(i) = abfqmmm_param%user_core(i)
   end do
  ! if neither above defined then use user defined qm list as center list
  else
   abfqmmm_param%n_center = abfqmmm_param%n_user_qm
   allocate(abfqmmm_param%center(abfqmmm_param%n_center), stat=ier)
   REQUIRE(ier==0)
   do i=1, abfqmmm_param%n_center
    abfqmmm_param%center(i) = abfqmmm_param%user_qm(i)
   end do
  end if

end subroutine abfqmmm_setup


subroutine abfqmmm_residue_labelling(res, label)

  implicit none

  integer, intent(in) :: res
  integer, intent(in) :: label

  integer :: i, j
  integer :: resstart, resstop

  if(abfqmmm_param%res_atom_number(res) > abfqmmm_param%solvent_atom_number) return
  resstart=abfqmmm_param%res_pointers(res)
  resstop=abfqmmm_param%res_pointers(res+1)-1
  do i=resstart, resstop
   if(abfqmmm_param%id(i) == label) then
    do j=resstart, resstop
     abfqmmm_param%id(j) = label
    end do
    exit
   end if
  end do

end subroutine abfqmmm_residue_labelling


recursive subroutine abfqmmm_recursive_labelling(atom, label)

  implicit none

  integer, intent(in) :: atom
  integer, intent(in) :: label

  integer :: i, j

  if(abfqmmm_param%n_recursive > abfqmmm_param%n_max_recursive) then
   write(6,*) 'ERROR: number of recursive labelling reached the maximum number: ', abfqmmm_param%n_max_recursive
   stop
  end if

  do i=1, abfqmmm_param%numbond(atom)
   j=abfqmmm_param%listbond(atom,i)
   if(abfqmmm_param%id(j) > label) then
    abfqmmm_param%n_recursive = abfqmmm_param%n_recursive + 1
    abfqmmm_param%id(j) = label
    call abfqmmm_recursive_labelling(j, label)
   end if
  enddo

  return

end subroutine abfqmmm_recursive_labelling


recursive subroutine abfqmmm_recursive_find_label_heavy_atom(atom, label)

  implicit none

  integer, intent(in) :: atom
  integer, intent(out) :: label

  integer :: i, j

  if(abfqmmm_param%n_recursive > abfqmmm_param%n_max_recursive) then
   write(6,*) 'ERROR: number of recursive labelling of nonheavy atom reached the maximum number: ', abfqmmm_param%n_max_recursive
   stop
  end if

  do i=1, abfqmmm_param%numbond(atom)
   j=abfqmmm_param%listbond(atom,i)
   if(abfqmmm_param%mass(j) < abfqmmm_param%min_heavy_mass) then
    call abfqmmm_recursive_find_label_heavy_atom(j,label)
   else
    label = abfqmmm_param%id(j)
    exit
   end if
  end do

  return

end subroutine abfqmmm_recursive_find_label_heavy_atom


subroutine abfqmmm_update_qmatoms(crd)

  implicit none

#include "box.h"

  _REAL_, intent(inout) :: crd(3*abfqmmm_param%natom)

  _REAL_ :: mindist
  _REAL_ :: center(3), radius
  _REAL_ :: c(3), r

  integer :: n,i,j,k,l, ierr

  _REAL_ :: box_center(3)

  ! save coordinates for printing out trajectory/pdb 
  abfqmmm_param%x(:) = crd(:)

  if(ntb > 0) call iwrap2(abfqmmm_param%n_user_qm, abfqmmm_param%user_qm, crd, box_center)

  ! First: Build core region
  if(abfqmmm_param%n_user_core > 0) then

   ! calculate center and radius if necessary
   select case(abfqmmm_param%selection_type)
    ! atom-atom distance selection
    case(1)
     radius = 0.0d0
    ! flexible sphere selection
    case(2)
     ! calculate the center of center list and radius of sphere of user defined core list
     call abfqmmm_center_and_radius_of_atom_list(abfqmmm_param%center_type, abfqmmm_param%n_center, abfqmmm_param%center, &
                                                 abfqmmm_param%n_user_core, abfqmmm_param%user_core, crd, center, radius)
    ! fixed sphere com selection
    case(3)
     ! calculate the center of center list and radius of sphere of center list
     call abfqmmm_center_and_radius_of_atom_list(abfqmmm_param%center_type, abfqmmm_param%n_center, abfqmmm_param%center, &
                                                 abfqmmm_param%n_center, abfqmmm_param%center, crd, center, radius)
   end select

   if(abfqmmm_param%fix_atom_list == 0) then

    ! => breakable bond based atom selection
    do n=1,abfqmmm_param%n_subset_core
     i = abfqmmm_param%subset_core(n)
     ! cycle if it is already user defined core or hydrogen
     if( (abfqmmm_param%id(i) < 2) .or. (abfqmmm_param%mass(i) < abfqmmm_param%min_heavy_mass) ) cycle
     ! select selection type
     select case(abfqmmm_param%selection_type)
      ! atom-atom distance selection
      case(1)
       ! calculate the shortest distance between atom and user defined core list
       mindist=abfqmmm_atom_group_min_distance(i, abfqmmm_param%n_user_core, abfqmmm_param%user_core, crd)
      ! flexible sphere selection
      case(2)
       ! calculate distance between atom and center of center list
       mindist=abfqmmm_point_atom_dist(center, i, crd)
      ! fixed sphere selection
      case(3)
       ! calculate distance between atom and center of center list
       mindist=abfqmmm_point_atom_dist(center, i, crd)
     end select
     ! label the atom according to distance critera
     if( (abfqmmm_param%read_idrst_file == '') .and. (abfqmmm_param%qmstep == 1) ) then
      ! if first step and not abfqmmm restart calculation
      select case(abfqmmm_param%initial_selection_type)
       ! inner sphere selection
       case(-1)
        call abfqmmm_atom_labelling_by_cutoff(i, 2, mindist, abfqmmm_param%r_core_in+radius, abfqmmm_param%r_core_in+radius)
       ! middle sphere selection
       case(0)
        call abfqmmm_atom_labelling_by_cutoff(i, 2, mindist, (abfqmmm_param%r_core_in+abfqmmm_param%r_core_out)/2.0d0+radius, &
                                                             (abfqmmm_param%r_core_in+abfqmmm_param%r_core_out)/2.0d0+radius)
       ! outer sphere selection
       case(1)
        call abfqmmm_atom_labelling_by_cutoff(i, 2, mindist, abfqmmm_param%r_core_out+radius, abfqmmm_param%r_core_out+radius)
      end select
     else
       call abfqmmm_atom_labelling_by_cutoff(i, 2, mindist, abfqmmm_param%r_core_in+radius, abfqmmm_param%r_core_out+radius) 
     end if
    end do

    ! if labeled as core, recursively label its environment
    do i=1,abfqmmm_param%natom
     if(abfqmmm_param%id(i) /= 2) cycle
     ! cycle if it is hydrogen
     if(abfqmmm_param%mass(i) < abfqmmm_param%min_heavy_mass) cycle
     abfqmmm_param%n_recursive = 0
     call abfqmmm_recursive_labelling(i, 2)
    end do

    ! label hydrogen atoms 
    do i=1,abfqmmm_param%natom
     ! cycle if it is not hydrogen
     if(abfqmmm_param%mass(i) >= abfqmmm_param%min_heavy_mass) cycle
     abfqmmm_param%n_recursive = 0
     call abfqmmm_recursive_find_label_heavy_atom(i, j)
     abfqmmm_param%id(i) = j
    end do

   else

    ! => fixed atom list based selection 
    do n=1,abfqmmm_param%n_subset_core
     i = abfqmmm_param%subset_core(n)
     ! select selection type
     select case(abfqmmm_param%selection_type)
      ! atom-atom distance selection
      case(1)
       ! calculate the shortest distance between res and user defined core list
       mindist=abfqmmm_residue_group_min_distance(i, 2, abfqmmm_param%n_user_core, abfqmmm_param%user_core, crd)
      ! flexible sphere selection
      case(2)
       ! calculate distance between res and center of center list
       mindist=abfqmmm_point_res_dist(center, i, 2, crd)
      ! fixed sphere selection
      case(3)
       ! calculate distance between res and center of center list
       mindist=abfqmmm_point_res_dist(center, i, 2, crd)
      end select
     ! cycle if already user defined core labeled residue
     if(mindist < 0) cycle
     ! label the whole residue according to distance criteria
     if( (abfqmmm_param%read_idrst_file == '') .and. (abfqmmm_param%qmstep == 1) ) then
      ! if first step and not abfqmmm restart calculation
      select case(abfqmmm_param%initial_selection_type)
       ! inner sphere selection
       case(-1)
        call abfqmmm_res_labelling_by_cutoff(i, 2, mindist, abfqmmm_param%r_core_in+radius, abfqmmm_param%r_core_in+radius)
       ! middle sphere selection
       case(0)
        call abfqmmm_res_labelling_by_cutoff(i, 2, mindist, (abfqmmm_param%r_core_in+abfqmmm_param%r_core_out)/2.0d0+radius, &
                                                            (abfqmmm_param%r_core_in+abfqmmm_param%r_core_out)/2.0d0+radius)
       ! outer sphere selection
       case(1)
        call abfqmmm_res_labelling_by_cutoff(i, 2, mindist, abfqmmm_param%r_core_out+radius, abfqmmm_param%r_core_out+radius)
      end select
     else
       call abfqmmm_res_labelling_by_cutoff(i, 2, mindist, abfqmmm_param%r_core_in+radius, abfqmmm_param%r_core_out+radius)
     end if
    end do

   end if

   ! Allocate core array (core exteded) and qm array (core extended + qm user)

   if( associated(abfqmmm_param%core) ) then
     deallocate(abfqmmm_param%core, stat=ierr)
     REQUIRE (ierr == 0)
   end if
   if( associated(abfqmmm_param%qm) ) then
     deallocate(abfqmmm_param%qm, stat=ierr)
     REQUIRE (ierr == 0)
   end if
   abfqmmm_param%n_core=0
   abfqmmm_param%n_qm=0
   do i=1,abfqmmm_param%natom
    if(abfqmmm_param%id(i)<=2) abfqmmm_param%n_core=abfqmmm_param%n_core+1
    if(abfqmmm_param%id(i)<=3) abfqmmm_param%n_qm=abfqmmm_param%n_qm+1
   end do
   if(abfqmmm_param%n_core>0) allocate(abfqmmm_param%core(abfqmmm_param%n_core))
   if(abfqmmm_param%n_qm>0) allocate(abfqmmm_param%qm(abfqmmm_param%n_qm))
   j=1
   k=1
   do i=1,abfqmmm_param%natom
    if(abfqmmm_param%id(i)<=2) then
     abfqmmm_param%core(j)=i
     j=j+1
     abfqmmm_param%qm(k)=i
     k=k+1
    else if(abfqmmm_param%id(i)<=3) then
     abfqmmm_param%qm(k)=i
     k=k+1
    end if
   end do

  else

   if( associated(abfqmmm_param%qm) ) then
     deallocate(abfqmmm_param%qm, stat=ierr)
     REQUIRE (ierr == 0)
   end if
   abfqmmm_param%n_qm=0
   do i=1,abfqmmm_param%natom
    if(abfqmmm_param%id(i)<=3) abfqmmm_param%n_qm=abfqmmm_param%n_qm+1
   end do
   if(abfqmmm_param%n_qm>0) then
    allocate(abfqmmm_param%qm(abfqmmm_param%n_qm))
    k=1
    do i=1,abfqmmm_param%natom
     if(abfqmmm_param%id(i)<=3) then
      abfqmmm_param%qm(k)=i
      k=k+1
     end if
    end do
   end if

  end if

  ! Second: Build qm region

  ! calculate center and radius if necessary
  select case(abfqmmm_param%selection_type)
   ! atom-atom distance selection 
   case(1)
    radius = 0.0d0
   ! flexible sphere selection
   case(2)
    ! calculate com of center list and radius of sphere of (core extended + qm user) list
    call abfqmmm_center_and_radius_of_atom_list(abfqmmm_param%center_type, abfqmmm_param%n_center, abfqmmm_param%center, &
                                                abfqmmm_param%n_qm, abfqmmm_param%qm, crd, center, radius)
   ! fixed sphere selection
   case(3)
    ! calculate center of center list and radius of sphere of center list
    call abfqmmm_center_and_radius_of_atom_list(abfqmmm_param%center_type, abfqmmm_param%n_center, abfqmmm_param%center, &
                                                abfqmmm_param%n_center, abfqmmm_param%center, crd, center, radius) 
    ! if there is user defined core region then radius is elongated by the mean distance of r_core_in and r_core_out
    if(abfqmmm_param%n_user_core > 0) radius = radius + (abfqmmm_param%r_core_in + abfqmmm_param%r_core_out) / 2.0d0 
  end select

  if(abfqmmm_param%fix_atom_list == 0) then

   ! => breakable bond based atom selection
   do n=1,abfqmmm_param%n_subset_qm
    i = abfqmmm_param%subset_qm(n)
    ! cycle if it is already extended core or user defined qm or hydrogen
    if( (abfqmmm_param%id(i) < 4) .or. (abfqmmm_param%mass(i) < abfqmmm_param%min_heavy_mass) ) cycle
    ! select selection type
    select case(abfqmmm_param%selection_type)
     ! atom-atom distance selection
     case(1)
      ! calculate the shortest distance between atom and extended core + user defined qm list
      mindist=abfqmmm_atom_group_min_distance(i, abfqmmm_param%n_qm, abfqmmm_param%qm, crd)
     ! flexible sphere selection
     case(2)
      ! calculate distance between atom and center of center list
      mindist=abfqmmm_point_atom_dist(center, i, crd)
     ! fixed sphere selection
     case(3)
      ! calculate distance between atom and center of center list
      mindist=abfqmmm_point_atom_dist(center, i, crd)
    end select
    ! label the atom according to distance critera
    if( (abfqmmm_param%read_idrst_file == '') .and. (abfqmmm_param%qmstep == 1) ) then
     ! if first step and not abfqmmm restart calculation
     select case(abfqmmm_param%initial_selection_type)
      ! inner sphere selection
      case(-1)
       call abfqmmm_atom_labelling_by_cutoff(i, 4, mindist, abfqmmm_param%r_qm_in+radius, abfqmmm_param%r_qm_in+radius)
      ! middle sphere selection
      case(0)
       call abfqmmm_atom_labelling_by_cutoff(i, 4, mindist, (abfqmmm_param%r_qm_in+abfqmmm_param%r_qm_out)/2.0d0+radius, &
                                                            (abfqmmm_param%r_qm_in+abfqmmm_param%r_qm_out)/2.0d0+radius)
      ! outer sphere selection
      case(1)
       call abfqmmm_atom_labelling_by_cutoff(i, 4, mindist, abfqmmm_param%r_qm_out+radius, abfqmmm_param%r_qm_out+radius)
     end select
    else
      call abfqmmm_atom_labelling_by_cutoff(i, 4, mindist, abfqmmm_param%r_qm_in+radius, abfqmmm_param%r_qm_out+radius)
    end if
   end do

   ! if labeled as qm, recursively label its environment
   do i=1,abfqmmm_param%natom
    if(abfqmmm_param%id(i) /= 4) cycle
    ! cycle if it is hydrogen
    if(abfqmmm_param%mass(i) < abfqmmm_param%min_heavy_mass) cycle
    abfqmmm_param%n_recursive = 0
    call abfqmmm_recursive_labelling(i, 4)
   end do 

   ! label hydrogen atoms 
   do i=1,abfqmmm_param%natom
    ! cycle if it is not hydrogen
    if(abfqmmm_param%mass(i) >= abfqmmm_param%min_heavy_mass) cycle
    abfqmmm_param%n_recursive = 0
    call abfqmmm_recursive_find_label_heavy_atom(i, j)
    abfqmmm_param%id(i) = j
   end do

  else

   ! => fixed atom list based selection 
   do n=1,abfqmmm_param%n_subset_qm
    i = abfqmmm_param%subset_qm(n)
    select case(abfqmmm_param%selection_type)
     ! atom-atom distance selection
     case(1)
      ! calculate the shortest distance between res and extended core + user defined qm list
      mindist=abfqmmm_residue_group_min_distance(i, 4, abfqmmm_param%n_qm, abfqmmm_param%qm, crd)
     ! flexible sphere selection
     case(2)
      ! calculate distance between res and center of center list
      mindist=abfqmmm_point_res_dist(center, i, 4, crd)
     ! fixed sphere selection
     case(3)
      ! calculate distance between res and center of center list
      mindist=abfqmmm_point_res_dist(center, i, 4, crd)
    end select
    ! cycle if already extended core or user defined qm labeled residue
    if(mindist < 0) cycle
    ! label the whole residue according to distance criteria
    if( (abfqmmm_param%read_idrst_file == '') .and. (abfqmmm_param%qmstep == 1) ) then
     ! if first step and not abfqmmm restart calculation
     select case(abfqmmm_param%initial_selection_type)
      ! inner sphere selection
      case(-1)
       call abfqmmm_res_labelling_by_cutoff(i, 4, mindist, abfqmmm_param%r_qm_in+radius, abfqmmm_param%r_qm_in+radius)
      ! middle sphere selection
      case(0)
       call abfqmmm_res_labelling_by_cutoff(i, 4, mindist, (abfqmmm_param%r_qm_in+abfqmmm_param%r_qm_out)/2.0d0+radius, &
                                                           (abfqmmm_param%r_qm_in+abfqmmm_param%r_qm_out)/2.0d0+radius) 
      ! outer sphere selection
      case(1)
       call abfqmmm_res_labelling_by_cutoff(i, 4, mindist, abfqmmm_param%r_qm_out+radius, abfqmmm_param%r_qm_out+radius)
     end select
    else
      call abfqmmm_res_labelling_by_cutoff(i, 4, mindist, abfqmmm_param%r_qm_in+radius, abfqmmm_param%r_qm_out+radius)
    end if
   end do

  end if

  ! Allocate qm array (core extended + qm extended)

  if( associated(abfqmmm_param%qm) ) then
    deallocate(abfqmmm_param%qm, stat=ierr)
    REQUIRE (ierr == 0)
  end if
  abfqmmm_param%n_qm=0
  do i=1,abfqmmm_param%natom
   if(abfqmmm_param%id(i)<=4) abfqmmm_param%n_qm=abfqmmm_param%n_qm+1
  end do
  if(abfqmmm_param%n_qm>0) then
   allocate(abfqmmm_param%qm(abfqmmm_param%n_qm))
   j=1
   do i=1,abfqmmm_param%natom
    if(abfqmmm_param%id(i)<=4) then
     abfqmmm_param%qm(j)=i
     j=j+1
    end if
   end do
  end if

  ! Third: Build buffer region

  ! calculate center and radius if necessary
  select case(abfqmmm_param%selection_type)
   ! atom-atom distance selection 
   case(1)
    radius = 0
   ! flexible sphere selection
   case(2)
    ! calculate center of center list and radius of sphere of (core extended + qm extended) list
    call abfqmmm_center_and_radius_of_atom_list(abfqmmm_param%center_type, abfqmmm_param%n_center, abfqmmm_param%center, &
                                                abfqmmm_param%n_qm, abfqmmm_param%qm, crd, center, radius)
   ! fixed sphere selection
   case(3)
    ! calculate center of center list and radius of sphere of center list
    call abfqmmm_center_and_radius_of_atom_list(abfqmmm_param%center_type, abfqmmm_param%n_center, abfqmmm_param%center, &
                                                abfqmmm_param%n_center, abfqmmm_param%center, crd, center, radius)
    ! if there is user defined core region then radius is elongated by the mean distance of r_core_in and r_core_out
    if(abfqmmm_param%n_user_core > 0) radius = radius + (abfqmmm_param%r_core_in + abfqmmm_param%r_core_out) / 2.0d0
    ! radius is elongated by the mean distance of r_qm_in and r_qm_out
    radius = radius + (abfqmmm_param%r_qm_in + abfqmmm_param%r_qm_out) / 2.0d0
  end select

  if(abfqmmm_param%fix_atom_list == 0) then

  ! => breakable bond based atom selection
   do n=1,abfqmmm_param%n_subset_buffer
    i = abfqmmm_param%subset_buffer(n)
    ! cycle if it is already extended core or extended qm or user defined buffer hydrogen
    if( (abfqmmm_param%id(i) < 6) .or. (abfqmmm_param%mass(i) < abfqmmm_param%min_heavy_mass) ) cycle
    ! select selection type
    select case(abfqmmm_param%selection_type)
     ! atom-atom distance selection
     case(1)
      ! calculate the shortest distance between atom and extended core + extended qm list
      mindist=abfqmmm_atom_group_min_distance(i, abfqmmm_param%n_qm, abfqmmm_param%qm, crd)
     ! flexible sphere selection
     case(2)
      ! calculate distance between atom and center of center list
      mindist=abfqmmm_point_atom_dist(center, i, crd)
     ! fixed sphere com selection
     case(3)
      ! calculate distance between atom and center of center list
      mindist=abfqmmm_point_atom_dist(center, i, crd)
    end select
    ! label the atom according to distance critera
    if( (abfqmmm_param%read_idrst_file == '') .and. (abfqmmm_param%qmstep == 1) ) then
     ! if first step and not abfqmmm restart calculation
     select case(abfqmmm_param%initial_selection_type)
      ! inner sphere selection
      case(-1)
       call abfqmmm_atom_labelling_by_cutoff(i, 6, mindist, abfqmmm_param%r_buffer_in+radius, abfqmmm_param%r_buffer_in+radius)
      ! middle sphere selection
      case(0)
       call abfqmmm_atom_labelling_by_cutoff(i, 6, mindist, (abfqmmm_param%r_buffer_in+abfqmmm_param%r_buffer_out)/2.0d0+radius, &
                                                            (abfqmmm_param%r_buffer_in+abfqmmm_param%r_buffer_out)/2.0d0+radius)
      ! outer sphere selection
      case(1)
       call abfqmmm_atom_labelling_by_cutoff(i, 6, mindist, abfqmmm_param%r_buffer_out+radius, abfqmmm_param%r_buffer_out+radius)
     end select
    else
      call abfqmmm_atom_labelling_by_cutoff(i, 6, mindist, abfqmmm_param%r_buffer_in+radius, abfqmmm_param%r_buffer_out+radius)
    end if
   end do

   ! if labeled as buffer, recursively label its environment
   do i=1,abfqmmm_param%natom
    if(abfqmmm_param%id(i) /= 6) cycle
    ! cycle if it is hydrogen
    if(abfqmmm_param%mass(i) < abfqmmm_param%min_heavy_mass) cycle
    abfqmmm_param%n_recursive = 0
    call abfqmmm_recursive_labelling(i, 6)
   end do

   ! label hydrogen atoms 
   do i=1,abfqmmm_param%natom
    ! cycle if it is not hydrogen
    if(abfqmmm_param%mass(i) >= abfqmmm_param%min_heavy_mass) cycle
    abfqmmm_param%n_recursive = 0
    call abfqmmm_recursive_find_label_heavy_atom(i, j)
    abfqmmm_param%id(i) = j
   end do

  else

   ! => fixed atom list based selection 
   do n=1,abfqmmm_param%n_subset_buffer
    i = abfqmmm_param%subset_buffer(n)
    ! select selection type
    select case(abfqmmm_param%selection_type)
     ! atom-atom distance selection
     case(1)
      ! calculate the shortest distance between res and extended core + extended qm list
      mindist=abfqmmm_residue_group_min_distance(i, 6, abfqmmm_param%n_qm, abfqmmm_param%qm, crd)
     ! flexible sphere selection
     case(2)
      ! calculate distance between res and center of center list
      mindist=abfqmmm_point_res_dist(center, i, 6, crd)
     ! fixed sphere selection
     case(3)
      ! calculate distance between res and center of center list
      mindist=abfqmmm_point_res_dist(center, i, 6, crd)
    end select
    ! cycle if already extended core or extended qm or user defined buffer labeled residue
    if(mindist < 0) cycle
    ! label the whole residue according to distance criteria
    if( (abfqmmm_param%read_idrst_file == '') .and. (abfqmmm_param%qmstep == 1) ) then
     ! if first step and not abfqmmm restart calculation
     select case(abfqmmm_param%initial_selection_type)
      ! inner sphere selection
      case(-1)
       call abfqmmm_res_labelling_by_cutoff(i, 6, mindist, abfqmmm_param%r_buffer_in+radius, abfqmmm_param%r_buffer_in+radius)
      ! middle sphere selection
      case(0)
       call abfqmmm_res_labelling_by_cutoff(i, 6, mindist, (abfqmmm_param%r_buffer_in+abfqmmm_param%r_buffer_out)/2.0d0+radius, &
                                                           (abfqmmm_param%r_buffer_in+abfqmmm_param%r_buffer_out)/2.0d0+radius)
      ! outer sphere selection
      case(1)
       call abfqmmm_res_labelling_by_cutoff(i, 6, mindist, abfqmmm_param%r_buffer_out+radius, abfqmmm_param%r_buffer_out+radius)
     end select
    else
      call abfqmmm_res_labelling_by_cutoff(i, 6, mindist, abfqmmm_param%r_buffer_in+radius, abfqmmm_param%r_buffer_out+radius)
    end if
   end do

  end if

  ! Allocate buffer array (buffer extended)

  if( associated(abfqmmm_param%buffer) ) then
    deallocate(abfqmmm_param%buffer, stat=ierr)
    REQUIRE (ierr == 0)
  end if
  abfqmmm_param%n_buffer=0
  do i=1,abfqmmm_param%natom
   if( (abfqmmm_param%id(i) > 4) .and. (abfqmmm_param%id(i) < 7) ) abfqmmm_param%n_buffer=abfqmmm_param%n_buffer+1
  end do
  if(abfqmmm_param%n_buffer>0) then
   allocate(abfqmmm_param%buffer(abfqmmm_param%n_buffer))
   j=1
   do i=1,abfqmmm_param%natom
    if(abfqmmm_param%id(i)>4 .and. abfqmmm_param%id(i)<7) then
     abfqmmm_param%buffer(j)=i
     j=j+1
    end if
   end do
  end if

  if( abfqmmm_param%hot_spot == 1 ) then
    abfqmmm_param%r_hot_spot_in = radius
    abfqmmm_param%r_hot_spot_out = abfqmmm_param%r_buffer_out+radius 
    if( associated(abfqmmm_param%r_hot_spot) ) then
      deallocate(abfqmmm_param%r_hot_spot, stat=ierr)
      REQUIRE (ierr == 0)
    end if
    allocate(abfqmmm_param%r_hot_spot(abfqmmm_param%nres))
    do j=1,abfqmmm_param%nres
       k=abfqmmm_param%res_pointers(j)
       l=abfqmmm_param%res_pointers(j)+abfqmmm_param%res_atom_number(j)-1
       call abfqmmm_center_and_radius_of_atom_list(abfqmmm_param%center_type, &
               abfqmmm_param%res_atom_number(j), (/ (i, i=k, l) /), &
               abfqmmm_param%res_atom_number(j), (/ (i, i=k, l) /), crd, c, r)
       abfqmmm_param%r_hot_spot(j) = abfqmmm_point_point_dist(center, c) 
    end do
  end if

end subroutine abfqmmm_update_qmatoms


function abfqmmm_atom_group_min_distance(atom,num_group,list_group,crd) result(mindist)

  implicit none

  integer, intent(in) :: atom
  integer, intent(in) :: num_group
  integer, intent(in) :: list_group(*)
  _REAL_, intent(in) :: crd(*)

  _REAL_ :: dist, mindist

  integer :: j, k
  logical :: first

  first=.true.
  mindist = 0.d0
  do j=1,num_group
   k=list_group(j)
   ! if hydrogen then cycle
   if(abfqmmm_param%mass(k) < abfqmmm_param%min_heavy_mass) cycle
   dist=abfqmmm_atom_atom_dist(atom,k,crd)
   if(first) then
    mindist=dist
    first=.false.
   end if
   if(mindist>dist) mindist=dist
  end do

end function abfqmmm_atom_group_min_distance

function abfqmmm_residue_group_min_distance(res,label,num_group,list_group,crd) result(mindist)

  implicit none

  integer, intent(in) :: res
  integer, intent(in) :: label
  integer, intent(in) :: num_group
  integer, intent(in) :: list_group(*)
  _REAL_, intent(in) :: crd(*)

  _REAL_ :: dist, mindist

  integer :: i, j, k
  integer :: resstart, resstop
  logical :: first

  resstart=abfqmmm_param%res_pointers(res)
  resstop=abfqmmm_param%res_pointers(res+1)-1
  first=.true.
  mindist = 0.d0

  if(abfqmmm_param%res_atom_number(res) > abfqmmm_param%solvent_atom_number .or. abfqmmm_param%id(resstart) < label) then
   mindist = -1.d0
   return 
  end if

  do i=resstart, resstop
   ! if hydrogen then cycle
   if(abfqmmm_param%mass(i) < abfqmmm_param%min_heavy_mass) cycle
   do j=1,num_group
    k=list_group(j)
    ! if hydrogen then cycle
    if(abfqmmm_param%mass(k) < abfqmmm_param%min_heavy_mass) cycle
    dist=abfqmmm_atom_atom_dist(i,k,crd)
    if(first) then
     mindist=dist
     first=.false.
    endif
    if(mindist>dist) mindist=dist
   end do
  enddo

end function abfqmmm_residue_group_min_distance


subroutine abfqmmm_atom_labelling_by_cutoff(atom, label, dist, min_cutoff, max_cutoff)

  implicit none

  integer, intent(in) :: atom
  integer, intent(in) :: label
  _REAL_, intent(in) :: dist
  _REAL_, intent(in) :: min_cutoff
  _REAL_, intent(in) :: max_cutoff

  if(dist<=min_cutoff) then
   abfqmmm_param%id(atom) = label
  else if(dist>max_cutoff .and. abfqmmm_param%id(atom) == label) then
   abfqmmm_param%id(atom) = abfqmmm_param%id_orig(atom)
  end if

end subroutine abfqmmm_atom_labelling_by_cutoff


subroutine abfqmmm_res_labelling_by_cutoff(res, label, dist, min_cutoff, max_cutoff)

  implicit none

  integer, intent(in) :: res
  integer, intent(in) :: label
  _REAL_, intent(in) :: dist
  _REAL_, intent(in) :: min_cutoff
  _REAL_, intent(in) :: max_cutoff

  integer :: resstart, resstop

  integer :: i

  resstart=abfqmmm_param%res_pointers(res)
  resstop=abfqmmm_param%res_pointers(res+1)-1

  if(dist<=min_cutoff) then
   do i=resstart, resstop
    abfqmmm_param%id(i) = label
   end do
  else if(dist>max_cutoff .and. abfqmmm_param%id(resstart) == label) then
   do i=resstart, resstop
    abfqmmm_param%id(i) = abfqmmm_param%id_orig(i)
   end do
  end if

end subroutine abfqmmm_res_labelling_by_cutoff


function abfqmmm_atom_atom_dist(atom1, atom2, crd) result(dist)

  implicit none

  integer, intent(in) :: atom1
  integer, intent(in) :: atom2
  _REAL_, intent(in) :: crd(*)

  integer :: i, j, k
  _REAL_ :: d(3)
  _REAL_ :: dist

  j = (atom1-1)*3
  k = (atom2-1)*3

  do i=1,3
   d(i)=crd(j+i)-crd(k+i)
  end do

  dist=sqrt(d(1)*d(1)+d(2)*d(2)+d(3)*d(3))

end function abfqmmm_atom_atom_dist


function abfqmmm_point_point_dist(point1, point2) result(dist)

  implicit none

  _REAL_, intent(in) :: point1(3)
  _REAL_, intent(in) :: point2(3)

  integer :: i
  _REAL_ :: d(3)
  _REAL_ :: dist

  do i=1,3
   d(i)=point1(i)-point2(i)
  end do

  dist=sqrt(d(1)*d(1)+d(2)*d(2)+d(3)*d(3))

end function abfqmmm_point_point_dist


function abfqmmm_point_atom_dist(point, atom, crd) result(dist)

  implicit none

  _REAL_, intent(in) :: point(3)
  integer, intent(in) :: atom
  _REAL_, intent(in) :: crd(*)

  integer :: i, j
  _REAL_ :: d(3)
  _REAL_ :: dist

  j = (atom-1)*3

  do i=1,3
   d(i)=point(i)-crd(j+i)
  end do

  dist=sqrt(d(1)*d(1)+d(2)*d(2)+d(3)*d(3))

end function abfqmmm_point_atom_dist


function abfqmmm_point_res_dist(point, res, label, crd) result(mindist)

  implicit none

  _REAL_, intent(in) :: point(3)
  integer, intent(in) :: res
  integer, intent(in) :: label
  _REAL_, intent(in) :: crd(*)

  integer :: i
  _REAL_ :: dist
  _REAL_ :: mindist

  integer :: resstart, resstop
  logical :: first

  resstart=abfqmmm_param%res_pointers(res)
  resstop=abfqmmm_param%res_pointers(res+1)-1
  first=.true.

  if(abfqmmm_param%res_atom_number(res) > abfqmmm_param%solvent_atom_number .or. abfqmmm_param%id(resstart) < label) then
   mindist = -1.0d0
   return
  end if

  mindist = 0.d0
  do i=resstart, resstop
   ! if hydrogen then cycle
   if(abfqmmm_param%mass(i) < abfqmmm_param%min_heavy_mass) cycle
   dist=abfqmmm_point_atom_dist(point,i,crd)
   if(first) then
    mindist=dist
    first=.false.
   endif
   if(mindist>dist) mindist=dist
  enddo

end function abfqmmm_point_res_dist


subroutine abfqmmm_center_and_radius_of_atom_list(selection_type,c_num_group,c_list_group,num_group,list_group,crd,center,radius)

  implicit none

  integer, intent(in) :: selection_type
  integer, intent(in) :: c_num_group
  integer, intent(in) :: c_list_group(*)
  integer, intent(in) :: num_group
  integer, intent(in) :: list_group(*)
  _REAL_, intent(in) :: crd(*)
  _REAL_, intent(out) :: center(3)
  _REAL_, intent(out) :: radius
  
  integer :: i, j, k, l
  _REAL_ :: totalmass
  _REAL_ :: dist

  center(:) = 0.0d0
  totalmass = 0.0d0
  radius = 0.0d0

  ! calculate com or geometric center of center group
  do j=1,c_num_group
   k=c_list_group(j)
   ! if hydrogen then cycle
   if(abfqmmm_param%mass(k) < abfqmmm_param%min_heavy_mass) cycle
   i = (k-1)*3
   select case(selection_type)
    ! com selection
    case(1)
     do l=1,3
      center(l) = center(l) + crd(i+l)*abfqmmm_param%mass(k)
     end do
     totalmass = totalmass + abfqmmm_param%mass(k)
    ! geometric center selection
    case(2)
     do l=1,3
      center(l) = center(l) + crd(i+l)
     end do
     totalmass = totalmass + 1.0d0
   end select
  end do

  center(:) = center(:) / totalmass

  ! calculate radius of sphere
  do j=1,num_group
   k=list_group(j)
   ! if hydrogen then cycle
   if(abfqmmm_param%mass(k) < abfqmmm_param%min_heavy_mass) cycle
   dist=abfqmmm_point_atom_dist(center,k,crd)
   if(dist > radius) radius = dist
  end do

end subroutine abfqmmm_center_and_radius_of_atom_list


subroutine abfqmmm_select_system_qmatoms(natom)

  use qmmm_module, only : qmmm_nml

  implicit none

  integer, intent(in) :: natom
  integer :: i
  integer :: charge

  if(abfqmmm_param%system == 1) then
   charge = 0
   do i=1,natom
    if(abfqmmm_param%id(i) <= 6) then
     abfqmmm_param%isqm(i)=1
     ! for not user defined atoms calculate the extra charges based on oxidation numbers
     if(abfqmmm_param%fix_atom_list == 0) then
      if(abfqmmm_param%id_orig(i) /= 1 .and. abfqmmm_param%id_orig(i) /= 3 .and. abfqmmm_param%id_orig(i) /= 5) then
       charge = charge + abfqmmm_param%oxidation_number(i)
      end if
     end if
    else
     abfqmmm_param%isqm(i)=0
    end if
   end do
   abfqmmm_param%abfcharge = abfqmmm_param%corecharge + abfqmmm_param%qmcharge + abfqmmm_param%buffercharge + charge
   qmmm_nml%ifqnt = .true.
  else if(abfqmmm_param%system == 2) then
   ! do full MM if there are no qm and buffer atoms AND an unbuffered adaptive QM/MM is performed (i.e. only core radii are defined) 
   if( (abfqmmm_param%n_user_core == abfqmmm_param%n_user_qm) .and. (abfqmmm_param%n_user_buffer == 0) .and. &
       (abfqmmm_param%r_qm_out == 0.0d0) .and. (abfqmmm_param%r_buffer_out == 0.0d0) ) then
    do i=1,natom
     abfqmmm_param%isqm(i)=0
    enddo
    qmmm_nml%ifqnt = .false.
   ! do full MM if there are no core atoms
   else if (abfqmmm_param%n_user_core == 0) then
    do i=1,natom
     abfqmmm_param%isqm(i)=0
    enddo
    qmmm_nml%ifqnt = .false.
   else
    charge = 0
    do i=1,natom
     if(abfqmmm_param%id(i) <= 2) then
      abfqmmm_param%isqm(i)=1
      ! for not user defined atoms calculate the extra charges based on oxidation numbers
      if(abfqmmm_param%fix_atom_list == 0) then
       if(abfqmmm_param%id_orig(i) /= 1) then
        charge = charge + abfqmmm_param%oxidation_number(i)
       end if
      end if
     else
      abfqmmm_param%isqm(i)=0
     end if
    end do
    abfqmmm_param%abfcharge = abfqmmm_param%corecharge + charge
    qmmm_nml%ifqnt = .true.
   end if
  end if

end subroutine abfqmmm_select_system_qmatoms


subroutine abfqmmm_allocate_arrays_of_parameters(numbnd, nbonh, nbona, ntheth, ntheta, nphih, nphia)

  implicit none

  integer, intent(in) :: numbnd
  integer, intent(in) :: nbonh
  integer, intent(in) :: nbona
  integer, intent(in) :: ntheth
  integer, intent(in) :: ntheta
  integer, intent(in) :: nphih
  integer, intent(in) :: nphia

  integer :: ierr

  abfqmmm_param%numbnd = numbnd
  abfqmmm_param%nbonh  = nbonh
  abfqmmm_param%nbona  = nbona
  abfqmmm_param%ntheth = ntheth
  abfqmmm_param%ntheta = ntheta
  abfqmmm_param%nphih  = nphih
  abfqmmm_param%nphia  = nphia

    ! mm charges
    if ( .not. associated(abfqmmm_param%charge) ) then
       allocate (abfqmmm_param%charge(abfqmmm_param%natom), stat=ierr)
       REQUIRE(ierr == 0)
    end if

    ! general bond parameters
    if ( .not. associated(abfqmmm_param%rk) ) then
       allocate (abfqmmm_param%rk(numbnd), stat=ierr)
       REQUIRE(ierr == 0)
    end if 

    if ( .not. associated(abfqmmm_param%req) ) then
       allocate (abfqmmm_param%req(numbnd), stat=ierr)
       REQUIRE(ierr == 0)
    end if

    ! bonds with Hydrogen
    if ( .not. associated(abfqmmm_param%iibh) ) then
       allocate (abfqmmm_param%iibh(nbonh), stat=ierr)
       REQUIRE(ierr == 0)
    end if

    if ( .not. associated(abfqmmm_param%ijbh) ) then
       allocate (abfqmmm_param%ijbh(nbonh), stat=ierr)
       REQUIRE(ierr == 0)
    end if

    if ( .not. associated(abfqmmm_param%icbh) ) then
       allocate (abfqmmm_param%icbh(nbonh), stat=ierr)
       REQUIRE(ierr == 0)
    end if

    ! bonds without hydrogen atoms
    if ( .not. associated(abfqmmm_param%iiba) ) then
       allocate (abfqmmm_param%iiba(nbona), stat=ierr)
       REQUIRE(ierr == 0)
    end if

    if ( .not. associated(abfqmmm_param%ijba) ) then
       allocate (abfqmmm_param%ijba(nbona), stat=ierr)
       REQUIRE(ierr == 0)
    end if

    if ( .not. associated(abfqmmm_param%icba) ) then
       allocate (abfqmmm_param%icba(nbona), stat=ierr)
       REQUIRE(ierr == 0)
    end if

    ! angles with hydrogen atoms
    if ( .not. associated(abfqmmm_param%iith) ) then
       allocate (abfqmmm_param%iith(ntheth), stat=ierr)
       REQUIRE(ierr == 0)
    end if

    if ( .not. associated(abfqmmm_param%ijth) ) then
       allocate (abfqmmm_param%ijth(ntheth), stat=ierr)
       REQUIRE(ierr == 0)
    end if

    if ( .not. associated(abfqmmm_param%ikth) ) then
       allocate (abfqmmm_param%ikth(ntheth), stat=ierr)
       REQUIRE(ierr == 0)
    end if

    if ( .not. associated(abfqmmm_param%icth) ) then
       allocate (abfqmmm_param%icth(ntheth), stat=ierr)
       REQUIRE(ierr == 0)
    end if

    ! angles without hydrogen atoms
    if ( .not. associated(abfqmmm_param%iita) ) then
       allocate (abfqmmm_param%iita(ntheta), stat=ierr)
       REQUIRE(ierr == 0)
    end if

    if ( .not. associated(abfqmmm_param%ijta) ) then
       allocate (abfqmmm_param%ijta(ntheta), stat=ierr)
       REQUIRE(ierr == 0)
    end if

    if ( .not. associated(abfqmmm_param%ikta) ) then
       allocate (abfqmmm_param%ikta(ntheta), stat=ierr)
       REQUIRE(ierr == 0)
    end if

    if ( .not. associated(abfqmmm_param%icta) ) then
       allocate (abfqmmm_param%icta(ntheta), stat=ierr)
       REQUIRE(ierr == 0)
    end if

    ! dihedrals with hydrogen atoms
    if ( .not. associated(abfqmmm_param%iiph) ) then
       allocate (abfqmmm_param%iiph(nphih), stat=ierr)
       REQUIRE(ierr == 0)
    end if

    if ( .not. associated(abfqmmm_param%ijph) ) then
       allocate (abfqmmm_param%ijph(nphih), stat=ierr)
       REQUIRE(ierr == 0)
    end if

    if ( .not. associated(abfqmmm_param%ikph) ) then
       allocate (abfqmmm_param%ikph(nphih), stat=ierr)
       REQUIRE(ierr == 0)
    end if

    if ( .not. associated(abfqmmm_param%ilph) ) then
       allocate (abfqmmm_param%ilph(nphih), stat=ierr)
       REQUIRE(ierr == 0)
    end if

    if ( .not. associated(abfqmmm_param%icph) ) then
       allocate (abfqmmm_param%icph(nphih), stat=ierr)
       REQUIRE(ierr == 0)
    end if

    ! dihedrals without hydrogen atoms
    if ( .not. associated(abfqmmm_param%iipa) ) then
       allocate (abfqmmm_param%iipa(nphia), stat=ierr)
       REQUIRE(ierr == 0)
    end if

    if ( .not. associated(abfqmmm_param%ijpa) ) then
       allocate (abfqmmm_param%ijpa(nphia), stat=ierr)
       REQUIRE(ierr == 0)
    end if

    if ( .not. associated(abfqmmm_param%ikpa) ) then
       allocate (abfqmmm_param%ikpa(nphia), stat=ierr)
       REQUIRE(ierr == 0)
    end if

    if ( .not. associated(abfqmmm_param%ilpa) ) then
       allocate (abfqmmm_param%ilpa(nphia), stat=ierr)
       REQUIRE(ierr == 0)
    end if

    if ( .not. associated(abfqmmm_param%icpa) ) then
       allocate (abfqmmm_param%icpa(nphia), stat=ierr)
       REQUIRE(ierr == 0)
    end if

    ! charmm
    if(charmm_active) then

       ! Urey-Bradley
       abfqmmm_param%charmm_nub = charmm_nub
       if ( .not. associated(abfqmmm_param%charmm_ang_ub) ) then
          allocate (abfqmmm_param%charmm_ang_ub(charmm_nub), stat=ierr)
          REQUIRE(ierr == 0)
       end if

       ! improper
       abfqmmm_param%charmm_nimphi = charmm_nimphi
       if ( .not. associated(abfqmmm_param%charmm_imp) ) then
          allocate (abfqmmm_param%charmm_imp(charmm_nimphi), stat=ierr)
          REQUIRE(ierr == 0)
       end if

       ! cmap
       abfqmmm_param%cmap_term_count = cmap_term_count
       if ( .not. associated(abfqmmm_param%cmap_index) ) then
          allocate (abfqmmm_param%cmap_index(6,cmap_term_count), stat=ierr)
          REQUIRE(ierr == 0)
       end if

    end if
  
end subroutine abfqmmm_allocate_arrays_of_parameters 


subroutine abfqmmm_store_parameters(iibh, ijbh, icbh, &
                                    iiba, ijba, icba, &
                                    iith, ijth, ikth, icth, &
                                    iita, ijta, ikta, icta, &
                                    iiph, ijph, ikph, ilph, icph, &
                                    iipa, ijpa, ikpa, ilpa, icpa, &
                                    charge, rk, req)

  implicit none

    integer, intent(in) :: iibh(abfqmmm_param%nbonh), ijbh(abfqmmm_param%nbonh), icbh(abfqmmm_param%nbonh)
    integer, intent(in) :: iiba(abfqmmm_param%nbona), ijba(abfqmmm_param%nbona), icba(abfqmmm_param%nbona)
    integer, intent(in) :: iith(abfqmmm_param%ntheth), ijth(abfqmmm_param%ntheth), &
                           ikth(abfqmmm_param%ntheth), icth(abfqmmm_param%ntheth)
    integer, intent(in) :: iita(abfqmmm_param%ntheta), ijta(abfqmmm_param%ntheta), &
                           ikta(abfqmmm_param%ntheta), icta(abfqmmm_param%ntheta)
    integer, intent(in) :: iiph(abfqmmm_param%nphih), ijph(abfqmmm_param%nphih), ikph(abfqmmm_param%nphih), &
                           ilph(abfqmmm_param%nphih), icph(abfqmmm_param%nphih)
    integer, intent(in) :: iipa(abfqmmm_param%nphia), ijpa(abfqmmm_param%nphia), ikpa(abfqmmm_param%nphia), &
                           ilpa(abfqmmm_param%nphia), icpa(abfqmmm_param%nphia)

    _REAL_, intent(in) :: charge(abfqmmm_param%natom)

    _REAL_, intent(in) :: rk(abfqmmm_param%numbnd)
    _REAL_, intent(in) :: req(abfqmmm_param%numbnd)

    abfqmmm_param%iibh = iibh
    abfqmmm_param%ijbh = ijbh
    abfqmmm_param%icbh = icbh

    abfqmmm_param%iiba = iiba
    abfqmmm_param%ijba = ijba
    abfqmmm_param%icba = icba

    abfqmmm_param%iith = iith
    abfqmmm_param%ijth = ijth
    abfqmmm_param%ikth = ikth
    abfqmmm_param%icth = icth

    abfqmmm_param%iita = iita
    abfqmmm_param%ijta = ijta
    abfqmmm_param%ikta = ikta
    abfqmmm_param%icta = icta

    abfqmmm_param%iiph = iiph
    abfqmmm_param%ijph = ijph
    abfqmmm_param%ikph = ikph
    abfqmmm_param%ilph = ilph
    abfqmmm_param%icph = icph

    abfqmmm_param%iipa = iipa
    abfqmmm_param%ijpa = ijpa
    abfqmmm_param%ikpa = ikpa
    abfqmmm_param%ilpa = ilpa
    abfqmmm_param%icpa = icpa

    abfqmmm_param%charge = charge

    abfqmmm_param%rk = rk
    abfqmmm_param%req = req

    if(charmm_active) then

      abfqmmm_param%charmm_ang_ub = charmm_ang_ub
      abfqmmm_param%charmm_imp = charmm_imp
      abfqmmm_param%cmap_index = cmap_index

    end if

end subroutine abfqmmm_store_parameters


subroutine abfqmmm_set_parameters(numbnd, nbonh, nbona, ntheth, ntheta, nphih, nphia, &
                                  iibh, ijbh, icbh, &
                                  iiba, ijba, icba, &
                                  iith, ijth, ikth, icth, &
                                  iita, ijta, ikta, icta, &
                                  iiph, ijph, ikph, ilph, icph, &
                                  iipa, ijpa, ikpa, ilpa, icpa, &
                                  charge, rk, req)


  implicit none

  integer, intent(inout) :: numbnd
  integer, intent(inout) :: nbonh
  integer, intent(inout) :: nbona
  integer, intent(inout) :: ntheth
  integer, intent(inout) :: ntheta
  integer, intent(inout) :: nphih
  integer, intent(inout) :: nphia

  integer, intent(inout) :: iibh(abfqmmm_param%nbonh), ijbh(abfqmmm_param%nbonh), icbh(abfqmmm_param%nbonh)
  integer, intent(inout) :: iiba(abfqmmm_param%nbona), ijba(abfqmmm_param%nbona), icba(abfqmmm_param%nbona)
  integer, intent(inout) :: iith(abfqmmm_param%ntheth), ijth(abfqmmm_param%ntheth), &
                            ikth(abfqmmm_param%ntheth), icth(abfqmmm_param%ntheth)
  integer, intent(inout) :: iita(abfqmmm_param%ntheta), ijta(abfqmmm_param%ntheta), &
                            ikta(abfqmmm_param%ntheta), icta(abfqmmm_param%ntheta)
  integer, intent(inout) :: iiph(abfqmmm_param%nphih), ijph(abfqmmm_param%nphih), ikph(abfqmmm_param%nphih), &
                            ilph(abfqmmm_param%nphih), icph(abfqmmm_param%nphih)
  integer, intent(inout) :: iipa(abfqmmm_param%nphia), ijpa(abfqmmm_param%nphia), ikpa(abfqmmm_param%nphia), &
                            ilpa(abfqmmm_param%nphia), icpa(abfqmmm_param%nphia)

  _REAL_, intent(inout) :: charge(abfqmmm_param%natom)

  _REAL_, intent(inout) :: rk(abfqmmm_param%numbnd)
  _REAL_, intent(inout) :: req(abfqmmm_param%numbnd)

    numbnd = abfqmmm_param%numbnd

    nbonh  = abfqmmm_param%nbonh    
    nbona  = abfqmmm_param%nbona  
    ntheth = abfqmmm_param%ntheth 
    ntheta = abfqmmm_param%ntheta 
    nphih  = abfqmmm_param%nphih  
    nphia  = abfqmmm_param%nphia  

    iibh  = abfqmmm_param%iibh
    ijbh  = abfqmmm_param%ijbh
    icbh  = abfqmmm_param%icbh
            
    iiba  = abfqmmm_param%iiba
    ijba  = abfqmmm_param%ijba
    icba  = abfqmmm_param%icba
            
    iith  = abfqmmm_param%iith
    ijth  = abfqmmm_param%ijth
    ikth  = abfqmmm_param%ikth
    icth  = abfqmmm_param%icth
            
    iita  = abfqmmm_param%iita
    ijta  = abfqmmm_param%ijta
    ikta  = abfqmmm_param%ikta
    icta  = abfqmmm_param%icta
            
    iiph  = abfqmmm_param%iiph
    ijph  = abfqmmm_param%ijph
    ikph  = abfqmmm_param%ikph
    ilph  = abfqmmm_param%ilph
    icph  = abfqmmm_param%icph
            
    iipa  = abfqmmm_param%iipa
    ijpa  = abfqmmm_param%ijpa
    ikpa  = abfqmmm_param%ikpa
    ilpa  = abfqmmm_param%ilpa
    icpa  = abfqmmm_param%icpa

    charge = abfqmmm_param%charge

    rk = abfqmmm_param%rk
    req = abfqmmm_param%req

    if(charmm_active) then

       charmm_nub = abfqmmm_param%charmm_nub
       if(allocated(charmm_ang_ub)) deallocate(charmm_ang_ub)
       allocate(charmm_ang_ub(charmm_nub))
       charmm_ang_ub = abfqmmm_param%charmm_ang_ub

       charmm_nimphi = abfqmmm_param%charmm_nimphi
       if(allocated(charmm_imp)) deallocate(charmm_imp)
       allocate(charmm_imp(charmm_nimphi))
       charmm_imp = abfqmmm_param%charmm_imp

       cmap_term_count = abfqmmm_param%cmap_term_count
       if(allocated(cmap_index)) deallocate(cmap_index)
       allocate(cmap_index(6,cmap_term_count))
       cmap_index = abfqmmm_param%cmap_index

    end if

end subroutine abfqmmm_set_parameters


subroutine abfqmmm_vel_verlet1(x, v, invmass)

  implicit none

#include "../include/md.h"

  _REAL_, intent(inout) :: x(*)

  _REAL_, intent(inout) :: v(*)

  _REAL_, intent(in) :: invmass(*)

  integer i, j

  _REAL_ dtx, halfdt, halfdtperm

  ! this happens before force calculation: calculate new positions [x(t+dt)] and velocities at half step [v(t+1/2dt)]
  ! if first loop then skip
  ! combined forces [abfqmmm_param%f] are already calculated after the second system in abfqmmm_vel_verlet2 

  if(abfqmmm_param%qmstep == 1) return

  ! calculate dt and 1/2dt in unit compatible with kcal/mol

  dtx=dt*20.455d+00
  halfdt=0.5*dtx

  ! calculate v(t+1/2dt) using v(t) and f(t): v(t+1/2dt) = v(t) + 1/2dt*f/m

  do i=1,abfqmmm_param%natom
   j=3*(i-1)
     
   halfdtperm=halfdt*invmass(i)

   v(j+1)=v(j+1)+halfdtperm*abfqmmm_param%f(j+1)
   v(j+2)=v(j+2)+halfdtperm*abfqmmm_param%f(j+2)
   v(j+3)=v(j+3)+halfdtperm*abfqmmm_param%f(j+3)
  end do

  ! calculate x(t+dt) using v(t+1/2dt): x(t+dt) = x(t) + dt*v(t+1/2dt)

  do i=1,abfqmmm_param%natom
   j=3*(i-1)
    
   x(j+1)=x(j+1)+dtx*v(j+1)
   x(j+2)=x(j+2)+dtx*v(j+2)
   x(j+3)=x(j+3)+dtx*v(j+3)
  end do

end subroutine abfqmmm_vel_verlet1


subroutine abfqmmm_combine_forces()

  implicit none

  integer   :: i, j, k, id
  _REAL_    :: force(3), absforce(3), absacc(3)
  _REAL_    :: mass, num
  _REAL_    :: smoothing

  force(:) = 0
  absforce(:) = 0
  absacc(:) = 0
  num = 0
  mass = 0

  ! if pure adaptive QM/MM is performed then keep forces only from the QM/MM ('extended') calculation
  if(abfqmmm_param%r_qm_out == 0.0d0 .and. abfqmmm_param%r_buffer_out == 0.0d0) then
   abfqmmm_param%f(:) = abfqmmm_param%f1(:)
   return
  end if

  ! chose threshold id for momentum conservation
  id = 0
  select case(abfqmmm_param%mom_cons_region)
   ! core region
   case(0)
    id = 2
   ! core+qm region
   case(1)
    id = 4
   ! core+qm+buffer region
   case(2)
    id = 6
   ! all atom
   case(3)
    id = 7
   case default
    id = 0
    REQUIRE( .false. )
  end select

  do i=1,abfqmmm_param%natom
   j=3*(i-1)

   if(abfqmmm_param%id(i)<=4) then
    abfqmmm_param%f(j+1) = abfqmmm_param%f1(j+1)
    abfqmmm_param%f(j+2) = abfqmmm_param%f1(j+2)
    abfqmmm_param%f(j+3) = abfqmmm_param%f1(j+3)
   else
    abfqmmm_param%f(j+1) = abfqmmm_param%f2(j+1)
    abfqmmm_param%f(j+2) = abfqmmm_param%f2(j+2)
    abfqmmm_param%f(j+3) = abfqmmm_param%f2(j+3)
   end if

   if(abfqmmm_param%hot_spot == 1) then
     if( (abfqmmm_param%id(i)==5) .or. (abfqmmm_param%id(i)==6) ) then
        k = abfqmmm_param%res_id_of_atom(i)
        if(abfqmmm_param%r_hot_spot(k) > abfqmmm_param%r_hot_spot_out) then
         smoothing = 0.0d0
        else if(abfqmmm_param%r_hot_spot(k) <= abfqmmm_param%r_hot_spot_in) then
         smoothing = 1.0d0
        else
         smoothing = (abfqmmm_param%r_hot_spot_out**2-abfqmmm_param%r_hot_spot(k)**2)**2 &
                   * (abfqmmm_param%r_hot_spot_out**2+2.0d0*abfqmmm_param%r_hot_spot(k)**2-3.0d0*abfqmmm_param%r_hot_spot_in**2) &
                   / (abfqmmm_param%r_hot_spot_out**2-abfqmmm_param%r_hot_spot_in**2)**3
        end if
        abfqmmm_param%f(j+1) = smoothing * abfqmmm_param%f1(j+1) + (1.0d0 - smoothing) * abfqmmm_param%f2(j+1)
        abfqmmm_param%f(j+2) = smoothing * abfqmmm_param%f1(j+2) + (1.0d0 - smoothing) * abfqmmm_param%f2(j+2)
        abfqmmm_param%f(j+3) = smoothing * abfqmmm_param%f1(j+3) + (1.0d0 - smoothing) * abfqmmm_param%f2(j+3)
     end if
   end if

   force(1) = force(1) + abfqmmm_param%f(j+1)
   force(2) = force(2) + abfqmmm_param%f(j+2)
   force(3) = force(3) + abfqmmm_param%f(j+3)

   select case(abfqmmm_param%mom_cons_type)
    ! no momentum conservation
    case(0)
    ! equal acceleration
    case(1)
     if(abfqmmm_param%id(i)<=id) then
      mass = mass + abfqmmm_param%mass(i)
     end if
    ! equal force
    case(2)
     if(abfqmmm_param%id(i)<=id) then
      num = num + 1.0d0
     end if
    ! proportional acceleration
    case(-1)
     if(abfqmmm_param%id(i)<=id) then
      mass = mass + abfqmmm_param%mass(i)
      absacc(1) = absacc(1) + abs(abfqmmm_param%f(j+1)) / abfqmmm_param%mass(i)
      absacc(2) = absacc(2) + abs(abfqmmm_param%f(j+2)) / abfqmmm_param%mass(i)
      absacc(3) = absacc(3) + abs(abfqmmm_param%f(j+3)) / abfqmmm_param%mass(i)
      absforce(1) = absforce(1) + abs(abfqmmm_param%f(j+1))
      absforce(2) = absforce(2) + abs(abfqmmm_param%f(j+2))
      absforce(3) = absforce(3) + abs(abfqmmm_param%f(j+3))
     end if
    ! proportional force
    case(-2)
     if(abfqmmm_param%id(i)<=id) then
      absforce(1) = absforce(1) + abs(abfqmmm_param%f(j+1))
      absforce(2) = absforce(2) + abs(abfqmmm_param%f(j+2))
      absforce(3) = absforce(3) + abs(abfqmmm_param%f(j+3))
     end if
   end select
  end do

  do i=1,abfqmmm_param%natom
   j=3*(i-1)

   select case(abfqmmm_param%mom_cons_type)
    ! no momentum conservation
    case(0)
    ! equal acceleration
    case(1)
     if(abfqmmm_param%id(i)<=id) then
      abfqmmm_param%f(j+1) = abfqmmm_param%f(j+1) - force(1) * abfqmmm_param%mass(i) / mass
      abfqmmm_param%f(j+2) = abfqmmm_param%f(j+2) - force(2) * abfqmmm_param%mass(i) / mass
      abfqmmm_param%f(j+3) = abfqmmm_param%f(j+3) - force(3) * abfqmmm_param%mass(i) / mass
     end if
    ! equal force
    case(2)
     if(abfqmmm_param%id(i)<=id) then
      abfqmmm_param%f(j+1) = abfqmmm_param%f(j+1) - force(1) / num
      abfqmmm_param%f(j+2) = abfqmmm_param%f(j+2) - force(2) / num
      abfqmmm_param%f(j+3) = abfqmmm_param%f(j+3) - force(3) / num
     end if
    ! proportional acceleration
    case(-1)
     if(abfqmmm_param%id(i)<=id) then
      abfqmmm_param%f(j+1) = abfqmmm_param%f(j+1) - force(1) * abs(abfqmmm_param%f(j+1)) / absacc(1) / mass
      abfqmmm_param%f(j+2) = abfqmmm_param%f(j+2) - force(2) * abs(abfqmmm_param%f(j+2)) / absacc(2) / mass
      abfqmmm_param%f(j+3) = abfqmmm_param%f(j+3) - force(3) * abs(abfqmmm_param%f(j+3)) / absacc(3) / mass
     end if
    ! proportional force
    case(-2)
     if(abfqmmm_param%id(i)<=id) then
      abfqmmm_param%f(j+1) = abfqmmm_param%f(j+1) - abs(abfqmmm_param%f(j+1)) * force(1) / absforce(1)
      abfqmmm_param%f(j+2) = abfqmmm_param%f(j+2) - abs(abfqmmm_param%f(j+2)) * force(2) / absforce(2)
      abfqmmm_param%f(j+3) = abfqmmm_param%f(j+3) - abs(abfqmmm_param%f(j+3)) * force(3) / absforce(3)
     end if
   end select
  end do

end subroutine abfqmmm_combine_forces


subroutine abfqmmm_vel_verlet2(v, invmass)

  implicit none

#include "../include/md.h"

  _REAL_, intent(inout) :: v(*)

  _REAL_, intent(in) :: invmass(*)

  integer i, j

  _REAL_ dtx, halfdt, halfdtperm

  ! this happens after the force caclulation: combine the forces and calculate velocities [v(t+dt)]
  ! if first loop than just get the forces for the next step
  ! if qmstep > 1 than do the velocity calculation

  ! combine the forces

  do i=1,abfqmmm_param%natom
   j=3*(i-1)

   if(abfqmmm_param%id(i)<=2) then
    abfqmmm_param%f(j+1)=abfqmmm_param%f1(j+1)
    abfqmmm_param%f(j+2)=abfqmmm_param%f1(j+2)
    abfqmmm_param%f(j+3)=abfqmmm_param%f1(j+3)
   else
    abfqmmm_param%f(j+1)=abfqmmm_param%f2(j+1)
    abfqmmm_param%f(j+2)=abfqmmm_param%f2(j+2)
    abfqmmm_param%f(j+3)=abfqmmm_param%f2(j+3)
   end if
  end do

  if(abfqmmm_param%qmstep == 1) return

  ! calculate dt and 1/2dt in unit compatible with kcal/mol

  dtx=dt*20.455d+00
  halfdt=0.5*dtx

  ! calculate v(t+dt) using v(t+1/2dt) and f(t+dt): v(t+dt) = v(t+1/2dt) + 1/2dt*f/m

  do i=1,abfqmmm_param%natom
   j=3*(i-1)
     
   halfdtperm=halfdt*invmass(i)
 
   v(j+1)=v(j+1)+halfdtperm*abfqmmm_param%f(j+1)
   v(j+2)=v(j+2)+halfdtperm*abfqmmm_param%f(j+2)
   v(j+3)=v(j+3)+halfdtperm*abfqmmm_param%f(j+3)
  end do

end subroutine abfqmmm_vel_verlet2


subroutine abfqmmm_next_step()

  implicit none

#include "../include/md.h"

  if (abfqmmm_param%system == 1) then               
    abfqmmm_param%system = 2                         
  else
    abfqmmm_param%system = 1                         
    abfqmmm_param%qmstep = abfqmmm_param%qmstep + 1
  end if

end subroutine abfqmmm_next_step


subroutine abfqmmm_read_idrst()

  implicit none

  integer :: i, j
  integer :: ios

  write(6,*) ''
  if (abfqmmm_param%read_idrst_file /= '') then
    open(unit=1989,file=abfqmmm_param%read_idrst_file,status='old',iostat=ios)

    if(ios /=0) then
      write(6,*) 'ERROR: cannot open read_idrst_file: ', trim(abfqmmm_param%read_idrst_file)
      stop
    end if

    write(6,'(a)') 'IDRST FILE WAS SPECIFIED - USER QM/CORE/BUFFER ID SPECIFICATIONS IN INPUT FILE ARE IGNORED'
    write(6,'(a)') 'WARNING: PARAMETERS OTHER THAN USED IN PREVIOUS RUN CAN CAUSE SIGNIFICANT PROBLEMS!' 
    write(6,'(a)') 'INITIAL ATOM ASSIGNMENT IS READ FROM IDRST FILE'
  
    do i=1,abfqmmm_param%natom
     read(1989,'(i6,1x,a6,1x,a3,1x,i7,1x,i6,1x,i6)',iostat=ios) j, abfqmmm_param%atom_name(i), &
                                abfqmmm_param%res_name_of_atom(i), abfqmmm_param%res_id_of_atom(i), &
                                abfqmmm_param%id_orig(i), abfqmmm_param%id(i)
     if(ios /=0) then
      write(6,*) 'ERROR: cannot read id info for atom: ', i
      stop
     end if
    end do
    close(1989)
  else
    write(6,'(a)') 'NO IDRST FILE WAS SPECIFIED'
    write(6,'(a)') 'ATOM SELECTION IS BASED ON USER SPECIFICATION'
  end if

end subroutine abfqmmm_read_idrst


subroutine abfqmmm_write_idrst()

  implicit none

  integer :: i
  integer :: ios

  if (abfqmmm_param%ntwidrst == 0) then
    if (abfqmmm_param%qmstep /= abfqmmm_param%maxqmstep) return
  else
    if (mod(abfqmmm_param%qmstep,abfqmmm_param%ntwidrst) /= 0 .and. abfqmmm_param%qmstep /= abfqmmm_param%maxqmstep) return
  end if

  open(unit=1991,file=abfqmmm_param%write_idrst_file,status='unknown',iostat=ios)

  if(ios /=0) then
    write(6,*) 'ERROR: cannot open write_idrst_file: ', trim(abfqmmm_param%write_idrst_file)
    stop
  end if

  do i=1,abfqmmm_param%natom
     write(1991,'(i6,1x,a6,1x,a3,1x,i7,1x,i6,1x,i6)',iostat=ios) i, abfqmmm_param%atom_name(i), &
                                abfqmmm_param%res_name_of_atom(i), abfqmmm_param%res_id_of_atom(i), &
                                abfqmmm_param%id_orig(i), abfqmmm_param%id(i)
  end do

  close(1991)

end subroutine abfqmmm_write_idrst


subroutine abfqmmm_write_pdb(crd,nsp)

use molecule, only: n_iwrap_mask_atoms, iwrap_mask_atoms

  implicit none

#include "../include/md.h"
#include "box.h"

  _REAL_, intent(in) :: crd(3*abfqmmm_param%natom)
  integer, intent(in) :: nsp(*)

  _REAL_ :: tmp_crd(3*abfqmmm_param%natom)
  integer :: i,j,k
  integer :: ios
  character(len=11) :: tag
  _REAL_ :: box_center(3)

  if (abfqmmm_param%ntwpdb == 0) return
  if ((mod(abfqmmm_param%qmstep,abfqmmm_param%ntwpdb) /= 0) .and. (abfqmmm_param%ntwpdb > 0)) return

  if ((abfqmmm_param%qmstep == abfqmmm_param%ntwpdb) .or. (abfqmmm_param%ntwpdb < 0)) then
    open(unit=1987,file=abfqmmm_param%pdb_file,status='unknown',iostat=ios)

    if(ios /=0) then
      write(6,*) 'ERROR: cannot open pdb_file: ', trim(abfqmmm_param%pdb_file)
      stop
    end if
  end if

  tmp_crd(:)=crd(:)

  if(ntb > 0) then
   if(iwrap == 1) then
    call wrap_molecules(nspm,nsp,tmp_crd)
    if (ifbox == 2) call wrap_to(nspm,nsp,tmp_crd,box)
   else if (iwrap == 2) then
    call iwrap2(n_iwrap_mask_atoms,iwrap_mask_atoms,tmp_crd,box_center)
   end if
  end if

  if(ntb > 0) then
    write(1987,'(a6,f9.3,f9.3,f9.3,f7.2,f7.2,f7.2,a2,i4)') 'CRYST1', box(1), box(2), box(3), 90.0, 90.0, 90.0, 'P', 1
  end if  
  do i=1,abfqmmm_param%natom
   j=(i-1)*3
   select case (abfqmmm_param%id(i))
    case (1)
     tag='CORE_USER'
    case (2)
     tag='CORE_EXT'
    case (3)
     tag='QM_USER'
    case (4)
     tag='QM_EXT'
    case (5)
     tag='BUFFER_USER'
    case (6)
     tag='BUFFER_EXT'
    case (7)
     tag='MM'
   end select
   write(1987,advance='NO',fmt='(a5,i6,1x,a4,1x,a3,1x,i5,4x,f8.3,f8.3,f8.3,i6,i5,2x,a11)') 'ATOM ', i, abfqmmm_param%atom_name(i), &
                                            abfqmmm_param%res_name_of_atom(i), abfqmmm_param%res_id_of_atom(i), &
                                            tmp_crd(j+1), tmp_crd(j+2), tmp_crd(j+3), abfqmmm_param%oxidation_number(i), &
                                            abfqmmm_param%id(i), tag
   if(abfqmmm_param%id(i)<7) then
    do k=1,abfqmmm_param%cutnumbond(i)
     if(abfqmmm_param%id(abfqmmm_param%cutlistbond(i,k)) == 7) write(1987,advance='NO',fmt='(1x,i6)') abfqmmm_param%cutlistbond(i,k)
    end do
   end if
   write(1987,*)
  end do 
  write(1987,'(a3)') 'END'

  if ( (abfqmmm_param%ntwpdb > abfqmmm_param%maxqmstep - abfqmmm_param%qmstep) .or. (abfqmmm_param%ntwpdb < 0) ) then
    close(1987)
  end if

  if(abfqmmm_param%ntwpdb < 0) then
   write(6,*) 'INFO: ntwpdb < 0 so this run was used only for generating pdb file!'
   write(6,*) 'INFO: printing is done, job has finished.'
  end if

end subroutine abfqmmm_write_pdb


subroutine abfqmmm_diffusion_restraint(crd)

  implicit none

  _REAL_, intent(in) :: crd(3*abfqmmm_param%natom)

  _REAL_ :: center(3), radius, com(3)
  _REAL_ :: r_in, r_out, r
  _REAL_ :: vrest, frest
  _REAL_ :: totalmass
  integer :: i, j, k
  integer :: res, resstart, resstop

   if(abfqmmm_param%r_diff_in == 0.0d0) then
    ! calculate center of center list and radius of sphere of (core extended + qm extended) list
    call abfqmmm_center_and_radius_of_atom_list(abfqmmm_param%center_type, abfqmmm_param%n_center, abfqmmm_param%center, &
                                                abfqmmm_param%n_qm, abfqmmm_param%qm, crd, center, radius)
    r_in = radius
   else
    r_in = abfqmmm_param%r_diff_in
   end if

   if(abfqmmm_param%r_diff_out == 0.0d0) then
    ! calculate center of center list and radius of sphere of buffer list
    call abfqmmm_center_and_radius_of_atom_list(abfqmmm_param%center_type, abfqmmm_param%n_center, abfqmmm_param%center, &
                                                abfqmmm_param%n_buffer, abfqmmm_param%buffer, crd, center, radius)

    r_out = radius
   else
    r_out = abfqmmm_param%r_diff_out
   end if

   vrest = 0.0d0
   do res=1,abfqmmm_param%nres
    if(abfqmmm_param%res_atom_number(res) > abfqmmm_param%solvent_atom_number) cycle
    if(abfqmmm_param%id(abfqmmm_param%res_pointers(res)) < 5) cycle
    resstart=abfqmmm_param%res_pointers(res)
    resstop=abfqmmm_param%res_pointers(res+1)-1
    ! calculate center of mass of solvent res
    com(:) = 0.0d0
    totalmass = 0.0d0
    do i=resstart,resstop
     totalmass = totalmass + abfqmmm_param%mass(i)
     j=3*(i-1)
     do k=1,3
      com(k) = com(k) + abfqmmm_param%mass(i)*crd(j+k)
     end do
    end do
    com(:) = com(:) / totalmass
    ! calculate potential and -grad of potential
    r = abfqmmm_point_point_dist(center, com)
!   vrest = vrest + abfqmmm_diffusion_potential(r, r_out)
    frest = abfqmmm_diffusion_force(r, r_out)
    ! calculate the forces on each atom of solvent res
    do i=resstart,resstop
     j=3*(i-1)
     do k=1,3
      abfqmmm_param%f(j+k) = abfqmmm_param%f(j+k) + abfqmmm_param%mass(i) / totalmass * frest * (com(k)-center(k)) / r
     end do
    end do
   end do

   abfqmmm_param%diff_potential = vrest

end subroutine abfqmmm_diffusion_restraint


!function abfqmmm_diffusion_potential(r, r_in, r_out) result(potential)
!
!use constants, only : PI
!
!  implicit none
!
!  _REAL_, intent(in) :: r
!  _REAL_, intent(in) :: r_in
!  _REAL_, intent(in) :: r_out
!
!  _REAL_ :: potential
!
!  potential = 0.0d0 
!! if( r < r_in ) then
!!  potential = 0.0d0
!! else if( r > r_out ) then
!!  potential = abfqmmm_param%diff_k
!! else
!!  potential = abfqmmm_param%diff_k * 0.50d0 * ( 1.d0 - cos( PI*(r-r_in)/(r_out-r_in) ) )
!! end if
!  
!end function abfqmmm_diffusion_potential


function abfqmmm_diffusion_force(r, r_out) result(force)

use constants, only : PI

  implicit none

  _REAL_, intent(in) :: r
! _REAL_, intent(in) :: r_in
  _REAL_, intent(in) :: r_out

  _REAL_ :: force

  if( r < r_out ) then
   force = -abfqmmm_param%diff_k
  else
   force = 0.0d0
  end if
! if( r < r_in ) then
!  force = 0.0d0
! else if( r > r_out ) then
!  force = 0.0d0
! else
!  force = -abfqmmm_param%diff_k * 0.50d0 * sin( PI*(r-r_in)/(r_out-r_in) ) * PI / (r_out-r_in)
! end if

end function abfqmmm_diffusion_force


subroutine abfqmmm_calc_diff_coeff(crd)

  implicit none

#include "box.h"
#include "../include/md.h"

  _REAL_, intent(in) :: crd(3*abfqmmm_param%natom)

  integer :: res, resstart, resstop
  integer :: n, i, j, k
  integer :: ier=0

  _REAL_ :: d(3)
  _REAL_ :: a, b, x, y, xx, xy

  if(abfqmmm_param%qmstep == 1) then
    ! find O atoms in water molecules and allocate arrays
    n=0
    do res=1, abfqmmm_param%nres
      if(abfqmmm_param%res_atom_number(res) > abfqmmm_param%solvent_atom_number) cycle
      resstart=abfqmmm_param%res_pointers(res)
      resstop=abfqmmm_param%res_pointers(res+1)-1
      do i=resstart, resstop
        if(abfqmmm_param%mass(i) < 4.0d0) cycle
        n=n+1
      end do
    end do
    
    abfqmmm_param%water_o_n = n
    
    allocate(abfqmmm_param%water_o_id(abfqmmm_param%water_o_n), stat=ier)
    REQUIRE(ier==0)

    allocate(abfqmmm_param%water_o_is_qm(abfqmmm_param%water_o_n), stat=ier)
    REQUIRE(ier==0)

    allocate(abfqmmm_param%water_o_time(abfqmmm_param%water_o_n), stat=ier)
    REQUIRE(ier==0)

    n=0
    do res=1, abfqmmm_param%nres
      resstart=abfqmmm_param%res_pointers(res)
      resstop=abfqmmm_param%res_pointers(res+1)-1
      if(abfqmmm_param%res_atom_number(res) > abfqmmm_param%solvent_atom_number) cycle
      do i=resstart, resstop 
        if(abfqmmm_param%mass(i) < 4.0d0) cycle
        n=n+1
        abfqmmm_param%water_o_id(n) = i
        abfqmmm_param%water_o_time(n) = 0
        if(abfqmmm_param%id(i) <= 4) then
          abfqmmm_param%water_o_is_qm(n) = 1
        else
          abfqmmm_param%water_o_is_qm(n) = 0
        end if
      end do
    end do

    allocate(abfqmmm_param%x0(3*abfqmmm_param%natom), stat=ier)
    REQUIRE(ier==0)
    abfqmmm_param%x0(:) = crd(:)

    allocate(abfqmmm_param%dx(3*abfqmmm_param%natom), stat=ier)
    REQUIRE(ier==0)
    abfqmmm_param%dx(:) = 0.0d0

    allocate(abfqmmm_param%ndiffmm(abfqmmm_param%maxqmstep), stat=ier)
    REQUIRE(ier==0)

    allocate(abfqmmm_param%diffmm(abfqmmm_param%maxqmstep), stat=ier)
    REQUIRE(ier==0)

    allocate(abfqmmm_param%ndiffqm(abfqmmm_param%maxqmstep), stat=ier)
    REQUIRE(ier==0)

    allocate(abfqmmm_param%diffqm(abfqmmm_param%maxqmstep), stat=ier)
    REQUIRE(ier==0)

    abfqmmm_param%ndiffmm = 0
    abfqmmm_param%ndiffqm = 0

    abfqmmm_param%diffmm = 0.0d0
    abfqmmm_param%diffqm = 0.0d0

    return
  end if

  do n=1,abfqmmm_param%water_o_n
     i = abfqmmm_param%water_o_id(n)
     j = (i-1)*3
     if(abfqmmm_param%id(i) <= 4) then
       if(abfqmmm_param%water_o_is_qm(n) == 1) then
         abfqmmm_param%water_o_time(n) = abfqmmm_param%water_o_time(n) + 1
         abfqmmm_param%ndiffqm(abfqmmm_param%water_o_time(n)) = abfqmmm_param%ndiffqm(abfqmmm_param%water_o_time(n)) + 1
         do k=1,3
           d(k) = crd(j+k) - abfqmmm_param%x0(j+k)
           d(k) = d(k) - nint(d(k)/box(k))*box(k)
           abfqmmm_param%dx(j+k) = abfqmmm_param%dx(j+k) + d(k)
           abfqmmm_param%diffqm(abfqmmm_param%water_o_time(n)) = abfqmmm_param%diffqm(abfqmmm_param%water_o_time(n)) + &
                                                                 abfqmmm_param%dx(j+k)*abfqmmm_param%dx(j+k)
         end do
       else
         abfqmmm_param%water_o_is_qm(n) = 1
         abfqmmm_param%water_o_time(n) = 0
         do k=1,3
            abfqmmm_param%dx(j+k) = 0.0d0
         end do
       end if
     else
       if(abfqmmm_param%water_o_is_qm(n) == 0) then
         abfqmmm_param%water_o_time(n) = abfqmmm_param%water_o_time(n) + 1
         abfqmmm_param%ndiffmm(abfqmmm_param%water_o_time(n)) = abfqmmm_param%ndiffmm(abfqmmm_param%water_o_time(n)) + 1
         do k=1,3
           d(k) = crd(j+k) - abfqmmm_param%x0(j+k)
           d(k) = d(k) - nint(d(k)/box(k))*box(k)
           abfqmmm_param%dx(j+k) = abfqmmm_param%dx(j+k) + d(k)
           abfqmmm_param%diffmm(abfqmmm_param%water_o_time(n)) = abfqmmm_param%diffmm(abfqmmm_param%water_o_time(n)) + &
                                                                 abfqmmm_param%dx(j+k)*abfqmmm_param%dx(j+k)
         end do
       else
         abfqmmm_param%water_o_is_qm(n) = 0
         abfqmmm_param%water_o_time(n) = 0
         do k=1,3
            abfqmmm_param%dx(j+k) = 0.0d0
         end do
       end if
     end if
     do k=1,3
        abfqmmm_param%x0(j+k) = crd(j+k)
     end do
  end do

  n = 0
  x = 0.0d0
  y = 0.0d0
  xx = 0.0d0
  xy = 0.0d0
  do i=1,abfqmmm_param%maxqmstep
    if(abfqmmm_param%ndiffqm(i) == 0) exit
    n = n+1
    a = real(i)*dt
    b = abfqmmm_param%diffqm(i)/real(abfqmmm_param%ndiffqm(i))/6.0d0/dt/real(i)
    x = x + a
    y = y + b
    xx = xx + a*a
    xy = xy + a*b
  end do

  n = 0
  x = 0.0d0
  y = 0.0d0
  xx = 0.0d0
  xy = 0.0d0
  do i=1,abfqmmm_param%maxqmstep
    if(abfqmmm_param%ndiffmm(i) == 0) exit
    n = n+1
    a = real(i)*dt
    b = abfqmmm_param%diffmm(i)/real(abfqmmm_param%ndiffmm(i))/6.0d0/dt/real(i)
    x = x + a
    y = y + b
    xx = xx + a*a
    xy = xy + a*b
  end do

end subroutine abfqmmm_calc_diff_coeff

! Does not seem to be used currently
!subroutine abfqmmm_write_outputs(x,v,nsp)
!
!use file_io_dat
!use molecule, only: n_iwrap_mask_atoms, iwrap_mask_atoms
!use stack
!
!  implicit none
!
!#include "../include/md.h"
!#include "box.h"
!#include "../include/memory.h"
!
!  _REAL_, intent(inout) :: x(*)
!  _REAL_, intent(inout) :: v(*)
!  integer, intent(inout) :: nsp(*)
!  
!  integer :: m, nrx, nr, nr3, l_temp
!  logical :: itdump, ivdump, ifdump, ixdump, loutfm
!  character(kind=1,len=5) :: routine="runmd"
!  _REAL_ :: box_center(3)
!
!  loutfm = ioutfm <= 0
!  nr=nrp
!  nr3=3*natom
!  nrx=nr3
!  if(ntwprt > 0) nrx = ntwprt*3
!
!  itdump = .false.             ! Write coordinates this step?
!  ivdump = .false.             ! Write velocities this step?
!  ifdump = .false.             ! Write forces this step?
!  ixdump = .false.             ! Write restart this step?
!
!  if (ntwx>0) itdump = mod(abfqmmm_param%qmstep,ntwx) == 0 ! Trajectory coords
!  if (ntwv>0) ivdump = mod(abfqmmm_param%qmstep,ntwv) == 0 ! Velocity
!  if (ntwf>0) ifdump = mod(abfqmmm_param%qmstep,ntwf) == 0 ! Force
!  if (ntwr /= 0) then
!     if ( mod(abfqmmm_param%qmstep,ntwr) == 0 ) ixdump = .true. ! Restart
!  endif
!  if (abfqmmm_param%qmstep==abfqmmm_param%maxqmstep) ixdump = .true. ! Final restart
!  if (ntwv == -1 .and. itdump) ivdump = .true. !Combined crdvel file
!
!  ! rst archive
!  if(ixdump) then
!   if(iwrap == 0) then
!    call mdwrit(abfqmmm_param%qmstep,nr,ntxo,ntb,x,v,t,temp0)
!   else if(iwrap == 1) then
!   ! --- use temp. array to hold coords. so that the master's values
!   !     are always identical to those on all other nodes:
!    call get_stack(l_temp,nr3,routine)
!    if(.not. rstack_ok)then
!     deallocate(r_stack)
!     allocate(r_stack(1:lastrst),stat=alloc_ier)
!     call reassign_rstack(routine)
!    endif
!    REQUIRE(rstack_ok)
!    do m=1,nr3
!     r_stack(l_temp+m-1) = x(m)
!    end do
!    call wrap_molecules(nspm,nsp,r_stack(l_temp)) 
!    if (ifbox == 2) call wrap_to(nspm,nsp,r_stack(l_temp),box)
!    call mdwrit(abfqmmm_param%qmstep,nr,ntxo,ntb,r_stack(l_temp),v,t,temp0)
!   else if (iwrap == 2) then
!    ! GMS ------------------------------------------
!    ! We are wrapping around a pre-determined mask
!    ! Need to center it on the mask COM first, then
!    ! wrap it normally as it happens on the iwrap=1 
!    ! case.
!    ! GMS ------------------------------------------
!    call get_stack(l_temp,nr3,routine)
!    if(.not. rstack_ok)then
!     deallocate(r_stack)
!     allocate(r_stack(1:lastrst),stat=alloc_ier)
!     call reassign_rstack(routine)
!    endif
!    REQUIRE(rstack_ok)
!    do m=1,nr3
!     r_stack(l_temp+m-1) = x(m)
!    end do
!    ! Now, wrap the coordinates around the iwrap_mask:
!    call iwrap2(n_iwrap_mask_atoms,iwrap_mask_atoms,r_stack(l_temp),box_center)
!    call mdwrit(abfqmmm_param%qmstep,nr,ntxo,ntb,r_stack(l_temp),v,t,temp0)
!    call free_stack(l_temp,routine)
!   end if
!  end if
!
!   ! crd archive
!   if(itdump) then
!    if(iwrap == 0) then
!     call corpac(x,1,nrx,MDCRD_UNIT,loutfm)
!     if(ntb > 0)  call corpac(box,1,3,MDCRD_UNIT,loutfm)    
!    else if(iwrap == 1) then
!     call get_stack(l_temp,nr3,routine)
!     if(.not. rstack_ok)then
!      deallocate(r_stack)
!      allocate(r_stack(1:lastrst),stat=alloc_ier)
!      call reassign_rstack(routine)
!     endif
!     REQUIRE(rstack_ok)
!     do m=1,nr3
!      r_stack(l_temp+m-1) = x(m)
!     end do
!     call wrap_molecules(nspm,nsp,r_stack(l_temp))
!     if (ifbox == 2) call wrap_to(nspm,nsp,r_stack(l_temp),box)
!     call corpac(r_stack(l_temp),1,nrx,MDCRD_UNIT,loutfm)
!     call corpac(box,1,3,MDCRD_UNIT,loutfm)
!     call free_stack(l_temp,routine)
!    else if (iwrap == 2) then
!    ! GMS ------------------------------------------
!    ! We are wrapping around a pre-determined mask
!    ! Need to center it on the mask COM first, then
!    ! wrap it normally as it happens on the iwrap=1 
!    ! case.
!    ! GMS ------------------------------------------
!     call get_stack(l_temp,nr3,routine)
!     if(.not. rstack_ok)then
!      deallocate(r_stack)
!      allocate(r_stack(1:lastrst),stat=alloc_ier)
!      call reassign_rstack(routine)
!     endif
!     REQUIRE(rstack_ok)
!     do m=1,nr3
!      r_stack(l_temp+m-1) = x(m)
!     end do
!     call iwrap2(n_iwrap_mask_atoms,iwrap_mask_atoms, r_stack(l_temp),box_center)
!     call corpac(r_stack(l_temp),1,nrx,MDCRD_UNIT,loutfm)
!     call corpac(box,1,3,MDCRD_UNIT,loutfm)
!     call free_stack(l_temp,routine)
!    end if
!   end if
!
!   ! vel archive
!   if(ivdump) then
!    call corpac(v,1,nrx,MDVEL_UNIT,loutfm)
!   end if
!
!   ! frc archive
!   if(ifdump) then
!    call corpac(abfqmmm_param%f,1,nrx,MDFRC_UNIT,loutfm)
!   end if
!
!end subroutine abfqmmm_write_outputs


end module abfqmmm_module
