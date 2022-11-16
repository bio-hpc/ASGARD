module amoeba_parm
implicit none

private
integer, save :: num_atoms = 0
integer, save :: num_xyz_atoms = 0
integer, save :: num_dyn_atoms = 0
integer,save  :: num_res = 0
integer, save :: num_molecules = 0
integer, save :: num_12 = 0
integer, save :: num_13 = 0
integer, save :: num_14 = 0
integer, save :: num_15 = 0
integer, save :: num_polarizable = 0
integer, save :: num_polargroup = 0

integer, save, allocatable :: atomtype(:),atomclass(:),atomic_number(:), &
              atomvalence(:),num_atoms_in_molecule(:)
double precision, save, allocatable ::atomic_weight(:) 
double precision, save, allocatable ::polarizability(:)

double precision, save, allocatable :: xyz_crd(:,:) 

logical,save :: dyn_read = .FALSE.
double precision, save, allocatable :: dyn_crd(:,:) 
double precision, save, allocatable :: dyn_vel(:,:) 
double precision, save, allocatable :: dyn_accel(:,:) 
double precision, save, allocatable :: dyn_old_accel(:,:) 
double precision, save :: dyn_cell_length(3) 
double precision, save :: dyn_cell_angles(3) 
double precision, parameter :: sander_to_tinker_time_convert = 20.455d0
double precision, parameter :: tinker_to_sander_time_convert = 1.d0 / 20.455d0

character(len=4),allocatable :: atomname(:)
character(len=4),allocatable :: atomsymbol(:)
character(len=4),allocatable :: reslabel(:)
integer,allocatable :: startres(:)
character(len=80)title

integer,save  :: num_bonds = 0
integer,save, allocatable :: bond_list(:,:)
integer,save :: num_bond_params
double precision,save, allocatable ::bond_param(:,:)

integer,save, allocatable :: vdw_atom_parent(:)
double precision,save, allocatable :: vdw_parent_weight(:)
integer,save, allocatable :: atom_vdw_type(:)
integer,save :: num_vdw_types=0
double precision, parameter :: vdw_buffer_delta=0.07d0
double precision, parameter :: vdw_buffer_gamma=0.12d0
double precision,save, allocatable :: vdw_mix_radius(:,:),vdw_mix_epsilon(:,:)


integer,save :: num_urey_bradley = 0
integer,save, allocatable :: urey_bradley_list(:,:)
integer, save :: num_ureyb_params
double precision,save, allocatable ::urey_bradley_param(:,:)

integer,save :: num_tot_angles
integer,save :: num_regular_angles = 0,num_regular_angle_params
integer,save :: num_trigonal_angles = 0,num_trigonal_angle_params
integer, allocatable,save :: regular_angle_list(:,:)
double precision,save, allocatable ::regular_angle_param(:,:) 
integer, allocatable :: trigonal_angle_list(:,:),trigonal_angle_parm_ptr(:)
double precision, allocatable ::trigonal_angle_param(:,:) 

integer,save :: num_out_of_plane_bends = 0,num_opbend_params
integer,save, allocatable :: outofplane_list(:,:)
double precision,save, allocatable ::outofplane_param(:)

integer,save :: num_torsions = 0,num_torsion_params
integer,save, allocatable :: torsion_list(:,:)
double precision,save, allocatable ::torsion_param(:,:)

integer,save :: num_str_torsions = 0,num_str_torsion_params
integer,save, allocatable :: str_torsion_list(:,:)
double precision,save, allocatable ::str_torsion_param(:,:)

integer,save :: num_pitorsions = 0,num_pitorsion_params
integer,save, allocatable :: pitorsion_list(:,:)
double precision,save, allocatable ::pitorsion_param(:,:)

integer,save :: num_stretch_bends = 0,num_strbend_params
integer,save, allocatable :: stretch_bend_list(:,:)
double precision,save, allocatable ::stretch_bend_param(:,:)

type :: angle_angle_functable
   integer :: dim1 = 0,dim2 = 0
   double precision,pointer :: angle1(:) => null()
   double precision,pointer :: angle2(:) => null()
   double precision, pointer :: func(:,:) => null()
   double precision, pointer :: dfunc_dangle1(:,:) => null()
   double precision, pointer :: dfunc_dangle2(:,:) => null()
   double precision, pointer :: d2func_dangle1_dangle2(:,:) => null()
end type  angle_angle_functable
integer, save :: num_tor_tors,num_tortor_tables
integer,save, allocatable :: tortor_list(:,:)
type(angle_angle_functable),save, allocatable :: tortor_table(:)
logical, save :: get_tortors_from_analout

integer,save, allocatable :: bond_nghbors(:),offset_nghbors(:)
integer,save, allocatable :: nghb_12(:,:),nghb_13(:,:),nghb_14(:,:),nghb_15(:,:)
integer,save, allocatable :: nghb_polargroup(:),num_in_polar_group(:), &
                             offset_polar_group(:)
integer,save, allocatable :: num_excluded(:),excluded_atomlist(:)
integer,save, allocatable :: adjust_list(:,:)
integer,save :: size_excluded_list
integer, parameter :: bit_12=1,bit_13=2,bit_14=3,bit_15=4,bit_polgrp=0

double precision,parameter  ::  pi= 3.14159265358979323846d0
double precision,parameter :: radians_to_degrees = 180.d0/pi
double precision,parameter :: degrees_to_radians = pi / 180.d0

type :: frame_list_entry
   integer :: frame_index ! which atomic frame does this refer to
   integer :: frame_point_number ! which frame point (1 or 2)
   integer :: vector_tail_index ! unit vector to add into point def
   integer :: vector_head_index
   integer :: num_vectors ! number of unit vec contribs to frame def point
end type frame_list_entry
integer, save :: num_frame_list
type(frame_list_entry),allocatable,save :: frame_list(:)
integer, save :: num_multipoles
double precision, save, allocatable :: local_multipole(:,:)

character(len=80),save :: realformat = '(5E16.8)'
character(len=80),save :: outrealformat = '%FORMAT(5E16.8)'
type :: chiral_frame
  integer :: frame_index
  integer :: fourth_atom
  integer :: chirality
end type chiral_frame
integer,save :: num_chiral_frame_list
type(chiral_frame), save, allocatable :: chiral_frame_list(:)

double precision, save :: boxa = 0.d0,boxb = 0.d0,boxc = 0.d0, &
               box_alpha = 0.d0,box_beta = 0.d0,box_gamma = 0.d0

public AM_PARM_get_atom_params,AM_PARM_get_bond_params,  &
       AM_PARM_get_angle_params,AM_PARM_get_strbend_params,  &
       AM_PARM_get_UreyB_params,AM_PARM_get_outofplane_params, &
       AM_PARM_get_torsion_params,AM_PARM_bond_nghb_list, &
       AM_PARM_find_vdw_parent,AM_PARM_get_atom_vdw_params, &
       AM_PARM_get_pitorsion_params,AM_PARM_get_tortor_params, &
       AM_PARM_deallocate,AM_PARM_get_xyz,AM_PARM_get_pdb, &
       AM_PARM_write_top,AM_PARM_write_inpcrd,AM_PARM_get_mol_list,  &
       AM_PARM_write_bond_params,AM_PARM_get_polar_params, &
       AM_PARM_get_nghb_12131415_lists,AM_PARM_get_excluded_atom_list, &
       AM_PARM_write_UreyB_params,AM_PARM_write_angle_params, &
       AM_PARM_write_outofplane_params,AM_PARM_write_torsion_params, &
       AM_PARM_write_pitorsion_params,AM_PARM_write_strbend_params, &
       AM_PARM_write_tortor_params,AM_PARM_write_multipole, &
       AM_PARM_get_multipole,AM_PARM_process_key_file, &
       AM_PARM_write_adjust,AM_PARM_write_polar,AM_PARM_write_vdw,title, &
       AM_PARM_get_dyn,AM_PARM_check_atom_num, &
       AM_PARM_get_stretch_tor_params,AM_PARM_write_str_tor_params, &
       AM_PARM_set_real_format
contains
!------------------------------------------------------------------
subroutine AM_PARM_set_real_format()
      realformat = '(4E20.12)'
      outrealformat = '%FORMAT(4E20.12)'
end subroutine AM_PARM_set_real_format
!------------------------------------------------------------------
subroutine AM_PARM_write_top(prmtop_unit)
  integer,intent(in) :: prmtop_unit

  double precision :: r_zero = 0.d0,r_one = 1.d0,beta = 90.d0
  integer :: n,one=1,zero=0,nttyp,ntype,nspsol = 0,iptres = 0

  integer :: natom,ntypes,nbonh,mbona,ntheth,mtheta,nphih,mphia, &
         nhparm,nparm,nnb,nres,nbona,ntheta,nphia, &
         numbnd,numang,nptra,natyp,nphb,ifpert,nbper,ngper,ndper, &
         mbper,mgper,mdper,ifbox,nmxrs,ifcap,numextra,ncopy
  integer :: values(8)
  character(len=12) :: date,time,zone
  character(len=8)word,word1
  character(len=16)word2

  call date_and_time(date,time,zone,values)
  write(prmtop_unit,'(a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2)')  &
           '%VERSION  VERSION_STAMP = V0001.000  DATE = ', &
           values(2),'/',values(3),'/',values(1)-2000, &
           '  ',values(5),':',values(6),':',values(7)
  write(prmtop_unit,'(a)')'%FLAG TITLE'
  write(prmtop_unit,'(a)')'%FORMAT(a)'
  write(prmtop_unit,'(a)')title(1:len_trim(title))

  natom = num_atoms
  ntypes = 1
  nbonh = 1
  mbona = 1
  ntheth = 1
  mtheta = 1
  nphih = 1
  mphia = 1
  nhparm = 0
  nparm = 0
  nnb = size_excluded_list 
  nres = num_res
  nbona = mbona
  ntheta = mtheta
  nphia = mphia
  numbnd = 1
  numang = 1
  nptra = 1
  natyp = 1
  nphb = 1
! perturb
  ifpert = 0
  nbper = 0
  ngper = 0
  ndper = 0
  mbper = 0
  mgper = 0
  mdper = 0
! PBC
  ifbox = 1
  nmxrs = 0
  ifcap = 0
  numextra = 0
  ncopy = 0

  nttyp = ntypes*(ntypes+1)/2
  ntype = ntypes*ntypes
  write(prmtop_unit,'(a)')'%FLAG POINTERS'
  write(prmtop_unit,'(a)')'%FORMAT(10I8)'
  write(prmtop_unit,'(10I8)')  &
         natom,ntypes,nbonh,mbona,ntheth,mtheta,nphih,mphia, &
         nhparm,nparm,nnb,nres,nbona,ntheta,nphia, &
         numbnd,numang,nptra,natyp,nphb,ifpert,nbper,ngper,ndper, &
         mbper,mgper,mdper,ifbox,nmxrs,ifcap,numextra,ncopy

  write(word,'(I8)')num_atoms
  word1 = adjustl(word)
  word2 = '('//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%FLAG ATOM_NAME'
  !write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')'%FORMAT(20a4)'
  write(prmtop_unit,'(20a4)')(atomname(n),n=1,num_atoms)

  write(prmtop_unit,'(a)')'%FLAG CHARGE'
  !write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat)(r_zero,n=1,num_atoms)

  write(prmtop_unit,'(a)')'%FLAG MASS'
  !write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat)(atomic_weight(n),n=1,num_atoms)

  write(prmtop_unit,'(a)')'%FLAG ATOM_TYPE_INDEX'
  !write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')'%FORMAT(10I8)'
  write(prmtop_unit,'(10I8)')(one,n=1,num_atoms)

  write(prmtop_unit,'(a)')'%FLAG NUMBER_EXCLUDED_ATOMS'
  !write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')'%FORMAT(10I8)'
  write(prmtop_unit,'(10I8)')(num_excluded(n),n=1,num_atoms)

  write(word,'(I8)')ntype
  word1 = adjustl(word)
  word2 = '('//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%FLAG NONBONDED_PARM_INDEX'
  !write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')'%FORMAT(10I8)'
  write(prmtop_unit,'(10I8)')(one,n=1,ntype)

  write(word,'(I8)')num_res
  word1 = adjustl(word)
  word2 = '('//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%FLAG RESIDUE_LABEL'
  !write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')'%FORMAT(20a4)'
  write(prmtop_unit,'(20a4)')(reslabel(n),n=1,num_res)

  write(prmtop_unit,'(a)')'%FLAG RESIDUE_POINTER'
  !write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')'%FORMAT(10I8)'
  write(prmtop_unit,'(10I8)')(startres(n),n=1,num_res)

  write(word,'(I8)')numbnd
  word1 = adjustl(word)
  word2 = '('//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%FLAG BOND_FORCE_CONSTANT'
  !write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat)(r_zero,n=1,numbnd)

  write(prmtop_unit,'(a)')'%FLAG BOND_EQUIL_VALUE'
  !write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat)(r_zero,n=1,numbnd)

  write(word,'(I8)')numang
  word1 = adjustl(word)
  word2 = '('//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%FLAG ANGLE_FORCE_CONSTANT'
  !write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat)(r_zero,n=1,numang)

  write(prmtop_unit,'(a)')'%FLAG ANGLE_EQUIL_VALUE'
  !write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat)(r_zero,n=1,numang)

  write(word,'(I8)')nptra
  word1 = adjustl(word)
  word2 = '('//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%FLAG DIHEDRAL_FORCE_CONSTANT'
  !write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat)(r_zero,n=1,nptra)

  write(prmtop_unit,'(a)')'%FLAG DIHEDRAL_PERIODICITY'
  !write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat)(r_zero,n=1,nptra)

  write(prmtop_unit,'(a)')'%FLAG DIHEDRAL_PHASE'
  !write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat)(r_zero,n=1,nptra)

  write(word,'(I8)')natyp
  word1 = adjustl(word)
  word2 = '('//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%FLAG SOLTY'
  !write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat)(r_zero,n=1,natyp)

  write(word,'(I8)')nttyp
  word1 = adjustl(word)
  word2 = '('//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%FLAG LENNARD_JONES_ACOEF'
  !write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat)(r_one,n=1,nttyp)

  write(prmtop_unit,'(a)')'%FLAG LENNARD_JONES_BCOEF'
  !write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat)(r_one,n=1,nttyp)

  write(word,'(I8)')nbonh
  word1 = adjustl(word)
  word2 = '(3,'//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%FLAG BONDS_INC_HYDROGEN'
  !write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')'%FORMAT(10I8)'
  write(prmtop_unit,'(10I8)')(one,n=1,3)

  write(word,'(I8)')mbona
  word1 = adjustl(word)
  word2 = '(3,'//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%FLAG BONDS_WITHOUT_HYDROGEN'
  !write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')'%FORMAT(10I8)'
  write(prmtop_unit,'(10I8)')(one,n=1,3)

  write(word,'(I8)')ntheth
  word1 = adjustl(word)
  word2 = '(4,'//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%FLAG ANGLES_INC_HYDROGEN'
  !write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')'%FORMAT(10I8)'
  write(prmtop_unit,'(10I8)')(one,n=1,4)

  write(word,'(I8)')mtheta
  word1 = adjustl(word)
  word2 = '(4,'//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%FLAG ANGLES_WITHOUT_HYDROGEN'
  !write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')'%FORMAT(10I8)'
  write(prmtop_unit,'(10I8)')(one,n=1,4)

  write(word,'(I8)')nphih
  word1 = adjustl(word)
  word2 = '(5,'//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%FLAG DIHEDRALS_INC_HYDROGEN'
  !write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')'%FORMAT(10I8)'
  write(prmtop_unit,'(10I8)')(one,n=1,5)

  write(word,'(I8)')mphia
  word1 = adjustl(word)
  word2 = '(5,'//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%FLAG DIHEDRALS_WITHOUT_HYDROGEN'
  !write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')'%FORMAT(10I8)'
  write(prmtop_unit,'(10I8)')(one,n=1,5)

  write(word,'(I8)')nnb
  word1 = adjustl(word)
  word2 = '('//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%FLAG EXCLUDED_ATOMS_LIST'
  !write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')'%FORMAT(10I8)'
  write(prmtop_unit,'(10I8)')(excluded_atomlist(n),n=1,size_excluded_list)

  write(word,'(I8)')nphb
  word1 = adjustl(word)
  word2 = '('//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%FLAG HBOND_ACOEF'
  !write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat)(r_zero,n=1,nphb)

  write(prmtop_unit,'(a)')'%FLAG HBOND_BCOEF'
  !write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat)(r_zero,n=1,nphb)

  write(prmtop_unit,'(a)')'%FLAG HBCUT'
  !write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat)(r_one,n=1,nphb)

  write(word,'(I8)')natom
  word1 = adjustl(word)
  word2 = '('//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%FLAG AMBER_ATOM_TYPE'
  !write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')'%FORMAT(20a4)'
  write(prmtop_unit,'(20a4)')(atomsymbol(n),n=1,natom)

  write(prmtop_unit,'(a)')'%FLAG TREE_CHAIN_CLASSIFICATION'
  !write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')'%FORMAT(20a4)'
  write(prmtop_unit,'(20a4)')('BLA ',n=1,natom)

  write(prmtop_unit,'(a)')'%FLAG JOIN_ARRAY'
  !write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')'%FORMAT(10I8)'
  write(prmtop_unit,'(10I8)')(zero,n=1,natom)

  write(prmtop_unit,'(a)')'%FLAG IROTAT'
  !write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')'%FORMAT(10I8)'
  write(prmtop_unit,'(10I8)')(zero,n=1,natom)

  write(prmtop_unit,'(a)')'%FLAG SOLVENT_POINTERS'
  write(prmtop_unit,'(a)')'%FORMAT(10I8)'
  write(prmtop_unit,'(10I8)')iptres,num_molecules,nspsol

  write(word,'(I8)')num_molecules
  word1 = adjustl(word)
  word2 = '('//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%FLAG ATOMS_PER_MOLECULE'
  !write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')'%FORMAT(10I8)'
  write(prmtop_unit,'(10I8)')(num_atoms_in_molecule(n),n=1,num_molecules)

  write(prmtop_unit,'(a)')'%FLAG BOX_DIMENSIONS'
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat)beta,r_one,r_one,r_one

  write(prmtop_unit,'(a)')'%FLAG AMOEBA_FORCEFIELD'
  write(prmtop_unit,'(a,a)')'%COMMENT This indicates that ', &
                           'this parm file is specific to amoeba'
  write(prmtop_unit,'(a,a)')'%COMMENT This must be present if', &
                          ' do_amoeba(in mdin) is 1'
  write(prmtop_unit,'(a)')'%COMMENT This must NOT be present if do_amoeba is 0'
  write(prmtop_unit,'(a)')'%FORMAT(i5)'
  write(prmtop_unit,'(a)')'    1'

  write(word,'(I8)')num_atoms
  word1 = adjustl(word)
  word2 = '('//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_ATOM_TYPE_INDEX'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')'%FORMAT(10I8)'
  write(prmtop_unit,'(10I8)')(atomtype(n),n=1,num_atoms)

  write(word,'(I8)')num_atoms
  word1 = adjustl(word)
  word2 = '('//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_ATOMIC_NUMBER'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')'%FORMAT(10I8)'
  write(prmtop_unit,'(10I8)')(atomic_number(n),n=1,num_atoms)

  write(word,'(I8)')num_atoms
  word1 = adjustl(word)
  word2 = '('//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_ATOM_CLASS_INDEX'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')'%FORMAT(10I8)'
  write(prmtop_unit,'(10I8)')(atomclass(n),n=1,num_atoms)

end subroutine AM_PARM_write_top
!------------------------------------------------------------------
subroutine AM_PARM_write_inpcrd(inpcrd_unit)
  integer, intent(in) :: inpcrd_unit

  integer :: j,n
  double precision :: zero=0.d0
  integer :: values(8)
  character(len=12) :: date,time,zone
  character(len=8)word,word1
  character(len=16)word2

#ifdef BEEMAN_FORMAT
  call date_and_time(date,time,zone,values)
  write(inpcrd_unit,'(a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a)')  &
           '%VERSION  VERSION_STAMP = V0001.000  DATE = ', &
           values(2),'/',values(3),'/',values(1)-2000, &
           '  ',values(5),':',values(6),':',values(7), ' '
  write(inpcrd_unit,'(a)')'%FLAG TITLE'
  write(inpcrd_unit,'(a)')'%FORMAT(a)'
  write(inpcrd_unit,'(a)')title(1:len_trim(title))
  write(inpcrd_unit,'(a)')'%FLAG ATOMIC_COORDS_SIMULATION_TIME'
  write(inpcrd_unit,'(a)')'%FORMAT(E16.8)'
  write(inpcrd_unit,'(E16.8)')zero
!-----------------------------------------------------------
  write(inpcrd_unit,'(a)')'%FLAG ATOMIC_COORDS_NUM_LIST'
  write(inpcrd_unit,'(a)')'%FORMAT(i8)'
  write(inpcrd_unit,'(I8)')num_atoms
  write(word,'(I8)')num_atoms
  word1 = adjustl(word)
  word2 = '(3,'//word1(1:len_trim(word1))//')'
  write(inpcrd_unit,'(a)')'%FLAG ATOMIC_COORDS_LIST'
  write(inpcrd_unit,'(a)')'%COMMENT   dimension = '//word2
  write(inpcrd_unit,'(a)')'%FORMAT(3e20.12)'
!-----------------------------------------------------------
  if ( dyn_read )then
     write(inpcrd_unit,'(3e20.12)')((dyn_crd(j,n),j=1,3),n=1,num_atoms)
  else
     write(inpcrd_unit,'(3e20.12)')((xyz_crd(j,n),j=1,3),n=1,num_atoms)
  endif
  ! velocities, accelerations
  if ( dyn_read )then
     write(inpcrd_unit,'(a)')'%FLAG ATOMIC_VELOCITIES_NUM_LIST'
     write(inpcrd_unit,'(a)')'%FORMAT(i8)'
     write(inpcrd_unit,'(I8)')num_atoms
     write(inpcrd_unit,'(a)')'%FLAG ATOMIC_VELOCITIES_LIST'
     write(inpcrd_unit,'(a)')'%COMMENT   dimension = '//word2
     write(inpcrd_unit,'(a)')'%FORMAT(3e20.12)'
     write(inpcrd_unit,'(3e20.12)') &
       ((tinker_to_sander_time_convert*dyn_vel(j,n),j=1,3),n=1,num_atoms)
     write(inpcrd_unit,'(a)')'%FLAG ATOMIC_ACCELERATIONS_NUM_LIST'
     write(inpcrd_unit,'(a)')'%FORMAT(i8)'
     write(inpcrd_unit,'(I8)')num_atoms
     write(inpcrd_unit,'(a)')'%FLAG ATOMIC_ACCELERATIONS_LIST'
     write(inpcrd_unit,'(a)')'%COMMENT   dimension = '//word2
     write(inpcrd_unit,'(a)')'%FORMAT(3e20.12)'
     write(inpcrd_unit,'(3e20.12)') &
       ((tinker_to_sander_time_convert**2*dyn_accel(j,n),j=1,3),n=1,num_atoms)
     write(inpcrd_unit,'(a)')'%FLAG OLD_ATOMIC_ACCELERATIONS_NUM_LIST'
     write(inpcrd_unit,'(a)')'%FORMAT(i8)'
     write(inpcrd_unit,'(I8)')num_atoms
     write(inpcrd_unit,'(a)')'%FLAG OLD_ATOMIC_ACCELERATIONS_LIST'
     write(inpcrd_unit,'(a)')'%COMMENT   dimension = '//word2
     write(inpcrd_unit,'(a)')'%FORMAT(3e20.12)'
     write(inpcrd_unit,'(3e20.12)') &
       ((tinker_to_sander_time_convert**2*dyn_old_accel(j,n),j=1,3), &
                                                             n=1,num_atoms)
  endif
!-----------------------------------------------------------
  write(inpcrd_unit,'(a)')'%FLAG UNIT_CELL_PARAMETERS'
  write(inpcrd_unit,'(a,a)')'%COMMENT lengths a,b,c; ', &
                          'then angles alpha,beta,gamma'
  write(inpcrd_unit,'(a)')'%FORMAT(3e20.12)'
  if ( dyn_read )then
     write(inpcrd_unit,'(3e20.12)') &
         (dyn_cell_length(j),j=1,3),(dyn_cell_angles(j),j=1,3)
  else
     write(inpcrd_unit,'(3e20.12)')boxa,boxb,boxc,box_alpha,box_beta,box_gamma
  endif
#else
  write(inpcrd_unit,'(a)')title(1:len_trim(title))
  write(inpcrd_unit,'(i6)') num_atoms
  write(inpcrd_unit,'(6f12.7)') ((xyz_crd(j,n),j=1,3),n=1,num_atoms)
  write(inpcrd_unit,'(6f12.7)') boxa,boxb,boxc,box_alpha,box_beta,box_gamma
  close(inpcrd_unit)
#endif
  return
!-----------------------------------------------------------
end subroutine AM_PARM_write_inpcrd
!------------------------------------------------------------------
subroutine AM_PARM_check_atom_num()
  if ( num_atoms /= num_xyz_atoms )then
      write(6,*)'mismatch between num_atoms,num_xyz_atoms: ', &
                 num_atoms,num_xyz_atoms
      stop
  endif
  if ( dyn_read )then
     if ( num_dyn_atoms /= num_xyz_atoms )then
         write(6,*)'mismatch between num_dyn_atoms,num_xyz_atoms: ', &
                 num_dyn_atoms,num_xyz_atoms
         stop
     endif
  endif
end subroutine AM_PARM_check_atom_num
!------------------------------------------------------------------
subroutine AM_PARM_get_xyz(xyz_unit)
  integer, intent(in) :: xyz_unit

  character(len=4) :: aname
  character(len=1024) :: line
  integer :: k,n,ier,n_lines
  logical :: box_read = .false.
  double precision :: a, b, c, alpha, beta, gamma
  read(xyz_unit,*)num_xyz_atoms
  write(6,*)'num_xyz_atoms = ',num_xyz_atoms
  allocate(xyz_crd(3,num_xyz_atoms), stat=ier)
  if ( ier /= 0 )then
    write(6,*)' AM_PARM_get_xyz: problem allocating xyz_crd'
    stop
  endif
  n = 0
  n_lines = 1
  do while (n < num_xyz_atoms)
    read(xyz_unit,'(A)') line
    n_lines = n_lines + 1
    if (n_lines == 2 .and. .not. box_read) then
      ! Try to read in the box
      read(line,*,err=10) a, b, c, alpha, beta, gamma
      ! If we didn't jump to 10, we read the box and now overwrite the main box
      ! parameters with the ones we read in here.
      box_read = .true.
      boxa = a
      boxb = b
      boxc = c
      box_alpha = alpha
      box_beta = beta
      box_gamma = gamma
      cycle
10    continue
      ! If we are here, the box reading failed, so just go on
    end if
    n = n + 1
    read(line,*,iostat=ier)k,aname,xyz_crd(1,n),xyz_crd(2,n),xyz_crd(3,n)
    if (ier /= 0) then
      write(0,*) 'Could not understand line ', n_lines, ' in xyz file'
      write(0,'(a)') '   ', trim(line)
      stop
    end if
  enddo
  ! Make sure we have box coordinates set here (it used to be an error when
  ! reading the keyword file)
  if (boxa == 0.d0 .and. boxb == 0.d0 .and. boxc == 0.d0) then
    write(0,*) 'No box information in XYZ or .key file'
    stop
  end if
end subroutine AM_PARM_get_xyz
!------------------------------------------------------------------
subroutine AM_PARM_get_dyn(dyn_unit)
  integer, intent(in) :: dyn_unit

  character(len=120)line
  integer j,n,ier1,ier2,ier3,ier4
  read(dyn_unit,'(a)')line
  read(dyn_unit,*)num_dyn_atoms
  read(dyn_unit,'(a)')line
  read(dyn_unit,*)(dyn_cell_length(j),j=1,3)
  read(dyn_unit,*)(dyn_cell_angles(j),j=1,3)

  allocate(dyn_crd(3,num_xyz_atoms),stat=ier1)
  allocate(dyn_vel(3,num_xyz_atoms),stat=ier1)
  allocate(dyn_accel(3,num_xyz_atoms),stat=ier1)
  allocate(dyn_old_accel(3,num_xyz_atoms),stat=ier1)
  read(dyn_unit,'(a)')line
  do n = 1,num_dyn_atoms
     read(dyn_unit,*)(dyn_crd(j,n),j=1,3)
  enddo
  read(dyn_unit,'(a)')line
  do n = 1,num_dyn_atoms
     read(dyn_unit,*)(dyn_vel(j,n),j=1,3)
  enddo
  read(dyn_unit,'(a)')line
  do n = 1,num_dyn_atoms
     read(dyn_unit,*)(dyn_accel(j,n),j=1,3)
  enddo
  read(dyn_unit,'(a)')line
  do n = 1,num_dyn_atoms
     read(dyn_unit,*)(dyn_old_accel(j,n),j=1,3)
  enddo
  dyn_read = .TRUE.
end subroutine AM_PARM_get_dyn
!------------------------------------------------------------------
subroutine AM_PARM_get_pdb(pdb_unit)
  integer, intent(in) :: pdb_unit

  character(len=120)line,word
  character(len=4),allocatable :: tmpatomname(:)
  character(len=3),allocatable :: tmpresname(:)
  character(len=1),allocatable :: tmpchainname(:)
  integer,allocatable :: tmpresnum(:)
  double precision,allocatable :: tmpx(:),tmpy(:),tmpz(:)
  integer,allocatable :: xyz_to_pdb(:),pdb_to_xyz(:)
  character(len=6) :: string
  integer :: num_pdb_atoms = 0,num,numres,curr_res,m,n,p
  integer :: ios,ier1,ier2,ier3,ier4,ier5,ier6,ier7,ier8,anum
  double precision :: x,y,z
  logical :: warned = .false.

  rewind(pdb_unit)
  do
    read(pdb_unit,'(A)',iostat=ios)line
    if ( ios /= 0 )exit
    if ( line(1:6) == 'ATOM  ' .or. line(1:6) == 'HETATM')then
      num_pdb_atoms = num_pdb_atoms + 1
    endif
  enddo
  if ( num_pdb_atoms /= num_xyz_atoms )then
    write(6,*)'num pdb atoms /= num_xyz_atoms !!'
    stop
  endif
  allocate(tmpatomname(num_pdb_atoms),stat=ier1)
  allocate(atomname(num_pdb_atoms),stat=ier2)
  allocate(tmpresname(num_pdb_atoms),stat=ier3)
  allocate(tmpchainname(num_pdb_atoms),stat=ier4)
  allocate(tmpx(num_pdb_atoms),stat=ier5)
  allocate(tmpy(num_pdb_atoms),stat=ier6)
  allocate(tmpz(num_pdb_atoms),stat=ier7)
  allocate(tmpresnum(num_pdb_atoms),stat=ier8)
  if ( (ier1 /= 0) .or. (ier2 /= 0) .or. (ier3 /= 0) .or. (ier4 /= 0) .or. &
       (ier5 /= 0) .or. (ier6 /= 0) .or. (ier7 /= 0) .or. (ier8 /= 0) )then
    write(6,*)'AM_PARM_get_pdb: prob allocating'
    stop
  endif
  rewind(pdb_unit)
  num = 0
  do
    read(pdb_unit,'(A)',iostat=ios)line
    if ( ios /= 0 )exit
    if ( line(1:6) == 'ATOM  ' .or. line(1:6) == 'HETATM')then
      num = num + 1
      read(line,'(a6,i5,1x,a4,1x,a3,1x,a1,i4,4x,3f8.3)') &
                 string,anum,tmpatomname(num),tmpresname(num), &
                 tmpchainname(num),tmpresnum(num), &
                 tmpx(num),tmpy(num),tmpz(num)
    endif
  enddo
  ! get residues info
  numres = 0
  curr_res = 0
  do num = 1,num_pdb_atoms
    if ( tmpresnum(num) /= curr_res )then
      curr_res = tmpresnum(num)
      numres = numres + 1
    endif
  enddo
  write(6,*)'numres = ',numres
  num_res = numres
  ! allocate
  allocate(reslabel(numres),stat=ier1)
  allocate(startres(numres+1),stat=ier2)
  allocate(xyz_to_pdb(num_xyz_atoms),stat=ier3)
  allocate(pdb_to_xyz(num_xyz_atoms),stat=ier4)
  if ( (ier1 /= 0) .or. (ier2 /= 0) .or. (ier3 /= 0) .or. (ier4 /= 0) )then
    write(6,*)'problem allocating residue info'
    stop
  endif
  numres = 0
  curr_res = 0
  do num = 1,num_pdb_atoms
    if ( tmpresnum(num) /= curr_res )then
      curr_res = tmpresnum(num)
      numres = numres + 1
      startres(numres) = num
      reslabel(numres)(1:3) = tmpresname(num)
      reslabel(numres)(4:4) = ' '
    endif
  enddo
  startres(numres+1) = num_pdb_atoms+1
  ! match the atom names for xyz atoms
  do num = 1,num_pdb_atoms
    xyz_to_pdb(num) = 0
    pdb_to_xyz(num) = 0
  enddo
  do n = 1,numres
    do m = startres(n),startres(n+1)-1
      x = xyz_crd(1,m)
      y = xyz_crd(2,m)
      z = xyz_crd(3,m)
      call TPDB_match(x,y,z,tmpx,tmpy,tmpz,startres(n),startres(n+1)-1,p)
      if ( p == 0 )then
        ! In this case, the conformations must be different. So just assume both
        ! files have the same ordering
        if (.not. warned) then
          write(0,*)'TPDB_read_file: xyz atom ',m,' fails to match pdb atom'
          write(0,*)'Assuming atom ordering is identical in PDB and XYZ'
          warned = .true.
        end if
        p = m
      end if
        ! rotate hydrogens
      atomname(m)(1:1) = tmpatomname(p)(2:2)
      if ( tmpatomname(p)(3:3) == ' ' )then
        atomname(m)(2:2) = tmpatomname(p)(1:1)
        atomname(m)(3:4) = '  '
      else
        atomname(m)(2:2) = tmpatomname(p)(3:3)
        if ( tmpatomname(p)(4:4) == ' ' )then
          atomname(m)(3:3) = tmpatomname(p)(1:1)
          atomname(m)(4:4) = ' '
        else
          atomname(m)(3:3) = tmpatomname(p)(4:4)
          atomname(m)(4:4) = tmpatomname(p)(1:1)
        endif
      endif
      xyz_to_pdb(m) = p
      pdb_to_xyz(p) = m
    enddo
  enddo
  if (warned) then
    do m = 1,num_pdb_atoms
      pdb_to_xyz(m) = m
      xyz_to_pdb(m) = m
      atomname(m) = adjustl(tmpatomname(m))
    end do
  end if
  ! sanity check---for double matches
  do m = 1,num_pdb_atoms
    if ( pdb_to_xyz(m) == 0 )then
      write(6,*)'no xyz atom matched to pdb atom ',m
      stop
    endif
    if ( xyz_to_pdb(m) == 0 )then
      write(6,*)'no pdb atom matched to xyz atom ',m
      stop
    endif
  enddo
  
end subroutine AM_PARM_get_pdb
!------------------------------------------------------------------
subroutine TPDB_match(x,y,z,tmpx,tmpy,tmpz,startlist,endlist,res)
  double precision,intent(in) :: x,y,z,tmpx(*),tmpy(*),tmpz(*)
  integer,intent(in) :: startlist,endlist
  integer,intent(out) :: res

  integer p
  double precision tol,dx,dy,dz,dis,mindis
  tol = 1.d-3

  res = 0
  mindis = 100.d0
  do p = startlist,endlist
    dx = tmpx(p) - x
    dy = tmpy(p) - y
    dz = tmpz(p) - z
    dis = sqrt(dx*dx+dy*dy+dz*dz)
    if ( dis < mindis )then
      res = p
      mindis = dis
    endif
  enddo
end subroutine TPDB_match
!------------------------------------------------------------------
subroutine AM_PARM_get_atom_params(analout_unit)
  integer, intent(in) :: analout_unit

  integer :: success,num
  integer :: ier1,ier2,ier3,ier4,ier5,ier6
  character(len=120) :: line
  integer :: ios,n,k
  character(len=4)atomname

  rewind(analout_unit)
  call AM_PARM_get_one_num(analout_unit,'Atoms in System',num,success)
  if ( success > 0 )then
    num_atoms = num
    write(6,*)'num_atoms = ',num_atoms
  else
    write(6,*)'no success finding Atoms in System!'
    stop
  endif
  if ( num_atoms /= num_xyz_atoms )then
    write(6,*)'number mismatch between analout and xyz files!!'
    stop
  endif
  ! allocate
  allocate(atomtype(num_atoms),stat=ier1)
  allocate(atomclass(num_atoms),stat=ier2)
  allocate(atomic_number(num_atoms),stat=ier3)
  allocate(atomvalence(num_atoms),stat=ier4)
  allocate(atomic_weight(num_atoms),stat=ier5)
  allocate(atomsymbol(num_atoms),stat=ier6)
  if ( (ier1 /= 0) .or. (ier2 /= 0) .or. (ier3 /= 0) .or.  &
       (ier4 /= 0) .or. (ier5 /= 0) .or. (ier6 /= 0)) then
    write(6,*)'atomtype allocation prob!!'
    stop
  endif
  call AM_PARM_get_section(analout_unit, &
           'Atom Type Definition Parameters',success)
  if ( success == 0 )then
    write(6,*)'no success getting Atom Type Definition Parameters'
    stop
  endif
  ! read 3 lines
  read(analout_unit,'(A)',iostat=ios)line
  read(analout_unit,'(A)',iostat=ios)line
  read(analout_unit,'(A)',iostat=ios)line
  do n = 1,num_atoms
    read(analout_unit,*)k,atomsymbol(k),atomtype(k),atomclass(k), &
            atomic_number(k),atomic_weight(k),atomvalence(k)
  enddo
  
end subroutine AM_PARM_get_atom_params
!------------------------------------------------------------------
subroutine AM_PARM_get_atom_vdw_params(analout_unit)
  integer, intent(in) :: analout_unit

  double precision,allocatable :: atom_vdw(:,:),atom_newvdw(:,:)
  integer :: i,j,k,n,ier1,ier2,ier3,ier4,ios,success
  character(len=120) :: line
  integer :: dim1
  double precision :: wt,radi,radj,epsi,epsj
  

  allocate(atom_vdw(2,num_atoms),stat=ier1)
  allocate(vdw_parent_weight(num_atoms),stat=ier2)
  allocate(atom_newvdw(2,num_atoms),stat=ier3)
  allocate(atom_vdw_type(num_atoms),stat=ier4)
  if ( (ier1 /= 0) .or. (ier2 /= 0) .or. (ier3 /= 0) .or. (ier4 /= 0) ) then
    write(6,*)'vdw allocation prob!!'
    stop
  endif
  rewind(analout_unit)
  call AM_PARM_get_section(analout_unit, &
           'Van der Waals Parameters',success)
  if ( success == 0 )then
    write(6,*)'no success getting Van der Waals Parameters'
    stop
  endif
  ! read 3 lines
  read(analout_unit,'(A)',iostat=ios)line
  read(analout_unit,'(A)',iostat=ios)line
  read(analout_unit,'(A)',iostat=ios)line
  do n = 1,num_atoms
     read(analout_unit,'(A)',iostat=ios)line
     if ( ios < 0 )then
        write(6,*)'unexpected end of file in vdw!'
        stop
     endif
     read(line,*,iostat=ios)j,k,atom_vdw(1,n),atom_vdw(2,n),wt
     if ( ios /= 0 )then !missing last value
        wt = 0.d0
     endif
     vdw_parent_weight(n) = 1.d0 - wt
  enddo
  dim1 = 2
  call AM_PARM_compress_params(dim1,num_atoms,atom_vdw, &
                      num_vdw_types,atom_newvdw,atom_vdw_type)
  write(6,*)'vdw: num_vdw_types = ',num_vdw_types
  allocate(vdw_mix_radius(num_vdw_types,num_vdw_types),stat=ier1)
  allocate(vdw_mix_epsilon(num_vdw_types,num_vdw_types),stat=ier2)
  if ( (ier1 /= 0) .or. (ier2 /= 0) ) then
    write(6,*)'vdw allocation prob!!'
    stop
  endif
  do j = 1,num_vdw_types
    radj = atom_newvdw(1,j)
    epsj = atom_newvdw(2,j)
    do i = 1,num_vdw_types
       radi = atom_newvdw(1,i)
       epsi = atom_newvdw(2,i)
       vdw_mix_radius(i,j) = 2.d0*(radi**3+radj**3)/(radi**2+radj**2)
       vdw_mix_epsilon(i,j) = 4.d0*epsi*epsj / (sqrt(epsi)+sqrt(epsj))**2
    enddo
  enddo
  deallocate(atom_vdw)
  deallocate(atom_newvdw)
  
end subroutine AM_PARM_get_atom_vdw_params
!------------------------------------------------------------------
subroutine AM_PARM_get_bond_params(analout_unit)
  integer, intent(in) :: analout_unit

  integer :: success,num,ier1,ier2,ier3,ier4
  character(len=120) :: line
  integer :: ios,n,j,k,dim1
  integer, allocatable :: bond_parm_ptr(:)
  double precision, allocatable ::old_param(:,:) 

  rewind(analout_unit)
  call AM_PARM_get_one_num(analout_unit,'Bond Stretches',num,success)
  if ( success > 0 )then
    num_bonds = num
    write(6,*)'num_bonds = ',num_bonds
  else
    write(6,*)'no success finding Bond Stretches!'
    num_bonds = 0
    return
  endif
  ! allocate
  allocate(bond_list(3,num_bonds),stat=ier1)
  allocate(bond_parm_ptr(num_bonds),stat=ier2)
  allocate(bond_param(2,num_bonds),stat=ier3)
  allocate(old_param(2,num_bonds),stat=ier4)
  if ( (ier1 /= 0) .or. (ier2 /= 0) .or. (ier3 /= 0) .or. (ier4 /= 0) ) then
    write(6,*)'bond allocation prob!!'
    stop
  endif
  call AM_PARM_get_section(analout_unit, &
           'Bond Stretching Parameters',success)
  if ( success == 0 )then
    write(6,*)'no success getting bond Stretching Parameters'
    stop
  endif
  ! read 3 lines
  read(analout_unit,'(A)',iostat=ios)line
  read(analout_unit,'(A)',iostat=ios)line
  read(analout_unit,'(A)',iostat=ios)line
  do n = 1,num_bonds
    read(analout_unit,*)k,(bond_list(j,k),j=1,2), &
                        (old_param(j,k),j=1,2)
  enddo
  dim1 = 2
  call AM_PARM_compress_params(dim1,num_bonds,old_param, &
                      num_bond_params,bond_param,bond_parm_ptr)
  do n = 1,num_bonds
    bond_list(3,n) = bond_parm_ptr(n)
  enddo
  ! deallocate
  if ( allocated(bond_parm_ptr) )deallocate(bond_parm_ptr)
  if ( allocated(old_param) )deallocate(old_param)
end subroutine AM_PARM_get_bond_params
!------------------------------------------------------------------
subroutine AM_PARM_write_bond_params(prmtop_unit)
  integer, intent(in) :: prmtop_unit

  ! functable:
  integer :: degree = 4
  character(len=8)word,word1
  character(len=16)word2
  double precision :: coeff(0:4) = (/0.d0,0.d0,1.d0,-2.55d0,3.793125d0/)
  integer :: j,k

  if ( num_bonds == 0 )return !didn't find any
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_REGULAR_BOND_NUM_LIST'
  write(prmtop_unit,'(a)')'%FORMAT(I8)'
  write(prmtop_unit,'(I8)')num_bonds
  write(word,'(I8)')num_bonds
  word1 = adjustl(word)
  word2 = '(3,'//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_REGULAR_BOND_LIST'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')'%FORMAT(10I8)'
  write(prmtop_unit,'(10I8)')((bond_list(j,k),j=1,3),k=1,num_bonds)
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_REGULAR_BOND_NUM_PARAMS'
  write(prmtop_unit,'(a)')'%FORMAT(I8)'
  write(prmtop_unit,'(I8)')num_bond_params
  write(word,'(I8)')num_bond_params
  word1 = adjustl(word)
  word2 = '('//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_REGULAR_BOND_FORCE_CONSTANT'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat)(bond_param(1,k),k=1,num_bond_params)
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_REGULAR_BOND_EQUIL_VALUE'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat)(bond_param(2,k),k=1,num_bond_params)
  ! bond ftab nex
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_REGULAR_BOND_FTAB_DEGREE'
  write(prmtop_unit,'(a)')'%FORMAT(I8)'
  write(prmtop_unit,'(2I8)')degree
  write(word,'(I8)')degree
  word1 = adjustl(word)
  word2 = '(0:'//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_REGULAR_BOND_FTAB_COEFFS'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat)(coeff(j),j=0,4)
end subroutine AM_PARM_write_bond_params
!------------------------------------------------------------------
subroutine AM_PARM_get_UreyB_params(analout_unit)
  integer, intent(in) :: analout_unit

  integer :: success,num,ier1,ier2,ier3,ier4
  character(len=120) :: line
  integer :: ios,n,j,k,dim1,idum,ntok
  integer, allocatable :: urey_bradley_parm_ptr(:)
  double precision, allocatable ::old_param(:,:)

  rewind(analout_unit)
  call AM_PARM_get_one_num(analout_unit,'Urey-Bradley',num,success)
  if ( success > 0 )then
    num_urey_bradley = num
    write(6,*)'num_urey_bradley = ',num_urey_bradley
  else
    write(6,*)'no success finding Urey-Bradley!'
    num_urey_bradley = 0
    return
  endif
  ! allocate
  allocate(urey_bradley_list(3,num_urey_bradley),stat=ier1)
  allocate(urey_bradley_parm_ptr(num_urey_bradley),stat=ier2)
  allocate(urey_bradley_param(2,num_urey_bradley),stat=ier3)
  allocate(old_param(2,num_urey_bradley),stat=ier4)
  if ( (ier1 /= 0) .or. (ier2 /= 0) .or. (ier3 /= 0) .or. (ier4 /= 0) ) then
    write(6,*)'UreyBradley allocation prob!!'
    stop
  endif
  call AM_PARM_get_section(analout_unit, &
           'Urey-Bradley Parameters',success)
  if ( success == 0 )then
    write(6,*)'no success getting bond Stretching Parameters'
    stop
  endif
  ! read 3 lines
  read(analout_unit,'(A)',iostat=ios)line
  read(analout_unit,'(A)',iostat=ios)line
  read(analout_unit,'(A)',iostat=ios)line
  do n = 1,num_urey_bradley
    read(analout_unit,'(a)') line
    call get_num_tokens(line, ntok)
    if (ntok == 6) then
      read(line, *) k,urey_bradley_list(1,k),idum,urey_bradley_list(2,k), &
                    (old_param(j,k),j=1,2)
    else if (ntok == 5) then
      read(line, *) k,urey_bradley_list(1,k),urey_bradley_list(2,k), &
                    (old_param(j,k),j=1,2)
    else
      write(0,*) 'Could not parse [[', trim(line), ']] in Urey-Bradley terms.'
      stop
    end if
  enddo
  dim1 = 2
  call AM_PARM_compress_params(dim1,num_urey_bradley,old_param, &
                      num_ureyb_params,urey_bradley_param, &
                      urey_bradley_parm_ptr)
  do n = 1,num_urey_bradley
    urey_bradley_list(3,n) = urey_bradley_parm_ptr(n)
  enddo
  ! deallocate
  if ( allocated(urey_bradley_parm_ptr) )deallocate(urey_bradley_parm_ptr)
  if ( allocated(old_param) )deallocate(old_param)
end subroutine AM_PARM_get_UreyB_params
!------------------------------------------------------------------
subroutine AM_PARM_write_UreyB_params(prmtop_unit)
  integer,intent(in) :: prmtop_unit

  ! functable:
  integer :: degree = 2
  double precision :: coeff(0:2) = (/0.d0,0.d0,1.d0/)
  character(len=8)word,word1
  character(len=16)word2
  integer j,k

  if ( num_urey_bradley == 0 )return
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_UREY_BRADLEY_BOND_NUM_LIST'
  write(prmtop_unit,'(a)')'%FORMAT(I8)'
  write(prmtop_unit,'(I8)')num_urey_bradley
  write(word,'(I8)')num_urey_bradley
  word1 = adjustl(word)
  word2 = '(3,'//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_UREY_BRADLEY_BOND_LIST' 
  write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')'%FORMAT(10I8)'
  write(prmtop_unit,'(10I8)')((urey_bradley_list(j,k),j=1,3),  &
                                  k=1,num_urey_bradley)
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_UREY_BRADLEY_BOND_NUM_PARAMS'
  write(prmtop_unit,'(a)')'%FORMAT(I8)'
  write(prmtop_unit,'(I8)')num_ureyb_params
  write(word,'(I8)')num_ureyb_params
  word1 = adjustl(word)
  word2 = '('//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_UREY_BRADLEY_BOND_FORCE_CONSTANT'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat)(urey_bradley_param(1,k),k=1,num_ureyb_params)
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_UREY_BRADLEY_BOND_EQUIL_VALUE'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat)(urey_bradley_param(2,k),k=1,num_ureyb_params)
  ! bond ftab next
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_UREY_BRADLEY_BOND_FTAB_DEGREE'
  write(prmtop_unit,'(a)')'%FORMAT(I8)'
  write(prmtop_unit,'(2I8)')degree
  write(word,'(I8)')degree
  word1 = adjustl(word)
  word2 = '(0:'//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_UREY_BRADLEY_BOND_FTAB_COEFFS'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat)(coeff(j),j=0,2)
end subroutine AM_PARM_write_UreyB_params
!------------------------------------------------------------------
subroutine AM_PARM_get_angle_params(analout_unit)
  integer, intent(in) :: analout_unit

  integer :: success,num
  character(len=120) :: line
  integer :: ios,n,j,k,dim1
  integer :: ier1,ier2,ier3,ier4,ier5,ier6,ier7,ier8
  integer, allocatable :: regular_angle_parm_ptr(:)
  double precision, allocatable ::regular_angle_old_param(:,:) 
  integer, allocatable :: trigonal_angle_parm_ptr(:)
  double precision, allocatable ::trigonal_angle_old_param(:,:) 

  rewind(analout_unit)
  call AM_PARM_get_one_num(analout_unit,'Angle Bends',num,success)
  if ( success > 0 )then
    num_tot_angles = num
    write(6,*)'num_tot_angles = ',num_tot_angles
  else
    write(6,*)'no success finding Angle Bending Parameters!'
    num_tot_angles = 0
    return
  endif
  ! allocate
  allocate(regular_angle_list(4,num_tot_angles),stat=ier1)
  allocate(regular_angle_parm_ptr(num_tot_angles),stat=ier2)
  allocate(regular_angle_param(2,num_tot_angles),stat=ier3)
  allocate(regular_angle_old_param(2,num_tot_angles),stat=ier4)
  ! need a 4th atom for trigonal!!
  allocate(trigonal_angle_list(5,num_tot_angles),stat=ier5)
  allocate(trigonal_angle_parm_ptr(num_tot_angles),stat=ier6)
  allocate(trigonal_angle_param(2,num_tot_angles),stat=ier7)
  allocate(trigonal_angle_old_param(2,num_tot_angles),stat=ier8)
  if ( (ier1 /= 0) .or. (ier2 /= 0) .or. (ier3 /= 0) .or.  &
       (ier4 /= 0) .or. (ier5 /= 0) .or. (ier6 /= 0) .or.  & 
       (ier7 /= 0) .or. (ier8 /= 0) ) then
    write(6,*)'angle allocation prob!!'
    stop
  endif
  ! read in data
  call AM_PARM_get_section(analout_unit, &
           'Angle Bending Parameters',success)
  if ( success == 0 )then
    write(6,*)'no success getting angle Bending Parameters'
    stop
  endif
  ! read 3 lines
  read(analout_unit,'(A)',iostat=ios)line
  read(analout_unit,'(A)',iostat=ios)line
  read(analout_unit,'(A)',iostat=ios)line
  num_trigonal_angles = 0
  num_regular_angles = 0
  do n = 1,num_tot_angles
    read(analout_unit,'(A)',iostat=ios)line
    if ( index(line,'In-Plane') > 0 )then
      num_trigonal_angles = num_trigonal_angles + 1
      read(line,*)k,(trigonal_angle_list(j,num_trigonal_angles),j=1,3),  &
                    (trigonal_angle_old_param(j,num_trigonal_angles),j=1,2)
      call AM_PARM_find_trig_4th_atom(  &
            trigonal_angle_list(:,num_trigonal_angles))
    else
      num_regular_angles = num_regular_angles + 1
      read(line,*)k,(regular_angle_list(j,num_regular_angles),j=1,3),  &
                    (regular_angle_old_param(j,num_regular_angles),j=1,2)
    endif
  enddo
  write(6,*)'num regular, trigonal angles = ',num_regular_angles, &
                 num_trigonal_angles
  dim1 = 2
  call AM_PARM_compress_params(dim1,num_regular_angles,  &
                      regular_angle_old_param, &
                      num_regular_angle_params,regular_angle_param,  &
                      regular_angle_parm_ptr)
  do n = 1,num_regular_angles
    regular_angle_list(4,n) = regular_angle_parm_ptr(n)
  enddo
  call AM_PARM_compress_params(dim1,num_trigonal_angles,  &
                      trigonal_angle_old_param, &
                      num_trigonal_angle_params,trigonal_angle_param,  &
                      trigonal_angle_parm_ptr)
  do n = 1,num_trigonal_angles
    trigonal_angle_list(5,n) = trigonal_angle_parm_ptr(n)
  enddo
  ! deallocate
  if ( allocated(regular_angle_parm_ptr) )deallocate(regular_angle_parm_ptr)
  if ( allocated(regular_angle_old_param) )deallocate(regular_angle_old_param)
  if ( allocated(trigonal_angle_parm_ptr) )deallocate(trigonal_angle_parm_ptr)
  if ( allocated(trigonal_angle_old_param) )deallocate(trigonal_angle_old_param)
end subroutine AM_PARM_get_angle_params
!------------------------------------------------------------------
subroutine AM_PARM_write_angle_params(prmtop_unit)
  integer, intent(in) :: prmtop_unit

  ! functable:
  integer :: degree = 6
  double precision :: coeff(0:6) =   &
          (/0.d0,0.d0,1.d0,-0.014d0, 0.000056d0,-0.0000007d0,0.000000022d0/)
  character(len=8)word,word1
  character(len=16)word2
  integer j,k
  
  if ( num_tot_angles == 0 )return
  ! dump out REGULAR angles
  if ( num_regular_angles > 0 )then
    write(prmtop_unit,'(a)')'%FLAG AMOEBA_REGULAR_ANGLE_NUM_LIST'
    write(prmtop_unit,'(a)')'%FORMAT(I8)'
    write(prmtop_unit,'(I8)')num_regular_angles
    write(word,'(I8)')num_regular_angles
    word1 = adjustl(word)
    word2 = '(4,'//word1(1:len_trim(word1))//')'
    write(prmtop_unit,'(a)')'%FLAG AMOEBA_REGULAR_ANGLE_LIST'
    write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
    write(prmtop_unit,'(a)')'%FORMAT(10I8)'
    write(prmtop_unit,'(10I8)')((regular_angle_list(j,k),j=1,4),  &
                               k=1,num_regular_angles)
    write(prmtop_unit,'(a)')'%FLAG AMOEBA_REGULAR_ANGLE_NUM_PARAMS'
    write(prmtop_unit,'(a)')'%FORMAT(I8)'
    write(prmtop_unit,'(I8)')num_regular_angle_params
    write(word,'(I8)')num_regular_angle_params
    word1 = adjustl(word)
    word2 = '('//word1(1:len_trim(word1))//')'
    write(prmtop_unit,'(a)')'%FLAG AMOEBA_REGULAR_ANGLE_FORCE_CONSTANT'
    write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
    write(prmtop_unit,'(a)')outrealformat
    write(prmtop_unit,realformat)(regular_angle_param(1,k),  &
                                  k=1,num_regular_angle_params)
    write(prmtop_unit,'(a)')'%FLAG AMOEBA_REGULAR_ANGLE_EQUIL_VALUE'
    write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
    write(prmtop_unit,'(a)')outrealformat
    write(prmtop_unit,realformat)(regular_angle_param(2,k),  &
                                 k=1,num_regular_angle_params)
  ! angle ftab next
    write(prmtop_unit,'(a)')'%FLAG AMOEBA_REGULAR_ANGLE_FTAB_DEGREE'
    write(prmtop_unit,'(a)')'%FORMAT(I8)'
    write(prmtop_unit,'(2I8)')degree
    write(word,'(I8)')degree
    word1 = adjustl(word)
    word2 = '(0:'//word1(1:len_trim(word1))//')'
    write(prmtop_unit,'(a)')'%FLAG AMOEBA_REGULAR_ANGLE_FTAB_COEFFS'
    write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
    write(prmtop_unit,'(a)')outrealformat
    write(prmtop_unit,realformat)(coeff(j),j=0,6)
  endif ! num_regular_angles>0
  ! dump out TRIGONAL angles
  if ( num_trigonal_angles > 0 )then
    write(prmtop_unit,'(a)')'%FLAG AMOEBA_TRIGONAL_ANGLE_NUM_LIST'
    write(prmtop_unit,'(a)')'%FORMAT(I8)'
    write(prmtop_unit,'(I8)')num_trigonal_angles
    write(word,'(I8)')num_trigonal_angles
    word1 = adjustl(word)
    word2 = '(5,'//word1(1:len_trim(word1))//')'
    write(prmtop_unit,'(a)')'%FLAG AMOEBA_TRIGONAL_ANGLE_LIST'
    write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
    write(prmtop_unit,'(a)')'%FORMAT(10I8)'
    write(prmtop_unit,'(10I8)')((trigonal_angle_list(j,k),j=1,5),  &
                               k=1,num_trigonal_angles)
    write(prmtop_unit,'(a)')'%FLAG AMOEBA_TRIGONAL_ANGLE_NUM_PARAMS'
    write(prmtop_unit,'(a)')'%FORMAT(I8)'
    write(prmtop_unit,'(I8)')num_trigonal_angle_params
    write(word,'(I8)')num_trigonal_angle_params
    word1 = adjustl(word)
    word2 = '('//word1(1:len_trim(word1))//')'
    write(prmtop_unit,'(a)')'%FLAG AMOEBA_TRIGONAL_ANGLE_FORCE_CONSTANT'
    write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
    write(prmtop_unit,'(a)')outrealformat
    write(prmtop_unit,realformat)(trigonal_angle_param(1,k),  &
                                 k=1,num_trigonal_angle_params)
    write(prmtop_unit,'(a)')'%FLAG AMOEBA_TRIGONAL_ANGLE_EQUIL_VALUE'
    write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
    write(prmtop_unit,'(a)')outrealformat
    write(prmtop_unit,realformat)(trigonal_angle_param(2,k),  &
                                 k=1,num_trigonal_angle_params)
    ! angle ftab next
    write(prmtop_unit,'(a)')'%FLAG AMOEBA_TRIGONAL_ANGLE_FTAB_DEGREE'
    write(prmtop_unit,'(a)')'%FORMAT(I8)'
    write(prmtop_unit,'(2I8)')degree
    write(word,'(I8)')degree
    word1 = adjustl(word)
    word2 = '(0:'//word1(1:len_trim(word1))//')'
    write(prmtop_unit,'(a)')'%FLAG AMOEBA_TRIGONAL_ANGLE_FTAB_COEFFS'
    write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
    write(prmtop_unit,'(a)')outrealformat
    write(prmtop_unit,realformat)(coeff(j),j=0,6)
  endif ! num_trigonal_angles > 0
end subroutine AM_PARM_write_angle_params
!------------------------------------------------------------------
subroutine AM_PARM_get_outofplane_params(analout_unit)
  integer, intent(in) :: analout_unit

  integer, allocatable :: outofplane_parm_ptr(:)
  double precision, allocatable ::old_param(:)

  integer :: success,num,dim1,ier1,ier2,ier3,ier4
  character(len=120) :: line
  integer :: ios,n,j,k

  rewind(analout_unit)
  call AM_PARM_get_one_num(analout_unit,'Out-of-Plane Bends',num,success)
  if ( success > 0 )then
    num_out_of_plane_bends = num
    write(6,*)'num_out_of_plane_bends = ',num_out_of_plane_bends
  else
    write(6,*)'no success finding Out-of-Plane Bends'
    num_out_of_plane_bends = 0
    return
  endif
  ! allocate
  allocate(outofplane_list(5,num_out_of_plane_bends),stat=ier1)
  allocate(outofplane_parm_ptr(num_out_of_plane_bends),stat=ier2)
  allocate(outofplane_param(num_out_of_plane_bends),stat=ier3)
  allocate(old_param(num_out_of_plane_bends),stat=ier4)
  if ( (ier1 /= 0) .or. (ier2 /= 0) .or. (ier3 /= 0) .or. (ier4 /= 0) ) then
    write(6,*)'Out-of-Plane allocation prob!!'
    stop
  endif
  call AM_PARM_get_section(analout_unit, &
           'Out-of-Plane Bending Parameters',success)
  if ( success == 0 )then
    write(6,*)'no success getting Out-of-Plane Parameters'
    stop
  endif
  ! read 3 lines
  read(analout_unit,'(A)',iostat=ios)line
  read(analout_unit,'(A)',iostat=ios)line
  read(analout_unit,'(A)',iostat=ios)line
  do n = 1,num_out_of_plane_bends
    ! note the scrambled order---amoeba_valence code treats this atom list
    ! the same as the trigonal angles---should match up 1-to-1
    ! thus the 4th atom is the "wagging" one, while atom 2 is angle center
    ! which gets projected
    read(analout_unit,*)k,outofplane_list(4,k),outofplane_list(2,k), &
                        outofplane_list(1,k),outofplane_list(3,k), &
                        old_param(k)
  enddo
  dim1 = 1
  call AM_PARM_compress_params(dim1,num_out_of_plane_bends,  &
                      old_param, &
                      num_opbend_params,outofplane_param,  &
                      outofplane_parm_ptr)
  do n = 1,num_out_of_plane_bends
    outofplane_list(5,n) = outofplane_parm_ptr(n)
  enddo
  if ( allocated(outofplane_parm_ptr) )deallocate(outofplane_parm_ptr)
  if ( allocated(old_param) )deallocate(old_param)
end subroutine AM_PARM_get_outofplane_params
!------------------------------------------------------------------
subroutine AM_PARM_write_outofplane_params(prmtop_unit)
  integer, intent(in) :: prmtop_unit
  character(len=8)word,word1
  character(len=16)word2
  ! functable:
  integer :: degree = 6
  double precision :: coeff(0:6) =   &
          (/0.d0,0.d0,1.d0,-0.014d0, 0.000056d0,-0.0000007d0,0.000000022d0/)
  integer j,k
  ! dump out angles
  if ( num_out_of_plane_bends == 0 )return
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_OPBEND_ANGLE_NUM_LIST'
  write(prmtop_unit,'(a)')'%FORMAT(I8)'
  write(prmtop_unit,'(I8)')num_out_of_plane_bends
  write(word,'(I8)')num_out_of_plane_bends
  word1 = adjustl(word)
  word2 = '(5,'//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_OPBEND_ANGLE_LIST'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')'%FORMAT(10I8)'
  write(prmtop_unit,'(10I8)')((outofplane_list(j,k),j=1,5),  &
                               k=1,num_out_of_plane_bends)
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_OPBEND_ANGLE_NUM_PARAMS'
  write(prmtop_unit,'(a)')'%FORMAT(I8)'
  write(prmtop_unit,'(I8)')num_opbend_params
  write(word,'(I8)')num_opbend_params
  word1 = adjustl(word)
  word2 = '('//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_OPBEND_ANGLE_FORCE_CONSTANT'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat)(outofplane_param(k),k=1,num_opbend_params)
  ! angle ftab next
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_OPBEND_ANGLE_FTAB_DEGREE'
  write(prmtop_unit,'(a)')'%FORMAT(I8)'
  write(prmtop_unit,'(I8)')degree
  write(word,'(I8)')degree
  word1 = adjustl(word)
  word2 = '(0:'//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_OPBEND_ANGLE_FTAB_COEFFS'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat)(coeff(j),j=0,6)
end subroutine AM_PARM_write_outofplane_params
!------------------------------------------------------------------
subroutine AM_PARM_get_torsion_params(analout_unit)
  integer, intent(in) :: analout_unit

  integer, allocatable :: torsion_parm_ptr(:)
  double precision, allocatable ::old_param(:,:)
  integer :: success,num,numtor,numterms,num1,ier1,ier2,ier3,ier4
  character(len=120) :: line
  integer :: ios,n,j,k,i1,i2,i3,i4,phase(10),period(10)
  integer start,find,temp,dim1
  double precision :: coeff(10),theta

  rewind(analout_unit)
  call AM_PARM_get_one_num(analout_unit,'Torsional Angles',num,success)
  if ( success > 0 )then
    numtor = num
    write(6,*)'numtor = ',numtor
  else
    write(6,*)'no success finding Torsional Angles'
    num_torsions = 0
    return
  endif
  call AM_PARM_get_section(analout_unit, &
           'Torsional Angle Parameters',success)
  if ( success == 0 )then
    write(6,*)'no success getting Torsional Angle Parameters'
    stop
  endif
  ! read 3 lines
  read(analout_unit,'(A)',iostat=ios)line
  read(analout_unit,'(A)',iostat=ios)line
  read(analout_unit,'(A)',iostat=ios)line
  num_torsions = 0
  num1 = 0
  do n = 1,numtor
    read(analout_unit,'(A)',iostat=ios)line
    call AM_PARM_get_numtorterms(line,numterms)
    if ( numterms > 0 )num1 = num1 + 1
    num_torsions = num_torsions + numterms
  enddo
  write(6,*)'num1,num_torsions = ',num1,num_torsions
  ! allocate
  allocate(torsion_list(5,num_torsions),stat=ier1)! 4 atoms plus parm ptr
  allocate(torsion_parm_ptr(num_torsions),stat=ier2)
  allocate(torsion_param(3,num_torsions),stat=ier3)
  allocate(old_param(3,num_torsions),stat=ier4)
  if ( (ier1 /= 0) .or. (ier2 /= 0) .or. (ier3 /= 0) .or. (ier4 /= 0) ) then
    write(6,*)'torsion allocation prob!!'
    stop
  endif
  rewind(analout_unit)
  call AM_PARM_get_section(analout_unit, &
           'Torsional Angle Parameters',success)
  if ( success == 0 )then
    write(6,*)'no success getting Torsional Angle Parameters'
    stop
  endif
  ! read 3 lines
  read(analout_unit,'(A)',iostat=ios)line
  read(analout_unit,'(A)',iostat=ios)line
  read(analout_unit,'(A)',iostat=ios)line
  num_torsions = 0
  do n = 1,numtor
    read(analout_unit,'(A)',iostat=ios)line
    call AM_PARM_get_numtorterms(line,numterms)
    if ( numterms > 0 )then
      start = 1
      find = index(line(start:),"/")
      read(line(start:find-1),*)j,i1,i2,i3,i4,coeff(1),phase(1)
      read(line(find+1:),*)period(1)
      num_torsions = num_torsions + 1
      torsion_list(1,num_torsions) = i1
      torsion_list(2,num_torsions) = i2
      torsion_list(3,num_torsions) = i3
      torsion_list(4,num_torsions) = i4
      old_param(1,num_torsions) = coeff(1)
      old_param(2,num_torsions) = dble(period(1))
      theta = degrees_to_radians*phase(1)
      old_param(3,num_torsions) = theta
      do k = 2,numterms
        start = start + find
        find = index(line(start:),"/")
        read(line(start:start+find-2),*)temp,coeff(k),phase(k)
        read(line(start+find:),*)period(k)
        num_torsions = num_torsions + 1
        torsion_list(1,num_torsions) = i1
        torsion_list(2,num_torsions) = i2
        torsion_list(3,num_torsions) = i3
        torsion_list(4,num_torsions) = i4
        old_param(1,num_torsions) = coeff(k)
        old_param(2,num_torsions) = dble(period(k))
        theta = degrees_to_radians*phase(k)
        old_param(3,num_torsions) = theta
      enddo
    endif
  enddo 
  dim1 = 3
  call AM_PARM_compress_params(dim1,num_torsions,  &
                      old_param, &
                      num_torsion_params,torsion_param,  &
                      torsion_parm_ptr)
  do n = 1,num_torsions
    torsion_list(5,n) = torsion_parm_ptr(n)
  enddo
  if ( allocated(torsion_parm_ptr) )deallocate(torsion_parm_ptr)
  if ( allocated(old_param) )deallocate(old_param)
end subroutine AM_PARM_get_torsion_params
!------------------------------------------------------------------
subroutine AM_PARM_write_torsion_params(prmtop_unit)
  integer, intent(in) :: prmtop_unit
  character(len=8)word,word1
  character(len=16)word2

  integer j,k
  ! dump out torsions
  if ( num_torsions == 0 )return
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_TORSION_NUM_LIST'
  write(prmtop_unit,'(a)')'%FORMAT(I8)'
  write(prmtop_unit,'(I8)')num_torsions
  write(word,'(I8)')num_torsions
  word1 = adjustl(word)
  word2 = '(5,'//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_TORSION_LIST'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')'%FORMAT(10I8)'
  write(prmtop_unit,'(10I8)')((torsion_list(j,k),j=1,5),  &
                               k=1,num_torsions)
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_TORSION_NUM_PARAMS'
  write(prmtop_unit,'(a)')'%FORMAT(I8)'
  write(prmtop_unit,'(I8)')num_torsion_params
  write(word,'(I8)')num_torsion_params
  word1 = adjustl(word)
  word2 = '('//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_TORSION_FORCE_CONSTANT'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat)(torsion_param(1,k),k=1,num_torsion_params)
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_TORSION_PERIODICITY'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat)(torsion_param(2,k),k=1,num_torsion_params)
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_TORSION_PHASE'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat)(torsion_param(3,k),k=1,num_torsion_params)
end subroutine AM_PARM_write_torsion_params
!------------------------------------------------------------------
subroutine AM_PARM_get_stretch_tor_params(analout_unit)
  integer, intent(in) :: analout_unit

  integer, allocatable :: str_torsion_parm_ptr(:)
  double precision, allocatable ::old_param(:,:)
  integer :: success,num,numstrtor,numterms,num1,ier1,ier2,ier3,ier4
  character(len=120) :: line
  integer :: ios,n,j,k,i1,i2,i3,i4,phase(10),period(10)
  integer start,find,temp,dim1
  double precision :: coeff(10),theta,bondlength

  rewind(analout_unit)
  call AM_PARM_get_one_num(analout_unit,'Stretch-Torsions',num,success)
  if ( success > 0 )then
    numstrtor = num
    write(6,*)'numstrtor = ',numstrtor
  else
    write(6,*)'no success finding Stretch-Torsions'
    num_str_torsions = 0
    return
  endif
  call AM_PARM_get_section(analout_unit, &
           'Stretch-Torsion Parameters',success)
  if ( success == 0 )then
    write(6,*)'no success getting Stretch-Torsion Parameters'
    stop
  endif
  ! read 3 lines
  read(analout_unit,'(A)',iostat=ios)line
  read(analout_unit,'(A)',iostat=ios)line
  read(analout_unit,'(A)',iostat=ios)line
  num_str_torsions = 0
  num1 = 0
  do n = 1,numstrtor
    read(analout_unit,'(A)',iostat=ios)line
    call AM_PARM_get_numtorterms(line,numterms)
    if ( numterms > 0 )num1 = num1 + 1
    num_str_torsions = num_str_torsions + numterms
  enddo
  write(6,*)'num1,num_str_torsions = ',num1,num_str_torsions
  ! allocate
  allocate(str_torsion_list(5,num_str_torsions),stat=ier1)! 4 atoms + parm ptr
  allocate(str_torsion_parm_ptr(num_str_torsions),stat=ier2)
  allocate(str_torsion_param(4,num_str_torsions),stat=ier3)
  allocate(old_param(4,num_str_torsions),stat=ier4)
  if ( (ier1 /= 0) .or. (ier2 /= 0) .or. (ier3 /= 0) .or. (ier4 /= 0) ) then
    write(6,*)'torsion allocation prob!!'
    stop
  endif
  rewind(analout_unit)
  call AM_PARM_get_section(analout_unit, &
           'Stretch-Torsion Parameters',success)
  if ( success == 0 )then
    write(6,*)'no success getting Stretch-Torsion Parameters'
    stop
  endif
  ! read 3 lines
  read(analout_unit,'(A)',iostat=ios)line
  read(analout_unit,'(A)',iostat=ios)line
  read(analout_unit,'(A)',iostat=ios)line
  num_str_torsions = 0
  do n = 1,numstrtor
    read(analout_unit,'(A)',iostat=ios)line
    call AM_PARM_get_numtorterms(line,numterms)
    if ( numterms > 0 )then
      start = 1
      find = index(line(start:),"/")
      read(line(start:find-1),*)j,i1,i2,i3,i4,bondlength,coeff(1),phase(1)
      read(line(find+1:),*)period(1)
      num_str_torsions = num_str_torsions + 1
      str_torsion_list(1,num_str_torsions) = i1
      str_torsion_list(2,num_str_torsions) = i2
      str_torsion_list(3,num_str_torsions) = i3
      str_torsion_list(4,num_str_torsions) = i4
      old_param(1,num_str_torsions) = coeff(1)
      old_param(2,num_str_torsions) = dble(period(1))
      theta = degrees_to_radians*phase(1)
      old_param(3,num_str_torsions) = theta
      old_param(4,num_str_torsions) = bondlength
      do k = 2,numterms
        start = start + find
        find = index(line(start:),"/")
        read(line(start:start+find-2),*)temp,coeff(k),phase(k)
        read(line(start+find:),*)period(k)
        num_str_torsions = num_str_torsions + 1
        str_torsion_list(1,num_str_torsions) = i1
        str_torsion_list(2,num_str_torsions) = i2
        str_torsion_list(3,num_str_torsions) = i3
        str_torsion_list(4,num_str_torsions) = i4
        old_param(1,num_str_torsions) = coeff(k)
        old_param(2,num_str_torsions) = dble(period(k))
        theta = degrees_to_radians*phase(k)
        old_param(3,num_str_torsions) = theta
        old_param(4,num_str_torsions) = bondlength
      enddo
    endif
  enddo 
  dim1 = 4
  call AM_PARM_compress_params(dim1,num_str_torsions,  &
                      old_param, &
                      num_str_torsion_params,str_torsion_param,  &
                      str_torsion_parm_ptr)
  do n = 1,num_str_torsions
    str_torsion_list(5,n) = str_torsion_parm_ptr(n)
  enddo
  if ( allocated(str_torsion_parm_ptr) )deallocate(str_torsion_parm_ptr)
  if ( allocated(old_param) )deallocate(old_param)
end subroutine AM_PARM_get_stretch_tor_params
!------------------------------------------------------------------
subroutine AM_PARM_write_str_tor_params(prmtop_unit)
  integer, intent(in) :: prmtop_unit
  character(len=8)word,word1
  character(len=16)word2

  integer j,k
  ! dump out torsions
  if ( num_str_torsions == 0 )return
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_STRETCH_TORSION_NUM_LIST'
  write(prmtop_unit,'(a)')'%FORMAT(I8)'
  write(prmtop_unit,'(I8)')num_str_torsions
  write(word,'(I8)')num_str_torsions
  word1 = adjustl(word)
  word2 = '(5,'//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_STRETCH_TORSION_LIST'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')'%FORMAT(10I8)'
  write(prmtop_unit,'(10I8)')((str_torsion_list(j,k),j=1,5),  &
                               k=1,num_str_torsions)
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_STRETCH_TORSION_NUM_PARAMS'
  write(prmtop_unit,'(a)')'%FORMAT(I8)'
  write(prmtop_unit,'(I8)')num_str_torsion_params
  write(word,'(I8)')num_str_torsion_params
  word1 = adjustl(word)
  word2 = '('//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_STRETCH_TORSION_FORCE_CONSTANT'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat) &
          (str_torsion_param(1,k),k=1,num_str_torsion_params)
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_STRETCH_TORSION_PERIODICITY'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat) &
          (str_torsion_param(2,k),k=1,num_str_torsion_params)
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_STRETCH_TORSION_PHASE'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat) &
          (str_torsion_param(3,k),k=1,num_str_torsion_params)
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_STRETCH_TORSION_BOND_EQUIL_VALUE'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat) &
          (str_torsion_param(4,k),k=1,num_str_torsion_params)
end subroutine AM_PARM_write_str_tor_params
!------------------------------------------------------------------
subroutine AM_PARM_get_pitorsion_params(analout_unit)
  integer, intent(in) :: analout_unit

  integer, allocatable :: pitorsion_parm_ptr(:)
  double precision, allocatable :: old_param(:,:)

  integer :: success,num,dim1,ier1,ier2,ier3,ier4
  character(len=120) :: line
  integer :: ios,n,j,k

  rewind(analout_unit)
  call AM_PARM_get_one_num(analout_unit,'Pi-Orbital Torsions',num,success)
  if ( success > 0 )then
    num_pitorsions = num
    write(6,*)'num_pitorsions = ',num_pitorsions
  else
    write(6,*)'no success finding Pi-Orbital Torsions'
    num_pitorsions = 0
    return
  endif
  ! allocate
  allocate(pitorsion_list(7,num_pitorsions),stat=ier1)! 6 atoms plus parm ptr
  allocate(pitorsion_parm_ptr(num_pitorsions),stat=ier2)
  allocate(pitorsion_param(3,num_pitorsions),stat=ier3)
  allocate(old_param(3,num_pitorsions),stat=ier4)
  if ( (ier1 /= 0) .or. (ier2 /= 0) .or. (ier3 /= 0) .or. (ier4 /= 0) ) then
    write(6,*)'PiTorsion allocation prob!!'
    stop
  endif
  call AM_PARM_get_section(analout_unit, &
           'Pi-Orbital Torsion Parameters',success)
  if ( success == 0 )then
    write(6,*)'no success getting Pi-Orbital Torsion Parameters'
    stop
  endif
  ! read 3 lines
  read(analout_unit,'(A)',iostat=ios)line
  read(analout_unit,'(A)',iostat=ios)line
  read(analout_unit,'(A)',iostat=ios)line
  do n = 1,num_pitorsions
    ! note just central pair listed--fill out rest ourselves
    read(analout_unit,*)k,pitorsion_list(3,k),pitorsion_list(4,k), &
                        old_param(1,k)
    call AM_PARM_fillout_pitorsion_6(pitorsion_list(:,k))
    old_param(2,k) = 2.d0 ! period is 2
    old_param(3,k) = pi ! phase is 180
  enddo
  dim1 = 3
  call AM_PARM_compress_params(dim1,num_pitorsions,  &
                      old_param, &
                      num_pitorsion_params,pitorsion_param,  &
                      pitorsion_parm_ptr)
  do n = 1,num_pitorsions
    pitorsion_list(7,n) = pitorsion_parm_ptr(n)
  enddo
  if ( allocated(pitorsion_parm_ptr) )deallocate(pitorsion_parm_ptr)
  if ( allocated(old_param) )deallocate(old_param)
end subroutine AM_PARM_get_pitorsion_params
!------------------------------------------------------------------
subroutine AM_PARM_write_pitorsion_params(prmtop_unit)
  integer, intent(in) :: prmtop_unit

  integer j,k
  character(len=8)word,word1
  character(len=16)word2
  ! dump out pitorsions
  if ( num_pitorsions == 0 )return
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_PI_TORSION_NUM_LIST'
  write(prmtop_unit,'(a)')'%FORMAT(I8)'
  write(prmtop_unit,'(I8)')num_pitorsions
  write(word,'(I8)')num_pitorsions
  word1 = adjustl(word)
  word2 = '(7,'//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_PI_TORSION_LIST'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')'%FORMAT(10I8)'
  write(prmtop_unit,'(10I8)')((pitorsion_list(j,k),j=1,7),  &
                               k=1,num_pitorsions)
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_PI_TORSION_NUM_PARAMS'
  write(prmtop_unit,'(a)')'%FORMAT(I8)'
  write(prmtop_unit,'(I8)')num_pitorsion_params
  write(word,'(I8)')num_pitorsion_params
  word1 = adjustl(word)
  word2 = '('//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_PI_TORSION_FORCE_CONSTANT'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat)(pitorsion_param(1,k),k=1,num_pitorsion_params)
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_PI_TORSION_PERIODICITY'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat)(pitorsion_param(2,k),k=1,num_pitorsion_params)
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_PI_TORSION_PHASE'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat)(pitorsion_param(3,k),k=1,num_pitorsion_params)
end subroutine AM_PARM_write_pitorsion_params
!------------------------------------------------------------------
subroutine AM_PARM_get_strbend_params(analout_unit)
  integer, intent(in) :: analout_unit

  integer :: success,num,ier1,ier2,ier3,ier4
  character(len=120) :: line
  integer :: ios,n,j,k,dim1
  integer, allocatable :: stretch_bend_parm_ptr(:)
  double precision, allocatable :: old_param(:,:)

  rewind(analout_unit)
  call AM_PARM_get_one_num(analout_unit,'Stretch-Bends',num,success)
  if ( success > 0 )then
    num_stretch_bends = num
    write(6,*)'num_stretch_bends = ',num_stretch_bends
  else
    write(6,*)'no success finding Stretch-Bends!'
    num_stretch_bends = 0
    return
  endif
  ! allocate
  allocate(stretch_bend_list(4,num_stretch_bends),stat=ier1)!3 atoms, 1 ptr
  allocate(stretch_bend_parm_ptr(num_stretch_bends),stat=ier2)
  allocate(stretch_bend_param(4,num_stretch_bends),stat=ier3)
  allocate(old_param(4,num_stretch_bends),stat=ier4)
  if ( (ier1 /= 0) .or. (ier2 /= 0) .or. (ier3 /= 0) .or. (ier4 /= 0) ) then
    write(6,*)'stretch_bend allocation prob!!'
    stop
  endif
  call AM_PARM_get_section(analout_unit, &
           'Stretch-Bend Parameters',success)
  if ( success == 0 )then
    write(6,*)'no success getting stretch_bend Parameters'
    stop
  endif
  ! read 3 lines
  read(analout_unit,'(A)',iostat=ios)line
  read(analout_unit,'(A)',iostat=ios)line
  read(analout_unit,'(A)',iostat=ios)line
  do n = 1,num_stretch_bends
    read(analout_unit,*)k,(stretch_bend_list(j,k),j=1,3), &
                        (old_param(j,k),j=1,4)
  enddo
  dim1 = 4
  call AM_PARM_compress_params(dim1,num_stretch_bends,  &
                      old_param, &
                      num_strbend_params,stretch_bend_param,  &
                      stretch_bend_parm_ptr)
  do n = 1,num_stretch_bends
    stretch_bend_list(4,n) = stretch_bend_parm_ptr(n)
  enddo
  if ( allocated(stretch_bend_parm_ptr) )deallocate(stretch_bend_parm_ptr)
  if ( allocated(old_param) )deallocate(old_param)
end subroutine AM_PARM_get_strbend_params
!------------------------------------------------------------------
subroutine AM_PARM_write_strbend_params(prmtop_unit)
  integer, intent(in) :: prmtop_unit

  character(len=8)word,word1
  character(len=16)word2
  integer j,k
  ! dump out stretch_bends
  if ( num_stretch_bends == 0 )return
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_STRETCH_BEND_NUM_LIST'
  write(prmtop_unit,'(a)')'%FORMAT(I8)'
  write(prmtop_unit,'(I8)')num_stretch_bends
  write(word,'(I8)')num_stretch_bends
  word1 = adjustl(word)
  word2 = '(4,'//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_STRETCH_BEND_LIST'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')'%FORMAT(10I8)'
  write(prmtop_unit,'(10I8)')((stretch_bend_list(j,k),j=1,4),  &
                               k=1,num_stretch_bends)
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_STRETCH_BEND_NUM_PARAMS'
  write(prmtop_unit,'(a)')'%FORMAT(I8)'
  write(prmtop_unit,'(I8)')num_strbend_params
  write(word,'(I8)')num_strbend_params
  word1 = adjustl(word)
  word2 = '('//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_STRETCH_BEND_FORCE_CONSTANT'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat)(stretch_bend_param(1,k),k=1,num_strbend_params)
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_STRETCH_BEND_ANGLE_EQUIL_VALUE'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat)(stretch_bend_param(2,k),k=1,num_strbend_params)
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_STRETCH_BEND_BOND1_EQUIL_VALUE'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat)(stretch_bend_param(3,k),k=1,num_strbend_params)
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_STRETCH_BEND_BOND2_EQUIL_VALUE'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat)(stretch_bend_param(4,k),k=1,num_strbend_params)
end subroutine AM_PARM_write_strbend_params
!------------------------------------------------------------------
subroutine AM_PARM_get_tortor_params(frcfield_unit,analout_unit)
  integer, intent(in) :: frcfield_unit,analout_unit

  integer :: n,nn,success,dim1,dim2,k1,k2,dmax,res,num
  integer ier1,ier2,ier3,ier4,ier5,ier6,ios,j,k,ntok
  integer,allocatable :: tortorclass(:,:)
  double precision,allocatable :: bs(:),cs(:),ds(:),tmp1(:),tmp2(:),tmp3(:)
  character(len=120) :: line

  ! get num_tor_tors
  rewind(analout_unit)
  call AM_PARM_get_one_num(analout_unit,'Torsion-Torsions',num,success)
  if ( success > 0 )then
    num_tor_tors = num
    write(6,*)'num_tor_tors = ',num_tor_tors
  else
    write(6,*)'no success finding Torsion-Torsions'
    num_tor_tors = 0
    return
  endif
  call AM_PARM_num_matching_lines(frcfield_unit,'tortors',num_tortor_tables)
  write(6,*)'num tortor tables = ',num_tortor_tables
  allocate(tortor_table(num_tortor_tables),stat=ier1)
  allocate(tortorclass(5,num_tortor_tables),stat=ier2)
  if ( (ier1 /= 0) .or. (ier2 /= 0) ) then
    write(6,*)'tortor table allocation prob!!'
    stop
  endif
  rewind(frcfield_unit)
  do n = 1,num_tortor_tables
    call AM_PARM_get_next_tortor_table(frcfield_unit,success,  &
                                         tortorclass(:,n),dim1,dim2)
    if ( success == 0 )then
      write(6,*)'problem with AM_PARM_get_next_tortor_table!'
      stop
    endif
    tortor_table(n)%dim1 = dim1
    tortor_table(n)%dim2 = dim2
    allocate(tortor_table(n)%angle1(dim1),stat=ier1)
    allocate(tortor_table(n)%angle2(dim2),stat=ier2)
    allocate(tortor_table(n)%func(dim1,dim2),stat=ier3)
    allocate(tortor_table(n)%dfunc_dangle1(dim1,dim2),stat=ier4)
    allocate(tortor_table(n)%dfunc_dangle2(dim1,dim2),stat=ier5)
    allocate(tortor_table(n)%d2func_dangle1_dangle2(dim1,dim2),stat=ier6)
    if ( (ier1 /= 0) .or. (ier2 /= 0) .or. (ier3 /= 0) .or.  &
         (ier4 /= 0) .or. (ier5 /= 0) .or. (ier6 /= 0) ) then
      write(6,*)'tortor table allocation prob!!'
      stop
    endif
    do k2 = 1,dim2
      do k1 = 1,dim1
        read(frcfield_unit,*)tortor_table(n)%angle1(k1), &
                             tortor_table(n)%angle2(k2),  &
                             tortor_table(n)%func(k1,k2)
      enddo
    enddo
    dmax = max(dim1,dim2)
    allocate(bs(dmax+1),stat=ier1)
    allocate(cs(dmax+1),stat=ier2)
    allocate(ds(dmax+1),stat=ier3)
    allocate(tmp1(dmax+1),stat=ier4)
    allocate(tmp2(dmax+1),stat=ier5)
    allocate(tmp3(5*dmax+1),stat=ier6)
    if ( (ier1 /= 0) .or. (ier2 /= 0) .or. (ier3 /= 0) .or.  &
         (ier4 /= 0) .or. (ier5 /= 0) .or. (ier6 /= 0) ) then
      write(6,*)'tortor table allocation prob!!'
      stop
    endif
    call AM_PARM_tortorspline( &
           dim1,dim2,tortor_table(n)%angle1,tortor_table(n)%angle2, &
           tortor_table(n)%func,tortor_table(n)%dfunc_dangle1, &
           tortor_table(n)%dfunc_dangle2, &
           tortor_table(n)%d2func_dangle1_dangle2,  &
           bs,cs,ds,tmp1,tmp2,tmp3) 
    deallocate(bs)
    deallocate(cs)
    deallocate(ds)
    deallocate(tmp1)
    deallocate(tmp2)
    deallocate(tmp3)
  enddo
  ! fill the tortors
  allocate(tortor_list(6,num_tor_tors),stat=ier1)! 5 atoms plus parm ptr
  if ( (ier1 /= 0) ) then
    write(6,*)'Torsion-Torsions allocation prob!!'
    stop
  endif
  ! read in atom list
  call AM_PARM_get_section(analout_unit, &
           'Torsion-Torsion Parameters',success)
  if ( success == 0 )then
    write(6,*)'no success getting Torsion-Torsion Parameters'
    stop
  endif
  ! read 3 lines
  read(analout_unit,'(A)',iostat=ios)line
  read(analout_unit,'(A)',iostat=ios)line
  read(analout_unit,'(A)',iostat=ios)line
  ! the next line is the first we will actually parse
  read(analout_unit,'(A)',iostat=ios)line
  do n = 1,num_tor_tors
    read(line,*)k,(tortor_list(j,k),j=1,5),dim1,dim2
    read(analout_unit, '(A)') line
    ! See if the table definition is supplied here. If there are only 3 words on
    ! this line, we have a table definition that we need to skip over
    call get_num_tokens(line, ntok)
    if (ntok == 3) then
      ! Now read over the table definition and ignore them. This is currently
      ! extracted from the frcfield file. It can be changed to be read out of
      ! this section instead, but for now keep things as-is so we have something
      ! to compare to.
      do nn = 1, dim1*dim2
        read(analout_unit,'(A)') line
      end do
    end if
    call AM_PARM_find_tortor_table(tortor_list(:,k),num_tortor_tables, &
             tortorclass,res)
    if ( res == 0 )then
      write(6,*)'trouble finding table for tortor num ',n
    else
      tortor_list(6,k) = res
    endif
  enddo
  ! deallocate
  deallocate(tortorclass)
end subroutine AM_PARM_get_tortor_params
!------------------------------------------------------------------
subroutine AM_PARM_write_tortor_params(prmtop_unit)
  integer, intent(in) :: prmtop_unit

  character(len=2)word
  character(len=8)word1,word2,word3
  character(len=25)word4
  integer j,k,n
  ! output tortor info
  if ( num_tor_tors == 0 )return
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_TORSION_TORSION_NUM_LIST'
  write(prmtop_unit,'(a)')'%FORMAT(I8)'
  write(prmtop_unit,'(I8)')num_tor_tors
  write(word1,'(I8)')num_tor_tors
  word2 = adjustl(word1)
  word3 = '(6,'//word2(1:len_trim(word2))//')'
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_TORSION_TORSION_LIST'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word3
  write(prmtop_unit,'(a)')'%FORMAT(10I8)'
  write(prmtop_unit,'(10I8)')((tortor_list(j,k),j=1,6),  &
                               k=1,num_tor_tors)
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_TORSION_TORSION_NUM_PARAMS'
  write(prmtop_unit,'(a)')'%FORMAT(I8)'
  write(prmtop_unit,'(I8)')num_tortor_tables
  do n = 1,num_tortor_tables
    write(word,'(i2.2)')n
    write(prmtop_unit,'(a)')  &
       '%FLAG AMOEBA_TORSION_TORSION_TORTOR_TABLE_'//word//'_DIMS'
    write(prmtop_unit,'(a)')'%COMMENT   dimension = (2)'
    write(prmtop_unit,'(a)')'%FORMAT(2I8)'
    write(prmtop_unit,'(2I8)')tortor_table(n)%dim1,tortor_table(n)%dim2
    write(word1,'(I8)')tortor_table(n)%dim1
    word2 = adjustl(word1)
    write(word1,'(I8)')tortor_table(n)%dim2
    word3 = adjustl(word1)
    word4 = '('//word2(1:len_trim(word2))//')'
    write(prmtop_unit,'(a)')  &
       '%FLAG AMOEBA_TORSION_TORSION_TORTOR_TABLE_'//word//'_ANGLE1'
    write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word4
    write(prmtop_unit,'(a)')outrealformat
    write(prmtop_unit,realformat)tortor_table(n)%angle1
    write(prmtop_unit,'(a)')  &
       '%FLAG AMOEBA_TORSION_TORSION_TORTOR_TABLE_'//word//'_ANGLE2'
    word4 = '('//word3(1:len_trim(word3))//')'
    write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word4
    write(prmtop_unit,'(a)')outrealformat
    write(prmtop_unit,realformat)tortor_table(n)%angle2
    write(prmtop_unit,'(a)') &
       '%FLAG AMOEBA_TORSION_TORSION_TORTOR_TABLE_'//word//'_FUNC'
    word4 = '('//word2(1:len_trim(word2))//','//word3(1:len_trim(word3))//')'
    write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word4
    write(prmtop_unit,'(a)')outrealformat
    write(prmtop_unit,realformat)tortor_table(n)%func
    write(prmtop_unit,'(a)')  &
       '%FLAG AMOEBA_TORSION_TORSION_TORTOR_TABLE_'//word//'_DFUNC_DANGLE1'
    write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word4
    write(prmtop_unit,'(a)')outrealformat
    write(prmtop_unit,realformat)tortor_table(n)%dfunc_dangle1
    write(prmtop_unit,'(a)')  &
       '%FLAG AMOEBA_TORSION_TORSION_TORTOR_TABLE_'//word//'_DFUNC_DANGLE2'
    write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word4
    write(prmtop_unit,'(a)')outrealformat
    write(prmtop_unit,realformat)tortor_table(n)%dfunc_dangle2
    write(prmtop_unit,'(a)')  &
  '%FLAG AMOEBA_TORSION_TORSION_TORTOR_TABLE_'//word//'_D2FUNC_DANGLE1_DANGLE2'
    write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word4
    write(prmtop_unit,'(a)')outrealformat
    write(prmtop_unit,realformat)tortor_table(n)%d2func_dangle1_dangle2
  enddo
end subroutine AM_PARM_write_tortor_params
!------------------------------------------------------------------
! utility subroutines
!------------------------------------------------------------------
subroutine AM_PARM_find_tortor_table(atomlist,num_tab,tortorclass,res)
  integer,intent(in) :: atomlist(5),num_tab
  integer,intent(in) :: tortorclass(5,num_tab)
  integer,intent(out) :: res

  integer :: n
  do n = 1,num_tab
    if ( ( atomclass(atomlist(1)) == tortorclass(1,n) ) .and.  &
         ( atomclass(atomlist(2)) == tortorclass(2,n) ) .and.  &
         ( atomclass(atomlist(3)) == tortorclass(3,n) ) .and.  &
         ( atomclass(atomlist(4)) == tortorclass(4,n) ) .and.  &
         ( atomclass(atomlist(5)) == tortorclass(5,n) ) ) then
      res = n
      return
    endif
    if ( ( atomclass(atomlist(5)) == tortorclass(1,n) ) .and.  &
         ( atomclass(atomlist(4)) == tortorclass(2,n) ) .and.  &
         ( atomclass(atomlist(3)) == tortorclass(3,n) ) .and.  &
         ( atomclass(atomlist(2)) == tortorclass(4,n) ) .and.  &
         ( atomclass(atomlist(1)) == tortorclass(5,n) ) ) then
      res = n
      return
    endif
  enddo
  res = 0 
end subroutine AM_PARM_find_tortor_table
!------------------------------------------------------------------
subroutine AM_PARM_tortorspline(nphi1,nphi2,phi1,phi2,func, &
                               grad1_func,grad2_func,grad12_func, &
                               bs,cs,ds,tmp1,tmp2,tmp3)
  integer :: nphi1,nphi2
  double precision,intent(in) :: phi1(nphi1),phi2(nphi2),func(nphi1,nphi2)
  double precision,intent(out) ::   grad1_func(nphi1,nphi2), &
                                    grad2_func(nphi1,nphi2),  &
                                    grad12_func(nphi1,nphi2)
  double precision,intent(out) :: bs(*),cs(*),ds(*),tmp1(*),tmp2(*),tmp3(*)

  integer nx,ny,m,j,k

  ! first get grad1; fit each row to a 1d periodic spline
  nx = nphi1 - 1
  do k = 1,nphi1
    tmp1(k) = phi1(k)
  enddo
  do j = 1,nphi2
    do k = 1,nphi1
      tmp2(k) = func(k,j)
    enddo
    call isplpe(nx,tmp1,tmp2,1,bs,cs,ds,tmp3(1),tmp3(nx+2), &
                tmp3(2*nx+2),tmp3(3*nx+2),tmp3(4*nx+2))
    do k = 1,nphi1
      grad1_func(k,j) = bs(k)
    enddo
  enddo
  ! next get grad2; fit each column to a 1d periodic spline
  ny = nphi2 - 1
  do k = 1,nphi2
    tmp1(k) = phi2(k)
  enddo
  do j = 1,nphi1
    do k = 1,nphi2
!      tmp2(k) = func(m+(k-1)*nphi1)
      tmp2(k) = func(j,k)
    enddo
    call isplpe(ny,tmp1,tmp2,1,bs,cs,ds,tmp3(1),tmp3(ny+2), &
                tmp3(2*ny+2),tmp3(3*ny+2),tmp3(4*ny+2))
    do k = 1,nphi2
      grad2_func(j,k) = bs(k)
    enddo
  enddo
  ! finally get grad2 of grad1--fit each column to a 1d periodic spline
  do k = 1,nphi2
    tmp1(k) = phi2(k)
  enddo
  do j = 1,nphi1
    do k = 1,nphi2
      tmp2(k) = grad1_func(j,k)
    enddo
    call isplpe(ny,tmp1,tmp2,1,bs,cs,ds,tmp3(1),tmp3(ny+2), &
                tmp3(2*ny+2),tmp3(3*ny+2),tmp3(4*ny+2))
    do k = 1,nphi2
      grad12_func(j,k) = bs(k)
    enddo
  enddo
  return
end subroutine AM_PARM_tortorspline
!------------------------------------------------------------------
subroutine AM_PARM_num_matching_lines(in_unit,search_string,num)
  integer, intent(in) :: in_unit
  character(len=*), intent(in) :: search_string
  integer, intent(out) :: num

  character(len=120) :: line
  integer :: ios,success
  num = 0
  rewind(in_unit)
  do
    read(in_unit,'(A)',iostat=ios)line
    if ( ios /= 0 )return
    success = index(line,search_string)
    if ( success > 0 )then
      num = num + 1
    endif
  enddo
end subroutine AM_PARM_num_matching_lines
!------------------------------------------------------------------
subroutine AM_PARM_get_one_num(analout_unit,search_string,num,success)
  integer, intent(in) :: analout_unit
  character(len=*), intent(in) :: search_string
  integer, intent(out) :: num,success

  character(len=120) :: line
  integer :: ios
  ios = 0
  success = 0
  do while ( ios == 0 )
    read(analout_unit,'(A)',iostat=ios)line
    success = index(line,search_string)
    num = -1
    if ( success > 0 )then
      read(line(success+len(search_string)+1:),*)num
      return
    endif
  enddo
end subroutine AM_PARM_get_one_num
!------------------------------------------------------------------
subroutine AM_PARM_get_section(analout_unit,search_string,success)
  integer, intent(in) :: analout_unit
  character(len=*), intent(in) :: search_string
  integer, intent(out) :: success

  character(len=130) :: line,lline
  integer :: ios
  ios = 0
  success = 0
  do while ( ios == 0 )
    read(analout_unit,'(A)',iostat=ios)line
    lline = adjustl(line)
    success = index(lline,search_string)
    if ( success == 1 )then
      return
    endif
  enddo
end subroutine AM_PARM_get_section
!------------------------------------------------------------------
subroutine AM_PARM_get_numtorterms(line,numtor)
  character(len=*), intent(in) :: line
  integer, intent(out) :: numtor
  
  integer :: start,res
  start = 1
  numtor = 0
  do
    res = index(line(start:),"/") 
    if ( res == 0 )return
    numtor = numtor + 1
    start = start + res
  enddo
end subroutine AM_PARM_get_numtorterms
!------------------------------------------------------------------
subroutine AM_PARM_compress_params(dim1,num_old_params,old_params, &
                      num_new_params,new_params,paramptr)
  integer, intent(in) :: dim1,num_old_params
  double precision, intent(in) :: old_params(dim1,num_old_params)
  integer, intent(out) :: num_new_params,  &
                          paramptr(num_old_params)
  double precision, intent(out) :: new_params(dim1,num_old_params)

  integer :: j,k,n,res
  double precision diff,small
  ! copy 1st one to start things
  num_new_params = 1
  small = 1.d-6
  do j = 1,dim1
    new_params(j,1) = old_params(j,1)
  enddo
  paramptr(1) = 1
  do n = 2,num_old_params
    res = -1
    do k = 1,num_new_params
      diff = 0.d0
      do j = 1,dim1
        diff = diff + abs(new_params(j,k)-old_params(j,n))
      enddo
      if ( diff < small )then
        res = k
        exit
      endif
    enddo
    if ( res < 0 )then
      num_new_params = num_new_params + 1
      do j = 1,dim1
        new_params(j,num_new_params) = old_params(j,n)
      enddo
      paramptr(n) = num_new_params
    else
      paramptr(n) = res
    endif
  enddo
end subroutine AM_PARM_compress_params
!------------------------------------------------------------------
subroutine AM_PARM_bond_nghb_list()
  integer :: n,k1,k2,sizlist,off

  integer, allocatable :: numnghbors(:)
  integer :: ier,siz_list
  ! allocate offset
  allocate(offset_nghbors(num_atoms),stat=ier)
  if ( ier /= 0 )then
    write(6,*)'bond_nghb_list allocation prob!!'
    stop
  endif
  allocate(numnghbors(num_atoms),stat=ier)
  if ( ier /= 0 )then
    write(6,*)'bond_nghb_list allocation prob!!'
    stop
  endif
  offset_nghbors(1) = 0
  do n = 2,num_atoms
    offset_nghbors(n) = offset_nghbors(n-1) + atomvalence(n-1)
  enddo
  siz_list = offset_nghbors(num_atoms) + atomvalence(num_atoms)
  allocate(bond_nghbors(siz_list),stat=ier)
  if ( ier /= 0 )then
    write(6,*)'bond_nghb_list allocation prob!!'
    stop
  endif
  do n = 1,num_atoms
    numnghbors(n) = 0
  enddo
  do n = 1,num_bonds
    k1 = bond_list(1,n)
    k2 = bond_list(2,n)
    numnghbors(k1) = numnghbors(k1) + 1
    if ( numnghbors(k1) > atomvalence(k1) )then
      write(6,*)'too many nghbors for ',k1
      stop
    endif
    if ( numnghbors(k2) > atomvalence(k2) )then
      write(6,*)'too many nghbors for ',k2
      stop
    endif
    numnghbors(k2) = numnghbors(k2) + 1
    bond_nghbors(offset_nghbors(k1)+numnghbors(k1)) = k2
    bond_nghbors(offset_nghbors(k2)+numnghbors(k2)) = k1
  enddo 
  do n = 1,num_atoms
    if ( numnghbors(n) /= atomvalence(n) )then
      write(6,*)'num nghbors not equal to valence for ',n
      stop
    endif
  enddo
  deallocate(numnghbors)
end subroutine AM_PARM_bond_nghb_list
!------------------------------------------------------------------
subroutine AM_PARM_find_vdw_parent()

  integer :: n,ier

  allocate(vdw_atom_parent(num_atoms),stat=ier)
  if ( ier /= 0 )then
    write(6,*)'vdw_atom_parent allocation prob!!'
    stop
  endif
  do n = 1,num_atoms
     if ( atomic_number(n) == 1 )then
        vdw_atom_parent(n) = bond_nghbors(offset_nghbors(n)+1)
     else ! atom is its own vdw parent
        vdw_atom_parent(n) = n
     endif
  enddo
end subroutine AM_PARM_find_vdw_parent
!------------------------------------------------------------------
subroutine AM_PARM_get_mol_list()
  integer,allocatable :: atom_flag(:),start_mol(:)
  integer :: k,m,n,off,num,base,ier

  allocate(atom_flag(num_atoms),stat=ier)
  if ( ier /= 0 )then
    write(6,*)'AM_PARM_get_mol_list: problem allocating atom_flag'
    stop
  endif
  ! first pass, get num_molecules
  do n = 1,num_atoms
    atom_flag(n) = 0
  enddo
  num_molecules = 0
  do n = 1,num_atoms
    if ( atom_flag(n) == 0 )then
      num_molecules = num_molecules + 1
      base = n
    endif
    off = offset_nghbors(n)
    num = atomvalence(n)
    do k = 1,num
      m = bond_nghbors(off+k)
      if ( m > n )then
        atom_flag(m) = base
      endif
    enddo
  enddo
  write(6,*)'num_molecules = ',num_molecules
  allocate(num_atoms_in_molecule(num_molecules),stat=ier)
  if ( ier /= 0 )then
    write(6,*)'AM_PARM_get_mol_list: problem allocating num_atoms_in_molecule'
    stop
  endif
  allocate(start_mol(num_molecules),stat=ier)
  if ( ier /= 0 )then
    write(6,*)'AM_PARM_get_mol_list: problem allocating start_mol'
    stop
  endif
  ! 2nd pass, get num_atoms_in_molecule
  m = 0
  do n = 1,num_atoms
    if ( atom_flag(n) == 0 )then
      m = m + 1
      start_mol(m) = n
    endif 
  enddo
  do m = 1,num_molecules-1
    num_atoms_in_molecule(m) = start_mol(m+1) - start_mol(m)
  enddo
  num_atoms_in_molecule(num_molecules) = num_atoms+1 - start_mol(num_molecules)
  deallocate(atom_flag)
  deallocate(start_mol)
end subroutine AM_PARM_get_mol_list
!------------------------------------------------------------------
subroutine AM_PARM_get_multipole(analout_unit)
  integer,intent(in) :: analout_unit

  character(len=256) :: line
  integer success,num,ios,ier1,ier2,ier3,numf,j,k,n,i0,iz,ix,iy
  double precision cg,Mux,Muy,Muz,THxx,THyy,THzz,THxy,THxz,THyz
  character(len=12)ftype
  ! get num framelists
  rewind(analout_unit)
  call AM_PARM_get_one_num(analout_unit,'Polarizable Multipoles',num,success)
  if ( success > 0 )then
    num_multipoles = num
    write(6,*)'num_multipoles = ',num_multipoles
  else
    write(6,*)'no success finding Polarizable Multipoles!'
    num_multipoles = 0
    return
  endif
  num_frame_list = 3*num_multipoles ! over-estimate
  ! each frame has 2 or 3 frame def list entries
  allocate(frame_list(num_frame_list),stat=ier1)
  allocate(chiral_frame_list(num_multipoles),stat=ier2)
  allocate(local_multipole(10,num_multipoles),stat=ier3)
  if ( (ier1/= 0) .or. (ier2/= 0) .or. (ier3/= 0) )then
    write(6,*)'problem allocating multipole lists'
    stop
  endif
  call AM_PARM_get_section(analout_unit, &
           'Atomic Multipole Parameters',success)
  if ( success == 0 )then
    write(6,*)'no success getting Atomic Multipole Parameters'
    stop
  endif
  ! read 3 lines
  read(analout_unit,'(A)',iostat=ios)line
  read(analout_unit,'(A)',iostat=ios)line
  read(analout_unit,'(A)',iostat=ios)line
  num_frame_list = 0
  num_chiral_frame_list = 0
  do n = 1,num_multipoles
    read(analout_unit,'(a)',iostat=ios)line
    call AM_PARM_num_fields(line,numf)
    if ( numf == 6 )then
      read(line,*)k,i0,iz,ix,ftype,cg
    elseif (numf == 7 )then
      read(line,*)k,i0,iz,ix,iy,ftype,cg
      num_chiral_frame_list = num_chiral_frame_list + 1
      if ( iy > 0 )then
        chiral_frame_list(num_chiral_frame_list)%frame_index = n
        chiral_frame_list(num_chiral_frame_list)%fourth_atom = iy
        chiral_frame_list(num_chiral_frame_list)%chirality = 1
      else
        chiral_frame_list(num_chiral_frame_list)%frame_index = n
        chiral_frame_list(num_chiral_frame_list)%fourth_atom = abs(iy)
        chiral_frame_list(num_chiral_frame_list)%chirality = -1
      endif
    else
      write(6,*)'wrong number of fields in multipole line: numf = ',numf
      stop
    endif
    if ( ftype(1:8) == 'Z-then-X' )then
      num_frame_list = num_frame_list + 1 ! 1 for def point 1
      frame_list(num_frame_list)%frame_index = n
      frame_list(num_frame_list)%frame_point_number = 1
      frame_list(num_frame_list)%vector_tail_index = i0
      frame_list(num_frame_list)%vector_head_index = iz
      frame_list(num_frame_list)%num_vectors = 1
      num_frame_list = num_frame_list + 1 ! 1 for def point 2
      frame_list(num_frame_list)%frame_index = n
      frame_list(num_frame_list)%frame_point_number = 2
      frame_list(num_frame_list)%vector_tail_index = i0
      frame_list(num_frame_list)%vector_head_index = ix
      frame_list(num_frame_list)%num_vectors = 1
    elseif (  ftype(1:8) == 'Bisector' )then
      num_frame_list = num_frame_list + 1 ! 2 for def point 1
      frame_list(num_frame_list)%frame_index = n
      frame_list(num_frame_list)%frame_point_number = 1
      frame_list(num_frame_list)%vector_tail_index = i0
      frame_list(num_frame_list)%vector_head_index = iz
      frame_list(num_frame_list)%num_vectors = 2
      num_frame_list = num_frame_list + 1 ! 2 for def point 1
      frame_list(num_frame_list)%frame_index = n
      frame_list(num_frame_list)%frame_point_number = 1
      frame_list(num_frame_list)%vector_tail_index = i0
      frame_list(num_frame_list)%vector_head_index = ix
      frame_list(num_frame_list)%num_vectors = 2
      num_frame_list = num_frame_list + 1 ! 1 for def point 2
      frame_list(num_frame_list)%frame_index = n
      frame_list(num_frame_list)%frame_point_number = 2
      frame_list(num_frame_list)%vector_tail_index = i0
      frame_list(num_frame_list)%vector_head_index = ix
      frame_list(num_frame_list)%num_vectors = 1
    else
      write(6,*)'unknown frame type: ftype = ',ftype
      stop
    endif
    read(analout_unit,*)Mux,Muy,Muz
    read(analout_unit,*)THxx
    read(analout_unit,*)THxy,THyy
    read(analout_unit,*)THxz,THyz,THzz
    local_multipole(1,n) = cg
    local_multipole(2,n) = Mux
    local_multipole(3,n) = Muy
    local_multipole(4,n) = Muz
    local_multipole(5,n) = 0.5d0*THxx
    local_multipole(6,n) = 0.5d0*THyy
    local_multipole(7,n) = 0.5d0*THzz
    local_multipole(8,n) = THxy
    local_multipole(9,n) = THxz
    local_multipole(10,n) = THyz
  enddo
  write(6,*)'number of frame def list = ',num_frame_list

end subroutine AM_PARM_get_multipole
!------------------------------------------------------------------
subroutine AM_PARM_write_polar(prmtop_unit)
  integer,intent(in) :: prmtop_unit

  character(len=8)word,word1
  character(len=17)word2
  integer :: n
  
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_POLARIZABILITY_NUM_LIST'
  write(prmtop_unit,'(a)')'%FORMAT(I8)'
  write(prmtop_unit,'(I8)')num_atoms
  write(word,'(I8)')num_atoms
  word1 = adjustl(word)
  word2 = '('//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_POLARIZABILITY_LIST'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat)(polarizability(n),n=1,num_atoms)
end subroutine AM_PARM_write_polar
!------------------------------------------------------------------
subroutine AM_PARM_write_adjust(prmtop_unit)
  integer,intent(in) :: prmtop_unit

  character(len=8)word,word1
  character(len=17)word2
  integer :: j,n
  double precision vdw_wt(9),mpole_wt(9),polar_wt(9),direct_wt(9),mutual_wt(9)
  
  if ( size_excluded_list == 0 )return
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_ADJUST_NUM_LIST'
  write(prmtop_unit,'(a)')'%FORMAT(I8)'
  write(prmtop_unit,'(I8)')size_excluded_list
  write(word,'(I8)')size_excluded_list
  word1 = adjustl(word)
  word2 = '(3,'//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_ADJUST_LIST'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')'%FORMAT(10I8)'
  write(prmtop_unit,'(10I8)')((adjust_list(j,n),j=1,3),n=1,size_excluded_list)
  ! fill adjust weight arrays
  ! recall indices of 1-2,1-3,1-4 and 1-5 are 1,2,3,4 and polar group adds 5-9
  ! 5 means in polar group and 1-2, 6 means in polar group and 1-3, 
  ! 7 means in polar group and 1-4, 8 means in polar group and 1-5
  ! and 9 means in polar group but not in 1-2 through 1-5
  vdw_wt(1) = 0.d0;  vdw_wt(2) = 0.d0; vdw_wt(3) = 1.d0; vdw_wt(4) = 1.d0
  vdw_wt(5) = 0.d0;  vdw_wt(6) = 0.d0; vdw_wt(7) = 1.d0; vdw_wt(8) = 1.d0
  vdw_wt(9) = 1.d0
  mpole_wt(1) = 0.0d0; mpole_wt(2) = 0.0d0; mpole_wt(3) = 0.4d0
  mpole_wt(4) = 0.8d0; mpole_wt(5) = 0.0d0; mpole_wt(6) = 0.0d0
  mpole_wt(7) = 0.4d0; mpole_wt(8) = 0.8d0; mpole_wt(9) = 1.0d0
  direct_wt(1) = 1.d0; direct_wt(2) = 1.d0; direct_wt(3) = 1.d0
  direct_wt(4) = 1.d0; direct_wt(5) = 0.d0; direct_wt(6) = 0.d0
  direct_wt(7) = 0.d0; direct_wt(8) = 0.d0; direct_wt(9) = 0.d0
  polar_wt(1) = 0.d0; polar_wt(2) = 0.d0; polar_wt(3) = 1.d0
  polar_wt(4) = 1.d0; polar_wt(5) = 0.d0; polar_wt(6) = 0.d0
  ! note polar wt of 1/2 for pair in polar group and 1-4
  polar_wt(7) = 0.5d0; polar_wt(8) = 1.d0; polar_wt(9) = 1.d0
  mutual_wt(1) = 1.d0; mutual_wt(2) = 1.d0; mutual_wt(3) = 1.d0
  mutual_wt(4) = 1.d0; mutual_wt(5) = 1.d0; mutual_wt(6) = 1.d0
  mutual_wt(7) = 1.d0; mutual_wt(8) = 1.d0; mutual_wt(9) = 1.d0
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_ADJUST_VDW_WEIGHTS_LIST'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = (9)'
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat)(vdw_wt(n),n=1,9)
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_ADJUST_MPOLE_WEIGHTS_LIST'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = (9)'
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat)(mpole_wt(n),n=1,9)
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_ADJUST_DIRECT_WEIGHTS_LIST'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = (9)'
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat)(direct_wt(n),n=1,9)
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_ADJUST_POLAR_WEIGHTS_LIST'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = (9)'
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat)(polar_wt(n),n=1,9)
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_ADJUST_MUTUAL_WEIGHTS_LIST'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = (9)'
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat)(mutual_wt(n),n=1,9)
end subroutine AM_PARM_write_adjust
!------------------------------------------------------------------
subroutine AM_PARM_write_vdw(prmtop_unit)
  integer,intent(in) :: prmtop_unit

  character(len=8)word,word1
  character(len=20)word2
  integer :: j,k,n

  if ( num_vdw_types == 0 )return
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_VDW_ATOM_TYPES_NUM_LIST'
  write(prmtop_unit,'(a)')'%FORMAT(I8)'
  write(prmtop_unit,'(I8)')num_atoms
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_VDW_ATOM_TYPES_LIST'
  write(word,'(I8)')num_atoms
  word1 = adjustl(word)
  word2 = '(1,'//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')'%FORMAT(10I8)'
  write(prmtop_unit,'(10I8)')(atom_vdw_type(n),n=1,num_atoms)
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_VDW_ATOM_PARENT_LIST'
  write(word,'(I8)')num_atoms
  word1 = adjustl(word)
  word2 = '(1,'//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')'%FORMAT(10I8)'
  write(prmtop_unit,'(10I8)')(vdw_atom_parent(n),n=1,num_atoms)
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_VDW_PARENT_COORD_WEIGHT_LIST'
  write(word,'(I8)')num_atoms
  word1 = adjustl(word)
  word2 = '(1,'//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat)(vdw_parent_weight(n),n=1,num_atoms)
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_VDW_BUFFER_DELTA'
  write(prmtop_unit,'(a)')'%FORMAT(E16.8)'
  write(prmtop_unit,'(E16.8)')vdw_buffer_delta
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_VDW_BUFFER_GAMMA'
  write(prmtop_unit,'(a)')'%FORMAT(E16.8)'
  write(prmtop_unit,'(E16.8)')vdw_buffer_gamma
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_VDW_PARAMS_NUM_LIST'
  write(prmtop_unit,'(a)')'%FORMAT(I8)'
  write(prmtop_unit,'(I8)')num_vdw_types
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_VDW_MIXED_RADII_LIST'
  write(word,'(I8)')num_vdw_types
  word1 = adjustl(word)
  word2 = '('//word1(1:len_trim(word1))//','//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat) &
         ((vdw_mix_radius(j,k),j=1,num_vdw_types),k=1,num_vdw_types)
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_VDW_MIXED_EPSILONS_LIST'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat) &
         ((vdw_mix_epsilon(j,k),j=1,num_vdw_types),k=1,num_vdw_types)
end subroutine AM_PARM_write_vdw
!------------------------------------------------------------------
subroutine AM_PARM_write_multipole(prmtop_unit)
  integer,intent(in) :: prmtop_unit

  character(len=8)word,word1
  character(len=17)word2
  integer :: j,k,n
  !write out params
  if ( num_multipoles == 0 )return
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_LOCAL_FRAME_MULTIPOLES_NUM_LIST'
  write(prmtop_unit,'(a)')'%FORMAT(I8)'
  write(prmtop_unit,'(I8)')num_multipoles

  write(word,'(I8)')num_multipoles
  word1 = adjustl(word)
  word2 = '(10,'//word1(1:len_trim(word1))//')'
  write(prmtop_unit,'(a)')'%FLAG AMOEBA_LOCAL_FRAME_MULTIPOLES_LIST'
  write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
  write(prmtop_unit,'(a)')outrealformat
  write(prmtop_unit,realformat) &
      ((local_multipole(j,n),j=1,10),n=1,num_multipoles)

  if ( num_chiral_frame_list > 0 )then
    write(prmtop_unit,'(a)')'%FLAG AMOEBA_CHIRAL_FRAME_NUM_LIST'
    write(prmtop_unit,'(a)')'%FORMAT(I8)'
    write(prmtop_unit,'(I8)')num_chiral_frame_list

    write(word,'(I8)')num_chiral_frame_list
    word1 = adjustl(word)
    word2 = '(3,'//word1(1:len_trim(word1))//')'
    write(prmtop_unit,'(a)')'%FLAG AMOEBA_CHIRAL_FRAME_LIST'
    write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
    write(prmtop_unit,'(a)')'%FORMAT(10I8)'
    write(prmtop_unit,'(10I8)') &
                            (chiral_frame_list(n)%frame_index,  &
                             chiral_frame_list(n)%fourth_atom,  &
                             chiral_frame_list(n)%chirality,  &
                             n=1,num_chiral_frame_list)
  endif ! num_chiral_frame_list > 0
  if ( num_frame_list > 0 )then
    write(prmtop_unit,'(a)')'%FLAG AMOEBA_FRAME_DEF_NUM_LIST'
    write(prmtop_unit,'(a)')'%FORMAT(I8)'
    write(prmtop_unit,'(I8)')num_frame_list

    write(word,'(I8)')num_frame_list
    word1 = adjustl(word)
    word2 = '(5,'//word1(1:len_trim(word1))//')'
    write(prmtop_unit,'(a)')'%FLAG AMOEBA_FRAME_DEF_LIST'
    write(prmtop_unit,'(a)')'%COMMENT   dimension = '//word2
    write(prmtop_unit,'(a)')'%FORMAT(10I8)'
    write(prmtop_unit,'(10I8)')(frame_list(n)%frame_index, &
                              frame_list(n)%frame_point_number, &
                              frame_list(n)%vector_tail_index, &
                              frame_list(n)%vector_head_index, &
                              frame_list(n)%num_vectors, &
                              n=1,num_frame_list)
  endif ! num_frame_list > 0
end subroutine AM_PARM_write_multipole
!------------------------------------------------------------------
subroutine AM_PARM_get_polar_params(analout_unit)
  integer,intent(in) :: analout_unit

  character(len=256) :: line
  integer :: pgroup(256)
  integer :: i,j,k,n,num,success,numf, &
             numpol,ier1,ier2,ier3,ios,siz,off

  ! allocate
  allocate(polarizability(num_atoms),stat=ier1)
  allocate(num_in_polar_group(num_atoms),stat=ier2)
  allocate(offset_polar_group(num_atoms),stat=ier3)
  if ( (ier1/=0) .or. (ier2/=0) .or. (ier3/=0) )then
    write(6,*)'AM_PARM_get_polar_params: trouble allocating'
    stop
  endif
  ! first pass, count polar group
  call AM_PARM_get_section(analout_unit, &
           'Dipole Polarizability Parameters',success)
  if ( success == 0 )then
    write(6,*)'no success getting Dipole Polarizability Parameters'
    stop
  endif
  ! read 3 lines
  read(analout_unit,'(A)',iostat=ios)line
  read(analout_unit,'(A)',iostat=ios)line
  read(analout_unit,'(A)',iostat=ios)line
  num_in_polar_group = 0
  do n = 1,num_atoms
    read(analout_unit,'(a)',iostat=ios)line
    call AM_PARM_num_fields(line,numf)
    read(line,*)j,k
    if ( k /= n )then
      write(6,*)'mismatch in line!!'
    endif
    num_in_polar_group(k) = numf - 3
  enddo
  offset_polar_group(1) = 0
  do n = 2,num_atoms
    offset_polar_group(n) = offset_polar_group(n-1) + num_in_polar_group(n-1)
  enddo
  siz = offset_polar_group(num_atoms) + num_in_polar_group(num_atoms)
  allocate(nghb_polargroup(siz),stat=ier1)
  if ( (ier1/=0) )then
    write(6,*)'AM_PARM_get_polar_params: trouble allocating nghb_polargroup'
    stop
  endif
  ! pass 2 fill polar parms
  rewind(analout_unit)
  call AM_PARM_get_section(analout_unit, &
           'Dipole Polarizability Parameters',success)
  if ( success == 0 )then
    write(6,*)'no success getting Dipole Polarizability Parameters'
    stop
  endif
  ! read 3 lines
  read(analout_unit,'(A)',iostat=ios)line
  read(analout_unit,'(A)',iostat=ios)line
  read(analout_unit,'(A)',iostat=ios)line
  do n = 1,num_atoms
    num = num_in_polar_group(n)
    read(analout_unit,*)j,k,polarizability(n),(pgroup(i),i=1,num)
    off = offset_polar_group(n)
    do i = 1,num
      nghb_polargroup(off+i) = pgroup(i)
    enddo
  enddo
end subroutine AM_PARM_get_polar_params
!------------------------------------------------------------------
subroutine AM_PARM_get_nghb_12131415_lists(analout_unit)
  integer,intent(in) :: analout_unit

  integer n,num,ios,ier2,ier3,ier4,ier5,success
  character(len=120) :: line

  rewind(analout_unit)
  call AM_PARM_get_one_num(analout_unit,'Number of 1-2 Pairs',num,success)
  if ( success > 0 )then
    num_12 = num
    write(6,*)'num_12 = ',num_12
  else
    write(6,*)'no success finding Number of 1-2 Pairs!'
    num_12 = 0
  endif
  call AM_PARM_get_one_num(analout_unit,'Number of 1-3 Pairs',num,success)
  if ( success > 0 )then
    num_13 = num
    write(6,*)'num_13 = ',num_13
  else
    write(6,*)'no success finding Number of 1-3 Pairs!'
    num_13 = 0
  endif
  call AM_PARM_get_one_num(analout_unit,'Number of 1-4 Pairs',num,success)
  if ( success > 0 )then
    num_14 = num
    write(6,*)'num_14 = ',num_14
  else
    write(6,*)'no success finding Number of 1-4 Pairs!'
    num_14 = 0
  endif
  call AM_PARM_get_one_num(analout_unit,'Number of 1-5 Pairs',num,success)
  if ( success > 0 )then
    num_15 = num
    write(6,*)'num_15 = ',num_15
  else
    write(6,*)'no success finding Number of 1-5 Pairs!'
    num_15 = 0
  endif
  ! allocate
  ier2 = 0; ier3 = 0; ier4 = 0; ier5 = 0
  if ( num_12>0)allocate(nghb_12(2,num_12),stat=ier2)
  if ( num_13>0)allocate(nghb_13(2,num_13),stat=ier3)
  if ( num_14>0)allocate(nghb_14(2,num_14),stat=ier4)
  if ( num_15>0)allocate(nghb_15(2,num_15),stat=ier5)
  if ( (ier2/=0) .or. (ier3/=0) .or. (ier4/=0) .or. (ier5/=0) )then
    write(6,*)'AM_PARM_get_nb_modify_lists: trouble allocating nghb_lists'
    stop
  endif
  rewind(analout_unit)
  if ( num_12 > 0 )then
    call AM_PARM_get_section(analout_unit, &
           'List of 1-2 Connected Atomic Interactions',success)
    if ( success == 0 )then
      write(6,*)'no success getting List of 1-2 Connected Atomic Interactions'
      stop
    endif
    ! read a lines
    read(analout_unit,'(A)',iostat=ios)line
    do n = 1,num_12
      read(analout_unit,*)nghb_12(1,n),nghb_12(2,n)
    enddo
  endif

  if ( num_13 > 0 )then
    call AM_PARM_get_section(analout_unit, &
           'List of 1-3 Connected Atomic Interactions',success)
    if ( success == 0 )then
      write(6,*)'no success getting List of 1-3 Connected Atomic Interactions'
      stop
    endif
    ! read a lines
    read(analout_unit,'(A)',iostat=ios)line
    do n = 1,num_13
      read(analout_unit,*)nghb_13(1,n),nghb_13(2,n)
    enddo
  endif

  if ( num_14 > 0 )then
    call AM_PARM_get_section(analout_unit, &
           'List of 1-4 Connected Atomic Interactions',success)
    if ( success == 0 )then
      write(6,*)'no success getting List of 1-4 Connected Atomic Interactions'
      stop
    endif
    ! read a lines
    read(analout_unit,'(A)',iostat=ios)line
    do n = 1,num_14
      read(analout_unit,*)nghb_14(1,n),nghb_14(2,n)
    enddo
  endif

  if ( num_15 > 0 )then
    call AM_PARM_get_section(analout_unit, &
           'List of 1-5 Connected Atomic Interactions',success)
    if ( success == 0 )then
      write(6,*)'no success getting List of 1-5 Connected Atomic Interactions'
      stop
    endif
    ! read a lines
    read(analout_unit,'(A)',iostat=ios)line
    do n = 1,num_15
      read(analout_unit,*)nghb_15(1,n),nghb_15(2,n)
    enddo
  endif
end subroutine AM_PARM_get_nghb_12131415_lists
!------------------------------------------------------------------
subroutine AM_PARM_get_excluded_atom_list(analout_unit)
  integer,intent(in) :: analout_unit

  integer,allocatable :: num_pre_exclude(:),pre_exclude_list(:), &
                         off_pre_exclude(:),flag(:),off_excluded(:)
  integer, allocatable :: pre_mask(:),mask_flag(:)
  integer :: j,k,n,maxflag,ier1,ier2,ier3,ier4, &
             siz,off1,off2,offj,offk,num,num_ex
  allocate(num_pre_exclude(num_atoms),stat=ier1)
  allocate(off_pre_exclude(num_atoms),stat=ier2)
  allocate(flag(num_atoms),stat=ier3)
  allocate(mask_flag(num_atoms),stat=ier4)
  if ( (ier1/=0) .or. (ier2/=0) .or. (ier3/=0) .or. (ier4/=0) )then
    write(6,*)'AM_PARM_get_excluded_atom_list: trouble allocating'
    stop
  endif
  num_pre_exclude = 0
  flag = 0
  ! count them
  ! start with polar group
  do n = 1,num_atoms
    num_pre_exclude(n) = num_pre_exclude(n) + num_in_polar_group(n)
  enddo
  ! next 1-2
  do n = 1,num_12
    j = nghb_12(1,n)
    k = nghb_12(2,n)
    num_pre_exclude(j) = num_pre_exclude(j) + 1
    num_pre_exclude(k) = num_pre_exclude(k) + 1
  enddo 
  ! next 1-3
  do n = 1,num_13
    j = nghb_13(1,n)
    k = nghb_13(2,n)
    num_pre_exclude(j) = num_pre_exclude(j) + 1
    num_pre_exclude(k) = num_pre_exclude(k) + 1
  enddo 
  ! next 1-4
  do n = 1,num_14
    j = nghb_14(1,n)
    k = nghb_14(2,n)
    num_pre_exclude(j) = num_pre_exclude(j) + 1
    num_pre_exclude(k) = num_pre_exclude(k) + 1
  enddo 
  ! next 1-5
  do n = 1,num_15
    j = nghb_15(1,n)
    k = nghb_15(2,n)
    num_pre_exclude(j) = num_pre_exclude(j) + 1
    num_pre_exclude(k) = num_pre_exclude(k) + 1
  enddo 
  off_pre_exclude(1) = 0
  do n = 2,num_atoms
    off_pre_exclude(n) = off_pre_exclude(n-1) + num_pre_exclude(n-1)
  enddo
  siz = off_pre_exclude(num_atoms) + num_pre_exclude(num_atoms)
  allocate(pre_exclude_list(siz),stat=ier1)
  if ( ier1 /= 0 )then
     write(6,*)'problem allocating pre_exclude_list'
     stop
  endif
  allocate(pre_mask(siz),stat=ier1)
  if ( ier1 /= 0 )then
     write(6,*)'problem allocating pre_mask'
     stop
  endif
 !2nd pass . fill in
  num_pre_exclude = 0
  pre_mask = 0
  ! start with polar group
  do n = 1,num_atoms
    off1 = offset_polar_group(n)
    off2 = off_pre_exclude(n)
    do j = 1,num_in_polar_group(n)
      pre_exclude_list(off2+j) = nghb_polargroup(off1+j)
      pre_mask(off2+j) = ibset(pre_mask(off2+j),bit_polgrp)
    enddo
    num_pre_exclude(n) = num_pre_exclude(n) + num_in_polar_group(n)
  enddo
  ! next 1-2
  do n = 1,num_12
    j = nghb_12(1,n)
    k = nghb_12(2,n)
    offj = off_pre_exclude(j)
    offk = off_pre_exclude(k)
    num_pre_exclude(j) = num_pre_exclude(j) + 1
    num_pre_exclude(k) = num_pre_exclude(k) + 1
    pre_exclude_list(offj+num_pre_exclude(j)) = k
    pre_exclude_list(offk+num_pre_exclude(k)) = j
    pre_mask(offj+num_pre_exclude(j)) =  &
             ibset(pre_mask(offj+num_pre_exclude(j)),bit_12)
    pre_mask(offk+num_pre_exclude(k)) =  &
             ibset(pre_mask(offk+num_pre_exclude(k)),bit_12)
  enddo 
  ! next 1-3
  do n = 1,num_13
    j = nghb_13(1,n)
    k = nghb_13(2,n)
    offj = off_pre_exclude(j)
    offk = off_pre_exclude(k)
    num_pre_exclude(j) = num_pre_exclude(j) + 1
    num_pre_exclude(k) = num_pre_exclude(k) + 1
    pre_exclude_list(offj+num_pre_exclude(j)) = k
    pre_exclude_list(offk+num_pre_exclude(k)) = j
    pre_mask(offj+num_pre_exclude(j)) =  &
             ibset(pre_mask(offj+num_pre_exclude(j)),bit_13)
    pre_mask(offk+num_pre_exclude(k)) =  &
             ibset(pre_mask(offk+num_pre_exclude(k)),bit_13)
  enddo 
  ! next 1-4
  do n = 1,num_14
    j = nghb_14(1,n)
    k = nghb_14(2,n)
    offj = off_pre_exclude(j)
    offk = off_pre_exclude(k)
    num_pre_exclude(j) = num_pre_exclude(j) + 1
    num_pre_exclude(k) = num_pre_exclude(k) + 1
    pre_exclude_list(offj+num_pre_exclude(j)) = k
    pre_exclude_list(offk+num_pre_exclude(k)) = j
    pre_mask(offj+num_pre_exclude(j)) =  &
             ibset(pre_mask(offj+num_pre_exclude(j)),bit_14)
    pre_mask(offk+num_pre_exclude(k)) =  &
             ibset(pre_mask(offk+num_pre_exclude(k)),bit_14)
  enddo 
  ! next 1-5
  do n = 1,num_15
    j = nghb_15(1,n)
    k = nghb_15(2,n)
    offj = off_pre_exclude(j)
    offk = off_pre_exclude(k)
    num_pre_exclude(j) = num_pre_exclude(j) + 1
    num_pre_exclude(k) = num_pre_exclude(k) + 1
    pre_exclude_list(offj+num_pre_exclude(j)) = k
    pre_exclude_list(offk+num_pre_exclude(k)) = j
    pre_mask(offj+num_pre_exclude(j)) =  &
             ibset(pre_mask(offj+num_pre_exclude(j)),bit_15)
    pre_mask(offk+num_pre_exclude(k)) =  &
             ibset(pre_mask(offk+num_pre_exclude(k)),bit_15)
  enddo 
  ! now for the real excluded---remove redundencies
  allocate(num_excluded(num_atoms),stat=ier1)
  allocate(off_excluded(num_atoms),stat=ier2)
  if ( (ier1/=0)  .or. (ier2/=0) )then
     write(6,*)'problem allocating exclude_list'
     stop
  endif
  ! first pass--count
  do n = 1,num_atoms
    off1 = off_pre_exclude(n)
    num = num_pre_exclude(n)
    maxflag = n
    do j = 1,num
      k = pre_exclude_list(off1+j)
      if ( k > n )then
         flag(k) = n
         if ( k > maxflag )maxflag = k
      endif
    enddo
    num_ex = 0
    do j = n+1,maxflag
      if ( flag(j) == n )num_ex = num_ex + 1
    enddo
    num_excluded(n) = num_ex
  enddo
  off_excluded(1) = 0
  do n = 2,num_atoms
    off_excluded(n) = off_excluded(n-1) + num_excluded(n-1)
  enddo
  size_excluded_list = off_excluded(num_atoms) + num_excluded(num_atoms)
  allocate(excluded_atomlist(size_excluded_list),stat=ier1)
  if ( ier1 /= 0 )then
    write(6,*)'problems with exclusion list'
    stop
  endif
  allocate(adjust_list(3,size_excluded_list),stat=ier1)
  if ( ier1 /= 0 )then
    write(6,*)'problems with exclusion list'
    stop
  endif
  ! 2nd pass--fill list
  do n = 1,num_atoms
    off1 = off_pre_exclude(n)
    num = num_pre_exclude(n)
    maxflag = n
    do j = 1,num
      k = pre_exclude_list(off1+j)
      if ( k > n )then
         flag(k) = n
         mask_flag(k) = 0
         if ( k > maxflag )maxflag = k
      endif
    enddo
    ! repeat to fill the mask_flag
    maxflag = n
    do j = 1,num
      k = pre_exclude_list(off1+j)
      if ( k > n )then
         flag(k) = n
         mask_flag(k) = ior(mask_flag(k),pre_mask(off1+j))
         if ( k > maxflag )maxflag = k
      endif
    enddo
    num_ex = 0
    off2 = off_excluded(n)
    do j = n+1,maxflag
      if ( flag(j) == n )then
         num_ex = num_ex + 1
         excluded_atomlist(off2+num_ex) = j
         adjust_list(1,off2+num_ex) = n
         adjust_list(2,off2+num_ex) = j
         adjust_list(3,off2+num_ex) = mask_flag(j)
      endif
    enddo
    num_excluded(n) = num_ex
  enddo
  ! compress the values of adjust_list(3,n)
  ! check for sensible values
  do n = 1,size_excluded_list
    if ( adjust_list(3,n) == 2 )then 
      adjust_list(3,n) = 1 
    elseif ( adjust_list(3,n) == 4 )then 
      adjust_list(3,n) = 2 
    elseif ( adjust_list(3,n) == 8 )then 
      adjust_list(3,n) = 3 
    elseif ( adjust_list(3,n) == 16 )then 
      adjust_list(3,n) = 4 
    elseif ( adjust_list(3,n) == 3 )then 
      adjust_list(3,n) = 5 
    elseif ( adjust_list(3,n) == 5 )then 
      adjust_list(3,n) = 6 
    elseif ( adjust_list(3,n) == 9 )then 
      adjust_list(3,n) = 7 
    elseif ( adjust_list(3,n) == 17 )then 
      adjust_list(3,n) = 8 
    elseif ( adjust_list(3,n) == 1 )then 
      adjust_list(3,n) = 9 
    else
      ! only one of bit_12 through bit_15 should be set
      ! thus the above are the only sensible possibilities
      write(6,*)'adjust list: bad value n,adjust_list(k,n) = ', &
               n,adjust_list(1,n),adjust_list(2,n),adjust_list(3,n)
    endif
  enddo
  deallocate(num_pre_exclude)
  deallocate(off_pre_exclude)
  deallocate(pre_exclude_list)
  deallocate(pre_mask)
  deallocate(flag)
  deallocate(mask_flag)
  deallocate(off_excluded)
end subroutine AM_PARM_get_excluded_atom_list
!------------------------------------------------------------------
subroutine AM_PARM_find_trig_4th_atom(trigonal_angle)
  integer,intent(inout) :: trigonal_angle(4)

  integer n,k1,k2,j,k
  k1 = trigonal_angle(1)
  n = trigonal_angle(2)
  k2 = trigonal_angle(3)
  if ( atomvalence(n) /= 3 )then
    write(6,*)'trigonal angle screwup: ',n
    stop
  endif
  do j = 1,3
    k = bond_nghbors(offset_nghbors(n)+j)
    if ( (k /= k1) .and. (k /= k2) )then
      trigonal_angle(4) = k
    endif
  enddo
end subroutine AM_PARM_find_trig_4th_atom
!------------------------------------------------------------------
subroutine AM_PARM_fillout_pitorsion_6(atom)
  integer,intent(inout) :: atom(6)

  integer k1,k2,k3

  if ( atomvalence(atom(3)) /= 3 )then
    write(6,*)'pitorsion wrong valence for ',atom(3)
    stop
  endif
  k1 = bond_nghbors(offset_nghbors(atom(3))+1)
  k2 = bond_nghbors(offset_nghbors(atom(3))+2)
  k3 = bond_nghbors(offset_nghbors(atom(3))+3)
  if ( k1 == atom(4) )then
    atom(1) = k2
    atom(2) = k3
  elseif ( k2 == atom(4) )then
    atom(1) = k1
    atom(2) = k3
  elseif ( k3 == atom(4) )then
    atom(1) = k1
    atom(2) = k2
  else
    write(6,*)'pitorsion: problem with ',atom(3),atom(4)
    stop
  endif
  if ( atomvalence(atom(4)) /= 3 )then
    write(6,*)'pitorsion wrong valence for ',atom(4)
    stop
  endif
  k1 = bond_nghbors(offset_nghbors(atom(4))+1)
  k2 = bond_nghbors(offset_nghbors(atom(4))+2)
  k3 = bond_nghbors(offset_nghbors(atom(4))+3)
  if ( k1 == atom(3) )then
    atom(5) = k2
    atom(6) = k3
  elseif ( k2 == atom(3) )then
    atom(5) = k1
    atom(6) = k3
  elseif ( k3 == atom(3) )then
    atom(5) = k1
    atom(6) = k2
  else
    write(6,*)'pitorsion: problem with ',atom(3),atom(4)
    stop
  endif

end subroutine AM_PARM_fillout_pitorsion_6
!------------------------------------------------------------------
subroutine AM_PARM_get_next_tortor_table(frcfield_unit,success,  &
                                         atomtypes,dim1,dim2)
  integer,intent(in) :: frcfield_unit
  integer,intent(out) :: success,atomtypes(5),dim1,dim2

  character(len=120) :: line
  integer :: ios,j
  do
    read(frcfield_unit,'(A)',iostat=ios)line
    if ( ios /= 0 )then
      success = 0
      return
    endif
    success = index(line,'tortors')
    if ( success > 0 )then
      read(line(success+len('tortors'):),*)(atomtypes(j),j=1,5),dim1,dim2
      return
    endif
  enddo
end subroutine AM_PARM_get_next_tortor_table
!------------------------------------------------------------------
subroutine AM_PARM_num_fields(string,num_fields)
  integer,intent(out) :: num_fields
  character(len=*),intent(in) :: string

  character(len=256) :: string1,string2
  integer :: num,ptr,lstr,lstr1

  lstr1 = len_trim(string)
  if ( lstr1 == 0 )then
    num_fields = 0
    return
  endif
  string1(1:lstr1) = string(1:lstr1)
  num = 0
  do
    string2(1:lstr1) = adjustl(string1(1:lstr1))
    lstr = len_trim(string2(1:lstr1))
    if ( lstr > 0 )then
      num = num+1
      ptr = index(string2(1:lstr),' ') 
      if ( ptr == 0 )exit
      lstr1 = lstr-ptr+1
      string1(1:lstr1) = string2(ptr:lstr)
    else
      exit
    endif
  enddo
  num_fields = num
end subroutine AM_PARM_num_fields
!------------------------------------------------------------------
subroutine AM_PARM_deallocate()

  integer :: n
  if ( allocated(atomtype) )deallocate(atomtype)
  if ( allocated(atomclass) )deallocate(atomclass)
  if ( allocated(atomic_number) )deallocate(atomic_number)
  if ( allocated(atomvalence) )deallocate(atomvalence)
  if ( allocated(atomic_weight) )deallocate(atomic_weight)
  if ( allocated(bond_list) )deallocate(bond_list)
  if ( allocated(bond_param) )deallocate(bond_param)
  if ( allocated(urey_bradley_list) )deallocate(urey_bradley_list)
  if ( allocated(urey_bradley_param) )deallocate(urey_bradley_param)
  if ( allocated(regular_angle_list) )deallocate(regular_angle_list)
  if ( allocated(regular_angle_param) )deallocate(regular_angle_param)
  if ( allocated(trigonal_angle_list) )deallocate(trigonal_angle_list)
  if ( allocated(trigonal_angle_param) )deallocate(trigonal_angle_param)
  if ( allocated(outofplane_list) )deallocate(outofplane_list)
  if ( allocated(outofplane_param) )deallocate(outofplane_param)
  if ( allocated(torsion_list) )deallocate(torsion_list)
  if ( allocated(torsion_param) )deallocate(torsion_param)
  if ( allocated(str_torsion_list) )deallocate(str_torsion_list)
  if ( allocated(str_torsion_param) )deallocate(str_torsion_param)
  if ( allocated(pitorsion_list) )deallocate(pitorsion_list)
  if ( allocated(pitorsion_param) )deallocate(pitorsion_param)
  if ( allocated(stretch_bend_list) )deallocate(stretch_bend_list)
  if ( allocated(stretch_bend_param) )deallocate(stretch_bend_param)
  if ( allocated(tortor_list) )deallocate(tortor_list)
  if ( allocated(tortor_table) )then
    do n = 1,num_tortor_tables
      deallocate(tortor_table(n)%angle1)
      deallocate(tortor_table(n)%angle2)
      deallocate(tortor_table(n)%func)
      deallocate(tortor_table(n)%dfunc_dangle1)
      deallocate(tortor_table(n)%dfunc_dangle2)
      deallocate(tortor_table(n)%d2func_dangle1_dangle2)
    enddo
    deallocate(tortor_table)
  endif
  if ( allocated(bond_nghbors) )deallocate(bond_nghbors)
  if ( allocated(offset_nghbors) )deallocate(offset_nghbors)
  if ( allocated(vdw_atom_parent) )deallocate(vdw_atom_parent)
  if ( allocated(xyz_crd) )deallocate(xyz_crd)
  if ( allocated(atomname) )deallocate(atomname)
  if ( allocated(reslabel) )deallocate(reslabel)
  if ( allocated(startres) )deallocate(startres)
  if ( allocated(atomsymbol) )deallocate(atomsymbol)
  if ( allocated(num_atoms_in_molecule) )deallocate(num_atoms_in_molecule)
  if ( allocated(polarizability) )deallocate(polarizability)
  if ( allocated(num_in_polar_group) )deallocate(num_in_polar_group)
  if ( allocated(nghb_polargroup) )deallocate(nghb_polargroup)
  if ( allocated(offset_polar_group) )deallocate(offset_polar_group)
  if ( allocated(nghb_12) )deallocate(nghb_12)
  if ( allocated(nghb_13) )deallocate(nghb_13)
  if ( allocated(nghb_14) )deallocate(nghb_14)
  if ( allocated(nghb_15) )deallocate(nghb_15)
  if ( allocated(num_excluded) )deallocate(num_excluded)
  if ( allocated(excluded_atomlist) )deallocate(excluded_atomlist)
  if ( allocated(adjust_list) )deallocate(adjust_list)
  if ( allocated(frame_list) )deallocate(frame_list)
  if ( allocated(chiral_frame_list) )deallocate(chiral_frame_list)
  if ( allocated(local_multipole) )deallocate(local_multipole)
  if ( allocated(vdw_atom_parent) )deallocate(vdw_atom_parent)
  if ( allocated(vdw_parent_weight) )deallocate(vdw_parent_weight)
  if ( allocated(atom_vdw_type) )deallocate(atom_vdw_type)
  if ( allocated(vdw_mix_radius) )deallocate(vdw_mix_radius)
  if ( allocated(vdw_mix_epsilon) )deallocate(vdw_mix_epsilon)
  if ( allocated(dyn_crd) )deallocate(dyn_crd)
  if ( allocated(dyn_vel) )deallocate(dyn_vel)
  if ( allocated(dyn_accel) )deallocate(dyn_accel)
  if ( allocated(dyn_old_accel) )deallocate(dyn_old_accel)
end subroutine AM_PARM_deallocate
!------------------------------------------------------------------
subroutine AM_PARM_process_key_file(key_unit,frcfieldfile)
  integer,intent(in) :: key_unit
  character(len=*) :: frcfieldfile

  double precision :: boxmax
  character(len=120) :: line,lline,word,tempf
  integer :: ios,success
  rewind(key_unit)
  success = 0
  word = ' '
  do
    read(key_unit,'(A)',iostat=ios)line
    if ( ios /= 0 )then
      write(6,*)'failed to get forcefield file from key file!'
      stop
    endif
    lline = adjustl(line)
    success = index(lline,'parameters')
    if ( success == 1 )then
      success = index(line(success:),' ')
      if ( success > 0 )then
        read(line(success:),'(a)')tempf
        tempf = adjustl(tempf)
        frcfieldfile = tempf(1:len_trim(tempf)) // '.prm'
        write(6,*)'frcfieldfile = ',frcfieldfile
        exit
      else
        write(6,*)'failed to get forcefield file from key file!'
        stop
      endif
    endif
  enddo
  rewind(key_unit)
  do
    read(key_unit,'(A)',iostat=ios)line
    if ( ios /= 0 )then
      exit
    endif
    if ( len_trim(line) > 0 )then
       read(line,*)word
       if ( word(1:6) == 'a-axis' )then
         success = index(line,'a-axis')
         success = index(line(success:),' ')
         if ( success > 0 )read(line(success:),*)boxa
         write(6,*)'found a-axis = ',boxa
       elseif( word(1:6) == 'b-axis' )then
         success = index(line,'b-axis')
         success = index(line(success:),' ')
         if ( success > 0 )read(line(success:),*)boxb
         write(6,*)'found b-axis = ',boxb
       elseif( word(1:6) == 'c-axis' )then
         success = index(line,'c-axis')
         success = index(line(success:),' ')
         if ( success > 0 )read(line(success:),*)boxc
         write(6,*)'found c-axis = ',boxc
       elseif( word(1:5) == 'alpha' )then
         success = index(line,'alpha')
         success = index(line(success:),' ')
         if ( success > 0 )read(line(success:),*)box_alpha
         write(6,*)'found alpha = ',box_alpha
       elseif( word(1:4) == 'beta' )then
         success = index(line,'beta')
         success = index(line(success:),' ')
         if ( success > 0 )read(line(success:),*)box_beta
         write(6,*)'found beta = ',box_beta
       elseif( word(1:5) == 'gamma' )then
         success = index(line,'gamma')
         success = index(line(success:),' ')
         if ( success > 0 )read(line(success:),*)box_gamma
         write(6,*)'found gamma = ',box_gamma
       endif
    endif
  enddo
  boxmax = max(boxa,boxb,boxc)
  ! This is no longer fatal, since we can read the box information from the XYZ
  ! file instead. Any box info in the XYZ file will overwrite the box info read
  ! from the .key file
! if ( boxmax == 0.d0 )then
!   write(6,*)'failed to get unitcell params from key file'
!   stop
! endif
  if ( boxa == 0.d0 )boxa = boxmax
  if ( boxb == 0.d0 )boxb = boxmax
  if ( boxc == 0.d0 )boxc = boxmax
  if ( box_alpha == 0.d0 )box_alpha = 90.d0
  if ( box_beta == 0.d0 )box_beta = 90.d0
  if ( box_gamma == 0.d0 )box_gamma = 90.d0
end subroutine AM_PARM_process_key_file

subroutine get_num_tokens(string, token_num)

  implicit none

! Passed arguments

  character(*), intent(in) :: string

  integer, intent(out)     :: token_num

! Local variables

  integer :: string_loc  ! our location in the string
  integer :: iend        ! last non-whitespace character location

  string_loc = 1
  iend = len_trim(string)
  token_num = 0

  do while (string_loc .le. iend)

    if ( string(string_loc:string_loc) .le. ' ' ) then
      string_loc = string_loc + 1
    else

      do while ( string(string_loc:string_loc) .gt. ' ' )
        string_loc = string_loc + 1
      end do

      token_num = token_num + 1
    end if
  end do

end subroutine get_num_tokens
!------------------------------------------------------------------
end module amoeba_parm
!------------------------------------------------------------------
program amoeba_parm_main
  use amoeba_parm
  implicit none
  integer :: argc,arg,IArgC
  character(len=256) :: argv
  character(len=80)  :: analoutfile,prmtopfile,frcfieldfile,key_file, &
                        xyz_file,dyn_file,pdb_file,inpcrd_file,filehead
  integer :: analout_unit,prmtop_unit,frcfield_unit,xyz_unit,pdb_unit, &
             dyn_unit,inpcrd_unit,key_unit
  integer j,ios
  logical dyn_exists

  analoutfile = ''
  prmtopfile = ''
  frcfieldfile = ''
  xyz_file = ''
  dyn_file = ''
  pdb_file = ''
  title = ''
  inpcrd_file = ''
  filehead = ''
  key_file = ''
  arg = 1
  argc = IArgC()
  do while (arg <= argc)
    call GetArg(arg,argv)
    if ( argv == '-analout')then
      arg=arg+1
      call GetArg(arg,analoutfile)
    elseif ( argv == '-prmtop')then
      arg=arg+1
      call GetArg(arg,prmtopfile)
    elseif ( argv == '-xyz')then
      arg=arg+1
      call GetArg(arg,xyz_file)
    elseif ( argv == '-dyn')then
      arg=arg+1
      call GetArg(arg,dyn_file)
    elseif ( argv == '-pdb')then
      arg=arg+1
      call GetArg(arg,pdb_file)
    elseif ( argv == '-inpcrd')then
      arg=arg+1
      call GetArg(arg,inpcrd_file)
    elseif ( argv == '-title')then
      arg=arg+1
      call GetArg(arg,title)
    elseif ( argv == '-name')then
      arg=arg+1
      call GetArg(arg,filehead)
    elseif ( argv == '-key')then
      arg=arg+1
      call GetArg(arg,key_file)
    elseif ( argv == '-precise')then
      arg=arg+1
      call AM_PARM_set_real_format()
    endif
    arg = arg + 1
  enddo
  if ( filehead /= '' )then
    ! assign some unassigned using filehead
    if ( analoutfile == '' )analoutfile = trim(filehead)//'.analout'
    if ( xyz_file == '' )   xyz_file = trim(filehead)//'.xyz'
    if ( dyn_file == '' )   dyn_file = trim(filehead)//'.dyn'
    if ( pdb_file == '' )   pdb_file = trim(filehead)//'.pdb'
    if ( prmtopfile == '' ) prmtopfile = trim(filehead)//'.prmtop'
    if ( inpcrd_file == '' )inpcrd_file = trim(filehead)//'.inpcrd'
    if ( key_file == '' )   key_file = trim(filehead)//'.key'
  endif
  if ( analoutfile == '' )then
    write(6,*)'analoutfile: '
    read(5,*)analoutfile
  endif
  if ( xyz_file == '' )then
    write(6,*)'xyz_file: '
    read(5,*)xyz_file
  endif
  if ( pdb_file == '' )then
    write(6,*)'pdb_file:'
    read(5,*)pdb_file
  endif
  if ( title == '' )then
    write(6,*)'title: '
    read(5,'(a)')title
  endif
  if ( prmtopfile == '' )then
    write(6,*)'prmtopfile (output): '
    read(5,*)prmtopfile
  endif
  if ( inpcrd_file == '' )then
    write(6,*)'inpcrd_file (output): '
    read(5,*)inpcrd_file
  endif
  key_unit = 7
  open(unit=key_unit,file=key_file,status='old',iostat=ios)
  if ( ios /= 0 )then
    write(6,*)'trouble opening key_file ',key_file
    stop
  endif
  call AM_PARM_process_key_file(key_unit,frcfieldfile)
  analout_unit=8
  open(unit=analout_unit,file=analoutfile,status='old',iostat=ios)
  if ( ios /= 0 )then
    write(6,*)'trouble opening analoutfile ',analoutfile
    stop
  endif
  prmtop_unit=9
  open(unit=prmtop_unit,file=prmtopfile,status='new',iostat=ios)
  if ( ios /= 0 )then
    write(6,*)'trouble opening prmtopfile ',prmtopfile
    stop
  endif
  frcfield_unit=10
  open(unit=frcfield_unit,file=frcfieldfile,status='old',iostat=ios)
  if ( ios /= 0 )then
    write(6,*)'trouble opening frcfieldfile ',frcfieldfile
    stop
  endif
  xyz_unit=11
  open(unit=xyz_unit,file=xyz_file,status='old',iostat=ios)
  if ( ios /= 0 )then
    write(6,*)'trouble opening xyz_file ',xyz_file
    stop
  endif
  pdb_unit=12
  open(unit=pdb_unit,file=pdb_file,status='old',iostat=ios)
  if ( ios /= 0 )then
    write(6,*)'trouble opening pdb_file ',pdb_file
    stop
  endif
  inpcrd_unit=13
  open(unit =inpcrd_unit,file=inpcrd_file,status='new',iostat=ios)
  if ( ios /= 0 )then
    write(6,*)'trouble opening inpcrd_unit ',inpcrd_unit
    stop
  endif
  inquire(file=dyn_file,exist=dyn_exists)
  
  call AM_PARM_get_xyz(xyz_unit)
  call AM_PARM_get_pdb(pdb_unit)
  if ( dyn_exists )then
     dyn_unit=14
     open(unit =dyn_unit,file=dyn_file,status='old',iostat=ios)
     if ( ios /= 0 )then
        write(6,*)'trouble opening dyn_unit ',dyn_unit
        stop
     endif
     call AM_PARM_get_dyn(dyn_unit)
  endif
  call AM_PARM_get_atom_params(analout_unit)
  call AM_PARM_check_atom_num()
  call AM_PARM_get_atom_vdw_params(analout_unit)
  call AM_PARM_get_bond_params(analout_unit)
  call AM_PARM_get_UreyB_params(analout_unit)
  call AM_PARM_bond_nghb_list()
  call AM_PARM_find_vdw_parent()
  call AM_PARM_get_angle_params(analout_unit)
  call AM_PARM_get_outofplane_params(analout_unit)
  call AM_PARM_get_torsion_params(analout_unit)
  call AM_PARM_get_stretch_tor_params(analout_unit)
  call AM_PARM_get_pitorsion_params(analout_unit)
  call AM_PARM_get_strbend_params(analout_unit)
  call AM_PARM_get_tortor_params(frcfield_unit,analout_unit)
  call AM_PARM_get_mol_list()
  call AM_PARM_get_multipole(analout_unit)
  call AM_PARM_get_polar_params(analout_unit)
  call AM_PARM_get_nghb_12131415_lists(analout_unit)
  call AM_PARM_get_excluded_atom_list(analout_unit)
  ! write out params
  call AM_PARM_write_inpcrd(inpcrd_unit)
  call AM_PARM_write_top(prmtop_unit)
  call AM_PARM_write_bond_params(prmtop_unit)
  call AM_PARM_write_UreyB_params(prmtop_unit)
  call AM_PARM_write_angle_params(prmtop_unit)
  call AM_PARM_write_outofplane_params(prmtop_unit)
  call AM_PARM_write_torsion_params(prmtop_unit)
  call AM_PARM_write_str_tor_params(prmtop_unit)
  call AM_PARM_write_pitorsion_params(prmtop_unit)
  call AM_PARM_write_strbend_params(prmtop_unit)
  call AM_PARM_write_tortor_params(prmtop_unit)
  call AM_PARM_write_vdw(prmtop_unit)
  call AM_PARM_write_multipole(prmtop_unit)
  call AM_PARM_write_adjust(prmtop_unit)
  call AM_PARM_write_polar(prmtop_unit)
  call AM_PARM_deallocate()

end program amoeba_parm_main
