! <compile=optimized>
#include "../include/assert.fh"

! This is a module to encapsulate all of the AMBER common blocks and
! other global junk into a module namespace.
module xray_common_module
use file_io_dat
end module xray_common_module

module xray_interface_module
   ! Calls from main SANDER code:
   !  dynlib.f:   if (xray_active) call xray_write_md_state(stdout=6)
   !  printe.f:   if (xray_active) call xray_write_min_state(stdout=6)
   !  force.f:    if (xray_active) call xray_get_derivative(coor=x,deriv=f)
   !  mdread.f:   call xray_read_mdin(mdin_lun=5) !<== Also does the global init
   !  sander.f:   call xray_read_parm(prmtop_lun=8,stdout=6)
   !  sander.f:   call xray_init()
   !  sander.f:   call xray_fini()

   use xray_globals_module
   implicit none
   private

   namelist /xray/ &
         pdb_infile, pdb_outfile, &
         pdb_read_coordinates, &
         pdb_use_segid, &
         pdb_wrap_names, &
         spacegroup_name, &
         reflection_infile, &
         reflection_infile_format, &
         reflection_fobs, &
         reflection_sigfobs, &
         resolution_low, &
         resolution_high, &
         xray_weight, &
         solvent_mask_probe_radius, &
         solvent_mask_expand, &
         solvent_mask_outfile, &
         solvent_mask_reflection_outfile, &
         solvent_mask_update_interval, &
         solvent_scale, &
         solvent_bfactor, &
         fft_method, &
         fft_grid_size, &
         fft_grid_spacing, &
         fft_bfactor_sharpen, &
         fft_density_tolerance, &
         fft_reflection_tolerance, &
         fft_radius_min, fft_radius_max, &
         bfactor_min, bfactor_max, &
         bfactor_refinement_interval, &
         atom_selection_mask

   ! Common public entities all have an xray_ prefix.
   ! Others assume localization by the module scope.
   public :: xray_init_globals, xray_init, xray_fini
   public :: xray_active, xray_get_derivative, xray_write_md_state
   public :: xray_read_mdin, xray_read_parm, xray_read_pdb
   public :: xray_write_options, xray_write_min_state

   !-------------------------------------------------------------------
contains

   subroutine xray_read_mdin(mdin_lun)
      implicit none
      integer, intent(in) :: mdin_lun
      integer :: stat
      call xray_init_globals()
      if (.not.xray_active) return
      rewind(mdin_lun)
      read(unit=mdin_lun,nml=xray,iostat=stat)
      if (stat /= 0) then
         write(stdout,'(A)') 'Error reading namelist &xray.'
         call mexit(stdout,1)
      end if
      !write(unit=6,nml=xray)
   end subroutine xray_read_mdin

   subroutine xray_write_options()
      write(stdout,'(/,A)') 'X-ray Refinement Parameters:'
      write(stdout,'(5X,2A)') 'PDB InFile: ',trim(pdb_infile)
      write(stdout,'(5X,2A)') 'PDB OutFile:',trim(pdb_outfile)
      write(stdout,'(5X,A,L1)') 'PDB Read Coordinates: ',pdb_read_coordinates
      write(stdout,'(5X,A,L1)') 'PDB Use SegID: ',pdb_use_segid
      write(stdout,'(5X,A,L1)') 'PDB Wrap Names: ',pdb_wrap_names
      write(stdout,'(5X,2A)') 'Spacegroup: ',trim(spacegroup_name)
      write(stdout,'(5X,2A)') 'Reflection InFile: ',trim(reflection_infile)
      !write(stdout,'(5X,A)') reflection_infile_format
      !write(stdout,'(5X,A)') reflection_fobs
      !write(stdout,'(5X,A)') reflection_sigfobs
      write(stdout,'(5X,2(A,F8.3))') 'Resolution Range: ',resolution_low,',',resolution_high
      write(stdout,'(5X,A,E10.3)') 'X-ray weight: ',xray_weight
      write(stdout,'(5X,A,F8.3)') 'Solvent mask probe radius: ',solvent_mask_probe_radius
      write(stdout,'(5X,A,F8.3)') 'Solvent mask expand: ',solvent_mask_expand
      write(stdout,'(5X,2A)') 'Solvent Mask OutFile:',trim(solvent_mask_outfile)
      write(stdout,'(5X,2A)') 'Solvent Mask Reflection OutFile:',trim(solvent_mask_reflection_outfile)
      write(stdout,'(5X,A,I4)') 'Solvent Mask Update Interval: ',solvent_mask_update_interval
      write(stdout,'(5X,2(A,F8.3))') 'Solvent scale:',solvent_scale,', B-factor:', solvent_bfactor
      write(stdout,'(5X,A,I2)')   'FFT method: ',fft_method
      write(stdout,'(5X,A,3(5X,I5))') 'FFT Grid Size: ',fft_grid_size
      write(stdout,'(5X,A,F9.5)') 'FFT Grid Spaceing: ',fft_grid_spacing
      write(stdout,'(5X,A,F8.3)') 'FFT B-factor Sharpen: ',fft_bfactor_sharpen
      write(stdout,'(5X,A,E10.3)') 'FFT Densty Toleranec: ',fft_density_tolerance
      write(stdout,'(5X,A,E10.3)') 'FFT Reflection Tolerance: ',fft_reflection_tolerance
      write(stdout,'(5X,2(A,F8.3))') 'FFT Radius Min:',fft_radius_min,', Max: ',fft_radius_max
      write(stdout,'(5X,2(A,F8.3))') 'B-Factor Min:',bfactor_min,', Max: ',bfactor_max
      write(stdout,'(5X,A,I4)') 'B-factor Refinement Interval: ',bfactor_refinement_interval
      write(stdout,'(5X,2A)') 'Atom Selection Mask: ',trim(atom_selection_mask)
   end subroutine xray_write_options

   ! Read X-ray data from the PRMTOP file, and also save pointers to global
   ! PRMTOP data.
   subroutine xray_read_parm(prmtop_lun,out_lun)
      use memory_module, only: natom, nres
      use xray_reciprocal_space_module, only: SYMM_TRICLINIC
      implicit none
      integer, intent(in) :: prmtop_lun, out_lun
      ! local
      character(len=32) :: fmt
      integer :: alloc_status, ierr

      if (pdb_outfile /= '') then
         num_atoms = natom
         num_residues = nres

         allocate(atom_bfactor(natom), atom_occupancy(natom), &
               atom_selection(natom), residue_chainid(nres), residue_icode(nres), &
               atom_element(natom), atom_altloc(natom), residue_number(nres), &
               stat=alloc_status)
         REQUIRE(alloc_status==0)

         call nxtsec(prmtop_lun,STDOUT,0,'*','RESIDUE_NUMBER',fmt,ierr)
         read(prmtop_lun,fmt) residue_number
         call nxtsec(prmtop_lun,STDOUT,0,'*','RESIDUE_CHAINID',fmt,ierr)
         read(prmtop_lun,fmt) residue_chainid

         call nxtsec(prmtop_lun,STDOUT,1,'*','RESIDUE_ICODE',fmt,ierr)
         if (ierr==0) then
            read(prmtop_lun,fmt) residue_icode
         else
            residue_icode=' '
         end if

         call nxtsec(prmtop_lun,STDOUT,1,'*','ATOM_ALTLOC',fmt,ierr)
         if (ierr==0) then
            read(prmtop_lun,fmt) atom_altloc
         else
            atom_altloc=' '
         end if

         call nxtsec(prmtop_lun,STDOUT,0,'*','ATOM_ELEMENT',fmt,ierr)
         read(prmtop_lun,fmt) atom_element
      end if

      if (reflection_infile == '') return

      call nxtsec(prmtop_lun,out_lun,0,'*','XRAY_NUM_SCATTER_TYPES',fmt,ierr)
      if (fmt=='*') then
         write(stdout,'(A)') &
            'ERROR: XRAY_NUM_SCATTER_TYPES not found in PRMTOP file.'
         call mexit(stdout,1)
      end if
      read(prmtop_lun,fmt) num_scatter_types

      allocate(atom_scatter_type(natom),  &
            scatter_coefficients(2,scatter_ncoeffs,num_scatter_types), &
            stat=alloc_status)
      REQUIRE(alloc_status==0)

      call nxtsec(prmtop_lun,out_lun,0,'*','XRAY_ATOM_SCATTER_TYPE_INDEX',fmt,ierr)
      read(prmtop_lun,fmt) atom_scatter_type
      call nxtsec(prmtop_lun,out_lun,0,'*','XRAY_SCATTER_COEFFICIENTS',fmt,ierr)
      read(prmtop_lun,fmt) scatter_coefficients
      call nxtsec(prmtop_lun,out_lun,1,'*','XRAY_SYMMETRY_TYPE',fmt,ierr)
      if (ierr==-2) then
         write(STDOUT,*) &
               'XRAY_SYMMETRY_TYPE not found in PRMTOP file; assuming P1'
         num_symmops = 1
         spacegroup_number = 1
         spacegroup_name = 'P 1'
         au_type = 1
         system = SYMM_TRICLINIC
         symmop(:,:,1) = reshape((/1,0,0, 0,1,0, 0,0,1, 0,0,0/),(/3,4/))
         symmop_inv(:,:,1) = symmop(:,:,1)
      else
         stop 'ONLY P1 SUPPORTED FOR NOW'
         read(prmtop_lun,fmt) num_symmops, spacegroup_number, au_type, system
         call nxtsec(prmtop_lun,out_lun,1,'*', &
               'XRAY_SYMMETRY_OPERATORS',fmt,ierr)
         ! ...
      end if
   end subroutine xray_read_parm

   subroutine xray_read_pdb(filename)
      use memory_module, only: residue_label,atom_name,coordinate
      use xray_utils_module, only: allocate_lun
      implicit none
      character(len=*), intent(in) :: filename
      ! locals
      character(len=4) :: name,resName,segID,element,altLoc,chainID,iCode
      integer :: serial,resSeq
      real(real_kind) :: xyz(3),occupancy,tempFactor
      character(len=80) :: line
      integer :: unit, iostat, iatom, ires, i, j, ndup, nmiss
      real(real_kind), parameter :: MISSING = -999.0_rk_
      ! begin
      atom_occupancy(:)=MISSING
      call amopen(allocate_lun(unit),filename,'O','F','R')
      ndup=0
      iatom=1
      ires=1
      do
         read(unit,'(A)',iostat=iostat) line
         if (iostat/=0) exit
         if (line(1:6)=='END   ') exit
         if (line(1:6)=='ATOM  ' .or. line(1:6)=='HETATM') then
            read(line,'(6X,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,6X,2A4)') &
                  serial,name,altLoc,resName,chainID,resSeq,iCode, &
                  xyz,occupancy,tempFactor,segID,element
            i = find_atom(name,resName,chainID,resSeq,iCode)
            if (i<0) then
               write(stdout,'(2A)') 'Atom not found: ',trim(line)
               stop
            end if
            if (atom_occupancy(i) >= 0) then
               ndup=ndup+1
               if (ndup<10) then
                  write(stdout,'(3(A,1X),A,I4,A)') 'PDB: Duplicate ATOM:', &
                        name,resName,chainID(1:1),resSeq,iCode(1:1)
               end if
            end if
            if (pdb_read_coordinates) coordinate(1:3,i) = xyz
            atom_bfactor(i) = tempFactor
            atom_occupancy(i) = occupancy
         end if
      end do
      nmiss = count(atom_occupancy==MISSING)
      if (nmiss>0) then
         write(stdout,'(A,I4,A)') 'PDB: missing data for ',nmiss,' atoms.'
         j=0
         do i=1,num_atoms
            if (atom_occupancy(i)==MISSING) then
               atom_occupancy(i)=0
               j=j+1
               if (j<=10) then
                  write(stdout,'(3(A,1X),A,I4,A)') 'PDB: Missing ATOM:', &
                        atom_name(i),residue_label(i),residue_chainID(i)(1:1),&
                        residue_number(i),residue_iCode(i)(1:1)
               end if
            end if
         end do
      end if
      if (nmiss==0 .and. ndup==0) then
         write(stdout,'(A)') 'PDB: All atoms read successfully.'
      end if
      close(unit)
      return
   end subroutine xray_read_pdb

   function find_atom(name,resName,chainID,resSeq,iCode) result(atom_serial)
      use memory_module, only: residue_pointer,residue_label,atom_name
      implicit none
      integer :: atom_serial
      character(len=4), intent(in) :: name, resName, chainID, iCode
      integer, intent(in) :: resSeq
      ! locals
      character(len=4) :: lname
      integer, save :: ires = 1
      integer :: i,j
      lname = adjustl(name)
      ! Unwrap PDB v2 hydrogen names
      if (name(1:1)>='0' .and. name(2:2)<='9') then
         lname = trim(lname(2:))//lname(1:1)
      end if
      ! first find the matching residue:
      do i=1,num_residues
         if (resSeq==residue_number(ires) &
               .and. chainID==residue_chainid(ires) &
               .and. iCode==residue_icode(ires) &
               .and. resName==residue_label(ires)) then
            ! then find the matching atom name:
            do j = residue_pointer(ires),residue_pointer(ires+1)-1
               if (lname==atom_name(j)) then
                  atom_serial = j
                  return
               end if
            end do
            ! Continue searching, just in case there is a residue
            ! that has been split into two parts.
         end if
         ires = ires + 1
      end do
      atom_serial = -1
      return
   end function find_atom

   subroutine xray_write_pdb(filename)
      use xray_common_module, only: owrite, title, title1
      use xray_utils_module, only: allocate_lun
      use memory_module, only: &
            residue_pointer,residue_label,atom_name,coordinate, &
            num_bonds=>nbona, &
            bond_atom1=>bonds_without_hydrogen_1, &
            bond_atom2=>bonds_without_hydrogen_2
      implicit none
      character(len=*), intent(in) :: filename
      ! locals
      integer :: unit, iatom, ires, ierr
      integer :: first1, last1, first2, last2, ibond
      character(len=4) :: name
      character(len=8) :: date
      character(len=10) :: time
      logical, allocatable :: linked(:)
      ! character(len=1) :: altloc
      ! character(len=4) :: segid
      character(len=3) :: resName
      logical          :: isStandardRes
      character(len=3), parameter :: standard_pdb_residues(28) = (/ &
            "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE", &
            "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL", &
            " DG"," DA"," DT"," DC","  G","  A","  U","  C" /)
      ! Amber modres types: "CYX","HID","HIE","HIP",
      character(len=*), parameter :: pdbfmt_MODRES = &
            '("MODRES",1X,A4,1X,A3,1X,A1,1X,I4,A1,1X,A3,2X,A41)'
      ! GMS: Fix for pgf90 compiler
      character(len=4) :: this_residue_chainid
      ! begin
      allocate(linked(num_residues),stat=ierr)
      if (ierr==0) then
         do ires = 1,num_residues-1
            first1 = 3*(residue_pointer(ires))
            last1  = 3*(residue_pointer(ires+1)-1)
            first2 = 3*(residue_pointer(ires+1))
            last2  = 3*(residue_pointer(ires+2)-1)
            do ibond=1,num_bonds
               if (bond_atom1(ibond) >= first1 &
                     .and. bond_atom1(ibond) <= last1 &
                     .and. bond_atom2(ibond) >= first2 &
                     .and. bond_atom2(ibond) <= last2) then
                  linked(ires)=.true.
                  exit
               end if
            end do
         end do
      end if

      call amopen(allocate_lun(unit),filename,owrite,'F','R')
      call date_and_time(date,time)
      if (title/='') write(unit,'(2A)') 'REMARK  ', title
      if (title1/='') write(unit,'(2A)') 'REMARK  ', title1
      write(unit,'(12A)') 'REMARK  Written by Amber 10, SANDER, ', &
            date(1:4),'.',date(5:6),'.',date(7:8),'  ', &
            time(1:2),':',time(3:4),':',time(5:6)

      ! Actually '(6A,3F9.3A9,3F7.2,1X,A11,I4)', with last value = Z
      write(unit,'(A6,3F9.3,3F7.2,1X,A11)') &
            'CRYST1',unit_cell,spacegroup_name
#if 0
      do ires = 1,num_residues
         if (residue_chainid(ires)=='*') cycle
         if (residue_label(ires)=='HID') then
            write(unit,pdbfmt_MODRES) &
                  '----','HID',
            residue_chainid(ires)(1:1), &
                  residue_number(ires),residue_icode(ires)(1:1), &
                  'HIS','HE2 ATOM REMOVED'
         else if (residue_label(ires)=='HIE') then
            write(unit,pdbfmt_MODRES) &
                  '----','HIE',
            residue_chainid(ires)(1:1), &
                  residue_number(ires),residue_icode(ires)(1:1), &
                  'HIS','HD1 ATOM REMOVED'
         else if (residue_label(ires)=='HIP') then
            write(unit,pdbfmt_MODRES) &
                  '----','HIP',
            residue_chainid(ires)(1:1), &
                  residue_number(ires),residue_icode(ires)(1:1), &
                  'HIS','HD1 AND HE2 ATOMS REMOVED'
         end if
      end do
#endif

      do ires = 1,num_residues
         if (residue_chainid(ires)=='*') cycle
         do iatom = residue_pointer(ires), residue_pointer(ires+1)-1
            ! ***NOTE***
            ! This code only adds a leading space to give element-alignment
            ! where possible. It is impossible to follow the PDB version 3
            ! "remediated" format correctly, because it has no alignement rules.
            ! Instead, it assumes you have a complete database of all known
            ! residues, and any other residue names are a fatal error.
            name = atom_name(iatom)
            if (atom_element(iatom)(1:1)==' ' &
                  .and. name(1:1) == atom_element(iatom)(2:2)) then
               if (len_trim(name) < 4 .or. pdb_wrap_names) then
                  name = name(4:4)//name(1:3)
               end if
            end if
            resName=residue_label(ires)(1:3)
            resName=adjustr(resName)
            ! GMS: Fix for pgf90 compiler
            this_residue_chainid = residue_chainid(ires)
            ! DRR: PGI does not seem to like any() inside merge() intrinsic.
            isStandardRes = any(resName==standard_pdb_residues)
            write(unit,'(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,6X,2A4)')&
                  merge('ATOM  ', 'HETATM', isStandardRes), &
                  iatom,name,atom_altloc(iatom)(1:1), &
                  resName,residue_chainid(ires)(1:1), &
                  residue_number(ires),residue_icode(ires)(1:1), &
                  coordinate(1:3,iatom), &
                  atom_occupancy(iatom), &
                  atom_bfactor(iatom), &
                  ! GMS: Fix for pgf90 compiler
                  !merge(residue_chainid(ires),'    ',pdb_use_segid), &
                  merge(this_residue_chainid,'    ',pdb_use_segid), &
                  atom_element(iatom)
         end do
      end do
      write(unit,'(A)') 'END'
      close(unit)
      return
   end subroutine xray_write_pdb

   subroutine xray_init()
      use nblist, only: a, b, c, alpha, beta, gamma
      use xray_utils_module, only: allocate_lun
      use xray_reciprocal_space_module, only: derive_cell_info
      !use xray_FFT_module, only: FFT_setup, get_mss4
      use xray_fourier_module, only: get_mss4
      use findmask, only: atommask
      use memory_module, only: natom,nres,ih,m02,m04,m06,ix,i02,x,lcrd
      implicit none
      ! local
      integer :: hkl_lun, i, alloc_status

      if (pdb_infile /= '') call xray_read_pdb(trim(pdb_infile))

      if (reflection_infile == '') xray_active = .false.

      if (.not.xray_active) then
         unit_cell = (/a, b, c, alpha, beta, gamma/)
         spacegroup_name = 'P 1'
         return
      end if

      write(stdout,'(A,3F9.3,3F7.2)') &
            'XRAY: UNIT CELL= ',a, b, c, alpha, beta, gamma
      call derive_cell_info(a, b, c, alpha, beta, gamma)
      ! Ewald/X-ray equivalences:
      !     XRAY                EWALD
      ! orth_to_frac    == transpose(recip)
      ! frac_to_orth    == ucell
      ! volume          == volume
      ! unit_cell(1:3)  == dirlng
      ! recip_cell(1:3) == 1.0/(reclng)
      !
      ! NOTE: orth_to_frac and ewald:recip are the same as PDB SCALEn records

      !--------------------------------------------------------------
      ! Read reflection data
      call amopen(allocate_lun(hkl_lun),reflection_infile,'O','F','R')
      read(hkl_lun,*,end=1,err=1) num_hkl
      allocate(hkl_index(3,num_hkl),Fobs(num_hkl),sigFobs(num_hkl), &
            mSS4(num_hkl),test_flag(num_hkl),stat=alloc_status)
      REQUIRE(alloc_status==0)
      do i = 1,num_hkl
         read(hkl_lun,*,end=1,err=1) &
               hkl_index(1:3,i),Fobs(i),sigFobs(i),test_flag(i)
      end do
      close(hkl_lun)
      !call hkl_reduce()

      !--------------------------------------------------------------
      ! call FFT_setup()

      !--------------------------------------------------------------
      if (atom_selection_mask/='') then
         call atommask(natom=natom,nres=nres,prnlev=0, &
               igraph=ih(m04),isymbl=ih(m06),ipres=ix(i02), &
               lbres=ih(m02),crd=x(lcrd), &
               maskstr=atom_selection_mask,mask=atom_selection)
      end if
      call get_mss4(num_hkl, hkl_index, mSS4 )

      return
      1 continue
      write(stdout,'(A)') 'Error reading HKL file.'
      call mexit(stdout,1)
   end subroutine xray_init

   subroutine xray_init_globals()
      pdb_infile = ''
      pdb_outfile = ''
      pdb_read_coordinates = .false.
      pdb_use_segid = .false.
      pdb_wrap_names = .false.
      spacegroup_name = 'P 1'
      reflection_infile = ''
      reflection_infile_format = 'raw'
      reflection_fobs = '4'
      reflection_sigfobs = '5'
      resolution_low = 50
      resolution_high = 0
      xray_weight = 3e+5 ! typical for Residual target
      solvent_mask_probe_radius = 1.0
      solvent_mask_expand = 0.8
      solvent_mask_reflection_outfile = ''
      solvent_mask_outfile = ''
      solvent_mask_update_interval = 0
      solvent_scale = -1
      solvent_bfactor = -1
      fft_method = 1
      fft_grid_size = (/0,0,0/)
      fft_grid_spacing = 0.33_rk_
      fft_bfactor_sharpen = 20
      fft_density_tolerance = 1e-4_rk_
      fft_reflection_tolerance = 1e-4_rk_
      fft_radius_min = 1.0
      fft_radius_max = 4.0
      bfactor_min = 1.0
      bfactor_max = 999.0
      bfactor_refinement_interval = 0
      atom_selection_mask = '!@H='
   end subroutine xray_init_globals

   ! Write X-ray output files and deallocate. Bond info is included
   ! here only to check for places to insert TER in PDB output files.
   subroutine xray_fini()
      implicit none
#     include "extra.h"
      ! local
      integer :: dealloc_status
      if (master .and. pdb_outfile /= '') then
         call xray_write_pdb(trim(pdb_outfile))
      end if

      if (.not.xray_active) return

      deallocate(atom_bfactor,atom_occupancy,atom_scatter_type, &
            atom_selection,residue_chainid,residue_icode, &
            atom_element,atom_altloc,residue_number, &
            scatter_coefficients, &
            hkl_index,Fobs,sigFobs,mSS4,test_flag,density_map, &
            stat=dealloc_status)
      REQUIRE(dealloc_status==0)
   end subroutine xray_fini

   ! NOTE: CNS supports user-defined functions for:
   !
   ! TARGET function (expression to minimize)
   ! DERIVATIVE function (derivative of TARGET)
   ! SELECTION expression to define R-work set
   ! CV-SELECTION expression to define R-free set
   ! MONITOR function (working R-factor)

   subroutine xray_get_derivative(xyz,dxyz)
      use xray_fourier_module
      use xray_utils_module, only: pack_index
      use xray_reciprocal_space_module, only: scale_data
      real(real_kind), intent(in) :: xyz(3,num_atoms)
      real(real_kind), intent(inout) :: dxyz(3,num_atoms)
      ! local
      integer, allocatable :: sel_index(:)
      real(real_kind), allocatable :: frac_xyz(:,:)
      real(real_kind), allocatable :: xray_dxyz(:,:)
      complex(real_kind), allocatable :: Fcalc(:)
      real(real_kind), allocatable, target :: Fcalc_amplitude(:)
      complex(real_kind), allocatable :: dF(:)
      type(hkl_data_scale_type) :: data(2)
      integer :: status, alloc_status, num_selected, dealloc_status

      allocate(sel_index(num_atoms),stat=alloc_status)
      REQUIRE(alloc_status==0)
      call pack_index(atom_selection(:)==1 .and. atom_scatter_type(:)>0, sel_index, num_selected)

      allocate(frac_xyz(3,num_selected),dF(num_hkl),Fcalc(num_hkl), &
            Fcalc_amplitude(num_hkl), &
            xray_dxyz(3,num_selected),stat=alloc_status)
      REQUIRE(alloc_status==0)

      frac_xyz=modulo(matmul(transpose(orth_to_frac),xyz(:,sel_index(1:num_selected))),1.0_rk_)
      
!  dac: try calling fourier_Fcalc() instead of FFT_Fcalc(), since Joe has
!       not yet implemented the latter:
#if 0
      call FFT_Fcalc(num_hkl,Fcalc,test_flag-1, &
            num_selected,frac_xyz, &
            atom_bfactor(sel_index(1:num_selected)), &
            atom_occupancy(sel_index(1:num_selected)), &
            atom_scatter_type(sel_index(1:num_selected)), &
            num_scatter_types,scatter_ncoeffs,scatter_coefficients)
#else
      ! dac note: this is probably in the wrong place
      ! call get_mss4(num_hkl, hkl_index, mSS4 )
      call fourier_Fcalc(num_hkl,hkl_index,Fcalc,mSS4,test_flag, &
            num_selected,frac_xyz, &
            atom_bfactor(sel_index(1:num_selected)), &
            atom_occupancy(sel_index(1:num_selected)), &
            atom_scatter_type(sel_index(1:num_selected)) )
#endif

      Fcalc_amplitude = abs(Fcalc)

      data(1)%f => Fobs
      data(1)%sigma => sigFobs
      data(1)%name = "Fobs"
      data(1)%refine_bfactor = .false.

      data(2)%f => Fcalc_amplitude
      data(2)%name = "Fcalc"
      data(2)%scale = -1
      data(2)%refine_scale = .false.

      call scale_data(num_hkl,hkl_index, & ! selection=ALL
            num_sets=2,data=data, &
            scale_min=1e-4_rk_, scale_max=1e+4_rk_, &
            bfactor_min=1.0_rk_,bfactor_max=100.0_rk_, &
            max_cycles=20,tolerance=1e-4_rk_, &
            reference_set=2,status=status,print=.true.)

      call dTarget_dF(num_hkl, Fobs,Fcalc,selected=test_flag-1,residual=r_free)
      call dTarget_dF(num_hkl, Fobs,Fcalc,selected=test_flag,deriv=dF, &
           residual=r_work, xray_energy=xray_energy)
      xray_energy = xray_weight * xray_energy
      dF = xray_weight * dF

#if 0
! hkl_index and mSS4 are not used
      call FFT_dXYZBQ_dF(num_hkl,dF,test_flag-1, &
            num_selected,frac_xyz, &
            atom_bfactor(sel_index(1:num_selected)), &
            atom_occupancy(sel_index(1:num_selected)), &
            atom_scatter_type(sel_index(1:num_selected)), &
            num_scatter_types,scatter_ncoeffs,scatter_coefficients,xray_dxyz)
#else
      call fourier_dXYZBQ_dF(num_hkl,hkl_index,dF,mSS4,test_flag, &
            num_selected,frac_xyz, &
            atom_bfactor(sel_index(1:num_selected)), &
            atom_scatter_type(sel_index(1:num_selected)), &
            dxyz=xray_dxyz )
#endif
      ! Convert xray_dxyz() back to orthogonal coordinates: 
      xray_dxyz(:,:) = matmul(orth_to_frac,xray_dxyz(:,:))

      dxyz(:,sel_index(1:num_selected)) = dxyz(:,sel_index(1:num_selected)) &
          - xray_dxyz(:,:)

      deallocate(frac_xyz,dF,Fcalc, &
            Fcalc_amplitude, xray_dxyz,stat=dealloc_status)
      REQUIRE(dealloc_status==0)
      deallocate(sel_index,stat=dealloc_status)
      REQUIRE(dealloc_status==0)

   end subroutine xray_get_derivative

   subroutine xray_write_md_state(unit)
      integer, intent(in) :: unit
      write(unit,'(3(1x,A,f14.4))') &
        'E(XRAY)= ',xray_energy,' R(WORK) = ',r_work,' R(FREE)   = ',r_free
   end subroutine xray_write_md_state

   subroutine xray_write_min_state(unit)
      integer, intent(in) :: unit
      write(unit,'(3(1x,A,f13.4))') &
        'E(XRAY) = ',xray_energy,' R(WORK) = ',r_work,' R(FREE)    = ',r_free
   end subroutine xray_write_min_state

end module xray_interface_module
