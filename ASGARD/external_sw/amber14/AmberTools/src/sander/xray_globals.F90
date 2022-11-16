! <compile=optimized>
#include "../include/dprec.fh"
module xray_globals_module
   use file_io_dat, only : MAX_FN_LEN
   implicit none
   public

   _REAL_, private :: r_dummy
   integer, parameter :: real_kind = kind(r_dummy)
   integer, parameter :: rk_ = real_kind

   ! HP Fortran uses unit 7 for stderr. Most others use 0.
   ! F2003 gives the actual values in intrinsic module ISO_FORTRAN_ENV.
   ! SANDER does not use stderr, so it is assigned to stdout here,
   ! which is really mdout.
   integer, parameter :: STDERR=6, STDIN=5, STDOUT=6

   logical, save :: xray_active = .false.

   !-------------------------------------------------------------------
   ! NameList Input Parameters
   integer, parameter :: REFL_LABEL_MAXLEN=32

   character(len=MAX_FN_LEN), save :: pdb_infile, pdb_outfile

   ! If true, PDB coordinates will overwrite INPCRD coordinates.
   ! NOTE: the cell still comes from the INPCRD!
   logical, save :: pdb_read_coordinates

   ! If TRUE, write the full 4-character ChainID to the SegID field.
   logical, save :: pdb_use_segid, pdb_wrap_names

   ! Standard condensed spacegroup name, or integer spacegroup number
   character(len=16), save :: spacegroup_name

   ! Filename for reflection input file.
   character(len=MAX_FN_LEN), save :: reflection_infile

   ! raw: unlabelled columns, starting with H,K,L.
   ! Input labels are the column number.
   character(len=16), save :: reflection_infile_format

   ! Labels for input Fobs, sigFobs. More labels in the future.
   character(len=REFL_LABEL_MAXLEN), save :: &
         reflection_Fobs, reflection_sigFobs

   ! Resolution limits for all X-ray calculations
   real(real_kind), save :: resolution_low, resolution_high

   ! Weight term for X-ray force:
   real(real_kind), save :: xray_weight

   ! Solvent mask generation parameters
   real(real_kind), save :: solvent_mask_probe_radius, solvent_mask_expand

   ! Output file for bulk-solvent reflections (Fbulk) and mask
   character(len=MAX_FN_LEN), save :: solvent_mask_reflection_outfile
   character(len=MAX_FN_LEN), save :: solvent_mask_outfile

   ! Number of cycles between updates of the bulk-solvent mask
   integer, save :: solvent_mask_update_interval

   _REAL_, save :: solvent_scale, solvent_bfactor

   ! 0=direct, fft=1. This is currently is a boolean option, but
   ! may have other fft options for spacegroup-optimized fft.
   integer, save :: fft_method

   ! Define specific grid size for FFT
   integer, save :: fft_grid_size(3)

   ! Maximum grid spacing for automatic grid-size determination
   real(real_kind), save :: fft_grid_spacing

   ! Fourier sharpening factor (value subtracted from atomic B-factor)
   real(real_kind), save :: fft_bfactor_sharpen

   ! Maximum density outside of the cutoff-radius for atomic density gridding
   real(real_kind), save :: fft_density_tolerance

   ! Maximum reflection intensity at resolution_max; enforces a minimum radius
   real(real_kind), save :: fft_reflection_tolerance

   ! Maximum and minimum radius for atomic gridding
   real(real_kind), save :: fft_radius_min, fft_radius_max

   ! Maximum and minimum B-factor
   real(real_kind), save :: bfactor_min, bfactor_max

   integer, save :: bfactor_refinement_interval

   character(len=132) :: atom_selection_mask

   !----------------------------------------------------------------------------
   ! GLOBALS:

   !----------------------------------------------------------------------------
   ! Atom data:
   integer, save :: num_atoms, num_residues
   real(real_kind), allocatable, save :: atom_bfactor(:), atom_occupancy(:)
   integer, allocatable, save :: atom_scatter_type(:)
   integer, allocatable, save :: atom_selection(:)

   ! Residue and atom data that may become SANDER globals:
   character(len=4), allocatable, save :: residue_chainid(:), residue_icode(:)
   character(len=4), allocatable, save :: atom_element(:), atom_altloc(:)
   integer, allocatable, save :: residue_number(:)

   !----------------------------------------------------------------------------
   ! Reflection data:

   integer, save :: num_hkl
   integer, allocatable, save :: hkl_index(:,:) ! (3,num_hkl)

   real(real_kind), allocatable, target, save :: Fobs(:), sigFobs(:)
   real(real_kind), allocatable, save :: mSS4(:)
   integer, allocatable, save :: test_flag(:)
   integer, allocatable, save :: test_selection(:), work_selection(:)

   !----------------------------------------------------------------------------
   ! Symmetry, unit cell, and transformations:

   real(real_kind), save :: cell_volume
   real(real_kind), save :: unit_cell(6), recip_cell(6)
   real(real_kind), dimension(3,3), save :: orth_to_frac, frac_to_orth
   integer, parameter :: MAX_SYMMOPS = 16 ! actually 96
   real(real_kind), save :: symmop(3,4,MAX_SYMMOPS), symmop_inv(3,4,MAX_SYMMOPS)
   integer, save :: num_symmops

   integer, save :: spacegroup_number
   integer, save :: au_type ! Laue code index
   integer, save :: system  ! i.e. SYMM_TRICLINIC

   character(len=1), save :: lattice_type

   real(real_kind), save :: fft_normalization
   real(real_kind), save :: xray_energy, r_work, r_free

   !----------------------------------------------------------------------------
   ! Electron-density and FFT configuration and conversions:

   integer, save :: grid_size(3), recip_size(3)
   real(real_kind), dimension(3,3), save :: orth_to_grid, grid_to_orth
   real(real_kind), dimension(3), save :: frac_to_grid, angstrom_to_grid

   !dimension(0:grid_size(1)-1,0:grid_size(2)-1,0:grid_size(3)-1)
   real(real_kind), allocatable, target, save :: density_map(:,:,:)
   ! NOTE: the reciprocal map is complex, but an in-place complex-to-real
   ! FFT is used. The complex map is an alias of the real-space map,
   ! but F90 won't let us equivalence allocatable or pointer arrays.
   ! To access the array as complex numbers, it must be passed via a
   ! non-prototyped external procedure.

   !----------------------------------------------------------------------------
   ! Fourier coefficients:

   ! number of gaussian coefficients, including the constant term
   integer, parameter :: scatter_ncoeffs = 5
   integer, save :: num_scatter_types
   real(real_kind), allocatable, save :: scatter_coefficients(:,:,:)
                                                       !(2,ncoeffs,ntypes)
#ifdef QMMM_OMP
   integer, save :: fft3_max_threads = 2
#endif

end module xray_globals_module
