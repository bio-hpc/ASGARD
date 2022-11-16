! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Definition for EVB modules                                             |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#include "../include/assert.fh"
#include "../include/dprec.fh"

!  +---------------------------------------------------------------------------+
!  |  EVB permitted keywords                                                   |
!  +---------------------------------------------------------------------------+

   module evb_allowed

   implicit none

   save

   integer, parameter :: ndia_type = 2, nxch_type = 6, nevb_dyn = 9

   character(512), dimension(ndia_type), parameter :: &
      dia_type_key = (/ "ab_initio  " &
                      , "force_field" &
                     /)

   character(512), dimension(nxch_type), parameter :: &
      xch_type_key = (/ "constant      " &
                      , "exp           " &
                      , "gauss         " &
                      , "chang_miller  " &
                      , "minichino_voth" &
                      , "dist_gauss    " &
                     /)

   character(512), dimension(nevb_dyn), parameter :: &
      evb_dyn_key = (/ "groundstate  " &
                     , "evb_map      " &
                     , "egap_umb     " &
                     , "bond_umb     " &
                     , "dbonds_umb   " &
                     , "qi_bond_pmf  " &
                     , "qi_dbonds_pmf" &
                     , "qi_bond_dyn  " &
                     , "qi_dbonds_dyn" &
                    /)

   end module evb_allowed

!  +---------------------------------------------------------------------------+
!  |  EVB namelist derived types                                               |
!  +---------------------------------------------------------------------------+

   module evb_nml

   implicit none

   save
!  .............................................................................
!  :  maximum size of EVB namelist array variables                             :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   integer, parameter :: xdg_size = 80
   integer, parameter :: dia_size = 20
   integer, parameter :: xch_size = dia_size * ( dia_size - 1 ) / 2
   integer, parameter :: umb_size = 40
   integer, parameter :: nmorse_size   = 200
   integer, parameter :: nmodvdw_size  = 200
   integer, parameter :: nUFF_size  = 200
   integer, parameter :: nAdihed_size  = 200
   integer, parameter :: nAangle_size  = 200

!  .............................................................................
!  :  EVB diabatic state energy shift                                          :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   type nml_dia_shift
      integer :: st                    ! EVB diabatic state index
      _REAL_  :: nrg_offset            ! energy offset
   endtype

!  .............................................................................
!  :  constant exchange element                                                :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   type nml_xch_cnst
      integer :: ist                   ! EVB diabatic state index
      integer :: jst                   ! EVB diabatic state index
      _REAL_  :: xcnst                 ! coupling between states ist & jst
   endtype
!  .............................................................................
!  :  Warshel's exponential exchange element                                   :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   type nml_xch_exp
      integer :: ist                   ! EVB diabatic state index
      integer :: jst                   ! EVB diatabic state index
      integer :: iatom                 ! atom index involved in bond length
      integer :: jatom                 ! atom index involved in bond length
      _REAL_  :: A                     ! constant prefactor
      _REAL_  :: u                     ! constant prefactor inside exponential
      _REAL_  :: r0                    ! distance between atoms iatom & jatom
   endtype
!  .............................................................................
!  :  gaussian exchange element                                                :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   type nml_xch_gauss
      integer :: ist                   ! EVB diabatic state index
      integer :: jst                   ! EVB diabatic state index
      integer :: iatom                 ! atom index involved in bond length
      integer :: jatom                 ! atom index involved in bond length
      _REAL_  :: A                     ! constant prefactor
      _REAL_  :: sigma                 ! standard deviation
      _REAL_  :: r0                    ! distance between atoms iatom & jatom;
   endtype
!  .............................................................................
!  :  "morsify" specified standard harmonic bonded interactions                :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   type nml_morsify
      integer :: iatom                 ! atom index involved in bond length
      integer :: jatom                 ! atom index involved in bond length
      _REAL_  :: D                     ! well-depth in morse potential
      _REAL_  :: a                     ! alpha in exponential
      _REAL_  :: r0                    ! distance between atoms iatom & jatom;
   endtype
   type( nml_morsify    ) :: morsify(nmorse_size)
!  .............................................................................
!  :  exclude specified van der waals interactions                             :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   type nml_modvdw
      integer :: iatom                 ! atom index of non-bonded interaction
      integer :: jatom                 ! atom index of non-bonded interaction
   endtype
   type( nml_modvdw ) :: modvdw(nmodvdw_size)
!  .............................................................................
!  :  Schlegel DG-EVB: include UFF VDW term in V_ii                            :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   type nml_UFF
      integer :: iatom                 ! UFF VDW pair index
      integer :: jatom                 ! UFF VDW pair index
   endtype
   type( nml_UFF ) :: UFF(nUFF_size)
!  .............................................................................
!  :  Schlegel DG-EVB: include Amber torsion term in V_ii                      :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   type nml_Amber_dihed
      integer :: iatom                 ! dihedral atom index
      integer :: jatom                 ! dihedral atom index
      integer :: katom                 ! dihedral atom index
      integer :: latom                 ! dihedral atom index
      _REAL_  :: V                     ! dihedral force constant
      _REAL_  :: period                ! dihedral periodicity
      _REAL_  :: phase                 ! dihedral phase shift
   endtype
   type( nml_Amber_dihed ) :: Adihed(nUFF_size)

!  .............................................................................
!  :  Schlegel DG-EVB: include Amber angle term in V_ii                        :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   type nml_Amber_angle
      integer :: iatom                 ! angle atom index
      integer :: jatom                 ! angle atom index
      integer :: katom                 ! angle atom index
      _REAL_  :: K                     ! angle force constant
      _REAL_  :: theta                 ! equilibrium angle
   endtype
   type( nml_Amber_angle ) :: Aangle(nUFF_size)

!  .............................................................................
!  :  EVB mapping potential                                                    :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   type nml_emap
      integer :: ist                   ! EVB diabatic state index
      integer :: jst                   ! EVB diabatic state index
      _REAL_  :: lambda                ! progress parameter between ist & jst
   endtype
!  .............................................................................
!  :  EVB harmonic umbrella sampling along an energy gap RC                    :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   type nml_egap_umb
      integer :: ist                   ! EVB diabatic state index
      integer :: jst                   ! EVB diabatic state index
      _REAL_  :: k                     ! EVB energy gap umbrella force constant
      _REAL_  :: ezero                 ! EVB energy gap anchor point
   endtype
!  .............................................................................
!  :  EVB harmonic umbrella sampling along a difference of 2 distances RC      :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   type nml_dbonds_umb
      integer :: iatom                 ! atom index of donor
      integer :: jatom                 ! atom index of transferring particle
      integer :: katom                 ! atom index of acceptor
      _REAL_  :: k                     ! umbrella force constant
      _REAL_  :: ezero                 ! umbrella anchor point
   endtype
!  .............................................................................
!  :  EVB harmonic umbrella sampling along a distance RC                       :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   type nml_bond_umb
      integer :: iatom                 ! atom index involved in distance
      integer :: jatom                 ! atom index involved in distance
      _REAL_  :: k                     ! umbrella force constant
      _REAL_  :: ezero                 ! umbrella anchor point
   endtype
!  .............................................................................
!  :  EVB reactive flux                                                        :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   type nml_react_flux
      integer :: nbias_ndx             !! biasing potential index
      integer :: p_space               !! number of initial momenta sampling
      integer :: minstep               !! minimum no. of steps before checking stopping criteria
      integer :: nsample               !! no. of RC points used for computing RC_rmsd
      _REAL_  :: s0_surface            !! location of dividing surface
      _REAL_  :: s0_toler              !! tolerance for location of dividing surface
      _REAL_  :: RC_rmsd_toler         !! tolerance for determining commitment to a FE well
      _REAL_  :: RC_RS_toler           !! tolerance for determining commitment to the RS well
      _REAL_  :: RC_PS_toler           !! tolerance for determining commitment to the PS well
      character(512) :: RC_type
   endtype

   end module evb_nml

!  +---------------------------------------------------------------------------+
!  |  AMBER external variables that interface with EVB in force                |
!  +---------------------------------------------------------------------------+

   module evb_amber

   use evb_nml, only: nml_morsify, nml_modvdw

   implicit none

   _REAL_ , allocatable :: xnrg(:), xq(:), xf(:,:)

   logical :: evb_amber_ready = .false.  !  already allocated (ready == .true.)

!  .............................................................................
!  :  replace Amber harmonic interaction with Morse                            :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   type( nml_morsify ), allocatable :: morse(:)  ! contains Morse parameters

   _REAL_ , allocatable :: k_harm(:), r0_harm(:) ! Amber harmonic parameters

   logical :: morsify_initialized = .false. 

!  .............................................................................
!  :  modify VDW interactions                                                  :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   type( nml_modvdw  ), allocatable :: vdw_mod(:)

   integer, allocatable :: iac(:), ico(:), cn1(:), cn2(:)

   logical :: modvdw_initialized = .false.

   end module evb_amber

!  +---------------------------------------------------------------------------+
!  |  Global variables used in EVB module                                      |
!  +---------------------------------------------------------------------------+

   module evb_parm

   use evb_nml, only: nml_xch_exp, nml_xch_gauss

   implicit none

   save

!  .............................................................................
!  :  EVB array dimensions and io unit                                         :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   integer :: nevb                     ! No. of EVB diabatic states
   integer :: nxch                     ! No. of coupling terms: nevb*(nevb+1)/2
   integer :: nbias                    ! No. of biased potentials
   integer :: nmorse                   ! No. of "morsifications"
   integer :: nmodvdw                  ! No. of modifications to VDW
   integer :: ntw_evb                  ! Output to evbout every ntw_evb
   integer :: ioe                      ! IO unit for EVB

   logical :: inc_dbonds_RC   = .false.
   logical :: inc_bond_RC     = .false.
   logical :: out_RCdot       = .false.
   logical :: egapRC          = .false.

!  .............................................................................
!  :  EVB diabatic type, exchange type and dynamics type                       :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   character(512) :: dia_type, xch_type, evb_dyn

!  .............................................................................
!  :  EVB states and exchange decryption files                                 :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   integer, allocatable :: evb_dcrypt(:,:)
!dEVB   integer, allocatable :: xch_dcrypt(:)

!  .............................................................................
!  :  biased sampling                                                          :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   integer, allocatable :: bias_ndx(:,:)

   _REAL_ , allocatable :: lambda(:), k_umb(:), r0_umb(:), f(:,:)

   type atm_index
      integer :: iatom                   ! atom index
      integer :: jatom                   ! atom index
      integer :: katom                   ! atom index
   endtype
 
   type( atm_index ), allocatable :: dbonds_RC(:)
   type( atm_index ), allocatable :: bond_RC(:)

!  .............................................................................
!  :  constant energy shift and coupling                                       :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   _REAL_ , allocatable :: C_evb(:), C_xch(:)

!  .............................................................................
!  :  parameters for the exponential and gaussian exchange functional forms    :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   type( nml_xch_exp   ), allocatable :: xch_expdat(:)
   type( nml_xch_gauss ), allocatable :: xch_gaussdat(:)

!  .............................................................................
!  :  Tells if arrays have already been allocated                              :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   logical :: input_ready = .false.     ! already allocated (ready == .true.)

   end module evb_parm

!  +---------------------------------------------------------------------------+
!  |  Exchange element functional form types                                   |
!  +---------------------------------------------------------------------------+

   module evb_xchff

   implicit none

   save

!  .............................................................................
!  :  Warshel's exponential functional form                                    :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   type xch_warshel_type 
      _REAL_ :: xch                    ! coupling element
      _REAL_ , pointer :: gxch(:)      ! gradient of xch 
   endtype

!  .............................................................................
!  :  gaussian functional form                                                 :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   type xch_gauss_type
      _REAL_ :: xch                    ! coupling element
      _REAL_ , pointer :: gxch(:)      ! gradient of xch 
   endtype

   type( xch_warshel_type ), allocatable :: xch_warshel(:)
   type( xch_gauss_type   ), allocatable :: xch_gauss(:)

   logical :: xchff_warshel_ready = .false.   ! (ready == .true.)
   logical :: xchff_gauss_ready   = .false.   ! (ready == .true.)

   end module evb_xchff

!  +---------------------------------------------------------------------------+
!  |  Schlegel-Sonnenberg distributed gaussian EVB                             |
!  +---------------------------------------------------------------------------+

   module schlegel

   implicit none
   save

!  .............................................................................
!  :  conversion factors                                                       :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   _REAL_  , parameter :: Bohr2Angstrom   =   0.529177208319
   _REAL_  , parameter :: Hartree2Kcalmol = 627.50947169

!  .............................................................................
!  :  dimensions                                                               :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   integer :: ndg                      ! number of distributed Gaussians
   integer :: ncoord                   ! number of coordinates 
   integer :: natm                     ! number of atoms
   integer :: dgdim                    ! 1 + ncoord + ncoord^2
   logical :: dg_ready = .false.       ! already allocated (ready == .true.)
   _REAL_  :: diis_tol = 1.0d-9        ! Convergence criterium for DIIS
   _REAL_  :: svd_cond = 1.0d-9        ! Singular value decomposition tolerance

!dEVB   _REAL_ , parameter :: gmat_tol = 1.0d-9 !! define zero as less than this number

   character(512) :: full_solve = "berny_solve"

!  .............................................................................
!  :  Redundant internal coordinate definitions                                :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   integer :: nbond                    ! No. of bonds in internals
   integer :: nangle                   ! No. of angles in internals
   integer :: ndihed                   ! No. of dihedrals in internals
   integer :: ncart                    ! No. of Cartsians
   integer :: nselect                  ! No. of selected coordinates for DG-EVB
   integer :: rdgdim                   ! No. of data pts. in each DG
   integer, allocatable :: ibond(:,:)  ! atom indices for bonds
   integer, allocatable :: iangle(:,:) ! atom indices for angles
   integer, allocatable :: idihed(:,:) ! atom indices for dihedrals
   integer, allocatable :: scoord(:)   ! decryption for selected coordinates

   logical :: use_cartesian = .false. 
   logical :: init_ready = .false.     !! already allocated (ready == .true.)

!  .............................................................................
!  :  type for external DG data                                                :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   type berny_type
      _REAL_ , pointer :: q(  :)       ! coordinate
      _REAL_ , pointer :: d(  :)       ! gradient
      _REAL_ , pointer :: k(:,:)       ! hessian
      _REAL_           :: v            ! energy
      _REAL_ , pointer :: qcart(:)     ! Cartesian coordinate
      character(512)   :: filename     ! external data filename
   endtype
!  .............................................................................
!  :  External ab initio data for DG method                                    :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   type( berny_type ), allocatable :: xdat_xdg(:)   ! all DG points
   type( berny_type ), allocatable :: xdat_min(:)   ! DG at minimum geometry
   type( berny_type ), allocatable :: xdat_ts (:)   ! DG at TS geometry

!  .............................................................................
!  :  type for distributed gaussian                                            :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   type dg_type
      _REAL_         :: stol           ! coordinate selection tolerance
      character(512) :: stype          ! coordinate selection type
      character(512) :: lin_solve      ! linear solutions method
      character(512) :: xfile_type     ! DG external file type
   endtype

   type( dg_type ) :: dist_gauss

!  .............................................................................
!  :  Input parameters                                                         :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   _REAL_ , allocatable :: alpha(:)         ! Gaussian function parameter
   _REAL_ , allocatable :: dg_weight(:)     ! DG weights
   _REAL_ :: gbasis_weight(3)               ! weight on each basis type

!  .............................................................................
!  :  Intermediate arrays involving diabatic state                             :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   _REAL_ , pointer ::   vmm(    :)         ! diabatic state energy 
   _REAL_ , pointer ::  dvmm(  :,:)         ! 1st derivative of energy
   _REAL_ , pointer :: ddvmm(:,:,:)         ! 2nd derivative of energy

!  .............................................................................
!  :  Intermediate arrays involving coupling term                              :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   _REAL_ :: vkl                            ! coupling term 
   _REAL_ , pointer ::  dvkl(  :)           ! 1st derivative of coupling
   _REAL_ , pointer :: ddvkl(:,:)           ! 2nd derivative of coupling

!  .............................................................................
!  :  UFF VDW interaction in V_ii                                              :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   integer :: nUFF                          ! No. of UFF VDW pairs
   integer, allocatable :: nUFF_dcrypt(:,:) ! atom indices in VDW pairs
   integer, allocatable :: atomic_numbers(:)! atomic numbers

   _REAL_ :: UFF_scale                      ! scale factor for UFF interaction
   _REAL_ , pointer :: xvdw(:)              ! VDW bond length parameters
   _REAL_ , pointer :: dvdw(:)              ! VDW well depth parameters
   _REAL_ , pointer :: zvdw(:)              ! VDW shape parameters
   _REAL_ , pointer :: Avdw(:)              ! E_vdw = A e(-Bx) - C_6/x^6
   _REAL_ , pointer :: Bvdw(:)              ! E_vdw = A e(-Bx) - C_6/x^6

!  .............................................................................
!  :  Amber torsion interaction in V_ii                                        :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   integer :: nAdihed                         ! No. of Amber dihedrals
   integer, allocatable :: nAdihed_dcrypt(:,:)! atom indices in dihedral

   _REAL_ , allocatable :: V(:)               ! dihedral force constant
   _REAL_ , allocatable :: period(:)          ! periodicity of dihedral
   _REAL_ , allocatable :: phase(:)           ! phase shift of dihedral

!  .............................................................................
!  :  Amber angle interaction in V_ii                                          :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   integer :: nAangle                         ! No. of Amber angles 
   integer, allocatable :: nAangle_dcrypt(:,:)! atom indices in angle

   _REAL_ , allocatable :: K_angle(:)         ! angle force constant
   _REAL_ , allocatable :: theta(:)           ! equilibrium angle

!  .............................................................................
!  :  solution and Cartesian <-> internal transformation                       :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   _REAL_ , pointer :: b(:)                 ! solution of linear equations
   _REAL_ , pointer :: wdcmat(:,:)          ! Wilson-Decius-Cross trans. matrix


   end module schlegel

!  +---------------------------------------------------------------------------+
!  |  Internal arrays for the Miller quantum instanton rate                    |
!  +---------------------------------------------------------------------------+

   module miller

   implicit none
   save

   integer, parameter :: ndiv = 2      ! No. of dividing surfaces
   integer :: i_qi = 0                 ! ( i_qi > 0 ) => perform QI
   integer :: div_ndx(ndiv)            ! bead index for dividing surface

   _REAL_ , pointer :: gradRC(:,:)     ! gradient of RC
   _REAL_ , pointer :: gradRC_norm(:)  ! mass-weighted norm of gradient of RC

   _REAL_ :: qi_fact1                  ! m ( i p / (2 hbar B) )^2
   _REAL_ :: qi_fact2                  ! - m p / (hbar^2 B^2)
   _REAL_ :: qi_fact3                  ! 2 d p / B^2
   _REAL_ :: qi_fact4                  ! 4 m p / (hbar^2 B^3)


!! :::: these variables are mirrored in the nml_reactive_flux ::::::::::::::
   integer :: nbias_ndx                ! biasing potential index
   integer :: p_space                  ! number of initial momenta sampling
   integer :: minstep                  ! minimum No. of steps before checking stopping criteria
   integer :: nsample                  ! No. of RC points used for computing RC_armsd
   _REAL_  :: s0_surface               ! location of dividing surface
   _REAL_  :: s0_toler                 ! tolerance for location of dividing surface
   _REAL_  :: RC_rmsd_toler            ! tolerance for location of dividing surface
   _REAL_  :: RC_RS_toler              ! tolerance for determining commitment to the RS well
   _REAL_  :: RC_PS_toler              ! tolerance for determining commitment to the PS well
   character(512) :: RC_type

!! :::: these variables are mirrored in the nml_reactive_flux :::::::::::::::


!  type corr_type
!     _REAL_           :: RC_s0        !! actual location of dividing surface
!     _REAL_ , pointer :: gradRC(:)    !! gradient of RC
!     _REAL_ , pointer :: coord_DS(:)  !! coordinates at dividing surface
!  endtype

!  type( corr_type ), allocatable :: corr(:)   !! data required for flux-flux and delta-delta
                                       ! correlation functions

   integer :: pm1                      ! p - 1 slice
   integer :: php1                     ! p/2 + 1 slice
   integer :: phm1                     ! p/2 - 1 slice

   _REAL_ :: dof                       ! degrees of freedom

!  .............................................................................
!  :  Thermodynamic integration by mass                                        :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   type ti_mass_type
      _REAL_  :: lambda                ! (1-lambda) mass_i + lambda mass_f
      _REAL_  :: init                  ! initial isotope mass
      _REAL_  :: final                 ! final isotop mass
      integer :: atm_ndx               ! atomic index for the TI
   endtype

   type( ti_mass_type ) :: ti_mass(2)  ! Allow for secondary isotope effect
   integer :: nti_mass                 ! No. of mass modifications
   logical :: do_ti_mass = .false.     ! perform TI by mass

   _REAL_  :: dlnQ_dl                  ! (d/dl) ln Q
   _REAL_  :: dlnCdd_dl                ! (d/dl) ln C_dd

   end module miller

!  +---------------------------------------------------------------------------+
!  |  EVB data                                                                 |
!  +---------------------------------------------------------------------------+

   module evb_data

   implicit none

!  .............................................................................
!  :  EVB data averages                                                        :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   _REAL_ :: evb_nrg(3)                 ! {V_tot, V_EVB, V_umb}
   _REAL_ :: evb_nrg_ave(3) = 0.0d0     ! Averages for above energies
   _REAL_ :: evb_nrg_rms(3) = 0.0d0     ! RMS for above averages
   _REAL_ :: evb_nrg_tmp(3) = 0.0d0     ! temporary array
   _REAL_ :: evb_nrg_old(3) = 0.0d0     ! temporary array
   _REAL_ :: evb_nrg_tmp2(3) = 0.0d0    ! temporary array
   _REAL_ :: evb_nrg_old2(3) = 0.0d0    ! temporary array
!  .............................................................................
!  :  EVB matrix arrays                                                        :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   type evb_mat_type
      _REAL_ , pointer :: evb_mat  (:,:)    ! EVB matrix
      _REAL_ , pointer :: evb_vec0 (  :)    ! EVB ground-state vector
      _REAL_           :: evb_nrg           ! EVB ground-state energy
      _REAL_ , pointer :: grad_xch (  :)    ! gradient of exchange term
      _REAL_ , pointer :: b1_kl    (  :)    ! gradient of V_ij^2
      logical :: ready = .false.            ! (ready == .true.)
   endtype 

   type( evb_mat_type ),save :: evb_Hmat

   logical :: evb_frc_ready = .false.       ! (ready == .true.)

!  .............................................................................
!  :  EVB force arrays                                                         :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   type evb_frc_type
      _REAL_ , pointer :: evb_f(  :)        ! -gradient of EVB potential
      _REAL_ :: evb_nrg                     ! EVB energy
      _REAL_ :: evb_nrg_ave(2)  = 0.0d0     ! Average EVB energy
      _REAL_ :: evb_nrg_rms(2)  = 0.0d0     ! RMS of EVB energy
      _REAL_ :: evb_nrg_tmp(2)  = 0.0d0
      _REAL_ :: evb_nrg_tmp2(2) = 0.0d0
      _REAL_ :: evb_nrg_old(2)  = 0.0d0
      _REAL_ :: evb_nrg_old2(2) = 0.0d0
      logical :: evb_ave = .false. 
      logical :: evb_rms = .false.
   endtype

   type( evb_frc_type  ),save :: evb_frc, evb_vel0, fevb_tmp

!  .............................................................................
!  :  EVB biased sampling arrays                                               :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   type evb_bias_type
      _REAL_ , pointer :: RC       (  :)    ! RC in umbrella sampling
      _REAL_ , pointer :: nrg_bias (  :)    ! energy of H = H_0 + H_biased
      _REAL_ , pointer :: evb_fbias(:,:)    ! -gradient of H = H_0 + H_biased
      _REAL_ :: nrg_bias_ave  = 0.0d0
      _REAL_ :: nrg_bias_rms  = 0.0d0
      _REAL_ :: nrg_bias_tmp  = 0.0d0
      _REAL_ :: nrg_bias_tmp2 = 0.0d0
      _REAL_ :: nrg_bias_old  = 0.0d0
      _REAL_ :: nrg_bias_old2 = 0.0d0
   endtype

   type( evb_bias_type  ),save :: evb_bias

!  .............................................................................
!  :  EVB output                                                               :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   _REAL_ :: nrg_frc(3), nrg_vel0(3)

!  .............................................................................
!  :  TST dynamical frequency factor                                           :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   _REAL_ :: RCdot

!  .............................................................................
!  :  QI dynamical factors                                                     :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   _REAL_ :: f_v, F_QI, G_QI


   end module evb_data

!  +---------------------------------------------------------------------------+
!  |  Reactive flux                                                            |
!  +---------------------------------------------------------------------------+

   module wigner

   implicit none

!  .............................................................................
!  :  Reactive flux dynamics for determining the kappa factor                  :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

!! :::: these variables are mirrored in the nml_reactive_flux ::::::::::::
   integer :: nbias_ndx                ! biasing potential index
   integer :: p_space                  ! number of initial momenta sampling
   integer :: minstep                  ! minimum No. of steps before stopping 
   integer :: nsample                  ! No. of pts. used for computing RC_rmsd
   _REAL_  :: s0_surface               ! location of dividing surface
   _REAL_  :: s0_toler                 ! tolerance for locating s0 
   _REAL_  :: RC_rmsd_toler            ! rmsd for RC 
   _REAL_  :: RC_RS_toler              ! tolerance for commitment to the RS well
   _REAL_  :: RC_PS_toler              ! tolerance for commitment to the PS well
   character(512) :: RC_type
!! :::: these variables are mirrored in the nml_reactive_flux ::::::::::::

   integer :: nreact_RC = 30000        ! array limit for RC_sampled
   integer :: RC_ndx = 0               ! index for sampled RC
   integer :: p_ndx  = 0               ! index for sampled momentum distribution

   logical :: init_rflux = .true.      ! initialize reactive flux
   logical :: back_prop  = .true.      ! perform backwards propagation
   logical :: new_prop   = .true.      ! perform reactive flux from new config

   type rflux_type
      _REAL_           :: RC_s0             ! location of dividing surface
      _REAL_ , pointer :: normal_vect(:)    ! unit vector normal to s0 
      _REAL_ , pointer :: normal_velv(:)    ! velocity vector normal to the s0 
      _REAL_ , pointer :: RC_sampled (:)    ! RC during backward & forward paths
      _REAL_ , pointer :: coord_TS   (:)    ! initial coordinates at s0 
      _REAL_ , pointer :: force_TS   (:)    ! initial forces at s0
      _REAL_ , pointer :: xi         (:)    ! recrossing factor
      _REAL_ , pointer :: RC_trj     (:)    ! Sampled RC (temporary storage)
   endtype

   type( rflux_type ), save :: rflux

   end module wigner

!  +---------------------------------------------------------------------------+
!  |  For debugging EVB                                                        |
!  +---------------------------------------------------------------------------+

   module evb_check

   implicit none

   integer, parameter :: what_size = 50

   _REAL_  :: fdeps = 1.0d-9           ! finite difference step size
   _REAL_  :: debug_toler = 1.0d-42    ! tolerance for numer. vs anal. check

   logical, save :: full_evb_debug = .false.
   logical, save :: schlegel_debug = .false.
   logical, save :: gradRC_debug   = .false.
   logical, save :: morsify_debug  = .false.
   logical, save :: mod_vdw_debug  = .false.
   logical, save :: xwarshel_debug = .false.
   logical, save :: xgauss_debug   = .false.
   logical, save :: dbonds_debug   = .false.
   logical, save :: bond_debug     = .false.
   logical, save :: evb_pimd_debug = .false.

   type evb_debug_type
      _REAL_  :: schlegel_toler = 1.0d-5    ! tolerance for debugging DG
      character(512) :: what(what_size)     ! Which EVB part to debug
      character(512) :: fmathnb             ! Mathematica notebook output file
   endtype

   type( evb_debug_type ), save :: evb_debug

   character(512), save :: ab_initio_gridfile     ! ab initio relaxed scan grid 
   
   contains

!  .............................................................................
!  :  Compute rms deviation                                                    :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   function deviation ( a, b )

   implicit none

   intrinsic :: size
   _REAL_ , intent(in), dimension(:) :: a, b
   _REAL_ :: deviation(2), largest
   intrinsic :: sqrt, dot_product, maxval, abs

   !  ..........................................................................

   _REAL_ , dimension( size(a) ) :: c

   if( size( a ) /= size( b ) ) then
      write(6,'(A)') 'ERROR: Cannot compare arrays of two different sizes'
      write(6,'(A,I8,A)') 'a(', size(a), ')'
      write(6,'(A,I8,A)') 'b(', size(b), ')'
      stop 6
   endif

   c(:) = a(:) - (b)
   largest = maxval( abs( c(:) ) )

   if( size(a) < 2 ) then
      deviation(1) = abs( c(1) )
   else
      deviation(1) = sqrt( dot_product( c, c ) )
   endif

   deviation(2) = largest

   end function deviation

   end module evb_check

