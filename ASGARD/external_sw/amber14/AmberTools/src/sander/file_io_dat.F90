!<compile=optimized>

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Module file_io_dat contains all of the information regarding file IO objects
module file_io_dat

!+ Specification and control of Amber's Input/Output

! Size of the file names
integer, parameter :: MAX_FN_LEN = 256
integer, parameter :: MAX_LINE_BUF_LEN = 256

! File names
character(len=4096) , save :: groupbuffer    ! Buffer for groupfile lines
character(len=MAX_FN_LEN), save :: mdin      ! Input file
character(len=MAX_FN_LEN), save :: mdout     ! Output file
character(len=MAX_FN_LEN), save :: inpcrd    ! Input coordinate file
character(len=MAX_FN_LEN), save :: parm      ! Topology file
character(len=MAX_FN_LEN), save :: restrt    ! Restart file (output)
character(len=MAX_FN_LEN), save :: refc      ! Reference coordinate file (input)
character(len=MAX_FN_LEN), save :: mdvel     ! MD velocity dump file (output)
character(len=MAX_FN_LEN), save :: mdfrc     ! MD force dump file (output)
character(len=MAX_FN_LEN), save :: mden      ! MD energy dump file (output)
character(len=MAX_FN_LEN), save :: mdcrd     ! MD trajectory file (output)
character(len=MAX_FN_LEN), save :: mdinfo    ! MD info file (output)
character(len=MAX_FN_LEN), save :: mtmd      ! Multiple targetted MD file
character(len=MAX_FN_LEN), save :: vecs      ! Eigenvector file? ??
character(len=MAX_FN_LEN), save :: radii     ! PB Radii file? ??
character(len=MAX_FN_LEN), save :: freqe     ! ??
character(len=MAX_FN_LEN), save :: redir(9)  ! NMR redirection file
character(len=MAX_FN_LEN), save :: rstdip    ! Dipole restart file (output)
character(len=MAX_FN_LEN), save :: mddip     ! MD dipole trajectory (output)
character(len=MAX_FN_LEN), save :: inpdip    ! Input dipole file (input)
character(len=MAX_FN_LEN), save :: groups    ! Group file name (input)
character(len=MAX_FN_LEN), save :: gpes      ! Proc. allocation file (unused)
character(len=MAX_FN_LEN), save :: cpin      ! Constant pH input file
character(len=MAX_FN_LEN), save :: cpout     ! Constant pH output file
character(len=MAX_FN_LEN), save :: cprestrt  ! Constant pH restart file
character(len=MAX_FN_LEN), save :: evbin     ! EVB input file
character(len=MAX_FN_LEN), save :: evbout    ! EVB output file
character(len=MAX_FN_LEN), save :: inptraj   ! Input trajectory file
character(len=MAX_FN_LEN), save :: pimdout   ! PIMD output file
character(len=MAX_FN_LEN), save :: amdlog    ! Log file for accelerated MD
character(len=MAX_FN_LEN), save :: scaledMDlog    ! Log file for scaled MD
character(len=MAX_FN_LEN), save :: cph_dump  ! dump of CpHMD statistics
character(len=MAX_FN_LEN), save :: sechgname ! dump of atomic charges in SEBOMD
#ifdef RISMSANDER
character(len=MAX_FN_LEN), save :: rismcrdfil
character(len=MAX_FN_LEN), save :: rismfrcfil
character(len=MAX_FN_LEN), save :: rismcrdrstfil
character(len=MAX_FN_LEN), save :: rismfrcrstfil
#endif /* RISMSANDER */
#ifdef MPI
character(len=MAX_FN_LEN), save :: remlog        ! REM log file
character(len=MAX_FN_LEN), save :: remtype       ! REM type file
character(len=MAX_FN_LEN), save :: remstripcoord ! Stripped coords (hybrid REMD)
character(len=MAX_FN_LEN), save :: saveenefile   ! Energies of reservoir structs
character(len=MAX_FN_LEN), save :: clusterinfofile ! self-explanatory
character(len=MAX_FN_LEN), save :: reservoirname  ! Name of reservoir file
character(len=MAX_FN_LEN), save :: remd_dimension_file ! Multi-D REMD input file
#endif /* MPI */

character, save :: owrite ! overwrite: N [no owrite] R [replace] or U [append]
character, save :: facc   ! file access: W [write] A [append] R [read]

integer, save :: numgroup  ! Number of groups (-ng CL flag)
integer, save :: nslice    ! EVB slices

! PLUMED variables
integer,                   save :: plumed
character(len=MAX_FN_LEN), save :: plumedfile

! File units
! An I/O Unit resource manager does not exist.
integer, parameter :: MDINFO_UNIT = 7
integer, parameter :: INPCRD_UNIT = 9
integer, parameter :: MDCRD_UNIT = 12
integer, parameter :: MDVEL_UNIT = 13
integer, parameter :: MDFRC_UNIT = 14
integer, parameter :: MDEN_UNIT = 15
integer, parameter :: CNSTPH_UNIT = 18
integer, parameter :: CPOUT_UNIT = 19
integer, parameter :: CPH_DUMP_UNIT = 21
integer, parameter :: INPTRAJ_UNIT = 24
integer, parameter :: EVB_UNIT = 75
integer, parameter :: AMDLOG_UNIT = 77
integer, parameter :: scaledMDLOG_UNIT = 78
integer, parameter :: SCHLEGEL_UNIT = 80 ! EVB
integer, parameter :: PIMD_UNIT = 277
integer, parameter :: SECHGUNIT = 31
#ifdef RISMSANDER
integer, parameter :: RISMCRD_UNIT = 90
integer, parameter :: RISMFRC_UNIT = 91
integer, parameter :: RISMCRDRST_UNIT = 92
integer, parameter :: RISMFRCRST_UNIT = 93
#endif
#ifdef MPI
integer, parameter :: REMIN_UNIT = 32
integer, parameter :: RESERVOIR_UNIT = 39
integer, parameter :: REMSTRIPCOORD_UNIT = 97
integer, parameter :: REMTYPE_UNIT = 98
integer, parameter :: REMLOG_UNIT = 99
integer, parameter :: REMD_DIM_UNIT = 100
#endif /* MPI */

! File related controls and options
character(len=80), save :: title
character(len=80), save :: title1

! Presence of certain namelists in the mdin file
logical, save :: mdin_ewald
logical, save :: mdin_pb
logical, save :: mdin_amoeba
#ifdef APBS
logical, save :: mdin_apbs
logical, save :: sp_apbs
#endif /* APBS */

! File writing input flags
integer, save :: ntpr   ! How often to write energies to mdout
integer, save :: ntwr   ! How often to write restarts
integer, save :: ntwx   ! How often to write snapshots to trajectory
integer, save :: ntwv   ! How often to write velocities to mdvel
integer, save :: ntwf   ! How often to write forces to mdfrc
integer, save :: ntwe   ! How often to write energies to mden
integer, save :: ntpp   ! ??
integer, save :: ioutfm ! Write NetCDF trajectories? (0=no, 1=yes)
integer, save :: ntwprt ! Number of atoms for which to write crds in mdcrd
integer, save :: ntave  ! Interval of steps to collect data for averaging

!      NMRRDR : Contains information about input/output file redirection
!               REDIR and IREDIR contain information regarding
!               LISTIN, LISTOUT, READNMR, NOESY, SHIFTS, DUMPAVE,
!               PCSHIFT and DIPOLE respectively. If IREDIR(I) > 0,
!               then that input/output has been redirected.

integer, save :: iredir(9)

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Provides default file names for the various files
subroutine initialize_fnames

   implicit none

   mdin = 'mdin'
   mdout = 'mdout'
   inpcrd ='inpcrd'
   parm = 'prmtop'
   restrt = 'restrt'
   refc = 'refc'
   mdvel = 'mdvel'
   mden = 'mden'
   mdfrc = 'mdfrc'
   mdcrd = 'mdcrd'
   mdinfo = 'mdinfo'
   mtmd = 'mtmd'
   vecs = 'vecs'
   radii = 'radii'
   freqe = 'dummy'
   redir(:) = ' '
   rstdip = 'rstdip'
   mddip = 'mddip'
   inpdip = 'inpdip'
   groups = ' '
   gpes = ' '
   cpin = 'cpin'
   cpout = 'cpout'
   cprestrt = 'cprestrt'
   evbin = 'evbin'
   evbout = 'evbout'
   sechgname = 'sebomd.chg'
   inptraj = 'inptraj'
   pimdout = 'pimdout'
#ifdef RISMSANDER
   rismcrdfil = ''
   rismfrcfil = ''
   rismcrdrstfil = ''
   rismfrcrstfil = ''
#endif /* RISMSANDER */
#ifdef MPI
   remlog = 'rem.log'
   remtype = 'rem.type'
   remstripcoord = ' '
   saveenefile = 'saveene'
   clusterinfofile = 'cluster.info'
   reservoirname = 'reserv/frame'
   remd_dimension_file = ' '
#endif /* MPI */

end subroutine initialize_fnames

end module file_io_dat
