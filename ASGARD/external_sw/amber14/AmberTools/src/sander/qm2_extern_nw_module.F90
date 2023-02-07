#include "../include/dprec.fh"
module qm2_extern_nw_module
! ----------------------------------------------------------------
! Interface for NWChem based QM MD 
!
! Currently supports:
! pure QM
!
! Initial version:
! Mark J. Williamson (Unilever Centre, Cambridge)
! Date: February 2012
!
! Extensions:
! Matthew A. Clark (SDSC)
! Andreas W. Goetz (SDSC)
! Date: June 2012
! 
! Still TODO:
!   Templates
!   Dipole support
!   PIMD
!   Check energy conservation and possibly let users adjust
!    integral neglect thresholds / XC quadrature grid params etc
! 
! ----------------------------------------------------------------

  use qm2_extern_util_module, only: debug_enter_function, debug_exit_function

  implicit none

  private
  public :: get_nw_forces

  character(len=*), parameter :: module_name = "qm2_extern_nw_module"

  type nw_nml_type
     character(len=20) :: method
     character(len=20) :: basis
     character(len=20) :: guess
     integer :: scf_conv
     integer :: ntpr
     integer :: num_threads
     integer :: debug
     
     ! Deprecated
     integer :: charge
     integer :: spinmult
  end type nw_nml_type

contains

  ! --------------------------------------------
  ! Get QM energy and forces from NWChem
  ! --------------------------------------------
  subroutine get_nw_forces( do_grad, ntpr_default, id, &
       nqmatoms, qmcoords, qmtypes, nclatoms, clcoords, escf, dxyzqm, dxyzcl, &
       charge, spinmult )

    use qm2_extern_util_module, only: print_results, check_installation
    use constants, only: CODATA08_AU_TO_KCAL, CODATA08_A_TO_BOHRS, ZERO
    use file_io_dat

    implicit none

    logical, intent(in) :: do_grad              ! Return gradient/not
    integer, intent(in) :: ntpr_default         ! frequency of printing
    character(len=3), intent(in) :: id          ! ID number for PIMD or REMD
    integer, intent(in) :: nqmatoms             ! Number of QM atoms
    _REAL_,  intent(in) :: qmcoords(3,nqmatoms) ! QM atom coordinates
    integer, intent(in) :: qmtypes(nqmatoms)    ! QM atom types (nuclear charge in au)
    integer, intent(in) :: nclatoms             ! Number of MM atoms
    _REAL_,  intent(in) :: clcoords(4,nclatoms) ! MM atom coordinates and charges in au
    _REAL_, intent(out) :: escf                 ! SCF energy
    _REAL_, intent(out) :: dxyzqm(3,nqmatoms)   ! SCF QM force
    _REAL_, intent(out) :: dxyzcl(3,nclatoms)   ! SCF MM force
    integer, intent(in) :: charge, spinmult     ! Charge and spin multiplicity

    type(nw_nml_type), save      :: nw_nml
    logical, save                :: first_call = .true.
    integer                      :: i
    character (len=512)          :: call_buffer
    character(len=*), parameter  :: program='nwchem'
    character(len=*), parameter  :: basename = 'nwchem'
    character(len=*), parameter  :: inpext   = '.nw'
    character(len=*), parameter  :: logext   = '.log'
    character(len=*), parameter  :: bqfext   = '.bqforce.dat'
    character(len=*), parameter  :: movecext = '.movecs' ! MO vector file
    character(len=*), parameter  :: dbext    = '.db'     ! Database file; essentially current state of calculation
    character(len=25)            :: inpfile, movecfile, logfile, bqforcefile, dbfile
    ! Need to prepend subdirectory if doing REMD, PIMD or multi-region QM/MM. 
    !   This is triggered if 'id' is defined (not empty). 
    character(len=25)            :: subdir 

    ! for system calls
    integer :: system
    integer :: stat

    ! assemble input - / output data filenames
    inpfile     = basename//trim(id)//inpext
    logfile     = basename//trim(id)//logext
    movecfile   = basename//trim(id)//movecext
    dbfile      = basename//trim(id)//dbext
    bqforcefile = basename//trim(id)//bqfext

    ! Setup on first call
    if ( first_call ) then
      first_call = .false.
      write (6,'(/,a,/)') '  >>> Running calculations with NWChem <<<'
      call get_namelist( ntpr_default, nw_nml )
      call print_namelist( nw_nml ) 
      call check_installation( program, id, .true., nw_nml%debug  )

      write (6,'(80a)') ('-', i=1,80)
      write (6,'(a)') '   4.  RESULTS'
      write (6,'(80a)') ('-', i=1,80)
      ! Remove old inpfile, logfile, database file, and movecs file at the 
      ! beginning of a run so only the latest run is stored.
      stat = system('rm -f '//inpfile//' '//dbfile//' '//logfile//' '//movecfile)

      if ( stat /= 0 ) then
        call sander_bomb('get_nw_forces (qm2_extern_nw_module)', & 
          'Error with system call (removing files)', &
          'Will quit now.')
      end if
    end if
    
    call write_inpfile( trim(inpfile), &
         nqmatoms, qmcoords, qmtypes, nclatoms, clcoords, &
         nw_nml, do_grad, charge, spinmult )

    if ( nw_nml%debug > 1 ) then
      write(6,'(a)') ' NWChem input file successfully written; calling NWChem...'
    end if

    ! Run NWChem
    ! Separate runs into different directories if we are doing PIMD or REMD
    subdir = ''
    call_buffer = ''
    if (trim(id)/='') then 
      subdir = './'//trim(id)//'/'
      call_buffer = ' mkdir -p '//trim(subdir)//'; cd '//trim(subdir)//&
        '; mv ../'//inpfile//' .;'
    end if
    ! remove restart files if we don't want to read an MO guess
    if ( trim(nw_nml%guess) /= 'read' ) then
       call_buffer = trim(call_buffer)//'rm -f '//movecfile//' '//dbfile//';'
    end if
    ! call nwchem
    call_buffer = trim(call_buffer)//' '//program//' '//inpfile//'>'//logfile
    stat = system(trim(call_buffer))
    if ( stat /= 0 ) then
      call sander_bomb('get_nw_forces (qm2_extern_nw_module)', & 
        'Error with system call (executing NWChem)', &
        'Will quit now.')
    end if

    if ( nw_nml%debug > 0 ) then    
      write(6,'(a)') ' NWChem execution success; Processing NWChem results...'
    end if    

    ! Call read_results - retrieve data from NWChem .log file
    ! Will output data to escf and dxyqm for pure QM or MM runs
    ! For QM/MM runs will also return dxyzcl containing electric field strength
    
    ! Search in subdir for logfile if doing PIMD
    ! Otherwise, search current directory
    call read_results( nw_nml, trim(subdir)//trim(logfile), trim(subdir)//trim(bqforcefile), &
      nqmatoms, escf, dxyzqm, nclatoms, dxyzcl, do_grad, nw_nml%debug )


    ! Save copy of last input and log files 
    stat = system('mv '//trim(subdir)//inpfile//' '//trim(subdir)//'old.'//inpfile)
    stat = stat + system('mv '//trim(subdir)//logfile//' '//trim(subdir)//'old.'//logfile)
    if ( do_grad .and. ( nclatoms > 0 ) ) then
       stat = stat + system('mv '//trim(subdir)//bqforcefile//' '//trim(subdir)//'old.'//bqforcefile)
    end if
    if ( stat /= 0 ) then
      call sander_bomb('get_nw_forces (qm2_extern_nw_module)', & 
        'Error with system call (moving / removing files)', &
        'Will quit now.')
    end if

    ! Convert gradient from au to kcal/(mol*A)
    if ( do_grad ) then
      dxyzqm(:,:) = dxyzqm(:,:) * CODATA08_AU_TO_KCAL * CODATA08_A_TO_BOHRS
      if ( nclatoms > 0 ) then
        dxyzcl(:,:) = dxyzcl(:,:) * CODATA08_AU_TO_KCAL * CODATA08_A_TO_BOHRS
      end if
    else
      dxyzqm = ZERO
      if ( nclatoms > 0 ) dxyzcl = ZERO
    end if
    
    escf = escf * CODATA08_AU_TO_KCAL

    call print_results( 'qm2_extern_nw_module', escf, nqmatoms, dxyzqm,&
      nw_nml%debug, nclatoms, dxyzcl )

  end subroutine get_nw_forces

  ! ---------------------------------------------
  ! Read NWChem namelist values from file mdin,
  ! use default values if none are present.
  ! ---------------------------------------------
    
 
  subroutine get_namelist(ntpr_default, nw_nml)

    implicit none
    integer, intent(in) :: ntpr_default
    type(nw_nml_type), intent(out) :: nw_nml

    character(len=20) :: method, basis, guess
    integer :: debug
    integer :: scf_conv, ntpr, num_threads
    integer :: charge, spinmult ! deprecated
    namelist /nw/ method, basis, guess, scf_conv, ntpr, &
      num_threads, debug,&
      charge, spinmult

    integer :: ierr

    ! Set default values for nw namelist values
    method       = 'BLYP'
    basis        = '6-31G*'
    guess        = 'read'
    scf_conv     = 8 
    ntpr         = ntpr_default
    num_threads  = 1
    debug        = 0

    ! These are now deprecated and should be specified in the &qmmmm namelist
    charge   = -351
    spinmult = -351

    ! Read namelist
    rewind 5
    read(5,nml=nw,iostat=ierr)

    if ( ierr > 0 ) then
      call sander_bomb('get_namelist (qm2_extern_nw_module)', &
        '&nw namelist read error', &
        'Please check your input.')
    else if ( ierr < 0 ) then
      write(6,'(a,/,a)') '&nw namelist read encountered end of file', &
        'Please check your input if the calculation encounters a problem'
    end if

    if ( charge /= -351 .or. spinmult /= -351 ) then
      call sander_bomb('get_namelist (qm2_extern_nw_module)', &
        'The charge and spin keywords are deprecated', &
        'Please specify charge (qmcharge) and spin multiplicity (spin) in the &qmmm namelist.')
    end if

    ! Assign namelist values to nw_nml data type
    nw_nml%method       = method
    nw_nml%basis        = basis
    nw_nml%guess        = guess
    nw_nml%scf_conv     = scf_conv
    nw_nml%ntpr         = ntpr
    nw_nml%num_threads  = num_threads
    nw_nml%debug        = debug

  end subroutine get_namelist

  ! --------------------------------
  ! Print NWChem namelist settings
  ! --------------------------------
  subroutine print_namelist(nw_nml)

    use UtilitiesModule, only: Upcase

    implicit none
    type(nw_nml_type), intent(in) :: nw_nml

    write(6, '(a)')       '| &nw'
    write(6, '(2a)')      '|   method       = ', nw_nml%method
    write(6, '(2a)')      '|   basis        = ', nw_nml%basis
    write(6, '(2a)')      '|   guess        = ', nw_nml%guess
    write(6, '(a,i2)')    '|   scf_conv     = ', nw_nml%scf_conv
    write(6, '(a,i0)')    '|   ntpr         = ', nw_nml%ntpr
    write(6, '(a,i2)')    '|   num_threads  = ', nw_nml%num_threads
    write(6, '(a,i2)')    '|   debug        = ', nw_nml%debug
    write(6,'(a)')        '| /'

  end subroutine print_namelist

  ! -----------------------------
  ! Write input file for NWChem
  ! -----------------------------

  subroutine write_inpfile( inpfile, nqmatoms, qmcoords, qmtypes, &
       nclatoms, clcoords, nw_nml, do_grad, charge, spinmult )

    use ElementOrbitalIndex, only : elementSymbol

    implicit none

    character(len=*), intent(in)   :: inpfile
    integer, intent(in)            :: nqmatoms
    _REAL_,  intent(in)            :: qmcoords(:,:)
    integer, intent(in)            :: qmtypes(:)
    integer, intent(in)            :: nclatoms
    _REAL_,  intent(in)            :: clcoords(:,:)
    type(nw_nml_type), intent(in)  :: nw_nml
    logical, intent(in)            :: do_grad
    integer, intent(in)            :: charge, spinmult

    integer, parameter :: iunit = 351, tplunit = 352
    integer            :: i, ierr
    logical, save      :: first_call = .true.

    call debug_enter_function( 'write_inpfile', module_name, nw_nml%debug )

    open(iunit, file=inpfile, iostat=ierr)
    if ( ierr /= 0 ) then
      call sander_bomb('write_inpfile (qm2_extern_nw_module)', &
        'Error opening NWChem inpfile '//inpfile//' for writing', &
        'Will quit now.')
    end if

    ! Write title
    write (iunit,'(a)') 'title "NWChem run using SANDER external interface."'

    ! Charge
    write (iunit,'(a,i0,/)'),'charge ',charge

    ! Write QM atoms and coordinates
    ! More info here: http://www.nwchem-sw.org/index.php/Release61:Geometry
    write(iunit,'(a)') 'geometry units an noautoz nocenter noautosym noprint'
    do i = 1, nqmatoms
      write(iunit,'(a2,x,3f25.16)') elementSymbol(qmtypes(i)), qmcoords(1:3,i)
    end do
    write(iunit,'(a,/)') 'end'

    ! If QM/MM, write the external point charges
    if ( nclatoms > 0 ) then
       write ( iunit, '(a)' ) 'bq'
       if ( do_grad ) then
          write ( iunit, '(a)' ) 'force'
       end if
       do i = 1, nclatoms
          write ( iunit, '(4(x,f19.9))' ) clcoords(1:4,i)
       end do
       write ( iunit, '(a/)' ) 'end'
    end if

    ! Basis
    write(iunit,'(a)')       'basis'
    write(iunit,'(2a)')      '  * library ', trim(nw_nml%basis)
    write(iunit,'(a,/)')     'end'

    ! Method
    ! TODO; this is far too fagile

    ! HF Method
    ! More info here: http://www.nwchem-sw.org/index.php/Release61:Hartree-Fock_Theory_for_Molecules
    if (nw_nml%method == "hf") then
      write(iunit,'(a)')     'scf'
      ! Convergence
      write(iunit,'(a,i0)'), '  thresh 1e-',nw_nml%scf_conv
      ! Calculate all integrals "on-the-fly"
      write(iunit,'(a)')     '  direct '
      ! Do not print the MO vector coefficients; just too much data.
      write(iunit,'(a)')     '  noprint "final vectors analysis"'
      write(iunit,'(a,/)')   'end'
    end if

    ! DFT Method
    ! More info here: http://www.nwchem-sw.org/index.php/Density_Functional_Theory_for_Molecules
    if (nw_nml%method == "BLYP") then
      write(iunit,'(a)')     'dft'
      ! Convergence
      write(iunit,'(a,i1)')  '  convergence energy 1e-',nw_nml%scf_conv
      ! Grid
      write(iunit,'(a)'),    '  grid fine'
      write(iunit,'(a,a,a)') '  xc',' becke88', ' lyp'
      write(iunit,'(a,i1)'), '  mult ', spinmult
      write(iunit,'(a)'),    '  noio '
      ! Calculate all integrals "on-the-fly"
      write(iunit,'(a)')     '  direct '
      ! Do not print the MO vector coefficients; just too much data.
      write(iunit,'(a)')     '  noprint "final vectors analysis"'
      write(iunit,'(a,/)')   'end'
    end if

    ! One will need at least one of these task directives

    if (nw_nml%method == "hf") then
       if ( do_grad ) then
          write(iunit,'(a)') 'task scf gradient'
       else
          ! Don't calculate gradient anymore 
          write(iunit,'(a)') 'task scf energy'
       endif
    end if

    ! TODO, this is just not good enough.
    ! We  needs a generic list of DFT methods here
    ! i.e. if (nw_nml%method == "BLYP" | "B3LYP" | ) then
    ! or perhaps even change the "nw_nml_type" to include
    ! a functional entry?

    if (nw_nml%method == "BLYP") then
       if ( do_grad ) then
          write(iunit,'(a)') 'task dft gradient'
       else
          ! Don't calculate gradient anymore
          write(iunit,'(a)') 'task dft energy'
       endif
    end if
    
    close(iunit, iostat=ierr)

    if ( ierr /= 0 ) then
       call sander_bomb('write_inpfile (qm2_extern_nw_module)', &
            'Error closing NWChem runfile after writing', &
            'Will quit now.')
    end if

    first_call = .false.

   call debug_exit_function( 'write_inpfile', module_name, nw_nml%debug )

  end subroutine write_inpfile

  ! Parse the output of the nwchem.log file and the nwchem.bqforce.dat files
  !   * QM Energy
  !   * Gradient on QM atoms
  !   * Gradient on MM atoms
  subroutine read_results( nw_nml, datfile, bqforcefile, nqmatoms, escf, dxyzqm,&
       nclatoms, dxyzcl, do_grad, debug )

    implicit none

    type(nw_nml_type), intent(in) :: nw_nml
    character(len=*), intent(in)  :: datfile, bqforcefile
    integer, intent(in)           :: nqmatoms, nclatoms
    _REAL_, intent(out)           :: escf, dxyzqm(3,nqmatoms), &
                                     dxyzcl(3,nclatoms)
    logical, intent(in)           :: do_grad
    integer, intent(in)           :: debug

    integer :: ios, i
    integer :: nulli
    _REAL_  :: nullr
    integer, parameter :: iunit = 351
    character(len=256) :: read_buffer

    call debug_enter_function( 'read_results', module_name, debug )

    open(iunit, file=datfile, status='old', iostat=ios)
    if ( ios /= 0 ) then
      call sander_bomb('read_results (qm2_extern_nw_module)', &
        'Error opening NWChem log file '//datfile//' (expected in same dir as input file).', &
        'Will quit now')
    end if

    do
      read (iunit, '(a)', iostat = ios) read_buffer
      ! End of file; nothing left to read
      if (ios < 0) then
        exit
      end if

      if (nw_nml%method == "hf") then
        !         Total SCF energy =   -243.826232163838
        !      One-electron energy =   -684.067100443940
        !      Two-electron energy =    263.613349226485
        ! Nuclear repulsion energy =    176.627519053618

        if ( read_buffer(1:25) == '         Total SCF energy' ) then
          ! Read value after "=" sign
          read(read_buffer(index(read_buffer,'=')+1:),*) escf
        end if

      end if

      if (nw_nml%method == "BLYP") then
        !         Total DFT energy =      -74.827214054279
        !      One electron energy =     -115.784049745746
        !           Coulomb energy =       43.899915115726
        !    Exchange-Corr. energy =       -8.837093294400
        ! Nuclear repulsion energy =        5.894013870142

        if ( read_buffer(1:25) == '         Total DFT energy' ) then
          ! Read value after "=" sign
          read(read_buffer(index(read_buffer,'=')+1:),*) escf
        end if

      end if

      if( do_grad ) then
        ! Gradients
        !                         DFT ENERGY GRADIENTS
        !
        !    atom               coordinates                        gradient
        !                 x          y          z           x          y          z
        !   1 O       0.000000   0.000000   0.400871    0.000000   0.000000   0.151975
        !   2 H       2.004357   0.000000  -1.603486    0.073745   0.000000  -0.075987
        !   3 H      -2.004357   0.000000  -1.603486   -0.073745   0.000000  -0.075987

        if ( read_buffer(1:80) == '    atom               coordinates                        gradient' ) then
          ! Skip over this line
          read(iunit, '(a)') read_buffer
          ! Skip over next line (x,y,z heading)
          read(iunit, '(a)') read_buffer

          if ( nw_nml%debug > 1 ) then
            write (6,'(a)') 'read_results() - read in gradients:'
            write (6,'(a)') 'QM region:'
          end if

          do i = 1, nqmatoms
            ! Format from nwchem-6.0/src/gradients/grad_force.F line 1054
            ! format(1X,I3,1X,A4,2(1X,3(1X,F10.6)))

            ! MJW TODO; this is too dirty
            read (read_buffer, '(1X,I3,1X,A4,2(1X,3(1X,F10.6)))', iostat = ios) &
                                 nulli, nulli, nullr, nullr, nullr, dxyzqm(1:3,i)
            if ( nw_nml%debug > 1 ) then
              write(6,*) dxyzqm(1:3,i)
            end if

            ! Next line
            read(iunit, '(a)') read_buffer
          end do
        end if

      end if
    end do

    close(iunit)

    ! MM forces
    if ( do_grad .and. ( nclatoms > 0 ) ) then

       open(iunit, file=bqforcefile, status='old', iostat=ios)
       if ( ios /= 0 ) then
          call sander_bomb('read_results (qm2_extern_nw_module)', &
               'Error opening NWChem bqforce file '//bqforcefile//' (expected in same dir as input file).', &
        'Will quit now')
       end if

       read ( iunit, '(a)' ) read_buffer

       do i = 1, nclatoms
          read ( iunit, * ) dxyzcl(:,i)
       end do
       
    end if
     
    call debug_exit_function( 'read_results', module_name, debug )
 
  end subroutine read_results


end module qm2_extern_nw_module
