#include "../include/dprec.fh"
module qm2_extern_qc_module
! ----------------------------------------------------------------
! Interface for Q-Chem based QM MD 
!
! Currently supports:
! pure QM
! QM/MM with cutoff for QM-MM electrostatics 
!
! Implementation by
! Balachandar Kesavan, Roger Ouyang, Stephen Clark
! (SDSC REHS interns)
! under supervision of
! Andreas Goetz (SDSC)
! 
! Date: August 2012, August 2013
!
! Implementation extended by
! Evan Wildenhain, Pietro Sette, Prathyush Katukojwala
! (SDSC REHS interns)
! under the supervision of
! Andreas Goetz (SDSC)
!
! Date: July 2014
!
! ----------------------------------------------------------------
  use qm2_extern_util_module, only: debug_enter_function, debug_exit_function

  implicit none

  private
  public :: get_qc_forces

  character(len=*), parameter :: module_name = "qm2_extern_qc_module"

  type qc_nml_type
     character(len=20) :: method
     character(len=20) :: basis
     character(len=20) :: aux_basis
     character(len=20) :: exchange
     character(len=20) :: correlation
     character(len=20) :: guess
     integer :: scf_conv
     integer :: debug
     integer :: ntpr
     integer :: num_threads 
     logical :: dipole
     logical :: use_template
  end type qc_nml_type

contains

  ! --------------------------------------------
  ! Get QM energy and forces from Q-Chem
  ! --------------------------------------------
  subroutine get_qc_forces(do_grad, nstep, ntpr_default, id, nqmatoms, qmcoords,&
    qmtypes, nclatoms, clcoords, escf, dxyzqm, dxyzcl, charge, spinmult)

    use qm2_extern_util_module, only: check_installation, print_results, write_dipole
    use constants, only: CODATA08_AU_TO_KCAL, CODATA08_A_TO_BOHRS, ZERO
    use file_io_dat

    implicit none
    
    logical, intent(in) :: do_grad              ! Return gradient/not
    integer, intent(in) :: nstep                ! MD step number
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
    _REAL_              :: dipxyz(3), dipole    ! Dipole moment
    integer, intent(in) :: charge, spinmult     ! Charge and spin multiplicity

    type(qc_nml_type), save      :: qc_nml
    logical, save                :: first_call = .true.
    integer                      :: i
    integer                      :: printed =-1 ! Used to tell if we have printed this step yet 
                                                ! since the same step may be called multiple times
    character (len=150)          :: call_buffer
    character(len=6), save       :: program  = 'qchem'
    character(len=*), parameter  :: basename = 'qc_job'
    character(len=*), parameter  :: inpext = '.inp'
    character(len=*), parameter  :: logext = '.log'
    character(len=*), parameter  :: rstext = '.chk' ! Restart from checkpoint files
    character(len=*), parameter  :: dipext = '.dip'
    character(len=*), parameter  :: tplext = '.tpl'
    character(len=*), parameter  :: savext = '.sav'
    character(len=14)            :: inpfile, logfile, fortfile, dipfile, tplfile, datfile, savfile
    ! Need to prepend subdirectory if doing REMD, PIMD
    character(len=25)            :: subdir 

    ! for system call
    integer :: system
    integer :: stat

    ! assemble input - / output data filenames
    inpfile  = basename//trim(id)//inpext
    logfile  = basename//trim(id)//logext
    dipfile  = basename//trim(id)//dipext
    savfile  = basename//trim(id)//savext
    tplfile  = basename//tplext
    fortfile = 'fort.7'
    datfile  = 'efield.dat' 

    ! Setup on first call
    if ( first_call ) then
       first_call = .false.
       write (6,'(/,a,/)') '  >>> Running calculation with Q-Chem <<<'
       call get_namelist( ntpr_default, qc_nml )
       call print_namelist( qc_nml ) 
       ! Check for version of Q-Chem to use; store as 'program'
       call check_installation( program, id, .true., qc_nml%debug )

       write (6,'(80a)') ('-', i=1,80)
       write (6,'(a)') '   4.  RESULTS'
       write (6,'(80a)') ('-', i=1,80)
       ! Remove old inpfile, logfile, fort.7, and dipfile at the 
       ! beginning of a run so only the latest run is stored.
       stat = system('rm -f '//inpfile//' '//logfile//' '//fortfile//' '//dipfile//' '//datfile)
       if ( stat /= 0 ) then
          call sander_bomb('get_qc_forces (qm2_extern_qc_module)', & 
               'Error with system call (removing files)', &
               'Will quit now.')
       end if
    end if
    
    call write_inpfile( trim(inpfile), trim(tplfile), &
      nqmatoms, qmcoords, qmtypes, nclatoms, clcoords, qc_nml, do_grad, charge, spinmult )

    ! Run Q-Chem
    ! Search in subdir for logfile if doing PIMD
    subdir=''
    call_buffer=''
    if(trim(id)/='') then 
      subdir='./'//trim(id)//'/'
      call_buffer=' mkdir -p '//trim(subdir)//'; cd '//trim(subdir)//'; mv ../'//inpfile//' .;'
    end if

    call_buffer = trim(call_buffer)//' '//program

    ! parallelization
    if (qc_nml%num_threads > 1) then
       write(call_buffer, '(2a,i0)') &
            trim(call_buffer), ' -np ', qc_nml%num_threads
    end if

    ! input/output/scratch files
    call_buffer = trim(call_buffer)//' '//inpfile//' '//logfile//' '//savfile
       
    if ( qc_nml%debug > 0 ) then
       write (6,'(2a)') 'call_buffer=', trim(call_buffer)
    end if
    
    stat = system(trim(call_buffer))

    if ( stat /= 0 ) then
       call sander_bomb('get_qc_forces (qm2_extern_qc_module)', & 
            'Error with system call (executing Q-Chem)', &
            'Will quit now.')
    end if

    ! Retrieve data from Q-Chem .log and efield.dat files
    ! Will output data to escf and dxyqm for pure QM or mechanical embedding QM/MM
    ! Search in subdir for logfile if doind PIMD
    ! Otherwise, search current directory
    call read_results( trim(subdir)//trim(logfile), trim(subdir)//trim(datfile), & 
      nqmatoms, escf, dxyzqm, nclatoms, dxyzcl, dipxyz, dipole, qc_nml%dipole)

    ! Write dipole moment to dipfile
    if ( qc_nml%ntpr > 0 .and. mod(nstep, qc_nml%ntpr) == 0 ) then
       if ( printed /= nstep .and. qc_nml%dipole ) then
          call write_dipole(trim(dipfile), dipxyz, dipole, qc_nml%debug)
          printed = nstep
       end if
    end if

    ! F = E*q to get gradients
    do i = 1, nclatoms
       dxyzcl(:,i) = -dxyzcl(:,i) * clcoords(4,i)
    end do
 
    if ( do_grad ) then
       ! Convert Hartree/Bohr -> kcal/(mol*A)  
       dxyzqm(:,:) = dxyzqm(:,:) * CODATA08_AU_TO_KCAL * CODATA08_A_TO_BOHRS
       if ( nclatoms > 0 ) then
          dxyzcl(:,:) = dxyzcl(:,:) * CODATA08_AU_TO_KCAL * CODATA08_A_TO_BOHRS
       end if
    else
      dxyzqm = ZERO
      if ( nclatoms > 0 ) dxyzcl = ZERO
    end if
    escf = escf * CODATA08_AU_TO_KCAL 
    
    call print_results( 'qm2_extern_qc_module', escf, nqmatoms, dxyzqm,&
      qc_nml%debug, nclatoms, dxyzcl )
 
    ! Save copy of last input and log files
    stat = system('mv '//trim(subdir)//inpfile//' '//trim(subdir)//'old.'//inpfile)
    stat = stat + system('mv '//trim(subdir)//logfile//' '//trim(subdir)//'old.'//logfile)
    stat = stat + system('mv '//trim(subdir)//datfile//' '//trim(subdir)//'old.'//datfile)
    if ( stat /= 0 ) then
      call sander_bomb('get_qc_forces (qm2_extern_qc_module)', & 
        'Error with system call (moving / removing files)', &
        'Will quit now.')
    end if

  end subroutine get_qc_forces

  ! ---------------------------------------------
  ! Read Q-Chem namelist values from file mdin,
  ! use default values if none are present.
  ! ---------------------------------------------
 
  subroutine get_namelist(ntpr_default, qc_nml)

    use UtilitiesModule, only : Upcase

    implicit none
    integer, intent(in) :: ntpr_default
    type(qc_nml_type), intent(out) :: qc_nml

    character(len=20) :: method, basis, exchange, correlation, aux_basis, guess
    character(len=20) :: default_method, default_basis, default_mp2_basis, default_aux_basis
    integer :: debug
    integer :: scf_conv, ntpr, num_threads, dipole, use_template
    namelist /qc/ method, basis, exchange, correlation, aux_basis, guess, scf_conv, ntpr, &
      num_threads, debug, dipole, use_template

    integer :: ierr

    call debug_enter_function( 'get_namelist', module_name, qc_nml%debug )

    ! Set default values for qc namelist values
    method            = ''
    default_method    = 'BLYP'
    basis             = ''
    default_basis     = '6-31G*'
    default_mp2_basis = 'cc-pVDZ'
    aux_basis         = ''
    default_aux_basis = 'rimp2-cc-pVDZ'
    exchange          = ''
    correlation       = ''
    guess             = 'read' ! have qchem read restart files by default
    scf_conv          = 6 
    ntpr              = ntpr_default
    num_threads       = 1
    dipole            = 0
    use_template      = 0
    debug             = 0

    ! Read namelist
    rewind 5
    read(5,nml=qc,iostat=ierr)

    if ( ierr > 0 ) then
       call sander_bomb('get_namelist (qm2_extern_qc_module)', &
            '&qc namelist read error', &
            'Please check your input.')
    else if ( ierr < 0 ) then
       write(6,'(a/a)') '&qc namelist read encountered end of file', &
            'Please check your input if the calculation encounters a problem'
    end if

    ! Assign namelist values to qc_nml data type
    ! QM method
    if ( (method == '') .and. (exchange == '') .and. (correlation == '') &
         .and. (use_template == 0) ) then
       qc_nml%method = default_method
    else
       qc_nml%method = method
    end if

    if ( (method /= '') .and. ( (exchange /= '') .or. (correlation /= '') ) ) then
       call sander_bomb('get_namelist (qm2_extern_qc_module)', &
            'You cannot specify "method" and "exchange", "correlation" at the same time.', &
            'Will quit now.')
    end if
    qc_nml%exchange     = exchange
    qc_nml%correlation  = correlation

    ! Basis set
    if ( basis == '' ) then
       if ( (Upcase(qc_nml%method) == 'MP2') .or. &
            (Upcase(qc_nml%method) == 'RIMP2') .or. &
            (Upcase(qc_nml%method) == 'RI-MP2') )  then
          qc_nml%basis = default_mp2_basis
       else
          qc_nml%basis = default_basis
       end if
    else
       qc_nml%basis        = basis
    end if

    ! Aux basis set
    if ( (aux_basis == '') .and. &
         ( (Upcase(qc_nml%method) == 'RIMP2') .or. &
           (Upcase(qc_nml%method) == 'RI-MP2') ) ) then
       qc_nml%aux_basis = default_aux_basis
    else
       qc_nml%aux_basis = aux_basis
    end if

    qc_nml%guess        = guess
    qc_nml%scf_conv     = scf_conv
    qc_nml%ntpr         = ntpr
    qc_nml%num_threads  = num_threads
    qc_nml%debug        = debug

    if ( dipole == 0 ) then
       qc_nml%dipole = .false.
    else if ( dipole == 1 ) then
    !  qc_nml%dipole = .true. <---- not supported by Q-Chem interface as of 2014-07-21
       call sander_bomb('get_namelist (qm2_extern_qc module)', &
            '&qc dipole value not allowed', &
            'Q-Chem interface does not currently support reading dipole.')
    else
       call sander_bomb('get_namelist (qm2_extern_qc_module)', &
            '&qc dipole value not allowed', &
            'Please check your input. dipole can only be 0 or 1.')
    end if

    if ( use_template == 0 ) then
       qc_nml%use_template = .false.
    else if ( use_template == 1 ) then
       qc_nml%use_template = .true.
    else
       call sander_bomb('get_namelist (qm2_extern_qc_module)', &
            '&qc use_template value not allowed', &
            'Please check your input. use_template can only be 0 or 1.')
    end if

    call debug_exit_function( 'get_namelist', module_name, qc_nml%debug )

  end subroutine get_namelist

  ! --------------------------------
  ! Print Q-Chem namelist settings
  ! --------------------------------
  subroutine print_namelist(qc_nml)

    implicit none
    type(qc_nml_type), intent(in) :: qc_nml

    write(6, '(a)')       '| &qc'
    write(6, '(2a)')      '|   method       = ', qc_nml%method
    write(6, '(2a)')      '|   basis        = ', qc_nml%basis
    write(6, '(2a)')      '|   exchange     = ', qc_nml%exchange
    write(6, '(2a)')      '|   correlation  = ', qc_nml%correlation
    write(6, '(2a)')      '|   aux_basis    = ', qc_nml%aux_basis
    write(6, '(2a)')      '|   scf_guess    = ', qc_nml%guess   
    write(6, '(a,i2)')    '|   scf_conv     = ', qc_nml%scf_conv
    write(6, '(a,i0)')    '|   ntpr         = ', qc_nml%ntpr
    write(6, '(a,i2)')    '|   num_threads  = ', qc_nml%num_threads
    write(6, '(a,i2)')    '|   debug        = ', qc_nml%debug
    write(6, '(a,l)')     '|   dipole       = ', qc_nml%dipole
    write(6, '(a,l)')     '|   use_template = ', qc_nml%use_template
    write(6,'(a)')        '| /'

  end subroutine print_namelist

  ! -----------------------------
  ! Write input file for Q-Chem
  ! -----------------------------
  subroutine write_inpfile( inpfile, tplfile, nqmatoms, qmcoords,&
    qmtypes, nclatoms, clcoords, qc_nml, do_grad, charge, spinmult )

    use ElementOrbitalIndex, only : elementSymbol
    use UtilitiesModule, only     : Upcase
    implicit none

    character(len=*), intent(in)   :: inpfile, tplfile
    integer, intent(in)            :: nqmatoms
    _REAL_,  intent(in)            :: qmcoords(:,:)
    integer, intent(in)            :: qmtypes(:)
    integer, intent(in)            :: nclatoms
    _REAL_,  intent(in)            :: clcoords(:,:)
    type(qc_nml_type), intent(in) :: qc_nml
    logical, intent(in)            :: do_grad
    integer, intent(in)            :: charge, spinmult

    integer, parameter :: iunit = 351, tplunit = 352
    integer            :: i, ierr
    logical, save      :: first_call = .true.

    ! for system call
    integer :: system
    integer :: stat

    call debug_enter_function( 'write_inpfile', module_name, qc_nml%debug )

    open(iunit, file=inpfile, iostat=ierr)
    if ( ierr /= 0 ) then
       call sander_bomb('write_inpfile (qm2_extern_qc_module)', &
            'Error opening Q-Chem inpfile '//inpfile//' for writing', &
            'Will quit now.')
    end if

    ! Start writing molecule
    write(iunit,'(a)') '$molecule'
    
    ! Write charge and spin multiplicity
    write (iunit,'(i0,a,i0)') charge,' ', spinmult
    
    ! Write QM atoms and coordinates
    do i = 1, nqmatoms
       write(iunit,'(a2,1x,3f25.16)') elementSymbol(qmtypes(i)), qmcoords(1:3,i)
    end do
    
    write(iunit,'(a)') '$end'
    write(iunit,'(a)')
    
    ! When electrostatic embedding QM/MM is in use 
    ! Write MM coordinates with point charges
    if ( nclatoms > 0 ) then
       write(iunit,'(a)') '$external_charges'
       do i = 1, nclatoms
          write(iunit,'(4f21.16)') clcoords(:,i) 
       end do
       write(iunit,'(a)') '$end'
       write(iunit,'(a)')
    end if
    
    ! Write REM options
    write(iunit,'(a)')    '$rem'
    write(iunit,'(a)')    'SYMMETRY false' 
    write(iunit,'(a)')    'SYM_IGNORE true'
    if (do_grad) then
       write(iunit,'(2a)')   'JOBTYPE ', 'force' 
    else
       write(iunit,'(2a)')   'JOBTYPE ', 'sp' 
    end if
    write(iunit,'(a)') 'QM_MM true'

    if ( qc_nml%use_template ) then 
       ! Template file for QM method
       close(iunit)
       stat = system('cat '//trim(tplfile)//' >>  '//inpfile)
       if ( stat /= 0 ) then
          call sander_bomb('write_inpfile (qm2_extern_qc_module)', & 
               'Error with system call (cat-ing files)', &
               'Will quit now.')
       end if
       open(iunit, file=inpfile, access='append', iostat=ierr)
       if ( ierr /= 0 ) then
          call sander_bomb('write_inpfile (qm2_extern_qc_module)', &
               'Error opening Q-Chem inpfile '//inpfile//' for appending', &
               'Will quit now.')
       end if
    else  
       ! QM Method
       if (Upcase(qc_nml%method) == 'MP2') then   ! Default values for MP2
          write(iunit,'(a)')   'EXCHANGE HF'
          write(iunit,'(a)')   'CORRELATION MP2'
       else if (Upcase(qc_nml%method) == 'RIMP2' .or. &
                Upcase(qc_nml%method) == 'RI-MP2') then ! Default values for RIMP2
          write(iunit,'(a)')   'EXCHANGE HF'
          write(iunit,'(a)')   'CORRELATION RIMP2'
       else if(Upcase(qc_nml%method) == 'BLYP') then ! Default settings for BLYP
          write(iunit,'(a)')   'EXCHANGE BECKE'
          write(iunit,'(a)')   'CORRELATION LYP'
       else
          if (qc_nml%exchange /= '') then
             write(iunit,'(2a)')   'EXCHANGE ', qc_nml%exchange
          else if (qc_nml%method /= '') then
             ! We will assume this is an xc functional known to Q-Chem
             write(iunit,'(2a)')   'EXCHANGE ', qc_nml%method
          end if
          if(qc_nml%correlation /= '') then
             write(iunit,'(2a)')   'CORRELATION ', qc_nml%correlation
          end if
       end if

       ! Basis set
       write(iunit,'(2a)')   'BASIS ', qc_nml%basis
       if ( qc_nml%aux_basis /= '' ) then
          write(iunit,'(2a)')   'AUX_BASIS ', qc_nml%aux_basis
       end if

       write(iunit,'(a,i0)') 'SCF_CONVERGENCE ', qc_nml%scf_conv
    end if

    if ( .not. first_call .and. (Upcase(trim(qc_nml%guess)) == 'READ') ) then
       write(iunit,'(2a)') 'SCF_GUESS ', qc_nml%guess
    end if

    write(iunit,'(a)')    '$end'
    
    close(iunit, iostat=ierr)
    if ( ierr /= 0 ) then
       call sander_bomb('write_inpfile (qm2_extern_qc_module)', &
            'Error closing Q-Chem runfile after writing', &
            'Will quit now.')
    end if

    first_call = .false.

    call debug_exit_function( 'write_inpfile', module_name, qc_nml%debug )

  end subroutine write_inpfile
  
  ! QM or QM/MM results
  subroutine read_results( logfile, datfile, nqmatoms, escf, dxyzqm,&
    nclatoms, dxyzcl, dipxyz, dipole, do_dipole)

    use constants, only: CODATA08_AU_TO_KCAL, CODATA08_A_TO_BOHRS
    use UtilitiesModule, only     : Upcase
    implicit none
    
    character(len=*), intent(in) :: logfile, datfile
    integer, intent(in)          :: nqmatoms, nclatoms
    _REAL_, intent(out)          :: escf, dxyzqm(3,nqmatoms), & 
                                    dxyzcl(3,nclatoms) ! dxyzcl will return containing the electric field at x,y,z
    _REAL_, intent(out)          :: dipxyz(3), dipole ! return in Debye
    logical, intent(in)          :: do_dipole

    
    integer :: ios, i
    integer, parameter :: iunit = 351
    character(len=120) :: read_buffer
    _REAL_  :: echg ! Charge-Charge Energy

    open(iunit, file=logfile, status='old', iostat=ios)
    if ( ios /= 0 ) then
      call sander_bomb('read_results (qm2_extern_qc_module)', &
        'Error opening Q-Chem log file '//logfile//' (expected in same dir as input file).', &
        'Will quit now')
    end if

    escf = 0.0d0
    echg = 0.0d0

    do
       read (iunit, '(a)', iostat = ios) read_buffer
       ! End of file; nothing left to read
       if (ios < 0) then
          ! Subtract charge-charge energy at the end
          exit
       end if

       ! Point charge self energy
       if ( read_buffer(1:21) == ' Charge-charge energy' ) then
          read(read_buffer(index(read_buffer,'=')+1:),'(f17.10)') echg
       end if
       ! Find SCF energy line
       if ( read_buffer(1:13) == ' SCF   energy' ) then
          ! Read value after the first "=" sign
          read(read_buffer(index(read_buffer,'=')+2:),*) escf
       end if
       ! Find RIMP2 energy line if it exists
       if ( read_buffer(1:19) == 'RI-MP2 TOTAL ENERGY' ) then
          read(read_buffer(index(read_buffer,'=')+1:),'(f27.10)') escf
       end if
       ! Find MP2 energy line if it exists
       if ( read_buffer(1:32) == '        MP2         total energy' ) then
          read(read_buffer(index(read_buffer,'=')+1:),'(f21.10)') escf
       end if

       if ( do_dipole ) then
          ! Dipole moment
          if (index(read_buffer,'Dipole Moment') > 0) then
             read (iunit, '(a)') read_buffer
             read(read_buffer(11:),'(f13.4)') dipxyz(1) ! already in Debye
             read(read_buffer(31:),'(f13.4)') dipxyz(2)
             read(read_buffer(51:),'(f13.4)') dipxyz(3)
             read (iunit, '(a)') read_buffer
             read(read_buffer(11:),'(f13.4)') dipole
          end if
       end if        
    end do
    close(iunit)

    ! Remove charge self-energy from total energy
    escf = escf - echg
    
    ! Open file for reading gradients
    open(iunit, file=datfile, status='old', iostat=ios)
    if ( ios /= 0 ) then
       call sander_bomb('read_results (qm2_extern_qc_module)', &
            'Error opening Q-Chem log file "efield.dat" (expected in same dir as input file).', &
            'Will quit now')
    end if

    ! QM/MM with electrostatic embedding: E-field at point charges
    if( nclatoms > 0 ) then       
       do i = 1, nclatoms
          read (iunit, '(3(f22.16))', iostat = ios) dxyzcl(:,i)
          if (ios < 0) then
             call sander_bomb('read_results (qm2_extern_qc_module)', &
                  'Error reading e-field at point charge positions.',&
                  'Will quit now')
             exit
          end if
       end do
    end if

    ! read QM gradients
    do i = 1, nqmatoms
       read(iunit, '(3(f25.20))', iostat = ios) dxyzqm(:,i)
       if ( ios < 0 ) then
          call sander_bomb('read_results (qm2_extern_qc_module)',&
               'Error reading QM gradients.',&
               'Will quit now')
       end if
    end do
    
    close (iunit)

  end subroutine read_results
  
end module qm2_extern_qc_module
