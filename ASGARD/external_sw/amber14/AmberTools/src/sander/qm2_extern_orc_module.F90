#include "../include/dprec.fh"
module qm2_extern_orc_module
! ----------------------------------------------------------------
! Interface for Orca based QM and QM/MM MD 
!
! Currently supports:
! pure QM
! QM/MM with cutoff for QM-MM electrostatics under periodic
! boundary conditions
!
! Author: Matthew Clark
! Based on qm2_extern_tc_module.f
!
! Date: July 2011
!
! ----------------------------------------------------------------

  use qm2_extern_util_module, only: debug_enter_function, debug_exit_function

  implicit none

  private
  public :: get_orc_forces

  character(len=*), parameter :: module_name = "qm2_extern_orc_module"

  type orc_nml_type
     character(len=20) :: basis
     character(len=20) :: cbasis
     character(len=20) :: jbasis
     character(len=20) :: method
     character(len=20) :: convkey
     character(len=20) :: guess
     integer :: scfconv
     integer :: grid
     integer :: finalgrid
     integer :: maxiter
     integer :: maxcore
     integer :: ntpr
     integer :: num_threads
     integer :: debug
     logical :: dipole
     logical :: use_template
     
     ! Deprecated
     integer :: charge
     integer :: spinmult
  end type orc_nml_type


contains

  ! --------------------------------------
  ! Get QM energy and forces from Orca
  ! --------------------------------------
  subroutine get_orc_forces( do_grad, nstep, ntpr_default, id, &
       nqmatoms, qmcoords, qmtypes, nclatoms, clcoords, escf, dxyzqm, dxyzcl, &
       charge, spinmult )

    use qm2_extern_util_module, only: print_results, check_installation, write_dipole
    use constants,              only: CODATA08_AU_TO_KCAL, CODATA08_A_TO_BOHRS, ZERO

    logical, intent(in) :: do_grad              ! Return gradient/not
    integer, intent(in) :: nstep                ! MD step number
    integer, intent(in) :: ntpr_default         ! frequency of printing
    character(len=3), intent(in) :: id          ! ID number for PIMD or REMD
    integer, intent(in) :: nqmatoms             ! Number of QM atoms
    _REAL_,  intent(in) :: qmcoords(3,nqmatoms) ! QM coordinates
    integer, intent(in) :: qmtypes(nqmatoms)    ! QM atom types (nuclear charge in au)
    integer, intent(in) :: nclatoms             ! Number of MM atoms
    _REAL_,  intent(in) :: clcoords(4,nqmatoms) ! MM coordinates and charges in au
    _REAL_, intent(out) :: escf                 ! SCF energy
    _REAL_, intent(out) :: dxyzqm(3,nqmatoms)   ! SCF QM force
    _REAL_, intent(out) :: dxyzcl(3,nclatoms)   ! SCF MM force
    _REAL_              :: dipxyz(3), dipole    ! Dipole moment
    integer, intent(in) :: charge, spinmult     ! Charge and spin multiplicity

    type(orc_nml_type), save :: orc_nml
    logical, save :: first_call = .true.
    integer :: i
    integer :: printed =-1 ! Used to tell if we have printed this step yet 
                           ! since the same step may be called multiple times
    character(len=512) :: call_buffer
    character(len=*), parameter  :: program='orca'
    character(len=*),  parameter :: basename = 'orc_job'
    character(len=*),  parameter :: runext = '.inp'
    character(len=*),  parameter :: datext = '.dat'
    character(len=*),  parameter :: gbwext = '.gbw'
    character(len=*),  parameter :: dipext = '.dip'
    character(len=*),  parameter :: tplext = '.tpl'
    character(len=14) :: inpfile, datfile, gbwfile, dipfile, crdfile, chgfile, tplfile
    ! Need to prepend subdirectory if doing REMD, PIMD or multi-region QM/MM. 
    !   This is triggered if 'id' is defined (not empty). 
    character(len=25) :: subdir 
    character(len=256), save :: orca_path
    
    ! assemble input - / output data filenames
    inpfile = basename//trim(id)//runext 
    datfile = basename//trim(id)//datext
    gbwfile = basename//trim(id)//gbwext
    dipfile = basename//trim(id)//dipext 
    crdfile = 'inpfile'//trim(id)//'.xyz'
    chgfile = 'ptchrg'//trim(id)//'.xyz'
    tplfile = basename//tplext 
    
    ! Setup on first program call
    if ( first_call ) then
      first_call = .false.
      write (6,'(/,a,/)') '   >>> Running QM calculation with Orca <<<'
      call get_namelist( ntpr_default, orc_nml )
      call check_installation( program, id, .true., orc_nml%debug, path=orca_path )
      call print_namelist( orc_nml )

      write (6,'(80a)') ('-', i=1,80)
      write (6,'(a)') '   4.  RESULTS'
      write (6,'(80a)') ('-', i=1,80)
      call system('rm -f '//dipfile//' '//gbwfile)
    end if

    call system('rm -f '//inpfile)
    call write_inpfile( trim(inpfile), trim(crdfile), trim(chgfile), trim(tplfile), &
         nqmatoms, qmcoords, qmtypes, nclatoms, clcoords, &
         orc_nml, do_grad, charge, spinmult )

    call_buffer=''
    subdir=''
    if (trim(id)/='') then
      subdir='./'//trim(id)//'/'
      call_buffer=' mkdir -p '//trim(subdir)//&
                  '; cd '//trim(subdir)//&
                  '; mv ../'//inpfile//' .'//&
                  '; mv ../'//crdfile//' .;'
      if (nclatoms > 0) then
        call_buffer=trim(call_buffer)//' mv ../'//chgfile//' .;'
      end if
    end if

    ! remove gbw file if we don't want to re-use MOs
    if ( trim(orc_nml%guess) /= 'read' ) then
       call_buffer = trim(call_buffer)//' rm -f '//trim(gbwfile)//';'
    end if

    ! Run Orca with file inpfile
    call_buffer = trim(call_buffer)//' rm -f '//datfile//';'
    write(call_buffer,'(2a)') &
         trim(call_buffer), trim(orca_path)//'/'//program//' '//trim(inpfile)//' > '//datfile

    if ( orc_nml%debug > 0 ) then
       write (6,'(2a)') 'call_buffer=', trim(call_buffer)
    end if
    
    call system(trim(call_buffer))

    ! If working in a subdirectory, move everything back for reading
    if (trim(id)/='') then
      call system('mv '//trim(subdir)//'* .;')
    end if

    ! Read Orca results
    call read_results(trim(datfile), basename//trim(id), nqmatoms, nclatoms, &
    escf, dxyzqm, dxyzcl, dipxyz, dipole, do_grad, orc_nml%dipole, orc_nml%debug )
 
    ! Write dipole moment to dipfile
    if ( orc_nml%ntpr > 0 .and. mod(nstep, orc_nml%ntpr) == 0 ) then
      if ( printed /= nstep .and. orc_nml%dipole ) then
        call write_dipole(trim(dipfile), dipxyz, dipole, orc_nml%debug)
        printed = nstep
      end if
    end if

    ! Backup old input and datfiles (if in subdirectory, copied out earlier)
    call system('mv '//trim(inpfile)//' '//'old.'//inpfile)
    call system('mv '//trim(datfile)//' '//'old.'//datfile)
    call system('mv '//trim(crdfile)//' '//'old.'//crdfile)
    if (nclatoms > 0) then
       call system('mv '//trim(chgfile)//' '//'old.'//chgfile)
    end if

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

    call print_results( 'qm2_extern_orc_module', escf, nqmatoms, dxyzqm,&
      orc_nml%debug, nclatoms, dxyzcl )

  end subroutine get_orc_forces

  ! -----------------------------------------------
  ! Read Orca orc namelist values from file mdin,
  ! use default values if none are present.
  ! -----------------------------------------------
  subroutine get_namelist( ntpr_default, orc_nml )

    use UtilitiesModule, only: Upcase
    implicit none

    integer, intent(in) :: ntpr_default
    type(orc_nml_type), intent(out) :: orc_nml
    character(len=20):: basis, cbasis, jbasis, method, convkey, guess
    integer :: scfconv, grid, finalgrid, maxiter, maxcore, &
                 ntpr,num_threads, debug, dipole, use_template
    integer :: charge, spinmult ! deprecated
    namelist /orc/ basis, cbasis, jbasis, method, convkey, guess, &
      scfconv, grid, finalgrid, &
      maxiter, maxcore, ntpr, num_threads, debug, &
      dipole, use_template,&
      charge, spinmult
    integer :: ierr

    ! Default values
    basis        = 'SV(P)'
    cbasis       = 'NONE'
    jbasis       = 'NONE'
    method       = 'blyp'
    convkey      = 'verytightscf'
    guess        = 'read'
    scfconv      = -1
    grid         = 4
    finalgrid    = 6
    maxiter      = 100
    maxcore      = 1024
    ntpr         = ntpr_default
    num_threads  = 1
    debug        = 0
    dipole       = 0
    use_template = 0

    ! These are now deprecated and should be specified in the &qmmmm namelist
    charge   = -351
    spinmult = -351

    ! Read namelist
    rewind 5
    read(5,nml=orc,iostat=ierr)

    if ( ierr > 0 ) then
       call sander_bomb('get_namelist (qm2_extern_orc_module)', &
            '&orc namelist read error', &
            'Please check your input.')
    else if ( ierr < 0 ) then
       write(6,'(a,/,a)') '&orc namelist read encountered end of file', &
            'Please check your input if the calculation encounters a problem'
    end if

    if ( charge /= -351 .or. spinmult /= -351 ) then
      call sander_bomb('get_namelist (qm2_extern_orc_module)', &
        'The charge and spin keywords are deprecated', &
        'Please specify charge (qmcharge) and spin multiplicity (spin) in the &qmmm namelist.')
    end if

    ! Assign namelist values to orc_nml data type
    orc_nml%basis       = Upcase(basis)
    orc_nml%cbasis      = Upcase(cbasis)
    orc_nml%jbasis      = Upcase(jbasis)
    orc_nml%method      = method
    orc_nml%convkey     = convkey
    orc_nml%guess       = guess
    orc_nml%scfconv     = scfconv
    orc_nml%grid        = grid
    orc_nml%finalgrid   = finalgrid
    orc_nml%maxiter     = maxiter
    orc_nml%maxcore     = maxcore
    orc_nml%ntpr        = ntpr
    orc_nml%num_threads = num_threads
    orc_nml%debug       = debug

    if ( dipole == 0 ) then
       orc_nml%dipole = .false.
    else if ( dipole == 1 ) then
       orc_nml%dipole = .true.
    else
       call sander_bomb('get_namelist (qm2_extern_orc_module)', &
            '&orc dipole value not allowed', &
            'Please check your input. dipole can only be 0 or 1.')
    end if

    if ( use_template == 0 ) then
       orc_nml%use_template = .false.
    else if ( use_template == 1 ) then
       orc_nml%use_template = .true.
    else
       call sander_bomb('get_namelist (qm2_extern_orc_module)', &
            '&orc use_template value not allowed', &
            'Please check your input. use_template can only be 0 or 1.')
    end if


  end subroutine get_namelist

  ! --------------------------------
  ! Print Orca namelist settings
  ! --------------------------------
  subroutine print_namelist( orc_nml )

    implicit none
    type(orc_nml_type), intent(in) :: orc_nml

    write(6, '(/,a)')      '| "NONE" or "-1" indicates ORCA default values'
    write(6, '(/,a)')      '| &orc'
    write(6, '(2a)')       '|   basis        = ', orc_nml%basis
    write(6, '(2a)')       '|   cbasis       = ', orc_nml%cbasis
    write(6, '(2a)')       '|   jbasis       = ', orc_nml%jbasis
    write(6, '(2a)')       '|   method       = ', orc_nml%method
    write(6, '(2a)')       '|   convkey      = ', orc_nml%convkey
    write(6, '(2a)')       '|   guess        = ', orc_nml%guess
    write(6, '(a,i4)')     '|   scfconv      = ', orc_nml%scfconv
    write(6, '(a,i4)')     '|   grid         = ', orc_nml%grid
    write(6, '(a,i4)')     '|   finalgrid    = ', orc_nml%finalgrid
    write(6, '(a,i4)')     '|   maxiter      = ', orc_nml%maxiter
    write(6, '(a,i0)')     '|   maxcore      = ', orc_nml%maxcore
    write(6, '(a,i0)')     '|   ntpr         = ', orc_nml%ntpr
    write(6, '(a,i0)')     '|   num_threads  = ', orc_nml%num_threads
    write(6, '(a,i2)')     '|   debug        = ', orc_nml%debug
    write(6, '(a,l)')      '|   dipole       = ', orc_nml%dipole
    write(6, '(a,l)')      '|   use_template = ', orc_nml%use_template
    write(6,'(a)')         '| /'

  end subroutine print_namelist

  ! ---------------------------------
  ! Write the input file for Orca
  ! ---------------------------------
  subroutine write_inpfile( inpfile, crdfile, chgfile, tplfile, &
       nqmatoms, qmcoords, qmtypes, nclatoms, clcoords, &
       orc_nml, do_grad, charge, spinmult )

    use ElementOrbitalIndex,    only: elementSymbol
    use qm2_extern_util_module, only: write_chgfile
    use UtilitiesModule,        only: Upcase

    implicit none

    character (len=*) , intent(in) :: inpfile
    character (len=*) , intent(in) :: crdfile
    character (len=*) , intent(in) :: chgfile
    character (len=*) , intent(in) :: tplfile
    integer           , intent(in) :: nqmatoms
    _REAL_            , intent(in) :: qmcoords(:,:)
    integer           , intent(in) :: qmtypes(:)
    integer           , intent(in) :: nclatoms
    _REAL_            , intent(in) :: clcoords(:,:)
    type(orc_nml_type), intent(in) :: orc_nml
    logical           , intent(in) :: do_grad
    integer           , intent(in) :: charge, spinmult

    character(len=256) :: simple_input
    integer            :: i, ios
    integer, parameter :: iurun = 10
    logical, save :: first_call = .true.
    logical, save :: dft

    call debug_enter_function( 'write_inpfile', module_name, orc_nml%debug )

    if ( first_call ) then
       ! We assume that we are doing DFT if we do neither HF or MP2
       dft = (index(Upcase(orc_nml%method), 'MP2') == 0) .and. ( index(Upcase(orc_nml%method), 'HF') == 0 )
       first_call = .false.
    end if

    ! Are we using a template input file?
    if ( orc_nml%use_template ) then
      call system('cp '//trim(tplfile)//' '//inpfile)
    end if

    open(iurun, file=inpfile, position='append', iostat=ios)
    if ( ios > 0 ) then
      call sander_bomb('write_inpfile (qm2_extern_orc_module)', &
           'Error opening Orca input file '//inpfile//' for writing', &
           'Will quit now')
    end if

    ! If no template is used, write appropriate keywords to input file
    if ( .not. orc_nml%use_template ) then 

       write(iurun, '(a,/,a)')  '# Run using SANDER file-based interface for Orca','#'
       ! Number of processors for parallel runs
       write(iurun, '(a,i0,a)')      '%pal nprocs ',orc_nml%num_threads,' end'

       ! Method, basis set, SCF convergence via simplified input
       write(simple_input,'(2a)')    '!',trim(orc_nml%method)//' '//trim(orc_nml%basis)
       if ( orc_nml%cbasis /= 'NONE' ) then
          write(simple_input,'(a)')   trim(simple_input)//' '//trim(orc_nml%cbasis)
       endif
       if ( orc_nml%jbasis /= 'NONE' ) then
          write(simple_input,'(a)')   trim(simple_input)//' '//trim(orc_nml%jbasis)
       endif
       if ( orc_nml%convkey /= 'NONE' ) then
          write(simple_input,'(a)')   trim(simple_input)//' '//trim(orc_nml%convkey)
       endif
       if ( orc_nml%scfconv > 4 ) then
          if ( orc_nml%scfconv < 10 ) then
             write(simple_input,'(a,i1)')  trim(simple_input)//' SCFCONV',orc_nml%scfconv
          else
             write(simple_input,'(a,i2)')  trim(simple_input)//' SCFCONV',orc_nml%scfconv
          endif
       endif
       write(iurun, '(a)')           trim(simple_input)

       if ( dft ) then
          ! DFT grid info
          if ( (orc_nml%grid > 0) .or. (orc_nml%finalgrid > 0) ) then
             write(iurun, '(a)')           '%method'
             if ( orc_nml%grid > 0 ) then
                write(iurun, '(a,i1)')     '  grid ', orc_nml%grid
             end if
             if ( orc_nml%finalgrid > 0 ) then
                write(iurun, '(a,i1)')     '  finalgrid ', orc_nml%finalgrid
             end if
             write(iurun, '(a)')           'end'
          end if
       end if

       ! SCF iterations
       write(iurun, '(a)')           '%scf'
       write(iurun, '(a,i0)')        '  maxiter ',orc_nml%maxiter
       write(iurun, '(a)')           'end'

       ! Maximum RAM
       write(iurun, '(a,i0)')        '%MaxCore ',orc_nml%maxcore

    end if

    ! Gradient (default) or single point (for post processing of trajectories)
    if ( do_grad ) then
       write(iurun, '(2a)')        '! ENGRAD'
    else
       write(iurun, '(2a)')        '! ENERGY'
    end if

    ! Never use symmetry
    write(iurun, '(a)')           '! Angs NoUseSym'

    ! Point charges for QM/MM
    if ( nclatoms > 0 ) then
       write(iurun, '(3a)')        '%pointcharges "',chgfile,'"'
    end if

    ! Charge / Spin / Coordinates
    write(iurun, '(a,i0,x,i0,x,a)') '*xyzfile ',charge,&
         spinmult, crdfile

    close(iurun)

    ! Now write xyz file with geometry
    open(iurun, file=crdfile, iostat=ios)
    if ( ios > 0 ) then
      call sander_bomb('write_inpfile (qm2_extern_orc_module)', &
           'Error opening coordinate file '//crdfile//'for writing', &
           'Will quit now')
    end if
    write(iurun,'(i5,/)') nqmatoms
    do i = 1, nqmatoms
       write(iurun,'(a2,1x,3f21.16)') &
          elementSymbol(qmtypes(i)), &
          qmcoords(:,i)
    enddo
    close(iurun)

    if ( nclatoms > 0 ) then
      call write_chgfile( chgfile, nclatoms, clcoords, orc_nml%debug )
    end if

    call debug_exit_function( 'write_inpfile', module_name, orc_nml%debug )

  end subroutine write_inpfile

  ! ---------------------
  ! Read Orca results
  ! ---------------------
  subroutine read_results( datfile, basename, nqmatoms, nclatoms, &
    escf, dxyzqm, dxyzcl, dipxyz, dipole, do_grad, do_dipole, debug )

    use constants, only:  CODATA08_AU_TO_DEBYE

    implicit none

    character(len=*), intent(in) :: datfile
    ! Used to get the .engrad and .pcgrad files
    character(len=*), intent(in) :: basename
    integer, intent(in) :: nqmatoms, nclatoms 
    _REAL_, intent(out) :: escf, dxyzqm(3,nqmatoms), dxyzcl(3,nclatoms)
    logical, intent(in) :: do_grad, do_dipole
    integer, intent(in) :: debug

    integer :: ios, i
    integer, parameter :: iunit = 351
    _REAL_, intent(out) :: dipxyz(3), dipole ! return in Debye
    character(len=256) :: read_buffer

    call debug_enter_function( 'read_results', module_name, debug )

    if ( do_grad ) then

      ! Energy and gradient; file will be empty if not do_grad
      open(iunit, file=basename//'.engrad', iostat=ios)
      if ( ios /= 0 ) then
        call sander_bomb('read_results (qm2_extern_orc_module)', &
          'Error opening Orca output file '//basename//'.engrad for reading', &
          'Will quit now')
      end if
      do
        read (iunit, '(a)', iostat = ios) read_buffer
        ! End of file; data not found
        if (ios < 0) then
          call sander_bomb('read_results (qm2_extern_orc_module)', &
            'Error reading Orca output from file '//basename//'.engrad', &
            '(Current total energy or gradient not found or unsupported units.)')
        end if
        if ( index(read_buffer,'# The current total energy in Eh') > 0 ) then
          ! Skip a line, read energy
          read (iunit, '(a)') read_buffer
          read (iunit,'(f22.12)') escf
        else if ( index(read_buffer,'# The current gradient in Eh/bohr') > 0 ) then
          ! Skip a line, read QM gradient data
          read (iunit, '(a)') read_buffer
          do i = 1, nqmatoms
            read (iunit, '(f16.10)') dxyzqm(1,i)
            read (iunit, '(f16.10)') dxyzqm(2,i)
            read (iunit, '(f16.10)') dxyzqm(3,i)
          end do
          exit
        end if
      end do
      close(iunit)

      ! MM gradients (none if do_grad)
      if ( nclatoms > 0 ) then
        open(iunit, file=basename//'.pcgrad', iostat=ios)
        if ( ios /= 0 ) then
          call sander_bomb('read_results (qm2_extern_orc_module)', &
            'Error opening Orca output file '//basename//'.pcgrad for reading', &
            'Will quit now')
        end if
        do
          read (iunit, '(a)', iostat = ios) read_buffer !First line is nclatoms
          ! End of file; data not found
          if (ios < 0) then
            call sander_bomb('read_results (qm2_extern_orc_module)', &
              'Error reading Orca output from file '//basename//'.pcgrad', &
              '(Empty file)')
          end if
          ! Read MM gradient data
          do i = 1, nclatoms
            read (iunit, '(3(f17.12))') dxyzcl(:,i)
          end do
          exit
        end do
        close(iunit)
      end if
    end if

    if ( do_dipole .or. (.not. do_grad) ) then
      ! Dipole moment
      open(iunit, file=datfile, iostat=ios)
      if ( ios /= 0 ) then
        call sander_bomb('read_results (qm2_extern_orc_module)', &
          'Error opening Orca output file '//datfile//' for reading', &
          'Will quit now')
      end if
      do
        read (iunit, '(a)', iostat = ios) read_buffer
        ! End of file; data not found
        if (ios < 0) then
          call sander_bomb('read_results (qm2_extern_orc_module)', &
            'Error reading Orca output from file '//datfile, &
            '(Dipole moment not found.)')
        end if
        if ( .not. do_grad .and. index(read_buffer, 'FINAL SINGLE POINT')>0 ) then
          read(read_buffer(29:),'(f22.12)') escf
        end if
        if (index(read_buffer,'Total Dipole Moment') > 0) then
          read(read_buffer(29:),*) dipxyz(:)
          read (iunit, '(a)') read_buffer ! Skip
          read (iunit, '(a)') read_buffer
          read(read_buffer(29:),*) dipole
          ! convert au to Debye
          dipxyz(:) = dipxyz(:) * CODATA08_AU_TO_DEBYE
          dipole = dipole * CODATA08_AU_TO_DEBYE
          exit
        end if
      end do
      close(iunit)
    end if

    call debug_exit_function( 'read_results', module_name, debug )

  end subroutine read_results

end module qm2_extern_orc_module
