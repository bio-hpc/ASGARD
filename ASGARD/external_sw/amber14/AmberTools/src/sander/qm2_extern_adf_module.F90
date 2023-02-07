#include "../include/dprec.fh"
module qm2_extern_adf_module
! ----------------------------------------------------------------
! Interface for ADF based QM MD 
!
! Currently supports:
! pure QM
!
! Initial implementation by
! Matthew Clark and Prithvi Undavalli
! (SDSC LSSI summer highschool students)
! under supervision of
! Andreas Goetz and Ross Walker (SDSC)
!
! Date: August 2010
!
! Extensions by Andreas Goetz (agoetz@sdsc.edu)
!
! Date: November 2010
!
! Input options extended by
! Prathyush Katukojwala, Evan Wildenhain and Pietro Sette
! (REHS 14 interns)
! under supervision of
! Andreas Goetz
! 
! Date: July 2014
! ----------------------------------------------------------------

  use qm2_extern_util_module, only: debug_enter_function, debug_exit_function

  implicit none

  private
  public :: get_adf_forces

  character(len=*), parameter :: module_name = "qm2_extern_adf_module"

  type adf_nml_type
     character(len=20) :: xc
     character(len=20) :: basis
     character(len=20) :: core
     character(len=20) :: fit_type
     character(len=20) :: guess
     _REAL_ :: integration
     _REAL_ :: scf_conv
     integer :: scf_iter
     integer :: ntpr
     integer :: num_threads
     integer :: linear_scaling
     integer :: debug
     logical :: use_dftb
     logical :: oldgradients
     logical :: dipole
     logical :: exactdensity
     logical :: use_template
     character(len=20) :: beckegrid
     character(len=20) :: zlmfit
     ! Deprecated
     integer :: charge
     integer :: spin
  end type adf_nml_type

contains

  ! -------------------------------------------------------
  ! Get QM energy and forces from ADF (adf or dftb program)
  ! -------------------------------------------------------
  subroutine get_adf_forces( do_grad, nstep, ntpr_default, id, natoms, &
       coords, qmtypes, escf, dxyzqm, charge, spinmult )

    use qm2_extern_util_module, only: print_results, check_installation, write_dipole
    use constants, only: CODATA08_AU_TO_KCAL, CODATA08_A_TO_BOHRS, CODATA08_AU_TO_DEBYE, ZERO
    use file_io_dat

    implicit none

    logical, intent(in)  :: do_grad            ! Return gradient/not
    integer, intent(in)  :: nstep              ! MD step number
    integer, intent(in)  :: ntpr_default       ! frequency of printing
    character(len=3), intent(in) :: id         ! ID number for PIMD or REMD
    integer, intent(in)  :: natoms             ! Total number of atoms
    _REAL_,  intent(in)  :: coords(3,natoms)   ! QM atom coordinates
    integer, intent(in)  :: qmtypes(natoms)    ! QM atom types (nuclear charge in au)
    _REAL_,  intent(out) :: escf               ! SCF energy
    _REAL_,  intent(out) :: dxyzqm(3,natoms)   ! SCF QM force
    integer, intent(in)  :: charge, spinmult   ! Charge and spin multiplicity

    _REAL_               :: dipxyz(3), diptot  ! Dipole Moment

    type(adf_nml_type), save :: adf_nml
    logical, save :: first_call = .true.
    integer :: i
    integer :: printed =-1 ! Used to tell if we have printed this step yet 
                           ! since the same step may be called multiple times

    character(len=150)           :: call_buffer
    character(len=80), save      :: program
    character(len=*),  parameter :: basename = 'adf_job'
    character(len=*),  parameter :: inpext = '.inp'
    character(len=*),  parameter :: outext = '.out'
    character(len=*),  parameter :: dipext = '.dip'
    character(len=*),  parameter :: tplext = '.tpl'
    character(len=14)            :: inpfile, outfile, dipfile, keyfile, tplfile
    ! Need to prepend subdirectory if doing REMD, PIMD or multi-region QM/MM. 
    !   This is triggered if 'id' is defined (not empty). 
    character(len=25)            :: subdir 

    ! assemble input - / output data filenames
    inpfile = basename//trim(id)//inpext 
    outfile = basename//trim(id)//outext
    dipfile = basename//trim(id)//dipext
    tplfile = basename//tplext

    ! Setup on first call
    if ( first_call ) then
      first_call = .false.
      write (6,'(/,a,/)') '  >>> Running QM calculation with ADF <<<'
      call get_namelist(ntpr_default, adf_nml)
      call print_namelist( adf_nml ) 
      if ( adf_nml%use_dftb ) then
        program = 'dftb'
      else
        program = 'adf'
      end if
      call check_installation( program, id, .true., adf_nml%debug )

      write (6,'(80a)') ('-', i=1,80)
      write (6,'(a)') '   4.  RESULTS'
      write (6,'(80a)') ('-', i=1,80)
    end if
    if ( adf_nml%use_dftb ) then
      keyfile='DFTB.kf'
    else
      keyfile='TAPE21'
    end if
    ! Remove the logfile at the beginning of a run so only the latest run is
    ! stored. ADF will also throw a false 'error detected' if it finds a TAPE13 
    ! file in the working directory; this is to prevent confusion.
    call system('rm -f logfile TAPE13')

    call system('rm -f '//inpfile)
    call write_inpfile( trim(inpfile), trim(tplfile), coords, qmtypes, natoms, &
      adf_nml, do_grad, charge, spinmult )

    if ( adf_nml%debug > 1 ) then
      write(6,*) ' ADF input file successfully written; calling ADF...'
    end if

    ! Run ADF
    ! Separate runs into different directories if we are doing PIMD
    subdir=''
    call_buffer=''
    if ( trim(id)/='' ) then
      subdir='./'//trim(id)//'/'
      call_buffer=' mkdir -p '//trim(subdir)//'; cd '//trim(subdir)//'; mv ../'//inpfile//' .;'
    end if

    write(call_buffer,'(a)') trim(call_buffer)//trim(program)
    ! Inputting 0 will cause ADF to use default number of threads
    ! Otherwise, we add what the user specified
    if ( adf_nml%num_threads /= 0 ) then
      write(call_buffer,'(a,i0)') trim(call_buffer)//' -n ',adf_nml%num_threads
    end if

    write(call_buffer,'(a)') trim(call_buffer)//' < ./'//trim(inpfile)//' > '//outfile 

    call system(trim(call_buffer))
 
    if ( adf_nml%debug > 1 ) then
      write(6,*) ' ADF execution success; Processing ADF results...'
    end if

    ! Call read_adf_results - a function in C program (qm2_read_adf_results.c) 
    ! to read TAPE21 file; will output the data to escf and dxyzqm
    call read_adf_results( natoms, escf, dxyzqm, adf_nml%use_dftb, dipxyz, &
      trim(subdir)//trim(keyfile)//char(0), do_grad )

    ! Call write_dipole - dipfile, dipxyz, and magnitude of dipxyz;
    ! will write output to dipfile
    if ( adf_nml%ntpr > 0 .and. mod(nstep, adf_nml%ntpr) == 0 ) then
      if ( printed /= nstep .and. adf_nml%dipole ) then
         dipxyz = dipxyz * CODATA08_AU_TO_DEBYE
         diptot = dsqrt( dipxyz(1)**2 + dipxyz(2)**2 + dipxyz(3)**2 )
        call write_dipole(dipfile, dipxyz, diptot, adf_nml%debug)
        printed = nstep
      end if
    end if

    if ( adf_nml%use_dftb ) then
      call system('mv '//trim(subdir)//keyfile//' '//trim(subdir)//'kf.DFTB')
    else
      ! Save old TAPE21 file to restart from in next calculation
      call system('mv '//trim(subdir)//keyfile//' '//trim(subdir)//'adf.t21')  
    end if
    call system('mv '//trim(subdir)//trim(inpfile)//' '//trim(subdir)//'old.'//inpfile)
    call system('mv '//trim(subdir)//trim(outfile)//' '//trim(subdir)//'old.'//outfile)

    ! Convert gradient from au to kcal/(mol*A)
    if ( do_grad ) then
      dxyzqm(:,:) = dxyzqm(:,:) * CODATA08_AU_TO_KCAL * CODATA08_A_TO_BOHRS
    else
      dxyzqm = ZERO
    end if
    
    escf = escf * CODATA08_AU_TO_KCAL

    call print_results( 'qm2_extern_adf_module', escf, natoms, dxyzqm, &
      adf_nml%debug )

  end subroutine get_adf_forces

  ! ----------------------------------------
  ! Read ADF namelist values from file mdin,
  ! use default values if none are present.
  ! ----------------------------------------
  subroutine get_namelist(ntpr_default, adf_nml)
    
    implicit none
    integer, intent(in) :: ntpr_default
    type(adf_nml_type), intent(out) :: adf_nml

    character(len=20) :: xc, basis, core, fit_type, guess, beckegrid, zlmfit, &
                         beckegrid_default, zlmfit_default
    integer :: scf_iter, ntpr, num_threads, linear_scaling, debug, use_dftb, &
    oldgradients, dipole, exactdensity, use_template
    integer :: charge, spin ! deprecated
    _REAL_ :: integration, scf_conv, scf_conv_dftb, scf_conv_adf
    namelist /adf/ xc, basis, core, fit_type, guess, &
         integration, scf_conv, scf_iter, ntpr, &
         num_threads, linear_scaling, debug, use_dftb, oldgradients, &
         dipole, exactdensity, use_template,&
         beckegrid, zlmfit, charge, spin 

    integer :: ierr

    ! Set default values for adf namelist values
    xc           = 'GGA BLYP'
    basis        = 'DZP'
    core         = 'None'
    fit_type     = ''
    guess        = 'read'
    integration  = -1.0D0
    scf_conv     = 1.0D30    ! dummy value
    scf_conv_dftb= 1.0d-12   ! dftb default value
    scf_conv_adf = 1.0d-06   ! adf default value
    ntpr         = ntpr_default
    scf_iter     = 50
    num_threads  = 0         ! 0 = ADF default (all available)
    linear_scaling =-1
    debug        = 0
    use_dftb     = 0
    oldgradients = 0
    dipole       = 0
    exactdensity = 0
    use_template = 0
    beckegrid    = ''
    beckegrid_default = 'good'
    zlmfit       = ''
    zlmfit_default = 'good'

    ! These are now deprecated and should be specified in the &qmmmm namelist
    charge = -351
    spin   = -351
    
    ! Read namelist
    rewind 5
    read(5,nml=adf,iostat=ierr)

    if ( ierr > 0 ) then
      call sander_bomb('get_namelist (qm2_extern_adf_module)', &
        '&adf namelist read error', &
        'Please check your input.')
    else if ( ierr < 0 ) then
       write(6,'(a,/,a)') '&adf namelist read encountered end of file', &
         'Please check your input if the calculation encounters a problem'
    end if

    if ( charge /= -351 .or. spin /= -351 ) then
      call sander_bomb('get_namelist (qm2_extern_adf_module)', &
        'The charge and spin keywords are deprecated', &
        'Please specify charge (qmcharge) and spin multiplicity (spin) in the &qmmm namelist.')
    end if

    ! Assign namelist values to adf_nml data type
    adf_nml%xc           = xc
    adf_nml%basis        = basis
    adf_nml%core         = core
    adf_nml%guess        = guess

    adf_nml%fit_type = fit_type
    adf_nml%zlmfit   = zlmfit
    if ( (fit_type /= '') .and. (zlmfit /= '') ) then
       call sander_bomb('get_namelist (qm2_extern_adf_module)',&
            '&adf: you cannot use pair fit and Zlm fit at the same time.',&
            'Please check your input (fit_type, zlmfit).')
    else if ( (fit_type == '') .and. (zlmfit == '') ) then
       ! use zlmfit by default
       adf_nml%zlmfit = zlmfit_default
    end if

    adf_nml%integration = integration
    adf_nml%beckegrid   = beckegrid
    if ( (integration > 0.0d0) .and. (beckegrid /= '') ) then
       call sander_bomb('get_namelist (qm2_extern_adf_module)',&
            '&adf: you cannot use the teVelde Baerends and Becke grid at the same time.',&
            'Please check your input (integration, beckegrid).')
    else if ( (integration < 0.0d0) .and. (beckegrid == '') ) then
       ! use beckegrid by default
       adf_nml%beckegrid = beckegrid_default
    end if

    if ( use_dftb > 0 ) then
       adf_nml%use_dftb = .true.
    else
       adf_nml%use_dftb = .false.
    end if
    if ( scf_conv > 1.0D29 ) then
       adf_nml%scf_conv = scf_conv_dftb ! use default value
    else
       adf_nml%scf_conv = scf_conv
    end if
    if ( scf_conv > 1.0D29 ) then
       adf_nml%scf_conv = scf_conv_adf  ! use default value
    else
       adf_nml%scf_conv = scf_conv
    end if
    adf_nml%scf_iter     = scf_iter
    adf_nml%ntpr         = ntpr
    adf_nml%num_threads  = num_threads
    adf_nml%linear_scaling = linear_scaling
    adf_nml%debug        = debug
    if ( oldgradients == 1 ) then
       adf_nml%oldgradients = .true.
    else if ( oldgradients == 0 ) then
       adf_nml%oldgradients = .false.
    else
      call sander_bomb('get_namelist (qm2_extern_adf_module)', &
        '&adf oldgradients value not allowed', &
        'Please check your input. oldgradients can only be 0 or 1.')
    end if

    if ( dipole == 1 ) then
       adf_nml%dipole = .true.
    else if ( dipole == 0 ) then
       adf_nml%dipole = .false.
    else
      call sander_bomb('get_namelist (qm2_extern_adf_module)', &
        '&adf dipole value not allowed', &
        'Please check your input. dipole can only be 0 or 1.')
    end if

    if ( exactdensity == 1 ) then
       adf_nml%exactdensity = .true.
    else if ( exactdensity == 0 ) then
       adf_nml%exactdensity = .false.
    else
      call sander_bomb('get_namelist (qm2_extern_adf_module)', &
        '&adf exactdensity value not allowed', &
        'Please check your input. exactdensity can only be 0 or 1.')
    end if

    if ( use_template == 0 ) then
       adf_nml%use_template = .false.
    else if ( use_template == 1 ) then
       adf_nml%use_template = .true.
    else
      call sander_bomb('get_namelist (qm2_extern_adf_module)', &
        '&adf use_template value not allowed', &
        'Please check your input. use_template can only be 0 or 1.')
    end if

  end subroutine get_namelist
  
  ! ---------------------------
  ! Print ADF namelist settings
  ! ---------------------------
  subroutine print_namelist( adf_nml )

    implicit none
    type(adf_nml_type), intent(in) :: adf_nml

    write(6, '(a)')          '| &adf'
    if ( .not. adf_nml%use_dftb ) then
      write(6, '(2a)')       '|   xc             = ', adf_nml%xc
      write(6, '(2a)')       '|   basis          = ', adf_nml%basis
      write(6, '(2a)')       '|   core           = ', adf_nml%core
      write(6, '(2a)')       '|   fit_type       = ', adf_nml%fit_type
      write(6, '(2a)')       '|   guess          = ', adf_nml%guess
      write(6, '(a,E22.16)') '|   integration    = ', adf_nml%integration
    end if
    write(6, '(a,es10.2)')   '|   scf_conv       = ', adf_nml%scf_conv
    write(6, '(a,i5)')       '|   scf_iter       = ', adf_nml%scf_iter
    write(6, '(a,i0)')       '|   ntpr           = ', adf_nml%ntpr
    write(6, '(a,i3)')       '|   num_threads    = ', adf_nml%num_threads
    write(6, '(a,i3)')       '|   linear_scaling = ', adf_nml%linear_scaling
    write(6, '(a,l)')        '|   use_dftb       = ', adf_nml%use_dftb
    if ( .not. adf_nml%use_dftb ) then
      write(6, '(a,l)')      '|   oldgradients   = ', adf_nml%oldgradients
      write(6, '(a,l)')      '|   dipole         = ', adf_nml%dipole
      write(6, '(a,l)')      '|   exactdensity   = ', adf_nml%exactdensity
    end if
    write(6, '(a,l)')        '|   use_template   = ', adf_nml%use_template
    write(6, '(2a)')        '|   beckegrid   = ', adf_nml%beckegrid
    write(6, '(2a)')        '|   zlmfit   = ', adf_nml%zlmfit
    write(6,'(a)')           '| /'

  end subroutine print_namelist

  ! ------------------------
  ! Write input file for ADF
  ! ------------------------
  subroutine write_inpfile( inpfile, tplfile, coords, qmtypes, natoms, &
       adf_nml, do_grad, charge, spinmult )

    use ElementOrbitalIndex, only : elementSymbol

    implicit none

    character(len=*), intent(in) :: inpfile, tplfile
    _REAL_,  intent(in) :: coords(:,:)
    integer, intent(in) :: qmtypes(:)
    integer, intent(in) :: natoms
    type(adf_nml_type), intent(in) :: adf_nml
    logical, intent(in) :: do_grad
    integer, intent(in) :: charge, spinmult

    character(len=20) :: bakfile

    integer :: temp_atoms(natoms), i, j, ierr, tempvar
    integer, parameter :: iunit = 351
    character(len=20) :: tmp_buffer=''
    logical, save :: first_call = .true.

    ! 'temp_atoms' is used to write only the unique atoms in the fragments keyword
    temp_atoms = qmtypes

    call debug_enter_function( 'write_inpfile', module_name, adf_nml%debug )

    if ( adf_nml%use_template ) then
      bakfile = tplfile//'.bak'
      if ( first_call ) then
        call copy_template( tplfile, trim(bakfile), do_grad, adf_nml%debug )
        call system('cp '//tplfile//' '//inpfile)
      else
        call system('cp '//trim(bakfile)//' '//inpfile)
      end if
    end if

    open(iunit, file=inpfile, iostat=ierr, position='append')
    if ( ierr /= 0 ) then
      call sander_bomb('write_inpfile (qm2_extern_adf_module)', &
        'Error opening ADF inpfile '//inpfile//' for writing', &
        'Will quit now.')
    end if

    ! ATOMS keyword
    write (iunit,'(a)') 'Atoms'
    do i = 1, natoms
      write(iunit,'(a2,1x,3f25.16)') elementSymbol(qmtypes(i)), coords(1:3,i)
    end do
    write(iunit,'(a,/)') 'End'
    
    ! Don't include BASIS/FRAGMENTS etc if we are on the first call
    if ( first_call .and.  adf_nml%use_template ) then
      first_call=.false.
      if ( do_grad ) write (iunit,'(a)') 'GRADIENT'
      return
    end if

    if ( .not. adf_nml%use_dftb ) then

      ! ---------------------------
      ! ADF PROGRAM INPUT SPECIFICS
      ! ---------------------------

      ! BASIS/RESTART/FRAGMENTS keyword
      if ( first_call ) then

        first_call = .false.
        
        ! write this only in the first call, not needed for restarts
        write (iunit,'(a,/,2(a,a,/))') &
          'Basis'                      , &      
          ' type ',trim(adf_nml%basis) , &
          ' core ',trim(adf_nml%core)
        if ( adf_nml%fit_type/='standard' .and. adf_nml%fit_type /= '') then
           tmp_buffer=' FitType '// trim(adf_nml%fit_type)
           write(iunit, '(a,/)') tmp_buffer
        end if
        write(iunit, '(2(a,/))') &
             ' createoutput None', &
             'End'
     else
        ! Use last t21 file as restart file
        if ( trim(adf_nml%guess) == 'read' ) then
           write(iunit,'(a,/,a,/,a,/)') 'Restart adf.t21 &', 'nogeo', 'END'
        end if
        ! Set any duplicate atomic numbers in our array to zero
        do i = 1, natoms
          tempvar = qmtypes(i)
          do j = i+1, natoms
            if ( temp_atoms(j) == tempvar ) then
              temp_atoms(j)=0
            end if
          end do
        end do
        ! Begin writing fragments
        write (iunit,'(a)') 'Fragments'
        do i = 1, natoms
          if ( temp_atoms(i) /= 0 ) then
            write (iunit,'(a,a,a)') elementSymbol(temp_atoms(i)), ' t21.', &
              elementSymbol(temp_atoms(i))
          end if
        end do
        write (iunit,'(a,/)') 'End'
      end if

      if ( adf_nml%use_template ) then
        return
      end if

      ! ELECTRIC field (for external point charges for QM/MM)
      ! AWG: This is disabled at the moment
      ! AWG: We need to add QM region extraction and link atom setup first
      ! if ( natoms > num_atoms ) then
      !   write(iunit,'(a)') 'Efield'
      !   do i = 1, natoms
      !     if ( qmmm_struct%qm_xcrd(4,i) /= 0 ) then
      !       write(iunit, '(f0.12,1x, f0.12, 1x, f0.12, 1x, f0.12)') qmmm_struct%qm_xcrd(1:4,i)
      !     end if
      !   end do
      !     write(iunit,'(a,/)') 'End' 
      ! end if
       
      ! XC keyword
      write(iunit,'(a,/,a,/,a,/)')&
           'XC'             , &
           trim(adf_nml%xc) , &
           'END'
      
      ! SCF keyword
      write (iunit,'(a,/,a,i0,/,a,E22.16,/,a,/)') &
           'SCF '                         , &
           'iterations ', adf_nml%scf_iter, &
           'converge '  , adf_nml%scf_conv, &
           'END'

      ! INTEGRATION, BECKEGRID, ZLMFIT keywords
      if(adf_nml%integration > 0.0d0) then 
         write (iunit,'(a,E22.16,/)') 'INTEGRATION ', adf_nml%integration
      end if

      if ( adf_nml%beckegrid /= '' ) then
         write (iunit,'(a,/,2a,/,a,/)') &
           'BECKEGRID'                         , &
           ' quality ',trim(adf_nml%beckegrid) , &
           'END'
      end if

      if ( adf_nml%zlmfit /= '' ) then
         write (iunit,'(a,/,2a,/,a,/)') &
           'ZLMFIT'                         , &
           ' quality ',trim(adf_nml%zlmfit) , &
           'END'
      end if

      ! GRADIENT, SYMMETRY, EXACTDENSITY keywords!
      if ( adf_nml%oldgradients ) then
        write(iunit,'(a)')'OLDGRADIENTS'
      end if
      ! Need to retreive force data
      if ( do_grad ) write (iunit,'(a)') 'GRADIENT'
      ! Run without symmetry, we won't need it for MD
      write (iunit,'(a)') 'SYMMETRY NOSYM'
      ! Use exact density for XC calculations if specified
      if ( adf_nml%exactdensity ) then
        write (iunit,'(a)') 'EXACTDENSITY'
      end if

      ! LINEARSCALING Keyword
      if ( adf_nml%linear_scaling/=-1 ) then
        write (iunit,'(a,i3)') 'LINEARSCALING ',adf_nml%linear_scaling
      end if

      ! CHARGE keyword
      ! Note ADF wants the number of unpaired electrons (spinmult-1)
      write(iunit,'(a,i2,i2,/)') 'CHARGE ', charge, spinmult-1
      if ( spinmult > 1 ) then
        write(iunit,'(a,/)') 'UNRESTRICTED'
      end if

      ! SAVE files (Need TAPE21 to extract sander data and restart)
      write(iunit,'(a,/)')'SAVE TAPE21'

    else   

      ! ----------------------------
      ! DFTB PROGRAM INPUT SPECIFICS
      ! ----------------------------

      write(iunit,'(a,i2,/)') 'CHARGE ', charge

      ! SCF Keywords
      write (iunit,'(a,/,a,E22.16,/,a,/)') &
           'SCF '                       , &
           'converge ', adf_nml%scf_conv, &
           'END'

    end if

    ! End writing inpfile for ADF
    close(iunit, iostat=ierr)
    if ( ierr /= 0 ) then
      call sander_bomb('write_inpfile (qm2_extern_adf_module)', &
          'Error closing ADF inpfile for writing', &
          'Will quit now.')
    end if

    call debug_exit_function( 'write_inpfile', module_name, adf_nml%debug )

  end subroutine write_inpfile
  

  subroutine copy_template( tplfile, bakfile, do_grad, debug )

    use UtilitiesModule, only: Upcase

    implicit none
    character(len=*), intent(in) :: tplfile, bakfile
    logical, intent(in) :: do_grad
    integer, intent(in) :: debug

    integer, parameter :: tplunit = 351, bakunit=352
    character(len=256) :: read_buffer
    integer :: tplerr, bakerr, ios
    logical :: in_basis = .false.

    call debug_enter_function( 'copy_template', module_name, debug )

    open(tplunit, file=tplfile, iostat=tplerr )
    open(bakunit, file=bakfile, iostat=bakerr )
    if ( tplerr /= 0 ) then
      call sander_bomb('copy_template (qm2_extern_adf_module)', &
        'Error opening ADF template file '//tplfile//' for reading', &
        'Will quit now.')
    end if
    if ( bakerr /= 0 ) then
      call sander_bomb('copy_template (qm2_extern_adf_module)', &
        'Error opening ADF template backup file '//bakfile//' for writing', &
        'Will quit now.')
    end if

    ! Write tplfile (without basis key) to bakfile
    do
       read (tplunit, '(a)', iostat = ios) read_buffer
       ! End of file; stop writing
       if ( ios < 0 ) then
         exit
       end if
       if ( index(Upcase(read_buffer), 'BASIS') > 0 ) then
          ! Stop writing until past basis
          in_basis=.true.
       end if
       if ( .not. in_basis ) then
         write (bakunit, '(a)') read_buffer
       else if ( in_basis .and. index(Upcase(read_buffer), 'END') > 0 ) then
         in_basis=.false.
       end if
    end do

    if ( do_grad ) write(bakunit,'(a)') 'GRADIENT'
      
    close(tplunit)
    close(bakunit)

    call debug_exit_function( 'copy_template', module_name, debug )

  end subroutine copy_template

end module qm2_extern_adf_module
