#include "../include/dprec.fh"
module qm2_extern_gau_module
! ----------------------------------------------------------------
! Interface for Gaussian based QM MD 
!
! Currently supports:
! pure QM
!
! Initial implementation by
! Matthew Clark
! under supervision of
! Andreas Goetz and Ross Walker (SDSC)
! 
! Date: February 2011
!
! Extensions by Andreas Goetz (SDSC)
!
! ----------------------------------------------------------------

  use qm2_extern_util_module, only: debug_enter_function, debug_exit_function

  implicit none

  private
  public :: get_gau_forces

  character(len=*), parameter :: module_name = "qm2_extern_gau_module"

  type gau_nml_type
     character(len=20) :: method
     character(len=20) :: basis
     character(len=50) :: mem
     character(len=20) :: guess
     integer :: scf_conv
     integer :: ntpr
     integer :: num_threads 
     integer :: debug
     logical :: dipole
     logical :: use_template
     
     ! Deprecated
     integer :: charge
     integer :: spinmult
  end type gau_nml_type

contains

  ! --------------------------------------------
  ! Get QM energy and forces from Gaussian
  ! --------------------------------------------
  subroutine get_gau_forces( do_grad, nstep, ntpr_default, id, &
       nqmatoms, qmcoords, qmtypes, nclatoms, clcoords, &
       escf, dxyzqm, dxyzcl, charge, spinmult )

    use qm2_extern_util_module, only: print_results, check_installation, write_dipole
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

    type(gau_nml_type), save     :: gau_nml
    logical, save                :: first_call = .true.
    logical                      :: exist
    integer                      :: i
    integer                      :: printed =-1 ! Used to tell if we have printed this step yet 
                                                ! since the same step may be called multiple times
    character(len=1024)          :: call_buffer
    character(len=80), save      :: program     = 'g09'
    character(len=*), parameter  :: alt_program = 'g03'
    character(len=*), parameter  :: basename = 'gau_job'
    character(len=*), parameter  :: inpext = '.inp'
    character(len=*), parameter  :: logext = '.log'
    character(len=*), parameter  :: rstext = '.chk' ! Restart from checkpoint files
    character(len=*), parameter  :: frstext= '.fchk'! Formatted checkpoint file
    character(len=*), parameter  :: dipext = '.dip'
    character(len=*), parameter  :: tplext = '.tpl'
    character(len=14)            :: inpfile, rstfile, frstfile, logfile, dipfile, tplfile
    ! Need to prepend subdirectory if doing REMD, PIMD or multi-region QM/MM. 
    !   This is triggered if 'id' is defined (not empty). 
    character(len=25)            :: subdir 

    ! for system call
    integer :: system
    integer :: stat

    ! assemble input - / output data filenames
    inpfile = basename//trim(id)//inpext
    logfile = basename//trim(id)//logext
    rstfile = basename//trim(id)//rstext
    frstfile= basename//trim(id)//frstext
    dipfile = basename//trim(id)//dipext
    tplfile = basename//tplext

    ! Setup on first call
    if ( first_call ) then
      first_call = .false.
      write (6,'(/,a,/)') '  >>> Running calculations with Gaussian <<<'
      call get_namelist( ntpr_default, gau_nml )
      call print_namelist( gau_nml ) 

      ! Check for version of Gaussian to use; store as 'program'
      call check_installation( trim(program), id, .false., gau_nml%debug, found=exist )
      ! If we did not find g09, try looking for g03 instead
      if ( exist .eqv. .false. ) then
        program = alt_program
        ! Will quit if this is not found
        call check_installation( trim(program), id, .true., gau_nml%debug )
      end if

      write (6,'(80a)') ('-', i=1,80)
      write (6,'(a)') '   4.  RESULTS'
      write (6,'(80a)') ('-', i=1,80)
      ! Remove old inpfile, logfile and dipfile at the 
      ! beginning of a run so only the latest run is stored.
      stat = system('rm -f '//inpfile//' '//logfile//' '//rstfile//' '//dipfile)
      if ( stat /= 0 ) then
        call sander_bomb('get_gau_forces (qm2_extern_gau_module)', & 
          'Error with system call (removing files)', &
          'Will quit now.')
      end if
    end if
    
    call write_inpfile( trim(inpfile), trim(rstfile), trim(tplfile), &
         nqmatoms, qmcoords, qmtypes, nclatoms, clcoords, &
         gau_nml, do_grad, charge, spinmult)

    if ( gau_nml%debug > 0 ) then
      write(6,'(a)') ' Input file written successfully; calling Gaussian...'
    end if

    ! Run g09/g03
    ! Separate runs into different directories if we are doing PIMD or REMD
    subdir=''
    call_buffer=''
    if (trim(id)/='') then 
      subdir='./'//trim(id)//'/'
      call_buffer=' mkdir -p '//trim(subdir)//'; cd '//trim(subdir)//'; mv ../'//inpfile//' .;'
    end if
    call_buffer = trim(call_buffer)//' '//trim(program)//' '//inpfile
    call_buffer = trim(call_buffer)//'; formchk '//trim(rstfile)//' '//trim(frstfile)//' >/dev/null'
    stat = system(call_buffer) 
    if ( stat /= 0 ) then
      call sander_bomb('get_gau_forces (qm2_extern_gau_module)', & 
        'Error with system call (executing Gaussian)', &
        'Will quit now.')
    end if

    if ( gau_nml%debug > 0 ) then    
      write(6,'(a)') ' Gaussian execution success; Processing Gaussian results...'
    end if    

    ! Call read_results - retrieve data from Gaussian .log file
    ! Will output data to escf and dxyqm for pure QM or MM runs
    ! For QM/MM runs will also return dxyzcl containing electric field strength
    
    ! Search in subdir for logfile if doing PIMD
    ! Otherwise, search current directory
    call read_results( trim(subdir)//trim(logfile), trim(subdir)//trim(frstfile), &
         nqmatoms, escf, dxyzqm, nclatoms, dxyzcl, dipxyz, dipole, &
         do_grad, gau_nml%debug )

    ! Call write_dipole with dipfile, dipxyz, and magnitude of dipxyz;
    ! will write output to dipfile
    if ( gau_nml%ntpr > 0 .and. mod(nstep, gau_nml%ntpr) == 0 ) then
      if ( printed /= nstep .and. gau_nml%dipole ) then
        call write_dipole( trim(dipfile), dipxyz, dipole, gau_nml%debug )
        printed = nstep
      end if
    end if

    ! Save copy of last input and log files
    stat = system('mv '//trim(subdir)//inpfile//' '//trim(subdir)//'old.'//inpfile)
    stat = stat + system('mv '//trim(subdir)//logfile//' '//trim(subdir)//'old.'//logfile)
    if ( stat /= 0 ) then
      call sander_bomb('get_gau_forces (qm2_extern_gau_module)', & 
        'Error with system call (moving / removing files)', &
        'Will quit now.')
    end if
    

    ! F = E*q to get gradients
    ! Note dxyzcl is currently holding the electric field strength
    do i = 1, nclatoms
      dxyzcl(:,i)= -dxyzcl(:,i) * clcoords(4,i)
    end do

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

    call print_results( 'qm2_extern_gau_module', escf, nqmatoms, dxyzqm, &
      gau_nml%debug, nclatoms, dxyzcl )

  end subroutine get_gau_forces

  ! ---------------------------------------------
  ! Read Gaussian namelist values from file mdin,
  ! use default values if none are present.
  ! ---------------------------------------------
    
 
  subroutine get_namelist(ntpr_default, gau_nml)

    implicit none
    integer, intent(in) :: ntpr_default
    type(gau_nml_type), intent(out) :: gau_nml

    character(len=20) :: method, basis, guess
    character(len=50) :: mem
    integer :: debug
    integer :: scf_conv, ntpr, num_threads, dipole, use_template
    integer :: charge, spinmult ! deprecated
    namelist /gau/ method, basis, mem, guess, scf_conv, ntpr, &
      num_threads, debug, dipole, use_template,&
      charge, spinmult

    integer :: ierr

    ! Set default values for gau namelist values
    method       = 'BLYP'
    basis        = '6-31G*'
    mem          = '256MB'
    guess        = 'read'
    scf_conv     = 8 
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
    read(5,nml=gau,iostat=ierr)

    if ( ierr > 0 ) then
       call sander_bomb('get_namelist (qm2_extern_gau_module)', &
            '&gau namelist read error', &
            'Please check your input.')
    else if ( ierr < 0 ) then
       write(6,'(a,/,a)') '&gau namelist read encountered end of file', &
            'Please check your input if the calculation encounters a problem'
    end if

    if ( charge /= -351 .or. spinmult /= -351 ) then
      call sander_bomb('get_namelist (qm2_extern_gau_module)', &
        'The charge and spin keywords are deprecated', &
        'Please specify charge (qmcharge) and spin multiplicity (spin) in the &qmmm namelist.')
    end if

    ! Assign namelist values to gau_nml data type
    gau_nml%method       = method
    gau_nml%basis        = basis
    gau_nml%mem          = mem
    gau_nml%guess        = guess
    gau_nml%scf_conv     = scf_conv
    gau_nml%ntpr         = ntpr
    gau_nml%num_threads  = num_threads
    gau_nml%debug        = debug
    if ( dipole == 0 ) then
       gau_nml%dipole = .false.
    else if ( dipole == 1 ) then
       gau_nml%dipole = .true.
    else
       call sander_bomb('get_namelist (qm2_extern_gau_module)', &
            '&gau dipole value not allowed', &
            'Please check your input. dipole can only be 0 or 1.')
    end if

    if ( use_template == 0 ) then
       gau_nml%use_template = .false.
    else if ( use_template == 1 ) then
       gau_nml%use_template = .true.
    else
       call sander_bomb('get_namelist (qm2_extern_gau_module)', &
            '&gau use_template value not allowed', &
            'Please check your input. use_template can only be 0 or 1.')
    end if

  end subroutine get_namelist

  ! --------------------------------
  ! Print Gaussian namelist settings
  ! --------------------------------
  subroutine print_namelist(gau_nml)

    implicit none
    type(gau_nml_type), intent(in) :: gau_nml

    write(6, '(a)')       '| &gau'
    write(6, '(2a)')      '|   method       = ', gau_nml%method
    write(6, '(2a)')      '|   basis        = ', gau_nml%basis
    write(6, '(2a)')      '|   guess        = ', gau_nml%guess
    write(6, '(2a)')      '|   mem          = ', trim(gau_nml%mem)
    write(6, '(a,i2)')    '|   scf_conv     = ', gau_nml%scf_conv
    write(6, '(a,i0)')    '|   ntpr         = ', gau_nml%ntpr
    write(6, '(a,i2)')    '|   num_threads  = ', gau_nml%num_threads
    write(6, '(a,i2)')    '|   debug        = ', gau_nml%debug
    write(6, '(a,l)')     '|   dipole       = ', gau_nml%dipole
    write(6, '(a,l)')     '|   use_template = ', gau_nml%use_template
    write(6,'(a)')        '| /'

  end subroutine print_namelist

  ! -----------------------------
  ! Write input file for Gaussian
  ! -----------------------------

  subroutine write_inpfile( inpfile, rstfile, tplfile, &
       nqmatoms, qmcoords, qmtypes, nclatoms, clcoords, &
       gau_nml, do_grad, charge, spinmult )

    use ElementOrbitalIndex, only : elementSymbol

    implicit none

    character(len=*), intent(in)   :: inpfile, rstfile, tplfile
    integer, intent(in)            :: nqmatoms
    _REAL_,  intent(in)            :: qmcoords(:,:)
    integer, intent(in)            :: qmtypes(:)
    integer, intent(in)            :: nclatoms
    _REAL_,  intent(in)            :: clcoords(:,:)
    type(gau_nml_type), intent(in) :: gau_nml
    logical, intent(in)            :: do_grad
    integer, intent(in)            :: charge, spinmult

    integer, parameter :: iunit = 351, tplunit = 352
    integer            :: i, ierr
    integer            :: tplerr, ios
    character(len=256) :: route
    character(len=256) :: read_buffer
    logical, save      :: first_call = .true.

    call debug_enter_function( 'write_inpfile', module_name, gau_nml%debug )

    open(iunit, file=inpfile, iostat=ierr)
    if ( ierr /= 0 ) then
      call sander_bomb('write_inpfile (qm2_extern_gau_module)', &
        'Error opening Gaussian inpfile '//inpfile//' for writing', &
        'Will quit now.')
    end if

    ! Write link option
    ! chkfile so we can restart from previous runs
    write (iunit,'(a)') '%chk='//rstfile

    ! Shared memory parallelism
    write (iunit, '(a,i0)') '%NProcShared=',gau_nml%num_threads

    ! RAM to be used
    write (iunit, '(2a)') '%mem=',gau_nml%mem

    ! Assemble route
    ! If using template, write route information in template file to route card
    route=''
    if ( gau_nml%use_template) then

      open(tplunit, file=tplfile, iostat=tplerr)
      if ( tplerr /= 0 ) then
        call sander_bomb('write_inpfile (qm2_extern_gau_module)', &
          'Error opening Gaussian template file '//tplfile//' for reading', &
          'Will quit now.')
      end if
      read (tplunit, '(a)', iostat = ios) read_buffer
      if (ios < 0) then
        call sander_bomb('write_inpfile (qm2_extern_gau_module)', &
          'Error reading Gaussian template file '//tplfile, &
          'Will quit now.')
      end if
      close(tplunit)
      write(route,'(a)') trim(route)//' '//trim(read_buffer)//' '

    else

      ! Method/Basis
      route = '#P '//trim(gau_nml%method)//'/'//trim(gau_nml%basis)
      ! SCF convergence setting
      write(route,'(a,i0,a)') trim(route)//' SCF=(Conver=', gau_nml%scf_conv, ')'
    end if
    write (route,'(a)') trim(route)//' NoSymm'
    ! Gradient or single point
    if ( do_grad ) then
      route = trim(route)//' Force'
    end if
    ! If we are not on our first run, restart from .chk file 
    if ( .not. first_call .and. (trim(gau_nml%guess) == 'read') ) then
      route = trim(route)//' Guess=Read'
    end if
    ! If doing electrostatic embedding QM/MM, 
    ! read external point charges, print to .log
    if ( nclatoms > 0 ) then
      route = trim(route)//' Charge Prop=(Field,Read)'
      if ( index(gau_nml%method, 'MP2') > 0 ) then
         route = trim(route)//' Density=MP2'
      end if
    end if

    write(iunit,'(a)') trim(route)
    write(iunit,'(a)')

    ! Write comment line
    if ( gau_nml%use_template) then
      write (iunit,'(a,/)') 'Gaussian run using SANDER external interface, route from template '//tplfile
    else
      write (iunit,'(a,/)') 'Gaussian run using SANDER external interface.'
    end if

    ! Write charge and spin
    write (iunit,'(i0,a,i0)') charge,' ', spinmult

    ! Write QM atoms and coordinates
    do i = 1, nqmatoms
      write(iunit,'(a2,1x,3f25.16)') elementSymbol(qmtypes(i)), qmcoords(1:3,i)
    end do
    write(iunit,'()')

    ! When electrostatic embadding QM/MM is in use 
    ! write MM coordinates with point charges
    if ( nclatoms > 0 ) then
      do i = 1, nclatoms
        write(iunit,'(4f21.16)') clcoords(:,i)
      end do
      write(iunit,'()')

      ! Write a second time without charges
      do i = 1, nclatoms
        write(iunit,'(4f21.16)') clcoords(:3,i)
      end do

      close(iunit)
    end if

   write(iunit,'()')
   close(iunit, iostat=ierr)

   if ( ierr /= 0 ) then
     call sander_bomb('write_inpfile (qm2_extern_gau_module)', &
       'Error closing Gaussian runfile after writing', &
       'Will quit now.')
   end if
   first_call = .false.

   call debug_exit_function( 'write_inpfile', module_name, gau_nml%debug )

  end subroutine write_inpfile

  subroutine read_results( datfile, fchkfile, nqmatoms, escf, dxyzqm,&
       nclatoms, dxyzcl, dipxyz, dipole, do_grad, debug )

    implicit none

    character(len=*), intent(in)  :: datfile
    character(len=*), intent(in)  :: fchkfile
    integer, intent(in)           :: nqmatoms, nclatoms
    _REAL_, intent(out)           :: escf, dxyzqm(3,nqmatoms), & 
                                     dxyzcl(3,nclatoms) ! dxyzcl returns containing the electric field at x,y,z
    _REAL_, intent(out)           :: dipxyz(3), dipole
    logical, intent(in)           :: do_grad
    integer, intent(in)           :: debug

    _REAL_ :: self_energy = 0 ! Temporary variable to hold self energy of point charges
    integer :: ios, i
    integer, parameter :: iunit = 351
    character(len=256) :: read_buffer

    call debug_enter_function( 'read_results', module_name, debug )

    open(iunit, file=fchkfile, status='old', iostat=ios)
    if ( ios /= 0 ) then
      call sander_bomb('read_results (qm2_extern_gau_module)', &
        'Error opening Gaussian formatted checkpoint file '//datfile//' (expected in same dir as input file).', &
        'Will quit now')
    end if

    do
      read (iunit, '(a)', iostat = ios) read_buffer
      ! End of file; nothing left to read
      if (ios < 0) then
        exit
      end if 
    
      ! Store SCF energy to escf 
      if ( read_buffer(1:12) == 'Total Energy' ) then
        ! Read value after the first "R"
        read(read_buffer(index(read_buffer,'R')+1:),*) escf
      end if
       
      ! Read forces on QM atoms
      if ( do_grad ) then
         if ( read_buffer(1:18) == 'Cartesian Gradient' ) then
            read (iunit, *) dxyzqm(:,:)
         end if
      end if

      ! Read dipole moment to vars dipxyz and dipole
      ! Note: This dipole moment is in units of Debye
      if (read_buffer(1:13) == 'Dipole Moment') then
         read (iunit, *) dipxyz(:)
      end if

    end do
    
    close(iunit)

    dipole = dsqrt ( dipxyz(1)*dipxyz(1) + dipxyz(2)*dipxyz(2) + dipxyz(3)*dipxyz(3) )

    ! QM/MM electrostatic embedding
    ! Required data is not on gaussian checkpoint file
    if ( nclatoms > 0 ) then

       open(iunit, file=datfile, status='old', iostat=ios)
       if ( ios /= 0 ) then
          call sander_bomb('read_results (qm2_extern_gau_module)', &
               'Error opening Gaussian log file '//datfile//' (expected in same dir as input file).', &
               'Will quit now')
       end if

       do
          read (iunit, '(a)', iostat = ios) read_buffer
          ! End of file; nothing left to read
          if (ios < 0) then
             exit
          end if
    
          ! Retrieve self energy of charges (must add to escf)
          if (read_buffer(1:30) == ' Self energy of the charges =') then
             read(read_buffer(index(read_buffer,'=')+1:),*) self_energy
          end if

          if ( do_grad ) then
             ! Read efield at CL atoms into dxyzcl
             if ( read_buffer(1:10) == '    1 Atom') then
                ! Skip QM atoms
                do  i = 1, nqmatoms-1
                   read(iunit, '(a)') read_buffer
                end do
                ! Read into dxyzcl
                do i = 1, nclatoms
                   read (iunit, '(a)', iostat = ios) read_buffer
                   read(read_buffer(25:),*) dxyzcl(:,i ) 
                end do
             end if
          end if

       end do
    
       close(iunit)
     
       escf = escf - self_energy

    end if

    call debug_exit_function( 'read_results', module_name, debug )

  end subroutine read_results

end module qm2_extern_gau_module
