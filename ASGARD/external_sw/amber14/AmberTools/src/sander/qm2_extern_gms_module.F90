#include "../include/dprec.fh"
module qm2_extern_gms_module
! ----------------------------------------------------------------
! Interface for GAMESS based QM MD 
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
!       January  2011
!
! ----------------------------------------------------------------

  use qm2_extern_util_module, only: debug_enter_function, debug_exit_function

  implicit none

  private
  public :: get_gms_forces

  character(len=*), parameter :: module_name = "qm2_extern_gms_module"

  type gms_nml_type
     character(len=20) :: basis
     character(len=20) :: method
     character(len=20) :: gms_version
     _REAL_ :: scf_conv
     integer :: nrad
     integer :: nleb
     integer :: maxit
     integer :: num_threads
     integer :: ntpr
     integer :: mwords
     integer :: debug
     logical :: chelpg
     logical :: dipole
     logical :: use_template
     
     ! Deprecated
     integer :: charge
     integer :: spinmult
  end type gms_nml_type

contains

  ! ------------------------------------
  ! Get QM energy and forces from GAMESS
  ! ------------------------------------
  subroutine get_gms_forces( do_grad, nstep, ntpr_default, id, &
       natoms, coords, qmtypes, escf, dxyzqm, charge, spinmult )

    use qm2_extern_util_module, only: print_results, check_installation, write_dipole, write_charges
    use constants, only: CODATA08_AU_TO_KCAL, CODATA08_A_TO_BOHRS, ZERO

    logical, intent(in) :: do_grad            ! Return gradient/not
    integer, intent(in) :: nstep              ! MD step number
    integer, intent(in) :: ntpr_default       ! frequency of printing
    character(len=3), intent(in) :: id        ! ID number for PIMD or REMD
    integer, intent(in) :: qmtypes(natoms)    ! QM atom types (nuclear charge in au
    integer, intent(in) :: natoms             ! Total number of atoms
    _REAL_,  intent(in) :: coords(3,natoms)   ! Sander coordinate data
    _REAL_, intent(out) :: escf               ! SCF energy
    _REAL_, intent(out) :: dxyzqm(3,natoms)   ! SCF QM force
    integer, intent(in) :: charge, spinmult   ! Charge and spin multiplicity

    _REAL_ :: dipxyz(3), dipole, charges(natoms) ! Dipole and charges from GAMESS

    type(gms_nml_type), save :: gms_nml
    logical, save :: first_call = .true.
    integer :: i
    integer :: printed =-1 ! Used to tell if we have printed this step yet 
                           ! since the same step may be called multiple times

    character(len=512)          :: call_buffer
    character(len=*), parameter :: program= 'rungms'
    character(len=*), parameter :: basename = 'gms_job'
    character(len=*), parameter :: inpext = '.inp'
    character(len=*), parameter :: outext = '.out'
    character(len=*), parameter :: datext = '.dat'
    character(len=*), parameter :: logext = '.log'
    character(len=*), parameter :: dipext = '.dip'
    character(len=*), parameter :: chgext = '.chg'
    character(len=*), parameter :: tplext = '.tpl'
    character(len=14) :: inpfile, outfile, datfile, logfile, dipfile, chgfile, tplfile
    ! Need to prepend subdirectory if doing REMD, PIMD or multi-region QM/MM. 
    !   This is triggered if 'id' is defined (not empty). 
    character(len=25)           :: subdir 

    ! assemble input - / output data filenames
    inpfile = basename//trim(id)//inpext
    outfile = basename//trim(id)//outext
    datfile = basename//trim(id)//datext
    logfile = basename//trim(id)//logext
    dipfile = basename//trim(id)//dipext
    chgfile = basename//trim(id)//chgext
    tplfile = basename//tplext

    ! Setup on first program call
    if ( first_call ) then
      first_call = .false.
      write (6,'(/,a,/)') '   >>> Running QM calculation with GAMESS <<<'
      call system('rm -f '//inpfile//' '//outfile//' '//datfile//' '//logfile)
      call system('rm -f '//dipfile//' '//chgfile)

      call check_installation( program, id, .true., gms_nml%debug )
      call get_namelist( ntpr_default, gms_nml )
      call print_namelist( gms_nml ) 

      write (6,'(80a)') ('-', i=1,80)
      write (6,'(a)') '   4.  RESULTS'
      write (6,'(80a)') ('-', i=1,80)
    end if

    call write_inpfile(trim(inpfile), trim(tplfile), coords, qmtypes, natoms, &
                       gms_nml, do_grad, charge, spinmult )

    ! Run GAMESS with file inpfile, version gms_version, using num_threads threads

    ! BPK-Separate runs into different directories if we are doing PIMD or REMD
    subdir = ''
    call_buffer = ''
    if (trim(id)/='') then 
      subdir = './'//trim(id)//'/'
      call_buffer = ' mkdir -p '//trim(subdir)//'; cd '//trim(subdir)//&
        '; mv ../'//inpfile//' .;'
    end if
    
    ! We assume that the datfile will be written to this same directory
    ! First remove a datfile that may potentially have been left over
    call_buffer = trim(call_buffer)//' rm -f '//datfile//';'

    ! Add program call to call_buffer
    write(call_buffer,'(a,a,i0,a)') &
         trim(call_buffer), &
         program//' '//trim(inpfile)//' '//gms_nml%gms_version//' ', &
         gms_nml%num_threads, &
         ' > '//trim(outfile)//' 2> '//trim(logfile)
 
    call system(trim(call_buffer))


    ! Call read_results - retrieve data from GAMESS .dat and .out files
    call read_results( trim(subdir)//trim(datfile), trim(subdir)//trim(outfile), gms_nml%chelpg, gms_nml%dipole,&
      natoms, gms_nml%method, escf, dxyzqm, dipxyz, dipole, charges, do_grad, gms_nml%debug )

    if ( gms_nml%ntpr > 0 .and. mod(nstep, gms_nml%ntpr) == 0 ) then
      if ( printed /= nstep ) then
        if ( gms_nml%dipole ) then
          ! Call write_dipole  - will write output to dipfile
          call write_dipole( trim(dipfile), dipxyz, dipole, gms_nml%debug )
        end if
        if ( gms_nml%chelpg ) then
          ! Call write_charges - will write output to chgfile
          call write_charges( trim(chgfile), charges, gms_nml%debug )
        end if
        printed = nstep
      end if
    end if

    call system('mv '//trim(subdir)//inpfile//' '//trim(subdir)//'old.'//inpfile)
    call system('mv '//trim(subdir)//datfile//' '//trim(subdir)//'old.'//datfile)
    call system('mv '//trim(subdir)//outfile//' '//trim(subdir)//'old.'//outfile)
    call system('mv '//trim(subdir)//logfile//' '//trim(subdir)//'old.'//logfile)

    ! Convert Hartree/Bohr -> kcal/(mol*A)
    if ( do_grad ) then
      dxyzqm(:,:) = dxyzqm(:,:) * CODATA08_AU_TO_KCAL * CODATA08_A_TO_BOHRS
    else
      dxyzqm = ZERO
    end if

    ! Convert Hartree -> kcal/mol
    escf = escf * CODATA08_AU_TO_KCAL

    call print_results( 'qm2_extern_gms_module', escf, natoms, dxyzqm,&
      gms_nml%debug )

  end subroutine get_gms_forces

  ! -----------------------------------------------
  ! Read GAMESS gms namelist values from file mdin,
  ! use default values if none are present.
  ! -----------------------------------------------
  subroutine get_namelist(ntpr_default, gms_nml)
    
    implicit none

    integer, intent(in) :: ntpr_default
    type(gms_nml_type), intent(out) :: gms_nml

    character(len=20) :: basis, method, gms_version
    _REAL_ :: scf_conv
    integer :: nrad, nleb, maxit, num_threads, ntpr, mwords, &
      debug, chelpg, dipole, use_template
    integer :: charge, spinmult ! deprecated
    namelist /gms/ basis, method, gms_version, scf_conv, &
      nrad, nleb, maxit, num_threads, ntpr, mwords, &
      debug, chelpg, dipole, use_template,&
      charge, spinmult

    integer :: ios

    ! default values
    basis        = '6-31G*'
    method       = 'BP86'
    gms_version  = '00'    ! This is the two number version chosen during GAMESS compilation
    scf_conv     = 1.0D-06 ! SCF (1.0D-05 is GAMESS default)
    nrad         = 96      ! radial grid points in DFT grid (96 is GAMESS default)
    nleb         = 590     ! Angular Lebedev DFT grid points (302 is GAMESS default)
    maxit        = 50      ! Maximum number of SCF iterations
    num_threads  = 1
    ntpr         = ntpr_default
    mwords       = 50
    debug        = 0
    chelpg       = 0
    dipole       = 0
    use_template = 0
    
    ! These are now deprecated and should be specified in the &qmmmm namelist
    charge   = -351
    spinmult = -351

    ! Read namelist
    rewind 5
    read (5, nml=gms, iostat=ios)

    if ( ios > 0 ) then
       call sander_bomb('get_namelist (qm2_extern_gms_module)', &
            '&gms namelist read error', &
            'Please check your input.')
    else if ( ios < 0 ) then
       write(6,'(a,/,a)') '&gms namelist read encountered end of file', &
            'Please check your input if the calculation encounters a problem'
    end if

    if ( charge /= -351 .or. spinmult /= -351 ) then
      call sander_bomb('get_namelist (qm2_extern_gms_module)', &
        'The charge and spin keywords are deprecated', &
        'Please specify charge (qmcharge) and spin multiplicity (spin) in the &qmmm namelist.')
    end if

    ! Assign namelist values to gms_nml data type
    gms_nml%basis        = basis   ! Will be split and parsed later in write_inpfile
    gms_nml%method       = method
    gms_nml%gms_version  = gms_version
    gms_nml%scf_conv     = scf_conv
    gms_nml%nrad         = nrad
    gms_nml%nleb         = nleb
    gms_nml%maxit        = maxit
    gms_nml%num_threads  = num_threads
    gms_nml%ntpr         = ntpr
    gms_nml%mwords       = mwords
    gms_nml%debug        = debug

    if ( chelpg == 0 ) then
      gms_nml%chelpg = .false.
    else if ( chelpg == 1 ) then
      gms_nml%chelpg = .true.
    else
      call sander_bomb('get_namelist (qm2_extern_gms_module)', &
        '&gms chelpg value not allowed', &
        'Please check your input. chelpg can only be 0 or 1.')
    end if

    if ( dipole == 0 ) then
      gms_nml%dipole = .false.
    else if ( dipole == 1 ) then
      gms_nml%dipole = .true.
    else
      call sander_bomb('get_namelist (qm2_extern_gms_module)', &
        '&gms dipole value not allowed', &
        'Please check your input. dipole can only be 0 or 1.')
    end if

    if ( use_template == 0 ) then
      gms_nml%use_template = .false.
    else if ( use_template == 1 ) then
      gms_nml%use_template = .true.
    else
      call sander_bomb('get_namelist (qm2_extern_gms_module)', &
        '&gms use_template value not allowed', &
        'Please check your input. use_template can only be 0 or 1.')
    end if
  end subroutine get_namelist

  ! ------------------------------
  ! Print GAMESS namelist settings
  ! ------------------------------
  subroutine print_namelist(gms_nml)

    implicit none
    type(gms_nml_type), intent(in) :: gms_nml

    write(6, '(/,a)')     "| &gms"
    write(6, '(2a)')      "|   basis        = ", gms_nml%basis
    write(6, '(2a)')      "|   method       = ", gms_nml%method
    write(6, '(2a)')      "|   gms_version  = ", gms_nml%gms_version
    write(6, '(a,es10.2)')"|   scf_conv     = ", gms_nml%scf_conv
    write(6, '(a,i3)')    "|   nrad         = ", gms_nml%nrad
    write(6, '(a,i3)')    "|   nleb         = ", gms_nml%nleb
    write(6, '(a,i3)')    "|   maxit        = ", gms_nml%maxit
    write(6, '(a,i3)')    "|   num_threads  = ", gms_nml%num_threads
    write(6, '(a,i0)')    "|   ntpr         = ", gms_nml%ntpr
    write(6, '(a,i0)')    "|   mwords       = ", gms_nml%mwords
    write(6, '(a,l)')     "|   chelpg       = ", gms_nml%chelpg
    write(6, '(a,l)')     "|   dipole       = ", gms_nml%dipole
    write(6, '(a,l)')     "|   use_template = ", gms_nml%use_template
    write(6,'(a)')        "| /"

  end subroutine print_namelist
  
  ! -------------------------------
  ! Write the input file for GAMESS
  ! -------------------------------
  subroutine write_inpfile(inpfile, tplfile, coords, qmtypes, natoms, &
                           gms_nml, do_grad, charge, spinmult )

    use UtilitiesModule, only: Upcase
    use ElementOrbitalIndex, only : elementSymbol

    implicit none

    character (len=*) , intent(in) :: inpfile, tplfile
    _REAL_            , intent(in) :: coords(:,:)
    integer           , intent(in) :: qmtypes(:)
    integer           , intent(in) :: natoms
    type(gms_nml_type), intent(in) :: gms_nml
    logical, intent(in)            :: do_grad
    integer           , intent(in) :: charge, spinmult

    character(len=256) :: gbasis, read_buffer
    character(len=3)   :: scftyp
    integer            :: ngauss, nstars, nplus, i, j, ios
!   integer            :: norb
    integer, parameter :: iurun = 10, iudat = 11
    logical            :: pople, dft
!   logical            :: found
!   logical, save      :: first_call = .true.
    character(len=10)  :: runtyp = 'GRADIENT' ! Normal run type

    ! We assume that we are doing DFT if we do neither HF or MP2
    dft = (index(Upcase(gms_nml%method), 'MP2') == 0) .and. ( index(Upcase(gms_nml%method), 'HF') == 0 )

    call debug_enter_function( 'write_inpfile', module_name, gms_nml%debug )
    
    if ( gms_nml%use_template ) then
      call system('cp '//tplfile//' '//inpfile)
    end if

    open(iurun, file=inpfile, position='append', iostat=ios)
    if ( ios > 0 ) then
       call sander_bomb('write_inpfile (qm2_extern_gms_module)', &
            'Error opening GAMESS input file '//inpfile//' for writing', &
            'Will quit now')
    end if

    if (.not. gms_nml%use_template) then
      write (iurun, '(a,/,a)') '! Run using SANDER interface for GAMESS', '!'

      ! $SYSTEM card
      write (iurun,'(a,/,a, i0,/,a,/)') &
           ' $SYSTEM',  &
           'MWORDS=',gms_nml%mwords, &
           ' $END'

      ! $CONTRL card
      if ( spinmult > 1 ) then
         scftyp = 'UHF'
      else
         scftyp = 'RHF'
      end if

      if ( .not. do_grad ) runtyp = 'ENERGY'

      write (iurun,'(a,/,2a,/,a,/,a,i0,/,a,i0,/,a,i0,/,a,/,a,/,a)') &
           ' $CONTRL'        ,&
           'SCFTYP=',scftyp  ,&
           'RUNTYP='//runtyp ,&
           'ICHARG='         ,charge,&
           'MULT='           ,spinmult,&
           'MAXIT='          ,gms_nml%maxit,&
           'COORD=UNIQUE'    ,&
           'UNITS=ANGS'      ,&
           'ISPHER=1'
      if ( Upcase(gms_nml%method) == 'MP2' ) then
         write (iurun,'(a,/)') &
              'MPLEVL=2'
      else if ( dft ) then
         ! We are assuming that if we are not doing MP2 or HF, then we do DFT
         write (iurun,'(2a)') &
              'DFTTYP=', trim(gms_nml%method)
      end if
      write (iurun,'(a,/)') &
           ' $END'

      ! $BASIS card
      pople = .true.
      
      ! First, count and filter out '+' and '*' characters and store result in read_buffer
      nstars = 0
      nplus = 0
      j=0
      read_buffer=''
        ! Search through each character of our basis name and write to 
        !  read_buffer if it does not contain a star or plus
        do i=1, len(gms_nml%basis)
        if ( gms_nml%basis(i:i)=='*') then
          nstars=nstars+1
        else if (gms_nml%basis(i:i)=='+') then
          nplus=nplus+1
        else
          j=j+1
          read_buffer(j:j)=gms_nml%basis(i:i)
        end if
      end do
      
      ! Now, see if we fall under a pople case:
      select case (Upcase(trim(read_buffer)))
      case ('STO-3G')
         gbasis = 'STO'
         ngauss = 3
      case ('STO-6G')
         gbasis = 'STO'
         ngauss = 6
      case ('3-21G')
         gbasis = 'N21'
         ngauss = 3
      case ('6-31G')
         gbasis = 'N31'
         ngauss = 6
      case ('6-311G')
         gbasis = 'N311'
         ngauss = 6
      case default
         pople = .false.
         nstars = 0
         nplus = 0
         gbasis = gms_nml%basis
      end select

      write(iurun,'(a,/,2a)') &
           ' $BASIS', &
           'GBASIS=', gbasis
      if ( pople ) then
         write(iurun,'(a,i1)') 'NGAUSS=', ngauss
         if ( nstars > 0 ) then
            write(iurun,'(a,i0)') 'NDFUNC=1'
            if ( nstars > 1 ) then
               write(iurun,'(a,i0)') 'NPFUNC=1'
            end if
         end if
         if ( nplus > 0 ) then
            write(iurun,'(a,i0)') 'DIFFSP=.TRUE.'
            if ( nplus > 1 ) then
               write(iurun,'(a,i0)') 'DIFFS=.TRUE.'
            end if
         end if
      end if

      write(iurun,'(a,/)')' $END'

      ! $DFT card
      if ( dft ) then
         write(iurun,'(a,/,a,/,a,i3,/,a,i3,/,a,/)') &
              ' $DFT'      , &
              'METHOD=GRID', &
              'NRAD=', gms_nml%nrad, &
              'NLEB=', gms_nml%nleb, &
              ' $END'
      end if

      ! $SCF card
      write(iurun,'(a,/,a,/,a,E22.16,/,a,/)') &
           ' $SCF'        , &
           'DIRSCF=.TRUE.', &
           'CONV='        , gms_nml%scf_conv, &
           ' $END'

      ! $PDC card
      if ( gms_nml%chelpg ) then
         write(iurun, '(a,/,a,/,a,/,a,/,a,/,a,/,a,/,a,/,a,/)') &
              ' $ELPOT', &
              'IEPOT = 1', &
              'WHERE = PDC', &
              ' $END', &
              ' ', &
              ' $PDC', &
              'PTSEL = CHELPG', &
              'CONSTR = CHARGE', &
              ' $END'
      end if
    end if
    ! $DATA card
    write(iurun,'(a,/,a,/,a)') &
         ' $DATA', &
         'Atom cards are inserted below; C1 = no symmetry', &
         'C1'
    do i = 1, natoms
       write(iurun,'(a2,1x,i0,1x,3f21.16)') &
            elementSymbol(qmtypes(i)), &
            qmtypes(i)             , &
            coords(1:3,i)
    end do
    write(iurun,'(a,/)') ' $END'

    ! $VEC card and $GUESS cards (for restart)
    !if (first_call) then
    !
    !   first_call = .false.
    !
    !else
    !
    !   open(unit=iudat, file=datfile, status='old', iostat=ios)
    !   if ( ios /= 0 ) then
    !      call sander_bomb('write_inpfile (qm2_extern_gms_module)', &
    !           'Error opening GAMESS punch file '//trim(datfile)//' from previous step', &
    !           'Will quit now')
    !   end if
    !   
    !   found = .false.
    !   do
    !      read (iudat, '(a)', iostat = ios) read_buffer
    !      ! End of file; data not found
    !      if (ios < 0) then
    !         call sander_bomb('write_inpfile (qm2_extern_gms_module)', &
    !              'Error reading GAMESS data from file '//trim(datfile)//' ($VEC card not found)', &
    !              'Will quit now')
    !      end if
    !      if (trim(read_buffer) == ' $VEC') then
    !         found = .true.
    !         write (iurun, '(a)') read_buffer
    !         do
    !            read (iudat, '(a)', iostat = ios) read_buffer
    !            ! End of file; End block not found
    !            if (ios < 0) then
    !               call sander_bomb('write_inpfile (qm2_extern_gms_module)', &
    !                    'Error reading GAMESS data from file '//trim(datfile)//' ($END of $VEC card not found)', &
    !                    'Will quit now')
    !            end if
    !            if (trim(read_buffer) == ' $END') then
    !               exit
    !            else
    !               write (iurun, '(a)') read_buffer
    !               read (read_buffer, *) norb
    !            end if
    !         end do
    !      end if
    !      if (found) exit
    !   end do
    !
    !   close(iudat)
    !   
    !   write (iurun, '(2(a,/),a,i0,/,a,/)') &
    !        ' $GUESS'       , &
    !        'GUESS = MOREAD', &
    !        'NORB=',norb    , &
    !        ' $END'
    !
    !end if
    
    close(iurun)

    call debug_exit_function( 'write_inpfile', module_name, gms_nml%debug )

  end subroutine write_inpfile

  ! -------------------
  ! Read GAMESS results
  ! -------------------
  subroutine read_results( datfile, outfile, chelpg, do_dipole, natoms, method, &
  escf, dxyzqm, dipxyz, dipole, charges, do_grad, debug )

    implicit none

    character(len=*), intent(in) :: datfile
    character(len=*), intent(in) :: outfile
    logical, intent(in) :: chelpg, do_dipole
    integer, intent(in) :: natoms
    character(len=*), intent(in) :: method
    _REAL_, intent(out) :: escf, dxyzqm(3,natoms)
    _REAL_, intent(out) :: dipxyz(3), dipole, charges(natoms)
    logical, intent(in) :: do_grad
    integer, intent(in) :: debug

    integer :: ios, i, ifound
    integer, parameter :: iunit = 351, junit = 352
    character(len=256) :: read_buffer

    call debug_enter_function( 'read_results', module_name, debug )

    ! read energy and gradient
    open(iunit, file=datfile, status='old', iostat=ios)
    if ( ios /= 0 ) then
      call sander_bomb('read_results (qm2_extern_gms_module)', &
        'Error opening GAMESS punch file '//datfile//' (expected in same dir as input file).', &
        'Will quit now')
    end if

    do
      read (iunit, '(a)', iostat = ios) read_buffer
      ! End of file; data not found
      if (ios < 0) then
        call sander_bomb('read_results (qm2_extern_gms_module)', &
          'Error reading GAMESS data from file '//datfile//' ($GRAD card not found)', &
          'Will quit now - please check the GAMESS output file.')
      end if
      if ( do_grad ) then
        if (trim(read_buffer) == ' $GRAD') then
          ! Read energy line
          read (iunit, '(2x,f22.10)') escf
          ! Read gradient data
          do i = 1, natoms
            read (iunit, '(15x,3f20.16)') dxyzqm(:,i)
          end do
          exit
        end if
      ! No gradients; just get energy
      else if (index(read_buffer,'--- CLOSED SHELL ORBITALS ---')>0) then
        read (iunit, '(a)') read_buffer ! Don't need this line
        read (iunit, '(15x,f22.10)') escf
        exit
      end if
    end do

    close(iunit)

    ! open output file for reading dipole moment and charges
    open(junit, file=outfile, status='old', iostat=ios)
    if ( ios /= 0 ) then
      call sander_bomb('read_results (qm2_extern_gms_module)', &
        'Error opening GAMESS punch file '//outfile//' (expected in same dir as input file).', &
        'Will quit now')
    end if

    if ( do_dipole ) then
      ! read dipole moment
      ! the first dipole moment is for the reference wave function
      if ( method == 'MP2' ) then
        ifound = 0
      else
        ifound = 1
      end if
      do
        read (junit, '(a)', iostat = ios) read_buffer
        if (ios < 0) then
          ! End of file; data not found
          call sander_bomb('read_results (qm2_extern_gms_module)', &
            'Error reading GAMESS data from file '//outfile//' (Dipole moment DX etc not found)', &
            'Will quit now')
         end if
         if (index(read_buffer,'DX          DY          DZ         /D/  (DEBYE)') > 0) then
           if ( ifound == 0 ) then
             ifound = 1
             cycle
           else
             read (junit, '(x,4(x,f11.6))') dipxyz(:), dipole
             exit
           end if
         end if
      end do
    end if

    ! read CHELPG charges
    if (chelpg) then
      do
        read (junit, '(a)', iostat = ios) read_buffer
        if (ios < 0) then
          ! End of file; data not found
          call sander_bomb('read_results (qm2_extern_gms_module)', &
            'Error reading GAMESS data from file '//outfile//' (NET CHARGES not found)', &
            'Will quit now')
        end if
        if (trim(read_buffer) == ' NET CHARGES:') then
          do i = 1, 3
            read(junit,'(a)') read_buffer
          end do
          do i = 1, natoms
            read (junit,'(20x,f7.4)') charges(i)
          end do
          exit
        end if
      end do
    end if
       
    ! close output file after reading dipole moment and charges
    close(junit)

    call debug_exit_function( 'read_results', module_name, debug )

  end subroutine read_results

end module qm2_extern_gms_module
