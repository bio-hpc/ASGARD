#include "../include/dprec.fh"
module qm2_extern_tc_module
! ----------------------------------------------------------------
! Interface for TeraChem based QM and QM/MM MD 
!
! Currently supports:
! pure QM
! QM/MM with cutoff for QM-MM electrostatics under periodic
! boundary conditions
!
! Author: Andreas Goetz (agoetz@sdsc.edu)
!
! Date: November 2010
!
! ----------------------------------------------------------------

  use qm2_extern_util_module, only: debug_enter_function, debug_exit_function

  implicit none

  private
  public :: get_tc_forces, tc_finalize
  logical, save :: do_mpi = .false.  ! Used in finalize subroutine

  character(len=*), parameter :: module_name = "qm2_extern_tc_module"

  type tc_nml_type
     character(len=20) :: basis
     character(len=20) :: method
     character(len=20) :: precision
     character(len=20) :: executable
     character(len=20) :: dftd
     character(len=20) :: guess
     character(len=20) :: cis
     character(len=20) :: charge_analysis
     _REAL_ :: threall
     _REAL_ :: convthre
     integer :: maxit
     integer :: dftgrid
     integer :: ngpus
     integer, dimension(:), pointer :: gpuids => null()
     integer :: cisnumstates
     integer :: cistarget
     integer :: mpi
     integer :: ntpr
     integer :: debug
     logical :: dipole
     logical :: use_template
     
     ! Deprecated
     integer :: charge
     integer :: spinmult
  end type tc_nml_type

  integer, save         :: newcomm ! Initialized in mpi_init subroutine

contains

  ! --------------------------------------
  ! Get QM energy and forces from TeraChem
  ! --------------------------------------
  subroutine get_tc_forces( do_grad, nstep, ntpr_default, id, nqmatoms, qmcoords,&
       qmtypes, nclatoms, clcoords, escf, dxyzqm, dxyzcl, charge, spinmult )

    use qm2_extern_util_module, only: print_results, check_installation, write_dipole, write_charges
    use constants, only: CODATA08_AU_TO_KCAL, CODATA08_A_TO_BOHRS, ZERO

    logical, intent(in) :: do_grad              ! Return gradient/not
    integer, intent(in) :: nstep                ! MD step number
    integer, intent(in) :: ntpr_default         ! frequency of printing
    character(len=3), intent(in) :: id          ! ID number for PIMD or REMD
    integer, intent(in) :: nqmatoms             ! Number of QM atoms
    _REAL_,  intent(in) :: qmcoords(3,nqmatoms) ! QM coordinates
    integer, intent(in) :: qmtypes(nqmatoms)    ! QM atom types (nuclear charge in au)
    integer, intent(in) :: nclatoms             ! Number of MM atoms
    _REAL_,  intent(in) :: clcoords(4,nclatoms) ! MM coordinates and charges in au
    _REAL_, intent(out) :: escf                 ! SCF energy
    _REAL_, intent(out) :: dxyzqm(3,nqmatoms)   ! SCF QM force
    _REAL_, intent(out) :: dxyzcl(3,nclatoms)   ! SCF MM force
    integer, intent(in) :: charge, spinmult     ! Charge and spin multiplicity

    _REAL_              :: dipmom(4,3)          ! Dipole moment {x, y, z, |D|}, {QM, MM, TOT}
    _REAL_              :: qmcharges(nqmatoms)  ! QM charges from population analysis

    type(tc_nml_type), save :: tc_nml
    logical, save :: first_call = .true.
    integer :: i
    integer :: printed =-1 ! Used to tell if we have printed this step yet 
                           ! since the same step may be called multiple times
    character(len=150) :: call_buffer
    character(len=*),  parameter :: basename = 'tc_job'
    character(len=*),  parameter :: runext = '.inp'
    character(len=*),  parameter :: datext = '.dat'
    character(len=*),  parameter :: dipext = '.dip'
    character(len=*),  parameter :: chgext = '.chg'
    character(len=*),  parameter :: tplext = '.tpl'
    character(len=14) :: inpfile, datfile, dipfile, chgfile, crdfile, ptcfile, tplfile
    ! Need to prepend subdirectory if doing REMD, PIMD or multi-region QM/MM. 
    !   This is triggered if 'id' is defined (not empty). 
    character(len=25)            :: subdir 
    
    ! assemble input - / output data filenames
    inpfile = basename//trim(id)//runext 
    datfile = basename//trim(id)//datext
    dipfile = basename//trim(id)//dipext 
    chgfile = basename//trim(id)//chgext 
    crdfile = 'inpfile'//trim(id)//'.xyz' 
    ptcfile = 'ptchrg'//trim(id)//'.xyz'
    tplfile = basename//tplext 
    
    ! Setup on first program call
    if ( first_call ) then
      first_call = .false.
      write (6,'(/,a,/)') '   >>> Running QM calculation with TeraChem <<<'
      call get_namelist( ntpr_default, tc_nml )
      call check_installation( trim(tc_nml%executable), id, .true., tc_nml%debug )
      call print_namelist(tc_nml)

      write (6,'(80a)') ('-', i=1,80)
      write (6,'(a)') '   4.  RESULTS'
      write (6,'(80a)') ('-', i=1,80)
      call system('rm -f '//dipfile//' '//chgfile)
    end if

#ifdef MPI
# ifndef MPI_1
    if (tc_nml%mpi==1 ) then ! Do mpi (forced to 0 ifndef MPI)
      call mpi_hook( trim(tplfile), nqmatoms, qmcoords, qmtypes, nclatoms, clcoords,&
        tc_nml, escf, dxyzqm, dxyzcl, dipmom, qmcharges, do_grad, id, charge, spinmult )
    else
# else
    ! If we are using MPI 1.x the code will not compile since
    ! MPI_LOOKUP_NAME is part of the MPI 2 standard, so  just quit
    if (tc_nml%mpi==1 ) then 
    call sander_bomb('(qm2_extern_tc_module)', &
      '&unsupported MPI version', &
      'Will quit now.')
    else
# endif
#endif
      
       call system('rm -f '//inpfile)
       call write_inpfile( trim(inpfile), trim(crdfile), trim(ptcfile), trim(tplfile), &
            nqmatoms, qmcoords, qmtypes, nclatoms, clcoords, &
            tc_nml, do_grad, charge, spinmult )

       call_buffer=''
       subdir=''
       if ( trim(id)/='' ) then
         subdir='./'//trim(id)//'/'
         call_buffer=' mkdir -p '//trim(subdir)//&
                     '; cd '//trim(subdir)//&
                     '; mv ../'//inpfile//' .;'//&
                      ' mv ../'//crdfile//' .;'
         if ( nclatoms > 0 ) then
           call_buffer=trim(call_buffer)//' mv ../'//ptcfile//' .;'
         end if
       end if

      ! Run TeraChem with file inpfile
       call system('rm -f '//datfile)
       write(call_buffer,'(2a)') &
         trim(call_buffer),trim(tc_nml%executable)//' '//trim(inpfile)//' > '//datfile
      
       call system(trim(call_buffer))

       ! If working in a subdirectory, move datfile back for reading
       if ( trim(id)/='' ) then
         call system('mv '//trim(subdir)//datfile//' .;')
       end if

       ! Read TeraChem results
       call read_results( trim(datfile), nqmatoms, nclatoms, escf, dxyzqm, dxyzcl, &
            dipmom, qmcharges, tc_nml%charge_analysis, do_grad, tc_nml%debug )
 
       call system( 'mv '//trim(subdir)//trim(inpfile)//' '//trim(subdir)//'old.'//inpfile )
       call system( 'mv '//trim(datfile)//' '//trim(subdir)//'old.'//datfile )
#ifdef MPI
    end if
#endif

    ! Write dipole and charges to file
    if ( tc_nml%ntpr > 0 .and. mod(nstep, tc_nml%ntpr) == 0 ) then
       if ( printed /= nstep ) then
          printed = nstep
          if ( tc_nml%dipole ) then
             ! using util module's write dipole, writing only QM dipole moment
             call write_dipole( trim(dipfile), dipmom(1:3,1), dipmom(4,1), tc_nml%debug )
          end if
          if ( trim(tc_nml%charge_analysis) /= 'NONE' ) then
             call write_charges( trim(chgfile), qmcharges, tc_nml%debug )
          end if
       end if
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

    call print_results( 'qm2_extern_tc_module', escf, nqmatoms, dxyzqm,&
      tc_nml%debug, nclatoms, dxyzcl )


  end subroutine get_tc_forces

  ! -----------------------------------------------
  ! Read TeraChem tc namelist values from file mdin,
  ! use default values if none are present.
  ! -----------------------------------------------
  subroutine get_namelist( ntpr_default, tc_nml)

    use UtilitiesModule, only: Upcase
    implicit none

    integer, intent(in) :: ntpr_default
    type(tc_nml_type), intent(out) :: tc_nml
    character(len=20):: basis, method, dftd, precision, executable, guess, cis, charge_analysis
    _REAL_ :: threall, convthre
    integer, parameter :: maxgpus = 64
    integer :: maxit, dftgrid, ngpus, gpuids(maxgpus),  &
         cisnumstates, cistarget, mpi, ntpr, debug, dipole, use_template
    integer :: charge, spinmult ! deprecated
    namelist /tc/ basis, method, dftd, precision, executable, guess, cis, charge_analysis, &
      threall, convthre, &
      maxit, dftgrid, ngpus, gpuids, cisnumstates, cistarget, &
      mpi, ntpr, debug, dipole, use_template,&
      charge, spinmult
    integer :: i, ierr

    ! Default values
    basis           = '6-31g'
    method          = 'blyp'
    dftd            = 'no'
    precision       = 'mixed'
    executable      = 'terachem'
    guess           = 'scr/c0'
    cis             = 'no'
    charge_analysis = 'none'
    threall         = 1.0d-11
    convthre        = 3.0d-05
    maxit           = 100
    dftgrid         = 1
    ngpus           = 0 ! Use all available GPUs
    do i = 1, maxgpus
       gpuids(i) = i-1
    end do
    cisnumstates = 1
    cistarget    = 1
    mpi          = 1 ! Default to using MPI if available
    ntpr         = ntpr_default
    debug        = 0
    dipole       = 0
    use_template = 0

    ! These are now deprecated and should be specified in the &qmmmm namelist
    charge   = -351
    spinmult = -351

    ! Read namelist
    rewind 5
    read(5,nml=tc,iostat=ierr)

    if ( ierr > 0 ) then
       call sander_bomb('get_namelist (qm2_extern_tc_module)', &
            '&tc namelist read error', &
            'Please check your input.')
    else if ( ierr < 0 ) then
       write(6,'(a/a)') '&tc namelist read encountered end of file', &
            'Please check your input if the calculation encounters a problem'
    end if

    if ( charge /= -351 .or. spinmult /= -351 ) then
      call sander_bomb('get_namelist (qm2_extern_tc_module)', &
        'The charge and spin keywords are deprecated', &
        'Please specify charge (qmcharge) and spin multiplicity (spin) in the &qmmm namelist.')
    end if

    charge_analysis = Upcase(charge_analysis)
    if ( (trim(charge_analysis) /= 'NONE') .and. (trim(charge_analysis) /= 'MULLIKEN') ) then
       call sander_bomb('get_namelist (qm2_extern_tc_module)', &
        'Only Mulliken charge analysis is supported.', &
        'Please correct charge_analysis in the &qmmm namelist.')
    end if

    ! Assign namelist values to tc_nml data type
    tc_nml%basis           = Upcase(basis)
    tc_nml%method          = method
    tc_nml%dftd            = dftd
    tc_nml%precision       = precision
    tc_nml%executable      = executable
    tc_nml%guess           = guess
    tc_nml%cis             = cis
    tc_nml%charge_analysis = charge_analysis
    tc_nml%threall         = threall
    tc_nml%convthre        = convthre
    tc_nml%maxit           = maxit
    tc_nml%dftgrid         = dftgrid
    tc_nml%ngpus           = ngpus
    tc_nml%cisnumstates    = cisnumstates
    tc_nml%cistarget       = cistarget
    if ( ngpus > 0 ) then
       allocate ( tc_nml%gpuids(ngpus), stat = ierr )
       if ( ierr /= 0 ) then
          call sander_bomb('get_namelist (qm2_extern_tc_module)', &
               'Allocation error for gpuids(:)', &
               'Will quit now')
       end if
       tc_nml%gpuids(:) = gpuids(:ngpus)
    end if
#ifndef MPI
        if ( tc_nml%mpi == 1 ) then
          write(6,'(a)') '| Warning: mpi=1 selected but sander was not compiled with MPI support.'
          write(6,'(a)') '| Continuing with mpi=0'
        end if
        tc_nml%mpi         = 0 ! Can't pick MPI if not available 
#else
        tc_nml%mpi         = mpi
#endif

    ! Need this variable so we don't call MPI_Send in the finalize subroutine
    if (mpi==1 ) then
      do_mpi=.true.
    end if

    tc_nml%ntpr      = ntpr
    tc_nml%debug     = debug

    if ( dipole == 0 ) then
       tc_nml%dipole = .false.
    else if ( dipole == 1 ) then
       tc_nml%dipole = .true.
    else
       call sander_bomb('get_namelist (qm2_extern_tc_module)', &
            '&tc dipole value not allowed', &
            'Please check your input. dipole can only be 0 or 1.')
    end if

    if ( use_template == 0 ) then
       tc_nml%use_template = .false.
    else if ( use_template == 1 ) then
       tc_nml%use_template = .true.
    else
       call sander_bomb('get_namelist (qm2_extern_tc_module)', &
            '&tc use_template value not allowed', &
            'Please check your input. use_template can only be 0 or 1.')
    end if


  end subroutine get_namelist

  ! --------------------------------
  ! Print TeraChem namelist settings
  ! --------------------------------
  subroutine print_namelist( tc_nml )

    implicit none
    type(tc_nml_type), intent(in) :: tc_nml

    integer :: i, j, jstart, jend
    integer, parameter :: jstep = 10
    character(len=30) :: tmpstr

    write(6, '(/,a)')      '| &tc'
    write(6, '(2a)')       '|   basis           = ', tc_nml%basis
    write(6, '(2a)')       '|   method          = ', tc_nml%method
    write(6, '(2a)')       '|   dftd            = ', tc_nml%dftd
    write(6, '(2a)')       '|   precision       = ', tc_nml%precision
    write(6, '(2a)')       '|   executable      = ', tc_nml%executable
    write(6, '(2a)')       '|   guess           = ', tc_nml%guess
    write(6, '(2a)')       '|   cis             = ', tc_nml%cis
    write(6, '(2a)')       '|   charge_analysis = ', tc_nml%charge_analysis
    write(6, '(a,es10.2)') '|   threall         = ', tc_nml%threall
    write(6, '(a,es10.2)') '|   convthre        = ', tc_nml%convthre
    write(6, '(a,i4)')     '|   maxit           = ', tc_nml%maxit
    write(6, '(a,i4)')     '|   dftgrid         = ', tc_nml%dftgrid
    write(6, '(a,i4)')     '|   ngpus           = ', tc_nml%ngpus
    write(6, '(a,i4)')     '|   cisnumstates    = ', tc_nml%cisnumstates
    write(6, '(a,i4)')     '|   cistarget       = ', tc_nml%cistarget
    if ( tc_nml%ngpus > 0 ) then
       jstart = 1
       do i = 1, tc_nml%ngpus / jstep + 1
          if ( i == 1 ) then
             tmpstr =      '|   gpuids          = '
          else
             tmpstr =      '                  '
          end if
          jend = min ( (jstart + jstep - 1), tc_nml%ngpus )
          write(6,'(a,9999(i5))') tmpstr, (tc_nml%gpuids(j), j = jstart, jend)
          jstart = jstart + jstep
       end do
    end if
    write(6, '(a,i1)')     '|   mpi             = ', tc_nml%mpi
    write(6, '(a,i0)')     '|   ntpr            = ', tc_nml%ntpr
    write(6, '(a,i2)')     '|   debug           = ', tc_nml%debug
    write(6, '(a,l)')      '|   dipole          = ', tc_nml%dipole
    write(6, '(a,l)')      '|   use_template    = ', tc_nml%use_template
    write(6,'(a)')         '| /'

  end subroutine print_namelist

#if defined(MPI) && !defined(MPI_1)
  ! Perform MPI communications with terachem. Requires MPI 2.0 or above to use
  subroutine mpi_hook( tplfile, nqmatoms, qmcoords, qmtypes, nclatoms, clcoords,&
       tc_nml, escf, dxyzqm, dxyzcl, dipmom, qmcharges, do_grad, id, charge, spinmult )
    
    use ElementOrbitalIndex, only : elementSymbol
    
    implicit none
    include 'mpif.h'

    character(len=*), intent(in)  :: tplfile
    integer, intent(in) :: nqmatoms
    _REAL_,  intent(in) :: qmcoords(3,nqmatoms) 
    integer, intent(in) :: qmtypes(nqmatoms)
    integer, intent(in) :: nclatoms
    _REAL_,  intent(in) :: clcoords(4,nqmatoms)
    type(tc_nml_type), intent(in) :: tc_nml
    _REAL_, intent(out) :: escf
    _REAL_, intent(out) :: dxyzqm(3,nqmatoms)
    _REAL_, intent(out) :: dxyzcl(3,nclatoms)
    _REAL_, intent(out) :: dipmom(4,3)
    _REAL_, intent(out) :: qmcharges(nqmatoms)
    logical, intent(in) :: do_grad
    character(len=3), intent(in) :: id
    integer         , intent(in) :: charge, spinmult

    character(len=2)    :: atom_types(nqmatoms)
    _REAL_              :: coords(3,nqmatoms+nclatoms)
    _REAL_              :: charges(nclatoms)
    _REAL_              :: dxyz_all(3,nclatoms+nqmatoms)

    logical,save        :: first_call=.true.
    integer             :: i, status(MPI_STATUS_SIZE)
    integer             :: ierr

    call debug_enter_function( 'mpi_hook', module_name, tc_nml%debug )

    ! Determine atom types
    ! TeraChem needs those both for initialization and later during the MD run
    do i = 1, nqmatoms
      atom_types(i)=elementSymbol(qmtypes(i))
   end do

    ! ---------------------------------------------------
    ! Initialization: Connect to "terachem_port", set    
    ! newcomm (global), send relevant namelist variables.
    ! ---------------------------------------------------
    if (first_call) then 
      first_call=.false.
      call connect_to_terachem( tplfile, tc_nml, nqmatoms, atom_types, do_grad, id, charge, spinmult )
    end if

    ! -----------------------------------------
    ! Begin sending data each step to terachem
    ! -----------------------------------------
    if ( tc_nml%debug > 1 ) then
       write(6,'(a)') 'Sending data to TeraChem'
       call flush(6)
    end if

    ! Send nqmatoms and the type of each qmatom
    if ( tc_nml%debug > 2 ) then
       write(6,'(/, a, i0)') 'Sending nqmatoms = ', nqmatoms
       call flush(6)
    end if
    call MPI_Send( nqmatoms, 1, MPI_INTEGER, 0, 2, newcomm, ierr )

    if ( tc_nml%debug > 2 ) then
       write(6,'(/,a)') 'Sending QM atom types: '
       do i = 1, nqmatoms
          write(6,'(a)') atom_types(i)
          call flush(6)
       end do
    end if
    call MPI_Send( atom_types, 2*size(atom_types), MPI_CHARACTER, 0, 2, newcomm, ierr )

    ! Send QM coordinate array
    if ( tc_nml%debug > 2 ) then
       write(6,'(a)') 'Sending QM coords: '
    end if
    do i=1, nqmatoms
       if ( tc_nml%debug > 2 ) then
          write(6,*) 'Atom ',i,': ',qmcoords(:,i)
          call flush(6)
       end if
    end do 
    call MPI_Send( qmcoords, 3*nqmatoms, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr ) 

    ! Send nclatoms and the charge of each atom
    if ( tc_nml%debug > 2 ) then
       write(6,'(a, i0)') 'Sending nclatoms = ', nclatoms
       call flush(6)
    end if
    call MPI_Send( nclatoms, 1, MPI_INTEGER, 0, 2, newcomm, ierr ) 

    if ( tc_nml%debug > 2 ) then
       write(6,'(a)') 'Sending charges: '
    end if
    do i=1, nclatoms
      charges(i)=clcoords(4,i)
      if ( tc_nml%debug > 2 ) then
         write(6,*) 'Charge ',i,':',charges(i)
         call flush(6)
      end if
    end do
    call MPI_Send( charges, nclatoms, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr ) 

    ! Send MM point charge coordinate array
    if ( tc_nml%debug > 2 ) then
       write(6,'(a)') 'Sending CL coords: '
    end if
    do i=1, nclatoms
      coords(:,i)=clcoords(:3,i)
      if ( tc_nml%debug > 2 ) then
         write(6,*) 'Atom ',i,': ',coords(:,i)
         call flush(6)
      end if
    end do 
    call MPI_Send( coords, 3*nclatoms, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr ) 

    ! -----------------------------------
    ! Begin receiving data from terachem
    ! -----------------------------------

    ! Energy
    if ( tc_nml%debug > 2 ) then
       write(6,'(a)') 'Waiting to receive scf energy from TeraChem...'
       call flush(6)
    end if
    call MPI_Recv( escf, 1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr )
    if ( tc_nml%debug > 1 ) then
       write(6,'(a,es15.6)') 'Received scf energy from server:', escf
       call flush(6)
    end if

    ! Charges (Mulliken or other)
    if ( tc_nml%debug > 2 ) then
       write(6,'(a)') 'Waiting to receive charges...'
    end if
    call MPI_Recv( qmcharges(:), nqmatoms, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr )
    if ( tc_nml%debug > 2 ) then
       write(6,'(a)') 'Received the following charges from server:'
       do i=1, nqmatoms
          write(6,*) 'Atom ',i, ': ', qmcharges(i)
       end do
       call flush(6)
    end if

    ! Dipole moment
    if ( tc_nml%debug > 2 ) then
       write(6,'(a)') 'Waiting to receive dipole moment...'
    end if
    ! QM dipole moment
    call MPI_Recv( dipmom(:,1), 4, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr )
    if ( tc_nml%debug > 1 ) then
       write(6,'(a,4es15.6)') 'Received QM  dipole moment from server:', dipmom(:,1)
       call flush(6)
    end if
    ! MM dipole moment
    call MPI_Recv( dipmom(:,2), 4, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr )
    if ( tc_nml%debug > 1 ) then
       write(6,'(a,4es15.6)') 'Received MM  dipole moment from server:', dipmom(:,2)
       call flush(6)
    end if
    ! TOT dipole moment
    call MPI_Recv( dipmom(:,3), 4, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr )
    if ( tc_nml%debug > 1 ) then
       write(6,'(a,4es15.6)') 'Received TOT dipole moment from server:', dipmom(:,3)
       call flush(6)
    end if
    
    ! QM gradients
    if ( do_grad ) then
       if ( tc_nml%debug > 2 ) then
          write(6,'(a)') 'Waiting to receive gradients...'
       end if
       call MPI_Recv( dxyz_all, 3*(nqmatoms+nclatoms), MPI_DOUBLE_PRECISION, &
            MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr )
       if ( tc_nml%debug > 2 ) then
          write(6,'(a)') 'Received the following gradients from server:'
          do i=1, nqmatoms+nclatoms
             write(6,*) 'Atom ',i, ': ',dxyz_all(:,i)
          end do
          call flush(6)
       end if
       
       ! Poplulate our output arrays with gradients from terachem
       do i=1, nqmatoms
          dxyzqm(:,i)=dxyz_all(:,i)
       end do
       do i=1, nclatoms
          dxyzcl(:,i)=dxyz_all(:,i+nqmatoms)
       end do

    end if

    call debug_exit_function( 'connect_to_terachem', module_name, tc_nml%debug )

  end subroutine mpi_hook

  ! -------------------------------------------------
  ! Search for name published by TeraChem and connect
  ! (this step initializes newcomm)
  ! Send relevant namelist variables to terachem
  ! -------------------------------------------------
  subroutine connect_to_terachem( tplfile, tc_nml, nqmatoms, atom_types, do_grad, id, charge, spinmult )

    implicit none
    include 'mpif.h'

    character(len=*) , intent(in) :: tplfile
    type(tc_nml_type), intent(in) :: tc_nml
    integer          , intent(in) :: nqmatoms
    character(len=2) , intent(in) :: atom_types(nqmatoms)
    logical          , intent(in) :: do_grad
    character(len=3) , intent(in) :: id
    integer          , intent(in) :: charge, spinmult

    character(len=17) :: server_name="terachem_port"
    integer, parameter  :: clen=128 ! Length of character strings we are using
    character(255) :: port_name
    character(len=clen) :: dbuffer(2,32)
    _REAL_          :: timer
    integer         :: ierr, i, j, irow
    logical         :: done=.false.

    call debug_enter_function( 'connect_to_terachem', module_name, tc_nml%debug )

    ! -----------------------------------
    ! Look for server_name, get port name
    ! After 60 seconds, exit if not found
    ! -----------------------------------
    if ( trim(id) /= '' ) then
      server_name = trim(server_name)//'.'//trim(id)
    end if
    if ( tc_nml%debug > 1 ) then
      write(6,'(2a)') 'Looking up server under name:', trim(server_name)
      call flush(6)
    end if
    timer = MPI_WTIME(ierr)
    do while (done .eqv. .false.)

      call MPI_LOOKUP_NAME(trim(server_name), MPI_INFO_NULL, port_name, ierr)
      if (ierr == MPI_SUCCESS) then
        if ( tc_nml%debug > 1 ) then
          write(6,'(2a)') 'Found port: ', trim(port_name)
          call flush(6)
        end if
        done=.true.

      end if

      if ( (MPI_WTIME(ierr)-timer) > 60 ) then ! Time out after 60 seconds
        call sander_bomb('connect_to_terachem() ('//module_name//')', &
          '"'//trim(server_name)//'" not found. Timed out after 60 seconds.', &
          'Will quit now')
      end if

    end do

    ! ----------------------------------------
    ! Establish new communicator via port name
    ! ----------------------------------------
    call MPI_COMM_CONNECT(port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF, newcomm, ierr)
    if ( tc_nml%debug > 1 ) then
      write(6,'(a,i0)') 'Established new communicator:', newcomm
      call flush(6)
    end if

    ! --------------------------------
    ! Send job information to terachem
    ! --------------------------------
    dbuffer(:,:) = ''
    if ( tc_nml%use_template ) then
       ! using template input file for specs of QM method
      call read_template(tplfile, tc_nml%debug, dbuffer, irow )
    else
      ! specs of QM method from AMBER input file 
      irow = 1
      write(dbuffer(:,irow),'(a,/,a)') 'basis',      tc_nml%basis
      irow = irow + 1
      write(dbuffer(:,irow),'(a,/,a)') 'method',     tc_nml%method
      irow = irow + 1
      write(dbuffer(:,irow),'(a,/,a)') 'dftd',       tc_nml%dftd
      irow = irow + 1
      write(dbuffer(:,irow),'(a,/,a)') 'precision',  tc_nml%precision
      irow = irow + 1
      write(dbuffer(:,irow),'(a,/,E22.16)') 'threall', tc_nml%threall
      irow = irow + 1
      write(dbuffer(:,irow),'(a,/,E22.16)') 'convthre', tc_nml%convthre
      irow = irow + 1
      write(dbuffer(:,irow),'(a,/,i0)') 'maxit', tc_nml%maxit
      irow = irow + 1
      write(dbuffer(:,irow),'(a,/,i0)') 'dftgrid', tc_nml%dftgrid
      irow = irow + 1
      write(dbuffer(:,irow),'(a,/,a)') 'cis',        tc_nml%cis
      if ( trim(tc_nml%cis) == 'yes' ) then
        irow = irow + 1
        write(dbuffer(:,irow),'(a,/,i0)') 'cisnumstates', tc_nml%cisnumstates
        irow = irow + 1
        write(dbuffer(:,irow),'(a,/,i0)') 'cistarget', tc_nml%cistarget
      end if
    end if

    ! common data provided both for template / AMBER input
    irow = irow + 1
    write(dbuffer(:,irow),'(a,/,i0)') 'charge', charge
    irow = irow + 1
    write(dbuffer(:,irow),'(a,/,i0)') 'spinmult', spinmult
    irow = irow + 1
    if ( do_grad ) then
       write(dbuffer(:,irow),'(a,/,a)') 'run', 'gradient'
    else
       write(dbuffer(:,irow),'(a,/,a)') 'run', 'energy'
    end if
    irow = irow + 1
    write(dbuffer(:,irow),'(a,/,a)') 'guess',       tc_nml%guess
    ! This is set to 1 because this instructs terachem to skip 
    ! calculating the self-energy of the charges
    irow = irow + 1
    write(dbuffer(:,irow),'(a,/,a)') 'amber', 'yes'
    ! Write gpus
    if ( tc_nml%ngpus > 0 ) then
       irow = irow + 1
       write(dbuffer(:,irow), '(a,/,9999(i3))') 'gpus      ', tc_nml%ngpus, (tc_nml%gpuids(i), i = 1, tc_nml%ngpus)
    end if
    ! Finish writing - send 'end'
    irow = irow + 1
    write(dbuffer(:,irow),'(a)') 'end', ''

    if ( tc_nml%debug > 2 ) then
      write(6,'(a)') '(debug) sending namelist data:'
      do j=1, 32
        write(6,*) trim(dbuffer(1,j)), ' = ', trim(dbuffer(2,j))
      end do
      call flush(6)
    end if

    call MPI_Send( dbuffer, 2*clen*size(dbuffer,2), MPI_CHARACTER, 0, 2, newcomm, ierr )

    ! -----------------------------------------
    ! Send nqmatoms and the type of each qmatom
    ! TeraChem needs this information to correctly initialize the GPUs
    ! (depending on whether d orbitals are in use or not)
    ! -----------------------------------------
    if ( tc_nml%debug > 2 ) then
      write(6,'(/, a, i0)') 'Sending nqmatoms = ', nqmatoms
      call flush(6)
    end if
    call MPI_Send( nqmatoms, 1, MPI_INTEGER, 0, 2, newcomm, ierr )

    if ( tc_nml%debug > 2 ) then
      write(6,'(/,a)') 'Sending QM atom types: '
      do i = 1, nqmatoms
        write(6,'(a)') atom_types(i)
        call flush(6)
      end do
    end if
    call MPI_Send( atom_types, 2*size(atom_types), MPI_CHARACTER, 0, 2, newcomm, ierr )


    call debug_exit_function( 'connect_to_terachem', module_name, tc_nml%debug )

  end subroutine connect_to_terachem

#endif

  ! ---------------------------------
  ! Write the input file for TeraChem
  ! ---------------------------------
  subroutine write_inpfile( inpfile, crdfile, ptcfile, tplfile, &
       nqmatoms, qmcoords, qmtypes, nclatoms, clcoords, &
       tc_nml, do_grad, charge, spinmult )

    use ElementOrbitalIndex, only : elementSymbol
    use qm2_extern_util_module, only: write_chgfile

    implicit none

    character (len=*) , intent(in) :: inpfile
    character (len=*) , intent(in) :: crdfile
    character (len=*) , intent(in) :: ptcfile
    character (len=*) , intent(in) :: tplfile
    integer           , intent(in) :: nqmatoms
    _REAL_            , intent(in) :: qmcoords(:,:)
    integer           , intent(in) :: qmtypes(:)
    integer           , intent(in) :: nclatoms
    _REAL_            , intent(in) :: clcoords(:,:)
    type(tc_nml_type) , intent(in) :: tc_nml
    logical, intent(in) :: do_grad
    integer           , intent(in) :: charge, spinmult

    character(len=20)  :: bakfile
    integer            :: i, ios
    integer, parameter :: iurun = 10
    logical, save :: first_call = .true.

    call debug_enter_function( 'write_inpfile', module_name, tc_nml%debug )

    if ( tc_nml%use_template ) then
      bakfile = tplfile//'.bak'
      if ( first_call ) then
        call copy_template(tplfile, trim(bakfile), tc_nml%debug)
      end if
      call system('cp '//trim(bakfile)//' '//inpfile)
    end if

    open(iurun, file=inpfile, position='append', iostat=ios)
    if ( ios > 0 ) then
      call sander_bomb('write_inpfile (qm2_extern_tc_module)', &
           'Error opening TeraChem input file '//inpfile//' for writing', &
           'Will quit now')
    end if

    ! Following should be in template input file, if used
    if ( .not. tc_nml%use_template ) then
      write(iurun, '(a,/,a)')'# Run using SANDER file-based interface for TeraChem','#'
      write(iurun, '(2a)')       'basis        ', trim(tc_nml%basis)
      write(iurun, '(2a)')       'method       ', trim(tc_nml%method)
      write(iurun, '(2a)')       'precision    ', trim(tc_nml%precision)
      write(iurun, '(a,E22.16)') 'threall      ', tc_nml%threall
      write(iurun, '(a,E22.16)') 'convthre     ', tc_nml%convthre
      write(iurun, '(2a)')       'dftd         ', trim(tc_nml%dftd)
      write(iurun, '(a,i4)')     'maxit        ', tc_nml%maxit
      write(iurun, '(a,i3)')     'dftgrid      ', tc_nml%dftgrid
      write(iurun, '(2a)')       'cis          ', trim(tc_nml%cis)
      if ( trim(tc_nml%cis) == 'yes' ) then
        write(iurun, '(a,i3)')     'cisnumstates ', tc_nml%cisnumstates
        write(iurun, '(a,i3)')     'cistarget    ', tc_nml%cistarget
      end if
    end if

    ! common to templete / input via AMBER &tc namelist
    write(iurun, '(a,i3)')     'charge       ', charge
    write(iurun, '(a,i3)')     'spinmult     ', spinmult
    if ( tc_nml%ngpus > 0 ) then
       write(iurun, '(a,65(i3))')    'gpus      ', tc_nml%ngpus, (tc_nml%gpuids(i), i = 1, tc_nml%ngpus)
    end if

    if ( .not. first_call ) then
      write(iurun, '(2a)') 'guess       ', tc_nml%guess
    end if

    if ( do_grad ) then
      write(iurun, '(a)') 'run         gradient'
    else
      write(iurun, '(a)') 'run         energy'
    end if
    write(iurun, '(a)') 'coordinates '//crdfile

    if ( nclatoms > 0 ) then
      write(iurun, '(a)') 'pointcharges '//ptcfile
      write(iurun, '(a)') 'amber       yes'
    end if

    write(iurun, '(a)') 'end'

    close(iurun)
    
    ! Now write xyz file with geometry
    open(iurun, file=crdfile, iostat=ios)
    if ( ios > 0 ) then
      call sander_bomb('write_inpfile (qm2_extern_tc_module)', &
        'Error opening coordinate file '//crdfile//' for writing', &
        'Will quit now')
    end if
    write(iurun,'(i5,/)') nqmatoms
    do i = 1, nqmatoms
      write(iurun,'(a2,1x,3f21.16)') &
        elementSymbol(qmtypes(i)), &
        qmcoords(:,i)
    enddo
    close(iurun)

    ! Now write xyz file with point charges
    if ( nclatoms > 0 ) then
      call write_chgfile( ptcfile, nclatoms, clcoords, tc_nml%debug )
    end if
    
    call debug_exit_function( 'write_inpfile', module_name, tc_nml%debug )

    if (first_call) then
      first_call = .false.
    end if

  end subroutine write_inpfile

  ! ---------------------
  ! Read TeraChem results
  ! ---------------------
  subroutine read_results( datfile, nqmatoms, nclatoms, escf, dxyzqm, dxyzcl,&
       dipmom, qmcharges, charge_analysis, do_grad, debug )

    implicit none

    character(len=*), intent(in) :: datfile
    integer, intent(in) :: nqmatoms, nclatoms
    _REAL_, intent(out) :: escf, dxyzqm(3,nqmatoms), dxyzcl(3,nclatoms)
    logical, intent(in) :: do_grad
    integer, intent(in) :: debug

    integer :: ios, i, itmp
    integer, parameter :: iunit = 351
    _REAL_, intent(out) :: dipmom(4,3)
    _REAL_, intent(out) :: qmcharges(nqmatoms)
    character(len=*), intent(in) :: charge_analysis
    character(len=256) :: read_buffer

    character(len=19), parameter :: mulliken_charge_file = 'scr/charge_mull.xls'
    character, parameter :: tab = achar(9)

    call debug_enter_function( 'read_results', module_name, debug )

    open(iunit, file=datfile, iostat=ios)
    if ( ios /= 0 ) then
       call sander_bomb('read_results (qm2_extern_tc_module)', &
            'Error opening TeraChem output file '//datfile//' for reading', &
            'Will quit now')
    end if

    do
      read (iunit, '(a)', iostat = ios) read_buffer
      ! End of file; data not found
      if (ios < 0 ) then
        call sander_bomb('read_results (qm2_extern_tc_module)', &
          'Error reading TeraChem output from file '//datfile, &
          '("FINAL ENERGY" or "dE/dX" or "MM / Point charge part" not found)')
      end if

      if ( read_buffer(1:14) == 'FINAL ENERGY: ' ) then
        ! Read energy
        itmp = len_trim(read_buffer)
        read (read_buffer(15:itmp-5),*) escf
        if ( .not. do_grad ) exit
      else if ( read_buffer(9:13) == 'dE/dX' ) then
        ! Read QM gradient data
        do i = 1, nqmatoms
           read (iunit, '(3(f16.10,1x))') dxyzqm(:,i)
        end do
        if ( nclatoms == 0 ) exit
      else if ( read_buffer(9:31) == 'MM / Point charge part' ) then
        ! Read MM atoms gradient data
        do i = 1, nclatoms
          read (iunit, '(3(f16.10,1x))') dxyzcl(:,i)
        end do
        exit
      end if

      if (index(read_buffer,'QM  DIPOLE MOMENT: ') > 0 ) then
        read(read_buffer(index(read_buffer,'{')+1:index(read_buffer,'}')-1),*) dipmom(1:3,1)
        read(read_buffer(index(read_buffer,'|D| = ')+5:index(read_buffer,') DEBYE')-1),*) dipmom(4,1)
      end if
      if (index(read_buffer,'MM  DIPOLE MOMENT: ') > 0 ) then
        read(read_buffer(index(read_buffer,'{')+1:index(read_buffer,'}')-1),*) dipmom(1:3,2)
        read(read_buffer(index(read_buffer,'|D| = ')+5:index(read_buffer,') DEBYE')-1),*) dipmom(4,2)
      end if
      if (index(read_buffer,'TOT DIPOLE MOMENT: ') > 0 ) then
        read(read_buffer(index(read_buffer,'{')+1:index(read_buffer,'}')-1),*) dipmom(1:3,3)
        read(read_buffer(index(read_buffer,'|D| = ')+5:index(read_buffer,') DEBYE')-1),*) dipmom(4,3)
      end if

      if (index(read_buffer(1:15),'DIPOLE MOMENT: ') > 0 ) then
        read(read_buffer(index(read_buffer,'{')+1:index(read_buffer,'}')-1),*) dipmom(1:3,1)
        read(read_buffer(index(read_buffer,'|D| = ')+5:index(read_buffer,') DEBYE')-1),*) dipmom(4,1)
      end if

    end do

    close(iunit)

    ! Read Mulliken charges if requested
    if ( trim(charge_analysis) == 'MULLIKEN' ) then

       open(iunit, file=mulliken_charge_file, iostat=ios)
       if ( ios /= 0 ) then
          call sander_bomb('read_results (qm2_extern_tc_module)', &
               'Error opening TeraChem Mulliken charge file '//mulliken_charge_file//' for reading', &
               'Will quit now')
       end if

       do i = 1, nqmatoms

          read (iunit, '(a)', iostat=ios) read_buffer
          ! End of file; data not found
          if (ios < 0 ) then
             call sander_bomb('read_results (qm2_extern_tc_module)', &
                  'Error reading TeraChem output from Mulliken charge file '//mulliken_charge_file, &
                  '(not enough entries)')
          end if

          itmp = index (read_buffer, tab)
          read (read_buffer(itmp:), *) qmcharges(i)

       end do
       
       close(iunit)
       
    end if

    call debug_exit_function( 'read_results', module_name, debug )

  end subroutine read_results


  subroutine tc_finalize()

    implicit none
#ifdef MPI
    include 'mpif.h'

    integer :: ierr
    _REAL_  :: empty
    if (do_mpi) then

      call MPI_Send( empty, 1, MPI_DOUBLE_PRECISION, 0, 0, newcomm, ierr )

    end if
#endif

  end subroutine tc_finalize
      

  subroutine copy_template( tplfile, bakfile, debug )

    use UtilitiesModule, only: Upcase
    implicit none
    character(len=*), intent(in) :: tplfile, bakfile
    integer         , intent(in) :: debug

    integer, parameter :: tplunit = 351, bakunit=352
    character(len=256) :: read_buffer
    integer :: tplerr, bakerr, ios

    call debug_enter_function( 'copy_template', module_name, debug )

    open(tplunit, file=tplfile, iostat=tplerr )
    open(bakunit, file=bakfile, iostat=bakerr )
    if ( tplerr /= 0 ) then
      call sander_bomb('copy_template (qm2_extern_tc_module)', &
        'Error opening ADF template file '//tplfile//' for reading', &
        'Will quit now.')
    end if
    if ( bakerr /= 0 ) then
      call sander_bomb('copy_template (qm2_extern_tc_module)', &
        'Error opening ADF template backup file '//bakfile//' for writing', &
        'Will quit now.')
    end if

    ! Write tplfile (without end key) to bakfile
    do
       read (tplunit, '(a)', iostat = ios) read_buffer
       ! End of file; stop writing
       if (ios < 0 ) then
         exit
       end if
       if ( index(Upcase(read_buffer), 'END') > 0 ) then
          exit
       end if
       write(bakunit, '(a)') trim(read_buffer)
    end do

    close(tplunit)
    close(bakunit)

    call debug_exit_function( 'copy_template', module_name, debug )

  end subroutine copy_template

  subroutine read_template( tplfile, debug, dbuffer, irow )

    use UtilitiesModule, only: Upcase
    implicit none
    character(len=*), intent(in) :: tplfile
    integer         , intent(in) :: debug
    character(len=*), intent(out):: dbuffer(:,:)
    integer         , intent(out):: irow
    integer, parameter :: tplunit = 351
    character(len=256) :: read_buffer
    integer :: tplerr, ios, ispace

    call debug_enter_function( 'read_template', module_name, debug )

    dbuffer=''

    open(tplunit, file=tplfile, iostat=tplerr )
    if ( tplerr /= 0 ) then
      call sander_bomb('read_template (qm2_extern_tc_module)', &
        'Error opening TeraChem template file '//tplfile//' for reading', &
        'Will quit now.')
    end if

    irow = 0
    do
       read (tplunit, '(a)', iostat = ios) read_buffer
       ! End of file; stop reading
       if (ios < 0 ) then
         exit
       end if
       if ( index(Upcase(read_buffer), 'END') > 0 ) then
          exit
       end if
       ispace=index(read_buffer,' ')
       irow = irow + 1
       write(dbuffer(:,irow), '(a,/,a)') read_buffer(1:ispace-1), trim(adjustl(read_buffer(ispace:)))
    end do

    close(tplunit)

    call debug_exit_function( 'read_template', module_name, debug )

  end subroutine read_template

end module qm2_extern_tc_module
