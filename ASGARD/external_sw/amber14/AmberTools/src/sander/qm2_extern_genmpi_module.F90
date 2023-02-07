#include "../include/dprec.fh"
module qm2_extern_genmpi_module
! ----------------------------------------------------------------
! Interface for QM and QM/MM MD via MPI interface
!
! Currently supports:
! pure QM
! QM/MM with electronic embedding of point charges
!
! Author: Andreas Goetz (agoetz@sdsc.edu)
!
! Date: August 2013
!
! ----------------------------------------------------------------

  use qm2_extern_util_module, only: debug_enter_function, debug_exit_function

  implicit none

  private
#ifdef MPI
  public :: get_genmpi_forces, genmpi_finalize

  character(len=*), parameter :: module_name = "qm2_extern_mpi_module"

  type genmpi_nml_type
     character(len=20) :: method
     character(len=20) :: basis
     character(len=20) :: jbasis
     character(len=20) :: cbasis
     character(len=20) :: scfconv
     character(len=20) :: scfiter
     character(len=20) :: guess
     character(len=20) :: grid
     integer           :: ntpr
     integer           :: debug
     logical           :: dipole
     logical           :: use_template
  end type genmpi_nml_type

  integer, save         :: newcomm ! Initialized in mpi_init subroutine

contains

  ! ------------------------------------------
  ! Get QM energy and forces via MPI interface
  ! ------------------------------------------
  subroutine get_genmpi_forces( do_grad, nstep, ntpr_default, id, nqmatoms, qmcoords,&
       qmtypes, nclatoms, clcoords, escf, dxyzqm, dxyzcl, charge, spinmult )

    use qm2_extern_util_module, only: print_results, write_dipole
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

    _REAL_              :: dipmom(4)          ! Dipole moment {x, y, z, |D|}

    type(genmpi_nml_type), save :: self
    logical, save :: first_call = .true.
    integer :: i
    integer :: printed =-1 ! Used to tell if we have printed this step yet 
                           ! since the same step may be called multiple times
    
    character(len=*),  parameter :: basename = 'mpi_job'
    character(len=*),  parameter :: dipext = '.dip'
    character(len=*),  parameter :: tplext = '.tpl'
    character(len=14) :: dipfile, tplfile
    
    ! assemble output data filenames
    dipfile = basename//trim(id)//dipext 
    tplfile = basename//tplext 

    ! Setup on first program call
    if ( first_call ) then
      first_call = .false.
      write (6,'(/,a,/)') '   >>> Running QM calculation with generic MPI interface <<<'
      call get_namelist( ntpr_default, self )
      call print_namelist(self)
      write (6,'(80a)') ('-', i=1,80)
      write (6,'(a)') '   4.  RESULTS'
      write (6,'(80a)') ('-', i=1,80)
      call system('rm -f '//dipfile)
    end if


    call mpi_hook( trim(tplfile), nqmatoms, qmcoords, qmtypes, nclatoms, clcoords,&
         self, escf, dxyzqm, dxyzcl, dipmom, do_grad, id, charge, spinmult )


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

    call print_results( module_name, escf, nqmatoms, dxyzqm,&
         self%debug, nclatoms, dxyzcl )

    ! Write dipole moment to file
    if ( self%ntpr > 0 .and. mod(nstep, self%ntpr) == 0 ) then
       if ( printed /= nstep .and. self%dipole ) then
          ! using util module's write dipole, writing only QM dipole moment
          call write_dipole( trim(dipfile), dipmom(1:3), dipmom(4), self%debug )
          printed = nstep
       end if
    end if

  end subroutine get_genmpi_forces

  ! ---------------------------------------
  ! Read mpi namelist values from file mdin
  ! use default values if none are present
  ! ---------------------------------------
  subroutine get_namelist( ntpr_default, self )

    use UtilitiesModule, only: Upcase
    implicit none

    type(genmpi_nml_type), intent(out) :: self
    integer, intent(in) :: ntpr_default
    character(len=20):: method, basis, jbasis, cbasis, scfconv, scfiter, guess, grid
    integer :: ntpr, debug, dipole, use_template
    namelist /genmpi/ method, basis, jbasis, cbasis, scfconv, scfiter, guess, grid, &
         use_template, ntpr, debug, dipole

    integer :: ierr

    ! Default values
    method   = 'blyp'
    basis    = '6-31g*'
    jbasis   = 'none'
    cbasis   = 'none'
    scfconv  = '1E-08'
    scfiter  = '100'
    guess    = 'read'
    grid     = 'none'
    ntpr     = ntpr_default
    debug    = 0
    dipole   = 0
    use_template = 0

    ! Read namelist
    rewind 5
    read(5,nml=genmpi,iostat=ierr)

    if ( ierr > 0 ) then
       call sander_bomb('get_namelist ('//module_name//')', &
            '&mpi namelist read error', &
            'Please check your input.')
    else if ( ierr < 0 ) then
       write(6,'(a/a)') '&mpi namelist read encountered end of file', &
            'Please check your input if the calculation encounters a problem'
    end if

    ! Assign namelist values to self data type

    self%method   = method  
    self%basis    = basis   
    self%jbasis   = jbasis  
    self%cbasis   = cbasis  
    self%scfconv  = scfconv 
    self%scfiter  = scfiter 
    self%guess    = guess   
    self%grid     = grid    
    self%ntpr     = ntpr
    self%debug    = debug

    if ( dipole == 0 ) then
       self%dipole = .false.
    else if ( dipole == 1 ) then
       self%dipole = .true.
    else
       call sander_bomb('get_namelist ('//module_name//')', &
            '&mpi dipole value not allowed', &
            'Please check your input. dipole can only be 0 or 1.')
    end if

    if ( use_template == 0 ) then
       self%use_template = .false.
    else if ( use_template == 1 ) then
       self%use_template = .true.
    else
       call sander_bomb('get_namelist ('//module_name//')', &
            '&mpi use_template value not allowed', &
            'Please check your input. use_template can only be 0 or 1.')
    end if


  end subroutine get_namelist

  ! -----------------------
  ! Print namelist settings
  ! -----------------------
  subroutine print_namelist( self )

    implicit none
    type(genmpi_nml_type), intent(in) :: self

    write(6, '(/,a)')      '| &mpi'
    write(6, '(2a)')       '|   method       = ', self%method
    write(6, '(2a)')       '|   basis        = ', self%basis
    write(6, '(2a)')       '|   jbasis       = ', self%jbasis
    write(6, '(2a)')       '|   cbasis       = ', self%cbasis
    write(6, '(2a)')       '|   scfconv      = ', self%scfconv
    write(6, '(2a)')       '|   scfiter      = ', self%scfiter
    write(6, '(2a)')       '|   guess        = ', self%guess
    write(6, '(2a)')       '|   grid         = ', self%grid
    write(6, '(a,i0)')     '|   ntpr         = ', self%ntpr
    write(6, '(a,i2)')     '|   debug        = ', self%debug
    write(6, '(a,l)')      '|   dipole       = ', self%dipole
    write(6, '(a,l)')      '|   use_template = ', self%use_template
    write(6,'(a)')         '| /'

  end subroutine print_namelist


#  ifndef MPI_1
  ! Perform MPI communications. Requires MPI 2.0 or above to use
  subroutine mpi_hook( tplfile, nqmatoms, qmcoords, qmtypes, nclatoms, clcoords,&
       self, escf, dxyzqm, dxyzcl, dipmom, do_grad, id, charge, spinmult )
    
    use ElementOrbitalIndex, only : elementSymbol
    
    implicit none
    include 'mpif.h'

    character(len=*), intent(in)  :: tplfile
    integer, intent(in) :: nqmatoms
    _REAL_,  intent(in) :: qmcoords(3,nqmatoms) 
    integer, intent(in) :: qmtypes(nqmatoms)
    integer, intent(in) :: nclatoms
    _REAL_,  intent(in) :: clcoords(4,nqmatoms)
    type(genmpi_nml_type), intent(in) :: self
    _REAL_, intent(out) :: escf
    _REAL_, intent(out) :: dxyzqm(3,nqmatoms)
    _REAL_, intent(out) :: dxyzcl(3,nclatoms)
    _REAL_, intent(out) :: dipmom(4)
    logical, intent(in) :: do_grad
    character(len=3), intent(in) :: id
    integer         , intent(in) :: charge, spinmult

    character(len=2)    :: atom_types(nqmatoms)
    _REAL_              :: coords(3,nclatoms)
    _REAL_              :: charges(nclatoms)

    logical,save        :: first_call=.true.
    integer             :: i, status(MPI_STATUS_SIZE)
    integer             :: ierr

    ! AWG FIXME: write charges if requested
    ! for now just receive them locally to make the interface work
    _REAL_ :: qmcharges(nqmatoms)

    call debug_enter_function( 'mpi_hook', module_name, self%debug )

    ! Determine atom types since we are sending element symbols, not nuclear charges
    do i = 1, nqmatoms
       atom_types(i) = elementSymbol(qmtypes(i))
    end do
    
    ! Determine classical charges and coords to be sent
    do i = 1, nclatoms
       charges(i)  = clcoords(4,i)
       coords(:,i) = clcoords(:3,i)
    end do

    ! ---------------------------------------------------
    ! Initialization: Connect to "qc_program_port", set    
    ! newcomm (global), send relevant namelist variables.
    ! ---------------------------------------------------
    if (first_call) then 
      first_call=.false.
      call connect( self, id )
      call send_job_info( tplfile, self, do_grad )
    end if


    ! ----------------------------
    ! Begin sending data each step
    ! ----------------------------
    if ( self%debug > 1 ) then
       write(6,'(a)') 'Sending data to QC program'
       call flush(6)
    end if

    ! Send charge and spin
    if ( self%debug > 2 ) then
       write(6,'(/a,i0, /a,i0)') &
            'Sending charge   = ', charge, &
            'Sending spinmult = ', spinmult
       call flush(6)
    end if
    call MPI_Send( charge,   1, MPI_INTEGER, 0, 1, newcomm, ierr )
    call MPI_Send( spinmult, 1, MPI_INTEGER, 0, 1, newcomm, ierr )

    ! Send nqmatoms and the type of each qmatom
    if ( self%debug > 2 ) then
       write(6,'(/, a, i0)') 'Sending nqmatoms = ', nqmatoms
       call flush(6)
    end if
    call MPI_Send( nqmatoms, 1, MPI_INTEGER, 0, 1, newcomm, ierr )

    if ( self%debug > 2 ) then
       write(6,'(/,a)') 'Sending QM atom types: '
       do i = 1, nqmatoms
          write(6,'(a)') atom_types(i)
          call flush(6)
       end do
    end if
    call MPI_Send( atom_types, 2*size(atom_types), MPI_CHARACTER, 0, 1, newcomm, ierr )

    ! Send QM coordinate array
    if ( self%debug > 2 ) then
       write(6,'(a)') 'Sending QM coords: '
       do i=1, nqmatoms
          write(6,*) 'Atom ',i,': ',qmcoords(:,i)
          call flush(6)
       end do
    end if
    call MPI_Send( qmcoords, 3*nqmatoms, MPI_DOUBLE_PRECISION, 0, 1, newcomm, ierr ) 

    ! Send nclatoms and the charge of each atom
    if ( self%debug > 2 ) then
       write(6,'(a, i0)') 'Sending nclatoms = ', nclatoms
       call flush(6)
    end if
    call MPI_Send( nclatoms, 1, MPI_INTEGER, 0, 1, newcomm, ierr ) 

    if ( self%debug > 2 ) then
       write(6,'(a)') 'Sending charges: '
       do i=1, nclatoms
          write(6,*) 'Charge ',i,':',charges(i)
          call flush(6)
       end do
    end if
    call MPI_Send( charges, nclatoms, MPI_DOUBLE_PRECISION, 0, 1, newcomm, ierr ) 

    ! Send MM point charge coordinate array
    if ( self%debug > 2 ) then
       write(6,'(a)') 'Sending charge coords: '
       do i=1, nclatoms
          write(6,*) 'Charge ',i,': ',coords(:,i)
         call flush(6)
      end do
    end if
    call MPI_Send( coords, 3*nclatoms, MPI_DOUBLE_PRECISION, 0, 1, newcomm, ierr ) 

    ! ------------------------------------
    ! Begin receiving data from QC program
    ! ------------------------------------

    ! Energy
    if ( self%debug > 2 ) then
       write(6,'(a)') 'Waiting to receive scf energy from QM program...'
       call flush(6)
    end if
    call MPI_Recv( escf, 1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr )
    if ( self%debug > 1 ) then
       write(6,'(a,es15.6)') 'Received scf energy from QM program:', escf
       call flush(6)
    end if

    ! Charges (Mulliken or other)
    call MPI_Recv( qmcharges(:), nqmatoms, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr )
    if ( self%debug > 2 ) then
       write(6,'(a)') 'Received the following charges from QM program:'
       do i=1, nqmatoms
          write(6,*) 'Atom ',i, ': ', qmcharges(i)
       end do
       call flush(6)
    end if

    ! Dipole moment
    call MPI_Recv( dipmom(:), 4, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr )
    if ( self%debug > 1 ) then
       write(6,'(a,4es15.6)') 'Received QM  dipole moment from QM program:', dipmom(:)
       call flush(6)
    end if
    
    ! QM gradients
    if ( do_grad ) then

       call MPI_Recv( dxyzqm, 3*nqmatoms, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr )
       if ( self%debug > 2 ) then
          write(6,'(a)') 'Received the following QM gradients from the QM program:'
          do i=1, nqmatoms
             write(6,*) 'Atom ',i, ': ',dxyzqm(:,i)
          end do
          call flush(6)
       end if

       call MPI_Recv( dxyzcl, 3*nclatoms, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr )
       if ( self%debug > 2 ) then
          write(6,'(a)') 'Received the following MM gradients from the QM program:'
          do i=1, nclatoms
             write(6,*) 'Atom ',i, ': ',dxyzcl(:,i)
          end do
          call flush(6)
       end if
       
    end if

    call debug_exit_function( 'mpi_hook', module_name, self%debug )

  end subroutine mpi_hook

  ! ---------------------------------------------------
  ! Search for name published by QC program and connect
  ! (this step initializes newcomm)
  ! ---------------------------------------------------
  subroutine connect( self, id )

    implicit none
    include 'mpif.h'

    type(genmpi_nml_type), intent(in) :: self
    character(len=3)     , intent(in) :: id

    character(len=19) :: server_name="qc_program_port" ! Allow connection to 999 copies
    character(255) :: port_name
    _REAL_          :: timer
    integer         :: ierr
    logical         :: done=.false.

    call debug_enter_function( 'connect', module_name, self%debug )

    ! -----------------------------------
    ! Look for server_name, get port name
    ! After 60 seconds, exit if not found
    ! -----------------------------------
    if ( trim(id) /= '' ) then
      server_name = trim(server_name)//'.'//trim(id)
    end if

    if ( self%debug > 1 ) then
      write(6,'(2a)') 'Looking up server under name:', trim(server_name)
      call flush(6)
    end if

    timer = MPI_WTIME(ierr)

    do while ( done .eqv. .false. )

      call MPI_LOOKUP_NAME(trim(server_name), MPI_INFO_NULL, port_name, ierr)

      if ( ierr == MPI_SUCCESS ) then
        if ( self%debug > 1 ) then
          write(6,'(2a)') 'Found port: ', trim(port_name)
          call flush(6)
        end if
        done=.true.
      end if

      if ( (MPI_WTIME(ierr)-timer) > 60 ) then ! Time out after 60 seconds
        call sander_bomb('connect() ('//module_name//')', &
             '"'//trim(server_name)//'" not found. Timed out after 60 seconds.', &
             'Will quit now')
      end if

    end do

    ! ----------------------------------------
    ! Establish new communicator via port name
    ! ----------------------------------------
    call MPI_COMM_CONNECT(port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF, newcomm, ierr)
    if ( self%debug > 1 ) then
      write(6,'(a,i0)') 'Established new communicator:', newcomm
      call flush(6)
    end if

    call debug_exit_function( 'connect', module_name, self%debug )

  end subroutine connect


  ! ----------------------------
  ! Send initial job information
  ! ----------------------------
  subroutine send_job_info( tplfile, self, do_grad )

    implicit none
    include 'mpif.h'

    character(len=*)     , intent(in) :: tplfile
    type(genmpi_nml_type), intent(in) :: self
    logical              , intent(in) :: do_grad

    integer, parameter  :: clen=256 ! Length of character strings we are using
    character(len=clen) :: dbuffer(128)
    integer             :: i, ierr, irow

    call debug_enter_function( 'send_job_info', module_name, self%debug )


    ! Fill bugger with data to send
    dbuffer(:) = ''

    if ( .not. self%use_template ) then

      ! Send namelist data
      irow = 1
      write(dbuffer(irow),'(a)') 'method' // self%method
      irow = irow + 1                       
      write(dbuffer(irow),'(a)') 'basis'  // self%basis
      irow = irow + 1                       
      write(dbuffer(irow),'(a)') 'jbasis' // self%jbasis
      irow = irow + 1                       
      write(dbuffer(irow),'(a)') 'cbasis' // self%cbasis
      irow = irow + 1                       
      write(dbuffer(irow),'(a)') 'scfconv'// self%scfconv
      irow = irow + 1                       
      write(dbuffer(irow),'(a)') 'scfiter'// self%scfiter
      irow = irow + 1                       
      write(dbuffer(irow),'(a)') 'guess'  // self%guess
      irow = irow + 1                       
      write(dbuffer(irow),'(a)') 'grid'   // self%grid
      irow = irow + 1
      if ( do_grad ) then
        write(dbuffer(irow),'(a)') 'gradient true'
      else
        write(dbuffer(irow),'(a)') 'gradient false'
      end if

    else

      call read_template(tplfile, self%debug, dbuffer )

    end if

    if ( self%debug > 2 ) then
      write(6,'(a)') '(debug) sending namelist data:'
      do i = 1, size(dbuffer)
        write(6,*) trim(dbuffer(i))
      end do
      call flush(6)
    end if

    ! Now send the data
    call MPI_Send( dbuffer, clen*size(dbuffer), MPI_CHARACTER, 0, 1, newcomm, ierr )

    call debug_exit_function( 'send_job_info', module_name, self%debug )

  end subroutine send_job_info

#  endif /* MPI_1 */


  ! Send a final message with tag 0 to tell the QC program that we are done
  ! (this is sent in place of qmcharge)
  subroutine genmpi_finalize()

    implicit none
    include 'mpif.h'

    integer :: ierr, empty

    call MPI_Send( empty, 1, MPI_INTEGER, 0, 0, newcomm, ierr )

  end subroutine genmpi_finalize
      
  ! Read template input file into buffer to be sent
  subroutine read_template( tplfile, debug, dbuffer )

    use UtilitiesModule, only: Upcase

    implicit none

    character(len=*), intent(in) :: tplfile
    integer         , intent(in) :: debug
    character(len=*), intent(out):: dbuffer(:)

    integer, parameter :: tplunit = 351
    character(len=256) :: read_buffer
    integer :: tplerr, ios, i

    call debug_enter_function( 'read_template', module_name, debug )

    dbuffer=''

    open( tplunit, file=tplfile, iostat=tplerr )
    if ( tplerr /= 0 ) then
      call sander_bomb('read_template() ('//module_name//')', &
           'Error opening template file '//tplfile//' for reading', &
           'Will quit now.')
    end if

    i = 1
    do
       read (tplunit, '(a)', iostat = ios) read_buffer
       ! End of file; stop reading
       if (ios < 0 ) then
         exit
       end if
       write(dbuffer(i), '(a)') read_buffer
       i=i+1
    end do

    close(tplunit)

    call debug_exit_function( 'read_template', module_name, debug )

  end subroutine read_template
#endif /* MPI */

end module qm2_extern_genmpi_module
