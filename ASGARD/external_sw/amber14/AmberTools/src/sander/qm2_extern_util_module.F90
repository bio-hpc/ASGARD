#include "../include/dprec.fh"
module qm2_extern_util_module
! ----------------------------------------------------------------
! ----------------------------------------------------------------

  implicit none

  private
  public :: print_results, check_installation, write_dipole, write_charges, &
       write_chgfile, debug_enter_function, debug_exit_function
!    au_to_kcal, a_to_bohr, set_zero
  
  character(len=*), parameter :: module_name = "qm2_extern_util_module"

  contains

  ! -----------------------------------------------------------------------
  ! Print out escf and QM + MM forces at the end of get_forces
  ! Note that specifiying MM forces is optional
  ! -----------------------------------------------------------------------
  subroutine print_results( extern_routine, escf, nqmatoms, dxyzqm, debug, nclatoms, dxyzcl )
 
    character(len=*), intent(in)  :: extern_routine  ! Name of calling module
    _REAL_          , intent(in)  :: escf            ! SCF energy
    integer         , intent(in)  :: nqmatoms        ! Number of QM atoms
    _REAL_          , intent(in)  :: dxyzqm(:,:)     ! SCF QM force
    integer         , intent(in)  :: debug
    integer,optional, intent(in)  :: nclatoms        ! Number of MM atoms
    _REAL_, optional, intent(in)  :: dxyzcl(:,:)     ! SCF MM force

    integer :: i

    call debug_enter_function( 'print_results', module_name, debug )

    if ( debug > 1 ) then
      write (6,'(a,/)') 'print_results - final energy in kcal and gradient(s) in kcal/(mol*Å):'
      write (6,'(a)') extern_routine//' - final energy:'
      write(6,'(f20.8)') escf
  
      write(6,'(a)') extern_routine//' - final gradient(s):'
      write(6,'(a)') 'QM region:'
      do i = 1, nqmatoms
         write(6,'(3(x,f16.10))') dxyzqm(:, i)
      end do
      ! If the user specified MM forces, print these as well
      if ( present( nclatoms ) .and. present( dxyzcl ) ) then
        if ( nclatoms > 0 ) then
          write(6,'(a)') 'MM region:'
          do i = 1, nclatoms
            write(6,'(3(x,f16.10))') dxyzcl(:, i)
          end do
        end if
      end if
    end if

    call debug_exit_function( 'print_results', module_name, debug )
 
  end subroutine print_results
 
  ! -----------------------------------------------------------------------
  ! Check whether the program at 'path' is properly installed.
  ! If exist is specified, return whether the program was found.
  ! If required is true, call sander bomb if the program is missing.
  ! -----------------------------------------------------------------------
  subroutine check_installation( program, id, required, debug, found, path )
 
    implicit none
 
    character(len=*) , intent(in)    :: program
    character(len=*) , intent(in)    :: id
    logical          , intent(in)    :: required
    integer          , intent(in)    :: debug
    logical, optional, intent(out)   :: found
    character(len=*), optional, intent(out) :: path

    character(len=80) :: read_buffer
    character(len=80) :: call_buffer
    character(len=80) :: filename
    integer :: iunit = 77
    integer :: stat
    integer :: system
 
    call debug_enter_function( 'check_installation', module_name, debug )

    filename = 'extern_location'//trim(id)

    ! Search for executable
    call_buffer = 'which '//trim(program)//' > '//trim(filename)
    stat = system(trim(call_buffer))

    if ( stat == 0 ) then
      ! Found program
      write(6,'(a)') '| Program '//trim(program)//' found!'
      if ( present( found ) ) then
        found = .true.
      end if
    else if ( required ) then
      ! Required program not found
      call sander_bomb('check_installation (qm2_extern_util_module)', &
        'Executable "'//trim(program)//'" not found', &
        'Please check that "'//trim(program)//'" is in your path')
    else
      ! Program not found, but don't bomb since we will look for another
      ! Instead, set output and leave the function
      write(6,'(a)') '| Program '//trim(program)//' not found!'
      if ( present( found ) ) then
        found = .false.
      end if
      return
    end if

    ! If we found the program, print the full path 
    open (unit=iunit, file=trim(filename), form='formatted', iostat=stat)
    if ( stat /= 0 ) then
      call sander_bomb('check_installation (qm2_extern_util_module)', &
        'Internal error opening file with path location of executable.', &
        'Quitting now.')
    end if
    read (iunit, '(a)', iostat=stat) read_buffer
    if ( stat /= 0 ) then
      call sander_bomb('check_installation (qm2_extern_util_module)', &
        'Internal error reading from file with path location of executable.', &
        'Quitting now.')
    end if
    close (unit=iunit, status='delete', iostat=stat)
    if ( stat /= 0 ) then
      call sander_bomb('check_installation (qm2_extern_util_module)', &
        'Internal error closing and deleting file with path location of executable.', &
        'Quitting now.')
    end if

    write(6,'(2a,/)') '| Executable location: ', trim(read_buffer)

    if (present(path)) then
       path=read_buffer(1:index(read_buffer, '/', .true.)-1)
    end if

    if ( debug > 0) then
      call debug_exit_function( 'check_installation', module_name, debug )
    end if
 
  end subroutine check_installation
 
  ! -----------------------------------------------------------------------
  ! Write dipole data to the dipole moment property file 'dipfile'
  ! -----------------------------------------------------------------------
  subroutine write_dipole( dipfile, dipxyz, dipole, debug )

    use constants, only: CODATA08_AU_TO_DEBYE

    implicit none

    character(len=*), intent(in) :: dipfile
    _REAL_          , intent(in) :: dipxyz(:), dipole ! in Debye
    integer         , intent(in) :: debug
    integer :: system

    integer                      :: iunit = 351, ios, stat
    logical, save                :: first_call = .true.

    call debug_enter_function( 'write_dipole', module_name, debug )

    ! Remove any existing dipole file on first run
    if ( first_call ) then ! This is set to false later
      stat = system('rm -f '//dipfile)
    end if

    ! write dipole moment to dipole moment property file
    open (iunit, file=dipfile, position='append', iostat=ios)
    if ( ios /= 0 ) then
       call sander_bomb('write_dipole (qm2_extern_util_module)', &
            'Error opening file '//dipfile//' for appending.', &
            'Will quit now')
    end if

    ! Write out information on the first call
    if ( first_call ) then
      first_call=.false.
      write(iunit,'(a)') "| Dipole moment in Debye: {x, y, z}, |D|"
    end if
    ! Write dipole moment
    write(iunit,'(4f15.6)') dipxyz(:), dipole

    close(iunit)

    call debug_exit_function( 'write_dipole', module_name, debug )

  end subroutine write_dipole


  ! -------------------------
  ! Write charges to chgfile
  ! -------------------------
  subroutine write_charges(chgfile, charges, debug)

    implicit none

    character(len=*), intent(in) :: chgfile
    _REAL_, intent(in)  :: charges(:) ! charges from GAMESS
    integer, intent(in) :: debug

    integer        :: iunit = 351, ios, stat
    integer        :: system
    logical, save  :: first_call = .true.

    call debug_enter_function( 'write_charges', module_name, debug )

    ! Remove any existing charge file on first run
    if ( first_call ) then ! This is set to false later
      stat = system('rm -f '//chgfile)
      first_call = .false.
    end if

    ! write charges to charges property file
    open(iunit, file=chgfile, position='append', iostat=ios)
    if ( ios /= 0 ) then
      call sander_bomb('write_charges(qm2_extern_gms_module)', &
        'Error opening file '//chgfile//' for appending.', &
        'Will quit now')
    end if

    write(iunit,'(10f8.4)') charges(:)

    close(iunit)

    call debug_exit_function( 'write_charges', module_name, debug )

  end subroutine write_charges


  ! Now write point charges to 'chgfile'
  subroutine write_chgfile( chgfile, nclatoms, clcoords, debug )

    implicit none

    character(len=*), intent(in) :: chgfile
    integer         , intent(in) :: nclatoms
    _REAL_          , intent(in) :: clcoords(:,:)
    integer         , intent(in) :: debug

    integer                      :: iunit = 351, ios, i

    call debug_enter_function( 'write_chgfile', module_name, debug )

    if ( nclatoms > 0 ) then
      open(iunit, file=chgfile, iostat=ios)
      if ( ios > 0 ) then
        call sander_bomb('write_inpfile (qm2_extern_tc_module)', &
          'Error opening point charge coordinate file '//chgfile//' for writing', &
          'Will quit now')
      end if
      write(iunit,'(i6,/)') nclatoms
      do i = 1, nclatoms
         write(iunit,'(4f21.16)') clcoords(4,i), clcoords(:3,i)
      enddo
      close(iunit)
    end if

    call debug_exit_function( 'write_chgfile', module_name, debug )
  
  end subroutine write_chgfile

  subroutine debug_enter_function( funcname, modname, debug )

    implicit none    

    character(len=*), intent(in) :: funcname, modname
    integer         , intent(in) :: debug

    if ( debug > 0 ) then
      write (6, '(a)') '>>>>> Entered '//funcname//' ('//modname//')'
      call flush(6)
    end if

  end subroutine debug_enter_function
 
  subroutine debug_exit_function( funcname, modname, debug )

    implicit none    

    character(len=*), intent(in) :: funcname, modname
    integer         , intent(in) :: debug

    if ( debug > 0 ) then
      write (6, '(a)') '<<<<< Left '//funcname//' ('//modname//')'
      call flush(6)
    end if

  end subroutine debug_exit_function


end module qm2_extern_util_module
