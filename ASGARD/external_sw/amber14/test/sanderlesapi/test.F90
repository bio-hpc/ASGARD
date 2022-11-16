! This program tests the Fortran 90 sander interface to compute forces and
! energies

subroutine compare(computed, regression, failed, desc)
    
    implicit none

    double precision, intent(in) :: computed
    double precision, intent(in) :: regression
    character(len=*), intent(in) :: desc
    logical, intent(in out) :: failed

    ! Compare to 4 decimal places
    if (abs(computed - regression) .gt. 2.0d-4) then
        write(6, '(a,2(a,f15.4))') &
                desc, ' failed: Expected ', regression, ' got ', computed
        failed = .true.
    end if

end subroutine compare

program test_program

    use sanderles_api, only: sander_input, gas_sander_input, pme_sander_input,&
                             sander_setup, potential_energy_rec, energy_forces,&
                             sander_cleanup, sander_natom, read_inpcrd_file, &
                             get_inpcrd_natom

    implicit none

    double precision, allocatable, dimension(:) :: forces 
    double precision, allocatable, dimension(:) :: coordinates
    double precision, dimension(6)              :: box
    type(potential_energy_rec) :: energies
    type(sander_input)         :: options
    integer                    :: alloc_failed
    integer                    :: natom
    integer                    :: ierr

    logical :: failed = .false.
    logical :: failed2

    write(6, '(a)') 'Testing GB sander interface (diffcoords w/ RDT)'
    call gas_sander_input(options, 7)
    options%cut = 9999.d0
    options%rgbmax = 100.0d0;
    options%rdt = 0.01d0;

    call get_inpcrd_natom("../LES_GB/les.diffcoords.r", natom)
    allocate(coordinates(3*natom))
    call read_inpcrd_file("../LES_GB/les.diffcoords.r", coordinates, box, ierr)
    call sander_setup("../LES_GB/les.prm", coordinates, box, options, ierr=ierr)
    call sander_natom(natom)
    allocate(forces(3*natom), stat=alloc_failed)
    if (alloc_failed .ne. 0) then
        write(0,*) 'Failed allocating force array'
        stop 1
    end if

    call energy_forces(energies, forces)

    ! Compare the energies to the output from the relevant test:
! Etot   =       -13.2705  EKtot   =         0.0000  EPtot      =       -13.2705
! BOND   =        16.5749  ANGLE   =        21.5250  DIHED      =        35.5749
! 1-4 NB =         6.4411  1-4 EEL =       140.5502  VDWAALS    =        -4.6590
! EELEC  =      -198.7892  EGB     =       -30.4884  RESTRAINT  =         0.0000
    call compare(energies%bond, 16.5749d0, failed, 'Bond')
    call compare(energies%angle, 21.5250d0, failed, 'Angle')
    call compare(energies%dihedral, 35.5749d0, failed, 'Dihedral')
    call compare(energies%vdw_14, 6.4411d0, failed, '1-4 vdW')
    call compare(energies%elec_14, 140.5502d0, failed, '1-4 Elec')
    call compare(energies%vdw, -4.6590d0, failed, 'van der Waals')
    call compare(energies%elec, -198.7892d0, failed, 'Electrostatic')
    call compare(energies%gb, -30.4884d0, failed, 'EGB')
    call compare(energies%surf, 0.d0, failed, 'SASA (GBSA)')

    call sander_cleanup
    deallocate(coordinates, forces)

    if (failed) then
        write(6, '(a)') 'Possible FAILURE'
    else
        write(6, '(a)') 'PASSED'
    end if
    write(6, '(62("="))')

    ! Now try explicit solvent with PME

    write(6, '(a)') 'Testing GB sander interface (samecoords w/out RDT)'
    failed2 = failed
    failed = .false.
    options%rdt = 0.d0

    call get_inpcrd_natom("../LES_GB/les.samecoords.r", natom)
    allocate(coordinates(3*natom))
    call read_inpcrd_file("../LES_GB/les.samecoords.r", coordinates, box, ierr)
    call sander_setup("../LES_GB/les.alt.prm", coordinates, box, options, ierr=ierr)
    call sander_natom(natom)
    allocate(forces(3*natom))

    call energy_forces(energies, forces)

    ! Compare the energies to the output from the relevant test:
! NSTEP =        0   TIME(PS) = 1000010.000  TEMP(K) =     0.00  PRESS =     0.0
! Etot   =       -27.4464  EKtot   =         0.0000  EPtot      =       -27.4464
! BOND   =         5.8375  ANGLE   =        19.0846  DIHED      =        32.7197
! 1-4 NB =         7.1039  1-4 EEL =       141.3377  VDWAALS    =        -3.0346
! EELEC  =      -202.2822  EGB     =       -28.2130  RESTRAINT  =         0.0000
    call compare(energies%bond, 5.8375d0, failed, 'Bond')
    call compare(energies%angle, 19.0846d0, failed, 'Angle')
    call compare(energies%dihedral, 32.7197d0, failed, 'Dihedral')
    call compare(energies%vdw_14, 7.1039d0, failed, '1-4 vdW')
    call compare(energies%elec_14, 141.3377d0, failed, '1-4 Elec')
    call compare(energies%vdw, -3.0346d0, failed, 'van der Waals')
    call compare(energies%elec, -202.2822d0, failed, 'Electrostatic')
    call compare(energies%gb, -28.2130d0, failed, 'EGB')

    call sander_cleanup
    deallocate(forces, coordinates)

    if (failed) then
        write(6, '(a)') 'Possible FAILURE'
    else
        write(6, '(a)') 'PASSED'
    end if
    write(6, '(62("="))')

    write(6, '(a)') 'Testing PME sander interface'
    failed2 = failed
    failed = .false.
    call pme_sander_input(options)
    options%cut = 8.d0

    call get_inpcrd_natom("../LES/md.LES.x", natom)
    allocate(coordinates(3*natom))
    call read_inpcrd_file("../LES/md.LES.x", coordinates, box, ierr)
    call sander_setup("../LES/LES.prmtop.save", coordinates, box, options, ierr=ierr)
    call sander_natom(natom)
    allocate(forces(3*natom), stat=alloc_failed)
    if (alloc_failed .ne. 0) then
        write(0, *) 'Failed allocating force array'
        stop 1
    end if

    call energy_forces(energies, forces)

    ! Compare the energies to the output from the relevant test:
! NSTEP =        1   TIME(PS) =      10.002  TEMP(K) =   289.54  PRESS =     0.0
! Etot   =     -3114.6572  EKtot   =       963.1779  EPtot      =     -4077.8351
! BOND   =        14.7095  ANGLE   =        34.6208  DIHED      =        35.3483
! 1-4 NB =        13.0097  1-4 EEL =       274.1453  VDWAALS    =       545.6397
! EELEC  =     -4995.3084  EHBOND  =         0.0000  RESTRAINT  =         0.0000
    
    call compare(energies%bond, 14.7095d0, failed, 'Bond')
    call compare(energies%angle, 34.6208d0, failed, 'Angle')
    call compare(energies%dihedral, 35.3483d0, failed, 'Dihedral')
    call compare(energies%vdw_14, 13.0097d0, failed, '1-4 vdW')
    call compare(energies%elec_14, 274.1453d0, failed, '1-4 Elec')
    call compare(energies%vdw, 545.6397d0, failed, 'van der Waals')
    call compare(energies%elec, -4995.3084d0, failed, 'Electrostatic')

    deallocate(forces, coordinates)
    call sander_cleanup

    if (failed) then
        write(6, '(a)') 'Possible FAILURE'
    else
        write(6, '(a)') 'PASSED'
    end if
    write(6, '(62("="))')

    ! Exit without success if either step failed
    if (failed .or. failed2) stop 1

end program test_program
