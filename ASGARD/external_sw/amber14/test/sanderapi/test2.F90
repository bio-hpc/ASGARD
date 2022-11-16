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

    use sander_api, only: sander_input, gas_sander_input, pme_sander_input, &
                          sander_setup2, potential_energy_rec, energy_forces, &
                          sander_cleanup, qmmm_input_options, qm_sander_input, &
                          read_inpcrd_file, get_inpcrd_natom, set_positions, &
                          prmtop_struct, read_prmtop_file, destroy_prmtop_struct

    implicit none

    double precision, allocatable, dimension(:) :: forces 
    double precision, allocatable, dimension(:) :: coordinates
    double precision, dimension(6)              :: box
    type(potential_energy_rec) :: energies
    type(sander_input)         :: options
    type(qmmm_input_options)   :: qmmm_options
    type(prmtop_struct)        :: parm
    integer                    :: alloc_failed
    integer                    :: ierr
    integer                    :: natom

    logical :: failed = .false.
    logical :: failed2

    write(6, '(a)') 'Testing GB sander interface'
    call gas_sander_input(options, 7)
    options%cut = 9999.d0
    options%saltcon = 0.2d0
    options%gbsa = 1
    call read_prmtop_file("../gb7_trx/prmtop_an", parm, ierr)
    allocate(forces(3*parm%natom), coordinates(3*parm%natom))
    call read_inpcrd_file("../gb7_trx/trxox.2.4ns.x", coordinates, box, ierr)
    call sander_setup2(parm, coordinates, box, options, ierr=ierr)
    call energy_forces(energies, forces)

    ! Compare the energies to the output from the relevant test:
! BOND    =      631.8993  ANGLE   =      898.2543  DIHED      =      566.4453
! VDWAALS =     -768.3629  EEL     =    -7874.4913  EGB        =    -1943.0838
! 1-4 VDW =      348.8246  1-4 EEL =     5980.5047  RESTRAINT  =        0.0000
! ESURF   =       33.8338
    call compare(energies%bond, 631.8993d0, failed, 'Bond')
    call compare(energies%angle, 898.2543d0, failed, 'Angle')
    call compare(energies%dihedral, 566.4453d0, failed, 'Dihedral')
    call compare(energies%vdw_14, 348.8246d0, failed, '1-4 vdW')
    call compare(energies%elec_14, 5980.5047d0, failed, '1-4 Elec')
    call compare(energies%vdw, -768.3629d0, failed, 'van der Waals')
    call compare(energies%elec, -7874.4913d0, failed, 'Electrostatic')
    call compare(energies%gb, -1943.0838d0, failed, 'EGB')
    call compare(energies%surf, 33.8338d0, failed, 'SASA (GBSA)')

    call sander_cleanup

    if (failed) then
        write(6, '(a)') 'Possible FAILURE'
    else
        write(6, '(a)') 'PASSED'
    end if
    write(6, '(62("="))')

    ! Now try explicit solvent with PME

    write(6, '(a)') 'Testing PME sander interface'
    failed2 = failed
    failed = .false.
    call pme_sander_input(options)
    options%cut = 8.d0

    deallocate(forces, coordinates)

    call destroy_prmtop_struct(parm)
    call read_prmtop_file("../4096wat/prmtop", parm, ierr)
    call get_inpcrd_natom("../4096wat/eq1.x", natom)
    allocate(coordinates(3*natom))
    allocate(forces(3*parm%natom))
    call read_inpcrd_file("../4096wat/eq1.x", coordinates, box, ierr)
    call sander_setup2(parm, coordinates, box, options, ierr=ierr)

    call energy_forces(energies, forces)

    ! Compare the energies to the output from the relevant test:
! NSTEP =        1   TIME(PS) =       1.001  TEMP(K) =   298.28  PRESS =     0.0
! Etot   =    -32059.8471  EKtot   =      7282.8008  EPtot      =    -39342.6479
! BOND   =         0.0000  ANGLE   =         0.0000  DIHED      =         0.0000
! 1-4 NB =         0.0000  1-4 EEL =         0.0000  VDWAALS    =      6028.9517
! EELEC  =    -45371.5995  EHBOND  =         0.0000  RESTRAINT  =         0.0000
    call compare(energies%bond, 0.d0, failed, 'Bond')
    call compare(energies%angle, 0.d0, failed, 'Angle')
    call compare(energies%dihedral, 0.d0, failed, 'Dihedral')
    call compare(energies%vdw_14, 0.d0, failed, '1-4 vdW')
    call compare(energies%elec_14, 0.d0, failed, '1-4 Elec')
    call compare(energies%vdw, 6028.9517d0, failed, 'van der Waals')
    call compare(energies%elec, -45371.5995d0, failed, 'Electrostatic')

    deallocate(forces, coordinates)
    call sander_cleanup

    if (failed) then
        write(6, '(a)') 'Possible FAILURE'
    else
        write(6, '(a)') 'PASSED'
    end if
    write(6, '(62("="))')

    ! Now test the various QM/MM capabilities

    failed2 = failed2 .or. failed
    failed = .false.
    write(6, '(a)') 'Testing the QM/MM non-periodic interface'

    call gas_sander_input(options, 1)
    options%cut = 99.d0
    options%ifqnt = 1

    call qm_sander_input(qmmm_options)
    qmmm_options%iqmatoms(1:3) = (/ 8, 9, 10 /)
    qmmm_options%qm_theory = 'PM3'
    qmmm_options%qmcharge = 0
    qmmm_options%qmgb = 2
    qmmm_options%adjust_q = 0

    call destroy_prmtop_struct(parm)
    call read_prmtop_file("../qmmm2/lysine_PM3_qmgb2/prmtop", parm, ierr)
    call get_inpcrd_natom("../qmmm2/lysine_PM3_qmgb2/lysine.crd", natom)
    allocate(coordinates(3*natom))
    allocate(forces(3*parm%natom))
    call read_inpcrd_file("../qmmm2/lysine_PM3_qmgb2/lysine.crd", coordinates, &
                          box, ierr)
    call sander_setup2(parm, coordinates, box, options, qmmm_options, ierr=ierr)

    call energy_forces(energies, forces)

    call compare(energies%bond, 0.0016d0, failed, 'Bond')
    call compare(energies%angle, 0.3736d0, failed, 'Angle')
    call compare(energies%dihedral, 0.0026d0, failed, 'Dihedral')
    call compare(energies%vdw_14, 3.7051d0, failed, '1-4 vdW')
    call compare(energies%elec_14, 65.9137d0, failed, '1-4 Elec')
    call compare(energies%vdw, 0.1908d0, failed, 'van der Waals')
    call compare(energies%elec, -4.1241d0, failed, 'Electrostatic')
    call compare(energies%gb, -80.1406d0, failed, 'EGB')
    call compare(energies%scf, -11.9100d0, failed, 'QM Escf')

    deallocate(forces, coordinates)
    call sander_cleanup

    if (failed) then
        write(6, '(a)') 'Possible FAILURE'
    else
        write(6, '(a)') 'PASSED'
    end if
    write(6, '(62("="))')

    failed2 = failed2 .or. failed
    failed = .false.
    write(6, '(a)') 'Testing the QM/MM periodic interface (PM3-PDDG)'

    call pme_sander_input(options)
    options%cut = 8.d0
    options%ifqnt = 1
    options%jfastw = 4

    call qm_sander_input(qmmm_options)
    qmmm_options%qmmask = ':1-2'
    qmmm_options%qm_theory = 'PDDG-PM3'
    qmmm_options%qmcharge = 0
    qmmm_options%scfconv = 1.d-10
    qmmm_options%tight_p_conv = 1
    qmmm_options%qmmm_int = 5 ! mechanical embedding

    call destroy_prmtop_struct(parm)
    call read_prmtop_file("../qmmm2/MechEm_nma-spcfwbox/prmtop", parm, ierr)
    call get_inpcrd_natom("../qmmm2/MechEm_nma-spcfwbox/inpcrd", natom)
    allocate(coordinates(3*natom))
    allocate(forces(3*parm%natom))
    call read_inpcrd_file("../qmmm2/MechEm_nma-spcfwbox/inpcrd", coordinates, &
                          box, ierr)
    call sander_setup2(parm, coordinates, box, options, qmmm_options, ierr=ierr)

    call energy_forces(energies, forces)

    call compare(energies%bond, 605.7349d0, failed, 'Bond')
    call compare(energies%angle, 331.7679d0, failed, 'Angle')
    call compare(energies%dihedral, 0.0000d0, failed, 'Dihedral')
    call compare(energies%vdw_14, 0.0000d0, failed, '1-4 vdW')
    call compare(energies%elec_14, 0.0000d0, failed, '1-4 Elec')
    call compare(energies%vdw, 1281.8450d0, failed, 'van der Waals')
    call compare(energies%elec, -7409.7167d0, failed, 'Electrostatic')
    call compare(energies%scf, -37.1277d0, failed, 'QM Escf')

    ! Don't deallocate forces, since we want to change the Hamiltonian slightly
    ! but keep the same system
    call sander_cleanup
    call destroy_prmtop_struct(parm)

    if (failed) then
        write(6, '(a)') 'Possible FAILURE'
    else
        write(6, '(a)') 'PASSED'
    end if
    write(6, '(62("="))')

end program test_program
