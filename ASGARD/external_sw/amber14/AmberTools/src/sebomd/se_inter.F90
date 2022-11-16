! module to define between two atoms whether the interaction is intermolecular
! or intramolecular
module se_inter
  implicit none
  private
  ! subroutines
  public :: se_inter_init
  public :: se_interres
  ! arrays
  integer, allocatable, dimension(:) :: molnum
  ! variables

contains
  subroutine se_inter_init(natoms, nmol, pmolid, symbol2)
    use sebomd_module, only: sebomd_obj
    implicit none
    integer :: natoms ! number of atoms
    integer :: nmol   ! number of molecules
    integer, dimension(nmol) :: pmolid ! number of atoms per molecule
    character*4, dimension(natoms) :: symbol2

    integer :: i, mol, ind, nmolfound, nsolventfound
    logical :: solute

    ! do this only for PIF2 and PIF3 (ctype = 1 or 2)
    if ((sebomd_obj%ctype /= 1) .and. (sebomd_obj%ctype /= 2)) return

    ! nmol is the number of molecules as found in the topology file
    ! (atoms_per_molecule flag). if the flag is not found, then nmol=1

    ! in the case of nmol = 1, we will check for water molecules using 
    ! the symbol2 array (atom type array)

    ! molnum array is defined as follow:
    ! for each atom index i, molnum(i) = + the molecule index if solute
    !                                  = - the molecule index if solvent

    if (.not.allocated(molnum)) then
      allocate(molnum(natoms))

      if (nmol.eq.1) then
        mol = 1
        nmolfound=1
        nsolventfound=0
        solute = .false.
        do i = 1, natoms
          if (symbol2(i).eq.'OW  ') then
            mol = mol+1
            nsolventfound = nsolventfound + 1
            nmolfound = nmolfound + 1
            molnum(i) = -mol ! solvent case (= OW or HW)
          else if (symbol2(i).eq.'HW  ') then
            molnum(i) = -mol ! solvent case (= OW or HW)
          else
            molnum(i) = mol  ! solute case
            solute = .true.
          endif
        end do
        if (.not.solute) nmolfound = nmolfound - 1 ! no solute thus nmolfound
                                                   ! should have been initialized to zero

      else ! more than one molecule: atoms_per_molecule is set
        ind = 0
        nmolfound=0
        nsolventfound=0
        do mol = 1, nmol
          nmolfound = nmolfound + 1
          do i = 1, pmolid(mol)
            if (symbol2(ind+i).eq.'OW  ') then
              molnum(ind+i) = -mol ! solvent case (= OW or HW)
              nsolventfound = nsolventfound + 1
            else if (symbol2(ind+i).eq.'HW  ') then
              molnum(ind+i) = -mol ! solvent case (= OW or HW)
            else
              molnum(ind+i) = mol  ! solute case
            endif
          end do
          ind = ind+pmolid(mol)
        end do
      endif
    endif
    write(6,'("")')
    write(6,'("SEBOMD: PIF info: ",i5," molecules found, including ",i5," solvent molecules")') nmolfound, nsolventfound
!   write(6,'("se_inter_init",20i4)') (molnum(i),i=1,natoms)
    return
  end subroutine se_inter_init
!--------------------------------------------------------------------------------
  subroutine se_interres(i,j,resinter)
    implicit none
! determines whether the interaction between atom i and atom j is an
! inter-residue interaction or an intra-residue interaction

    integer :: resinter
    integer :: i
    integer :: j
    integer :: imol
    integer :: jmol

    resinter = 0      ! = 0 -> intra-residue interaction
                         ! = 1 -> inter-residue interaction
                         !        water-water
                         ! = 2 -> inter-residue interaction
                         !        solute-water
    imol = molnum(i)
    jmol = molnum(j)

    ! if imol /= jmol then
    !              +-----------+
    !              |    imol   |
    !              +-----+-----+
    !              | < 0 | > 0 |
    !   +----+-----+-----+-----+
    !   |    |     |     |     |
    !   | j  | < 0 |  1  |  2  |
    !   | m  |     |     |     |
    !   | o  +-----+-----+-----+
    !   | l  |     |     |     |
    !   |    | > 0 |  2  |error|
    !   |    |     |     |     |
    !   +----+-----+-----+-----+

    if (imol == jmol) then
      resinter = 0 ! same molecule = intramolecular
    else if (imol*jmol < 0) then
      resinter = 2 ! one negative / one positive = inter solute/solvent
    else
      if ((imol < 0).and.(jmol < 0)) then
        resinter = 1 ! both negative = inter solvent/solvent
      else
        write(6,'("SEBOMD: PIF Error")')
        write(6,'("SEBOMD: solute-solute intermolecular interaction between atom ",i4," and atom ",i4)') i,j
        write(6,'("SEBOMD: this is not defined in the PIF scheme")')
        write(6,'("SEBOMD: exiting...")')
        call mexit(6,1)
      endif
    endif
  end subroutine se_interres
end module se_inter
