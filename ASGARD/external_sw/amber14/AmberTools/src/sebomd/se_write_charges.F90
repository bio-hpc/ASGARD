subroutine se_write_charges(natoms, edc, iatnum, xyz, atchg, atchg2, atchg3, dipole_atm)
  use file_io_dat, only : sechgunit
  implicit none
  integer :: natoms
  integer, dimension(natoms) :: iatnum
  double precision :: edc
  double precision, dimension (3, natoms) :: xyz
  double precision, dimension (3, natoms) :: dipole_atm
  double precision, dimension (natoms) :: atchg
  double precision, dimension (natoms) :: atchg2
  double precision, dimension (natoms) :: atchg3

  integer :: i

  write(sechgunit,'(i6)') natoms
  write(sechgunit,'(a,f18.8)') "ESEBOMD = ",edc
  do i = 1, natoms
    write(sechgunit,'(i3,1x,10f12.7)')  &
        iatnum(i) &
       ,xyz(1,i) &
       ,xyz(2,i) &
       ,xyz(3,i) &
       ,atchg(i) &
       ,dipole_atm(1,i) &
       ,dipole_atm(2,i) &
       ,dipole_atm(3,i) &
       ,atchg2(i)  &
       ,atchg3(i)
  end do
end subroutine se_write_charges
