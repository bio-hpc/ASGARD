subroutine se_init(nat, boxsander, xyzsander, gradsander, chgsander,&
                   maxatm, natoms, xyz, chgpme, dbox, dhalf, boxvol)
  implicit none
  integer :: nat      ! from sander
  double precision :: boxsander(3)
  double precision :: xyzsander(3*nat)
  double precision :: gradsander(3*nat)
  double precision :: chgsander(nat)
  integer :: natoms   ! to sebomd
  integer :: maxatm
  double precision :: xyz(3,maxatm)
  double precision :: chgpme(maxatm)
  double precision :: dbox(3), dhalf(3), boxvol

  integer :: i

  natoms = nat

  if (natoms.gt.maxatm) then
    write(6,'(" SEBOMD ERROR: Too many atoms: increase maxatm value in", &
              "sebomd.dim")')
    call mexit(6,1)
  endif

  ! gradient to zero
  do i = 1, 3*nat
    gradsander(i) = 0.0d0
  end do

  ! transfer xyz coordinates and atomic charges
  do i = 1, natoms
    xyz(1,i) = xyzsander(i*3-2)
    xyz(2,i) = xyzsander(i*3-1)
    xyz(3,i) = xyzsander(i*3-0)
    chgpme(i) = chgsander(i)
  end do

  ! transfer box coordinates
  dbox(1) = boxsander(1)
  dbox(2) = boxsander(2)
  dbox(3) = boxsander(3)
  dhalf(1) = 0.5d0*dbox(1)
  dhalf(2) = 0.5d0*dbox(2)
  dhalf(3) = 0.5d0*dbox(3)
  boxvol = dbox(1)*dbox(2)*dbox(3)

  return
end subroutine se_init
