subroutine se_backinbox(natoms, nres, xyz, dbox, irpnt)
  ! put all atoms back in the box (PBC)
  ! works for orthorombic cells
  ! box is centered on (0,0,0)
  implicit none
  integer :: natoms
  integer :: nres
  integer :: irpnt(nres+1)
  double precision :: dbox(3)
  double precision :: dhalf(3)
  double precision :: xyz(3, natoms)

  integer :: i, ires, i1, i2
  double precision :: boxinv(3)
  double precision :: x1, y1, z1
  double precision :: xi, yi, zi
  double precision :: dxij, dyij, dzij

  dhalf(1) = 0.5d0*dbox(1)
  dhalf(2) = 0.5d0*dbox(2)
  dhalf(3) = 0.5d0*dbox(3)

  ! step 1: all back in the box (atom based)
  boxinv(1) = 1/dbox(1)
  boxinv(2) = 1/dbox(2)
  boxinv(3) = 1/dbox(3)

  do i = 1, natoms
    xyz(1, i) = xyz(1,i)-anint(xyz(1,i)*boxinv(1))*dbox(1)
    xyz(2, i) = xyz(2,i)-anint(xyz(2,i)*boxinv(2))*dbox(2)
    xyz(3, i) = xyz(3,i)-anint(xyz(3,i)*boxinv(3))*dbox(3)
  end do

  x1 = xyz(1,1)
  y1 = xyz(2,1)
  z1 = xyz(3,1)
  do i = 1, natoms
    xi = xyz(1,i)
    yi = xyz(2,i)
    zi = xyz(3,i)
    dxij = x1-xi
    dyij = y1-yi
    dzij = z1-zi
    if(abs(dxij).gt.dhalf(1)) xi = xi + sign(dbox(1),dxij)
    if(abs(dyij).gt.dhalf(2)) yi = yi + sign(dbox(2),dyij)
    if(abs(dzij).gt.dhalf(3)) zi = zi + sign(dbox(3),dzij)
    xyz(1,i) = xi
    xyz(2,i) = yi
    xyz(3,i) = zi
  end do

  ! step 2: correct for residues based
  ! (reference = first atom of the residue)
  do ires = 1, nres
    i1 = irpnt(ires)
    i2 = irpnt(ires+1)-1
    x1 = xyz(1,i1)
    y1 = xyz(2,i1)
    z1 = xyz(3,i1)
    do i = i1+1, i2
!     call se_pbcxyz(i1,i,xi,yi,zi)
      ! apply PBC for i with i1 as a reference (same code as pbcxyz)
      xi = xyz(1,i)
      yi = xyz(2,i)
      zi = xyz(3,i)
      dxij = x1-xi
      dyij = y1-yi
      dzij = z1-zi
      if(abs(dxij).gt.dhalf(1)) xi = xi + sign(dbox(1),dxij)
      if(abs(dyij).gt.dhalf(2)) yi = yi + sign(dbox(2),dyij)
      if(abs(dzij).gt.dhalf(3)) zi = zi + sign(dbox(3),dzij)
      xyz(1,i) = xi
      xyz(2,i) = yi
      xyz(3,i) = zi
    end do
  end do
  return
end subroutine se_backinbox
