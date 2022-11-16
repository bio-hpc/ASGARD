subroutine se_deeprec(energy,eecrt)
!
! compute machine precision for current energy
!
  use sebomd_module, only : sebomd_obj
  implicit none

! subroutine parameters
  double precision :: energy, eecrt

! local variables
  double precision :: a, b, c
  integer :: iprec

  iprec = sebomd_obj%iprec
  a = energy
  b = 0.1d0
  c = energy
  do while ( (a+b).ne.(c) )
    b = b*0.1d0
  end do
  eecrt = b*(10.0d0**iprec)
  return
end subroutine se_deeprec
