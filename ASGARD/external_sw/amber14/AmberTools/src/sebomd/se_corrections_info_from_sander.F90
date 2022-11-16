subroutine se_corrections_info_from_sander(mm_pep, ctype, k_pep)
!
! transfer information from sander to sebomd
! here: data from SE corrections (core-core functions)
!

! modules
  use sebomd_module, only: sebomd_obj
  implicit none

! subroutine parameters
  logical :: mm_pep
  integer :: ctype
  double precision :: k_pep

  mm_pep = .false.
  ctype = sebomd_obj%ctype
  if (sebomd_obj%peptcorr.ne.0) then
    mm_pep =  .true.
  endif
  k_pep = sebomd_obj%peptk
end subroutine se_corrections_info_from_sander
