subroutine se_mpi_vars(nproc, myid, commsebomd)
!
! transfer MPI values from sander to sebomd
!
  implicit none

! subroutine parameters
  integer :: nproc
  integer :: myid
  integer :: commsebomd

! external variables
#include "../sander/parallel.h"

  nproc = sandersize
  myid = sanderrank
  commsebomd = commsander
  return
end subroutine se_mpi_vars
