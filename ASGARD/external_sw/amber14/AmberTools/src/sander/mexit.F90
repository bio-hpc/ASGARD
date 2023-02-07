#include "copyright.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Platform independent exit; returns an exit status to the OS.
subroutine mexit(output_unit, status)
   
   !  mexit() - machine-dependent exit() procedure, designed to return an
   !            appropriate (success/failure) value to the operating system.
   
#ifdef PUPIL_SUPPORT
   use pupildata
#endif

   implicit none
   integer output_unit  ! close this unit if greater than zero, non-MPI
   integer status       ! exit status; error if non-zero

#ifdef MPI
   include 'mpif.h'
   integer ierr
#  include "parallel.h"
   
   !       ...status .gt. 0 implies an error condition, therefore
   !       kill all the nodes.  mpi_abort on the world communicator
   !       should do this, but it does not on some implementations.
   !       some MPI's have an environmental sledge hammer that kills
   !       every MPI process if one dies: mpiexec -kill
   
   if (status /= 0) then
      call amflsh(output_unit)
      call mpi_abort(MPI_COMM_WORLD, status, ierr)
   else
      call mpi_finalize(ierr)
   end if
#endif

#ifdef PUPIL_SUPPORT
!jtc ========================= PUPIL INTERFACE =========================
!     Terminate the PUPIL CORBA interface, only if such an interface
!     exists.
      if (pupactive) then
         puperror = 0
         call killcorbaintfc(puperror)
         if(puperror /= 0) write(6,*) 'Error ending PUPIL CORBA interface.'
      end if
!jtc ========================= PUPIL INTERFACE =========================
#endif

   if (output_unit > 0 .and. status/=0) then
      close(unit=output_unit)
   end if

#ifdef XLF90
   if (status /= 0) then
      stop 1
   else
      stop 0
   end if
#else
   call exit(status)
#endif
end subroutine mexit 
