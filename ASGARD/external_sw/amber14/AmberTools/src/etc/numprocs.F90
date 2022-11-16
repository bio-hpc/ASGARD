! This program is just a quick program designed to print out how many
! processors are assigned to the world communicator. This is a *better* way
! of counting the # of lines printed by $DO_PARALLEL echo "I'm here" or trying
! to parse $DO_PARALLEL for -n or -np (since this isn't required on all
! platforms). Indeed, the best way of determining how many processors are
! assigned to an MPI program is to actually count them and return it.  That's
! what's done here.

program count_procs

   implicit none
   include 'mpif.h'

   integer  :: rank
   integer  :: numprocs
   logical  :: master
   integer  :: ier

   call MPI_Init(ier)

   call MPI_Comm_size(MPI_COMM_WORLD, numprocs, ier)

   call MPI_Comm_rank(MPI_COMM_WORLD, rank, ier)

   master = rank .eq. 0

   ! I can't imagine you'd be using more than 999 processors for the
   ! test...
   if (master) write(6, '(I3)') numprocs

   call MPI_Finalize(ier)

end program count_procs

