      subroutine se_mpi_allreduce_dble(values,ncount)

! routine to perform mpi_allreduce for double precision numbers
! size of the buffer is gtmp_mpi(3,maxatm) as defined in bal_mpi.h
! input: values(ncount)

      implicit none
#include "sebomd.dim"
#include "se_parproc.h"
#include "se_bal_mpi.h"
#ifdef MPI
#include "mpif.h"
#endif
      double precision values(*)
      integer ncount
#ifdef MPI
      integer ier
      integer i
      do i = 1, ncount
        tmp_mpi(i) = 0.0d0
      end do
      call mpi_allreduce(values,tmp_mpi,ncount,
     .                   MPI_DOUBLE_PRECISION,MPI_SUM,
     .                   commsebomd,ier)
      do i = 1, ncount
        values(i) = tmp_mpi(i)
      end do
#else
      return
#endif
      end
