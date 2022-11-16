      subroutine se_etimer(t)

!     returns the elapsed cpu time in seconds

      implicit none
#ifdef MPI
#include "mpif.h"
#endif

      double precision t

      double precision tnow, t0

      logical first
      data first /.true./
      save t0,first

#ifdef MPI 
      tnow = mpi_wtime()
#else
#  ifdef __PGI
      call cpu_time(tnow)
#  else
#     if (__GNUC__ == 4 && __GNUC_MINOR__ <= 1)
      tnow = second()
#     else
      tnow = 0.01*mclock()
#     endif
#  endif
#endif
      if(first)then
         t0 = tnow
         first = .false.
      endif
      t = tnow - t0
      return
      end
