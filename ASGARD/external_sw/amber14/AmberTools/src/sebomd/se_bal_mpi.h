
!     my_subs is an array of susbsytems (by number) that each PE in a
!     parallel job will own. These can be mixed and arranged across
!     PE-s in order to load balance.
!     my_numsubs is the number of subsystems this PE owns
!     time_subs is an array of times each susbsytem took to complete
!     in the mosub routine - this should be filled in by subsytem index
!     i.e. the time for subsystem 5 should go in time_subs(5) regardless
!     of which PE gets that subsystem

#ifdef MPI 
      double precision gtmp_mpi, virtmp_mpi, TMP_MPI, time_subs, watch
      integer my_numsubs,my_subs
      COMMON/se_bal_mpi/
     +               gtmp_mpi(3,maxatm),virtmp_mpi(4),
     +               TMP_MPI(MXDIAT),time_subs(MAXSUB),
     +               watch(24),my_subs(MAXSUB),my_numsubs
#endif
