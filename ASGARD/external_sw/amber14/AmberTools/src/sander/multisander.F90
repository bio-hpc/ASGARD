#include "copyright.h"
#include "../include/dprec.fh"
#include "../include/assert.fh"
#include "ncsu-config.h"
!
!        SANDER, version 14
!
! The Molecular Dynamics/NMR Refinement/Modeling Module of the AMBER Package.
!     sander is a subroutine.  See below.
!
!     MULTISANDER: multiple sander module
!     Originally by Guanglei Cui (Simmerling lab), Jan. 2002
!
!     Modified by tec3 to allow runtime specification of the number of
!     copies to run (-ng) and to allow completely different input files.
!     Also added the capability to allocate the processors in groups
!     sequentially or dispersed (to do dispersed use the -ng-nonsequential
!     argument).  Added comments.  April 2002.

!     Changed the names of the communicators to be more descriptive.
!     Allowed specification of input/output files for each group to be specified
!     in a file (tec3, Jan 03).  Spliced the code into the main distribution
!     such that sander is minimally modified and acts like normal as long as
!     group related runtime/command line arguments are not specified.
!     TODO: Modify so that each group can run using a different number of
!           processors.
!
! Developers notes:
!
! o  All MPI communications should be done using the new communicators rather
!     than MPI_COMM_WORLD.  A number of new communicators are defined:
!     CommSander  -- communications within a given sander job
!                    (replaces MPI_COMM_WORLD)
!     CommWorld   -- communications to ALL processors across multiple sander
!                    jobs
!     CommMaster  -- communications to the master node of each separate sander
!                    job each has corresponding size and rank,
!                    i.e. MasterRank, MasterSize
!
! o  mdfil() has been modified to allow reading command line input from a
!       file or string and some new arguments have been added.
!    -ng N
!       ...set the number of groups (multiple instances of sander) to N.
!    -groupfile FILE
!       ...reads command line input from a file for each of the N groups using
!          mdfil.  Note the command line for the groupfile FILE for each group
!          must be on a single line and must be less than or equal to 512
!          characters.
!          If this option is not present, assume that the I/O file names
!          have 000, 001, ... appended to the names for each group.
!    -ng-nonsequential
!       ...allocate the processors per group via arithmetic modulo numgroup
!          rather than sequentially.
!    -gpes FILE (future option)
!       ...specify the number processors for this particular group/sander run.
!
! o  If you add files to files.h or file_io_dat, you may need to update the code 
!    here to broadcast the names appropriately when -groupfile is not specified.
!
! TODO:
! o  generalize the error handling to allow one job dying not to kill the
!     rest (on error) and the opposite (if one job dies for any reason, kill
!     all the others) based on specification of a command line argument.
! o  modify the code to allow differing numbers of processors per sander
!     instance.


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Main program
!-----------------------------------------------------------------------
!     --- MULTISANDER ---
!-----------------------------------------------------------------------
!  Setup MPI and file handling. 
!  Call sander to perform calculations.

program multisander

#ifdef MPI
   use remd, only : rem, rremd, repnum, numreps
#endif
   use commandline_module, only : mdfil
   use file_io_dat
#ifdef BINTRAJ
   use AmberNetcdf_mod, only: NC_setupAmberNetcdf
#endif

#ifndef DISABLE_NCSU
   use ncsu_sander_hooks, only: ncsu_on_multisander_exit => on_multisander_exit
#endif /* DISABLE_NCSU */

#ifdef MPI
#  ifdef LES
   use evb_pimd, only: evb_begin, evb_end, lpimd_dcrypt, nslice_per_group &
                     , nslice_per_node, lpimd_group
#  endif /* LES */
#endif 

   implicit none

   ! Update this when the version changes! (Make sure to update the len if
   ! necessary)
   ! Also update call to NC_setupAmberNetcdf below
   character (len=4), parameter :: VERSION_STRING = "14.0"
   logical :: version_requested

#include "parallel.h"

#ifdef MPI

#  include "ew_parallel.h"
   include 'mpif.h'
#  include "../include/md.h"

#ifdef LES
   integer :: m, n, nn, ndx, lpimd_id, alloc_error
#endif /* LES */

   integer ierror,  masterid, i

   ! set default value for rem here, not in mdfil()
   ! if it is in mdfil, -rem needs to be in the groupfile too
   ! this way we can call mdfil again and it won't reset rem value
   ! rem value must be set on the command line, not on the group file (like -ng)
   rem=0

   ! Set default for RREMD. This tells us the type of reservoir.
   ! This is set as a command line flag to sander: 0 - No reservoir (std. REMD),
   ! 1 - Boltzmann weighted, 2 - 1/N weighted, 3 - user defined weights
   rremd=0

   call mpi_init(ierror)

   !   Set up the world communicators and associated variables:

   CommWorld = MPI_COMM_WORLD
   call mpi_comm_rank( CommWorld, worldrank,  ierror )
   call mpi_comm_size( CommWorld, worldsize,  ierror )
   call mpi_barrier( CommWorld, ierror )

   mytaskid = worldrank
   numtasks = worldsize

   ! Call mdfil to set the names of the input and output files and
   ! process other command line arguments.

   call initialize_fnames

   if (worldrank == 0) then
      ng_sequential = .true.
      numgroup = 1
      num_recip = -1
#ifdef LES
      nslice = 0
#endif
      ! (mdfil may modify these depending on the command line options)
      version_requested = .false.
      call mdfil(VERSION_STRING, version_requested)
   end if

   ! Everyone should know if we're stopping
   call mpi_bcast(version_requested, 1, mpi_logical, 0, commworld, ierror)

   if (version_requested) call mexit(6,0)

#ifdef LES
   call mpi_bcast(nslice, 1, mpi_integer, 0, commworld, ierror)
#endif

   ! rem should be known to each processor 
   call mpi_bcast(rem, 1, mpi_integer, 0, commworld, ierror)
   
   ! rremd too
   call mpi_bcast(rremd, 1, mpi_integer, 0, commworld, ierror)

   !  Broadcast and validate the number of groups:

   call mpi_bcast(numgroup, 1, mpi_integer, 0, commworld, ierror)

   if (numgroup > worldsize) then
      if (worldrank == 0) then
         write(6,*) 'Error: specified more groups (', numgroup, &
               ') than the number of processors (', numtasks, ') !'
      end if
      call mexit(6,1)
   end if

   if ( mod(worldsize, numgroup ) /= 0 ) then
      if (worldrank == 0) then
         write(6,*) 'Error: the number of processors ', &
               'is not a multiple of the number of groups!'
      end if
      call mexit(6,1)
   end if

   !  Broadcast the group file information:

   call mpi_bcast(groups, len(groups), mpi_character, 0, commworld, ierror)

   !  Fix up the NUM_RECIP semantics:

   call mpi_bcast(num_recip, 1, mpi_integer, 0, commworld, ierror)
   if (num_recip == -1) then
      num_recip = numtasks/numgroup
      num_direct= numtasks/numgroup
   else
      num_direct=(numtasks/numgroup)-num_recip
   end if

   !  Processor allocation:

   call mpi_bcast(ng_sequential, 1, mpi_logical, 0, commworld, ierror)
   if ( ng_sequential ) then
      nodeid = worldrank / (worldsize / numgroup)
   else
      nodeid = mod(worldrank, numgroup)
   end if
   call mpi_barrier(commworld, ierror)

   !  Create a communicator for each group of -ng NumGroup processors:

   commsander = mpi_comm_world
   sandersize = worldsize
   sanderrank = worldrank
   if (numgroup > 1) then
      commsander = mpi_comm_null
      call mpi_comm_split(commworld, nodeid, worldrank, &
            commsander, ierror)
      if (commsander == mpi_comm_null) then
         if (worldrank == 0) then
            write(6,'(a,i5,a,i5)') 'Error: NULL Communicator on PE ', &
                  worldrank, ' from group ', nodeid
         end if
         call mexit(6,1)
      end if
      call mpi_comm_size(commsander, sandersize, ierror)
      call mpi_comm_rank(commsander, sanderrank, ierror)
   end if

   !  Define a communicator (CommMaster) that only talks between the local
   !  "master" in each group.  This is equivalent to a SanderRank .eq. 0

   masterid = 0
   masterrank = MPI_UNDEFINED
   mastersize = 0
   if (numgroup > 1) then
      commmaster = mpi_comm_null
      if(sanderrank /= 0) then
         masterid = MPI_UNDEFINED
      end if

      call mpi_comm_split(commworld, masterid, worldrank, &
            commmaster, ierror)
      ! will this be emitted when using the default MPI error handler ?
      if (ierror /= MPI_SUCCESS) then
         write(6,*) 'Error: MPI_COMM_SPLIT error ', ierror, &
               ' on PE ', worldrank
      end if

      if(commmaster /= mpi_comm_null) then
         call mpi_comm_size(commmaster, mastersize, ierror)
         call mpi_comm_rank(commmaster, masterrank, ierror)
      end if

   end if

   !  +-----------------------------------------------------------------------+
   !  | Define a communicator for multisander-EVB/LES-PIMD & assign bead      |
   !  |  tasks to PEs via evb_begin & evb_end                                 |
   !  +-----------------------------------------------------------------------+

#ifdef LES

   nslice_per_group = nslice / numgroup
   lpimd_id = worldrank / ( worldsize / numgroup )
   lpimd_group = lpimd_id

   lpimd_size = 0
   lpimd_rank = MPI_UNDEFINED

   if( numgroup > 1 .and. nslice > 0 ) then
      comm_lpimd = mpi_comm_null
      if( sanderrank > nslice_per_group - 1 ) lpimd_id = MPI_UNDEFINED
      call mpi_comm_split(commworld, lpimd_id, worldrank, &
            comm_lpimd, ierror)
      if( ierror /= MPI_SUCCESS ) then
         write(6,'(2(A,2X,I8))') 'COMM_LPIMD Error: MPI_COMM_SPLIT error ' &
                , ierror, ' on PE ', worldrank
      end if

      if( comm_lpimd /= MPI_COMM_NULL ) then

         call mpi_comm_size(comm_lpimd, lpimd_size, ierror)
         call mpi_comm_rank(comm_lpimd, lpimd_rank, ierror)


         if( nslice_per_group > lpimd_size ) then
            nslice_per_node = nslice_per_group / lpimd_size
         else
            nslice_per_node = 1
         endif

         allocate( evb_begin(lpimd_size), evb_end(lpimd_size) &
                                        , stat = alloc_error )
         REQUIRE( alloc_error == 0 )
         do n = 1, lpimd_size
            ndx = lpimd_group * nslice_per_group
            evb_begin(n) = ( n - 1 ) * nslice_per_node + 1 + ndx
            evb_end  (n) = n * nslice_per_node + ndx
         enddo

         allocate( lpimd_dcrypt(lpimd_size,nslice_per_node) &
                                          , stat = alloc_error )
         REQUIRE( alloc_error == 0 )
         do n = 1, lpimd_size
            nn = 0
            do m = evb_begin(n), evb_end(n)
               nn = nn + 1
               lpimd_dcrypt(n,nn) = m
            enddo
         enddo

!        if( worldrank == 0 ) then
!           write(6,'( A, 34(I4, X) )') '| evb_begin = ', evb_begin(:)
!           write(6,'( A, 34(I4, X) )') '| evb_end   = ', evb_end(:)
!           do n = 1, nslice_per_node
!              write(6,'( A, I2, A, 34(I4, X) )') &
!                 '| lpimd_dcrypt(:,', n, ') = ', lpimd_dcrypt(:,n)
!           enddo
!        endif
      endif
   endif

#endif /* LES */

   !  Setup mytaskid and numtasks for each group:

   mytaskid = sanderrank
   numtasks = sandersize

   !  Determine and communicate the file information to the masters:

   if (numgroup > 1) then

      !   Use the -groupfile command line option to specify the input file
      !   names for each processor.  This is done assuming *each* master has
      !  access to the groupfile.

      if (sanderrank == 0) then
         call amopen(15,groups,'O','F','R')
         i = 0
         do while (i < numgroup )
            read(15,'(a512)') groupbuffer

            if ( groupbuffer(1:1) /= '#' .and. &
               groupbuffer(1:1) /= '/' .and. &
               len_trim(groupbuffer) > 0) then

               ! Each master should read its own line in the groupfile
               if (i == masterrank) call mdfil

               i = i + 1
            end if
         end do
         close( 15 )
      end if

      ! broadcast rem again in case user put it in the groupfile
      ! Also broadcast rremd. 

      call mpi_bcast(rem, 1, mpi_integer, 0, commworld, ierror)
      call mpi_bcast(rremd, 1, mpi_integer, 0, commworld, ierror)

      call mpi_barrier(commworld, ierror)

      !   Print summary of multisander run:

      if (worldrank == 0) then
         write(6, '(a)') ''
         write(6, '(a)') ' Running multisander version of sander Amber14'
         write(6, '(a,i5)') '    Total processors = ', worldsize
         write(6, '(a,i5)') '    Number of groups = ', numgroup
         if ( .not. ng_sequential) then
            write(6, '(a)') '    Allocation of processors is non-sequential.'
         end if
         write(6, '(a)') ' '
      end if

#if 0
      if (debugremd) then
         if (worldrank == 0) then
         write(6,*) '    Looping over processors:'
            write(6,*) '       WorldRank is the global PE rank'
            write(6,*) '       NodeID is the local PE rank in current group'
            write(6,*) ''
         end if
         do i = 1, worldsize
            if (worldrank == i-1) then
               if (sanderrank == 0) then
                  write(6,*) '       Group     = ', masterrank
               end if
               write(6,*) '       WorldRank = ', worldrank
               write(6,*) '       NodeID    = ', sanderrank
               write(6,*) ''
            end if
            call mpi_barrier(commworld, ierror)
         end do
      end if
#endif

   end if  ! (numgroup > 1)

   ! ------====== REMD INITIALIZATION ======------
   if (rem /= 0) then

      ! Set replica # 
      repnum = nodeid + 1
      numreps = numgroup

   endif ! rem /= 0
   ! ------====== END REMD INITIALIZATION ======------

#else /* non-MPI section below: */

   numgroup = 1
   mytaskid = 0
   numtasks = 1
   version_requested = .false.
   call mdfil(VERSION_STRING, version_requested)
   if (version_requested) call mexit(6,0)

#endif /* MPI */
#  ifdef BINTRAJ
   ! Activate NetCDF interface.
   call NC_setupAmberNetcdf(6, "sander", VERSION_STRING)
#  endif
   call sander()

   !  Clean up and exit

#ifdef MPI

   if(numgroup > 1 .and. commmaster /= mpi_comm_null) then
      call mpi_comm_free(commmaster, ierror)
      commmaster = mpi_comm_null
   end if
   if(numgroup > 1 .and. commsander /= mpi_comm_null) then
      call mpi_comm_free(commsander, ierror)
      commsander = mpi_comm_null
   end if

#endif

#ifndef DISABLE_NCSU
   call ncsu_on_multisander_exit()
#endif /* DISABLE_NCSU */

   call mexit(-8,0)

end program multisander
