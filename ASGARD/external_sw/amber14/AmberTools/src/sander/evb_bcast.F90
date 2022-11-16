! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  MPI_BCAST EVB data to all PE                                           |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#  include "../include/dprec.fh"

   subroutine evb_bcast

   use evb_parm, only: nevb, nxch, xch_type 
   use schlegel, only: ndg, ncoord, ncart, natm, xdat_min, xdat_ts, xdat_xdg &
                     , nbond, nangle, ndihed, ibond, iangle, idihed &
                     , atomic_numbers, use_cartesian

   implicit none

#  include "parallel.h"
   include 'mpif.h'

   !  ..........................................................................

   integer :: n, ierr

!  +---------------------------------------------------------------------------+
!  |  Broadcast array size variables for DG-EVB                                |
!  +---------------------------------------------------------------------------+

! write(6,*) 'evb_cast 1'

   if( trim( adjustl( xch_type ) ) == "dist_gauss" ) then
      call mpi_bcast (  nbond, 1, MPI_INTEGER, 0, commworld, ierr )
      call mpi_bcast ( nangle, 1, MPI_INTEGER, 0, commworld, ierr )
      call mpi_bcast ( ndihed, 1, MPI_INTEGER, 0, commworld, ierr )
      call mpi_bcast ( ncoord, 1, MPI_INTEGER, 0, commworld, ierr )
      call mpi_bcast (  ncart, 1, MPI_INTEGER, 0, commworld, ierr )
      call mpi_bcast (   natm, 1, MPI_INTEGER, 0, commworld, ierr )
      call mpi_bcast ( use_cartesian, 1, MPI_LOGICAL, 0, commworld, ierr )
   endif

!dEVB   call mpi_barrier ( commworld, ierr )

!  +---------------------------------------------------------------------------+
!  |  Allocate memory for DG-EVB                                               |
!  +---------------------------------------------------------------------------+

! write(6,*) 'evb_cast 2'

   if( worldrank /= 0 ) call evb_init_alloc

!  +---------------------------------------------------------------------------+
!  |  Broadcast coordinate, energy, gradient and hessian for DG-EVB            |
!  +---------------------------------------------------------------------------+

   select case( trim( adjustl( xch_type ) ) )

      case( "dist_gauss" )

         call mpi_bcast (  ibond,  nbond*2, MPI_INTEGER, 0, commworld, ierr )
         call mpi_bcast ( iangle, nangle*3, MPI_INTEGER, 0, commworld, ierr )
         call mpi_bcast ( idihed, ndihed*4, MPI_INTEGER, 0, commworld, ierr )
         call mpi_bcast ( atomic_numbers, natm, MPI_INTEGER, 0, commworld, ierr )

         do n = 1, nevb
            call mpi_bcast ( xdat_min(n)%q, ncoord, MPI_DOUBLE_PRECISION &
                           , 0, commworld, ierr )
            call mpi_bcast ( xdat_min(n)%d, ncoord, MPI_DOUBLE_PRECISION &
                           , 0, commworld, ierr )
            call mpi_bcast ( xdat_min(n)%k, ncoord*ncoord, MPI_DOUBLE_PRECISION &
                           , 0, commworld, ierr )
            call mpi_bcast ( xdat_min(n)%v, 1, MPI_DOUBLE_PRECISION &
                           , 0, commworld, ierr )
            call mpi_bcast ( xdat_min(n)%qcart, natm*3, MPI_DOUBLE_PRECISION &
                           , 0, commworld, ierr )
            call mpi_bcast ( xdat_min(n)%filename, 512, MPI_CHARACTER &
                           , 0, commworld, ierr )
         enddo

         do n = 1, nxch
            call mpi_bcast ( xdat_ts(n)%q, ncoord, MPI_DOUBLE_PRECISION &
                           , 0, commworld, ierr )
            call mpi_bcast ( xdat_ts(n)%d, ncoord, MPI_DOUBLE_PRECISION &
                           , 0, commworld, ierr )
            call mpi_bcast ( xdat_ts(n)%k, ncoord*ncoord, MPI_DOUBLE_PRECISION &
                           , 0, commworld, ierr )
            call mpi_bcast ( xdat_ts(n)%v, 1, MPI_DOUBLE_PRECISION &
                           , 0, commworld, ierr )
            call mpi_bcast ( xdat_ts(n)%qcart, natm*3, MPI_DOUBLE_PRECISION &
                           , 0, commworld, ierr )
            call mpi_bcast ( xdat_ts(n)%filename, 512, MPI_CHARACTER &
                           , 0, commworld, ierr )
         enddo

         do n = 1, ndg
            call mpi_bcast ( xdat_xdg(n)%q, ncoord, MPI_DOUBLE_PRECISION &
                           , 0, commworld, ierr )
            call mpi_bcast ( xdat_xdg(n)%d, ncoord, MPI_DOUBLE_PRECISION &
                           , 0, commworld, ierr )
            call mpi_bcast ( xdat_xdg(n)%k, ncoord*ncoord, MPI_DOUBLE_PRECISION &
                           , 0, commworld, ierr )
            call mpi_bcast ( xdat_xdg(n)%v, 1, MPI_DOUBLE_PRECISION &
                           , 0, commworld, ierr )
            call mpi_bcast ( xdat_xdg(n)%qcart, natm*3, MPI_DOUBLE_PRECISION &
                           , 0, commworld, ierr )
            call mpi_bcast ( xdat_xdg(n)%filename, 512, MPI_CHARACTER &
                           , 0, commworld, ierr )
         enddo

   end select 

! write(6,*) 'evb_cast 3'

!dEVB   call mpi_barrier ( commworld, ierr )

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine evb_bcast


