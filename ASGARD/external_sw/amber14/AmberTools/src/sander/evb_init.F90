! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  EVB Initialization                                                     |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#include "../include/assert.fh"
#include "../include/dprec.fh"

   subroutine evb_init

   use evb_parm, only: nevb, nxch, xch_type, ioe 
   use schlegel, only: dist_gauss, ndg, ncoord, natm, xdat_min, xdat_ts &
                     , xdat_xdg, init_ready, nbond, nangle, ndihed &
                     , use_cartesian

   implicit none

#  include "parallel.h"

   !  ..........................................................................

   integer :: n, ifound, ios
   logical :: lread
   character(  1) :: ctype
   character( 40) :: label
   character(512) :: cline 

!  +---------------------------------------------------------------------------+
!  |  Obtain system dimension for Schlegel's DG-EVB                            |
!  +---------------------------------------------------------------------------+

   if( worldrank == 0 ) then

      select case( trim( adjustl( xch_type ) ) )
         case( "dist_gauss" )
!dEVB       open( ioe, file = trim( adjustl( xdat_min(1)%filename ) ) )
            open( ioe, file = trim( adjustl( xdat_xdg(1)%filename ) ) )
            lread = .true.
         case default
            lread = .false.
      end select

      if( lread ) then

         ifound   = 0

         select case( trim( adjustl( dist_gauss%xfile_type ) ) )

!  +---------------------------------------------------------------------------+
!  |  Gaussian formatted checkpoint file                                       |
!  +---------------------------------------------------------------------------+
            case( "gaussian_fchk" )
               do
                  read( ioe, '(A)', iostat = ios ) cline
!  .............................................................................
!  :  Check for EOF and read errors                                            :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
                  select case( ios )
                     case(:-1)          ! end of file encountered
                       write(6,'(A)') "EOF but did not find requested data in " &
                                   //  trim( adjustl( xdat_xdg(1)%filename ) )
                       call mexit(6,1)
                     case(1:)           ! error during read
                        write(6,'(A)') "ERROR encountered while reading " &
                                    //  trim( adjustl( xdat_xdg(1)%filename ) )
                        call mexit(6,1)
                  end select
!  .............................................................................
!  :  Read No. of atoms and internal coordinates                               :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
                  select case( trim( adjustl( cline(1:40) ) ) )
                     case( 'Number of atoms' )
                        backspace( ioe )
                        read( ioe, 100 ) label, ctype, natm
                        ifound = ifound + 1
                     case( 'Redundant internal dimensions' )
                        backspace( ioe )
                        read( ioe, 100 ) label, ctype, n
                        read( ioe, 2000 ) ncoord, nbond, nangle, ndihed
                        ifound = ifound + 1
                  end select
                  if( ifound  == 2 ) then
                     write(6,200) 'No. of redundant internal bonds     = ', nbond
                     write(6,200) 'No. of redundant internal angles    = ', nangle
                     write(6,200) 'No. of redundant internal dihedrals = ', ndihed
                     write(6,200) 'No. of internal coordinates         = ', ncoord
                     exit
                  endif
               enddo

!  .............................................................................
!  :  .EVB formatted file                                                      :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
            case( "old_EVB" )
               do
                  read( ioe, '(A)', iostat = ios ) cline
!  .............................................................................
!  :  Check for EOF and read errors                                            :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
                  select case( ios )
                     case(:-1)          ! end of file encountered
                       write(6,'(A)') "EOF but did not find requested data in " &
                                   //  trim( adjustl( xdat_xdg(1)%filename ) )
                       call mexit(6,1)
                     case(1:)           ! error during read
                        write(6,'(A)') 'Error encountered while reading ' &
                                    //  trim( adjustl( xdat_xdg(1)%filename ) )
                        call mexit(6,1)
                  end select
!  .............................................................................
!  :  Read coordinate type, No. of atoms & internal coordinates                :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
                  select case( trim( adjustl( cline ) ) )
                     case( '[coordinate type]' )
                        read(ioe,'(A)') cline
                        if( trim( adjustl( cline ) ) == "use_cartesian" ) &
                           use_cartesian = .true. 
                     case( '[external evb data dimension]' )
                        read( ioe, '(5I12)' ) ncoord, natm, nbond, nangle, ndihed
                        ifound = ifound + 1 
                  end select
                  if( ifound == 1 ) exit
               enddo

         end select

         close( ioe )

      endif

!  +---------------------------------------------------------------------------+
!  |  Allocate memory & read coordinate, energy, gradient and hessian          |
!  +---------------------------------------------------------------------------+

      if( .not. init_ready ) then
         call evb_init_alloc
         init_ready = .true.
      endif

      if( trim( adjustl( xch_type ) ) == "dist_gauss" ) then

         select case( trim( adjustl( dist_gauss%xfile_type ) ) )

!  .............................................................................
!  :  Gaussian formatted checkpoint file                                       :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
            case( "gaussian_fchk" )
               do n = 1, nevb
                  call evb_gfchk ( xdat_min(n), ioe )
               enddo
               do n = 1, nxch
                  call evb_gfchk ( xdat_ts(n), ioe )
               enddo
               do n = 1, ndg
                  call evb_gfchk ( xdat_xdg(n), ioe )
               enddo

!  .............................................................................
!  :  .EVB formatted file                                                      :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
            case( "EVB", "old_EVB" )
               do n = 1, nevb
                  call evb_io ( xdat_min(n), ioe )
               enddo
               do n = 1, nxch
                  call evb_io ( xdat_ts(n), ioe )
               enddo
               do n = 1, ndg
                  call evb_io ( xdat_xdg(n), ioe )
               enddo

         end select

      endif

   endif 

 100  format( A40, 3X, A1, 5X, I12 )
 200  format( A, I12 )
 2000 format( 4I12 )
 
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine evb_init


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Allocate storage for DG-EVB                                            |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   subroutine evb_init_alloc

   use evb_parm, only: nevb, nxch, xch_type
   use schlegel, only: ndg, xdat_min, xdat_ts, xdat_xdg, natm, atomic_numbers &
                     , nbond, nangle, ndihed, ibond, iangle, idihed

   implicit none

   !  ..........................................................................

   integer :: n, alloc_error


   select case( trim( adjustl( xch_type ) ) )

      case( "dist_gauss" )

         allocate( atomic_numbers(natm), stat = alloc_error )
         REQUIRE( alloc_error == 0 )

         if( nbond > 0 ) then
            allocate( ibond(nbond,2), stat = alloc_error )
            REQUIRE( alloc_error == 0 )
         endif
         if( nangle > 0 ) then
            allocate( iangle(nangle,3), stat = alloc_error )
            REQUIRE( alloc_error == 0 )
         endif
         if( ndihed > 0 ) then
            allocate( idihed(ndihed,4), stat = alloc_error )
            REQUIRE( alloc_error == 0 )
         endif

         do n = 1, nevb
            call evb_gfchk_alloc ( xdat_min(n) )
         enddo
         do n = 1, nxch
            call evb_gfchk_alloc ( xdat_ts(n) )
         enddo
         do n = 1, ndg
            call evb_gfchk_alloc ( xdat_xdg(n) )
         enddo

   end select

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine evb_init_alloc


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Deallocate storage for DG-EVB                                          |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   subroutine evb_init_dealloc

   use evb_parm, only: nevb, nxch, xch_type
   use schlegel, only: ndg, xdat_min, xdat_ts, xdat_xdg, atomic_numbers &
                     , nbond, nangle, ndihed, ibond, iangle, idihed

   implicit none

   !  ..........................................................................

   integer :: n, dealloc_error


   select case( trim( adjustl( xch_type ) ) )

      case( "dist_gauss" )

         deallocate( atomic_numbers, stat = dealloc_error )

         if( nbond > 0 ) then
            deallocate( ibond, stat = dealloc_error )
            REQUIRE( dealloc_error == 0 )
         endif
         if( nangle > 0 ) then
            deallocate( iangle, stat = dealloc_error )
            REQUIRE( dealloc_error == 0 )
         endif
         if( ndihed > 0 ) then
            deallocate( idihed, stat = dealloc_error )
            REQUIRE( dealloc_error == 0 )
         endif

         do n = 1, nevb
            call evb_gfchk_dealloc ( xdat_min(n) )
         enddo
         do n = 1, nxch
            call evb_gfchk_dealloc ( xdat_ts(n) )
         enddo
         do n = 1, ndg
            call evb_gfchk_dealloc ( xdat_xdg(n) )
         enddo

   end select

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine evb_init_dealloc

