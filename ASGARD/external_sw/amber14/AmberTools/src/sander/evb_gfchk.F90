! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Read ab initio energy, gradient and hessian from Gaussian formatted    |#
! #|  checkpoint file                                                        |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#include "../include/assert.fh"
#include "../include/dprec.fh"

   subroutine evb_gfchk ( xdat, ioe )

   use schlegel, only: berny_type, ncoord, natm, atomic_numbers, ibond &
                     , iangle, idihed 

   implicit none

#  include "parallel.h"

   integer, intent(in) :: ioe
   type( berny_type ), intent(inout) :: xdat

   !  ..........................................................................

   integer :: ios, m, n, nn, mbond, mangle, mdihed, alloc_error, dealloc_error
   integer :: ndx(ncoord*4), ndx1, ndx2, ndx3, ndx4
#ifdef DEBUG_EVB
   _REAL_  :: fdeps
#endif
   _REAL_, allocatable :: grad_cart(:), hess_cart(:,:), th(:), th_cart(:) 
   _REAL_, allocatable :: bmat(:,:), binv(:,:), dbmat(:,:), mass(:)

   logical :: read_qint = .false.

   character(  1) :: ctype
   character( 40) :: label
   character(512) :: cline

!  .............................................................................
!  :  Allocate temporary arrays                                                :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   allocate( grad_cart(natm*3), hess_cart(natm*3,natm*3), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( th(ncoord*(ncoord+1)/2), th_cart(natm*3*(natm*3+1)/2), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( bmat(ncoord,natm*3), binv(natm*3,ncoord), stat = alloc_error ) 
   REQUIRE( alloc_error == 0 )
   allocate( dbmat(ncoord,natm*3*natm*3), mass(natm), stat = alloc_error ) 
   REQUIRE( alloc_error == 0 )

!  +---------------------------------------------------------------------------+
!  |  Read EVB data from external files                                        |
!  +---------------------------------------------------------------------------+

   if( worldrank == 0 ) then 

      write(6,'(A)') "DG::  Opening file: " // trim( adjustl( xdat%filename ) )

      open( ioe, file = trim( adjustl( xdat%filename ) ) )

      do 

         read( ioe, '(A)', iostat = ios ) cline
!  .............................................................................
!  :  Check for EOF and read errors                                            :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
         select case( ios ) 
            case(:-1)          ! end of file encountered
               exit 
            case(1:)           ! error during read 
               write(6,'(A)') "Error encountered while reading " &
                           //  trim( adjustl( xdat%filename ) )
               call mexit(6,1)
         end select 

         select case( trim( adjustl( cline(1:40) ) ) )

!  +---------------------------------------------------------------------------+
!  |  Read energy                                                              |
!  +---------------------------------------------------------------------------+
            case( "Total Energy" )
               backspace( ioe )
               read( ioe, 200 ) label, ctype, xdat%v

!  +---------------------------------------------------------------------------+
!  |  Read atomic numbers                                                      |
!  +---------------------------------------------------------------------------+
            case( "Atomic numbers" )
               backspace( ioe )
               read( ioe, 100 ) label, ctype, nn
               read( ioe, 2000 ) ( atomic_numbers(n), n = 1, nn )

!  +---------------------------------------------------------------------------+
!  |  Read redundant internal coordinate definitions                           |
!  :...........................................................................:
!  |  1 -- 2                3                  1                               |
!  |                      /   \                  \                             |
!  |                    1       2                  2 -- 3                      |
!  |                                                      \                    |
!  |                                                        4                  |
!  +---------------------------------------------------------------------------+
            case( "Redundant internal coordinate indices" )
               backspace( ioe )
               read( ioe, 100 ) label, ctype, nn
               if( nn /= ncoord * 4 ) then
                  write(6,'(A)') "ERROR: the number of internal coordinate " &
                              // "indices is inconsistent with ncoord * 4"
                  write(6,'(A,I12)') "ncoord = ", ncoord
                  write(6,'(A,I12)') "# read = ", nn/4 
                  call mexit(6,1)
               endif
               read( ioe, 2000 ) ( ndx(n), n = 1, nn )

!  .............................................................................
!  :  Populate internal coordinate definitions                                 :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
               mbond  = 0
               mangle = 0
               mdihed = 0
               do n = 1, ncoord
                  ndx1 = ndx( (n-1)*4 + 1 )
                  ndx2 = ndx( (n-1)*4 + 2 )
                  ndx3 = ndx( (n-1)*4 + 3 )
                  ndx4 = ndx( (n-1)*4 + 4 )
                  if( ndx3 == 0 ) then
                     mbond = mbond + 1
                     ibond(mbond,1) = ndx1
                     ibond(mbond,2) = ndx2
                  else if( ndx4 == 0 ) then
                     mangle = mangle + 1
                     iangle(mangle,1) = ndx1
                     iangle(mangle,2) = ndx3
                     iangle(mangle,3) = ndx2
                  else
                     mdihed = mdihed + 1
                     idihed(mdihed,1) = ndx1
                     idihed(mdihed,2) = ndx2
                     idihed(mdihed,3) = ndx3
                     idihed(mdihed,4) = ndx4
                  endif
               enddo

!  +---------------------------------------------------------------------------+
!  |  Read redundant internal coordinates                                      |
!  +---------------------------------------------------------------------------+
!dEVB       case( 'Redundant internal coordinates', 'Redundant Internal Coordinates' )
            case( "Redundant internal coordinates" )
               backspace( ioe )
               read( ioe, 100 ) label, ctype, nn
               read( ioe, 1000 ) ( xdat%q(n), n = 1, nn )
               read_qint = .true.

!  +---------------------------------------------------------------------------+
!  |  Read Cartesian coordinate                                                |
!  +---------------------------------------------------------------------------+
            case( "Current cartesian coordinates" )
               backspace( ioe )               
               read( ioe, 100 ) label, ctype, nn
               read( ioe, 1000 ) ( xdat%qcart(n), n = 1, nn )

!  +---------------------------------------------------------------------------+
!  |  Read atomic masses                                                       |
!  +---------------------------------------------------------------------------+
            case( "Real atomic weights" )
               backspace( ioe )
               read( ioe, 100 ) label, ctype, nn
               read( ioe, 1000 ) ( mass(n), n = 1, nn )

!  +---------------------------------------------------------------------------+
!  |  Read redundant internal forces                                           |
!  +---------------------------------------------------------------------------+
            case( "Internal Forces" )
               backspace( ioe )
               read( ioe, 100 ) label, ctype, nn
               read( ioe, 1000 ) ( xdat%d(n), n = 1, nn )
               xdat%d(:) = - xdat%d(:)

!  +---------------------------------------------------------------------------+
!  |  Read Cartesian coordinate gradient                                       |
!  +---------------------------------------------------------------------------+
            case( "Cartesian Gradient" )
               backspace( ioe )
               read( ioe, 100 ) label, ctype, nn
               read( ioe, 1000 ) ( grad_cart(n), n = 1, nn )

!  +---------------------------------------------------------------------------+
!  |  Read (lower triangle) internal coordinate hessian                        |
!  +---------------------------------------------------------------------------+
            case( "Internal Force Constants" )
               backspace( ioe )
               read( ioe, 100 ) label, ctype, nn
               read( ioe, 1000 ) ( th(n), n = 1, nn )

               if( nn /= ncoord * ( ncoord + 1 ) / 2 ) then
                  write(6,'(A)') "ERROR: the number of lower-T internal hessian " &
                              // "elements read is inconsistent with ncoord"
                  write(6,'(A,I12)') "ncoord = ", ncoord
                  write(6,'(A,I12,A)') "expect ", ncoord*(ncoord+1)/2, "elements"
                  write(6,'(A,I12,A)') "... but read ", nn, "elements"
                  call mexit(6,1)
               endif 

               nn = 0
               do m = 1, ncoord
                  do n = 1, m 
                     nn = nn + 1
                     xdat%k(m,n) = th(nn)
                     xdat%k(n,m) = xdat%k(m,n)
                  enddo
               enddo

!  +---------------------------------------------------------------------------+
!  |  Read (lower triangle) Cartesian coordinate hessian                       |
!  +---------------------------------------------------------------------------+
            case( "Cartesian Force Constants" )
               backspace( ioe )
               read( ioe, 100 ) label, ctype, nn
               read( ioe, 1000 ) ( th_cart(n), n = 1, nn )

               if( nn /= natm*3 * ( natm*3 + 1 ) / 2 ) then
                  write(6,'(A)') "ERROR: the number of lower-T Cartesian hessian " &
                              // "elements read is inconsistent with natm*3"
                  write(6,'(A,I12)') "natm*3 = ", natm * 3
                  write(6,'(A,I12,A)') "expect ", natm*3*(natm*3+1)/2, "elements"
                  write(6,'(A,I12,A)') "... but read ", nn, "elements"

                  call mexit(6,1)
               endif

               nn = 0
               do m = 1, natm*3
                  do n = 1, m
                     nn = nn + 1
                     hess_cart(m,n) = th_cart(nn)
                     hess_cart(n,m) = hess_cart(m,n)
                  enddo
               enddo

         end select

      enddo 

      close( ioe )

   endif 

!  +---------------------------------------------------------------------------+
!  |  Transform Cartesian data to redundant internals                          |
!  +---------------------------------------------------------------------------+

   if( .not. read_qint ) call cart2internal ( xdat%qcart, xdat%q )


#ifdef DEBUG_EVB

   call wdc_bmat ( xdat%qcart, bmat )

   call bmat_inv (  bmat, binv, natm, ncoord )

   fdeps = 0.001
   call bmat_grad ( xdat%qcart, dbmat, natm, fdeps )

   call evb_io_debug ( xdat%qcart, grad_cart, hess_cart, mass, bmat, binv &
                     , dbmat, xdat, ioe )
#endif

!  .............................................................................
!  :  Deallocate temporary arrays                                              :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   deallocate( grad_cart, hess_cart, th, th_cart, stat = dealloc_error )
   REQUIRE( dealloc_error == 0 )
   deallocate( bmat, binv, dbmat, mass, stat = dealloc_error )
   REQUIRE( dealloc_error == 0 )

  100 format( A40, 3X, A1, 5X, I12 )
  200 format( A40, 3X, A1, 5X, E22.15 )
 1000 format( 5( 1PE16.8 ) )
 2000 format( 6I12 )

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine evb_gfchk


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Allocate storage for DG ab initio data                                 |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   subroutine evb_gfchk_alloc ( xdat )

   use schlegel, only: berny_type, ncoord, natm

   implicit none

   type( berny_type ) :: xdat

   !  ..........................................................................

   integer :: alloc_error


   allocate( xdat%q(ncoord), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

   allocate( xdat%d(ncoord), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

   allocate( xdat%k(ncoord,ncoord), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

   allocate( xdat%qcart(natm*3), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine evb_gfchk_alloc


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Deallocate storage for xdat                                            |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   subroutine evb_gfchk_dealloc ( xdat )

   use schlegel, only: berny_type

   implicit none

   type( berny_type ) :: xdat

   !  ..........................................................................

   integer :: dealloc_error


   deallocate( xdat%q, xdat%d, xdat%k, xdat%qcart, stat = dealloc_error )
   REQUIRE( dealloc_error == 0 )

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine evb_gfchk_dealloc


