! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Read ab initio energy, gradient and hessian from .EVB file.            |#
! #|  xdat% is allocated in evb_gfchk_alloc.                                 |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#include "../include/assert.fh"
#include "../include/dprec.fh"

   subroutine evb_io ( xdat, ioe )

   use schlegel, only: berny_type, ncoord, natm, nbond, nangle, ndihed, ncart &
                     , ibond, iangle, idihed, atomic_numbers, use_cartesian 

   implicit none

#  include "parallel.h"

   integer, intent(in) :: ioe
   type( berny_type ), intent(inout) :: xdat

   !  ..........................................................................

   integer :: ios, m, n, nn, mbond, mangle, mdihed, ndx1, ndx2, ndx3, ndx4
   integer :: alloc_error, dealloc_error
   _REAL_  :: fdeps
   _REAL_, allocatable :: mass(:), qcart(:) 
   _REAL_, allocatable :: ih(:), ch(:), grad_cart(:), hess_cart(:,:) 
   _REAL_, allocatable :: bmat(:,:), binv(:,:), dbmat(:,:), hess_tmp(:,:)
   logical :: read_qint = .false.
   character(512) :: cline

!  .............................................................................
!  :  Allocate temporary arrays                                                :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   allocate( ih(ncoord*(ncoord+1)/2), ch(natm*3*(natm*3+1)/2), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( mass(natm), qcart(natm*3), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( grad_cart(natm*3), hess_cart(natm*3,natm*3), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( bmat(ncoord,natm*3), binv(natm*3,ncoord), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( dbmat(ncoord,natm*3*natm*3), hess_tmp(natm*3,natm*3), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

!  +---------------------------------------------------------------------------+
!  |  Read EVB data from external files                                        |
!  +---------------------------------------------------------------------------+

   mbond  = 0
   mangle = 0
   mdihed = 0

   if( worldrank == 0 ) then 

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

         select case( trim( adjustl( cline ) ) )

!  +---------------------------------------------------------------------------+
!  |  Read ab initio energy                                                    |
!  +---------------------------------------------------------------------------+
            case( "[electronic energy]" )
               read( ioe, * ) xdat%v
               read( ioe, * )

!  +---------------------------------------------------------------------------+
!  |  Read atomic numbers                                                      |
!  +---------------------------------------------------------------------------+
            case( "atomic numbers" )
               read( ioe, 2000 ) ( atomic_numbers(n), n = 1, natm )
               read( ioe, * )

!  +---------------------------------------------------------------------------+
!  |  Read redundant internal coordinate definitions                           |
!  :''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''':
!  |  1 -- 2              3                  1                                 |
!  |                    /   \                  \                               |
!  |                  1       2                  2 -- 3                        |
!  |                                                    \                      |
!  |                                                      4                    |
!  +---------------------------------------------------------------------------+
            case( "[redundant internal indices]" )
               do n = 1, nbond + nangle + ndihed
                  read( ioe, 3000 ) ndx1, ndx2, ndx3, ndx4
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
!  |  Cartesian                                                                |
!  +---------------------------------------------------------------------------+
!  .............................................................................
!  :  Read cartesian coordinate                                                :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
            case( "[cartesian coordinates]" )
               read( ioe, 1000 ) ( xdat%qcart(n), n = 1, natm*3)
               read( ioe, * )

!  .............................................................................
!  :  Read cartesian gradient                                                  :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
            case( "[cartesian gradient]" )
               read( ioe, 1000 ) ( grad_cart(n), n = 1, natm*3)
               read( ioe, * )

!  .............................................................................
!  :  Read (lower triangle) cartesian hessian                                  :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
            case( "[cartesian hessian]" )
               read( ioe, 1000 ) ( ch(n), n = 1, natm*3*(natm*3+1)/2 )
               read( ioe, * )
               nn = 0
               do m = 1, natm * 3
                  do n = 1, m
                     nn = nn + 1
                     hess_cart(m,n) = ch(nn)
                     hess_cart(n,m) = hess_cart(m,n)
                  enddo
               enddo

!  +---------------------------------------------------------------------------+
!  |  Internals                                                                |
!  +---------------------------------------------------------------------------+
!  .............................................................................
!  :  Read converted redundant internal coordinates                            :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
            case( "[redundant internal coordinates]" )
               read( ioe, 1000 ) ( xdat%q(n), n = 1, ncoord )
               read( ioe, * )
               read_qint = .true.

!  .............................................................................
!  :  Read redundant internal gradient                                         :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
            case( "[redundant internal gradient]" )
               read( ioe, 1000 ) ( xdat%d(n), n = 1, ncoord )
               read( ioe, * )

!  .............................................................................
!  :  Read redundant internal hessian                                          :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
            case( "[redundant internal hessian]" )
               read( ioe, 1000 ) ( ih(n), n = 1, ncoord*(ncoord+1)/2 )
               read( ioe, * )
               nn = 0
               do m = 1, ncoord
                  do n = 1, m
                     nn = nn + 1
                     xdat%k(m,n) = ih(nn)
                     xdat%k(n,m) = xdat%k(m,n)
                  enddo
               enddo

         end select

      enddo 

      close( ioe )

!  +---------------------------------------------------------------------------+
!  |  Use Cartesian coordinates for DG fit                                     |
!  +---------------------------------------------------------------------------+

      if( use_cartesian ) then
         nbond  = 0
         nangle = 0
         ndihed = 0
         ncart  = natm * 3
         ncoord = natm * 3
         xdat%q(  :) = xdat%qcart(  :)
         xdat%d(  :) =  grad_cart(  :)
         xdat%k(:,:) =  hess_cart(:,:)
      endif

   endif 

!  +---------------------------------------------------------------------------+
!  |  Convert Cartesians to internals, so as to perform DG fit in internals    |
!  +---------------------------------------------------------------------------+

   if( .not. use_cartesian .and. .not. read_qint ) then

      call cart2internal ( xdat%qcart, xdat%q )

!  .............................................................................
!  :  Transform gradients from Cartesian to internal frame                     :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
!  :  [g] = [B^-1]^t [G]                                                       :
!  .............................................................................

      call wdc_bmat ( xdat%qcart, bmat )
      call bmat_inv (  bmat, binv, natm, ncoord )

      xdat%d(:) = matmul( transpose(binv), grad_cart )

!  .............................................................................
!  :  Transform Hessian from Cartesian to internal frame                       :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
!  :  [h] = [(B^-1)^t] ( [H] - [(B')^t] [g] ) [B^-1]                           :
!  .............................................................................

      fdeps = 0.0001
      call bmat_grad ( xdat%qcart, dbmat, natm, fdeps )

      hess_tmp(:,:) =  hess_cart(:,:) - reshape( matmul( &
                    transpose(dbmat), xdat%d ), (/ natm*3, natm*3 /) )

      xdat%k(:,:) = matmul( transpose(binv), matmul( hess_tmp, binv ) )


   endif

#ifdef DEBUG_EVB
   call evb_io_debug2 ( xdat%qcart, grad_cart, hess_cart, mass, xdat, ioe )
#endif 

!  .............................................................................
!  :  Deallocate temporary arrays                                              :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   deallocate( mass, qcart, ih, ch, grad_cart, hess_cart, stat = dealloc_error )
   REQUIRE( dealloc_error == 0 )
   deallocate( bmat, binv, dbmat, hess_tmp, stat = dealloc_error )
   REQUIRE( dealloc_error == 0 )

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

 1000 format( 5( 1PE16.8 ) )
 2000 format( 6I12 )
 3000 format( 4I12 )

   end subroutine evb_io


