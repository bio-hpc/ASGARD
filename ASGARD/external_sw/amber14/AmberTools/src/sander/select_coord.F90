! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Select internal coordinates for DG representation of V_12              |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#include "../include/assert.fh"
#include "../include/dprec.fh"

   subroutine select_coord ( stype, stol, scoord, nselect )

   use schlegel, only: ncoord, ncart, nbond, nangle, ndihed, xdat_min, xdat_ts
   use evb_math, only: dqint

   implicit none

   integer, intent(out) :: scoord(ncoord), nselect
   _REAL_ , intent( in) :: stol
   character(512), intent(in) :: stype

   !............................................................................

   integer :: n 
   _REAL_  :: diff_RP(ncoord), diff_RTS(ncoord), diff_PTS(ncoord)

!  +---------------------------------------------------------------------------+
!  |  Assume the internal coordinates are ordered as bonds, angles, dihedrals  |
!  +---------------------------------------------------------------------------+

   nselect = 0
   scoord(:) = 0

   select case( trim( adjustl( stype ) ) )

      case( "all_coords" )
         do n = 1, ncoord
            scoord(n) = n
         enddo
         nselect = ncoord
      case( "bonds_only" )
         do n = 1, nbond
            scoord(n) = n
         enddo 
         nselect = nbond
      case( "no_dihedrals" )
         do n = 1, nbond + nangle
            scoord(n) = n
         enddo
         nselect = nbond + nangle
      case( "react-product" )
         diff_RP(:) = abs( dqint ( xdat_min(1)%q, xdat_min(2)%q, nbond &
                                 , nangle, ndihed, ncoord, ncart ) )
         do n = 1, ncoord
            if( diff_RP(n) > stol ) then
               nselect = nselect + 1
               scoord(nselect) = n
            endif         
         enddo

      case( "react-ts-product" )
         diff_RP(:) = abs( dqint ( xdat_min(1)%q, xdat_min(2)%q, nbond &
                                 , nangle, ndihed, ncoord, ncart ) )
         diff_RTS(:) = abs( dqint ( xdat_min(1)%q, xdat_ts(1)%q, nbond &
                                  , nangle, ndihed, ncoord, ncart ) )
         diff_PTS(:) = abs( dqint ( xdat_min(2)%q, xdat_ts(1)%q, nbond &
                                  , nangle, ndihed, ncoord, ncart ) )
         do n = 1, ncoord
            if( diff_RP(n)  > stol .and. &
                diff_RTS(n) > stol .and. diff_PTS(n) > stol ) then
               nselect = nselect + 1
               scoord(nselect) = n
            endif
         enddo

   end select

   write(6,'(A)') 'DG::  subspace = ' // trim( adjustl( stype ) )
   if( trim( adjustl( stype ) ) == 'react-product' ) &
      write(6,'(A,1PE16.8)') 'DG:: selection tolerance = ', stol
   write(6,'(A,2X,I8)') 'DG::  nselect = ', nselect
   write(6,'(A/,(5(2X,I8)))') 'DG:: scoord(:) = ', scoord(1:nselect)

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine select_coord

