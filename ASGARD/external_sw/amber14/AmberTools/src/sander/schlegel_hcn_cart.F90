! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Check EVB implementation against published HCN system                  |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#  include "../include/assert.fh"
#  include "../include/dprec.fh"

   subroutine schlegel_hcn_cart

   use schlegel, only: ncoord, xdat_min, Bohr2Angstrom
   use file_io_dat

   implicit none

#  include "parallel.h"
#  include "extra.h"

   !  ..........................................................................

   integer , parameter :: xndata = 60
   integer , parameter :: zndata = 60
   integer :: xn, zn, io, alloc_error
   _REAL_  , allocatable :: contour(:,:)
   _REAL_  :: xbegin, xend, xincr, zbegin, zend, zincr, toBohr
   _REAL_  :: V_EVB, dV_EVB(ncoord), ddV_EVB(ncoord*ncoord), q(ncoord), nrg_M1

!  +---------------------------------------------------------------------------+
!  |  Generate data for Mathematica visualization of PES contours              |
!  +---------------------------------------------------------------------------+

   io = schlegel_unit
   toBohr = 1.0d0 / Bohr2Angstrom
   nrg_M1 = xdat_min(1)%v

   allocate ( contour(xndata,zndata), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

   xbegin =  0.0d0 * toBohr
   xend   =  2.0d0 * toBohr
   zbegin = -1.4d0 * toBohr
   zend   =  2.6d0 * toBohr

   xincr = ( xend - xbegin ) / dble( xndata )
   zincr = ( zend - zbegin ) / dble( zndata )

   q(:) = 0.00d0 * toBohr
   q(6) = 1.16d0 * toBohr

   do xn = 1, xndata
      q(7) = xbegin + dble(xn) * xincr 
      do zn = 1, zndata
         q(9) = zbegin + dble(zn) * zincr 
         call schlegel_evb( q, V_EVB, dV_EVB, ddV_EVB )
         contour(xn,zn) = V_EVB
      enddo
   enddo

!  +---------------------------------------------------------------------------+
!  |  Output data compatible with Mathematica ListContourPlot                  |
!  +---------------------------------------------------------------------------+

   if( worldrank == 0 ) then
      open( io, file = "EVB_hcn_cart.nb" )

      write(io,'(A/)')
      write(io,'(A)') "Off[General::spell]"
      write(io,'(A/)')
      write(io,'(A)') ' "For visualization of published HCN EVB PES in Mathematica" '
      write(io,'(A/)')
      write(io,'(A)',advance="no") "pes = {"
      do xn = 1, xndata
         write(io,1000,advance="no") "{", contour(xn,1)
         write(io,2000,advance="no") contour(xn,2:zndata-1)
         write(io,3000,advance="no") contour(xn,zndata)
         if( xn == xndata ) then
            write(io,'(A)') "}};"
         else
            write(io,'(A)') "},"
         endif
      enddo
      write(io,'(A/)')
      write(io,'(A)') "MatrixForm[pes];"
      write(io,'(A/)')
      write(io,'(A,F14.8)') "energyMin1 = ", nrg_M1
      write(io,'(A/)')
      write(io,'(A)') "contourvals=0.031872{-2,-1,0,1,2,3,4,5,6,7,8,9,10};"
      write(io,'(A)') "pesplot = ListContourPlot[pes - energyMin1, Axes -> None, " &
                   // "Frame -> True, ContourShading -> False, " &
                   // "PlotRange -> {{1,60},{1,60},Automatic}, " &
                  // "Contours -> contourvals, AspectRatio -> .5, FrameTicks -> None]"

      close( io )
   endif

 1000 format( A, 2X, F14.8 )
 2000 format( ( 5(',', 2X, F14.8) ) )
 3000 format( 2X, F14.8 )

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine schlegel_hcn_cart


