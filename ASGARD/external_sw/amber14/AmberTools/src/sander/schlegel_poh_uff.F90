! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Check EVB implementation against published hydroxypyridine/pyridone    |#
! #|  system with the inclusion of UFF VDW interactions to V_ii              |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#  include "../include/assert.fh"
#  include "../include/dprec.fh"

   subroutine schlegel_poh_uff

   use schlegel,  only: xdat_xdg, ndg, ncoord, natm, b, rdgdim 
   use evb_check, only: schlegel_debug, ab_initio_gridfile
   use file_io_dat

   implicit none

#  include "parallel.h"
#  include "extra.h"

   !............................................................................

   integer :: io, ios, alloc_error
   integer :: k, m, n, mm, nn, xn, ntmp, ngrid_pts, ngrid 
   integer , allocatable :: ngeom_opt(:)
   intrinsic :: int, max
   _REAL_  , allocatable :: qcart_grid(:,:,:), qint_grid(:,:), energy_grid(:)
   _REAL_  , allocatable :: EVB_grid(:), dEVB_grid(:,:), ddEVB_grid(:,:,:)
   _REAL_  , allocatable :: V12SQ_grid(:), dV12SQ_grid(:,:), ddV12SQ_grid(:,:,:)
   _REAL_  , allocatable :: mEVB_grid(:,:), mV12SQ_grid(:,:), menergy_grid(:,:)
   _REAL_  , allocatable :: tmp(:)
   _REAL_  :: vkl_sq, nrg_TS, nrg_M1, nrg_M2, nrg_tmp(3)
   _REAL_  :: ddvkl_sqDG(1+ncoord+ncoord*ncoord) 
   intrinsic :: dble, sqrt 
   character(  1) :: ctype
   character( 40) :: label
   character(512) :: cline

!  +---------------------------------------------------------------------------+
!  |  Read energies from Gaussian fchk file for TS, RS, and PS (assumes        |
!  |  xdat_xdg(1:3) are in that order)                                         |
!  +---------------------------------------------------------------------------+

   io = schlegel_unit
   if( master ) then
      do m = 1, 3
         open( io, file = trim( adjustl( xdat_xdg(m)%filename ) ) )
         do
            read( io, '(A)', iostat = ios ) cline
!  .............................................................................
!  :  Check for eof and read errors                                            :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
            select case( ios )
               case(:-1)          ! end of file encountered
                  exit
               case(1:)           ! error during read
                  write(6,'(A)') 'Error encountered while reading ' &
                              //  trim( adjustl( xdat_xdg(m)%filename ) )
                  call mexit(6,1)
            end select
!  .............................................................................
!  :  Grab the ab initio energies                                              :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
            select case( trim( adjustl( cline(1:40) ) ) )
               case( 'Total Energy' )
                  backspace( io )
                  read( io, 200 ) label, ctype, nrg_tmp(m)
                  exit
            end select
         enddo
         close( io )
      enddo
   endif

   nrg_TS = nrg_tmp(1)
   nrg_M1 = nrg_tmp(2)
   nrg_M2 = nrg_tmp(3)

!  +---------------------------------------------------------------------------+
!  |  Read grid of Cartesian coordinates from external file.  For pyridone,    |
!  |  the 2D scan uses the O7-H8 and N1-H8 distances as the progress variables |
!  +---------------------------------------------------------------------------+

   if( master ) then
      open( io, file = trim( adjustl( ab_initio_gridfile ) ) )
      do
         read( io, '(A)', iostat = ios ) cline
!  .............................................................................
!  :  Check for eof and read errors                                            :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
         select case( ios )
            case(:-1)          ! end of file encountered
               exit
            case(1:)           ! error during read
               write(6,'(A)') 'Error encountered while reading ' &
                            // trim( adjustl( ab_initio_gridfile ) )
               call mexit(6,1)
         end select
!  .............................................................................
!  :  First pass grabs the # of grid points and the # of geom opt per grid     :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
         select case( trim( adjustl( cline(1:40) ) ) )
            case( 'Optimization Number of geometries' )
               backspace( io )
               read( io, 100 ) label, ctype, ngrid_pts
               allocate( ngeom_opt(ngrid_pts), stat = alloc_error )
               REQUIRE( alloc_error == 0 )
               read( io, 2000 ) ( ngeom_opt(n), n = 1, ngrid_pts )
               ngrid = int( sqrt( dble( ngrid_pts ) ) )
               exit
         end select
      enddo

      ntmp = 0
      do n = 1, ngrid_pts
         ntmp = max( ntmp, ngeom_opt(n) )
      enddo
      ntmp = ntmp * 3 * natm
      allocate( tmp(ntmp), stat = alloc_error )
      REQUIRE( alloc_error == 0 )
      allocate( energy_grid(ngrid_pts), stat = alloc_error )
      REQUIRE( alloc_error == 0 )
      allocate( qcart_grid(3,natm,ngrid_pts), stat = alloc_error )
      REQUIRE( alloc_error == 0 )
      rewind( io )

      k = 0
      m = 0   
      do
         read( io, '(A)', iostat = ios ) cline
!  .............................................................................
!  :  Check for eof and read errors                                            :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
         select case( ios )
            case(:-1)          ! end of file encountered
               exit
            case(1:)           ! error during read
               write(6,'(A)') 'Error encountered while reading ' &
                           //  trim( adjustl( xdat_xdg(m)%filename ) )
               call mexit(6,1)
         end select
!  .............................................................................
!  :  2nd pass grabs the energies & Cartesian geometries for the grids         :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
         select case( trim( adjustl( cline(19:40) ) ) )
            case( 'Results for each geome' )
               backspace( io )
               read( io, 100 ) label, ctype, nn
               read( io, 1000 ) ( tmp(n), n = 1, nn )
               k = k + 1
               energy_grid(k) = tmp(nn-1)
            case( 'Geometries' )
               backspace( io )
               read( io, 100 ) label, ctype, nn 
               read( io, 1000 ) ( tmp(n), n = 1, nn )
               m = m + 1
               mm = ( ngeom_opt(m) - 1 ) * 3 * natm + 1
               qcart_grid(:,:,m) = reshape( tmp(mm:nn), (/ 3, natm /) )
         end select
      enddo
      close( io )
      if( m /= ngrid ) write(6,'(A)') "Warning:: the # of geometries read is " &
                    // "different from the number of specified grid points." 
   endif

!  +---------------------------------------------------------------------------+
!  |  Convert grid to redundant internals                                      |
!  +---------------------------------------------------------------------------+

   allocate( qint_grid(ncoord,ngrid_pts), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

   do m = 1, ngrid_pts
      call cart2internal ( qcart_grid(:,:,m), qint_grid(:,m) )
   enddo

   if( schlegel_debug ) then
      open ( io, file = "qint_grid.debug" )
      write( io, '( A )' ) '[Data Size]'
      write( io, '( I20 )' ) ncoord * ngrid_pts
      write( io, '( A )' ) '[qint_grid]'
      write( io, '(5(1PE16.8))' ) ( qint_grid(:,:)  )
      close( io )
   endif

!  +---------------------------------------------------------------------------+
!  |  Form EVB surface                                                         |
!  +---------------------------------------------------------------------------+

   allocate( EVB_grid(ngrid_pts), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( dEVB_grid(ncoord,ngrid_pts), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( ddEVB_grid(ncoord,ncoord,ngrid_pts), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
 
   do m = 1, ngrid_pts
      call schlegel_EVBR ( qint_grid(:,m), vkl_sq, EVB_grid(m), dEVB_grid(:,m) &
                         , ddEVB_grid(:,:,m) )
   enddo

   if( schlegel_debug ) then
      open ( io, file = "EVB_grid.debug" )
      write( io, '( A )' ) '[Data Size]'
      write( io, '( I20 )' ) ngrid_pts
      write( io, '( A )' ) '[EVB_grid]'
      write( io, '(5(1PE16.8))' ) ( EVB_grid(n), n = 1, ngrid_pts )
      close( io )

      open ( io, file = "dEVB_grid.debug" )
      write( io, '( A )' ) '[Data Size]'
      write( io, '( I20 )' ) ncoord * ngrid_pts
      write( io, '( A )' ) '[dEVB_grid]'
      write( io, '(5(1PE16.8))' ) ( dEVB_grid(:,:)  )
      close( io )

      open ( io, file = "ddEVB_grid.debug" )
      write( io, '( A )' ) '[Data Size]'
      write( io, '( I20 )' ) ncoord * ncoord * ngrid_pts
      write( io, '( A )' ) '[ddEVB_grid]'
      write( io, '(5(1PE16.8))' ) ( ddEVB_grid(:,:,:)  )
      close( io )
   endif

!  +---------------------------------------------------------------------------+
!  |  Form V12 surface                                                         |
!  +---------------------------------------------------------------------------+

   allocate( V12SQ_grid(ngrid_pts), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( dV12SQ_grid(ncoord,ngrid_pts), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( ddV12SQ_grid(ncoord,ncoord,ngrid_pts), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

   do m = 1, ngrid_pts
      call ddv12sqDG ( qint_grid(:,m), b(1:rdgdim*ndg), ddvkl_sqDG )
        V12SQ_grid(    m) = ddvkl_sqDG(1)
       dV12SQ_grid(  :,m) = ddvkl_sqDG(2:ncoord+1)
      ddV12SQ_grid(:,:,m) = reshape( ddvkl_sqDG(ncoord+2:ncoord+1 &
                              +ncoord*ncoord), (/ ncoord, ncoord /) )
   enddo

   if( schlegel_debug ) then
      open ( io, file = "V12SQ_grid.debug" )
      write( io, '( A )' ) '[Data Size]'
      write( io, '( I20 )' ) ngrid_pts
      write( io, '( A )' ) '[V12SQ_grid]'
      write( io, '(5(1PE16.8))' ) ( V12SQ_grid(n), n = 1, ngrid_pts )
      close( io )

      open ( io, file = "dV12SQ_grid.debug" )
      write( io, '( A )' ) '[Data Size]'
      write( io, '( I20 )' ) ncoord * ngrid_pts
      write( io, '( A )' ) '[dV12SQ_grid]'
      write( io, '(5(1PE16.8))' ) ( dV12SQ_grid(:,:)  )
      close( io )

      open ( io, file = "ddV12SQ_grid.debug" )
      write( io, '( A )' ) '[Data Size]'
      write( io, '( I20 )' ) ncoord * ncoord * ngrid_pts
      write( io, '( A )' ) '[ddV12SQ_grid]'
      write( io, '(5(1PE16.8))' ) ( ddV12SQ_grid(:,:,:)  )
      close( io )
   endif

!  +---------------------------------------------------------------------------+
!  |  Output data compatible with Mathematica ListContourPlot                  |
!  +---------------------------------------------------------------------------+

!  .............................................................................
!  :  PES from grid file                                                       :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   if( worldrank == 0 ) then

      allocate( menergy_grid(ngrid,ngrid), stat = alloc_error )
      REQUIRE( alloc_error == 0 )

      mm = 0
      do m = 1, ngrid
         do nn = 1, ngrid
            if( mod(m,2) /= 0 ) then
               n = nn
            else
               n = ngrid + 1 - nn
            endif
            mm = mm + 1
            menergy_grid(m,n) = energy_grid(mm)
         enddo
      enddo 

      open( io, file = "HFgrid_poh.nb" )

      write(io,'(A/)')
      write(io,'(A)') "Off[General::spell]"
      write(io,'(A/)')
      write(io,'(A)') ' "For visualization of published POH/PO EVB PES in'&
                   // ' Mathematica" '
      write(io,'(A/)')
      write(io,'(A)',advance="no") "HFgrid = {"
      do xn = 1, ngrid
         write(io,6000,advance="no") "{", menergy_grid(xn,1)
         write(io,7000,advance="no") menergy_grid(xn,2:ngrid-1)
         write(io,8000,advance="no") menergy_grid(xn,ngrid)
         if( xn == ngrid ) then
            write(io,'(A)') "}};"
         else
            write(io,'(A)') "},"
         endif
      enddo
      write(io,'(A/)')
      write(io,'(A)') 'MatrixForm[HFgrid];'
      write(io,'(A/)')
      write(io,'(A)') 'AbInitioPlotRange = {Min[HFgrid],Min[HFgrid] + 0.15};'
      write(io,'(A/)')
      write(io,'(A)') 'GridMeshRange = {{0.8,2.6}, {0.8,2.6}};'
      write(io,'(A/)')
      write(io,'(A)') 'Slices = 15;'
      write(io,'(A/)')
      write(io,'(A)') 'GridLabelX = "R(N-H) / Angstroms"; '
      write(io,'(A)') 'GridLabelY = "R(O-H) / Angstroms"; '
!  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
!  :  2D EVB surface                                                           :
!  `````````````````````````````````````````````````````````````````````````````
      write(io,'(A/)')
      write(io,'(A)') 'pes2Dplot = ListContourPlot[HFgrid, ' &
                   // 'ContourShading -> False, Contours -> Slices, ' &
                   // 'PlotRange -> AbInitioPlotRange, ' &
                   // 'MeshRange -> GridMeshRange, ' &
                   // 'FrameLabel -> {GridLabelX, GridLabelY} ] ' 
!  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
!  :  3D EVB surface                                                           :
!  `````````````````````````````````````````````````````````````````````````````
      write(io,'(A/)')
      write(io,'(A)') 'GridViewPoint = {0.8,-2.4,1.5};'
      write(io,'(A/)')
      write(io,'(A)') 'GridImageSize = 450;'
      write(io,'(A/)')
      write(io,'(A)') 'pes3Dplot = ListPlot3D[HFgrid, ' &
                   // 'PlotRange -> AbInitioPlotRange, ' &
                   // 'ClipFill -> Automatic, ViewPoint -> GridViewPoint, ' &
                   // 'AspectRatio -> 1, ' &
                   // 'AxesEdge -> {Automatic, {1,-1}, Automatic}, ' &
                   // 'AxesLabel -> {GridLabelX, GridLabelY, "Energy in au  "}, ' &
                   // 'MeshRange -> GridMeshRange, ' &
                   // 'ImageSize -> GridImageSize ]'
      close( io )

!  .............................................................................
!  :  EVB surface                                                              :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

      allocate( mEVB_grid(ngrid,ngrid), stat = alloc_error )
      REQUIRE( alloc_error == 0 )
      allocate( mV12SQ_grid(ngrid,ngrid), stat = alloc_error )
      REQUIRE( alloc_error == 0 )

      mm = 0
      do m = 1, ngrid
         do nn = 1, ngrid
            if( mod(m,2) /= 0 ) then
               n = nn
            else
               n = ngrid + 1 - nn
            endif
            mm = mm + 1
            mEVB_grid(m,n) = EVB_grid(mm)
            mV12SQ_grid(m,n) = V12SQ_grid(mm)
         enddo
      enddo

      open( io, file = "EVB_poh_uff.nb" )

      write(io,'(A/)')
      write(io,'(A)') "Off[General::spell]"
      write(io,'(A/)')
      write(io,'(A)') ' "For visualization of published POH/PO EVB PES in'&
                   // ' Mathematica" '
      write(io,'(A/)')
      write(io,'(A)',advance="no") "EVBgrid = {"
      do xn = 1, ngrid
         write(io,6000,advance="no") "{", mEVB_grid(xn,1)
         write(io,7000,advance="no") mEVB_grid(xn,2:ngrid-1)
         write(io,8000,advance="no") mEVB_grid(xn,ngrid)
         if( xn == ngrid ) then
            write(io,'(A)') "}};"
         else
            write(io,'(A)') "},"
         endif
      enddo
      write(io,'(A/)')
      write(io,'(A,F14.8)') "nrgM1 = ", nrg_M1
      write(io,'(A,F14.8)') "nrgM2 = ", nrg_M2
      write(io,'(A,F14.8)') "nrgTS = ", nrg_TS
      write(io,'(A)') 'GridZeroEnergy = Min[nrgM1,nrgM2]'
      write(io,'(A/)')
      write(io,'(A)') 'MatrixForm[EVBgrid];'
      write(io,'(A/)')
      write(io,'(A)') 'GridPlotRange = {0,nrgTS-GridZeroEnergy+0.3};'
      write(io,'(A/)')
      write(io,'(A)') 'GridMeshRange = {{0.8,2.6}, {0.8,2.6}};'
      write(io,'(A/)')
      write(io,'(A)') 'Slices = 15;'
      write(io,'(A/)')
      write(io,'(A)') 'GridLabelX = "R(N-H) / Angstroms"; '
      write(io,'(A)') 'GridLabelY = "R(O-H) / Angstroms"; '
!  .............................................................................
!  :  2D EVB surface                                                           :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
      write(io,'(A/)')
      write(io,'(A)') 'pes2Dplot = ListContourPlot[EVBgrid-GridZeroEnergy, ' &
                   // 'ContourShading -> False, Contours -> Slices, ' &
                   // 'PlotRange -> GridPlotRange, ' &
                   // 'FrameLabel -> {GridLabelX, GridLabelY}, ' &
                   // 'MeshRange -> GridMeshRange ] '
!  .............................................................................
!  :  3D EVB surface                                                           :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
      write(io,'(A/)')
      write(io,'(A)') 'GridViewPoint = {0.8,-2.4,1.5};'
      write(io,'(A/)')
      write(io,'(A)') 'GridImageSize = 450;'
      write(io,'(A/)')
      write(io,'(A)') 'pes3Dplot = ListPlot3D[EVBgrid-GridZeroEnergy, ' &
                   // 'PlotRange -> GridPlotRange, ' &
                   // 'ClipFill -> Automatic, ViewPoint -> GridViewPoint, ' &
                   // 'AspectRatio -> 1, ' &
                   // 'AxesEdge -> {Automatic, {1,-1}, Automatic}, ' &
                   // 'AxesLabel -> {GridLabelX, GridLabelY, "Energy in au  "}, ' &
                   // 'MeshRange -> GridMeshRange, ' &
                   // 'ImageSize -> GridImageSize ]'
!  .............................................................................
!  :  2D V12 surface                                                           :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
      write(io,'(A/)')
      write(io,'(A)',advance="no") "V12SQgrid = {"
      do xn = 1, ngrid
         write(io,6000,advance="no") "{", mV12SQ_grid(xn,1)
         write(io,7000,advance="no") mV12SQ_grid(xn,2:ngrid-1)
         write(io,8000,advance="no") mV12SQ_grid(xn,ngrid)
         if( xn == ngrid ) then
            write(io,'(A)') "}};"
         else
            write(io,'(A)') "},"
         endif
      enddo
      write(io,'(A/)')
      write(io,'(A)') 'vij2Dplot = ListContourPlot[V12SQgrid, ' &
                   // 'ContourShading -> False, Contours -> Slices, ' &
                   // 'PlotRange -> {Min[V12SQgrid], Max[V12SQgrid]}, ' &
                   // 'FrameLabel -> {GridLabelX,GridLabelY}, ' &
                   // 'MeshRange -> GridMeshRange ]' 
!  .............................................................................
!  :  3D V12 surface                                                           :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
      write(io,'(A/)')
      write(io,'(A)') 'vij3Dplot = ListPlot3D[V12SQgrid, ' &
                   // 'ClipFill -> Automatic, ViewPoint -> GridViewPoint, ' &
                   // 'PlotRange -> {Min[V12SQgrid], Max[V12SQgrid]}, ' &
                   // 'AspectRatio -> 1, ' &
                   // 'AxesEdge -> {Automatic, {1,-1}, Automatic}, ' &
                   // 'AxesLabel -> {GridLabelX, GridLabelY, "Energy in au^2  "}, ' &
                   // 'MeshRange -> GridMeshRange, ' &
                   // 'ImageSize -> GridImageSize ]'
      close( io )

   endif

  100 format( A40, 3X, A1, 5X, I12 )
  200 format( A40, 3X, A1, 5X, E22.15 )
 1000 format( 5( 1PE16.8 ) )
 2000 format( 6I12 )
 6000 format( A, 2X, F14.8 )
 7000 format( ( 5(',', 2X, F14.8) ) )
 8000 format( 2X, F14.8 )

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine schlegel_poh_uff





























