! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Check EVB implementation against published hydroxypyridine/pyridone    |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#  include "../include/assert.fh"
#  include "../include/dprec.fh"

   subroutine schlegel_poh_irc

   use schlegel,  only: xdat_min, xdat_ts, ndg, ncoord, natm, b, rdgdim, ndg &
                      , Bohr2Angstrom
   use evb_math,  only: dist
   use evb_check, only: schlegel_debug
   use file_io_dat

   implicit none

#  include "parallel.h"
#  include "extra.h"

   !  ..........................................................................

   integer :: xn, io, ios, alloc_error
   integer :: j, k, m, n, nn, jend, kend, nroots
   intrinsic :: int
   _REAL_  :: qcart(3,natm,ndg), qcartTS(3,natm), qcartM1(3,natm) &
            , qcartM2(3,natm), qcartS(3,natm), dqM2_M1(3,natm), dqS_TS(3,natm)
   _REAL_  :: dha_dist(ndg,2), qtmp(natm*3)
   _REAL_  , allocatable :: roots(:,:,:), roots_del(:,:), rootsa(:), rootsb(:) &
                          , cart_grid(:,:,:,:), int_grid(:,:,:)
   _REAL_  , allocatable :: EVB_grid(:,:), dEVB_grid(:,:,:), ddEVB_grid(:,:,:,:)
   _REAL_  , allocatable :: V12SQ_grid(:,:), dV12SQ_grid(:,:,:) &
                          , ddV12SQ_grid(:,:,:,:) 
   _REAL_  :: lame_begin, lame_end, lame_step, lame_tol, jdel, kdel, odel, del &
            , f1, f2, y1, y2 
   _REAL_  :: vkl_sq, nrg_TS, nrg_M1, nrg_M2, nrg_tmp(3)
   _REAL_  :: ddvkl_sqDG(1+ncoord+ncoord*ncoord)
   intrinsic :: dble, sqrt 
   logical :: read_roots
   character(  1) :: ctype
   character( 40) :: label
   character(512) :: cline


   read_roots = .true.
   nroots     = 37

!  +---------------------------------------------------------------------------+
!  |  Read Cartesian coordinates from external files                           |
!  +---------------------------------------------------------------------------+

   io = schlegel_unit

 goto 888 
   if( master ) then
      do m = 1, 3
!        open( io, file = trim( adjustl( xdat_xdg(m)%filename ) ) )
         do
            read( io, '(A)', iostat = ios ) cline
!  .............................................................................
!  :  Check for eof and read errors                                            :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
            select case( ios )
               case(:-1)          ! end of file encountered
                  exit
               case(1:)           ! error during read
!                 write(6,'(A)') 'Error encountered while reading ' &
!                             //  trim( adjustl( xdat_xdg(m)%filename ) )
                  call mexit(6,1)
            end select
!  .............................................................................
!  :  Grab the ab initio energies                                              :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
            select case( trim( adjustl( cline(1:40) ) ) )
               case( 'Total Energy' )
                  backspace( io )
                  read( io, 200 ) label, ctype, nrg_tmp(m)
               case( 'Current cartesian coordinates' )
                  backspace( io )
                  read( io, 100 ) label, ctype, nn
                  read( io, 1000 ) ( qtmp(n), n = 1, nn )
                  qcart(:,:,m) = reshape( qtmp(:), (/ 3, natm /) )
                  exit
            end select
         enddo
         close( io )
      enddo
   endif

 888 continue

   nrg_TS = xdat_ts (1)%v
   nrg_M1 = xdat_min(1)%v
   nrg_M2 = xdat_min(2)%v

   qcartTS(:,:) = reshape( xdat_ts (1)%qcart(:), (/ 3, natm /) )
   qcartM1(:,:) = reshape( xdat_min(1)%qcart(:), (/ 3, natm /) )
   qcartM2(:,:) = reshape( xdat_min(2)%qcart(:), (/ 3, natm /) )

!  +---------------------------------------------------------------------------+
!  |  For pyridone, use the O7-H8 and N1-H8 distances                          |
!  +---------------------------------------------------------------------------+

   do n = 1, ndg
      dha_dist(n,1) = dist( qcart(:,:,n), 7, 8 )   
      dha_dist(n,2) = dist( qcart(:,:,n), 8, 1 )   
   enddo

!  +---------------------------------------------------------------------------+
!  |  Construct additional points by displacing H8 along the C1-H8 vector      |
!  +---------------------------------------------------------------------------+

   do n = 1, natm
      if( n == 8 ) then
         qcartS(:,n) = qcartTS(:,n) - 0.50d0 * ( qcartTS(:,8) - qcartTS(:,2) )
      else
         qcartS(:,n) = qcartTS(:,n)
      endif
   enddo

!  +---------------------------------------------------------------------------+
!  |  Determine the number of roots from the # of lines in rootsa.dat; this    |
!  |  overrides the default nroots defined at the top                          |
!  +---------------------------------------------------------------------------+

   if( read_roots .and. master ) then
      nroots = 0
      open( io, file = "rootsa.dat" )
      do
         nroots = nroots + 1
         read( io, '(A)', iostat = ios ) cline
!  .............................................................................
!  :  Check for eof and read errors                                            :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
         select case( ios )
            case(:-1)          ! end of file encountered
               nroots = nroots - 1
               nroots = int( sqrt( dble( nroots ) ) )
               exit
            case(1:)           ! error during read
               write(6,'(A)') 'Error encountered while reading rootsa.dat'
               call mexit(6,1)
         end select
      enddo
      close( io )
   endif

   write(6,'(A,2X,I8)') 'nroots = ', nroots

   allocate( roots(nroots,nroots,2), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( roots_del(nroots,nroots), stat = alloc_error ) 
   REQUIRE( alloc_error == 0 )
   allocate( rootsa(nroots*nroots), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( rootsb(nroots*nroots), stat = alloc_error ) 
   REQUIRE( alloc_error == 0 )

!  +---------------------------------------------------------------------------+
!  |  Choose a linear combination of O-H and N-H distances from a uniform      |
!  |  grid from 0.8 to 2.6 Angstrom                                            |
!  +---------------------------------------------------------------------------+

   dqM2_M1(:,:) = qcartM2(:,:) - qcartM1(:,:)
   dqS_TS(:,:)  =  qcartS(:,:) - qcartTS(:,:)

   if( read_roots ) then
      open( io, file = "rootsa.dat" )
      do n = 1, nroots * nroots
         read(io, * ) rootsa(n)
      enddo
      close( io )
      open( io, file = "rootsb.dat" )
      do n = 1, nroots * nroots
         read(io, * ) rootsb(n)
      enddo
      close( io )
      roots(:,:,1) = reshape( rootsa, (/ nroots, nroots /) )
      roots(:,:,2) = reshape( rootsb, (/ nroots, nroots /) )
   else
      lame_tol = 1.0d-3
      lame_begin = -1.5d0
      lame_end   =  1.5d0
      lame_step  = 1.0d-2
      kend = int(( lame_end - lame_begin ) / lame_step)
      jend = kend
      roots(:,:,:) = 0.0d0
      do n = 1, nroots
         y1 = ( 0.75d0 + 0.05d0 * dble(n) ) / Bohr2Angstrom
         do m = 1, nroots
            y2 = ( 0.75d0 + 0.050d0 * dble(m) ) / Bohr2Angstrom
            odel = 1.0d9
   lame:    do k = 1, kend
               kdel = lame_begin + dble(k) * lame_step
               do j = 1, jend
                  jdel = lame_begin + dble(j) * lame_step
                  f1 = dist( qcartTS(:,:) + jdel * dqM2_M1(:,:) &
                     + kdel * dqS_TS(:,:), 8, 7 )
                  f2 = dist( qcartTS(:,:) + jdel * dqM2_M1(:,:) &
                     + kdel * dqS_TS(:,:), 8, 1 )
                  del = abs( f1 - y1 ) + abs( f2 - y2 )
                  if( del < odel ) then
                     roots(m,n,1) = jdel
                     roots(m,n,2) = kdel
                     roots_del(m,n) = del 
                     odel = del
                     if( del < lame_tol ) exit lame
                  endif
               enddo
            enddo lame
         enddo
      enddo
   endif

!  +---------------------------------------------------------------------------+
!  |  Form the Cartesian grid points                                           |
!  +---------------------------------------------------------------------------+

   allocate( cart_grid(3,natm,nroots,nroots), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( int_grid(ncoord,nroots,nroots), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

   do n = 1, nroots
      do m = 1, nroots
         cart_grid(:,:,m,n) = qcartTS(:,:) + roots(m,n,1) * dqM2_M1(:,:) &
                                           + roots(m,n,2) *  dqS_TS(:,:)
      enddo
   enddo

!  +---------------------------------------------------------------------------+
!  |  Convert grid to redundant internals                                      |
!  +---------------------------------------------------------------------------+

   do n = 1, nroots
      do m = 1, nroots
         call cart2internal ( cart_grid(:,:,m,n), int_grid(:,m,n) )
      enddo
   enddo

   if( schlegel_debug ) then
      open ( io, file = "int_grid.debug" )
      write( io, '( A )' ) '[Data Size]'
      write( io, '( I20 )' ) nroots * nroots * ncoord
      write( io, '( A )' ) '[int_grid]'
      write( io, '(5(1PE16.8))' ) ( int_grid(:,:,:)  )
      close( io )
   endif

!  +---------------------------------------------------------------------------+
!  |  Form EVB surface                                                         |
!  +---------------------------------------------------------------------------+

   allocate( EVB_grid(nroots,nroots), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( dEVB_grid(ncoord,nroots,nroots), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( ddEVB_grid(ncoord,ncoord,nroots,nroots), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

   do n = 1, nroots
      do m = 1, nroots
         call schlegel_EVBR ( int_grid(:,m,n), vkl_sq, EVB_grid(m,n) &
                            , dEVB_grid(:,m,n), ddEVB_grid(:,:,m,n) )
      enddo
   enddo

   if( schlegel_debug ) then
      open ( io, file = "EVB_grid.debug" )
      write( io, '( A )' ) '[Data Size]'
      write( io, '( I20 )' ) nroots * nroots 
      write( io, '( A )' ) '[EVB_grid]'
      write( io, '(5(1PE16.8))' ) ( EVB_grid(:,n), n = 1, nroots  )
      close( io )

      open ( io, file = "dEVB_grid.debug" )
      write( io, '( A )' ) '[Data Size]'
      write( io, '( I20 )' ) ncoord * nroots * nroots 
      write( io, '( A )' ) '[dEVB_grid]'
      write( io, '(5(1PE16.8))' ) ( dEVB_grid(:,:,:)  )
      close( io )

      open ( io, file = "ddEVB_grid.debug" )
      write( io, '( A )' ) '[Data Size]'
      write( io, '( I20 )' ) ncoord * ncoord * nroots * nroots
      write( io, '( A )' ) '[ddEVB_grid]'
      write( io, '(5(1PE16.8))' ) ( ddEVB_grid(:,:,:,:)  )
      close( io )
   endif

!  +---------------------------------------------------------------------------+
!  |  Form V12 surface                                                         |
!  +---------------------------------------------------------------------------+

   allocate( V12SQ_grid(nroots,nroots), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( dV12SQ_grid(ncoord,nroots,nroots), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( ddV12SQ_grid(ncoord,ncoord,nroots,nroots), stat = alloc_error )

   do n = 1, nroots
      do m = 1, nroots
         call ddv12sqDG ( int_grid(:,m,n), b(1:rdgdim*ndg), ddvkl_sqDG )
           V12SQ_grid(    m,n) = ddvkl_sqDG(1)
          dV12SQ_grid(  :,m,n) = ddvkl_sqDG(2:ncoord+1)
         ddV12SQ_grid(:,:,m,n) = reshape( ddvkl_sqDG(ncoord+2:ncoord+1 &
                                   +ncoord*ncoord), (/ ncoord, ncoord /) )
      enddo
   enddo

   if( schlegel_debug ) then
      open ( io, file = "V12SQ_grid.debug" )
      write( io, '( A )' ) '[Data Size]'
      write( io, '( I20 )' ) nroots * nroots
      write( io, '( A )' ) '[V12SQ_grid]'
      write( io, '(5(1PE16.8))' ) ( V12SQ_grid(:,n), n = 1, nroots  )
      close( io )

      open ( io, file = "dV12SQ_grid.debug" )
      write( io, '( A )' ) '[Data Size]'
      write( io, '( I20 )' ) ncoord * nroots * nroots
      write( io, '( A )' ) '[dV12SQ_grid]'
      write( io, '(5(1PE16.8))' ) ( dV12SQ_grid(:,:,:)  )
      close( io )

      open ( io, file = "ddV12SQ_grid.debug" )
      write( io, '( A )' ) '[Data Size]'
      write( io, '( I20 )' ) ncoord * ncoord * nroots * nroots
      write( io, '( A )' ) '[ddV12SQ_grid]'
      write( io, '(5(1PE16.8))' ) ( ddV12SQ_grid(:,:,:,:)  )
      close( io )
   endif

!  +---------------------------------------------------------------------------+
!  |  Output data compatible with Mathematica ListContourPlot                  |
!  +---------------------------------------------------------------------------+

   if( worldrank == 0 ) then

      EVB_grid(:,:) = transpose( EVB_grid(:,:) )
      V12SQ_grid(:,:) = transpose( V12SQ_grid(:,:) )

      open( io, file = "EVB_poh_irc.nb" )

      write(io,'(A/)')
      write(io,'(A)') "Off[General::spell]"
      write(io,'(A/)')
      write(io,'(A)') ' "For visualization of published POH/PO EVB PES in'&
                   // ' Mathematica" '
      write(io,'(A/)')
      write(io,'(A)',advance="no") "EVBgrid = {"
      do xn = 1, nroots
         write(io,6000,advance="no") "{", EVB_grid(xn,1)
         write(io,7000,advance="no") EVB_grid(xn,2:nroots-1)
         write(io,8000,advance="no") EVB_grid(xn,nroots)
         if( xn == nroots ) then
            write(io,'(A)') "}};"
         else
            write(io,'(A)') "},"
         endif
      enddo
      write(io,'(A/)')
      write(io,'(A)') "MatrixForm[EVBgrid];"
      write(io,'(A/)')
      write(io,'(A,F14.8)') "nrgM1 = ", nrg_M1
      write(io,'(A,F14.8)') "nrgTS = ", nrg_TS
!  .............................................................................
!  :  2D EVB surface                                                           :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
      write(io,'(A/)')
      write(io,'(A)') 'pes2Dplot = ListContourPlot[EVBgrid, ' &
                   // 'ContourShading -> False, Contours -> 12, ' &
                   // 'PlotRange -> {nrgM1-0.01,nrgTS+0.108}, ' &
                   // 'FrameLabel -> {"R(O-H)", "R(N-H)"}, ' &
                   // 'FrameTicks -> {{{5,"1.0"}, {15,"1.5"}, {25, "2.0"}, ' &
                   // '{35,"2.5"}}, {{5,"1.0"},{15,"1.5"}, {25,"2.0"}, ' &
                   // '{35,"2.5"}}, {}, {}}]'
!  .............................................................................
!  :  3D EVB surface                                                           :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
      write(io,'(A/)')
      write(io,'(A)') 'pes3Dplot = ListPlot3D[EVBgrid, ' &
                   // 'PlotRange -> {nrgM1-0.01,nrgTS+0.11}, ' &
                   // 'ClipFill -> Automatic, ViewPoint -> {0.8,-2.4,1.5}, ' &
                   // 'AspectRatio -> 1, ' &
                   // 'AxesEdge -> {Automatic, {1,-1}, Automatic}, ' &
                   // 'AxesLabel -> {"R(O-H)", "R(N-H)", "Energy   "}, ' &
                   // 'Ticks -> {{{5,"1.0"}, {15,"1.5"}, {25, "2.0"}, ' &
                   // '{35,"2.5"}}, {{5,"1.0"},{15,"1.5"}, {25,"2.0"}, ' &
                   // '{35,"2.5"}}, {{-322.7053,"0"}, {-322.6575,"30"}, ' &
                   // '{-322.6097,"60"},{-322.5619,"90"}}}]'
!  .............................................................................
!  :  2D V12 surface                                                           :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
      write(io,'(A/)')
      write(io,'(A)',advance="no") "V12SQgrid = {"
      do xn = 1, nroots
         write(io,6000,advance="no") "{", V12SQ_grid(xn,1)
         write(io,7000,advance="no") V12SQ_grid(xn,2:nroots-1)
         write(io,8000,advance="no") V12SQ_grid(xn,nroots)
         if( xn == nroots ) then
            write(io,'(A)') "}};"
         else
            write(io,'(A)') "},"
         endif
      enddo
      write(io,'(A/)')
      write(io,'(A)') 'vij2Dplot = ListContourPlot[V12SQgrid, ' &
                   // 'ContourShading -> False, Contours -> 12, ' &
                   // 'FrameLabel -> {"R(O-H)", "R(N-H)"}, ' &
                   // 'FrameTicks -> {{{5,"1.0"}, {15,"1.5"}, {25, "2.0"}, ' &
                   // '{35,"2.5"}}, {{5,"1.0"},{15,"1.5"}, {25,"2.0"}, ' &
                   // '{35,"2.5"}}, {}, {}}]'
!  .............................................................................
!  :  3D V12 surface                                                           :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
      write(io,'(A/)')
      write(io,'(A)') 'vij3Dplot = ListPlot3D[V12SQgrid, ' &
                   // 'ClipFill -> Automatic, ViewPoint -> {0.8,-2.4,1.5}, ' &
                   // 'AspectRatio -> 1, ' &
                   // 'AxesEdge -> {Automatic, {1,-1}, Automatic}, ' &
                   // 'AxesLabel -> {"R(O-H)", "R(N-H)", "V12SQ    "}, ' &
                   // 'Ticks -> {{{5,"1.0"}, {15,"1.5"}, {25, "2.0"}, ' &
                   // '{35,"2.5"}}, {{5,"1.0"},{15,"1.5"}, {25,"2.0"}, ' &
                   // '{35,"2.5"}}, Automatic}]'

      close( io )

   endif

  100 format( A40, 3X, A1, 5X, I12 )
  200 format( A40, 3X, A1, 5X, E22.15 )
 1000 format( 5( 1PE16.8 ) )
 6000 format( A, 2X, F14.8 )
 7000 format( ( 5(',', 2X, F14.8) ) )
 8000 format( ',', 2X, F14.8 )

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine schlegel_poh_irc





























