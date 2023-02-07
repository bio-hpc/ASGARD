! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Derivative of Wilson B matrix via finite difference                    |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#include "../include/dprec.fh"

   subroutine bmat_grad ( qcart, bmat_grad_anal, natm, fdeps ) 

   use schlegel,  only: ncoord

   implicit none

   integer, intent(in   ) :: natm 
!  _REAL_ , intent(in   ) :: bmat_grad_anal(ncoord,natm*3,natm*3)
   _REAL_ , intent(in   ) :: fdeps
   _REAL_ , intent(inout) :: qcart(3,natm)

   !............................................................................

   integer :: m, n, mdx
   _REAL_  :: bmat_p(ncoord,natm*3), bmat_m(ncoord,natm*3), &
              fdeps_inv, fdeps2
   _REAL_  :: bmat_grad_num (ncoord,natm*3,natm*3) &
            , bmat_grad_anal(ncoord,natm*3,natm*3)

!  +---------------------------------------------------------------------------+
!  |  Analytical                                                               |
!  +---------------------------------------------------------------------------+
!  |  [F] = [B]^t [f]                                                          |
!  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  |  Numerical                                                                |
!  +---------------------------------------------------------------------------+
!  |  dB_x = [ B(x+e/2) - B(x-e/2) ] / e                                       |
!  |  dB_y = [ B(y+e/2) - B(y-e/2) ] / e                                       |
!  |  dB_z = [ B(z+e/2) - B(z-e/2) ] / e                                       |
!  +---------------------------------------------------------------------------+

   fdeps2 = fdeps * 0.50d0 
   fdeps_inv = 1.0d0 / fdeps

   mdx = 0
   do n = 1, natm
      do m = 1, 3
         mdx = mdx + 1

         qcart(m,n) = qcart(m,n) + fdeps2 
         call wdc_bmat ( qcart, bmat_p )

         qcart(m,n) = qcart(m,n) - fdeps
         call wdc_bmat ( qcart, bmat_m )

         bmat_grad_num(:,:,mdx) = ( bmat_p(:,:) - bmat_m(:,:) ) * fdeps_inv 
         qcart(m,n) = qcart(m,n) + fdeps2
      enddo 
   enddo

   bmat_grad_anal(:,:,:) = bmat_grad_num(:,:,:)

!  +---------------------------------------------------------------------------+
!  |  Compare numerical and analytical forces                                  |
!  +---------------------------------------------------------------------------+

!  write(6,'(2(/))')
!  write(6,'(A)') '|^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^' &
!              // '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^|'
!  write(6,'(A)') '| ++++++++ Checking analytical derivative of B ' &
!              // 'matrix vs. numerical ++++++++++ |'
!  write(6,'(A)') '|,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,' &
!              // ',,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,|'

!  rms(:) = deviation( bmat_grad_anal, bmat_grad_num )
!  if( rms(1) > debug_toler ) then
!     write(6,'(/)')
!     write(6,'(A)') '|ERROR: gradient of B  >> numerical result ' &
!                 // 'differs from the analytical'
!     write(6,1000) 'B_A       = ', bmat_grad_anal
!     write(6,1000) 'B_#       = ', bmat_grad_num
!     write(6,4000) '|Deviation; largest component = ', rms(:)
!  else if( rms(1) == 0.0d0 ) then
!     write(6,'(/)')
!     write(6,'(A)') '|WARNING: gradient of B  >> numerical result ' &
!                 // 'equals EXACTLY the analytical'
!     write(6,1000) 'B_A       = ', bmat_grad_anal
!     write(6,1000) 'B_#       = ', bmat_grad_num
!  endif

!1000 format( A/, (5(2X,F14.8)) )
!4000 format( A ,  2(2X,E14.8)  )

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine bmat_grad

