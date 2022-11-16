! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Check Cartesian gradients from Schlegel-Sonnenberg (redundant          |#
! #|  internals) EVB potential against finite difference                     |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#include "../include/dprec.fh"

   subroutine dg_grad_anal2num ( qcart, dEVB_anal, natm ) 

   use schlegel,  only: ncoord
   use evb_check, only: debug_toler, deviation, fdeps
   implicit none

   integer, intent(in   ) :: natm 
   _REAL_ , intent(in   ) :: dEVB_anal(natm*3)
   _REAL_ , intent(inout) :: qcart(3,natm)

   !  ..........................................................................

   integer :: m, n, mdx
   _REAL_  :: Ep, Em, fdeps_inv, fdeps2, rms(2)
   _REAL_  :: V_EVB, dV_EVB(ncoord), ddV_EVB(ncoord*ncoord) &
            , q(ncoord)
   _REAL_  :: EVB_num(natm*3)

!  +---------------------------------------------------------------------------+
!  |  Analytical                                                               |
!  :...........................................................................:
!  |  [G] = [B]^t [g]                                                          |
!  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  |  Numerical  (Note: this is the gradient, not the force)                   |
!  :...........................................................................:
!  |  G_x = [ E(x+e/2) - E(x-e/2) ] / e                                        |
!  |  G_y = [ E(y+e/2) - E(y-e/2) ] / e                                        |
!  |  G_z = [ E(z+e/2) - E(z-e/2) ] / e                                        |
!  +---------------------------------------------------------------------------+

   fdeps2 = fdeps * 0.50d0 
   fdeps_inv = 1.0d0 / fdeps

   mdx = 0
   do n = 1, natm
      do m = 1, 3
         mdx = mdx + 1

         qcart(m,n) = qcart(m,n) + fdeps2 
         call cart2internal ( qcart, q )
         call schlegel_evb( q, V_EVB, dV_EVB, ddV_EVB )
         Ep = V_EVB

         qcart(m,n) = qcart(m,n) - fdeps
         call cart2internal ( qcart, q )
         call schlegel_evb( q, V_EVB, dV_EVB, ddV_EVB )
         Em = V_EVB

         EVB_num(mdx) = ( Ep - Em ) * fdeps_inv 
         qcart(m,n) = qcart(m,n) + fdeps2

      enddo 
   enddo

!  +---------------------------------------------------------------------------+
!  |  Compare numerical and analytical forces                                  |
!  +---------------------------------------------------------------------------+

   write(6,'(2(/))')
   write(6,'(A)') '|^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^' &
               // '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^|'
   write(6,'(A)') '| ++++++ Checking analytical forces from DG EVB' &
               // ' potential vs. numerical ++++++ |'
   write(6,'(A)') '|,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,' &
               // ',,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,|'

   rms(:) = deviation( dEVB_anal, EVB_num )
   write(6,*) 'rms(1), debug_toler = ', rms(1), debug_toler
   if( rms(1) > debug_toler ) then
      write(6,'(/)')
      write(6,'(A)') '|ERROR: DG EVB forces  >> numerical result ' &
                  // 'differs from the analytical'
      write(6,1000) 'FM_A       = ', dEVB_anal
      write(6,1000) 'FM_#       = ', EVB_num
      write(6,4000) '|Deviation; largest component = ', rms(:)
   else if( rms(1) == 0.0d0 ) then
      write(6,'(/)')
      write(6,'(A)') '|WARNING: DG EVB forces  >> numerical result ' &
                  // 'equals EXACTLY the analytical'
      write(6,1000) 'FM_A       = ', dEVB_anal
      write(6,1000) 'FM_#       = ', EVB_num
   endif

 1000 format( A/, (5(2X,F14.8)) )
 4000 format( A ,  2(2X,E14.8)  )

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine dg_grad_anal2num

