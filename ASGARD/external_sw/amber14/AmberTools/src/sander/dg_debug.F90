! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Check EVB data against DG data                                         |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#  include "../include/dprec.fh"

   subroutine evb_vs_abinitio

   use schlegel,  only: ndg, ncoord, xdat_xdg
   use evb_check, only: deviation, evb_debug

   implicit none

   !  ..........................................................................

   integer :: n
   _REAL_  :: q(ncoord), grad(ncoord), hess(ncoord*ncoord), nrg(1)
   _REAL_  :: V_EVB(1), dV_EVB(ncoord), ddV_EVB(ncoord*ncoord)
   _REAL_  :: rms(2)
   intrinsic :: reshape

!  +---------------------------------------------------------------------------+
!  |  Compare EVB derived energy, gradient, and hessian to DG counterpart      |
!  +---------------------------------------------------------------------------+

   do n = 1, ndg
      q(:) = xdat_xdg(n)%q(:)
      call schlegel_evb( q, V_EVB, dV_EVB, ddV_EVB )
      nrg(1)  = xdat_xdg(n)%v
      grad(:) = xdat_xdg(n)%d(:)
      hess(:) = reshape( xdat_xdg(n)%k(:,:), (/ ncoord*ncoord /) )

      rms(:) = deviation( V_EVB, nrg )
      if( rms(1) > evb_debug%schlegel_toler ) then
         write(6,'(/)')
         write(6,'(A)') '|ERROR: EVB energy - ab initio energy > toler ' 
         write(6,1000) 'EVB       = ', V_EVB(1)
         write(6,1000) 'ab initio = ', nrg(1)
         write(6,4000) '|Deviation                    = ', rms(1)
      endif

      rms(:) = deviation( dV_EVB, grad )
      if( rms(1) > evb_debug%schlegel_toler ) then
         write(6,'(/)')
         write(6,'(A)') '|ERROR: EVB gradient - ab initio gradient > toler '
         write(6,1000) 'EVB       = ', dV_EVB(:)
         write(6,1000) 'ab initio = ', grad(:)
         write(6,4000) '|Deviation; largest component = ', rms(:)
      endif

      rms(:) = deviation( ddV_EVB, hess )
      if( rms(1) > evb_debug%schlegel_toler ) then
         write(6,'(/)')
         write(6,'(A)') '|ERROR: EVB hessian - ab initio hessian > toler '
         write(6,1000) 'EVB       = ', ddV_EVB(:)
         write(6,1000) 'ab initio = ', hess(:)
         write(6,4000) '|Deviation; largest component = ', rms(:)
      endif

   enddo

 1000 format( A/, (5(2X,F14.8)) )
 4000 format( A ,  2(2X,E14.8)  )

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine evb_vs_abinitio

