! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Check gradient of gaussian coupling against finite difference          |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#include "../include/dprec.fh"

   subroutine xgauss_anal2num ( q, ndim ) 

   use evb_parm,  only: nxch, xch_gaussdat
   use evb_xchff, only: xch_gauss
   use evb_check, only: debug_toler, deviation, fdeps

   implicit none

   integer, intent(in) :: ndim
   _REAL_ , intent(in) :: q(ndim)

   !............................................................................

   integer :: m, n, idx, jdx 
   _REAL_  :: A, sigma, r0, pdt, mdt, Ep, Em, rij, fdeps2, fdeps_inv &
            , drSQ(3), dtSQ(3), rms(2)
   _REAL_  :: fgauss_anal(ndim), fgauss_num(ndim)
   intrinsic :: sqrt, exp

!  +---------------------------------------------------------------------------+
!  |  Analytical                                                               |
!  +---------------------------------------------------------------------------+
!  |  V_kl = A_kl * exp [ - ( R_ij - R_ij^(0,kl) )^2 / sigma^2 ]               |
!  |  F_Ri = V_kl * ( - 2.0/sigma^2 * ( R_ij - R_ij^(0,kl) )                   |
!  |       / R_ij * ( R_i - R_j )                                              |
!  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  |  Numerical (here, we just want the gradient ... not the force             |
!  |  ... hence, the positive sign in the finite diff)                         |
!  +---------------------------------------------------------------------------+
!  |  G_x = +[ E(x+e/2) - E(x-e/2) ] / e                                       |
!  |  G_y = +[ E(y+e/2) - E(y-e/2) ] / e                                       |
!  |  G_z = +[ E(z+e/2) - E(z-e/2) ] / e                                       |
!  +---------------------------------------------------------------------------+

   fdeps2 = fdeps * 0.50d0 
   fdeps_inv = 1.0d0 / fdeps
   fgauss_num(:) = 0.0d0

   do n = 1, nxch

      fgauss_anal(:) = xch_gauss(n)%gxch(:)

      A      = xch_gaussdat(n)%A
      sigma  = xch_gaussdat(n)%sigma
      r0     = xch_gaussdat(n)%r0

      idx = ( xch_gaussdat(n)%iatom - 1 ) * 3
      jdx = ( xch_gaussdat(n)%jatom - 1 ) * 3

      drSQ(1) = ( q(idx+1) - q(jdx+1) )**2
      drSQ(2) = ( q(idx+2) - q(jdx+2) )**2
      drSQ(3) = ( q(idx+3) - q(jdx+3) )**2

      do m = 1, 3

         dtSQ(:) = drSQ(:)

         pdt = q(idx+m) + fdeps2 - q(jdx+m)
         dtSQ(m) = pdt * pdt
         rij = sqrt( dtSQ(1) + dtSQ(2) + dtSQ(3) )
         Ep = A * exp( - ( rij - r0 )**2 / sigma**2 ) 

         mdt = q(idx+m) - fdeps2 - q(jdx+m)
         dtSQ(m) = mdt * mdt
         rij = sqrt( dtSQ(1) + dtSQ(2) + dtSQ(3) )
         Em = A * exp( - ( rij - r0 )**2 / sigma**2 ) 

         fgauss_num(idx+m) = fgauss_num(idx+m) + ( Ep - Em ) * fdeps_inv

      enddo 

      do m = 1, 3

         dtSQ(:) = drSQ(:)

         pdt = q(idx+m) - ( q(jdx+m) + fdeps2 )
         dtSQ(m) = pdt * pdt
         rij = sqrt( dtSQ(1) + dtSQ(2) + dtSQ(3) )
         Ep = A * exp( - ( rij - r0 )**2 / sigma**2 )

         mdt = q(idx+m) - ( q(jdx+m) - fdeps2 )
         dtSQ(m) = mdt * mdt
         rij = sqrt( dtSQ(1) + dtSQ(2) + dtSQ(3) )
         Em = A * exp( - ( rij - r0 )**2 / sigma**2 )

         fgauss_num(jdx+m) = fgauss_num(jdx+m) + ( Ep - Em ) * fdeps_inv

      enddo

!  +---------------------------------------------------------------------------+
!  |  Compare numerical and analytical gradients                               |
!  +---------------------------------------------------------------------------+

      write(6,'(2(/))')
      write(6,'(A)') '|^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^' &
                  // '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^|'
      write(6,'(A)') '| +++++ Checking analytical gradients from gaus' &
                  // 'sian coupling vs. numerical ++++|'
      write(6,'(A)') '|,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,' &
                  // ',,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,|'

      rms(:) = deviation( fgauss_anal, fgauss_num )
      if( rms(1) > debug_toler ) then
         write(6,'(/)')
         write(6,'(A)') '|ERROR: Gauss gradient  >> numerical result ' &
                     // 'differs from the analytical'
         write(6,1000) 'FM_A       = ', fgauss_anal
         write(6,1000) 'FM_#       = ', fgauss_num
         write(6,4000) '|Deviation; largest component = ', rms(:)
      else if( rms(1) == 0.0d0 ) then
         write(6,'(/)')
         write(6,'(A)') '|WARNING: Gauss gradient  >> numerical result ' &
                     // 'equals EXACTLY the analytical'
         write(6,1000) 'FM_A       = ', fgauss_anal
         write(6,1000) 'FM_#       = ', fgauss_num
      endif

   enddo
 
 1000 format( A/, (5(2X,F14.8)) )
 4000 format( A ,  2(2X,E14.8)  )

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine xgauss_anal2num


