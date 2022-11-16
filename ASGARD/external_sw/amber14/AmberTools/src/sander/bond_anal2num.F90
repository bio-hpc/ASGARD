! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Check forces from bond harmonic umbrella term                          |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#include "../include/dprec.fh"

   subroutine bond_anal2num ( q, fbond_anal, ndim )

   use evb_parm,  only: nbias, k_umb, r0_umb, bond_RC
   use evb_check, only: debug_toler, deviation, fdeps

   implicit none

   integer, intent(in) :: ndim
   _REAL_ , intent(in) :: q(ndim), fbond_anal(ndim)

   !............................................................................

   integer :: m, n, idx, jdx 

   _REAL_  :: pdt, mdt, Ep, Em, rij, fdeps_inv, fdeps2 &
            , drSQ(3), dtSQ(3), rms(2)
   _REAL_  :: fbond_num(ndim)
   intrinsic :: sqrt

!  +---------------------------------------------------------------------------+
!  |  Analytical                                                               |
!  +---------------------------------------------------------------------------+
!  |  /\ = R_ij                                                                |
!  |                                                                           |
!  |  V = 0.5 k_evb * ( /\ - eta )^2                                           |
!  |   F_Ri = - k_evb * ( /\ - eta ) / R_ij * ( R_i - R_j )                    |
!  |   F_Rj = + k_evb * ( /\ - eta ) / R_ij * ( R_i - R_j )                    |
!  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  |  Numerical                                                                |
!  +---------------------------------------------------------------------------+
!  |  F_x = -[ E(x+e/2) - E(x-e/2) ] / e                                       |
!  |  F_y = -[ E(y+e/2) - E(y-e/2) ] / e                                       |
!  |  F_z = -[ E(z+e/2) - E(z-e/2) ] / e                                       |
!  +---------------------------------------------------------------------------+

   fdeps2 = fdeps * 0.50d0
   fdeps_inv = 1.0d0 / fdeps 
   fbond_num(:) = 0.0d0 

   do n = 1, nbias

      idx = ( bond_RC(n)%iatom - 1 ) * 3
      jdx = ( bond_RC(n)%jatom - 1 ) * 3

      drSQ(1) = ( q(idx+1) - q(jdx+1) )**2
      drSQ(2) = ( q(idx+2) - q(jdx+2) )**2
      drSQ(3) = ( q(idx+3) - q(jdx+3) )**2

!  +---------------------------------------------------------------------------+
!  |  /\ = R_ij                                                                |
!  +---------------------------------------------------------------------------+

      do m = 1, 3

         dtSQ(:) = drSQ(:) 

         pdt = q(idx+m) + fdeps2 - q(jdx+m)
         dtSQ(m) = pdt * pdt
         rij = sqrt( dtSQ(1) + dtSQ(2) + dtSQ(3) )
         Ep = 0.50d0 * k_umb(n) * ( rij - r0_umb(n) )**2
         
         mdt = q(idx+m) - fdeps2 - q(jdx+m)
         dtSQ(m) = mdt * mdt
         rij = sqrt( dtSQ(1) + dtSQ(2) + dtSQ(3) )
         Em = 0.50d0 * k_umb(n) * ( rij - r0_umb(n) )**2

         fbond_num(idx+m) = fbond_num(idx+m) - ( Ep - Em ) * fdeps_inv

      enddo 

      do m = 1, 3

         dtSQ(:) = drSQ(:)

         pdt = q(idx+m) - ( q(jdx+m) + fdeps2 )
         dtSQ(m) = pdt * pdt
         rij = sqrt( dtSQ(1) + dtSQ(2) + dtSQ(3) )
         Ep = 0.50d0 * k_umb(n) * ( rij - r0_umb(n) )**2

         mdt = q(idx+m) - ( q(jdx+m) - fdeps2 )
         dtSQ(m) = mdt * mdt
         rij = sqrt( dtSQ(1) + dtSQ(2) + dtSQ(3) )
         Em = 0.50d0 * k_umb(n) * ( rij - r0_umb(n) )**2

         fbond_num(jdx+m) = fbond_num(jdx+m) - ( Ep - Em ) * fdeps_inv

      enddo

   enddo

!  +---------------------------------------------------------------------------+
!  |  Compare numerical and analytical forces                                  |
!  +---------------------------------------------------------------------------+

   write(6,'(2(/))')
   write(6,'(A)') '|^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^' &
               // '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^|'
   write(6,'(A)') '| +++++++ Checking analytical forces from bond ' &
               // 'potential vs. numerical +++++++ |'
   write(6,'(A)') '|,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,' &
               // ',,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,|'

   rms(:) = deviation( fbond_anal, fbond_num )
   if( rms(1) > debug_toler ) then
      write(6,'(/)')
      write(6,'(A)') '|ERROR: bond forces  >> numerical result ' &
                  // 'differs from the analytical'
      write(6,1000) 'FM_A       = ', fbond_anal
      write(6,1000) 'FM_#       = ', fbond_num
      write(6,4000) '|Deviation; largest component = ', rms(:)
   else if( rms(1) == 0.0d0 ) then
      write(6,'(/)')
      write(6,'(A)') '|WARNING: bond forces  >> numerical result ' &
                  // 'equals EXACTLY the analytical'
      write(6,1000) 'FM_A       = ', fbond_anal
      write(6,1000) 'FM_#       = ', fbond_num
   endif
 
 1000 format( A/, (5(2X,F14.8)) )
 4000 format( A ,  2(2X,E14.8)  )

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine bond_anal2num


