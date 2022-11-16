! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Check forces from Lennard-Jones potential against finite difference    |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#include "../include/dprec.fh"

   subroutine vdw_anal2num ( q, fvdw_anal, ndim ) 

   use parms,     only: cn1, cn2
   use evb_parm,  only: nmodvdw
   use evb_amber, only: vdw_mod, iac, ico
   use evb_check, only: debug_toler, deviation, fdeps

   implicit none

#  include "../include/memory.h"

   integer, intent(in) :: ndim
   _REAL_ , intent(in) :: q(ndim), fvdw_anal(ndim)

   !............................................................................

   integer :: m, n, i, j, idx, jdx, ic 
   _REAL_  :: pdt, mdt, Ep, Em, fdeps2, fdeps_inv, r6 &
            , drSQ(3), dtSQ(3), rms(2)
   _REAL_  :: fvdw_num(ndim)

!  +---------------------------------------------------------------------------+
!  |  Analytical                                                               |
!  +---------------------------------------------------------------------------+
!  |  V = A_ij / R_ij^12 - B_ij / R_ij^6                                       |
!  |  F_Ri = (12 * A_ij / R_ij^14 - 6 * B_ij / R_ij^8 )                        |
!  |       * ( R_i - R_j )                                                     |
!  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  |  Numerical                                                                |
!  +---------------------------------------------------------------------------+
!  |  F_x = -[ E(x+e/2) - E(x-e/2) ] / e                                       |
!  |  F_y = -[ E(y+e/2) - E(y-e/2) ] / e                                       |
!  |  F_z = -[ E(z+e/2) - E(z-e/2) ] / e                                       |
!  +---------------------------------------------------------------------------+

   fdeps2 = fdeps * 0.50d0 
   fdeps_inv = 1.0d0 / fdeps
   fvdw_num(:) = 0.0d0

   do n = 1, nmodvdw

      i = vdw_mod(n)%iatom
      j = vdw_mod(n)%jatom

      ic = ico( ntypes * ( iac(i) - 1 ) + iac(j) )

      idx = ( i - 1 ) * 3
      jdx = ( j - 1 ) * 3

      drSQ(1) = ( q(idx+1) - q(jdx+1) )**2
      drSQ(2) = ( q(idx+2) - q(jdx+2) )**2
      drSQ(3) = ( q(idx+3) - q(jdx+3) )**2

      do m = 1, 3

         dtSQ(:) = drSQ(:)
         pdt = q(idx+m) + fdeps2 - q(jdx+m)
         dtSQ(m) = pdt * pdt
         r6 = 1.0d0 / ( dtSQ(1) + dtSQ(2) + dtSQ(3) )**3
         Ep = cn1(ic) * r6 * r6 - cn2(ic) * r6

         mdt = q(idx+m) - fdeps2 - q(jdx+m)
         dtSQ(m) = mdt * mdt
         r6 = 1.0d0 / ( dtSQ(1) + dtSQ(2) + dtSQ(3) )**3
         Em = cn1(ic) * r6 * r6 - cn2(ic) * r6

         fvdw_num(idx+m) = fvdw_num(idx+m) - ( Ep - Em ) * fdeps_inv

      enddo 

      do m = 1, 3

         dtSQ(:) = drSQ(:)
         pdt = q(idx+m) - ( q(jdx+m) + fdeps2 )
         dtSQ(m) = pdt * pdt
         r6 = 1.0d0 / ( dtSQ(1) + dtSQ(2) + dtSQ(3) )**3
         Ep = cn1(ic) * r6 * r6 - cn2(ic) * r6

         mdt = q(idx+m) - ( q(jdx+m) - fdeps2 )
         dtSQ(m) = mdt * mdt
         r6 = 1.0d0 / ( dtSQ(1) + dtSQ(2) + dtSQ(3) )**3
         Em = cn1(ic) * r6 * r6 - cn2(ic) * r6

         fvdw_num(jdx+m) = fvdw_num(jdx+m) - ( Ep - Em ) * fdeps_inv

      enddo

   enddo

!  +---------------------------------------------------------------------------+
!  |  Compare numerical and analytical forces                                  |
!  +---------------------------------------------------------------------------+

   write(6,'(2(/))')
   write(6,'(A)') '|^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^' &
               // '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^|'
   write(6,'(A)') '| +++ Checking analytical forces from Lennard-J' &
               // 'ones potential vs. numerical +++|'
   write(6,'(A)') '|,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,' &
               // ',,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,|'

   rms(:) = deviation( fvdw_anal, fvdw_num )
   if( rms(1) > debug_toler ) then
      write(6,'(/)')
      write(6,'(A)') '|ERROR: Lennard-Jones forces  >> numerical result ' &
                  // 'differs from the analytical'
      write(6,1000) 'FM_A       = ', fvdw_anal
      write(6,1000) 'FM_#       = ', fvdw_num
      write(6,4000) '|Deviation; largest component = ', rms(:)
   else if( rms(1) == 0.0d0 ) then
      write(6,'(/)')
      write(6,'(A)') '|WARNING: Lennard-Jones forces  >> numerical result ' &
                  // 'equals EXACTLY the analytical'
      write(6,1000) 'FM_A       = ', fvdw_anal
      write(6,1000) 'FM_#       = ', fvdw_num
   endif
 
 1000 format( A/, (5(2X,F14.8)) )
 4000 format( A ,  2(2X,E14.8)  )

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine vdw_anal2num


