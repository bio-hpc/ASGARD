! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Check forces from dbonds harmonic umbrella term                        |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#include "../include/dprec.fh"

   subroutine dbonds_anal2num ( q, fdbonds_anal, ndim )

   use evb_parm,  only: nbias, k_umb, r0_umb, dbonds_RC
   use evb_check, only: debug_toler, deviation, fdeps

   implicit none

   integer, intent(in) :: ndim
   _REAL_ , intent(in) :: q(ndim), fdbonds_anal(ndim)

   !............................................................................

   integer :: m, n, nn, idx, jdx, kdx
   _REAL_  :: pdt, mdt, Ep, Em, rij, rkj, RC, fdeps_inv, fdeps2 &
            , drSQ(3), srSQ(3), dtSQ(3), rms(2)
   _REAL_  :: fdbonds_num(ndim)

!  +---------------------------------------------------------------------------+
!  |  Analytical                                                               |
!  +---------------------------------------------------------------------------+
!  |  /\ = R_ij - R_kj                                                         |
!  |                                                                           |
!  |  V = 0.5 k_evb * ( /\ - eta )^2                                           |
!  |   F_Ri = - k_evb * ( /\ - eta ) / R_ij * ( R_i - R_j )                    |
!  |   F_Rk = + k_evb * ( /\ - eta ) / R_kj * ( R_k - R_j )                    |
!  |   F_Rj = + k_evb * ( /\ - eta ) / R_ij * ( R_i - R_j )                    |
!  |          - k_evb * ( /\ - eta ) / R_kj * ( R_k - R_j )                    |
!  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  |  Numerical                                                                |
!  +---------------------------------------------------------------------------+
!  |  F_x = -[ E(x+e/2) - E(x-e/2) ] / e                                       |
!  |  F_y = -[ E(y+e/2) - E(y-e/2) ] / e                                       |
!  |  F_z = -[ E(z+e/2) - E(z-e/2) ] / e                                       |
!  +---------------------------------------------------------------------------+

   fdeps2 = fdeps * 0.50d0
   fdeps_inv = 1.0d0 / fdeps
   fdbonds_num(:) = 0.0d0 

   do n = 1, nbias

      idx = ( dbonds_RC(n)%iatom - 1 ) * 3
      jdx = ( dbonds_RC(n)%jatom - 1 ) * 3
      kdx = ( dbonds_RC(n)%katom - 1 ) * 3

      do nn = 1, 3
         drSQ(nn) = ( q(idx+nn) - q(jdx+nn) )**2
         srSQ(nn) = ( q(kdx+nn) - q(jdx+nn) )**2
      enddo

!  +---------------------------------------------------------------------------+
!  |  /\ = R_ij - R_kj    (operate derivative on first term)                   |
!  +---------------------------------------------------------------------------+

      rkj = sqrt( srSQ(1) + srSQ(2) + srSQ(3) )

      do m = 1, 3

         dtSQ(:) = drSQ(:) 

         pdt = q(idx+m) + fdeps2 - q(jdx+m)
         dtSQ(m) = pdt * pdt
         rij = sqrt( dtSQ(1) + dtSQ(2) + dtSQ(3) )
         RC = rij - rkj
         Ep = 0.50d0 * k_umb(n) * ( RC - r0_umb(n) )**2

         mdt = q(idx+m) - fdeps2 - q(jdx+m)
         dtSQ(m) = mdt * mdt
         rij = sqrt( dtSQ(1) + dtSQ(2) + dtSQ(3) )
         RC = rij - rkj
         Em = 0.50d0 * k_umb(n) * ( RC - r0_umb(n) )**2

         fdbonds_num(idx+m) = fdbonds_num(idx+m) - ( Ep - Em ) * fdeps_inv

      enddo 

      do m = 1, 3

         dtSQ(:) = drSQ(:)

         pdt = q(idx+m) - ( q(jdx+m) + fdeps2 )
         dtSQ(m) = pdt * pdt
         rij = sqrt( dtSQ(1) + dtSQ(2) + dtSQ(3) )
         RC = rij - rkj
         Ep = 0.50d0 * k_umb(n) * ( RC - r0_umb(n) )**2

         mdt = q(idx+m) - ( q(jdx+m) - fdeps2 )
         dtSQ(m) = mdt * mdt
         rij = sqrt( dtSQ(1) + dtSQ(2) + dtSQ(3) )
         RC = rij - rkj
         Em = 0.50d0 * k_umb(n) * ( RC - r0_umb(n) )**2

         fdbonds_num(jdx+m) = fdbonds_num(jdx+m) - ( Ep - Em ) * fdeps_inv

      enddo

!  +---------------------------------------------------------------------------+
!  |  /\ = R_ij - R_kj    (operate derivative on second term)                  |
!  +---------------------------------------------------------------------------+

      rij = sqrt( drSQ(1) + drSQ(2) + drSQ(3) )

      do m = 1, 3

         dtSQ(:) = srSQ(:)

         pdt = q(kdx+m) + fdeps2 - q(jdx+m)
         dtSQ(m) = pdt * pdt
         rkj = sqrt( dtSQ(1) + dtSQ(2) + dtSQ(3) )
         RC = rij - rkj
         Ep = 0.50d0 * k_umb(n) * ( RC - r0_umb(n) )**2

         mdt = q(kdx+m) - fdeps2 - q(jdx+m)
         dtSQ(m) = mdt * mdt
         rkj = sqrt( dtSQ(1) + dtSQ(2) + dtSQ(3) )
         RC = rij - rkj
         Em = 0.50d0 * k_umb(n) * ( RC - r0_umb(n) )**2

         fdbonds_num(kdx+m) = fdbonds_num(kdx+m) - ( Ep - Em ) * fdeps_inv

      enddo

      do m = 1, 3

         dtSQ(:) = srSQ(:)

         pdt = q(kdx+m) - ( q(jdx+m) + fdeps2 )
         dtSQ(m) = pdt * pdt
         rkj = sqrt( dtSQ(1) + dtSQ(2) + dtSQ(3) )
         RC = rij - rkj
         Ep = 0.50d0 * k_umb(n) * ( RC - r0_umb(n) )**2

         mdt = q(kdx+m) - ( q(jdx+m) - fdeps2 )
         dtSQ(m) = mdt * mdt
         rkj = sqrt( dtSQ(1) + dtSQ(2) + dtSQ(3) )
         RC = rij - rkj
         Em = 0.50d0 * k_umb(n) * ( RC - r0_umb(n) )**2

         fdbonds_num(jdx+m) = fdbonds_num(jdx+m) - ( Ep - Em ) * fdeps_inv

      enddo

   enddo

!  +---------------------------------------------------------------------------+
!  |  Compare numerical and analytical forces                                  |
!  +---------------------------------------------------------------------------+

   write(6,'(2(/))')
   write(6,'(A)') '|^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^' &
               // '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^|'
   write(6,'(A)') '| ++++++ Checking analytical forces from dbonds' &
               // ' potential vs. numerical ++++++ |'
   write(6,'(A)') '|,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,' &
               // ',,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,|'

   rms(:) = deviation( fdbonds_anal, fdbonds_num )
   if( rms(1) > debug_toler ) then
      write(6,'(/)')
      write(6,'(A)') '|ERROR: dbonds forces  >> numerical result ' &
                  // 'differs from the analytical'
      write(6,1000) 'FM_A       = ', fdbonds_anal
      write(6,1000) 'FM_#       = ', fdbonds_num
      write(6,4000) '|Deviation; largest component = ', rms(:)
   else if( rms(1) == 0.0d0 ) then
      write(6,'(/)')
      write(6,'(A)') '|WARNING: dbonds forces  >> numerical result ' &
                  // 'equals EXACTLY the analytical'
      write(6,1000) 'FM_A       = ', fdbonds_anal
      write(6,1000) 'FM_#       = ', fdbonds_num
   endif
 
 1000 format( A/, (5(2X,F14.8)) )
 4000 format( A ,  2(2X,E14.8)  )

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine dbonds_anal2num


