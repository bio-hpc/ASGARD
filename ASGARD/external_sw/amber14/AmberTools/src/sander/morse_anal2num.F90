! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Check forces from morse potential against finite difference            |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#include "../include/dprec.fh"

   subroutine morse_anal2num ( q, fmorse_anal, ndim ) 

   use evb_parm,  only: nmorse
   use evb_amber, only: morse
   use evb_check, only: debug_toler, deviation, fdeps

   implicit none

   integer, intent(in) :: ndim
   _REAL_ , intent(in) :: q(ndim), fmorse_anal(ndim)

   !............................................................................

   integer :: m, n, idx, jdx 
   _REAL_  :: a, D, r0, pdt, mdt, Ep, Em, rij, fdeps2, fdeps_inv &
            , drSQ(3), dtSQ(3), rms(2)
   _REAL_  :: fmorse_num(ndim)
   intrinsic :: sqrt, exp

!  +---------------------------------------------------------------------------+
!  |  Analytical                                                               |
!  +---------------------------------------------------------------------------+
!  |  V = D * { 1 - exp[ - alpha * ( R_ij - r0 ) ] }^2                         |
!  |  F_Ri = -2.0 * D * { 1 - exp[ - alpha * ( R_ij - r0 ) ] }                 |
!  |                * alpha * exp[ - alpha * ( R_ij - r0 ) ]                   |
!  |                                / R_ij * ( R_i - R_j )                     |
!  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  |  Numerical                                                                |
!  +---------------------------------------------------------------------------+
!  |  F_x = -[ E(x+e/2) - E(x-e/2) ] / e                                       |
!  |  F_y = -[ E(y+e/2) - E(y-e/2) ] / e                                       |
!  |  F_z = -[ E(z+e/2) - E(z-e/2) ] / e                                       |
!  +---------------------------------------------------------------------------+

   fdeps2 = fdeps * 0.50d0 
   fdeps_inv = 1.0d0 / fdeps
   fmorse_num(:) = 0.0d0

   do n = 1, nmorse

      D  = morse(n)%D
      a  = morse(n)%a
      r0 = morse(n)%r0

      idx = ( morse(n)%iatom - 1 ) * 3
      jdx = ( morse(n)%jatom - 1 ) * 3

      drSQ(1) = ( q(idx+1) - q(jdx+1) )**2
      drSQ(2) = ( q(idx+2) - q(jdx+2) )**2
      drSQ(3) = ( q(idx+3) - q(jdx+3) )**2

      do m = 1, 3

         dtSQ(:) = drSQ(:)

         pdt = q(idx+m) + fdeps2 - q(jdx+m)
         dtSQ(m) = pdt * pdt
         rij = sqrt( dtSQ(1) + dtSQ(2) + dtSQ(3) )
         Ep = D * ( 1.0d0 - exp( - a * ( rij - r0 ) ) )**2

         mdt = q(idx+m) - fdeps2 - q(jdx+m)
         dtSQ(m) = mdt * mdt
         rij = sqrt( dtSQ(1) + dtSQ(2) + dtSQ(3) )
         Em = D * ( 1.0d0 - exp( - a * ( rij - r0 ) ) )**2

         fmorse_num(idx+m) = fmorse_num(idx+m) - ( Ep - Em ) * fdeps_inv

      enddo 

      do m = 1, 3

         dtSQ(:) = drSQ(:)

         pdt = q(idx+m) - ( q(jdx+m) + fdeps2 )
         dtSQ(m) = pdt * pdt
         rij = sqrt( dtSQ(1) + dtSQ(2) + dtSQ(3) )
         Ep = D * ( 1.0d0 - exp( - a * ( rij - r0 ) ) )**2

         mdt = q(idx+m) - ( q(jdx+m) - fdeps2 )
         dtSQ(m) = mdt * mdt
         rij = sqrt( dtSQ(1) + dtSQ(2) + dtSQ(3) )
         Em = D * ( 1.0d0 - exp( - a * ( rij - r0 ) ) )**2

         fmorse_num(jdx+m) = fmorse_num(jdx+m) - ( Ep - Em ) * fdeps_inv

      enddo

   enddo

!  +---------------------------------------------------------------------------+
!  |  Compare numerical and analytical forces                                  |
!  +---------------------------------------------------------------------------+

   write(6,'(2(/))')
   write(6,'(A)') '|^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^' &
               // '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^|'
   write(6,'(A)') '| ++++++ Checking analytical forces from Morse ' &
               // 'potential vs. numerical +++++++ |'
   write(6,'(A)') '|,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,' &
               // ',,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,|'

   rms(:) = deviation( fmorse_anal, fmorse_num )
   if( rms(1) > debug_toler ) then
      write(6,'(/)')
      write(6,'(A)') '|ERROR: Morse force  >> numerical result ' &
                  // 'differs from the analytical'
      write(6,1000) 'FM_A       = ', fmorse_anal
      write(6,1000) 'FM_#       = ', fmorse_num
      write(6,4000) '|Deviation; largest component = ', rms(:)
   else if( rms(1) == 0.0d0 ) then
      write(6,'(/)')
      write(6,'(A)') '|WARNING: Morse force  >> numerical result ' &
                  // 'equals EXACTLY the analytical'
      write(6,1000) 'FM_A       = ', fmorse_anal
      write(6,1000) 'FM_#       = ', fmorse_num
   endif
 
 1000 format( A/, (5(2X,F14.8)) )
 4000 format( A ,  2(2X,E14.8)  )

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine morse_anal2num


