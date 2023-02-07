! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Replace specified Amber harmonic bond interaction with Morse potential |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#include "../include/dprec.fh"

   subroutine morsify ( q, nrg_bnd, f, ndim, nbead ) 

   use evb_parm,  only: nmorse 
   use evb_amber, only: morse, k_harm, r0_harm
   use les_data,  only: cnum
#ifdef DEBUG_EVB
   use evb_check, only: full_evb_debug, morsify_debug
#endif

   implicit none

#  include "parallel.h"
#  include "../include/memory.h"

   integer, intent(in   ) :: ndim, nbead
   _REAL_ , intent(in   ) :: q(ndim)
   _REAL_ , intent(out  ) :: nrg_bnd(nbead)
   _REAL_ , intent(inout) :: f(ndim)

   !...........................................................................

   integer :: n, i, j, idx, jdx, bead_dx 
   _REAL_  :: dx, dy, dz, fx, fy, fz, rij, ff
   _REAL_  :: vharm(nbead), vmorse(nbead), fharm(ndim), fmorse(ndim) &
            , a, D, r0, exp_part
   intrinsic :: sqrt, exp
   logical :: lmorsify_error

!  +---------------------------------------------------------------------------+
!  |  V = K * ( rij - r0 )^2    /* K is already scaled by 0.5 */               |
!  |  F_Ri = - 2.0 * K * ( rij - r0 ) / rij * ( R_i - R_j )                    |
!  +---------------------------------------------------------------------------+

   lmorsify_error = .false. 

   vharm  (:) = 0.0d0
   vmorse (:) = 0.0d0

   fharm (:) = 0.0d0
   fmorse(:) = 0.0d0

   do n = 1, nmorse

      i = morse(n)%iatom
      j = morse(n)%jatom

      idx = ( i - 1 ) * 3 
      jdx = ( j - 1 ) * 3

      dx = q(idx+1) - q(jdx+1) 
      dy = q(idx+2) - q(jdx+2) 
      dz = q(idx+3) - q(jdx+3) 

      rij = sqrt( dx**2 + dy**2 + dz**2 )

#ifdef LES
      bead_dx = cnum(i) 
      if( bead_dx /= cnum(j) ) then 
         if( bead_dx == 0 ) then
            bead_dx = cnum(j)
         else
            lmorsify_error = .true.
         endif
      endif

#ifdef DEBUG_EVB
      if( full_evb_debug .or. morsify_debug ) then
         write(6,'(A,4(I10,2X))') '         ', i, j, cnum(i), cnum(j)
      endif
#endif

#else
      bead_dx = 1
#endif

      vharm(bead_dx) = vharm(bead_dx) + k_harm(n) * ( rij - r0_harm(n) )**2

      ff = 2.0 * k_harm(n) * ( rij - r0_harm(n) ) / rij

      fx = ff * dx
      fy = ff * dy
      fz = ff * dz
 
      fharm(idx+1) = fharm(idx+1) - fx
      fharm(idx+2) = fharm(idx+2) - fy
      fharm(idx+3) = fharm(idx+3) - fz

      fharm(jdx+1) = fharm(jdx+1) + fx
      fharm(jdx+2) = fharm(jdx+2) + fy
      fharm(jdx+3) = fharm(jdx+3) + fz

!  +---------------------------------------------------------------------------+
!  |  V = D * { 1 - exp[ - alpha * ( rij - r0 ) ] }^2                          |
!  |  F_Ri = -2.0 * D * { 1 - exp[ - alpha * ( rij - r0 ) ] }                  |
!  |       * alpha * exp[ - alpha * ( rij - r0 ) ] / rij * ( R_i - R_j )       |
!  +---------------------------------------------------------------------------+

      D  = morse(n)%D
      a  = morse(n)%a
      r0 = morse(n)%r0

      exp_part = exp( - a * ( rij - r0 ) )
      vmorse(bead_dx) = vmorse(bead_dx) + D * ( 1.0d0 - exp_part )**2

      ff = 2.0 * D * ( 1.0d0 - exp_part ) * a * exp_part / rij

      fx = ff * dx
      fy = ff * dy
      fz = ff * dz

      fmorse(idx+1) = fmorse(idx+1) - fx
      fmorse(idx+2) = fmorse(idx+2) - fy
      fmorse(idx+3) = fmorse(idx+3) - fz

      fmorse(jdx+1) = fmorse(jdx+1) + fx
      fmorse(jdx+2) = fmorse(jdx+2) + fy
      fmorse(jdx+3) = fmorse(jdx+3) + fz

   enddo

   if( lmorsify_error ) then
      do n = 1, nmorse
         i = morse(n)%iatom
         j = morse(n)%jatom
         write(6,'(A,4(2X,I8))') 'morsify between atoms ', i, j, cnum(i), cnum(j)
      enddo
      call mexit(6,0) 
   endif

!  +---------------------------------------------------------------------------+
!  |  "Morsify" -- Bill Miller                                                 |
!  +---------------------------------------------------------------------------+

   f(:) = f(:) - fharm(:) + fmorse(:)
   nrg_bnd(:) = - vharm(:) + vmorse(:) 

#ifdef DEBUG_EVB
   if( full_evb_debug .or. morsify_debug ) call morse_anal2num ( q, fmorse, ndim )
#endif

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine morsify


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Read Amber harmonic bond parameters and prepare to morsify             |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   subroutine morsify_init ( ix ) 

   use parms,     only: rk, req
   use evb_parm,  only: nmorse
   use evb_amber, only: morse, k_harm, r0_harm, morsify_initialized

   implicit none

#  include "../include/memory.h"

   integer :: ix(*)

   !............................................................................

   integer :: m, n, ii, jj, idx, jdx
   integer :: istart, jstart, cstart 


!  +---------------------------------------------------------------------------+
!  |  Store the Amber harmonic bond parameters for EVB usage.  The index       |
!  |  pointers iibh, ijbh, iicbh point to the location of the ith and jth      |
!  |  atoms forming bond i with the associated harmonic constant and           |
!  |  equilibrium distance at location iicbh                                   |
!  +---------------------------------------------------------------------------+

    k_harm(:) = 9999.0d0
   r0_harm(:) = 9999.0d0 

!! '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
!! ;  Do bonds involving hydrogen                                              ;
!! ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

   istart =  iibh - 1
   jstart =  ijbh - 1
   cstart = iicbh - 1

   do n = 1, nbonh
      idx = ix( istart + n ) / 3 + 1
      jdx = ix( jstart + n ) / 3 + 1
      do m = 1, nmorse 
         ii = morse(m)%iatom
         jj = morse(m)%jatom
         if( idx == ii .and. jdx == jj .or. &
             idx == jj .and. jdx == ii        ) then 
             k_harm(m) =  rk( ix( cstart + n ) )
            r0_harm(m) = req( ix( cstart + n ) )
         endif
      enddo
   enddo 

!! '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
!! ;  Do all other bonds                                                       ;
!! ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

   istart =  iiba - 1
   jstart =  ijba - 1
   cstart = iicba - 1

   do n = 1, nbona
      idx = ix( istart + n ) / 3 + 1
      jdx = ix( jstart + n ) / 3 + 1
      do m = 1, nmorse
         ii = morse(m)%iatom 
         jj = morse(m)%jatom 
         if( idx == ii .and. jdx == jj .or. &
             idx == jj .and. jdx == ii        ) then
             k_harm(m) =  rk( ix( cstart + n ) )
            r0_harm(m) = req( ix( cstart + n ) )
         endif
      enddo
   enddo 

   do n = 1, nmorse
      if(  k_harm(n) == 9999.0d0 .or. &
          r0_harm(n) == 9999.0d0        ) then
         write(6,'(A)') 'ERROR: User requested morsification for the bond below'
         write(6,'(A)') '       but Amber harmonic parameters are not defined.'
         write(6,'(A,I8)' ) '       iatom = ', morse(n)%iatom
         write(6,'(A,I8)' ) '       jatom = ', morse(n)%jatom
         write(6,'(A)')
         write(6,'(A)') '... possible reason could be that these atoms do ' &
                     // 'not form a bond.'
         call mexit(6,1)
      endif
   enddo 

   morsify_initialized = .true. 

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine morsify_init


