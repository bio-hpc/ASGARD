! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  EVB_2STDEBUG                                                           |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#include "../include/dprec.fh"

   subroutine egap_umb_2stdebug ( xf, v0, ndim )

   use evb_parm,  only: nevb, nbias, bias_ndx, k_umb, r0_umb 
   use evb_data,  only: evb_Hmat, evb_frc
   use evb_check, only: debug_toler, deviation

   implicit none

   integer, intent(in) :: ndim
   _REAL_ , intent(in) :: xf(ndim,nevb), v0

   !  ..........................................................................

   integer :: n, ni, nf
   _REAL_  :: evb_f(ndim), f_tmp(ndim), V_umb(1), nrg_tmp(1) &
            , evb_del, umb_pot, dumb_dLAM, rms(2)
   intrinsic :: sqrt             


!  +---------------------------------------------------------------------------+
!  |  Sampling according to Case's umbrella biasing potential                  |
!  |                                                                           |
!  |  V_0 = 0.5 ( V_11 + V_22) - [ ( V_11 - V_22 )^2 / 4 + V_12^2 ]^0.5        |
!  |                                                                           |
!  |  /\ = V_i - V_f                                                           |
!  |                                                                           |
!  |  V_umb(eta) = V_0 + 0.5 k_evb * ( /\ - eta )^2                            |
!  |                                                                           |
!  |  dV_umb / dR = ( dV_umb / d/\ ) * ( d/\ / dR )                            |
!  |    = 0.50 * ( force_i + force_f ) - [ 0.25 * ( V_i - V_f )                |
!  |           * [ ( V_i - V_f )^2/4 + V_if^2 ]^-0.5                           |
!  |           + k_evb ( /\ - eta ) ] ( force_i - force_f )                    |
!  +---------------------------------------------------------------------------+

   write(6,'(2(/))')
   write(6,'(A)') '|^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^' &
               // '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^|'
   write(6,'(A)') '| ++++ 2-state NRG_GAP UMBRELLA@ checking numerical ' &
               // 'against analytical results ++++ |'
   write(6,'(A)') '|,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,' &
               // ',,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,|'
   write(6,3000) '| ++++++++ Flag [Error] if RMS > ', debug_toler

   do n = 1, nbias

      ni = bias_ndx(n,1)
      nf = bias_ndx(n,2)

      evb_del = evb_Hmat%evb_mat(ni,ni) - evb_Hmat%evb_mat(nf,nf)
      umb_pot = 0.50d0 * k_umb(n) * ( evb_del - r0_umb(n) )**2
      V_umb = v0 + umb_pot
      dumb_dLAM = -0.25d0 * evb_del / SQRT( 0.25d0 * evb_del**2 &
                + evb_Hmat%evb_mat(ni,nf)**2 ) + k_umb(n) * ( evb_del - r0_umb(n) )

      evb_f(:) = ( 0.50d0 + dumb_dLAM ) * xf(:,ni) &
               + ( 0.50d0 - dumb_dLAM ) * xf(:,nf)

   enddo

   nrg_tmp = evb_frc%evb_nrg
   rms(:) = deviation( nrg_tmp, V_umb )
   if( rms(1) > debug_toler ) then
      write(6,'(/)')
      write(6,'(A)') '|ERROR: 2-state NRG_GAP UMBRELLA ENERGY >> numerical result ' &
                  // 'differs from the analytical'
      write(6,2000) 'E_#       = ', nrg_tmp
      write(6,2000) 'E_A       = ', V_umb
      write(6,3000) '|Deviation                    = ', rms(1)
   else if( rms(1) == 0.0d0 ) then
      write(6,'(/)')
      write(6,'(A)') '|WARNING: 2-state NRG_GAP UMBRELLA ENERGY >> numerical result ' &
                  // 'equals EXACTLY the analytical'
      write(6,2000) 'E_#       = ', nrg_tmp
      write(6,2000) 'E_A       = ', V_umb
   endif

   f_tmp(:) = evb_frc%evb_f(:)
   rms(:) = deviation( f_tmp, evb_f )
   if( rms(1) > debug_toler ) then
      write(6,'(/)')
      write(6,'(A)') '|ERROR: 2-state NRG_GAP UMBRELLA FORCE  >> numerical result ' &
                  // 'differs from the analytical'
      write(6,1000) 'F_#       = ', f_tmp
      write(6,1000) 'F_A       = ', evb_f
      write(6,4000) '|Deviation; largest component = ', rms(:)
   else if( rms(1) == 0.0d0 ) then 
      write(6,'(/)')
      write(6,'(A)') '|WARNING: 2-state NRG_GAP UMBRELLA FORCE  >> numerical result ' &
                  // 'equals EXACTLY the analytical'
      write(6,1000) 'F_#       = ', f_tmp
      write(6,1000) 'F_A       = ', evb_f
   endif

   write(6,'(2(/))')
   write(6,'(A)') '|"""""""""""""""""""""""""""""""""""""""""""""""""""' &
               // '""""""""""""""""""""""""""""""""|'
   write(6,'(A)') '| ---- Done checking 2-state NRG_GAP UMBRELLA numerical ' &
               // 'vs. analytical results ----- |'
   write(6,'(A)') '|vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv' &
               // 'vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv|'

 1000 format( A/, (5(2X,F14.8)) )
 2000 format( A ,    2X,F14.8   )
 3000 format( A ,    2X,E14.8   )
 4000 format( A ,  2(2X,E14.8)  )

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine egap_umb_2stdebug

