! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  EVB_2STDEBUG                                                           |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#include "../include/dprec.fh"

   subroutine evb_2stdebug ( evb_vec, xf, ndim )

   use evb_parm,  only: nevb 
   use evb_data,  only: evb_mat_type, evb_frc_type, evb_Hmat, evb_frc
   use evb_check, only: debug_toler, deviation

   implicit none

   integer, intent(in) :: ndim
   _REAL_ , intent(in) :: xf(ndim,nevb), evb_vec(nevb,nevb)

   !  ..........................................................................

   _REAL_  :: nrg_evb(1), evb_f(ndim), c_11, c_21, wnorm &
            , nrg_tmp(1), f_tmp(ndim), vec_tmp(2), rms(2)
   _REAL_  :: s1_kl(ndim), d1_kl(ndim), b0_kl, s0_kl, d0_kl

!  +---------------------------------------------------------------------------+
!  |  EVB lowest state potential                                               |
!  |                                                                           |
!  |  V_0 = 0.5 ( V_11 + V_22) - [ ( V_11 - V_22 )^2 / 4 + V_12^2 ]^0.5        |
!  |                                                                           |
!  |  e_1 = ( -V_12    V_11-V_0 )^T                                            |
!  |                                                                           |
!  |  V_0^1(:) = S_kl^1(:) - [ D_kl^0 D_kl^0 + B_kl ]^-0.5                     |
!  |           * [ D_kl^0 D_kl^1(:) + 0.5 B_kl^1(:)                            |
!  |                                                                           |
!  |  where                                                                    |
!  |                                                                           |
!  |  B_kl^1(:) = B_kl^0 [ A1_kl(:) - A2_kl(:,:) * dq_kl(:)                    |
!  +---------------------------------------------------------------------------+

   write(6,'(2(/))')
   write(6,'(A)') '|^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^' &
               // '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^|'
   write(6,'(A)') '| ++++++++ 2-state EVB@ checking numerical against ' &
               // 'analytical results ++++++++ |'
   write(6,'(A)') '|,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,' &
               // ',,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,|'
   write(6,3000) '| ++++++++ Flag [Error] if RMS > ', debug_toler

   b0_kl = evb_Hmat%evb_mat(1,2)**2

   s0_kl      =  0.50d0 * ( evb_Hmat%evb_mat(1,1) + evb_Hmat%evb_mat(2,2) )
   s1_kl(:)   = -0.50d0 * ( xf(:,1) + xf(:,2) )

   d0_kl      =  0.50d0 * ( evb_Hmat%evb_mat(1,1) - evb_Hmat%evb_mat(2,2) )
   d1_kl(:)   = -0.50d0 * ( xf(:,1) - xf(:,2) )

   nrg_evb = s0_kl - sqrt( d0_kl * d0_kl + b0_kl )

   evb_f(:) = -s1_kl(:) + ( d0_kl * d1_kl(:) + 0.50d0 * evb_Hmat%b1_kl(:) ) &
                        / sqrt( d0_kl * d0_kl + b0_kl )

   c_11 = - evb_Hmat%evb_mat(1,2)
   c_21 =   evb_Hmat%evb_mat(1,1) - evb_frc%evb_nrg

   wnorm = sqrt( c_11 * c_11 + c_21 * c_21 )

   c_11 = c_11 / wnorm
   c_21 = c_21 / wnorm

   vec_tmp(1) = c_11 
   vec_tmp(2) = c_21 
   rms(:) = deviation( abs( vec_tmp ), abs( evb_vec(:,1) ) )
   if( rms(1) > debug_toler ) then
      write(6,'(/)')
      write(6,'(A)') '|ERROR: 2-state EVB AMPLITUDE >> numerical result ' &
                  // 'differs from the analytical'
      write(6,1000) 'C_#       = ', evb_vec(:,1)
      write(6,1000) 'C_A       = ', vec_tmp(:  )
      write(6,4000) '|Deviation; largest component = ', rms(:)
   else if( rms(1) == 0.0d0 ) then 
      write(6,'(/)')
      write(6,'(A)') '|WARNING: 2-state EVB AMPLITUDE >> numerical result ' &
                  // 'equals EXACTLY the analytical'
      write(6,1000) 'C_#       = ', evb_vec(:,1)
      write(6,1000) 'C_A       = ', vec_tmp(:  )
   endif

   nrg_tmp = evb_frc%evb_nrg
   rms(:) = deviation( nrg_tmp, nrg_evb )
   if( rms(1) > debug_toler ) then
      write(6,'(/)')
      write(6,'(A)') '|ERROR: 2-state EVB ENERGY >> numerical result ' &
                  // 'differs from the analytical'
      write(6,2000) 'E_#       = ', nrg_tmp
      write(6,2000) 'E_A       = ', nrg_evb
      write(6,3000) '|Deviation                    = ', rms(1)
   else if( rms(1) == 0.0d0 ) then 
      write(6,'(/)')
      write(6,'(A)') '|WARNING: 2-state EVB ENERGY >> numerical result ' &
                  // 'equals EXACTLY the analytical'
      write(6,2000) 'E_#       = ', nrg_tmp
      write(6,2000) 'E_A       = ', nrg_evb
   endif 

   f_tmp(:) = evb_frc%evb_f(:)
   rms(:) = deviation( f_tmp, evb_f )
   if( rms(1) > debug_toler ) then
      write(6,'(/)')
      write(6,'(A)') '|ERROR: 2-state EVB FORCE  >> numerical result ' &
                  // 'differs from the analytical'
      write(6,1000) 'F_#       = ', f_tmp
      write(6,1000) 'F_A       = ', evb_f
      write(6,4000) '|Deviation; largest component = ', rms(:)
   else if( rms(1) == 0.0d0 ) then 
      write(6,'(/)')
      write(6,'(A)') '|WARNING: 2-state EVB FORCE  >> numerical result ' &
                  // 'equals EXACTLY the analytical'
      write(6,1000) 'F_#       = ', f_tmp
      write(6,1000) 'F_A       = ', evb_f
   endif

   write(6,'(2(/))')
   write(6,'(A)') '|"""""""""""""""""""""""""""""""""""""""""""""""""""' &
               // '"""""""""""""""""""""""""""|'
   write(6,'(A)') '| ------- Done checking 2-state EVB numerical  vs. ' &
               // 'analytical results -------- |'
   write(6,'(A)') '|vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv' &
               // 'vvvvvvvvvvvvvvvvvvvvvvvvvvv|'

 1000 format( A/, (5(2X,F14.8)) )
 2000 format( A ,    2X,F20.8   )
 3000 format( A ,    2X,E14.8   )
 4000 format( A ,  2(2X,E14.8)  )

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine evb_2stdebug

