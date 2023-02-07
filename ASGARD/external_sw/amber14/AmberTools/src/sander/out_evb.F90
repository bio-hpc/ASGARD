! <compile=optimized>
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  OUT_EVB                                                        |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#include "../include/dprec.fh"

   subroutine out_evb ( nstep )

   use evb_parm,  only: nbias, evb_dyn, ntw_evb, inc_dbonds_RC &
                      , out_RCdot
   use evb_data,  only: evb_bias, RCdot
   use file_io_dat
#ifdef LES
   use evb_data,  only: nrg_vel0, nrg_frc, G_QI, f_v, F_QI
   use pimd_vars, only: nbead, itimass
   use evb_pimd,  only: evb_mat_bead, evb_vec0_bead, nbead_inv
   use miller,    only: dlnQ_dl, dlnCdd_dl, i_qi, gradRC_norm, ndiv
   use evb_parm,  only: nevb
#else
   use evb_data,  only: evb_Hmat
#endif

   implicit none

#  include "parallel.h"
#  include "../include/md.h"

   integer, intent(in) :: nstep

   !  +---------------------------------------------------------------+

   integer :: n
#ifdef LES
   _REAL_  :: evb_matrix(nevb,nevb), evb_pop(nevb)
#endif

   if( worldrank == 0 ) then

!  +---------------------------------------------------------------+
!  |  Output EVB base info                                         |
!  +---------------------------------------------------------------+

      if( mod(nstep,ntw_evb) == 0 ) then

         write(evb_unit,'(/)')
         write(evb_unit,888) '{NSTEP}: ', nstep, '{TIME}: ', nstep*dt
#ifdef LES
         evb_matrix(:,:) = 0.0d0
         evb_pop(:) = 0.0d0
         do n = 1, nbead 
            evb_matrix(:,:) = evb_matrix(:,:) + evb_mat_bead(:,:,n)
            evb_pop(:) = evb_pop(:) + evb_vec0_bead(:,n)**2
         enddo

         write(evb_unit,'(A)')
         write(evb_unit,1000) '  [EVB MATRIX]', evb_matrix(:,:)
         write(evb_unit,'(A)')
         write(evb_unit,1000) '  [EVB VEC_0^2]', evb_pop(:) * nbead_inv

!FOR DEBUGGING
!        write(evb_unit,'(A)')
!        write(evb_unit,1000) '  [VEL_0 BEAD]', vel0_bead(:)
!FOR DEBUGGING

#else
         write(evb_unit,'(A)')
         write(evb_unit,1000) '  [EVB MATRIX]', evb_Hmat%evb_mat(:,:)
         write(evb_unit,'(A)')
         write(evb_unit,1000) '  [EVB VEC_0]', evb_Hmat%evb_vec0(:)
#endif

         select case( trim( adjustl( evb_dyn) ) )

            case( "groundstate" )
#ifdef LES
               if( inc_dbonds_RC ) then
                  do n = 1, nbias
                     write(evb_unit,'(A)')
                     write(evb_unit,1000) '  [RC EVB]', evb_bias%RC(n), evb_bias%RC(n)
                  enddo
               endif

               write(evb_unit,'(A)')
               write(evb_unit,999) '{VEL0_PIMD}: ', ( nrg_frc(n), n = 1, 3 )

#else
               if( inc_dbonds_RC ) then
                  do n = 1, nbias
                     write(evb_unit,'(A)')
                     write(evb_unit,1000) '  [RC EVB]', evb_bias%RC(n), evb_bias%RC(n)
                  enddo
               endif

#endif

!  +---------------------------------------------------------------+
!  |  Sampling according to Warshel's mapping potential            |
!  |                                                               |
!  |  /\ = V_i - V_f                                               |
!  |                                                               |
!  |  V_map(n) = V_i + lambda(n) * ( V_f - V_i )                   |
!  +---------------------------------------------------------------+

            case( "evb_map" )
#ifdef LES
               write(evb_unit,'(A)')
               write(evb_unit,999) '{EMAP_PIMD}: ', ( nrg_frc(n), n = 1, 3 )
               write(evb_unit,'(A)')
               write(evb_unit,999) '{VEL0_PIMD}: ', ( nrg_vel0(n), n = 1, 3 )
#else
               if( inc_dbonds_RC ) then
                  do n = 1, nbias
                     write(evb_unit,'(A)')
                     write(evb_unit,1000) '  [RC EVB]', evb_bias%RC(n), evb_bias%RC(n)
                  enddo
               endif
#endif 

!  +---------------------------------------------------------------+
!  |  Sampling according to Case's umbrella biasing potential      |
!  |                                                               |
!  |  /\ = V_i - V_f                                               |
!  |                                                               |
!  |  V_umb(eta) = V_0 + 0.5 k_evb * ( /\ - eta )^2                |
!  |                                                               |
!  |  dV_umb / dR = ( dV_umb / d/\ ) * ( d/\ / dR )                |
!  |    = dV_0 / dR + k_evb * ( /\ - eta ) * d/dR ( V_i - V_f )    |
!  +---------------------------------------------------------------+

            case( "egap_umb" )
#ifdef LES
               write(evb_unit,'(A)')
               write(evb_unit,999) '{VTOT_PIMD}: ', ( nrg_frc(n), n = 1, 3 )
               write(evb_unit,'(A)')
               write(evb_unit,999) '{VEL0_PIMD}: ', ( nrg_vel0(n), n = 1, 3 )
#endif 

!  +---------------------------------------------------------------+
!  |  Harmonic umbrella biased sampling                            |
!  |                                                               |
!  |  (1) /\ = r_ij - r_kj                                         |
!  |  (2) /\ = r_ij                                                |
!  +---------------------------------------------------------------+

            case( "dbonds_umb", "bond_umb",  &
                  "qi_dbonds_pmf", "qi_bond_pmf", "qi_dbonds_dyn", "qi_bond_dyn" )
               do n = 1, nbias
                  write(evb_unit,'(A)')
                  write(evb_unit,1000) '  [RC EVB]', evb_bias%RC(n), evb_bias%RC(n)
               enddo
#ifdef LES
               write(evb_unit,'(A)')
               write(evb_unit,999) '{VEL0_PIMD}: ', ( nrg_frc(n), n = 1, 3 )
#endif

         end select

!  +---------------------------------------------------------------+
!  |  Thermodynamic integration by mass                            |
!  +---------------------------------------------------------------+

#ifdef LES
         if( itimass > 0 ) then
            write(evb_unit,'(A)')
            write(evb_unit,999) '{TI MASS: (d/dl) V_eff}: ', dlnQ_dl

            if( i_qi > 0 ) then
               write(evb_unit,'(A)')
               write(evb_unit,999) '{TI MASS: -KT * (d/dl) C_dd}: ', dlnCdd_dl
            endif
         endif
#endif

!  +---------------------------------------------------------------+
!  |  TST dynamical frequency factor                               |
!  +---------------------------------------------------------------+

         if( out_RCdot ) then
            write(evb_unit,'(A)')
            write(evb_unit,1000) '{TST: (d/dt) RC}: ', RCdot
         endif

!  +---------------------------------------------------------------+
!  |  QI dynamical factors                                         |
!  +---------------------------------------------------------------+

#ifdef LES
         if( i_qi == 2 ) then
            write(evb_unit,'(A)')
            write(evb_unit,1000) '{QI rate: f_v, F, G}: ', f_v, F_QI, G_QI
         endif
         if( i_qi > 0 ) then
            write(evb_unit,'(A)')
            write(evb_unit,1000) '{QI rate: gradRC_norm}: ', ( gradRC_norm(n), n = 1, ndiv )
         endif
#endif

      endif 

   endif 

  888 format( A, 2X, I10, 2X, A, 2X, F20.8 )
#ifdef LES
  999 format( A/, 3(2X, F20.8) )
#endif
 1000 format( A/, (5(2X,F20.8)) )

!  +---------------------------------------------------------------+
!  +---------------------------------------------------------------+

   end subroutine out_evb


