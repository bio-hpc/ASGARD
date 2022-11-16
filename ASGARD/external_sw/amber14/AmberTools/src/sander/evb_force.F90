! <compile=optimized>
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  EVB_FORCE                                                              |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#include "../include/assert.fh"
#include "../include/dprec.fh"

#ifdef LES
   subroutine evb_force ( xf, xq, ndim )
#else
   subroutine evb_force ( xf, xq, mass, ndim )
#endif

   use evb_parm,  only: evb_dcrypt, xch_type, nevb, nxch
   use evb_data,  only: evb_Hmat, evb_frc, evb_vel0
   use constants, only: A_TO_BOHRS, AU_TO_KCAL, kB
   use schlegel,  only: ncoord, vmm
   use file_io_dat
#ifdef DEBUG_EVB
   use evb_check, only: full_evb_debug
#endif
#ifdef LES
   use pimd_vars, only: nbead_inv
#else
   use evb_data, only : RCdot, evb_bias
   use evb_parm, only : out_RCdot, r0_umb, nbias, lambda, k_umb, inc_bond_RC &
                      , evb_dyn, inc_dbonds_RC, egapRC, bias_ndx, bond_RC &
                      , dbonds_RC
#endif

   implicit none

#  include "parallel.h"
#  include "../include/md.h"

   integer, intent(in) :: ndim
   _REAL_ , intent(in) :: xf(ndim,nevb), xq(ndim)

   !  ..........................................................................

   integer :: k, l, n, info
   _REAL_  :: evb_vec(nevb,nevb), nrg_evb(nevb)
   _REAL_  :: work(26*nevb) 
   _REAL_  :: bmat(ncoord,ndim), dEVB(ncoord), ddEVB(ncoord,ncoord), vkl_sq &
            , EVB_NRG, qint(ncoord), vsub
#ifndef LES
   _REAL_  :: dr(3), RC, pi, rij, rkj, mass(ndim/3), gradRC(ndim)
   integer :: nn, ni, nf, jdx, kdx, idx
#endif

!  +---------------------------------------------------------------------------+
!  |  Analytical solution to a 2x2 EVB matrix                                  |
!  |                                                                           |
!  |  e_0 = 0.50 ( V11 + V22 ) - sqrt{ [ 0.50 ( V11 - V22 ) ]^2 + V12^2}       |
!  |  c_01 = ( e_0 - V22 ) / ||c_0||           c_02 = V12 / ||c_0||            |
!  +---------------------------------------------------------------------------+

   if( trim( adjustl( xch_type ) ) == "dist_gauss" ) then

!   write(111,*) 'AU_TO_KCAL = ', AU_TO_KCAL
!#ifdef LES
!   write(111,*) 'nbead_inv = ', nbead_inv
!#endif
!   write(111,*) 'ndim = ', ndim
!   write(111,*) xq(:)
!   write(111,*)

      call cart2internal ( xq * A_TO_BOHRS, qint )
      call schlegel_EVBR ( qint, vkl_sq, EVB_NRG, dEVB(:), ddEVB(:,:) )
      call wdc_bmat ( xq * A_TO_BOHRS, bmat )

!  write(222,*) qint(:)

      n = 1
      k = evb_dcrypt(n,1)
      l = evb_dcrypt(n,2)

!  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
!  |  Patch the case where the argument of the sqrt is negative                |
!  :```````````````````````````````````````````````````````````````````````````:
!  |  In schlegel_EVBR, we set sqrt{ [ 0.50 ( V11 - V22 ) ]^2 + V12^2} = 0     |
!  |    ==> V12^2 = [ 0.50 ( V11 - V22 ) ]^2                                   |
!  |    ... let's simply take the positive root.                               |
!  |  Note that V12^2 can be negative even if the argument of sqrt{} > 0.      |
!  |  That is why we use sign(sqrt(abs(*)),*) below ... i.e., the EVB          |
!  |  amplitudes can have an arbitrary phase which does not impact the         |
!  |  populations.  In general, care must be taken to ensure that this         |
!  |  arbitrary phase does not change between time steps but since our forces  |
!  |  are computed analytically (without explicitly resorting to V12 as in     |
!  |  the Hellman-Feynman force from matrix elements), we are OK here.         |
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

      vsub = 0.50d0 * ( vmm(k) - vmm(l) )
      if( vsub**2 + vkl_sq < 0.0d0 ) vkl_sq = vsub**2

      evb_vec(1,1) = EVB_NRG - vmm(2)
      evb_vec(2,1) = sign( sqrt( abs(vkl_sq ) ), vkl_sq )
      evb_vec(:,1) = evb_vec(:,1) / sqrt( evb_vec(1,1)**2 + evb_vec(2,1)**2 )
      evb_Hmat%evb_vec0(:) = evb_vec(:,1)

      evb_Hmat%evb_mat(k,l) = vkl_sq * AU_TO_KCAL * AU_TO_KCAL
      evb_Hmat%evb_mat(l,k) = evb_Hmat%evb_mat(k,l)
      evb_Hmat%evb_mat(k,k) = vmm(k) * AU_TO_KCAL
      evb_Hmat%evb_mat(l,l) = vmm(l) * AU_TO_KCAL

      nrg_evb(1) = EVB_NRG * AU_TO_KCAL
      evb_vel0%evb_f(:) = - matmul( dEVB(:), bmat(:,:) ) * AU_TO_KCAL * A_TO_BOHRS
#ifdef LES
      evb_vel0%evb_f(:) = evb_vel0%evb_f(:) * nbead_inv
      evb_Hmat%evb_mat(:,:) = evb_Hmat%evb_mat(:,:) * nbead_inv
      nrg_evb(1) = nrg_evb(1) * nbead_inv
#endif

!    write(666,*) evb_Hmat%evb_mat(:,:)

   else 

!  +---------------------------------------------------------------------------+
!  |  Analytical solution of a 2X2 EVB matrix                                  |
!  |                                                                           |
!  |  e_0 = 0.50 ( V11 + V22 ) - sqrt{ [ 0.50 ( V11 - V22 ) ]^2 +V12^2}        |
!  |  c_01 = ( e_0 - V22 ) / ||c_0||           c_02 = V12 / ||c_0||            |
!  +---------------------------------------------------------------------------+

      if( nevb == 2 ) then

         nrg_evb(1) = 0.50d0 * ( evb_Hmat%evb_mat(1,1) + evb_Hmat%evb_mat(2,2) ) &
            - sqrt( 0.25d0 * ( evb_Hmat%evb_mat(1,1) - evb_Hmat%evb_mat(2,2) )**2 &
            + evb_Hmat%evb_mat(1,2)**2 )

         evb_vec(1,1) = nrg_evb(1) - evb_Hmat%evb_mat(2,2)
         evb_vec(2,1) = evb_Hmat%evb_mat(1,2)
         evb_vec(:,1) = evb_vec(:,1) / sqrt( evb_vec(1,1)**2 + evb_vec(2,1)**2 )

      else
!  +---------------------------------------------------------------------------+
!  |  Diagonalize EVB matrix using LAPACK and form the Hellmann-Feynman force  |
!  |                                                                           |
!  |  F(:) = - sum(k,l) C_k1 * C_l1 * grad( V_kl )                             |
!  +---------------------------------------------------------------------------+

         evb_vec(:,:) = evb_Hmat%evb_mat(:,:)

         call D_OR_S()syev ( 'V', 'L', nevb, evb_vec, nevb, nrg_evb, work &
                    , 26*nevb, info )

         if( info /= 0 ) then 
            write(6,'(A)') 'Problem diagonalizing the EVB matrix using LAPACK DSYEV'
            write(6,'(A, I8)') 'info = ', info 
            call mexit(6,1)
         endif 

      endif 
 
      evb_Hmat%evb_vec0(:) = evb_vec(:,1)
      evb_vel0%evb_f(:) = 0.0d0 

      do n = 1, nevb 
         evb_vel0%evb_f(:) = evb_vel0%evb_f(:) + evb_vec(n,1)**2 * xf(:,n) 
      enddo 

      do n = 1, nxch 

         k = evb_dcrypt(n,1)
         l = evb_dcrypt(n,2) 

         evb_vel0%evb_f(:) = evb_vel0%evb_f(:) - 2.0d0 &
                           * evb_vec(k,1) * evb_vec(l,1) * evb_Hmat%grad_xch(:)

      enddo 

   endif 

   evb_vel0%evb_nrg = nrg_evb(1) 
   evb_frc%evb_f(:) = evb_vel0%evb_f(:)
   evb_frc%evb_nrg  = nrg_evb(1) 

#ifdef DEBUG_EVB
!  +---------------------------------------------------------------------------+
!  |  Debugging: analytical vs. numerical results for 2-state EVB              |
!  +---------------------------------------------------------------------------+

   if( nevb == 2 .and. full_evb_debug ) call evb_2stdebug ( evb_vec, xf, ndim )
#endif

!  +---------------------------------------------------------------------------+
!  |  The EVB/PIMD umbrella potential is done in evb_umb inside runmd.  Note   |
!  |  that PI-QTST is formulated in terms of RCs that are linear in the        |
!  |  centroid coordinate.  A PI-QTST theory based on the energy gap RC has    |
!  |  not been worked out yet.                                                 |
!  +---------------------------------------------------------------------------+

#ifndef LES

!  +---------------------------------------------------------------------------+
!  |  12062008: For outputting RC when doing groundstate dynamics              |
!  +---------------------------------------------------------------------------+

         do n = 1, nbias

            if( inc_dbonds_RC ) then

               idx = ( dbonds_RC(n)%iatom - 1 ) * 3
               jdx = ( dbonds_RC(n)%jatom - 1 ) * 3
               kdx = ( dbonds_RC(n)%katom - 1 ) * 3

               do nn = 1, 3
                  dr(nn) = xq(idx+nn) - xq(jdx+nn)
               enddo
               rij = sqrt( dr(1)**2 + dr(2)**2 + dr(3)**2 )

               do nn = 1, 3
                  dr(nn)  = xq(kdx+nn) - xq(jdx+nn)
               enddo
               rkj = sqrt( dr(1)**2 + dr(2)**2 + dr(3)**2 )

               RC = rij - rkj
               evb_bias%RC(n) = RC

            endif

            if( inc_bond_RC ) then

               idx = ( bond_RC(n)%iatom - 1 ) * 3
               jdx = ( bond_RC(n)%jatom - 1 ) * 3

               do nn = 1, 3
                  dr(nn) = xq(idx+nn) - xq(jdx+nn)
               enddo
               rij = sqrt( dr(1)**2 + dr(2)**2 + dr(3)**2 )

               RC = rij
               evb_bias%RC(n) = RC

            endif

         enddo


   select case( trim( adjustl( evb_dyn) ) )

!  +---------------------------------------------------------------------------+
!  |  Sampling according to Warshel's mapping potential                        |
!  |                                                                           |
!  |  /\ = V_i - V_f                                                           |
!  |                                                                           |
!  |  V_map(n) = V_i + lambda(n) * ( V_f - V_i )                               |
!  +---------------------------------------------------------------------------+

      case( "evb_map" )

         egapRC = .true.

         do n = 1, nbias

            ni = bias_ndx(n,1) 
            nf = bias_ndx(n,2) 

            evb_bias%RC(n) = evb_Hmat%evb_mat(ni,ni) - evb_Hmat%evb_mat(nf,nf)

            evb_bias%nrg_bias(n) = evb_Hmat%evb_mat(ni,ni) + lambda(n) &
               * ( evb_Hmat%evb_mat(nf,nf) - evb_Hmat%evb_mat(ni,ni) ) 

            evb_bias%evb_fbias(:,n) = ( 1.0d0 - lambda(n) ) * xf(:,ni) &
                           + lambda(n) * xf(:,nf)

            evb_frc%evb_f(:) = evb_bias%evb_fbias(:,n)
            evb_frc%evb_nrg  = evb_bias%nrg_bias(n)

            if( inc_dbonds_RC ) then

               idx = ( dbonds_RC(n)%iatom - 1 ) * 3
               jdx = ( dbonds_RC(n)%jatom - 1 ) * 3
               kdx = ( dbonds_RC(n)%katom - 1 ) * 3

               do nn = 1, 3
                  dr(nn) = xq(idx+nn) - xq(jdx+nn)
               enddo 
               rij = sqrt( dr(1)**2 + dr(2)**2 + dr(3)**2 )

               do nn = 1, 3 
                  dr(nn)  = xq(kdx+nn) - xq(jdx+nn) 
               enddo
               rkj = sqrt( dr(1)**2 + dr(2)**2 + dr(3)**2 )
         
               RC = rij - rkj
               evb_bias%RC(n) = RC

            endif

            if( inc_bond_RC ) then

               idx = ( bond_RC(n)%iatom - 1 ) * 3
               jdx = ( bond_RC(n)%jatom - 1 ) * 3

               do nn = 1, 3
                  dr(nn) = xq(idx+nn) - xq(jdx+nn)
               enddo
               rij = sqrt( dr(1)**2 + dr(2)**2 + dr(3)**2 )

               RC = rij 
               evb_bias%RC(n) = RC

            endif
 
         enddo 

!  +---------------------------------------------------------------------------+
!  |  Sampling according to Case's umbrella biasing potential                  |
!  |                                                                           |
!  |  /\ = V_i - V_f                                                           |
!  |                                                                           |
!  |  V_umb(eta) = V_0 + 0.5 k_evb * ( /\ - eta )^2                            |
!  |                                                                           |
!  |  dV_umb / dR = ( dV_umb / d/\ ) * ( d/\ / dR )                            |
!  |    = dV_0 / dR + k_evb * ( /\ - eta ) * d/dR ( V_i - V_f )                |
!  +---------------------------------------------------------------------------+

      case( "egap_umb" )

         egapRC = .true.

         do n = 1, nbias

            ni = bias_ndx(n,1)
            nf = bias_ndx(n,2)

            RC = evb_Hmat%evb_mat(ni,ni) - evb_Hmat%evb_mat(nf,nf)

            evb_bias%RC(n) = RC
            evb_bias%nrg_bias(n) = 0.50d0 * k_umb(n) * ( RC - r0_umb(n) )**2
            evb_bias%evb_fbias(:,n) = k_umb(n) * ( RC - r0_umb(n) ) &
                                    * ( xf(:,ni) - xf(:,nf) )
         enddo

         do n = 1, nbias
            evb_frc%evb_f(:) = evb_frc%evb_f(:) + evb_bias%evb_fbias(:,n)
            evb_frc%evb_nrg  = evb_frc%evb_nrg + evb_bias%nrg_bias(n)
         enddo

#ifdef DEBUG_EVB
!  +---------------------------------------------------------------------------+
!  |  Debugging: analytical vs. numerical results for 2-state EVB umbrella     |
!  |  sampling along an energy gap RC                                          |
!  +---------------------------------------------------------------------------+

         if( nevb == 2 .and. full_evb_debug ) &
            call egap_umb_2stdebug ( xf, nrg_evb(1), ndim )
#endif

   end select

   if( out_RCdot .and. egapRC ) then

      ni = bias_ndx(1,1)
      nf = bias_ndx(1,2)
      gradRC(:) = xf(:,ni) - xf(:,nf)

      pi = acos( -1.0d0 )
      RCdot = 0.0d0
      do n = 1, ndim/3
         idx = ( n - 1 ) * 3
         RCdot = RCdot + ( gradRC(idx+1)**2 + gradRC(idx+2)**2 &
                         + gradRC(idx+3)**2 ) / mass(n)
      enddo

      RCdot = sqrt( 2.0d0 * kB * temp0 * RCdot / pi )

   endif


#endif /* ! LES */

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

end subroutine evb_force


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Allocate storage for evb_frc_type                                      |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   subroutine evb_frc_type_alloc ( ndim )

   use evb_parm, only: nbias, evb_dyn 
   use evb_data, only: evb_frc, evb_vel0, evb_bias

   implicit none

   integer, intent(in) :: ndim

   !  +------------------------------------------------------------------------+

   integer :: alloc_error

   allocate( evb_frc%evb_f(ndim), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

   allocate( evb_vel0%evb_f(ndim), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

   allocate( evb_bias%RC(nbias), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

   select case( trim( adjustl( evb_dyn ) ) )

      case( "evb_map", "egap_umb", "dbonds_umb", "bond_umb",  &
            "qi_dbonds_pmf", "qi_bond_pmf", "qi_dbonds_dyn", "qi_bond_dyn" )

         allocate( evb_bias%nrg_bias(nbias), stat = alloc_error )
         REQUIRE( alloc_error == 0 )

         allocate( evb_bias%evb_fbias(ndim,nbias), stat = alloc_error )
         REQUIRE( alloc_error == 0 )

   end select

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine evb_frc_type_alloc


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Deallocate storage for evb_frc_type                                    |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   subroutine evb_frc_type_dealloc

   use evb_parm, only: evb_dyn
   use evb_data, only: evb_frc, evb_vel0, evb_bias

   implicit none

   integer :: dealloc_error

   deallocate( evb_frc%evb_f, evb_vel0%evb_f, stat = dealloc_error )
   REQUIRE( dealloc_error == 0 )

   select case( trim( adjustl( evb_dyn ) ) )

      case( "evb_map", "egap_umb", "dbonds_umb", "bond_umb",  &
            "qi_dbonds_pmf", "qi_bond_pmf", "qi_dbonds_dyn", "qi_bond_dyn" )
 
         deallocate( evb_bias%RC, evb_bias%nrg_bias, evb_bias%evb_fbias, stat = dealloc_error )
         REQUIRE( dealloc_error == 0 )

   end select

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine evb_frc_type_dealloc
