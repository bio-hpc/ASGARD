! <compile=optimized>

! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  REACT_FLUX_INIT: Initialize reactive flux                      |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#include "../include/assert.fh"
#include "../include/dprec.fh"

   subroutine react_flux_init ( xf, xq, ndim )

   use evb_parm,  only: bias_ndx, nevb, dbonds_RC, bond_RC 
   use evb_data,  only: evb_Hmat 
   use wigner,    only: rflux, RC_type, nbias_ndx, s0_surface, s0_toler
#ifdef DEBUG_EVB
   use evb_check, only: full_evb_debug, gradRC_debug
#endif

   implicit none

   integer, intent(in) :: ndim
   _REAL_ , intent(in) :: xf(ndim,nevb), xq(ndim)

   !  +---------------------------------------------------------------+

   integer :: n, nn, ni, nf, idx, jdx, kdx
   _REAL_  :: fharm(ndim), dr(3), fr(3), rij, rkj, rij_inv, rkj_inv


   call react_flux_alloc ( ndim )

   select case( trim( adjustl( RC_type ) ) )

!  +---------------------------------------------------------------+
!  |  Sampling according to Warshel's mapping potential &          |
!  |  Case's umbrella sampling along an energy gap /\              |
!  |                                                               |
!  |  /\ = V_ii - V_ff                                             |
!  |                                                               |
!  | - d/dR /\  = F_ii - F_ff                                      |
!  +---------------------------------------------------------------+

      case( "egap" )

         n = nbias_ndx

         ni = bias_ndx(n,1)
         nf = bias_ndx(n,2)

         rflux%normal_vect(:) = - ( xf(:,ni) - xf(:,nf) )

         rflux%RC_s0 = evb_Hmat%evb_mat(ni,ni) - evb_Hmat%evb_mat(nf,nf)

!  +---------------------------------------------------------------+
!  |  Difference of 2 bonds RC harmonic umbrella sampling          |
!  |                                                               |
!  |  /\ = r_ij - r_kj                                             |
!  |                                                               |
!  |  - d/dR_i /\ = - ( R_i - R_j ) / rij                          |
!  |  - d/dR_k /\ = + ( R_k - R_j ) / rkj                          |
!  |  - d/dR_j /\ = + ( R_i - R_j ) / rij - ( R_k - R_j ) / rkj    |
!  +---------------------------------------------------------------+

      case( "dbonds" )

         fharm(:) = 0.0d0

         n = nbias_ndx

         idx = ( dbonds_RC(n)%iatom - 1 ) * 3
         jdx = ( dbonds_RC(n)%jatom - 1 ) * 3
         kdx = ( dbonds_RC(n)%katom - 1 ) * 3

         do nn = 1, 3
            dr(nn) = xq(idx+nn) - xq(jdx+nn)
         enddo 

         rij = sqrt( dr(1)**2 + dr(2)**2 + dr(3)**2 )
         rij_inv = 1.0d0 / rij

         do nn = 1, 3
            fr(nn) = dr(nn) * rij_inv
            fharm(idx+nn) = fharm(idx+nn) - fr(nn)
            fharm(jdx+nn) = fharm(jdx+nn) + fr(nn)
         enddo 

         do nn = 1, 3 
            dr(nn)  = xq(kdx+nn) - xq(jdx+nn)
         enddo

         rkj = sqrt( dr(1)**2 + dr(2)**2 + dr(3)**2 )
         rkj_inv = 1.0d0 / rkj

         do nn = 1, 3
            fr(nn) = dr(nn) * rkj_inv
            fharm(kdx+nn) = fharm(kdx+nn) + fr(nn)
            fharm(jdx+nn) = fharm(jdx+nn) - fr(nn)
         enddo

#ifdef DEBUG_EVB
         if( full_evb_debug .or. gradRC_debug ) &
            call RCdbonds_anal2num ( xq, fharm, n, ndim )
#endif

         rflux%normal_vect(:) =  - fharm(:)

         rflux%RC_s0 = rij - rkj

!  +---------------------------------------------------------------+
!  |  Bond RC harmonic umbrella sampling                           |
!  |                                                               |
!  |  /\ = r_ij                                                    |
!  |                                                               |
!  |  - d/dR_i /\ = - ( R_i - R_j ) / rij                          |
!  |  - d/dR_j /\ = + ( R_i - R_j ) / rij                          |
!  +---------------------------------------------------------------+

      case( "bond" )

         fharm(:) = 0.0d0

         n = nbias_ndx

         idx = ( bond_RC(n)%iatom - 1 ) * 3
         jdx = ( bond_RC(n)%jatom - 1 ) * 3

         do nn = 1, 3
            dr(nn) = xq(idx+nn) - xq(jdx+nn)
         enddo

         rij = sqrt( dr(1)**2 + dr(2)**2 + dr(3)**2 )
         rij_inv = 1.0d0 / rij

         do nn = 1, 3
            fr(nn) = dr(nn) * rij_inv
            fharm(idx+nn) = fharm(idx+nn) - fr(nn)
            fharm(jdx+nn) = fharm(jdx+nn) + fr(nn)
         enddo

#ifdef DEBUG_EVB
         if( gradRC_debug ) call RCbond_anal2num ( xq, fharm, n, ndim )
#endif

         rflux%normal_vect(:) = - fharm(:)

         rflux%RC_s0 = rij 

   end select

   if( abs( rflux%RC_s0 - s0_surface ) > s0_toler ) then 
      write(6,'(A)') '|WARNING: RC value for config is outside of ' &
                  // 'tolerance for placement of s0_surface '

      write(6,'(A, F14.4)') 'Initial RC = ', rflux%RC_s0
      write(6,'(A, F14.4)') 's0_surface = ', s0_surface
      write(6,'(A, F14.4)') 's0_toler   = ', s0_toler

   endif 

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine react_flux_init


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Allocate storage space for react_flux                          |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   subroutine react_flux_alloc (ndim )

   use wigner, only: rflux, nreact_RC, p_space

   implicit none

   integer, intent(in) :: ndim

   !  +---------------------------------------------------------------+

   integer :: alloc_error

   allocate( rflux%normal_vect(ndim), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

   allocate( rflux%normal_velv(ndim), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

   allocate( rflux%RC_sampled(nreact_RC), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

   allocate( rflux%RC_trj(nreact_RC), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

   allocate( rflux%coord_TS(ndim), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

   allocate( rflux%force_TS(ndim), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

   allocate( rflux%xi(p_space), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

   end subroutine react_flux_alloc


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Deallocate storage space for react_flux                        |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   subroutine react_flux_dealloc

   use wigner, only: rflux

   implicit none

   integer :: dealloc_error

   deallocate( rflux%normal_vect, rflux%normal_velv &
             , rflux%RC_sampled, rflux%coord_TS, rflux%force_TS &
             , rflux%xi, rflux%RC_trj, stat = dealloc_error )
   REQUIRE( dealloc_error == 0 )

   end subroutine react_flux_dealloc
