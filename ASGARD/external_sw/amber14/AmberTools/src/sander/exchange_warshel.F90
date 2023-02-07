! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Compute Warshel's exponential functional form for the coupling         |#
! #|  J. Am. Chem. Soc., 110, 5297-5311 (1988)                               |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#include "../include/assert.fh"
#include "../include/dprec.fh"

   subroutine exchange_warshel ( q, ndim ) 

   use evb_parm,  only: nxch, xch_expdat 
   use evb_xchff, only: xch_warshel
#ifdef DEBUG_EVB
   use evb_check, only: xwarshel_debug
#endif

   implicit none

   integer, intent(in) :: ndim
   _REAL_ , intent(in) :: q(ndim)

   !  ..........................................................................

   integer :: n, idx, jdx 
   _REAL_  :: dx, dy, dz, rij, ff, A, u, r0 

!  +---------------------------------------------------------------------------+
!  |  Each coupling term is of the form                                        |
!  |                                                                           |
!  |  V_kl = A_kl * exp [ - u_kl ( r_ij - r_kl^0 ) ]                           |
!  +---------------------------------------------------------------------------+

   do n = 1, nxch

      A  = xch_expdat(n)%A
      u  = xch_expdat(n)%u
      r0 = xch_expdat(n)%r0

      idx = ( xch_expdat(n)%iatom - 1 ) * 3 
      jdx = ( xch_expdat(n)%jatom - 1 ) * 3

      dx = q(idx+1) - q(jdx+1) 
      dy = q(idx+2) - q(jdx+2) 
      dz = q(idx+3) - q(jdx+3) 

      rij = sqrt( dx**2 + dy**2 + dz**2 ) 
      xch_warshel(n)%xch = A * exp( - u * ( rij - r0 ) ) 

!  +---------------------------------------------------------------------------+
!  |  Compute gradient of the coupling matrix element                          |
!  |                                                                           |
!  |  grad_Ri V_kl = V_kl * [ -u_kl / rij * ( R_i - R_j ) ]                    |
!  |  grad_Rj V_kl = - grad_Ri V_kl                                            |
!  +---------------------------------------------------------------------------+

      xch_warshel(n)%gxch(:) = 0.0d0

      ff = - u * xch_warshel(n)%xch / rij

      xch_warshel(n)%gxch(idx+1) = ff * dx
      xch_warshel(n)%gxch(idx+2) = ff * dy
      xch_warshel(n)%gxch(idx+3) = ff * dz

      xch_warshel(n)%gxch(jdx+1) = - xch_warshel(n)%gxch(idx+1)
      xch_warshel(n)%gxch(jdx+2) = - xch_warshel(n)%gxch(idx+2)
      xch_warshel(n)%gxch(jdx+3) = - xch_warshel(n)%gxch(idx+3)

   enddo

#ifdef DEBUG_EVB
   if( xwarshel_debug ) call xwarshel_anal2num ( q, ndim )
#endif

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine exchange_warshel


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Allocate storage for xch_warshel                                       |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   subroutine xch_warshel_alloc ( ndim )

   use evb_parm,  only: nxch 
   use evb_xchff, only: xch_warshel

   implicit none 

   integer, intent(in) :: ndim

   !  ..........................................................................

   integer :: n, alloc_error

   allocate( xch_warshel(nxch), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

   do n = 1, nxch
      allocate( xch_warshel(n)%gxch(ndim), stat = alloc_error )
      REQUIRE( alloc_error == 0 )
   enddo

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine xch_warshel_alloc


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Deallocate storage for xch_warshel                                     |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   subroutine xch_warshel_dealloc

   use evb_parm,  only: nxch
   use evb_xchff, only: xch_warshel

   implicit none

   !  ..........................................................................

   integer :: n, dealloc_error

   do n = 1, nxch 
      deallocate( xch_warshel(n)%gxch, stat = dealloc_error )
      REQUIRE( dealloc_error == 0 )
   enddo 

   deallocate( xch_warshel, stat = dealloc_error )
   REQUIRE( dealloc_error == 0 )

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine xch_warshel_dealloc


