! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Compute gaussian functional form for the coupling                      |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#include "../include/assert.fh"
#include "../include/dprec.fh"

   subroutine exchange_gauss ( q, ndim ) 

   use evb_parm,  only: nxch, xch_gaussdat 
   use evb_xchff, only: xch_gauss
#ifdef DEBUG_EVB
   use evb_check, only: xgauss_debug
#endif

   implicit none

   integer, intent(in) :: ndim
   _REAL_ , intent(in) :: q(ndim)

   !  ..........................................................................

   integer :: n, idx, jdx 
   _REAL_  :: dx, dy, dz, rij, ff, A, sigma, r0 

!  +---------------------------------------------------------------------------+
!  |  Each coupling term is of the form                                        |
!  |                                                                           |
!  |  V_kl = A_kl * exp [ - ( r_ij - r_kl^0 )^2 / sigma^2 ]                    |
!  +---------------------------------------------------------------------------+

   do n = 1, nxch

      A      = xch_gaussdat(n)%A
      sigma  = xch_gaussdat(n)%sigma
      r0     = xch_gaussdat(n)%r0

      idx = ( xch_gaussdat(n)%iatom - 1 ) * 3 
      jdx = ( xch_gaussdat(n)%jatom - 1 ) * 3

      dx = q(idx+1) - q(jdx+1) 
      dy = q(idx+2) - q(jdx+2) 
      dz = q(idx+3) - q(jdx+3) 

      rij = sqrt( dx**2 + dy**2 + dz**2 )
      xch_gauss(n)%xch = A * exp( - ( rij - r0 )**2 / sigma**2  ) 

!  +---------------------------------------------------------------------------+
!  |  Compute gradient of the coupling matrix element                          |
!  |                                                                           |
!  |  grad_Ri V_kl = V_kl * [ - 2 * ( r_ij - r_kl^0 ) / sigma**2 / rij         |
!  |               * ( R_i - R_j ) ]                                           |
!  |  grad_Rj V_kl = - grad_Ri V_kl                                            |
!  +---------------------------------------------------------------------------+

      xch_gauss(n)%gxch(:) = 0.0d0

      ff = - 2.0d0 * ( rij - r0 ) / sigma**2 / rij * xch_gauss(n)%xch

      xch_gauss(n)%gxch(idx+1) = ff * dx 
      xch_gauss(n)%gxch(idx+2) = ff * dy 
      xch_gauss(n)%gxch(idx+3) = ff * dz 

      xch_gauss(n)%gxch(jdx+1) = - xch_gauss(n)%gxch(idx+1)
      xch_gauss(n)%gxch(jdx+2) = - xch_gauss(n)%gxch(idx+2)
      xch_gauss(n)%gxch(jdx+3) = - xch_gauss(n)%gxch(idx+3)

   enddo

#ifdef DEBUG_EVB
   if( xgauss_debug ) call xgauss_anal2num ( q, ndim )
#endif

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine exchange_gauss


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Allocate storage for exchange_gauss                                    |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   subroutine xch_gauss_alloc ( ndim )

   use evb_parm,  only: nxch 
   use evb_xchff, only: xch_gauss

   implicit none 

   integer, intent(in) :: ndim

   !  ..........................................................................

   integer :: n, alloc_error

   allocate( xch_gauss(nxch), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

   do n = 1, nxch
      allocate( xch_gauss(n)%gxch(ndim), stat = alloc_error )
      REQUIRE( alloc_error == 0 )
   enddo

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine xch_gauss_alloc


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Deallocate storage for exchange_gauss                                  |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   subroutine xch_gauss_dealloc

   use evb_parm,  only: nxch
   use evb_xchff, only: xch_gauss

   implicit none

   !  ..........................................................................

   integer :: n, dealloc_error

   do n = 1, nxch 
      deallocate( xch_gauss(n)%gxch, stat = dealloc_error )
      REQUIRE( dealloc_error == 0 )
   enddo 

   deallocate( xch_gauss, stat = dealloc_error )
   REQUIRE( dealloc_error == 0 )

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine xch_gauss_dealloc

