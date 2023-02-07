! <compile=optimized>
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Populate the EVB matrix                                                |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#include "../include/assert.fh"
#include "../include/dprec.fh"

   subroutine evb_matrix ( xnrg, q, ndim )

   use evb_parm,  only: xch_type, evb_dcrypt, C_evb, C_xch, nevb, nxch
   use evb_xchff, only: xch_warshel, xch_gauss
   use evb_data,  only: evb_Hmat

#if defined(LES)
   use evb_pimd,  only: nbead
#endif

   implicit none

   integer, intent(in ) :: ndim
   _REAL_ , intent(in ) :: xnrg(nevb), q(ndim)

   !  ..........................................................................

   integer :: n, k, l

!  +---------------------------------------------------------------------------+
!  |  Populate diagonal terms                                                  |
!  +---------------------------------------------------------------------------+

   evb_Hmat%evb_mat(:,:) = 0.0d0

   do n = 1, nevb 
      evb_Hmat%evb_mat(n,n) = xnrg(n) + C_evb(n)
   enddo

!  +---------------------------------------------------------------------------+
!  |  Compute coupling terms                                                   |
!  +---------------------------------------------------------------------------+

   select case( trim( adjustl( xch_type ) ) )

      case( "constant" )

         do n = 1, nxch

            k = evb_dcrypt(n,1)
            l = evb_dcrypt(n,2)

            evb_Hmat%evb_mat(k,l) = C_xch(n)
            evb_Hmat%evb_mat(l,k) = evb_Hmat%evb_mat(k,l)

            evb_Hmat%grad_xch(:) = 0.0d0 
            evb_Hmat%b1_kl(:) = 0.0d0

         enddo

      case( "exp" )

         call exchange_warshel ( q, ndim )

         do n = 1, nxch

            k = evb_dcrypt(n,1)
            l = evb_dcrypt(n,2)

            evb_Hmat%evb_mat(k,l) = xch_warshel(n)%xch
            evb_Hmat%evb_mat(l,k) = evb_Hmat%evb_mat(k,l)

            evb_Hmat%grad_xch(:) = xch_warshel(n)%gxch(:)

            evb_Hmat%b1_kl(:) = 2.0d0 * xch_warshel(n)%xch &
                              * evb_Hmat%grad_xch(:)

         enddo

      case( "gauss" )

         call exchange_gauss ( q, ndim )

         do n = 1, nxch

            k = evb_dcrypt(n,1)
            l = evb_dcrypt(n,2)

            evb_Hmat%evb_mat(k,l) = xch_gauss(n)%xch
            evb_Hmat%evb_mat(l,k) = evb_Hmat%evb_mat(k,l)

            evb_Hmat%grad_xch(:) = xch_gauss(n)%gxch(:)

            evb_Hmat%b1_kl(:) = 2.0d0 * xch_gauss(n)%xch &
                              * evb_Hmat%grad_xch(:)

         enddo

      case( "polynom_gaussian" )

!        call exchange_PG ( q )

!  +---------------------------------------------------------------------------+
!  |  Important note: data from Gaussian are in atomic units ... so the DG     |
!  |  must be done in a.u. and then the energies and forces are converted      |
!  |  back to Amber kcal/mol and Angstroms                                     |
!  +---------------------------------------------------------------------------+

      case( "dist_gauss" )

   end select 

!  +---------------------------------------------------------------------------+
!  |  Scale exchange terms for PIMD                                            |
!  +---------------------------------------------------------------------------+

#if defined(LES)

   do n = 1, nxch 

      k = evb_dcrypt(n,1)
      l = evb_dcrypt(n,2)

      evb_Hmat%b1_kl(:) = evb_Hmat%b1_kl(:) / ( nbead**2 )

      evb_Hmat%evb_mat(k,l) = evb_Hmat%evb_mat(k,l) / nbead
      evb_Hmat%evb_mat(l,k) = evb_Hmat%evb_mat(k,l)

      evb_Hmat%grad_xch(:) = 0.50d0 / evb_Hmat%evb_mat(k,l) &
                                     * evb_Hmat%b1_kl(:) 

   enddo

#endif /* LES */

!  +---------------------------------------------------------------------------+
!  +---------------------------------------------------------------------------+

   end subroutine evb_matrix


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Allocate storage for evb_matrix                                        |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   subroutine evb_mat_type_alloc ( ndim )

   use evb_parm, only: nevb 
   use evb_data, only: evb_Hmat 

   implicit none

   integer, intent(in) :: ndim

   !  ..........................................................................

   integer :: alloc_error

   allocate( evb_Hmat%evb_mat(nevb,nevb), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

   allocate( evb_Hmat%grad_xch(ndim), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

   allocate( evb_Hmat%b1_kl(ndim), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

   allocate( evb_Hmat%evb_vec0(nevb), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

   end subroutine evb_mat_type_alloc


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Deallocate storage for evb_matrix                                      |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   subroutine evb_mat_type_dealloc

   use evb_data, only: evb_Hmat

   implicit none

   !  ..........................................................................

   integer :: dealloc_error

   deallocate( evb_Hmat%evb_mat, evb_Hmat%grad_xch, evb_Hmat%b1_kl &
             , stat = dealloc_error )
   REQUIRE( dealloc_error == 0 )

   end subroutine evb_mat_type_dealloc


