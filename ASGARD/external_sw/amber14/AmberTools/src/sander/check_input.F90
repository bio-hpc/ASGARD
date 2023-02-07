! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Check EVB variable dependencies                                        |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   subroutine check_input ( xch_check, dyn_check, dia_check, xch_type &
                          , evb_dyn ) 

   use file_io_dat
   implicit none

#  include "parallel.h"

   integer, intent(in) :: xch_check(4), dyn_check(4), dia_check(1) 
   character(512), intent(in) :: xch_type, evb_dyn

!  +---------------------------------------------------------------------------+
!  |  Check exchange type inputs                                               |
!  +---------------------------------------------------------------------------+

   select case( trim( adjustl( xch_type ) ) )

      case( "constant" ) 
         if( xch_check(1) /= 1 ) then 
            write(6,'(A)')
            write(6,'(A)') 'xch_type = ' // trim( adjustl( xch_type ) ) &
                        // ' but incomplete data specified in ' &
                        // trim( adjustl( evbin ) ) 
            write(6,'(A)') '==> need [constant exchange] block'
            call mexit(6,1)
         endif 

      case( "exp" )
         if( xch_check(2) /= 1 ) then
            write(6,'(A)')
            write(6,'(A)') 'xch_type = ' // trim( adjustl( xch_type ) ) &
                        // ' but incomplete data specified in ' &
                        // trim( adjustl( evbin ) )
            write(6,'(A)') '==> need [exponential exchange] block'
            call mexit(6,1)
         endif

      case( "gauss" )
         if( xch_check(3) /= 1 ) then
            write(6,'(A)')
            write(6,'(A)') 'xch_type = ' // trim( adjustl( xch_type ) ) &
                        // ' but incomplete data specified in ' &
                        // trim( adjustl( evbin ) )
            write(6,'(A)') '==> need [gaussian exchange] block'
            call mexit(6,1)
         endif

      case( "chang_miller", "minichino_voth" )
         if( dia_check(1) /= 1 ) then
            write(6,'(A)')
            write(6,'(A)') 'xch_type = ' // trim( adjustl( xch_type ) ) &
                        // ' but incomplete data specified in ' &
                        // trim( adjustl( evbin ) )
            write(6,'(A)') '==> need [diabatic state external files] block'
            call mexit(6,1)
         endif 
         if( xch_check(4) /= 1 ) then
            write(6,'(A)')
            write(6,'(A)') 'xch_type = ' // trim( adjustl( xch_type ) ) &
                        // ' but incomplete data specified in ' &
                        // trim( adjustl( evbin ) )
            write(6,'(A)') '==> need [exchange element external files] block'
            call mexit(6,1)
         endif

   end select

!  +---------------------------------------------------------------------------+
!  |  Check dynamics type inputs                                               |
!  +---------------------------------------------------------------------------+

   select case( trim( adjustl( evb_dyn ) ) )

      case( "evb_map" )
         if( dyn_check(1) /= 1 ) then
            write(6,'(A)')
            write(6,'(A)') 'evb_dyn = ' // trim( adjustl( evb_dyn ) ) &
                        // ' but incomplete data specified in ' &
                        // trim( adjustl( evbin ) )
            write(6,'(A)') '==> need [mapping potential] block'
            call mexit(6,1)
         endif 

      case( "egap_umb" )
         if( dyn_check(2) /= 1 ) then
            write(6,'(A)')
            write(6,'(A)') 'evb_dyn = ' // trim( adjustl( evb_dyn ) ) &
                        // ' but incomplete data specified in ' &
                        // trim( adjustl( evbin ) )
            write(6,'(A)') '==> need [energy gap RC umbrella] block'
            call mexit(6,1)
         endif 

      case( "dbonds_umb", "qi_dbonds_pmf", "qi_dbonds_dyn" )
         if( dyn_check(3) /= 1 ) then
            write(6,'(A)')
            write(6,'(A)') 'evb_dyn = ' // trim( adjustl( evb_dyn ) ) &
                        // ' but incomplete data specified in ' &
                        // trim( adjustl( evbin ) )
            write(6,'(A)') '==> need [difference of 2 bonds RC umbrella] block'
            call mexit(6,1)
         endif

      case( "bond_umb", "qi_bond_pmf", "qi_bond_dyn" )
         if( dyn_check(4) /= 1 ) then
            write(6,'(A)')
            write(6,'(A)') 'evb_dyn = ' // trim( adjustl( evb_dyn ) ) &
                        // ' but incomplete data specified in ' &
                        // trim( adjustl( evbin ) )
            write(6,'(A)') '==> need [difference of 2 bonds RC umbrella] block'
            call mexit(6,1)
         endif

   end select

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine check_input

