! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Read inputs for EVB calculation                                        |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#include "../include/assert.fh"
#include "../include/dprec.fh"

   subroutine evb_outoptions

   use evb_parm,  only: evb_dcrypt, nevb, nxch, nbias, nmorse, lambda &
                      , xch_type, evb_dyn, ntw_evb, C_evb, C_xch, bias_ndx &
                      , k_umb, r0_umb, dbonds_RC, bond_RC, xch_expdat &
                      , xch_gaussdat
   use evb_amber, only: morse
   use miller,    only: ti_mass, nti_mass, do_ti_mass
   use file_io_dat


   implicit none

#  include "../include/memory.h"

   !  ..........................................................................

   integer :: n

!  +---------------------------------------------------------------------------+
!  |  Output EVB options                                                       |
!  +---------------------------------------------------------------------------+

   write(6,'(/a)') 'EVB options:'
   write(6,100) 'nevb = ', nevb, ', nbias  = ', nbias &
             ,', nmorse = ', nmorse, ', ntw_evb = ', ntw_evb

   write(6,'(5x,a)' ) 'xch_type = ' // trim( adjustl( xch_type ) )
   write(6,'(5x,a)' ) 'evb_dyn  = ' // trim( adjustl( evb_dyn  ) )

   write(6,200) ( 'dia_shift(', n, ') = ', C_evb(n), n = 1, nevb )

   if( nmorse > 0 ) &
      write(6,300) ('morsify(',morse(n)%iatom,',',morse(n)%jatom,') ::', &
          'd = ',morse(n)%d,', a = ',morse(n)%a,', r0 = ',morse(n)%r0, &
          n=1,nmorse)

   if( nbias > 0 ) then
      select case( trim( adjustl( evb_dyn) ) )
         case( "evb_map" )
            write(6,400) ('emap(',bias_ndx(n,1),',',bias_ndx(n,2),') ::', &
                'lambda = ',lambda,n=1,nbias)
         case( "egap_umb" )
            write(6,500) ('egap_umb(',bias_ndx(n,1),',',bias_ndx(n,2),') ::', &
                'k = ',k_umb(n),'ezero= ',r0_umb(n),n=1,nbias)
         case( "dbonds_umb", "qi_dbonds_pmf", "qi_dbonds_dyn" )
            write(6,600) ('dbonds_umb(',dbonds_RC(n)%iatom,',', &
                dbonds_RC(n)%jatom,',',dbonds_RC(n)%katom,') ::','k = ', &
                k_umb(n),'ezero= ',r0_umb(n),n=1,nbias)
         case( "bond_umb", "qi_bond_pmf", "qi_bond_dyn" )
            write(6,700) ('bond_umb(',bond_RC(n)%iatom,',',bond_RC(n)%jatom, &
                ') ::','k = ',k_umb(n),'ezero= ',r0_umb(n),n=1,nbias)
      end select
   endif

   select case( trim( adjustl( xch_type ) ) )
      case( "constant" )
         write(6,800) ('xch_cnst(',evb_dcrypt(n,1),',',evb_dcrypt(n,2), &
             ') = ',C_xch(n),n=1,nxch)
      case( "exp" )
         write(6,900) ('xch_exp(',evb_dcrypt(n,1),',',evb_dcrypt(n,2), &
             ') ::','iatom = ',xch_expdat(n)%iatom,', jatom = ', &
             xch_expdat(n)%jatom,', a = ',xch_expdat(n)%A,', u = ', &
             xch_expdat(n)%u,', r0 = ', xch_expdat(n)%r0,n = 1,nxch)
      case( "gauss" )
         write(6,900) ('xch_gauss(',evb_dcrypt(n,1),',',evb_dcrypt(n,2), &
             ') ::','iatom = ',xch_gaussdat(n)%iatom,', jatom = ', &
             xch_gaussdat(n)%jatom,', a = ',xch_gaussdat(n)%A,', sigma = ', &
             xch_gaussdat(n)%sigma,', r0 = ',xch_gaussdat(n)%r0,n=1,nxch)
   end select

   if( do_ti_mass) then
      write(6,1000) ('ti_mass(',ti_mass(n)%atm_ndx,') ::','lambda = ', &
          ti_mass(n)%lambda,', M_i = ',ti_mass(n)%init,', M_f = ', &
          ti_mass(n)%final,n=1,nti_mass)
   endif

   write(6,'(/)')

 100  format( ( 5x, 4 ( a, 2x, i4 ) ) )
 200  format( ( 5x, a, i2, a, f10.5 ) )
 300  format( ( 5x, 2 ( a, i10 ), a, 2x, 3 ( a, f10.5 ) ) )
 400  format( ( 5x, 2 ( a, i4 ), a, 2x, a, f10.5 ) )
 500  format( ( 5x, 2 ( a, i4 ), a, 2x, 2 ( 2x, a, f10.5 ) ) )
 600  format( ( 5x, 3 ( a, i10 ), a, 2x, 2 ( 2x, a, f10.5, 2x ) ) )
 700  format( ( 5x, 2 ( a, i10 ), a, 2x, 2 ( 2x, a, f10.5, 2x ) ) )
 800  format( ( 5x, 2 ( a, i4 ), a, f10.5 ) )
 900  format( ( 5x, 2 ( a, i4 ), a, 2x, 2 ( a, i10 ), 3 ( a, f10.5 ) ) )
1000  format( ( 5x, ( a, i10 ), a, 2x, 3 ( a, f10.5 ) ) )


!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine evb_outoptions


