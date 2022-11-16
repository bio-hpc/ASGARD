! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  EVB_ALLOC                                                              |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#include "../include/dprec.fh"

   subroutine evb_alloc

   use evb_parm,  only: nbias, nevb, xch_type, bias_ndx, lambda, evb_dyn &
                      , ntw_evb, dbonds_RC, bond_RC, k_umb, r0_umb
   use evb_amber, only: evb_amber_ready 
   use evb_xchff, only: xchff_warshel_ready, xchff_gauss_ready
   use evb_data,  only: evb_Hmat, evb_frc_ready
   use schlegel,  only: ncoord, dgdim, dg_ready, dist_gauss
   use file_io_dat
!  use miller,    only: ti_mass, nti_mass, do_ti_mass
#if defined(LES)
   use pimd_vars, only: itimass
   use evb_pimd,  only: natomCL
#endif

   implicit none

#  include "parallel.h"
#  include "../include/memory.h"
#  include "../include/md.h"

   !  ..........................................................................

   integer :: n
!KFW   integer :: ti_mass_ndx(nti_mass*ncopy)


!  +---------------------------------------------------------------------------+
!  |  Exchange type                                                            |
!  +---------------------------------------------------------------------------+

   select case( trim( adjustl( xch_type ) ) )

!  :,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,:
!  :  Warshel's exponential functional form for the exchange term              :
!  :''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''':

      case( "exp" ) 
         if( .not. xchff_warshel_ready ) then
#if defined(LES)
            call xch_warshel_alloc ( natomCL * 3 )
#else
            call xch_warshel_alloc ( natom * 3 )
#endif
            xchff_warshel_ready = .true.
         endif

!  :,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,:
!  :  Gaussian functional form for the exchange term                           :
!  :''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''':

      case( "gauss" )
         if( .not. xchff_gauss_ready ) then
#if defined(LES)
            call xch_gauss_alloc ( natomCL * 3 )
#else
            call xch_gauss_alloc ( natom * 3 )
#endif
            xchff_gauss_ready = .true.
         endif

!  :,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,:
!  :  Berny's distributed gaussian PES fitting                                 :
!  :''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''':

      case( "dist_gauss" )
         if( trim( adjustl( dist_gauss%lin_solve ) ) == "diis" ) then
            dgdim = 1 + ncoord + ncoord * ( ncoord + 1 ) / 2
         else if( trim( adjustl( dist_gauss%lin_solve ) ) == "full" ) then
            dgdim = 1 + ncoord + ncoord * ncoord
         endif
         if( .not. dg_ready ) then
            call schlegel_dg_alloc
            dg_ready = .true.
         endif

   end select

!  +---------------------------------------------------------------------------+
!  |  Allocate memory for evb_Hmat in evb_matrix                               |
!  +---------------------------------------------------------------------------+

   if( .not. evb_Hmat%ready ) then
#if defined(LES)
      call evb_mat_type_alloc ( natomCL * 3 )
#else
      call evb_mat_type_alloc ( natom * 3 )
#endif
      evb_Hmat%ready = .true.
   endif 

!  +---------------------------------------------------------------------------+
!  |  Allocate memory for evb_frc in evb_force                                 |
!  +---------------------------------------------------------------------------+

   if( .not. evb_frc_ready ) then
#if defined(LES)
      call evb_frc_type_alloc ( natomCL * 3 )
#else
      call evb_frc_type_alloc ( natom * 3 )
#endif
      evb_frc_ready = .true.
   endif

!  +---------------------------------------------------------------------------+
!  |  Allocate memory for xnrg, xf, and xq in force                            |
!  +---------------------------------------------------------------------------+

   if( .not. evb_amber_ready ) then
      call evb_amber_alloc
      evb_amber_ready = .true.
   endif

!  +---------------------------------------------------------------------------+
!  |  Open EVB output file unit                                                |
!  +---------------------------------------------------------------------------+

   if( worldrank == 0 ) then

      write(6,'(A)') '| EVB data will be written to '//trim( adjustl( evbout ) )

      call amopen ( evb_unit, evbout, 'R', 'F', 'W' )

      write(evb_unit,'(/)')
#ifdef LES
      write(evb_unit, 888) '  [DYNAMICS TYPE]: ', 'pimd_'//trim( adjustl(evb_dyn) )
#else
      write(evb_unit, 888) '  [DYNAMICS TYPE]: ', trim( adjustl(evb_dyn) )
#endif
      write(evb_unit,'(/)')
      write(evb_unit, 999) '  [NBEAD]: ', ncopy
      write(evb_unit,'(/)')
      write(evb_unit,1000) '  [NEVB]   [NBIAS]   [NTW_EVB]   [DT]: ' &
                         , nevb, nbias, ntw_evb, dt
      write(evb_unit,'(/)')

      select case( trim( adjustl( evb_dyn) ) )

         case( "evb_map" ) 
            do n = 1, nbias
               write(evb_unit,2000) '  [MAPPING POTENTIAL]:   ni, nf, lambda ' &
                                    , bias_ndx(n,1), bias_ndx(n,2), lambda(n)
            enddo
         case( "egap_umb" ) 
            do n = 1, nbias
               write(evb_unit,3000) '  [NRG_GAP UMBRELLA]:   ni, nf, k, ezero ' &
                                    , bias_ndx(n,1), bias_ndx(n,2), k_umb(n) &
                                    , r0_umb(n)
            enddo

         case( "dbonds_umb", "qi_dbonds_pmf", "qi_dbonds_dyn" )
            do n = 1, nbias
               write(evb_unit,4000) '  [DBONDS UMBRELLA]:   iatom, jatom, katom, k, ezero ' &
                                    ,  dbonds_RC(n)%iatom,  dbonds_RC(n)%jatom &
                                    ,  dbonds_RC(n)%katom, k_umb(n), r0_umb(n)
            enddo

         case( "bond_umb", "qi_bond_pmf", "qi_bond_dyn" )
            do n = 1, nbias
               write(evb_unit,3000) '  [BOND UMBRELLA]:   iatom, jatom, k, ezero ' &
                                    ,  bond_RC(n)%iatom,  bond_RC(n)%jatom &
                                    ,  k_umb(n), r0_umb(n)
            enddo

      end select

#ifdef LES

!     if( do_ti_mass ) then

!        mm = 0
!        do n = 1, nti_mass
!           nn = ti_mass(n)%atm_ndx
!           do m = 1, ncopy
!              mm = mm + 1
!              ti_mass_ndx(mm) = bead_dcrypt(nn,m)
!           enddo
!        enddo

!        do n = 1, nti_mass
!           write(evb_unit,'(/)')
!           write(evb_unit,5000) '  [TI BY MASS]:   base atm_ndx, lambda, mass_init, ' &
!                               //  'mass_final', ti_mass(n)%atm_ndx, ti_mass(n)%lambda &
!                                , ti_mass(n)%init, ti_mass(n)%final
!        enddo

!        write(evb_unit,'(/)')
!        write(evb_unit,6000) '  [TI BY MASS ATM_NDX]', ti_mass_ndx(:)

!     endif

      if( itimass > 0 ) then
         write(evb_unit,'(/)')
         write(evb_unit,5000) '  [TI BY MASS]:   lambda ', clambda 
      endif

#endif /* LES */

   endif

  888 format( A/, A )
  999 format( A/, I8 )
 1000 format( A/, 3(2X,I8), 2X, F14.8 )
 2000 format( A/, 2I8, F14.8 )
 3000 format( A/, 2I8, F14.8, 2X, F14.8 )
 4000 format( A/, 3I8, F14.8, 2X, F14.8 )
#ifdef LES
 5000 format( A/, 2X, F14.8 )
!6000 format( A/, ( 6(2X,I10) ) )
#endif

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine evb_alloc



