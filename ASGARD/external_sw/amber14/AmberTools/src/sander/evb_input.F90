! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Read inputs for EVB calculation                                        |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#include "../include/assert.fh"
#include "../include/dprec.fh"

   subroutine evb_input 

   use evb_parm,  only: evb_dcrypt, C_evb, C_xch, xch_expdat &
                      , xch_gaussdat, bias_ndx, lambda, k_umb, r0_umb &
                      , nevb, nxch, nbias, nmorse, nmodvdw, ioe, dia_type &
                      , xch_type, evb_dyn, input_ready, ntw_evb, dbonds_RC &
                      , bond_RC, inc_dbonds_RC, out_RCdot 
!dEVB   use evb_parm,  only: xch_dcrypt, inc_bond_RC
   use evb_nml,   only: dia_size, xch_size, umb_size, nmorse_size, nmodvdw_size &
                      , nUFF_size, nAdihed_size, nAangle_size, xdg_size &
                      , nml_dia_shift, nml_xch_cnst &
                      , nml_xch_exp, nml_xch_gauss, nml_emap, nml_egap_umb &
                      , nml_dbonds_umb, nml_bond_umb, morsify, modvdw, UFF &
                      , Adihed, Aangle, nml_react_flux 
   use evb_amber, only: morse, vdw_mod
   use schlegel,  only: ndg, xdat_xdg, xdat_ts, xdat_min, alpha, dg_weight &
                      , gbasis_weight, svd_cond, dist_gauss, full_solve &
                      , nUFF, nUFF_dcrypt, nAdihed, nAdihed_dcrypt &
                      , nAangle, nAangle_dcrypt &
                      , UFF_scale &
                      , K_angle, theta &
                      , V, period, phase 
   use miller,    only: ti_mass, nti_mass, do_ti_mass
   use wigner,    only: nbias_ndx, p_space, minstep, nsample &
                      , s0_surface, s0_toler, RC_rmsd_toler &
                      , RC_RS_toler, RC_PS_toler, RC_type
   use evb_check, only: fdeps, debug_toler, evb_debug, what_size &
                      , ab_initio_gridfile, full_evb_debug, schlegel_debug &
                      , gradRC_debug, morsify_debug, mod_vdw_debug, xwarshel_debug &
                      , xgauss_debug, dbonds_debug, bond_debug, evb_pimd_debug

   use file_io_dat

   implicit none

#  include "../include/memory.h"

   type( nml_dia_shift  ) ::  dia_shift(dia_size)
   type( nml_xch_cnst   ) ::   xch_cnst(xch_size)
   type( nml_xch_exp    ) ::    xch_exp(xch_size)
   type( nml_xch_gauss  ) ::  xch_gauss(xch_size)
   type( nml_emap       ) ::       emap(umb_size)
   type( nml_egap_umb   ) ::   egap_umb(umb_size)
   type( nml_dbonds_umb ) :: dbonds_umb(umb_size)
   type( nml_bond_umb   ) ::   bond_umb(umb_size)
   type( nml_react_flux ) :: react_flux

   character(512) :: min_xfile(dia_size), ts_xfile(xch_size), xdg_xfile(xdg_size)  

   !  ..........................................................................

   integer :: n, nn, ios, alloc_error

   integer :: nmin_xfile, nts_xfile, ndg_xfile
   integer :: xch_check(4), dyn_check(4), dia_check(1)
   _REAL_  :: dgpt_alpha(xdg_size), dgpt_weight(xdg_size)
   _REAL_  :: gauss_weight(3)

   namelist / evb / nevb, nbias, nmorse, nmodvdw, nUFF, nAdihed, nAangle, dia_type, xch_type &
                  , evb_dyn, debug_toler, ioe, dia_shift, xch_cnst, xch_exp &
                  , xch_gauss, morsify, modvdw, UFF, UFF_scale, Adihed, Aangle &
                  , emap, egap_umb, dbonds_umb &
                  , bond_umb, ntw_evb, react_flux, min_xfile, ts_xfile &
                  , xdg_xfile, dist_gauss, dgpt_alpha, dgpt_weight &
                  , gauss_weight, svd_cond, ti_mass, evb_debug, fdeps &
                  , ab_initio_gridfile, full_solve, out_RCdot


   xch_check(:) = 0
   dyn_check(:) = 0
   dia_check(:) = 0

!  +---------------------------------------------------------------------------+
!  |  Initialize EVB input arrays                                              |
!  +---------------------------------------------------------------------------+

   dia_shift(:)%st  = 9999
   morsify(:)%iatom = 9999
   modvdw(:)%iatom  = 9999
   xch_cnst(:)%ist  = 9999
   xch_exp(:)%ist   = 9999
   xch_gauss(:)%ist = 9999

   ti_mass(:)%init = 0.0d0

   dgpt_alpha(:)    = 9999.9d0
   dgpt_weight(:)   = 9999.9d0
   gauss_weight(:)  = 9999.9d0
   gbasis_weight(:) = 9999.9d0

   min_xfile(:) = 'UNDEFINED'
   ts_xfile(:)  = 'UNDEFINED'
   xdg_xfile(:) = 'UNDEFINED'

   evb_debug%what(:) = 'XXX'

   emap(:)%ist         = 9999
   egap_umb(:)%ist     = 9999
   dbonds_umb(:)%iatom = 9999
   bond_umb(:)%iatom   = 9999

   nevb     = 0
   nxch     = 0
   nbias    = 0
   nmorse   = 0 
   nmodvdw  = 0 
   nUFF     = 0
   nAdihed  = 0
   nAangle  = 0
   nti_mass = 0
   ntw_evb  = 50
   UFF_scale = 1.0d0

   dia_type = 'force_field'
   xch_type = 'constant'
   evb_dyn  = 'groundstate'

   dist_gauss%stol = 1.0d-1
   dist_gauss%stype = 'react-product'
   dist_gauss%lin_solve = 'diis'
   dist_gauss%xfile_type = 'gaussian_fchk'

!  +---------------------------------------------------------------------------+
!  |  Read EVB input file                                                      |
!  +---------------------------------------------------------------------------+

   write(6,'(/)')
   write(6,'(A)') 'Reading EVB input file from ' // trim( adjustl( evbin ) )

   call amopen ( evb_unit, evbin, 'O', 'F', 'R' )
   read( evb_unit, nml = evb, IOSTAT = ios )
!  .............................................................................
!  :  Check for eof and read errors                                            :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   select case( ios )
      case(:-1)          ! end of file encountered
         write(6,'(A)')
         write(6,'(A)') 'ERROR: ' // trim( adjustl( evbin ) ) // ' does not exist'
         call mexit(6,1)
      case(1:)           ! error during read
         write(6,'(A)')
         write(6,'(A)') 'Error encountered while reading ' & 
                     // trim( adjustl( evbin ) )
         call mexit(6,1)
   end select

   close( evb_unit )

   call evb_keywrd( dia_type, xch_type, evb_dyn )

!  +---------------------------------------------------------------------------+
!  |  Maximum number of exchange terms in a symmetric matrix                   |
!  +---------------------------------------------------------------------------+

   nxch = nevb * ( nevb - 1 ) / 2

!  +---------------------------------------------------------------------------+
!  |  Allocate memory                                                          |
!  +---------------------------------------------------------------------------+

   if( .not. input_ready ) then
      if( dbonds_umb(nbias)%iatom /= 9999 ) inc_dbonds_RC = .true.
      call evb_input_alloc
      input_ready = .true.
   endif

!  +---------------------------------------------------------------------------+
!  |  Morsify Amber harmonic bond types                                        |
!  +---------------------------------------------------------------------------+

   do n = 1, nmorse
      morse(n)%iatom = morsify(n)%iatom
      morse(n)%jatom = morsify(n)%jatom
#ifdef LES
      morse(n)%D  = morsify(n)%D / dble( ncopy )
#else
      morse(n)%D  = morsify(n)%D 
#endif
      morse(n)%a  = morsify(n)%a 
      morse(n)%r0 = morsify(n)%r0
   enddo

!  +---------------------------------------------------------------------------+
!  |  Exclude certain van der waals interactions                               |
!  +---------------------------------------------------------------------------+

   do n = 1, nmodvdw
      vdw_mod(n)%iatom = modvdw(n)%iatom
      vdw_mod(n)%jatom = modvdw(n)%jatom
   enddo

!  +---------------------------------------------------------------------------+
!  |  Schlegel DG-EVB: include UFF VDW interaction in V_ii                     |
!  +---------------------------------------------------------------------------+

   do n = 1, nUFF
      nUFF_dcrypt(n,1) = UFF(n)%iatom
      nUFF_dcrypt(n,2) = UFF(n)%jatom
   enddo

!  +---------------------------------------------------------------------------+
!  |  Schlegel DG-EVB: include Amber torsion interaction in V_ii               |
!  +---------------------------------------------------------------------------+

   do n = 1, nAdihed
      nAdihed_dcrypt(n,1) = Adihed(n)%iatom
      nAdihed_dcrypt(n,2) = Adihed(n)%jatom
      nAdihed_dcrypt(n,3) = Adihed(n)%katom
      nAdihed_dcrypt(n,4) = Adihed(n)%latom

      V(n)      = Adihed(n)%V
      period(n) = Adihed(n)%period
      phase(n)  = Adihed(n)%phase
   enddo

!  +---------------------------------------------------------------------------+
!  |  Schlegel DG-EVB: include Amber angle interaction in V_ii                 |
!  +---------------------------------------------------------------------------+

   do n = 1, nAangle
      nAangle_dcrypt(n,1) = Aangle(n)%iatom
      nAangle_dcrypt(n,2) = Aangle(n)%jatom
      nAangle_dcrypt(n,3) = Aangle(n)%katom

      K_angle(n) = Aangle(n)%K
      theta(n) = Aangle(n)%theta
   enddo

!  +---------------------------------------------------------------------------+
!  |  Diabatic state relative energy shifts                                    |
!  +---------------------------------------------------------------------------+

   do n = 1, nevb
      nn = dia_shift(n)%st
      C_evb(nn) = dia_shift(n)%nrg_offset
   enddo

!  +---------------------------------------------------------------------------+
!  |  Exchange type                                                            |
!  +---------------------------------------------------------------------------+

   select case( trim( adjustl( xch_type ) ) )

!  .............................................................................
!  :  Constant exchange term                                                   :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
      case( "constant" ) 
         do n = 1, nxch
            evb_dcrypt(n,1) = xch_cnst(n)%ist 
            evb_dcrypt(n,2) = xch_cnst(n)%jst 
            C_xch(n) = xch_cnst(n)%xcnst
!dEVB       xch_dcrypt( evb_dcrypt(n,1) ) = n
!dEVB       xch_dcrypt( evb_dcrypt(n,2) ) = n
            xch_check(1) = 1
         enddo

!  .............................................................................
!  :  Warshel's exponential functional form for the exchange term              :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
      case( "exp" ) 
         do n = 1, nxch
            evb_dcrypt(n,1) = xch_exp(n)%ist
            evb_dcrypt(n,2) = xch_exp(n)%jst
            xch_expdat(n) = xch_exp(n)
!dEVB       xch_dcrypt( evb_dcrypt(n,1) ) = n
!dEVB       xch_dcrypt( evb_dcrypt(n,2) ) = n
            xch_check(2) = 1
         enddo

!  .............................................................................
!  :  Gaussian functional form for the exchange term                           :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
      case( "gauss" ) 
         do n = 1, nxch
            evb_dcrypt(n,1) = xch_gauss(n)%ist
            evb_dcrypt(n,2) = xch_gauss(n)%jst
            xch_gaussdat(n) = xch_gauss(n)
!dEVB       xch_dcrypt( evb_dcrypt(n,1) ) = n
!dEVB       xch_dcrypt( evb_dcrypt(n,2) ) = n
            xch_check(3) = 1
         enddo

!  .............................................................................
!  :  DG-EVB: determine the number of external files                           :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
      case( "dist_gauss" )
         nts_xfile = 0
         do n = 1, dia_size
            if( trim( adjustl( ts_xfile(n) ) ) /= 'UNDEFINED' ) &
            nts_xfile = nts_xfile + 1
         enddo
         if( nts_xfile /= nxch ) then
            write(6,'(2(/))')
            write(6,'(2(A,I8,A))') 'ERROR: the no. of  TS external ' &
                                // 'files (nts_xfile = ', nts_xfile, ')' &
                                 , ' is not equal to nxch = ', nxch, '.'
            call mexit(6,1)
         endif
         nmin_xfile = 0
         do n = 1, dia_size
            if( trim( adjustl( min_xfile(n) ) ) /= 'UNDEFINED' ) &
            nmin_xfile = nmin_xfile + 1
         enddo
         if( nmin_xfile /= nevb ) then
            write(6,'(2(/))')
            write(6,'(2(A,I8,A))') 'ERROR: the no. of  PES minimum external ' &
                                // 'files (nmin_xfile = ', nmin_xfile, ')' &
                                 , ' is not equal to nevb = ', nevb, '.'
            call mexit(6,1)
         endif
         ndg_xfile = 0
         do n = 1, xdg_size
            if( trim( adjustl( xdg_xfile(n) ) ) /= 'UNDEFINED' ) &
            ndg_xfile = ndg_xfile + 1
         enddo

!  .............................................................................
!  :  Order the DG files as TS, MIN, additional points                         :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
         allocate( xdat_min(nevb), stat = alloc_error )
         REQUIRE( alloc_error == 0 )
         allocate( xdat_ts(nxch), stat = alloc_error )
         REQUIRE( alloc_error == 0 )
         allocate( xdat_xdg(nts_xfile+nevb+ndg_xfile), stat = alloc_error )
         REQUIRE( alloc_error == 0 )

         if(.not. allocated(evb_dcrypt) ) then
            allocate( evb_dcrypt(1,2), stat = alloc_error )
            REQUIRE( alloc_error == 0 )
         endif

         nn = 0
         do n = 1, nxch
            nn = nn + 1
            xdat_ts(n)%filename   = ts_xfile(n)
            xdat_xdg(nn)%filename = ts_xfile(n)
         enddo
         do n = 1, nevb
            nn = nn + 1
            xdat_min(n)%filename  = min_xfile(n)
            xdat_xdg(nn)%filename = min_xfile(n)
            evb_dcrypt(1,n) = n
         enddo
         do n = 1, ndg_xfile
            nn = nn + 1
            xdat_xdg(nn)%filename = xdg_xfile(n)
         enddo

         ndg = nts_xfile + nmin_xfile + ndg_xfile

         write(6,'(A,I8)') 'No. of DG at TS      = ', nts_xfile
         write(6,'(A,I8)') 'No. of DG at minima  = ', nevb
         write(6,'(A,I8)') 'No. of additional DG = ', ndg_xfile
         write(6,'(A,I8)') 'Total DG points      = ', ndg

         allocate( alpha(ndg), stat = alloc_error )
         REQUIRE( alloc_error == 0 )

         alpha(:) = dgpt_alpha(1)
         do n = 1, ndg
            if( dgpt_alpha(n) /= 9999.9d0 ) alpha(n) = dgpt_alpha(n)
         enddo

         allocate( dg_weight(ndg), stat = alloc_error )
         REQUIRE( alloc_error == 0 )

         dg_weight(:) = dgpt_weight(1)
         do n = 1, ndg
            if( dgpt_weight(n) /= 9999.9d0 ) then
               dg_weight(n) = dgpt_weight(n)
            else
               dg_weight(n) = 1.0d0
            endif
         enddo

         do n = 1, 3
            if( gauss_weight(n) /= 9999.9d0 ) gbasis_weight(n) = gauss_weight(n) 
         enddo

   end select 

!  +---------------------------------------------------------------------------+
!  |  Dynamics type                                                            |
!  +---------------------------------------------------------------------------+

   select case( trim( adjustl( evb_dyn ) ) )

! 12062008: for outputting RC value when performing groundstate dynamics

!  .............................................................................
!  :  Umbrella sampling: difference of 2 bonds RC                              :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

!        allocate( dbonds_RC(nbias), stat = alloc_error )
!        REQUIRE( alloc_error == 0 )
      case( "groundstate" )
!        write(6,*) 'nbias = ', nbias
!        write(6,*) "evb_input: this allocation will break others"
         if( inc_dbonds_RC ) then
            do n = 1, nbias
               dbonds_RC(n)%iatom = dbonds_umb(n)%iatom
               dbonds_RC(n)%jatom = dbonds_umb(n)%jatom
               dbonds_RC(n)%katom = dbonds_umb(n)%katom
            enddo
         endif

!  .............................................................................
!  :  Reactive flux                                                            :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
      case( "react_flux" )
         nbias_ndx     = react_flux%nbias_ndx
         p_space       = react_flux%p_space
         minstep       = react_flux%minstep
         nsample       = react_flux%nsample
         s0_surface    = react_flux%s0_surface
         s0_toler      = react_flux%s0_toler
         RC_rmsd_toler = react_flux%RC_rmsd_toler
         RC_RS_toler   = react_flux%RC_RS_toler
         RC_PS_toler   = react_flux%RC_PS_toler
         RC_type       = react_flux%RC_type

         do n = 1, nbias
!  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
!  [] Mapping potential                                                        :
!  `````````````````````````````````````````````````````````````````````````````
            nn = emap(n)%ist
            if( nn /= 9999 ) then
               bias_ndx(n,1) = emap(n)%ist
               bias_ndx(n,2) = emap(n)%jst
            endif
!  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
!  [] Energy gap RC harmonic umbrella                                          :
!  `````````````````````````````````````````````````````````````````````````````
            nn = egap_umb(n)%ist
            if( nn /= 9999 ) then
               bias_ndx(n,1) = egap_umb(n)%ist
               bias_ndx(n,2) = egap_umb(n)%jst
            endif
!  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
!  [] Difference of 2 bonds RC harmonic umbrella sampling                      :
!  `````````````````````````````````````````````````````````````````````````````
            nn = dbonds_umb(n)%iatom
            if( nn /= 9999 ) then
               dbonds_RC(n)%iatom = dbonds_umb(n)%iatom
               dbonds_RC(n)%jatom = dbonds_umb(n)%jatom
               dbonds_RC(n)%katom = dbonds_umb(n)%katom
            endif
!  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
!  [] Bond RC harmonic umbrella sampling                                       :
!  `````````````````````````````````````````````````````````````````````````````
            nn = bond_umb(n)%iatom
            if( nn /= 9999 ) then
               bond_RC(n)%iatom = bond_umb(n)%iatom
               bond_RC(n)%jatom = bond_umb(n)%jatom
            endif
         enddo

!  .............................................................................
!  :  Warshel's mapping potential                                              :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
      case( "evb_map" ) 
         do n = 1, nbias
            bias_ndx(n,1) = emap(n)%ist
            bias_ndx(n,2) = emap(n)%jst
            lambda(n) = emap(n)%lambda
         enddo
!        do n = 1, nbias
!           dbonds_RC(n)%iatom = dbonds_umb(n)%iatom
!           dbonds_RC(n)%jatom = dbonds_umb(n)%jatom
!           dbonds_RC(n)%katom = dbonds_umb(n)%katom
!        enddo
!        if( dbonds_umb(nbias)%iatom /= 9999 ) inc_dbonds_RC = .true.
!        do n = 1, nbias
!           bond_RC(n)%iatom = bond_umb(n)%iatom
!           bond_RC(n)%jatom = bond_umb(n)%jatom
!        enddo
!        if( bond_umb(nbias)%iatom /= 9999 ) inc_bond_RC = .true.
         dyn_check(1) = 1

!  .............................................................................
!  :  Energy gap RC harmonic umbrella                                          :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
      case( "egap_umb" ) 
         do n = 1, nbias
            bias_ndx(n,1) = egap_umb(n)%ist
            bias_ndx(n,2) = egap_umb(n)%jst
                 k_umb(n) = egap_umb(n)%k
                r0_umb(n) = egap_umb(n)%ezero
         enddo
         dyn_check(2) = 1

!  .............................................................................
!  :  Umbrella sampling: difference of 2 bonds RC                              :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
      case( "dbonds_umb", "qi_dbonds_pmf", "qi_dbonds_dyn" )
         do n = 1, nbias
            dbonds_RC(n)%iatom = dbonds_umb(n)%iatom
            dbonds_RC(n)%jatom = dbonds_umb(n)%jatom
            dbonds_RC(n)%katom = dbonds_umb(n)%katom
                      k_umb(n) = dbonds_umb(n)%k
                     r0_umb(n) = dbonds_umb(n)%ezero
         enddo
         dyn_check(3) = 1

!  .............................................................................
!  :  Umbrella sampling: distance RC                                           :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
      case( "bond_umb", "qi_bond_pmf", "qi_bond_dyn" )
         do n = 1, nbias
            bond_RC(n)%iatom = bond_umb(n)%iatom
            bond_RC(n)%jatom = bond_umb(n)%jatom
                    k_umb(n) = bond_umb(n)%k
                   r0_umb(n) = bond_umb(n)%ezero
         enddo
         dyn_check(4) = 1

   end select

!  +---------------------------------------------------------------------------+
!  |  Check for TI by mass                                                     |
!  +---------------------------------------------------------------------------+

   if( ti_mass(1)%init /= 0.0d0 ) then
      nti_mass = 1
      do_ti_mass = .true.
      if( ti_mass(2)%init /= 0.0d0 ) nti_mass = 2
   endif

!  +---------------------------------------------------------------------------+
!  |  Determine debugging options                                              |
!  +---------------------------------------------------------------------------+

   do n = 1, what_size
      select case( trim( adjustl( evb_debug%what(n) ) ) )
         case( "full"     ) 
            full_evb_debug = .true.
         case( "schlegel" ) 
            schlegel_debug = .true.
         case( "gradRC"   ) 
            gradRC_debug   = .true.
         case( "morsify"  ) 
            morsify_debug  = .true.
         case( "mod_vdw"  ) 
            mod_vdw_debug  = .true.
         case( "xwarshel" ) 
            xwarshel_debug = .true.
         case( "xgauss"   ) 
            xgauss_debug   = .true.
         case( "dbonds"   ) 
            dbonds_debug   = .true.
         case( "bond"     ) 
            bond_debug     = .true.
         case( "evb_pimd" ) 
            evb_pimd_debug = .true.
      end select
   enddo

!  +---------------------------------------------------------------------------+
!  |  Check for EVB variable dependences                                       |
!  +---------------------------------------------------------------------------+

   call check_input ( xch_check, dyn_check, dia_check, xch_type, evb_dyn )

!  +---------------------------------------------------------------------------+
!  |  Output EVB options                                                       |
!  +---------------------------------------------------------------------------+

   call evb_outoptions


!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine evb_input


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Allocate storage space for evb_input                                   |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   subroutine evb_input_alloc 

   use evb_parm,  only: evb_dcrypt, C_evb, C_xch &
                      , xch_expdat, xch_gaussdat, bias_ndx, lambda &
                      , k_umb, r0_umb, nevb, nxch &
                      , nmorse, nmodvdw, nbias, xch_type, evb_dyn &
                      , dbonds_RC, bond_RC, inc_dbonds_RC
!dEVB   use evb_parm,  only: xch_dcrypt
   use evb_amber, only: morse, k_harm, r0_harm, vdw_mod
   use schlegel,  only: nUFF, nUFF_dcrypt, nAdihed, nAdihed_dcrypt &
                      , nAangle, nAangle_dcrypt, K_angle, theta, V, period, phase


   implicit none


   !...........................................................................

   integer :: alloc_error

!  +---------------------------------------------------------------------------+
!  |  Morse parameters                                                         |
!  +---------------------------------------------------------------------------+

   if( nmorse > 0 ) then
      allocate( morse(nmorse), stat = alloc_error )
      REQUIRE( alloc_error == 0 )
      allocate( k_harm(nmorse), stat = alloc_error )
      REQUIRE( alloc_error == 0 )
      allocate( r0_harm(nmorse), stat = alloc_error )
      REQUIRE( alloc_error == 0 )
   endif

!  +---------------------------------------------------------------------------+
!  |  VDW exclusion                                                            |
!  +---------------------------------------------------------------------------+

   if( nmodvdw > 0 ) then
      allocate( vdw_mod(nmodvdw), stat = alloc_error )
      REQUIRE( alloc_error == 0 )
   endif

!  +---------------------------------------------------------------------------+
!  |  Schlegel DG-EVB: include UFF VDW interaction in V_ii                     |
!  +---------------------------------------------------------------------------+

   if( nUFF > 0 ) then
      allocate( nUFF_dcrypt(nUFF,3), stat = alloc_error )
      REQUIRE( alloc_error == 0 )
   endif

!  +---------------------------------------------------------------------------+
!  |  Schlegel DG-EVB: include Amber torsion interaction in V_ii               |
!  +---------------------------------------------------------------------------+

   if( nAdihed > 0 ) then
      allocate( nAdihed_dcrypt(nAdihed,5), stat = alloc_error )
      REQUIRE( alloc_error == 0 )
      allocate( V(nAdihed), period(nAdihed), phase(nAdihed), stat = alloc_error )
      REQUIRE( alloc_error == 0 )
   endif

!  +---------------------------------------------------------------------------+
!  |  Schlegel DG-EVB: include Amber angle interaction in V_ii                 |
!  +---------------------------------------------------------------------------+

   if( nAangle > 0 ) then
      allocate( nAangle_dcrypt(nAangle,4), stat = alloc_error )
      REQUIRE( alloc_error == 0 )
      allocate( K_angle(nAangle), theta(nAangle), stat = alloc_error )
      REQUIRE( alloc_error == 0 )
   endif

!  +---------------------------------------------------------------------------+
!  |  Diabatic state relative energy shifts                                    |
!  +---------------------------------------------------------------------------+

   allocate( C_evb(nevb), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   C_evb(:) = 0.d0

!  +---------------------------------------------------------------------------+
!  |  Exchange type                                                            |
!  +---------------------------------------------------------------------------+

   select case( trim( adjustl( xch_type ) ) )

!  :,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,:
!  :  Constant exchange term                                                   :
!  :''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''':

      case( "constant" ) 
         allocate( evb_dcrypt(nxch,2), stat = alloc_error )
         REQUIRE( alloc_error == 0 )
!dEVB    allocate( xch_dcrypt(nevb), stat = alloc_error )
!dEVB    REQUIRE( alloc_error == 0 )
         allocate( C_xch(nxch), stat = alloc_error )
         REQUIRE( alloc_error == 0 )
         C_xch(:) = 0.d0

!  :,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,:
!  :  Exponential functional form for exchange term                            :
!  :''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''':

      case( "exp" ) 
         allocate( evb_dcrypt(nxch,2), stat = alloc_error )
         REQUIRE( alloc_error == 0 )
!dEVB    allocate( xch_dcrypt(nevb), stat = alloc_error )
!dEVB    REQUIRE( alloc_error == 0 )
         allocate( xch_expdat(nxch), stat = alloc_error )
         REQUIRE( alloc_error == 0 )

!  :,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,:
!  :  Gaussian functional form for exchange term                               :
!  :''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''':

      case( "gauss" ) 
         allocate( evb_dcrypt(nxch,2), stat = alloc_error )
         REQUIRE( alloc_error == 0 )
!dEVB    allocate( xch_dcrypt(nevb), stat = alloc_error )
!dEVB    REQUIRE( alloc_error == 0 )
         allocate( xch_gaussdat(nxch), stat = alloc_error )
         REQUIRE( alloc_error == 0 )

   end select

!  +---------------------------------------------------------------------------+
!  |  Dynamics type                                                            |
!  +---------------------------------------------------------------------------+

   select case( trim( adjustl( evb_dyn ) ) )

!  :,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,:
!  :  "groundstate but with RC output                              :
!  :''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''':

      case( "groundstate" )
         if( inc_dbonds_RC ) then
            allocate( dbonds_RC(nbias), stat = alloc_error )
            REQUIRE( alloc_error == 0 )
         endif

!  :,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,:
!  :  Reactive flux                                                            :
!  :''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''':

      case( "react_flux" )
         allocate( bias_ndx(nbias,2), stat = alloc_error )
         REQUIRE( alloc_error == 0 )
         allocate( dbonds_RC(nbias), stat = alloc_error )
         REQUIRE( alloc_error == 0 )
         allocate( bond_RC(nbias), stat = alloc_error )
         REQUIRE( alloc_error == 0 )

!  :,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,:
!  :  Mapping potential parameter                                              :
!  :''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''':

      case( "evb_map" ) 
         allocate( lambda(nbias), stat = alloc_error )
         REQUIRE( alloc_error == 0 )
         allocate( bias_ndx(nbias,2), stat = alloc_error )
         REQUIRE( alloc_error == 0 )
         allocate( dbonds_RC(nbias), stat = alloc_error )
         REQUIRE( alloc_error == 0 )

!  :,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,:
!  :  Energy gap RC harmonic umbrella parameter                                :
!  :''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''':

      case( "egap_umb" ) 
         allocate( bias_ndx(nbias,2), stat = alloc_error )
         REQUIRE( alloc_error == 0 )
         allocate( k_umb(nbias), stat = alloc_error )
         REQUIRE( alloc_error == 0 )
         allocate( r0_umb(nbias), stat = alloc_error )
         REQUIRE( alloc_error == 0 )

!  :,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,:
!  :  Difference of 2 bonds RC harmonic umbrella parameter         :
!  :''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''':

      case( "dbonds_umb", "qi_dbonds_pmf", "qi_dbonds_dyn" ) 
         allocate( dbonds_RC(nbias), stat = alloc_error )
         REQUIRE( alloc_error == 0 )
         allocate( k_umb(nbias), stat = alloc_error )
         REQUIRE( alloc_error == 0 )
         allocate( r0_umb(nbias), stat = alloc_error )
         REQUIRE( alloc_error == 0 )

!  :,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,:
!  :  Bond RC harmonic umbrella parameter                                      :
!  :''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''':

      case( "bond_umb", "qi_bond_pmf", "qi_bond_dyn" )
         allocate( bond_RC(nbias), stat = alloc_error )
         REQUIRE( alloc_error == 0 )
         allocate( k_umb(nbias), stat = alloc_error )
         REQUIRE( alloc_error == 0 )
         allocate( r0_umb(nbias), stat = alloc_error )
         REQUIRE( alloc_error == 0 )

   end select

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine evb_input_alloc


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Deallocate storage space for evb_input                                 |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   subroutine evb_input_dealloc

   use evb_parm,  only: evb_dcrypt, C_evb, C_xch &
                      , xch_expdat, xch_gaussdat, bias_ndx, lambda &
                      , k_umb, r0_umb &
                      , xch_type, evb_dyn &
                      , nmorse &
                      , dbonds_RC, bond_RC

!dEVB   use evb_parm,  only: xch_dcrypt

   use evb_amber, only: morse, k_harm, r0_harm

   implicit none

   !...........................................................................

   integer :: dealloc_error


!  +---------------------------------------------------------------+
!  |  Morse potential                                              |
!  +---------------------------------------------------------------+

   if( nmorse > 0 ) then 

      deallocate( morse, k_harm, r0_harm, stat = dealloc_error )
      REQUIRE( dealloc_error == 0 )

   endif

!  +---------------------------------------------------------------------------+
!  |  Diabatic state relative energy shifts                                    |
!  +---------------------------------------------------------------------------+

   deallocate( C_evb, stat = dealloc_error )
   REQUIRE( dealloc_error == 0 )

!  +---------------------------------------------------------------------------+
!  |  Constant exchange term                                                   |
!  +---------------------------------------------------------------------------+

   if( trim( adjustl( xch_type ) ) == 'constant' ) then

      deallocate( evb_dcrypt, C_xch, stat = dealloc_error )
      REQUIRE( dealloc_error == 0 )

!dEVB deallocate( xch_dcrypt, stat = dealloc_error )
!dEVB REQUIRE( dealloc_error == 0 )

   endif

!  +---------------------------------------------------------------------------+
!  |  Exponential functional form for exchange term                            |
!  +---------------------------------------------------------------------------+

   if( trim( adjustl( xch_type ) ) == 'exp' ) then

      deallocate( evb_dcrypt, xch_expdat, stat = dealloc_error )
      REQUIRE( dealloc_error == 0 )

   endif

!  +---------------------------------------------------------------------------+
!  |  Gaussian functional form for exchange term                               |
!  +---------------------------------------------------------------------------+

   if( trim( adjustl( xch_type ) ) == 'gauss' ) then

      deallocate( evb_dcrypt, xch_gaussdat, stat = dealloc_error )
      REQUIRE( dealloc_error == 0 )

   endif

   select case( trim( adjustl( evb_dyn ) ) )

!  +---------------------------------------------------------------------------+
!  |  Mapping potential parameter                                              |
!  +---------------------------------------------------------------------------+

      case( "evb_map" ) 

         deallocate( lambda, bias_ndx, dbonds_RC, stat = dealloc_error )
         REQUIRE( dealloc_error == 0 )

!  +---------------------------------------------------------------------------+
!  |  Energy gap RC harmonic umbrella sampling                                 |
!  +---------------------------------------------------------------------------+

      case( "egap_umb" ) 

         deallocate( bias_ndx, k_umb, r0_umb, stat = dealloc_error )
         REQUIRE( dealloc_error == 0 )

!  +---------------------------------------------------------------------------+
!  |  Difference of 2 bonds RC harmonic umbrella parameter                     |
!  +---------------------------------------------------------------------------+

      case( "dbonds_umb", "qi_dbonds_pmf", "qi_dbonds_dyn" )

         deallocate( dbonds_RC, k_umb, r0_umb, stat = dealloc_error )
         REQUIRE( dealloc_error == 0 )

!  +---------------------------------------------------------------------------+
!  |  Bond RC harmonic umbrella parameter                                      |
!  +---------------------------------------------------------------------------+

      case( "bond_umb", "qi_bond_pmf", "qi_bond_dyn" )

         deallocate( bond_RC, k_umb, r0_umb, stat = dealloc_error )
         REQUIRE( dealloc_error == 0 )

   end select 

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine evb_input_dealloc


