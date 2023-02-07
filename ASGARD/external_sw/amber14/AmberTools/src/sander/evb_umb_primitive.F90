! <compile=optimized>
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  EVB_UMB                                                                |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#include "../include/assert.fh"
#include "../include/dprec.fh"

   subroutine evb_umb_primitive ( f, q, mass, natom, istart, iend )

   use constants, only: kB
   use evb_parm,  only: k_umb, r0_umb, evb_dyn, nbias, dbonds_RC, bond_RC &
                      , out_RCdot, egapRC
   use evb_data,  only: evb_frc, evb_bias, RCdot
#ifdef DEBUG_EVB
   use evb_check, only: full_evb_debug, dbonds_debug, bond_debug
#endif /* DEBUG_EVB */
   use file_io_dat
#ifdef LES
   use pimd_vars, only: nbead, nbead_inv
   use evb_pimd,  only: bead_dcrypt, natomCL
!  use miller,    only: i_qi, div_ndx
#endif

   implicit none

#  include "parallel.h"
#  include "../include/md.h"

   integer, intent(in   ) :: natom, istart, iend
   _REAL_ , intent(in   ) :: q(3,natom)
   _REAL_ , intent(in   ) :: mass(natom)
   _REAL_ , intent(inout) :: f(3,natom)

   !  ..........................................................................

   integer :: m, n, nn
   integer :: idx, jdx, kdx
#ifdef LES
   _REAL_  :: q_centroid(3,natomCL), fharm(3,natomCL), evb_fbias(3,natomCL,nbias) &
            , f_les(3,natom)
   integer :: mm
#else
   _REAL_  :: fharm(3,natom), evb_fbias(3,natom,nbias) 
#endif /* LES */
   _REAL_  :: RC, dr(3), fr(3), rij, rkj, rij_inv, rkj_inv, pi
   intrinsic :: sqrt, acos

!  +---------------------------------------------------------------------------+
!  |  Apply umbrella on the centroid coordinates.  Construct the centroid      |
!  |  from the Cartesians of the beads.                                        |
!  +---------------------------------------------------------------------------+

#ifdef LES
!  write(6,*) '>>> nbead, natomCL, natom = ', nbead, natomCL, natom
!  write(6,*) '>>> nbead_inv = ', nbead_inv
   q_centroid(:,:) = 0.0d0
   do n = 1, nbead 
      do m = 1, natomCL
         mm = bead_dcrypt(m,n)
         q_centroid(:,m) = q_centroid(:,m) + q(:,mm)

!        write(6,*) '>>>', m, n, mm

      enddo
   enddo
   q_centroid(:,:) = q_centroid(:,:) * nbead_inv

!  write(100,'((6F12.7))') q_centroid(:,:)
!  write(200,'((6F12.7))') q(:,:)

! +++++++++++++++++++++
!  do n = 1, nbead
!     q_bead(:,:) = 0.0d0
!     do m = 1, natomCL
!        mm = bead_dcrypt(m,n)
!        q_bead(:,m) = q(:,mm)
!     enddo
!     write(200,*) 'bead = ', n 
!     write(200,'((6F12.7))') q_bead(:,:)
!  enddo

#endif /* LES */

   select case( trim( adjustl( evb_dyn) ) )

!  +---------------------------------------------------------------------------+
!  |  Difference of 2 bonds RC harmonic umbrella sampling                      |
!  |                                                                           |
!  |  /\ = r_ij - r_kj                                                         |
!  |                                                                           |
!  |  V_umb(eta) = V_0 + 0.5 k_evb * ( /\ - eta )^2                            |
!  |                                                                           |
!  |  dV_umb / dR = ( dV_umb / d/\ ) * ( d/\ / dR )                            |
!  |              = dV_0 / dR + k_evb * ( /\ - eta ) * d/dR ( r_ij - r_kj )    |
!  +---------------------------------------------------------------------------+

      case( "dbonds_umb", "qi_dbonds_pmf", "qi_dbonds_dyn" )

         fharm(:,:) = 0.0d0

         do n = 1, nbias

            idx = dbonds_RC(n)%iatom
            jdx = dbonds_RC(n)%jatom
            kdx = dbonds_RC(n)%katom

! write(6,*) 'idx, jdx, kdx = ', idx, jdx, kdx

#ifdef LES
            do nn = 1, 3
               dr(nn) = q_centroid(nn,idx) - q_centroid(nn,jdx)
            enddo 
#else
            do nn = 1, 3
               dr(nn) = q(nn,idx) - q(nn,jdx)
            enddo
#endif /* LES */
            rij = sqrt( dr(1)**2 + dr(2)**2 + dr(3)**2 )
            rij_inv = 1.0d0 / rij

            do nn = 1, 3
               fr(nn) = dr(nn) * rij_inv
               fharm(nn,idx) = fharm(nn,idx) - fr(nn)
               fharm(nn,jdx) = fharm(nn,jdx) + fr(nn)
            enddo 

#ifdef LES
            do nn = 1, 3 
               dr(nn)  = q_centroid(nn,kdx) - q_centroid(nn,jdx)
            enddo
#else
            do nn = 1, 3
               dr(nn)  = q(nn,kdx) - q(nn,jdx)
            enddo
#endif /* LES */
            rkj = sqrt( dr(1)**2 + dr(2)**2 + dr(3)**2 )
            rkj_inv = 1.0d0 / rkj

            do nn = 1, 3
               fr(nn) = dr(nn) * rkj_inv
               fharm(nn,kdx) = fharm(nn,kdx) + fr(nn)
               fharm(nn,jdx) = fharm(nn,jdx) - fr(nn)
            enddo

#ifdef LES
!           if( i_qi > 0 ) then
!              do m = 1, natomCL
!                 mm = bead_dcrypt( m, div_ndx(n) )
!                 gradRC(m,n) = fharm(mm) 
!              enddo
!           endif
#endif
            RC = rij - rkj

! write(6,*) 'rij = ', rij
! write(6,*) 'rkj = ', rkj
! write(6,*) 'RC  = ', RC
  

            evb_bias%RC(n) = RC
            evb_bias%nrg_bias(n) = 0.50d0 * k_umb(n) * ( RC - r0_umb(n) )**2
#ifdef LES
            evb_fbias(:,:,n) = k_umb(n) * ( RC - r0_umb(n) ) * fharm(:,:) &
                             * nbead_inv
#else
            evb_fbias(:,:,n) = k_umb(n) * ( RC - r0_umb(n) ) * fharm(:,:)
#endif /* LES */

#ifdef DEBUG_EVB
            if( dbonds_debug .or. full_evb_debug ) &
               call dbonds_anal2num ( q, evb_fbias(:,:,n), natom*3 )
#endif

         enddo

!  +---------------------------------------------------------------------------+
!  |  Distribute the umbrella force onto the Cartesian beads                   |
!  +---------------------------------------------------------------------------+

         do n = 1, nbias
            evb_frc%evb_nrg  = evb_frc%evb_nrg + evb_bias%nrg_bias(n)
#ifdef LES
            f_les(:,:) = 0.0d0
            do nn = 1, nbead
               do m = 1, natomCL
                  mm = bead_dcrypt(m,nn)
                  f_les(:,mm) = f_les(:,mm) + evb_fbias(:,m,n)
               enddo
            enddo
            do m = istart, iend
               f(:,m) = f(:,m) + f_les(:,m)
            enddo
#else
            do m = istart, iend
               f(:,m) = f(:,m) + evb_fbias(:,m,n)
            enddo
#endif /* LES */
         enddo

!  +---------------------------------------------------------------------------+
!  |  Bond RC harmonic umbrella sampling                                       |
!  |                                                                           |
!  |  /\ = r_ij                                                                |
!  |                                                                           |
!  |  V_umb(eta) = V_0 + 0.5 k_evb * ( /\ - eta )^2                            |
!  |                                                                           |
!  |  dV_umb / dR = ( dV_umb / d/\ ) * ( d/\ / dR )                            |
!  |              = dV_0 / dR + k_evb * ( /\ - eta ) * d/dR ( r_ij )           |
!  +---------------------------------------------------------------------------+

      case( "bond_umb", "qi_bond_pmf", "qi_bond_dyn" )

         fharm(:,:) = 0.0d0

         do n = 1, nbias

            idx = bond_RC(n)%iatom
            jdx = bond_RC(n)%jatom

#ifdef LES
            do nn = 1, 3
               dr(nn) = q_centroid(nn,idx) - q_centroid(nn,jdx)
            enddo

! write(6,*) '>>> idx, jdx = ', idx, jdx
! write(6,*) '>>> dr = ', dr(:)

#else
            do nn = 1, 3
               dr(nn) = q(nn,idx) - q(nn,jdx)
            enddo
#endif

            rij = sqrt( dr(1)**2 + dr(2)**2 + dr(3)**2 )
            rij_inv = 1.0d0 / rij

            do nn = 1, 3
               fr(nn) = dr(nn) * rij_inv
               fharm(nn,idx) = fharm(nn,idx) - fr(nn)
               fharm(nn,jdx) = fharm(nn,jdx) + fr(nn)
            enddo

#ifdef LES
!           if( i_qi > 0 ) then
!              do m = 1, natomCL
!                 mm = bead_dcrypt( m, div_ndx(n) )
!                 gradRC(m,n) = fharm(mm)
!              enddo
!           endif
#endif

            RC = rij 

!write(6,*) '>>> RC = ', RC

            evb_bias%RC(n) = RC
            evb_bias%nrg_bias(n) = 0.50d0 * k_umb(n) * ( RC - r0_umb(n) )**2
#ifdef LES
            evb_fbias(:,:,n) = k_umb(n) * ( RC - r0_umb(n) ) * fharm(:,:) &
                             * nbead_inv
#else
            evb_fbias(:,:,n) = k_umb(n) * ( RC - r0_umb(n) ) * fharm(:,:) 
#endif /* LES */

#ifdef DEBUG_EVB
            if( bond_debug .or. full_evb_debug ) &
               call bond_anal2num ( q, evb_fbias(:,:,n), natom*3 )
#endif

         enddo

!  +---------------------------------------------------------------------------+
!  |  Distribute the umbrella force onto the Cartesian beads                   |
!  +---------------------------------------------------------------------------+

         do n = 1, nbias
            evb_frc%evb_nrg  = evb_frc%evb_nrg + evb_bias%nrg_bias(n)
#ifdef LES
            f_les(:,:) = 0.0d0
            do nn = 1, nbead
               do m = 1, natomCL
                  mm = bead_dcrypt(m,nn)
                  f_les(:,mm) = f_les(:,mm) + evb_fbias(:,m,n)
               enddo
            enddo
            do m = istart, iend
               f(:,m) = f(:,m) + f_les(:,m)
            enddo
#else
            do m = istart, iend
               f(:,m) = f(:,m) + evb_fbias(:,m,n)
            enddo
#endif /* LES */
         enddo

   end select

!  write(6,*) 'exiting now!'
!  call mexit(6,1) 

   if( out_RCdot .and. .not. egapRC ) then

      pi = acos( -1.0d0 )
      RCdot = 0.0d0
#ifdef LES
      do n = 1, natomCL
         nn = bead_dcrypt(n,1)
         RCdot = RCdot + ( fharm(1,n)**2 + fharm(2,n)**2 &
                         + fharm(3,n)**2 ) / mass(nn)
      enddo
#else
      do n = 1, natom
         RCdot = RCdot + ( fharm(1,n)**2 + fharm(2,n)**2 &
                         + fharm(3,n)**2 ) / mass(n)
      enddo
#endif
      RCdot = sqrt( 2.0d0 * kB * temp0 * RCdot / pi )

   endif

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine evb_umb_primitive


