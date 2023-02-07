! <compile=optimized>
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  EVB_UMB                                                                |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#  include "../include/dprec.fh"

   subroutine evb_umb ( f, q, mass, natom, istart3, iend3 )

   use constants, only: kB
   use evb_parm,  only: k_umb, r0_umb, evb_dyn, nbias, dbonds_RC, bond_RC &
                      , out_RCdot, inc_dbonds_RC, inc_bond_RC
   use evb_data,  only: evb_frc, evb_bias, RCdot
#ifdef DEBUG_EVB
   use evb_check, only: full_evb_debug, dbonds_debug, bond_debug 
#endif /* DEBUG_EVB */
#ifdef LES
   use evb_pimd,  only: bead_dcrypt, natomCL
   use miller,    only: gradRC, gradRC_norm, i_qi, div_ndx
#endif

   implicit none

#include "../include/md.h"

   integer, intent(   in) :: natom, istart3, iend3
   _REAL_ , intent(   in) :: q(natom*3)
   _REAL_ , intent(   in) :: mass(natom)
   _REAL_ , intent(inout) :: f(natom*3)

   !  ..........................................................................

   integer :: m, n, nn
   integer :: idx, jdx, kdx
   _REAL_  :: fharm(natom*3), dr(3), fr(3), rij, rkj, rij_inv, rkj_inv
   _REAL_  :: evb_fbias(natom*3,nbias), RC, pi
#ifdef LES
   integer :: mm_dx, mm, m_dx
   _REAL_  :: norm_gradRC
#endif

   select case( trim( adjustl( evb_dyn) ) )

!  +---------------------------------------------------------------------------+
!  |  Difference of 2 bonds RC harmonic umbrella sampling                      |
!  |                                                                           |
!  |  /\ = r_ij - r_kj                                                         |
!  |                                                                           |
!  |  V'(/\) = V_0 + 0.5 k_evb * ( /\ - /\_0 )^2                               |
!  |                                                                           |
!  |  dV' / dR = ( dV' / d/\ ) * ( d/\ / dR )                                  |
!  |           = dV_0 / dR + k_evb * ( /\ - /\_0 ) * d/dR ( r_ij - r_kj )      |
!  +---------------------------------------------------------------------------+

!  +---------------------------------------------------------------------------+
!  |  12312008: For outputting RC when doing groundstate dynamics              |
!  +---------------------------------------------------------------------------+

      case( "groundstate" )

         do n = 1, nbias

            if( inc_dbonds_RC ) then

               idx = ( dbonds_RC(n)%iatom - 1 ) * 3
               jdx = ( dbonds_RC(n)%jatom - 1 ) * 3
               kdx = ( dbonds_RC(n)%katom - 1 ) * 3

               do nn = 1, 3
                  dr(nn) = q(idx+nn) - q(jdx+nn)
               enddo
               rij = sqrt( dr(1)**2 + dr(2)**2 + dr(3)**2 )

               do nn = 1, 3
                  dr(nn)  = q(kdx+nn) - q(jdx+nn)
               enddo
               rkj = sqrt( dr(1)**2 + dr(2)**2 + dr(3)**2 )

               RC = rij - rkj
               evb_bias%RC(n) = RC

            endif

            if( inc_bond_RC ) then

               idx = ( bond_RC(n)%iatom - 1 ) * 3
               jdx = ( bond_RC(n)%jatom - 1 ) * 3

               do nn = 1, 3
                  dr(nn) = q(idx+nn) - q(jdx+nn)
               enddo
               rij = sqrt( dr(1)**2 + dr(2)**2 + dr(3)**2 )

               RC = rij
               evb_bias%RC(n) = RC

            endif

         enddo

!  |  12312008: For outputting RC when doing groundstate dynamics              |

      case( "dbonds_umb", "qi_dbonds_pmf", "qi_dbonds_dyn" )

!        fharm(:) = 0.0d0

         do n = 1, nbias

            fharm(:) = 0.0d0

            idx = ( dbonds_RC(n)%iatom - 1 ) * 3
            jdx = ( dbonds_RC(n)%jatom - 1 ) * 3
            kdx = ( dbonds_RC(n)%katom - 1 ) * 3

! write(6,*) 'nbias cmd >>> ', n, idx, jdx, kdx
! write(6,*) 'q(i) = ', q(idx+1), q(idx+2), q(idx+3)
! write(6,*) 'q(j) = ', q(jdx+1), q(jdx+2), q(jdx+3)
! write(6,*) 'q(k) = ', q(kdx+1), q(kdx+2), q(kdx+3)

            do nn = 1, 3
               dr(nn) = q(idx+nn) - q(jdx+nn)
            enddo 

            rij = sqrt( dr(1)**2 + dr(2)**2 + dr(3)**2 )
            rij_inv = 1.0d0 / rij

            do nn = 1, 3
               fr(nn) = dr(nn) * rij_inv
               fharm(idx+nn) = fharm(idx+nn) - fr(nn)
               fharm(jdx+nn) = fharm(jdx+nn) + fr(nn)
            enddo 

            do nn = 1, 3 
               dr(nn)  = q(kdx+nn) - q(jdx+nn)
            enddo

            rkj = sqrt( dr(1)**2 + dr(2)**2 + dr(3)**2 )
            rkj_inv = 1.0d0 / rkj

            do nn = 1, 3
               fr(nn) = dr(nn) * rkj_inv
               fharm(kdx+nn) = fharm(kdx+nn) + fr(nn)
               fharm(jdx+nn) = fharm(jdx+nn) - fr(nn)
            enddo

#ifdef LES
            if( i_qi > 0 ) then
               do m = 1, natomCL
                  mm = bead_dcrypt( m, div_ndx(n) )
                   m_dx = (  m - 1 ) * 3
                  mm_dx = ( mm - 1 ) * 3
                  gradRC(m_dx+1,n) = fharm(mm_dx+1)
                  gradRC(m_dx+2,n) = fharm(mm_dx+2)
                  gradRC(m_dx+3,n) = fharm(mm_dx+3)
               enddo
            endif
#endif
            RC = rij - rkj

! write(6,*) 'cmd rij = ', rij
! write(6,*) 'cmd rkj = ', rkj
! write(6,*) 'cmd RC  = ', RC 

            evb_bias%RC(n) = RC
            evb_bias%nrg_bias(n) = 0.50d0 * k_umb(n) * ( RC - r0_umb(n) )**2
            evb_fbias(:,n) = k_umb(n) * ( RC - r0_umb(n) ) * fharm(:)
#ifdef DEBUG_EVB
            if( dbonds_debug .or. full_evb_debug ) &
               call dbonds_anal2num ( q, evb_fbias(:,n), natom*3 )
#endif

         enddo

         do n = 1, nbias
            do m = istart3, iend3
               f(m) = f(m) + evb_fbias(m,n)
            enddo
            evb_frc%evb_nrg  = evb_frc%evb_nrg + evb_bias%nrg_bias(n)
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

!        fharm(:) = 0.0d0

         do n = 1, nbias

            fharm(:) = 0.0d0

            idx = ( bond_RC(n)%iatom - 1 ) * 3
            jdx = ( bond_RC(n)%jatom - 1 ) * 3

            do nn = 1, 3
               dr(nn) = q(idx+nn) - q(jdx+nn)
            enddo

            rij = sqrt( dr(1)**2 + dr(2)**2 + dr(3)**2 )
            rij_inv = 1.0d0 / rij

            do nn = 1, 3
               fr(nn) = dr(nn) * rij_inv
               fharm(idx+nn) = fharm(idx+nn) - fr(nn)
               fharm(jdx+nn) = fharm(jdx+nn) + fr(nn)
            enddo

#ifdef LES
            if( i_qi > 0 ) then
               do m = 1, natomCL
                  mm = bead_dcrypt( m, div_ndx(n) )
                   m_dx = (  m - 1 ) * 3
                  mm_dx = ( mm - 1 ) * 3
                  gradRC(m_dx+1,n) = fharm(mm_dx+1)
                  gradRC(m_dx+2,n) = fharm(mm_dx+2)
                  gradRC(m_dx+3,n) = fharm(mm_dx+3)
               enddo
            endif
#endif

            RC = rij 

            evb_bias%RC(n) = RC
            evb_bias%nrg_bias(n) = 0.50d0 * k_umb(n) * ( RC - r0_umb(n) )**2
            evb_fbias(:,n) = k_umb(n) * ( RC - r0_umb(n) ) * fharm(:)
#ifdef DEBUG_EVB
            if( bond_debug .or. full_evb_debug ) &
               call bond_anal2num ( q, evb_fbias(:,n), natom*3 )
#endif

         enddo

         do n = 1, nbias
            do m = istart3, iend3
               f(m) = f(m) + evb_fbias(m,n)
            enddo
            evb_frc%evb_nrg  = evb_frc%evb_nrg + evb_bias%nrg_bias(n)
         enddo

   end select

! Output gradRC_norm for QI

#ifdef LES
   if( i_qi > 0 ) then
      do n = 1, nbias
         norm_gradRC = 0.0d0
         do m = 1, natomCL
            m_dx = (  m - 1 ) * 3
            mm = bead_dcrypt( m, div_ndx(n) )
            norm_gradRC = norm_gradRC + ( gradRC(m_dx+1,n)**2 + gradRC(m_dx+2,n)**2 &
                                        + gradRC(m_dx+3,n)**2 ) / mass(mm)

!write(6,'(A,I8, 4F12.8)') '>>>', m, gradRC(m_dx+1,n), gradRC(m_dx+2,n), gradRC(m_dx+3,n), mass(mm)

         enddo
         gradRC_norm(n) = sqrt( norm_gradRC )

!write(6,*) 'gradRC_norm = ', n, gradRC_norm(n)

      enddo
   endif
#endif

   if( out_RCdot ) then

      pi = acos( -1.0d0 )
      RCdot = 0.0d0
      do n = 1, natom
         idx = ( n - 1 ) * 3
         RCdot = RCdot + ( fharm(idx+1)**2 + fharm(idx+2)**2 &
                         + fharm(idx+3)**2 ) / mass(n)

! write(6,*) 'RCdot >>> ', n, mass(n)

      enddo

      RCdot = sqrt( 2.0d0 * kB * temp0 * RCdot / pi )

   endif

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine evb_umb


