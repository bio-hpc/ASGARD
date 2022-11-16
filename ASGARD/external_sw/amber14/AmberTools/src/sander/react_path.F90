#include "../include/dprec.fh"

! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  REACT_TRJ: sample a reaction path by going downhill            |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   subroutine react_path ( xq, ndim, RC_ndx )  

   use evb_parm,  only: bias_ndx, dbonds_RC, bond_RC 
   use evb_data,  only: evb_Hmat, evb_bias 
   use wigner,    only: rflux, RC_type, nbias_ndx
   use file_io_dat

   implicit none

   integer, intent(in) :: ndim, RC_ndx
   _REAL_ , intent(in) :: xq(ndim)

#  include "parallel.h"
#  include "../include/md.h"

   !  +---------------------------------------------------------------+

   integer :: n, nn, ni, nf, idx, jdx, kdx
   _REAL_  :: dr(3), rij, rkj


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

         rflux%RC_trj(RC_ndx) &
            = evb_Hmat%evb_mat(ni,ni) - evb_Hmat%evb_mat(nf,nf)

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

         n = nbias_ndx

         idx = ( dbonds_RC(n)%iatom - 1 ) * 3
         jdx = ( dbonds_RC(n)%jatom - 1 ) * 3
         kdx = ( dbonds_RC(n)%katom - 1 ) * 3

         do nn = 1, 3
            dr(nn) = xq(idx+nn) - xq(jdx+nn)
         enddo 

         rij = sqrt( dr(1)**2 + dr(2)**2 + dr(3)**2 )

         do nn = 1, 3 
            dr(nn)  = xq(kdx+nn) - xq(jdx+nn)
         enddo

         rkj = sqrt( dr(1)**2 + dr(2)**2 + dr(3)**2 )

         rflux%RC_trj(RC_ndx) = rij - rkj

!  +---------------------------------------------------------------+
!  |  Bond RC harmonic umbrella sampling                           |
!  |                                                               |
!  |  /\ = r_ij                                                    |
!  |                                                               |
!  |  - d/dR_i /\ = - ( R_i - R_j ) / rij                          |
!  |  - d/dR_j /\ = + ( R_i - R_j ) / rij                          |
!  +---------------------------------------------------------------+

      case( "bond" )

         n = nbias_ndx

         idx = ( bond_RC(n)%iatom - 1 ) * 3
         jdx = ( bond_RC(n)%jatom - 1 ) * 3

         do nn = 1, 3
            dr(nn) = xq(idx+nn) - xq(jdx+nn)
         enddo

         rij = sqrt( dr(1)**2 + dr(2)**2 + dr(3)**2 )

         rflux%RC_trj(RC_ndx) = rij 

   end select

   evb_bias%RC(1) = rflux%RC_trj(RC_ndx)

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine react_path

