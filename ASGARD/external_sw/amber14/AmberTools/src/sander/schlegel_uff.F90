! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Schlegel DG-EVB: include UFF VDW interaction to quadratic              |#
! #|  approximation for V_ii                                                 |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#include "../include/dprec.fh"

   module dg_uff

   implicit none

   contains 

! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  UFF combination rules for VDW parameters                               |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   subroutine vdw_mix 

   use constants, only: A_TO_BOHRS, AU_TO_KCAL
   use schlegel,  only: nUFF, nUFF_dcrypt, xvdw, dvdw, zvdw, Avdw, Bvdw &
                      , atomic_numbers, Bohr2Angstrom, Hartree2Kcalmol &
                      , UFF_scale

   implicit none

   ! ...........................................................................

   integer :: i, j, n, ii, jj
   intrinsic :: sqrt, exp

#include "schlegel_uff_parm.h"

!  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
!  :  Geometric mean form for VDW bond length, well depth & shape parameter:   :
!  :                                                                           :
!  :  x_ij = sqrt( x_i x_j )                                                   :
!  :  D_ij = sqrt( D_i D_j )                                                   :
!  :  z_ij = sqrt( z_i z_j )                                                   :
!  :                                                                           :
!  :  Note: (1) xvdw_i & dvdw_i have been converted to Bohrs & Hartrees        :
!  :        (2) diagonal parameters are defined in schlegel_uff_parm.h         :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   do n = 1, nUFF

      i = nUFF_dcrypt(n,1)
      j = nUFF_dcrypt(n,2)

      ii = atomic_numbers(i)
      jj = atomic_numbers(j)

      xvdw(n) = sqrt( UFFvdw(ii)%x * UFFvdw(jj)%x )
      dvdw(n) = sqrt( UFFvdw(ii)%d * UFFvdw(jj)%d )
      zvdw(n) = sqrt( UFFvdw(ii)%z * UFFvdw(jj)%z )

   enddo

!  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
!  :  E_vdw = A exp( -B q ) - C_6 / q^6                                        :
!  :                                                                           :
!  :  A_ij = D_ij ( 6 / (z_ij - 6) ) exp( z_ij )                               :
!  :  B_ij = z_ij / x_ij                                                       :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   Avdw(:) = UFF_scale * dvdw(:) * ( 6.0d0 / ( zvdw(:) - 6.0d0 ) ) * exp( zvdw(:) )
   Bvdw(:) = zvdw(:) / xvdw(:) 

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine vdw_mix


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  UFF VDW initialization: form atom pair decryption file                 |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   subroutine vdw_repulsive_init

   use schlegel, only: nUFF, nUFF_dcrypt, nbond, ibond

   implicit none

   !............................................................................

   integer :: i, j, k, ii, jj, n

   do n = 1, nUFF

      i = nUFF_dcrypt(n,1)
      j = nUFF_dcrypt(n,2)

      do k = 1, nbond

         ii = ibond(k,1)
         jj = ibond(k,2)

         if( i == ii .and. j == jj .or. i == jj .and. j == ii ) &
            nUFF_dcrypt(n,3) = k

      enddo

   enddo

   do n = 1, nUFF
      write(6,'(A,4I8)') 'DG:: nUFF_D = ', n, nUFF_dcrypt(n,1) &
                          , nUFF_dcrypt(n,2), nUFF_dcrypt(n,3)
   enddo

   call vdw_mix

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine vdw_repulsive_init


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  UFF VDW repulsive exponential B term                                   |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   function vdw_repulsive ( q, ncoord )

   use schlegel, only: nUFF, nUFF_dcrypt, Avdw, Bvdw

   implicit none

   integer, intent(in) :: ncoord
   _REAL_ , intent(in) :: q(ncoord)

   !............................................................................

   integer :: k, n
   _REAL_  :: vdw_repulsive
   intrinsic :: exp

!  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
!  :  E_vdw = A exp( -B q ) - C_6 / q^6                                        :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   vdw_repulsive = 0.0d0

   do n = 1, nUFF

      k = nUFF_dcrypt(n,3)
      vdw_repulsive = vdw_repulsive + Avdw(n) * exp( - Bvdw(n) * q(k) )

   enddo

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end function vdw_repulsive


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  UFF VDW repulsive exponential B term (1st derivative)                  |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   function dvdw_repulsive ( q, ncoord )

   use schlegel, only: nUFF, nUFF_dcrypt, Avdw, Bvdw

   implicit none

   integer, intent(in) :: ncoord
   _REAL_ , intent(in) :: q(ncoord)

   !............................................................................

   integer :: k, n
   _REAL_  :: dvdw_repulsive( ncoord )
   intrinsic :: exp

!  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
!  :  (d/q) vdw = - A B exp( -B q ) + 6 C_6 / q^7                              :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   dvdw_repulsive(:) = 0.0d0

   do n = 1, nUFF

      k = nUFF_dcrypt(n,3)

      dvdw_repulsive(k) = dvdw_repulsive(k) - Avdw(n) * Bvdw(n) &
                                            * exp( - Bvdw(n) * q(k) ) 

   enddo

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end function dvdw_repulsive


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  UFF VDW repulsive exponential B term (2nd derivative)                  |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   function ddvdw_repulsive ( q, ncoord )

   use schlegel, only: nUFF, nUFF_dcrypt, Avdw, Bvdw

   implicit none

   integer, intent(in) :: ncoord
   _REAL_ , intent(in) :: q(ncoord)

   !............................................................................

   integer :: k, n
   _REAL_  :: ddvdw_repulsive(ncoord,ncoord)
   intrinsic :: exp

!  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
!  :  (d/dq_i) (d/dq_i) vdw = A B^2 exp( -B q ) - 42 C_6 / q^8                 :
!  :                                                                           :
!  :  q is a set of internals ... so (d/dq_i) (d/dq_j) vdw = 0                 :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   ddvdw_repulsive(:,:) = 0.0d0

   do n = 1, nUFF

      k = nUFF_dcrypt(n,3)
      ddvdw_repulsive(k,k) = ddvdw_repulsive(k,k) + Avdw(n) * Bvdw(n)**2 &
                           * exp( - Bvdw(n) * q(k) ) 

   enddo

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end function ddvdw_repulsive

   end module dg_uff
