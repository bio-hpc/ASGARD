! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Schlegel DG-EVB: include Amber angle interaction to quadratic          |#
! #|  approximation for V_ii                                                 |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#include "../include/assert.fh"
#include "../include/dprec.fh"

   module dg_angle

   implicit none

   contains 

! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Amber angle initialization: form atom triplet decryption file          |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   subroutine amber_angle_init

   use constants, only: AU_TO_KCAL 
   use schlegel,  only: nAangle, nAangle_dcrypt, iangle, nbond, nangle, ndihed, K_angle

   implicit none

   !............................................................................

   integer :: i, j, k, m, n, ii, jj, kk, offset

!  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
!  :  For internal coordinate derivatives                                      :
!  :''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''':
!  |                      3                                                    |
!  |                    /   \                                                  |
!  |                  1       2                                                |
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   write(6,*) 'nbond = ', nbond
   write(6,*) 'nangle = ', nangle
   write(6,*) 'ndihed = ', ndihed
   write(6,*) 'K_angle(:) = ', K_angle(:)

   do m = 1, nangle
      write(6,*) m, iangle(m,1), iangle(m,3), iangle(m,2)
   enddo
 
   offset = nbond 

   do n = 1, nAangle

      i = nAangle_dcrypt(n,1)
      j = nAangle_dcrypt(n,2)
      k = nAangle_dcrypt(n,3)

      do m = 1, nangle

         ii = iangle(m,1)
         jj = iangle(m,2)
         kk = iangle(m,3)

  write(6,*) i, j, k
  write(6,*) ii, kk, jj
  write(6,*)

         if( i == ii .and. j == kk .and. k == jj .or. &
             i == jj .and. j == kk .and. k == ii        ) &
            nAangle_dcrypt(n,4) = m + offset

      enddo

   enddo

   do n = 1, nAangle
      write(6,'(A,5I8)') "DG:: nAangle_D = ", n, nAangle_dcrypt(n,1) &
             , nAangle_dcrypt(n,2), nAangle_dcrypt(n,3), nAangle_dcrypt(n,4) 
   enddo

!  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
!  :  Convert angle force constant from kcal/mol to hartree                    :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   K_angle(:) = K_angle(:) / AU_TO_KCAL


!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine amber_angle_init


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Amber angle term                                                       |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   function amber_angle ( q, ncoord )

   use schlegel, only: nAangle, nAangle_dcrypt, K_angle, theta

   implicit none

   integer, intent(in) :: ncoord
   _REAL_ , intent(in) :: q(ncoord)

   !............................................................................

   integer :: k, n
   _REAL_  :: amber_angle

!  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
!  :  E_angle = K_theta * ( theta - theta_0 )**2                               :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   amber_angle = 0.0d0

   do n = 1, nAangle

      k = nAangle_dcrypt(n,4)

      amber_angle = amber_angle + K_angle(n) * ( q(k) - theta(n) )**2

   enddo

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end function amber_angle


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Amber angle term (1st derivative)                                      |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   function damber_angle ( q, ncoord )

   use schlegel, only: nAangle, nAangle_dcrypt, K_angle, theta

   implicit none

   integer, intent(in) :: ncoord
   _REAL_ , intent(in) :: q(ncoord)

   !............................................................................

   integer :: k, n
   _REAL_  :: damber_angle( ncoord )

!  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
!  :  (d/q) E_angle = 2.0 * K_theta * ( q(k) - theta(n) )                      :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   damber_angle(:) = 0.0d0

   do n = 1, nAangle

      k = nAangle_dcrypt(n,4)

      damber_angle(k) = damber_angle(k) + 2.0d0 * K_angle(n) * (q(k) - theta(n) )

   enddo

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end function damber_angle


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Amber angle term (2nd derivative)                                      |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   function ddamber_angle ( ncoord )

   use schlegel, only: nAangle, nAangle_dcrypt, K_angle

   implicit none

   integer, intent(in) :: ncoord

   !............................................................................

   integer :: k, n
   _REAL_  :: ddamber_angle(ncoord,ncoord)

!  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
!  :  (d/q_i) (d/dq_i) E_angle =  2.0 * K_theta                                :
!  :                                                                           :
!  :  q is a set of internals ... so (d/dq_i) (d/dq_j) E_angle = 0             :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   ddamber_angle(:,:) = 0.0d0

   do n = 1, nAangle

      k = nAangle_dcrypt(n,4)
      ddamber_angle(k,k) = ddamber_angle(k,k) + 2.0d0 * K_angle(n)

   enddo

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end function ddamber_angle

   end module dg_angle
