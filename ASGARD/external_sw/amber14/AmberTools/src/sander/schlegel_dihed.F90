! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Schlegel DG-EVB: include Amber torsion interaction to quadratic        |#
! #|  approximation for V_ii                                                 |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#include "../include/assert.fh"
#include "../include/dprec.fh"


   module dg_dihed

   implicit none

   contains 

! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Amber dihedral initialization: form atom quartet decryption file       |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   subroutine amber_dihed_init

   use constants, only: AU_TO_KCAL 
   use schlegel,  only: nAdihed, nAdihed_dcrypt, ndihed, idihed, nbond, nangle, V

   implicit none

   !............................................................................

   integer :: i, j, k, l, m, n, ii, jj, kk, ll, offset

!  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
!  :  For internal coordinate derivatives                                      :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   offset = nbond + nangle

   do n = 1, nAdihed

      i = nAdihed_dcrypt(n,1)
      j = nAdihed_dcrypt(n,2)
      k = nAdihed_dcrypt(n,3)
      l = nAdihed_dcrypt(n,4)

      do m = 1, ndihed

         ii = idihed(m,1)
         jj = idihed(m,2)
         kk = idihed(m,3)
         ll = idihed(m,4)

         if( i == ii .and. j == jj .and. k == kk .and. l == ll .or. &
             i == ll .and. j == jj .and. k == kk .and. l == ii        ) &
            nAdihed_dcrypt(n,5) = m + offset

      enddo

   enddo

   do n = 1, nAdihed
      write(6,'(A,6I8)') "DG:: nAdihed_D = ", n, nAdihed_dcrypt(n,1) &
                           , nAdihed_dcrypt(n,2), nAdihed_dcrypt(n,3) &
                           , nAdihed_dcrypt(n,4), nAdihed_dcrypt(n,5)
   enddo

!  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
!  :  Convert dihedral force constant from kcal/mol to hartree                 :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   V(:) = V(:) / AU_TO_KCAL


!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine amber_dihed_init


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Amber dihedral term                                                    |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   function amber_dihed ( q, ncoord )

   use schlegel, only: nAdihed, nAdihed_dcrypt, V, period, phase

   implicit none

   integer, intent(in) :: ncoord
   _REAL_ , intent(in) :: q(ncoord)

   !............................................................................

   integer :: k, n
   _REAL_  :: amber_dihed

!  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
!  :  E_dihed = V_n / 2 * [ 1 + cos( n phi - y ) ]                             :
!  :                                                                           :
!  :  Note: V(n) should already be halved                                      :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   amber_dihed = 0.0d0

   do n = 1, nAdihed

      k = nAdihed_dcrypt(n,5)
      amber_dihed = amber_dihed + V(n) &
                  * ( 1.0d0 + cos( period(n) * q(k) - phase(n) ) )

   enddo

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end function amber_dihed


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Amber dihedral term (1st derivative)                                   |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   function damber_dihed ( q, ncoord )

   use schlegel, only: nAdihed, nAdihed_dcrypt, V, period, phase

   implicit none

   integer, intent(in) :: ncoord
   _REAL_ , intent(in) :: q(ncoord)

   !............................................................................

   integer :: k, n
   _REAL_  :: damber_dihed( ncoord )

!  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
!  :  (d/q) E_dihed = - V_n / 2 * sin( n phi - y ) * n                         :
!  :                                                                           :
!  :  Note: V(n) should already be halved                                      :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   damber_dihed(:) = 0.0d0

   do n = 1, nAdihed

      k = nAdihed_dcrypt(n,5)

      damber_dihed(k) = damber_dihed(k) - V(n) &
                      * sin( period(n) * q(k) - phase(n) ) * period(n)

   enddo

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end function damber_dihed


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Amber dihedral term (2nd derivative)                                   |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   function ddamber_dihed ( q, ncoord )

   use schlegel, only: nAdihed, nAdihed_dcrypt, V, period, phase

   implicit none

   integer, intent(in) :: ncoord
   _REAL_ , intent(in) :: q(ncoord)

   !............................................................................

   integer :: k, n
   _REAL_  :: ddamber_dihed(ncoord,ncoord)

!  ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
!  :  (d/q_i) (d/dq_i) E_dihed = - V_n / 2 * cos( n phi - y ) * n^2            :
!  :                                                                           :
!  :  Note: V(n) should already be halved                                      :
!  :  q is a set of internals ... so (d/dq_i) (d/dq_j) E_dihed = 0             :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   ddamber_dihed(:,:) = 0.0d0

   do n = 1, nAdihed

      k = nAdihed_dcrypt(n,5)
      ddamber_dihed(k,k) = ddamber_dihed(k,k) - V(n) &
                      * cos( period(n) * q(k) - phase(n) ) * period(n)**2

   enddo

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end function ddamber_dihed

   end module dg_dihed
