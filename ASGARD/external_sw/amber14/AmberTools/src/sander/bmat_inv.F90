! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Compute the inverse of the B matrix                                    |#
! #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#+#+#+#+#

#include "../include/dprec.fh"

   subroutine bmat_inv ( bmat, binv, natm, ncoord )

   use schlegel, only: svd_cond

   implicit none

   integer, intent( in) :: natm, ncoord 
   _REAL_ , intent( in) :: bmat(ncoord,natm*3) 
   _REAL_ , intent(out) :: binv(natm*3,ncoord)

   !  ..........................................................................

   integer :: n, info 
   _REAL_  :: s(natm*3), s_inv(natm*3,ncoord), u(ncoord,ncoord) &
            , vt(natm*3,natm*3), work(ncoord*128)
   intrinsic :: matmul, transpose
   _REAL_ :: bmat_tmp(ncoord,natm*3)

!  +---------------------------------------------------------------------------+
!  |  Solve for B^-1 using singular value decomposition                        |
!  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  |  For a given B, B^-1 = V S^-1 U^T                                         |
!  +---------------------------------------------------------------------------+

! necessary because ifort permits dgevd to change bmat despite the intent(in)
   bmat_tmp(:,:) = bmat(:,:)  
! necessary because ifort permits dgevd to change bmat despite the intent(in)
   call dgesvd( 'A', 'A', ncoord, natm*3, bmat_tmp, ncoord, s, u &
              , ncoord, vt, natm*3, work, ncoord*128, info )

   if( info /= 0 ) then
      write(6,'(A)') "ERROR encountered during SVD for solving B^-1"
      write(6,'(A,I8)') 'info = ', info
      call mexit(6, 1)
   endif

   s_inv(:,:) = 0.0d0
   do n = 1, natm*3
      if( s(n) > svd_cond ) s_inv(n,n) = 1.0d0 / s(n)
   enddo

   binv = matmul( transpose( vt ), matmul( s_inv, transpose( u ) ) )

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine bmat_inv

