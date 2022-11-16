! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Form the Schlegel-Sonnenberg D matrix and F vector and perform SVD     |#
! #|  for determining the B vector ( D B = F )                               |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#  include "../include/assert.fh"
#  include "../include/dprec.fh"

   subroutine schlegel_full ( f, d )

   use schlegel, only: ndg, ncoord, dgdim, b, dg_weight, gbasis_weight &
                     , svd_cond, full_solve

   implicit none

   _REAL_   , intent(in) :: f(dgdim*ndg), d(dgdim*ndg,dgdim*ndg)

   !  ..........................................................................

   integer :: k, n, ndx, info, lwork, rank, alloc_error, dealloc_error
   integer, allocatable :: iwork(:), jpvt(:)
   _REAL_  :: wnrg, wgrad, whess
   intrinsic :: matmul, transpose
   _REAL_,  allocatable :: w(:), wd(:,:), a(:,:), a_inv(:,:)
   _REAL_,  allocatable :: s(:), s_inv(:,:), u(:,:), vt(:,:) 
   _REAL_,  allocatable :: work(:)
   intrinsic :: trim, adjustl

!  .............................................................................
!  :  Allocate temporary arrays                                                :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   allocate( iwork( dgdim*ndg*128 ), jpvt(dgdim*ndg), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( w(dgdim*ndg), wd(dgdim*ndg,dgdim*ndg), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( a(dgdim*ndg,dgdim*ndg), a_inv(dgdim*ndg,dgdim*ndg), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( s(dgdim*ndg), s_inv(dgdim*ndg,dgdim*ndg), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( u(dgdim*ndg,dgdim*ndg), vt(dgdim*ndg,dgdim*ndg), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( work(dgdim*ndg*2*128), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

!  +---------------------------------------------------------------------------+
!  |  Form diagonal weighting matrix                                           |
!  +---------------------------------------------------------------------------+

   wnrg  = gbasis_weight(1)
   wgrad = gbasis_weight(2)
   whess = gbasis_weight(3)

   do n = 1, ndg
       ndx = dgdim * ( n - 1 ) + 1
       w(ndx) = dg_weight(n) * wnrg
       w(ndx+1:ndx+ncoord) = dg_weight(n) * wgrad
       w(ndx+ncoord+1:ndx+ncoord+ncoord*ncoord) = dg_weight(n) * whess
   enddo 

!  +--------------------------------------------------------------------------+
!  |  Solve for the B vector                                                  |
!  :..........................................................................:
!  |        D B =       F                                                     |
!  |  D^T W D B = D^T W F                                                     |
!  |          B = [D^T W D]^-1 D^T W F                                        |
!  |                                                                          |
!  |  solve for [D^T W D]^-1 = def A^-1 using singular value decomposition,   |
!  |     where A = U S V^T ==> A^-1 = V S^-1 U^T                              |
!  +--------------------------------------------------------------------------+

   do n = 1, dgdim * ndg
     wd(:,n) = w(:) * d(:,n)
   enddo

   a = matmul( transpose( d ), wd )

   select case( trim( adjustl( full_solve ) ) )

      case( "berny_solve" )

         call dgesvd( 'A', 'A', dgdim*ndg, dgdim*ndg, a, dgdim*ndg, s, u &
                    , dgdim*ndg, vt, dgdim*ndg, work, dgdim*ndg*10, info ) 

         if( info /= 0 ) then 
            write(6,'(A)') "ERROR encountered during SVD in schlegel_dg"
            write(6,'(A,I8)') "info = ", info
            call mexit(6, 1)
         endif 

         s_inv(:,:) = 0.0d0
         k = 0
         do n = 1, dgdim * ndg
            if( s(n) > svd_cond ) s_inv(n,n) = 1.0d0 / s(n)
            if( s_inv(n,n) /= 0.0d0 ) k = k + 1
         enddo

         write(6,100) ">>> Inverted ", k, "singular values in schlegel_full"

         a_inv = matmul( transpose( vt ), matmul( s_inv, transpose( u ) ) )
         b = matmul( a_inv, matmul( transpose( d ), w(:) * f(:) ) ) 

      case( "lapack_dgelsd" )

         b(:) = matmul( transpose( d ), w(:) * f(:) )

         lwork = dgdim * ndg * 128 * 2
         call dgelsd ( dgdim*ndg, dgdim*ndg, 1, a, dgdim*ndg, b, dgdim*ndg &
                     , s, svd_cond, rank, work, lwork, iwork, info )

         write(6,'(A,I8)') "rank of matrix a = ", rank
         if( info /= 0 ) then
            write(6,'(A)') "ERROR encountered during dgelsd in schlegel_full"
            write(6,'(A,I8)') 'info = ', info
            write(6,'(A,I8)') 'rank of matrix a = ', rank
            call mexit(6, 1)
         endif

      case( "lapack_dgelsy" )

         b(:) = matmul( transpose( d ), w(:) * f(:) )

         lwork = dgdim * ndg * 128 * 2

         jpvt(:) = 0
         call dgelsy ( dgdim*ndg, dgdim*ndg, 1, a, dgdim*ndg, b, dgdim*ndg &
                     , jpvt, svd_cond, rank, work, lwork, info )

         write(6,'(A,I8)') "rank of matrix a = ", rank
         if( info /= 0 ) then
            write(6,'(A)') "ERROR encountered during dgelsy in schlegel_full"
            write(6,'(A,I8)') "info = ", info
            call mexit(6, 1)
         endif

   end select

!  .............................................................................
!  :  Deallocate temporary arrays                                              :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   deallocate( iwork, jpvt, stat = dealloc_error )
   REQUIRE( dealloc_error == 0 )
   deallocate( w, wd, a, a_inv, stat = dealloc_error )
   REQUIRE( dealloc_error == 0 ) 
   deallocate( s, s_inv, u, vt, work, stat = dealloc_error )
   REQUIRE( dealloc_error == 0 )

 100  format( A, I8, A )

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine schlegel_full


