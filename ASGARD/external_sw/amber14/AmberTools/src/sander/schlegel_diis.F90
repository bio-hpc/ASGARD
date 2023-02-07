! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Form the Schlegel-Sonnenberg D matrix and F vector and perform DIIS    |#
! #|  for determining the B vector: ( D B = F )                              |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#  include "../include/assert.fh"
#  include "../include/dprec.fh"

   subroutine schlegel_diis ( f, b, ndiis_max, ndiis_dim )

   use schlegel,  only: ndg, dgdim, diis_tol 

   implicit none

   integer, intent( in) :: ndiis_max, ndiis_dim
   _REAL_ , intent( in) :: f(dgdim*ndg)
   _REAL_ , intent(out) :: b(dgdim*ndg)

   !  ..........................................................................

   integer :: i, j, k, l, kp1, nvec, info, alloc_error, dealloc_error
   integer :: ipiv(ndiis_dim)
   intrinsic :: mod, min, max
   _REAL_  :: rms
   intrinsic :: dot_product, transpose, matmul, sqrt
   _REAL_, allocatable :: Basis(:,:), DBasis(:,:), TB(:,:), TDB(:,:)
   _REAL_, allocatable :: DIISmat(:,:), DIISvec(:), DBmat(:,:), coeff(:)
   _REAL_, allocatable :: Db(:), shit(:), work(:)

!  .............................................................................
!  :  Allocate temporary arrays                                                :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   allocate( Basis(dgdim*ndg,ndiis_dim), stat = alloc_error )
   REQUIRE( alloc_error == 0 ) 
   allocate( DBasis(dgdim*ndg,ndiis_dim), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( TB(ndiis_dim,dgdim*ndg), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( TDB(ndiis_dim,dgdim*ndg), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( DIISmat(ndiis_dim,ndiis_dim), stat = alloc_error )
   REQUIRE( alloc_error == 0 ) 
   allocate( DIISvec(ndiis_dim), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( DBmat(ndiis_dim,ndiis_dim), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( coeff(ndiis_dim), Db(dgdim*ndg), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( shit(dgdim*ndg), work(ndiis_dim*2*128), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   
!  +---------------------------------------------------------------------------+
!  |  Solve the linear equations [D] [B] = [F] using GMRES                     |
!  :...........................................................................:
!  |  D B = F                                                                  |
!  |                                                                           |
!  |  B = sum_i c_i b_i                                                        |
!  |                                                                           |
!  |  R^2 = 0.5 ( D B - F )^2 = 0.5 ( D sum_i c_i b_i - F )^2                  |
!  |      = 0.50 [ ( D sum_i c_i b_i )^T ( D sum_j c_j b_j )                   |
!  |      - 2 ( D sum_i c_i b_i )^T F + F^2 ]                                  |
!  |                                                                           |
!  |  MINIMIZE R^2 with respect to c_i:                                        |
!  |  d/dc_i R^2 = ( D sum_i b_i )^T ( D sum_j c_j b_j )                       |
!  |             - ( D sum_i b_i )^T F = 0                                     |
!  |                                                                           |
!  |  DIISmat_ij = ( D b_i )^T ( D b_j )                                       |
!  |  DIISvec_i  = ( D b_i )^T F                                               |
!  |                                                                           |
!  |  B_i+1 = B_i + F - D B_i                                                  |
!  +---------------------------------------------------------------------------+

   Basis(:,1) = f(:)

   do i = 1, ndiis_max 

      k    = mod( i - 1, ndiis_dim ) + 1
      kp1  = mod( i, ndiis_dim ) + 1
      nvec = min( i, ndiis_dim )

      call schlegel_DB ( Basis(:,k), DBasis(:,k) )

      DIISvec(k) = dot_product( f, DBasis(:,k) )

      do j = max( 1, i - ndiis_dim + 1 ), i

         l = mod( j - 1, ndiis_dim ) + 1
         DIISmat(k,l) = dot_product( DBasis(:,k), DBasis(:,l) )
         DIISmat(l,k) = DIISmat(k,l)

      enddo

      coeff(  :) = DIISvec(  :)
      DBmat(:,:) = DIISmat(:,:)

!  +---------------------------------------------------------------------------+
!  |  Solve the linear system of equations D B = F, where D is symmetric       |
!  +---------------------------------------------------------------------------+

      call dsysv ( 'L', nvec, 1, DBmat(1:nvec,1:nvec), nvec, ipiv(1:nvec) &
                 , coeff(1:nvec), nvec, work, nvec*128, info )

       TB(:,:) = transpose(  Basis(:,:) )
      TDB(:,:) = transpose( DBasis(:,:) )

      Basis(:,kp1) = matmul( coeff(1:nvec),  TB(1:nvec,1:dgdim*ndg) ) &
                   - matmul( coeff(1:nvec), TDB(1:nvec,1:dgdim*ndg) ) &
                   + f(:)

      shit(:) = Basis(:,kp1) - Basis(:,k)
      rms = sqrt ( dot_product( shit, shit ) )

      write(6,'(A, I8, 1PE16.8)') "DIIS iter, rms = ", i, rms

      if( rms <= diis_tol ) exit
 
   enddo

   if( rms <= diis_tol ) then
      write(6,100) "DIIS converged to ", diis_tol, " in ", i, " iterations"
      b(:) = Basis(:,kp1) 
      call schlegel_DB ( b(:), Db(:) )
      rms = sqrt ( dot_product( Db(:) - f(:), Db(:) - f(:) ) )
      write(6,'(A,1PE16.8)') ">>> FINAL DEVIATION = ", rms
   else
      write(6,100) "DIIS failed to converged to ", diis_tol, " after " &
                 , i, " iterations"
   endif

!  .............................................................................
!  :  Deallocate temporary arrays                                              :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   deallocate( Basis, DBasis, TB, TDB, stat = dealloc_error )
   REQUIRE( dealloc_error == 0 )
   deallocate( DIISmat, DIISvec, DBmat, coeff, stat = dealloc_error )
   REQUIRE( dealloc_error == 0 )
   deallocate( Db, shit, work, stat = dealloc_error )
   REQUIRE( dealloc_error == 0 )


100 format( A, E12.5, A, I8, A )

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine schlegel_diis


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Form the F vector :: D B = F                                           |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   subroutine schlegel_DB ( b, f )

   use schlegel, only: ndg, xdat_xdg, rdgdim 

   implicit none

   _REAL_ , intent( in) :: b(rdgdim*ndg)
   _REAL_ , intent(out) :: f(rdgdim,ndg)

   !............................................................................

   integer :: n 
   intrinsic :: reshape

!  +---------------------------------------------------------------------------+
!  |  Form F vector                                                            |
!  +---------------------------------------------------------------------------+

   do n = 1, ndg
      call ddv12sqDGR ( xdat_xdg(n)%q(:), b, f(:,n) )
   enddo

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine schlegel_DB


