! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  GMRES for iterative solution of the B vector :: D B = F                |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#  include "../include/assert.fh"
#  include "../include/dprec.fh"

   subroutine schlegel_gmres ( b0, f, ndim, ndiis_dim )

   use schlegel,  only: b, diis_tol, alpha 
   use file_io_dat
#ifdef DEBUG_EVB
   use evb_check, only: schlegel_debug
#endif

   implicit none

   integer, intent(in) :: ndim, ndiis_dim
   _REAL_ , intent(in) :: b0(ndim), f(ndim)

   !  ..........................................................................

   integer :: j, k, ipiv(ndiis_dim), info
#ifdef DEBUG_EVB
   integer :: io
#endif
   intrinsic :: min
   _REAL_  :: e(ndim,ndiis_dim), De(ndim), r0(ndim), r0_nrm
   _REAL_  :: uhess(ndiis_dim,ndiis_dim), qr(ndiis_dim,ndiis_dim) &
            , HTH(ndiis_dim,ndiis_dim), coeff(ndiis_dim)
   _REAL_  :: b_gmres(ndim,ndiis_dim), resid(ndiis_dim), v(ndiis_dim)
   _REAL_  :: tau(ndiis_dim), work(ndiis_dim*2*128) 
   intrinsic :: sqrt, dot_product, transpose, matmul, abs

!  +---------------------------------------------------------------------------+
!  |  e_1 = ( F - D B_0 ) / || ( F - D B_0 ) ||                                |
!  |  h_ik = < e_i | D e_k >      (i=1,k; k=1,m)                               |
!  |  e_(k+1) = ( D e_k - sum(i->k) h_ik e_i ) / h_(k+1,k)                     |
!  |  h_(k+1,k) = || D e_k - sum(i->k) h_ik e_i ||                             |
!  +---------------------------------------------------------------------------+

   call schlegel_DB ( b0(:), De(:) )

#ifdef DEBUG_EVB
   if( schlegel_debug ) then
      io = schlegel_unit
      open ( io, file = "DB_vector.debug" )
      write( io, '( A )' ) "[Data Size]"
      write( io, '( I20 )' ) ndim
      write( io, '( A )' ) "[DB Vector]"
      write( io, '(5(1PE16.8))' ) ( De(j), j = 1, ndim )
      close( io )
   endif
#endif

   r0(:) = f(:) - De(:)
   r0_nrm = sqrt ( dot_product( r0, r0 ) )

   b_gmres(:,1) = b0(:)
  
   v(:) = 0.0d0
   v(1) = 1.0d0

   write(6,'(A,2X,I8)') "DG::  dimension = ", ndim
   write(6,'(A,1PE16.8)') "DG::  diss_tol = ", diis_tol
   write(6,'(A/,(5(F10.5,2X)))') "DG::  Current alpha values = ", alpha(:)
   write(6,'(A,2X,F10.5)') "DG::  Norm of r0 = ", r0_nrm

   if( r0_nrm /= 0.0d0 ) then

      e(:,1) = r0(:) / r0_nrm
      resid(1) = 1.0d0
      coeff(:) = 0.0d0
      k = 0

      do while ( k <= min( ndim-1, ndiis_dim-1 ) .and. resid(k+1) > diis_tol )

!  .............................................................................
!  :  Arnoldi: apply modified Gram-Schmidt orthogonalization to Krylov e_i     :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

         k = k + 1

         call schlegel_DB ( e(:,k), De(:) )
 
         do j = 1, k
            uhess(j,k) = dot_product ( e(:,j), De(:) )
            De(:) = De(:) - uhess(j,k) * e(:,j)
         enddo

         uhess(k+1,k) = sqrt( dot_product( De, De ) ) 

         if( uhess(k+1,k) >= 1.0d-16 .and. k < ndim ) then
            e(:,k+1) = De(:) / uhess(k+1,k)
         endif

!  +---------------------------------------------------------------------------+
!  |  QR decomposition of H: (Q_k)^T R_k = H_k                                 |
!  +---------------------------------------------------------------------------+
      
         qr(1:k+1,1:k) = uhess(1:k+1,1:k)

         call dgeqrf ( k+1, k, qr(1:k+1,1:k), k+1, tau(1:k), work, k, info )

         if( info /= 0 ) then
            write(6,'(A)') "ERROR encountered during dgeqrf in schlegel_gmres"
            write(6,'(A,I8)') "info = ", info
            call mexit(6,1)
         endif

         call dorgqr ( k+1, k, k, qr(1:k+1,1:k), k+1, tau(1:k), work, k, info )

         if( info /= 0 ) then
            write(6,'(A)') "ERROR encountered during dorgqr in schlegel_gmres"
            write(6,'(A,I8)') "info = ", info
            call mexit(6,1)
         endif

!  +---------------------------------------------------------------------------+
!  |  resid_k = Q_k || F - D B_0 || v                                          |
!  +---------------------------------------------------------------------------+

         resid(k+1) = abs( dot_product ( qr(:,k), r0_nrm * v(1:k+1) ) )

         write(6,'(A,I4,A,1PE16.8)') "| DG::  residual(", k+1, ") = ", resid(k+1)

      enddo

      if( resid(k+1) > diis_tol ) then
         write(6,'(A,I8,A)') "| DG:: GMRES failed to converge in ", ndiis_dim &
                           , " steps."
         call mexit(6,1)   
      else
         write(6,'(A,I8,A)') "| DG:: GMRES converged in ", k+1, " steps."
         write(6,'(A,I4,A,1PE16.8)') "| DG::  final residual(", k+1, ") = " &
                                   , resid(k+1)
      endif

!  +---------------------------------------------------------------------------+
!  |  ( H^T H ) C = ( H^T || F - D B_0 || v )                                  |
!  +---------------------------------------------------------------------------+

      HTH(1:k,1:k) = matmul( transpose( uhess(1:k+1,1:k) ), uhess(1:k+1,1:k) ) 

      coeff(1:k) = matmul( transpose( uhess(1:k+1,1:k) ), r0_nrm * v(1:k+1) )

      call dsysv ( 'L', k, 1, HTH(1:k,1:k), k, ipiv(1:k), coeff(1:k), k, work &
                 , k*128, info )

      if( info /= 0 ) then
         write(6,'(A)') "ERROR encountered during dsysv in schlegel_gmres"
         write(6,'(A,I8)') "info = ", info
         call mexit(6,1)
      endif

!  +---------------------------------------------------------------------------+
!  |  B = B_0 + E^T C                                                          |
!  +---------------------------------------------------------------------------+

      b_gmres(:,k+1) = b0(:) + matmul( e(:,1:k), coeff(1:k) )
      b(:) = b_gmres(:,k+1)

#ifdef DEBUG_EVB
      if( schlegel_debug ) then
         open ( io, file = "GMRES_coeff.debug" )
         write( io, '( A )' ) "[Data Size]"
         write( io, '( I20 )' ) k
         write( io, '( A )' ) "[Coeff Vector]"
         write( io, '(5(1PE16.8))' ) ( coeff(j), j = 1, k )
         close( io )
      endif  
#endif

   endif 

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine schlegel_gmres


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  EVB potential and derivatives (using a partial set of coordinates)     |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   subroutine schlegel_EVBR ( q, vkl_sq, V_EVB, dV_EVB, ddV_EVB )

   use schlegel, only: ndg, ncoord, rdgdim, b, vmm, dvmm, ddvmm 
   use evb_parm, only: evb_dcrypt
   use evb_math, only: outer

   implicit none

#  include "parallel.h"
   include 'mpif.h'
#  include "extra.h"

   _REAL_ , intent( in) :: q(ncoord)
   _REAL_ , intent(out) :: vkl_sq, V_EVB, dV_EVB(ncoord), ddV_EVB(ncoord,ncoord)

   !  ..........................................................................

   integer :: k, l, n, alloc_error, dealloc_error
   _REAL_  :: vsub, vsqrt
   intrinsic :: sqrt, reshape
   _REAL_,  allocatable :: dv(:), ddv(:,:)
   _REAL_,  allocatable :: dvkl_sq(:), ddvkl_sq(:,:), ddvkl_sqDG(:)
   _REAL_,  allocatable :: dvsub(:), ddvsub(:,:), dvsqrt(:), ddvsqrt(:,:)

!  .............................................................................
!  :  Allocate temporary arrays                                                :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   allocate( dv(ncoord), ddv(ncoord,ncoord), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( dvkl_sq(ncoord), ddvkl_sq(ncoord,ncoord), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( ddvkl_sqDG(1+ncoord+ncoord*ncoord), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( dvsub(ncoord), ddvsub(ncoord,ncoord), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( dvsqrt(ncoord), ddvsqrt(ncoord,ncoord), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

!  +---------------------------------------------------------------------------+
!  |  Coupling term and derivatives                                            |
!  +---------------------------------------------------------------------------+

   call ddv12sqDG ( q, b(1:rdgdim*ndg), ddvkl_sqDG )

     vkl_sq      = ddvkl_sqDG(1)
    dvkl_sq(  :) = ddvkl_sqDG(2:ncoord+1)
   ddvkl_sq(:,:) = reshape( ddvkl_sqDG(ncoord+2:ncoord+1+ncoord*ncoord) &
                 , (/ ncoord, ncoord /) )

!  +---------------------------------------------------------------------------+
!  |  EVB potential and derivatives                                            |
!  :...........................................................................:
!  |  V_EVB = s_kl - sqrt( d_kl * d_kl + v_kl^2 )                              |
!  |                                                                           |
!  |  where s_kl = 0.50 * ( v_kk + v_ll )                                      |
!  |        d_kl = 0.50 * ( v_kk - v_ll )                                      |
!  +---------------------------------------------------------------------------+

   n = 1
   k = evb_dcrypt(n,1)
   l = evb_dcrypt(n,2)

!  call schlegel_vmm ( q, v, dv, ddv )

!  if( master ) then 
!     call mpi_allgather ( v, 1, MPI_DOUBLE_PRECISION, vmm, 1 &
!                        , MPI_DOUBLE_PRECISION, commmaster, ierr )
!     call mpi_allgather ( dv, ncoord, MPI_DOUBLE_PRECISION, dvmm, ncoord &
!                        , MPI_DOUBLE_PRECISION, commmaster, ierr )
!     call mpi_allgather ( ddv, ncoord*ncoord, MPI_DOUBLE_PRECISION, ddvmm &
!                        , ncoord*ncoord, MPI_DOUBLE_PRECISION, commmaster, ierr )
!  endif

     vsub      = 0.50d0 * (   vmm(    k) -   vmm(    l) )
    dvsub(  :) = 0.50d0 * (  dvmm(  :,k) -  dvmm(  :,l) )
   ddvsub(:,:) = 0.50d0 * ( ddvmm(:,:,k) - ddvmm(:,:,l) )

!  .............................................................................
!  :  Patch the case where the argument of the sqrt is negative                :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   if( vsub**2 + vkl_sq < 0.0d0 ) then

      write(6,'(A)') "WARNING:: about to take SQRT(-#) in DG GMRES"
      write(6,'(A)') "Setting vsqrt and higher derivatives to 0.0"
      write(6,'(A)') "Check routine schlegel_gmres for details"
        vsqrt      = 0.0d0
       dvsqrt(  :) = 0.0d0
      ddvsqrt(:,:) = 0.0d0

!  .............................................................................
!  :  Otherwise, perform the analytical expression                             :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   else

        vsqrt      = sqrt( vsub**2 + vkl_sq )
       dvsqrt(  :) = ( dvsub(:) * vsub +  0.50d0 * dvkl_sq(:) ) / vsqrt
      ddvsqrt(:,:) = ( ddvsub(:,:) * vsub + outer( dvsub(:), dvsub(:)  ) &
                   + 0.50d0 * ddvkl_sq(:,:) ) / vsqrt &
                   - outer( dvsqrt(:), dvsqrt(:) ) / vsqrt
   endif

     V_EVB      = 0.50d0 * (   vmm(    k) +   vmm(    l) ) -   vsqrt
    dV_EVB(  :) = 0.50d0 * (  dvmm(  :,k) +  dvmm(  :,l) ) -  dvsqrt(  :)
   ddV_EVB(:,:) = 0.50d0 * ( ddvmm(:,:,k) + ddvmm(:,:,l) ) - ddvsqrt(:,:)

!  .............................................................................
!  :  Deallocate temporary arrays                                              :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   deallocate( dv, ddv, dvkl_sq, ddvkl_sq, ddvkl_sqDG, stat = dealloc_error )
   REQUIRE( dealloc_error == 0 )
   deallocate( dvsub, ddvsub, dvsqrt, ddvsqrt, stat = dealloc_error )
   REQUIRE( dealloc_error == 0 )

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine schlegel_EVBR


