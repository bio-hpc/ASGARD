! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Form the Schlegel-Sonnenberg D matrix and F vector and solve for B:    |#
! #|                                                                         |#
! #|                                 D B = F                                 |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#  include "../include/assert.fh"
#  include "../include/dprec.fh"

   subroutine schlegel_dg 

   use schlegel,  only: ndg, ncoord, dgdim, vkl, dvkl, ddvkl, xdat_xdg &
                      , alpha, scoord, nselect, rdgdim, dist_gauss, nUFF &
                      , nAdihed, nAangle, vmm, dvmm, ddvmm 
   use dg_uff,    only: vdw_repulsive_init
   use dg_dihed,   only: amber_dihed_init
   use dg_angle,   only: amber_angle_init
   use file_io_dat
#ifdef DEBUG_EVB
   use evb_check, only: evb_debug, schlegel_debug
   use schlegel, only : b
#endif

   implicit none

#  include "parallel.h"
   include 'mpif.h'
#  include "extra.h"

   !  ..........................................................................

   integer :: k, l, n, ndx, nstart, nend, alloc_error, dealloc_error, ierr
#ifdef DEBUG_EVB
   integer :: io
#endif
   _REAL_  :: v
   intrinsic :: reshape
   _REAL_, allocatable :: q(:), dv(:), ddv(:,:)
   _REAL_, allocatable :: gcomp(:), dgcomp(:), ddgcomp(:)
   _REAL_, allocatable :: f(:), b0(:), d(:,:)

!  .............................................................................
!  :  Allocate temporary arrays                                                :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   allocate( q(ncoord), dv(ncoord), ddv(ncoord,ncoord), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate(   gcomp(dgdim*ndg), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate(  dgcomp(dgdim*ndg), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( ddgcomp(dgdim*ndg), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( f(dgdim*ndg), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( b0(dgdim*ndg), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( d(dgdim*ndg,dgdim*ndg), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

   write(6,'(A)') "DG:: external file name and associated alpha exponent"
   do n = 1, ndg
      write(6,100) n, trim( adjustl( xdat_xdg(n)%filename ) ), alpha(n)
   enddo

!  +---------------------------------------------------------------------------+
!  |  Select coordinates for DG procedure                                      |
!  +---------------------------------------------------------------------------+

   allocate( scoord(ncoord), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

   call select_coord ( dist_gauss%stype, dist_gauss%stol, scoord, nselect )

   if( trim( adjustl( dist_gauss%lin_solve ) ) == "diis" ) then
      rdgdim = 1 + nselect + nselect * ( nselect + 1 ) / 2
   else if( trim( adjustl( dist_gauss%lin_solve ) ) == "full" ) then
      rdgdim = 1 + nselect + nselect * nselect
   endif

   write(6,'(A,I8)') "DG:: No. of DG ab initio configuration points = ", ndg
   write(6,'(A,I8)') "DG:: No. of DG data points within each " &
                  // "configuration point = ", rdgdim

!  +---------------------------------------------------------------------------+
!  |  Initialize UFF VDW parameters (UFF VDW term is added to V_ii)            |
!  +---------------------------------------------------------------------------+

   if( nUFF > 0 ) call vdw_repulsive_init

!  +---------------------------------------------------------------------------+
!  |  Initialize Amber torsion parameters (Amber dihed term is added to V_ii)  |
!  +---------------------------------------------------------------------------+

   if( nAdihed > 0 ) call amber_dihed_init

!  +---------------------------------------------------------------------------+
!  |  Initialize Amber torsion parameters (Amber dihed term is added to V_ii)  |
!  +---------------------------------------------------------------------------+

   if( nAangle > 0 ) call amber_angle_init

!  +---------------------------------------------------------------------------+
!  |  Form F vector:: { V_12^2; (d/dq) V_12^2; (d^2/dq^2) V_12^2 }             |
!  +---------------------------------------------------------------------------+

   do n = 1, ndg

      q(:) = xdat_xdg(n)%q(:)


   call schlegel_vmm ( q, v, dv, ddv  )

   if( master ) then
      call mpi_allgather ( v, 1, MPI_DOUBLE_PRECISION, vmm, 1 &
                         , MPI_DOUBLE_PRECISION, commmaster, ierr )
      call mpi_allgather ( dv, ncoord, MPI_DOUBLE_PRECISION, dvmm, ncoord &
                         , MPI_DOUBLE_PRECISION, commmaster, ierr )
      call mpi_allgather ( ddv, ncoord*ncoord, MPI_DOUBLE_PRECISION, ddvmm &
                         , ncoord*ncoord, MPI_DOUBLE_PRECISION, commmaster, ierr )
   endif


      call schlegel_vkl ( n )
      ndx = rdgdim * ( n - 1 ) + 1
      f(ndx) = vkl

      do k = 1, nselect
         ndx = ndx + 1
         f(ndx) = dvkl( scoord(k) )
      enddo 

      if( trim( adjustl( dist_gauss%lin_solve ) ) == "diis" ) then
         do k = 1, nselect
            do l = 1, k
               ndx = ndx + 1
               f(ndx) = ddvkl( scoord(k),scoord(l) )
            enddo
         enddo
      else if( trim( adjustl( dist_gauss%lin_solve ) ) == "full" ) then
         nstart = ndx + 1 
         nend   = nstart + ncoord * ncoord - 1
         f(nstart:nend) = reshape( ddvkl(:,:), (/ ncoord*ncoord /) )
      endif

   enddo

#ifdef DEBUG_EVB
   if( schlegel_debug ) then
      io = schlegel_unit
      open ( io, file = "F_vector.debug" )
      write( io, '( A )' ) "[Data Size]"
      write( io, '( I20 )' ) rdgdim * ndg
      write( io, '( A )' ) "[F Vector]"
      write( io, '( 5(1PE16.8) )' ) ( f(n), n = 1, rdgdim * ndg )
      close( io )
   endif
#endif

!  +---------------------------------------------------------------------------+
!  |  Solve B using the Krylov subspace approach                               |
!  +---------------------------------------------------------------------------+

   if( trim( adjustl( dist_gauss%lin_solve ) ) == "diis" ) then

!     call schlegel_diis ( f, b, 100, 20 )
      b0(:) = f(:)
      call schlegel_gmres ( b0, f, rdgdim*ndg, 200 )

!  +---------------------------------------------------------------------------+
!  |  Form D matrix from Gaussian centers/derivatives & solve for B via SVD    |
!  +---------------------------------------------------------------------------+

   else if( trim( adjustl( dist_gauss%lin_solve ) ) == "full" ) then

      do n = 1, ndg
         q(:) =xdat_xdg(n)%q(:)
         call schlegel_gauss ( q, gcomp )
         ndx = dgdim * ( n - 1 ) + 1
         d(ndx,:) = gcomp(:)
         do k = 1, ncoord
            call schlegel_dgauss ( q, dgcomp, k )
            d(ndx+k,:) = dgcomp(:)
            do l = 1, ncoord
               call schlegel_ddgauss ( q, ddgcomp, k, l )
               d(ndx+ncoord+ncoord*(k-1)+l,:) = ddgcomp(:)
            enddo
         enddo
      enddo

#ifdef DEBUG_EVB
      if( schlegel_debug ) then
         io = schlegel_unit
         open ( io, file = "D_matrix.debug" )
         write( io, '( A )' ) "[Data Size]"
         write( io, '( I20 )' ) rdgdim * ndg * rdgdim * ndg
         write( io, '( A )' ) "[D Matrix]"
         write( io, '( 5(1PE16.8) )' ) ( d(:,n), n = 1, rdgdim * ndg )
         close( io )
      endif
#endif

      call schlegel_full ( f, d )

   endif

#ifdef DEBUG_EVB
   if( schlegel_debug ) then
      io = schlegel_unit
      open ( io, file = "B_vector.debug" )
      write( io, '( A )' ) "[Data Size]"
      write( io, '( I20 )' ) rdgdim * ndg
      write( io, '( A )' ) "[B Vector]"
      write( io, '( 5(1PE16.8) )' ) ( b(n), n = 1, rdgdim * ndg )
      close( io )
   endif

!  +---------------------------------------------------------------------------+
!  |  Debugging against Schlegel-Sonnenberg Mathematica files                  |
!  +---------------------------------------------------------------------------+

   if( schlegel_debug ) then
!     call evb_vs_abinitio
      if( master ) then
         select case( trim( adjustl( evb_debug%fmathnb ) ) )
            case( "EVB_hcn_int_irc.nb" )
               call schlegel_hcn_int_irc
               write(6,'(A)') "|  DONE generating EVB PES for HCN system " &
                           // "|  (Internal IRC case)."
               write(6,'(A)') "|  Data points are written to EVB_hcn_int_irc.nb"
            case( "EVB_hcn_cart_irc.nb" )
               call schlegel_hcn_cart_irc
               write(6,'(A)') "|  DONE generating EVB PES for HCN system " &
                           // "|  (Cartesian IRC case)."
               write(6,'(A)') "|  Data points are written to EVB_hcn_cart_irc.nb"
            case( "EVB_hcn_cart.nb" )
               call schlegel_hcn_cart
               write(6,'(A)') "|  DONE generating EVB PES for HCN system " &
                           // "|  (Cartesian case)."
               write(6,'(A)') "|  Data points are written to EVB_hcn_cart.nb"
            case( "EVB_poh_irc.nb" )
               call schlegel_poh_irc
               write(6,'(A)') "|  DONE generating EVB PES for POH system " &
                           // "|  (Internal IRC case)."
               write(6,'(A)') "|  Data points are written to EVB_poh_irc.nb"
            case( "EVB_poh_uff.nb" )
               call schlegel_poh_uff
               write(6,'(A)') "|  DONE generating EVB PES for POH system " &
                           // "|  (Internal IRC case with UFF VDW)."
               write(6,'(A)') "|  Data points are written to EVB_poh_uff.nb"
               write(6,'(A)') "|  ab initio PES written to HFgrid_poh.nb"
         end select
      endif
      call mexit(6, 0)
   endif 

#endif

!  .............................................................................
!  :  Deallocate temporary arrays                                              :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   deallocate( q, dv, ddv, stat = dealloc_error )
   REQUIRE( dealloc_error == 0 )
   deallocate( gcomp, dgcomp, ddgcomp, stat = dealloc_error )
   REQUIRE( dealloc_error == 0 )
   deallocate( f, b0, d, stat = dealloc_error )
   REQUIRE( dealloc_error == 0 )

 100 format( I8, 2X, A45, 2X, F12.5 )

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine schlegel_dg


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  DG-EVB potential and derivatives                                       |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   subroutine schlegel_evb ( q, V_EVB, dV_EVB, ddV_EVB )

   use evb_parm, only: evb_dcrypt, nxch
   use schlegel, only: ndg, ncoord, dgdim, b, vmm, dvmm, ddvmm 
   use evb_math, only: outer

   implicit none

#  include "parallel.h"
   include 'mpif.h'
#  include "extra.h"

   _REAL_ , intent( in) :: q(ncoord)
   _REAL_ , intent(out) :: V_EVB, dV_EVB(ncoord), ddV_EVB(ncoord,ncoord)

   !  ..........................................................................

   integer :: k, l, n, alloc_error, dealloc_error
   _REAL_  :: vkl_sq, vsub, vsqrt
   intrinsic :: dot_product, sqrt
   _REAL_, allocatable :: dv(:), ddv(:,:), dvkl_sq(:), ddvkl_sq(:,:)
   _REAL_, allocatable :: dvsub(:), ddvsub(:,:), dvsqrt(:), ddvsqrt(:,:)
   _REAL_, allocatable :: gcomp(:), dgcomp(:), ddgcomp(:)

!  .............................................................................
!  :  Allocate temporary arrays                                                :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   allocate( dv(ncoord), ddv(ncoord,ncoord), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( dvkl_sq(ncoord), ddvkl_sq(ncoord,ncoord), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( dvsub(ncoord), ddvsub(ncoord,ncoord), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( dvsqrt(ncoord), ddvsqrt(ncoord,ncoord), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate(   gcomp(dgdim*ndg), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate(  dgcomp(dgdim*ndg), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( ddgcomp(dgdim*ndg), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
 
!  +---------------------------------------------------------------------------+
!  |  Coupling term and derivatives                                            |
!  |                                                                           |
!  |             V_12^2(q) = sum(K;i>=j>=0) B_ijK g(q,q_K,i,j,a)               |
!  |      (d/dq) V_12^2(q) = sum(K;i>=j>=0) B_ijK (d/dq) g(q,q_K,i,j,a)        |
!  |  (d^2/dq^2) V_12^2(q) = sum(K;i>=j>=0) B_ijK (d^2/dq^2) g(q,q_K,i,j,a)    |
!  +---------------------------------------------------------------------------+

   call schlegel_gauss ( q, gcomp )
   vkl_sq = dot_product( b(:), gcomp(:) )   
  
   do k = 1, ncoord 
      call schlegel_dgauss ( q, dgcomp, k )
      dvkl_sq(k) = dot_product( b(:), dgcomp(:) )
   enddo

   do k = 1, ncoord
      do l = 1, ncoord 
         call schlegel_ddgauss ( q, ddgcomp, k, l )
         ddvkl_sq(k,l) = dot_product( b(:), ddgcomp(:) )
      enddo
   enddo

!  +---------------------------------------------------------------------------+
!  |  EVB potential and derivatives                                            |
!  :...........................................................................:
!  |  V_EVB = s_kl - sqrt( d_kl * d_kl - v_kl^2 )                              |
!  |                                                                           |
!  |  where s_kl = 0.50 * ( v_kk + v_ll )                                      |
!  |        d_kl = 0.50 * ( v_kk - v_ll )                                      |
!  +---------------------------------------------------------------------------+

   do n = 1, nxch

      k = evb_dcrypt(n,1)
      l = evb_dcrypt(n,2)

!     call schlegel_vmm ( q, v, dv, ddv  )

!  if( master ) then 
!     call mpi_allgather ( v, 1, MPI_DOUBLE_PRECISION, vmm, 1 &
!                        , MPI_DOUBLE_PRECISION, commmaster, ierr )
!     call mpi_allgather ( dv, ncoord, MPI_DOUBLE_PRECISION, dvmm, ncoord &
!                        , MPI_DOUBLE_PRECISION, commmaster, ierr )
!     call mpi_allgather ( ddv, ncoord*ncoord, MPI_DOUBLE_PRECISION, ddvmm &
!                        , ncoord*ncoord, MPI_DOUBLE_PRECISION, commmaster, ierr )
!  endif

      vsub = 0.50d0 * ( vmm(k) - vmm(l) )
      dvsub(:) = 0.50d0 * ( dvmm(:,k) - dvmm(:,l) )
      ddvsub(:,:) = 0.50d0 * ( ddvmm(:,:,k) - ddvmm(:,:,l) )

      vsqrt = sqrt( vsub**2 + vkl_sq )
      dvsqrt(:) = (dvsub(:) * vsub +  0.50d0 * dvkl_sq(:) ) / vsqrt

      ddvsqrt(:,:) = ( ddvsub(:,:) * vsub & 
                   + outer( dvsub(:), dvsub(:)  ) &
                   + 0.50d0 * ddvkl_sq(:,:) ) / vsqrt &
                   - outer( dvsqrt(:), dvsqrt(:) ) / vsqrt
 
      V_EVB = 0.50d0 * ( vmm(k) + vmm(l) ) - vsqrt
      dV_EVB(:) = 0.50d0 * ( dvmm(:,k) + dvmm(:,l) ) - dvsqrt(:)
      ddV_EVB(:,:) = 0.50d0 * ( ddvmm(:,:,k) + ddvmm(:,:,l) ) - ddvsqrt(:,:)

   enddo 

!  .............................................................................
!  :  Deallocate temporary arrays                                              :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   deallocate( dv, ddv, dvkl_sq, ddvkl_sq, stat = dealloc_error )
   REQUIRE( dealloc_error == 0 )
   deallocate( dvsub, ddvsub, dvsqrt, ddvsqrt, stat = dealloc_error )
   REQUIRE( dealloc_error == 0 )
   deallocate( gcomp, dgcomp, ddgcomp, stat = dealloc_error )
   REQUIRE( dealloc_error == 0 )

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine schlegel_evb


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  EVB coupling & derivatives based on a polynomial times a Gaussian form |#
! #|                                                                         |#
! #|  V_12^2(q) = A[1 + B^T*dq + 0.5 dq^T*(C+aI)*dq] exp(-0.5 a dq^2)        |#
! #|                                                                         |#
! #|  A = [V_11(q) - V(q)][V_22(q) - V(q)]                                   |#
! #|  B = D_1/(V_11-V) + D_2/(V_22-V)      D_i = (d/dq) V_ii - (d/dq) V      |#
! #|  C = (D_1 D_2^T + D_2 D_1^T)/A + K_1/(V11-V) + K_2/(V_22-V)             |#
! #|                                  K_i = (d^2/dq^2) V_ii - (d^2/dq^2) V   |#
! #|                                                                         |#
! #|  Note that for the components of the F vector, dq = 0                   |#
! #|  ... simplifying V_12^2 = A                                             |#
! #|           (d/dq) V_12^2 = D_1 (V_22-V) + D_2 (V_11-V)                   |#
! #|       (d^2/dq^2) V_12^2 = K_1 (V_22-V) + D_1 D_2^T + K_2 (V_11-V)       |#
! #|                                          D_2 D_1^T                      |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   subroutine schlegel_vkl ( dg_pt )

   use evb_parm, only: evb_dcrypt, nxch
   use schlegel, only: ncoord, vkl, dvkl, ddvkl, vmm, dvmm, ddvmm, xdat_xdg
   use evb_math, only: outer

   implicit none

#  include "parallel.h"
   include 'mpif.h'
#  include "extra.h"

   integer, intent(in) :: dg_pt

   !............................................................................

   integer :: k, l, n, alloc_error, dealloc_error
   _REAL_  :: nrg, dEk, dEl
   _REAL_, allocatable :: dv(:), ddv(:,:), grad(:), hess(:,:), dGk(:), dGl(:)

!  .............................................................................
!  :  Allocate temporary arrays                                                :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   allocate( dv(ncoord), ddv(ncoord,ncoord), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( grad(ncoord), hess(ncoord,ncoord), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( dGk(ncoord), dGl(ncoord), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

!  +---------------------------------------------------------------------------+
!  |  Harmonic expansion of potential about ab initio minimum                  |
!  +---------------------------------------------------------------------------+

!  call schlegel_vmm ( q, v, dv, ddv  )

!  if( master ) then 
!     call mpi_allgather ( v, 1, MPI_DOUBLE_PRECISION, vmm, 1 &
!                        , MPI_DOUBLE_PRECISION, commmaster, ierr )
!     call mpi_allgather ( dv, ncoord, MPI_DOUBLE_PRECISION, dvmm, ncoord &
!                        , MPI_DOUBLE_PRECISION, commmaster, ierr )
!     call mpi_allgather ( ddv, ncoord*ncoord, MPI_DOUBLE_PRECISION, ddvmm &
!                        , ncoord*ncoord, MPI_DOUBLE_PRECISION, commmaster, ierr )
!  endif

!  +---------------------------------------------------------------------------+
!  |  Coupling and derivatives for 2-state EVB                                 |
!  +---------------------------------------------------------------------------+

   do n = 1, nxch

      k = evb_dcrypt(n,1)
      l = evb_dcrypt(n,2)

      nrg       = xdat_xdg(dg_pt)%v
      grad(  :) = xdat_xdg(dg_pt)%d(  :)
      hess(:,:) = xdat_xdg(dg_pt)%k(:,:)

      dEk = vmm(k) - nrg
      dEl = vmm(l) - nrg

      dGk(:) = dvmm(:,k) - grad(:)
      dGl(:) = dvmm(:,l) - grad(:)

      vkl = dEk * dEl

      dvkl(:) = dGk(:) * dEl + dGl(:) * dEk

      ddvkl(:,:) = outer( dGk(:), dGl(:) ) + outer( dGl(:), dGk(:) ) &
                 + ( ddvmm(:,:,k) - hess(:,:) ) * dEl &
                 + ( ddvmm(:,:,l) - hess(:,:) ) * dEk

   enddo

!  .............................................................................
!  :  Deallocate temporary arrays                                              :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   deallocate( dv, ddv, stat = dealloc_error )
   REQUIRE( dealloc_error == 0 )
   deallocate( grad, hess, stat = dealloc_error )
   REQUIRE( dealloc_error == 0 )
   deallocate( dGk, dGl, stat = dealloc_error )
   REQUIRE( dealloc_error == 0 )

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine schlegel_vkl


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Quadratic expansion of potential about ab initio minimum               |#
! #|                                                                         |#
! #|  V_mm(q) = V_0 + G*(q-q_0) + 0.5 (q-q_0)^T*H*(q-q_0)                    |#
! #|  G_mm(q) = G + (q-q_0)^T*H                                              |#
! #|  H_mm(q) = H                                                            |#
! #|                                                                         |#
! #|  ... where          G = gradient      H = Hessian                       |# 
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   subroutine schlegel_vmm ( q, v, dv, ddv )

   use schlegel, only: nbond, nangle, ndihed, ncoord, ncart, xdat_min, nUFF &
                     , nAdihed, nAangle
   use dg_uff,   only: vdw_repulsive, dvdw_repulsive, ddvdw_repulsive
   use dg_dihed,   only: amber_dihed, damber_dihed, ddamber_dihed
   use dg_angle,   only: amber_angle, damber_angle, ddamber_angle
   use evb_math, only: dqint

   implicit none

#  include "parallel.h"

   _REAL_ , intent( in) :: q(ncoord)
   _REAL_ , intent(out) :: v, dv(ncoord), ddv(ncoord,ncoord)

   !  ..........................................................................

   integer :: n, alloc_error, dealloc_error
   _REAL_  :: nrg
   intrinsic :: dot_product, matmul
   _REAL_, allocatable :: grad(:), hess(:,:), dq(:)

!  .............................................................................
!  :  Allocate temporary arrays                                                :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   allocate( grad(ncoord), hess(ncoord,ncoord), dq(ncoord), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

   n = masterrank + 1

   dq(:) = dqint ( q, xdat_min(n)%q, nbond, nangle, ndihed, ncoord, ncart )

   nrg       = xdat_min(n)%v
   grad(  :) = xdat_min(n)%d(  :)
   hess(:,:) = xdat_min(n)%k(:,:)

   v = nrg + dot_product( grad(:), dq(:) ) &
           + 0.50d0 * dot_product( dq(:), matmul( hess(:,:), dq(:) ) )

   dv(:) = grad(:) + matmul( hess(:,:), dq(:) )

   ddv(:,:) = hess(:,:)

   if( nUFF > 0 ) then

      v = v + vdw_repulsive(q,ncoord) - vdw_repulsive(xdat_min(n)%q,ncoord)

!     write(6,*) 'n = ', n
!     write(6,*) 'vdw_repulsive = ', vdw_repulsive(q,ncoord) 
!     write(6,*) 'D vdw_repulsive = ', vdw_repulsive(xdat_min(n)%q,ncoord) 
!     write(6,*) 'v = ', v

      dv(:) = dv(:) + dvdw_repulsive(q,ncoord) 
!- dvdw_repulsive(xdat_min(n)%q,ncoord)

      ddv(:,:) = ddv(:,:) + ddvdw_repulsive(q,ncoord) 
!- ddvdw_repulsive(xdat_min(n)%q,ncoord)

!     write(100,*) 'nUFF = ', nUFF
!     write(100,*) dvdw_repulsive(q,ncoord)
!     write(100,*)

!     write(200,*) 'nUFF = ', nUFF
!     write(200,*) dvdw_repulsive(xdat_min(n)%q,ncoord)
!     write(200,*)

   endif

   if( nAdihed > 0 ) then
      v = v + amber_dihed(q,ncoord) - amber_dihed(xdat_min(n)%q,ncoord)
      dv(:) = dv(:) + damber_dihed(q,ncoord)
      ddv(:,:) = ddv(:,:) + ddamber_dihed(q,ncoord)
   endif

   if( nAangle > 0 ) then
      v = v + amber_angle(q,ncoord) - amber_angle(xdat_min(n)%q,ncoord)
      dv(:) = dv(:) + damber_angle(q,ncoord)
      ddv(:,:) = ddv(:,:) + ddamber_angle(ncoord)
   endif

!  .............................................................................
!  :  Deallocate temporary arrays                                              :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
   deallocate( grad, hess, dq, stat = dealloc_error )
   REQUIRE( dealloc_error == 0 )

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine schlegel_vmm


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Allocate storage space for schlegel_dg                                 |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   subroutine schlegel_dg_alloc

   use evb_parm, only: nevb
   use schlegel, only: ndg, dgdim, ncoord, natm, nUFF, vmm, dvmm, ddvmm, dvkl &
                     , ddvkl, b, gbasis_weight, wdcmat &
                     , xvdw, dvdw, zvdw, Avdw, Bvdw

   implicit none

   !  ..........................................................................

   integer :: alloc_error

   allocate( vmm(nevb), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( dvmm(ncoord,nevb), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( ddvmm(ncoord,ncoord,nevb), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

   allocate( dvkl(ncoord), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( ddvkl(ncoord,ncoord), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

   allocate( b(dgdim*ndg), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

   allocate( wdcmat(ncoord,natm*3), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

!  allocate( alpha(ndg), stat = alloc_error )
!  REQUIRE( alloc_error == 0 )

!  alpha(:) = 0.50d0
!  alpha(:) = 0.80d0

   if( gbasis_weight(1) == 9999.9d0 ) then
      gbasis_weight(1) = 1.0d0 
      gbasis_weight(2) = 1.0d0 / dble( ncoord )
      gbasis_weight(3) = 1.0d0 / dble( ncoord * ncoord )
   endif

   if( nUFF > 0 ) then
      allocate( xvdw(nUFF), stat = alloc_error )
      REQUIRE( alloc_error == 0 )
      allocate( dvdw(nUFF), stat = alloc_error )
      REQUIRE( alloc_error == 0 )
      allocate( zvdw(nUFF), stat = alloc_error )
      REQUIRE( alloc_error == 0 )
      allocate( Avdw(nUFF), stat = alloc_error )
      REQUIRE( alloc_error == 0 )
      allocate( Bvdw(nUFF), stat = alloc_error )
      REQUIRE( alloc_error == 0 )
   endif

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
 
   end subroutine schlegel_dg_alloc 


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Deallocate storage space for schlegel_dg                               |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   subroutine schlegel_dg_dealloc

   use schlegel, only: vmm, dvmm, ddvmm, dvkl, ddvkl, b, wdcmat, alpha, dg_weight

   implicit none

   !  ..........................................................................

   integer :: dealloc_error

   if (allocated(alpha)) deallocate(alpha)
   deallocate( vmm, dvmm, ddvmm, dvkl, ddvkl, b, wdcmat, &
               dg_weight, stat = dealloc_error )
   REQUIRE( dealloc_error == 0 )

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine schlegel_dg_dealloc


