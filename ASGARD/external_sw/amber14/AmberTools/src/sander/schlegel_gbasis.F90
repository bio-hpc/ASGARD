! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Distributed Gaussian components                                        |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#  include "../include/dprec.fh"

   subroutine schlegel_gauss ( q, gcomp )

   use schlegel, only: ndg, nbond, nangle, ndihed, ncoord, ncart, dgdim, alpha &
                     , xdat_xdg
   use evb_math, only: dqint

   implicit none

   _REAL_ , intent( in) :: q(ncoord)
   _REAL_ , intent(out) :: gcomp(dgdim*ndg)

   !  ..........................................................................

   integer :: i, j, n, ndx
   _REAL_  :: dq(ncoord), dqSQ, a
   _REAL_  , external  :: gaussf 
   intrinsic :: dot_product 

   ndx = 0

   do n = 1, ndg

      a = alpha(n)
      dq(:) = dqint ( q, xdat_xdg(n)%q, nbond, nangle, ndihed, ncoord, ncart )

      dqSQ = dot_product( dq(:), dq(:) )

      ndx = ndx + 1
      gcomp(ndx) = gaussf ( dq, dqSQ, a, 0, 0 )

      do i = 1, ncoord
         ndx = ndx + 1
         gcomp(ndx) = gaussf ( dq, dqSQ, a, i, 0 )
      enddo
 
      do i = 1, ncoord
         do j = 1, ncoord
            ndx = ndx + 1
            gcomp(ndx) = 0.50d0 * gaussf ( dq, dqSQ, a, i, j )
         enddo
      enddo 

   enddo

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine schlegel_gauss


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  1st derivative of distributed Gaussian components                      |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   subroutine schlegel_dgauss ( q, dgcomp, k )

   use schlegel, only: ndg, nbond, nangle, ndihed, ncoord, ncart, dgdim, alpha &
                     , xdat_xdg
   use evb_math, only: dqint

   implicit none

   integer, intent( in) :: k
   _REAL_ , intent( in) :: q(ncoord)
   _REAL_ , intent(out) :: dgcomp(dgdim*ndg)

   !  ..........................................................................

   integer :: i, j, n, ndx
   _REAL_  :: dq(ncoord), dqSQ, a
   _REAL_  , external  :: dgaussf
   intrinsic :: dot_product

   ndx = 0

   do n = 1, ndg

      a = alpha(n)
      dq(:) = dqint ( q, xdat_xdg(n)%q, nbond, nangle, ndihed, ncoord, ncart )

      dqSQ = dot_product( dq(:), dq(:) )

      ndx = ndx + 1
      dgcomp(ndx) = dgaussf ( dq, dqSQ, a, 0, 0, k )

      do i = 1, ncoord
         ndx = ndx + 1
         dgcomp(ndx) = dgaussf ( dq, dqSQ, a, i, 0, k )
      enddo

      do i = 1, ncoord
         do j = 1, ncoord
            ndx = ndx + 1
            dgcomp(ndx) = 0.50d0 * dgaussf ( dq, dqSQ, a, i, j, k )
         enddo
      enddo

   enddo

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine schlegel_dgauss


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  2nd derivative of distributed Gaussian components                      |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   subroutine schlegel_ddgauss ( q, ddgcomp, k, l )

   use schlegel, only: ndg, nbond, nangle, ndihed, ncoord, ncart, dgdim, alpha &
                     , xdat_xdg
   use evb_math, only: dqint

   implicit none

   integer, intent( in) :: k, l
   _REAL_ , intent( in) :: q(ncoord)
   _REAL_ , intent(out) :: ddgcomp(dgdim*ndg)

   !  ..........................................................................

   integer :: i, j, n, ndx
   _REAL_  :: dq(ncoord), dqSQ, a
   _REAL_  , external  :: ddgaussf
   intrinsic :: dot_product

   ndx = 0

   do n = 1, ndg

      a = alpha(n)
      dq(:) = dqint ( q, xdat_xdg(n)%q, nbond, nangle, ndihed, ncoord, ncart )

      dqSQ = dot_product( dq(:), dq(:) )

      ndx = ndx + 1
      ddgcomp(ndx) = ddgaussf ( dq, dqSQ, a, 0, 0, k, l )

      do i = 1, ncoord
         ndx = ndx + 1
         ddgcomp(ndx) = ddgaussf ( dq, dqSQ, a, i, 0, k, l )
      enddo

      do i = 1, ncoord
         do j = 1, ncoord
            ndx = ndx + 1
            ddgcomp(ndx) = 0.50d0 * ddgaussf ( dq, dqSQ, a, i, j, k, l )
         enddo
      enddo

   enddo

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine schlegel_ddgauss


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Compute s-, p-, d-type gaussian functions                              |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

!  .............................................................................
!  :  g(DQ,a,0,0) = exp ( - 0.5 * a * | DQ |^2 )                               :
!  :  g(DQ,a,i,0) = g(DQ,a,0,0) * DQ_i                                         :
!  :  g(DQ,a,i,j) = g(DQ,a,0,0) * DQ_i * DQ_j                                  :
!  :  DQ =  q - q_K                                                            :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   function gaussf ( dq, dqSQ, a, i, j ) 

   implicit none

   integer, intent(in) :: i, j 
   _REAL_ , intent(in) :: dq(*), dqSQ, a

   !  ..........................................................................

   _REAL_ :: gaussf, dq_i, dq_j
   intrinsic :: exp 

   if( i == 0 ) then
      dq_i = 1.0d0
   else
      dq_i = dq(i)
   endif 

   if( j == 0 ) then
      dq_j = 1.0d0
   else
      dq_j = dq(j)
   endif

   gaussf = dq_i * dq_j * exp( - 0.50d0 * a * dqSQ )

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end function gaussf


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Compute 1st derivatives of s-, p-, d-type gaussian functions           |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

!  .............................................................................
!  :  (d/dDQ_k) g(DQ,a,i,j) = - a * DQ_k * g(DQ,a,i,j)                         :
!  :  (d/dDQ_{k=i}) g(DQ,a,i,j) = (d/dDQ_k) g(DQ,a,i,j) + g(DQ,a,0,j)          :
!  :  (d/dDQ_{k=j}) g(DQ,a,i,j) = (d/dDQ_k) g(DQ,a,i,j) + g(DQ,a,i,0)          :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   function dgaussf ( dq, dqSQ, a, i, j, k )

   implicit none

   integer, intent(in) :: i, j, k 
   _REAL_ , intent(in) :: dq(*), dqSQ, a

   !  ..........................................................................

   _REAL_ :: dgaussf
   _REAL_ , external :: gaussf

   dgaussf = - a * dq(k) * gaussf ( dq, dqSQ, a, i, j )

   if( k == i ) then 
      dgaussf = dgaussf + gaussf ( dq, dqSQ, a, 0, j )
   endif 

   if( k == j ) then
      dgaussf = dgaussf + gaussf ( dq, dqSQ, a, i, 0 )
   endif 

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end function dgaussf


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Compute 2nd derivatives of s-, p-, d-type gaussian functions           |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

!  .............................................................................
!  :  (d/dDQ_k) (d/dDQ_k) g(DQ,a,i,j) = - a * g(DQ,a,i,j)                      :
!  :                                    - a * DQ_k dg(DQ,a,i,j,l)              :
!  :  (d/dDQ{k=i}) (d/dDQ_k) g(DQ,i,0,a) = (d/dDQ_k) (d/dDQ_k)                 :
!  :                         g(DQ,a,i,j) + dg(dQ,a,0,j,l)                      :
!  :  (d/dDQ{k=j}) (d/dDQ_k) g(DQ,i,0,a) = (d/dDQ_k) (d/dDQ_k)                 :
!  :                         g(DQ,a,i,j) + dg(DQ,a,0,j,l)                      :
!  :                                     + dg(DQ,a,i,0,l)                      :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   function ddgaussf ( dq, dqSQ, a, i, j, k, l )

   implicit none

   integer, intent(in) :: i, j, k, l 
   _REAL_ , intent(in) :: dq(*), dqSQ, a

   !  ..........................................................................

   _REAL_ :: ddgaussf
   _REAL_ , external :: gaussf, dgaussf

   ddgaussf = - a * dq(k) * dgaussf ( dq, dqSQ, a, i, j, l )

   if( l == k ) ddgaussf = ddgaussf - a * gaussf ( dq, dqSQ, a, i, j ) 

   if( k == i ) ddgaussf = ddgaussf +  dgaussf ( dq, dqSQ, a, 0, j, l )

   if( k == j ) ddgaussf = ddgaussf +  dgaussf ( dq, dqSQ, a, i, 0, l )

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end function ddgaussf


