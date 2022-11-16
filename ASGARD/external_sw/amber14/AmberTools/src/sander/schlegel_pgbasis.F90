! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Polynomial times a Gaussian Distributed Gaussian components            |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#  include "../include/dprec.fh"

   subroutine v12sqDG ( q, b, f )

   use schlegel, only: ndg, ncoord, ncart, nbond, nangle, ndihed, alpha &
                     , xdat_xdg, nselect, scoord, rdgdim
   use evb_math, only: dqint

   implicit none

   _REAL_ , intent( in) :: q(ncoord), b(rdgdim,ndg)
   _REAL_ , intent(out) :: f

   !  ..........................................................................

   integer :: m, n
   _REAL_  :: v12sqF, dq(nselect), dq_tmp(ncoord) 

   f = 0.0d0
   do n = 1, ndg

      dq_tmp(:) = dqint ( q, xdat_xdg(n)%q, nbond, nangle, ndihed, ncoord, ncart )
      do m = 1, nselect
         dq(m) = dq_tmp( scoord(m) )
      enddo

      call v12sqI ( v12sqF, dq, b(:,n), alpha(n), rdgdim, nselect )

      f = f + v12sqF

   enddo

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine v12sqDG


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  0 & 1st derivative of distributed Gaussian components (full internal   |#
! #|  coordinate representation)                                             |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   subroutine dv12sqDG ( q, b, f )

   use schlegel, only: ndg, ncoord, ncart, nbond, nangle, ndihed, alpha &
                     , xdat_xdg, nselect, scoord, rdgdim
   use evb_math, only: dqint
   implicit none

   _REAL_ , intent( in) :: q(ncoord), b(rdgdim,ndg)
   _REAL_ , intent(out) :: f(ncoord+1)

   !  ..........................................................................

   integer :: m, n
   _REAL_  :: dv12sqF(nselect+1), f_tmp(nselect+1), dq(nselect), dq_tmp(ncoord)

   f_tmp(:) = 0.0d0
   do n = 1, ndg

      dq_tmp(:) = dqint ( q, xdat_xdg(n)%q, nbond, nangle, ndihed, ncoord, ncart )
      do m = 1, nselect
         dq(m) = dq_tmp( scoord(m) )
      enddo

      call dv12sqI ( dv12sqF, dq, b(:,n), alpha(n), rdgdim, nselect )

      f_tmp(:) = f_tmp(:) + dv12sqF(:)

   enddo

   f(:) = 0.0d0
   f(1) = f_tmp(1)
   do n = 1, nselect
      f( scoord(n)+1 ) = f_tmp(n+1)
   enddo

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine dv12sqDG


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  0, 1st & 2nd derivatives of distributed Gaussian components (full      |#
! #|  internal coordinate representation)                                    |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   subroutine ddv12sqDG ( q, b, f )

   use schlegel, only: ndg, ncoord, ncart, nbond, nangle, ndihed, alpha &
                     , xdat_xdg, nselect, scoord, rdgdim 
   use evb_math, only: lt_ndx, dqint

   implicit none

   _REAL_ , intent( in) :: q(ncoord), b(rdgdim,ndg)
   _REAL_ , intent(out) :: f(1+ncoord+ncoord*ncoord) 

   !  ..........................................................................

   integer :: m, n, mm, nn, mdx 
   _REAL_  :: ddv12sqF(rdgdim), f_tmp(rdgdim) 
   _REAL_  :: dq(nselect), h_tmp(ncoord,ncoord), dq_tmp(ncoord)


   f_tmp(:) = 0.0d0
   do n = 1, ndg

      dq_tmp(:) = dqint ( q, xdat_xdg(n)%q, nbond, nangle, ndihed, ncoord, ncart )
      do m = 1, nselect
         dq(m) = dq_tmp( scoord(m) )
      enddo

      call ddv12sqI ( ddv12sqF, dq, b(:,n), alpha(n), rdgdim, nselect )

      f_tmp(:) = f_tmp(:) + ddv12sqF(:)

   enddo

!  .............................................................................
!  :  Accumulate components in f(:)                                            :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   f(:) = 0.0d0
   f(1) = f_tmp(1)
   do n = 1, nselect
      f( scoord(n)+1 ) = f_tmp(n+1)
   enddo

   h_tmp(:,:) = 0.0d0
   do n = 1, nselect
      nn = scoord(n)
      do m = 1, nselect
         mm = scoord(m)
         mdx = lt_ndx(  m,  n )
         h_tmp(mm,nn) = f_tmp( 1 + nselect + mdx )
         h_tmp(nn,mm) = h_tmp(mm,nn)
      enddo
   enddo

   f( 1+ncoord+1:1+ncoord+ncoord*ncoord ) = reshape( h_tmp(:,:) &
                                          , (/ ncoord * ncoord /) )

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
 
   end subroutine ddv12sqDG


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  0 & 1st derivative of distributed Gaussian components (partial         |#
! #|  internal coordinate representation)                                    |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   subroutine dv12sqDGR ( q, b, f )

   use schlegel, only: ndg, ncoord, ncart, nbond, nangle, ndihed, alpha &
                     , xdat_xdg, nselect, scoord, rdgdim
   use evb_math, only: dqint

   implicit none

   _REAL_ , intent( in) :: q(ncoord), b(rdgdim,ndg)
   _REAL_ , intent(out) :: f(nselect+1)

   !  ..........................................................................

   integer :: m, n
   _REAL_  :: dv12sqF(nselect+1), dq(nselect), dq_tmp(ncoord) 


   f(:) = 0.0d0
   do n = 1, ndg

      dq_tmp(:) = dqint ( q, xdat_xdg(n)%q, nbond, nangle, ndihed, ncoord, ncart )
      do m = 1, nselect
         dq(m) = dq_tmp( scoord(m) )
      enddo

      call dv12sqI ( dv12sqF, dq, b(:,n), alpha(n), rdgdim, nselect )

      f(:) = f(:) + dv12sqF(:) 

   enddo

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine dv12sqDGR


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  0, 1st & 2nd derivative of distributed Gaussian components (partial    |#
! #|  internal coordinate representation)                                    |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   subroutine ddv12sqDGR ( q, b, f )

   use schlegel, only: ndg, ncoord, ncart, nbond, nangle, ndihed, alpha &
                     , xdat_xdg, nselect, scoord, rdgdim 
   use evb_math, only: dqint

   implicit none

   _REAL_ , intent( in) :: q(ncoord), b(rdgdim,ndg)
   _REAL_ , intent(out) :: f(rdgdim) 

   !  ..........................................................................

   integer :: m, n
   _REAL_  :: ddv12sqF(rdgdim), dq(nselect), dq_tmp(ncoord)


   f(:) = 0.0d0
   do n = 1, ndg

      dq_tmp(:) = dqint ( q, xdat_xdg(n)%q, nbond, nangle, ndihed, ncoord, ncart )
      do m = 1, nselect
         dq(m) = dq_tmp( scoord(m) )
      enddo
 
      call ddv12sqI ( ddv12sqF, dq, b(:,n), alpha(n), rdgdim, nselect )
      
      f(:) = f(:) + ddv12sqF(:)

   enddo

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine ddv12sqDGR


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Polynomial times a Gaussian with weight B                              |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

!  .............................................................................
!  :  V_12^2(q) = sum(k) sum(i>=j>=0) B_ijK g(q,q_K,i,j,a_k)                   :
!  :                                                                           :
!  :  g(DQ,a,0,0) = ( 1 + 0.50 a_K DQ^2 ) exp(-0.5 a_K DQ^2)                   :
!  :  g(DQ,a,i,0) = DQ_i exp(-0.5 a_K DQ^2)                                    :
!  :  g(DQ,a,i,j) = DQ_i * DQ_j exp(-0.5 a_K DQ^2)                             :
!  :                                                                           :
!  :  DQ =  q - q_K      dgdim = 1 + ncoord + ncoord * ( ncoord + 1 ) / 2      :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   subroutine v12sqI ( v12sqF, dq, b, a, dgdim, ncoord ) 

   implicit none

   integer, intent( in) :: dgdim, ncoord
   _REAL_ , intent( in) :: dq(ncoord), b(dgdim), a
   _REAL_ , intent(out) :: v12sqF

   !  ..........................................................................

   integer :: i, j, ii
   _REAL_  :: dqSQ, poly, tmp
   intrinsic :: dot_product, exp


   dqSQ = dot_product( dq(:), dq(:) )

   poly = b(1) * ( 1.0d0 + 0.50d0 * a * dqSQ ) &
        + dot_product( b(2:ncoord+1), dq(:) )

   do i = 2, ncoord
      do j = 1, i - 1
         ii = i * ( i - 1 ) / 2 + j + ncoord + 1
         poly = poly + b(ii) * dq(i) * dq(j)
      enddo
   enddo

   tmp = 0.0d0
   do i = 1, ncoord
      ii = i * ( i + 1 ) / 2 + ncoord + 1
      tmp = tmp + b(ii) * dq(i) * dq(i)
   enddo
   poly = poly + 0.50d0 * tmp

   v12sqF = poly * exp( - 0.50d0 * a * dqSQ )

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine v12sqI 


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Function & 1st derivative of polynomial times a Gaussian with weight B |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

!  .............................................................................
!  :  (d/dq_i) V_12^2(q) = sum(k) sum(i>=j>=0) B_ijK (d/dq_i) g(q,q_K,i,j,a_k) :
!  :                     = sum(k) sum(i>=j>=0) B_ijK                           :
!  :                     [ (d/dq_i) poly exp(-0.5 a_K DQ^2)                    :
!  :                     - poly a_K DQ_i exp(-0.5 a_K DQ^2) ]                  :
!  :                                                                           :
!  :  poly = { ( 1 + 0.50 a_K DQ^2 ); DQ_i; DQ_i DQ_j }                        :
!  :  (d/dq_i) poly = [ { a_K DQ_i; 1; DQ_j } ]_i                              :
!  :                                                                           :
!  :  DQ =  q - q_K      dgdim = 1 + ncoord + ncoord * ( ncoord + 1 ) / 2      :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   subroutine dv12sqI ( dv12sqF, dq, b, a, dgdim, ncoord )

   implicit none

   integer, intent( in) :: dgdim, ncoord
   _REAL_ , intent( in) :: dq(ncoord), b(dgdim), a
   _REAL_ , intent(out) :: dv12sqF(ncoord+1)

   !  ..........................................................................

   integer :: i, j, ii, jj 
   _REAL_  :: dpoly(ncoord), dqSQ, poly, gauss0, tmp
   intrinsic :: dot_product, exp


   dqSQ = dot_product( dq(:), dq(:) )

   poly = b(1) * ( 1.0d0 + 0.50d0 * a * dqSQ ) &
        + dot_product( b(2:ncoord+1), dq(:) )

   do i = 2, ncoord
      do j = 1, i - 1
         ii = i * ( i - 1 ) / 2 + j + ncoord + 1
         poly = poly + b(ii) * dq(i) * dq(j)
      enddo
   enddo

   tmp = 0.0d0
   do i = 1, ncoord
      ii = i * ( i + 1 ) / 2 + ncoord + 1
      tmp = tmp + b(ii) * dq(i) * dq(i)
   enddo
   poly = poly + 0.50d0 * tmp

   dpoly(:) = b(1) * a * dq(:) + b(2:ncoord+1)

   ii = 0
   do i = 1, ncoord
      ii = ii + 1
      do j = 1, i
         jj = i * ( i - 1 ) / 2 + j + ncoord + 1
         dpoly(ii) = dpoly(ii) + b(jj) * dq(j)
      enddo
      do j = i + 1, ncoord
         jj = j * ( j - 1 ) / 2 + i + ncoord + 1
         dpoly(ii) = dpoly(ii) + b(jj) * dq(j)
      enddo
   enddo

   gauss0 = exp( - 0.50d0 * a * dqSQ )

   dv12sqF(         1) = poly * gauss0
   dv12sqF(2:ncoord+1) = ( dpoly(:) - poly * a * dq(:) ) * gauss0

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine dv12sqI


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Function, 1st & 2nd derivatives of polynomial times a Gaussian times B |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

!  .............................................................................
!  :  (d/dq_i) (d/dq_j) V_12^2(q)                                              :
!  :          = sum(k) sum(i>=j>=0) B_ijK (d/dq_i) (d/dq_j) g(q,q_K,i,j,a_k)   :
!  :          = sum(k) sum(i>=j>=0) B_ijK [ (d/dq_i) (d/dq_j) poly             :
!  :            exp(-0.5 a_K DQ^2) - [ (d/dq_j) poly a_K DQ_i                  :
!  :          - poly a_K ] exp(-0.5 a_K DQ^2)                                  :
!  :          - [ (d/dq_i) poly a_K DQ_j + poly (a_K DQ_i) ( a_K DQ_j) ]       :
!  :            exp(-0.5 a_K DQ^2) ]                                           :
!  :                                                                           :
!  :  poly = { ( 1 + 0.50 a_K DQ^2 ); DQ_i; DQ_i DQ_j }                        :
!  :  (d/dq_i) poly = [ { a_K DQ_i; 1; DQ_j } ]_i                              :
!  :  (d/dq_i) (d/dq_j) poly = [{ a_K, for j=i; 0; 1 } ]_i                     :
!  :                                                                           :
!  :  DQ =  q - q_K      dgdim = 1 + ncoord + ncoord * ( ncoord + 1 ) / 2      :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   subroutine ddv12sqI ( ddv12sqF, dq, b, a, dgdim, ncoord )

   implicit none

   integer, intent( in) :: dgdim, ncoord
   _REAL_ , intent( in) :: dq(ncoord), b(dgdim), a
   _REAL_ , intent(out) :: ddv12sqF(dgdim)

   !  ..........................................................................

   integer :: i, j, ii, jj
   _REAL_  :: ddpoly(ncoord*(ncoord+1)/2), ddpoly_tmp(ncoord*(ncoord+1)/2) &
            , dpoly(ncoord), poly, dqSQ, gauss0, tmp
   intrinsic :: dot_product, exp


   dqSQ = dot_product( dq(:), dq(:) )

   poly = b(1) * ( 1.0d0 + 0.50d0 * a * dqSQ ) & 
        + dot_product( b(2:ncoord+1), dq(:) )

   do i = 2, ncoord
      do j = 1, i - 1
         ii = i * ( i - 1 ) / 2 + j + ncoord + 1
         poly = poly + b(ii) * dq(i) * dq(j)
      enddo
   enddo

   tmp = 0.0d0
   do i = 1, ncoord
      ii = i * ( i + 1 ) / 2 + ncoord + 1
      tmp = tmp + b(ii) * dq(i) * dq(i)
   enddo
   poly = poly + 0.50d0 * tmp

   dpoly(:) = b(1) * a * dq(:) + b(2:ncoord+1)

   ii = 0
   do i = 1, ncoord
      ii = ii + 1
      do j = 1, i
         jj = i * ( i - 1 ) / 2 + j + ncoord + 1
         dpoly(ii) = dpoly(ii) + b(jj) * dq(j)
      enddo
      do j = i + 1, ncoord
         jj = j * ( j - 1 ) / 2 + i + ncoord + 1
         dpoly(ii) = dpoly(ii) + b(jj) * dq(j)
      enddo
   enddo

   jj = 0
   ddpoly(:) = 0.0d0
   do i = 1, ncoord
      do j = 1, i
         jj = jj + 1
         if( i == j ) ddpoly(jj) = ddpoly(jj) + b(1) * a
         ii = i * ( i - 1 ) / 2 + j + ncoord + 1
         ddpoly(jj) = ddpoly(jj) + b(ii)
      enddo
   enddo

   jj = 0
   ddpoly_tmp(:) = 0.0d0
   do i = 1, ncoord
      do j = 1, i
         jj = jj + 1
         if( i == j ) ddpoly_tmp(jj) = ddpoly_tmp(jj) + poly
         ddpoly_tmp(jj) = ddpoly_tmp(jj) + dq(i) * dpoly(j) &
                         + dq(j) * dpoly(i) - a * poly * dq(i) * dq(j) 
      enddo
   enddo

   gauss0 = exp( - 0.50d0 * a * dqSQ )

   ii = ncoord * ( ncoord + 1 ) / 2
   ddv12sqF(                   1) = poly * gauss0
   ddv12sqF(          2:ncoord+1) = (  dpoly(:) - poly * a * dq(:) ) * gauss0
   ddv12sqF(ncoord+2:ncoord+1+ii) = ( ddpoly(:) - ddpoly_tmp(:) * a ) * gauss0

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine ddv12sqI



