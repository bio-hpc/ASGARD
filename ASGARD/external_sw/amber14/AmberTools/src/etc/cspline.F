c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2003  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine isplpe  --  cubic periodic interpolating spline  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "isplpe" computes the coefficients for a cubic periodic
c     interpolating spline
c
c     literature reference:
c
c     G. Engeln-Mullges and F. Uhlig, Numerical Algorithms with Fortran,
c     Springer Verlag, 1996, Section 10.1.2
c
c
      subroutine isplpe (ng,xn,fn,mrep,b,c,d,h,du,dm,rc,rs)
      implicit none
      integer i,ng,mrep,iflag
      real*8 eps,average
      real*8 temp1,temp2
      real*8 xn(0:ng),fn(0:ng)
      real*8 b(0:ng),c(0:ng)
      real*8 d(0:ng),h(0:ng)
      real*8 du(ng),dm(ng)
      real*8 rc(ng),rs(ng)
c
c
c     check the periodicity of fn, and for subsequent call
c
      eps = 0.00001d0
      if (abs(fn(ng)-fn(0)) .gt. eps) then
         write (6,10)
   10    format (' ISPLPE  --  Warning, Input Values to Spline',
     &              ' are Not Periodic')
      end if
      average = 0.5d0 * (fn(0) + fn(ng))
      fn(0) = average
      fn(ng) = average
      if (mrep.ne.1 .and. mrep.ne.2)  return
c
c     get auxiliary variables and matrix elements on first call
c
      if (mrep .eq. 1) then
         do i = 0, ng-1
            h(i) = xn(i+1) - xn(i)
         end do
         h(ng) = h(0)
         do i = 1, ng-1
            du(i) = h(i)
         end do
         du(ng) = h(0)
         do i = 1, ng
            dm(i) = 2.0d0 * (h(i-1)+h(i))
         end do
      end if
c
c     compute the right hand side
c
      temp1 = (fn(1)-fn(0)) / h(0)
      do i = 1, ng-1, 1
         temp2 = (fn(i+1)-fn(i)) / h(i)
         rs(i)  = 3.0d0 * (temp2-temp1)
         temp1 = temp2
      end do
      rs(ng) = 3.0d0 * ((fn(1)-fn(0))/h(0)-temp1)
c
c     solve the linear system, with factorization on first call
c
      if (mrep .eq. 1) then
         call cytsy (ng,dm,du,rc,rs,c(1),iflag)
         if (iflag .ne. 1)  return
c
c     proceed without factorization on subsequent calls
c
      else
         call cytsys (ng,dm,du,rc,rs,c(1))
      end if
c
c     compute remaining spline coefficients
c
      c(0) = c(ng)
      do i = 0, ng-1
         b(i) = (fn(i+1)-fn(i))/h(i) - h(i)/3.0d0*(c(i+1)+2.0d0*c(i))
         d(i) = (c(i+1)-c(i)) / (3.0d0*h(i))
      end do
      b(ng) = (fn(1)-fn(ng))/h(ng) - h(ng)/3.0d0*(c(1)+2.0d0*c(ng))
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine cytsy  --  solve cyclic tridiagonal system  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "cytsy" solves a system of linear equations for a cyclically
c     tridiagonal, symmetric, positive definite matrix
c
c     literature reference:
c
c     G. Engeln-Mullges and F. Uhlig, Numerical Algorithms with Fortran,
c     Springer Verlag, 1996, Section 4.11.2
c
c
      subroutine cytsy (ng,dm,du,cr,rs,x,iflag)
      implicit none
      integer ng,iflag
      real*8 dm(1:ng),du(1:ng)
      real*8 cr(1:ng),rs(1:ng)
      real*8 x(1:ng)
c
c
c     factorization of the input matrix
c
      iflag = -2
      if (ng .lt. 3)  return
      call cytsyp (ng,dm,du,cr,iflag)
c
c     update and backsubstitute as necessary
c
      if (iflag .eq. 1) then
         call cytsys (ng,dm,du,cr,rs,x)
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine cytsyp  --  tridiagonal Cholesky factorization  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "cytsyp" finds the Cholesky factors of a cyclically tridiagonal
c     symmetric, positive definite matrix given by two vectors
c
c     literature reference:
c
c     G. Engeln-Mullges and F. Uhlig, Numerical Algorithms with Fortran,
c     Springer Verlag, 1996, Section 4.11.2
c
c
      subroutine cytsyp (ng,dm,du,cr,iflag)
      implicit none
      integer i,ng,iflag
      real*8 eps,row,d
      real*8 temp1,temp2
      real*8 dm(1:ng),du(1:ng)
      real*8 cr(1:ng)
c
c
c     set error bound and test for condition ng > 2
c
      eps = 0.00000001d0
      iflag = -2
      if (ng .lt. 3)  return
c
c     checking for positive definite matrix a and for strong
c     nonsingularity of a for ng=1 and if a is
c
      row = abs(dm(1)) + abs(du(1)) + abs(du(ng))
      if (row .eq. 0.0d0) then
         iflag = 0
         return
      end if
      d = 1.0d0 / row
      if (dm(1) .lt. 0.0d0) then
         iflag = -1
         return
      else if (abs(dm(1))*d .le. eps) then
         iflag = 0
         return
      end if
c
c     factoring a while checking for a positive definite and strong
c     nonsingular matrix a
c
      temp1 = du(1)
      du(1) = du(1) / dm(1)
      cr(1) = du(ng) / dm(1)
      do i = 2, ng-1
         row = abs(dm(i)) + abs(du(i)) + abs(temp1)
         if (row .eq. 0.0d0) then
            iflag = 0
            return
         end if
         d = 1.0d0 / row
         dm(i) = dm(i) - temp1*du(i-1)
         if (dm(i) .lt. 0.0d0) then
            iflag = -1
            return
         else if (abs(dm(i))*d .le. eps) then
            iflag = 0
            return
         end if
         if (i .lt. (ng-1)) then
            cr(i) = -temp1 * cr(i-1) / dm(i)
            temp1 = du(i)
            du(i) = du(i) / dm(i)
         else
            temp2 = du(i)
            du(i) = (du(i) - temp1*cr(i-1)) / dm(i)
         end if
      end do
      row = abs(du(ng)) + abs(dm(ng)) + abs(temp2)
      if (row .eq. 0.0d0) then
         iflag = 0
         return
      end if
      d = 1.0d0 / row
      dm(ng) = dm(ng) - dm(ng-1)*du(ng-1)*du(ng-1)
      temp1 = 0.0d0
      do i = 1, ng-2
         temp1 = temp1 + dm(i)*cr(i)*cr(i)
      end do
      dm(ng) = dm(ng) - temp1
      if (dm(ng) .lt. 0) then
         iflag = -1
         return
      else if (abs(dm(ng))*d .le. eps) then
         iflag = 0
         return
      end if
      iflag = 1
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine cytsys  --  tridiagonal solution from factors  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "cytsys" solves a cyclically tridiagonal linear system
c     given the Cholesky factors
c
c     literature reference:
c
c     G. Engeln-Mullges and F. Uhlig, Numerical Algorithms with Fortran,
c     Springer Verlag, 1996, Section 4.11.2
c
c
      subroutine cytsys (ng,dm,du,cr,rs,x)
      implicit none
      integer i,ng
      real*8 sum,temp
      real*8 dm(1:ng),du(1:ng)
      real*8 cr(1:ng),rs(1:ng)
      real*8 x(1:ng)
c
c
c     updating phase
c
      temp = rs(1)
      rs(1) = temp / dm(1)
      sum = cr(1) * temp
      do i = 2, ng-1
         temp = rs(i) - du(i-1)*temp
         rs(i) = temp / dm(i)
         if (i .ne. (ng-1))  sum = sum + cr(i)*temp
      end do
      temp = rs(ng) - du(ng-1)*temp
      temp = temp - sum
      rs(ng) = temp / dm(ng)
c
c     backsubstitution phase
c
      x(ng) = rs(ng)
      x(ng-1) = rs(ng-1) - du(ng-1)*x(ng)
      do i = ng-2, 1, -1
         x(i) = rs(i) - du(i)*x(i+1) - cr(i)*x(ng)
      end do
      return
      end
