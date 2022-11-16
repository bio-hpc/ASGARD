
!************************************************************************
!                              AMBER                                   **
!                                                                      **
!                  Copyright (c) 1986, 1991, 1995                      **
!             Regents of the University of California                  **
!                       All Rights Reserved.                           **
!                                                                      **
!  This software provided pursuant to a license agreement containing   **
!  restrictions on its disclosure, duplication, and use. This software **
!  contains confidential and proprietary information, and may not be   **
!  extracted or distributed, in whole or in part, for any purpose      **
!  whatsoever, without the express written permission of the authors.  **
!  This notice, and the associated author list, must be attached to    **
!  all copies, or extracts, of this software. Any additional           **
!  restrictions set forth in the license agreement also apply to this  **
!  software.                                                           **
!************************************************************************


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine movecm here]
subroutine movecm(c,wt,nat)
   implicit double precision (a-h,o-z)
   
   !     ----- ROUTINE TO MOVE MOLECULE SO THAT ORIGIN IS AT COM
   
   dimension c(2),wt(2),vect1(3)
   
   data zero/0.0d+00/,one/1.0e+00/
   
   write(6,9008)
   sumx = zero
   sumy = zero
   sumz = zero
   sum  = zero
   i3 = 0
   do 10 i =1,nat
      wti = wt(i)
      sumx = sumx+c(i3+1)*wti
      sumy = sumy+c(i3+2)*wti
      sumz = sumz+c(i3+3)*wti
      sum = sum+wti
      i3 = i3+3
   10 continue
   xc = sumx/sum
   yc = sumy/sum
   zc = sumz/sum
   i3 = 0
   do 20 i = 1,nat
      c(i3+1) = c(i3+1)-xc
      c(i3+2) = c(i3+2)-yc
      c(i3+3) = c(i3+3)-zc
      i3 = i3+3
   20 continue
   return
   9008 format(5x,'Coordinate origin is now the molecular center of mass')
end subroutine movecm 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine movecm2 here]
subroutine movecm2(c,nat,ndim,c1)
   implicit double precision (a-h,o-z)
   
   !     ----- ROUTINE TO ORIENT THE MOLECULE ALONG THE PRINCIPAL
   !           AXIS ASSUMING THAT ALL THE ATOMIC MASSES ARE 1.00D0
   
   character(len=1) uplo,jobz
   data uplo,jobz / 'U','V' /
   dimension t(6),vec(9),eig(3),b(30)
   dimension c(*),c1(*),vect1(3)
   dimension ititle(20),ititle1(20)
   
   data zero/0.0d+00/,one/1.0e+00/
   
   sumx = zero
   sumy = zero
   sumz = zero
   sum  = zero
   i3 = 0
   do 100 i =1,nat
      sumx = sumx+c(i3+1)
      sumy = sumy+c(i3+2)
      sumz = sumz+c(i3+3)
      sum = sum+ 1
      i3 = i3+3
   100 continue
   xc = sumx/sum
   yc = sumy/sum
   zc = sumz/sum
   i3 = 0
   do 120 i = 1,nat
      c(i3+1) = c(i3+1)-xc
      c(i3+2) = c(i3+2)-yc
      c(i3+3) = c(i3+3)-zc
      c1(i3+1) = c1(i3+1)-xc
      c1(i3+2) = c1(i3+2)-yc
      c1(i3+3) = c1(i3+3)-zc
      i3 = i3+3
   120 continue
   
   !     ----- NOW REORIENT THE MOLECULE ALONG THE PRINCIPAL AXES -----
   
   sxx = zero
   syy = zero
   szz = zero
   sxy = zero
   sxz = zero
   syz = zero
   i3 = 0
   do 140 i = 1,nat
      x = c(i3+1)
      y = c(i3+2)
      z = c(i3+3)
      xx = x*x
      yy = y*y
      zz = z*z
      xy = x*y
      yz = y*z
      xz = x*z
      sxx = sxx+(yy+zz)
      syy = syy+(xx+zz)
      szz = szz+(xx+yy)
      sxy = sxy+xy
      sxz = sxz+xz
      syz = syz+yz
      i3 = i3+3
   140 continue
   
   t(1) = sxx
   t(2) = -sxy
   t(3) = syy
   t(4) = -sxz
   t(5) = -syz
   t(6) = szz
   
   !     ----- CALCULATE PRINCIPAL AXES OF INERTIA -----
   
   call dspev(jobz,uplo,3,t,eig,vec,3,b,ier)
   if(ier /= 0) write(6,9018) ier
   
   !-- check coord. and reinvert if nec.
   
   vect1(1) = vec(2)*vec(6) - vec(3)*vec(5)
   vect1(2) = vec(3)*vec(4) - vec(1)*vec(6)
   vect1(3) = vec(1)*vec(5) - vec(2)*vec(4)
   dot = (vect1(1)*vec(7)+vect1(2)*vec(8)+vect1(3) &
         *vec(9))
   if(dot < 0) then
      do 150 i=7,9
      150 vec(i) = -vec(i)
   end if
   
   i3 = 0
   do 160 i = 1,nat
      x = c(i3+1)
      y = c(i3+2)
      z = c(i3+3)
      xt = vec(1)*x+vec(2)*y+vec(3)*z
      yt = vec(4)*x+vec(5)*y+vec(6)*z
      zt = vec(7)*x+vec(8)*y+vec(9)*z
      c(i3+1) = xt
      c(i3+2) = yt
      c(i3+3) = zt
      x1 = c1(i3+1)
      y1 = c1(i3+2)
      z1 = c1(i3+3)
      xt1 = vec(1)*x1+vec(2)*y1+vec(3)*z1
      yt1 = vec(4)*x1+vec(5)*y1+vec(6)*z1
      zt1 = vec(7)*x1+vec(8)*y1+vec(9)*z1
      c1(i3+1) = xt1
      c1(i3+2) = yt1
      c1(i3+3) = zt1
      i3 = i3+3
   160 continue
   return
   9018 format(5x,'DSPEV RETURNS ERROR')
end subroutine movecm2 
