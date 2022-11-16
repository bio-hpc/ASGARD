#include "../include/dprec.fh"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine project here]
subroutine project(l,m,n,h,hx,hy,hz,i,j,k,x,y,z,phi, x1,y1,z1)
  implicit none
   
   !       **************************************************************
   !       *                                                            *
   !       *  project find the projection of the interface of a given   *
   !       *  grid point (x_i,y_j,z_k)                                  *
   !       *                                                            *
   !       **************************************************************

   !common /lmn/l, m, n, nirreg
   integer l,m,n,i,j,k
   !common /hxyz/hx,hy,hz,hmax
  
  !passed variables
   _REAL_ h,hx,hy,hz,x1,y1,z1
   _REAL_ x(0:l+1),y(0:m+1),z(0:n+1)
   _REAL_ phi(0:l+1,0:m+1,0:n+1)
  
  !local variables
   _REAL_ hmax,phx,phy,phz,phxx,phyy,phzz,phxy,phyz,phxz,phn,phn1,&
          px,py,pz,temp1,temp2,temp3,a0,a1,a2,r1,r2,r3,r4,r5,drop,&
          dproj,xx,yy,zz
   integer info

   ! ----- Compute the deepest descent direction (-grad phi(i,j,k)), then
   !       locates the control points on the level set

   hmax=h
   !ph0 = phi(i,j,k)
   !       write(100,*) phi
   call phidv1(l,m,n,h,hx,hy,hz,i,j,k,x,y,z,phi, phx,phy,phz)
   call phidv2(l,m,n,h,hx,hy,hz,i,j,k,x,y,z,phi, phxx,phyy,phzz,phxy,phxz,phyz)
 
   
   phn = phx*phx + phy*phy + phz*phz
   phn1 = sqrt(phn)

   if (abs(phn1) < 1.0e-20) then
      write(*,*) "   Error: phn1 = 0.0 in project()!"
      stop
   end if

   px = phx/phn1
   py = phy/phn1
   pz = phz/phn1
   temp1 = phxx*px + phxy*py + phxz*pz
   temp2 = phxy*px + phyy*py + phyz*pz
   temp3 = phxz*px + phyz*py + phzz*pz

   a0 = phi(i,j,k)
   a1 = -phn1
   a2 =( px*temp1+py*temp2+pz*temp3)*0.5

   call rootp2(a0,a1,a2, r1,r2,info)

   ! ----- If there is at least one root

   if (info > 0) then
      r3 = min(r1,r2)
      r4 = max(r1,r2)

      if (phi(i,j,k) >= 0.0) then
         if (r3 >= 0.0) then
            r5 = r3
         else
            r5 = r4
         end if
      else
         if (r4 <= 0.0) then
            r5 = r4
         else
            r5 = r3
         end if
      end if

      x1 = x(i) - r5*px
      y1 = y(j) - r5*py
      z1 = z(k) - r5*pz

      dproj =  (x1-x(i))*(x1-x(i)) + (y1-y(j))*(y1-y(j)) &
            + (z1-z(k))*(z1-z(k))

   end if

   ! ----- if no real root is found, or it is too far away (hmax),
   !       function may be too bad, choose the intersection between
   !       the interface and the grid lines

   if (info <= 0 .or. sqrt(dproj) > 2.0*hmax) then
      if (phi(i,j,k)*phi(i+1,j,k) < 0.0) then
         if (i == 0) then
            info = 0
            call update(x(i),x(i+1),x(i+2), phi(i,j,k), &
                  phi(i+1,j,k),phi(i+2,j,k),info,xx)
         else
            info = 1
            call update(x(i-1),x(i),x(i+1), phi(i-1,j,k), &
                  phi(i,j,k),phi(i+1,j,k),info,xx)
         end if
         if (info == 1) then
            x1 = xx
         else
            x1 = x(i) + hx*phi(i,j,k)/(phi(i,j,k)-phi(i+1,j,k))
         end if
         y1 = y(j)
         z1 = z(k)
         goto 777
      end if

      if (phi(i,j,k)*phi(i-1,j,k) < 0.0) then
         if (i-1 == 0) then
            info = 0
            call update(x(i-1),x(i),x(i+1), phi(i-1,j,k), &
                  phi(i,j,k),phi(i+1,j,k),info,xx)
         else
            info = 1
            call update(x(i-2),x(i-1),x(i), phi(i-2,j,k), &
                  phi(i-1,j,k),phi(i,j,k),info,xx)
         end if
         if (info == 1) then
            x1 = xx
         else
            x1 = x(i) - hx*phi(i,j,k)/(phi(i,j,k)-phi(i-1,j,k))
         end if
         y1 = y(j)
         z1 = z(k)
         goto 777
      end if

      if (phi(i,j,k)*phi(i,j+1,k) < 0.0) then
         if (j == 0) then
            info = 0
            call update(y(j),y(j+1),y(j+2), phi(i,j,k), &
                  phi(i,j+1,k),phi(i,j+2,k),info,yy)
         else
            info = 1
            call update(y(j-1),y(j),y(j+1), phi(i,j-1,k), &
                  phi(i,j,k),phi(i,j+1,k),info,yy)
         end if
         if (info == 1) then
            y1 = yy
         else
            y1 = y(j) + hy*phi(i,j,k)/(phi(i,j,k)-phi(i,j+1,k))
         end if
         x1 = x(i)
         z1 = z(k)
         goto 777
      end if

      if (phi(i,j,k)*phi(i,j-1,k) < 0.0) then
         if (j-1 == 0) then
            info = 0
            call update(y(j-1),y(j),y(j+1), phi(i,j-1,k), &
                  phi(i,j,k),phi(i,j+1,k),info,yy)
         else
            info = 1
            call update(y(j-2),y(j-1),y(j), phi(i,j-2,k), &
                  phi(i,j-1,k),phi(i,j,k),info,yy)
         end if
         if (info == 1) then
            y1 = yy
         else
            y1 = y(j) - hy*phi(i,j,k)/(phi(i,j,k)-phi(i,j-1,k))
         end if
         x1 = x(i)
         z1 = z(k)
         goto 777
      end if

      if (phi(i,j,k)*phi(i,j,k+1) < 0.0) then
         if (k == 0) then
            info = 0
            call update(z(k),z(k+1),z(k+2), phi(i,j,k), &
                  phi(i,j,k+1),phi(i,j,k+2),info,zz)
         else
            info = 1
            call update(z(k-1),z(k),z(k+1), phi(i,j,k-1), &
                  phi(i,j,k),phi(i,j,k+1),info,zz)
         end if
         if (info == 1) then
            z1 = zz
         else
            z1 = z(k) + hz*phi(i,j,k)/(phi(i,j,k)-phi(i,j,k+1))
         end if
         x1 = x(i)
         y1 = y(j)
         goto 777
      end if

      if (phi(i,j,k)*phi(i,j,k-1) < 0.0) then
         if (k-1 == 0) then
            info = 0
            call update(z(k-1),z(k),z(k+1), phi(i,j,k-1), &
                  phi(i,j,k),phi(i,j,k+1),info,zz)
         else
            info = 1
            call update(z(k-2),z(k-1),z(k), phi(i,j,k-2), &
                  phi(i,j,k-1),phi(i,j,k),info,zz)
         end if
         if (info == 1) then
            z1 = zz
         else
            z1 = z(k) - hz*phi(i,j,k)/(phi(i,j,k)-phi(i,j,k-1))
         end if
         x1 = x(i)
         y1 = y(j)
         goto 777
      end if

   end if  ! (info <= 0 .or. sqrt(dproj) > 2.0*hmax)

   777 continue

   if (abs(x1-x(i)) >= hx .or. abs(y1-y(j)) >= hy &
         .or. abs(z1-z(k)) >= hz ) then
      if (phi(i,j,k)*phi(i+1,j,k) < 0.0) then
         x1 = x(i) + hx*phi(i,j,k)/(phi(i,j,k)-phi(i+1,j,k))
         y1 = y(j)
         z1 = z(k)
         return
      end if

      if (phi(i,j,k)*phi(i-1,j,k) < 0.0) then
         x1 = x(i) - hx*phi(i,j,k)/(phi(i,j,k)-phi(i-1,j,k))
         y1 = y(j)
         z1 = z(k)
         return
      end if

      if (phi(i,j,k)*phi(i,j+1,k) < 0.0) then
         x1 = x(i)
         y1 = y(j) + hy*phi(i,j,k)/(phi(i,j,k)-phi(i,j+1,k))
         z1 = z(k)
         return
      end if

      if (phi(i,j,k)*phi(i,j-1,k) < 0.0) then
         x1 = x(i)
         y1 = y(j) - hy*phi(i,j,k)/(phi(i,j,k)-phi(i,j-1,k))
         z1 = z(k)
         return
      end if

      if (phi(i,j,k)*phi(i,j,k+1) < 0.0) then
         x1 = x(i)
         y1 = y(j)
         z1 = z(k) + hz*phi(i,j,k)/(phi(i,j,k)-phi(i,j,k+1))
         return
      end if

      if (phi(i,j,k)*phi(i,j,k-1) < 0.0) then
         x1 = x(i)
         y1 = y(j)
         z1 = z(k) - hz*phi(i,j,k)/(phi(i,j,k)-phi(i,j,k-1))
         return
      end if

   end if  !  (abs(x1-x(i)) >= hx .or. abs(y1-y(j)) >= hy

   !       call GrToPr(x,y,z,x1,y1,z1,i0,j0,k0,phi,1,0,0,
   !     1              q0,qx,qy,qz,qxx,qyy,qzz,qxy,qxz,qyz)

   !        if(abs(q0) .lt. 1.0e-12) return

   ! ----- Call newsol() to update (x1,y1,z1) so that it is indeed on the interface

   !       call newsol(x,y,z,x1,y1,z1,phi)

  ! write(1010,*) x1,y1,z1,px,py,pz !used to collect projection info
   return
end subroutine project 




!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine newsol here]
!XP: newsol not called..
subroutine newsol(l,m,n,h,hx,hy,hz,x,y,z,x1,y1,z1,phi)
  implicit _REAL_ (a-h,o-z) !Not used, not optimized.
  !implicit none
   
   !       **************************************************************
   !       *                                                            *
   !       *  newsol uses Newton method to update (x1,y1,z1) so that it *
   !       *  is indeed on the interface.                               *
   !       *                                                            *
   !       **************************************************************

   !common /lmn/l, m, n, nirreg
   integer l,m,n
   !common /hxyz/hx,hy,hz,hmax

   _REAL_ x(0:l+1),y(0:m+1),z(0:n+1)
   _REAL_ phi(0:l+1,0:m+1,0:n+1)

   iter = 0
   call grtopr(l,m,n,x,y,z,h,hx,hy,hz,x1,y1,z1,phi,1,1,0, &
         q0,qx,qy,qz,qxx,qyy,qzz,qxy,qxz,qyz)
   111 if (qx /= 0.0) then
      x1 = x1 - q0/qx
      call grtopr(l,m,n,x,y,z,h,hx,hy,hz,x1,y1,z1,phi,1,1,0, &
            q0,qx,qy,qz,qxx,qyy,qzz,qxy,qxz,qyz)
      if(abs(q0) < 1.0e-14) goto 777
      iter = iter + 1
      if (iter < 101) goto 111
   end if

   iter = 0
   222 if (qy /= 0.0) then
      y1 = y1 - q0/qy
      call grtopr(l,m,n,x,y,z,h,hx,hy,hz,x1,y1,z1,phi,1,1,0, &
            q0,qx,qy,qz,qxx,qyy,qzz,qxy,qxz,qyz)
      if(abs(q0) < 1.0e-14) goto 777
      iter = iter + 1
      if (iter < 101) goto 222
   end if

   iter = 0
   333 if (qz /= 0.0) then
      z1 = z1 - q0/qz
      call grtopr(l,m,n,x,y,z,h,hx,hy,hz,x1,y1,z1,phi,1,1,0, &
            q0,qx,qy,qz,qxx,qyy,qzz,qxy,qxz,qyz)
      if(abs(q0) < 1.0e-14) goto 777
      iter = iter + 1
      if (iter < 101) goto 333
   end if

   write(*,*) "   Warning: Some projection point NOT on the interface!"

   777 return
end subroutine newsol 




!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine newint here]
subroutine newint(x0,x1,x2,y0,y1,y2, a0, a1, a2)
   implicit none
   
   !       **************************************************************
   !       *                                                            *
   !       *  newint uses Newton interpolation to determine the 2-nd    *
   !       *  polynomial using the three points (x0,y0), (x1,y1),       *
   !       *  (x2,y2).                                                  *
   !       *                                                            *
   !       *  The resulted polynomial is p(x)=a2*x^2+a1*x+a0            *
   !       *                                                            *
   !       **************************************************************
  
   !Passed variables
   _REAL_ x0,x1,x2,y0,y1,y2,a0,a1,a2

   !Local variables
   _REAL_ f10,f21,f210

   f10  = (y1 - y0)/(x1 - x0)
   f21  = (y2 - y1)/(x2 - x1)
   f210 = (f21 - f10)/(x2-x0)

   a2 = f210
   a1 = -(x0+x1)*f210 + f10
   a0 = x0*x1*f210 - x0*f10 + y0

   return
end subroutine newint 



!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine update here]
subroutine update(x0,x1,x2,y0,y1,y2,info,xx)
   implicit none
   
   !       **************************************************************
   !       *                                                            *
   !       *  update finds one of two zeros of p(x)=a2*x^2+a1*x+a0      *
   !       *  which is in [x0,x1] if info=0 or in [x1,x2] if info=1     *
   !       *                                                            *
   !       **************************************************************
  
   !Passed variables:
   _REAL_ x0,x1,x2,y0,y1,y2,xx
   integer info
  
   !Local variables:
   _REAL_ a0,a1,a2,r1,r2
   integer info1

   call newint(x0,x1,x2,y0,y1,y2, a0, a1, a2)
   call rootp2(a0,a1,a2, r1,r2,info1)

   if (info1 <= 0 ) then
      info = 0
   else
      if (info .eq. 0) then
         if (r1 >= x0-1.0e-15 .and. r1 <= x1+1.0e-15) then
            xx = r1
            info = 1
            return
         else if (r2 >= x0-1.0e-15 .and. r2 <= x1+1.0e-15) then
            xx = r2
            info = 1
            return
         else
            info = 0
            return
         end if
      else
         if (r1 >= x1-1.0e-15 .and. r1 <= x2+1.0e-15) then
            xx = r1
            info = 1
            return
         else if (r2 >= x1-1.0e-15 .and. r2 <= x2+1.0e-15) then
            xx = r2
            info = 1
            return
         else
            info = 0
            return
         end if
      end if
   end if
   return
end subroutine update 







