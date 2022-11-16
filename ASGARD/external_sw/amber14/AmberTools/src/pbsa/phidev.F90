#include "../include/dprec.fh"

!       ------------------------------------------
!       Computer derivatives of level set function
!       ------------------------------------------


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine phidv1 here]
subroutine phidv1(l,m,n,h,hx,hy,hz,i,j,k,x,y,z,phi, phx,phy,phz)
   implicit none
   
   !       **************************************************************
   !       *                                                            *
   !       *  phidv1  find the  first  derivatives  of a  level set     *
   !       *  phi(x,y,z)  at a grid  point (x_i,y_j,z_k)  using the     *
   !       *  central deifference.                                    *
   !       *                                                            *
   !       *  assumption: we don''t know the exact expression of the     *
   !       *  function phi(x,y,z), but its values at each grid point    *
   !       *                                                            *
   !       **************************************************************
   
   !Passed variables
   integer l,m,n,i,j,k
   _REAL_ h,hx,hy,hz,phx,phy,phz
   _REAL_ x(0:l+1),y(0:m+1),z(0:n+1)
   _REAL_ phi(0:l+1,0:m+1,0:n+1)

   if (i == l+1) then
      phx = (3*phi(l+1,j,k)+phi(l-1,j,k)-4*phi(l,j,k))/(2.0*hx)
   else
      if (i == 0) then
         phx = (4*phi(1,j,k)-3*phi(0,j,k)-phi(2,j,k))/(2.0*hx)
      else
         phx = (phi(i+1,j,k) - phi(i-1,j,k))/(2.0*hx)
      end if
   end if

   if (j == m+1) then
      phy = (3*phi(i,m+1,k)+phi(i,m-1,k)-4*phi(i,m,k))/(2.0*hy)
   else
      if (j == 0) then
         phy = (4*phi(i,1,k)-3*phi(i,0,k)-phi(i,2,k))/(2.0*hy)
      else
         phy = (phi(i,j+1,k) - phi(i,j-1,k))/(2.0*hy)
      end if
   end if

   if (k == n+1) then
      phz = (3*phi(i,j,n+1)+phi(i,j,n-1)-4*phi(i,j,n))/(2.0*hz)
   else
      if (k == 0) then
         phz = (4*phi(i,j,1)-3*phi(i,j,0)-phi(i,j,2))/(2.0*hz)
      else
         phz = (phi(i,j,k+1) - phi(i,j,k-1))/(2.0*hz)
      end if
   end if

   return
end subroutine phidv1 




!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine phidv2 here]
subroutine phidv2(l,m,n,h,hx,hy,hz,i,j,k,x,y,z,phi, &
      phxx,phyy,phzz,phxy,phxz,phyz)
  implicit none
   
   !       **************************************************************
   !       *                                                            *
   !       *  phidv2  find the second  derivatives  of a  level set     *
   !       *  phi(x,y,z)  at a grid  point (x_i,y_j,z_k)  using the     *
   !       *  central deifference.                                      *
   !       *                                                            *
   !       *  assumption: we don''t know the exact expression of the     *
   !       *  function phi(x,y,z), but its values at each grid point    *
   !       *                                                            *
   !       **************************************************************

   !Passed variables
   integer l,m,n,i,j,k
   _REAL_ hx,hy,hz,h,phxx,phyy,phzz,phxy,phxz,phyz
   _REAL_ x(0:l+1),y(0:m+1),z(0:n+1)
   _REAL_ phi(0:l+1,0:m+1,0:n+1)

   if (i == 0 .or. i == l+1) return
   if (j == 0 .or. j == m+1) return
   if (k == 0 .or. k == n+1) return

   if (i == 0) then
      phxx=(phi(i,j,k)+phi(i+2,j,k)-2.0*phi(i+1,j,k))/hx/hx
   else if (i == l+1) then
      phxx=(phi(i-2,j,k)+phi(i,j,k)-2.0*phi(i-1,j,k))/hx/hx
   else
      phxx=(phi(i+1,j,k)+phi(i-1,j,k)-2.0*phi(i,j,k))/hx/hx
   end if

   if (j == 0) then
      phyy=(phi(i,j,k)+phi(i,j+2,k)-2.0*phi(i,j+1,k))/hy/hy
   else if (j == m+1) then
      phyy=(phi(i,j-2,k)+phi(i,j,k)-2.0*phi(i,j-1,k))/hy/hy
   else
      phyy=(phi(i,j+1,k)+phi(i,j-1,k)-2.0*phi(i,j,k))/hy/hy
   end if

   if (k == 0) then
      phzz=(phi(i,j,k)+phi(i,j,k+2)-2.0*phi(i,j,k+1))/hz/hz
   else if (k == n+1) then
      phzz=(phi(i,j,k-2)+phi(i,j,k)-2.0*phi(i,j,k-1))/hz/hz
   else
      phzz=(phi(i,j,k+1)+phi(i,j,k-1)-2.0*phi(i,j,k))/hz/hz
   end if

   phxy = (phi(i+1,j+1,k)+phi(i-1,j-1,k)-phi(i+1,j-1,k) &
         -phi(i-1,j+1,k))/(4.0*hx*hy)
   phxz = (phi(i+1,j,k+1)+phi(i-1,j,k-1)-phi(i+1,j,k-1) &
         -phi(i-1,j,k+1))/(4.0*hx*hz)
   phyz = (phi(i,j+1,k+1)+phi(i,j-1,k-1)-phi(i,j+1,k-1) &
         -phi(i,j-1,k+1))/(4.0*hy*hz)

   return 
 end subroutine phidv2 



