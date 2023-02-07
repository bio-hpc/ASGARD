
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
!+ [Enter a one-line description of subroutine bond here]
subroutine bond (nb,   ib,    jb,      icb, &
      rk,   req,   ntb,     x, &
      f,    omega, eb,      ndrv,  nbel,   natsys)
   
   !     ----- routine to get bond energy and forces for the potential
   !           of cb*(b-b0)**2
   
   implicit double precision (a-h,o-z)
   dimension ib(nb),jb(nb),icb(nb),rk(2),req(2),x(2), &
         f(2),omega(2),nbel(2)
   dimension xij(3)
   
   eb = 0.0d+00
   do 40 n = 1,nb
      i3 = ib(n)
      j3 = jb(n)
      ipc = icb(n)
      rij2 = 0.e0
      do 10 m = 1,3
         xij(m) = x(i3+m)-x(j3+m)
         rij2 = rij2+xij(m)**2
      10 continue
      if(ntb /= 0) call percon(rij2,xij)
      rij = sqrt(rij2)
      db = rij-req(ipc)
      df = rk(ipc)*db
      ebh = df*db
      if (ebh > 100.0) then
         write(6,101) i3/3+1,j3/3+1,rij,req(ipc),rk(ipc),ebh
         101 format('Bad bond:',2i5,4e12.5)
      end if
      eb = eb+ebh
      if(ndrv == 2) then
         dv1 = 2.0e0*df
         dv2 = 2.0e0*rk(ipc)
         call difbon(omega,dv1,dv2,xij,rij,i3,j3,nbel)
      end if
      df = 2.0e0*df/rij
      do 20 m = 1,3
         xh = xij(m)*df
         f(i3+m) = f(i3+m)-xh
         f(j3+m) = f(j3+m)+xh
      20 continue
      
   40 continue
   return
end subroutine bond 
