
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
!+ [Enter a one-line description of subroutine dihed here]
subroutine dihed (x,igraph,idih)
   
   !     calculates angle for a dihedral i-j-k-l
   !       the coordinates for atom i are assumed to be present
   !       in positions x([i-1]*3+1), x([i-1]*3+2),x([i-1]*3+2)
   
   implicit double precision (a-h,o-z)
   dimension x(*),igraph(*),idih(*)
   pi = 3.14159d0
   
   if(idih(1) == 0) return
   
   write(6,100)
   do 10 i=1,400,4
      if(idih(i) == 0) return
      
      i3 = (idih(i)-1)*3
      j3 = (idih(i+1)-1)*3
      k3 = (idih(i+2)-1)*3
      l3 = (idih(i+3)-1)*3
      
      xij = x(i3+1)-x(j3+1)
      yij = x(i3+2)-x(j3+2)
      zij = x(i3+3)-x(j3+3)
      xkj = x(k3+1)-x(j3+1)
      ykj = x(k3+2)-x(j3+2)
      zkj = x(k3+3)-x(j3+3)
      xkl = x(k3+1)-x(l3+1)
      ykl = x(k3+2)-x(l3+2)
      zkl = x(k3+3)-x(l3+3)
      
      dx = yij*zkj - zij*ykj
      dy = zij*xkj - xij*zkj
      dz = xij*ykj - yij*xkj
      gx = zkj*ykl - ykj*zkl
      gy = xkj*zkl - zkj*xkl
      gz = ykj*xkl - xkj*ykl
      bi = dx**2 + dy**2 + dz**2
      bk = gx**2 + gy**2 + gz**2
      ct = dx*gx + dy*gy + dz*gz
      ct = ct/sqrt(bi*bk)
      if (ct < -1.) ct=-1.
      if (ct > 1.)  ct= 1.
      ap = acos(ct)
      d  = xkj*(dz*gy-dy*gz) + ykj*(dx*gz-dz*gx) + zkj*(dy*gx-dx*gy)
      if (d < 0) ap=-ap
      ap = pi-ap
      g = ap
      
      app = 180.00d0*ap/pi
      if (app > 180.d0) app = app -360.d0
      write (6,110) (igraph(idih(j)),idih(j), j=i,i+3),app
      
   10 continue
   
   return
   100 format(5x,' Dihedral angles ')
   110 format(5x,3(a5,'(',i3,') - '),a5,'(',i3,') = ',f10.5)
end subroutine dihed 
