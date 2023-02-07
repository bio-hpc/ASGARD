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
!+ [Enter a one-line description of subroutine bmat here]
subroutine bmat(b,iat,jat,kat,lat,f,nint,ncart,x,avint)
   
   ! ---  sets up vibrational b-matrix of dimension n x n
   !       first index is internal coordinate, second is cartesian.
   !       order of the internal coordinates is same as in force list
   
   ! --- also sets up arrays iat,jat,kat,lat identifying each of the
   !       nint internal coordinates by its atom numbers, and f array,
   !       specifying the force constant along each internal coordinate
   
   implicit double precision(a-h,o-z)
#  include "sizes2.h"
#  include "bad.h"
   
   common/minpar/ntrun,maxcyc,ncyc,iopt,nvect,izzz,dxm,dele,drms
   common/conver/ipt(maxint)
   dimension b(ncart,nint),iat(nint),jat(nint),kat(nint), &
         lat(nint),f(nint),avint(nint)
   dimension x(*),ibs(100),jbs(100)
   
   ! --- zero out b matrix
   
   do 6 i=1,nint
      do 5 j=1,ncart
         b(j,i) = 0.00d0
      5 continue
   6 continue
   ratodeg = 180./3.14159
   
   kb = 0
#ifdef Z_AXIS_FLUC
   
   !---here is trial code to compute fluctuations of angle this
   !       vector makes with the z-axis
   !       Let f array contain the equilibrium angles!
   !      first internal is theta, second is phi.
   
   if(ntrun == -1) then
      do mm=1,nint/2
         i = 3*(ibs(mm)) - 3
         j = 3*(jbs(mm)) - 3
         dx = x(j+1) - x(i+1)
         dy = x(j+2) - x(i+2)
         dz = x(j+3) - x(i+3)
         s = dx*dx + dy*dy + dz*dz
         s = 1./sqrt(s)
         cost = dz*s
         sint = sqrt(1. - cost*cost)
         den = s**3/sint
         kb = kb + 1
         b(i+1,kb) = dx*dz*den
         b(i+2,kb) = dy*dz*den
         b(i+3,kb) = -(dx*dx+dy*dy)*den
         b(j+1,kb) = -dx*dz*den
         b(j+2,kb) = -dy*dz*den
         b(j+3,kb) = (dx*dx+dy*dy)*den
         iat(kb) = i/3 + 1
         jat(kb) = j/3 + 1
         kat(kb) = 0
         lat(kb) = 0
         f(kb) = dacosd(cost)
         ipt(kb) = 2
         
         sp = dx*dx + dy*dy
         dp = 1./sqrt(sp)
         cosp = dx*dp
         sinp = dy*dp
         den = s**3/sinp
         kb = kb + 1
         b(i+1,kb) = -dy*dy*den
         b(i+2,kb) = dx*dy*den
         b(j+1,kb) = dy*dy*den
         b(j+2,kb) = -dx*dy*den
         iat(kb) = i/3 + 1
         jat(kb) = j/3 + 1
         kat(kb) = 0
         lat(kb) = 0
         f(kb) = datan2d(dy,dx)
         ipt(kb) = 2
         
      end do  !  mm=1,nint/2
      return
   end if  ! (ntrun == -1)
#endif
   
   if(ntrun == 1 .or. ntrun == 0) then
      nbond = nbonh + nbona
      do 20 mm = 1,nbond
         if (mm <= nbonh) then
            i = ibh(mm)
            j = jbh(mm)
            ic = icbh(mm)
         else
            i = iba(mm-nbonh)
            j = jba(mm-nbonh)
            ic = icba(mm-nbonh)
         end if
         dx = x(i+1) - x(j+1)
         dy = x(i+2) - x(j+2)
         dz = x(i+3) - x(j+3)
         s = dx*dx + dy*dy + dz*dz
         s = 1./sqrt(s)
         dx = dx*s
         dy = dy*s
         dz = dz*s
         kb = kb + 1
         b(i+1,kb) = dx
         b(i+2,kb) = dy
         b(i+3,kb) = dz
         b(j+1,kb) = -dx
         b(j+2,kb) = -dy
         b(j+3,kb) = -dz
         iat(kb) = i/3 + 1
         jat(kb) = j/3 + 1
         kat(kb) = 0
         lat(kb) = 0
         avint(kb) = 1./s
         if (ic > 0) then
            f(kb) = rk(ic)
         else
            f(kb) = 1.0
         end if
         ipt(kb) = 1
      20 continue
   end if  !  20 mm = 1,nbond
   
   if(ntrun == 0) return
   if(ntrun <= 2) then
      nthets = ntheth + ntheta
      do 30 ithe=1,nthets
         if (ithe <= ntheth) then
            i3=ith(ithe)
            j3=jth(ithe)
            k3=kth(ithe)
            ic = icth(ithe)
         else
            kthe = ithe - ntheth
            i3 = ita(kthe)
            j3 = jta(kthe)
            k3 = kta(kthe)
            ic = icta(kthe)
         end if
         dx=x(i3+1)-x(j3+1)
         dy=x(i3+2)-x(j3+2)
         dz=x(i3+3)-x(j3+3)
         gx=x(k3+1)-x(j3+1)
         gy=x(k3+2)-x(j3+2)
         gz=x(k3+3)-x(j3+3)
         bij=dx*dx+dy*dy+dz*dz
         bkj=gx*gx+gy*gy+gz*gz
         cst=dx*gx+dy*gy+dz*gz
         
         boi2=1./bij
         boj2=1./bkj
         bij=sqrt(bij)
         bkj=sqrt(bkj)
         a=1./bij
         b1=1./bkj
         bioj=bij*b1
         bjoi=bkj*a
         cst=cst*a*b1
         if (cst > 1. ) cst= 1.
         if (cst < (-1.)) cst=-1.
         at=acos(cst)
         st=1.-cst*cst
         !      DERIVATIVE UNDEFINED AT THETA=180.
         if(abs(st) < 0.000001) st=0.000001
         st=-1.0/sqrt(st)
         boi2=st*boi2
         boj2=st*boj2
         baa=boi2*bioj
         bab=boi2*cst
         bdd=boj2*bjoi
         bde=boj2*cst
         dt1=baa*gx-bab*dx
         dt2=baa*gy-bab*dy
         dt3=baa*gz-bab*dz
         dt7=bdd*dx-bde*gx
         dt8=bdd*dy-bde*gy
         dt9=bdd*dz-bde*gz
         dt4=-dt1-dt7
         dt5=-dt2-dt8
         dt6=-dt3-dt9
         
         kb = kb + 1
         b(i3+1,kb) = dt1
         b(i3+2,kb) = dt2
         b(i3+3,kb) = dt3
         b(j3+1,kb) = dt4
         b(j3+2,kb) = dt5
         b(j3+3,kb) = dt6
         b(k3+1,kb) = dt7
         b(k3+2,kb) = dt8
         b(k3+3,kb) = dt9
         iat(kb) = i3/3 + 1
         jat(kb) = j3/3 + 1
         kat(kb) = k3/3 + 1
         lat(kb) = 0
         avint(kb) = at*ratodeg
         if (ic > 0) then
            f(kb) = tk(ic)
         else
            f(kb) = 1.0
         end if
         ipt(kb) = 2
      30 continue
   end if  !  30 ithe=1,nthets
   
   if(ntrun <= 3) then
      nphi = nphih + nphia
      do 10 iphi=1,nphi
         if (iphi <= nphih) then
            i3=iph(iphi)
            j3=jph(iphi)
            k3=kph(iphi)
            l3=lph(iphi)
            k3=iabs(k3)
            l3=iabs(l3)
            ic = icph(iphi)
         else
            kphi = iphi - nphih
            i3 = ipa(kphi)
            j3 = jpa(kphi)
            k3 = iabs(kpa(kphi))
            l3 = iabs(lpa(kphi))
            ic = icpa(kphi)
         end if
         
         xij=x(i3+1)-x(j3+1)
         yij=x(i3+2)-x(j3+2)
         zij=x(i3+3)-x(j3+3)
         xkj=x(k3+1)-x(j3+1)
         ykj=x(k3+2)-x(j3+2)
         zkj=x(k3+3)-x(j3+3)
         xkl=x(k3+1)-x(l3+1)
         ykl=x(k3+2)-x(l3+2)
         zkl=x(k3+3)-x(l3+3)
         
         dx=yij*zkj-zij*ykj
         dy=zij*xkj-xij*zkj
         dz=xij*ykj-yij*xkj
         gx=zkj*ykl-ykj*zkl
         gy=xkj*zkl-zkj*xkl
         gz=ykj*xkl-xkj*ykl
         
         bi=dx*dx+dy*dy+dz*dz
         bk=gx*gx+gy*gy+gz*gz
         ct=dx*gx+dy*gy+dz*gz
         
         if(bi < 0.01 .or. bk < 0.01) goto 10
         boi2=1./bi
         boj2=1./bk
         bi  =sqrt(bi)
         bk  =sqrt(bk)
         z1  =1./bi
         z2  =1./bk
         bioj=bi*z2
         bjoi=bk*z1
         ct  =ct*z1*z2
         if (ct > 1. ) ct= 1.
         if (ct < (-1.)) ct=-1.
         ap  =acos(ct)
         
         s = xkj*(dz*gy - dy*gz)  +  ykj*(dx*gz - dz*gx) &
               +zkj*(dy*gx - dx*gy)
         
         if (s < 0.) ap=-ap
         pi = 2.*asin(1.0)
         ap=pi-ap
         
         
         !  NOW FOR THE MISERABLE GEOMETRIC TRANSFORMATION TO CARTESIANS.
         !  FOR NEARLY-PLANAR ANGLES, TAKE THE SHORT-CUT TO WILSON'S FORMULA.
         
         if(abs(ct) <= 0.99999) then
            
            st  =1.-ct*ct
            if (st < 1.e-6) st=1.e-6
            st  =sqrt(st)
            if (s > 0.) st=-st
            st  =-1.0/st
            boi2=boi2*st
            boj2=boj2*st
            ba=boi2*bioj
            bb=boi2*ct
            ga=boj2*bjoi
            gb=boj2*ct
            dc1=ba*gx-bb*dx
            dc2=ba*gy-bb*dy
            dc3=ba*gz-bb*dz
            dc4=ga*dx-gb*gx
            dc5=ga*dy-gb*gy
            dc6=ga*dz-gb*gz
            
            xik=x(i3+1)-x(k3+1)
            yik=x(i3+2)-x(k3+2)
            zik=x(i3+3)-x(k3+3)
            xlj=x(l3+1)-x(j3+1)
            ylj=x(l3+2)-x(j3+2)
            zlj=x(l3+3)-x(j3+3)
            
            dp1 =ykj*dc3 - zkj*dc2
            dp2 =zkj*dc1 - xkj*dc3
            dp3 =xkj*dc2 - ykj*dc1
            dp4 =yik*dc3 - zik*dc2 + ykl*dc6 - zkl*dc5
            dp5 =zik*dc1 - xik*dc3 + zkl*dc4 - xkl*dc6
            dp6 =xik*dc2 - yik*dc1 + xkl*dc5 - ykl*dc4
            dp7 =zij*dc2 - yij*dc3 + ylj*dc6 - zlj*dc5
            dp8 =xij*dc3 - zij*dc1 + zlj*dc4 - xlj*dc6
            dp9 =yij*dc1 - xij*dc2 + xlj*dc5 - ylj*dc4
            dp10=zkj*dc5 - ykj*dc6
            dp11=xkj*dc6 - zkj*dc4
            dp12=ykj*dc4 - xkj*dc5
            
         else
            
            bij=sqrt(xij*xij+yij*yij+zij*zij)
            bjk=sqrt(xkj*xkj+ykj*ykj+zkj*zkj)
            bkl=sqrt(xkl*xkl+ykl*ykl+zkl*zkl)
            ctijk=xij*xkj+yij*ykj+zij*zkj
            ctjkl=xkj*xkl+ykj*ykl+zkj*zkl
            
            ctijk=ctijk/(bij*bjk)
            ctjkl=ctjkl/(bjk*bkl)
            stijk=sqrt(abs(1.-ctijk**2))
            stjkl=sqrt(abs(1.-ctjkl**2))
            if(stijk < 1.e-4 .or. stjkl < 1.e-4) goto 10
            dna=1.0/(bi*bij*stijk)
            dnd=1.0/(bk*bkl*stjkl)
            dnb=ctijk/(bi*bjk*stijk)
            dnc=ctjkl/(bk*bjk*stjkl)
            
            dp1 =dx*dna
            dp2 =dy*dna
            dp3 =dz*dna
            dp10=gx*dnd
            dp11=gy*dnd
            dp12=gz*dnd
            difx=dx*dnb-gx*dnc
            dify=dy*dnb-gy*dnc
            difz=dz*dnb-gz*dnc
            dp4 =-dp1 +difx
            dp5 =-dp2 +dify
            dp6 =-dp3 +difz
            dp7 =-dp10-difx
            dp8 =-dp11-dify
            dp9 =-dp12-difz
            
         end if  ! (abs(ct) <= 0.99999)
         
         kb = kb + 1
         b(i3+1,kb) = dp1
         b(i3+2,kb) = dp2
         b(i3+3,kb) = dp3
         b(j3+1,kb) = dp4
         b(j3+2,kb) = dp5
         b(j3+3,kb) = dp6
         b(k3+1,kb) = dp7
         b(k3+2,kb) = dp8
         b(k3+3,kb) = dp9
         b(l3+1,kb) = dp10
         b(l3+2,kb) = dp11
         b(l3+3,kb) = dp12
         iat(kb) = i3/3 + 1
         jat(kb) = j3/3 + 1
         kat(kb) = k3/3 + 1
         lat(kb) = l3/3 + 1
         if (ap < 3.14159) then
            avint(kb) = ap*ratodeg
         else
            avint(kb) = ap*ratodeg - 360.
         end if
         if (ic > 0) then
            f(kb) = 0.5*pn(ic)*pn(ic)*pk(ic)
         else
            f(kb) = 1.0
         end if
         ipt(kb) = 3
         
      10 continue
   end if  !  10 iphi=1,nphi
   return
end subroutine bmat 
