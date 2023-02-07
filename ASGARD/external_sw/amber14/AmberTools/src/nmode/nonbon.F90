
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
!+ [Enter a one-line description of subroutine nonbon here]
subroutine nonbon(natom,npair,iar1,iar2,iac,ico,x,f,h, &
      cn1,cn2,asol,bsol,hbcut,cg,xchrg,enb,ehb,eel, &
      idiel,ntypes,ndrv,nbel,natsys)
   implicit double precision (a-h,o-z)
   logical sysi, sysj, vsb, fori, forj
   
   common/setbox/box(3),beta,boxh(3),boxhm,boxhm2,boxoh,boxoq, &
         cosb,cosb2,ntm,ntb
   
   dimension iar1(*),iar2(*),iac(*),ico(*),x(*),f(*),cg(*)
   dimension h(*)
   dimension xchrg(*),cn1(*),cn2(*),asol(*),bsol(*),hbcut(*)
   dimension xij(3),xijp(3),fac2(2),nbel(natom)
   data fac2/2.0d+00,6.0d+00/
   
   enb = 0.0d+00
   eel = 0.0d+00
   ehb = 0.0d+00
   lpair = 0
   kdiel = 1
   if(idiel == 0) kdiel = 2
   
   !     ----- transform cartesian to oblique if necessary -----
   
   if(ntm /= 0) call traco(natom,0,x,beta,1)
   
   do 100 i = 1,natom
      npr = iar1(i)
      if (npr /= 0) then
         i3 = 3*i-3
         cgi = cg(i)
         iaci = iac(i)
         jaci = ntypes*(iaci-1)
         do 50 jj = 1,npr
            lpair = lpair+1
            j = iar2(lpair)
            j3 = 3*j-3
            do 10 m = 1,3
               xij(m) = x(i3+m)-x(j3+m)
            10 continue
            s = xij(1)**2 + xij(2)**2 + xij(3)**2
            if (ntb /= 0) then
               do 20 m = 1,3
                  if(xij(m) >= boxh(m)) then
                     xij(m) = xij(m)-box(m)
                  else if(xij(m) < -boxh(m)) then
                     xij(m) = xij(m)+box(m)
                  end if
               20 continue
               s = xij(1)**2+xij(2)**2+xij(3)**2
               if(ntb < 0) then
                  s = s+dmin1(0.0d0,boxoq-dabs(xij(1))-dabs(xij(2)) &
                        -dabs(xij(3)))*box(1)
               else if(ntm /= 0) then
                  s = s+cosb2*xij(1)*xij(3)
               end if
            end if
            
            iacj = iac(j)
            index = jaci+iacj
            ic = ico(index)
            if(ntm /= 0) call traco(1,0,xij,beta,-1)
            
            !     ----- claculate the electrostaic energy -----
            
            r2 = 1.0e0/s
            if (idiel == 0) then
               
               !     ----- distance dependent dielctric -----
               
               g = cg(i)*cg(j)*r2
               eel = eel+g
               df2 = -(g+g)
               
               !     ----- constant dielectric -----
               
            else if (idiel == 1) then
               r1 = sqrt(r2)
               g = cg(i)*cg(j)*r1
               eel = eel+g
               df2 = -g
               
               !     ----- sigmoidal dielectric,  d= 78., s=0.3 hardwired
               
            else if (idiel == 2) then
               r1 = 1./sqrt(r2)
               sr = 0.3*r1
               expsr = exp(-sr)
               srsq = sr*sr
               epsr = r1*(79.0 - 39.*expsr*(srsq + sr + sr + 2))
               g = cg(i)*cg(j)/epsr
               df2 = -g*(1.0 + 1.053*expsr/(r2*r2*epsr))
               eel = eel+g
            end if
            
            !     ----- vdw energy ... branch to h-bond if ic .lt. 0 -----
            
            if (ic > 0) then
               r6 = r2*r2*r2
               r12 = r6*r6
               f1 = cn1(ic)*r12
               f2 = cn2(ic)*r6
               e = f1-f2
               enb = enb + e
               df1 = (-12.*f1 + 6.*f2)
               df = (df1+df2)*r2
               if (ndrv == 2) then
                  bl = sqrt(s)
                  dv1 = df*bl
                  if (idiel < 2) then
                     dv2 = r2*(156.0e0*f1-42.0e0*f2)+fac2(kdiel)*g*r2
                  else
                     dv2 = r2*(156.0e0*f1-42.0e0*f2) + &
                           g*(2.*r2 + 2.217618*expsr*expsr/(r6*epsr*epsr) &
                           +0.3159*expsr*(bl**3)/epsr)
                  end if
                  call difbon(h,dv1,dv2,xij,bl,i3,j3,nbel)
               end if
            else
               
               !     ----- h-bond pairs 10-12 potential -----
               
               ic = iabs(ic)
               df1 = 0.0e0
               !                      ***note: amber 3.0 no longer uses hbcut!
               !             if (s.le.hbcut(ic)) then
               r10 = r2**5
               f1 = asol(ic)*r10*r2
               f2 = bsol(ic)*r10
               ehb = ehb+f1-f2
               df1 = (-12.0e0*f1+10.0e0*f2)
               !             end if
               df = (df1+df2)*r2
               
               !     ----- calculate the second derivative if needed -----
               
               if (ndrv == 2) then
                  bl = sqrt(s)
                  dv1 = df*bl
                  if (idiel < 2) then
                     dv2 = r2*(156.0*f1-110.0*f2)+fac2(kdiel)*g*r2
                  else
                     dv2 = r2*(156.0*f1-110.0*f2) + &
                           g*(2.*r2 + 2.217618*expsr*expsr/(r2*r2*r2*epsr*epsr) &
                           +0.3159*expsr*(bl**3)/epsr)
                  end if
                  call difbon(h,dv1,dv2,xij,bl,i3,j3,nbel)
               end if
            end if  ! (ic > 0)
            
            !     ----- update the force array -----
            
            xa = xij(1)*df
            ya = xij(2)*df
            za = xij(3)*df
            f(i3+1) = f(i3+1)-xa
            f(i3+2) = f(i3+2)-ya
            f(i3+3) = f(i3+3)-za
            f(j3+1) = f(j3+1)+xa
            f(j3+2) = f(j3+2)+ya
            f(j3+3) = f(j3+3)+za
            
         50 continue
      end if  ! (ntb /= 0)
   100 continue
   
   !     ----- transform the oblique coordinates to cartesian if needed
   
   if (ntm /= 0) call traco(natom,0,x,beta,-1)
   return
end subroutine nonbon 
