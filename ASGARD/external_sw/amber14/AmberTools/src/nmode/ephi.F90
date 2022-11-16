
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
!+ [Enter a one-line description of subroutine ephi here]
subroutine ephi(nphi,ip,jp,kp,lp,icp,pk,pn,phase,cn1,cn2,ntypes, &
      cg,no,iac,ntb,x,f,dd, &
      ep,enbp,eelp,snb,see, &
      mphi,ecn,idiel,ndrv,nbel,natsys)
   implicit double precision (a-h,o-z)
   
   dimension nbel(*)
   dimension ip(nphi),jp(nphi),kp(nphi),lp(nphi),icp(nphi), &
         pk(*),pn(*),phase(*)

   dimension cn1(*),cn2(*),cg(*),no(*),iac(*),x(*),f(*),dd(*)

   dimension dc(6),t(6),dr(12),dtx(6,12),ddc(6,6),ddr(12,12)
   dimension xxij(3),xxkj(3),xxkl(3),xxil(3),fac2(2)
   
   ia(i) = i*(i-1)/2
   ibel(m3) = nbel(m3/3+1)
   data pi/3.1415926d0/
   data fac2/2.0d+00,6.0d+00/
   
   ep = 0.0d+00
   enbp = 0.0d+00
   eelp = 0.0d+00
   ecn = 0.0d+00
   scnb = snb
   scee = see
   kdiel = 1
   if(idiel == 0) kdiel = 2
   
   do 200 iphi = 1,nphi
      i3 = ip(iphi)
      j3 = jp(iphi)
      k3 = kp(iphi)
      l3 = lp(iphi)
      ic = icp(iphi)
      
      !     ----- l3.lt.0 for improper torsion angles
      !           k3.lt.0 for torsions with end atoms that should not
      !           have nonbonded terms calculated between them -----
      
      kdiv = 1
      if(l3 < 0) kdiv = 0
      if(k3 < 0) kdiv = 0
      k3 = iabs(k3)
      l3 = iabs(l3)
      
      !     ----- calculation of the dihedral -----
      
      do m = 1,3
         xj = x(j3+m)
         xk = x(k3+m)
         xxij(m) = x(i3+m)-xj
         xxkj(m) = xk-xj
         xxkl(m) = xk-x(l3+m)
      end do

      if(ntb /= 0) then
      
         !     ----- apply periodic boundary condition -----
      
         rij2 = xxij(1)**2+xxij(2)**2+xxij(3)**2
         rkj2 = xxkj(1)**2+xxkj(2)**2+xxkj(3)**2
         rkl2 = xxkl(1)**2+xxkl(2)**2+xxkl(3)**2
         call percon(rij2,xxij)
         call percon(rkj2,xxkj)
         call percon(rkl2,xxkl)
      end if

      xij = xxij(1)
      yij = xxij(2)
      zij = xxij(3)
      xkj = xxkj(1)
      ykj = xxkj(2)
      zkj = xxkj(3)
      xkl = xxkl(1)
      ykl = xxkl(2)
      zkl = xxkl(3)
      
      dx = yij*zkj-zij*ykj
      dy = zij*xkj-xij*zkj
      dz = xij*ykj-yij*xkj
      gx = zkj*ykl-ykj*zkl
      gy = xkj*zkl-zkj*xkl
      gz = ykj*xkl-xkj*ykl
      
      bi = dx*dx+dy*dy+dz*dz
      bk = gx*gx+gy*gy+gz*gz
      ct = dx*gx+dy*gy+dz*gz
      
      !     ----- branch if linear dihedral -----
      
      if(bk < 0.01.or.bi < 0.01) goto 301
      
      boi2 = 1./bi
      boj2 = 1./bk
      bi = sqrt(bi)
      bk = sqrt(bk)
      z1 = 1./bi
      z2 = 1./bk
      bioj = bi*z2
      bjoi = bk*z1
      ct = ct*z1*z2
      if(ct > 1. ) ct =  1.
      if(ct < (-1.)) ct = -1.
      ap = acos(ct)
      
      s = xkj*(dz*gy-dy*gz)+ykj*(dx*gz-dz*gx)+zkj*(dy*gx-dx*gy)
      
      if(s < 0.0d0) ap = -ap
      ap = pi-ap
      cphi = cos(ap)
      sphi = sin(ap)
      
      !     ----- ct value above is actually -cosphi; here we change
      !           its sign -----
      
      ct = -ct
      
      !     ----- calculate the energy and derivatives -----
      
      
      !     ----- get df = first der. of potential w/respect to cosphi; and
      !           ddf = second der. of potntial w/respect to cosphi -----
      
      !           the torsional potential is assumed to have the form:
      !            e = pk(ic) * (1.0+phase*cos(pn(ic)*phi)
      !            where phase = 1.0 or -1.0, and pn = 1,2,3,4, or 6
      
      !     ----- energy terms for dihedrals are expressed in terms of
      !           cosphi in order to eliminate problems for planar angles ----
      
      iper = int(dabs(pn(ic))+0.0001)
      if(iper == int(dabs(pn(ic))-0.0001)) goto 400
      if(iper > 6 .or. iper < 1) goto 400
      goto (30,35,40,45,50,55),iper
      
      !     ----- here for iper = 1 -----
      
      30 e = ct
      df = 1.0
      ddf = 0.0
      goto 60
      
      !     ----- here for iper = 2 -----
      
      35 e = 2.0*ct*ct-1.0
      df = 4.0*ct
      ddf = 4.0
      goto 60
      
      !     ----- here for iper = 3 -----
      
      40 ct2 = ct*ct
      e = ct*(4.0*ct2-3.0)
      df  = 12.0*ct2-3.0
      ddf = 24.0*ct
      goto 60
      
      !     ----- here for iper = 4 -----
      
      45 ct2 = ct*ct
      e = 1.0+ct2*8.0*(ct2-1.0)
      df = 16.0*ct*(ct2+ct2-1.0)
      ddf = 16.0*(6.0*ct2-1.0)
      goto 60
      
      !     ----- here for iper = 5 -----
      
      50 ct2 = ct*ct
      e = 1.0+ 16.*ct2*ct2*ct - 20.*ct2*ct + 5.*ct
      df = 80.*ct2*ct2 - 60.*ct2 + 5.
      ddf = 320.*ct2*ct - 120.*ct
      goto 60
      !     ----- here for iper = 6 -----
      
      55 ct2 = ct*ct
      e = ct2*(ct2*(ct2*32.0-48.0)+18.0)-1.0
      df = ct*(ct2*(ct2*192.0-192.0)+36.0)
      ddf = ct2*(ct2*960.0-576.0)+36.0
      
      60 arg = pk(ic)
      if(phase(ic) == 0.0) goto 65
      arg = -arg
      
      !     ----- if phase angle is other than 0 or pi then assume this angle
      !           has the old energy form e = pk*( 1.0+cos(pn*phi-phase) -----
      
      if(dabs(3.14159d0-phase(ic)) < 0.01) goto 65
      arg = -arg
      apn = abs(pn(ic))
      argum = apn*ap-phase(ic)
      e = cos(argum)
      dfp = -apn*sin(argum)
      ddfp = -apn*apn*e
      df = -dfp/sphi
      ddf = (ddfp-cphi*dfp/sphi)/(sphi*sphi)
      write(6,402) iphi,i3,j3,k3,l3
      
      65 e = pk(ic)+arg*e
      df = df*arg
      ddf = ddf*arg
      
      ep = ep+e
      if(iphi > mphi) ecn = ecn+e
      
      !     ----- skip the derivative sections if only energy needed -----
      
      if(ndrv <= 0) goto 301
      
      !     ----- now do torsional first and second derivatives -----
      
      !     ----- first, set up t array -----
      
      t(1) = dx
      t(2) = dy
      t(3) = dz
      t(4) = -gx
      t(5) = -gy
      t(6) = -gz
      
      !     ----- now, set up array dc = first der. of cosphi w/respect
      !           to the cartesian differences t -----
      
      z11 = z1*z1
      z12 = z1*z2
      z22 = z2*z2
      do i = 1,3
         dc(i) = t(i+3)*z12-cphi*t(i)*z11
         dc(i+3) = t(i)*z12-cphi*t(i+3)*z22
      end do
      
      !     ----- subroutine difang will now create array ddc which is second
      !           derivative of cosphi with respect to the t s -----
      
      if(ndrv >= 2) call difang(cphi,t,dc,ddc)
      
      !     ----- now set up array s, given on page 118 of cff book -----
      
      s1 = xij
      s2 = yij
      s3 = zij
      s4 = xkj
      s5 = ykj
      s6 = zkj
      s7 = -xkj
      s8 = -ykj
      s9 = -zkj
      s10 = -xkl
      s11 = -ykl
      s12 = -zkl
      
      !     ----- set up dtx(i,j) = derivative of t(i) w/respect to x(j)
      !           see p. 120 of cff book -----
      
      do i = 1,6
         do j = 1,12
            dtx(i,j) = 0.0
         end do
      end do
      dtx(1,2) = s6
      dtx(1,3) = -s5
      dtx(1,5) = s3-s6
      dtx(1,6) = s5-s2
      dtx(1,8) = -s3
      dtx(1,9) = s2
      dtx(2,1) = -s6
      dtx(2,3) = s4
      dtx(2,4) = s6-s3
      dtx(2,6) = s1-s4
      dtx(2,7) = s3
      dtx(2,9) = -s1
      dtx(3,1) = s5
      dtx(3,2) = -s4
      dtx(3,4) = s2-s5
      dtx(3,5) = s4-s1
      dtx(3,7) = -s2
      dtx(3,8) = s1
      dtx(4,5) = s12
      dtx(4,6) = -s11
      dtx(4,8) = s9-s12
      dtx(4,9) = s11-s8
      dtx(4,11) = -s9
      dtx(4,12) = s8
      dtx(5,4) = -s12
      dtx(5,6) = s10
      dtx(5,7) = s12-s9
      dtx(5,9) = s7-s10
      dtx(5,10) = s9
      dtx(5,12) = -s7
      dtx(6,4) = s11
      dtx(6,5) = -s10
      dtx(6,7) = s8-s11
      dtx(6,8) = s10-s7
      dtx(6,10) = -s8
      dtx(6,11) = s7
      
      !     ----- set up dr array, containing -first derivative of cosphi with
      !           respect to cartesians -----
      
      do i = 1,12
         dum = 0.0d0
         do j = 1,6
            dum = dum+dc(j)*dtx(j,i)
         end do
         dr(i) = -dum
      end do
      
      !     ----- update the force array -----
      
      f(i3+1) = f(i3+1)+df*dr(1)
      f(i3+2) = f(i3+2)+df*dr(2)
      f(i3+3) = f(i3+3)+df*dr(3)
      f(j3+1) = f(j3+1)+df*dr(4)
      f(j3+2) = f(j3+2)+df*dr(5)
      f(j3+3) = f(j3+3)+df*dr(6)
      f(k3+1) = f(k3+1)+df*dr(7)
      f(k3+2) = f(k3+2)+df*dr(8)
      f(k3+3) = f(k3+3)+df*dr(9)
      f(l3+1) = f(l3+1)+df*dr(10)
      f(l3+2) = f(l3+2)+df*dr(11)
      f(l3+3) = f(l3+3)+df*dr(12)
      
      if(ndrv <= 1) goto 301
      
      !     ----- now set up the ddr array = second der. of cosphi w/respect
      !           to cartesians; first we take the first term of last formula
      !           on p. 113 of cff book -----
      
      do i = 1,12
         do j = i,12
            ddr(i,j) = 0.0
            do k = 1,6
               if(dtx(k,i) == 0.0) continue
               do l = 1,6
                  ddr(i,j) = ddr(i,j)+ddc(k,l)*dtx(k,i)*dtx(l,j)
               end do
            end do
         end do
      end do
      
      !     ----- now do the second term of this equation -----
      
      ddr(2,9) = ddr(2,9)+dc(1)
      ddr(3,5) = ddr(3,5)+dc(1)
      ddr(6,8) = ddr(6,8)+dc(1)
      ddr(2,6) = ddr(2,6)-dc(1)
      ddr(3,8) = ddr(3,8)-dc(1)
      ddr(5,9) = ddr(5,9)-dc(1)
      ddr(1,6) = ddr(1,6)+dc(2)
      ddr(3,7) = ddr(3,7)+dc(2)
      ddr(4,9) = ddr(4,9)+dc(2)
      ddr(1,9) = ddr(1,9)-dc(2)
      ddr(3,4) = ddr(3,4)-dc(2)
      ddr(6,7) = ddr(6,7)-dc(2)
      ddr(1,8) = ddr(1,8)+dc(3)
      ddr(2,4) = ddr(2,4)+dc(3)
      ddr(5,7) = ddr(5,7)+dc(3)
      ddr(1,5) = ddr(1,5)-dc(3)
      ddr(2,7) = ddr(2,7)-dc(3)
      ddr(4,8) = ddr(4,8)-dc(3)
      ddr(5,12) = ddr(5,12)+dc(4)
      ddr(6,8) = ddr(6,8)+dc(4)
      ddr(9,11) = ddr(9,11)+dc(4)
      ddr(5,9) = ddr(5,9)-dc(4)
      ddr(6,11) = ddr(6,11)-dc(4)
      ddr(8,12) = ddr(8,12)-dc(4)
      ddr(4,9) = ddr(4,9)+dc(5)
      ddr(6,10) = ddr(6,10)+dc(5)
      ddr(7,12) = ddr(7,12)+dc(5)
      ddr(4,12) = ddr(4,12)-dc(5)
      ddr(6,7) = ddr(6,7)-dc(5)
      ddr(9,10) = ddr(9,10)-dc(5)
      ddr(4,11) = ddr(4,11)+dc(6)
      ddr(5,7) = ddr(5,7)+dc(6)
      ddr(8,10) = ddr(8,10)+dc(6)
      ddr(4,8) = ddr(4,8)-dc(6)
      ddr(5,10) = ddr(5,10)-dc(6)
      ddr(7,11) = ddr(7,11)-dc(6)
      
      !     ----- NOW FORM THE SECOND DERIVATIVE MATRIX -----
      
      do 260 i = 1,12
         ist = ibel(l3)
         if(i <= 9) ist = ibel(k3)
         if(i <= 6) ist = ibel(j3)
         if(i <= 3) ist = ibel(i3)
         if (ist /= -1) then
            iof = mod(i,3)
            if(iof == 0) iof = 3
            inew = ist+iof
            
            do 255 j = i,12
               jst = ibel(l3)
               if(j <= 9) jst = ibel(k3)
               if(j <= 6) jst = ibel(j3)
               if(j <= 3) jst = ibel(i3)
               if (jst /= -1) then
                  jof = mod(j,3)
                  if(jof == 0) jof = 3
                  jnew = jst+jof
                  
                  istore = ia(jnew)+inew
                  if(jnew < inew) istore = ia(inew)+jnew
                  dd(istore) = dd(istore)+ddf*dr(i)*dr(j)+df*ddr(i,j)
               end if
            255 continue
         end if
      260 continue
      
      !     ----- end of torsional second derivative calculations -----
      
      301 if(pn(ic) < 0) goto 200
      
      !     ----- compute nonbonded interactions -----
      
      if(kdiv == 0) goto 200
      
      ril2 = 0.0d0
      do 510 m = 1,3
         xxil(m) = x(i3+m)-x(l3+m)
         ril2 = ril2+xxil(m)**2
      510 continue
      if(ntb == 0) goto 520
      
      !     ----- apply periodic boundary condition -----
      
      call percon(ril2,xxil)
      520 continue
      
      xd = xxil(1)
      yd = xxil(2)
      zd = xxil(3)
      ii = (i3+3)/3
      jj = (l3+3)/3
      ia1 = iac(ii)
      ia2 = iac(jj)
      index = ntypes*(ia1-1)+ia2
      ic = no(index)
      
      !     ----- calculate the 14-eel energy -----
      
      r2 = 1.0d0/ril2
      if (idiel == 1) then
         
         !     ----- constant dielectric -----
         
         r1 = sqrt(r2)
         g = cg(ii)*cg(jj)*r1
         eelp = eelp+g
         df2 = -g/scee
         
         !     ----- distance dependent dielectric -----
         
      else if (idiel == 0) then
         g = cg(ii)*cg(jj)*r2
         eelp = eelp+g
         df2 = -(g+g)/scee
         
         !     ----- sigmoidal dielectric,  d= 78., s=0.3 hardwired
         
      else if (idiel == 2) then
         r1 = 1./sqrt(r2)
         sr = 0.3*r1
         expsr = exp(-sr)
         srsq = sr*sr
         epsr = r1*(79.0 - 39.*expsr*(srsq + sr + sr + 2))
         g = cg(ii)*cg(jj)/epsr
         df2 = -g*(1.0 + 1.053*expsr/(r2*r2*epsr))/scee
         eelp = eelp+g
#ifdef debug
         write(6,*) ii,jj,g,eelp
#endif
      end if

      
      !     ----- calculate regular 6-12 vdw for 1-4 h-bonds -----
      
      if(ic < 0) then
         ibig = max0(ia1,ia2)
         isml = min0(ia1,ia2)
         ic = ibig*(ibig-1)/2+isml
      end if
      r6 = r2*r2*r2
      r12 = r6*r6
      f1 = cn1(ic)*r12
      f2 = cn2(ic)*r6
      e = f1-f2
      enbp = enbp+e
      
      !     ----- scale the force appropriately -----
      
      df1 = (-12.*f1+6.*f2)/scnb
      
      !     ----- skip if derivatives not needed -----
      
      if(ndrv <= 0) goto 200
      
      !     ----- separate scale factor for nb and ee -----
      
      df = (df1+df2)*r2
      xa = xd*df
      ya = yd*df
      za = zd*df
      f(i3+1) = f(i3+1)-xa
      f(i3+2) = f(i3+2)-ya
      f(i3+3) = f(i3+3)-za
      f(l3+1) = f(l3+1)+xa
      f(l3+2) = f(l3+2)+ya
      f(l3+3) = f(l3+3)+za
      
      !     ----- calculate the second derivative if needed -----
      
      if(ndrv <= 1) goto 200
      
      ril = sqrt(ril2)
      dv1 = df*ril
      if (idiel < 2) then
         dv2 = r2*(156.0d0*f1-42.0d0*f2)/scnb + fac2(kdiel)*g*r2/scee
      else
         dv2 = r2*(156.0d0*f1-42.0d0*f2)/scnb + &
               g*(2.*r2 + 2.217618*expsr*expsr/(r6*epsr*epsr) &
               +0.3159*expsr*(ril**3)/epsr)/scee
      end if
      call difbon(dd, dv1,dv2,xxil,ril,i3,l3,nbel)
      
   200 continue
   
   enbp = enbp/scnb
   eelp = eelp/scee
   return
   400 write(6,402) iphi,i3,j3,k3,l3,iper, phase(ic)
   call mexit(6,1)
   402 format (' DIHEDRAL ANGLE ERROR: ',6i6,f10.4)
   403 format(' DIHEDRAL ',4i5,'  HAS ENERGY = ',f10.5)
end subroutine ephi 
