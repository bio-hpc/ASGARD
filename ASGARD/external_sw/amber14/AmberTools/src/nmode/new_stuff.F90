
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine politr here]
subroutine politr(natom,nres,c,q,pol,ipair,jpair, &
      eold,enew,eper,epol,xwij,r2,ipres, &
      aveper,aveind,avetot,emtot,nstep, &
      iptsol,n14,ni14,iarx,jpw,imgslt, &
      iptatm,scr1,scr2)
   !************************************************************************
   !                              AMBER                                   **
   !                                                                      **
   !                 COPYRIGHT (c) 1986, 1991,1995, 1997                  **
   !               regents of the university of california                **
   !                       all rights reserved.                           **
   !                                                                      **
   !  this software provided pursuant to a license agreement containing   **
   !  restrictions on its disclosure, duplication, and use. this software **
   !  contains confidential and proprietary information, and may not be   **
   !  extracted or distributed, in whole or in part, for any purpose      **
   !  whatsoever, without the express written permission of the authors.  **
   !  this notice, and the associated author list, must be attached to    **
   !  all copies, or extracts, of this software. any additional           **
   !  restrictions set forth in the license agreement also apply to this  **
   !  software.                                                           **
   !************************************************************************
   implicit real*8 (a-h,o-z)
   logical iskip
   
   !     ----- routine to calculate the polarization energy of the the system -----
   
   !     on entry:   c = coordinates
   !                 q = charges
   !                 pol = polarizabilities
   !                 ipair, jpair = pairwise interactions
   
   
   !     on exit:    epol = polarizability energy
   !                 eper = field due to fixed charges.
   !                 eold = final polarization field
   !                 enew = "new" induced moments
   
   common/heavy/tmas,ave_thet
   common/setbox/box(3),beta,boxh(3),boxhm,boxhm2,boxoh,boxoq, &
         cosb,cosb2,ntm,ntb

   dimension c(3,*),eold(3,*),enew(3,*),eper(3,*), &
         q(*),pol(*),scr1(*),scr2(*)
   dimension ipair(*),jpair(*),xwij(*),r2(*),ipres(*)
   dimension n14(*),ni14(15,*),iarx(*),jpw(*)

   zero = 0.0d0
   tol  = 0.01d0
   maxit = 20d0
   iskip = .false.

   !     ----- evaluate the initial field at each atoms -----

   do 20 i = 1,natom
      do 20 j = 1,3
         eper(j,i) = zero
   20 continue

   call efield(natom,c,eper,q,ipair,jpair,xwij,r2, &
         n14,ni14,iarx,jpw,imgslt,iptatm)

   do 120 i = 1,natom
      do 120 j =1,3
         eold(j,i) = eper(j,i)
   120 continue
   iter = 0
   
   100 continue
   iter = iter + 1

   !   ----- evaluate the total field and the rms deviation -----

   rms = zero
   call indip(natom,c,eold,enew,eper,q,pol,ipair,jpair, &
         rms,xwij,r2,scr1,scr2, &
         n14,ni14,iarx,jpw,imgslt,iptatm)
   write(6,1000)iter,rms
   1000 format(t2,'at iteration ',i2,' the rms =',f15.8)

   !   ------ eold and enew contain the same field at self consistency

   if(rms <= tol .or. iter == maxit) goto 140
   goto 100
   140 continue

   !     ----- calculate the induced moments for the energy calc

   do i = 1,natom
      enew(1,i) = eold(1,i)*pol(i)
      enew(2,i) = eold(2,i)*pol(i)
      enew(3,i) = eold(3,i)*pol(i)
   end do


   nv = 0
   pper = 0.0d0
   pind = 0.0d0
   ptot = 0.0d0
   emx = 0.0d0
   emy = 0.0d0
   emz = 0.0d0
   tot_thet = 0.0d0
   pp_x = 0.0d0
   pp_y = 0.0d0
   pp_z = 0.0d0
   pi_x = 0.0d0
   pi_y = 0.0d0
   pi_z = 0.0d0
   do 11 j = 1,nres
      px = 0.0d0
      py = 0.0d0
      pz = 0.0d0
      pindx = 0.0d0
      pindy = 0.0d0
      pindz = 0.0d0
      llim = ipres(j)
      iulim = ipres(j+1)-1
      if(iulim -llim > 0 ) then
         do i = llim,iulim
            px = px + q(i)*c(1,i)
            py = py + q(i)*c(2,i)
            pz = pz + q(i)*c(3,i)
            pindx = pindx + enew(1,i)
            pindy = pindy + enew(2,i)
            pindz = pindz + enew(3,i)
         end do
         emx = emx + px+pindx
         emy = emy + py+pindy
         emz = emz + pz+pindz

         !  sum over "solute"

         if(j <= iptsol) then
            pp_x = pp_x + px
            pp_y = pp_y + py
            pp_z = pp_z + pz
            pi_x = pi_x + pindx
            pi_y = pi_y + pindy
            pi_z = pi_z + pindz
         end if

         !  sum over "solvent"

         if(j > iptsol) then
            nv = nv + 1
            pi_mag = sqrt(pindx**2 + pindy**2 + pindz**2)
            pp_mag = sqrt(px**2 + py**2 + pz**2)
            pt_mag = sqrt( (px+pindx)**2 + (py+pindy)**2 &
                  + (pz+pindz)**2)

            pper = pper + pp_mag
            pind = pind + pi_mag
            ptot = ptot + pt_mag

            xtheta1 = (px*pindx + py*pindy + pz*pindz)
            if(pi_mag > 1.d-10.and.pp_mag > 1.d-10) then
               xtheta = (xtheta1)/(pi_mag*pp_mag)
               if(xtheta > 1.0d0)xtheta=1.0d0
               if(xtheta < -1.0d0)xtheta=-1.0d0
               theta = acos(xtheta)*180.0d0/3.14159d0
               tot_thet = tot_thet + theta
            end if
         end if

      end if  ! (iulim -llim > 0 )
   11 continue
   term1 = 4.8d0/18.2223d0
   emsq = ((emx*term1)**2 + (emy*term1)**2 + (emz*term1)**2)
   emtot = emtot + emsq/nres
   if(nv > 0) then
      ave_thet = tot_thet/(nv)
      aveper = (pper / nv) * term1
      aveind = (pind / nv) * term1
      avetot = (ptot / nv) * term1
   end if

   epol1 = 0.0d0
   do i = 1,natom
      do j = 1,3
         epol1 = epol1 + enew(j,i)*eper(j,i)
      end do
   end do
   epol = -epol1/2.0d0

   iskip = .true.
   return
end subroutine politr 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine efield here]
subroutine efield(natom,c,eper,q,ipol_pair,jpol_pair,xwij,r2, &
      n14,ni14,iarx,jpw,imgslt,iptatm)
   !************************************************************************
   !                            amber                                     **
   !                                                                      **
   ! copyright (c) 1986, 1991 regents of the university of california     **
   !                       all rights reserved.                           **
   !                                                                      **
   !  this software provided pursuant to a license agreement containing   **
   !  restrictions on its disclosure, duplication, and use. this software **
   !  contains confidential and proprietary information, and may not be   **
   !  extracted or distributed, in whole or in part, for any purpose      **
   !  whatsoever, without the express written permission of the authors.  **
   !  this notice, and the associated author list, must be attached to    **
   !  all copies, or extracts, of this software. any additional           **
   !  restrictions set forth in the license agreement also apply to this  **
   !  software.                                                           **
   !************************************************************************
   implicit real*8 (a-h,o-z)
   logical trnoct,sltimg,obliq
   
   common/polcon/tol
   common/setbox/box(3),beta,boxh(3),boxhm,boxhm2,boxoh,boxoq, &
         cosb,cosb2,ntm,ntb
   common/nbterm/cut,scnb,scee,idum,iddd,dielc,nbuck,numpk,nbit
   
   dimension c(3,*),eper(3,*),q(*),ipol_pair(*),jpol_pair(*), &
         xwij(3,*),r2(*)
   dimension n14(*),ni14(15,*),iarx(*),jpw(*)
   
   zero = 0.0d0
   one = 1.0d0
   trnoct = (ntb < 0)
   obliq =  (ntm /= 0)
   sltimg = (imgslt <= 0)
   ipair = 0
   do 400 i = 1,natom-1
      npack = 1
      dumx = zero
      dumy = zero
      dumz = zero
      qi = q(i)
      npr1 = (ipol_pair(i))
      if(npr1 > 0) then
         do 110 jn = 1,npr1
            iarx(jn) = jpol_pair(ipair+jn)
         110 continue
      end if
      npr2 = 0
      do 112 jn = 1,n14(i)
         iarx(npr1+jn) = ni14(jn,i)
         npr2 = npr2 + 1
      112 continue
      npr = npr1 + npr2
      if(npr <= 0) goto 400

      !     ----- loop over surrounding atoms -----

      do 140 jn = 1,npr
         j = iarx(jn)
         xwij(1,jn) = c(1,i) - c(1,j)
         xwij(2,jn) = c(2,i) - c(2,j)
         xwij(3,jn) = c(3,i) - c(3,j)
      140 continue

      !-----------periodic boundry conditions

      if(ntb == 0) then
         do jn = 1,npr
            r2(jn) = &
                  one/(xwij(1,jn)**2+xwij(2,jn)**2+xwij(3,jn)**2)
         end do
      else
         if (i <= iptatm) then
            if (sltimg) then
               do 200 jn = 1,npr
                  if(iarx(jn) > iptatm) then
                     if (abs(xwij(1,jn)) > boxh(1)) &
                           xwij(1,jn) = xwij(1,jn)-sign(box(1),xwij(1,jn))
                     if (abs(xwij(2,jn)) > boxh(2)) &
                           xwij(2,jn) = xwij(2,jn)-sign(box(2),xwij(2,jn))
                     if (abs(xwij(3,jn)) > boxh(3)) &
                           xwij(3,jn) = xwij(3,jn)-sign(box(3),xwij(3,jn))
                  end if
                  r2(jn) = &
                        one/(xwij(1,jn)**2+xwij(2,jn)**2+xwij(3,jn)**2)
               200 continue
            else
               do 201 jn = 1,npr
                  r2(jn) = &
                        one/(xwij(1,jn)**2+xwij(2,jn)**2+xwij(3,jn)**2)
               201 continue
            end if
         else
            do 210 jn = 1,npr
               if (abs(xwij(1,jn)) > boxh(1)) &
                     xwij(1,jn) = xwij(1,jn)-sign(box(1),xwij(1,jn))
               if (abs(xwij(2,jn)) > boxh(2)) &
                     xwij(2,jn) = xwij(2,jn)-sign(box(2),xwij(2,jn))
               if (abs(xwij(3,jn)) > boxh(3)) &
                     xwij(3,jn) = xwij(3,jn)-sign(box(3),xwij(3,jn))
               r2(jn) = &
                     one/(xwij(1,jn)**2+xwij(2,jn)**2+xwij(3,jn)**2)
            210 continue
         end if
         
         !     --- transform to cartesian if needed ---
         
      end if  !  200 jn = 1,npr
      300 continue
      do 310 jn = 1,npr
         j = iarx(jn)
         r3 = r2(jn)*sqrt(r2(jn))
         r3qi = qi*r3
         eper(1,j) = eper(1,j) - xwij(1,jn)*r3qi
         eper(2,j) = eper(2,j) - xwij(2,jn)*r3qi
         eper(3,j) = eper(3,j) - xwij(3,jn)*r3qi
         r3qj = q(j)*r3
         dumx = dumx + xwij(1,jn)*r3qj
         dumy = dumy + xwij(2,jn)*r3qj
         dumz = dumz + xwij(3,jn)*r3qj
      310 continue
      eper(1,i) = eper(1,i) + dumx
      eper(2,i) = eper(2,i) + dumy
      eper(3,i) = eper(3,i) + dumz
      ipair = ipair + npr1
   400 continue
   if(ntm /= 0)call traco(natom,0,c,beta,-1)
   return
end subroutine efield 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine indip here]
subroutine indip(natom,c,eold,enew,eper,q,pol,ipol_pair,jpol_pair &
      ,rms,xwij,r2,rsq,ri, &
      n14,ni14,iarx,jpw,imgslt,iptatm)
   !************************************************************************
   !                            AMBER                                     **
   !                                                                      **
   !                Copyright (c) 1986, 1991, 1995, 1997                  **
   !                Regents of the University of California               **
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
   implicit real*8 (a-h,o-z)
   
#define cforcevector c$dif no_recurrence
   
   
   logical obliq,sltimg,trnoct
   
   common/polcon/tol
   common/setbox/box(3),beta,boxh(3),boxhm,boxhm2,boxoh,boxoq, &
         cosb,cosb2,ntm,ntb
   common/nbterm/cut,scnb,scee,idum,iddd,dielc,nbuck,numpk,nbit
   
   dimension n14(*),ni14(15,*),iarx(*),jpw(*)
   dimension c(3,*),eold(3,*),enew(3,*),q(*),pol(*),eper(3,*)
   dimension ipol_pair(*),jpol_pair(*), &
         xwij(3,*),r2(*),ri(*),rsq(*)

   !     ----- calculate the induced polarization term -----

   trnoct = (ntm /= 0)
   obliq = (ntb < 0)
   sltimg = (imgslt == 1)
   zero = 0.0d0
   one = 1.0d0
   three = 3.0d0
   nat3 = 3*natom
   do 100 i = 1, natom
      do 100 j = 1,3
   100 enew(j,i) = eper(j,i)

   ipair = 0
   do 400 i = 1,natom-1
      poli = pol(i)
      sumx = zero
      sumy = zero
      sumz = zero
      npr1 = (ipol_pair(i))
      if(npr1 > 0) then
         do jn = 1,npr1
            iarx(jn) = jpol_pair(ipair+jn)
         end do
      end if
      npr2 = 0
      do jn = 1,n14(i)
         iarx(npr1+jn) = ni14(jn,i)
         npr2 = npr2 + 1
      end do
      npr = npr1 + npr2

      if (npr == 0) goto 400

      do jn = 1,npr
         j = iarx(jn)
         xwij(1,jn) = c(1,i) - c(1,j)
         xwij(2,jn) = c(2,i) - c(2,j)
         xwij(3,jn) = c(3,i) - c(3,j)
      end do

      if(ntb == 0) then
         do jn = 1,npr
            r2(jn) = &
                  one/(xwij(1,jn)**2+xwij(2,jn)**2+xwij(3,jn)**2)
         end do
      else
         if (i <= iptatm) then
            if (sltimg) then
               do 198 jn = 1,npr
                  if(iarx(jn) > iptatm) then
                     if (abs(xwij(1,jn)) > boxh(1)) &
                           xwij(1,jn) = xwij(1,jn)-sign(box(1),xwij(1,jn))
                     if (abs(xwij(2,jn)) > boxh(2)) &
                           xwij(2,jn) = xwij(2,jn)-sign(box(2),xwij(2,jn))
                     if (abs(xwij(3,jn)) > boxh(3)) &
                           xwij(3,jn) = xwij(3,jn)-sign(box(3),xwij(3,jn))
                  end if
               198 continue
            end if
         else
            do 210 jn = 1,npr
               if (abs(xwij(1,jn)) > boxh(1)) &
                     xwij(1,jn) = xwij(1,jn)-sign(box(1),xwij(1,jn))
               if (abs(xwij(2,jn)) > boxh(2)) &
                     xwij(2,jn) = xwij(2,jn)-sign(box(2),xwij(2,jn))
               if (abs(xwij(3,jn)) > boxh(3)) &
                     xwij(3,jn) = xwij(3,jn)-sign(box(3),xwij(3,jn))
            210 continue
         end if
         do 215 jn = 1,npr
            r2(jn) = one/(xwij(1,jn)**2+xwij(2,jn)**2+xwij(3,jn)**2)
         215 continue
         
         !     ----- transform to cartesian if needed ---
         
      end if
      300 continue
      do jn = 1,npr
         ri(jn)= 1.d0/r2(jn)
      end do
      do jn = 1,npr
         rsq(jn) = sqrt(r2(jn))
      end do
      !forcevector
      do 310 jn = 1,npr
         j = iarx(jn)
         polj = pol(j)
         xd3 = three*xwij(1,jn)
         yd3 = three*xwij(2,jn)
         zd3 = three*xwij(3,jn)
         
         xxd = xd3*xwij(1,jn)
         xyd = xd3*xwij(2,jn)
         xzd = xd3*xwij(3,jn)
         
         yyd = yd3*xwij(2,jn)
         yzd = yd3*xwij(3,jn)
         
         zzd = zd3*xwij(3,jn)
         
         r5 = rsq(jn)*r2(jn)*r2(jn)
         r5ai = poli*r5
         r5aj = polj*r5
         
         sumx = sumx+ &
               ((xxd-ri(jn))*eold(1,j)+xyd*eold(2,j)+xzd*eold(3,j))*r5aj
         sumy = sumy+ &
               (xyd*eold(1,j)+(yyd-ri(jn))*eold(2,j)+yzd*eold(3,j))*r5aj
         sumz = sumz+ &
               (xzd*eold(1,j)+yzd*eold(2,j)+(zzd-ri(jn))*eold(3,j))*r5aj
         
         dumxa = ((xxd-ri(jn))*eold(1,i)+xyd*eold(2,i)+xzd*eold(3,i)) &
               *r5ai
         dumya = (xyd*eold(1,i)+(yyd-ri(jn))*eold(2,i)+yzd*eold(3,i)) &
               *r5ai
         dumza = (xzd*eold(1,i)+yzd*eold(2,i)+(zzd-ri(jn))*eold(3,i)) &
               *r5ai
         
         enew(1,j) = enew(1,j) + dumxa
         enew(2,j) = enew(2,j) + dumya
         enew(3,j) = enew(3,j) + dumza
         
      310 continue
      enew(1,i) = enew(1,i) + sumx
      enew(2,i) = enew(2,i) + sumy
      enew(3,i) = enew(3,i) + sumz
      ipair = ipair + npr1
   400 continue
   
   !     ----- calculate the rms deviation -----
   
   dum = zero
   do 500 i = 1,natom
      duma = enew(1,i) - eold(1,i)
      dumb = enew(2,i) - eold(2,i)
      dumc = enew(3,i) - eold(3,i)
      dum = dum + duma*duma + dumb*dumb + dumc*dumc
      eold(1,i) = enew(1,i)
      eold(2,i) = enew(2,i)
      eold(3,i) = enew(3,i)
   500 continue
   rms = sqrt(dum/float(natom))
   return
end subroutine indip 

!************************************************************************
!                            amber                                     **
!                                                                      **
! copyright (c) 1986, 1991 regents of the university of california     **
!                       all rights reserved.                           **
!                                                                      **
!  this software provided pursuant to a license agreement containing   **
!  restrictions on its disclosure, duplication, and use. this software **
!  contains confidential and proprietary information, and may not be   **
!  extracted or distributed, in whole or in part, for any purpose      **
!  whatsoever, without the express written permission of the authors.  **
!  this notice, and the associated author list, must be attached to    **
!  all copies, or extracts, of this software. any additional           **
!  restrictions set forth in the license agreement also apply to this  **
!  software.                                                           **
!************************************************************************

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine pol2der here]
subroutine pol2der(natom,x,f,p,q,h,ipair,jpair,xrc,xij,r2,fw,vt, &
      n14,ni14,iarx,jpw,imgslt,iptatm,ndrv)
   implicit double precision (a-h,o-z)
   integer imat(6,6)
   logical trnoct,obliq,sltimg
   common/setbox/box(3),beta,boxh(3),boxhm,boxhm2,boxoh,boxoq, &
         cosb,cosb2,ntm,ntb
   common/nbterm/cut,scnb,scee,idiel,iddd,dielc,nbuck,numpk,nbit
   dimension x(3,*),f(3,*),p(3,*),q(*),ipair(*),jpair(*) &
         ,fw(3,*),xij(3,*),r2(*),vt(4),xrc(3,*) &
         ,n14(*),ni14(15,*),iarx(*),jpw(*) &
         ,h(*)

#define cforcevector c$dif no_recurrence

   !    calculate the 1st and 2nd derivatives for
   !    scf polariztion induced dipoles

   !    1st deriv code: jw caldwell and lx dang  (dec-jan 1992)

   !    2nd deriv code: jw caldwell (jan 1994)

   !     entry:
   !            x = coordinates (3,natom)
   !            p = scf polariz    "
   !            q = charges
   !            f = first deriv
   !            h = second deriv
   !         ipair,jpair pointers
   !         xij,r2,fw scratch

   !     exit:
   !            f = new 1st derivatives
   !            h = new 2nd deriv matrix
   !           vt = virials

   sltimg = (imgslt == 1)
   trnoct = (ntm /= 0)
   obliq = (ntb < 0)
   zero = 0.0d0
   one  = 1.0d0
   three = 3.0d0
   five = 5.0d0

   vt(1) = zero
   vt(2) = zero
   vt(3) = zero
   lpair = 0

   !  --- start here

   lim = natom - 1
   do 500 i = 1,lim
      npack = 1
      qi = q(i)
      npr1 = (ipair(i))
      if(npr1 > 0) then
         do 110 jn = 1,npr1
            iarx(jn) = jpair(lpair+jn)
         110 continue
      end if
      npr2 = 0
      do 112 jn = 1,n14(i)
         iarx(npr1+jn) = ni14(jn,i)
         npr2 = npr2 + 1
      112 continue
      npr = npr1 + npr2
      if (npr <= 0) goto 480
      do 140 jn = 1,npr
         j = iarx(jn)
         xij(1,jn) = x(1,i) - x(1,j)
         xij(2,jn) = x(2,i) - x(2,j)
         xij(3,jn) = x(3,i) - x(3,j)
      140 continue

      ! ------put in the p.b.c. -----------

      if(ntb == 0) then
         do 190 jn = 1,npr
            r2(jn) = &
                  one/(xij(1,jn)**2+xij(2,jn)**2+xij(3,jn)**2)
         190 continue
      else
         if (i <= iptatm) then
            if (sltimg) then
               do 200 jn = 1,npr
                  if(iarx(jn) > iptatm) then
                     if (abs(xij(1,jn)) > boxh(1)) &
                           xij(1,jn) = xij(1,jn)-sign(box(1),xij(1,jn))
                     if (abs(xij(2,jn)) > boxh(2)) &
                           xij(2,jn) = xij(2,jn)-sign(box(2),xij(2,jn))
                     if (abs(xij(3,jn)) > boxh(3)) &
                           xij(3,jn) = xij(3,jn)-sign(box(3),xij(3,jn))
                  end if
                  r2(jn) = &
                        one/(xij(1,jn)**2+xij(2,jn)**2+xij(3,jn)**2)
               200 continue
            else
               do 201 jn = 1,npr
                  r2(jn) = &
                        one/(xij(1,jn)**2+xij(2,jn)**2+xij(3,jn)**2)
               201 continue
            end if
         else
            do 210 jn = 1,npr
               if (abs(xij(1,jn)) > boxh(1)) &
                     xij(1,jn) = xij(1,jn)-sign(box(1),xij(1,jn))
               if (abs(xij(2,jn)) > boxh(2)) &
                     xij(2,jn) = xij(2,jn)-sign(box(2),xij(2,jn))
               if (abs(xij(3,jn)) > boxh(3)) &
                     xij(3,jn) = xij(3,jn)-sign(box(3),xij(3,jn))
               r2(jn) = &
                     one/(xij(1,jn)**2+xij(2,jn)**2+xij(3,jn)**2)
            210 continue
         end if


      end if

      do 390 jn = 1,npr
      390 r2(jn) = one/sqrt(r2(jn))

      !forcevector
      do 400 jn = 1,npr
         j = iarx(jn)
         rij   = r2(jn)
         rij2  = rij*rij
         rij4  = rij2*rij2
         rij5  = rij4*rij
         qj    = q(j)
         
         tx1j  =  ( &
               p(1,i) * (rij2-3.0d0*xij(1,jn)**2) &
               - 3.d0*p(2,i) * xij(1,jn)*xij(2,jn) &
               - 3.d0*p(3,i) * xij(1,jn)*xij(3,jn) &
               )



         ty1j  =  ( &
               p(2,i) * (rij2-3.0d0*xij(2,jn)**2) &
               - 3.d0*p(1,i) * xij(1,jn)*xij(2,jn) &
               - 3.d0*p(3,i) * xij(2,jn)*xij(3,jn) &
               )



         tz1j  =  ( &
               p(3,i) * (rij2-3.0d0*xij(3,jn)**2) &
               - 3.d0*p(1,i) * xij(1,jn)*xij(3,jn) &
               - 3.d0*p(2,i) * xij(2,jn)*xij(3,jn) &
               )



         tx1i  =  ( &
               p(1,j) * (rij2-3.0d0*xij(1,jn)**2) &
               - 3.d0*p(2,j) * xij(1,jn)*xij(2,jn) &
               - 3.d0*p(3,j) * xij(1,jn)*xij(3,jn) &
               )



         ty1i  =  ( &
               p(2,j) * (rij2-3.0d0*xij(2,jn)**2) &
               - 3.d0*p(1,j) * xij(1,jn)*xij(2,jn) &
               - 3.d0*p(3,j) * xij(2,jn)*xij(3,jn) &
               )



         tz1i  =  ( &
               p(3,j) * (rij2-3.0d0*xij(3,jn)**2) &
               - 3.d0*p(1,j) * xij(1,jn)*xij(3,jn) &
               - 3.d0*p(2,j) * xij(2,jn)*xij(3,jn) &
               )



         t2x1  = &
               p(1,j)*p(1,i)                            *xij(1,jn)**2 &
               + (p(1,j)*p(2,i) +  p(2,j)*p(1,i))*xij(2,jn)*xij(1,jn) &
               + (p(1,j)*p(3,i) +  p(3,j)*p(1,i))*xij(3,jn)*xij(1,jn)



         t2x2 =    ( &
               p(2,j)*p(2,i)                  *xij(2,jn)**2 &
               + (p(2,j)*p(3,i) + p(3,j)*p(2,i)) *xij(3,jn)*xij(2,jn) &
               +  p(3,j)*p(3,i)                  *xij(3,jn)**2 &
               )



         t2x3 = -(3.0d0*p(1,j)*p(1,i) &
               + p(2,j)*p(2,i) &
               + p(3,j)*p(3,i))



         t2x4 = &
               - ( (p(2,j)*p(1,i) + p(1,j)*p(2,i)) *xij(2,jn) &
               +   (p(3,j)*p(1,i) + p(3,i)*p(1,j)) *xij(3,jn) )



         t2y1 = &
               (p(1,j)*p(2,i) + p(2,j)*p(1,i))*xij(1,jn)*xij(2,jn) &
               +  p(2,j)*p(2,i)                           *xij(2,jn)**2 &
               + (p(2,j)*p(3,i) + p(3,j)*p(2,i))*xij(3,jn)*xij(2,jn)



         t2y2 = &
               (p(1,j)*p(1,i)                 *xij(1,jn)**2 &
               + (p(1,j)*p(3,i) + p(3,j)*p(1,i))*xij(1,jn)*xij(3,jn) &
               +  p(3,j)*p(3,i)                 *xij(3,jn)**2 &
               )



         t2y3 = -(3.0d0* p(2,j)*p(2,i) &
               + p(1,j)*p(1,i) &
               + p(3,j)*p(3,i))



         t2y4 = &
               -  ( (p(1,j)*p(2,i) + p(2,j)*p(1,i)) *xij(1,jn) &
               +(p(3,j)*p(2,i) + p(2,j)*p(3,i)) *xij(3,jn))



         t2z1 = &
               (p(1,j)*p(3,i) + p(3,j)*p(1,i))*xij(1,jn)*xij(3,jn) &
               + (p(3,j)*p(2,i) + p(2,j)*p(3,i))*xij(2,jn)*xij(3,jn) &
               +  p(3,j)*p(3,i)                           *xij(3,jn)**2



         t2z2 =        ( &
               p(1,j)*p(1,i)                 *xij(1,jn)**2 &
               + (p(1,j)*p(2,i) + p(2,j)*p(1,i))*xij(1,jn)*xij(2,jn) &
               +  p(2,j)*p(2,i)                 *xij(2,jn)**2 &
               )



         t2z3 = -(3.0d0*p(3,j)*p(3,i) &
               + p(1,j)*p(1,i) &
               + p(2,j)*p(2,i))



         t2z4 = &
               - ( (p(1,j)*p(3,i) + p(3,j)*p(1,i)) *xij(1,jn) &
               +   (p(2,j)*p(3,i) + p(3,j)*p(2,i)) *xij(2,jn))



         ! assemble the 1st and 2nd first derivative terms

         fx1i  =  (qi/rij5)*tx1i
         fy1i  =  (qi/rij5)*ty1i
         fz1i  =  (qi/rij5)*tz1i

         fx1j  =  (qj/rij5)*tx1j
         fy1j  =  (qj/rij5)*ty1j
         fz1j  =  (qj/rij5)*tz1j

         ! assemble the 3rd first derivative term

         fx2 = (3.0d0/rij5)* (5.0d0/rij2 * t2x1 * xij(1,jn) &
               + 5.0d0/rij2 * t2x2 * xij(1,jn) &
               + t2x3 * xij(1,jn) &
               + t2x4 )

         fy2 = (3.0d0/rij5)* (5.0d0/rij2 * t2y1 * xij(2,jn) &
               + 5.0d0/rij2 * t2y2 * xij(2,jn) &
               + t2y3 * xij(2,jn) &
               + t2y4)

         fz2 = (3.0d0/rij5)* (5.0d0/rij2 * t2z1 * xij(3,jn) &
               + 5.0d0/rij2 * t2z2 * xij(3,jn) &
               + t2z3 * xij(3,jn) &
               + t2z4)

         
         ! make the 1st derivative

         fw(1,jn) =  fx1i - fx1j + fx2
         fw(2,jn) =  fy1i - fy1j + fy2
         fw(3,jn) =  fz1i - fz1j + fz2

         f(1,j) = f(1,j) + fw(1,jn)
         f(2,j) = f(2,j) + fw(2,jn)
         f(3,j) = f(3,j) + fw(3,jn)

         f(1,i) = f(1,i) - fw(1,jn)
         f(2,i) = f(2,i) - fw(2,jn)
         f(3,i) = f(3,i) - fw(3,jn)

         if(ndrv < 2) goto 400

         !   now the make 2nd derivative terms

         dfx_dx = &
               - 5.d0*xij(1,jn)*fw(1,jn)/rij2 &

               - (qi/rij5)*(4.d0*p(1,j)*xij(1,jn) &
               + 3.d0*p(2,j)*xij(2,jn) &
               + 3.d0*p(3,j)*xij(3,jn) ) &

               + (qj/rij5)*(4.d0*p(1,i)*xij(1,jn) &
               + 3.d0*p(2,i)*xij(2,jn) &
               + 3.d0*p(3,i)*xij(3,jn) ) &

               + (3.d0/rij5)* ( &
               -  (10.d0*xij(1,jn)/rij4)*(t2x1+t2x2)*xij(1,jn)      &! #3&4
               +  t2x3                                              &! #5

               + (5.d0/rij2)*  ( &
               3.0d0*p(1,i)*p(1,j)            *xij(1,jn)**2     &! d3/dx
               +     2.d0*( p(1,i)*p(2,j)+p(2,i)*p(1,j) )             &!  "
               *xij(2,jn)*xij(1,jn)        &!  "
               +     2.d0*( p(1,i)*p(3,j)+p(3,i)*p(1,j) )             &!  "
               *xij(3,jn)*xij(1,jn)        &!  "
               +     t2x2)                                            &! d4/dx
               )

         !----------------------d2e/dxdy

         dfx_dy = &
               - 5.d0*xij(2,jn)*fx1i/rij2                      &! 1a

               +    (qi/rij5)*(2.d0*p(1,j)*xij(2,jn)           &! 1b
               - 3.d0*p(2,j)*xij(1,jn)) &

               + 5.d0*xij(2,jn)*fx1j/rij2                      &! 2a

               -    (qj/rij5)*(2.d0*p(1,i)*xij(2,jn)           &! 2b
               - 3.d0*p(2,i)*xij(1,jn)) &

               - 5.d0*xij(2,jn)*fx2/rij2                       &! 3.1

               + 3.d0/rij5 * ( &
               -(10.d0*xij(2,jn)/rij4)*(t2x1+t2x2)*xij(1,jn)                &! 3.3&5

               + (5.d0/rij2)*   (                                           &!
               ( p(1,i)*p(2,j) + p(2,i)*p(1,j) )          *xij(1,jn)**2  &! 3.2
               +  2.d0*p(2,i)*p(2,j)       *xij(2,jn)*xij(1,jn)     &! 3.4
               +( p(2,i)*p(3,j) + p(3,i)*p(2,j) )*xij(3,jn)*xij(1,jn)     &! "
               ) &

               - ( p(2,j)*p(1,i)+p(1,j)*p(2,j) )                            &! 3.6
               )

         !----------------------d2e/dxdz

         dfx_dz = &
               - 5.d0*xij(3,jn)*fx1i/rij2 &

               + (qi/rij5)*(2.d0*p(1,j)*xij(3,jn) - 3.d0*p(3,j)*xij(1,jn)) &

               + 5.d0*xij(3,jn)*fx1j/rij2 &

               - (qj/rij5)*(2.d0*p(1,i)*xij(3,jn) - 3.d0*p(3,i)*xij(1,jn)) &

               - 5.d0*xij(3,jn)*fx2/rij2 &

               + 3.d0/rij5 * ( &
               - (10.d0*xij(3,jn)/rij4)*(t2x1+t2x2)*xij(1,jn) &
               + (5.d0/rij2)*  ( &
               ( p(1,i)*p(3,j)+p(3,i)*p(1,j))         *xij(1,jn)**2 &
               +         2.d0*p(3,i)*p(3,j)      *xij(3,jn)*xij(1,jn) &
               +    (p(2,i)*p(3,j)+p(3,i)*p(2,j))*xij(2,jn)*xij(1,jn) &
               ) &

               - (p(3,j)*p(1,i)+p(1,j)*p(3,i)) &
               )

         ! --------------d2e/dydx

         dfy_dx = &
               -5.d0*xij(1,jn)*fy1i/rij2 &

               +(qi/rij5)*( 2.d0*p(2,j)*xij(1,jn) - 3.d0*p(1,j)*xij(2,jn)) &

               +5.d0*xij(1,jn)*fy1j/rij2 &

               -(qj/rij5)*( 2.d0*p(2,i)*xij(1,jn) - 3.d0*p(1,i)*xij(2,jn)) &

               - 5.d0*xij(1,jn)*fy2/rij2 &

               + 3.d0/rij5 * ( &

               - (10.d0*xij(1,jn)/rij4) *(t2y1+t2y2)*xij(2,jn) &

               +  (5.d0/rij2)*( &
               ( p(1,i)*p(2,j)+p(2,i)*p(1,j) )          *xij(2,jn)**2 &
               + 2.d0*p(1,i)*p(1,j)   *xij(1,jn)*xij(2,jn) &
               +( p(1,j)*p(3,i)+p(3,j)*p(1,i) )*xij(3,jn)*xij(2,jn) &
               ) &
               -  (p(1,j)*p(2,i)+p(2,j)*p(1,i)) &
               )

         ! --------------d2e/dy2

         dfy_dy = &
               -5.d0*xij(2,jn)*fy1i/rij2 &
               - (qi/rij5)*(4.d0*p(2,j)*xij(2,jn) &
               +  3.d0*p(1,j)*xij(1,jn) &
               +  3.d0*p(3,j)*xij(3,jn) ) &

               +5.d0*xij(2,jn)*fy1j/rij2 &
               + (qj/rij5)*(4.d0*p(2,i)*xij(2,jn) &
               +  3.d0*p(1,i)*xij(1,jn) &
               +  3.d0*p(3,i)*xij(3,jn) ) &

               -5.d0*xij(2,jn)*fy2/rij2 &

               +  3.d0/rij5* ( &
               -(10.d0*xij(2,jn)/rij4)*(t2y1+t2y2)*xij(2,jn) &
               +(5.d0/rij2)*( &

               3.0d0*p(2,i)*p(2,j)    *xij(2,jn)**2 &

               + 2.d0*(p(1,j)*p(2,i)+p(2,j)*p(1,i)) &
               *xij(1,jn)*xij(2,jn) &

               + 2.d0*(p(2,j)*p(3,i)+p(3,j)*p(2,i)) &
               *xij(3,jn)*xij(2,jn) &
               + t2y2) + t2y3 &
               )

         ! -------------d2e/dydz

         dfy_dz = &
               - 5.d0*xij(3,jn)*fy1i/rij2 &

               + (qi/rij5)*(2.d0*p(2,j)*xij(3,jn) - 3.d0*p(3,j)*xij(2,jn)) &

               + 5.d0*xij(3,jn)*fy1j/rij2 &

               - (qj/rij5)*(2.d0*p(2,i)*xij(3,jn) - 3.d0*p(3,i)*xij(2,jn)) &

               - (5.d0*xij(3,jn)/rij2)*fy2 &
               + 3.d0/rij5 * ( &
               - (10.d0*xij(3,jn)/rij4)*(t2y1+t2y2)*xij(2,jn) &
               +  (5.d0/rij2)*  ( &
               (p(2,i)*p(3,j)+p(3,i)*p(2,j))*xij(2,jn)**2 &
               + 2.d0*p(3,i)*p(3,j)     *xij(3,jn)*xij(2,jn) &
               +(p(1,i)*p(3,j)+p(3,i)*p(1,j))*xij(1,jn)*xij(2,jn) &
               ) &

               -  (p(2,j)*p(3,i)+p(2,i)*p(3,j)) &
               )

         !-------------- d2e/dzdx

         dfz_dx  = &
               - 5.d0*xij(1,jn)*fz1i/rij2 &
               +  (qi/rij5)*( 2.d0*p(3,j)*xij(1,jn) -3.d0*p(1,j)*xij(3,jn)) &

               + 5.d0*xij(1,jn)*fz1j/rij2 &
               -  (qj/rij5)*( 2.d0*p(3,i)*xij(1,jn) -3.d0*p(1,i)*xij(3,jn)) &

               -  5.d0*xij(1,jn)*fz2/rij2 &
               + (3.d0/rij5) *   ( &
               -(10.d0*xij(1,jn)/rij4)*(t2z1+t2z2)*xij(3,jn) &
               + (5.d0/rij2)*( &
               ( p(1,i)*p(3,j)+p(3,i)*p(1,j) )         *xij(3,jn)**2 &

               + 2.d0*p(1,i)*p(1,j)         *xij(1,jn)*xij(3,jn) &

               +( p(2,i)*p(1,j)+p(1,i)*p(2,j) ) *xij(2,jn)*xij(3,jn) &
               ) &
               - (p(3,j)*p(1,i)+p(1,j)*p(3,i)) &
               )

         ! ----------------d2e/dzdy

         dfz_dy = &
               -5.d0*xij(2,jn)*fz1i/rij2 &
               +  qi/rij5*(2.d0*p(3,j)*xij(2,jn) - 3.d0*p(2,j)*xij(3,jn)) &

               +5.d0*xij(2,jn)*fz1j/rij2 &
               -  qj/rij5*(2.d0*p(3,i)*xij(2,jn) - 3.d0*p(2,i)*xij(3,jn)) &

               -  5.d0*xij(2,jn)*fz2/rij2 &

               + 3.d0/rij5 * ( &
               - (10.d0*xij(2,jn)/rij4)*(t2z1+t2z2)*xij(3,jn) &
               + (5.d0/rij2)*  ( &

               ( p(2,i)*p(3,j)+p(3,i)*p(2,j) )          *xij(3,jn)**2 &
               +2.d0*p(2,j)*p(2,i)        *xij(2,jn)*xij(3,jn) &
               +( p(1,j)*p(2,i)+p(2,j)*p(1,i) )*xij(1,jn)*xij(3,jn) &
               ) &
               -  (p(2,j)*p(3,i)+p(3,j)*p(2,i)) &
               )

         !------------------d2e/dz2

         dfz_dz = &
               -5.d0*xij(3,jn)*fz1i/rij2 &
               -(qi/rij5)*(4.d0*p(3,j)*xij(3,jn) &
               + 3.d0*p(1,j)*xij(1,jn) &
               + 3.d0*p(2,j)*xij(2,jn)) &

               +5.d0*xij(3,jn)*fz1j/rij2 &

               +(qj/rij5)*(4.d0*p(3,i)*xij(3,jn) &
               + 3.d0*p(1,i)*xij(1,jn) &
               + 3.d0*p(2,i)*xij(2,jn)) &

               -5.d0*xij(3,jn)*fz2/rij2 &

               + 3.d0/rij5* ( &
               -  (10.d0*xij(3,jn)/rij4)*(t2z1+t2z2)*xij(3,jn) &
               +  (5.d0/rij2)*( &
               3.0d0*p(3,i)*p(3,j)                 *xij(3,jn)**2 &
               +   2.d0*(p(1,i)*p(3,j)+p(3,i)*p(1,j)) &
               *xij(1,jn)*xij(3,jn) &
               +   2.d0*(p(2,i)*p(3,j)+p(3,i)*p(2,j)) &
               *xij(2,jn)*xij(3,jn) &
               +  t2z2 )  + t2z3 &
               )


         ! calculate the indices to load 2nd deriv matrix

         i3 = 3*(i-1)
         j3 = 3*(j-1)
         call loadit(imat,i3,j3)

         !  d2e/xixi (1)

         index1 = imat(1,1)
         h(index1) = h(index1) + dfx_dx

         !  d2e/xixj (1)

         index2= imat(1,4)
         h(index2) = h(index2) - dfx_dx

         !  d2e/xjxj (1)

         index3= imat(4,4)
         h(index3) = h(index3) + dfx_dx

         !-------------------------
         !  d2e/yiyi (2)

         index1 = imat(2,2)
         h(index1) = h(index1) + dfy_dy

         !  d2e/yiyj (2)

         index2 = imat(2,5)
         h(index2) = h(index2) - dfy_dy

         !  d2e/yjyj (2)

         index3 = imat(5,5)
         h(index3) = h(index3) + dfy_dy

         !-------------------------
         !  d2e/zizi (3)

         index1 = imat(3,3)
         h(index1) = h(index1) + dfz_dz

         !  d2e/zizj (3)

         index2 = imat(3,6)
         h(index2) = h(index2) - dfz_dz

         !  d2e/zjzj (3)

         index3 = imat(6,6)
         h(index3) = h(index3) + dfz_dz

         !-------------------------
         ! d2e/xiyi (4)

         index1 = imat(1,2)
         h(index1) = h(index1) + dfx_dy

         ! d2e/xiyj (4)

         index2 = imat(1,5)
         h(index2) = h(index2) - dfx_dy

         ! d2e/xjyi (4)

         index3 = imat(2,4)
         h(index3) = h(index3) - dfx_dy

         ! d2e/xjyj (4)

         index4 = imat(4,5)
         h(index4) = h(index4) + dfx_dy

         !-------------------------
         ! d2e/xizi (5)

         index1 = imat(1,3)
         h(index1) = h(index1) + dfx_dz

         ! d2e/xizj (5)

         index2 = imat(1,6)
         h(index2) = h(index2) - dfx_dz

         ! d2e/xjzi (5)

         index3 = imat(3,4)
         h(index3) = h(index3) - dfx_dz

         ! d2e/xjzj (5)

         index4 = imat(4,6)
         h(index4) = h(index4) + dfx_dz

         !-------------------------
         ! d2e/yizi (6)

         index1 = imat(2,3)
         h(index1) = h(index1) + dfy_dz

         ! d2e/yizj (6)

         index2 = imat(2,6)
         h(index2) = h(index2) - dfy_dz

         ! d2e/yjzi  (6)

         index3 = imat(3,5)
         h(index3) = h(index3) - dfy_dz

         ! d2e/yjzj  (6)

         index4 = imat(5,6)
         h(index4) = h(index4) + dfy_dz

         !--------------
      400 continue

      !  calculate the contribution to the virial

      do 410 jn = 1,npr
         j = iarx(jn)
         vt(1) = vt(1) + fw(1,jn)*(xij(1,jn) - xrc(1,i)+xrc(1,j))
         vt(2) = vt(2) + fw(2,jn)*(xij(2,jn) - xrc(2,i)+xrc(2,j))
         vt(3) = vt(3) + fw(3,jn)*(xij(3,jn) - xrc(3,i)+xrc(3,j))
      410 continue
      480 continue
      lpair = lpair + npr1
   500 continue

   vt(1) = vt(1)*0.5d0
   vt(2) = vt(2)*0.5d0
   vt(3) = vt(3)*0.5d0

   !      call matout(h,3*natom)
   return
end subroutine pol2der 

!************************************************************************
!                            amber                                     **
!                                                                      **
! copyright (c) 1986, 1991 regents of the university of california     **
!                       all rights reserved.                           **
!                                                                      **
!  this software provided pursuant to a license agreement containing   **
!  restrictions on its disclosure, duplication, and use. this software **
!  contains confidential and proprietary information, and may not be   **
!  extracted or distributed, in whole or in part, for any purpose      **
!  whatsoever, without the express written permission of the authors.  **
!  this notice, and the associated author list, must be attached to    **
!  all copies, or extracts, of this software. any additional           **
!  restrictions set forth in the license agreement also apply to this  **
!  software.                                                           **
!************************************************************************


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine loadit here]
subroutine loadit(istore,i3,j3)
   dimension istore(6,6)
   logical sw
   
   !     ----- form the index array and merge df and ddf -----
   
   sw = .false.
   if(i3 > j3) sw = .true.
   do i = 1,6
      do j = 1,6
         istore(i,j) = 0
      end do
   end do
   if(sw) then
      item = j3
      j3 = i3
      i3 = item
   end if
   istore(1,1) = i3*(i3+1)/2+i3+1
   istore(1,2) = istore(1,1)+i3+1
   istore(2,2) = istore(1,2)+1
   istore(1,3) = istore(1,2)+i3+2
   istore(2,3) = istore(1,3)+1
   istore(3,3) = istore(2,3)+1
   istore(1,4) = j3*(j3+1)/2+i3+1
   istore(2,4) = istore(1,4)+1
   istore(3,4) = istore(2,4)+1
   istore(1,5) = istore(1,4)+j3+1
   istore(2,5) = istore(1,5)+1
   istore(3,5) = istore(2,5)+1
   istore(1,6) = istore(1,5)+j3+2
   istore(2,6) = istore(1,6)+1
   istore(3,6) = istore(2,6)+1
   istore(4,4) = istore(1,4)-i3+j3
   istore(4,5) = istore(4,4)+j3+1
   istore(5,5) = istore(4,5)+1
   istore(4,6) = istore(4,5)+j3+2
   istore(5,6) = istore(4,6)+1
   istore(6,6) = istore(5,6)+1
   if(sw) then

      ! yeah, i know this is ugly but so is the above.

      itemp = istore(1,1)
      istore(1,1) = istore(4,4)
      istore(4,4) = itemp

      itemp = istore(2,2)
      istore(2,2) = istore(5,5)
      istore(5,5) = itemp

      itemp = istore(3,3)
      istore(3,3) = istore(6,6)
      istore(6,6) = itemp

      itemp = istore(1,2)
      istore(1,2) = istore(4,5)
      istore(4,5) = itemp

      itemp = istore(1,3)
      istore(1,3) = istore(4,6)
      istore(4,6) = itemp

      itemp = istore(2,3)
      istore(2,3) = istore(5,6)
      istore(5,6) = itemp


      itemp = istore(1,5)
      istore(1,5) = istore(2,4)
      istore(2,4) = itemp
      itemp = istore(1,6)
      istore(1,6) = istore(3,4)
      istore(3,4) = itemp
      itemp = istore(2,6)
      istore(2,6) = istore(3,5)
      istore(3,5) = itemp

   end if

   return
end subroutine loadit 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine tripl here]
subroutine tripl(natom,iar1,iar2,ityp,itrp)
   implicit double precision(a-h,o-z)
   common/thrbod/iat1(5),iat2(5),acon(5),beta3b(5),gama3b(5) &
         ,i3b(5),n3b,nion

   !    routine to generate triplets of atoms

   !    note: the ions must be at the beginning of the
   !          atom list

   !    for an interaction i,j,k there must exist pairs
   !    i,j and i,k and j,k

   !    note that the order of the pair atom numbers will not
   !    necessarily be in ascending order due to the seperation
   !    of 6-12 and 10-12 lists


   dimension iar1(*),iar2(*),itrp(3,*),ityp(*)
   dimension ix1(5000),ix2(5000)
   
   !     iar1 = number of pairs for each atom (iar1(1,i) regular pairs
   !                                          (iar1(2,i) h-b pairs
   !     iar2 = second index of i,j
   !     itrp will contain the triplet( i(ion)  j,k (water oxygens))

   !     alph = three body interactions parameters

   jt_cnt = 0
   zero = 0.0d0
   one  = 1.0d0
   itot = 0
   jbig = 0
   do 110 jj = 1,n3b
      ioff = 0
      it_cnt = 0
      if(nion == 0) then
         limi = natom - 1
      else
         limi = nion
      end if
      
      ! outer loop over i (i-j-k triplet)
      
      lpack1 = 1
      lpack2 = 1
      
      do 100 i = 1,limi
         ipr = iar1(i)
         if(ipr == 0) goto 100
         
         do 210 kk = 1,ipr
            ix1(kk) = abs(iar2(kk+ioff))
         210 continue
         do 80 jn = 1,ipr
            j = ix1(jn)
            if(ityp(j) == iat2(jj))then
               jpr = iar1(j)
               npack = 1
               if(jpr == 0) goto 80
               
               ! now for i,j pairs (joff will be the number of pairs
               ! before the current j)
               
               joff = 0
               lpack2 = 1
               do 30 jo = 1,j-1
                  jtemp =  iar1(jo)
                  joff = joff + jtemp
                  lpack2 = lpack2 + jtemp + 1
               30 continue
               
               do 310 kj = 1,jpr
                  ix2(kj) = iar2(kj+joff)
               310 continue
               do 50 kn = 1,jpr
                  k = abs(ix2(kn))
                  if(ityp(k) == iat2(jj)) then
                     
                     !   we have a-b and b-c
                     !   now test if a-c exists
                     
                     do 48 ln = 1,ipr
                        l = ix1(ln)
                        if(l == k) then
                           it_cnt  = it_cnt + one
                           if(it_cnt > 1000) then
                              print *,'too many triplets max = 1000...calc =',it_cnt
                              call exit
                           end if
                           jt_cnt  = jt_cnt + 1
                           itrp(1,jt_cnt ) = i
                           itrp(2,jt_cnt ) = j
                           itrp(3,jt_cnt ) = k
                        end if
                     48 continue
                  end if
               50 continue
            end if
         80 continue
         ioff = ioff + ipr
      100 continue
      i3b(jj) = it_cnt
      itot = itot + it_cnt
   110 continue
   
   !  end of triplet generation
   
   write (6,'(t2,''total 3bods = '',i5)')itot

   999 continue
   return
end subroutine tripl 
!------------------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine ugly here]
subroutine ugly(n3b,itrp)
   implicit double precision(a-h,o-z)
   dimension itrp(3,*)
   common/thrbod/iat1(5),iat2(5),acon(5),beta3b(5),gama3b(5) &
         ,i3b(5),idum,jdum

   i3b(1) = 0
   if(n3b == -1)then

      !  specific for ion-ethylene

      itrp(1,1) = 1
      itrp(2,1) = 2
      itrp(3,1) = 5
      i3b(1) = 1

   else if(n3b == -6) then

      ! specific for ion-benzene

      itrp(1,1) = 1
      itrp(2,1) = 2
      itrp(3,1) = 4
      itrp(1,2) = 1
      itrp(2,2) = 4
      itrp(3,2) = 6
      itrp(1,3) = 1
      itrp(2,3) = 6
      itrp(3,3) = 8
      itrp(1,4) = 1
      itrp(2,4) = 8
      itrp(3,4) = 10
      itrp(1,5) = 1
      itrp(2,5) = 10
      itrp(3,5) = 12
      itrp(1,6) = 1
      itrp(2,6) = 12
      itrp(3,6) = 2
      i3b(1) = 6
   else
      write(6,'(t2,''n3b value not allowed: '',i3)')n3b
      call mexit(6,0)
   end if
   n3b = 1
   write (6,'(t2,''ugly 3bods = '',i5)')i3b(1)
   return
end subroutine ugly 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine threeb here]
subroutine threeb(x,f,itrp,e3bod,natom)
   implicit double precision(a-h,o-z)

   !     ----- routine to evaluate the energy and gradient due to
   !           the three body non-bonded interactions

   common/setbox/box(3),beta,boxh(3),boxhm,boxhm2,boxoh,boxoq, &
         cosb,cosb2,ntm,ntb

   common/thrbod/iat1(5),iat2(5),acon(5),beta3b(5),gama3b(5) &
         ,i3b(5),n3b,nion

   dimension x(3,*),f(3,*),itrp(3,*)

   !     ----- main loop over the ions -----

   icount = 0
   do 900 j = 1,abs(n3b)
      do 900 i = 1,i3b(j)

         icount = icount + 1
         ii = itrp(1,icount)
         jj = itrp(2,icount)
         kk = itrp(3,icount)

         x12 = x(1,ii)-x(1,jj)
         y12 = x(2,ii)-x(2,jj)
         z12 = x(3,ii)-x(3,jj)

         if(ntb > 0) then
            if(abs(x12) > boxh(1)) x12 = x12-sign(box(1),x12)
            if(abs(y12) > boxh(2)) y12 = y12-sign(box(2),y12)
            if(abs(z12) > boxh(3)) z12 = z12-sign(box(3),z12)
         end if
         r12 = sqrt(x12**2+y12**2+z12**2)
         xe12 = x12/r12
         ye12 = y12/r12
         ze12 = z12/r12

         !  indices are "backwards"

         x13 = x(1,ii)-x(1,kk)
         y13 = x(2,ii)-x(2,kk)
         z13 = x(3,ii)-x(3,kk)
         if(ntb > 0) then
            if(abs(x13) > boxh(1)) x13 = x13-sign(box(1),x13)
            if(abs(y13) > boxh(2)) y13 = y13-sign(box(2),y13)
            if(abs(z13) > boxh(3)) z13 = z13-sign(box(3),z13)
         end if
         r13 = sqrt(x13**2+y13**2+z13**2)
         xe13 = x13/r13
         ye13 = y13/r13
         ze13 = z13/r13

         ! distance between first and second water (atoms jj and kk)

         x23 = x(1,jj)-x(1,kk)
         y23 = x(2,jj)-x(2,kk)
         z23 = x(3,jj)-x(3,kk)
         if(ntb > 0) then
            if(abs(x23) > boxh(1)) x23 = x23-sign(box(1),x23)
            if(abs(y23) > boxh(2)) y23 = y23-sign(box(2),y23)
            if(abs(z23) > boxh(3)) z23 = z23-sign(box(3),z23)
         end if
         r23 = sqrt(x23**2+y23**2+z23**2)
         xe23 = x23/r23
         ye23 = y23/r23
         ze23 = z23/r23
         
         ! calculate and sum up the potential
         
         e12   = exp(-beta3b(j)*r12)
         e13   = exp(-beta3b(j)*r13)
         e23   = exp( -gama3b(j)*r23)
         e123  = acon(j)*e12*e13*e23
         e3bod = e3bod+e123

         ! calculate and sum up the gradient ----- ion

         term1 = beta3b(j)*(e123*xe12 + e123*xe13)
         term2 = beta3b(j)*(e123*ye12 + e123*ye13)
         term3 = beta3b(j)*(e123*ze12 + e123*ze13)
         f(1,ii) = f(1,ii) + term1
         f(2,ii) = f(2,ii) + term2
         f(3,ii) = f(3,ii) + term3

         ! calculate and sum up the gradient ----- first water

         term1= (-beta3b(j)*e123*xe12 + gama3b(j)*e123*xe23)
         term2= (-beta3b(j)*e123*ye12 + gama3b(j)*e123*ye23)
         term3= (-beta3b(j)*e123*ze12 + gama3b(j)*e123*ze23)
         f(1,jj) = f(1,jj) + term1
         f(2,jj) = f(2,jj) + term2
         f(3,jj) = f(3,jj) + term3

         ! calculate and sum up the gradient ----- second water

         term1= (-beta3b(j)*e123*xe13 - gama3b(j)*e123*xe23)
         term2= (-beta3b(j)*e123*ye13 - gama3b(j)*e123*ye23)
         term3= (-beta3b(j)*e123*ze13 - gama3b(j)*e123*ze23)
         f(1,kk) = f(1,kk) + term1
         f(2,kk) = f(2,kk) + term2
         f(3,kk) = f(3,kk) + term3
   900 continue
   
   return
end subroutine threeb 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine thr_2nd here]
subroutine thr_2nd(x,h,i,j,k,e,beta,gama)
   implicit double precision(a-h,o-z)
   dimension x(3,*),h(*)
   dimension index(6,6),jndex(6,6),kndex(6,6)

   ! -- 2nd deriviative over "3" body terms.


   xij = x(1,i)-x(1,j)
   yij = x(2,i)-x(2,j)
   zij = x(3,i)-x(3,j)

   xik = x(1,i)-x(1,k)
   yik = x(2,i)-x(2,k)
   zik = x(3,i)-x(3,k)

   xjk = (x(1,j)-x(1,k))
   yjk = (x(2,j)-x(2,k))
   zjk = (x(3,j)-x(3,k))

   rij = sqrt( xij**2 + yij**2 + zij**2)
   rik = sqrt( xik**2 + yik**2 + zik**2)
   rjk = sqrt( xjk**2 + yjk**2 + zjk**2)

   !----------- atom 1

   d2e_dx1_2 =  e *beta * ( beta* (xij/rij + xik/rik)**2 &
         +  (xij**2/rij**3  - 1.d0/rij + xik**2/rik**3 - 1.d0/rik ))

   d2e_dy1_2 =  e *beta * ( beta* (yij/rij + yik/rik)**2 &
         +  (yij**2/rij**3  - 1.d0/rij + yik**2/rik**3 - 1.d0/rik ))

   d2e_dz1_2 =  e *beta * ( beta* (zij/rij + zik/rik)**2 &
         +  (zij**2/rij**3  - 1.d0/rij + zik**2/rik**3 - 1.d0/rik ))

   !----------- atom 2

   d2e_dx2_2 =    e *((beta*xij/rij - gama*xjk/rjk)**2 &
         +    (beta*(-1.d0/rij +   xij**2/rij**3) &
         +     gama*(-1.d0/rjk +   xjk**2/rjk**3)))

   d2e_dy2_2 =    e *((beta*yij/rij - gama*yjk/rjk)**2 &
         +    (beta*(-1.d0/rij +   yij**2/rij**3) &
         +     gama*(-1.d0/rjk +   yjk**2/rjk**3)))

   d2e_dz2_2 =    e *((beta*zij/rij - gama*zjk/rjk)**2 &
         +    (beta*(-1.d0/rij +   zij**2/rij**3) &
         +     gama*(-1.d0/rjk +   zjk**2/rjk**3)))

   !----------- atom 3

   d2e_dx3_2 =    e *((beta*xik/rik + gama*xjk/rjk)**2 &
         +  (beta*  (-1.d0/rik +   xik**2/rik**3) &
         +   gama*  ( 1.d0/rjk -   xjk**2/rjk**3)))

   d2e_dy3_2 =    e *((beta*yik/rik + gama*yjk/rjk)**2 &
         +  (beta*  (-1.d0/rik +   yik**2/rik**3) &
         +   gama*  ( 1.d0/rjk -   yjk**2/rjk**3)))

   d2e_dz3_2 =    e *((beta*zik/rik + gama*zjk/rjk)**2 &
         +  (beta*  (-1.d0/rik +   zik**2/rik**3) &
         +   gama*  ( 1.d0/rjk -   zjk**2/rjk**3)))


   !----------- 1-2

   d2e_dx2_dx1 = -beta * e * ( &
         (xij/rij + xik/rik) * (beta*xij/rij - gama*xjk/rjk) &
         -1.d0/rij + xij**2/rij**3 )

   d2e_dy2_dx1 = -beta * e * ( &
         (xij/rij + xik/rik) * (beta*yij/rij - gama*yjk/rjk) &
         + xij*yij/rij**3 )

   d2e_dz2_dx1 = -beta * e * ( &
         (xij/rij + xik/rik) * (beta*zij/rij - gama*zjk/rjk) &
         + xij*zij/rij**3 )

   d2e_dx2_dy1 = -beta * e * ( &
         (yij/rij + yik/rik) * (beta*xij/rij - gama*xjk/rjk) &
         + xij*yij/rij**3 )

   d2e_dy2_dy1 = -beta * e *( &
         (yij/rij + yik/rik) * (beta*yij/rij - gama*yjk/rjk) &
         -1.d0/rij + yij**2/rij**3 )

   d2e_dz2_dy1 = -beta * e *( &
         (yij/rij + yik/rik) * (beta*zij/rij - gama*zjk/rjk) &
         + yij*zij/rij**3 )


   d2e_dx2_dz1 = -beta * e *( &
         (zij/rij + zik/rik) * (beta*xij/rij - gama*xjk/rjk) &
         + xij*zij/rij**3 )

   d2e_dy2_dz1 = -beta * e *( &
         (zij/rij + zik/rik) * (beta*yij/rij - gama*yjk/rjk) &
         + yij*zij/rij**3 )

   d2e_dz2_dz1 = -beta * e *( &
         (zij/rij + zik/rik) * (beta*zij/rij - gama*zjk/rjk) &
         -1.d0/rij + zij**2/rij**3 )

   !----------- 1-3

   d2e_dx3_dx1 = -beta * e * ( &
         (xij/rij + xik/rik) * (beta*xik/rik + gama*xjk/rjk) &
         - 1.d0/rik + xik**2/rik**3)

   d2e_dy3_dx1 = -beta * e * ( &
         (xij/rij + xik/rik) * (beta*yik/rik + gama*yjk/rjk) &
         + xik*yik/rik**3 )

   d2e_dz3_dx1 = -beta * e * ( &
         (xij/rij + xik/rik) * (beta*zik/rik + gama*zjk/rjk) &
         + xik*zik/rik**3 )

   d2e_dx3_dy1 = -beta * e * ( &
         (yik/rik + yij/rij)*(beta*xik/rik + gama*xjk/rjk) &
         + yik*xik/rik**3 )

   d2e_dy3_dy1 = -beta * e * ( &
         (yik/rik + yij/rij)*(beta*yik/rik + gama*yjk/rjk) &
         - 1.d0/rik + yik**2/rik**3 )

   d2e_dz3_dy1 = -beta * e * ( &
         (yik/rik + yij/rij)*(beta*zik/rik + gama*zjk/rjk) &
         + yik*zik/rik**3 )

   d2e_dx3_dz1 = -beta * e * ( &
         (zik/rik + zij/rij)*(beta*xik/rik + gama*xjk/rjk) &
         + zik*xik/rik**3 )

   d2e_dy3_dz1 = -beta * e * ( &
         (zik/rik + zij/rij)*(beta*yik/rik + gama*yjk/rjk) &
         + yik*zik/rik**3 )

   d2e_dz3_dz1 = -beta * e * ( &
         (zik/rik + zij/rij)*(beta*zik/rik + gama*zjk/rjk) &
         -1.d0/rik + zik**2/rik**3 )

   !----------- 2-3

   d2e_dx3_dx2 =  e * ( &
         (beta*xij/rij - gama*xjk/rjk)*(beta*xik/rik + gama*xjk/rjk) &
         + gama * (1.d0/rjk - xjk**2/rjk**3) )

   d2e_dy3_dx2 =  e * ( &
         (beta*xij/rij - gama*xjk/rjk)*(beta*yik/rik + gama*yjk/rjk) &
         - gama *             xjk*yjk/rjk**3 )

   d2e_dz3_dx2 =  e * ( &
         (beta*xij/rij - gama*xjk/rjk)*(beta*zik/rik + gama*zjk/rjk) &
         - gama *             xjk*zjk/rjk**3 )


   d2e_dx3_dy2 =  e * ( &
         (beta*yij/rij - gama*yjk/rjk)*(beta*xik/rik + gama*xjk/rjk) &
         -gama *             xjk*yjk/rjk**3 )

   d2e_dy3_dy2 =  e * ( &
         (beta*yij/rij - gama*yjk/rjk)*(beta*yik/rik + gama*yjk/rjk) &
         +gama * (1.d0/rjk - yjk**2/rjk**3))

   d2e_dz3_dy2 =  e * ( &
         (beta*yij/rij - gama*yjk/rjk)*(beta*zik/rik + gama*zjk/rjk) &
         -gama *             yjk*zjk/rjk**3 )


   d2e_dx3_dz2 = e * ( &
         (beta*zij/rij - gama*zjk/rjk)*(beta*xik/rik + gama*xjk/rjk) &
         -gama *             xjk*zjk/rjk**3)

   d2e_dy3_dz2 = e * ( &
         (beta*zij/rij - gama*zjk/rjk)*(beta*yik/rik + gama*yjk/rjk) &
         -gama *             yjk*zjk/rjk**3 )

   d2e_dz3_dz2 = e * ( &
         (beta*zij/rij - gama*zjk/rjk)*(beta*zik/rik + gama*zjk/rjk) &
         +gama * (1.d0/rjk - zjk**2/rjk**3))

   !----------- atom 1 off diagonal

   d2e_dy1_dx1 = beta * e * &
         (beta* (xij/rij + xik/rik)*(yij/rij + yik/rik) &
         + (xik*yik/rik**3 + xij*yij/rij**3))

   d2e_dz1_dx1 = beta * e * &
         (beta* (xij/rij + xik/rik)*(zij/rij + zik/rik) &
         + (xik*zik/rik**3 + xij*zij/rij**3))

   d2e_dz1_dy1 = beta * e * &
         (beta* (yij/rij + yik/rik)*(zij/rij + zik/rik) &
         + (yik*zik/rik**3 + yij*zij/rij**3))

   !----------- atom 2 off diagonal

   d2e_dy2_dx2 = e * ( &
         (beta*xij/rij - gama*xjk/rjk)*(beta*yij/rij - gama*yjk/rjk) &
         + (beta*xij*yij/rij**3 + gama*xjk*yjk/rjk**3))

   d2e_dz2_dx2 = e * ( &
         (beta*xij/rij - gama*xjk/rjk)*(beta*zij/rij - gama*zjk/rjk) &
         + (beta*xij*zij/rij**3 + gama*xjk*zjk/rjk**3))

   d2e_dz2_dy2 = e *( &
         (beta*yij/rij - gama*yjk/rjk)*(beta*zij/rij - gama*zjk/rjk) &
         + (beta*yij*zij/rij**3 + gama*yjk*zjk/rjk**3))

   !----------- atom 3 off diagonal

   d2e_dy3_dx3 = e * ( &
         (beta*xik/rik + gama*xjk/rjk)*(beta*yik/rik + gama*yjk/rjk) &
         +  (beta*yik*xik/rik**3 + gama*yjk*xjk/rjk**3))

   d2e_dz3_dx3 = e * ( &
         (beta*xik/rik + gama*xjk/rjk)*(beta*zik/rik + gama*zjk/rjk) &
         +  (beta*xik*zik/rik**3 + gama*xjk*zjk/rjk**3))

   d2e_dz3_dy3 = e * ( &
         (beta*yik/rik + gama*yjk/rjk)*(beta*zik/rik + gama*zjk/rjk) &
         +  (beta*yik*zik/rik**3 + gama*yjk*zjk/rjk**3))

   !----------- 45


   call loadit(index,3*(i-1),3*(j-1))
   call loadit(jndex,3*(j-1),3*(k-1))
   call loadit(kndex,3*(i-1),3*(k-1))


   ! diagonal terms

   h(index(1,1)) = h(index(1,1)) + d2e_dx1_2
   h(index(2,2)) = h(index(2,2)) + d2e_dy1_2
   h(index(3,3)) = h(index(3,3)) + d2e_dz1_2

   h(index(4,4)) = h(index(4,4)) + d2e_dx2_2
   h(index(5,5)) = h(index(5,5)) + d2e_dy2_2
   h(index(6,6)) = h(index(6,6)) + d2e_dz2_2

   h(kndex(4,4)) = h(kndex(4,4)) + d2e_dx3_2
   h(kndex(5,5)) = h(kndex(5,5)) + d2e_dy3_2
   h(kndex(6,6)) = h(kndex(6,6)) + d2e_dz3_2

   ! one center term

   h(index(1,2)) = h(index(1,2)) + d2e_dy1_dx1
   h(index(1,3)) = h(index(1,3)) + d2e_dz1_dx1
   h(index(2,3)) = h(index(2,3)) + d2e_dz1_dy1

   h(index(4,5)) = h(index(4,5)) + d2e_dy2_dx2
   h(index(4,6)) = h(index(4,6)) + d2e_dz2_dx2
   h(index(5,6)) = h(index(5,6)) + d2e_dz2_dy2

   h(kndex(4,5)) = h(kndex(4,5)) + d2e_dy3_dx3
   h(kndex(4,6)) = h(kndex(4,6)) + d2e_dz3_dx3
   h(kndex(5,6)) = h(kndex(5,6)) + d2e_dz3_dy3

   ! two center term

   h(index(1,4)) = h(index(1,4)) + d2e_dx2_dx1
   h(index(1,5)) = h(index(1,5)) + d2e_dy2_dx1
   h(index(1,6)) = h(index(1,6)) + d2e_dz2_dx1

   h(index(2,4)) = h(index(2,4)) + d2e_dx2_dy1
   h(index(2,5)) = h(index(2,5)) + d2e_dy2_dy1
   h(index(2,6)) = h(index(2,6)) + d2e_dz2_dy1

   h(index(3,4)) = h(index(3,4)) + d2e_dx2_dz1
   h(index(3,5)) = h(index(3,5)) + d2e_dy2_dz1
   h(index(3,6)) = h(index(3,6)) + d2e_dz2_dz1

   !2
   h(kndex(1,4)) = h(kndex(1,4)) + d2e_dx3_dx1
   h(kndex(1,5)) = h(kndex(1,5)) + d2e_dy3_dx1
   h(kndex(1,6)) = h(kndex(1,6)) + d2e_dz3_dx1

   h(kndex(2,4)) = h(kndex(2,4)) + d2e_dx3_dy1
   h(kndex(2,5)) = h(kndex(2,5)) + d2e_dy3_dy1
   h(kndex(2,6)) = h(kndex(2,6)) + d2e_dz3_dy1

   h(kndex(3,4)) = h(kndex(3,4)) + d2e_dx3_dz1
   h(kndex(3,5)) = h(kndex(3,5)) + d2e_dy3_dz1
   h(kndex(3,6)) = h(kndex(3,6)) + d2e_dz3_dz1

   !3
   if(k < j) then
      call loadit(jndex,3*(k-1),3*(j-1))
   end if
   h(jndex(1,4)) = h(jndex(1,4)) + d2e_dx3_dx2
   h(jndex(1,5)) = h(jndex(1,5)) + d2e_dy3_dx2
   h(jndex(1,6)) = h(jndex(1,6)) + d2e_dz3_dx2

   h(jndex(2,4)) = h(jndex(2,4)) + d2e_dx3_dy2
   h(jndex(2,5)) = h(jndex(2,5)) + d2e_dy3_dy2
   h(jndex(2,6)) = h(jndex(2,6)) + d2e_dz3_dy2

   h(jndex(3,4)) = h(jndex(3,4)) + d2e_dx3_dz2
   h(jndex(3,5)) = h(jndex(3,5)) + d2e_dy3_dz2
   h(jndex(3,6)) = h(jndex(3,6)) + d2e_dz3_dz2

   !4
   return
end subroutine thr_2nd 
