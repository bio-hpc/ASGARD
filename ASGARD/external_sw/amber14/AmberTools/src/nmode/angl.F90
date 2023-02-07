
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
!+ [Enter a one-line description of subroutine angl here]
subroutine angl(nba, it, jt, kt, ict, tk, teq, fk, fpk, qeq, ntb, x, &
      f, dd, eba, ndrv, nbel, natsys)
   
   implicit double precision (a-h,o-z)
   
   !     ----- routine to get the angle energies and forces for the
   !           potential of the type ct*(t-t0)**2
   
   dimension it(nba),jt(nba),kt(nba),ict(2),tk(2), &
         teq(2),fk(2),fpk(2)
   dimension x(2),f(2),xxik(3),qeq(2)
   dimension xxij(3),xxkj(3),s(6),dc(6),ddc(6,6),dr(9),ddr(9,9)
   dimension dd(*), chi(2), add(9,9),nbel(*)
   equivalence (xxij(1),s(1)),(xxkj(1),s(4)),(xxik(1),xxij(1))
   ia(i) = i*(i-1)/2
   ibel(m3) = nbel(m3/3+1)
   
   eba = 0.0d0
   eub = 0.0d0
   
   do n = 1,nba
      
      i3 = it(n)
      j3 = jt(n)
      k3 = kt(n)
      ic = ict(n)
      rij2 = 0.d0
      rkj2 = 0.d0
      do m = 1,3
         xj = x(j3+m)
         xxij(m) = x(i3+m)-xj
         xxkj(m) = x(k3+m)-xj
         rij2 = rij2+xxij(m)**2
         rkj2 = rkj2+xxkj(m)**2
      end do
      if(ntb /= 0) then
         call percon(rij2,xxij)
         call percon(rkj2,xxkj)
      end if
      rij = sqrt(rij2)
      rkj = sqrt(rkj2)
      rrik = rij*rkj
      cst = (s(1)*s(4)+s(2)*s(5)+s(3)*s(6))/rrik
      if(cst > 1.0d0) cst = 1.0d0
      if(cst < -1.0d0) cst = -1.0d0
      at = acos(cst)
      
      !     ----- calculation of the energy -----
      
      da = at-teq(ic)
      df = tk(ic)*da
      ebah = df*da
      eba = eba+ebah
      if (ebah > 100.0d0) then
         i3p = (i3+3)/3
         j3p = (j3+3)/3
         k3p = (k3+3)/3
         write(6,101) i3/3+1,j3/3+1,k3/3+1,at,teq(ic),tk(ic),ebah
         101 format(' Bad angle:',3i5,4e12.5)
      end if
      
      !     ----- calculation of the force -----
      
      if(ndrv > 0) then
         st = sin(at)
         if(dabs(st) < 1.0d-08) st = 1.0d-08
         df= df+df
         ddf = tk(ic)+tk(ic)
         
         !     ----- now set up first and second derivatives -----
         
         !     ----- set up dc = first derivative of costheta with respect to
         !           cartesian differences -----
         
         do i=1,3
            dc(i) = (s(i+3)/rkj-cst*s(i)/rij)/rij
            dc(i+3) = (s(i)/rij-cst*s(i+3)/rkj)/rkj
         end do
      end if
      if(ndrv > 1) then
         call difang(cst,s,dc,ddc,nbel)
         
         !     ----- now change ddc to be the second derivative of theta
         !           w/respect to cartesian differences
         
         st2c = cst/(st*st)
         do i=1,6
            do j=i,6
               ddc(i,j) = -(ddc(i,j)+dc(i)*dc(j)*st2c)/st
               ddc(j,i) = ddc(i,j)
            end do
         end do
      end if
      
      !     ----- now make dc to be the first derivative of theta w/respect
      !           to cartesian differences -----
      
      if (ndrv >= 1) then
         do i=1,6
            dc(i) = -dc(i)/st
         end do
         
         !     ----- set up dr = -first derivatives of theta with respect
         !           to cartesians -----
         
         do i=1,3
            dr(i) = -dc(i)
            dr(i+6) = -dc(i+3)
            dr(i+3) = dc(i)+dc(i+3)
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
      end if
      if (ndrv == 2) then
         
         !     ----- set up ddr = second derivatives of theta w/resp.
         !           to cartesians -----
         
         do i=1,3
            do j=i,3
               ddr(i,j) = ddc(i,j)
               ddr(i+6,j+6) = ddc(i+3,j+3)
            end do
            do j=1,3
               ddr(i,j+6) = ddc(i,j+3)
               ddr(i,j+3) = -ddc(i,j)-ddc(i,j+3)
               ddr(i+3,j+6) = -ddc(i+3,j)-ddc(i+3,j+3)
            end do
         end do
         
         do i=1,3
            do j=i,3
               ddr(i+3,j+3) = -ddr(i,j+3)-ddr(i+3,j+6)
            end do
         end do
         
         !     ----- FORM THE DD ARRAY -----
         
#ifdef NEC_SX
         !         -- Novector required for NEC compiler problem
         !vdir novector
#endif
         do i=1,9
            ist = ibel(k3)
            if(i <= 6) ist = ibel(j3)
            if(i <= 3) ist = ibel(i3)
            if(ist /= -1) then
               iof = mod(i,3)
               if(iof == 0) iof = 3
               inew = ist+iof
               
#ifdef NEC_SX
               !             -- Novector required for NEC compiler problem
               !vdir novector
#endif
               do j=i,9
                  jst = ibel(k3)
                  if(j <= 6) jst = ibel(j3)
                  if(j <= 3) jst = ibel(i3)
                  if(jst /= -1) then
                     jof = mod(j,3)
                     if(jof == 0) jof = 3
                     jnew = jst+jof
                     
                     istore = ia(jnew)+inew
                     if(jnew < inew) istore = ia(inew)+jnew
                     dd(istore) = dd(istore) + ddf*dr(i)*dr(j) + df*ddr(i,j)
                  end if
               end do
            end if
         end do
      end if
#ifdef CHARMM
      !     ---- compute Urey-Bradley terms:
      
      i3 = it(n)
      j3 = kt(n)
      rij2 = 0.d0
      do m = 1,3
         xxij(m) = x(i3+m)-x(j3+m)
         rij2 = rij2+xxij(m)**2
      end do
      if(ntb /= 0) call percon(rij2,xxij)
      rij = sqrt(rij2)
      db = rij-qeq(ic)
      df = fk(ic)*db
      ebh = df*db
      eub = eub+ebh
      if(ndrv == 2) then
         dv1 = 2.0d0*df
         dv2 = 2.0d0*fk(ic)
         call difbon(dd,dv1,dv2,xxij,rij,i3,j3,nbel)
      end if
      df = 2.0d0*df/rij
      do m = 1,3
         xh = xxij(m)*df
         f(i3+m) = f(i3+m)-xh
         f(j3+m) = f(j3+m)+xh
      end do
#endif
      
      !     ----- finish loop over angles -----
      
   end do
   
#ifdef CHARMM
   !     write(6,*) 'Urey-Bradley energy: ',eub
   eba = eba + eub
#endif
   return
end subroutine angl 
