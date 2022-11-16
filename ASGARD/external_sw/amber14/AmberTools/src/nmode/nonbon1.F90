
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
!+ [Enter a one-line description of subroutine nonbon1 here]
subroutine nonbon1(natom,npair,iar1,iar2,iac,ico,x,f,dd,cn1, &
      cn2,asol,bsol,hbcut,cg,xchrg,enb,ehb,eel, &
      idiel,ntypes,ndrv)
   
   !     this routine calculates the electrostatic and either the
   !     van der waals -or- hbond interactions of each atom i of the
   !     system with all non-bonded atoms j within a specified sphere
   !     of influence.  this version of nonbon runs fastest on a
   !     cray xmp-4 with hardware gather/scatter capabilities. it will
   !     run faster than the original nonbon on any cray, however.
   
   !     variable      function
   !     ------------------------------------------------------------------
   !     iar1()        number of nonbonded atoms j for each atom i
   !     iar2()        nonbonded pairlist. iar2(1) to iar2(iar1(1)) are
   !                     all the atoms j interacting with atom no. 1
   !     iac()         amber atom type of atom i
   !     x(,)          coordinates of the atoms x(1,i)=x, x(2,i)=y, x(3,i)=z
   !     f(,)          force on atoms, x,y,z as x(,)
   !     cn1()         vdw a parameter
   !     cn2()         vdw b parameter
   !     asol()        hbond a parameter
   !     bsol()        hbond b parameter
   !     cg()          charge on atom i
   !     p()           pointer to atoms in current nb pairlist (j)
   !     hbptr()       pointer to hbonding atoms in current nb pairlist (j)
   !     vdwptr()      pointer to atoms with vdw interactions in current pairlist
   !     iparmp(,)     pointer to parameters.  indices are atom types.  set in
   !                     subroutine  pload
   !     ------------------------------------------------------------------
   implicit double precision (a-h,o-z)
   logical dield
   
   common/crayon/ iparmp(50,50)
   
   dimension iar1(*),iar2(*),iac(*),ico(*),cg(*)
   dimension xchrg(*),cn1(*),cn2(*),asol(*),bsol(*),hbcut(*)
   
   dimension x(3,*),f(3,*)
   
   !     -- scratch arrays --
   
   !        MAXPR is the maximum number of non-bonded partners for any atom:
   parameter (maxpr=2800)
   
   !        MAXAT is the maximum number of atoms in the system:
   parameter (maxat=7500)
   
   dimension ic(maxpr),xij(maxpr),yij(maxpr),zij(maxpr), &
         s(maxpr),r2(maxpr),r6(maxpr),r10(maxpr),r12(maxpr), &
         g(maxpr),df(maxpr),df1(maxpr),df2(maxpr),hf1(maxpr), &
         hf2(maxpr),vf1(maxpr),vf2(maxpr)
   dimension icx(0:maxpr)
   dimension tempxi(0:maxpr), tempyi(0:maxpr), tempzi(0:maxpr)
   double precision dx(0:maxat), dy(0:maxat), dz(0:maxat)
   double precision dasol(0:200),dbsol(0:200), &
         dcn1(0:1830),dcn2(0:1830)
   
   !     -- pointer arrays --
   
   integer p(maxpr),hbptr(maxpr),vdwptr(maxpr)
   
   if (natom > maxat) then
      write(6,1) natom
      1 format('Stop in nonbon1, natom=',i6)
      call mexit(6, 1)
   end if
   do 10 j = 1, natom
      dx(j) = x(1,j)
      dy(j) = x(2,j)
      dz(j) = x(3,j)
   10 continue

   dasol(0) = 0.0
   dbsol(0) = 0.0
   do 20 j = 1, 200
      dasol(j) = asol(j)
      dbsol(j) = bsol(j)
   20 continue

   dcn1(0) = 0.0
   dcn2(0) = 0.0
   do 30 j = 1, 1830
      dcn1(j) = cn1(j)
      dcn2(j) = cn2(j)
   30 continue

   enb = 0.0d+00
   eel = 0.0d+00
   ehb = 0.0d+00
   lpair = 1
   lim = natom-1
   
   !     loop over all atoms i in the system
   
   do 1000 i=1,lim
      
      npr = iar1(i)
      if (npr == 0) goto 1000
      if (npr > maxpr) then
         write(6,2) npr
         2 format('Stop in nonbon1, npr=',i6)
         call mexit(6, 1)
      end if
      
      !     load atom pointer array p()
      
      kount = 0
      do 100 j = lpair, lpair + npr - 1
         kount = kount + 1
         p(kount) = iar2(j)
      100 continue
      lpair = lpair + npr
      
      !     load parameter pointer array ic()
      
      lowic = 0
      do 150 j = 1, npr
         ic(j) = iparmp(iac(i),iac(p(j)))
         lowic = min (lowic, ic(j))
      150 continue
      
      !     calculate interatomic distances and reciprocals
      
      do 200 j = 1, npr
         xij(j) = dx(i) - dx(p(j))
         yij(j) = dy(i) - dy(p(j))
         zij(j) = dz(i) - dz(p(j))
         s(j) = xij(j)**2 + yij(j)**2 + zij(j)**2
         r2(j) = 1.0d0 / s(j)
         r6(j) = r2(j)**3
         r12(j) = r6(j) * r6(j)
      200 continue
      
      if (idiel == 1) then
         
         !     -- constant dielectric --
         
         dsum = 0.0d0
         do 250 j = 1, npr
            g(j) = cg(i) * cg(p(j)) * sqrt(r2(j))
            df2(j) = -g(j)
            dsum = dsum + g(j)
         250 continue
         
         eel = eel + dsum
      else if (idiel == 0) then
         
         !     -- distance dependent dielectric --
         
         320 continue
         dsum = 0.0d0
         do 300 j = 1, npr
            g(j) = cg(i) * cg(p(j)) * r2(j)
            df2(j) = -(g(j) + g(j))
            dsum = dsum + g(j)
         300 continue
         
         eel = eel + dsum
      else if (idiel == 2) then
         
         !         ---- sigmoidal dielectric ------
         !                 (d = 78, s=0.3 hardwired)
         
         dsum = 0.0
         do 324 j=1,npr
            rwq = 1./sqrt(r2(j))
            sr = 0.3*rwq
            expsr = exp(-sr)
            srsq = sr*sr
            epsr = rwq*(79.0 - 39.*expsr*(srsq + sr + sr + 2))
            g(j) = cg(i)*cg(p(j))/epsr
            df2(j) = -g(j)*(1.0 + 1.053*expsr/(r2(j)*r2(j)*epsr))
            dsum = dsum + g(j)
         324 continue
         
         eel = eel + dsum
      end if
      
      !     -- check for no hbonds in current pairlist
      
      330 if (lowic < 0) goto 340
      
      !     -- do vdw energy this way when there are no hbonds
      
      do 332 j = 1, npr
         vf1(j) = dcn1(ic(j)) * r12(j)
         vf2(j) = dcn2(ic(j)) * r6(j)
         df1(j) = -12.0d0 * vf1(j) + 6.0d0 * vf2(j)
         df(j) = (df1(j) + df2(j)) * r2(j)
         if (vf1(j) > 50.) write(6,331) i,p(j),s(j),vf1(j)
         331 format('Bad nonbon:',2i5,2e12.5)
      332 continue
      
      dsum = 0.0d0
      do 334 j = 1, npr
         dsum = dsum + vf1(j)
      334 continue
      enb = enb + dsum
      dsum = 0.0d0
      do 335 j = 1, npr
         dsum = dsum + vf2(j)
      335 continue
      enb = enb - dsum
      goto 575
      
      
      340 continue
      
      !     make hbond parameter pointers positive
      !       and throw out h-bonds longer than hbcut
      
      do 342 j = 1,npr
         if (ic(j) < 0) then
            icx(j) = -ic(j)
            !                   ***note: amber 3.0 no longer uses hbcut!
            !         if (s(j).gt.hbcut(icx(j))) icx(j) = 0
         else
            icx(j) = 0
         end if
      342 continue

      do 343 j = 1, npr
         if (ic(j) < 0) ic(j) = 0
      343 continue

      
      !     -- hbond pairs 10 - 12 potential --
      
      do 344 j = 1, npr
         r10(j) = r6(j) * r2(j) * r2(j)
         hf1(j) = dasol(icx(j)) * r12(j)
         hf2(j) = dbsol(icx(j)) * r10(j)
         vf1(j) = dcn1(ic(j)) * r12(j)
         vf2(j) = dcn2(ic(j)) * r6(j)
         df1(j) = -12.0d0 * (hf1(j) + vf1(j)) + 10.0d0 * hf2(j) + &
               6.0d0 * vf2(j)
         df(j) = (df1(j) + df2(j)) * r2(j)
         if(vf1(j) > 50.) write(6,331) i,p(j),s(j),vf1(j)
      344 continue
      
      !     ehb = ehb + dsum(nhbond,hf1,1) - dsum(nhbond,hf2,1)
      dsum = 0.0d0
      do 460 j = 1, npr
         dsum = dsum + hf1(j)
      460 continue
      ehb = ehb + dsum
      dsum = 0.0d0
      do 470 j = 1, npr
         dsum = dsum + hf2(j)
      470 continue
      ehb = ehb - dsum
      
      !     enb = enb + dsum(nvdw,vf1,1) - dsum(nvdw,vf2,1)
      dsum = 0.0d0
      do 510 j = 1, npr
         dsum = dsum + vf1(j)
      510 continue
      enb = enb + dsum
      dsum = 0.0d0
      do 520 j = 1, npr
         dsum = dsum + vf2(j)
      520 continue
      enb = enb - dsum
      
      
      575 if (ndrv <= 0) goto 1000
      
      !     -- update the force array --
      
      !$dir no_recurrence
      do 600 j = 1, npr
         tempxi(j) = xij(j) * df(j)
         tempyi(j) = yij(j) * df(j)
         tempzi(j) = zij(j) * df(j)
         f(1,i) = f(1,i) - tempxi(j)
         f(2,i) = f(2,i) - tempyi(j)
         f(3,i) = f(3,i) - tempzi(j)
      600 continue
      !$dir no_recurrence
      do 700 j = 1, npr
         f(1,p(j)) = f(1,p(j)) + tempxi(j)
         f(2,p(j)) = f(2,p(j)) + tempyi(j)
         f(3,p(j)) = f(3,p(j)) + tempzi(j)
      700 continue
   1000 continue
   return
end subroutine nonbon1 
