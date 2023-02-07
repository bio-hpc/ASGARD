
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
!+ [Enter a one-line description of subroutine putvar here]
subroutine putvar(nat,igrp,x,xp)
   real*8 x,xp
   
   !     ----- ROUTINE TO STUFF THE VARYING COORDINATES INTO THE
   !           COORDINATE ARRAY FOR FORCE -----
   
   dimension igrp(*),x(*),xp(*)
   
   j3 = 0
   do 10 i = 1,nat
      if (igrp(i) > 0) then
         i3 = 3*i-3
         xp(i3+1) = x(j3+1)
         xp(i3+2) = x(j3+2)
         xp(i3+3) = x(j3+3)
         j3 = j3+3
      end if
   10 continue
   return
end subroutine putvar 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine setvar here]
subroutine setvar
   implicit double precision (a-h,o-z)
   
   !     ----- ROUTINE TO DO THE NECESSARY ACCOMODATIONS FOR PROTEIN
   !           BELLY MINIMISATIONS -----
   
#  include "sizes2.h"
#  include "bad.h"
#  include "carts.h"
#  include "infoa.h"
#  include "optf.h"
#  include "opta.h"
   common/belly/ibelly,natbel,igroup(maxatom),isymbl(maxatom), &
         itree(maxatom),igres(600)
   
   common/consnb/npair,ntypes,idiel,iyyy,dielc,scnb,scee
   common/runhed/ihead(20),ihead1(20)
   common/ibel/nbel(3000),mbel(3000)

   
   !     ----- DELETE BONDS WHICH ARE IN THE BELLY ALONE -----
   
   if (ibelly /= 0) then
      mbonh = nbonh
      call setbon(nbonh,mbonh,ibh,jbh,icbh,igroup)
      call setbon(nbona,mbona,iba,jba,icba,igroup)
      
      !     ----- DELETE THE BONDS WHICH ARE IN THE BELLY ALONE -----
      
      mtheth = ntheth
      call setang(ntheth,mtheth,ith,jth,kth,icth,igroup)
      call setang(ntheta,mtheta,ita,jta,kta,icta,igroup)
      
      !     ----- DELETE THE DIHEDRALS -----
      
      mphih = nphih
      call setdih(nphih,mphih,iph,jph,kph,lph,icph,igroup)
      call setdih(nphia,mphia,ipa,jpa,kpa,lpa,icpa,igroup)
   end if
   
   !     ----- FIND THE TOTAL NUMBER OF ACTIVE ATOMS AND RESIDUES -----
   
   call setatm(natom,nres,natbel,ipres,igroup,igres)
   
   !     ----- MOVE THE ACTIVE COORDINATES and masses TO THE WORKING ARRAY
   !           FOR MINIMISER -----
   
   call setcor(natom,igroup,x,xbel)
   call getm(natom,igroup,amass)
   
   !     ----- set up NBEL array: NBEL(I) = new pointer array for atom i
   
   if (lat /= 0) return
   j = 0
   do 10 i=1,natom
      if (igroup(i) == 0) then
         nbel(i) = -1
      else
         j = j + 1
         nbel(i) = 3*j - 3
      end if
   10 continue
   return
end subroutine setvar 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine setbon here]
subroutine setbon(nb,mb,ib,jb,icb,igrp)
   dimension ib(*),jb(*),icb(*),igrp(*)
   
   nba = 0
   nca = 0
   do 10 i = 1,nb
      iat = ib(i)/3+1
      jat = jb(i)/3+1
      if (igrp(iat) > 0 .or. igrp(jat) > 0) then
         if (i > mb) nca = nca+1
         nba = nba+1
         ib(nba) = ib(i)
         jb(nba) = jb(i)
         icb(nba) = icb(i)
      end if
   10 continue
   nb = nba
   mb = nba-nca
   return
end subroutine setbon 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine setang here]
subroutine setang(nt,mt,it,jt,kt,ict,igrp)
   dimension it(*),jt(*),kt(*),ict(*),igrp(*)
   
   nta = 0
   mta = 0
   do 10 i = 1,nt
      iat = it(i)/3+1
      jat = jt(i)/3+1
      kat = kt(i)/3+1
      if (igrp(iat) > 0 .or. igrp(jat) > 0 .or. igrp(kat) > 0) &
            then
         if (i > mt) mta = mta+1
         nta = nta+1
         it(nta) = it(i)
         jt(nta) = jt(i)
         kt(nta) = kt(i)
         ict(nta) = ict(i)
      end if
   10 continue
   nt = nta
   mt = nta-mta
   return
end subroutine setang 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine setdih here]
subroutine setdih(np,mp,ip,jp,kp,lp,icp,igrp)
   dimension ip(*),jp(*),kp(*),lp(*),icp(*),igrp(*)
   
   npa = 0
   mpa = 0
   do 10 i = 1,np
      iat = ip(i)/3+1
      jat = jp(i)/3+1
      kat = iabs(kp(i))/3+1
      lat = iabs(lp(i))/3+1
      if (igrp(iat) > 0 .or. igrp(jat) > 0 .or. igrp(kat) > 0 &
            .or. igrp(lat) > 0) then
         if (i > mp) mpa = mpa+1
         npa = npa+1
         ip(npa) = ip(i)
         jp(npa) = jp(i)
         kp(npa) = kp(i)
         lp(npa) = lp(i)
         icp(npa) = icp(i)
      end if
   10 continue
   np = npa
   mp = npa-mpa
   return
end subroutine setdih 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine setatm here]
subroutine setatm(nat,nres,natb,ipres,igrp,igres)
   logical active
   dimension ipres(600),igrp(3000),igres(600)
   
   idum = 0
   do 10 i = 1,nat
      if (igrp(i) > 0) idum = idum+1
   10 continue
   natb = idum
   
   do 20 i = 1,nres
      i1 = ipres(i)
      i2 = ipres(i+1)-1
      igres(i) = 0
      active = .false.
      do 15 j = i1,i2
         if (igrp(j) > 0) active = .true.
      15 continue
      if (active) igres(i) = 1
   20 continue
   return
end subroutine setatm 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine setcor here]
subroutine setcor(nat,igrp,xp,x)
   double precision xp,x
   dimension xp(*),x(*),igrp(*)
   j3 = 0
   do 10 i = 1,nat
      if (igrp(i) > 0) then
         i3 = 3*i-3
         x(j3+1) = xp(i3+1)
         x(j3+2) = xp(i3+2)
         x(j3+3) = xp(i3+3)
         j3 = j3+3
      end if
   10 continue
   return
end subroutine setcor 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine getf here]
subroutine getf(nat,igrp,f)
   real*8 f
   
   !     ----- ROUTINE TO EXTRACT THE FORCES OF VARYING
   !           ATOMS FOR THE MINIMISER -----
   
   dimension igrp(*),f(*)
   
   j3 = 0
   do 10 i = 1,nat
      if (igrp(i) > 0) then
         i3 = 3*i-3
         f(j3+1) = f(i3+1)
         f(j3+2) = f(i3+2)
         f(j3+3) = f(i3+3)
         j3 = j3+3
      end if
   10 continue
   return
end subroutine getf 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine getm here]
subroutine getm(nat,igrp,f)
   double precision f
   
   !     ----- ROUTINE TO EXTRACT THE masses OF VARYING
   !           ATOMS FOR THE MINIMISER -----
   
   dimension igrp(*),f(*)
   
   j = 1
   do 10 i = 1,nat
      if (igrp(i) > 0) then
         f(j) = f(i)
         j = j + 1
      end if
   10 continue
   return
end subroutine getm 
