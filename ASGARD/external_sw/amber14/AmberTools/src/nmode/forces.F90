
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
!+ [Enter a one-line description of subroutine forces here]
subroutine forces (x, ix, ih, ipair, ndrv)
   
   implicit double precision (a-h, o-z)
   dimension  x(*)
   dimension ix(*)
   character(len=4) ih(*)
   
   !     ----- Input parameters:
   !     -----   ipair   = 1 read pairlist from file PRLIST
   !     -----           = 2 make fresh pairlist now
   !     -----           = 3 make pairlist now and write it
   !                           to file PRLIST
   !     -----           = 4 use pairlist already in memory
   !     -----           = -2 or -3 make pairlist and convert to indices
   !     -----           =          corresponding to rearranged atomic list
   
   !     -----   ndrv    = 0 energy only
   !     -----           = 1 energy and force
   !     -----           = 2 energy, force and second derivatives
   

#  include "pointer.h"
#  include "inpdat.h"
#  include "epot.h"
#  include "pol.h"
#  include "sizes2.h"

   dimension scr1(3*maxatom),enew(3*maxatom),scr2(3*maxatom), &
         scr3(3*maxatom),scr4(3*maxatom)
   dimension scr5(maxatom),scr7(maxatom),scr8(maxatom)
   dimension iscr1(maxatom),iscr2(maxatom),iscr3(maxatom)
   dimension vt(4)
   dimension itrp(3,maxatom)
   imgslt = 0
   iptatm = 0

   
   lat = 0
   ntb = 0
   if(natom > maxatom) then
      write(6,'(a,i7,i7)') 'too many atoms: ',natom,maxatom
      write(6,'(a)') 'change MAXATOM in sizes2.h'
      call exit
   end if

   call putvar(natom,ix(mgroup),x(mxbel),x(mx))
   
   iapair = iabs(ipair)
   if(iapair == 1) &
         call pairr (natom,     npair,     kpair,     nhb, &
         ix(miar1), ix(miar2))
   if(iapair == 2 .or. iapair == 3) &
         call bresna (natom,      nres,      npair,     nhb, &
         ix(mipres),  ix(miac),  ix(mico),  ix(miblo), &
         ix(minb),    ix(miar1), ix(miar2), x(mx), &
         ix(migres),  ntb,       cut,       ntypes, &
         nbmax)
   if(iapair == 3) &
         call pairw (natom,     npair,     kpair,     nhb, &
         ix(miar1), ix(miar2))
   if(iapair <= 3) then
      write (6,1) npair, nhb
      1 format (/4x, 'Number of non-bonded pairs = ', i10, &
            /4x, 'Number of H-bonded pairs   = ', i10)
      if(npair > nbmax) call errstp(2,nbmax)
   end if
   
   if(n3b /= 0.and.itriple == 0) then
      if(n3b > 0) then
         call tripl(natom,ix(miar1),ix(miar2),ih(msymbl),itrp)
      else
         call ugly(n3b,itrp)
      end if
      itriple = 1
   end if
   
   call zerout (x(mf), nr3)
   if(ndrv > 1) then
      ntemp = ns3*(ns3+1)/2
      call zerout (x(mh), ntemp)
   end if
   
   if(ndrv == 2) then
      call nonbon (natom,     npair,    ix(miar1), ix(miar2), &
            ix(miac),  ix(mico), x(mx),     x(mf), &
            x(mh),     x(mcn1), &
            x(mcn2),   x(masol), x(mbsol),  x(mhbcut), &
            x(mchrg),  x(mchrg), enba,      enbh, &
            eel,       idiel,    ntypes,    ndrv, &
            ix(mnbel),  natsys)
   else
      call nonbon1 (natom,    npair,     ix(miar1), ix(miar2), &
            ix(miac), ix(mico),  x(mx),     x(mf), &
            dummy,    x(mcn1),   x(mcn2),   x(masol), &
            x(mbsol), x(mhbcut), x(mchrg),  x(mchrg), &
            enba,     enbh,      eel,       idiel, &
            ntypes,   ndrv)
   end if
   if(ipol /= 0) then
      !        subroutine politr(natom,nres,c,q,pol,ipair,jpair,
      !    -                  eold,enew,eper,epol,xwij,r2,ipres,
      !    -                  aveper,aveind,avetot,emtot,nstep,
      !    -                  iptsol,n14,ni14,iarx,jpw,imgslt
      !    -                  ,iptatm,scr1,scr2)
      
      call politr(natom,nres,x(mx),x(mchrg),x(mpol),ix(miar1), &
            ix(miar2), &
            scr1,enew,scr3,epol,scr4,scr5,ix(mipres), &
            aveper,aveind,avetot,emtot,nstep, &
            iptsol,n14,ni14,iscr1,iscr2,imgslt, &
            iptatm,scr7,scr8)
      
      !     subroutine pol2der(natom,x,f,p,q,h,ipair,
      !    -                  jpair,xrc,xij,r2,fw,vt,
      !    -                  n14,ni14,iarx,jpw,imgslt,iptatm,ndrv)
      

      call pol2der (natom,x(mx),x(mf),enew,x(mchrg),x(mh), &
            ix(miar1),ix(miar2),scr1,scr3,scr4,scr5,vt, &
            n14,ni14,iscr1,iscr3,imgslt,iptatm,ndrv)
   end if

   if(n3b /= 0) then
      e3bod = 0.0d0
      call threeb(x(mx),x(mf),itrp,e3bod,natom)
      do jn = 1,n3b
         do kn = 1,i3b(jn)
            ii = itrp(1,kn)
            jj = itrp(2,kn)
            kk = itrp(3,kn)

            call thr_2nd(x(mx),x(mh),ii,jj,kk, &
                  e3bod,beta3b(jn),gama3b(jn),ih(mgraph))
         end do
      end do
   end if


   call bond (nbonh,     ix(mibh),  ix(mjbh), ix(micbh), &
         x(mrk),    x(mreq),   ntb,      x(mx), &
         x(mf),     x(mh),     ebh,      ndrv,      ix(mnbel), &
         natsys)

   call bond (nbona,     ix(miba),  ix(mjba), ix(micba), &
         x(mrk),    x(mreq),   ntb,      x(mx), &
         x(mf),     x(mh), &
         eba,       ndrv,      ix(mnbel), &
         natsys)

   call angl (ntheth,    ix(mith),  ix(mjth), ix(mkth), &
         ix(micth), x(mtk),    x(mteq),  x(mfk), &
         x(mfpk),   x(mqeq),   ntb,      x(mx), &
         x(mf),     x(mh), &
         eah,       ndrv,      ix(mnbel), &
         natsys)

   call angl (ntheta,    ix(mita),  ix(mjta), ix(mkta), &
         ix(micta), x(mtk),    x(mteq),  x(mfk), &
         x(mfpk),   x(mqeq),   ntb,      x(mx), &
         x(mf),     x(mh), &
         eaa,       ndrv,      ix(mnbel), &
         natsys)

   call ephi (nphih,     ix(miph),  ix(mjph), ix(mkph), &
         ix(mlph),  ix(micph), x(mpk),   x(mpn), &
#ifdef CHARMM
         x(mphase), x(mcn114), x(mcn214),ntypes, &
#else
         x(mphase), x(mcn1),   x(mcn2),  ntypes, &
#endif
         x(mchrg),  ix(mico),   ix(miac),  ntb, &
         x(mx),     x(mf),     x(mh), &
         edihh,     enb14h,    eel14h, &
         scnb,      scee,      nphih,    epcnh, &
         idiel,     ndrv,      ix(mnbel), natsys)
   

   call ephi (nphia,     ix(mipa),  ix(mjpa), ix(mkpa), &
         ix(mlpa),  ix(micpa), x(mpk),   x(mpn), &
#ifdef CHARMM
         x(mphase), x(mcn114), x(mcn214),ntypes, &
#else
         x(mphase), x(mcn1),   x(mcn2),  ntypes, &
#endif
         x(mchrg),  ix(mico),   ix(miac),  ntb, &
         x(mx),     x(mf),     x(mh), &
         ediha,     enb14a,    eel14a, &
         scnb,      scee,      nphia,    epcna, &
         idiel,     ndrv,      ix(mnbel), natsys)

   if(icons == 1) call xconst(econs,x(mx),x(mf),x(mh),ndrv, &
         x(mxref),x(mwref),natom)
   !!
   !!     ----- put the belly force vector into the first part of array f

   call getf(natom,ix(mgroup),x(mf))

   eb    = ebh    + eba
   ea    = eah    + eaa
   edih  = edihh  + ediha
   enb   = enba

   enb14 = enb14h + enb14a
   eel14 = eel14h + eel14a

   etot  = eb + ea + edih + enb + eel + enb14 + eel14 + enbh + econs &
         + e3bod + epol

   return
end subroutine forces 
