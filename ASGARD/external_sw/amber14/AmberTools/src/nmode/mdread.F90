
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
!+ [Enter a one-line description of subroutine mdread here]
subroutine mdread(x)
   
   !------reads in main input data for the program
   
   implicit double precision (a-h,o-z)
   character(len=8) intname
   integer blank
#  include "sizes2.h"
#  include "anal.h"
#  include "infoa.h"
#  include "optf.h"
#  include "opta.h"
#  include "files.h"
   common/consnb/npair,ntypes,idiel,iyyy,dielc,scnb,scee
   common/nbterm/cut,nbuck,iblo(maxatom),inb(9000)
   common/minpar/ntrun,maxcyc,ncyc,iopt,nvect,izzz,dxm,dele,drms
   common/runhed/ihead(20),ihead1(20)
   common/belly/ibelly,natbel,igroup(maxatom),isymbl(maxatom), &
         itree(maxatom),igres(600)
   common/ntrun5/nrvec,nrat,iat,jat,imov
   common/tcor/ tmax, tintvl, intname(maxint)
   common/ired/ vtemp, iorder, npairs, langevin
   dimension xc(2),weit(2),x(*)
   
   namelist/data/ntrun,nvect,pcut,ibelly,ntx,ibeg,iend,nrgrp, &
         irms,nrvec,nrat,iat,jat,imov,ifluc,ipro,kform,ieff,bose, &
         ivform,natom,ihsful,first,last,iskip,tmax,tintvl, &
         iorder,npairs,langevin,vtemp
   
   data blank/4h    /

   call amopen(6,nmdout,owrite,'F','W')
   call amopen(4,nmdin,'O','F','R')
   write(6,40)
   
   !   ---- set the default values for namelist variables
   
   ntrun = 1
   nvect = 50
   pcut = 0.02
   ibelly = 0
   ntx = 1
   ibeg = 7
   iend = 50
   vtemp = 300.
   nrgrp = 0
   irms = 0
   nrvec = 0
   nrat = 0
   iat = 1
   jat = 0
   imov = 0
   ifluc = 0
   ipro = 1
   kform = 2
   ieff = 0
   bose = .false.
   ivform = 1
   natom = 1
   ihsful = 1
   first = 1
   last = 9999
   iskip = 1
   tmax = 0.0
   tintvl = 1.0
   iorder = 2
   npairs = 1
   langevin = 0
   
   !     ----- READ DATA CHARACTERIZING THE RUN -----
   
   read(4,70) ihead
   
   ! ---- read in namelist variables
   
   read(4,data)
   if (natom > maxatom) then
      write(6,*) 'natom = ',natom,' is too big; max is ',maxatom
      call mexit(6,1)
   end if
   
   !     ----- READ MOLECULAR TOPOLOGY FILE, PRINT IHEAD -----
   
   if (ntrun > -1 .and. ntrun < 7) then
      call rdparm2()
      
      write(6,'(/ /,a,/ /)') &
            '   1.  M O L E C U L A R   T O P O L O G Y   F I L E :'
      write(6,60) ihead1
   else
      do i=1,natom
         igrph(i) = blank
      end do
   end if
   
   !  --- set default values for some variables
   
   ntf = 1
   ntn = 1
   ntb = 0
   npm = 1
   nrp = natom
   nsm = 0
   nram = 0
   nrpt = npm*nrp
   if(jat == 0) jat=natom
   
   !---load the belly atoms if required
   
   if(ibelly == 0) then
      do 20 i=1,natom
      20 igroup(i) = 1
      natbel = natom
      write(6,30) natbel
      30 format('No belly:',i4,' atoms in moving group')
   else
      ngrp = 0
      write(6,80)
      call rgroup(natom,natbel,nres,ngrp,ipres,lbres,igrph, &
            isymbl,itree,igroup,weit,xc,.false.,.true.,4)
   end if
   if(ibeg < 1) ibeg = 1
   if(iend > (3*natbel)) iend = 3*natbel
   
   ! --- write out data characterizing the run
   
   write(6,data)
   
   return
   
   40 format(/ /10x,54(1h*),/10x,'INITIATE THE NMANAL MODULE',/10x,54(1h*),/ /)
   60 format(2x,20a4)
   70 format(20a4)
   80 format('0Loading belly atoms into igroup array')
end subroutine mdread 
