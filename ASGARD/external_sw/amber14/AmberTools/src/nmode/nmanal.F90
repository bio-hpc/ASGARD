
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
!+ [Enter a one-line description of program nmanal here]
program nmanal
   
   implicit double precision (a-h,o-z)
#  include "sizes2.h"
   character(len=8) intname,type
   character(len=4) star
   character(len=40) title
   
#  include "anal.h"
#  include "bad.h"
#  include "files.h"
#  include "infoa.h"
#  include "carts.h"
#  include "optf.h"
#  include "opta.h"
   common/consnb/npair,ntypes,idiel,iyyy,dielc,scnb,scee
   common/minpar/ntrun,maxcyc,ncyc,iopt,nvect,izzz,dxm,dele,drms
   common/deriv/h(memdrv),ih(maxint)
   common/runhed/ihead(20),ihead1(20)
   common/belly/ibelly,natbel,igroup(maxatom),isymbl(maxatom), &
         itree(maxatom),igres(600)
   common/ntrun5/nrvec,nrat,iat,jat,imov
   common/mass/atmas(maxatom)
   common/coords/c(3*maxatom),c1(3*maxatom)
   common/tcor/ tmax, tintvl, intname(maxint)
   common/ired/ vtemp, iorder, npairs, langevin
   
   !     ----- READ THE NECESSARY DATA TO INITIATE THE RUN -----
   
   call nmdfil2
   
   call mdread(x)
   
   !     ----- EVALUATE SOME CONSTANTS -----
   
   ncart =  3*natbel
   nrcar = 3*nrat
   if (ntrun < -1 .or. ntrun > 10) then
      write(6,70) ntrun
      call mexit(6,1)
   end if
   if ((ntrun == 9 .or. ntrun == 10) .and. &
         (iorder < 0 .or. iorder > 2)) then
      write(6,75) iorder
      call mexit(6,1)
   end if
   if (ntrun == -1) then
      nint=maxint
   else if (ntrun == 0)  then
      nint = nbonh + nbona
   else if (ntrun == 1) then
      nint = nbonh + nbona + ntheth + ntheta + nphih + nphia
   else if (ntrun == 2) then
      nint = ntheth + ntheta + nphih + nphia
   else if (ntrun == 3) then
      nint =  nphia + nphih
   else if (ntrun == 4) then
      nint = nrgrp*6
   else
      nint = 0
   end if
   
   !     ----- PARTITION THE MEMORY -----
   
   !     ----- LOCATIONS OF DIFFERENT ARRAYS -----
   
   !           B      ... H(I10)
   !           iat    ... H(I20)
   !           jat    ... H(I30)
   !           kat    ... H(I40)
   !           lat    ... H(I50)
   !           fint   ... H(I60)
   !           avint  ... H(I65)
   !           vect   ... H(I70)
   !           freq   ... h(i80)
   !           rvec   ... h(i90)
   !           rfreq  ... h(i95)
   !           --- partitions below are used for
   !           ---   (ired) correlation functions
   !           cftmp  ... h(i100) (real array of complex numbers)
   !           p2cftmp... h(i102) (real array of complex numbers, only for ntrun=10)
   !           rcftmp ... h(i104) ('')
   !           cf     ... h(i110)
   !           p2cf   ... h(i112) (only for ntrun=10)
   !           rcf    ... h(i114) ('')
   !           data   ... h(i120) (real array of complex numbers)
   !           table  ... h(i125)
   !           rave   ... h(i130) (only for ntrun=10)
   !           r3iave ... h(i140) ('')
   !           r6iave ... h(i150) ('')
   !           avcrd  ... h(i160) ('')
   !           --- integer heap
   !           pair1  ...ih(j10)
   !           pair2  ...ih(j20)
   
   l1 = ncart*nint
   i10 = 1
   i20 = i10 + l1
   i30 = i20 + nint
   i40 = i30 + nint
   i50 = i40 + nint
   i60 = i50 + nint
   i65 = i60 + nint
   i70 = i65 + nint
   if( langevin == 0 ) then
      i80 = i70 + (ncart*nvect)
      i90 = i80 + nvect
   else
      i80 = i70
      i90 = i80
   end if
   i95 = i90 + (nrcar*nrvec)
   iwrkspc = i95 + nrvec
   if(ntrun == 9 .or. ntrun == 10) then
      i100 = iwrkspc
      ndata = 2**(int(log(dble(4 * int((last-first+1)/iskip) - 1)) / &
            log(2.0d0)) + 1)
      if(ntrun == 9) then
         i102 = i100
         i104 = i100
         i110 = i100 + 2 * &
               int((last - first + 1) / iskip) * &
               (2 * iorder + 1) * &
               (iend - ibeg + 1)
         i112 = i110
         i114 = i110
         i120 = i110 + (iend - ibeg + 1) * &
               (int(tmax/tintvl) + 1)
         i125 = i120 + ndata
         i130 = i125 + 2 * ndata + 15
         i140 = i130
         i150 = i130
         i160 = i130
         iwrkspc = i130
      else if(ntrun == 10) then
         i102 = i100 + 2 * &
               int((last - first + 1) / iskip) * &
               (2 * iorder + 1) * npairs
         i104 = i102 + 2 * &
               int((last - first + 1) / iskip) * &
               (2 * iorder + 1) * npairs
         i110 = i104 + 2 * &
               int((last - first + 1) / iskip) * &
               npairs
         i112 = i110 + npairs * (int(tmax/tintvl) + 1)
         i114 = i112 + npairs * (int(tmax/tintvl) + 1)
         i120 = i114 + npairs * (int(tmax/tintvl) + 1)
         i125 = i120 + ndata
         i130 = i125 + 2 * ndata + 15
         i140 = i130 + npairs
         i150 = i140 + npairs
         i160 = i150 + npairs
         iwrkspc = i160 + 3 * npairs
      end if  ! (ntrun == 9)
   end if  ! (ntrun == 9 .or. ntrun == 10)
   if(iwrkspc > memdrv) then
      write(6,60) iwrkspc
      call mexit(6, 1)
   else
      write(6,90) iwrkspc
   end if
   
   if(ntrun == 9 .or. ntrun == 10) then
      j10 = 1
      j20 = j10 + npairs
      iwrkspc = j20 + npairs
      if(iwrkspc > maxint) then
         write(6,65) iwrkspc
         call mexit(6, 1)
      else
         write(6,95) iwrkspc
      end if
   end if
   
   if (ntrun == -1) then
      nbonh  = 0
      nbona  = 0
      ntheth = 0
      ntheta = 0
      nphih  = 0
      nphia  = 0
      kb = 0
      10 continue
      it = 0
      read(4,*,end=20) type, ia, ja, ka, la
      kb = kb + 1
      intname(kb) = type
      if (ka <= 0) then
         nbonh = nbonh + 1
         ibh(nbonh) = 3*(ia-1)
         jbh(nbonh) = 3*(ja-1)
         icbh(nbonh) = it
      else if (la <= 0) then
         ntheth = ntheth + 1
         ith(ntheth) = 3*(ia-1)
         jth(ntheth) = 3*(ja-1)
         kth(ntheth) = 3*(ka-1)
         icth(ntheth) = it
      else
         nphih = nphih + 1
         iph(nphih) = 3*(ia-1)
         jph(nphih) = 3*(ja-1)
         kph(nphih) = 3*(ka-1)
         lph(nphih) = 3*(la-1)
         icph(nphih) = it
      end if
      goto 10
      20 continue
      nint = nbonh + ntheth + nphih
      if (nint > maxint) then
         write (6,*) ' Not dimensioned to handle so many time'
         write (6,*) '        correlation functions'
         call mexit(6, 1)
      end if
   else if(ntrun /= 9 .and. ntrun /= 10) then
      tmax = 0.0
      tintvl = 1.0
   end if  ! (ntrun == -1)
   
   if(ntrun /= 10) then
      if (ivform == 0) then
         call amopen(50,vecs,'O','U','R')
      else
         call amopen(50,vecs,'O','F','R')
      end if
      if(ntrun == 9) then
         call rdvect(h(i70),h(i80),nvect,npairs,50,x,ieff,ivform,langevin)
      else
         call rdvect(h(i70),h(i80),nvect,ncart,50,x,ieff,ivform,langevin)
      end if
      if (iend > nvect) iend = nvect
   end if
   if (ntrun == 5) then
      do 30 i=1,3*natom
         c(i) = x(i)
      30 continue
      ncary = ncart
      if (ivform == 0) then
         call amopen(52,rvecs,'O','U','R')
      else
         call amopen(52,rvecs,'O','F','R')
      end if
      call rdvect(h(i90),h(i95),nrvec,nrcar,52,c1,ieff,ivform,langevin)
      call excise(h(i70),ncart,nvect)
   end if
   call setvar
   
   !   --- here calculate the b - matrix and analyze the normal modes
   !       bmat   -- b matrix in internal coord.
   !       bgroup -- b matrix in rigid groups
   
   if(ntrun < 4) then
      call bmat(h(i10),h(i20),h(i30),h(i40),h(i50),h(i60), &
            nint,ncart,xbel,h(i65))
      call proj(h(i10),h(i20),h(i30),h(i40),h(i50),h(i60), &
            h(i70),h(i80),nint,ncart,nvect,h(i65))
   end if
   if(ntrun == 4) then
      call group(h(i10),nrgrp,nint,ncart)
      call grpprj (h(i10),h(i70),h(i80),nint,ncart,nvect)
   end if
   if(ntrun == 5) then
      call dotref(h(i70),h(i90),nvect,nrcar,ncart,ibeg, &
            iend,ncary)
   end if
   if(ntrun == 6) call rmsf(h(i70),h(i80),nvect,ncart)
   if (ntrun == 7) then
      
      !     --- read in atom pairs and compute dipole-dipole
      !         correlation functions:
      
      write(6,'(a,l1,a,f6.2)') 'Dipole-dipole correction factors: bose=', &
         bose,', vtemp=',vtemp
      write(6,*) '---------------------------------------------------'
      write(6,*) '   i    j   mode    freq    1-S**2    P2(cum)'
      40 read(4,*,end=50) kup,lup
      if( langevin == 0 ) then
         call corf(ncart,nvect,h(i80),h(i70),x,kup,lup,bose,ibeg,iend,vtemp)
      else
         nat6 = 2*ncart
         np = int(tmax/tintvl) + 1
         call corfl(ncart, nat6, nvect, kup, lup, 2, tmax, np, bose, x)
      end if
      goto 40
   end if
   if (ntrun == 8) call quasip(h(i70),x,c1,ncart)
   if (ntrun == 9 .or. ntrun == 10) then
      
      !     --- read in atom pairs
      
      npai = 0
      45 read(4,*,end=47) ih(j10+npai), ih(j20+npai)
      npai = npai + 1
      if(npai > npairs) then
         write(6,80) npai, npairs
         call mexit(6,1)
      end if
      goto 45
      47 continue
      
      !     --- calc correlation functions
      
      call corffft(ntrun, nvect,h(i70),ih(j10),ih(j20), &
            h(i100),h(i102),h(i104), &
            h(i110),h(i112),h(i114), &
            ndata,h(i120),h(i125), &
            h(i130),h(i140),h(i150),h(i160))
   end if
   
   50 call mexit(6, 0)
   60 format('.... not enough memory for H...need:',i10)
   65 format('.... not enough memory for IH...need:',i10)
   70 format('.... ntrun ( = ',i5,') is not in range...stop')
   75 format('.... iorder ( = ',i5,') is not in range...stop')
   80 format('.... number of pairs read (=',i6, &
         ') exceeds npairs (=',i6,')...stop')
   90 format('.... H core used =',i10)
   95 format('.... IH core used =',i10)
end program nmanal 
