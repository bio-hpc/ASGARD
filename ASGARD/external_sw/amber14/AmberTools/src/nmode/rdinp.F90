
!************************************************************************
!                              AMBER                                   **
!                                                                      **
!                  Copyright (c) 1986, 1991, 1995, 1997                **
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
!+ [Enter a one-line description of subroutine rdinp here]
subroutine rdinp
   implicit double precision (a-h, o-z)
   character(len=80) header
   
#  include "inpdat.h"
   
   common/thrbod/iat1(5),iat2(5),acon(5),beta3b(5),gama3b(5) &
         ,i3b(5),n3b,nion
   namelist/data/ntrun,ibelly,maxcyc,drms,nvect,kform, &
         cut,scnb,scee,dielc,idiel,nsave,dfpred, &
         ndiag,smx,emx,alpha,bdwnhl,iprr,iprw, &
         lat,nbel,mbel,irxn,sdstp,stplim,nprint,mindim, &
         icons, eta, ioseen, hrmax, &
         istart,ivect,iflag,isdir,idir,hnot, &
         isw,buphl,idih,ilevel,ivform,ntx,ntxo,t,patm, &
         ipol,i3bod,ismem,idecomp
   
   !   ---- set the default values for namelist variables
   
   ntrun = 1
   ismem = 0
   ibelly = 0
   maxcyc = 100
   drms = 1.0d-5
   nvect = 10
   kform = 1
   cut = 99.0
   scnb = 2.0
   scee = 2.0
   dielc = 1.0
   idiel = 0
   nsave = 20
   dfpred = 0.01
   ndiag = 10
   smx = 0.08
   emx = 0.30
   alpha = 0.8
   bdwnhl = 0.1
   iprr = 0
   iprw = 0
   lat = 0
   irxn = 0
   sdstp = 0.011
   stplim = 0.5
   nprint = 1
   mindim = 5
   icons = 0
   eta = 0.2
   ioseen = 0
   hrmax = 0.77
   istart = 0
   ivect = 2
   isdir = 1
   idir = 1
   hnot = 0.1
   isw = 40
   do 8 i=1,101
      idih(i) = 0
   8 continue
   iflag = 0
   buphl = -0.1
   ipol = 0
   i3bod = 0
   ilevel = 0
   ivform = 1
   ntx = 1
   ntxo = 1
   t = 298.15d+00
   patm = 1.0
   idecomp = 0
   
   write(6,1)
   
   !     ----- READ DATA CHARACTERIZING THE RUN -----
   
   read(5,'(a80)') header
   write(6,'(2x,a80)') header
#ifdef NMLEQ
   read(5,nml=data)
#else
   read(5,data)
#endif
   
   !     the "small memory model" (ismem=1) does not allow for vectors:
   
   if( ismem == 1 ) nvect = 0
   
   eta = eta * 2.945
   !               ... in internal units
   
   !     ----- PRINT DATA CHARACTERIZING THE RUN -----
   
   write(6,42)
   write(6,43) ntrun,maxcyc,ibelly,drms
   write(6,74)
   write(6,76) cut,scnb,scee,dielc,idiel
   write(6,77)
   write(6,78) nsave,dfpred,bdwnhl,smx,emx,alpha,ndiag
   write(6,'(t2,''ipol = '',i2)')ipol
   write(6,'(t2,''i3bod = '',i2)')i3bod
   write(6,'(t2,''idecomp = '',i2)')idecomp
   write(6,'(t2,''nvect = '',i4)')nvect

   
   if (ntrun == 2) then
      write(6,10)
      write(6,20)
      write(6,21) istart,ivect,isdir,idir,isw,iflag
      write(6,44)
      write(6,41) hnot,buphl
   end if
   if (ntrun == 5) write (6, &
         '(/,5x, a6, f5.2,1x,a2, 5x, a9, i1)' &
         ) 'eta = ', eta/2.945, 'cp', 'ioseen = ', ioseen
   if (ntx == 0) write(6,'(5x,a)') &
         'Binary format used for input coords.'
   if (ntxo == 0) write(6,'(5x,a)') &
         'Binary format used for output coords.'

   if(i3bod > 0) then
      write(6,'(''Reading 3-bod information:'')')
      read(5,'(2i5)')n3b,nion
      if(nion == 0 ) then
         write(6,'(t2,''##########################nion = 0'')')
         call exit
      end if
      write(6,'(/t2,''Number of triplets = '',i5,'' nion = '',i3,/, &
            & t2,''sing pair    acons     beta     gamma  '')')n3b,nion
      itemp = n3b
      if(n3b < 0)itemp=1
      do 55 jn = 1,itemp
         read(5,'(a4,a4,2x,3e10.3)')    iat1(jn),iat2(jn) &
               ,acon(jn),beta3b(jn),gama3b(jn)
         write(6,'(t2,a4,1x,a4,1x,3e12.4)')iat1(jn),iat2(jn) &
               ,acon(jn),beta3b(jn),gama3b(jn)
      55 continue
   end if
   
   
   !     check input
   
   inerr = 0

   if (ntrun < 1 .or. ntrun > 5) then
      write(6, '(a,i3,a)') 'NTRUN (',ntrun,') must be in 1..5'
      inerr = 1
   end if
   if (ibelly /= 0 .and. ibelly /= 1) then
      write(6, '(a,i3,a)') 'IBELLY (',ibelly,') must be 0 or 1'
      inerr = 1
   end if
   if (icons /= 0 .and. icons /= 1) then
      write(6, '(a,i3,a)') 'ICONS (',icons,') must be 0 or 1'
      inerr = 1
   end if
   if (idiel /= 0 .and. idiel /= 1) then
      write(6, '(a,i3,a)') 'IDIEL (',idiel,') must be 0 or 1.'
      inerr = 1
   end if
   if (iprr /= 0 .and. iprr /= 1) then
      write(6, '(a,i3,a)') 'IPRR (',iprr,') must be 0 or 1.'
      inerr = 1
   end if
   if (iprw /= 0 .and. iprw /= 1) then
      write(6, '(a,i3,a)') 'IPRW (',iprw,') must be 0 or 1.'
      inerr = 1
   end if
   if (ioseen < 0 .or. ioseen > 2) then
      write(6, '(a,i3,a)') 'IOSEEN (',ioseen,') must be 0, 1 or 2.'
      inerr = 1
   end if
   if (istart /= 0 .and. istart /= 1) then
      write(6, '(a,i3,a)') 'ISTART (',istart,') must be 0 or 1.'
      inerr = 1
   end if
   if (iflag < -1 .or. iflag > 1) then
      write(6, '(a,i3,a)') 'IFLAG (',iflag,') must be 0, 1 or -1.'
      inerr = 1
   end if
   if (idir /= 1 .and. idir /= -1) then
      write(6, '(a,i3,a)') 'IDIR (',idir,') must be 1 or -1.'
      inerr = 1
   end if
   if (ilevel /= 0) write(6,*) 'translational and rotational ', &
         'modes will be shifted up'
   if (ivform == 0) write(6,*) 'unformatted eigenvectors will be', &
         ' written out.'
   if (inerr == 1) call mexit(6,1)
   return
   1 format(/ /10x,55(1h*),/10x,'Initiate the NMODE module of ', &
         'AMBER 8  ',/10x,55(1h*),/ /)
   10 format ( / /,'   Options for transition state search:'/)
   20 format('    istart  ivect  isdir   idir   isw   iflag')
   21 format(1x,7i7)
   41 format(1x,2f10.5)
   42 format('       ntrun    maxcyc    ibelly       drms')
   43 format(' ',3i10,e10.2)
   44 format('       hnot      buphl')
   74 format('       rcut      scnb      scee      dielc   idiel')
   76 format(3x,4f10.5,i6)
   77 format('      nsave   dfpred    bdwnhl     smx       emx    ' &
         ,'  alpha       ndiag ')
   78 format(1x,i10,5f10.5,i10)
end subroutine rdinp 
