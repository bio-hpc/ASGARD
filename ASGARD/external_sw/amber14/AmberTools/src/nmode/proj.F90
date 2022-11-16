
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
!+ [Enter a one-line description of subroutine proj here]
subroutine proj(b,iat,jat,kat,lat,fint,vect,freq, &
      nint,ncart,nvect,avint)
   implicit double precision (a-h,o-z)
   
#  include "sizes2.h"
#  include "anal.h"
#  include "infoa.h"
   
   common/belly/ibelly,natbel,igroup(maxatom),isymbl(maxatom), &
         itree(maxatom),igres(600)
   dimension b(ncart,nint),iat(*),jat(*),kat(*),lat(*)
   dimension iat1(maxint),jat1(maxint),kat1(maxint),lat1(maxint)
   dimension fint(*),vect(ncart,nvect),freq(*),avint(*)
   dimension prjct(maxint,maxvec),prjt1(maxint),igraph(maxatom)
   dimension ped(maxint),ped1(maxint)
   
   if (iend > maxvec) then
      write(6,*) 'Too many vectors for subroutine proj'
      write(6,*) '  -- change MAXVEC to be at least ',iend
      call mexit(6, 1)
   end if
   if (nint > maxint) then
      write(6,*) 'Too many internal coordinates for subroutine proj'
      write(6,*) '  -- change MXINT to be at least ',nint
      call mexit(6, 1)
   end if
   
   ! -- repack the graph name of atoms
   
   call setgrp(natom,igroup,igrph,igraph)
   
   ! --- here start big loop over all eigenvectors
   
   do 100  ivec= ibeg,iend
      do 10 izero = 1,nint
         prjct(izero,ivec)  = 0.00d0
         prjt1(izero) = 0.00d0
         ped(izero)    = 0.00d0
         ped1(izero)   = 0.00d0
      10 continue
      sped = 0.0d0
      
      ! --- calculate the projection
      
      do 30 l=1,nint
         do 20 k = 1,ncart
         20 prjct(l,ivec) = prjct(l,ivec) + vect(k,ivec)*b(k,l)
         ped(l) = fint(l)*prjct(l,ivec)*prjct(l,ivec)
         sped = sped + ped(l)
      30 continue
      do 40 ln = 1,nint
         if(sped /= 0.0) ped(ln) = ped(ln)/sped
      40 continue
      if(ipro == 0) goto 90
      
      ! --- check for cut off
      
      jkl = 0
      write(6,120) ivec,freq(ivec)
      write(6,110)
      do 50 i=1,nint
         if (ped(i) >= pcut) then
            !          if (abs(prjct(i,ivec)).ge.pcut) then
            jkl=jkl+1
            iat1(jkl)=iat(i)
            jat1(jkl)=jat(i)
            kat1(jkl)=kat(i)
            lat1(jkl)=lat(i)
            ped1(jkl)= ped(i)
            prjt1(jkl)=prjct(i,ivec)*57.296
         end if
      50 continue
      
      ! -- selection sort
      
      if(jkl == 0) goto 90
      jkk = jkl - 1
      do 70 j = 1,jkk
         l =  j
         jj = j + 1
         do 60 i= jj,jkl
            !           if (abs(prjt1(l)).gt.abs(prjt1(i))) go to 210
            if (ped1(l) > ped1(i)) goto 60
            l = i
         60 continue
         tped = ped1(l)
         ped1(l) = ped1(j)
         ped1(j) = tped
         t = prjt1(l)
         prjt1(l) = prjt1(j)
         prjt1(j) = t
         i1 = iat1(l)
         iat1(l) = iat1(j)
         iat1(j) = i1
         j1 = jat1(l)
         jat1(l) = jat1(j)
         jat1(j) = j1
         k1 = kat1(l)
         kat1(l) = kat1(j)
         kat1(j) = k1
         l1 = lat1(l)
         lat1(l) = lat1(j)
         lat1(j) = l1
      70 continue
      
      imaxm = min(jkl,10)
      do 80 imax = 1 , imaxm
         !         if(ped1(imax).eq.0.0000d0) go to 230
         if(lat1(imax) /= 0) then
            write(6,130)  igraph(iat1(imax)),iat1(imax), &
                  igraph(jat1(imax)),jat1(imax), &
                  igraph(kat1(imax)),kat1(imax), &
                  igraph(lat1(imax)),lat1(imax), &
                  prjt1(imax),ped1(imax)
            goto 80
         end if
         if(kat1(imax) /= 0) then
            write(6,140)  igraph(iat1(imax)),iat1(imax), &
                  igraph(jat1(imax)),jat1(imax), &
                  igraph(kat1(imax)),kat1(imax), &
                  prjt1(imax),ped1(imax)
            goto 80
         end if
         write(6,150)  igraph(iat1(imax)),iat1(imax), &
               igraph(jat1(imax)),jat1(imax), &
               prjt1(imax)/57.296,ped1(imax)
         
      80 continue
      90 continue
   100 continue
   if(ifluc >= 1) call fluct(prjct,iat,jat,kat,lat,nint,igraph, &
         freq,fint,avint)
   return
   
   110 format(58x,'  proj.      PE%')
   120 format(/ /' NORMAL MODE ',i3,'  WITH FREQUENCY =',f12.4/)
   130 format(' dih.   ',3(a5,'(',i3,') - '),a5,'(',i3,')',2f10.5)
   140 format(' angle  ',2(a5,'(',i3,') - '),a5,'(',i3,')',13x,2f10.5)
   150 format(' bond   ',a5,'(',i3,') - ',a5,'(',i3,')',26x,2f10.5)
end subroutine proj 
