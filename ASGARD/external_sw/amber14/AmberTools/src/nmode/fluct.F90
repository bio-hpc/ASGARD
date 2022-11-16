
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
!+ [Enter a one-line description of subroutine fluct here]
subroutine fluct(prjct,iat,jat,kat,lat,nint,igraph, &
      freq,fint,avint)
   implicit double precision(a-h,o-z)
   character(len=8) intname
#  include "sizes2.h"
#  include "anal.h"
   common/conver/ipt(maxint)
   common/tcor/tmax, tintvl, intname(maxint)
   common/minpar/ntrun,maxcyc,ncyc,iopt,nvect,izzz,dxm,dele,drms
   dimension prjct(maxint,*),igraph(*),iat(*),jat(*),kat(*), &
         lat(*),freq(*),rms(maxint),fac(3),fint(*),avint(*)
   
   parameter (const= 7.03398e-13)
   !                  ----CONST = kTN/(2*pi*c)**2 in amu, with T=300K
   !                       (use for classical, Boltzmann stats.)
   parameter (consq = 2.39805e-3)
   !                  ----CONSQ = hc/2kT in cm, with T=300K
   !                       (use for quantum, Bose statistics)
   ratodeg = 180/3.14159
   ibw = iend-ibeg+1
   write(6,90) ibw
   write(6,100)
   
   fac(1) = 1.e8
   fac(2) = ratodeg*1.e8
   fac(3) = ratodeg*1.e8
   
   !     ----- loop over all internal coordinates-----
   
   do 40 i=1,nint
      !                ---comment out (for now) code that does cross-correl.
      !     do 30 k=i,nint
      k = i
      
      t = 0.0
20    continue
      ! do 20 t = 0.0, tmax, tintvl
         sum = 0.0
         do 10 j=ibeg,iend
            prj = prjct(k,j)*prjct(i,j)
            fre = freq(j)*freq(j)
            if (bose) then
               argq = consq*freq(j)
               fre = fre*tanh(argq)/argq
            end if
            sum = sum + (prj/fre) * cos(freq(j)*t)
         10 continue
         
         !    --- get rms fluctuation and convert to Angstroms or degrees:
         
         sum = const*sum
         if (sum >= 0.0) rms(i) = (sqrt(sum))*fac(ipt(i))
         corf = sum*fac(ipt(i))*fac(ipt(i))
         if (t == 0.0) then
            if(k == i) then
               !                         --print diagonal correlations
               if(lat(i) /= 0) then
                  write(6,50)  'dihed:  ',igraph(iat(i)),iat(i), &
                        igraph(jat(i)),jat(i), &
                        igraph(kat(i)),kat(i), &
                        igraph(lat(i)),lat(i), &
                        avint(i),rms(i)
               else if(kat(i) /= 0) then
                  write(6,60)  'angle:  ',igraph(iat(i)),iat(i), &
                        igraph(jat(i)),jat(i), &
                        igraph(kat(i)),kat(i), &
                        avint(i),rms(i)
               else
                  write(6,70)  'bond:   ',igraph(iat(i)),iat(i), &
                        igraph(jat(i)),jat(i), &
                        avint(i),rms(i)
               end if
               
            else
               !                                      ---print cross correlations
               write(6,80) i,k,corf
            end if
            
         else
            !                                       ---print time correlations
            write (6, '(50x,f10.2,8x,e12.5)') t,corf
         end if
         
      t = t + tintvl
      if (t <= tmax ) go to 20
      !20 continue

      !                    ---only do cross correlations when ntrun=-1
      !     if (ntrun.ne.-1) go to 40
      !  30 continue
   40 continue
   return
   
   50 format(a8,3(a4,1x,i4,' - '),a4,1x,i4,2f12.5)
   60 format(a8,2(a4,1x,i4,' - '),a4,1x,i4,12x,2f12.5)
   70 format(a8,a4,1x,i4,' - ',a4,1x,i4,24x,2f12.5)
   80 format(' cross correlation:', 2i5,35x,f12.5)
   90 format(/ /5x,'INTERNAL COORDINATES RMS FLUCTUATIONS BASED ON ', &
         i5,' NORMAL MODES',/ /)
   100 format(55x,'average     rms fluct.'/)
end subroutine fluct 
