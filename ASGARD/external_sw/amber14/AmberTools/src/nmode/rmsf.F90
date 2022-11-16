
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
!+ [Enter a one-line description of subroutine rmsf here]
subroutine rmsf(vect,freq,nvect,n3)
   
   !     ---calculates rms fluctuations for each atom---
   
   implicit double precision(a-h,o-z)
#  include "sizes2.h"
#  include "anal.h"
#  include "infoa.h"
   common/belly/ibelly,natbel,igroup(maxatom),isymbl(maxatom), &
         itree(maxatom),igres(600)
   dimension vect(n3,nvect),freq(*),igraph(maxatom)
   dimension rms(maxatom),rmsx(maxatom),rmsy(maxatom),rmsz(maxatom)
   parameter (consq = 2.39805e-3)
   !                  ----CONSQ = hc/2kT in cm, with T=300K
   !                       (use for quantum, Bose statistics)
   
   !       these inline functions were substituted in below for
   !       the Decstation compiler
   !     vx(i,j) = vect(3*(i-1)+1,j)
   !     vy(i,j) = vect(3*(i-1)+2,j)
   !     vz(i,j) = vect(3*(i-1)+3,j)
   
   call setgrp(natom,igroup,igrph,igraph)
   tkbc2 = 0.46105e-34    !  kT/c**2 in cgs units (grams)
   avo = 6.023e23         !  converts amu to grams
   cnst = tkbc2*avo
   cmtoa = 1.000e08       !  cm to Angstroms
   twopi = 6.2832         !  convert Hz to angular frequency, omega
   cont  = cmtoa/twopi
   
   !     ---loop over all atoms----
   
   do 500 i=1,natbel
      sumx = 0.0
      sumy = 0.0
      sumz = 0.0
      nvec = 0
      do 300 j=1,nvect
         
         !         --- here loop over the eigenvectors - skip all those associated
         !              with zero or negative eigenvalues ---
         
         if(freq(j) >= 0.5)then
            nvec = nvec + 1
            distx = vect(3*(i-1)+1,j) **2
            disty = vect(3*(i-1)+2,j) **2
            distz = vect(3*(i-1)+3,j) **2
            fre = freq(j)*freq(j)
            if (bose) then
               argq = consq*freq(j)
               fre = fre*tanh(argq)/argq
            end if
            sumx = sumx + (distx/fre)
            sumy = sumy + (disty/fre)
            sumz = sumz + (distz/fre)
         end if
      300 continue
      sumx = cnst*sumx
      sumy = cnst*sumy
      sumz = cnst*sumz
      
      rmsx(i) = (sqrt(sumx))*cont
      rmsy(i) = (sqrt(sumy))*cont
      rmsz(i) = (sqrt(sumz))*cont
      rms(i) = (sqrt(sumx + sumy + sumz))*cont
   500 continue
   write(6,155) nvec
   write(6,156)
   do 350 i=1,natom
      write(6,157) igraph(i),i,rms(i),rmsx(i),rmsy(i),rmsz(i)
   350 continue
   return
   155 format(10x,'RMS FLUCTUATIONS BASED ON ',i5,' NORMAL MODES')
   156 format(1x,'atom number',4x,'rms fluctuation',3x,'x-component', &
         3x,'y-component',3x,'z-component')
   157 format(1x,a4,'(',i5,')',7x,4(f12.5,2x))
end subroutine rmsf 
