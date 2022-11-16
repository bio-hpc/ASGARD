
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
!+ [Enter a one-line description of subroutine setgam here]
subroutine setgam (natom, a, amass, exp, eta, gamma, x, &
      ioseen, hrmax,igraph)
   
   !     ----- to set up the friction matrix gamma
   
   implicit double precision (a-h,o-z)
   parameter (maxat=3000)
   dimension expr(maxat)
   character(len=4) label
   dimension a(6*natom,6*natom), amass(natom), exp(natom), gamma(*)
   dimension x(3*natom)
   character(len=4) igraph(natom)
#  include "files.h"
   namelist /expos/ expr
#if UXP || IBM3090
   logical there
#endif
   nr3 = 3*natom
   nr6 = 6*natom
   
   !     ----- set up hydrodynamic radii
   
   if (natom > maxat) then
      write(6,*) 'MAXAT is too small in subroutine setgam!'
      call mexit(6, 1)
   end if
   do 10 i = 1, natom
      expr(i) = 0.0
   10 continue
   
   !     err=40 prevents amopen() from being used here
#if UXP || IBM3090
   inquire (file = expfil, exist = there)
   if (.not. there) goto 40
#endif
   open (10, file = expfil, status = 'old', err=40)
#ifdef NMLEQ
   read (10,nml=expos)
#else
   read (10, expos)
#endif
   close (10)
   do 11 i=1,natom
      exp(i) = expr(i)
   11 continue
   
   big = exp(1)
   do 20 i = 2, natom
      if (exp(i) > big) big = exp(i)
   20 continue
   
   factor = hrmax / sqrt(big)
   do 30 i = 1, natom
      exp(i) = sqrt(exp(i)) * factor
      write(label,'(a4)') igraph(i)
      if (label(1:1) == 'H') exp(i) = min(0.2d0,exp(i))
   30 continue
   goto 50
   
   40 continue
   write (6,*)
   write (6,'(t1,a)') ' expfile is not present. Hydrodynamic radius'
   write (6,'(t1,a,f6.4,a)') &
         '   of ',hrmax,' angs. is used for all atoms'
   write (6,'(t1,a)') '   except hydrogen, where max. value is 0.2'
   write (6,*)
   do 45 i = 1, natom
      exp(i) = hrmax
      write(label,'(a4)') igraph(i)
      if (label(1:1) == 'H') exp(i) = min(0.2d0,exp(i))
   45 continue
   
   50 continue
   
   !     ----- set up friction matrix gamma
   
   if (ioseen == 0) then
      
      !       ----- diagonal gamma
      
      sxpita = 6.0 * 3.14159 * eta
      ij = 0
      do 70 j = 1, 3 * natom
         do 60 i = 1, j-1
            ij = ij + 1
            gamma(ij) =  0.0
         60 continue
         ij = ij + 1
         jat = (j-1)/3 + 1
         gamma(ij) = sxpita * exp(jat) / amass(jat)
      70 continue
      
   else
      
      !       ----- gamma with hydrodynamic interaction
      
      call oseen (ioseen, natom, nr3, nr6, a, eta, exp, x, &
            gamma, amass)
      do 90 j = 1, nr6
         do 80 i = 1, nr6
            a(i,j) = 0.0
         80 continue
      90 continue
      
   end if
   
   !     ----- copy gamma into matrix a
   
   do 110 j = 1, nr3
      do 100 i = 1, j
         index = j*(j-1)/2 + i
         a(i+nr3,j+nr3) = - gamma(index)
         a(j+nr3,i+nr3) = - gamma(index)
      100 continue
   110 continue
   
   return
end subroutine setgam 
