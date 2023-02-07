
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
!+ [Enter a one-line description of subroutine eprnt here]
subroutine eprnt(ith,nat,f)
   
   !     ----- to print out the potential energies
   
   implicit double precision (a-h, o-z)
   dimension f(*)
   
#  include "epot.h"
   
   gmax = grdmax(f,rms,3*nat)
   write(6,5)
   write(6,10) ith
   write(6,15) etot,gmax,rms
   write(6,20)
   write(6,25) enb,eel,enbh,eb
   write(6,21)
   write(6,25) ea,edih,enb14,eel14
   write(6,23)
   write(6,25) epol,e3bod
   if (econs /= 0.d0) then
      write(6,22)
      write(6,25) econs
   end if
   call amflsh(6)
   
   5 format(/ /,3(' *****************       '))
   10 format(2x,' step = ',i10)
   15 format(5x,'F = ',e14.6,4x,'GRDMAX = ',e12.6,4x,'GNORM = ',e12.6)
   20 format('       E-NONB              E-ELE        ', &
         '       E-HBOND             E-BOND')
   21 format('       E-ANGLE             E-DIHED      ', &
         '       E-NB14              E-EEL14')
   23 format('       E-POL               E-3BOD ')
   22 format('       CONSTRAINT')
   25 format(4(4x,e12.5,4x))
   
   return
end subroutine eprnt 
