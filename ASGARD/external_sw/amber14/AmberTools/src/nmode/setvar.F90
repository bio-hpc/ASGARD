
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
!+ [Enter a one-line description of subroutine setvar here]
subroutine setvar(x,ix,igroup,nbel)
   implicit double precision (a-h,o-z)
   dimension ix(*),igroup(*),nbel(*)
   dimension x(*)
#  include "pointer.h"
#  include "inpdat.h"
   
   !     ----- ROUTINE TO DO THE NECESSARY ACCOMODATIONS FOR PROTEIN
   !           BELLY MINIMISATIONS -----
   
   
   !     ----- DELETE BONDS WHICH ARE IN THE BELLY ALONE -----
   
   if (ibelly /= 0) then
      mbonh = nbonh
      call setbon(nbonh,mbonh,ix(mibh),ix(mjbh),ix(micbh),ix(mgroup))
      call setbon(nbona,mbona,ix(miba),ix(mjba),ix(micba),ix(mgroup))
      
      !     ----- DELETE THE BONDS WHICH ARE IN THE BELLY ALONE -----
      
      mtheth = ntheth
      call setang(ntheth,mtheth,ix(mith),ix(mjth),ix(mkth),ix(micth), &
            ix(mgroup))
      call setang(ntheta,mtheta,ix(mita),ix(mjta),ix(mkta),ix(micta), &
            ix(mgroup))
      
      !     ----- DELETE THE DIHEDRALS -----
      
      mphih = nphih
      call setdih(nphih,mphih,ix(miph),ix(mjph),ix(mkph),ix(mlph), &
            ix(micph),ix(mgroup))
      call setdih(nphia,mphia,ix(mipa),ix(mjpa),ix(mkpa),ix(mlpa), &
            ix(micpa),ix(mgroup))
      
   end if
   
   !     ----- move the active coordinates and masses to the working array
   
   call setcor(natom,ix(mgroup),x(mx),x(mxbel))

   call getm(natom,ix(mgroup),x(mamass))
   
   !     ----- set up NBEL array: NBEL(I) = new pointer array for atom i
   
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
