
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

!---------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine pload here]
subroutine pload(natom,ntypes,iac,ico)
   
   !     EXPANDS PARAMETER POINTER ICO() INTO 2 DIM ARRAY IPARMP(,)
   !     INDICES OF IPARMP(,) ARE THE ATOM TYPES.  VALUES OF IPARMP(,)
   !     POINT TO VDW A AND B PARAMETERS IN CN1() AND CN2() IF POSITIVE,
   !     IF THEY ARE NEGATIVE, THE ABSOLUTE VALUES POINT TO HBOND
   !     PARAMETERS ASOL() AND BSOL().
   !          THIS ROUTINE MUST BE CALLED ONCE AT THE BEGINNING OF ANY
   !     AMBER/CRAY RUN IN WHICH THE VECTOR NONBONDED ROUTINE IS USED.
   
   common/crayon/iparmp(50,50)

   dimension iac(*),ico(*)
   
   !     CALCULATE INDEX OF ICO ARRAY (IX(I06) IS STARTING ADDRESS)
   !     AND LOAD VALUES IN IPARMP(,)
   
   do 200 i = 1, ntypes
      jaci = ntypes * (i - 1)
      do 150 j = 1, natom
         index = jaci + iac(j)
         iparmp(i,iac(j)) = ico(index)
      150 continue
   200 continue
   
   return
end subroutine pload 
