
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
!+ [Enter a one-line description of subroutine setgr here]
subroutine setgr (natom,      natc,        nres,        ngroup, &
      ibelly,     natsys,      igres,       igroup, &
      ipres,      lbres,       igraph,      isymbl, &
      itree,      icons,       igrp2,       wref, &
      xref)
   
   character(len=4) lbres, igraph, isymbl, itree
   dimension  igres(*),  igroup(*),  ipres(*),  lbres(*), &
         igraph(*), isymbl(*), itree(*),   igrp2(*), &
         wref(*),   xref(*)
   
   if(ibelly /= 0) then
      call rgroup (natom,     natc,        nres,        ngroup, &
            ipres,     lbres,       igraph,      isymbl, &
            itree,     igroup,      dummy,       dummy, .false., .true., 5)
   else
      do 10 i = 1, natom
         igroup(i) = 1
      10 continue
   end if
   
   !     ----- FIND THE TOTAL NUMBER OF ACTIVE ATOMS AND RESIDUES -----
   
   call setatm(natom,nres,natsys,ipres,igroup,igres)
   
   if (icons /= 0) then
      call rgroup (natom,     natc,        nres,        ngroup, &
            ipres,     lbres,       igraph,      isymbl, &
            itree,     igrp2,       wref,        xref, .true., .false., 5)
   end if
   return
end subroutine setgr 
