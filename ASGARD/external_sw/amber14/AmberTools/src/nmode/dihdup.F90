
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
!+ [Enter a one-line description of subroutine dihdup here]
subroutine dihdup(mphia,nphia,ip,jp,kp,lp,icp,pn,maxdup)

   !     duplicates pointers to multi-term dihedrals for vector ephi.
   !     H-atom diheds are duplicated at the end of the h-atom lists,
   !     but heavy atom diheds must be inserted between the end of the
   !     heavy-atom diheds, but before the constraint diheds if they are
   !     present.  In order to use this code, extra space MUST be allocated
   !     in LOCMEM (parameter MAXDUP set in main)
   
   
   implicit double precision (a-h,o-z)

!     ----- The following are pointers to the first element of
!     ----- various arrays. The pointer name is derived by prefixing
!     ----- array name by m. For example amass begins at x(mamass)
!     ----- Exceptions are the following arrays with six character
!     ----- names:
!     -----          igraph, isymbl, igroup
!     ----- In these cases m replaces the first character

common /pointr/   mchrg,  mamass,  mrk,    mreq,   mtk, &
      mteq,   mfk,     mfpk,   mqeq,   mpk, &
      mpn,    mphase,  msolty, mcn1,   mcn2, &
      masol,  mbsol,   mhbcut, mxref,  msf, &
      momega, mgraph,  miac,   miblo,  mico, &
      mlbres, mipres,  mibh,   mjbh,   micbh, &
      miba,   mjba,    micba,  mith,   mjth, &
      mkth,   micth,   mita,   mjta,   mkta, &
      micta,  miph,    mjph,   mkph,   mlph, &
      micph,  mipa,    mjpa,   mkpa,   mlpa, &
      micpa,  minb,    msymbl, mitree, mgroup, &
      migres, miar1,   mx,     mf,     mh, &
      mnbel,  mxbel,   miar2,  mcscr,  mcval, &
      mcvec,  mdd,     mb,     mroots, mvect, &
      mxdir,  mwref,   migrp2, ma,     mwr, &
      mwi,    mz,      mfv1,   miv1,   mhrad, &
      mgam,   mwinv,   mkpvt,  mxinit, mcn114, &
      mcn214, mjoin,   mrotat, mpol
   
   !     INPUT:
   
   integer mphia
   !        ... number of real dihedrals
   integer nphia
   !        ... number of real dihedrals + number of constraint dihedrals
   dimension ip(*), jp(*), kp(*), lp(*)
   !        ... pointers to atoms of dihedrals
   dimension icp(*)
   !        ... pointers to dihedral parameters
   dimension pn(*)
   !        ... periodicity; negative if multi-term, read until + encountered
   integer MAXDUP
   !        ... max number of duplicated dihedrals

   integer ndup
   logical const
   const = (nphia > mphia)
   ndup = 0
   ierr = 0

   do 200 i = 1, mphia
      ic = icp(i)
      100 continue
      if (pn(ic) < 0.0) then
         ndup = ndup + 1
         if (ndup > maxdup) then
            ierr = 1
         else
            !               --- move constraint dihed to end, if necessary ---
            if (const) then
               ip(nphia+ndup)  = ip(mphia+ndup)
               jp(nphia+ndup)  = jp(mphia+ndup)
               kp(nphia+ndup)  = kp(mphia+ndup)
               lp(nphia+ndup)  = lp(mphia+ndup)
               icp(nphia+ndup) = icp(mphia+ndup)
            end if
            !                 --- duplicate pointers of multi-term dihedrals ---
            ip(mphia+ndup)  = ip(i)
            jp(mphia+ndup)  = jp(i)
            kp(mphia+ndup)  = kp(i)
            lp(mphia+ndup)  = lp(i)
            !                 --- but use the NEXT parameter pointer ---
            icp(mphia+ndup) = ic + 1
            !                 --- check for a third or higher term ---
         end if
         ic = ic + 1
         goto 100
      end if
   200 continue
   
   if (ierr == 1) then
      write(6,'(/,5x,a,i5,a)') 'MAXDUP =',maxdup,' exceeded in dihdup'
      write(6,'(/,5x,a,i5,a)') 'set MAXDUP =',ndup,' in nmode.f'
      call mexit(6, 1)
   end if
   mphia = mphia + ndup
   nphia = nphia + ndup
   write(6,'(/,5x,a,i5,a)') 'Duplicated',ndup,' dihedrals'
   return
end subroutine dihdup 
