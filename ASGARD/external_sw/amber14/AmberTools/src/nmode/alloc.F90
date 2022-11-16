
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
!+ [Enter a one-line description of subroutine alloc here]
subroutine alloc(maxdup, memusd_x, memusd_i, memusd_h)
   
   !     ----- to allocate memory.
   
   implicit double precision (a-h, o-z)
   
#  include "pointer.h"
#  include "inpdat.h"
   
   ntype  = ntypes * ntypes
   nttyp  = ntypes * (ntypes+1) / 2
   
   !     ----- pointers to real arrays
   
   mchrg  = 1
   mamass = mchrg  + natom
   mrk    = mamass + natom
   mreq   = mrk    + numbnd
   mtk    = mreq   + numbnd
   mteq   = mtk    + numang
   mfk    = mteq   + numang
   mfpk   = mfk    + numang
   mqeq   = mfpk   + numang
   mpk    = mqeq   + numang
   mpn    = mpk    + nptra
   mphase = mpn    + nptra
   msolty = mphase + nptra
   mcn1   = msolty + natyp
   mcn2   = mcn1   + nttyp
#ifdef CHARMM
   mcn114 = mcn2   + nttyp
   mcn214 = mcn114 + nttyp
   masol  = mcn214 + nttyp
#else
   masol  = mcn2   + nttyp
#endif
   mbsol  = masol  + nphb
   mhbcut = mbsol  + nphb

   mxref  = mhbcut + nphb
   mwref  = mxref  + nr3
   msf    = mwref  + natom
   mpol   = msf    + nr3
   mx     = mpol   + natom

   mf     = mx     + nr3
   mwinv  = mf     + nr3
   mh     = mwinv  + nr3

   if(ntrun == 1) then
      memusd_x  = mh     + (ns3*(ns3+1)/2)
   end if
   
   if(ntrun == 2) then
      !                                 ---space for tstate searches
      mdd    = mh     + (ns3*(ns3+1)/2)
      mvect  = mdd    + (ns3*(ns3+1)/2)
      mroots = mvect  + ns3*ivect
      mb     = mroots + ns3
      mxdir  = mb     + 9*ns3
      mxinit = mxdir  + ns3
      memusd_x = mxinit + ns3

   end if
   
   if(ntrun == 3) then
      !                                 ---no 2nd deriv; scratch for zxcgr
      mcscr = mh
      memusd_x = mcscr   + 6*ns3

   end if
   
   if(ntrun == 4) then
      if(ndiag /= 0) then
         !                                 ---extra space for routine findl
         mxdir  = mh     + (ns3*(ns3+1)/2)
         mdd    = mxdir  + ns3
         mb     = mdd    + (ns3*(ns3+1)/2)
         mroots = mb     + 9*ns3
         mvect  = mroots + ns3
         memusd_x  = mvect  + ns3
      else
         !                                 ---just space for xdir
         mxdir  = mh     + (ns3*(ns3+1)/2)
         memusd_x  = mxdir  + ns3
      end if
   end if
   if(ntrun == 5) then
      !                                --- space for langevin modes
      !                                --- need 2nd der and some more
      nat6   = 6 * natom
      nat3   = 3 * natom
      ma     = mh     + (ns3*(ns3+1)/2)
      mwr    = ma     + nat6*nat6
      mwi    = mwr    + nat6
      mz     = mwi    + nat6
      mfv1   = mz     + nat6*nat6
      mhrad  = mfv1   + nat6
      mgam   = mhrad  + natom
      memusd_x = mgam   + nat3*(nat3+1)/2
   end if

   mxbel  = memusd_x  + natom   + 1
   memusd_x  = mxbel  + ns3

   if(ntrun == 1) then
      mcscr  = mxbel
      if(mod(mcscr,2) /= 1) mcscr = mcscr + 1
      mcval  = mcscr  + ns3*9
      mcvec  = mcval  + ns3
      if( ismem == 0 ) then
         memusd_x = mcvec  + ns3*ns3
      else
         memusd_x = mcvec  + ns3*nvect
      end if
   end if
   
   !=============================================================================
   
   !       ----- begin: pointers to integer arrays

   miac   = 1
   miblo  = miac   + natom   + 1
   mico   = miblo  + natom   + 1
   mipres = mico   + ntype   + 1

   !       ----- The following three arrays (i.e. ibh, jbh and icbh)
   !       ----- should always appear together with no other
   !       ----- array in between them, because rdparm considers
   !       ----- them as one array

   mibh   = mipres + nres    + 2
   mjbh   = mibh   + nbonh   + 1
   micbh  = mjbh   + nbonh   + 1

   !       ----- Another pack of three arrys

   miba   = micbh  + nbonh   + 1
   mjba   = miba   + nbona   + 1
   micba  = mjba   + nbona   + 1

   !       ----- Now a pack of four arrays

   mith   = micba  + nbona   + 1
   mjth   = mith   + ntheth  + 1
   mkth   = mjth   + ntheth  + 1
   micth  = mkth   + ntheth  + 1

   !       -----     ... and so on

   mita   = micth  + ntheth  + 1
   mjta   = mita   + ntheta  + 1
   mkta   = mjta   + ntheta  + 1
   micta  = mkta   + ntheta  + 1

   miph   = micta  + ntheta  + 1
   mjph   = miph   + (nphia+nphih+maxdup)+ 1
   mkph   = mjph   + (nphia+nphih+maxdup)+ 1
   mlph   = mkph   + (nphia+nphih+maxdup)+ 1
   micph  = mlph   + (nphia+nphih+maxdup)+ 1

   mipa   = micph  + (nphia+nphih+maxdup)+ 1
   mjpa   = mipa   + (nphia+nphih+maxdup)+ 1
   mkpa   = mjpa   + (nphia+nphih+maxdup)+ 1
   mlpa   = mkpa   + (nphia+nphih+maxdup)+ 1
   micpa  = mlpa   + (nphia+nphih+maxdup)+ 1

   !       ----- Arrays being together no longer important

   minb   = micpa  + (nphia+nphih+maxdup)+ 1
   mjoin  = minb   + nnb     + 1
   mrotat = mjoin  + natom   + 1
   mgroup = mrotat + natom   + 1
   migrp2 = mgroup + natom   + 1
   migres = migrp2 + natom   + 1
   miar1  = migres + natom   + 1
   memusd_i    = miar1  + natom   + 1
   
   ! end of "regular" integer arrays
   
   if(ntrun == 2) then

      mkpvt  =  memusd_i
      memusd_i = mkpvt + ns3

   end if

   if(ntrun == 5) then
      !                                --- space for langevin modes
      !                                --- need 2nd der and some more
      nat6   = 6 * natom
      nat3   = 3 * natom
      miv1 = memusd_i
      memusd_i = miv1   + nat6
   end if
   
   mnbel  = memusd_i
   miar2 = mnbel  + natom   + 1
   memusd_i = miar2 + natom*(natom-1)/2  ! miar2 begins the pairlist
   nbmax  = (memusd_i - miar2 - 1)
   
   !===========================================================================
   
   !       ----- begin: pointers to character arrays

   mgraph = 1
   mlbres = mgraph + natom + 1
   msymbl = mlbres + nres + 1
   mitree = msymbl + natom + 1
   memusd_h = mitree + natom + 1


#if 0
        write (6,1)
    1   format ('   Memory allocation : ',/)
        write (6,2)    mchrg,  mamass,  mrk,    mreq,   mtk,    &
                       mteq,   mfk,     mfpk,   mqeq,   mpk,    &
                       mpn,    mphase,  msolty, mcn1,   mcn2,    &
                       masol,  mbsol,   mhbcut, mxref,  msf,    &
                       miac,   miblo,   mico,   mipres, mibh,  &
                       mjbh,   micbh,   miba,   mjba,   micba, &
                       mith,   mjth,    &
                       mkth,   micth,   mita,   mjta,   mkta,    &
                       micta,  miph,    mjph,   mkph,   mlph,    &
                       micph,  mipa,    mjpa,   mkpa,   mlpa,    &
                       micpa,  minb,    mgroup,    &
                       migres, miar1,   mx,     mf,     mh,    &
                       mnbel,  mxbel,   miar2,  mcscr,  mcval,    &
                       mcvec,  mgam,    mcn114, mcn214, mjoin,    &
                       mrotat, mpol,    mkpvt, mgraph, mlbres,   &
                       msymbl, mitree
   2    format (5i12)
#endif
   
   write (6,3) memusd_x
   3 format (/, 'Total memory required : ', i12,' real words')
   write (6,7) memusd_i
   7 format (/, 'Total memory required : ', i12,' integer words')
   write (6,8) memusd_h
   8 format (/, 'Total memory required : ', i12,' 4-character words')
   write(6,'(/,''Maximum nonbond pairs'', i12)')nbmax
   
   return
end subroutine alloc 
