
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
!+ [Enter a one-line description of subroutine bresna here]
subroutine bresna(natom,nres,npair,nhb,ipres,iac,ico,iblo,inb, &
      iar1,iar2,x,igres,ntb,cut,ntypes,nbmax)
   implicit double precision (a-h,o-z)
   logical skip
   logical frozi,frozj
   
   dimension iar1(*),iar2(*),iac(*),ico(*)
   dimension iblo(*),inb(*),ipres(*),x(*),igres(*)
   parameter (maxres=2500)
   dimension ifj(100),skip(maxres)
   
   !     ----- calculate nb pairs -----
   
   lhb = 0
   lpair = 0
   ji = 0
   if (nres > maxres) then
      write(6,*) 'too many residues! maximum is ',maxres
      write(6,*) 'edit bresna.f'
      call mexit(6, 1)
   end if
   
   do 300 ii = 1,nres
      ip1 = ipres(ii)
      ip2 = ipres(ii+1)-1
      frozi = igres(ii) <= 0
      
      !       ----- generate the nb residues to be included for the
      !             current residue ii -----
      
      do 280 jj = ii,nres
         skip(jj) = .false.
         frozj = igres(jj) <= 0
         if(frozi .or. frozj) skip(jj) = .true.
         if(skip(jj)) goto 280
         
         jp1 = ipres(jj)
         jp2 = ipres(jj+1)-1
         call resdis(ip1,ip2,jp1,jp2,x,val)
         if(val > cut) skip(jj) = .true.
      280 continue
      
      do 200 i = ip1,ip2
         lpr = 0
         nx = iblo(i)
         
         do 220 j = 1,nx
            ifj(j) = inb(ji+j)
         220 continue
         
         ifj(nx+1) = 0
         ji = ji+nx
         lj = 1
         ityi = iac(i)
         
         !         ----- if iac = 0 then exit the loop and pick up next i -----
         
         if(ityi == 0) goto 180
         iaci = ntypes*(ityi-1)
         
         !         ----- inner do loop to pick up atom j -----
         
         do 210 jj = ii,nres
            jp1 = ipres(jj)
            jp2 = ipres(jj+1)-1
            
            !           ----- skip the residue if not within cutoff -----
            
            if(ii == jj) jp1 = i+1
            
            do 160 j = jp1,jp2
               
               if(j /= ifj(lj)) goto 140
               lj = lj+1
               goto 160
               
               140 continue
               if(skip(jj)) goto 160
               iacj = iac(j)
               if(iacj == 0) goto 160
               index = iaci+iacj
               ic = ico(index)
               
               !             ----- check whether the pair is a nb or hb -----
               
               if(ic < 0) lhb = lhb+1
               lpr = lpr+1
               lpair = lpair+1
               if( lpair >= nbmax ) then
                  write(6,*) 'too many pairs:', nbmax, i
                  call mexit(6,1)
               end if
               iar2(lpair) = j
            160 continue
         210 continue
         180 continue
         iar1(i) = lpr
      200 continue
   300 continue
   nhb = lhb
   npair = lpair
   return
end subroutine bresna 
