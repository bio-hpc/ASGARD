
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
!+ [Enter a one-line description of subroutine spsort here]
subroutine spsort (n, dbl, wr, wi, ind, nreal)
   
   !     ----- a special sort routine called by lmode
   
   implicit double precision (a-h,o-z)
   dimension dbl(n), ind(n), wr(n), wi(n)
   
   !     ----- First make ind in natural order
   
   do 10 i = 1, n
      ind(i) = i
      dbl(i) = abs(wi(i))
   10 continue
   
   !     ----- sort keeping track of which element goes where
   
   do 30 i = 1, n-1
      
      !       ----- get pointer to minimum value of remaining entries
      
      ip = i
      do 20 j = i+1, n
         if (dbl(j) < dbl(ip)) ip = j
      20 continue
      
      !       ----- interchange ip-th and i-th dbl
      
      dt = dbl(i)
      dbl(i) = dbl(ip)
      dbl(ip) = dt
      
      !       ----- interchange ip-th and i-th ind
      
      it = ind(i)
      ind(i) = ind(ip)
      ind(ip) = it
      
   30 continue
   
   !     ---- oops! Bring back complex conjugates in original order
   !     ---- i.e., (a + ib) followed by (a - ib)
   
   i = 1
#ifdef ENDDO
   do while (i <= n)
#else
   40 if (i > n) goto 45
#endif
      it = ind(i)
      if (wi(it) /= 0.0 .and. wr(it) == wr(ind(i+1))) then
         if (wi(it) < 0.0) then
            ind(i) = ind(i+1)
            ind(i+1) = it
         end if
         i = i + 2
      else
         i = i + 1
      end if
#ifdef ENDDO
   end do
#else
      goto 40
      45 continue
#endif
   
   !     ----- count reals
   
   nreal = 0
#ifdef ENDDO
   do while (abs(wi(ind(nreal+1))) < 1.0e-3 .and. nreal < n)
#else
   50 if (.not. (abs(wi(ind(nreal+1))) < 1.0e-3 .and. nreal < n)) &
         goto 55
#endif
      nreal = nreal + 1
      dbl(nreal) = abs(wr(ind(nreal)))
#ifdef ENDDO
   end do
#else
      goto 50
      55 continue
#endif
   
   !     ----- sort reals
   
   do 70 i = 1, nreal-1
      ip = i
      do 60 j = i+1, nreal
         if (dbl(j) < dbl(ip)) ip = j
      60 continue
      dt = dbl(i)
      dbl(i) = dbl(ip)
      dbl(ip) = dt
      it = ind(i)
      ind(i) = ind(ip)
      ind(ip) = it
   70 continue
   
   return
end subroutine spsort 
