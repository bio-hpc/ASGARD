
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
!+ [Enter a one-line description of subroutine makel here]
subroutine makel (neig, nm, n, z, wr, wi)
   
   !     ----- eq.(4.1), j. chem. phys., 85, 7334 (1986)
   
   implicit double precision (a-h,o-z)
   complex el, pi, sql, lambda, minusi, zero
   dimension z(nm, n), wr(n), wi(n)
   
   minusi = (0.0, -1.0)
   zero   = (0.0,  0.0)
   
   j = 1
#ifdef ENDDO
   do while (j <= neig)
#else
   10 if (j > neig) goto 60
#endif
      lambda = cmplx ( wr(j), wi(j) )
      sql = sqrt(lambda)
      if (lambda == zero) then
         do 20 i = 1, n/2
            z(i,j) = 0.0
         20 continue
         j = j + 1
      else if (wi(j) == 0.0) then
         if (wr(j) > 0.0) then
            do 30 i = 1, n/2
               pi = z(i,j)
               el = pi*sql
               z(i,j) = real(el)
            30 continue
         else
            do 40 i = 1, n/2
               pi = z(i,j)
               el = pi*sql
               z(i,j) = real(minusi*el)
            40 continue
         end if
         j = j + 1
      else
         do 50 i = 1, n/2
            pi = cmplx ( z(i,j), z(i,j+1) )
            el = pi*sql
            z(i,j) = real(el)
            z(i,j+1) = real(minusi*el)
         50 continue
         j = j + 2
      end if
#ifdef ENDDO
   end do
#else
      goto 10
      60 continue
#endif
   
   return
end subroutine makel 
