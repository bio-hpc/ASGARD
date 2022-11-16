
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
!+ [Enter a one-line description of subroutine getcor here]
subroutine getcor (natom, ntx, x)
   implicit double precision (a-h,o-z)
   character(len=80) line
   
   dimension x(*)
   character(len=4) ititl(20)
   if (ntx == 0) then
      
      !       ----- read title etc.
      
      read (53) ititl
      read (53) natomd
      if (natomd /= natom) then
         write(6,*) 'Number of atoms in -p and -c files do not agree!'
         call mexit(6, 1)
      end if
      
      !         ----- Read X
      
      nr3 = natom*3
      read (53) (x(i), i=1,nr3)
      
   else
      
      !       ----- read title etc.
      
      read (53, '(20a4)') ititl
      read (53, '(a)' ) line
      if( line(6:6) == ' ' ) then
         read( line, '(i5)' ) natomd
      else
         read( line, '(i6)' ) natomd
      end if
      if (natomd /= natom) then
         write(6,*) 'Number of atoms in -p and -c files do not agree!'
         call mexit(6, 1)
      end if
      
      !         ----- Read X
      
      nr3 = natom*3
      read (53,3) (x(i), i=1,nr3)
      3 format(6f12.7)
   end if  ! (ntx == 0)
   
   write (6,*)
   write(6,4) ititl
   4 format(' Getting coordinates from file with title:'/5x,20a4)
   close(53)
   return
end subroutine getcor 
