
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
!+ [Enter a one-line description of subroutine eigout here]
subroutine eigout (ns3, neig, val, vec, x, mass, ivform)
   implicit double precision (a-h,o-z)
   double precision mass(*)
   dimension val(ns3), vec(ns3,ns3), x(ns3)
   parameter (maxat=2800)
   dimension xw(3*maxat)
   
   do i=1,ns3
      if (val(i) > 0.0) then
         val(i) = 108.587 * sqrt(val(i))
      end if
   end do
   !     write (6,3) (val(i),i=1,ns3)
   
   if (neig == 0) return
   if (ivform == 0) then
      
      !   --- unformatted, single precision eigenvector output file:
      
      if (ns3 > 3*maxat) then
         write(6,*) 'MAXAT is not big enough in eigout: ',ns3,maxat
         call mexit(6,1)
      end if
      write(61) ns3
      do i=1,ns3
         xw(i) = x(i)
      end do
      write(61) (xw(i),i=1,ns3)
      do j=1,neig
         val4 = val(j)
         write(61) j,val4
         do i=1,ns3
            xw(i) = vec(i,j)
         end do
         write(61) (xw(i),i=1,ns3)
      end do
   else if (ivform == 1 ) then
      
      !   --- (amber) formatted eigenvector output file:
      
      if(idecomp == 0) then
         write(61,1)
         write(61,2) ns3
         write(61,3) x
         do j = 1,neig
            write (61,'(1x,a)') '****'
            write (61,2) j, val(j)
            write (61,3) (vec(i,j),i=1,ns3)
         end do
      end if
      
   else if (ivform == 2) then
      
      !   --- (mkl) formatted eigenvector output file, for molekel:
      
      write(61,4)
      write(61,5)
      natom = ns3/3
      do i=1,natom
         ino = 1
         if( abs(mass(i)-12.d0) < 0.25 ) ino=6
         if( abs(mass(i)-13.01d0) < 0.05 ) ino=6  ! united atom CH
         if( abs(mass(i)-15.03d0) < 0.05 ) ino=6  ! united atom CH3
         if( abs(mass(i)-14.d0) < 0.25 ) ino=7
         if( abs(mass(i)-16.d0) < 0.25 ) ino=8
         if( abs(mass(i)-19.d0) < 0.25 ) ino=9
         if( abs(mass(i)-31.d0) < 0.25 ) ino=15
         if( abs(mass(i)-32.d0) < 0.25 ) ino=16
         if( abs(mass(i)-35.45d0) < 0.25 ) ino=17
         write(61,7) ino, x(3*i-2), x(3*i-1), x(3*i)
      end do
      write(61,6)
      
      write(61,12)
      write(61,13)
      write(61,6)
      
      write(61,8)
      do j=1,neig,3
         write(61,9)
         write(61,10) val(j),val(j+1),val(j+2)
         do i=1,natom
            f= sqrt(mass(i))
            write(61,11)f*vec(3*i-2,j),f*vec(3*i-1,j),f*vec(3*i,j), &
                  f*vec(3*i-2,j+1),f*vec(3*i-1,j+1),f*vec(3*i,j+1), &
                  f*vec(3*i-2,j+2),f*vec(3*i-1,j+2),f*vec(3*i,j+2)
         end do
      end do
      write(61,6)
      
   end if  ! (ivform == 0)
   close (61)
   
   1 format ('NORMAL COORDINATE FILE')
   2 format (i5,f12.5)
   3 format (7f11.5)
   4 format ('$MKL')
   5 format ('$COORD')
   6 format ('$END')
   7 format (i5,3f15.5)
   8 format ('$FREQ')
   9 format ('  A1                      A1                      A1')
   10 format (f8.2,16x,f8.2,16x,f8.2)
   11 format (9f8.4)
   12 format ('$CHAR_MULT')
   13 format ('0 1')
   
   return
end subroutine eigout 
