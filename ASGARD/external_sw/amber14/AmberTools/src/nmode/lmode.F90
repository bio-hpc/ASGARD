
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
!+ [Enter a one-line description of subroutine lmode here]
subroutine lmode(natom, a, wr, wi, z, fv1, iv1, dd, amass, x, &
      expos, eta, gamma, winv, nvect, ioseen, hrmax,igraph)
   
   !      ----- Langevin mode analysis
   !      ----- G. Lamm and A. Szabo, J. Chem. Phys. 85, 7334 (1986)
   
   !                              j kottalam 23 sept 87
   
   implicit double precision (a-h,o-z)
   parameter(lwork=24*450)
   character(len=1) jobvl,jobvr
   dimension work(lwork)
   dimension a(6*natom,6*natom), wr(6*natom), wi(6*natom), &
         z(6*natom,6*natom), fv1(6*natom),iv1(6*natom), &
         dd(3*natom*(3*natom+1)/2), amass(natom), gamma(*), &
         winv(3*natom), x(3*natom),igraph(natom)
   common /blocks/ are, prohi, bited
   common /timer/ timtot, tlast, elapst
#  include "files.h"
   
   call timit (timtot, tlast, elapst)
   nr6 = 6*natom
   nr3 = 3*natom
   mn  = nr6
   neig = min (nvect ,6*natom)
   
   !     ----- set up the friction matrix and place it
   !     ----- in the lower right hand corner of A
   
   call setgam (natom, a, amass, expos, eta, gamma, x, ioseen, &
         hrmax, igraph)
   
   !     ----- set double derivatives to lower left of A
   
   ij = 0
   do 20 j = 1,nr3
      do 10 i = nr3+1,nr3+j
         ij = ij + 1
         a(i,j) = - dd(ij)
         a(nr3+j,i-nr3) = -dd(ij)
      10 continue
   20 continue
   
   !     ----- set upper right to unit matrix
   !     ----- and leave upper left zero
   
   do 30 i = 1, nr3
      a(i,i+nr3) = 1.0
   30 continue
   
   !     ---- now, A is ready to be diagonalized, but first
   !     ---- make sure output file can be opened and written to
   
   call amopen(11, lmod, owrite, 'F', 'W')
   write (11,'(2i12)') natom, neig
   write (11,'(a)') ' Langevin frequencies '
   
   !     ----- now diagonalize the matrix
   
   if (nvect > 0) then
      matz  = 1
   else
      matz = 0
   end if
   call timit (timtot, tlast, elapst)
   write(6,'(a,f8.4)') '| Time to set up A matrix: ',elapst
   
   !     --diagonalize with lapack routine:
   
   jobvl = 'N'
   if (nvect > 0) then
      jobvr = 'V'
   else
      jobvr = 'N'
   end if
   ldum = 1
   call dgeev(jobvl,jobvr,nr6,a,nr6,wr,wi,dum,ldum,z,nr6, &
         work,lwork,info)
   if (info /= 0) then
      write(6,*) ' error code returned from dgeev :', info
   end if
   
   !     --following for historical interest only, using rg:
   
   !     call rg (mn, nr6, a, wr, wi, matz, z, iv1, fv1, ierr)
   !     if (ierr.ne.0) then
   !       write(6,*) ' error code returned from rg :', ierr
   !     end if
   
   call timit (timtot, tlast, elapst)
   write(6,'(a,f8.4)') '| Time to diagonalize A matrix: ',elapst
   
   !     ----- Matrix A has been utilized, its space can now be used.
   !     ----- Sort the eigenvalues according to the imaginary parts.
   
   call spsort (nr6, a(1,2), wr, wi, a(1,3), nreal)
   write (6,'(/,a,i8,a,/)') ' There are ', nreal, ' real modes'
   
   if (neig > 0) then
      
      !       ----- compute normalized L matrix
      
      call norml (neig, mn, nr6, z, wr, wi, gamma, a, a(1,3))
      
      !       ----- remove translational and rotational modes
      
      if (nreal > 12) then
         
         call movecm (x, amass, natom)
         do 40 k = 1, 6*natom
            a(k,4) = 0.0
         40 continue
         if (nreal > neig) nreal = neig
         call remrot (natom, nreal, a(1,5), x, z, wr, wi, a(1,3), &
               a(1,4), amass, a(1,2), a)
         
      end if
      
      !       ----- mass "unweight" L
      
      call mweit (mn, neig, a(1,3), z, amass, winv)

   end if
   call timit (timtot, tlast, elapst)
   write(6,'(a,f12.4)')'| Time to normalize, sort, mass weight: ', &
         elapst
   
   !     ----- output langevin modes
   
   call langout (neig, a(1,3), mn, nr6, wr, wi, z, a, a(1,4))
   
   return
end subroutine lmode 
