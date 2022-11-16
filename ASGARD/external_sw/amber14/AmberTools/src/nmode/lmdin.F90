
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine lmin here]
subroutine lmin (nat3, nat6, nvect, wr, wi, z, iclass)
   
   !                            j. kottalam
   
   implicit double precision (a-h,o-z)
#  include "files.h"
   dimension wr(*), wi(*), z(nat3,nvect), iclass(*)
   
   call amopen (10, lmod, 'O', 'F', 'R')
   
   read (10,'(3i12)') natom, nvecs, neigv
   
   if (neigv == 0) neigv = nat6
   if (3*natom /= nat3) then
      write(6,*) 'Number of atoms in lmode and inpcrd do not match.',  &
           natom,nat3
      call mexit(6, 1)
   end if
   if (nvecs < nvect) then
      write(6,*) 'Only ', nvecs,' vectors found in lmode file'
      nvect = nvecs
   end if
   
   read (10,*)
   
   !     ----- read eigenvalues
   
   do 10 i = 1, neigv
      read (10,'(1x, i4, 3x, f10.0, 2x, f10.0)') &
            idummy, wr(i), wi(i)
      
      !     ---- convert from cm**-1 to internal units:
      
      wr(i) = wr(i) / 108.59
      wi(i) = wi(i) / 108.59
   10 continue
   
   !     ----- read eigenvectors
   
   do 30 j = 1, nvect
      read (10,*)
      read (10,*)
      do 20 i = 1, nat3, 6
         read (10,'(6(2x,e11.0))') (z(i+l,j),l=0,5)
      20 continue
      read (10,'(2x,i10)') iclass(j)
   30 continue
   
   close (10)
   return
end subroutine lmin 
