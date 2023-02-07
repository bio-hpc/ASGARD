
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine distr here]
subroutine distr (mn, n, wr, wi, ntrun)
   
   !                         j. kottalam
   
   !     --distributions of langevin modes:
   !         ntrun=4: print distribution of frequencies (i.e. imaginary
   !                     parts of lambda
   !         ntrun=5: print distribution of decay times (i.e. real
   !                     parts of lambda
   !         ntrun=6: same as 4, but only for overdamped modes
   
   implicit double precision (a-h,o-z)
   parameter (mkount=500)
#  include "files.h"
   dimension wr(mn), wi(mn), kount(mkount)
   
   neig = 0
   call  lmin (mn, n, neig, wr, wi, dummy, idummy)
   
   if (ntrun == 4) then
      idel = 10
      del = 10.0/108.59
      !                 ----- t0 inverse to wave number
      call count (wi, n, kount, mkount, k, del)
   else
      if (ntrun == 6) then
         n = 0
         !         do while (wi(n+1).eq.0.0)
         do 10 nnn=1,999999
            if (wi(n+1) /= 0.0) goto 11
            n = n + 1
         10 continue
         11 continue
      end if
      do 20 i = 13, n
         wr(i) = 1.0 / wr(i)
      20 continue
      idel = 1
      del = 20.455
      !                 ----- t0 to ps
      call dsort (n-12, wr(13))
      call count (wr, n, kount, mkount, k, del)
   end if
   
   call amopen (19, 'DENS', owrite, 'F', 'W')
   write (19,*) ' &data'
   write (19,*) ' x(1,1) = '

   call amopen (20, 'CUMU', owrite, 'F', 'W')
   write (20,*) ' &data'
   write (20,*) ' x(1,1) = '
   
   do 30 i = 1, k
      write (19,'(1x,f10.4,a1)') idel*(i-0.5), ','
      write (20,'(1x,f10.4,a1)') idel*(i-0.5), ','
   30 continue
   
   write (19,*) ' y(1,1) = '
   write (20,*) ' y(1,1) = '
   
   sum = 0.0
   do 40 i = 1, k
      term = float(kount(i))
      sum = sum + term
      write (19,'(1x,f10.4,a1)') term, ','
      write (20,'(1x,f10.4,a1)') sum, ','
   40 continue
   
   write (19,*) ' hscale = 0.85, box = .false., yfact = 0.5,'
   write (19,*) ' histgm = .true., bwidth = ', &
         float(idel), ', ybase = 0.0,'
   write (19,*) ' fillin = 30,'
   write (19,*) ' &end'
   write (19,*) 'density of modes'
   if (ntrun == 4) then
      write (19,*) 'wave number'
   else
      write (19,*) 'decay time (ps)'
   end if
   write (19,*) 'number of modes'
   
   write (20,*) ' hscale = 0.85, box = .false., yfact = 0.5,'
   write (20,*) ' &end'
   write (20,*) 'cumulative distribution'
   if (ntrun == 4) then
      write (20,*) 'wave number'
   else
      write (20,*) 'decay time (ps)'
   end if
   write (20,*) 'number of modes'
   
   return
end subroutine distr 
