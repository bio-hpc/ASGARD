      subroutine getunit(iounit)
      integer   iounit, low, high, next, range
      integer namel,i,level
#include "UNITS.h"
      parameter (low = 26, high = 99, range=high-low+1)
      logical used(low:high)
      character*7 name
      data next /low/
      data used/range*.false./
      save next, used

      name  = 'getunit'
      namel = 7
c Loop over possibilities, and find a unit that
c is not attached to an open file (just in case someone
c used a unit without checking with this subroutine).
      do 1 i=low,high
         if (next .gt. high) next = low
         if (.not.used(next)) goto 100
         next = next + 1
1     continue

20    continue
c !!!!!!!!!put error message here, since no unit can be assigned!
      level = 1
      call alert('Unit cannot be open',19)
      iounit = -1
      return

c success---a unit is found, recorded, and returned.
100   continue
      iounit = next
      next = next+1
      used(iounit) = .true.
      write(*,101)iounit
101   format(1x,' The following unit number was assigned ',i5)
      return


C frees a unit number if no longer needed
      entry freeunit(iounit)
      used(iounit) = .false.
      return


      entry closeall
c  close all files that might have been opened
      do 500 i=low, high
         if (used(i)) then
            close(i)
            used(i) = .false.
         end if
500   continue
      return

      end

