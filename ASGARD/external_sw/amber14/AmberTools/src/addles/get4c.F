      subroutine get4c(c4,empty)
c
c get4c return the first four characters of the first expr (of line)
c which is not blank
c

      logical empty
      character*4 c4
      character*1 blank

      integer i,k,l

#include "LINE.h"
#include "UNITS.h"
#include "DEBUG.h"

      data blank/' '/

      c4 = ' '
      empty = .true.
      do 1 i=0,nexp-1
      	if (i.eq.0) then
      		k=1
      		l=point(1)
      	else
      		k=point(i)+1
      		l=point(i+1)
      	end if
      	 if (line(k:k).ne.blank) then
      		if (l-k+1 .ge. 4) then
      		 c4(1:4) = line(k:k+3)
      		else if (l-k+1 .lt. 4) then
      		 c4(1:l-k+1) = line(k:l)
      		 c4(l-k+2:4) = blank
      		end if
      		line(k:l) = ' '
      		empty = .false.
      		return
      	 end if
1	continue
      return
      end
