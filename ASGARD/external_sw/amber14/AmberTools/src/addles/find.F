      logical function find(expr)
c check if command line (line) includes a key word (expr). If yes
c set "find" to .true. otherwise find=.false.
c
      character*4 expr
#include "LINE.h"
#include "UNITS.h"
#include "DEBUG.h"
c
c local
c
      integer i,j,k
      find = .false.
c
c verify that the line exists ...
c
      if (point(100).le.0 .or. nexp.le.0 ) return
        j=1
      do 2 i=1,nexp
      	k = point(i) - j + 1
      	if (k.gt.4) k = 4
      	if (k.lt.4 .and. expr(k+1:4).ne.' ') go to 1
      	if(expr(1:k).eq.line(j:j+k-1)) then
      	 find = .true.
      	 line(j:point(i))=' '
      	 return
      	end if
1		continue
      	j=point(i)+1
2	continue
      find = .false.
      return
      end
