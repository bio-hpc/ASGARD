      logical function nxt(expr,length)
c
c check if the next expression (that is not blank) includes
c a key word (expr). If yes
c set "nxt" to .true. otherwise nxt=.false.
c
      character*4 expr
      integer length
#include "LINE.h"
#include "UNITS.h"
#include "DEBUG.h"
c
c local
c
      integer i,j,k
      nxt = .false.
c
c verify that the line exists ...
c
      if (point(100).le.0 .or. nexp.le.0 ) return
        j=1
      i=1
2	continue
      	k = point(i) - j + 1
      	if (k.gt.4) k = 4
      	if(expr(1:length).eq.line(j:j+length-1)) then
      	 nxt = .true.
      	 line(j:point(i))=' '
      	 return
      	end if
1		continue
      	if (line(j:j).ne.' ') return
      	if (i.eq.nexp) return
      	j = point(i)+1
      	i = i+1
      go to 2
      end
