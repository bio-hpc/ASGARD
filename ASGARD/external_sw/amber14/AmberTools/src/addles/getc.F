      character*(*) function getc(expr,def,lcpext)
c check if command line (line) includes a key word (expr). If yes
c get getc, if not set getc to def
c
      integer lcpext
      character*4 expr
      character*1 BLANK
      character*80 def
#include "LINE.h"
#include "DEBUG.h"
#include "UNITS.h"

c
c local
c
      integer i,j,k,lcp
      DATA BLANK/' '/
c
c verify that the line exists ...
c
      if (point(100).le.0 .or. nexp.le.0 ) then
      	call alert('getc',4,'No input line',13,0)
      	return
      end if
        j=1
      do 2 i=1,nexp
      	if(expr(1:4).eq.line(j:j+3)) then
      	 do 1 k=j,point(i)
c
c check if this is a substitution command
c
      	  if (line(k:k).eq.'=') then
                     getc(1:point(i)-k-2)=line(k+2:point(i)-1)
      	     line(j:point(i))=BLANK
           	     lcp=point(i)-k-2
      	     if (lcp.lt.lcpext) getc(lcp+1:lcpext)=BLANK
      	     if (lcpext.eq.0) lcpext = lcp
      	     return
                  end if
1		 continue
      	 call alert('getc',4,'Missing =',9,0)
      	 getc(1:lcpext)=def(1:lcpext)
      	 return
      	end if
      	j=point(i)+1
2	continue
      getc(1:lcpext)=def(1:lcpext)
      return
      end
