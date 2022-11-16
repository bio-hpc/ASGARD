      integer function geti(expr,def)
c check if command line (line) includes a key word (expr). If yes
c get geti, if not set geti to def
c
      character*4 expr
      integer def
#include "LINE.h"
#include "UNITS.h"
#include "DEBUG.h"
c
c local
c
      integer i,j,k
c
c verify that the line exists ...
c
      if (point(100).le.0 .or. nexp.le.0 ) return
        j=1
      do 2 i=1,nexp
      	if(expr(1:4).eq.line(j:j+3)) then
      	 do 1 k=j,point(i)
c
c check if this is a substitution command
c
      	  if (line(k:k).eq.'=') then
      	   rewind jnkf
      	   write(jnkf,*)line(k+1:point(i))
      	   rewind jnkf
      	   read(jnkf,*,err=11)geti
      	   line(j:point(i)) = ' '
      	   return
                  end if
1		 continue
      	 write(*,*)' **missing = in expression. Nothing done'
      	 write(*,*)' **variable was set to default'
      	 geti=def
      	 return
      	end if
      	j=point(i)+1
2	continue
      geti=def
      return
11	continue
      write(*,*)' *** error while reading variable in line '
      write(*,*)line
      write(*,*)' *** geti ignored'
      return
      end
