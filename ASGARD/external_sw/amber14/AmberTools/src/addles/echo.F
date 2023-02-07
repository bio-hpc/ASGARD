      subroutine echo(sbr,sbrl,line,linel)
c
c echoing the command line
c
#include "UNITS.h"
#include "DEBUG.h"
      character*80 prline
        character*(*) sbr
      character*(*) line
      integer sbrl,linel,istart,lpr1,lpr2,length

c get subroutine name
c
      prline(1:sbrl)=sbr(1:sbrl)
      prline(sbrl+1:sbrl+2)='> '
c insert the command line to the character if too long split
c to several printout
c
      istart=1
1	continue
c length is how many characters still left
c
      length=linel-istart+sbrl+3
c if length less than 80 print all of it and return
c
      if (length.lt.80) then
      	prline(sbrl+3:length)=line(istart:linel)
      	write(*,'(A)')prline(1:length)
      	return
      end if
c if length larger than 80 need to split the line. 
c
      lpr1=istart+80-sbrl-3
      lpr2=80
2	continue
c since we dont want to break an expression, test that
c we break the line at a space. otherwise print less..
c
      if (line(lpr1:lpr1).ne.' ') then
      	lpr1=lpr1-1
      	lpr2=lpr2-1
      	go to 2
      end if
      prline(sbrl+3:lpr2)=line(istart:lpr1)
      write(*,'(A)')prline(1:lpr2)
      istart=lpr1+1
      go to 1
      end
