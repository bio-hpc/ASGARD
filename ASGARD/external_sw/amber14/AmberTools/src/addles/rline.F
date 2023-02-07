      subroutine rline(name,namel,unit)
c
c       input:
c       name - the name of the subroutine which called "rline"
c		(character)
c	namel - the length of the character name
c	unit  - unit number of stream file
c
c       output:
c	line - character of length 300 which includes the command
c		line
c	nexp - the number of experssions found in line,  maximum 99
c	point(i) - pointer to the position in line in which 
c		exppression i starts

c	maximum 99 expressions are allowed in a line and it can
c	have at most 300 characters.
c	point(100) is the end of the string.
c	line is compressed on output to include no separators (spaces).
c
c
c	this subroutine reads a line and separate it to expressions
c
#include "LINE.h"
#include "UNITS.h"
#include "DEBUG.h"
      character*80 name
      character*80 lseg
      integer unit,namel,empty
      logical space
c
c local
c
      integer i,istart,iend
      integer jend,j,k,level
c
c initialize the pointer
c
      point(100)=0
      empty = 0
c
c	find the end of the line
c
      jend=0
      iend=0
1	continue
      istart=jend+1
      jend=istart+80-1
      if (iend.gt.300) then
      	write(*,*)' line exceeds 300 characters'
      	write(*,*)' too long... command ignored'
      end if

11	read(unit,2,end=99,err=999)lseg
2	format(a)
      do 22 i=istart,jend
       if (lseg(i:i).eq.'~') then
c
c this is a remark line. Find end of line
c
        do 21 k=80,1,-1
         if (lseg(k:k).ne.' ') then
          call echo(name,namel,lseg,k)
          go to 11
         end if
21	  continue
       end if
22	continue
      line(istart:jend)=lseg(1:80)

      iend=jend
      do 3 i=iend,istart,-1
         if (line(i:i).ne.' ') then
      	point(100)=i
c
c		now test if this line has continuation
c
      	if (line(i:i).eq.'-') then
c
c the "-" for the continuation is ignored
c
      	 jend=i-1
      	 go to 1
      	end if
      	go to 4
       end if
3	continue
4	continue
c
c	If the line is empty simply print empty line and read another
c
      if (point(100).eq.0) then
       write(*,100)
100	 format(/)
       go to 11
      end if

      call echo(name,namel,line,point(100))
c
c       and now find the different words
c
      nexp=0
      i=0
      space = .true.
6	continue
      i=i+1
c
c separator (space) ?
c
        if(line(i:i).eq.' ') then
c
c if the space is not at the beginning of the line
c
      	if (i.ne.1) then
c
c test that the number of expressions-nexp is less than 99
c
      	if (.not.space) nexp=nexp+1
      	if (nexp.gt.99) then
      		write(*,*)' Line with more than 99 words !'
      		write(*,*)' Cannot read line, ignored'
      		point(100)=0
      		return
      	end if
c
c set the pointer to expression number nexp.
c
      		point(nexp)=i-1
      	end if
c
c compress the string by eliminating the space
c
      	do 5 j=i+1,point(100)
      		line(j-1:j-1)=line(j:j)
5		continue
c we have to check the point i now since we may have move space
c to i so...
c
      	i=i-1
c decrease the string length after space compression
c
      	point(100)=point(100)-1
         	space = .true.
        else
      	space = .false.
        end if
c loop check if end of stricg was reached
c
      if (i.lt.point(100)) go to 6
c
c since the end of the line does not include space there is
c one more expression at the end
c
      nexp=nexp+1
c
c the end of the last expression is in point(100)
c
      point(nexp)=point(100)
      return
99	continue
      empty = empty + 1
      if (empty.lt.50) go to 11
      write(*,*)
      write(*,*)
      write(*,*) 'End of File '
      write(*,*) 'last called by routine: ',name(1:namel)
      level = 1
      call alert(' File too short for read required ',34)
999	continue
      write(*,*)
      write(*,*)
      write(*,*) 'Error during read'
      write(*,*) 'last called by routine: ',name(1:namel)
      level = 1
      call alert(' Input problems ',17)
      return
      end
