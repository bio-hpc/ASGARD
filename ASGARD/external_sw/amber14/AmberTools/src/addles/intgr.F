      function intgr()
c
c integer return the first integer number in a line (separated
c by spaces from previous expressions) and erase the number from the
c line
c
      integer intgr
      character*11 legal
c
c a floating point number can start with one of the characters below
c

      integer i,j,k,l,level
      logical success

#include "UNITS.h"
#include "LINE.h"

      data legal/'1234567890-'/

      k = 1
      l = point(1)
      success = .false.
      do 2 i=0,nexp
      	if (i.eq.0) then
      		k=1
      		l=point(1)
      	else
      		k=point(i)+1
      		l=point(i+1)
      	end if
         	do 1 j=1,11
      	 if (line(k:k).eq.legal(j:j)) then
      		rewind jnkf
      		write(jnkf,*)line(k:l)
      		rewind jnkf
      		read(jnkf,*,err=999)intgr
      		success = .true.
      		line(k:l) = ' '
      		go to 3
      	 end if
1		continue
2	continue
      if (.not. success) then
       level = 1
       call alert('no integer',10)
       return
      end if
3	return
999	continue
      level = 1
      call alert('Format error in integer number',30)
      return
      end
