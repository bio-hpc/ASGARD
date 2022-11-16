      subroutine pick(ipick,ngrp,origpt)
c
c pick particles. On return ipick has integer values for all
c the atoms which were selected. If the group keyword
c was used ipick will have integer values from 1 up according
c to the group found. If only one group exists ipick is set to
c 1 for picked atoms. If no selection is made all particles are
c picked. It is assumed that rline was called BEFORE pick
c
#include "SIZE.h"
#include "TOP.h"
#include "UNITS.h"
#include "DEBUG.h"
      integer*2 ipick1(maxnatom)
      integer ngrp,ipick(*),origpt(*)
      logical find,nxt,first

      integer level,nptpick

      integer i,j,k
c
      nptpick=natom
c
c initialize ipick to default value
      do 1 i=1,nptpick
       ipick(i)  = 1
       ipick1(i) = 0
1	continue
c 
c chang for CARLOS LES
c
c	ngrp = 0
c
      first = .true.
      if (.not.find('pick')) go to 6

c start loop on pick-up line
c Go from left to right and pick up logical expressions
2	continue
      if (nxt('done',4)) go to 6

c .and.(&) or first
      if (nxt('&',1) .or.first) then
       first = .false.
       call pick2(ipick1,ngrp,origpt)
       do 3 i=1,nptpick
        if(ipick(i).gt.0 .and. ipick1(i).gt.0) then
         ipick(i)=ipick1(i)
        else
         ipick(i)=0
        end if
3	 continue
       go to 2

c .or.(|)
      else if (nxt('|',1)) then
       call pick2(ipick1,ngrp,origpt)
       do 4 i=1,nptpick
        if (ipick1(i).gt.0)  ipick(i) = ipick1(i)
4	 continue
       go to 2

c .not. (!=)
      else if (nxt('!=',2)) then
       call pick2(ipick1,ngrp,origpt)
       do 5 i=1,nptpick
        if (ipick1(i).gt.0) ipick(i) = 0
5	 continue
       go to 2
      end if
      level = 1
      write(*,100)
100	format(1x,'Illegal logical pick syntax')
6	continue
      j = 0
      do 7 i=1,nptpick
       if (ipick(i).ne.0) j = j + 1
7	continue
      return
      end
c----------------------------------------------------------


c-----------------------------------------------------------
      subroutine pick2(ipick1,ngrp,origpt)
c
c pick atoms in a logical segment. I.e. expression enclosed
c between & (.and.) or | (.or.) signs.
c
c we need origpt(i) ONLY when looking for atoms by number, and we want the
c original number so user can pick by number even after some new atoms 
c have been defined
c
#include "SIZE.h"
#include "TOP.h"
#include "UNITS.h"
#include "DEBUG.h"
      integer ngrp
      integer intgr,origpt(*)
      integer*2 ipick1(maxnatom)
      logical nxt

      integer level,i,j,k,i1,i2,i5

      integer i3
      logical empty
      character*4 char1

c initialize ipick1
      do i=1,natom
        ipick1(i) = 0
      end do

c change for CARLOS LES
c
c check if groups exist
c	if (nxt('grou',4)) then
c	  ngrp = ngrp + 1
c	  j = intgr()
c	else
c	  ngrp = 1
c	  j    = 1
c	end if
c
c j is number of copies
c
      j=ngrp

c start processing
      if (nxt('#prt',4)) then
        i1 = intgr()
        i2 = intgr()
        do 3 i=i1,i2
c
c change for CARLOS LES
c
      	do 300 i5=1,natom
      		if (origpt(i5).eq.i) then
          			ipick1(i5) = j
      		endif
300		continue
3	  continue
c
c carlos's add: side chain pick by mono # (all but c,o,n,h,ca,ha)
c
      else if (nxt('#sid',4)) then
        i1 = intgr()
        i2 = intgr()
c
c CHANGE FOR CARLOS LES
c
        if (i2.ge.nres) then
         i3=natom
        else
         i3 = ipres(i2+1)-1
        end if
        do 4 i=ipres(i1),i3
          if (igraph(i).ne.'C   '.and.igraph(i).ne.'O   '.and.
     &		igraph(i).ne.'N   '.and.igraph(i).ne.'HN  '.and.
     &          igraph(i).ne.'H   '.and.
     &		igraph(i).ne.'CA  '.and.igraph(i).ne.'HA  ') then
          		ipick1(i) = j
          endif
4	  continue
c
c pick by mono #
c
      else if (nxt('#mon',4)) then
        i1 = intgr()
        i2 = intgr()
c
c CHANGE FOR CARLOS LES
c
        if (i2.ge.nres) then
         i3=natom
        else
         i3 = ipres(i2+1)-1
        end if
        do 144 i=ipres(i1),i3
          ipick1(i) = j
144	  continue
c
c carlos add:
c pick mono by mono #, BUT PUT MONO CUT POINT AT CA-C (include CO of prev
c mono, exclude CO of this mono)
c
      else if (nxt('#cca',4)) then
        i1 = intgr()
        i2 = intgr()
        if (i2.ge.nres) then
         i3=natom
        else
         i3 = ipres(i2+1)-1
        end if
c
c add prev mono CO if there is a prev mono
c
        if (ipres(i1).ne.1) then
      	do 104 i=ipres(i1-1),ipres(i1)-1
c
c compare actual name of this atom
c
      		if (igraph(i).eq.'C   '.or.igraph(i)
     &			.eq.'O   ') then
          			ipick1(i) = j
      		endif
104		continue
        endif
c
c add the mono
c
        do 154 i=ipres(i1),i3
          		ipick1(i) = j
154	  continue
c
c IMPORTANT!!!!!
c remove the final mono's CO (unless it is the last mono)
c but some judgement here is required by the user!
c
        if (i3.ne.natom) then
            do 164 i=ipres(i2),i3
      	if (igraph(i).eq.'C   '.or.igraph(i)
     &		.eq.'O   ') then
          		ipick1(i) = 0
      	endif
164	    continue
        endif
c 
c chemical name picks
c
      else if (nxt('chem',4)) then
        call get4c(char1,empty)
        if (empty) then
         level = 1
         call alert('Illegal chem pickup',19)
        end if
        if (char1.eq.'prtc') then
         call get4c(char1,empty)
         if (empty) then
          level = 1
          call alert('missing particle name',21)
         end if
         do 5 i=1,natom
c
c check name match as well as possible wildcards
c
          if (igraph(i).eq.char1 .or. char1(1:1).eq.'*') ipick1(i) = j
          if (igraph(i)(1:1).eq.char1(1:1).and.char1(2:2).eq.'*') 
     &		ipick1(i) = j
          if (igraph(i)(1:2).eq.char1(1:2).and.char1(3:3).eq.'*') 
     &		ipick1(i) = j
          if (igraph(i)(1:3).eq.char1(1:3).and.char1(4:4).eq.'*') 
     &		ipick1(i) = j
5	   continue
        else if (char1.eq.'mono') then
         call get4c(char1,empty)
         if (empty) then
          level = 1
          call alert('Illegal chem pickup',19)
         end if
         if (char1(1:4).eq.'none') then
          level = 1
          call alert('Illegal chem pickup',19)
         end if
         do 7 i=1,nres
          if ((labres(i).eq.char1 .or. char1(1:1).eq.'*').or.
     &	    (labres(i)(1:1).eq.char1(1:1).and.char1(2:2).eq.'*').or.
     &	    (labres(i)(1:2).eq.char1(1:2).and.char1(3:3).eq.'*').or.
     &	    (labres(i)(1:3).eq.char1(1:3).and.char1(4:4).eq.'*')) then
c
c CHANGE FOR CARLOS LES
c
           if (i.ge.nres) then
            i2 = natom
           else
            i2 = ipres(i+1)-1
           end if
           do 6 k=ipres(i),i2
            ipick1(k) = j
6	     continue
          end if
7	   continue
        end if
       end if
       return
       end

