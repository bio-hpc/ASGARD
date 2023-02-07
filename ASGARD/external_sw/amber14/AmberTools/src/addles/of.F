      function of()
c
c Open File: a routine to open a file, does the following
c (i) extract a file name and check that file exists
c (ii) extract a unit number and return the unit number as "of"
c
      implicit none
#include "SIZE.h"
#include "MISC.h"
#include "LINE.h"
#include "DEBUG.h"
      character*80 namef, def
      character*80 getc
      integer lc, level, packdef
      logical exist
c get file name
c
      namef = ' '
      lc = 80
      def = 'none'
      namef = getc('name',def,lc)
      packdef = npack
      npack = geti('pack',packdef)
      if (debug) then
       write(*,*)' name of file ',namef(1:lc)
       write(*,*)' length of name ',lc
      end if
      inquire(file=namef(1:lc),exist=exist)
      if (find('read')) then
       of = geti('unit',-1)
       if (of.eq.5) then
        write(*,100)
        return
       end if
       call getunit(of)
       if (of.eq.-1) then
        of = stdi
        return
       else if (of.eq.stdi) then
        write(*,100)
100     format(/,1x,'Input taken from standard file',/)
        return
       else if (find('bina')) then
        if (.not.exist) then
         level = 1
         call alert('READ file missing',17)
         return
        end if
        open(unit=of,file=namef(1:lc),form='unformatted',
     1       status='old')
        rewind of
        return
       else
        if (.not.exist) then
         level = 1
         call alert('READ file missing',17)
         return
        end if
        open (unit=of,file=namef(1:lc),status='old')
        rewind of
        return
       end if
      else if (find('wovr')) then
       call getunit(of)
       if (find('bina')) then
        open(unit=of,file=namef(1:lc),form='unformatted',
     1       status='unknown')
        rewind of
        return
       else
        open (unit=of,file=namef(1:lc),status='unknown')
        rewind of
        return
       end if
      else if (find('writ')) then
       call getunit(of)
       if (exist) then
        write(*,101)
101     format(1x,/,' *** ILLEGAL ATTEMPT TO OPEN FILE ',/,
     1         1x,'*** IF EXISTS MUST SPECIFY read OR wovr ',/)
        level = 1
        call alert(' Illegal Open File statement',28)
       end if
       if (find('bina')) then
        open(unit=of,file=namef(1:lc),form='unformatted',status='new')
        rewind of
       else
        open (unit=of,file=namef(1:lc),status='new')
        rewind of
       end if
      end if
       
      return
      end
