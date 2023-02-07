!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     contains all the routines to convert strings to numbers
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine se_rdnum(string,istart,value,ierror)
!
!     extracts a double precision floating point number from a character
!     string.  the field of search starts at string(istart:istart) or
!     after the first equals sign following the istart position.  the
!     number is returned in value.  if an error is encountered, ierror
!     is set to one.  this routine expects that there are no blank spaces
!     embedded anywhere within the numerical field.
!
      implicit none
      integer istart, ierror
      double precision value
      character string*(*),efield(4)*1
      data efield /'E','e','D','d'/
      save efield
      integer ibeg, istop, iend, i, ieq

      ierror = 0
      ibeg = istart
      istop = len(string)
      iend = istop
      do i=istart,istop
        if(string(i:i).eq.' ')then
          iend = i-1
          exit
        endif
      enddo
      if(iend.lt.ibeg)then
        ierror = 1
      else
        ieq = index(string(ibeg:iend),'=')
        if(ieq.ne.0) ibeg = ibeg + ieq
        call se_getnum(string,ibeg,iend,value,ierror)
      endif
      end
!--------------------------------------------------------------------------------
      subroutine se_getnum(string,ibeg,iend,value,ierror)
      implicit none

      integer ibeg, iend, ierror
      double precision value

      character string*(*),char*1,efield(4)*1
      data efield /'E','e','D','d'/
      save efield

      double precision asign, xleft, xright, esign, expart
      integer idecml, i1, i2, i, ie


      value = 0.0d0
!
!     check for algebraic sign.
!
      char = string(ibeg:ibeg)
      if(char.eq.'-')then
        asign = -1.0d0
        ibeg = ibeg + 1
      elseif(char.eq.'+')then
        asign = 1.0d0
        ibeg = ibeg + 1
      else
        asign = 1.0d0
      endif
      if(ibeg.gt.iend)then
        ierror = 1
        return
      endif
!
!     first determine the whole number equivalent of whatever is
!     to the left of any decimal point.
!
      idecml = index(string(ibeg:iend),'.')
      if(idecml.eq.1)then
        if(ibeg.eq.iend)then
!
!         number is just a decimal point.  assume a value of zero.
!
          value = 0.0d0
          return
        endif
        xleft = 0.0d0
        ibeg = ibeg+1
      else
        i1 = ibeg
        if(idecml.eq.0)then
          i2 = iend
        else
          i2 = ibeg+idecml-2
        endif
        call se_whole(string,i1,i2,xleft,ierror)
        if(ierror.ne.0) return
        value = xleft*asign
        if(idecml.eq.0.or.i2.eq.(iend-1)) return
        ibeg = i2+2
      endif
!
!     determine the whole number equivalent of whatever is to the
!     right of the decimal point.  account for e or d field format.
!
      do i=1,4
        ie = index(string(ibeg:iend),efield(i))
        if(ie.ne.0) exit
      enddo
      if(ie.eq.1)then
        value = xleft*asign
        ibeg = ibeg + 1
      else
        i1 = ibeg
        if(ie.eq.0)then
          i2 = iend
        else
          i2 = ibeg+ie-2
        endif
        call se_whole(string,i1,i2,xright,ierror)
        if(ierror.ne.0) return
        xright = xright*10.0d0**(i1-i2-1)
        value = value + xright*asign
        if(ie.eq.0.or.i2.eq.(iend-1)) return
        ibeg = i2+2
      endif
!
!     get the exponential portion.
!
      char = string(ibeg:ibeg)
      if(char.eq.'-')then
        esign = -1.0d0
        ibeg = ibeg + 1
      elseif(char.eq.'+')then
        esign = 1.0d0
        ibeg = ibeg + 1
      else
        esign = 1.0d0
      endif
      if(ibeg.gt.iend) return
      i1 = ibeg
      i2 = iend
      call se_whole(string,i1,i2,expart,ierror)
      if(ierror.ne.0) return
      value = value*10.0d0**(esign*expart)
      return
      end
!--------------------------------------------------------------------------------
      subroutine se_whole(string,ibeg,iend,value,ierror)
!
!     returns the whole number in the field string(ibeg:iend).  only
!     the numbers 0-9 are allowed to be present.
!
      implicit none
      integer ibeg, iend, ierror
      double precision value
      character string*(*)
      integer i, idigit, ichar0

      ierror = 0
      value = 0.0d0
      ichar0 = ichar('0')
      do i=ibeg,iend
        idigit = ichar(string(i:i)) - ichar0
        if(idigit.lt.0.or.idigit.gt.9)then
          ierror = 1
          return
        endif
        value = 10.0d0*value + idigit
      enddo
      return
      end
!--------------------------------------------------------------------------------
      subroutine se_rdinum(string,istart,ivalue,ierror)

      implicit none
      integer istart, ivalue, ierror
      character string*(*)

      integer ibeg, istop, iend, i, ieq

      ierror = 0
      ibeg = istart
      istop = len(string)
      iend = istop
      do i=istart,istop
        if(string(i:i).eq.' ')then
          iend = i-1
          exit
        endif
      enddo
      if(iend.lt.ibeg)then
        ierror = 1
        return
      endif
      ieq = index(string(ibeg:iend),'=')
      if(ieq.ne.0) ibeg = ibeg + ieq
      call se_getinum(string,ibeg,iend,ivalue,ierror)
      return
      end
!--------------------------------------------------------------------------------
      subroutine se_getinum(string,ib,ie,ivalue,ierror)

      implicit none
      integer ib, ie, ivalue, ierror
      character string*(*)
      
      integer isign

      ierror = 0
      ivalue = 0
      if (string(ib:ib) .eq. '-') then
         ib = ib + 1
         isign = -1
      elseif (string(ib:ib) .eq. '+') then
         ib = ib + 1
         isign = 1
      else
         isign = 1
      endif

      call se_iwhole(string,ib,ie,ivalue,ierror)
      ivalue = ivalue * isign

      return
      end
!--------------------------------------------------------------------------------
      subroutine se_iwhole(string,ibeg,iend,ivalue,ierror)
!
!     returns the whole number in the field string(ibeg:iend).  only
!     the numbers 0-9 are allowed to be present.
!
      implicit none
      integer ibeg, iend, ivalue, ierror
      character string*(*)

      integer i, idigit, ichar0

      ivalue = 0
      ichar0 = ichar('0')
      do i=ibeg,iend
        idigit = ichar(string(i:i)) - ichar0
        if(idigit.lt.0.or.idigit.gt.9)then
          ierror = 1
          return
        endif
        ivalue = 10*ivalue + idigit
      enddo
      return
      end
!--------------------------------------------------------------------------------
      subroutine se_iatoi(str, istart, lstr, integ, ierror)

      implicit none
      integer istart, lstr, integ, ierror

      character str*(*), ch
      logical int, min

      integer i, izero, nstr

      integ = 0
      izero = ichar('0')
      nstr = len(str)
      do i=istart,nstr
         ch = str(i:i)
         call se_whatis2(ch, int, min)
         if (.not.int) exit
      enddo
      lstr = i-1
      if (lstr.eq.0) return
      call se_getinum(str,istart,lstr,integ,ierror)

      end
!--------------------------------------------------------------------------------
      subroutine se_iatoimp(str, istart, lstr, integ, ierror)

      implicit none
      integer istart, lstr, integ, ierror
      character str*(*), ch
      logical int
      integer i, izero, nstr

      integ = 0
      izero = ichar('0')
      nstr = len(str)
      do i=istart,nstr
         ch = str(i:i)
         call se_whatis1i(ch, int)
         if (.not.int) exit
      enddo
      lstr = i-1
      if (lstr.eq.0) return
      call se_getinum(str,istart,lstr,integ,ierror)

      end
