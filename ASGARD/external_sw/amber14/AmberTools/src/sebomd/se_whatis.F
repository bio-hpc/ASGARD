      subroutine se_whatis1(this, float)

      implicit none

      character this
      logical float

      integer ithis, i0, i9

      float = .false.

      ithis = ichar(this)
      i0 = ichar('0')
      i9 = ichar('9')
      if ((ithis.ge.i0).and.(ithis.le.i9)) then
         float = .true.
      elseif (this.eq.'.') then
         float = .true.
      elseif (this.eq.'-') then
         float = .true.
      elseif (this.eq.'+') then
         float = .true.
      elseif (this.eq.'E') then
         float = .true.
      elseif (this.eq.'e') then
         float = .true.
      elseif (this.eq.'D') then
         float = .true.
      elseif (this.eq.'d') then
         float = .true.
      endif

      end
!--------------------------------------------------------------------------------
      subroutine se_whatis1i(this, int)

      implicit none

      character this
      logical int

      integer ithis, i0, i9

      int = .false.

      if (this.eq.'-') then
         int = .true.
      elseif (this.eq.'+') then
         int = .true.
      else
         ithis = ichar(this)
         i0 = ichar('0')
         i9 = ichar('9')
         if ((ithis.ge.i0).and.(ithis.le.i9)) int = .true.
      endif

      end
!--------------------------------------------------------------------------------
      subroutine se_whatis2(this, int, min)

      implicit none

      character this
      logical int, min

      integer ithis, i0, i9

      int = .false.
      min = .false.

      if (this.eq.'-') then
         min = .true.
      else
         ithis = ichar(this)
         i0 = ichar('0')
         i9 = ichar('9')
         if ((ithis.ge.i0).and.(ithis.le.i9)) int = .true.
      endif
      end
!--------------------------------------------------------------------------------
      subroutine se_whatis7(this,char,num,parl,parr,comma,eq,white)
      implicit none

      character this
      logical char,num,parl,parr,comma,eq,white

      integer ithis, i0, i9, ia, iz, iaa, izz

      char = .false.
      num = .false.
      parl = .false.
      parr = .false.
      comma = .false.
      eq = .false.
      white = .false.
      if (this.eq.' ') then
         white = .true.
      elseif (this.eq.',') then
         comma = .true.
      elseif (this.eq.'=') then
         eq = .true.
      elseif (this.eq.'(') then
         parl = .true.
      elseif (this.eq.')') then
         parr = .true.
      elseif (this.eq.'/') then
         eq = .true.
      elseif (this.eq.'+') then
         num = .true.
      elseif (this.eq.'.') then
         num = .true.
      elseif (this.eq.'-') then
         num = .true.
      else
         ithis = ichar(this)
         ia = ichar('a')
         iz = ichar('z')
         iaa = ichar('A')
         izz = ichar('Z')
         i0 = ichar('0')
         i9 = ichar('9')
         if (((ithis.ge.ia).and.(ithis.le.iz)).or.
     &        ((ithis.ge.iaa).and.(ithis.le.izz))) then
            char = .true.
         elseif ((ithis.ge.i0).and.(ithis.le.i9)) then
            num = .true.
         else
            white = .true.
         endif
      endif
      end
