
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine quasip here]
subroutine quasip(vect,x,proj,ncart)
   
   !     ---- project dyanmics onto modes
   
   implicit double precision(a-h,o-z)
#  include "sizes2.h"
#  include "anal.h"
#  include "files.h"
   real crdr(3*maxatom)
   integer ftyp
   character(len=40) fmt2
   character(len=3) chno
   dimension vect(ncart,*),x(*),proj(*),atmas(maxatom), &
         atmas3(3*maxatom)
   
   !  --- set up formats and  start the output file
   
   nch = iend - ibeg
   if (nch < 9) then
      write(chno,'(i1)') nch+1
   else if (nch < 99) then
      write(chno,'(i2)') nch+1
   else if (nch < 999) then
      write(chno,'(i3)') nch+1
   else
      write(0,*) 'nch is too big: ', nch
      call mexit(6,0)
   end if
   fmt2 = '(' // chno // 'f9.2)'
   call amopen(15,prjfil,owrite,'F','W')
   
   !  ---   read the mass file
   
   natom = ncart/3
   call amopen(16,masfil,'O','F','R')
   read(16,*) (atmas(i),i=1,natom)
   close(16)
   k = 0
   do i=1,natom
      atmas3(k+1) = atmas(i)
      atmas3(k+2) = atmas(i)
      atmas3(k+3) = atmas(i)
      k = k + 3
   end do
   
   !  ---   open the input file
   
   ftyp = 5
   call openinp(ftyp)
   
   !  ---   read input file
   
   ierr=0
   ieof=0
   isnap = 0
   do isp=1,9999999
      call readfile(natom,nread,ener,crdr,vel,box,ftyp,ierr,ieof)
      if (ierr == 1) then
         write(0,*) 'Error during read: natomal=',natom,' nread=',nread
         stop
      end if
      if (ieof == 1) goto 30
      if (isp > last) goto 30
      if (isp >= first .and. mod(isp,iskip) == 0) then
         isnap = isnap + 1
         do imode=ibeg,iend
            proj(imode) = 0.0
            do i=1,3*natom
               proj(imode) = proj(imode) + (crdr(i)-x(i))* &
                     atmas3(i)*vect(i,imode)
               
               !    ----------debug: check that modes are orthonormal:---------------
               !             proj(imode) = proj(imode) + atmas3(i)*vect(i,imode)*
               !    .                        vect(i,isp)
               !    -----------------------------------------------------------------
               
            end do
         end do
         write(15,fmt2) (proj(i),i=ibeg,iend)
      end if
   end do
   
   30 write(6,*) 'Projected ',isnap,' coordinate sets'
   return
end subroutine quasip 
