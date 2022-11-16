!     md file i/o routines

!     input file types:
!     1 - .crd2 file
!     2 - Amber coord file made with NTWX
!     3 - traj  file
!     4 - .xyz file
!     5 - binpos file
!     6 - .crd2 file without velocities


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine openio here]
subroutine openio(ftyp)
   
   !     used for opening input and output coordinate files
   
   integer ftyp
   character(len=80) title,header
   
   if (ftyp == 2) then
      read(5,1000) title
      write(6,1000) title
      write(0,*) 'Title in input file: ',title
   end if
   
   if (ftyp == 4) then
      read(5,1000) title
      write(6,1000) title
      read(5,1000) header
      write(6,1000) header
   end if
   
   if (ftyp == 5) then
      call startbinpos()
      call openbinpos()
   end if
   
   1000 format(a80)
   
   return
end subroutine openio 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine openinp here]
subroutine openinp(ftyp)
   
   !     used for opening input coordinate file only
   
   integer ftyp
   character(len=80) title,header
   
   if (ftyp == 2 .or. ftyp == 3) then
      read(5,1000) title
      write(6,*) 'Title in input file: ',title
   end if
   
   if (ftyp == 4) then
      read(5,1000) title
      read(5,1000) header
   end if
   
   if (ftyp == 5) then
      call openbinpos()
   end if
   
   1000 format(a80)
   
   return
end subroutine openinp 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine readfile here]
subroutine readfile(ntot,nread,ener,crd,vel,box,ftyp, &
      hasbox,ierr,ieof)
   
#  include "fileio.h"
   dimension crd(3*mxatom),vel(3*mxatom),ener(36),box(3)
   integer ntot,nread,ftyp,ierr,ieof
   logical hasbox
   
   nread=0
   ierr=0
   ieof=0
   
   if (ftyp == 1 .or. ftyp == 6) &
         read(5,1010,end=10,err=20) (ener(i),i=1,36)
   if (ftyp /= 5) read(5,1020,end=10,err=20) (crd(i),i=1,3*ntot)
   if (ftyp == 2 .and. hasbox) &
         read(5,1020,end=10,err=20) boxx,boxy,boxz
   if (ftyp == 1) read(5,1020,end=10,err=20) (vel(i),i=1,3*ntot)
   if (ftyp == 4 .and. hasbox) &
         read(5,1000,end=10,err=20) (box(i),i=1,3)
   
   if (ftyp == 5) then
      call readbinpos(nread,crd,ieof)
      if (ieof == 1) goto 10
      if (nread /= ntot) goto 20
   end if
   
   goto 30
   
   10 ieof=1
   goto 30
   
   20 ierr=1
   30 return
   
   1000 format(6f12.7)
   1010 format(6e12.5)
   1020 format(10f8.3)
   
end subroutine readfile 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine writefile here]
subroutine writefile(ntot,ener,crd,vel,box,ftyp,ierr)
   
#  include "fileio.h"
   
   dimension crd(3*mxatom),vel(3*mxatom),ener(36),box(3)
   integer ntot,ftyp,ierr
   
   ierr=0
   if (ftyp == 1 .or. ftyp == 6) &
         write(6,1010,err=10) (ener(i),i=1,36)
   if (ftyp /= 5) write(6,1020,err=10) (crd(i),i=1,3*ntot)
   if (ftyp == 1) write(6,1020,err=10) (vel(i),i=1,3*ntot)
   if (ftyp == 4) write(6,1000,err=10) (box(i),i=1,3)
   if (ftyp == 5) call writebinpos(ntot,crd)
   
   goto 20
   
   10 ierr=1
   20 return
   
   1000 format(6f12.7)
   1010 format(6e12.5)
   1020 format(10f8.3)
   
end subroutine writefile 
