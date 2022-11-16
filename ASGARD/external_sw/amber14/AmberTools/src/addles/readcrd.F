      subroutine readcrd(urcrd)
      implicit none
#include "SIZE.h"
#include "TOP.h"
#include "MISC.h"
      integer tempint, ipack
      character*80 line_test

c
c read title
c
      do ipack = 1, npack
         read(urcrd,105,err=1010)crdline
105      format (a)
         write (6,105) crdline
c
c try to read time, if error then read natom only
c
          time=0.
          read(urcrd,'(a80)') line_test
          if( line_test(6:6).eq.' ' ) then ! this is an old, i5 file
             read(line_test,'(i5,e15.7)',err=18) tempint,time
          else                   ! assume a new, i6 file
             read(line_test,'(i6,e15.7)',err=18) tempint,time
          end if
          havetime=.true.
          go to 19
c
18        write (6,*) "No time index present"
          havetime=.false.
          read(urcrd,'(a80)') line_test
          if( line_test(6:6).eq.' ' ) then ! this is an old, i5 file
             read(line_test,'(i5,e15.7)') tempint
          else                   ! assume a new, i6 file
             read(line_test,'(i6,e15.7)') tempint
          end if
c
19        continue
c
          if (tempint.ne.natom) then
             write(6,*) 'prm/crd natom mismatch, stopping'
             stop
          endif
c
          tempint = (ipack-1)*natom
          write(6,*) 'Reading coordinates from input file'
          read (urcrd,1400,end=1000,err=1010) (x(i),y(i),z(i), 
     &         i=tempint+1, tempint+natom)
1400      format (6(f12.7))
c
c read velos if requested
c
          if (rcvdok.or.rcvbok) then
             write(6,*) 'Reading velocities from input file'
             read (urcrd,1400,end=1000,err=1010) 
     &           (vx(i),vy(i),vz(i),i=tempint+1,tempint+natom)
          endif

          if (rcbdok.or.rcvbok) then
             write(6,*) 'Reading box coords from input file'
c
c try to read PME box angles if present
c
             read (urcrd,1400,end=22,err=22) 
     &          (crdbox(i),i=1,3),ew_alpha,ew_beta,ew_gamma
                ewald_box=.true.
             go to 23

22           backspace (urcrd)
             read (urcrd,1400,end=1000,err=1010) 
     &          (crdbox(i),i=1,3)
             ewald_box=.false.

23           continue
          endif
       end do
 
      return
c
1000   write(6,*) 'crd file end'
       stop        
c
1010   write(6,*) 'crd file error'
       stop
       end


