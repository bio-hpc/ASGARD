       subroutine writcrd(uwcrd)

#include "SIZE.h"
#include "TOP.h"
#include "MISC.h"
       integer tempint

       write(uwcrd,105,err=1010)crdline
105	format (a)
c
c write time info if it was read
c
16     format (i6,1x,e14.7)

        if (havetime) then
            write (uwcrd,16,err=1010) natom,time
        else
            write (uwcrd,16,err=1010) natom
        endif
c
       write (*,*) 'Writing coordinates to output file'  
       write (uwcrd,1400,err=1010) (x(i),y(i),z(i),i=1,natom)
1400    format (6(f12.7))
c
       if (rcvdok.or.rcvbok) then
           write (*,*) 'Writing velocities to output file'
           write (uwcrd,1400,err=1010) (vx(i),vy(i),vz(i),i=1,natom)
        endif
        if (rcvbok.or.rcbdok) then
         if (ewald_box) then
           write (*,*) 'Writing box coords and angles to output file'
           write (uwcrd,1400,err=1010) (crdbox(i),i=1,3),
     &     ew_alpha,ew_beta,ew_gamma
         else
           write (*,*) 'Writing box coords to output file'
           write (uwcrd,1400,err=1010) (crdbox(i),i=1,3)
         endif
        endif

       return
c
1010	write (*,*) 'coordinate file write error'
       stop	
       end


