      subroutine se_dodmxinter()

      implicit none

#include "sebomd.dim"
#include "se_bpairs.h"
#include "se_fulsys.h"
#include "se_glbmat1.h"
#include "se_glbmat2.h"
#include "se_glbmat4.h"

      double precision xiia(3),yiia(3)
      double precision x
      integer n
      integer i, ii, iatm, ijpair, jatm, ijpair3, ijpair2, ijpair1

       n = 3 ! interpolation order
       xiia(1)=1*1.0d0
       xiia(2)=2*1.0d0
       xiia(3)=3*1.0d0
       x=4.0d0

       ! interpolation of the diagonal block
       do i=1,iimat(natoms+1)-1
         yiia(1)=pdiag3(i)
         yiia(2)=pdiag2(i)
         yiia(3)=pdiag1(i)
         call se_polint(xiia,yiia,n,x,pdiag(i))
       enddo

       ! interpolation of the non diagonal part
       do iatm = 2, natoms
         do ijpair = ip1(iatm), ip1(iatm+1)-1
           jatm = ipair(ijpair)
           ! ijpair = interaction between atom iatm and atom jatm
           ! -> does it exists in previous pairs?

           ! in ipair3
           ijpair3 = -1
           do i = ip13(iatm), ip13(iatm+1)-1
              if (ipair3(i).eq.jatm) then
                ijpair3 = i
                exit
              endif
           end do

           ! in ipair2
           ijpair2 = -1
           do i = ip12(iatm), ip12(iatm+1)-1
              if (ipair2(i).eq.jatm) then
                ijpair2 = i
                exit
              endif
           end do

           ! in ipair1
           ijpair1 = -1
           do i = ip11(iatm), ip11(iatm+1)-1
              if (ipair1(i).eq.jatm) then
                ijpair1 = i
                exit
              endif
           end do

           do ii = 0,ijmat(ijpair+1)-ijmat(ijpair)-1
             if (ijpair3.gt.0) then
               yiia(1)=pdiat3(ijmat3(ijpair3)+ii)
             else
               yiia(1) = 0.0d0
             endif
             if (ijpair2.gt.0) then
               yiia(2)=pdiat2(ijmat2(ijpair2)+ii)
             else
               yiia(2) = 0.0d0
             endif
             if (ijpair1.gt.0) then
               yiia(3)=pdiat1(ijmat1(ijpair1)+ii)
             else
               yiia(3) = 0.0d0
             endif
             call se_polint(xiia,yiia,n,x,pdiat(ijmat(ijpair)+ii))
           end do
         end do
       end do
       END

!_______________________________________________________________________

      subroutine se_polint(xa,ya,n,x,y)

      implicit none
      integer i,j,n
      double precision xa(n),ya(n),x,y,yi
     

       y=0.0d0
      do i=1,n
        yi=ya(i)
        do j=1,n 
          if (j.ne.i) then
          yi=yi*(x-xa(j))/(xa(i)-xa(j)) 
          endif  
        enddo
       y=y+yi
      enddo

      end

!_______________________________________________________________________
      subroutine se_updatedmxinter()
      implicit none
#include "sebomd.dim"
#include "se_bpairs.h"
#include "se_fulsys.h"
#include "se_glbmat1.h"
#include "se_glbmat2.h"
#include "se_glbmat4.h"

      integer ii
      ! save ip1
      do ii=1,natoms+1
        ip13(ii) = ip12(ii)
      end do
      ! save ipair and ijmat
      do ii=1,ip13(natoms+1)
        ipair3(ii) = ipair2(ii)
        ijmat3(ii) = ijmat2(ii)
      end do

      do ii=1,natoms+1
        ip12(ii) = ip11(ii)
      end do
      do ii=1,ip12(natoms+1)
        ipair2(ii) = ipair1(ii)
        ijmat2(ii) = ijmat1(ii)
      end do

      do ii=1,natoms+1
        ip11(ii) = ip1(ii)
      end do
      do ii=1,ip11(natoms+1)
        ipair1(ii) = ipair(ii)
        ijmat1(ii) = ijmat(ii)
      end do

      ! save pdiag
      do ii=1,iimat(natoms+1)
        pdiag3(ii) = pdiag2(ii)
        pdiag2(ii) = pdiag1(ii)
        pdiag1(ii) = pdiag(ii)
      enddo

      ! save pdiat
      do ii = 1, ijmat3(ip13(natoms+1))
        pdiat3(ii) = pdiat2(ii)
      end do
      do ii = 1, ijmat2(ip12(natoms+1))
        pdiat2(ii) = pdiat1(ii)
      end do
      do ii = 1, ijmat1(ip11(natoms+1))
        pdiat1(ii) = pdiat(ii)
      end do
      return
      end
