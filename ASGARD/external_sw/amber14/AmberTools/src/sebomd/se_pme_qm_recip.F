      subroutine se_pme_qm_recip(erec,derec,iat,jat)

      implicit none
#include "sebomd.dim"
#include "se_logic.h"
#include "se_pbcbox.h"
#include "se_pmeee.h"
#include "se_xyzcrd.h"
#include "se_constants.h"
      double precision erec
      double precision derec(3)
      integer iat, jat

      double precision eta

      double precision xi, yi, zi
      double precision xj, yj, zj
      double precision fact
      double precision uk1, uk2, uk3
      double precision scal, uksqr, expo
      double precision rsum_tol, rsum
      double precision a1, a2, a3, a4, a5, a6
      double precision uksqr_inv
      integer m1, m2, m3

      eta = betapme
      derec(1)=0.0d0
      derec(2)=0.0d0
      derec(3)=0.0d0

      xi = xyz(1,iat)
      yi = xyz(2,iat)
      zi = xyz(3,iat)
      if (pbc) then
        call se_pbcxyz(iat,jat,xj,yj,zj)
      else
        xj = xyz(1,jat)
        yj = xyz(2,jat)
        zj = xyz(3,jat)
      endif

      erec=0.0d0
      fact=4.0d0*pi/boxvol

      a1 = 2.0d0*pi/dbox(1)
      a2 = 2.0d0*pi/dbox(2)
      a3 = 2.0d0*pi/dbox(3)
      a4 = 1.0d0/(4.0d0*eta*eta)
      do m1=-mmax,mmax
        do m2=-mmax,mmax
          do m3=-mmax,mmax
            uk1=m1*a1
            uk2=m2*a2
            uk3=m3*a3
            if (m1.ne.0.or.m2.ne.0.or.m3.ne.0) then

              scal = uk1*(xi-xj)+uk2*(yi-yj)+uk3*(zi-zj)
              uksqr=uk1*uk1+uk2*uk2+uk3*uk3
              uksqr_inv = 1.0d0/uksqr
              expo=dexp(-uksqr*a4)
              a5 = fact*expo*uksqr_inv
              erec=erec + a5*dcos(scal)
              a6 = a5*dsin(scal)
!             if (iat.eq.jat) then
!               derec(1)=0.0d0
!               derec(2)=0.0d0
!               derec(3)=0.0d0
!             else
                derec(1)=derec(1) - uk1*a6
                derec(2)=derec(2) - uk2*a6
                derec(3)=derec(3) - uk3*a6
!             endif
            endif  
          enddo
        enddo
      enddo

      erec = elec2eV*erec
      derec(1) = elec2eV*derec(1)
      derec(2) = elec2eV*derec(2)
      derec(3) = elec2eV*derec(3)

      end
