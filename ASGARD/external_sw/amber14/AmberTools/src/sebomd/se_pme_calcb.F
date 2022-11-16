      subroutine se_pme_calcb

      implicit none
#include "sebomd.dim"
#include "se_pmeee.h"

!     this routine calculates |1/bi|**2, i=x,y,z as given by [Eq.4.4]
!     of Essmann et.al., J.Chem.Phys. 103 8577

!     written by Arjan van der Vaart

      double precision smn(maxspline)
      double precision pi2k1pme, pi2k2pme, pi2k3pme
      double precision arg, argx, argy, argz
      double precision bxr, bxi
      double precision byr, byi
      double precision bzr, bzi

      integer m3, m2,  m1, k

      ! some constants
      pi2k3pme = pi2/dk3pme
      pi2k2pme = pi2/dk2pme
      pi2k1pme = pi2/dk1pme

      ! do the splining
      call se_mkbspline(smn)

      ! |bz|**2
      do m3=0,k3pmem1
        argz = pi2k3pme*dble(m3)
        bzr = smn(1)
        bzi = 0.0
        do k=1,nspline-2
          arg = argz*dble(k)
          bzr = bzr + smn(k+1)*cos(arg)
          bzi = bzi + smn(k+1)*sin(arg)
        end do
        bfaczpme(m3) = bzr*bzr + bzi*bzi
      end do
 
      ! |by|**2
      do m2=0,k2pmem1
        argy = pi2k2pme*dble(m2)
        byr = smn(1)
        byi = 0.0
        do k=1,nspline-2
          arg = argy*dble(k)
          byr = byr + smn(k+1)*cos(arg)
          byi = byi + smn(k+1)*sin(arg)
        end do
        bfacypme(m2) = byr*byr + byi*byi
      end do

      ! |bx|**2
      do m1=0,k1pmem1
        argx = pi2k1pme*dble(m1)
        bxr = smn(1)
        bxi = 0.0
        do k=1,nspline-2
          arg = argx*dble(k)
          bxr = bxr + smn(k+1)*cos(arg)
          bxi = bxi + smn(k+1)*sin(arg)
        end do
        bfacxpme(m1) = bxr*bxr + bxi*bxi
      end do

      end
