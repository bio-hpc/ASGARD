      subroutine se_pme_calctheta

      implicit none
#include "sebomd.dim"
#include "se_pmeee.h"
#include "se_pbcbox.h"

!     this routine calculates the array thetapme = F(C.B)
!     as defined in J.Chem.Phys. 103 ('95) 8577
!
!     thetapme is needed in the calculation of the reciprocal
!     energy:
!     .   Erec = 0.5 Sum{m} [ Q Conv(thetapme . Q) ]    [Eq.4.7]
!
!     thetapme is stored in wrapped-around order.
!
!     this routine only needs to be called once at the start of
!     the simulation if k1pme,k2pme,k3pme are not dependent on
!     the boxdimension (gridsize varies).
!
!     Written by Arjan van der Vaart, Oct. '97
!
!
      double precision
     &     vkx(3,0:maxkpme), vky(3,0:maxkpme),
     &     vkz(3,0:maxkpme), vks(3)

      double precision pib, piv, x, y, z, vk, vk2, v, b
      integer k1h, k2h, k3h, iw, iwz
      integer mx, my, mz, k, iz, iy

!     some constants
      pib = -pisqr/betapme2
      piv = pi*boxvol
      k3h = k3pme/2
      k2h = k2pme/2
      k1h = k1pme/2

      do mz=0,k3h
         do k=1,3
            vkz(k,mz) = dble(mz)*recip3(k)
         enddo
      enddo
      do mz=k3h+1,k3pmem1
         do k=1,3
            vkz(k,mz) = dble(mz-k3pme)*recip3(k)
         enddo
      enddo
      do my=0,k2h
         do k=1,3
            vky(k,my) = dble(my)*recip2(k)
         enddo
      enddo
      do my=k2h+1,k2pmem1
         do k=1,3
            vky(k,my) = dble(my-k2pme)*recip2(k)
         enddo
      enddo
      do mx=0,k1h
         do k=1,3
            vkx(k,mx) = dble(mx)*recip1(k)
         enddo
      enddo
      do mx=k1h+1,k1pmem1
         do k=1,3
            vkx(k,mx) = dble(mx-k1pme)*recip1(k)
         enddo
      enddo

!     (0,0,0)
      thetapme(1) = 0.0

!     (mx,0,0)
      z = bfaczpme(0)
      y = piv*z*bfacypme(0)
      do mx=1,k1pmem1
         vk2 = 0.0
         do k=1,3
            vk2 = vk2 + vkx(k,mx)*vkx(k,mx)
         enddo
         v = exp(pib*vk2)
         b = vk2*y*bfacxpme(mx)
         thetapme(mx+1) = v/b
      enddo

!     (mx,my,0)
      do my=1,k2pmem1
         y = piv*z*bfacypme(my)
         iw = my*k1pme+1
         do mx=1,k1pmem1
            vk2 = 0.0
            do k=1,3
               vk = vky(k,my) + vkx(k,mx)
               vk2 = vk2 + vk*vk
            enddo
            v = exp(pib*vk2)
            b = vk2*y*bfacxpme(mx)
            thetapme(mx+iw) = v/b
         enddo
      enddo

!     (0,0,mz), (mx,0,mz), (0,my,mz), (mx,my,mz)
      do mz=1,k3pmem1
         z = piv*bfaczpme(mz)
         iwz = mz*k1pmek2pme + 1
         do my=0,k2pmem1
            y = z*bfacypme(my)
            iw = iwz + my*k1pme
            do k=1,3
               vks(k) = vkz(k,mz) + vky(k,my)
            enddo
            do mx=0,k1pmem1
               vk2 = 0.0
               do k=1,3
                  vk = vks(k) + vkx(k,mx)
                  vk2 = vk2 + vk*vk
               enddo
               v = exp(pib*vk2)
               b = vk2*y*bfacxpme(mx)
               thetapme(mx+iw) = v/b
            enddo
         enddo
      enddo

!     (0,my,0)
      x = piv*bfaczpme(0)*bfacxpme(0)
      do my=1,k2pmem1
         vk2 = 0.0
         do k=1,3
            vk2 = vk2 + vky(k,my)*vky(k,my)
         enddo
         v = exp(pib*vk2)
         b = vk2*x*bfacypme(my)
         thetapme(my*k1pme+1) = v/b
      enddo

!     The PME-approximation fails when 2|mi|=Kipme, so set the value
!     of thetapme at these mi values to zero

!     (k1pme/2,y,z)
      k1h = k1h + 1
      do mz=0,k3pmem1
         iz = mz*k1pmek2pme + k1h
         do my=0,k2pmem1
            iy = iz + my*k1pme
            thetapme(iy) = 0.0
         enddo
      enddo

!     (x,k2pme/2,z)
      k2h = k2h*k1pme + 1
      do mz=0,k3pmem1
         iz = mz*k1pmek2pme + k2h
         do mx=0,k1pmem1
            thetapme(iz+mx) = 0.0
         enddo
      enddo

!     (x,y,k3pme/2)
      k3h = k3h*k1pmek2pme + 1
      do my=0,k2pmem1
         iy = my*k1pme + k3h
         do mx=0,k1pmem1
            thetapme(iy+mx) = 0.0
         enddo
      enddo

      end


