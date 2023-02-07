      subroutine se_pme_derec(derec)

      implicit none
#include "sebomd.dim"
#include "se_pmeee.h"
#include "se_fulsys.h"
#include "se_xyzcrd.h"
#include "se_pbcbox.h"
#include "se_esave.h"
#include "se_constants.h"

!     this routine calculates the derivative of the reciprocal energy
!     to the coordinates.
!
!     in order to save memory, dQ(k1,k2,k3)/dr is not calculated in
!     pme_calcq, but calculated here where Conv(theta,Q) is known.
!
!
!     Written by Arjan van der Vaart, Dec. '97
!
!     parameters:
!
!     *) derec, a wrapped-around real array, will be filled with
!     .         the derivative dErec/dr
!
!     Note that the performance of this algorithm is very high
!     by demanding that nspline <= K_i, i=1,2,3 and by
!     shifting the scaled fractional coordinates upme such that
!     0 <= upme_i <= K_i, i=1,2,3.
!
!     derec needs to be recalculated whenever the coordinates or
!     charges change.
!

      double precision derec(maxpar)

      logical ly1,ly2,lx1,lx2,lx,lyx1,lyx2

      double precision xmn(maxspline2p1),
     .                 ymn(maxspline2p1),
     .                 zmn(maxspline2p1),
     .                 dxmn(maxsplinep1),
     .                 dymn(maxsplinep1),
     .                 dzmn(maxsplinep1)
      double precision f1
      double precision f2
      double precision f3
      integer iat
      integer i1
      integer i2
      integer i3
      integer j
      integer iupmex
      integer iupmey
      integer iupmez
      integer iw, ixp1, ixpk1pmem1
      integer kxe1, kxb1, kxb2, ixp, kx, ix
      integer kye1, kyb1, kyb2, iyp, ky, iy, iwy
      integer kze1, kzb1, kzb2, izp, kz, iz, iwz
      double precision x, upmex
      double precision y, upmey
      double precision z, upmez
      double precision dz0, dy, dz
      

!----------------------------------------------------------

      f1 = elec2kcal*dk1pme*recip1(1)
      f2 = elec2kcal*dk2pme*recip2(2)
      f3 = elec2kcal*dk3pme*recip3(3)

      do iat=1,natoms

         i3 = 3*iat
         i2 = i3-1
         i1 = i2-1

!     .  calculate scaled fractional coordinates
         x = 0.0
         y = 0.0
         z = 0.0
         do j=1,3
            x = x + recip1(j)*(xyz(j,iat)-xyzmin(j))
            y = y + recip2(j)*(xyz(j,iat)-xyzmin(j))
            z = z + recip3(j)*(xyz(j,iat)-xyzmin(j))
         enddo
         upmex = dk1pme*x
         iupmex = int(upmex-nspline)
         upmey = dk2pme*y
         iupmey = int(upmey-nspline)
         upmez = dk3pme*z
         iupmez = int(upmez-nspline)

!     .  calculate the cardinal B-Splines Mn(upme-k) and derivatives.
!     .  Only one call is needed to fill these in the x,y and z-direction
!     .  for all the values of k for which Mn(upme-k) is NOT zero.

         call se_mbspline(upmex,upmey,upmez,xmn,ymn,zmn,dxmn,dymn,dzmn)

!     .  calculate the k values for which Mn(u-k-nK) is NOT zero.
!     .  n = 0 for k=kxb1,kxe1 (same in y and z direction)
!     .  n = 1 for k=kxb2,kxe2 (same in y and z direction)
!     .  These are the only possible values of n (since 0 <= upme <= K
!     .  and nspline <= K).

         kxe1 = int(upmex)
         kxb2 = kxe1+k1pmenspl
         kxe1 = min(kxe1,k1pmem1)
         kxb1 = max(izero,iupmex)
         ixp = int(nspline-upmex)
         ixp1 = ixp+1
         ixpk1pmem1 = ixp-k1pmem1

         kye1 = int(upmey)
         kyb2 = kye1+k2pmenspl
         kye1 = min(kye1,k2pmem1)
         kyb1 = max(izero,iupmey)
         iyp = int(nspline-upmey)

         kze1 = int(upmez)
         kzb2 = kze1+k3pmenspl
         kze1 = min(kze1,k3pmem1)
         kzb1 = max(izero,iupmez)
         izp = int(nspline-upmez)

!     .  computation time can be safed if the program doesn't
!     .  have to go in triple loops of which the inner loops
!     .  are empty

         lx1 = kxb1.le.kxe1
         lx2 = kxb2.le.k1pmem1
         lx = lx1.or.lx2

         ly1 = kyb1.le.kye1
         ly2 = kyb2.le.k2pmem1
         lyx1 = ly1.and.lx
         lyx2 = ly2.and.lx

!     .  in the following the derivatives dQ(kx,ky,kz)/di, i=x,y,z
!     .  are performed for Mn_i .ne. 0 (i=x,y,z)

         if (lyx1) then
            do kz=kzb1,kze1
               iz = izp+1+kz
               z = chgpme(iat)*zmn(iz)
               dz0 = chgpme(iat)*dzmn(iz)
               iwz = kz*k1pmek2pme+1
               do ky=kyb1,kye1
                  iy =  iyp+1+ky
                  y = z*ymn(iy)
                  dy = z*dymn(iy)
                  dz = dz0*ymn(iy)
                  iwy = iwz + ky*k1pme
                  do kx=kxb1,kxe1
                     iw = 2*(iwy + kx) - 1
                     ix = ixp1+kx
                     derec(i1) = derec(i1) + qpmec(iw)*y*dxmn(ix)
                     derec(i2) = derec(i2) + qpmec(iw)*dy*xmn(ix)
                     derec(i3) = derec(i3) + qpmec(iw)*dz*xmn(ix)
                  enddo
                  do kx=kxb2,k1pmem1
                     iw = 2*(iwy + kx) - 1
                     ix = ixpk1pmem1+kx
                     derec(i1) = derec(i1) + qpmec(iw)*y*dxmn(ix)
                     derec(i2) = derec(i2) + qpmec(iw)*dy*xmn(ix)
                     derec(i3) = derec(i3) + qpmec(iw)*dz*xmn(ix)
                  enddo
               enddo
            enddo
         endif

         if (lyx2) then
            do kz=kzb1,kze1
               iz = izp+1+kz
               z = chgpme(iat)*zmn(iz)
               dz0 = chgpme(iat)*dzmn(iz)
               iwz = kz*k1pmek2pme+1
               do ky=kyb2,k2pmem1
                  iy = iyp+ky-k2pmem1
                  y = z*ymn(iy)
                  dy = z*dymn(iy)
                  dz = dz0*ymn(iy)
                  iwy = iwz + ky*k1pme
                  do kx=kxb1,kxe1
                     iw = 2*(iwy + kx) - 1
                     ix = ixp1+kx
                     derec(i1) = derec(i1) + qpmec(iw)*y*dxmn(ix)
                     derec(i2) = derec(i2) + qpmec(iw)*dy*xmn(ix)
                     derec(i3) = derec(i3) + qpmec(iw)*dz*xmn(ix)
                  enddo
                  do kx=kxb2,k1pmem1
                     iw = 2*(iwy + kx) - 1
                     ix = ixpk1pmem1+kx
                     derec(i1) = derec(i1) + qpmec(iw)*y*dxmn(ix)
                     derec(i2) = derec(i2) + qpmec(iw)*dy*xmn(ix)
                     derec(i3) = derec(i3) + qpmec(iw)*dz*xmn(ix)
                  enddo
               enddo
            enddo
         endif

         if (lyx1) then
            do kz=kzb2,k3pmem1
               iz = izp+kz-k3pmem1
               z = chgpme(iat)*zmn(iz)
               dz0 = chgpme(iat)*dzmn(iz)
               iwz = kz*k1pmek2pme+1
               do ky=kyb1,kye1
                  iy = iyp+1+ky
                  y = z*ymn(iy)
                  dy = z*dymn(iy)
                  dz = dz0*ymn(iy)
                  iwy = iwz + ky*k1pme
                  do kx=kxb1,kxe1
                     iw = 2*(iwy + kx) - 1
                     ix = ixp1+kx
                     derec(i1) = derec(i1) + qpmec(iw)*y*dxmn(ix)
                     derec(i2) = derec(i2) + qpmec(iw)*dy*xmn(ix)
                     derec(i3) = derec(i3) + qpmec(iw)*dz*xmn(ix)
                  enddo
                  do kx=kxb2,k1pmem1
                     iw = 2*(iwy + kx) - 1
                     ix = ixpk1pmem1+kx
                     derec(i1) = derec(i1) + qpmec(iw)*y*dxmn(ix)
                     derec(i2) = derec(i2) + qpmec(iw)*dy*xmn(ix)
                     derec(i3) = derec(i3) + qpmec(iw)*dz*xmn(ix)
                  enddo
               enddo
            enddo
         endif

         if (lyx2) then
            do kz=kzb2,k3pmem1
               iz = izp+kz-k3pmem1
               z = chgpme(iat)*zmn(iz)
               dz0 = chgpme(iat)*dzmn(iz)
               iwz = kz*k1pmek2pme+1
               do ky=kyb2,k2pmem1
                  iy = iyp+ky-k2pmem1
                  y = z*ymn(iy)
                  dy = z*dymn(iy)
                  dz = dz0*ymn(iy)
                  iwy = iwz + ky*k1pme
                  do kx=kxb1,kxe1
                     iw = 2*(iwy + kx) - 1
                     ix = ixp1+kx
                     derec(i1) = derec(i1) + qpmec(iw)*y*dxmn(ix)
                     derec(i2) = derec(i2) + qpmec(iw)*dy*xmn(ix)
                     derec(i3) = derec(i3) + qpmec(iw)*dz*xmn(ix)
                  enddo
                  do kx=kxb2,k1pmem1
                     iw = 2*(iwy + kx)-1
                     ix = ixpk1pmem1+kx
                     derec(i1) = derec(i1) + qpmec(iw)*y*dxmn(ix)
                     derec(i2) = derec(i2) + qpmec(iw)*dy*xmn(ix)
                     derec(i3) = derec(i3) + qpmec(iw)*dz*xmn(ix)
                  enddo
               enddo
            enddo
         endif

         derec(i1) = f1*derec(i1)
         derec(i2) = f2*derec(i2)
         derec(i3) = f3*derec(i3)

      enddo

      return
      end
