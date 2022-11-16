      subroutine se_pme_calcq

      implicit none
#include "sebomd.dim"
#include "se_pmeee.h"
#include "se_fulsys.h"
#include "se_xyzcrd.h"
#include "se_esave.h"
#include "se_pbcbox.h"

!     this routine calculates Q as defined in Eq.4.6 of
!     J.Chem.Phys. 103 ('95) 8577.
!
!     Written by Arjan van der Vaart, Dec. '97
!
!
!     Note that the performance of this algorithm is very high
!     by demanding that nspline <= K_i, i=1,2,3 and by
!     shifting the scaled fractional coordinates upme such that
!     0 <= upme_i <= K_i, i=1,2,3.
!
!     Q needs to be recalculated whenever the coordinates or
!     charges change.
!
!

      logical ly1,ly2,lx1,lx2,lx,lyx1,lyx2

      double precision xmn(maxspline2p1),
     &     ymn(maxspline2p1),
     &     zmn(maxspline2p1)

      integer i, iat, ix, iw, iwy, iz, iy, iwz
      integer j
      double precision x, y, z
      double precision upmex, upmey, upmez
      integer kx, ky, kz
      integer kxe1, kxb1, kxb2, ixp, ixp1, ixpk1pmem1
      integer kye1, kyb1, kyb2, iyp
      integer kze1, kzb1, kzb2, izp
      integer iupmex, iupmey, iupmez

!     qpme(maxkpme3)   /pme/  |  ww(msorb2)    /work/   | ww(maxovlp2)
!

!     no sharing of memory of the pme arrays with the /work/ arrays

!     initialize
      do i=1,k123pme
         qpme(i) = 0.0
      enddo

      do iat=1,natoms

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

!     .  calculate the cardinal B-Splines Mn(upme-k). Only one call
!     .  is needed to fill these in the x,y and z-direction for
!     .  all the values of k for which Mn(upme-k) is NOT zero.

         call se_mbspline1(upmex,upmey,upmez,xmn,ymn,zmn)

!     .  calculate the k values for which Mn(u-k-nK) is NOT zero.
!     .  n = 0  for k=kxb1,kxe1 (same in y and z direction)
!     .  n = -1 for k=kxb2,kxe2 (same in y and z direction)
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

!     .  in the following the summation Q = q_i*Mn_x*Mn_y*Mn_z
!     .  is performed for Mn_i .ne. 0 (i=x,y,z)

         if (lyx1) then
            do kz=kzb1,kze1
               iz = izp+1+kz
               z = chgpme(iat)*zmn(iz)
               iwz = kz*k1pmek2pme+1
               do ky=kyb1,kye1
                  iy =  iyp+1+ky
                  y = z*ymn(iy)
                  iwy = iwz + ky*k1pme
                  do kx=kxb1,kxe1
                     iw = iwy + kx
                     ix = ixp1+kx
                     qpme(iw) = qpme(iw) + y*xmn(ix)
                  enddo
                  do kx=kxb2,k1pmem1
                     iw = iwy + kx
                     ix = ixpk1pmem1+kx
                     qpme(iw) = qpme(iw) + y*xmn(ix)
                  enddo
               enddo
            enddo
         endif

         if (lyx2) then
            do kz=kzb1,kze1
               iz = izp+1+kz
               z = chgpme(iat)*zmn(iz)
               iwz = kz*k1pmek2pme+1
               do ky=kyb2,k2pmem1
                  iy = iyp+ky-k2pmem1
                  y = z*ymn(iy)
                  iwy = iwz + ky*k1pme
                  do kx=kxb1,kxe1
                     iw = iwy + kx
                     ix = ixp1+kx
                     qpme(iw) = qpme(iw) + y*xmn(ix)
                  enddo
                  do kx=kxb2,k1pmem1
                     iw = iwy + kx
                     ix = ixpk1pmem1+kx
                     qpme(iw) = qpme(iw) + y*xmn(ix)
                  enddo
               enddo
            enddo
         endif

         if (lyx1) then
            do kz=kzb2,k3pmem1
               iz = izp+kz-k3pmem1
               z = chgpme(iat)*zmn(iz)
               iwz = kz*k1pmek2pme+1
               do ky=kyb1,kye1
                  iy = iyp+1+ky
                  y = z*ymn(iy)
                  iwy = iwz + ky*k1pme
                  do kx=kxb1,kxe1
                     iw = iwy + kx
                     ix = ixp1+kx
                     qpme(iw) = qpme(iw) + y*xmn(ix)
                  enddo
                  do kx=kxb2,k1pmem1
                     iw = iwy + kx
                     ix = ixpk1pmem1+kx
                     qpme(iw) = qpme(iw) + y*xmn(ix)
                  enddo
               enddo
            enddo
         endif

         if (lyx2) then
            do kz=kzb2,k3pmem1
               iz = izp+kz-k3pmem1
               z = chgpme(iat)*zmn(iz)
               iwz = kz*k1pmek2pme+1
               do ky=kyb2,k2pmem1
                  iy = iyp+ky-k2pmem1
                  y = z*ymn(iy)
                  iwy = iwz + ky*k1pme
                  do kx=kxb1,kxe1
                     iw = iwy + kx
                     ix = ixp1+kx
                     qpme(iw) = qpme(iw) + y*xmn(ix)
                  enddo
                  do kx=kxb2,k1pmem1
                     iw = iwy + kx
                     ix = ixpk1pmem1+kx
                     qpme(iw) = qpme(iw) + y*xmn(ix)
                  enddo
               enddo
            enddo
         endif
      enddo
      end
