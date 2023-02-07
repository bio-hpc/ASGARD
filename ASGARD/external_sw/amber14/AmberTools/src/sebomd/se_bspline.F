!------------------------------------------------------------------------
!     This file contains all the Cardinal B-spline functions needed
!     by the PME routines:
!
!     mbspline1    massive B-Spline algorithm, returns splines
!     mbspline     massive B-Spline algorithm, returns splines + derivatives
!     mkbspline    massive B-Spline algorithm for integers

      subroutine se_mbspline1(upmex,upmey,upmez,xmn,ymn,zmn)

!     massive bspline algorithm
!     returns the Cardinal B-spline Mn(w) at point w
!     for w=u+int(n-u), u+int(n-u-1), .. , u, u-1,
!     .                  u-2, ... , u-n
!     in the x, y and z-direction.
!     [note that n=nspline is the order of the B-Spline]
!
!     Arjan van der Vaart, Oct. 1997
!
!     parameters:
!     *) iat is the atom number. This is needed for the scaled fractional
!     .  coordinate defined by
!     .  u(i) = K(i)*a(i)*r(i), i=1,3
!     .    with K integer representing the grid size
!     .         a the reciprocal box vector
!     .         r the (Cartesian) coordinate of the particle
!     *) xmn will contain the Cardinal B-spline of order n in for
!     .  the x coordinate;
!     .    xmn(1) = Mn(u+int(n-u))
!     .    xmn(int(u)+int(n-u)+1) = Mn(u-n)
!     .  xmn should be n+int(u)+int(n-u)+1-2 <= 2n-1 elements long,
!     .  this extra length of smn is used as workspace
!     *) ymn, zmn see xmn
!
!     note that nspline (= n) is the order of the B-spline,
!     nspline > 2, nspline must be even

      implicit none
#include "sebomd.dim"
#include "se_pmeee.h"

      double precision upmex
      double precision upmey
      double precision upmez
      double precision xmn(*)
      double precision ymn(*)
      double precision zmn(*)

      integer iplusx, iminx, ikx, mx
      double precision wx
      integer iplusy, iminy, iky, my
      double precision wy
      integer iplusz, iminz, ikz, mz
      double precision wz

      double precision x, y, z
      double precision di, dim1

      integer m, i, j

!     x coordinate
      iplusx = int(nspline-upmex)
      iminx = int(upmex)
      ikx = iminx+iplusx
      mx = ikx+1
      wx = upmex+iplusx

!     y coordinate
      iplusy = int(nspline-upmey)
      iminy = int(upmey)
      iky = iminy+iplusy
      my = iky+1
      wy = upmey+iplusy

!     z coordinate
      iplusz = int(nspline-upmez)
      iminz = int(upmez)
      ikz = iminz+iplusz
      mz = ikz+1
      wz = upmez+iplusz

      m = max(mx,my,mz)

!     fill M2
      do i=1,nspline+m-2
         xmn(i) = 0.0
         ymn(i) = 0.0
         zmn(i) = 0.0
      enddo
      xmn(ikx+1) = upmex-iminx
      xmn(ikx) = 1.0-xmn(ikx+1)
      ymn(iky+1) = upmey-iminy
      ymn(iky) = 1.0-ymn(iky+1)
      zmn(ikz+1) = upmez-iminz
      zmn(ikz) = 1.0-zmn(ikz+1)

!     fill M3 through M(nspline)
      do i=3,nspline
         di = dble(i)
         dim1 = di-1.0
         do j=1,nspline-i+mx
            x = wx-(j-1)
            xmn(j) = (x*xmn(j) + (di-x)*xmn(j+1))/dim1
         enddo
         do j=1,nspline-i+my
            y = wy-(j-1)
            ymn(j) = (y*ymn(j) + (di-y)*ymn(j+1))/dim1
         enddo
         do j=1,nspline-i+mz
            z = wz-(j-1)
            zmn(j) = (z*zmn(j) + (di-z)*zmn(j+1))/dim1
         enddo
      enddo

      end
!--------------------------------------------------------------------------------
      subroutine se_mbspline(upmex,upmey,upmez,xmn,ymn,zmn,dxmn,dymn,
     .                       dzmn)

!     massive bspline algorithm
!     returns the derivative of the Cardinal B-spline at point w
!     (d Mn(w)/dw) for w=u+int(n-u), u+int(n-u-1), .. , u, u-1,
!     .                  u-2, ... , u-n
!     in the x, y and z-direction.
!     [note that n=nspline is the order of the B-Spline]
!
!     Arjan van der Vaart, Oct. 1997
!
!     parameters:
!     *) iat is the atom number. This is needed for the scaled fractional
!     .  coordinate defined by
!     .  u(i) = K(i)*a(i)*r(i), i=1,3
!     .    with K integer representing the grid size
!     .         a the reciprocal box vector
!     .         r the (Cartesian) coordinate of the particle
!     *) xmn will contain the Cardinal B-spline of order n in for
!     .  the x coordinate;
!     .    xmn(1) = Mn(u+int(n-u))
!     .    xmn(int(u)+int(n-u)+1) = Mn(u-n)
!     .  xmn should be n+int(u)+int(n-u)+1-2 <= 2n-1 elements long,
!     .  this extra length of smn is used as workspace
!     *) ymn, zmn see xmn
!     *) dxmn contains the derivatives dMn(w)/dw in the x direction;
!     .    dxmn(1) = dMn(u+int(n-u))/d(u+int(n-u))
!     .    dxmn(int(u)+int(n-u)+1) = dMn(u-n)/d(u-n)
!     .  dxmn should be int(u)+int(n-u)+1 <= n+1 elements long.
!     .  analogous for dymn and dzmn
!
!     note that nspline (= n) is the order of the B-spline,
!     nspline > 2, nspline must be even

      implicit none
#include "sebomd.dim"
#include "se_pmeee.h"

!     parameters:
      double precision upmex
      double precision upmey
      double precision upmez
      double precision xmn(*)
      double precision ymn(*)
      double precision zmn(*)
      double precision dxmn(*)
      double precision dymn(*)
      double precision dzmn(*)


      integer iplusx, iminx, ikx, mx
      double precision wx, xnmw
      integer iplusy, iminy, iky, my
      double precision wy, ynmw
      integer iplusz, iminz, ikz, mz
      double precision wz, znmw

      double precision x, y, z
      double precision di, dim1, dnm1

      integer m, i, j, im1

!     x coordinate
      iplusx = int(nspline-upmex)
      iminx = int(upmex)
      ikx = iminx+iplusx
      mx = ikx+1
      wx = upmex+iplusx

!     y coordinate
      iplusy = int(nspline-upmey)
      iminy = int(upmey)
      iky = iminy+iplusy
      my = iky+1
      wy = upmey+iplusy

!     z coordinate
      iplusz = int(nspline-upmez)
      iminz = int(upmez)
      ikz = iminz+iplusz
      mz = ikz+1
      wz = upmez+iplusz

      m = max(mx,my,mz)

!     fill M2
      do i=1,nspline+m-2
         xmn(i) = 0.0
         ymn(i) = 0.0
         zmn(i) = 0.0
      enddo
      xmn(ikx+1) = upmex-iminx
      xmn(ikx) = 1.0-xmn(ikx+1)
      ymn(iky+1) = upmey-iminy
      ymn(iky) = 1.0-ymn(iky+1)
      zmn(ikz+1) = upmez-iminz
      zmn(ikz) = 1.0-zmn(ikz+1)

!     fill M3 through M(nspline-1)
      do i=3,nspline-1
         di = dble(i)
         dim1 = di-1.0
         do j=1,nspline-i+mx
            x = wx-(j-1)
            xmn(j) = (x*xmn(j) + (di-x)*xmn(j+1))/dim1
         enddo
         do j=1,nspline-i+my
            y = wy-(j-1)
            ymn(j) = (y*ymn(j) + (di-y)*ymn(j+1))/dim1
         enddo
         do j=1,nspline-i+mz
            z = wz-(j-1)
            zmn(j) = (z*zmn(j) + (di-z)*zmn(j+1))/dim1
         enddo
      enddo

!     fill dMn and Mn
      dnm1 = dnspline-1.0
      xnmw = dnspline-wx
      do i=1,mx
         im1 = i-1
         dxmn(i) = xmn(i)-xmn(i+1)
         xmn(i) = ((wx-im1)*xmn(i)+(xnmw+im1)*xmn(i+1))/dnm1
      enddo
      ynmw = dnspline-wy
      do i=1,my
         im1 = i-1
         dymn(i) = ymn(i)-ymn(i+1)
         ymn(i) = ((wy-im1)*ymn(i)+(ynmw+im1)*ymn(i+1))/dnm1
      enddo
      znmw = dnspline-wz
      do i=1,mz
         im1 = i-1
         dzmn(i) = zmn(i)-zmn(i+1)
         zmn(i) = ((wz-im1)*zmn(i)+(znmw+im1)*zmn(i+1))/dnm1
      enddo
      end
!--------------------------------------------------------------------------------
      subroutine se_mkbspline(smn)

!     massive bspline algorithm
!     returns the Cardinal B-spline Mn(w) at point k+1 for
!     k=0, 1, .. , n-2
!
!     Arjan van der Vaart, Oct. 1997
!
!     parameters:
!     *) smn will contain the Cardinal B-spline;
!     .    smn(1) = Mn(1)
!     .    smn(k) = Mn(k+1)
!
!     nspline is the order of the B-Spline
!
      implicit none
#include "sebomd.dim"
#include "se_pmeee.h"

      double precision smn(*), wmn(maxspline2m3)

      integer m, n2m1, n2m3, i, nj, j
      double precision xi, xim1, x

      m = nspline-1
      n2m1 = m+nspline
      n2m3 = n2m1-2

!     fill Mn2
      do i=1,n2m3
         wmn(i) = 0.0
      enddo
      wmn(m) = 1.0

!     fill Mn3 through Mn(nspline-1)
      do i=3,m
         xi = dble(i)
         xim1 = xi-1.0
         do j=1,n2m1-i
            x = dnspline-j
            wmn(j) = (x*wmn(j) + (xi-x)*wmn(j+1))/xim1
         enddo
      enddo

!     fill Mn(nspline)
      xim1 = dnspline-1.0
      do j=1,m
         nj = nspline-j
         x = dble(nj)
         smn(nj) = (x*wmn(j) + (dnspline-x)*wmn(j+1))/xim1
      enddo

      end
