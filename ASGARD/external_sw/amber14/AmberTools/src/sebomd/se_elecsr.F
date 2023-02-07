      subroutine se_elecsr(eelecsr,deelecsr,iat,jat)

      implicit none
#include "sebomd.dim"
#include "se_logic.h"
#include "se_xyzcrd.h"
#include "se_constants.h"
      double precision deelecsr(3) 
      double precision eelecsr
      integer iat, jat

      integer i, iat3, jat3
      double precision xi, yi, zi
      double precision xj, yj, zj
      double precision x, y, z
      double precision dx, dy, dz
      double precision rsqr, rr, rr3

      eelecsr=0.0d0
      do i=1,3
        deelecsr(i)=0.0d0
      enddo
      if (iat.ne.jat) then
        iat3 = 3*iat-2
        xi = xyz(1,iat)
        yi = xyz(2,iat)
        zi = xyz(3,iat)
        jat3 = 3*jat-2
        if (pbc) then
          call se_pbcxyz(iat,jat,xj,yj,zj)
        else
          xj = xyz(1,jat)
          yj = xyz(2,jat)
          zj = xyz(3,jat)
        endif
        x = xi-xj
        y = yi-yj
        z = zi-zj

        rsqr = x*x + y*y + z*z
        rr = dsqrt(rsqr)
        rr3 = rr**3

        dx = x/rr3
        dy = y/rr3
        dz = z/rr3
!       eelecsr= - 14.39964521/rr
!       deelecsr(1)=14.39964521*dx
!       deelecsr(2)=14.39964521*dy
!       deelecsr(3)=14.39964521*dz
        eelecsr= - elec2eV/rr
        deelecsr(1)=elec2eV*dx
        deelecsr(2)=elec2eV*dy
        deelecsr(3)=elec2eV*dz
      endif
      end
