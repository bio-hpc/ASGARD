      subroutine se_pme_qm_direct(edir,eself,eclmb,dedir,iat,jat)

      implicit none
#include "sebomd.dim"
#include "se_logic.h"
#include "se_pmeee.h"
#include "se_pbcbox.h"
#include "se_xyzcrd.h"
#include "se_constants.h"

      double precision edir, eself, eclmb, dedir(3)
      integer iat, jat

      double precision qp, a1, a2, a3, a4, a5, eta
      parameter (qp=0.3275911d0)
      parameter (a1=0.254829592d0)
      parameter (a2=-0.284496736d0)
      parameter (a3=1.421413741d00)
      parameter (a4=-1.453152027d0)
      parameter (a5=1.061405429d00)
!     parameter (eta=0.70d0)

      integer i, iat3, jat3, nn1, nn2, nn3
      double precision xi, yi, zi, chi
      double precision xj, yj, zj
      double precision x, y, z
      double precision dx, dy, dz
      double precision rsqr, rr, rr3, ch, chk
      double precision etar, expar2, qt, erfc, erf, dfe

!
!     this routine returns the direct energy, the selfenergy and
!     the short-range (classical) Coulomb energy and gradients.
!
!
!     Error Function ( erf(x)= 1-erfc(x)) { erfc(x) is obtained from the
!     polynomial expression}. The term expar2 is the expotential appearing
!     in the Integral appearing in the exact expression for erf.
!     { erf(x) = (2/sqrt(pi))*Integ(exp(-t*t)dt)}
!     The polynomial expression for erfc(x) makes use of the same exponental.
!     However, the derivative is obtained from the actual expression, not the
!     polynomial.
!
!     d(erf(ax)/x) = -erf(ax) + (2/sqrt(pi)) {ax*exp(-ax*ax)}
!     ------------   ----------------------------------------
!     dx                               x*x
!
!

        
!        eta = betapme/6.0d0
      eta = betapme
!     units of eclmb, edir and eself are eV
      eclmb = 0.0
      edir = 0.0
      eself = 0.0
      do i=1,3
        dedir(i)=0.d0
      enddo

!     unit of gradient (dedir) is kcal/A
!     tokcal = 332.0704
!     bpisqrt = betapme/pisqrt
                
      iat3 = 3*iat-2
      xi = xyz(1,iat)
      yi = xyz(2,iat)
      zi = xyz(3,iat)
      chi = 1.0d0
    

!     selfenergy (is a constant, it doesn't contribute to the gradient)
!     eself = eself - bpisqrt*chi*chi

       
      jat3 = 3*jat-2
      if (pbc) then
        call se_pbcxyz(iat,jat,xj,yj,zj)
      else
        xj = xyz(1,jat)
        yj = xyz(2,jat)
        zj = xyz(3,jat)
      endif
!     do nn1=-1,1
!       do nn2=-1,1
!         do nn3=-1,1
! loops on neighbours are not needed since betapme is calculated
! to cancel out interactions with with them
      nn1 = 0
      nn2 = 0
      nn3 = 0
            if ((iat.eq.jat).and.
     .          ((abs(nn1)+abs(nn2)+abs(nn3)).eq.0)) then
              edir=edir - 2.0d0*eta/dsqrt(pi)
            else
              x = xi-xj - nn1*dbox(1)
              y = yi-yj - nn2*dbox(2)
              z = zi-zj - nn3*dbox(3)
              rsqr = x*x + y*y + z*z
              rr = dsqrt(rsqr)
              rr3 = rr**3
              ch = 1.0d0
              chk = 1.0d0
!             short range Coulomb energy
!             eclmb = eclmb + ch/rr
!             calculate contribution to gradient of short range Coulomb
              dx = x/rr3
              dy = y/rr3
              dz = z/rr3
              etar = eta*rr
              expar2 = exp(-etar*etar)
              qt     = 1.0/(1.0+qp*etar)
              erfc = ((((a5*qt+a4)*qt+a3)*qt+a2)*qt+a1)*qt*expar2
              erf = 1.00D0 - erfc

!             direct energy
              edir = edir + erfc/rr

!             calculate contribution to gradient of direct energy;
!             dfe is gradient/|r| of direct energy
              dfe = (-twopisqrt*etar*expar2+erf-1)/(rr*rsqr)

!             add direct and subtract coulomb contributions to gradient:
!             dEclmb/dr_i = SUM(j) -q_i*q_j*r_ij/|r_ij|^3
!             but here dEclmb/dr_i is used as a correction,
!             so change the sign.

              dedir(1) = dedir(1) + x*dfe
              dedir(2) = dedir(2) + y*dfe
              dedir(3) = dedir(3) + z*dfe
              if ((abs(nn1)+abs(nn2)+abs(nn3)).eq.0) then
                dedir(1) = dedir(1) + dx
                dedir(2) = dedir(2) + dy
                dedir(3) = dedir(3) + dz
                edir = edir - 1/rr
              endif
            endif
!         enddo
!       enddo
!     enddo

!     dedir(4) = - x*dfe - dx
!     dedir(5) = - y*dfe - dy
!     dedir(6) = - z*dfe - dz

!     energy conversion: until now energies are in e*e/A
!     multiply this by 1/(4 pi eps0) and then convert to V

      edir = elec2eV*edir
      dedir(1) = elec2eV*dedir(1)
      dedir(2) = elec2eV*dedir(2)
      dedir(3) = elec2eV*dedir(3)
!     edir = 14.39964521*edir
! 128 eself = 14.39964521*eself
!     eclmb = 14.39964521*eclmb
!     do i=1,3*natoms
!       dedir(i)=14.39964521* dedir(i)
!     enddo
      end
