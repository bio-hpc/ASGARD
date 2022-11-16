      subroutine se_pme_direct2(edir,eself,eclmb,dedir)

      implicit none
#include "sebomd.dim"
#include "se_logic.h"
#include "se_esave.h"
#include "se_pmeee.h"
#include "se_fulsys.h"
#include "se_xyzcrd.h"
#include "se_pbcbox.h"
#include "se_constants.h"

      double precision edir, eself, eclmb, dedir(maxpar)

      double precision qp, a1, a2, a3, a4, a5
      double precision eta
      integer i, iat, iat3
      integer    jat, jat3
      double precision xi, yi, zi, chi
      double precision xj, yj, zj
      double precision x, y, z
      double precision dx, dy, dz
      double precision rsqr, rr, rr3, ch
      double precision etar, expar2, qt, erf, erfc, dfe
      integer nn1, nn2, nn3

      parameter (qp=0.3275911d0)
      parameter (a1=0.254829592d0)
      parameter (a2=-0.284496736d0)
      parameter (a3=1.421413741d00)
      parameter (a4=-1.453152027d0)
      parameter (a5=1.061405429d00)
!     parameter (eta=0.22d0)

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

!     eta = betapme/6.0d0
!     units of eclmb, edir and eself are eV
      call se_debugmsg0('in pme_direct2')
      eclmb = 0.0
      edir = 0.0
      eself = 0.0
      do i=1,3*natoms
        dedir(i)=0.d0
      enddo

!     unit of gradient (dedir) is kcal/A

!     tokcal = 332.0704
!     tokcal = 332.0806806d0
!     bpisqrt = betapme/pisqrt
      eta=betapme                        
      do iat=1,natoms
        iat3 = 3*iat-2
        xi = xyz(1,iat)
        yi = xyz(2,iat)
        zi = xyz(3,iat)
        chi = chgpme(iat)
    

!     . selfenergy (is a constant, it doesn't contribute to the gradient)
!       eself = eself - bpisqrt*chi*chi

        do jat=1,natoms
       
          jat3 = 3*jat-2
          if (pbc) then
            call se_pbcxyz(iat,jat,xj,yj,zj)
          else
            xj = xyz(1,jat)
            yj = xyz(2,jat)
            zj = xyz(3,jat)
          endif
!         do nn1=-1,1
!           do nn2=-1,1
!             do nn3=-1,1
          nn1=0
          nn2=0
          nn3=0
                if ((iat.eq.jat).and.
     .             ((abs(nn1)+abs(nn2)+abs(nn3)).eq.0)) then
                  eself=eself - 2.0d0*chi*chi*eta/dsqrt(pi)
!                 goto 128
!                 edir = 0.0
                else
                  x = xi-xj - nn1*dbox(1)
                  y = yi-yj - nn2*dbox(2)
                  z = zi-zj - nn3*dbox(3)
                  rsqr = x*x + y*y + z*z
                  rr = dsqrt(rsqr)
                  rr3 = rr**3
                  ch = chi*chgpme(jat)
!     .           short range Coulomb energy
!                 eclmb = eclmb + ch/rr

!     .           calculate contribution to gradient of short range Coulomb
                  dx = (ch*x)/rr3
                  dy = (ch*y)/rr3
                  dz = (ch*z)/rr3
                  etar = eta*rr
                  expar2 = exp(-etar*etar)
                  qt     = 1.0/(1.0+qp*etar)
                  erfc = ((((a5*qt+a4)*qt+a3)*qt+a2)*qt+a1)*qt*expar2
                  erf = 1.00D0 - erfc
                 
!     .           direct energy
                  edir = edir + ch*erfc/rr


!     .           calculate contribution to gradient of direct energy;
!     .           dfe is gradient/|r| of direct energy
                  dfe = ch*(-twopisqrt*etar*expar2+erf-1)/(rr*rsqr)
                 
                 
!     .           add direct and subtract coulomb contributions to gradient:
!     .           dEclmb/dr_i = SUM(j) -q_i*q_j*r_ij/|r_ij|^3
!     .           but here dEclmb/dr_i is used as a correction,
!     .           so change the sign.
                 
                  dedir(iat3+0) = dedir(iat3+0) + x*dfe
                  dedir(iat3+1) = dedir(iat3+1) + y*dfe
                  dedir(iat3+2) = dedir(iat3+2) + z*dfe
                 
                  dedir(jat3+0) = dedir(jat3+0) - x*dfe
                  dedir(jat3+1) = dedir(jat3+1) - y*dfe
                  dedir(jat3+2) = dedir(jat3+2) - z*dfe

                  if ((abs(nn1)+abs(nn2)+abs(nn3)).eq.0) then
                    dedir(iat3+0) = dedir(iat3+0) + dx
                    dedir(iat3+1) = dedir(iat3+1) + dy
                    dedir(iat3+2) = dedir(iat3+2) + dz
                  
                    dedir(jat3+0) = dedir(jat3+0) - dx
                    dedir(jat3+1) = dedir(jat3+1) - dy
                    dedir(jat3+2) = dedir(jat3+2) - dz    
                  
!                   edir = edir - ch/rr
                    eclmb= eclmb + ch/rr
                  endif
                endif
!             enddo
!           enddo
!         enddo
        enddo
      enddo
!     dedir(4) = - x*dfe - dx
!     dedir(5) = - y*dfe - dy
!     dedir(6) = - z*dfe - dz

!     energy conversion: until now energies are in e*e/A
!     multiply this by 1/(4 pi eps0) and then convert to V

      edir = elec2eV*edir*0.5d0
      eself = elec2eV*eself*0.5d0
      eclmb = elec2eV*eclmb*0.5d0
      do i=1,3*natoms
        dedir(i) = elec2eV*dedir(i)*eV2kcal*0.5d0
      enddo
!     edir = 14.39964521*edir
! 128 eself = 14.39964521*eself
!     eclmb = 14.39964521*eclmb
!     do i=1,3*natoms
!       dedir(i)=14.39964521* dedir(i)
!     enddo
      call se_debugmsg0('out pme_direct2')
      return
      end


