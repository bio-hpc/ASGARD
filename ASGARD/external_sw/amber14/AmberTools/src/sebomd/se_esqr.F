      SUBROUTINE se_ESQR(K,NORBSK,EVEC1)

!     determines the coefficients evecsq for subsystem k:

!     evecsq(l0+l) = d11*evec1(1,l)**2 + d22*evec1(2,l)**2 + . . .

!     here l0 is just a global pointer for subsystem k, and dii is
!     the normalization factor that accounts for subsystem overlap.

!     once all of the entries of evecsq have been determined, then
!     the number of electrons can be computed for any set of fermi
!     occupation numbers:

!     no. of electrons = fermi(1)*evecsq(1) + fermi(2)*evecsq(2) + . . .

      implicit none
#include "sebomd.dim"
#include "se_electr.h"
#include "se_elemnt.h"
#include "se_fulsys.h"
#include "se_subsys.h"
      double precision evec1
      dimension evec1(norbsk,*)
      integer k, norbsk

      double precision suml, dii, sumi
      integer l0, l, i, iatm, norbsi, nsubsi, i0, iorb

      l0 = iorbpt(k)-1
      do l=1,norbsk
        suml = 0.0d0
        do i=iatom1(k),iatom1(k+1)-1
          iatm = iatoms(i)
!         skip contribution if this is a buffer atom.
          if(iabuff(i).eq.0) then
            norbsi = natorb(iatnum(iatm))
            nsubsi = isub1(iatm+1) - isub1(iatm)
            dii = 2.0d0/nsubsi
            i0 = iorb1(i)-1
            sumi = 0.0d0
            do iorb=1,norbsi
              sumi = sumi + evec1(i0+iorb,l)**2
            end do
            suml = suml + dii*sumi
          endif
        end do
        evecsq(l0+l) = suml
      end do
      return
      end
