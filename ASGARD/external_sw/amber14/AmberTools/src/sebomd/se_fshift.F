      subroutine se_fshift(idiag,eshift)

!     routine to carry out dynamic or static level shift to prevent
!     auto-polarization of charge distribution.

      implicit none
#include "sebomd.dim"
#include "se_bpairs.h"
#include "se_elemnt.h"
#include "se_fulsys.h"
#include "se_glbmat1.h"
#include "se_glbmat2.h"
#include "se_logic.h"
      integer idiag
      double precision eshift

c     local:
      double precision oldshf, eshft2
      integer iimax, ijmax, ii, i, norbsi, iorb, ij

      logical first
      data first /.true./
      save first
      save oldshf

      if(first)then
        first = .false.
      endif

c     no shifting if density matrix has diagonal form.

      if(idiag.eq.1)then
        oldshf = eshift
        return
      endif

      if(stand)then
        eshift = 0.5d0*(eshift + oldshf)
        oldshf = eshift
      endif
      eshft2 = eshift*2.0d0
      iimax = iimat(natoms+1)-1
      ijmax = ijmat(ip1(natoms+1))-1
      do ii=1,iimax
        fdiag(ii) = fdiag(ii) + eshift*pdiag(ii)
      end do
      ii = 0
      do i=1,natoms
        norbsi = natorb(iatnum(i))
        if(norbsi.ne.0) then
          do iorb=1,norbsi
            ii = ii + iorb
            fdiag(ii) = fdiag(ii) - eshft2
          end do
        end if
      end do
      do ij=1,ijmax
        fdiat(ij) = fdiat(ij) + eshift*pdiat(ij)
      end do
      return
      end
