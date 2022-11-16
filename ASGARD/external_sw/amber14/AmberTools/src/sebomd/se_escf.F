      subroutine se_escf(eelect1)

!     determines the electronic energy for the current scf iteration.

      implicit none
#include "sebomd.dim"
#include "se_bpairs.h"
#include "se_elemnt.h"
#include "se_esave.h"
#include "se_fulsys.h"
#include "se_glbmat1.h"
#include "se_glbmat2.h"
#include "se_logic.h"
#include "se_parproc.h"
#include "se_pmeqm.h"
#include "se_xyzcrd.h"

      double precision eelect1, eii, epmeqm

      integer iai, norbsi, iorb, ijphi
      integer iimax, ijmax, i, j, ii, ij


      eelect1 = 0.0d0
      iimax = iimat(natoms+1)-1
      ijmax = ijmat(ip1(natoms+1))-1
      do i=1,iimax
        eelect1 = eelect1 + (hdiag(i) + fdiag(i))*pdiag(i)
      end do
      eii = 0.0d0
      ii = 0
      do i=1,natoms
        iai = iatnum(i)
        norbsi = natorb(iai)

!       skip atom if it's a dummy or a sparkle.
        if(norbsi.ne.0) then
          do iorb=1,norbsi
            ii = ii + iorb
            eii = eii + (hdiag(ii) + fdiag(ii))*pdiag(ii)
          end do
        endif
      end do

      eelect1 = eelect1 - 0.5d0*eii

      do ij=1,ijmax
        eelect1 = eelect1 + (hdiat(ij) + fdiat(ij))*pdiat(ij)
      end do

      if (pmeqm) then
        ijphi=0
        epmeqm=0.0d0

        do i=1,natoms
        iai = iatnum(i)
          do j=1,natoms
            ijphi=ijphi+1
            if (chewald) then
              epmeqm = epmeqm + phipme(ijphi)*chgpme(j)*chgpme(i)
            elseif (mullewald) then
              epmeqm = epmeqm + phipme(ijphi)*atchg(j)*atchg(i)
            elseif (cmewald) then
              epmeqm = epmeqm + phipme(ijphi)*atchg2(j)*atchg2(i)
            else 
              write(iout,*)'definition of the set of charges (escf1)'
              stop
            endif
!           epmeqm = epmeqm + phipme(ijphi)*atchg(j)*zchg(iai)
          end do
        enddo
        
        eelect1 = eelect1 + 0.5d0*epmeqm
        virlrqm = 0.5d0*epmeqm
        ijphi=0
        do i=1,natoms
           do j=1,natoms
              ijphi=ijphi+1
              if (chewald) then
                virlrqm = virlrqm-0.5d0*phisr(ijphi)*chgpme(j)*
     &          chgpme(i)
              elseif (mullewald) then
                virlrqm = virlrqm-0.5d0*phisr(ijphi)*atchg(j)*atchg(i)
              elseif (cmewald) then
                virlrqm = virlrqm - 0.5d0*phisr(ijphi)*atchg2(j)*
     &          atchg2(i)
              else 
                write(iout,*)'definition of the set of charges (escf2)'
                stop
              endif
!           epmeqm = epmeqm + phipme(ijphi)*atchg(j)*zchg(iai)
          end do
        enddo
      endif
      return
      end
