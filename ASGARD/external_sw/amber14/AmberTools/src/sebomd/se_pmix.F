      subroutine se_pmix(iter)
!
!     uses the improved iteration scheme (iis) to accelerate convergence
!     on the density matrix.  for details see:  p. badziag and f. solms,
!     comput. chem., vol 12, pp 233-236 (1988).  the iteration counter
!     corresponds to the scf cycle, i.e., iter=1 for the first iteration,
!     iter=2 for the second, etc.  this differs from badziag and solms
!     convention where iter=0 corresponds to the first cycle.
!
      implicit none
#include "sebomd.dim"
#include "se_bpairs.h"
#include "se_convrg.h"
#include "se_elemnt.h"
#include "se_fulsys.h"
#include "se_glbmat1.h"
#include "se_glbmat2.h"

      integer iter

      integer iimax, ijmax, ii, ij, i, norbsi, ip, jp
      double precision pold, pnew, dd12, dd11, dii, d11tmp, d12tmp
      double precision denom
      double precision a1, a2, xlold, xlnew, xltest

      logical mixed
      save a1,a2,xlold,xlnew,mixed

      if(iter.eq.1)then
        a1 = 1.0d0
        a2 = 1.0d0
        xlnew = 1.0d0
        mixed = .false.
      endif
      iimax = iimat(natoms+1)-1
      ijmax = ijmat(ip1(natoms+1))-1
      if(iter.eq.2)then
        do ii=1,iimax
          dpii(ii) = pdiag(ii) - piiold(ii)
        end do
      endif
      if(iter.lt.3) return
!
!     do the first mixing step if we've completed two iterations.
!
      if(iter.ge.3)then
        pold = 1.0d0 - a2
        pnew = a2
        do ii=1,iimax
          pdiag(ii) = pnew*pdiag(ii) + pold*piiold(ii)
        end do
        do ij=1,ijmax
          pdiat(ij) = pnew*pdiat(ij) + pold*pijold(ij)
        end do
      endif
!
!     get new lambda parameter (xlnew) and form new difference matrix.
!
      xlold = xlnew
      dd12 = 0.0d0
      dd11 = 0.0d0
      do i=1,natoms
        norbsi = natorb(iatnum(i))
        if(norbsi.ne.0) then
          ii = iimat(i)
          do ip=1,norbsi
            do jp=1,ip
              dii = dpii(ii)
              d11tmp = dii*dii
              dd11 = dd11 + d11tmp
              dpii(ii) = pdiag(ii) - piiold(ii)
              d12tmp = dii*dpii(ii)
              dd12 = dd12 + d12tmp
              if(ip.eq.jp)then
                dd11 = dd11 - 0.5d0*d11tmp
                dd12 = dd12 - 0.5d0*d12tmp
              endif
              ii = ii + 1
            end do
          end do
        endif
      end do
      xlnew = dd12/dd11
      xltest = abs(xlnew-xlold)/max(abs(xlold),abs(xlnew))
      if(.not.mixed)then
        xltest = abs(xlnew-xlold)/max(abs(xlold),abs(xlnew))
        if(xltest.le.0.1d0)then
          mixed = .true.
        else
          return
        endif
      endif
      denom = 1.0d0 - xlnew
      if(abs(denom).lt.0.01d0) denom = sign(0.01d0,denom)
      a1 = 1.0d0/denom
      if((a1.gt.0.0d0.and.a1.le.1.0d0))then
        pnew = a1
        pold = 1.0d0 - pnew
        do ii=1,iimax
          pdiag(ii) = pnew*pdiag(ii) + pold*piiold(ii)
          dpii(ii) = pnew*dpii(ii)
        end do
        do ij=1,ijmax
          pdiat(ij) = pnew*pdiat(ij) + pold*pijold(ij)
        end do
        xlnew = pnew*xlnew
      endif
      a2 = a1*a2
      if(a2.le.0.0d0.or.a2.gt.1.0d0)then
        denom = 1.0d0 + a2
        if(abs(denom).lt.0.01d0) denom = sign(0.01d0,denom)
        a2 = 2.0d0*a2/denom
      endif
      if(abs(a2).gt.1.5d0) a2 = sign(1.5d0,a2)
      return
      end
