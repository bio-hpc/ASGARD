
!     new code: not parallelized!!!!!

      subroutine se_doferm(ierror)
!
!     uses an iterative bisection technique to determine the fermi
!     energy efermi.  error flag is set if more than 100 iterations
!     are necessary.
!
      implicit none
#include "sebomd.dim"
#include "se_electr.h"
#include "se_electr1.h"
#include "se_subsys.h"

      integer ierror
      double precision betah, etoler, emin, emax
      integer maxit, iorbmx
      integer imin, imax, niter, i, k
      double precision elecs, etotal, deltae, arg, diff

!
!     boltzmann constant in ev/kelvin:
!
      double precision boltz
      data boltz /8.617335408d-5/
      save boltz
!
      betah = 1.0d0/(boltz*tempk)
      etoler = 1.0d-8
      maxit = 100

      iorbmx = iorbpt(nsub+1)-1
!
!     determine the observed range of eigenvalues for all subsystems.
!
      emin = eval(1)
      emax = eval(iorbpt(2)-1)
      if(nsub.gt.1)then
        do k=1,nsub
          imin = iorbpt(k)
          imax = iorbpt(k+1)-1
          emin = min(emin,eval(imin))
          emax = max(emax,eval(imax))
        enddo
      endif

!
!     use a bisection technique to determine the fermi energy.
!
      elecs = nelec
      niter = 0
 100  efermi(1) = (emax + emin)/2.0d0
!
!     set fermi occupation numbers and use the sum of weighted squared
!     eigenvector coefficients evecsq to get the number of electrons.
!
      niter = niter + 1
      etotal = 0.0d0
      do i=1,iorbmx
        deltae = eval(i) - efermi(1)
        arg = betah*deltae
        if(arg.gt.50.0d0)then
          fermi(i) = 0.0d0
        elseif(arg.lt.-50.0d0)then
          fermi(i) = 1.0d0
        else
          fermi(i) = 1.0d0/(1.0d0 + exp(arg))
        endif
        etotal = etotal + fermi(i)*evecsq(i)
      enddo
!
!     we're done if the computed number of electrons is close enough
!     to the target.
!
      diff = abs(etotal - elecs)

      if(diff.lt.etoler) return
!
!     throw away lower or upper half of interval, depending upon whether
!     the computed number of electrons is too low or too high.  note that
!     raising the fermi energy increases the number of electrons because
!     more mo's become fully occupied.
!

      if(etotal.lt.elecs)then
        emin = efermi(1)
      else
        emax = efermi(1)
      endif
      if(niter.eq.maxit)then
        ierror = 1
        return
      endif
      go to 100

      end
