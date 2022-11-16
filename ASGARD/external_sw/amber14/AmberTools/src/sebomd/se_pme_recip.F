!------------------------------------------------------------------------
      subroutine se_pme_recip(erec,derec)

      implicit none
#include "sebomd.dim"
#include "se_pmeee.h"
#include "se_constants.h"

!     this routine calculates the reciprocal enery and derivatives
!     by the PME method of J.Chem.Phys. 103 ('95) 8577.
!
!     written by Arjan van der Vaart, Dec. '97
!
!     parameters:
!     *) erec            the reciprocal energy
!     *) derec           the derivative of erec to the coordinates
!

!     parameters:
      double precision derec(maxpar), erec

!     local:
      integer ndim, ndimn(3)
      integer i, j, k, kk

      call se_debugmsg0('in pme_recip')
!     initialization for the Fourier transforms
      ndim = 3
      ndimn(1) = k1pme
      ndimn(2) = k2pme
      ndimn(3) = k3pme

!     calculate the function Q
      call se_pme_calcq

!     copy Q to qpme, since qpmec will be overwritten with F[Q]
      j = 1
      do i=1,k123pme
         qpmec(j) = qpme(i)
         j = j + 1
         qpmec(j) = 0.0
         j = j + 1
      end do

!     Fourier transform Q and multiply with theta
      call se_fourn(qpmec,ndimn,ndim,iminone)
      j = 1
      do i=1,k123pme
         qpmec(j) = qpmec(j)*thetapme(i)
         j = j + 1
         qpmec(j) = qpmec(j)*thetapme(i)
         j = j + 1
      end do

!     after the following transformation, qpmec will contain Conv(theta,Q)
      call se_fourn(qpmec,ndimn,ndim,ione)

!     calculate the gradient
      call se_pme_derec(derec)

!     calculate the reciprocal energy in kcal/mol
      erec = 0.0
      kk = 1
      do k=1,k123pme
         erec = erec + qpmec(kk)*qpme(k)
         kk = kk + 2
      end do
     
      erec = elec2ev*erec*0.5d0
!      do i=1,3*natoms
!         derec(i)=derec(i)*7.199822605
!      enddo
   
      call se_debugmsg0('out pme_recip')
      return
      end
