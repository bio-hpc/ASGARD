      subroutine se_pme_setup
!
! setup some variables for Constant Charge PME
!

      implicit none
#include "sebomd.dim"
#include "se_pbcbox.h"
#include "se_pmeee.h"
      integer i
      double precision cutoff
      double precision skinnb
      double precision dsum_tol
      double precision rsum, rsum_tol
      double precision uk1, uk2, uk3, uksqr, expo, fact

!     define betapme

      cutoff=0.5d0*min(dbox(1), dbox(2), dbox(3))
      skinnb = 2.0d0
      dsum_tol = 1.0d-5
      call find_ewaldcof(cutoff-skinnb,dsum_tol,betapme)
      betapme2=betapme*betapme

      call se_debugmsg1f("pme_setup - betapme=",betapme)

!     define k1pme, k2pme, k3pme
      k1pme = 2
      do while (k1pme < dbox(1))
        k1pme = k1pme*2
      end do
      k2pme = 2
      do while (k2pme < dbox(2))
        k2pme = k2pme*2
      end do
      k3pme = 2
      do while (k3pme < dbox(3))
        k3pme = k3pme*2
      end do

      call se_debugmsg3i("pme_setup - k[123]pme = ", k1pme,k2pme,k3pme)
      ! check values
      if (k1pme > maxkpme) then
        write(6,'("ERROR: in pme_setup: k1pme (",i4,
     .            ") is larger than MAXKPME (",i4")")') k1pme, maxkpme
        write(6,'("Increase MAXKPME value in sebomd.dim")')
        call mexit(6,1)
      endif
      if (k2pme > maxkpme) then
        write(6,'("ERROR: in pme_setup: k2pme (",i4,
     .            ") is larger than MAXKPME (",i4")")') k2pme, maxkpme
        write(6,'("Increase MAXKPME value in sebomd.dim")')
        call mexit(6,1)
      endif
      if (k3pme > maxkpme) then
        write(6,'("ERROR: in pme_setup: k3pme (",i4,
     .            ") is larger than MAXKPME (",i4")")') k3pme, maxkpme
        write(6,'("Increase MAXKPME value in sebomd.dim")')
        call mexit(6,1)
      endif

!     define nspline

      nspline = max(4,max(k1pme,k2pme,k3pme)/2)

      call se_debugmsg1i("pme_setup - nspline = ", nspline)
      if (nspline > maxspline) then
        write(6,'("ERROR: in pme_setup: nspline (",i4,
     .       ") is greater than MAXSPLINE (",i4")")') nspline, maxspline
        write(6,'("Increase MAXSPLINE value in sebomd.dim")')
        call mexit(6,1)
      endif

!     search for mmax
      rsum_tol = 1.d-8
      rsum = 1.0d0
      mmax = 1
      fact=4.0d0*pi/boxvol
      do while (abs(rsum)>rsum_tol)
        uk1=mmax*2.0d0*pi/dbox(1)
        uk2=mmax*2.0d0*pi/dbox(2)
        uk3=mmax*2.0d0*pi/dbox(3)
        uksqr=uk1*uk1+uk2*uk2+uk3*uk3
        expo=dexp(-uksqr/(4.0d0*betapme**2))
        rsum = fact*expo/uksqr
        mmax = mmax+1
      end do
      call se_debugmsg1i("pme_setup: mmax = ",mmax)

      k123pme = k1pme*k2pme*k3pme
      k123pme2 = 2*k123pme

      dk1pme = dble(k1pme)
      dk2pme = dble(k2pme)
      dk3pme = dble(k3pme)

      k1pmem1 = k1pme-1
      k2pmem1 = k2pme-1
      k3pmem1 = k3pme-1

      k1pmek2pme = k1pme*k2pme
      k3pmenspl = k3pme-nspline+1
      k2pmenspl = k2pme-nspline+1
      k1pmenspl = k1pme-nspline+1

      dnspline = dble(nspline)

!     we always have a rectangular box
      do i=1,3
         recip1(i) = 0.0
         recip2(i) = 0.0
         recip3(i) = 0.0
      enddo
      recip1(1) = 1.0/dbox(1)
      recip2(2) = 1.0/dbox(2)
      recip3(3) = 1.0/dbox(3)

      end
