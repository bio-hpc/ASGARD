      subroutine se_atmchg

      implicit none
#include "sebomd.dim"
#include "se_xyzcrd.h"
#include "se_fulsys.h"
#include "se_bpairs.h"
#include "se_logic.h"
#include "se_glbmat1.h"
#include "se_glbmat2.h"
#include "se_elemnt.h"

      double precision sp2(maxatm)

      integer ii, i, iai, norbsi, iatk, iat
      integer j, norbsj, jatk, jat
      integer norbsij, ijstart, ijend, ij, npairs, iiorb, ijaddr
      double precision zi, psumii
      double precision bk, bnc1, bnc2, bos, bno, bhk
      double precision p2
      logical calccm1, calccm2

      calccm1 = .true.
      calccm2 = .true.

!     get the Mulliken charges
      ii = 0
      do i=1,natoms
        iai = iatnum(i)
        if(iai.eq.0) cycle
        norbsi = natorb(iai)
        zi = zchg(iai)
        psumii = 0.0d0
        if(norbsi.gt.0)then
          do iiorb=1,norbsi
            ii = ii + iiorb
            psumii = psumii + pdiag(ii)
          enddo
        endif
        atchg(i) = zi - psumii
        atchg2(i) = atchg(i)
        atchg3(i) = atchg(i)
      enddo

!     for now atchg refers to Mulliken charges, atchg2 to CM1 charges,
!     atchg3 to CM2 charges (swap this if needed at the end of the
!     routine)

      if (calccm1.or.calccm2) then

!     .  there are two efficient ways to calculate the CM1 charges:
!     .  the faster one is a loop from iat=2,natoms with
!     .  an inner loop from jat=1,iat-1, calculate and store the
!     .  B(iat,jat) = B(jat, iat) coefficients and then
!     .  do another iat=2,natoms - jat=1,iat-1 loop for
!     .  the normalization using the stored B(iat,jat) = B(jat,iat)
!     .  coefficients.
!     .  The second slower one loops from iat=1,natoms in which the
!     .  inner loop jat=1,natoms calculates the B(iat,jat) coefficients.
!     .  After this jat-loop a partial normalization is performed.
!     .  Method one is faster because it only performs the
!     .  density matrix summation needed for the B(iat,jat) calculation
!     .  N^2/2 times (method 2 does this N^2 times),
!     .  but requires a lot more memory: N^2/2 elements to
!     .  store the B(iat,jat) coefficients. Method 2 only needs to
!     .  store the N intermediate B(iat,jat) coefficients within the
!     .  iat-loop (N is the number of atoms).
!     .  Here memory considerations lead me to implement the second
!     .  (slower) method.
!     .  CM2 charges can be cheaply obtained in the same sweep.

         if (am1) then

!     .     calculate the CM1A / CM2A charges

           do iat=1,natoms
             iatk = iatnum(iat)
             if (iatk.eq.0) cycle

             bk = 0.0
             bnc1 = 0.0
             bnc2 = 0.0
             bos = 0.0
             bno = 0.0
             bhk = 0.0

             norbsi = natorb(iatk)
             npairs = ip1(natoms+1)-1
             do jat=1,natoms
               sp2(jat) = 0.0
               if (jat.eq.iat) cycle
               jatk = iatnum(jat)
               if (jatk.ne.0) then
                 i = max(iat,jat)
                 j = min(iat,jat)
                 call se_ijfind(npairs,i,j,ijaddr)
                 if (ijaddr.ne.0) then
                   norbsj = natorb(jatk)
                   norbsij = norbsi*norbsj - 1
                   ijstart = ijmat(ijaddr)
                   ijend = ijstart + norbsij
                   do ij=ijstart,ijend
                     p2 = pdiat(ij)*pdiat(ij)
                     sp2(jat) = sp2(jat) + p2
                     bk = bk + p2
                   enddo
                   if (calccm1)
     &               call se_getcm1a(iatk,jat,jatk,bhk,
     &                 bnc1,bnc2,bno,bos,sp2)
                   if (calccm2)
     &               call se_getcm2a(iat,iatk,jat,jatk,sp2)
                 endif
               endif
             enddo
             if (calccm1)
     &         call se_calccm1a(iat,iatk,bhk,bnc1,bnc2,
     &           bno,bos,bk,sp2)
           enddo
         elseif (pm3) then
!     .    calculate CM1P / CM2P charges
           do iat=1,natoms
             iatk = iatnum(iat)
             if (iatk.eq.0) cycle
             bk = 0.0
             bhk = 0.0
             norbsi = natorb(iatk)
             npairs = ip1(natoms+1)-1
             do jat=1,natoms
               sp2(jat) = 0.0
               if (jat.eq.iat) cycle
               jatk = iatnum(jat)
               if (jatk.ne.0) then
                 i = max(iat,jat)
                 j = min(iat,jat)
                 call se_ijfind(npairs,i,j,ijaddr)
                 if (ijaddr.ne.0) then
                   norbsj = natorb(jatk)
                   norbsij = norbsi*norbsj - 1
                   ijstart = ijmat(ijaddr)
                   ijend = ijstart + norbsij
                   do ij=ijstart,ijend
                     p2 = pdiat(ij)*pdiat(ij)
                     sp2(jat) = sp2(jat) + p2
                     bk = bk + p2
                   enddo
                   if (calccm1)
     &               call se_getcm1p(iatk,jat,jatk,bhk,sp2)
                   if (calccm2)
     &               call se_getcm2p(iat,iatk,jat,jatk,sp2)
                 endif
               endif
             enddo
             if (calccm1)
     &         call se_calccm1p(iat,iatk,bhk,bk,sp2)
           enddo
         endif
      endif
      end
!
!     subroutines for CM1/AM1 charge model
!
      subroutine se_getcm1a(iatk,jat,jatk,bhk,bnc1,bnc2,bno,bos,sp2)

      implicit none
#include "sebomd.dim"

      integer iatk, jat, jatk
      double precision bhk, bnc1, bnc2, bno, bos, sp2(maxatm)

      double precision d, t

      if (iatk.eq.1) then
!     .  H atom
         if (jatk.eq.7) then
!     .     overlap with N atom
            d = 0.0850
         elseif (jatk.eq.8) then
!     .     overlap with O atom
            d = 0.1447
         elseif (jatk.eq.14) then
!     .     overlap with Si atom
            d = 0.0640
         else
            d = 0.0
         endif
         bhk = bhk + d*sp2(jat)
      elseif (iatk.eq.7) then
!     .   N atom
         if (jatk.eq.6) then
!     .     overlap with C atom
            t = 10.0*(sp2(jat)-2.3)
            t = 0.5*tanh(t)+0.5
            bnc1 = bnc1 + t
            bnc2 = bnc2 - 0.0880*t
         elseif (jatk.eq.8) then
!     .     overlap with O atom
            bno = bno - 0.0630*sp2(jat)
         endif
      elseif ((iatk.eq.8).and.(jatk.eq.16)) then
!     .  overlap S atom with O atom
         bos = bos - 0.060*sp2(jat)
      endif
      end
!--------------------------------------------------------------------------------
      subroutine se_calccm1a(iat,iatk,bhk,bnc1,bnc2,bno,bos,bk,sp2)
      implicit none
#include "sebomd.dim"
#include "se_fulsys.h"
#include "se_xyzcrd.h"

      integer iat, iatk
      double precision bhk, bnc1, bnc2, bno, bos, bk, sp2(maxatm)

      integer jat
      double precision c, d, dq

      bnc1 = -0.3846*bnc1

      if (iatk.eq.7) then
!     .  N atom
         c = 0.3846
         d = 0.0
      elseif (iatk.eq.8) then
!     .  O atom
         c = 0.0
         d = -0.0283
      elseif (iatk.eq.9) then
!     .  F atom
         c = 0.1468
         d = 0.0399
      elseif (iatk.eq.16) then
!     .  S atom
         c = -0.1311
         d = -0.0956
      elseif (iatk.eq.17) then
!     .  Cl atom
         c = 0.0405
         d = -0.0276
      elseif (iatk.eq.35) then
!     .  Br atom
         c = 0.1761
         d = -0.0802
      elseif (iatk.eq.53) then
!     .  I atom
         c = 0.2380
         d = -0.1819
      else
         c = 0.0
         d = 0.0
      endif

      dq = (c + bnc1)*atchg(iat) + d + bhk + bos +
     &     bnc2 + bno
      atchg2(iat) = atchg2(iat) + dq*bk
!     partial normalization:
      do jat=1,natoms
         atchg2(jat) = atchg2(jat) - sp2(jat)*dq
      enddo

      end
!
!     subroutines for CM1/PM3 charge model
!
      subroutine se_getcm1p(iatk,jat,jatk,bhk,sp2)

      implicit none

#include "sebomd.dim"
      integer iatk, jat, jatk
      double precision bhk, sp2(maxatm)


      double precision d

      if (iatk.eq.1) then
!     .  H atom
         if (jatk.eq.7) then
!     .     overlap with N atom
            d = 0.1854
         elseif (jatk.eq.8) then
!     .     overlap with O atom
            d = 0.1434
         elseif (jatk.eq.14) then
!     .     overlap with Si atom
            d = -0.1004
         else
            d = 0.0
         endif
         bhk = bhk + d*sp2(jat)
      endif

      end
!--------------------------------------------------------------------
      subroutine se_calccm1p(iat,iatk,bhk,bk,sp2)

      implicit none
#include "sebomd.dim"
#include "se_fulsys.h"
#include "se_xyzcrd.h"

      integer iat, iatk
      double precision bhk, bk, sp2(maxatm)

      integer jat
      double precision c, d, dq

      if (iatk.eq.7) then
!     .  N atom
         c = 0.0
         d = -0.0909
      elseif (iatk.eq.8) then
!     .  O atom
         c = 0.0
         d = -0.0449
      elseif (iatk.eq.9) then
!     .  F atom
         c = 0.3381
         d = 0.0148
      elseif (iatk.eq.16) then
!     .  S atom
         c = -0.0834
         d = -0.0848
      elseif (iatk.eq.17) then
!     .  Cl atom
         c = -0.1080
         d = -0.1168
      elseif (iatk.eq.35) then
!     .  Br atom
         c = -0.0116
         d = -0.0338
      elseif (iatk.eq.53) then
!     .  I atom
         c = -0.3213
         d = -0.0636
      else
         c = 0.0
         d = 0.0
      endif

      dq = c*atchg(iat) + d + bhk
      atchg2(iat) = atchg2(iat) + dq*bk
!     partial normalization:
      do jat=1,natoms
         atchg2(jat) = atchg2(jat) - sp2(jat)*dq
      enddo

      end
!
!     subroutines for CM2/AM1 charge model
!
      subroutine se_getcm2a(iat,iatk,jat,jatk,sp2)

      implicit none
#include "sebomd.dim"
#include "se_xyzcrd.h"
      integer iat, iatk, jat, jatk
      double precision sp2(maxatm)

      integer i, j
      double precision c, d, t

      i = min(iatk,jatk)
      j = max(iatk,jatk)

      c = 0.0
      d = 0.0
      if (i.eq.1) then
!     .  H atom
         if (j.eq.6) then
!     .     overlap with C atom
            c = -0.020
         elseif (j.eq.7) then
!     .     overlap with N atom
            c = 0.207
         elseif (j.eq.8) then
!     .     overlap with O atom
            c = 0.177
         elseif (j.eq.14) then
!     .     overlap with Si atom
            c = -0.083
         elseif (j.eq.15) then
!     .     overlap with P atom
            d = 0.103
         elseif (j.eq.16) then
!     .     overlap with S atom
            c = 0.038
         endif
      elseif (i.eq.6) then
!     .  C atom
         if (j.eq.7) then
!     .     overlap with N atom
            c = 0.008
            d = 0.086
         elseif (j.eq.8) then
!     .     overlap with O atom
            c = 0.026
            d = 0.016
         elseif (j.eq.9) then
!     .     overlap with F atom
            d = 0.019
         elseif (j.eq.14) then
!     .     overlap with Si atom
            c = 0.062
         elseif (j.eq.15) then
!     .     overlap with P atom
            d = -0.019
         elseif (j.eq.16) then
!     .     overlap with S atom
            c = -0.059
            d = 0.171
         elseif (j.eq.17) then
!     .     overlap with Cl atom
            d = 0.027
         elseif (j.eq.35) then
!     .     overlap with Br atom
            d = 0.081
         elseif (j.eq.53) then
!     .     overlap with I atom
            d = 0.147
         endif
      elseif (i.eq.7) then
!     .  N atom
         if (j.eq.8) then
!     .     overlap with O atom
            c = -0.197
            d = 0.134
         endif
      elseif (i.eq.8) then
!     .  O atom
         if (j.eq.15) then
!     .     overlap with P atom
            d = 0.088
         endif
      elseif (i.eq.9) then
!     .  F atom
         if (j.eq.15) then
!     .     overlap with P atom
            d = 0.252
         endif
      elseif (i.eq.15) then
!     .  P atom
         if (j.eq.16) then
!     .     overlap with S atom
            d = -0.080
         endif
      endif

      t = sp2(jat)*(d + c*sp2(jat))
      if (i.ne.iatk) t = -t
      atchg3(iat) = atchg3(iat) + t

      end
!--------------------------------------------------------------------------------
      subroutine se_getcm2p(iat,iatk,jat,jatk,sp2)

      implicit none
#include "sebomd.dim"
#include "se_xyzcrd.h"

      integer iat, iatk, jat, jatk
      double precision sp2(maxatm)

      integer i, j
      double precision c, d, t

      i = min(iatk,jatk)
      j = max(iatk,jatk)

      c = 0.0
      d = 0.0
      if (i.eq.1) then
!     .  H atom
         if (j.eq.6) then
!     .     overlap with C atom
            c = 0.003
         elseif (j.eq.7) then
!     .     overlap with N atom
            c = 0.274
         elseif (j.eq.8) then
!     .     overlap with O atom
            c = 0.185
         elseif (j.eq.14) then
!     .     overlap with Si atom
            c = -0.021
         elseif (j.eq.15) then
!     .     overlap with P atom
            d = 0.253
         elseif (j.eq.16) then
!     .     overlap with S atom
            c = 0.089
         endif
      elseif (i.eq.6) then
!     .  C atom
         if (j.eq.7) then
!     .     overlap with N atom
            c = -0.022
            d = 0.156
         elseif (j.eq.8) then
!     .     overlap with O atom
            c = 0.025
            d = 0.016
         elseif (j.eq.9) then
!     .     overlap with F atom
            d = 0.025
         elseif (j.eq.14) then
!     .     overlap with Si atom
            c = -0.107
         elseif (j.eq.15) then
!     .     overlap with P atom
            d = 0.082
         elseif (j.eq.16) then
!     .     overlap with S atom
            c = -0.033
            d = 0.112
         elseif (j.eq.17) then
!     .     overlap with Cl atom
            d = 0.117
         elseif (j.eq.35) then
!     .     overlap with Br atom
            d = 0.040
         elseif (j.eq.53) then
!     .     overlap with I atom
            d = -0.032
         endif
      elseif (i.eq.7) then
!     .  N atom
         if (j.eq.8) then
!     .     overlap with O atom
            c = -0.030
            d = -0.043
         endif
      elseif (i.eq.8) then
!     .  O atom
         if (j.eq.15) then
!     .     overlap with P atom
            d = 0.181
         elseif (j.eq.16) then
!     .     overlap with S atom
            d = 0.056
         endif
      elseif (i.eq.9) then
!     .  F atom
         if (j.eq.15) then
!     .     overlap with P atom
            d = 0.244
         endif
      elseif (i.eq.15) then
!     .  P atom
         if (j.eq.16) then
!     .     overlap with S atom
            d = -0.087
         endif
      endif

      t = sp2(jat)*(d + c*sp2(jat))
      if (i.ne.iatk) t = -t
      atchg3(iat) = atchg3(iat) + t

      end
