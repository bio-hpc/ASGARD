      subroutine se_pme_calccm1

      implicit none
#include "sebomd.dim"
#include "se_xyzcrd.h"
#include "se_glbmat1.h"
#include "se_glbmat2.h"
#include "se_bpairs.h"
#include "se_fulsys.h"
#include "se_elemnt.h"
#include "se_parproc.h"
#ifdef MPI
#include "mpif.h"
#endif

      double precision sp2(maxatm)
      double precision bksave(maxatm)
      double precision dqsave(maxatm)
      double precision cksave(maxatm)
 
      logical calccm1, calccm2

      integer iimax, ipmax, ijmax
      integer ii, ij
      integer i, iai, norbsi, iatk, iat, iorb, iiorb, iorbi
      integer norbsj, jatk, jat, iorbj
      integer npairs, norbsij, ijstart, ijend, ijdiag, ijbond, ij1
      integer ijaddr
      double precision zi, psumii, bk, bhk, ck, dq
      double precision p2, dhki, dhkj

      iimax = iimat(natoms+1)-1
      ipmax = ip1(natoms+1)
      ijmax = ijmat(ipmax)-1

      do i=1,mxpair
       bkk1(i)=0.0d0
      enddo
      do ii=1,iimax
        fpmediag(ii) = 0.0d0
      enddo
      do ij=1,ijmax
        fpmediat1(ij) = 0.0d0
        fpmediat2(ij) = 0.0d0
      enddo  
      calccm1 = .true.
      calccm2 = .false.


!     get the mulliken charges
      ii = 0
      do i=1,natoms
        iai = iatnum(i)

!       skip dummy atom.
        if(iai.ne.0) then
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
        endif
      enddo

!     for now atchg refers to Mulliken charges, atchg2 to CM1 charges,
!     atchg3 to CM2 charges (swap this if needed at the end of the
!     routine)

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



!     .     calculate CM1P / CM2P charges

      do iat=1,natoms
        iatk = iatnum(iat)
        if (iatk.ne.0) then
          bk = 0.0
          bhk = 0.0
         
          norbsi = natorb(iatk)
          npairs = ip1(natoms+1)-1
         
          do jat=1,natoms
            sp2(jat) = 0.0
            if (jat.ne.iat) then
              jatk = iatnum(jat)
              if (jatk.ne.0) then
!               i = max(iat,jat)
!               j = min(iat,jat)
                if (iat.gt.jat) then
                  call se_ijfind(npairs,iat,jat,ijaddr)
                else
                  call se_ijfind(npairs,jat,iat,ijaddr)
                endif
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
                  call se_pme_getcm1p(iatk,jat,jatk,bhk,sp2)
                endif
              endif
            endif
          enddo
         
          call se_pme_calccm1p(iat,iatk,bhk,bk,sp2,ck,dq)
         
          dqsave(iat)=dq
          bksave(iat)=bk
          cksave(iat)=ck
        endif
      enddo

#ifdef MPI 
      do  iat=myid+1,natoms,nproc
#else
      do  iat=1,natoms
#endif
        iatk = iatnum(iat)
        if (iatk.ne.0) then
          norbsi = natorb(iatk)
          ijdiag = iimat(iat)-1
          do iorb=1,norbsi
            ijdiag=ijdiag+iorb
            fpmediag(ijdiag)= -1 -bksave(iat)*cksave(iat)
          enddo  
        endif
      enddo
           
      ijbond=1
      do iat=2,natoms
        iatk = iatnum(iat)
        if (iatk.ne.0) then
          norbsi = natorb(iatk)
          do jat=1,iat-1
            jatk = iatnum(jat)
            if (jatk.ne.0) then
              call se_ijfind(npairs,iat,jat,ijaddr)
              if (ijaddr.ne.0) then
                norbsj = natorb(jatk)
                ij1= ijmat(ijaddr)
                ij=ijmat(ijaddr)
                call se_pme_getdhk(iatk,jatk,dhkj)
                call se_pme_getdhk(jatk,iatk,dhki)
                bkk1(ijaddr)=0.0d0
                do iorbi=1,norbsi  
                  do iorbj=1,norbsj  
                    bkk1(ijaddr) = bkk1(ijaddr) + pdiat(ij1)*pdiat(ij1)
                    ij1=ij1+1
                  enddo
                enddo
                do iorbi=1,norbsi  
                  do iorbj=1,norbsj  
                    fpmediat1(ij) = bksave(iat)*dhkj*pdiat(ij) 
     &                            - bkk1(ijaddr)*dhki*pdiat(ij)
     &                            - pdiat(ij)*dqsave(jat)
     &                            + pdiat(ij)*dqsave(iat)
                    fpmediat2(ij) = bksave(jat)*dhki*pdiat(ij)
     &                            - bkk1(ijaddr)*dhkj*pdiat(ij)
     &                            - pdiat(ij)*dqsave(iat)
     &                            + pdiat(ij)*dqsave(jat)
         
                    ij=ij+1
                  enddo
                enddo
              endif
            endif  
            ijbond=ijbond+1
          enddo
        endif
      enddo
      end
!--------------------------------------------------------------------------------
      subroutine se_pme_getcm1p(iatk,jat,jatk,bhk,sp2)
!
!     subroutines for CM1/PM3 charge model
!
      implicit none
#include "sebomd.dim"
      integer iatk, jat, jatk
      double precision bhk, sp2(maxatm)
      double precision d

      if (iatk.eq.1) then
!       H atom
        if (jatk.eq.7) then
!         overlap with N atom
          d = 0.1854d0
        elseif (jatk.eq.8) then
!         overlap with O atom
          d = 0.1434d0
        elseif (jatk.eq.14) then
!         overlap with Si atom
          d = -0.1004d0
        else
          d = 0.0d0
        endif
        bhk = bhk + d*sp2(jat)
      endif

      end
!--------------------------------------------------------------------------------
      subroutine se_pme_calccm1p(iat,iatk,bhk,bk,sp2,c,dq)

      implicit none
#include "sebomd.dim"
#include "se_fulsys.h"
#include "se_xyzcrd.h"
      integer iat, iatk
      double precision bhk, bk, sp2(maxatm), c, dq
      double precision d
      integer jat

      if (iatk.eq.7) then
!       N atom
        c = 0.0d0
        d = -0.0909d0
      elseif (iatk.eq.8) then
!       O atom
        c = 0.0d0
        d = -0.0449d0
      elseif (iatk.eq.9) then
!       F atom
        c = 0.3381d0
        d = 0.0148d0
      elseif (iatk.eq.16) then
!       S atom
        c = -0.0834d0
        d = -0.0848d0
      elseif (iatk.eq.17) then
!       Cl atom
        c = -0.1080d0
        d = -0.1168d0
      elseif (iatk.eq.35) then
!       Br atom
        c = -0.0116d0
        d = -0.0338d0
      elseif (iatk.eq.53) then
!       I atom
        c = -0.3213d0
        d = -0.0636d0
      else
        c = 0.0d0
        d = 0.0d0
      endif
       
      dq = c*atchg(iat) + d + bhk
      atchg2(iat) = atchg2(iat) + dq*bk
!     partial normalization:
      do jat=1,natoms
        atchg2(jat) = atchg2(jat) - sp2(jat)*dq
      enddo
      end
!--------------------------------------------------------------------------------
      subroutine se_pme_getdhk(iatk,jatk,dhk)

      implicit none
      integer iatk, jatk
      double precision dhk

      double precision d

      if (iatk.eq.1) then
!       H atom
        if (jatk.eq.7) then
!         overlap with N atom
          d = 0.1854d0
        elseif (jatk.eq.8) then
!         overlap with O atom
          d = 0.1434d0
        elseif (jatk.eq.14) then
!         overlap with Si atom
          d = -0.1004d0
        else
          d = 0.0d0
        endif
      else
        d=0.0d0 
      endif
       
      dhk=d
      end
