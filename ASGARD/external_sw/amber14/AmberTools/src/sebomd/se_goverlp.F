      subroutine se_goverlp(ia,norbsa,npqa,xisa,xipa,xida,xa,ya,za,
     .                   ib,norbsb,npqb,xisb,xipb,xidb,xb,yb,zb,
     .                   rab,sabdia,dxsabdia,dysabdia,dzsabdia)
!
!     routine to calculate diatomic overlaps for the s-p-d shells of
!     of atoms a and b.  returns s, p, and d overlaps in sabdia.
!
!     input variables:
!
!     ia,ib     = atomic numbers for atoms a and b.
!
!     norbsa,
!     norbsb    = numbers of atomic orbitals centered on atoms a and b.
!
!     npqa,npqb = principal quantum numbers.
!
!     xisa,xisb = orbital exponents for s-type atomic orbitals.
!
!     xipa,xipb = orbital exponents for p-type atomic orbitals.
!
!     xida,xidb = orbital exponents for d-type atomic orbitals.
!
!     xa,ya,za,
!     xb,yb,zb  = cartesian coordinate (angstroms).
!
!     rab       = interatomic distance (angstroms).
!
!
!
!     returned:
!
!     sabdia = diatomic overlap matrix for atoms a and b.  the format of
!              sabdia is shown below.  for example, sabdia(6,4) would
!              contain the overlap between the dxz orbital of atom a and
!              the pz orbital of atom b.
!
!
!     format of sabdia:
!
!              s(b) px(b) py(b) pz(b) dzz(b) dxz(b) dyz(b) dxx-yy(b) dxy(b)
!            --------------------------------------------------------------
!  s(a)      |                    .
!  px(a)     |                    .
!  py(a)     |                    .
!  pz(a)     |                    .
!  dzz(a)    |                    .
!  dxz(a)    | . . . . . . . . .(6,4) = <dxz(a)|pz(b)>
!  dyz(a)    |
!  dxx-yy(a) |
!  dxy(a)    |
!
!
!
!
!     utilizes the following additional subprograms:  sab
!                                                     slimit
!                                                     rtrans
!
!
      implicit none
#include "sebomd.dim"
#include "se_ovraux.h"
#include "se_constants.h"
      integer ia, norbsa, npqa
      integer ib, norbsb, npqb
      double precision xisa, xipa, xida, xa, ya, za
      double precision xisb, xipb, xidb, xb, yb, zb
      double precision rab
      double precision sabdia(9,9)
      double precision dxsabdia(9,9),dysabdia(9,9),dzsabdia(9,9)
!
!     auxiliary inforation for efficient calculation of overlaps:
!
!     exhalf(la,ia) = sqrt(0.5)*xila**(npqa+0.5),
!                     where xila = orbital exponent corresponding to la.
!
!
!     local arrays:
!
      double precision sabloc(9,9),sctemp(9,9),c(9,9)
      integer icode(18)
      double precision dxsctemp(9,9),dysctemp(9,9),dzsctemp(9,9)
      integer nijsc(6),isc(6,81),jsc(6,81),nkijsc(6,81),kijsc(6,81,9),
     .          nkc(6,9),kc(6,9,9)
      double precision dsabloc(9,9),dxsabloc(9,9),dysabloc(9,9)
      double precision dzsabloc(9,9)
      double precision dxc(9,9),dyc(9,9),dzc(9,9)
      logical newexp
!
      logical first
      data first /.true./
!
!     code for the six possible values of norbsa+norbsb:
!
      data icode /0,1,0,0,2,0,0,3,0,4,0,0,5,0,0,0,0,6/
!
      save first,icode,nijsc,isc,jsc,nkijsc,kijsc,nkc,kc

      double precision xmin, xtest, sij, stemp, rnab
      integer i, iswap, itemp, ntemp, norbs
      integer j, ij, k, m
      integer la, nptemp, ncentr
      integer lb, iab, kindx
      double precision xtemp, ytemp, ztemp
      double precision exhsa, exhpa, exhda
      double precision exhsb, exhpb, exhdb
      double precision r, tmpij
      double precision dxtmpij, dxsij, dxstemp
      double precision dytmpij, dysij, dystemp
      double precision dztmpij, dzsij, dzstemp
      double precision, external :: sab
!
!     on the first call, set limits to take advantage of the sparsity
!     of sabloc and c.
!
      if(first)then
        first = .false.
        call se_gslimit(nijsc,isc,jsc,nkijsc,kijsc,nkc,kc)
      endif
!
!     before doing any overlap calculations, check to see if the orbital
!     exponents and interatomic distance are such that the overlaps would
!     be negligible.
!
      xmin = xisa + xisb
      if(norbsa.gt.1) xmin = min(xmin,xipa+xisb)
      if(norbsa.gt.4) xmin = min(xmin,xida+xisb)
      if(norbsb.gt.1) xmin = min(xmin,xisa+xipb)
      if(norbsb.gt.1) xmin = min(xmin,xisa+xidb)
      if(norbsa.gt.1.and.norbsb.gt.1)then
        xmin = min(xmin,xipa+xipb)
        if(norbsa.gt.4) xmin = min(xmin,xida+xipb)
        if(norbsb.gt.4) xmin = min(xmin,xipa+xidb)
        if(norbsa.gt.4.and.norbsb.gt.4) xmin = min(xmin,xida+xidb)
      endif
      xtest = rab*xmin
      if(xtest.gt.25.0d0)then
        do i=1,norbsa
          do j=1,norbsb
            sabdia(i,j) = 0.0d0
             dxsabdia(i,j) = 0.0d0
             dysabdia(i,j) = 0.0d0
             dzsabdia(i,j) = 0.0d0
            sabloc(i,j) = 0.0d0
             dsabloc(i,j) = 0.0d0
             dxsabloc(i,j) = 0.0d0
             dysabloc(i,j) = 0.0d0
             dzsabloc(i,j) = 0.0d0
          enddo
        enddo
        return
      endif
!
!
!     rotation of local overlaps to molecular frame requires that norbsa
!     be greater than or equal to norbsb.  if this isn't the case, then
!     temporarily swap everything.
!
      iswap = 0
      if(norbsb.gt.norbsa)then
        iswap = 1
        itemp = ia
        ia = ib
        ib = itemp
        xtemp = xisa
        xisa = xisb
        xisb = xtemp
        xtemp = xipa
        xipa = xipb
        xipb = xtemp
        xtemp = xida
        xida = xidb
        xidb = xtemp
        xtemp = xa
        xa = xb
        xb = xtemp
        ytemp = ya
        ya = yb
        yb = ytemp
        ztemp = za
        za = zb
        zb = ztemp
        nptemp = npqa
        npqa = npqb
        npqb = nptemp
        ntemp = norbsa
        norbsa = norbsb
        norbsb = ntemp
      endif
!
!     now compute the nonzero local overlaps sabloc.  these are defined
!     in a frame wherein atom b is directly above atom a (i.e., xa=xb,
!     ya=yb, and zb-za=rab).  the local frame overlaps are then rotated
!     to yield the molecular frame overlaps sabdia.
!

      r = rab/bohr2ang
      ncentr = 2
      exhsa = exhalf(0,ia)
      exhpa = exhalf(1,ia)
      exhda = exhalf(2,ia)
      exhsb = exhalf(0,ib)
      exhpb = exhalf(1,ib)
      exhdb = exhalf(2,ib)
      rnab = r**(npqa+npqb+1)
!
!     <s(a)|s(b)>:
!
      la = 0
      lb = 0
      m = 0
      newexp = .true.

       call se_gfctsab(ncentr,npqa,la,xisa,exhsa,
     .                         npqb,lb,xisb,exhsb,
     .                         newexp,m,r,rnab,sabloc(1,1),dsabloc(1,1))
!
!     <pz(a)|s(b)>:
!
      if(norbsa.gt.1)then
        la = 1
        lb = 0
        m = 0
        newexp = .true.

        call se_gfctsab(ncentr,npqa,la,xipa,exhpa,
     .                      npqb,lb,xisb,exhsb,
     .                      newexp,m,r,rnab,sabloc(4,1),dsabloc(4,1))
      endif
!
!     <dzz(a)|s(b)>:
!
      if(norbsa.gt.4)then
        la = 2
        lb = 0
        m = 0
        newexp = .true.

        call se_gfctsab(ncentr,npqa,la,xida,exhda,
     .                         npqb,lb,xisb,exhsb,
     .                         newexp,m,r,rnab,sabloc(5,1),dsabloc(5,1))
      endif
!
!     <dxz(a)|px(b)>, <dyz(a)|py(b)>, and <dzz(a)|pz(b)>:
!
      if(norbsa.gt.4.and.norbsb.gt.1)then
        la = 2
        lb = 1
        m = 1
        newexp = .true.
        call se_gfctsab(ncentr,npqa,la,xida,exhda,
     .                         npqb,lb,xipb,exhpb,
     .                         newexp,m,r,rnab,sabloc(6,2),dsabloc(6,2))
        sabloc(7,3) = sabloc(6,2)
         dsabloc(7,3) = dsabloc(6,2)
        m = 0
        newexp = .false.

        call se_gfctsab(ncentr,npqa,la,xida,exhda,
     .                         npqb,lb,xipb,exhpb,
     .                         newexp,m,r,rnab,sabloc(5,4),dsabloc(5,4))
          sabloc(5,4) = -sabloc(5,4)
           dsabloc(5,4) = -dsabloc(5,4)
      endif
!
!     <s(a)|pz(b)>:
!
      if(norbsb.gt.1)then
        if(ia.eq.ib)then
          sabloc(1,4) = -sabloc(4,1)
           dsabloc(1,4) = -dsabloc(4,1)
        else
          la = 0
          lb = 1
          m = 0
          newexp = .true.
        call se_gfctsab(ncentr,npqa,la,xisa,exhsa,
     .                         npqb,lb,xipb,exhpb,
     .                         newexp,m,r,rnab,sabloc(1,4),dsabloc(1,4))
              sabloc(1,4) = -sabloc(1,4)
               dsabloc(1,4) = -dsabloc(1,4)
        endif
      endif
!
!     <s(a)|dzz(b)>:
!
      if(norbsb.gt.4)then
        if(ia.eq.ib)then
          sabloc(1,5) = sabloc(5,1)
           dsabloc(1,5) = dsabloc(5,1)
        else
          la = 0
          lb = 2
          m = 0
          newexp = .true.
        call se_gfctsab(ncentr,npqa,la,xisa,exhsa,
     .                        npqb,lb,xidb,exhdb,
     .                        newexp,m,r,rnab,sabloc(1,5),dsabloc(1,5))
        endif
      endif
!
!     <px(a)|dxz(b)>, <py(a)|dyz(b)>, and <pz(a)|dzz(b)>:
!
      if(norbsa.gt.1.and.norbsb.gt.4)then
        if(ia.eq.ib)then
          sabloc(2,6) = -sabloc(6,2)
           dsabloc(2,6) = -dsabloc(6,2)
          sabloc(4,5) = -sabloc(5,4)
           dsabloc(4,5) = -dsabloc(5,4)
        else
          la = 1
          lb = 2
          m = 1
          newexp = .true.
        call se_gfctsab(ncentr,npqa,la,xipa,exhpa,
     .                         npqb,lb,xidb,exhdb,
     .                         newexp,m,r,rnab,sabloc(2,6),dsabloc(2,6))
             sabloc(2,6) = -sabloc(2,6)
              dsabloc(2,6) = -dsabloc(2,6)
          m = 0
          newexp = .false.
        call se_gfctsab(ncentr,npqa,la,xipa,exhpa,
     .                        npqb,lb,xidb,exhdb,
     .                        newexp,m,r,rnab,sabloc(4,5),dsabloc(4,5))
        endif
        sabloc(3,7) = sabloc(2,6)
         dsabloc(3,7) = dsabloc(2,6)
      endif
!
!     <px(a)|px(b)>, <py(a)|py(b)>, <pz(a)|pz(b)>:
!
      if(norbsa.gt.1.and.norbsb.gt.1)then
        la = 1
        lb = 1
        m = 1
        newexp = .true.
        call se_gfctsab(ncentr,npqa,la,xipa,exhpa,
     .                         npqb,lb,xipb,exhpb,
     .                         newexp,m,r,rnab,sabloc(2,2),dsabloc(2,2))
        newexp = .false.
        sabloc(3,3) = sabloc(2,2)
         dsabloc(3,3) = dsabloc(2,2)
        m = 0
        call se_gfctsab(ncentr,npqa,la,xipa,exhpa,
     .                         npqb,lb,xipb,exhpb,
     .                         newexp,m,r,rnab,sabloc(4,4),dsabloc(4,4))
            sabloc(4,4) = -sabloc(4,4)
             dsabloc(4,4) = -dsabloc(4,4)
      endif
!
!     the five <d(a)|d(b)> overlaps:
!
      if(norbsa.gt.4.and.norbsb.gt.4)then
        la = 2
        lb = 2
        m = 0
        newexp = .true.

        call se_gfctsab(ncentr,npqa,la,xida,exhda,
     .                        npqb,lb,xidb,exhdb,
     .                        newexp,m,r,rnab,sabloc(5,5),dsabloc(5,5))

        m = 1
        newexp = .false.
        call se_gfctsab(ncentr,npqa,la,xida,exhda,
     .                        npqb,lb,xidb,exhdb,
     .                        newexp,m,r,rnab,sabloc(6,6),dsabloc(6,6))
           sabloc(6,6) = -sabloc(6,6)
            dsabloc(6,6) = -dsabloc(6,6)

        sabloc(7,7) = sabloc(6,6)
         dsabloc(7,7) = dsabloc(6,6)
        m = 2
        call se_gfctsab(ncentr,npqa,la,xida,exhda,
     .                        npqb,lb,xidb,exhdb,
     .                        newexp,m,r,rnab,sabloc(8,8),dsabloc(8,8))
        sabloc(9,9) = sabloc(8,8)
         dsabloc(9,9) = dsabloc(8,8)
      endif
!
!       d        d
!       --  -->  --  
!       dr       dx
!
        if(iswap.eq.1)then
         do i=1,9
          do j=1,9
          dxsabloc(j,i)=-dsabloc(j,i)*(xa-xb)/rab/bohr2ang
          dysabloc(j,i)=-dsabloc(j,i)*(ya-yb)/rab/bohr2ang
          dzsabloc(j,i)=-dsabloc(j,i)*(za-zb)/rab/bohr2ang
          enddo
         enddo
        else
         do i=1,9
          do j=1,9
          dxsabloc(j,i)=dsabloc(j,i)*(xa-xb)/rab/bohr2ang
          dysabloc(j,i)=dsabloc(j,i)*(ya-yb)/rab/bohr2ang
          dzsabloc(j,i)=dsabloc(j,i)*(za-zb)/rab/bohr2ang
          enddo
         enddo
        endif
!
!     now determine a transformation matrix c that relates the local frame
!     overlaps to the true diatomic overlaps.  the matrix c is computed
!     from a rotation matrix that relates the true coordinate system to
!     a coordinate system which is consistent with the local frame.
!
      norbs = max(norbsa,norbsb)
      if(norbs.eq.1)then
        sabdia(1,1) = sabloc(1,1)
         dxsabloc(1,1)=dsabloc(1,1)*(xa-xb)/rab/bohr2ang
         dysabloc(1,1)=dsabloc(1,1)*(ya-yb)/rab/bohr2ang
         dzsabloc(1,1)=dsabloc(1,1)*(za-zb)/rab/bohr2ang
         dxsabdia(1,1)=dxsabloc(1,1)
         dysabdia(1,1)=dysabloc(1,1)
         dzsabdia(1,1)=dzsabloc(1,1)
      else
        call se_grtrans(xa,ya,za,xb,yb,zb,rab,norbs,c,dxc,dyc,dzc,iswap)
!
!       if there are d orbitals then compute the molecular frame overlaps
!       using the general formula,
!
!                        sabdia = c'*sabloc*c.
!
!       if we only have s and p orbitals, then do rotations by hand.
!
        if(norbsa.gt.4.or.norbsb.gt.4)then

!
!         use general formula.
!
          iab = icode(norbsa+norbsb)
          do ij=1,nijsc(iab)
            i = isc(iab,ij)
            j = jsc(iab,ij)
            tmpij = 0.0d0
             dxtmpij = 0.0d0
             dytmpij = 0.0d0
             dztmpij = 0.0d0
            do k=1,nkijsc(iab,ij)
              kindx = kijsc(iab,ij,k)
              tmpij = tmpij + sabloc(i,kindx)*c(kindx,j)
               dxtmpij = dxtmpij + dxsabloc(i,kindx)*c(kindx,j)
     .                           + sabloc(i,kindx)*dxc(kindx,j)
               dytmpij = dytmpij + dysabloc(i,kindx)*c(kindx,j)
     .                           + sabloc(i,kindx)*dyc(kindx,j)
               dztmpij = dztmpij + dzsabloc(i,kindx)*c(kindx,j)
     .                           + sabloc(i,kindx)*dzc(kindx,j)
               
            enddo
            sctemp(i,j) = tmpij
             dxsctemp(i,j) = dxtmpij
             dysctemp(i,j) = dytmpij
             dzsctemp(i,j) = dztmpij
          enddo
          do i=1,norbsa
            do j=1,norbsb
              sij = 0.0d0
               dxsij = 0.0d0
               dysij = 0.0d0
               dzsij = 0.0d0
              do k=1,nkc(iab,i)
                kindx = kc(iab,i,k)
                sij = sij + c(kindx,i)*sctemp(kindx,j)
                 dxsij = dxsij + dxc(kindx,i)*sctemp(kindx,j)
     .                     + c(kindx,i)*dxsctemp(kindx,j)
                 dysij = dysij + dyc(kindx,i)*sctemp(kindx,j)
     .                     + c(kindx,i)*dysctemp(kindx,j)
                 dzsij = dzsij + dzc(kindx,i)*sctemp(kindx,j)
     .                     + c(kindx,i)*dzsctemp(kindx,j)
              enddo
              sabdia(i,j) = sij
               dxsabdia(i,j) = dxsij
               dysabdia(i,j) = dysij
               dzsabdia(i,j) = dzsij
            enddo
          enddo
        else



!
!         manual rotations.
!
          sabdia(1,1) = sabloc(1,1)
           dxsabdia(1,1) = dxsabloc(1,1)
           dysabdia(1,1) = dysabloc(1,1)
           dzsabdia(1,1) = dzsabloc(1,1)
          if(norbsb.gt.1)then
           sabdia(1,2) = sabloc(1,4)*c(4,2)
            dxsabdia(1,2) = dxsabloc(1,4)*c(4,2) + sabloc(1,4)*dxc(4,2)
            dysabdia(1,2) = dysabloc(1,4)*c(4,2) + sabloc(1,4)*dyc(4,2)
            dzsabdia(1,2) = dzsabloc(1,4)*c(4,2) + sabloc(1,4)*dzc(4,2)
           sabdia(1,3) = sabloc(1,4)*c(4,3)
            dxsabdia(1,3) = dxsabloc(1,4)*c(4,3) + sabloc(1,4)*dxc(4,3)
            dysabdia(1,3) = dysabloc(1,4)*c(4,3) + sabloc(1,4)*dyc(4,3)
            dzsabdia(1,3) = dzsabloc(1,4)*c(4,3) + sabloc(1,4)*dzc(4,3)
           sabdia(1,4) = sabloc(1,4)*c(4,4)
            dxsabdia(1,4) = dxsabloc(1,4)*c(4,4) + sabloc(1,4)*dxc(4,4)
            dysabdia(1,4) = dysabloc(1,4)*c(4,4) + sabloc(1,4)*dyc(4,4)
            dzsabdia(1,4) = dzsabloc(1,4)*c(4,4) + sabloc(1,4)*dzc(4,4)
          endif
          if(norbsa.gt.1)then
            sabdia(2,1) = sabloc(4,1)*c(4,2)
             dxsabdia(2,1) = dxsabloc(4,1)*c(4,2) + sabloc(4,1)*dxc(4,2)
             dysabdia(2,1) = dysabloc(4,1)*c(4,2) + sabloc(4,1)*dyc(4,2)
             dzsabdia(2,1) = dzsabloc(4,1)*c(4,2) + sabloc(4,1)*dzc(4,2)
            sabdia(3,1) = sabloc(4,1)*c(4,3)
             dxsabdia(3,1) = dxsabloc(4,1)*c(4,3) + sabloc(4,1)*dxc(4,3)
             dysabdia(3,1) = dysabloc(4,1)*c(4,3) + sabloc(4,1)*dyc(4,3)
             dzsabdia(3,1) = dzsabloc(4,1)*c(4,3) + sabloc(4,1)*dzc(4,3)
            sabdia(4,1) = sabloc(4,1)*c(4,4)
             dxsabdia(4,1) = dxsabloc(4,1)*c(4,4) + sabloc(4,1)*dxc(4,4)
             dysabdia(4,1) = dysabloc(4,1)*c(4,4) + sabloc(4,1)*dyc(4,4)
             dzsabdia(4,1) = dzsabloc(4,1)*c(4,4) + sabloc(4,1)*dzc(4,4)
          endif
          if(norbsa.gt.1.and.norbsb.gt.1)then
            sabdia(2,2) = sabloc(2,2)*c(2,2)*c(2,2)
     .                  + sabloc(3,3)*c(3,2)*c(3,2)
     .                  + sabloc(4,4)*c(4,2)*c(4,2)

             dxsabdia(2,2) = dxsabloc(2,2)*c(2,2)*c(2,2)
     .                   + sabloc(2,2)*dxc(2,2)*c(2,2)
     .                   + sabloc(2,2)*c(2,2)*dxc(2,2)
     .                   + dxsabloc(3,3)*c(3,2)*c(3,2)
     .                   + sabloc(3,3)*dxc(3,2)*c(3,2)
     .                   + sabloc(3,3)*c(3,2)*dxc(3,2)
     .                   + dxsabloc(4,4)*c(4,2)*c(4,2)
     .                   + sabloc(4,4)*dxc(4,2)*c(4,2)
     .                   + sabloc(4,4)*c(4,2)*dxc(4,2)

             dysabdia(2,2) = dysabloc(2,2)*c(2,2)*c(2,2)
     .                   + sabloc(2,2)*dyc(2,2)*c(2,2)
     .                   + sabloc(2,2)*c(2,2)*dyc(2,2)
     .                   + dysabloc(3,3)*c(3,2)*c(3,2)
     .                   + sabloc(3,3)*dyc(3,2)*c(3,2)
     .                   + sabloc(3,3)*c(3,2)*dyc(3,2)
     .                   + dysabloc(4,4)*c(4,2)*c(4,2)
     .                   + sabloc(4,4)*dyc(4,2)*c(4,2)
     .                   + sabloc(4,4)*c(4,2)*dyc(4,2)

             dzsabdia(2,2) = dzsabloc(2,2)*c(2,2)*c(2,2)
     .                   + sabloc(2,2)*dzc(2,2)*c(2,2)
     .                   + sabloc(2,2)*c(2,2)*dzc(2,2)
     .                   + dzsabloc(3,3)*c(3,2)*c(3,2)
     .                   + sabloc(3,3)*dzc(3,2)*c(3,2)
     .                   + sabloc(3,3)*c(3,2)*dzc(3,2)
     .                   + dzsabloc(4,4)*c(4,2)*c(4,2)
     .                   + sabloc(4,4)*dzc(4,2)*c(4,2)
     .                   + sabloc(4,4)*c(4,2)*dzc(4,2)

!
            sabdia(2,3) = sabloc(2,2)*c(2,2)*c(2,3)
     .                  + sabloc(3,3)*c(3,2)*c(3,3)
     .                  + sabloc(4,4)*c(4,2)*c(4,3)

            dxsabdia(2,3) = dxsabloc(2,2)*c(2,2)*c(2,3)
     .                  + sabloc(2,2)*dxc(2,2)*c(2,3)
     .                  + sabloc(2,2)*c(2,2)*dxc(2,3)
     .                  + dxsabloc(3,3)*c(3,2)*c(3,3)
     .                  + sabloc(3,3)*dxc(3,2)*c(3,3)
     .                  + sabloc(3,3)*c(3,2)*dxc(3,3)
     .                  + dxsabloc(4,4)*c(4,2)*c(4,3)
     .                  + sabloc(4,4)*dxc(4,2)*c(4,3)
     .                  + sabloc(4,4)*c(4,2)*dxc(4,3)

            dysabdia(2,3) = dysabloc(2,2)*c(2,2)*c(2,3)
     .                  + sabloc(2,2)*dyc(2,2)*c(2,3)
     .                  + sabloc(2,2)*c(2,2)*dyc(2,3)
     .                  + dysabloc(3,3)*c(3,2)*c(3,3)
     .                  + sabloc(3,3)*dyc(3,2)*c(3,3)
     .                  + sabloc(3,3)*c(3,2)*dyc(3,3)
     .                  + dysabloc(4,4)*c(4,2)*c(4,3)
     .                  + sabloc(4,4)*dyc(4,2)*c(4,3)
     .                  + sabloc(4,4)*c(4,2)*dyc(4,3)

            dzsabdia(2,3) = dzsabloc(2,2)*c(2,2)*c(2,3)
     .                  + sabloc(2,2)*dzc(2,2)*c(2,3)
     .                  + sabloc(2,2)*c(2,2)*dzc(2,3)
     .                  + dzsabloc(3,3)*c(3,2)*c(3,3)
     .                  + sabloc(3,3)*dzc(3,2)*c(3,3)
     .                  + sabloc(3,3)*c(3,2)*dzc(3,3)
     .                  + dzsabloc(4,4)*c(4,2)*c(4,3)
     .                  + sabloc(4,4)*dzc(4,2)*c(4,3)
     .                  + sabloc(4,4)*c(4,2)*dzc(4,3)
!
            sabdia(2,4) = sabloc(2,2)*c(2,2)*c(2,4)
     .                  + sabloc(3,3)*c(3,2)*c(3,4)
     .                  + sabloc(4,4)*c(4,2)*c(4,4)

            dxsabdia(2,4) = dxsabloc(2,2)*c(2,2)*c(2,4)
     .                  + sabloc(2,2)*dxc(2,2)*c(2,4)
     .                  + sabloc(2,2)*c(2,2)*dxc(2,4)
     .                  + dxsabloc(3,3)*c(3,2)*c(3,4)
     .                  + sabloc(3,3)*dxc(3,2)*c(3,4)
     .                  + sabloc(3,3)*c(3,2)*dxc(3,4)
     .                  + dxsabloc(4,4)*c(4,2)*c(4,4)
     .                  + sabloc(4,4)*dxc(4,2)*c(4,4)
     .                  + sabloc(4,4)*c(4,2)*dxc(4,4)

            dysabdia(2,4) = dysabloc(2,2)*c(2,2)*c(2,4)
     .                  + sabloc(2,2)*dyc(2,2)*c(2,4)
     .                  + sabloc(2,2)*c(2,2)*dyc(2,4)
     .                  + dysabloc(3,3)*c(3,2)*c(3,4)
     .                  + sabloc(3,3)*dyc(3,2)*c(3,4)
     .                  + sabloc(3,3)*c(3,2)*dyc(3,4)
     .                  + dysabloc(4,4)*c(4,2)*c(4,4)
     .                  + sabloc(4,4)*dyc(4,2)*c(4,4)
     .                  + sabloc(4,4)*c(4,2)*dyc(4,4)

            dzsabdia(2,4) = dzsabloc(2,2)*c(2,2)*c(2,4)
     .                  + sabloc(2,2)*dzc(2,2)*c(2,4)
     .                  + sabloc(2,2)*c(2,2)*dzc(2,4)
     .                  + dzsabloc(3,3)*c(3,2)*c(3,4)
     .                  + sabloc(3,3)*dzc(3,2)*c(3,4)
     .                  + sabloc(3,3)*c(3,2)*dzc(3,4)
     .                  + dzsabloc(4,4)*c(4,2)*c(4,4)
     .                  + sabloc(4,4)*dzc(4,2)*c(4,4)
     .                  + sabloc(4,4)*c(4,2)*dzc(4,4)


!
            sabdia(3,2) = sabdia(2,3)
              dxsabdia(3,2) = dxsabdia(2,3)
              dysabdia(3,2) = dysabdia(2,3)
              dzsabdia(3,2) = dzsabdia(2,3)
!
            sabdia(3,3) = sabloc(2,2)*c(2,3)*c(2,3)
     .                  + sabloc(3,3)*c(3,3)*c(3,3)
     .                  + sabloc(4,4)*c(4,3)*c(4,3)

            dxsabdia(3,3) = dxsabloc(2,2)*c(2,3)*c(2,3)
     .                  + sabloc(2,2)*dxc(2,3)*c(2,3)
     .                  + sabloc(2,2)*c(2,3)*dxc(2,3)
     .                  + dxsabloc(3,3)*c(3,3)*c(3,3)
     .                  + sabloc(3,3)*dxc(3,3)*c(3,3)
     .                  + sabloc(3,3)*c(3,3)*dxc(3,3)
     .                  + dxsabloc(4,4)*c(4,3)*c(4,3)
     .                  + sabloc(4,4)*dxc(4,3)*c(4,3)
     .                  + sabloc(4,4)*c(4,3)*dxc(4,3)

            dysabdia(3,3) = dysabloc(2,2)*c(2,3)*c(2,3)
     .                  + sabloc(2,2)*dyc(2,3)*c(2,3)
     .                  + sabloc(2,2)*c(2,3)*dyc(2,3)
     .                  + dysabloc(3,3)*c(3,3)*c(3,3)
     .                  + sabloc(3,3)*dyc(3,3)*c(3,3)
     .                  + sabloc(3,3)*c(3,3)*dyc(3,3)
     .                  + dysabloc(4,4)*c(4,3)*c(4,3)
     .                  + sabloc(4,4)*dyc(4,3)*c(4,3)
     .                  + sabloc(4,4)*c(4,3)*dyc(4,3)

            dzsabdia(3,3) = dzsabloc(2,2)*c(2,3)*c(2,3)
     .                  + sabloc(2,2)*dzc(2,3)*c(2,3)
     .                  + sabloc(2,2)*c(2,3)*dzc(2,3)
     .                  + dzsabloc(3,3)*c(3,3)*c(3,3)
     .                  + sabloc(3,3)*dzc(3,3)*c(3,3)
     .                  + sabloc(3,3)*c(3,3)*dzc(3,3)
     .                  + dzsabloc(4,4)*c(4,3)*c(4,3)
     .                  + sabloc(4,4)*dzc(4,3)*c(4,3)
     .                  + sabloc(4,4)*c(4,3)*dzc(4,3)
!
            sabdia(3,4) = sabloc(2,2)*c(2,3)*c(2,4)
     .                  + sabloc(3,3)*c(3,3)*c(3,4)
     .                  + sabloc(4,4)*c(4,3)*c(4,4)

            dxsabdia(3,4) = dxsabloc(2,2)*c(2,3)*c(2,4)
     .                  + sabloc(2,2)*dxc(2,3)*c(2,4)
     .                  + sabloc(2,2)*c(2,3)*dxc(2,4)
     .                  + dxsabloc(3,3)*c(3,3)*c(3,4)
     .                  + sabloc(3,3)*dxc(3,3)*c(3,4)
     .                  + sabloc(3,3)*c(3,3)*dxc(3,4)
     .                  + dxsabloc(4,4)*c(4,3)*c(4,4)
     .                  + sabloc(4,4)*dxc(4,3)*c(4,4)
     .                  + sabloc(4,4)*c(4,3)*dxc(4,4)

            dysabdia(3,4) = dysabloc(2,2)*c(2,3)*c(2,4)
     .                  + sabloc(2,2)*dyc(2,3)*c(2,4)
     .                  + sabloc(2,2)*c(2,3)*dyc(2,4)
     .                  + dysabloc(3,3)*c(3,3)*c(3,4)
     .                  + sabloc(3,3)*dyc(3,3)*c(3,4)
     .                  + sabloc(3,3)*c(3,3)*dyc(3,4)
     .                  + dysabloc(4,4)*c(4,3)*c(4,4)
     .                  + sabloc(4,4)*dyc(4,3)*c(4,4)
     .                  + sabloc(4,4)*c(4,3)*dyc(4,4)

            dzsabdia(3,4) = dzsabloc(2,2)*c(2,3)*c(2,4)
     .                  + sabloc(2,2)*dzc(2,3)*c(2,4)
     .                  + sabloc(2,2)*c(2,3)*dzc(2,4)
     .                  + dzsabloc(3,3)*c(3,3)*c(3,4)
     .                  + sabloc(3,3)*dzc(3,3)*c(3,4)
     .                  + sabloc(3,3)*c(3,3)*dzc(3,4)
     .                  + dzsabloc(4,4)*c(4,3)*c(4,4)
     .                  + sabloc(4,4)*dzc(4,3)*c(4,4)
     .                  + sabloc(4,4)*c(4,3)*dzc(4,4)
!
            sabdia(4,2) = sabdia(2,4)
             dxsabdia(4,2) = dxsabdia(2,4)
             dysabdia(4,2) = dysabdia(2,4)
             dzsabdia(4,2) = dzsabdia(2,4)

            sabdia(4,3) = sabdia(3,4)
             dxsabdia(4,3) = dxsabdia(3,4)
             dysabdia(4,3) = dysabdia(3,4)
             dzsabdia(4,3) = dzsabdia(3,4)
!
            sabdia(4,4) = sabloc(2,2)*c(2,4)*c(2,4)
     .                  + sabloc(3,3)*c(3,4)*c(3,4)
     .                  + sabloc(4,4)*c(4,4)*c(4,4)

            dxsabdia(4,4) = dxsabloc(2,2)*c(2,4)*c(2,4)
     .                  + sabloc(2,2)*dxc(2,4)*c(2,4)
     .                  + sabloc(2,2)*c(2,4)*dxc(2,4)
     .                  + dxsabloc(3,3)*c(3,4)*c(3,4)
     .                  + sabloc(3,3)*dxc(3,4)*c(3,4)
     .                  + sabloc(3,3)*c(3,4)*dxc(3,4)
     .                  + dxsabloc(4,4)*c(4,4)*c(4,4)
     .                  + sabloc(4,4)*dxc(4,4)*c(4,4)
     .                  + sabloc(4,4)*c(4,4)*dxc(4,4)

            dysabdia(4,4) = dysabloc(2,2)*c(2,4)*c(2,4)
     .                  + sabloc(2,2)*dyc(2,4)*c(2,4)
     .                  + sabloc(2,2)*c(2,4)*dyc(2,4)
     .                  + dysabloc(3,3)*c(3,4)*c(3,4)
     .                  + sabloc(3,3)*dyc(3,4)*c(3,4)
     .                  + sabloc(3,3)*c(3,4)*dyc(3,4)
     .                  + dysabloc(4,4)*c(4,4)*c(4,4)
     .                  + sabloc(4,4)*dyc(4,4)*c(4,4)
     .                  + sabloc(4,4)*c(4,4)*dyc(4,4)

            dzsabdia(4,4) = dzsabloc(2,2)*c(2,4)*c(2,4)
     .                  + sabloc(2,2)*dzc(2,4)*c(2,4)
     .                  + sabloc(2,2)*c(2,4)*dzc(2,4)
     .                  + dzsabloc(3,3)*c(3,4)*c(3,4)
     .                  + sabloc(3,3)*dzc(3,4)*c(3,4)
     .                  + sabloc(3,3)*c(3,4)*dzc(3,4)
     .                  + dzsabloc(4,4)*c(4,4)*c(4,4)
     .                  + sabloc(4,4)*dzc(4,4)*c(4,4)
     .                  + sabloc(4,4)*c(4,4)*dzc(4,4)
          endif
       endif
      endif
!
!     if swap was made previously, return input variables to original
!     values and transpose the diatomic overlap matrix.
!
      if(iswap.eq.1)then
        itemp = ia
        ia = ib
        ib = itemp
        xtemp = xisa
        xisa = xisb
        xisb = xtemp
        xtemp = xipa
        xipa = xipb
        xipb = xtemp
        xtemp = xida
        xida = xidb
        xidb = xtemp
        xtemp = xa
        xa = xb
        xb = xtemp
        ytemp = ya
        ya = yb
        yb = ytemp
        ztemp = za
        za = zb
        zb = ztemp
        nptemp = npqa
        npqa = npqb
        npqb = nptemp
        ntemp = norbsa
        norbsa = norbsb
        norbsb = ntemp
        do i=2,norbsb
          do j=1,i-1
            stemp = sabdia(i,j)
             dxstemp = dxsabdia(i,j)
             dystemp = dysabdia(i,j)
             dzstemp = dzsabdia(i,j)
            sabdia(i,j) = sabdia(j,i)
             dxsabdia(i,j) = dxsabdia(j,i)
             dysabdia(i,j) = dysabdia(j,i)
             dzsabdia(i,j) = dzsabdia(j,i)
            sabdia(j,i) = stemp
             dxsabdia(j,i) = dxstemp 
             dysabdia(j,i) = dystemp 
             dzsabdia(j,i) = dzstemp 
          enddo
        enddo
      endif
      return
      end
!--------------------------------------------------------------------------------
       subroutine se_gfctsab(ncentr,na,la,xia,exhla,
     .                   nb,lb,xib,exhlb,
     .                   newexp,m,rab,rnab,sab,dsab)
!
!     portable overlap function.  computes the overlap integral between
!     slater-type orbitals a and b by transformation to prolate spheroidal
!     coordinates.  in terms of cartesian coordinates, this means that
!     the xyz coordinate systems of centers a and b are mirror images,
!     with the z-axis of a pointing toward b and vice versa.
!
!     handles overlap for orbitals of type s, p, d, f, ..., etc.
!
!
!     input:
!
!     ncentr = the number of centers involved (ncentr=2 means orbitals
!              a and b are centered on different atoms, ncentr=1 means
!              they are centered on the same atom.  default=2).
!     na,nb = principal quantum numbers for the orbitals.
!     la,lb = angular momentum quantum numbers for the orbitals.
!     xia,xib = exponents for the radial terms in the sto.
!     exhla = sqrt(0.5)*xia**(na+0.5)
!     exhlb = sqrt(0.5)*xib**(na+0.5)
!     newexp = .true. if either xia or xib has changed from previous
!              call for a given atom, .false. otherwise.
!     m = magnetic quantum number (must be the same for both orbitals).
!     rab = distance between the atomic centers in atomic units
!           (atomic units = angtroms/bohr2ang)
! ken changed back to original mopac value:
!            (atomic units = angtroms/0.529167)
!     rnab = rab**(na+nb+1)
!
!
!     returned:
!
!     sab = overlap between orbitals a and b.
!
!
!
!     sto's are of the form:
!
!     sto(n,l,m,xi,r,theta,phi) = norm*r**(n-1)*exp(-xi*r)*y(l,m,theta,phi);
!
!     where,
!
!     n,l,m = quantum numbers;
!     xi = orbital exponent;
!     r,theta,phi = spherical polar coordinates;
!     norm = radial function normalization factor;
!     y(l,m,theta,phi) = normalized, complex spherical harmonic function.
!
!
!     the formulation used is taken directly from:
!     j. p. p. stewart, journal of computer-aided molecular design,
!     vol. 4, pp 1-105 (1990).  the final formula for the overlap
!     appears on page 28 of this article.  one correction was necessary
!     for the exponent iexp -- see below in code.
!
!     programmed by s. l. dixon, oct., 1991.
!
!
      implicit none
      integer ncentr, m
      integer na, la
      integer nb, lb
      double precision xia, exhla
      double precision xib, exhlb
      double precision rab, rnab
      logical newexp
      double precision fact(0:15),cterm(0:4,0:4,0:2),binom(0:10,0:10),
     .          aterm(0:15),bterm(0:15),abterm(0:15,0:15),neg1(-15:15),
     .          daterm(0:15),dbterm(0:15),dabterm(0:15,0:15)
      save rplus,rminus,abterm,neg1
!
!     initialize factorial values (up to 15!).
!
      data fact /1.0d0,1.0d0,2.0d0,6.0d0,24.0d0,120.0d0,720.0d0,
     .           5040.0d0,40320.0d0,362880.0d0,36288.0d2,399168.0d2,
     .           4790016.0d2,62270208.0d2,871782912.0d2,1307674368.0d3/
!
!     initialize spherical harmonic normalization factors -- up to
!     an angular momentum quantum number of 4.  corresponds to the
!     term c(l,m,j) of stewart reference.
!
      data cterm(0,0,0) / 1.000000000000d0/
      data cterm(1,0,0) / 1.000000000000d0/
      data cterm(1,1,0) / 0.707106781187d0/
      data cterm(2,0,0) / 1.500000000000d0/
      data cterm(2,0,1) /-0.500000000000d0/
      data cterm(2,1,0) / 1.224744871392d0/
      data cterm(2,2,0) / 0.612372435696d0/
      data cterm(3,0,0) / 2.500000000000d0/
      data cterm(3,0,1) /-1.500000000000d0/
      data cterm(3,1,0) / 2.165063509461d0/
      data cterm(3,1,1) /-0.433012701892d0/
      data cterm(3,2,0) / 1.369306393763d0/
      data cterm(3,3,0) / 0.559016994375d0/
      data cterm(4,0,0) / 4.375000000000d0/
      data cterm(4,0,1) /-3.750000000000d0/
      data cterm(4,0,2) / 0.375000000000d0/
      data cterm(4,1,0) / 3.913118960625d0/
      data cterm(4,1,1) /-1.677050983125d0/
      data cterm(4,2,0) / 2.766992952647d0/
      data cterm(4,2,1) /-0.395284707521d0/
      data cterm(4,3,0) / 1.479019945775d0/
      data cterm(4,4,0) / 0.522912516584d0/
!
!     initialize binomial coefficients:
!
!           binom(n1,n2) = (n1!)/((n1-n2)!(n2!))  for n1.ge.n2
!
!     assign values up to n1=10.
!
      data binom /11*1.0d0,0.0d0,1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0,
     .            7.0d0,8.0d0,9.0d0,10.0d0,2*0.0d0,1.0d0,3.0d0,6.0d0,
     .            10.0d0,15.0d0,21.0d0,28.0d0,36.0d0,45.0d0,3*0.0d0,
     .            1.0d0,4.0d0,10.0d0,20.0d0,35.0d0,56.0d0,84.0d0,
     .            120.0d0,4*0.0d0,1.0d0,5.0d0,15.0d0,35.0d0,70.0d0,
     .            126.0d0,210.0d0,5*0.0d0,1.0d0,6.0d0,21.0d0,56.0d0,
     .            126.0d0,252.0d0,6*0.0d0,1.0d0,7.0d0,28.0d0,84.0d0,
     .            210.0d0,7*0.0d0,1.0d0,8.0d0,36.0d0,120.0d0,8*0.0d0,
     .            1.0d0,9.0d0,45.0d0,9*0.0d0,1.0d0,10.0d0,10*0.0d0,
     .            1.0d0/
!
!     initialize (-1)**i for i=-15,-14,...,0,...15.
!
      data neg1 /-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,
     .           -1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1/
!
      double precision sqrtl(0:7,0:7),sqrtf(0:7,0:7)
      logical first
      data first /.true./
      save fact,cterm,sqrtl,sqrtf,binom,first

      integer ia, ja, ipa, ja2, ipamax
      integer ib, jb, ipb, jb2, ipbmax
      integer msave, nab, i, j
      integer nasub, nbsub, nbsub0
      integer iexp, iexp0, iexp2, iexp3, iexp4, iexp6
      integer iqa, ka, jamax, iqamax, ipam0, iqam0
      integer iqb, kb, jbmax, iqbmax, ipbm0, iqbm0
      integer nasub3, nasub5, nasub6, nasub7
      integer nbsub1, nbsub2, nbsub4, nbsub5, nbsub6, nbsub7
      double precision terml, termf, rplus, rminus, r, explus
      double precision termxi, termnl
      double precision sum1, sum2, t45, t345, ab, radial, t2345
      double precision term1, term2, term3, term4, term5, t12345
      double precision sumja, sumqa, sumka
      double precision sumjb, sumqb
      double precision dexplus, dsum1, dsum2, sab, dsab, dab
      double precision dsumja, dsumqa, dsumqb, dsumka, dsumjb

      if(first)then
        first = .false.
        do ia=0,7
          do ib=0,7
            terml = (2*ia+1)*(2*ib+1)
            sqrtl(ia,ib) = dsqrt(terml)
            termf = fact(2*ia)*fact(2*ib)
            sqrtf(ia,ib) = dsqrt(termf)
          enddo
        enddo
      endif
!
!     if the calling routine has supplied a negative magnetic quantum
!     number, temporarily change it to the equivalent positive number.
!
      if(ncentr.ne.1.and.ncentr.ne.2) ncentr = 2
      msave = m
      if(m.lt.0) m = -m
      r = rab
!
!     some frequently used quantities:
!
      if(newexp)then
        rplus = r*(xia + xib)*0.5d0
        rminus = r*(xia - xib)*0.5d0
      endif
!
!     if the orbitals are on different centers, compute the overlap
!     according to j.p.p. stewart's formula.  note that term1, term2,
!     term3, term4, and term5 correspond, in order, to the first 5
!     terms in the innermost sum on p. 28 of the reference.  the term
!     (-1)**iexp has been corrected and differs from stewart's formula.
!     aterm and bterm correspond to the a and b integrals on p. 28.
!
      if(ncentr.eq.2)then
        if(na.eq.1.and.nb.eq.1.and.abs(rminus).lt.1.0d-6)then
!
!         use simple formula for hydrogen-hydrogen overlaps.
!
          explus = exp(-rplus)
           dexplus = -exp(-rplus)*(xia + xib)*0.5d0
          sum1 = 1.0d0/(3.0d0*rplus)
           dsum1 = -1.0d0/(3.0d0*rplus**2)*(xia + xib)*0.5d0
          sum2 = (rplus + 1.0d0)/rplus**3
           dsum2 = -(xia + xib)*0.5d0*(2.0d0/rplus**3+3.0d0/rplus**4)
          sab = explus*(sum1 + sum2)*(xia**3)*rnab
           dsab = dexplus*(sum1 + sum2)*(xia**3)*rnab
     .          + explus*(dsum1 + dsum2)*(xia**3)*rnab
     .          + explus*(sum1 + sum2)*(xia**3)*(na+nb+1)*rnab/r
        else
!
!         at least one non-hydrogen involved.
!
!
!         assign a-type and b-type integrals if xia or xib has changed.
!
          nab = na + nb
          if(newexp)then
            call se_gabintg(nab,rplus,rminus,aterm,bterm,daterm,dbterm)
            do i=0,nab
              do j=0,nab
                abterm(i,j) = aterm(i)*bterm(j)
                 dabterm(i,j) = daterm(i)*bterm(j)*(xia + xib)*0.5d0
     .                        + aterm(i)*dbterm(j)*(xia - xib)*0.5d0
              enddo
            enddo
          endif
!
!         compute sums.
!
          if(la.eq.0.and.lb.eq.0)then
!
!           s(a)-s(b) overlap.
!
            sumja = 0.0d0
            dsumja = 0.0d0
            do ipa=0,na
              t45 = binom(na,ipa)
              do ipb=0,nb
                t345 = binom(nb,ipb)*t45
                nasub = ipa + ipb
                nbsub = nab - nasub
                iexp = nb - ipb
                ab = abterm(nasub,nbsub)
                 dab = dabterm(nasub,nbsub)
                sumja = sumja + t345*ab*neg1(iexp)
                 dsumja = dsumja + t345*dab*neg1(iexp)
              enddo
            enddo

          elseif(la.eq.1.and.lb.eq.0)then
!
!           p(a)-s(b) overlap.
!
            sumja = 0.0d0
             dsumja = 0.0d0
            nbsub0 = nab - 1
            iexp0 = nb
            do ipa=0,na-1
              t45 = binom(na-1,ipa)
              nasub5 = ipa
              nbsub5 = nbsub0 - ipa
              do ipb=0,nb
                t345 = binom(nb,ipb)*t45
                nasub6 = nasub5 + ipb
                nbsub6 = nbsub5 - ipb
                iexp6 = iexp0 - ipb
                sumqa = 0.0d0
                 dsumqa = 0.0d0
                do iqa=0,1
                  nasub = nasub6 + iqa
                  nbsub = nbsub6 + iqa
                  iexp = iexp6
                  ab = abterm(nasub,nbsub)
                   dab = dabterm(nasub,nbsub)
                  sumqa = sumqa + ab*neg1(iexp)
                   dsumqa = dsumqa + dab*neg1(iexp)
                enddo
                sumja = sumja + t345*sumqa
                 dsumja = dsumja + t345*dsumqa
              enddo
            enddo
          elseif(la.eq.0.and.lb.eq.1)then
!
!           s(a)-p(b) overlap.
!
            sumja = 0.0d0
             dsumja = 0.0d0
            nbsub0 = nab - 1
            iexp0 = nb - 1
            do ipa=0,na
              t45 = binom(na,ipa)
              nasub5 = ipa
              nbsub5 = nbsub0 - ipa
              do ipb=0,nb-1
                t345 = binom(nb-1,ipb)*t45
                nasub6 = nasub5 + ipb
                nbsub6 = nbsub5 - ipb
                iexp6 = iexp0 - ipb
                sumqb = 0.0d0
                 dsumqb = 0.0d0
                do iqb=0,1
                  nasub = nasub6 + iqb
                  nbsub = nbsub6 + iqb
                  iexp = iexp6 + iqb
                  ab = abterm(nasub,nbsub)
                   dab = dabterm(nasub,nbsub)
                  sumqb = sumqb + ab*neg1(iexp)
                   dsumqb = dsumqb + dab*neg1(iexp)
                enddo
                sumja = sumja + sumqb*t345
                 dsumja = dsumja + dsumqb*t345
              enddo
            enddo
          elseif(la.eq.1.and.lb.eq.1)then
!
!           p(a)-p(b) overlap.
!
            nbsub0 = nab - 2
            iexp0 = m + nb - 1
            sumka = 0.0d0
             dsumka = 0.0d0
            do ka=0,m
              nasub3 = 2*ka
              iexp3 = iexp0 - ka
              do kb=0,m
                nbsub4 = nbsub0 + 2*kb
                iexp4 = iexp3 + kb
                do ipa=0,na-1
                  t45 = binom(na-1,ipa)
                  nasub5 = nasub3 + ipa
                  nbsub5 = nbsub4 - ipa
                  do ipb=0,nb-1
                    t345 = binom(nb-1,ipb)*t45
                    nasub6 = nasub5 + ipb
                    nbsub6 = nbsub5 - ipb
                    iexp6 = iexp4 - ipb
                    sumqa = 0.0d0
                     dsumqa = 0.0d0
                    do iqa=0,la-m
                      nasub7 = nasub6 + iqa
                      nbsub7 = nbsub6 + iqa
                      do iqb=0,lb-m
                        nasub = nasub7 + iqb
                        nbsub = nbsub7 + iqb
                        iexp = iexp6 + iqb
                        ab = abterm(nasub,nbsub)
                         dab = dabterm(nasub,nbsub)
                        sumqa = sumqa + ab*neg1(iexp)
                         dsumqa = dsumqa + dab*neg1(iexp)
                      enddo
                    enddo
                    sumka = sumka + sumqa*t345
                     dsumka = dsumka + dsumqa*t345
                  enddo
                enddo
              enddo
            enddo
            sumja = sumka*cterm(1,m,0)**2
             dsumja = dsumka*cterm(1,m,0)**2
          else
!
!           use general formula from stewart.
!
            jamax = (la-m)/2
            jbmax = (lb-m)/2
            sumja = 0.0d0
             dsumja = 0.0d0
            nbsub0 = na + nb - la - lb
            iexp0 = m + nb - lb
            ipam0 = na - la
            ipbm0 = nb - lb
            iqam0 = la - m
            iqbm0 = lb - m
            do ja=0,jamax
              ja2 = 2*ja
              ipamax = ipam0 + ja2
              iqamax = iqam0 - ja2
              nbsub1 = nbsub0 + ja2
              sumjb = 0.0d0
               dsumjb = 0.0d0
              do jb=0,jbmax
                jb2 = 2*jb
                ipbmax = ipbm0 + jb2
                iqbmax = iqbm0 - jb2
                nbsub2 = nbsub1 + jb2
                iexp2 = iexp0 + jb2
                sumka = 0.0d0
                 dsumka = 0.0d0
                do ka=0,m
                  nasub3 = 2*ka
                  iexp3 = iexp2 - ka
                  do kb=0,m
                    term5 = binom(m,ka)*binom(m,kb)
                    nbsub4 = nbsub2 + 2*kb
                    iexp4 = iexp3 + kb
                    do ipa=0,ipamax
                      term4 = binom(ipamax,ipa)
                      t45 = term4*term5
                      nasub5 = nasub3 + ipa
                      nbsub5 = nbsub4 - ipa
                      do ipb=0,ipbmax
                        term3 = binom(ipbmax,ipb)
                        t345 = term3*t45
                        nasub6 = nasub5 + ipb
                        nbsub6 = nbsub5 - ipb
                        iexp6 = iexp4 - ipb
                        do iqa=0,iqamax
                          term2 = binom(iqamax,iqa)
                          t2345 = term2*t345
                          nasub7 = nasub6 + iqa
                          nbsub7 = nbsub6 + iqa
                          do iqb=0,iqbmax
                            term1 = binom(iqbmax,iqb)
                            t12345 = term1*t2345
                            nasub = nasub7 + iqb
                            nbsub = nbsub7 + iqb
                            iexp = iexp6 + iqb
!
!                           note overall quantities:
!
!                           nasub = 2*ka + ipa + ipb + iqa + iqb
!
!                           nbsub = 2*kb + na - la + 2*ja + nb - lb
!                                 + 2*jb - ipa - ipb + iqa + iqb
!
!                           iexp = m - ka + kb + nb - lb + 2*jb - ipb + iqb
!
                            ab = abterm(nasub,nbsub)
                             dab = dabterm(nasub,nbsub)
                            sumka = sumka + t12345*ab*neg1(iexp)
                             dsumka = dsumka + t12345*dab*neg1(iexp)
                          enddo
                        enddo
                      enddo
                    enddo
                  enddo
                enddo
                sumjb = sumjb + sumka*cterm(lb,m,jb)
                 dsumjb = dsumjb + dsumka*cterm(lb,m,jb)
              enddo
              sumja = sumja + sumjb*cterm(la,m,ja)
               dsumja = dsumja + dsumjb*cterm(la,m,ja)
            enddo
          endif
          termxi = exhla*exhlb
          termnl = sqrtl(la,lb)/sqrtf(na,nb)
          sab = termxi*termnl*rnab*sumja
           dsab = termxi*termnl*rnab*sumja*(na+nb+1)/r
     .          + termxi*termnl*rnab*dsumja
        endif
      else
!
!       the orbitals are centered on the same atom.  since sto's do
!       not have radial nodes, the radial functions are not orthogonal,
!       and therefore the orbitals may not be orthogonal.
!
        if(la.ne.lb)then
!
!         orbitals are orthogonal due to their orthogonal spherical
!         harmonics.
!
          sab = 0.0d0
           dsab = 0.0d0
        else
!
!         orbitals have the same total angular momentum.  compute the
!         overlap based on an integration over the radial portions of
!         the sto's.
!
          radial = fact(na+nb)/((xia+xib)**(na+nb+1))
          termxi = (2**(na+1))*exhla*(2**(nb+1))*exhlb
!
!         take into account the fact that the spherical harmonics could
!         be 180 degrees out of phase since the coordinate systems are
!         mirror images.
!
          sab = neg1(la+m)*radial*termxi/sqrtf(na,nb)
           dsab = 0.0d0
        endif
      endif
!
!     return the magnetic quantum number back to its original value
!     if it was changed.
!
      m = msave
      return
      end
!--------------------------------------------------------------------------------
      subroutine se_gabintg(n,xa,xb,aterm,bterm,daterm,dbterm)
!
!     computes a-type and b-type integrals.
!
      implicit none
      integer n
      double precision xa, xb
      double precision aterm(0:15),bterm(0:15)
      double precision daterm(0:15),dbterm(0:15)
!
!     local variables:
!
      double precision fact(0:15)
!
!     initialize factorial values (up to 15!).
!
      data fact /1.0d0,1.0d0,2.0d0,6.0d0,24.0d0,120.0d0,720.0d0,
     .           5040.0d0,40320.0d0,362880.0d0,36288.0d2,399168.0d2,
     .           4790016.0d2,62270208.0d2,871782912.0d2,1307674368.0d3/
      save fact
!
!     on first call, initialize some quantities used repetitively in the
!     b-type integrals.
!
      double precision xbexp(15),bfactr(0:7,0:15),bterm0(0:15)
      integer neg1(0:16)
      double precision dxbexp(15)
      logical first
      data first /.true./
!
!     initialize (-1)**i for i=0,1,...,15
!
      data neg1 /1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1/
!
      save first,bfactr,bterm0,neg1

      integer i, j, ij1, jlast
      double precision xainv, axb, xbinv, exb, exbinv, bi
      double precision dxbinv, dexb, dexbinv, dbi

      if(first)then
        first = .false.
        do i=0,7
          do j=0,15
            ij1 = i+j+1
            bfactr(i,j) = 2*mod(ij1,2)/(fact(j)*ij1)
          enddo
        enddo
        do i=0,15
          bterm0(i) = (1.0d0 - dfloat(neg1(i+1)))/dfloat(i+1)
        enddo
      endif
!
      xainv = 1.0d0/xa
      aterm(0) = xainv*exp(-xa)
       daterm(0) = -xainv*exp(-xa) - exp(-xa)*xainv**2
      do i=1,n
        aterm(i) = aterm(0) + i*aterm(i-1)*xainv
         daterm(i) = daterm(0) + i*daterm(i-1)*xainv
     .                         - i*aterm(i-1)*xainv**2
      enddo
!
!     the direct computation of bterm can be unstable, depending upon
!     the size of xb.  if abs(xb) is large or very small, then there
!     is no problem.  however, for intermediate sized values of abs(xb),
!     a taylor expansion for the integral must be used.
!
      axb = abs(xb)
      if(axb.lt.1.0d-6)then
        do i=0,n
          bterm(i) = bterm0(i)
           dbterm(i) = 0.0d0
        enddo
      elseif(axb.gt.3.0d0)then
        xbinv = 1.0d0/xb
         dxbinv = -1.0d0/xb**2
        exb = exp(xb)
         dexb = exp(xb)
        exbinv = 1.0d0/exb
         dexbinv = -exbinv  
        bterm(0) = xbinv*(exb - exbinv)
         dbterm(0) = dxbinv*(exb - exbinv) + xbinv*(dexb - dexbinv)
        do i=1,n
          bterm(i) = xbinv*(i*bterm(i-1) + neg1(i)*exb - exbinv)
           dbterm(i) = dxbinv*(i*bterm(i-1) + neg1(i)*exb - exbinv)
     .               + xbinv*(i*dbterm(i-1) + neg1(i)*dexb - dexbinv)
        enddo
      else
        if(axb.gt.2.0d0)then
          jlast = 15
        elseif(axb.gt.1.0d0)then
          jlast = 12
        elseif(axb.gt.0.5d0)then
          jlast = 7
        else
          jlast = 6
        endif
        xbexp(1) = -xb
         dxbexp(1) = -1.0d0
        do j=2,jlast
          xbexp(j) = -xb*xbexp(j-1)
           dxbexp(j) = -xb*dxbexp(j-1) - xbexp(j-1)
        enddo
        do i=0,n
          bi = bfactr(i,0)
           dbi = 0.0d0
          do j=1,jlast
            bi = bi + xbexp(j)*bfactr(i,j)
             dbi = dbi + dxbexp(j)*bfactr(i,j)
          enddo
          bterm(i) = bi
           dbterm(i) = dbi
        enddo
      endif
      return
      end
!--------------------------------------------------------------------------------
      subroutine se_gslimit(nijsc,isc,jsc,nkijsc,kijsc,nkc,kc)
!
!     sets up array limits to take advantage of sparsity of matrices
!     used in diatomic overlap calculation.  the ultimate goal is
!     to compute c'*sabloc*c using the minimum number of floating
!     point operations.  the degree of sparsity depends depends on
!     how many atomic orbitals are on each atom.  for norbsa.ge.norbsb,
!     the combinations are coded according to the variable nab:
!
!                        norbsa  norbsb  nab
!                           1       1     1
!                           4       1     2
!                           4       4     3
!                           9       1     4
!                           9       4     5
!                           9       9     6
!
!     once the value of nab is established, the sparsity variables are
!     defined:
!
!     nijsc(nab) = the number of nonzero elements in sabloc*c.
!
!     *************************************************************************
!     *                                                                       *
!     * for ij=1,nijsc(nab)                                                   *
!     *                                                                       *
!     *   isc(nab,ij) = row number of ijth nonzero element of sabloc*c        *
!     *                                                                       *
!     *   jsc(nab,ij) = column number of ijth nonzero element of sabloc*c     *
!     *                                                                       *
!     *   nkijsc(nab,ij) = number multiplications required to compute the     *
!     *                    (isc(nab,ij),jsc(nab,ij)) entry of sabloc*c        *
!     *                                                                       *
!     *   kijsc(nab,ij,k) = inner index for elements required to compute      *
!     *                     the (isc(nab,ij),jsc(nab,ij)) entry of            *
!     *                     sabloc*c, i.e.,                                   *
!     *                                                                       *
!     *                           i=isc(nab,ij)                               *
!     *                           j=jsc(nab,ij)                               *
!     *                           scij = 0.0                                  *
!     *                           do k=1,nkijsc(nab,ij)                       *
!     *                             kindx = kijsc(nab,ij,k)                   *
!     *                             scij = scij + sabloc(i,kindx)*c(kindx,j)  *
!     *                           enddo                                       *
!     *                           (i,j) elment of sabloc*c = scij             *
!     *                                                                       *
!     *************************************************************************
!
!     nkc(nab,j) = the number of nonzero elements in column j of c.
!     (kc(nab,j,k), k=1,nkc(nab,j)) = list of row numbers in column j
!                                     of c that have nonzero entries.
!
      implicit none
 
      integer nijsc(6),isc(6,81),jsc(6,81),nkijsc(6,81),kijsc(6,81,9),
     .          nkc(6,9),kc(6,9,9)
!
!     local arrays:
!
      integer is(9,9),ic(9,9),norbs(3)
      integer iab, na, nb, ij, norbsb, norbsa, i, j, ijflag, k
!
!     flags for non-zero entries of local frame overlap matrix sabloc:
!
      data is /1,0,0,1,1,0,0,0,0,
     .         0,1,0,0,0,1,0,0,0,
     .         0,0,1,0,0,0,1,0,0,
     .         1,0,0,1,1,0,0,0,0,
     .         1,0,0,1,1,0,0,0,0,
     .         0,1,0,0,0,1,0,0,0,
     .         0,0,1,0,0,0,1,0,0,
     .         0,0,0,0,0,0,0,1,0,
     .         0,0,0,0,0,0,0,0,1/
!
!     flags for non-zero entries of transformation matrix c:
!
      data ic /1,0,0,0,0,0,0,0,0,
     .         0,1,1,1,0,0,0,0,0,
     .         0,1,1,1,0,0,0,0,0,
     .         0,1,1,1,0,0,0,0,0,
     .         0,0,0,0,1,1,1,1,1,
     .         0,0,0,0,1,1,1,1,1,
     .         0,0,0,0,1,1,1,1,1,
     .         0,0,0,0,1,1,1,1,1,
     .         0,0,0,0,1,1,1,1,1/
!
!     possible numbers of orbitals on an atom:
!
      data norbs /1,4,9/
!
      save is,ic,norbs
!
!     determine which entries of sabloc*c need to be computed based on
!     the non-zero entries of sabloc and c.  consider each possible
!     combination of norbsa and norbsb with norbsa.ge.norbsb.
!
      iab = 0
      do na=1,3
        norbsa = norbs(na)
        do nb=1,na
          iab = iab + 1
          ij = 0
          norbsb = norbs(nb)
          do i=1,norbsa
            nkc(iab,i) = 0
            do j=1,norbsa
              if(ic(j,i).eq.1)then
                nkc(iab,i) = nkc(iab,i) + 1
                kc(iab,i,nkc(iab,i)) = j
              endif
              ijflag = 0
              do k=1,norbsb
                if(is(i,k).ne.0.and.ic(k,j).ne.0)then
!
!                 (i,j) entry of c*sabloc is nonzero.
!
                  if(ijflag.eq.0)then
                    ijflag = 1
                    ij = ij + 1
                    nkijsc(iab,ij) = 0
                  endif
                  nkijsc(iab,ij) = nkijsc(iab,ij) + 1
                  kijsc(iab,ij,nkijsc(iab,ij)) = k
                endif
              enddo
              if(ijflag.eq.1)then
                isc(iab,ij) = i
                jsc(iab,ij) = j
              endif
            enddo
          enddo
          nijsc(iab) = ij
        enddo
      enddo
      return
      end
!--------------------------------------------------------------------------------
      subroutine se_grtrans(xa,ya,za,xb,yb,zb,rab,norbs,
     .                            c,dxc,dyc,dzc,iswap)
!
!     computes a matrix c that transforms local frame overlaps sabloc
!     to the molecular frame overlaps sabdia:
!
!     sabdia = c'*sabloc*c.
!
!     c will be either 4 by 4 or 9 by 9, depending on whether norbs=4
!     or norbs=9.
!
      implicit none
      integer norbs, iswap
      double precision xa, ya, za
      double precision xb, yb, zb
      double precision rab
      double precision c(9,9)
      double precision dxc(9,9),dyc(9,9),dzc(9,9)
!
!     local arrays:
!
      double precision rot(3,3),d(5,3,3)
      double precision sq2, sq6
      integer i, j, k, k4
      double precision dxd(5,3,3),dyd(5,3,3),dzd(5,3,3)
      double precision dxrot(3,3),dyrot(3,3),dzrot(3,3)
!
!     sqrt(1/2) and sqrt(1/6):
!
      data sq2,sq6 /0.707106769084930420d0, 0.408248305320739746d0/
      save sq2,sq6
!
!     coefficient for s-type orbitals not effected by directionality:
!
      c(1,1) = 1.0d0
       dxc(1,1) = 0.0d0
       dyc(1,1) = 0.0d0
       dzc(1,1) = 0.0d0
!
!     compute a rotation matrix rot that relates the x, y, and z axes
!     of atoms a and b to a coordinate system (x*,y*,z*), where z* points
!     from atom a to atom b.  the directions of x* and y* can be defined
!     arbitrarily, as long as the system (x*,y*,z*) is right-handed.
!
      call se_grmatrx(xa,ya,za,xb,yb,zb,rab,rot,dxrot,dyrot,dzrot,iswap)
!
!     now assign coefficients for p-type orbitals using the fact that
!     the px, py, and pz orbitals transform just like the x, y, and z
!     axes.
!
      do i=1,3
        do j=1,3
          c(i+1,j+1) = rot(i,j)
           dxc(i+1,j+1) = dxrot(i,j)
           dyc(i+1,j+1) = dyrot(i,j)
           dzc(i+1,j+1) = dzrot(i,j)
        enddo
      enddo
!
!     if atoms a and b have no d orbitals we are done.
!
      if(norbs.le.4) return
!
!     now form five symmetric tensors to represent the d(z**2), d(xz),
!     d(yz), d(x**2-y**2), and d(xy) orbitals in the (x*,y*,z*) coordinate
!     system.  this will involve taking dyad products of column vectors.
!     for column vectors x and y, the dyad product is defined as:
!
!        x dyad y  =  x*y'  =  outer product = tensor matrix.
!
      do i=1,3
        do j=1,3
!
!         d(z**2) = (3*(z dyad z) - i)/sqrt(6)  (i = identity tensor):
!
          d(1,i,j) = rot(i,3)*rot(j,3)*3.0d0*sq6
           dxd(1,i,j) = dxrot(i,3)*rot(j,3)*3.0d0*sq6
     .                + rot(i,3)*dxrot(j,3)*3.0d0*sq6
           dyd(1,i,j) = dyrot(i,3)*rot(j,3)*3.0d0*sq6
     .                + rot(i,3)*dyrot(j,3)*3.0d0*sq6
           dzd(1,i,j) = dzrot(i,3)*rot(j,3)*3.0d0*sq6
     .                + rot(i,3)*dzrot(j,3)*3.0d0*sq6
          if(i.eq.j) then
           d(1,i,j) = d(1,i,j) - sq6
            dxd(1,i,j) = dxd(1,i,j)
            dyd(1,i,j) = dyd(1,i,j)
            dzd(1,i,j) = dzd(1,i,j)
          endif
!
!         d(xz) = (x dyad z  +  z dyad x)/sqrt(2):
!
          d(2,i,j) = (rot(i,1)*rot(j,3) + rot(i,3)*rot(j,1))*sq2
           dxd(2,i,j) = (dxrot(i,1)*rot(j,3) + dxrot(i,3)*rot(j,1))*sq2
     .              + (rot(i,1)*dxrot(j,3) + rot(i,3)*dxrot(j,1))*sq2
           dyd(2,i,j) = (dyrot(i,1)*rot(j,3) + dyrot(i,3)*rot(j,1))*sq2
     .              + (rot(i,1)*dyrot(j,3) + rot(i,3)*dyrot(j,1))*sq2
           dzd(2,i,j) = (dzrot(i,1)*rot(j,3) + dzrot(i,3)*rot(j,1))*sq2
     .              + (rot(i,1)*dzrot(j,3) + rot(i,3)*dzrot(j,1))*sq2
!
!         d(yz) = (y dyad z  +  z dyad y)/sqrt(2):
!
          d(3,i,j) = (rot(i,2)*rot(j,3) + rot(i,3)*rot(j,2))*sq2
           dxd(3,i,j) = (dxrot(i,2)*rot(j,3) + dxrot(i,3)*rot(j,2))*sq2
     .              + (rot(i,2)*dxrot(j,3) + rot(i,3)*dxrot(j,2))*sq2
           dyd(3,i,j) = (dyrot(i,2)*rot(j,3) + dyrot(i,3)*rot(j,2))*sq2
     .              + (rot(i,2)*dyrot(j,3) + rot(i,3)*dyrot(j,2))*sq2
           dzd(3,i,j) = (dzrot(i,2)*rot(j,3) + dzrot(i,3)*rot(j,2))*sq2
     .              + (rot(i,2)*dzrot(j,3) + rot(i,3)*dzrot(j,2))*sq2

!
!         d(x**2-y**2) = (x dyad x  -  y dyad y)/sqrt(2):
!
          d(4,i,j) = (rot(i,1)*rot(j,1) - rot(i,2)*rot(j,2))*sq2
           dxd(4,i,j) = (dxrot(i,1)*rot(j,1) - dxrot(i,2)*rot(j,2))*sq2
     .              + (rot(i,1)*dxrot(j,1) - rot(i,2)*dxrot(j,2))*sq2
           dyd(4,i,j) = (dyrot(i,1)*rot(j,1) - dyrot(i,2)*rot(j,2))*sq2
     .              + (rot(i,1)*dyrot(j,1) - rot(i,2)*dyrot(j,2))*sq2
           dzd(4,i,j) = (dzrot(i,1)*rot(j,1) - dzrot(i,2)*rot(j,2))*sq2
     .              + (rot(i,1)*dzrot(j,1) - rot(i,2)*dzrot(j,2))*sq2
!
!         d(xy) = (x dyad y  +  y dyad x)/sqrt(2):
!
          d(5,i,j) = (rot(i,1)*rot(j,2) + rot(i,2)*rot(j,1))*sq2
           dxd(5,i,j) = (dxrot(i,1)*rot(j,2) + dxrot(i,2)*rot(j,1))*sq2
     .             + (rot(i,1)*dxrot(j,2) + rot(i,2)*dxrot(j,1))*sq2
           dyd(5,i,j) = (dyrot(i,1)*rot(j,2) + dyrot(i,2)*rot(j,1))*sq2
     .             + (rot(i,1)*dyrot(j,2) + rot(i,2)*dyrot(j,1))*sq2
           dzd(5,i,j) = (dzrot(i,1)*rot(j,2) + dzrot(i,2)*rot(j,1))*sq2
     .             + (rot(i,1)*dzrot(j,2) + rot(i,2)*dzrot(j,1))*sq2
        enddo
      enddo
!
!     now get the projections of the above tensors on the corresponding
!     tensors formed by using the standard basis:
!
!                     1          0          0
!                x* = 0     y* = 1     z* = 0
!                     0          0          1
!
!     this requires a scalar product of two tensors.  for tensors a and
!     b, the scalar product (a|b) is defined as:
!
!            (a|b) = (b|a) = sum(i=1,3;j=1,3: a(i,j)*b(i,j))
!
      do k=1,5
        k4 = k+4
      c(5,k4) = -d(k,1,1)*sq6 - d(k,2,2)*sq6 + d(k,3,3)*2.0*sq6
       dxc(5,k4) = -dxd(k,1,1)*sq6 - dxd(k,2,2)*sq6 + dxd(k,3,3)*2.0*sq6
       dyc(5,k4) = -dyd(k,1,1)*sq6 - dyd(k,2,2)*sq6 + dyd(k,3,3)*2.0*sq6
       dzc(5,k4) = -dzd(k,1,1)*sq6 - dzd(k,2,2)*sq6 + dzd(k,3,3)*2.0*sq6
      c(6,k4) = d(k,1,3)*sq2 + d(k,3,1)*sq2
       dxc(6,k4) = dxd(k,1,3)*sq2 + dxd(k,3,1)*sq2
       dyc(6,k4) = dyd(k,1,3)*sq2 + dyd(k,3,1)*sq2
       dzc(6,k4) = dzd(k,1,3)*sq2 + dzd(k,3,1)*sq2
      c(7,k4) = d(k,2,3)*sq2 + d(k,3,2)*sq2
       dxc(7,k4) = dxd(k,2,3)*sq2 + dxd(k,3,2)*sq2
       dyc(7,k4) = dyd(k,2,3)*sq2 + dyd(k,3,2)*sq2
       dzc(7,k4) = dzd(k,2,3)*sq2 + dzd(k,3,2)*sq2
      c(8,k4) = d(k,1,1)*sq2 - d(k,2,2)*sq2
       dxc(8,k4) = dxd(k,1,1)*sq2 - dxd(k,2,2)*sq2
       dyc(8,k4) = dyd(k,1,1)*sq2 - dyd(k,2,2)*sq2
       dzc(8,k4) = dzd(k,1,1)*sq2 - dzd(k,2,2)*sq2
      c(9,k4) = d(k,1,2)*sq2 + d(k,2,1)*sq2
       dxc(9,k4) = dxd(k,1,2)*sq2 + dxd(k,2,1)*sq2 
       dyc(9,k4) = dyd(k,1,2)*sq2 + dyd(k,2,1)*sq2 
       dzc(9,k4) = dzd(k,1,2)*sq2 + dzd(k,2,1)*sq2 
      enddo
      return
      end
!--------------------------------------------------------------------------------
      subroutine se_grmatrx(xa,ya,za,xb,yb,zb,rab,rot,dxrot,dyrot,dzrot
     &                  ,iswap)
!
!     constructs a rotation matrix rot that aligns the z axis with the
!     the vector (xb-xa,yb-ya,zb-za).
!
      implicit none
      double precision xa, ya, za
      double precision xb, yb, zb
      double precision rab
      double precision rot(3,3)
      double precision dxrot(3,3),dyrot(3,3),dzrot(3,3)
      double precision dxx(3),dxy(3),dxz(3),dyx(3),dyy(3),
     &          dyz(3),dzx(3),dzy(3),dzz(3)

      double precision x(3),y(3),z(3)
      double precision zinv, z12sqr, z12, yinv
      double precision dxzinv, dyzinv, dzzinv
      double precision dxyinv, dyyinv, dzyinv
      double precision dxz12sqr, dyz12sqr, dzz12sqr
      double precision dxz12, dyz12, dzz12
      integer j, iswap
    
       z(1)= -(xa - xb)
       z(2)= -(ya - yb)
       z(3)= -(za - zb)
       zinv = 1/rab
       if (iswap.ne.1) then
       dxzinv = -(xa-xb)/rab**3
       dyzinv = -(ya-yb)/rab**3
       dzzinv = -(za-zb)/rab**3
      z(1) = z(1)*zinv
       dxz(1)= -(xa - xb)*dxzinv - zinv
       dyz(1)= -(xa - xb)*dyzinv
       dzz(1)= -(xa - xb)*dzzinv
      z(2) = z(2)*zinv
       dxz(2)= -(ya - yb)*dxzinv
       dyz(2)= -(ya - yb)*dyzinv - zinv
       dzz(2)= -(ya - yb)*dzzinv
      z(3) = z(3)*zinv
       dxz(3)= -(za - zb)*dxzinv
       dyz(3)= -(za - zb)*dyzinv
       dzz(3)= -(za - zb)*dzzinv - zinv
       else
       dxzinv = (xa-xb)/rab**3
       dyzinv = (ya-yb)/rab**3
       dzzinv = (za-zb)/rab**3
      z(1) = z(1)*zinv
       dxz(1)= -(xa - xb)*dxzinv + zinv
       dyz(1)= -(xa - xb)*dyzinv
       dzz(1)= -(xa - xb)*dzzinv
      z(2) = z(2)*zinv
       dxz(2)= -(ya - yb)*dxzinv
       dyz(2)= -(ya - yb)*dyzinv + zinv
       dzz(2)= -(ya - yb)*dzzinv
      z(3) = z(3)*zinv
       dxz(3)= -(za - zb)*dxzinv
       dyz(3)= -(za - zb)*dyzinv
       dzz(3)= -(za - zb)*dzzinv + zinv
       endif

      if(abs(z(3)).gt.0.99999999d0)then
!
!       atoms a and b lie along the z axis.  construct a zero or 180 degree
!       rotation about the x axis.
!             
        z(1) = 0.0d0
         dxz(1) = 0.0d0
         dyz(1) = 0.0d0
         dzz(1) = 0.0d0
        z(2) = 0.0d0
         dxz(2) = 0.0d0
         dyz(2) = 0.0d0
         dzz(2) = 0.0d0
        z(3) = sign(1.0d0,z(3))
         dxz(3) = 0.0d0
         dyz(3) = 0.0d0
         dzz(3) = 0.0d0
        x(1) = 1.0d0
         dxx(1) = 0.0d0
         dyx(1) = 0.0d0
         dzx(1) = 0.0d0
        x(2) = 0.0d0
         dxx(2) = 0.0d0
         dyx(2) = 0.0d0
         dzx(2) = 0.0d0
        x(3) = 0.0d0
         dxx(3) = 0.0d0
         dyx(3) = 0.0d0
         dzx(3) = 0.0d0
        y(1) = 0.0d0
         dxy(1) = 0.0d0
         dyy(1) = 0.0d0
         dzy(1) = 0.0d0
        y(2) = z(3)
         dxy(2) = 0.0d0
         dyy(2) = 0.0d0
         dzy(2) = 0.0d0
        y(3) = 0.0d0
         dxy(3) = 0.0d0
         dyy(3) = 0.0d0
         dzy(3) = 0.0d0
      else
        z12sqr = z(1)**2 + z(2)**2
         dxz12sqr = 2.0d0*(z(1)*dxz(1)+z(2)*dxz(2))
         dyz12sqr = 2.0d0*(z(1)*dyz(1)+z(2)*dyz(2))
         dzz12sqr = 2.0d0*(z(1)*dzz(1)+z(2)*dzz(2))
        z12 = dsqrt(z12sqr)
         dxz12 = dxz12sqr/(2.0d0*z12)
         dyz12 = dyz12sqr/(2.0d0*z12)
         dzz12 = dzz12sqr/(2.0d0*z12)
        x(1) = -z(2)/z12
         dxx(1) = -dxz(2)/z12 + z(2)*dxz12/z12**2
         dyx(1) = -dyz(2)/z12 + z(2)*dyz12/z12**2
         dzx(1) = -dzz(2)/z12 + z(2)*dzz12/z12**2
        x(2) = z(1)/z12
         dxx(2) = dxz(1)/z12 - z(1)*dxz12/z12**2
         dyx(2) = dyz(1)/z12 - z(1)*dyz12/z12**2
         dzx(2) = dzz(1)/z12 - z(1)*dzz12/z12**2
        x(3) = 0.0d0
         dxx(3) = 0.0d0
         dyx(3) = 0.0d0
         dzx(3) = 0.0d0
        if(abs(z(3)).lt.1.0d-8)then
          y(1) = 0.0d0
           dxy(1) = 0.0d0
           dyy(1) = 0.0d0
           dzy(1) = 0.0d0
          y(2) = 0.0d0
           dxy(2) = 0.0d0
           dyy(2) = 0.0d0
           dzy(2) = 0.0d0
          y(3) = 1.0d0
           dxy(3) = 0.0d0
           dyy(3) = 0.0d0
           dzy(3) = 0.0d0

        else
          yinv = 1.0d0/(z12*dsqrt(1.0d0 + z12sqr/z(3)**2))
             
           dxyinv = -(dxz12*dsqrt(1.0d0 + z12sqr/z(3)**2)
     &                     +z12*0.5d0*(-2.0d0*z12sqr*dxz(3)/z(3)**3
     &           + dxz12sqr/z(3)**2)/dsqrt(1.0d0 + z12sqr/z(3)**2))
     &             /(z12*dsqrt(1.0d0 + z12sqr/z(3)**2))**2

           dyyinv = -(dyz12*dsqrt(1.0d0 + z12sqr/z(3)**2)
     &                     +z12*0.5d0*(-2.0d0*z12sqr*dyz(3)/z(3)**3
     &           + dyz12sqr/z(3)**2)/dsqrt(1.0d0 + z12sqr/z(3)**2))
     &             /(z12*dsqrt(1.0d0 + z12sqr/z(3)**2))**2

           dzyinv = -(dzz12*dsqrt(1.0d0 + z12sqr/z(3)**2)
     &                     +z12*0.5d0*(-2.0d0*z12sqr*dzz(3)/z(3)**3
     &           + dzz12sqr/z(3)**2)/dsqrt(1.0d0 + z12sqr/z(3)**2))
     &             /(z12*dsqrt(1.0d0 + z12sqr/z(3)**2))**2


          y(1) = z(1)*yinv
          dxy(1) = dxz(1)*yinv + z(1)*dxyinv
          dyy(1) = dyz(1)*yinv + z(1)*dyyinv
          dzy(1) = dzz(1)*yinv + z(1)*dzyinv   
          y(2) = z(2)*yinv
          dxy(2) = dxz(2)*yinv + z(2)*dxyinv
          dyy(2) = dyz(2)*yinv + z(2)*dyyinv
          dzy(2) = dzz(2)*yinv + z(2)*dzyinv  
          y(3) = -(z12sqr/z(3))*yinv
          dxy(3) = -(dxz12sqr/z(3))*yinv
     &              +(z12sqr*dxz(3)/z(3)**2)*yinv
     &              - (z12sqr/z(3))*dxyinv
          dyy(3) = -(dyz12sqr/z(3))*yinv
     &              +(z12sqr*dyz(3)/z(3)**2)*yinv
     &              - (z12sqr/z(3))*dyyinv
          dzy(3) = -(dzz12sqr/z(3))*yinv
     &              +(z12sqr*dzz(3)/z(3)**2)*yinv
     &              - (z12sqr/z(3))*dzyinv


        endif
      endif
      do j=1,3
        rot(1,j) = x(j)
        rot(2,j) = y(j)
        rot(3,j) = z(j)
        dxrot(1,j) = dxx(j)
        dxrot(2,j) = dxy(j)
        dxrot(3,j) = dxz(j)

        dyrot(1,j) = dyx(j)
        dyrot(2,j) = dyy(j)
        dyrot(3,j) = dyz(j)

        dzrot(1,j) = dzx(j)
        dzrot(2,j) = dzy(j)
        dzrot(3,j) = dzz(j)
      enddo
      return
      end
