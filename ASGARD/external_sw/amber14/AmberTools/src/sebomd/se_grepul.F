      subroutine se_grepul(iai,norbsi,ai0,ai1,ai2,di1,di2,xi,yi,zi,
     .                 iaj,norbsj,aj0,aj1,aj2,dj1,dj2,xj,yj,zj,
     .                 rij,repij,dxrepij,dyrepij,dzrepij)

      implicit none
#include "sebomd.dim"
#include "se_constants.h"
!
!     routine to calculate repulsion integrals for a pair of atoms i,j.
!
!     input:
!
!     iai,iaj     = atomic numbers for atoms i and j.
!
!     norbsi,
!     norbsj,     = numbers of atomic orbitals centered on atoms i and j.
!
!     ai0,ai1,ai2,
!     aj0,aj1,aj2 = monopole, dipole, and quadrupole additive terms.
!
!     di1,di2,
!     dj1,dj2,    = dipole and quadrupole charge separations.
!
!     xi,yi,zi,
!     xj,yj,zj,   = cartesian coordinates (angstroms).
!
!     rij         = interatomic distance (angstroms).
!
!
!     returned:
!
!     repij = 10x10 array of repulsion integrals organized according to:
!
!             s_s  s_px  s_py  s_pz  px_px  px_py  px_pz  py_py  py_pz  pz_pz
!     s_s    ----------------------------------------------------------------
!     s_px   |                  .
!     s_py   |                  .
!     s_pz   |                  .
!     px_px  |                  .
!     px_py  |................(6,4) = (px(i)*py(i)|s(j)*pz(j))
!     px_pz  |
!     py_py  |
!     py_pz  |
!     pz_pz  |
!
      integer iai, norbsi
      integer iaj, norbsj
      double precision ai0, ai1, ai2, di1, di2, xi, yi, zi
      double precision aj0, aj1, aj2, dj1, dj2, xj, yj, zj
      double precision rij, repij(10,10)
      double precision dxrepij(10,10),dyrepij(10,10),dzrepij(10,10)
!
!     local variables:
!
      double precision reploc(22),sqroot(100),dreploc(22)
      double precision dxreploc(22),dyreploc(22),dzreploc(22)
      double precision ee
      double precision eehalf
      double precision ee4th
      double precision ee8th
      double precision ee16th

      integer i, j
      double precision r
      double precision ra
      double precision rpsqr
      double precision rsqr
      double precision rp
      double precision rm
      double precision rma
      double precision rmm
      double precision rpa
      double precision rpp
      double precision rpm
      double precision rmp, rmj, rmi
      double precision rpisq, rmisq
      double precision rmsqr, rpjsq, rmjsq

      double precision a0
      double precision a00
      double precision a10
      double precision a01
      double precision a11
      double precision a20
      double precision a02
      double precision a22
      double precision a12
      double precision a21
      double precision aii2sq, ajj2sq

      double precision ri0, ri1, ri2, rpi
      double precision rj0, rj1, rj2, rpj


      double precision di22
      double precision dj22
      double precision di22sq
      double precision dj22sq
      double precision dipj, dipjsq
      double precision dimj, dimjsq
      double precision di2a
      double precision dj2a
      double precision dpsqa
      double precision dmsqa

      double precision qq,     diffqq
      double precision qdz,    diffqdz
      double precision qxxq,   diffqxxq
      double precision qqxx,   diffqqxx
      double precision qqzz,   diffqqzz
      double precision qzzq,   diffqzzq
      double precision qxyqxy, diffqxyqxy
      double precision qxzqxz, diffqxzqxz
      double precision qxxqyy, diffqxxqyy
      double precision qxxqzz, diffqxxqzz
      double precision qzzqxx, diffqzzqxx
      double precision qxxqxx, diffqxxqxx
      double precision qzzqzz, diffqzzqzz
      double precision qqqqxx, diffqqqqxx
      double precision qqqqzz, diffqqqqzz

      double precision dp
      double precision dm
      double precision dzqxx, diffdzqxx
      double precision dzqzz, diffdzqzz
      double precision qxxdz, diffqxxdz
      double precision qzzdz, diffqzzdz
      double precision qxzdx, diffqxzdx
      double precision dxqxz, diffdxqxz
      double precision dzq, diffdzq
      double precision dxdx, diffdxdx
      double precision dzdz, diffdzdz
!
!     ee = conversion factor from atomic units to ev.
!
      data ee     /27.210000d0/
      data eehalf /13.605000d0/
      data ee4th  / 6.802500d0/
      data ee8th  / 3.401250d0/
      data ee16th / 1.700625d0/
      save ee,eehalf,ee4th,ee8th,ee16th
!

      do 20 i=1,10
        do 10 j=1,10
          repij(i,j) = 0.0d0
         dxrepij(i,j) = 0.0d0
         dyrepij(i,j) = 0.0d0
         dzrepij(i,j) = 0.0d0
 10     continue
 20   continue
!
!     convert rij to atomic units.
!
      a0 = bohr2ang
      r = rij/a0
      rsqr = r*r
!
!     assign semiempirical "additive" terms that determine one-center
!     limits for two-center integrals.
!
      ri0 = 0.5d0/ai0
      rj0 = 0.5d0/aj0
      if(norbsi.gt.1)then
        ri1 = 0.5d0/ai1
        ri2 = 0.5d0/ai2
      endif
      if(norbsj.gt.1)then
        rj1 = 0.5d0/aj1
        rj2 = 0.5d0/aj2
      endif
!
!     compute the nonzero multipole-multipole interactions for a local
!     frame with xi=xj, yi=yj, and zj-zi=r, i.e.,
!
!                           z
!                           |
!                           |
!                           |
!                   atom i  |_______ y
!                          /.
!                         / .
!                        x  .
!                           .
!                           .
!                           z
!                           |
!                           |
!                           |
!                   atom j  |_______ y
!                          /
!                         /
!                        x
!
!
!     [q(i),q(i)]:
!
       a00 = (ri0 + rj0)**2
       sqroot(1) = dsqrt(rsqr + a00)
       qq = ee/sqroot(1)
       diffqq = -ee*r/sqroot(1)**3
      if(norbsi.gt.1)then
      
!
!       [dz(i),q(j)]:
!
        a10 = (ri1 + rj0)**2
        sqroot(2) = dsqrt((r+di1)**2 + a10)
        sqroot(3) = dsqrt((r-di1)**2 + a10)
        dzq = eehalf/sqroot(2) - eehalf/sqroot(3)
        diffdzq = -eehalf*(r+di1)/sqroot(2)**3 
     &           + eehalf*(r-di1)/sqroot(3)**3
!
!       [qzz(i),q(j)]:
!
        a20 = (ri2 + rj0)**2
        di22 = 2.0d0*di2
        sqroot(4) = dsqrt((r+di22)**2 + a20)
        sqroot(5) = dsqrt((r-di22)**2 + a20)
        sqroot(6) = dsqrt(rsqr + a20)
        qzzq = ee4th/sqroot(4) + ee4th/sqroot(5) - eehalf/sqroot(6)
        diffqzzq = -ee4th*(r+di22)/sqroot(4)**3
     &            - ee4th*(r-di22)/sqroot(5)**3
     &            + eehalf*r/sqroot(6)**3
!
!       [qxx(i),q(j)] or (qyy(i),q(j)]:
!
        di22sq = di22*di22
        sqroot(7) = dsqrt(rsqr + a20 + di22sq)
        sqroot(8) = sqroot(6)
        qxxq = eehalf/sqroot(7) - eehalf/sqroot(8)
        diffqxxq = -eehalf*r/sqroot(7)**3 
     &            + eehalf*r/sqroot(8)**3 
!
        if(norbsj.gt.1)then


!
!         [qxz(i),dx(j)] or [qyz(i),dy(j)]:
!
          a21 = (ri2 + rj1)**2
          rp = (r + di2)**2
          rm = (r - di2)**2
          rpa = rp + a21
          rma = rm + a21
          dp = (di2 + dj1)**2
          dm = (di2 - dj1)**2
          sqroot(9) = dsqrt(rpa + dm)
          sqroot(10) = dsqrt(rma + dm)
          sqroot(11) = dsqrt(rma + dp)
          sqroot(12) = dsqrt(rpa + dp)
          qxzdx = ee4th/sqroot(9) - ee4th/sqroot(10)
     .          + ee4th/sqroot(11) - ee4th/sqroot(12)
          diffqxzdx = -ee4th*(r + di2)/sqroot(9)**3
     &               + ee4th*(r - di2)/sqroot(10)**3
     &               - ee4th*(r - di2)/sqroot(11)**3
     &               + ee4th*(r + di2)/sqroot(12)**3
!
!         [qzz(i),dz(j)]:
!
          rp = r + dj1
          rpsqr = rp*rp
          rm = r - dj1
          rmsqr = rm*rm
          rpp = (rp + di22)**2
          rmm = (rm - di22)**2
          rpm = (rp - di22)**2
          rmp = (rm + di22)**2
          sqroot(13) = dsqrt(rmm + a21)
          sqroot(14) = dsqrt(rpp + a21)
          sqroot(15) = dsqrt(rmp + a21)
          sqroot(16) = dsqrt(rpm + a21)
          sqroot(17) = dsqrt(rmsqr + a21)
          sqroot(18) = dsqrt(rpsqr + a21)
          qzzdz = ee8th/sqroot(13) - ee8th/sqroot(14) + ee8th/sqroot(15)
     .          - ee8th/sqroot(16) - ee4th/sqroot(17) + ee4th/sqroot(18)
          diffqzzdz = -ee8th*(rm - di22)/sqroot(13)**3
     &               + ee8th*(rp + di22)/sqroot(14)**3
     &               - ee8th*(rm + di22)/sqroot(15)**3
     &               + ee8th*(rp - di22)/sqroot(16)**3
     &               + ee4th*rm/sqroot(17)**3
     &               - ee4th*rp/sqroot(18)**3
!
!         [qxx(i),dz(j)] or [qyy(i),dz(j)]:
!
          di2a = di22sq + a21
          sqroot(19) = dsqrt(rmsqr + di2a)
          sqroot(20) = dsqrt(rpsqr + di2a)
          sqroot(21) = sqroot(17)
          sqroot(22) = sqroot(18)
          qxxdz = ee4th/sqroot(19) - ee4th/sqroot(20)
     .          - ee4th/sqroot(21) + ee4th/sqroot(22)
          diffqxxdz = -ee4th*rm/sqroot(19)**3
     &               + ee4th*rp/sqroot(20)**3
     &               + ee4th*rm/sqroot(21)**3
     &               - ee4th*rp/sqroot(22)**3
         
        endif
      endif
      if(norbsj.gt.1)then
          
! 
!       [q(i),dz(j)]:
!
        if(iai.eq.iaj)then
          qdz = -dzq
          diffqdz = -diffdzq
        else
          a01 = (ri0 + rj1)**2
          sqroot(23) = dsqrt((r+dj1)**2 + a01)
          sqroot(24) = dsqrt((r-dj1)**2 + a01)
          qdz = -eehalf/sqroot(23) + eehalf/sqroot(24)
          diffqdz = eehalf*(r+dj1)/sqroot(23)**3
     &            - eehalf*(r-dj1)/sqroot(24)**3
        endif
!
!       [q(i),qzz(j)]:
!
        dj22 = 2.0d0*dj2
        if(iai.eq.iaj)then
          qqzz = qzzq
          diffqqzz = diffqzzq
        else
          a02 = (ri0 + rj2)**2
          sqroot(25) = dsqrt((r+dj22)**2 + a02)
          sqroot(26) = dsqrt((r-dj22)**2 + a02)
          sqroot(27) = dsqrt(rsqr + a02)
          qqzz = ee4th/sqroot(25) + ee4th/sqroot(26) - eehalf/sqroot(27)
          diffqqzz = -ee4th*(r+dj22)/sqroot(25)**3
     &              - ee4th*(r-dj22)/sqroot(26)**3
     &              + eehalf*r/sqroot(27)**3
        endif
!
!       [q(i),qxx(j)]:
!
        dj22sq = dj22*dj22
        if(iai.eq.iaj)then
          qqxx = qxxq
          diffqqxx = diffqxxq
        else
          sqroot(28) = dsqrt(rsqr + a02 + dj22sq)
          sqroot(29) = sqroot(27)
          qqxx = eehalf/sqroot(28) - eehalf/sqroot(29)
          diffqqxx = -eehalf*r/sqroot(28)**3
     &          + eehalf*r/sqroot(29)**3
        endif
        if(norbsi.gt.1)then
    
!
!         [dx(i),qxz(j)] or [dy(i),qyz(j)]:
!
          if(iai.eq.iaj)then
            dxqxz = -qxzdx
            diffdxqxz = -diffqxzdx
          else
            a12 = (ri1 + rj2)**2
            rp = (r + dj2)**2
            rm = (r - dj2)**2
            rpa = rp + a12
            rma = rm + a12
            dp = (dj2+di1)**2
            dm = (dj2-di1)**2
            sqroot(30) = dsqrt(rpa + dm)
            sqroot(31) = dsqrt(rma + dm)
            sqroot(32) = dsqrt(rma + dp)
            sqroot(33) = dsqrt(rpa + dp)
            dxqxz = -ee4th/sqroot(30) + ee4th/sqroot(31)
     .              -ee4th/sqroot(32) + ee4th/sqroot(33)
            diffdxqxz = ee4th*(r + dj2)/sqroot(30)**3
     &            - ee4th*(r - dj2)/sqroot(31)**3
     &            + ee4th*(r - dj2)/sqroot(32)**3
     &            - ee4th*(r + dj2)/sqroot(33)**3
          endif
!
!         [dz(i),qzz(j)]:
!
          if(iai.eq.iaj)then
            dzqzz = -qzzdz
            diffdzqzz = -diffqzzdz
          else
            rp = r + di1
            rpsqr = rp*rp
            rm = r - di1
            rmsqr = rm*rm
            rpp = (rp + dj22)**2
            rmm = (rm - dj22)**2
            rpm = (rp - dj22)**2
            rmp = (rm + dj22)**2
            sqroot(34) = dsqrt(rmm + a12)
            sqroot(35) = dsqrt(rpp + a12)
            sqroot(36) = dsqrt(rmp + a12)
            sqroot(37) = dsqrt(rpm + a12)
            sqroot(38) = dsqrt(rmsqr + a12)
            sqroot(39) = dsqrt(rpsqr + a12)
            dzqzz = -ee8th/sqroot(34) + ee8th/sqroot(35)
     .              -ee8th/sqroot(36) + ee8th/sqroot(37)
     .              +ee4th/sqroot(38) - ee4th/sqroot(39)
            diffdzqzz = ee8th*(rm - dj22)/sqroot(34)**3 
     &                - ee8th*(rp + dj22)/sqroot(35)**3
     &                + ee8th*(rm + dj22)/sqroot(36)**3
     &                - ee8th*(rp - dj22)/sqroot(37)**3
     &                - ee4th*rm/sqroot(38)**3
     &                + ee4th*rp/sqroot(39)**3
          endif
!
!         [dz(i),qxx(j)] or [dz(i),qyy(j)]:
!
          if(iai.eq.iaj)then
            dzqxx = -qxxdz
            diffdzqxx = -diffqxxdz
          else
            dj2a = dj22sq + a12
            sqroot(40) = dsqrt(rmsqr + dj2a)
            sqroot(41) = dsqrt(rpsqr + dj2a)
            sqroot(42) = sqroot(38)
            sqroot(43) = sqroot(39)
            dzqxx = -ee4th/sqroot(40) + ee4th/sqroot(41)
     .              +ee4th/sqroot(42) - ee4th/sqroot(43)
            diffdzqxx = ee4th*rm/sqroot(40)**3 
     &                - ee4th*rp/sqroot(41)**3
     &                - ee4th*rm/sqroot(42)**3 
     &                + ee4th*rp/sqroot(43)**3
          endif
!
!         [dx(i),dx(j)] or [dy(i),dy(j)]:
!
          a11 = (ri1 + rj1)**2
          ra = rsqr + a11
          dipj = di1 + dj1
          dimj = di1 - dj1
          sqroot(44) = dsqrt(ra + dimj*dimj)
          sqroot(45) = dsqrt(ra + dipj*dipj)
          dxdx = eehalf/sqroot(44) - eehalf/sqroot(45)
          diffdxdx = -eehalf*r/sqroot(44)**3 
     &              + eehalf*r/sqroot(45)**3
!
!         [dz(i),dz(j)]:
!
          sqroot(46) = dsqrt((r+dimj)**2 + a11)
          sqroot(47) = dsqrt((r+dipj)**2 + a11)
          sqroot(48) = dsqrt((r-dipj)**2 + a11)
          sqroot(49) = dsqrt((r-dimj)**2 + a11)
          dzdz = ee4th/sqroot(46) - ee4th/sqroot(47)
     .         - ee4th/sqroot(48) + ee4th/sqroot(49)
          diffdzdz = -ee4th*(r+dimj)/sqroot(46)**3 
     &              + ee4th*(r+dipj)/sqroot(47)**3 
     &              + ee4th*(r-dipj)/sqroot(48)**3 
     &              - ee4th*(r-dimj)/sqroot(49)**3 


!
!         [qzz(i),qzz(j)]:
!
          a22 = (ri2 + rj2)**2
          ra = rsqr + a22
          rpi = r + di22
          rmi = r - di22
          rpj = r + dj22
          rmj = r - dj22
          sqroot(50) = dsqrt((rpi - dj22)**2 + a22)
          sqroot(51) = dsqrt((rpi + dj22)**2 + a22)
          sqroot(52) = dsqrt((rmi - dj22)**2 + a22)
          sqroot(53) = dsqrt((rmi + dj22)**2 + a22)
          sqroot(54) = dsqrt(rpi*rpi + a22)
          sqroot(55) = dsqrt(rmi*rmi + a22)
          sqroot(56) = dsqrt(rpj*rpj + a22)
          sqroot(57) = dsqrt(rmj*rmj + a22)
          sqroot(58) = dsqrt(ra)
          qzzqzz = ee16th/sqroot(50) + ee16th/sqroot(51)
     .           + ee16th/sqroot(52) + ee16th/sqroot(53)
     .           - ee8th/sqroot(54) - ee8th/sqroot(55)
     .           - ee8th/sqroot(56) - ee8th/sqroot(57)
     .           + ee4th/sqroot(58)
         
          diffqzzqzz = -ee16th*(rpi - dj22)/sqroot(50)**3
     &            - ee16th*(rpi + dj22)/sqroot(51)**3
     &            - ee16th*(rmi - dj22)/sqroot(52)**3 
     &            - ee16th*(rmi + dj22)/sqroot(53)**3
     &            + ee8th*rpi/sqroot(54)**3
     &            + ee8th*rmi/sqroot(55)**3
     &            + ee8th*rpj/sqroot(56)**3
     &            + ee8th*rmj/sqroot(57)**3
     &            - ee4th*r/sqroot(58)**3
!
!         [qxx(i),qxx(j)] or [qyy(i),qyy(i)]:
!
          sqroot(59) = dsqrt(ra + (di22-dj22)**2)
          sqroot(60) = dsqrt(ra + (di22+dj22)**2)
          sqroot(61) = dsqrt(ra + di22sq)
          sqroot(62) = dsqrt(ra + dj22sq)
          sqroot(63) = sqroot(58)
          qxxqxx = ee8th/sqroot(59) + ee8th/sqroot(60)
     .           - ee4th/sqroot(61) - ee4th/sqroot(62)
     .           + ee4th/sqroot(63)
          diffqxxqxx = -ee8th*r/sqroot(59)**3 
     &                - ee8th*r/sqroot(60)**3
     &                + ee4th*r/sqroot(61)**3 
     &                + ee4th*r/sqroot(62)**3
     &                - ee4th*r/sqroot(63)**3
!
!         [qzz(i),qxx(j)] or [qzz(i),qyy(i)]:
!
          ajj2sq = a22 + dj22sq
          rpisq = rpi*rpi
          rmisq = rmi*rmi
          sqroot(64) = dsqrt(rmisq + ajj2sq)
          sqroot(65) = dsqrt(rpisq + ajj2sq)
          sqroot(66) = sqroot(55)
          sqroot(67) = sqroot(54)
          sqroot(68) = sqroot(62)
          sqroot(69) = sqroot(58)
          qzzqxx = ee8th/sqroot(64) + ee8th/sqroot(65)
     .           - ee8th/sqroot(66) - ee8th/sqroot(67)
     .           - ee4th/sqroot(68) + ee4th/sqroot(69)
          diffqzzqxx = -ee8th*rmi/sqroot(64)**3 
     &                - ee8th*rpi/sqroot(65)**3
     &                + ee8th*rmi/sqroot(66)**3
     &                + ee8th*rpi/sqroot(67)**3
     &                + ee4th*r/sqroot(68)**3
     &                - ee4th*r/sqroot(69)**3
!
!         [qxx(i),qzz(j)] or [qyy(i),qzz(j)]
!
          if(iai.eq.iaj)then
            qxxqzz = qzzqxx
            diffqxxqzz = diffqzzqxx
          else
            aii2sq = a22 + di22sq
            rpjsq = rpj*rpj
            rmjsq = rmj*rmj
            sqroot(70) = dsqrt(rmjsq + aii2sq)
            sqroot(71) = dsqrt(rpjsq + aii2sq)
            sqroot(72) = sqroot(57)
            sqroot(73) = sqroot(56)
            sqroot(74) = sqroot(61)
            sqroot(75) = sqroot(58)
            qxxqzz = ee8th/sqroot(70) + ee8th/sqroot(71)
     .             - ee8th/sqroot(72) - ee8th/sqroot(73)
     .             - ee4th/sqroot(74) + ee4th/sqroot(75)
            diffqxxqzz = -ee8th*rmj/sqroot(70)**3 
     &                  - ee8th*rpj/sqroot(71)**3 
     &                  + ee8th*rmj/sqroot(72)**3  
     &                  + ee8th*rpj/sqroot(73)**3 
     &                  + ee4th*r/sqroot(74)**3  
     &                  - ee4th*r/sqroot(75)**3 
          endif
!
!         [qxx(i),qyy(j)] or [qyy(i),qxx(j)]:
!
          sqroot(76) = dsqrt(ra + di22sq + dj22sq)
          sqroot(77) = sqroot(61)
          sqroot(78) = sqroot(62)
          sqroot(79) = sqroot(58)
          qxxqyy = ee4th/sqroot(76) - ee4th/sqroot(77)
     .           - ee4th/sqroot(78) + ee4th/sqroot(79)
          diffqxxqyy = -ee4th*r/sqroot(76)**3
     &                + ee4th*r/sqroot(77)**3
     &                + ee4th*r/sqroot(78)**3
     &                - ee4th*r/sqroot(79)**3
!
!         [qxz(i),qxz(j)] or [qyz(i),qyz(j)]:
!

          dipj = di2 + dj2
          dimj = di2 - dj2
          dipjsq = dipj*dipj
          dimjsq = dimj*dimj
          rpp = (r+dipj)**2
          rmp = (r-dipj)**2
          rpm = (r+dimj)**2
          rmm = (r-dimj)**2
          dpsqa = dipjsq + a22
          dmsqa = dimjsq + a22
          sqroot(80) = dsqrt(rpm + dmsqa)
          sqroot(81) = dsqrt(rpm + dpsqa)
          sqroot(82) = dsqrt(rpp + dmsqa)
          sqroot(83) = dsqrt(rpp + dpsqa)
          sqroot(84) = dsqrt(rmp + dmsqa)
          sqroot(85) = dsqrt(rmp + dpsqa)
          sqroot(86) = dsqrt(rmm + dmsqa)
          sqroot(87) = dsqrt(rmm + dpsqa)
          qxzqxz = ee8th/sqroot(80) - ee8th/sqroot(81)
     .         - ee8th/sqroot(82) + ee8th/sqroot(83)
     .         - ee8th/sqroot(84) + ee8th/sqroot(85)
     .         + ee8th/sqroot(86) - ee8th/sqroot(87)
          diffqxzqxz = -ee8th*(r+dimj)/sqroot(80)**3
     &                + ee8th*(r+dimj)/sqroot(81)**3
     &                + ee8th*(r+dipj)/sqroot(82)**3 
     &                - ee8th*(r+dipj)/sqroot(83)**3
     &                + ee8th*(r-dipj)/sqroot(84)**3
     &                - ee8th*(r-dipj)/sqroot(85)**3
     &                - ee8th*(r-dimj)/sqroot(86)**3
     &                + ee8th*(r-dimj)/sqroot(87)**3
!
!     [qxy(i),qxy(j)] or  [qxx-yy(i),qxx-yy(j)]:
!
          qxyqxy = 0.5d0*(qxxqxx - qxxqyy)
          diffqxyqxy = 0.5d0*(diffqxxqxx - diffqxxqyy)

        endif
      endif
!
!     now use the local multipole-multipole interactions to construct
!     the 22 nonzero local repulsion integrals.
!
!     (s s | s s):
! 
    

      reploc(1) = qq
      dreploc(1) = diffqq
      if(norbsi.gt.1)then
!
!       (s pz | s s):
!
        reploc(2) = dzq
        dreploc(2) = diffdzq
!
!       (pz pz | s s):
!
        reploc(3) = qzzq + qq
        dreploc(3) = diffqzzq + diffqq
!
!       (px px | s s) or (py py | s s):
!
        reploc(4) = qxxq + qq
        dreploc(4) = diffqxxq + diffqq
      endif
      if(norbsj.gt.1)then
!
!       (s s | s pz):
!
        reploc(5) = qdz
        dreploc(5) = diffqdz
        if(norbsi.gt.1)then
!
!         (s pz | s pz):
!
          reploc(6) = dzdz
          dreploc(6) = diffdzdz
!
!         (s px | s px) or (s py | s py):
!
          reploc(7) = dxdx
          dreploc(7) = diffdxdx
!
!         (pz pz | s pz):
!
          reploc(8) = qzzdz + qdz
          dreploc(8) = diffqzzdz + diffqdz
!
!         (px px | s pz) or (py py | s pz):
!
          reploc(9) = qxxdz + qdz
          dreploc(9) = diffqxxdz + diffqdz
!
!         (px pz | s px):
!
          reploc(10) = qxzdx
          dreploc(10) = diffqxzdx
        endif
!
!       (s s | pz pz):
!
        reploc(11) = qqzz + qq
        dreploc(11) = diffqqzz + diffqq
!
!       (s s | px px) or (s s | py py):
!
        reploc(12) = qqxx + qq
        dreploc(12) = diffqqxx + diffqq
        if(norbsi.gt.1)then
!
!         (s pz | pz pz):
!
          reploc(13) = dzqzz + dzq
          dreploc(13) = diffdzqzz + diffdzq
!
!         (s pz | px px) or (s pz | py py):
!
          reploc(14) = dzqxx + dzq
          dreploc(14) = diffdzqxx + diffdzq
!
!         (s px | px pz) or (s py | py pz):
!
          reploc(15) = dxqxz
          dreploc(15) = diffdxqxz
          qqqqzz = qq + qqzz
          diffqqqqzz = diffqq + diffqqzz
!
!         (pz pz | pz pz):
!
          reploc(16) = qqqqzz + qzzq + qzzqzz
          dreploc(16) = diffqqqqzz + diffqzzq + diffqzzqzz
!
!         (px px | pz pz) or (py py | pz pz):
!
          reploc(17) = qqqqzz + qxxq + qxxqzz
          dreploc(17) = diffqqqqzz + diffqxxq + diffqxxqzz
          qqqqxx = qq + qqxx
          diffqqqqxx = diffqq + diffqqxx
!
!         (pz pz | px px) or (pz pz | py py):
!
          reploc(18) = qqqqxx + qzzq + qzzqxx
          dreploc(18) = diffqqqqxx + diffqzzq + diffqzzqxx
!
!         (px px | px px) or (py py | py py):
!
          reploc(19) = qqqqxx + qxxq + qxxqxx
          dreploc(19) = diffqqqqxx + diffqxxq + diffqxxqxx
!
!         (px pz | px pz) or (py pz | py pz):
!
          reploc(20) = qxzqxz
          dreploc(20) = diffqxzqxz
!
!         (px px | py py) or (py py | px px):
!
          reploc(21) = qqqqxx + qxxq +  qxxqyy
          dreploc(21) = diffqqqqxx + diffqxxq +  diffqxxqyy
!
!         (px py | px py):
!
          reploc(22) = qxyqxy
          dreploc(22) = diffqxyqxy
        endif
      endif
!
!       d        d
!       --  -->  --  
!       dr       dx
!
         do i=1,22
          dxreploc(i)=dreploc(i)*(xi-xj)/rij/bohr2ang
          dyreploc(i)=dreploc(i)*(yi-yj)/rij/bohr2ang
          dzreploc(i)=dreploc(i)*(zi-zj)/rij/bohr2ang
         enddo

!
!     rotate local repulsion integrals to molecular frame.
!
      if(norbsi.eq.1.and.norbsj.eq.1)then
        repij(1,1) = reploc(1)
         dxrepij(1,1) = dreploc(1)*(xi-xj)/rij/bohr2ang
         dyrepij(1,1) = dreploc(1)*(yi-yj)/rij/bohr2ang
         dzrepij(1,1) = dreploc(1)*(zi-zj)/rij/bohr2ang
      else
        call se_grrep(norbsi,norbsj,xi,yi,zi,xj,yj,zj,rij,reploc,repij,
     &            dxrepij,dyrepij,dzrepij,dxreploc,dyreploc,dzreploc)
      endif

      return
      end
!--------------------------------------------------------------------------------
      subroutine se_grrep(norbsi,norbsj,xi,yi,zi,xj,yj,zj,rij,reploc,
     &         repij,dxrepij,dyrepij,dzrepij,dxreploc,dyreploc,dzreploc)

      implicit none
      integer norbsi
      integer norbsj
      double precision xi, yi, zi
      double precision xj, yj, zj
      double precision rij
      double precision reploc(22)
      double precision repij(10,10)
      double precision dxreploc(22)
      double precision dxrepij(10,10)
      double precision dyreploc(22)
      double precision dyrepij(10,10)
      double precision dzreploc(22)
      double precision dzrepij(10,10)

      double precision z12, z12sqr, yinv
      double precision dxz12, dxz12sqr, dxyinv
      double precision dyz12, dyz12sqr, dyyinv
      double precision dzz12, dzz12sqr, dzyinv

      double precision r11r11, r11r12, r11r13 ! 1-1
      double precision         r12r12, r12r13
      double precision                 r13r13

      double precision dxr11r11, dxr11r12, dxr11r13
      double precision           dxr12r12, dxr12r13
      double precision                     dxr13r13

      double precision dyr11r11, dyr11r12, dyr11r13
      double precision           dyr12r12, dyr12r13
      double precision                     dyr13r13

      double precision dzr11r11, dzr11r12, dzr11r13
      double precision           dzr12r12, dzr12r13
      double precision                     dzr13r13

      double precision r11r21, r11r22, r11r23 ! 1-2
      double precision         r12r22, r12r23
      double precision                 r13r23

      double precision dxr11r21, dxr11r22, dxr11r23
      double precision           dxr12r22, dxr12r23
      double precision                     dxr13r23

      double precision dyr11r21, dyr11r22, dyr11r23
      double precision           dyr12r22, dyr12r23
      double precision                     dyr13r23

      double precision dzr11r21, dzr11r22, dzr11r23
      double precision           dzr12r22, dzr12r23
      double precision                     dzr13r23

      double precision r11r31, r11r32, r11r33 ! 1-3
      double precision r12r31, r12r32, r12r33
      double precision                 r13r33

      double precision dxr11r31, dxr11r32, dxr11r33
      double precision dxr12r31, dxr12r32, dxr12r33
      double precision                     dxr13r33

      double precision dyr11r31, dyr11r32, dyr11r33
      double precision dyr12r31, dyr12r32, dyr12r33
      double precision                     dyr13r33

      double precision dzr11r31, dzr11r32, dzr11r33
      double precision dzr12r31, dzr12r32, dzr12r33
      double precision                     dzr13r33

      double precision r21r11, r21r12, r21r13 ! 2-1
      double precision         r22r12, r22r13
      double precision                 r23r13

      double precision dxr21r11, dxr21r12, dxr21r13
      double precision           dxr22r12, dxr22r13
      double precision                     dxr23r13

      double precision dyr21r11, dyr21r12, dyr21r13
      double precision           dyr22r12, dyr22r13
      double precision                     dyr23r13

      double precision dzr21r11, dzr21r12, dzr21r13
      double precision           dzr22r12, dzr22r13
      double precision                     dzr23r13

      double precision r21r21, r21r22, r21r23 ! 2-2
      double precision         r22r22, r22r23
      double precision                 r23r23

      double precision dxr21r21, dxr21r22, dxr21r23
      double precision           dxr22r22, dxr22r23
      double precision                     dxr23r23

      double precision dyr21r21, dyr21r22, dyr21r23
      double precision           dyr22r22, dyr22r23
      double precision                     dyr23r23

      double precision dzr21r21, dzr21r22, dzr21r23
      double precision           dzr22r22, dzr22r23
      double precision                     dzr23r23

      double precision r21r31, r21r32, r21r33 ! 2-3
      double precision r22r31, r22r32, r22r33
      double precision                 r23r33

      double precision dxr21r31, dxr21r32, dxr21r33
      double precision dxr22r31, dxr22r32, dxr22r33
      double precision                     dxr23r33

      double precision dyr21r31, dyr21r32, dyr21r33
      double precision dyr22r31, dyr22r32, dyr22r33
      double precision                     dyr23r33

      double precision dzr21r31, dzr21r32, dzr21r33
      double precision dzr22r31, dzr22r32, dzr22r33
      double precision                     dzr23r33

      double precision r31r11, r31r12, r31r13 ! 3-1
      double precision         r32r12, r32r13
      double precision                 r33r13

      double precision dxr31r11, dxr31r12, dxr31r13
      double precision           dxr32r12, dxr32r13
      double precision                     dxr33r13

      double precision dyr31r11, dyr31r12, dyr31r13
      double precision           dyr32r12, dyr32r13
      double precision                     dyr33r13

      double precision dzr31r11, dzr31r12, dzr31r13
      double precision           dzr32r12, dzr32r13
      double precision                     dzr33r13

      double precision r31r21, r31r22, r31r23 ! 3-2
      double precision         r32r22, r32r23
      double precision                 r33r23

      double precision dxr31r21, dxr31r22, dxr31r23
      double precision           dxr32r22, dxr32r23
      double precision                     dxr33r23

      double precision dyr31r21, dyr31r22, dyr31r23
      double precision           dyr32r22, dyr32r23
      double precision                     dyr33r23

      double precision dzr31r21, dzr31r22, dzr31r23
      double precision           dzr32r22, dzr32r23
      double precision                     dzr33r23

      double precision r31r31, r31r32, r31r33 ! 3-3
      double precision         r32r32, r32r33
      double precision                 r33r33

      double precision dxr31r31, dxr31r32, dxr31r33
      double precision           dxr32r32, dxr32r33
      double precision                     dxr33r33

      double precision dyr31r31, dyr31r32, dyr31r33
      double precision           dyr32r32, dyr32r33
      double precision                     dyr33r33

      double precision dzr31r31, dzr31r32, dzr31r33
      double precision           dzr32r32, dzr32r33
      double precision                     dzr33r33

      double precision radd01, dxradd01, dyradd01, dzradd01
      double precision radd02, dxradd02, dyradd02, dzradd02
      double precision radd03, dxradd03, dyradd03, dzradd03
      double precision radd04, dxradd04, dyradd04, dzradd04
      double precision radd05, dxradd05, dyradd05, dzradd05
      double precision radd06, dxradd06, dyradd06, dzradd06
      double precision radd07, dxradd07, dyradd07, dzradd07
      double precision radd08, dxradd08, dyradd08, dzradd08
      double precision radd09, dxradd09, dyradd09, dzradd09
      double precision radd10, dxradd10, dyradd10, dzradd10
      double precision radd11, dxradd11, dyradd11, dzradd11
      double precision radd12, dxradd12, dyradd12, dzradd12
      double precision radd13, dxradd13, dyradd13, dzradd13
      double precision radd14, dxradd14, dyradd14, dzradd14
      double precision radd15, dxradd15, dyradd15, dzradd15
      double precision radd16, dxradd16, dyradd16, dzradd16
      double precision radd17, dxradd17, dyradd17, dzradd17
!
!     rotates the local frame repulsion integrals to the molecular frame.
!
      double precision r(3,3),dxr(3,3),dyr(3,3),dzr(3,3)
!
!     generate rotation matrix.
!
      r(3,1) = (xi-xj)/rij
       dxr(3,1) = (1/rij - (xi-xj)*(xi-xj)/rij**3)
       dyr(3,1) = -(xi-xj)*(yi-yj)/rij**3
       dzr(3,1) = -(xi-xj)*(zi-zj)/rij**3
      r(3,2) = (yi-yj)/rij
       dxr(3,2) = -(yi-yj)*(xi-xj)/rij**3
       dyr(3,2) = (1/rij - (yi-yj)*(yi-yj)/rij**3)
       dzr(3,2) = -(yi-yj)*(zi-zj)/rij**3
      r(3,3) = (zi-zj)/rij
       dxr(3,3) = -(zi-zj)*(xi-xj)/rij**3
       dyr(3,3) = -(zi-zj)*(yi-yj)/rij**3
       dzr(3,3) = (1/rij - (zi-zj)*(zi-zj)/rij**3)
      if(abs(r(3,3)).gt.0.99999999d0)then
        r(3,1) = 0.0d0
         dxr(3,1) = 0.0d0
         dyr(3,1) = 0.0d0
         dzr(3,1) = 0.0d0
        r(3,2) = 0.0d0
         dxr(3,2) = 0.0d0
         dyr(3,2) = 0.0d0
         dzr(3,2) = 0.0d0
        r(3,3) = sign(1.0d0,r(3,3))
         dxr(3,3) = 0.0d0
         dyr(3,3) = 0.0d0
         dzr(3,3) = 0.0d0
        r(1,1) = 1.0d0
         dxr(1,1) = 0.0d0
         dyr(1,1) = 0.0d0
         dzr(1,1) = 0.0d0
        r(1,2) = 0.0d0
         dxr(1,2) = 0.0d0
         dyr(1,2) = 0.0d0
         dzr(1,2) = 0.0d0
        r(1,3) = 0.0d0
         dxr(1,3) = 0.0d0
         dyr(1,3) = 0.0d0
         dzr(1,3) = 0.0d0
        r(2,1) = 0.0d0
         dxr(2,1) = 0.0d0
         dyr(2,1) = 0.0d0
         dzr(2,1) = 0.0d0
        r(2,2) = r(3,3)
         dxr(2,2) = 0.0d0
         dyr(2,2) = 0.0d0
         dzr(2,2) = 0.0d0      
        r(2,3) = 0.0d0
         dxr(2,3) = 0.0d0
         dyr(2,3) = 0.0d0
         dzr(2,3) = 0.0d0
      else
        z12sqr = r(3,1)**2 + r(3,2)**2
         dxz12sqr = 2.0d0*(r(3,1)*dxr(3,1)+r(3,2)*dxr(3,2))
         dyz12sqr = 2.0d0*(r(3,1)*dyr(3,1)+r(3,2)*dyr(3,2))
         dzz12sqr = 2.0d0*(r(3,1)*dzr(3,1)+r(3,2)*dzr(3,2))
        z12 = dsqrt(z12sqr)
         dxz12 = dxz12sqr/(2.0d0*z12)
         dyz12 = dyz12sqr/(2.0d0*z12)
         dzz12 = dzz12sqr/(2.0d0*z12)
        r(1,1) = -r(3,2)/z12
         dxr(1,1) = (-dxr(3,2)/z12 + r(3,2)*dxz12/z12sqr)
         dyr(1,1) = (-dyr(3,2)/z12 + r(3,2)*dyz12/z12sqr)
         dzr(1,1) = (-dzr(3,2)/z12 + r(3,2)*dzz12/z12sqr)
        r(1,2) = r(3,1)/z12
         dxr(1,2) = (dxr(3,1)/z12 - r(3,1)*dxz12/z12sqr)
         dyr(1,2) = (dyr(3,1)/z12 - r(3,1)*dyz12/z12sqr)
         dzr(1,2) = (dzr(3,1)/z12 - r(3,1)*dzz12/z12sqr)
        r(1,3) = 0.0d0
         dxr(1,3) = 0.0d0
         dyr(1,3) = 0.0d0
         dzr(1,3) = 0.0d0
        if(abs(r(3,3)).lt.1.0d-8)then
          r(2,1) = 0.0d0
           dxr(2,1) = 0.0d0
           dyr(2,1) = 0.0d0
           dzr(2,1) = 0.0d0
          r(2,2) = 0.0d0
           dxr(2,2) = 0.0d0
           dyr(2,2) = 0.0d0
           dzr(2,2) = 0.0d0
          r(2,3) = 1.0d0
           dxr(2,3) = 0.0d0
           dyr(2,3) = 0.0d0
           dzr(2,3) = 0.0d0
        else
          yinv = 1.0d0/(z12*dsqrt(1.0d0 + z12sqr/r(3,3)**2))
           dxyinv = -(dxz12*dsqrt(1.0d0 + z12sqr/r(3,3)**2)
     &                     +z12*0.5d0*(-2.0d0*z12sqr*dxr(3,3)/r(3,3)**3
     &           + dxz12sqr/r(3,3)**2)/dsqrt(1.0d0 + z12sqr/r(3,3)**2))
     &             /(z12*dsqrt(1.0d0 + z12sqr/r(3,3)**2))**2

           dyyinv = -(dyz12*dsqrt(1.0d0 + z12sqr/r(3,3)**2)
     &                     +z12*0.5d0*(-2.0d0*z12sqr*dyr(3,3)/r(3,3)**3
     &           + dyz12sqr/r(3,3)**2)/dsqrt(1.0d0 + z12sqr/r(3,3)**2))
     &             /(z12*dsqrt(1.0d0 + z12sqr/r(3,3)**2))**2


           dzyinv = -(dzz12*dsqrt(1.0d0 + z12sqr/r(3,3)**2)
     &                     +z12*0.5d0*(-2.0d0*z12sqr*dzr(3,3)/r(3,3)**3
     &           + dzz12sqr/r(3,3)**2)/dsqrt(1.0d0 + z12sqr/r(3,3)**2))
     &             /(z12*dsqrt(1.0d0 + z12sqr/r(3,3)**2))**2



          r(2,1) = r(3,1)*yinv
           dxr(2,1) = (dxr(3,1)*yinv + r(3,1)*dxyinv)
           dyr(2,1) = (dyr(3,1)*yinv + r(3,1)*dyyinv)
           dzr(2,1) = (dzr(3,1)*yinv + r(3,1)*dzyinv)       
          r(2,2) = r(3,2)*yinv
           dxr(2,2) = (dxr(3,2)*yinv + r(3,2)*dxyinv)
           dyr(2,2) = (dyr(3,2)*yinv + r(3,2)*dyyinv)
           dzr(2,2) = (dzr(3,2)*yinv + r(3,2)*dzyinv)
          r(2,3) = -(z12sqr/r(3,3))*yinv
           dxr(2,3) = (-(dxz12sqr/r(3,3))*yinv 
     &              + (z12sqr*dxr(3,3)/r(3,3)**2)*yinv 
     &              - (z12sqr/r(3,3))*dxyinv)
           dyr(2,3) = (-(dyz12sqr/r(3,3))*yinv 
     &              + (z12sqr*dyr(3,3)/r(3,3)**2)*yinv 
     &              - (z12sqr/r(3,3))*dyyinv)
           dzr(2,3) = (-(dzz12sqr/r(3,3))*yinv 
     &              + (z12sqr*dzr(3,3)/r(3,3)**2)*yinv 
     &              - (z12sqr/r(3,3))*dzyinv)
        endif
      endif

      r11r11 = r(1,1)*r(1,1)
       dxr11r11 = dxr(1,1)*r(1,1) + r(1,1)*dxr(1,1)
       dyr11r11 = dyr(1,1)*r(1,1) + r(1,1)*dyr(1,1)
       dzr11r11 = dzr(1,1)*r(1,1) + r(1,1)*dzr(1,1)
      r11r12 = r(1,1)*r(1,2)
       dxr11r12 = dxr(1,1)*r(1,2) + r(1,1)*dxr(1,2)
       dyr11r12 = dyr(1,1)*r(1,2) + r(1,1)*dyr(1,2)
       dzr11r12 = dzr(1,1)*r(1,2) + r(1,1)*dzr(1,2)
      r11r13 = r(1,1)*r(1,3)
       dxr11r13 = dxr(1,1)*r(1,3) + r(1,1)*dxr(1,3)
       dyr11r13 = dyr(1,1)*r(1,3) + r(1,1)*dyr(1,3)
       dzr11r13 = dzr(1,1)*r(1,3) + r(1,1)*dzr(1,3)
      r12r12 = r(1,2)*r(1,2)
       dxr12r12 = dxr(1,2)*r(1,2) + r(1,2)*dxr(1,2)
       dyr12r12 = dyr(1,2)*r(1,2) + r(1,2)*dyr(1,2)
       dzr12r12 = dzr(1,2)*r(1,2) + r(1,2)*dzr(1,2)
      r12r13 = r(1,2)*r(1,3) 
       dxr12r13 = dxr(1,2)*r(1,3) + r(1,2)*dxr(1,3)
       dyr12r13 = dyr(1,2)*r(1,3) + r(1,2)*dyr(1,3)
       dzr12r13 = dzr(1,2)*r(1,3) + r(1,2)*dzr(1,3)
      r13r13 = r(1,3)*r(1,3) 
       dxr13r13 = dxr(1,3)*r(1,3) + r(1,3)*dxr(1,3)
       dyr13r13 = dyr(1,3)*r(1,3) + r(1,3)*dyr(1,3)
       dzr13r13 = dzr(1,3)*r(1,3) + r(1,3)*dzr(1,3)
      r21r21 = r(2,1)*r(2,1)
       dxr21r21 = dxr(2,1)*r(2,1) + r(2,1)*dxr(2,1)
       dyr21r21 = dyr(2,1)*r(2,1) + r(2,1)*dyr(2,1)
       dzr21r21 = dzr(2,1)*r(2,1) + r(2,1)*dzr(2,1)
      r21r22 = r(2,1)*r(2,2)
       dxr21r22 = dxr(2,1)*r(2,2) + r(2,1)*dxr(2,2)
       dyr21r22 = dyr(2,1)*r(2,2) + r(2,1)*dyr(2,2)
       dzr21r22 = dzr(2,1)*r(2,2) + r(2,1)*dzr(2,2)
      r21r23 = r(2,1)*r(2,3)
       dxr21r23 = dxr(2,1)*r(2,3) + r(2,1)*dxr(2,3)
       dyr21r23 = dyr(2,1)*r(2,3) + r(2,1)*dyr(2,3)
       dzr21r23 = dzr(2,1)*r(2,3) + r(2,1)*dzr(2,3)
      r22r22 = r(2,2)*r(2,2) 
       dxr22r22 = dxr(2,2)*r(2,2) + r(2,2)*dxr(2,2)
       dyr22r22 = dyr(2,2)*r(2,2) + r(2,2)*dyr(2,2)
       dzr22r22 = dzr(2,2)*r(2,2) + r(2,2)*dzr(2,2)
      r22r23 = r(2,2)*r(2,3) 
       dxr22r23 = dxr(2,2)*r(2,3) + r(2,2)*dxr(2,3)
       dyr22r23 = dyr(2,2)*r(2,3) + r(2,2)*dyr(2,3)
       dzr22r23 = dzr(2,2)*r(2,3) + r(2,2)*dzr(2,3)
      r23r23 = r(2,3)*r(2,3)
       dxr23r23 = dxr(2,3)*r(2,3) + r(2,3)*dxr(2,3)
       dyr23r23 = dyr(2,3)*r(2,3) + r(2,3)*dyr(2,3)
       dzr23r23 = dzr(2,3)*r(2,3) + r(2,3)*dzr(2,3)  
      r31r31 = r(3,1)*r(3,1)
       dxr31r31 = dxr(3,1)*r(3,1) + r(3,1)*dxr(3,1)
       dyr31r31 = dyr(3,1)*r(3,1) + r(3,1)*dyr(3,1)
       dzr31r31 = dzr(3,1)*r(3,1) + r(3,1)*dzr(3,1)
      r31r32 = r(3,1)*r(3,2) 
       dxr31r32 = dxr(3,1)*r(3,2) + r(3,1)*dxr(3,2)
       dyr31r32 = dyr(3,1)*r(3,2) + r(3,1)*dyr(3,2)
       dzr31r32 = dzr(3,1)*r(3,2) + r(3,1)*dzr(3,2)
      r31r33 = r(3,1)*r(3,3)  
       dxr31r33 = dxr(3,1)*r(3,3) + r(3,1)*dxr(3,3)
       dyr31r33 = dyr(3,1)*r(3,3) + r(3,1)*dyr(3,3)
       dzr31r33 = dzr(3,1)*r(3,3) + r(3,1)*dzr(3,3)
      r32r32 = r(3,2)*r(3,2)
       dxr32r32 = dxr(3,2)*r(3,2) + r(3,2)*dxr(3,2)
       dyr32r32 = dyr(3,2)*r(3,2) + r(3,2)*dyr(3,2)
       dzr32r32 = dzr(3,2)*r(3,2) + r(3,2)*dzr(3,2)
      r32r33 = r(3,2)*r(3,3)
       dxr32r33 = dxr(3,2)*r(3,3) + r(3,2)*dxr(3,3)
       dyr32r33 = dyr(3,2)*r(3,3) + r(3,2)*dyr(3,3)
       dzr32r33 = dzr(3,2)*r(3,3) + r(3,2)*dzr(3,3)
      r33r33 = r(3,3)*r(3,3)
       dxr33r33 = dxr(3,3)*r(3,3) + r(3,3)*dxr(3,3)
       dyr33r33 = dyr(3,3)*r(3,3) + r(3,3)*dyr(3,3)
       dzr33r33 = dzr(3,3)*r(3,3) + r(3,3)*dzr(3,3)
      if(norbsi.gt.1.or.norbsj.gt.1)then
        radd01 = r11r11+r21r21
         dxradd01 = dxr11r11+dxr21r21
         dyradd01 = dyr11r11+dyr21r21
         dzradd01 = dzr11r11+dzr21r21
        radd02 = r11r12+r21r22
         dxradd02 = dxr11r12+dxr21r22
         dyradd02 = dyr11r12+dyr21r22
         dzradd02 = dzr11r12+dzr21r22
        radd03 = r11r13+r21r23
         dxradd03 = dxr11r13+dxr21r23
         dyradd03 = dyr11r13+dyr21r23
         dzradd03 = dzr11r13+dzr21r23
        radd04 = r12r12+r22r22
         dxradd04 = dxr12r12+dxr22r22
         dyradd04 = dyr12r12+dyr22r22
         dzradd04 = dzr12r12+dzr22r22
        radd05 = r12r13+r22r23
         dxradd05 = dxr12r13+dxr22r23
         dyradd05 = dyr12r13+dyr22r23
         dzradd05 = dzr12r13+dzr22r23
        radd06 = r13r13+r23r23
         dxradd06 = dxr13r13+dxr23r23
         dyradd06 = dyr13r13+dyr23r23
         dzradd06 = dzr13r13+dzr23r23
      endif
      if(norbsi.gt.1.and.norbsj.gt.1)then
        r11r21 = r(1,1)*r(2,1)*2.0d0
         dxr11r21 = (dxr(1,1)*r(2,1) + r(1,1)*dxr(2,1))*2.0d0
         dyr11r21 = (dyr(1,1)*r(2,1) + r(1,1)*dyr(2,1))*2.0d0
         dzr11r21 = (dzr(1,1)*r(2,1) + r(1,1)*dzr(2,1))*2.0d0
        r11r22 = r(1,1)*r(2,2)
         dxr11r22 = dxr(1,1)*r(2,2) + r(1,1)*dxr(2,2)
         dyr11r22 = dyr(1,1)*r(2,2) + r(1,1)*dyr(2,2)
         dzr11r22 = dzr(1,1)*r(2,2) + r(1,1)*dzr(2,2)
        r11r23 = r(1,1)*r(2,3)
         dxr11r23 = dxr(1,1)*r(2,3) + r(1,1)*dxr(2,3)
         dyr11r23 = dyr(1,1)*r(2,3) + r(1,1)*dyr(2,3)
         dzr11r23 = dzr(1,1)*r(2,3) + r(1,1)*dzr(2,3)
        r11r31 = r(1,1)*r(3,1)*2.0d0
         dxr11r31 = (dxr(1,1)*r(3,1) + r(1,1)*dxr(3,1))*2.0d0
         dyr11r31 = (dyr(1,1)*r(3,1) + r(1,1)*dyr(3,1))*2.0d0
         dzr11r31 = (dzr(1,1)*r(3,1) + r(1,1)*dzr(3,1))*2.0d0
        r11r32 = r(1,1)*r(3,2)
         dxr11r32 = dxr(1,1)*r(3,2) + r(1,1)*dxr(3,2)
         dyr11r32 = dyr(1,1)*r(3,2) + r(1,1)*dyr(3,2)
         dzr11r32 = dzr(1,1)*r(3,2) + r(1,1)*dzr(3,2)
        r11r33 = r(1,1)*r(3,3)
         dxr11r33 = dxr(1,1)*r(3,3) + r(1,1)*dxr(3,3)
         dyr11r33 = dyr(1,1)*r(3,3) + r(1,1)*dyr(3,3)
         dzr11r33 = dzr(1,1)*r(3,3) + r(1,1)*dzr(3,3)
        r12r22 = r(1,2)*r(2,2)*2.0d0
         dxr12r22 = (dxr(1,2)*r(2,2) + r(1,2)*dxr(2,2))*2.0d0
         dyr12r22 = (dyr(1,2)*r(2,2) + r(1,2)*dyr(2,2))*2.0d0
         dzr12r22 = (dzr(1,2)*r(2,2) + r(1,2)*dzr(2,2))*2.0d0
        r12r23 = r(1,2)*r(2,3)
         dxr12r23 = dxr(1,2)*r(2,3) + r(1,2)*dxr(2,3)
         dyr12r23 = dyr(1,2)*r(2,3) + r(1,2)*dyr(2,3)
         dzr12r23 = dzr(1,2)*r(2,3) + r(1,2)*dzr(2,3)
        r12r31 = r(1,2)*r(3,1)
         dxr12r31 = dxr(1,2)*r(3,1) + r(1,2)*dxr(3,1)
         dyr12r31 = dyr(1,2)*r(3,1) + r(1,2)*dyr(3,1)
         dzr12r31 = dzr(1,2)*r(3,1) + r(1,2)*dzr(3,1)
        r12r32 = r(1,2)*r(3,2)*2.0d0
         dxr12r32 = (dxr(1,2)*r(3,2) + r(1,2)*dxr(3,2))*2.0d0
         dyr12r32 = (dyr(1,2)*r(3,2) + r(1,2)*dyr(3,2))*2.0d0
         dzr12r32 = (dzr(1,2)*r(3,2) + r(1,2)*dzr(3,2))*2.0d0
        r12r33 = r(1,2)*r(3,3)
         dxr12r33 = dxr(1,2)*r(3,3) + r(1,2)*dxr(3,3)
         dyr12r33 = dyr(1,2)*r(3,3) + r(1,2)*dyr(3,3)
         dzr12r33 = dzr(1,2)*r(3,3) + r(1,2)*dzr(3,3)
        r13r23 = r(1,3)*r(2,3)*2.0d0
         dxr13r23 = (dxr(1,3)*r(2,3) + r(1,3)*dxr(2,3))*2.0d0
         dyr13r23 = (dyr(1,3)*r(2,3) + r(1,3)*dyr(2,3))*2.0d0
         dzr13r23 = (dzr(1,3)*r(2,3) + r(1,3)*dzr(2,3))*2.0d0
        r13r33 = r(1,3)*r(3,3)*2.0d0
         dxr13r33 = (dxr(1,3)*r(3,3) + r(1,3)*dxr(3,3))*2.0d0
         dyr13r33 = (dyr(1,3)*r(3,3) + r(1,3)*dyr(3,3))*2.0d0
         dzr13r33 = (dzr(1,3)*r(3,3) + r(1,3)*dzr(3,3))*2.0d0
        r21r11 = r11r21
         dxr21r11 = dxr11r21 
         dyr21r11 = dyr11r21 
         dzr21r11 = dzr11r21 
        r21r12 = r(2,1)*r(1,2)
         dxr21r12 = dxr(2,1)*r(1,2) + r(2,1)*dxr(1,2)
         dyr21r12 = dyr(2,1)*r(1,2) + r(2,1)*dyr(1,2)
         dzr21r12 = dzr(2,1)*r(1,2) + r(2,1)*dzr(1,2)
        r21r13 = r(2,1)*r(1,3)
         dxr21r13 = dxr(2,1)*r(1,3) + r(2,1)*dxr(1,3)
         dyr21r13 = dyr(2,1)*r(1,3) + r(2,1)*dyr(1,3)
         dzr21r13 = dzr(2,1)*r(1,3) + r(2,1)*dzr(1,3)
        r21r21 = r(2,1)*r(2,1)
         dxr21r21 = dxr(2,1)*r(2,1) + r(2,1)*dxr(2,1)
         dyr21r21 = dyr(2,1)*r(2,1) + r(2,1)*dyr(2,1)
         dzr21r21 = dzr(2,1)*r(2,1) + r(2,1)*dzr(2,1)
        r21r31 = r(2,1)*r(3,1)*2.0d0
         dxr21r31 = (dxr(2,1)*r(3,1) + r(2,1)*dxr(3,1))*2.0d0
         dyr21r31 = (dyr(2,1)*r(3,1) + r(2,1)*dyr(3,1))*2.0d0
         dzr21r31 = (dzr(2,1)*r(3,1) + r(2,1)*dzr(3,1))*2.0d0
        r21r32 = r(2,1)*r(3,2)
         dxr21r32 = dxr(2,1)*r(3,2) + r(2,1)*dxr(3,2)
         dyr21r32 = dyr(2,1)*r(3,2) + r(2,1)*dyr(3,2)
         dzr21r32 = dzr(2,1)*r(3,2) + r(2,1)*dzr(3,2)
        r21r33 = r(2,1)*r(3,3)
         dxr21r33 = dxr(2,1)*r(3,3) + r(2,1)*dxr(3,3)
         dyr21r33 = dyr(2,1)*r(3,3) + r(2,1)*dyr(3,3)
         dzr21r33 = dzr(2,1)*r(3,3) + r(2,1)*dzr(3,3)
        r22r12 = r12r22
         dxr22r12 = dxr12r22 
         dyr22r12 = dyr12r22 
         dzr22r12 = dzr12r22 
        r22r13 = r(2,2)*r(1,3)
         dxr22r13 = dxr(2,2)*r(1,3) + r(2,2)*dxr(1,3)
         dyr22r13 = dyr(2,2)*r(1,3) + r(2,2)*dyr(1,3)
         dzr22r13 = dzr(2,2)*r(1,3) + r(2,2)*dzr(1,3)
        r22r31 = r(2,2)*r(3,1)
         dxr22r31 = dxr(2,2)*r(3,1) + r(2,2)*dxr(3,1)
         dyr22r31 = dyr(2,2)*r(3,1) + r(2,2)*dyr(3,1)
         dzr22r31 = dzr(2,2)*r(3,1) + r(2,2)*dzr(3,1)
        r22r32 = r(2,2)*r(3,2)*2.0d0
         dxr22r32 = (dxr(2,2)*r(3,2) + r(2,2)*dxr(3,2))*2.0d0
         dyr22r32 = (dyr(2,2)*r(3,2) + r(2,2)*dyr(3,2))*2.0d0
         dzr22r32 = (dzr(2,2)*r(3,2) + r(2,2)*dzr(3,2))*2.0d0
        r22r33 = r(2,2)*r(3,3)
         dxr22r33 = dxr(2,2)*r(3,3) + r(2,2)*dxr(3,3)
         dyr22r33 = dyr(2,2)*r(3,3) + r(2,2)*dyr(3,3)
         dzr22r33 = dzr(2,2)*r(3,3) + r(2,2)*dzr(3,3)
        r23r13 = r13r23
         dxr23r13 = dxr13r23 
         dyr23r13 = dyr13r23 
         dzr23r13 = dzr13r23 
        r23r33 = r(2,3)*r(3,3)*2.0d0
         dxr23r33 = (dxr(2,3)*r(3,3) + r(2,3)*dxr(3,3))*2.0d0
         dyr23r33 = (dyr(2,3)*r(3,3) + r(2,3)*dyr(3,3))*2.0d0
         dzr23r33 = (dzr(2,3)*r(3,3) + r(2,3)*dzr(3,3))*2.0d0
        r31r11 = r11r31
         dxr31r11 = dxr11r31
         dyr31r11 = dyr11r31
         dzr31r11 = dzr11r31
        r31r12 = r(3,1)*r(1,2)
         dxr31r12 = dxr(3,1)*r(1,2) + r(3,1)*dxr(1,2)
         dyr31r12 = dyr(3,1)*r(1,2) + r(3,1)*dyr(1,2)
         dzr31r12 = dzr(3,1)*r(1,2) + r(3,1)*dzr(1,2)
        r31r13 = r(3,1)*r(1,3)
         dxr31r13 = dxr(3,1)*r(1,3) + r(3,1)*dxr(1,3)
         dyr31r13 = dyr(3,1)*r(1,3) + r(3,1)*dyr(1,3)
         dzr31r13 = dzr(3,1)*r(1,3) + r(3,1)*dzr(1,3)
        r31r21 = r21r31
         dxr31r21 = dxr21r31 
         dyr31r21 = dyr21r31 
         dzr31r21 = dzr21r31 
        r31r22 = r22r31
         dxr31r22 = dxr22r31 
         dyr31r22 = dyr22r31 
         dzr31r22 = dzr22r31 
        r31r23 = r(3,1)*r(2,3)
         dxr31r23 = dxr(3,1)*r(2,3) + r(3,1)*dxr(2,3)
         dyr31r23 = dyr(3,1)*r(2,3) + r(3,1)*dyr(2,3)
         dzr31r23 = dzr(3,1)*r(2,3) + r(3,1)*dzr(2,3)
        r32r12 = r12r32
         dxr32r12 = dxr12r32 
         dyr32r12 = dyr12r32 
         dzr32r12 = dzr12r32 
        r32r13 = r(3,2)*r(1,3)
         dxr32r13 = dxr(3,2)*r(1,3) + r(3,2)*dxr(1,3)
         dyr32r13 = dyr(3,2)*r(1,3) + r(3,2)*dyr(1,3)
         dzr32r13 = dzr(3,2)*r(1,3) + r(3,2)*dzr(1,3)
        r32r22 = r22r32
         dxr32r22 = dxr22r32 
         dyr32r22 = dyr22r32 
         dzr32r22 = dzr22r32 
        r32r23 = r(3,2)*r(2,3)
         dxr32r23 = dxr(3,2)*r(2,3) + r(3,2)*dxr(2,3)
         dyr32r23 = dyr(3,2)*r(2,3) + r(3,2)*dyr(2,3)
         dzr32r23 = dzr(3,2)*r(2,3) + r(3,2)*dzr(2,3)
        r33r13 = r13r33
         dxr33r13 = dxr13r33 
         dyr33r13 = dyr13r33 
         dzr33r13 = dzr13r33 
        r33r23 = r23r33
         dxr33r23 = dxr23r33 
         dyr33r23 = dyr23r33 
         dzr33r23 = dzr23r33 
        radd07 = r11r32+r12r31
         dxradd07 = dxr11r32+dxr12r31
         dyradd07 = dyr11r32+dyr12r31
         dzradd07 = dzr11r32+dzr12r31
        radd08 = r21r32+r22r31
         dxradd08 = dxr21r32+dxr22r31
         dyradd08 = dyr21r32+dyr22r31
         dzradd08 = dzr21r32+dzr22r31
        radd09 = r11r33+r31r13
         dxradd09 = dxr11r33+dxr31r13
         dyradd09 = dyr11r33+dyr31r13
         dzradd09 = dzr11r33+dzr31r13
        radd10 = r21r33+r31r23
         dxradd10 = dxr21r33+dxr31r23
         dyradd10 = dyr21r33+dyr31r23
         dzradd10 = dzr21r33+dzr31r23
        radd11 = r11r32+r31r12
         dxradd11 = dxr11r32+dxr31r12
         dyradd11 = dyr11r32+dyr31r12 
         dzradd11 = dzr11r32+dzr31r12 
        radd12 = r21r32+r31r22
         dxradd12 = dxr21r32+dxr31r22 
         dyradd12 = dyr21r32+dyr31r22 
         dzradd12 = dzr21r32+dzr31r22 
        radd13 = r12r33+r32r13
         dxradd13 = dxr12r33+dxr32r13 
         dyradd13 = dyr12r33+dyr32r13 
         dzradd13 = dzr12r33+dzr32r13 
        radd14 = r22r33+r32r23
         dxradd14 = dxr22r33+dxr32r23 
         dyradd14 = dyr22r33+dyr32r23 
         dzradd14 = dzr22r33+dzr32r23 
        radd15 = r11r23+r21r13
         dxradd15 = dxr11r23+dxr21r13 
         dyradd15 = dyr11r23+dyr21r13 
         dzradd15 = dzr11r23+dzr21r13 
        radd16 = r12r23+r22r13
         dxradd16 = dxr12r23+dxr22r13 
         dyradd16 = dyr12r23+dyr22r13 
         dzradd16 = dzr12r23+dzr22r13 
        radd17 = r11r22+r21r12
         dxradd17 = dxr11r22+dxr21r12 
         dyradd17 = dyr11r22+dyr21r12 
         dzradd17 = dzr11r22+dzr21r12 
      endif
!
!     (s s | s s):
!
      repij(1,1) = reploc(1)
       dxrepij(1,1) = dxreploc(1)
       dyrepij(1,1) = dyreploc(1)
       dzrepij(1,1) = dzreploc(1)
      if(norbsj.gt.1)then
!
!       (s s | s px):
!
        repij(1,2) = reploc(5)*r(3,1)
         dxrepij(1,2) = dxreploc(5)*r(3,1) + reploc(5)*dxr(3,1)
         dyrepij(1,2) = dyreploc(5)*r(3,1) + reploc(5)*dyr(3,1)
         dzrepij(1,2) = dzreploc(5)*r(3,1) + reploc(5)*dzr(3,1)
!
!       (s s | s py):
!
        repij(1,3) = reploc(5)*r(3,2)
         dxrepij(1,3) = dxreploc(5)*r(3,2) + reploc(5)*dxr(3,2)
         dyrepij(1,3) = dyreploc(5)*r(3,2) + reploc(5)*dyr(3,2)
         dzrepij(1,3) = dzreploc(5)*r(3,2) + reploc(5)*dzr(3,2)
!
!       (s s | s pz):
!
        repij(1,4) = reploc(5)*r(3,3)
         dxrepij(1,4) = dxreploc(5)*r(3,3) + reploc(5)*dxr(3,3)
         dyrepij(1,4) = dyreploc(5)*r(3,3) + reploc(5)*dyr(3,3)
         dzrepij(1,4) = dzreploc(5)*r(3,3) + reploc(5)*dzr(3,3)
!
!       (s s | px px):
!
        repij(1,5) = reploc(11)*r31r31
     .             + reploc(12)*radd01
         dxrepij(1,5) = dxreploc(11)*r31r31 + reploc(11)*dxr31r31
     .                + dxreploc(12)*radd01 + reploc(12)*dxradd01
         dyrepij(1,5) = dyreploc(11)*r31r31 + reploc(11)*dyr31r31
     .                + dyreploc(12)*radd01 + reploc(12)*dyradd01
         dzrepij(1,5) = dzreploc(11)*r31r31 + reploc(11)*dzr31r31
     .                + dzreploc(12)*radd01 + reploc(12)*dzradd01
!
!       (s s | px py):
!
        repij(1,6) = reploc(11)*r31r32
     .             + reploc(12)*radd02
         dxrepij(1,6) = dxreploc(11)*r31r32 + reploc(11)*dxr31r32
     .                + dxreploc(12)*radd02 + reploc(12)*dxradd02
         dyrepij(1,6) = dyreploc(11)*r31r32 + reploc(11)*dyr31r32
     .                + dyreploc(12)*radd02 + reploc(12)*dyradd02
         dzrepij(1,6) = dzreploc(11)*r31r32 + reploc(11)*dzr31r32
     .                + dzreploc(12)*radd02 + reploc(12)*dzradd02
!
!       (s s | px pz):
!
        repij(1,7) = reploc(11)*r31r33
     .             + reploc(12)*radd03
         dxrepij(1,7) = dxreploc(11)*r31r33 + reploc(11)*dxr31r33
     .                + dxreploc(12)*radd03 + reploc(12)*dxradd03
         dyrepij(1,7) = dyreploc(11)*r31r33 + reploc(11)*dyr31r33
     .                + dyreploc(12)*radd03 + reploc(12)*dyradd03
         dzrepij(1,7) = dzreploc(11)*r31r33 + reploc(11)*dzr31r33
     .                + dzreploc(12)*radd03 + reploc(12)*dzradd03
!
!       (s s | py py):
!
        repij(1,8) = reploc(11)*r32r32
     .             + reploc(12)*radd04
         dxrepij(1,8) = dxreploc(11)*r32r32 + reploc(11)*dxr32r32
     .                + dxreploc(12)*radd04 + reploc(12)*dxradd04
         dyrepij(1,8) = dyreploc(11)*r32r32 + reploc(11)*dyr32r32
     .                + dyreploc(12)*radd04 + reploc(12)*dyradd04
         dzrepij(1,8) = dzreploc(11)*r32r32 + reploc(11)*dzr32r32
     .                + dzreploc(12)*radd04 + reploc(12)*dzradd04
!
!       (s s | py pz):
!
        repij(1,9) = reploc(11)*r32r33
     .             + reploc(12)*radd05
         dxrepij(1,9) = dxreploc(11)*r32r33 + reploc(11)*dxr32r33
     .                + dxreploc(12)*radd05 + reploc(12)*dxradd05
         dyrepij(1,9) = dyreploc(11)*r32r33 + reploc(11)*dyr32r33
     .                + dyreploc(12)*radd05 + reploc(12)*dyradd05
         dzrepij(1,9) = dzreploc(11)*r32r33 + reploc(11)*dzr32r33
     .                + dzreploc(12)*radd05 + reploc(12)*dzradd05
!
!       (s s | pz pz):
!
        repij(1,10) = reploc(11)*r33r33
     .              + reploc(12)*radd06
         dxrepij(1,10) = dxreploc(11)*r33r33 + reploc(11)*dxr33r33
     .                 + dxreploc(12)*radd06 + reploc(12)*dxradd06
         dyrepij(1,10) = dyreploc(11)*r33r33 + reploc(11)*dyr33r33
     .                 + dyreploc(12)*radd06 + reploc(12)*dyradd06
         dzrepij(1,10) = dzreploc(11)*r33r33 + reploc(11)*dzr33r33
     .                 + dzreploc(12)*radd06 + reploc(12)*dzradd06
         
      endif
      if(norbsi.gt.1)then
!
!       (s px | s s):
!
        repij(2,1) = reploc(2)*r(3,1)
         dxrepij(2,1) = dxreploc(2)*r(3,1) + reploc(2)*dxr(3,1)
         dyrepij(2,1) = dyreploc(2)*r(3,1) + reploc(2)*dyr(3,1)
         dzrepij(2,1) = dzreploc(2)*r(3,1) + reploc(2)*dzr(3,1)
!
!       (s py | s s):
!
        repij(3,1) = reploc(2)*r(3,2)
         dxrepij(3,1) = dxreploc(2)*r(3,2) + reploc(2)*dxr(3,2)
         dyrepij(3,1) = dyreploc(2)*r(3,2) + reploc(2)*dyr(3,2)
         dzrepij(3,1) = dzreploc(2)*r(3,2) + reploc(2)*dzr(3,2)
!
!       (s pz | s s):
!
        repij(4,1) = reploc(2)*r(3,3)
         dxrepij(4,1) = dxreploc(2)*r(3,3) + reploc(2)*dxr(3,3)
         dyrepij(4,1) = dyreploc(2)*r(3,3) + reploc(2)*dyr(3,3)
         dzrepij(4,1) = dzreploc(2)*r(3,3) + reploc(2)*dzr(3,3)
!
!       (px px | s s):
!
        repij(5,1) = reploc(3)*r31r31 
     .             + reploc(4)*radd01 
         dxrepij(5,1) = dxreploc(3)*r31r31 + reploc(3)*dxr31r31
     .                + dxreploc(4)*radd01 + reploc(4)*dxradd01
         dyrepij(5,1) = dyreploc(3)*r31r31 + reploc(3)*dyr31r31
     .                + dyreploc(4)*radd01 + reploc(4)*dyradd01
         dzrepij(5,1) = dzreploc(3)*r31r31 + reploc(3)*dzr31r31
     .                + dzreploc(4)*radd01 + reploc(4)*dzradd01
!
!       (px py | s s):
!
        repij(6,1) = reploc(3)*r31r32
     .             + reploc(4)*radd02
         dxrepij(6,1) = dxreploc(3)*r31r32 + reploc(3)*dxr31r32
     .                + dxreploc(4)*radd02 + reploc(4)*dxradd02
         dyrepij(6,1) = dyreploc(3)*r31r32 + reploc(3)*dyr31r32
     .                + dyreploc(4)*radd02 + reploc(4)*dyradd02
         dzrepij(6,1) = dzreploc(3)*r31r32 + reploc(3)*dzr31r32
     .                + dzreploc(4)*radd02 + reploc(4)*dzradd02
!
!       (px pz | s s):
!
        repij(7,1) = reploc(3)*r31r33
     .             + reploc(4)*radd03
         dxrepij(7,1) = dxreploc(3)*r31r33 + reploc(3)*dxr31r33
     .                + dxreploc(4)*radd03 + reploc(4)*dxradd03
         dyrepij(7,1) = dyreploc(3)*r31r33 + reploc(3)*dyr31r33
     .                + dyreploc(4)*radd03 + reploc(4)*dyradd03
         dzrepij(7,1) = dzreploc(3)*r31r33 + reploc(3)*dzr31r33
     .                + dzreploc(4)*radd03 + reploc(4)*dzradd03
!
!       (py py | s s):
!
        repij(8,1) = reploc(3)*r32r32
     .             + reploc(4)*radd04
         dxrepij(8,1) = dxreploc(3)*r32r32 + reploc(3)*dxr32r32
     .                + dxreploc(4)*radd04 + reploc(4)*dxradd04
         dyrepij(8,1) = dyreploc(3)*r32r32 + reploc(3)*dyr32r32
     .                + dyreploc(4)*radd04 + reploc(4)*dyradd04
         dzrepij(8,1) = dzreploc(3)*r32r32 + reploc(3)*dzr32r32
     .                + dzreploc(4)*radd04 + reploc(4)*dzradd04
!
!       (py pz | s s);
!
        repij(9,1) = reploc(3)*r32r33
     .             + reploc(4)*radd05
         dxrepij(9,1) = dxreploc(3)*r32r33 + reploc(3)*dxr32r33
     .                + dxreploc(4)*radd05 + reploc(4)*dxradd05
         dyrepij(9,1) = dyreploc(3)*r32r33 + reploc(3)*dyr32r33
     .                + dyreploc(4)*radd05 + reploc(4)*dyradd05
         dzrepij(9,1) = dzreploc(3)*r32r33 + reploc(3)*dzr32r33
     .                + dzreploc(4)*radd05 + reploc(4)*dzradd05
!
!       (pz pz | s s):
!
        repij(10,1) = reploc(3)*r33r33
     .              + reploc(4)*radd06
         dxrepij(10,1) = dxreploc(3)*r33r33 + reploc(3)*dxr33r33
     .                 + dxreploc(4)*radd06 + reploc(4)*dxradd06
         dyrepij(10,1) = dyreploc(3)*r33r33 + reploc(3)*dyr33r33
     .                 + dyreploc(4)*radd06 + reploc(4)*dyradd06
         dzrepij(10,1) = dzreploc(3)*r33r33 + reploc(3)*dzr33r33
     .                 + dzreploc(4)*radd06 + reploc(4)*dzradd06
      endif
      if(norbsi.gt.1.and.norbsj.gt.1)then
!
!       (s px | s px):
!
        repij(2,2) = reploc(6)*r31r31 
     .             + reploc(7)*radd01 
         dxrepij(2,2) = dxreploc(6)*r31r31 + reploc(6)*dxr31r31
     .                + dxreploc(7)*radd01 + reploc(7)*dxradd01
         dyrepij(2,2) = dyreploc(6)*r31r31 + reploc(6)*dyr31r31
     .                + dyreploc(7)*radd01 + reploc(7)*dyradd01
         dzrepij(2,2) = dzreploc(6)*r31r31 + reploc(6)*dzr31r31
     .                + dzreploc(7)*radd01 + reploc(7)*dzradd01
!
!       (s px | s py):
!
        repij(2,3) = reploc(6)*r31r32
     .             + reploc(7)*radd02
         dxrepij(2,3) = dxreploc(6)*r31r32 + reploc(6)*dxr31r32
     .                + dxreploc(7)*radd02 + reploc(7)*dxradd02
         dyrepij(2,3) = dyreploc(6)*r31r32 + reploc(6)*dyr31r32
     .                + dyreploc(7)*radd02 + reploc(7)*dyradd02
         dzrepij(2,3) = dzreploc(6)*r31r32 + reploc(6)*dzr31r32
     .                + dzreploc(7)*radd02 + reploc(7)*dzradd02
! 
!       (s py | s px):
!
        repij(3,2) = repij(2,3)
         dxrepij(3,2) = dxrepij(2,3) 
         dyrepij(3,2) = dyrepij(2,3)
         dzrepij(3,2) = dzrepij(2,3)
!
!       (s px | s pz):
!
        repij(2,4) = reploc(6)*r31r33
     .             + reploc(7)*radd03
         dxrepij(2,4) = dxreploc(6)*r31r33 + reploc(6)*dxr31r33
     .                + dxreploc(7)*radd03 + reploc(7)*dxradd03
         dyrepij(2,4) = dyreploc(6)*r31r33 + reploc(6)*dyr31r33
     .                + dyreploc(7)*radd03 + reploc(7)*dyradd03
         dzrepij(2,4) = dzreploc(6)*r31r33 + reploc(6)*dzr31r33
     .                + dzreploc(7)*radd03 + reploc(7)*dzradd03
!
!       (s pz | s px):
!
        repij(4,2) = repij(2,4)
         dxrepij(4,2) = dxrepij(2,4)
         dyrepij(4,2) = dyrepij(2,4)
         dzrepij(4,2) = dzrepij(2,4)
!
!       (s px | px px):
!
        repij(2,5) = reploc(13)*r(3,1)*r31r31 
     .             + reploc(14)*r(3,1)*radd01
     .             + reploc(15)*(r(1,1)*r11r31
     .                         + r(2,1)*r21r31)
         dxrepij(2,5) = dxreploc(13)*r(3,1)*r31r31
     .               + reploc(13)*dxr(3,1)*r31r31
     .               + reploc(13)*r(3,1)*dxr31r31
     .              + dxreploc(14)*r(3,1)*radd01
     .               + reploc(14)*dxr(3,1)*radd01
     .               + reploc(14)*r(3,1)*dxradd01
     .              + dxreploc(15)*(r(1,1)*r11r31
     .                          + r(2,1)*r21r31)
     .               + reploc(15)*(dxr(1,1)*r11r31
     .                           + dxr(2,1)*r21r31)
     .               + reploc(15)*(r(1,1)*dxr11r31
     .                           + r(2,1)*dxr21r31)

         dyrepij(2,5) = dyreploc(13)*r(3,1)*r31r31
     .               + reploc(13)*dyr(3,1)*r31r31  
     .               + reploc(13)*r(3,1)*dyr31r31  
     .              + dyreploc(14)*r(3,1)*radd01
     .               + reploc(14)*dyr(3,1)*radd01
     .               + reploc(14)*r(3,1)*dyradd01
     .              + dyreploc(15)*(r(1,1)*r11r31
     .                          + r(2,1)*r21r31)
     .               + reploc(15)*(dyr(1,1)*r11r31
     .                           + dyr(2,1)*r21r31)
     .               + reploc(15)*(r(1,1)*dyr11r31
     .                           + r(2,1)*dyr21r31)

         dzrepij(2,5) = dzreploc(13)*r(3,1)*r31r31
     .               + reploc(13)*dzr(3,1)*r31r31
     .               + reploc(13)*r(3,1)*dzr31r31
     .              + dzreploc(14)*r(3,1)*radd01
     .               + reploc(14)*dzr(3,1)*radd01
     .               + reploc(14)*r(3,1)*dzradd01
     .              + dzreploc(15)*(r(1,1)*r11r31
     .                          + r(2,1)*r21r31)
     .               + reploc(15)*(dzr(1,1)*r11r31
     .                           + dzr(2,1)*r21r31)
     .               + reploc(15)*(r(1,1)*dzr11r31
     .                           + r(2,1)*dzr21r31)
!
!       (px px | s px):
!
        repij(5,2) = reploc(8)*r(3,1)*r31r31
     .             + reploc(9)*r(3,1)*radd01
     .             + reploc(10)*(r(1,1)*r11r31
     .                         + r(2,1)*r21r31)

         dxrepij(5,2) = dxreploc(8)*r(3,1)*r31r31
     .               + reploc(8)*dxr(3,1)*r31r31
     .               + reploc(8)*r(3,1)*dxr31r31
     .              + dxreploc(9)*r(3,1)*radd01
     .               + reploc(9)*dxr(3,1)*radd01
     .               + reploc(9)*r(3,1)*dxradd01
     .              + dxreploc(10)*(r(1,1)*r11r31
     .                          + r(2,1)*r21r31)
     .               + reploc(10)*(dxr(1,1)*r11r31
     .                           + dxr(2,1)*r21r31)
     .               + reploc(10)*(r(1,1)*dxr11r31
     .                           + r(2,1)*dxr21r31)

         dyrepij(5,2) = dyreploc(8)*r(3,1)*r31r31
     .               + reploc(8)*dyr(3,1)*r31r31
     .               + reploc(8)*r(3,1)*dyr31r31
     .              + dyreploc(9)*r(3,1)*radd01
     .               + reploc(9)*dyr(3,1)*radd01
     .               + reploc(9)*r(3,1)*dyradd01
     .              + dyreploc(10)*(r(1,1)*r11r31
     .                          + r(2,1)*r21r31)
     .               + reploc(10)*(dyr(1,1)*r11r31
     .                           + dyr(2,1)*r21r31)
     .               + reploc(10)*(r(1,1)*dyr11r31
     .                           + r(2,1)*dyr21r31)

         dzrepij(5,2) = dzreploc(8)*r(3,1)*r31r31
     .               + reploc(8)*dzr(3,1)*r31r31
     .               + reploc(8)*r(3,1)*dzr31r31
     .              + dzreploc(9)*r(3,1)*radd01
     .               + reploc(9)*dzr(3,1)*radd01
     .               + reploc(9)*r(3,1)*dzradd01
     .              + dzreploc(10)*(r(1,1)*r11r31
     .                          + r(2,1)*r21r31)
     .               + reploc(10)*(dzr(1,1)*r11r31
     .                           + dzr(2,1)*r21r31)
     .               + reploc(10)*(r(1,1)*dzr11r31
     .                           + r(2,1)*dzr21r31)
!
!       (s px | px py):
!
        repij(2,6) = reploc(13)*r(3,1)*r31r32
     .              + reploc(14)*r(3,1)*radd02
     .              + reploc(15)*(r(1,1)*radd07
     .                          + r(2,1)*radd08)

         dxrepij(2,6) = dxreploc(13)*r(3,1)*r31r32
     .               + reploc(13)*dxr(3,1)*r31r32
     .               + reploc(13)*r(3,1)*dxr31r32
     .              + dxreploc(14)*r(3,1)*radd02
     .               + reploc(14)*dxr(3,1)*radd02
     .               + reploc(14)*r(3,1)*dxradd02
     .              + dxreploc(15)*(r(1,1)*radd07
     .                          + r(2,1)*radd08)
     .               + reploc(15)*(dxr(1,1)*radd07
     .                           + dxr(2,1)*radd08)
     .               + reploc(15)*(r(1,1)*dxradd07
     .                           + r(2,1)*dxradd08)

         dyrepij(2,6) = dyreploc(13)*r(3,1)*r31r32
     .               + reploc(13)*dyr(3,1)*r31r32
     .               + reploc(13)*r(3,1)*dyr31r32
     .              + dyreploc(14)*r(3,1)*radd02
     .               + reploc(14)*dyr(3,1)*radd02
     .               + reploc(14)*r(3,1)*dyradd02
     .              + dyreploc(15)*(r(1,1)*radd07
     .                          + r(2,1)*radd08)
     .               + reploc(15)*(dyr(1,1)*radd07
     .                           + dyr(2,1)*radd08)
     .               + reploc(15)*(r(1,1)*dyradd07
     .                           + r(2,1)*dyradd08)

         dzrepij(2,6) = dzreploc(13)*r(3,1)*r31r32
     .               + reploc(13)*dzr(3,1)*r31r32
     .               + reploc(13)*r(3,1)*dzr31r32
     .              + dzreploc(14)*r(3,1)*radd02
     .               + reploc(14)*dzr(3,1)*radd02
     .               + reploc(14)*r(3,1)*dzradd02
     .              + dzreploc(15)*(r(1,1)*radd07
     .                          + r(2,1)*radd08)
     .               + reploc(15)*(dzr(1,1)*radd07
     .                           + dzr(2,1)*radd08)
     .               + reploc(15)*(r(1,1)*dzradd07
     .                           + r(2,1)*dzradd08)
!
!       (px py | s px):
!
        repij(6,2) = reploc(8)*r(3,1)*r31r32
     .             + reploc(9)*r(3,1)*radd02
     .             + reploc(10)*(r(1,1)*radd07
     .                         + r(2,1)*radd08)

         dxrepij(6,2) = dxreploc(8)*r(3,1)*r31r32
     .               + reploc(8)*dxr(3,1)*r31r32
     .               + reploc(8)*r(3,1)*dxr31r32
     .              + dxreploc(9)*r(3,1)*radd02
     .               + reploc(9)*dxr(3,1)*radd02
     .               + reploc(9)*r(3,1)*dxradd02
     .              + dxreploc(10)*(r(1,1)*radd07
     .                          + r(2,1)*radd08)
     .               + reploc(10)*(dxr(1,1)*radd07
     .                           + dxr(2,1)*radd08)
     .               + reploc(10)*(r(1,1)*dxradd07
     .                           + r(2,1)*dxradd08)

         dyrepij(6,2) = dyreploc(8)*r(3,1)*r31r32
     .               + reploc(8)*dyr(3,1)*r31r32
     .               + reploc(8)*r(3,1)*dyr31r32
     .              + dyreploc(9)*r(3,1)*radd02
     .               + reploc(9)*dyr(3,1)*radd02
     .               + reploc(9)*r(3,1)*dyradd02
     .              + dyreploc(10)*(r(1,1)*radd07
     .                          + r(2,1)*radd08)
     .               + reploc(10)*(dyr(1,1)*radd07
     .                           + dyr(2,1)*radd08)
     .               + reploc(10)*(r(1,1)*dyradd07
     .                           + r(2,1)*dyradd08)

         dzrepij(6,2) = dzreploc(8)*r(3,1)*r31r32
     .               + reploc(8)*dzr(3,1)*r31r32
     .               + reploc(8)*r(3,1)*dzr31r32
     .              + dzreploc(9)*r(3,1)*radd02
     .               + reploc(9)*dzr(3,1)*radd02
     .               + reploc(9)*r(3,1)*dzradd02
     .              + dzreploc(10)*(r(1,1)*radd07
     .                          + r(2,1)*radd08)
     .               + reploc(10)*(dzr(1,1)*radd07
     .                           + dzr(2,1)*radd08)
     .               + reploc(10)*(r(1,1)*dzradd07
     .                           + r(2,1)*dzradd08)
!
!       (s px | px pz):
!
        repij(2,7) = reploc(13)*r(3,1)*r31r33
     .             + reploc(14)*r(3,1)*radd03
     .             + reploc(15)*(r(1,1)*radd09
     .                         + r(2,1)*radd10)

        dxrepij(2,7) = dxreploc(13)*r(3,1)*r31r33
     .              + reploc(13)*dxr(3,1)*r31r33
     .              + reploc(13)*r(3,1)*dxr31r33
     .             + dxreploc(14)*r(3,1)*radd03
     .              + reploc(14)*dxr(3,1)*radd03
     .              + reploc(14)*r(3,1)*dxradd03
     .             + dxreploc(15)*(r(1,1)*radd09
     .                         + r(2,1)*radd10)
     .              + reploc(15)*(dxr(1,1)*radd09
     .                          + dxr(2,1)*radd10)
     .              + reploc(15)*(r(1,1)*dxradd09
     .                          + r(2,1)*dxradd10)

        dyrepij(2,7) = dyreploc(13)*r(3,1)*r31r33
     .              + reploc(13)*dyr(3,1)*r31r33
     .              + reploc(13)*r(3,1)*dyr31r33
     .             + dyreploc(14)*r(3,1)*radd03
     .              + reploc(14)*dyr(3,1)*radd03
     .              + reploc(14)*r(3,1)*dyradd03
     .             + dyreploc(15)*(r(1,1)*radd09
     .                         + r(2,1)*radd10)
     .              + reploc(15)*(dyr(1,1)*radd09
     .                          + dyr(2,1)*radd10)
     .              + reploc(15)*(r(1,1)*dyradd09
     .                          + r(2,1)*dyradd10)

        dzrepij(2,7) = dzreploc(13)*r(3,1)*r31r33
     .              + reploc(13)*dzr(3,1)*r31r33
     .              + reploc(13)*r(3,1)*dzr31r33
     .             + dzreploc(14)*r(3,1)*radd03
     .              + reploc(14)*dzr(3,1)*radd03
     .              + reploc(14)*r(3,1)*dzradd03
     .             + dzreploc(15)*(r(1,1)*radd09
     .                         + r(2,1)*radd10)
     .              + reploc(15)*(dzr(1,1)*radd09
     .                          + dzr(2,1)*radd10)
     .              + reploc(15)*(r(1,1)*dzradd09
     .                          + r(2,1)*dzradd10)

!
!       (px pz | s px):
!
        repij(7,2) = reploc(8)*r(3,1)*r31r33
     .             + reploc(9)*(r(3,1)*r11r13
     .                        + r(3,1)*r21r23)
     .             + reploc(10)*(r(1,1)*radd09
     .                         + r(2,1)*radd10)

        dxrepij(7,2) = dxreploc(8)*r(3,1)*r31r33
     .              + reploc(8)*dxr(3,1)*r31r33 
     .              + reploc(8)*r(3,1)*dxr31r33 
     .             + dxreploc(9)*(r(3,1)*r11r13
     .                        + r(3,1)*r21r23)
     .              + reploc(9)*(dxr(3,1)*r11r13
     .                         + dxr(3,1)*r21r23)
     .              + reploc(9)*(r(3,1)*dxr11r13
     .                         + r(3,1)*dxr21r23)
     .             + dxreploc(10)*(r(1,1)*radd09
     .                        + r(2,1)*radd10)
     .              + reploc(10)*(dxr(1,1)*radd09
     .                         + dxr(2,1)*radd10)
     .              + reploc(10)*(r(1,1)*dxradd09
     .                         + r(2,1)*dxradd10)

        dyrepij(7,2) = dyreploc(8)*r(3,1)*r31r33
     .              + reploc(8)*dyr(3,1)*r31r33 
     .              + reploc(8)*r(3,1)*dyr31r33 
     .             + dyreploc(9)*(r(3,1)*r11r13
     .                        + r(3,1)*r21r23)
     .              + reploc(9)*(dyr(3,1)*r11r13
     .                         + dyr(3,1)*r21r23)
     .              + reploc(9)*(r(3,1)*dyr11r13
     .                         + r(3,1)*dyr21r23)
     .             + dyreploc(10)*(r(1,1)*radd09
     .                         + r(2,1)*radd10)
     .             + reploc(10)*(dyr(1,1)*radd09
     .                         + dyr(2,1)*radd10)
     .             + reploc(10)*(r(1,1)*dyradd09
     .                         + r(2,1)*dyradd10)

        dzrepij(7,2) = dzreploc(8)*r(3,1)*r31r33
     .              + reploc(8)*dzr(3,1)*r31r33 
     .              + reploc(8)*r(3,1)*dzr31r33 
     .             + dzreploc(9)*(r(3,1)*r11r13
     .                        + r(3,1)*r21r23)
     .              + reploc(9)*(dzr(3,1)*r11r13
     .                         + dzr(3,1)*r21r23)
     .              + reploc(9)*(r(3,1)*dzr11r13
     .                         + r(3,1)*dzr21r23)
     .             + dzreploc(10)*(r(1,1)*radd09
     .                         + r(2,1)*radd10)
     .             + reploc(10)*(dzr(1,1)*radd09
     .                         + dzr(2,1)*radd10)
     .             + reploc(10)*(r(1,1)*dzradd09
     .                         + r(2,1)*dzradd10)

!
!       (s px | py py):
!
        repij(2,8) = reploc(13)*r(3,1)*r32r32
     .             + reploc(14)*r(3,1)*radd04
     .             + reploc(15)*(r(1,1)*r12r32
     .                         + r(2,1)*r22r32)

        dxrepij(2,8) = dxreploc(13)*r(3,1)*r32r32
     .             +  reploc(13)*dxr(3,1)*r32r32
     .             +  reploc(13)*r(3,1)*dxr32r32
     .             + dxreploc(14)*r(3,1)*radd04
     .              + reploc(14)*dxr(3,1)*radd04
     .              + reploc(14)*r(3,1)*dxradd04
     .             + dxreploc(15)*(r(1,1)*r12r32
     .                         + r(2,1)*r22r32)
     .              + reploc(15)*(dxr(1,1)*r12r32
     .                          + dxr(2,1)*r22r32)
     .              + reploc(15)*(r(1,1)*dxr12r32
     .                          + r(2,1)*dxr22r32)

        dyrepij(2,8) = dyreploc(13)*r(3,1)*r32r32
     .             +  reploc(13)*dyr(3,1)*r32r32
     .             +  reploc(13)*r(3,1)*dyr32r32
     .             + dyreploc(14)*r(3,1)*radd04
     .              + reploc(14)*dyr(3,1)*radd04
     .              + reploc(14)*r(3,1)*dyradd04
     .             + dyreploc(15)*(r(1,1)*r12r32
     .                         + r(2,1)*r22r32)
     .              + reploc(15)*(dyr(1,1)*r12r32
     .                          + dyr(2,1)*r22r32)
     .              + reploc(15)*(r(1,1)*dyr12r32
     .                          + r(2,1)*dyr22r32)

        dzrepij(2,8) = dzreploc(13)*r(3,1)*r32r32
     .             +  reploc(13)*dzr(3,1)*r32r32
     .             +  reploc(13)*r(3,1)*dzr32r32
     .             + dzreploc(14)*r(3,1)*radd04
     .              + reploc(14)*dzr(3,1)*radd04
     .              + reploc(14)*r(3,1)*dzradd04
     .             + dzreploc(15)*(r(1,1)*r12r32
     .                         + r(2,1)*r22r32)
     .              + reploc(15)*(dzr(1,1)*r12r32
     .                          + dzr(2,1)*r22r32)
     .              + reploc(15)*(r(1,1)*dzr12r32
     .                          + r(2,1)*dzr22r32)
!
!       (py py | s px):
!
        repij(8,2) = reploc(8)*r(3,1)*r32r32
     .             + reploc(9)*r(3,1)*radd04
     .             + reploc(10)*(r(1,1)*r12r32
     .                         + r(2,1)*r22r32)

        dxrepij(8,2) = dxreploc(8)*r(3,1)*r32r32
     .              + reploc(8)*dxr(3,1)*r32r32
     .              + reploc(8)*r(3,1)*dxr32r32
     .             + dxreploc(9)*r(3,1)*radd04
     .             + reploc(9)*dxr(3,1)*radd04
     .             + reploc(9)*r(3,1)*dxradd04
     .             + dxreploc(10)*(r(1,1)*r12r32
     .                         + r(2,1)*r22r32)
     .             + reploc(10)*(dxr(1,1)*r12r32
     .                         + dxr(2,1)*r22r32)
     .             + reploc(10)*(r(1,1)*dxr12r32
     .                         + r(2,1)*dxr22r32)

        dyrepij(8,2) = dyreploc(8)*r(3,1)*r32r32
     .              + reploc(8)*dyr(3,1)*r32r32
     .              + reploc(8)*r(3,1)*dyr32r32
     .             + dyreploc(9)*r(3,1)*radd04
     .             + reploc(9)*dyr(3,1)*radd04
     .             + reploc(9)*r(3,1)*dyradd04
     .             + dyreploc(10)*(r(1,1)*r12r32
     .                         + r(2,1)*r22r32)
     .             + reploc(10)*(dyr(1,1)*r12r32
     .                         + dyr(2,1)*r22r32)
     .             + reploc(10)*(r(1,1)*dyr12r32
     .                         + r(2,1)*dyr22r32)

        dzrepij(8,2) = dzreploc(8)*r(3,1)*r32r32
     .              + reploc(8)*dzr(3,1)*r32r32
     .              + reploc(8)*r(3,1)*dzr32r32
     .             + dzreploc(9)*r(3,1)*radd04
     .             + reploc(9)*dzr(3,1)*radd04
     .             + reploc(9)*r(3,1)*dzradd04
     .             + dzreploc(10)*(r(1,1)*r12r32
     .                         + r(2,1)*r22r32)
     .             + reploc(10)*(dzr(1,1)*r12r32
     .                         + dzr(2,1)*r22r32)
     .             + reploc(10)*(r(1,1)*dzr12r32
     .                         + r(2,1)*dzr22r32)
!
!       (s px | py pz):
!
!
        repij(2,9) = reploc(13)*r(3,1)*r32r33
     .             + reploc(14)*r(3,1)*radd05
     .             + reploc(15)*(r(1,1)*radd13
     .                         + r(2,1)*radd14)

        dxrepij(2,9) = dxreploc(13)*r(3,1)*r32r33
     .             + reploc(13)*dxr(3,1)*r32r33
     .             + reploc(13)*r(3,1)*dxr32r33
     .             + dxreploc(14)*r(3,1)*radd05
     .             + reploc(14)*dxr(3,1)*radd05
     .             + reploc(14)*r(3,1)*dxradd05
     .             + dxreploc(15)*(r(1,1)*radd13
     .                         + r(2,1)*radd14)
     .             + reploc(15)*(dxr(1,1)*radd13
     .                         + dxr(2,1)*radd14)
     .             + reploc(15)*(r(1,1)*dxradd13
     .                         + r(2,1)*dxradd14)

        dyrepij(2,9) = dyreploc(13)*r(3,1)*r32r33
     .             + reploc(13)*dyr(3,1)*r32r33
     .             + reploc(13)*r(3,1)*dyr32r33
     .             + dyreploc(14)*r(3,1)*radd05
     .             + reploc(14)*dyr(3,1)*radd05
     .             + reploc(14)*r(3,1)*dyradd05
     .             + dyreploc(15)*(r(1,1)*radd13
     .                         + r(2,1)*radd14)
     .             + reploc(15)*(dyr(1,1)*radd13
     .                         + dyr(2,1)*radd14)
     .             + reploc(15)*(r(1,1)*dyradd13
     .                         + r(2,1)*dyradd14)

        dzrepij(2,9) = dzreploc(13)*r(3,1)*r32r33
     .             + reploc(13)*dzr(3,1)*r32r33
     .             + reploc(13)*r(3,1)*dzr32r33
     .             + dzreploc(14)*r(3,1)*radd05
     .             + reploc(14)*dzr(3,1)*radd05
     .             + reploc(14)*r(3,1)*dzradd05
     .             + dzreploc(15)*(r(1,1)*radd13
     .                         + r(2,1)*radd14)
     .             + reploc(15)*(dzr(1,1)*radd13
     .                         + dzr(2,1)*radd14)
     .             + reploc(15)*(r(1,1)*dzradd13
     .                         + r(2,1)*dzradd14)

!
!       (py pz | s px):
!
        repij(9,2) = reploc(8)*r(3,1)*r32r33
     .             + reploc(9)*r(3,1)*radd05
     .             + reploc(10)*(r(1,1)*radd13
     .                         + r(2,1)*radd14)

        dxrepij(9,2) = dxreploc(8)*r(3,1)*r32r33
     .             + reploc(8)*dxr(3,1)*r32r33
     .             + reploc(8)*r(3,1)*dxr32r33
     .             + dxreploc(9)*r(3,1)*radd05
     .             + reploc(9)*dxr(3,1)*radd05
     .             + reploc(9)*r(3,1)*dxradd05
     .             + dxreploc(10)*(r(1,1)*radd13
     .                         + r(2,1)*radd14)
     .             + reploc(10)*(dxr(1,1)*radd13
     .                         + dxr(2,1)*radd14)
     .             + reploc(10)*(r(1,1)*dxradd13
     .                         + r(2,1)*dxradd14)

        dyrepij(9,2) = dyreploc(8)*r(3,1)*r32r33
     .             + reploc(8)*dyr(3,1)*r32r33
     .             + reploc(8)*r(3,1)*dyr32r33
     .             + dyreploc(9)*r(3,1)*radd05
     .             + reploc(9)*dyr(3,1)*radd05
     .             + reploc(9)*r(3,1)*dyradd05
     .             + dyreploc(10)*(r(1,1)*radd13
     .                         + r(2,1)*radd14)
     .             + reploc(10)*(dyr(1,1)*radd13
     .                         + dyr(2,1)*radd14)
     .             + reploc(10)*(r(1,1)*dyradd13
     .                         + r(2,1)*dyradd14)

        dzrepij(9,2) = dzreploc(8)*r(3,1)*r32r33
     .             + reploc(8)*dzr(3,1)*r32r33
     .             + reploc(8)*r(3,1)*dzr32r33
     .             + dzreploc(9)*r(3,1)*radd05
     .             + reploc(9)*dzr(3,1)*radd05
     .             + reploc(9)*r(3,1)*dzradd05
     .             + dzreploc(10)*(r(1,1)*radd13
     .                         + r(2,1)*radd14)
     .             + reploc(10)*(dzr(1,1)*radd13
     .                         + dzr(2,1)*radd14)
     .             + reploc(10)*(r(1,1)*dzradd13
     .                         + r(2,1)*dzradd14)

!
!       (s px | pz pz):
!
        repij(2,10) = reploc(13)*r(3,1)*r33r33
     .              + reploc(14)*r(3,1)*radd06
     .              + reploc(15)*(r(1,1)*r13r33
     .                          + r(2,1)*r23r33)

        dxrepij(2,10) = dxreploc(13)*r(3,1)*r33r33
     .              + reploc(13)*dxr(3,1)*r33r33
     .              + reploc(13)*r(3,1)*dxr33r33
     .              + dxreploc(14)*r(3,1)*radd06
     .              + reploc(14)*dxr(3,1)*radd06
     .              + reploc(14)*r(3,1)*dxradd06
     .              + dxreploc(15)*(r(1,1)*r13r33
     .                          + r(2,1)*r23r33)
     .              + reploc(15)*(dxr(1,1)*r13r33
     .                          + dxr(2,1)*r23r33)
     .              + reploc(15)*(r(1,1)*dxr13r33
     .                          + r(2,1)*dxr23r33)

        dyrepij(2,10) = dyreploc(13)*r(3,1)*r33r33
     .              + reploc(13)*dyr(3,1)*r33r33
     .              + reploc(13)*r(3,1)*dyr33r33
     .              + dyreploc(14)*r(3,1)*radd06
     .              + reploc(14)*dyr(3,1)*radd06
     .              + reploc(14)*r(3,1)*dyradd06
     .              + dyreploc(15)*(r(1,1)*r13r33
     .                          + r(2,1)*r23r33)
     .              + reploc(15)*(dyr(1,1)*r13r33
     .                          + dyr(2,1)*r23r33)
     .              + reploc(15)*(r(1,1)*dyr13r33
     .                          + r(2,1)*dyr23r33)

        dzrepij(2,10) = dzreploc(13)*r(3,1)*r33r33
     .              + reploc(13)*dzr(3,1)*r33r33
     .              + reploc(13)*r(3,1)*dzr33r33
     .              + dzreploc(14)*r(3,1)*radd06
     .              + reploc(14)*dzr(3,1)*radd06
     .              + reploc(14)*r(3,1)*dzradd06
     .              + dzreploc(15)*(r(1,1)*r13r33
     .                          + r(2,1)*r23r33)
     .              + reploc(15)*(dzr(1,1)*r13r33
     .                          + dzr(2,1)*r23r33)
     .              + reploc(15)*(r(1,1)*dzr13r33
     .                          + r(2,1)*dzr23r33)
!
!       (pz pz | s px):
!
        repij(10,2) = reploc(8)*r(3,1)*r33r33
     .              + reploc(9)*r(3,1)*radd06
     .              + reploc(10)*(r(1,1)*r13r33
     .                          + r(2,1)*r23r33)

        dxrepij(10,2) = dxreploc(8)*r(3,1)*r33r33
     .              + reploc(8)*dxr(3,1)*r33r33
     .              + reploc(8)*r(3,1)*dxr33r33
     .              + dxreploc(9)*r(3,1)*radd06
     .              + reploc(9)*dxr(3,1)*radd06
     .              + reploc(9)*r(3,1)*dxradd06
     .              + dxreploc(10)*(r(1,1)*r13r33
     .                          + r(2,1)*r23r33)
     .              + reploc(10)*(dxr(1,1)*r13r33
     .                          + dxr(2,1)*r23r33)
     .              + reploc(10)*(r(1,1)*dxr13r33
     .                          + r(2,1)*dxr23r33)

        dyrepij(10,2) = dyreploc(8)*r(3,1)*r33r33
     .              + reploc(8)*dyr(3,1)*r33r33
     .              + reploc(8)*r(3,1)*dyr33r33
     .              + dyreploc(9)*r(3,1)*radd06
     .              + reploc(9)*dyr(3,1)*radd06
     .              + reploc(9)*r(3,1)*dyradd06
     .              + dyreploc(10)*(r(1,1)*r13r33
     .                          + r(2,1)*r23r33)
     .              + reploc(10)*(dyr(1,1)*r13r33
     .                          + dyr(2,1)*r23r33)
     .              + reploc(10)*(r(1,1)*dyr13r33
     .                          + r(2,1)*dyr23r33)

        dzrepij(10,2) = dzreploc(8)*r(3,1)*r33r33
     .              + reploc(8)*dzr(3,1)*r33r33
     .              + reploc(8)*r(3,1)*dzr33r33
     .              + dzreploc(9)*r(3,1)*radd06
     .              + reploc(9)*dzr(3,1)*radd06
     .              + reploc(9)*r(3,1)*dzradd06
     .              + dzreploc(10)*(r(1,1)*r13r33
     .                          + r(2,1)*r23r33)
     .              + reploc(10)*(dzr(1,1)*r13r33
     .                          + dzr(2,1)*r23r33)
     .              + reploc(10)*(r(1,1)*dzr13r33
     .                          + r(2,1)*dzr23r33)
!
!       (s py | s py):
!
        repij(3,3) = reploc(6)*r32r32
     .             + reploc(7)*radd04

         dxrepij(3,3) = dxreploc(6)*r32r32 + reploc(6)*dxr32r32
     .              + dxreploc(7)*radd04 + reploc(7)*dxradd04

         dyrepij(3,3) = dyreploc(6)*r32r32 + reploc(6)*dyr32r32
     .              + dyreploc(7)*radd04 + reploc(7)*dyradd04

         dzrepij(3,3) = dzreploc(6)*r32r32 + reploc(6)*dzr32r32
     .              + dzreploc(7)*radd04 + reploc(7)*dzradd04
!
!       (s py | s pz):
!
        repij(3,4) = reploc(6)*r32r33 
     .             + reploc(7)*radd05

        dxrepij(3,4) = dxreploc(6)*r32r33 + reploc(6)*dxr32r33
     .             + dxreploc(7)*radd05 + reploc(7)*dxradd05

        dyrepij(3,4) = dyreploc(6)*r32r33 + reploc(6)*dyr32r33
     .             + dyreploc(7)*radd05 + reploc(7)*dyradd05

        dzrepij(3,4) = dzreploc(6)*r32r33 + reploc(6)*dzr32r33
     .             + dzreploc(7)*radd05 + reploc(7)*dzradd05
!
!       (s pz | s py):
!
        repij(4,3) = repij(3,4)
         dxrepij(4,3) = dxrepij(3,4)
         dyrepij(4,3) = dyrepij(3,4)
         dzrepij(4,3) = dzrepij(3,4)
!
!       (s py | px px):
!
        repij(3,5) = reploc(13)*r(3,2)*r31r31
     .             + reploc(14)*r(3,2)*radd01
     .             + reploc(15)*(r(1,2)*r11r31
     .                         + r(2,2)*r21r31)

        dxrepij(3,5) = dxreploc(13)*r(3,2)*r31r31
     .             + reploc(13)*dxr(3,2)*r31r31
     .             + reploc(13)*r(3,2)*dxr31r31
     .             + dxreploc(14)*r(3,2)*radd01
     .             + reploc(14)*dxr(3,2)*radd01
     .             + reploc(14)*r(3,2)*dxradd01
     .             + dxreploc(15)*(r(1,2)*r11r31
     .                         + r(2,2)*r21r31)
     .             + reploc(15)*(dxr(1,2)*r11r31
     .                         + dxr(2,2)*r21r31)
     .             + reploc(15)*(r(1,2)*dxr11r31
     .                         + r(2,2)*dxr21r31)

        dyrepij(3,5) = dyreploc(13)*r(3,2)*r31r31
     .             + reploc(13)*dyr(3,2)*r31r31
     .             + reploc(13)*r(3,2)*dyr31r31
     .             + dyreploc(14)*r(3,2)*radd01
     .             + reploc(14)*dyr(3,2)*radd01
     .             + reploc(14)*r(3,2)*dyradd01
     .             + dyreploc(15)*(r(1,2)*r11r31
     .                         + r(2,2)*r21r31)
     .             + reploc(15)*(dyr(1,2)*r11r31
     .                         + dyr(2,2)*r21r31)
     .             + reploc(15)*(r(1,2)*dyr11r31
     .                         + r(2,2)*dyr21r31)

        dzrepij(3,5) = dzreploc(13)*r(3,2)*r31r31
     .             + reploc(13)*dzr(3,2)*r31r31
     .             + reploc(13)*r(3,2)*dzr31r31
     .             + dzreploc(14)*r(3,2)*radd01
     .             + reploc(14)*dzr(3,2)*radd01
     .             + reploc(14)*r(3,2)*dzradd01
     .             + dzreploc(15)*(r(1,2)*r11r31
     .                         + r(2,2)*r21r31)
     .             + reploc(15)*(dzr(1,2)*r11r31
     .                         + dzr(2,2)*r21r31)
     .             + reploc(15)*(r(1,2)*dzr11r31
     .                         + r(2,2)*dzr21r31)
!
!       (px px | s py):
!
        repij(5,3) = reploc(8)*r(3,2)*r31r31
     .             + reploc(9)*r(3,2)*radd01
     .             + reploc(10)*(r(1,2)*r11r31
     .                         + r(2,2)*r21r31)

        dxrepij(5,3) = dxreploc(8)*r(3,2)*r31r31
     .             + reploc(8)*dxr(3,2)*r31r31
     .             + reploc(8)*r(3,2)*dxr31r31
     .             + dxreploc(9)*r(3,2)*radd01
     .             + reploc(9)*dxr(3,2)*radd01
     .             + reploc(9)*r(3,2)*dxradd01
     .             + dxreploc(10)*(r(1,2)*r11r31
     .                         + r(2,2)*r21r31)
     .             + reploc(10)*(dxr(1,2)*r11r31
     .                         + dxr(2,2)*r21r31)
     .             + reploc(10)*(r(1,2)*dxr11r31
     .                         + r(2,2)*dxr21r31)

        dyrepij(5,3) = dyreploc(8)*r(3,2)*r31r31
     .             + reploc(8)*dyr(3,2)*r31r31
     .             + reploc(8)*r(3,2)*dyr31r31
     .             + dyreploc(9)*r(3,2)*radd01
     .             + reploc(9)*dyr(3,2)*radd01
     .             + reploc(9)*r(3,2)*dyradd01
     .             + dyreploc(10)*(r(1,2)*r11r31
     .                         + r(2,2)*r21r31)
     .             + reploc(10)*(dyr(1,2)*r11r31
     .                         + dyr(2,2)*r21r31)
     .             + reploc(10)*(r(1,2)*dyr11r31
     .                         + r(2,2)*dyr21r31)

        dzrepij(5,3) = dzreploc(8)*r(3,2)*r31r31
     .             + reploc(8)*dzr(3,2)*r31r31
     .             + reploc(8)*r(3,2)*dzr31r31
     .             + dzreploc(9)*r(3,2)*radd01
     .             + reploc(9)*dzr(3,2)*radd01
     .             + reploc(9)*r(3,2)*dzradd01
     .             + dzreploc(10)*(r(1,2)*r11r31
     .                         + r(2,2)*r21r31)
     .             + reploc(10)*(dzr(1,2)*r11r31
     .                         + dzr(2,2)*r21r31)
     .             + reploc(10)*(r(1,2)*dzr11r31
     .                         + r(2,2)*dzr21r31)
!
!       (s py | px py):
!
        repij(3,6) = reploc(13)*r(3,2)*r31r32
     .             + reploc(14)*r(3,2)*radd02
     .             + reploc(15)*(r(1,2)*radd11
     .                         + r(2,2)*radd12)

        dxrepij(3,6) = dxreploc(13)*r(3,2)*r31r32
     .             + reploc(13)*dxr(3,2)*r31r32
     .             + reploc(13)*r(3,2)*dxr31r32
     .             + dxreploc(14)*r(3,2)*radd02
     .             + reploc(14)*dxr(3,2)*radd02
     .             + reploc(14)*r(3,2)*dxradd02
     .             + dxreploc(15)*(r(1,2)*radd11
     .                         + r(2,2)*radd12)
     .             + reploc(15)*(dxr(1,2)*radd11
     .                         + dxr(2,2)*radd12)
     .             + reploc(15)*(r(1,2)*dxradd11
     .                         + r(2,2)*dxradd12)

        dyrepij(3,6) = dyreploc(13)*r(3,2)*r31r32
     .             + reploc(13)*dyr(3,2)*r31r32
     .             + reploc(13)*r(3,2)*dyr31r32
     .             + dyreploc(14)*r(3,2)*radd02
     .             + reploc(14)*dyr(3,2)*radd02
     .             + reploc(14)*r(3,2)*dyradd02
     .             + dyreploc(15)*(r(1,2)*radd11
     .                         + r(2,2)*radd12)
     .             + reploc(15)*(dyr(1,2)*radd11
     .                         + dyr(2,2)*radd12)
     .             + reploc(15)*(r(1,2)*dyradd11
     .                         + r(2,2)*dyradd12)

        dzrepij(3,6) = dzreploc(13)*r(3,2)*r31r32
     .             + reploc(13)*dzr(3,2)*r31r32
     .             + reploc(13)*r(3,2)*dzr31r32
     .             + dzreploc(14)*r(3,2)*radd02
     .             + reploc(14)*dzr(3,2)*radd02
     .             + reploc(14)*r(3,2)*dzradd02
     .             + dzreploc(15)*(r(1,2)*radd11
     .                         + r(2,2)*radd12)
     .             + reploc(15)*(dzr(1,2)*radd11
     .                         + dzr(2,2)*radd12)
     .             + reploc(15)*(r(1,2)*dzradd11
     .                         + r(2,2)*dzradd12)

!
!       (px py | s py):
!
        repij(6,3) = reploc(8)*r(3,2)*r31r32
     .             + reploc(9)*r(3,2)*radd02
     .             + reploc(10)*(r(1,2)*radd11
     .                         + r(2,2)*radd12)

        dxrepij(6,3) = dxreploc(8)*r(3,2)*r31r32
     .             + reploc(8)*dxr(3,2)*r31r32
     .             + reploc(8)*r(3,2)*dxr31r32
     .             + dxreploc(9)*r(3,2)*radd02
     .             + reploc(9)*dxr(3,2)*radd02
     .             + reploc(9)*r(3,2)*dxradd02
     .             + dxreploc(10)*(r(1,2)*radd11
     .                         + r(2,2)*radd12)
     .             + reploc(10)*(dxr(1,2)*radd11
     .                         + dxr(2,2)*radd12)
     .             + reploc(10)*(r(1,2)*dxradd11
     .                         + r(2,2)*dxradd12)

        dyrepij(6,3) = dyreploc(8)*r(3,2)*r31r32
     .             + reploc(8)*dyr(3,2)*r31r32
     .             + reploc(8)*r(3,2)*dyr31r32
     .             + dyreploc(9)*r(3,2)*radd02
     .             + reploc(9)*dyr(3,2)*radd02
     .             + reploc(9)*r(3,2)*dyradd02
     .             + dyreploc(10)*(r(1,2)*radd11
     .                         + r(2,2)*radd12)
     .             + reploc(10)*(dyr(1,2)*radd11
     .                         + dyr(2,2)*radd12)
     .             + reploc(10)*(r(1,2)*dyradd11
     .                         + r(2,2)*dyradd12)

        dzrepij(6,3) = dzreploc(8)*r(3,2)*r31r32
     .             + reploc(8)*dzr(3,2)*r31r32
     .             + reploc(8)*r(3,2)*dzr31r32
     .             + dzreploc(9)*r(3,2)*radd02
     .             + reploc(9)*dzr(3,2)*radd02
     .             + reploc(9)*r(3,2)*dzradd02
     .             + dzreploc(10)*(r(1,2)*radd11
     .                         + r(2,2)*radd12)
     .             + reploc(10)*(dzr(1,2)*radd11
     .                         + dzr(2,2)*radd12)
     .             + reploc(10)*(r(1,2)*dzradd11
     .                         + r(2,2)*dzradd12)
!
!       (s py | px pz):
!
        repij(3,7) = reploc(13)*r(3,2)*r31r33
     .             + reploc(14)*r(3,2)*radd03
     .             + reploc(15)*(r(1,2)*radd09
     .                         + r(2,2)*radd10)

        dxrepij(3,7) = dxreploc(13)*r(3,2)*r31r33
     .             + reploc(13)*dxr(3,2)*r31r33
     .             + reploc(13)*r(3,2)*dxr31r33
     .             + dxreploc(14)*r(3,2)*radd03
     .             + reploc(14)*dxr(3,2)*radd03
     .             + reploc(14)*r(3,2)*dxradd03
     .             + dxreploc(15)*(r(1,2)*radd09
     .                         + r(2,2)*radd10)
     .             + reploc(15)*(dxr(1,2)*radd09
     .                         + dxr(2,2)*radd10)
     .             + reploc(15)*(r(1,2)*dxradd09
     .                         + r(2,2)*dxradd10)

        dyrepij(3,7) = dyreploc(13)*r(3,2)*r31r33
     .             + reploc(13)*dyr(3,2)*r31r33
     .             + reploc(13)*r(3,2)*dyr31r33
     .             + dyreploc(14)*r(3,2)*radd03
     .             + reploc(14)*dyr(3,2)*radd03
     .             + reploc(14)*r(3,2)*dyradd03
     .             + dyreploc(15)*(r(1,2)*radd09
     .                         + r(2,2)*radd10)
     .             + reploc(15)*(dyr(1,2)*radd09
     .                         + dyr(2,2)*radd10)
     .             + reploc(15)*(r(1,2)*dyradd09
     .                         + r(2,2)*dyradd10)

        dzrepij(3,7) = dzreploc(13)*r(3,2)*r31r33
     .             + reploc(13)*dzr(3,2)*r31r33
     .             + reploc(13)*r(3,2)*dzr31r33
     .             + dzreploc(14)*r(3,2)*radd03
     .             + reploc(14)*dzr(3,2)*radd03
     .             + reploc(14)*r(3,2)*dzradd03
     .             + dzreploc(15)*(r(1,2)*radd09
     .                         + r(2,2)*radd10)
     .             + reploc(15)*(dzr(1,2)*radd09
     .                         + dzr(2,2)*radd10)
     .             + reploc(15)*(r(1,2)*dzradd09
     .                         + r(2,2)*dzradd10)
!
!       (px pz | s py):
!
        repij(7,3) = reploc(8)*r(3,2)*r31r33
     .             + reploc(9)*r(3,2)*radd03
     .             + reploc(10)*(r(1,2)*radd09
     .                         + r(2,2)*radd10)

        dxrepij(7,3) = dxreploc(8)*r(3,2)*r31r33
     .             + reploc(8)*dxr(3,2)*r31r33
     .             + reploc(8)*r(3,2)*dxr31r33
     .             + dxreploc(9)*r(3,2)*radd03
     .             + reploc(9)*dxr(3,2)*radd03
     .             + reploc(9)*r(3,2)*dxradd03
     .             + dxreploc(10)*(r(1,2)*radd09
     .                         + r(2,2)*radd10)
     .             + reploc(10)*(dxr(1,2)*radd09
     .                         + dxr(2,2)*radd10)
     .             + reploc(10)*(r(1,2)*dxradd09
     .                         + r(2,2)*dxradd10)

        dyrepij(7,3) = dyreploc(8)*r(3,2)*r31r33
     .             + reploc(8)*dyr(3,2)*r31r33
     .             + reploc(8)*r(3,2)*dyr31r33
     .             + dyreploc(9)*r(3,2)*radd03
     .             + reploc(9)*dyr(3,2)*radd03
     .             + reploc(9)*r(3,2)*dyradd03
     .             + dyreploc(10)*(r(1,2)*radd09
     .                         + r(2,2)*radd10)
     .             + reploc(10)*(dyr(1,2)*radd09
     .                         + dyr(2,2)*radd10)
     .             + reploc(10)*(r(1,2)*dyradd09
     .                         + r(2,2)*dyradd10)

        dzrepij(7,3) = dzreploc(8)*r(3,2)*r31r33
     .             + reploc(8)*dzr(3,2)*r31r33
     .             + reploc(8)*r(3,2)*dzr31r33
     .             + dzreploc(9)*r(3,2)*radd03
     .             + reploc(9)*dzr(3,2)*radd03
     .             + reploc(9)*r(3,2)*dzradd03
     .             + dzreploc(10)*(r(1,2)*radd09
     .                         + r(2,2)*radd10)
     .             + reploc(10)*(dzr(1,2)*radd09
     .                         + dzr(2,2)*radd10)
     .             + reploc(10)*(r(1,2)*dzradd09
     .                         + r(2,2)*dzradd10)
!
!       (s py | py py):
!
        repij(3,8) = reploc(13)*r(3,2)*r32r32
     .             + reploc(14)*r(3,2)*radd04
     .             + reploc(15)*(r(1,2)*r12r32
     .                         + r(2,2)*r22r32)

        dxrepij(3,8) = dxreploc(13)*r(3,2)*r32r32
     .             + reploc(13)*dxr(3,2)*r32r32
     .             + reploc(13)*r(3,2)*dxr32r32
     .             + dxreploc(14)*r(3,2)*radd04
     .             + reploc(14)*dxr(3,2)*radd04
     .             + reploc(14)*r(3,2)*dxradd04
     .             + dxreploc(15)*(r(1,2)*r12r32
     .                         + r(2,2)*r22r32)
     .             + reploc(15)*(dxr(1,2)*r12r32
     .                         + dxr(2,2)*r22r32)
     .             + reploc(15)*(r(1,2)*dxr12r32
     .                         + r(2,2)*dxr22r32)

        dyrepij(3,8) = dyreploc(13)*r(3,2)*r32r32
     .             + reploc(13)*dyr(3,2)*r32r32
     .             + reploc(13)*r(3,2)*dyr32r32
     .             + dyreploc(14)*r(3,2)*radd04
     .             + reploc(14)*dyr(3,2)*radd04
     .             + reploc(14)*r(3,2)*dyradd04
     .             + dyreploc(15)*(r(1,2)*r12r32
     .                         + r(2,2)*r22r32)
     .             + reploc(15)*(dyr(1,2)*r12r32
     .                         + dyr(2,2)*r22r32)
     .             + reploc(15)*(r(1,2)*dyr12r32
     .                         + r(2,2)*dyr22r32)

        dzrepij(3,8) = dzreploc(13)*r(3,2)*r32r32
     .             + reploc(13)*dzr(3,2)*r32r32
     .             + reploc(13)*r(3,2)*dzr32r32
     .             + dzreploc(14)*r(3,2)*radd04
     .             + reploc(14)*dzr(3,2)*radd04
     .             + reploc(14)*r(3,2)*dzradd04
     .             + dzreploc(15)*(r(1,2)*r12r32
     .                         + r(2,2)*r22r32)
     .             + reploc(15)*(dzr(1,2)*r12r32
     .                         + dzr(2,2)*r22r32)
     .             + reploc(15)*(r(1,2)*dzr12r32
     .                         + r(2,2)*dzr22r32)
!
!       (py py | s py):
!
        repij(8,3) = reploc(8)*r(3,2)*r32r32
     .             + reploc(9)*r(3,2)*radd04
     .             + reploc(10)*(r(1,2)*r12r32
     .                         + r(2,2)*r22r32)

        dxrepij(8,3) = dxreploc(8)*r(3,2)*r32r32
     .             + reploc(8)*dxr(3,2)*r32r32
     .             + reploc(8)*r(3,2)*dxr32r32
     .             + dxreploc(9)*r(3,2)*radd04
     .             + reploc(9)*dxr(3,2)*radd04
     .             + reploc(9)*r(3,2)*dxradd04
     .             + dxreploc(10)*(r(1,2)*r12r32
     .                         + r(2,2)*r22r32)
     .             + reploc(10)*(dxr(1,2)*r12r32
     .                         + dxr(2,2)*r22r32)
     .             + reploc(10)*(r(1,2)*dxr12r32
     .                         + r(2,2)*dxr22r32)

        dyrepij(8,3) = dyreploc(8)*r(3,2)*r32r32
     .             + reploc(8)*dyr(3,2)*r32r32
     .             + reploc(8)*r(3,2)*dyr32r32
     .             + dyreploc(9)*r(3,2)*radd04
     .             + reploc(9)*dyr(3,2)*radd04
     .             + reploc(9)*r(3,2)*dyradd04
     .             + dyreploc(10)*(r(1,2)*r12r32
     .                         + r(2,2)*r22r32)
     .             + reploc(10)*(dyr(1,2)*r12r32
     .                         + dyr(2,2)*r22r32)
     .             + reploc(10)*(r(1,2)*dyr12r32
     .                         + r(2,2)*dyr22r32)

        dzrepij(8,3) = dzreploc(8)*r(3,2)*r32r32
     .             + reploc(8)*dzr(3,2)*r32r32
     .             + reploc(8)*r(3,2)*dzr32r32
     .             + dzreploc(9)*r(3,2)*radd04
     .             + reploc(9)*dzr(3,2)*radd04
     .             + reploc(9)*r(3,2)*dzradd04
     .             + dzreploc(10)*(r(1,2)*r12r32
     .                         + r(2,2)*r22r32)
     .             + reploc(10)*(dzr(1,2)*r12r32
     .                         + dzr(2,2)*r22r32)
     .             + reploc(10)*(r(1,2)*dzr12r32
     .                         + r(2,2)*dzr22r32)
!
!       (s py | py pz):
!
        repij(3,9) = reploc(13)*r(3,2)*r32r33
     .             + reploc(14)*r(3,2)*radd05
     .             + reploc(15)*(r(1,2)*radd13
     .                         + r(2,2)*radd14)

        dxrepij(3,9) = dxreploc(13)*r(3,2)*r32r33
     .             + reploc(13)*dxr(3,2)*r32r33
     .             + reploc(13)*r(3,2)*dxr32r33
     .             + dxreploc(14)*r(3,2)*radd05
     .             + reploc(14)*dxr(3,2)*radd05
     .             + reploc(14)*r(3,2)*dxradd05
     .             + dxreploc(15)*(r(1,2)*radd13
     .                         + r(2,2)*radd14)
     .             + reploc(15)*(dxr(1,2)*radd13
     .                         + dxr(2,2)*radd14)
     .             + reploc(15)*(r(1,2)*dxradd13
     .                         + r(2,2)*dxradd14)

        dyrepij(3,9) = dyreploc(13)*r(3,2)*r32r33
     .             + reploc(13)*dyr(3,2)*r32r33
     .             + reploc(13)*r(3,2)*dyr32r33
     .             + dyreploc(14)*r(3,2)*radd05
     .             + reploc(14)*dyr(3,2)*radd05
     .             + reploc(14)*r(3,2)*dyradd05
     .             + dyreploc(15)*(r(1,2)*radd13
     .                         + r(2,2)*radd14)
     .             + reploc(15)*(dyr(1,2)*radd13
     .                         + dyr(2,2)*radd14)
     .             + reploc(15)*(r(1,2)*dyradd13
     .                         + r(2,2)*dyradd14)

        dzrepij(3,9) = dzreploc(13)*r(3,2)*r32r33
     .             + reploc(13)*dzr(3,2)*r32r33
     .             + reploc(13)*r(3,2)*dzr32r33
     .             + dzreploc(14)*r(3,2)*radd05
     .             + reploc(14)*dzr(3,2)*radd05
     .             + reploc(14)*r(3,2)*dzradd05
     .             + dzreploc(15)*(r(1,2)*radd13
     .                         + r(2,2)*radd14)
     .             + reploc(15)*(dzr(1,2)*radd13
     .                         + dzr(2,2)*radd14)
     .             + reploc(15)*(r(1,2)*dzradd13
     .                         + r(2,2)*dzradd14)

!
!       (py pz | s py):
!
        repij(9,3) = reploc(8)*r(3,2)*r32r33
     .             + reploc(9)*r(3,2)*radd05
     .             + reploc(10)*(r(1,2)*radd13
     .                         + r(2,2)*radd14)

        dxrepij(9,3) = dxreploc(8)*r(3,2)*r32r33
     .             + reploc(8)*dxr(3,2)*r32r33
     .             + reploc(8)*r(3,2)*dxr32r33
     .             + dxreploc(9)*r(3,2)*radd05
     .             + reploc(9)*dxr(3,2)*radd05
     .             + reploc(9)*r(3,2)*dxradd05
     .             + dxreploc(10)*(r(1,2)*radd13
     .                         + r(2,2)*radd14)
     .             + reploc(10)*(dxr(1,2)*radd13
     .                         + dxr(2,2)*radd14)
     .             + reploc(10)*(r(1,2)*dxradd13
     .                         + r(2,2)*dxradd14)

        dyrepij(9,3) = dyreploc(8)*r(3,2)*r32r33
     .             + reploc(8)*dyr(3,2)*r32r33
     .             + reploc(8)*r(3,2)*dyr32r33
     .             + dyreploc(9)*r(3,2)*radd05
     .             + reploc(9)*dyr(3,2)*radd05
     .             + reploc(9)*r(3,2)*dyradd05
     .             + dyreploc(10)*(r(1,2)*radd13
     .                         + r(2,2)*radd14)
     .             + reploc(10)*(dyr(1,2)*radd13
     .                         + dyr(2,2)*radd14)
     .             + reploc(10)*(r(1,2)*dyradd13
     .                         + r(2,2)*dyradd14)

        dzrepij(9,3) = dzreploc(8)*r(3,2)*r32r33
     .             + reploc(8)*dzr(3,2)*r32r33
     .             + reploc(8)*r(3,2)*dzr32r33
     .             + dzreploc(9)*r(3,2)*radd05
     .             + reploc(9)*dzr(3,2)*radd05
     .             + reploc(9)*r(3,2)*dzradd05
     .             + dzreploc(10)*(r(1,2)*radd13
     .                         + r(2,2)*radd14)
     .             + reploc(10)*(dzr(1,2)*radd13
     .                         + dzr(2,2)*radd14)
     .             + reploc(10)*(r(1,2)*dzradd13
     .                         + r(2,2)*dzradd14)

!
!       (s py | pz pz):
!
        repij(3,10) = reploc(13)*r(3,2)*r33r33
     .              + reploc(14)*r(3,2)*radd06
     .              + reploc(15)*(r(1,2)*r13r33
     .                          + r(2,2)*r23r33)

        dxrepij(3,10) = dxreploc(13)*r(3,2)*r33r33
     .              + reploc(13)*dxr(3,2)*r33r33
     .              + reploc(13)*r(3,2)*dxr33r33
     .              + dxreploc(14)*r(3,2)*radd06
     .              + reploc(14)*dxr(3,2)*radd06
     .              + reploc(14)*r(3,2)*dxradd06
     .              + dxreploc(15)*(r(1,2)*r13r33
     .                          + r(2,2)*r23r33)
     .              + reploc(15)*(dxr(1,2)*r13r33
     .                          + dxr(2,2)*r23r33)
     .              + reploc(15)*(r(1,2)*dxr13r33
     .                          + r(2,2)*dxr23r33)

        dyrepij(3,10) = dyreploc(13)*r(3,2)*r33r33
     .              + reploc(13)*dyr(3,2)*r33r33
     .              + reploc(13)*r(3,2)*dyr33r33
     .              + dyreploc(14)*r(3,2)*radd06
     .              + reploc(14)*dyr(3,2)*radd06
     .              + reploc(14)*r(3,2)*dyradd06
     .              + dyreploc(15)*(r(1,2)*r13r33
     .                          + r(2,2)*r23r33)
     .              + reploc(15)*(dyr(1,2)*r13r33
     .                          + dyr(2,2)*r23r33)
     .              + reploc(15)*(r(1,2)*dyr13r33
     .                          + r(2,2)*dyr23r33)

        dzrepij(3,10) = dzreploc(13)*r(3,2)*r33r33
     .              + reploc(13)*dzr(3,2)*r33r33
     .              + reploc(13)*r(3,2)*dzr33r33
     .              + dzreploc(14)*r(3,2)*radd06
     .              + reploc(14)*dzr(3,2)*radd06
     .              + reploc(14)*r(3,2)*dzradd06
     .              + dzreploc(15)*(r(1,2)*r13r33
     .                          + r(2,2)*r23r33)
     .              + reploc(15)*(dzr(1,2)*r13r33
     .                          + dzr(2,2)*r23r33)
     .              + reploc(15)*(r(1,2)*dzr13r33
     .                          + r(2,2)*dzr23r33)
!
!       (pz pz | s py):
!
        repij(10,3) = reploc(8)*r(3,2)*r33r33
     .              + reploc(9)*r(3,2)*radd06
     .              + reploc(10)*(r(1,2)*r13r33
     .                          + r(2,2)*r23r33)

        dxrepij(10,3) = dxreploc(8)*r(3,2)*r33r33
     .              + reploc(8)*dxr(3,2)*r33r33
     .              + reploc(8)*r(3,2)*dxr33r33
     .              + dxreploc(9)*r(3,2)*radd06
     .              + reploc(9)*dxr(3,2)*radd06
     .              + reploc(9)*r(3,2)*dxradd06
     .              + dxreploc(10)*(r(1,2)*r13r33
     .                          + r(2,2)*r23r33)
     .              + reploc(10)*(dxr(1,2)*r13r33
     .                          + dxr(2,2)*r23r33)
     .              + reploc(10)*(r(1,2)*dxr13r33
     .                          + r(2,2)*dxr23r33)

        dyrepij(10,3) = dyreploc(8)*r(3,2)*r33r33
     .              + reploc(8)*dyr(3,2)*r33r33
     .              + reploc(8)*r(3,2)*dyr33r33
     .              + dyreploc(9)*r(3,2)*radd06
     .              + reploc(9)*dyr(3,2)*radd06
     .              + reploc(9)*r(3,2)*dyradd06
     .              + dyreploc(10)*(r(1,2)*r13r33
     .                          + r(2,2)*r23r33)
     .              + reploc(10)*(dyr(1,2)*r13r33
     .                          + dyr(2,2)*r23r33)
     .              + reploc(10)*(r(1,2)*dyr13r33
     .                          + r(2,2)*dyr23r33)

        dzrepij(10,3) = dzreploc(8)*r(3,2)*r33r33
     .              + reploc(8)*dzr(3,2)*r33r33
     .              + reploc(8)*r(3,2)*dzr33r33
     .              + dzreploc(9)*r(3,2)*radd06
     .              + reploc(9)*dzr(3,2)*radd06
     .              + reploc(9)*r(3,2)*dzradd06
     .              + dzreploc(10)*(r(1,2)*r13r33
     .                          + r(2,2)*r23r33)
     .              + reploc(10)*(dzr(1,2)*r13r33
     .                          + dzr(2,2)*r23r33)
     .              + reploc(10)*(r(1,2)*dzr13r33
     .                          + r(2,2)*dzr23r33)

!
!       (s pz | s pz):
!
        repij(4,4) = reploc(6)*r33r33
     .             + reploc(7)*radd06

        dxrepij(4,4) = dxreploc(6)*r33r33 + reploc(6)*dxr33r33
     .             + dxreploc(7)*radd06 + reploc(7)*dxradd06

        dyrepij(4,4) = dyreploc(6)*r33r33 + reploc(6)*dyr33r33
     .             + dyreploc(7)*radd06 + reploc(7)*dyradd06

        dzrepij(4,4) = dzreploc(6)*r33r33 + reploc(6)*dzr33r33
     .             + dzreploc(7)*radd06 + reploc(7)*dzradd06
!
!       (s pz | px px):
!
        repij(4,5) = reploc(13)*r(3,3)*r31r31
     .             + reploc(14)*r(3,3)*radd01
     .             + reploc(15)*(r(1,3)*r11r31
     .                         + r(2,3)*r21r31)

        dxrepij(4,5) = dxreploc(13)*r(3,3)*r31r31
     .             + reploc(13)*dxr(3,3)*r31r31
     .             + reploc(13)*r(3,3)*dxr31r31
     .             + dxreploc(14)*r(3,3)*radd01
     .             + reploc(14)*dxr(3,3)*radd01
     .             + reploc(14)*r(3,3)*dxradd01
     .             + dxreploc(15)*(r(1,3)*r11r31
     .                         + r(2,3)*r21r31)
     .             + reploc(15)*(dxr(1,3)*r11r31
     .                         + dxr(2,3)*r21r31)
     .             + reploc(15)*(r(1,3)*dxr11r31
     .                         + r(2,3)*dxr21r31)

        dyrepij(4,5) = dyreploc(13)*r(3,3)*r31r31
     .             + reploc(13)*dyr(3,3)*r31r31
     .             + reploc(13)*r(3,3)*dyr31r31
     .             + dyreploc(14)*r(3,3)*radd01
     .             + reploc(14)*dyr(3,3)*radd01
     .             + reploc(14)*r(3,3)*dyradd01
     .             + dyreploc(15)*(r(1,3)*r11r31
     .                         + r(2,3)*r21r31)
     .             + reploc(15)*(dyr(1,3)*r11r31
     .                         + dyr(2,3)*r21r31)
     .             + reploc(15)*(r(1,3)*dyr11r31
     .                         + r(2,3)*dyr21r31)

        dzrepij(4,5) = dzreploc(13)*r(3,3)*r31r31
     .             + reploc(13)*dzr(3,3)*r31r31
     .             + reploc(13)*r(3,3)*dzr31r31
     .             + dzreploc(14)*r(3,3)*radd01
     .             + reploc(14)*dzr(3,3)*radd01
     .             + reploc(14)*r(3,3)*dzradd01
     .             + dzreploc(15)*(r(1,3)*r11r31
     .                         + r(2,3)*r21r31)
     .             + reploc(15)*(dzr(1,3)*r11r31
     .                         + dzr(2,3)*r21r31)
     .             + reploc(15)*(r(1,3)*dzr11r31
     .                         + r(2,3)*dzr21r31)
!
!       (px px | s pz):
!
        repij(5,4) = reploc(8)*r(3,3)*r31r31
     .             + reploc(9)*r(3,3)*radd01
     .             + reploc(10)*(r(1,3)*r11r31
     .                         + r(2,3)*r21r31)

        dxrepij(5,4) = dxreploc(8)*r(3,3)*r31r31
     .             + reploc(8)*dxr(3,3)*r31r31
     .             + reploc(8)*r(3,3)*dxr31r31
     .             + dxreploc(9)*r(3,3)*radd01
     .             + reploc(9)*dxr(3,3)*radd01
     .             + reploc(9)*r(3,3)*dxradd01
     .             + dxreploc(10)*(r(1,3)*r11r31
     .                         + r(2,3)*r21r31)
     .             + reploc(10)*(dxr(1,3)*r11r31
     .                         + dxr(2,3)*r21r31)
     .             + reploc(10)*(r(1,3)*dxr11r31
     .                         + r(2,3)*dxr21r31)

        dyrepij(5,4) = dyreploc(8)*r(3,3)*r31r31
     .             + reploc(8)*dyr(3,3)*r31r31
     .             + reploc(8)*r(3,3)*dyr31r31
     .             + dyreploc(9)*r(3,3)*radd01
     .             + reploc(9)*dyr(3,3)*radd01
     .             + reploc(9)*r(3,3)*dyradd01
     .             + dyreploc(10)*(r(1,3)*r11r31
     .                         + r(2,3)*r21r31)
     .             + reploc(10)*(dyr(1,3)*r11r31
     .                         + dyr(2,3)*r21r31)
     .             + reploc(10)*(r(1,3)*dyr11r31
     .                         + r(2,3)*dyr21r31)

        dzrepij(5,4) = dzreploc(8)*r(3,3)*r31r31
     .             + reploc(8)*dzr(3,3)*r31r31
     .             + reploc(8)*r(3,3)*dzr31r31
     .             + dzreploc(9)*r(3,3)*radd01
     .             + reploc(9)*dzr(3,3)*radd01
     .             + reploc(9)*r(3,3)*dzradd01
     .             + dzreploc(10)*(r(1,3)*r11r31
     .                         + r(2,3)*r21r31)
     .             + reploc(10)*(dzr(1,3)*r11r31
     .                         + dzr(2,3)*r21r31)
     .             + reploc(10)*(r(1,3)*dzr11r31
     .                         + r(2,3)*dzr21r31)
!
!       (s pz | px py):
!
        repij(4,6) = reploc(13)*r(3,3)*r31r32
     .             + reploc(14)*r(3,3)*radd02
     .             + reploc(15)*(r(1,3)*radd11
     .                         + r(2,3)*radd12)

        dxrepij(4,6) = dxreploc(13)*r(3,3)*r31r32
     .             + reploc(13)*dxr(3,3)*r31r32
     .             + reploc(13)*r(3,3)*dxr31r32
     .             + dxreploc(14)*r(3,3)*radd02
     .             + reploc(14)*dxr(3,3)*radd02
     .             + reploc(14)*r(3,3)*dxradd02
     .             + dxreploc(15)*(r(1,3)*radd11
     .                         + r(2,3)*radd12)
     .             + reploc(15)*(dxr(1,3)*radd11
     .                         + dxr(2,3)*radd12)
     .             + reploc(15)*(r(1,3)*dxradd11
     .                         + r(2,3)*dxradd12)

        dyrepij(4,6) = dyreploc(13)*r(3,3)*r31r32
     .             + reploc(13)*dyr(3,3)*r31r32
     .             + reploc(13)*r(3,3)*dyr31r32
     .             + dyreploc(14)*r(3,3)*radd02
     .             + reploc(14)*dyr(3,3)*radd02
     .             + reploc(14)*r(3,3)*dyradd02
     .             + dyreploc(15)*(r(1,3)*radd11
     .                         + r(2,3)*radd12)
     .             + reploc(15)*(dyr(1,3)*radd11
     .                         + dyr(2,3)*radd12)
     .             + reploc(15)*(r(1,3)*dyradd11
     .                         + r(2,3)*dyradd12)

        dzrepij(4,6) = dzreploc(13)*r(3,3)*r31r32
     .             + reploc(13)*dzr(3,3)*r31r32
     .             + reploc(13)*r(3,3)*dzr31r32
     .             + dzreploc(14)*r(3,3)*radd02
     .             + reploc(14)*dzr(3,3)*radd02
     .             + reploc(14)*r(3,3)*dzradd02
     .             + dzreploc(15)*(r(1,3)*radd11
     .                         + r(2,3)*radd12)
     .             + reploc(15)*(dzr(1,3)*radd11
     .                         + dzr(2,3)*radd12)
     .             + reploc(15)*(r(1,3)*dzradd11
     .                         + r(2,3)*dzradd12)
!
!       (px py | s pz):
!
        repij(6,4) = reploc(8)*r(3,3)*r31r32
     .             + reploc(9)*r(3,3)*radd02
     .             + reploc(10)*(r(1,3)*radd11
     .                         + r(2,3)*radd12)

        dxrepij(6,4) = dxreploc(8)*r(3,3)*r31r32
     .             + reploc(8)*dxr(3,3)*r31r32
     .             + reploc(8)*r(3,3)*dxr31r32
     .             + dxreploc(9)*r(3,3)*radd02
     .             + reploc(9)*dxr(3,3)*radd02
     .             + reploc(9)*r(3,3)*dxradd02
     .             + dxreploc(10)*(r(1,3)*radd11
     .                         + r(2,3)*radd12)
     .             + reploc(10)*(dxr(1,3)*radd11
     .                         + dxr(2,3)*radd12)
     .             + reploc(10)*(r(1,3)*dxradd11
     .                         + r(2,3)*dxradd12)

        dyrepij(6,4) = dyreploc(8)*r(3,3)*r31r32
     .             + reploc(8)*dyr(3,3)*r31r32
     .             + reploc(8)*r(3,3)*dyr31r32
     .             + dyreploc(9)*r(3,3)*radd02
     .             + reploc(9)*dyr(3,3)*radd02
     .             + reploc(9)*r(3,3)*dyradd02
     .             + dyreploc(10)*(r(1,3)*radd11
     .                         + r(2,3)*radd12)
     .             + reploc(10)*(dyr(1,3)*radd11
     .                         + dyr(2,3)*radd12)
     .             + reploc(10)*(r(1,3)*dyradd11
     .                         + r(2,3)*dyradd12)

        dzrepij(6,4) = dzreploc(8)*r(3,3)*r31r32
     .             + reploc(8)*dzr(3,3)*r31r32
     .             + reploc(8)*r(3,3)*dzr31r32
     .             + dzreploc(9)*r(3,3)*radd02
     .             + reploc(9)*dzr(3,3)*radd02
     .             + reploc(9)*r(3,3)*dzradd02
     .             + dzreploc(10)*(r(1,3)*radd11
     .                         + r(2,3)*radd12)
     .             + reploc(10)*(dzr(1,3)*radd11
     .                         + dzr(2,3)*radd12)
     .             + reploc(10)*(r(1,3)*dzradd11
     .                         + r(2,3)*dzradd12)
!
!       (s pz | px pz):
!
        repij(4,7) = reploc(13)*r(3,3)*r31r33
     .             + reploc(14)*r(3,3)*radd03
     .             + reploc(15)*(r(1,3)*radd09
     .                         + r(2,3)*radd10)

        dxrepij(4,7) = dxreploc(13)*r(3,3)*r31r33
     .             + reploc(13)*dxr(3,3)*r31r33
     .             + reploc(13)*r(3,3)*dxr31r33
     .             + dxreploc(14)*r(3,3)*radd03
     .             + reploc(14)*dxr(3,3)*radd03
     .             + reploc(14)*r(3,3)*dxradd03
     .             + dxreploc(15)*(r(1,3)*radd09
     .                         + r(2,3)*radd10)
     .             + reploc(15)*(dxr(1,3)*radd09
     .                         + dxr(2,3)*radd10)
     .             + reploc(15)*(r(1,3)*dxradd09
     .                         + r(2,3)*dxradd10)

        dyrepij(4,7) = dyreploc(13)*r(3,3)*r31r33
     .             + reploc(13)*dyr(3,3)*r31r33
     .             + reploc(13)*r(3,3)*dyr31r33
     .             + dyreploc(14)*r(3,3)*radd03
     .             + reploc(14)*dyr(3,3)*radd03
     .             + reploc(14)*r(3,3)*dyradd03
     .             + dyreploc(15)*(r(1,3)*radd09
     .                         + r(2,3)*radd10)
     .             + reploc(15)*(dyr(1,3)*radd09
     .                         + dyr(2,3)*radd10)
     .             + reploc(15)*(r(1,3)*dyradd09
     .                         + r(2,3)*dyradd10)

        dzrepij(4,7) = dzreploc(13)*r(3,3)*r31r33
     .             + reploc(13)*dzr(3,3)*r31r33
     .             + reploc(13)*r(3,3)*dzr31r33
     .             + dzreploc(14)*r(3,3)*radd03
     .             + reploc(14)*dzr(3,3)*radd03
     .             + reploc(14)*r(3,3)*dzradd03
     .             + dzreploc(15)*(r(1,3)*radd09
     .                         + r(2,3)*radd10)
     .             + reploc(15)*(dzr(1,3)*radd09
     .                         + dzr(2,3)*radd10)
     .             + reploc(15)*(r(1,3)*dzradd09
     .                         + r(2,3)*dzradd10)
!
!       (px pz | s pz):
!
        repij(7,4) = reploc(8)*r(3,3)*r31r33
     .             + reploc(9)*r(3,3)*radd03
     .             + reploc(10)*(r(1,3)*radd09
     .                         + r(2,3)*radd10)

        dxrepij(7,4) = dxreploc(8)*r(3,3)*r31r33
     .             + reploc(8)*dxr(3,3)*r31r33
     .             + reploc(8)*r(3,3)*dxr31r33
     .             + dxreploc(9)*r(3,3)*radd03
     .             + reploc(9)*dxr(3,3)*radd03
     .             + reploc(9)*r(3,3)*dxradd03
     .             + dxreploc(10)*(r(1,3)*radd09
     .                         + r(2,3)*radd10)
     .             + reploc(10)*(dxr(1,3)*radd09
     .                         + dxr(2,3)*radd10)
     .             + reploc(10)*(r(1,3)*dxradd09
     .                         + r(2,3)*dxradd10)

        dyrepij(7,4) = dyreploc(8)*r(3,3)*r31r33
     .             + reploc(8)*dyr(3,3)*r31r33
     .             + reploc(8)*r(3,3)*dyr31r33
     .             + dyreploc(9)*r(3,3)*radd03
     .             + reploc(9)*dyr(3,3)*radd03
     .             + reploc(9)*r(3,3)*dyradd03
     .             + dyreploc(10)*(r(1,3)*radd09
     .                         + r(2,3)*radd10)
     .             + reploc(10)*(dyr(1,3)*radd09
     .                         + dyr(2,3)*radd10)
     .             + reploc(10)*(r(1,3)*dyradd09
     .                         + r(2,3)*dyradd10)

        dzrepij(7,4) = dzreploc(8)*r(3,3)*r31r33
     .             + reploc(8)*dzr(3,3)*r31r33
     .             + reploc(8)*r(3,3)*dzr31r33
     .             + dzreploc(9)*r(3,3)*radd03
     .             + reploc(9)*dzr(3,3)*radd03
     .             + reploc(9)*r(3,3)*dzradd03
     .             + dzreploc(10)*(r(1,3)*radd09
     .                         + r(2,3)*radd10)
     .             + reploc(10)*(dzr(1,3)*radd09
     .                         + dzr(2,3)*radd10)
     .             + reploc(10)*(r(1,3)*dzradd09
     .                         + r(2,3)*dzradd10)
!
!       (s pz | py py):
!
        repij(4,8) = reploc(13)*r(3,3)*r32r32
     .             + reploc(14)*r(3,3)*radd04
     .             + reploc(15)*(r(1,3)*r12r32
     .                         + r(2,3)*r22r32)

        dxrepij(4,8) = dxreploc(13)*r(3,3)*r32r32
     .             + reploc(13)*dxr(3,3)*r32r32
     .             + reploc(13)*r(3,3)*dxr32r32
     .             + dxreploc(14)*r(3,3)*radd04
     .             + reploc(14)*dxr(3,3)*radd04
     .             + reploc(14)*r(3,3)*dxradd04
     .             + dxreploc(15)*(r(1,3)*r12r32
     .                         + r(2,3)*r22r32)
     .             + reploc(15)*(dxr(1,3)*r12r32
     .                         + dxr(2,3)*r22r32)
     .             + reploc(15)*(r(1,3)*dxr12r32
     .                         + r(2,3)*dxr22r32)

        dyrepij(4,8) = dyreploc(13)*r(3,3)*r32r32
     .             + reploc(13)*dyr(3,3)*r32r32
     .             + reploc(13)*r(3,3)*dyr32r32
     .             + dyreploc(14)*r(3,3)*radd04
     .             + reploc(14)*dyr(3,3)*radd04
     .             + reploc(14)*r(3,3)*dyradd04
     .             + dyreploc(15)*(r(1,3)*r12r32
     .                         + r(2,3)*r22r32)
     .             + reploc(15)*(dyr(1,3)*r12r32
     .                         + dyr(2,3)*r22r32)
     .             + reploc(15)*(r(1,3)*dyr12r32
     .                         + r(2,3)*dyr22r32)

        dzrepij(4,8) = dzreploc(13)*r(3,3)*r32r32
     .             + reploc(13)*dzr(3,3)*r32r32
     .             + reploc(13)*r(3,3)*dzr32r32
     .             + dzreploc(14)*r(3,3)*radd04
     .             + reploc(14)*dzr(3,3)*radd04
     .             + reploc(14)*r(3,3)*dzradd04
     .             + dzreploc(15)*(r(1,3)*r12r32
     .                         + r(2,3)*r22r32)
     .             + reploc(15)*(dzr(1,3)*r12r32
     .                         + dzr(2,3)*r22r32)
     .             + reploc(15)*(r(1,3)*dzr12r32
     .                         + r(2,3)*dzr22r32)
!
!       (py py | s pz):
!
        repij(8,4) = reploc(8)*r(3,3)*r32r32
     .             + reploc(9)*r(3,3)*radd04
     .             + reploc(10)*(r(1,3)*r12r32
     .                         + r(2,3)*r22r32)

        dxrepij(8,4) = dxreploc(8)*r(3,3)*r32r32
     .             + reploc(8)*dxr(3,3)*r32r32
     .             + reploc(8)*r(3,3)*dxr32r32
     .             + dxreploc(9)*r(3,3)*radd04
     .             + reploc(9)*dxr(3,3)*radd04
     .             + reploc(9)*r(3,3)*dxradd04
     .             + dxreploc(10)*(r(1,3)*r12r32
     .                         + r(2,3)*r22r32)
     .             + reploc(10)*(dxr(1,3)*r12r32
     .                         + dxr(2,3)*r22r32)
     .             + reploc(10)*(r(1,3)*dxr12r32
     .                         + r(2,3)*dxr22r32)

        dyrepij(8,4) = dyreploc(8)*r(3,3)*r32r32
     .             + reploc(8)*dyr(3,3)*r32r32
     .             + reploc(8)*r(3,3)*dyr32r32
     .             + dyreploc(9)*r(3,3)*radd04
     .             + reploc(9)*dyr(3,3)*radd04
     .             + reploc(9)*r(3,3)*dyradd04
     .             + dyreploc(10)*(r(1,3)*r12r32
     .                         + r(2,3)*r22r32)
     .             + reploc(10)*(dyr(1,3)*r12r32
     .                         + dyr(2,3)*r22r32)
     .             + reploc(10)*(r(1,3)*dyr12r32
     .                         + r(2,3)*dyr22r32)

        dzrepij(8,4) = dzreploc(8)*r(3,3)*r32r32
     .             + reploc(8)*dzr(3,3)*r32r32
     .             + reploc(8)*r(3,3)*dzr32r32
     .             + dzreploc(9)*r(3,3)*radd04
     .             + reploc(9)*dzr(3,3)*radd04
     .             + reploc(9)*r(3,3)*dzradd04
     .             + dzreploc(10)*(r(1,3)*r12r32
     .                         + r(2,3)*r22r32)
     .             + reploc(10)*(dzr(1,3)*r12r32
     .                         + dzr(2,3)*r22r32)
     .             + reploc(10)*(r(1,3)*dzr12r32
     .                         + r(2,3)*dzr22r32)

!
!       (s pz | py pz):
!
        repij(4,9) = reploc(13)*r(3,3)*r32r33
     .             + reploc(14)*r(3,3)*radd05
     .             + reploc(15)*(r(1,3)*radd13
     .                         + r(2,3)*radd14)

        dxrepij(4,9) = dxreploc(13)*r(3,3)*r32r33
     .             + reploc(13)*dxr(3,3)*r32r33
     .             + reploc(13)*r(3,3)*dxr32r33
     .             + dxreploc(14)*r(3,3)*radd05
     .             + reploc(14)*dxr(3,3)*radd05
     .             + reploc(14)*r(3,3)*dxradd05
     .             + dxreploc(15)*(r(1,3)*radd13
     .                         + r(2,3)*radd14)
     .             + reploc(15)*(dxr(1,3)*radd13
     .                         + dxr(2,3)*radd14)
     .             + reploc(15)*(r(1,3)*dxradd13
     .                         + r(2,3)*dxradd14)

        dyrepij(4,9) = dyreploc(13)*r(3,3)*r32r33
     .             + reploc(13)*dyr(3,3)*r32r33
     .             + reploc(13)*r(3,3)*dyr32r33
     .             + dyreploc(14)*r(3,3)*radd05
     .             + reploc(14)*dyr(3,3)*radd05
     .             + reploc(14)*r(3,3)*dyradd05
     .             + dyreploc(15)*(r(1,3)*radd13
     .                         + r(2,3)*radd14)
     .             + reploc(15)*(dyr(1,3)*radd13
     .                         + dyr(2,3)*radd14)
     .             + reploc(15)*(r(1,3)*dyradd13
     .                         + r(2,3)*dyradd14)

        dzrepij(4,9) = dzreploc(13)*r(3,3)*r32r33
     .             + reploc(13)*dzr(3,3)*r32r33
     .             + reploc(13)*r(3,3)*dzr32r33
     .             + dzreploc(14)*r(3,3)*radd05
     .             + reploc(14)*dzr(3,3)*radd05
     .             + reploc(14)*r(3,3)*dzradd05
     .             + dzreploc(15)*(r(1,3)*radd13
     .                         + r(2,3)*radd14)
     .             + reploc(15)*(dzr(1,3)*radd13
     .                         + dzr(2,3)*radd14)
     .             + reploc(15)*(r(1,3)*dzradd13
     .                         + r(2,3)*dzradd14)
!
!       (py pz | s pz):
!
        repij(9,4) = reploc(8)*r(3,3)*r32r33
     .             + reploc(9)*r(3,3)*radd05
     .             + reploc(10)*(r(1,3)*radd13
     .                         + r(2,3)*radd14)

        dxrepij(9,4) = dxreploc(8)*r(3,3)*r32r33
     .             + reploc(8)*dxr(3,3)*r32r33
     .             + reploc(8)*r(3,3)*dxr32r33
     .             + dxreploc(9)*r(3,3)*radd05
     .             + reploc(9)*dxr(3,3)*radd05
     .             + reploc(9)*r(3,3)*dxradd05
     .             + dxreploc(10)*(r(1,3)*radd13
     .                         + r(2,3)*radd14)
     .             + reploc(10)*(dxr(1,3)*radd13
     .                         + dxr(2,3)*radd14)
     .             + reploc(10)*(r(1,3)*dxradd13
     .                         + r(2,3)*dxradd14)

        dyrepij(9,4) = dyreploc(8)*r(3,3)*r32r33
     .             + reploc(8)*dyr(3,3)*r32r33
     .             + reploc(8)*r(3,3)*dyr32r33
     .             + dyreploc(9)*r(3,3)*radd05
     .             + reploc(9)*dyr(3,3)*radd05
     .             + reploc(9)*r(3,3)*dyradd05
     .             + dyreploc(10)*(r(1,3)*radd13
     .                         + r(2,3)*radd14)
     .             + reploc(10)*(dyr(1,3)*radd13
     .                         + dyr(2,3)*radd14)
     .             + reploc(10)*(r(1,3)*dyradd13
     .                         + r(2,3)*dyradd14)

        dzrepij(9,4) = dzreploc(8)*r(3,3)*r32r33
     .             + reploc(8)*dzr(3,3)*r32r33
     .             + reploc(8)*r(3,3)*dzr32r33
     .             + dzreploc(9)*r(3,3)*radd05
     .             + reploc(9)*dzr(3,3)*radd05
     .             + reploc(9)*r(3,3)*dzradd05
     .             + dzreploc(10)*(r(1,3)*radd13
     .                         + r(2,3)*radd14)
     .             + reploc(10)*(dzr(1,3)*radd13
     .                         + dzr(2,3)*radd14)
     .             + reploc(10)*(r(1,3)*dzradd13
     .                         + r(2,3)*dzradd14)
!
!       (s pz | pz pz):
!
        repij(4,10) = reploc(13)*r(3,3)*r33r33
     .              + reploc(14)*r(3,3)*radd06
     .              + reploc(15)*(r(1,3)*r13r33
     .                          + r(2,3)*r23r33)

        dxrepij(4,10) = dxreploc(13)*r(3,3)*r33r33
     .              + reploc(13)*dxr(3,3)*r33r33
     .              + reploc(13)*r(3,3)*dxr33r33
     .              + dxreploc(14)*r(3,3)*radd06
     .              + reploc(14)*dxr(3,3)*radd06
     .              + reploc(14)*r(3,3)*dxradd06
     .              + dxreploc(15)*(r(1,3)*r13r33
     .                          + r(2,3)*r23r33)
     .              + reploc(15)*(dxr(1,3)*r13r33
     .                          + dxr(2,3)*r23r33)
     .              + reploc(15)*(r(1,3)*dxr13r33
     .                          + r(2,3)*dxr23r33)

        dyrepij(4,10) = dyreploc(13)*r(3,3)*r33r33
     .              + reploc(13)*dyr(3,3)*r33r33
     .              + reploc(13)*r(3,3)*dyr33r33
     .              + dyreploc(14)*r(3,3)*radd06
     .              + reploc(14)*dyr(3,3)*radd06
     .              + reploc(14)*r(3,3)*dyradd06
     .              + dyreploc(15)*(r(1,3)*r13r33
     .                          + r(2,3)*r23r33)
     .              + reploc(15)*(dyr(1,3)*r13r33
     .                          + dyr(2,3)*r23r33)
     .              + reploc(15)*(r(1,3)*dyr13r33
     .                          + r(2,3)*dyr23r33)

        dzrepij(4,10) = dzreploc(13)*r(3,3)*r33r33
     .              + reploc(13)*dzr(3,3)*r33r33
     .              + reploc(13)*r(3,3)*dzr33r33
     .              + dzreploc(14)*r(3,3)*radd06
     .              + reploc(14)*dzr(3,3)*radd06
     .              + reploc(14)*r(3,3)*dzradd06
     .              + dzreploc(15)*(r(1,3)*r13r33
     .                          + r(2,3)*r23r33)
     .              + reploc(15)*(dzr(1,3)*r13r33
     .                          + dzr(2,3)*r23r33)
     .              + reploc(15)*(r(1,3)*dzr13r33
     .                          + r(2,3)*dzr23r33)

!
!       (pz pz | s pz):
!
        repij(10,4) = reploc(8)*r(3,3)*r33r33
     .              + reploc(9)*r(3,3)*radd06
     .              + reploc(10)*(r(1,3)*r13r33
     .                          + r(2,3)*r23r33)

        dxrepij(10,4) = dxreploc(8)*r(3,3)*r33r33
     .              + reploc(8)*dxr(3,3)*r33r33
     .              + reploc(8)*r(3,3)*dxr33r33
     .              + dxreploc(9)*r(3,3)*radd06
     .              + reploc(9)*dxr(3,3)*radd06
     .              + reploc(9)*r(3,3)*dxradd06
     .              + dxreploc(10)*(r(1,3)*r13r33
     .                          + r(2,3)*r23r33)
     .              + reploc(10)*(dxr(1,3)*r13r33
     .                          + dxr(2,3)*r23r33)
     .              + reploc(10)*(r(1,3)*dxr13r33
     .                          + r(2,3)*dxr23r33)

        dyrepij(10,4) = dyreploc(8)*r(3,3)*r33r33
     .              + reploc(8)*dyr(3,3)*r33r33
     .              + reploc(8)*r(3,3)*dyr33r33
     .              + dyreploc(9)*r(3,3)*radd06
     .              + reploc(9)*dyr(3,3)*radd06
     .              + reploc(9)*r(3,3)*dyradd06
     .              + dyreploc(10)*(r(1,3)*r13r33
     .                          + r(2,3)*r23r33)
     .              + reploc(10)*(dyr(1,3)*r13r33
     .                          + dyr(2,3)*r23r33)
     .              + reploc(10)*(r(1,3)*dyr13r33
     .                          + r(2,3)*dyr23r33)

        dzrepij(10,4) = dzreploc(8)*r(3,3)*r33r33
     .              + reploc(8)*dzr(3,3)*r33r33
     .              + reploc(8)*r(3,3)*dzr33r33
     .              + dzreploc(9)*r(3,3)*radd06
     .              + reploc(9)*dzr(3,3)*radd06
     .              + reploc(9)*r(3,3)*dzradd06
     .              + dzreploc(10)*(r(1,3)*r13r33
     .                          + r(2,3)*r23r33)
     .              + reploc(10)*(dzr(1,3)*r13r33
     .                          + dzr(2,3)*r23r33)
     .              + reploc(10)*(r(1,3)*dzr13r33
     .                          + r(2,3)*dzr23r33)
!
!       (px px | px px):
!
        repij(5,5) = reploc(16)*r31r31*r31r31
     .             + reploc(17)*radd01*r31r31
     .             + reploc(18)*radd01*r31r31
     .             + reploc(19)*(r11r11*r11r11
     .                         + r21r21*r21r21)
     .             + reploc(20)*(r11r31*r11r31
     .                         + r21r31*r21r31)
     .             + reploc(21)*r11r11*r21r21*2.0d0
     .             + reploc(22)*r11r21*r11r21

        dxrepij(5,5) = dxreploc(16)*r31r31*r31r31
     .             + reploc(16)*dxr31r31*r31r31
     .             + reploc(16)*r31r31*dxr31r31
     .             + dxreploc(17)*radd01*r31r31
     .             + reploc(17)*dxradd01*r31r31
     .             + reploc(17)*radd01*dxr31r31
     .             + dxreploc(18)*radd01*r31r31
     .             + reploc(18)*dxradd01*r31r31
     .             + reploc(18)*radd01*dxr31r31
     .             + dxreploc(19)*(r11r11*r11r11
     .                         + r21r21*r21r21)
     .             + reploc(19)*(dxr11r11*r11r11
     .                         + dxr21r21*r21r21)
     .             + reploc(19)*(r11r11*dxr11r11
     .                         + r21r21*dxr21r21)
        dxrepij(5,5) = dxrepij(5,5)
     .             + dxreploc(20)*(r11r31*r11r31
     .                         + r21r31*r21r31)
     .             + reploc(20)*(dxr11r31*r11r31
     .                         + dxr21r31*r21r31)
     .             + reploc(20)*(r11r31*dxr11r31
     .                         + r21r31*dxr21r31)
     .             + dxreploc(21)*r11r11*r21r21*2.0d0
     .             + reploc(21)*dxr11r11*r21r21*2.0d0
     .             + reploc(21)*r11r11*dxr21r21*2.0d0
     .             + dxreploc(22)*r11r21*r11r21
     .             + reploc(22)*dxr11r21*r11r21
     .             + reploc(22)*r11r21*dxr11r21

        dyrepij(5,5) = dyreploc(16)*r31r31*r31r31
     .             + reploc(16)*dyr31r31*r31r31
     .             + reploc(16)*r31r31*dyr31r31
     .             + dyreploc(17)*radd01*r31r31
     .             + reploc(17)*dyradd01*r31r31
     .             + reploc(17)*radd01*dyr31r31
     .             + dyreploc(18)*radd01*r31r31
     .             + reploc(18)*dyradd01*r31r31
     .             + reploc(18)*radd01*dyr31r31
     .             + dyreploc(19)*(r11r11*r11r11
     .                         + r21r21*r21r21)
     .             + reploc(19)*(dyr11r11*r11r11
     .                         + dyr21r21*r21r21)
     .             + reploc(19)*(r11r11*dyr11r11
     .                         + r21r21*dyr21r21)
         dyrepij(5,5) = dyrepij(5,5)
     .             + dyreploc(20)*(r11r31*r11r31
     .                         + r21r31*r21r31)
     .             + reploc(20)*(dyr11r31*r11r31
     .                         + dyr21r31*r21r31)
     .             + reploc(20)*(r11r31*dyr11r31
     .                         + r21r31*dyr21r31)
     .             + dyreploc(21)*r11r11*r21r21*2.0d0
     .             + reploc(21)*dyr11r11*r21r21*2.0d0
     .             + reploc(21)*r11r11*dyr21r21*2.0d0
     .             + dyreploc(22)*r11r21*r11r21
     .             + reploc(22)*dyr11r21*r11r21
     .             + reploc(22)*r11r21*dyr11r21

        dzrepij(5,5) = dzreploc(16)*r31r31*r31r31
     .             + reploc(16)*dzr31r31*r31r31
     .             + reploc(16)*r31r31*dzr31r31
     .             + dzreploc(17)*radd01*r31r31
     .             + reploc(17)*dzradd01*r31r31
     .             + reploc(17)*radd01*dzr31r31
     .             + dzreploc(18)*radd01*r31r31
     .             + reploc(18)*dzradd01*r31r31
     .             + reploc(18)*radd01*dzr31r31
     .             + dzreploc(19)*(r11r11*r11r11
     .                         + r21r21*r21r21)
     .             + reploc(19)*(dzr11r11*r11r11
     .                         + dzr21r21*r21r21)
     .             + reploc(19)*(r11r11*dzr11r11
     .                         + r21r21*dzr21r21)
        dzrepij(5,5) = dzrepij(5,5)
     .             + dzreploc(20)*(r11r31*r11r31
     .                         + r21r31*r21r31)
     .             + reploc(20)*(dzr11r31*r11r31
     .                         + dzr21r31*r21r31)
     .             + reploc(20)*(r11r31*dzr11r31
     .                         + r21r31*dzr21r31)
     .             + dzreploc(21)*r11r11*r21r21*2.0d0
     .             + reploc(21)*dzr11r11*r21r21*2.0d0
     .             + reploc(21)*r11r11*dzr21r21*2.0d0
     .             + dzreploc(22)*r11r21*r11r21
     .             + reploc(22)*dzr11r21*r11r21
     .             + reploc(22)*r11r21*dzr11r21
!
!       (px px | px py):
!
        repij(5,6) = reploc(16)*r31r31*r31r32
     .             + reploc(17)*radd01*r31r32
     .             + reploc(18)*radd02*r31r31
     .             + reploc(19)*(r11r11*r11r12
     .                         + r21r21*r21r22)
     .             + reploc(20)*(r11r31*r11r32
     .                         + r31r11*r31r12
     .                         + r21r31*r21r32
     .                         + r31r21*r31r22)
     .             + reploc(21)*(r11r11*r21r22
     .                         + r21r21*r11r12)
     .             + reploc(22)*(r11r21*r11r22
     .                         + r21r11*r21r12)

        dxrepij(5,6) = dxreploc(16)*r31r31*r31r32
     .             + reploc(16)*dxr31r31*r31r32
     .             + reploc(16)*r31r31*dxr31r32
     .             + dxreploc(17)*radd01*r31r32
     .             + reploc(17)*dxradd01*r31r32
     .             + reploc(17)*radd01*dxr31r32
     .             + dxreploc(18)*radd02*r31r31
     .             + reploc(18)*dxradd02*r31r31
     .             + reploc(18)*radd02*dxr31r31
     .             + dxreploc(19)*(r11r11*r11r12
     .                         + r21r21*r21r22)
     .             + reploc(19)*(dxr11r11*r11r12
     .                         + dxr21r21*r21r22)
     .             + reploc(19)*(r11r11*dxr11r12
     .                         + r21r21*dxr21r22)
        dxrepij(5,6) = dxrepij(5,6)
     .             + dxreploc(20)*(r11r31*r11r32
     .                         + r31r11*r31r12
     .                         + r21r31*r21r32
     .                         + r31r21*r31r22)
     .             + reploc(20)*(dxr11r31*r11r32
     .                         + dxr31r11*r31r12
     .                         + dxr21r31*r21r32
     .                         + dxr31r21*r31r22)
     .             + reploc(20)*(r11r31*dxr11r32
     .                         + r31r11*dxr31r12
     .                         + r21r31*dxr21r32
     .                         + r31r21*dxr31r22)
     .             + dxreploc(21)*(r11r11*r21r22
     .                         + r21r21*r11r12)
     .             + reploc(21)*(dxr11r11*r21r22
     .                         + dxr21r21*r11r12)
     .             + reploc(21)*(r11r11*dxr21r22
     .                         + r21r21*dxr11r12)
     .             + dxreploc(22)*(r11r21*r11r22
     .                         + r21r11*r21r12)
     .             + reploc(22)*(dxr11r21*r11r22
     .                         + dxr21r11*r21r12)
     .             + reploc(22)*(r11r21*dxr11r22
     .                         + r21r11*dxr21r12)

        dyrepij(5,6) = dyreploc(16)*r31r31*r31r32
     .             + reploc(16)*dyr31r31*r31r32
     .             + reploc(16)*r31r31*dyr31r32
     .             + dyreploc(17)*radd01*r31r32
     .             + reploc(17)*dyradd01*r31r32
     .             + reploc(17)*radd01*dyr31r32
     .             + dyreploc(18)*radd02*r31r31
     .             + reploc(18)*dyradd02*r31r31
     .             + reploc(18)*radd02*dyr31r31
     .             + dyreploc(19)*(r11r11*r11r12
     .                         + r21r21*r21r22)
     .             + reploc(19)*(dyr11r11*r11r12
     .                         + dyr21r21*r21r22)
     .             + reploc(19)*(r11r11*dyr11r12
     .                         + r21r21*dyr21r22)
        dyrepij(5,6) = dyrepij(5,6)
     .             + dyreploc(20)*(r11r31*r11r32
     .                         + r31r11*r31r12
     .                         + r21r31*r21r32
     .                         + r31r21*r31r22)
     .             + reploc(20)*(dyr11r31*r11r32
     .                         + dyr31r11*r31r12
     .                         + dyr21r31*r21r32
     .                         + dyr31r21*r31r22)
     .             + reploc(20)*(r11r31*dyr11r32
     .                         + r31r11*dyr31r12
     .                         + r21r31*dyr21r32
     .                         + r31r21*dyr31r22)
     .             + dyreploc(21)*(r11r11*r21r22
     .                         + r21r21*r11r12)
     .             + reploc(21)*(dyr11r11*r21r22
     .                         + dyr21r21*r11r12)
     .             + reploc(21)*(r11r11*dyr21r22
     .                         + r21r21*dyr11r12)
     .             + dyreploc(22)*(r11r21*r11r22
     .                         + r21r11*r21r12)
     .             + reploc(22)*(dyr11r21*r11r22
     .                         + dyr21r11*r21r12)
     .             + reploc(22)*(r11r21*dyr11r22
     .                         + r21r11*dyr21r12)

        dzrepij(5,6) = dzreploc(16)*r31r31*r31r32
     .             + reploc(16)*dzr31r31*r31r32
     .             + reploc(16)*r31r31*dzr31r32
     .             + dzreploc(17)*radd01*r31r32
     .             + reploc(17)*dzradd01*r31r32
     .             + reploc(17)*radd01*dzr31r32
     .             + dzreploc(18)*radd02*r31r31
     .             + reploc(18)*dzradd02*r31r31
     .             + reploc(18)*radd02*dzr31r31
     .             + dzreploc(19)*(r11r11*r11r12
     .                         + r21r21*r21r22)
     .             + reploc(19)*(dzr11r11*r11r12
     .                         + dzr21r21*r21r22)
     .             + reploc(19)*(r11r11*dzr11r12
     .                         + r21r21*dzr21r22)
     .             + dzreploc(20)*(r11r31*r11r32
     .                         + r31r11*r31r12
     .                         + r21r31*r21r32
     .                         + r31r21*r31r22)
         dzrepij(5,6) =dzrepij(5,6) 
     .             + reploc(20)*(dzr11r31*r11r32
     .                         + dzr31r11*r31r12
     .                         + dzr21r31*r21r32
     .                         + dzr31r21*r31r22)
     .             + reploc(20)*(r11r31*dzr11r32
     .                         + r31r11*dzr31r12
     .                         + r21r31*dzr21r32
     .                         + r31r21*dzr31r22)
     .             + dzreploc(21)*(r11r11*r21r22
     .                         + r21r21*r11r12)
     .             + reploc(21)*(dzr11r11*r21r22
     .                         + dzr21r21*r11r12)
     .             + reploc(21)*(r11r11*dzr21r22
     .                         + r21r21*dzr11r12)
     .             + dzreploc(22)*(r11r21*r11r22
     .                         + r21r11*r21r12)
     .             + reploc(22)*(dzr11r21*r11r22
     .                         + dzr21r11*r21r12)
     .             + reploc(22)*(r11r21*dzr11r22
     .                         + r21r11*dzr21r12)

     .     
!
!       (px py | px px):
!
        repij(6,5) = reploc(16)*r31r32*r31r31
     .             + reploc(17)*radd02*r31r31
     .             + reploc(18)*radd01*r31r32
     .             + reploc(19)*(r11r12*r11r11
     .                         + r21r22*r21r21)
     .             + reploc(20)*(radd11*r11r31
     .                         + radd12*r21r31)
     .             + reploc(21)*(r11r12*r21r21
     .                         + r21r22*r11r11)
     .             + reploc(22)*radd17*r11r21

        dxrepij(6,5) = dxreploc(16)*r31r32*r31r31
     .             + reploc(16)*dxr31r32*r31r31
     .             + reploc(16)*r31r32*dxr31r31
     .             + dxreploc(17)*radd02*r31r31
     .             + reploc(17)*dxradd02*r31r31
     .             + reploc(17)*radd02*dxr31r31
     .             + dxreploc(18)*radd01*r31r32
     .             + reploc(18)*dxradd01*r31r32
     .             + reploc(18)*radd01*dxr31r32
     .             + dxreploc(19)*(r11r12*r11r11
     .                         + r21r22*r21r21)
     .             + reploc(19)*(dxr11r12*r11r11
     .                         + dxr21r22*r21r21)
     .             + reploc(19)*(r11r12*dxr11r11
     .                         + r21r22*dxr21r21)
     .             + dxreploc(20)*(radd11*r11r31
     .                         + radd12*r21r31)
     .             + reploc(20)*(dxradd11*r11r31
     .                         + dxradd12*r21r31)
          dxrepij(6,5)=dxrepij(6,5)
     .             + reploc(20)*(radd11*dxr11r31
     .                         + radd12*dxr21r31)
     .             + dxreploc(21)*(r11r12*r21r21
     .                         + r21r22*r11r11)
     .             + reploc(21)*(dxr11r12*r21r21
     .                         + dxr21r22*r11r11)
     .             + reploc(21)*(r11r12*dxr21r21
     .                         + r21r22*dxr11r11)
     .             + dxreploc(22)*radd17*r11r21
     .             + reploc(22)*dxradd17*r11r21
     .             + reploc(22)*radd17*dxr11r21

        dyrepij(6,5) = dyreploc(16)*r31r32*r31r31
     .             + reploc(16)*dyr31r32*r31r31
     .             + reploc(16)*r31r32*dyr31r31
     .             + dyreploc(17)*radd02*r31r31
     .             + reploc(17)*dyradd02*r31r31
     .             + reploc(17)*radd02*dyr31r31
     .             + dyreploc(18)*radd01*r31r32
     .             + reploc(18)*dyradd01*r31r32
     .             + reploc(18)*radd01*dyr31r32
     .             + dyreploc(19)*(r11r12*r11r11
     .                         + r21r22*r21r21)
     .             + reploc(19)*(dyr11r12*r11r11
     .                         + dyr21r22*r21r21)
     .             + reploc(19)*(r11r12*dyr11r11
     .                         + r21r22*dyr21r21)
     .             + dyreploc(20)*(radd11*r11r31
     .                         + radd12*r21r31)
     .             + reploc(20)*(dyradd11*r11r31
     .                         + dyradd12*r21r31)
         dyrepij(6,5)=dyrepij(6,5)
     .             + reploc(20)*(radd11*dyr11r31
     .                         + radd12*dyr21r31)
     .             + dyreploc(21)*(r11r12*r21r21
     .                         + r21r22*r11r11)
     .             + reploc(21)*(dyr11r12*r21r21
     .                         + dyr21r22*r11r11)
     .             + reploc(21)*(r11r12*dyr21r21
     .                         + r21r22*dyr11r11)
     .             + dyreploc(22)*radd17*r11r21
     .             + reploc(22)*dyradd17*r11r21
     .             + reploc(22)*radd17*dyr11r21

        dzrepij(6,5) = dzreploc(16)*r31r32*r31r31
     .             + reploc(16)*dzr31r32*r31r31
     .             + reploc(16)*r31r32*dzr31r31
     .             + dzreploc(17)*radd02*r31r31
     .             + reploc(17)*dzradd02*r31r31
     .             + reploc(17)*radd02*dzr31r31
     .             + dzreploc(18)*radd01*r31r32
     .             + reploc(18)*dzradd01*r31r32
     .             + reploc(18)*radd01*dzr31r32
     .             + dzreploc(19)*(r11r12*r11r11
     .                         + r21r22*r21r21)
     .             + reploc(19)*(dzr11r12*r11r11
     .                         + dzr21r22*r21r21)
     .             + reploc(19)*(r11r12*dzr11r11
     .                         + r21r22*dzr21r21)
     .             + dzreploc(20)*(radd11*r11r31
     .                         + radd12*r21r31)
     .             + reploc(20)*(dzradd11*r11r31
     .                         + dzradd12*r21r31)
         dzrepij(6,5)= dzrepij(6,5)
     .             + reploc(20)*(radd11*dzr11r31
     .                         + radd12*dzr21r31)
     .             + dzreploc(21)*(r11r12*r21r21
     .                         + r21r22*r11r11)
     .             + reploc(21)*(dzr11r12*r21r21
     .                         + dzr21r22*r11r11)
     .             + reploc(21)*(r11r12*dzr21r21
     .                         + r21r22*dzr11r11)
     .             + dzreploc(22)*radd17*r11r21
     .             + reploc(22)*dzradd17*r11r21
     .             + reploc(22)*radd17*dzr11r21
!
!       (px px | px pz):
!
        repij(5,7) = reploc(16)*r31r31*r31r33
     .             + reploc(17)*radd01*r31r33
     .             + reploc(18)*radd03*r31r31
     .             + reploc(19)*(r11r11*r11r13
     .                         + r21r21*r21r23)
     .             + reploc(20)*(r11r31*r11r33
     .                         + r31r11*r31r13
     .                         + r21r31*r21r33
     .                         + r31r21*r31r23)
     .             + reploc(21)*(r11r11*r21r23
     .                         + r21r21*r11r13)
     .             + reploc(22)*(r11r21*r11r23
     .                         + r21r11*r21r13)

        dxrepij(5,7) = dxreploc(16)*r31r31*r31r33
     .             + reploc(16)*dxr31r31*r31r33
     .             + reploc(16)*r31r31*dxr31r33
     .             + dxreploc(17)*radd01*r31r33
     .             + reploc(17)*dxradd01*r31r33
     .             + reploc(17)*radd01*dxr31r33
     .             + dxreploc(18)*radd03*r31r31
     .             + reploc(18)*dxradd03*r31r31
     .             + reploc(18)*radd03*dxr31r31
     .             + dxreploc(19)*(r11r11*r11r13
     .                         + r21r21*r21r23)
     .             + reploc(19)*(dxr11r11*r11r13
     .                         + dxr21r21*r21r23)
     .             + reploc(19)*(r11r11*dxr11r13
     .                         + r21r21*dxr21r23)
     .             + dxreploc(20)*(r11r31*r11r33
     .                         + r31r11*r31r13
     .                         + r21r31*r21r33
     .                         + r31r21*r31r23)
             dxrepij(5,7)=dxrepij(5,7)
     .             + reploc(20)*(dxr11r31*r11r33
     .                         + dxr31r11*r31r13
     .                         + dxr21r31*r21r33
     .                         + dxr31r21*r31r23)
     .             + reploc(20)*(r11r31*dxr11r33
     .                         + r31r11*dxr31r13
     .                         + r21r31*dxr21r33
     .                         + r31r21*dxr31r23)
     .             + dxreploc(21)*(r11r11*r21r23
     .                         + r21r21*r11r13)
     .             + reploc(21)*(dxr11r11*r21r23
     .                         + dxr21r21*r11r13)
     .             + reploc(21)*(r11r11*dxr21r23
     .                         + r21r21*dxr11r13)
     .             + dxreploc(22)*(r11r21*r11r23
     .                         + r21r11*r21r13)
     .             + reploc(22)*(dxr11r21*r11r23
     .                         + dxr21r11*r21r13)
     .             + reploc(22)*(r11r21*dxr11r23
     .                         + r21r11*dxr21r13)

        dyrepij(5,7) = dyreploc(16)*r31r31*r31r33
     .             + reploc(16)*dyr31r31*r31r33
     .             + reploc(16)*r31r31*dyr31r33
     .             + dyreploc(17)*radd01*r31r33
     .             + reploc(17)*dyradd01*r31r33
     .             + reploc(17)*radd01*dyr31r33
     .             + dyreploc(18)*radd03*r31r31
     .             + reploc(18)*dyradd03*r31r31
     .             + reploc(18)*radd03*dyr31r31
     .             + dyreploc(19)*(r11r11*r11r13
     .                         + r21r21*r21r23)
     .             + reploc(19)*(dyr11r11*r11r13
     .                         + dyr21r21*r21r23)
     .             + reploc(19)*(r11r11*dyr11r13
     .                         + r21r21*dyr21r23)
     .             + dyreploc(20)*(r11r31*r11r33
     .                         + r31r11*r31r13
     .                         + r21r31*r21r33
     .                         + r31r21*r31r23)
          dyrepij(5,7)=dyrepij(5,7)
     .             + reploc(20)*(dyr11r31*r11r33
     .                         + dyr31r11*r31r13
     .                         + dyr21r31*r21r33
     .                         + dyr31r21*r31r23)
     .             + reploc(20)*(r11r31*dyr11r33
     .                         + r31r11*dyr31r13
     .                         + r21r31*dyr21r33
     .                         + r31r21*dyr31r23)
     .             + dyreploc(21)*(r11r11*r21r23
     .                         + r21r21*r11r13)
     .             + reploc(21)*(dyr11r11*r21r23
     .                         + dyr21r21*r11r13)
     .             + reploc(21)*(r11r11*dyr21r23
     .                         + r21r21*dyr11r13)
     .             + dyreploc(22)*(r11r21*r11r23
     .                         + r21r11*r21r13)
     .             + reploc(22)*(dyr11r21*r11r23
     .                         + dyr21r11*r21r13)
     .             + reploc(22)*(r11r21*dyr11r23
     .                         + r21r11*dyr21r13)

        dzrepij(5,7) = dzreploc(16)*r31r31*r31r33
     .             + reploc(16)*dzr31r31*r31r33
     .             + reploc(16)*r31r31*dzr31r33
     .             + dzreploc(17)*radd01*r31r33
     .             + reploc(17)*dzradd01*r31r33
     .             + reploc(17)*radd01*dzr31r33
     .             + dzreploc(18)*radd03*r31r31
     .             + reploc(18)*dzradd03*r31r31
     .             + reploc(18)*radd03*dzr31r31
     .             + dzreploc(19)*(r11r11*r11r13
     .                         + r21r21*r21r23)
     .             + reploc(19)*(dzr11r11*r11r13
     .                         + dzr21r21*r21r23)
     .             + reploc(19)*(r11r11*dzr11r13
     .                         + r21r21*dzr21r23)
     .             + dzreploc(20)*(r11r31*r11r33
     .                         + r31r11*r31r13
     .                         + r21r31*r21r33
     .                         + r31r21*r31r23)
            dzrepij(5,7)= dzrepij(5,7)
     .             + reploc(20)*(dzr11r31*r11r33
     .                         + dzr31r11*r31r13
     .                         + dzr21r31*r21r33
     .                         + dzr31r21*r31r23)
     .             + reploc(20)*(r11r31*dzr11r33
     .                         + r31r11*dzr31r13
     .                         + r21r31*dzr21r33
     .                         + r31r21*dzr31r23)
     .             + dzreploc(21)*(r11r11*r21r23
     .                         + r21r21*r11r13)
     .             + reploc(21)*(dzr11r11*r21r23
     .                         + dzr21r21*r11r13)
     .             + reploc(21)*(r11r11*dzr21r23
     .                         + r21r21*dzr11r13)
     .             + dzreploc(22)*(r11r21*r11r23
     .                         + r21r11*r21r13)
     .             + reploc(22)*(dzr11r21*r11r23
     .                         + dzr21r11*r21r13)
     .             + reploc(22)*(r11r21*dzr11r23
     .                         + r21r11*dzr21r13)
!
!       (px pz | px px):
!
        repij(7,5) = reploc(16)*r31r33*r31r31
     .             + reploc(17)*radd03*r31r31
     .             + reploc(18)*radd01*r31r33
     .             + reploc(19)*(r11r13*r11r11
     .                         + r21r23*r21r21)
     .             + reploc(20)*(radd09*r11r31
     .                          +radd10*r21r31)
     .             + reploc(21)*(r11r13*r21r21
     .                         + r21r23*r11r11)
     .             + reploc(22)*radd15*r11r21

        dxrepij(7,5) = dxreploc(16)*r31r33*r31r31
     .             + reploc(16)*dxr31r33*r31r31
     .             + reploc(16)*r31r33*dxr31r31
     .             + dxreploc(17)*radd03*r31r31
     .             + reploc(17)*dxradd03*r31r31
     .             + reploc(17)*radd03*dxr31r31
     .             + dxreploc(18)*radd01*r31r33
     .             + reploc(18)*dxradd01*r31r33
     .             + reploc(18)*radd01*dxr31r33
     .             + dxreploc(19)*(r11r13*r11r11
     .                         + r21r23*r21r21)
     .             + reploc(19)*(dxr11r13*r11r11
     .                         + dxr21r23*r21r21)
     .             + reploc(19)*(r11r13*dxr11r11
     .                         + r21r23*dxr21r21)
            dxrepij(7,5)=dxrepij(7,5)
     .             + dxreploc(20)*(radd09*r11r31
     .                          +radd10*r21r31)
     .             + reploc(20)*(dxradd09*r11r31
     .                          +dxradd10*r21r31)
     .             + reploc(20)*(radd09*dxr11r31
     .                          +radd10*dxr21r31)
     .             + dxreploc(21)*(r11r13*r21r21
     .                         + r21r23*r11r11)
     .             + reploc(21)*(dxr11r13*r21r21
     .                         + dxr21r23*r11r11)
     .             + reploc(21)*(r11r13*dxr21r21
     .                         + r21r23*dxr11r11)
     .             + dxreploc(22)*radd15*r11r21
     .             + reploc(22)*dxradd15*r11r21
     .             + reploc(22)*radd15*dxr11r21

        dyrepij(7,5) = dyreploc(16)*r31r33*r31r31
     .             + reploc(16)*dyr31r33*r31r31
     .             + reploc(16)*r31r33*dyr31r31
     .             + dyreploc(17)*radd03*r31r31
     .             + reploc(17)*dyradd03*r31r31
     .             + reploc(17)*radd03*dyr31r31
     .             + dyreploc(18)*radd01*r31r33
     .             + reploc(18)*dyradd01*r31r33
     .             + reploc(18)*radd01*dyr31r33
     .             + dyreploc(19)*(r11r13*r11r11
     .                         + r21r23*r21r21)
     .             + reploc(19)*(dyr11r13*r11r11
     .                         + dyr21r23*r21r21)
     .             + reploc(19)*(r11r13*dyr11r11
     .                         + r21r23*dyr21r21)
           dyrepij(7,5)=dyrepij(7,5)
     .             + dyreploc(20)*(radd09*r11r31
     .                          +radd10*r21r31)
     .             + reploc(20)*(dyradd09*r11r31
     .                          +dyradd10*r21r31)
     .             + reploc(20)*(radd09*dyr11r31
     .                          +radd10*dyr21r31)
     .             + dyreploc(21)*(r11r13*r21r21
     .                         + r21r23*r11r11)
     .             + reploc(21)*(dyr11r13*r21r21
     .                         + dyr21r23*r11r11)
     .             + reploc(21)*(r11r13*dyr21r21
     .                         + r21r23*dyr11r11)
     .             + dyreploc(22)*radd15*r11r21
     .             + reploc(22)*dyradd15*r11r21
     .             + reploc(22)*radd15*dyr11r21

        dzrepij(7,5) = dzreploc(16)*r31r33*r31r31
     .             + reploc(16)*dzr31r33*r31r31
     .             + reploc(16)*r31r33*dzr31r31
     .             + dzreploc(17)*radd03*r31r31
     .             + reploc(17)*dzradd03*r31r31
     .             + reploc(17)*radd03*dzr31r31
     .             + dzreploc(18)*radd01*r31r33
     .             + reploc(18)*dzradd01*r31r33
     .             + reploc(18)*radd01*dzr31r33
     .             + dzreploc(19)*(r11r13*r11r11
     .                         + r21r23*r21r21)
     .             + reploc(19)*(dzr11r13*r11r11
     .                         + dzr21r23*r21r21)
     .             + reploc(19)*(r11r13*dzr11r11
     .                         + r21r23*dzr21r21)
            dzrepij(7,5)=dzrepij(7,5)
     .             + dzreploc(20)*(radd09*r11r31
     .                          +radd10*r21r31)
     .             + reploc(20)*(dzradd09*r11r31
     .                          +dzradd10*r21r31)
     .             + reploc(20)*(radd09*dzr11r31
     .                          +radd10*dzr21r31)
     .             + dzreploc(21)*(r11r13*r21r21
     .                         + r21r23*r11r11)
     .             + reploc(21)*(dzr11r13*r21r21
     .                         + dzr21r23*r11r11)
     .             + reploc(21)*(r11r13*dzr21r21
     .                         + r21r23*dzr11r11)
     .             + dzreploc(22)*radd15*r11r21
     .             + reploc(22)*dzradd15*r11r21
     .             + reploc(22)*radd15*dzr11r21
 
!
!       (px px | py py):
!
        repij(5,8) = reploc(16)*r31r31*r32r32
     .             + reploc(17)*radd01*r32r32
     .             + reploc(18)*radd04*r31r31
     .             + reploc(19)*(r11r11*r12r12
     .                         + r21r21*r22r22)
     .             + reploc(20)*(r11r31*r12r32
     .                         + r21r31*r22r32)
     .             + reploc(21)*(r11r11*r22r22
     .                         + r21r21*r12r12)
     .             + reploc(22)*r11r21*r12r22

        dxrepij(5,8) = dxreploc(16)*r31r31*r32r32
     .             + reploc(16)*dxr31r31*r32r32
     .             + reploc(16)*r31r31*dxr32r32
     .             + dxreploc(17)*radd01*r32r32
     .             + reploc(17)*dxradd01*r32r32
     .             + reploc(17)*radd01*dxr32r32
     .             + dxreploc(18)*radd04*r31r31
     .             + reploc(18)*dxradd04*r31r31
     .             + reploc(18)*radd04*dxr31r31
     .             + dxreploc(19)*(r11r11*r12r12
     .                         + r21r21*r22r22)
     .             + reploc(19)*(dxr11r11*r12r12
     .                         + dxr21r21*r22r22)
     .             + reploc(19)*(r11r11*dxr12r12
     .                         + r21r21*dxr22r22)
               dxrepij(5,8)=dxrepij(5,8)
     .             + dxreploc(20)*(r11r31*r12r32
     .                         + r21r31*r22r32)
     .             + reploc(20)*(dxr11r31*r12r32
     .                         + dxr21r31*r22r32)
     .             + reploc(20)*(r11r31*dxr12r32
     .                         + r21r31*dxr22r32)
     .             + dxreploc(21)*(r11r11*r22r22
     .                         + r21r21*r12r12)
     .             + reploc(21)*(dxr11r11*r22r22
     .                         + dxr21r21*r12r12)
     .             + reploc(21)*(r11r11*dxr22r22
     .                         + r21r21*dxr12r12)
     .             + dxreploc(22)*r11r21*r12r22
     .             + reploc(22)*dxr11r21*r12r22
     .             + reploc(22)*r11r21*dxr12r22

        dyrepij(5,8) = dyreploc(16)*r31r31*r32r32
     .             + reploc(16)*dyr31r31*r32r32
     .             + reploc(16)*r31r31*dyr32r32
     .             + dyreploc(17)*radd01*r32r32
     .             + reploc(17)*dyradd01*r32r32
     .             + reploc(17)*radd01*dyr32r32
     .             + dyreploc(18)*radd04*r31r31
     .             + reploc(18)*dyradd04*r31r31
     .             + reploc(18)*radd04*dyr31r31
     .             + dyreploc(19)*(r11r11*r12r12
     .                         + r21r21*r22r22)
     .             + reploc(19)*(dyr11r11*r12r12
     .                         + dyr21r21*r22r22)
     .             + reploc(19)*(r11r11*dyr12r12
     .                         + r21r21*dyr22r22)
              dyrepij(5,8)=dyrepij(5,8)
     .             + dyreploc(20)*(r11r31*r12r32
     .                         + r21r31*r22r32)
     .             + reploc(20)*(dyr11r31*r12r32
     .                         + dyr21r31*r22r32)
     .             + reploc(20)*(r11r31*dyr12r32
     .                         + r21r31*dyr22r32)
     .             + dyreploc(21)*(r11r11*r22r22
     .                         + r21r21*r12r12)
     .             + reploc(21)*(dyr11r11*r22r22
     .                         + dyr21r21*r12r12)
     .             + reploc(21)*(r11r11*dyr22r22
     .                         + r21r21*dyr12r12)
     .             + dyreploc(22)*r11r21*r12r22
     .             + reploc(22)*dyr11r21*r12r22
     .             + reploc(22)*r11r21*dyr12r22

        dzrepij(5,8) = dzreploc(16)*r31r31*r32r32
     .             + reploc(16)*dzr31r31*r32r32
     .             + reploc(16)*r31r31*dzr32r32
     .             + dzreploc(17)*radd01*r32r32
     .             + reploc(17)*dzradd01*r32r32
     .             + reploc(17)*radd01*dzr32r32
     .             + dzreploc(18)*radd04*r31r31
     .             + reploc(18)*dzradd04*r31r31
     .             + reploc(18)*radd04*dzr31r31
     .             + dzreploc(19)*(r11r11*r12r12
     .                         + r21r21*r22r22)
     .             + reploc(19)*(dzr11r11*r12r12
     .                         + dzr21r21*r22r22)
     .             + reploc(19)*(r11r11*dzr12r12
     .                         + r21r21*dzr22r22)
            dzrepij(5,8)=dzrepij(5,8)
     .             + dzreploc(20)*(r11r31*r12r32
     .                         + r21r31*r22r32)
     .             + reploc(20)*(dzr11r31*r12r32
     .                         + dzr21r31*r22r32)
     .             + reploc(20)*(r11r31*dzr12r32
     .                         + r21r31*dzr22r32)
     .             + dzreploc(21)*(r11r11*r22r22
     .                         + r21r21*r12r12)
     .             + reploc(21)*(dzr11r11*r22r22
     .                         + dzr21r21*r12r12)
     .             + reploc(21)*(r11r11*dzr22r22
     .                         + r21r21*dzr12r12)
     .             + dzreploc(22)*r11r21*r12r22
     .             + reploc(22)*dzr11r21*r12r22
     .             + reploc(22)*r11r21*dzr12r22
!
!     (py py | px px):
!
        repij(8,5) = reploc(16)*r32r32*r31r31
     .             + reploc(17)*radd04*r31r31
     .             + reploc(18)*radd01*r32r32
     .             + reploc(19)*(r12r12*r11r11
     .                         + r22r22*r21r21)
     .             + reploc(20)*(r12r32*r11r31
     .                         + r22r32*r21r31)
     .             + reploc(21)*(r12r12*r21r21
     .                         + r22r22*r11r11)
     .             + reploc(22)*r12r22*r11r21

        dxrepij(8,5) = dxreploc(16)*r32r32*r31r31
     .             + reploc(16)*dxr32r32*r31r31
     .             + reploc(16)*r32r32*dxr31r31
     .             + dxreploc(17)*radd04*r31r31
     .             + reploc(17)*dxradd04*r31r31
     .             + reploc(17)*radd04*dxr31r31
     .             + dxreploc(18)*radd01*r32r32
     .             + reploc(18)*dxradd01*r32r32
     .             + reploc(18)*radd01*dxr32r32
     .             + dxreploc(19)*(r12r12*r11r11
     .                         + r22r22*r21r21)
     .             + reploc(19)*(dxr12r12*r11r11
     .                         + dxr22r22*r21r21)
     .             + reploc(19)*(r12r12*dxr11r11
     .                         + r22r22*dxr21r21)
             dxrepij(8,5)=dxrepij(8,5)
     .             + dxreploc(20)*(r12r32*r11r31
     .                         + r22r32*r21r31)
     .             + reploc(20)*(dxr12r32*r11r31
     .                         + dxr22r32*r21r31)
     .             + reploc(20)*(r12r32*dxr11r31
     .                         + r22r32*dxr21r31)
     .             + dxreploc(21)*(r12r12*r21r21
     .                         + r22r22*r11r11)
     .             + reploc(21)*(dxr12r12*r21r21
     .                         + dxr22r22*r11r11)
     .             + reploc(21)*(r12r12*dxr21r21
     .                         + r22r22*dxr11r11)
     .             + dxreploc(22)*r12r22*r11r21
     .             + reploc(22)*dxr12r22*r11r21
     .             + reploc(22)*r12r22*dxr11r21

        dyrepij(8,5) = dyreploc(16)*r32r32*r31r31
     .             + reploc(16)*dyr32r32*r31r31
     .             + reploc(16)*r32r32*dyr31r31
     .             + dyreploc(17)*radd04*r31r31
     .             + reploc(17)*dyradd04*r31r31
     .             + reploc(17)*radd04*dyr31r31
     .             + dyreploc(18)*radd01*r32r32
     .             + reploc(18)*dyradd01*r32r32
     .             + reploc(18)*radd01*dyr32r32
     .             + dyreploc(19)*(r12r12*r11r11
     .                         + r22r22*r21r21)
     .             + reploc(19)*(dyr12r12*r11r11
     .                         + dyr22r22*r21r21)
     .             + reploc(19)*(r12r12*dyr11r11
     .                         + r22r22*dyr21r21)
            dyrepij(8,5)=dyrepij(8,5)
     .             + dyreploc(20)*(r12r32*r11r31
     .                         + r22r32*r21r31)
     .             + reploc(20)*(dyr12r32*r11r31
     .                         + dyr22r32*r21r31)
     .             + reploc(20)*(r12r32*dyr11r31
     .                         + r22r32*dyr21r31)
     .             + dyreploc(21)*(r12r12*r21r21
     .                         + r22r22*r11r11)
     .             + reploc(21)*(dyr12r12*r21r21
     .                         + dyr22r22*r11r11)
     .             + reploc(21)*(r12r12*dyr21r21
     .                         + r22r22*dyr11r11)
     .             + dyreploc(22)*r12r22*r11r21
     .             + reploc(22)*dyr12r22*r11r21
     .             + reploc(22)*r12r22*dyr11r21


        dzrepij(8,5) = dzreploc(16)*r32r32*r31r31
     .             + reploc(16)*dzr32r32*r31r31
     .             + reploc(16)*r32r32*dzr31r31
     .             + dzreploc(17)*radd04*r31r31
     .             + reploc(17)*dzradd04*r31r31
     .             + reploc(17)*radd04*dzr31r31
     .             + dzreploc(18)*radd01*r32r32
     .             + reploc(18)*dzradd01*r32r32
     .             + reploc(18)*radd01*dzr32r32
     .             + dzreploc(19)*(r12r12*r11r11
     .                         + r22r22*r21r21)
     .             + reploc(19)*(dzr12r12*r11r11
     .                         + dzr22r22*r21r21)
     .             + reploc(19)*(r12r12*dzr11r11
     .                         + r22r22*dzr21r21)
               dzrepij(8,5)= dzrepij(8,5)
     .             + dzreploc(20)*(r12r32*r11r31
     .                         + r22r32*r21r31)
     .             + reploc(20)*(dzr12r32*r11r31
     .                         + dzr22r32*r21r31)
     .             + reploc(20)*(r12r32*dzr11r31
     .                         + r22r32*dzr21r31)
     .             + dzreploc(21)*(r12r12*r21r21
     .                         + r22r22*r11r11)
     .             + reploc(21)*(dzr12r12*r21r21
     .                         + dzr22r22*r11r11)
     .             + reploc(21)*(r12r12*dzr21r21
     .                         + r22r22*dzr11r11)
     .             + dzreploc(22)*r12r22*r11r21
     .             + reploc(22)*dzr12r22*r11r21
     .             + reploc(22)*r12r22*dzr11r21
!
!       (px px | py pz):
!
        repij(5,9) = reploc(16)*r31r31*r32r33
     .             + reploc(17)*radd01*r32r33
     .             + reploc(18)*radd05*r31r31
     .             + reploc(19)*(r11r11*r12r13
     .                         + r21r21*r22r23)
     .             + reploc(20)*(r11r31*r12r33
     .                         + r31r11*r32r13
     .                         + r21r31*r22r33
     .                         + r31r21*r32r23)
     .             + reploc(21)*(r11r11*r22r23
     .                         + r21r21*r12r13)
     .             + reploc(22)*(r11r21*r12r23
     .                         + r21r11*r22r13)

        dxrepij(5,9) = dxreploc(16)*r31r31*r32r33
     .             + reploc(16)*dxr31r31*r32r33
     .             + reploc(16)*r31r31*dxr32r33
     .             + dxreploc(17)*radd01*r32r33
     .             + reploc(17)*dxradd01*r32r33
     .             + reploc(17)*radd01*dxr32r33
     .             + dxreploc(18)*radd05*r31r31
     .             + reploc(18)*dxradd05*r31r31
     .             + reploc(18)*radd05*dxr31r31
     .             + dxreploc(19)*(r11r11*r12r13
     .                         + r21r21*r22r23)
     .             + reploc(19)*(dxr11r11*r12r13
     .                         + dxr21r21*r22r23)
     .             + reploc(19)*(r11r11*dxr12r13
     .                         + r21r21*dxr22r23)
     .             + dxreploc(20)*(r11r31*r12r33
     .                         + r31r11*r32r13
     .                         + r21r31*r22r33
     .                         + r31r21*r32r23)
            dxrepij(5,9)=dxrepij(5,9)
     .             + reploc(20)*(dxr11r31*r12r33
     .                         + dxr31r11*r32r13
     .                         + dxr21r31*r22r33
     .                         + dxr31r21*r32r23)
     .             + reploc(20)*(r11r31*dxr12r33
     .                         + r31r11*dxr32r13
     .                         + r21r31*dxr22r33
     .                         + r31r21*dxr32r23)
     .             + dxreploc(21)*(r11r11*r22r23
     .                         + r21r21*r12r13)
     .             + reploc(21)*(dxr11r11*r22r23
     .                         + dxr21r21*r12r13)
     .             + reploc(21)*(r11r11*dxr22r23
     .                         + r21r21*dxr12r13)
     .             + dxreploc(22)*(r11r21*r12r23
     .                         + r21r11*r22r13)
     .             + reploc(22)*(dxr11r21*r12r23
     .                         + dxr21r11*r22r13)
     .             + reploc(22)*(r11r21*dxr12r23
     .                         + r21r11*dxr22r13)

        dyrepij(5,9) = dyreploc(16)*r31r31*r32r33
     .             + reploc(16)*dyr31r31*r32r33
     .             + reploc(16)*r31r31*dyr32r33
     .             + dyreploc(17)*radd01*r32r33
     .             + reploc(17)*dyradd01*r32r33
     .             + reploc(17)*radd01*dyr32r33
     .             + dyreploc(18)*radd05*r31r31
     .             + reploc(18)*dyradd05*r31r31
     .             + reploc(18)*radd05*dyr31r31
     .             + dyreploc(19)*(r11r11*r12r13
     .                         + r21r21*r22r23)
     .             + reploc(19)*(dyr11r11*r12r13
     .                         + dyr21r21*r22r23)
     .             + reploc(19)*(r11r11*dyr12r13
     .                         + r21r21*dyr22r23)
     .             + dyreploc(20)*(r11r31*r12r33
     .                         + r31r11*r32r13
     .                         + r21r31*r22r33
     .                         + r31r21*r32r23)
             dyrepij(5,9)=dyrepij(5,9)
     .             + reploc(20)*(dyr11r31*r12r33
     .                         + dyr31r11*r32r13
     .                         + dyr21r31*r22r33
     .                         + dyr31r21*r32r23)
     .             + reploc(20)*(r11r31*dyr12r33
     .                         + r31r11*dyr32r13
     .                         + r21r31*dyr22r33
     .                         + r31r21*dyr32r23)
     .             + dyreploc(21)*(r11r11*r22r23
     .                         + r21r21*r12r13)
     .             + reploc(21)*(dyr11r11*r22r23
     .                         + dyr21r21*r12r13)
     .             + reploc(21)*(r11r11*dyr22r23
     .                         + r21r21*dyr12r13)
     .             + dyreploc(22)*(r11r21*r12r23
     .                         + r21r11*r22r13)
     .             + reploc(22)*(dyr11r21*r12r23
     .                         + dyr21r11*r22r13)
     .             + reploc(22)*(r11r21*dyr12r23
     .                         + r21r11*dyr22r13)

        dzrepij(5,9) = dzreploc(16)*r31r31*r32r33
     .             + reploc(16)*dzr31r31*r32r33
     .             + reploc(16)*r31r31*dzr32r33
     .             + dzreploc(17)*radd01*r32r33
     .             + reploc(17)*dzradd01*r32r33
     .             + reploc(17)*radd01*dzr32r33
     .             + dzreploc(18)*radd05*r31r31
     .             + reploc(18)*dzradd05*r31r31
     .             + reploc(18)*radd05*dzr31r31
     .             + dzreploc(19)*(r11r11*r12r13
     .                         + r21r21*r22r23)
     .             + reploc(19)*(dzr11r11*r12r13
     .                         + dzr21r21*r22r23)
     .             + reploc(19)*(r11r11*dzr12r13
     .                         + r21r21*dzr22r23)
     .             + dzreploc(20)*(r11r31*r12r33
     .                         + r31r11*r32r13
     .                         + r21r31*r22r33
     .                         + r31r21*r32r23)
             dzrepij(5,9)=dzrepij(5,9)
     .             + reploc(20)*(dzr11r31*r12r33
     .                         + dzr31r11*r32r13
     .                         + dzr21r31*r22r33
     .                         + dzr31r21*r32r23)
     .             + reploc(20)*(r11r31*dzr12r33
     .                         + r31r11*dzr32r13
     .                         + r21r31*dzr22r33
     .                         + r31r21*dzr32r23)
     .             + dzreploc(21)*(r11r11*r22r23
     .                         + r21r21*r12r13)
     .             + reploc(21)*(dzr11r11*r22r23
     .                         + dzr21r21*r12r13)
     .             + reploc(21)*(r11r11*dzr22r23
     .                         + r21r21*dzr12r13)
     .             + dzreploc(22)*(r11r21*r12r23
     .                         + r21r11*r22r13)
     .             + reploc(22)*(dzr11r21*r12r23
     .                         + dzr21r11*r22r13)
     .             + reploc(22)*(r11r21*dzr12r23
     .                         + r21r11*dzr22r13)
!
!       (py pz | px px):
!
        repij(9,5) = reploc(16)*r32r33*r31r31
     .             + reploc(17)*radd05*r31r31
     .             + reploc(18)*radd01*r32r33
     .             + reploc(19)*(r12r13*r11r11
     .                         + r22r23*r21r21)
     .             + reploc(20)*(radd13*r11r31
     .                         + radd14*r21r31)
     .             + reploc(21)*(r12r13*r21r21
     .                         + r22r23*r11r11)
     .             + reploc(22)*radd16*r11r21

        dxrepij(9,5) = dxreploc(16)*r32r33*r31r31
     .             + reploc(16)*dxr32r33*r31r31
     .             + reploc(16)*r32r33*dxr31r31
     .             + dxreploc(17)*radd05*r31r31
     .             + reploc(17)*dxradd05*r31r31
     .             + reploc(17)*radd05*dxr31r31
     .             + dxreploc(18)*radd01*r32r33
     .             + reploc(18)*dxradd01*r32r33
     .             + reploc(18)*radd01*dxr32r33
     .             + dxreploc(19)*(r12r13*r11r11
     .                         + r22r23*r21r21)
     .             + reploc(19)*(dxr12r13*r11r11
     .                         + dxr22r23*r21r21)
     .             + reploc(19)*(r12r13*dxr11r11
     .                         + r22r23*dxr21r21)
                dxrepij(9,5)=dxrepij(9,5)
     .             + dxreploc(20)*(radd13*r11r31
     .                         + radd14*r21r31)
     .             + reploc(20)*(dxradd13*r11r31
     .                         + dxradd14*r21r31)
     .             + reploc(20)*(radd13*dxr11r31
     .                         + radd14*dxr21r31)
     .             + dxreploc(21)*(r12r13*r21r21
     .                         + r22r23*r11r11)
     .             + reploc(21)*(dxr12r13*r21r21
     .                         + dxr22r23*r11r11)
     .             + reploc(21)*(r12r13*dxr21r21
     .                         + r22r23*dxr11r11)
     .             + dxreploc(22)*radd16*r11r21
     .             + reploc(22)*dxradd16*r11r21
     .             + reploc(22)*radd16*dxr11r21

        dyrepij(9,5) = dyreploc(16)*r32r33*r31r31
     .             + reploc(16)*dyr32r33*r31r31
     .             + reploc(16)*r32r33*dyr31r31
     .             + dyreploc(17)*radd05*r31r31
     .             + reploc(17)*dyradd05*r31r31
     .             + reploc(17)*radd05*dyr31r31
     .             + dyreploc(18)*radd01*r32r33
     .             + reploc(18)*dyradd01*r32r33
     .             + reploc(18)*radd01*dyr32r33
     .             + dyreploc(19)*(r12r13*r11r11
     .                         + r22r23*r21r21)
     .             + reploc(19)*(dyr12r13*r11r11
     .                         + dyr22r23*r21r21)
     .             + reploc(19)*(r12r13*dyr11r11
     .                         + r22r23*dyr21r21)
             dyrepij(9,5)=dyrepij(9,5)
     .             + dyreploc(20)*(radd13*r11r31
     .                         + radd14*r21r31)
     .             + reploc(20)*(dyradd13*r11r31
     .                         + dyradd14*r21r31)
     .             + reploc(20)*(radd13*dyr11r31
     .                         + radd14*dyr21r31)
     .             + dyreploc(21)*(r12r13*r21r21
     .                         + r22r23*r11r11)
     .             + reploc(21)*(dyr12r13*r21r21
     .                         + dyr22r23*r11r11)
     .             + reploc(21)*(r12r13*dyr21r21
     .                         + r22r23*dyr11r11)
     .             + dyreploc(22)*radd16*r11r21
     .             + reploc(22)*dyradd16*r11r21
     .             + reploc(22)*radd16*dyr11r21

        dzrepij(9,5) = dzreploc(16)*r32r33*r31r31
     .             + reploc(16)*dzr32r33*r31r31
     .             + reploc(16)*r32r33*dzr31r31
     .             + dzreploc(17)*radd05*r31r31
     .             + reploc(17)*dzradd05*r31r31
     .             + reploc(17)*radd05*dzr31r31
     .             + dzreploc(18)*radd01*r32r33
     .             + reploc(18)*dzradd01*r32r33
     .             + reploc(18)*radd01*dzr32r33
     .             + dzreploc(19)*(r12r13*r11r11
     .                         + r22r23*r21r21)
     .             + reploc(19)*(dzr12r13*r11r11
     .                         + dzr22r23*r21r21)
     .             + reploc(19)*(r12r13*dzr11r11
     .                         + r22r23*dzr21r21)
               dzrepij(9,5)=dzrepij(9,5)
     .             + dzreploc(20)*(radd13*r11r31
     .                         + radd14*r21r31)
     .             + reploc(20)*(dzradd13*r11r31
     .                         + dzradd14*r21r31)
     .             + reploc(20)*(radd13*dzr11r31
     .                         + radd14*dzr21r31)
     .             + dzreploc(21)*(r12r13*r21r21
     .                         + r22r23*r11r11)
     .             + reploc(21)*(dzr12r13*r21r21
     .                         + dzr22r23*r11r11)
     .             + reploc(21)*(r12r13*dzr21r21
     .                         + r22r23*dzr11r11)
     .             + dzreploc(22)*radd16*r11r21
     .             + reploc(22)*dzradd16*r11r21
     .             + reploc(22)*radd16*dzr11r21
!
!       (px px | pz pz):
!
        repij(5,10) = reploc(16)*r31r31*r33r33
     .              + reploc(17)*radd01*r33r33
     .              + reploc(18)*radd06*r31r31
     .              + reploc(19)*(r11r11*r13r13
     .                          + r21r21*r23r23)
     .              + reploc(20)*(r11r31*r13r33
     .                          + r21r31*r23r33)
     .              + reploc(21)*(r11r11*r23r23
     .                          + r21r21*r13r13)
     .              + reploc(22)*r11r21*r13r23

        dxrepij(5,10) = dxreploc(16)*r31r31*r33r33
     .              + reploc(16)*dxr31r31*r33r33
     .              + reploc(16)*r31r31*dxr33r33
     .              + dxreploc(17)*radd01*r33r33
     .              + reploc(17)*dxradd01*r33r33
     .              + reploc(17)*radd01*dxr33r33
     .              + dxreploc(18)*radd06*r31r31
     .              + reploc(18)*dxradd06*r31r31
     .              + reploc(18)*radd06*dxr31r31
     .              + dxreploc(19)*(r11r11*r13r13
     .                          + r21r21*r23r23)
     .              + reploc(19)*(dxr11r11*r13r13
     .                          + dxr21r21*r23r23)
     .              + reploc(19)*(r11r11*dxr13r13
     .                          + r21r21*dxr23r23)
               dxrepij(5,10)=dxrepij(5,10)
     .              + dxreploc(20)*(r11r31*r13r33
     .                          + r21r31*r23r33)
     .              + reploc(20)*(dxr11r31*r13r33
     .                          + dxr21r31*r23r33)
     .              + reploc(20)*(r11r31*dxr13r33
     .                          + r21r31*dxr23r33)
     .              + dxreploc(21)*(r11r11*r23r23
     .                          + r21r21*r13r13)
     .              + reploc(21)*(dxr11r11*r23r23
     .                          + dxr21r21*r13r13)
     .              + reploc(21)*(r11r11*dxr23r23
     .                          + r21r21*dxr13r13)
     .              + dxreploc(22)*r11r21*r13r23
     .              + reploc(22)*dxr11r21*r13r23
     .              + reploc(22)*r11r21*dxr13r23

        dyrepij(5,10) = dyreploc(16)*r31r31*r33r33
     .              + reploc(16)*dyr31r31*r33r33
     .              + reploc(16)*r31r31*dyr33r33
     .              + dyreploc(17)*radd01*r33r33
     .              + reploc(17)*dyradd01*r33r33
     .              + reploc(17)*radd01*dyr33r33
     .              + dyreploc(18)*radd06*r31r31
     .              + reploc(18)*dyradd06*r31r31
     .              + reploc(18)*radd06*dyr31r31
     .              + dyreploc(19)*(r11r11*r13r13
     .                          + r21r21*r23r23)
     .              + reploc(19)*(dyr11r11*r13r13
     .                          + dyr21r21*r23r23)
     .              + reploc(19)*(r11r11*dyr13r13
     .                          + r21r21*dyr23r23)
                 dyrepij(5,10)=dyrepij(5,10)
     .              + dyreploc(20)*(r11r31*r13r33
     .                          + r21r31*r23r33)
     .              + reploc(20)*(dyr11r31*r13r33
     .                          + dyr21r31*r23r33)
     .              + reploc(20)*(r11r31*dyr13r33
     .                          + r21r31*dyr23r33)
     .              + dyreploc(21)*(r11r11*r23r23
     .                          + r21r21*r13r13)
     .              + reploc(21)*(dyr11r11*r23r23
     .                          + dyr21r21*r13r13)
     .              + reploc(21)*(r11r11*dyr23r23
     .                          + r21r21*dyr13r13)
     .              + dyreploc(22)*r11r21*r13r23
     .              + reploc(22)*dyr11r21*r13r23
     .              + reploc(22)*r11r21*dyr13r23

        dzrepij(5,10) = dzreploc(16)*r31r31*r33r33
     .              + reploc(16)*dzr31r31*r33r33
     .              + reploc(16)*r31r31*dzr33r33
     .              + dzreploc(17)*radd01*r33r33
     .              + reploc(17)*dzradd01*r33r33
     .              + reploc(17)*radd01*dzr33r33
     .              + dzreploc(18)*radd06*r31r31
     .              + reploc(18)*dzradd06*r31r31
     .              + reploc(18)*radd06*dzr31r31
     .              + dzreploc(19)*(r11r11*r13r13
     .                          + r21r21*r23r23)
     .              + reploc(19)*(dzr11r11*r13r13
     .                          + dzr21r21*r23r23)
     .              + reploc(19)*(r11r11*dzr13r13
     .                          + r21r21*dzr23r23)
               dzrepij(5,10)=dzrepij(5,10)
     .              + dzreploc(20)*(r11r31*r13r33
     .                          + r21r31*r23r33)
     .              + reploc(20)*(dzr11r31*r13r33
     .                          + dzr21r31*r23r33)
     .              + reploc(20)*(r11r31*dzr13r33
     .                          + r21r31*dzr23r33)
     .              + dzreploc(21)*(r11r11*r23r23
     .                          + r21r21*r13r13)
     .              + reploc(21)*(dzr11r11*r23r23
     .                          + dzr21r21*r13r13)
     .              + reploc(21)*(r11r11*dzr23r23
     .                          + r21r21*dzr13r13)
     .              + dzreploc(22)*r11r21*r13r23
     .              + reploc(22)*dzr11r21*r13r23
     .              + reploc(22)*r11r21*dzr13r23
!
!       (pz pz | px px):
!
        repij(10,5) = reploc(16)*r33r33*r31r31
     .              + reploc(17)*radd06*r31r31
     .              + reploc(18)*radd01*r33r33
     .              + reploc(19)*(r13r13*r11r11
     .                          + r23r23*r21r21)
     .              + reploc(20)*(r13r33*r11r31
     .                          + r23r33*r21r31)
     .              + reploc(21)*(r13r13*r21r21
     .                          + r23r23*r11r11)
     .              + reploc(22)*r13r23*r11r21

        dxrepij(10,5) = dxreploc(16)*r33r33*r31r31
     .              + reploc(16)*dxr33r33*r31r31
     .              + reploc(16)*r33r33*dxr31r31
     .              + dxreploc(17)*radd06*r31r31
     .              + reploc(17)*dxradd06*r31r31
     .              + reploc(17)*radd06*dxr31r31
     .              + dxreploc(18)*radd01*r33r33
     .              + reploc(18)*dxradd01*r33r33
     .              + reploc(18)*radd01*dxr33r33
     .              + dxreploc(19)*(r13r13*r11r11
     .                          + r23r23*r21r21)
     .              + reploc(19)*(dxr13r13*r11r11
     .                          + dxr23r23*r21r21)
     .              + reploc(19)*(r13r13*dxr11r11
     .                          + r23r23*dxr21r21)
                 dxrepij(10,5)=dxrepij(10,5)
     .              + dxreploc(20)*(r13r33*r11r31
     .                          + r23r33*r21r31)
     .              + reploc(20)*(dxr13r33*r11r31
     .                          + dxr23r33*r21r31)
     .              + reploc(20)*(r13r33*dxr11r31
     .                          + r23r33*dxr21r31)
     .              + dxreploc(21)*(r13r13*r21r21
     .                          + r23r23*r11r11)
     .              + reploc(21)*(dxr13r13*r21r21
     .                          + dxr23r23*r11r11)
     .              + reploc(21)*(r13r13*dxr21r21
     .                          + r23r23*dxr11r11)
     .              + dxreploc(22)*r13r23*r11r21
     .              + reploc(22)*dxr13r23*r11r21
     .              + reploc(22)*r13r23*dxr11r21

        dyrepij(10,5) = dyreploc(16)*r33r33*r31r31
     .              + reploc(16)*dyr33r33*r31r31
     .              + reploc(16)*r33r33*dyr31r31
     .              + dyreploc(17)*radd06*r31r31
     .              + reploc(17)*dyradd06*r31r31
     .              + reploc(17)*radd06*dyr31r31
     .              + dyreploc(18)*radd01*r33r33
     .              + reploc(18)*dyradd01*r33r33
     .              + reploc(18)*radd01*dyr33r33
     .              + dyreploc(19)*(r13r13*r11r11
     .                          + r23r23*r21r21)
     .              + reploc(19)*(dyr13r13*r11r11
     .                          + dyr23r23*r21r21)
     .              + reploc(19)*(r13r13*dyr11r11
     .                          + r23r23*dyr21r21)
                dyrepij(10,5)=dyrepij(10,5)
     .              + dyreploc(20)*(r13r33*r11r31
     .                          + r23r33*r21r31)
     .              + reploc(20)*(dyr13r33*r11r31
     .                          + dyr23r33*r21r31)
     .              + reploc(20)*(r13r33*dyr11r31
     .                          + r23r33*dyr21r31)
     .              + dyreploc(21)*(r13r13*r21r21
     .                          + r23r23*r11r11)
     .              + reploc(21)*(dyr13r13*r21r21
     .                          + dyr23r23*r11r11)
     .              + reploc(21)*(r13r13*dyr21r21
     .                          + r23r23*dyr11r11)
     .              + dyreploc(22)*r13r23*r11r21
     .              + reploc(22)*dyr13r23*r11r21
     .              + reploc(22)*r13r23*dyr11r21

        dzrepij(10,5) = dzreploc(16)*r33r33*r31r31
     .              + reploc(16)*dzr33r33*r31r31
     .              + reploc(16)*r33r33*dzr31r31
     .              + dzreploc(17)*radd06*r31r31
     .              + reploc(17)*dzradd06*r31r31
     .              + reploc(17)*radd06*dzr31r31
     .              + dzreploc(18)*radd01*r33r33
     .              + reploc(18)*dzradd01*r33r33
     .              + reploc(18)*radd01*dzr33r33
     .              + dzreploc(19)*(r13r13*r11r11
     .                          + r23r23*r21r21)
     .              + reploc(19)*(dzr13r13*r11r11
     .                          + dzr23r23*r21r21)
     .              + reploc(19)*(r13r13*dzr11r11
     .                          + r23r23*dzr21r21)
                dzrepij(10,5)=dzrepij(10,5)
     .              + dzreploc(20)*(r13r33*r11r31
     .                          + r23r33*r21r31)
     .              + reploc(20)*(dzr13r33*r11r31
     .                          + dzr23r33*r21r31)
     .              + reploc(20)*(r13r33*dzr11r31
     .                          + r23r33*dzr21r31)
     .              + dzreploc(21)*(r13r13*r21r21
     .                          + r23r23*r11r11)
     .              + reploc(21)*(dzr13r13*r21r21
     .                          + dzr23r23*r11r11)
     .              + reploc(21)*(r13r13*dzr21r21
     .                          + r23r23*dzr11r11)
     .              + dzreploc(22)*r13r23*r11r21
     .              + reploc(22)*dzr13r23*r11r21
     .              + reploc(22)*r13r23*dzr11r21
!
!       (px py | px py):
!
        repij(6,6) = reploc(16)*r31r32*r31r32
     .             + reploc(17)*radd02*r31r32
     .             + reploc(18)*radd02*r31r32
     .             + reploc(19)*(r11r12*r11r12
     .                         + r21r22*r21r22)
     .             + reploc(20)*(radd11**2
     .                         + radd12**2)
     .             + reploc(21)*2.0d0*r11r12*r21r22
     .             + reploc(22)*radd17**2

        dxrepij(6,6) = dxreploc(16)*r31r32*r31r32
     .             + reploc(16)*dxr31r32*r31r32
     .             + reploc(16)*r31r32*dxr31r32
     .             + dxreploc(17)*radd02*r31r32
     .             + reploc(17)*dxradd02*r31r32
     .             + reploc(17)*radd02*dxr31r32
     .             + dxreploc(18)*radd02*r31r32
     .             + reploc(18)*dxradd02*r31r32
     .             + reploc(18)*radd02*dxr31r32
     .             + dxreploc(19)*(r11r12*r11r12
     .                         + r21r22*r21r22)
     .             + reploc(19)*(dxr11r12*r11r12
     .                         + dxr21r22*r21r22)
               dxrepij(6,6)=dxrepij(6,6)
     .             + reploc(19)*(r11r12*dxr11r12
     .                         + r21r22*dxr21r22)
     .             + dxreploc(20)*(radd11**2
     .                         + radd12**2)
     .             + reploc(20)*2.0d0*(radd11*dxradd11
     .                         + radd12*dxradd12)
     .             + dxreploc(21)*2.0d0*r11r12*r21r22
     .             + reploc(21)*2.0d0*dxr11r12*r21r22
     .             + reploc(21)*2.0d0*r11r12*dxr21r22
     .             + dxreploc(22)*radd17**2
     .             + reploc(22)*radd17*dxradd17*2.0d0

        dyrepij(6,6) = dyreploc(16)*r31r32*r31r32
     .             + reploc(16)*dyr31r32*r31r32
     .             + reploc(16)*r31r32*dyr31r32
     .             + dyreploc(17)*radd02*r31r32
     .             + reploc(17)*dyradd02*r31r32
     .             + reploc(17)*radd02*dyr31r32
     .             + dyreploc(18)*radd02*r31r32
     .             + reploc(18)*dyradd02*r31r32
     .             + reploc(18)*radd02*dyr31r32
     .             + dyreploc(19)*(r11r12*r11r12
     .                         + r21r22*r21r22)
     .             + reploc(19)*(dyr11r12*r11r12
     .                         + dyr21r22*r21r22)
                 dyrepij(6,6)=dyrepij(6,6)
     .             + reploc(19)*(r11r12*dyr11r12
     .                         + r21r22*dyr21r22)
     .             + dyreploc(20)*(radd11**2
     .                         + radd12**2)
     .             + reploc(20)*2.0d0*(radd11*dyradd11
     .                         + radd12*dyradd12)
     .             + dyreploc(21)*2.0d0*r11r12*r21r22
     .             + reploc(21)*2.0d0*dyr11r12*r21r22
     .             + reploc(21)*2.0d0*r11r12*dyr21r22
     .             + dyreploc(22)*radd17**2
     .             + reploc(22)*2.0d0*radd17*dyradd17
     .        

        dzrepij(6,6) = dzreploc(16)*r31r32*r31r32
     .             + reploc(16)*dzr31r32*r31r32
     .             + reploc(16)*r31r32*dzr31r32
     .             + dzreploc(17)*radd02*r31r32
     .             + reploc(17)*dzradd02*r31r32
     .             + reploc(17)*radd02*dzr31r32
     .             + dzreploc(18)*radd02*r31r32
     .             + reploc(18)*dzradd02*r31r32
     .             + reploc(18)*radd02*dzr31r32
     .             + dzreploc(19)*(r11r12*r11r12
     .                         + r21r22*r21r22)
     .             + reploc(19)*(dzr11r12*r11r12
     .                         + dzr21r22*r21r22)
                dzrepij(6,6)=dzrepij(6,6)
     .             + reploc(19)*(r11r12*dzr11r12
     .                         + r21r22*dzr21r22)
     .             + dzreploc(20)*(radd11**2
     .                         + radd12**2)
     .             + reploc(20)*2.0d0*(radd11*dzradd11
     .                         + radd12*dzradd12)
     .             + dzreploc(21)*2.0d0*r11r12*r21r22
     .             + reploc(21)*2.0d0*dzr11r12*r21r22
     .             + reploc(21)*2.0d0*r11r12*dzr21r22
     .             + dzreploc(22)*radd17**2
     .             + reploc(22)*2.0d0*radd17*dzradd17
!
!       (px py | px pz):
!
        repij(6,7) = reploc(16)*r31r32*r31r33
     .             + reploc(17)*radd02*r31r33
     .             + reploc(18)*r31r32*radd03
     .             + reploc(19)*(r11r12*r11r13
     .                         + r21r22*r21r23)
     .             + reploc(20)*(radd11*radd09
     .                         + radd12*radd10)
     .             + reploc(21)*(r11r12*r21r23
     .                         + r21r22*r11r13)
     .             + reploc(22)*radd17*(r11r23+r21r13)

        dxrepij(6,7) = dxreploc(16)*r31r32*r31r33
     .             + reploc(16)*dxr31r32*r31r33
     .             + reploc(16)*r31r32*dxr31r33
     .             + dxreploc(17)*radd02*r31r33
     .             + reploc(17)*dxradd02*r31r33
     .             + reploc(17)*radd02*dxr31r33
     .             + dxreploc(18)*r31r32*radd03
     .             + reploc(18)*dxr31r32*radd03
     .             + reploc(18)*r31r32*dxradd03
     .             + dxreploc(19)*(r11r12*r11r13
     .                         + r21r22*r21r23)
     .             + reploc(19)*(dxr11r12*r11r13
     .                         + dxr21r22*r21r23)
     .             + reploc(19)*(r11r12*dxr11r13
     .                         + r21r22*dxr21r23)
                dxrepij(6,7)=dxrepij(6,7)
     .             + dxreploc(20)*(radd11*radd09
     .                         + radd12*radd10)
     .             + reploc(20)*(dxradd11*radd09
     .                         + dxradd12*radd10)
     .             + reploc(20)*(radd11*dxradd09
     .                         + radd12*dxradd10)
     .             + dxreploc(21)*(r11r12*r21r23
     .                         + r21r22*r11r13)
     .             + reploc(21)*(dxr11r12*r21r23
     .                         + dxr21r22*r11r13)
     .             + reploc(21)*(r11r12*dxr21r23
     .                         + r21r22*dxr11r13)
     .             + dxreploc(22)*radd17*(r11r23+r21r13)
     .             + reploc(22)*dxradd17*(r11r23+r21r13)
     .             + reploc(22)*radd17*(dxr11r23+dxr21r13)

        dyrepij(6,7) = dyreploc(16)*r31r32*r31r33
     .             + reploc(16)*dyr31r32*r31r33
     .             + reploc(16)*r31r32*dyr31r33
     .             + dyreploc(17)*radd02*r31r33
     .             + reploc(17)*dyradd02*r31r33
     .             + reploc(17)*radd02*dyr31r33
     .             + dyreploc(18)*r31r32*radd03
     .             + reploc(18)*dyr31r32*radd03
     .             + reploc(18)*r31r32*dyradd03
     .             + dyreploc(19)*(r11r12*r11r13
     .                         + r21r22*r21r23)
     .             + reploc(19)*(dyr11r12*r11r13
     .                         + dyr21r22*r21r23)
     .             + reploc(19)*(r11r12*dyr11r13
     .                         + r21r22*dyr21r23)
                dyrepij(6,7)=dyrepij(6,7) 
     .             + dyreploc(20)*(radd11*radd09
     .                         + radd12*radd10)
     .             + reploc(20)*(dyradd11*radd09
     .                         + dyradd12*radd10)
     .             + reploc(20)*(radd11*dyradd09
     .                         + radd12*dyradd10)
     .             + dyreploc(21)*(r11r12*r21r23
     .                         + r21r22*r11r13)
     .             + reploc(21)*(dyr11r12*r21r23
     .                         + dyr21r22*r11r13)
     .             + reploc(21)*(r11r12*dyr21r23
     .                         + r21r22*dyr11r13)
     .             + dyreploc(22)*radd17*(r11r23+r21r13)
     .             + reploc(22)*dyradd17*(r11r23+r21r13)
     .             + reploc(22)*radd17*(dyr11r23+dyr21r13)

        dzrepij(6,7) = dzreploc(16)*r31r32*r31r33
     .             + reploc(16)*dzr31r32*r31r33
     .             + reploc(16)*r31r32*dzr31r33
     .             + dzreploc(17)*radd02*r31r33
     .             + reploc(17)*dzradd02*r31r33
     .             + reploc(17)*radd02*dzr31r33
     .             + dzreploc(18)*r31r32*radd03
     .             + reploc(18)*dzr31r32*radd03
     .             + reploc(18)*r31r32*dzradd03
     .             + dzreploc(19)*(r11r12*r11r13
     .                         + r21r22*r21r23)
     .             + reploc(19)*(dzr11r12*r11r13
     .                         + dzr21r22*r21r23)
     .             + reploc(19)*(r11r12*dzr11r13
     .                         + r21r22*dzr21r23)
               dzrepij(6,7)=dzrepij(6,7)
     .             + dzreploc(20)*(radd11*radd09
     .                         + radd12*radd10)
     .             + reploc(20)*(dzradd11*radd09
     .                         + dzradd12*radd10)
     .             + reploc(20)*(radd11*dzradd09
     .                         + radd12*dzradd10)
     .             + dzreploc(21)*(r11r12*r21r23
     .                         + r21r22*r11r13)
     .             + reploc(21)*(dzr11r12*r21r23
     .                         + dzr21r22*r11r13)
     .             + reploc(21)*(r11r12*dzr21r23
     .                         + r21r22*dzr11r13)
     .             + dzreploc(22)*radd17*(r11r23+r21r13)
     .             + reploc(22)*dzradd17*(r11r23+r21r13)
     .             + reploc(22)*radd17*(dzr11r23+dzr21r13)
!
!       (px pz | px py):
!
        repij(7,6) = repij(6,7)
         dxrepij(7,6) = dxrepij(6,7)
         dyrepij(7,6) = dyrepij(6,7)
         dzrepij(7,6) = dzrepij(6,7)
        
!
!       (px py | py py):
!
        repij(6,8) = reploc(16)*r31r32*r32r32
     .             + reploc(17)*radd02*r32r32
     .             + reploc(18)*radd04*r31r32
     .             + reploc(19)*(r11r12*r12r12
     .                         + r21r22*r22r22)
     .             + reploc(20)*(radd11*r12r32
     .                         + radd12*r22r32)
     .             + reploc(21)*(r11r12*r22r22
     .                         + r21r22*r12r12)
     .             + reploc(22)*radd17*r12r22

        dxrepij(6,8) = dxreploc(16)*r31r32*r32r32
     .             + reploc(16)*dxr31r32*r32r32
     .             + reploc(16)*r31r32*dxr32r32
     .             + dxreploc(17)*radd02*r32r32
     .             + reploc(17)*dxradd02*r32r32
     .             + reploc(17)*radd02*dxr32r32
     .             + dxreploc(18)*radd04*r31r32
     .             + reploc(18)*dxradd04*r31r32
     .             + reploc(18)*radd04*dxr31r32
     .             + dxreploc(19)*(r11r12*r12r12
     .                         + r21r22*r22r22)
     .             + reploc(19)*(dxr11r12*r12r12
     .                         + dxr21r22*r22r22)
     .             + reploc(19)*(r11r12*dxr12r12
     .                         + r21r22*dxr22r22)
                  dxrepij(6,8)=dxrepij(6,8)
     .             + dxreploc(20)*(radd11*r12r32
     .                         + radd12*r22r32)
     .             + reploc(20)*(dxradd11*r12r32
     .                         + dxradd12*r22r32)
     .             + reploc(20)*(radd11*dxr12r32
     .                         + radd12*dxr22r32)
     .             + dxreploc(21)*(r11r12*r22r22
     .                         + r21r22*r12r12)
     .             + reploc(21)*(dxr11r12*r22r22
     .                         + dxr21r22*r12r12)
     .             + reploc(21)*(r11r12*dxr22r22
     .                         + r21r22*dxr12r12)
     .             + dxreploc(22)*radd17*r12r22
     .             + reploc(22)*dxradd17*r12r22
     .             + reploc(22)*radd17*dxr12r22


        dyrepij(6,8) = dyreploc(16)*r31r32*r32r32
     .             + reploc(16)*dyr31r32*r32r32
     .             + reploc(16)*r31r32*dyr32r32
     .             + dyreploc(17)*radd02*r32r32
     .             + reploc(17)*dyradd02*r32r32
     .             + reploc(17)*radd02*dyr32r32
     .             + dyreploc(18)*radd04*r31r32
     .             + reploc(18)*dyradd04*r31r32
     .             + reploc(18)*radd04*dyr31r32
     .             + dyreploc(19)*(r11r12*r12r12
     .                         + r21r22*r22r22)
     .             + reploc(19)*(dyr11r12*r12r12
     .                         + dyr21r22*r22r22)
     .             + reploc(19)*(r11r12*dyr12r12
     .                         + r21r22*dyr22r22)
               dyrepij(6,8)=dyrepij(6,8)
     .             + dyreploc(20)*(radd11*r12r32
     .                         + radd12*r22r32)
     .             + reploc(20)*(dyradd11*r12r32
     .                         + dyradd12*r22r32)
     .             + reploc(20)*(radd11*dyr12r32
     .                         + radd12*dyr22r32)
     .             + dyreploc(21)*(r11r12*r22r22
     .                         + r21r22*r12r12)
     .             + reploc(21)*(dyr11r12*r22r22
     .                         + dyr21r22*r12r12)
     .             + reploc(21)*(r11r12*dyr22r22
     .                         + r21r22*dyr12r12)
     .             + dyreploc(22)*radd17*r12r22
     .             + reploc(22)*dyradd17*r12r22
     .             + reploc(22)*radd17*dyr12r22

        dzrepij(6,8) = dzreploc(16)*r31r32*r32r32
     .             + reploc(16)*dzr31r32*r32r32
     .             + reploc(16)*r31r32*dzr32r32
     .             + dzreploc(17)*radd02*r32r32
     .             + reploc(17)*dzradd02*r32r32
     .             + reploc(17)*radd02*dzr32r32
     .             + dzreploc(18)*radd04*r31r32
     .             + reploc(18)*dzradd04*r31r32
     .             + reploc(18)*radd04*dzr31r32
     .             + dzreploc(19)*(r11r12*r12r12
     .                         + r21r22*r22r22)
     .             + reploc(19)*(dzr11r12*r12r12
     .                         + dzr21r22*r22r22)
     .             + reploc(19)*(r11r12*dzr12r12
     .             + r21r22*dzr22r22)
                 dzrepij(6,8)=dzrepij(6,8)
     .             + dzreploc(20)*(radd11*r12r32
     .                         + radd12*r22r32)
     .             + reploc(20)*(dzradd11*r12r32
     .                         + dzradd12*r22r32)
     .             + reploc(20)*(radd11*dzr12r32
     .                         + radd12*dzr22r32)
     .             + dzreploc(21)*(r11r12*r22r22
     .                         + r21r22*r12r12)
     .             + reploc(21)*(dzr11r12*r22r22
     .                         + dzr21r22*r12r12)
     .             + reploc(21)*(r11r12*dzr22r22
     .                         + r21r22*dzr12r12)
     .             + dzreploc(22)*radd17*r12r22
     .             + reploc(22)*dzradd17*r12r22
     .             + reploc(22)*radd17*dzr12r22
!
!       (py py | px py):
!
        repij(8,6) = reploc(16)*r32r32*r31r32
     .             + reploc(17)*radd04*r31r32
     .             + reploc(18)*radd02*r32r32
     .             + reploc(19)*(r12r12*r11r12
     .                         + r22r22*r21r22)
     .             + reploc(20)*(r12r32*r11r32
     .                         + r32r12*r31r12
     .                         + r22r32*r21r32
     .                         + r32r22*r31r22)
     .             + reploc(21)*(r12r12*r21r22
     .                         + r22r22*r11r12)
     .             + reploc(22)*(r12r22*r11r22
     .                         + r22r12*r21r12)

        dxrepij(8,6) = dxreploc(16)*r32r32*r31r32
     .             + reploc(16)*dxr32r32*r31r32
     .             + reploc(16)*r32r32*dxr31r32
     .             + dxreploc(17)*radd04*r31r32
     .             + reploc(17)*dxradd04*r31r32
     .             + reploc(17)*radd04*dxr31r32
     .             + dxreploc(18)*radd02*r32r32
     .             + reploc(18)*dxradd02*r32r32
     .             + reploc(18)*radd02*dxr32r32
     .             + dxreploc(19)*(r12r12*r11r12
     .                         + r22r22*r21r22)
     .             + reploc(19)*(dxr12r12*r11r12
     .                         + dxr22r22*r21r22)
     .             + reploc(19)*(r12r12*dxr11r12
     .                         + r22r22*dxr21r22)
     .             + dxreploc(20)*(r12r32*r11r32
     .                         + r32r12*r31r12
     .                         + r22r32*r21r32
     .                         + r32r22*r31r22)
                    dxrepij(8,6)=dxrepij(8,6)
     .             + reploc(20)*(dxr12r32*r11r32
     .                         + dxr32r12*r31r12
     .                         + dxr22r32*r21r32
     .                         + dxr32r22*r31r22)
     .             + reploc(20)*(r12r32*dxr11r32
     .                         + r32r12*dxr31r12
     .                         + r22r32*dxr21r32
     .                         + r32r22*dxr31r22)
     .             + dxreploc(21)*(r12r12*r21r22
     .                         + r22r22*r11r12)
     .             + reploc(21)*(dxr12r12*r21r22
     .                         + dxr22r22*r11r12)
     .             + reploc(21)*(r12r12*dxr21r22
     .                         + r22r22*dxr11r12)
     .             + dxreploc(22)*(r12r22*r11r22+ r22r12*r21r12)
     .             + reploc(22)*(dxr12r22*r11r22+ dxr22r12*r21r12)
     .             + reploc(22)*(r12r22*dxr11r22+ r22r12*dxr21r12)

        dyrepij(8,6) = dyreploc(16)*r32r32*r31r32
     .             + reploc(16)*dyr32r32*r31r32
     .             + reploc(16)*r32r32*dyr31r32
     .             + dyreploc(17)*radd04*r31r32
     .             + reploc(17)*dyradd04*r31r32
     .             + reploc(17)*radd04*dyr31r32
     .             + dyreploc(18)*radd02*r32r32
     .             + reploc(18)*dyradd02*r32r32
     .             + reploc(18)*radd02*dyr32r32
     .             + dyreploc(19)*(r12r12*r11r12
     .                         + r22r22*r21r22)
     .             + reploc(19)*(dyr12r12*r11r12
     .                         + dyr22r22*r21r22)
     .             + reploc(19)*(r12r12*dyr11r12
     .                         + r22r22*dyr21r22)
     .             + dyreploc(20)*(r12r32*r11r32
     .                         + r32r12*r31r12
     .                         + r22r32*r21r32
     .                         + r32r22*r31r22)
                    dyrepij(8,6)=dyrepij(8,6)
     .             + reploc(20)*(dyr12r32*r11r32
     .                         + dyr32r12*r31r12
     .                         + dyr22r32*r21r32
     .                         + dyr32r22*r31r22)
     .             + reploc(20)*(r12r32*dyr11r32
     .                         + r32r12*dyr31r12
     .                         + r22r32*dyr21r32
     .                         + r32r22*dyr31r22)
     .             + dyreploc(21)*(r12r12*r21r22
     .                         + r22r22*r11r12)
     .             + reploc(21)*(dyr12r12*r21r22
     .                         + dyr22r22*r11r12)
     .             + reploc(21)*(r12r12*dyr21r22
     .                         + r22r22*dyr11r12)
     .             + dyreploc(22)*(r12r22*r11r22+ r22r12*r21r12)
     .             + reploc(22)*(dyr12r22*r11r22+ dyr22r12*r21r12)
     .             + reploc(22)*(r12r22*dyr11r22+ r22r12*dyr21r12)

        dzrepij(8,6) = dzreploc(16)*r32r32*r31r32
     .             + reploc(16)*dzr32r32*r31r32
     .             + reploc(16)*r32r32*dzr31r32
     .             + dzreploc(17)*radd04*r31r32
     .             + reploc(17)*dzradd04*r31r32
     .             + reploc(17)*radd04*dzr31r32
     .             + dzreploc(18)*radd02*r32r32
     .             + reploc(18)*dzradd02*r32r32
     .             + reploc(18)*radd02*dzr32r32
     .             + dzreploc(19)*(r12r12*r11r12
     .                         + r22r22*r21r22)
     .             + reploc(19)*(dzr12r12*r11r12
     .                         + dzr22r22*r21r22)
     .             + reploc(19)*(r12r12*dzr11r12
     .                         + r22r22*dzr21r22)
     .             + dzreploc(20)*(r12r32*r11r32
     .                         + r32r12*r31r12
     .                         + r22r32*r21r32
     .                         + r32r22*r31r22)
                 dzrepij(8,6)= dzrepij(8,6)
     .             + reploc(20)*(dzr12r32*r11r32
     .                         + dzr32r12*r31r12
     .                         + dzr22r32*r21r32
     .                         + dzr32r22*r31r22)
     .             + reploc(20)*(r12r32*dzr11r32
     .                         + r32r12*dzr31r12
     .                         + r22r32*dzr21r32
     .                         + r32r22*dzr31r22)
     .             + dzreploc(21)*(r12r12*r21r22
     .                         + r22r22*r11r12)
     .             + reploc(21)*(dzr12r12*r21r22
     .                         + dzr22r22*r11r12)
     .             + reploc(21)*(r12r12*dzr21r22
     .                         + r22r22*dzr11r12)
     .             + dzreploc(22)*(r12r22*r11r22+ r22r12*r21r12)
     .             + reploc(22)*(dzr12r22*r11r22+ dzr22r12*r21r12)
     .             + reploc(22)*(r12r22*dzr11r22+ r22r12*dzr21r12)
!
!       (px py | py pz):
!
        repij(6,9) = reploc(16)*r31r32*r32r33
     .             + reploc(17)*radd02*r32r33
     .             + reploc(18)*r31r32*radd05
     .             + reploc(19)*(r11r12*r12r13
     .                         + r21r22*r22r23)
     .             + reploc(20)*(radd11*radd13
     .                         + radd12*radd14)
     .             + reploc(21)*(r11r12*r22r23
     .                         + r21r22*r12r13)
     .             + reploc(22)*radd17*(r12r23+r22r13)

        dxrepij(6,9) = dxreploc(16)*r31r32*r32r33
     .             + reploc(16)*dxr31r32*r32r33
     .             + reploc(16)*r31r32*dxr32r33
     .             + dxreploc(17)*radd02*r32r33
     .             + reploc(17)*dxradd02*r32r33
     .             + reploc(17)*radd02*dxr32r33
     .             + dxreploc(18)*r31r32*radd05
     .             + reploc(18)*dxr31r32*radd05
     .             + reploc(18)*r31r32*dxradd05
     .             + dxreploc(19)*(r11r12*r12r13
     .                         + r21r22*r22r23)
     .             + reploc(19)*(dxr11r12*r12r13
     .                         + dxr21r22*r22r23)
     .             + reploc(19)*(r11r12*dxr12r13
     .                         + r21r22*dxr22r23)
                   dxrepij(6,9)=dxrepij(6,9)
     .             + dxreploc(20)*(radd11*radd13
     .                         + radd12*radd14)
     .             + reploc(20)*(dxradd11*radd13
     .                         + dxradd12*radd14)
     .             + reploc(20)*(radd11*dxradd13
     .                         + radd12*dxradd14)
     .             + dxreploc(21)*(r11r12*r22r23
     .                         + r21r22*r12r13)
     .             + reploc(21)*(dxr11r12*r22r23
     .                         + dxr21r22*r12r13)
     .             + reploc(21)*(r11r12*dxr22r23
     .                         + r21r22*dxr12r13)
     .             + dxreploc(22)*radd17*(r12r23+r22r13)
     .             + reploc(22)*dxradd17*(r12r23+r22r13)
     .             + reploc(22)*radd17*(dxr12r23+dxr22r13)

        dyrepij(6,9) = dyreploc(16)*r31r32*r32r33
     .             + reploc(16)*dyr31r32*r32r33
     .             + reploc(16)*r31r32*dyr32r33
     .             + dyreploc(17)*radd02*r32r33
     .             + reploc(17)*dyradd02*r32r33
     .             + reploc(17)*radd02*dyr32r33
     .             + dyreploc(18)*r31r32*radd05
     .             + reploc(18)*dyr31r32*radd05
     .             + reploc(18)*r31r32*dyradd05
     .             + dyreploc(19)*(r11r12*r12r13
     .                         + r21r22*r22r23)
     .             + reploc(19)*(dyr11r12*r12r13
     .                         + dyr21r22*r22r23)
     .             + reploc(19)*(r11r12*dyr12r13
     .                         + r21r22*dyr22r23)
                  dyrepij(6,9)=dyrepij(6,9)
     .             + dyreploc(20)*(radd11*radd13
     .                         + radd12*radd14)
     .             + reploc(20)*(dyradd11*radd13
     .                         + dyradd12*radd14)
     .             + reploc(20)*(radd11*dyradd13
     .                         + radd12*dyradd14)
     .             + dyreploc(21)*(r11r12*r22r23
     .                         + r21r22*r12r13)
     .             + reploc(21)*(dyr11r12*r22r23
     .                         + dyr21r22*r12r13)
     .             + reploc(21)*(r11r12*dyr22r23
     .                         + r21r22*dyr12r13)
     .             + dyreploc(22)*radd17*(r12r23+r22r13)
     .             + reploc(22)*dyradd17*(r12r23+r22r13)
     .             + reploc(22)*radd17*(dyr12r23+dyr22r13)

        dzrepij(6,9) = dzreploc(16)*r31r32*r32r33
     .             + reploc(16)*dzr31r32*r32r33
     .             + reploc(16)*r31r32*dzr32r33
     .             + dzreploc(17)*radd02*r32r33
     .             + reploc(17)*dzradd02*r32r33
     .             + reploc(17)*radd02*dzr32r33
     .             + dzreploc(18)*r31r32*radd05
     .             + reploc(18)*dzr31r32*radd05
     .             + reploc(18)*r31r32*dzradd05
     .             + dzreploc(19)*(r11r12*r12r13
     .                         + r21r22*r22r23)
     .             + reploc(19)*(dzr11r12*r12r13
     .                         + dzr21r22*r22r23)
     .             + reploc(19)*(r11r12*dzr12r13
     .                         + r21r22*dzr22r23)
               dzrepij(6,9)=dzrepij(6,9)
     .             + dzreploc(20)*(radd11*radd13
     .                         + radd12*radd14)
     .             + reploc(20)*(dzradd11*radd13
     .                         + dzradd12*radd14)
     .             + reploc(20)*(radd11*dzradd13
     .                         + radd12*dzradd14)
     .             + dzreploc(21)*(r11r12*r22r23
     .                         + r21r22*r12r13)
     .             + reploc(21)*(dzr11r12*r22r23
     .                         + dzr21r22*r12r13)
     .             + reploc(21)*(r11r12*dzr22r23
     .                         + r21r22*dzr12r13)
     .             + dzreploc(22)*radd17*(r12r23+r22r13)
     .             + reploc(22)*dzradd17*(r12r23+r22r13)
     .             + reploc(22)*radd17*(dzr12r23+dzr22r13)
!
!       (py pz | px py):
!
        repij(9,6) = repij(6,9)
                 dxrepij(9,6) = dxrepij(6,9)
                 dyrepij(9,6) = dyrepij(6,9)
                 dzrepij(9,6) = dzrepij(6,9)

!
!       (px py | pz pz):
!
        repij(6,10) = reploc(16)*r31r32*r33r33
     .              + reploc(17)*radd02*r33r33
     .              + reploc(18)*radd06*r31r32
     .              + reploc(19)*(r11r12*r13r13
     .                          + r21r22*r23r23)
     .              + reploc(20)*(radd11*r13r33
     .                          + radd12*r23r33)
     .              + reploc(21)*(r11r12*r23r23
     .                          + r21r22*r13r13)
     .              + reploc(22)*radd17*r13r23

        dxrepij(6,10) = dxreploc(16)*r31r32*r33r33
     .              + reploc(16)*dxr31r32*r33r33
     .              + reploc(16)*r31r32*dxr33r33
     .              + dxreploc(17)*radd02*r33r33
     .              + reploc(17)*dxradd02*r33r33
     .              + reploc(17)*radd02*dxr33r33
     .              + dxreploc(18)*radd06*r31r32
     .              + reploc(18)*dxradd06*r31r32
     .              + reploc(18)*radd06*dxr31r32
     .              + dxreploc(19)*(r11r12*r13r13
     .                          + r21r22*r23r23)
     .              + reploc(19)*(dxr11r12*r13r13
     .                          + dxr21r22*r23r23)
     .              + reploc(19)*(r11r12*dxr13r13
     .                          + r21r22*dxr23r23)
     .              + dxreploc(20)*(radd11*r13r33
     .                          + radd12*r23r33)
     .              + reploc(20)*(dxradd11*r13r33
     .                          + dxradd12*r23r33)
     .              + reploc(20)*(radd11*dxr13r33
     .                          + radd12*dxr23r33)
     .              + dxreploc(21)*(r11r12*r23r23
     .                          + r21r22*r13r13)
     .              + reploc(21)*(dxr11r12*r23r23
     .                          + dxr21r22*r13r13)
     .              + reploc(21)*(r11r12*dxr23r23
     .                          + r21r22*dxr13r13)
     .              + dxreploc(22)*radd17*r13r23
     .              + reploc(22)*dxradd17*r13r23
     .              + reploc(22)*radd17*dxr13r23

        dyrepij(6,10) = dyreploc(16)*r31r32*r33r33
     .              + reploc(16)*dyr31r32*r33r33
     .              + reploc(16)*r31r32*dyr33r33
     .              + dyreploc(17)*radd02*r33r33
     .              + reploc(17)*dyradd02*r33r33
     .              + reploc(17)*radd02*dyr33r33
     .              + dyreploc(18)*radd06*r31r32
     .              + reploc(18)*dyradd06*r31r32
     .              + reploc(18)*radd06*dyr31r32
     .              + dyreploc(19)*(r11r12*r13r13
     .                          + r21r22*r23r23)
     .              + reploc(19)*(dyr11r12*r13r13
     .                          + dyr21r22*r23r23)
     .              + reploc(19)*(r11r12*dyr13r13
     .                          + r21r22*dyr23r23)
     .              + dyreploc(20)*(radd11*r13r33
     .                          + radd12*r23r33)
     .              + reploc(20)*(dyradd11*r13r33
     .                          + dyradd12*r23r33)
     .              + reploc(20)*(radd11*dyr13r33
     .                          + radd12*dyr23r33)
     .              + dyreploc(21)*(r11r12*r23r23
     .                          + r21r22*r13r13)
     .              + reploc(21)*(dyr11r12*r23r23
     .                          + dyr21r22*r13r13)
     .              + reploc(21)*(r11r12*dyr23r23
     .                          + r21r22*dyr13r13)
     .              + dyreploc(22)*radd17*r13r23
     .              + reploc(22)*dyradd17*r13r23
     .              + reploc(22)*radd17*dyr13r23

        dzrepij(6,10) = dzreploc(16)*r31r32*r33r33
     .              + reploc(16)*dzr31r32*r33r33
     .              + reploc(16)*r31r32*dzr33r33
     .              + dzreploc(17)*radd02*r33r33
     .              + reploc(17)*dzradd02*r33r33
     .              + reploc(17)*radd02*dzr33r33
     .              + dzreploc(18)*radd06*r31r32
     .              + reploc(18)*dzradd06*r31r32
     .              + reploc(18)*radd06*dzr31r32
     .              + dzreploc(19)*(r11r12*r13r13
     .                          + r21r22*r23r23)
     .              + reploc(19)*(dzr11r12*r13r13
     .                          + dzr21r22*r23r23)
     .              + reploc(19)*(r11r12*dzr13r13
     .                          + r21r22*dzr23r23)
     .              + dzreploc(20)*(radd11*r13r33
     .                          + radd12*r23r33)
     .              + reploc(20)*(dzradd11*r13r33
     .                          + dzradd12*r23r33)
     .              + reploc(20)*(radd11*dzr13r33
     .                          + radd12*dzr23r33)
     .              + dzreploc(21)*(r11r12*r23r23
     .                          + r21r22*r13r13)
     .              + reploc(21)*(dzr11r12*r23r23
     .                          + dzr21r22*r13r13)
     .              + reploc(21)*(r11r12*dzr23r23
     .                          + r21r22*dzr13r13)
     .              + dzreploc(22)*radd17*r13r23
     .              + reploc(22)*dzradd17*r13r23
     .              + reploc(22)*radd17*dzr13r23
!
!       (pz pz | px py):
!
        repij(10,6) = reploc(16)*r33r33*r31r32
     .              + reploc(17)*radd06*r31r32
     .              + reploc(18)*radd02*r33r33
     .              + reploc(19)*(r13r13*r11r12
     .                          + r23r23*r21r22)
     .              + reploc(20)*(r13r33*r11r32
     .                          + r33r13*r31r12
     .                          + r23r33*r21r32
     .                          + r33r23*r31r22)
     .              + reploc(21)*(r13r13*r21r22
     .                          + r23r23*r11r12)
     .              + reploc(22)*(r13r23*r11r22
     .                          + r23r13*r21r12)

        dxrepij(10,6) = dxreploc(16)*r33r33*r31r32
     .              + reploc(16)*dxr33r33*r31r32
     .              + reploc(16)*r33r33*dxr31r32
     .              + dxreploc(17)*radd06*r31r32
     .              + reploc(17)*dxradd06*r31r32
     .              + reploc(17)*radd06*dxr31r32
     .              + dxreploc(18)*radd02*r33r33
     .              + reploc(18)*dxradd02*r33r33
     .              + reploc(18)*radd02*dxr33r33
     .              + dxreploc(19)*(r13r13*r11r12
     .                          + r23r23*r21r22)
     .              + reploc(19)*(dxr13r13*r11r12
     .                          + dxr23r23*r21r22)
     .              + reploc(19)*(r13r13*dxr11r12
     .                          + r23r23*dxr21r22)
                 dxrepij(10,6)= dxrepij(10,6)
     .              + dxreploc(20)*(r13r33*r11r32
     .                          + r33r13*r31r12
     .                          + r23r33*r21r32
     .                          + r33r23*r31r22)
     .              + reploc(20)*(dxr13r33*r11r32
     .                          + dxr33r13*r31r12
     .                          + dxr23r33*r21r32
     .                          + dxr33r23*r31r22)
     .              + reploc(20)*(r13r33*dxr11r32
     .                          + r33r13*dxr31r12
     .                          + r23r33*dxr21r32
     .                          + r33r23*dxr31r22)
     .              + dxreploc(21)*(r13r13*r21r22
     .                          + r23r23*r11r12)
     .              + reploc(21)*(dxr13r13*r21r22
     .                          + dxr23r23*r11r12)
     .              + reploc(21)*(r13r13*dxr21r22
     .                          + r23r23*dxr11r12)
     .              + dxreploc(22)*(r13r23*r11r22
     .                          + r23r13*r21r12)
     .              + reploc(22)*(dxr13r23*r11r22
     .                          + dxr23r13*r21r12)
     .              + reploc(22)*(r13r23*dxr11r22
     .                          + r23r13*dxr21r12)

        dyrepij(10,6) = dyreploc(16)*r33r33*r31r32
     .              + reploc(16)*dyr33r33*r31r32
     .              + reploc(16)*r33r33*dyr31r32
     .              + dyreploc(17)*radd06*r31r32
     .              + reploc(17)*dyradd06*r31r32
     .              + reploc(17)*radd06*dyr31r32
     .              + dyreploc(18)*radd02*r33r33
     .              + reploc(18)*dyradd02*r33r33
     .              + reploc(18)*radd02*dyr33r33
     .              + dyreploc(19)*(r13r13*r11r12
     .                          + r23r23*r21r22)
     .              + reploc(19)*(dyr13r13*r11r12
     .                          + dyr23r23*r21r22)
     .              + reploc(19)*(r13r13*dyr11r12
     .                          + r23r23*dyr21r22)
           dyrepij(10,6)=dyrepij(10,6)
     .              + dyreploc(20)*(r13r33*r11r32
     .                          + r33r13*r31r12
     .                          + r23r33*r21r32
     .                          + r33r23*r31r22)
     .              + reploc(20)*(dyr13r33*r11r32
     .                          + dyr33r13*r31r12
     .                          + dyr23r33*r21r32
     .                          + dyr33r23*r31r22)
     .              + reploc(20)*(r13r33*dyr11r32
     .                          + r33r13*dyr31r12
     .                          + r23r33*dyr21r32
     .                          + r33r23*dyr31r22)
     .              + dyreploc(21)*(r13r13*r21r22
     .                          + r23r23*r11r12)
     .              + reploc(21)*(dyr13r13*r21r22
     .                          + dyr23r23*r11r12)
     .              + reploc(21)*(r13r13*dyr21r22
     .                          + r23r23*dyr11r12)
     .              + dyreploc(22)*(r13r23*r11r22
     .                          + r23r13*r21r12)
     .              + reploc(22)*(dyr13r23*r11r22
     .                          + dyr23r13*r21r12)
     .              + reploc(22)*(r13r23*dyr11r22
     .                          + r23r13*dyr21r12)

        dzrepij(10,6) = dzreploc(16)*r33r33*r31r32
     .              + reploc(16)*dzr33r33*r31r32
     .              + reploc(16)*r33r33*dzr31r32
     .              + dzreploc(17)*radd06*r31r32
     .              + reploc(17)*dzradd06*r31r32
     .              + reploc(17)*radd06*dzr31r32
     .              + dzreploc(18)*radd02*r33r33
     .              + reploc(18)*dzradd02*r33r33
     .              + reploc(18)*radd02*dzr33r33
     .              + dzreploc(19)*(r13r13*r11r12
     .                          + r23r23*r21r22)
     .              + reploc(19)*(dzr13r13*r11r12
     .                          + dzr23r23*r21r22)
     .              + reploc(19)*(r13r13*dzr11r12
     .                          + r23r23*dzr21r22)
              dzrepij(10,6)=dzrepij(10,6)
     .              + dzreploc(20)*(r13r33*r11r32
     .                          + r33r13*r31r12
     .                          + r23r33*r21r32
     .                          + r33r23*r31r22)
     .              + reploc(20)*(dzr13r33*r11r32
     .                          + dzr33r13*r31r12
     .                          + dzr23r33*r21r32
     .                          + dzr33r23*r31r22)
     .              + reploc(20)*(r13r33*dzr11r32
     .                          + r33r13*dzr31r12
     .                          + r23r33*dzr21r32
     .                          + r33r23*dzr31r22)
     .              + dzreploc(21)*(r13r13*r21r22
     .                          + r23r23*r11r12)
     .              + reploc(21)*(dzr13r13*r21r22
     .                          + dzr23r23*r11r12)
     .              + reploc(21)*(r13r13*dzr21r22
     .                          + r23r23*dzr11r12)
     .              + dzreploc(22)*(r13r23*r11r22
     .                          + r23r13*r21r12)
     .              + reploc(22)*(dzr13r23*r11r22
     .                          + dzr23r13*r21r12)
     .              + reploc(22)*(r13r23*dzr11r22
     .                          + r23r13*dzr21r12)

!
!       (px pz | px pz):
!
        repij(7,7) = reploc(16)*r31r33*r31r33
     .             + reploc(17)*radd03*r31r33
     .             + reploc(18)*radd03*r31r33
     .             + reploc(19)*(r11r13*r11r13
     .                         + r21r23*r21r23)
     .             + reploc(20)*(radd09**2
     .                         + radd10**2)
     .             + reploc(21)*r11r13*r21r23*2.0d0
     .             + reploc(22)*radd15**2

 

        dxrepij(7,7) = dxreploc(16)*r31r33*r31r33
     .             + reploc(16)*dxr31r33*r31r33
     .             + reploc(16)*r31r33*dxr31r33
     .             + dxreploc(17)*radd03*r31r33
     .             + reploc(17)*dxradd03*r31r33
     .             + reploc(17)*radd03*dxr31r33
     .             + dxreploc(18)*radd03*r31r33
     .             + reploc(18)*dxradd03*r31r33
     .             + reploc(18)*radd03*dxr31r33
     .             + dxreploc(19)*(r11r13*r11r13
     .                         + r21r23*r21r23)
     .             + reploc(19)*(dxr11r13*r11r13
     .                         + dxr21r23*r21r23)
     .             + reploc(19)*(r11r13*dxr11r13
     .                         + r21r23*dxr21r23)
     .             + dxreploc(20)*(radd09**2
     .                         + radd10**2)
     .             + reploc(20)*2.0d0*(radd09*dxradd09
     .                               + radd10*dxradd10)
     .             + dxreploc(21)*r11r13*r21r23*2.0d0
     .             + reploc(21)*dxr11r13*r21r23*2.0d0
     .             + reploc(21)*r11r13*dxr21r23*2.0d0
     .             + dxreploc(22)*radd15**2
     .             + reploc(22)*radd15*dxradd15*2.0d0

        dyrepij(7,7) = dyreploc(16)*r31r33*r31r33
     .             + reploc(16)*dyr31r33*r31r33
     .             + reploc(16)*r31r33*dyr31r33
     .             + dyreploc(17)*radd03*r31r33
     .             + reploc(17)*dyradd03*r31r33
     .             + reploc(17)*radd03*dyr31r33
     .             + dyreploc(18)*radd03*r31r33
     .             + reploc(18)*dyradd03*r31r33
     .             + reploc(18)*radd03*dyr31r33
     .             + dyreploc(19)*(r11r13*r11r13
     .                         + r21r23*r21r23)
     .             + reploc(19)*(dyr11r13*r11r13
     .                         + dyr21r23*r21r23)
     .             + reploc(19)*(r11r13*dyr11r13
     .                         + r21r23*dyr21r23)
     .             + dyreploc(20)*(radd09**2
     .                         + radd10**2)
     .             + reploc(20)*2.0d0*(radd09*dyradd09
     .                         + radd10*dyradd10)
     .             + dyreploc(21)*r11r13*r21r23*2.0d0
     .             + reploc(21)*dyr11r13*r21r23*2.0d0
     .             + reploc(21)*r11r13*dyr21r23*2.0d0
     .             + dyreploc(22)*radd15**2
     .             + reploc(22)*radd15*dyradd15*2.0d0

        dzrepij(7,7) = dzreploc(16)*r31r33*r31r33
     .             + reploc(16)*dzr31r33*r31r33
     .             + reploc(16)*r31r33*dzr31r33
     .             + dzreploc(17)*radd03*r31r33
     .             + reploc(17)*dzradd03*r31r33
     .             + reploc(17)*radd03*dzr31r33
     .             + dzreploc(18)*radd03*r31r33
     .             + reploc(18)*dzradd03*r31r33
     .             + reploc(18)*radd03*dzr31r33
     .             + dzreploc(19)*(r11r13*r11r13
     .                         + r21r23*r21r23)
     .             + reploc(19)*(dzr11r13*r11r13
     .                         + dzr21r23*r21r23)
     .             + reploc(19)*(r11r13*dzr11r13
     .                         + r21r23*dzr21r23)
     .             + dzreploc(20)*(radd09**2
     .                         + radd10**2)
     .             + reploc(20)*2.0d0*(radd09*dzradd09
     .                               + radd10*dzradd10)
     .             + dzreploc(21)*r11r13*r21r23*2.0d0
     .             + reploc(21)*dzr11r13*r21r23*2.0d0
     .             + reploc(21)*r11r13*dzr21r23*2.0d0
     .             + dzreploc(22)*radd15**2
     .             + reploc(22)*radd15*dzradd15*2.0d0
!
!       (px pz | py py):
!
        repij(7,8) = reploc(16)*r31r33*r32r32
     .             + reploc(17)*radd03*r32r32
     .             + reploc(18)*radd04*r31r33
     .             + reploc(19)*(r11r13*r12r12
     .                         + r21r23*r22r22)
     .             + reploc(20)*(radd09*r12r32
     .                         + radd10*r22r32)
     .             + reploc(21)*(r11r13*r22r22
     .                         + r21r23*r12r12)
     .             + reploc(22)*radd15*r12r22

        dxrepij(7,8) = dxreploc(16)*r31r33*r32r32
     .             + reploc(16)*dxr31r33*r32r32
     .             + reploc(16)*r31r33*dxr32r32
     .             + dxreploc(17)*radd03*r32r32
     .             + reploc(17)*dxradd03*r32r32
     .             + reploc(17)*radd03*dxr32r32
     .             + dxreploc(18)*radd04*r31r33
     .             + reploc(18)*dxradd04*r31r33
     .             + reploc(18)*radd04*dxr31r33
     .             + dxreploc(19)*(r11r13*r12r12
     .                         + r21r23*r22r22)
     .             + reploc(19)*(dxr11r13*r12r12
     .                         + dxr21r23*r22r22)
     .             + reploc(19)*(r11r13*dxr12r12
     .                         + r21r23*dxr22r22)
                dxrepij(7,8)= dxrepij(7,8)
     .             + dxreploc(20)*(radd09*r12r32
     .                         + radd10*r22r32)
     .             + reploc(20)*(dxradd09*r12r32
     .                         + dxradd10*r22r32)
     .             + reploc(20)*(radd09*dxr12r32
     .                         + radd10*dxr22r32)
     .             + dxreploc(21)*(r11r13*r22r22
     .                         + r21r23*r12r12)
     .             + reploc(21)*(dxr11r13*r22r22
     .                         + dxr21r23*r12r12)
     .             + reploc(21)*(r11r13*dxr22r22
     .                         + r21r23*dxr12r12)
     .             + dxreploc(22)*radd15*r12r22
     .             + reploc(22)*dxradd15*r12r22
     .             + reploc(22)*radd15*dxr12r22

        dyrepij(7,8) = dyreploc(16)*r31r33*r32r32
     .             + reploc(16)*dyr31r33*r32r32
     .             + reploc(16)*r31r33*dyr32r32
     .             + dyreploc(17)*radd03*r32r32
     .             + reploc(17)*dyradd03*r32r32
     .             + reploc(17)*radd03*dyr32r32
     .             + dyreploc(18)*radd04*r31r33
     .             + reploc(18)*dyradd04*r31r33
     .             + reploc(18)*radd04*dyr31r33
     .             + dyreploc(19)*(r11r13*r12r12
     .                         + r21r23*r22r22)
     .             + reploc(19)*(dyr11r13*r12r12
     .                         + dyr21r23*r22r22)
     .             + reploc(19)*(r11r13*dyr12r12
     .                         + r21r23*dyr22r22)
             dyrepij(7,8)=dyrepij(7,8)
     .             + dyreploc(20)*(radd09*r12r32
     .                         + radd10*r22r32)
     .             + reploc(20)*(dyradd09*r12r32
     .                         + dyradd10*r22r32)
     .             + reploc(20)*(radd09*dyr12r32
     .                         + radd10*dyr22r32)
     .             + dyreploc(21)*(r11r13*r22r22
     .                         + r21r23*r12r12)
     .             + reploc(21)*(dyr11r13*r22r22
     .                         + dyr21r23*r12r12)
     .             + reploc(21)*(r11r13*dyr22r22
     .                         + r21r23*dyr12r12)
     .             + dyreploc(22)*radd15*r12r22
     .             + reploc(22)*dyradd15*r12r22
     .             + reploc(22)*radd15*dyr12r22

        dzrepij(7,8) = dzreploc(16)*r31r33*r32r32
     .             + reploc(16)*dzr31r33*r32r32
     .             + reploc(16)*r31r33*dzr32r32
     .             + dzreploc(17)*radd03*r32r32
     .             + reploc(17)*dzradd03*r32r32
     .             + reploc(17)*radd03*dzr32r32
     .             + dzreploc(18)*radd04*r31r33
     .             + reploc(18)*dzradd04*r31r33
     .             + reploc(18)*radd04*dzr31r33
     .             + dzreploc(19)*(r11r13*r12r12
     .                         + r21r23*r22r22)
     .             + reploc(19)*(dzr11r13*r12r12
     .                         + dzr21r23*r22r22)
     .             + reploc(19)*(r11r13*dzr12r12
     .                         + r21r23*dzr22r22)
            dzrepij(7,8)=dzrepij(7,8)
     .             + dzreploc(20)*(radd09*r12r32
     .                         + radd10*r22r32)
     .             + reploc(20)*(dzradd09*r12r32
     .                         + dzradd10*r22r32)
     .             + reploc(20)*(radd09*dzr12r32
     .                         + radd10*dzr22r32)
     .             + dzreploc(21)*(r11r13*r22r22
     .                         + r21r23*r12r12)
     .             + reploc(21)*(dzr11r13*r22r22
     .                         + dzr21r23*r12r12)
     .             + reploc(21)*(r11r13*dzr22r22
     .                         + r21r23*dzr12r12)
     .             + dzreploc(22)*radd15*r12r22
     .             + reploc(22)*dzradd15*r12r22
     .             + reploc(22)*radd15*dzr12r22
!
!       (py py | px pz):
!
        repij(8,7) = reploc(16)*r32r32*r31r33
     .             + reploc(17)*radd04*r31r33
     .             + reploc(18)*radd03*r32r32
     .             + reploc(19)*(r12r12*r11r13
     .                         + r22r22*r21r23)
     .             + reploc(20)*(r12r32*r11r33
     .                         + r32r12*r31r13
     .                         + r22r32*r21r33
     .                         + r32r22*r31r23)
     .             + reploc(21)*(r12r12*r21r23
     .                         + r22r22*r11r13)
     .             + reploc(22)*(r12r22*r11r23
     .                         + r22r12*r21r13)

        dxrepij(8,7) = dxreploc(16)*r32r32*r31r33
     .             + reploc(16)*dxr32r32*r31r33
     .             + reploc(16)*r32r32*dxr31r33
     .             + dxreploc(17)*radd04*r31r33
     .             + reploc(17)*dxradd04*r31r33
     .             + reploc(17)*radd04*dxr31r33
     .             + dxreploc(18)*radd03*r32r32
     .             + reploc(18)*dxradd03*r32r32
     .             + reploc(18)*radd03*dxr32r32
     .             + dxreploc(19)*(r12r12*r11r13
     .                         + r22r22*r21r23)
     .             + reploc(19)*(dxr12r12*r11r13
     .                         + dxr22r22*r21r23)
     .             + reploc(19)*(r12r12*dxr11r13
     .                         + r22r22*dxr21r23)
     .             + dxreploc(20)*(r12r32*r11r33
     .                         + r32r12*r31r13
     .                         + r22r32*r21r33
     .                         + r32r22*r31r23)
               dxrepij(8,7)=dxrepij(8,7)
     .             + reploc(20)*(dxr12r32*r11r33
     .                         + dxr32r12*r31r13
     .                         + dxr22r32*r21r33
     .                         + dxr32r22*r31r23)
     .             + reploc(20)*(r12r32*dxr11r33
     .                         + r32r12*dxr31r13
     .                         + r22r32*dxr21r33
     .                         + r32r22*dxr31r23)
     .             + dxreploc(21)*(r12r12*r21r23
     .                         + r22r22*r11r13)
     .             + reploc(21)*(dxr12r12*r21r23
     .                         + dxr22r22*r11r13)
     .             + reploc(21)*(r12r12*dxr21r23
     .                         + r22r22*dxr11r13)
     .             + dxreploc(22)*(r12r22*r11r23
     .                         + r22r12*r21r13)
     .             + reploc(22)*(dxr12r22*r11r23
     .                         + dxr22r12*r21r13)
     .             + reploc(22)*(r12r22*dxr11r23
     .                         + r22r12*dxr21r13)

        dyrepij(8,7) = dyreploc(16)*r32r32*r31r33
     .             + reploc(16)*dyr32r32*r31r33
     .             + reploc(16)*r32r32*dyr31r33
     .             + dyreploc(17)*radd04*r31r33
     .             + reploc(17)*dyradd04*r31r33
     .             + reploc(17)*radd04*dyr31r33
     .             + dyreploc(18)*radd03*r32r32
     .             + reploc(18)*dyradd03*r32r32
     .             + reploc(18)*radd03*dyr32r32
     .             + dyreploc(19)*(r12r12*r11r13
     .                         + r22r22*r21r23)
     .             + reploc(19)*(dyr12r12*r11r13
     .                         + dyr22r22*r21r23)
     .             + reploc(19)*(r12r12*dyr11r13
     .                         + r22r22*dyr21r23)
     .             + dyreploc(20)*(r12r32*r11r33
     .                         + r32r12*r31r13
     .                         + r22r32*r21r33
     .                         + r32r22*r31r23)
               dyrepij(8,7)=dyrepij(8,7)
     .             + reploc(20)*(dyr12r32*r11r33
     .                         + dyr32r12*r31r13
     .                         + dyr22r32*r21r33
     .                         + dyr32r22*r31r23)
     .             + reploc(20)*(r12r32*dyr11r33
     .                         + r32r12*dyr31r13
     .                         + r22r32*dyr21r33
     .                         + r32r22*dyr31r23)
     .             + dyreploc(21)*(r12r12*r21r23
     .                         + r22r22*r11r13)
     .             + reploc(21)*(dyr12r12*r21r23
     .                         + dyr22r22*r11r13)
     .             + reploc(21)*(r12r12*dyr21r23
     .                         + r22r22*dyr11r13)
     .             + dyreploc(22)*(r12r22*r11r23
     .                         + r22r12*r21r13)
     .             + reploc(22)*(dyr12r22*r11r23
     .                         + dyr22r12*r21r13)
     .             + reploc(22)*(r12r22*dyr11r23
     .                         + r22r12*dyr21r13)

        dzrepij(8,7) = dzreploc(16)*r32r32*r31r33
     .             + reploc(16)*dzr32r32*r31r33
     .             + reploc(16)*r32r32*dzr31r33
     .             + dzreploc(17)*radd04*r31r33
     .             + reploc(17)*dzradd04*r31r33
     .             + reploc(17)*radd04*dzr31r33
     .             + dzreploc(18)*radd03*r32r32
     .             + reploc(18)*dzradd03*r32r32
     .             + reploc(18)*radd03*dzr32r32
     .             + dzreploc(19)*(r12r12*r11r13
     .                         + r22r22*r21r23)
     .             + reploc(19)*(dzr12r12*r11r13
     .                         + dzr22r22*r21r23)
     .             + reploc(19)*(r12r12*dzr11r13
     .                         + r22r22*dzr21r23)
     .             + dzreploc(20)*(r12r32*r11r33
     .                         + r32r12*r31r13
     .                         + r22r32*r21r33
     .                         + r32r22*r31r23)
                 dzrepij(8,7)= dzrepij(8,7)
     .             + reploc(20)*(dzr12r32*r11r33
     .                         + dzr32r12*r31r13
     .                         + dzr22r32*r21r33
     .                         + dzr32r22*r31r23)
     .             + reploc(20)*(r12r32*dzr11r33
     .                         + r32r12*dzr31r13
     .                         + r22r32*dzr21r33
     .                         + r32r22*dzr31r23)
     .             + dzreploc(21)*(r12r12*r21r23
     .                         + r22r22*r11r13)
     .             + reploc(21)*(dzr12r12*r21r23
     .                         + dzr22r22*r11r13)
     .             + reploc(21)*(r12r12*dzr21r23
     .                         + r22r22*dzr11r13)
     .             + dzreploc(22)*(r12r22*r11r23
     .                         + r22r12*r21r13)
     .             + reploc(22)*(dzr12r22*r11r23
     .                         + dzr22r12*r21r13)
     .             + reploc(22)*(r12r22*dzr11r23
     .                         + r22r12*dzr21r13)
!
!       (px pz | py pz):
!
        repij(7,9) = reploc(16)*r31r33*r32r33
     .             + reploc(17)*radd03*r32r33
     .             + reploc(18)*r31r33*radd05
     .             + reploc(19)*(r11r13*r12r13
     .                         + r21r23*r22r23)
     .             + reploc(20)*(radd09*radd13
     .                         + radd10*r22r33
     .                         + radd10*r32r23)
     .             + reploc(21)*(r11r13*r22r23
     .                         + r21r23*r12r13)
     .             + reploc(22)*radd15*radd16

        dxrepij(7,9) = dxreploc(16)*r31r33*r32r33
     .             + reploc(16)*dxr31r33*r32r33
     .             + reploc(16)*r31r33*dxr32r33
     .             + dxreploc(17)*radd03*r32r33
     .             + reploc(17)*dxradd03*r32r33
     .             + reploc(17)*radd03*dxr32r33
     .             + dxreploc(18)*r31r33*radd05
     .             + reploc(18)*dxr31r33*radd05
     .             + reploc(18)*r31r33*dxradd05
     .             + dxreploc(19)*(r11r13*r12r13
     .                         + r21r23*r22r23)
     .             + reploc(19)*(dxr11r13*r12r13
     .                         + dxr21r23*r22r23)
     .             + reploc(19)*(r11r13*dxr12r13
     .                         + r21r23*dxr22r23)
     .             + dxreploc(20)*(radd09*radd13
     .                         + radd10*r22r33
     .                         + radd10*r32r23)
              dxrepij(7,9)=dxrepij(7,9)
     .             + reploc(20)*(dxradd09*radd13
     .                         + dxradd10*r22r33
     .                         + dxradd10*r32r23)
     .             + reploc(20)*(radd09*dxradd13
     .                         + radd10*dxr22r33
     .                         + radd10*dxr32r23)
     .             + dxreploc(21)*(r11r13*r22r23
     .                         + r21r23*r12r13)
     .             + reploc(21)*(dxr11r13*r22r23
     .                         + dxr21r23*r12r13)
     .             + reploc(21)*(r11r13*dxr22r23
     .                         + r21r23*dxr12r13)
     .             + dxreploc(22)*radd15*radd16
     .             + reploc(22)*dxradd15*radd16
     .             + reploc(22)*radd15*dxradd16

        dyrepij(7,9) = dyreploc(16)*r31r33*r32r33
     .             + reploc(16)*dyr31r33*r32r33
     .             + reploc(16)*r31r33*dyr32r33
     .             + dyreploc(17)*radd03*r32r33
     .             + reploc(17)*dyradd03*r32r33
     .             + reploc(17)*radd03*dyr32r33
     .             + dyreploc(18)*r31r33*radd05
     .             + reploc(18)*dyr31r33*radd05
     .             + reploc(18)*r31r33*dyradd05
     .             + dyreploc(19)*(r11r13*r12r13
     .                         + r21r23*r22r23)
     .             + reploc(19)*(dyr11r13*r12r13
     .                         + dyr21r23*r22r23)
     .             + reploc(19)*(r11r13*dyr12r13
     .                         + r21r23*dyr22r23)
     .             + dyreploc(20)*(radd09*radd13
     .                         + radd10*r22r33
     .                         + radd10*r32r23)
               dyrepij(7,9)=dyrepij(7,9)
     .             + reploc(20)*(dyradd09*radd13
     .                         + dyradd10*r22r33
     .                         + dyradd10*r32r23)
     .             + reploc(20)*(radd09*dyradd13
     .                         + radd10*dyr22r33
     .                         + radd10*dyr32r23)
     .             + dyreploc(21)*(r11r13*r22r23
     .                         + r21r23*r12r13)
     .             + reploc(21)*(dyr11r13*r22r23
     .                         + dyr21r23*r12r13)
     .             + reploc(21)*(r11r13*dyr22r23
     .                         + r21r23*dyr12r13)
     .             + dyreploc(22)*radd15*radd16
     .             + reploc(22)*dyradd15*radd16
     .             + reploc(22)*radd15*dyradd16

        dzrepij(7,9) = dzreploc(16)*r31r33*r32r33
     .             + reploc(16)*dzr31r33*r32r33
     .             + reploc(16)*r31r33*dzr32r33
     .             + dzreploc(17)*radd03*r32r33
     .             + reploc(17)*dzradd03*r32r33
     .             + reploc(17)*radd03*dzr32r33
     .             + dzreploc(18)*r31r33*radd05
     .             + reploc(18)*dzr31r33*radd05
     .             + reploc(18)*r31r33*dzradd05
     .             + dzreploc(19)*(r11r13*r12r13
     .                         + r21r23*r22r23)
     .             + reploc(19)*(dzr11r13*r12r13
     .                         + dzr21r23*r22r23)
     .             + reploc(19)*(r11r13*dzr12r13
     .                         + r21r23*dzr22r23)
     .             + dzreploc(20)*(radd09*radd13
     .                         + radd10*r22r33
     .                         + radd10*r32r23)
     .             + reploc(20)*(dzradd09*radd13
     .                         + dzradd10*r22r33
     .                         + dzradd10*r32r23)
     .             + reploc(20)*(radd09*dzradd13
     .                         + radd10*dzr22r33
     .                         + radd10*dzr32r23)
     .             + dzreploc(21)*(r11r13*r22r23
     .                         + r21r23*r12r13)
     .             + reploc(21)*(dzr11r13*r22r23
     .                         + dzr21r23*r12r13)
     .             + reploc(21)*(r11r13*dzr22r23
     .                         + r21r23*dzr12r13)
     .             + dzreploc(22)*radd15*radd16
     .             + reploc(22)*dzradd15*radd16
     .             + reploc(22)*radd15*dzradd16
!
!       (py pz | px pz):
!
        repij(9,7) = repij(7,9)
         dxrepij(9,7) = dxrepij(7,9)
         dyrepij(9,7) = dyrepij(7,9)
         dzrepij(9,7) = dzrepij(7,9)
!
!       (px pz | pz pz):
!
        repij(7,10) = reploc(16)*r31r33*r33r33
     .              + reploc(17)*radd03*r33r33
     .              + reploc(18)*radd06*r31r33
     .              + reploc(19)*(r11r13*r13r13
     .                          + r21r23*r23r23)
     .              + reploc(20)*(radd09*r13r33
     .                          + radd10*r23r33)
     .              + reploc(21)*(r11r13*r23r23
     .                          + r21r23*r13r13)
     .              + reploc(22)*radd15*r13r23

        dxrepij(7,10) = dxreploc(16)*r31r33*r33r33
     .              + reploc(16)*dxr31r33*r33r33
     .              + reploc(16)*r31r33*dxr33r33
     .              + dxreploc(17)*radd03*r33r33
     .              + reploc(17)*dxradd03*r33r33
     .              + reploc(17)*radd03*dxr33r33
     .              + dxreploc(18)*radd06*r31r33
     .              + reploc(18)*dxradd06*r31r33
     .              + reploc(18)*radd06*dxr31r33
     .              + dxreploc(19)*(r11r13*r13r13
     .                          + r21r23*r23r23)
     .              + reploc(19)*(dxr11r13*r13r13
     .                          + dxr21r23*r23r23)
     .              + reploc(19)*(r11r13*dxr13r13
     .                          + r21r23*dxr23r23)
               dxrepij(7,10)=dxrepij(7,10)
     .              + dxreploc(20)*(radd09*r13r33
     .                          + radd10*r23r33)
     .              + reploc(20)*(dxradd09*r13r33
     .                          + dxradd10*r23r33)
     .              + reploc(20)*(radd09*dxr13r33
     .                          + radd10*dxr23r33)
     .              + dxreploc(21)*(r11r13*r23r23
     .                          + r21r23*r13r13)
     .              + reploc(21)*(dxr11r13*r23r23
     .                          + dxr21r23*r13r13)
     .              + reploc(21)*(r11r13*dxr23r23
     .                          + r21r23*dxr13r13)
     .              + dxreploc(22)*radd15*r13r23
     .              + reploc(22)*dxradd15*r13r23
     .              + reploc(22)*radd15*dxr13r23

        dyrepij(7,10) = dyreploc(16)*r31r33*r33r33
     .              + reploc(16)*dyr31r33*r33r33
     .              + reploc(16)*r31r33*dyr33r33
     .              + dyreploc(17)*radd03*r33r33
     .              + reploc(17)*dyradd03*r33r33
     .              + reploc(17)*radd03*dyr33r33
     .              + dyreploc(18)*radd06*r31r33
     .              + reploc(18)*dyradd06*r31r33
     .              + reploc(18)*radd06*dyr31r33
     .              + dyreploc(19)*(r11r13*r13r13
     .                          + r21r23*r23r23)
     .              + reploc(19)*(dyr11r13*r13r13
     .                          + dyr21r23*r23r23)
     .              + reploc(19)*(r11r13*dyr13r13
     .                          + r21r23*dyr23r23)
               dyrepij(7,10)=dyrepij(7,10)
     .              + dyreploc(20)*(radd09*r13r33
     .                          + radd10*r23r33)
     .              + reploc(20)*(dyradd09*r13r33
     .                          + dyradd10*r23r33)
     .              + reploc(20)*(radd09*dyr13r33
     .                          + radd10*dyr23r33)
     .              + dyreploc(21)*(r11r13*r23r23
     .                          + r21r23*r13r13)
     .              + reploc(21)*(dyr11r13*r23r23
     .                          + dyr21r23*r13r13)
     .              + reploc(21)*(r11r13*dyr23r23
     .                          + r21r23*dyr13r13)
     .              + dyreploc(22)*radd15*r13r23
     .              + reploc(22)*dyradd15*r13r23
     .              + reploc(22)*radd15*dyr13r23

        dzrepij(7,10) = dzreploc(16)*r31r33*r33r33
     .              + reploc(16)*dzr31r33*r33r33
     .              + reploc(16)*r31r33*dzr33r33
     .              + dzreploc(17)*radd03*r33r33
     .              + reploc(17)*dzradd03*r33r33
     .              + reploc(17)*radd03*dzr33r33
     .              + dzreploc(18)*radd06*r31r33
     .              + reploc(18)*dzradd06*r31r33
     .              + reploc(18)*radd06*dzr31r33
     .              + dzreploc(19)*(r11r13*r13r13
     .                          + r21r23*r23r23)
     .              + reploc(19)*(dzr11r13*r13r13
     .                          + dzr21r23*r23r23)
     .              + reploc(19)*(r11r13*dzr13r13
     .                          + r21r23*dzr23r23)
                dzrepij(7,10)=dzrepij(7,10)
     .              + dzreploc(20)*(radd09*r13r33
     .                          + radd10*r23r33)
     .              + reploc(20)*(dzradd09*r13r33
     .                          + dzradd10*r23r33)
     .              + reploc(20)*(radd09*dzr13r33
     .                          + radd10*dzr23r33)
     .              + dzreploc(21)*(r11r13*r23r23
     .                          + r21r23*r13r13)
     .              + reploc(21)*(dzr11r13*r23r23
     .                          + dzr21r23*r13r13)
     .              + reploc(21)*(r11r13*dzr23r23
     .                          + r21r23*dzr13r13)
     .              + dzreploc(22)*radd15*r13r23
     .              + reploc(22)*dzradd15*r13r23
     .              + reploc(22)*radd15*dzr13r23
!
!       (pz pz | px pz):
!
        repij(10,7) = reploc(16)*r33r33*r31r33
     .              + reploc(17)*radd06*r31r33
     .              + reploc(18)*radd03*r33r33
     .              + reploc(19)*(r13r13*r11r13
     .                          + r23r23*r21r23)
     .              + reploc(20)*(r13r33*r11r33
     .                          + r33r13*r31r13
     .                          + r23r33*r21r33
     .                          + r33r23*r31r23)
     .              + reploc(21)*(r13r13*r21r23
     .                          + r23r23*r11r13)
     .              + reploc(22)*(r13r23*r11r23
     .                          + r23r13*r21r13)

        dxrepij(10,7) = dxreploc(16)*r33r33*r31r33
     .              + reploc(16)*dxr33r33*r31r33
     .              + reploc(16)*r33r33*dxr31r33
     .              + dxreploc(17)*radd06*r31r33
     .              + reploc(17)*dxradd06*r31r33
     .              + reploc(17)*radd06*dxr31r33
     .              + dxreploc(18)*radd03*r33r33
     .              + reploc(18)*dxradd03*r33r33
     .              + reploc(18)*radd03*dxr33r33
     .              + dxreploc(19)*(r13r13*r11r13
     .                          + r23r23*r21r23)
     .              + reploc(19)*(dxr13r13*r11r13
     .                          + dxr23r23*r21r23)
     .              + reploc(19)*(r13r13*dxr11r13
     .                          + r23r23*dxr21r23)
     .              + dxreploc(20)*(r13r33*r11r33
     .                          + r33r13*r31r13
     .                          + r23r33*r21r33
     .                          + r33r23*r31r23)
                dxrepij(10,7)=dxrepij(10,7)
     .              + reploc(20)*(dxr13r33*r11r33
     .                          + dxr33r13*r31r13
     .                          + dxr23r33*r21r33
     .                          + dxr33r23*r31r23)
     .              + reploc(20)*(r13r33*dxr11r33
     .                          + r33r13*dxr31r13
     .                          + r23r33*dxr21r33
     .                          + r33r23*dxr31r23)
     .              + dxreploc(21)*(r13r13*r21r23
     .                          + r23r23*r11r13)
     .              + reploc(21)*(dxr13r13*r21r23
     .                          + dxr23r23*r11r13)
     .              + reploc(21)*(r13r13*dxr21r23
     .                          + r23r23*dxr11r13)
     .              + dxreploc(22)*(r13r23*r11r23
     .                          + r23r13*r21r13)
     .              + reploc(22)*(dxr13r23*r11r23
     .                          + dxr23r13*r21r13)
     .              + reploc(22)*(r13r23*dxr11r23
     .                          + r23r13*dxr21r13)

        dyrepij(10,7) = dyreploc(16)*r33r33*r31r33
     .              + reploc(16)*dyr33r33*r31r33
     .              + reploc(16)*r33r33*dyr31r33
     .              + dyreploc(17)*radd06*r31r33
     .              + reploc(17)*dyradd06*r31r33
     .              + reploc(17)*radd06*dyr31r33
     .              + dyreploc(18)*radd03*r33r33
     .              + reploc(18)*dyradd03*r33r33
     .              + reploc(18)*radd03*dyr33r33
     .              + dyreploc(19)*(r13r13*r11r13
     .                          + r23r23*r21r23)
     .              + reploc(19)*(dyr13r13*r11r13
     .                          + dyr23r23*r21r23)
     .              + reploc(19)*(r13r13*dyr11r13
     .                          + r23r23*dyr21r23)
     .              + dyreploc(20)*(r13r33*r11r33
     .                          + r33r13*r31r13
     .                          + r23r33*r21r33
     .                          + r33r23*r31r23)
     .              + reploc(20)*(dyr13r33*r11r33
     .                          + dyr33r13*r31r13
     .                          + dyr23r33*r21r33
     .                          + dyr33r23*r31r23)
                 dyrepij(10,7)=dyrepij(10,7)
     .              + reploc(20)*(r13r33*dyr11r33
     .                          + r33r13*dyr31r13
     .                          + r23r33*dyr21r33
     .                          + r33r23*dyr31r23)
     .              + dyreploc(21)*(r13r13*r21r23
     .                          + r23r23*r11r13)
     .              + reploc(21)*(dyr13r13*r21r23
     .                          + dyr23r23*r11r13)
     .              + reploc(21)*(r13r13*dyr21r23
     .                          + r23r23*dyr11r13)
     .              + dyreploc(22)*(r13r23*r11r23
     .                          + r23r13*r21r13)
     .              + reploc(22)*(dyr13r23*r11r23
     .                          + dyr23r13*r21r13)
     .              + reploc(22)*(r13r23*dyr11r23
     .                          + r23r13*dyr21r13)

        dzrepij(10,7) = dzreploc(16)*r33r33*r31r33
     .              + reploc(16)*dzr33r33*r31r33
     .              + reploc(16)*r33r33*dzr31r33
     .              + dzreploc(17)*radd06*r31r33
     .              + reploc(17)*dzradd06*r31r33
     .              + reploc(17)*radd06*dzr31r33
     .              + dzreploc(18)*radd03*r33r33
     .              + reploc(18)*dzradd03*r33r33
     .              + reploc(18)*radd03*dzr33r33
     .              + dzreploc(19)*(r13r13*r11r13
     .                          + r23r23*r21r23)
     .              + reploc(19)*(dzr13r13*r11r13
     .                          + dzr23r23*r21r23)
     .              + reploc(19)*(r13r13*dzr11r13
     .                          + r23r23*dzr21r23)
     .              + dzreploc(20)*(r13r33*r11r33
     .                          + r33r13*r31r13
     .                          + r23r33*r21r33
     .                          + r33r23*r31r23)
                 dzrepij(10,7)= dzrepij(10,7)
     .              + reploc(20)*(dzr13r33*r11r33
     .                          + dzr33r13*r31r13
     .                          + dzr23r33*r21r33
     .                          + dzr33r23*r31r23)
     .              + reploc(20)*(r13r33*dzr11r33
     .                          + r33r13*dzr31r13
     .                          + r23r33*dzr21r33
     .                          + r33r23*dzr31r23)
     .              + dzreploc(21)*(r13r13*r21r23
     .                          + r23r23*r11r13)
     .              + reploc(21)*(dzr13r13*r21r23
     .                          + dzr23r23*r11r13)
     .              + reploc(21)*(r13r13*dzr21r23
     .                          + r23r23*dzr11r13)
     .              + dzreploc(22)*(r13r23*r11r23
     .                          + r23r13*r21r13)
     .              + reploc(22)*(dzr13r23*r11r23
     .                          + dzr23r13*r21r13)
     .              + reploc(22)*(r13r23*dzr11r23
     .                          + r23r13*dzr21r13)

!
!       (py py | py py):
!
        repij(8,8) = reploc(16)*r32r32*r32r32
     .             + reploc(17)*radd04*r32r32
     .             + reploc(18)*radd04*r32r32
     .             + reploc(19)*(r12r12*r12r12
     .                         + r22r22*r22r22)
     .             + reploc(20)*(r12r32*r12r32
     .                         + r22r32*r22r32)
     .             + reploc(21)*r12r12*r22r22*2.0d0
     .             + reploc(22)*r12r22*r12r22

        dxrepij(8,8) = dxreploc(16)*r32r32*r32r32
     .             + reploc(16)*dxr32r32*r32r32
     .             + reploc(16)*r32r32*dxr32r32
     .             + dxreploc(17)*radd04*r32r32
     .             + reploc(17)*dxradd04*r32r32
     .             + reploc(17)*radd04*dxr32r32
     .             + dxreploc(18)*radd04*r32r32
     .             + reploc(18)*dxradd04*r32r32
     .             + reploc(18)*radd04*dxr32r32
     .             + dxreploc(19)*(r12r12*r12r12
     .                         + r22r22*r22r22)
     .             + reploc(19)*(dxr12r12*r12r12
     .                         + dxr22r22*r22r22)
     .             + reploc(19)*(r12r12*dxr12r12
     .                         + r22r22*dxr22r22)
              dxrepij(8,8)=dxrepij(8,8)
     .             + dxreploc(20)*(r12r32*r12r32
     .                         + r22r32*r22r32)
     .             + reploc(20)*(dxr12r32*r12r32
     .                         + dxr22r32*r22r32)
     .             + reploc(20)*(r12r32*dxr12r32
     .                         + r22r32*dxr22r32)
     .             + dxreploc(21)*r12r12*r22r22*2.0d0
     .             + reploc(21)*dxr12r12*r22r22*2.0d0
     .             + reploc(21)*r12r12*dxr22r22*2.0d0
     .             + dxreploc(22)*r12r22*r12r22
     .             + reploc(22)*dxr12r22*r12r22
     .             + reploc(22)*r12r22*dxr12r22

        dyrepij(8,8) = dyreploc(16)*r32r32*r32r32
     .             + reploc(16)*dyr32r32*r32r32
     .             + reploc(16)*r32r32*dyr32r32
     .             + dyreploc(17)*radd04*r32r32
     .             + reploc(17)*dyradd04*r32r32
     .             + reploc(17)*radd04*dyr32r32
     .             + dyreploc(18)*radd04*r32r32
     .             + reploc(18)*dyradd04*r32r32
     .             + reploc(18)*radd04*dyr32r32
     .             + dyreploc(19)*(r12r12*r12r12
     .                         + r22r22*r22r22)
     .             + reploc(19)*(dyr12r12*r12r12
     .                         + dyr22r22*r22r22)
     .             + reploc(19)*(r12r12*dyr12r12
     .                         + r22r22*dyr22r22)
                  dyrepij(8,8)=dyrepij(8,8)
     .             + dyreploc(20)*(r12r32*r12r32
     .                         + r22r32*r22r32)
     .             + reploc(20)*(dyr12r32*r12r32
     .                         + dyr22r32*r22r32)
     .             + reploc(20)*(r12r32*dyr12r32
     .                         + r22r32*dyr22r32)
     .             + dyreploc(21)*r12r12*r22r22*2.0d0
     .             + reploc(21)*dyr12r12*r22r22*2.0d0
     .             + reploc(21)*r12r12*dyr22r22*2.0d0
     .             + dyreploc(22)*r12r22*r12r22
     .             + reploc(22)*dyr12r22*r12r22
     .             + reploc(22)*r12r22*dyr12r22

        dzrepij(8,8) = dzreploc(16)*r32r32*r32r32
     .             + reploc(16)*dzr32r32*r32r32
     .             + reploc(16)*r32r32*dzr32r32
     .             + dzreploc(17)*radd04*r32r32
     .             + reploc(17)*dzradd04*r32r32
     .             + reploc(17)*radd04*dzr32r32
     .             + dzreploc(18)*radd04*r32r32
     .             + reploc(18)*dzradd04*r32r32
     .             + reploc(18)*radd04*dzr32r32
     .             + dzreploc(19)*(r12r12*r12r12
     .                         + r22r22*r22r22)
     .             + reploc(19)*(dzr12r12*r12r12
     .                         + dzr22r22*r22r22)
     .             + reploc(19)*(r12r12*dzr12r12
     .                         + r22r22*dzr22r22)
                dzrepij(8,8)=dzrepij(8,8)
     .             + dzreploc(20)*(r12r32*r12r32
     .                         + r22r32*r22r32)
     .             + reploc(20)*(dzr12r32*r12r32
     .                         + dzr22r32*r22r32)
     .             + reploc(20)*(r12r32*dzr12r32
     .                         + r22r32*dzr22r32)
     .             + dzreploc(21)*r12r12*r22r22*2.0d0
     .             + reploc(21)*dzr12r12*r22r22*2.0d0
     .             + reploc(21)*r12r12*dzr22r22*2.0d0
     .             + dzreploc(22)*r12r22*r12r22
     .             + reploc(22)*dzr12r22*r12r22
     .             + reploc(22)*r12r22*dzr12r22

!
!       (py py | py pz):
!
        repij(8,9) = reploc(16)*r32r32*r32r33
     .             + reploc(17)*radd04*r32r33
     .             + reploc(18)*radd05*r32r32
     .             + reploc(19)*(r12r12*r12r13
     .                         + r22r22*r22r23)
     .             + reploc(20)*(r12r32*r12r33
     .                         + r32r12*r32r13
     .                         + r22r32*r22r33
     .                         + r32r22*r32r23)
     .             + reploc(21)*(r12r12*r22r23
     .                         + r22r22*r12r13)
     .             + reploc(22)*(r12r22*r12r23
     .                         + r22r12*r22r13)

        dxrepij(8,9) = dxreploc(16)*r32r32*r32r33
     .             + reploc(16)*dxr32r32*r32r33
     .             + reploc(16)*r32r32*dxr32r33
     .             + dxreploc(17)*radd04*r32r33
     .             + reploc(17)*dxradd04*r32r33
     .             + reploc(17)*radd04*dxr32r33
     .             + dxreploc(18)*radd05*r32r32
     .             + reploc(18)*dxradd05*r32r32
     .             + reploc(18)*radd05*dxr32r32
     .             + dxreploc(19)*(r12r12*r12r13
     .                         + r22r22*r22r23)
     .             + reploc(19)*(dxr12r12*r12r13
     .                         + dxr22r22*r22r23)
     .             + reploc(19)*(r12r12*dxr12r13
     .                         + r22r22*dxr22r23)
     .             + dxreploc(20)*(r12r32*r12r33
     .                         + r32r12*r32r13
     .                         + r22r32*r22r33
     .                         + r32r22*r32r23)
                dxrepij(8,9)=dxrepij(8,9)
     .             + reploc(20)*(dxr12r32*r12r33
     .                         + dxr32r12*r32r13
     .                         + dxr22r32*r22r33
     .                         + dxr32r22*r32r23)
     .             + reploc(20)*(r12r32*dxr12r33
     .                         + r32r12*dxr32r13
     .                         + r22r32*dxr22r33
     .                         + r32r22*dxr32r23)
     .             + dxreploc(21)*(r12r12*r22r23
     .                         + r22r22*r12r13)
     .             + reploc(21)*(dxr12r12*r22r23
     .                         + dxr22r22*r12r13)
     .             + reploc(21)*(r12r12*dxr22r23
     .                         + r22r22*dxr12r13)
     .             + dxreploc(22)*(r12r22*r12r23
     .                         + r22r12*r22r13)
     .             + reploc(22)*(dxr12r22*r12r23
     .                         + dxr22r12*r22r13)
     .             + reploc(22)*(r12r22*dxr12r23
     .                         + r22r12*dxr22r13)

        dyrepij(8,9) = dyreploc(16)*r32r32*r32r33
     .             + reploc(16)*dyr32r32*r32r33
     .             + reploc(16)*r32r32*dyr32r33
     .             + dyreploc(17)*radd04*r32r33
     .             + reploc(17)*dyradd04*r32r33
     .             + reploc(17)*radd04*dyr32r33
     .             + dyreploc(18)*radd05*r32r32
     .             + reploc(18)*dyradd05*r32r32
     .             + reploc(18)*radd05*dyr32r32
     .             + dyreploc(19)*(r12r12*r12r13
     .                         + r22r22*r22r23)
     .             + reploc(19)*(dyr12r12*r12r13
     .                         + dyr22r22*r22r23)
     .             + reploc(19)*(r12r12*dyr12r13
     .                         + r22r22*dyr22r23)
     .             + dyreploc(20)*(r12r32*r12r33
     .                         + r32r12*r32r13
     .                         + r22r32*r22r33
     .                         + r32r22*r32r23)
                dyrepij(8,9)=dyrepij(8,9)
     .             + reploc(20)*(dyr12r32*r12r33
     .                         + dyr32r12*r32r13
     .                         + dyr22r32*r22r33
     .                         + dyr32r22*r32r23)
     .             + reploc(20)*(r12r32*dyr12r33
     .                         + r32r12*dyr32r13
     .                         + r22r32*dyr22r33
     .                         + r32r22*dyr32r23)
     .             + dyreploc(21)*(r12r12*r22r23
     .                         + r22r22*r12r13)
     .             + reploc(21)*(dyr12r12*r22r23
     .                         + dyr22r22*r12r13)
     .             + reploc(21)*(r12r12*dyr22r23
     .                         + r22r22*dyr12r13)
     .             + dyreploc(22)*(r12r22*r12r23
     .                         + r22r12*r22r13)
     .             + reploc(22)*(dyr12r22*r12r23
     .                         + dyr22r12*r22r13)
     .             + reploc(22)*(r12r22*dyr12r23
     .                         + r22r12*dyr22r13)

        dzrepij(8,9) = dzreploc(16)*r32r32*r32r33
     .             + reploc(16)*dzr32r32*r32r33
     .             + reploc(16)*r32r32*dzr32r33
     .             + dzreploc(17)*radd04*r32r33
     .             + reploc(17)*dzradd04*r32r33
     .             + reploc(17)*radd04*dzr32r33
     .             + dzreploc(18)*radd05*r32r32
     .             + reploc(18)*dzradd05*r32r32
     .             + reploc(18)*radd05*dzr32r32
     .             + dzreploc(19)*(r12r12*r12r13
     .                         + r22r22*r22r23)
     .             + reploc(19)*(dzr12r12*r12r13
     .                         + dzr22r22*r22r23)
     .             + reploc(19)*(r12r12*dzr12r13
     .                         + r22r22*dzr22r23)
     .             + dzreploc(20)*(r12r32*r12r33
     .                         + r32r12*r32r13
     .                         + r22r32*r22r33
     .                         + r32r22*r32r23)
                dzrepij(8,9)=dzrepij(8,9)
     .             + reploc(20)*(dzr12r32*r12r33
     .                         + dzr32r12*r32r13
     .                         + dzr22r32*r22r33
     .                         + dzr32r22*r32r23)
     .             + reploc(20)*(r12r32*dzr12r33
     .                         + r32r12*dzr32r13
     .                         + r22r32*dzr22r33
     .                         + r32r22*dzr32r23)
     .             + dzreploc(21)*(r12r12*r22r23
     .                         + r22r22*r12r13)
     .             + reploc(21)*(dzr12r12*r22r23
     .                         + dzr22r22*r12r13)
     .             + reploc(21)*(r12r12*dzr22r23
     .                         + r22r22*dzr12r13)
     .             + dzreploc(22)*(r12r22*r12r23
     .                         + r22r12*r22r13)
     .             + reploc(22)*(dzr12r22*r12r23
     .                         + dzr22r12*r22r13)
     .             + reploc(22)*(r12r22*dzr12r23
     .                         + r22r12*dzr22r13)
!
!       (py pz | py py):
!

        repij(9,8) = reploc(16)*r32r33*r32r32
     .             + reploc(17)*radd05*r32r32
     .             + reploc(18)*radd04*r32r33
     .             + reploc(19)*(r12r13*r12r12
     .                         + r22r23*r22r22)
     .             + reploc(20)*(radd13*r12r32
     .                         + radd14*r22r32)
     .             + reploc(21)*(r12r13*r22r22
     .                         + r22r23*r12r12)
     .             + reploc(22)*radd16*r12r22

        dxrepij(9,8) = dxreploc(16)*r32r33*r32r32
     .             + reploc(16)*dxr32r33*r32r32
     .             + reploc(16)*r32r33*dxr32r32
     .             + dxreploc(17)*radd05*r32r32
     .             + reploc(17)*dxradd05*r32r32
     .             + reploc(17)*radd05*dxr32r32
     .             + dxreploc(18)*radd04*r32r33
     .             + reploc(18)*dxradd04*r32r33
     .             + reploc(18)*radd04*dxr32r33
     .             + dxreploc(19)*(r12r13*r12r12
     .                         + r22r23*r22r22)
     .             + reploc(19)*(dxr12r13*r12r12
     .                         + dxr22r23*r22r22)
     .             + reploc(19)*(r12r13*dxr12r12
     .                         + r22r23*dxr22r22)
              dxrepij(9,8)=dxrepij(9,8)
     .             + dxreploc(20)*(radd13*r12r32
     .                         + radd14*r22r32)
     .             + reploc(20)*(dxradd13*r12r32
     .                         + dxradd14*r22r32)
     .             + reploc(20)*(radd13*dxr12r32
     .                         + radd14*dxr22r32)
     .             + dxreploc(21)*(r12r13*r22r22
     .                         + r22r23*r12r12)
     .             + reploc(21)*(dxr12r13*r22r22
     .                         + dxr22r23*r12r12)
     .             + reploc(21)*(r12r13*dxr22r22
     .                         + r22r23*dxr12r12)
     .             + dxreploc(22)*radd16*r12r22
     .             + reploc(22)*dxradd16*r12r22
     .             + reploc(22)*radd16*dxr12r22

        dyrepij(9,8) = dyreploc(16)*r32r33*r32r32
     .             + reploc(16)*dyr32r33*r32r32
     .             + reploc(16)*r32r33*dyr32r32
     .             + dyreploc(17)*radd05*r32r32
     .             + reploc(17)*dyradd05*r32r32
     .             + reploc(17)*radd05*dyr32r32
     .             + dyreploc(18)*radd04*r32r33
     .             + reploc(18)*dyradd04*r32r33
     .             + reploc(18)*radd04*dyr32r33
     .             + dyreploc(19)*(r12r13*r12r12
     .                         + r22r23*r22r22)
     .             + reploc(19)*(dyr12r13*r12r12
     .                         + dyr22r23*r22r22)
     .             + reploc(19)*(r12r13*dyr12r12
     .                         + r22r23*dyr22r22)
              dyrepij(9,8)=dyrepij(9,8)
     .             + dyreploc(20)*(radd13*r12r32
     .                         + radd14*r22r32)
     .             + reploc(20)*(dyradd13*r12r32
     .                         + dyradd14*r22r32)
     .             + reploc(20)*(radd13*dyr12r32
     .                         + radd14*dyr22r32)
     .             + dyreploc(21)*(r12r13*r22r22
     .                         + r22r23*r12r12)
     .             + reploc(21)*(dyr12r13*r22r22
     .                         + dyr22r23*r12r12)
     .             + reploc(21)*(r12r13*dyr22r22
     .                         + r22r23*dyr12r12)
     .             + dyreploc(22)*radd16*r12r22
     .             + reploc(22)*dyradd16*r12r22
     .             + reploc(22)*radd16*dyr12r22

        dzrepij(9,8) = dzreploc(16)*r32r33*r32r32
     .             + reploc(16)*dzr32r33*r32r32
     .             + reploc(16)*r32r33*dzr32r32
     .             + dzreploc(17)*radd05*r32r32
     .             + reploc(17)*dzradd05*r32r32
     .             + reploc(17)*radd05*dzr32r32
     .             + dzreploc(18)*radd04*r32r33
     .             + reploc(18)*dzradd04*r32r33
     .             + reploc(18)*radd04*dzr32r33
     .             + dzreploc(19)*(r12r13*r12r12
     .                         + r22r23*r22r22)
     .             + reploc(19)*(dzr12r13*r12r12
     .                         + dzr22r23*r22r22)
     .             + reploc(19)*(r12r13*dzr12r12
     .                         + r22r23*dzr22r22)
              dzrepij(9,8)=dzrepij(9,8)
     .             + dzreploc(20)*(radd13*r12r32
     .                         + radd14*r22r32)
     .             + reploc(20)*(dzradd13*r12r32
     .                         + dzradd14*r22r32)
     .             + reploc(20)*(radd13*dzr12r32
     .                         + radd14*dzr22r32)
     .             + dzreploc(21)*(r12r13*r22r22
     .                         + r22r23*r12r12)
     .             + reploc(21)*(dzr12r13*r22r22
     .                         + dzr22r23*r12r12)
     .             + reploc(21)*(r12r13*dzr22r22
     .                         + r22r23*dzr12r12)
     .             + dzreploc(22)*radd16*r12r22
     .             + reploc(22)*dzradd16*r12r22
     .             + reploc(22)*radd16*dzr12r22


!
!       (py py | pz pz):
!
        repij(8,10) = reploc(16)*r32r32*r33r33
     .              + reploc(17)*radd04*r33r33
     .              + reploc(18)*radd06*r32r32
     .              + reploc(19)*(r12r12*r13r13
     .                          + r22r22*r23r23)
     .              + reploc(20)*(r12r32*r13r33
     .                          + r22r32*r23r33)
     .              + reploc(21)*(r12r12*r23r23
     .                          + r22r22*r13r13)
     .              + reploc(22)*r12r22*r13r23

        dxrepij(8,10) = dxreploc(16)*r32r32*r33r33
     .              + reploc(16)*dxr32r32*r33r33
     .              + reploc(16)*r32r32*dxr33r33
     .              + dxreploc(17)*radd04*r33r33
     .              + reploc(17)*dxradd04*r33r33
     .              + reploc(17)*radd04*dxr33r33
     .              + dxreploc(18)*radd06*r32r32
     .              + reploc(18)*dxradd06*r32r32
     .              + reploc(18)*radd06*dxr32r32
     .              + dxreploc(19)*(r12r12*r13r13
     .                          + r22r22*r23r23)
     .              + reploc(19)*(dxr12r12*r13r13
     .                          + dxr22r22*r23r23)
     .              + reploc(19)*(r12r12*dxr13r13
     .                          + r22r22*dxr23r23)
               dxrepij(8,10)=dxrepij(8,10)
     .              + dxreploc(20)*(r12r32*r13r33
     .                          + r22r32*r23r33)
     .              + reploc(20)*(dxr12r32*r13r33
     .                          + dxr22r32*r23r33)
     .              + reploc(20)*(r12r32*dxr13r33
     .                          + r22r32*dxr23r33)
     .              + dxreploc(21)*(r12r12*r23r23
     .                          + r22r22*r13r13)
     .              + reploc(21)*(dxr12r12*r23r23
     .                          + dxr22r22*r13r13)
     .              + reploc(21)*(r12r12*dxr23r23
     .                          + r22r22*dxr13r13)
     .              + dxreploc(22)*r12r22*r13r23
     .              + reploc(22)*dxr12r22*r13r23
     .              + reploc(22)*r12r22*dxr13r23

        dyrepij(8,10) = dyreploc(16)*r32r32*r33r33
     .              + reploc(16)*dyr32r32*r33r33
     .              + reploc(16)*r32r32*dyr33r33
     .              + dyreploc(17)*radd04*r33r33
     .              + reploc(17)*dyradd04*r33r33
     .              + reploc(17)*radd04*dyr33r33
     .              + dyreploc(18)*radd06*r32r32
     .              + reploc(18)*dyradd06*r32r32
     .              + reploc(18)*radd06*dyr32r32
     .              + dyreploc(19)*(r12r12*r13r13
     .                          + r22r22*r23r23)
     .              + reploc(19)*(dyr12r12*r13r13
     .                          + dyr22r22*r23r23)
     .              + reploc(19)*(r12r12*dyr13r13
     .                          + r22r22*dyr23r23)
               dyrepij(8,10)=dyrepij(8,10)
     .              + dyreploc(20)*(r12r32*r13r33
     .                          + r22r32*r23r33)
     .              + reploc(20)*(dyr12r32*r13r33
     .                          + dyr22r32*r23r33)
     .              + reploc(20)*(r12r32*dyr13r33
     .                          + r22r32*dyr23r33)
     .              + dyreploc(21)*(r12r12*r23r23
     .                          + r22r22*r13r13)
     .              + reploc(21)*(dyr12r12*r23r23
     .                          + dyr22r22*r13r13)
     .              + reploc(21)*(r12r12*dyr23r23
     .                          + r22r22*dyr13r13)
     .              + dyreploc(22)*r12r22*r13r23
     .              + reploc(22)*dyr12r22*r13r23
     .              + reploc(22)*r12r22*dyr13r23

        dzrepij(8,10) = dzreploc(16)*r32r32*r33r33
     .              + reploc(16)*dzr32r32*r33r33
     .              + reploc(16)*r32r32*dzr33r33
     .              + dzreploc(17)*radd04*r33r33
     .              + reploc(17)*dzradd04*r33r33
     .              + reploc(17)*radd04*dzr33r33
     .              + dzreploc(18)*radd06*r32r32
     .              + reploc(18)*dzradd06*r32r32
     .              + reploc(18)*radd06*dzr32r32
     .              + dzreploc(19)*(r12r12*r13r13
     .                          + r22r22*r23r23)
     .              + reploc(19)*(dzr12r12*r13r13
     .                          + dzr22r22*r23r23)
     .              + reploc(19)*(r12r12*dzr13r13
     .                          + r22r22*dzr23r23)
               dzrepij(8,10)=dzrepij(8,10)
     .              + dzreploc(20)*(r12r32*r13r33
     .                          + r22r32*r23r33)
     .              + reploc(20)*(dzr12r32*r13r33
     .                          + dzr22r32*r23r33)
     .              + reploc(20)*(r12r32*dzr13r33
     .                          + r22r32*dzr23r33)
     .              + dzreploc(21)*(r12r12*r23r23
     .                          + r22r22*r13r13)
     .              + reploc(21)*(dzr12r12*r23r23
     .                          + dzr22r22*r13r13)
     .              + reploc(21)*(r12r12*dzr23r23
     .                          + r22r22*dzr13r13)
     .              + dzreploc(22)*r12r22*r13r23
     .              + reploc(22)*dzr12r22*r13r23
     .              + reploc(22)*r12r22*dzr13r23
!
!       (pz pz | py py):
!

        repij(10,8) = reploc(16)*r33r33*r32r32
     .              + reploc(17)*radd06*r32r32
     .              + reploc(18)*radd04*r33r33
     .              + reploc(19)*(r13r13*r12r12
     .                          + r23r23*r22r22)
     .              + reploc(20)*(r13r33*r12r32
     .                          + r23r33*r22r32)
     .              + reploc(21)*(r13r13*r22r22
     .                          + r23r23*r12r12)
     .              + reploc(22)*r13r23*r12r22

        dxrepij(10,8) = dxreploc(16)*r33r33*r32r32
     .              + reploc(16)*dxr33r33*r32r32
     .              + reploc(16)*r33r33*dxr32r32
     .              + dxreploc(17)*radd06*r32r32
     .              + reploc(17)*dxradd06*r32r32
     .              + reploc(17)*radd06*dxr32r32
     .              + dxreploc(18)*radd04*r33r33
     .              + reploc(18)*dxradd04*r33r33
     .              + reploc(18)*radd04*dxr33r33
     .              + dxreploc(19)*(r13r13*r12r12
     .                          + r23r23*r22r22)
     .              + reploc(19)*(dxr13r13*r12r12
     .                          + dxr23r23*r22r22)
     .              + reploc(19)*(r13r13*dxr12r12
     .                          + r23r23*dxr22r22)
                   dxrepij(10,8)=dxrepij(10,8)
     .              + dxreploc(20)*(r13r33*r12r32
     .                          + r23r33*r22r32)
     .              + reploc(20)*(dxr13r33*r12r32
     .                          + dxr23r33*r22r32)
     .              + reploc(20)*(r13r33*dxr12r32
     .                          + r23r33*dxr22r32)
     .              + dxreploc(21)*(r13r13*r22r22
     .                          + r23r23*r12r12)
     .              + reploc(21)*(dxr13r13*r22r22
     .                          + dxr23r23*r12r12)
     .              + reploc(21)*(r13r13*dxr22r22
     .                          + r23r23*dxr12r12)
     .              + dxreploc(22)*r13r23*r12r22
     .              + reploc(22)*dxr13r23*r12r22
     .              + reploc(22)*r13r23*dxr12r22

        dyrepij(10,8) = dyreploc(16)*r33r33*r32r32
     .              + reploc(16)*dyr33r33*r32r32
     .              + reploc(16)*r33r33*dyr32r32
     .              + dyreploc(17)*radd06*r32r32
     .              + reploc(17)*dyradd06*r32r32
     .              + reploc(17)*radd06*dyr32r32
     .              + dyreploc(18)*radd04*r33r33
     .              + reploc(18)*dyradd04*r33r33
     .              + reploc(18)*radd04*dyr33r33
     .              + dyreploc(19)*(r13r13*r12r12
     .                          + r23r23*r22r22)
     .              + reploc(19)*(dyr13r13*r12r12
     .                          + dyr23r23*r22r22)
     .              + reploc(19)*(r13r13*dyr12r12
     .                          + r23r23*dyr22r22)
               dyrepij(10,8)=dyrepij(10,8)
     .              + dyreploc(20)*(r13r33*r12r32
     .                          + r23r33*r22r32)
     .              + reploc(20)*(dyr13r33*r12r32
     .                          + dyr23r33*r22r32)
     .              + reploc(20)*(r13r33*dyr12r32
     .                          + r23r33*dyr22r32)
     .              + dyreploc(21)*(r13r13*r22r22
     .                          + r23r23*r12r12)
     .              + reploc(21)*(dyr13r13*r22r22
     .                          + dyr23r23*r12r12)
     .              + reploc(21)*(r13r13*dyr22r22
     .                          + r23r23*dyr12r12)
     .              + dyreploc(22)*r13r23*r12r22
     .              + reploc(22)*dyr13r23*r12r22
     .              + reploc(22)*r13r23*dyr12r22


        dzrepij(10,8) = dzreploc(16)*r33r33*r32r32
     .              + reploc(16)*dzr33r33*r32r32
     .              + reploc(16)*r33r33*dzr32r32
     .              + dzreploc(17)*radd06*r32r32
     .              + reploc(17)*dzradd06*r32r32
     .              + reploc(17)*radd06*dzr32r32
     .              + dzreploc(18)*radd04*r33r33
     .              + reploc(18)*dzradd04*r33r33
     .              + reploc(18)*radd04*dzr33r33
     .              + dzreploc(19)*(r13r13*r12r12
     .                          + r23r23*r22r22)
     .              + reploc(19)*(dzr13r13*r12r12
     .                          + dzr23r23*r22r22)
     .              + reploc(19)*(r13r13*dzr12r12
     .                          + r23r23*dzr22r22)
                 dzrepij(10,8)=dzrepij(10,8)
     .              + dzreploc(20)*(r13r33*r12r32
     .                          + r23r33*r22r32)
     .              + reploc(20)*(dzr13r33*r12r32
     .                          + dzr23r33*r22r32)
     .              + reploc(20)*(r13r33*dzr12r32
     .                          + r23r33*dzr22r32)
     .              + dzreploc(21)*(r13r13*r22r22
     .                          + r23r23*r12r12)
     .              + reploc(21)*(dzr13r13*r22r22
     .                          + dzr23r23*r12r12)
     .              + reploc(21)*(r13r13*dzr22r22
     .                          + r23r23*dzr12r12)
     .              + dzreploc(22)*r13r23*r12r22
     .              + reploc(22)*dzr13r23*r12r22
     .              + reploc(22)*r13r23*dzr12r22


!
!       (py pz | py pz):
!
        repij(9,9) = reploc(16)*r32r33*r32r33
     .             + reploc(17)*radd05*r32r33
     .             + reploc(18)*radd05*r32r33
     .             + reploc(19)*(r12r13*r12r13
     .                         + r22r23*r22r23)
     .             + reploc(20)*(radd13**2
     .                         + radd14**2)
     .             + reploc(21)*(r12r13*r22r23)*2.0d0
     .             + reploc(22)*radd16**2

        dxrepij(9,9) = dxreploc(16)*r32r33*r32r33
     .             + reploc(16)*dxr32r33*r32r33
     .             + reploc(16)*r32r33*dxr32r33
     .             + dxreploc(17)*radd05*r32r33
     .             + reploc(17)*dxradd05*r32r33
     .             + reploc(17)*radd05*dxr32r33
     .             + dxreploc(18)*radd05*r32r33
     .             + reploc(18)*dxradd05*r32r33
     .             + reploc(18)*radd05*dxr32r33
     .             + dxreploc(19)*(r12r13*r12r13
     .                         + r22r23*r22r23)
     .             + reploc(19)*(dxr12r13*r12r13
     .                         + dxr22r23*r22r23)
     .             + reploc(19)*(r12r13*dxr12r13
     .                         + r22r23*dxr22r23)
                 dxrepij(9,9)=dxrepij(9,9)
     .             + dxreploc(20)*(radd13**2
     .                         + radd14**2)
     .             + reploc(20)*2.0d0*(radd13*dxradd13
     .                               + radd14*dxradd14)
     .             + dxreploc(21)*(r12r13*r22r23)*2.0d0
     .             + reploc(21)*(dxr12r13*r22r23)*2.0d0
     .             + reploc(21)*(r12r13*dxr22r23)*2.0d0
     .             + dxreploc(22)*radd16**2
     .             + reploc(22)*radd16*dxradd16*2.0d0

        dyrepij(9,9) = dyreploc(16)*r32r33*r32r33
     .             + reploc(16)*dyr32r33*r32r33
     .             + reploc(16)*r32r33*dyr32r33
     .             + dyreploc(17)*radd05*r32r33
     .             + reploc(17)*dyradd05*r32r33
     .             + reploc(17)*radd05*dyr32r33
     .             + dyreploc(18)*radd05*r32r33
     .             + reploc(18)*dyradd05*r32r33
     .             + reploc(18)*radd05*dyr32r33
     .             + dyreploc(19)*(r12r13*r12r13
     .                         + r22r23*r22r23)
     .             + reploc(19)*(dyr12r13*r12r13
     .                         + dyr22r23*r22r23)
     .             + reploc(19)*(r12r13*dyr12r13
     .                         + r22r23*dyr22r23)
                 dyrepij(9,9)=dyrepij(9,9)
     .             + dyreploc(20)*(radd13**2
     .                         + radd14**2)
     .             + reploc(20)*2.0d0*(radd13*dyradd13
     .                         + radd14*dyradd14)
     .             + dyreploc(21)*(r12r13*r22r23)*2.0d0
     .             + reploc(21)*(dyr12r13*r22r23)*2.0d0
     .             + reploc(21)*(r12r13*dyr22r23)*2.0d0
     .             + dyreploc(22)*radd16**2
     .             + reploc(22)*radd16*dyradd16*2.0d0

        dzrepij(9,9) = dzreploc(16)*r32r33*r32r33
     .             + reploc(16)*dzr32r33*r32r33
     .             + reploc(16)*r32r33*dzr32r33
     .             + dzreploc(17)*radd05*r32r33
     .             + reploc(17)*dzradd05*r32r33
     .             + reploc(17)*radd05*dzr32r33
     .             + dzreploc(18)*radd05*r32r33
     .             + reploc(18)*dzradd05*r32r33
     .             + reploc(18)*radd05*dzr32r33
     .             + dzreploc(19)*(r12r13*r12r13
     .                         + r22r23*r22r23)
     .             + reploc(19)*(dzr12r13*r12r13
     .                         + dzr22r23*r22r23)
     .             + reploc(19)*(r12r13*dzr12r13
     .                         + r22r23*dzr22r23)
                 dzrepij(9,9)=dzrepij(9,9)
     .             + dzreploc(20)*(radd13**2
     .                         + radd14**2)
     .             + reploc(20)*2.0d0*(radd13*dzradd13
     .                         + radd14*dzradd14)
     .             + dzreploc(21)*(r12r13*r22r23)*2.0d0
     .             + reploc(21)*(dzr12r13*r22r23)*2.0d0
     .             + reploc(21)*(r12r13*dzr22r23)*2.0d0
     .             + dzreploc(22)*radd16**2
     .             + reploc(22)*radd16*dzradd16*2.0d0
!
!       (py pz | pz pz):
!
        repij(9,10) = reploc(16)*r32r33*r33r33
     .              + reploc(17)*radd05*r33r33
     .              + reploc(18)*radd06*r32r33
     .              + reploc(19)*(r12r13*r13r13
     .                          + r22r23*r23r23)
     .              + reploc(20)*(radd13*r13r33
     .                          + radd14*r23r33)
     .              + reploc(21)*(r12r13*r23r23
     .                          + r22r23*r13r13)
     .              + reploc(22)*radd16*r13r23

        dxrepij(9,10) = dxreploc(16)*r32r33*r33r33
     .              + reploc(16)*dxr32r33*r33r33
     .              + reploc(16)*r32r33*dxr33r33
     .              + dxreploc(17)*radd05*r33r33
     .              + reploc(17)*dxradd05*r33r33
     .              + reploc(17)*radd05*dxr33r33
     .              + dxreploc(18)*radd06*r32r33
     .              + reploc(18)*dxradd06*r32r33
     .              + reploc(18)*radd06*dxr32r33
     .              + dxreploc(19)*(r12r13*r13r13
     .                          + r22r23*r23r23)
     .              + reploc(19)*(dxr12r13*r13r13
     .                          + dxr22r23*r23r23)
     .              + reploc(19)*(r12r13*dxr13r13
     .                          + r22r23*dxr23r23)
               dxrepij(9,10)=dxrepij(9,10)
     .              + dxreploc(20)*(radd13*r13r33
     .                          + radd14*r23r33)
     .              + reploc(20)*(dxradd13*r13r33
     .                          + dxradd14*r23r33)
     .              + reploc(20)*(radd13*dxr13r33
     .                          + radd14*dxr23r33)
     .              + dxreploc(21)*(r12r13*r23r23
     .                          + r22r23*r13r13)
     .              + reploc(21)*(dxr12r13*r23r23
     .                          + dxr22r23*r13r13)
     .              + reploc(21)*(r12r13*dxr23r23
     .                          + r22r23*dxr13r13)
     .              + dxreploc(22)*radd16*r13r23
     .              + reploc(22)*dxradd16*r13r23
     .              + reploc(22)*radd16*dxr13r23

        dyrepij(9,10) = dyreploc(16)*r32r33*r33r33
     .              + reploc(16)*dyr32r33*r33r33
     .              + reploc(16)*r32r33*dyr33r33
     .              + dyreploc(17)*radd05*r33r33
     .              + reploc(17)*dyradd05*r33r33
     .              + reploc(17)*radd05*dyr33r33
     .              + dyreploc(18)*radd06*r32r33
     .              + reploc(18)*dyradd06*r32r33
     .              + reploc(18)*radd06*dyr32r33
     .              + dyreploc(19)*(r12r13*r13r13
     .                          + r22r23*r23r23)
     .              + reploc(19)*(dyr12r13*r13r13
     .                          + dyr22r23*r23r23)
     .              + reploc(19)*(r12r13*dyr13r13
     .                          + r22r23*dyr23r23)
               dyrepij(9,10)=dyrepij(9,10)
     .              + dyreploc(20)*(radd13*r13r33
     .                          + radd14*r23r33)
     .              + reploc(20)*(dyradd13*r13r33
     .                          + dyradd14*r23r33)
     .              + reploc(20)*(radd13*dyr13r33
     .                          + radd14*dyr23r33)
     .              + dyreploc(21)*(r12r13*r23r23
     .                          + r22r23*r13r13)
     .              + reploc(21)*(dyr12r13*r23r23
     .                          + dyr22r23*r13r13)
     .              + reploc(21)*(r12r13*dyr23r23
     .                          + r22r23*dyr13r13)
     .              + dyreploc(22)*radd16*r13r23
     .              + reploc(22)*dyradd16*r13r23
     .              + reploc(22)*radd16*dyr13r23

        dzrepij(9,10) = dzreploc(16)*r32r33*r33r33
     .              + reploc(16)*dzr32r33*r33r33
     .              + reploc(16)*r32r33*dzr33r33
     .              + dzreploc(17)*radd05*r33r33
     .              + reploc(17)*dzradd05*r33r33
     .              + reploc(17)*radd05*dzr33r33
     .              + dzreploc(18)*radd06*r32r33
     .              + reploc(18)*dzradd06*r32r33
     .              + reploc(18)*radd06*dzr32r33
     .              + dzreploc(19)*(r12r13*r13r13
     .                          + r22r23*r23r23)
     .              + reploc(19)*(dzr12r13*r13r13
     .                          + dzr22r23*r23r23)
     .              + reploc(19)*(r12r13*dzr13r13
     .                          + r22r23*dzr23r23)
               dzrepij(9,10)=dzrepij(9,10)
     .              + dzreploc(20)*(radd13*r13r33
     .                          + radd14*r23r33)
     .              + reploc(20)*(dzradd13*r13r33
     .                          + dzradd14*r23r33)
     .              + reploc(20)*(radd13*dzr13r33
     .                          + radd14*dzr23r33)
     .              + dzreploc(21)*(r12r13*r23r23
     .                          + r22r23*r13r13)
     .              + reploc(21)*(dzr12r13*r23r23
     .                          + dzr22r23*r13r13)
     .              + reploc(21)*(r12r13*dzr23r23
     .                          + r22r23*dzr13r13)
     .              + dzreploc(22)*radd16*r13r23
     .              + reploc(22)*dzradd16*r13r23
     .              + reploc(22)*radd16*dzr13r23
!
!       (pz pz | py pz):
!
        repij(10,9) = reploc(16)*r33r33*r32r33
     .              + reploc(17)*radd06*r32r33
     .              + reploc(18)*radd05*r33r33
     .              + reploc(19)*(r13r13*r12r13
     .                          + r23r23*r22r23)
     .              + reploc(20)*(r13r33*r12r33
     .                          + r33r13*r32r13
     .                          + r23r33*r22r33
     .                          + r33r23*r32r23)
     .              + reploc(21)*(r13r13*r22r23
     .                          + r23r23*r12r13)
     .              + reploc(22)*(r13r23*r12r23
     .                          + r23r13*r22r13)

        dxrepij(10,9) = dxreploc(16)*r33r33*r32r33
     .              + reploc(16)*dxr33r33*r32r33
     .              + reploc(16)*r33r33*dxr32r33
     .              + dxreploc(17)*radd06*r32r33
     .              + reploc(17)*dxradd06*r32r33
     .              + reploc(17)*radd06*dxr32r33
     .              + dxreploc(18)*radd05*r33r33
     .              + reploc(18)*dxradd05*r33r33
     .              + reploc(18)*radd05*dxr33r33
     .              + dxreploc(19)*(r13r13*r12r13
     .                          + r23r23*r22r23)
     .              + reploc(19)*(dxr13r13*r12r13
     .                          + dxr23r23*r22r23)
     .              + reploc(19)*(r13r13*dxr12r13
     .                          + r23r23*dxr22r23)
              dxrepij(10,9)=dxrepij(10,9)
     .              + dxreploc(20)*(r13r33*r12r33
     .                          + r33r13*r32r13
     .                          + r23r33*r22r33
     .                          + r33r23*r32r23)
     .              + reploc(20)*(dxr13r33*r12r33
     .                          + dxr33r13*r32r13
     .                          + dxr23r33*r22r33
     .                          + dxr33r23*r32r23)
     .              + reploc(20)*(r13r33*dxr12r33
     .                          + r33r13*dxr32r13
     .                          + r23r33*dxr22r33
     .                          + r33r23*dxr32r23)
     .              + dxreploc(21)*(r13r13*r22r23
     .                          + r23r23*r12r13)
     .              + reploc(21)*(dxr13r13*r22r23
     .                          + dxr23r23*r12r13)
     .              + reploc(21)*(r13r13*dxr22r23
     .                          + r23r23*dxr12r13)
     .              + dxreploc(22)*(r13r23*r12r23
     .                          + r23r13*r22r13)
     .              + reploc(22)*(dxr13r23*r12r23
     .                          + dxr23r13*r22r13)
     .              + reploc(22)*(r13r23*dxr12r23
     .                          + r23r13*dxr22r13)

        dyrepij(10,9) = dyreploc(16)*r33r33*r32r33
     .              + reploc(16)*dyr33r33*r32r33
     .              + reploc(16)*r33r33*dyr32r33
     .              + dyreploc(17)*radd06*r32r33
     .              + reploc(17)*dyradd06*r32r33
     .              + reploc(17)*radd06*dyr32r33
     .              + dyreploc(18)*radd05*r33r33
     .              + reploc(18)*dyradd05*r33r33
     .              + reploc(18)*radd05*dyr33r33
     .              + dyreploc(19)*(r13r13*r12r13
     .                          + r23r23*r22r23)
     .              + reploc(19)*(dyr13r13*r12r13
     .                          + dyr23r23*r22r23)
     .              + reploc(19)*(r13r13*dyr12r13
     .                          + r23r23*dyr22r23)
              dyrepij(10,9)=dyrepij(10,9)
     .              + dyreploc(20)*(r13r33*r12r33
     .                          + r33r13*r32r13
     .                          + r23r33*r22r33
     .                          + r33r23*r32r23)
     .              + reploc(20)*(dyr13r33*r12r33
     .                          + dyr33r13*r32r13
     .                          + dyr23r33*r22r33
     .                          + dyr33r23*r32r23)
     .              + reploc(20)*(r13r33*dyr12r33
     .                          + r33r13*dyr32r13
     .                          + r23r33*dyr22r33
     .                          + r33r23*dyr32r23)
     .              + dyreploc(21)*(r13r13*r22r23
     .                          + r23r23*r12r13)
     .              + reploc(21)*(dyr13r13*r22r23
     .                          + dyr23r23*r12r13)
     .              + reploc(21)*(r13r13*dyr22r23
     .                          + r23r23*dyr12r13)
     .              + dyreploc(22)*(r13r23*r12r23
     .                          + r23r13*r22r13)
     .              + reploc(22)*(dyr13r23*r12r23
     .                          + dyr23r13*r22r13)
     .              + reploc(22)*(r13r23*dyr12r23
     .                          + r23r13*dyr22r13)

        dzrepij(10,9) = dzreploc(16)*r33r33*r32r33
     .              + reploc(16)*dzr33r33*r32r33
     .              + reploc(16)*r33r33*dzr32r33
     .              + dzreploc(17)*radd06*r32r33
     .              + reploc(17)*dzradd06*r32r33
     .              + reploc(17)*radd06*dzr32r33
     .              + dzreploc(18)*radd05*r33r33
     .              + reploc(18)*dzradd05*r33r33
     .              + reploc(18)*radd05*dzr33r33
     .              + dzreploc(19)*(r13r13*r12r13
     .                          + r23r23*r22r23)
     .              + reploc(19)*(dzr13r13*r12r13
     .                          + dzr23r23*r22r23)
     .              + reploc(19)*(r13r13*dzr12r13
     .                          + r23r23*dzr22r23)
               dzrepij(10,9)=dzrepij(10,9)
     .              + dzreploc(20)*(r13r33*r12r33
     .                          + r33r13*r32r13
     .                          + r23r33*r22r33
     .                          + r33r23*r32r23)
     .              + reploc(20)*(dzr13r33*r12r33
     .                          + dzr33r13*r32r13
     .                          + dzr23r33*r22r33
     .                          + dzr33r23*r32r23)
     .              + reploc(20)*(r13r33*dzr12r33
     .                          + r33r13*dzr32r13
     .                          + r23r33*dzr22r33
     .                          + r33r23*dzr32r23)
     .              + dzreploc(21)*(r13r13*r22r23
     .                          + r23r23*r12r13)
     .              + reploc(21)*(dzr13r13*r22r23
     .                          + dzr23r23*r12r13)
     .              + reploc(21)*(r13r13*dzr22r23
     .                          + r23r23*dzr12r13)
     .              + dzreploc(22)*(r13r23*r12r23
     .                          + r23r13*r22r13)
     .              + reploc(22)*(dzr13r23*r12r23
     .                          + dzr23r13*r22r13)
     .              + reploc(22)*(r13r23*dzr12r23
     .                          + r23r13*dzr22r13)

!
!       (pz pz | pz pz):
!
        repij(10,10) = reploc(16)*r33r33*r33r33
     .               + reploc(17)*radd06*r33r33
     .               + reploc(18)*radd06*r33r33
     .               + reploc(19)*(r13r13*r13r13
     .                           + r23r23*r23r23)
     .               + reploc(20)*(r13r33*r13r33
     .                           + r23r33*r23r33)
     .               + reploc(21)*r13r13*r23r23*2.0d0
     .               + reploc(22)*r13r23*r13r23

        dxrepij(10,10) = dxreploc(16)*r33r33*r33r33
     .               + reploc(16)*dxr33r33*r33r33
     .               + reploc(16)*r33r33*dxr33r33
     .               + dxreploc(17)*radd06*r33r33
     .               + reploc(17)*dxradd06*r33r33
     .               + reploc(17)*radd06*dxr33r33
     .               + dxreploc(18)*radd06*r33r33
     .               + reploc(18)*dxradd06*r33r33
     .               + reploc(18)*radd06*dxr33r33
     .               + dxreploc(19)*(r13r13*r13r13
     .                           + r23r23*r23r23)
     .               + reploc(19)*(dxr13r13*r13r13
     .                           + dxr23r23*r23r23)
     .               + reploc(19)*(r13r13*dxr13r13
     .                           + r23r23*dxr23r23)
              dxrepij(10,10) =dxrepij(10,10) 
     .               + dxreploc(20)*(r13r33*r13r33
     .                           + r23r33*r23r33)
     .               + reploc(20)*(dxr13r33*r13r33
     .                           + dxr23r33*r23r33)
     .               + reploc(20)*(r13r33*dxr13r33
     .                           + r23r33*dxr23r33)
     .               + dxreploc(21)*r13r13*r23r23*2.0d0
     .               + reploc(21)*dxr13r13*r23r23*2.0d0
     .               + reploc(21)*r13r13*dxr23r23*2.0d0
     .               + dxreploc(22)*r13r23*r13r23
     .               + reploc(22)*dxr13r23*r13r23
     .               + reploc(22)*r13r23*dxr13r23


        dyrepij(10,10) = dyreploc(16)*r33r33*r33r33
     .               + reploc(16)*dyr33r33*r33r33
     .               + reploc(16)*r33r33*dyr33r33
     .               + dyreploc(17)*radd06*r33r33
     .               + reploc(17)*dyradd06*r33r33
     .               + reploc(17)*radd06*dyr33r33
     .               + dyreploc(18)*radd06*r33r33
     .               + reploc(18)*dyradd06*r33r33
     .               + reploc(18)*radd06*dyr33r33
     .               + dyreploc(19)*(r13r13*r13r13
     .                           + r23r23*r23r23)
     .               + reploc(19)*(dyr13r13*r13r13
     .                           + dyr23r23*r23r23)
     .               + reploc(19)*(r13r13*dyr13r13
     .                           + r23r23*dyr23r23)
              dyrepij(10,10) =dyrepij(10,10) 
     .               + dyreploc(20)*(r13r33*r13r33
     .                           + r23r33*r23r33)
     .               + reploc(20)*(dyr13r33*r13r33
     .                           + dyr23r33*r23r33)
     .               + reploc(20)*(r13r33*dyr13r33
     .                           + r23r33*dyr23r33)
     .               + dyreploc(21)*r13r13*r23r23*2.0d0
     .               + reploc(21)*dyr13r13*r23r23*2.0d0
     .               + reploc(21)*r13r13*dyr23r23*2.0d0
     .               + dyreploc(22)*r13r23*r13r23
     .               + reploc(22)*dyr13r23*r13r23
     .               + reploc(22)*r13r23*dyr13r23

        dzrepij(10,10) = dzreploc(16)*r33r33*r33r33
     .               + reploc(16)*dzr33r33*r33r33
     .               + reploc(16)*r33r33*dzr33r33
     .               + dzreploc(17)*radd06*r33r33
     .               + reploc(17)*dzradd06*r33r33
     .               + reploc(17)*radd06*dzr33r33
     .               + dzreploc(18)*radd06*r33r33
     .               + reploc(18)*dzradd06*r33r33
     .               + reploc(18)*radd06*dzr33r33
     .               + dzreploc(19)*(r13r13*r13r13
     .                           + r23r23*r23r23)
     .               + reploc(19)*(dzr13r13*r13r13
     .                           + dzr23r23*r23r23)
     .               + reploc(19)*(r13r13*dzr13r13
     .                           + r23r23*dzr23r23)
              dzrepij(10,10) =dzrepij(10,10) 
     .               + dzreploc(20)*(r13r33*r13r33
     .                           + r23r33*r23r33)
     .               + reploc(20)*(dzr13r33*r13r33
     .                           + dzr23r33*r23r33)
     .               + reploc(20)*(r13r33*dzr13r33
     .                           + r23r33*dzr23r33)
     .               + dzreploc(21)*r13r13*r23r23*2.0d0
     .               + reploc(21)*dzr13r13*r23r23*2.0d0
     .               + reploc(21)*r13r13*dzr23r23*2.0d0
     .               + dzreploc(22)*r13r23*r13r23
     .               + reploc(22)*dzr13r23*r13r23
     .               + reploc(22)*r13r23*dzr13r23
      endif
      return
      end
